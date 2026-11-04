"""Shared loading for the three CD20 replication fields (Zenodo 7795826):
EDF10/rep2, Figure4/untreated, Figure4/RTX treated.

Same study and instrument class as the Figure2 and EDF10/rep1 fields, so
the frozen conventions apply verbatim: pixel 130 nm (asserted from each
YAML), lp <= 8 nm with lp = 130*sqrt((lpx^2+lpy^2)/2), event linkage 10 nm,
16 quality strata of 8 nm^2, q = {0.10, 0.15} nm^-1, four angles.

ONE GENERALIZATION, disclosed in design_expectations.json before analysis:
panels are per-round frame quarters with 20-frame guards computed from each
round's own YAML frame count, because both Figure4 fields have
round-unequal frame counts (untreated R2: 34,579 vs 40,000; RTX R3/R4:
33,458/34,841 vs 40,000). Panels are within-round time bins — the
estimator compares panels within a round's windows and never aligns
absolute frames across rounds — so this is the natural extension of the
shared-panel convention, not a new degree of freedom.

Group centers follow the SDC/EDF10 convention verbatim: union-find
single-linkage at 70 nm over the study's own merged RESI reconstruction,
components with >= 5 averaged sources, localizations assigned within
120 nm. On membrane-receptor fields this percolates into density-selected
neighborhoods (EDF10/rep1 finding, expected to repeat, worst on RTX);
the 250-nm grid variant in analyze_field.py is the robustness check.
All code and oracle contrasts are invariant to per-group translation.

Data are read from the unrar extraction (CRC-verified; the manifest
records dual-toolchain SHA-256 agreement for every member).
"""
from pathlib import Path
import os,re
import numpy as np
import h5py

P=Path(__file__).resolve().parent
SCRATCH=Path(os.environ.get('SCRATCH',
    str(P.parents[1]/'scratch')))/'replcd20'
FIELDS={'edf10rep2':'EDF10/rep2/',
        'fig4untreated':'Figure4/untreated/',
        'fig4rtx':'Figure4/RTX treated/'}
PX=130.0
LP_MAX=8.0
ANGLES=np.arange(4)*np.pi/4
CLUSTER_NM=70.
MIN_SOURCES=5
RADIUS=120.

def paths(tag):
    root=SCRATCH/tag/'extract_unrar'/FIELDS[tag].rstrip('/')
    return dict(rounds=root/'RESI/RESI rounds',
        recon=root/'RESI/RESI results/RESI_R1,R2,R3,R4_multi.hdf5')

def frames_from_yaml(tag,r):
    txt=(paths(tag)['rounds']/f'RESI_R{r}_apicked.yaml').read_text()
    m=re.search(r'^Frames:\s*(\d+)',txt,re.M)
    px=re.search(r'^Pixelsize:\s*(\d+)',txt,re.M)
    assert px and float(px.group(1))==PX,('pixel size mismatch',tag,r)
    return int(m.group(1))

def make_panels(nframes,guard=10):
    b=[round(nframes*k/4) for k in (1,2,3)]
    return [(0,b[0]-guard),(b[0]+guard,b[1]-guard),
            (b[1]+guard,b[2]-guard),(b[2]+guard,nframes)]

def round_panels(tag):
    """Per-round frame-quarter panels from each round's own YAML count."""
    return {r:make_panels(frames_from_yaml(tag,r)) for r in (1,2,3,4)}

def centers_from_recon(tag):
    """Union-find single-linkage at CLUSTER_NM on merged-reconstruction
    positions; centers of components with >= MIN_SOURCES averaged sources.
    Returns (centers, sizes, n_recon_sources, frac_in_components)."""
    with h5py.File(paths(tag)['recon']) as f:a=f['locs'][:]
    xy=np.column_stack([a['x'],a['y']])*PX
    lab=np.arange(len(xy))
    def find(i):
        while lab[i]!=i:lab[i]=lab[lab[i]];i=lab[i]
        return i
    cell=(xy//CLUSTER_NM).astype(int)
    buckets={}
    for i,c in enumerate(map(tuple,cell)):buckets.setdefault(c,[]).append(i)
    for i,c in enumerate(map(tuple,cell)):
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for j in buckets.get((c[0]+dx,c[1]+dy),[]):
                    if j>i and ((xy[i]-xy[j])**2).sum()<CLUSTER_NM**2:
                        lab[find(i)]=find(j)
    roots=np.array([find(i) for i in range(len(xy))])
    _,inv=np.unique(roots,return_inverse=True)
    centers=[];sizes=[]
    for s in range(inv.max()+1):
        m=inv==s
        if m.sum()>=MIN_SOURCES:centers.append(xy[m].mean(0));sizes.append(int(m.sum()))
    sizes=np.array(sizes)
    return np.array(centers),sizes,len(xy),float(sizes.sum()/len(xy))

def _grid_map(centers,radius):
    m={}
    cc=np.floor(centers/radius).astype(np.int64)
    for i,(cx,cy) in enumerate(cc):
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                m.setdefault((cx+dx,cy+dy),[]).append(i)
    return {k:np.array(v) for k,v in m.items()}

def assign(xy,centers,gmap=None,radius=RADIUS):
    if gmap is None:gmap=_grid_map(centers,radius)
    best=np.full(len(xy),-1,np.int64)
    cells=np.floor(xy/radius).astype(np.int64)
    code=cells[:,0]*(1<<32)+cells[:,1]
    order=np.argsort(code,kind='stable');code_s=code[order]
    starts=np.flatnonzero(np.concatenate([[True],code_s[1:]!=code_s[:-1]]))
    starts=np.append(starts,len(code_s))
    for a,b in zip(starts[:-1],starts[1:]):
        rows=order[a:b]
        cx,cy=cells[rows[0]]
        cand=gmap.get((cx,cy))
        if cand is None:continue
        d=((xy[rows,None,:]-centers[cand][None,:,:])**2).sum(-1)
        j=np.argmin(d,axis=1);dm=d[np.arange(len(rows)),j]
        ok=dm<radius**2
        best[rows[ok]]=cand[j[ok]]
    return best

def load_rounds(tag,centers):
    """Filtered per-round arrays with pick = assigned center id.
    Localizations beyond RADIUS of every center are dropped (count kept)."""
    gmap=_grid_map(centers,RADIUS)
    rounds={}
    for r in (1,2,3,4):
        with h5py.File(paths(tag)['rounds']/f'RESI_R{r}_apicked.hdf5') as f:
            a=f['locs'][:]
        lp=PX*np.sqrt((a['lpx']**2+a['lpy']**2)/2)
        keep=np.isfinite(lp)&(lp<=LP_MAX)
        a=a[keep];lp=lp[keep]
        xy=np.column_stack([a['x'],a['y']])*PX
        pick=assign(xy,centers,gmap)
        m=pick>=0
        rounds[r]=dict(frame=a['frame'][m].astype(np.int64),
            pick=pick[m],xy=xy[m],
            lpx=(a['lpx'][m]*PX).astype(np.float64),
            lpy=(a['lpy'][m]*PX).astype(np.float64),lp=lp[m],
            z=(a['z'][m].astype(np.float64) if 'z' in a.dtype.names else None),
            rows_filtered=int(keep.sum()),rows_assigned=int(m.sum()),
            rows_total=int(len(keep)))
    return rounds

def load_rounds_raw(tag):
    """lp-filtered per-round arrays WITHOUT center assignment."""
    rounds={}
    for r in (1,2,3,4):
        with h5py.File(paths(tag)['rounds']/f'RESI_R{r}_apicked.hdf5') as f:
            a=f['locs'][:]
        lp=PX*np.sqrt((a['lpx']**2+a['lpy']**2)/2)
        keep=np.isfinite(lp)&(lp<=LP_MAX)
        a=a[keep];lp=lp[keep]
        rounds[r]=dict(frame=a['frame'].astype(np.int64),
            xy=np.column_stack([a['x'],a['y']])*PX,
            lpx=(a['lpx']*PX).astype(np.float64),
            lpy=(a['lpy']*PX).astype(np.float64),lp=lp)
    return rounds

def tiles(centers,nx=4,ny=3):
    lo=centers.min(0);hi=centers.max(0)
    ix=np.minimum(((centers[:,0]-lo[0])/(hi[0]-lo[0])*nx).astype(int),nx-1)
    iy=np.minimum(((centers[:,1]-lo[1])/(hi[1]-lo[1])*ny).astype(int),ny-1)
    return ix*ny+iy
