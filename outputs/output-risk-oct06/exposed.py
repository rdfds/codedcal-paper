"""Exposed-data harness: prediction core never receives native calls or offsets."""
import argparse,gc,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent.parent.parent;HERE=Path(__file__).resolve().parent
STAR=ROOT/'outputs/starfinder-public-confirmation-oct06';BAR=ROOT/'outputs/barseq3-native-repair-oct05'
sys.path[:0]=[str(STAR),str(ROOT/'tool/src')]
import adapter as A
import tissue as T
import probe as P
import quality_catalog as Q
import local_trace as L
import volume_native as V
import recovery as R
from codedcal import output_risk as O

def read(p):return json.loads(Path(p).read_text())
def write(p,d):A.write(p,d)
def sha(p):return A.sha(p)
def as3(m):
    m=np.rint(m).astype(int);return np.c_[np.zeros(len(m),int),m] if m.shape[1]==2 else m
CASES=[(s,k,t) for s,fields,kinds in [('bar',[0,1],['zero','common','color']),('culture',[4],['zero','common','color']),('tissue',[4],['zero','common','color','outclass'])] for k in fields for t in kinds]
def name(case):return '_'.join(map(str,case))

def bank(case):
    source,index,kind=case;maps=[];labels=[];ordinary=[];base=[];provenance=[]
    def load(p):provenance.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)));return read(p)
    def add(label,m,conventional=False,incumbent=False):
        m=as3(m);relative=m-m[:1]
        if np.max(np.abs(relative[:,0]))>6 or np.max(np.abs(relative[:,1:]))>20:return
        token=relative.tobytes();keys=[(a-a[:1]).tobytes() for a in maps]
        if token in keys:i=keys.index(token);labels[i].append(label)
        else:i=len(maps);maps.append(m);labels.append([label])
        if conventional and i not in ordinary:ordinary.append(i)
        if incumbent and i not in base:base.append(i)
    if source=='bar':
        model=load(STAR/'bar-model'/f'{index}_{kind}-geometry.json');add('author_selected',model['maps'],True,True)
        for row in model['candidates']:add('author_'+row['name'],row['maps'],row['name'].endswith('_refined'))
        init=load(STAR/'initialization'/f'barseq3_{index}_{kind}_correlation.json')['diagnostics']
        for row in init['starts']:add('image_'+row['start'],row['maps'],True,True)
        old=load(BAR/'recovery'/f'field-{index}'/f'{kind}-geometry.json')
        for p,m in old['maps'].items():
            if p not in ['oracle']:add('old_'+p,m,True,p in ['barcode_model','population_template'])
        global_geo=load(BAR/'global-repair'/f'field-{index}'/f'{kind}-geometry.json');group=global_geo['search'];add('group_selected',global_geo['maps']['global_group'])
        seed=as3(model['maps'])
    elif source=='culture':
        geom=load(STAR/'probe'/f'{kind}-geometry.json');group=geom['diagnostics']['global']
        for p,m in geom['maps'].items():
            if p=='oracle':continue
            add(p,m,p!='global_group',p in ['author_multistart','population_template'])
        for row in geom['diagnostics']['author_models']:
            if 'maps' in row:add('author_'+row['name'],row['maps'],True)
        init=load(STAR/'initialization'/f'starfinder_{kind}_0_correlation.json')['diagnostics']
        for row in init['starts']:add('image_'+row['start'],row['maps'],True,True)
        three=load(STAR/'source-comparators'/f'{kind}__frame_template_3d__correlation.json');add('template3d',three['diagnostics']['maps_zyx'],True,True)
        seed=as3(geom['maps']['author_multistart'])
    else:
        geom=load(STAR/'fullfield-tissue'/f'{kind}-geometry.json');group=geom['diagnostics']['global_group']
        for p,m in geom['maps_zyx'].items():
            if p=='oracle':continue
            add(p,m,p!='global_group',p in ['author_multistart','template_2d_multistart','template_3d_multistart'])
        for row in geom['diagnostics']['template_2d_multistart']['starts']:add('image_'+row['start'],row['maps'],True)
        seed=as3(geom['maps_zyx']['author_multistart'])
        if kind=='color':
            vol=load(STAR/'volumetric-model/geometry.json');add('author_volume_selected',vol['maps'],True,True);seed=as3(vol['maps'])
            for row in vol['candidates']:add('author_volume_'+row['name'],row['maps'],row['name'].endswith('_refined'))
    common=np.asarray(group['common_maps'],int)
    for fam in group['families']:
        initial=np.zeros_like(common) if fam['name']=='zero' else common
        for partition in range(3):
            row=min((r for r in fam['diagnostics']['grid'] if r['partition']==partition),key=lambda r:(r['loss'],sum(t*t for t in r['shift_yx']),r['shift_yx']))
            sign=np.where(np.isin(np.arange(len(seed))%4,row['plus_colors']),1,-1)
            add(f"group_{fam['name']}_{partition}",initial+sign[:,None]*row['shift_yx'])
    for c in range(4):
        for axis in [1,2]:
            for sign in [-1,1]:
                m=seed.copy();m[np.arange(len(m))%4==c,axis]+=sign;add(f'local_color{c}_axis{axis}_{sign}',m)
    for r in range(1,len(seed)//4):
        for axis in [1,2]:
            for sign in [-1,1]:
                m=seed.copy();m[r*4:r*4+4,axis]+=4*sign;add(f'local_round{r}_axis{axis}_{sign}',m)
    assert base and ordinary;return np.stack(maps),labels,ordinary,base,provenance

def inputs(case):
    source,index,kind=case
    if source=='bar':
        base,support=R.inputs(index);base=base[:,None];support=support[:,None];book,negative,names=R.qualify.codebook();unique=list(dict.fromkeys(names));ids=np.array([unique.index(n) for n in names]);imposed=R.offsets(kind);refdir=BAR/'recovery'/f'field-{index}'
        refs={d:read(refdir/f'reference-{d}.json') for d in R.DECODERS};inputsha=[read(BAR/f'field-{index}'/f'round-{r}-demons.json')['output_sha256'] for r in range(10)]
    else:
        if source=='culture':
            path=A.CACHE/'local-trace-Position351/crop-4.npz';book,negative,names,_=A.book();imposed=P.shifts(kind);refs=read(STAR/'probe/references.json')
        else:
            path=T.CACHE/'fullfield-conditioned-crop-4.npz';book,_=T.book();negative=None;imposed=T.offsets(kind);refs={d:read(STAR/'fullfield-tissue/native'/f'4__{d}.json')['catalog'] for d in A.DECODERS}
        with np.load(path) as d:base=d['images'];support=d['support']
        expected=read(STAR/'fullfield-tissue/crop-4.json')['array_sha256'] if source=='tissue' else read(STAR/'quality-catalog/native/crop-4__correlation.json')['input_sha256'];assert sha(path)==expected
        ids=np.arange(book.shape[1]);inputsha=[sha(path)]
    observed=P.warp(base,imposed);alpha=P.warp(support,imposed);physical=book[:,:998] if source=='culture' else book
    return observed,alpha,book,physical,negative,ids,imposed,refs,inputsha

def subset(points,words,source):
    points=np.asarray(points,float).reshape(-1,3);words=np.asarray(words,int)
    keep=(points[:,0]>=40)&(points[:,0]<88)&(points[:,1]>=40)&(points[:,1]<88)
    if source!='bar':keep&=(points[:,2]>=8)&(points[:,2]<22)
    return dict(xyz=points[keep].tolist(),words=words[keep].tolist())

def decode(observed,alpha,m,book,source,decoder,ids,negative):
    out=T.move3d(observed,m);sup=T.move3d(alpha,m)
    if source=='bar':
        xy,w,runtime=R.decode(out[:,0],book,decoder);points=np.c_[xy,np.zeros(len(xy))];raw=dict(xyz=points.tolist(),words=w.tolist())
    else:
        if source=='culture':
            cuts={d:r['selected']['threshold'] for d,r in read(STAR/'quality-catalog/state.json')['groups'].items()};raw=P.catalog(out,sup,book,negative,decoder,cuts[decoder])
        else:raw=T.catalog(out,sup,book,decoder)
        points=np.asarray(raw['xyz']).reshape(-1,3);w=np.asarray(raw['words'],int)
    canonical=points-m[0][::-1];cat=subset(canonical,ids[w],source)
    negcount=None if negative is None else int(negative[w].sum());coding=len(w) if negative is None else len(w)-negcount
    proxy=None if negative is None else negcount/max(coding,1)*int((~negative).sum())/int(negative.sum())
    return cat,dict(full_catalog_calls=len(w),negative=negcount,blank_proxy=proxy,source_abstentions=raw.get('abstained'),source_unsupported=raw.get('unsupported_candidates'))

def decoder_confidence(observed,m,book,source,decoder,cat):
    out=T.move3d(observed,m)
    if source=='bar':
        data=out[:,0,None,16:112,16:112]
        if decoder=='correlation':density=R.qualify.bardensr.spot_calling.estimate_density_singleshot(data,book,.05)
        else:density,_=R.qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),book,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
        density=np.asarray(density).reshape(96,96,-1);points=np.asarray(cat['xyz']).reshape(-1,3)+m[0][::-1];xy=np.rint(points[:,:2]).astype(int)-16;valid=np.all((xy>=0)&(xy<96),axis=1);values=density[xy[valid,1],xy[valid,0]]
        ordered=np.sort(values,axis=1);margin=(ordered[:,-1]-ordered[:,-2])/np.maximum(ordered[:,-1],1e-20)
    else:
        points=np.rint(np.asarray(cat['xyz']).reshape(-1,3)+m[0][::-1]).astype(int);rounds=len(book)//4
        if not len(points):return None
        vals=np.stack([out[:,z-1:z+2,y-2:y+3,x-2:x+3].sum(axis=(1,2,3)) for x,y,z in points],axis=1);cube=vals.T.reshape(-1,rounds,4);cube/=np.sqrt((cube*cube).sum(axis=2,keepdims=True))+1e-6
        if decoder=='source_spot':ordered=np.sort(cube,axis=2);margin=(ordered[:,:,-1]-ordered[:,:,-2]).min(axis=1)
        else:
            data=cube.reshape(-1,len(book)).T[:,None,None,:]
            if decoder=='correlation':density=R.qualify.bardensr.spot_calling.estimate_density_singleshot(data,book,.05)
            else:density,_=R.qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),book,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
            ordered=np.sort(np.asarray(density).reshape(len(points),-1),axis=1);margin=(ordered[:,-1]-ordered[:,-2])/np.maximum(ordered[:,-1],1e-20)
    return float(np.median(margin)) if len(margin) else None

def guard():
    L.guard();T.guard();freeze=read(HERE/'freeze.json');assert read(HERE/'controls.json')['passed']
    for rel,digest in freeze['files'].items():assert sha(ROOT/rel)==digest;assert subprocess.check_output(['git','show','HEAD:'+rel],cwd=ROOT)==(ROOT/rel).read_bytes(),rel
    return sha(HERE/'freeze.json')

def run(case):
    freeze=guard();source,index,kind=case;dest=HERE/'cases'/name(case);dest.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    if (dest/'evaluation.json').exists():assert read(dest/'evaluation.json')['freeze_sha256']==freeze;print('RETAINED',name(case),flush=True);return
    observed,alpha,book,physical,negative,ids,imposed,refs,inputsha=inputs(case);maps,labels,ordinary,bases,provenance=bank(case)
    region=((0,1),(40,88),(40,88)) if source=='bar' else ((8,20),(40,88),(40,88));block=(1,12,12) if source=='bar' else (4,12,12)
    profile=O.profile_losses(observed,alpha,physical,maps,region,block);comp=O.compatibility(profile);write(dest/'profile.json',dict(freeze_sha256=freeze,input_sha256=inputsha,geometry_sources=provenance,maps_zyx=maps.tolist(),labels=labels,conventional_indices=ordinary,incumbent_indices=bases,common_coverage=profile['coverage'],residual=profile['residual'].tolist(),energy=profile['energy'].tolist(),block_counts=profile['block_counts'].tolist(),compatibility=comp))
    if comp['abstained']:
        write(dest/'evaluation.json',dict(freeze_sha256=freeze,stage='abstained',reason=comp['reason'],development=True,independent_confirmation=False));print('ABSTAIN',name(case),comp['reason'],flush=True);return
    incumbent=min(bases,key=lambda k:(comp['losses'][k],float(np.square(maps[k]).sum()),k));decoders=['correlation'] if source=='bar' else ['source_spot','correlation'];rows=[]
    for decoder in decoders:
        catalogs=[];quality=[]
        for k,m in enumerate(maps):
            path=dest/f'catalog-{decoder}-{k}.json'
            if path.exists():row=read(path);assert row['freeze_sha256']==freeze;cat=row['catalog'];q=row['quality']
            else:
                cat,q=decode(observed,alpha,m,book,source,decoder,ids,negative);write(path,dict(freeze_sha256=freeze,candidate=k,decoder=decoder,catalog=cat,quality=q))
            cat['negative_ids']=[] if negative is None else np.unique(ids[np.flatnonzero(negative)]).tolist()
            central_negative=int(np.isin(cat['words'],cat['negative_ids']).sum());central_coding=len(cat['words'])-central_negative
            q=q|dict(central_coding=central_coding,central_negative=central_negative,blank_proxy=None if negative is None else central_negative/max(central_coding,1)*int((~negative).sum())/int(negative.sum()))
            catalogs.append(cat);quality.append(q)
        risk=O.output_risk(catalogs,comp,incumbent,ordinary);sens={str(t):O.output_risk(catalogs,O.compatibility(profile,tolerance=t),incumbent,ordinary) for t in [0.,.005]};confidence=decoder_confidence(observed,maps[incumbent],book,source,decoder,catalogs[incumbent]);prediction=dict(decoder=decoder,incumbent=incumbent,incumbent_labels=labels[incumbent],risk=risk,tolerance_sensitivities=sens,baseline_fit_loss=comp['losses'][incumbent],baseline_fit_gap=comp['losses'][incumbent]-min(comp['losses']),baseline_decoder_confidence=confidence,baseline_quality=quality[incumbent],freeze_sha256=freeze);write(dest/f'prediction-{decoder}.json',prediction)
        # Evaluation-only truth appears only AFTER prediction is durably written.
        ref=refs[decoder];points=np.asarray(ref['xy' if source=='bar' else 'xyz']).reshape(-1,2 if source=='bar' else 3)
        if source=='bar':points=np.c_[points,np.zeros(len(points))]
        points=points+np.r_[imposed[0][::-1],0];reference=subset(points,ids[np.asarray(ref['words'],int)],source);reference['negative_ids']=catalogs[incumbent]['negative_ids'];actual=O.catalog_distance(reference,catalogs[incumbent]);failure=actual['identity_loss']>.05 or actual['composition_tv']>.05
        rows.append(dict(decoder=decoder,reference_calls=len(reference['words']),incumbent_calls=len(catalogs[incumbent]['words']),actual_consistency_error=actual,failed_consistency=failure,risk=risk,baseline_fit_loss=comp['losses'][incumbent],baseline_fit_gap=comp['losses'][incumbent]-min(comp['losses']),baseline_decoder_confidence=confidence,baseline_quality=quality[incumbent]))
        print('RISK',name(case),decoder,'calls',len(catalogs[incumbent]['words']),'flag',risk.get('flagged'),'truth_failure',failure,'loss',actual['identity_loss'],'risk',risk.get('identity_risk'),flush=True);gc.collect()
    write(dest/'evaluation.json',dict(stage='complete',case=name(case),source=source,field=index,kind=kind,freeze_sha256=freeze,candidates=len(maps),incumbent=incumbent,rows=rows,runtime_s=time.monotonic()-started,development=True,independent_confirmation=False,independent_RNA_truth=False))
if __name__=='__main__':
    cli=argparse.ArgumentParser();cli.add_argument('--case',type=int,choices=range(13));args=cli.parse_args()
    cases=CASES if args.case is None else [CASES[args.case]]
    for case in cases:run(case)
