"""Render a focused diagram and complete retained-data displays; no fitting."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BLUE, ORANGE, INK = '#2675a9', '#c75b38', '#1e2935'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
    'axes.titlesize': 12, 'axes.labelsize': 11, 'xtick.labelsize': 10,
    'ytick.labelsize': 10, 'legend.fontsize': 9, 'pdf.fonttype': 42,
    'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})


def save(fig, out, name):
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(out/f'{name}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
        if ext=='svg':
            target=out/f'{name}.{ext}'
            target.write_text('\n'.join(line.rstrip() for line in target.read_text().splitlines())+'\n')
    plt.close(fig)


def construction(out):
    fig = plt.figure(figsize=(10.2, 6.5), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[1.1, 1], hspace=.12, wspace=.12)
    ax = fig.add_subplot(grid[0, :]);ax.set(xlim=(0, 12), ylim=(0, 4));ax.axis('off')
    ax.set_title('a  Permanent source codes, unequal counts, disjoint panels', loc='left', weight='bold')
    rng = np.random.default_rng(72)  # Layout only; not simulated scientific outcomes.
    colors = [BLUE, ORANGE, BLUE, '#26906e']  # A collision is allowed and shown.
    counts = [[3+(s+p)%4 for s in range(4)] for p in range(4)]
    for p, label in enumerate(('a', 'b', 'c', 'd')):
        x0 = p*3+.1
        ax.add_patch(Rectangle((x0, .9), 2.6, 2.3, facecolor='#f3f6f8', edgecolor='#b8c2cc'))
        ax.text(x0+1.3, 3.45, f'Window {"A" if p<2 else "B"}: panel {label}', ha='center', fontsize=10)
        for s, (x,y) in enumerate(((.55,1.4),(1.8,1.35),(.65,2.5),(1.8,2.55))):
            n = counts[p][s]
            ax.scatter(x0+x+rng.normal(0,.055,n), y+rng.normal(0,.055,n), s=18, color=colors[s])
    ax.text(.1,.35,'One code per persistent source; different sources may share a code.', fontsize=11)
    ax = fig.add_subplot(grid[1, 0]);ax.axis('off')
    ax.set_title('b  Match the count product', loc='left', weight='bold')
    lines = [r'$R_A=E\,\mathrm{Re}\,H_3(Z_a,\overline{Z_b},n_c)$',
        r'$C_A=\frac{1}{2}E\,\mathrm{Re}\,[H_3(Z_a,n_b,\overline{Z_c})$',
        r'$\qquad\qquad+H_3(n_a,Z_b,\overline{Z_c})]$',
        r'$w_s=n_{sa}n_{sb}n_{sc}$ in both moments',
        r'$C_A=a(q)R_A$ under the stated model']
    for i, line in enumerate(lines):ax.text(0,.88-i*.17,line,va='top',fontsize=11)
    ax = fig.add_subplot(grid[1, 1]);ax.axis('off')
    ax.set_title('c  Separate normalization can retain bias', loc='left', weight='bold')
    lines = [r'$R_2=\sum_s u_s k_s,\quad C_2=a(q)\sum_s v_s k_s$',
        r'$C_2/R_2=a(q)\,\frac{\sum_s v_s k_s}{\sum_s u_s k_s}$',
        r'Kernel invariance for every admissible $k_s>0$',
        r'holds exactly when $u_s=v_s$ for all sources.',
        'A population identity, not a physical certificate.']
    for i, line in enumerate(lines):ax.text(0,.88-i*.18,line,va='top',fontsize=10.5)
    save(fig,out,'F1')
    return dict(kind='schematic and algebra; not experimental data', layout_seed=72,
                source_code_indices=[0,1,0,2], illustrative_counts=counts)


def reference(data, out):
    rows = json.loads((data/'figure2-all-views.json').read_text())['views']
    labels = ['NPC','E1','E2','F4','R1','R2','R3','SDC c','SDC r']
    fig, axs = plt.subplots(1,3,figsize=(11.7,5.5),sharey=True,
        gridspec_kw={'width_ratios':[1.1,1.35,.8]}, layout='constrained')
    y = np.arange(len(rows))
    for i,r in enumerate(rows):
        a,b = r['signed_full_gap'],r['signed_selected_gap']
        axs[0].plot([a,b],[i,i],color='#b5bdc4',lw=1.1,zorder=1)
        axs[0].scatter(a,i,color='#777777',marker='x',s=35,zorder=3)
        axs[0].scatter(b,i,color=BLUE,s=30,zorder=3)
        for j,d in enumerate(r['frequency_differences']):
            v = d['code_minus_diagonal']; interval = d['conditional_interval']
            yy = i+(.13 if j else -.13);color = (BLUE,ORANGE)[j]
            if interval is not None:
                axs[1].hlines(yy,*interval,color=color,lw=1.25)
                axs[1].scatter(v,yy,color=color,s=21,zorder=3)
            else:axs[1].scatter(0 if v is None else v,yy,color=color,marker='x',s=35)
        pct = 100*r['code_input_coverage']
        axs[2].barh(i,pct,color=BLUE,height=.57)
        axs[2].text(pct+.8,i,f'{pct:.1f}',va='center',fontsize=10)
    axs[0].set(yticks=y,yticklabels=labels,ylim=(8.8,-.8),xlabel='Fit minus diagonal\n(nm$^2$-equivalent)')
    axs[1].set(xlabel='Same-row difference with\ncentral 95% valid-draw range')
    axs[2].set(xlabel='Code-input rows\nretained (%)',xlim=(0,58))
    for a,title in zip(axs,('a  All nine views','b  Uncertainty and residual differences','c  Selection coverage')):
        a.set_title(title,loc='left',weight='bold',fontsize=11)
        a.grid(axis='y',color='#edf0f3',lw=.7);a.set_axisbelow(True)
    for a in axs[:2]:a.axvline(0,color='#777777',lw=1,ls=':')
    axs[0].legend(handles=[Line2D([],[],color='#777777',marker='x',ls='',label='Original'),
        Line2D([],[],color=BLUE,marker='o',ls='',label='Matched rows')],loc='lower right',frameon=False)
    axs[1].legend(handles=[Line2D([],[],color=c,marker='o',ls='-',label=f'q={q}')
        for c,q in ((BLUE,'0.10'),(ORANGE,'0.15'))],loc='lower left',frameon=False)
    save(fig,out,'F2')
    return dict(kind='retained data; no new fitting', views=[r['view'] for r in rows],
                point_values='data/figure2-all-views.json', omitted_views=0)


def covariance(data, out):
    d=json.loads((data/'figure4-complete.json').read_text());s=d['summary']
    ratios=[r['ratio'] for r in d['component_variance_ratios']]
    fig,axs=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    hist=axs[0].hist(ratios,bins=35,color=BLUE,edgecolor='white',linewidth=.4)
    axs[0].axvline(1,color=INK,ls='--',lw=1.1)
    axs[0].set(xlabel='Monte Carlo / exact variance',ylabel='Statistic components')
    for variant,color,offset in [('arithmetic',BLUE,-.15),('geometric',ORANGE,.15)]:
        rows=[r for r in s['summary'] if r['variant']==variant]
        axs[1].plot([r['K'] for r in rows],[r['median_predicted_variance_factor'] for r in rows],
            'o-',color=color,label=variant,lw=1.5)
        axs[2].bar(np.arange(3)+offset,[100*r['usable_covers']/r['total'] for r in rows],
            width=.30,color=color,label=variant)
    axs[1].axhline(1,color=INK,ls='--',lw=1.1)
    axs[1].set(xlabel='Source codes K',ylabel='Median first-order variance / K=4',xticks=[4,8,16])
    axs[1].legend(frameon=False,fontsize=9)
    axs[2].set(xlabel='Codes / independent groups',ylabel='Usable and covering / all calls (%)',
        ylim=(0,100),xticks=range(3),xticklabels=['4 / 64','8 / 32','16 / 16'])
    for ax,title in zip(axs,('a  Coding covariance check','b  Stipulated equal budget','c  Finite measurement support')):
        ax.set_title(title,loc='left',weight='bold',fontsize=11)
    save(fig,out,'F4')
    return dict(kind='retained summaries; no simulations rerun', component_count=len(ratios),
                histogram_counts=hist[0].astype(int).tolist(),histogram_edges=hist[1].tolist(),
                summary=s['summary'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--selection-python',default=sys.executable,
        help='Python with the recorded historical Matplotlib environment for Figure 3')
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    records={'1':construction(out),'2':reference(args.data,out),'4':covariance(args.data,out)}
    subprocess.run([args.selection_python,str(Path(__file__).with_name('selection_figure.py')),'--output',str(out)],check=True)
    records['3']=dict(kind='complete regenerated selection display',settings=81,views=9)
    files=[dict(path=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
           for p in sorted(out.iterdir())]
    (out/'manifest.json').write_text(json.dumps(dict(figures=records,files=files,
        new_fits=0,new_outcome_study=False),indent=2)+'\n')


if __name__=='__main__':main()
