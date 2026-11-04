"""Complete selection-sensitivity display; all 81 retained settings."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
O=Path(__file__).resolve().parents[1]/'outputs'
TAGS=['npc','edf10rep1','edf10rep2','fig4untreated','rtx_rep1','rtx_rep2','rtx_rep3','sdc_corrected','sdc_raw']
LABELS=['NPC','EDF10 1','EDF10 2','Untreated','RTX 1','RTX 2','RTX 3','SDC corrected','SDC raw']
plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
    'legend.fontsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,
    'svg.fonttype':'none','savefig.bbox':'tight','axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e5e7eb','grid.linewidth':.6,
    'axes.axisbelow':True,'figure.dpi':110})

def read(p):return json.loads((O/p).read_text())

def title(ax,l,s):ax.set_title(f'{l}  {s}',loc='left',fontweight='bold',pad=10)

def selection(save):
    results=[read(f'reconciliation-sensitivity-oct23/result-{t}.json') for t in TAGS]
    gain=np.array([[a['full_to_diagonal_fit_gap']-a['selected_to_diagonal_fit_gap'] for a in r['arms']] for r in results])
    cov=np.array([[100*a['localization_coverage'] for a in r['arms']] for r in results])
    fig=plt.figure(figsize=(12,9));gs=fig.add_gridspec(2,2,height_ratios=[2,1],hspace=.52,wspace=.38)
    names=[a['name'].replace('-','\n') for a in results[0]['arms']]
    extent=max(abs(gain).max(),1e-6)
    for k,(values,cmap,norm,heading) in enumerate([(gain,'RdBu',TwoSlopeNorm(0,-extent,extent),'Matching can narrow or widen a gap'),(cov,'viridis',None,'Reference coverage changes with selection')]):
        ax=fig.add_subplot(gs[0,k]);im=ax.imshow(values,cmap=cmap,norm=norm,**({} if norm else {'vmin':0,'vmax':100}))
        ax.grid(False)
        ax.set(xticks=range(9),xticklabels=names,yticks=range(9),yticklabels=LABELS)
        ax.tick_params(labelsize=6.5);title(ax,'ab'[k],heading)
        for y in range(9):
            for x in range(9):
                a=results[y]['arms'][x]
                supported=all(a[z]['uncertainty']['all_frequencies_usable'] for z in ('selected_code','diagonal'))
                s=f'{values[y,x]:.2f}' if k==0 else f'{values[y,x]:.1f}'
                if not supported:s+='*'
                ax.text(x,y,s,ha='center',va='center',fontsize=6,color='white' if (abs(values[y,x])>extent*.55 if k==0 else values[y,x]<45) else 'black')
        ax.axvline(2.5,color='white',lw=1);ax.axvline(5.5,color='white',lw=1)
        ax.set_xlabel('Radius factor × eligibility (original: r1-original)',fontsize=7)
        fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label='Full gap − matched gap (nm²-equiv.)' if k==0 else 'Code-input rows retained (%)')
    ax=fig.add_subplot(gs[1,0]);d=read('reconciliation-sensitivity-oct23/singleton-diagnostic.json')['totals']
    vals=np.array([[68,13,0],[d['narrowed'],d['worsened'],d['unavailable']]])
    left=np.zeros(2)
    for i,(label,c) in enumerate([('Narrower','#2777b5'),('Wider','#d95c36'),('Unsupported','#999999')]):
        ax.barh([1,0],vals[:,i],left=left,color=c,label=label,height=.5)
        for y,x,w in zip([1,0],left,vals[:,i]):
            if w:ax.text(x+w/2,y,str(w),ha='center',va='center',fontsize=8,color='white' if i<2 else 'black')
        left+=vals[:,i]
    ax.set(yticks=[1,0],yticklabels=['All groups','Multi-source only'],xlim=(0,81),xlabel='Settings (complete 81-setting grid)')
    title(ax,'c','Exclude guaranteed singleton agreement')
    ax.legend(frameon=False,fontsize=7,ncol=3,loc='upper center',bbox_to_anchor=(.5,-.25))
    ax=fig.add_subplot(gs[1,1]);ax.axis('off');title(ax,'d','Agreement has a limited interpretation')
    for i,s in enumerate(['One inferred source per coarse group:', 'code and diagonal terms coincide algebraically.', '', '78/81 settings: individual Fieller sets usable.', '3/81: twelve all-real or disjoint sets retained.', '', 'Selected identities remain unvalidated;', 'zero-spanning intervals do not show equivalence.']):
        ax.text(0,.94-i*.12,s,va='top',fontsize=8,fontweight='bold' if i==0 else 'normal')
    fig.suptitle('All eight acquisitions / nine views; post-hoc conditional selection sensitivity',fontsize=12,y=.98)
    fig.subplots_adjust(top=.91,bottom=.09,left=.11,right=.95)
    save(fig,'F3')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    def save(fig,name):
        for ext in ('png','pdf','svg'):
            f=a.output/f'{name}.{ext}'
            if f.exists():raise ValueError('Refusing existing figure: '+str(f))
            fig.savefig(f,dpi=220)
            if ext=='svg':f.write_text('\n'.join(line.rstrip() for line in f.read_text().splitlines())+'\n')
        plt.close(fig)
    selection(save)
