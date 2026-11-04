"""Show all fixed external conditions and alternatives, including the zero control."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from condition import HERE,read

def main():
    s=read(HERE/'recovery/summary.json');assert s['full_cohort_comparison_supported']
    policies=['none','common_iterative','population_template','barcode_model','oracle']
    labels=['No added correction','Strong round-only','Per-frame template','Author model','Known-offset oracle']
    colors=['#888888','#ac6441','#287899','#745399','#373a3c']
    fig,axes=plt.subplots(2,2,figsize=(10.4,7.1),sharex=True,sharey=True)
    for field in [0,1]:
        for col,decoder in enumerate(['correlation','regression']):
            ax=axes[field,col];r=s['fields'][str(field)]['decoders'][decoder]
            for j,(p,label,c) in enumerate(zip(policies,labels,colors)):
                y=[case['scores'][p]['5'] for case in r['cases']]
                x=np.arange(3)+(j-2)*.045
                ax.plot(x,y,marker=['o','s','^','D','v'][j],color=c,lw=1.2,ms=5,label=label,alpha=.9)
            chosen=[case['scores'][case['recommend']]['5'] for case in r['cases']]
            ax.scatter(np.arange(3),chosen,s=90,facecolors='none',edgecolors='#111111',linewidths=1.2,zorder=8,label='Frozen ≥2 px rule')
            ax.set_title(f"{'Development' if field==0 else 'Spatial holdout'} · {decoder}",fontsize=11)
            ax.set(ylim=(-.025,1.08),xticks=[0,1,2],xticklabels=['Zero','By round','By color'])
            ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
            ax.text(.025,.025,f"Harmful recommendations: {r['harmful_recommendations']}/3",transform=ax.transAxes,fontsize=8)
    for ax in axes[:,0]:ax.set_ylabel('Native-output identity consistency F1')
    for ax in axes[1]:ax.set_xlabel('Error grouping (round/color distribution identical)')
    handles,legend=axes[0,0].get_legend_handles_labels();fig.legend(handles,legend,loc='lower center',ncol=3,fontsize=9,bbox_to_anchor=(.5,.055))
    fig.suptitle('Matched full error distributions across calibration and decoder families',fontsize=13)
    fig.text(.5,.02,'60/60 fixed arms · one external BARseq3 preparation · algorithm reference, not independent RNA truth · integer corrections',ha='center',fontsize=8)
    fig.subplots_adjust(left=.08,right=.98,top=.89,bottom=.23,hspace=.28,wspace=.14)
    for ext in ['png','pdf']:fig.savefig(HERE/'recovery'/('comparison.'+ext),dpi=180)

if __name__=='__main__':main()
