"""Research summary from retained results only; does not read acquisition arrays."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent
cal=json.loads((BASE/'calibration-feasibility/result.json').read_text())
sim=json.loads((BASE/'simulation-results.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
phi=np.linspace(0,np.pi/2,200)
for a,color in [(1.,'#1565c0'),(.8,'#bb4b00')]:
    ax[0].plot(phi,a*np.cos(phi),color=color,label=f'Actual |m|={a:g}')
    ax[0].axhline(a,color=color,ls=':',lw=1)
ax[0].axhline(.96,color='.5',ls='--',lw=1)
ax[0].set(xlabel='Source-specific phase change (radians)',ylabel='Population coherence G',
          ylim=(0,1.05),title='A  Conservative bound, not equality')
ax[0].legend(frameon=False,fontsize=9)
names=['Original','Prior cubic','New quartic'];code=[];scatter=[]
for key in ['original_groups','prior_cubic_calibration_groups','new_quartic_calibration_groups']:
    m=cal['configurations'][key]['methods'][0]
    code.append(m['conditional_code_witness_variance']);scatter.append(m['scatter_proxy'])
x=np.arange(3)
ax[1].bar(x-.18,code,width=.36,label='Conditional code variance',color='#1565c0')
ax[1].bar(x+.18,scatter,width=.36,label='Calibration scatter proxy',color='#8099a6')
ax[1].set(yscale='log',xticks=x,xticklabels=names,ylabel='Variance in common scaled units',
          title='B  Calibration-only design improvement')
ax[1].legend(frameon=False,fontsize=8,loc='upper right')
for i,(key,label) in enumerate([('ordinary_H4','H4'),('source_resolved_oracle','Resolved oracle'),('collision_aware','Known collision flag')]):
    rows=[next(r for r in sim['results'] if r['scenario']==name) for name in ['stable','large_synthetic_stable']]
    rates=[r['methods'][key]['e']['rate'] for r in rows]
    ax[2].bar(np.arange(2)+(i-1)*.24,rates,width=.24,label=label)
ax[2].set(xticks=[0,1],xticklabels=['139 groups\n1500 repetitions','4096 groups\n200 repetitions'],
          ylim=(0,1.15),ylabel='Finite mixture test call fraction',
          title='C  Favorable synthetic control, g=.96')
ax[2].legend(frameon=False,fontsize=8,loc='upper left')
fig.savefig(BASE/'summary.png',dpi=160)
fig.savefig(BASE/'summary.svg')
