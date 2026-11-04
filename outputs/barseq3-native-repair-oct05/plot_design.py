"""Scientific schematic of the frozen design; contains no outcome data."""
import ast
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from condition import HERE
ns={'np':np}
node=next(n for n in ast.parse((HERE/'recovery.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='offsets')
exec(compile(ast.Module(body=[node],type_ignores=[]),'frozen design','exec'),ns)
fig,axes=plt.subplots(1,3,figsize=(9.2,4.9),gridspec_kw={'width_ratios':[1,1,1.2]})
for ax,kind,title in zip(axes[:2],['common','color'],['Offsets grouped by round','Offsets grouped by color']):
    values=ns['offsets'](kind).reshape(10,4,2)
    ax.imshow(values[:,:,0]>0,cmap=ListedColormap(['#b75938','#287899']),vmin=0,vmax=1,aspect='equal')
    ax.set_xticks(range(4),['G','T','C','A']);ax.set_yticks(range(10),[str(i+1) for i in range(10)])
    ax.set(xlabel='Readout color',ylabel='Round',title=title)
    ax.set_xticks(np.arange(-.5,4,1),minor=True);ax.set_yticks(np.arange(-.5,10,1),minor=True)
    ax.grid(which='minor',color='white',lw=1);ax.tick_params(which='minor',bottom=False,left=False)
axes[2].axis('off')
axes[2].text(0,.98,'Identical full error distribution',fontsize=11,fontweight='bold',va='top')
axes[2].text(0,.85,'20 frames: +(4, 2) px\n20 frames: −(4, 2) px',fontsize=11,linespacing=1.8,va='top')
axes[2].text(0,.63,'Mean = (0, 0)\nRMS = √20 = 4.47 px\nCovariance = [[16, 8], [8, 4]]',fontsize=10,linespacing=1.8,va='top')
axes[2].text(0,.34,'Five calibration policies\nTwo actual decoder families\nTwo fixed external fields\n60 planned arms',fontsize=10,linespacing=1.6,va='top')
fig.suptitle('Equal error budgets, different readout structure',fontsize=14)
fig.text(.5,.04,'Frozen before perturbed decode outcomes · one BARseq3 preparation · controlled recovery of native outputs',ha='center',fontsize=9)
fig.subplots_adjust(left=.07,right=.98,bottom=.16,top=.83,wspace=.32)
for extension in ['png','pdf']:fig.savefig(HERE/('recovery-design.'+extension),dpi=180)
