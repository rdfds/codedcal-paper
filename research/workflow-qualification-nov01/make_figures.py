"""Descriptive scientific figures only, from stopped/adjudicated evidence."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
H=Path(__file__).resolve().parent;d=json.loads((H/'adjudicated-summary.json').read_text())
methods=['conventional','mixed_conventional','balanced_random','readout_aware'];names=['Round/color','Mixed conventional','Balanced random','Readout-aware']
rows=[r for r in d['rows'] if r['source_qualified'] is True]
fig=plt.figure(figsize=(12,max(8,.31*len(rows)+3)));gs=fig.add_gridspec(2,2,height_ratios=[max(5,.3*len(rows)),2.4],width_ratios=[1.5,1],hspace=.47,wspace=.35)
ax=fig.add_subplot(gs[0,0]);value={'unsupported':0,'accept':1,'reject':2};matrix=np.array([[value[r[k]['decision']] for k in ['truth']+methods] for r in rows])
ax.imshow(matrix,aspect='auto',cmap=ListedColormap(['#B9BDC4','#357C61','#BC4B40']),vmin=0,vmax=2)
ax.set_xticks(range(5),['Held-out\nbank','Round/\ncolor','Mixed','Random','Readout-\naware'],fontsize=9)
ax.set_yticks(range(len(rows)),[r['name']+' / '+r['decoder']+' / '+r['policy']+(' *' if not r['native_workflow_qualified'] else '') for r in rows],fontsize=8)
ax.set_title('Qualified-source workflows; finite digital challenges',fontsize=11)
ax.legend(handles=[Patch(facecolor=c,label=t) for c,t in [('#357C61','Adequate / accept'),('#BC4B40','Failed / reject'),('#B9BDC4','Unsupported')]],loc='lower center',bbox_to_anchor=(.5,-.22),ncol=3,fontsize=8,frameon=False)
bx=fig.add_subplot(gs[0,1]);palette=['#5B5B5B','#8B6BB1','#CC9850','#2865A8'];plotted=0
for method,name,color in zip(methods,names,palette):
 c=d['confusion'][method]
 if c['false_alert_rate'] is None or c['failed_detection_rate'] is None:continue
 bx.scatter(c['false_alert_rate'],c['failed_detection_rate'],s=90,color=color,label=name,alpha=.8);plotted+=1
bx.set(xlim=(-.04,1.04),ylim=(-.04,1.04),xlabel='False-alert rate on adequate workflows',ylabel='Detection rate on failed workflows')
bx.grid(alpha=.2);bx.set_title(f"Adequate n={d['adequate_workflows']}; failed n={d['failed_workflows']}",fontsize=11)
if plotted:bx.legend(loc='lower right',fontsize=8,frameon=False)
else:bx.text(.5,.5,'Rates not jointly estimable\nDo not interpret reject-all as utility',ha='center',va='center',transform=bx.transAxes,fontsize=10)
cx=fig.add_subplot(gs[1,:]);policies=list(d['per_policy_measured_union_runtime']);cal=[];dec=[];other=[]
for p in policies:
 r=d['per_policy_measured_union_runtime'][p];cal.append(r['calibration_seconds']);dec.append(r['decoder_seconds']);other.append(max(0,r['wall_seconds']-r['calibration_seconds']-r['decoder_seconds']))
x=np.arange(len(policies));cx.bar(x,cal,label='Calibration',color='#2865A8');cx.bar(x,dec,bottom=cal,label='Decoding',color='#4E927E');cx.bar(x,other,bottom=np.array(cal)+dec,label='Other adapter work / exact cache lookup',color='#B9BDC4')
cx.set_xticks(x,policies);cx.set_ylabel('Measured seconds');cx.set_title('Measured union of queries; shared cache reuse is charged once',fontsize=11);cx.legend(fontsize=8,frameon=False)
status='Completed development matrix' if d['complete'] else 'Stopped / incomplete development matrix'
fig.suptitle(status+' — four challenge pairs per strategy',fontsize=14,y=.995)
fig.text(.02,.015,'Descriptive workflow cases, not independent biological replicates. Eight source-gate and eight compute exclusions are retained in the full report.\n* Native workflow gate failed; held-out truth unsupported. Invocation matched; identical wall-clock cost is not asserted. No independent comparison: utility gate failed.',fontsize=8)
fig.subplots_adjust(left=.32,right=.98,top=.94,bottom=.13)
for suffix in ['png','svg']:fig.savefig(H/('qualification-development.'+suffix),dpi=180,bbox_inches='tight')
plt.close(fig)
