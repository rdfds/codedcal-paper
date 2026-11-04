"""Descriptive scientific display; all applicable worlds retained."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('--prefix',required=True);p.add_argument('--report',required=True);args=p.parse_args()
prefix=Path(args.prefix);report=Path(args.report)
s=json.loads((report/'summary.json').read_text())
ratio=[]
for f in prefix.glob('*.json'):
    d=json.loads(f.read_text())
    if d.get('theorem_applicable'):
        v=d['verification'];ratio.extend(np.diag(v['mc_second_moment'])/np.diag(v['predicted_covariance']))
fig,axs=plt.subplots(1,3,figsize=(12,3.4),layout='constrained')
axs[0].hist(ratio,bins=35,color='#4477aa',edgecolor='white',linewidth=.3)
axs[0].axvline(1,color='black',ls='--',lw=1)
axs[0].set(xlabel='Monte Carlo / exact conditional variance',ylabel='Statistic components',title='A  Conditional covariance check')
for variant,color,offset in [('arithmetic','#4477aa',-.12),('geometric','#cc6677',.12)]:
    rows=[r for r in s['summary'] if r['variant']==variant]
    x=np.array([r['K'] for r in rows]);v=[r['median_predicted_variance_factor'] for r in rows]
    axs[1].plot(x,v,'o-',color=color,label=variant)
    axs[2].bar(np.arange(3)+offset,[100*r['usable_covers']/r['total'] for r in rows],width=.24,color=color,label=variant)
axs[1].axhline(1,color='black',ls='--',lw=1)
axs[1].set(xlabel='Permanent source codes K',ylabel='Median first-order variance / four codes',title='B  Equal stipulated budget')
axs[1].set_xticks([4,8,16]);axs[1].legend(frameon=False,fontsize=8)
axs[2].set(xlabel='Codes / independent groups',ylabel='Usable and covering calls (%)',title='C  Finite measurement support',ylim=(0,100))
axs[2].set_xticks(range(3),['4 / 64','8 / 32','16 / 16'])
for ax in axs:
    ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8)
fig.savefig(report/'covariance-design.png',dpi=180)
fig.savefig(report/'covariance-design.svg')
(report/'FIGURE_CAPTION.md').write_text('''# Conditional covariance and finite-budget design

A: All 864 component variance comparisons from 108 exogenous worlds. Each
Monte Carlo estimate uses 512 independent assignments to fixed source
observations; the histogram retains every discrepancy. The dashed line is
agreement, not a confidence bound. B: Across all 216 supported world/frequency
conditions per variant, median predicted first-order variance relative to
four codes. C: Descriptive usable-and-covering frequency, including every
planned call; no biological-replication or nominal-coverage claim. Error bars
are omitted because pooling overlapping worlds does not yield independent
biological uncertainty. All panels use exposed observations. B and C assume
budget B=K*G=256 with counts/kernels/photons unchanged and observations drawn
uniformly from a fixed sixteen-group prefix. The physical and code variance
terms use that same law. This is conditional computational verification, not
physical acquisition validation. Lower variance does not ensure more usable
measurements. See full tables, failures and scope in RESULTS.md.
''')
