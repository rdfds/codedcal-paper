"""Report complete verified comparisons; no estimator or criterion selection."""
import json
import gzip
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import HERE, TAGS, sha, write


def main():
    verification=json.loads((HERE/'verification.json').read_text());assert verification['all_pass']
    allrows=json.loads((HERE/'summary.json').read_text())['summaries']
    rows=[r for r in allrows if r['scope']=='observed' and r['method']=='coded']
    table=['| Stress | Variant | Available cases | Mean joint usable coverage | Median absolute sample-mean ratio bias | Errors >0.01 / quiet-D errors |',
           '|---|---|---:|---:|---:|---:|']
    overview={}
    for world in ('noise_shift','source_noise_shift','code_offset'):
        overview[world]={}
        for variant in ('arithmetic','geometric'):
            a=[r for r in rows if r['world']==world and r['variant']==variant]
            stat=dict(cases=sum(r['records'] for r in a),available_cases=round(sum(r['available_fraction']*r['records'] for r in a)),
                      usable_joint_coverage=float(np.mean([r['usable_joint_coverage'] for r in a])),
                      median_abs_mean_bias=float(np.median([abs(r['ratio_bias']) for r in a if r['ratio_bias'] is not None])),
                      bad_errors=sum(r['bad_physical_errors'] for r in a),quiet_bad=sum(r['bad_errors_quiet_D'] for r in a))
            overview[world][variant]=stat
            table.append(f"| {world} | {variant} | {stat['available_cases']}/{stat['cases']} | {stat['usable_joint_coverage']:.1%} | {stat['median_abs_mean_bias']:.4f} | {stat['bad_errors']} / {stat['quiet_bad']} |")
    primary=['| Method | Variant | All inherited gates /144 | Above / below nominal bootstrap failures | Bias / usability / joint-coverage failures |',
             '|---|---|---:|---:|---:|']
    for method in ('coded','known_source'):
        for variant in ('arithmetic','geometric'):
            a=[r for r in allrows if r['scope']=='observed' and r['method']==method and r['variant']==variant and r['primary']]
            failures=Counter(k for r in a for k in r['failed_gates']);directions=Counter(r['coverage_bootstrap_nominal_failure_direction'] for r in a)
            primary.append(f"| {method} | {variant} | {sum(r['qualified'] for r in a)}/144 | {directions['above_nominal']} / {directions['below_nominal']} | {failures['ratio_bias']} / {failures['usable_fraction']} / {failures['usable_coverage']} |")
    matched={}
    for world in ('noise_shift','source_noise_shift','code_offset'):
        biases={v:[] for v in ('arithmetic','geometric')};common=0;total=0
        for file in (HERE/'worlds').glob('*-'+world+'.jsonl.gz'):
            buckets={}
            with gzip.open(file,'rt') as stream:
                for line in stream:
                    row=json.loads(line)
                    if row['scope']=='observed' and row['method']=='coded':
                        buckets.setdefault(row['q'],{}).setdefault((row['rep'],row['assignment']),{})[row['variant']]=row
            for pairs in buckets.values():
                available=[pair for pair in pairs.values() if all(pair[v]['available'] for v in biases)]
                total+=len(pairs);common+=len(available)
                for variant in biases:
                    if available:biases[variant].append(abs(np.mean([pair[variant]['ratio']-pair[variant]['physical_target'] for pair in available])))
        matched[world]=dict(common_available=common,total=total,median_absolute_sample_mean_bias={v:float(np.median(a)) for v,a in biases.items()})
    write(HERE/'availability-matched.json',dict(post_hoc=True,scope='observed coded; identical available paired draws',conditions=matched))
    full=['| View | Noise | q | World | Variant | Bias | Usable | Joint coverage | Width | Correct coverage + positive-motion resolution | Gates |',
          '|---|---|---:|---|---|---:|---:|---:|---:|---:|---|']
    for r in rows:
        if not r['primary']:continue
        width='undefined' if r['median_usable_interval_width'] is None else f"{r['median_usable_interval_width']:.4f}"
        bias='undefined' if r['ratio_bias'] is None else f"{r['ratio_bias']:.4f}"
        full.append(f"| {r['tag']} | {r['model']} | {r['q']:.2f} | {r['world']} | {r['variant']} | {bias} | {r['usable_fraction']:.1%} | {r['usable_joint_coverage']:.1%} | {width} | {r['usable_covers_and_excludes_zero_fraction']:.1%} | {','.join(r['failed_gates']) or 'all pass'} |")
    text='''# Existing temporal variants: robustness, support and resolution

Completed 17 October 2026. A separately specified retained-observation comparison,
not an independent experiment. All 126 Source Data 70 worlds, 3,528 simulated
physical fields, 84,672 variant cases and 2,016 condition summaries are retained.
No new realization or candidate image was generated or acquired.

## Finding and limits

The existing geometric construction cancels a source-uniform multiplicative
kernel change under the declared common-motion law. The arithmetic construction
does not generally cancel it. On the specified uniform-noise stress, median
absolute Monte Carlo mean bias across 36 coded view/noise/frequency conditions
falls from 0.0292 to 0.0043; mean joint usable coverage rises from 72.6% to 96.9%.
This is a conditional robustness benefit, not physical instrument validation.

Its conservative joint sets are wider and often lose small-motion resolution.
For NPC Gaussian T=3, correct coverage with positive-motion resolution falls
from 25.0% to 0% at q=.10 and from 78.8% to 11.2% at q=.15. EDF10 replicate 2
retains poor support. Source-dependent kernel and code-geometry violations remain;
quiet-D failures are not repaired. No data-selected switch or safe automatic
recommendation is validated.

Arithmetic uses a pointwise 95% Fieller set. Geometric uses two 97.5% window
Fieller sets propagated through their positive product, yielding a conservative
nominal 95% joint rectangle when marginal approximations are calibrated. No
independence between windows is assumed; both require independent groups.
Concentrated contributions can defeat either approximation. More coverage alone
is not more information. Gaussian-equivalent variances are conditional on the
motion model, not total molecular-distance errors.

## All declared violations

Each stress has 864 observed-count coded cases. Coverage averages and median
absolute *sample-mean bias* are across conditions, not independent preparations.
Available ratios and errors do not imply usable confidence support.

'''+ '\n'.join(table)+'''

## Availability-matched diagnostic (post hoc)

Point availability differs between variants. Recounting identical available
paired draws retains 862/864 uniform-kernel stress cases and the same median
absolute sample-mean biases, 0.0292 versus 0.0043. Source-dependent stress
retains 858/864, with 0.1267 versus 0.0415; code geometry retains 862/864,
with 0.0038 versus 0.0036. This post-hoc fairness diagnostic is not a replacement
gate or an estimate over excluded draws. Unconditional coverage denominators
remain all 864 cases. Full records are in availability-matched.json.

## Inherited gate accounting, not a ranking of conservative intervals

The original four gates are retained literally. The nominal-in-bootstrap gate
asks whether coverage is compatible with exactly 95%. That is not an appropriate
accuracy ranking for deliberately conservative joint rectangles: 100% observed
coverage can fail it from above. Neither those failures nor all-success bootstrap
degeneracy demonstrate undercoverage. Do not advertise geometric pass counts as
evidence that it is universally worse, or its greater coverage as a universal win.
Every failure, including below-nominal coverage and bias, remains below and in CSV.

'''+ '\n'.join(primary)+'''

## Targets and independent recount

Each variant is scored against its own baseline physical functional, and its
actual-diagonal observable target is recorded separately. Under heterogeneous
source motion the geometric and arithmetic targets can differ: changing the
weighted functional is not recovering the original arithmetic target. summary.csv
retains the original-target error and target change as well. The known-source
comparator has privileged synthetic labels and is not a deployable alternative.

All 42,336 arithmetic cases and 1,008 summaries reproduce Source Data 70,
including its 115/144 coded qualification verdict. Independent frozen-formula
recounts check both geometric windows in 84,672 comparisons. The initial
verification stop is retained: raw versus normalized float64 summation differs
by 1.44e-12 relatively for an already-unusable, nearly cancelled reference.
verify_v2 uses the pre-existing arithmetic point-replay tolerance of 1e-10;
observations, confidence sets, summary numbers and scientific gates did not change.
The original verifier, failed log and precision-only amendment remain.

Reproduction (repository root, pinned analysis Python): controls.py, run.py freeze
in a new working copy before outcomes; run.py execute; run.py summarize;
verify_v2.py; report.py. The existing freeze and completed outcomes refuse silent
replacement. All source hashes and compressed per-case sets are retained.

## Complete observed-count coded primary table

'''+ '\n'.join(full)+'\n'
    (HERE/'RESULTS.md').write_text(text)
    supplementary='# Supplementary Table 4: existing temporal variants\n\nSource Data 71 reuses all 126 worlds and 3,528 simulated fields; 84,672 cases\nare independently recounted. This is conditional computation, not physical\nvalidation. Original Source Data 70 results and failed gates remain unchanged.\n\n'+text.split('## Finding and limits',1)[1]
    (HERE/'table4.md').write_text(supplementary)
    (HERE.parents[1]/'manuscript/temporal-variants.md').write_text(supplementary)
    write(HERE/'overview.json',overview)
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'pdf.fonttype':42,
                         'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(10,7.1),gridspec_kw={'hspace':.47,'wspace':.38})
    index={(r['tag'],r['model'],r['world'],r['q'],r['variant']):r for r in rows}
    pairs=lambda world:[(index[tag,model,world,q,'arithmetic'],index[tag,model,world,q,'geometric'])
                       for tag in TAGS for model in ('gaussian','residual') for q in (.1,.15)]
    ax=axes[0,0]
    for world,color,marker,label in [('noise_shift','#2777b5','o','Uniform kernel'),
            ('source_noise_shift','#d95c36','^','Source-dependent kernel'),('code_offset','#777777','x','Code geometry')]:
        xy=[(abs(a['ratio_bias']),abs(g['ratio_bias'])) for a,g in pairs(world) if a['ratio_bias'] is not None and g['ratio_bias'] is not None]
        ax.scatter(*np.array(xy).T,s=19,c=color,marker=marker,label=label,alpha=.8)
    ax.set_xscale('symlog',linthresh=.01);ax.set_yscale('symlog',linthresh=.01)
    ax.plot([0,1],[0,1],c='#888',ls=':',lw=.8)
    ax.set(xlabel='Arithmetic: absolute sample-mean ratio bias',ylabel='Geometric: absolute sample-mean ratio bias',title='a  Robustness depends on the violation')
    ax.legend(frameon=False,fontsize=7)
    ax=axes[0,1]
    for world in ('zero','homogeneous','group_common','heterogeneous'):
        for a,g in pairs(world):
            color='#d95c36' if a['tag']=='edf10rep2' else '#2777b5'
            ax.scatter(a['usable_joint_coverage'],g['usable_joint_coverage'],s=16,c=color,alpha=.6)
    ax.plot([0,1],[0,1],c='#888',ls=':',lw=.8);ax.axhline(.90,c='#ccc',lw=.7);ax.axvline(.90,c='#ccc',lw=.7)
    ax.set(xlabel='Arithmetic joint usable coverage',ylabel='Geometric joint usable coverage',title='b  Coverage retains unsupported cases',xlim=(-.03,1.04),ylim=(-.03,1.04))
    ax.text(.04,.96,'144 primary conditions; orange: EDF10 replicate 2\nConservative rectangles are not exact-95% tests',transform=ax.transAxes,va='top',fontsize=7)
    ax=axes[1,0]
    for a,g in pairs('homogeneous'):
        ax.scatter(a['usable_covers_and_excludes_zero_fraction'],g['usable_covers_and_excludes_zero_fraction'],
                   c='#2777b5' if a['model']=='gaussian' else '#d95c36',marker='o' if a['q']==.1 else '^',s=23)
    ax.plot([0,1],[0,1],c='#888',ls=':',lw=.8)
    ax.set(xlabel='Arithmetic: correct coverage + positive-motion resolution',ylabel='Geometric: same endpoint',title='c  Wider sets often lose small-motion resolution',xlim=(-.03,1.04),ylim=(-.03,1.04))
    ax.text(.04,.84,'T=3; 36 conditions, not biological replicates\nBlue/orange: Gaussian/residual\nCircle/triangle: q=.10/.15',transform=ax.transAxes,va='top',fontsize=7)
    ax=axes[1,1];worlds=['noise_shift','source_noise_shift','code_offset'];x=np.arange(3)
    for shift,variant,color in [(-.18,'arithmetic','#2777b5'),(.18,'geometric','#d95c36')]:
        values=[overview[w][variant]['quiet_bad'] for w in worlds]
        ax.bar(x+shift,values,width=.34,color=color,label=variant)
        for xx,y in zip(x+shift,values):ax.text(xx,y+6,str(y),ha='center',fontsize=8)
    ax.set(xticks=x,xticklabels=['Uniform\nkernel','Source-dependent\nkernel','Code\ngeometry'],ylabel='Errors >0.01 with quiet D /864 cases',title='d  Quiet D remains insufficient',ylim=(0,570))
    ax.legend(frameon=False)
    for extension in ('png','pdf','svg'):fig.savefig(HERE/f'variants.{extension}',dpi=200,bbox_inches='tight')
    svg=HERE/'variants.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)
    files=['RESULTS.md','table4.md','availability-matched.json','overview.json','variants.png','variants.pdf','variants.svg','report.py','summary.json','summary.csv','verification.json','verification-amendment.json','verify_v2.py']
    write(HERE/'report-manifest.json',{name:sha(HERE/name) for name in files})
    print('Reported complete verified comparison',flush=True)


if __name__=='__main__':main()
