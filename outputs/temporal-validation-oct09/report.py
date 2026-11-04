"""Complete-cohort interpretation and an exportable scientific summary."""
import json
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import validation as V
from common import save, sha


def pct(x):return f'{100*x:.1f}%'


def main():
    s=json.loads((V.P/'summary.json').read_text());assert s['complete']
    verify=json.loads((V.P/'verification.json').read_text());assert verify['all_pass']
    primary=[r for r in s['records'] if r['primary'] and r['scope']=='observed']
    coded=[r for r in primary if r['method']=='coded'];known=[r for r in primary if r['method']=='known_source']
    stress=[r for r in s['records'] if not r['primary'] and r['scope']=='observed' and r['method']=='coded']
    n_code=sum(r['qualified'] for r in coded);n_known=sum(r['qualified'] for r in known)
    both=sum(c['qualified'] and k['qualified'] for c in coded for k in known if all(c[x]==k[x] for x in ('tag','model','world','q')))
    fail=Counter(g for r in coded for g in r['failed_gates'])
    direction=Counter(r['coverage_bootstrap_nominal_failure_direction'] for r in coded)
    lines=['# Preparation-derived temporal validation', '',
           'Completed 10 October 2026; settings frozen 9 October. Complete, frozen **conditional simulation** from eight public acquisitions / nine raw-corrected views. No new acquisition, physical randomized-readout intervention or independent physical source truth.', '',
           f"All 126 prescribed worlds completed: 72 primary and 54 stress worlds, {s['physical_fields']:,} simulated physical-field realizations, {s['coded_fields']:,} coded observations, and {s['coded_group_draws']:,} group/code draws. Paired codes share a physical field; they are clustered in the bootstrap. All {verify['case_recounts']:,} pointwise cases, 126 analytic truth tables, nine source catalogs and retained artifact hashes passed recount. This audit pass is not a scientific qualification pass.", '',
           '## Main finding', '',
           f'At the observed group counts, **{n_code}/{len(coded)} coded primary world/frequency conditions pass all frozen gates**, compared with **{n_known}/{len(known)}** for the source-label-informed Fourier comparator; **{both}** pass for both. Gates require absolute mean ratio bias ≤0.01, ≥90% finite positive confidence calls, ≥90% joint usable physical-target coverage, and a whole-field bootstrap interval including nominal 95%. These are pointwise study criteria, not established field standards or multiplicity-adjusted certification.', '',
           'Fieller sets use a group-mean normal approximation with a Student critical value; they are nominal confidence sets, not exact finite-sample coverage guarantees. Highly concentrated count/kernel weights can leave far fewer effectively contributing groups than the raw group count. Coverage failures in either method remain failures. More gate passes do not imply better precision: broader intervals may cover more often.', '',
           'The group-common world shares motion within one synthetic pore/group and draws it independently across groups. A motion component shared by an entire acquisition is outside this independent-group uncertainty test. The persistent-window law is imposed, not fitted from real frame lags.', '',
           'The source-informed comparator has privileged synthetic labels and is an information benchmark. It does not demonstrate an implementable stronger pipeline, and this experiment does not establish a benefit over Shaw–Veatch, which measures a different reference-relative lag quantity. No optimizer boundary or quiet diagnostic certifies a measurement.', '',
           f"Coded primary gate failures (overlap allowed): {dict(fail)}. Failed criteria remain failures; uncertainty intervals and extreme ratio tails are retained rather than discarded.", '',
           f'Nominal bootstrap failures by direction: {dict(direction)}. A percentile bootstrap at all-success coverage is degenerate; above-nominal failures are not evidence of undercoverage. A separately labeled post-hoc Chernoff interval for the mean bounded physical-field coverage is included for every condition, preserving arbitrary dependence between its two code assignments. It is a finite-sample sensitivity and does not replace any frozen gate.', '',
           'EDF10 replicate 2 is the strongest failure case: the coded finite-ratio sample mean has maximum absolute bias 0.407, minimum joint usable coverage 12.5%, and a ratio tail reaching 23.7. The known-source comparator reaches only 5% joint coverage in its worst heterogeneous residual condition. Median agreement or a correct extraction identity cannot erase these failures. Finite-ratio mean/RMSE diagnostics are sample summaries; finite population ratio moments are not assumed.', '',
           '## All views and noise laws', '',
           '| View | Noise | Coded qualified /8 | Known-source /8 | Coded usable range | Coded joint coverage range | Worst absolute ratio bias |',
           '|---|---|---:|---:|---:|---:|---:|']
    for tag in V.TAGS:
        for model in V.MODELS:
            c=[r for r in coded if r['tag']==tag and r['model']==model];k=[r for r in known if r['tag']==tag and r['model']==model]
            lines.append(f"| {tag} | {model} | {sum(r['qualified'] for r in c)}/8 | {sum(r['qualified'] for r in k)}/8 | {pct(min(r['usable_fraction'] for r in c))}–{pct(max(r['usable_fraction'] for r in c))} | {pct(min(r['usable_joint_coverage'] for r in c))}–{pct(max(r['usable_joint_coverage'] for r in c))} | {max((abs(r['ratio_bias']) for r in c if r['ratio_bias'] is not None),default=float('nan')):.4f} |")
    lines += ['', 'These are repeated simulations conditioned on exposed preparations, not eighteen independent acquisitions. SDC raw/corrected share one acquisition; the Figure4 RTX alias is counted once. The 512-group sensitivity uses the same draws and is not independent confirmation.', '',
              '## Measurement resolution', '',
              'Finite positive Fieller uncertainty does not imply useful precision for small motion. The table reports pointwise exclusion of zero motion (attenuation ratio=1), plus correct target coverage; these descriptive endpoints were added after early outcomes and do not change frozen qualification.', '',
              '| View | Noise | Homogeneous T=3: correct coverage and zero exclusion q=.10 / .15 | Median interval width q=.10 / .15 | Known-source correct zero exclusion q=.10 / .15 |',
              '|---|---|---:|---:|---:|']
    for tag in V.TAGS:
        for model in V.MODELS:
            c=sorted([r for r in coded if r['tag']==tag and r['model']==model and r['world']=='homogeneous'],key=lambda r:r['q'])
            k=sorted([r for r in known if r['tag']==tag and r['model']==model and r['world']=='homogeneous'],key=lambda r:r['q'])
            lines.append(f"| {tag} | {model} | {' / '.join(pct(r['usable_covers_and_excludes_zero_fraction']) for r in c)} | {' / '.join((format(r['median_usable_interval_width'],'.4f') if r['median_usable_interval_width'] is not None else 'unavailable') for r in c)} | {' / '.join(pct(r['usable_covers_and_excludes_zero_fraction']) for r in k)} |")
    lines += ['', '## Physical interpretation and diagnostic misses', '',
              'Violations test the frozen arithmetic ratio. The existing split geometric construction can cancel source-uniform multiplicative kernel changes under its own assumptions; these results must not be described as failures of every temporal variant. No new strongest-variant comparison was frozen here.', '',
              'The physical attenuation target and observable source-diagonal ratio are recorded separately in every world. Window-B noise changes can preserve randomized source extraction while changing the ratio’s interpretation. Static code-dependent positions can break extraction itself. Heterogeneous source motion targets the exact count/kernel-weighted attenuation, not arithmetic motion variance.', '',
              '| Violation | Coded cases | Errors >0.01 | Such errors with quiet D | Invalid D | Alerts / valid D | Observable-to-physical target difference range |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for world in V.STRESS:
        r=[a for a in stress if a['world']==world];cases=sum(a['records'] for a in r);bad=sum(a['bad_physical_errors'] for a in r);miss=sum(a['bad_errors_quiet_D'] for a in r)
        inv=sum(a['D_invalid'] for a in r);alert=sum(a['D_alerts'] for a in r);valid=cases-inv
        diff=[abs(a['observable_target']-a['physical_target']) for a in r]
        lines.append(f"| {world} | {cases} | {bad} | {miss} | {inv} | {pct(alert/valid) if valid else 'unavailable'} | {min(diff):.4f}–{max(diff):.4f} |")
    d_null=[r for r in coded if r['common_motion_D_null_declared']]
    false=sum(r['D_alerts'] for r in d_null);invalid=sum(r['D_invalid'] for r in d_null);n=sum(r['records'] for r in d_null)
    lines += ['', f'Across common-motion null coded cases, D produces {false}/{n-invalid} valid-test alerts and {invalid} invalid tests. These pooled counts are descriptive; acquisitions, frequencies and paired codes are not independent biological replicates. Heterogeneous source motion can produce different weighted window targets despite stable kernels; its alerts are not false alarms under a common-motion null. Both analytic window D targets and per-condition rates are in `summary.json`.', '',
              '## Sampled reference coverage', '',
              'Only homogeneous T=3 uses the frozen first-eight-fields, first-sixteen-groups clustering sample. Historical radii/QC remain fixed. Mixed assignments earn no known-source coverage; fragments receive their captured count-triple fraction. This diagnostic tests inferred synthetic identities, not independent biological labeling.', '',
              '| View | Noise | Passing clusters | Passing impurity | Localization coverage | Pure stationary kernel mass retained |',
              '|---|---|---:|---:|---:|---:|']
    for w in s['worlds']:
        c=w['clustering']
        if c is not None:
            lines.append(f"| {w['tag']} | {w['model']} | {c['passing_clusters']} | {pct(c['impure_passing_fraction']) if c['impure_passing_fraction'] is not None else 'unavailable'} | {pct(c['localization_coverage'])} | {pct(c['pure_kernel_mass_coverage'])} |")
    lines += ['', '## What this changes', '',
              'This closes the missing preparation-derived finite-observation test under its declared laws. It identifies exactly which conditions support useful uncertainty and which do not, while retaining all source/QC failures. It does **not** establish that the acquired source assignments or persistent-window law describe the instrument, validate an unperformed randomized MERFISH intervention, or resolve the independent molecular-benefit gap. Passing conditions support a conditional measurement claim; failed conditions limit it.', '',
              'No favorable subset becomes the headline for unrestricted accuracy. The complete table contains every preparation, noise model, frequency, scope and comparator; the summary plot shows all primary observed-count conditions. All original failed S3/ISS/molecular results and density qualifications remain unchanged.', '',
              '## Retained evidence', '',
              '`PROTOCOL.md` and `freeze.json` predate outcomes. `controls.json` records analytic extraction, kernel, persistent-motion and Fieller checks; `IMPLEMENTATION.md` records the pre-outcome type-label repair. `catalog-*.npz` holds all retained inferred sources; hashed empirical inputs hold event residuals. Each `worlds/*/truth.npz` holds physical/observable expected moments; every `field-*.npz` retains group observations and template IDs, with seeds/uncertainty/diagnostics/QC in its JSON. `summary.csv` is the full condition table; `verification.json` and `artifact-manifest.json` provide recount and hashes. No outcome-driven additional simulation or threshold selection was performed.', '']
    (V.P/'RESULTS.md').write_text('\n'.join(lines))
    info=json.loads((V.P/'population-information.json').read_text())['records']
    lines += ['', '## What population is measured', '',
              'All retained candidate clusters enter generation, but sources lacking three-panel support have zero diagonal weight. The observable is count/kernel weighted, not the mean budget of every candidate or raw localization. Noise fallback is inferred from the preparation; it is not measured physical truth. Effective weight counts below quantify concentration, not independent biological replication.', '',
              '| View | Noise | Candidates | Triple support | Effective sources q=.10 / .15 | Effective groups q=.10 / .15 | Heterogeneous apparent T q=.10 / .15 (source mean=3) |',
              '|---|---|---:|---:|---:|---:|---:|']
    for tag in V.TAGS:
        for model in V.MODELS:
            r=sorted([a for a in info if a['tag']==tag and a['model']==model],key=lambda a:a['q'])
            lines.append(f"| {tag} | {model} | {r[0]['inferred_sources']} | {r[0]['triple_supported_sources']} | {' / '.join(format(a['effective_weighted_sources'],'.1f') for a in r)} | {' / '.join(format(a['effective_weighted_groups'],'.1f') for a in r)} | {' / '.join(format(a['heterogeneous_finite_q_apparent_T'],'.2f') for a in r)} |")
    (V.P/'RESULTS.md').write_text('\n'.join(lines))
    rowkeys=list(__import__('itertools').product(V.TAGS,V.MODELS));columns=list(__import__('itertools').product(V.PRIMARY,V.Q))
    fig,axes=plt.subplots(1,2,figsize=(13,8),sharey=True)
    for ax,method,title in zip(axes,['coded','known_source'],['Aggregate random-code estimator','Known-source Fourier comparator']):
        a=np.zeros((len(rowkeys),len(columns)));coverage=np.zeros_like(a)
        for i,(tag,model) in enumerate(rowkeys):
            for j,(world,q) in enumerate(columns):
                r=next(r for r in primary if r['tag']==tag and r['model']==model and r['world']==world and r['q']==float(q) and r['method']==method)
                a[i,j]=1 if r['qualified'] else (2 if r['failed_gates']==['nominal_coverage_in_bootstrap'] and r['coverage_bootstrap_nominal_failure_direction']=='above_nominal' else 0);coverage[i,j]=100*r['usable_joint_coverage']
        ax.imshow(a,cmap=ListedColormap(['#f2d2c9','#d0e8da','#f1e0b4']),vmin=0,vmax=2,aspect='auto')
        for i in range(len(rowkeys)):
            for j in range(len(columns)):ax.text(j,i,f'{coverage[i,j]:.0f}',ha='center',va='center',fontsize=7,color='#202020')
        ax.set_title(title,fontsize=11)
        ax.set_xticks(range(8),[f'{w.replace("group_common","common").replace("heterogeneous","hetero")}\nq={q:.2f}' for w,q in columns],fontsize=7,rotation=30,ha='right')
        ax.set_yticks(range(18),[f'{tag} / {model[:3]}' for tag,model in rowkeys],fontsize=8)
        ax.set_xticks(np.arange(-.5,8,1),minor=True);ax.set_yticks(np.arange(-.5,18,1),minor=True)
        ax.grid(which='minor',color='white',linewidth=1);ax.tick_params(which='minor',length=0)
    fig.suptitle('Conditional temporal measurement validation: all preparations and primary worlds',fontsize=12)
    fig.text(.5,.012,'Numbers: joint usable nominal 95% coverage (%). Green: all frozen gates pass. Red: a bias, usability or coverage gate fails.\nAmber: only the above-nominal bootstrap gate fails; the frozen verdict remains failed.\nObserved group counts; 40 simulated fields (two coded assignments each). Inferred synthetic labels; no physical validation.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.15,right=.98,top=.93,bottom=.14,wspace=.06)
    fig.savefig(V.P/'validation.png',dpi=180);fig.savefig(V.P/'validation.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    save(V.P/'report-manifest.json',{name:sha(V.P/name) for name in ('RESULTS.md','validation.png','validation.pdf','report.py','integrate.py','UNCERTAINTY.md')})
    print(f'coded {n_code}/{len(coded)}; known source {n_known}/{len(known)}; joint {both}',flush=True)

if __name__=='__main__':main()
