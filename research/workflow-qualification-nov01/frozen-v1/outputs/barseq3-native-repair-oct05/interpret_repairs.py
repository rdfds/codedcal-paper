"""Synthesize every repair outcome and the earned development conclusion."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from condition import HERE,ROOT,read,write,sha

def main():
    verified=read(HERE/'repair-verification.json');assert verified['verified_arms']==60
    rows=[]
    policies=[('recovery','barcode_model','Original author'),('structured-repair','structured','Pooling'),
              ('structured-repair','structured_author','Pooling + author'),('global-repair','global_group','Global search'),
              ('global-repair','global_author','Global + author')]
    for index in [0,1]:
        for kind in ['zero','common','color']:
            for cohort,policy,label in policies:
                for decoder in ['correlation','regression']:
                    p=HERE/cohort/f'field-{index}'/f'{kind}__{policy}__{decoder}.json';a=read(p)
                    assert a['status']=='completed'
                    row=dict(index=index,kind=kind,policy=policy,label=label,decoder=decoder,metrics=a['metrics'],quality=a['quality'],
                             residual_added_rms=a['residual_added_rms'],path=str(p.relative_to(ROOT)),sha256=sha(p))
                    row['development_target_met']=a['metrics']['5']['gene_class_identity_f1']>=.95 and (kind=='zero' or a['metrics']['5']['composition_tv']<=.05)
                    rows.append(row)
    global_rows=[r for r in rows if r['policy']=='global_group'];assert len(global_rows)==12
    assert all(r['development_target_met'] for r in global_rows)
    assert all(r['metrics']['5']['gene_class_identity_f1']==1 and r['metrics']['5']['composition_tv']==0 for r in global_rows)
    summaries=[]
    for cohort,policy,label in policies:
        selected=[r for r in rows if r['policy']==policy]
        summaries.append(dict(policy=policy,label=label,records=len(selected),targets_met=sum(r['development_target_met'] for r in selected),
            mean_nonzero_identity_f1=float(np.mean([r['metrics']['5']['gene_class_identity_f1'] for r in selected if r['kind']!='zero']))))
    runtime=[]
    for index in [0,1]:
        for kind in ['zero','common','color']:
            g=read(HERE/'global-repair'/f'field-{index}'/f'{kind}-geometry.json')
            runtime.append(dict(index=index,kind=kind,calibration_search_s=g['search']['runtime_s'],
                candidate_entries=sum(len(f['diagnostics']['grid']) for f in g['search']['families']),
                selected_round_initialization=g['search']['selected_round_initialization'],
                selected_partition=g['search']['selected_partition'],selected_shift_yx=g['search']['selected_shift_yx'],
                optional_author_refinement_s=None,note='Refinement time was not separately instrumented; no speed ratio.'))
    summary=dict(rows=rows,policies=summaries,global_group_records=12,global_group_exact_native_recovery=12,
        global_group_development_target_met=True,planned_global_arms=24,completed_global_arms=24,failed_global_arms=0,
        planned_pooling_arms=24,completed_pooling_arms=24,failed_pooling_arms=0,oracle_start_cases=6,
        synthetic_global_controls=22,synthetic_global_controls_exact=22,independent_transfer_established=False,
        temporal_intervention_validated=False,runtime=runtime,script_sha256=sha(__file__),verification_sha256=sha(HERE/'repair-verification.json'))
    write(HERE/'REPAIR_SUMMARY.json',summary)
    lines=['# Registration failure repaired — 6 October 2026','',
      '**The bounded global group search restores native outputs exactly in all twelve field/case/decoder combinations. This closes the exposed computational recovery failure, not independent validation or the Nature Methods acceptance gap.**','',
      'All global-search primary outcomes have identity F1 1.000, composition TV 0, the original yields (field 0: 70 correlation / 141 regression; field 1: 125 / 66), zero observed blank calls and complete scoring support. The complete 24-record global/refinement study and preceding 24-record pooling repair are retained. Every global primary decoder input is byte-identical to its native reference; calls are explicitly reused rather than counted as independent solver replications.','',
      '## What caused the failure','',
      'On field one’s color case, the centered failed geometry has padded author loss 0.496903 versus 0.385586 at the known inverse. Starting the unchanged optimizer at the inverse restores both native call lists exactly. The nine-point profile reveals a barrier in the padded objective along the line from the failed geometry to the inverse. This supports an initialization/local-basin limitation; it does not prove that overlaps or a single-code model mismatch caused it.','',
      'A separate limitation remains: field zero’s common-round oracle start drifts by one pixel in one of forty frames. Primary correlation/regression F1 becomes 0.9787/0.8114; label-constrained agreement and composition remain separate. A lower calibration loss can therefore change native outputs without proving improved or worsened physical RNA accuracy.','',
      '## What changed','',
      'Pooling rounds and colors repairs development initialization but still misses field one. Synthetic controls then expose destructive round-first alignment. The final search compares both zero-round and fitted-round starts, all three balanced color partitions, and fixed bounded shifts; it composes maps without discarding photons and ranks candidates using the original codebook objective. It receives no imposed maps, assignment labels, references or gene-call outcomes. Its accelerated grouped projections agree with actual author TensorFlow loss and independent direct float64 projections.','',
      'The search has 540 candidate entries per case, including repeated central refinement points. Rounds are bounded to ten pixels and opposite color-group shifts to nine; twenty-pixel padding preserves photons. It models balanced, static two-group color offsets plus a tested common-round initialization. It cannot represent arbitrary per-readout/nonrigid errors or certify that the model is sufficient.','',
      '## Complete color-case comparison','',
      '| Field | Decoder | Original author | Pooling | Pooling + author | Global search | Global + author |','|---|---|---:|---:|---:|---:|---:|']
    for index in [0,1]:
        for decoder in ['correlation','regression']:
            vals=[next(r for r in rows if r['index']==index and r['decoder']==decoder and r['kind']=='color' and r['policy']==p)['metrics']['5']['gene_class_identity_f1'] for _,p,_ in policies]
            lines.append(f"| {index} | {decoder} | "+' | '.join(f'{v:.4f}' for v in vals)+' |')
    lines+=['','All zero/common global-search outcomes are also exactly 1.000, preventing a color-only success claim. Optional author refinement preserves 11/12 combinations but damages field zero’s common case; it is not silently selected as the final policy. Original frozen selector and interior-fit failures remain. The corrected global controls pass 22/22; the prior six failures, failed phase initialization and all eight pooling controls are preserved.','',
      '## Verification and remaining test','',
      f"Retained calls, three radii, whole-domain maximum-cardinality/minimum-distance matching, blank counts, native-support masks, exact reuse, twenty array hashes and six independent global losses verify. Equal-cost assignments differ in {verified['equal_cost_identity_differences']}/180 radius-specific comparisons (largest observed F1 difference {verified['largest_observed_equal_cost_f1_change']:.5f}, not an exhaustive tie bound). Global primary exact-image/call recovery is unaffected by that ambiguity.",'',
      'These are post-hoc repairs on two fields of one exposed BARseq3 preparation, using controlled integer residuals and native algorithm references. The exposed amplitude lies on the coarse grid; odd-amplitude synthetic controls test refinement, not physical transfer. No independent RNA truth, prospective confirmation, robust universal recommendation, temporal T attribution or randomized intervention follows. The new global method is a stronger conventional alternative; unique early-calibration superiority is not inferred.','',
      'The appropriate next test freezes global_group without optional local refinement on a different public preparation. The already-known whole-cortex BARseq release supplies two metadata-qualified, previously unused D076 acquisitions, but its download host currently times out. D078-1L was used earlier in Source Data 23 and is excluded from fresh confirmation; the first two metadata records also have swapped sample identities. Exact metadata/access evidence and the calibration core freeze are preserved in outputs/barseq-public-transfer-audit-oct06/. No D076 candidate images have been examined. The adapter and full independent study are still pending; reserved kidney fields remain unopened.','',
      'Search runtimes and all policy/component outcomes are in REPAIR_SUMMARY.json. Optional refinement timing was not separately instrumented; reused kernel runtimes are not actual rerun times. No controlled speed claim is made. Failed S3/ISS, two temporal-estimator modality classes, fiducial-density qualifications, simulation-only randomized MERFISH, and all original strong alternatives remain unchanged.']
    (HERE/'REPAIR_RESULTS.md').write_text('\n'.join(lines)+'\n')
    # Complete panel: all conditions/policies, rather than only the favorable case.
    fig=plt.figure(figsize=(12,6),layout='constrained');grid=fig.add_gridspec(2,2,width_ratios=[1,1.45])
    ax=fig.add_subplot(grid[0,0]);d=read(HERE/'objective-diagnosis/field-1/color.json')
    for key,label in [('padded','Author padded objective'),('fixed_support_geometry_denominator','Fixed observed interior')]:
        ax.plot([p['fraction_toward_oracle'] for p in d['profile']],[p['objective'][key] for p in d['profile']],marker='o',ms=3,label=label)
    ax.set(xlabel='Fraction from failed geometry toward known inverse',ylabel='Single-code unexplained energy',title='a  Failed field: better basin beyond local fit');ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=fig.add_subplot(grid[1,0]);labels=[];old=[]
    for i in [0,1]:
        for decoder in ['correlation','regression']:
            labels.append(f"{i} / {'corr.' if decoder=='correlation' else 'reg.'}")
            old.append(next(r for r in rows if r['index']==i and r['kind']=='color' and r['decoder']==decoder and r['policy']=='barcode_model')['metrics']['5']['gene_class_identity_f1'])
    x=np.arange(4);ax.bar(x-.18,old,.36,label='Original author');ax.bar(x+.18,np.ones(4),.36,label='Global search',color='#238b45')
    ax.set(xticks=x,xticklabels=labels,ylim=(0,1.06),ylabel='Native identity F1');ax.set_title('b  Both fields and both decoder families',pad=34)
    ax.legend(fontsize=8,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.01));ax.grid(axis='y',alpha=.2)
    ax=fig.add_subplot(grid[:,1]);values=[];labels=[]
    for i in [0,1]:
        for kind in ['zero','common','color']:
            for decoder in ['correlation','regression']:
                labels.append(f"{i} / {kind} / {'corr.' if decoder=='correlation' else 'reg.'}")
                values.append([next(r for r in rows if r['index']==i and r['kind']==kind and r['decoder']==decoder and r['policy']==p)['metrics']['5']['gene_class_identity_f1'] for _,p,_ in policies])
    values=np.array(values);ax.imshow(values,vmin=0,vmax=1,cmap='YlGn',aspect='auto')
    ax.set(xticks=np.arange(5),xticklabels=['Original\nauthor','Pooling','Pooling\n+ author','Global\nsearch','Global\n+ author'],yticks=np.arange(12),yticklabels=labels,title='c  All cases; optional refinement harm retained')
    ax.tick_params(labelsize=8)
    for y in range(12):
        for x in range(5):ax.text(x,y,f'{values[y,x]:.3f}',ha='center',va='center',fontsize=8,color='white' if values[y,x]>.65 else 'black')
    ax.axhline(5.5,color='white',lw=2)
    fig.suptitle('Exposed-data repair: exact group recovery; independent transfer remains untested',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(HERE/f'repair-comparison.{ext}',dpi=200)
    plt.close(fig);print('Retained every repair outcome; global primary passes 12/12 development targets',flush=True)

if __name__=='__main__':main()
