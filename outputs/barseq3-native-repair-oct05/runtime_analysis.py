"""Account for both model starts and both decision fits; no portable speed claim."""
from condition import HERE,read,write,sha

s=read(HERE/'recovery/summary.json');assert s['full_cohort_comparison_supported']
rows=[];native=[]
for i in [0,1]:
    geometry=sum(read(HERE/f'field-{i}/round-{r}-geometry.json')['elapsed_s'] for r in range(10))
    conditioning=sum(read(HERE/f'field-{i}/round-{r}-demons.json')['elapsed_s'] for r in range(10))
    native.append(dict(field=i,full_field_geometry_s=geometry,recorded_conditioning_s=conditioning,
        scope='Shared native preparation; recorded stage times do not include all download/profile/loading overhead.'))
    folder=HERE/'recovery'/f'field-{i}'
    for kind in ['zero','common','color']:
        geo=read(folder/f'{kind}-geometry.json');g=geo['runtime']
        costs=dict(none=0.,oracle=0.,common_iterative=g['shared_phase_s']+g['common_s'],
            population_template=g['shared_phase_s']+g['frame_s'],
            barcode_model=g['shared_phase_s']+g['frame_s']+g['model_s'],
            selected=g['shared_phase_s']+g['common_s']+g['frame_s'])
        for d in ['correlation','regression']:
            for p in ['none','common_iterative','population_template','barcode_model','oracle','selected']:
                chosen=geo['recommend'] if p=='selected' else p
                arm=read(folder/f'{kind}__{chosen}__{d}.json');r=arm['runtime']
                reused=r.get('exact_reference_kernel_reused',False)
                # Reference cost remains explicit even when the execution was
                # skipped after exact input equality. Never rank the oracle by
                # its zero cached execution time.
                measured=r.get('actual_execution_kernel_s',r['kernel_s'])+(0. if reused else r['peaks_s'])+r['warp_s']
                standalone=costs[p]+r['kernel_s']+r['peaks_s']+r['warp_s']
                rows.append(dict(field=i,kind=kind,decoder=d,policy=p,chosen=chosen,
                    required_fit_s=costs[p],recorded_execution_with_fit_s=costs[p]+measured,
                    standalone_component_cost_s=standalone,reference_kernel_s=r['kernel_s'] if reused else None,
                    exact_reference_reused=reused,quality_pass=arm['quality']['qualified'],
                    note='Selected requires both maps; model includes both starts plus its per-frame initialization. Standalone totals omit shared native preparation/score and are descriptive component accounting.'))
write(HERE/'recovery/runtime-accounting.json',dict(rows=rows,native_shared=native,
    original_scientific_arms=60,policy_accounting_rows=len(rows),script_sha256=sha(__file__),
    conclusion='Local runtime bookkeeping, not controlled speed superiority. Cached oracle/reference arms are not independent solver executions.'))
print('Accounted for 60 arms and 12 selected-policy rows; both model starts and both decision maps included',flush=True)
