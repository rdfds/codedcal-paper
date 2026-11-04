"""Complete-case interpretation and sensitivities; never refits the frozen rule."""
import analyze as study
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

assert study.read(study.HERE/'comparison-summary.json')['successful_arms']==312
freeze=study.read(study.HERE/'preflight-freeze.json')
summary={}
fig,axs=plt.subplots(2,2,figsize=(9,7),constrained_layout=True)
colors={'round':'#3979b5','color':'#d36935','zero':'#777777'}
for split in ['development','validation']:
    records=study.records(split)
    summary[split]={}
    for col,decoder in enumerate(study.DECODERS):
        threshold=freeze['rules'][decoder]['selected']['threshold']
        cases=[]
        for obj in records:
            kind=obj['condition'].split('-')[0]
            arm=obj['arms'];get=lambda pol,r=5:arm[pol+'__'+decoder]['metrics'][str(r)]
            feature=study.feature(obj)
            selected='population_template' if feature>=threshold else 'common_iterative'
            metrics={pol:get(pol) for pol in study.POLICIES}
            gain=metrics['population_template']['identity_f1']-metrics['common_iterative']['identity_f1']
            if kind!='zero':
                cases.append(dict(condition=obj['condition'],kind=kind,feature=feature,selected=selected,
                    identity_f1={p:m['identity_f1'] for p,m in metrics.items()},
                    composition_tv={p:m['composition_tv'] for p,m in metrics.items()},
                    calls={p:m['calls'] for p,m in metrics.items()},gain=gain,
                    radius_gains={str(r):get('population_template',r)['identity_f1']-get('common_iterative',r)['identity_f1'] for r in [1.5,3,5]},
                    constrained_gain=get('population_template')['label_constrained_f1']-get('common_iterative')['label_constrained_f1'],
                    decision_core_s=sum(arm[p+'__'+decoder]['runtime']['calibration_s'] for p in ['population_template','common_iterative'])+
                        sum(arm[selected+'__'+decoder]['runtime'][k] for k in ['preprocessing_s','decoder_s','warp_s']),
                    always_template_core_s=sum(arm['population_template__'+decoder]['runtime'][k] for k in ['preprocessing_s','decoder_s','warp_s','calibration_s']),
                    model_core_s=sum(arm['barcode_model__'+decoder]['runtime'][k] for k in ['preprocessing_s','decoder_s','warp_s','calibration_s'])))
            axs[0,col].scatter(feature,gain,color=colors[kind],marker='o' if split=='validation' else 'x',alpha=.8,
                label=split+' '+kind if not any(line.get_label()==split+' '+kind for line in axs[0,col].collections) else None)
            if split=='validation' and kind!='zero':
                axs[1,col].scatter(metrics['common_iterative']['identity_f1'],metrics[selected]['identity_f1'],color=colors[kind])
        predictions=[c['selected']=='population_template' for c in cases]
        by_radius={str(r):study.confusion([c['radius_gains'][str(r)]>=.05 for c in cases],predictions) for r in [1.5,3,5]}
        means={p:float(np.mean([c['identity_f1'][p] for c in cases])) for p in study.POLICIES}
        selected_mean=float(np.mean([c['identity_f1'][c['selected']] for c in cases]))
        summary[split][decoder]=dict(cases=cases,threshold=threshold,confusion=by_radius['5'],
            radius_sensitivity=by_radius,constrained_matching_confusion=study.confusion([c['constrained_gain']>=.05 for c in cases],predictions),
            mean_identity_f1=means|{'selected':selected_mean},
            selected_minus_model=selected_mean-means['barcode_model'],
            harms_vs_common=sum(c['identity_f1'][c['selected']]-c['identity_f1']['common_iterative']<=-.05 for c in cases),
            mean_workflow_core_s={name:float(np.mean([c[name] for c in cases])) for name in ['decision_core_s','always_template_core_s','model_core_s']},
            native={p:next(o for o in records if o['condition'].startswith('zero'))['arms'][p+'__'+decoder]['metrics']['5'] for p in study.POLICIES})
for col,decoder in enumerate(study.DECODERS):
    ax=axs[0,col];ax.axhline(.05,color='black',ls=':');ax.axvline(freeze['rules'][decoder]['selected']['threshold'],color='black',ls='--')
    ax.set(xlabel='Observed correction disagreement (pixels RMS)',ylabel='Per-frame minus strong common-round\nidentity F1',title=decoder.capitalize())
    ax.legend(fontsize=7)
    ax=axs[1,col];ax.plot([0,1],[0,1],color='#888888',ls='--');ax.set(xlim=(0,1),ylim=(0,1),xlabel='Strong common-round identity F1',ylabel='Frozen recommendation identity F1',title='Annotated validation crop')
for ext in ['pdf','png']:fig.savefig(study.HERE/('decision-comparison.'+ext),dpi=180)
study.write(study.HERE/'interpretation.json',dict(studies=summary,
    independence='Two author-annotated crops; independent animal/preparation provenance not established.',
    runtime_scope='Calibration, preprocessing, warp and kernel only; excludes network, peak/sensitivity analysis; concurrent local hardware. Decision includes both calibration fits, conservatively double counts shared phase initialization.',
    interpretation='Controlled digital offsets and author annotation agreement; not prospective guidance, independent physical molecular truth or randomized-readout validation.'))
for split,value in summary.items():
    for decoder,entry in value.items():print(split,decoder,entry['confusion'],entry['mean_identity_f1'],flush=True)
