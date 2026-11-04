"""Post-hoc common-canvas repair on exposed inputs; preserves original study."""
import argparse,gc,json,sys,time
from pathlib import Path
import numpy as np
import exposed as E
from codedcal import output_risk as O
from codedcal import output_risk_weighted as W
HERE=E.HERE/'common-domain'

def decode(images,support,maps,book,source,decoder,ids,negative):
    # Canonicalize BEFORE resampling, so all nonlinear kernels share one canvas.
    maps=O.canonical_maps(np.asarray(maps)[None])[0]
    out=E.T.move3d(images,maps);sup=E.T.move3d(support,maps)
    if source=='bar':
        xy,w,runtime=E.R.decode(out[:,0],book,decoder);points=np.c_[xy,np.zeros(len(xy))];raw=dict(runtime=runtime)
    elif source=='culture':
        cut=E.read(E.STAR/'quality-catalog/state.json')['groups'][decoder]['selected']['threshold'];raw=E.P.catalog(out,sup,book,negative,decoder,cut);points=np.asarray(raw['xyz']).reshape(-1,3);w=np.asarray(raw['words'],int)
    else:
        raw=E.T.catalog(out,sup,book,decoder);points=np.asarray(raw['xyz']).reshape(-1,3);w=np.asarray(raw['words'],int)
    cat=dict(xyz=points.tolist(),words=ids[w].tolist(),negative_ids=[] if negative is None else np.unique(ids[np.flatnonzero(negative)]).tolist())
    n=int(np.isin(cat['words'],cat['negative_ids']).sum());coding=len(w)-n
    quality=dict(full_catalog_calls=len(w),central_coding=coding,central_negative=n,blank_proxy=None if negative is None else n/max(coding,1)*int((~negative).sum())/int(negative.sum()),source_abstentions=raw.get('abstained'),source_unsupported=raw.get('unsupported_candidates'),runtime=raw.get('runtime'))
    return cat,quality

def native_input(case):
    source,index,kind=case
    if source=='bar':a,s=E.R.inputs(index);return a[:,None],s[:,None]
    path=E.A.CACHE/'local-trace-Position351/crop-4.npz' if source=='culture' else E.T.CACHE/'fullfield-conditioned-crop-4.npz'
    with np.load(path) as d:return d['images'],d['support']

def controls():
    E.L.guard();E.T.guard();records=[]
    for case in E.CASES:
        source,index,kind=case;obs,sup,book,physical,negative,ids,imposed,refs,sha=E.inputs(case);base,alpha=native_input(case);g=E.as3(imposed)[0];oracle=-E.as3(imposed);canonical=O.canonical_maps(oracle[None])[0]
        corrected=E.T.move3d(obs,canonical);correcteds=E.T.move3d(sup,canonical);expected=E.T.move3d(base,np.tile(g,(len(base),1)));expecteds=E.T.move3d(alpha,np.tile(g,(len(base),1)))
        region=(slice(None),slice(3,27),slice(16,112),slice(16,112)) if source!='bar' else (slice(None),slice(None),slice(16,112),slice(16,112))
        assert np.array_equal(corrected[region],expected[region]);assert np.array_equal(correcteds[region],expecteds[region]);assert np.array_equal(canonical,O.canonical_maps((oracle+np.array([1,2,-3]))[None])[0])
        decoders=['correlation'] if source=='bar' else ['source_spot','correlation']
        for dec in decoders:
            a,q=decode(obs,sup,oracle,book,source,dec,ids,negative);b,r=decode(expected,expecteds,np.zeros_like(oracle),book,source,dec,ids,negative)
            assert a==b,(case,dec,O.catalog_distance(a,b));assert len(a['words'])>=10,(case,dec,len(a['words']))
            records.append(dict(case=E.name(case),decoder=dec,calls=len(a['words']),coding=q['central_coding'],oracle_exact=True));print('CANVAS CONTROL',E.name(case),dec,len(a['words']),flush=True)
        gc.collect()
    E.write(HERE/'controls.json',dict(passed=True,records=records,exposed_development=True,independent_RNA_truth=False,nonnegative_backend_requalification=False,input_tensor_equivalence_in_full_kernel_context=True))

def guard():
    f=E.read(HERE/'freeze.json')
    for rel,digest in f['files'].items():assert E.sha(E.ROOT/rel)==digest,rel
    assert E.read(HERE/'controls.json')['passed'];E.L.guard();E.T.guard();return E.sha(HERE/'freeze.json')

def run(case):
    freeze=guard();source,index,kind=case;dest=HERE/'cases'/E.name(case);dest.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    if (dest/'evaluation.json').exists():assert E.read(dest/'evaluation.json')['freeze_sha256']==freeze;return
    observed,alpha,book,physical,negative,ids,imposed,refs,inputsha=E.inputs(case);maps,labels,ordinary,bases,provenance=E.bank(case);maps=O.canonical_maps(maps)
    region=((0,1),(24,104),(24,104)) if source=='bar' else ((4,26),(24,104),(24,104));block=(1,16,16) if source=='bar' else (2,16,16)
    profile=O.profile_losses(observed,alpha,physical,maps,region,block);comp=O.compatibility(profile);selection=W.selection_weights(profile)
    E.write(dest/'profile.json',dict(freeze_sha256=freeze,input_sha256=inputsha,geometry_sources=provenance,maps_zyx=maps.tolist(),labels=labels,conventional_indices=ordinary,incumbent_indices=bases,common_coverage=profile['coverage'],residual=profile['residual'].tolist(),energy=profile['energy'].tolist(),block_counts=profile['block_counts'].tolist(),compatibility=comp,selection=selection))
    if comp['abstained']:E.write(dest/'evaluation.json',dict(freeze_sha256=freeze,stage='abstained',reason=comp['reason'],development=True));print('DOMAIN ABSTAIN',E.name(case),comp['reason'],flush=True);return
    incumbent=min(bases,key=lambda k:(comp['losses'][k],float(np.square(maps[k]).sum()),k));decoders=['correlation'] if source=='bar' else ['source_spot','correlation'];rows=[]
    for dec in decoders:
        cats=[];quality=[]
        for k,m in enumerate(maps):
            path=dest/f'catalog-{dec}-{k}.json'
            if path.exists():d=E.read(path);assert d['freeze_sha256']==freeze;cat,q=d['catalog'],d['quality']
            else:cat,q=decode(observed,alpha,m,book,source,dec,ids,negative);E.write(path,dict(freeze_sha256=freeze,candidate=k,decoder=dec,catalog=cat,quality=q))
            cats.append(cat);quality.append(q)
        risk=O.output_risk(cats,comp,incumbent,ordinary);weighted=W.weighted_output_risk(cats,selection,incumbent,maps)
        if not weighted['abstained']:
            for key in ['ordinary_multistart_identity_disagreement','ordinary_multistart_composition_disagreement','unfiltered_bank_identity_disagreement','unfiltered_bank_composition_disagreement']:weighted[key]=risk[key]
        confidence=E.decoder_confidence(observed,maps[incumbent],book,source,dec,cats[incumbent]);pred=dict(decoder=dec,incumbent=incumbent,risk=risk,weighted_risk=weighted,baseline_fit_loss=comp['losses'][incumbent],baseline_fit_gap=comp['losses'][incumbent]-min(comp['losses']),baseline_decoder_confidence=confidence,baseline_quality=quality[incumbent],freeze_sha256=freeze);E.write(dest/f'prediction-{dec}.json',pred)
        # Evaluate only after saving observable forecasts, on exactly this canvas.
        base,sup=native_input(case);g=E.as3(imposed)[0];reference,rq=decode(E.T.move3d(base,np.tile(g,(len(base),1))),E.T.move3d(sup,np.tile(g,(len(base),1))),np.zeros_like(maps[0]),book,source,dec,ids,negative)
        actual=O.catalog_distance(reference,cats[incumbent]);E.write(dest/f'reference-{dec}.json',dict(catalog=reference,quality=rq,evaluation_only=True,known_global_gauge=g.tolist()))
        rows.append(pred|dict(reference_calls=len(reference['words']),incumbent_calls=len(cats[incumbent]['words']),actual_consistency_error=actual,failed_consistency=actual['identity_loss']>.05 or actual['composition_tv']>.05));print('DOMAIN RISK',E.name(case),dec,'calls',len(cats[incumbent]['words']),'actual',actual['identity_loss'],'max',risk.get('identity_risk'),'weighted',weighted.get('identity_risk'),flush=True);gc.collect()
    E.write(dest/'evaluation.json',dict(stage='complete',case=E.name(case),source=source,field=index,kind=kind,freeze_sha256=freeze,rows=rows,runtime_s=time.monotonic()-started,development=True,independent_confirmation=False,independent_RNA_truth=False))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true');args=parser.parse_args()
    if args.controls:controls()
    else:
        for case in E.CASES:run(case)
