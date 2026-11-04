"""Declared codebook-only stratification; preserve every planned recovery arm."""
import sys
import numpy as np
from condition import HERE,ROOT,read,write,sha
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.evaluation import spatial_match

def main():
    summary=read(HERE/'recovery/summary.json')
    assert summary['full_cohort_comparison_supported'],'No complete-cohort mechanism claim from partial results'
    classes=read(HERE/'codebook-classes.json')
    labels=np.array(classes['names']);blank=np.array(classes['blank'],bool)
    k=np.array([sum(c in 'GC' for c in s) for s in classes['sequences']])
    balanced=np.maximum(k,10-k)<=6
    ledger=[]
    for field in [0,1]:
        folder=HERE/'recovery'/f'field-{field}'
        for decoder in ['correlation','regression']:
            ref=read(folder/f'reference-{decoder}.json');a=np.array(ref['xy'],float).reshape(-1,2);aw=np.array(ref['words'],int)
            for kind in ['zero','common','color']:
                for policy in ['none','common_iterative','population_template','barcode_model','oracle']:
                    path=folder/f'{kind}__{policy}__{decoder}.json';arm=read(path)
                    assert arm['status']=='completed'
                    b=np.array(arm['xy'],float).reshape(-1,2);bw=np.array(arm['words'],int)
                    radius_records={}
                    for radius in [1.5,3,5]:
                        pairs,diag=spatial_match(a,b,radius)
                        assert diag['skipped_components']==0
                        correct=labels[aw[pairs[:,0]]]==labels[bw[pairs[:,1]]]
                        seqcorrect=aw[pairs[:,0]]==bw[pairs[:,1]]
                        counts={}
                        for name,mask in [('balanced',balanced),('unbalanced',~balanced)]:
                            native=mask[aw]&~blank[aw];predicted=mask[bw]&~blank[bw]
                            matched_group=native[pairs[:,0]]
                            nref=int(native.sum());ngene=int((matched_group&correct).sum());nseq=int((matched_group&seqcorrect).sum())
                            counts[name]=dict(native_reference_calls=nref,predicted_nonblank_calls=int(predicted.sum()),
                                gene_consistent_matches=ngene,sequence_consistent_matches=nseq,
                                missed_or_gene_inconsistent_reference_calls=nref-ngene,
                                native_gene_consistency_recall=ngene/nref if nref else None,
                                native_sequence_consistency_recall=nseq/nref if nref else None)
                        radius_records[str(radius)]=counts
                    ledger.append(dict(field=field,decoder=decoder,kind=kind,policy=policy,
                        groups=radius_records,quality=arm['quality'],arm_sha256=sha(path)))
    assert len(ledger)==60
    write(HERE/'recovery/readout-stratification.json',dict(arms=ledger,complete_arms=60,
        group_definition='Public codebook only: balanced=max(GC_count,10-GC_count)<=6; other classes unbalanced. Blank classes excluded from group recall but preserved in arm quality.',
        reference_scope='Recall of native algorithm outputs, not independent molecular truth. Group membership follows reference sequence; gene consistency pools synonymous sequence classes.',
        specification_sha256=sha(HERE/'READOUT_STRUCTURE_ANALYSIS.md'),script_sha256=sha(__file__),
        codebook_sha256=classes['codebook_sha256']))
    lines=['# Declared readout-structure diagnostic','',
        'Native-reference consistency, not biological sensitivity or true FDR. Groups were fixed from the public codebook before perturbed recovery outcomes.','',
        '| Field | Decoder | Error grouping | Policy | Balanced reference recall | Unbalanced reference recall | Blank proxy |',
        '|---|---|---|---|---:|---:|---:|']
    def fmt(x):return 'undefined' if x is None else f'{x:.4f}'
    for r in ledger:
        if r['kind']=='zero':continue
        g=r['groups']['5']
        lines.append(f"| {r['field']} | {r['decoder']} | {r['kind']} | {r['policy']} | {fmt(g['balanced']['native_gene_consistency_recall'])} | {fmt(g['unbalanced']['native_gene_consistency_recall'])} | {fmt(r['quality']['proxy'])} |")
    lines+=['','All 60 arms, including zero cases, all group counts and all three radii are retained in readout-stratification.json. A low or zero blank proxy is not an independent label check. Restricted-support bounds need not describe actual pixels with PSF tails, overlap and background. Synonymous barcodes may cross the camera-balance groups, so gene and exact-sequence recalls are both retained.']
    (HERE/'recovery/READOUT_RESULTS.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
