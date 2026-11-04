"""Constructive codebook counterexample, not an estimate of empirical frequency."""
import csv
from itertools import combinations
import numpy as np
from condition import HERE,sha,write

rows=list(csv.DictReader((HERE.parent/'independent-comparison-oct05/barseq3-adult-codebook.csv').open()))
groups={}
for r in rows:groups.setdefault(r['Barcode'],set()).add(r['Name'])
codes=list(groups);names=['|'.join(sorted(groups[s])) for s in codes]
blank=np.array([all('no_gene' in n.lower() for n in groups[s]) for s in codes])
B=np.zeros((40,len(codes)),bool)
for j,s in enumerate(codes):
    for r,c in enumerate(s):B[r*4+'GTCA'.index(c),j]=True
camera=np.tile([True,False,True,False],10)
partial=B&camera[:,None];k=partial.sum(axis=0);eligible=(k>=4)&(k<=6)&~blank
found=None
for target in np.flatnonzero((k>=7)&~blank):
    # Parents may contribute only the target's positive-camera support.
    candidates=np.flatnonzero(eligible&~np.any(partial&~partial[:,target,None],axis=0))
    for a,b in combinations(candidates,2):
        if groups[codes[target]]&(groups[codes[a]]|groups[codes[b]]):continue
        if groups[codes[a]]&groups[codes[b]]:continue
        x=.5*(partial[:,a].astype(float)+partial[:,b])
        evidence=float(x@B[:,target]/((np.linalg.norm(x)+.05)*np.sqrt(10)))
        if evidence<=.775:continue
        union=B[:,a]|B[:,b]
        supported=np.flatnonzero(~np.any(B&~union[:,None],axis=0))
        # With known gains and noiseless complete nonnegative signal, no other
        # class can explain the zero frames: the two parents are identifiable.
        if set(supported)!={int(a),int(b)}:continue
        found=dict(parent_indices=[int(a),int(b)],parent_names=[names[a],names[b]],parent_codes=[codes[a],codes[b]],
            absent_target_index=int(target),absent_target_name=names[target],absent_target_code=codes[target],
            predicted_false_target_evidence=evidence,isolated_point_signal=x.tolist(),
            complete_point_signal=(.5*(B[:,a].astype(float)+B[:,b])).tolist(),
            complete_support_permitted_classes=supported.tolist(),parent_camera_support=[int(k[a]),int(k[b])],
            target_camera_support=int(k[target]))
        break
    if found:break
assert found is not None,'No qualifying constructive example found'
found.update(protocol='Deterministic first lexicographic witness in public codebook order satisfying declared analytic constraints. This is a selected counterexample, not a sampled performance benchmark.',
    gain_assumption='Known gains/noiseless positive mixture for complete-signal identifiability; not a claim that the gain-learning native solver always recovers it.',
    all_blank_codes_camera_counts=k[blank].astype(int).tolist(),all_blanks_below_correlation_cut_for_any_positive_camera_only_signal=bool(np.all(np.sqrt(np.maximum(k[blank],10-k[blank])/10)<.775)),
    codebook_sha256=sha(HERE.parent/'independent-comparison-oct05/barseq3-adult-codebook.csv'),script_sha256=sha(__file__),before_actual_kernel_validation=True)
write(HERE/'counterexample-specification.json',found)
print(found['parent_names'],'-> absent',found['absent_target_name'],'predicted evidence',found['predicted_false_target_evidence'],flush=True)
