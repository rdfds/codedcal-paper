"""Pinned case reconstruction and original-reference scoring; no optimization."""
from pathlib import Path
import json
import sys
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[1]
P=R/'research/workflow-validation-completion-nov03'
sys.path[:0]=[str(R/'tool/src'),str(R/'research/workflow-validation-controller-nov02')]
from common import sha,cache_key,blank_bound
from codedcal.workflow_validation import (ImageCase,Catalog,WorkflowResult,warp,
    canonical,pad_for_challenges,assess,compare,write_json)


def load(path):return json.loads(Path(path).read_text())


def observed(case,challenge):
    if challenge is None:return case.for_adapter()
    maps=np.asarray(challenge['maps'],float)
    assert np.array_equal(maps,np.round(maps)),'Integer challenges only'
    return case.replace_images(warp(case.images,maps),warp(case.support.astype(np.float32),maps)).for_adapter()


def retained(entry):
    row=load(R/entry['source_path'])
    return row if entry['challenge'] is None else row['runs'][0]


def score(result,reference,case,challenge,native_maps):
    reference=Catalog.from_dict(reference['catalog'])
    r=WorkflowResult(Catalog.from_dict(result['catalog']),np.asarray(result['corrections']),result['coverage'])
    assessment=assess(reference,r,case.blank_classes,minimum_reference=10,minimum_coverage=.99,radius=3.,tolerance=.05)
    if challenge is not None:
        assessment['matching_sensitivity']={str(rad):compare(reference,r.catalog,case.blank_classes,rad) for rad in [1.5,5.]}
        residual=canonical(r.corrections+np.asarray(challenge['maps'])-np.asarray(native_maps))
        assessment['known_probe_map_rms']=float(np.sqrt(np.mean(np.sum(residual**2,axis=1))))
    return assessment


def verify_manifest(path):
    rows=load(path)['files']
    for row in rows:
        p=R/row['path']
        assert p.is_file() and not p.is_symlink() and p.stat().st_size==row['bytes'],row['path']
        assert sha(p)==row['sha256'],row['path']
    return len(rows)
