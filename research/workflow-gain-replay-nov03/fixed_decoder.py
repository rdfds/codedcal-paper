"""Actual source inference with retained parameters; no image optimization."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

H=Path(__file__).resolve().parent;R=H.parents[1]
sys.path[:0]=[str(R/'tool/src'),str(R.parent/'serval-upstream.nosync')]
from codedcal.workflow_validation import Catalog,apply_corrections
from codedcal.workflow_gain import anchor_relative_gains
from serval.codebook import Codebook
from serval.image import ImageStack
from serval.decode.pixel import NearestNeigbourPixelDecoder,CosineOptimizedPixelDecoder


def catalogs_reproduced(actual,expected,tolerance=1e-12):
    a,b=Catalog.from_dict(actual),Catalog.from_dict(expected)
    return bool(np.array_equal(a.classes,b.classes) and a.points.shape==b.points.shape
                and np.allclose(a.points,b.points,rtol=0,atol=tolerance))


class HistoricalGate:
    def __init__(self,expected):
        self.expected=set(expected);self.seen=set();self.failed=False
    def record(self,case_id,matched):
        if case_id not in self.expected or case_id in self.seen:
            raise ValueError('Unexpected or repeated historical case')
        self.seen.add(case_id);self.failed|=not matched
    def require_complete(self):
        if self.failed or self.seen!=self.expected:
            raise ValueError('Complete historical reproduction required before anchored outcomes')


class FixedDecoder:
    """A fixed-codebook NN index plus the real CoSiD prediction method.

    The source fit method is never called. Its prediction method uses only the
    assigned scaling_factors and fixed underlying NN decoder. Original catalog
    extraction and scalar operation order are retained.
    """
    def __init__(self,case,scales):
        self.scales=np.asarray(scales,float)
        if self.scales.shape!=(len(case.images),) or not np.isfinite(self.scales).all() or np.any(self.scales<=0):
            raise ValueError('Invalid frozen source scales')
        self.codebook=case.codebook.copy();self.class_ids=case.class_ids.copy()
        book=Codebook(pd.DataFrame(case.codebook.T,index=np.arange(case.codebook.shape[1])))
        nn=NearestNeigbourPixelDecoder(book,max_dist=.5167,min_norm=1)
        self.decoder=CosineOptimizedPixelDecoder(book,nn,init_scaling_factors=np.ones(len(case.images)),
            min_area=5,max_area=12,max_size=100000,penalty_entropy=0,penalty_l2=0)
        self.reference=self.scales/self.scales.mean()
        self.predictions=0;self.fixed_codebook_index_builds=1

    def run(self,case,maps,gains,mode):
        if mode not in {'historical','anchored'}:raise ValueError('Unknown fixed-fit inference mode')
        if not np.array_equal(case.codebook,self.codebook) or not np.array_equal(case.class_ids,self.class_ids):
            raise ValueError('Frozen codebook/class mapping changed')
        gains=np.asarray(gains,float)
        if gains.shape!=self.scales.shape or not np.isfinite(gains).all() or np.any(gains<=0):
            raise ValueError('Invalid retained gains')
        used=gains.copy();gauge=None
        if mode=='anchored':used,gauge=anchor_relative_gains(gains,self.reference)
        aligned,canonical_maps=apply_corrections(case.for_adapter(),maps)
        z,y,x=aligned.images.shape[1:]
        context=tuple(slice(*r) for r in aligned.metadata.get('context',[[0,z],[16,y-16],[16,x-16]]))
        data=aligned.images[(slice(None),)+context]
        if data.shape[1]!=1:raise ValueError('Frozen source inference is single-plane')
        raw=ImageStack(data[:,0].astype(float)*self.scales[:,None,None],fov=0,z=7)
        self.decoder.scaling_factors=used
        self.predictions+=1
        table=self.decoder.predict(raw).spots
        xy=table[['y','x']].to_numpy(float)
        points=np.c_[np.zeros(len(xy)),xy]+np.array([s.start or 0 for s in context])
        words=table['barcode_id'].to_numpy(int)
        region=np.array(aligned.metadata.get('scoring_region',[[0,z],[24,y-24],[24,x-24]]))
        take=np.all((points>=region[:,0])&(points<region[:,1]),axis=1)
        points=points[take];words=words[take];order=np.lexsort((points[:,2],points[:,1],words))
        catalog=Catalog(points[order],aligned.class_ids[words[order]])
        coverage_region=tuple(slice(*r) for r in aligned.metadata.get('coverage_region',[[s.start,s.stop] for s in context]))
        coverage=float(np.all(aligned.support[(slice(None),)+coverage_region]>=1-1e-6,axis=0).mean())
        return dict(catalog=catalog.as_dict(),coverage=coverage,corrections=canonical_maps,
            gains_used=used,gain_gauge=gauge,decoder_only_executions=1,
            calibration_optimizations=0,gain_optimizations=0,
            original_source_max_distance=.5167,original_source_min_norm=1)
