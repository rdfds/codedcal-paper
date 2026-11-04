"""Declared strong hybrid: native geometric model plus actual source gain fit."""
from pathlib import Path
import sys,random
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path.insert(0,str(H))
from native_serval_adapter import SourceAdapter
from codedcal.workflow_validation import ImageCase,Catalog
from serval.codebook import Codebook
from serval.image import ImageStack
from serval.decode.pixel import NearestNeigbourPixelDecoder,CosineOptimizedPixelDecoder

class FullSourceAdapter(SourceAdapter):
 def physical_case(self,case):
    keep=~np.isin(case.class_ids,case.blank_classes)
    if not keep.any():raise ValueError('No physical codewords for geometric calibration')
    return ImageCase(case.images,case.support,case.codebook[:,keep],case.rounds,case.colors,case.class_ids[keep],np.array([],int),case.name,case.preparation,case.metadata)
 def calibrate(self,case,policy):
    fit_case=self.physical_case(case);maps,info=super().calibrate(fit_case,policy)
    return maps,dict(info,calibration_physical_codewords=fit_case.codebook.shape[1],named_blank_codewords_excluded=case.codebook.shape[1]-fit_case.codebook.shape[1])
 def objective(self,images,book):return super().objective(images,self._physical_book)
 def run(self,case,policy):
    self._physical_book=self.physical_case(case).codebook
    return super().run(case,policy)
 def decode(self,case):
    scales=np.asarray(self.config['source_native_scales'],float)
    if scales.shape!=(len(case.images),) or not np.isfinite(scales).all() or np.any(scales<=0):raise ValueError('Invalid source-unit scales')
    context=self.context(case);a=case.images[(slice(None),)+context]
    if a.shape[1]!=1:raise ValueError('Declared source workflow is single-plane')
    raw=ImageStack(a[:,0].astype(float)*scales[:,None,None],fov=0,z=7)
    book=Codebook(pd.DataFrame(case.codebook.T,index=np.arange(case.codebook.shape[1])))
    fitting_nn=NearestNeigbourPixelDecoder(book,max_dist=.5167,min_norm=0)
    decoder=CosineOptimizedPixelDecoder(book,fitting_nn,init_scaling_factors=scales/scales.mean(),min_area=5,max_area=12,max_size=100000,penalty_entropy=0,penalty_l2=0)
    random.seed(110126)
    for _ in range(10):decoder.fit([raw])
    decoder.decoder=NearestNeigbourPixelDecoder(book,max_dist=.5167,min_norm=1)
    table=decoder.predict(raw).spots
    xy=table[['y','x']].to_numpy(float);points=np.c_[np.zeros(len(xy)),xy]+np.array([s.start or 0 for s in context]);words=table['barcode_id'].to_numpy(int)
    region=np.array(self.score_region(case));take=np.all((points>=region[:,0])&(points<region[:,1]),axis=1)
    points=points[take];words=words[take];order=np.lexsort((points[:,2],points[:,1],words))
    return Catalog(points[order],case.class_ids[words[order]]),dict(native_kernel='actual Serval CoSiD10 plus native NN/connected components',hybrid_geometric_workflow=True,
        fitted_source_gains=decoder.scaling_factors,source_max_distance=.5167,source_final_min_norm=1,added_output_area_cutoff=None,
        source_fit_iterations=10,source_fit_min_area=5,source_fit_max_area=12,source_fit_max_size=100000,
        mean_decode_distance=float(table['mean_distance'].mean()) if len(table) else None)

def make_adapter(config=None):return FullSourceAdapter(config)
