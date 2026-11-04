"""Actual source NN/connected-component decoder; no Ray/Dask or new thresholds."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(R/'tool/src'),str(R.parent/'serval-upstream.nosync')]
from codedcal.workflow_bardensr import BardensrAdapter
from codedcal.workflow_validation import Catalog
from serval.codebook import Codebook
from serval.image import ImageStack
from serval.decode.pixel import NearestNeigbourPixelDecoder
class SourceAdapter(BardensrAdapter):
 def __init__(self,config=None):
    self.config=dict(config or {});self.decoder='source_nn';self.cache={};self.geometry_cache={};self._native=None
 @property
 def native(self):
    if self._native is None:
        import bardensr
        self._native=bardensr
    return self._native
 def decode(self,case):
    context=self.context(case);a=case.images[(slice(None),)+context]
    if a.shape[1]!=1:raise ValueError('Source nearest adapter is qualified for the fixed single-plane target')
    book=Codebook(pd.DataFrame(case.codebook.T,index=np.arange(case.codebook.shape[1])))
    decoder=NearestNeigbourPixelDecoder(book,max_dist=.5167,min_norm=1)
    result=decoder.predict(ImageStack(a[:,0],fov=0,z=7));table=result.spots
    xyz=table[['y','x']].to_numpy(float);coords=np.c_[np.zeros(len(xyz)),xyz]+np.array([s.start or 0 for s in context]);words=table['barcode_id'].to_numpy(int)
    region=np.array(self.score_region(case));take=np.all((coords>=region[:,0])&(coords<region[:,1]),axis=1)
    coords=coords[take];words=words[take];order=np.lexsort((coords[:,2],coords[:,1],words))
    return Catalog(coords[order],case.class_ids[words[order]]),dict(native_kernel='Serval actual nearest-neighbor and native connected components',source_max_distance=.5167,source_min_norm=1,added_area_cutoff=None,
        native_components=int(len(table)),mean_decode_distance=float(table['mean_distance'].mean()) if len(table) else None)
def make_adapter(config=None):return SourceAdapter(config)
