"""Pinned, read-only historical object adapters plus actual native image kernels."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT.parent/'bardensr-upstream.nosync'),str(ROOT/'outputs/starfinder-public-confirmation-oct06'),str(ROOT/'tool/src')]
from codedcal.workflow_bardensr import BardensrAdapter
from codedcal.workflow_validation import Catalog
import adapter as A
import probe as P
import quality_catalog as Q
import tissue as T
import global_repair as G

class NativeAdapter(BardensrAdapter):
    def decode(self,case):
        kind=case.metadata.get('object_adapter')
        if kind is None:return super().decode(case)
        region=tuple(slice(*r) for r in case.metadata.get('calibration_region',[[0,case.images.shape[1]],[0,case.images.shape[2]],[0,case.images.shape[3]]]))
        images=case.images[(slice(None),)+region];support=case.support[(slice(None),)+region]>=1-1e-6
        if kind=='culture':raw=P.catalog(images,support,case.codebook,np.isin(case.class_ids,case.blank_classes),self.decoder,0.)
        elif kind=='tissue':raw=T.catalog(images,support,case.codebook,self.decoder)
        else:raise ValueError('Unknown frozen source object target')
        xyz=np.asarray(raw['xyz'],float).reshape(-1,3);words=np.asarray(raw['words'],int)
        unsupported=raw.get('unsupported_candidates',raw.get('runtime',{}).get('unsupported_candidates',0));candidates=raw['candidate_objects']
        coverage=candidates/(candidates+unsupported) if candidates+unsupported else 0.
        return Catalog(xyz[:,[2,1,0]]+np.array([s.start or 0 for s in region]),case.class_ids[words]),dict(native_kernel='BarDensr '+self.decoder,
            object_target=kind,source_abstentions=raw.get('abstained'),source_unsupported=raw.get('unsupported_candidates'),
            candidate_support_coverage=coverage,supported_source_candidates=candidates,unsupported_source_candidates=unsupported,
            original_whole_interface_failures_preserved=True)

    def group_fit(self,case,common):
        region=tuple(slice(*r) for r in case.metadata.get('calibration_region',[[0,case.images.shape[1]],[0,case.images.shape[2]],[0,case.images.shape[3]]]))
        image=case.images[(slice(None),)+region];alpha=case.support[(slice(None),)+region]
        projection=np.where(alpha>=1-1e-6,image,0).max(axis=1);support=alpha.max(axis=1)
        maps,info=G.fit(projection,support,case.codebook)
        return np.c_[np.zeros(len(maps)),maps],dict(calibration_fit_source='retained_balanced_group_search',integer_only_policy=True,group_loss=info['actual_author_loss'])


def make_adapter(config=None):return NativeAdapter(config)
