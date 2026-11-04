"""Four tiny actual decoder controls, plus pure gate/validation tests."""
import json
from pathlib import Path
import sys
from unittest.mock import patch
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fixed_decoder import FixedDecoder,HistoricalGate,catalogs_reproduced
from codedcal.workflow_validation import ImageCase, Catalog
from serval.codebook import Codebook
from serval.image import ImageStack
from serval.decode.pixel import NearestNeigbourPixelDecoder,CosineOptimizedPixelDecoder
import pandas as pd


def test_actual_source_predictions_match_replay_and_anchored_scale_invariance():
    receipt=Path(__file__).resolve().parent/'software-control-receipt.json'
    assert not receipt.exists(),'Retain controls; do not silently re-execute native predictions'
    receipt.write_text(json.dumps(dict(stage='started',maximum_native_decoder_predictions=4))+'\n')
    images=np.zeros((2,1,12,12),np.float32)
    images[0,0,5,4]=8;images[0,0,5,5]=4;images[1,0,7,8]=12
    case=ImageCase(images,np.ones_like(images),np.eye(2),[0,0],[0,1],[10,11],[],
        metadata=dict(context=[[0,1],[0,12],[0,12]],scoring_region=[[0,1],[0,12],[0,12]],coverage_region=[[0,1],[0,12],[0,12]]))
    gains=np.array([2.,4.]);maps=np.zeros((2,3));scales=np.ones(2)
    # Actual source prediction with pre-set parameters, no CoSiD optimizer.
    source=CosineOptimizedPixelDecoder(Codebook(pd.DataFrame(np.eye(2))),
        NearestNeigbourPixelDecoder(Codebook(pd.DataFrame(np.eye(2))),max_dist=.5167,min_norm=1),init_scaling_factors=gains.copy())
    with patch.object(CosineOptimizedPixelDecoder,'fit',side_effect=AssertionError('No fit allowed')):
        table=source.predict(ImageStack(images[:,0].astype(float),fov=0,z=7)).spots
        expected=Catalog(np.c_[np.zeros(len(table)),table[['y','x']].to_numpy(float)],case.class_ids[table['barcode_id'].to_numpy(int)])
        engine=FixedDecoder(case,scales)
        historical=engine.run(case,maps,gains,'historical')
        anchored=engine.run(case,maps,gains,'anchored')
        multiplied=engine.run(case,maps,100*gains,'anchored')
    assert catalogs_reproduced(historical['catalog'],expected.as_dict(),0)
    np.testing.assert_array_equal(expected.classes,[10,11])
    np.testing.assert_allclose(expected.points,[[0,5,26/6],[0,7,8]],rtol=0,atol=1e-14)
    assert catalogs_reproduced(anchored['catalog'],multiplied['catalog'])
    np.testing.assert_array_equal(gains,[2.,4.])
    assert engine.predictions==3
    receipt.write_text(json.dumps(dict(passed=True,actual_native_decoder_predictions=4,
        registration_optimizations=0,gain_optimizations=0,source_predict_path_reproduced=True),indent=2)+'\n')


def test_incomplete_or_failed_historical_gate_blocks_new_outcomes():
    gate=HistoricalGate(['one','two'])
    gate.record('one',True)
    with pytest.raises(ValueError):gate.require_complete()
    gate.record('two',False)
    with pytest.raises(ValueError):gate.require_complete()


def test_duplicate_unknown_cases_rejected_and_complete_gate_opens():
    gate=HistoricalGate(['one','two'])
    gate.record('one',True)
    with pytest.raises(ValueError):gate.record('one',True)
    with pytest.raises(ValueError):gate.record('three',True)
    gate.record('two',True);gate.require_complete()


def test_catalog_gate_cannot_hide_changed_labels_or_coordinate_errors():
    original=dict(points=[[0,1,2]],classes=[5])
    assert catalogs_reproduced(original,original)
    assert not catalogs_reproduced(dict(points=[[0,1,2]],classes=[6]),original)
    assert not catalogs_reproduced(dict(points=[[0,1,2.0000001]],classes=[5]),original)
    assert not catalogs_reproduced(dict(points=[],classes=[]),original)
