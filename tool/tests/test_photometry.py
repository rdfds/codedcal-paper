import numpy as np
import pytest
from codedcal.photometry import fit_scales, score_traces, select_calls, sample_traces
W = np.array([[1,1,1,1,0,0,0,0], [0,0,0,0,1,1,1,1]])

def test_known_missing_bit_ties_and_label_symmetric_selection():
    x = np.array([W[0], W[1], [1,1,1,0,0,0,0,0], np.ones(8), W[0]], float)
    score = score_traces(x, W, np.ones(8))
    assert .51 < score['distance'][2] < .53 and score['tie'][3]
    xy = np.array([[0,0],[10,0],[20,0],[30,0],[.1,.1]])
    ids, counts = select_calls(xy, score, np.arange(5), 2, .6, .25)
    np.testing.assert_array_equal(ids, [0,1,2]);np.testing.assert_array_equal(counts, [2,1])
    _, switched = select_calls(xy, score_traces(x, W[::-1], np.ones(8)), np.arange(5), 2, .6, .25)
    np.testing.assert_array_equal(switched, counts[::-1])
    perm = [4,3,2,1,0]
    inds, _ = select_calls(xy[perm], score_traces(x[perm], W, np.ones(8)), np.array(perm), 2, .6, .25)
    np.testing.assert_array_equal(np.array(perm)[inds], ids)

def test_training_only_scaling_and_unsupported_abstention():
    x = np.vstack([W[0]*2., W[1]*3., W[0]*1000.]);ids = np.array([0,1,0]);train = np.array([True,True,False])
    scales, support = fit_scales(x, ids, W, train, min_support=1)
    np.testing.assert_array_equal(scales, [2]*4+[3]*4);np.testing.assert_array_equal(support, np.ones(8))
    with pytest.raises(ValueError, match='unsupported scale'):fit_scales(x, ids, W, train, min_support=2)

def test_inverse_affine_sampling_geometry_and_stream_length():
    image = np.zeros((17,17), np.float32);image[8,8] = 100
    C = np.array([[1,0],[0,1],[5,-3.]])
    traces, good = sample_traces(iter([image]), np.array([[13,5],[5,-3]]), C[None])
    assert traces[0,0] > 0 and good.tolist() == [True,False]
    shifted, _ = sample_traces([image], np.array([[8,8]]), np.array([[[1,0],[0,1],[0,0]]]))
    assert traces[0,0] == shifted[0,0]
    with pytest.raises(ValueError, match='fewer images'):sample_traces([], np.array([[8,8]]), C[None])
    with pytest.raises(ValueError, match='singular'):sample_traces([image], np.array([[8,8]]), np.zeros((1,3,2)))

@pytest.mark.parametrize('bad', [np.zeros(8), np.ones(7), np.full(8,np.nan)])
def test_invalid_scales_cannot_silently_decode(bad):
    with pytest.raises(ValueError):score_traces(W, W, bad)

def test_invalid_book_and_zero_trace_abstention():
    with pytest.raises(ValueError, match='unique'):score_traces(W, np.vstack([W[0],W[0]]), np.ones(8))
    with pytest.raises(ValueError, match='equal'):score_traces(W, np.array([[1,0]*4,[1]*8]), np.ones(8))
    assert score_traces(np.zeros((1,8)), W, np.ones(8))['tie'][0]
