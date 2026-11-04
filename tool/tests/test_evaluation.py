import numpy as np
import pytest
from codedcal.evaluation import evaluate_calls, spatial_match, warp_readouts, identity_affines


def test_balanced_swaps_do_not_pass_identity_despite_preserved_composition():
    xy=np.array([[0.,0.],[10.,0.]])
    r=evaluate_calls(xy,np.array([0,1]),xy,np.array([1,0]),np.array([False,False,True]))
    assert r['identity_f1']==0 and r['composition_tv']==0 and r['blank_proxy']==0


def test_matching_maximizes_cardinality_without_label_information():
    p,d=spatial_match(np.array([[0.,0.],[2.,0.]]),np.array([[.9,0.],[-1.1,0.]]),1.2)
    assert set(map(tuple,p))=={(0,1),(1,0)}
    assert d['ambiguous_prediction']==1


def test_empty_calls_are_missing_signal_not_zero_risk():
    r=evaluate_calls(np.array([[0.,0.]]),np.array([0]),np.empty((0,2)),np.array([],int),np.array([False,True]))
    assert r['identity_f1']==0 and r['blank_proxy'] is None and r['identity_precision'] is None


def test_oversized_components_remain_in_denominator():
    xy=np.zeros((3,2))
    pairs,d=spatial_match(xy,xy,max_component=2)
    assert len(pairs)==0 and d['skipped_components']==1


def test_real_reference_does_not_manufacture_truth():
    r=evaluate_calls(np.array([[0.,0.]]),np.array([0]),np.array([[0.,0.]]),np.array([0]),np.array([False,True]),reference_kind='unperturbed_call_reference')
    assert r['identity_f1']==1 and r['false_call_fraction'] is None


def test_affine_identity_and_singular_rejection():
    image=np.arange(128).reshape(2,8,8)
    assert np.array_equal(warp_readouts(image,identity_affines(2),(100,200)),image)
    with pytest.raises(ValueError):warp_readouts(image,np.zeros((2,3,2)))


def test_spatial_assignment_agrees_with_exhaustive_small_graphs():
    from itertools import combinations,permutations
    rng=np.random.default_rng(73)
    for _ in range(30):
        a=rng.uniform(0,3,(4,2));b=rng.uniform(0,3,(4,2));radius=1.5
        dist=np.linalg.norm(a[:,None]-b[None],axis=2)
        best=(0,0.)
        for size in range(1,5):
            for rows in combinations(range(4),size):
                for cols in permutations(range(4),size):
                    ds=dist[rows,cols]
                    if np.all(ds<=radius):
                        candidate=(size,-float(ds.sum()))
                        if candidate>best:best=candidate
        pairs,_=spatial_match(a,b,radius)
        assert len(pairs)==best[0]
        assert -dist[pairs[:,0],pairs[:,1]].sum()==pytest.approx(best[1])


def test_affine_composition_matches_direct_coordinates_and_image_ramp():
    from codedcal.evaluation import compose_affines
    a=identity_affines(1);a[0,:2]=[[1.1,.02],[-.03,.9]];a[0,2]=[.4,-.6]
    b=identity_affines(1);b[0,:2]=[[.97,-.02],[.01,1.04]];b[0,2]=[-.5,.7]
    c=compose_affines(a,b)
    xy=np.array([[10.,20.],[14.,22.]])
    assert np.allclose((xy@a[0,:2]+a[0,2])@b[0,:2]+b[0,2],xy@c[0,:2]+c[0,2])
    yy,xx=np.mgrid[:32,:32];origin=np.array([100.,200.])
    image=(2*(xx+origin[0])+3*(yy+origin[1]))[None].astype(float)
    transform=identity_affines(1);transform[0,2]=[1.2,-.7]
    result=warp_readouts(image,transform,origin)
    expected=2*(xx+origin[0]-1.2)+3*(yy+origin[1]+.7)
    assert np.allclose(result[0,3:-3,3:-3],expected[3:-3,3:-3],atol=1e-4)
