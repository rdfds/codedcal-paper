import numpy as np
from prepare_source_fullfield import stream_bit
from serval.image import ImageStack
from serval.transform import HighPassImageTransform,DeconvoleImageTransform,LowPassImageTransform


def test_streamed_native_filters_match_author_batch_for_both_dtype_modes():
    data=np.random.default_rng(3).integers(0,300,(3,64,64),dtype=np.uint16)
    transforms=[HighPassImageTransform(sigma=2),DeconvoleImageTransform(filter_size=9,num_iters=10,sigma=1.2),LowPassImageTransform(sigma=1.)]
    for preserve in [True,False]:
        batch=ImageStack(data)
        for t in transforms:batch=t.transform(batch,preserve_dtype=preserve)
        streamed=np.stack([stream_bit(a,i,transforms,preserve) for i,a in enumerate(data)])
        np.testing.assert_array_equal(streamed,batch.imgs)
