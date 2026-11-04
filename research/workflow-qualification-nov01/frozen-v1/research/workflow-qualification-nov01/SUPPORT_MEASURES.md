# Match support to the declared computational target

Pre-comparison native v2 showed a measurement mismatch: the historical tissue
object workflow has99.425% eligible source candidates, but98.060% supported voxels
in the bounding box. The new wrapper incorrectly substituted the latter for the
former. Object catalogs are conditioned on source candidates; pixel catalogs are
scored over their voxel domain. Native v3 therefore reports BOTH measurements
and applies the unchanged99% gate to the declared target's measure. It never
calls the entire tissue volume99%-supported. Culture candidate coverage is also
measured, including its previously discarded unsupported candidates; this can
remove a culture backend that the erroneous voxel-only gate admitted.

Fractional challenges form a separate panel. Their analytic composed-operator
control must pass, but fractional support loss does not erase a qualified integer
panel. Record fractional support as its own outcome; never call interpolation an
exact inverse. v2 failures remain intact. These are pre-comparison adapter and
target-measure repairs; no method/challenge comparative outcome informed them.
