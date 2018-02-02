dtw_c
====

Dynamic Time Warping (DTW) implementation in C for Python

Perform DTW of source features (org) into target features (trg) with mel-cepstral distortion (mcd) distance

Install
----

pip install dtw_c

Usage of dtw_org_to_trg
----

from dtw_c import dtw_c as dtw

dtw_org, twf_mat, mcd, mcd_mat = dtw.dtw_org_to_trg(org, trg, sdim, ldim, shiftm, startm, endm)

Variable desc.
----

**dtw_org**: result of dtw-ed source features: trg_frame x dim

**twf_mat**: time warping data used to compute dtw_org: trg_frame x 2

**mcd**: average mcd over all frames: double

**mcd_mat**: mcd per frame: trg_frame x 1

**org**: source features data: org_frame x dim

**trg**: target features data: trg_frame x dim

**sdim**: starting dimension to compute distance, e.g., from first dimension would be 0: default=0

**ldim**: last dimension to compute distance, e.g., for 35 dimensional features, and want to compute until the last dimension, then ldim is 34: default=until last dimension

**shiftm**: frame shift value in miliseconds (ms), e.g, 5.0: default = 5.0

**startm**: how many ms would be regarded as starting point, e.g., 0.0 means starting point would only be the 1st frame, if it is 10.0, then possible starting points are the 1st and 2nd, because 10.0/5.0 = 2, and if it is 5.0, it would be also 1st frame only: default = 0.0

**endm**: how many ms would be regarded as ending point, e.g., if 0.0 means ending point would only be the last frame: default = 0.0

Usage of calc_mcd
----

from dtw_c import dtw_c as dtw

mcd, mcd_mat = dtw.calc_mcd(trg_mat, src_mat)

Variable desc.
----

**mcd**: average mcd over all frames: double

**mcd_mat**: mcd per frame: n_frm x 1

**trg_mat**: target feature vector sequence: n_frm x dim

**src_mat**: source feature vector sequence: n_frm x dim

To-do:
----

- other distance measures, e.g., rmse
- function for performing dtw of both source and target sides
- docs
