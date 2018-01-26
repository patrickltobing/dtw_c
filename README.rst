mic Time Warping (DTW) implementation in C for Python

Perform DTW of source features (org) into target features (trg) with mel-cepstral distortion (mcd) distance

<<<<<<< HEAD
Install
----

pip install dtw_c

=======
Perform DTW of source features (org) into target features (trg) with mel-cepstral distortion (mcd) distance

Install
----

pip install dtw_c

>>>>>>> 602cf27a99a0cb4f7edad05b5ca2956c5d7b7964
Usage
----

from dtw_c import dtw_c as dtw

dtw_org, twf_mat, mcd, mcd_mat = dtw.dtw_org_to_trg(trg, org, sdim, ldim, shiftm, startm, endm)

Variable desc.
----
<<<<<<< HEAD

**dtw_org**: result of dtw-ed source features

**twf_mat**: time warping data used to compute dtw_org

**mcd**: average mcd over all frames

**mcd_mat**: mcd per frame

**org**: source features data, org_frame x dim

**trg**: target features data, trg_frame x dim

**sdim**: starting dimension to compute distance, e.g., from first dimension would be 0, default=0

**ldim**: last dimension to compute distance, e.g., for 35 dimensional features, and want to compute until the last dimension, then ldim is 34, default=24

**shiftm**: frame shift value in miliseconds (ms), e.g, 5.0, default = 5.0

**startm**: how many ms would be regarded as starting point, e.g., 0.0 means starting point would only be the 1st frame, if it is 10.0, then possible starting points are the 1st and 2nd, because 10.0/5.0 = 2, and if it is 5.0, it would be also 1st frame only, default = 0.0

**endm**: how many ms would be regarded as ending point, e.g., if 0.0 means ending point would only be the last frame, default = 0.0
=======

**dtw_org**: result of dtw-ed source features

**twf_mat**: time warping data used to compute dtw_org

**mcd**: average mcd over all frames

**mcd_mat**: mcd per frame

**sdim**: starting dimension to compute distance, e.g., from first dimension would be 0

**ldim**: last dimension to compute distance, e.g., for 35 dimensional features, and want to compute until the last dimension, then ldim is 34

**shiftm**: frame shift value in miliseconds (ms), e.g, 5.0

**startm**: how many ms would be regarded as starting point, e.g., 0.0 means starting point would only be the 1st frame, if it is 10.0, then possible starting points are the 1st and 2nd, because 10.0/5.0 = 2, and if it is 5.0, it would be also 1st frame only

**endm**: how many ms would be regarded as ending point, e.g., if 0.0 means ending point would only be the last frame
>>>>>>> 602cf27a99a0cb4f7edad05b5ca2956c5d7b7964

To-do:
----

<<<<<<< HEAD
- other distance measures, e.g., rmse
- function for performing dtw of both source and target sides
- docs
=======
- assert for input
- set default values
- other distance measures, e.g., rmse
- function for performing dtw of both source and target sides
>>>>>>> 602cf27a99a0cb4f7edad05b5ca2956c5d7b7964
