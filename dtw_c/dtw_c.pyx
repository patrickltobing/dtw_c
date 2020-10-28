# By: Patrick Lumban Tobing (Nagoya University, August 2018 - October 2020)

## Based on WORLD codec implementation by Jeremy Hsu
## https://github.com/JeremyCCHsu/Python-Wrapper-for-World-Vocoder

import cython
import numpy as np
cimport numpy as np

cdef extern from "dtw_sub_c.h":
    void c_calc_distmat(const double* const * y, int row_y, const double* const * x, int row_x, int col, double** distmat, int mcd)
    double c_calc_mcd(const double* const * y, const double* const * x, int row, int col, double* mcdmat, int mcd)
    void c_calc_twfunc_asym(const double* const * distmat, int winl, int tar_row, int org_row, double** sumdistmat, int** pathmat1, int** pathmat2, int** twfunc)
    void c_calc_dtwmat(const double* const * x, const int* const * twfunc, int tar_row, int org_row, int col, double** dtwmat)
    
def calc_distmat(np.ndarray[double, ndim=2, mode="c"] tar_mat not None, np.ndarray[double, ndim=2, mode="c"] org_mat not None, int mcd=1):
    cdef int tar_row, org_row, col
    
    tar_row, org_row, col = tar_mat.shape[0], org_mat.shape[0], tar_mat.shape[1]
    
    cdef double[:, ::1] tar_data = tar_mat
    cdef double[:, ::1] org_data = org_mat
    cdef double[:, ::1] distmat_data = np.zeros((tar_row, org_row), dtype=np.dtype('float64'))
    
    cdef np.intp_t[:] tmp1 = np.zeros(tar_row, dtype=np.intp)
    cdef np.intp_t[:] tmp2 = np.zeros(org_row, dtype=np.intp)
    cdef np.intp_t[:] tmp3 = np.zeros(tar_row, dtype=np.intp)
    
    cdef double** cpp_tar_mat = <double**> (<void*> &tmp1[0])
    cdef double** cpp_org_mat = <double**> (<void*> &tmp2[0])
    cdef double** cpp_distmat = <double**> (<void*> &tmp3[0])
    
    cdef np.intp_t i
    
    for i in range(tar_row):
        cpp_tar_mat[i] = &tar_data[i, 0]
        cpp_distmat[i] = &distmat_data[i, 0]
    
    for i in range(org_row):
        cpp_org_mat[i] = &org_data[i, 0]
    
    c_calc_distmat(cpp_tar_mat, tar_row, cpp_org_mat, org_row, col, cpp_distmat, mcd)
    
    return np.array(distmat_data, dtype=np.float64)

def calc_twfunc_asym(np.ndarray[double, ndim=2, mode="c"] distmat not None, int winl):
    cdef int tar_row, org_row
    
    tar_row, org_row = distmat.shape[0], distmat.shape[1]
    
    cdef double[:, ::1] distmat_data = distmat
    cdef double[:, ::1] sumdistmat_data = np.zeros((tar_row, org_row), dtype=np.dtype('float64'))
    cdef int[:, ::1] pathmat1_data = np.zeros((tar_row, org_row), dtype=np.dtype('int32'))
    cdef int[:, ::1] pathmat2_data = np.zeros((tar_row, org_row), dtype=np.dtype('int32'))
    cdef int[:, ::1] twfunc_data = np.zeros((tar_row, 2), dtype=np.dtype('int32'))
    
    cdef np.intp_t[:] tmp1 = np.zeros(tar_row, dtype=np.intp)
    cdef np.intp_t[:] tmp2 = np.zeros(tar_row, dtype=np.intp)
    cdef np.intp_t[:] tmp3 = np.zeros(tar_row, dtype=np.intp)
    cdef np.intp_t[:] tmp4 = np.zeros(tar_row, dtype=np.intp)
    cdef np.intp_t[:] tmp5 = np.zeros(tar_row, dtype=np.intp)
    
    cdef double** cpp_distmat = <double**> (<void*> &tmp1[0])
    cdef double** cpp_sumdistmat = <double**> (<void*> &tmp2[0])
    cdef int** cpp_pathmat1 = <int**> (<void*> &tmp3[0])
    cdef int** cpp_pathmat2 = <int**> (<void*> &tmp4[0])
    cdef int** cpp_twfunc = <int**> (<void*> &tmp5[0])
    
    cdef np.intp_t i
    
    for i in range(tar_row):
        cpp_distmat[i] = &distmat_data[i, 0]
        cpp_sumdistmat[i] = &sumdistmat_data[i, 0]
        cpp_pathmat1[i] = &pathmat1_data[i, 0]
        cpp_pathmat2[i] = &pathmat2_data[i, 0]
        cpp_twfunc[i] = &twfunc_data[i, 0]
    
    c_calc_twfunc_asym(cpp_distmat, winl, tar_row, org_row, cpp_sumdistmat, cpp_pathmat1, cpp_pathmat2, cpp_twfunc)
    
    cdef int[:, :, ::1] pathmats_data = np.zeros((2, tar_row, org_row), dtype=np.dtype('int32'))
    pathmats_data[0] = pathmat1_data
    pathmats_data[1] = pathmat2_data
    
    return np.array(twfunc_data, dtype=np.int32)

def calc_dtwmat(np.ndarray[double, ndim=2, mode="c"] org_mat not None, np.ndarray[int, ndim=2, mode="c"] twf_mat not None):
    cdef int org_row, tar_row, col
    
    org_mat = np.array(org_mat, dtype=np.float64)
    twf_mat = np.array(twf_mat, dtype=np.int32)
    
    org_row, tar_row, col = org_mat.shape[0], twf_mat.shape[0], org_mat.shape[1]
    
    cdef double[:, ::1] orgmat_data = org_mat
    cdef int[:, ::1] twfmat_data = twf_mat
    cdef double[:, ::1] dtwmat_data = np.zeros((tar_row, col), dtype=np.dtype('float64'))
    
    cdef np.intp_t[:] tmp1 = np.zeros(org_row, dtype=np.intp)
    cdef np.intp_t[:] tmp2 = np.zeros(tar_row, dtype=np.intp)
    cdef np.intp_t[:] tmp3 = np.zeros(tar_row, dtype=np.intp)
    
    cdef double** cpp_orgmat = <double**> (<void*> &tmp1[0])
    cdef int** cpp_twfmat = <int**> (<void*> &tmp2[0])
    cdef double** cpp_dtwmat = <double**> (<void*> &tmp3[0])
    
    cdef np.intp_t i
    
    for i in range(tar_row):
        cpp_twfmat[i] = &twfmat_data[i, 0]
        cpp_dtwmat[i] = &dtwmat_data[i, 0]
    
    for i in range(org_row):
        cpp_orgmat[i] = &orgmat_data[i, 0]
    
    c_calc_dtwmat(cpp_orgmat, cpp_twfmat, tar_row, org_row, col, cpp_dtwmat)
    
    return np.array(dtwmat_data, dtype=np.float64)

def calc_mcd(np.ndarray[double, ndim=2, mode="c"] tar_mat not None, np.ndarray[double, ndim=2, mode="c"] org_mat not None, int mcd=1):
    assert org_mat.ndim == 2, "org.ndim = %d != 2" % org_mat.ndim
    assert tar_mat.ndim == 2, "trg.ndim = %d != 2" % tar_mat.ndim
    assert org_mat.shape[0] == tar_mat.shape[0], "org.shape[0] = %d != %d = trg.shape[0]" % (org_mat.shape[0], tar_mat.shape[0])
    assert org_mat.shape[1] == tar_mat.shape[1], "org.shape[1] = %d != %d = trg.shape[1]" % (org_mat.shape[1], tar_mat.shape[1])
    
    cdef int row, col
    
    if mcd < 0:
        tar_mat = np.array(np.clip(tar_mat, a_min=1e-16, a_max=None), dtype=np.float64)
        org_mat = np.array(np.clip(org_mat, a_min=1e-16, a_max=None), dtype=np.float64)
    else:
        tar_mat = np.array(tar_mat, dtype=np.float64)
        org_mat = np.array(org_mat, dtype=np.float64)
    
    row, col = tar_mat.shape[0], tar_mat.shape[1]
    
    cdef double[:, ::1] tarmat_data = tar_mat
    cdef double[:, ::1] orgmat_data = org_mat
    cdef double[::1] mcdmat_data = np.zeros(row, dtype=np.dtype('float64'))
    
    cdef np.intp_t[:] tmp1 = np.zeros(row, dtype=np.intp)
    cdef np.intp_t[:] tmp2 = np.zeros(row, dtype=np.intp)
    cdef np.intp_t[:] tmp3 = np.zeros(row, dtype=np.intp)
    
    cdef double** cpp_tarmat = <double**> (<void*> &tmp1[0])
    cdef double** cpp_orgmat = <double**> (<void*> &tmp2[0])
    cdef double* cpp_mcdmat = <double*> (<void*> &tmp3[0])
    
    cdef np.intp_t i
    
    for i in range(row):
        cpp_tarmat[i] = &tarmat_data[i, 0]
        cpp_orgmat[i] = &orgmat_data[i, 0]
    
    cpp_mcdmat = &mcdmat_data[0]
    
    mcd_res = c_calc_mcd(cpp_tarmat, cpp_orgmat, row, col, cpp_mcdmat, mcd)
    
    return mcd_res, np.array(mcdmat_data, dtype=np.float64)

def dtw_org_to_trg(np.ndarray[double, ndim=2, mode="c"] org_mat not None, np.ndarray[double, ndim=2, mode="c"] tar_mat not None, int sdim=0, int ldim=-1, double shiftm=5.0, double winlenm=100.0, int mcd=1):
    assert org_mat.ndim == 2, "org.ndim = %d != 2" % org_mat.ndim
    assert tar_mat.ndim == 2, "trg.ndim = %d != 2" % tar_mat.ndim
    assert org_mat.shape[1] == tar_mat.shape[1], "org.shape[1] = %d != %d = trg.shape[1]" % (org_mat.shape[1], tar_mat.shape[1])
    if ldim < 0:
        ldim = tar_mat.shape[1]-1
    assert sdim >= 0 and ldim >= 0 and sdim <= ldim, "sdim = %d, ldim = %d, sdim and ldim should be >= 0, sdim should be <= ldim" % (sdim, ldim)
    assert sdim < org_mat.shape[1], "sdim = %d >= %d = org.shape[1]" % (sdim, org_mat.shape[1])
    assert ldim < org_mat.shape[1], "ldim = %d >= %d = org.shape[1]" % (ldim, org_mat.shape[1])
    assert sdim < tar_mat.shape[1], "sdim = %d >= %d = trg.shape[1]" % (sdim, tar_mat.shape[1])
    assert ldim < tar_mat.shape[1], "ldim = %d >= %d = trg.shape[1]" % (ldim, tar_mat.shape[1])
    assert winlenm >= 0.0, "shiftm = %d < 0.0" % winlenm
    winl = (int)(winlenm/shiftm)
    if winl+1 > org_mat.shape[0]:
        winl = org_mat.shape[1]-1
    
    cdef int tar_row = tar_mat.shape[0]
    ldim += 1
    
    if mcd < 0:
        tar_mat = np.array(np.clip(tar_mat, a_min=1e-16, a_max=None), dtype=np.float64)
        org_mat = np.array(np.clip(org_mat, a_min=1e-16, a_max=None), dtype=np.float64)
    else:
        tar_mat = np.array(tar_mat, dtype=np.float64)
        org_mat = np.array(org_mat, dtype=np.float64)
    dist_mat = calc_distmat(np.array(tar_mat[:,sdim:ldim], dtype=np.float64), np.array(org_mat[:,sdim:ldim], dtype=np.float64), mcd)
    twf_func = calc_twfunc_asym(dist_mat, winl)
    dtw_mat = calc_dtwmat(org_mat, twf_func)
    mcd_res, mcd_mat = calc_mcd(tar_mat, dtw_mat, mcd)
    
    return dtw_mat, twf_func, mcd_res, mcd_mat
