import cython
import numpy as np
cimport numpy as np

cdef extern from "dtw_sub_c.h":
	void c_calc_distmat(const double* const * y, int row_y, const double* const * x, int row_x, int col, double** distmat)
	double c_calc_mcd(const double* const * y, const double* const * x, int row, int col, double* mcdmat)
	void c_calc_sumdistmat_asym(const double* const * distmat, int startl, int tar_row, int org_row, double** sumdistmat)
	void c_calc_pathmats_asym(const double* const * sumdistmat, int startl, int tar_row, int org_row, int** pathmat1, int** pathmat2)
	void c_calc_bestpath_asym(const double* const * sumdistmat, const int* const * pathmat, int endl, int tar_row, int org_row, int* pathres)
	void c_calc_twfunc_asym(const double* const * sumdistmat, const int* const * pathmat, int ci, int tar_row, int org_row, int** twfunc)
	void c_calc_dtwmat(const double* const * x, const int* const * twfunc, int tar_row, int org_row, int col, double** dtwmat)

def calc_distmat(np.ndarray[double, ndim=2, mode="c"] tar_mat not None, np.ndarray[double, ndim=2, mode="c"] org_mat not None):
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

	c_calc_distmat(cpp_tar_mat, tar_row, cpp_org_mat, org_row, col, cpp_distmat)

	return np.array(distmat_data, dtype=np.float64)

def calc_sumdistmat_asym(np.ndarray[double, ndim=2, mode="c"] distmat not None, int startl):
	cdef int tar_row, org_row

	tar_row, org_row = distmat.shape[0], distmat.shape[1]

	cdef double[:, ::1] distmat_data = distmat
	cdef double[:, ::1] sumdistmat_data = np.zeros((tar_row, org_row), dtype=np.dtype('float64'))

	cdef np.intp_t[:] tmp1 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp2 = np.zeros(tar_row, dtype=np.intp)

	cdef double** cpp_distmat = <double**> (<void*> &tmp1[0])
	cdef double** cpp_sumdistmat = <double**> (<void*> &tmp2[0])

	cdef np.intp_t i

	for i in range(tar_row):
		cpp_distmat[i] = &distmat_data[i, 0]
		cpp_sumdistmat[i] = &sumdistmat_data[i, 0]

	c_calc_sumdistmat_asym(cpp_distmat, startl, tar_row, org_row, cpp_sumdistmat)

	return np.array(sumdistmat_data, dtype=np.float64)

def calc_pathmats_asym(np.ndarray[double, ndim=2, mode="c"] sumdistmat not None, int startl):
	#print("pathmats_init")
	cdef int tar_row, org_row

	tar_row, org_row = sumdistmat.shape[0], sumdistmat.shape[1]

	cdef double[:, ::1] sumdistmat_data = sumdistmat
	cdef int[:, ::1] pathmat1_data = np.zeros((tar_row, org_row), dtype=np.dtype('int32'))
	cdef int[:, ::1] pathmat2_data = np.zeros((tar_row, org_row), dtype=np.dtype('int32'))

	cdef np.intp_t[:] tmp1 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp2 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp3 = np.zeros(tar_row, dtype=np.intp)

	cdef double** cpp_sumdistmat = <double**> (<void*> &tmp1[0])
	cdef int** cpp_pathmat1 = <int**> (<void*> &tmp2[0])
	cdef int** cpp_pathmat2 = <int**> (<void*> &tmp3[0])

	#print("pathmats_init2")
	cdef np.intp_t i, j

	for i in range(tar_row):
		cpp_sumdistmat[i] = &sumdistmat_data[i, 0]
		cpp_pathmat1[i] = &pathmat1_data[i, 0]
		cpp_pathmat2[i] = &pathmat2_data[i, 0]

	#print("pathmats_start")
	c_calc_pathmats_asym(cpp_sumdistmat, startl, tar_row, org_row, cpp_pathmat1, cpp_pathmat2)
	#print("pathmats_end")

	cdef int[:, :, ::1] pathmats_data = np.zeros((2, tar_row, org_row), dtype=np.dtype('int32'))
	pathmats_data[0] = pathmat1_data
	pathmats_data[1] = pathmat1_data

	return np.array(pathmats_data, dtype=np.int32)

def calc_bestpath_asym(np.ndarray[double, ndim=2, mode="c"] sumdistmat not None, np.ndarray[int, ndim=2, mode="c"] pathmat not None, int endl):
	cdef int tar_row, org_row

	tar_row, org_row = sumdistmat.shape[0], sumdistmat.shape[1]

	cdef double[:, ::1] sumdistmat_data = sumdistmat
	cdef int[:, ::1] pathmat_data = pathmat
	cdef int[::1] pathres_data = np.zeros(2, dtype=np.dtype('int32'))

	cdef np.intp_t[:] tmp1 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp2 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp3 = np.zeros(2, dtype=np.intp)

	cdef double** cpp_sumdistmat = <double**> (<void*> &tmp1[0])
	cdef int** cpp_pathmat = <int**> (<void*> &tmp2[0])
	cdef int* cpp_pathres = <int*> (<void*> &tmp3[0])

	cdef np.intp_t i, j

	for i in range(tar_row):
		cpp_sumdistmat[i] = &sumdistmat_data[i, 0]
		cpp_pathmat[i] = &pathmat_data[i, 0]

	cpp_pathres = &pathres_data[0]

	c_calc_bestpath_asym(cpp_sumdistmat, cpp_pathmat, endl, tar_row, org_row, cpp_pathres)

	return np.array(pathres_data, dtype=np.int32)

def calc_twfunc_asym(np.ndarray[double, ndim=2, mode="c"] sumdistmat not None, np.ndarray[int, ndim=2, mode="c"] pathmat not None, int ci):
	cdef int tar_row, org_row

	tar_row, org_row = sumdistmat.shape[0], sumdistmat.shape[1]

	cdef double[:, ::1] sumdistmat_data = sumdistmat
	cdef int[:, ::1] pathmat_data = pathmat
	cdef int[:, ::1] twfunc_data = np.zeros((tar_row, 2), dtype=np.dtype('int32'))

	cdef np.intp_t[:] tmp1 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp2 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp3 = np.zeros(tar_row, dtype=np.intp)

	cdef double** cpp_sumdistmat = <double**> (<void*> &tmp1[0])
	cdef int** cpp_pathmat = <int**> (<void*> &tmp2[0])
	cdef int** cpp_twfunc = <int**> (<void*> &tmp3[0])

	cdef np.intp_t i

	for i in range(tar_row):
		cpp_sumdistmat[i] = &sumdistmat_data[i, 0]
		cpp_pathmat[i] = &pathmat_data[i, 0]
		cpp_twfunc[i] = &twfunc_data[i, 0]

	c_calc_twfunc_asym(cpp_sumdistmat, cpp_pathmat, ci, tar_row, org_row, cpp_twfunc)

	return np.array(twfunc_data, dtype=np.int32)

def dtw_body_asym(np.ndarray[double, ndim=2, mode="c"] distmat not None, int shiftm, int startm, int endm):
	cdef int tar_row, org_row, ci
	cdef double mindist

	tar_row, org_row = distmat.shape[0], distmat.shape[1]
	startl = (int)(startm/shiftm)
	endl = (int)(endm/shiftm)
	if ((startl + endl) > tar_row or (startl + endl) > org_row):
		print("error: dtw_body_asym: startl=%d, endl=%d, tar_row=%d, org_row=%d" % (startl, endl, tar_row, org_row))
		exit()
	startl = np.maximum(startl, 1)	

	# calc acc dist mat
	sumdistmat = calc_sumdistmat_asym(distmat, startl)
	# calc backtrack path mat
	pathmats = calc_pathmats_asym(sumdistmat, startl)
	# calc best path and dist
	pathres = calc_bestpath_asym(sumdistmat, pathmats[1], endl)
	mindist = pathres[0]
	if (mindist == 10E16):
		print("error: dtw_body_asym: can't reach and end range")
		print("error: dtw_body_asym: mindist=%lf, shiftm=%d, startm=%d, endm=%d" % (mindist, shiftm, startm, endm))
		exit()
	ci = pathres[1]
	print("#sum_distance [%d][%d]: %lf" % (tar_row-1, ci, sumdistmat[tar_row-1][ci]))
	# calc twf function
	twfunc = calc_twfunc_asym(sumdistmat, pathmats[0], ci)

	return twfunc 

def calc_dtwmat(np.ndarray[double, ndim=2, mode="c"] org_mat not None, np.ndarray[int, ndim=2, mode="c"] twf_mat not None, int tar_row):
	cdef int org_row, col

	org_row, col = org_mat.shape[0], org_mat.shape[1]

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

def dtw_org_to_trg(np.ndarray[double, ndim=2, mode="c"] tar_mat not None, np.ndarray[double, ndim=2, mode="c"] org_mat not None, int sdim, int ldim, int shiftm, int startm, int endm):
	cdef int tar_row, col
	
	tar_row, col = tar_mat.shape[0], tar_mat.shape[1]
	ldim += 1

	dist_mat = calc_distmat(np.array(tar_mat[:,sdim:ldim], dtype=np.float64), np.array(org_mat[:,sdim:ldim], dtype=np.float64))
	twf_func = dtw_body_asym(dist_mat, shiftm, startm, endm)
	dtw_mat = calc_dtwmat(org_mat, twf_func, tar_row)

	cdef double[:, ::1] tarmat_data = tar_mat
	cdef double[:, ::1] dtwmat_data = dtw_mat
	cdef double[::1] mcdmat_data = np.zeros(tar_row, dtype=np.dtype('float64'))

	cdef np.intp_t[:] tmp1 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp2 = np.zeros(tar_row, dtype=np.intp)
	cdef np.intp_t[:] tmp3 = np.zeros(tar_row, dtype=np.intp)

	cdef double** cpp_tarmat = <double**> (<void*> &tmp1[0])
	cdef double** cpp_dtwmat = <double**> (<void*> &tmp2[0])
	cdef double* cpp_mcdmat = <double*> (<void*> &tmp3[0])

	cdef np.intp_t i

	for i in range(tar_row):
		cpp_tarmat[i] = &tarmat_data[i, 0]
		cpp_dtwmat[i] = &dtwmat_data[i, 0]
	
	cpp_mcdmat = &mcdmat_data[0]
	
	mcd = c_calc_mcd(cpp_tarmat, cpp_dtwmat, tar_row, col, cpp_mcdmat)

	return dtw_mat, twf_func, mcd, np.array(mcdmat_data, dtype=np.float64)
