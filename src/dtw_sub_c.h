void c_calc_distmat(const double* const * y, int row_y, const double* const * x, int row_x, int col, double** distmat);
double c_calc_mcd(const double* const * y, const double* const * x, int row, int col, double* mcdmat);
void c_calc_sumdistmat_asym(const double* const * distmat, int startl, int tar_row, int org_row, double** sumdistmat);
void c_calc_pathmats_asym(const double* const * sumdistmat, int startl, int tar_row, int org_row, int** pathmat1, int** pathmat2);
void c_calc_bestpath_asym(const double* const * sumdistmat, const int* const * pathmat, int endl, int tar_row, int org_row, int* pathres);
void c_calc_twfunc_asym(const double* const * sumdistmat, const int* const * pathmat, int ci, int tar_row, int org_row, int** twfunc);
void c_calc_dtwmat(const double* const * x, const int* const * twfunc, int tar_row, int org_row, int col, double** dtwmat);
