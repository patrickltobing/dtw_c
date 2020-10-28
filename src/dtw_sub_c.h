// By: Patrick Lumban Tobing (Nagoya University, August 2018 - October 2020)

void c_calc_distmat(const double* const * y, int row_y, const double* const * x, int row_x, int col, double** distmat, int mcd);
double c_calc_mcd(const double* const * y, const double* const * x, int row, int col, double* mcdmat, int mcd);
void c_calc_twfunc_asym(const double* const * distmat, int winl, int tar_row, int org_row, double** sumdistmat, int** pathmat1, int** pathmat2, int** twfunc);
void c_calc_dtwmat(const double* const * x, const int* const * twfunc, int tar_row, int org_row, int col, double** dtwmat);
