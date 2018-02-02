#include "dtw_sub_c.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) < (Y) ? (X) : (Y))

void c_calc_distmat(const double* const * y, int row_y, const double* const * x, int row_x, int col, double** distmat) {
	int i, j, k;
	double sumcdk;
	double log_fact = 10.0 / log(10.0);

	// calculate mcd in each frame-pair as distance
	for (i = 0; i < row_y; i++) {
		for (j = 0; j < row_x; j++) {
			for (k = 0, sumcdk = 0.0; k < col; k++) {
				sumcdk += pow(y[i][k]-x[j][k],2);
			} 
			distmat[i][j] = log_fact * sqrt(2.0 * sumcdk);
		}
	}

	return;
}

double c_calc_mcd(const double* const * y, const double* const * x, int row, int col, double* mcdmat) {
	int i, j;
	double sumcdk, summcd;
	double log_fact = 10.0 / log(10.0);

	// calculate mcd in each frame-pair as distance
	for (i = 0, summcd = 0.0; i < row; i++) {
		for (j = 0, sumcdk = 0.0; j < col; j++) {
			sumcdk += pow(y[i][j]-x[i][j],2);
		} 
		mcdmat[i] = log_fact * sqrt(2.0 * sumcdk);
		summcd += mcdmat[i];
	}

	return summcd/(double)row;
}

void c_calc_sumdistmat_asym(const double* const * distmat, int startl, int tar_row, int org_row, double** sumdistmat) {
	int ci, ri, cie;
	double dist, diadist = -1.0, underdist = -1.0, leftdist = -1.0;

	// init first row for sumdistmat and pathmats1 and patmaths2
	for (ci = 0; ci < org_row; ci++) {
		if (ci < startl) {
			sumdistmat[0][ci] = distmat[0][ci];
		} else {
			sumdistmat[0][ci] = -1.0;	
		}
	}

	// iterate from 1st row of target
	for (ri = 1; ri < tar_row; ri++) {
		cie = MAX(startl + 2 * ri, org_row);
		// iterate for every row of source
		for (ci = 0; ci < org_row; ci++) {
			if (ci < cie) {
				// diagonal
				if (ci-1 >= 0) {
					if (sumdistmat[ri-1][ci-1] != -1.0)
						diadist = sumdistmat[ri-1][ci-1] + distmat[ri][ci];
					else diadist = -1.0;
				} else diadist = -1.0;
				// under
				if (sumdistmat[ri-1][ci] != -1.0)
					underdist = sumdistmat[ri-1][ci] + distmat[ri][ci];
				else underdist = -1.0;
				// left
				if (ci-2 >= 0) {
					if (sumdistmat[ri-1][ci-2] != -1.0)
						leftdist = sumdistmat[ri-1][ci-2] + distmat[ri][ci];
					else leftdist = -1.0;
				} else leftdist = -1.0;
				if (diadist != -1.0) {
					dist = diadist;
					if (underdist != -1.0)
						if (underdist < dist)
							dist = underdist;
					if (leftdist != -1.0)
						if (leftdist < dist)
							dist = leftdist;
				} else if (underdist != -1.0) {
					dist = underdist;
					if (leftdist != -1.0)
						if (leftdist < dist)
							dist = leftdist;
				} else dist = leftdist;
				if (dist == -1.0) {
					fprintf(stderr, "error: calc_sumdistmat_asym: dist = -1.0 at [%d][%d]\n", ri, ci);
					exit(1);
				}
				sumdistmat[ri][ci] = dist;
			} else {
				sumdistmat[ri][ci] = -1.0;
			}
		}
	}

	return;
}

void c_calc_pathmats_asym(const double* const * sumdistmat, int startl, int tar_row, int org_row, int** pathmat1, int** pathmat2) {
	int ci, ri, cie;
	double dist, diadist = -1.0, underdist = -1.0, leftdist = -1.0;

	// init first row for sumdistmat and pathmats1 and patmaths2
	for (ci = 0; ci < org_row; ci++) {
		if (ci < startl) {
			pathmat1[0][ci] = 0;
			pathmat2[0][ci] = -ci;
		} else {
			pathmat1[0][ci] = -1;
			pathmat2[0][ci] = -ci;
		}
	}

	// iterate from 1st row of target
	for (ri = 1; ri < tar_row; ri++) {
		cie = MAX(startl + 2 * ri, org_row);
		// iterate for every row of source
		for (ci = 0; ci < org_row; ci++) {
			if (ci < cie) {
				// diagonal
				if (ci-1 >= 0) {
					if (sumdistmat[ri-1][ci-1] != -1.0)
						diadist = sumdistmat[ri-1][ci-1];
					else diadist = -1.0;
				} else diadist = -1.0;
				// under
				if (sumdistmat[ri-1][ci] != -1.0)
					underdist = sumdistmat[ri-1][ci];
				else underdist = -1.0;
				// left
				if (ci-2 >= 0) {
					if (sumdistmat[ri-1][ci-2] != -1.0)
						leftdist = sumdistmat[ri-1][ci-2];
					else leftdist = -1.0;
				} else leftdist = -1.0;
				// select the best path (diagonal, under, or left) w/ minimum dist
				if (diadist != -1.0) {
					dist = diadist;
					if (underdist != -1.0)
						if (underdist < dist)
							dist = underdist;
					if (leftdist != -1.0)
						if (leftdist < dist)
							dist = leftdist;
				} else if (underdist != -1.0) {
					dist = underdist;
					if (leftdist != -1.0)
						if (leftdist < dist)
							dist = leftdist;
				} else dist = leftdist;
				if (dist == -1.0) {
					fprintf(stderr, "error: calc_pathmats_asym: dist = -1.0 at [%d][%d]\n", ri, ci);
					exit(1);
				}
				// set backtrack point for the best path
				if (dist == diadist) {
					pathmat1[ri][ci] = 1;
					pathmat2[ri][ci] = pathmat2[ri-1][ci-1];
				} else if (dist == underdist) {
					pathmat1[ri][ci] = 2;
					pathmat2[ri][ci] = pathmat2[ri-1][ci];
				} else if (dist == leftdist) {
					pathmat1[ri][ci] = 3;
					pathmat2[ri][ci] = pathmat2[ri-1][ci-3];
				} else {
					fprintf(stderr, "error: calc_pathmats_asym: dist=%lf, diadist=%lf, underdist=%lf, leftdist=%lf\n", dist, diadist, underdist, leftdist);
					exit(1);
				}
			} else {
				pathmat1[ri][ci] = -1;
				pathmat2[ri][ci] = -1;
			}
		}
	}

	return;
}

void c_calc_bestpath_asym(const double* const * sumdistmat, const int* const * pathmat, int endl, int tar_row, int org_row, int* pathres) {
	int ri, ci, ris, cis, org_row_e;
	double meandist;

	ris = tar_row-1;
	org_row_e = org_row-endl-1;
	for (cis = org_row_e; cis >= 0; cis--) {
		if (pathmat[ris][cis] < 0) {
			ri = 0;
			ci = -1 * pathmat[ris][cis];
		} else {
			ri = ci = 0;
		}

		meandist = sumdistmat[ris][cis] / (ris - ri + cis - ci + 2);

		if (sumdistmat[ris][cis] < 0.0)
			meandist = 10E16;

		if (cis == org_row_e) {
			pathres[0] = meandist;
			pathres[1] = cis;
		} else {
			pathres[0] = MIN(pathres[0], meandist);
			if (pathres[0] == meandist)
				pathres[1] = cis;
		}
	}

	return;
}

void c_calc_twfunc_asym(const double* const * sumdistmat, const int* const * pathmat, int ci, int tar_row, int org_row, int** twfunc) {
	int ri, cis, tar_row_e = tar_row-1;

	for (ri = tar_row_e, cis = ci; ri > 0; ri--) {
		twfunc[ri][0] = ci;
		twfunc[ri][1] = 1;
		if (pathmat[ri][ci] == 2) {
			// under
		} else if (pathmat[ri][ci] == 3) {
			// left
			ci -= 2;
		} else if (pathmat[ri][ci] == 1) {
			// diagonal
			ci--;
		} else {
			if (pathmat[ri][ci] < 1) {
				fprintf(stderr, "error: c_calc_twfunc_asym: pathmat[%d][%d] = %d < 1\n", ri, ci, pathmat[ri][ci]);
				exit(1);
			}
		}
	}
	twfunc[ri][0] = ci;
	twfunc[ri][1] = 1;

	//fprintf(stderr, "normalized distance %lf\n", sumdistmat[tar_row_e][cis]/(double)(tar_row));

	return;	
}

void c_calc_dtwmat(const double* const * x, const int* const * twfunc, int tar_row, int org_row, int col, double** dtwmat) {
	int ri, ci, dtwri;

	for (ri = 0; ri < tar_row; ri++)
		if (twfunc[ri][0] >= 0) break;

	for(; ri < tar_row; ri++) {
		if (twfunc[ri][0] >= 0)
			dtwri = twfunc[ri][0] + twfunc[ri][1] - 1;
		if (dtwri >= org_row) {
			fprintf(stderr, "error: c_calc_dtw_mat: dtwri=%d >= org_row=%d at [%d], twfunc[%d][0]=%d, twfunc[%d][1]\n",
						dtwri, org_row, ri, ri, twfunc[ri][0], twfunc[ri][1]);
			exit(1);
		}
		if (twfunc[ri][1] == 1) {
			for (ci = 0; ci < col; ci++)
				dtwmat[ri][ci] = x[dtwri][ci];
		} else {
			for (ci = 0; ci < col; ci++)
				dtwmat[ri][ci] = (x[dtwri][ci] + x[twfunc[ri][0]][ci]) / 2;
		}
	}

	return;
}
