// By: Patrick Lumban Tobing (Nagoya University, August 2018 - October 2020)

/* Based on a GMM-based VC implementation by Tomoki Toda (Nara Institute of Science and Technology)
   "Voice conversion based on maximum-likelihood estimation of spectral parameter trajectory, 2007" */

#include "dtw_sub_c.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

void c_calc_distmat(const double* const * y, int row_y, const double* const * x, int row_x, int col, double** distmat, int mcd) {
	int i, j, k;
	double sumcdk, sumprod, sumx, sumy;
	double log_fact = 10.0 / log(10.0);

    if (mcd > 0) { // use mcd distance
	    for (i = 0; i < row_y; i++) {
	    	for (j = 0; j < row_x; j++) {
	    		for (k = 0, sumcdk = 0.0; k < col; k++) {
	    			sumcdk += pow(y[i][k]-x[j][k],2);
	    		}
	    		distmat[i][j] = log_fact * sqrt(2.0 * sumcdk);
	    	}
	    }
    } else if (mcd < 0) { // use log-spectral distance of magnitude
	    for (i = 0; i < row_y; i++) {
	    	for (j = 0; j < row_x; j++) {
	    		for (k = 0, sumcdk = 0.0; k < col; k++) {
	    			sumcdk += pow(20*(log10(y[i][k])-log10(x[j][k])),2);
	    		}
	    		distmat[i][j] = sqrt(sumcdk/k);
	    	}
	    }
    } else { // use reversed cosine similarity distance (because the algorithm takes minimum distance)
	    for (i = 0; i < row_y; i++) {
	    	for (j = 0; j < row_x; j++) {
	    		for (k = 0, sumprod = 0.0, sumx = 0.0, sumy = 0.0; k < col; k++) {
	    			sumprod += x[j][k]*y[i][k];
	    			sumx += pow(x[j][k],2);
	    			sumy += pow(y[i][k],2);
	    		}
	    		distmat[i][j] = (sumprod/(sqrt(sumx)*sqrt(sumy)))*(-1);
	    	}
	    }
    }

	return;
}

double c_calc_mcd(const double* const * y, const double* const * x, int row, int col, double* mcdmat, int mcd) {
	int i, j;
	double sumcdk, summcd, sumprod, sumx, sumy;
	double log_fact = 10.0 / log(10.0);

    if (mcd > 0) { // calculate frame-pair mcd
	    for (i = 0, summcd = 0.0; i < row; i++) {
	    	for (j = 0, sumcdk = 0.0; j < col; j++) {
	    		sumcdk += pow(y[i][j]-x[i][j],2);
	    	} 
	    	mcdmat[i] = log_fact * sqrt(2.0 * sumcdk);
	    	summcd += mcdmat[i];
	    }
    } else if (mcd < 0) { // calculate frame-pair lsd
	    for (i = 0, summcd = 0.0; i < row; i++) {
	    	for (j = 0, sumcdk = 0.0; j < col; j++) {
	    	    sumcdk += pow(20*(log10(y[i][j])-log10(x[i][j])),2);
	    	} 
	    	mcdmat[i] = sqrt(sumcdk/j);
	    	summcd += mcdmat[i];
	    }
    } else { // calculate frame-pair cosine similarity
	    for (i = 0, summcd = 0.0; i < row; i++) {
	    	for (j = 0, sumprod = 0.0, sumx = 0.0, sumy = 0.0; j < col; j++) {
	    		sumprod += x[i][j]*y[i][j];
	    		sumx += pow(x[i][j],2);
	    		sumy += pow(y[i][j],2);
	    	} 
	    	mcdmat[i] = sumprod/(sqrt(sumx)*sqrt(sumy));
	    	summcd += mcdmat[i];
	    }
    }

	return summcd/(double)row;
}

void c_calc_twfunc_asym(const double* const * distmat, int winl, int tar_row, int org_row, double** sumdistmat, int** pathmat1, int** pathmat2, int** twfunc) {
	double dist, diadist = -1.0, underdist = -1.0, leftdist = -1.0;
	double underdiadist = -1.0, dialeftdist = -1.0;
	double factor, bound, meandist, minmeandist = 1E6;
	int ri, ci, cis;
	int i = 1, j = 0, k, winl_e = winl + 1;
	int idx_under, idx_underdia, idx_dia, idx_dialeft, idx_left;
	int last_j, last_left, last_right, left, right, last_j_diff, j_diff, count;

	for (ci = 0; ci < org_row; ci++) { // init first target row
		if (ci < winl_e) {
			sumdistmat[0][ci] = distmat[0][ci];
			pathmat1[0][ci] = 0;
			pathmat2[0][ci] = -ci;
			//fprintf(stderr, "ci:%d;  dist=%lf;  sumdist=%lf\n", ci, distmat[0][ci], sumdistmat[0][ci]);
		} else {
			sumdistmat[0][ci] = -1.0;
			pathmat1[0][ci] = -1;
			pathmat2[0][ci] = -ci;
			//fprintf(stderr, "~ ci:%d;  dist=%lf;  sumdist=%lf\n", ci, distmat[0][ci], sumdistmat[0][ci]);
		}
	}
	for (ri = 1; ri < tar_row; ri++)
		for (ci = 0; ci < org_row; ci++) {
			sumdistmat[ri][ci] = -1.0;
			pathmat1[ri][ci] = -1;
			pathmat2[ri][ci] = -1;
		}
	//fprintf(stderr, "dist_under=%lf\n", sumdistmat[0][0]);
	//fprintf(stderr, "~dist_under=%lf\n", sumdistmat[i-1][0]);

	if (tar_row > org_row) { // target is longer
		factor = (double)tar_row/org_row;
		count = 1;
		bound = factor * count;
		left = j-winl;
		if (left < 0) left = 0;
		right = j+winl_e;
		if (right > org_row) right = org_row;
		for (; i < tar_row; i++) {
			last_j = j;
			last_left = last_j-winl;
			if (last_left < 0) last_left = 0;
			last_right = last_j+winl_e;
			if (last_right > org_row) last_right = org_row;
			if ((float)i >= bound) {
				j = count;
				left = j-winl;
				if (left < 0) left = 0;
				right = j+winl_e;
				if (right > org_row) right = org_row;
				count++;
				bound = factor * count;
			}
			last_j_diff = j-last_j;
			for (k = left; k < right; k++) {
				j_diff = j-k;
				idx_under = last_j+(last_j_diff-j_diff);
				idx_underdia = idx_under-1;
				idx_dia = idx_underdia-1;
				idx_dialeft = idx_dia-1;
				idx_left = idx_dialeft-1;
				// take minimum accumulated distance
				if (idx_under < last_right && idx_under >= last_left) { // record under dist
					underdist = sumdistmat[i-1][idx_under] + distmat[i][k];
				} else underdist = -1;
				if (idx_underdia < last_right && idx_underdia >= last_left) { // record under_dia dist
					underdiadist = sumdistmat[i-1][idx_underdia] + distmat[i][k];
				} else underdiadist = -1;
				if (idx_dia < last_right && idx_dia >= last_left) { // record dia dist
					diadist = sumdistmat[i-1][idx_dia] + distmat[i][k];
				} else diadist = -1;
				if (idx_dialeft < last_right && idx_dialeft >= last_left) { // record dia_left dist
					dialeftdist = sumdistmat[i-1][idx_dialeft] + distmat[i][k];
				} else dialeftdist = -1;
				if (idx_left < last_right && idx_left >= last_left) { // record left dist
					leftdist = sumdistmat[i-1][idx_left] + distmat[i][k];
				} else leftdist = -1;
                if (underdist != -1) { // compare under dist
					dist = underdist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist <= dist)
							dist = dialeftdist;
					if (underdiadist != - 1)
						if (underdiadist <= dist)
							dist = underdiadist;
					if (diadist != - 1)
						if (diadist <= dist)
							dist = diadist;
				} else if (underdiadist != -1) { // compare under_dia dist
					dist = underdiadist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist < dist)
							dist = dialeftdist;
					if (diadist != - 1)
						if (diadist <= dist)
							dist = diadist;
				} else if (diadist != -1) { // compare dia dist
					dist = diadist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist < dist)
							dist = dialeftdist;
				} else if (dialeftdist != -1) { // compare dia_left dist
					dist = dialeftdist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
				} else dist = leftdist;
				if (dist == -1.0) {
					fprintf(stderr, "a error: calc_sumdistmat_asym: dist = -1.0 at [%d][%d]\n", i, k);
					exit(1);
				}
				sumdistmat[i][k] = dist;
				// set backtrack point for the best path
				if (dist == diadist) {
					pathmat1[i][k] = 1;
					pathmat2[i][k] = pathmat2[i-1][k-2];
				} else if (dist == underdiadist) {
					pathmat1[i][k] = 2;
					pathmat2[i][k] = pathmat2[i-1][k-1];
				} else if (dist == dialeftdist) {
					pathmat1[i][k] = 3;
					pathmat2[i][k] = pathmat2[i-1][k-3];
				} else if (dist == underdist) {
					pathmat1[i][k] = 4;
					pathmat2[i][k] = pathmat2[i-1][k];
				} else if (dist == leftdist) {
					pathmat1[i][k] = 5;
					pathmat2[i][k] = pathmat2[i-1][k-4];
				} else {
					fprintf(stderr, "a error: calc_pathmats_asym: dist=%lf, diadist=%lf, underdiadist=%lf, dialeftdist=%lf, underdist=%lf, leftdist=%lf\n", dist, diadist, underdiadist, dialeftdist, underdist, leftdist);
					exit(1);
				}
			}
		}
	} else if (tar_row < org_row) { // source is longer
		factor = (double)org_row/tar_row;
		//fprintf(stderr, "b i-1=%d;  dist_under=%lf\n", i-1, sumdistmat[i-1][0]);
		for (; i < tar_row; i++) {
			last_j = j;
			last_left = last_j-winl;
			if (last_left < 0) last_left = 0;
			last_right = last_j+winl_e;
			if (last_right > org_row) last_right = org_row;
			j = (int)floor(i*factor);
			if (j-last_j > 4) j = last_j+4;
			left = j-winl;
			if (left < 0) left = 0;
			right = j+winl_e;
			if (right > org_row) right = org_row;
			last_j_diff = j-last_j;
			//fprintf(stderr, "b i=%d;  last_j=%d;  j=%d;  last_j_diff=%d;  last_left=%d; last_right=%d;  left=%d;  right=%d\n", i, last_j, j, last_j_diff, last_left, last_right, left, right);
			for (k = left; k < right; k++) {
				j_diff = j-k;
				idx_under = last_j+(last_j_diff-j_diff);
				idx_underdia = idx_under-1;
				idx_dia = idx_underdia-1;
				idx_dialeft = idx_dia-1;
				idx_left = idx_dialeft-1;
				//fprintf(stderr, "b k=%d;  dist=%lf;  j_diff=%d;  idx_under=%d;  idx_underdia=%d;  idx_dia=%d;  idx_dialeft=%d;  idx_left=%d\n", k, distmat[i][k], j_diff, idx_under, idx_underdia, idx_dia, idx_dialeft, idx_left);
				// take minimum accumulated distance
				if (idx_under < last_right && idx_under >= last_left) { // record under dist
					underdist = sumdistmat[i-1][idx_under] + distmat[i][k];
					//fprintf(stderr, "b i-1=%d;  idx_under=%d;  underdist=%lf;  dist_under=%lf\n", i-1, idx_under, underdist, sumdistmat[0][idx_under]);
				} else underdist = -1;
				if (idx_underdia < last_right && idx_underdia >= last_left) { // record under_dia dist
					underdiadist = sumdistmat[i-1][idx_underdia] + distmat[i][k];
					//fprintf(stderr, "b dist_underdia=%lf\n", sumdistmat[i-1][idx_underdia]);
				} else underdiadist = -1;
				if (idx_dia < last_right && idx_dia >= last_left) { // record dia dist
					diadist = sumdistmat[i-1][idx_dia] + distmat[i][k];
					//fprintf(stderr, "b dist_dia=%lf\n", sumdistmat[i-1][idx_dia]);
				} else diadist = -1;
				if (idx_dialeft < last_right && idx_dialeft >= last_left) { // record dia_left dist
					dialeftdist = sumdistmat[i-1][idx_dialeft] + distmat[i][k];
					//fprintf(stderr, "b dist_dialeft=%lf\n", sumdistmat[i-1][idx_dialeft]);
				} else dialeftdist = -1;
				if (idx_left < last_right && idx_left >= last_left) { // record left dist
					leftdist = sumdistmat[i-1][idx_left] + distmat[i][k];
					//fprintf(stderr, "b dist_left=%lf\n", sumdistmat[i-1][idx_left]);
				} else leftdist = -1;
                if (underdist != -1) { // compare under dist
					dist = underdist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist <= dist)
							dist = dialeftdist;
					if (underdiadist != - 1)
						if (underdiadist <= dist)
							dist = underdiadist;
					if (diadist != - 1)
						if (diadist <= dist)
							dist = diadist;
				} else if (underdiadist != -1) { // compare under_dia dist
					dist = underdiadist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist < dist)
							dist = dialeftdist;
					if (diadist != - 1)
						if (diadist <= dist)
							dist = diadist;
				} else if (diadist != -1) { // compare dia dist
					dist = diadist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist < dist)
							dist = dialeftdist;
				} else if (dialeftdist != -1) { // compare dia_left dist
					dist = dialeftdist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
				} else dist = leftdist;
				//fprintf(stderr, "b dist=%lf, diadist=%lf, underdiadist=%lf, dialeftdist=%lf, underdist=%lf, leftdist=%lf\n", dist, diadist, underdiadist, dialeftdist, underdist, leftdist);
				if (dist == -1.0) {
					fprintf(stderr, "b error: calc_sumdistmat_asym: dist = -1.0 at [%d][%d]\n", i, k);
					exit(1);
				}
				sumdistmat[i][k] = dist;
				// set backtrack point for the best path
				if (dist == diadist) {
					pathmat1[i][k] = 1;
					pathmat2[i][k] = pathmat2[i-1][k-2];
				} else if (dist == underdiadist) {
					pathmat1[i][k] = 2;
					pathmat2[i][k] = pathmat2[i-1][k-1];
				} else if (dist == dialeftdist) {
					pathmat1[i][k] = 3;
					pathmat2[i][k] = pathmat2[i-1][k-3];
				} else if (dist == underdist) {
					pathmat1[i][k] = 4;
					pathmat2[i][k] = pathmat2[i-1][k];
				} else if (dist == leftdist) {
					pathmat1[i][k] = 5;
					pathmat2[i][k] = pathmat2[i-1][k-4];
				} else {
					fprintf(stderr, "b error: calc_pathmats_asym: dist=%lf, diadist=%lf, underdiadist=%lf, dialeftdist=%lf, underdist=%lf, leftdist=%lf\n", dist, diadist, underdiadist, dialeftdist, underdist, leftdist);
					exit(1);
				}
			}
		}
	} else { // same length
		for (; i < tar_row; i++) {
			last_j = j;
			last_left = last_j-winl;
			if (last_left < 0) last_left = 0;
			last_right = last_j+winl_e;
			if (last_right > org_row) last_right = org_row;
			j = i;
			left = j-winl;
			if (left < 0) left = 0;
			right = j+winl_e;
			if (right > org_row) right = org_row;
			last_j_diff = j-last_j;
			for (k = left; k < right; k++) {
				j_diff = j-k;
				idx_under = last_j+(last_j_diff-j_diff);
				idx_underdia = idx_under-1;
				idx_dia = idx_underdia-1;
				idx_dialeft = idx_dia-1;
				idx_left = idx_dialeft-1;
				// take minimum accumulated distance
				if (idx_under < last_right && idx_under >= last_left) { // record under dist
					underdist = sumdistmat[i-1][idx_under] + distmat[i][k];
				} else underdist = -1;
				if (idx_underdia < last_right && idx_underdia >= last_left) { // record under_dia dist
					underdiadist = sumdistmat[i-1][idx_underdia] + distmat[i][k];
				} else underdiadist = -1;
				if (idx_dia < last_right && idx_dia >= last_left) { // record dia dist
					diadist = sumdistmat[i-1][idx_dia] + distmat[i][k];
				} else diadist = -1;
				if (idx_dialeft < last_right && idx_dialeft >= last_left) { // record dia_left dist
					dialeftdist = sumdistmat[i-1][idx_dialeft] + distmat[i][k];
				} else dialeftdist = -1;
				if (idx_left < last_right && idx_left >= last_left) { // record left dist
					leftdist = sumdistmat[i-1][idx_left] + distmat[i][k];
				} else leftdist = -1;
                if (underdist != -1) { // compare under dist
					dist = underdist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist <= dist)
							dist = dialeftdist;
					if (underdiadist != - 1)
						if (underdiadist <= dist)
							dist = underdiadist;
					if (diadist != - 1)
						if (diadist <= dist)
							dist = diadist;
				} else if (underdiadist != -1) { // compare under_dia dist
					dist = underdiadist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist < dist)
							dist = dialeftdist;
					if (diadist != - 1)
						if (diadist <= dist)
							dist = diadist;
				} else if (diadist != -1) { // compare dia dist
					dist = diadist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
					if (dialeftdist != - 1)
						if (dialeftdist < dist)
							dist = dialeftdist;
				} else if (dialeftdist != -1) { // compare dia_left dist
					dist = dialeftdist;
					if (leftdist != - 1)
						if (leftdist < dist)
							dist = leftdist;
				} else dist = leftdist;
				if (dist == -1.0) {
					fprintf(stderr, "c error: calc_sumdistmat_asym: dist = -1.0 at [%d][%d]\n", i, k);
					exit(1);
				}
				sumdistmat[i][k] = dist;
				// set backtrack point for the best path
				if (dist == diadist) {
					pathmat1[i][k] = 1;
					pathmat2[i][k] = pathmat2[i-1][k-2];
				} else if (dist == underdiadist) {
					pathmat1[i][k] = 2;
					pathmat2[i][k] = pathmat2[i-1][k-1];
				} else if (dist == dialeftdist) {
					pathmat1[i][k] = 3;
					pathmat2[i][k] = pathmat2[i-1][k-3];
				} else if (dist == underdist) {
					pathmat1[i][k] = 4;
					pathmat2[i][k] = pathmat2[i-1][k];
				} else if (dist == leftdist) {
					pathmat1[i][k] = 5;
					pathmat2[i][k] = pathmat2[i-1][k-4];
				} else {
					fprintf(stderr, "c error: calc_pathmats_asym: dist=%lf, diadist=%lf, underdiadist=%lf, dialeftdist=%lf, underdist=%lf, leftdist=%lf\n", dist, diadist, underdiadist, dialeftdist, underdist, leftdist);
					exit(1);
				}
			}
		}
	}

	left = j-winl;
	if (left < 0) left = 0;
	right = j+winl_e-1;
	if (right >= org_row) right = org_row-1;
	// take the minimum accumulated distance and its optimum path
	for (k = right, cis = k, ri = i-1; k >= left; k--) {
		//fprintf(stderr, "b i=%d;  j=%d;  k=%d;  left=%d;  right=%d\n", i, j, k, left, right);
		//fprintf(stderr, "b i=%d;  j=%d;  k=%d;  left=%d;  right=%d;  sumdist=%lf\n", i, j, k, left, right, sumdistmat[i][k]);
		//fprintf(stderr, "b i=%d;  j=%d;  k=%d;  left=%d;  right=%d;  sumdist=%lf;  pathmat2=%d\n", i, j, k, left, right, sumdistmat[i][k], pathmat2[i][k]);
		if (pathmat2[ri][k] < 0)
			ci = -1 * pathmat2[ri][k];
		else ci = 0;
		meandist = sumdistmat[ri][k] / (ri + k - ci + 2);
		if (k < right) {
			if (meandist < minmeandist) {
				minmeandist = meandist;
				cis = k;
			}
		} else minmeandist = meandist;
	}

	// compute twf function
	for (ri = tar_row-1, ci = cis; ri > 0; ri--) {
		twfunc[ri][0] = ci;
		twfunc[ri][1] = 1;
		//fprintf(stderr, "b ri=%d;  ci=%d\n", ri, ci);
		if (pathmat1[ri][ci] == 1) {
			// dia
			ci -= 2;
		} else if (pathmat1[ri][ci] == 2) {
			// under_dia
			ci--;
		} else if (pathmat1[ri][ci] == 3) {
			// dia_left
			ci -= 3;
		} else if (pathmat1[ri][ci] == 5) {
			// left
			ci -= 4;
		} else {
			if (pathmat1[ri][ci] < 1) {
				fprintf(stderr, "error: c_calc_twfunc_asym: pathmat[%d][%d] = %d < 1, %d, %lf\n", ri, ci, pathmat1[ri][ci], pathmat2[ri][ci], sumdistmat[ri][ci]);
				exit(1);
			}
		}
		//fprintf(stderr, "b ri=%d; ===> ci=%d\n", ri, ci);
	}
	twfunc[ri][0] = ci;
	twfunc[ri][1] = 1;

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
