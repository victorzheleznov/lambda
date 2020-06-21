#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#include "lambda.h"

void lambda(int dim, int num_cands, double* trian_matr, double* diag_matr, double* est, int* found, candidates* cands)
{
	// input:
	// dim - dimension of estimate
	// num_cands - requested number of candidates
	// trian_matr - array of size (dim * dim) with unit lower triangular matrix of L^TDL decomposition of the covariance matrix of the float ambiguities
	// diag_matr - array of size (dim) with diagonal matrix of L^TDL decomposition of the covariance matrix of the float ambiguities
	// est - array of size (dim) with float ambiguities
	//
	// output:
	// found - number of candidates found by algorithm
	// cands - array of structures (size (num_cands)) with integer ambiguities (size (dim)) and norm values
	//
	// trian_matr, diag_mat, est will be overwritten!
	//
	// references:
	// 1. De Jonge P, Tiberius C (1996) The LAMBDA method of intger ambiguity estimation: implementation aspects

	int i = 0;
	int j = 0;
	int k = 0;

	int* transf_matr = NULL;
	int* temp = NULL;

	// memory allocation

	transf_matr = (int*) malloc(dim * dim * sizeof(int));
	temp = (int*) malloc(dim * sizeof(int));

	// main algorithm

	transformation(dim, trian_matr, diag_matr, est, transf_matr);
	search(dim, num_cands, trian_matr, diag_matr, est, found, cands);

	// back transformation

	for (k = 0; k < *found; k++)
	{
		for (i = 0; i < dim; i++)
		{
			temp[i] = 0;
			for (j = 0; j < dim; j++) temp[i] += transf_matr[i * dim + j] * cands[k].est[j];
		}

		for (i = 0; i < dim; i++) cands[k].est[i] = temp[i];
	}

	// clean up memory

	free(transf_matr);
	free(temp);

	return;
}
 
void transformation(int dim, double* trian_matr, double* diag_matr, double* est, int* transf_matr)
{
	// input:
	// dim - dimension of estimate
	// trian_matr - array of size (dim * dim) with unit lower triangular matrix of L^TDL decomposition of the covariance matrix of the float ambiguities
	// diag_matr - array of size (dim) with diagonal matrix of L^TDL decomposition of the covariance matrix of the float ambiguities
	// est - array of size (dim) with float ambiguities
	//
	// output:
	// transf_matr - transpose and inverse matrix to integer Z-transformation matrix
	//
	// trian_matr, diag_mat, est will be transformed to correspond to new ambiguities!
	//
	// references:
	// 1. De Jonge P, Tiberius C (1996) The LAMBDA method of intger ambiguity estimation: implementation aspects

	int swap = 1;
	int i = 0;
	int j = 0;
	int k = 0;
	int i_last = dim - 2;
	
	// buffer variables
	
	int mu_i = 0;
	int temp_i = 0;

	double mu = 0.0;
	double delta = 0.0;
	double eta = 0.0;
	double lambda_1 = 0.0;
	double lambda_2 = 0.0;
	double lambda_3 = 0.0;
	double temp = 0.0;

	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
		{
			if (i == j) transf_matr[i * dim + j] = 1;
			else transf_matr[i * dim + j] = 0;
		}
	}

	while (swap == 1)
	{
		i = dim - 1;
		swap = 0;

		while (swap == 0 && i > 0)
		{
			i--;

			if (i <= i_last)
			{
				// integer Gauss transformation to column i 

				for (j = i + 1; j < dim; j++)
				{
					mu = round(trian_matr[j * dim + i]);
					mu_i = (int) mu;

					if (mu_i != 0)
					{
						for (k = 0; k < dim; k++)
						{
							transf_matr[k * dim + j] += mu_i * transf_matr[k * dim + i];
							if (k >= j) trian_matr[k * dim + i] -= mu * trian_matr[k * dim + j];
						}

						est[i] -= mu * est[j];
					}
				}
			}

			delta = diag_matr[i] + pow(trian_matr[(i + 1) * dim + i], 2.0) * diag_matr[i + 1];

			if (delta < diag_matr[i + 1])
			{
				lambda_3 = diag_matr[i + 1] * trian_matr[(i + 1) * dim + i] / delta;
				eta = diag_matr[i] / delta;

				diag_matr[i] = eta * diag_matr[i + 1];
				diag_matr[i + 1] = delta;

				for (j = 0; j < i; j++)
				{
					lambda_1 = trian_matr[i * dim + j];
					lambda_2 = trian_matr[(i + 1) * dim + j];

					trian_matr[i * dim + j] = -trian_matr[(i + 1) * dim + i] * lambda_1 + lambda_2;
					trian_matr[(i + 1) * dim + j] = eta * lambda_1 + lambda_2 * lambda_3;
				}

				trian_matr[(i + 1) * dim + i] = lambda_3;

				for (k = 0; k < dim; k++)
				{
					temp_i = transf_matr[k * dim + i];
					transf_matr[k * dim + i] = transf_matr[k * dim + (i + 1)];
					transf_matr[k * dim + (i + 1)] = temp_i;

					if (k >= i + 2)
					{
						temp = trian_matr[k * dim + i];
						trian_matr[k * dim + i] = trian_matr[k * dim + (i + 1)];
						trian_matr[k * dim + (i + 1)] = temp;
					}
				}

				temp = est[i];
				est[i] = est[i + 1];
				est[i + 1] = temp;

				i_last = i;
				swap = 1;
			}
		}
	}
}

void search(int dim, int num_cands, double* trian_matr, double* diag_matr, double* est, int* found, candidates* cands)
{
	// counters

	int i = 0;
	int j = 0;
	int k = 0;
	
	// flags

	int end_search = 0;
	int stop_backtrack = 0;
	int last_coord_change = 0;

	// other variables
	
	int i_old = 0;
	int num_cand = 0;
	int pos_max = 0;

	double ellipsoid_bound = 0.0;
	double sum = 0.0;
	double reach = 0.0;
	double delta = 0.0;
	double norm = 0.0;
	double norm_max = 0.0;
	double temp = 0.0;

	double* trian_inv = NULL;
	double* diag_inv = NULL;
	double* lef = NULL;
	double* left = NULL;
	double* right = NULL;
	double* end = NULL;
	double* dist = NULL;

	// memory allocation

	trian_inv = (double*) malloc(dim * dim * sizeof(double));
	diag_inv = (double*) malloc(dim * sizeof(double));
	lef = (double*) malloc(dim * sizeof(double));
	left = (double*) malloc((dim + 1) * sizeof(double));
	right = (double*) malloc((dim + 1) * sizeof(double));
	end = (double*) malloc(dim * sizeof(double));
	dist = (double*) malloc(dim * sizeof(double));

	// inversion to obtain decomposition of the inverse covariance matrix

	for (i = 0; i < dim; i++)
	{
		diag_inv[i] = 1.0 / diag_matr[i];
		trian_inv[i * dim + i] = 1.0 / trian_matr[i * dim + i];

		for (j = i + 1; j < dim; j++)
		{
			sum = 0.0;
			for (k = i; k < j; k++) sum += trian_matr[j * dim + k] * trian_inv[k * dim + i];

			trian_inv[j * dim + i] = -sum / trian_matr[j * dim + j];
		}
	}

	// size of the search ellipsoid

	ellipsoid_bound = set_ellipsoid_bound(dim, num_cands, 1.5, trian_inv, diag_inv, est);
	
	// initial conditions

	i = dim;
	i_old = i;
	right[dim] = ellipsoid_bound;
	left[dim] = 0.0;

	// main loop

	while (end_search == 0)
	{
		i--;

		if (i_old <= i) lef[i] = lef[i] + trian_inv[(i + 1) * dim + i];
		else
		{
			lef[i] = 0.0;
			for (j = i + 1; j < dim; j++) lef[i] += trian_inv[j * dim + i] * dist[j];
		}
		i_old = i;

		if (i != dim - 1) right[i] = (right[i + 1] - left[i + 1]) * diag_inv[i + 1] / diag_inv[i];
		else right[i] = (right[i + 1] - left[i + 1]) / diag_inv[i];

		reach = sqrt(right[i]);
		delta = est[i] - reach - lef[i];
		dist[i] = ceil(delta) - est[i];

		if (dist[i] > reach - lef[i])
		{
			// backtrack

			stop_backtrack = 0;
			last_coord_change = 0;

			while (!stop_backtrack && i < dim - 1)
			{
				i++;

				if (dist[i] < end[i])
				{
					dist[i] += 1.0;
					left[i] = dist[i] + lef[i];
					left[i] = pow(left[i], 2.0);

					stop_backtrack = 1;
					if (i == dim - 1) last_coord_change = 1;
				}
			}

			if (i == dim - 1 && !last_coord_change) end_search = 1;
		}
		else
		{
			// set the right border

			end[i] = reach - lef[i] - 1.0;
			left[i] = dist[i] + lef[i];
			left[i] = pow(left[i], 2.0);
		}

		if (i == 0)
		{
			// collect

			norm = ellipsoid_bound - (right[0] - left[0]) * diag_inv[0];
			end[0] += 1.0;

			while (dist[0] <= end[0])
			{
				if (num_cand < num_cands)
				{
					cands[num_cand].norm = norm;
					for (j = 0; j < dim; j++) 
					{
						temp = dist[j] + est[j];
						cands[num_cand].est[j] = (int) round(temp);
					}

					num_cand++;
				}
				else
				{
					norm_max = cands[0].norm;
					pos_max = 0;

					for (j = 1; j < num_cands; j++)
					{
						if (cands[j].norm > norm_max) 
						{
							norm_max = cands[j].norm;
							pos_max = j;
						}
					}

					if (norm < norm_max)
					{
						cands[pos_max].norm = norm;
						for (j = 0; j < dim; j++) 
						{
							temp = dist[j] + est[j];
							cands[pos_max].est[j] = (int) round(temp);
						}
					}
				}

				norm += (2 * (dist[0] + lef[0]) + 1.0) * diag_inv[0];
				dist[0] += 1.0;
			}

			// backtrack

			stop_backtrack = 0;
			last_coord_change = 0;

			while (!stop_backtrack && i < dim - 1)
			{
				i++;

				if (dist[i] < end[i])
				{
					dist[i] += 1.0;
					left[i] = dist[i] + lef[i];
					left[i] = pow(left[i], 2.0);

					stop_backtrack = 1;
					if (i == dim - 1) last_coord_change = 1;
				}
			}

			if (i == dim - 1 && !last_coord_change) end_search = 1;
		}
	}

	// clean up memeory

	free(trian_inv);
	free(diag_inv);
	free(lef);
	free(left);
	free(right);
	free(end);
	free(dist);

	// check for errors

	if (num_cand < num_cands) 
	{
		printf("only %d of %d candidates were found!\n", num_cand, num_cands);

		*found = num_cand;
	}
	else *found = num_cands;

	// sorting according to the norm

	qsort(cands, *found, sizeof(candidates), compare_candidates);

	return;
}

int compare_candidates(const void* first, const void* second)
{
	const candidates* first_cand = (const candidates*) first;
	const candidates* second_cand = (const candidates*) second;

	if (first_cand->norm < second_cand->norm) return -1;
	else if (first_cand->norm > second_cand->norm) return 1;
	else return 0;
}

double set_ellipsoid_bound(int dim, int num_cands, double factor, double* trian_inv, double* diag_inv, double* est)
{
	int i = 0;
	int j = 0;
	int k = 0;
	int num_bounds = 0;

	double sign = 0.0;
	double volume = 0.0;
	double det = 0.0;
	double power = 0.0;
	double ellipsoid_bound = 0.0;
	
	const double pi = 3.14159265358979323846;

	double* cand = NULL;
	double* sum = NULL;
	double* bounds = NULL;

	if (num_cands <= dim + 1)
	{
		// memory allocation

		cand = (double*) malloc(dim * sizeof(double));
		sum = (double*) malloc(dim * sizeof(double));
		bounds = (double*) malloc((dim + 1) * sizeof(double));

		// norm of the first candidate (rounding each coordinate to the nearest integer)

		for (i = 0; i < dim; i++) cand[i] = round(est[i]);

		bounds[0] = 0.0;
		for (i = 0; i < dim; i++)
		{
			sum[i] = 0.0;
			for (j = i; j < dim; j++) sum[i] += trian_inv[j * dim + i] * (cand[j] - est[j]);

			bounds[0] += diag_inv[i] * pow(sum[i], 2.0);
		}

		for (k = 0; k < dim; k++)
		{
			if (est[k] - cand[k] > 0.0) sign = 1.0;
			else if (est[k] - cand[k] < 0.0) sign = -1.0;
			else sign = 0.0;

			bounds[k + 1] = bounds[0];
			for (i = 0; i <= k; i++) bounds[k + 1] += diag_inv[i] * (2.0 * sum[i] * trian_inv[k * dim + i] * sign + pow(trian_inv[k * dim + i], 2.0) * pow(sign, 2.0));
		}

		free(cand);
		free(sum);

		// sorting

		num_bounds = dim + 1;
		qsort(bounds, num_bounds, sizeof(double), compare_norms);

		ellipsoid_bound = bounds[num_cands - 1];
		ellipsoid_bound += 1e-6;	// to make sure there is no boundary problem

		free(bounds);
	}
	else
	{
		// volume of an n-ball

		if (dim % 2 == 0) 
		{
			volume = pi;
			for (k = 4; k <= dim; k = k + 2) volume = 2.0 * pi * volume / ((double) k);
		}
		else 
		{
			volume = 2.0;
			for (k = 3; k <= dim; k = k + 2) volume = 2.0 * pi * volume / ((double) k);
		}

		// determinant of diagonal matrix 

		det = 1.0;
		for (i = 0; i < dim; i++) det *= diag_inv[i];

		// ellipsoid bound calculation
		
		power = 2.0 / ((double) dim);
		ellipsoid_bound = num_cands * sqrt(det) / volume;	// volume -> sqrt(volume) to expand search ellipsoid
		ellipsoid_bound = pow(ellipsoid_bound, power);
		ellipsoid_bound *= factor;
	}

	return ellipsoid_bound;
}

int compare_norms(const void* first, const void* second)
{
	const double* first_norm = (const double*) first;
	const double* second_norm = (const double*) second;

	if (*first_norm < *second_norm) return -1;
	else if (*first_norm > *second_norm) return 1;
	else return 0;
}

void decomposition(int dim, double* matr, double* trian_matr, double* diag_matr)
{
	// L^TDL

	int i = 0;
	int j = 0;
	int k = 0;

	for (i = dim - 1; i >= 0; i--)
	{
		diag_matr[i] = matr[i * dim + i];
		for (k = 0; k <= i; k++) trian_matr[i * dim + k] = matr[i * dim + k] / sqrt(matr[i * dim + i]);

		for (j = 0; j < i; j++) 
		{
			for (k = 0; k <= j; k++) matr[j * dim + k] -= trian_matr[i * dim + k] * trian_matr[i * dim + j];
		}

		for (k = 0; k <= i; k++) trian_matr[i * dim + k] /= trian_matr[i * dim + i];
	}
}