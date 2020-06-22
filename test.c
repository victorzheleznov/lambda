#include "crossplatform.h"

#include <conio.h>
#include <stdio.h>
#include <stdlib.h>

#include "test.h"
#include "lambda.h"


int main(void)
{
	// counters

	int i = 0;
	int j = 0;

	// parameters

	int dim = 0;
	int num_cands = 0;

	// flag

	int perform_zero_test = 0;

	// input data

	double* cov_matr = NULL;
	double* est = NULL;
	int* true_amb = NULL;

	// matrices

	double* trian_matr = NULL;
	double* diag_matr = NULL;

	// candidates

	candidates* cands = {0};
	int found = 0;

	// input files

	FILE* cov_matr_txt = NULL;
	FILE* est_txt = NULL;
	FILE* true_amb_txt = NULL;

	// output_files

	FILE* transformation_matrix_txt = NULL;
	FILE* candidates_txt = NULL;
	FILE* sq_norm_txt = NULL;
	FILE* zero_test_txt = NULL;

	char cwd[BUFFER_SIZE];


	// dimension and number of candidates

	printf("Enter dimension: ");
	scanf("%d", &dim);

	printf("Enter desired number of candidates: ");
	scanf("%d", &num_cands);

	// open input files

	getcwd(cwd, BUFFER_SIZE);
	mkdir(INPUT_FOLDER);
	chdir(INPUT_FOLDER);

	cov_matr_txt = fopen(COVARIANCE_MATR_FILE, "r");
	est_txt = fopen(ESTIMATE_FILE, "r");
	true_amb_txt = fopen(TRUE_AMBIGUITIES_FILE, "r");

	if (true_amb_txt != NULL) 
		perform_zero_test = 1;

	chdir(cwd);

	// memory allocation for input data

	cov_matr = (double*) malloc(dim * dim * sizeof(double));
	est = (double*) malloc(dim * sizeof(double));

	if (perform_zero_test) 
		true_amb = (int*) malloc(dim * sizeof(int));

	// memory allocation for decomposition

	trian_matr = (double*) malloc(dim * dim * sizeof(double));
	diag_matr = (double*) malloc(dim * sizeof(double));

	// memory allocation for candidates

	cands = (candidates*) malloc(num_cands * sizeof(candidates));

	for (int j = 0; j < num_cands; j++)
	{
		cands[j].est = NULL;
		cands[j].norm = 0.0;

		cands[j].est = (int*) malloc(dim * sizeof(int));
	}

	// read input data

	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
		{
			fscanf(cov_matr_txt, "%lf", &cov_matr[i * dim + j]);
		}

		fscanf(est_txt, "%lf", &est[i]);

		if (perform_zero_test)
			fscanf(true_amb_txt, "%d", &true_amb[i]);
	}

	// LAMBDA

	decomposition(dim, cov_matr, trian_matr, diag_matr);
	lambda(dim, num_cands, trian_matr, diag_matr, est, &found, cands);

	// output files

	mkdir(OUTPUT_FOLDER);
	chdir(OUTPUT_FOLDER);

	candidates_txt = fopen(CANDIDATES_FILE, "w");
	sq_norm_txt = fopen(SQUARED_NORM_FILE, "w");
	
	if (perform_zero_test) 
		zero_test_txt = fopen(ZERO_TEST_FILE, "w");

	chdir(cwd);

	// output to files

	for (int i = 0; i < found; i++)
	{
		fprintf(sq_norm_txt, "%-16.8lf", cands[i].norm);

		for (int j = 0; j < dim; j++)
		{
			fprintf(candidates_txt, "%-5d", cands[i].est[j]);

			if (perform_zero_test)
				fprintf(zero_test_txt, "%-5d", cands[i].est[j] - true_amb[j]);
		}

		fprintf(candidates_txt, "\n");

		if (perform_zero_test)
			fprintf(zero_test_txt, "\n");
	}

	// cleaning up memory

	free(cov_matr);
	free(est);

	if (perform_zero_test)
		free(true_amb);

	free(trian_matr);
	free(diag_matr);

	for (int i = 0; i < num_cands; i++)
	{
		free(cands[i].est);
	}
	free(cands);

	// closing files

	fclose(cov_matr_txt);
	fclose(est_txt);
	fclose(candidates_txt);
	fclose(sq_norm_txt);

	if (perform_zero_test)
	{
		fclose(true_amb_txt);
		fclose(zero_test_txt);
	}

	return 1;
}