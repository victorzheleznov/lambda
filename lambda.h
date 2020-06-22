#ifndef LAMBDA_H
#define LAMBDA_H

// structures

typedef struct {
	double norm;	// squared norm!
	int* est; 
} candidates;

// functions

void lambda(int dim, int num_cands, double* trian_matr, double* diag_matr, double* est, int* found, candidates* cands);
void decomposition(int dim, double* matr, double* trian_matr, double* diag_matr);
void transformation(int dim, double* trian_matr, double* diag_matr, double* est, int* transf_matr);
void search(int dim, int num_cands, double* trian_matr, double* diag_matr, double* est, int* found, candidates* cands);
double set_ellipsoid_bound(int dim, int num_cands, double factor, double* trian_inv, double* diag_inv, double* est);
int compare_candidates(const void* first, const void* second);
int compare_norms(const void* first, const void* second);

#endif //LAMBDA_H