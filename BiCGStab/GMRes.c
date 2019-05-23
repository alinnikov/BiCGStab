#include "GMRes.h"
#include "Multiplicate.h"
#include <math.h> 
#include <omp.h>
#include <string.h> 

int GMRes(struct CSR_matrix *A, double *b, double *x, double tol, int max_iter, int m) {

	int k;
	int iter;
	double beta;
	double gamma;
	double h;
	double h_tilda;
	
	double *r = (double*)malloc((A->num_rows) * sizeof(double));
	double *Ax = (double*)malloc((A->num_rows) * sizeof(double));
	
	double *w = (double*)malloc((A->num_rows) * sizeof(double));

	double **Hess = (double **)malloc(m * sizeof(double));
	for (int i = 0; i < m; i++) {
		Hess[i] = (double *)malloc(m * sizeof(double));
	}

	double *V = (double *)malloc((A->num_rows) * m * sizeof(double));
	double *Z = (double *)malloc((A->num_rows) * m * sizeof(double));
					    
	double *c = (double *)malloc(m * sizeof(double));
	double *s = (double *)malloc(m * sizeof(double));
	double *g = (double *)malloc(m * sizeof(double));

	iter = 0;
	
Start:
	
	// initial guess
	Ax = spmv(*A, x);
	for (int i = 0; i < A->num_rows; i++) {
		r[i] = b[i] - Ax[i];
		w[i] = r[i];
	}
	beta = sqrt(dot_product(r, r, A->num_rows));
	gamma = beta;
	k = m;

	// main loop
	for (int j = 0; j < m; j++) {
		
		for (int i = 0; i < A->num_rows; i++) {
			V[j * A->num_rows + i] = w[i] / beta;
			// preconditioning
			Z[j * A->num_rows + i] = V[j * A->num_rows + i];
		}
		
		// Arnoldi orthogonalization
		w = spmv(*A, &Z[j * A->num_rows]);
		
		for (int i = 0; i <= j; i++) {
			Hess[i][j] = dot_product(w, &V[i * A->num_rows], A->num_rows);
			for (int l = 0; l < A->num_rows; l++) {
				w[l] -= Hess[i][j] * V[i * A->num_rows + l];
			}	
		}

		// Givens rotations
		for (int i = 0; i < j; i++) {
			h_tilda = Hess[i][j];
			h = Hess[i + 1][j];
			Hess[i][j] = c[i] * h_tilda + s[i] * h;
			Hess[i + 1][j] = -s[i] * h_tilda + c[i] * h;
		}
		h_tilda = Hess[j][j];
		beta = sqrt(dot_product(w, w, A->num_rows));
		Hess[j][j] = sqrt(h_tilda * h_tilda + beta * beta);
		s[j] = beta / Hess[j][j];
		c[j] = h_tilda / Hess[j][j];
		g[j] = gamma * c[j];	
		gamma = -gamma * s[j];
		if (fabs(gamma) < tol) {
			k = j + 1;
			j = m + 1;
		}

	}

	iter += k;

	//printf("m = %d\n", m);
	//printf("krylov_subspace_dim = %d\n", k);
	//printf("gamma_residual_norm = %e\n", gamma);

	double *y = (double *)malloc(k * sizeof(double));

	// y = H^(-1) * g
	for (int i = k - 1; i >= 0; i--) {
		y[i] = g[i];
		for (int j = i + 1; j < k; j++) {
			y[i] -= Hess[i][j] * y[j];
		}
		y[i] /= Hess[i][i];
	}

	// x = x + Z * y
	for (int i = 0; i < A->num_rows; i++) {
		for (int j = 0; j < k; j++) {
			x[i] += Z[j * A->num_rows + i] * y[j];
		}
	}

	free(y);

	/*
	Ax = spmv(*A, x);
	for (int i = 0; i < A->num_rows; i++) {
		r[i] = b[i] - Ax[i];
	}
	printf("true__residual_norm = %e\n", sqrt(dot_product(r, r, A->num_rows)));
	*/

	if ((fabs(gamma) > tol) && (iter <= max_iter)) {
		goto Start;
	}

	printf("iter = %d, residual_norm = %e\n", iter, gamma);

	for (int i = 0; i < m; i++) {
		free(Hess[i]);
	}
	free(Hess);
	free(r);
	free(Ax);
	free(w);
	free(V);
	free(Z);
	free(c);
	free(s);
	free(g);
	
	return 0;
}