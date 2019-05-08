#include "GMRes.h"
#include "Multiplicate.h"
#include <math.h> 
#include <omp.h>
#include <string.h> 

int GMRes(struct CSR_matrix *m, double *b, double *x_n, double eps, int max_iterations) {

	int number = 10;
	double *r_j = (double*)malloc((m->num_rows) * sizeof(double));
	double *w = (double*)malloc((m->num_rows) * sizeof(double));
	double beta;
	double gamma;
	double k;
	double h_tilda;
	double h;
	double *z_j = (double*)malloc((m->num_rows) * sizeof(double));
	double *omega = (double*)malloc((m->num_rows) * sizeof(double));
	
	double *result_vector = (double*)malloc((m->num_rows) * sizeof(double));
	double *result_vector2 = (double*)malloc((m->num_rows) * sizeof(double));

	double *H = (double*)malloc((m->num_rows) * number * sizeof(double));
	double *V = (double*)malloc((m->num_rows) * number * sizeof(double));
	double *Z = (double*)malloc((m->num_rows) * number * sizeof(double));
	double *G = (double*)malloc((m->num_rows) * number * sizeof(double));
	
	double *y_m = (double*)malloc(number * sizeof(double));
	double *Ax_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *c = (double*)malloc(number * sizeof(double));
	double *s = (double*)malloc(number * sizeof(double));
	double *g = (double*)malloc(number * sizeof(double));

	for (int i = 0; i < number; i++) {
		c[i] = 1;
		s[i] = 1;
	}


	//Начальное приближение
	Ax_n = spnv_pointer(*m, x_n);

	for (int i = 0; i < m->num_rows; i++) {
		r_j[i] = b[i] - Ax_n[i];
	}

	beta = sqrt(dot_product(r_j, r_j, m->num_rows));
	gamma = beta;
	double antibeta = 1 / beta;

	for (int i = 0; i < m->num_rows; i++) {
		w[i] = r_j[i];
	}

	k = number;

	for (int j = 1; j < number; j++) {
		
		for (int i = 0; i < m->num_rows; i++) {
				V[i*j + i] = w[i] * antibeta;
				//Without Preconditioner
				Z[i*j + i] = V[i*j + i];
			}
		
		w = spnv_pointer(*m, &Z[j*m->num_rows]);

		for (int i = 0; i <= j; i++) {
			H[i*j + i] = dot_product(w, &V[i*m->num_rows], m->num_rows);
			for (int d = 0; d < m->num_rows; d++) {
				w[d] = w[d] + H[i*j + i] * V[i*m->num_rows + d];
			}
		}

		for (int i = 0; i <= (j - 1); i++) {

			h_tilda = H[i*j + i];
			h = H[i*j + i + 1];
			H[i*j + i] = c[i] * h_tilda + s[i] * h;
			H[i*j + i + 1] = -s[i] * h_tilda + c[i] * h;
		}

		h_tilda = H[j*number + j];
		beta = dot_product(w, w, m->num_rows);

		H[j*number + j] = sqrt(h_tilda*h_tilda + beta * beta);

		c[j] = beta / H[j*number + j];
		s[j] = h_tilda / H[j*number + j];

		g[j] = gamma * c[j];
		
		gamma = -gamma * s[j];

		if (gamma*gamma < eps*eps) {
			k = j;
			j = m + 1;
		}

	}

	for (int i = 0; i < m->num_rows; i++)
	{
		result_vector[i] = 0.0;
		for (int j = 0; j < number; j++)
		{
			result_vector[i] += H[i*j+j] * g[j];
		}
	}

	for (int i = 0; i < m->num_rows; i++)
	{
		result_vector2[i] = 0.0;
		for (int j = 0; j < number; j++)
		{
			result_vector2[i] += Z[i*j + j] * result_vector[j];
		}
	}

	for (int i = 0; i < m->num_rows; i++) {
		x_n[i] = x_n[i] + result_vector2[i];
	}

	printf("x[num_rows-1] = %.40lf\n", x_n[m->num_rows - 1]);
}