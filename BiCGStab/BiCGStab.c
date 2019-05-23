#include "BiCGStab.h"
#include "Multiplicate.h"
#include <math.h> 
#include <omp.h>


double alpha(struct CSR_matrix *m, double *r_0_help, double r0help_rn, double *p_n, double *Ap_n) {
	
	

	double numerator = r0help_rn;
	double denumerator = dot_product(Ap_n, r_0_help, m->num_rows);
	return numerator / denumerator;
	
}

int s(struct CSR_matrix *m, double *r_n, double *p_n, double alpha_n, double* s_n, double *Ap_n) {
	

	for (int i = 0; i < m->num_rows; i++) {
		s_n[i] = r_n[i] - alpha_n * Ap_n[i];
	}
	return 0;
}

double beta(struct CSR_matrix *m, double r0help_rn, double *r_n, double * r_0_help, double alpha_n, double omega_n) {

	double numerator = alpha_n * dot_product(r_n, r_0_help , m->num_rows);
	double denumerator = omega_n * r0help_rn;

	return numerator / denumerator;
}

double omega(struct CSR_matrix *m, double *s_n, double *As_n) {

	double numerator = dot_product(As_n, s_n, m->num_rows);
	double denumerator = dot_product(As_n, As_n, m->num_rows);

	return numerator / denumerator;

}


int BiCGStab(struct CSR_matrix *m, double *b, double *x_n, double tol, int max_iter) {
	//Создаём необходимые массивы
	double *r_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *r_0_help = (double*)malloc((m->num_rows) * sizeof(double));
	double *p_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *Ax_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *s_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *As_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *Ap_n = (double*)malloc((m->num_rows) * sizeof(double));
	double alpha_n;
	double omega_n;
	double beta_n;
	double r0help_rn;

	double L2_norm = 1.0;
	int number_of_iterations = 0;
	// Задаём A*x0

	Ax_n = spmv(*m, x_n);
	
	//Задаём r0*, r0 и p0

	for (int i = 0; i < m->num_rows; i++) {
		r_n[i] = b[i] - Ax_n[i];
		r_0_help[i] = r_n[i];
		p_n[i] = r_n[i];
	}

	//Цикл
	for (int num = 0; num < max_iter && L2_norm>tol * tol; num++) {

		Ap_n = spmv(*m, p_n);

		r0help_rn = dot_product(r_0_help, r_n, m->num_rows);

		alpha_n = alpha(m, r_0_help, r0help_rn, p_n, Ap_n);

		s(m, r_n, p_n, alpha_n, s_n, Ap_n);

		As_n = spmv(*m, s_n);

		omega_n = omega(m, s_n, As_n);

		for (int i = 0; i < m->num_rows; i++) {
			x_n[i] = x_n[i] + alpha_n * p_n[i] + omega_n * s_n[i];
		}


		for (int i = 0; i < m->num_rows; i++) {
			r_n[i] = s_n[i] - omega_n * As_n[i];
		}
	
		beta_n = beta(m, r0help_rn, r_0_help, r_n, alpha_n, omega_n);

		for (int i = 0; i < m->num_rows; i++) {
			p_n[i] = r_n[i] - beta_n * omega_n * Ap_n[i] + beta_n * p_n[i];
		}

		if (number_of_iterations % 100 == 0) {
			printf("%.40lf\n", x_n[m->num_rows-1]);
			printf("L2_norm=%.40lf\n", L2_norm);
		}

		L2_norm = dot_product(r_n,r_n,m->num_rows);
		number_of_iterations++;
		free(Ap_n);
		free(As_n);
	}
	
	printf("Number of iterations = %d\n", number_of_iterations);
	printf("L2_norm=%.60lf\n", L2_norm);
	printf("x[num_rows-1] = %.40lf\n", x_n[m->num_rows-1]);
	return 0;
}