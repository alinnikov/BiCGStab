#include "BiCGStab.h"
#include "Multiplicate.h"
#include <math.h> 
#include <omp.h>


double alpha(struct CSR_matrix *m, double *r_0_help, double r0help_rn, double *p_n, double *Ap_n) {
	
	

	double numerator = r0help_rn;//scalar_multiplicate(r_0_help, r_n, m->num_rows);
	double denumerator = dot_product(Ap_n, r_0_help, m->num_rows);
	//free(Ap_n);
	return numerator / denumerator;
	
}

int s(struct CSR_matrix *m, double *r_n, double *p_n, double alpha_n, double* s_n, double *Ap_n) {
	
	//multiplicate(*m, p_n, Ap_n);

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

	//multiplicate(*m, s_n, As_n);
	double numerator = dot_product(As_n, s_n, m->num_rows);
	double denumerator = dot_product(As_n, As_n, m->num_rows);

	return numerator / denumerator;

}


int BiCGStab(struct CSR_matrix *m, double *b){//, double *x_n), double eps, int numbner_of_iterations) {
//#pragma omp parallel
	//Создаём необходимые массивы
	double *x_n = (double*)malloc((m->num_rows) * sizeof(double));
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

	double L2_norm = 1;
	double eps = 0.000000000001;
	int number_of_iterations=0;
	//Задаём начальный вектор
	for (int i = 0; i < m->num_rows; i++) {
		x_n[i] = 10.0;
	}
	// Задаём A*x0

	Ax_n = spnv_pointer(*m, x_n);
	
	//Задаём r0*, r0 и p0

	for (int i = 0; i < m->num_rows; i++) {
		r_n[i] = b[i] - Ax_n[i];
		r_0_help[i] = r_n[i];
		p_n[i] = r_n[i];
	}
	//int num;
	//Цикл
//#pragma omp parallel private(x_n,x_n_plus, r_n, r_n_plus, r_n_help,p_n,p_n_plus,AX_n,s_n,As_n,AP_n,alpha_n,omega_n, beta_n, m)
	for (int num = 0; num < 30000000 && L2_norm>eps*eps; num++) {

		Ap_n = spnv_pointer(*m, p_n);

		r0help_rn = dot_product(r_0_help, r_n, m->num_rows);

		alpha_n = alpha(m, r_0_help, r0help_rn, p_n, Ap_n);

		s(m, r_n, p_n, alpha_n, s_n, Ap_n);

		As_n = spnv_pointer(*m, s_n);

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

		if (number_of_iterations % 1000 == 0) {
		//	printf("%.40lf\n", x_n[m->num_rows-1]);
		//	printf("L2_norm=%.40lf\n", L2_norm);
		}

		//double abs = 0;
		//for (int i = 0; i < m->num_rows; i++) {
		//	abs += fabs(r_n[i])*fabs(r_n[i]);
		//}
		L2_norm = dot_product(r_n,r_n,m->num_rows);


		number_of_iterations++;

		free(Ap_n);
		free(As_n);
	}
	
	//printf("%lf\n", r_n[m->num_rows]);
	printf("Number of iterations = %d\n", number_of_iterations);
	//printf("L2_norm=%.40lf\n", L2_norm);
	//printf("x[num_rows-1] = %lf\n", x_n[m->num_rows-1]);
	return 0;
}