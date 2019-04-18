#include "BiCGStab.h"
#include "Multiplicate.h"
#include <math.h> 

double alpha(struct CSR_matrix **m, double *r_0_help, double *r_n, double *p_n, double *Ap_n) {
	
	//double *result = (double*)malloc((m->num_rows) * sizeof(double));
	//double *Ap_n = (double*)malloc(((*m)->num_rows) * sizeof(double));
	multiplicate(*m, p_n, Ap_n);

	double numerator = scalar_multiplicate(r_0_help, r_n, (*m)->num_rows);
	double denumerator = scalar_multiplicate(Ap_n, r_0_help, (*m)->num_rows);
	return numerator / denumerator;
}

double s(struct CSR_matrix **m, double *r_n, double *p_n, double alpha_n, double* s_n, double *Ap_n) {
	
	//double *Ap_n = (double*)malloc(((*m)->num_rows) * sizeof(double));
	multiplicate(*m, p_n, Ap_n);

	for (int i = 0; i <= (*m)->num_rows; i++) {
		s_n[i] = r_n[i] - alpha_n * Ap_n[i];
	}

}

double beta(struct CSR_matrix **m, double * r_n_plus, double * r_0_help, double *r_n, double alpha_n, double omega_n) {

	double numerator = alpha_n * scalar_multiplicate(r_n_plus, r_0_help , (*m)->num_rows);
	double denumerator = omega_n * scalar_multiplicate(r_n_plus, r_0_help, (*m)->num_rows);

	return numerator / denumerator;
}

double omega(struct CSR_matrix **m, double *s_n, double *As_n) {
	//double *As_n = (double*)malloc(((*m)->num_rows) * sizeof(double));

	multiplicate(*m, s_n, As_n);
	double numerator = scalar_multiplicate(As_n, s_n, (*m)->num_rows);
	double denumerator = scalar_multiplicate(As_n, As_n, (*m)->num_rows);

	return numerator / denumerator;

}


int BiCGStab(struct CSR_matrix *m, double *b) {
	//Создаём необходимые массивы
	double *x_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *x_n_plus = (double*)malloc((m->num_rows) * sizeof(double));
	double *r_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *r_n_plus = (double*)malloc((m->num_rows) * sizeof(double));
	double *r_0_help = (double*)malloc((m->num_rows) * sizeof(double));
	double *p_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *p_n_plus = (double*)malloc((m->num_rows) * sizeof(double));
	double *multi_result = (double*)malloc((m->num_rows) * sizeof(double));
	double *s_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *As_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *Ap_n = (double*)malloc((m->num_rows) * sizeof(double));
	double alpha_n;
	double omega_n;
	double beta_n;

	double L2_norm = 1;
	double eps = 0.00000001;

	//Задаём начальный вектор
	for (int i = 0; i <= m->num_rows; i++) {
		x_n[i] = 1000;
	}
	// Задаём A*x0

	multiplicate(m, x_n, multi_result);
	
	//Задаём r0*, r0 и p0

	for (int i = 0; i <= m->num_rows; i++) {
		r_n[i] = b[i] - multi_result[i];
		r_0_help[i] = r_n[i];
		p_n[i] = r_n[i];
	}
	//Цикл
	for (int num = 0; num < 300000 && L2_norm>eps; num++) {
		alpha_n = alpha(&m, r_0_help, r_n, p_n, Ap_n);
		
		s(&m, r_n, p_n, alpha_n, s_n, Ap_n);

		omega_n = omega(&m, s_n, As_n);

		for (int i = 0; i <= m->num_rows; i++) {
			x_n_plus[i] = x_n[i] + alpha_n * p_n[i] + omega_n * s_n[i];
		}

		multiplicate(m, s_n, As_n);
		for (int i = 0; i <= m->num_rows; i++) {
			r_n_plus[i] = s_n[i] - omega_n * As_n[i];
		}

		beta_n = beta(&m, r_n_plus, r_0_help, r_n, alpha_n, omega_n);

		multiplicate(m, s_n, Ap_n);
		for (int i = 0; i <= m->num_rows; i++) {
			p_n_plus[i] = r_n_plus[i] - beta_n * omega_n * Ap_n[i] + beta_n * p_n[i];
		}

		for (int i = 0; i <= m->num_rows; i++) {
			p_n[i] = p_n_plus[i];
			r_n[i] = r_n_plus[i];
			x_n[i] = x_n_plus[i];
		}
		//printf("%lf\n", x_n[0]);

		double abs = 0;
		for (int i = 0; i <= m->num_rows; i++) {
		abs += fabs(r_n[i])*fabs(r_n[i]);
		}
		L2_norm = sqrt(abs);
	}
	//printf("%lf\n", r_n[m->num_rows]);
	printf("\n%lf", L2_norm);
	printf("%lf\n", x_n[0]);
	return 0;
}