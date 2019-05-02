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


int BiCGStab(struct CSR_matrix *m, double *b, double *x_n, double eps, int max_iterations) {
	//Структура с матрицей M
	struct CSR_matrix LU;
	LU.num_rows = m->num_rows;
	LU.num_values = m->num_values;
	LU.array_rows = (int*)malloc((m->num_rows + 1) * sizeof(int));
	LU.array_columns = (int*)malloc((m->num_values) * sizeof(int));
	LU.array_values = (double*)malloc((m->num_values) * sizeof(double));
	memcpy(LU.array_rows, m->array_rows, (m->num_rows+1) * sizeof(int));
	memcpy(LU.array_columns, m->array_columns, m->num_values * sizeof(int));
	memcpy(LU.array_values, m->array_values, m->num_values * sizeof(double));
	//Переменные для ILU
	double *lu_values = (double*)malloc((m->num_values) * sizeof(double));
	double *result = (double*)malloc((m->num_rows) * sizeof(double));
	double *y_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *z_n = (double*)malloc((m->num_rows) * sizeof(double));
	
	//ILU
	ILU0_fast(&LU, lu_values);
	memcpy(LU.array_values, lu_values, m->num_values * sizeof(double));

	//Создаём необходимые массивы
	double *r_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *r_0_help = (double*)malloc((m->num_rows) * sizeof(double));
	double *p_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *Ax_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *s_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *Az_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *Ay_n = (double*)malloc((m->num_rows) * sizeof(double));
	double alpha_n;
	double omega_n;
	double beta_n;
	double r0help_rn;
	double L2_norm = 1.0;
	int number_of_iterations = 0;
	
	// Задаём A*x0
	Ax_n = spnv_pointer(*m, x_n);
	
	//Задаём r0*, r0 и p0
	for (int i = 0; i < m->num_rows; i++) {
		r_n[i] = b[i] - Ax_n[i];
		r_0_help[i] = r_n[i];
		p_n[i] = r_n[i];
	}

	//Считаем диагональный вектор
	int *iptr = (int*)malloc((m->num_rows) * sizeof(int));
	int iptr_count = 0;
	for (int i = 0; i < m->num_rows; i++) {
		for (int j = m->array_rows[i]; j < m->array_rows[i + 1]; j++) {
			if (i == m->array_columns[j]) {
				iptr[iptr_count] = j;
				iptr_count++;
			}

		}
	}


	//Цикл
	for (int num = 0; num < max_iterations && L2_norm>eps * eps; num++) {
		GaussSolve(&LU, y_n, p_n, iptr);


		Ay_n = spnv_pointer(*m, y_n);

		r0help_rn = dot_product(r_0_help, r_n, m->num_rows);

		alpha_n = alpha(m, r_0_help, r0help_rn, y_n, Ay_n);

		s(m, r_n, p_n, alpha_n, s_n, Ay_n);

		GaussSolve(&LU, z_n, s_n,iptr);

		Az_n = spnv_pointer(*m, z_n);

		omega_n = omega(m, s_n, Az_n);

		for (int i = 0; i < m->num_rows; i++) {
			x_n[i] = x_n[i] + alpha_n * y_n[i] + omega_n * z_n[i];
		}


		for (int i = 0; i < m->num_rows; i++) {
			r_n[i] = s_n[i] - omega_n * Az_n[i];
		}
	
		beta_n = beta(m, r0help_rn, r_0_help, r_n, alpha_n, omega_n);

		for (int i = 0; i < m->num_rows; i++) {
			p_n[i] = r_n[i] - beta_n * omega_n * Ay_n[i] + beta_n * p_n[i];
		}

		if (number_of_iterations % 100 == 0) {
			printf("%.40lf\n", x_n[m->num_rows-1]);
			printf("L2_norm=%.40lf\n", L2_norm);
		}

		L2_norm = dot_product(r_n,r_n,m->num_rows);
		number_of_iterations++;
		
		free(Ay_n);
		free(Az_n);
	}
	printf("Number of iterations = %d\n", number_of_iterations);
	printf("L2_norm=%.60lf\n", sqrt(L2_norm));
	printf("x[num_rows-1] = %.40lf\n", x_n[m->num_rows-1]);
	return 0;
}