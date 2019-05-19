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
	int k;
	double h_tilda;
	double h;
	double *z_j = (double*)malloc((m->num_rows) * sizeof(double));
	double *omega = (double*)malloc((m->num_rows) * sizeof(double));
	
	double *H = (double*)malloc(number * (number + 1) * sizeof(double));

	double Hess[10][11];
	double anti_Hess[11][10];
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 11; j++) {
			Hess[i][j] = 0.0;
			anti_Hess[i][j] = 0.0;
			H[i*number + j] = 0.0;

		}
	}


	double *anti_Hess_vector = (double*)malloc(number * sizeof(double));
	for (int i = 0; i < number; i++) {
		anti_Hess_vector[i] = 0.0;
	}
	
	double *V = (double*)malloc((m->num_rows) * number * sizeof(double));
	double *Z = (double*)malloc((m->num_rows) * number * sizeof(double));
	double *G = (double*)malloc((m->num_rows) * number * sizeof(double));
	double *v_j = (double*)malloc((m->num_rows) * sizeof(double));


	//double *y_m = (double*)malloc(number * sizeof(double));
	double *Ax_n = (double*)malloc((m->num_rows) * sizeof(double));
	double *c = (double*)malloc(number * sizeof(double));
	double *s = (double*)malloc(number * sizeof(double));
	
	double *g = (double*)malloc(number * sizeof(double));
	for (int i = 0; i < number; i++) {
		g[i] = 0.0;
	}

	for (int i = 0; i < number; i++) {
		c[i] = 1;
		s[i] = 1;
	}


	
Start:
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

	for (int j = 0; j < number; j++) {
		
		for (int i = 0; i < m->num_rows; i++) {

			//Получается слишком маленькое Z
				V[j*m->num_rows+i] = w[i] * antibeta;
				//Without Preconditioner
				Z[j*m->num_rows+i] = V[j*m->num_rows + i];
			}

		w = spnv_pointer(*m, &Z[j*m->num_rows]);

		for (int i = 0; i <= j; i++) {
			Hess[i][j] = dot_product(w, &V[i*m->num_rows], m->num_rows);
			for (int d = 0; d < m->num_rows; d++) {
				w[d] = w[d] + Hess[i][j] * V[i*m->num_rows + d];
			}
		}

		for (int i = 0; i <= j - 1; i++) {

			h_tilda = Hess[i][j];
			h = Hess[i+1][j];
			Hess[j][i] = c[i] * h_tilda + s[i] * h;
			Hess[j][i+1] = -s[i] * h_tilda + c[i] * h;
		}

		h_tilda = Hess[j][j];
		beta = sqrt(dot_product(w, w, m->num_rows));

		Hess[j][j] = sqrt(h_tilda*h_tilda + beta * beta);
		
		
		for (int l = 0; l < number; l++) {
			H[j*number + l] = Hess[j][l];
		}



		c[j] = beta / Hess[j][j];
		s[j] = h_tilda / Hess[j][j];

		g[j] = gamma * c[j];
		
		gamma = -gamma * s[j];

		if (gamma*gamma < eps*eps) {
			k = j;
			j = number + 1;
		}

	}

	//Вспомогательный вектор для общей суммы произведения двух матриц на g_k
	double *summ = (double*)malloc((m->num_rows) * sizeof(double));
	//Получаем вектор из умножения обратной матрицы Hess на g
	Gauss(H, g, anti_Hess_vector, number, number);
	//Находим произведение матрицы Z на вектор anti_Hess_vector. Значение summ[i] для каждого члена слишком мало
	MatrMultiply(number, m->num_rows, Z, anti_Hess_vector, summ);
	for (int i = 0; i < m->num_rows; i++) {
		x_n[i] = x_n[i] + summ[i];
	}

	printf("x[num_rows-1] = %.40lf\n", x_n[0]);

	if (gamma*gamma >= eps * eps) {
		goto Start;
	}


}