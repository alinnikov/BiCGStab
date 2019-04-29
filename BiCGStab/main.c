#include "GetMatrix.h"
#include "multiplicate.h"
#include <time.h>
#include <dos.h> 
#include "Preconditioners.h"
#include <stdio.h>  

int main(void)
{

#ifdef _OPENMP
	printf("OpenMP is supported!\n");
#endif 

	

	struct CSR_matrix m; 
	struct CSR_matrix LU;
	int count = 1;
	double eps = 0.000000000001;
	int max_iterations = 10000;
	char name[1024];
	printf("Print name of matrix\n");
	gets(name);
	read_csr_matrix(&m, name);

	create_csr_matrix(&LU, m.array_rows, m.array_columns, m.array_values);

	double *x_n = (double*)malloc((m.num_rows) * sizeof(double));
	double *x = (double*)malloc((m.num_rows) * sizeof(double));
	double *b;
	double *lu_values = (double*)malloc((m.num_values) * sizeof(double));
	double *result = (double*)malloc((m.num_rows) * sizeof(double));

	for (int i = 0; i < m.num_rows; i++) {
		x[i] = 9.0;
		x_n[i] = 10.0;
	}

	unsigned long start = GetTickCount();
	for (int i = 0; i < count; i++) {

		for (int i = 0; i < m.num_rows; i++) {
			x[i] = 9.0;
			x_n[i] = 10.0;
		}

		b = spnv_pointer(m, x);

		BiCGStab(&m, b, x_n, eps, max_iterations);
	
	}

	unsigned long delta = GetTickCount() - start;
	printf("used %lu milliseconds in average\n", delta/count);
	getchar();
	return 0;
}