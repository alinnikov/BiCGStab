#include "GetMatrix.h"
#include "multiplicate.h"
#include <time.h>
#include <dos.h> 


int main(void)
{

	#ifdef _OPENMP
		printf("OpenMP is supported!\n");
	#endif 

	struct CSR_matrix A; 

	int count = 1;

	double tol = 1e-12;
	int max_iter = 10000;
	int m = 60;
	char name[1024];
	
	printf("Print name of matrix: ");
	gets(name);
	
	read_csr_matrix(&A, name);

	double *x = (double*)malloc((A.num_rows) * sizeof(double));
	double *x_ref = (double*)malloc((A.num_rows) * sizeof(double));
	double *rhs;

	unsigned long start = GetTickCount();
	for (int i = 0; i < count; i++) {

		for (int i = 0; i < A.num_rows; i++) {
			x_ref[i] = 1.0;
			x[i] = 0.0;
		}

		rhs = spmv(A, x_ref);

		//GMRes(&A, rhs, x, tol, max_iter, m);
		BiCGStab(&A, rhs, x, tol, max_iter);
	
	}
	unsigned long delta = GetTickCount() - start;
	printf("used %lu milliseconds in average\n", delta/count);
	
	/*
	for (int i = 0; i < m.num_rows; i++) {
		printf("%e\n", x[i] - x_ref[i]);
		getchar();
	}
	*/
	
	getchar();
	return 0;
}