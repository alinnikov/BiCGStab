#include "GetMatrix.h"
#include "multiplicate.h"
#include <time.h>
#include <dos.h> 


int main(void)
{

#ifdef _OPENMP
	printf("OpenMP is supported!\n");
#endif 

	unsigned long start = GetTickCount();

	struct CSR_matrix m; 

	int count = 1;

	double eps = 0.00000001;
	int max_iterations = 3000000;

	char name[] = "ex35.rb";
	read_csr_matrix(&m, name);

	double *x_n = (double*)malloc((m.num_rows) * sizeof(double));
	double *x = (double*)malloc((m.num_rows) * sizeof(double));
	double *b;

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