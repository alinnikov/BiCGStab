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

	int count = 100;

	char name[] = "matrix.rb";
	read_csr_matrix(&m, name);

	double *x = (double*)malloc((m.num_rows) * sizeof(double));
	double *b;

	for (int i = 0; i < count; i++) {

		for (int i = 0; i < m.num_rows; i++) {
			x[i] = 3.0;
		}
		//printf("%lf", x[m.num_rows]);


		b = spnv_pointer(m, x);

		BiCGStab(&m, b);
	
	}
	unsigned long delta = GetTickCount() - start;
	printf("used %lu milliseconds in average\n", delta/count);

	return 0;
}