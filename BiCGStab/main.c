#include "GetMatrix.h"

int main(void)
{
	struct CSR_matrix m;
	char name[] = "matrix.rb";
	get_csr_matrix(&m, name);

	double *x = (double*)malloc((m.num_rows+1) * sizeof(double));
	double *b = (double*)malloc((m.num_rows+1) * sizeof(double));


	for (int i = 0; i < m.num_rows; i++) {
		x[i] = 1.0;
	}

	multiplicate(&m, x, b);

	//printf("%d", m.num_rows);
	//getchar();
	return 0;
}