#include "GetMatrix.h"

int main(void)
{
	struct CSR_matrix m;
	get_csr_matrix(&m);

	float *x = (float*)malloc((m.num_rows+1) * sizeof(float));
	float *b = (float*)malloc((m.num_rows+1) * sizeof(float));


	for (int i = 0; i < m.num_rows; i++) {
		x[i] = 1.0;
	}

	multiplicate(&m, x, b);

	//printf("%d", m.num_rows);
	//getchar();
	return 0;
}