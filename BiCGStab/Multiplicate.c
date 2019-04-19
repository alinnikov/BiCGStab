#include "Multiplicate.h"



int multiplicate(struct CSR_matrix Matrix, double *x, double *b) {
#pragma omp parallel for num_threads(1)

	for (int i = 0; i < Matrix.num_rows; i++)
	{
		b[i] = 0.0;
		for (int j = Matrix.array_rows[i]; j < Matrix.array_rows[i + 1]; j++) {

			b[i] += Matrix.array_values[j] * x[Matrix.array_columns[j]];

		}
	}
	return 0;
}

// spmv
double *spnv_pointer(struct CSR_matrix Matrix, double *x) {
	double *b = calloc((Matrix.num_rows), sizeof(double));
#pragma omp parallel for num_threads(1)
	for (int i = 0; i < Matrix.num_rows; i++)
	{
		for (int j = Matrix.array_rows[i]; j < Matrix.array_rows[i + 1]; j++) {

			b[i] += Matrix.array_values[j] * x[Matrix.array_columns[j]];

		}
	}
	return b;
}


// dot product
double dot_product(double *a, double *b, int num) {
	double result=0;
#pragma omp parallel shared(result)  num_threads(2)
#pragma omp for reduction(+:result)
	for (int i = 0; i < num; i++) {
		result += a[i] * b[i];
	}
	
	return result;
}