#include "Multiplicate.h"


int multiplicate(struct CSR_matrix *Matrix, double *x, double *b) {

	for (int i = 0; i <= Matrix->num_rows; i++)
	{
		b[i] = 0.0;
		for (int j = Matrix->array_rows[i]; j < Matrix->array_rows[i + 1]; j++) {

			b[i] += Matrix->array_values[j] * x[Matrix->array_columns[j]];

		}
	}
	return 0;
}

double scalar_multiplicate(double *a, double *b, int num) {
	double result=0;
	
	for (int i=0; i < num; i++) {
		result += a[i] * b[i];
	}
	return result;
}