#include "Multiplicate.h"


int multiplicate(struct CSR_matrix *Matrix, float *x, float *b) {

	for (int i = 0; i < Matrix->num_rows; i++)
	{
		b[i] = 0.0;
		for (int j = Matrix->array_rows[i]; j < Matrix->array_rows[i + 1]; j++) {
			//printf("%d\n", Matrix->array_columns[j-1]);
			b[i] += Matrix->array_values[j-1] * x[Matrix->array_columns[j-1]-1];
			
		}
	//printf("%f\n", b[i]);
	}

	return 0;
}