#include "Multiplicate.h"


int multiplicate(struct CSR_matrix *Matrix, double *x, double *b) {

	for (int i = 0; i < Matrix->num_rows; i++)
	{
		b[i] = 0.0;
		for (int j = Matrix->array_rows[i]; j < Matrix->array_rows[i + 1]; j++) {
			//printf("%d\n", Matrix->array_columns[j-1]);
			b[i] += Matrix->array_values[j] * x[Matrix->array_columns[j]];
			
		}
		//printf("%d\n", i);
		//printf("%f\n", b[i]);
	}
	//printf("%f\n", b[Matrix->num_rows]);
	return 0;
}