#include "Preconditioners.h"

int ILU0(struct CSR_matrix *Matrix, double *lu_values) {
	
	int *iptr = (int*)malloc((Matrix->num_rows) * sizeof(int));
	int s1, s2;
	int y1, y2, y_end1, y_end2;
	int pr1, pr2;
	int k,b;


	//»щем диагональные элементы
	int iptr_count = 0;
	for (int i = 0; i < Matrix->num_rows; i++) {
		for (int j = Matrix->array_rows[i]; j < Matrix->array_rows[i + 1]; j++) {
			if (i == Matrix->array_columns[j]) {
				iptr[iptr_count] = j;
			//	printf("%d\n", j-Matrix->array_rows[i]);
				iptr_count++;
			}

		}
	}

	memcpy(lu_values, Matrix->array_values, Matrix->num_values * sizeof(double));

	for (int i = 1; i < Matrix->num_rows; i++) {
		s1 = Matrix->array_rows[i];
		pr1 = 1;
		while (pr1 == 1) {
			k = Matrix->array_columns[s1];
			if (k >= i) {
				break;
				//pr1 = 0;
			}
			lu_values[s1] = lu_values[s1] / lu_values[iptr[k]];
			//printf("After = %lf\n", lu_values[s1]);
			s2 = s1;
			s1++;
			y1 = s1;
			y_end1 = Matrix->array_rows[i + 1];
			y2 = iptr[k] + 1;
			y_end2 = Matrix->array_rows[k + 1];

			if ((y_end1 > y1) && (y_end2 > y2)) {
				continue;
			}
			pr2 = 1;
			while (pr2 == 1) {
				if (Matrix->array_columns[y1] == Matrix->array_columns[y2]) {
					lu_values[y1] = lu_values[y1] - lu_values[s2] * lu_values[y2];
					y1++;
					y2++;
				}
				if (Matrix->array_columns[y1] > Matrix->array_columns[y2]) {
					y2++;
				}
				if (Matrix->array_columns[y1] < Matrix->array_columns[y2]) {
					y1++;
				}
				if ((Matrix->array_rows[k+1]<y2) || (Matrix->array_rows[i + 1] < y2)) {
					pr2 = 0;
				}
			}



		}


	}
	
	return 0;
}



void GaussSolve(struct CSR_matrix *Matrix, double* y, double* z)
{
	int s1, s2,i;

	int *iptr = (int*)malloc((Matrix->num_rows) * sizeof(int));
	int iptr_count = 0;
	for (int i = 0; i < Matrix->num_rows; i++) {
		for (int j = Matrix->array_rows[i]; j < Matrix->array_rows[i + 1]; j++) {
			if (i == Matrix->array_columns[j]) {
				iptr[iptr_count] = j;
				iptr_count++;
			}

		}
	}



	for (int i = 0; i < Matrix->num_rows; i++) {
		s1 = Matrix->array_rows[i];
		s2 = Matrix->array_rows[i+1];
		y[i] = z[i];
		for (int k = s1; k < s2; k++) {
			if (Matrix->array_columns[k] < i) {
				y[i] = y[i] - Matrix->array_values[k] * y[Matrix->array_columns[k]];
				printf("Lower = %lf\n", y[i]);
			}

		}
	}

	for (int i = Matrix->num_rows-1; i >= 0; i--) {
		s1 = Matrix->array_rows[i];
		s2 = Matrix->array_rows[i + 1];
		for (int k = s1; k < s2; k++) {
			if (Matrix->array_columns[k] > i) {
				y[i] = y[i] - Matrix->array_values[k] * y[Matrix->array_columns[k]];
				printf("Upper = %lf\n", y[i]);
			}

		}
		y[i] = z[i] / Matrix->array_values[iptr[i]];
	}
	free(iptr);
}