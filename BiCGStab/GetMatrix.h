#include <stdio.h>  
#include <io.h> 
#include <malloc.h>
#pragma once
#ifndef GET_MATRIX
#define GET_MATRIX

struct CSR_matrix {
	int num_rows;
	int num_values;
	int *array_rows;
	int *array_columns;
	double *array_values;
};
int read_csr_matrix(struct CSR_matrix *m, char *name);
int get_num_rows(char *name);
int get_num_values(char *name);
int create_csr_matrix(struct CSR_matrix *m, int *rowind, int *col, double *val);
#endif

