#include <stdio.h>  
#include <io.h> 
#pragma once
#ifndef GET_MATRIX
#define GET_MATRIX

struct CSR_matrix {
	int num_rows;
	int num_values;
	int *array_rows;
	int *array_columns;
	float *array_values;
};
void get_csr_matrix(CSR_matrix);
int get_num_rows();
int get_num_values();
#endif

