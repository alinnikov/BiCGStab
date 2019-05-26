#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"
#include <omp.h>
#pragma once
#ifndef MULTIPLICATE	
#define MULTIPLICATE
int multiplicate(struct CSR_matrix Matrix, double *x, double *b);
double dot_product(double *a, double *b, int num);
double *spmv(struct CSR_matrix Matrix, double *x);
#endif
