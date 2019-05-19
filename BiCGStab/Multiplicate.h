#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"
#include <omp.h>
#pragma once
#ifndef MULTIPLICATE	
#define MULTIPLICATE
int multiplicate(struct CSR_matrix Matrix, double *x, double *b);
double dot_product(double *a, double *b, int num);
double *spnv_pointer(struct CSR_matrix Matrix, double *x);
void Gauss(double *Hess, double* b, double* x, int n, int m);
void MatrMultiply(int n, int m, double *matrix, double *vektor, double *res);
#endif
