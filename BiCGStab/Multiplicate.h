#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"
#include <omp.h>
#pragma once
#ifndef MULTIPLICATE	
#define MULTIPLICATE
int multiplicate(struct CSR_matrix Matrix, double *x, double *b);
double scalar_multiplicate(double *a, double *b, int num);
#endif
