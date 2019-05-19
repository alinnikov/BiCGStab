#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"

#pragma once
#ifndef FGMRES
#define FGMRES
int GMRes(struct CSR_matrix *m, double *b, double *x_n, double eps, int max_iterations);

#endif