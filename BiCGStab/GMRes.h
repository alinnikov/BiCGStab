#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"

#pragma once
#ifndef FGMRES
#define FGMRES
int GMRes(struct CSR_matrix *A, double *b, double *x, double tol, int max_iter, int m);

#endif