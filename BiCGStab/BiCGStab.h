#pragma once
#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"

#ifndef BISGSTAB	
#define BISGSTAB
int BiCGStab(struct CSR_matrix *Matrix, double *b, double *x_n, double tol, int max_iter);
//double alpha(struct CSR_matrix *Matrix, double *r_0_help, double *r_n, double *p_n);
#endif