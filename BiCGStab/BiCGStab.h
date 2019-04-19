#pragma once
#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"

#pragma once
#ifndef BISGSTAB	
#define BISGSTAB
int BiCGStab(struct CSR_matrix *Matrix, double *b, double *x_n, double eps, int max_iterations);
//double alpha(struct CSR_matrix *Matrix, double *r_0_help, double *r_n, double *p_n);
#endif