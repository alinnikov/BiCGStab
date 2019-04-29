#pragma once
#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"
#include "Preconditioners.h"

#pragma once
#ifndef BISGSTAB	
#define BISGSTAB
int BiCGStab(struct CSR_matrix *Matrix, double *b, double *x_n, double eps, int max_iterations);
#endif