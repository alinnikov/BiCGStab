#include "GetMatrix.h"
#include "multiplicate.h"
#include <time.h>
#include <dos.h> 
#include "malloc.h"
#include <stdio.h>  
#include <io.h> 

#pragma once
#ifndef PRECONDITIONERS	
#define PRECONDITIONERS	
int ILU0(struct CSR_matrix *Matrix, double *luvalues);
void GaussSolve(struct CSR_matrix *Matrix, double* b, double* result);
#endif
