#include <stdio.h>  
#include <io.h> 
#include "GetMatrix.h"
#pragma once
#ifndef MULTIPLICATE	
#define MULTIPLICATE
int multiplicate(struct CSR_matrix *Matrix, float *x, float *b);
#endif
