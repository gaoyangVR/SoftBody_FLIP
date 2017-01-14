#ifndef  MATHMATRIX_H
#define  MATHMATRIX_H

#include<vector_types.h>
#include "utility.h"
#include "mymesh.h"
#include<math.h>
#include <stdio.h>

float oneNorm(matrix3x3 m);
float infNorm(matrix3x3 m);

void eigenDecomposition(matrix3x3 &A, matrix3x3 &eigenVecs, float3 &eigenVals);
void jacobiRotate(matrix3x3 &A, matrix3x3 &R, int p, int q);

float Apq(matrix3x3 A,int p, int q);
void setApq(matrix3x3 &A, int p, int q, float tmp);

float vec9(int a, int b);

#endif // ! MATHMATRIX_H
