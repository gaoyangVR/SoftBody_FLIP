#include "mathmatrix.h"

/** Return the one norm of the matrix.
*/
float oneNorm(matrix3x3 A)
{
	const float sum1 = fabs(A.x00) + fabs(A.x10) + fabs(A.x20);
	const float sum2 = fabs(A.x01) + fabs(A.x11) + fabs(A.x21);
	const float sum3 = fabs(A.x02) + fabs(A.x12) + fabs(A.x22);
	float maxSum = sum1;
	if (sum2 > maxSum)
		maxSum = sum2;
	if (sum3 > maxSum)
		maxSum = sum3;
	return maxSum;

}

/** Return the inf norm of the matrix.
*/
float infNorm(matrix3x3 A)
{
	const float sum1 = fabs(A.x00) + fabs(A.x01) + fabs(A.x02);
	const float sum2 = fabs(A.x10) + fabs(A.x11) + fabs(A.x12);
	const float sum3 = fabs(A.x20) + fabs(A.x21) + fabs(A.x22);
	float maxSum = sum1;
	if (sum2 > maxSum)
		maxSum = sum2;
	if (sum3 > maxSum)
		maxSum = sum3;
	return maxSum;

}

void jacobiRotate(matrix3x3 &A, matrix3x3 &R, int p, int q)
{
	// rotates A through phi in pq-plane to set A(p,q) = 0
	// rotation stored in R whose columns are eigenvectors of A
	if (Apq(A,p,q)==0.0)
		return;

	float d = (Apq(A, p, p) - Apq(A, q, q)) / (2.0*Apq(A,p, q));
	float t = 1.0 / (fabs(d) + sqrt(d*d + 1.0));
	if (d < 0.0) t = -t;
	float c = 1.0 / sqrt(t*t + 1);
	float s = t*c;
	float tmp=0.;
	tmp = Apq(A, p, p)+ t*Apq(A, p, q);
	setApq(A, p, p,tmp);
	tmp= Apq(A, q, q)-t*Apq(A, p, q);
	setApq(A, q, q, tmp);
	setApq(A, p, q, 0);setApq(A, q, p, 0);
	// transform A
	int k;
	for (k = 0; k < 3; k++) {
		if (k != p && k != q) {
			float Akp = c*Apq(A,k, p) + s*Apq(A,k, q);
			float Akq = -s*Apq(A,k, p) + c*Apq(A,k, q);
			setApq(A, k, p, Akp); setApq(A, p, k, Akp);
			setApq(A, k, q, Akq); setApq(A, q, k,Akq);
		}
	}
	// store rotation in R
	for (k = 0; k < 3; k++) {
		float Rkp = c*Apq(R,k, p) + s*Apq(R,k, q);
		float Rkq = -s*Apq(R, k, p) + c*Apq(R,k, q);
		setApq(R,k, p, Rkp);
		setApq(R,k, q, Rkq);
	}
}

void eigenDecomposition(matrix3x3 &A, matrix3x3 &eigenVecs, float3 &eigenVals)
{
	const int numJacobiIterations = 10;
	const float epsilon = 1e-15;

	matrix3x3 D = A;

	// only for symmetric matrices!
	eigenVecs.x00 = eigenVecs.x01 = eigenVecs.x02 = eigenVecs.x10 = eigenVecs.x11 = eigenVecs.x12 = eigenVecs.x20 = eigenVecs.x21 = eigenVecs.x22 = 0.; eigenVecs.x00 = eigenVecs.x11 = eigenVecs.x22 = 1.;	// unit matrix
	int iter = 0;
	while (iter < numJacobiIterations) {	// 3 off diagonal elements
		// find off diagonal element with maximum modulus
		int p, q;
		float a, max;
		max = fabs(D.x01);
		p = 0; q = 1;
		a = fabs(D.x02);
		if (a > max) { p = 0; q = 2; max = a; }
		a = fabs(D.x12);
		if (a > max) { p = 1; q = 2; max = a; }
		// all small enough -> done
		if (max < epsilon) break;
		// rotate matrix with respect to that element
		jacobiRotate(D, eigenVecs, p, q);
		iter++;
	}
	eigenVals.x = D.x00;
	eigenVals.y = D.x11;
	eigenVals.z = D.x22;

}

float Apq(matrix3x3 A,int p, int q)
{
	float tmp;
	if (p == 0 && q == 0) tmp = A.x00;
	if (p == 0 && q == 1) tmp = A.x01;
	if (p == 0 && q == 2) tmp = A.x02;
	if (p == 1 && q == 0) tmp = A.x10;
	if (p == 1 && q == 1) tmp = A.x11;
	if (p == 1 && q == 2) tmp = A.x12;
	if (p == 2 && q == 0) tmp = A.x20;
	if (p == 2 && q == 1) tmp = A.x21;
	if (p == 2 && q == 2) tmp = A.x22;

	return tmp;

}

void setApq(matrix3x3 &A, int p, int q, float tmp)
{
	if (p == 0 && q == 0) A.x00 = tmp;
	if (p == 0 && q == 1)  A.x01 = tmp;
	if (p == 0 && q == 2)  A.x02 = tmp;
	if (p == 1 && q == 0)  A.x10 = tmp;
	if (p == 1 && q == 1)  A.x11 = tmp;
	if (p == 1 && q == 2)  A.x12 = tmp;
	if (p == 2 && q == 0)  A.x20 = tmp;
	if (p == 2 && q == 1)  A.x21 = tmp;
	if (p == 2 && q == 2)  A.x22 = tmp;

}
