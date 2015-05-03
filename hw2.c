/*
Author: Seoyun Lee						Class: CPSC 445
Homework Assignment: 2

Uses steepest descent algorithm to find the solution to a linear system
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps);
void matvec(double *matrix, double *vector, int n, double *result);
void multiply(double *a, double *b, int n, double *ab);
void subtract(double *left, double *right, int n, double *result);
void transpose(double *matrix, int n, double *trans);
double dotproduct(double *vec1, double *vec2, int n);
double step(double *a, double *axb, int n);
void gradient(double *a, double *axy, int n, double *g);

int main() {
	return 0;
}

void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps) {
	double *xn;
	xn = malloc(n * sizeof(double));

	double *ax;
	ax = malloc(n * sizeof(double));

	double *axy;
	axy = malloc(n * sizeof(double));

	double *g;
	g = malloc(n * sizeof(double));

	double stepsize;
	stepsize = 0.0;

	/* x0 = 0 */
	int i;
	for(i = 0; i < n; i++) {
		xn[i] = 0.0;
	}

	double approx;
	approx = 0.0;


	for(*niter = 0; *niter < numit; (*niter)++) {
		matvec(a, xn, n, ax);
		subtract(ax, y, n, axy);

		approx = dotproduct(axy, axy, n);

		discreps[*niter] = approx;

		if(approx < eps) {
			break;
		}

		gradient(a, axy, n, g);
		stepsize = step(a, axy, n);

		for(i = 0; i < n; i++) {
			g[i] = stepsize * g[i];
		}
		subtract(xn, g, n, xn);
	}

	for(i = 0; i < n; i++) {
		x[i] = xn[i];
	}


	free(g);
	free(axy);
	free(ax);
	free(xn);
}

/*
	multiplies a n x n matrix with a n vector
*/
void matvec(double *matrix, double *vector, int n, double *result) {
	int i;
	int j;

	double sum;

	for(i = 0; i < n; i++) {
		sum = 0.0;
		for(j = 0; j < n; j++) {
			sum = sum + (matrix[i * n + j] * vector[j]);
		}
		result[i] = sum;
	}
}

/*
	multiplies two matrices together
*/
void multiply(double *a, double *b, int n, double *ab) {
	int i;
	int j;
	int k;
	double sum;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			sum = 0.0;
			for(k = 0; k < n; k++) {
				sum = sum + a[i * n + k] * b[k * n + j];
			}
			ab[i * n + j] = sum;
		}
	}
}

/*
	left vector - right vector = result vector
*/
void subtract(double *left, double *right, int n, double *result) {
	int i;

	for(i = 0; i < n; i++) {
		result[i] = left[i] - right[i];
	}
}

/*
	transposes an n x n matrix
*/
void transpose(double *matrix, int n, double *trans) {
	int i;
	int j;

	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			trans[i * n + j] = matrix[j * n + i];
		}
	}
}

/*
	returns the dot product of 2 vectors
*/
double dotproduct(double *vec1, double *vec2, int n) {
	double result;
	result = 0.0;

	int i;

	for(i = 0; i < n; i++) {
		result = result + vec1[i] * vec2[i];
	}

	return result;
}

/*
	returns step size for steepest descent
*/
double step(double *a, double *axy, int n) {

	double *at;
	at = malloc(n * n * sizeof(double));
	/* at now contains the transpose of a */
	transpose(a, n, at);

	double *ataxy;
	ataxy = malloc(n * sizeof(double));
	/* multiply at by axb and store in ataxb */
	matvec(at, axy, n, ataxy);

	double num;
	num = dotproduct(ataxy, ataxy, n);

	double *aat;
	aat = malloc(n * n * sizeof(double));
	/* multiply a and at */
	multiply(a, at, n, aat);

	double *aataxy;
	aataxy = malloc(n * sizeof(double));
	/* a * at * (ax - b) */
	matvec(aat, axy, n, aataxy);

	double denom;
	denom = dotproduct(aataxy, aataxy, n);

	free(aataxy);
	free(aat);
	free(ataxy);
	free(at);

	double result;
	result = (num / (2.0 * denom));

	return result;
}

/*
	finds gradient
*/
void gradient(double *a, double *axy, int n, double *g) {

	double *at;
	at = malloc(n * n * sizeof(double));
	/* at now contains the transpose of a */
	transpose(a, n, at);

	int i;
	for(i = 0; i < (n * n); i ++) {
		at[i] = 2 * at[i];
	}

	/* calculate 2 * at * ax-y and store in gradient */
	matvec(at, axy, n, g);

	free(at);
}



