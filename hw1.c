/*
Author: Seoyun Lee						Class: CPSC 445
Homework Assignment: 1

Finds the SVD of an n x n matrix using the Jacobi algorithm.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
	Finds alpha/beta in the AT*A matrix
	*a: pointer to matrix A
	col: column which you are computing for
	n: size of matrix A
*/
double alb(double *a, int col, int n) {
	double al;
	int k;

	al = 0;
	for(k = 0; k < n; k++) {
		al = al + a[k * n + col] * a[k * n + col];
	}

	return al;
}

/*
	Finds gamma in the AT*A matrix
	*a: pointer to matrix A
	i: column i
	j: column j
	n: size of matrix A
*/
double gamma1(double *a, int i, int j, int n) {
	double g;
	int k;

	g = 0;
	for(k = 0; k < n; k++) {
		g = g + a[k * n + i] * a[k * n + j];
	}

	return g;
}

/*
	z: number to find the signum of
*/
int signum(double z) {
	if(z > 0) {
		return 1;
	}
	else if (z == 0) {
		return 0;
	}
	else {
		return -1;
	}
}

/*
	After doing the roatations of AT*A, we update U and V in the following way
	*matrix: pointer to matrix U or V
	n: size of matrix
	i: col i
	j: col j
	c: calculated cosine
	s: calculated sine
*/
void updatematrix(double *matrix, int n, int i, int j, double c, double s) {
	int k;
	double t;

	for(k = 0; k < n; k++) {
		t = matrix[k * n + i];
		matrix[k * n + i] = c * t - s * matrix[k * n + j];
		matrix[k * n + j] = s * t + c * matrix[k * n + j];
	}

}

/*
	Makes an n x n identity matrix
	*id: pointer to n x n user allocated memory
	n: size of matrix
*/
void identity(double *id, int n) {
	int i;
	int j;

	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			if(i == j) {
				id[i * n + j] = 1;
			}
			else {
				id[i * n + j] = 0;
			}
		}
	}
}

/*
	Finds the magnitude of the column i in matrix A
	*a: pointer to matrix 
	i: column i
	n: size of matrix
*/
double mag(double *a, int i, int n) {
	double norm2;

	norm2 = 0;
	int k;

	for(k = 0; k < n; k++) {
		norm2 = norm2 + a[k * n + i] * a[k * n + i];
	}

	double norm;
	norm = sqrt(norm2);

	return norm;
}

/*
	a: number you want to find the absolute value of
*/
double absolute(double a) {
	if(a >= 0) {
		return a;
	}
	else {
		return -a;
	}
}

/*
	Makes the U matrix by normalizing each column
	*a: pointer to matrix A
	n: size of matrix
	*u: pointer to user allocated matrix U
*/
void makeu(double *a, int n, double *u) {
	int i;
	int j;

	double norm;

	for(j = 0; j < n; j++) {
		norm = mag(a, j, n);
		for(i = 0; i < n; i++) {
			if(norm == 0.0) {
				u[i * n + j] = a[i * n + j];
			}
			else {
				u[i * n + j] = a[i * n + j] / norm;
			}
		}
	}
}

/*
	Makes the matrix S which is the magnitude of the cols of A
	*a: pointer to matrix A
	n: size of matrix
	*s: pointer to user allocated length n array S
*/
void makes(double *a, int n, double *s) {
	int i;
	int j;

	for(j = 0; j < n; j++) {
		double norm;
		norm = mag(a, j, n);
		s[j] = norm;
	}
}

/*
	Multiplies two matrices together, not used in jacobi alg
	*a: pointer to matrix A
	*b: pointer to matrix B
	n: size of matrix
	*ab: pointer to n x n space for A x B
*/
void multiply(double *a, double *b, int n, double *ab) {
	int i;
	int j;
	int k;
	double sum;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			sum = 0;
			for(k = 0; k < n; k++) {
				sum = sum + a[i * n + k] * b[k * n + j];
			}
			ab[i * n + j] = sum;
		}
	}
}

/*
	Finds the transpose of a matrix A, not used in algorithm
	*a: pointer to matrix A
	*at: pointer to space for AT
	n: size of matrix
*/
void transpose(double *a, double *at, int n) {
	int i;
	int j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			at[i * n + j] = a[j * n + i];
		}
	}
}

/*
	Swaps columns i and i+1 in matrix U
	*u: pointer to matrix
	n: size of matrix
	i: column i
*/
void swapcols(double *u, int n, int i) {
	double temp;

	int j;

	for(j = 0; j < n; j++) {
		temp = u[j * n + i];
		u[j * n + i] = u[j * n + (i + 1)];
		u[j * n + (i + 1)] = temp; 
	}
}

/*
	Uses bubble sort to sort S into decreasing order and
	swaps columns of U and V in the process
	*u: pointer to U
	*s: pointer to S
	*v: pointer to V
	n: size of matrix
*/
void sort(double *u, double *s, double *v, int n) {
	double swap;
	int c;
	int i;

	for(c = 0; c < (n - 1); c++) {
		for(i = 0; i < n - c - 1; i++) {
			if(s[i] < s[(i + 1)]) {
				swap = s[i];
				s[i] = s[(i + 1)];
				s[(i + 1)] = swap;
				swapcols(u, n, i);
				swapcols(v, n, i);
			}
		}
	}
}

/*
	The jacobi algorithm
	*a: pointer to the matrix A that we use the algorithm on
	n: size of A (n x n)
	*s: pointer to user allocated memory for array S which contains
		the singular values of A
	*u: pointer to user allocated memory for matrix U
	*v: pointer to user allocated memory for matrix V
*/
void jacobi(double *a, int n, double *s, double *u, double *v) {
	/* start v as the identity matrix */
	identity(v, n);

	double converge;
	converge = 0;

	int i;
	int j;

	while (converge != 1.0) {
		converge = 1.0;
		for(j = 1; j < n; j++) {
			for(i = 0; i < j; i++) {
				double al;
				double b;
				double g;

				al = alb(a, i, n);
				b = alb(a, j, n);
				g = gamma1(a, i, j, n);

				double zeta;
				double t;
				double c;
				double s;

				zeta = (b - al)/(2*g);
				t = signum(zeta) / (absolute(zeta) + sqrt(1 + zeta * zeta));
				c = 1 / sqrt(1 + t * t);
				s = c * t;

				converge = c;

				updatematrix(a, n, i, j, c, s);
				updatematrix(v, n, i, j, c, s);
			}
		}
	}


	makeu(a, n, u);
	makes(a, n, s);
	sort(u, s, v, n);
}

int main() {
	return 0;
}
