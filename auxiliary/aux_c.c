#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#define M_PI acos(-1.0)



SEXP bignum2func(SEXP yx, SEXP yy, SEXP cx, SEXP cy, SEXP h){
	
	int nC = length(cx);
	int nY = length(yy);
	double hh = asReal(h);

	double *pyx, *pyy, *pcx, *pcy;
	pyx = REAL(yx);
	pyy = REAL(yy);
	pcx = REAL(cx);
	pcy = REAL(cy);
	double bignum = 0;
	double bignum2 = 0;
	double temp1, temp2, t1, t2;
	double h2 = 2*hh*hh;
	double h3 = -nY*log((M_PI*h2));

	for(int i=0; i<nY; i++){
		t1 = pyx[i];
		t2 = pyy[i];
		for(int j=0; j<nC; j++){
			temp1= t1 - pcx[j];
			temp2= t2 - pcy[j];	
			bignum += exp(-(temp1*temp1 + temp2*temp2)/(h2));
		}
		bignum2 += log(bignum);
		bignum = 0;

	}
	bignum2 = h3+bignum2;

	SEXP out = PROTECT(allocVector(REALSXP,1));
	REAL(out)[0] = bignum2;
	UNPROTECT(1);

	return out;
}



SEXP gethsamp(SEXP yx, SEXP yy, SEXP cx, SEXP cy){
	int nC = length(cx);
	int nY = length(yy);

	double *pyx, *pyy, *pcx, *pcy, *pout;
	pyx = REAL(yx);
	pyy = REAL(yy);
	pcx = REAL(cx);
	pcy = REAL(cy);

	SEXP out = PROTECT(allocVector(REALSXP,nY));
	pout = REAL(out);


	double temp1, temp2, temp3;
	double temp = LONG_MAX;

	for(int i=0; i<nY; i++){
		for(int j=0; j<nC; j++){
			temp1 = pyx[i] - pcx[j];
			temp2 = pyy[i] - pcy[j];
			temp3 = temp1*temp1 + temp2*temp2;
			if (temp3 < temp){
				temp = temp3;
			}
		}

		pout[i] = sqrt(temp);
		temp = LONG_MAX;
		
	}

	UNPROTECT(1);

	return out;
}