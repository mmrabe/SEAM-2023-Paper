
#ifndef __GAMMALIKE__

#define __CGAMMALIKE__

#include <complex.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>

#include "gammalike.c"
#include "logsumexp.c"
#include "array.c"


double cgammaengine2(double x, double * a, double * b, int N, int mode) {
	if(N < 1) error(1, "Convolution must have at least one component!");
	else if(N == 1) return gammaloglike(x, a[1], b[1]);
	const int max_iter = 100000;
	if(x < 0.0) return -INFINITY;
	double b1 = b[1];
	int i, j;
	double rho = a[1];
	for(i = 2; i <= N; i++) {
		if(b1 < b[i]) b1 = b[i];
		rho += a[i];
	}
	double complex log_bigc = 0.0;
	for(i = 1; i <= N; i++) {
		log_bigc += a[i]*clog(b[i]/b1);
	}
	int k = 0;
	double complex new_v, log_deltak2;
	double complex * deltaks = vector(double complex, max_iter);
	double complex * v = vector(double complex, max_iter);
	for(k = 0; k == 0 || (k < max_iter && creal(new_v) >= log(DBL_EPSILON)) ; k++) {
		if(k == 0) {
			log_deltak2 = 0.0;
		} else {
			double complex * dvec = vector(double complex, k);
			for(i = 1; i <= k; i++) {
				double complex * log_gamma_k_i = vector(double complex, N);
				for(j = 1; j <= N; j++) {
					log_gamma_k_i[j] = clog(a[j]) + i * clog((double complex) (1.0-b[j]/b1)) - clog(i);
				}
				dvec[i] = clog(i) + clogsumexp(log_gamma_k_i, N);
				free_vector(double complex, log_gamma_k_i);
				if(i < k) dvec[i] += deltaks[k-i];
			}
			for(i=1;i<=k;i++) {
				//printf("dvec[k=%d,i=%d]=%lf+%lfi\n", k, i, creal(v[i]), cimag(v[i]));
			}
			log_deltak2 = deltaks[k] = -clog(k) + clogsumexp(dvec, k);
			free_vector(double complex, dvec);
		}
		double p;
		if(mode == 1) {
			p = gammaloglike(x,rho+k,b1);
		} else if(mode == 2) {
			p = gammalogcdf(x,rho+k,b1);
		}
		//printf("p(%g,%g,%g)=%g\n", x, rho+k, b1, p);
		new_v = log_deltak2 + p;
		//printf("v[%d]=%g+%gi + %g+%gi = %g+%gi\n", k+1, creal(log_deltak2), cimag(log_deltak2), creal(p), cimag(p), creal(new_v), cimag(new_v));
		v[k+1] = new_v;
	}
	free_vector(double complex, deltaks);
	double ret = creal(log_bigc + clogsumexp(v, k));
	free_vector(double complex, v);
	return ret;
}

double cgammaengine(double x, double * a, double * b, int N, int mode) {
	double * new_a = vector(double, N);
	double * new_b = vector(double, N);
	int n=0, i;
	for(i=1; i<=N; i++) {
		if(a[i] > DBL_EPSILON && b[i] > DBL_EPSILON && !isinf(b[i])) {
			n++;
			new_a[n] = a[i];
			new_b[n] = b[i];
		}
	}
	double ret = cgammaengine2(x, new_a, new_b, n, mode);
	free_vector(double, new_a);
	free_vector(double, new_b);
	return ret;
}

double cgammaloglike(double x, double * a, double * b, int N) {
	return cgammaengine(x, a, b, N, 1);
}

double cgammalike(double x, double * a, double * b, int N) {
	return exp(cgammaloglike(x, a, b, N));
}


double cgammalogpdf(double x, double * a, double * b, int N) {
	return cgammaengine(x, a, b, N, 2);
}

double cgammapdf(double x, double * a, double * b, int N) {
	return exp(cgammalogpdf(x, a, b, N));
}


#endif


