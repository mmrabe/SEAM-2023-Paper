
#include<stdio.h>
#include<math.h>

double integral_trapezoidal(double f(double x, void * params), void * params, double a, double b, int n){
	double h, sum=0;
	int i;
	h = fabs(b-a) / n;
	for(i=1;i<n;i++){
		sum += f(a+i*h);
	}
	return (h/2)*(f(a)+f(b)+2*sum);
}


