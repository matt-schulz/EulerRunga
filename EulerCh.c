#include <stdio.h>
#include <math.h>

double accel = 4;
double t0 = 0, tf = 100, dt = 0.1, xint = 16, vint = -7.3, accuracy;
double k = 4, m = 1, w;
double y[2];
double yexact[2];
double dydx[2];
int ndim = sizeof(y);
double  t;

void Dydx(double, double [], double []);
void Step(double, double [], double [], int, double,
    void (*derivs)(double, double [], double []));
double function(double);
double absval(double);

int main(){
	t = t0;
	w = pow((k/m),0.5);
	y[0] = xint;
	y[1] = vint;
	yexact[0] = xint + vint*(t-t0) + 0.5*accel*(t-t0)*(t-t0);
	yexact[1] = vint + accel*(t-t0) ;
	FILE *printer;
	printer = fopen("EulerCh.txt", "w");
	
	while(t <= tf){
		accuracy = (1-absval(1-absval(yexact[0]/y[0])))*100;
		printf("(t,ycalc):(%f,%f), yex = %f, accuracy = %f\n",t, y[0], yexact[0], accuracy);
		fprintf(printer, "%f\t%f\t%f\t%f\t%f\n",t, y[0], y[1], yexact[0], yexact[1]);
		t += dt;
		Step(t,y,dydx,ndim,dt,Dydx);
		accel = 4;//(-1)*(pow(w,2))*y[0];
		yexact[0] = xint + vint*(t-t0) + 0.5*accel*(t-t0)*(t-t0);
		yexact[1] = vint + accel*(t-t0) ;
//		yexact[0] = cos(w*t);
//		yexact[1] = (-1)*w*sin(w*t);
	}
	fclose(printer);
	return 0;
}

void Dydx(double x, double y[],double dydx[])
{
	int j ;
	dydx[0] = y[1] ;
	dydx[1] = accel ;
}

void Step(double x, double y[], double dydx[], int n, double dx,
	void (*derivs)(double, double [], double [])){
		int j ;
		(*derivs)(x,y,dydx) ;
		for(j = 0 ; j < n ; j++) y[j] += dydx[j]*dx ;
}

double absval(double input){
	double returner = (sqrt(pow(input,2)));
	return returner;
}
