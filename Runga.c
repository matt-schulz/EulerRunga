#include <stdio.h>
#include <math.h>

double t0 = 0, tf = 100, dt = 0.01, xint = 1, vint = 0, accuracy, accel, accelt;
double PEapp, KEapp, PEexa, KEexa, Eapp, Eexa;
double k = 4, m = 3, w;
double y[2];
double yt[2];
double yexact[2];
int ndim = sizeof(y);
double  t;

void Dydx(double, double [], double []);
void Step(double, double [], double [], int, double,
    void (*derivs)(double, double [], double []));
void Runga(double [], double [], double, double);
void accelerator(double*, double*, double [], double []);
void functio(double [], double);
void energy(double*, double*, double*, double*, double*, double*, double [], double []);
double absval(double);

int main(){
	t = t0;
	w = pow((k/m),0.5);
	y[0] = xint;
	y[1] = vint;
	functio(yexact,t);
	Runga(y,yt,dt,accel);
	accelerator(&accel,&accelt,y,yt);
	FILE *printer;
	printer = fopen("Runga.txt", "w");
	
	while(t <= tf){
		energy(&PEapp,&KEapp,&PEexa,&KEexa,&Eapp,&Eexa,y,yexact);
		accuracy = (1-absval(1-absval(Eapp/Eexa)))*100;
		printf("(t,ycalc):(%f,%f), yex = %f, accuracy = %f\n",t, y[0], yexact[0], accuracy);
		fprintf(printer, "%f\t%f\t%f\t%f\t%f\t%f\n",t, y[0], y[1], yexact[0], yexact[1], accuracy);
		t += dt;
		Runga(y,yt,dt,accel);
		accelerator(&accel,&accelt,y,yt);
		Step(t,y,yt,ndim,dt,Dydx);
		functio(yexact,t);
	}
	fclose(printer);
	return 0;
}

void accelerator(double* accel, double* accelt, double y[], double yt[]){
	*accel = (-1)*(pow(w,2))*y[0];
	*accelt = (-1)*(pow(w,2))*yt[0];
}

void functio(double yexact[], double t){
	yexact[0] = cos(w*t);
	yexact[1] = (-1)*w*sin(w*t);
}

void Runga(double y[], double yt[], double dt, double accel){
	yt[0] = y[0]+((dt/2)*y[1]);
	yt[1] = y[1]+((dt/2)*accel);
}

void Dydx(double x, double yt[],double dydx[]){
	int j ;
	dydx[0] = yt[1] ;
	dydx[1] = accelt ;
}

void energy(double* PEapp, double* KEapp, double* PEexa, double* KEexa, double* Eapp, double* Eexa, double y[], double yexact[]){
		*PEapp = (0.5)*k*(pow(y[0],2));
		*KEapp = (0.5)*m*(pow(y[1],2));
		*PEexa = (0.5)*k*(pow(yexact[0],2));
		*KEexa = (0.5)*m*(pow(yexact[1],2));
		*Eapp = *PEapp+*KEapp;
		*Eexa = *PEexa+*KEexa;
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
