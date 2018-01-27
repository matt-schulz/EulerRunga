#include <stdio.h>
#include <math.h>

int ndim = 2;
double t0 = 0, tf = 10, dt, accuracy;
double x0 = 1, v0 = 0;
double k = 4, m = 1, w;
double y[2], yexact[2], dydx[2];
double Eapp[3], Eexa[3];
double  dts = 0.001, t, dd;

void Dydx(double [], double []);
void Step(double, double [], int, double, double [],
	void(*derivs)(double [], double []));
void functio(double [], double);
void energy(double [], double [], double [], double []);
double absval(double);

int main(){
	w = pow((k/m),0.5);
	FILE *printer;
	printer = fopen("RungaErrorTester.txt", "w");
	
	for(dt = 0.001; dt < (0.1+dts); dt += 0.001){
		t = t0;
		y[0] = x0;
		y[1] = v0;
		while(t < tf){
			Step(t,y,ndim,dt,dydx,Dydx);
			t += dt;
			functio(yexact,t);
			energy(Eapp,Eexa,y,yexact);
		}
		accuracy = absval(Eapp[2]-Eexa[2]);
		dd = pow(dt,(3.0));
		printf("dt = %.10f, accuracy = %f\n",dd, accuracy);
		fprintf(printer, "%.10f\t%f\n",dd, accuracy);
	}
	fclose(printer);
	return 0;
}

void Dydx(double y[], double dydx[]){
	double a, yt0, yt1, at;
	//Calculate the acceleration at t = t
	a = (-1)*(pow(w,2))*y[0];
	//Calculate Runga-Kutta Intermediates (t = t + dt/2)
	yt0 = y[0]+((dt/2)*y[1]);
	yt1 = y[1]+((dt/2)*a);
	//Calculate dv/dt (acceleration)
	at = (-1)*(pow(w,2))*yt0;
	//Set corresponding Runga-Kutta intermediates into dydx array
	dydx[0] = yt1;
	dydx[1] = at;
}

void Step(double t, double y[], int ndim, double dt, double dydx[],
	void(*derivs)(double [], double [])){
	//Compute final Runga-Kutta step
	int j;
	(*derivs)(y,dydx);
	for(j = 0 ; j < ndim ; j++) y[j] += dydx[j]*dt ;
}

void functio(double yexact[], double t){
	yexact[0] = x0*cos(w*t);
	yexact[1] = (-1)*x0*w*sin(w*t);
}

void energy(double Eapp [], double Eexa [], double y[], double yexact[]){
	Eapp[0] = (0.5)*k*(pow(y[0],2));
	Eapp[1] = (0.5)*m*(pow(y[1],2));
	Eexa[0] = (0.5)*k*(pow(yexact[0],2));
	Eexa[1] = (0.5)*m*(pow(yexact[1],2));
	Eapp[2] = Eapp[0]+Eapp[1];
	Eexa[2] = Eexa[0]+Eexa[1];
}

double absval(double input){
	double returner = (sqrt(pow(input,2)));
	return returner;
}
