#include <stdio.h>
#include <math.h>

int ndim = 2, osci = 7;
//ndim is the number of elements in y
double t0 = 0, tf = 10, dt = 0.01, accuracy;
//t0 is initial time, tf is final time, dt is the timestep
double x0 = 1, v0 = 0;
//x0 is the intial x (position) value, v0 is the intial v (velocity) value
double k = 4, m = 1, w, Eint;
//(For oscillators) k is the spring constant, m is mass, and w is omega
//Eint is the initial energy
double y[2], dydx[2];
/*y[] stores the x value in y[0], and the v value in y[1];
dydx[] stores the derivatives of y[], with their respective indecies*/
double Eapp[3];
//Eapp[] stores the energy of the approximation
//Exxx[0] = PE, Exxx[1] = KE, Exxx[2] = Total E
double  t;
//t is time

void Dydx(double [], double [], int*);
/*Dydx is the derivative function. Despite the name, the function actually
calculates the derivative at t = t+ dt/2, per the Second-Order Runga-Kutta Method*/
void Step(double, double [], int, double, double [],
	void(*derivs)(double [], double [], int*));
//Step uses the dydx[] array to move forward a full dt step in time
void energy(double [], double [], int*);
//Calculates all the neede values for energy
double absval(double);
//Calculates the absolute value of a number

int main(){
	t = t0; //initializes the time
	w = pow((k/m),0.5); //defines omega
	y[0] = x0; //initializes x
	y[1] = v0; //initializes v
	energy(Eapp,y,&osci); //initializes energy
	Eint = Eapp[2];
	accuracy = (1-absval(1-absval(Eapp[2]/Eint)))*100; //calculates initial accuracy
	FILE *printer;
	printer = fopen("RungaOscillatorsOsci7.txt", "w");
	printf("(t,ycalc):(%f,%f) accuracy = %f\n",t, y[0], accuracy);
	fprintf(printer, "%f\t%f\t%f\t%f\t%f\n",t, y[0], y[1], x0, -x0);
	
	while(t < tf){
		Step(t,y,ndim,dt,dydx,Dydx);
		t += dt;
		/*Since Step moves the model ahead dt units of time, he exact values
		should be calculated at t+dt*/
		energy(Eapp,y,&osci);
		accuracy = (1-absval(1-absval(Eapp[2]/Eint)))*100;
		printf("(t,ycalc):(%f,%f) accuracy = %f\n",t, y[0], accuracy);
		fprintf(printer, "%f\t%f\t%f\t%f\t%f\n",t, y[0], y[1], x0, -x0);
	}
	fclose(printer);
	return 0;
}

void Dydx(double y[], double dydx[], int* osci){
	double a, yt0, yt1, at;
	//Calculate the acceleration at t = t
	a = (-1)*(pow(w,2))*pow(y[0],(2**osci-1));
	//Calculate Runga-Kutta Intermediates (t = t + dt/2)
	yt0 = y[0]+((dt/2)*y[1]);
	yt1 = y[1]+((dt/2)*a);
	//Calculate dv/dt for Runga-Kutta (acceleration)
	at = (-1)*(pow(w,2))*pow(yt0,(2**osci-1));
	//Set corresponding Runga-Kutta intermediates into dydx array
	dydx[0] = yt1;
	dydx[1] = at;
}

void Step(double t, double y[], int ndim, double dt, double dydx[],
	void(*derivs)(double [], double [], int*)){
	//Compute final Runga-Kutta step
	int j;
	(*derivs)(y,dydx,&osci);
	for(j = 0 ; j < ndim ; j++) y[j] += dydx[j]*dt ;
}

void energy(double Eapp[], double y[], int* osci){
	Eapp[0] = (1.0/(2.0**osci))*k*(pow(y[0],(2**osci)));
	Eapp[1] = (0.5)*m*(pow(y[1],2));
	Eapp[2] = Eapp[0]+Eapp[1];
}

double absval(double input){
	double returner = (sqrt(pow(input,2)));
	return returner;
}
