#include <stdio.h>
#include <math.h>

int ndim = 2;
//ndim is the number of elements in y
double t0 = , tf = , dt = , accuracy;
//t0 is initial time, tf is final time, dt is the timestep
double x0 = , v0 =;
//x0 is the intial x (position) value, v0 is the intial v (velocity) value
double k = , m = , w;
//(For oscillators) k is the spring constant, m is mass, and w is omega
double y[2], yexact[2], dydx[2];
/*y[] stores the x value in y[0], and the v value in y[1];
yexact[] stores the exact x value in y[0], and the exact v value in y[1];
dydx[] stores the derivatives of y[], with their respective indecies*/
double Eapp[3], Eexa[3];
//Eapp[] stores the energy of the approximation
//Eexa[] stores the energy of the exact values
//Exxx[0] = PE, Exxx[1] = KE, Exxx[2] = Total E
double  t;
//t is time

void Dydx(double [], double []);
/*Dydx is the derivative function. Despite the name, the function actually
calculates the derivative at t = t+ dt/2, per the Second-Order Runga-Kutta Method*/
void Step(double, double [], int, double, double [],
	void(*derivs)(double [], double []));
//Step uses the dydx[] array to move forward a full dt step in time
void functio(double [], double);
//Functio calculates the exact values, to which the program is modeling
void energy(double [], double [], double [], double []);
//Calculates all the neede values for energy
double absval(double);
//Calculates the absolute value of a number

int main(){
	t = t0; //initializes the time
	w = pow((k/m),0.5); //defines omega
	y[0] = x0; //intializes x
	y[1] = v0; //intializes v
	functio(yexact,t); //calculates the exact values which the program is trying to model
	accuracy = (1-absval(1-absval(Eapp[2]/Eexa[2])))*100; //calculates initial accuracy
	FILE *printer;
	printer = fopen("Runga.txt", "w");
	printf("(t,ycalc):(%f,%f), yex = %f, accuracy = %f\n",t, y[0], yexact[0], accuracy);
	fprintf(printer, "%f\t%f\t%f\t%f\t%f\n",t, y[0], y[1], yexact[0], yexact[1]);
	
	while(t < tf){
		Step(t,y,ndim,dt,dydx,Dydx);
		t += dt;
		/*Since Step moves the model ahead dt units of time, he exact values
		should be calculated at t+dt*/
		functio(yexact,t);
		energy(Eapp,Eexa,y,yexact);
		accuracy = (1-absval(1-absval(Eapp[2]/Eexa[2])))*100;
		printf("(t,ycalc):(%f,%f), yex = %f, accuracy = %f\n",t, y[0], yexact[0], accuracy);
		fprintf(printer, "%f\t%f\t%f\t%f\t%f\n",t, y[0], y[1], yexact[0], yexact[1]);
	}
	fclose(printer);
	return 0;
}

void Dydx(double y[], double dydx[]){
	double a, yt0, yt1, at;
	//Calculate the acceleration at t = t
	a = ;
	//Calculate Runga-Kutta Intermediates (t = t + dt/2)
	yt0 = y[0]+((dt/2)*y[1]);
	yt1 = y[1]+((dt/2)*a);
	//Calculate dv/dt for Runga-Kutta (acceleration)
	at = ;
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
	yexact[0] = ;
	yexact[1] = ;
}

void energy(double Eapp [], double Eexa [], double y[], double yexact[]){
	Eapp[0] = ;
	Eapp[1] = (0.5)*m*(pow(y[1],2));
	Eexa[0] = ;
	Eexa[1] = (0.5)*m*(pow(yexact[1],2));
	Eapp[2] = Eapp[0]+Eapp[1];
	Eexa[2] = Eexa[0]+Eexa[1];
}

double absval(double input){
	double returner = (sqrt(pow(input,2)));
	return returner;
}
