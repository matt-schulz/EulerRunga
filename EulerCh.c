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
	/*accel is the acceleration of the particle in question
	t0 is the intial time
	tf is the final time
	dt is the time step
	xint is the initial x (position) value
	vint is the initial v (velocity) value
	k is the spring constant, for oscillators
	m is the mass of the object in question
	w is omega, for oscillators
	y[] is the array that stores the intermediate x and v
		values in y[0] and y[1], respectively
	yexact[] is the array that stores the exact intermediate
		values for x and y (same indexing as y[])
	dydx[] is the array that stores the derivatives of y[]'s
		values (y[1] for y[0], and accel for y[1])
	ndim is the number of dimensions in y[]
	t is time*/

void Dydx(double, double [], double []);
	/*This function casts the dydx indecies to their
	respective values*/
void Step(double, double [], double [], int, double,
    void (*derivs)(double, double [], double []));
    /*This function completes Euler's Method*/
double absval(double);

int main(){
	t = t0; //setting t at the initial time
	w = pow((k/m),0.5); //calculating omega for oscillators
	y[0] = xint; //seting x initial
	y[1] = vint; //setting v initial
	yexact[0] = xint + vint*(t-t0) + 0.5*accel*(t-t0)*(t-t0); //exact x function
	yexact[1] = vint + accel*(t-t0) ; //exact v function
	FILE *printer;
	printer = fopen("EulerCh.txt", "w");
	
	while(t <= tf){
		accuracy = (1-absval(1-absval(yexact[0]/y[0])))*100; //method to calculate accuracy
		printf("(t,ycalc):(%f,%f), yex = %f, accuracy = %f\n",t, y[0], yexact[0], accuracy);
		fprintf(printer, "%f\t%f\t%f\t%f\t%f\n",t, y[0], y[1], yexact[0], yexact[1]);
		t += dt;
		Step(t,y,dydx,ndim,dt,Dydx); //completes the step in Eulers method
		accel = 4; //recasting of acceleration for jerking systems
		yexact[0] = xint + vint*(t-t0) + 0.5*accel*(t-t0)*(t-t0);
		yexact[1] = vint + accel*(t-t0) ;
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
