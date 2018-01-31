#include <stdio.h>
#include <math.h>

int ndim = 3;
//ndim is the number of differentiable elements in y
double t0 = 0, tf = 50, dt = 0.001, accuracy;
//t0 is initial time, tf is final time, dt is the timestep
double sigma = 10, beta = 8/3, rho = 13;
//lorentz constants
double xx0 = 1, xy0 = 1, xz0 = 1;
//initializes the poisition values
double y[6], dydx[6];
/*y[] stores the position values in y[0]-y[2], and the velocity
values in y[3]-y[5]; dydx[] stores the derivatives of y[], with
their respective indecies*/
double  t;
//t is time

void Dydx(double [], double []);
//Dydx is the derivative function
void Step(double, double [], int, double, double [],
	void(*derivs)(double [], double []));
//Step uses the dydx[] array to move forward a full dt step in time
double absval(double);
//Calculates the absolute value of a number

int main(){
	t = t0; //initializes the time
	y[0] = xx0; //initializes x
	y[1] = xy0; //initializes y
	y[2] = xz0; //initializes z
	FILE *printer;
	printer = fopen("RungaLorentz.txt", "w");
	//printf("t: %f | x = %.5f | y = %.5f | z = %.5f\n",t, y[0], y[1], y[2]);
	fprintf(printer, "%f\t%f\t%f\t%f\n",t, y[0], y[1], y[2]);
	
	while(t < tf){
		Step(t,y,ndim,dt,dydx,Dydx);
		t += dt;
		//printf("t: %f | x = %.5f | y = %.5f | z = %.5f\n",t, y[0], y[1], y[2]);
		fprintf(printer, "%f\t%f\t%f\t%f\n",t, y[0], y[1], y[2]);
	}
	fclose(printer);
	return 0;
}

void Dydx(double y[], double dydx[]){
	y[3] = sigma*(y[1]-y[0]);
	y[4] = y[0]*(rho-y[2])-y[1];
	y[5] = (y[0]*y[1])-beta*y[2];
	dydx[0] = y[3];
	dydx[1] = y[4];
	dydx[2] = y[5];
}

void Step(double t, double y[], int ndim, double dt, double dydx[],
	void(*derivs)(double [], double [])){
	//Compute final step
	int j;
	(*derivs)(y,dydx);
	for(j = 0 ; j < ndim ; j++) y[j] += dydx[j]*dt ;
}

double absval(double input){
	double returner = (sqrt(pow(input,2)));
	return returner;
}
