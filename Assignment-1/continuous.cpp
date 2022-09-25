#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <cmath>
using namespace std;

/*
    Continuous simuilation of a projectile using Runge-kutta technique
*/

Eigen::VectorXd calculate(double t,const Eigen::VectorXd &in)
{
	Eigen::VectorXd out;
    const double g = -9.8; //Acceleration due to gravity

	Eigen::MatrixXd A(4, 4);
    Eigen::MatrixXd k(4, 1);
	A << 0,1,0,0,
         0,0,0,0,
         0,0,0,1,
         0,0,0,0;
    
    k << 0,0,0,g;

	out= A*in + k; // Implementation of state space equation X'(t) = A*x(t) +Bu(t) + k || B = 0
	return(out);
}

// Implemenation of the Runge-Kutta technique
Eigen::VectorXd Runge(double t, const Eigen::VectorXd &in, double h)
{
        Eigen::VectorXd k1, k2, k3, k4;
        k1 = calculate(t, in);
        k2 = calculate(t + h / 2, in + (h / 2) * k1);
        k3 = calculate(t + h / 2, in + (h / 2) * k2);
        k4 = calculate(t + h, in + h * k3);

        return (in + h* (k1 + 2 * k2 + 2 * k3 + k4) / 6);
}

void solver(const Eigen::VectorXd &x0,double h,double time)
{
	ofstream myfile;
	myfile.open ("Continuous_Data");

	double  t; // Keep track of time in while loop
	t=0;
	Eigen::VectorXd x=x0;
    myfile<<"\t"<<t<<"\t"<<x(0,0)<<"\t"<<x(1,0)<<"\t"<<x(2,0)<<"\t"<<x(3,0)<<endl;
    t=t+h;
	while(t<time){
		x=Runge(t,x,h);
        myfile<<"\t"<<t<<"\t"<<x(0,0)<<"\t"<<x(1,0)<<"\t"<<x(2,0)<<"\t"<<x(3,0)<<endl; //Formatting 
		t=t+h;
	}
}

int main(void)
{
	Eigen::VectorXd X(4,1);
    double h, velocity, total_time, t,alpha;
    velocity = 20; // Initial velocity = 20 m/s
    alpha = 60; // Launch angle = 60 deg
    h = 0.01;
    total_time = (velocity*sin(alpha*M_PI/180))/4.9; // Time projectile is in air

	X << 0,velocity*cos(alpha*M_PI/180),0,velocity*sin(alpha*M_PI/180); //Initial vector [pos_x,vel_x,pos_y,vel_y]
    solver(X,h,total_time);
    
    return 0;
}