#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <cmath>

using namespace std;

// Continuous representation using Runge-kutta technique

double g = -9.8;

Eigen::VectorXd u_func(double t,const Eigen::VectorXd &in)
{

	Eigen::VectorXd out;

	return(out);
}

Eigen::VectorXd m_func(double t,const Eigen::VectorXd &in)
{
	Eigen::VectorXd out;


//	Eigen::VectorXd u(1);
	Eigen::MatrixXd A(4, 4);
    Eigen::MatrixXd k(4, 1);
	A << 0,1,0,0,
         0,0,0,0,
         0,0,0,1,
         0,0,0,0;
    
    k << 0,0,0,g;

	out= A*in + k;
	return(out);
}

Eigen::VectorXd run(double t, const Eigen::VectorXd &in, double h)
{
        Eigen::VectorXd k1, k2, k3, k4;
        k1 = m_func(t, in);
        k2 = m_func(t + h / 2, in + (h / 2) * k1);
        k3 = m_func(t + h / 2, in + (h / 2) * k2);
        k4 = m_func(t + h, in + h * k3);

        return in + h* (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

void solve(const Eigen::VectorXd &x0,double h,double time)
{
	ofstream myfile;
	myfile.open ("Continuous_Data");

	double  t;
	t=0;
	Eigen::VectorXd x=x0;
    myfile<<"\t"<<t<<"\t"<<x(0,0)<<"\t"<<x(1,0)<<"\t"<<x(2,0)<<"\t"<<x(3,0)<<endl;
    t=t+h;
	while(t<time){
		x=run(t,x,h);
        myfile<<"\t"<<t<<"\t"<<x(0,0)<<"\t"<<x(1,0)<<"\t"<<x(2,0)<<"\t"<<x(3,0)<<endl;
		t=t+h;
	}
}

int main(void)
{
        
	Eigen::VectorXd X(4,1);
    Eigen::VectorXd Y(4,1);
    double h, velocity, total_time, t,alpha;
    velocity = 20;
    alpha = 60;
    h = 0.01;
    total_time = (velocity*sin(alpha*M_PI/180))/4.9;
    

	X << 0,velocity*cos(alpha*M_PI/180),0,velocity*sin(alpha*M_PI/180);
    solve(X,h,total_time);
    return 0;
}