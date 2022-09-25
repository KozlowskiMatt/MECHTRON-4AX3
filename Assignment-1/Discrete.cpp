#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <cmath>

// Discrete representation of a projectile

using namespace std;

int main()
{
    double h, t, total_time, velocity, g;
    int alpha = 60; //initial launch angle

    Eigen::MatrixXd A(4,4); 
    Eigen::MatrixXd C(2,4);
    Eigen::MatrixXd k(4,1);
    Eigen::MatrixXd X (4,1); // State vector for pos_x, vel_x, pos_y, vel_y
    Eigen::MatrixXd Y (2,1); // To store appropriate x and y coordinates

    ofstream position_file;
    position_file.open("discrete_data_1");

    h = 0.01; // Initialize step size
    t = 0; // Keep track of time
    g = -9.8; // Acceleration due to gravity
    velocity = 20;

    A << 1,h,0,0,
         0,1,0,0,
         0,0,1,h,
         0,0,0,1;

    C << 1,0,0,0,
         0,0,1,0;

    k << 0,0,0,g*h;

    X << 0,velocity*cos(alpha*M_PI/180),0,velocity*sin(alpha*M_PI/180);
    Y <<0,0;

    total_time = (velocity*sin(alpha*M_PI/180))/4.9;
    while (t<total_time)
    {
        position_file<<"\t"<<t<<"\t"<<Y(0,0)<<"\t"<<Y(1,0)<<endl;
        X = A*X +k;
        Y = C*X;
        t = t+h; //increment t by step size
    }
    return 0;
}