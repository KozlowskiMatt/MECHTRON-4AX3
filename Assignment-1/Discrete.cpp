#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <cmath>

// Discrete representation

using namespace std;

int main()
{
    double h, t, total_time, velocity, g;
    int alpha = 60; //initial launch angle

    Eigen::MatrixXd A(4,4); 
    Eigen::MatrixXd C(2,4);
    Eigen::MatrixXd k(4,1);
    Eigen::MatrixXd X (4,1); // State vector of pos_x, vel_x, pos_y, vel_y
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
    position_file<<"\tTime\tPos_x\tPos_y"<<endl;
    position_file<<"\t"<<t<<"\t"<<Y(0,0)<<"\t"<<Y(1,0)<<endl;

    while (t<total_time)
    {
        X = A*X +k;
        Y = C*X;
        myfile<<"\t"<<t<<"\t"<<X(0,0)<<"\t"<<X(1,0)<<"\t"<<X(2,0)<<"\t"<<X(3,0)<<endl;
        t = t+h; //increment t by step size
        position_file<<"\t"<<t<<"\t"<<Y(0,0)<<"\t"<<Y(1,0)<<endl;
    }
    return 0;
}