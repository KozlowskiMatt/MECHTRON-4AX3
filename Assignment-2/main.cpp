#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

#define INPUT_DATA "DATA.txt"
#define OUTPUT_DATA "OUTPUT_4.txt"
#define SIZE_OF_DATA 2048

using namespace std;

/*
    comp_coeff calculates and returns the necessary coefficients for 
    a n degree polynomial using the LLS method. 
*/

Eigen::MatrixXd comp_coeff(Eigen::MatrixXd phi,Eigen::MatrixXd Y, int m, int n)
{
    Eigen::MatrixXd theta (n+1,1); // Vector that consists of coefficients
    theta = (phi.transpose()*phi).inverse() * phi.transpose() * Y;
    return theta;
}

/*
    coeff_filter applies the caluclated LLS coefficients and produces the filtered data
    point (Y_tilda)
*/
double coeff_filter(int m, int n,Eigen::MatrixXd phi, Eigen::MatrixXd theta,int value_wanted)
{
    Eigen::MatrixXd Y_tilda(m,1);
    Y_tilda = phi*theta;
    return Y_tilda(value_wanted,0);
}

int main()
{
    int m,n; //inputs from user
    int j,q,i,counter; //Used for indexing amongst for and while loops
    int m_maninpulated;
    
    double l; //To read data from the dataset
    double filtered_data_point; //To store a specific filtered data point

    cout << "Please give your desired frame length (m): ";
    cin >> m;
    cout << "Please give your desired power of polynomial (n): ";
    cin >> n;

    ifstream my_file;
    ofstream write_file;
    my_file.open(INPUT_DATA);
    write_file.open(OUTPUT_DATA);

    Eigen::MatrixXd Y (m,1); // size of vector should be (m,1)
    Eigen::MatrixXd phi(m,n+1);
    Eigen::MatrixXd theta(n+1,1);
    Eigen::MatrixXd temp(m,1); //Used as a temp variable to store filtered data

    // Get the first m data points from the dataset
    for (i=0;i<m;i++)
    {
        my_file >> l;
        Y(i,0) = l;
    }

    // Making the phi vector
    m_maninpulated = -(m-1)/2; // to establish bounds centered around 0, from [-(m-1)/2, (m-1)/2]  m= 5 --> [-2,-1,0,1,2]
    for (j = 0; j<m; j++)
    {
        for (q=0; q<=n; q++)
        {
            phi(j,q) = pow(m_maninpulated,q);
        }
        m_maninpulated ++;
    }
   // m_maninpulated = -(m-1)/2;

    theta = (phi.transpose()*phi).inverse() * phi.transpose() * Y; //   Compute the first set of coefficients
    temp = phi*theta; //The first set of m filtered data points
    j=0;
    for (i=(-(m-1)/2); i <1;i++) //Write the first (m+1)/2 (range --> [-(m-1)/2, 0]) fitlered data points into the text file
    {
        filtered_data_point = temp(j,0);
        write_file << filtered_data_point<<endl;
        j++;
    }
    counter = 0;
    temp = Y;
    while (my_file >> l)
    {
        temp(m-1,0) = l;
        for (i = 0; i<(m-1); i++)
        {
            temp(i,0) = Y(i+1,0);
        }
        Y = temp;
        theta = comp_coeff(phi,Y,m,n);
        filtered_data_point = coeff_filter(m,n,phi,theta,(m-1)/2); // Gather the 'x=0' term from the calculated matrix 
        write_file << filtered_data_point<<endl;
        if (SIZE_OF_DATA - counter == (m+1))
            break;
        counter++;
    }

    // To write the last set of data points once we have reached eof
    for (i=j; i<m;i++)
        write_file << coeff_filter(m,n,phi,theta,i)<<endl;

    cout<<"\n\nCOMPLETE"<<endl;
    return 0;
}