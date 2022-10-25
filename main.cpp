
#include <iostream>
#include <Eigen/Dense>
#include "gnc.h"
using namespace std;
using namespace Eigen;

int main()
{
    
    MatrixXi A = MatrixXi::Random(3,3);
    Matrix3i B = Matrix3i::Zero(3,3);
    MatrixXi M  {{1,2,3,4,5,6},{2,4,6,8,10,12},{1,3,5,7,9,11},{12,13,14,15,16,17},{21,22,23,24,25,26},{33,35,36,37,45,67}};
    MatrixXf Ma = MatrixXf::Ones(6,6);
  
    Vector3f w (1.0,2.0,3.0);
    VectorXf nu(6);
    nu<< 1.,
     23., 4., 5., 6., 4.; 
    VectorXf eta(6);
    eta<<0.0,0.0,0.0,0.0,0.0,0.0;
    float phi {0};
    float theta {0};
    float psi {0};
    float sampleTime =0.02;
    cout<<eta<<endl;
    VectorXf T = attitudeEuler(eta,nu,sampleTime);
    cout<<"-----"<<endl;
    cout<<M.block<2,2>(0,0)<<endl;
 


    return 0;
}
