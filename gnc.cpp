#include "gnc.h"
#include <stdexcept>
#include <iostream>

using namespace std;

Matrix3f  Smtrx(const Vector3f &a)
{

    Matrix3f S  {{0,-a[2],a[1]}
                ,{a[2],0,-a[0]}
                ,{-a[1],a[0],0}};
    return S;  
}

MatrixXf  Hmtrx(const VectorXf &a)
{
    MatrixXf H = MatrixXf::Identity(6,6);
    return H.block(0,3,3,3) = Smtrx(a).transpose();
}

Matrix3f Rzyx(float phi,float theta, float psi)
{
    float cphi = cos(phi);
    float sphi = cos(sphi);
    float cth = cos(theta);
    float sth = sin(theta);
    float cpsi = cos(psi);
    float spsi = sin(psi);

    Matrix3f R  {
                  { cpsi*cth, -spsi*cphi+cpsi*sth*sphi, spsi*sphi+cpsi*cphi*sth },
                  { spsi*cth,  cpsi*cphi+sphi*sth*spsi, -cpsi*sphi+sth*spsi*cphi },
                  { -sth,              cth*sphi,                 cth*cphi }
                };
 return R;
}

Matrix3f Tzyx( float phi, float theta)
{
    float cphi = cos(phi);
    float sphi = sin(phi);
    float cth  = cos(theta);
    float sth  = sin(theta);

    try
        {
               Matrix3f   T  {
            { 1,  sphi*sth/cth,  cphi*sth/cth },
            { 0,  cphi,          -sphi},
            { 0,  sphi/cth,      cphi/cth} };
            
            return T;
        }

    catch(...) {
            std::cout<<"Tzyx is singular for theta = +-90 degrees."<<std::endl;
               } 

}

VectorXf attitudeEuler(VectorXf &eta, VectorXf &nu,float sampleTime)

{   
    
    VectorXf p_dot = Rzyx(eta[3],eta[4],eta[5])*nu(seqN(0,3));
    VectorXf v_dot = Tzyx(eta[3],eta[4])*nu(seqN(3,3));
    // // Forward Euler Integration
    // eta(seqN(0,3)) = eta(seqN(0,3))+sampleTime*p_dot;
    // eta(seqN(3,3)) = eta(seqN(3,3))+sampleTime*v_dot;
    Vector3f FirstEta= eta(seqN(0,3))+sampleTime*p_dot;
    Vector3f LastEta = eta(seqN(3,3)) +sampleTime*v_dot;
    eta.block<3,1>(0,0) = FirstEta ;
    eta.block<3,1>(3,0) = LastEta ;
    cout  <<eta<<endl;
    return eta;
}

// Matrix3f m2c( MatrixXf &M, VectorXf &eta)

// {
//     M = 0.5 * (M + M.transpose()); // systematization of the inertia matrix
//     // float M11 = 
//     // float
//     // float
//     // float
//     return M 
// }


// float Hoerner (const float B, float T)
// {
//     Vector2f DATA1 {{ 0.0109,0.1766,0.3530,0.4519,0.4728,0.4929,0.4933,0.5585,0.6464,0.8336,0.9880,1.3081,1.6392,1.8600,2.3129,2.6000,3.0088,3.4508, 3.7379,4.0031 }};
//     Vector2f DATA2 {{1.9661,1.9657,1.8976,1.7872,1.5837,1.2786,1.2108,1.0836,0.9986,0.8796,0.8284,0.7599,0.6914,0.6571,0.6307,0.5962,0.5868,0.5859,0.5599,0.5593}};
    
//     return 0.0
// }

MatrixXf crossFlowDrag(const float L, const float B, float T, VectorXf nu_r)
{
    int rho = 1026;  // Density of water
    int n = 20;       // number of strip 
    
    float dx = L /20;
    // float Cd_2d = Hoerner(B,T); // 2D drag coefficient based on Hoerner's curve
    float Cd_2d = 1.004184470989761;
    float Yh = 0.0;
    float Nh = 0.0;
    float xL = -L/2;

    for(int i = 0;i<=n;i++)
    {
        float v_r = nu_r[1];  // relative sway velocity
        float r = nu_r[5];    // yaw rate
        float Ucf = abs(v_r + xL * r)*(v_r + xL * r );
        Yh = Yh - 0.5 * rho * T * Cd_2d * Ucf * dx ;  // sway force
        Nh =  Nh - 0.5 * rho * T * Cd_2d * xL * Ucf * dx ; // yaw moment 
        xL = xL + dx;   
    };

    VectorXf tau_crossflow(6);
    tau_crossflow<< 0,Yh,0,0,0,Nh;

    return tau_crossflow;
}