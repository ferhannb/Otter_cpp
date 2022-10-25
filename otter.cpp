#include "otter.h"
#include "gnc.h"
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

Otter::Otter(string controlSystem ,
             float r,
             float V_current,
             float beta_current,
             int tau_x )
    // member init part
    : nu {{0, 0, 0, 0, 0, 0}}
    , u_actual {{0, 0}}
    , u_control{{0, 0}}
    , current_eta {{0, 0, 0, 0, 0, 0}}
    , speed(0)
    , C(MatrixXf::Zero(6,6))
    , ref(r)
    , V_c(V_current)
    , beta_c(beta_current)
    , tauX(tau_x)

{
    if (controlSystem == "headingAutopilot")
    {
        controlDescription = "Heading autopilot, psi_d = " + std::to_string(r) + " deg";
    }
    else
    {
        controlDescription = "Step inputs for n1 and n2";
        controlSystem = "stepInput";
    }

    controlMode = controlSystem;

    // controls[0] = "Left propeller shaft speed (rad/s)";
    // controls[1] = "Right propeller shaft speed (rad/s)";

    name = "Otter USV (see otter.py for more details)";
    // Initialize the Otter USV model
    T_n = 0.32;
    L = 2.0;
    B = 1.08;
    dimU = 2;
    float g = 9.81;      // acceleration of gravity (m/s^2)
    int rho = 1025;      // density of water
    float m = 55.0;      // mass (kg)
    float mp = 25.0;     // payload (kg)
    m_total = m+mp;
    Vector3f rp = {0,0,0};              // location of payload (m)
    Vector3f rg = {0.2,0,-0.2};         // CG for hull only (m)
    rg = (m*rg+mp*rp)/(m+mp);  // CG corrected for payload

    S_rg = Smtrx(rg);
    H_rg = Hmtrx(rg);
    S_rp = Smtrx(rp);
    
    float R44 = 0.4*B;
    float R55 = 0.25*L;
    float R66 = 0.25*L;

    float T_yaw = 1.0; // time constant in yaw (s) 
    float Umax = 6 * 0.5144; // max forward speed (m/s)

    // Data for one pontoon

    B_pont = 0.25;
    float y_pont = 0.395;
    float Cw_pont = 0.75;
    float Cb_pont = 0.4;

    // Inertia dyadic, volume displacement and draft

    float nabla = (m + mp ) / rho ;
    T = nabla / (2*Cb_pont*B_pont*L); // Draft 
    Matrix3f diagM;
    diagM << pow(R44,2),0,0,
                0,pow(R55,2),0,
                0,0,pow(R66,2);

    Matrix3f Ig_CG = m * diagM;
    Ig = Ig_CG - m* S_rg*S_rg - mp * S_rg * S_rg;
     
    // Experimental propeller data including lever arms

    l1 = -y_pont;
    l2 = y_pont;
    k_pos = 0.02216 / 2;
    k_neg = 0.01289 / 2;
    n_max = sqrt((0.5*24.4*g) / k_pos);
    n_min = sqrt((0.5 * 13.6 * g) / k_neg);


    Matrix3f Mident;
    Mident<< 1,0,0,
             0,1,0,
             0,0,1 ;

    MatrixXf MRB_CG = MatrixXf::Zero(6,6);
    MRB_CG.block<3,3>(0,0) = (m + mp) * Mident;
    MRB_CG.block<3,3>(3,3) = Ig;
    MatrixXf MRB = H_rg.transpose() * MRB_CG *H_rg ;

    // Hydrodynamic added mass (best practice)
    float Xudot = -0.1 * m;
    float Yvdot = -1.5 * m;
    float Zwdot = -1.0 * m;
    float Kpdot = -0.2 * Ig(0,0);
    float Mqdot = -0.8 * Ig(1,1);
    float Nrdot = -1.7 * Ig(2,2);

    MA << -Xudot,0,0,0,0,0,
            0,-Yvdot,0,0,0,0,
            0,0,-Zwdot,0,0,0,
            0,0,0,-Kpdot,0,0,0,
            0,0,0,0, -Mqdot,0,
            0,0,0,0,0,-Nrdot;


    // System mass Matrix

    M = MRB + MA ;
    
}

