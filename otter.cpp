#include "otter.h"
#include <Eigen/Dense>

using namespace std;

Otter::Otter(string controlSystem = "stepInput",
             float r = 0,
             float V_current = 0,
             float beta_current = 0,
             int tau_x = 120)
    // member init part
    : nu{0, 0, 0, 0, 0, 0}
    , u_actual{0, 0}
    , u_control{0, 0}
    , current_eta{0, 0, 0, 0, 0, 0}
    , speed(0)
    , C{{{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}}}
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
    Vector3f rg = (m*rg+mp*rp)/(m+mp);  // CG corrected for payload

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
    T = nabla / (2*Cb_pont*B_pont*L) // Draft 
    Matrix3f Ig_CG = m * 


    float Xudot = -0.1*m;
    float Yvdot;
    float Zwdot;
    float Kpdot;
    float Mqdot;
    float Nrdot;
}

