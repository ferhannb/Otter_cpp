#include "otter.h"
#include "gnc.h"
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;


Otter::Otter(string controlSystem,
             float r,
             float V_current,
             float beta_current,
             int tau_x)
    // member init part
    : nu{{0, 0, 0, 0, 0, 0}}, u_actual{{0, 0}}, u_control{{0, 0}}, current_eta{{0, 0, 0, 0, 0, 0}}, speed(0), C(MatrixXf::Zero(6, 6)), ref(r), V_c(V_current), beta_c(beta_current), tauX(tau_x)

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
    float g = 9.81;  // acceleration of gravity (m/s^2)
    int rho = 1025;  // density of water
    float m = 55.0;  // mass (kg)
    float mp = 25.0; // payload (kg)
    m_total = m + mp;
    Vector3f rp = {0, 0, 0};            // location of payload (m)
    Vector3f rg = {0.2, 0, -0.2};       // CG for hull only (m)
    rg = (m * rg + mp * rp) / (m + mp); // CG corrected for payload

    S_rg = Smtrx(rg);
    H_rg = Hmtrx(rg);
    S_rp = Smtrx(rp);

    float R44 = 0.4 * B;
    float R55 = 0.25 * L;
    float R66 = 0.25 * L;

    float T_yaw = 1.0;       // time constant in yaw (s)
    float Umax = 6 * 0.5144; // max forward speed (m/s)

    // Data for one pontoon

    B_pont = 0.25;
    float y_pont = 0.395;
    float Cw_pont = 0.75;
    float Cb_pont = 0.4;

    // Inertia dyadic, volume displacement and draft

    float nabla = (m + mp) / rho;
    T = nabla / (2 * Cb_pont * B_pont * L); // Draft
    Matrix3f diagM;
    diagM << pow(R44, 2), 0, 0,
        0, pow(R55, 2), 0,
        0, 0, pow(R66, 2);

    Matrix3f Ig_CG = m * diagM;
    Ig = Ig_CG - m * S_rg * S_rg - mp * S_rg * S_rg;

    // Experimental propeller data including lever arms

    l1 = -y_pont;
    l2 = y_pont;
    k_pos = 0.02216 / 2;
    k_neg = 0.01289 / 2;
    n_max = sqrt((0.5 * 24.4 * g) / k_pos);
    n_min = sqrt((0.5 * 13.6 * g) / k_neg);

    Matrix3f Mident;
    Mident << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;

    MatrixXf MRB_CG = MatrixXf::Zero(6, 6);
    MRB_CG.block<3, 3>(0, 0) = (m + mp) * Mident;
    MRB_CG.block<3, 3>(3, 3) = Ig;
    MatrixXf MRB = H_rg.transpose() * MRB_CG * H_rg;

    // Hydrodynamic added mass (best practice)
    float Xudot = -0.1 * m;
    float Yvdot = -1.5 * m;
    float Zwdot = -1.0 * m;
    float Kpdot = -0.2 * Ig(0, 0);
    float Mqdot = -0.8 * Ig(1, 1);
    float Nrdot = -1.7 * Ig(2, 2);

    MA << -Xudot, 0, 0, 0, 0, 0,
        0, -Yvdot, 0, 0, 0, 0,
        0, 0, -Zwdot, 0, 0, 0,
        0, 0, 0, -Kpdot, 0, 0, 0,
        0, 0, 0, 0, -Mqdot, 0,
        0, 0, 0, 0, 0, -Nrdot;

    // System mass Matrix
    M = MRB + MA;
    Minv = M.inverse();

    // Hydrostatic quantities (Fossen 2021, Chapter 4)
    float Aw_pont = Cw_pont * L * B_pont;

    float I_T = (
            2
            * (1.0 / 12)
            * L
            * pow(B_pont, 3)
            * (6 * pow(Cw_pont,3) / ((1 + Cw_pont) * (1 + 2 * Cw_pont)))
            + 2 * Aw_pont * pow(y_pont, 2)
        );

    float I_L = 0.8 * 2 * (1.0 / 12) * B_pont * pow(L, 3);

    float KB = (1.0 / 3) * (5 * T / 2 - 0.5 * nabla / (L * B_pont));
    float BM_T = I_T / nabla; // BM values
    float BM_L = I_L / nabla;
    float KM_T = KB + BM_T;// KM values
    float KM_L = KB + BM_L;
    float KG = T - rg[2];
    float GM_T = KM_T - KG;// GM values
    float GM_L = KM_L - KG;

    float G33 = rho * g * (2 * Aw_pont);// spring stiffness
    float G44 = rho * g * nabla * GM_T;
    float G55 = rho * g * nabla * GM_L;
    // spring stiff. matrix in CF
    Matrix<float, 6,6> G_CF;
    G_CF << 0, 0,   0,   0,   0, 0,
            0, G33, 0,   0,   0, 0,
                    G44, 0,   0, 0,
                         G55, 0, 0,
                         0,   0, 0;
    float LCF = -0.2;
    MatrixXf H = Hmtrx(Vector3f {{LCF,0,0}});  // transform G_CF from CF to CO

    G = H.transpose() * G_CF * H;

    // Natural frequencies

    float w3 = sqrt(G33 / M(2,2));
    float w4 = sqrt(G44 / M(3,3));
    float w5 = sqrt(G55 / M(4,4));

    // Linear damping terms (hydroddynamic derivatives)

    float Xu = -24.4 * g / Umax;
    float Yv = 0;
    float Zw = -2 * 0.3 * w3 * M(2,2);
    float Kp = -2 * 0.2 * w4 * M(3,3);
    float Mq = -2 * 0.4 * w5 * M(4,4);
    float Nr = M(5,5) / T_yaw;

    
    D<<-Xu,0,0,0,0,0,
        0,-Yv,0,0,0,0,
        0,0,-Zw,0,0,0,
        0,0,0,-Kp,0,0,
        0,0,0,0,-Mq,0,
        0,0,0,0,0,-Nr;

    // Trim: theta = -7.5 deg corresponds to 13.5 cm less height aft

    trim_moment = 0;
    trim_setpoint = 280;

    // Propeller configuration/input matrix 
    
    Matrix2f B = Matrix2f{{1,1},{-l1,-l2}}*k_pos;
    Matrix2f Binv = B.inverse();

    // Heading autoplot

    e_int = 0; // integral state
    wn = 1.2; // PID pole placement 
    zeta = 0.8; 

    // Referance model 


    r_max = 10*M_PI/180; // maximum yaw rate
    psi_d = 0;
    r_d = 0;
    a_d = 0;
    wn_d = wn / 5;
    zeta_d = 1 ;

}


Otter::dynamic_func_output Otter::dynamics(VectorXf eta, VectorXf nu, Vector2f u_actual,Vector2f u_control, float SampleTime)

{
    // input vector
    Vector2f n {u_actual(0),u_actual(1)};

    // Current Velocities
    float u_c = V_c * cos(beta_c-eta(5)); // current surge vel.
    float v_c = V_c * sin(beta_c-eta(5)); // current sway vel.

    VectorXf nu_c = VectorXf {{u_c,v_c,0,0,0,0}};   // current velocity vector 
    VectorXf nu_r = nu - nu_c;                      // relative velocity vector

    // Rigid body and added mass Coriolis and centripetal matrices
    // CRB_CG = [ (m+mp) * Smtrx(nu2)          O3   (Fossen 2021, Chapter 6)
    //              O3                   -Smtrx(Ig*nu2)  ]

    MatrixXf CRB_CG = MatrixXf::Zero(6,6);
    CRB_CG.block<3,3>(0,0) = m_total * Smtrx(nu(seqN(0,3)));
    CRB_CG.block<3,3>(3,3) = -Smtrx(Ig*nu(seqN(3,3)));
    VectorXf CRB = H_rg.transpose() * CRB_CG * H_rg;  // transform CRB from CG to CO

    MatrixXf CA = m2c(MA,nu_r); 
    CA(5,0) = 0;// assume that yhe munk moment in yaw can be neglected 
    CA(5,1) = 0;// if nonzero, must be balanced by adding nonlinear damping 

    MatrixXf C = CRB + CA;

    // Ballast 
    VectorXf g_0 {{0,0,0,0,trim_moment,0}};

    // Control forces and moments - with propeller revolution saturation

    Vector2f thrust = Vector2f::Zero(2,0);
    for(int i=0;i<2;i++)
    {
        n[i] = sat(n(i),n_min,n_max);

        if (n(i)>0)
            {
            thrust(i) = k_pos * n(i) * abs(n(i));
            }
        else 
            {
            thrust(i) = k_neg * n(i) * abs(n(i));
            }
    }


    // Control forces and moments

    VectorXf tau {{thrust(0)+thrust(1),0,0,0,0,-l1 * thrust(0)-l2 * thrust(1)}};

    // Hydrodynamic linear damping + nonlinear yaw damping

    VectorXf tau_damp = -D*nu_r;
    tau_damp(5) = tau_damp(5) - 10 * D(5,5) * abs(nu_r(5) * nu_r(5));
    VectorXf tau_crossflow = crossFlowDrag(L,B_pont,T,nu_r);
    VectorXf sum_tau = tau + tau_damp + tau_crossflow - (C * nu_r) - (G * eta) - g_0;


    VectorXf nu_dot = Minv * sum_tau ; // USV dynamics
    Vector2f n_dot = (u_control - n) / T_n ; // propeller dynamics
    
    float trim_momment  = trim_moment + SampleTime * trim

}
