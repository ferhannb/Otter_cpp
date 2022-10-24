#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

using namespace Eigen;

class Otter
{
public:
    Otter(std::string controlSystem = "stepInput",
          float r = 0,
          float v_current = 0,
          float beta_current = 0,
          int tau_x = 120);

    void initialize_otter();

private:

    VectorXf nu = {0.0,0.0,0.0,0.0,0.0,0.0};
    VectorXf current_eta = {0.0,0.0,0.0,0.0,0.0,0.0};
    Vector2f u_actual = {0,0};
    Vector2f u_control = {0,0};
    int speed;
    std::string controlDescription;
    MatrixXf C;
    float ref;
    float V_c;
    float beta_c;
    int tauX;
    std::string controlMode;

    Vector2f controls={0,0};
    std::string name;
    
    float T_n;
    float L;
    float B;
    float m_total;
    int dimU;
    Matrix3f S_rg;
    MatrixXf H_rg;
    Matrix3f S_rp;

    float T_yaw; 
    float Umax;
    // Data for one pontoon
    float B_pont;

    
    // Inertia dyadic, volume displacement and draft

    float T;
    Matrix3f Ig;

    
    // Experimental propeller data including lever arms
    
    float l1;
    float l2;
    float k_pos;
    float k_neg;
    float n_max;
    float n_min;
    

    MatrixXf MA;
    MatrixXf D;
    int trim_moment = 0;
    int trim_setpoint = 280;

    Matrix2f Binv;

    // Heading autopilot
    
    int e_int;
    float wn ;
    float zeta ;

    // Reference model

    float r_max;
    float psi_d;
    float r_d;
    float a_d;
    float wn_d;
    float zeta_d;
};
