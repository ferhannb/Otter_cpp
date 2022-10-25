#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

Matrix3f  Smtrx(const Vector3f &a);

MatrixXf  Hmtrx(const VectorXf &a);

VectorXf attitudeEuler(VectorXf &eta,VectorXf &nu,float sampleTime);

Matrix3f Rzyx( float phi, float theta, float psi);

Matrix3f Tzyx(float phi, float theta);

// Matrix3f m2c( MatrixXf &M, VectorXf &eta);

float Hoerner (const float B, float T); 

MatrixXf crossFlowDrag(const float L, const float B, float T, VectorXf nu_r);

