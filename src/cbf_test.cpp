#include <iostream>
#include <fstream>
#include <string>
#include <osqp.h>
#include "mpc/low_level_class.hpp"
#include "mpc/segway_example.hpp"
// Global declarations

using namespace std;
double offset_angle_ = 0.138324;
double dt_;
bool lowLevelActive_ = true;

using namespace ControlBarrierFunction;

CBF *cbf; // Define MPC object as a pointer
// state =[x, y, theta, v, thetaDot, psi, psiDot]

// Initialize constant parameters
const int nx = 7;
const int nu = 2;
double X[nx]     = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};     // true state
double Xn[nx]    = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};    // nominal state
double X_IC[nx]  = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};     // true state initial condition to CBF QP
double Xn_IC[nx] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};    // nominal state initial condition to CBF QP
double uMPC[nu]  = {0.0, 0.0};
double delay_ms_;
double Alinear[nx*nx] = {};
double Blinear[nx*nu] = {};
double Clinear[nx] = {};

double xmax        = {};
double ymax        = {};
double thetamax    = {};
double vmax        = {};
double thetaDotmax = {};
double psimax      = {};
double psiDotmax   = {};

int flagInputPub_ = 0;
int flagNominalStatePub_ = 0;
std::vector<double> inputTotBuffer1_;
std::vector<double> inputMpcBuffer1_;
std::vector<double> inputTotBuffer2_;
std::vector<double> inputMpcBuffer2_;

//// Main Definition
int main (int argc, char *argv[])
{
	double H_x[3] = {10000,10000,1};
	cbf = new CBF(nx_, nu_, true, H_x,x_eq_, fullDynamics,safetySet,clf);

	Xn[5] = offset_angle_;
	X[5]  = offset_angle_;

	cbf->setIC(X, Xn);

	cbf->setMatrices(Alinear, Blinear, Clinear);
	double uMPC[2] = {0.0,0.0};

	cbf->evaluateCBFqpConstraintMatrices(uMPC, 0);
	cbf->setUpOSQP(0);
    cbf->solveQP(10);

	cout << "[LOW LEVEL] cbf->uCBF[0]: " << cbf->uCBF[0] << endl;
	cout << "[LOW LEVEL] cbf->uCBF[1]: " << cbf->uCBF[1] << endl;

	Xn[5] = offset_angle_;
	X[5]  = offset_angle_ - 0.001;
	X[0] = 0.001;
	X[1] = 0.0005;
	X[4] = 0.0005;
	
	cbf->setIC(X, Xn);

	cbf->setMatrices(Alinear, Blinear, Clinear);
	cbf->evaluateCBFqpConstraintMatrices(uMPC, 0);
	cbf->setUpOSQP(0);
    cbf->solveQP(10);

	cout << "[LOW LEVEL] cbf->uCBF[0]: " << cbf->uCBF[0] << endl;
	cout << "[LOW LEVEL] cbf->uCBF[1]: " << cbf->uCBF[1] << endl;


	Xn[0] = offset_angle_;
	X[0]  = offset_angle_ + 0.01;
	X[1] = 0.001;
	X[2] = 0.0021;
	Alinear[2] = 0.001;
	Blinear[3] = 0.004;
	Clinear[4] = 0.003;

	X[0] = 0.001;
	X[1] = 0.0005;
	X[4] = 0.0005;

	uMPC[0] = 0.1;

	cbf->setIC(X, Xn);

	cbf->setMatrices(Alinear, Blinear, Clinear);
	cbf->evaluateCBFqpConstraintMatrices(uMPC, 0);
	cbf->setUpOSQP(0);
    cbf->solveQP(10);

	cout << "[LOW LEVEL] cbf->uCBF[0]: " << cbf->uCBF[0] << endl;
	cout << "[LOW LEVEL] cbf->uCBF[1]: " << cbf->uCBF[1] << endl;


}