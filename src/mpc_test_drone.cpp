#include <iostream>
#include <fstream>
#include <string>
#include <osqp.h>
#include "mpc/mpc.hpp"
#include "mpc/drone_example.hpp"

// Global declarations
std::vector<double> inputBuffer_u1_;
std::vector<double> inputBuffer_u2_;

using namespace std;
using namespace ModelPredictiveControllerValFun;
double dt_;
double dtDelay_;
double delay_ms_;
double x_start;
double y_start;
double linearization_IC_[nx_];
double vMax = {0.0};

double error_max_[nx_]        = {};
double enlarge    = {};
bool lowLevelActive;

int flag_state_measurement = 0;
int flag_goalSetAndState_measurement = 0;

// MPC *mpc; // Define MPC object as a pointer
MPCValFun *mpcValFun; // Define MPC object as a pointer

//// Main Definition
int main (int argc, char *argv[])
{
	dt_ = 0.05;

	// Read Launch File
	delay_ms_ = 0;
	linearization_IC_[0] = .5;
	linearization_IC_[1] = 6.5;
	dtDelay_ = delay_ms_/1000.0;
	
	enlarge = 0;
	lowLevelActive = false;

	// Initialize QP parameters. 
	// static const uint32_t constraint = 0; // 0 = no printing, 1 = only x_t, 2 = minimal, 3 = minimal + initialization, 4 = all

	static const int N  = 15; 	  // Horizon length

	static const uint32_t printLevel =  0; // 0 = no printing, 1 = only x_t, 2 = minimal, 3 = minimal + initialization, 4 = all

	double goalSetAndStateVector[11] = {0.0};
	// States are:        x,   y, theta,   v, thetaDot,            psi, psiDot

	// Initial condition and goal in relative reference frame
	double x_IC[nx_] = {  0.5, 6.5,   0.0,    0.0};  // initial condition to
	double x_g[nx_]  = {  0.5, 6.5,   0.0,    0.0};  // initial condition to

	goalSetAndStateVector[0]  = x_g[0];
	goalSetAndStateVector[1]  = x_g[1];
	goalSetAndStateVector[2]  = -10;
	goalSetAndStateVector[3]  = 10;
	goalSetAndStateVector[4]  = -10;
	goalSetAndStateVector[5]  = 10;
	goalSetAndStateVector[6]  = 0;
	goalSetAndStateVector[7]  = -10;
	goalSetAndStateVector[8]  = 10;
	goalSetAndStateVector[9]  = -10;
	goalSetAndStateVector[10] = 10;


	cout << "========== START Initializing MPC Object" << endl;
	string matrix_prefix_path = "/home/drew/multirate-mpc-cbf";
	mpcValFun = new MPCValFun(nx_, nu_, N, dt_, dtDelay_, printLevel, x_eq_, error_max_, enlarge, lowLevelActive, linearization_IC_, matrix_prefix_path, linearizedDroneDynamics);
	mpcValFun->setGoalState(x_g);
	mpcValFun->setIC(x_IC);  // Solve QP to check the everything works
	mpcValFun->readAB();
	mpcValFun->readCost();
	
	mpcValFun->initiConstrVector();

	mpcValFun->linearize();

	mpcValFun->updateGoalSetAndState(goalSetAndStateVector);
	
	mpcValFun->buildConstrMatrix();
	mpcValFun->buildCost();
	mpcValFun->buildConstrVector();
	
	mpcValFun->setUpOSQP(1); // Initialize QP with verbose = 1

	// mpcValFun->setUpOSQP(0); // Set Verbose option to 0
	cout <<"========== DONE Initializing MPC Object" << endl;

	mpcValFun->solveQP();  // Solve QP to check the everything works

	mpcValFun->linearize();			
	mpcValFun->buildCost();
	mpcValFun->buildConstrMatrix();
	mpcValFun->buildConstrVector();
	mpcValFun->solveQP();  // Solve QP to check the everything works

	cout << "Solver flag: " << mpcValFun->solverFlag << endl;
	cout << "OPTIMAL States New:"<< std::endl;
	for (int i = 0; i< nx_; i++){
		for (int j = 0; j < N+1; ++j)
			{cout << mpcValFun->xPred[j*nx_+i] <<",";}
		cout << endl;
	}	
	cout << endl;

	cout << "OPTIMAL Inputs:"<< std::endl;
	for (int i = 0; i< nu_; i++){
		for (int j = 0; j < N; ++j)
			cout << mpcValFun->uPred[j*nu_+i] <<",";
	cout << endl;
	}	
	cout << endl;

	if (delay_ms_>0){
		inputBuffer_u1_.resize( int(dt_/(delay_ms_/1000)) );
		inputBuffer_u2_.resize( int(dt_/(delay_ms_/1000)) );
		cout << "===== MPC delay: " << delay_ms_ << ". Initialize inputBuffer_u1_: ";
		for (int i = 0; i < int(dt_/(delay_ms_/1000)); ++i){
			cout << inputBuffer_u1_[i] << ", ";
		}
		cout << endl;
	}

	// Take it for a spin

}