#ifndef SEGWAY_EXAMPLE_H
#define SEGWAY_EXAMPLE_H

#include <cstring>
#include <string>

const int nu_ = 2;
const int nx_ = 4;
double x_eq_[nx_] = {0,0,0,0};

// Initialize Parameters

std::string prefix_path = "/home/drew/mpc";

void linearizedDroneDynamics(const double* x, 
				   const double* x_eq,
				   const double* u,
	               double* A,
	               double* B,
				   double* C,
				   const int idxN)
{
	int i = 0;
}


#endif