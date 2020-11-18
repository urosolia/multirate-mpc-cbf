#ifndef SEGWAY_EXAMPLE_H
#define SEGWAY_EXAMPLE_H

#include <cstring>
#include <string>

const int nu_ = 2;
const int nx_ = 7;

// Initialize Parameters
double mb = 44.798;
double mw = 2.485;
double Jw = 0.055936595310797;
double a2 = -0.02322718759275;
double c2 = 0.166845864363019;
double A2 = 3.604960049044268;
double B2 = 3.836289730154863;
double C2 = 1.069672194414735;
double K  = 1.261650363363571;
double r  = 0.195;
double L  = 0.5;
double gGravity = 9.81;
double FricCoeffViscous = 0.;
double velEps = 1.0e-3;
double FricCoeff = 1.225479467549329;

double x_max_[nx_] = {0.04,0.04,0.1,0.3,0.3,0.1,0.3};
double x_eq_[nx_] = {0,0,0,0,0,.138324423615,0};
double x_cost_[nx_] = {100,100,100,100,10000,10000,100};

std::string prefix_path = "/home/drew/mpc";

void linearizedDynamics(const double* x, 
				   const double* x_eq,
				   const double* u,
	               double* A,
	               double* B,
				   double* C,
				   const int idxN)
{
			double theta_t    = x[idxN*nx_ + 2];
			double v_t        = x[idxN*nx_ + 3];
			double thetaDot_t = x[idxN*nx_ + 4];
			double psi_t      = x[idxN*nx_ + 5] + x_eq_[5];
			double psiDot_t   = x[idxN*nx_ + 6];

			double u1_t       = u[idxN*nu_ + 0];
			double u2_t       = u[idxN*nu_ + 1];

			double fx = {};
			double gx[2] = {};
			double dfdx_5[nx_] = {};
			double dgdx_5_1[nx_] = {};
			double dgdx_5_2[nx_] = {};
			double dgdu_5[nu_];

			// Compute linearization x = v_t*cos(theta_t)
			A[0*nx_ + 2] = -v_t*sin(theta_t); // dpx/dtheta
			A[0*nx_ + 3] = cos(theta_t);     // dpx/dvt

			// Compute linearization y = v_t*sin(theta_t)
			A[1*nx_ + 2] =  v_t*cos(theta_t); // dpx/dtheta
			A[1*nx_ + 3] =  sin(theta_t);     // dpx/dvt

			// Compute affine term
			for (int row = 0; row < nx_; row++){
				C[idxN*nx_ + row] = 0.0;
			}
			C[idxN*nx_ + 0] = (v_t*cos(theta_t) + cos(theta_t)*(-v_t) - v_t*sin(theta_t)*(-theta_t) );
			C[idxN*nx_ + 1] = (v_t*sin(theta_t) + sin(theta_t)*(-v_t) + v_t*cos(theta_t)*(-theta_t) );
						
			fx =  pow(r,2)*thetaDot_t*((-2)*a2*mb*v_t*cos(psi_t)+(-4)*a2*c2*
			      mb*psiDot_t*cos(2*psi_t)+(-2)*(c2*mb*v_t+(A2+(-1)*C2+(-2)*
			      pow(a2,2)*mb+2*pow(c2,2)*mb)*psiDot_t*cos(psi_t))*sin(psi_t)) / (Jw*pow(L,2)+
			      pow(L,2)*mw*pow(r,2)+2*(C2+pow(a2,2)*mb)*pow(r,2)*pow(cos(psi_t),2)+2*(A2+pow(c2,2)*
			      mb)*pow(r,2)*pow(sin(psi_t),2)+2*a2*c2*mb*pow(r,2)*sin(2*psi_t));

			gx[0] =  (-1)*K*L/(Jw*pow(L,2)/r+pow(L,2)*mw*r+2*(C2+pow(a2,2)*mb)*r*pow(cos(psi_t),2)
					+2*(A2+pow(c2,2)*mb)*r*pow(sin(psi_t),2)+2*a2*c2*mb*r*sin(2*psi_t)),

			gx[1] = K*L/ ( Jw*pow(L,2)/r + pow(L,2)*mw*r + 2*(C2+pow(a2,2)*mb)*r*pow(cos(psi_t),2) 
					+ 2*(A2+pow(c2,2)*mb)*r*pow(sin(psi_t),2) + 2*a2*c2*mb*r*sin(2*psi_t));

    		dfdx_5[0]   = 0;
    		dfdx_5[1]   = 0;
    		dfdx_5[2]   = 0;
    		dfdx_5[3]   = (0.038025*thetaDot_t*(2.0811*cos(psi_t) - 14.949*sin(psi_t)))/(0.083187*pow(cos(psi_t),2) - 0.013203*sin(2.0*psi_t) + 0.369*pow(sin(psi_t), 2) + 0.037607);
    		dfdx_5[4]   = (0.038025*(0.69443*psiDot_t*cos(2.0*psi_t) + 2.0811*v_t*cos(psi_t) - 1.0*sin(psi_t)*(14.949*v_t + 9.9622*psiDot_t*cos(psi_t))))/(0.083187*pow(cos(psi_t),2) - 0.013203*sin(2.0*psi_t) + 0.369*pow(sin(psi_t), 2) + 0.037607);
    		dfdx_5[5]   = (0.038025*thetaDot_t*(0.026406*cos(2.0*psi_t) - 0.57162*cos(psi_t)*sin(psi_t))*(0.69443*psiDot_t*cos(2.0*psi_t) + 2.0811*v_t*cos(psi_t) - 1.0*sin(psi_t)*(14.949*v_t + 9.9622*psiDot_t*cos(psi_t))))/pow( (0.083187*pow(cos(psi_t), 2) - 0.013203*sin(2.0*psi_t) + 0.369*pow(sin(psi_t),2) + 0.037607), 2) - (0.038025*thetaDot_t*(1.3889*psiDot_t*sin(2.0*psi_t) - 9.9622*psiDot_t*pow(sin(psi_t),2) + 2.0811*v_t*sin(psi_t) + cos(psi_t)*(14.949*v_t + 9.9622*psiDot_t*cos(psi_t))))/(0.083187*pow(cos(psi_t),2) - 0.013203*sin(2.0*psi_t) + 0.369*pow(sin(psi_t),2) + 0.037607);
    		dfdx_5[6]   = (0.038025*thetaDot_t*(0.69443*cos(2.0*psi_t) - 9.9622*cos(psi_t)*sin(psi_t)))/(0.083187*pow(cos(psi_t),2) - 0.013203*sin(2.0*psi_t) + 0.369*pow(sin(psi_t),2) + 0.037607);

			dgdx_5_1[0] = 0;
			dgdx_5_1[1] = 0;
			dgdx_5_1[2] = 0;
			dgdx_5_1[3] = 0;
			dgdx_5_1[4] = 0;
			dgdx_5_1[5] = -(0.63083*(0.13541*cos(2.0*psi_t) - 2.9314*cos(psi_t)*sin(psi_t)))/ pow((0.4266*pow(cos(psi_t), 2) - 0.067707*sin(2.0*psi_t) + 1.8923*pow(sin(psi_t),2) + 0.19286), 2);
			dgdx_5_1[6] = 0;

			dgdx_5_2[0] = 0;
			dgdx_5_2[1] = 0;
			dgdx_5_2[2] = 0;
			dgdx_5_2[3] = 0;
			dgdx_5_2[4] = 0;
			dgdx_5_2[5] =  (0.63083*(0.13541*cos(2.0*psi_t) - 2.9314*cos(psi_t)*sin(psi_t)))/ pow((0.4266*pow(cos(psi_t), 2) - 0.067707*sin(2.0*psi_t) + 1.8923*pow(sin(psi_t),2) + 0.19286), 2);
			dgdx_5_2[6] = 0;
 
    		dgdu_5[0] = -0.63083/(0.4266*pow(cos(psi_t), 2) - 0.067707*sin(2.0*psi_t) + 1.8923*pow(sin(psi_t), 2) + 0.19286);
 		    dgdu_5[1] =  0.63083/(0.4266*pow(cos(psi_t), 2) - 0.067707*sin(2.0*psi_t) + 1.8923*pow(sin(psi_t), 2) + 0.19286);

 		    // Now update Alinear and Blinear
 		    for (int i = 0; i < nx_; i++){
 		    	// if (idxN == 0) cout << "dfdx_5: " << dfdx_5[i] << ",dgdx_5_1: " << dgdx_5_1[i]*u1_t + dgdx_5_2[i]*u2_t << endl;
 		    	A[4*nx_ + i] = dfdx_5[i] + dgdx_5_1[i]*u1_t + dgdx_5_2[i]*u2_t;
 		    	// A[4*nx_ + i] = dfdx_5[i] + dgdx_5_1[i]*u1_t + dgdx_5_1[i]*u2_t;
 		    }

 		    for (int i = 0; i < nu_; ++i){
 		    	B[4*nu_ + i] = dgdu_5[i];
 		    }

		    C[idxN*nx_ + 4] = fx + gx[0]*u1_t + gx[1]*u2_t;
			for (int i = 0; i < nx_; i++){
 		    	C[idxN*nx_ + 4] = C[idxN*nx_ + 4] + A[4*nx_ + i] * (-x[idxN*nx_ + i]);
 		    }
			for (int i = 0; i < nu_; i++){
 		    	C[idxN*nx_ + 4] = C[idxN*nx_ + 4] + B[4*nu_ + i] * (-u[idxN*nu_ + i]);
 		    }
}

void fullDynamics(const double *X, double *f, double *g)
{
	// Extract States
    double theta    = X[2]; 
    double v        = X[3]; 
    double thetaDot = X[4]; 
    double psi      = X[5]; 
    double psiDot   = X[6]; 

    double g1[nx_];
    double g2[nx_];

    double Fric = FricCoeff*tanh((v-psiDot*r)/velEps) + FricCoeffViscous*(v-psiDot*r);

    f[0] = v*cos(theta); 
    f[1] = v*sin(theta);
    f[2] = thetaDot;
    f[3] = (1/2)*r*(1/(4*B2*Jw+4*pow(a2, 2)*Jw*mb+4*pow(c2, 2)*Jw*mb+2*B2*mb*pow(r, 2)+pow(a2, 2)*pow(mb, 2)*pow(r, 2)+pow(c2, 2)*pow(mb, 2)*pow(r, 2)+4*B2*mw*pow(r, 2)+4*pow(a2, 2)*mb*mw*pow(r, 2)
    		+4*pow(c2, 2)*mb*mw*pow(r, 2)+(pow(a2, 2)+(-1)*pow(c2, 2))*pow(mb, 2)*pow(r, 2)*cos(2*psi)+2*a2*c2*pow(mb, 2)*pow(r, 2)*sin(2*psi)))*((-8)*B2*Fric+(-8)*pow(a2, 2)*Fric*mb
	  		+(-8)*pow(c2, 2)*Fric*mb+mb*r*((-8)*c2*Fric+a2*((-1)*A2+C2)*pow(thetaDot, 2)+4*a2*B2*(pow(psiDot, 2)+pow(thetaDot, 2))
	  		+pow(a2, 3)*mb*(4*pow(psiDot, 2)+3*pow(thetaDot, 2))+a2*pow(c2, 2)*mb*(4*pow(psiDot, 2)+3*pow(thetaDot, 2)))*cos(psi)+(-4)*a2*c2*gGravity*pow(mb, 2)*r*cos(2*psi)+a2*A2*mb*r*pow(thetaDot, 2)*cos(3*psi)
	  		+(-1)*a2*C2*mb*r*pow(thetaDot, 2)*cos(3*psi)+pow(a2, 3)*pow(mb, 2)*r*pow(thetaDot, 2)*cos(3*psi)+(-3)*a2*pow(c2, 2)*pow(mb, 2)*r*pow(thetaDot, 2)*cos(3*psi)+8*a2*Fric*mb*r*sin(psi)+4*B2*c2*mb*pow(psiDot, 2)*r*sin(psi)
	  		+4*pow(a2, 2)*c2*pow(mb, 2)*pow(psiDot, 2)*r*sin(psi)+4*pow(c2, 3)*pow(mb, 2)*pow(psiDot, 2)*r*sin(psi)+A2*c2*mb*r*pow(thetaDot, 2)*sin(psi)+4*B2*c2*mb*r*pow(thetaDot, 2)*sin(psi)+(-1)*c2*C2*mb*r*pow(thetaDot, 2)*sin(psi)
	  		+3*pow(a2, 2)*c2*pow(mb, 2)*r*pow(thetaDot, 2)*sin(psi)+3*pow(c2, 3)*pow(mb, 2)*r*pow(thetaDot, 2)*sin(psi)+ 2*pow(a2, 2)*gGravity*pow(mb, 2)*r*sin(2*psi)+(-2)*pow(c2, 2)*gGravity*pow(mb, 2)*r*sin(2*psi)+A2*c2*mb*r*pow(thetaDot, 2)*sin(3*psi)
	  		+(-1)*c2*C2*mb*r*pow(thetaDot, 2)*sin(3*psi)+3*pow(a2, 2)*c2*pow(mb, 2)*r*pow(thetaDot, 2)*sin(3*psi)+(-1)*pow(c2, 3)*pow(mb, 2)*r*pow(thetaDot, 2)*sin(3*psi));
    f[4] = pow(r, 2)*thetaDot*((-2)*a2*mb*v*cos(psi)+(-4)*a2*c2*mb*psiDot*cos(2*psi)+(-2)*(c2*mb*v+(A2+(-1)*C2+(-2)*pow(a2, 2)*mb+2*pow(c2, 2)*mb)*psiDot*cos(psi))*sin(psi))*(1/(Jw*pow(L, 2)+pow(L, 2)*mw*pow(r, 2)
    		+2*(C2+pow(a2, 2)*mb)*pow(r, 2)*pow(cos(psi), 2)+2*(A2+pow(c2, 2)*mb)*pow(r, 2)*pow(sin(psi), 2)+2*a2*c2*mb*pow(r, 2)*sin(2*psi)));
  	f[5] = psiDot;
	f[6] = (1/(4*B2*Jw+4*pow(a2, 2)*Jw*mb+4*pow(c2, 2)*Jw*mb+2*B2*mb*pow(r, 2)+pow(a2, 2)*pow(mb, 2)*pow(r, 2)+pow(c2, 2)*pow(mb, 2)*pow(r, 2)+4*B2*mw*pow(r, 2)+4*pow(a2, 2)*mb*mw*pow(r, 2)+4*pow(c2, 2)*mb*mw*pow(r, 2)+(pow(a2, 2)+(-1)*pow(c2, 2))*pow(mb, 2)*pow(r, 2)*cos(2*psi)
			+2*a2*c2*pow(mb, 2)*pow(r, 2)*sin(2*psi)))*(8*Fric*Jw+4*Fric*mb*pow(r, 2)+8*Fric*mw*pow(r, 2)+2*mb*(2*c2*Fric*r+a2*gGravity*(2*Jw+(mb+2*mw)*pow(r, 2)))*cos(psi)+(-2)*a2*c2*mb*(mb*pow(psiDot, 2)*pow(r, 2)
			+(-2)*(Jw+mw*pow(r, 2))*pow(thetaDot, 2))*cos(2*psi)+4*c2*gGravity*Jw*mb*sin(psi)+(-4)*a2*Fric*mb*r*sin(psi)+2*c2*gGravity*pow(mb, 2)*pow(r, 2)*sin(psi)+4*c2*gGravity*mb*mw*pow(r, 2)*sin(psi)
			+pow(a2, 2)*pow(mb, 2)*pow(psiDot, 2)*pow(r, 2)*sin(2*psi)+(-1)*pow(c2, 2)*pow(mb, 2)*pow(psiDot, 2)*pow(r, 2)*sin(2*psi)+(-2)*A2*Jw*pow(thetaDot, 2)*sin(2*psi)+2*C2*Jw*pow(thetaDot, 2)*sin(2*psi)+(-2)*pow(a2, 2)*Jw*mb*pow(thetaDot, 2)*sin(2*psi)
			+2*pow(c2, 2)*Jw*mb*pow(thetaDot, 2)*sin(2*psi)+(-1)*A2*mb*pow(r, 2)*pow(thetaDot, 2)*sin(2*psi)+C2*mb*pow(r, 2)*pow(thetaDot, 2)*sin(2*psi)+(-2)*A2*mw*pow(r, 2)*pow(thetaDot, 2)*sin(2*psi)+2*C2*mw*pow(r, 2)*pow(thetaDot, 2)*sin(2*psi)
			+(-2)*pow(a2, 2)*mb*mw*pow(r, 2)*pow(thetaDot, 2)*sin(2*psi)+2*pow(c2, 2)*mb*mw*pow(r, 2)*pow(thetaDot, 2)*sin(2*psi));

 	g1[0] = 0;
 	g2[0] = 0;
 	
 	g1[1] = 0;
 	g2[1] = 0;
 	
 	g1[2] = 0;
 	g2[2] = 0;
 	
 	g1[3] = K*r*(B2+pow(a2, 2)*mb+pow(c2, 2)*mb+c2*mb*r*cos(psi)+(-1)*a2*mb*r*sin(psi))*(1/(2*B2*Jw+2*pow(a2, 2)*Jw*mb+2*pow(c2, 2)*Jw*mb+B2*mb*pow(r, 2)
 			+pow(a2, 2)*pow(mb, 2)*pow(r, 2)+pow(c2, 2)*pow(mb, 2)*pow(r, 2)+2*B2*mw*pow(r, 2)+2*pow(a2, 2)*mb*mw*pow(r, 2)+2*pow(c2, 2)*mb*mw*pow(r, 2)
 			+(-1)*pow(c2, 2)*pow(mb, 2)*pow(r, 2)*pow(cos(psi), 2)+(-1)*pow(a2, 2)*pow(mb, 2)*pow(r, 2)*pow(sin(psi), 2)+a2*c2*pow(mb, 2)*pow(r, 2)*sin(2*psi)));
	g2[3] = K*r*(B2+pow(a2, 2)*mb+pow(c2, 2)*mb+c2*mb*r*cos(psi)+(-1)*a2*mb*r*sin(psi))*(1/(2*B2*Jw+2*pow(a2, 2)*Jw*mb+2*pow(c2, 2)*Jw*mb+B2*mb*pow(r, 2)
			+pow(a2, 2)*pow(mb, 2)*pow(r, 2)+pow(c2, 2)*pow(mb, 2)*pow(r, 2)+2*B2*mw*pow(r, 2)+2*pow(a2, 2)*mb*mw*pow(r, 2)+2*pow(c2, 2)*mb*mw*pow(r, 2)
			+(-1)*pow(c2, 2)*pow(mb, 2)*pow(r, 2)*pow(cos(psi), 2)+(-1)*pow(a2, 2)*pow(mb, 2)*pow(r, 2)*pow(sin(psi), 2)+a2*c2*pow(mb, 2)*pow(r, 2)*sin(2*psi)));

 	g1[4] = (-1)*K*L*(1/(Jw*pow(L, 2)*1/r+pow(L, 2)*mw*r+2*(C2+pow(a2, 2)*mb)*r*pow(cos(psi), 2)+2*(A2+pow(c2, 2)*mb)*r*pow(sin(psi), 2)+2*a2*c2*mb*r*sin(2*psi)));
	g2[4] = K*L*(1/(Jw*pow(L, 2)*1/r+pow(L, 2)*mw*r+2*(C2+pow(a2, 2)*mb)*r*pow(cos(psi), 2)+2*(A2+pow(c2, 2)*mb)*r*pow(sin(psi), 2)+2*a2*c2*mb*r*sin(2*psi)));

 	g1[5] = 0;
	g2[5] = 0;

 	g1[6] = (-2)*K*(2*Jw+mb*pow(r, 2)+2*mw*pow(r, 2)+c2*mb*r*cos(psi)+(-1)*a2*mb*r*sin(psi))*(1/(4*B2*Jw+4*pow(a2, 2)*Jw*mb+4*pow(c2, 2)*Jw*mb+2*B2*mb*pow(r, 2)
 			+pow(a2, 2)*pow(mb, 2)*pow(r, 2)+pow(c2, 2)*pow(mb, 2)*pow(r, 2)+4*B2*mw*pow(r, 2)+4*pow(a2, 2)*mb*mw*pow(r, 2)
 			+4*pow(c2, 2)*mb*mw*pow(r, 2)+(pow(a2, 2)+(-1)*pow(c2, 2))*pow(mb, 2)*pow(r, 2)*cos(2*psi)+2*a2*c2*pow(mb, 2)*pow(r, 2)*sin(2*psi)));
	g2[6] = (-2)*K*(2*Jw+mb*pow(r, 2)+2*mw*pow(r, 2)+c2*mb*r*cos(psi)+(-1)*a2*mb*r*sin(psi))*(1/(4*B2*Jw+4*pow(a2, 2)*Jw*mb+4*pow(c2, 2)*Jw*mb+2*B2*mb*pow(r, 2)
			+pow(a2, 2)*pow(mb, 2)*pow(r, 2)+pow(c2, 2)*pow(mb, 2)*pow(r, 2)+4*B2*mw*pow(r, 2)+4*pow(a2, 2)*mb*mw*pow(r, 2)+4*pow(c2, 2)*mb*mw*pow(r, 2)
  			+(pow(a2, 2)+(-1)*pow(c2, 2))*pow(mb, 2)*pow(r, 2)*cos(2*psi)+2*a2*c2*pow(mb, 2)*pow(r, 2)*sin(2*psi)));

	memcpy(g, g1, sizeof(double)*nx_);
	memcpy(g+nx_, g2, sizeof(double)*nx_);
}

void safetySet(const double X[], const double Xn[], double *h, double *Dh) 
{
	double x        = X[0]; 
    double y        = X[1]; 
    double theta    = X[2]; 
    double v        = X[3]; 
    double thetaDot = X[4]; 
    double psi      = X[5]; 
    double psiDot   = X[6]; 

    double xn        = Xn[0]; 
    double yn        = Xn[1]; 
    double thetan    = Xn[2]; 
    double vn        = Xn[3]; 
    double thetaDotn = Xn[4]; 
    double psin      = Xn[5]; 
    double psiDotn   = Xn[6]; 
	h[0] = 1 - pow((x - xn),2)/pow(x_max_[0],2) - pow((y - yn),2)/pow(x_max_[1],2)  - pow((theta - thetan),2)/pow(x_max_[2],2)  
              - pow((v - vn),2)/pow(x_max_[3],2) - pow((thetaDot - thetaDotn),2)/pow(x_max_[4],2)
        	  - pow((psi - psin),2)/pow(x_max_[5],2)- pow((psiDot - psiDotn),2)/pow(x_max_[6],2);


	Dh[0] = -(2*x - 2*xn)/pow(x_max_[0],2);
    Dh[1] = -(2*y - 2*yn)/pow(x_max_[1],2);
	Dh[2] = -(2*theta - 2*thetan)/pow(x_max_[2],2);
    Dh[3] = -(2*v - 2*vn)/pow(x_max_[3],2);
    Dh[4] = -(2*thetaDot - 2*thetaDotn)/pow(x_max_[4],2);
    Dh[5] = -(2*psi - 2*psin)/pow(x_max_[5],2);
    Dh[6] = -(2*psiDot - 2*psiDotn)/pow(x_max_[6],2);
	Dh[7] = (2*x - 2*xn)/pow(x_max_[0],2);
    Dh[8] = (2*y - 2*yn)/pow(x_max_[1],2);
	Dh[9] = (2*theta - 2*thetan)/pow(x_max_[2],2);
    Dh[10] = (2*v - 2*vn)/pow(x_max_[3],2);
    Dh[11] = (2*thetaDot - 2*thetaDotn)/pow(x_max_[4],2);
    Dh[12] = (2*psi - 2*psin)/pow(x_max_[5],2);
    Dh[13] = (2*psiDot - 2*psiDotn)/pow(x_max_[6],2);
}
void clf(const double X[], const double Xn[], double *V, double *Dv)
{
	double x        = X[0]; 
	double y        = X[1]; 
	double theta    = X[2]; 
	double v        = X[3]; 
	double thetaDot = X[4]; 
	double psi      = X[5]; 
	double psiDot   = X[6]; 

	double xn        = Xn[0]; 
	double yn        = Xn[1]; 
	double thetan    = Xn[2]; 
	double vn        = Xn[3]; 
	double thetaDotn = Xn[4]; 
	double psin      = Xn[5]; 
	double psiDotn   = Xn[6]; 

    V[0] =  x_cost_[0]*pow((x - xn),2) + x_cost_[1]*pow((y - yn),2) + x_cost_[2]*pow((theta - thetan),2) 
		+ x_cost_[3]*pow((v - vn),2) + x_cost_[4]*pow((thetaDot - thetaDotn),2) 
		+ x_cost_[5]*pow((psi - psin),2) + x_cost_[6]*pow((psiDot - psiDotn),2) ;

	
         
    Dv[0] = x_cost_[0]*(2*x - 2*xn),
    Dv[1] = x_cost_[1]*(2*y - 2*yn),
   	Dv[2] = x_cost_[2]*(2*theta - 2*thetan),
    Dv[3] = x_cost_[3]*(2*v - 2*vn),
	Dv[4] = x_cost_[4]*(2*thetaDot - 2*thetaDotn),
    Dv[5] = x_cost_[5]*(2*psi - 2*psin),
    Dv[6] = x_cost_[6]*(2*psiDot - 2*psiDotn),
    Dv[7] = -x_cost_[0]*(2*x - 2*xn),
    Dv[8] = -x_cost_[1]*(2*y - 2*yn),
   	Dv[9] = -x_cost_[2]*(2*theta - 2*thetan),
    Dv[10] = -x_cost_[3]*(2*v - 2*vn),
	Dv[11] = -x_cost_[4]*(2*thetaDot - 2*thetaDotn),
    Dv[12] = -x_cost_[5]*(2*psi - 2*psin),
    Dv[13] = -x_cost_[6]*(2*psiDot - 2*psiDotn);
}

#endif