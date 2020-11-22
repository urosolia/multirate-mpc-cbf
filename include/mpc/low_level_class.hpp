#ifndef LLC_H
#define LLC_H

#include <math.h>       /* pow */
#include <fstream>
#include <cstring>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/iterator/n_step_iterator.hpp>
enum class STATUS : uint8_t
{
	FAILURE = 0,
	RUNNING = 1
};

namespace ControlBarrierFunction
{	
	///////////////////////////// CBF

	class CBF
	{

	public:
		CBF(const int nx, const int nu, const bool lowLevelActive, c_float Hx[], const double x_eq_[], const double x_max_[],
			std::function<void(const double* /*x*/,
					       double* /*A*/,
					       double* /*B*/)> fullDynamics,
			std::function<void(const double* /*x*/,
							   const double* /*x_eq*/,
						       const double* /*x_max*/,
						       double* /*A*/,
						       double* /*B*/)> safetySet,
			std::function<void(const double* /*x*/,
							   const double* /*x_eq*/,
						       double* /*A*/,
						       double* /*B*/)> clf);

		~CBF();

		void evaluateCBFqpConstraintMatrices(double uMPC[], int printLevel);
		void setIC(double X_in[], double Xn_in[]);
		void setMatrices(double A_in[], double B_in[], double C_in[]);
		int32_t setUpOSQP(int verbose);
		void solveQP(int printLevel);

		double *uCBF;
		double gamma;
		double h[1];
		double V[1];

		bool clfActive;

		double h_alpha[1] = {0.1};
		double V_alpha[1] = {10.0};

		double uDes[2] = {0.0, 0.0};


		double *X;
		double *Xn;

		double flagQP;
	protected:
		const int nx_;
		const int nu_;
		const bool lowLevelActive_;
		double *f;
		double *g;
		double *Alinear_;
		double *Blinear_;
		double *Clinear_;
		c_float *Hx_;//[nv_]  = {0.0};
		const double *x_eq_;
		const double *x_max_;
		std::function<void(const double* /*x*/,
	                             double* /*f*/,
								 double* /*g*/)> fullDynamics_;
		std::function<void(const double* /*x*/,
								 const double* /*x*/,
								 const double* /*x*/,
	                             double* /*f*/,
								 double* /*g*/)> safetySet_;
		std::function<void(const double* /*x*/,
						   const double* /*x*/,
	                             double* /*f*/,
								 double* /*g*/)> clf_;		

		// QP Constraint Matrices (the structure is fixed upper triangular)
		int nv_ = (int) nu_+1;
		c_float *Fx_;//[2*nu_+1] = {0.0};
		c_int *Fr_;//[2*nu_+1] = {0}; // CSC representation where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]]
		c_int *Fc_;//[nv_+1] = {0};
		c_int F_nnz_ = 2*nu_+1;

		// QP lower and upper bounds
		c_float lb_x_[2] = {0.0, 0.0};
		c_float ub_x_[2] = {0.0, 0.0};

		// QP Cost Matrices (There are constants)
		c_int   *Hr_;//[nv_]  = {0.0};
		c_int   *Hc_;//[nv_+1]  = {0.0};
		c_int H_nnz_ = nv_;
		c_float *qv_;//[nv_]  = {0.0}; //Linear cost vector
		
		// OSQP settings
		OSQPSettings *settings_;
		OSQPWorkspace *work_;
		OSQPData *data_;

	}; //end CBF class

	// Define constructor
	CBF::CBF(const int nx, const int nu, const bool lowLevelActive, c_float Hx[], const double x_eq[], const double x_max[],
		std::function<void(const double* /*x*/,
					       double* /*A*/,
					       double* /*B*/)> fullDynamics,
		std::function<void(const double* /*x*/,
						   const double* /*x_eq*/,
						   const double* /*x_max*/,
					       double* /*A*/,
					       double* /*B*/)> safetySet,
		std::function<void(const double* /*x*/,
						   const double* /*x_eq*/,
					       double* /*A*/,
					       double* /*B*/)> clf):
	nx_(nx), nu_(nu), lowLevelActive_(lowLevelActive), Hx_(Hx), x_eq_(x_eq), x_max_(x_max), fullDynamics_(fullDynamics), safetySet_(safetySet), clf_(clf)
	{
		Fx_ = new c_float[2*nu_+1]{};
		Fr_ = new c_int[2*nu_+1]{};
		Fc_ = new c_int[nv_+1]{};

		Hr_ = new c_int[nv_]{};
		Hc_ = new c_int[nv_+1]{};
		qv_ = new c_float[nv_]{};

		uCBF = new double[nu_]{};

		//Initialize Fr
		for (int i = 0; i < nu_;i++) {
			Fr_[2*i+1] = 1;
		}

		//Initialize Fc
		for (int i = 0; i < nv_; i++) {
			Fc_[i] = 2*i;
		}
		Fc_[nv_] = nu_*2+1;

		//Initialize Hr
		for (int i = 0; i < nv_; i++) {
			Hr_[i] = i;
		}

		//Initialize Hc
		for (int i = 0; i < (nv_ + 1); i++) {
			Hc_[i] = i;
		}

 		// Initialize OSQP
		settings_ = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
		data_ = (OSQPData *)c_malloc(sizeof(OSQPData));


		// Initialize Matrices to planning system
		f  = new double[nx_]{};
		g = new double[nx_*nu_]{};

		X = new double[nx_]{};
		Xn = new double[nx_]{};

		Alinear_  = new double[nx_*nx_]{}; 
		Blinear_  = new double[nx_*nu_]{}; 
		Clinear_  = new double[nx_]{};

		std::cout << "Low level x_max_" << std::endl;
		for (int i=0; i < nx; i++){
			std::cout << "i: "<< i << ", x_max_[i]: " << x_max_[i] << std::endl;
		}
 	}

 	// Deallocation Memory
	CBF::~CBF(void)
	{
		osqp_cleanup(work_);
		if (data_)
		{
			c_free(data_->A);
			c_free(data_->P);
			c_free(data_);
		}
		if (settings_)  c_free(settings_);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
 	/////////////////////////////////////// Define Functions /////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void CBF::setIC(double X_in[], double Xn_in[])
	{
		for (int i = 0; i < nx_; i++ ){
			X[i] = X_in[i];
			Xn[i] = Xn_in[i];
		}
	}

	void CBF::setMatrices(double A_in[], double B_in[], double C_in[])
	{
		for (int i = 0; i < nx_*nx_; i++ ){
			Alinear_[i] = A_in[i];
		}

		for (int i = 0; i < nx_*nu_; i++ ){
			Blinear_[i] = B_in[i];
		}

		for (int i = 0; i < nx_; i++ ){
			Clinear_[i] = C_in[i];
		}
	}
	int32_t CBF::setUpOSQP(int verbose)
	{
		std::cout << "START Setting-up Solver" << std::endl;
		// Setting default settings
		osqp_set_default_settings(settings_);
		settings_->max_iter = 2000;
		// settings_->polish  = 1;
		// settings_->verbose = verbose;

		// Fill up osqp data
		data_->n = 3; // two varaibles
		data_->m = 2; // two constraints
		data_->P = csc_matrix(data_->n, data_->n, H_nnz_, Hx_, Hr_, Hc_);
		data_->q = qv_;
		data_->A = csc_matrix(data_->m, data_->n, F_nnz_, Fx_, Fr_, Fc_);
		data_->l = lb_x_;
		data_->u = ub_x_;

		std::cout << "DONE Setting-up Solver" << std::endl;
		
		return osqp_setup(&work_, data_, settings_);
	}

	void CBF::solveQP(int printLevel)
	{
		// Update constraints
		data_->l = lb_x_;
		data_->u = ub_x_;
		data_->A = csc_matrix(data_->m, data_->n, F_nnz_, Fx_, Fr_, Fc_);

		// osqp_setup(&work_, data_, settings_);
		osqp_update_lower_bound(work_,lb_x_);
		osqp_update_upper_bound(work_,ub_x_);
		for (int i = 0; i < nu_; i++){
			qv_[i] = -2 * Hx_[i] * uDes[i];
		}	
		osqp_update_lin_cost(work_,qv_);

		osqp_update_A(work_, Fx_, OSQP_NULL, 2*nu_+1);

		if (printLevel >=1 ) std::cout << "Solving QP... "<< std::endl;
		osqp_solve(work_);
		if ((work_->info->status_val != OSQP_SOLVED) and (lowLevelActive_ > 0.5)){
			std::cout << "CBF QP not optimal, solver flag: " << work_->info->status_val << std::endl;
		}
		flagQP = work_->info->status_val;

		uCBF[0]  = static_cast<double>(work_->solution->x[0]);
		uCBF[1]  = static_cast<double>(work_->solution->x[1]);
		gamma = static_cast<double>(work_->solution->x[2]);

		if (printLevel >=1 ) {
			std::cout << " OPTIMAL SOLUTION " << std::endl;
			std::cout << "uCBF: " << uCBF[0]<< ", " << uCBF[1] << ", gamma: "<< gamma << std::endl;			
		}
	}

	void CBF::evaluateCBFqpConstraintMatrices(double uMPC[], int printLevel)
	{
		// evaluate dynamics
		fullDynamics_(X,f,g);

        double dxdt_realsys_constaTerm[nx_];
        double dxdt_realsys_linearTerm[nx_*nu_];

        for (int i=0; i<=nx_; ++i){
        	dxdt_realsys_constaTerm[i] = f[i];
        	for (int j = 0; j < nu_; j++){
        		dxdt_realsys_constaTerm[i] += g[i + j*nx_]*uMPC[j];
        		dxdt_realsys_linearTerm[i + j*nx_] = g[i + j*nx_];
        	}
        	// std::cout << "g1[i]: " << g1[i] << ", i: " << i << std::endl;
        	// std::cout << "g2[i]: " << g2[i] << ", i: " << i << std::endl;
        }

        double dxdt_plansys_constaTerm[nx_];
        double dxdt_plansys_linearTerm[nx_*nu_] = {0.0};
 
        for (int i=0; i < nx_; ++i){
        	dxdt_plansys_constaTerm[i] = 0.0;
        	// std::cout << "Alinear_: ";
        	for (int j=0; j < nx_; ++j)
        	{
        		dxdt_plansys_constaTerm[i] = dxdt_plansys_constaTerm[i] + Alinear_[i*nx_ + j] *(Xn[j] - x_eq_[j]);
        		// std::cout << Alinear_[i*nx_ + j] << ",";
        	}
        	// std::cout << " i: " << i << std::endl;
        	dxdt_plansys_constaTerm[i] = dxdt_plansys_constaTerm[i] + Blinear_[i*nu_ + 0] * uMPC[0] + Blinear_[i*nu_ + 1] * uMPC[1] + Clinear_[i];
        	// std::cout << "dxdt_plansys_constaTerm: " << dxdt_plansys_constaTerm[i] << " i: " << i << std::endl;
        	// std::cout << "Blinear_: " << Blinear_[i*nx_ + 0] << ", " << Blinear_[i*nx_ + 0] << " i: " << i << std::endl;
        	// std::cout << "Clinear_: " << Clinear_[i]  << " i: " << i << std::endl;
        	// std::cout << "uMPC: " << uMPC[0] << ", " << uMPC[1] << " i: " << i << std::endl;
        	
        }

        // Aug Sys dynamcis 
        double dxdt_augmsys_constaTerm[2*nx_];
        double dxdt_augmsys_linearTerm[2*nx_*nu_];
        memcpy(dxdt_augmsys_constaTerm,dxdt_realsys_constaTerm,sizeof(double)*nx_);
        memcpy(dxdt_augmsys_constaTerm+nx_,dxdt_plansys_constaTerm,sizeof(double)*nx_);
		memcpy(dxdt_augmsys_linearTerm,dxdt_realsys_linearTerm,sizeof(double)*nx_*nu_);
        memcpy(dxdt_augmsys_linearTerm+(nx_*nu_),dxdt_plansys_linearTerm,sizeof(double)*nx_*nu_);
        
        // Barrier
        double Dh[2*nx_];
        double Dv[2*nx_];
        safetySet_(X, Xn, x_max_, h, Dh);
        clf_(X, Xn, V, Dv);

        // constraints
		// write upper bounds [-10.0*self.V -np.dot(dvdx, dxdt_augmsys_constaTerm) ,  np.inf];
        // write lower bounds [-np.inf                                             , -10.0*self.h -np.dot(dhdx, dxdt_augmsys_constaTerm)];
		double temp1 = 0.0;
		double temp2 = 0.0;
		for (int i=0; i<2*nx_; ++i){
			// std::cout << "dxdt_augmsys_constaTerm " << dxdt_augmsys_constaTerm[i] << std::endl;
			temp1 = temp1 + Dv[i] * dxdt_augmsys_constaTerm[i];
			temp2 = temp2 + Dh[i] * dxdt_augmsys_constaTerm[i];
		}
		if (V_alpha[0]<10000.0){
			ub_x_[0] = - V_alpha[0] *V[0] - temp1; 			
		}else{
			ub_x_[0] = OSQP_INFTY; 			
		}
		ub_x_[1] = OSQP_INFTY;
		lb_x_[0] = -OSQP_INFTY;
		lb_x_[1] = - h_alpha[0] *h[0] - temp2; 

		// write constraint matrix
        // self.A = csc_matrix(np.array([[np.dot(dvdx, dxdt_augmsys_linearTerm1), np.dot(dvdx, dxdt_augmsys_linearTerm2), -1],
        //                               [np.dot(dhdx, dxdt_augmsys_linearTerm1), np.dot(dhdx, dxdt_augmsys_linearTerm1),  0]]));
		for (int i = 0; i < 2*nu_+1; i++)
			Fx_[i] = 0;

		for (int j=0; j<nu_; ++j){
			for (int  i=0; i<nx_; i++) {
				Fx_[2*j] += Dv[i] * dxdt_augmsys_linearTerm[i+j*nx_];
				Fx_[2*j+1] += Dh[i] * dxdt_augmsys_linearTerm[i+j*nx_];
			}
		}
		Fx_[2*nu_] = -1;

		if (printLevel >= 1){
			std::cout << "ub: " << ub_x_[0] << ", " << ub_x_[1] << std::endl;
			std::cout << "lb: " << lb_x_[0] << ", " << lb_x_[1] << std::endl;
			for (int i = 0; i < 2*nu_+1; i++)
				std::cout << "Fx_[" << i << "]: " << Fx_[i] << std::endl;
		}

	}

	
}

#endif