#ifndef MPC_H
#define MPC_H

#include <osqp.h>
#include <fstream>
#include <functional>
#include <cstring>
#include <string>
#include <vector>

enum class STATUS : uint8_t
{
	FAILURE = 0,
	RUNNING = 1
};

using namespace std;


namespace ModelPredictiveControllerValFun
{	
	///////////////////////////// MPC
	class MPCValFun
	{

	public:
		MPCValFun(const int nx,
			const int nu,
			const int N,
			const double dt,
			const double dtDelay,
			const uint32_t printFlag,
			const double x_eq[],
			const double error_max[],
			const double enlarge,
			const bool lowLevelActive,
			const double linearization_IC[],
			const string matrix_prefix_path,
			std::function<void(const double* /*x*/,
							   const double* /*x_eq*/,
							   const double* /*u*/,
		                             double* /*A*/,
		                             double* /*B*/,
									 double* /*C*/,
									 const int /*idxN*/)> linearizedDynamics);

		~MPCValFun(void);

		void matrixProd(double Aout[], double A1[], double A2[]);
		void setIC(double xt[]);
		void oneStepPrediction(double ut[]);
		void oneStepPredictionDelay(double ut[]);
		void solveQP();
		void readAB(void);
		int setUpOSQP(int verbose);
		void solveQP(double xt[]);
		void initiConstrVector(void);
		void buildConstrVector(void);
		void buildConstrMatrix(void);
		void readCost(void);
		void buildCost(void);
		void linearize(void);
		void updateHorizon(void);
		void resetHorizon(void);
		void setGoalState(double xt[]);
		void updateGoalSetAndState(double goalSetAndState[]);

		double *xPred;
		double *uPred;
		double solverFlag = 0.0;
		double *AlinearOut;
		double *BlinearOut;
		double *ClinearOut;
		double *x_IC_;

	protected:
		double *xLin_;
		double *uLin_;

		const int nx_;
		const int nu_;
		const int  N_;
		const double dt_;
		const double dtDelay_;
		const uint32_t  printLevel_;
		const double *x_eq_;
	    const double *max_error_;
		const double enlarge_ = {};
	    const bool lowLevelActive_ = {};
		const double *linearization_IC_;
		const string matrix_prefix_path_;
		std::function<void(const double* /*x*/,
						   const double* /*x_eq*/,
						   const double* /*u*/,
	                             double* /*A*/,
	                             double* /*B*/,
								 double* /*C*/,
						   const int /*idxN*/)> linearizedDynamics_;
		int nv_;
		int nc_;
		double *x_g_;
		double *Alinear2_;
		double *Alinear3_;
		double *Alinear4_;


		double *Alin_;
		double *Blin_;

		double *Alinear_;
		double *Blinear_;
		double *Clinear_;
		double *AlinDelay_;
		double *BlinDelay_;
		double *xConstrLB_;
		double *xConstrUB_;
		double *xConstrTermLB_;
		double *xConstrTermUB_;
		double *uConstr_;

		c_float *Hx_;
		c_int *Hr_;
		c_int *Hc_;
		c_int H_nnz_;
		c_float *qv_; //Linear cost vector

		c_float *Fx_;
		c_int *Fr_;
		c_int *Fc_;
		c_int F_znn_;
		int Nterm_;


		c_float *lb_x_;
		c_float *ub_x_;

		double highLevTime = 0.0;

		// States are:                  x,    y, theta,    v, thetaDot,  psi,  psiDot
		double *tightening;//{0.02, 0.02,  0.01, 0.01,     0.01, 0.01,    0.01};
		double *oldInput;
		
		float *xCost_;
		float *xCostFinal_;
		float *uCost_;
		float *Qlin_;
		float *Qstate_rate;

		float *Qerr ;
		float Qrate_ = 0.0;
		float Qdiff_ = 0.0;

		// OSQP settings
		OSQPSettings *settings_;
		OSQPWorkspace *work_;
		OSQPData *data_;

	}; //end MPC class

	// Define constructor
	MPCValFun::MPCValFun(const int nx,  
						 const int nu, 
						 const int N, 
						 const double dt, 
						 const double dtDelay, 
						 const uint32_t printFlag, 
						 const double x_eq[], 
						 const double max_error[], 
						 const double enlarge, 
						 const bool lowLevelActive, 
						 const double linearization_IC[],
						 const string matrix_prefix_path, 
						 std::function<void(const double* /*x*/,
										    const double* /*x_eq*/,
										    const double* /*u*/,
					                              double* /*A*/,
					                              double* /*B*/,
											  	  double* /*C*/,
									 			  const int /*idxN*/)> linearizedDynamics):
	nx_(nx),nu_(nu), N_(N), dt_(dt), dtDelay_(dtDelay), printLevel_(printFlag), x_eq_(x_eq), max_error_(max_error), enlarge_(enlarge), lowLevelActive_(lowLevelActive), linearization_IC_(linearization_IC), matrix_prefix_path_(matrix_prefix_path), linearizedDynamics_(linearizedDynamics)
	{
					
		// Initialize horizon
		Nterm_ = N_;

		// Tightening paramters
		cout << "[MPC] Low Level Active: " << lowLevelActive_ << endl;
		cout << "[MPC] Enlarge: " << enlarge << endl;

		cout << "[MPC] linearization_IC_: " << linearization_IC_[2] << endl;

		tightening = new double[nx_]{};
		Qerr = new float[nx_]{};
		oldInput = new double[nu_]{};

		if (lowLevelActive_ == true){
			cout << " Low Level Active --> Use Robust MPC with tightening:" << endl;
			for (int i = 0; i < nx_; i++){
				tightening[i] = max_error_[i];	
				cout << "i :" << i << " : " << tightening[i] << endl; 
			}
			ifstream myfile;
			myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/Qerr.txt"), std::ios::app);
			for (int i = 0; i < nx_; i++){
				myfile >> Qerr[i];
			}
			myfile.close();
		}else{
			for (int i = 0; i < nx_; i++){
				tightening[i] = 0;
				Qerr[i]       = 0;	
				cout << "i :" << i << " : " << tightening[i] << endl; 
			}
		}
		cout << "Tightening Paramters: " << endl;
		cout << "error x: "        << tightening[0] << endl;
		cout << "error y: "        << tightening[1] << endl;
		cout << "error theta: "    << tightening[2] << endl;
		cout << "error v: "        << tightening[3] << endl;
		for (int i = 0; i < nx_; i++){
			cout << "error cost: "     << Qerr[i] << " i: " << i << endl;
		}

		// Initialize Constraint
		nv_ = nx_*(N_+1) + nu_*N_; // Total number of variable
		nc_ = nx_*(N_+1) + nx_*(N_+1) + nu_*N_; // Total number of constraints: nx_*(N+1) equality constraints + (nx_*(N_+1) + nu_*N_) box constraints
		
		F_znn_ = nx*(N_+1) + nx*nx*N_ + nx*nu*N_ + nx_*(N_+1) + nu_*N_;
		Fx_    = new c_float[F_znn_];
		Fr_    = new c_int[F_znn_]; // Note that dim(Hr_) = dim(FxOld_) by definition
		Fc_    = new c_int[nv_ + 1];
		lb_x_  = new c_float[nc_];
		ub_x_  = new c_float[nc_];

		// Initialize Hessian
		// nx_*N_         = diagonal state cost
		// nu_*N_         = diagonal input cost
		// nu_*(N_ - 1)*2 = off diagonal rate cost (u1[k] - u1[k+1])^2 + (u2[k] - u2[k+1])^2
		// nx_*(N_)*2     = off diagonal rate state cost (x1[k] - x1[k+1])^2 etc .. 
		// nu_*N          = off diagonal difference cost (u1[k] - u2[k])^2
		H_nnz_ = nx_*(N_+1) + nu_*N_ + nu_*(N_ - 1) + nx_*N_ + N_ ;
		Hx_    = new c_float[H_nnz_];
		Hr_    = new c_int[H_nnz_]; // Note that dim(Hr_) = dim(Hx_) by definition
		Hc_    = new c_int[nv_ + 1]; // Note dim Hc_ is number of non-zero columns + 1
		qv_    = new c_float[nv_];

		x_IC_ = new double[nx_];
		x_g_  = new double[nx_];

		xConstrLB_     = new double[nx_];
		xConstrUB_     = new double[nx_];
		xConstrTermLB_ = new double[nx_];
		xConstrTermUB_ = new double[nx_];
		uConstr_       = new double[nu_];


		Alinear_  = new double[nx_*nx_]; // This is overwritten after Alin_ is computed
		Blinear_  = new double[nx_*nu_]; // This is overwritten after Blin_ is computed
		Clinear_  = new double[nx_*N_];  // This is NOT overwritten because we simply have that Clin = Clinear * dt_
				
		Alinear2_ = new double[nx_*nx_];
		Alinear3_ = new double[nx_*nx_];
		Alinear4_ = new double[nx_*nx_];

		
		Blin_ = new double[N_*nx_*nu_];

		AlinDelay_ = new double[nx_*nx_];
		BlinDelay_ = new double[nx_*nu_];

		Alin_ = new double[N_*nx_*nx_];
		
		AlinearOut = new double[nx_*nx_];
		BlinearOut = new double[nx_*nu_];
		ClinearOut = new double[nx_];
		
		xPred = new double[nx_*(N_+1)];
		uPred = new double[nu_*(N_)];

		xLin_ = new double[nx_*(N_+1)];
		uLin_ = new double[nu_*(N_)];

		for (int i =0; i<(N_+1); i++){
			for (int j = 0; j<nx_; j++){
				xLin_[i*nx_ + j] = 0.0;
			}
		}
		for (int i = 0; i<(N_+1); i++){
			for (int j = 0; j<nx_; j++){
				xLin_[i*nx_ + j] = linearization_IC_[j];
			}
		}
		
		for (int i =0; i<N_; i++){
			for (int j = 0; j<nu_; j++){
				uLin_[i*nu_ + j] = 0.0;
			}
		}			

		xCost_      = new float[nx_];
		xCostFinal_ = new float[nx_];
		uCost_      = new float[nu_];

		Qlin_        = new float[nx_];
		Qstate_rate  = new float[nx_];

 		// Initialize OSQP
		settings_ = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
		data_ = (OSQPData *)c_malloc(sizeof(OSQPData));

 	}

 	// Deallocation Memory
	MPCValFun::~MPCValFun(void)
	{
		osqp_cleanup(work_);
		if (data_)
		{
			c_free(data_->A);
			c_free(data_->P);
			c_free(data_);
		}
		if (settings_)  c_free(settings_);

		delete[] Hx_;
		delete[] Hr_;
		delete[] Hc_;
		delete[] lb_x_;
		delete[] ub_x_;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
 	/////////////////////////////////////// Define Functions /////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void MPCValFun::updateHorizon(){
		if (Nterm_ > 2){
			Nterm_ = Nterm_ - 1;
		}
		// cout << "Horizon Length: " << Nterm_ << endl;
	}

	void MPCValFun::resetHorizon(){
		Nterm_ = N_;
	}
	

	void MPCValFun::updateGoalSetAndState(double goalSetAndState[])
	{
		x_g_[0] = goalSetAndState[0];
		x_g_[1] = goalSetAndState[1];
		
		xConstrLB_[0] = goalSetAndState[2] - enlarge_;
		xConstrUB_[0] = goalSetAndState[3] + enlarge_;
		xConstrLB_[1] = goalSetAndState[4] - enlarge_;
		xConstrUB_[1] = goalSetAndState[5] + enlarge_;

		highLevTime   = goalSetAndState[6];

		xConstrTermLB_[0] = goalSetAndState[7];
		xConstrTermUB_[0] = goalSetAndState[8];
		xConstrTermLB_[1] = goalSetAndState[9];
		xConstrTermUB_[1] = goalSetAndState[10];

		// cout << "xConstrLB_[0]: " << xConstrLB_[0] << " xConstrUB_[0]: " << xConstrUB_[0] << endl;
		// cout << "xConstrLB_[1]: " << xConstrLB_[1] << " xConstrUB_[1]: " << xConstrUB_[1] << endl;

	}

	void MPCValFun::setGoalState(double x_goal[])
	{
		for (int i = 0; i < nx_; i++){
			x_g_[i] = x_goal[i];
		}
	}


		
	void MPCValFun::readCost()
	{
		// Read Cost Matrices
		ifstream myfile;
		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/Qrate_Qdiff.txt"), std::ios::app);
		myfile >> Qrate_;
		myfile >> Qdiff_;
		myfile.close();

		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/Qstate_rate.txt"), std::ios::app);
		for (int i = 0; i<nx_; i++){
			myfile >> Qstate_rate[i];
		}
		myfile.close();

		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/Qlin.txt"), std::ios::app);
		for (int i = 0; i<nx_; i++){
			myfile >> Qlin_[i];
		}
		myfile.close();
		
		// Read Cost Matrices
		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/Q.txt"), std::ios::app);
		for (int j = 0; j < nx_; j++) {
			myfile >> xCost_[j];
		}
		myfile.close();

		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/Qf.txt"), std::ios::app);
		for (int j = 0; j < nx_; j++) {
			myfile >> xCostFinal_[j];
		}
		myfile.close();

		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/R.txt"), std::ios::app);
		for (int j = 0; j < nu_; j++) {
			myfile >> uCost_[j];
		}
		myfile.close();
	}

	void MPCValFun::buildCost()
	{
		// IMPORTANT: the structure is fixed upper triangular!!!

		int Hc_counter = 0;
		int Hx_counter = 0;

		// Initialize cost for states
		Hc_[Hc_counter] = 0;
		Hc_counter += 1;
		// cout << " Writing H" << endl;
		for (int i = 0; i< (N_+1); i++){
			for (int j = 0; j< nx_; j++){
				// Rate cost off diagonal
				if (i > 0){
					Hx_[Hx_counter]     = -Qstate_rate[j];
					Hr_[Hx_counter]     = i*nx_ + j - nx_;
					Hx_counter+=1;
				}
				// Diagonal cost
				if (i == 0){
					Hx_[Hx_counter] = xCost_[j]      + Qlin_[j] + Qstate_rate[j] + Qerr[j];
					// cout << Hx_[Hx_counter] << endl;
				}else if ( i < Nterm_ ){
					Hx_[Hx_counter] = xCost_[j]      + Qlin_[j] + 2*Qstate_rate[j];
				}else if ( ( i >= Nterm_ ) and (i < N_ ) ){
					Hx_[Hx_counter] = xCostFinal_[j] + Qlin_[j] + 2*Qstate_rate[j];
				}else{				
					Hx_[Hx_counter] = xCostFinal_[j] + Qlin_[j] + Qstate_rate[j] ;
				}
				Hr_[Hx_counter]     = i*nx_ + j;
				Hx_counter+=1;

				Hc_[Hc_counter] = Hx_counter;
				Hc_counter += 1;

			}	
		}

		// Initialize cost for inputs
		for (int i = 0; i< N_; i++){
			for (int j = 0; j< nu_; j++){
				// Rate cost
				if (i > 0){
					Hx_[Hx_counter]     = -Qrate_;
					Hr_[Hx_counter]     = (N_ + 1)*nx_ + i*nu_ + j - nu_;
					Hx_counter+=1;
				}

				// Diff cost
				if (j==1){
					Hx_[Hx_counter]     = -Qdiff_;
					Hr_[Hx_counter]     = (N_ + 1)*nx_ + i*nu_ + j - 1;
					Hx_counter+=1;
				}

				// Diagonal cost
				if ( i==0 ) {
					Hx_[Hx_counter]     = uCost_[j] + Qdiff_ + 2*Qrate_;
				}else if ( i == N_ - 1){
					Hx_[Hx_counter]     = uCost_[j] + Qdiff_ + Qrate_;
				}else{
					Hx_[Hx_counter]     = uCost_[j] + Qdiff_ + 2*Qrate_;
				}
				Hr_[Hx_counter]     = (N_ + 1)*nx_ + i*nu_ + j;
				Hx_counter+=1;

				Hc_[Hc_counter] = Hx_counter;
				Hc_counter += 1;
			}	
		}

		// Hessian --> multiply by 2
		for (int i = 0; i < Hx_counter; i++){
			// cout << "Hx_ : " << Hx_[i] << endl; 
			Hx_[i] = 2*Hx_[i];
		}

		// cout << " Writing qv" << endl;
		// Set linear cost
		for (int i = 0; i < (N_ + 1); i++){
			for (int j = 0; j < nx_; j++ ){
				if ( i==0 ){
					qv_[i*nx_ + j] = -2*xCost_[j]*x_g_[j]      - 2*Qlin_[j]*xLin_[i*nx_ + j] - 2*Qerr[j]*x_IC_[j];
					// cout << qv_[i*nx_ + j] << endl;
				}else if (i < Nterm_ ){
					qv_[i*nx_ + j] = -2*xCost_[j]*x_g_[j]      - 2*Qlin_[j]*xLin_[i*nx_ + j];
				}else{				
					qv_[i*nx_ + j] = -2*xCostFinal_[j]*x_g_[j] - 2*Qlin_[j]*xLin_[i*nx_ + j];
				}
			} 
		}
		// cout << "Nterm_: " << Nterm_ << endl;
		// cout << "goal x: ";
		// for (int i = 0; i< nx_; i++){
		// 	cout << x_g_[i] << ", ";
		// }
		// cout << endl;

		// cout << "goal xCostFinal_: ";
		// for (int i = 0; i< nx_; i++){
		// 	cout << xCostFinal_[i] << ", ";
		// }
		// cout << endl;

		for (int i = 0; i < N_; i++){
			for (int j = 0; j < nu_; j++ ){
				if (i ==0){
					qv_[(N_+1)*nx_ + i*nu_ + j] = -2*Qrate_*oldInput[j];
				}else{
					qv_[(N_+1)*nx_ + i*nu_ + j] = 0.0;
				}
			} 
		}
		// for (int i = 0; i < nv_; i++){
		// 	cout << "qv_ : " << qv_[i] << endl; 
		// }

		// cout << "Initial Condition: ";
		// for (int i = 0; i<nx_; i++){
		// 	cout << x_IC_[i] << ", ";
		// }
		// cout << endl;
	}

	void MPCValFun::initiConstrVector()
	{
		// Read constraint vector
		ifstream myfile;
		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/x_constr.txt"), std::ios::app);
		// cout << "state constraint: ";
		for (int j = 0; j < nx_; j++) {
			myfile >> xConstrUB_[j];
			xConstrLB_[j] = -xConstrUB_[j];

			// cout << xConstrLB_[j];

			xConstrTermUB_[j] =  xConstrUB_[j];
			xConstrTermLB_[j] =  xConstrLB_[j];
		}
		// cout << endl;
		myfile.close();

		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/u_constr.txt"), std::ios::app);
		for (int j = 0; j < nu_; j++) {
			myfile >> uConstr_[j];
		}
		myfile.close();
	}

	void MPCValFun::buildConstrVector()
	{

		for (int i = 0; i< nx_; i++){
			lb_x_[i] = 0.0;
			ub_x_[i] = 0.0;
		}
		// Affine dynamics bound
		for (int i = 0; i< nx_*N_; i++){
			lb_x_[nx_ + i] = Clinear_[i]*dt_;
			ub_x_[nx_ + i] = Clinear_[i]*dt_;
		}

		// State box constraints
		for (int i = 0; i< (N_+1); i++){
			for (int j = 0; j< nx_; j++){

				if (i == 0){
					lb_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrLB_[j]-10; // No constraint on first predicted state
					ub_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrUB_[j]+10; // No constraint on fisrt predicted state					
				}
				else if (i < Nterm_){
					lb_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrLB_[j] + tightening[j];
					ub_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrUB_[j] - tightening[j];
				}else{
					lb_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrTermLB_[j] + tightening[j];
					ub_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrTermUB_[j] - tightening[j];					
				}
			}	
		}

		// Input box constraints
		for (int i = 0; i< N_; i++){
			for (int j = 0; j< nu_; j++){
				lb_x_[nx_*(N_+1) + nx_*(N_+1) + i*nu_ + j] = -uConstr_[j];
				ub_x_[nx_*(N_+1) + nx_*(N_+1) + i*nu_ + j] =  uConstr_[j];
			}	
		}

	}

	void MPCValFun::buildConstrMatrix()
	{
		int Fx_counter = 0;
		int Fc_counter = 0;
		int col = 0;
		int row = 0;

		// Note the matrix Fx = [F_equality_MPC_sparse_formulation, Identify_for_constraint_dim(nx_*N + nu_*N)]

		// Populate matrix related to state
		Fc_[Fc_counter] = 0;
		Fc_counter += 1;
		for (int i = 0; i < (N_ + 1); i++){
			for (int j = 0; j < nx_; j++){
				// Pick row col associated with identiy matrix
				col = i*nx_ + j;
				row = i*nx_ + j;

				// Add identity
				Fx_[Fx_counter] = 1;
				Fr_[Fx_counter] = row;
				Fx_counter += 1;		
				
				// Now matrix A
				if (i < N_){
					for (int ii = 0; ii < nx_; ii++){
						row = (i+1)*nx_ + ii;
						Fx_[Fx_counter] = -Alin_[i*(nx_*nx_) + ii*nx_ + j];
						Fr_[Fx_counter] = row;
						Fx_counter += 1;
					}					
				}
				// Add Identity for state constraints
				row = nx_*(N_ + 1) + col;
				// Add identity
				Fx_[Fx_counter] = 1;
				Fr_[Fx_counter] = row;
				Fx_counter += 1;

				Fc_[Fc_counter] = Fx_counter;
				Fc_counter += 1;
			}
		}

		// Populate matrix related to inputs
		for (int i = 0; i < N_; i++){
			for (int j = 0; j < nu_; j++){
				col = nx_*(N_ + 1) + i*nu_ + j;
				for (int ii = 0; ii < nx_; ii++){
					row = (i+1)*nx_ + ii;
					Fx_[Fx_counter] = -Blin_[i*(nx_*nu_) + ii*nu_ + j];
					Fr_[Fx_counter] = row;
					Fx_counter += 1;
				}

				// Add Identity for input constraints
				row = nx_*(N_+1) +  col;
				// Add identity
				Fx_[Fx_counter] = 1;
				Fr_[Fx_counter] = row;
				Fx_counter += 1;


				Fc_[Fc_counter] = Fx_counter;		
				Fc_counter += 1;
			}
		}

		if (printLevel_ >= 6){
			cout << "Fx: " << endl;
			int countern = 0;
			for (int i = 0; i < Fx_counter; i++){
				if (i > Fc_[countern+1]-1){
					countern = countern + 1;
				}
				cout << Fx_[i] << ", Row: " << Fr_[i] << ", Col: "<< countern << endl;
			}			
		}
	}


	void MPCValFun::readAB()
	{	

		// Read Matrix A
		ifstream myfile;
		myfile.open((matrix_prefix_path_+"/qp_matrices/planningModelMatrices/Alinear.txt"), std::ios::app);
	
		if (printLevel_ >= 1) cout << "START Reading Matrix A" << endl;
		for (int i = 0; i < nx_; i++) {
			for (int j = 0; j < nx_; j++) {
				myfile >> Alinear_[i*nx_ + j];
				if (printLevel_ >= 2) cout << Alinear_[i*nx_ + j] << ", ";
			}
			if (printLevel_ >= 2) cout << endl;
		}
		if (printLevel_ >= 1) cout << "DONE Reading Matrix A" << endl;
		myfile.close();

		// Read Matrix B (Simple Euler (integral of matrix exponential hard due to singularity))
		myfile.open((matrix_prefix_path_+"/qp_matrices/planningModelMatrices/Blinear.txt"), std::ios::app);
		
		if (printLevel_ >= 1) cout << "START Reading Matrix B" << endl;
		for (int i = 0; i < nx_; i++) {
			for (int j = 0; j < nu_; j++) {
				myfile >> Blinear_[i*nu_ + j];
				if (printLevel_ >= 2) cout << Blinear_[i*nu_ + j] << ", ";
			}
			if (printLevel_ >= 2) cout << endl;
		}
		myfile.close();
		if (printLevel_ >= 1) cout << "DONE Reading Matrix B" << endl;
	}
	
	void MPCValFun::linearize()
	{
		for (int idxN = 0; idxN < N_; idxN ++){
			// Get linearized dynamics
			linearizedDynamics_(xLin_,x_eq_, uLin_,Alinear_,Blinear_,Clinear_,idxN);

			if (idxN == 0){
 		    	for (int i = 0; i<nx_*nx_; i++){
 		    		AlinearOut[i] = Alinear_[i];
 		    	}
 		    	for (int i = 0; i<nx_*nu_; i++){
 		    		BlinearOut[i] = Blinear_[i];
 		    	}
 		    	for (int i = 0; i<nx_; i++){
 		    		ClinearOut[i] = Clinear_[i];
 		    	}
 		    }

			// Linearize x y
			matrixProd(Alinear2_, Alinear_, Alinear_);
			matrixProd(Alinear3_, Alinear2_, Alinear_);
			matrixProd(Alinear4_, Alinear3_, Alinear_);

			// Discretize Matrix A (approximate matrix Explonential)
			if (printLevel_ >= 4) cout << "Alin Approximation with dt_: "<< dt_ << " idxN: " << idxN << endl;
			for (int row = 0; row < nx_; row ++){
				for (int col = 0; col < nx_; col ++){
					if (row == col){
						Alin_[idxN*(nx_*nx_) + row*nx_ + col] = 1.0 + Alinear_[row*nx_ + col]*dt_ + 0.5*Alinear2_[row*nx_ + col]*pow(dt_,2) + (1.0/6.0)*Alinear3_[row*nx_ + col]*pow(dt_,3) + (1.0/24.0)*Alinear4_[row*nx_ + col]*pow(dt_,4);
					}else{
						Alin_[idxN*(nx_*nx_) + row*nx_ + col] =       Alinear_[row*nx_ + col]*dt_ + 0.5*Alinear2_[row*nx_ + col]*pow(dt_,2) + (1.0/6.0)*Alinear3_[row*nx_ + col]*pow(dt_,3) + (1.0/24.0)*Alinear4_[row*nx_ + col]*pow(dt_,4);
					}
					if (printLevel_ >= 4) cout << Alin_[idxN*(nx_*nx_) + row*nx_ + col] << " ";
				}
				if (printLevel_ >= 4) cout << endl;
			}

			if ((printLevel_ >= 3) and (idxN == 0)) cout << "Alin Approximation with dt_: "<< dt_ << " idxN: " << idxN << endl;
			for (int row = 0; row < nx_; row ++){
				for (int col = 0; col < nx_; col ++){
					if ((printLevel_ >= 3) and (idxN == 0)) cout << Alin_[idxN*(nx_*nx_) + row*nx_ + col] << " ";
				}
				if ((printLevel_ >= 3) and (idxN == 0)) cout << endl;
			}

			// Discretize Matrix B
			if (printLevel_ >= 4) cout << "Blin Approximation with dt_: "<< dt_ << " idxN: " << idxN << endl;
			for (int row = 0; row < nx_; row ++){
				for (int col = 0; col < nu_; col ++){
					Blin_[idxN*(nx_*nu_) + row*nu_ + col] = Blinear_[row*nu_ + col] * dt_;
					if (printLevel_ >= 4) cout << Blin_[idxN*(nx_*nu_) + row*nu_ + col] << " ";
				}
				if (printLevel_ >= 4) cout << endl;
			}

			if ((printLevel_ >= 3) and (idxN == 0)) cout << "Blin Approximation with dt_: "<< dt_ << " idxN: " << idxN << endl;
			for (int row = 0; row < nx_; row ++){
				for (int col = 0; col < nu_; col ++){
					if ((printLevel_ >= 3) and (idxN == 0)) cout << Blin_[idxN*(nx_*nu_) + row*nu_ + col] << " ";
				}
				if ((printLevel_ >= 3) and (idxN == 0)) cout << endl;
			}

			// if idxN == 0 ---> take into account delay
			if (idxN == 0){
				// Discretize Dealy Matrix propagation
				if (printLevel_ >= 4) cout << "AlinDelay_ Approximation with dtDelay_: "<< dtDelay_ << endl;
				for (int row = 0; row < nx_; row ++){
					for (int col = 0; col < nx_; col ++){
						if (row == col){
							AlinDelay_[row*nx_ + col] = 1.0 + Alinear_[row*nx_ + col]*dtDelay_ + 0.5*Alinear2_[row*nx_ + col]*pow(dtDelay_,2) + (1.0/6.0)*Alinear3_[row*nx_ + col]*pow(dtDelay_,3) + (1.0/24.0)*Alinear4_[row*nx_ + col]*pow(dtDelay_,4);
						}else{
							AlinDelay_[row*nx_ + col] =       Alinear_[row*nx_ + col]*dtDelay_ + 0.5*Alinear2_[row*nx_ + col]*pow(dtDelay_,2) + (1.0/6.0)*Alinear3_[row*nx_ + col]*pow(dtDelay_,3) + (1.0/24.0)*Alinear4_[row*nx_ + col]*pow(dtDelay_,4);
						}
						if (printLevel_ >= 4) cout << AlinDelay_[row*nx_ + col] << " ";
					}
					if (printLevel_ >= 4) cout << endl;
				}

				// Discretize Matrix B delay
				if (printLevel_ >= 4) cout << "BlinDelay_ Approximation with dtDelay_: "<< dtDelay_ << endl;
				for (int row = 0; row < nx_; row ++){
					for (int col = 0; col < nu_; col ++){
						BlinDelay_[row*nu_ + col] = Blinear_[row*nu_ + col] * dtDelay_;
						if (printLevel_ >= 4) cout << BlinDelay_[row*nx_ + col] << " ";
					}
					if (printLevel_ >= 4) cout << endl;
				}
			}
		}
	}

	void MPCValFun::matrixProd(double Aout[], double A1[], double A2[])
	{	
		// compute x_IC = Alin * xCurr + Blin * uCurr
		if (printLevel_ >= 4) cout << "==================================== matrixProd" << endl;
		for (int i=0; i < nx_; i++) {
			for (int j = 0; j < nx_; ++j){
				Aout[i*nx_ + j] = 0.0;
				for (int ii = 0; ii < nx_; ii ++)
				{
					Aout[i*nx_ + j] = Aout[i*nx_ + j] + A1[i*nx_ + ii] * A2[ii*nx_ + j];
				}
				if (printLevel_ >= 4) cout << Aout[i*nx_ + j] << ", ";
			}
			if (printLevel_ >= 4) cout << endl;
		}

		if (printLevel_>=4){
			cout << "A matrix" << endl;
			for (int i = 0; i<nx_; i++){
				for (int j = 0; j <nx_; j++){
					cout << A1[i*nx_ + j] << ",";
				}
				cout << endl;
			}

			cout << "A matrix" << endl;
			for (int i = 0; i<nx_; i++){
				for (int j = 0; j <nx_; j++){
					cout << A2[i*nx_ + j] << ",";
				}
				cout << endl;
			}			
		}

	}
	
	int MPCValFun::setUpOSQP(int verbose)
	{
		if (printLevel_ >= 1) cout << "START Setting-up Solver" << endl;
		cout << "START Setting-up Solver" << endl;
		// Setting default settings
		osqp_set_default_settings(settings_);
		settings_->max_iter = 2000;
		// settings_->polish = 1;
		// settings_->verbose = verbose;

		// Fill up osqp data
		data_->n = nv_;
		data_->m = nc_;
		data_->P = csc_matrix(data_->n, data_->n, H_nnz_, Hx_, Hr_, Hc_);
		data_->q = qv_;
		// data_->A = csc_matrix(data_->m, data_->n, FOld_znn_, FxOld_, FrOld_, FcOld_);
		data_->A = csc_matrix(data_->m, data_->n, F_znn_, Fx_, Fr_, Fc_);
		data_->l = lb_x_;
		data_->u = ub_x_;

		if (printLevel_ >= 1) cout << "DONE Setting-up Solver" << endl;
	
		osqp_setup(&work_, data_, settings_);

		cout << "HERE !!! DONE Setting-up Solver" << endl;

		return osqp_setup(&work_, data_, settings_);
	}

	void MPCValFun::setIC(double xt[])
	{
		for (int i = 0; i < nx_; i++){
			x_IC_[i] = xt[i];
		}
	}

	void MPCValFun::oneStepPrediction(double ut[])
	{	
		double xDummy[nx_] = {};
		// compute x_IC = Alin * xCurr + Blin * uCurr + Clinear *dt_
		for (int i=0; i < nx_; i++) {
			xDummy[i] = Clinear_[0*nx_ + i] * dt_;
			for (int j = 0; j < nx_; ++j){
				xDummy[i] = xDummy[i] + Alin_[0*(nx_*nx_) + i*nx_ + j] * x_IC_[j];
			}
			for (int j = 0; j < nu_; ++j){
				xDummy[i] = xDummy[i] + Blin_[0*(nx_*nu_) + i*nu_ + j] * ut[j];	
			}
		}

		for (int i = 0; i < nx_; i++){
			x_IC_[i] = xDummy[i];
		}		
	}

	void MPCValFun::oneStepPredictionDelay(double ut[])
	{	
		double xDummy[nx_] = {};
		// compute x_IC = Alin * xCurr + Blin * uCurr + Clinear *dtDelay_
		for (int i=0; i < nx_; i++) {
			xDummy[i] = Clinear_[0*nx_ + i] * dtDelay_;
			for (int j = 0; j < nx_; ++j){
				xDummy[i] = xDummy[i] + AlinDelay_[i*nx_ + j] * x_IC_[j];
			}
			for (int j = 0; j < nu_; ++j){
				xDummy[i] = xDummy[i] + BlinDelay_[i*nu_ + j] * ut[j];	
			}
		}

		for (int i = 0; i < nx_; i++){
			x_IC_[i] = xDummy[i];
		}		
	}

	void MPCValFun::solveQP()
	{
		// Update initial condition (x0 is not a variable so need to constraint x1 ---> Ax0 in lb_x_ and ub_x_)
		if (printLevel_ >= 3) cout << "Compute Ax0 to to update lb_x_ and ub_x_"<< nx_ << std::endl;
		for (int i = 0; i < nx_; i++) {
			lb_x_[i] = x_IC_[i]; - tightening[i];
			ub_x_[i] = x_IC_[i]; + tightening[i]; 
		}

		if (printLevel_ >= 3){
			for (int i = 0; i < nx_; i++) {
				cout << " Udpated Initial condition lb_x_[i] " << lb_x_[i] << std::endl;
				cout << " Udpated Initial condition ub_x_[i] " << ub_x_[i] << std::endl;
			}			
		}

		// data_->n = nv_;
		// data_->m = nc_;
		// data_->P = csc_matrix(data_->n, data_->n, H_nnz_, Hx_, Hr_, Hc_);
		// data_->q = qv_;
		// data_->A = csc_matrix(data_->m, data_->n, F_znn_, Fx_, Fr_, Fc_);

		// data_->l = lb_x_;
		// data_->u = ub_x_;
		// osqp_setup(&work_, data_, settings_);

		osqp_update_lower_bound(work_,lb_x_);
		osqp_update_upper_bound(work_,ub_x_);	
		osqp_update_lin_cost(work_,qv_);	
		osqp_update_A(work_, Fx_, OSQP_NULL, F_znn_);
		osqp_update_P(work_, Hx_, OSQP_NULL, H_nnz_);
		
		// Now solve QP
		osqp_solve(work_);
		if (printLevel_ >= 1) cout << "QP solved optimally: "<< work_->info->status_val << std::endl;
		solverFlag = work_->info->status_val;

		if (work_->info->status_val != OSQP_SOLVED){
			cout << "MPC not optimal, solver flag: " << work_->info->status_val << endl;
		}

		

		// Update xLin_
		if ( (work_->info->status_val == 1) or (work_->info->status_val == 2) ){	
			// Store optimal solution 
			for (int i = 0; i<(N_ + 1); i++){
				for (int j = 0; j < nx_; ++j){
					xPred[i*nx_+j] = static_cast<double>(work_->solution->x[i*nx_+j]);
				}
			}

			for (int i = 0; i<N_; i++){
				for (int j = 0; j < nu_; ++j){
					uPred[i*nu_+j] = static_cast<double>(work_->solution->x[nx_*(N_+1) + i*nu_+j]);
					if (i==0){
						oldInput[j] = uPred[i*nu_+j];
					}
				}
			}		


			for (int i = 0; i<(N_+1); i++){
				for (int j = 0; j < nx_; ++j){
					if (i < N_ ){
						xLin_[i*nx_ + j] = xPred[(i+1)*nx_+j];
					}else{
						xLin_[i*nx_ + j] = xPred[i*nx_+j];
					}
				}
			}

			// Update xLin_
			for (int i = 0; i<N_; i++){
				for (int j = 0; j < nu_; ++j){
					if (i < N_ - 1){
						uLin_[i*nu_ + j] = uPred[(i+1)*nu_+j];
					}else{
						uLin_[i*nu_ + j] = uPred[(i)*nu_+j];
					}
				}
			}
		}
		else{
			// Store optimal solution 
			for (int i = 0; i<(N_ + 1); i++){
				for (int j = 0; j < nx_; ++j){
					if (j < N_ -1){				
						xPred[i*nx_+j] = xPred[i*nx_+j+1];
					}else{
						xPred[i*nx_+j] = xPred[i*nx_+j];
					}
				}
			}

			for (int i = 0; i<N_; i++){
				for (int j = 0; j < nu_; ++j){
					if (i < N_-1){
						uPred[i*nu_+j] = uPred[i*nu_+j+1];
					}else{
						uPred[i*nu_+j] = uPred[i*nu_+j];
					}

				}
			}					
		}

		// cout << "Optimal state: ";
		// for (int i = 0; i<nx_; i++){
		// 	cout << xPred[i] << ", ";
		// }
		// cout << endl;
		// cout << "Initial Condition: ";
		// for (int i = 0; i<nx_; i++){
		// 	cout << x_IC_[i] << ", ";
		// }
		// cout << endl;
		// Print optimal solution

		if ((printLevel_ >= 3)){// or (work_->info->status_val >= 2)){
			cout << "OPTIMAL States:"<< std::endl;
			for (int i = 0; i< nx_; i++){
				for (int j = 0; j < N_+1; ++j){
					cout << xPred[j*nx_+i] <<",";
				}
				cout << endl;
			}	

			cout << "OPTIMAL Inputs:"<< std::endl;
			for (int i = 0; i< nu_; i++){
				for (int j = 0; j < N_; ++j){
					cout << uPred[j*nu_+i] <<",";
				}
				cout << endl;
			}	

			for (int i = 0; i < 3; i++){
				cout << static_cast<double>(work_->solution->x[nx_*N_ + nu_*N_ + i]) << ", i: " << i << endl;
			}	
		}

	}

}
#endif
