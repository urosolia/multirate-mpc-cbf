#ifndef MPC_H
#define MPC_H

#include <osqp.h>
#include <fstream>
#include <functional>
#include <cstring>
#include <string>
#include <vector>
#include <array>

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

enum class STATUS : uint8_t
{
	FAILURE = 0,
	RUNNING = 1
};

namespace ModelPredictiveControllerValFun
{
template<int nx_, int nu_, int N_>
class MPCValFunAbstract
{

 public:
  MPCValFunAbstract() {
    x_IC_ = new double[nx_];
    AlinearOut = new double[nx_*nx_];
    BlinearOut = new double[nx_*nu_];
    ClinearOut = new double[nx_];

    xPred = new double[nx_*(N_+1)];
    uPred = new double[nu_*(N_)];
  }

  ~MPCValFunAbstract() {
    delete [] x_IC_;
    delete [] AlinearOut;
    delete [] BlinearOut;
    delete [] ClinearOut;
    delete [] xPred;
    delete [] uPred;
  }

  virtual void setIC(double xt[]) = 0;
  virtual void oneStepPrediction(double ut[]) = 0;
  virtual void oneStepPredictionDelay(double ut[]) = 0;
  virtual void solveQP() = 0;
  virtual void readAB(void) = 0;
  virtual int setUpOSQP(int verbose) = 0;
  virtual void initiConstrVector(void) = 0;
  virtual void buildConstrVector(void) = 0;
  virtual void buildConstrMatrix(void) = 0;
  virtual void readCost(void) = 0;
  virtual void buildCost(void) = 0;
  virtual void linearize(void) = 0;
  virtual void updateHorizon(void) = 0;
  virtual void resetHorizon(void) = 0;
  virtual void setGoalState(double xt[]) = 0;
  virtual void updateGoalSetAndState(const std::array<double, 12> &goalSetAndState) = 0;

  double *xPred;
  double *uPred;
  double solverFlag = 0.0;
  double *AlinearOut;
  double *BlinearOut;
  double *ClinearOut;
  double *x_IC_;
};
	///////////////////////////// MPC                                              S
	using LinDyn = std::function<void(const double* /*x*/,
                                        const double* /*x_eq*/,
                                        const double* /*u*/,
                                        double* /*A*/,
                                        double* /*B*/,
                                        double* /*C*/,
                                        const int /*idxN*/)>;
  template<int nx_, int nu_, int N_, typename Linearizer = LinDyn>
	class MPCValFun : public MPCValFunAbstract<nx_, nu_, N_>
	{

	public:

		MPCValFun(
			const double dt,
			const double dtDelay,
			const uint32_t printFlag,
			const double x_eq[],
			const double max_error[],
			const double enlarge,
			const double lowLevelActive,
			const double linearization_IC[],
			const std::string matrix_prefix_path,
			Linearizer linearizedDynamics);

		~MPCValFun(void);

		void matrixProd(double Aout[], double A1[], double A2[]);
		virtual void setIC(double xt[]) override;
		virtual void oneStepPrediction(double ut[]) override;
		virtual void oneStepPredictionDelay(double ut[]) override;
		virtual void solveQP() override;
		virtual void readAB(void) override;
		virtual int setUpOSQP(int verbose) override;
		virtual void initiConstrVector(void) override;
		virtual void buildConstrVector(void) override;
		virtual void buildConstrMatrix(void) override;
		virtual void readCost(void) override;
		virtual void buildCost(void) override;
		virtual void linearize(void) override;
		virtual void updateHorizon(void) override;
		virtual void resetHorizon(void) override;
		virtual void setGoalState(double xt[]) override;
		virtual void updateGoalSetAndState(const std::array<double, 12> &goalSetAndState) override;

	protected:
		double *xLin_;
		double *uLin_;

		const double dt_;
		const double dtDelay_;
		const uint32_t  printLevel_;
		const double *x_eq_;
	    const double *max_error_;
		const double enlarge_ = {};
	    const double lowLevelActive_ = {};
		const double *linearization_IC_;
		const std::string matrix_prefix_path_;
    Linearizer linearizedDynamics_;
		int nv_;
		int nc_;
		double *x_g_;


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
using LinDynCov = std::function<void(const double * /*x*/,
                                     const double * /*x_eq*/,
                                     const double * /*u*/,
                                     double * /*A*/,
                                     double * /*B*/,
                                     double * /*C*/,
                                     double * /*K*/,
                                     const int /*idxN*/)>;

template<int nx_, int nu_, int N_, typename LINEARIZER=LinDynCov>
class MPCValFunGP : public MPCValFun<nx_, nu_, N_, LINEARIZER> {
 protected:
  double * Klinear_; // This is overwritten after Alin_ is computed
  double * Klin_;
 public:
  using PARENT = MPCValFun<nx_, nu_, N_, LINEARIZER>;
  MPCValFunGP(
    const double dt,
    const double dtDelay,
    const uint32_t printFlag,
    const double x_eq[],
    const double max_error[],
    const double enlarge,
    const double lowLevelActive,
    const double linearization_IC[],
    const std::string matrix_prefix_path,
    LINEARIZER linearizer) :
    PARENT(dt, dtDelay, printFlag, x_eq, max_error, enlarge,
           lowLevelActive, linearization_IC, matrix_prefix_path, linearizer) {
    Klinear_  = new double[nx_*nx_]; // This is overwritten after Alin_ is computed
    Klin_ = new double[N_*N_*nx_*nx_];
  }

  ~MPCValFunGP(){
    delete [] Klinear_;
    delete [] Klin_;
  }
  void buildConstrVector(void) override {
    PARENT ::buildConstrVector();
    // Affine dynamics bound without dt
    for (int i = 0; i< nx_*N_; i++){

      this->lb_x_[nx_ + i] = this->Clinear_[i];
      this->ub_x_[nx_ + i] = this->Clinear_[i];
    }
  }
  virtual void linearize(void) override {
    using AMat = Eigen::Matrix<double, nx_, nx_>;
    using BMat = Eigen::Matrix<double, nx_, nu_>;
    using CMat = Eigen::Matrix<double, nx_, 1>;
    using KMat = Eigen::Matrix<double, nx_ + nu_, nx_ + nu_>;

    for (int idxN = 0; idxN < N_; idxN ++){
      if constexpr (std::is_same<LINEARIZER, LinDynCov>::value){
        // Get linearized dynamics
        this->linearizedDynamics_(this->xLin_,
                                  this->x_eq_,
                                  this->uLin_,
                                  this->Alinear_,
                                  this->Blinear_,
                                  this->Clinear_,
                                  this->Klinear_,
                                  idxN);
      } else {
        return;
      }

      auto AlinearMat = Eigen::Map<AMat>(this->Alinear_);
      auto BlinearMat = Eigen::Map<BMat>(this->Blinear_);
      auto ClinearMat = Eigen::Map<CMat>(this->Clinear_);
      auto KlinearMat = Eigen::Map<KMat>(this->Klinear_);

      auto AlinearOutMat = Eigen::Map<AMat>(this->AlinearOut);
      auto BlinearOutMat = Eigen::Map<BMat>(this->BlinearOut);
      auto ClinearOutMat = Eigen::Map<CMat>(this->ClinearOut);

      if (idxN == 0){
        AlinearOutMat = AlinearMat;
        BlinearOutMat = BlinearMat;
        ClinearOutMat = ClinearMat;
      }

      // Discretize Matrix A (approximate matrix Explonential)
      auto AlinIdxNmat = Eigen::Map<AMat>(this->Alin_ + idxN*(nx_*nx_));
      AlinIdxNmat = AlinearMat;

      // Discretize Matrix B
      auto BlinIdxNmat = Eigen::Map<BMat>(this->Blin_ + idxN*(nx_*nu_));
      BlinIdxNmat = BlinearMat;

      auto KlinIdxMat = Eigen::Map<KMat>(this->Klin_ + idxN*idxN*nx_*nx_);
      // if idxN == 0 ---> take into account delay
      //TODO: No computationally cheap way to handle this case with GP
      //with GPs we get the exp(adt)
      if (this->printLevel_ >= 2) std::cout << "AlinearMat IdxN (" << idxN <<"):"  << std::endl << AlinearMat << std::endl;
      if (this->printLevel_ >= 2) std::cout << "BlinearMat IdxN (" << idxN <<"):" << std::endl << BlinearMat << std::endl;
      if (this->printLevel_ >= 2) std::cout << "ClinearMat IdxN (" << idxN <<"):" << std::endl << ClinearMat << std::endl;
      if (this->printLevel_ >= 2) std::cout << "KlinearMat IdxN (" << idxN <<"):" << std::endl << KlinIdxMat << std::endl;
    }
  }
};

using BatchLinDynCov = std::function<void(const double * /*x*/,
                                          const double * /*x_eq*/,
                                          const double * /*u*/,
                                          double * /*A*/,
                                          double * /*B*/,
                                          double * /*C*/,
                                          double * /*K*/)>;

template<int nx_, int nu_, int N_, typename LINEARIZER=BatchLinDynCov>
class MPCValFunGPBatch : public MPCValFunGP<nx_, nu_, N_, LINEARIZER> {
  using PARENT = MPCValFunGP<nx_, nu_, N_, LINEARIZER>;
 public:
  MPCValFunGPBatch(
    const double dt,
    const double dtDelay,
    const uint32_t printFlag,
    const double x_eq[],
    const double max_error[],
    const double enlarge,
    const double lowLevelActive,
    const double linearization_IC[],
    const std::string matrix_prefix_path,
    LINEARIZER linearizer) :
    PARENT(dt, dtDelay, printFlag, x_eq, max_error, enlarge,
           lowLevelActive, linearization_IC, matrix_prefix_path, linearizer) {}

  virtual void linearize(void) override {
    using AMat = Eigen::Matrix<double, nx_, nx_>;
    using BMat = Eigen::Matrix<double, nx_, nu_>;
    using CMat = Eigen::Matrix<double, nx_, 1>;
    using KMat = Eigen::Matrix<double, nx_ + nu_, nx_ + nu_>;

    this->linearizedDynamics_(this->xLin_,
                              this->x_eq_,
                              this->uLin_,
                              this->Alin_,
                              this->Blin_,
                              this->Clinear_,
                              this->Klin_);
    for (int idxN = 0; idxN < N_; idxN ++){

      // Matrix A
      auto AlinIdxNmat = Eigen::Map<AMat>(this->Alin_ + idxN*(nx_*nx_));
      // Matrix B
      auto BlinIdxNmat = Eigen::Map<BMat>(this->Blin_ + idxN*(nx_*nu_));
      // Matrix C
      auto ClinIdxMat = Eigen::Map<CMat>(this->Clinear_ + idxN*nx_);
      // Matrix K
      auto KlinIdxMat = Eigen::Map<KMat>(this->Klin_ + idxN*(nx_ + nu_)*(nx_ + nu_));

      if (this->printLevel_ >= 2) std::cout << "AlinearMat IdxN (" << idxN <<"):"  << std::endl << AlinIdxNmat << std::endl;
      if (this->printLevel_ >= 2) std::cout << "BlinearMat IdxN (" << idxN <<"):" << std::endl << BlinIdxNmat << std::endl;
      if (this->printLevel_ >= 2) std::cout << "ClinearMat IdxN (" << idxN <<"):" << std::endl << ClinIdxMat << std::endl;
      if (this->printLevel_ >= 2) std::cout << "KlinearMat IdxN (" << idxN <<"):" << std::endl << KlinIdxMat << std::endl;
    }
  }
};

	// Define constructor
  template<int nx_, int nu_, int N_, typename Linearizer>
	MPCValFun<nx_,nu_,N_, Linearizer>::MPCValFun(
    const double dt,
    const double dtDelay,
    const uint32_t printFlag,
    const double x_eq[],
    const double max_error[],
    const double enlarge,
    const double lowLevelActive,
    const double linearization_IC[],
    const std::string matrix_prefix_path,
    Linearizer linearizedDynamics):
    dt_(dt), dtDelay_(dtDelay), printLevel_(printFlag), x_eq_(x_eq),
    max_error_(max_error), enlarge_(enlarge), lowLevelActive_(lowLevelActive), linearization_IC_(linearization_IC),
    matrix_prefix_path_(matrix_prefix_path), linearizedDynamics_(linearizedDynamics)
	{
					
		// Initialize horizon
		Nterm_ = N_;

		// Tightening paramters
    std::cout << "[MPC] Low Level Active: " << lowLevelActive_ << std::endl;
		std::cout << "[MPC] Enlarge: " << enlarge << std::endl;

    std::cout << "[MPC] linearization_IC_: " << linearization_IC_[2] << std::endl;

		tightening = new double[nx_]{};
		Qerr = new float[nx_]{};
		oldInput = new double[nu_]{};

		if (lowLevelActive_ > 0.5){
			for (int i = 0; i < nx_; i++)
				tightening[i] = max_error_[i];
      std::ifstream myfile;
			myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/Qerr.txt"), std::ios::app);
			for (int i = 0; i < nx_; i++){
				myfile >> Qerr[i];
			}
			myfile.close();
		}
		std::cout << "Tightening Paramters: " << std::endl;
		std::cout << "error x: "        << tightening[0] << std::endl;
		std::cout << "error y: "        << tightening[1] << std::endl;
		std::cout << "error theta: "    << tightening[2] << std::endl;
		std::cout << "error v: "        << tightening[3] << std::endl;
		std::cout << "error thetaDot: " << tightening[4] << std::endl;
		std::cout << "error psi: "      << tightening[5] << std::endl;
		std::cout << "error psiDot: "   << tightening[6] << std::endl;
		for (int i = 0; i < nx_; i++){
			std::cout << "error cost: "     << Qerr[i] << " i: " << i << std::endl;
		}

		// Initialize Constraint
		nv_ = nx_*(N_+1) + nu_*N_; // Total number of variable
		nc_ = nx_*(N_+1) + nx_*(N_+1) + nu_*N_; // Total number of constraints: nx_*(N+1) equality constraints + (nx_*(N_+1) + nu_*N_) box constraints
		
		F_znn_ = nx_*(N_+1) + nx_*nx_*N_ + nx_*nu_*N_ + nx_*(N_+1) + nu_*N_;
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

		x_g_  = new double[nx_];

		xConstrLB_     = new double[nx_];
		xConstrUB_     = new double[nx_];
		xConstrTermLB_ = new double[nx_];
		xConstrTermUB_ = new double[nx_];
		uConstr_       = new double[nu_];


		Alinear_  = new double[nx_*nx_]; // This is overwritten after Alin_ is computed
		Blinear_  = new double[nx_*nu_]; // This is overwritten after Blin_ is computed
		Clinear_  = new double[nx_*N_];  // This is NOT overwritten because we simply have that Clin = Clinear * dt_

		
		Blin_ = new double[N_*nx_*nu_];

		AlinDelay_ = new double[nx_*nx_];
		BlinDelay_ = new double[nx_*nu_];

		Alin_ = new double[N_*nx_*nx_];

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

		Qlin_        = new float[nx_]; // Nnx*N
		Qstate_rate  = new float[nx_];

 		// Initialize OSQP
		settings_ = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
		data_ = (OSQPData *)c_malloc(sizeof(OSQPData));

 	}

 	// Deallocation Memory
  template<int nx_, int nu_, int N_, typename Linearizer>
  MPCValFun<nx_,nu_,N_, Linearizer>::~MPCValFun(void)
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
  template<int nx_, int nu_, int N_, typename Linearizer>
  void MPCValFun<nx_,nu_,N_, Linearizer>::updateHorizon(){
		if (Nterm_ > 2){
			Nterm_ = Nterm_ - 1;
		}
		//std::cout << "Horizon Length: " << Nterm_ <<  std::endl;
	}
template<int nx_, int nu_, int N_, typename Linearizer>
 void MPCValFun<nx_,nu_,N_, Linearizer>::resetHorizon(){
		Nterm_ = N_;
	}

template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::updateGoalSetAndState(const std::array<double, 12> &goalSetAndState)
	{
		x_g_[0] = goalSetAndState[0];
		x_g_[1] = goalSetAndState[1];
		x_g_[2] = goalSetAndState[2];

		xConstrLB_[0] = goalSetAndState[3] - enlarge_;
		xConstrUB_[0] = goalSetAndState[4] + enlarge_;
		xConstrLB_[1] = goalSetAndState[5] - enlarge_;
		xConstrUB_[1] = goalSetAndState[6] + enlarge_;

		highLevTime   = goalSetAndState[7];

		xConstrTermLB_[0] = goalSetAndState[8];
		xConstrTermUB_[0] = goalSetAndState[9];
		xConstrTermLB_[1] = goalSetAndState[10];
		xConstrTermUB_[1] = goalSetAndState[11];

		//std::cout << "xConstrLB_[0]: " << xConstrLB_[0] << " xConstrUB_[0]: " << xConstrUB_[0] <<  std::endl;
		//std::cout << "xConstrLB_[1]: " << xConstrLB_[1] << " xConstrUB_[1]: " << xConstrUB_[1] <<  std::endl;

	}
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::setGoalState(double x_goal[])
	{
		for (int i = 0; i < nx_; i++){
			x_g_[i] = x_goal[i];
		}
	}


template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::readCost()
	{
		// Read Cost Matrices
		std::ifstream myfile;
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
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::buildCost()
	{
		// IMPORTANT: the structure is fixed upper triangular!!!

		int Hc_counter = 0;
		int Hx_counter = 0;

		// Initialize cost for states
		Hc_[Hc_counter] = 0;
		Hc_counter += 1;
		//std::cout << " Writing H" <<  std::endl;
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
					//std::cout << Hx_[Hx_counter] <<  std::endl;
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
			//std::cout << "Hx_ : " << Hx_[i] <<  std::endl;
			Hx_[i] = 2*Hx_[i];
		}

		//std::cout << " Writing qv" <<  std::endl;
		// Set linear cost
		for (int i = 0; i < (N_ + 1); i++){
			for (int j = 0; j < nx_; j++ ){
				if ( i==0 ){
					qv_[i*nx_ + j] = -2*xCost_[j]*x_g_[j]      - 2*Qlin_[j]*xLin_[i*nx_ + j] - 2*Qerr[j]*this->x_IC_[j];
					//std::cout << qv_[i*nx_ + j] <<  std::endl;
				}else if (i < Nterm_ ){
					qv_[i*nx_ + j] = -2*xCost_[j]*x_g_[j]      - 2*Qlin_[j]*xLin_[i*nx_ + j];
				}else{				
					qv_[i*nx_ + j] = -2*xCostFinal_[j]*x_g_[j] - 2*Qlin_[j]*xLin_[i*nx_ + j];
				}
			} 
		}
		//std::cout << "Nterm_: " << Nterm_ <<  std::endl;
		//std::cout << "goal x: ";
		// for (int i = 0; i< nx_; i++){
		// 	std::cout << x_g_[i] << ", ";
		// }
		//std::cout <<  std::endl;

		//std::cout << "goal xCostFinal_: ";
		// for (int i = 0; i< nx_; i++){
		// 	std::cout << xCostFinal_[i] << ", ";
		// }
		//std::cout <<  std::endl;

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
		// 	std::cout << "qv_ : " << qv_[i] <<  std::endl;
		// }

		//std::cout << "Initial Condition: ";
		// for (int i = 0; i<nx_; i++){
		// 	std::cout << x_IC_[i] << ", ";
		// }
		//std::cout <<  std::endl;
	}
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::initiConstrVector()
	{
		// Read constraint vector
		std::ifstream myfile;
		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/x_constr.txt"), std::ios::app);
		for (int j = 0; j < nx_; j++) {
			myfile >> xConstrUB_[j];
			xConstrLB_[j] = -xConstrUB_[j];

			xConstrTermUB_[j] =  xConstrUB_[j];
			xConstrTermLB_[j] =  xConstrLB_[j];
		}
		myfile.close();

		myfile.open((matrix_prefix_path_+"/qp_matrices/tuning/u_constr.txt"), std::ios::app);
		for (int j = 0; j < nu_; j++) {
			myfile >> uConstr_[j];
		}
		myfile.close();
	}
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::buildConstrVector()
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
					lb_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrLB_[j]-10;
					ub_x_[nx_*(N_+1) + i*nx_ + j] =  xConstrUB_[j]+10;					
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
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::buildConstrMatrix()
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
				Fr_[Fx_counter] = row;				Fx_counter += 1;
				
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
		 	std::cout << "Fx: " <<  std::endl;
			int countern = 0;
			for (int i = 0; i < Fx_counter; i++){
				if (i > Fc_[countern+1]-1){
					countern = countern + 1;
				}
			 	std::cout << Fx_[i] << ", Row: " << Fr_[i] << ", Col: "<< countern <<  std::endl;
			}			
		}
	}

template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::readAB()
	{	

		// Read Matrix A
		std::ifstream myfile;
		myfile.open((matrix_prefix_path_+"/qp_matrices/planningModelMatrices/Alinear.txt"), std::ios::app);
	
		if (printLevel_ >= 1)std::cout << "START Reading Matrix A" <<  std::endl;
		for (int i = 0; i < nx_; i++) {
			for (int j = 0; j < nx_; j++) {
				myfile >> Alinear_[i*nx_ + j];
				if (printLevel_ >= 2)std::cout << Alinear_[i*nx_ + j] << ", ";
			}
			if (printLevel_ >= 2)std::cout <<  std::endl;
		}
		if (printLevel_ >= 1)std::cout << "DONE Reading Matrix A" <<  std::endl;
		myfile.close();

		// Read Matrix B (Simple Euler (integral of matrix exponential hard due to singularity))
		myfile.open((matrix_prefix_path_+"/qp_matrices/planningModelMatrices/Blinear.txt"), std::ios::app);
		
		if (printLevel_ >= 1) std::cout << "START Reading Matrix B" <<  std::endl;
		for (int i = 0; i < nx_; i++) {
			for (int j = 0; j < nu_; j++) {
				myfile >> Blinear_[i*nu_ + j];
				if (printLevel_ >= 2) std::cout << Blinear_[i*nu_ + j] << ", ";
			}
			if (printLevel_ >= 2) std::cout <<  std::endl;
		}
		myfile.close();
		if (printLevel_ >= 1) std::cout << "DONE Reading Matrix B" <<  std::endl;
	}

template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::linearize()
  {
		for (int idxN = 0; idxN < N_; idxN ++){
			// Get linearized dynamics
			if constexpr (std::is_same<Linearizer, LinDyn>::value){
        linearizedDynamics_(xLin_, x_eq_, uLin_, Alinear_, Blinear_, Clinear_, idxN);
			} else {
			  return;
			}

      using AMat = Eigen::Matrix<double, nx_, nx_>;
      using BMat = Eigen::Matrix<double, nx_, nu_>;
      using CMat = Eigen::Matrix<double, nx_, 1>;

      auto AlinearMat = Eigen::Map<AMat>(Alinear_);
      auto BlinearMat = Eigen::Map<BMat>(Blinear_);
      auto ClinearMat = Eigen::Map<CMat>(Clinear_);

      auto AlinearOutMat = Eigen::Map<AMat>(this->AlinearOut);
      auto BlinearOutMat = Eigen::Map<BMat>(this->BlinearOut);
      auto ClinearOutMat = Eigen::Map<CMat>(this->ClinearOut);

			if (idxN == 0){
			    AlinearOutMat = AlinearMat;
			    BlinearOutMat = BlinearMat;
			    ClinearOutMat = ClinearMat;
 		    }

			// Discretize Matrix A (approximate matrix Explonential)
			auto AlinIdxNmat = Eigen::Map<AMat>(Alin_ + idxN*(nx_*nx_));
			AlinIdxNmat = (AlinearMat * dt_).exp();

			// Discretize Matrix B
			auto BlinIdxNmat = Eigen::Map<BMat>(Blin_ + idxN*(nx_*nu_));
			BlinIdxNmat = BlinearMat * dt_;

      if (printLevel_ >= 2) std::cout << "AlinearMat IdxN (" << idxN <<"):"  << std::endl << AlinearMat << std::endl;
      if (printLevel_ >= 2) std::cout << "BlinearMat IdxN (" << idxN <<"):" << std::endl << BlinearMat << std::endl;
      if (printLevel_ >= 2) std::cout << "ClinearMat IdxN (" << idxN <<"):" << std::endl << ClinearMat << std::endl;

      // if idxN == 0 ---> take into account delay
      if (idxN == 0){
        // Discretize Dealy Matrix propagation
        if (printLevel_ >= 4) std::cout << "AlinDelay_ Approximation with dtDelay_: "<< dtDelay_ <<  std::endl;
        auto AlinDelayMat = Eigen::Map<AMat>(AlinDelay_);
        AlinDelayMat = (AlinearMat * dtDelay_).exp();

        // Discretize Matrix B delay
        if (printLevel_ >= 4) std::cout << "BlinDelay_ Approximation with dtDelay_: "<< dtDelay_ <<  std::endl;
        auto BlinDelayMat = Eigen::Map<BMat>(BlinDelay_);
        BlinDelayMat = BlinearMat * dtDelay_;

      }
		}
	}

template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::matrixProd(double Aout[], double A1[], double A2[])
	{	
		// compute x_IC = Alin * xCurr + Blin * uCurr
		if (printLevel_ >= 4) std::cout << "==================================== matrixProd" <<  std::endl;
		for (int i=0; i < nx_; i++) {
			for (int j = 0; j < nx_; ++j){
				Aout[i*nx_ + j] = 0.0;
				for (int ii = 0; ii < nx_; ii ++)
				{
					Aout[i*nx_ + j] = Aout[i*nx_ + j] + A1[i*nx_ + ii] * A2[ii*nx_ + j];
				}
				if (printLevel_ >= 4) std::cout << Aout[i*nx_ + j] << ", ";
			}
			if (printLevel_ >= 4) std::cout <<  std::endl;
		}

		if (printLevel_>=4){
		 	std::cout << "A matrix" <<  std::endl;
			for (int i = 0; i<nx_; i++){
				for (int j = 0; j <nx_; j++){
				 	std::cout << A1[i*nx_ + j] << ",";
				}
			 	std::cout <<  std::endl;
			}

		 	std::cout << "A matrix" <<  std::endl;
			for (int i = 0; i<nx_; i++){
				for (int j = 0; j <nx_; j++){
				 	std::cout << A2[i*nx_ + j] << ",";
				}
			 	std::cout <<  std::endl;
			}			
		}

	}
template<int nx_, int nu_, int N_, typename Linearizer>
int MPCValFun<nx_,nu_,N_, Linearizer>::setUpOSQP(int verbose)
	{
		if (printLevel_ >= 1) std::cout << "START Setting-up Solver" <<  std::endl;
	 	std::cout << "START Setting-up Solver" <<  std::endl;
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

		if (printLevel_ >= 1) std::cout << "DONE Setting-up Solver" <<  std::endl;
	
		osqp_setup(&work_, data_, settings_);

	 	std::cout << "HERE !!! DONE Setting-up Solver" <<  std::endl;

		return osqp_setup(&work_, data_, settings_);
	}
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::setIC(double xt[])
	{
		for (int i = 0; i < nx_; i++){
			this->x_IC_[i] = xt[i];
		}
	}
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::oneStepPrediction(double ut[])
	{	
		double xDummy[nx_] = {};
		// compute x_IC = Alin * xCurr + Blin * uCurr + Clinear *dt_
		for (int i=0; i < nx_; i++) {
			xDummy[i] = Clinear_[0*nx_ + i] * dt_;
			for (int j = 0; j < nx_; ++j){
				xDummy[i] = xDummy[i] + Alin_[0*(nx_*nx_) + i*nx_ + j] * this->x_IC_[j];
			}
			for (int j = 0; j < nu_; ++j){
				xDummy[i] = xDummy[i] + Blin_[0*(nx_*nu_) + i*nu_ + j] * ut[j];	
			}
		}

		for (int i = 0; i < nx_; i++){
			this->x_IC_[i] = xDummy[i];
		}		
	}
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::oneStepPredictionDelay(double ut[])
	{	
		double xDummy[nx_] = {};
		// compute x_IC = Alin * xCurr + Blin * uCurr + Clinear *dtDelay_
		for (int i=0; i < nx_; i++) {
			xDummy[i] = Clinear_[0*nx_ + i] * dtDelay_;
			for (int j = 0; j < nx_; ++j){
				xDummy[i] = xDummy[i] + AlinDelay_[i*nx_ + j] * this->x_IC_[j];
			}
			for (int j = 0; j < nu_; ++j){
				xDummy[i] = xDummy[i] + BlinDelay_[i*nu_ + j] * ut[j];	
			}
		}

		for (int i = 0; i < nx_; i++){
			this->x_IC_[i] = xDummy[i];
		}		
	}
template<int nx_, int nu_, int N_, typename Linearizer>
void MPCValFun<nx_,nu_,N_, Linearizer>::solveQP()
	{
		// Update initial condition (x0 is not a variable so need to constraint x1 ---> Ax0 in lb_x_ and ub_x_)
		if (printLevel_ >= 3) std::cout << "Compute Ax0 to to update lb_x_ and ub_x_"<< nx_ << std::endl;
		for (int i = 0; i < nx_; i++) {
			lb_x_[i] = this->x_IC_[i];
			ub_x_[i] = this->x_IC_[i];
		}

		if (printLevel_ >= 3){
			for (int i = 0; i < nx_; i++) {
			 	std::cout << " Udpated Initial condition lb_x_[i] " << lb_x_[i] << std::endl;
			 	std::cout << " Udpated Initial condition ub_x_[i] " << ub_x_[i] << std::endl;
			}			
		}

		// for (int j = 0; j< nu_; j++){
		// 	lb_x_[nx_*(N_+1) + nx_*(N_+1) + j] = oldInput[j] - 3.0;
		// 	ub_x_[nx_*(N_+1) + nx_*(N_+1) + j] = oldInput[j] + 3.0;
		// }	

		data_->n = nv_;
		data_->m = nc_;
		data_->P = csc_matrix(data_->n, data_->n, H_nnz_, Hx_, Hr_, Hc_);
		data_->q = qv_;
		// data_->A = csc_matrix(data_->m, data_->n, FOld_znn_, FxOld_, FrOld_, FcOld_);
		data_->A = csc_matrix(data_->m, data_->n, F_znn_, Fx_, Fr_, Fc_);

		data_->l = lb_x_;
		data_->u = ub_x_;
		c_int exit_flag = osqp_setup(&work_, data_, settings_);

		assert(exit_flag == 0); // Something is wrong with the initialization if this is the case.
		// Now solve QP
		osqp_solve(work_);
		if (printLevel_ >= 1) std::cout << "QP solved optimally: "<< work_->info->status_val << std::endl;
		this->solverFlag = work_->info->status_val;

		if (work_->info->status_val != OSQP_SOLVED){
		 	std::cout << "MPC not optimal, solver flag: " << work_->info->status_val <<  std::endl;
		}

		

		// Update xLin_
		if ( (work_->info->status_val == 1) or (work_->info->status_val == 2) ){	
			// Store optimal solution 
			for (int i = 0; i<(N_ + 1); i++){
				for (int j = 0; j < nx_; ++j){
					this->xPred[i*nx_+j] = static_cast<double>(work_->solution->x[i*nx_+j]);
				}
			}

			for (int i = 0; i<N_; i++){
				for (int j = 0; j < nu_; ++j){
					this->uPred[i*nu_+j] = static_cast<double>(work_->solution->x[nx_*(N_+1) + i*nu_+j]);
					if (i==0){
						oldInput[j] = this->uPred[i*nu_+j];
					}
				}
			}		


			for (int i = 0; i<(N_+1); i++){
				for (int j = 0; j < nx_; ++j){
					if (i < N_ ){
						xLin_[i*nx_ + j] = this->xPred[(i+1)*nx_+j];
					}else{
						xLin_[i*nx_ + j] = this->xPred[i*nx_+j];
					}
				}
			}

			// Update xLin_
			for (int i = 0; i<N_; i++){
				for (int j = 0; j < nu_; ++j){
					if (i < N_ - 1){
						uLin_[i*nu_ + j] = this->uPred[(i+1)*nu_+j];
					}else{
						uLin_[i*nu_ + j] = this->uPred[(i)*nu_+j];
					}
				}
			}
		}
		else{
			// Store optimal solution 
			for (int i = 0; i<(N_ + 1); i++){
				for (int j = 0; j < nx_; ++j){
					if (j < N_ -1){				
						this->xPred[i*nx_+j] = this->xPred[i*nx_+j+1];
					}else{
						this->xPred[i*nx_+j] = this->xPred[i*nx_+j];
					}
				}
			}

			for (int i = 0; i<N_; i++){
				for (int j = 0; j < nu_; ++j){
					if (i < N_-1){
						this->uPred[i*nu_+j] = this->uPred[i*nu_+j+1];
					}else{
						this->uPred[i*nu_+j] = this->uPred[i*nu_+j];
					}

				}
			}					
		}

		// std::cout << "Optimal state: ";
		// for (int i = 0; i<nx_; i++){
		// 	std::cout << xPred[i] << ", ";
		// }
		// std::cout <<  std::endl;
		// std::cout << "Initial Condition: ";
		// for (int i = 0; i<nx_; i++){
		// 	std::cout << x_IC_[i] << ", ";
		// }
		// std::cout <<  std::endl;
		// Print optimal solution

		if ((printLevel_ >= 3)){// or (work_->info->status_val >= 2)){
		 	std::cout << "OPTIMAL States:"<< std::endl;
			for (int i = 0; i< nx_; i++){
				for (int j = 0; j < N_+1; ++j){
				 	std::cout << this->xPred[j*nx_+i] <<",";
				}
			 	std::cout <<  std::endl;
			}	

		 	std::cout << "OPTIMAL Inputs:"<< std::endl;
			for (int i = 0; i< nu_; i++){
				for (int j = 0; j < N_; ++j){
				 	std::cout << this->uPred[j*nu_+i] <<",";
				}
			 	std::cout <<  std::endl;
			}	

			for (int i = 0; i < 3; i++){
			 	std::cout << static_cast<double>(work_->solution->x[nx_*N_ + nu_*N_ + i]) << ", i: " << i <<  std::endl;
			}	
		}

	}

}
#endif
