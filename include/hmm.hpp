#pragma once

#include "utils.hpp"

using namespace std;

namespace cdHMM {

template <typename obs_type> class multiHMM; //forward declaration so hmm members can access protected

template <typename obs_type> class HMM{
friend class multiHMM<obs_type>;
public:

	int M; /*!< Number of emission model parameters */
	int N; /*!< Number of Hidden States */
	int T; /*!< Number of Observations */
	int maxIters; /*!< maximum allowed number of EM iterations */
	int minIters; /*!< minimum allowed number of EM iterations */
	 
	bool print_iter; /*!< Show Likelihood info every EM iteration */
	
	//fix all or some of the transition/emission/starting probabilities
	bool fixA; 			/*!< Fix all values of transmission matrix A*/
	vector<bool> fixArow;	/*!< Fix values of transmission matrix A where element == true*/
	bool fixB; 				/*!< Fix all values of emission matrix B*/
	vector<bool> fixBrow;	/*!< Fix values of emission matrix B where element == true*/
	bool fixpi; 			/*!< Fix all values of starting vector pi*/
	vector<bool> fixpirow;	/*!< Fix values of starting vector pi where element == true*/
	
	default_random_engine generator; /*!< Random number generator */

	//Type of emission model
	string type;				/*!< Model type ID string*/

	//Best Parameters
	vector< vector<double> > maxA; /*!< Most likely transmission matrix */ 
	vector< vector<double> > maxB; /*!< Most likely emission matrix */ 
	vector<double> maxpi;			/*!< Most likely starting vector */ 
	double maxlhood;				/*!< Max likelihood  */ 

	//best state sequence
	vector<unsigned> state;		/*!< Most likely state sequence */ 

	//Some distributions are bounded below
	vector<double> lb;	/*!<  Emission lower bound */
	
protected:	
	//Parameters
	vector< vector<double> > A; /*!< Transmission matrix */ 
	vector< vector<double> > B; /*!< Emission matrix */
	vector<double> pi;			/*!< Starting vector */
		
	vector< vector<double> > logA; /*!< Log transmission matrix */ 
	vector< vector<double> > logB; /*!< Log emission matrix */ 
	vector<double> logpi;			/*!< Log starting vector */ 

	//Forward-Backward variables
	vector< vector<double> > alpha;	/*!< Baum-Welch forward probability */
	vector< vector<double> > beta;	/*!< Baum-Welch backward probability */
	vector< vector<double> > gamma;	/*!< Baum-Welch lhood at t */
	vector< vector<double> > delta; /*!< Viterbi lhood at t*/
	vector< vector<double> > delta2; /*!< Viterbi bookkeeping*/
	vector< double > sumgamma; /*!< Baum-Welch normalization */
	vector< double > sumgammap; /*!< Baum-Welch sumgammap = sumgamma - last term */
	vector< double > c; /*!< Baum-Welch forward-backward normalization */
	vector< vector< vector<double> > > digamma; /*!<  Baum-Welch i,j transmission at t */

	//model likelihood
	double lhood;	/*!<  Model likelihood */
	
	//log of observations
	bool set_logO;	/*!<  Already calculated log of observations */
	bool no_logO;   /*!<  Cannot calculate log of observations */
	vector<double> logO;	/*!< Log of observations */

	double logsmall; /*!< Log of a small number (avoid log(0)) */

	max_lhood_params mlp; /*!<  Parameters for max likelihood estimate solvers */

public:
/*!
 * Sets emission distribution lower bound (if applicable),
 * @param[in]  in_lb  The vector of lower bounds (one per hidden state).
 */ void set_lower_bound(vector<double> &in_lb){ lb = in_lb; }
/*!
 * Sets max likelihood parameter estimate precision 
 * @param[in]  eps solver precision.
 */	void set_optimise_precision(double eps){ mlp.eps = eps; }
/*!
 * Sets max likelihood parameter estimate pmax iteration 
 * @param[in]  max_iter solver max iteration count.
 */ void set_optimise_max_iter(double max_iter){ mlp.max_iter = max_iter; }
/*!
 * Sets max likelihood parameter estimate starting vector 
 * @param[in]  start solver starting point.
 */ void set_optimise_start(vector<double> start){ mlp.x_start = start; }
	
	
/*!
 * Sets transmission matrix A 
 * @param[in]  inA Transmission matrix A.
 */	void setA(vector<vector<double> > &inA){ 
		if( inA.size() != N ){ cerr << "Setting A with mis-sized vector" << endl; exit(1); }
		for(int i=0; i<N; ++i){
			if( inA[i].size() != N ){ cerr << "Setting A with mis-sized vector " << endl; exit(1); }
		}
		A = inA; 
    }
/*!
 * Gets transmission matrix A 
 */ vector< vector<double> > getA(){ return A; }
/*!
 * Gets max lhood transmission matrix A 
 */ vector< vector<double> > getmaxA(){ return maxA; }
 //Virtual so that generic HMM has a setB but can overwrite in multiHMM
/*!
 * Sets emission matrix B
 * @param[in]  inB Emission matrix B.
 */ virtual void setB(vector<vector<double> > &inB){ 
		if( inB.size() != N ){ cerr << "Setting B with mis-sized vector" << endl; exit(1); }
		for(int i=0; i<N; ++i){
			if( inB[i].size() != M ){ cerr << "Setting B with mis-sized vector " << endl; exit(1); }
		}
		B = inB; 
	}
/*!
 * Gets emission matrix B
 */vector< vector<double> > getB(){ return B; }
/*!
 * Gets max lhood emission matrix A 
 */vector< vector<double> > getmaxB(){ return maxB; }
    
/*!
 * Sets starting vector pi
 * @param[in]  starting vector inpi.
 */	void setpi( vector<double> &inpi){ 
		if( inpi.size() != N ){ cerr << "Setting A with mis-sized vector" << endl; exit(1); }
		pi = inpi; 
	}
/*!
 * Gets starting vector pi
 */    vector<double> getpi(){ return pi; }
/*!
 * Gets max lhood starting vector pi 
 */    vector<double> getmaxpi(){ return maxpi; }
    
//set current params to best params
/*!
 * Set current transmission, emission and starting vector to max lhood values.
 */ void set_to_max(){
		A = maxA;
		B = maxB;
		pi = maxpi;
		lhood = maxlhood;
	}
/*!
 * Clear max lhood values.
 */ void clear_max(){
		maxA = A;
		maxB = B;
		maxpi = pi;
		maxlhood = -numeric_limits<double>::infinity();	
	}
/*!
 * Print info about this model to error console
 */ virtual void info() = 0;
	
protected:
/*!
 * Initialize 
 * @param[in] N_ Number of hidden states.
 * @param[in] M_ Number of emission parameters.
 */ void setsize(int N_, int M_){
		M = M_;
		N = N_;
		
		//fixing parameters
		fixArow = vector<bool>(N, false);
		fixBrow = vector<bool>(N, false);
		fixpirow = vector<bool>(N, false);
		
		//parameters
		A = vector< vector<double> >( N, vector<double>(N, 1 ) );
		B = vector< vector<double> >( N, vector<double>(M, 1 ) );
		pi = vector<double>(N, 1);
		//log parameters
		logA = vector< vector<double> >( N, vector<double>(N, 1 ) );
		logB = vector< vector<double> >( N, vector<double>(M, 1 ) );
		logpi = vector<double>(N, 1);
	}
	
/*!
 * Set minimum and maximum iterations
 * @param[in] min_ Min number EM iterations.
 * @param[in] max_ Max number EM iterations.
 */	void setIters(int min_, int max_){
		maxIters = max_;
		minIters = min_;
		print_iter = false;
		maxlhood = -numeric_limits<double>::infinity();
		//default optimisation parameters
		mlp.dim=1;
		mlp.x_start = {1e-10, 1e5};
		mlp.eps=1e-5;
		mlp.max_iter=1000;
		//fixing parameters
		fixA = false;
		fixB = false;
		fixpi = false;
		//log obs info		
		set_logO = false;
		no_logO = false;
		logO.resize(0);
		//log(small)
		logsmall = -1e12;
	}

/*!
 * Initialise emission matrix
 */ virtual void initB() = 0;
		
/*!
 * Calculate log of A
 */	void calc_logA(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){ logA[i][j] = log(A[i][j]);}
		}
	}
/*!
 * Calculate log of B
 */	void calc_logB(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<M; ++j){ logB[i][j] = log(B[i][j]);}
		}
	}
/*!
 * Calculate log of pi
 */	void calc_logpi(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){ logpi[i] = log(pi[i]);}
		}
	}
/*!
 * Calculate log of model params
 */ void calc_logP(){
		calc_logA();
		calc_logB();
		calc_logpi();
	}
/*!
 * Calculate log of the observation
 * @param[in] O vector of observations.
 */	void calc_logO(vector<obs_type> &O){
		logO.resize( O.size() );
		for(int i=0; i<O.size(); ++i){ 
			if( O[i] <= 0 ){ 
				//cerr << "<=0 observation at " << i << " = " << O[i] << "\tCan't compute log" << endl;  
				no_logO = true; 
				break;
			}
			logO[i] = log(O[i]); 
		}
		set_logO = true;
	}
public:
	//set A, B, pi to arbitrary initial values
/*!
 * Set parameters to random initial starting values
 */ void init(){	
		
		set_logO = false;
		no_logO = false;
		logO.resize(0);
		
		double npi = 0;
		for(unsigned i=0; i<N; ++i){
			double norm = 0;
			for(unsigned j=0; j<N; ++j){ 
				A[i][j] = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
				norm += A[i][j]; 
			} 
			for(int j=0; j<N; ++j){ A[i][j] /= norm; }
			
			pi[i] = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
			npi += pi[i]; 
		}
		for(unsigned i=0; i<N; ++i){ pi[i] /= npi; }
		initB();
		
		lhood = -numeric_limits<double>::infinity();
	
		calc_logP();
	}
		
	//Default
/*!
 * Default construction
 */	HMM(){
		setsize(0,0);
		setIters(0,0);
	}
	
/*!
 * Constructor
 * @param[in] N_ number of hidden states
 * @param[in] M_ number of emission parameters
 * @param[in] min_ minimum number of EM iterations
 * @param[in] max_ maximum number of EM iterations
 */	HMM(int N_, int M_, int min_, int max_){
		setsize(N_,M_);
		setIters(min_,max_);
	}
	
protected:	
/*!
 * Re-order states according to no_obs.second
 * @param[in] no_obs order vector.
 * @param[in] fA Transmission matrix.
 * @param[in] fB Emission matrix.
 * @param[in] fpi Starting vector.
 */ void sortparams( vector< pair<double, int> > &no_obs, 
	vector<vector<double> > &fA, 
	vector<vector<double> > &fB, 
	vector<double> &fpi ){
		
		vector< vector<double> > At(N, vector<double>(N) );
		vector< vector<double> > Bt(N, vector<double>(M) );
		vector<double> pit(N);
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){
				At[i][j] = fA[no_obs[i].second][no_obs[j].second];
			} 
			for(int j=0; j<M; ++j){
				Bt[i][j] = fB[no_obs[i].second][j];
			} 
			pit[i] = fpi[no_obs[i].second];
		}
		fA = At;
		fB = Bt;
		fpi = pit;
	}
public:

/*!
 * re-order states so B[0][0] > B[1][0] > B[2][0] > ...
 */	void sortparams(){ 
		
		vector< pair<double, int> > no_obs(N);
		for(int i=0; i<N; ++i){ 
				no_obs[i].first = B[i][0];
				no_obs[i].second = i;
		}
		sort(no_obs.begin(), no_obs.end(), pair_sort_first<double, int>);
		sortparams(no_obs, A, B, pi);

		if( maxB.size() == N ){
			for(int i=0; i<N; ++i){ 
					no_obs[i].first = maxB[i][0];
					no_obs[i].second = i;
			}
			sort(no_obs.begin(), no_obs.end(), pair_sort_first<double, int>);
			sortparams(no_obs, maxA, maxB, maxpi);
		}
	}
protected:
	

/*!
 * reserve space for viterbi
 * @param[in] T number of observations
 */ inline void setup_delta(int T_){	
		T = T_;
		delta = vector<vector<double> >(N, vector<double>(T,0)); 
		delta2 = vector<vector<double> >(N, vector<double>(T,1));
		
		state = vector<unsigned>(T,0);
	}
	
/*!
 * reserve space for Baum-Welch
 * @param[in] T number of observations
 */	inline void setup(int T_){	
		
		T = T_;
		
		alpha = vector<vector<double> >(N, vector<double>(T,0) ); 
		
		beta = vector<vector<double> >(N, vector<double>(T,1) ); 
		gamma = vector<vector<double> >(N, vector<double>(T,0) ); 
		sumgamma = vector<double>(N,0); 
		sumgammap = vector<double>(N,0); 
		
		state = vector<unsigned>(T,0);
		
		digamma = vector< vector< vector<double> > >( N, vector<vector<double> >(N, vector<double>(T-1,0) ) );
		c = vector<double>(T,0);
		
	}

/*!
 * reserve space for Baum-Welch
 * @param[in] T vector of observations
 */ inline void setup(vector<obs_type> &O){	
		setup( O.size() );
	}
public:
/*!
* Probability of emitting O in state i
@param[in] i HMM state
@param[in] O observation
*/ virtual double pB(int i, obs_type O) = 0;
/*!
* Log probability of emitting O in state i
@param[in] i HMM state
@param[in] O observation
*/ 	virtual double log_pB(int i, obs_type O) = 0;	
	
protected:	

/*!
 * First forward probability
 * @param[in] O observation
 */inline void compute_alpha0(obs_type O){
		c[0] = 0;
		for(int i=0; i<N; ++i){
			alpha[i][0] = pi[i] * pB(i,O); 
			c[0] += alpha[i][0];
		}
		c[0] = 1.0/c[0];
		for(int i=0; i<N; ++i){ 
			alpha[i][0] *= c[0]; 
		}
	}
/*!
 * Updates forward probabilities. assumes alpha is big enough! 
 * @param[in] O observation
 * @param[in] t time
 */	inline void update_alpha(obs_type O, int t){
		double sum = 0;
		for(int i=0; i<N; ++i){
			alpha[i][t] = alpha[0][t-1] * A[0][i]; 
			for(int j=1; j<N; ++j){
				alpha[i][t] += alpha[j][t-1] * A[j][i]; 
			}
			alpha[i][t] *=  pB(i,O); 
			sum += alpha[i][t];
		}
		c[t] = 1.0/sum;
		for(int i=0; i<N; ++i){ alpha[i][t] *= c[t]; }
	}
/*!
 * Fill alpha from observation vector
 * @param[in] O observation
 */	inline void compute_alpha(vector<obs_type> &O){	
		compute_alpha0(O[0]); 
		for(int t=1; t<T; ++t){ update_alpha(O[t], t);	} //exit(1);
	}
/*!
 * Fill beta (backward probabilities) from observation vector	
 * @param[in] O observation
 */	inline void compute_beta(vector<obs_type> &O){
		for(int i=0; i<N; ++i){ beta[i][ T-1] = c[T-1];}
		for(int t=T-1; t>0; --t){			
			for(int i=0; i<N; ++i){
				
				beta[i][t-1] = A[i][0] * pB( 0, O[t] ) * beta[0][t];
				for(int j=1; j<N; ++j){
					beta[i][t-1] += A[i][j] * pB( j, O[t] ) * beta[j][t];
				}
				beta[i][t-1] *=  c[t-1];
			}
		}
	}
/*!
 * fill gamma and digamma (state probabilities) from observation vector		
 * @param[in] O observation
 */	inline void compute_gamma(vector<obs_type> &O){
		
		double sden, sdenl;
		for(int t=0; t<T-1; ++t){
			
			sden = 0;
			for(int i=0; i<N; ++i){
					
				gamma[i][t] = 0;
				sdenl = 0;
				for(int j=0; j<N; ++j){
					digamma[i][j][t] = alpha[i][t] * A[i][j] * pB( j, O[t+1] ) * beta[j][ t+1 ];
					gamma[i][t] += digamma[i][j][t];
					sdenl += alpha[j][t] * A[j][i]; 
				}
				sdenl *= pB( i, O[t+1] );
				sden += sdenl * beta[i][ t+1 ];			
			}
			for(int i=0; i<N; ++i){
				for(int j=0; j<N; ++j){ digamma[i][j][t] /= sden; }
				gamma[i][t] /= sden;
			}
	
		}
		sden = 0;
		for(int i=0; i<N; ++i){
			gamma[i][ T-1] = alpha[i][ T-1] * beta[i][ T-1];
			sden += gamma[i][T-1];
		}
		for(int i=0; i<N; ++i){ gamma[i][ T-1] /= sden; }
		
	}
/*!
 * compute gamma norm	
 */	inline void compute_sumgamma(){ 
		for(int i=0; i<N; ++i){
			sumgamma[i] = 0;
			for(int t=0; t<T; ++t){
				sumgamma[i] += gamma[i][t];
			}
		}
	}
/*!
 * First log forward probability		
 * @param[in] O observation
 */	inline void compute_log_alpha0(obs_type O){
		for(int i=0; i<N; ++i){
			alpha[i][0] = logpi[i] + log_pB(i,O); 
		}
	}
/*!
 * Log forward probability		
 * @param[in] O observation
 * @param[in] t time
 */	inline void update_log_alpha(obs_type O, int t){
		double sum;
		for(int i=0; i<N; ++i){
			
			alpha[i][t] = alpha[0][t-1] + logA[0][i]; 
			for(int j=1; j<N; ++j){
				lsum(alpha[i][t], alpha[j][t-1] + logA[j][i]);
			}
			alpha[i][t] += log_pB(i,O); 
		}
	}
/*!
 * fill forward log probabilities from observation vector		
 * @param[in] O observation
 */	inline void compute_log_alpha(vector<obs_type> &O){	
		compute_log_alpha0(O[0]);  		
		for(int t=1; t<T; ++t){ update_log_alpha(O[t], t);	} 
	}
/*!
 * fill backwards log probabilities from observation vector		
 * @param[in] O observation
 */	inline void compute_log_beta(vector<obs_type> &O){
		for(int i=0; i<N; ++i){ beta[i][ T-1] = 0; } 
		
		for(int t=T-1; t>0; --t){			
			for(int i=0; i<N; ++i){
				
				beta[i][t-1] = logA[i][0] + log_pB( 0, O[t] ) + beta[0][t];
				for(int j=1; j<N; ++j){
					lsum( beta[i][t-1] , (logA[i][j] + log_pB( j, O[t] ) + beta[j][t]) );
				}
			}
		}
	}
/*!
 * Fill state transition probabilities from observation vector	
 * @param[in] O observation
 */	inline void compute_log_gamma(vector<obs_type> &O){
		
		double sden, sdenl;
		for(int t=0; t<T-1; ++t){
			
			sden = 0;
			for(int i=0; i<N; ++i){
					
				sdenl = 0;
				for(int j=0; j<N; ++j){
					digamma[i][j][t] = alpha[i][t] + logA[i][j] + log_pB( j, O[t+1] ) + beta[j][ t+1 ];
					if(j==0){ 
						gamma[i][t] = digamma[i][j][t];
						sdenl = alpha[j][t] + logA[j][i]; 
					} else { 
						lsum(gamma[i][t] , digamma[i][j][t]);
						lsum(sdenl , (alpha[j][t] + logA[j][i])); 
					}

				}
				sdenl += log_pB( i, O[t+1] );
				if(i==0){ sden = sdenl + beta[i][ t+1 ]; } else { lsum(sden, sdenl + beta[i][ t+1 ]); }			
			}
			sden *= -1.0;
			for(int i=0; i<N; ++i){
				for(int j=0; j<N; ++j){ 
					digamma[i][j][t] += sden; 
				}
				gamma[i][t] += sden;
			}
	
		}
		for(int i=0; i<N; ++i){
			gamma[i][ T-1] = alpha[i][ T-1] +  beta[i][ T-1];
			if(i==0){ sden = gamma[i][T-1]; } else { lsum(sden, gamma[i][T-1]); } 
		}
		sden *= -1.0;
		for(int i=0; i<N; ++i){ 
			gamma[i][ T-1] += sden; 
		}
		
	}
	
/*!
 * Compute norm of transition probabilities
 */	inline void compute_log_sumgamma(){ 
		for(int i=0; i<N; ++i){
			for(int t=0; t<T; ++t){
				if( t==0 ){ sumgamma[i] = gamma[i][t]; } else { lsum(sumgamma[i] , gamma[i][t]); }
				if( t == T-2 ){ sumgammap[i] = sumgamma[i]; }
			}
		}
	}
	

		
/*!
 * Reestimate starting vector	
 */	inline void reestimate_pi(){
		for(int i=0; i<N; ++i){ 
			if( !fixpirow[i] ){
				pi[i] = gamma[i][0]; 
			}
		}
	}
/*!
 * Reestimate log starting vector	
 */	inline void reestimate_log_pi(){
		for(int i=0; i<N; ++i){
			if( !fixpirow[i] ){
				pi[i] = exp(gamma[i][0]); 
				logpi[i] = gamma[i][0]; 
			}
		}
	}
	
/*!
 * Reestimate transition matrix
 */	inline void reestimate_A(){
		double sum = 0; 
		for( int i=0; i<N; ++i ){ 
			if( !fixArow[i] ){
				for( int j=0; j<N; ++j ) { 
					sum = 0;
					for(int t=0; t<T-1; ++t){
						sum += digamma[i][j][t];
					}
					A[i][j] = sum/(sumgamma[i] - gamma[i][T-1]); 
				}
			}
		}
	}	
/*!
 * Reestimate log transition matrix
 */	inline void reestimate_log_A(){
		double sum = 0; 
		for( int i=0; i<N; ++i ){
			if( !fixArow[i] ){
				for( int j=0; j<N; ++j ) { 
					for(int t=0; t<T-1; ++t){
						if(t==0){ sum = digamma[i][j][t]; } else { lsum(sum, digamma[i][j][t]); }
					}
					double denom = sumgammap[i]; 
					logA[i][j] = sum - denom;
					A[i][j] = exp( logA[i][j] );
				}
			}
		}
	}
	
/*!
 * Reestimate emission matrix	
 */	virtual void reestimate_B(vector<obs_type> &O) = 0;
/*!
 * Reestimate log emission matrix	
 */	virtual void reestimate_log_B(vector<obs_type> &O) = 0;
	
/*!
 * Reestimate model
 */	inline void reestimate(vector<obs_type> &O){ 
		if( !fixpi ){ reestimate_pi(); }
		if( !fixA ){ reestimate_A(); }
		if( !fixB ){ reestimate_B(O); }
	}
/*!
 * Reestimate log model
 */	inline void reestimate_log(vector<obs_type> &O){ 
		if( !fixpi ){ reestimate_log_pi(); }
		if( !fixA ){ reestimate_log_A(); }
		if( !fixB ){ reestimate_log_B(O); }
	}
/*!
 * calc model likelihood
 */	void evaluate(){
		lhood = 0;
		for(int t=0; t<T; ++t){ lhood -= log( c[t] ); }
	}
/*!
 * calc model likelihood from log parameters
 */	void evaluate_log(){
		lhood = 0;
		for(int i=0; i<N; ++i){ 
			if(i==0){ lhood = alpha[i][T-1]; } else { lsum(lhood, alpha[i][T-1]); }; 
		}
	}
public:	
/*!
* calc model likelihood for observation seq O
* @param[in] O observations
*/ void evaluate(vector<obs_type> &O){
		
		T = O.size();
		alpha = vector<vector<double> >(N, vector<double>(T,0) ); 	
		c = vector<double>(T,0); 
		compute_alpha(O);
		evaluate();
		
	}
protected:

/*!
 * Start dynamic programming matrix
 * @param[in] O observation
 */	inline void compute_delta0(obs_type O){
		for(int i=0; i<N; ++i){ 
			delta[i][0] = logpi[i] + log_pB(i,O);
			delta2[i][0] = 0;
		}	
	 }
/*!
 * Update dynamic programming matrix
 * @param[in] O observation
 * @param[in] t time
 */	inline void update_delta(obs_type O, int t){
		double tmp, tmp2;
		for(int i=0; i<N; ++i){
			tmp = delta[0][t-1] + logA[0][i] + log_pB(i,O);
			delta2[i][t] = 0;
			for(int j=1; j<N; ++j){ 
				tmp2 = delta[j][t-1] + logA[j][i] + log_pB(i,O);
				if( tmp2 > tmp ){ 
					tmp = tmp2; 
					delta2[i][t] = j;
				}
			}
			delta[i][t] = tmp;
		}
	}
/*!
 * Dynamic programming state sequence
 * @param[in] last Final observation time
 */	inline void backtrack_viterbi(int last){
		double best = delta[0][last-1];
		state[ last-1 ] = 0;
		for(int i=1; i<N; ++i){ 
				if( delta[i][last-1] > best ){ best = delta[i][last-1]; state[ last-1 ] = i; }
		}
		for(int t=last-1; t>0; --t){
			state[t-1] = delta2[state[t]][ t]; 
		}
	}
public:
/*!
* calc most likely state sequence for observation sequence O using Viterbi algorithm. Call set_max first to use max lhood estimates
* for A, B, pi.
* @param[in] O observations
*/ void set_state_viterbi(vector<obs_type> &O){
		
		setup_delta(O.size());
		compute_delta0(O[0]);	
		for(int t=1; t<O.size(); ++t){ 
			update_delta(O[t], t);
		}
		backtrack_viterbi(O.size());
		
	}

	

protected:
/*!
* Choose most likely state at each time
 */	void set_state(){
		for(int t=0; t<T; ++t){ 
			state[t] = 0;
			double best = alpha[0][t];
			for(int i=1; i<N; ++i){
				if(alpha[i][t] > best){ best = alpha[i][t]; state[t] = i; }			
			}
		}
	}
public:
/*!
* calc most likely state sequence for observation sequence O using DP algorithm. Call set_max first to use max lhood estimates
* for A, B, pi.
* @param[in] O observations
*/ void set_state(vector<obs_type> &O){
		setup(O);
		compute_alpha(O);
		set_state();
	}
	
	
/*!
* print A to stdout
*/ 	void printA(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){
				cout << A[i][j] << " ";
			} cout << "\n";
		}
	}
/*!
* print A to file
* @param[in] ofile file stream
*/  void printA(ofstream &ofile){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){
				ofile << A[i][j] << " ";
			} ofile << "\n";
		}
	}
	
/*!
* print B to stdout
*/	void printB(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<B[i].size(); ++j){
				cout << B[i][j] << " ";
			} cout << "\n";
		}
	}
/*!
* print B to file
* @param[in] ofile file stream
*/  void printB(ofstream &ofile){
		for(int i=0; i<N; ++i){
			for(int j=0; j<B[i].size(); ++j){
				ofile << B[i][j] << " ";
			} ofile << "\n";
		}
	}

/*!
* print pi to stdout
*/	void printpi(){
		for(int i=0; i<N; ++i){
				cout << pi[i] << " ";
		} cout << "\n";
	}
/*!
* print pi to file
* @param[in] ofile file stream
*/  void printpi(ofstream &ofile){
		for(int i=0; i<N; ++i){
				ofile << pi[i] << " ";
		} ofile << "\n";
	}

	
/*!
* Fit a HMM to observation sequence. Returns { num_iters, max_log_lhood }
* @param[in] O observation sequence
* @param[in] eps EM stopping precision
* @param[in] uselog If false try to avoid using log probabilities, can have numerical issues!
*/  vector<double> fit(vector<obs_type> &O, double eps, bool uselog=true){
		
		setup(O);
		
		double old_lP = -numeric_limits<double>::infinity();
		bool done = false;
		
		vector< vector<double> > tmpA;
		vector< vector<double> > tmpB;
		vector<double> tmppi;
		
		int iters = 0;
		do{
			if( uselog ){

				compute_log_alpha(O);
				compute_log_beta(O); 
				compute_log_gamma(O); 
				
				//re-estimate	
				tmpA = A;
				tmpB = B;
				tmppi = pi;
			
				compute_log_sumgamma();
				for(int i=0; i<N; ++i){
					if( sumgamma[i] == 0 ){
						cerr << "problem with state " << i << endl;
						for(int t=0; t<T; ++t){
								cout << t << " " << i << " " << 
								alpha[i][t] << " " << beta[i][t] << " " << gamma[i][t] << " " << sumgamma[i] << endl;
						}
						cout << "A "; printA();
						cout << "pi "; printpi();
						cout << "B "; printB();
						cout << endl;
						exit(1);
					}
				}

				reestimate_log(O);
				
				evaluate_log();
								
			} else {

				compute_alpha(O); 
				compute_beta(O);
				compute_gamma(O);
				
				//re-estimate	
				tmpA = A;
				tmpB = B;
				tmppi = pi;
			
				compute_sumgamma();		
				

				
						
				reestimate(O);
								
				evaluate();
								
			}
			
			++iters;
		    if( print_iter ){ 
				cerr << "Iter: " << iters << " Lhood: " << lhood << endl; 
			}

			if( lhood != lhood ){ 
				cerr << "lhood = nan!" << endl; 
				for(int t=0; t<T; ++t){
					for(int i=0; i<N; ++i){
						cout << t << " " << i << " " << 
						alpha[i][t] << " " << beta[i][t] << " " << gamma[i][t] << " " << sumgamma[i] << endl;
					}
				}
				cout << "A "; printA();
				cout << "pi "; printpi();
				cout << "B "; printB();
				cout << endl;
			exit(1); }


			if(lhood > maxlhood){
				maxlhood = lhood;
				maxA = A;
				maxB = B;
				maxpi = pi;
			}

			double dA = -1;	//max change in A matrix
			double dB = -1; //max change in B matrix
			double dpi = -1;  //max change in pi matrix
			for(unsigned i=0; i<N; ++i){ 
				for(unsigned j=0; j<N; ++j){ 
					double diff = fabs(A[i][j] - tmpA[i][j]);
					if( diff > dA){ dA = diff; }
				} 
				for(unsigned j=0; j<M; ++j){ 
					double diff = fabs(B[i][j] - tmpB[i][j]);
					if( diff > dB){ dB = diff; }
				} 
				double diff = fabs(pi[i] - tmppi[i]);
				if( diff > dpi){ dpi = diff; }
			}
						
			if( iters < minIters ){	//at least minIters done
				old_lP = lhood;
			} else if(!done && 
			iters < maxIters && //at most maxIters
			fabs(lhood - old_lP) > eps && //lhood change
			(dA > eps || dB > eps || dpi > eps ) //param change
			){
				old_lP = lhood;
			} else {
				done = true;
			}
			
		} while (!done);

		
		//most likely state
		set_state();
		
		vector<double> ret = { (double) iters, old_lP };
		return ret;
	}
protected:
/*!
* Random initial state
*/  discrete_distribution<int> setup_distpi(){
		discrete_distribution<int> pi_dist(pi.begin(), pi.end());
		return pi_dist;
	}				
/*!
* Random transition
*/  vector< discrete_distribution<int> > setup_distA(){
		vector< discrete_distribution<int> > A_dist;
		for(int i=0; i<N; ++i){
			vector<double> a_b(N);
			for(int j=0; j<N; ++j){
				a_b[j] = A[i][j];
			}
			discrete_distribution<int> b_dist(a_b.begin(), a_b.end());
			A_dist.push_back( b_dist );
		}
		return A_dist;
	}
	
/*!
* Random emission
* @param[in] state Emission from this state
*/  virtual obs_type gen_obs(int state_) = 0;

public:	
/*!
* Generate observation sequence from HMM
* @param[out] O observation sequence
* @param[in] T number of observations to generate
* @param[in] print If true print O to stdout
*/ 	void generate_seq(vector<obs_type> &O, int T, bool print = false){	

		O.resize(T);
		//initial state
		discrete_distribution<int> pi_dist = setup_distpi(); 
		//Transitions
		vector< discrete_distribution<int> > A_dist = setup_distA();

		//run
		int cur_state = pi_dist(generator);
		for(int t=0; t<T; ++t){
			//generate observation
			O[t] = gen_obs(cur_state);
			if(print){ cout << O[t] << "\n"; }
			//transition
			cur_state = A_dist[cur_state](generator);
		}
  
	}
	
};

//crap to make the linker happy
template class HMM<double>;
template class HMM<float>;
template class HMM<int>;

}
