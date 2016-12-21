#pragma once

#include "utils.hpp"

using namespace std;

namespace cdHMM {

template <typename obs_type> class HMM{
public:

	int M;
	int N;
	int T;
	int maxIters;
	int minIters;
	
	bool print_iter;
	
	//fix all or some of the transition/emission/starting probabilities
	bool fixA; vector<bool> fixArow;
	bool fixB; vector<bool> fixBrow;
	bool fixpi; vector<bool> fixpirow;
	
	default_random_engine generator;

	//Type of emission model
	string type;
	
	//Parameters
	vector< vector<double> > A, B;
	vector<double> pi;

	//Best Parameters
	vector< vector<double> > maxA, maxB;
	vector<double> maxpi;
	double maxlhood;
		
	//Parameters
	vector< vector<double> > logA, logB;
	vector<double> logpi;
	
	//best state sequence
	vector<unsigned> state;

	//Forward-Backward variables
	vector< vector<double> > alpha, beta, gamma;	//Baum-Welch variables
	vector< vector<double> > delta, delta2; //Viterbi variables
	vector< double > sumgamma, c, sumgammap; //sumgammap = sumgamma - last term
	vector< vector< vector<double> > > digamma;

	//model likelihood
	double lhood;
	
	//log of observations
	bool set_logO;
	bool no_logO;
	vector<double> logO;

	//Some distributions are bounded below
	vector<double> lb;

	//Some distributions need to get max likelihood parameters numerically
	max_lhood_params mlp;
	
	//So generic HMM has a setB method without needing to reimplement it in ecery realization
	virtual void setB(vector<vector<double> > &inB){ B = inB; }

	//Initialize Num states = N and num observations = M
	void setsize(int N_, int M_){
		M = M_;
		N = N_;
	}
	
	//Set minimum double of iterations and maximum double
	void setIters(int min_, int max_){
		maxIters = max_;
		minIters = min_;
		print_iter = false;
		maxlhood = -numeric_limits<double>::infinity();
		//default optimisation parameters
		mlp.dim=1;
		mlp.x_start = {1e-10, 1e5};
		mlp.eps=1e-5;
		mlp.max_iter=1000;
	}

	virtual void info() = 0;
	virtual void initB() = 0;
		
	void calc_logA(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){ logA[i][j] = log(A[i][j]);}
		}
	}
	void calc_logB(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<M; ++j){ logB[i][j] = log(B[i][j]);}
		}
	}
	void calc_logpi(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){ logpi[i] = log(pi[i]);}
		}
	}
	void calc_logP(){
		logA = vector< vector<double> >( N, vector<double>(N, 1 ) );
		logB = vector< vector<double> >( N, vector<double>(M, 1 ) );
		logpi = vector<double>(N, 1);
		calc_logA();
		calc_logB();
		calc_logpi();
	}
	void calc_logO(vector<obs_type> &O){
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
	
	//set A, B, pi to arbitrary initial values
	void init(){	
		fixA = false;
		fixB = false;
		fixpi = false;
		
		fixArow = vector<bool>(N, false);
		fixBrow = vector<bool>(N, false);
		fixpirow = vector<bool>(N, false);
		
		set_logO = false;
		no_logO = false;
		logO.resize(0);
		
		A = vector< vector<double> >( N, vector<double>(N, 1 ) );
		pi = vector<double>(N, 1);

		double npi = 0;
		for(unsigned i=0; i<N; ++i){
			double norm = 0;
			for(unsigned j=0; j<N; ++j){ 
				A[i][j] += -1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/2));
				norm += A[i][j]; 
			} 
			for(int j=0; j<N; ++j){ A[i][j] /= norm; }
			
			pi[i] += -1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/2));
			npi += pi[i]; 
		}
		for(unsigned i=0; i<N; ++i){ pi[i] /= npi; }
		initB();
		
		lhood = -numeric_limits<double>::infinity();
		//don't zero maxA, maxB, ... , to retain between multiple restarts
	

		calc_logP();
	}
		
	//Default
	HMM(){
		setsize(0,0);
		setIters(0,0);
	}
	
	//Constructor
	HMM(int N_, int M_, int min_, int max_){
		setsize(N_,M_);
		setIters(min_,max_);
	}
	
	//re-order states according to no_obs.second
	void sortparams( vector< pair<double, int> > &no_obs, 
	vector<vector<double> > &fA, vector<vector<double> > &fB, vector<double> &fpi ){
		
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
	
	//re-order states so B[0][0] > B[1][0] > B[2][0] > ...
	void sortparams(){ 
		
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

	
	//reserve space for viterbi
	inline void setup_delta(int T_){	
		T = T_;
		delta = vector<vector<double> >(N, vector<double>(T,0)); 
		delta2 = vector<vector<double> >(N, vector<double>(T,1));
		
		state = vector<unsigned>(T,0);
	}
	
	//reserve space for HMM forward-backwards algorithm
	inline void setup(int T_){	
		
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
	
	//reserve space for HMM forward-backwards algorithm
	inline void setup(vector<obs_type> &O){	
		setup( O.size() );
	}
	
	//probability of emitting O in state i
	virtual double pB(int i, obs_type O) = 0;
	//log probability of emitting O in state i
	virtual double log_pB(int i, obs_type O) = 0;	
	
	
	///////////////////
	//	FB algorthm  //
	///////////////////
	//Initial probability
	inline void compute_alpha0(obs_type O){
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
	//Updates forward probabilities. assumes alpha is big enough! 
	inline void update_alpha(obs_type O, int t){
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
	//fill alpha from observation vector
	inline void compute_alpha(vector<obs_type> &O){	
		compute_alpha0(O[0]); 
		for(int t=1; t<T; ++t){ update_alpha(O[t], t);	} //exit(1);
	}
	//fill beta (backward probabilities) from observation vector	
	inline void compute_beta(vector<obs_type> &O){
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
	//fill gamma and digamma (state probabilities) from observation vector	
	inline void compute_gamma(vector<obs_type> &O){
		
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
	
	//compute gamma norm
	inline void compute_sumgamma(){ 
		for(int i=0; i<N; ++i){
			sumgamma[i] = 0;
			for(int t=0; t<T; ++t){
				sumgamma[i] += gamma[i][t];
			}
		}
	}
	
	//////////////////////////////////
	//	FB algorithm (Log probs)	//
	//////////////////////////////////
	inline void compute_log_alpha0(obs_type O){
		for(int i=0; i<N; ++i){
			alpha[i][0] = logpi[i] + log_pB(i,O); 
		}
	}
	
	inline void update_log_alpha(obs_type O, int t){
		double sum;
		for(int i=0; i<N; ++i){
			
			alpha[i][t] = alpha[0][t-1] + logA[0][i]; 
			for(int j=1; j<N; ++j){
				lsum(alpha[i][t], alpha[j][t-1] + logA[j][i]);
			}
			alpha[i][t] += log_pB(i,O); 
		}
	}

	inline void compute_log_alpha(vector<obs_type> &O){	
		compute_log_alpha0(O[0]);  		
		for(int t=1; t<T; ++t){ update_log_alpha(O[t], t);	} 
	}
	

	inline void compute_log_beta(vector<obs_type> &O){
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
		
	inline void compute_log_gamma(vector<obs_type> &O){
		
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
	
	inline void compute_log_sumgamma(){ 
		for(int i=0; i<N; ++i){
			for(int t=0; t<T; ++t){
				if( t==0 ){ sumgamma[i] = gamma[i][t]; } else { lsum(sumgamma[i] , gamma[i][t]); }
				if( t == T-2 ){ sumgammap[i] = sumgamma[i]; }
			}
		}
	}
	
	//////////////////////
	//		re-estimate //
	//////////////////////
		
	//re-estimate pi from model
	inline void reestimate_pi(){
		for(int i=0; i<N; ++i){ 
			if( !fixpirow[i] ){
				pi[i] = gamma[i][0]; 
			}
		}
	}
	inline void reestimate_log_pi(){
		for(int i=0; i<N; ++i){
			if( !fixpirow[i] ){
				pi[i] = exp(gamma[i][0]); 
				logpi[i] = gamma[i][0]; 
			}
		}
	}
	
	//re-estimate A from model
	inline void reestimate_A(){
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
	inline void reestimate_log_A(){
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
	
	//re-estimate B from model
	virtual void reestimate_B(vector<obs_type> &O) = 0;
	virtual void reestimate_log_B(vector<obs_type> &O) = 0;
	
	//Re-estimate HMM
	inline void reestimate(vector<obs_type> &O){ 
		if( !fixpi ){ reestimate_pi(); }
		if( !fixA ){ reestimate_A(); }
		if( !fixB ){ reestimate_B(O); }
	}
	inline void reestimate_log(vector<obs_type> &O){ 
		if( !fixpi ){ reestimate_log_pi(); }
		if( !fixA ){ reestimate_log_A(); }
		if( !fixB ){ reestimate_log_B(O); }
	}
	
	//calc model likelihood
	void evaluate(){
		lhood = 0;
		for(int t=0; t<T; ++t){ lhood -= log( c[t] ); }
	}
	void evaluate_log(){
		lhood = 0;
		for(int i=0; i<N; ++i){ 
			if(i==0){ lhood = alpha[i][T-1]; } else { lsum(lhood, alpha[i][T-1]); }; 
		}
	}
	
	//calc model likelihood for observation seq O
	void evaluate(vector<obs_type> &O){
		
		T = O.size();
		alpha = vector<vector<double> >(N, vector<double>(T,0) ); 	
		c = vector<double>(T,0); 
		compute_alpha(O);
		evaluate();
		
	}
	
	//most likely state sequence in sense of dynamic programing
	/* Possible sequences
	 * aaa
	 * aab
	 * aba
	 * abb
	 * baa
	 * bab
	 * ...
	 * 
	 * Choose one with highest prob
	 * */
	 //Update dynamic programming matrix
	 inline void compute_delta0(obs_type O){
		for(int i=0; i<N; ++i){ 
			delta[i][0] = logpi[i] + log_pB(i,O);
			delta2[i][0] = 0;
		}	
	 }
	 inline void update_delta(obs_type O, int t){
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
	
	inline void backtrack_viterbi(int last){
		double best = delta[0][last-1];
		state[ last-1 ] = 0;
		for(int i=1; i<N; ++i){ 
				if( delta[i][last-1] > best ){ best = delta[i][last-1]; state[ last-1 ] = i; }
		}
		for(int t=last-1; t>0; --t){
			state[t-1] = delta2[state[t]][ t]; 
		}
	}
	
	void set_state_viterbi(vector<obs_type> &O){
		
		setup_delta(O.size());
		compute_delta0(O[0]);	
		for(int t=1; t<O.size(); ++t){ 
			update_delta(O[t], t);
		}
		backtrack_viterbi(O.size());
		
	}

	
	//most likely state sequence in E-M sense
	/* Possible sequences
	 * aaa
	 * aab
	 * aba
	 * abb
	 * baa
	 * bab
	 * ...
	 * gives
	 * P(a,0) P(a,1) P(a,2) ...
	 * P(b,0) P(a,1) P(a,2) ...
	 * Choose most likely state at each time
	 * */
	void set_state(){
		for(int t=0; t<T; ++t){ 
			state[t] = 0;
			double best = alpha[0][t];
			for(int i=1; i<N; ++i){
				if(alpha[i][t] > best){ best = alpha[i][t]; state[t] = i; }			
			}
		}
	}
	//Most likely state given the model to explain O.
	void set_state(vector<obs_type> &O){
		setup(O);
		compute_alpha(O);
		set_state();
	}
	
	
	void printA(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){
				cout << A[i][j] << " ";
			} cout << "\n";
		}
	}
	void printA(ofstream &ofile){
		for(int i=0; i<N; ++i){
			for(int j=0; j<N; ++j){
				ofile << A[i][j] << " ";
			} ofile << "\n";
		}
	}
	
	void printB(){
		for(int i=0; i<N; ++i){
			for(int j=0; j<B[i].size(); ++j){
				cout << B[i][j] << " ";
			} cout << "\n";
		}
	}
	void printB(ofstream &ofile){
		for(int i=0; i<N; ++i){
			for(int j=0; j<B[i].size(); ++j){
				ofile << B[i][j] << " ";
			} ofile << "\n";
		}
	}

	void printpi(){
		for(int i=0; i<N; ++i){
				cout << pi[i] << " ";
		} cout << "\n";
	}
	void printpi(ofstream &ofile){
		for(int i=0; i<N; ++i){
				ofile << pi[i] << " ";
		} ofile << "\n";
	}

	
	//Fit HMM to obs seq O
	vector<double> fit(vector<obs_type> &O, double eps, bool uselog=true){
		
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
			} else if( 
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
	
	///////////////////////////////
	//  Use HMM to generate seq  //
	///////////////////////////////
	//Initial state
	discrete_distribution<int> setup_distpi(){
		discrete_distribution<int> pi_dist(pi.begin(), pi.end());
		return pi_dist;
	}
					
	//Transitions
	vector< discrete_distribution<int> > setup_distA(){
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
	
	//generate observations
	virtual obs_type gen_obs(int state) = 0;
	
	//Generate observation sequence from HMM
	void generate_seq(vector<obs_type> &O, int T, bool print = false){	

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
