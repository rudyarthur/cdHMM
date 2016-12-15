#pragma once

#include "utils.hpp"
#include "hmm.hpp"
#include "weibull.hpp"

using namespace std;

//HMM with weibull emission
//prob( Emit O | state i ) = ( k / lambda ) ( O / lambda )^(k-1)  exp( -(O / lambda )^k ) 
//k_i = B[i][0]
//lambda_i = B[i][1]
template <typename obs_type> class weibullHMM : public HMM<obs_type> {
public:

	//Default
	weibullHMM(){
		this->setsize(0,2);
		this->setIters(0,0);
	}
	
	//Constructor
	weibullHMM(int N_, int min_, int max_){
		this->setsize(N_,2);
		this->setIters(min_,max_);
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		for(int i=0; i<this->N; ++i){
			this->B[i][0] = 0.5 + ( static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)) );
			this->B[i][1] = 1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
		}
	}
	
	inline double pB(int i, obs_type O){
		double lrat = ( log(O) - log( this->B[i][1] ) );
		return (this->B[i][0] / this->B[i][1]) *exp( (this->B[i][0] - 1) * lrat - exp( this->B[i][0] * lrat ) );
	}
	inline double log_pB(int i, obs_type O){		
		double lrat = ( log(O) - log( this->B[i][1] ) );
		return log(this->B[i][0] / this->B[i][1]) + ( (this->B[i][0] - 1) * lrat - exp( this->B[i][0] * lrat ) );
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			
			this->B[i][0] = weibull_solve(O, this->gamma, this->sumgamma, i, 1e-10, 100000, 100, 1e-5);
			
			this->B[i][1] = log( this->gamma[i][0] ) + log( O[0] );
			for(int t=1; t<this->T; ++t){ 
				lsum( this->B[i][1], log( this->gamma[i][t] ) + log( O[t] ) );  //log( E[ O^k ] ) 
			}
			this->B[i][1] -= log(this->sumgamma[i]); //log( E[ O^k ] ) 
			this->B[i][1] = exp( this->B[i][1] / this->B[i][0] ); 

		}}
		this->calc_logB();

	}
	void reestimate_log_B(vector<obs_type> &O){ 

		for(int i=0; i<this->N; ++i){
			for( int t=0; t<this->T; ++t){
				this->gamma[i][t] = exp(  this->gamma[i][t] );
			}
			this->sumgamma[i] = exp( this->sumgamma[i] );
		}
		this->reestimate_B(O);
	}	

	obs_type gen_obs(int state){
		weibull_distribution<double> b_dist(this->B[state][0], this->B[state][1]);
		return b_dist(this->generator);
	}
	
	
	
	
};

