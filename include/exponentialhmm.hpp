#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

//HMM with exponential emission
//prob( Emit O | state i ) = lambda exp( -lambda_i O ) 
//lambda_i = B[i][0]
template <typename obs_type> class exponentialHMM : public HMM<obs_type> {
public:

vector<double> lb;

	//Default
	exponentialHMM(){
		this->setsize(0,1);
		this->setIters(0,0);
	}
	
	//Constructor
	exponentialHMM(int N_, int min_, int max_){
		this->setsize(N_,1);
		this->setIters(min_,max_);
		lb = vector<double>(N_,0);
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
	}
	
	void initB(vector<obs_type> &O){
		this->initB();
		double av = 0;
		int nm = 0;
		for(int t=0; t<O.size(); ++t){
			if(O[t] != 0){ av += O[t]-lb[0]; ++nm; }
		}
		av /= nm;
		for(int i=0; i<this->N; ++i){
			this->B[i][0] = 1.0/av;
		}
		this->calc_logB();
	}

	
	inline double pB(int i, obs_type O){
		return this->B[i][0] * exp( - this->B[i][0] * (O-lb[i]) );	
	}
	inline double log_pB(int i, obs_type O){
		return this->logB[i][0] + ( - this->B[i][0] * (O-lb[i]) );	
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			this->B[i][0] = 0; 
			for( int t=0; t<this->T; ++t){
				this->B[i][0] += this->gamma[i][t] * (O[t]-lb[i]);
			}
			this->B[i][0] = this->sumgamma[i]/this->B[i][0]; 
		}}
		this->calc_logB();

	}
	void reestimate_log_B(vector<obs_type> &O){ 

		if( !this->set_logO){ this->calc_logO(O); }
		if( this->no_logO ){ 
			for(int i=0; i<this->N; ++i){
				for( int t=0; t<this->T; ++t){
					this->gamma[i][t] = exp(  this->gamma[i][t] );
				}
				this->sumgamma[i] = exp( this->sumgamma[i] );
			}
			this->reestimate_B(O);
		} else {
			
			for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
				this->B[i][0] = 0;  
				
				this->B[i][0] = this->gamma[i][0] + log(O[0] - lb[i]); //+ this->logO[0];
				for( int t=1; t<this->T; ++t){
					lsum( this->B[i][0] , this->gamma[i][t] + log(O[t] - lb[i]) ); //this->logO[t] ); //E[O]
				} 

				this->B[i][0] -= this->sumgamma[i]; //log[ E[O] ] 
				this->B[i][0] = exp( -this->B[i][0] ); //1 / E[O]
			}}
			this->calc_logB();
			
		}
	}	

	obs_type gen_obs(int state){
		exponential_distribution<double> b_dist(this->B[state][0]);
		return b_dist(this->generator) + lb[state];
	}
	
	
	
	
};

