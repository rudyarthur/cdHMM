#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

//HMM with laplace emission
//prob( Emit O | state i ) = (1 / 2b ) exp( -|x - mu|/b )
//mu_i = B[i][0]
//b_i = B[i][1]
template <typename obs_type> class laplaceHMM : public HMM<obs_type> {
public:

	//Default
	laplaceHMM(){
		this->setsize(0,2);
		this->setIters(0,0);
		
		this->type = "Laplace";
	}
	
	//Constructor
	laplaceHMM(int N_, int min_, int max_){
		this->setsize(N_,2);
		this->setIters(min_,max_);
		
		this->type = "Laplace";
	}
		
	void info(){
		cerr << this->type << " HMM" << endl;
		cerr << "prob( Emit O | state=i ) = (1 / 2b_i ) exp( -|O - mu_i|/b_i )\n";
		cerr << "mu_i = B[i][0]\n";
		cerr << "b_i = B[i][1]" << endl;
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
	}

	
	inline double pB(int i, obs_type O){
		return exp( - fabs(O - this->B[i][0]) / this->B[i][1] ) / (2*this->B[i][1]);	
	}
	inline double log_pB(int i, obs_type O){
		return ( - fabs(O - this->B[i][0]) / this->B[i][1] ) - log(2*this->B[i][1]);	
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 
		
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			
			//These aren't Max Likelihood Estimates!
			this->B[i][0] = 0;
			for( int t=0; t<this->T; ++t){
				this->B[i][0] += this->gamma[i][t] * O[t];
			}
			this->B[i][0] /= this->sumgamma[i];
			
			this->B[i][1] = 0;
			for( int t=0; t<this->T; ++t){
				this->B[i][1] += this->gamma[i][t] * ( O[t] - this->B[i][0] ) * ( O[t] - this->B[i][0] );
			}
			this->B[i][1] /= this->sumgamma[i];
			this->B[i][1] = sqrt( 0.5*this->B[i][1] );
			 
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
		exponential_distribution<double> b_dist(1/this->B[state][1]); 
		return   (rand()%2) ? this->B[state][0] + b_dist(this->generator) : this->B[state][0] - b_dist(this->generator);
	}
	
	
	
	
};

