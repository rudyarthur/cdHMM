#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

//HMM with gaussian emission
//prob( Emit O | state i ) = ( 1 / sigma sqrt(2 pi) )  exp( -(O - mu)*(O - mu)/2sigma^2 ) 
//mu_i = B[i][0]
//sigma_i = B[i][1]
template <typename obs_type> class gaussianHMM : public HMM<obs_type> {
public:

//Avoid recomputing normalization
vector<double> factor;
vector<double> log_factor;

	//Default
	gaussianHMM(){
		this->setsize(0,2);
		this->setIters(0,0);
		this->type = "Gaussian";
	}
	
	//Constructor
	gaussianHMM(int N_, int min_, int max_){
		this->setsize(N_,2);
		this->setIters(min_,max_);
		factor.resize(N_);
		log_factor.resize(N_);
		this->type = "Gaussian";
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;
		cerr << "prob( Emit O | state=i ) = ( 1 / (sigma_i sqrt(2 pi)) ) * exp( -(O - mu_i)^2/( 2 sigma_i^2) )\n";
		cerr << "mu_i = B[i][0]\n";
		cerr << "sigma_i = B[i][1]" << endl;
	}
	
	void calc_factor(){
		for(int i=0; i<this->N; ++i){
		
			factor[i] = (1.0 / ( this->B[i][1] * sqrt(2*M_PI) ) );
			log_factor[i] = -log( this->B[i][1] * sqrt(2*M_PI) );
			
		}
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		for(int i=0; i<this->N; ++i){
			this->B[i][0] = 1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
			this->B[i][1] = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
		}
		calc_factor();
	}
	
	//Emission probability
	inline double pB(int i, obs_type O){
		return factor[i] * exp( -(O - this->B[i][0])*(O - this->B[i][0])/( 2*this->B[i][1]*this->B[i][1] ) );			  
	}
	inline double log_pB(int i, obs_type O){
		return log_factor[i] - (O - this->B[i][0])*(O - this->B[i][0])/( 2*this->B[i][1]*this->B[i][1] );			    
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 
		
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			this->B[i][0] = 0;  
			this->B[i][1] = 0;
			
			for( int t=0; t<this->T; ++t){
				this->B[i][0] += this->gamma[i][t] * O[t]; //E[O]
			}
			this->B[i][0] /= this->sumgamma[i]; //E[O] = mu
			for( int t=0; t<this->T; ++t){
				this->B[i][1] += this->gamma[i][t] * ( O[t] - this->B[i][0] )*( O[t] - this->B[i][0] ); //E[ (O-mu)^2 ]
			}
			this->B[i][1] = sqrt(this->B[i][1] / this->sumgamma[i]); //sigma = sqrt( E[ (O-mu)^2 ] )
		}}
		this->calc_logB();
		calc_factor();
	}
	
	void reestimate_log_B(vector<obs_type> &O){ 
				
		//O - mu can be negative so can't calculate sum as log without safe sum
		for(int i=0; i<this->N; ++i){
			for( int t=0; t<this->T; ++t){
				this->gamma[i][t] = exp(  this->gamma[i][t] );
			}
			this->sumgamma[i] = exp( this->sumgamma[i] );
		}
		this->reestimate_B(O);

	}	
	
	//Emission probability
	obs_type gen_obs(int state){
		normal_distribution<double> b_dist(this->B[state][0], this->B[state][1]);
		return b_dist(this->generator);
	}
	
	
};

