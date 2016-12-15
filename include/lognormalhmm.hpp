#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

//HMM with log normal emission
//prob( Emit O | state i ) = ( 1 / O sigma sqrt(2 pi) )  exp( -(ln(O) - mu)*(ln(O) - mu)/2sigma^2 ) 
//mu_i = B[i][0]
//sigma_i = B[i][1]
template <typename obs_type> class lognormalHMM : public HMM<obs_type> {
public:

vector<double> factor;
vector<double> log_factor;

	//Default
	lognormalHMM(){
		this->setsize(0,2);
		this->setIters(0,0);
	}
	
	//Constructor
	lognormalHMM(int N_, int min_, int max_){
		this->setsize(N_,2);
		this->setIters(min_,max_);
		factor.resize(N_);
		log_factor.resize(N_);
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
	
	void initB(vector<obs_type> &O){
		
		initB();
		double av = 0;
		int nm = O.size();
		for(int t=0; t<O.size(); ++t){ av += O[t]; }
		av /= nm; 
		double sig = 0;
		for(int t=0; t<O.size(); ++t){ sig += (O[t] - av)*(O[t] - av); }
		sig = sqrt( sig / nm );

		for(int i=0; i<this->N; ++i){
			this->B[i][0] = av; 
			this->B[i][1] = sig;
		}
		
		this->calc_logB();
		calc_factor();
	}
	
	//Blows up for 0 observations
	inline double pB(int i, obs_type O){
		return factor[i] * (1.0/O) * exp( -(log(O) - this->B[i][0])*(log(O) - this->B[i][0])/( 2*this->B[i][1]*this->B[i][1] ) );			  
	}
	inline double log_pB(int i, obs_type O){
		return log_factor[i] - (log(O) - this->B[i][0])*(log(O) - this->B[i][0])/( 2*this->B[i][1]*this->B[i][1] ) - log(O);			    
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 
		
		if( !this->set_logO){ this->calc_logO(O); }
		if( this->no_logO ){ 
			cerr << "Log Normal function requires O > 0" << endl; exit(1); 
		}
		
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			this->B[i][0] = 0;  
			this->B[i][1] = 0;
			
			for( int t=0; t<this->T; ++t){
				this->B[i][0] += this->gamma[i][t] * log( O[t] ); //E[O]
			}
			this->B[i][0] /= this->sumgamma[i]; //E[O] = mu
			for( int t=0; t<this->T; ++t){
				this->B[i][1] += this->gamma[i][t] * ( log( O[t] ) - this->B[i][0] )*( log( O[t] ) - this->B[i][0] ); //E[ (O-mu)^2 ]
			}
			this->B[i][1] = sqrt(this->B[i][1] / this->sumgamma[i]); //sigma = sqrt( E[ (O-mu)^2 ] )
		}}
		this->calc_logB();
		calc_factor();
	}
	
	void reestimate_log_B(vector<obs_type> &O){ 
		
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			for( int t=0; t<this->T; ++t){
				this->gamma[i][t] = exp(  this->gamma[i][t] );
			}
			this->sumgamma[i] = exp( this->sumgamma[i] );
		}}
		this->reestimate_B(O); //avoid log(log(O))
			
	}	
	
	obs_type gen_obs(int state){
		lognormal_distribution<double> b_dist(this->B[state][0], this->B[state][1]);
		return b_dist(this->generator);
	}
	
	
};
