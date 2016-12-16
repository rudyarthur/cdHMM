#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

//HMM with poisson emission
//prob( Emit O | state i ) = exp( -lambda_i ) O^\lambda_i / O!
//alpha_i = B[i][0]
template <typename obs_type> class poissonHMM : public HMM<obs_type> {
public:

	//Default
	poissonHMM(){
		this->setsize(0,1);
		this->setIters(0,0);
		
		this->type = "Poisson";
	}
	
	//Constructor
	poissonHMM(int N_, int min_, int max_){
		this->setsize(N_,1);
		this->setIters(min_,max_);
		
		this->type = "Poisson";
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;
		cerr << "prob( Emit O | state=i ) = exp( -lambda_i ) O^lambda_i / O!\n";
		cerr << "lambda_i = B[i][0]" << endl;
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
	}
	
	
	inline double pB(int i, obs_type O){
		return (O == 0) ? exp( -this->B[i][0] ) : exp( -this->B[i][0] + O * this->logB[i][0] - lgamma(O) );	//why is there no std::factorial ?!
	}
	inline double log_pB(int i, obs_type O){
		return (O == 0) ? ( -this->B[i][0] ) : ( -this->B[i][0] + O * this->logB[i][0] - lgamma(O) );	
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			this->B[i][0] = 0; 
			for( int t=0; t<this->T; ++t){
				this->B[i][0] += this->gamma[i][t] * O[t];
			}
			this->B[i][0] /= this->sumgamma[i];
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
				
				this->B[i][0] = this->gamma[i][0] + this->logO[0];
				for( int t=1; t<this->T; ++t){
					lsum( this->B[i][0] , this->gamma[i][t] + this->logO[t] ); //E[O]
				} 
				this->B[i][0] -= this->sumgamma[i]; //log[ E[O] ] 
				this->B[i][0] = exp( this->B[i][0] ); //E[O]
			}}
			this->calc_logB();
		}
	}	
	
	obs_type gen_obs(int state){
		poisson_distribution<int> b_dist(this->B[state][0]);
		return b_dist(this->generator);
	}
	
	
	
};

