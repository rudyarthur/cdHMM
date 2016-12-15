#pragma once

#include "utils.hpp"
#include "hmm.hpp"
#include "gamma.hpp"
#include <random>


using namespace std;


//HMM with gamma emission
//prob( Emit O | state i ) = ( beta_i^alpha_i / Gamma( alpha_i ) )   O^{alpha_i-1}   exp(-beta_i O) 
//alpha_i = B[i][0]
//beta_i = B[i][1]
//requires x > 0 strictly
template <typename obs_type> class gammaHMM : public HMM<obs_type> {
public:

vector<double> factor;
vector<double> log_factor;

	//Default
	gammaHMM(){
		this->mlp.dim = 1;
		this->mlp.x_start = {1e-10, 10000};
		this->mlp.max_iter = 100;
		
		this->setsize(0,2);
		this->setIters(0,0);
		
		this->type = "Gamma";
	}
	
	//Constructor
	gammaHMM(int N_, int min_, int max_){
		this->mlp.dim = 1;
		this->mlp.x_start = {1e-10, 10000};
		this->mlp.max_iter = 100;

		this->setsize(N_,2);
		this->setIters(min_,max_);
		factor.resize(N_);
		log_factor.resize(N_);

		this->type = "Gamma";
	}
	
	void calc_factor(){
		for(int i=0; i<this->N; ++i){

			
			factor[i] = exp( this->B[i][0] * log( this->B[i][1] ) - lgamma( this->B[i][0] ) );
			log_factor[i] = ( this->B[i][0] * log( this->B[i][1] ) - lgamma( this->B[i][0] ) );
			
		}
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;		
		cerr << "prob( Emit O | state=i ) = ( beta_i^alpha_i / Gamma( alpha_i ) )   O^{alpha_i-1}   exp(-beta_i O)\n";
		cerr << "alpha_i = B[i][0]\n";
		cerr << "beta_i = B[i][1]\n";
		cerr << "requires x > 0 strictly" << endl;
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		for(int i=0; i<this->N; ++i){
			this->B[i][0] = 1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
			this->B[i][1] = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
		}
		calc_factor();
	}

	//Blows up for 0 observations
	inline double pB(int i, obs_type O){	   
		return factor[i] * exp( (this->B[i][0] - 1) * log(O) - this->B[i][1] * O );			  
	}
	inline double log_pB(int i, obs_type O){
		return log_factor[i] + (this->B[i][0] - 1) * log(O) - this->B[i][1] * O;
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 
		
		if( !this->set_logO){ this->calc_logO(O); }
		if( this->no_logO ){ cerr << "Gamma function requires O > 0" << endl; exit(1); }
		
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){

			this->B[i][0] = 0;  this->B[i][1] = 0;
			
			for( int t=0; t<this->T; ++t){
				this->B[i][0] += this->gamma[i][t] * O[t]; //E[O]
				this->B[i][1] += this->gamma[i][t] * this->logO[t]; //E[ log(O) ] Use for Max Likelihood Estimate, no 0 obs!
			}
			this->B[i][0] /= this->sumgamma[i]; //E[O] = alpha/beta
			this->B[i][1] /= this->sumgamma[i]; //E[ log(O) ] = ...

			double exlx = ( log(this->B[i][0]) - this->B[i][1]);
			if( exlx != exlx ){
				cerr << "nan probabilities try to fit with uselog = true" << endl; exit(1);
			}
			double alpha = gamma_solve( exlx , this->mlp); 

			double beta = alpha / this->B[i][0] ;

		}}
		calc_factor();
		this->calc_logB();
	}
	
	void reestimate_log_B(vector<obs_type> &O){ 
		
		if( !this->set_logO){ this->calc_logO(O); }
		if( this->no_logO ){ cerr << "Gamma function requires O > 0" << endl; exit(1); }
		
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			
			double old_a = this->B[i][0];
			this->B[i][0] = this->gamma[i][0] + this->logO[0];
			this->B[i][1] = exp( this->gamma[i][0] )*this->logO[0];
			for( int t=1; t<this->T; ++t){
				lsum( this->B[i][0] , this->gamma[i][t] + this->logO[t] ); //E[O]
				this->B[i][1] += exp( this->gamma[i][t] )*this->logO[t]; //E[ log(O) ] avoid log(log())
			} 

			this->B[i][0] -= this->sumgamma[i]; //log[ E[O] ] 
			this->B[i][1]  /= exp( this->sumgamma[i] ); //E[ log[O] ] 			
			
			double exlx = ( this->B[i][0] - this->B[i][1] ); 
			double alpha = gamma_solve( exlx , this->mlp); 
			double beta = alpha / exp( this->B[i][0] ) ;
			
			this->B[i][0] = alpha;
			this->B[i][1] = beta;
			
		}}
		calc_factor();	
		this->calc_logB();
	}	
	
	obs_type gen_obs(int state){
		gamma_distribution<double> b_dist(this->B[state][0], 1.0/this->B[state][1]);
		return b_dist(this->generator);
	}
	
};

