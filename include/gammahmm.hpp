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

double solve_x_lo;
double solve_x_hi;
double solve_prec;

vector<double> factor;
vector<double> log_factor;

	//Default
	gammaHMM(){
		solve_x_lo = 1e-10;
		solve_x_hi = 10000;
		solve_prec  =1e-5;

		this->setsize(0,2);
		this->setIters(0,0);
	}
	
	//Constructor
	gammaHMM(int N_, int min_, int max_){
		solve_x_lo = 1e-10;
		solve_x_hi = 10000;
		solve_prec  =1e-5;
		
		this->setsize(N_,2);
		this->setIters(min_,max_);
		factor.resize(N_);
		log_factor.resize(N_);
	}
	
	void calc_factor(){
		for(int i=0; i<this->N; ++i){

			
			factor[i] = exp( this->B[i][0] * log( this->B[i][1] ) - lgamma( this->B[i][0] ) );
			log_factor[i] = ( this->B[i][0] * log( this->B[i][1] ) - lgamma( this->B[i][0] ) );
			
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
	/*
	void initB(vector<obs_type> &O){
		initB();
		double av = 0;
		double avl = 0;
		int nm = 0;
		for(int t=0; t<O.size(); ++t){
			if(O[t] != 0){ av += O[t]; avl += log( O[t] ); ++nm; }
		}
		av /= nm; 
		double lav = log( av );
		avl /= nm;

		cout << "here " << lav - avl << endl;
		double alpha = gamma_solve( (lav - avl) , 0.1, 100, 1e-5);
		double beta = av / alpha;
		cout << alpha << " " << beta << endl;
		for(int i=0; i<this->N; ++i){
			this->B[i][1] = beta; 
			this->B[i][0] = alpha;
		}
		
		
		calc_factor();
		this->calc_logB();
	}*/
	
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
			double alpha = gamma_solve( exlx , solve_x_lo, solve_x_hi, 100, solve_prec);

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
			double alpha = gamma_solve( exlx , solve_x_lo, solve_x_hi, 100, solve_prec);
			//double alpha = gamma_solve1( exlx , solve_x_lo, 100, solve_prec);

			//cerr << alpha << " " << alpha1 << endl;

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

