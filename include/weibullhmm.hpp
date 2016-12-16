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
		
		this->type = "Weibull";
	}
	
	//Constructor
	weibullHMM(int N_, int min_, int max_){
		this->setsize(N_,2);
		this->setIters(min_,max_);
		
		this->type = "Weibull";
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;
		cerr << "prob( Emit O | state=i ) = ( k_i / lambda_i ) ( O / lambda_i )^(k_i-1)  exp( -(O / lambda_i )^k_i )\n";
		cerr << "k_i = B[i][0]" << endl;
		cerr << "lambda_i = B[i][1]" << endl;
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

		if( !this->set_logO){ this->calc_logO(O); }
		if( this->no_logO ){ 
			cerr << "Weibull requires O > 0" << endl; exit(1); 
		}
		fit_params fp;
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			
			fp = weibull_solve(O, this->gamma, this->sumgamma, i, this->mlp);
			if( fp.iter == this->mlp.max_iter ){
				cerr << "Max Likelihood estimate for weibull params failed to converge in " << fp.iter << " iterations" << endl;
				cerr << "Residual = " << fp.residual << endl;
			}
			this->B[i][0] = fp.root;
			
			this->B[i][1] = log( this->gamma[i][0] ) + this->logO[0];
			for(int t=1; t<this->T; ++t){ 
				lsum( this->B[i][1], log( this->gamma[i][t] ) + this->logO[t]  );  //log( E[ O^k ] ) 
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

