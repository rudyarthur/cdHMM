#pragma once

#include "utils.hpp"
#include "hmm.hpp"
#include "gumbel.hpp"

using namespace std;

namespace cdHMM {
	
/*! HMM with gumbel emission \n
//prob( Emit O | state i ) = ( 1 / sigma )  exp( -(O - mu)/sigma - exp( - (O - mu)/sigma ) ) \n 
//mu_i = B[i][0] \n
//sigma_i = B[i][1]*/
template <typename obs_type> class gumbelHMM : public HMM<obs_type> {
public:

	//Default
	gumbelHMM(){
		this->setsize(0,2);
		this->setIters(0,0);
		
		this->type = "Gumbel";
	}
	
	//Constructor
	gumbelHMM(int N_, int min_, int max_){
		this->setsize(N_,2);
		this->setIters(min_,max_);

		this->type = "Gumbel";
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;
		cerr << "prob( Emit O | state=i ) = ( 1 / sigma_i )  exp( -(O - mu_i)/sigma_i - exp( - (O - mu_i)/sigma_i ) ) \n";
		cerr << "mu_i = B[i][0]" << endl;
		cerr << "sigma_i = B[i][1]" << endl;
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		for(int i=0; i<this->N; ++i){
			this->B[i][0] = -1 + 2* ( static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)) );
			this->B[i][1] = 1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
		}
	}
	
	inline double pB(int i, obs_type O){
		double t = ((O - this->B[i][0])/this->B[i][1]);
		return exp ( -t - exp(-t) ) / this->B[i][1];	
	}
	inline double log_pB(int i, obs_type O){		
		double t = ((O - this->B[i][0])/this->B[i][1]);
		return ( -t - exp(-t) ) - log( this->B[i][1] );	
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		fit_params fp;
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			this->B[i][0] = 0;  
			this->B[i][1] = 0;
			
			
			fp = gumbel_solve(O, this->gamma, this->sumgamma, i, this->mlp);
			if( fp.iter == this->mlp.max_iter ){
				cerr << "Max Likelihood estimate for gumbel params failed to converge in " << fp.iter << " iterations" << endl;
				cerr << "Residual = " << fp.residual << endl;
			}
			this->B[i][1] = fp.root;
			
			this->B[i][0] = log(this->gamma[i][0]) - (O[0] / this->B[i][1]);
			for( int t=0; t<this->T; ++t){
				lsum( this->B[i][0] , log(this->gamma[i][t]) - (O[t] / this->B[i][1]) );
			}
			this->B[i][0] -= log(this->sumgamma[i]);
			this->B[i][0] *= -this->B[i][1];
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
		extreme_value_distribution<double> b_dist(this->B[state][0], this->B[state][1]);
		return b_dist(this->generator);
	}
	
	
	
	
};

}
