#pragma once

#include "utils.hpp"
#include "hmm.hpp"
#include "extremevalue.hpp"

using namespace std;

namespace cdHMM {

/*! HMM with GEV emission \n
//prob( Emit O | state i ) = (1/sigma) * t^(kappa - 1) exp( - t ) \n
//t = ( 1 - kappa ((O - mu)/sigma) )^1/kappa \n
//if kappa = 0: t = exp( -(O - mu)/sigma ) \n
//wikipedia notation \n
//kappa = -xi \n
//mu_i = B[i][0] \n
//sigma_i = B[i][1] \n
//kappa_i = B[i][2]*/
template <typename obs_type> class extremevalueHMM : public HMM<obs_type> {

public:
	//Default
	extremevalueHMM(){
		this->setsize(0,3);
		this->setIters(0,0);
		this->type = "Extremevalue";

	}
	
	//Constructor
	extremevalueHMM(int N_, int min_, int max_){
		this->setsize(N_,3);
		this->setIters(min_,max_);
		this->type = "Extremevalue";
		
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;
		cerr << "prob( Emit O | state=i ) = (1/sigma_i) * t_i^(kappa_i - 1) exp( - t_i ) \n";
		cerr << "t_i = ( 1 - kappa_i ((O - mu_i)/sigma_i) )^1/kappa_i\n";
		cerr << "if kappa_i = 0: t_i = exp( -(O - mu_i)/sigma_i )\n";
		cerr << "mu_i = B[i][0]" << endl;
		cerr << "sigma_i = B[i][1]" << endl;
		cerr << "kappa_i = B[i][2]" << endl;
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		for(int i=0; i<this->N; ++i){
			this->B[i][0] = 0.5 + ( static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)) );
			this->B[i][1] = 1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
			this->B[i][2] = 0;
		}
	}
	
	inline double pB(int i, obs_type O){
		if(this->B[i][2] == 0){
			double t = ((O - this->B[i][0])/this->B[i][1]);
			return exp ( -t - exp(-t) ) / this->B[i][1];		
		}
		cerr << "not implemented" << endl;
		exit(1);
		return 0;
	}
	inline double log_pB(int i, obs_type O){		
		if(this->B[i][2] == 0){
			double t = ((O - this->B[i][0])/this->B[i][1]);
			return ( -t - exp(-t) ) - log( this->B[i][1] );			
		}
		double t = 1 + this->B[i][2] * ((O - this->B[i][0])/this->B[i][1]);

		if(t < 0){ //no support here, no contrib to lhood
			return this->logsmall; //-numeric_limits<double>::infinity();
		} 

		return -( (1 / this->B[i][2]) + 1) * log(t)  - pow( t, -(1 / this ->B[i][2]) ) - log( this->B[i][1] );
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		fit_params fp;
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
			if( this->sumgamma[i] == 0 ){ continue; }
			
			this->mlp.x_start = this->B[i];

			fp = extreme_solve(O, this->gamma, this->sumgamma, i, this->mlp);
			if( fp.iter == this->mlp.max_iter ){
				this->mlp.x_start = this->B[i]; this->mlp.x_start[2] = 0;
				fp = extreme_solve(O, this->gamma, this->sumgamma, i, this->mlp);
			}
			if( fp.iter == this->mlp.max_iter ){
				cerr << "Max Likelihood estimate for extremevalue params failed to converge in " << fp.iter << " iterations" << endl;
				cerr << "Residual = " << fp.residual << endl;
			}
			this->B[i] = fp.mroot;
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
		
		if(this->B[state][2] == 0) //exactly the gumbel case
		{
			extreme_value_distribution<double> b_dist(this->B[state][0], this->B[state][1]);
			return b_dist(this->generator);    
		}

		uniform_real_distribution<double> b_dist(0.0,1.0);
		double x = b_dist(this->generator);
		
		//(1/kappa) * [  mu kappa + sigma( -1 + (-lnF)^-kappa )  ]
        return (1.0/this->B[state][2]) * ( this->B[state][0] * this->B[state][2] + 
			this->B[state][1] * ( pow( -log(x), -this->B[state][2] ) - 1 ) );
    
	}
	
	
	
	
};

}
