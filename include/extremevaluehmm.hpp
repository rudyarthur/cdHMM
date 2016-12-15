#pragma once

#include "utils.hpp"
#include "hmm.hpp"
#include "polysolve.hpp"
#include "extremevalue.hpp"

using namespace std;

//HMM with GEV emission
//prob( Emit O | state i ) = (1/sigma) * t^(kappa - 1) exp( - t )
//t = ( 1 - kappa ((O - mu)/sigma) )^1/kappa
//if kappa = 0: t = exp( -(O - mu)/sigma )
//wikipedia notation
//kappa = -xi
//mu_i = B[i][0]
//sigma_i = B[i][1]
//kappa_i = B[i][2]
template <typename obs_type> class extremevalueHMM : public HMM<obs_type> {
public:

double logsmall;
double solver_eps;
int solver_max_iter;

	//Default
	extremevalueHMM(){
		this->setsize(0,3);
		this->setIters(0,0);
		logsmall = -1e6;
		solver_eps = 1e-4;
		solver_max_iter = 10000;
	}
	
	//Constructor
	extremevalueHMM(int N_, int min_, int max_){
		this->setsize(N_,3);
		this->setIters(min_,max_);
		logsmall = -1e6;
		solver_eps = 1e-4;
		solver_max_iter = 10000;
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
			return logsmall; //-numeric_limits<double>::infinity();
		} 

		return -( (1 / this->B[i][2]) + 1) * log(t)  - pow( t, -(1 / this ->B[i][2]) ) - log( this->B[i][1] );
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){
				
			vector<double> in(3);
			//if( this->B[i][2] == 0 ){
				//GEV_moments(in, O, this->gamma[i], this->sumgamma[i]);
			//} else {
				in = this->B[i];
			//}
			//in[2] = 0;
			
			//bool success = GEV_nrsolve(in, this->B[i], O, this->gamma[i], solver_max_iter, solver_eps);
			//if( !success ){
				int num_iters = multi_solve(this->B[i], O, this->gamma, this->sumgamma, i, in, solver_max_iter, solver_eps);
				if( num_iters == solver_max_iter ){
					cerr << "first iterate failed" << endl;
					in = this->B[i]; in[2] = 0;
					num_iters = multi_solve(this->B[i], O, this->gamma, this->sumgamma, i, in, solver_max_iter, solver_eps);
				}
				cerr << num_iters << " " << this->B[i][0] << " " << this->B[i][1] << " " << this->B[i][2] << endl;	
			//}
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

