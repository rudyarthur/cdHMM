#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

//HMM with exponential emission
//prob( Emit O | state i ) = a xm^a / x^(a+1) x > xm else 0
//xm_i = lb[i]
//a_i = B[i][0]
template <typename obs_type> class paretoHMM : public HMM<obs_type> {
public:

vector<double> lb;

	//Default
	paretoHMM(){
		this->setsize(0,1);
		this->setIters(0,0);
		
		this->type = "Pareto";
	}
	
	//Constructor
	paretoHMM(int N_, int min_, int max_){
		this->setsize(N_,1);
		this->setIters(min_,max_);
		lb = vector<double>(N_,0);
		
		this->type = "Pareto";
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		for(int i=0; i<this->N; ++i){
			this->B[i][0] = ( static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)) );
		}
	}
	
	
	inline double pB(int i, obs_type O){
		return (O < lb[i]) ? 0 : this->B[i][0] * pow( lb[i], this->B[i][0] ) / pow( O, this->B[i][0] + 1 );	
	}
	inline double log_pB(int i, obs_type O){
		return (O < lb[i]) ? -numeric_limits<double>::infinity() : log( this->B[i][0] ) + this->B[i][0] * log( lb[i] ) - (this->B[i][0] + 1)*log( O );	
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		if( !this->set_logO){ this->calc_logO(O); }
		if( this->no_logO ){ 
			cerr << "Pareto requires O > 0" << endl; exit(1); 
		}
		
		for(int i=0; i<this->N; ++i){ if(!this->fixBrow[i]){

			double av = 0;
			for(int t=0; t<O.size(); ++t){
				av += (log(O[t]) - log(lb[i])) * this->gamma[i][t];
			}
			this->B[i][0] = this->sumgamma[i]/av;

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
		exponential_distribution<double> b_dist(this->B[state][0]); return lb[state] * exp( b_dist(this->generator) );
		//exponential_distribution<double> b_dist(this->B[state][1]); return this->B[state][0] * exp( b_dist(this->generator) );
	}
	
	
	
	
};

