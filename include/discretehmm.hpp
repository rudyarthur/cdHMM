#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

namespace cdHMM {

//HMM with discrete emission
//prob( Emit O | state i ) = B[i][O]
template <typename obs_type> class discreteHMM : public HMM<obs_type> {
public:

	//Default
	discreteHMM(){
		this->setsize(0,0);
		this->setIters(0,0);
		
		this->type = "Discrete";
	}
	
	//Constructor
	discreteHMM(int N_, int M_, int min_, int max_){
		this->setsize(N_,M_);
		this->setIters(min_,max_);
		
		this->type = "Discrete";
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;
		cerr << "prob( Emit O | state=i ) = b_i\n";
		cerr << "b_i = B[i][0]" << endl;
	}
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		for(unsigned i=0; i<this->N; ++i){
			double norm = 0;
			for(unsigned j=0; j<this->M; ++j){ 
				this->B[i][j] += -1 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/2));
				norm += this->B[i][j]; 
			} 
			for(int j=0; j<this->M; ++j){ this->B[i][j] /= norm; }
		}
	}
	
	inline double pB(int i, obs_type O){
		return this->B[i][O];
	}
	inline double log_pB(int i, obs_type O){
		return this->logB[i][O];
	}
	
	//re-estimate B from model
	inline void reestimate_B(vector<obs_type> &O){ 
		
		for( int j=0; j<this->M; ++j ){
			for(int i=0; i<this->N; ++i){ this->B[i][j] = 0; }
			for( int t=0; t<this->T; ++t){
				if( O[t] == j ){ 
					for(int i=0; i<this->N; ++i){
						this->B[i][j] += this->gamma[i][t];
					}
				}
			}
			for(int i=0; i<this->N; ++i){ this->B[i][j] /= this->sumgamma[i]; }
		}
	}
	inline void reestimate_log_B(vector<obs_type> &O){ 
		for(int i=0; i<this->N; ++i){
			for( int t=0; t<this->T; ++t){
				this->gamma[i][t] = exp(  this->gamma[i][t] );
			}
			this->sumgamma[i] = exp( this->sumgamma[i] );
		}
		this->reestimate_B(O);
		this->calc_logB();
	}
	
	obs_type gen_obs(int state){
		discrete_distribution<int> b_dist(this->B[state].begin(), this->B[state].end());
		return b_dist(this->generator);
	}
	
};

}
