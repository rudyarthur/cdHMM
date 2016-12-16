#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;


//HMM with two emission types
template <typename obs_type> class jointHMM : public HMM<obs_type> {
public:

HMM<obs_type> *hmm1;
HMM<obs_type> *hmm2;
bool setup;

	//Default
	jointHMM(){
		this->setsize(0,0);
		this->setIters(0,0);
	}
	
	//Constructor
	jointHMM(HMM<obs_type> *hmm1_, HMM<obs_type> *hmm2_, int min_, int max_){
		hmm1 = hmm1_;
		hmm2 = hmm2_;
		this->N = hmm1->N + hmm2->N;
		this->M = (hmm1->M > hmm2->M) ? hmm1->M : hmm2->M;
		this->setsize(this->N,this->M);
		this->setIters(min_,max_);
		setup = false;
	}
	
	void info(){ cerr << "hi" << endl; }
	
	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 1) );
		hmm1->init();
		hmm2->init();
	}
	
	void setB(vector<vector<double> > &inB){
		this->B = inB;
		for(int i=0; i<hmm1->N; ++i){
			for(int j=0; j<hmm1->M; ++j){
				hmm1->B[i][j] = inB[i][j];
			}
		}
		for(int i=0; i<hmm2->N; ++i){ 
			for(int j=0; j<hmm2->M; ++j){
				hmm2->B[i][j] = inB[i+hmm1->N][j];
			}
		}
	}
	
	inline double pB(int i, obs_type O){
		return (i < hmm1->N) ? hmm1->pB(i,O) : hmm2->pB(i - hmm1->N,O);
	}
	inline double log_pB(int i, obs_type O){
		return (i < hmm1->N) ? hmm1->log_pB(i,O) : hmm2->log_pB(i - hmm1->N,O);
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		if( !setup ){
			hmm1->setup(O);
			hmm2->setup(O);
			setup = true;
		}
		for(int i=0; i<hmm1->N; ++i){
			hmm1->gamma[i] = this->gamma[i];
			hmm1->sumgamma[i] = this->sumgamma[i];
		}
		hmm1->reestimate_B(O);
		for(int i=hmm1->N; i<this->N; ++i){
			hmm2->gamma[i-hmm1->N] = this->gamma[i];
			hmm2->sumgamma[i-hmm1->N] = this->sumgamma[i];
		}
		hmm2->reestimate_B(O);

		for(int i=0; i<hmm1->N; ++i){ if(!this->fixBrow[i]){
			for(int j=0; j<hmm1->M; ++j){
				this->B[i][j] = hmm1->B[i][j];
			}
		}}
		for(int i=0; i<hmm2->N; ++i){ if(!this->fixBrow[i+hmm1->N]){
			for(int j=0; j<hmm2->M; ++j){
				this->B[i+hmm1->N][j] = hmm2->B[i][j];
			}
		}}
		
	}
	void reestimate_log_B(vector<obs_type> &O){ 

		if( !setup ){
			hmm1->setup(O);
			hmm2->setup(O);
			setup = true;
		}
		for(int i=0; i<hmm1->N; ++i){
			hmm1->gamma[i] = this->gamma[i];
			hmm1->sumgamma[i] = this->sumgamma[i];
		}
		hmm1->reestimate_log_B(O);
		for(int i=hmm1->N; i<this->N; ++i){
			hmm2->gamma[i-hmm1->N] = this->gamma[i];
			hmm2->sumgamma[i-hmm1->N] = this->sumgamma[i];
		}
		hmm2->reestimate_log_B(O);

		for(int i=0; i<hmm1->N; ++i){ if(!this->fixBrow[i]){
			for(int j=0; j<hmm1->M; ++j){
				this->B[i][j] = hmm1->B[i][j];
			}
		}}
		for(int i=0; i<hmm2->N; ++i){ if(!this->fixBrow[i+hmm1->N]){
			for(int j=0; j<hmm2->M; ++j){
				this->B[i+hmm1->N][j] = hmm2->B[i][j];
			}
		}}
		
	}	
	
	obs_type gen_obs(int state){
		return ( state < hmm1->N ) ? hmm1->gen_obs(state) : hmm2->gen_obs(state-hmm1->N);		
	}
	
};

