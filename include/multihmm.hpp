#pragma once

#include "utils.hpp"
#include "hmm.hpp"

using namespace std;

namespace cdHMM {

//HMM with two emission types
template <typename obs_type> class multiHMM : public HMM<obs_type> {
public:

vector< HMM<obs_type>* > hmm;
bool setup;

	//Default
	multiHMM(){
		this->setsize(0,0);
		this->setIters(0,0);
		
		this->type = "multi";
	}
	
	//Constructor
	multiHMM( vector< HMM<obs_type>* > hmm_, int min_, int max_){
		hmm = hmm_;
		this->N = 0;
		this->M = 0;
		for(int i=0; i<hmm.size(); ++i){
			this->N += hmm[i]->N;
			if(hmm[i]->M > this->M){ this->M = hmm[i]->M; }
		}
		this->setsize(this->N,this->M);
		this->setIters(min_,max_);
		setup = false;
		this->type = "multi";
	}
	
	void info(){
		cerr << this->type << " HMM" << endl;		
		cerr << "prob( Emit O | state=i )\n";
		int ct = 0;
		for(int i=0; i<hmm.size(); ++i){
			cerr << ct << " <= i < " << ct + hmm[i]->N << "\n";
			hmm[i]->info();
			ct += hmm[i]->N;
		}
	}
	

	
	void initB(){
		this->B = vector< vector<double> >( this->N, vector<double>(this->M, 0) );
		for(int i=0; i<hmm.size(); ++i){ hmm[i]->init(); }
	}
	
	void setB(vector<vector<double> > &inB){
		
		this->B = inB;
		double add = 0;
		for(int h=0; h<hmm.size(); ++h){
			for(int i=0; i<hmm[h]->N; ++i){
				for(int j=0; j<hmm[h]->M; ++j){
					hmm[h]->B[i][j] = inB[i+add][j];
				}
			}
			add += hmm[h]->N;
		}
		
	}
	void getB(){
		
		if(this->maxB.size() != this->N){ this->maxB = vector< vector<double> >( this->N, vector<double>(this->M, 0) ); }

		double add = 0;
		for(int h=0; h<hmm.size(); ++h){
			for(int i=0; i<hmm[h]->N; ++i){
				for(int j=0; j<hmm[h]->M; ++j){
					this->B[i+add][j] = hmm[h]->B[i][j];
					this->maxB[i+add][j] = hmm[h]->maxB[i][j];
				}
			}
			add += hmm[h]->N;
		}
		
	}
	void sort_params(){
		for(int i=0; i<hmm.size(); ++i){
			hmm[i]->sort_params();
		}
		getB();
	}
	
	inline double pB(int i, obs_type O){
		int cum = 0;
		int sub = 0;
		for(int h=0; h<hmm.size(); ++h){
			cum += hmm[h]->N;
			if( i < cum ){
				return hmm[h]->pB(i - sub, O);
			}
			sub += hmm[h]->N;
		}
	}
	inline double log_pB(int i, obs_type O){
		int cum = 0;
		int sub = 0;
		for(int h=0; h<hmm.size(); ++h){
			cum += hmm[h]->N;
			if( i < cum ){
				return hmm[h]->log_pB(i - sub, O);
			}
			sub += hmm[h]->N;
		}	
	}
	
	//re-estimate B from model
	void reestimate_B(vector<obs_type> &O){ 

		if( !setup ){
			for(int h=0; h<hmm.size(); ++h){ hmm[h]->setup(O); }
			setup = true;
		}
		int add = 0;
		for(int h=0; h<hmm.size(); ++h){
			for(int i=0; i<hmm[h]->N; ++i){
				hmm[h]->gamma[i] = this->gamma[i+add];
				hmm[h]->sumgamma[i] = this->sumgamma[i+add];
			}
			hmm[h]->reestimate_B(O);
			for(int i=0; i<hmm[h]->N; ++i){ if(!this->fixBrow[i+add]){
				for(int j=0; j<hmm[h]->M; ++j){
					this->B[i+add][j] = hmm[h]->B[i][j];
				}
			}}
			add += hmm[h]->N;
		}
		
	}
	void reestimate_log_B(vector<obs_type> &O){ 

		if( !setup ){
			for(int h=0; h<hmm.size(); ++h){ hmm[h]->setup(O); }
			setup = true;
		}
		int add = 0;
		for(int h=0; h<hmm.size(); ++h){
			for(int i=0; i<hmm[h]->N; ++i){
				hmm[h]->gamma[i] = this->gamma[i+add];
				hmm[h]->sumgamma[i] = this->sumgamma[i+add];
			}
			hmm[h]->reestimate_log_B(O);
			for(int i=0; i<hmm[h]->N; ++i){ if(!this->fixBrow[i+add]){
				for(int j=0; j<hmm[h]->M; ++j){
					this->B[i+add][j] = hmm[h]->B[i][j];
				}
			}}
			add += hmm[h]->N;
		}
		
	}	
	
	obs_type gen_obs(int state){
		int cum = 0;
		int sub = 0;
		for(int h=0; h<hmm.size(); ++h){
			cum += hmm[h]->N;
			if( state < cum ){
				return hmm[h]->gen_obs(state - sub);
			}
			sub += hmm[h]->N;
		}
	}
	
};

}
