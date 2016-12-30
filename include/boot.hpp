#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <iomanip>      
#include "hmm.hpp"

using namespace std;

namespace cdHMM {

class boot {
public:
	int Nboot;				/*!< Number of bootstrap samples */
	vector<double> data;	/*!< Resampled data */
	
	double av;				/*!< Sample Average */
	double err;				/*!< Error Width */
	double uerr;			/*!< Upper error */
	double lerr;			/*!< Lower error */
	bool shuffle;			/*!< Uncorrelation addition etc. */
	
	double confidence;		/*!< Confidence interval */
	
	/*! 
	 * Empty boot */
	boot(){
		Nboot = 100;
		data.resize(Nboot);
		shuffle = false;
		confidence = 0.6827;
	}
	
	/*! 
	 * Empty boot 
	 * @param[in] Nboot_ number of bootstrap samples*/
	boot(int Nboot_){
		Nboot = Nboot_;
		data.resize(Nboot);
		shuffle = false;
		confidence = 0.6827;
	}
	/*! 
	 * Empty boot 
	 * @param[in] data bootstrap samples*/	
	boot(vector<double> &data_){
		Nboot = data_.size();
		data = data_;
		shuffle = false;
		confidence = 0.6827;
	}
	/*! 
	 * Swap boot 
	 * @param[out] first 
	 * @param[out] second 
	 * */	
	friend void swap(boot& first, boot& second)
    {
        swap(first.Nboot, second.Nboot);
        swap(first.data, second.data);
    }
	/*! 
	 * Copy boot 
	 * @param[in] other 
	 * */		
    boot& operator=(boot other)
	{
		swap(*this, other); 
		return *this;
	}
	/*! 
	 * Sum equals bootstraps
	 * @param[in] rhs 
	 * */		
	boot& operator+=(const boot &rhs) {
		Nboot = rhs.Nboot;
		data.resize(Nboot);
		if( shuffle ){
			vector<double> shuf = rhs.data;
			random_shuffle ( shuf.begin(), shuf.end() );
			for(int i=0; i<Nboot; ++i){
				data[i] += shuf[i];
			}
		} else {
			for(int i=0; i<Nboot; ++i){
				data[i] += rhs.data[i];
			}
		}
		
		return *this;
    }
	/*! 
	 * Sum bootstraps
	 * @param[in] other 
	 * */	
	const boot operator+(const boot &other) const {
		return boot(*this) += other;
    }
    /*! 
	 * Times equals bootstraps
	 * @param[in] rhs 
	 * */	
    boot& operator*=(const boot &rhs) {
		Nboot = rhs.Nboot;
		data.resize(Nboot);
		if( shuffle ){
			vector<double> shuf = rhs.data;
			random_shuffle ( shuf.begin(), shuf.end() );
			for(int i=0; i<Nboot; ++i){
				data[i] *= shuf[i];
			}
		} else {
			for(int i=0; i<Nboot; ++i){
				data[i] *= rhs.data[i];
			}
		}
		return *this;
    }
    /*! 
	 * Multiply bootstraps
	 * @param[in] other 
	 * */	
	const boot operator*(const boot &other) const {
		return boot(*this) *= other;
    }
   	/*! 
	 * Calculate boot statistics
	 * @param[in] av_ sample average
	 * */	
	void calc(double av_){
		
		//sort average
		vector<double> sorted_data = data;
		sort( sorted_data.begin(), sorted_data.end() );

		//save average
		av = av_;
		
		//Find data before av and after av
		int lt=0, gt=0;
		for(int b=0; b<sorted_data.size(); ++b){ 
			if( sorted_data[b] < av ){ ++lt; }
			else{ ++gt; } 
		}
		double toss = 0.5*(1.0 - confidence);
		double lb = sorted_data[ (int)(toss*2*lt) ];
		double ub = sorted_data[ (int)(sorted_data.size() - 1 - toss*2*gt) ];

		lerr = (av-lb);
		uerr = (ub-av);
				
	}
	
	void pstream(stringstream& os) const {
		os << av << " +/- {" << uerr << "}/{" << lerr << "}";
	}

	friend ostream& operator<<(std::ostream& os, const boot& b){
		std::stringstream ss;
		b.pstream(ss);
		os << ss.str();
	  	return os;
	};
	
	
};

/*! 
 * Resample a HMM
 * @param[in] hmm_eval The model 
 * @param[in] bA 	   Resampled transmission matrix
 * @param[in] bB 	   Resampled emission matrix
 * @param[in] bpi 	   Resampled starting vector
 * @param[in] T 	   Number of observations
 * @param[in] Nboot    Number of bootstrap samples
 * @param[in] num_restarts    	Number of hmm restarts 
 * @param[in] eps    			hmm precision target
 * */	
template<typename number> void generate_bootstrap( 
	HMM<number>* hmm, vector<vector< boot> > &bA, vector<vector< boot> > &bB , 
	int T, int Nboot,
	int num_restarts, double eps
	){

	int N = hmm->N;
	int M = hmm->M;
		
	boot empty(Nboot);
	bA.resize(N); for(int i=0; i<N; ++i){ bA[i].resize(N); for(int j=0; j<N; ++j){ bA[i][j] = empty; } }
	bB.resize(N); for(int i=0; i<N; ++i){ bB[i].resize(M); for(int j=0; j<M; ++j){ bB[i][j] = empty; } }

	
	hmm->sortparams();

	vector<vector<number> > At = hmm->getmaxA();
	vector<vector<number> > Bt = hmm->getmaxB();
	vector<number> pit = hmm->getmaxpi();
    double maxlhood = hmm->maxlhood;
  
    //cerr << "Generating bootstraps" << endl;
	for(int b=0; b<Nboot; ++b){

		hmm->setA( At ); hmm->setB( Bt ); hmm->setpi( pit );
		
		vector<number> newO;
		hmm->generate_seq(newO, T, false);
		hmm->clear_max();
		
		for(int r=0;r<num_restarts;++r){
			hmm->init(); 	
			hmm->setA( At ); hmm->setB( Bt ); hmm->setpi( pit );
			vector<double> res = hmm->fit(newO, eps);
			hmm->sortparams();
		}	

		vector<vector<number> > Atmp = hmm->getmaxA();
		vector<vector<number> > Btmp = hmm->getmaxB();
	
		for(int i=0; i<N; ++i){ for(int j=0; j<N; ++j){ bA[i][j].data[b] = Atmp[i][j]; } }
		for(int i=0; i<N; ++i){ for(int j=0; j<M; ++j){ bB[i][j].data[b] = Btmp[i][j]; } }

	}

	for(int i=0; i<N; ++i){ for(int j=0; j<N; ++j){ bA[i][j].calc( At[i][j] ); } } 
    for(int i=0; i<N; ++i){ for(int j=0; j<M; ++j){ bB[i][j].calc( Bt[i][j] ); } } 
	hmm->maxA = At;
	hmm->maxB = Bt;
	hmm->maxpi = pit;
	hmm->maxlhood = maxlhood;
}

}
