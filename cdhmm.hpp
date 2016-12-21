#pragma once

#include "include/gaussianhmm.hpp" 
#include "include/exponentialhmm.hpp" 
#include "include/gammahmm.hpp" 
#include "include/lognormalhmm.hpp" 
#include "include/paretohmm.hpp"
#include "include/laplacehmm.hpp" 
#include "include/discretehmm.hpp" 
#include "include/poissonhmm.hpp" 
#include "include/gumbelhmm.hpp" 
#include "include/weibullhmm.hpp" 
#include "include/extremevaluehmm.hpp" 
#include "include/multihmm.hpp" 
#include "include/read.hpp" 

#include <string>

namespace cdHMM {

template <typename T> class HMMFactory {
public:
	static HMM<T>* generateHMM(string name, int N, int min_iter, int max_iter, int M=-1){
	
		string lname = name;
		transform(lname.begin(), lname.end(), lname.begin(), ::tolower);    

		if(lname == "gaussian"){
			gaussianHMM<T>* hmm = new gaussianHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "exponential"){
			exponentialHMM<T>* hmm = new exponentialHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "gamma"){
			gammaHMM<T>* hmm = new gammaHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "lognormal"){
			lognormalHMM<T>* hmm = new lognormalHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "pareto"){
			paretoHMM<T>* hmm = new paretoHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "laplace"){
			laplaceHMM<T>* hmm = new laplaceHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "discrete"){
			if( M <= 0 ){ cerr << "Specify number of observation synbols in last argumnet" << endl; exit(1); }
			discreteHMM<T>* hmm = new discreteHMM<T>(N, M, min_iter, max_iter);
			return hmm;
		} else if(lname == "poisson"){
			poissonHMM<T>* hmm = new poissonHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "gumbel"){
			gumbelHMM<T>* hmm = new gumbelHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "weibull"){
			weibullHMM<T>* hmm = new weibullHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else if(lname == "extremevalue"){
			extremevalueHMM<T>* hmm = new extremevalueHMM<T>(N, min_iter, max_iter);
			return hmm;
		} else {
			cerr << "unknown HMM type: " << name << endl;
			exit(1);
		}
	} 

};

}
