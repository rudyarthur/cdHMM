#pragma once

#include <math.h>       
#include "optimise/function.hpp"
#include "optimise/brent.hpp"

using namespace std;

namespace cdHMM {
/*!
 * Weibull probability density function
 * */
template <typename S> class Weibull : public Func1d<S>{
public:
/*!
 * Evaluate Weibull probability density function
 * */
	double eval(double k){
			
		  double ElnO = 0;
		  double EOk = 0;
		  int T = this->wp.O->size();
		  int i = this->wp.i;
		  int t = 0;
		  
		  vector<double> vlOk(T);
		  vector<double> vlnO(T);

		  for(int t=0; t<T; ++t){ 
			  
			  double lO = log( (double)(*this->wp.O) [t] );
			  
			  ElnO += ( (*this->wp.gamma)[i][t] * lO ); //E[ lnO ]

			  double lg = log( (*this->wp.gamma)[i][t] );
			  vlOk[t] = lg + k*lO;	//log( gamma_t O^k )
			  if(t == 0){ EOk = vlOk[t]; } else { lsum(EOk, vlOk[t]); } //log( E(O^k) ) 
			  vlnO[t] = lO;			//log(O);
			  
		  }
		  ElnO /= ( (*this->wp.sumgamma)[i] ); //E[ lnO ]
		  
		  pair<int, double> EOklnO = safe_log_sum(vlnO, vlOk);
		  //1/k = -E[ ln O ] + E[ O^k ln O ]/E[ O^k ]
		  return (1.0/k) + ElnO - exp( (EOklnO.first * EOklnO.second) - EOk );
	}

};

/*!
 * Max likelihood estimate of Weibull distribution parameters given observation and model weights
 * */
template <typename T> fit_params weibull_solve(vector<T> &O, vector< vector<double> > &gamma, vector<double> &sumgamma, int i,
max_lhood_params mlp){
	
	Weibull<T> my_f;
	my_f.wp.O = &O;
	my_f.wp.gamma = &gamma;
	my_f.wp.sumgamma = &sumgamma;
	my_f.wp.i = i;
	
	return brent(mlp.x_start, mlp.max_iter, mlp.eps, my_f);

	
}

}
