#pragma once

#include <math.h>       
#include "optimise/function.hpp"
#include "optimise/brent.hpp"

using namespace std;

class Weibull : public Func1d{
public:
	
	Weibull(){ 
		this->dim = 1;
	}
	
	double eval(double k){
			
		  double ElnO = 0;
		  double EOk = 0;
		  int T = wp.O->size();
		  int i = wp.i;
		  int t = 0;
		  
		  vector<double> vlOk(T);
		  vector<double> vlnO(T);

		  for(int t=0; t<T; ++t){ 
			  
			  double lO = log( (*wp.O) [t] );
			  
			  ElnO += ( (*wp.gamma)[i][t] * lO ); //E[ lnO ]

			  double lg = log( (*wp.gamma)[i][t] );
			  vlOk[t] = lg + k*lO;	//log( gamma_t O^k )
			  if(t == 0){ EOk = vlOk[t]; } else { lsum(EOk, vlOk[t]); } //log( E(O^k) ) 
			  vlnO[t] = lO;			//log(O);
			  
		  }
		  ElnO /= ( (*wp.sumgamma)[i] ); //E[ lnO ]
		  
		  pair<int, double> EOklnO = safe_log_sum(vlnO, vlOk);
		  //1/k = -E[ ln O ] + E[ O^k ln O ]/E[ O^k ]
		  return (1.0/k) + ElnO - exp( (EOklnO.first * EOklnO.second) - EOk );
	}

};

double weibull_solve(vector<double> &O, vector< vector<double> > &gamma, vector<double> &sumgamma, int i,
double start, double end, int max_iter, double eps){
	
	Weibull my_f;
	my_f.wp.O = &O;
	my_f.wp.gamma = &gamma;
	my_f.wp.sumgamma = &sumgamma;
	my_f.wp.i = i;
	
	vector<double> range = {start, end};
	return brent(range, max_iter, eps, my_f);
	return 0;
	
}
