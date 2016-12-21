#pragma once

#include <math.h>       
#include "optimise/function.hpp"
#include "optimise/simplex.hpp"

using namespace std;

namespace cdHMM {

template <typename S> class GEV : public FuncNd<S>{
public:
	
	GEV(){ 
		this->dim = 3;
	}
	//From GNU Octave
	double eval(vector<double> &v){
      
		  double mu = v[0]; 
		  double sigma = v[1]; 
		  double k = v[2]; 
		  
		  int T = this->wp.O->size();
		  int i = this->wp.i;
		  
		  vector<double> a(T); for(int t=0; t<T; ++t){ a[t] = ( (double)(*this->wp.O)[t] - mu)/sigma; } //a = (x-m)/s
		  vector<double> k_terms;
		  
		  double nlogL = 0;
		  if(k == 0){ //gumbel case
			for(int t=0; t<T; ++t){ 
				nlogL += (exp(-a[t]) + a[t] + log(sigma))* (*this->wp.gamma)[i][t]; 
			}
		  } else {
			vector<double> aa(T); //aa = k(x-m)/s
			double min_aa = numeric_limits<double>::infinity();
			double max_aa = -numeric_limits<double>::infinity();
			for(int t=0; t<T; ++t){ 
				aa[t] = k * a[t]; 
				if( fabs(aa[t]) < min_aa ){ min_aa = fabs(aa[t]); }
				if( fabs(aa[t]) > max_aa ){ max_aa = fabs(aa[t]); }
				
			}
			if( min_aa < 1e-3 && max_aa < 0.5 ){

			  k_terms.resize(T);  for(int t=0; t<T; ++t){ k_terms[t] = 1; }
			  double sgn = 1; 
			  double ti = 0;
			  while(true){ 
				//this computes \sum_i log(1 +  k(x_i-m)/s ) for small Y = k(x_i-m)/s
				//log(1 + Y) = Y - Y^2/2 + Y^3/3 - ...
				sgn = -sgn; ti++;
				double max_nt = -numeric_limits<double>::infinity();
				for(int t=0; t<T; ++t){ 
					double newterm = (sgn  / (ti + 1)) * pow(aa[t],ti); 
					if( fabs(newterm) > max_nt ){ max_nt = fabs(newterm); }
					k_terms[t] += newterm;
				}
				if( max_nt <= 1e-15 ){
				  break;
				}
			  }
			  //a[t] (k+1) k_terms[t] = k(x-m)/s
			  for(int t=0; t<T; ++t){ 
				  nlogL += (exp( - a[t] * k_terms[t] ) + a[t] * (k+1) * k_terms[t] + log(sigma))* (*this->wp.gamma)[i][t]; 
			  }
			} else {
			  vector<double> b(T); for(int t=0; t<T; ++t){ b[t] = 1 + aa[t]; if(b[t] <= 0 ){ return 1e200; } } //b = 1 + k(x-m)/s
			  for(int t=0; t<T; ++t){ 
				  nlogL += (pow( b[t], -1/k ) + (1+1/k) * log(b[t]) + log(sigma))* (*this->wp.gamma)[i][t]; 
			  }
			}
		  }
		  return nlogL/(*this->wp.sumgamma)[i];
	}

};

template<typename T> fit_params extreme_solve(vector<T> &O, vector< vector<double> > &gamma, vector<double> &sumgamma, int i,
max_lhood_params mlp){
	
	GEV<T> my_f;
	my_f.wp.O = &O;
	my_f.wp.gamma = &gamma;
	my_f.wp.sumgamma = &sumgamma;
	my_f.wp.i = i;
	

	return simplex(mlp.x_start, mlp.max_iter, mlp.eps, my_f);

}

}
