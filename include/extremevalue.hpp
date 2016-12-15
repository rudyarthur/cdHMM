#pragma once

#include <math.h>       
#include "optimise/function.hpp"
#include "optimise/brent.hpp"

using namespace std;

class GEV : public FuncNd{
public:
	
	GEV(vector<double> &p){ 
		this->dim = 3;
	}
	
	double eval(vector<double> &v){
      
		  double mu = v[0]; //gsl_vector_get(v, 0);
		  double sigma = v[1]; //gsl_vector_get(v, 1);
		  double k = v[2]; //gsl_vector_get(v, 2);
		  
		  int T = wp.O->size();
		  int i = wp.i;
		  
		  vector<double> a(T); for(int t=0; t<T; ++t){ a[t] = ( (*wp.O)[t] - mu)/sigma; } //a = (x-m)/s
		  vector<double> k_terms;
		  
		  double nlogL = 0;
		  if(k == 0){ //gumbel case
			for(int t=0; t<T; ++t){ 
				nlogL += (exp(-a[t]) + a[t] + log(sigma))* (*wp.gamma)[i][t]; 
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
				  nlogL += (exp( - a[t] * k_terms[t] ) + a[t] * (k+1) * k_terms[t] + log(sigma))* (*wp.gamma)[i][t]; 
			  }
			} else {
			  vector<double> b(T); for(int t=0; t<T; ++t){ b[t] = 1 + aa[t]; if(b[t] <= 0 ){ return 1e200; } } //b = 1 + k(x-m)/s
			  for(int t=0; t<T; ++t){ 
				  //cout << p->O[t] << " " << -(pow( b[t], -1/k ) + (1+1/k) * log(b[t]) + log(sigma)) << endl; exit(1);
				  nlogL += (pow( b[t], -1/k ) + (1+1/k) * log(b[t]) + log(sigma))* (*wp.gamma)[i][t]; 
			  }
			}
		  }
		  return nlogL/(*wp.sumgamma)[i];
	}

};

double gamma_solve1(double pm, double start, double end, int max_iter, double eps){
	
	/*vector<double> p = {pm};
	PsiGamma my_f(p);
	
	vector<double> range = {start, end};
	return brent(range, max_iter, eps, my_f);*/
	
}
