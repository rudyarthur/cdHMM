#pragma once

#include <math.h>       
#include "utils.hpp"
#include "optimise/function.hpp"
#include "optimise/brent.hpp"
#include "optimise/psi.hpp"

using namespace std;

class PsiGamma : public Func1d{
public:
	
	PsiGamma(vector<double> &p){ 
		if( p.size() != 1 ){
			cerr << "Gamma requires 1 parameter!" << endl; 
			exit(1);
		}
		this->params = p; 
		this->dim = 1;
	}
	
	double eval(double x){
			return psi (x) - log(x) + this->params[0];
	}

};

double gamma_solve(double pm, max_lhood_params mlp){ 
	
	vector<double> p = {pm};
	PsiGamma my_f(p);
	
	return brent(mlp.x_start, mlp.max_iter, mlp.eps, my_f);
	
}
