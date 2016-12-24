#pragma once

#include <math.h>       
#include "utils.hpp"
#include "optimise/function.hpp"
#include "optimise/brent.hpp"
#include "optimise/psi.hpp"

using namespace std;

namespace cdHMM {

/*!
 * digamma function, gamma(x)/gamma'(x)
 * */
class PsiGamma : public Func1d<double>{
public:
	
/*!
 * Construct digamma function
 * */
	PsiGamma(vector<double> &p){ 
		if( p.size() != 1 ){
			cerr << "Gamma requires 1 parameter!" << endl; 
			exit(1);
		}
		this->params = p; 
	}
/*!
 * Evaluate digamma function
 * */
	double eval(double x){
			return psi (x) - log(x) + this->params[0];
	}

};

/*!
 * Max likelihood estimate of gamma distribution parameters given pm = E[ log[O] ]
 * */
fit_params gamma_solve(double pm, max_lhood_params mlp){ 
	
	vector<double> p = {pm};
	PsiGamma my_f(p);
	
	return brent(mlp.x_start, mlp.max_iter, mlp.eps, my_f);
	
}

}
