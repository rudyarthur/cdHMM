#pragma once

#include <math.h>       
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

double gamma_solve(double pm, double start, double end, int max_iter, double eps){
	
	vector<double> p = {pm};
	PsiGamma my_f(p);
	
	vector<double> range = {start, end};
	return brent(range, max_iter, eps, my_f);
	
}
