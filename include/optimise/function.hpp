#pragma once

#include <vector>

using namespace std;

namespace cdHMM {

struct fit_params{
	int iter;
	double residual;
	double root;
	vector<double> mroot;
};

template <typename T> class weight_params
{
public:
	//parameters for evaluating the function
    vector<T>* O;				//observations
    vector< vector<double> >* gamma; //weights
    vector<double>* sumgamma;		//sum weights
    int i;							//state    
};

template <typename T> class Func1d{
public:

	int dim;
	weight_params<T> wp;
	vector<double> params;
	virtual double eval(double x) = 0;

};

template <typename T> class FuncNd{
public:

	int dim;
	weight_params<T> wp;
	vector<double> params;
	virtual double eval(vector<double> &x) = 0;

};

}
