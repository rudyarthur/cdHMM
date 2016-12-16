#pragma once

#include <vector>

using namespace std;

struct fit_params{
	int iter;
	double residual;
	double root;
	vector<double> mroot;
};

struct weight_params
{
	//parameters for evaluating the function
    vector<double>* O;				//observations
    vector< vector<double> >* gamma; //weights
    vector<double>* sumgamma;		//sum weights
    int i;							//state    
};

class Func1d{
public:

	int dim;
	weight_params wp;
	vector<double> params;
	virtual double eval(double x) = 0;

};

class FuncNd{
public:

	int dim;
	weight_params wp;
	vector<double> params;
	virtual double eval(vector<double> &x) = 0;

};
