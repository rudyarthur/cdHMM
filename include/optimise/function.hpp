#pragma once

#include <vector>

using namespace std;


namespace cdHMM {

/*!
 * Container for return values of optimisation functions
 * */
class fit_params{
public:
	int iter;  /*!< number of iterations*/
	double residual; /*!<error*/
	double root; /*!<max lhood value*/
	vector<double> mroot; /*!<max lhood values*/
};

/*!
 * Container for parameter values for probability functions
 * */
template <typename T> class weight_params
{
public:
    vector<T>* O;				/*!< observations*/
    vector< vector<double> >* gamma; /*!< weights*/
    vector<double>* sumgamma;		/*!< sum weights*/
    int i;							/*!< state    */
};

/*!
 * One dimensional function interface
 * */
template <typename T> class Func1d{
public:

	weight_params<T> wp;	/*!< observation dependent parameters*/
	vector<double> params; /*!< observation independent parameters*/
	/*! Function value
	 * @param[in] x current best parameter
	 * */
	virtual double eval(double x) = 0; 

};

/*!
 * N dimensional function interface
 * */
template <typename T> class FuncNd{
public:

	int dim; /*!< Function dimensions*/
	weight_params<T> wp;	/*!< observation dependent parameters*/
	vector<double> params; /*!< observation independent parameters*/
	/*! Function value
	* @param[in] x current best parameter
	* */
	virtual double eval(vector<double> &x) = 0;

};

}
