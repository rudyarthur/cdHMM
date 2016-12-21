#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <limits>
#include <math.h>
#include <stdlib.h>  
#include <random>
#include <algorithm>
#include <iomanip>  
#include <math.h>      
#include <map>      
#include <utility>

using namespace std;

#define EULER_MASCHERONI 0.57721566490153286060651209008240243104215933593992

namespace cdHMM {

struct max_lhood_params {
	
	int dim;
	//start vector (start range for brent or start point for simplex dim = x_start.size() ) 
	vector<double> x_start;
	//stopping criteria
	double eps;
	int max_iter;
	
};

template<typename S, typename T> bool pair_sort_first (pair<S,T> i, pair<S,T> j) { return (i.first>j.first); }
template<typename S, typename T> bool pair_sort_second (pair<S,T> i, pair<S,T> j) { return (i.second>j.second); }

//a = log x
//b = log y
// log( 1 + exp(b - a) ) ; 
// log( 1 + exp( log(y/x) ) )
// log( (x+y) ) - log(x)
// a = log( x+ y )
// if
// a = -inf = log(0)
// a = log(0 + y) = log(y) = b
// if
// b = -inf = log(0)
// a = log(x + 0) = log(x) = a
void lsum(double &a, double b){ 
	if( a == -numeric_limits<double>::infinity() ){ 
		a = b;
	} else if ( b != -numeric_limits<double>::infinity() ){ 
		if( a > b ){
			a += log ( 1.0  + exp(b - a) ); 
		} else {
			a = b + log ( 1.0  + exp(a - b) ); 
		}
	} 
}

void ldiff(double &a, double b){ 
	if( a > b ){
		a += log ( 1.0  - exp(b - a) ); 
	} else if (a == b ){
		a = -numeric_limits<double>::infinity();
	} else {
		cerr << "negative log error " << a << " " << b << endl; exit(1);
	}
}

/* v is a vector of positive and negative values
 * safe log sum returns
 *  ( sign(sum(v))   ,    | log( | sum(v, v>0) - sum(v, v<0) | ) | )
 * */
pair<int, double> safe_log_sum(vector<double> &v){
	int t = 0;
    double pos = (v[t] >= 0) ? log( v[t] ) : -numeric_limits<double>::infinity(); //log v
    double neg = (v[t] < 0) ? log( -v[t] ) : -numeric_limits<double>::infinity(); //log |v|
	for(t=0; t<v.size(); ++t){
		if( v[t] >= 0 ){
			lsum( pos , log( v[t] ) ); //log ( sum v[t] )
		} else {
			lsum( neg , log( -v[t] ) );//log ( sum |v[t]| )
		}
	}
	if( neg > pos ){	//log ( sum |v[t]| ) > log ( sum v[t] )
	  double res = neg;
	  ldiff(res, pos);  //log ( sum |v[t]| - sum v[t] )
	  return make_pair(-1,res); 
    } else { //log ( sum v[t] ) < log ( sum |v[t]| )
	  double res = pos;
	  ldiff(res, neg); //log ( sum v[t] - sum |v[t]| )
	  return make_pair(1,res);  
    }
    cerr << "safe_log_sum error " << neg << " " << pos << endl;
    return make_pair(1,0); 
}

/* v is a vector of positive and negative values
 * add is a vector of ( log(a0), log(a1), ... ) := log(a)
 * safe log sum returns
 *  ( sign(sum(v))   ,    | log( | sum(ai*vi, v>0) - sum(ai*vi, v<0) | ) | )
 * */
pair<int, double> safe_log_sum(vector<double> &v, vector<double> &add){
	int t = 0;
    double pos = (v[t] >= 0) ? add[t] + log( v[t] ) : -numeric_limits<double>::infinity(); //log v
    double neg = (v[t] < 0) ? add[t] + log( -v[t] ) : -numeric_limits<double>::infinity(); //log |v|
	for(t=0; t<v.size(); ++t){
		if( v[t] >= 0 ){
			lsum( pos , add[t] + log( v[t] ) ); //log ( sum v[t] )
		} else {
			lsum( neg , add[t] + log( -v[t] ) );//log ( sum |v[t]| )
		}
	}
	if( neg > pos ){	//log ( sum |v[t]| ) > log ( sum v[t] )
	  double res = neg;
	  ldiff(res, pos);  //log ( sum |v[t]| - sum v[t] )
	  return make_pair(-1,res); 
    } else { //log ( sum v[t] ) < log ( sum |v[t]| )
	  double res = pos;
	  ldiff(res, neg); //log ( sum v[t] - sum |v[t]| )
	  return make_pair(1,res);  
    }
    cerr << "safe_log_sum error " << neg << " " << pos << endl;
    return make_pair(1,0); 
}

}
