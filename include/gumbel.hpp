#pragma once

#include <math.h>       
#include "optimise/function.hpp"
#include "optimise/brent.hpp"

using namespace std;

class Gumbel : public Func1d{
public:
	
	Gumbel(){ 
		this->dim = 1;
	}
	
	double eval(double b){
			
		  double av = 0;
		  double Eexp = 0;
		  double EOexp = 0;
		  int T = wp.O->size();
		  int i = wp.i;
		  int t = 0;
		  av = ( (*wp.gamma)[i][t] * (*wp.O)[t] ); //E[O]
		  double tmp = (- (*wp.O)[t] / b);
		  Eexp = log( (*wp.gamma)[i][t] ) + ( tmp ); //E[ exp(-O/b) ]
		  
		  //x >= 0: log( sum x exp( -x/b ) )
		  //x < 0: log( sum |x| exp( -x/b ) )
		  double EOexp_pos = ( (*wp.O)[t] >= 0) ? log( (*wp.gamma)[i][t] ) + ( tmp ) + log( (*wp.O)[t] )  : -numeric_limits<double>::infinity(); //E[ O exp(-O/b) ]
		  double EOexp_neg = ( (*wp.O)[t] < 0)  ? log( (*wp.gamma)[i][t] ) + ( tmp ) + log( -(*wp.O)[t] ) : -numeric_limits<double>::infinity(); //E[ O exp(-O/b) ]
		  
		  for( int t=1; t<T; ++t){
			double logg = log( (*wp.gamma)[i][t] );
			av += (*wp.gamma)[i][t] * (*wp.O)[t]; //E[O]
			double tmp = (- (*wp.O)[t] / b);
			lsum( Eexp , logg + ( tmp ) ); //E[ exp(-O/b) ]
			if( (*wp.O)[t] >= 0 ){
				lsum( EOexp_pos , logg + ( tmp ) + log( (*wp.O)[t] ) );
			} else {
				lsum( EOexp_neg , logg + ( tmp ) + log( -(*wp.O)[t] ) );
			}
		  }
		  av /= (*wp.sumgamma)[i]; 

		  //log( sum x exp( -x/b ) ) < 0 => EOexp_neg > EOexp_pos
		  //ldiff( EOexp_neg - EOexp_pos ) = log( -sum x exp( -x/b ) ) = EOexp
		  //exp(EOexp - Eexp) = -sum( x e^{-x/b} ) / sum( e^{-x/b} )
		  //cout << "hi 1 " << (EOexp_neg) << " " << (EOexp_pos) << endl;
		  if( EOexp_neg > EOexp_pos ){
			  EOexp = EOexp_neg;
			  if( EOexp != EOexp ){
				  cerr << "error" << endl; exit(1);
			  }
			  ldiff(EOexp, EOexp_pos);
			  //cout << "n>p return " << b - exp(EOexp - Eexp) - (av) << endl; exit(1);
			  return b - exp(EOexp - Eexp) - (av); 
		  } else {
			  EOexp = EOexp_pos;
			  ldiff(EOexp, EOexp_neg);
			  //cout << "p>n return " << b + exp(EOexp - Eexp) - (av) << endl; exit(1);
			  return b + exp(EOexp - Eexp) - (av);  
		  }
		  return b - (av); 
	}

};

double gumbel_solve(vector<double> &O, vector< vector<double> > &gamma, vector<double> &sumgamma, int i,
double start, double end, int max_iter, double eps){
	
	Gumbel my_f;
	my_f.wp.O = &O;
	my_f.wp.gamma = &gamma;
	my_f.wp.sumgamma = &sumgamma;
	my_f.wp.i = i;
	
	vector<double> range = {start, end};
	return brent(range, max_iter, eps, my_f);
	return 0;
	
}