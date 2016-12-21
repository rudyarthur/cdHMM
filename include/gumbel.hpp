#pragma once

#include <math.h>       
#include "optimise/function.hpp"
#include "optimise/brent.hpp"

using namespace std;

namespace cdHMM {

template <typename S> class Gumbel : public Func1d<S>{
public:
	
	Gumbel(){ 
		this->dim = 1;
	}
	
	double eval(double b){
			
		  double av = 0;
		  double Eexp = 0;
		  double EOexp = 0;
		  int T = this->wp.O->size();
		  int i = this->wp.i;
		  int t = 0;
		  av = ( (*this->wp.gamma)[i][t] * (double)(*this->wp.O)[t] ); //E[O]
		  double tmp = (- (double)(*this->wp.O)[t] / b);
		  Eexp = log( (*this->wp.gamma)[i][t] ) + ( tmp ); //E[ exp(-O/b) ]
		  
		  //x >= 0: log( sum x exp( -x/b ) )
		  //x < 0: log( sum |x| exp( -x/b ) )
		  double EOexp_pos = ( (double)(*this->wp.O)[t] >= 0) ? log( (*this->wp.gamma)[i][t] ) + ( tmp ) + log( (double)(*this->wp.O)[t] )  : -numeric_limits<double>::infinity(); //E[ O exp(-O/b) ]
		  double EOexp_neg = ( (double)(*this->wp.O)[t] < 0)  ? log( (*this->wp.gamma)[i][t] ) + ( tmp ) + log( -(double)(*this->wp.O)[t] ) : -numeric_limits<double>::infinity(); //E[ O exp(-O/b) ]
		  
		  for( int t=1; t<T; ++t){
			double logg = log( (*this->wp.gamma)[i][t] );
			av += (*this->wp.gamma)[i][t] * (double)(*this->wp.O)[t]; //E[O]
			double tmp = (- (double)(*this->wp.O)[t] / b);
			lsum( Eexp , logg + ( tmp ) ); //E[ exp(-O/b) ]
			if( (double)(*this->wp.O)[t] >= 0 ){
				lsum( EOexp_pos , logg + ( tmp ) + log( (double)(*this->wp.O)[t] ) );
			} else {
				lsum( EOexp_neg , logg + ( tmp ) + log( -(double)(*this->wp.O)[t] ) );
			}
		  }
		  av /= (*this->wp.sumgamma)[i]; 

		  //log( sum x exp( -x/b ) ) < 0 => EOexp_neg > EOexp_pos
		  //ldiff( EOexp_neg - EOexp_pos ) = log( -sum x exp( -x/b ) ) = EOexp
		  //exp(EOexp - Eexp) = -sum( x e^{-x/b} ) / sum( e^{-x/b} )
		  if( EOexp_neg > EOexp_pos ){
			  EOexp = EOexp_neg;
			  if( EOexp != EOexp ){
				  cerr << "error" << endl; exit(1);
			  }
			  ldiff(EOexp, EOexp_pos);
			  return b - exp(EOexp - Eexp) - (av); 
		  } else {
			  EOexp = EOexp_pos;
			  ldiff(EOexp, EOexp_neg);
			  return b + exp(EOexp - Eexp) - (av);  
		  }
		  return b - (av); 
	}

};

template <typename T> fit_params gumbel_solve(vector<T> &O, vector< vector<double> > &gamma, vector<double> &sumgamma, int i,
max_lhood_params mlp){
	
	Gumbel<T> my_f;
	my_f.wp.O = &O;
	my_f.wp.gamma = &gamma;
	my_f.wp.sumgamma = &sumgamma;
	my_f.wp.i = i;
	

	return brent<T>(mlp.x_start, mlp.max_iter, mlp.eps, my_f);

}

}
