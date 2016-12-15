#pragma once

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include "utils.hpp"
#include "hmm.hpp"

using namespace std;



//Can construct anything we like with this
struct general_params
{
    vector<double> O;
    vector< vector<double> > gamma;
    vector<double> sumgamma;
    int i;
};

struct general_params_kterms
{
    vector<double> O;
    vector< vector<double> > gamma;
    vector<double> sumgamma;
    int i;
    vector<double> k_terms;
};


/**********************
 * 
 * GEV stuff
 * 
 * ********************/

double
extremevalue_lhood_f (const gsl_vector *v, void *params)
{

  struct general_params_kterms *p 
    = (struct general_params_kterms *) params;
      
  double mu = gsl_vector_get(v, 0);
  double sigma = gsl_vector_get(v, 1);
  double k = gsl_vector_get(v, 2);
  
  int T = p->O.size();
  int i = p->i;
  
  vector<double> a(T); for(int t=0; t<T; ++t){ a[t] = (p->O[t] - mu)/sigma; } //a = (x-m)/s

  double nlogL = 0;
  if(k == 0){ //gumbel case
	for(int t=0; t<T; ++t){ 
		nlogL += (exp(-a[t]) + a[t] + log(sigma))*p->gamma[i][t]; 
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

      p->k_terms.resize(T);  for(int t=0; t<T; ++t){ p->k_terms[t] = 1; }
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
			p->k_terms[t] += newterm;
		}
        if( max_nt <= 1e-15 ){
          break;
	    }
      }
      //a[t] (k+1) k_terms[t] = k(x-m)/s
      for(int t=0; t<T; ++t){ 
		  nlogL += (exp( - a[t] * p->k_terms[t] ) + a[t] * (k+1) * p->k_terms[t] + log(sigma))*p->gamma[i][t]; 
	  }
    } else {
      vector<double> b(T); for(int t=0; t<T; ++t){ b[t] = 1 + aa[t]; if(b[t] <= 0 ){ return 1e200; } } //b = 1 + k(x-m)/s
	  for(int t=0; t<T; ++t){ 
		  //cout << p->O[t] << " " << -(pow( b[t], -1/k ) + (1+1/k) * log(b[t]) + log(sigma)) << endl; exit(1);
		  nlogL += (pow( b[t], -1/k ) + (1+1/k) * log(b[t]) + log(sigma))*p->gamma[i][t]; 
	  }
    }
  }
  return nlogL/p->sumgamma[i];

}

void
extremevalue_lhood_df (const gsl_vector *v, void *params, 
       gsl_vector *df)
{
  struct general_params_kterms *p 
    = (struct general_params_kterms *) params;
      
  double m = gsl_vector_get(v, 0);
  double sigma = gsl_vector_get(v, 1);
  double k = gsl_vector_get(v, 2);
  
  
  int T = p->O.size();
  int i = p->i;
  

  vector<double> G(3,0); 
  vector<double> a(T); for(int t=0; t<T; ++t){ a[t] = (p->O[t] - m)/sigma; } //a = (x-m)/s

  if ( k == 0  ){ 
	vector<double> f(T); for(int t=0; t<T; ++t){ f[t] = exp(-a[t]) - 1; } 
	vector<double> g(T); 
	for(int t=0; t<T; ++t){ 
		G[2] += ( a[t] * (1.0 + a[t] * f[t] * 0.5) )*p->gamma[i][t]; ;
		G[1] += ( (a[t] * f [t] + 1.0) / sigma )*p->gamma[i][t]; ;
		G[0] += ( f[t] / sigma )*p->gamma[i][t]; ;
	} 
  } else {
	bool inf_val = false;
	vector<double> b(T), c(T); for(int t=0; t<T; ++t){ 
		b[t] = 1 + k * a[t]; 
		if( b[t] <= 0 ){
			G[0] = 0; G[1] = 0; G[2] = 0; inf_val = true; 
			break;
		}
		c[t] = log(b[t]);
	}
	if(!inf_val){
		vector<double> aa(T); //aa = k(x-m)/s
		double min_aa = numeric_limits<double>::infinity();
		double max_aa = -numeric_limits<double>::infinity();
		for(int t=0; t<T; ++t){ 
			aa[t] = k * a[t]; 
			if( fabs(aa[t]) < min_aa ){ min_aa = fabs(aa[t]); }
			if( fabs(aa[t]) > max_aa ){ max_aa = fabs(aa[t]); }
			
		}
		if( min_aa < 1e-3 && max_aa < 0.5 ){
		
		//if( p->k_terms.size() > 0 ){ //need more accurate log
			vector<double> kk_terms(T); for(int t=0; t<T; ++t){ kk_terms[t] = 0.5; }
			double sgn = 1; 
			int it = 0;
			while(true){
				sgn = -sgn; it++;
				double max_nt = -numeric_limits<double>::infinity();
				for(int t=0; t<T; ++t){ 
					double newterm = (sgn * (it + 1) / (it + 2)) * (  pow(k*a[t] , it) );
					if( fabs(newterm) > max_nt ){ max_nt = fabs(newterm); }
					kk_terms[t] += newterm;
				}

				if( max_nt <= 1e-15 ){
				  break;
				}
			}
			for(int t=0; t<T; ++t){ 
				double f = exp(-a[t] * p->k_terms[t]);
				G[2] += ( a[t] * ((a[t] * kk_terms[t]) * (f - 1 - k) + p->k_terms[t]) )*p->gamma[i][t]; 
				G[1] += ( (1 - a[t] * (a[t] * k * kk_terms[t] - p->k_terms[t]) * (f - k - 1)) / sigma )*p->gamma[i][t]; 
				G[0] += ( -(a[t] * k * kk_terms[t] - p->k_terms[t]) * (f - k - 1) / sigma )*p->gamma[i][t]; 
			}
		} else {
			double d = 1/k + 1;
			for(int t=0; t<T; ++t){ 
				G[2] += ( (c[t] / k - a[t] / b[t]) / (k * pow( b[t] , (1/k)) ) - c[t] / (k*k) + a[t] * d / b[t] )*p->gamma[i][t]; 
				G[1] += ( (a[t] * pow( b[t] , (-d)) - (k + 1) * a[t] / b[t] + 1) / sigma )*p->gamma[i][t]; 
				G[0] += ( (pow(b[t] , (-d)) - (k + 1) / b[t]) / sigma )*p->gamma[i][t]; 
			}
		}
	}
  }

  gsl_vector_set(df, 0, G[0]/p->sumgamma[i]);
  gsl_vector_set(df, 1, G[1]/p->sumgamma[i]);
  gsl_vector_set(df, 2, G[2]/p->sumgamma[i]);
  //return GSL_SUCCESS;
}

int 
extremevalue_lhood_ddf (const gsl_vector *v, void *params, 
       gsl_matrix * J)
{
  struct general_params_kterms *p 
    = (struct general_params_kterms *) params;
      
  double mu = gsl_vector_get(v, 0);
  double sigma = gsl_vector_get(v, 1);
  double k = gsl_vector_get(v, 2);
  
  
  int T = p->O.size();
  int i = p->i;
  
vector< vector<double> > ACOV(3, vector<double>(3,0) );

if(k == 0){ 
  
  for(int t=0; t<T; ++t){
	double a = (p->O[t] - mu)/sigma;
    double f = exp(-a);
    //k, k
    ACOV[2][2] += (a*a) * (a * (a/4.0 - 2.0/3.0) * f + (2.0/3.0) * a - 1.0);  

    //sigma, sigma
    ACOV[1][1] += (sigma * sigma) * (a * ((a - 2) * f + 2) - 1);

	//mu, mu
	ACOV[0][0] += (sigma * sigma) * f;

    //k, sigma
    ACOV[2][1] += (-a / sigma) * (a * (1 - a/2) * f - a + 1);
    ACOV[1][2] = ACOV[2][1];

	//k, mu
    ACOV[2][0] += (-1 / sigma) * (a * (1 - a/2) * f - a + 1);
	ACOV[0][2] = ACOV[2][0];
	
    //sigma, mu
    ACOV[1][0] += (1 + (a - 1) * f) / (sigma * sigma);
	ACOV[0][1] = ACOV[1][0];
  }
} else {
  double d = 1 / k + 1;
  for(int t=0; t<T; ++t){
	double a = (p->O[t] - mu)/sigma;
	double b = k * a + 1;
    if( b >= 0 ){
		double c = log(b);
		//k, k
		ACOV[2][2] += ((   pow((c / (k*k) ) - (a / (k * b)), 2)  ) / ( pow(b , (1 / k)) )) + 
  ((-2*c / (k*k*k)) + (2*a / (k*k*b)) + ( (a / b) * (a / b) / k)) / ( pow(b , (1 / k)) ) + 
  2*c / (k*k*k) - 
  (2*a / ( k*k*b)) - (d * (a / b) * (a / b) );

		//s,s
		ACOV[1][1] += (sigma * sigma ) * (
    -2*a * pow( b , (-d)) + 
    d * k * a * a * pow(b , (-d-1)) + 
    2 * d * k * a / b - 
    d * (k * a / b) * (k * a / b) - 1);

		//m,m
		ACOV[0][0] += (d / (sigma * sigma)) *  (
  k * pow(b , (-d-1)) - 
  (k / b) * (k / b) );

		//k, mu
		double dterm = ( pow(b , (-d)) * (c / k  - a / b) / k - 
a * pow(b , (-d-1)) + 
((1 / k) - d) / b +
a * k * d / (b * b )) / sigma;
		ACOV[2][0] += dterm;
		ACOV[0][2] = ACOV[2][0];

		//k, sigma
		ACOV[2][1] += a*dterm;
		ACOV[1][2] = ACOV[2][1];

		//sigma, mu
		ACOV[1][0] += ( -pow(b , (-d)) + 
a * k * d * pow(b , (-d-1)) + 
(d * k / b) - a * (k/b) * (k/b) * d) / (sigma * sigma);
		ACOV[0][1] = ACOV[1][0];
	
}

  }
}
/*cout << ACOV[2][2] << " " << ACOV[2][1] << " " << ACOV[2][0] << endl;
cout << ACOV[1][2] << " " << ACOV[1][1] << " " << ACOV[1][0] << endl;
cout << ACOV[0][2] << " " << ACOV[0][1] << " " << ACOV[0][0] << endl;*/
for(int i=0; i<3; ++i){ for(int j=0; j<3; ++j){ gsl_matrix_set (J, i, j, ACOV[i][j]); } }

  return GSL_SUCCESS;

}

/* Compute both f and df together. */
void 
extremevalue_lhood_fdf (const gsl_vector *x, void *params, 
        double *f, gsl_vector *df) 
{
  *f = extremevalue_lhood_f(x, params); 
  extremevalue_lhood_df(x, params, df);
}

int
extremevalue_lhood_2nd_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  extremevalue_lhood_df (x, params, f);
  extremevalue_lhood_ddf (x, params, J);

  return GSL_SUCCESS;
}

  /*gsl_matrix *J;
  gsl_matrix_alloc (3,3);
  cout << xstart[0] << " " << xstart[1] << " " << xstart[2] << endl;
  cout << extremevalue_lhood_f (x, &params) << endl;
  extremevalue_lhood_ddf (x, &params, J);*/
  
int multi_solve(vector<double> &out,
vector<double> &O, vector< vector<double> > &gamma, vector<double> &sumgamma, int i,
vector<double> &xstart, 
int max_iter, double eps){
  
  vector<double> kterms;
  struct general_params_kterms params = { O, gamma, sumgamma, i, kterms };
  /*
  // Gradient
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 3;
  my_func.f = extremevalue_lhood_f;
  my_func.df = extremevalue_lhood_df;
  my_func.fdf = extremevalue_lhood_fdf;
  my_func.params = &params;

  x = gsl_vector_alloc (3);
  gsl_vector_set (x, 0, xstart[0]);
  gsl_vector_set (x, 1, xstart[1]);
  gsl_vector_set (x, 2, xstart[2]);

  T = gsl_multimin_fdfminimizer_conjugate_pr;
  s = gsl_multimin_fdfminimizer_alloc (T, 3);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, eps, eps);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
	  cout << "iter " << iter << " " << (status == GSL_ENOPROG) << endl;
      if (status)
        break;

      status = gsl_multimin_test_gradient (s->gradient, eps);

      cout << iter << " x = " <<
              gsl_vector_get (s->x, 0) << " " <<
              gsl_vector_get (s->x, 1) << " " <<
              gsl_vector_get (s->x, 2) << " g = " <<
              gsl_vector_get (s->gradient, 0) << " " <<
              gsl_vector_get (s->gradient, 1) << " " <<
              gsl_vector_get (s->gradient, 2) << " " <<
              s->f << endl;

    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  if( iter == max_iter ){
	  cout << iter << " " <<
              gsl_vector_get (s->x, 0) << " " <<
              gsl_vector_get (s->x, 1) << " " <<
              gsl_vector_get (s->x, 2) << " " <<
              s->f << endl;
  }
  
  out.resize(3);
  out[0] = gsl_vector_get (s->x, 0);
  out[1] = gsl_vector_get (s->x, 1);
  out[2] = gsl_vector_get (s->x, 2);
  
  
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
  */
  

  //No Gradient
  size_t iter = 0;
  int status;
  double size;

  // Starting point 
  gsl_vector *x, *df;
  x = gsl_vector_alloc (3);
  df = gsl_vector_alloc (3);

  gsl_vector_set (x, 0, xstart[0]);
  gsl_vector_set (x, 1, xstart[1]);
  gsl_vector_set (x, 2, xstart[2]);
       
  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss;
  gsl_multimin_function minex_func;
  
  // Set initial step sizes to 1 
  ss = gsl_vector_alloc (3);
  gsl_vector_set_all (ss, 1.0);

  // Initialize method and iterate 
  minex_func.n = 3;
  minex_func.f = extremevalue_lhood_f;
  minex_func.params = &params;

  s = gsl_multimin_fminimizer_alloc (T, 3);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status){ 
		  iter == max_iter; 
		  cerr << "iterate fail!" << endl;
          break;
      }
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, eps);
      
      //cout << "iter " << iter << " " <<  gsl_vector_get (s->x, 0) << " " << gsl_vector_get (s->x, 1) << " " <<
      //        gsl_vector_get (s->x, 2) << " " << s->fval << " " << size << endl;
              
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_vector_free(ss); 

  //cout << "Iter= " << iter << endl;
  if( iter == max_iter ){
	  cerr << "Max_iter reached!" << endl;
	  cerr << "iter " << iter << " " << 
              gsl_vector_get (s->x, 0) << " " <<
              gsl_vector_get (s->x, 1) << " " <<
              gsl_vector_get (s->x, 2) << " " <<
              s->fval << " " << size << endl;
  } else {
	  out.resize(3);
	  out[0] = gsl_vector_get (s->x, 0);
	  out[1] = gsl_vector_get (s->x, 1);
	  out[2] = gsl_vector_get (s->x, 2);
  }
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);

  

/*
  int status;
  size_t iter = 0;

  const size_t n = 3;
  /*gsl_multiroot_function_fdf f;
  f.n = 3;
  f.f = extremevalue_lhood_df;
  f.df = extremevalue_lhood_ddf;
  f.fdf = extremevalue_lhood_2nd_fdf;
  f.params = &params;*/
  /*
  gsl_multiroot_function f = {&extremevalue_lhood_df, n, &params};
  
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, 9.1);
  gsl_vector_set (x, 1, 1.1);
  gsl_vector_set (x, 2, 0.1);


  /*const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  T = gsl_multiroot_fdfsolver_gnewton;
  s = gsl_multiroot_fdfsolver_alloc (T, n);
  gsl_multiroot_fdfsolver_set (s, &f, x);*/
/*
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 3);
  gsl_multiroot_fsolver_set (s, &f, x);

  do
    {
      //iter++;
      //status = gsl_multiroot_fdfsolver_iterate (s);
	  iter++;
      status = gsl_multiroot_fsolver_iterate (s);

	  cout << "iter " << iter << " " << 
              gsl_vector_get (s->x, 0) << " " <<
              gsl_vector_get (s->x, 1) << " " <<
              gsl_vector_get (s->x, 2) << " " << endl;
              
      if (status)
        break;

      status = gsl_multiroot_test_residual (s->f, 1e-5);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  exit(1);
  //gsl_multiroot_fdfsolver_free (s);
  //gsl_vector_free (x);
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);*/
  
  return iter;
}

