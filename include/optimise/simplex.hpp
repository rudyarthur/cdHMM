/* from MYGSL
/* multimin/simplex2.c
 * 
 * Copyright (C) 2007, 2008, 2009 Brian Gough
 * Copyright (C) 2002 Tuomo Keskitalo, Ivo Alxneit
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
   - Originally written by Tuomo Keskitalo <tuomo.keskitalo@iki.fi>
   - Corrections to nmsimplex_iterate and other functions 
     by Ivo Alxneit <ivo.alxneit@psi.ch>
   - Additional help by Brian Gough <bjg@network-theory.co.uk>
   - Optimisations added by Brian Gough <bjg@network-theory.co.uk>
         + use BLAS for frequently-called functions
         + keep track of the center to avoid unnecessary computation
         + compute size as RMS value, allowing linear update on each step
           instead of recomputing from all N+1 vectors.
*/

/* The Simplex method of Nelder and Mead, also known as the polytope
   search alogorithm.  Ref: Nelder, J.A., Mead, R., Computer Journal 7
   (1965) pp. 308-313.

   This implementation uses n+1 corner points in the simplex.
*/
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "function.hpp"

using namespace std;

#define MYGSL_SUCCESS 0
#define MYGSL_CONTINUE 2
#define MYGSL_EBADFUNC 9

static void dscal  (double alpha, vector<double> &X){
	for(int i=0; i<X.size(); ++i){ X[i] *= alpha; }
}
static void daxpy  (double alpha, const vector<double> &X, vector<double> &Y){
	for(int i=0; i<X.size(); ++i){ Y[i] += alpha * X[i]; }
}
static double ddot (vector<double> &X, vector<double> &Y){
	double result = 0; for(int i=0; i<X.size(); ++i){ result += Y[i] * X[i]; } return result;
}
static double dnrm2  (vector<double> &X){
	double result = 0; for(int i=0; i<X.size(); ++i){ result += X[i] * X[i]; } return sqrt(result);
}
static void set_zero  (vector<double> &X){
	for(int i=0; i<X.size(); ++i){ X[i] = 0; } 
}
static int min_index (vector<double> &v){
	double min = v[0]; int min_idx = 0;
	for(int i=0; i<v.size(); ++i){
		if( v[i] < min ){
			min = v[i]; min_idx = i;
		}
	}
	return min_idx;
}

/* minimization of non-differentiable functions */
typedef struct
{
  vector<vector<double> > x1;		/* simplex corner points */
  vector<double> y1;		/* function value at corner points */
  vector<double> ws1;		/* workspace 1 for algorithm */
  vector<double> ws2;		/* workspace 2 for algorithm */
  vector<double> center;		/* center of all points */
  vector<double> delta;		/* current step */
  vector<double> xmc;		/* x - center (workspace) */
  double S2;
  unsigned long count;
}
nmsimplex_state_t;

typedef struct 
{
  //const char *name;
  //size_t size;
  int alloc (nmsimplex_state_t &state, size_t n);
  int set (nmsimplex_state_t &state, FuncNd *f,
              vector<double> &x, 
              double * size,
              vector<double> &step_size);
  int iterate (nmsimplex_state_t &state, FuncNd *f, 
                  vector<double> &x, 
                  double * size,
                  double * fval);

}
gsl_multimin_fminimizer_type;



typedef struct 
{
  /* multi dimensional part */
  gsl_multimin_fminimizer_type type;
  FuncNd *f;

  double fval;
  vector<double> x;
  
  double size;

  nmsimplex_state_t state;
}
gsl_multimin_fminimizer;

gsl_multimin_fminimizer 
gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_type &T,
                               size_t n)
{
  int status;
  
  gsl_multimin_fminimizer s;
	//= (gsl_multimin_fminimizer *) malloc (sizeof (gsl_multimin_fminimizer));

 /* if (s == 0)
    {
	  cerr << "failed to allocate space for minimizer struct" << endl; exit(1);
    }*/

  s.type = T;

  s.x.resize(n);

  /*s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);
      cerr << "failed to allocate space for minimizer state" << endl; exit(1);
    }*/

  status = T.alloc(s.state, n);

  /*if (status != MYGSL_SUCCESS)
    {
      free (s->state);
      free (s);
      cerr << "failed to initialize minimizer state" << endl; exit(1);
    }*/

  return s;
}

int
gsl_multimin_fminimizer_set (gsl_multimin_fminimizer &s,
                             FuncNd &f,
                             vector<double> &x,
                             vector<double> &step_size)
{

  if (s.x.size() != f.dim)
    {
	  cerr << "function incompatible with solver size" << endl; exit(1);
    }

  if (x.size() != f.dim || step_size.size() != f.dim) 
    {
	  cerr << "vector length not compatible with function" << endl; exit(1);
    }  
    
  s.f = &f;

  (s.x) = (x); 
  
  return s.type.set(s.state, s.f, s.x, &(s.size), step_size);
}

int
gsl_multimin_fminimizer_iterate (gsl_multimin_fminimizer &s)
{
  return s.type.iterate(s.state, s.f, s.x, &(s.size), &(s.fval));
}

double
gsl_multimin_fminimizer_size (gsl_multimin_fminimizer &s)
{
  return s.size;
}

int
gsl_multimin_test_size (const double size, double epsabs)
{
  if (epsabs < 0.0)
    {
	  cerr << "absolute tolerance is negative" << endl; exit(1);
    }
  
  if (size < epsabs)
    {
      return MYGSL_SUCCESS;
    }

  return MYGSL_CONTINUE;
}



static int
compute_center (const nmsimplex_state_t &state, vector<double> &center);
static double
compute_size (nmsimplex_state_t &state, vector<double> &center);

static double
try_corner_move (const double coeff,
		 const nmsimplex_state_t &state,
		 size_t corner,
		 vector<double> &xc, FuncNd *f)
{
  /* moves a simplex corner scaled by coeff (negative value represents 
     mirroring by the middle point of the "other" corner points)
     and gives new corner in xc and function value at xc as a 
     return value 
   */
  
  const size_t P = state.x1.size();
  double newval;

  /* xc = (1-coeff)*(P/(P-1)) * center(all) + ((P*coeff-1)/(P-1))*x_corner */
  {
    double alpha = (1 - coeff) * P / (P - 1.0);
    double beta = (P * coeff - 1.0) / (P - 1.0);
    vector<double> row = state.x1[corner]; 
    
	xc = (state.center);
    dscal (alpha, xc);
    daxpy (beta, row, xc);
  }

  newval = f->eval(xc);

  return newval;
}


static void
update_point (nmsimplex_state_t &state, size_t i,
	      vector<double> &x, double val)
{
  const size_t P = state.x1.size();

  /* Compute delta = x - x_orig */
  state.delta = x;
  daxpy(-1.0, state.x1[i], state.delta);
  
  /* Compute xmc = x_orig - c */
  state.xmc = state.x1[i];
  daxpy(-1, state.center, state.xmc);

  /* Update size: S2' = S2 + (2/P) * (x_orig - c).delta + (P-1)*(delta/P)^2 */
  {
    double d = dnrm2 (state.delta);
    double xmcd;

    xmcd = ddot (state.xmc, state.delta);
    state.S2 += (2.0 / P) * xmcd + ((P - 1.0) / P) * (d * d / P);
  }

  /* Update center:  c' = c + (x - x_orig) / P */

  {
    double alpha = 1.0 / P;
    daxpy (-alpha, state.x1[i], state.center);
    daxpy (alpha, x, state.center);
  }

  state.x1[i] = x; 
  state.y1[i] = val;  
}

static int
contract_by_best (nmsimplex_state_t &state, size_t best,
		  vector<double> &xc, FuncNd *f)
{

  /* Function contracts the simplex in respect to best valued
     corner. That is, all corners besides the best corner are moved.
     (This function is rarely called in practice, since it is the last
     choice, hence not optimised - BJG)  */

  /* the xc vector is simply work space here */

  
  size_t i, j;
  double newval;

  int status = MYGSL_SUCCESS;

  for (i = 0; i < state.x1.size(); i++)
    {
      if (i != best)
	{
	  for (j = 0; j < state.x1[i].size(); j++)
	    {
	      state.x1[i][j] = 0.5 * (state.x1[i][j] + state.x1[best][j]); 
	    }

	  /* evaluate function in the new point */
	  state.y1[i] = f->eval( state.x1[i] );

	  /* notify caller that we found at least one bad function value.
	     we finish the contraction (and do not abort) to allow the user
	     to handle the situation */

	  if (!isfinite (newval))
	    {
	      status = MYGSL_EBADFUNC;
	    }
	}
    }

  
  /* We need to update the centre and size as well */
  compute_center (state, state.center);
  compute_size (state, state.center);

  return status;
}

static int
compute_center (const nmsimplex_state_t &state, vector<double> &center)
{
  /* calculates the center of the simplex and stores in center */  
  const size_t P = state.x1.size();
  size_t i;

  set_zero (center); 
  
  for (i = 0; i < P; i++)
    {
      daxpy(1.0, state.x1[i], center);
    }

  {
    const double alpha = 1.0 / P;
    dscal (alpha, center);
  }

  return MYGSL_SUCCESS;
}

static double
compute_size (nmsimplex_state_t &state, vector<double> &center)
{
  /* calculates simplex size as rms sum of length of vectors 
     from simplex center to corner points:     

     sqrt( sum ( || y - y_middlepoint ||^2 ) / n )
   */

  vector<vector<double> > x1 = (state.x1); 
  const size_t P = state.x1.size();
  size_t i;

  double ss = 0.0;

  for (i = 0; i < P; i++)
    {
      double t;
      state.ws1 = state.x1[i]; 
      daxpy(-1.0, center, state.ws1); 
      t = dnrm2(state.ws1); 
      ss += t * t;
    }

  /* Store squared size in the state */
  state.S2 = (ss / P);
  
  return sqrt (ss / P);
}

int
gsl_multimin_fminimizer_type::alloc (nmsimplex_state_t &state, size_t n)
{
  //nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  if (n == 0)
    {
      cerr << "invalid number of parameters specified" << endl;  exit(1);
    }

  state.x1 = vector< vector<double> >(n+1, vector<double>(n,0) );
  state.y1.resize(n+1);
  state.ws1.resize(n);
  state.ws2.resize(n);
  state.center.resize(n);
  state.delta.resize(n);
  state.xmc.resize(n);
  
  state.count = 0;

  return MYGSL_SUCCESS;
}

int
gsl_multimin_fminimizer_type::set (nmsimplex_state_t &state, FuncNd *f,
	       vector<double> &x,
	       double *size, vector<double> &step_size)
{
  int status;
  size_t i;
  double val;


  //nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  vector<double> xtemp = state.ws1; 
  if (xtemp.size() != x.size())
    {
      cerr << "incompatible size of x" << endl; exit(1);
    }

  if (xtemp.size() != step_size.size() )
    {
      cerr << "incompatible size of step_size" << endl; exit(1);
    }

  /* first point is the original x0 */
  val = f->eval( x );

  if (!isfinite (val))
    {
      cerr << "non-finite function value encountered" << endl; exit(1);
    }

  state.x1[0] = x;
  state.y1[0] = val;

  /* following points are initialized to x0 + step_size */

  for (i = 0; i < x.size(); i++)
    {
      xtemp = x; 

      {
	double xi = x[i];
	double si = step_size[i]; 

	xtemp[i] = xi + si; 
	val = f->eval(xtemp);
      }

      if (!isfinite (val))
	{
	  cerr << "non-finite function value encountered" << endl; exit(1);
	}
	  state.x1[i+1] = xtemp;
      state.y1[i+1] = val;
    }

  compute_center (state, state.center);

  /* Initialize simplex size */
  *size = compute_size (state, state.center);

  state.count++;

  return MYGSL_SUCCESS;
}

int
gsl_multimin_fminimizer_type::iterate (nmsimplex_state_t &state, FuncNd *f,
		   vector<double> &x, double *size, double *fval)
{

  /* Simplex iteration tries to minimize function f value */
  /* Includes corrections from Ivo Alxneit <ivo.alxneit@psi.ch> */

  //nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  /* xc and xc2 vectors store tried corner point coordinates */
  vector<double> xc = state.ws1; 
  vector<double> xc2 = state.ws2;  


  const size_t n = state.y1.size();
  size_t i;
  size_t hi, s_hi, lo;
  double dhi, ds_hi, dlo;
  int status;
  double val, val2;

  if (xc.size() != x.size() )
    {
      cerr << "incompatible size of x" << endl; exit(1);
    }

  /* get index of highest, second highest and lowest point */

  dhi = dlo = state.y1[0]; 
  hi = 0;
  lo = 0;

  ds_hi = state.y1[1]; 
  s_hi = 1;

  for (i = 1; i < n; i++)
    {
      val = state.y1[i]; 
      if (val < dlo)
	{
	  dlo = val;
	  lo = i;
	}
      else if (val > dhi)
	{
	  ds_hi = dhi;
	  s_hi = hi;
	  dhi = val;
	  hi = i;
	}
      else if (val > ds_hi)
	{
	  ds_hi = val;
	  s_hi = i;
	}
    }

  /* try reflecting the highest value point */
  val = try_corner_move (-1.0, state, hi, xc, f);
  
  if (isfinite (val) && val < state.y1[lo]) 
    {
      /* reflected point is lowest, try expansion */
      val2 = try_corner_move (-2.0, state, hi, xc2, f);

      if (isfinite (val2) && val2 < state.y1[lo])
	{
	  update_point (state, hi, xc2, val2);
	}
      else
	{
	  update_point (state, hi, xc, val);
	}
    }
  else if (!isfinite (val) || val > state.y1[s_hi]) 
    {
      /* reflection does not improve things enough, or we got a
         non-finite function value */

      if (isfinite (val) && val <= state.y1[hi]) 
	{
	  /* if trial point is better than highest point, replace
	     highest point */
	  update_point (state, hi, xc, val);
	}

      /* try one-dimensional contraction */
      val2 = try_corner_move (0.5, state, hi, xc2, f);

      if (isfinite (val2) && val2 <= state.y1[hi]) 
	{
	  update_point (state, hi, xc2, val2);
	}
      else
	{
	  /* contract the whole simplex about the best point */
	  status = contract_by_best (state, lo, xc, f);

	  if (status != MYGSL_SUCCESS)
	    {
	      cerr << "contraction failed" << endl; exit(1); 
	    }
	}
    }
  else
    {
      /* trial point is better than second highest point.  Replace
         highest point by it */
      update_point (state, hi, xc, val);
    }

  /* return lowest point of simplex as x */
  lo = min_index(state.y1);
  x = state.x1[lo];
  
  *fval = state.y1[lo]; 
  
  /* Update simplex size */
  {
    double S2 = state.S2;

    if (S2 > 0)
      {
	*size = sqrt (S2);
      }
    else
      {
	/* recompute if accumulated error has made size invalid */
	*size = compute_size (state, state.center);
      }
  }

  return MYGSL_SUCCESS;
}

static const gsl_multimin_fminimizer_type nmsimplex_type = 
{ //"nmsimplex2",	//name 
  //sizeof (nmsimplex_state_t),
  //&nmsimplex_alloc,
  //&nmsimplex_set,
  //&nmsimplex_iterate,
};

gsl_multimin_fminimizer_type gsl_multimin_fminimizer_nmsimplex2 = nmsimplex_type;

fit_params simplex(vector<double> &x, int max_iter, double eps, FuncNd &f){
	
  fit_params fp;
  fp.root = 0; //not used
  fp.iter = 0;
  
  int status;
  
  gsl_multimin_fminimizer_type T = 
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer s;// = NULL;


  // Set initial step sizes to 1 
  vector<double> ss( x.size() , 1.0 );

  // Initialize method and iterate 
  s = gsl_multimin_fminimizer_alloc (T, f.dim);
  gsl_multimin_fminimizer_set (s, f, x, ss);

  do
    {
      fp.iter++;
      status = gsl_multimin_fminimizer_iterate(s);
         
      if (status) 
        break;

      fp.residual = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (fp.residual, eps);
      
      
     /* if (status == MYGSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5d %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", 
              fp.iter,
              s.x[0], 
              s.x[1], 
              s.x[2], 
              s.fval, fp.residual);
       */       
    }
  while (status == MYGSL_CONTINUE && fp.iter < max_iter);
  fp.mroot = s.x;

  return fp;
}
