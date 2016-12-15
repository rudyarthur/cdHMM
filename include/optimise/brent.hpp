#pragma once

#include <iostream>
#include <vector>
#include "function.hpp" 
#include <math.h>

#define GSL_DBL_EPSILON        2.2204460492503131e-16

using namespace std;


/* From GSL */
/* roots/brent.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
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

/* brent.c -- brent root finding algorithm */
double brent(vector<double> &x, int max_iter, double eps, Func1d &f){
	
  double root;
  
  double a, b, c, d, e;
  double fa, fb, fc;

  root = 0.5 * (x[0] + x[1]) ;

  a = x[0];
  fa = f.eval(a);
  
  b = x[1];
  fb = f.eval(b);

  c = b;
  fc = fb;

  d = x[1] - x[0] ;
  e = d;

  if ((fa < 0.0 && fb < 0.0) || (fa > 0.0 && fb > 0.0))
  {
      cerr << "Brent::Endpoints do not straddle y=0" << endl;
      exit(1);
  }
  
  int iter = 0;
  
  double tolerance, old_val;
  bool more = true;
  
  double tol, m;

  do
  {
      
      int ac_equal = 0;

	if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
    {
      ac_equal = 1;
      c = a;
      fc = fa;
      d = b - a;
      e = b - a;
    }

  if (fabs (fc) < fabs (fb))
  {
      ac_equal = 1;
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
  }
  
  tol = 0.5 * GSL_DBL_EPSILON * fabs (b);
  m = 0.5 * (c - b);
  

  if (fb == 0)
    {
      root = b;
      x[0] = b;
      x[1] = b;
      
      return root;
    }
  
  if (fabs (m) <= tol)
    {
      root = b;

      if (b < c) 
        {
          x[0] = b;
          x[1] = c;
        }
      else
        {
          x[0] = c;
          x[1] = b;
        }

      return root;
    }
  
  if (fabs (e) < tol || fabs (fa) <= fabs (fb))
    {
      d = m;            /* use bisection */
      e = m;
    }
  else
    {

      double p, q, r;   /* use inverse cubic interpolation */
      double s = fb / fa;
      
      if (ac_equal)
        {
          p = 2 * m * s;
          q = 1 - s;
        }
      else
        {
          q = fa / fc;
          r = fb / fc;
          p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
          q = (q - 1) * (r - 1) * (s - 1);
        }
      
      if (p > 0)
        {
          q = -q;
        }
      else
        {
          p = -p;
        }
      
      if (2 * p < min (3 * m * q - fabs (tol * q), fabs (e * q)))
        {
          e = d;
          d = p / q;
        }
      else
        {
          /* interpolation failed, fall back to bisection */
          
          d = m;
          e = m;
        }
    }
  
  a = b;
  fa = fb;
  
  if (fabs (d) > tol)
    {
      b += d;
    }
  else
    {
      b += (m > 0 ? +tol : -tol);
    }
  
  fb = f.eval(b);

  
  /* Update the best estimate of the root and bounds on each
     iteration */
  
  root = b;
  
  double tc = c;
  if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) 
    {
      tc = a;
    }

  if (b < tc)
    {
      x[0] = b;
      x[1] = tc;
    }
  else
    {
      x[0] = tc;
      x[1] = b;
    }

    more = !(fabs(x[1] - x[0]) < eps * min( fabs(x[0]), fabs(x[1]) ));
    //cout << iter << " " << x[0] << " " << x[1] << endl;
    ++iter;
       
  } while (more && iter < max_iter);


  return root;
}
