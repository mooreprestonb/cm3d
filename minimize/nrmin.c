/* 
   Copyright (C) 1997 1998 1999 2000 Dr. Preston B. Moore

Dr. Preston B. Moore
Associate Director, Center for Molecular Modeling (CMM)
University of Pennsylvania, Department of Chemistry, Box 188 
231 S. 34th St. Philadelphia, PA 19104-6323 USA
EMAIL: moore@cmm.chem.upenn.edu  
WWW: http://www.cmm.upenn.edu/~moore

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#include <math.h>
#include <stdio.h>

#include "minimize.h"

#define NRANSI
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static double maxa1,maxa2;
#define MAX(a,b) (maxa1=(a),maxa2=(b),((maxa1)>(maxa2))?(maxa1):(maxa2))
#define SIGN(a,b) ((b)>0.0 ? fabs(a) : -fabs(b))

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, 
	    double *fc, double (*func)(double))
{

	mnbrak_new(ax, bx, cx, fa, fb, fc, func);

}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI


#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double golden(double ax, double bx, double cx, double (*f)(double), double tol,
	      double *xmin)
{

	return(golden_new(ax, bx, cx, f, tol, xmin));

}
#undef C
#undef R
#undef SHFT2
#undef SHFT3

#define NRANSI
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	     double *xmin)
{

	return(brent_new(ax, bx, cx, f, tol, xmin));


}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI

#define NRANSI
#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double dbrent(double ax, double bx, double cx, double (*f)(double),
	      double (*df)(double), double tol, double *xmin)
{

	return(dbrent_new(ax, bx, cx, f, df, tol, xmin));

}
#undef ITMAX
#undef ZEPS
#undef MOV3
#undef NRANSI
