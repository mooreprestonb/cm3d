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


#include <stdio.h>
#include <math.h>

#include "minimize.h"

#define NRANSI


static double sqra=0;
#define SQR(A) (((sqra=(A)) == 0.)?0.:sqra*sqra)


/*----------------------------------------------------------------------*/

void powell(double p[], double **xi, int n, double ftol, int *iter, 
	    double *fret,double (*func)(double []),int itmax)
{

	powell_new(p, xi, n, ftol, iter, fret, func, itmax);

}

#define TOL 2.0e-4

static int ncom;
static double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, 
	    double (*func)(double []))
{

	linmin_new(p, xi, n, fret, func);

}
#undef TOL

double f1dim(double x)
{

	return(f1dim_new(x));

}

#define EPS 1.0e-10

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []),
	    int itmax)
{

	frprmn(p, n, ftol, iter, fret, func, dfunc, itmax);

}

#undef EPS

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

double amotry(double **p, double y[], double psum[], int ndim,
	      double (*funk)(double []), int ihi, double fac);

void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk,int nmax)
{

	amoeba_new(p, y, ndim, ftol, funk, nfunk, nmax);

}
#undef SWAP

double amotry(double **p, double y[], double psum[], int ndim,
	      double (*funk)(double []), int ihi, double fac)
{

	return(amotry_new(p, y, psum, ndim, funk, ihi, fac));

}
#undef NRANSI
