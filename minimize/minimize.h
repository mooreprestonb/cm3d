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

void linmin(double p[], double xi[], int n, double *fret,
	    double (*func)(double []));
double brent(double ax, double bx, double cx,
	     double (*f)(double), double tol, double *xmin);
double f1dim(double x);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
	    double *fc, double (*func)(double));

void powell(double p[], double **xi, int n, double ftol, int *iter, 
	    double *fret,double (*func)(double []),int);

double amotry(double **p, double y[], double psum[], int ndim,
	      double (*funk)(double []), int ihi, double fac);
void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk,int);
void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []),int);



void linmin_new(double p[], double xi[], int n, double *fret,
	    double (*func)(double []));
double brent_new(double ax, double bx, double cx,
	     double (*f)(double), double tol, double *xmin);
double f1dim_new(double x);
void mnbrak_new(double *ax, double *bx, double *cx, double *fa, double *fb,
	    double *fc, double (*func)(double));

void powell_new(double p[], double **xi, int n, double ftol, int *iter, 
	    double *fret,double (*func)(double []),int);

double amotry_new(double **p, double y[], double psum[], int ndim,
	      double (*funk)(double []), int ihi, double fac);
void amoeba_new(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk,int);
void frprmn_new(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []),int);



void mal_verify_new(int);

double *dvector(int nl,int nh);
void free_dvector(double *v,int nl,int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);


