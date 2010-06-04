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

void powell_new(double p[], double **xi, int n, double ftol, int *iter, 
	    double *fret,double (*func)(double []),int itmax)
{
  int i,ibig,j;
  double del,fp,fptt,t,*pt,*ptt,*xit;
  
  pt=dvector(0,n-1);
  ptt=dvector(0,n-1);
  xit=dvector(0,n-1);
  *fret=(*func)(p);
  for (j=0;j<n;j++) pt[j]=p[j];
  for (*iter=1;*iter<=itmax;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0;
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(p,xit,n,fret,func);
      if (fabs(fptt-(*fret)) > del) {
	del=fabs(fptt-(*fret));
	ibig=i;
      }
    }
    /* if we have converged break out of loop */
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) break;
    /* 
    if (*iter == itmax) 
      fprintf(stderr,"powell exceeding maximum iterations.\n");
      */
    for (j=0;j<n;j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(ptt);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if (t < 0.0) {
	linmin(p,xit,n,fret,func);
	for (j=0;j<n;j++) {
	  xi[j][ibig]=xi[j][n-1];
	  xi[j][n-1]=xit[j];
	}
      }
    }
  }
  free_dvector(xit,0,n-1);
  free_dvector(ptt,0,n-1);
  free_dvector(pt,0,n-1);
}

#define TOL 2.0e-4

static int ncom;
static double *pcom,*xicom,(*nrfunc)(double []);

void linmin_new(double p[], double xi[], int n, double *fret, 
	    double (*func)(double []))
{
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  
  ncom=n;
  pcom=dvector(0,n-1);
  xicom=dvector(0,n-1);
  nrfunc=func;
  for (j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
  for (j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_dvector(xicom,0,n-1);
  free_dvector(pcom,0,n-1);
}
#undef TOL

double f1dim_new(double x)
{
  int j;
  double f,*xt;

  xt=dvector(0,ncom-1);
  for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free_dvector(xt,0,ncom-1);
  return f;
}

#define EPS 1.0e-10

void frprmn_new(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []),
	    int itmax)
{
  int j,its;
  double gg,gam,fp,dgg;
  double *g,*h,*xi;
  
  g=dvector(0,n-1);
  h=dvector(0,n-1);
  xi=dvector(0,n-1);
  fp=(*func)(p);
  (*dfunc)(p,xi);
  for (j=0;j<n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=1;its<=itmax;its++) {
    *iter=its;
    linmin(p,xi,n,fret,func);
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      free_dvector(xi,0,n-1);
      free_dvector(h,0,n-1);
      free_dvector(g,0,n-1);
      return;
    }
    fp=(*func)(p);
    (*dfunc)(p,xi);
    dgg=gg=0.0;
    for (j=0;j<n;j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if (gg == 0.0) {
      free_dvector(xi,0,n-1);
      free_dvector(h,0,n-1);
      free_dvector(g,0,n-1);
      return;
    }
    gam=dgg/gg;
    for (j=0;j<n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }
  /* fprintf(stderr,"Too many iterations in frprmn\n"); */
}

#undef EPS

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

double amotry_new(double **p, double y[], double psum[], int ndim,
	      double (*funk)(double []), int ihi, double fac);

void amoeba_new(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk,int nmax)
{
  int i,ihi,ilo,inhi,j,mpts=ndim+1;
  double rtol,sum,swap,ysave,ytry,*psum;
  
  psum=dvector(0,ndim-1);
  *nfunk=0;
  for (j=0;j<ndim;j++) {
    for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];psum[j]=sum;
  }

  for (;*nfunk<=nmax;) {
    ilo=0;
    ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
    for (i=0;i<mpts;i++) {
      if (y[i] <= y[ilo]) ilo=i;
      if (y[i] > y[ihi]) {
	inhi=ihi;
	ihi=i;
      } else if (y[i] > y[inhi] && i != ihi) inhi=i;
    }
    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if (rtol < ftol) {
      SWAP(y[0],y[ilo]);
      for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i]);
      break;
    }
    /* if (*nfunk >= nmax) fprintf(stderr,"NMAX exceeded\n"); */
    *nfunk += 2;
    ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
    if (ytry <= y[ilo]) ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
    else if (ytry >= y[inhi]) {
      ysave=y[ihi];
      ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
      if (ytry >= ysave) {
	for (i=0;i<mpts;i++) {
	  if (i != ilo) {
	    for (j=0;j<ndim;j++) p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
	    y[i]=(*funk)(psum);
	  }
	}
	*nfunk += ndim;
	for (j=0;j<ndim;j++) {
	  for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];psum[j]=sum;}
      }
    } else --(*nfunk);
  }
  free_dvector(psum,0,ndim-1);
}
#undef SWAP

double amotry_new(double **p, double y[], double psum[], int ndim,
	      double (*funk)(double []), int ihi, double fac)
{
  int j;
  double fac1,fac2,ytry,*ptry;
  
  ptry=dvector(0,ndim-1);
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  ytry=(*funk)(ptry);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=0;j<ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  free_dvector(ptry,0,ndim-1);
  return ytry;
}
#undef NRANSI
