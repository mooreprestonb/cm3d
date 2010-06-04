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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.41421356237309504880
#endif

/* #define SIMPLE */
double randme(void);
double gauss1(void);
void ranv(int,double *);
double dsum(int,double *,int);
void ggauss(int,double *);
int locate(int *,int,int);
int located(double *,int,double);
int insert(int *ia,int *n,int jpoint);
int exclij(int i,int j,int *nexcl,int **excl);
void splint(double [], double [], double [], int , double , double *);
void spline(double [], double [], int , double , double , double []);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

/*--------------------------------------------------------*/
void ggauss(int nran,double *gauss)
{
  int i;
  double phi,al;
  double *unif;

  unif = (double *)malloc(nran*sizeof(double));
  ranv(nran,unif);

#ifdef SIMPLE
  for(i=0;i<nran;i++){
    phi = 2.*M_PI*unif[i];
    al = sqrt(-log(unif[i]))*M_SQRT2;
    gauss[i] = al*cos(phi);
  }
#else
  if(nran%2 == 1) gauss[nran-1] = gauss1();
#ifdef SGI
#pragma ivdep
#endif
#ifdef CRAY
#pragma _CRI ivdep
#endif
  for(i=0;i<nran-1;i += 2){
    phi = 2.*M_PI*unif[i];  
    al = sqrt(-log(unif[i]))*M_SQRT2;
    gauss[i] = al*cos(phi); 
    gauss[i+1] = al*sin(phi);
  }
#endif

  free(unif);
}
/*--------------------------------------------------------*/
double gauss1(void)
{
  double phi,al;
  phi = 2.*M_PI*randme();
  al = sqrt(-log(randme()))*M_SQRT2;
  return (al*cos(phi));
}
/*--------------------------------------------------------*/
void ranv(int nran,double *unif)
{
  int i;

  for(i=0;i<nran;i++) unif[i] = randme();
}
/*--------------------------------------------------------*/
double dsum(int n,double *v,int inc)
{
  int i;
  double sum;

  sum = 0.0; for(i=0;i<n;i += inc) sum += v[i];
  return sum;
}
/*--------------------------------------------------------*/
/* search an order integer list j */
int locate(int *ia,int n,int jpoint)
{
  int ju,jm,jl;

  jl = 0;ju = n;
  if(n!=0 && jpoint>ia[n-1]) { return n;}
  while(ju-jl>1){
    jm=(ju+jl)>>1; /* midpoint */
    ((jpoint >= ia[jm])? (jl=jm): (ju=jm));
  }
  return jl;
}
/*--------------------------------------------------------*/
/* search an order double list j */
int located(double *a,int n,double point)
{
  int ju,jm,jl;

  jl = 0;ju = n;
  if(n!=0 && point>a[n-1]) { return n;}
  while(ju-jl>1){
    jm=(ju+jl)>>1; /* midpoint */
    ((point >= a[jm])? (jl=jm): (ju=jm));
  }
  return jl;
}
/*--------------------------------------------------------*/
int insert(int *ia,int *n,int jpoint)
{
  int i,np;

  np = locate(ia,*n,jpoint); /* locate midpoint */

  if(np == *n){              /* take care if first/last interaction */
    ia[np] = jpoint;
    (*n)++;
    return 1;
  }

  if(ia[np] != jpoint){     /* if not in list insert */
    if(ia[np] < jpoint) np++;  /* get correct np place */
    for(i=*n;i>np;i--){        /* move table down */
      ia[i] = ia[i-1];
    }
    ia[np] = jpoint;           /* insert */
    (*n)++;
    return 1;
  }
  return 0;
}
/*-----------------------------------------------------------------------*/
int exclij(int i,int j,int *nexcl,int **excl)
{
  int k;
  if(i>j){k=i;i=j;j=k;}

  k = locate(excl[i],nexcl[i],j);

  if(excl[i][k]==j) return 1;
  else return 0;
}
/*------------------------------------------------------------------------*/
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
  int klo,khi,k;
  double h,b,a;

  /* find point via bisection */
  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
#ifdef DEBUG
  if (h == 0.0) {
    fprintf(stderr,"Bad xa input to routine splint");
    exit(1);
  }
#endif
  /* calculate spline */
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
/*------------------------------------------------------------------------*/
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;

  u=(double *)malloc(n*sizeof(double));
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}
/*------------------------------------------------------------------------*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;

  dif=fabs(x-xa[1]);
  c=(double *)malloc((n+1)*sizeof(double));
  d=(double *)malloc((n+1)*sizeof(double));
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0){
	fprintf(stderr,"Error in routine polint");
	exit(1);
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free(d);
  free(c);
}
