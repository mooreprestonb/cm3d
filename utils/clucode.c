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

/* routines for solving linear algebra */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void invres(int,double **,double **);
void luerr(char *);
void lubksb(double **,int,int *,double []);
void ludcmp(double **,int,int *,double *);
void savgol(double [],int,int,int,int,int);
double **dmatrix(int,int,int,int);
void free_dmatrix(double **,int,int,int,int);

#ifndef MIN
#define MIN(A,B) ((A)<(B) ? (A) :(B))
#endif

#ifndef TINY
#define TINY 1.0e-20;
#endif

/* invert a matrix using ludcmp and lubksb */
void invres(int n,double **a,double **ai)
{
  int i,j,*indx;
  double *col,d;

  col  = (double *)malloc(n*sizeof(double));
  indx = (int *)malloc(n*sizeof(int));

  ludcmp(a,n,indx,&d);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++) col[j] = 0.;
    col[i] = 1.;
    lubksb(a,n,indx,col);
    for(j=0;j<n;j++) ai[j][i] = col[j];
  }

  free(col);free(indx);
}
/*------------------------------------------------------*/
void luerr(char *s){fprintf(stderr,"ERROR: %s\n",s); exit(1);}

/*-------------------------------------------------------*/
/* lower upper decomposition */
void ludcmp(double **a,int n,int *indx,double *d)
{
  int i,j,k,imax=0;
  double big,dum,sum,temp,*vv;

  vv=(double *)malloc(n*sizeof(double));
  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) luerr("Singular matrix in routine LUDCMP");
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}

/*--------------------------------------------------------*/

void lubksb(double **a,int n,int *indx,double b[])
{
  int i,j,ii=-1,ip;
  double sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii != -1) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }

  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

/*--------------------------------------------------------*/
/* savgol smoothing coeffecient function */

void savgol(double c[],int np,int nl,int nr,int ld,int m)
{
  int imj,ipj,j,k,kk,mm,*indx;
  double d,fac,sum,**a,*b;
  
  if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m){
    luerr("bad args in savgol");
  }

  indx=(int *)malloc((m+1)*sizeof(int));
  a=dmatrix(0,m+1,0,m+1);
  b = (double *)malloc((m+1)*sizeof(double));

  for (ipj=0;ipj<=(m << 1);ipj++) {
    sum=(ipj ? 0.0 : 1.0);
    for (k=1;k<=nr;k++) sum += pow((double)k,(double)ipj);
    for (k=1;k<=nl;k++) sum += pow((double)-k,(double)ipj);
    mm=MIN(ipj,2*m-ipj);
    for (imj = -mm;imj<=mm;imj+=2) a[(ipj+imj)/2][(ipj-imj)/2]=sum;
  }

  ludcmp(a,m+1,indx,&d);
  for (j=0;j<m+1;j++) b[j]=0.0;
  b[ld]=1.0;
  lubksb(a,m+1,indx,b);

  for (kk=0;kk<np;kk++) c[kk]=0.0;
  for (k = -nl;k<=nr;k++) {
    sum=b[0];
    fac=1.0;
    for (mm=0;mm<m;mm++) sum += b[mm+1]*(fac *= k);
    kk=((np-k)%np);
    c[kk]=sum;
  }
  free(b);  free(indx);  free_dmatrix(a,0,m+1,0,m+1);
}

