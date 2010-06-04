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


/* converted from numerical reciepies so that matrix 
   goes from 0 to n-1 NOT 1 to n */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

  /* #define DEBUG  out error when TQLI fails */

/* #define BLAS */
/* #define CRAY */

void ctred2(double **,int,double [],double []);
void ctred2v(double **,int,double [],double []);
void ctqli(double [],double [],int,double **);
int ctqliv(double [],double [],int,double **);
void ceigsrt(double [],double **,int);
void ceigsrtv(double [],double **,int);
extern void md_error(char *);
extern double **dmatrix(int,int,int,int);
extern void free_dmatrix(double **,int,int,int,int);

#ifdef ESSL
#define NEED_Z
#endif /* ESSL */

#ifdef CRAY
extern void RS(int*,int*,double*,double*,int*,double*,double*,double*,int*);
#define NEED_Z
#else
extern void lwdsyev_(int*,int*,int*,double*,int*,double*,double*,int*,int*);
#endif /* CRAY */

void rs_me(int n,double *eigval,double **amat,int iflag)
{
  int ijz,iul,lda,info,lwork;
  double *work;
#ifdef DEBUG
  double **tmat;
#endif
#ifdef ESSL
  double tmp;
#endif
#ifdef NEED_Z
  int i;
  double *z,*fv1,*fv2;

  fv1 = malloc(2*n*sizeof(double));
  fv2 = fv1+n;
  if(fv1 ==NULL){
    md_error("can't allocate work array (RS eigenvector)");
  }
  if(iflag!=0){
    if((z = malloc(n*n*sizeof(double)))==NULL){
      md_error("can't allocate extra matrix (RS eigenvector)\n");
      exit(1);
    }
  }
#endif /* NEED_Z */

  ijz = 0;
  lda = n;
  lwork = 3*n;

  if((work = malloc(lwork*sizeof(double)))==NULL){
    md_error("can't allocate work array for eigenvector solve\n");
    exit(1);
  }

  switch(iflag){
  case 0:  /* eigenvalues only */
    info = ijz = iul = 0;
    break;
  case 1:  /* eigenvalues AND eigenvectors */
    info = iul = 0;
    ijz = 1;
    break;
  }

#ifdef ESSL
    dspev(ijz,&amat[0][0],&eigval[0],z,n,n,&work[0],lwork);
    for(i=0;i<n*n;i++) amat[0][i] = z[i];
    /* transpose so that collum are the vectors instead of rows */
    for(ijz = 0;ijz<n-1;ijz++){
      for(iul = ijz+1;iul<n;iul++){
	tmp = amat[ijz][iul];
	amat[ijz][iul] = amat[iul][ijz];
	amat[iul][ijz] = tmp;
      }
    }
#else /* not ESSL */
#ifdef BLAS
#ifdef CRAY
    RS(&n,&n,&amat[0][0],&eigval[0],&iflag,&z[0],&fv1[0],&fv2[0],&info);
    /* we have to copy z back to a (use contiguous memory) */
    for(i=0;i<n*n;i++) amat[0][i] = z[i];
#else /* NOT CRAY */
#ifdef FORTRANUNDERSCORE
    lwdsyev_(&ijz,&iul,&n,&amat[0][0],&lda,&eigval[0],&work[0],&lwork,&info);
#else
    lwdsyev(&ijz,&iul,&n,&amat[0][0],&lda,&eigval[0],&work[0],&lwork,&info);
#endif /* FORTRANUNDERSCORE */
#endif /* CRAY */
    /* transpose so that collum are the vectors instead of rows */
    for(ijz = 0;ijz<n-1;ijz++){
      int tmp;
      for(iul = ijz+1;iul<n;iul++){
	tmp = amat[ijz][iul];
	amat[ijz][iul] = amat[iul][ijz];
	amat[iul][ijz] = tmp;
      }
    }
#else /* NOT BLAS */
    if(ijz==1){
#ifdef DEBUG
      tmat = dmatrix(0,n-1,0,n-1);
      for(ijz=0;ijz<n;ijz++) 
	for(iul=0;iul<n;iul++) 
	  tmat[ijz][iul]=amat[ijz][iul];
#endif
      ctred2v(amat,n,eigval,work);
      if(ctqliv(eigval,work,n,amat)){
	fprintf(stderr,"Too many iterations in NR's TQLI\n");
#ifdef DEBUG
	printf("n = %d\n",n);
	for(ijz=0;ijz<n;ijz++){
	  printf("\n%d ",ijz);
	  for(iul=0;iul<n;iul++){
	    printf("%g ",tmat[ijz][iul]);
	  }
	}
	printf("\n");
#endif
	exit(1);
      }
#ifdef DEBUG
      free_dmatrix(tmat,0,n-1,0,n-1);
#endif
    } else {
      ctred2(amat,n,eigval,work);
      ctqli(eigval,work,n,amat);
    }
#endif /* BLAS */
#endif /* ESSL */

  free(work);

#ifdef NEED_Z
  free(fv1);
  if(iflag!=0) free(z);
#endif

}

/*----------------------------------------------------------------*/
#define VECTOR
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

void ctred2(double **a,int n,double d[],double e[])
{
  int l,k,j,i;
  double scale,hh,h,g,f;
  
  for(i=n-1;i>=1;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=0;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = f>0 ? -sqrt(h) : sqrt(h);
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=0;j<=l;j++) {
	  g=0.0;
	  for (k=0;k<=j;k++) g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++) g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=0;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=0;k<=j;k++) a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  e[0]=0.0;
  for (i=0;i<n;i++) {
    d[i]=a[i][i];
  }
}

/*----------------------------------------------------------------*/
void ctred2v(double **a,int n,double d[],double e[])
{
  int l,k,j,i;
  double scale,hh,h,g,f;
  
  for(i=n-1;i>=1;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=0;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = f>0 ? -sqrt(h) : sqrt(h);
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=0;j<=l;j++) {
	  /* Next statement can be omitted if eigenvectors not wanted */
#ifdef VECTOR
	  a[j][i]=a[i][j]/h;
#endif
	  g=0.0;
	  for (k=0;k<=j;k++) g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++) g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=0;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=0;k<=j;k++) a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  /* Next statement can be omitted if eigenvectors not wanted */
#ifdef VECTOR
  d[0]=0.0;
#endif
  e[0]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i=0;i<n;i++) {
#ifdef VECTOR
    l=i-1;
    if (d[i]) {
      for (j=0;j<=l;j++) {
	g=0.0;
	for (k=0;k<=l;k++) g += a[i][k]*a[k][j];
	for (k=0;k<=l;k++) a[k][j] -= g*a[k][i];
      }
    }
#endif
    d[i]=a[i][i];
#ifdef VECTOR
    a[i][i]=1.0;
    for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0;
#endif
  }
}

/*----------------------------------------------------------------*/
/* without the vectors */
void ctqli(double d[],double e[],int n,double **z)
{
  int m,l,iter,i;
  double s,r,p,g,f,dd,c,b;
  
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-2;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if (fabs(e[m])+dd == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  fprintf(stderr,"Too many iterations in NR's TQLI\n");
	  exit(1);
	}
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=sqrt((g*g)+1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  if (fabs(f) >= fabs(g)) {
	    c=g/f;
	    r=sqrt((c*c)+1.0);
	    e[i+1]=f*r;
	    c *= (s=1.0/r);
	  } else {
	    s=f/g;
	    r=sqrt((s*s)+1.0);
	    e[i+1]=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  p=s*r;
	  d[i+1]=g+p;
	  g=c*r-b;
	}
	d[l]=d[l]-p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}

/*----------------------------------------------------------------*/

int ctqliv(double d[],double e[],int n,double **z)
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-2;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if (fabs(e[m])+dd == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  return(1);
	}
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=sqrt((g*g)+1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  if (fabs(f) >= fabs(g)) {
	    c=g/f;
	    r=sqrt((c*c)+1.0);
	    e[i+1]=f*r;
	    c *= (s=1.0/r);
	  } else {
	    s=f/g;
	    r=sqrt((s*s)+1.0);
	    e[i+1]=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  p=s*r;
	  d[i+1]=g+p;
	  g=c*r-b;
	  /* Next loop can be omitted if eigenvectors not wanted */
#ifdef VECTOR
	  for (k=0;k<n;k++) {
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
#endif
	}
	d[l]=d[l]-p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
  return 0;
}

/*-------------------------------------------------------------*/

void ceigsrt(double d[],double **v,int n)
{
  int k,j,i;
  double p;
  
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++) if (d[j] >= p) p=d[k=j];
    if (k != i) {d[k]=d[i]; d[i]=p;}
  }
}
/*-------------------------------------------------------------*/

void ceigsrtv(double d[],double **v,int n)
{
  int k,j,i;
  double p;
  
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++) if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];  d[i]=p;
#ifdef VECTOR
      for (j=0;j<n;j++) {p=v[j][i];v[j][i]=v[j][k];v[j][k]=p;}
#endif
    }
  }
}
/*-------------------------------------------------------------*/

void vsrtasnd(double d[],double **v,int m,int n)
{
  int k,j,i;
  double p;
  
  for (i=0;i<m-1;i++) {
    p=d[k=i];
    for (j=i+1;j<m;j++) if (d[j] <= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];  d[i]=p;
      for (j=0;j<n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;}
    }
  }
}
/*-------------------------------------------------------------------*/
