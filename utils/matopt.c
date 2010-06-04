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


/* 
   routine to perform matrix operations
   it calls libraries or fortran routines or does it itself 
   depending on what is defined
   */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef FORTRANUNDERSCORE
void bwdgemm_(int *,int *,int *,int *,int *,double *,void *,int *,
	      void *,int *,double *,void *,int *);
void bwdgemv_(int *,int *,int *,double *,double *,int *,double *,int *,
	      double *,double *,int *);
void bwdsytrf_(char *,int *,void *,int *,void *,void *,int *,int *);
void bwdsytri_(char *,int *,void *,int *,void *,void *,int *);
void bwdsytrs_(char *,int *,int *,void *,int *,void *,void *,int *,int *);
#else
void bwdgemm(int *,int *,int *,int *,int *,double *,void *,int *,
	      void *,int *,double *,void *,int *);
void bwdgemv(int *,int *,int *,double *,double *,int *,double *,int *,
	      double *,double *,int *);
void bwdsytrf(char *,int *,void *,int *,void *,void *,int *,int *);
void bwdsytri(char *,int *,void *,int *,void *,void *,int *);
void bwdsytrs(char *,int *,int *,void *,int *,void *,void *,int *,int *);
#endif


/*-----------------------------------------------------------------*/
void matmul(int n,double **a,double **b,double **c)
{
  int ita,itb,m,k,lda,ldb,ldc;
  double alpha,beta;

  m = k = n;
  lda = ldb = ldc = n;
  alpha = 1.;
  beta = 0.;

  ita = 0;
  itb = 0;

#ifdef BLAS
#ifdef CBLAS
  dgemm('N','N',m,n,k,alpha,&b[0][0],lda,&a[0][0],ldb,
	beta,&c[0][0],ldc);
#else
#ifdef FORTRANUNDERSCORE
  bwdgemm_(&ita,&itb,&m,&n,&k,&alpha,&b[0][0],&ldb,&a[0][0],&lda,
           &beta,&c[0][0],&ldc);
#else
  bwdgemm(&ita,&itb,&m,&n,&k,&alpha,&b[0][0],&ldb,&a[0][0],&lda,
	  &beta,&c[0][0],&ldc);
#endif /* FORTRANUNDERSCORE */
#endif /*CBLAS */
#else /* NO BLAS */
  for(ita=0;ita<n;ita++){
    for(itb=0;itb<n;itb++){
      beta = 0.;
      for(k=0;k<n;k++){
        beta += a[ita][k]*b[k][itb];
      }
      c[ita][itb] = beta;
    }
  }
#endif      /* BLAS */
}
/*-------------------------------------------------------------*/
void vecmat(int n,double **a,double *x,double *y)
{
  int ita,incx,incy,m,lda;
  double alpha,beta;

  ita = 0;
  incx = incy = 1;
  alpha = 1.;
  beta = 0.;
  m = lda = n;

#ifdef BLAS
#ifdef FORTRANUNDERSCORE
  bwdgemv_(&ita,&m,&n,&alpha,&a[0][0],&lda,&x[0],&incx,&beta,&y[0],&incy);
#else
  bwdgemv(&ita,&m,&n,&alpha,&a[0][0],&lda,&x[0],&incx,&beta,&y[0],&incy);
#endif
#else /* no blas */
  for(ita=0;ita<n;ita++){
    beta = 0.;
    for(m=0;m<n;m++){
      beta += a[m][ita]*x[m];
    }
    y[ita] = beta;
  }
#endif  
}
/*-------------------------------------------------------------*/
void invsymMatrix(int n,double **a)  /* invert symmetric matrix */
{
  char uplo;
  int lwork,lda,i,j;
  int *ipiv;
  double *work;
#ifdef BLAS
  int info;
#endif
  uplo = 'U';
  lwork = lda = n;

  ipiv = (int *)malloc(sizeof(int)*n);
  work = (double *)malloc(sizeof(double)*n);

#ifdef BLAS
#ifdef FORTRANUNDERSCORE /* have trouble w/ Preston's format */
  dsytrf_(&uplo,&n,&(a[0][0]),&lda,&(ipiv[0]),&(work[0]),&lwork,&info);
  if(info != 0){
    printf("dsytrf returned wrong value:%d\n",info);
    exit(2);
  }
  dsytri_(&uplo,&n,&(a[0][0]),&lda,&(ipiv[0]),&(work[0]),&info);
  if(info != 0){
    printf("dsytri returned wrong value\n");
    exit(2);
  }
#else
  bwdsytrf(&uplo,&n,&(a[0][0]),&lda,&(ipiv[0]),&(work[0]),&lwork,&info);
  if(info != 0){
    printf("dsytri returned wrong value\n");
    exit(2);
  }
  bwdsytri(&uplo,&n,&(a[0][0]),&lda,&(ipiv[0]),&(work[0]),&info);
  if(info != 0){
    printf("dsytri returned wrong value\n");
    exit(2);
  }
#endif
#else /* no blas */
  fprintf(stderr,"ERROR, invsymMatrix must be compilied with BLAS\n");
  exit(1);
#endif

  for(i=0;i<n;i++)/* fill in other half */
    for(j=0;j<n;j++) a[i][j]=a[j][i];

  
  free(ipiv);
  free(work);
  
  return;
}
/*-------------------------------------------------------------*/
/* solve a system of linear equations
   A*X = B with a real symmetric matrix A
   using the factorization A = U*D*U**T or
   A = L*D*L**T computed by DSYTRF  */
void sol4musymMat(int n,double **a,double *b)
{
  char uplo;
  int lwork,lda,ldb,nrhs;
  int *ipiv;
  double *work;
#ifdef BLAS
  int info,i,j;
#endif

  uplo = 'U';
  lwork = lda = ldb = n;
  nrhs = 1;

  ipiv = (int *)malloc(sizeof(int)*n);
  work = (double *)malloc(sizeof(double)*n);

  
#ifdef BLAS
#ifdef FORTRANUNDERSCORE /* have trouble w/ Preston's format */
  dsytrf_(&uplo,&n,&(a[0][0]),&lda,&(ipiv[0]),&(work[0]),&lwork,&info);
  if(info != 0){
    printf("dsytrf returned wrong value:%d\n",info);
    exit(2);
  }
  printf("Inside sol4musymMat \n");
  for(i=0;i<n/3;i+=3)
    for(j=0;j<n/3;j++)printf("%g %g %g\n",a[i][j],a[i+1][j],a[i+2][j]);
  printf("\n");

  dsytrs_(&uplo,&n,&nrhs,&(a[0][0]),&lda,&(ipiv[0]),&(b[0]),&ldb,&info);
  if(info != 0){
    printf("dsytrs returned wrong value\n");
    exit(2);
  }
#else
  bwdsytrf(&uplo,&n,&(a[0][0]),&lda,&(ipiv[0]),&(work[0]),&lwork,&info);
  if(info != 0){
    printf("dsytrf returned wrong value\n");
    exit(2);
  }
  bwdsytrs(&uplo,&n,&nrhs,&(a[0][0]),&lda,&(ipiv[0]),&(b[0]),&ldb,&info);
  if(info != 0){
    printf("dsytrs returned wrong value\n");
    exit(2);
  }
#endif
#else /* no blas */
  fprintf(stderr,"ERROR, sol4musymMat must be compilied with BLAS\n");
  exit(1);
#endif

  free(ipiv);
  free(work);
  
  return;
}
/*-------------------------------------------------------------*/



