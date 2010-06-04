/*   Copyright (C) 1997 1998 1999 2000 Dr. Preston B. Moore

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

#include "md.h"
/*
#define PRINT_DEFIELD
#define PRINT_TRIATIC
#define PRINT_T_DIPOLE
*/
/* #define HDEBUG */
/* #define EWALD *//* this doesn't appear to be working properly */
#define GETDIPDER
/* #define PRINT_DIPOLE */

void getdid(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **dmumat)
{
  int i,j;
  int nef,n3,n;
  int *ipiv;
  double a3;
  double *mu,*bvec,**amat,**mat2;

#ifdef BLAS
  char uplo;
  int lwork,lda,ldb,info,nrhs;
  double *work;

  n = simparms->natoms;
  n3 = 3*n;
  uplo = 'U';
  lwork = lda = ldb = n3;
  nrhs = 1;

  work = (double *)malloc(sizeof(double)*n3);
#endif
  
#ifdef EWALD
  nef = 0;
#else
  nef = 1;
#endif
  
  n = simparms->natoms;
  n3 = 3*n;
  
  mu = cmalloc(n3*sizeof(double));
  ipiv = cmalloc(sizeof(int)*n3);
  bvec = cmalloc(n3*sizeof(double));
  amat = dmatrix(0,n3-1,0,n3-1);
  mat2 = dmatrix(0,n3-1,0,n3-1);

  /* b vector and A matrix */
  /* set up the equation M x = b 
     M is the matrix with Tij's on the off diagonal block 
       and a^(-1) along the diagonal blocks
     x is the dipole that we want
     b is the electric field
  */
  getamatrix(simparms,coords,amat,nef,inter,bvec);

  /* get the Electric Field from Ewald sums */
  if(nef==0)getewef(simparms,coords,bvec);
  
/*   printf("Electric Field\n"); */
/*   for(i=0;i<n3;i++) printf("%g ",bvec[i]); */
/*   printf("\n"); */

  /* save electric field */
  for(i=0;i<n3;i++) mu[i] = bvec[i];

  /* now that the system is set up 
     solve the thing for the dipoles */
#ifdef BLAS
#ifdef FORTRANUNDERSCORE
  dsytrf_(&uplo,&n3,&(amat[0][0]),&lda,&(ipiv[0]),&(work[0]),&lwork,&info);
  if(info != 0){
    printf("dsytrf returned wrong value:%d\n",info);
    exit(2);
  }
  dsytrs_(&uplo,&n3,&nrhs,&(amat[0][0]),&lda,&(ipiv[0]),&(mu[0]),&ldb,&info);
  if(info != 0){
    printf("dsytrs returned wrong value\n");
    exit(2);
  }
#else
  bwdsytrf(&uplo,&n,&(amat[0][0]),&lda,&(ipiv[0]),&(work[0]),&lwork,&info);
  if(info != 0){
    printf("dsytrf returned wrong value\n");
    exit(2);
  }
  bwdsytrs(&uplo,&n,&nrhs,&(amat[0][0]),&lda,&(ipiv[0]),&(mu[0]),&ldb,&info);
  if(info != 0){
    printf("dsytrs returned wrong value\n");
    exit(2);
  }
#endif /* FORTRANUNDERSCORE*/
#else /* NO BLAS */
  ludcmp(amat,n3,ipiv,&a3);
  lubksb(amat,n3,ipiv,mu);
#endif
  
#ifdef PRINT_DIPOLE
  printf("solution\n");
  a3 = 0.;
  for(i=0;i<n3;i++){
    a3 += mu[i]*bvec[i];
    printf("%5d %10g %10g\n",i,mu[i],bvec[i]);
  }
  printf("V = %g\n",-.5*a3);
  exit(1);
#endif

  /* save dipole information */
  for(j=0,i=0;i<n;i++,j+=3){
    coords->ux[i] = mu[j];
    coords->uy[i] = mu[j+1];
    coords->uz[i] = mu[j+2];
  }

#ifdef GETDIPDER
  /* get the gradient of the dipole */
  /* set up matrix */
  for(i=0;i<n3;i++) for(j=0;j<n3;j++) mat2[i][j] = 0.0;  
  getdefield(simparms,coords,inter,mat2);
  getddid(simparms,coords,mat2);
  
  /* use the same inverse matrix */  
  for(i=0;i<n3;i++){
    /* load the vector */
    for(j=0;j<n3;j++) mu[j] = mat2[j][i];
#ifdef BLAS
#ifdef FORTRANUNDERSCORE
    dsytrs_(&uplo,&n3,&nrhs,&(amat[0][0]),&lda,&(ipiv[0]),&(mu[0]),&ldb,&info);
    if(info != 0){
      printf("dsytrs returned wrong value\n");
      exit(2);
    }
#else
    bwdsytrs(&uplo,&n,&nrhs,&(amat[0][0]),&lda,&(ipiv[0]),&(mu[0]),&ldb,&info);
    if(info != 0){
      printf("dsytrs returned wrong value\n");
      exit(2);
    }
#endif /* FORTRANUNDERSCORE */
#else
    lubksb(amat,n3,ipiv,mu);
#endif /* BLAS */
    /* unpack to mat2 */
    for(j=0;j<n3;j++) dmumat[j][i] = mu[j];    
  }
#endif /* GETDIPDER */
  
  free(ipiv); free(bvec); free(mu);

  free_dmatrix(amat,0,n3-1,0,n3-1);
  free_dmatrix(mat2,0,n3-1,0,n3-1);
}
/*----------------------------------------------------------------*/
void getddid(SIMPARMS *simparms,COORDS *coords,double **amat)
{
  int i,j,k,l,m,ii,jj,n,iperd,iensemble;
  double w[3];
  double tmat[3][3][3];
  double r,r2,r3,r5,r7,hmati[9];
  double *x,*y,*z,*hmat;
  double *ux,*uy,*uz;

  iperd = simparms->iperd;
  iensemble = simparms->iensemble;
  n = simparms->natoms;
  x = coords->px;
  y = coords->py;
  z = coords->pz;
  hmat = coords->hmat;
  ux = coords->ux;
  uy = coords->uy;
  uz = coords->uz;

  gethinv9(hmat,hmati);

  for(i=0;i<n;i++){
    ii = 3*i;
    for(j=0;j<n;j++){
      if(i!=j){
        jj = 3*j;
        w[0] = x[i] - x[j];
        w[1] = y[i] - y[j];
        w[2] = z[i] - z[j];
        period(1,&w[0],&w[1],&w[2],hmat,hmati,iperd,iensemble);
        r2 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
        r = sqrt(r2);
        r3 = 1./(r2*r);
        r5 = r3/r2;
        r7 = r5/r2;
        
        /* get triatic */
        for(m=0;m<3;m++){
          for(k=0;k<3;k++){
            for(l=0;l<3;l++){
              tmat[k][l][m] = 15.*r7*w[k]*w[l]*w[m];
              if(l==k) tmat[k][l][m] += -3.*w[m]*r5;
              if(l==m) tmat[k][l][m] += -3.*w[k]*r5;
              if(m==k) tmat[k][l][m] += -3.*w[l]*r5;
            }
          }
        }
        
#ifdef PRINT_TRIATIC
        printf("dipole = %d-%d %g %g %g %g %g %g\n",
           i,j,ux[i],uy[i],uz[i],ux[j],uy[j],uz[j]);
        for(k=0;k<3;k++) {
          for(l=0;l<3;l++){
            for(m=0;m<2;m++){
              printf("%g ",tmat[k][l][m]);
            }
            printf("%g\t",tmat[k][l][m]);
          }
          printf("\n");
        }
#endif
        /* i-j particle dipole */
        for(k=0;k<3;k++){
          for(l=0;l<3;l++){
            amat[ii+k][jj+l] += (ux[j]*tmat[k][0][l] +
               uy[j]*tmat[k][1][l] +
               uz[j]*tmat[k][2][l]);
          }
        }
        /* i-i particle dipole (negative of the off diagonal) */
        for(k=0;k<3;k++){
          for(l=0;l<3;l++){
            amat[ii+k][ii+l] -= (ux[j]*tmat[k][0][l] +
               uy[j]*tmat[k][1][l] +
               uz[j]*tmat[k][2][l]);
          }
        }
      } /* end if (i!=j) */
    }
  }
#ifdef PRINT_T_DIPOLE
  printf("Triatic\n");
  for(i=0;i<n*3;i++) {
    for(j=0;j<n*3;j++){
      printf("%g ",amat[i][j]);
    }
    printf("\n");
  }
  printf("\n");
#endif
}

/*----------------------------------------------------------------*/
void getdiatic(SIMPARMS *simparms,COORDS *coords,
	       int i_now,int j_now,double **diatic)
{
  int i,j,k,l,iperd,iensemble;
  double w[3],amat[3][3];
  double r,r2,r3,r5,hmati[9];
  double *x,*y,*z,*hmat;

  x = coords->px;
  y = coords->py;
  z = coords->pz;
  hmat = coords->hmat;
  iperd = simparms->iperd;
  iensemble = simparms->iensemble;
  gethinv9(hmat,hmati);

  for(i=0;i<3;i++) for(j=0;j<3;j++) amat[i][j] = 0.0;
  
  w[0] = x[i_now] - x[j_now];
  w[1] = y[i_now] - y[j_now];
  w[2] = z[i_now] - z[j_now];
  period(1,&w[0],&w[1],&w[2],hmat,hmati,iperd,iensemble);
  r2 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
  r = sqrt(r2);
  r3 = 1./(r2*r);
  r5 = r3/r2;
      
  /* get diatic */
  for(k=0;k<3;k++){
    for(l=0;l<3;l++){
      amat[k][l] = -3.*r5*w[k]*w[l];
      if(l==k) amat[k][l] += r3;
    }
  }

#ifdef PRINT_DIATIC
  for(i=0;i<3;i++) {
    for(j=0;j<2;j++){
      printf("%g ",amat[i][j]);
    }
    printf("%g\n",amat[i][j]);
  }
  printf("\n");
#endif

  /* get diatic */
  for(k=0;k<3;k++){
    for(l=0;l<3;l++){
      diatic[k][l] = amat[k][l];
    }
  }
}
/*----------------------------------------------------------------*/
void test_diatic(SIMPARMS *simparms,COORDS *coords)
{
  int i,j,n3;
  double *mat1,*mat2,**amat;
  double h,xo;

  n3 = simparms->natoms*3;
  h = 1.e-6;
  mat1 = malloc(n3*sizeof(double));
  mat2 = malloc(n3*sizeof(double));
  amat = dmatrix(0,n3-1,0,n3-1);
  for(i=0;i<n3;i++) for(j=0;j<n3;j++) amat[i][j] = 0.0;
  
  for(i=0;i<simparms->natoms;i++){
    xo = coords->px[i];
    coords->px[i] = xo - h; 
    getmu(simparms,coords,mat1);
    coords->px[i] = xo + h; 
    getmu(simparms,coords,mat2);
    coords->px[i] = xo; 
    
    printf("%d dx  ",i);
    for(j=0;j<n3;j++){
      xo = (mat2[j]-mat1[j])/(2*h);
      amat[i*3][j] = xo;
      printf("%g ",xo);
    }
    printf("\n");

    xo = coords->py[i];
    coords->py[i] = xo - h; 
    getmu(simparms,coords,mat1);
    coords->py[i] = xo + h; 
    getmu(simparms,coords,mat2);
    coords->py[i] = xo; 
    
    printf("%d dy  ",i);
    for(j=0;j<n3;j++){
      xo = (mat2[j]-mat1[j])/(2*h);
      amat[i*3+1][j] = xo;
      printf("%g ",xo);
    }
    printf("\n");

    xo = coords->pz[i];
    coords->pz[i] = xo - h; 
    getmu(simparms,coords,mat1);
    coords->pz[i] = xo + h; 
    getmu(simparms,coords,mat2);
    coords->pz[i] = xo; 
    
    printf("%d dz  ",i);
    for(j=0;j<n3;j++){
      xo = (mat2[j]-mat1[j])/(2*h);
      amat[i*3+2][j] = xo;
      printf("%g ",xo);
    }
    printf("\n");
  }

  exit(1);

  free(mat2); free(mat1);
  free_dmatrix(amat,0,n3-1,0,n3-1);
}  

/*----------------------------------------------------------------*/
void getefield(SIMPARMS *simparms,COORDS *coords,double *bfield)
{
  int i,j,ii,jj,n,iperd,iensemble;
  double a3,w[3];
  double r,r2,r3,hmati[9];
  double *x,*y,*z,*qch,*hmat;

  x = coords->px;
  y = coords->py;
  z = coords->pz;
  qch = coords->qch;
  hmat = coords->hmat;
  iperd = simparms->iperd;
  iensemble = simparms->iensemble;
  n = simparms->natoms;
  gethinv9(hmat,hmati);

  /* b is the electric field */

  for(i=0;i<3*n;i++) bfield[i] = 0.0;

  for(i=0;i<n-1;i++){
    ii = 3*i;
    for(j=i+1;j<n;j++){
      w[0] = x[i] - x[j];
      w[1] = y[i] - y[j];
      w[2] = z[i] - z[j];
      period(1,&w[0],&w[1],&w[2],hmat,hmati,iperd,iensemble);
      r2 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
      r = sqrt(r2);
      r3 = 1./(r2*r);
      
      a3 = r3*qch[j];
      bfield[ii  ] += w[0]*a3;
      bfield[ii+1] += w[1]*a3;
      bfield[ii+2] += w[2]*a3;
      
      jj = 3*j;
      a3 = r3*qch[i];
      bfield[jj  ] -= w[0]*a3;
      bfield[jj+1] -= w[1]*a3;
      bfield[jj+2] -= w[2]*a3;
    }
  }
}
/*---------------------------------------------------------------*/
void getdefield(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **amat)
{
  int i,j,k,l,ii,jj,n,iperd,iensemble;
  double a3,w[3];
  double tmat[3][3];
  double r,r2,r3,r5,hmati[9];
  double *x,*y,*z,*qch,*hmat;

  x = coords->px;
  y = coords->py;
  z = coords->pz;
  qch = coords->qch;
  hmat = coords->hmat;
  iperd = simparms->iperd;
  iensemble = simparms->iensemble;
  n = simparms->natoms;
  gethinv9(hmat,hmati);

  for(i=0;i<n-1;i++){
    ii = 3*i;
    for(j=i+1;j<n;j++){
      jj = 3*j;
      if(!exclij(i,j,inter->nexclude,inter->exclude)){
	w[0] = x[i] - x[j];
	w[1] = y[i] - y[j];
	w[2] = z[i] - z[j];
	period(1,&w[0],&w[1],&w[2],hmat,hmati,iperd,iensemble);
	r2 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
	r = sqrt(r2);
	r3 = 1./(r2*r);
	r5 = r3/r2;
	
	/* get diatic */
	for(k=0;k<3;k++){
	  for(l=0;l<3;l++){
	    tmat[k][l] = 3.*r5*w[k]*w[l];
	    if(l==k) tmat[k][l] -= r3;
	  }
	}
	
	/* i particle */
	a3 = qch[j];
	for(k=0;k<3;k++){
	  for(l=0;l<3;l++){
	    amat[ii+k][jj+l] += tmat[k][l]*a3;
	    amat[ii+k][ii+l] -= tmat[k][l]*a3;
	  }
	}
	/* j particle */
	a3 = qch[i];
	for(k=0;k<3;k++){
	  for(l=0;l<3;l++){
	    amat[jj+k][ii+l] += tmat[k][l]*a3;
	    amat[jj+k][jj+l] -= tmat[k][l]*a3;
	  }
	}
      } /* end exclusion */
    }
  }
#ifdef PRINT_DEFIELD
  printf("defield\n");
  for(i=0;i<n*3;i++) {
    for(j=0;j<n*3;j++){
      printf("%g ",amat[i][j]);
    }
    printf("\n");
  }
  printf("\n");
#endif
}
/*---------------------------------------------------------------*/
void getmu(SIMPARMS *simparms,COORDS *coords,double *mu)
{
  int i,j,k,l,ii,jj,n3,n,*indx,iperd,iensemble;
  double a3,w[3];
  double **amat;
  double r,r2,r3,r5,hmati[9];
  double *x,*y,*z,*qch,*alpha,*hmat;

  x = coords->px;
  y = coords->py;
  z = coords->pz;
  qch = coords->qch;
  alpha = coords->alpha;
  hmat = coords->hmat;
  iperd = simparms->iperd;
  iensemble = simparms->iensemble;
  n = simparms->natoms;
  n3 = 3*n;
  amat = dmatrix(0,n3-1,0,n3-1);
  gethinv9(hmat,hmati);
    
  /* b vector and A matrix */
  /* set up the equation M x = b 
     M is the matrix with Tij's on the off diagonal block 
       and a^(-1) along the diagonal blocks
     x is the dipole that we want
     b is the electric field
     */

  for(i=0;i<3*n;i++) mu[i] = 0.0;
  for(i=0;i<n3;i++) for(j=0;j<n3;j++) amat[i][j] = 0.0;

  /* diagonal 3x3 blocks are alpha^(-1) easy since alpha is diagonal */
  for(i=0;i<n;i++) {
    ii = 3*i;
    if(alpha[i]==0.0) {
      md_error("Code not set up for zero alpha... try setting alpha small");
    }
  }

  for(i=0;i<n-1;i++){
    ii = 3*i;
    for(j=i+1;j<n;j++){

      w[0] = x[i] - x[j];
      w[1] = y[i] - y[j];
      w[2] = z[i] - z[j];
      period(1,&w[0],&w[1],&w[2],hmat,hmati,iperd,iensemble);
      r2 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
      r = sqrt(r2);
      r3 = 1./(r2*r);
      r5 = r3/r2;
      
      a3 = r3*qch[j];
      mu[ii  ] += w[0]*a3;
      mu[ii+1] += w[1]*a3;
      mu[ii+2] += w[2]*a3;
      
      jj = 3*j;
      a3 = r3*qch[i];
      mu[jj  ] -= w[0]*a3;
      mu[jj+1] -= w[1]*a3;
      mu[jj+2] -= w[2]*a3;
      
      /* get diatic */
      for(k=0;k<3;k++){
        for(l=0;l<3;l++){
          amat[ii+k][jj+l] = -3.*r5*w[k]*w[l];
          if(l==k) amat[ii+k][jj+l] += r3;
        }
      }

      /* reflect matrix (ie real symetric ) */
      for(k=0;k<3;k++){
        for(l=0;l<3;l++){
          amat[jj+k][ii+l] = amat[ii+l][jj+k];
        }
      }
    }
  }

  /* now that the system is set up 
     solve the thing for the dipoles */

#ifdef BLAS
  sol4musymMat(n3,amat,mu); 
#else
  indx = cmalloc(n3*sizeof(int));
  ludcmp(amat,n3,indx,&r);
  lubksb(amat,n3,indx,mu);
  free(indx);
#endif
  
  free_dmatrix(amat,0,n3-1,0,n3-1);
}
/*----------------------------------------------------------------*/
/* get the Electric Field from Ewald sums */
void getewef(SIMPARMS *simparms,COORDS *coords,double *bvec)
{
  int i,n,n3,ii;
  double *fxl,*fyl,*fzl,*qch;

  n = simparms->natoms;
  n3 = 3*n;

  qch = coords->qch;
  fxl = coords->fxl;
  fyl = coords->fyl;
  fzl = coords->fzl;

  for(i=0;i<n;i++)coords->fxl[i] = coords->fyl[i] = coords->fzl[i] = 0.0;  
  fk_ewald(simparms,coords);
  ecorr(simparms,coords); 
  for(i=0;i<n;i++){
    ii = i*3;

#ifdef HDEBUG
    printf("charge=%lg,fxl=%lg fyl=%lg fzl=%lg\n",qch[i],fxl[i],fyl[i],fzl[i]);
    printf(" OLD ELECTRIC FIELD:");
    printf("%lg %lg %lg\n",bvec[ii],bvec[ii+1],bvec[ii+2]);
#endif
    
    bvec[ii  ] = coords->fxl[i]/qch[i];
    bvec[ii+1] = coords->fyl[i]/qch[i];
    bvec[ii+2] = coords->fzl[i]/qch[i];

#ifdef HDEBUG
    printf(" NEW ELECTRIC FIELD:");
    printf("%lg %lg %lg\n",bvec[ii],bvec[ii+1],bvec[ii+2]);
#endif
  }
}
