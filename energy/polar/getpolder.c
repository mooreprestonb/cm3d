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

/* paper Applequist 1972 JACS p2952
 * paper Applequist 1977 JCP p3455
 */

/* #define PRINT_T_DIPOLE */
/* #define HDEBUG2 */
#include "md.h"

/* calculate equation 1 in JCP paper 
 * final result is in getpolder's sum_matrix field */

/* two adjustable parameters */
double global_dadRn[2];

/* test routines */
void printmatrix(int n,double **a);

/*----------------------------------------*/
/* get A matrix which appeared in paper 2 JCP, equation 5
 * it's also the A matrix in equation 1 in JCP, paper 1 */
/* amat is the matrix with Tij's on the off diagonal block 
 * and alpha^(-1) along the diagonal blocks */
void getamatrix(SIMPARMS *simparms,COORDS *coords,double **amat,int nef,INTER *inter,double *bvec)
{
    int i,j,k,l,ii,jj,n3,n,iperd,iensemble;
    double a3,w[3];
    double r,r2,r3,r5, hmati[9];
    double *x,*y,*z,*qch,*alpha,*hmat;
    /* temporary cutoff changes */
    double p,sij,v,v3,v4,rdamp;
    p = 1./6.;
    /* for water- Mol. Phys. 1996, v87, 1333-1347 */
    rdamp = 1.662; /* HLA switch to input field */

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
    
    gethinv9(hmat,hmati);

    if(nef==1)
      for(i=0;i<n3;i++) bvec[i]=0.0;
    
    for(i=0;i<n3;i++)
      for(j=0;j<n3;j++)
        amat[i][j]=0.0;
    
    for(i=0;i<n;i++) {
      ii = 3*i;
      if(alpha[i]==0.0) {
        md_error("Code not set up for zero alpha... try setting alpha small");
      }
      else
        amat[ii][ii] = amat[ii+1][ii+1] = amat[ii+2][ii+2] = 1./alpha[i];
    }

    for(i=0;i<n-1;i++){
      ii = 3*i;
      for(j=i+1;j<n;j++){
        jj = 3*j;
        /* temporary cutoff changes */
        sij = rdamp * pow(alpha[i]*alpha[j],p);
        
        w[0] = x[i] - x[j];
        w[1] = y[i] - y[j];
        w[2] = z[i] - z[j];
        period(1,&w[0],&w[1],&w[2],hmat,hmati,iperd,iensemble);
        r2 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
        r = sqrt(r2);
        r3 = 1./(r2*r);
        r5 = r3/r2;
        v = (r < sij) ? r/sij : 1.;
        v3 = v * v * v;
        v4 = v3 * v;

        if(nef==1){  /* calculating Field w/o Ewald sums */
          if(!exclij(i,j,inter->nexclude,inter->exclude)){ 
            a3 = r3*qch[j]; 
            bvec[ii  ] += w[0]*a3; 
            bvec[ii+1] += w[1]*a3; 
            bvec[ii+2] += w[2]*a3; 
            
            a3 = r3*qch[i]; 
            bvec[jj  ] -= w[0]*a3; 
            bvec[jj+1] -= w[1]*a3; 
            bvec[jj+2] -= w[2]*a3; 
          }
        }
        
        for(k=0;k<3;k++){
          for(l=0;l<3;l++){
            amat[ii+k][jj+l] = -3.*r5*w[k]*w[l]*v4;
            if(l==k) amat[ii+k][jj+l] += r3*(4.*v3-3*v4);
          }
        }
        /* reflect matrix (ie real symmetric) since it's symmetric */
        for(k=0;k<3;k++){
          for(l=0;l<3;l++){
            amat[jj+k][ii+l] = amat[ii+l][jj+k];
          }
        }
      } /* end j */      
    } /* end i, amat has been constructed */
}

/* caculate Snk in equation (6), which is a unit vector from
 * atom j to atom i, this internal coordinate is related
 * with bond length*/
/* caculate the vector Snk which relate internal coordinate n
 * with atom k, put the vector in s */
/* actually, we only deal with two internal coordinates:
 * n = 0 is the first X-Y bond, n = 1 is the other X-Y bond */
void getSnkBond(COORDS *coords, int n, int k, double *s)
{
  int j=0;
  double normal = 0;
  double *x, *y, *z;
  
  x = coords->px;
  y = coords->py;
  z = coords->pz;
  
  /* loop over all atoms j bonded with atoms k */
  if( (n == 0 && k%3 == 2) || (n == 1 && k%3 == 1)){
    s[0] = 0.;
    s[1] = 0.;
    s[2] = 0.;
    return;
  }
  
  if(n == 0) /* n = 0, the first internal coordinate
              * then if k%3 ==0, k is X and Y is k+1
              * otherwise k is the Y atom */
    j = ((k%3 == 0) ? k+1 : k-1);
  else if(n == 1)
    j = ((k%3 == 0) ? k+2 : k-2);
  
  
  /* get the vector for atom j to atom i */
  s[0] = x[k] - x[j];
  s[1] = y[k] - y[j];
  s[2] = z[k] - z[j];
  
  normal = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
  if(normal == 0 ){
    printf("ERROR: exit");
    exit(1);
  }
  
  normal = sqrt(normal);
  s[0] /= normal;
  s[1] /= normal;
  s[2] /= normal; /* now s is the unit vector along j->i */
}

/* get r_k' in equation (6), which is derivative of atom k's
 * vector to normal mode imode , result is put in r*/
void getdRkdqt(COORDS *coords,double **fmat,int k,int imode,double *r)
{  
  int ii = k*3;
  double *amass = coords->amass;
  double w_mass;
  
  /* get the w1[] from force matrix, */
  w_mass = sqrt(amass[k]);
  
  r[0] = fmat[ii][imode]/w_mass;
  r[1] = fmat[ii+1][imode]/w_mass;
  r[2] = fmat[ii+2][imode]/w_mass; /* now w is the vector dr_k/dqimode */
}

/* alphap_der is actually the alpha_j^{-2}*alpha_j' in equation 2
 * this routine calculate all the atom's alphap_der of normal mode imode and store it in
 * alphap_der[] */
void  adjustAlphap_der(COORDS* coords,double **fmat,int n,int imode,double* alphap_der)
{
  int j;
  double r[3], s[3];
  
  double tmp;
  double *alpha = coords->alpha;
  
  /* this portion is for XY2, we have assumed that indice
   * 1,2 is for Y, 0 for X*/
  for(j=0;j<n;j++) { /* for each atom */
    /* tmp will hold a_j' */
    tmp = 0.0;
    
    if(j % 3 == 0){ /* this is X atom */
      getSnkBond(coords, 0, j, s);
      getdRkdqt(coords, fmat, j, imode, r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[0];
      
      getSnkBond(coords, 0, j+1, s);
      getdRkdqt(coords, fmat, j+1, imode, r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[0];
      
      getSnkBond(coords, 1, j, s);
      getdRkdqt(coords, fmat, j, imode, r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[0];
      
      getSnkBond(coords, 1, j+2, s);
      getdRkdqt(coords, fmat, j+2, imode,  r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[0];
      
    }
    else if(j % 3 == 1){ /* this is first Y atom */
      getSnkBond(coords, 0, j-1, s);
      getdRkdqt(coords, fmat, j-1, imode, r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[1];
      
      getSnkBond(coords, 0, j, s);
      getdRkdqt(coords, fmat, j, imode, r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[1];
      
    }
    else{ /* j%3 = 2, the second Y atom */
      getSnkBond(coords, 1, j-2, s);
      getdRkdqt(coords, fmat, j-2, imode, r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[1];
      
      getSnkBond(coords, 1, j, s);
      getdRkdqt(coords, fmat, j, imode, r);
      tmp += (s[0]*r[0] + s[1]*r[1] + s[2]*r[2]) * global_dadRn[1];
    }
    alphap_der[j] = tmp * pow(alpha[j%3], -2.) ;
    /* printf("%lg\t%lg\t%lg\n",tmp,alpha[j%3],alphap_der[j]); */
  }
}

/*----------------------------------------------------------------*/
/* get the A_{jk}' matrix appeared in paper JCP (paper 2)
 * A_{jj}'s is derivative to normal coordinate IMODE */
void getdadq(SIMPARMS *simparms,COORDS *coords,double **amat, double **fmat, int imode)
{
  int i,j,k,l,m,ii,jj,n,iperd,iensemble;
  double w[3],tmat[3][3][3];
  double r,r2,r3,r5,r7,hmati[9];
  double *x,*y,*z,*hmat;
  double *ux,*uy,*uz;
  double w1[3];  /* used to contain the force matrix element
                  * for caculating drdq.*/
  double *amass;
  double w_massi, w_massj;
#ifdef ADJUST
  double *alphap_der;
#endif
  amass = coords->amass;
  iperd = simparms->iperd;
  iensemble = simparms->iensemble;
  n = simparms->natoms;
  x = coords->px;  y = coords->py;  z = coords->pz;
  hmat = coords->hmat;
  ux = coords->ux;  uy = coords->uy;  uz = coords->uz;
  
#ifdef adjust  
  /* first get the diagonal element */
  alphap_der = (double *)malloc(n*sizeof(double));
  adjustAlphap_der(coords, fmat, n, imode, alphap_der);
  
  /* i-i particle polarizability */
  for(i=0;i<n;i++){
    printf("alphap_der=%lg\n",alphap_der[i]);
    ii = 3*i;
    for(l=0;l<3;l++)
      for(k=0;k<3;k++) {
        amat[ii+l][ii+k] = 0.0;
        if(l==k)amat[ii+l][ii+k] = -1.0 * alphap_der[i];
      }
  }
  free(alphap_der);
#else
  for(i=0;i<n;i++){
    ii = 3*i;
    for(l=0;l<3;l++)
      for(k=0;k<3;k++) 
        amat[ii+l][ii+k] = 0.0;

  }
#endif

  
  gethinv9(hmat,hmati);

  /* now the off-diagonal i!=j */
  for(i=0; i<n; i++){
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
        for(k=0;k<3;k++){
          for(l=0;l<3;l++){
            for(m=0;m<3;m++){
              tmat[k][l][m] = 15.*r7*w[k]*w[l]*w[m];
              if(l==k) tmat[k][l][m] += -3.*w[m]*r5;
              if(l==m) tmat[k][l][m] += -3.*w[k]*r5;
              if(m==k) tmat[k][l][m] += -3.*w[l]*r5;
            }
          }
        }
        
        /* get the w1[] from force matrix, */
        w_massi = 1.0/sqrt(amass[i]);
        w_massj = 1.0/sqrt(amass[j]);
        
        w1[0] = fmat[ii][imode]*w_massi - fmat[jj][imode]*w_massj;
        w1[1] = fmat[ii+1][imode]*w_massi - fmat[jj+1][imode]*w_massj;
        w1[2] = fmat[ii+2][imode]*w_massi - fmat[jj+2][imode]*w_massj;
        
        /* printf("%lg,%lg,%lg\n",w1[0],w1[1],w1[2]);  */
        
        /* i-j particle polarizability */
        for(k=0;k<3;k++){
          for(l=0;l<3;l++){
            amat[ii+k][jj+l] = (w1[0]*tmat[k][l][0] +
               w1[1]*tmat[k][l][1] +
               w1[2]*tmat[k][l][2]);
          }
        }
      } /* end if (i!=j) */
    }/* end j loop */
  }/* end i loop */
  
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

/*---------------------------------------------------------------*/
/* caculate alpha derivatives to normal coordinate IMODE
 * using equation 1 in JCP paper 1 */
void getpolder(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *nmodes)
{
    int imode,i,j,k,l,n3,n;
    double **akeep;
    double **amat,**arest;
    double **aniso, **sum_matrix;
    double *dvec; /*dummy electric field variable for getamatrix */

#ifndef BLAS
    int *indx;
    double r,*col,**ainv;
#endif
    n = simparms->natoms;
    n3 = 3*n;
    dvec = NULL;
#ifndef BLAS
    col = malloc(n3*sizeof(double));
    indx = malloc(n3*sizeof(int));
    ainv = dmatrix(0,n3-1,0,n3-1);
#endif
    amat  = dmatrix(0,n3-1,0,n3-1);
    akeep = dmatrix(0,n3-1,0,n3-1);
    arest = dmatrix(0,n3-1,0,n3-1);
    aniso = dmatrix(0,2,0,2);
    sum_matrix = dmatrix(0,2,0,2);

#ifdef OPT /* optimized derivatives */
    global_dadRn[0] = -0.392906;
    global_dadRn[1] = 0;
    printf("global0=%lg\tglobal1=%lg\n",global_dadRn[0],global_dadRn[1]);
#else
    global_dadRn[0] = 0;
    global_dadRn[1] = 0.0;
   /* printf("Warning: no optimization\n"); */
#endif

    /* get amat, which is A in equation 5 paper JACS
     * which is also the inverse B matrix in eqn 1 paper JCP */
#ifdef BLAS
    getamatrix(simparms,coords,amat,0,inter,dvec);
    invsymMatrix(n3,amat); /* get inverse of amat */
#endif
    
    for(imode = 0; imode < n3; imode++){    
      /* amat is constructed, lets get A_{jk}' in paper JCP
       * equation 1, which is stored in arest matrix */
      getdadq(simparms,coords,arest,nmodes->fcmat,imode);
#ifdef BLAS
      /* Calculate A_{ij}^-1 * A_{jk}' * A_{kl}^-1 eqn 1 paper JCP */
      matmul(n3,amat,arest,akeep);
      matmul(n3,akeep,amat,arest);
#else
      /* must refresh amat w/ each loop b/c ludcmp destroys it */
      getamatrix(simparms,coords,amat,0,inter,dvec);
      ludcmp(amat,n3,indx,&r);
      for(i=0;i<n3;i++){/* Calculate A_{ij}^-1 * A_{jk}' */
        for(j=0;j<n3;j++) col[j] = arest[i][j];
        lubksb(amat,n3,indx,col);
        for(j=0;j<n3;j++) akeep[j][i] = col[j];
      }
      for(i=0;i<n3;i++){/* Calculate A_{ik}* A_{kl}^-1 */
        for(j=0;j<n3;j++) col[j] = 0.0;
        col[i] = 1.0;
        lubksb(amat,n3,indx,col);
        for(j=0;j<n3;j++) ainv[j][i] = col[j];
      }
      matmul(n3,akeep,ainv,arest);
#endif
      /* now we have stored data in arest 3N x 3N matrix,
       * time to extract it and add them up to form (3 x 3) matrix. */
      for(i=0;i<3;i++)
        for(j=0;j<3;j++)sum_matrix[i][j]=0.0;
      
      for(i=0; i<3;i++)
        for(j=0;j<3;j++)
          for(k=i;k<n3;k+=3)
            for(l=j;l<n3;l+=3)
              sum_matrix[i][j] += arest[k][l];
      for(i=0; i<3;i++)
        for(j=0;j<3;j++)
          sum_matrix[i][j] = -sum_matrix[i][j];
      /* ok, now the polarizability derivative to normal coordinate
         imode has been stored in sum_matrix     */
      nmodes->iso_raman[imode] = 0.0;
      for(i=0; i<3;i++)
        nmodes->iso_raman[imode] += sum_matrix[i][i];
      nmodes->iso_raman[imode] /= 3.0;
      for(i=0; i<3;i++)
        for(j=0;j<3;j++){
          aniso[i][j] = sum_matrix[i][j];
          if(i==j)aniso[i][j] -= nmodes->iso_raman[imode];
        }
      matmul(3,aniso,aniso,sum_matrix);
      
      nmodes->aniso_raman[imode] = 0.0;
      for(i=0; i<3;i++)
        nmodes->aniso_raman[imode] += sum_matrix[i][i];
      nmodes->iso_raman[imode] *= nmodes->iso_raman[imode];
    }
/*     for(imode = 0; imode < n3; imode++) */
/*       printf("weight[%d]=%lf\n",imode,nmodes->aniso_raman[imode]); */
    
    free_dmatrix(amat,0,n3-1,0,n3-1);    free_dmatrix(akeep,0,n3-1,0,n3-1);
    free_dmatrix(arest,0,n3-1,0,n3-1);    free_dmatrix(sum_matrix,0,2,0,2);   
    free_dmatrix(aniso,0,2,0,2);
#ifndef BLAS
    free_dmatrix(ainv,0,n3-1,0,n3-1); 
    free(col);
#endif
}
/*----------------------------------------------------------------*/
void printmatrix(int n,double **a)
{
  int i,j;

  for(i=0;i<n;i++){
    for(j=0;j<n;j+=3) 
      printf("%-10.6g\t%-10.6g\t%-10.6g\n",a[j][i],a[j+1][i],a[j+2][i]); 
    printf("\n");
  }
}
