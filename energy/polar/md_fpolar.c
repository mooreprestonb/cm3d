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
  
#include "md.h"
#ifdef RETARDED
/*====================================================================*/
/* da/dr = B * dA/dr * B, where B = A^-1  (see Applequist refs.) */  
void getdadr(SIMPARMS *simparms,COORDS *coords,double dadr)
{
  int i,j,k,l,m,ii,jj;
  int n3,iperd,iensemble;
  double w[3],tmat[3][3][3];
  double r,r2,r3,r5,r7,hmati[9];
  double *x,*y,*z,*hmat;

  nef = NULL;
  amat  = dmatrix(0,n3-1,0,n3-1);
  amat2  = dmatrix(0,n3-1,0,n3-1);

  amass = coords->amass;
  iperd = simparms->iperd;
  iensemble = simparms->iensemble;
  n3 = simparms->natoms*3;
  x = coords->px;  y = coords->py;  z = coords->pz;
  hmat = coords->hmat;
  
  gethinv9(hmat,hmati);

  /* construct da/dr off-diagonal i!=j */
  /* diagonal elements are zero */
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
            } /* end m loop */
            amat[ii+k][jj+l] = (tmat[k][l][0] +
				tmat[k][l][1] +
				tmat[k][l][2]);
          } /* end l loop */
        } /* end k loop */
      } /* end if (i!=j) */
    }/* end j loop */
  }/* end i loop */


}
#endif
/*====================================================================*/
void force_dipole(int nlen,double *dx,double *dy,double *dz,double *qch,
		  double *alpha,double *fx,double *fy,double *fz)
{
  int i;
  double r,r2,r3,r5;
  double mu[3];

  for(i=0;i<nlen;i++){
      /* temporary cutoff changes */
      /* sij = rdamp * pow(alpha[i]*alpha[j],p); */
        

      r2 = dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i];
      r = sqrt(r2);
      r3 = 1./(r2*r);
      r5 = r3/r2;
/*       v = (r < sij) ? r/sij : 1.; */
/*       v3 = v * v * v; */
/*       v4 = v3 * v; */
      mu[0] = alpha[i]*qch[i]*dx[i]*r3;
      mu[1] = alpha[i]*qch[i]*dy[i]*r3;
      mu[2] = alpha[i]*qch[i]*dz[i]*r3;

      fx[i] -= 3.* mu[0] * dx[i] * dx[i] *r5;
      fy[i] -= 3.* mu[1] * dy[i] * dy[i] *r5;
      fz[i] -= 3.* mu[2] * dz[i] * dz[i] *r5;

      fx[i] += mu[0] * r3;
      fy[i] += mu[1] * r3;
      fz[i] += mu[2] * r3;

      fx[i] *= 2;
      fy[i] *= 2;
      fz[i] *= 2;

  }
}

  
/*====================================================================*/
double pot_dipole(int nlen,double *dx,double *dy,double *dz,double *qch,
  double *alpha)
{
  int i;
  double r2,r6;
  double pot;
  
  pot = 0.;
  for(i=0;i<nlen;i++){
    /* temporary cutoff changes */
    /* sij = rdamp * pow(alpha[i]*alpha[j],p); */
    
    r2 = dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i];
    r6 = r2 * r2 *r2;
    
    /*       v = (r < sij) ? r/sij : 1.; */
    /*       v3 = v * v * v; */
    /*       v4 = v3 * v; */
    pot += alpha[i]*qch[i]*dx[i]*dx[i]*r6;
    pot += alpha[i]*qch[i]*dy[i]*dy[i]*r6;
    pot += alpha[i]*qch[i]*dz[i]*dz[i]*r6;

  }

  return(pot);
}

  

