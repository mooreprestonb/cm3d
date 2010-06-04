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


/* zero a double vector */

#include <math.h>
#include <stdlib.h>

void md_error(char *);

void azero_new(double x[],const int n)
{
  int i;
  for (i=0;i<n;i++){x[i]=0.;}
}
/*-----------------------------------------------------------*/
/* take projected array and expand them back into the full space */
void expandztox_new(double *x,double *y,double *z,double *q,double **d2v,
		int m,int n)
{
  double *qs;
  int i,j,ii,n3;

  n3 = n*3;
  qs = (double *)malloc(n3*sizeof(double));
  azero(qs,n3);
  
  for(j=0;j<m;j++){
    for(i=0;i<n3;i++){
      qs[i] += d2v[i][j]*q[j];
    }
  }
  for(i=ii=0;i<n;i++,ii+=3){
    x[i]=qs[ii];
    y[i]=qs[ii+1];
    z[i]=qs[ii+2];
  }
  free(qs);
}

/*-----------------------------------------------------------*/
/* take unprojected array and project them back into the subspace */
void projectxtoz_new(double *x,double *y,double *z,double *q,double **d2v,
		 int m,int n)
{
  int i,j,ii,n3;
  double *qs;

  n3 = n*3;
  qs = (double *)malloc(n3*sizeof(double));

  /* Set up the reduced variable */
  
  for(i=ii=0;i<n;i++,ii+=3){
    qs[ii  ] = x[i];    
    qs[ii+1] = y[i];    
    qs[ii+2] = z[i];
  }

  azero(q,m);
  
  for(j=0;j<m;j++){    
    for(i=0;i<n3;i++){
      q[j] += d2v[i][j]*qs[i];
    }
  }
  free(qs);
}
/*-------------------------------------------------------------*/
void mass_weight_new(int n,double *x,double *y,double *z,double *amass,int iflag)
{
  int i;
  double tmass;

  if(iflag==1){
    for(i=0;i<n;i++){
      tmass = sqrt(amass[i]);
      x[i] *= tmass; y[i] *= tmass; z[i] *= tmass;
    }
  } 
  if(iflag==-1){
    for(i=0;i<n;i++){
      tmass = 1./sqrt(amass[i]);
      x[i] *= tmass; y[i] *= tmass; z[i] *= tmass;
    }
  }
  if(iflag!=1&&iflag!= -1)md_error("(in mass_weight) call technical support");
}
/*----------------------------------------------------------------*/
