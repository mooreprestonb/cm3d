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
#define DEBUG
#define PRINT_MATRIX

/* get_fcmat is different from fcmat which is in md_nmafcntr.c */
/*--------------------------------------------------------------*/
void get_fcmat(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **fcmat)
{
  int i,j,n3;
#ifdef PRINT_MATRIX
  FILE *fp;
#endif

  n3 = simparms->natoms*3;
  for(i=0;i<n3;i++) for(j=0;j<n3;j++) fcmat[i][j] = 0.0;

#ifdef NUMERIC
  md_stdout("Getting force constant matrix numerically");
  fcbondn(coords,fcmat);
  fcbendn(coords,fcmat);
  fctorsn(coords,fcmat);
  fconefourn(coords,fcmat);
  fcbondxn(coords,fcmat);
  fcintern(simparms,coords,inter,fcmat);
#else
  md_stdout("\n ****Getting force constant matrix*****"); 
  fcbond(coords,fcmat);
  fcbend(coords,fcmat);
  fctors(coords,fcmat);
  fconefour(coords,fcmat);
  fcbondx(coords,fcmat);
  fcinter(simparms,coords,inter,fcmat);
#endif

  massweight_fcmat(simparms->natoms,fcmat,coords->amass);
  reflect_fcmat(simparms->natoms,fcmat);

#ifdef PRINT_MATRIX
  fp = fopen("hess.dat","w");
  for(i=0;i<n3;i++) 
    for(j=0;j<n3;j++) fprintf(fp,"%d %d %.12g\n",i,j,fcmat[i][j]);
#endif
  
}
/*-----------------------------------------------------------*/
void massweight_fcmat(int natoms,double **fcmat,double *amass)
{
  int i,j,ii,jj;
  double rmass;

  for(i=0;i<natoms;i++){
    for(j=i;j<natoms;j++){
      rmass = 1./sqrt(amass[i]*amass[j]);
      
      ii=3*i-1; jj=3*j-1;

      fcmat[ii+1][jj+1] *= rmass;
      fcmat[ii+1][jj+2] *= rmass;
      fcmat[ii+1][jj+3] *= rmass;
      fcmat[ii+2][jj+1] *= rmass;
      fcmat[ii+2][jj+2] *= rmass;
      fcmat[ii+2][jj+3] *= rmass;
      fcmat[ii+3][jj+1] *= rmass;
      fcmat[ii+3][jj+2] *= rmass;
      fcmat[ii+3][jj+3] *= rmass;
    }
  }
}

/*---------------------------------------------------------*/

void reflect_fcmat(int natoms,double **fcmat)
{
  int i,j,ii,jj;
  /* fill in other half */

  for(i=0;i<natoms-1;i++){
    for(j=i+1;j<natoms;j++){
      ii=3*i-1;
      jj=3*j-1;
      
      fcmat[jj+1][ii+1] = fcmat[ii+1][jj+1];
      fcmat[jj+2][ii+1] = fcmat[ii+1][jj+2];
      fcmat[jj+3][ii+1] = fcmat[ii+1][jj+3];
      fcmat[jj+1][ii+2] = fcmat[ii+2][jj+1];
      fcmat[jj+2][ii+2] = fcmat[ii+2][jj+2];
      fcmat[jj+3][ii+2] = fcmat[ii+2][jj+3];
      fcmat[jj+1][ii+3] = fcmat[ii+3][jj+1];
      fcmat[jj+2][ii+3] = fcmat[ii+3][jj+2];
      fcmat[jj+3][ii+3] = fcmat[ii+3][jj+3];
    }
  }
}
/*--------------------------------------------------------------*/
