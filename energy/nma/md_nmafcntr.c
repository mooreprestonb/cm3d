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

/* #define NUMERIC  */
/* #define PRINT_MATRIX */
/* #define PRINT_MATRIX_MOL */
/*---------------------------------------------------------------*/

void fmat(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *nmodes)
{
  int i,j,imol,ispec,ioff,n3;
  double **fcmat;
#ifdef PRINT_MATRIX
  FILE *fp;
#endif

  n3 = simparms->natoms*3;
  fcmat = nmodes->fcmat;

  gethinv9(coords->hmat,coords->hmati);
  for(i=0;i<n3;i++) for(j=0;j<n3;j++) fcmat[i][j] = 0.;

#ifdef NUMERIC
  md_stdout("Getting force constant matrix numerically");

  fcbondn(coords,fcmat);
  fcbendn(coords,fcmat);
  fctorsn(coords,fcmat);
  fconefourn(coords,fcmat);
  fcbondxn(coords,fcmat); 

  /* if overlap copy internal forces to molmat */
  if(nmodes->flo){
    /* we have to generalize to loop over species and only include the
       intramolecular part of the force matrix */
    ioff = 0;
    for(ispec=0;ispec<nmodes->nspec;ispec++){
      for(imol=0;imol< nmodes->nmol[ispec]; imol++){
        for(i=0;i<nmodes->napm[ispec]*3;i++){
          for(j=0;j<nmodes->napm[ispec]*3;j++){
            nmodes->molmat[ioff][j] = 0.0;
          }
          ioff++;
        }
      }
    }
    ioff = 0;
    for(ispec=0;ispec<nmodes->nspec;ispec++){
      for(imol=0;imol< nmodes->nmol[ispec]; imol++){
        copy(3*nmodes->napm[ispec],3*nmodes->napm[ispec],fcmat,ioff,ioff,
           nmodes->molmat,ioff,0);
        ioff += 3*nmodes->napm[ispec];
      }
    }
  }
  fcintern(simparms,coords,inter,fcmat);
#else
  /* md_stdout("Getting force constant matrix"); */
  fcbond(coords,fcmat);
  fcbend(coords,fcmat);
  fctors(coords,fcmat);
  fconefour(coords,fcmat);
  fcbondx(coords,fcmat);
   /* if overlap copy internal forces to molmat */
  if(nmodes->flo){
    /* we have to generalize to loop over species and only include the
       intramolecular part of the force matrix */
    ioff = 0;
    for(ispec=0;ispec<nmodes->nspec;ispec++){
      for(imol=0;imol< nmodes->nmol[ispec]; imol++){
        for(i=0;i<nmodes->napm[ispec]*3;i++){
          for(j=0;j<nmodes->napm[ispec]*3;j++){
            nmodes->molmat[ioff][j] = 0.0;
          }
          ioff++;
        }
      }
    }
    ioff = 0;
    for(ispec=0;ispec<nmodes->nspec;ispec++){
      for(imol=0;imol< nmodes->nmol[ispec]; imol++){
        copy(3*nmodes->napm[ispec],3*nmodes->napm[ispec],fcmat,ioff,ioff,
           nmodes->molmat,ioff,0);
        ioff += 3*nmodes->napm[ispec];
      }
    }
  }
  fcinter(simparms,coords,inter,fcmat);
#endif

  massweight_fcmat(simparms->natoms,fcmat,coords->amass);
  reflect_fcmat(simparms->natoms,fcmat);
  
  if(nmodes->flo){
    ioff = 0;
    for(ispec=0;ispec<nmodes->nspec;ispec++){
      for(imol=0;imol< nmodes->nmol[ispec]; imol++){
        massweight_fcmat(nmodes->napm[ispec],&nmodes->molmat[ioff*3],
           &coords->amass[ioff]);
        reflect_fcmat(nmodes->napm[ispec],&nmodes->molmat[ioff*3]);
        ioff += nmodes->napm[ispec];
      }
    }
  }

#ifdef PRINT_MATRIX
  fp = fopen("hess.dat","w");
  for(i=0;i<n3;i++) 
    for(j=0;j<n3;j++) fprintf(fp,"%d %d %.12g\n",i,j,fcmat[i][j]);
  fclose(fp);
  
  /* Mathematic format
  printf("{");
  for(i=0;i<3*simparms->natoms-1;i++){
    printf("{");
    for(j=0;j<3*simparms->natoms-1;j++)
      printf("%lg,",fcmat[i][j]);
    printf("%lg},\n",fcmat[i][j]);
  }
  printf("{");
  for(j=0;j<3*simparms->natoms-1;j++)
    printf("%lg,",fcmat[i][j]);
  printf("%lg}}\n",fcmat[i][j]);
  */
#ifdef PRINT_MATRIX_MOL
  if(nmodes->flo){
    printf("\n{");
    for(ispec=0;ispec<nmodes->nspec;ispec++){
      for(imol=0;imol< nmodes->nmol[ispec]; imol++){
	for(i=0;i<nmodes->napm[ispec]*3;i++){
	  printf("{");
	  for(j=0;j<nmodes->napm[ispec]*3-1;j++)
	    printf("%lg,",molmat[i][j]);
	  printf("%lg},\n",molmat[i][j]);
	}
      }
    }
    printf("}\n",molmat[i][j]);
  }
#endif
#endif
}
/*------------------------------------------------------------------------*/
void copy(int nr,int nc,double **m1,int m1r,int m1c,
	  double **m2,int m2r,int m2c)
{
  int i,j;

  for(i=0;i<nr;i++) for(j=0;j<nc;j++) m2[i+m2r][j+m2c] = m1[i+m1r][j+m1c];
}
/*------------------------------------------------------------------*/
