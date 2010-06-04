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

/* subroutine to overlap force with eigenvectors */

#include "md.h"

void proj_force(int natoms,COORDS *coords,double **fcmat,double **molmat,
		int nspec,int *nmol,int *napm,
		double *suspx,double *suspy,double *suspz)
{
  int i,j,k,natoms3;

  natoms3 = natoms*3;

  /* zero all forces not involed in mode of choice */
  for(i=3;i<natoms;i++){
    suspx[i] = suspy[i] = suspz[i] = 0.;
  }
  /* project susceptibilites onto mode of choice */
  suspx[0] *= molmat[0][5];
  suspy[0] *= molmat[1][5];
  suspz[0] *= molmat[2][5];
  suspx[1] *= molmat[3][5];
  suspy[1] *= molmat[4][5];
  suspz[1] *= molmat[5][5];

  for(i=0;i<natoms3;i++){
    for(k=j=0;j<natoms3;j+=3,k++){
      fcmat[j  ][i] *= suspx[k];
      fcmat[j+1][i] *= suspy[k];
      fcmat[j+2][i] *= suspz[k];
    }
  }
}

/*---------------------------------------------------------------------*/
void get_fric(int natoms,int nspec,int *nmol,int *napm,
	      double **fcmat,double *fricmod)
{
  int i,j,natoms3;
  double sum;
  
  natoms3 = 3*natoms;
  for(i=0;i<natoms3;i++) fricmod[i] = 0.0;
  
  /*
  for(l=0;l<natoms3;l++){
    ioff = 0;
    for(ispec=0;ispec<nspec;ispec++){
      napm3 = napm[ispec]*3;
      for(imol=0;imol<nmol[ispec];imol++){
	for(k=0;k<napm3;k++){
	  sum = 0.;
	  j = ioff*3+napm3;
	  for(m=ioff*3;m<j;m++){sum += molmat[m][k]*molmat[m][k];}
	  fricmod[l] += sum;
	}
	ioff += napm[ispec];
      }
    }
  }
  */

  for(i=0;i<natoms3;i++){
    sum = 0.;
    for(j=0;j<natoms3;j++){
      sum += fcmat[j][i]*fcmat[j][i];
    }
    fricmod[i] += sum;
  }
}
/*--------------------------------------------------------------------*/
void get_suscep(SIMPARMS *simparms,COORDS *coords,INTER *inter,
		int nspec,int *nmol,int *napm,
		double *suspx,double *suspy,double *suspz)
{
  int i,j;
  double **fmat;
  /* construct df/dq = df/dr * dr/dq  */
  /* watch out for mass weighting not taken into consideration yet :-( */
  
  fmat = dmatrix(0,3*simparms->natoms-1,0,3*simparms->natoms-1);
  fcinter(simparms,coords,inter,fmat);

  for(j=i=0;i<simparms->natoms;i++,j+=3){
    suspx[i] = fmat[j  ][j  ];
    suspy[i] = fmat[j+1][j+1];
    suspz[i] = fmat[j+2][j+2];
  }
  free_dmatrix(fmat,0,3*simparms->natoms-1,0,3*simparms->natoms-1);
}

/*--------------------------------------------------------------------*/
/*
  for(i=0;i<natoms3;i++){
    ioff=0;
    for(ispec=0;ispec<nspec;ispec++){
      napm3 = napm[ispec]*3;
      for(imol=0;imol<nmol[ispec];imol++){
	for(k=0;k<napm3;k++){
	  j = ioff*3+napm3;
	  for(l=0,m=ioff*3;m<j;m+=3,l++){
	    
#ifdef glop
	    molmat[m  ][k] *= suspx[ioff+l]*(fcmat[m  ][i]-molmat[m  ][k]);
	    molmat[m+1][k] *= suspy[ioff+l]*(fcmat[m+1][i]-molmat[m+1][k]);
	    molmat[m+2][k] *= suspz[ioff+l]*(fcmat[m+2][i]-molmat[m+2][k]);
#endif
	    molmat[m  ][k] *= suspx[ioff+l]*fcmat[m  ][i];
	    molmat[m+1][k] *= suspy[ioff+l]*fcmat[m+1][i];
	    molmat[m+2][k] *= suspz[ioff+l]*fcmat[m+2][i];
	  }
	}
	ioff += napm[ispec];
      }
    }
  }  
  */
  
