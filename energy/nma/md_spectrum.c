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

/* subroutine to get then spectrum from the charges and eigenvalues */

#include "md.h"

/* #define DEBUG */

void get_spectrum(int natoms,double *amass,double *qch,double **fcmat,
		  double *spectrum,double **dipder)
{
  int i,j,k;
  double ex,ey,ez,qm;
  double **mat;

  /* allocate temporary vector */
  mat = dmatrix(0,3*natoms-1,0,3*natoms-1);

  for(i=0;i<3*natoms;i++){
    for(j=0;j<3*natoms;j++){
      qm = 0.;
      for(k=0;k<3*natoms;k++) qm += dipder[j][k]*fcmat[k][i]/sqrt(amass[k/3]);
      mat[j][i] = qm;
    }
  }

  for(i=0;i<natoms*3;i++){
    ex = ey = ez = 0.;
    for(j=0,k=0;k<natoms;k++,j+=3){
      ex += mat[j  ][i];
      ey += mat[j+1][i];
      ez += mat[j+2][i];
    }
    /* a(k) = 1/(12 cspeed epsilon_0) * [Sum (charge*vector/mass^(1/2))] */
    /* McQuarrie(delta function) or Atkins (units) or 
       Willson+Decises+Cross */
    spectrum[i] = (ex*ex+ey*ey+ez*ez)*M_PI/(3.*CSPEED);
#ifdef DEBUG
    printf("spectrum alpha %d = %g\n",i,spectrum[i]);
#endif
  }
  free_dmatrix(mat,0,3*natoms-1,0,3*natoms-1); /* free vector */
}

/*------------------------------------------------------------------*/
void get_specq(int natoms,double *amass,double *qch,double **fcmat,
                  double *spectrum)
{
  int i,j,k;
  double ex,ey,ez,qm;

  for(k=0;k<natoms*3;k++){
    ex = ey = ez = 0.;
    for(j=0,i=0;i<natoms;i++,j+=3){
      qm = qch[i]/sqrt(amass[i]);
      ex += qm*fcmat[j+0][k];
      ey += qm*fcmat[j+1][k];
      ez += qm*fcmat[j+2][k];
    }
    /* a(k) = 1/(12 cspeed epsilon_0) * [Sum (charge*vector/mass^(1/2))] */
    /* McQuarrie(delta function) or Atkins (units) or Willson+Decises */
    spectrum[k] = (ex*ex+ey*ey+ez*ez)*M_PI/(3.*CSPEED);
#ifdef DEBUG
    printf("spectrum charge %d = %g\n",k,spectrum[k]);
#endif
  }
}

/*------------------------------------------------------------------*/

void get_spect_full(int natoms,double *amass,double *qch,double **fcmat,
   double *spectrum,double **dipder)
{
  int i,j,k;
  double ex,ey,ez,qm;
  double **mat;
  
  /* allocate temporary vector */
  mat = dmatrix(0,3*natoms-1,0,3*natoms-1);

  for(i=0;i<3*natoms;i++){
    for(j=0;j<3*natoms;j++){
      qm = 0.;
      for(k=0;k<3*natoms;k++) qm += dipder[j][k]*fcmat[k][i]/sqrt(amass[k/3]);
      mat[j][i] = qm;
    }
  }

  for(i=0;i<natoms*3;i++){
    ex = ey = ez = 0.;
    for(j=0,k=0;k<natoms;k++,j+=3){
      qm = qch[k]/sqrt(amass[k]);
      ex += mat[j  ][i];
      ey += mat[j+1][i];
      ez += mat[j+2][i];

      ex += qm*fcmat[j+0][i];
      ey += qm*fcmat[j+1][i];
      ez += qm*fcmat[j+2][i];
    }
    /* a(k) = 1/(12 cspeed epsilon_0) * [Sum (charge*vector/mass^(1/2))] */
    /* McQuarrie(delta function) or Atkins (units) or 
       Willson+Decises+Cross */
    spectrum[i] = (ex*ex+ey*ey+ez*ez)*M_PI/(3.*CSPEED);
#ifdef DEBUG
    printf("spectrum full %d = %g\n",i,spectrum[i]);
#endif
  }
  free_dmatrix(mat,0,3*natoms-1,0,3*natoms-1); /* free vector */
}

/*------------------------------------------------------------------*/
void avg_spec(int natoms,NMODES *nmodes)
{
  int i,n3;
  double fr;

  n3 = natoms*3;
  /* average the strength of the spectrum for each mode */
  for(i=0;i<n3;i++){
    nmodes->avg_strn[i] += nmodes->fricmodq[i] + nmodes->fricmod[i];
  }

  /* average the frequencies for each mode sqrt of the eigenvalue */
  for(i=0;i<n3;i++){
    fr=nmodes->dn[i];
    if (fr<0.0) fr= sqrt(-fr);
    else fr=sqrt(fr);

    nmodes->avg_freq[i] += fr;
  }
}

/*------------------------------------------------------------------*/

