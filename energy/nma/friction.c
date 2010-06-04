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


/* subroutine to get the friction from the eigenvectors */

#include "md.h"

/* #define OMEGA2 */
#define BOND
/* #define FRIC */

void friction(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *nmodes)
{
  int i,j,natoms3;
  double *vec,*dp,vmax;
#ifdef FRIC
  int k;
  double wt;
#endif

  natoms3 = simparms->natoms*3;

  /* construct df/dq = df/dr * dr/dq  */
  /* watch out for mass weighting not taken into consideration yet :-( */
  
  vec  = (double *)cmalloc(natoms3*sizeof(double));

  /* project susceptibilites onto mode of choice 
     (along the bond stretching mode) */

  nmodes->fricmod[0] = nmodes->molmat[0][5];
  nmodes->fricmod[1] = nmodes->molmat[1][5];
  nmodes->fricmod[2] = nmodes->molmat[2][5];
  nmodes->fricmod[3] = nmodes->molmat[3][5];
  nmodes->fricmod[4] = nmodes->molmat[4][5];
  nmodes->fricmod[5] = nmodes->molmat[5][5];

  for(i=6;i<natoms3;i++) nmodes->fricmod[i] = 0.;

  /* get overlaps */
  vecmat(natoms3,nmodes->fcmat,nmodes->fricmod,vec);

#ifdef FRIC
  /* find bond "renomalized" frequency */
  wt = vmax = 0.;
  k = j = 0;
  for(i=0;i<natoms3;i++) {
    if(nmodes->dn[i]<0) k++;
    else wt += vec[i]*vec[i]/nmodes->dn[i];
    if(vec[i]*vec[i] > vmax){
      vmax = vec[i]*vec[i];
      j = i;
    } 
  }
  wt = (1./wt)*(1./(1.-(double)k/(double)natoms3));
  printf("Coordinate %d freq = %g (overlap = %g) renorm freq = %g %g %d\n",
	 j,sqrt(nmodes->dn[j]),vmax,sqrt(wt),
	 (1./(1.-(double)k/(double)natoms3)),k);

  /* square the overlap this takes care of the 1/w^2 as well) */
  dp = nmodes->fricmod;
  for(i=0;i<natoms3;i++) {
    vmax = vec[i]*(nmodes->dn[i]-wt);
#ifdef OMEGA2
    dp[i] = vmax*vmax/(nmodes->dn[i]);
#else
    dp[i] = vmax*vmax;
#endif
  }
#endif /* define fric */

#ifdef BOND
  /* square the overlap this takes care of the 1/w^2 as well) */
  dp = nmodes->fricmod;
  vmax = 0;
  j = 0;
  for(i=0;i<natoms3;i++) {
#ifdef OMEGA2
    dp[i] = vec[i]*vec[i]/(nmodes->dn[i]);
#else
    dp[i] = vec[i]*vec[i];
#endif
    if(vec[i]*vec[i] > vmax){
      vmax = vec[i]*vec[i];
      j = i;
    } 
  }
  printf("Coordinate %d freq = %g (overlap = %g) \n",
	 j,sqrt(nmodes->dn[j]),vmax);
#endif

  /* free temporary storage */
  free(vec);
}
