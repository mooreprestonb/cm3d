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

#define NHILL_ALLOC 1000
  
static int nhills_alloc=0;

void hills_init(){  
  /* if(restart){ */
  /* old hills from restart file containing:                         */
  /* time,hilldepth,pcolvar(1:ncolvar),pistcolvar(1:ncolvar),hillwidth(1:ncolvar)*/
}

void add_hills(COLVAR *colvar)
{
  int icol,mbytes,nhills;
  ++(colvar->nhills);
  nhills = (colvar->nhills);

  if(nhills > nhills_alloc){
    if(nhills_alloc ==0 ){
      nhills_alloc = NHILL_ALLOC;
      mbytes = nhills_alloc*sizeof(double);
      colvar->t_hilldepth =(double *)malloc(mbytes);
      for(icol=0;icol<colvar->ncolvar;++icol){
	colvar->t_hillwidth[icol] = (double *)malloc(mbytes);
	colvar->t_pcolvar[icol] =  (double *)malloc(mbytes);
      }
    } else {   
      nhills_alloc += NHILL_ALLOC;
      
      mbytes = nhills_alloc*sizeof(double);
      colvar->t_hilldepth =(double *)realloc(colvar->t_hilldepth,mbytes);
      for(icol=0;icol<colvar->ncolvar;++icol){
	colvar->t_hillwidth[icol] = 
	  (double *)realloc(colvar->t_hillwidth[icol],mbytes);
	colvar->t_pcolvar[icol] = 
	  (double *)realloc(colvar->t_pcolvar[icol],mbytes);
      }
    }
  }
  colvar->t_hilldepth[nhills-1] = colvar->hilldepth;
  for(icol=0;icol<colvar->ncolvar;++icol){
    colvar->t_hillwidth[icol][nhills-1]=colvar->hillwidth[icol];
    colvar->t_pcolvar[icol][nhills-1]=colvar->pcolvar[icol];
  }
}

/*-----------------------------------------------------------------------*/
int cv_displacement(COLVAR *colvar)
{
  int icol;
  double dv,sdv2;
  int nhills = colvar->nhills;
  sdv2=0;

  for(icol=0;icol<colvar->ncolvar;++icol){
    dv = colvar->pcolvar[icol] - colvar->t_pcolvar[icol][nhills-1];
    sdv2 += dv*dv;
  }
  if(sdv2 > colvar->mindtol2) {
    return 1;
  }
  return 0;
}

/*-----------------------------------------------------------------------*/
void hills(SIMPARMS *simparms,COORDS *coords)
{
  int icol;
  double *t_hilldepth  = coords->colvar.t_hilldepth;
  double **t_hillwidth  = coords->colvar.t_hillwidth;
  double **t_pcolvar    = coords->colvar.t_pcolvar;
  double *f_hills = coords->colvar.f_hills;
  double *pcolvar = coords->colvar.pcolvar;
  int ncolvar = coords->colvar.ncolvar;
  int nhills = coords->colvar.nhills;
  int *lhills = coords->colvar.lhills;
  int ihill,dt;
  double *dp, dp2;
  double ecolvar = 0;
  double gauss;

  /* shall we add a hill */
  ihill=0;
  dt = simparms->istep - coords->colvar.lasthill;
  if(dt > coords->colvar.maxstephill){
    ihill = 1;
  }
  else if(dt > coords->colvar.minstephill){
    ihill = cv_displacement(&(coords->colvar));
  }
  if(ihill==1) {
    add_hills(&(coords->colvar));
    coords->colvar.lasthill = simparms->istep;
  }

  /* compute forces from hills on colvars */
  /* todo:add parallel gauss correction and hilldepth tuning              */

  dp = (double *)malloc(ncolvar*sizeof(double));

  for(icol=0;icol<ncolvar;++icol){
    f_hills[icol]=0.0;
  }
  for(ihill=0;ihill<nhills;++ihill){
    dp2=0.;
    for(icol=0;icol<ncolvar;++icol){      
      if(lhills[icol]) {
	dp[icol]=(pcolvar[icol]-t_pcolvar[icol][ihill])/t_hillwidth[icol][ihill];
	dp2 += dp[icol]*dp[icol];
      }
    }
    gauss = -t_hilldepth[ihill]*exp(-0.5*dp2);
    for(icol=0;icol<ncolvar;++icol){
      if(lhills[icol]) {
	f_hills[icol] -= gauss*dp[icol]/t_hillwidth[icol][ihill];
      }
    }
    ecolvar += gauss;
  }
  
  for(icol=0;icol<ncolvar;++icol){
    coords->colvar.fcolvar[icol] += f_hills[icol];
  }
  
  free(dp);
  return;
}
