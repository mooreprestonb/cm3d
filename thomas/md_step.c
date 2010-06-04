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

/*------------------------------------------------------*/

void step_init(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies,
	       INTER *inter,NGBR *ngbr)
{
  char line[MAXLINELEN];

  /*  solvent - solvent potential cutoff  */
  if(simparms->rank==0){
    md_stdout("Calculating initial forces");
  }
  ngbr->update = -1;

  get_ngbr(simparms,coords,inter,ngbr);

#ifdef DEBUG
  check_distance(simparms,coords,inter);
#endif

  force(simparms,coords,inter,ngbr,0);
  force(simparms,coords,inter,ngbr,1);
  force(simparms,coords,inter,ngbr,2);

  if(simparms->rank==0){
    md_stdout("Calculating initial potential energies");
  }

  getvireal(simparms,coords,energies,inter,ngbr);
  geteng(simparms,coords,energies);
  gettngbr(ngbr);
  if(simparms->ivol) forceV(simparms,coords,energies);
}
/*------------------------------------------------------*/

void step(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies,
	  INTER *inter,NGBR *ngbr)
{
  int i,intr,intra,n3,k,iyosh;
  double dtl,dtl2,dts,dts2,dta2,dta,dtlnc2;
  double *pos,*vel,*acel;

  dtl  = simparms->dt;
  dts  = simparms->dt/(double)simparms->ninner;
  dta  = dts/(double)simparms->ninnra;
  dtl2 = (.5*dtl);
  dts2 = (.5*dts);
  dta2 = (.5*dta);

  n3 = simparms->natoms*3;
  pos = coords->px;
  vel = coords->vx;

  /* long range force  */
  acel= coords->fxl;
  if (simparms->iensemble==3) for(i=0;i<n3;i++) vel[i] = dtl2*acel[i];
  else for(i=0;i<n3;i++) vel[i] += dtl2*acel[i];
  
  /* short range force   */
  for(intr=0;intr<simparms->ninner;intr++){

    acel= coords->fxs;

    if (simparms->iensemble==3) for(i=0;i<n3;i++) vel[i] = dts2*acel[i];
    else for(i=0;i<n3;i++) vel[i] += dts2*acel[i];
  /* get half step */
  getvireal(simparms,coords,energies,inter,ngbr);
  if(simparms->ivol) {
   for (k=0;k<simparms->nhc;k++) {
   for (iyosh=0;iyosh<simparms->nyoshida;iyosh++) {
    forceV(simparms,coords,energies);
       dtlnc2=simparms->wyosh[iyosh]*dts2/((float) simparms->nhc);

#ifdef ISO
    coords->vvol += dtlnc2*coords->fvol;
#endif
    for(i=0;i<9;i++) coords->vvol9[i] += dtlnc2*coords->fvol9[i];
  }
  }
  } else {
    coords->vvol = 0.;
    for(i=0;i<9;i++) coords->vvol9[i] = 0.;
  }
    
    /* intra force */
    for(intra=0;intra<simparms->ninnra;intra++){
      acel = coords->fxa;
      if (simparms->iensemble==3) for(i=0;i<n3;i++) vel[i] = dta2*acel[i];
      else for(i=0;i<n3;i++) vel[i] += dta2*acel[i];
      
      for(i=0;i<n3;i++) pos[i] += dta2*vel[i];
      

      /* Volume coupled equations */

  
        if(simparms->ivol) intbar(simparms,coords,dta2,0);
        if(simparms->ntherm) intnhc(simparms,coords,dta2);
        if(simparms->ivol) intbar(simparms,coords,dta2,1);



      for(i=0;i<n3;i++) pos[i] += dta2*vel[i];
      
      force(simparms,coords,inter,ngbr,2);
      
      acel = coords->fxa;
      if (simparms->iensemble==3) for(i=0;i<n3;i++) vel[i] = dta2*acel[i]; 
      else for(i=0;i<n3;i++) vel[i] += dta2*acel[i];
    }

    nrefsh(simparms,coords,inter,ngbr);
    force(simparms,coords,inter,ngbr,0);
    
  getvireal(simparms,coords,energies,inter,ngbr);
  if(simparms->ivol){
   for (k=0;k<simparms->nhc;k++) {
   for (iyosh=0;iyosh<simparms->nyoshida;iyosh++) {
    forceV(simparms,coords,energies);
       dtlnc2=simparms->wyosh[iyosh]*dts2/((float) simparms->nhc);

#ifdef ISO
    coords->vvol += dtlnc2*coords->fvol;
#endif
    for(i=0;i<9;i++) coords->vvol9[i] += dtlnc2*coords->fvol9[i];
  }
  } 
  }
    acel = coords->fxs;
    if (simparms->iensemble==3) for(i=0;i<n3;i++) vel[i] = dts2*acel[i]; 
    else for(i=0;i<n3;i++)  vel[i] += dts2*acel[i];
  }
  
  /* long range force  */
  force(simparms,coords,inter,ngbr,1);

  acel = coords->fxl;
  if (simparms->iensemble==3) for(i=0;i<n3;i++) vel[i] = dtl2*acel[i];
  else for(i=0;i<n3;i++) vel[i] += dtl2*acel[i];
  gettngbr(ngbr);
}
/*-----------------------------------------------------------*/
