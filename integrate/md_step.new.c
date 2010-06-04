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

void step_init_new(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies,
	       INTER *inter,NGBR *ngbr)
{
  /*  solvent - solvent potential cutoff  */
  if(simparms->rank==0){
    md_stdout("Calculating initial forces");
  }
  ngbr->update = -1;

  get_ngbr(simparms,coords,inter,ngbr);

#ifdef DEBUG
  check_distance(simparms,coords,inter);
#endif

  force(simparms,coords,inter,ngbr,1);
  force(simparms,coords,inter,ngbr,2);
  force(simparms,coords,inter,ngbr,3);

  if(simparms->rank==0){
    md_stdout("Calculating initial potential energies");
  }

  getvireal(simparms,coords,energies,inter,ngbr,1);
  geteng(simparms,coords,energies);
  gettngbr(ngbr);
  if(simparms->ivol) forceV(simparms,coords,energies);
}
/*------------------------------------------------------*/

void step_new(SIMPARMS *simparms, COORDS *coords, ENERGIES *energies, INTER *inter, NGBR *ngbr) {

  int i,intr,intra,inrra,n3,ninter,ninter1,natoms;
  double dtl2,dts2,dta2,dtr2,dter,dter1;
  double *pos,*vel,*acel,*acel2;

  int ncolvar = coords->colvar.ncolvar;
  natoms = simparms->natoms;
  ninter = simparms->ninter;
  ninter1= ninter-1;
  dtl2 = .5*simparms->dt;
  dts2 = dtl2/(double)simparms->ninter;
  dtr2 = dts2/(double)simparms->ninrra;
  dta2 = dtr2/(double)simparms->ninnra;
  dter = (double)ninter*dts2;
  dter1= (double)(1-ninter)*dts2;

  n3 = 3*natoms;
  pos = coords->p;
  vel = coords->v;

  /* get half step */
  if(simparms->ivol) {
#ifdef ISO
    coords->vvol += dtl2*coords->fvol;
#endif
    for(i=0;i<9;i++) coords->vvol9[i] += dtl2*coords->fvol9[i];
  } else {
    coords->vvol = 0.;
    for(i=0;i<9;i++) coords->vvol9[i] = 0.;
  }
  
  for(intr=0;intr<ninter;intr++){
    /* long and short range intermolecular forces */
    if(intr==0){
      acel  = coords->fl;  acel2 = coords->fs;
      for(i=0;i<n3;i++) vel[i] += dter*acel[i]+dter1*acel2[i];
    } else {
      acel= coords->fs;
      for(i=0;i<n3;i++) vel[i] += dts2*acel[i];
    }

    /* intRa force (torsion one-four)*/
    for(inrra=0;inrra<simparms->ninrra;inrra++){
      acel= coords->fr;
      for(i=0;i<n3;i++) vel[i] += dtr2*acel[i];
    
      /* intra force (bonds bends xbonds) */
      for(intra=0;intra<simparms->ninnra;intra++){
	acel = coords->fa;

	for(i=0;i<n3;i++) vel[i] += dta2*acel[i];
	for(i=0;i<ncolvar;++i) coords->colvar.vcolvar[i] += dta2*coords->colvar.fcolvar[i];
	for(i=0;i<n3;i++) pos[i] += dta2*vel[i];
	
	for(i=0;i<ncolvar;++i) coords->colvar.pcolvar[i] += dta2*coords->colvar.vcolvar[i];
	/* Volume coupled equations */
	if(simparms->ivol) intbar(simparms,coords,dta2,0);
	if(simparms->ntherm) intnhc(simparms,coords,dta2);
	if(simparms->ivol) intbar(simparms,coords,dta2,1);
	
	for(i=0;i<ncolvar;++i) coords->colvar.pcolvar[i] += dta2*coords->colvar.vcolvar[i];
	for(i=0;i<n3;i++) pos[i] += dta2*vel[i];
	
	force(simparms,coords,inter,ngbr,2);	
	acel = coords->fa;
	for(i=0;i<n3;i++) vel[i] += dta2*acel[i];
	for(i=0;i<ncolvar;++i) coords->colvar.vcolvar[i] += dta2*coords->colvar.fcolvar[i];
      }
      force(simparms,coords,inter,ngbr,3);
      acel= coords->fr;
      for(i=0;i<n3;i++) vel[i] += dtr2*acel[i];
    }
    nrefsh(simparms,coords,inter,ngbr);
    
    if(intr==ninter1){
      force(simparms,coords,inter,ngbr,1);
      acel  = coords->fl;  acel2 = coords->fs;
      for(i=0;i<n3;i++) vel[i] += dter*acel[i]+dter1*acel2[i];
    } else {
      force(simparms,coords,inter,ngbr,0);
      acel= coords->fs;
      for(i=0;i<n3;i++) vel[i] += dts2*acel[i];
    }
  }
  
  /* long range force  */
  getvireal(simparms,coords,energies,inter,ngbr,1);
  if(simparms->ivol){
    forceV(simparms,coords,energies);
#ifdef ISO
    coords->vvol += dtl2*coords->fvol;
#endif
    for(i=0;i<9;i++) coords->vvol9[i] += dtl2*coords->fvol9[i];
  }
  gettngbr(ngbr);
}
/*-----------------------------------------------------------*/
