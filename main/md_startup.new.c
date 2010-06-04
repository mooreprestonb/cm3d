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

/* file that sets up the calculation */

#include "md.h"

void md_startup_new(FILENAMES *filenames,SIMPARMS *simparms,
		SUBSPACE *subspace,WRITE_STEP *write_step,
		INTER *inter,NGBR *ngbr,COORDS *coords)
{
  int i;
  STARTUP startup;
  SPECIES *spec_now,*spec_next;

  test_anint();
  /* set command and infile for argument or use default MDFILES */

  /* simulation data set up */
  read_sim_inputs(simparms,subspace,write_step,filenames,
		  inter,ngbr,coords,&startup);

  /* read in the set file which specifies molecule types */
  read_setfile(simparms,filenames->setfile,&startup.nspec,&startup.spec_root,
	       coords);

  /* read in the parameters */
  read_parmfile(simparms,coords,&startup,inter,ngbr);

  simparms->ndof = simparms->natoms*3; /* for alloc set ndof to 3*natoms */
  if(simparms->iensemble == 0 ){
    simparms->ndof = simparms->natoms*3-3;
    if(simparms->iperd==0) simparms->ndof = simparms->natoms*3-6;
    if(simparms->ndof<=0) simparms->ndof=1;
  }

  if(simparms->icalc_type == 0){
    set_baro(simparms,coords,startup.tau_vol);
    set_therm(simparms,coords,startup.nspec,startup.spec_root,startup.tau_nhc);
  }

  if(simparms->icalc_type == 1){
    set_baro(simparms,coords,startup.tau_vol);  
    set_therm(simparms,coords,startup.nspec,startup.spec_root,startup.tau_nhc);
    if(simparms->ntherm !=0 || simparms->nbar != 0){
      md_error("Subspace is only allowed in the nve ensemble");
    }
  }

  /* free spec_root children */
  spec_now = &startup.spec_root; spec_now = spec_now->next;
  for(i=0;i<startup.nspec;i++){
    spec_next = spec_now->next; free(spec_now); spec_now = spec_next;
  }

  if(simparms->icalc_type == 3){
    simparms->nchain  = 0;
    simparms->ntherm = 0;
    simparms->nbar = 1;
  }
  allocate_extend(simparms,coords);

  if(simparms->iextern == 1){
    allocate_exttab(simparms,coords,inter->ntable);
    read_extern_pot(filenames->extfile,simparms->ntypes,simparms->atom,
		    inter->ntable,coords);
  }
  if((simparms->icalc_type == 0 || simparms->icalc_type == 1 ||
      simparms->icalc_type == 4 ) &&
     simparms->rank==0){
    /* set the molecular names in md_io.c */
    openfiles(simparms,write_step,filenames);
  }
  if(simparms->icalc_type == 3){
    open_nma_in(filenames->configs);
    simparms->ntherm = 0;
  }


}

/*-------------------------------------------------------------------*/ 
