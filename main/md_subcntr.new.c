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


/* controling subroutine for subspace dynamcis 
   using nose thermostats (NVT) and extened system variables
   for the volume
   */

#include "md.h"
/* #define USER_FUNC */

void md_subcntr_new(char *command,FILENAMES *filenames,SIMPARMS *simparms,
		WRITE_STEP *write_step,COORDS *coords,
		SUBSPACE *subspace,INTER *inter,
		NGBR *ngbr)
{
  ENERGIES energies;

  /* read in the configuration */
  readpos(command,filenames->initfile,simparms,
	  &energies.estep,&energies,coords);

  /* calculate and setup initial tables and values */

  if(simparms->istart == 1){
    samvel(simparms,coords,1);
    if(simparms->iperd==0){
      rotate(simparms->natoms,coords->px,coords->py,coords->pz,
         coords->amass);
      zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
         coords->vx,coords->vy,coords->vz,coords->amass);
      zerotm(simparms->natoms,coords->vx,coords->vy,coords->vz,coords->amass);
      samvel(simparms,coords,0);
    }
    samveta(simparms,coords);
    if(simparms->ivol) samvbar(simparms,coords);
  }

  zero_averages(simparms,&energies);
  substep_init(simparms,coords,&energies,subspace,inter,ngbr);

  geteng(simparms,coords,&energies);

  if(simparms->iunits==0){
    getengn(simparms,coords,&energies.zken,&energies.potn);
    getengv(simparms,coords,&energies.zkev,&energies.potv);
  }

  /* output stuff the user might want to know */
  output_subparm(simparms,write_step,subspace,simparms->istep);
  
  output_initval(simparms,coords,ngbr,&energies);

  check_vals(9,energies.potra,energies.poter,energies.zke,
	     energies.potn,energies.zken,energies.potv,energies.zkev,
	     energies.tiout,energies.prsi);

  set_econi(&energies);

  /* initial energies */
  energies.acpu = 0.;
  simparms->fstep = simparms->istep;
  writescreen(simparms,ngbr,&energies,coords->hmat);

#ifdef USER_FUNC
  user_initial(simparms->istep,simparms->natoms,simparms->dt,
	       coords->px,coords->py,coords->pz,
	       coords->vx,coords->vy,coords->vz,
	       coords->amass,coords->qch,coords->hmat);
#endif

  /* 
  write_vectors(simparms->natoms,coords->px,coords->py,coords->pz,
		coords->amass,subspace,d2v,eigval,exp(coords->bar[0]));
		*/

  md_stdout("Begining Subspace Dynamics");
  
  /* **************  Commence Steps ******************** */
  energies.acpu = cputime();
  for(;simparms->istep<=simparms->nstep;simparms->istep++){

    sub_step(simparms,coords,&energies,inter,ngbr,subspace);
    geteng(simparms,coords,&energies);
    
    if(simparms->iunits==0){
      getengn(simparms,coords,&energies.zken,&energies.potn);
      getengv(simparms,coords,&energies.zkev,&energies.potv);
    } else {
      getengn_cheese(simparms,coords,&energies.zken,&energies.potn);
      getengv_cheese(simparms,coords,&energies.zkev,&energies.potv,
		     &energies.zken);
    }

    accumulate(coords->hmat,&energies);

    /* dump restart info and output averages every ndump passes */
    if(write_step->nscrn && !(simparms->istep%(write_step->nscrn))){
      energies.acpu += cputime();
      writescreen(simparms,ngbr,&energies,coords->hmat);
    }
    if(write_step->ndump && !(simparms->istep%(write_step->ndump))){
      writestart(filenames,simparms->istep,simparms,&energies,coords);
      check_distance(simparms,coords,inter);
    }

    if(write_step->rescalevel && !(simparms->istep%(write_step->rescalevel))){
      samvel(simparms,coords,0);
      if(simparms->iperd==0){
        zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
           coords->vx,coords->vy,coords->vz,coords->amass);
      }
      /* Set up new reduced coord velocities */
      mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
         coords->amass,1);
      projectxtoz(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
         subspace->nstate,simparms->natoms);
      expandztox(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
         subspace->nstate,simparms->natoms);
      mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
         coords->amass,-1);
    }

    if(write_step->resamvel && !(simparms->istep%(write_step->resamvel))){
      samvel(simparms,coords,1);
      if(simparms->iperd==0){
        zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
           coords->vx,coords->vy,coords->vz,coords->amass);
      }
      samveta(simparms,coords);
      if(simparms->ivol) samvbar(simparms,coords);
      mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
         coords->amass,1);
      projectxtoz(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
         subspace->nstate,simparms->natoms);
      expandztox(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
         subspace->nstate,simparms->natoms);
      mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
         coords->amass,-1);
    }
    
    if(write_step->nupdss && !(simparms->istep%write_step->nupdss)){
      update_subspace(simparms,coords,inter,subspace,ngbr);
    }

    /* -------  record instantaneous properties ------- */
    
    if(write_step->nconf && !(simparms->istep%(write_step->nconf)))
      save_conf(filenames->fconf,simparms->natoms,coords);
    
    if(write_step->nvel && !(simparms->istep%(write_step->nvel)))
      save_vel(filenames->fvel,simparms->natoms,coords);
    
    if(write_step->ninst && !(simparms->istep%(write_step->ninst)))
      save_inst(filenames->fham,filenames->feng,filenames->fext,
         simparms->istep,simparms->dt,&energies,coords->hmat);
    
#ifdef USER_FUNC
    user_function(simparms->istep,simparms->natoms,simparms->dt,
       coords->px,coords->py,coords->pz,
       coords->vx,coords->vy,coords->vz,
       coords->amass,coords->qch,coords->hmat);
#endif
  }
  md_stdout("Finnished Subspace Dynamics");
  finish(filenames,simparms,&energies,coords);
}

