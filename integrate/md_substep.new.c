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

void substep_init_new(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies,
		  SUBSPACE *subspace,INTER *inter,NGBR *ngbr)
{
  int i;

  md_stdout("Setting up usefull constants for the subspace integrator");

  zerotm(simparms->natoms,coords->vx,coords->vy,coords->vz,coords->amass);
  if(simparms->iperd==0){
    rotate(simparms->natoms,coords->px,coords->py,coords->pz,coords->amass);
    zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
	      coords->vx,coords->vy,coords->vz,coords->amass);
  }
  if(subspace->nstate_max==0){
    simparms->ndof = simparms->natoms*3;
  }

#ifdef DEBUG
  printf("SUBSPACE STUFF\n");
  printf("%d %d %d %d %d %d %d %d\n",
	 simparms->natoms,subspace->nstate,subspace->num_subst,
	 simparms->ninter,simparms->ninnra,simparms->ntherm,
	 simparms->nbar,simparms->ivol);
  for(i=0;i<subspace->num_subst;i++){
    printf("%d %d %d %d\n",i,subspace->num_vecsub_max[i],
	   subspace->num_real_max[i],subspace->num_imag_max[i]);
  }
#endif
  coords->xold = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  coords->yold = coords->xold+1*simparms->natoms;
  coords->zold = coords->xold+2*simparms->natoms;

  /* allocate maximum number of states possible */
  if(subspace->nstate_max == 0){ 
    subspace->nstate = simparms->natoms*3;
  } else {
    subspace->nstate = subspace->nstate_max;
  }

  /* allocate subspace memory */
  subspace->q0 = (double *)cmalloc(3*subspace->nstate*sizeof(double));
  subspace->q1 = subspace->q0+subspace->nstate;
  subspace->q2 = subspace->q1+subspace->nstate;
  subspace->d2v= NULL;

  md_stdout("Calculating initial forces");
  
  ngbr->update = -1;
  get_ngbr(simparms,coords,inter,ngbr);
#ifdef DEBUG
  check_distance(simparms,coords,inter);
#endif

  force(simparms,coords,inter,ngbr,1);
  force(simparms,coords,inter,ngbr,2);
  force(simparms,coords,inter,ngbr,3);
  getvireal(simparms,coords,energies,inter,ngbr,1);
  forceV(simparms,coords,energies);

  /* use the fact the memory is continous */
  for(i=0;i<3*simparms->natoms;i++){
    coords->fxt[i]=coords->fxl[i]+coords->fxr[i]+coords->fxa[i];
  }

  /* SET UP REDUCED SPACE USING EIGENVECTORS
     5/25/95 pbmbs if we are on the first step in new subspace run project
     positions velocities and accelerations
     set up the intial {x,y,z} displacement vectors */

  update_subspace(simparms,coords,inter,subspace,ngbr);

}
/*-----------------------------------------------------------------*/
void update_subspace_new(SIMPARMS *simparms,COORDS *coords,INTER *inter,
		     SUBSPACE *subspace,NGBR *ngbr)
{
  int i,n3;
  double *pos,*pold,ttemp,tiout;
  LINE line;

  /* use the fact the the memory is contiguous to store old positions*/
  n3 = simparms->natoms*3;
  pos = coords->px;  pold = coords->xold;
  for(i=0;i<n3;i++){ pold[i] = pos[i]; }

  /* recalculate forces for these positions 
     (don't really need to do this but might as well be safe) */
  force(simparms,coords,inter,ngbr,1);
  force(simparms,coords,inter,ngbr,2);
  force(simparms,coords,inter,ngbr,3);

  /* use the fact the memory is continous */
  for(i=0;i<3*simparms->natoms;i++){
    coords->fxt[i]=coords->fxs[i]+coords->fxl[i]+coords->fxa[i];
  }

  if(subspace->d2v!=NULL){
    free_dmatrix(subspace->d2v,0,simparms->natoms*3-1,0,subspace->nstate-1);
  }

  /* get kinetic energy of system NOW */
  tiout = get_temp(simparms->natoms,coords->vx,coords->vy,coords->vz,
		   coords->amass,simparms->ndof);

  sprintf(line,"Setting up subspace transform %d (3N) by %d",
	  simparms->natoms*3,subspace->nstate);
  md_stdout(line);
  get_subspace(simparms,coords,inter,coords->fxt,coords->fyt,coords->fzt,
	       subspace);
  simparms->ndof = subspace->nstate;

  if(subspace->stemp_update){
    /* expand velocities */
    expandztox(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
	       subspace->nstate,simparms->natoms);
    mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
		coords->amass,-1);
    
    /* set temp to current temp */
    ttemp = simparms->temp;  
    simparms->temp = tiout;  
    samvel(simparms,coords,0); /* rescale velocity and zero angular momentum */

    /* reproject velocities */
    mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
		coords->amass,1);
    projectxtoz(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
		subspace->nstate,simparms->natoms);
    expandztox(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
	       subspace->nstate,simparms->natoms);
    mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
		coords->amass,-1);

    simparms->temp = ttemp;  /* reset the desired temp */
  }
}
/*-----------------------------------------------------------------*/
void sub_step_new(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies,
	      INTER *inter,NGBR *ngbr,SUBSPACE *subspace)
{
  int i;
  double dt,dto2;
  double *pos,*vel,*acel,*pold;
  
  /* time constants */
  dt   = simparms->dt;
  dto2 = dt/2.;
  pos = subspace->q0;
  vel = subspace->q1;
  acel= subspace->q2;

  /* Predictor step */
  for(i=0;i<subspace->nstate;i++){vel[i] += dto2*acel[i];}

  for(i=0;i<subspace->nstate;i++){pos[i] += vel[i]*dt;}

  /* calculate actual positions to get forces */
  /* reform positions by adding initial position onto the displacement */ 

  expandztox(coords->px,coords->py,coords->pz,subspace->q0,subspace->d2v,
	     subspace->nstate,simparms->natoms);
  mass_weight(simparms->natoms,coords->px,coords->py,coords->pz,
	      coords->amass,-1);
  pos = coords->px;  pold = coords->xold;
  for(i=0;i<3*simparms->natoms;i++){ pos[i]+=pold[i]; }

  /* refesh neighbor list if nessesary */
  nrefsh(simparms,coords,inter,ngbr);

  /* calculate forces (-grad(v)) at predicted points */

  force(simparms,coords,inter,ngbr,1);
  force(simparms,coords,inter,ngbr,2);
  force(simparms,coords,inter,ngbr,3);
  getvireal(simparms,coords,energies,inter,ngbr,1);

  /* sum total force as short inter (fxs), long inter (fxl), and intra (fxa)*/
  for(i=0;i<3*simparms->natoms;i++){
    coords->fxt[i]=coords->fxl[i]+coords->fxr[i]+coords->fxa[i];
  }

  /* get forces on the q's  */
  mass_weight(simparms->natoms,coords->fxt,coords->fyt,coords->fzt,
	      coords->amass,1);
  projectxtoz(coords->fxt,coords->fyt,coords->fzt,subspace->q2,subspace->d2v,
	      subspace->nstate,simparms->natoms);

  /* final 1/2 step */
  for(i=0;i<subspace->nstate;i++){vel[i] += dto2*acel[i];}

  /* calculate actual velocities for energy and correlation functions
     may only want to do on occasion
     */

  expandztox(coords->vx,coords->vy,coords->vz,subspace->q1,subspace->d2v,
	     subspace->nstate,simparms->natoms);
  mass_weight(simparms->natoms,coords->vx,coords->vy,coords->vz,
	      coords->amass,-1);

}
/*----------------------------------------------------------------*/
