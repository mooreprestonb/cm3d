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

/*---------------------------------------------------------*/

void openfiles(SIMPARMS *simparms,WRITE_STEP *write_step,
	       FILENAMES *filenames)
{
  LINE line;
  /* open file pointers */

  filenames->fconf= NULL;
  filenames->fvel = NULL;
  filenames->fham = NULL;
  filenames->feng = NULL;
  filenames->fext = NULL;
  filenames->fcolv= NULL;

  /* check for filenames */
  if(write_step->nconf>0){
    if(strlen(filenames->configs)==0) 
      md_error("No position file name (pos_file)");
  }
  if(write_step->nvel>0){
    if(strlen(filenames->velfile)==0) 
      md_error("No velocity file name (vel_file)");
  }
  if(write_step->nforce>0){
    if(strlen(filenames->forcefile)==0) 
      md_error("No force file name (force_file)");
  }
  if(write_step->ncolv>0){
    if(strlen(filenames->colvfile)==0) 
      md_error("No colvar file name (colv_file)");
  }
  if(write_step->ninst>0){
    if(strlen(filenames->insteng)==0)
      md_error("No Energy file name(eng_file)");
    if(strlen(filenames->instham)==0)
      md_error("No Hamiltonian file name (ham_file)");
    if(strlen(filenames->instext)==0) 
      md_error("No ext value file name(eval_file)");
  }
  if (simparms->istart!=3){
    if(write_step->nconf>0)filenames->fconf = cfopenw(filenames->configs);    
    if(write_step->nvel>0) filenames->fvel  = cfopenw(filenames->velfile);
    if(write_step->nforce>0)filenames->fforce  = cfopenw(filenames->forcefile);
    if(write_step->ncolv>0) filenames->fcolv  = cfopenw(filenames->colvfile);
    if(write_step->ninst>0){
      filenames->fham  = cfopenw(filenames->instham);    
      filenames->feng  = cfopenw(filenames->insteng);
      filenames->fext  = cfopenw(filenames->instext);
    }
  } else {
    if(write_step->nconf>0)filenames->fconf = fopen(filenames->configs,"a");
    if(write_step->nvel>0) filenames->fvel  = fopen(filenames->velfile,"a");
    if(write_step->nforce>0)filenames->fforce=fopen(filenames->forcefile,"a");
    if(write_step->ncolv>0) filenames->fcolv  = fopen(filenames->colvfile,"a");
    if(write_step->ninst>0){
      filenames->fham  = fopen(filenames->instham,"a+");
      filenames->feng  = fopen(filenames->insteng,"a+");
      filenames->fext  = fopen(filenames->instext,"a+");
    }
  }
  /* configuration header */
  if(write_step->nconf > 0){
    if(filenames->fconf == NULL){
      sprintf(line,"can't open %s (config) datafile",filenames->configs);
      md_error(line);
    } else {
      if (simparms->istart!=3)
        fprintf(filenames->fconf,"# %d %d %g\n",
           simparms->natoms,write_step->nconf,simparms->dt);
    }
  }

  /* velocity header */
  if(write_step->nvel>0){
    if(filenames->fvel  == NULL ){
      sprintf(line,"can't open %s (velocity) datafile",filenames->velfile);
      md_error(line);
    } else {
      if (simparms->istart!=3)
        fprintf(filenames->fvel ,"# %d %d %g\n",
           simparms->natoms,write_step->nvel,simparms->dt);
    }
  }
  
  /* force header */
  if(write_step->nforce>0){
    if(filenames->fforce  == NULL ){
      sprintf(line,"can't open %s (force) datafile",filenames->forcefile);
      md_error(line);
    } else {
      if (simparms->istart!=3)
        fprintf(filenames->fforce ,"# %d %d %g\n",
           simparms->natoms,write_step->nvel,simparms->dt);
    }
  }

  /* collective variables header */
  if(write_step->ncolv>0){
    if(filenames->fcolv  == NULL ){
      sprintf(line,"can't open %s (collective variable) datafile",filenames->colvfile);
      md_error(line);
    } else {
      if (simparms->istart!=3)
        fprintf(filenames->fcolv ,"# %d %d %g\n",
           simparms->natoms,write_step->ncolv,simparms->dt);
    }
  }

  /* instant headers */
  if(write_step->ninst > 0){
    if(filenames->fham == NULL){
      sprintf(line,"can't open %s (hamiltonian) datafile",filenames->instham);
      md_error(line);
    } else {
      if (simparms->istart!=3){
        fprintf(filenames->fham,"# %d %d %d %d %d %g\n",simparms->natoms,
           simparms->ntherm,simparms->nchain,simparms->nbar,
           write_step->ninst,simparms->dt);
        fprintf(filenames->fham,
           "# time ham inter intra ke potn ken potv kev inter\n");
      }
    }
    
    if(filenames->feng == NULL){
      sprintf(line,"can't open %s (instant energy) datafile",
         filenames->insteng);
      md_error(line);
    } else {
      if (simparms->istart!=3){
        fprintf(filenames->feng,"# %d %d %9g\n",
           simparms->natoms,write_step->ninst,simparms->dt);
        fprintf(filenames->feng,
           "# time Ke inter bond bend tors onfo onfoe elec recip temp prs\n");
      }
    }
    if(filenames->fext == NULL){
      sprintf(line,"can't open %s (instant values) datafile",
         filenames->instext);
      md_error(line);
    } else {
      if (simparms->istart!=3){
        fprintf(filenames->fext,"# %d %d %9g\n",
           simparms->natoms,write_step->ninst,simparms->dt);
        fprintf(filenames->fext,
           "# time temp prs vol hmat pmat\n");
      }
    }
  }
}

/*---------------------------------------------------------*/
