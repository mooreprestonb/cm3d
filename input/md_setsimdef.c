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

#define DATA_BASE   "data_base"
#define VEC_CONF    "eigvec.his"
#define VEC_FILE    ""

void set_sim_default(SIMPARMS *simparms,WRITE_STEP *write_step,
		     SUBSPACE *subspace,FILENAMES *filenames,
		     INTER *inter,NGBR *ngbr,COORDS *coords,STARTUP *startup)
{
  simparms->icalc_type = 0;
  simparms->dt         = 1.0;
  simparms->temp       = 300.0;
  simparms->pext       = 0.0;
  simparms->iperd      = 3;
  simparms->iensemble  = 0;  /* nve */
  simparms->imin_type  = 0;  /* powell minimize structure */
  simparms->istart     = 1;
  simparms->nstep      = 1;
  simparms->ninter     = 1; /* short range steps */
  simparms->ninnra     = 1; /* intra steps */
  simparms->ninrra     = 1; /* torsion steps */
  simparms->nhc        = 1; /* nhc steps */
  simparms->nchain     = 0;
  simparms->ntherm     = 0;
  simparms->nbar       = 0;
  simparms->ivol       = 0;
  simparms->iunits     = 0;
  simparms->nlen       = 500; /* good vector length */
  simparms->npoints    = 0;
  simparms->nocell_all = 1; /* use cell to speed things up (usually) */
  simparms->max_exclude= 15; 
  simparms->iextern    = 0; /* no external potential */
  simparms->min_tol    = 1e-6;
  simparms->wallclock  = 1E12; /*make the time infinity*/
  simparms->vlrc       = 0;
  simparms->scaleeps   = 1.0;
  simparms->scalecharge = 1.0;
  simparms->scaletemp  = 1.0;
  simparms->nyoshida   = 3;
  simparms->ipolar     = 0;

  write_step->nscrn    = 0;
  write_step->ndump    = 0;
  write_step->nconf    = 0;
  write_step->nvel     = 0;
  write_step->nforce   = 0;
  write_step->ncolv    = 0;
  write_step->ninst    = 0;
  write_step->resamvel = 0;
  write_step->rescalevel = 0;
  write_step->nupdss   = 0;
  write_step->psysvec  = 0;
  write_step->peigval  = 0;

  subspace->stemp_update = 1; /* scale temp on updating subspace vectors */
  subspace->igetvec   = -1;   
  subspace->nstate    =  0;   /* number of total eigenvectors */
  subspace->num_subst =  1;   /* number of substructures */
  subspace->freq_min  = -500.0;
  subspace->freq_max  = 5000.0;
  strcpy(subspace->vecfile,VEC_FILE);
  strcpy(subspace->vecconf,VEC_CONF);

  ngbr->ilist       = 1; /* verlist */
  ngbr->ncells      = 7; 
  inter->skin_ter   = .5;
  inter->rheal      = .5;
  inter->ntable     = 500.; 

  coords->zmin_extern = 0;
  coords->zmax_extern = 50;

  coords->colvar.ncolvar = 0;
  coords->colvar.maxstephill = 1000;
  coords->colvar.minstephill = 100;
  coords->colvar.ltunehills = 0;

  startup->tau_nhc    = 1000;
  startup->tau_vol    = 1000;        
  startup->alp_ewald  = 0.15;
  startup->kmax       = 7.;           
  startup->ishift     = 0; 
  startup->scale_onfo = 1.;
  startup->scale_onfo_e = 1.;
  startup->rcute_max  = 10.;
  startup->rcute_min  = 1.;      
  startup->rcute_resp = 7.;

  strcpy(startup->inter_file,DATA_BASE);
  strcpy(startup->bond_file,DATA_BASE);
  strcpy(startup->bend_file,DATA_BASE);
  strcpy(startup->tors_file,DATA_BASE);
  strcpy(startup->onfo_file,DATA_BASE);
  strcpy(startup->nvecstore,"-1");
  strcpy(startup->nvecreal,"-1");
  strcpy(startup->nvecimag,"-1");

  strcpy(filenames->molvec,"data/mol_");
  strcpy(filenames->partrat,"partratio");
  strcpy(filenames->sysvec,"data/sys_");
  strcpy(filenames->nmfile,"dos.dat");
  strcpy(filenames->specfile,"spectra.dat");

  strcpy(filenames->configs,"");
  strcpy(filenames->velfile,"");
  strcpy(filenames->forcefile,"");
  strcpy(filenames->instham,"");
  strcpy(filenames->insteng,"");
  strcpy(filenames->instext,"");
  strcpy(filenames->extfile,"");
}


