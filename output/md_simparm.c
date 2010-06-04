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
void output_simparm(SIMPARMS *simparms,WRITE_STEP *write_step,int istep)
{
  LINE line;
  /* output stuff that the user might want to know */
  sprintf(line,"\ntype of calculation   = %d",simparms->istart);
  md_stdout(line);
  sprintf(line,"units                 = %d",simparms->iunits);
  md_stdout(line);
  sprintf(line,"periodicity           = %d",simparms->iperd);
  md_stdout(line);
  sprintf(line,"# of atoms            = %d",simparms->natoms);
  md_stdout(line);
  sprintf(line,"# of frozen atoms     = %d",simparms->nfreeze);
  md_stdout(line);
  sprintf(line,"MD time step          = %g ps",simparms->dt);
  md_stdout(line);
  sprintf(line,"# of inter time steps = %d (dt = %g)",
	  simparms->ninter,simparms->dt/simparms->ninter);
  md_stdout(line);
  sprintf(line,"# of Torsion steps    = %d (dt = %g)",simparms->ninrra,
	  simparms->dt/(simparms->ninter*simparms->ninrra));
  md_stdout(line);
  sprintf(line,"# of Bond/Bend steps  = %d (dt = %g)",simparms->ninnra,
	  simparms->dt/(simparms->ninter*simparms->ninrra*simparms->ninnra));
  md_stdout(line);
  sprintf(line,"# of thermostats      = %d",simparms->ntherm);
  md_stdout(line);
  sprintf(line,"# of nh chains links  = %d",simparms->nchain);
  md_stdout(line);
  sprintf(line,"# of barostats        = %d",simparms->nbar);
  md_stdout(line);
  if(simparms->ntherm){
    sprintf(line,"Temperature           = %g K",simparms->temp);
    md_stdout(line);
  } else {
    sprintf(line,"Temperature is not thermostated");
    md_stdout(line);
  }
  if(simparms->ivol){
    sprintf(line,"Pressure              = %g K/A^3",simparms->pext);
    md_stdout(line);
  } else {
    sprintf(line,"Volume will not fluctuate.");
    md_stdout(line);
  }

  sprintf(line,"\nAverages reported every       %d ts",write_step->nscrn);
  md_stdout(line);
  sprintf(line,"Restart info dumped every     %d ts",write_step->ndump);
  md_stdout(line);
  sprintf(line,"A config is save every        %d ts",write_step->nconf);
  md_stdout(line);
  sprintf(line,"A velocity file is save every %d ts",write_step->nvel);
  md_stdout(line);
  sprintf(line,"A force file is save every    %d ts",write_step->nforce);
  md_stdout(line);
  sprintf(line,"Instant quantaties recorded   %d ts",write_step->ninst);
  md_stdout(line);
  sprintf(line,"");
  md_stdout(line);
  sprintf(line,"Calcutaion time step %d to %d.",istep,simparms->nstep);
  md_stdout(line);
  sprintf(line,"total time of simulation %g ps.\n",
	  simparms->nstep*simparms->dt);
  md_stdout(line);
}
/*---------------------------------------------------------*/ 
