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

void save_inst(FILE *fham,FILE *feng,FILE *fext,int istep,double dt,
	       ENERGIES *energies,double *hmat)
{
  int i;
  fprintf(fham,"%6g %.10g %g %g %g %g %g %g %g %g\n",
	  istep*dt,energies->ham,energies->poth,energies->potra,energies->zke,
	  energies->potn,energies->zken,energies->potv,energies->zkev,
	  energies->poter);
  fprintf(feng,"%6g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n",
	  istep*dt,energies->zke,energies->pot_inter,energies->pot_bond,energies->pot_bend,
	  energies->pot_tors,energies->pot_onfo,energies->pot_onfo_e,
	  energies->pot_elec,energies->pot_recip);
  fprintf(fext,"%6g %10g %10g %10g ",
	  istep*dt,energies->tiout,energies->prsi,get_deth(hmat));
  for(i=0;i<9;i++) fprintf(fext,"%g ",hmat[i]);
  for(i=0;i<9;i++) fprintf(fext,"%g ",energies->prs_tensor[i]);
  fprintf(fext,"\n");
  fflush(fham); fflush(feng); fflush(fext);
}
