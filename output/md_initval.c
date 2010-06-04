/* 
   Copyright (C) Dr. Preston B. Moore

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
/* #define ENERGIES  */

/*---------------------------------------------------------*/
void output_initval(SIMPARMS *simparms,COORDS *coords,NGBR *ngbr,
		    ENERGIES *energies)
{
  double total,ratoms,tmass,volume;
  LINE line;

  ratoms = 1./(double)(simparms->natoms);
  total = (energies->poth+energies->potra+energies->zke+
	   energies->potn+energies->zken+energies->potv+energies->zkev);
  tmass = dsum1(simparms->natoms,coords->amass);
  volume = get_deth(coords->hmat);

  if (simparms->iunits==1){
    sprintf(line,"\n  *** initial report ***");
    md_stdout(line);
    sprintf(line,"number of pairs  = s=%d t=%d se=%d te=%d",
	    ngbr->nprs,ngbr->nprt,ngbr->nprse,ngbr->nprte);
    md_stdout(line);
    sprintf(line,"");
    md_stdout(line);
    sprintf(line,"econv            =     %15.6f",(total)/ECONV); 
    md_stdout(line);
    sprintf(line,"temperature      =     %15.6f",
	    energies->zke*ratoms*3./2.);
    md_stdout(line);
#ifdef ENERGIES
    sprintf(line,"part kin         =     %15.6f",energies->zke/(ECONV));
    md_stdout(line);
    sprintf(line,"part vter        =     %15.6f",energies->poth/(ECONV));
    md_stdout(line);
    sprintf(line,"part vtra        =     %15.6f",energies->potra/(ECONV));
    md_stdout(line);
    sprintf(line,"pressure         =     %15.6f",energies->prsi*PCONV);
    md_stdout(line);
    sprintf(line,"nhc ke           =     %15.6f",energies->zken*2.0/
	    (double)(simparms->nchain*(simparms->ntherm+1)));
    md_stdout(line);
    sprintf(line,"vol ke           =     %15.6f",energies->zkev*2.0);
    md_stdout(line);
#endif
    sprintf(line,"volume           =     %15.6f",
	    volume*LCONV*LCONV*LCONV);
    md_stdout(line);
  }else{
    sprintf(line,"\n  *** Initial Report ***");
    md_stdout(line);
    sprintf(line,"npairs = s=%d t=%d se=%d te=%d",
	    ngbr->nprs,ngbr->nprt,ngbr->nprse,ngbr->nprte);
    md_stdout(line);
#ifdef ENERGIES
    sprintf(line,"  Ura/N    = %g K",energies->potra*ratoms);
    md_stdout(line);
    sprintf(line,"  Uer/N    = %g K",energies->poter*ratoms);
    md_stdout(line);
    sprintf(line,"  K/N      = %g K",energies->zke*ratoms);
    md_stdout(line);
    if(simparms->ntherm){
      sprintf(line,"  U_eta/nc = %g K",energies->potn/
	      (double)(simparms->nchain*simparms->ntherm));
    md_stdout(line);
      sprintf(line,"  K_eta/nc = %g K",energies->zken/
	      (double)(simparms->nchain*simparms->ntherm));
    md_stdout(line);
    }
    if(simparms->ivol){
      sprintf(line,"  U_bar    = %g K",energies->potv/
	      ((double)simparms->nbar+1.));
    md_stdout(line);
      sprintf(line,"  K_bar    = %g K",energies->zkev/
	      ((double)simparms->nbar+1.));
    md_stdout(line);
    }
#endif
    sprintf(line,"  density  = %g g/cc",DCONV*tmass/volume);
    md_stdout(line);
    sprintf(line,"  Volume   = %g A^3 (cbrt = %g)",volume,cbrt(volume));
    md_stdout(line);
    sprintf(line,"  Temp     = %g K",energies->tiout);
    md_stdout(line);
    sprintf(line,"  Pressure = %g K/A^3 = %g ATM",energies->prsi,
	    energies->prsi*PCONV);
    md_stdout(line);
    sprintf(line,""); /* blank line */
    md_stdout(line);
  }
}
/*---------------------------------------------------------*/
