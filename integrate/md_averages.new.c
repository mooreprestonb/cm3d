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

/* routines that set and accumulate the averages */

#include "md.h"
/*------------------------------------------------------------------*/
void set_econi_new(ENERGIES *energies)
{
  energies->ham = energies->econi =
    (energies->poth+energies->potra+energies->zke+
     energies->potn+energies->zken+
     energies->potv+energies->zkev);
  
#ifdef DEBUG
  printf("energies->poth  = %9g\n",energies->poth);
  printf("energies->potra = %9g\n",energies->potra);
  printf("energies->zke   = %9g\n",energies->zke);
  printf("energies->potn  = %9g\n",energies->potn);
  printf("energies->zken  = %9g\n",energies->zken);
  printf("energies->potv  = %9g\n",energies->potv);
  printf("energies->zkev  = %9g\n",energies->zkev);
#endif
}

/*------------------------------------------------------------------*/
void zero_averages_new(SIMPARMS *simparms,ENERGIES *energies)
{
  int i;
  if (simparms->istart == 1|| simparms->istart == 2){
    md_stdout("***********************************");
    md_stdout("* istart=1 or 2 : Reset averages  *");
    md_stdout("*      and counters and deriv.    *");
    md_stdout("***********************************");

    simparms-> istep  = 1;
    energies -> avham   = 0.0;
    energies -> avpotra = 0.0;
    energies -> avpoter = 0.0;
    energies -> avzke   = 0.0;
    energies -> avpotn  = 0.0;
    energies -> avzken  = 0.0;
    energies -> avpotv  = 0.0;
    energies -> avzkev  = 0.0;
    energies -> avpot_inter = 0.0;
    energies -> avpot_bond  = 0.0;
    energies -> avpot_bend  = 0.0;
    energies -> avpot_tors  = 0.0;
    energies -> avpot_onfo  = 0.0;
    energies -> avpot_onfo_e= 0.0;
    energies -> avpot_elec  = 0.0;
    energies -> avpot_recip = 0.0;
    energies -> avpot_extern = 0.0;
    energies -> avpot_extern_e = 0.0;
    energies -> avtemp = 0.0;
    energies -> avprs  = 0.0;
    for(i=0;i<9;i++) energies -> avprs_tensor[i] = 0.0;
    energies -> avvol  = 0.0;
    energies -> econv  = 0.0;
    energies -> econ2  = 0.0;
  } else if (simparms->istart == 3){
    md_stdout("****************************************");
    md_stdout("* istart=3: Continue from previous run *");
    md_stdout("****************************************");
  } else {
    md_error("istart must be 1,2 or 3");
  }
}
/*------------------------------------------------------------------*/

void accumulate_new(double *hmat,ENERGIES *energies)
{
  double vol;
  int i;

  vol = get_deth(hmat);
  energies->ham = (energies->poth+energies->potra+energies->zke+
		   energies->potn+energies->zken+
		   energies->potv+energies->zkev); 

  energies->avpotra += energies->potra;
  energies->avpoter += energies->poter;
  energies->avzke   += energies->zke;
  energies->avpotn  += energies->potn;
  energies->avzken  += energies->zken;
  energies->avpotv  += energies->potv;
  energies->avzkev  += energies->zkev;
  energies->avham   += energies->ham;
  
  energies->avprs   += energies->prsi;
  for(i=0;i<9;i++) energies->avprs_tensor[i]  += energies->prs_tensor[i];
  energies->avtemp  += energies->tiout;
  energies->avvol   += vol;
  
  energies->avpot_inter  += energies->pot_inter;
  energies->avpot_bond   += energies->pot_bond;
  energies->avpot_bend   += energies->pot_bend;
  energies->avpot_tors   += energies->pot_tors;
  energies->avpot_onfo   += energies->pot_onfo;
  energies->avpot_onfo_e += energies->pot_onfo_e;
  energies->avpot_elec   += energies->pot_elec;
  energies->avpot_recip  += energies->pot_recip;
  energies->avpot_extern += energies->pot_extern;
  energies->avpot_extern_e += energies->pot_extern_e;
  
  energies->econv += energies->ham;
  energies->econ2 += energies->ham*energies->ham;
  
  (energies->estep)++;
  energies->deltae2 += (((energies->econi-energies->ham)*
			 (energies->econi-energies->ham))/
			(energies->econi*energies->econi));
}
/*------------------------------------------------------------------*/
