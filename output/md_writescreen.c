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

#define REPORT_ENERGIES
#define PRESSURE
/*--------------------------------------------------------*/
void writescreen(SIMPARMS *simparms,NGBR *ngbr,
		 ENERGIES *energies,double *hmat)
{     
  double update,vol,*apt,*pt;
  double stddev,unit_conv,step,ratoms;
  WORD sunit;
  LINE line;

  vol = get_deth(hmat);
  if(simparms->istep==0) { 
    step = 1.;
  } else {
    step = (double)simparms->istep;
  }
  ratoms = 1./(double)simparms->natoms;

  if(ngbr->update==0) {
    update=1.;
  } else {
    update = (double)ngbr->update;
  }
  
  sprintf(line,"time steps= %d, time= %g ps, <cpu>= %g, <cpu/ps>= %g",
	  simparms->istep,step*simparms->dt,
	  energies->acpu/
	  ((step-simparms->fstep)==0?1.:(step-simparms->fstep)),
	  energies->acpu/
	  (((step-simparms->fstep)==0?1.:(step-simparms->fstep))
	   *simparms->dt));
  md_stdout(line);
  sprintf(line,"nprs=%d, nprt=%d, nprse=%d, nprte=%d, <ts/update>=%g",
	  ngbr->nprst,ngbr->nprtt,ngbr->nprset,ngbr->nprtet,
	  (step-simparms->fstep)/update);
  md_stdout(line);

  switch(simparms->iunits) {
    /*-----------------------------------------------------*/
     default : case 0: /* Kelvin Angstroms picoseconds */
       unit_conv = 1.;  
       strcpy(sunit,"K"); 
       break;
     case 1: /* AU */
       unit_conv = 1./ECONV;
       strcpy(sunit,"AU"); 
       break;
     case 2:   /* Kcal Units */
       unit_conv = 1./KCAL;
       strcpy(sunit,"KCAL/mol"); 
       break;
  }

  if(energies->potra){
    if(isnan(unit_conv*energies->avpotra*ratoms/step) != 0 ||
       isnan(unit_conv*energies->potra*ratoms) != 0){
      sprintf(line,"###NAN ERROR###<Ura/N>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpotra*ratoms/step,
	      unit_conv*energies->potra*ratoms,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<Ura/N>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpotra*ratoms/step,
	      unit_conv*energies->potra*ratoms,sunit);
      md_stdout(line);
    }
  }
  if(energies->poter){
    if(isnan(unit_conv*energies->avpoter*ratoms/step) != 0 ||
       isnan(unit_conv*energies->poter*ratoms) !=0){
      sprintf(line,"###NAN ERROR###<Uer/N>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpoter*ratoms/step,
	      unit_conv*energies->poter*ratoms,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<Uer/N>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpoter*ratoms/step,
	      unit_conv*energies->poter*ratoms,sunit);
      md_stdout(line);
    }
  }
    if(energies->zke){
      if(isnan(unit_conv*energies->avzke*ratoms/step) != 0 ||
	 isnan(unit_conv*energies->zke*ratoms)!=0){
	sprintf(line,"###NAN ERROR###<Ke /N>  \t= %15.6f\t~ %15.6f %s",
		unit_conv*energies->avzke*ratoms/step,
		unit_conv*energies->zke*ratoms,sunit);
      md_error(line);
	exit(99);
      }else{
	sprintf(line,"<Ke /N>  \t= %15.6f\t~ %15.6f %s",
		unit_conv*energies->avzke*ratoms/step,
		unit_conv*energies->zke*ratoms,sunit);
	md_stdout(line);
      }
    }
#ifdef REPORT_ENERGIES
  if(energies->pot_extern){
    if(isnan(unit_conv*energies->avpot_extern/step) != 0 ||
       isnan(unit_conv*energies->pot_extern) != 0){
      sprintf(line,"###NAN ERROR###<extern>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_extern/step,
	      unit_conv*energies->pot_extern,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<extern>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_extern/step,
	      unit_conv*energies->pot_extern,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_extern_e){
    if(isnan(unit_conv*energies->avpot_extern_e/step) != 0 ||
       isnan(unit_conv*energies->pot_extern_e) != 0){
      sprintf(line,"###NAN ERROR###<extern_e>\t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_extern_e/step,
	      unit_conv*energies->pot_extern_e,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<extern_e>\t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_extern_e/step,
	      unit_conv*energies->pot_extern_e,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_inter){
    if(isnan(unit_conv*energies->avpot_inter/step) != 0 ||
       isnan(unit_conv*energies->pot_inter) !=0){
      sprintf(line,"###NAN ERROR <inter>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_inter/step,
	      unit_conv*energies->pot_inter,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<inter>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_inter/step,
	      unit_conv*energies->pot_inter,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_bond){
    if(isnan(unit_conv*energies->avpot_bond/step) != 0 ||
       isnan(unit_conv*energies->pot_bond) !=0){
      sprintf(line,"###NAN ERROR <bond>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_bond/step,
	      unit_conv*energies->pot_bond,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<bond>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_bond/step,
	      unit_conv*energies->pot_bond,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_bend){
    if(isnan(unit_conv*energies->avpot_bend/step) != 0 ||
       isnan(unit_conv*energies->pot_bend) !=0){
      sprintf(line,"###NAN ERROR <bend>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_bend/step,
	      unit_conv*energies->pot_bend,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<bend>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_bend/step,
	      unit_conv*energies->pot_bend,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_tors){
    if(isnan(unit_conv*energies->avpot_tors/step) != 0 ||
       isnan(unit_conv*energies->pot_tors) !=0){
      sprintf(line,"###NAN ERROR <tors>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_tors/step,
	      unit_conv*energies->pot_tors,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<tors>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_tors/step,
	      unit_conv*energies->pot_tors,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_onfo){
    if(isnan(unit_conv*energies->avpot_onfo/step) != 0 ||
       isnan(unit_conv*energies->pot_onfo) !=0){
      sprintf(line,"###NAN ERROR <onefour>\t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_onfo/step,
	      unit_conv*energies->pot_onfo,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<onefour>\t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_onfo/step,
	      unit_conv*energies->pot_onfo,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_onfo_e){
    if(isnan(unit_conv*energies->avpot_onfo_e/step) != 0 ||
       isnan(unit_conv*energies->pot_onfo_e) !=0){
      sprintf(line,"###NAN ERROR <onfo_e> \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_onfo_e/step,
	      unit_conv*energies->pot_onfo_e,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<onfo_e> \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_onfo_e/step,
	      unit_conv*energies->pot_onfo_e,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_elec){
    if(isnan(unit_conv*energies->avpot_elec/step) != 0 ||
       isnan(unit_conv*energies->pot_elec) !=0){
      sprintf(line,"###NAN ERROR <elec>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_elec/step,
	      unit_conv*energies->pot_elec,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<elec>   \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_elec/step,
	      unit_conv*energies->pot_elec,sunit);
      md_stdout(line);
    }
  }
  if(energies->pot_recip){
    if(isnan(unit_conv*energies->avpot_recip/step) != 0 ||
       isnan(unit_conv*energies->pot_recip) !=0){
      sprintf(line,"###NAN ERROR <recip>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_recip/step,
	      unit_conv*energies->pot_recip,sunit);
      md_error(line);
      exit(99);
    }else{
      sprintf(line,"<recip>  \t= %15.6f\t~ %15.6f %s",
	      unit_conv*energies->avpot_recip/step,
	      unit_conv*energies->pot_recip,sunit);
      md_stdout(line);
    }
  }
#endif
  if(simparms->ntherm){
    sprintf(line,"<Ueta/NC>\t= %15.6f\t~ %15.6f %s",
	    unit_conv*energies->avpotn/
	    (step*(double)simparms->nchain*(double)simparms->ntherm),
	    unit_conv*energies->potn/
	    ((double)simparms->nchain*(double)simparms->ntherm),sunit);
    md_stdout(line);
    sprintf(line,"<Keta/NC>\t= %15.6f\t~ %15.6f %s",
	    unit_conv*energies->avzken*2./
	    (step*(double)(simparms->nchain*(simparms->ntherm+1))),
	    unit_conv*energies->zken*2./
	    ((double)simparms->nchain*((double)simparms->ntherm+1.)),sunit);
    md_stdout(line);
  }
  if(simparms->ivol){
    sprintf(line,"<Ubar/NB>\t= %15.6f\t~ %15.6f %s",
	    unit_conv*energies->avpotv/(step*(double)(simparms->nbar+1)),
	    unit_conv*energies->potv/(double)(simparms->nbar+1),sunit);
    md_stdout(line);
    sprintf(line,"<Kbar/NB>\t= %15.6f\t~ %15.6f %s",
	    unit_conv*energies->avzkev/(step*(double)(simparms->nbar+1)),
	    unit_conv*energies->zkev/(double)(simparms->nbar+1),sunit);
    md_stdout(line);
    sprintf(line,"<vol>    \t= %15.6f\t~ %15.6f A^3 ",
	    energies->avvol/step,vol);
    md_stdout(line);
    sprintf(line,"\t cbrt(vol)=%15.6f A",cbrt(vol));
    md_stdout(line);
  }
  sprintf(line,"<Ham/N>  \t= %15.6f\t~ %15.6f %s",
	  unit_conv*energies->avham*ratoms/step,
	  unit_conv*energies->ham*ratoms,sunit);
  md_stdout(line);
  sprintf(line,"<temp>   \t= %15.6f\t~ %15.6f K",
	  energies->avtemp/step,energies->tiout);
  md_stdout(line);    
  sprintf(line,"<pres>   \t= %15.6f\t~ %15.6f %s/A^3",
	  unit_conv*energies->avprs/step,
	  unit_conv*energies->prsi,sunit);
  md_stdout(line);    
  sprintf(line,"<pres>   \t= %15.6f\t~ %15.6f atm",
	  energies->avprs*PCONV/(step),energies->prsi*PCONV);
  md_stdout(line);

#ifdef PRESSURE
  pt = energies->prs_tensor;
  apt = energies->avprs_tensor;
  
  sprintf(line,"Pressure Tensor");
  md_stdout(line);
  sprintf(line,"%10f %10f %10f ~ %10f %10f %10f atm",
     apt[0]*PCONV/(step), apt[1]*PCONV/(step), apt[2]*PCONV/(step),
	  pt[0]*PCONV,	  pt[1]*PCONV,	  pt[2]*PCONV);
  md_stdout(line);
  sprintf(line,"%10f %10f %10f ~ %10f %10f %10f atm",
     apt[3]*PCONV/(step),     apt[4]*PCONV/(step),     apt[5]*PCONV/(step),
	  pt[3]*PCONV,	  pt[4]*PCONV,	  pt[5]*PCONV);
  md_stdout(line);
  sprintf(line,"%10f %10f %10f ~ %10f %10f %10f atm",
     apt[6]*PCONV/(step),     apt[7]*PCONV/(step),     apt[8]*PCONV/(step),
	  pt[6]*PCONV,	  pt[7]*PCONV,	  pt[8]*PCONV);
  md_stdout(line);
#endif
  
#ifdef PRESSURE_1
  sprintf(line,"<pres xx>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[0]*PCONV/(step),
	  energies->prs_tensor[0]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres xy>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[1]*PCONV/(step),
	  energies->prs_tensor[1]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres xz>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[2]*PCONV/(step),
	  energies->prs_tensor[2]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres yx>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[3]*PCONV/(step),
	  energies->prs_tensor[3]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres yy>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[4]*PCONV/(step),
     energies->prs_tensor[4]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres yz>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[5]*PCONV/(step),
	  energies->prs_tensor[5]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres zx>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[6]*PCONV/(step),
	  energies->prs_tensor[6]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres zy>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[7]*PCONV/(step),
	  energies->prs_tensor[7]*PCONV);
  md_stdout(line);
  sprintf(line,"<pres zz>   \t= %15.6f\t~ %15.6f atm",
     energies->avprs_tensor[8]*PCONV/(step),
	  energies->prs_tensor[8]*PCONV);
  md_stdout(line);
#endif
  
  if(simparms->istep>1 && energies->estep!=0.0){
    stddev = fabs((step*energies->econ2-energies->econv*energies->econv)/
		  ((step-1.)*step));
    sprintf(line,"<dE>     \t= %15.6f\t~%15.6f",
	    sqrt(stddev)*step/fabs(energies->econv),
	    sqrt(energies->deltae2)/(double)(energies->estep));
  md_stdout(line);
  }
#ifdef GCODE
  if(simparms->ntherm){
    sprintf(line,"NCHAIN KE   \t= %15.6f\t~ %15.6f K",
	    energies->avzken*2./(step*(double)simparms->nchain*
				 ((double)simparms->ntherm+1.)),
	    energies->zken*2./((double)simparms->nchain*
			       ((double)simparms->ntherm+1.)));
  md_stdout(line);
    sprintf(line,"NCHAIN PE   \t= %15.6f\t~ %15.6f au",
	    energies->avpotn/(step*ECONV),energies->potn/ECONV);
  md_stdout(line);
  }
  if(simparms->ivol){
    sprintf(line,"VOL KE   \t= %15.6f\t~ %15.6f K",
	    2.0*energies->avzkev/step,2.0*energies->zkev);
  md_stdout(line);
    sprintf(line,"VOL PE   \t= %15.6f\t~ %15.6f au",
	    energies->avpotv/(step*ECONV),energies->potv/ECONV);
  md_stdout(line);
    sprintf(line,"VOLUME   \t= %15.6f\t~ %15.6f A^3",
	    energies->avvol*LCONV3/step,vol*LCONV3);
  md_stdout(line);
  }
#endif
  fflush(stdout); fflush(stderr);
}
