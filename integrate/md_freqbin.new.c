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

#define ZERO  /* remove zero frequency eigenvalues */
#define CHOP

/* #define DEBUG */

/* #define STEP */

/*------------------------------------------------------------------*/
void freqbin_new(SIMPARMS *simparms,int npoints,double freq_min,NMODES *nmodes)
{
  /* create histogram  */
  int i,j,k,zero[3],np;
  double fr,fr2,factor,kexp;
  LINE line;

  np = simparms->natoms*3;
  kexp = HBAR/(BOLTS*simparms->temp*FCONV);

  zero[0] = zero[1] = zero[2] = -1;
#ifdef ZERO
  fr = DBL_MAX;
  for(j=0;j<np;j++){
    if(fabs(nmodes->dn[j]) < fr){zero[0] = j;  fr = fabs(nmodes->dn[j]); }
  }
  fr = DBL_MAX;
  for(j=0;j<np;j++){
    if(fabs(nmodes->dn[j]) < fr && j != zero[0]){
      zero[1] = j;  
      fr = fabs(nmodes->dn[j]);
    }
  }
  fr = DBL_MAX;
  for(j=0;j<np;j++){
    if(fabs(nmodes->dn[j])<fr&&j!=zero[0]&&j!=zero[1]){
      zero[2] = j;
      fr=fabs(nmodes->dn[j]);
    }
  }
  if(nmodes->dn[zero[0]]>ERRMAX || nmodes->dn[zero[1]] > ERRMAX || 
     nmodes->dn[zero[2]] > ERRMAX){
    sprintf(line,"zero frequency exceeds tolerence of %g\n",ERRMAX);
    md_warning(line);
    sprintf(line,"%d %d %d %g %g %g\n",zero[0],zero[1],zero[2],
	    nmodes->dn[zero[0]],nmodes->dn[zero[1]],nmodes->dn[zero[2]]);
    md_warning(line);
  }
#endif

  for(i=0;i<np;i++){
#ifdef DEBUG
    printf("eigenvalue %d = %g\n",i,nmodes->dn[i]);
#endif
    if(i != zero[0] && i != zero[1] && i != zero[2]){
      fr=nmodes->dn[i];
      if (fr<0.0) fr= -sqrt(-fr);
      else fr=sqrt(fr);
      nmodes->dn[i] = fr;
    } else {
      nmodes->dn[i] = 0.;
    }
  }

  /*  convert eigenvalues in ps^-1 to cm-1 if you want before recording */
  for (i=0;i<np;i++){
    if(simparms->iunits==3) nmodes->dn[i] *= FCONV;
    if(nmodes->dn[i] < nmodes->minfreq) nmodes->minfreq = nmodes->dn[i];
    if(nmodes->dn[i] > nmodes->maxfreq) nmodes->maxfreq = nmodes->dn[i];
  }

  for(j=0;j<np;j++){
    if(j != zero[0] && j != zero[1] && j != zero[2]){
      fr=nmodes->dn[j];
      fr=(fr-freq_min)/nmodes->df;
      k=(int)floor(fr);
      k=MAX(k,0);
      k=MIN(k,npoints-2);
      fr -= k;

#ifdef CHOP
      if(k==0||k==(npoints-2)) fr = 0.;
      fr = MIN(fr,1.);
      fr = MAX(fr,0.);
#endif
      /* DOS ONLY */

#ifdef STEP
      nmodes->freq[k]++;
      for(i=0;i<nmodes->imodes*nmodes->nspec;i++)
        nmodes->decmod[i][k] += nmodes->permod[i][j];
#else
      nmodes->freq[k]   += 1.-fr;
      nmodes->freq[k+1] += fr;
      for(i=0;i<nmodes->imodes*nmodes->nspec;i++){
        nmodes->decmod[i][k]   += (1.-fr)*nmodes->permod[i][j];
        nmodes->decmod[i][k+1] +=    (fr)*nmodes->permod[i][j];
      }
#endif
      

      /* IR */
      if(simparms->nma_type == 0 || simparms->nma_type == 2){  
#ifdef STEP
        nmodes->fric[k] += nmodes->fricmod[j];
        nmodes->fpos[k] += nmodes->fricmodq[j];
        nmodes->ftot[k] += nmodes->fricmodcq[j];
#else
        nmodes->fric[k]   += (1.-fr)*nmodes->fricmod[j];
        nmodes->fric[k+1] += (fr)*nmodes->fricmod[j];
        nmodes->fpos[k]   += (1.-fr)*nmodes->fricmodq[j];
        nmodes->fpos[k+1] += (fr)*nmodes->fricmodq[j];
        nmodes->ftot[k]   += (1.-fr)*nmodes->fricmodcq[j];
        nmodes->ftot[k+1] += (fr)*nmodes->fricmodcq[j];
        nmodes->fpart[k]   += (1.-fr)*nmodes->part[j];
        nmodes->fpart[k+1] += (fr)*nmodes->part[j];
        
#endif
        for(i=0;i<nmodes->imodes*nmodes->nspec;i++){
#ifdef STEP
          nmodes->decmodd[i][k] += nmodes->permod[i][j]*nmodes->fricmod[j];
          nmodes->decmodq[i][k] += nmodes->permod[i][j]*nmodes->fricmodq[j];
          nmodes->decmoddq[i][k] += nmodes->permod[i][j]*nmodes->fricmodcq[j];
#else
          nmodes->decmodd[i][k]   += (1.-fr)*nmodes->permod[i][j]*nmodes->fricmod[j];
          nmodes->decmodd[i][k+1] +=    (fr)*nmodes->permod[i][j]*nmodes->fricmod[j];
          nmodes->decmodq[i][k]   += (1.-fr)*nmodes->permod[i][j]*nmodes->fricmodq[j];
          nmodes->decmodq[i][k+1] +=    (fr)*nmodes->permod[i][j]*nmodes->fricmodq[j];
          nmodes->decmoddq[i][k]   += (1.-fr)*nmodes->permod[i][j]*nmodes->fricmodcq[j];
          nmodes->decmoddq[i][k+1] +=    (fr)*nmodes->permod[i][j]*nmodes->fricmodcq[j];
#endif
        }
      }
      /* RAMAN */
      if(simparms->nma_type == 1 || simparms->nma_type == 2){
        fr2 = nmodes->dn[j];
        if (fr2<0.0)fr2 = -fr2;
        factor = 1./(fr2* (1. - exp(-kexp*fr2)));
#ifdef STEP
        nmodes->redram[k] += nmodes->aniso_raman[j];
        nmodes->raman[k] += nmodes->aniso_raman[j]*factor;
        nmodes->oke[k] += nmodes->aniso_raman[j]/fr2;
        nmodes->iso[k] += nmodes->iso_raman[j];
#else
        nmodes->redram[k]   += (1.-fr)*nmodes->aniso_raman[j];
        nmodes->redram[k+1] += (fr)*nmodes->aniso_raman[j];
        nmodes->raman[k]   += (1.-fr)*nmodes->aniso_raman[j]*factor;
        nmodes->raman[k+1] += (fr)*nmodes->aniso_raman[j]*factor;
        nmodes->oke[k]   += (1.-fr)*nmodes->aniso_raman[j]/fr2;
        nmodes->oke[k+1] += (fr)*nmodes->aniso_raman[j]/fr2;
        nmodes->iso[k] += (1.-fr)*nmodes->iso_raman[j];
        nmodes->iso[k+1] += (fr)*nmodes->iso_raman[j];
#endif
        for(i=0;i<nmodes->imodes*nmodes->nspec;i++){
#ifdef STEP
          nmodes->decmodrr[i][k] += nmodes->permod[i][j]*nmodes->aniso_raman[j];
          nmodes->decmod_isorr[i][k] += nmodes->permod[i][j]*nmodes->iso_raman[j];
#else
          nmodes->decmodrr[i][k]   += (1.-fr)*nmodes->permod[i][j]*nmodes->aniso_raman[j];
          nmodes->decmodrr[i][k+1] +=    (fr)*nmodes->permod[i][j]*nmodes->aniso_raman[j];
          nmodes->decmod_isorr[i][k]   += (1.-fr)*nmodes->permod[i][j]*nmodes->iso_raman[j];
          nmodes->decmod_isorr[i][k+1] +=    (fr)*nmodes->permod[i][j]*nmodes->iso_raman[j];
#endif
        }
      }
    }
  }
}
/*------------------------------------------------------------------*/


