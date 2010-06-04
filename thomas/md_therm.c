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

/* subroutines dealing with nose thermostats */

#include "md.h"


/*------------------------------------------------------------------*/
void set_therm(SIMPARMS *simparms,COORDS *coords,
	       int nspec,SPECIES spec_root,double tau_nhc)
{
  char line[MAXLINELEN];
  int i,j,k,ispec,ioff,*itherm;
  SPECIES *spec_now;

   /* set thermostats */

  simparms->ntherm = 1;
  coords->ithm = (int *)cmalloc(simparms->natoms*sizeof(int));

  /* temporary array */
  itherm = (int *)cmalloc(simparms->ntherm*sizeof(int));
  itherm[simparms->ntherm-1] = 0;

  ioff = 0;
  /* loop over the species and determine thermostats */
  for(ispec=0,spec_now=&spec_root;ispec<nspec;
      ispec++,spec_now=spec_now->next){
    
    if(strcasecmp(spec_now->thermopt,"none")==0){
      /* no thermostat for these atoms */
      for(k=0;k<spec_now->nmol;k++){
        for(j=0;j<spec_now->napm;j++){
          coords->ithm[ioff++] = -1;
        }
      }
    } else if (strcasecmp(spec_now->thermopt,"global")==0){
      /* attach these to the global thermostat */
      for(k=0;k<spec_now->nmol;k++){
        for(j=0;j<spec_now->napm;j++){
          coords->ithm[ioff++] = 0;
          itherm[0]++;
        }
      }
    } else if(strcasecmp(spec_now->thermopt,"glob_mol")==0){
      /* attach all molecule of this type to a thermostat */
      ++simparms->ntherm;
      itherm = (int *)realloc(itherm,simparms->ntherm*sizeof(int));
      itherm[simparms->ntherm-1] = 0;
      for(k=0;k<spec_now->nmol;k++){
        for(j=0;j<spec_now->napm;j++){
          coords->ithm[ioff++] = simparms->ntherm-1;
          itherm[simparms->ntherm-1]++;
        }
      }
    } else if(strcasecmp(spec_now->thermopt,"ind_mol")==0){
     /* attach each molecule to a thermostat */
      for(k=0;k<spec_now->nmol;k++){
        ++simparms->ntherm;
        itherm = (int *)realloc(itherm,simparms->ntherm*sizeof(int));
        itherm[simparms->ntherm-1] = 0;
        for(j=0;j<spec_now->napm;j++){
          coords->ithm[ioff++] = simparms->ntherm-1;
          itherm[simparms->ntherm-1]++;
        }
      }
    } else if(strcasecmp(spec_now->thermopt,"atm_mol")==0){
      /* attach each atom to a thermostat */
      for(k=0;k<spec_now->nmol;k++){
        for(j=0;j<spec_now->napm;j++){
          ++simparms->ntherm;
          itherm = (int *)realloc(itherm,simparms->ntherm*sizeof(int));
          itherm[simparms->ntherm-1] = 0;
          coords->ithm[ioff++] = simparms->ntherm-1;
          itherm[simparms->ntherm-1]++;
        }
      }
    } else {
      sprintf(line,"thermostat option must be none,global,glob_mol,ind_mol");
      sprintf(line,"%s or atm_mol\n\tnot \"%s\" (species = %d)",
	      line,spec_now->thermopt,ispec);
      md_error(line);
    }
  }
  if(ioff != simparms->natoms){
    sprintf(line,"while setting up thermostats! (set) %d != %d atoms",
	    ioff,simparms->natoms);
    md_error(line);
  }
  /* get rid of zero itherms */
  for(i=0;i<simparms->ntherm;i++){
    if(itherm[i]==0){
      for(j=0;j<simparms->natoms;j++){
        if(coords->ithm[j]>=i) coords->ithm[j]--;
      }
      --simparms->ntherm;
      for(j=i;j<simparms->ntherm;j++){
        itherm[i] = itherm[i+1];
      }
    }
  }

#ifdef DEBUG
  /*
  for(i=0;i<ioff;i++){
    printf("ithm[%d] = %d\n",i,coords->ithm[i]);
  } */
  for(i=0;i<simparms->ntherm;i++){
    printf("itherm[%d] = %d\n",i,itherm[i]);
  }
#endif

  /* check for consistancy */
  if(simparms->ntherm != 0 && 
     (simparms->iensemble == 0 || simparms->iensemble==3)){
    sprintf(line,"you have specified %d thermostats in the setfile\n and the NVE ensemble in the input file (fix one)!",simparms->ntherm);
    md_error(line);
  }
  /* check for consistancy */
  if(simparms->ntherm == 0 && 
     (simparms->iensemble == 1 || simparms->iensemble==2)){
    sprintf(line,"you have specified %d thermostats in the setfile\n and NOT the NVE ensemble in the input file (fix one)!",simparms->ntherm);
    md_error(line);
  }
  
  if(simparms->ntherm>0){
    /* allocate memory for static arrays to be used when thermostats are on*/
    coords->feta = dmatrix(0,simparms->ntherm-1,0,simparms->nchain-1);
    coords->gkt  = dmatrix(0,simparms->ntherm-1,0,simparms->nchain-1);
    coords->mnh  = dmatrix(0,simparms->ntherm-1,0,simparms->nchain-1);
    coords->p2mt = (double *)cmalloc(simparms->ntherm*sizeof(double));
  }
  for(j=0;j<simparms->ntherm;j++){
    coords->gkt[j][0] = simparms->temp*DIM*itherm[j];
    for(i=1;i<simparms->nchain;i++){
      coords->gkt[j][i] += simparms->temp;
    }
  }
  for(j=0;j<simparms->ntherm;j++){
    for(i=0;i<simparms->nchain;i++) 
      coords->mnh[j][i] = coords->gkt[j][i]*tau_nhc*tau_nhc;
  }

/* initialize the yoshidas */
if(strcasecmp(simparms->yoshidas,"on")==0){
   simparms->nyoshida=3;
   simparms->wyosh[0]= simparms->wyosh[2]=1.0/(2.0-pow(2.0,1./3.));
   simparms->wyosh[1]=1.0-2.0* simparms->wyosh[0];
} else if(strcasecmp(simparms->yoshidas,"off")==0){
   simparms->nyoshida=1;
   simparms->wyosh[0]= simparms->wyosh[2]=1.0;
   simparms->wyosh[1]=1.0;
} else {
      sprintf(line,"yoshida options on/off and not %s",simparms->yoshidas);
      md_error(line);
}

   

      /* no thermostat for these atoms */


  free(itherm);
}

/*------------------------------------------------------------------*/

void forcen(SIMPARMS *simparms,COORDS *coords)
{
  int i,j,*ithm;
  double *am,*vx,*vy,*vz,*p2m;
  
  ithm = coords->ithm;
  p2m = coords->p2mt;
  am  = coords->amass;
  vx  = coords->vx;
  vy  = coords->vy;
  vz  = coords->vz;

  for(i=0;i<simparms->ntherm;i++) p2m[i] = 0.;

  for(i=0;i<simparms->natoms;i++){
    p2m[ithm[i]] += am[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
  }

  for(j=0;j<simparms->ntherm;j++){
    coords->feta[j][0] = (p2m[j] - coords->gkt[j][0])/coords->mnh[j][0];

/*    for(i=1;i<simparms->nchain;i++){
      coords->feta[j][i]=(coords->mnh[j][i-1]*coords->veta[j][i-1]*
			  coords->veta[j][i-1]-coords->gkt[j][i])/
			    coords->mnh[j][i];
    }
*/
  }
}
/*----------------------------------------------------------------*/
void getengn(SIMPARMS *simparms,COORDS *coords,double *zken,double *potn)
{
  int i,j;
  double svel2,sqr;
  
  *zken = *potn = 0.;
  
  for(j=0;j<simparms->ntherm;j++){
    svel2 = 0.;
    for(i=0;i<simparms->nchain;i++) {
      svel2 += coords->gkt[j][i]*coords->eta[j][i];
    }
    *potn += svel2;
    
    svel2 = 0.;
    for(i=0;i<simparms->nchain;i++){
      sqr = coords->veta[j][i];
      svel2 += coords->mnh[j][i]*sqr*sqr;
    }
    *zken += .5*svel2;
  }
}

/*----------------------------------------------------------------*/
void getengn_cheese(SIMPARMS *simparms,COORDS *coords,
		    double *zken,double *potn)
{
  int i,j;
  double svel2,sqr;

  *zken = *potn = 0.;

  for(j=0;j<simparms->ntherm;j++){
    svel2 = 0.;
    for(i=0;i<simparms->nchain;i++){
      svel2 += coords->gkt[j][i]*coords->eta[j][i];
    }
    
    *potn += svel2;
    
    svel2 = 0.;
    for(i=0;i<simparms->nchain;i++) {
      sqr = coords->veta[j][i];
      svel2 += coords->mnh[j][i]*sqr*sqr;
    }
    *zken += .5*svel2;
  }
}

/*----------------------------------------------------------------*/
void samveta(SIMPARMS *simparms,COORDS *coords)
{
  int i,j;
  double fact,tm;

  for(j=0;j<simparms->ntherm;j++){
    for(i=0;i<simparms->nchain;i++) (coords->eta)[j][i]=0.0;

    ggauss(simparms->nchain,(coords->veta)[j]);
    
    if(simparms->nchain >1){
      fact = tm = 0.;

      for(i=0;i<simparms->nchain;i++) {
	fact += coords->mnh[j][i]*coords->veta[j][i];
	tm += coords->mnh[j][i];
      }
      
      fact /= tm;
      for(i=0;i<simparms->nchain;i++) coords->veta[j][i] -= fact;
    }
    fact = 0.;
    for(i=0;i<simparms->nchain;i++) {
      tm = coords->veta[j][i];
      fact += coords->mnh[j][i]*tm*tm;
    }
    
    if(simparms->nchain>0) fact = sqrt(((double)simparms->nchain)*
				    simparms->temp/(fact));
    for(i=0;i<simparms->nchain;i++) coords->veta[j][i] *= fact;
  }
}

/*---------------------------------------------------------------------*/
void intnhc(SIMPARMS *simparms,COORDS *coords,double dta2)
{
  int i,j,*ithm,nchain,k,iyosh;
  double expveta,dta,dta4,scale,forceeta,*scalevel,dtanc,dtanc2,dtanc4;
  double *vx,*vy,*vz,**eta,**veta,**feta,*p2m,*am;
  double wyosh[3],nc;

  
  scalevel = (double *)cmalloc(simparms->ntherm*sizeof(double)); 
  nchain=simparms->nchain;
  dta = 2.*dta2;
  dta4 = 0.5 * dta2;
  vx = coords->vx;
  vy = coords->vy;
  vz = coords->vz;
  eta = coords->eta;
  veta = coords->veta;
  feta = coords->feta;
  ithm = coords->ithm;
  p2m = coords->p2mt;
  am  = coords->amass;
 
  for(j=0;j<simparms->ntherm;j++) scalevel[j] = 1.0;
  veta[-1][0] = 0.;
  /* get funky velocities */

  forcen(simparms,coords);

  nc=(float) simparms->nhc;
  dtanc=dta/nc;
  dtanc2=dta2/nc;
  dtanc4=dta4/nc;

for(k=0;k<simparms->nhc;k++)  
{
for(iyosh=0;iyosh<simparms->nyoshida;iyosh++)
{
  dtanc =simparms->wyosh[iyosh]*dta/nc;
  dtanc2=simparms->wyosh[iyosh]*dta2/nc;
  dtanc4=simparms->wyosh[iyosh]*dta4/nc;

  for(j=0;j<simparms->ntherm;j++)
  { 
/* now do the outermost one */
    forceeta=(coords->mnh[j][simparms->nchain-2]*veta[j][simparms->nchain-2]*veta[j][simparms->nchain-2] 
             - coords->gkt[j][simparms->nchain-1])/coords->mnh[j][simparms->nchain-1];
    veta[j][simparms->nchain-1] += dtanc2*forceeta;

    for(i=simparms->nchain-2;i>=1;i--)
    { scale = exp(-dtanc4*veta[j][i+1]);
      forceeta=(coords->mnh[j][i-1]*veta[j][i-1]*veta[j][i-1] 
               - coords->gkt[j][i])/coords->mnh[j][i];
      veta[j][i] = veta[j][i] * scale * scale + dtanc2 * forceeta * scale; 
    }
/* take care of the innermost chain last */
      scale = exp(-dtanc4*veta[j][1]);

      forceeta = (p2m[j] -coords->gkt[j][0])/coords->mnh[j][0];
      veta[j][0] = veta[j][0] * scale * scale + dtanc2 * forceeta * scale; 

  }

  for(j=0;j<simparms->ntherm;j++)  
    for(i=0;i<simparms->nchain;i++)  eta[j][i]  += dtanc*veta[j][i];

  for(j=0;j<simparms->ntherm;j++)  
  {
      expveta = exp(-dtanc*veta[j][0]);
      p2m[j] = p2m[j]*expveta*expveta;
      scalevel[j] *= expveta;
  }


  for(j=0;j<simparms->ntherm;j++)
  { 
/* take care of the innermost chain first */
      scale = exp(-dtanc4*veta[j][1]);
      forceeta = (p2m[j] -coords->gkt[j][0])/coords->mnh[j][0];
      veta[j][0] = veta[j][0] * scale * scale + dtanc2 * forceeta * scale; 

    for(i=1;i<simparms->nchain-1;i++)
    { scale = exp(-dtanc4*veta[j][i+1]);
      forceeta=(coords->mnh[j][i-1]*veta[j][i-1]*veta[j][i-1] 
               - coords->gkt[j][i])/coords->mnh[j][i];
      veta[j][i] = veta[j][i] * scale * scale + dtanc2 * forceeta * scale; 
    }
/* now do the outermost one */
    forceeta=(coords->mnh[j][simparms->nchain-2]*veta[j][simparms->nchain-2]*veta[j][simparms->nchain-2] 
             - coords->gkt[j][simparms->nchain-1])/coords->mnh[j][simparms->nchain-1];
    veta[j][simparms->nchain-1] += dtanc2*forceeta;
  }
}
}
/* now scale the velocities */
    for(i=0;i<simparms->natoms;i++){
      vx[i] *= scalevel[ithm[i]];
      vy[i] *= scalevel[ithm[i]];
      vz[i] *= scalevel[ithm[i]];
    }
free(scalevel);
}
/*-----------------------------------------------------------*/

