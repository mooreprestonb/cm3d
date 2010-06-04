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

/* routine have to do with pressure and volume */

#include "md.h"
/* #define ISO */
/* #define ZERO_TM */
/*------------------------------------------------------------------*/
void set_baro(SIMPARMS *simparms,COORDS *coords,double tauv)
{
  int i;
  
  if(simparms->nbar<0){
    md_error("you can't have a negative number for the volume chains!");
  }
  if(simparms->nbar>0) {
    coords->mbar = (double *)cmalloc(simparms->nbar*sizeof(double));
  }

  coords->mvol = (simparms->ndof)*simparms->temp*tauv*tauv/3.;
#ifdef ISO
  md_stdout("ISOtropic cell fluctionation is Defined\n");
#endif

  for(i=0;i<simparms->nbar;i++) coords->mbar[i] = simparms->temp*tauv*tauv;
}
/*----------------------------------------------------------------*/
void forceV(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies)
{
  int i;
  double pv,vol,p2m,dnf,am,p2m_tensor[9];
  double *vx,*vy,*vz,*amass,*fvol9;

  vx = coords->vx;
  vy = coords->vy;
  vz = coords->vz;
  amass = coords->amass;
  vol = get_deth(coords->hmat);
  fvol9 = coords->fvol9;

  p2m = 0.;
  for(i=0;i<9;i++) p2m_tensor[i] = 0.0;
  dnf = 1./((double)simparms->ndof);
#include "vectorize.h"
  for(i=0;i<simparms->natoms;i++){
    am = amass[i];
    p2m += am*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);

    p2m_tensor[0] += am*(vx[i]*vx[i]);
    p2m_tensor[1] += am*(vy[i]*vx[i]);
    p2m_tensor[2] += am*(vz[i]*vx[i]);
    p2m_tensor[3] += am*(vx[i]*vy[i]);
    p2m_tensor[4] += am*(vy[i]*vy[i]);
    p2m_tensor[5] += am*(vz[i]*vy[i]);
    p2m_tensor[6] += am*(vx[i]*vz[i]);
    p2m_tensor[7] += am*(vy[i]*vz[i]);
    p2m_tensor[8] += am*(vz[i]*vz[i]);
  }
    
  pv = (1.+3.*dnf)*p2m+energies->W+energies->WI;
  coords->fvol = (pv-DIM*simparms->pext*vol)/(3.*coords->mvol);

  /* construct matrix */
  for(i=0;i<9;i++){
    fvol9[i]=p2m_tensor[i]+energies->Wtensor[i]+energies->WItensor[i];
  }

  /* add diagonal part */
  p2m = p2m*dnf-simparms->pext*vol;
  fvol9[0] += p2m;
  fvol9[4] += p2m; 
  fvol9[8] += p2m; 

  /* mass weight */
  pv = 1./coords->mvol;
  for(i=0;i<9;i++) fvol9[i] *= pv;
  sym33(fvol9); 

  switch(simparms->ivol){
  case 0:  /* no volume fluctuations (do nothing)*/
  case 3:  /* fully flexable (do nothing) */
    break;
  case 2:  /* only lenghts (zero off diagonal components)*/
    fvol9[1] = fvol9[3] = 0.;
    fvol9[2] = fvol9[6] = 0.;
    fvol9[5] = fvol9[7] = 0.;
    break;
  case 1: /* isotropic (zero off diagonal components and average trace*/
    fvol9[1] = fvol9[3] = 0.;
    fvol9[2] = fvol9[6] = 0.;
    fvol9[5] = fvol9[7] = 0.;
    pv = (fvol9[0] + fvol9[4] + fvol9[8])/3.;
    fvol9[0] = fvol9[4] = fvol9[8] = pv;
    break;
  }
}
/*----------------------------------------------------------------*/
void samvbar(SIMPARMS *simparms,COORDS *coords)
{
  int i;
  double fact,tm;
  
  /* the volume velocity */
  ggauss(1,&coords->vvol);
  ggauss(9,coords->vvol9);
  coords->vvol9[1] = coords->vvol9[3];
  coords->vvol9[2] = coords->vvol9[6];
  coords->vvol9[5] = coords->vvol9[7];

  switch(simparms->ivol){
  case 0:  /* no volume fluctuations (do nothing)*/
  case 3:  /* fully flexable (do nothing) */
    break;
  case 2:  /* only lenghts (zero off diagonal components) */
    coords->vvol9[1] = coords->vvol9[3] = 0.;
    coords->vvol9[2] = coords->vvol9[6] = 0.;
    coords->vvol9[5] = coords->vvol9[7] = 0.;
    break;
  case 1: /* isotropic (zero off diagonal components and average trace */
    coords->vvol9[1] = coords->vvol9[3] = 0.;
    coords->vvol9[2] = coords->vvol9[6] = 0.;
    coords->vvol9[5] = coords->vvol9[7] = 0.;
    coords->vvol9[0] = coords->vvol9[4] = coords->vvol9[8] = coords->vvol;
    break;
  }

  /* the +9 is for the volume (should be +1 for isotropic,
     but it doesn't hurt anything */


#ifdef ZERO_TM /* zero the momentum of the barostat */
#ifdef ISO
  fact = 3.*coords->mvol*coords->vvol;
  tm = 3.*coords->mvol;
#else
  for(i=0,fact = 0;i<9;i++) fact += coords->vvol9[i];
  fact *= coords->mvol;
  tm = 9.*coords->mvol;
#endif
  fact /= tm;
  coords->vvol -= fact;
  for(i=0;i<9;i++) coords->vvol9[i] -= fact;

  fact = 3.*coords->mvol*coords->vvol*coords->vvol;
  fact = sqrt(((double)(simparms->nbar+9))*simparms->temp/fact);
  coords->vvol *= fact;

  for(i=0,fact=0.;i<9;i++) fact += coords->vvol9[i]*coords->vvol9[i];
  fact *= coords->mvol;
  fact = sqrt(9.*simparms->temp/fact);
  for(i=0;i<9;i++) coords->vvol9[i] *= fact;
#else 
  coords->vvol = 0.;
  for(i=0;i<9;i++) coords->vvol9[i] = 0.;
#endif

  for(i=0;i<simparms->nbar;i++) coords->bar[i] = 0.; /* initial pos baro */

  if(simparms->nbar>0){
    
    ggauss(simparms->nbar,coords->vbar);
    if(simparms->nbar>1){
      fact = tm = 0.;
      for(i=0;i<simparms->nbar;i++) {
	fact += coords->mbar[i]*coords->vbar[i];
	tm += coords->mbar[i];
      }
      fact /= tm;
      for(i=0;i<simparms->nbar;i++) coords->vbar[i] -= fact;
    }

    for(i=0,fact=0.;i<simparms->nbar;i++) {
      tm = coords->vbar[i];
      fact += coords->mbar[i]*tm*tm;
    }
    fact = sqrt(((double)(simparms->nbar))*simparms->temp/fact);
    for(i=0;i<simparms->nbar;i++) coords->vbar[i] *= fact;
  }
}

/*----------------------------------------------------------------*/
void getengv(SIMPARMS *simparms,COORDS *coords,double *zkev,double *potv)
{
  int i;
  double svel2;

  if(simparms->ivol){
    svel2 = simparms->pext*get_deth(coords->hmat);
    for(i=0;i<simparms->nbar;i++) svel2 += simparms->temp*coords->bar[i];
    *potv = svel2;

#ifdef ISO
    svel2 = 3.*coords->mvol*coords->vvol*coords->vvol;
#else
    for(svel2=0.,i=0;i<9;i++) svel2 += coords->vvol9[i]*coords->vvol9[i];
    svel2 *= coords->mvol;
#endif
    for(i=0;i<simparms->nbar;i++){
      svel2 += coords->mbar[i]*coords->vbar[i]*coords->vbar[i];
    }
    *zkev = .5*svel2;
   } else {
    *potv = *zkev = 0.;
   }
}

/*--------------------------------------------------------------------*/
void getengv_cheese(SIMPARMS *simparms,COORDS *coords,double *zkev,
		    double *potv,double *zken)
{
  int i;
  double svel2;

  if (simparms->ivol){
    svel2 = simparms->pext*get_deth(coords->hmat);
    for(i=0;i<simparms->nbar;i++) svel2 += simparms->temp*coords->bar[i];
    *potv = svel2;

#ifdef ISO
    svel2 = 3.*coords->mvol*coords->vvol*coords->vvol;
#else
    for(svel2=0,i=0;i<9;i++) svel2 += coords->vvol9[i]*coords->vvol9[i];
    svel2 *= coords->mvol;
#endif
    *zkev = .5*svel2;
   } else {
    *potv = *zkev = 0.;
   }

  for(svel2=0.,i=0;i<simparms->nbar;i++){
    svel2 += coords->mbar[i]*coords->vbar[i]*coords->vbar[i];
  }
  *zken += .5*svel2;
}

/*--------------------------------------------------------------------*/
void intbar(SIMPARMS *simparms,COORDS *coords,double dta2,int iflag)
{
  int i,k,iyosh;
  double expvV,dnf;
  double dtanc2;

  dnf = 1./(double)simparms->ndof;
  
  switch(simparms->ivol){
  case 0:
  case 3:  /* fully flexable */
    break;
  case 2:  /* only lenghts */
    coords->vvol9[1] = coords->vvol9[3] = 0.;
    coords->vvol9[2] = coords->vvol9[6] = 0.;
    coords->vvol9[5] = coords->vvol9[7] = 0.;
    break;
  case 1: /* isotropic */
    coords->vvol9[1] = coords->vvol9[3] = 0.;
    coords->vvol9[2] = coords->vvol9[6] = 0.;
    coords->vvol9[5] = coords->vvol9[7] = 0.;
    expvV = (coords->vvol9[0] + coords->vvol9[4] + coords->vvol9[8])/3.;
    coords->vvol9[0] = coords->vvol9[4] = coords->vvol9[8] = expvV;
    break;
  }

  sym33(coords->vvol9);

  switch(iflag){
  case 0:
/* do the nhc and the yoshida loop */
   for (k=0;k<simparms->nhc;k++) {
   for (iyosh=0;iyosh<simparms->nyoshida;iyosh++) {
       dtanc2=simparms->wyosh[iyosh]*dta2/((float) simparms->nhc);
    app_ux2(simparms->natoms,dtanc2,coords);
    app_up2(simparms->natoms,dtanc2,dnf,coords);
    app_uh(dtanc2,coords);

    /* Volume Chain */
    if(simparms->nbar>0){
      expvV = exp(-dtanc2*coords->vbar[0]);
#ifdef ISO
      coords->vvol *= expvV;
#endif
      for(i=0;i<9;i++) coords->vvol9[i] *= expvV;
      for(i=0;i<simparms->nbar-1;i++){
	coords->vbar[i] *= exp(-dtanc2*coords->vbar[i+1]);
      }
      for(i=simparms->nbar-1;i>=1;i--) {
	coords->vbar[i] += dtanc2*(coords->mbar[i-1]*coords->vbar[i-1]*
				 coords->vbar[i-1]-simparms->temp)/
				   coords->mbar[i];
      }
#ifdef ISO
      coords->vbar[0] += (dtanc2*(3.*coords->mvol*coords->vvol*coords->vvol-
				simparms->temp)/coords->mbar[0]);
#else
      expvV = 0.;
      for(i=0;i<9;i++) expvV += coords->vvol9[i]*coords->vvol9[i];
      coords->vbar[0] += (dtanc2*(coords->mvol*expvV-simparms->temp)/
			  coords->mbar[0]);
#endif
      for(i=0;i<simparms->nbar;i++)coords->bar[i] += dtanc2*coords->vbar[i];
    }
   }}
    break;
  case 1:
   for (k=0;k<simparms->nhc;k++) {
   for (iyosh=0;iyosh<simparms->nyoshida;iyosh++) {
       dtanc2=simparms->wyosh[iyosh]*dta2/((float) simparms->nhc);
    if(simparms->nbar>0){
      for(i=0;i<simparms->nbar;i++)coords->bar[i]  += dtanc2*coords->vbar[i];
#ifdef ISO
      coords->vbar[0] += (dtanc2*(3.*coords->mvol*coords->vvol*coords->vvol-
				simparms->temp)/coords->mbar[0]);
#else
      expvV=0.;
      for(i=0;i<9;i++) expvV += coords->vvol9[i]*coords->vvol9[i];
      coords->vbar[0] += (dtanc2*(coords->mvol*expvV-simparms->temp)/
			  coords->mbar[0]);
#endif
      for(i=1;i<simparms->nbar;i++){
	coords->vbar[i] += dtanc2*(coords->mbar[i-1]*coords->vbar[i-1]*
				 coords->vbar[i-1]-simparms->temp)/
				   coords->mbar[i];
      }
      for(i=simparms->nbar-2;i>=0;i--){
	coords->vbar[i] *= exp(-dtanc2*coords->vbar[i+1]);
      }
      expvV = exp(-dtanc2*coords->vbar[0]);
      for(i=0;i<9;i++) coords->vvol9[i] *= expvV;
#ifdef ISO
      coords->vvol *= expvV;
#endif
    }
    sym33(coords->vvol9);
    app_uh(dtanc2,coords);
    app_up2(simparms->natoms,dtanc2,dnf,coords);
    app_ux2(simparms->natoms,dtanc2,coords);

    }
    }
    break;
  default:
    md_error("Call technical support trouble with int_baro!!!");
    break;
  }

  /* isotropic cell fluctuations */
#ifdef ISO
  expvV = exp(coords->pvol)/get_deth(coords->hmat);
  for(i=0;i<9;i++){coords->hmat[i] *= expvV;}
#endif
  gethinv9(coords->hmat,coords->hmati);
}
/*--------------------------------------------------------------------*/
void app_ux2(int natoms,double dta2,COORDS *coords)
{
  int i;
  double evol[3],**cvol,*vvol9;
  double *px,*py,*pz,rx,ry,rz,ev0,ev1,ev2;

  vvol9 = coords->vvol9;
  px = coords->px;
  py = coords->py;
  pz = coords->pz;
  
#ifdef ISO
  ev0 = exp(dta2*coords->vvol);
  for(i=0;i<natoms*3;i++) px[i] *= ev0;
#else 
  cvol = dmatrix(0,2,0,2);
  cvol[0][0] = vvol9[0];  cvol[0][1] = vvol9[1];  cvol[0][2] = vvol9[2];
  cvol[1][0] = vvol9[3];  cvol[1][1] = vvol9[4];  cvol[1][2] = vvol9[5];
  cvol[2][0] = vvol9[6];  cvol[2][1] = vvol9[7];  cvol[2][2] = vvol9[8];
  rs_me(3,evol,cvol,1);

  ev0 = exp(dta2*evol[0]);
  ev1 = exp(dta2*evol[1]);
  ev2 = exp(dta2*evol[2]);

  for(i=0;i<natoms;i++){
    rx = px[i]*cvol[0][0] + py[i]*cvol[1][0] + pz[i]*cvol[2][0];
    ry = px[i]*cvol[0][1] + py[i]*cvol[1][1] + pz[i]*cvol[2][1];
    rz = px[i]*cvol[0][2] + py[i]*cvol[1][2] + pz[i]*cvol[2][2];

    rx *= ev0;
    ry *= ev1;
    rz *= ev2;

    px[i] = rx*cvol[0][0] + ry*cvol[0][1] + rz*cvol[0][2];
    py[i] = rx*cvol[1][0] + ry*cvol[1][1] + rz*cvol[1][2];
    pz[i] = rx*cvol[2][0] + ry*cvol[2][1] + rz*cvol[2][2];
  }
  free_dmatrix(cvol,0,2,0,2);
#endif
}

/*----------------------------------------------------------------*/
void app_up2(int natoms,double dta2,double dnf,COORDS *coords)
{
  int i;
  double evol[3],**cvol,tr,*vvol9;
  double *vx,*vy,*vz,rx,ry,rz,ev0,ev1,ev2;

  vvol9 = coords->vvol9;
  vx = coords->vx;
  vy = coords->vy;
  vz = coords->vz;

#ifdef ISO
  ev0 = exp(-dta2*(1.+3.*dnf)*coords->vvol);
  for(i=0;i<natoms*3;i++) vx[i] *= ev0;
#else
  cvol = dmatrix(0,2,0,2);
  tr = (vvol9[0] + vvol9[4] + vvol9[8])*dnf; 
  cvol[0][0] = vvol9[0]+tr; cvol[0][1] = vvol9[1];    cvol[0][2] = vvol9[2];
  cvol[1][0] = vvol9[3];    cvol[1][1] = vvol9[4]+tr; cvol[1][2] = vvol9[5];
  cvol[2][0] = vvol9[6];    cvol[2][1] = vvol9[7];    cvol[2][2] = vvol9[8]+tr;
  rs_me(3,evol,cvol,1);

  ev0 = exp(-dta2*evol[0]);
  ev1 = exp(-dta2*evol[1]);
  ev2 = exp(-dta2*evol[2]);

  for(i=0;i<natoms;i++){
    rx = vx[i]*cvol[0][0] + vy[i]*cvol[1][0] + vz[i]*cvol[2][0];
    ry = vx[i]*cvol[0][1] + vy[i]*cvol[1][1] + vz[i]*cvol[2][1];
    rz = vx[i]*cvol[0][2] + vy[i]*cvol[1][2] + vz[i]*cvol[2][2];
    
    rx *= ev0;
    ry *= ev1;
    rz *= ev2;
    
    vx[i] = rx*cvol[0][0] + ry*cvol[0][1] + rz*cvol[0][2];
    vy[i] = rx*cvol[1][0] + ry*cvol[1][1] + rz*cvol[1][2];
    vz[i] = rx*cvol[2][0] + ry*cvol[2][1] + rz*cvol[2][2];
  }
  free_dmatrix(cvol,0,2,0,2);
#endif
}

/*----------------------------------------------------------------*/
void app_uh(double dta2,COORDS *coords)
{
  int i;
  double evol[3],**cvol,*vvol9;
  double rx,ry,rz,*hm,ev0,ev1,ev2;

#ifdef ISO
  coords->pvol = log(get_deth(coords->hmat)) + dta2*coords->vvol;
#else
  vvol9 = coords->vvol9;
  cvol = dmatrix(0,2,0,2);
  hm = coords->hmat;
  cvol[0][0] = vvol9[0];  cvol[0][1] = vvol9[1];  cvol[0][2] = vvol9[2];
  cvol[1][0] = vvol9[3];  cvol[1][1] = vvol9[4];  cvol[1][2] = vvol9[5];
  cvol[2][0] = vvol9[6];  cvol[2][1] = vvol9[7];  cvol[2][2] = vvol9[8];
  rs_me(3,evol,cvol,1);

  ev0 = exp(dta2*evol[0]);
  ev1 = exp(dta2*evol[1]);
  ev2 = exp(dta2*evol[2]);

  for(i=0;i<3;i++){
    rx = hm[i]*cvol[0][0] + hm[i+3]*cvol[1][0] + hm[i+6]*cvol[2][0];
    ry = hm[i]*cvol[0][1] + hm[i+3]*cvol[1][1] + hm[i+6]*cvol[2][1];
    rz = hm[i]*cvol[0][2] + hm[i+3]*cvol[1][2] + hm[i+6]*cvol[2][2];
    rx *= ev0;
    ry *= ev1;
    rz *= ev2;
    hm[i  ] = rx*cvol[0][0] + ry*cvol[0][1] + rz*cvol[0][2];
    hm[i+3] = rx*cvol[1][0] + ry*cvol[1][1] + rz*cvol[1][2];
    hm[i+6] = rx*cvol[2][0] + ry*cvol[2][1] + rz*cvol[2][2];
  }
  free_dmatrix(cvol,0,2,0,2);
#endif
}
