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

/* routines that zero total momentum, scale temperature, zero COM, etc.. */

#include "md.h"

/* #define DEBUG */
/*----------------------------------------------------------------------*/
void samvel(SIMPARMS *simparms,COORDS *coords,int irandom)
{
  /* take advantange of the fact that the velocities are contiguous */
  if(irandom){
    md_warning("Resampling velocities from a gaussian distribution");
    ggauss(simparms->natoms*3,coords->vx);
  }
  zerotm(simparms->natoms,coords->vx,coords->vy,coords->vz,coords->amass);  
  if(simparms->nfreeze>0) freeze_atom(simparms,coords);
  
  simparms->temp *= simparms->scaletemp;
  scale(simparms->natoms,coords->vx,coords->vy,coords->vz,coords->amass,
	simparms->temp,simparms->ndof);

}
/*--------------------------------------------------------------------*/
void zerotm(int natoms,double *vx,double *vy,double *vz,double *amass)
{
  int i;
  double xvel,yvel,zvel,mass,tmass;

  tmass = xvel = yvel = zvel = 0.0;

  if(natoms > 1){
    for(i=0;i<natoms;i++){
      mass = amass[i];    tmass += mass;
      xvel += vx[i]*mass; yvel += vy[i]*mass; zvel += vz[i]*mass;
    }
    xvel /= tmass; yvel /= tmass; zvel /= tmass;
    
    for(i=0;i<natoms;i++){
      vx[i] -= xvel; vy[i] -= yvel; vz[i] -= zvel;
    }
  }
}
/*---------------------------------------------------------------------*/
double get_temp(int natoms,double *vx,double *vy,double *vz,
		double *amass,int ndof)
{
  int i;
  double p2m;
  
  p2m = 0.;
  for(i=0;i<natoms;i++){
    p2m += amass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
  }
  return (p2m/(double)ndof);
}
/*---------------------------------------------------------------------*/
void scale(int natoms,double *vx,double *vy,double *vz,
	   double *amass,double temp,int ndof)
{
  int i;
  LINE line;
  double fact,svel2;
  
  svel2 = 0.0;
  for(i=0;i<natoms;i++)svel2 +=amass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);

  if(svel2<=0.0){
    sprintf(line,"You have all zero velocities!");
  } else {
    sprintf(line,"Scaled temperature from %g K to %g K",
	    svel2/((double)ndof),temp);
    
    fact = sqrt((double)ndof*temp/svel2);
    for(i=0;i<natoms;i++){
      vx[i] *= fact; vy[i] *= fact; vz[i] *= fact;
    }
  }
  md_warning(line);
  if(svel2 != svel2){
    md_error("A NaN was incountered while rescaling velocties!");
  }
}

/*-----------------------------------------------------------*/
void radial(int natoms,double *px,double *py,double *pz,
	    double *vx,double *vy,double *vz,double *amass)
{
  int i;
  double r,v;

  zerotm(natoms,vx,vy,vz,amass);
  zerocm(natoms,px,py,pz,amass);
  for(i=0;i<natoms;i++){
    v = sqrt(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    r = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);

    vx[i] = v*px[i]/r;
    vy[i] = v*py[i]/r;
    vz[i] = v*pz[i]/r;
  }
  zerotm(natoms,vx,vy,vz,amass);
}
/*-----------------------------------------------------------*/
void zerocm(int natoms,double *px,double *py,double *pz,double *amass)
{
  int j;
  double cmx,cmy,cmz,tm;
  cmx = cmy = cmz = tm = 0.;
  for(j=0;j<natoms;j++){
    cmx += px[j]*amass[j]; cmy += py[j]*amass[j];  cmz += pz[j]*amass[j];
    tm += amass[j];
  }
  cmx /= tm; cmy /= tm; cmz /= tm;
  for(j=0;j<natoms;j++){
    px[j] -= cmx; py[j] -= cmy; pz[j] -= cmz;
  }
}
/*--------------------------------------------------------------*/
void zero_angm(int natoms,double *px,double *py,double *pz,
	       double *vx,double *vy,double *vz,double *amass)
{
  int i;
  double dx,dy,dz,tmass;
  double xvel,yvel,zvel,amt;
  double **imat,**invmat,l[3],omega[3];
#ifdef DEBUG
  double dvx,dvy,dvz,cmx,cmy,cmz;
#endif

  md_warning("Zeroing total angular momentum");
 
  zerocm(natoms,px,py,pz,amass);

  if(natoms==1) return;

  imat   = dmatrix(0,2,0,2);
  invmat = dmatrix(0,2,0,2);

  for(i=0,tmass=0.;i<natoms;i++) tmass += amass[i];

  if(natoms==2){
    radial(natoms,px,py,pz,vx,vy,vz,amass);
  } else {
    
    /* zero total linear momentum */
    xvel = yvel = zvel = 0.;
    for(i=0;i<natoms;i++){
      amt = amass[i];
      xvel += vx[i]*amt; yvel += vy[i]*amt; zvel += vz[i]*amt;
    }
    xvel /= tmass; yvel /= tmass; zvel /= tmass;
    for(i=0;i<natoms;i++){
      vx[i] -= xvel; vy[i] -= yvel; vz[i] -= zvel;
    }
    
    /* calculate angular momentum */
    for(i=0;i<3;i++){
      imat[i][0] = imat[i][1] = imat[i][2] = 0.;
      invmat[i][0] = invmat[i][1] = invmat[i][2] = 0.;
    }
    l[0] = l[1] = l[2] = 0.;
    for(i=0;i<natoms;i++){
      amt = amass[i];
      dx = px[i];     dy = py[i];     dz = pz[i];
      xvel = vx[i];   yvel = vy[i];   zvel = vz[i];
      
      l[0] += amt*(dy*zvel-dz*yvel);
      l[1] += amt*(dz*xvel-dx*zvel);
      l[2] += amt*(dx*yvel-dy*xvel);
      
      imat[0][0] += amt*(dy*dy+dz*dz);
      imat[1][1] += amt*(dx*dx+dz*dz);
      imat[2][2] += amt*(dx*dx+dy*dy);
      
      imat[0][1] -= amt*dx*dy;
      imat[0][2] -= amt*dx*dz;
      imat[1][2] -= amt*dy*dz;
      
    }
    imat[1][0]=imat[0][1]; imat[2][0]=imat[0][2]; imat[2][1]=imat[1][2];

    invres(3,imat,invmat);
    
    omega[0] = l[0]*invmat[0][0] + l[1]*invmat[1][0] + l[2]*invmat[2][0];
    omega[1] = l[0]*invmat[0][1] + l[1]*invmat[1][1] + l[2]*invmat[2][1];
    omega[2] = l[0]*invmat[0][2] + l[1]*invmat[1][2] + l[2]*invmat[2][2];

    for(i=0;i<natoms;i++){
      xvel = omega[1]*pz[i] - omega[2]*py[i];
      yvel = omega[2]*px[i] - omega[0]*pz[i];
      zvel = omega[0]*py[i] - omega[1]*px[i];
    }

    for(i=0;i<natoms;i++){
      vx[i] -= pz[i]*omega[1] - py[i]*omega[2];
      vy[i] -= px[i]*omega[2] - pz[i]*omega[0];
      vz[i] -= py[i]*omega[0] - px[i]*omega[1];
    }
  }
#ifdef DEBUG
  /* recalculate center of mass, linear and angular momentum */
  l[0] = l[1] = l[2] = 0.;
  cmx = cmy = cmz = 0.;
  xvel = yvel = zvel = 0.;
  for(i=0;i<natoms;i++){
    amt = amass[i];
    dx = px[i]; dy = py[i]; dz = pz[i];
    dvx = vx[i]; dvy = vy[i]; dvz = vz[i];

    cmx += dx*amt;   cmy += dy*amt;   cmz += dz*amt;
    xvel += dvx*amt; yvel += dvy*amt; zvel += dvz*amt;

    l[0] += amt*(dy*dvz-dz*dvy);
    l[1] += amt*(dz*dvx-dx*dvz);
    l[2] += amt*(dx*dvy-dy*dvx);
  }
  cmx /= tmass; cmy /= tmass; cmz /= tmass;
  xvel /= tmass; yvel /= tmass; zvel /= tmass;

  printf("cm = %g %g %g\n",cmx,cmy,cmz);
  printf("vm = %g %g %g\n",xvel,yvel,zvel);
  printf("l  = %g %g %g\n",l[0],l[1],l[2]);
#endif

  free_dmatrix(imat,0,2,0,2);
  free_dmatrix(invmat,0,2,0,2);
}
/*------------------------------------------------------------------*/
void rotate(int natoms,double *px,double *py,double *pz,double *amass)
{
  int i;
  double amt,dx,dy,dz,**imat,d[3];

  md_stdout("Rotating so that axis are principle axis");

  imat = dmatrix(0,2,0,2);
  for(i=0;i<3;i++) {
    imat[i][0] = imat[i][1] = imat[i][2] = 0.;
  }

  for(i=0;i<natoms;i++){
    amt = amass[i];
    dx = px[i]; dy = py[i]; dz = pz[i];
    imat[0][0] += amt*(dy*dy+dz*dz);
    imat[1][1] += amt*(dx*dx+dz*dz);
    imat[2][2] += amt*(dx*dx+dy*dy);
    imat[0][1] -= amt*dx*dy;
    imat[0][2] -= amt*dx*dz;
    imat[1][2] -= amt*dy*dz;
  }
  imat[1][0] = imat[0][1];  imat[2][0] = imat[0][2];  imat[2][1] = imat[1][2];

  rs_me(3,d,imat,1);
  ceigsrtv(d,imat,3); /* we don't need them sorted, but it nice */
 
  /* rotate into priciple frame  */
  for(i=0;i<natoms;i++){
    dx = imat[0][0]*px[i]+imat[1][0]*py[i]+imat[2][0]*pz[i];
    dy = imat[0][1]*px[i]+imat[1][1]*py[i]+imat[2][1]*pz[i];
    dz = imat[0][2]*px[i]+imat[1][2]*py[i]+imat[2][2]*pz[i];
    px[i] = dx; py[i] = dy; pz[i] = dz;
  }
  free_dmatrix(imat,0,2,0,2);
}
/*------------------------------------------------------------------*/
