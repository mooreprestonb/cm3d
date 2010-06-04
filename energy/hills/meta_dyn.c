/* 
   Copyright (C) 1997 1998 1999 2000 Dr. Preston B. Moore

Dr. Preston B. Moore
Associate Director, Center for Molecular Modeling (CMM)
University of Pennsylvania, Department of Chemistry, Box 188 
231 S. 34th St. Philadelphia, PA 19104-6323 USA
EMAIL: moore@cmm.chem.upenn.edu  
WWW: http://www.cmm.upenn.edu/~moore
#ifdef PRINTFUNCTION

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

/* subroutines dealing with bonded interactions
   search_bond_base - searches the data base and sets up initial vectors
   fbond - gets the force of the bond
   getvbond - gets the potential energy of the bond
   bonded   - check to see if to atoms are bonded
*/

#define DEBUGOFF
#define PRINTFUNCTIONOFF
#define PRINTTORSIONOF

#include "md.h"

/* #define PRINT_COLVAR */

/* rescale colvar velocities to target temperature  */
/* if temperature larger than temp+thrmmass         */
void scale_colvar(COLVAR *colvar)
{
  int icol;
  double t_var;

  for(icol=0;icol<colvar->ncolvar;++icol){
    t_var = colvar->amass[icol]*colvar->vcolvar[icol]*colvar->vcolvar[icol];
    if(t_var > colvar->temp[icol]+colvar->thrmmass[icol]) {
      colvar->vcolvar[icol] = COPY_SIGN(sqrt(colvar->temp[icol]/colvar->amass[icol]),colvar->vcolvar[icol]);
    }
  }
}

/*---------------------------------------------------------------------*/
/* forces on atoms due to spring of cv icol                            */
double fcolvar_dzplane(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double dx,dy,dz,r,dr;
  double zcomi,zcomj,tmassi,tmassj,mass;
  double *amass = coords->amass;
  double func,dfunc;

  /* if(coords->colvar.ndim[icol] !=2) error(); */

  int nicolvar = coords->colvar.natom[icol][0];  /* i= com of lipid bilayer */
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];  /* j= object in lipid bilayer */
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;
  /*decomp1d((long)nbonds,simparms->size,simparms->rank,&ibegin,&iend);*/

  zcomi = 0.0;
  zcomj = 0.0;
  tmassi = tmassj = 0.0;

  /* center of mass of atoms from lipid bilayer */
  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    mass = amass[ii];
    zcomi += mass*pz[ii];
    tmassi += mass;
  }
  zcomi /= tmassi;

  /* center of mass of object in or above lipid bilayer */
  for(j=0;j<njcolvar;++j){
    jj = jcolvar[j];
    mass = amass[jj];
    dx=px[jj];
    dy=py[jj];
    dz=pz[jj];
    period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,3,0); 
    zcomj += mass*dz;
    tmassj += mass;
  }
  zcomj /= tmassj;
  
  /* heigth difference between centers of mass */ 
  r = zcomi - zcomj;
  dr = r - coords->colvar.pcolvar[icol];

  func = .5*fk*dr*dr;
  dfunc = fk*dr;
    
  /* change the force */
  if(iforce && (simparms->rank==0)){
    for(i=0;i<nicolvar;++i){
      ii = icolvar[i];
      fz[ii] -= dfunc*amass[ii]/tmassi;
    }
    for(i=0;i<njcolvar;++i){
      jj = jcolvar[i];
      fz[jj] += dfunc*amass[jj]/tmassj;
    }
  }

  coords->colvar.pistcolvar[icol] = r;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];
  return func;
}

/*---------------------------------------------------------------------*/
/* forces on atoms due to spring of cv icol                            */
double fcolvar_dzplaneabs(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double r,dr;
  double zcomi,zcomj,tmassi,tmassj,mass;
  double *amass = coords->amass;
  double func,dfunc;

  /* if(coords->colvar.ndim[icol] !=2) error(); */

  int nicolvar = coords->colvar.natom[icol][0];  /* i= com of lipid bilayer */
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];  /* j= object in lipid bilayer */
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];

  double *pz = coords->pz;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;
  /*decomp1d((long)nbonds,simparms->size,simparms->rank,&ibegin,&iend);*/

  zcomi = 0.0;
  zcomj = 0.0;
  tmassi = tmassj = 0.0;


  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    mass = amass[ii];
    zcomi += mass*pz[ii];
    tmassi += mass;
  }
  zcomi /= tmassi;

  
  /* heigth difference between centers of mass */ 
  r = zcomi - zcomj;
  dr = r - coords->colvar.pcolvar[icol];

  func = .5*fk*dr*dr;
  dfunc = fk*dr;
    
  /* change the force */
  if(iforce && (simparms->rank==0)){
    for(i=0;i<nicolvar;++i){
      ii = icolvar[i];
      fz[ii] -= dfunc*amass[ii]/tmassi;
    }
  }

  coords->colvar.pistcolvar[icol] = r;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];
  return func;
}


/*---------------------------------------------------------------------*/
/* forces on atoms due to spring of cv icol                            */
double fcolvar_bond(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double dx,dy,dz,r,fr,dr;
  double xcomi,ycomi,zcomi,xcomj,ycomj,zcomj,tmassi,tmassj,mass;
  double *amass = coords->amass;
  double func,dfunc;

  /* if(coords->colvar.ndim[icol] !=2) error(); */

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;
  /*decomp1d((long)nbonds,simparms->size,simparms->rank,&ibegin,&iend);*/

  xcomi = ycomi = zcomi = 0.0;
  xcomj = ycomj = zcomj = 0.0;
  tmassi = tmassj = 0.0;
  fr = 0.;

  /* center of mass on one end of the bond */
  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    mass = amass[ii];
    xcomi += mass*px[ii];
    ycomi += mass*py[ii];
    zcomi += mass*pz[ii];
    tmassi += mass;
  }
  xcomi /= tmassi;
  ycomi /= tmassi;
  zcomi /= tmassi;

  /* center of mass on other end of the bond */
  for(j=0;j<njcolvar;++j){
    jj = jcolvar[j];
    mass = amass[jj];
    xcomj += mass*px[jj];
    ycomj += mass*py[jj];
    zcomj += mass*pz[jj];
    tmassj += mass;
  }
  xcomj /= tmassj;
  ycomj /= tmassj;
  zcomj /= tmassj;
  
  /* distance between centers of mass */
  dx = xcomi - xcomj;  
  dy = ycomi - ycomj;  
  dz = zcomi - zcomj;

  period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,3,0); 
  r = sqrt(dx*dx+dy*dy+dz*dz);
  dr = r - coords->colvar.pcolvar[icol];

  func = .5*fk*dr*dr;
  dfunc = fk*dr;
    
  /* change the force */
  if(iforce && (simparms->rank==0)){
    for(i=0;i<nicolvar;++i){
      ii = icolvar[i];
      fx[ii] -= dfunc*dx*amass[ii]/(r*tmassi);
      fy[ii] -= dfunc*dy*amass[ii]/(r*tmassi);
      fz[ii] -= dfunc*dz*amass[ii]/(r*tmassi);
    }
    for(i=0;i<njcolvar;++i){
      jj = jcolvar[i];
      fx[jj] += dfunc*dx*amass[jj]/(r*tmassj);
      fy[jj] += dfunc*dy*amass[jj]/(r*tmassj);
      fz[jj] += dfunc*dz*amass[jj]/(r*tmassj);
    }
  }

  coords->colvar.pistcolvar[icol] = r;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];
  return func;
}

/*---------------------------------------------------------------------*/
/* forces due to tethering COM to origin                         */
double fcolvar_COMorigin(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double r,dr;
  double xcomi,ycomi,zcomi,xcomj,ycomj,zcomj,tmassi,tmassj,mass;
  double *amass = coords->amass;
  double func,dfunc;

  /* if(coords->colvar.ndim[icol] !=2) error(); */

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;

  xcomi = ycomi = zcomi = 0.0;
  xcomj = ycomj = zcomj = 0.0;
  tmassi = tmassj = 0.0;

  /* center of mass (j is a dummy variable here, only i matters.) */
  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    mass = amass[ii];
    xcomi += mass*px[ii];
    ycomi += mass*py[ii];
    zcomi += mass*pz[ii];
    tmassi += mass;
  }
  xcomi /= tmassi;
  ycomi /= tmassi;
  zcomi /= tmassi;

  r = sqrt(xcomi*xcomi + ycomi*ycomi + zcomi*zcomi);
  func = .5*fk*r*r;
  dfunc = fk*r;

  /* change the force */
  if(iforce && (simparms->rank==0)){
    for(i=0;i<nicolvar;++i){
      ii = icolvar[i];
      fx[ii] -= fk*xcomi/(double)nicolvar;
      fy[ii] -= fk*ycomi/(double)nicolvar;
      fz[ii] -= fk*zcomi/(double)nicolvar;
    }
  }

  coords->colvar.pistcolvar[icol] = r;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];

if(iforce){
#ifdef PARA
        coords->colvar.meanforce_cum += dfunc/((double)simparms->size);
#else
        coords->colvar.meanforce_inst += dfunc;
        coords->colvar.meanforce_cum += dfunc;
#endif
}

  return func;
}


/*---------------------------------------------------------------------*/

/* forces due to radius of gyration constraint                         */
double fcolvar_gyration(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double r,dr;
  double xcomi,ycomi,zcomi,xcomj,ycomj,zcomj,tmassi,tmassj,mass,Rg2;
  double *amass = coords->amass;
  double func,dfunc;

  /* if(coords->colvar.ndim[icol] !=2) error(); */

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;

  xcomi = ycomi = zcomi = 0.0;
  xcomj = ycomj = zcomj = 0.0;
  tmassi = tmassj = 0.0;

  /* center of mass (j is a dummy variable here, only i matters.) */
  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    mass = amass[ii];
    xcomi += mass*px[ii];
    ycomi += mass*py[ii];
    zcomi += mass*pz[ii];
    tmassi += mass;
  }
  xcomi /= tmassi;
  ycomi /= tmassi;
  zcomi /= tmassi;

  Rg2 = 0.0;      /* radius of gyration squared */           
  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    Rg2 += (px[ii] - xcomi)*(px[ii] - xcomi) + (py[ii] - ycomi)*(py[ii] - ycomi) + (pz[ii] - zcomi)*(pz[ii] - zcomi);
  }
  Rg2 /= (double)nicolvar;
  dr = Rg2 - coords->colvar.pcolvar[icol];
  func = .5*fk*dr*dr;
  dfunc = fk*dr;
    
  /* change the force */
  if(iforce && (simparms->rank==0)){
    for(i=0;i<nicolvar;++i){
      ii = icolvar[i];
      fx[ii] -= 2.0*dfunc*(px[ii]-xcomi)/(double)nicolvar;
      fy[ii] -= 2.0*dfunc*(py[ii]-ycomi)/(double)nicolvar;
      fz[ii] -= 2.0*dfunc*(pz[ii]-zcomi)/(double)nicolvar;
    }
  }

  coords->colvar.pistcolvar[icol] = Rg2;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];

if(iforce){
#ifdef PARA
        coords->colvar.meanforce_cum += dfunc/((double)simparms->size);
#else
        coords->colvar.meanforce_inst += dfunc;
        coords->colvar.meanforce_cum += dfunc;
#endif
}

  return func;
}


/*---------------------------------------------------------------------*/

/* forces due to cluster algorithm                         */
double fcolvar_cluster(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,k,ii,jj,flag,iclosest,jclosest,cluster_num,clusterloop;
  double r,dr,dx,dy,dz,rmin,dxmin,dymin,dzmin;
  double *amass = coords->amass;
  double func,dfunc;

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];
  int *cluster;
  cluster = (int*) malloc(sizeof(int)*nicolvar);

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;


/* I'm assuming that all particles are in list i, and I'm ignoring j */

  cluster[0] = icolvar[0]; /* pick the first member of the cluster */
 cluster_num = 1;   /* number of elements in the cluster  */
for(clusterloop=0;clusterloop<(nicolvar-1);++clusterloop){
 rmin=9999.99;
 for(i=0;i<cluster_num;++i){
  for(j=0;j<nicolvar;++j){
   jj = icolvar[j];
   flag = 0; 
   for(k=0;k<cluster_num;++k){
    if(jj == cluster[k]) flag = 1;
   }
  if(flag == 0){
   dx = px[cluster[i]] - px[jj];
   dy = py[cluster[i]] - py[jj];
   dz = pz[cluster[i]] - pz[jj];
   period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,3,0);
   r = sqrt(dx*dx+dy*dy+dz*dz);
   if(r < rmin ){
    rmin = r;
    dxmin = dx;
    dymin = dy;
    dzmin = dz;
    iclosest = cluster[i];
    jclosest = jj;
   }
  }
  }
 }
/* fprintf(stdout," rmin %g %d %d \n" ,rmin,iclosest,jclosest); */
 if(iforce && (simparms->rank==0)){
if(rmin > 12.3){
 fx[iclosest] -= 10000.0 * 3.0 * (rmin-12.3)*(rmin-12.3) * dxmin/rmin;
 fy[iclosest] -= 10000.0 * 3.0 * (rmin-12.3)*(rmin-12.3) * dymin/rmin;
 fz[iclosest] -= 10000.0 * 3.0 * (rmin-12.3)*(rmin-12.3) * dzmin/rmin;
 fx[jclosest] += 10000.0 * 3.0 * (rmin-12.3)*(rmin-12.3) * dxmin/rmin;
 fy[jclosest] += 10000.0 * 3.0 * (rmin-12.3)*(rmin-12.3) * dymin/rmin;
 fz[jclosest] += 10000.0 * 3.0 * (rmin-12.3)*(rmin-12.3) * dzmin/rmin;
}
}
 /* Now add this jclosest to the cluster and repeat */
cluster_num++;
cluster[cluster_num-1] = jclosest;
}
 
    
  func= 0.0;
  dfunc= 0.0;
  coords->colvar.pistcolvar[icol] = 0.0;
  f_spring[icol] = func;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];

if(iforce){
#ifdef PARA
        coords->colvar.meanforce_cum += dfunc/((double)simparms->size);
#else
        coords->colvar.meanforce_inst += dfunc;
        coords->colvar.meanforce_cum += dfunc;
#endif
}

  return func;
}
/*---------------------------------------------------------------------*/
/* forces on atoms due to spring of cv icol                            */
double fcolvar_clusteradd(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj,jclosest;
  /* long ibegin,iend; */
  double dx,dy,dz,r,dr,rmin,dxmin,dymin,dzmin;
  double *amass = coords->amass;
  double func,dfunc;

  /* assume the single particle is in 'i', and the existing cluster is in 'j' */

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;

 ii = icolvar[0];  /* the single particle */
 rmin = 9999.99;
  for(j=0;j<njcolvar;++j){
    jj = jcolvar[j];
    dx = px[ii] - px[jj];
    dy = py[ii] - py[jj];
    dz = pz[ii] - pz[jj];
   period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,3,0);
   r = sqrt(dx*dx+dy*dy+dz*dz);
   if(r < rmin ){
    rmin = r;
    dxmin = dx;
    dymin = dy;
    dzmin = dz;
    jclosest = jj;
   }
  }
  
  dr = rmin - coords->colvar.pcolvar[icol];

  func = .5*fk*dr*dr;
  dfunc = fk*dr;
    
  /* change the force */
  if(iforce && (simparms->rank==0)){
      ii = icolvar[0];
      fx[ii] -= dfunc*dxmin/rmin;
      fy[ii] -= dfunc*dymin/rmin;
      fz[ii] -= dfunc*dzmin/rmin;
      jj = jclosest;
      fx[jj] += dfunc*dxmin/rmin;
      fy[jj] += dfunc*dymin/rmin;
      fz[jj] += dfunc*dzmin/rmin;
  }

  coords->colvar.pistcolvar[icol] = rmin;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];
  return func;
}
/*---------------------------------------------------------------------*/
/* forces due to coordination number min. constraint                         */
double fcolvar_coordmin(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double r,dr,dx,dy,dz,t1,Ni;
  double *amass = coords->amass;
  double func,dfunc;

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];
  double rij[njcolvar];
  double xij[njcolvar];
  double yij[njcolvar];
  double zij[njcolvar];
  double dfrijdrij[njcolvar];

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;


/* I'm assuming that i is a single particle, and j is the rest of the cluster */

 Ni = 0.0;
 ii = icolvar[0];
  for(j=0;j<njcolvar;++j){
    jj = jcolvar[j];
    dx = px[ii] - px[jj];
    dy = py[ii] - py[jj];
    dz = pz[ii] - pz[jj];
    xij[j]=dx;
    yij[j]=dy;
    zij[j]=dz;
    period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,3,0);
    rij[j] = sqrt(dx*dx+dy*dy+dz*dz);
    t1 = exp(4.0*(rij[j]-13.0));
    Ni += 1.0/(t1+1.0);
    dfrijdrij[j] = -4.0 * t1 / ( (1.0+t1)*(1.0+t1) );
  }
 if(iforce && (simparms->rank==0)){
  for(j=0;j<njcolvar;++j){   /* force on particle i */
    jj = jcolvar[j];
    fx[ii] += 000.0 * 5.0 * exp(-5.0*(Ni-0.68)) * dfrijdrij[j] * xij[j]/rij[j];
    fy[ii] += 000.0 * 5.0 * exp(-5.0*(Ni-0.68)) * dfrijdrij[j] * yij[j]/rij[j];
    fz[ii] += 000.0 * 5.0 * exp(-5.0*(Ni-0.68)) * dfrijdrij[j] * zij[j]/rij[j];
  }
  for(j=0;j<njcolvar;++j){   /* force on particles j */
    jj = jcolvar[j];
    fx[jj] -= 000.0 * 5.0 * exp(-5.0*(Ni-0.68)) * dfrijdrij[j] * xij[j]/rij[j];
    fy[jj] -= 000.0 * 5.0 * exp(-5.0*(Ni-0.68)) * dfrijdrij[j] * yij[j]/rij[j];
    fz[jj] -= 000.0 * 5.0 * exp(-5.0*(Ni-0.68)) * dfrijdrij[j] * zij[j]/rij[j];
  }
 }
 
    
  func= 000.0 * exp(-5.0*(Ni-0.68));
  coords->colvar.pistcolvar[icol] = Ni;
  f_spring[icol] = func;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];

if(iforce){
#ifdef PARA
        coords->colvar.meanforce_cum += dfunc/((double)simparms->size);
#else
        coords->colvar.meanforce_inst += dfunc;
        coords->colvar.meanforce_cum += dfunc;
#endif
}

  return func;
}


/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/* forces due to tether to sphere surface                              */
double fcolvar_tether(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double dx,dy,dz,r,dr;
  double xcomi,ycomi,zcomi,xcomj,ycomj,zcomj,tmassi,tmassj,mass;
  double *amass = coords->amass;
  double func,dfunc;

  /* if(coords->colvar.ndim[icol] !=2) error(); */

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;
  /*decomp1d((long)nbonds,simparms->size,simparms->rank,&ibegin,&iend);*/

  xcomi = ycomi = zcomi = 0.0;
  xcomj = ycomj = zcomj = 0.0;
  tmassi = tmassj = 0.0;

  /* center of mass on one end of the bond */
  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    mass = amass[ii];
    xcomi += mass*px[ii];
    ycomi += mass*py[ii];
    zcomi += mass*pz[ii];
    tmassi += mass;
  }
  xcomi /= tmassi;
  ycomi /= tmassi;
  zcomi /= tmassi;

  /* center of mass on other end of the bond */
  for(j=0;j<njcolvar;++j){
    jj = jcolvar[j];
    mass = amass[jj];
    xcomj += mass*px[jj];
    ycomj += mass*py[jj];
    zcomj += mass*pz[jj];
    tmassj += mass;
  }
  xcomj /= tmassj;
  ycomj /= tmassj;
  zcomj /= tmassj;
  
  /* distance between centers of mass */
  dx = xcomi - xcomj;  
  dy = ycomi - ycomj;  
  dz = zcomi - zcomj;

  period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,3,0); 
  r = sqrt(dx*dx+dy*dy+dz*dz);
  dr = r - coords->colvar.pcolvar[0] - 3.0 - coords->colvar.pcolvar[icol];

  func = .5*fk*dr*dr;
  dfunc = fk*dr;
    
  /* change the force */
  if(iforce && (simparms->rank==0)){
    for(i=0;i<nicolvar;++i){
      ii = icolvar[i];
      fx[ii] -= dfunc*dx*amass[ii]/(r*tmassi);
      fy[ii] -= dfunc*dy*amass[ii]/(r*tmassi);
      fz[ii] -= dfunc*dz*amass[ii]/(r*tmassi);
    }
    for(i=0;i<njcolvar;++i){
      jj = jcolvar[i];
      fx[jj] += dfunc*dx*amass[jj]/(r*tmassj);
      fy[jj] += dfunc*dy*amass[jj]/(r*tmassj);
      fz[jj] += dfunc*dz*amass[jj]/(r*tmassj);
    }
  }

  if (coords->colvar.amass[0]!=0.0)
  coords->colvar.fcolvar[0] += dfunc/coords->colvar.amass[0];

  coords->colvar.pistcolvar[icol] = r - coords->colvar.pcolvar[0] - 3.0;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];

if(iforce){
#ifdef PARA
        coords->colvar.meanforce_cum += dfunc/((double)simparms->size);
#else
        coords->colvar.meanforce_inst += dfunc;
        coords->colvar.meanforce_cum += dfunc;
#endif
}

  return func;
}


/*---------------------------------------------------------------------*/
/* forces due to Morse tether to sphere surface                              */
double fcolvar_morsetether(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  int i,j,ii,jj;
  /* long ibegin,iend; */
  double dx,dy,dz,r,dr;
  double xcomi,ycomi,zcomi,xcomj,ycomj,zcomj,tmassi,tmassj,mass;
  double *amass = coords->amass;
  double func,dfunc;

  int nicolvar = coords->colvar.natom[icol][0];
  int *icolvar = coords->colvar.atom[icol][0];
  int njcolvar = coords->colvar.natom[icol][1];
  int *jcolvar = coords->colvar.atom[icol][1];
  double fk = coords->colvar.fk[icol];  /* this is the force constant of the harmonic fit to the Morse */
  double De = coords->colvar.mindwidth[icol];  /* define De in the mindwidth field  */
  double alpha = sqrt(fk/(2.0*De));    /* define alpha in terms of De and k  */

  double *px = coords->px;
  double *py = coords->py;
  double *pz = coords->pz;

  double *fx = coords->fxa;
  double *fy = coords->fya;
  double *fz = coords->fza;

  double *f_spring = coords->colvar.f_spring;
  /*decomp1d((long)nbonds,simparms->size,simparms->rank,&ibegin,&iend);*/

  xcomi = ycomi = zcomi = 0.0;
  xcomj = ycomj = zcomj = 0.0;
  tmassi = tmassj = 0.0;

  /* center of mass on one end of the bond */
  for(i=0;i<nicolvar;++i){
    ii = icolvar[i];
    mass = amass[ii];
    xcomi += mass*px[ii];
    ycomi += mass*py[ii];
    zcomi += mass*pz[ii];
    tmassi += mass;
  }
  xcomi /= tmassi;
  ycomi /= tmassi;
  zcomi /= tmassi;

  /* center of mass on other end of the bond */
  for(j=0;j<njcolvar;++j){
    jj = jcolvar[j];
    mass = amass[jj];
    xcomj += mass*px[jj];
    ycomj += mass*py[jj];
    zcomj += mass*pz[jj];
    tmassj += mass;
  }
  xcomj /= tmassj;
  ycomj /= tmassj;
  zcomj /= tmassj;
  
  /* distance between centers of mass */
  dx = xcomi - xcomj;  
  dy = ycomi - ycomj;  
  dz = zcomi - zcomj;

  period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,3,0); 
  r = sqrt(dx*dx+dy*dy+dz*dz);
  dr = r - coords->colvar.pcolvar[0] - 3.0 - coords->colvar.pcolvar[icol];

  func = De*(1.0-exp(-alpha*dr))*(1.0-exp(-alpha*dr));
  dfunc = 2.0*De*alpha*exp(-alpha*dr)*(1.0-exp(-alpha*dr));
    
  /* change the force */
  if(iforce && (simparms->rank==0)){
    for(i=0;i<nicolvar;++i){
      ii = icolvar[i];
      fx[ii] -= dfunc*dx*amass[ii]/(r*tmassi);
      fy[ii] -= dfunc*dy*amass[ii]/(r*tmassi);
      fz[ii] -= dfunc*dz*amass[ii]/(r*tmassi);
    }
    for(i=0;i<njcolvar;++i){
      jj = jcolvar[i];
      fx[jj] += dfunc*dx*amass[jj]/(r*tmassj);
      fy[jj] += dfunc*dy*amass[jj]/(r*tmassj);
      fz[jj] += dfunc*dz*amass[jj]/(r*tmassj);
    }
  }

  if (coords->colvar.amass[0]!=0.0)
  coords->colvar.fcolvar[0] += dfunc/coords->colvar.amass[0];

  coords->colvar.pistcolvar[icol] = r - coords->colvar.pcolvar[0] - 3.0;
  f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += f_spring[icol];
  coords->colvar.f_spring_cum[icol] += f_spring[icol];

if(iforce){
#ifdef PARA
        coords->colvar.meanforce_cum += dfunc/((double)simparms->size);
#else
        coords->colvar.meanforce_inst += dfunc;
        coords->colvar.meanforce_cum += dfunc;
#endif
}

  return func;
}

/*---------------------------------------------------------------------*/


double fcolvar_sphere_R(int icol,SIMPARMS *simparms,COORDS *coords,int iforce)
{
  /*int ncyl_idx = coords->colvar.natom[icol][0];*/
  int nsph_idx = coords->colvar.natom[icol][0];
  double fk = coords->colvar.fk[icol];
  double dr,func,dfunc;
  double r = coords->colvar.pcolvar[nsph_idx]; /* radius of sphere */

  dr = r - coords->colvar.pcolvar[icol];
 
  func = .5*fk*dr*dr;
  dfunc = fk*dr;

  /* change the force */
  if(iforce && (simparms->rank==0)){
    if (coords->colvar.amass[nsph_idx]!=0.0) /* if sphere mass set to be 0, mass = infinity */
      coords->colvar.fcolvar[nsph_idx] += -dfunc/coords->colvar.amass[nsph_idx];
  }

  coords->colvar.pistcolvar[icol] = r;
  coords->colvar.f_spring[icol] = dfunc;
  coords->colvar.fcolvar[icol] += coords->colvar.f_spring[icol];
  coords->colvar.f_spring_cum[icol] += coords->colvar.f_spring[icol];

  return func;
}

/*---------------------------------------------------------------------*/

void fcolvar_cntr(SIMPARMS *simparms,COORDS *coords,int iforce)
{              /* if iforce=0 forces on atoms are not applied    */
               /* for initialization called from input/colvar.c  */
  int ncolvar = coords->colvar.ncolvar;
  int *type = coords->colvar.type;
  int *lhills = coords->colvar.lhills;
  int *lconst_vel = coords->colvar.lconst_vel;
  int lthills = 0;
  int icol;
  double dxmax,dxmin,dx;
  char line[MAXLINELEN];

  if(ncolvar==0) return;
  /* get forces on atoms from harmonic springs with colvars */
  coords->colvar.ecolvar = 0.0;

  for(icol=0;icol<ncolvar;++icol){
    switch(type[icol]) {
    case 0: /* distance */
      coords->colvar.ecolvar += fcolvar_bond(icol,simparms,coords,iforce);
      break;
    case 5: /* dzplane */
      coords->colvar.ecolvar += fcolvar_dzplane(icol,simparms,coords,iforce);
      break;
    case 6: /* dzplaneabs */
      coords->colvar.ecolvar += fcolvar_dzplaneabs(icol,simparms,coords,iforce);
      break;
    case 11: /* sphere */
      coords->colvar.ecolvar = 0.0; /* the potential of sphere already exist*/ 
      break;
    case 12: /* sphere colvar*/
      coords->colvar.ecolvar += fcolvar_sphere_R(icol,simparms,coords,iforce);
      break;
    case 13: /* tether*/
      coords->colvar.ecolvar += fcolvar_tether(icol,simparms,coords,iforce);
      break;
    case 15: /* gyration*/
      coords->colvar.ecolvar += fcolvar_gyration(icol,simparms,coords,iforce);
      break;
    case 16: /* COMorigin*/
      coords->colvar.ecolvar += fcolvar_COMorigin(icol,simparms,coords,iforce);
      break;
    case 17: /* coordmin*/
      coords->colvar.ecolvar += fcolvar_coordmin(icol,simparms,coords,iforce);
      break;
    case 18: /* cluster*/
      coords->colvar.ecolvar += fcolvar_cluster(icol,simparms,coords,iforce);
      break;
    case 19: /* clusteradd*/
      coords->colvar.ecolvar += fcolvar_clusteradd(icol,simparms,coords,iforce);
      break;
    case 20: /* Morse tether*/
      coords->colvar.ecolvar += fcolvar_morsetether(icol,simparms,coords,iforce);
      break;
    default:
      md_error("in meta_dyn colvar type not known!");
    }      
    if(lhills[icol]){
      lthills=1;
      if(coords->colvar.ltunehills==2) coords->colvar.fcolvar[icol] = 0.0;
    }
    if(lconst_vel[icol]) coords->colvar.fcolvar[icol] = 0.0;
  }

  /* Reset pcolvar to pistcolvar and zero velocities and forces */
  if(iforce==0){
    for(icol=0;icol<ncolvar;++icol){
      if(coords->colvar.pcolvar[icol] == RARE){
	coords->colvar.pcolvar[icol] = coords->colvar.pistcolvar[icol];
      }
      if(coords->colvar.lconst_vel[icol]==0){
	coords->colvar.vcolvar[icol] = 0.0;
      }
      coords->colvar.fcolvar[icol] = 0.0;
      coords->colvar.f_hills[icol] = 0.0;
      coords->colvar.f_spring[icol] = 0.0;
      coords->colvar.f_spring_cum[icol] = 0.0;
    }
  }

  /* put hills and get forces on colvars from hills */

  if(lthills) hills(simparms,coords);
  /* scale velocities if we are too hot */
  scale_colvar(&coords->colvar);

  /* set reflective boundaries */

  for(icol=0;icol<ncolvar;++icol){
    dxmin = coords->colvar.pcolvar[icol] - coords->colvar.min_val[icol];
    dxmax = coords->colvar.pcolvar[icol] - coords->colvar.max_val[icol];
    dx = MIN(fabs(dxmax),fabs(dxmin));
    if( dx > (coords->colvar.max_val[icol]-coords->colvar.min_val[icol]) ){
      sprintf(line,"Boundaries for coll.var. %d to close together.",icol);
      md_warning(line);
      coords->colvar.pcolvar[icol]=coords->colvar.max_val[icol];
      coords->colvar.vcolvar[icol]=0.0;
      dxmin=0.0;
      dxmax=0.0;
    }
    if(dxmin < 0.0 ){
      coords->colvar.pcolvar[icol] -= 2.0*dxmin ;
      coords->colvar.vcolvar[icol] *= -1.0 ;

/* Adding a hill at the boundary is not a good idea if the colvar  */
/* diffuses slowly. This should be turned into an optional feature */
/*      if(lthills){
	for(icol=0;icol<ncolvar;++icol){
	  coords->colvar.pcolvarlast[icol] = coords->colvar.pcolvar[icol];
	}
	add_hills(&(coords->colvar),simparms->istep);
	coords->colvar.lasthill = simparms->istep;
      }
*/
    }
    else if( dxmax > 0.0 ){
      coords->colvar.pcolvar[icol] -= 2.0*dxmax ;
      coords->colvar.vcolvar[icol] *= -1.0 ;

/*
      if(lthills){
	for(icol=0;icol<ncolvar;++icol){
	  coords->colvar.pcolvarlast[icol] = coords->colvar.pcolvar[icol];
	}
	add_hills(&(coords->colvar),simparms->istep);
	coords->colvar.lasthill = simparms->istep;
      }
*/
    }
  }
  
  /* change force to acel */

  for(icol=0;icol<ncolvar;++icol){
    if(coords->colvar.type[icol]!=11)
      coords->colvar.fcolvar[icol] /= coords->colvar.amass[icol];
  }
}
/*---------------------------------------------------------------------*/


