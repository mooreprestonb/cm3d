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

/* subroutines dealing with bonded interactions
   search_xbond_base - searches the data base and sets up initial vectors
   fbondx - gets the force of the bond
   getvbondx - gets the potential energy of the bond
   bondxed   - check to see if to atoms are bonded
*/

#include "md.h"

static int ibonx,jbonx,kbonx;
static double eq1bonx,eq2bonx,fkbonx;
/*-----------------------------------------------------------------------*/

void fbondx(SIMPARMS *simparms,COORDS *coords)
{
  int i,j,k,itype,nb;
  long ibegin,iend;
  double ax,ay,az,al,bx,by,bz,bl,dr1,dr2,fr1,fr2;
  double *px,*py,*pz,*fx,*fy,*fz;
  int nbondxs,*idxbox,*ibondx,*jbondx,*kbondx;
  double *eq1bondx,*eq2bondx,*fkbondx;

  nbondxs = coords->bondxs.nbondxs;
  idxbox = coords->bondxs.idxbox;
  ibondx = coords->bondxs.ibondx;
  jbondx = coords->bondxs.jbondx;
  kbondx = coords->bondxs.kbondx;
  eq1bondx = coords->bondxs.eq1bondx;
  eq2bondx = coords->bondxs.eq2bondx;
  fkbondx  = coords->bondxs.fkbondx;

  px = coords->px;
  py = coords->py;
  pz = coords->pz;
  fx = coords->fxa;
  fy = coords->fya;
  fz = coords->fza;

  decomp1d((long)nbondxs,simparms->size,simparms->rank,&ibegin,&iend);
  for(nb=ibegin;nb<iend;nb++){
    i = ibondx[nb];
    j = jbondx[nb];
    k = kbondx[nb];
    itype = idxbox[nb];

    ax = px[i] - px[j];
    ay = py[i] - py[j];
    az = pz[i] - pz[j];
    al = sqrt(ax*ax+ay*ay+az*az);

    bx = px[j] - px[k];
    by = py[j] - py[k];
    bz = pz[j] - pz[k];
    bl = sqrt(bx*bx+by*by+bz*bz);

    dr1 = al-eq1bondx[itype];
    dr2 = bl-eq2bondx[itype];

    fr1 = fkbondx[itype]*dr2/al;
    fr2 = fkbondx[itype]*dr1/bl;

    fx[i] += -fr1*ax;
    fy[i] += -fr1*ay;
    fz[i] += -fr1*az;
    fx[j] += fr1*ax-fr2*bx;
    fy[j] += fr1*ay-fr2*by;
    fz[j] += fr1*az-fr2*bz;
    fx[k] += fr2*bx;
    fy[k] += fr2*by;
    fz[k] += fr2*bz;
  }
}
/*-------------------------------------------------------------------*/
void getvbondx(SIMPARMS *simparms,COORDS *coords,
	       double *potbox,double *prsbox,double ptensor[9])
{
  int i,j,k,itype,nb;
  long ibegin,iend;
  double ax,ay,az,al,bx,by,bz,bl,dr1,dr2;
  double spot=0.,sprs=0.,sprsl,spotl;
  double sprstensor[9],sprstl[9];
  int nbondxs,*idxbox,*ibondx,*jbondx,*kbondx;
  double *eq1bondx,*eq2bondx,*fkbondx;

  for(i=0;i<9;i++) sprstensor[i]=0.0;

  nbondxs = coords->bondxs.nbondxs;
  idxbox = coords->bondxs.idxbox;
  ibondx = coords->bondxs.ibondx;
  jbondx = coords->bondxs.jbondx;
  kbondx = coords->bondxs.kbondx;
  eq1bondx = coords->bondxs.eq1bondx;
  eq2bondx = coords->bondxs.eq2bondx;
  fkbondx  = coords->bondxs.fkbondx;

  decomp1d((long)nbondxs,simparms->size,simparms->rank,&ibegin,&iend);
  for(nb=ibegin;nb<iend;nb++){
    i = ibondx[nb];
    j = jbondx[nb];
    k = kbondx[nb];
    itype = idxbox[nb];

    ax = coords->px[i] - coords->px[j];
    ay = coords->py[i] - coords->py[j];
    az = coords->pz[i] - coords->pz[j];
    al = sqrt(ax*ax+ay*ay+az*az);

    bx = coords->px[j] - coords->px[k];
    by = coords->py[j] - coords->py[k];
    bz = coords->pz[j] - coords->pz[k];
    bl = sqrt(bx*bx+by*by+bz*bz);

    dr1 = al-eq1bondx[itype];
    dr2 = bl-eq2bondx[itype];

    spot += fkbondx[itype]*dr1*dr2;
    sprs -= fkbondx[itype]*(bl*dr1+al*dr2);

    sprstensor[0] -= fkbondx[itype]*(dr1*bx*bx/bl+dr2*ax*ax/al);
    sprstensor[1] -= fkbondx[itype]*(dr1*by*bx/bl+dr2*ay*ax/al);
    sprstensor[2] -= fkbondx[itype]*(dr1*bz*bx/bl+dr2*az*ax/al);

    sprstensor[3] -= fkbondx[itype]*(dr1*bx*by/bl+dr2*ax*ay/al);
    sprstensor[4] -= fkbondx[itype]*(dr1*by*by/bl+dr2*ay*ay/al);
    sprstensor[5] -= fkbondx[itype]*(dr1*bz*by/bl+dr2*az*ay/al);

    sprstensor[6] -= fkbondx[itype]*(dr1*bx*bz/bl+dr2*ax*az/al);
    sprstensor[7] -= fkbondx[itype]*(dr1*by*bz/bl+dr2*ay*az/al);
    sprstensor[8] -= fkbondx[itype]*(dr1*bz*bz/bl+dr2*az*az/al);
  }
#ifdef PARA
  MPI_Allreduce(&spot,&spotl,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&sprs,&sprsl,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(sprstensor,sprstl,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
  spotl = spot;
  sprsl = sprs;
  for(i=0;i<9;i++) sprstl[i] = sprstensor[i];
#endif
  *potbox += spotl;
  *prsbox += sprsl;
  for(i=0;i<9;i++) ptensor[i] += sprstl[i];
}
/*---------------------------------------------------------------*/
void fcbondx(COORDS *coords,double **fcmat)
{
  int i,j,k,ii,jj,nb,itype;
  double dr1,dr2,fk;
  double ax,ay,az,ra3,al,ral,bx,by,bz,rbl,rb3,bl;
  double dr1dx1,dr1dy1,dr1dz1,dr1dx2,dr1dy2,dr1dz2;
  double dr2dx2,dr2dy2,dr2dz2,dr2dx3,dr2dy3,dr2dz3;
  double dr1dxx,dr1dxy,dr1dxz,dr1dyx,dr1dyy,dr1dyz,dr1dzx,dr1dzy,dr1dzz;
  double dr2dxx,dr2dxy,dr2dxz,dr2dyx,dr2dyy,dr2dyz,dr2dzx,dr2dzy,dr2dzz;
  double dxx,dyy,dzz,dxy,dxz,dyz,dyx,dzx,dzy;
  int nbondxs,*idxbox,*ibondx,*jbondx,*kbondx;
  double *eq1bondx,*eq2bondx,*fkbondx;

  nbondxs = coords->bondxs.nbondxs;
  idxbox = coords->bondxs.idxbox;
  ibondx = coords->bondxs.ibondx;
  jbondx = coords->bondxs.jbondx;
  kbondx = coords->bondxs.kbondx;
  eq1bondx = coords->bondxs.eq1bondx;
  eq2bondx = coords->bondxs.eq2bondx;
  fkbondx  = coords->bondxs.fkbondx;
  
  for(nb=0;nb<nbondxs;nb++){

    i = ibondx[nb];
    j = jbondx[nb];
    k = kbondx[nb];
    itype = idxbox[nb];

    ax = coords->px[i] - coords->px[j];
    ay = coords->py[i] - coords->py[j];
    az = coords->pz[i] - coords->pz[j];
    al = sqrt(ax*ax+ay*ay+az*az);
    ral = 1./al;
    ra3 = ral*ral*ral;
  
    bx = coords->px[j] - coords->px[k];
    by = coords->py[j] - coords->py[k];
    bz = coords->pz[j] - coords->pz[k];
    bl = sqrt(bx*bx+by*by+bz*bz);
    rbl = 1./bl;
    rb3 = rbl*rbl*rbl;
    
    dr1 = al - eq1bondx[itype];
    dr2 = bl - eq2bondx[itype];
    fk = fkbondx[itype];

    dr1dx1 =  ax/al; dr1dy1 =  ay/al; dr1dz1 =  az/al;
    dr1dx2 = -ax/al; dr1dy2 = -ay/al; dr1dz2 = -az/al;
    dr2dx2 =  bx/bl; dr2dy2 =  by/bl; dr2dz2 =  bz/bl;
    dr2dx3 = -bx/bl; dr2dy3 = -by/bl; dr2dz3 = -bz/bl;

    dr1dxx = -ax*ax*ra3+ral; dr1dxy = -ax*ay*ra3;     dr1dxz = -ax*az*ra3;
    dr1dyx = -ay*ax*ra3;     dr1dyy = -ay*ay*ra3+ral; dr1dyz = -ay*az*ra3;
    dr1dzx = -az*ax*ra3;     dr1dzy = -az*ay*ra3;     dr1dzz = -az*az*ra3+ral;

    dr2dxx = -bx*bx*rb3+rbl; dr2dxy = -bx*by*rb3;     dr2dxz = -bx*bz*rb3;
    dr2dyx = -by*bx*rb3;     dr2dyy = -by*by*rb3+rbl; dr2dyz = -by*bz*rb3;
    dr2dzx = -bz*bx*rb3;     dr2dzy = -bz*by*rb3;     dr2dzz = -bz*bz*rb3+rbl;

    /* d f[x,y]h[x,y]/dxy=df/dx dh/dy+dh/dx df/dy+h df2/dxy+f dh2/dxy*/
  
    /*  calculate force constant matrix interact partical 1-1  */

    dxx = fk*(dr2*dr1dxx); dxy = fk*(dr2*dr1dxy);  dxz = fk*(dr2*dr1dxz);
    dyx = fk*(dr2*dr1dyx); dyy = fk*(dr2*dr1dyy);  dyz = fk*(dr2*dr1dyz);
    dzx = fk*(dr2*dr1dzx); dzy = fk*(dr2*dr1dzy);  dzz = fk*(dr2*dr1dzz);
    
    ii = 3*i-1;
    
    fcmat[ii+1][ii+1] += dxx;
    fcmat[ii+1][ii+2] += dxy;
    fcmat[ii+1][ii+3] += dxz;
    fcmat[ii+2][ii+1] += dyx;
    fcmat[ii+2][ii+2] += dyy;
    fcmat[ii+2][ii+3] += dyz;
    fcmat[ii+3][ii+1] += dzx;
    fcmat[ii+3][ii+2] += dzy;
    fcmat[ii+3][ii+3] += dzz;
    
    /*  calculate force constant matrix for interact 1-2 */    
    dxx = fk*(dr1dx1*dr2dx2-dr2*dr1dxx);
    dxy = fk*(dr1dx1*dr2dy2-dr2*dr1dxy);
    dxz = fk*(dr1dx1*dr2dz2-dr2*dr1dxz);
    
    dyx = fk*(dr1dy1*dr2dx2-dr2*dr1dyx);
    dyy = fk*(dr1dy1*dr2dy2-dr2*dr1dyy);
    dyz = fk*(dr1dy1*dr2dz2-dr2*dr1dyz);

    dzx = fk*(dr1dz1*dr2dx2-dr2*dr1dzx);
    dzy = fk*(dr1dz1*dr2dy2-dr2*dr1dzy);
    dzz = fk*(dr1dz1*dr2dz2-dr2*dr1dzz);
    
    ii = 3*i-1;
    jj = 3*j-1;
    
    fcmat[ii+1][jj+1] += dxx;
    fcmat[ii+1][jj+2] += dxy;
    fcmat[ii+1][jj+3] += dxz;
    fcmat[ii+2][jj+1] += dyx;
    fcmat[ii+2][jj+2] += dyy;
    fcmat[ii+2][jj+3] += dyz;
    fcmat[ii+3][jj+1] += dzx;
    fcmat[ii+3][jj+2] += dzy;
    fcmat[ii+3][jj+3] += dzz;

    fcmat[jj+1][ii+1] += dxx;
    fcmat[jj+2][ii+1] += dxy;
    fcmat[jj+3][ii+1] += dxz;
    fcmat[jj+1][ii+2] += dyx;
    fcmat[jj+2][ii+2] += dyy;
    fcmat[jj+3][ii+2] += dyz;
    fcmat[jj+1][ii+3] += dzx;
    fcmat[jj+2][ii+3] += dzy;
    fcmat[jj+3][ii+3] += dzz;

    /*  calculate force constant matrix for interact 1-3 */    
    dxx=fk*(dr1dx1*dr2dx3); dxy=fk*(dr1dx1*dr2dy3); dxz=fk*(dr1dx1*dr2dz3);
    dyx=fk*(dr1dy1*dr2dx3); dyy=fk*(dr1dy1*dr2dy3); dyz=fk*(dr1dy1*dr2dz3);
    dzx=fk*(dr1dz1*dr2dx3); dzy=fk*(dr1dz1*dr2dy3); dzz=fk*(dr1dz1*dr2dz3);
    
    ii = 3*i-1;
    jj = 3*k-1;
    
    fcmat[ii+1][jj+1] += dxx;
    fcmat[ii+1][jj+2] += dxy;
    fcmat[ii+1][jj+3] += dxz;
    fcmat[ii+2][jj+1] += dyx;
    fcmat[ii+2][jj+2] += dyy;
    fcmat[ii+2][jj+3] += dyz;
    fcmat[ii+3][jj+1] += dzx;
    fcmat[ii+3][jj+2] += dzy;
    fcmat[ii+3][jj+3] += dzz;

    fcmat[jj+1][ii+1] += dxx;
    fcmat[jj+2][ii+1] += dxy;
    fcmat[jj+3][ii+1] += dxz;
    fcmat[jj+1][ii+2] += dyx;
    fcmat[jj+2][ii+2] += dyy;
    fcmat[jj+3][ii+2] += dyz;
    fcmat[jj+1][ii+3] += dzx;
    fcmat[jj+2][ii+3] += dzy;
    fcmat[jj+3][ii+3] += dzz;
    
   /*  calculate force constant matrix interact partical 2-2  */

    dxx = fk*(dr1dx2*dr2dx2+dr1dx2*dr2dx2+dr2*dr1dxx+dr1*dr2dxx);
    dxy = fk*(dr1dx2*dr2dy2+dr1dy2*dr2dx2+dr2*dr1dxy+dr1*dr2dxy);
    dxz = fk*(dr1dx2*dr2dz2+dr1dz2*dr2dx2+dr2*dr1dxz+dr1*dr2dxz);

    dyx = fk*(dr1dy2*dr2dx2+dr1dx2*dr2dy2+dr2*dr1dyx+dr1*dr2dyx);
    dyy = fk*(dr1dy2*dr2dy2+dr1dy2*dr2dy2+dr2*dr1dyy+dr1*dr2dyy);
    dyz = fk*(dr1dy2*dr2dz2+dr1dz2*dr2dy2+dr2*dr1dyz+dr1*dr2dyz);

    dzx = fk*(dr1dz2*dr2dx2+dr1dx2*dr2dz2+dr2*dr1dzx+dr1*dr2dzx);
    dzy = fk*(dr1dz2*dr2dy2+dr1dy2*dr2dz2+dr2*dr1dzy+dr1*dr2dzy);
    dzz = fk*(dr1dz2*dr2dz2+dr1dz2*dr2dz2+dr2*dr1dzz+dr1*dr2dzz);
    
    ii = 3*j-1;
    
    fcmat[ii+1][ii+1] += dxx;
    fcmat[ii+1][ii+2] += dxy;
    fcmat[ii+1][ii+3] += dxz;
    fcmat[ii+2][ii+1] += dyx;
    fcmat[ii+2][ii+2] += dyy;
    fcmat[ii+2][ii+3] += dyz;
    fcmat[ii+3][ii+1] += dzx;
    fcmat[ii+3][ii+2] += dzy;
    fcmat[ii+3][ii+3] += dzz;

    /*  calculate force constant matrix for interact 2-3 */    

    dxx = fk*(dr1dx2*dr2dx3-dr1*dr2dxx);
    dxy = fk*(dr1dx2*dr2dy3-dr1*dr2dxy);
    dxz = fk*(dr1dx2*dr2dz3-dr1*dr2dxz);

    dyx = fk*(dr1dy2*dr2dx3-dr1*dr2dyx);
    dyy = fk*(dr1dy2*dr2dy3-dr1*dr2dyy);
    dyz = fk*(dr1dy2*dr2dz3-dr1*dr2dyz);

    dzx = fk*(dr1dz2*dr2dx3-dr1*dr2dzx);
    dzy = fk*(dr1dz2*dr2dy3-dr1*dr2dzy);
    dzz = fk*(dr1dz2*dr2dz3-dr1*dr2dzz);
   
    ii = 3*j-1;
    jj = 3*k-1;
    
    fcmat[ii+1][jj+1] += dxx;
    fcmat[ii+1][jj+2] += dxy;
    fcmat[ii+1][jj+3] += dxz;
    fcmat[ii+2][jj+1] += dyx;
    fcmat[ii+2][jj+2] += dyy;
    fcmat[ii+2][jj+3] += dyz;
    fcmat[ii+3][jj+1] += dzx;
    fcmat[ii+3][jj+2] += dzy;
    fcmat[ii+3][jj+3] += dzz;

    fcmat[jj+1][ii+1] += dxx;
    fcmat[jj+2][ii+1] += dxy;
    fcmat[jj+3][ii+1] += dxz;
    fcmat[jj+1][ii+2] += dyx;
    fcmat[jj+2][ii+2] += dyy;
    fcmat[jj+3][ii+2] += dyz;
    fcmat[jj+1][ii+3] += dzx;
    fcmat[jj+2][ii+3] += dzy;
    fcmat[jj+3][ii+3] += dzz;

     /*  calculate force constant matrix interact partical 3-3  */
    dxx = fk*(dr1*dr2dxx); dxy = fk*(dr1*dr2dxy);  dxz = fk*(dr1*dr2dxz);
    dyx = fk*(dr1*dr2dyx); dyy = fk*(dr1*dr2dyy);  dyz = fk*(dr1*dr2dyz);
    dzx = fk*(dr1*dr2dzx); dzy = fk*(dr1*dr2dzy);  dzz = fk*(dr1*dr2dzz);
    
    ii = 3*k-1;
    
    fcmat[ii+1][ii+1] += dxx;
    fcmat[ii+1][ii+2] += dxy;
    fcmat[ii+1][ii+3] += dxz;
    fcmat[ii+2][ii+1] += dyx;
    fcmat[ii+2][ii+2] += dyy;
    fcmat[ii+2][ii+3] += dyz;
    fcmat[ii+3][ii+1] += dzx;
    fcmat[ii+3][ii+2] += dzy;
    fcmat[ii+3][ii+3] += dzz;
  }
}
/*------------------------------------------------------------------------*/
void fcbondxn(COORDS *coords,double **fcmat)
{
  int nb,itype;
  for(nb=0;nb<coords->bondxs.nbondxs;nb++){
    
    ibonx = coords->bondxs.ibondx[nb];
    jbonx = coords->bondxs.jbondx[nb];
    kbonx = coords->bondxs.kbondx[nb];
    itype = coords->bondxs.idxbox[nb];

    eq1bonx = coords->bondxs.eq1bondx[itype];
    eq2bonx = coords->bondxs.eq2bondx[itype];
    fkbonx = coords->bondxs.fkbondx[itype];

    ngradv2(coords->px,coords->py,coords->pz,ibonx,ibonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ibonx,jbonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ibonx,kbonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jbonx,ibonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jbonx,jbonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jbonx,kbonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,kbonx,ibonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,kbonx,jbonx,potbondx,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,kbonx,kbonx,potbondx,fcmat);
  }
}
/*---------------------------------------------------------------------*/
double potbondx(double *px,double *py,double *pz)
{
  double ax,ay,az,al,bx,by,bz,bl,dr1,dr2,spot;

  ax = px[ibonx] - px[jbonx];
  ay = py[ibonx] - py[jbonx];
  az = pz[ibonx] - pz[jbonx];
  al = sqrt(ax*ax+ay*ay+az*az);

  bx = px[jbonx] - px[kbonx];
  by = py[jbonx] - py[kbonx];
  bz = pz[jbonx] - pz[kbonx];
  bl = sqrt(bx*bx+by*by+bz*bz);
  
  dr1 = al-eq1bonx;
  dr2 = bl-eq2bonx;
  spot = fkbonx*dr1*dr2;

  return spot;
}
/*---------------------------------------------------------------------*/
