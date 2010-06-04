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


/*---------------------------------------------------------------*/
/* subroutines dealing with torsed interactions
   search_bond_base - searches the data base and sets up initial vectors
   ftors - gets the force of the tors
   getvtors - gets the potential energy of the tors
   pair_tors - check to see if to atoms are involved in a tors
*/

#include "md.h"

#define PRESSURE /* calculate the torsional contribution to the pressure 
		    But this is zero if you don't have cross therms! */

/* #define DEBUG */

static int itor,jtor,ktor,ltor,idxtor;
static int *ind_pots,*ind_pot_nums;
static double *eqts,*fkts;
static VPHI *vphis;

/*---------------------------------------------------------------*/
int tors_pair(int iatom,int jatom,TORS *tors)
{
  int i,j,itemp,nt;

  if(iatom>jatom){
    itemp = iatom; iatom = jatom; jatom = itemp;
  }

  for(nt=0;nt<tors->ntorss;nt++){
    i = tors->itors[nt]; j = tors->jtors[nt];
    if(i>j){itemp = i; i = j; j = itemp;}
    if(iatom == i && jatom == j) return 1;

    i = tors->itors[nt]; j = tors->ktors[nt];
    if(i>j){itemp = i; i = j; j = itemp;}
    if(iatom == i && jatom == j) return 1;

    i = tors->itors[nt]; j = tors->ltors[nt];
    if(i>j){itemp = i; i = j; j = itemp;}
    if(iatom == i && jatom == j) return 1;

    i = tors->jtors[nt]; j = tors->ktors[nt];
    if(i>j){itemp = i; i = j; j = itemp;}
    if(iatom == i && jatom == j) return 1;

    i = tors->jtors[nt]; j = tors->ltors[nt];
    if(i>j){itemp = i; i = j; j = itemp;}
    if(iatom == i && jatom == j) return 1;

    i = tors->ktors[nt]; j = tors->ltors[nt];
    if(i>j){itemp = i; i = j; j = itemp;}
    if(iatom == i && jatom == j) return 1;
  }
  return 0;
}
/*-----------------------------------------------------------------------*/
void exclude_tors(int *nexclude,int **exclude,TORS *tors)
{
  int i,j,nb;
  
  for(nb=0;nb<tors->ntorss;nb++){
    if(tors->itors[nb]<tors->jtors[nb]){
      i=tors->itors[nb];j=tors->jtors[nb];
    }else {j=tors->itors[nb];i=tors->jtors[nb];}
    insert(exclude[i],&nexclude[i],j);

    if(tors->itors[nb]<tors->ktors[nb]){i=tors->itors[nb];j=tors->ktors[nb];
    }else {j=tors->itors[nb];i=tors->ktors[nb];}
    insert(exclude[i],&nexclude[i],j);

    if(tors->itors[nb]<tors->ltors[nb]){i=tors->itors[nb];j=tors->ltors[nb];
    }else {j=tors->itors[nb];i=tors->ltors[nb];}
    insert(exclude[i],&nexclude[i],j);

    if(tors->jtors[nb]<tors->ktors[nb]){i=tors->jtors[nb];j=tors->ktors[nb];
    }else {j=tors->jtors[nb];i=tors->ktors[nb];}
    insert(exclude[i],&nexclude[i],j);

    if(tors->jtors[nb]<tors->ltors[nb]){i=tors->jtors[nb];j=tors->ltors[nb];
    }else {j=tors->jtors[nb];i=tors->ltors[nb];}
    insert(exclude[i],&nexclude[i],j);

    if(tors->ktors[nb]<tors->ltors[nb]){i=tors->ktors[nb];j=tors->ltors[nb];
    }else {j=tors->ktors[nb];i=tors->ltors[nb];}
    insert(exclude[i],&nexclude[i],j);
  }
}
/*-----------------------------------------------------------------------*/
void ftors(SIMPARMS *simparms,COORDS *coords)
{
  int i,j,k,l,i_tors,nt,ii;
  long ibegin,iend;
  double ax,ay,az,bx,by,bz,cx,cy,cz;
  double dx,dy,dz,dl,dl2,ex,ey,ez,el,el2;
  double phi,cphi,dphi,dphidcos,dudphi,dudcos=0.;
  double g,dre,erd,r11;
  double dfx1,dfy1,dfz1,dfx2,dfy2,dfz2,dfx3,dfy3,dfz3,dfx4,dfy4,dfz4;
  double dgx1,dgy1,dgz1,dgx2,dgy2,dgz2,dgx3,dgy3,dgz3,dgx4,dgy4,dgz4;
  double *px,*py,*pz;
  double *fx,*fy,*fz;
  double *eqtors,*fktors;
  int ntorss,*itors,*jtors,*ktors,*ltors,*idxtors,*ind_pot,*ind_pot_num;
  VPHI *vphi;
  int ncos = (MAXPOWER_TORS-1)/2;

  px = coords->px;
  py = coords->py;
  pz = coords->pz;
#ifndef TWO
  fx = coords->fxr;
  fy = coords->fyr;
  fz = coords->fzr;
#else 
  fx = coords->fxa;
  fy = coords->fya;
  fz = coords->fza;
#endif

  ntorss = coords->tors.ntorss;
  itors = coords->tors.itors;
  jtors = coords->tors.jtors;
  ktors = coords->tors.ktors;
  ltors = coords->tors.ltors;
  idxtors= coords->tors.idxtors;
  ind_pot= coords->tors.ind_pot;
  ind_pot_num= coords->tors.ind_pot_num;
  eqtors = coords->tors.eqtors;
  fktors = coords->tors.fktors;
  vphi = coords->tors.vphi;

  decomp1d((long)ntorss,simparms->size,simparms->rank,&ibegin,&iend);
  for(nt=ibegin;nt<iend;nt++){
    i = itors[nt];
    j = jtors[nt];
    k = ktors[nt];
    l = ltors[nt];

    ax = px[i] - px[j];
    ay = py[i] - py[j];
    az = pz[i] - pz[j];
    
    bx = px[j] - px[k];
    by = py[j] - py[k];
    bz = pz[j] - pz[k];

    cx = px[k] - px[l];
    cy = py[k] - py[l];
    cz = pz[k] - pz[l];

    dx = ay*bz - az*by;
    dy = az*bx - ax*bz;
    dz = ax*by - ay*bx;
    dl2= dx*dx+dy*dy+dz*dz;
    dl = sqrt(dl2);

    ex=by*cz - bz*cy;
    ey=bz*cx - bx*cz;
    ez=bx*cy - by*cx;
    el2=ex*ex+ey*ey+ez*ez;
    el=sqrt(el2);
    
    cphi=(dx*ex+dy*ey+dz*ez)/(dl*el);

    i_tors = ind_pot_num[idxtors[nt]];

    switch(ind_pot[idxtors[nt]]){
    case 0: /* Harmonic */
      if(cphi>=1.0 || cphi <= -1.0) {
	phi      = 0.;
	dphi = phi - eqtors[i_tors];
	dudphi   = 0.;
	dphidcos = 0.;
	dudcos   = 0.;       
      } else {
	phi = acos(cphi);
	dphi = phi - eqtors[i_tors];
	dudphi = fktors[i_tors]*dphi;
	dphidcos = -1./sqrt(1.-cphi*cphi);
	dudcos  = dudphi*dphidcos;
      }
      break;
    case 1: /* cosine power series */
      /* 
	 tpot = vphi[i_tors][MAXPOWER_TORS-1];
	 for(ii=MAXPOWER_TORS-2;ii>=0;ii--){
	 tpot = tpot*cphi + vphi[i_tors][ii];
	 } */
      dudcos = ((double)MAXPOWER_TORS-1.)*vphi[i_tors][MAXPOWER_TORS-1];
      for(ii=MAXPOWER_TORS-2;ii>0;ii--){
	dudcos = dudcos*cphi + vphi[i_tors][ii]*(double)ii;
      }
      break;
    case 2: /* cosine with phase series */
      if(cphi>=1.0 || cphi <= -1.0) {
	dudcos = 0.;       
      } else {
	phi = acos(cphi); /* we have a problem if sin(x) ~ 0 */
	dudphi = 0;
	for(ii=1;ii<ncos;++ii){
	  dudphi -= ii*vphis[i_tors][2*ii-1]*sin(ii*phi+vphis[i_tors][2*ii]);
	}
	dphidcos = -1./sqrt(1.-cphi*cphi); /* -1/sin(x) */
	dudcos = dudphi*dphidcos; /* we have a problem if sin(x) ~ 0 */      
      }
      break;
    default:
      md_error("unkown torsional potential in ftors");
      exit(1);
    }

    /* Derivatives to f=d.e, g=dl*el  */

    g=dl*el;
    dre=dl/el;
    erd=el/dl;
    r11=1./g;

    dfx1 = by*ez - bz*ey;
    dfy1 = bz*ex - bx*ez;
    dfz1 = bx*ey - by*ex;
    
    dfx2 = (cy*dz-cz*dy+(az+bz)*ey-(ay+by)*ez);
    dfy2 = (cz*dx-cx*dz+(ax+bx)*ez-(az+bz)*ex);
    dfz2 = (cx*dy-cy*dx+(ay+by)*ex-(ax+bx)*ey);
    
    dfx3 = (ay*ez-az*ey+(bz+cz)*dy-(by+cy)*dz);
    dfy3 = (az*ex-ax*ez+(bx+cx)*dz-(bz+cz)*dx);
    dfz3 = (ax*ey-ay*ex+(by+cy)*dx-(bx+cx)*dy);
    
    dfx4 = by*dz - bz*dy;
    dfy4 = bz*dx - bx*dz;
    dfz4 = bx*dy - by*dx;

    dgx1 = (by*dz - bz*dy)*erd;
    dgy1 = (bz*dx - bx*dz)*erd;
    dgz1 = (bx*dy - by*dx)*erd;

    dgx2 = ((dy*(az+bz)-dz*(ay+by))*erd+(cy*ez-cz*ey)*dre);
    dgy2 = ((dz*(ax+bx)-dx*(az+bz))*erd+(cz*ex-cx*ez)*dre);
    dgz2 = ((dx*(ay+by)-dy*(ax+bx))*erd+(cx*ey-cy*ex)*dre);
    
    dgx3 = ((ay*dz-az*dy)*erd+(ey*(bz+cz)-ez*(by+cy))*dre);
    dgy3 = ((az*dx-ax*dz)*erd+(ez*(bx+cx)-ex*(bz+cz))*dre);
    dgz3 = ((ax*dy-ay*dx)*erd+(ex*(by+cy)-ey*(bx+cx))*dre);
    
    dgx4 = (by*ez - bz*ey)*dre;
    dgy4 = (bz*ex - bx*ez)*dre;
    dgz4 = (bx*ey - by*ex)*dre;
    
    dudcos *= r11;

    fx[i] -= dudcos*(dfx1 - cphi*dgx1);
    fy[i] -= dudcos*(dfy1 - cphi*dgy1);
    fz[i] -= dudcos*(dfz1 - cphi*dgz1);
    fx[j] -= dudcos*(dfx2 - cphi*dgx2);
    fy[j] -= dudcos*(dfy2 - cphi*dgy2);
    fz[j] -= dudcos*(dfz2 - cphi*dgz2);
    fx[k] -= dudcos*(dfx3 - cphi*dgx3);
    fy[k] -= dudcos*(dfy3 - cphi*dgy3);
    fz[k] -= dudcos*(dfz3 - cphi*dgz3);
    fx[l] -= dudcos*(dfx4 - cphi*dgx4);
    fy[l] -= dudcos*(dfy4 - cphi*dgy4);
    fz[l] -= dudcos*(dfz4 - cphi*dgz4);
  }
}

/*---------------------------------------------------------------*/
void getvtors(SIMPARMS *simparms,COORDS *coords,
	      double *UI,double *WI,double WItensor[9])
{
  int i,j,k,l,i_tors,nt,ii;
  long ibegin,iend;
  double ax,ay,az,bx,by,bz;
  double cx,cy,cz,dx,dy,dz,ex,ey,ez,dl,dl2,el,el2,dde;
  double phi,cphi,dphi=0.;
  double spot,tpot,spotl;
#ifdef PRESSURE
  double sprstensor[9],sprstl[9];
  double g,dre,erd,r11;
  double dphidcos,dudphi,dudcos=0.;
  double dfx1,dfy1,dfz1,dfx3,dfy3,dfz3,dfx4,dfy4,dfz4;
  double dgx1,dgy1,dgz1,dgx3,dgy3,dgz3,dgx4,dgy4,dgz4;
  /* double dfx2,dfy2,dfz2,dgx2,dgy2,dgz2; */
#endif
  double hpot,ppot;
  double *px,*py,*pz;
  double *eqtors,*fktors;
  int ntorss,*itors,*jtors,*ktors,*ltors,*idxtors,*ind_pot,*ind_pot_num;
  VPHI *vphi;
  int ncos =  (MAXPOWER_TORS-1)/2;

  px = coords->px;
  py = coords->py;
  pz = coords->pz;
  ntorss = coords->tors.ntorss;
  itors = coords->tors.itors;
  jtors = coords->tors.jtors;
  ktors = coords->tors.ktors;
  ltors = coords->tors.ltors;
  idxtors= coords->tors.idxtors;
  ind_pot= coords->tors.ind_pot;
  ind_pot_num= coords->tors.ind_pot_num;
  eqtors = coords->tors.eqtors;
  fktors = coords->tors.fktors;
  vphi = coords->tors.vphi;

  spot = hpot = ppot = 0.;
  for(i=0;i<9;i++) sprstensor[i]=0.0;

  decomp1d((long)ntorss,simparms->size,simparms->rank,&ibegin,&iend);
  for(nt=ibegin;nt<iend;nt++){
    i = itors[nt];
    j = jtors[nt];
    k = ktors[nt];
    l = ltors[nt];

    ax = px[i] - px[j];
    ay = py[i] - py[j];
    az = pz[i] - pz[j];
    
    bx = px[j] - px[k];
    by = py[j] - py[k];
    bz = pz[j] - pz[k];

    cx = px[k] - px[l];
    cy = py[k] - py[l];
    cz = pz[k] - pz[l];

    dx = ay*bz - az*by;
    dy = az*bx - ax*bz;
    dz = ax*by - ay*bx;
    dl2= dx*dx+dy*dy+dz*dz;
    dl = sqrt(dl2);

    ex=by*cz - bz*cy;
    ey=bz*cx - bx*cz;
    ez=bx*cy - by*cx;
    el2=ex*ex+ey*ey+ez*ez;
    el=sqrt(el2);
    
    dde=dx*ex+dy*ey+dz*ez;

    cphi=dde/(dl*el); 

    i_tors = ind_pot_num[idxtors[nt]];
    switch(ind_pot[idxtors[nt]]){
    case 0: /* harmonic */
      if(cphi>=1. || cphi <= -1.){
	phi=0.;
	dphi = phi - eqtors[i_tors];
#ifdef PRESSURE
	dudphi = 0.;
	dphidcos = 0.;
	dudcos  = 0.;       
#endif
      } else {
	phi = acos(cphi);
	dphi = phi - eqtors[i_tors];
#ifdef PRESSURE
	dudphi = fktors[i_tors]*dphi;
	dphidcos = -1./sqrt(1.-cphi*cphi);
	dudcos  = dudphi*dphidcos;
#endif
      }
      tpot = .5*fktors[i_tors]*dphi*dphi;
      hpot += tpot;
      spot += tpot;
      break;
    case 1: /* cosine series */
      tpot = vphi[i_tors][MAXPOWER_TORS-1];
      for(ii=MAXPOWER_TORS-2;ii>=0;ii--){
	tpot = tpot*cphi + vphi[i_tors][ii];
      }
      ppot += tpot;
      spot += tpot;
#ifdef PRESSURE
      dudcos = ((double)MAXPOWER_TORS-1.)*vphi[i_tors][MAXPOWER_TORS-1];
      for(ii=MAXPOWER_TORS-2;ii>0;ii--){
	dudcos = dudcos*cphi + vphi[i_tors][ii]*(double)ii;
      }
#endif
      break;
    case 2: /* cosine w/ phase */
      if(cphi>=1. || cphi <= -1.){
	phi=0.;
      } else {
	phi = acos(cphi);
      }
      tpot = vphis[i_tors][0];
      phi = acos(cphi);
      for(ii=1;ii<ncos;++ii){
	tpot += vphis[i_tors][2*ii-1]*cos(ii*phi+vphis[i_tors][2*ii]);
      }
      ppot += tpot;
      spot += tpot;
#ifdef PRESSURE /* if phi ~ 0 we have a problem */
      dudphi = 0;
      for(ii=1;ii<ncos;++ii){      
	dudphi -= ii*vphis[i_tors][2*ii-1]*sin(ii*phi+vphis[i_tors][2*ii]);
      }
      dphidcos = -1./sqrt(1.-cphi*cphi); /* -1/sin(x) */
      dudcos = dudphi*dphidcos;
#endif
    default:
      md_error("unknown torsional potential in getvtors");
      exit(1);
    }
#ifdef PRESSURE
    /* Derivatives to f=d.e, g=dl*el  */

    g=dl*el;
    dre=dl/el;
    erd=el/dl;
    r11=1./g;

    dfx1 = by*ez - bz*ey;
    dfy1 = bz*ex - bx*ez;
    dfz1 = bx*ey - by*ex;
    
    /* not needed since x2,y2,z2 is taken as the origin 
    dfx2 = (cy*dz-cz*dy+(az+bz)*ey-(ay+by)*ez);
    dfy2 = (cz*dx-cx*dz+(ax+bx)*ez-(az+bz)*ex);
    dfz2 = (cx*dy-cy*dx+(ay+by)*ex-(ax+bx)*ey);
    */

    dfx3 = (ay*ez-az*ey+(bz+cz)*dy-(by+cy)*dz);
    dfy3 = (az*ex-ax*ez+(bx+cx)*dz-(bz+cz)*dx);
    dfz3 = (ax*ey-ay*ex+(by+cy)*dx-(bx+cx)*dy);
    
    dfx4 = by*dz - bz*dy;
    dfy4 = bz*dx - bx*dz;
    dfz4 = bx*dy - by*dx;
    
    dgx1 = (by*dz - bz*dy)*erd;
    dgy1 = (bz*dx - bx*dz)*erd;
    dgz1 = (bx*dy - by*dx)*erd;

    /* not needed  since x2,y2,z2 is taken as the origin
    dgx2 = ((dy*(az+bz)-dz*(ay+by))*erd+(cy*ez-cz*ey)*dre);
    dgy2 = ((dz*(ax+bx)-dx*(az+bz))*erd+(cz*ex-cx*ez)*dre);
    dgz2 = ((dx*(ay+by)-dy*(ax+bx))*erd+(cx*ey-cy*ex)*dre);
    */

    dgx3 = ((ay*dz-az*dy)*erd+(ey*(bz+cz)-ez*(by+cy))*dre);
    dgy3 = ((az*dx-ax*dz)*erd+(ez*(bx+cx)-ex*(bz+cz))*dre);
    dgz3 = ((ax*dy-ay*dx)*erd+(ex*(by+cy)-ey*(bx+cx))*dre);
    
    dgx4 = (by*ez - bz*ey)*dre;
    dgy4 = (bz*ex - bx*ez)*dre;
    dgz4 = (bx*ey - by*ex)*dre;
    
    dudcos *= r11;

    /* tensor for i,j */
    sprstensor[0] -= ax*(dudcos*(dfx1 - cphi*dgx1));
    sprstensor[1] -= ax*(dudcos*(dfy1 - cphi*dgy1));
    sprstensor[2] -= ax*(dudcos*(dfz1 - cphi*dgz1));

    sprstensor[3] -= ay*(dudcos*(dfx1 - cphi*dgx1));
    sprstensor[4] -= ay*(dudcos*(dfy1 - cphi*dgy1));
    sprstensor[5] -= ay*(dudcos*(dfz1 - cphi*dgz1));

    sprstensor[6] -= az*(dudcos*(dfx1 - cphi*dgx1));
    sprstensor[7] -= az*(dudcos*(dfy1 - cphi*dgy1));
    sprstensor[8] -= az*(dudcos*(dfz1 - cphi*dgz1));

    /* tensor for j,k */
    sprstensor[0] += bx*(dudcos*(dfx3 - cphi*dgx3));
    sprstensor[1] += bx*(dudcos*(dfy3 - cphi*dgy3));
    sprstensor[2] += bx*(dudcos*(dfz3 - cphi*dgz3));

    sprstensor[3] += by*(dudcos*(dfx3 - cphi*dgx3));
    sprstensor[4] += by*(dudcos*(dfy3 - cphi*dgy3));
    sprstensor[5] += by*(dudcos*(dfz3 - cphi*dgz3));

    sprstensor[6] += bz*(dudcos*(dfx3 - cphi*dgx3));
    sprstensor[7] += bz*(dudcos*(dfy3 - cphi*dgy3));
    sprstensor[8] += bz*(dudcos*(dfz3 - cphi*dgz3));

    /* need to calculate the vector between j and l*/
    sprstensor[0] += (bx+cx)*dudcos*(dfx4 - cphi*dgx4);
    sprstensor[1] += (bx+cx)*dudcos*(dfy4 - cphi*dgy4);
    sprstensor[2] += (bx+cx)*dudcos*(dfz4 - cphi*dgz4);

    sprstensor[3] += (by+cy)*dudcos*(dfx4 - cphi*dgx4);
    sprstensor[4] += (by+cy)*dudcos*(dfy4 - cphi*dgy4);
    sprstensor[5] += (by+cy)*dudcos*(dfz4 - cphi*dgz4);

    sprstensor[6] += (bz+cz)*dudcos*(dfx4 - cphi*dgx4);
    sprstensor[7] += (bz+cz)*dudcos*(dfy4 - cphi*dgy4);
    sprstensor[8] += (bz+cz)*dudcos*(dfz4 - cphi*dgz4);

#endif
  }
#ifdef PARA
  MPI_Allreduce(&spot,&spotl,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#ifdef PRESSURE
  MPI_Allreduce(sprstensor,sprstl,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
#else
  spotl = spot;
#ifdef PRESSURE
  for(i=0;i<9;i++) sprstl[i] = sprstensor[i];
#endif
#endif
  *UI = spotl;
  *WI = sprstl[0]+sprstl[4]+sprstl[8];
#ifdef PRESSURE
  for(i=0;i<9;i++) WItensor[i] = sprstl[i];
#endif
#ifdef DEBUG
  fprintf(stdout,"hpot=%g ppot=%g UI=%g WI=%g %g %g\n",
	  hpot/KCAL,ppot/KCAL,spotl,sprstl[0],sprstl[4],sprstl[8]);
#endif
}
/*-----------------------------------------------------------------------*/
void fctors(COORDS *coords,double **fcmat)
{
  int i,j,k,l,ii,jj,nt,i_tors;
  double ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez;
  double dl,dl2,el,el2,r11,r22,r33,rdl2,rel2,dre3,erd3;
  double cphi,cpr11,phi,g,dre,erd,dde;
  double dphidcos,dphi,dudphi,dudcos=0.;
  double dfx1,dfy1,dfz1,dfx2,dfy2,dfz2,dfx3,dfy3,dfz3,dfx4,dfy4,dfz4;
  double dgx1,dgy1,dgz1,dgx2,dgy2,dgz2,dgx3,dgy3,dgz3,dgx4,dgy4,dgz4;
  double dhx1,dhy1,dhz1,dhx2,dhy2,dhz2,dhx3,dhy3,dhz3,dhx4,dhy4,dhz4;
  double dcdx1,dcdy1,dcdz1,dcdx2,dcdy2,dcdz2;
  double dcdx3,dcdy3,dcdz3,dcdx4,dcdy4,dcdz4;

  double du2dphi2,dphi2cos2,du2dcos2=0.;
  double dxx,dyy,dzz,dxy,dxz,dyz,dyx,dzx,dzy;
  double dc2dxx,dc2dyy,dc2dzz,dc2dxy,dc2dxz,dc2dyx,dc2dyz,dc2dzx,dc2dzy;
  
  double abx,aby,abz,bcx,bcy,bcz;

  double *px,*py,*pz;
  double *eqtors,*fktors;
  int ntorss,*itors,*jtors,*ktors,*ltors,*idxtors,*ind_pot,*ind_pot_num;
  VPHI *vphi;
  int ncos =  (MAXPOWER_TORS-1)/2;

  px = coords->px;
  py = coords->py;
  pz = coords->pz;
  ntorss = coords->tors.ntorss;
  itors = coords->tors.itors;
  jtors = coords->tors.jtors;
  ktors = coords->tors.ktors;
  ltors = coords->tors.ltors;
  idxtors= coords->tors.idxtors;
  ind_pot= coords->tors.ind_pot;
  ind_pot_num= coords->tors.ind_pot_num;
  eqtors = coords->tors.eqtors;
  fktors = coords->tors.fktors;
  vphi = coords->tors.vphi;

  for(nt=0;nt<ntorss;nt++){
    i = itors[nt];
    j = jtors[nt];
    k = ktors[nt];
    l = ltors[nt];

    ax = px[i] - px[j];
    ay = py[i] - py[j];
    az = pz[i] - pz[j];
    
    bx = px[j] - px[k];
    by = py[j] - py[k];
    bz = pz[j] - pz[k];

    cx = px[k] - px[l];
    cy = py[k] - py[l];
    cz = pz[k] - pz[l];

    dx=ay*bz-az*by;
    dy=az*bx-ax*bz;
    dz=ax*by-ay*bx;
    dl2=dx*dx+dy*dy+dz*dz;
    dl=sqrt(dl2);

    ex=by*cz-bz*cy;
    ey=bz*cx-bx*cz;
    ez=bx*cy-by*cx;
    el2=ex*ex+ey*ey+ez*ez;
    el=sqrt(el2);

    g=dl*el;
    dre=dl/el;
    erd=el/dl;
    r11=1./g;
    r22=r11*r11;
    r33=r22*r11;
    rdl2=1./dl2;
    rel2=1./el2;
    dre3=dre*rel2;
    erd3=erd*rdl2;
      
    dde=dx*ex+dy*ey+dz*ez;
    cphi=dde*r11;
    cpr11 = cphi*r11;

    i_tors = ind_pot_num[idxtors[nt]];

    switch(ind_pot[idxtors[nt]]){
    case 0:
      if(cphi>=1. || cphi <= -1.){
	phi=0.;
	dphidcos = 0.;
	dphi2cos2 = 0.;
      } else {
	  phi = acos(cphi);
	  dphidcos = -1./sqrt(1.-cphi*cphi);
	  dphi2cos2= -cphi/pow(1-cphi*cphi,3./2.);
      }
      dphi = phi - eqtors[i_tors];
      dudphi = fktors[i_tors]*dphi;
      du2dphi2 = fktors[i_tors];
      dudcos  = dudphi*dphidcos;
      du2dcos2= (du2dphi2*dphidcos*dphidcos+dudphi*dphi2cos2);
      break;
    case 1:
      dudcos = ((double)(MAXPOWER_TORS-1))*vphi[i_tors][MAXPOWER_TORS-1];
      for(ii=MAXPOWER_TORS-2;ii>0;ii--){
        dudcos = dudcos*cphi + vphi[i_tors][ii]*(double)ii;
      }
      du2dcos2 = ((double)((MAXPOWER_TORS-1)*(MAXPOWER_TORS-2)))*
	vphi[i_tors][MAXPOWER_TORS-1];
      for(ii=MAXPOWER_TORS-3;ii>0;ii--){
	du2dcos2 = du2dcos2*cphi + vphi[i_tors][ii+1]*(double)(ii*(ii+1));
      }
      break;
    case 2:
      if(cphi>=1. || cphi <= -1.){
	phi=0.;
	dphidcos = 0.;
	dphi2cos2 = 0.;
      } else {
	phi = acos(cphi);
	dphidcos = -1./sqrt(1.-cphi*cphi);
	dphi2cos2= -cphi/pow(1-cphi*cphi,3./2.);
      }
      dudphi = 0;
      for(ii=1;ii<ncos;++ii){      
	dudphi -= ii*vphis[i_tors][2*ii-1]*sin(ii*phi+vphis[i_tors][2*ii]);
      }
      du2dphi2 = 0;
      for(ii=1;ii<ncos;++ii){      
	du2dphi2 -= ii*ii*vphis[i_tors][2*ii-1]*cos(ii*phi+vphis[i_tors][2*ii]);
      }
      dudcos  = dudphi*dphidcos;
      du2dcos2= (du2dphi2*dphidcos*dphidcos+dudphi*dphi2cos2);
      break;
    default:
      md_error("unknown torsional potential in fctors");
      exit(1);
    }

    /* usefull constants */
    abx = ax+bx;
    aby = ay+by;
    abz = az+bz;
    bcx = bx+cx;
    bcy = by+cy;
    bcz = bz+cz;

    /* Derivatives to f=dde, g=dl*el  */ 
    
    dfx1 = by*ez-bz*ey;
    dfy1 = bz*ex-bx*ez;
    dfz1 = bx*ey-by*ex;
    dfx2 = cy*dz-cz*dy + abz*ey-aby*ez;
    dfy2 = cz*dx-cx*dz + abx*ez-abz*ex;
    dfz2 = cx*dy-cy*dx + aby*ex-abx*ey;
    dfx3 = ay*ez-az*ey + bcz*dy-bcy*dz;
    dfy3 = az*ex-ax*ez + bcx*dz-bcz*dx;
    dfz3 = ax*ey-ay*ex + bcy*dx-bcx*dy;
    dfx4 = by*dz-bz*dy;
    dfy4 = bz*dx-bx*dz;
    dfz4 = bx*dy-by*dx;
    
    dhx4 = dy*abz-dz*aby;
    dhy4 = dz*abx-dx*abz;
    dhz4 = dx*aby-dy*abx;
    dhx3 = cy*ez-cz*ey;
    dhy3 = cz*ex-cx*ez;
    dhz3 = cx*ey-cy*ex;
    dhx2 = ay*dz-az*dy;
    dhy2 = az*dx-ax*dz;
    dhz2 = ax*dy-ay*dx;
    dhx1 = ey*bcz-ez*bcy;
    dhy1 = ez*bcx-ex*bcz;
    dhz1 = ex*bcy-ey*bcx;

    dgx1 = dfx4*erd;
    dgy1 = dfy4*erd;
    dgz1 = dfz4*erd;
    dgx2 = dhx4*erd + dhx3*dre;
    dgy2 = dhy4*erd + dhy3*dre;
    dgz2 = dhz4*erd + dhz3*dre;

    dgx3 = dhx2*erd + dhx1*dre;
    dgy3 = dhy2*erd + dhy1*dre;
    dgz3 = dhz2*erd + dhz1*dre;
    dgx4 = dfx1*dre;
    dgy4 = dfy1*dre;
    dgz4 = dfz1*dre;
    
    /*  get the first derivatives to the cphi  */
    
    dcdx1 = r11*(dfx1 - cphi*dgx1);
    dcdy1 = r11*(dfy1 - cphi*dgy1);
    dcdz1 = r11*(dfz1 - cphi*dgz1);
    dcdx2 = r11*(dfx2 - cphi*dgx2);
    dcdy2 = r11*(dfy2 - cphi*dgy2);
    dcdz2 = r11*(dfz2 - cphi*dgz2);
    dcdx3 = r11*(dfx3 - cphi*dgx3);
    dcdy3 = r11*(dfy3 - cphi*dgy3);
    dcdz3 = r11*(dfz3 - cphi*dgz3);
    dcdx4 = r11*(dfx4 - cphi*dgx4);
    dcdy4 = r11*(dfy4 - cphi*dgy4);
    dcdz4 = r11*(dfz4 - cphi*dgz4);
    
    /*  calculate force constant matrix interact partical 1-1  */
    
    dc2dxx = (2.*dde*dgx1*dgx1-2.*g*dfx1*dgx1)*r33;
    dc2dxx += -((by*by+bz*bz)-dfx4*dfx4*rdl2)*erd*cpr11;
    dc2dyy = (2.*dde*dgy1*dgy1-2.*g*dfy1*dgy1)*r33;
    dc2dyy += -((bz*bz+bx*bx)-dfy4*dfy4*rdl2)*erd*cpr11;
    dc2dzz = (2.*dde*dgz1*dgz1-2.*g*dfz1*dgz1)*r33;
    dc2dzz += -((bx*bx+by*by)-dfz4*dfz4*rdl2)*erd*cpr11;
    
    dc2dxy = (2.*cphi*dgx1*dgy1-dgx1*dfy1-dfx1*dgy1)*r22;
    dc2dxy += ((bx*by)+dfy4*dfx4*rdl2)*erd*cpr11;
    dc2dxz = (2.*cphi*dgx1*dgz1-dgx1*dfz1-dfx1*dgz1)*r22;
    dc2dxz += -(-(bx*bz)-dfz4*dfx4*rdl2)*erd*cpr11;
    dc2dyz = (2.*cphi*dgy1*dgz1-dgy1*dfz1-dfy1*dgz1)*r22;
    dc2dyz += -(-(by*bz)-dfz4*dfy4*rdl2)*erd*cpr11;
    
    dc2dyx = dc2dxy;
    dc2dzx = dc2dxz;
    dc2dzy = dc2dyz;
    
    dxx = du2dcos2*dcdx1*dcdx1 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx1*dcdy1 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx1*dcdz1 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy1*dcdx1 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy1*dcdy1 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy1*dcdz1 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz1*dcdx1 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz1*dcdy1 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz1*dcdz1 + dudcos*dc2dzz;

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
    
    dc2dxx = (2.*cphi*dgx1*dgx2-dfx1*dgx2-dfx2*dgx1)*r22;
    dc2dxx += ((by*cy+bz*cz)*r11 -
	       (dfx4*dhx3*r11-(by*aby+bz*abz)*erd-dfx4*dhx4*erd3)*cpr11);
    dc2dxy = (2.*cphi*dgx1*dgy2-dfx1*dgy2-dfy2*dgx1)*r22;
    dc2dxy += ((bx*cy-2.*by*cx)*r11 -
	       (dfx4*dhy3*r11+(dz+by*abx)*erd-dfx4*dhy4*erd3)*cpr11);
    dc2dxz = (2.*cphi*dgx1*dgz2-dfx1*dgz2-dfz2*dgx1)*r22;
    dc2dxz += ((-ey-bz*cx)*r11 - 
	       (dfx4*dhz3*r11-(dy-bz*abx)*erd-dfx4*dhz4*erd3)*cpr11);
    dc2dyx = (2.*cphi*dgy1*dgx2-dfy1*dgx2-dfx2*dgy1)*r22;
    dc2dyx += ((-ez-bx*cy)*r11 -
	       (dfy4*dhx3*r11-(dz-bx*aby)*erd-dfy4*dhx4*erd3)*cpr11);
    dc2dyy = (2.*cphi*dgy1*dgy2-dfy1*dgy2-dfy2*dgy1)*r22;
    dc2dyy += ((bx*cx+bz*cz)*r11 -
	       (dfy4*dhy3*r11-(bz*abz+bx*abx)*erd-dfy4*dhy4*erd3)*cpr11);
    dc2dyz = (2.*cphi*dgy1*dgz2-dfy1*dgz2-dfz2*dgy1)*r22;
    dc2dyz += ((by*cz-2.*bz*cy)*r11 - 
	       (dfy4*dhz3*r11+(dx+bz*aby)*erd-dfy4*dhz4*erd3)*cpr11);
    dc2dzx = (2.*cphi*dgz1*dgx2-dfz1*dgx2-dfx2*dgz1)*r22;
    dc2dzx += ((bz*cx-2.*bx*cz)*r11 -
	       (dfz4*dhx3*r11+(dy+bx*abz)*erd-dfz4*dhx4*erd3)*cpr11);
    dc2dzy = (2.*cphi*dgz1*dgy2-dfz1*dgy2-dfy2*dgz1)*r22;
    dc2dzy += ((-ex-by*cz)*r11 -
	       (dfz4*dhy3*r11-(dx-by*abz)*erd-dfz4*dhy4*erd3)*cpr11);
    dc2dzz = (2.*cphi*dgz1*dgz2-dfz1*dgz2-dfz2*dgz1)*r22;
    dc2dzz += ((bx*cx+by*cy)*r11 -
	       (dfz4*dhz3*r11-(by*aby+bx*abx)*erd-dfz4*dhz4*erd3)*cpr11);
    
    dxx = du2dcos2*dcdx1*dcdx2 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx1*dcdy2 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx1*dcdz2 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy1*dcdx2 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy1*dcdy2 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy1*dcdz2 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz1*dcdx2 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz1*dcdy2 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz1*dcdz2 + dudcos*dc2dzz;

    ii=3*i-1;
    jj=3*j-1;
    
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

    dc2dxx = (2.*cphi*dgx1*dgx3-dfx1*dgx3-dfx3*dgx1)*r22;
    dc2dxx += (-(by*bcy+bz*bcz)*r11 -
	       (dfx4*dhx1*r11+(ay*by+az*bz)*erd-dhx2*dfx4*erd3)*cpr11);
    dc2dxy = (2.*cphi*dgx1*dgy3-dfx1*dgy3-dfy3*dgx1)*r22;
    dc2dxy += ((by*bcx-ez)*r11 -
	       (dfx4*dhy1*r11-(dz+ax*by)*erd-dhy2*dfx4*erd3)*cpr11);
    dc2dxz = (2.*cphi*dgx1*dgz3-dfx1*dgz3-dfz3*dgx1)*r22;
    dc2dxz += ((ey + bz*bcx)*r11 -
	       (dfx4*dhz1*r11+(dy-ax*bz)*erd-dhz2*dfx4*erd3)*cpr11);
    dc2dyx = (2.*cphi*dgy1*dgx3-dfy1*dgx3-dfx3*dgy1)*r22;
    dc2dyx += ((ez+bx*bcy)*r11 -
	       (dfy4*dhx1*r11+(dz-ay*bx)*erd-dhx2*dfy4*erd3)*cpr11);
    dc2dyy = (2.*cphi*dgy1*dgy3-dfy1*dgy3-dfy3*dgy1)*r22;
    dc2dyy += (-(bz*bcz+bx*bcx)*r11 -
	       (dfy4*dhy1*r11+(az*bz+ax*bx)*erd-dhy2*dfy4*erd3)*cpr11);
    dc2dyz = (2.*cphi*dgy1*dgz3-dfy1*dgz3-dfz3*dgy1)*r22;
    dc2dyz += ((bz*bcy-ex)*r11 -
	       (dfy4*dhz1*r11-(dx+ay*bz)*erd-dhz2*dfy4*erd3)*cpr11);
    dc2dzx = (2.*cphi*dgz1*dgx3-dfz1*dgx3-dfx3*dgz1)*r22;
    dc2dzx += ((bx*bcz-ey)*r11 -
	       (dfz4*dhx1*r11-(dy+az*bx)*erd-dhx2*dfz4*erd3)*cpr11);
    dc2dzy = (2.*cphi*dgz1*dgy3-dfz1*dgy3-dfy3*dgz1)*r22;
    dc2dzy += ((by*bcz+ex)*r11 -
	       (dfz4*dhy1*r11+(dx-az*by)*erd-dhy2*dfz4*erd3)*cpr11);
    dc2dzz = (2.*cphi*dgz1*dgz3-dfz1*dgz3-dfz3*dgz1)*r22;
    dc2dzz += (-(bx*bcx+by*bcy)*r11-
	       (dfz4*dhz1*r11+(ax*bx+ay*by)*erd-dhz2*dfz4*erd3)*cpr11);
    

    dxx = du2dcos2*dcdx1*dcdx3 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx1*dcdy3 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx1*dcdz3 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy1*dcdx3 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy1*dcdy3 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy1*dcdz3 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz1*dcdx3 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz1*dcdy3 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz1*dcdz3 + dudcos*dc2dzz;

    ii=3*i-1;
    jj=3*k-1;
    
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

    /*  calculate force constant matrix for interact 1-4 */
    
    dc2dxx = (2.*cphi*dgx1*dgx4-dfx1*dgx4-dfx4*dgx1)*r22;
    dc2dxx += ((by*by+bz*bz)*r11-(dfx4*dfx1*r11)*cpr11);
    dc2dxy = (2.*cphi*dgx1*dgy4-dfx1*dgy4-dfy4*dgx1)*r22;
    dc2dxy += ((-bx*by)*r11-(dfx4*dfy1*r11)*cpr11);
    dc2dxz = (2.*cphi*dgx1*dgz4-dfx1*dgz4-dfz4*dgx1)*r22;
    dc2dxz += ((-bx*bz)*r11-(dfx4*dfz1*r11)*cpr11);
    dc2dyx = (2.*cphi*dgy1*dgx4-dfy1*dgx4-dfx4*dgy1)*r22;
    dc2dyx += ((-bx*by)*r11-(dfy4*dfx1*r11)*cpr11);
    dc2dyy = (2.*cphi*dgy1*dgy4-dfy1*dgy4-dfy4*dgy1)*r22;
    dc2dyy += ((bx*bx+bz*bz)*r11-(dfy4*dfy1*r11)*cpr11);
    dc2dyz = (2.*cphi*dgy1*dgz4-dfy1*dgz4-dfz4*dgy1)*r22;
    dc2dyz += ((-by*bz)*r11-(dfy4*dfz1*r11)*cpr11);
    dc2dzx = (2.*cphi*dgz1*dgx4-dfz1*dgx4-dfx4*dgz1)*r22;
    dc2dzx += ((-bx*bz)*r11-(dfz4*dfx1*r11)*cpr11);
    dc2dzy = (2.*cphi*dgz1*dgy4-dfz1*dgy4-dfy4*dgz1)*r22;
    dc2dzy += ((-bz*by)*r11-(dfz4*dfy1*r11)*cpr11);
    dc2dzz = (2.*cphi*dgz1*dgz4-dfz1*dgz4-dfz4*dgz1)*r22;
    dc2dzz += ((bx*bx+by*by)*r11-(dfz4*dfz1*r11)*cpr11);
    
    dxx = du2dcos2*dcdx1*dcdx4 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx1*dcdy4 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx1*dcdz4 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy1*dcdx4 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy1*dcdy4 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy1*dcdz4 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz1*dcdx4 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz1*dcdy4 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz1*dcdz4 + dudcos*dc2dzz;

    ii=3*i-1;
    jj=3*l-1;
    
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

    dc2dxx = (2.*dde*dgx2*dgx2-2.*g*dfx2*dgx2)*r33;
    dc2dxx += ( -2.*(cy*aby+cz*abz)*r11 -
	       ((cy*cy+cz*cz)*dre+(abz*abz+aby*aby)*erd +
		dhx3*dhx4*2.*r11-dhx3*dhx3*dre3-dhx4*dhx4*erd3)*cpr11);
    dc2dyy = (2.*dde*dgy2*dgy2-2.*g*dfy2*dgy2)*r33;
    dc2dyy += ( -2.*(cz*abz+cx*abx)*r11 -
	       ((cx*cx+cz*cz)*dre+(abx*abx+abz*abz)*erd +
		dhy3*dhy4*2.*r11-dhy3*dhy3*dre3-dhy4*dhy4*erd3)*cpr11);
    dc2dzz = (2.*dde*dgz2*dgz2-2.*g*dfz2*dgz2)*r33;
    dc2dzz += ( -2.*(cx*abx+cy*aby)*r11 -
	       ((cx*cx+cy*cy)*dre+(abx*abx+aby*aby)*erd +
		dhz3*dhz4*2.*r11-dhz3*dhz3*dre3-dhz4*dhz4*erd3)*cpr11);
    dc2dxy = (2.*cphi*dgx2*dgy2-dgx2*dfy2-dfx2*dgy2)*r22;
    dc2dxy +=((cy*abx+cx*aby)*r11 - 
	      ( -cx*cy*dre-abx*aby*erd +
	       (dhx3*dhy4+dhy3*dhx4)*r11-dhy4*dhx4*erd3-dhy3*dhx3*dre3)*cpr11);
    dc2dxz = (2.*cphi*dgx2*dgz2-dgx2*dfz2-dfx2*dgz2)*r22;
    dc2dxz +=((cx*abz+cz*abx)*r11 -
	      ( -cx*cz*dre-abz*abx*erd +
	       (dhx3*dhz4+dhz3*dhx4)*r11-dhz4*dhx4*erd3-dhz3*dhx3*dre3)*cpr11);
    dc2dyz = (2.*cphi*dgy2*dgz2-dgy2*dfz2-dfy2*dgz2)*r22;
    dc2dyz +=((cz*aby+cy*abz)*r11 -
	      ( -cy*cz*dre-aby*abz*erd +
	       (dhz3*dhy4+dhy3*dhz4)*r11-dhy4*dhz4*erd3-dhz3*dhy3*dre3)*cpr11);

    dc2dyx = dc2dxy;
    dc2dzx = dc2dxz;
    dc2dzy = dc2dyz;

    dxx = du2dcos2*dcdx2*dcdx2 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx2*dcdy2 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx2*dcdz2 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy2*dcdx2 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy2*dcdy2 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy2*dcdz2 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz2*dcdx2 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz2*dcdy2 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz2*dcdz2 + dudcos*dc2dzz;

    ii=3*j-1;
    
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

    dc2dxx = (2.*cphi*dgx2*dgx3-dfx2*dgx3-dfx3*dgx2)*r22;
    dc2dxx +=((ay*cy+az*cz+abz*bcz+aby*bcy)*r11 -
	      (-(ay*aby+az*abz)*erd-(cz*bcz+cy*bcy)*dre+
	       (dhx2*dhx3+dhx4*dhx1)*r11-dhx2*dhx4*erd3-dhx3*dhx1*dre3)*cpr11);
    dc2dxy = (2.*cphi*dgx2*dgy3-dfx2*dgy3-dfy3*dgx2)*r22;
    dc2dxy +=((dz+ez-ax*cy-aby*bcx)*r11 -
	      ((dz+ax*aby)*erd+(ez+cy*bcx)*dre+		
	       (dhy2*dhx3+dhx4*dhy1)*r11-dhy2*dhx4*erd3-dhx3*dhy1*dre3)*cpr11);
    dc2dxz = (2.*cphi*dgx2*dgz3-dfx2*dgz3-dfz3*dgx2)*r22;
    dc2dxz +=((-dy-ey-ax*cz-bcx*abz)*r11 -
	      ((ax*abz-dy)*erd+(cz*bcx-ey)*dre+
	       (dhz2*dhx3+dhx4*dhz1)*r11-dhz2*dhx4*erd3-dhx3*dhz1*dre3)*cpr11);
    dc2dyx = (2.*cphi*dgy2*dgx3-dfy2*dgx3-dfx3*dgy2)*r22;
    dc2dyx +=((-dz-ez-ay*cx-bcy*abx)*r11 -
	      ((ay*abx-dz)*erd+(cx*bcy-ez)*dre+
	       (dhx2*dhy3+dhy4*dhx1)*r11-dhx2*dhy4*erd3-dhy3*dhx1*dre3)*cpr11);
    dc2dyy = (2.*cphi*dgy2*dgy3-dfy2*dgy3-dfy3*dgy2)*r22;
    dc2dyy +=((az*cz+ax*cx+abz*bcz+abx*bcx)*r11 -
	      (-(az*abz+ax*abx)*erd-(cz*bcz+cx*bcx)*dre+
	       (dhy2*dhy3+dhy4*dhy1)*r11-dhy2*dhy4*erd3-dhy3*dhy1*dre3)*cpr11);
    dc2dyz = (2.*cphi*dgy2*dgz3-dfy2*dgz3-dfz3*dgy2)*r22;
    dc2dyz +=((dx+ex-ay*cz-bcy*abz)*r11 -
	      ((ay*abz+dx)*erd+(cz*bcy+ex)*dre +
	       (dhz2*dhy3+dhy4*dhz1)*r11-dhz2*dhy4*erd3-dhy3*dhz1*dre3)*cpr11);
    dc2dzx = (2.*cphi*dgz2*dgx3-dfz2*dgx3-dfx3*dgz2)*r22;
    dc2dzx +=((dy+ey-az*cx-bcz*abx)*r11 -
	      ((az*abx+dy)*erd+(cx*bcz+ey)*dre +
	       (dhx2*dhz3+dhz4*dhx1)*r11-dhx2*dhz4*erd3-dhz3*dhx1*dre3)*cpr11);
    dc2dzy = (2.*cphi*dgz2*dgy3-dfz2*dgy3-dfy3*dgz2)*r22;
    dc2dzy +=((-dx-ex-az*cy-bcz*aby)*r11 -
	      ((az*aby-dx)*erd+(cy*bcz-ex)*dre +
	       (dhy2*dhz3+dhz4*dhy1)*r11-dhy2*dhz4*erd3-dhz3*dhy1*dre3)*cpr11);
    dc2dzz = (2.*cphi*dgz2*dgz3-dfz2*dgz3-dfz3*dgz2)*r22;
    dc2dzz +=((ax*cx+ay*cy+abx*bcx+aby*bcy)*r11 -
	      (-(ax*abx+ay*aby)*erd-(cx*bcx+cy*bcy)*dre+
	       (dhz2*dhz3+dhz4*dhz1)*r11-dhz2*dhz4*erd3-dhz3*dhz1*dre3)*cpr11);
    
    dxx = du2dcos2*dcdx2*dcdx3 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx2*dcdy3 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx2*dcdz3 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy2*dcdx3 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy2*dcdy3 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy2*dcdz3 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz2*dcdx3 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz2*dcdy3 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz2*dcdz3 + dudcos*dc2dzz;

    ii=3*j-1;
    jj=3*k-1;
    
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

    /*  calculate force constant matrix for interact 2-4 */
    
    dc2dxx = (2.*cphi*dgx2*dgx4-dfx2*dgx4-dfx4*dgx2)*r22;
    dc2dxx += (-(by*aby+bz*abz)*r11-
	       ((by*cy+bz*cz)*dre+dfx1*dhx4*r11-dfx1*dhx3*dre3)*cpr11);
    dc2dxy = (2.*cphi*dgx2*dgy4-dfx2*dgy4-dfy4*dgx2)*r22;
    dc2dxy += ((bx*aby-dz)*r11 -
	       (-(ez+bx*cy)*dre+dfy1*dhx4*r11-dfy1*dhx3*dre3)*cpr11);
    dc2dxz = (2.*cphi*dgx2*dgz4-dfx2*dgz4-dfz4*dgx2)*r22;
    dc2dxz += ((bx*abz+dy)*r11 -
	       ((ey-bx*cz)*dre+dfz1*dhx4*r11-dfz1*dhx3*dre3)*cpr11);
    dc2dyx = (2.*cphi*dgy2*dgx4-dfy2*dgx4-dfx4*dgy2)*r22;
    dc2dyx += ((dz+by*abx)*r11 -
	       ((ez-by*cx)*dre+dfx1*dhy4*r11-dfx1*dhy3*dre3)*cpr11);
    dc2dyy = (2.*cphi*dgy2*dgy4-dfy2*dgy4-dfy4*dgy2)*r22;
    dc2dyy += (-(bz*abz+bx*abx)*r11 -
	       ((bz*cz+bx*cx)*dre+dfy1*dhy4*r11-dfy1*dhy3*dre3)*cpr11);
    dc2dyz = (2.*cphi*dgy2*dgz4-dfy2*dgz4-dfz4*dgy2)*r22;
    dc2dyz += ((by*abz-dx)*r11 -
	       (-(ex+by*cz)*dre+dfz1*dhy4*r11 -
		dfz1*dhy3*dre3)*cpr11);
    dc2dzx = (2.*cphi*dgz2*dgx4-dfz2*dgx4-dfx4*dgz2)*r22;
    dc2dzx += ((bz*abx-dy)*r11 -
	       (-(ey+bz*cx)*dre+dfx1*dhz4*r11-dfx1*dhz3*dre3)*cpr11);
    dc2dzy = (2.*cphi*dgz2*dgy4-dfz2*dgy4-dfy4*dgz2)*r22;
    dc2dzy += ((bz*aby+dx)*r11 -
	       ((ex-bz*cy)*dre+dfy1*dhz4*r11-dfy1*dhz3*dre3)*cpr11);
    dc2dzz = (2.*cphi*dgz2*dgz4-dfz2*dgz4-dfz4*dgz2)*r22;
    dc2dzz += (-(bx*abx+by*aby)*r11 -
	       ((bx*cx+by*cy)*dre+dfz1*dhz4*r11-dfz1*dhz3*dre3)*cpr11);
    

    dxx = du2dcos2*dcdx2*dcdx4 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx2*dcdy4 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx2*dcdz4 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy2*dcdx4 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy2*dcdy4 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy2*dcdz4 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz2*dcdx4 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz2*dcdy4 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz2*dcdz4 + dudcos*dc2dzz;

    ii=3*j-1;
    jj=3*l-1;
    
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

    dc2dxx = (2.*dde*dgx3*dgx3-2.*g*dfx3*dgx3)*r33;
    dc2dxx += ( -2.*(ay*bcy+az*bcz)*r11 -
	       ((ay*ay+az*az)*erd+(bcy*bcy+bcz*bcz)*dre +
		dhx2*dhx1*2.*r11-dhx1*dhx1*dre3-dhx2*dhx2*erd3)*cpr11);
    dc2dyy = (2.*dde*dgy3*dgy3-2.*g*dfy3*dgy3)*r33;
    dc2dyy += ( -2.*(az*bcz+ax*bcx)*r11 -
	       ((az*az+ax*ax)*erd+(bcz*bcz+bcx*bcx)*dre +
		dhy2*dhy1*2.*r11-dhy1*dhy1*dre3-dhy2*dhy2*erd3)*cpr11);
    dc2dzz = (2.*dde*dgz3*dgz3-2.*g*dfz3*dgz3)*r33;
    dc2dzz += ( -2.*(ax*bcx+ay*bcy)*r11 -
	       ((ax*ax+ay*ay)*erd+(bcx*bcx+bcy*bcy)*dre +
		dhz2*dhz1*2.*r11-dhz1*dhz1*dre3-dhz2*dhz2*erd3)*cpr11);

    dc2dxy = (2.*cphi*dgx3*dgy3-dgx3*dfy3-dfx3*dgy3)*r22;
    dc2dxy += ((ay*bcx+ax*bcy)*r11 - 
	       (-ax*ay*erd-bcx*bcy*dre+(dhx2*dhy1+dhy2*dhx1)*r11-
		dhy1*dhx1*dre3-dhy2*dhx2*erd3)*cpr11);
    
    dc2dxz = (2.*cphi*dgx3*dgz3-dgx3*dfz3-dfx3*dgz3)*r22;
    dc2dxz += ((ax*bcz+az*bcx)*r11 - 
	       ( -ax*az*erd-bcx*bcz*dre+(dhx2*dhz1+dhz2*dhx1)*r11 -
		dhz1*dhx1*dre3-dhz2*dhx2*erd3)*cpr11);
    dc2dyz = (2.*cphi*dgy3*dgz3-dgy3*dfz3-dfy3*dgz3)*r22;
    dc2dyz += ((ay*bcz+az*bcy)*r11 - 
	       ( -ay*az*erd-bcy*bcz*dre+(dhz2*dhy1+dhy2*dhz1)*r11 -
		dhy1*dhz1*dre3-dhz2*dhy2*erd3)*cpr11);
    
    dc2dyx = dc2dxy;
    dc2dzx = dc2dxz;
    dc2dzy = dc2dyz;
    
    dxx = du2dcos2*dcdx3*dcdx3 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx3*dcdy3 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx3*dcdz3 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy3*dcdx3 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy3*dcdy3 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy3*dcdz3 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz3*dcdx3 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz3*dcdy3 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz3*dcdz3 + dudcos*dc2dzz;

    ii=3*k-1;
    
    fcmat[ii+1][ii+1] += dxx;
    fcmat[ii+1][ii+2] += dxy;
    fcmat[ii+1][ii+3] += dxz;
    fcmat[ii+2][ii+1] += dyx;
    fcmat[ii+2][ii+2] += dyy;
    fcmat[ii+2][ii+3] += dyz;
    fcmat[ii+3][ii+1] += dzx;
    fcmat[ii+3][ii+2] += dzy;
    fcmat[ii+3][ii+3] += dzz;

    /*  calculate force constant matrix for interact 3-4 */

    dc2dxx = (2.*cphi*dgx3*dgx4-dfx3*dgx4-dfx4*dgx3)*r22;
    dc2dxx += ((ay*by+az*bz)*r11 -
	       (-(by*bcy+bz*bcz)*dre+dhx2*dfx1*r11-dfx1*dhx1*dre3)*cpr11);
    dc2dxy = (2.*cphi*dgx3*dgy4-dfx3*dgy4-dfy4*dgx3)*r22;
    dc2dxy += ((ax*by-2.*ay*bx)*r11 -
	       ((ez+bx*bcy)*dre+dhx2*dfy1*r11-dfy1*dhx1*dre3)*cpr11);
    dc2dxz = (2.*cphi*dgx3*dgz4-dfx3*dgz4-dfz4*dgx3)*r22;
    dc2dxz += (-(dy+az*bx)*r11 -
	       ((bx*bcz-ey)*dre+dhx2*dfz1*r11-dfz1*dhx1*dre3)*cpr11);
    dc2dyx = (2.*cphi*dgy3*dgx4-dfy3*dgx4-dfx4*dgy3)*r22;
    dc2dyx += (-(dz+ax*by)*r11 -
	       ((by*bcx-ez)*dre+dhy2*dfx1*r11-dfx1*dhy1*dre3)*cpr11);
    dc2dyy = (2.*cphi*dgy3*dgy4-dfy3*dgy4-dfy4*dgy3)*r22;
    dc2dyy += ((az*bz+ax*bx)*r11 -
	       (-(bz*bcz+bx*bcx)*dre+dhy2*dfy1*r11-dfy1*dhy1*dre3)*cpr11);
    dc2dyz = (2.*cphi*dgy3*dgz4-dfy3*dgz4-dfz4*dgy3)*r22;
    dc2dyz += ((ay*bz-2.*az*by)*r11 -
	       ((by*bcz+ex)*dre+dhy2*dfz1*r11-dfz1*dhy1*dre3)*cpr11);
    dc2dzx = (2.*cphi*dgz3*dgx4-dfz3*dgx4-dfx4*dgz3)*r22;
    dc2dzx += ((az*bx-2.*ax*bz)*r11 -
	       ((bz*bcx+ey)*dre+dhz2*dfx1*r11-dfx1*dhz1*dre3)*cpr11);
    dc2dzy = (2.*cphi*dgz3*dgy4-dfz3*dgy4-dfy4*dgz3)*r22;
    dc2dzy += (-(ay*bz+dx)*r11 -
	       ((bz*bcy-ex)*dre+dhz2*dfy1*r11-dfy1*dhz1*dre3)*cpr11);
    
    dc2dzz = (2.*cphi*dgz3*dgz4-dfz3*dgz4-dfz4*dgz3)*r22;
    dc2dzz += ((ax*bx+ay*by)*r11 -
	       (-(bx*bcx+by*bcy)*dre+dhz2*dfz1*r11-dfz1*dhz1*dre3)*cpr11);
    
    dxx = du2dcos2*dcdx3*dcdx4 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx3*dcdy4 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx3*dcdz4 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy3*dcdx4 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy3*dcdy4 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy3*dcdz4 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz3*dcdx4 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz3*dcdy4 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz3*dcdz4 + dudcos*dc2dzz;

    ii=3*k-1;
    jj=3*l-1;
    
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

    /*  calculate force constant matrix interact partical 4-4  */

    dc2dxx = (2.*dde*dgx4*dgx4-2.*g*dfx4*dgx4)*r33;
    dc2dxx += -((by*by+bz*bz)*dre-dfx1*dfx1*dre3)*cpr11;
    dc2dyy = (2.*dde*dgy4*dgy4-2.*g*dfy4*dgy4)*r33;
    dc2dyy += -((bz*bz+bx*bx)*dre-dfy1*dfy1*dre3)*cpr11;
    dc2dzz = (2.*dde*dgz4*dgz4-2.*g*dfz4*dgz4)*r33;
    dc2dzz += -((bx*bx+by*by)*dre-dfz1*dfz1*dre3)*cpr11;
    
    dc2dxy = (2.*cphi*dgx4*dgy4-dgx4*dfy4-dfx4*dgy4)*r22;
    dc2dxy += (bx*by*dre+dfy1*dfx1*dre3)*cpr11;
    dc2dxz = (2.*cphi*dgx4*dgz4-dgx4*dfz4-dfx4*dgz4)*r22;
    dc2dxz += (bx*bz*dre+dfz1*dfx1*dre3)*cpr11;
    dc2dyz = (2.*cphi*dgy4*dgz4-dgy4*dfz4-dfy4*dgz4)*r22;
    dc2dyz += (by*bz*dre+dfz1*dfy1*dre3)*cpr11;
    
    dc2dyx = dc2dxy;
    dc2dzx = dc2dxz;
    dc2dzy = dc2dyz;
    
    dxx = du2dcos2*dcdx4*dcdx4 + dudcos*dc2dxx;
    dxy = du2dcos2*dcdx4*dcdy4 + dudcos*dc2dxy;
    dxz = du2dcos2*dcdx4*dcdz4 + dudcos*dc2dxz;
    dyx = du2dcos2*dcdy4*dcdx4 + dudcos*dc2dyx;
    dyy = du2dcos2*dcdy4*dcdy4 + dudcos*dc2dyy;
    dyz = du2dcos2*dcdy4*dcdz4 + dudcos*dc2dyz;
    dzx = du2dcos2*dcdz4*dcdx4 + dudcos*dc2dzx;
    dzy = du2dcos2*dcdz4*dcdy4 + dudcos*dc2dzy;
    dzz = du2dcos2*dcdz4*dcdz4 + dudcos*dc2dzz;

    ii=3*l-1;
    
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

/*----------------------------------------------------------------*/
void fctorsn(COORDS *coords,double **fcmat)
{
  int nt;

  vphis = coords->tors.vphi;
  fkts  = coords->tors.fktors;
  eqts  = coords->tors.eqtors;
  ind_pots = coords->tors.ind_pot;
  ind_pot_nums = coords->tors.ind_pot_num;

  for(nt=0;nt<coords->tors.ntorss;nt++){
    
    itor = coords->tors.itors[nt]; jtor = coords->tors.jtors[nt];  
    ktor = coords->tors.ktors[nt]; ltor = coords->tors.ltors[nt];
    idxtor = coords->tors.idxtors[nt];

    ngradv2(coords->px,coords->py,coords->pz,itor,itor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,itor,jtor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,itor,ktor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,itor,ltor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jtor,itor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jtor,jtor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jtor,ktor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jtor,ltor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ktor,itor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ktor,jtor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ktor,ktor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ktor,ltor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ltor,itor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ltor,jtor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ltor,ktor,pottors,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,ltor,ltor,pottors,fcmat);
  }
}
/*---------------------------------------------------------------------*/
double pottors(double *px,double *py,double *pz)
{
  int i_tors,ii,ncos;
  double ax,ay,az,bx,by,bz;
  double cx,cy,cz,dx,dy,dz,ex,ey,ez,dl,dl2,el,el2,dde;
  double phi,cphi;
  double spot,tpot;

  ncos = (MAXPOWER_TORS-1)/2;
  spot = 0.;

  ax = px[itor] - px[jtor];
  ay = py[itor] - py[jtor];
  az = pz[itor] - pz[jtor];

  bx = px[ktor] - px[jtor];
  by = py[ktor] - py[jtor];
  bz = pz[ktor] - pz[jtor];
  
  cx = px[ktor] - px[ltor];
  cy = py[ktor] - py[ltor];
  cz = pz[ktor] - pz[ltor];

  dx = ay*bz - az*by;
  dy = az*bx - ax*bz;
  dz = ax*by - ay*bx;
  dl2= dx*dx+dy*dy+dz*dz;
  dl = sqrt(dl2);

  ex=by*cz - bz*cy;
  ey=bz*cx - bx*cz;
  ez=bx*cy - by*cx;
  el2=ex*ex+ey*ey+ez*ez;
  el=sqrt(el2);
    
  dde=dx*ex+dy*ey+dz*ez;
  
  cphi=dde/(dl*el); 

  i_tors = ind_pot_nums[idxtor];
  switch(ind_pots[idxtor]){
  case 0: /* harmonic */
    phi = acos(cphi);
    phi = eqts[i_tors]-phi;
    tpot = .5*fkts[i_tors]*phi*phi;
    spot += tpot;
    break;
  case 1: /* cosine series */
    tpot = vphis[i_tors][MAXPOWER_TORS-1];
    for(ii=MAXPOWER_TORS-2;ii>=0;--ii){
      tpot = tpot*cphi + vphis[i_tors][ii];
    }
    spot += tpot;
    break;
  case 2: /* cosine w/ phase */
    tpot = vphis[i_tors][0];
    for(ii=1;ii<ncos;++ii){
      tpot += vphis[i_tors][2*ii-1]*cos(ii*phi+vphis[i_tors][2*ii]);
    }
    spot += tpot;
    break;
  default:
    md_error("unknown torsional potential in pottors");
  }

  return spot;
}
/*---------------------------------------------------------------------*/