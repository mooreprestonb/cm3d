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


/* subroutines included in the charge-charge interactions 
   (ie 1/r and ewald sums) */

#include "md.h"

/* #define CUBIC */
/* #define ANALYTIC */
/* #define CHECK_BKG */
/* #define PCORR */
/* #define DEBUG */

/* #define DEBUG */
#define TPI  (2.*M_PI)
#define HELPFULL

#ifdef ANALYTIC
extern double alpha_ewald;
#endif

/*-----------------------------------------------------------------------*/
void ewald_setup(SIMPARMS *simparms,COUL *coul,double *qch,int kmax,
		 double alpha)
{
  int i,j,k;
  int ksqmax,memk,ksq;
  double sqrt_pi,self,back;
  LINE line;

  sqrt_pi = sqrt(M_PI);
  coul->kappa = alpha;
  coul->kappa2 = coul->kappa*coul->kappa;
  coul->falp2 = 1./(4.*coul->kappa2);

  self = 0.;
  back = 0.;

  for(i=0;i<coul->ncharge;i++){
    j = coul->icharge[i];
    self += qch[j]*qch[j];
    back += qch[j];
  }
  coul->ewald_self = self;
  coul->ewald_bkg  = back;

  if(coul->ncharge == 0){
    coul->ktot = 0;
  } else {
    coul->ewald_self *= -coul->kappa/sqrt_pi;
    coul->ewald_bkg = -.5*coul->ewald_bkg*coul->ewald_bkg*M_PI/(coul->kappa2);
#ifdef CHECK_BKG
    if(ewald_bkg < ERRMAX) ewald_bkg = 0.;
#endif
    memk = kmax;
    coul->ka = (int *)cmalloc(memk*sizeof(int));
    coul->kb = (int *)cmalloc(memk*sizeof(int));
    coul->kc = (int *)cmalloc(memk*sizeof(int));
    coul->kup = (int *)cmalloc((memk+1)*sizeof(int));
    
    coul->ktot = 0;
    ksqmax = kmax*kmax;
    
    if(simparms->rank==0){
      sprintf(line,"Setting up k-vectors for ewald sum (kmax = %d)",kmax);
      md_stdout(line);
    }
    for(i=0;i<=kmax;i++){
      for(j = (i==0)?0:-kmax;j<=kmax;j++){
        coul->kup[coul->ktot] = 1;
        for(k= (i==0 && j==0)?1:-kmax;k<=kmax;k++){
          ksq = i*i+j*j+k*k;
          if(ksq<=ksqmax){
            coul->ka[coul->ktot] = i;
            coul->kb[coul->ktot] = j;
            coul->kc[coul->ktot] = k;
            coul->kup[coul->ktot+1] = 0;
            coul->ktot++;
            if(memk==coul->ktot){
              memk += kmax;
              coul->ka = (int *)realloc(coul->ka,memk*sizeof(int));
              coul->kb = (int *)realloc(coul->kb,memk*sizeof(int));
              coul->kc = (int *)realloc(coul->kc,memk*sizeof(int));
              coul->kup = (int *)realloc(coul->kup,(memk+1)*sizeof(int));
            }
          }
        }
      }
    }
    
    coul->ka = (int *)realloc(coul->ka,coul->ktot*sizeof(int));
    coul->kb = (int *)realloc(coul->kb,coul->ktot*sizeof(int));
    coul->kc = (int *)realloc(coul->kc,coul->ktot*sizeof(int));
    coul->kup = (int *)realloc(coul->kup,(coul->ktot+1)*sizeof(int));
    coul->cossc = (double *)malloc(4*coul->ncharge*sizeof(double));
    coul->sinsc = coul->cossc +   coul->ncharge;
    coul->helr  = coul->cossc + 2*coul->ncharge;
    coul->heli  = coul->cossc + 3*coul->ncharge;
    simparms->mem_bytes += 3*coul->ktot*sizeof(int);
    simparms->mem_bytes += (coul->ktot+1)*sizeof(int);
    simparms->mem_bytes += 4*coul->ncharge*sizeof(double);
  }

  if(simparms->rank==0){
    sprintf(line,"There are %d charged particles (out of a total of %d)",
	    coul->ncharge,simparms->natoms);
    md_stdout(line);
    sprintf(line,"There were %d k-vectors set up.",coul->ktot);
    md_stdout(line);
    sprintf(line,"Self correction term = %g K",coul->ewald_self);
    md_stdout(line);
    sprintf(line,"Back ground term = %g K = %g e",coul->ewald_bkg,back);
    md_stdout(line);
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
}

/*-------------------------------------------------------------*/
void k_ewald(COORDS *coords,double *vkspace)
{
  int i,j,jj;
  int *ka,*kb,*kc;
  double aka,akb,akc;
  double vol,vk,qi;
  double sumr,sumi,g2,preg,smag,rk;
  double *hmati,*px,*py,*pz;

  ka = coords->coul.ka;
  kb = coords->coul.kb;
  kc = coords->coul.kc;
  hmati = coords->hmati;
  px = coords->px;
  py = coords->py;
  pz = coords->pz;
  vol = get_deth(coords->hmat);
  vk = 0.;

  for(i=0;i<coords->coul.ktot;i++){
    aka = (ka[i]*hmati[0] + kb[i]*hmati[3] + kc[i]*hmati[6])*TPI; 
    akb = (ka[i]*hmati[1] + kb[i]*hmati[4] + kc[i]*hmati[7])*TPI;
    akc = (ka[i]*hmati[2] + kb[i]*hmati[5] + kc[i]*hmati[8])*TPI;

    g2 = aka*aka+akb*akb+akc*akc;
    preg = 2.*exp(-g2*coords->coul.falp2)/(g2*vol);
    sumi = sumr = 0.;
    for(j=0;j<coords->coul.ncharge;j++){
      jj = coords->coul.icharge[j];  
      qi = coords->qch[jj];
      rk = (px[jj]*aka + py[jj]*akb + pz[jj]*akc);
      sumr += qi*cos(rk);
      sumi += qi*sin(rk);
    }
    smag = sumr*sumr+sumi*sumi;
    vk += preg*smag;
  }
  *vkspace = vk;
}
/*----------------------------------------------------------*/
void fk_ewald(SIMPARMS *simparms,COORDS *coords)
{
  long ibegin,iend;
  int i,j,jj,*icharge,ncharge;
  int *ka,*kb,*kc,*kup;
  double aka,akb,akc,rvol,vk,vp,q;
  double sumr,sumi,g2,preg,prep,smag,rk;
  double srx,sry,srz,six,siy,siz;
  double *cossc,*sinsc,*px,*py,*pz,*qch,hmati[9];
  double vk_local,vp_local,vp_local_tensor[9],vp_tensor[9];
  double *fx,*fy,*fz,*tx,*ty,*tz,*tq;
  double *helr,*heli;
  double rhr,rhi,fact;
#ifdef PARA
  int nj;
#endif
  
  fact = sqrt(4.0*M_PI); /* to convert from 1/(4 PI E_0) to 1/(2 E_0) */
  ka = coords->coul.ka;  
  kb = coords->coul.kb;  
  kc = coords->coul.kc;
  kup = coords->coul.kup;
  icharge = coords->coul.icharge;
  px = coords->px;   py = coords->py;   pz = coords->pz;  qch= coords->qch;
  fx = coords->fxl;  fy = coords->fyl;  fz = coords->fzl;

  cossc = coords->coul.cossc;
  sinsc = coords->coul.sinsc;
  helr  = coords->coul.helr;
  heli  = coords->coul.heli;

  rvol = 1./gethinv9(coords->hmat,hmati);
  vp = vk = vp_local = vk_local=0.;
  for(i=0;i<9;i++) vp_tensor[i] = vp_local_tensor[i]=0.;
  ncharge = coords->coul.ncharge;

  tx = (double *)cmalloc(4*ncharge*sizeof(double));
  ty = tx +   ncharge;
  tz = tx + 2*ncharge;
  tq = tx + 3*ncharge;

#include "vectorize.h"
  for(j=0;j<ncharge;j++){
    jj = icharge[j];
    tx[j] = px[jj];    ty[j] = py[jj];    tz[j] = pz[jj];    
    tq[j] = qch[jj]*fact;
  }

#ifdef PARA
  for(i=0;i<simparms->natoms*3;i++) coords->scr_buf[i]=coords->scr_rec[i]=0.0;
#endif

  /* get sc compenents for recursion */
#include "vectorize.h"
  for(j=0;j<ncharge;j++){
    akc = (tx[j]*hmati[6] + ty[j]*hmati[7] + tz[j]*hmati[8])*TPI;
    cossc[j] = cos(akc);
    sinsc[j] = sin(akc);
  }

  decomp1d((long)coords->coul.ktot,simparms->size,simparms->rank,
	   &ibegin,&iend);
  for(i=ibegin;i<iend;i++){
    six = TPI*ka[i];    siy = TPI*kb[i];    siz = TPI*kc[i];
    aka = (six*hmati[0] + siy*hmati[3] + siz*hmati[6]); 
    akb = (six*hmati[1] + siy*hmati[4] + siz*hmati[7]);
    akc = (six*hmati[2] + siy*hmati[5] + siz*hmati[8]);

    g2 = aka*aka+akb*akb+akc*akc;
    preg = exp(-g2*coords->coul.falp2)*rvol/g2;
    /* preg = exp(-g2*coords->coul.falp2)/g2; */
    prep = -2.*preg*(g2*coords->coul.falp2+1.)/g2;

#ifndef HELPFULL
    kup[i] = 1; test to see that recursion works
#endif
    if(kup[i] || i==ibegin){
#include "vectorize.h" 
      for(j=0;j<ncharge;j++){
        rk = (tx[j]*aka + ty[j]*akb + tz[j]*akc);
        helr[j] = cos(rk); 
        heli[j] = sin(rk);
      }
#ifdef HELPFULL
    } else {
      /* update helpfull vectors with recursion 
         ie cos(a+b) = cos(a)cos(b)-sin(a)sin(b) etc..*/
#include "vectorize.h"
      for(j=0;j<ncharge;j++){
        rhr = helr[j];   rhi = heli[j];
        sumr = cossc[j]; sumi = sinsc[j];
        helr[j] = rhr*sumr-rhi*sumi;
        heli[j] = rhi*sumr+rhr*sumi;
      }
#endif
    }
    sumr = sumi = 0.;
    for(j=0;j<ncharge;j++){
      q = tq[j]; 
      sumr += q*helr[j]; 
      sumi += q*heli[j];
    }
    smag = sumr*sumr+sumi*sumi;

    vk_local += preg*smag;
    vp_local += prep*smag*g2;
    
    vp_local_tensor[0] += prep*smag*aka*aka;
    vp_local_tensor[1] += prep*smag*akb*aka;
    vp_local_tensor[2] += prep*smag*akc*aka;

    vp_local_tensor[3] += prep*smag*aka*akb;
    vp_local_tensor[4] += prep*smag*akb*akb;
    vp_local_tensor[5] += prep*smag*akc*akb;

    vp_local_tensor[6] += prep*smag*aka*akc;
    vp_local_tensor[7] += prep*smag*akb*akc;
    vp_local_tensor[8] += prep*smag*akc*akc;
    sumr *= 2.*preg; 
    sumi *= 2.*preg;

#include "vectorize.h"
     for(j=0;j<ncharge;j++){
      jj = icharge[j];      q = tq[j];
      rhr = helr[j];      rhi = heli[j];
      srx = aka*sumr*q; sry = akb*sumr*q; srz = akc*sumr*q;
      six = aka*sumi*q; siy = akb*sumi*q; siz = akc*sumi*q;

#ifdef PARA
      nj=3*jj;
      coords->scr_buf[nj  ] += srx*rhi - six*rhr;
      coords->scr_buf[nj+1] += sry*rhi - siy*rhr;
      coords->scr_buf[nj+2] += srz*rhi - siz*rhr;
#else
      fx[jj] += srx*rhi - six*rhr;
      fy[jj] += sry*rhi - siy*rhr;
      fz[jj] += srz*rhi - siz*rhr;
#endif
    }
  }

#ifdef PARA
  MPI_Allreduce(&vk_local,&vk,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&vp_local,&vp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(vp_local_tensor,vp_tensor,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else 
  vk = vk_local;
  vp = vp_local;
  for (i=0;i<9;i++) vp_tensor[i]=vp_local_tensor[i];
#endif

  six = coords->coul.ewald_bkg/get_deth(coords->hmat);
  coords->coul.vk_ewald = vk+six+coords->coul.ewald_self;
  coords->coul.vp_ewald = vp + 3.*(vk+six);
  
  coords->coul.vp_ewald_tensor[0]=vp_tensor[0]+vk+six;
  coords->coul.vp_ewald_tensor[1]=vp_tensor[1];
  coords->coul.vp_ewald_tensor[2]=vp_tensor[2];
  coords->coul.vp_ewald_tensor[3]=vp_tensor[3];
  coords->coul.vp_ewald_tensor[4]=vp_tensor[4]+vk+six;
  coords->coul.vp_ewald_tensor[5]=vp_tensor[5];
  coords->coul.vp_ewald_tensor[6]=vp_tensor[6];
  coords->coul.vp_ewald_tensor[7]=vp_tensor[7];
  coords->coul.vp_ewald_tensor[8]=vp_tensor[8]+vk+six;
 
#ifdef PARA
  MPI_Allreduce(coords->scr_buf,coords->scr_rec,simparms->natoms*3,
		MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(j=nj=0;j<simparms->natoms;j++,nj+=3){
    fx[j]+=coords->scr_rec[nj  ];
    fy[j]+=coords->scr_rec[nj+1];
    fz[j]+=coords->scr_rec[nj+2];
  }
#endif
  if (tx != NULL) free(tx);
  
}
/*------------------------------------------------------------*/

void getvewald(SIMPARMS *simparms,COORDS *coords,
	       double *pot_ewald,double *prs_ewald,double prs_ewald_tensor[9])
{
  int i;
  *pot_ewald += coords->coul.vk_ewald+coords->coul.pot_ecorr; 
  *prs_ewald += coords->coul.vp_ewald+coords->coul.prs_ecorr;
  for(i=0;i<9;i++) {
    prs_ewald_tensor[i] += (coords->coul.vp_ewald_tensor[i]+
			    coords->coul.prs_ecorr_tensor[i]); 
  }
#ifdef PCORR
  printf("ewald vk = %lg, ecorr = %lg\n",
     coords->coul.vk_ewald,coords->coul.pot_ecorr);
  printf("ewald vp = %lg, ecorr = %lg\n",
     coords->coul.vp_ewald,coords->coul.prs_ecorr);
#endif
}
/*--------------------------------------------------------------*/
/* routine to calculate ewald corrections to pairs of charged
   particles that don't interact via a coulomb term in the first image */

void ecorr(SIMPARMS *simparms,COORDS *coords)
{
  int i,ii,jj,ncorr,*icorr,*jcorr;
  double dvecorr,sprs,spot;
  double xdis,ydis,zdis,dis2,dis;
  double erfr,qijr,ralp,kappa,kappa2,sqrt_pi;
  double *fx,*fy,*fz,*px,*py,*pz,*qch;
  double sprs_tensor[9];

  spot = sprs = 0.;
  for(i=0;i<9;i++) sprs_tensor[i]=0.0;

  sqrt_pi = 1./sqrt(M_PI);
  kappa  = coords->coul.kappa;
  kappa2 = coords->coul.kappa2;
  fx = coords->fxl;  fy = coords->fyl;  fz = coords->fzl;
  ncorr = coords->coul.necorr;
  icorr = coords->coul.iecorr;
  jcorr = coords->coul.jecorr;
  px = coords->px;  py = coords->py;  pz = coords->pz;
  qch = coords->qch;
  
  for(i=0;i<ncorr;i++){
    ii = icorr[i];     jj = jcorr[i];
    xdis = px[ii] - px[jj]; 
    ydis = py[ii] - py[jj]; 
    zdis = pz[ii] - pz[jj];
    dis2 = xdis*xdis+ydis*ydis+zdis*zdis;
    dis = sqrt(dis2);
    qijr = qch[ii]*qch[jj]/dis;

#ifdef DEBUG
    {
      LINE line;
      if(dis == 0.0){
        sprintf(line,"distance error in ecorr: %d %d %d dis=%g\n",i,ii,jj,dis);
        md_error(line);
      }
    }
#endif
      
    ralp = dis*kappa;
    erfr = erf(ralp);
    spot += -qijr*erfr;
    dvecorr = -qijr*(-2.*kappa*exp(-kappa2*dis2)*sqrt_pi+erfr/dis);
    sprs += dis*dvecorr;

    dvecorr /= dis;

    sprs_tensor[0]+= dvecorr*xdis*xdis;
    sprs_tensor[1]+= dvecorr*xdis*ydis;
    sprs_tensor[2]+= dvecorr*xdis*zdis;

    sprs_tensor[3]+= dvecorr*ydis*xdis;
    sprs_tensor[4]+= dvecorr*ydis*ydis;
    sprs_tensor[5]+= dvecorr*ydis*zdis;

    sprs_tensor[6]+= dvecorr*zdis*xdis;
    sprs_tensor[7]+= dvecorr*zdis*ydis;
    sprs_tensor[8]+= dvecorr*zdis*zdis;
    
    xdis *= dvecorr; ydis *= dvecorr; zdis *= dvecorr;
    fx[ii] += xdis; fy[ii] += ydis; fz[ii] += zdis;
    fx[jj] -= xdis; fy[jj] -= ydis; fz[jj] -= zdis;
  }
  coords->coul.pot_ecorr = spot;
  coords->coul.prs_ecorr = sprs;
  for(i=0;i<9;i++) coords->coul.prs_ecorr_tensor[i] = sprs_tensor[i];
}
/*----------------------------------------------------------------------*/
void r_ewald(SIMPARMS *simparms,COORDS *coords,double *vrspace)
{
  int i,j,ii,jj;
  double qi,xi,yi,zi,dx,dy,dz,r,vr;

  vr = 0.0;

  for(i=0;i<coords->coul.ncharge-1;i++){
    ii = coords->coul.icharge[i]; 
    xi = coords->px[ii]; 
    yi = coords->py[ii]; 
    zi = coords->pz[ii]; 
    qi = coords->qch[ii];

    for(j=i+1;j<coords->coul.ncharge;j++){
      jj = coords->coul.icharge[j];  
      dx = xi - coords->px[jj]; 
      dy = yi - coords->py[jj]; 
      dz = zi - coords->pz[jj];

      period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,
	     simparms->iperd,simparms->ivol);

      r = sqrt(dx*dx+dy*dy+dz*dz);
      vr += qi*coords->qch[jj]*erfc(coords->coul.kappa*r)/r;
    }
  }
  *vrspace = vr;
}
/*------------------------------------------------------------------------*/
void fr_ewald(SIMPARMS *simparms,COORDS *coords)
{
  int i,j,ii,jj;
  double qi,xi,yi,zi,dx,dy,dz,r,r2,sqrt_pi;
  double expkr2,erfkr,vir,fxi,fyi,fzi;
  double *px,*py,*pz,*qch,*fx,*fy,*fz;

  px = coords->px;
  py = coords->py;
  pz = coords->pz;
  fx = coords->fxl;
  fy = coords->fyl;
  fz = coords->fzl;
  qch= coords->qch;
  sqrt_pi = sqrt(M_PI);
  
  for(i=0;i<coords->coul.ncharge-1;i++){
    ii = coords->coul.icharge[i]; 
    xi = px[ii]; 
    yi = py[ii]; 
    zi = pz[ii]; 
    qi = qch[ii];
    fxi = fyi = fzi = 0.;
    for(j=i+1;j<coords->coul.ncharge;j++){
      jj = coords->coul.icharge[j];  
      dx = xi - px[jj];  
      dy = yi - py[jj]; 
      dz = zi - pz[jj];

      period(1,&dx,&dy,&dz,coords->hmat,coords->hmati,
	     simparms->iperd,simparms->ivol);

      r2 = dx*dx+dy*dy+dz*dz;
      r = sqrt(r2);
      erfkr = erfc(coords->coul.kappa*r);
      expkr2 = exp(coords->coul.kappa2*r2);
      vir = qi*qch[jj]*(erfkr/r+2.*coords->coul.kappa/(sqrt_pi*expkr2));
      dx *= vir/r2; dy *= vir/r2;  dz *= vir/r2;
      fxi -= dx;    fyi -= dy;     fzi -= dz;
      fx[jj] += dx; 
      fy[jj] += dy;  
      fz[jj] += dz;
    }
    fx[ii] += fxi;
    fy[ii] += fyi;  
    fz[ii] += fzi;
  }
}
/*------------------------------------------------------------------------*/
