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

/* #define PRINT_ELEC */

static double spothe,spote,sprse,sprset[9];

/*-------------------------------------------------------------------*/
void zero_pote(void)
{
  int i;
  for(i=0;i<9;i++) sprset[i] = 0.0;
  spothe = spote = sprse = 0.;
}
/*-------------------------------------------------------------------*/

void force_e(int natoms,int nlen,int *indx,int *jndx,COORDS *coords,
	     INTER *inter,int iperd,int ivol,int iflg,int ipolar)
{
  int j,k,ni,nj,ierror,jer=0,ntable,natom2;
  double *x,*y,*z,*qch,*scr,*rec,*hmat,*qchij,*alpha,*alphi;
  double f,fm1,f0,f1,p,swv,dis2,fcc,qij,tablast,sphet,spet; 
  double *xdis,*ydis,*zdis,*fxds,*fyds,*fzds,*fx,*fy,*fz,hmati[9];
  double *sw,*dtab,*tab,rmin2,rmax2,rmins2,dx2;
  char line[120];
#ifdef CUBIC
  double f2;
#endif

  natom2 = natoms*2;
  ntable = inter->ntable;
  /* when iperd=3 tab & dtab are erfc fxs for ewald sums */
  dtab   = inter->dcoultab; tab = inter->coultab;
  dx2    = inter->dx2coul;  sw     = inter->switchse;
  rmins2 = inter->rmins_coul2;
  x = coords->px;  y = coords->py;  z = coords->pz;
  qch    = coords->qch;
  alpha  = coords->alpha;
  scr    = coords->scr_buf;  rec = coords->scr_rec;
  hmat   = coords->hmat;

  tablast = tab[ntable-4];
  sphet = spet = 0.;
  rmin2 = 0.;
  rmax2 = 0.;
  
  switch(iflg){
     case 0:
       rmax2 = inter->rmaxs_coul2;
       rmin2 = inter->rminl_coul2;
     case 1:
       rmax2 = rmin2 = inter->rmaxl_coul2;  
       break;
     default:
       sprintf(line,"in electonic inter-molecular force routine?!?!?"); 
       md_error(line);
       break;
  }
  
  ierror = 0;
  gethinv9(hmat,hmati);

  xdis = (double *)cmalloc(11*nlen*sizeof(double));
  ydis = xdis +   nlen;
  zdis = xdis + 2*nlen;
  fxds = xdis + 3*nlen;
  fyds = xdis + 4*nlen;
  fzds = xdis + 5*nlen;
  qchij = xdis + 6*nlen;
  alphi = xdis + 7*nlen;
  fx = xdis + 8*nlen;
  fy = xdis + 9*nlen;
  fz = xdis + 10*nlen;

#include "vectorize.h"
  for(j=0;j<nlen;j++){
    ni = indx[j]; nj = jndx[j];
    xdis[j] = x[ni] - x[nj];
    ydis[j] = y[ni] - y[nj];
    zdis[j] = z[ni] - z[nj];
    qchij[j] = qch[ni]*qch[nj];
    alphi[j] = alpha[ni];
  }

  period(nlen,xdis,ydis,zdis,hmat,hmati,iperd,ivol);
  
  if(ipolar == 1){
    force_dipole(nlen,xdis,ydis,zdis,qchij,alphi,fx,fy,fz);
  }

#include "vectorize.h"
  for(j=0;j<nlen;j++){
    dis2=xdis[j]*xdis[j]+ydis[j]*ydis[j]+zdis[j]*zdis[j];
    
    if(dis2<rmax2){
      p = (dis2-rmins2)*dx2;
      k = (int)p;  
      if(k<1){ierror = 1;jer = j;k=1;p=1;}
      p = p-(double)k;
      
      fm1= dtab[k-1];      f0 = dtab[k];      f1 = dtab[k+1];
#ifdef CUBIC
      f2 = dtab[k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else	
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      if(dis2>rmin2){
        fm1 = sw[k-1]; 	f0  = sw[k]; 	f1  = sw[k+1];
#ifdef CUBIC
        f2  = sw[k+2];
        swv = f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
           p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
        swv = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      } else {
        swv = 1.;
      }
      
      qij = qchij[j];
      fcc = qij*f;
      
#ifdef PRINT_ELEC
      printf("Electrostatic force (%d) %d-%d = %g\n",
         j,indx[j],jndx[j],fcc/KCAL);
#endif
      if(iflg==1){
        sprse -= dis2*fcc;
        sprset[0] -= xdis[j]*xdis[j]*fcc;
        sprset[1] -= ydis[j]*xdis[j]*fcc;
        sprset[2] -= zdis[j]*xdis[j]*fcc;
        sprset[4] -= ydis[j]*ydis[j]*fcc;
        sprset[5] -= zdis[j]*ydis[j]*fcc;
        sprset[8] -= zdis[j]*zdis[j]*fcc;
        
        fx[j] = xdis[j] * fcc;
        fy[j] = ydis[j] * fcc;
        fz[j] = zdis[j] * fcc;
        
        fxds[j] = fx[j]*swv;
        fyds[j] = fy[j]*swv;
        fzds[j] = fz[j]*swv;
        
        fm1 = tab[k-1];  f0  = tab[k];  f1  = tab[k+1];
#ifdef CUBIC
        f2  = tab[k+2];
        f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
           p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
        f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
        sphet += qij*(f - tablast);
        spet  += qij*f;
#ifdef PRINT_ELEC
        printf("Electrostatic energy (%d) %d-%d = %g, distance = %g\n",
           j,indx[j],jndx[j],qij*f/KCAL,sqrt(dis2));
#endif
      } else {
        fcc *= swv;
        fx[j] = xdis[j] * fcc;
        fy[j] = ydis[j] * fcc;
        fz[j] = zdis[j] * fcc;
      }
    } else {
      fx[j] = 0.;      fy[j] = 0.;      fz[j] = 0.;
      fxds[j] = 0.;    fyds[j] = 0.;    fzds[j] = 0.;
    }
  }

  spothe += sphet;
  spote  += spet;
  
  if(iflg==1){
    sprset[3] = sprset[1];
    sprset[6] = sprset[2];
    sprset[7] = sprset[5];
  }

  /* need break points to vectorize */
  for(j=0;j<nlen;j++){
    ni=indx[j];
    nj=jndx[j];
    scr[ni       ] -= fx[j];
    scr[ni+natoms] -= fy[j];
    scr[ni+natom2] -= fz[j];
    scr[nj       ] += fx[j];
    scr[nj+natoms] += fy[j];
    scr[nj+natom2] += fz[j];
    if(iflg==1){
      rec[ni       ] -= fxds[j];
      rec[ni+natoms] -= fyds[j];
      rec[ni+natom2] -= fzds[j];
      rec[nj       ] += fxds[j];
      rec[nj+natoms] += fyds[j];
      rec[nj+natom2] += fzds[j];
    }
  }

  free(xdis);
  
  if(ierror){
    sprintf(line,"distance out of range in finter coul\n");
    sprintf(line,"%s between atoms %d and %d",line,indx[jer],jndx[jer]);
    md_warning(line);
  }
}

/* ---------------------------------------------------*/
#ifdef ANALYTIC
	dis = sqrt(dis2);
	if(simparms->iperd == 3){
	  func_kerf(dis,&f0,&f,&f1,alpha_ewald);
	} else {
	  func_coul(dis,&f0,&f,&f1);
	}
	if(dis2>rminl_coul2 && dis2 < rmaxs_coul2){
	  br = (dis-sqrt(rminl_coul2))/rheal_res; 
	  sm = (1.+ br*br*(2.*br-3.));
	  if(iflg==0){
	    swv = sm; 
	  }else{
	    swv = 1.-sm;
	  }
	} else {
	  swv = 1.;
	}
	fcc = qch[indx[j]]*qch[jndx[j]]*f*swv;
	xdis[j] *= fcc; 
	ydis[j] *= fcc; 
	zdis[j] *= fcc;
#endif
/* -------------------------------------------------------------------*/
void vinter_e(int nlen,int *inp,int *jnp,COORDS *coords,INTER *inter,
	      int iperd,double *spelec,double *spelech,double *swelec,
	      double swelectensor[9],int ivol,int iflg,int ipolar)
{
  int i,k,ii,jj,ierror,jer,ier;
  int ntable;
  double *px,*py,*pz,*qch,*scr,*rec,*hmat,*alpha,*alphi;
  double *xdis,*ydis,*zdis,dis2,p,qijd;
  double f,fm1,f0,f1,f2,*qij,hmati[9];
  double pote,pothe,prse;
  double prsetensor[9];
  double *sw,*dtab,*tab,rmin2,rmax2,rmins2,dx2;
  char line[MAXLINELEN];
  
  ntable = inter->ntable;
  dtab   = inter->dcoultab; tab = inter->coultab;
  dx2    = inter->dx2coul;  sw  = inter->switchse;
  rmins2 = inter->rmins_coul2;
  scr    = coords->scr_buf; rec = coords->scr_rec;
  px = coords->px;  py = coords->py;  pz = coords->pz;
  qch    = coords->qch;
  alpha  = coords->alpha;
  hmat   = coords->hmat;
  alpha  = coords->alpha;
  gethinv9(hmat,hmati);
  pote = pothe = prse = 0.0;
  for(i=0;i<9;i++) prsetensor[i]=0.0;
  ier = jer = ierror = 0;

  xdis = (double *)malloc(5*nlen*sizeof(double));
  ydis = xdis +   nlen;
  zdis = xdis + 2*nlen;
  alphi = xdis + 3*nlen;
  qij = xdis + 4*nlen;

  rmin2 = rmax2 = 0.;
  switch(iflg){
     case 0:
       rmax2 = inter->rmaxs_coul2;
       rmin2 = inter->rminl_coul2;
     case 1:
       rmax2 = rmin2 = inter->rmaxl_coul2;  
       break;
     default:
       sprintf(line,"in electonic inter-molecular force routine?!?!?"); 
       md_error(line);
       break;
  }
  
#include "vectorize.h"
  for(i=0;i<nlen;i++){
    ii = inp[i]; jj = jnp[i];    
    xdis[i] = px[ii] - px[jj];
    ydis[i] = py[ii] - py[jj];
    zdis[i] = pz[ii] - pz[jj];
    qij[i] = qch[ii]*qch[jj];
    alphi[i] = alpha[ii];
  }

  period(nlen,xdis,ydis,zdis,hmat,hmati,iperd,ivol);
  
  if(ipolar == 1)
    pot_dipole(nlen,xdis,ydis,zdis,qij,alphi);
#include "vectorize.h"
  for(i=0;i<nlen;i++){
    dis2=xdis[i]*xdis[i]+ydis[i]*ydis[i]+zdis[i]*zdis[i];

    if(dis2 < rmax2){
      ii = inp[i]; jj = jnp[i];
      p = (dis2-rmin2)*dx2;
      k = (int)p;
      if(k<1) {ierror=1;jer=jj,ier=ii;k=1;p=0.0;}
      p = p-(double)k;

      fm1 = tab[k-1];  f0  = tab[k];  f1  = tab[k+1];
#ifdef CUBIC
      f2  = tab[k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      qijd = qij[i];
      pothe += qijd*(f - tab[ntable-4]);
      pote  += qijd*f;

      fm1= dtab[k-1]; f0 = dtab[k];   f1 = dtab[k+1];	  
#ifdef CUBIC
      f2 = dtab[k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else	
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      qijd *= f;
      prse -= dis2*qijd;
      prsetensor[0] -= xdis[i]*xdis[i]*qijd;
      prsetensor[1] -= ydis[i]*xdis[i]*qijd;
      prsetensor[2] -= zdis[i]*xdis[i]*qijd;

      prsetensor[3] -= xdis[i]*ydis[i]*qijd;
      prsetensor[4] -= ydis[i]*ydis[i]*qijd;
      prsetensor[5] -= zdis[i]*ydis[i]*qijd;

      prsetensor[6] -= xdis[i]*zdis[i]*qijd;
      prsetensor[7] -= ydis[i]*zdis[i]*qijd;
      prsetensor[8] -= zdis[i]*zdis[i]*qijd;
    }
  }
  *spelec  += pote;
  *spelech += pothe;
  *swelec  += prse;
  for(i=0;i<9;i++) swelectensor[i] += prsetensor[i];

  if(ierror){
    f0 = px[ier] - px[jer];
    f1 = py[ier] - py[jer];
    f2 = pz[ier] - pz[jer];
    period(1,&f0,&f1,&f2,hmat,hmati,iperd,ivol);
    dis2=f0*f0+f1*f1+f2*f2;

    sprintf(line,"distance out of range getvinter coul %g\n",sqrt(dis2));
    sprintf(line,"%sbetween atoms %d and %d",line,ier,jer);

    md_warning(line);
  }

  free(xdis);
}
/* -------------------------------------------------------------------*/
#ifdef ANALYTIC
    if(dis2<inter->rmax2){
      if(iperd == 3){
	func_kerf(sqrt(dis2),&f,&f0,&f1,alpha_ewald);
      } else {
	func_coul(sqrt(dis2),&f,&f0,&f1);
      }
      qijd = qch[ii]*qch[jj];
      *spelec  += qijd*f;
      *spelech += qijd*(f - inter->coultab[ntable-4]);
      *swelec  -= dis2*qijd*f0;
      *swelecx -= xdis[i]*xdis[i]*qijd*f0;
      *swelecy -= ydis[i]*ydis[i]*qijd*f0;
      *swelecz -= zdis[i]*zdis[i]*qijd*f0;
    }
#else
#endif


/* -------------------------------------------------------------------*/
void getsvirter_e(double *velec,double *welec,double welectensor[9])
{
  int i;
#ifdef PARA
  double sten[12],rten[12];

  for(i=0;i<9;i++) sten[i] = sprset[i];
  sten[9] = spote;
  sten[10] = sprse;
  sten[11] = spothe;
  MPI_Allreduce(sten,rten,12,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=0;i<9;i++) sprset[i] = rten[i];
  spote = rten[9];
  sprse = rten[10];
  spothe = rten[11];
#endif

  *velec += spote;
  *welec += sprse; 
  for(i=0;i<9;i++) welectensor[i] += sprset[i];
}
/* -------------------------------------------------------------------*/
