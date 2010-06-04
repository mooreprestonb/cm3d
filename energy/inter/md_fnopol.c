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


/* routine to get the forces from an interaction list */

#define NOERROR
#include "md.h"

static double spotl,spothl,sprsl,sprslt[9];

/* #define VERBOSE */
/* #define PRINT_VDW_ENERGIES */

void zero_potl(void)
{
  int i;
  for(i=0;i<9;i++) sprslt[i] = 0;
  spotl = spothl = sprsl = 0.;
}

/*-----------------------------------------------------------------------*/
void force_npol(SIMPARMS *simparms,COORDS *coords,int nlen,int *indx,int *jndx,
		int nbreak,int *ibreak,INTER *inter,int iflg)
{
  int natoms,iperd,ivol,natoms2;
  int j,k,in,ni,nj,jer,ierror,ntable,ntm4;
  int *ityp,*jtyp,*itype,**imap;
  double swv,f,p,fm1,f0,f1,dis2,fcc,spotlt,spothlt,df,d2f,dfr;
  double fm1q,f0q,f1q,f2q;
  double *xdis,*ydis,*zdis,*fxds,*fyds,*fzds,**sw,*dvt,*d2vt,*vt;
  double *x,*y,*z,*hmat,*scr,*rec;
  double *rmax2,*dx2tab,**dvtab,*rminh2,*rmint2,**vtab;
  double **d2vtab;
  double hmati[9];
  LINE line;
#ifdef CUBIC
  double f2;
#endif

  x = coords->px;  y = coords->py;  z = coords->pz; hmat = coords->hmat;
  scr = coords->scr_buf; rec = coords->scr_rec;
  natoms= simparms->natoms;
  iperd = simparms->iperd;
  ivol  = simparms->ivol;
  natoms2 = 2*natoms;
  spotlt = spothlt = 0.;
  rmax2 = NULL;
  ntable= inter->ntable;
  ntm4  = ntable-4;
  itype = inter->itype;
  sw    = inter->switchs;
  imap  = inter->map;
  rminh2= inter->rminl_cut2;
  /* rmaxh2= inter->rmaxs_cut2; if you want a long range force short cutoff */
  rmint2= inter->rmins_cut2;
  dx2tab= inter->dx2tab;
  dvtab = inter->dvtab;
  vtab  = inter->vtab;
  d2vtab = inter-> d2vtab;

  /* switch to long pair or short pair */
  switch(iflg){
  case 0:
    rmax2 = inter->rmaxs_cut2;  
    break;
  case 1:
    rmax2 = inter->rmaxl_cut2; 
    break;
  default:
    md_error("ERROR: in inter-molecular force routine?!?!?");
  }

  jer = -1; swv = 1.;  ierror=0;
  gethinv9(hmat,hmati);
  ityp = (int *)cmalloc(2*nlen*sizeof(int));
  jtyp = ityp +   nlen;
  xdis = (double *)cmalloc(6*nlen*sizeof(double));
  ydis = xdis +   nlen;
  zdis = xdis + 2*nlen;
  fxds = xdis + 3*nlen;  fyds = xdis + 4*nlen;  fzds = xdis + 5*nlen;

#ifdef VERBOSE
  for(j=0;j<nlen;j++)printf("interaction %d pair %d - %d\n",j,indx[j],jndx[j]);
#endif

#include "vectorize.h"
  for(j=0;j<nlen;j++){
    ityp[j] = itype[indx[j]]; 
    jtyp[j] = itype[jndx[j]];
  }
#include "vectorize.h"
  for(j=0;j<nlen;j++){
    ni = indx[j]; nj = jndx[j];
    xdis[j] = x[ni] - x[nj];
    ydis[j] = y[ni] - y[nj];
    zdis[j] = z[ni] - z[nj];
  }
  
  period(nlen,xdis,ydis,zdis,hmat,hmati,iperd,ivol);
  
#include "vectorize.h"
  for(j=0;j<nlen;j++){
    in = imap[ityp[j]][jtyp[j]];
    dis2 = xdis[j]*xdis[j]+ydis[j]*ydis[j]+zdis[j]*zdis[j];

    if(dis2<rmax2[in]){
      p = (dis2-rmint2[in])*dx2tab[in];
      k = (int)p;
      if(k<1){
#ifdef NOERROR
	ierror++;
#endif
	jer = j;k=1;p=1.;
      }
      
      p = p-(double)k;
      dvt = dvtab[in];
      fm1 = dvt[k-1]; f0=dvt[k]; f1=dvt[k+1];
#ifdef CUBIC
      f2  = dvt[k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      if(dis2>rminh2[in]){
        dvt = sw[in];
        fm1 = dvt[k-1];	f0=dvt[k];  f1=dvt[k+1];
#ifdef CUBIC
        f2  = dvt[k+2];
        swv = f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
           p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
        swv = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      } else {
        swv = 1.;
      }      
#ifdef VERBOSE 
      printf("force[%d] = %f (%f %d %f %d)\n",j,f,swv,k,p,in);
#endif
      if(iflg==1){
        sprsl -= dis2*f;
        sprslt[0] -= xdis[j]*xdis[j]*f;
        sprslt[1] -= ydis[j]*xdis[j]*f;
        sprslt[2] -= zdis[j]*xdis[j]*f;
        sprslt[4] -= ydis[j]*ydis[j]*f;
        sprslt[5] -= zdis[j]*ydis[j]*f;
        sprslt[8] -= zdis[j]*zdis[j]*f;
        
        xdis[j] *= f;
        ydis[j] *= f;  
        zdis[j] *= f;
        fxds[j] = xdis[j]*swv;
        fyds[j] = ydis[j]*swv;
        fzds[j] = zdis[j]*swv;
        
        vt = vtab[in];
        fm1=vt[k-1];	f0 =vt[k];	f1 =vt[k+1];
#ifdef CUBIC
        f2  = vt[k+2];
        fcc = f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
           p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
        fcc = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
        spotlt  += fcc;
        spothlt += fcc-vt[ntm4];
#ifdef VERBOSE 
        printf("potential[%d] = %f\n",j,fcc/KCAL);
#endif
      } else {
        f *= swv;
        xdis[j] *= f;
        ydis[j] *= f;  
        zdis[j] *= f;
      }
    } else {
      xdis[j] = 0.;      ydis[j] = 0.;      zdis[j] = 0.;
      fxds[j] = 0.;      fyds[j] = 0.;      fzds[j] = 0.;
    }
  }
  if(iflg==1){
    sprslt[3] = sprslt[1];
    sprslt[6] = sprslt[2];
    sprslt[7] = sprslt[5];
  }
  spotl += spotlt;
  spothl += spothlt;
  
  /* recursion use breakpoints (offsets) to vectorize inner loop */
  for(k=0;k<nbreak;k++){
#ifdef VERBOSE
    printf("break point %d from %d to %d\n",k,ibreak[k],ibreak[k+1]);
#endif
#include "vectorize.h"
    for(j=ibreak[k];j<ibreak[k+1];j++){
      ni=indx[j];
      nj=jndx[j];
#ifdef VERBOSE
      printf("adding force pairs %d %d-%d\n",j,indx[j],jndx[j]);
#endif
      scr[ni        ] -= xdis[j];
      scr[ni+natoms ] -= ydis[j];
      scr[ni+natoms2] -= zdis[j];
      scr[nj        ] += xdis[j];
      scr[nj+natoms ] += ydis[j];
      scr[nj+natoms2] += zdis[j];
      if(iflg==1){
        rec[ni        ] -= fxds[j];
        rec[ni+natoms ] -= fyds[j];
        rec[ni+natoms2] -= fzds[j];
        rec[nj        ] += fxds[j];
        rec[nj+natoms ] += fyds[j];
        rec[nj+natoms2] += fzds[j];
      }
    }
  }

  /* free temporary vectors */
  free(ityp);
  free(xdis);
  
#ifdef NOERROR
  if(ierror){
    sprintf(line,"%d distance out of range in force routine\n",ierror);
    sprintf(line,"%s (ie atoms %d and %d)",line,indx[jer],jndx[jer]);
    md_warning(line);
  }
#endif
}
/*-----------------------------------------------------------------------*/
void vinter_npol(int nlen,int *indx,int *jndx,int *itype,
		 double *x,double *y,double *z,double *hmat,int iperd,
		 int **imap,double *rmin_cut2,double *rmax_cut2,
		 double *dx2tab,int ntable,
		 double **vtab,double **dvtab,double *spot,double *spoth,
		 double *sprs,double sprst[9],int ivol)
{
  int i,j,k,*ityp,*jtyp,intr,ierror,jer,in,jn;
  double *xdis,*ydis,*zdis,dis2,p;
  double fm1,f0,f1,f,spotl,spothl,sprsl;
  double sprslt[9];
  double hmati[9];
  LINE line;
#ifdef CUBIC
  double f2;
#endif

  jer = -1;
  ierror = 0;
  gethinv9(hmat,hmati);
  ityp = (int *)cmalloc(2*nlen*sizeof(int));
  jtyp = ityp + nlen;
  xdis = (double *)cmalloc(3*nlen*sizeof(double));
  ydis = xdis +   nlen;
  zdis = xdis + 2*nlen;
  spotl = spothl = sprsl = 0.0;
  for(i=0;i<9;i++)  sprslt[i]=0.0;

#include "vectorize.h"
  for(j=0;j<nlen;j++){
    ityp[j] = itype[indx[j]]; 
    jtyp[j] = itype[jndx[j]];
  }
#include "vectorize.h"
  for(j=0;j<nlen;j++){
    in = indx[j] ; jn = jndx[j];
    xdis[j] = x[in] - x[jn];
    ydis[j] = y[in] - y[jn];
    zdis[j] = z[in] - z[jn];
  }
  
  period(nlen,xdis,ydis,zdis,hmat,hmati,iperd,ivol);
      
#include "vectorize.h"
  for(j=0;j<nlen;j++){
    intr = imap[ityp[j]][jtyp[j]];
    dis2 = xdis[j]*xdis[j]+ydis[j]*ydis[j]+zdis[j]*zdis[j];
    
    if(dis2<rmax_cut2[intr]){
      p = (dis2-rmin_cut2[intr])*dx2tab[intr];
      k = (int)p;
      if(k<1) {
#ifdef NOERROR
	ierror= 1; 
#endif
	jer = j;k = 1; p=1.;
      }
      
      p = p-(double)k;
      fm1=vtab[intr][k-1];
      f0 =vtab[intr][k];
      f1 =vtab[intr][k+1];
      
#ifdef CUBIC
      f2  = vtab[intr][k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      spotl  += f;
      spothl += f-vtab[intr][ntable-4];
      
#ifdef PRINT_VDW_ENERGIES
      printf("VDW %d-%d = %g\n",indx[j],jndx[j],f);
#endif
      fm1 = dvtab[intr][k-1];
      f0  = dvtab[intr][k];
      f1  = dvtab[intr][k+1];
#ifdef CUBIC
      f2  = dvtab[intr][k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      sprsl -= dis2*f;
      sprslt[0] -= xdis[j]*xdis[j]*f;
      sprslt[1] -= ydis[j]*xdis[j]*f;
      sprslt[2] -= zdis[j]*xdis[j]*f;

      sprslt[3] -= xdis[j]*ydis[j]*f;
      sprslt[4] -= ydis[j]*ydis[j]*f;
      sprslt[5] -= zdis[j]*ydis[j]*f;

      sprslt[6] -= xdis[j]*zdis[j]*f;
      sprslt[7] -= ydis[j]*zdis[j]*f;
      sprslt[8] -= zdis[j]*zdis[j]*f;

#ifdef VERBOSE 
      printf("distance %d %g %d %g %g\n",j,sqrt(dis2),intr,
	     sqrt(rmax_cut2[intr]),f);
#endif
    }
  }
#ifdef NOERROR
  if(ierror){
    sprintf(line,"distance out of range in energy routine\n");
    sprintf(line,"%s between atoms %d and %d)",line,indx[jer],jndx[jer]);
    md_warning(line);
  }
#endif

  *spot += spotl;
  *spoth += spothl;
  *sprs += sprsl;
  for(i=0;i<9;i++)  sprst[i] += sprslt[i];

  /* free temporary vectors */
  free(ityp);
  free(xdis);

#ifdef VERBOSE
  printf("energies %g %g %g\n",*spot,*spoth,*sprs);
#endif

}

/*-----------------------------------------------------------------------*/
void getsvirter(double *vhinter,double *vinter,double *winter,
		double wintert[9])
{
  int i;
#ifdef PARA
  double sprst[12],sprlt[12];

  for(i=0;i<9;i++) sprst[i] = sprslt[i];
  sprst[9] = spothl;
  sprst[10] = spotl;
  sprst[11] = sprsl;
  MPI_Allreduce(sprst,sprlt,12,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=0;i<9;i++) sprslt[i] = sprlt[i];
  spothl = sprlt[9];
  spotl = sprlt[10];
  sprsl = sprlt[11];
#endif

  *vhinter += spothl;
  *vinter  += spotl;
  *winter  += sprsl;
  for(i=0;i<9;i++) wintert[i] += sprslt[i];
}
