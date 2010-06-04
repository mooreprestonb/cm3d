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


#include "md.h"

/* #define PRINT_FORCE */

/* #define VERBOSE2 */
/* #define DEBUG */
/* #define KCALS */
/* #define CHECK */

void ptensorchk(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr,int iflg);
void ngradv(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr,int iflg);
void get_p2mt(int,double[],double[],double[],double[],double[9]);

enum funcs {BOND,BEND,XBOND};

/*------------------------------------------------------------------*/
void force(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr,int iflg)
{
  int i,n3;
  double *fxp;
#ifdef PARA
  double *rec_buf;
#endif

  n3 = simparms->natoms*3;

  switch(iflg){
  case 0:  /* short range intermolecular */
    if(simparms->iperd!=0)gethinv9(coords->hmat,coords->hmati);
    fxp = coords->fxs;
    for(i=0;i<n3;i++) fxp[i] = 0.;
    force_ter(simparms,coords,inter,ngbr,iflg);
    force_ter_e(simparms,coords,inter,ngbr,iflg);
    rmass(simparms->natoms,coords->fxs,coords->fys,coords->fzs,coords->amass);
    break;

  case 1: /* long range intermolecular */
    if(simparms->iperd!=0) gethinv9(coords->hmat,coords->hmati);
    fxp = coords->fxs;  for(i=0;i<n3;i++) fxp[i] = 0.;
    fxp = coords->fxl;  for(i=0;i<n3;i++) fxp[i] = 0.;
    force_ter(simparms,coords,inter,ngbr,iflg);
    force_ter_e(simparms,coords,inter,ngbr,iflg);
    if(simparms->iperd == 3){    /* ewald sums */
      fk_ewald(simparms,coords);
      ecorr(simparms,coords);
    }
    rmass(simparms->natoms,coords->fxl,coords->fyl,coords->fzl,coords->amass);
    rmass(simparms->natoms,coords->fxs,coords->fys,coords->fzs,coords->amass);
    break;

  case 2: /* intra molecular */
    fxp = coords->fxa;
    for(i=0;i<n3;i++) fxp[i] = 0.;
    fbond(simparms,coords);
    fbend(simparms,coords);
    fbondx(simparms,coords);
    /* should we zero the colvar forces */
    fcolvar_cntr(simparms,coords,1); /* 1 is to add force to particles */

#ifdef TWO
    ftors(simparms,coords);
    fonfo(simparms,coords); 
    if(simparms->iextern==1) force_extern(simparms,coords,inter->ntable);
#endif    

#ifdef PARA
    rec_buf = coords->scr_rec;
    for(i=0;i<n3;i++) rec_buf[i] = 0.;
    MPI_Allreduce(fxp,rec_buf,n3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<n3;i++) fxp[i] = rec_buf[i];
#endif
    rmass(simparms->natoms,coords->fxa,coords->fya,coords->fza,coords->amass);
    break;

  case 3: /* (torsion & 1-4) molecular */
    fxp = coords->fxr;
    for(i=0;i<n3;i++) fxp[i] = 0.;

#ifndef TWO
    ftors(simparms,coords);
    fonfo(simparms,coords); 
    if(simparms->iextern==1) force_extern(simparms,coords,inter->ntable);
#ifdef PARA
    rec_buf = coords->scr_rec;
    for(i=0;i<n3;i++) rec_buf[i] = 0.;
    MPI_Allreduce(fxp,rec_buf,n3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<n3;i++) fxp[i] = rec_buf[i];
#endif
#endif
    rmass(simparms->natoms,coords->fxr,coords->fyr,coords->fzr,coords->amass);
    break;
    
  default:
    fprintf(stderr,"ERROR: in forces!!!????\n");
    exit(1);
  }
  if(simparms->nfreeze>0) freeze_atom(simparms,coords);
  if(simparms->iensemble==3) {
    fxp = coords->vx;
    for(i=0;i<n3;i++) fxp[i] = 0.;
  }

#ifdef PRINT_FORCE
  printf("iflag = %d\n",iflg);
  for(i=0;i<n3;i++) printf("%d %g\n",i,fxp[i]);
#endif
}
/*----------------------------------------------------------------*/
void getvireal(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies,
	       INTER *inter,NGBR *ngbr,int iflg)
{
  int i;
  double rvol;

  energies->vinter = energies->vhinter = 0.0;
  energies->velec = energies->vewald = 0.0;
  energies->vbond = energies->vbend = energies->vtors = energies->vbondx = 0.;
  energies->vonfo = energies->vonfo_e = 0.0;
  energies->vextern = energies->vextern_e = 0.0;
  energies->winter = energies->welec = energies->wewald  = 0.0;
  energies->wbond = energies->wbend = energies->wtors = energies->wbondx = 0.;
  energies->wonfo = energies->wonfo_e = 0.0;
  energies->wextern = energies->wextern_e = 0.0;
  energies->US = energies->U = energies->UI = 0.0;
  energies->W = energies->WI = 0.0;
  for(i=0;i<9;i++) {
    energies->wbondxtensor[i]=0.0;
    energies->WItensor[i]= energies->Wtensor[i]=0.0;
    energies->wbondtensor[i]= energies->wbendtensor[i]=0.0;
    energies->wonfotensor[i]= energies->wonfo_etensor[i]=0.0;
    energies->wtorstensor[i]= energies->wintertensor[i]=0.0;
    energies->welectensor[i]= energies->wewaldtensor[i]=0.0;
  }

  if(iflg==1) {
    getsvirter(&energies->vhinter,&energies->vinter,&energies->winter,
	       energies->wintertensor);
    getsvirter_e(&energies->velec,&energies->welec,energies->welectensor);
  } else {
    getvirter(simparms,coords,inter,ngbr,&energies->vhinter,&energies->vinter,
	      &energies->winter,energies->wintertensor);
    getvirter_e(simparms,coords,inter,ngbr,&energies->velec,&energies->welec,
		energies->welectensor,iflg);
  }

  getvbond(simparms,coords,
	   &energies->vbond,&energies->wbond,energies->wbondtensor);
  getvbend(simparms,coords,
	   &energies->vbend,&energies->wbend,energies->wbendtensor);
  getvtors(simparms,coords,
	   &energies->vtors,&energies->wtors,energies->wtorstensor);
  getvonfo(simparms,coords,
	   &energies->vonfo,&energies->vonfo_e,&energies->wonfo,
	   &energies->wonfo_e,energies->wonfotensor,energies->wonfo_etensor);
  getvbondx(simparms,coords,
	    &energies->vbond,&energies->wbond,energies->wbondxtensor);
  
  if(simparms->iextern == 1) {
    getvextern(simparms,coords,&energies->vextern,&energies->wextern,
	       &energies->vextern_e,&energies->wextern_e,inter->ntable);
  }

  if(simparms->iperd == 3) {
    getvewald(simparms,coords,&energies->vewald,&energies->wewald,
	      energies->wewaldtensor);
  }

#ifdef CHECK
  ptensorchk(simparms,coords,inter,ngbr,iflg);
  ngradv(simparms,coords,inter,ngbr,iflg);
#endif

  energies->US = (energies->vhinter + energies->velec + energies->vewald +
		  energies->vextern + energies->vextern_e);
  energies->U  = (energies->vinter  + energies->velec + energies->vewald + 
		  energies->vextern + energies->vextern_e);
  energies->W  = (energies->winter  + energies->welec + energies->wewald +
		  energies->wextern + energies->wextern_e);
  energies->UI = (energies->vbond + energies->vbend + energies->vtors + 
		  energies->vonfo + energies->vonfo_e + energies->vbondx);
  energies->WI = (energies->wbond + energies->wbend + energies->wtors + 
		  energies->wonfo + energies->wonfo_e + energies->wbondx);

  for(i=0;i<9;i++) {    
    energies->WItensor[i]=(energies->wbondtensor[i]+energies->wbendtensor[i]+
			   energies->wtorstensor[i]+energies->wonfotensor[i]+
			   energies->wonfo_etensor[i]+
			   energies->wbondxtensor[i]);
    energies->Wtensor[i]=(energies->wintertensor[i]+energies->welectensor[i]+
			  energies->wewaldtensor[i]);
  }
  /* add external field to the z compent of the pressure tensor */
  energies->Wtensor[8] += energies->wextern + energies->wextern_e;

  /* add long range correction if iperd==3 */
  if(simparms->iperd == 3){
    rvol = 1./get_deth(coords->hmat);
    energies->W -= DIM*simparms->wlrc*rvol;
    energies->Wtensor[0] -= simparms->wlrc*rvol;
    energies->Wtensor[4] -= simparms->wlrc*rvol;
    energies->Wtensor[8] -= simparms->wlrc*rvol;
  }


#ifdef DEBUG
#ifdef KCALS
  printf("bonds = %g Kcal/mol\n",energies->vbond/KCAL);
  printf("Xbonds = %g Kcal/mol\n",energies->vbondx/KCAL);
  printf("bends = %g Kcal/mol\n",energies->vbend/KCAL);
  printf("tors  = %g Kcal/mol\n",energies->vtors/KCAL);
  printf("onfo  = %g Kcal/mol\n",energies->vonfo/KCAL);
  printf("onfo_e= %g Kcal/mol\n",energies->vonfo_e/KCAL);
  printf("inter = %g Kcal/mol\n",energies->vinter/KCAL);
  printf("elec  = %g Kcal/mol\n",energies->velec/KCAL);
  printf("recip = %g Kcal/mol\n",energies->vewald/KCAL);
#else
  printf("vhinter=%g vinter=%g velec=%g vrecip=%g\n",
	 energies->vhinter,energies->vinter,energies->velec,energies->vewald);
  printf("vbond=%g vbondx=%g vbend=%g vtors=%g vonfo=%g vonfo_e=%g\n",
         energies->vbond,energies->vbondx,energies->vbend,energies->vtors,
         energies->vonfo,energies->vonfo_e);
  printf("winter=%g welec=%g wrecip=%g\n",
	 energies->winter/ECONV,energies->welec/ECONV,energies->wewald/ECONV);
  printf("wbond=%g wbondx=%g wbend=%g wtors=%g wonfo = %g wonfo_e=%g\n",
         energies->wbond/ECONV,energies->wbondx/ECONV,energies->wbend/ECONV,energies->wtors/ECONV,
         energies->wonfo/ECONV,energies->wonfo_e/ECONV);
#endif
  exit(1);
#endif
}
/*----------------------------------------------------------------*/
void geteng(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies)
{
  int i;
  double p2m,vol,pv_tensor[9],p2m_tensor[9];

  for(i=0;i<9;i++)  pv_tensor[i]= p2m_tensor[i]=0.0;

  vol = get_deth(coords->hmat);

  energies->poth  = energies->US;
  energies->poter = energies->U;
  energies->potra = energies->UI;

  if(simparms->iperd == 3){
    energies->poter += simparms->vlrc/vol; 
    energies->poth  += simparms->vlrc/vol;
  }

  energies->pot_inter = energies->vinter;
  energies->pot_bond  = energies->vbond + energies->vbondx;
  energies->pot_bend  = energies->vbend;
  energies->pot_tors  = energies->vtors;
  energies->pot_onfo  = energies->vonfo;
  energies->pot_onfo_e= energies->vonfo_e;
  energies->pot_elec  = energies->velec;
  energies->pot_recip = energies->vewald;
  energies->pot_extern= energies->vextern;
  energies->pot_extern_e= energies->vextern_e;

  get_p2mt(simparms->natoms,coords->amass,coords->vx,coords->vy,coords->vz,
	   p2m_tensor);

  p2m = p2m_tensor[0]+p2m_tensor[4]+p2m_tensor[8];

  energies->prsi  = (p2m+energies->W+energies->WI)/(DIM*vol);
  for(i=0;i<9;i++){
    pv_tensor[i]=p2m_tensor[i]+energies->Wtensor[i]+energies->WItensor[i];
  }
  for(i=0;i<9;i++) energies->prs_tensor[i] = pv_tensor[i]/vol;

#ifdef DEBUG
  printf("pressure components %g %g %g %g,",p2m,energies->W,energies->WI,vol);
  printf(" pressure now = %g == %g\n",energies->prsi,
	 (pv_tensor[0]+pv_tensor[4]+pv_tensor[8])/(3.*vol));
#endif

  energies->zke   = .5*p2m;
  energies->tiout = p2m/(double)(simparms->ndof);

}
/*--------------------------------------------------------------------*/

void freeze_atom(SIMPARMS *simparms,COORDS *coords)
{
  int i,j;

  for(i=0;i<simparms->nfreeze;i++){
    j = coords->ifreeze[i];
    coords->vx[j] = coords->vy[j] = coords->vz[j] = 0.;

    coords->ffreeze[i*DIM] = coords->fxs[j]+coords->fxl[j]+coords->fxa[j]+coords->fxr[j]+coords->fxt[j];
    coords->ffreeze[i*DIM+1] = coords->fys[j]+coords->fyl[j]+coords->fya[j]+coords->fyr[j]+coords->fyt[j];
    coords->ffreeze[i*DIM+2] = coords->fzs[j]+coords->fzl[j]+coords->fza[j]+coords->fzr[j]+coords->fzt[j];


    coords->fxs[j] = coords->fys[j] = coords->fzs[j] = 0.;
    coords->fxl[j] = coords->fyl[j] = coords->fzl[j] = 0.;
    coords->fxa[j] = coords->fya[j] = coords->fza[j] = 0.;
    coords->fxr[j] = coords->fyr[j] = coords->fzr[j] = 0.;
    coords->fxt[j] = coords->fyt[j] = coords->fzt[j] = 0.;
  }
}
/*--------------------------------------------------------------------*/
void get_p2mt(int natoms,double amass[],double vx[],double vy[],double vz[],
	      double p2m_tensor[9])
{
  int i;
  double am;
  for(i=0;i<9;i++) p2m_tensor[i] = 0.;
#include "vectorize.h"
  for(i=0;i<natoms;i++){
    am = amass[i];

    p2m_tensor[0] += am*(vx[i]*vx[i]);
    p2m_tensor[4] += am*(vy[i]*vy[i]);
    p2m_tensor[8] += am*(vz[i]*vz[i]);

    p2m_tensor[1] += am*(vy[i]*vx[i]);
    p2m_tensor[2] += am*(vz[i]*vx[i]);
    p2m_tensor[5] += am*(vz[i]*vy[i]);
  }

  p2m_tensor[3] = p2m_tensor[1];
  p2m_tensor[6] = p2m_tensor[2];
  p2m_tensor[7] = p2m_tensor[5];
}

/*--------------------------------------------------------------------*/
