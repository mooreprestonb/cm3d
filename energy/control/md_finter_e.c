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

/* #define VERBOSE */
/* #define CUBIC */
/* #define ANALYTIC */

#ifdef ANALYTIC
extern double alpha_ewald;
#endif
int exclij(int i,int j,int *nexcl,int **excl);

/*-----------------------------------------------------------------------*/
void check_distance_e(SIMPARMS *simparms,COORDS *coords,INTER *inter)
{
  int i,j,ii,jj,*excl,nexcl;
  double xdis,ydis,zdis,rdis2;
  LINE line;

  for(i=0;i<coords->coul.ncharge-1;i++){
    ii = coords->coul.icharge[i];
    excl = &(coords->coul.exclude[ii][0]);
    nexcl = *(excl);
    for(j=i+1;j<coords->coul.ncharge;j++){
      jj = coords->coul.icharge[j];
      if(nexcl == jj){
	nexcl = *(++excl);
      } else {
	xdis = coords->px[ii] - coords->px[jj];
	ydis = coords->py[ii] - coords->py[jj];
	zdis = coords->pz[ii] - coords->pz[jj];

	period(1,&xdis,&ydis,&zdis,coords->hmat,coords->hmati,simparms->iperd,simparms->ivol);
	rdis2 = xdis*xdis+ydis*ydis+zdis*zdis;

        if(rdis2<inter->rmins_coul2){
          sprintf(line,"distance out of range in coulumb force\n\t");
	  sprintf(line,"%s distance = %g (between atoms %d and %d)\n\t",line,
		  sqrt(rdis2),ii,jj);
	  sprintf(line,"%s Distance should be > %g (interaction max=%g)",line,
		  sqrt(inter->rmins_coul2),sqrt(inter->rmaxl_coul2));
	  md_error(line);
        }
      }
    }
  }
}
/*--------------------------------------------------------------------*/
void force_ter_e(SIMPARMS *simparms,COORDS *coords,INTER *inter,
		 NGBR *ngbr,int iflg)
{
  int i,j,ii,jj,k,iperd,ivol,natoms,ipolar;
  int np,*inp,*jnp,nexcl,*excl,ncharge,icell,jcell;
  long ibegin,iend;
  double *fx,*fxs,*rbuf,*sbuf;
#ifdef CUBIC
  double f2;
#endif
#ifdef ANALYTIC
  double dis,sm,br,rheal_res;
  rheal_res = sqrt(inter->rminl_coul2)-sqrt(inter->rmaxs_coul2);
#endif

  natoms= simparms->natoms;
  iperd = simparms->iperd;
  ipolar = simparms->ipolar;
  ivol  = simparms->ivol;
  np = 3*natoms;
  for(i=0;i<np;i++) coords->scr_buf[i] = coords->scr_rec[i]=0.0;
  ncharge = coords->coul.ncharge;

  np = 0;
  inp = jnp = (int *)NULL;
  fx  = coords->fxl;
  fxs = coords->fxs;

  switch(iflg){
  case 0:
    np = ngbr->nprse; 
    inp = ngbr->inpse; 
    jnp = ngbr->jnpse; 
    fx = coords->fxs;
    break;
  case 1:
    zero_pote();
    np  = ngbr->nprte;
    inp = ngbr->inpte;
    jnp = ngbr->jnpte; 
    break;
  default:
    md_error("in electronic inter-molecular force routine?!?!?"); 
    break;
  }
  
  switch(ngbr->ilist){
     case 0:
       inp = coords->indx;
       jnp = coords->jndx;
       
       if(ncharge<=simparms->size) {
         ibegin = MIN(simparms->rank,ncharge-1);
         iend = MIN(simparms->rank+1,ncharge-1);
       } else {
         decomp_trig((long)ncharge,simparms->size,simparms->rank,&ibegin,&iend);
       }
       k=0;
       for(i=ibegin;i<iend;i++){
         ii = coords->coul.icharge[i];
         excl = &(coords->coul.exclude[ii][0]);
         nexcl = *(excl);
         for(j=i+1;j<ncharge;j++){
           jj = coords->coul.icharge[j];
           if(nexcl == jj){
             nexcl = *(++excl);
           } else {
             inp[k] = ii;
             jnp[k] = jj;
             k++;
             if(k==simparms->nlen){
               force_e(natoms,k,inp,jnp,coords,inter,iperd,ivol,iflg,ipolar);
               k=0;
             }
           }
         }
       }
       if(k>0) {
         force_e(natoms,k,inp,jnp,coords,inter,iperd,ivol,iflg,ipolar);
       }
       break;
     case 1:
       for(i=0;i<np;i+=simparms->nlen){
         k = MIN(simparms->nlen,np-i);
         force_e(natoms,k,&inp[i],&jnp[i],coords,inter,iperd,ivol,iflg,ipolar);
       }
       break;
     case 2:
       get_lnklist_e(simparms->natoms,ngbr->ncells,ngbr->hlist,ngbr->lnklist,
          iperd,coords->hmat,
          coords->px,coords->py,coords->pz,coords->qch);
       inp = coords->indx;
       jnp = coords->jndx;    
       
       k = 0;
       for(i=0;i<ngbr->npaircelle;i++){
         icell = ngbr->icelle[i];
         jcell = ngbr->jcelle[i];
         
         ii = ngbr->hlist[icell];
         
         while(ii!=-1){
           if(icell==jcell){
             jj = ngbr->lnklist[ii];
           } else {
             jj = ngbr->hlist[jcell];
           }
           while(jj!=-1){
             if(!exclij(ii,jj,coords->coul.nexclude,coords->coul.exclude)){ 
               inp[k] = ii;
               jnp[k] = jj;
               k++;
               
               if(k==simparms->nlen){
                 force_e(natoms,k,inp,jnp,coords,inter,iperd,ivol,iflg,ipolar);
                 k=0;
               }
             }
             jj = ngbr->lnklist[jj];
           }
           ii = ngbr->lnklist[ii];
         }
       }
       if(k>0) force_e(natoms,k,inp,jnp,coords,inter,iperd,ivol,iflg,ipolar);
       break;    
  }
  
  np = natoms*3;
#ifdef PARA
  rbuf = cmalloc(2*np*sizeof(double));
  sbuf = rbuf + np;

  MPI_Allreduce(coords->scr_buf,sbuf,np,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if(iflg==1)
    MPI_Allreduce(coords->scr_rec,rbuf,np,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
  sbuf = coords->scr_buf;
  rbuf = coords->scr_rec;
#endif

#include "vectorize.h"
  for(i=0;i<np;i++) fx[i] += sbuf[i];
  if(iflg==1){
#include "vectorize.h"
    for(i=0;i<np;i++) fxs[i] += rbuf[i];
  }

#ifdef PARA
  free(rbuf);
#endif
  
#ifdef VERBOSE
  printf("\nprinting forces\n");
  for(i=0;i<n3;i++)  printf("%d %g %g %g %g\n",i,fx[i],sbuf[i],fxs[i],rbuf[i]);
#endif
}
/*-----------------------------------------------------------------------*/
void getvirter_e(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr,
		 double *velec,double *welec,double welectensor[9],int iflg)
{
  int i,j,k,ii,jj,*inp,*jnp,nexcl,*excl,ncharge,icell,jcell;
  long ibegin,iend;
  double swelec,spelec,swelec_local,spelec_local,spelech_local;
  double *px,*py,*pz,*qch;
  double swelec_localtensor[9];
  double swelectensor[9];
#ifdef CUBIC
  double f2;
#endif

  swelec_local = spelec_local = spelech_local = 0.;
  for(i=0;i<9;i++) swelec_localtensor[i]=0.0;

  ncharge = coords->coul.ncharge;
  px = coords->px;
  py = coords->py;
  pz = coords->pz;
  qch= coords->qch; 

  switch(ngbr->ilist){
  case 0:
    inp = coords->indx;
    jnp = coords->jndx;

    if(ncharge<=simparms->size) {
      ibegin = MIN(simparms->rank,ncharge-1);
      iend = MIN(simparms->rank+1,ncharge-1);
    } else {
      decomp_trig((long)ncharge,simparms->size,simparms->rank,&ibegin,&iend);
    }

    k=0;
    for(i=ibegin;i<iend;i++){
      ii = coords->coul.icharge[i];
      excl = &(coords->coul.exclude[ii][0]);
      nexcl = *(excl);
      for(j=i+1;j<ncharge;j++){
        jj = coords->coul.icharge[j];
        if(nexcl == jj){ 
          nexcl = *(++excl);
        } else {
          inp[k] = ii;
          jnp[k] = jj;
          k++;
          if(k==simparms->nlen){
            vinter_e(k,inp,jnp,coords,inter,simparms->iperd,
		     &spelec_local,&spelech_local,&swelec_local,
		     swelec_localtensor,simparms->ivol,iflg,simparms->ipolar);
            k=0;
          }
        }
      }
    }
    if(k>0) {
      vinter_e(k,inp,jnp,coords,inter,simparms->iperd,
	       &spelec_local,&spelech_local,&swelec_local,
	       swelec_localtensor,simparms->ivol,iflg,simparms->ipolar);
    }
    break;
  case 1:
    for(i=0;i<ngbr->nprte;i+=simparms->nlen){
      k = MIN(simparms->nlen,ngbr->nprte-i);
      vinter_e(k,&ngbr->inpte[i],&ngbr->jnpte[i],coords,inter,
	       simparms->iperd,&spelec_local,&spelech_local,&swelec_local,
	       swelec_localtensor,simparms->ivol,iflg,simparms->ipolar);
    }
    break;
     case 2:
       get_lnklist_e(simparms->natoms,ngbr->ncells,ngbr->hlist,ngbr->lnklist,
          simparms->iperd,coords->hmat,
          coords->px,coords->py,coords->pz,coords->qch);
       inp = coords->indx;
       jnp = coords->jndx;    
       
       k = 0;
       for(i=0;i<ngbr->npaircelle;i++){
         icell = ngbr->icelle[i];
         jcell = ngbr->jcelle[i];
         
         ii = ngbr->hlist[icell];
         
         while(ii!=-1){
           if(icell==jcell){
             jj = ngbr->lnklist[ii];
           } else {
             jj = ngbr->hlist[jcell];
           }
           while(jj!=-1){
             if(!exclij(ii,jj,coords->coul.nexclude,coords->coul.exclude)){ 
               inp[k] = ii;
               jnp[k] = jj;
               k++;
               if(k==simparms->nlen){
		 vinter_e(k,inp,jnp,coords,inter,simparms->iperd,
			  &spelec_local,&spelech_local,&swelec_local,
			  swelec_localtensor,simparms->ivol,iflg,simparms->ipolar);
                 k=0;
               }
             }
             jj = ngbr->lnklist[jj];
           }
           ii = ngbr->lnklist[ii];
         }
       }
       if(k>0){
	 vinter_e(k,inp,jnp,coords,inter,simparms->iperd,
		  &spelec_local,&spelech_local,&swelec_local,
		  swelec_localtensor,simparms->ivol,iflg,simparms->ipolar);
       }
       break;
  }
  
#ifdef PARA
  MPI_Allreduce(&swelec_local ,&swelec ,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(swelec_localtensor ,swelectensor ,9,MPI_DOUBLE,MPI_SUM,
		MPI_COMM_WORLD);
  MPI_Allreduce(&spelec_local ,&spelec ,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  /*  MPI_Allreduce(&spelech_local,&spelech,1,
      MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
#else
  swelec = swelec_local;
  for(i=0;i<9;i++) swelectensor[i] = swelec_localtensor[i];

  spelec = spelec_local;
#endif
  *welec += swelec; 
  for(i=0;i<9;i++) welectensor[i] += swelectensor[i];
  *velec += spelec;
}
/*-----------------------------------------------------------------------*/
