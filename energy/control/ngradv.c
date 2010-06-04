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
#define DELTA (1.e-5)

/* calculate first derivative numericaly */

void ngradv(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr,int iflg)
{
  int i,j,n,r,c,np,iloop;
  double *sx,*sy,*sz;  /* stored positions */
  double *x,*y,*z,*fkx,*fky,*fkz;
  double v,v_e,w,w_e,vup,vlow;
  double *hmat,hinv[9],tensor[9],etensor[9];
       
  md_stdout("checking forces numerically");

  np = simparms->natoms;
  sx = (double *)calloc(6*np,sizeof(double));
  sy = sx +   np;
  sz = sx + 2*np;
  fkx = sx+ 3*np;
  fky = sx+ 4*np;
  fkz = sx+ 5*np;
  x = coords->px;  y = coords->py;  z = coords->pz;
  hmat = coords->hmat;
  
  for(i=0;i<9;i++) tensor[i] = etensor[i] = 0.;
  /* store particle posistions  and zero forces */
  for(n=0;n<np;n++) { 
    sx[n]= x[n];    sy[n]= y[n];    sz[n]= z[n];
  }

  gethinv9(hmat,hinv);
     
  for(iloop=0;iloop<10;iloop++){
    for(n=0;n<np;n++) fkx[n] = fky[n] = fkz[n] = 0.;
    for(c=0;c<np;c++){
      for(r=0;r<3;r++){
	switch (r){
	case 0:
	  x[c] = sx[c] + DELTA;
	  break;
	case 1:
	  y[c] = sy[c] + DELTA;
	  break;
	case 2:
	  z[c] = sz[c] + DELTA;
	  break;
	}
	
	/* get pot high */
	v = w = v_e = w_e = 0.;
	switch (iloop) {
	case 0:
	  getvbond(simparms,coords,&v,&w,tensor);
	  break;
	case 1:
	  getvbondx(simparms,coords,&v,&w,tensor);
	  break;
	case 2:
	  getvbend(simparms,coords,&v,&w,tensor);
	  break;
	case 3:
	  getvtors(simparms,coords,&v,&w,tensor);
	  break;
	case 4:
	  getvonfo(simparms,coords,&v,&v_e,&w,&w_e,tensor,etensor);
	  v += v_e; w += w_e;
	  break;
	case 5:
	  getvirter(simparms,coords,inter,ngbr,&v_e,&v,&w,tensor);
	  break;
	case 6:
	  getvirter_e(simparms,coords,inter,ngbr,&v,&w,tensor,iflg);
	  break;
	case 7:
	  if(simparms->iperd==3){
	    fk_ewald(simparms,coords);
	    v = coords->coul.vk_ewald;
	  }
	  break;
	case 8:
	  if(simparms->iperd==3){
	    ecorr(simparms,coords);
	    v = coords->coul.pot_ecorr;
	  }
	  break;
	case 9:
	  if(simparms->iextern){
	    getvextern(simparms,coords,&v,&w,&v_e,&w_e,inter->ntable);
	    v += v_e;
	  }
	  break;
	} 
	
	vup = v;
	
	switch (r){
	case 0:
	  x[c] = sx[c] - DELTA;
	  break;
	case 1:
	  y[c] = sy[c] - DELTA;
	  break;
	case 2:
	  z[c] = sz[c] - DELTA;
	  break;
	}
	/* get pot low */
	v = w = v_e = w_e = 0.;
	switch (iloop) {
	case 0:
	  getvbond(simparms,coords,&v,&w,tensor);
	  break;
	case 1:
	  getvbondx(simparms,coords,&v,&w,tensor);
	  break;
	case 2:
	  getvbend(simparms,coords,&v,&w,tensor);
	  break;
	case 3:
	  getvtors(simparms,coords,&v,&w,tensor);
	  break;
	case 4:
	  getvonfo(simparms,coords,&v,&v_e,&w,&w_e,tensor,etensor);
	  v += v_e; w += w_e;
	  break;
	case 5:
	  getvirter(simparms,coords,inter,ngbr,&v_e,&v,&w,tensor);
	  break;
	case 6:
	  getvirter_e(simparms,coords,inter,ngbr,&v,&w,tensor,iflg);
	  break;
	case 7:
	  if(simparms->iperd==3){
	    fk_ewald(simparms,coords);
	    v = coords->coul.vk_ewald;
	  }
	  break;
	case 8:
	  if(simparms->iperd==3){
	    ecorr(simparms,coords);
	    v = coords->coul.pot_ecorr;
	  }
	  break;
	case 9:
	  if(simparms->iextern){
	    getvextern(simparms,coords,&v,&w,&v_e,&w_e,inter->ntable);
	    v += v_e;
	  }
	  break;
	case 10:
	  if(simparms->iextern)
	    getvextern(simparms,coords,&v_e,&w_e,&v,&w,inter->ntable);
	  break;
	} 
      
	vlow = v;
	
	/* calculate -derivative and reset pos */
	switch (r){
	case 0:
	  fkx[c] = -(vup-vlow)/(2.*DELTA); x[c] = sx[c];
	  break;
	case 1:
	  fky[c] = -(vup-vlow)/(2.*DELTA); y[c] = sy[c];
	  break;
	case 2:
	  fkz[c] = -(vup-vlow)/(2.*DELTA); z[c] = sz[c];
	  break;
	}
      }
    }
    /* calculate analytic forces */
    for(n=0;n<np;n++) { 
      coords->fxs[n] = coords->fys[n] = coords->fzs[n] = 0.0;
      coords->fxl[n] = coords->fyl[n] = coords->fzl[n] = 0.0;
      coords->fxr[n] = coords->fyr[n] = coords->fzr[n] = 0.0;
      coords->fxa[n] = coords->fya[n] = coords->fza[n] = 0.0;
      coords->fxt[n] = coords->fyt[n] = coords->fzt[n] = 0.0;
    }
    
    w = v = 0;
    switch (iloop) {
    case 0:
      printf("\nfbond \n");
      fbond(simparms,coords);
      break;
    case 1:
      printf("\nfbondx \n");
      fbondx(simparms,coords);
      break;
    case 2:
      printf("\nfbend \n"); 
      fbend(simparms,coords);
      break;
    case 3:
      printf("\nftors \n"); 
      ftors(simparms,coords);
      break;
    case 4:
      printf("\nfonfo \n"); 
      fonfo(simparms,coords);
      break;
    case 5:
      printf("\ngetvirter \n"); 
      force_ter(simparms,coords,inter,ngbr,1);
      break;
    case 6:
      printf("\ngetvirter_e\n"); 
      force_ter_e(simparms,coords,inter,ngbr,1);
      break;
    case 7:
      if(simparms->iperd==3){
	printf("\nfk_ewald \n"); 
	fk_ewald(simparms,coords);
      }
      break;
    case 8:
      if(simparms->iperd==3){
	printf("\necorr \n"); 
	ecorr(simparms,coords);
      }
      break;
    case 9:
      if(simparms->iextern==1){
	printf("\nget_extern\n"); 
	force_extern(simparms,coords,inter->ntable);
      }
      break;
    } 
    /* fxl is the long + short range force */
    for(j=0;j<np;j++) {
      coords->fxt[j] = (coords->fxl[j]+coords->fxa[j]+coords->fxr[j]);
      coords->fyt[j] = (coords->fyl[j]+coords->fya[j]+coords->fyr[j]);
      coords->fzt[j] = (coords->fzl[j]+coords->fza[j]+coords->fzr[j]);
    }
      
    for(i=0;i<np;i++) {
      if(fkx[i]!=0.0 && coords->fxt[i]!=0.0){
	printf("%d fx=% 11.5g fxn=% 11.5g diff=% 11.5g\n",i,coords->fxt[i],fkx[i],
	       fabs((coords->fxt[i]-fkx[i])/fkx[i]));
      }
      if(fky[i]!=0.0 && coords->fyt[i]!=0.0){
	printf("%d fy=% 11.5g fyn=% 11.5g diff=% 11.5g\n",i,coords->fyt[i],fky[i],
	       fabs((coords->fyt[i]-fky[i])/fky[i]));
      }
      if(fkz[i]!=0.0 && coords->fzt[i]!=0.0){
	printf("%d fz=% 11.5g fzn=% 11.5g diff=% 11.5g\n",i,coords->fzt[i],fkz[i],
	       fabs((coords->fzt[i]-fkz[i])/fkz[i]));
      }
    }
  }
  
  free(sx);
  exit(1);
  return;
}
