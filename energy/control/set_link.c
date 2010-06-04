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

/* create the link list cell  */

#include "md.h"

void set_link(SIMPARMS *simparms,NGBR *ngbr,double hmat[],
	      double rmax,double rmaxe)
{
  int icell,ic,jc,kc,i,j,k,in,jn,kn,nmax;
  int ncell,ncells,ncell2,iperd,num_cells,mem_now;
  long ibegin,iend;
  double rcell,dx,dy,dz,sinangle,ddot,r;
  static int istart=0,ncell_now,ncell_nowe;
  LINE line;
  
  iperd = simparms->iperd;
  ngbr->npaircell  =0;
  ngbr->npaircelle =0;
  nmax = (ncell=ngbr->ncells)/2;
  rcell = 1./(double)ncell;
  ncell2 = ncell*ncell;
  ncells = ncell*ncell*ncell;
  mem_now = simparms->mem_bytes;
  
  decomp1d((long)ncells,simparms->size,simparms->rank,&ibegin,&iend);
  num_cells = iend-ibegin;

  if(istart==0){
    if(simparms->rank==0){
      sprintf(line,"Setting up link list neighbors (~%d cells/processor)",
	      num_cells);
      md_stdout(line);
      if(ncells<simparms->size){
	sprintf(line,"%d cells is less then %d prossesors",
		num_cells,simparms->size);
	md_error(line);
      }
    }
    
    istart=1;
    count_link_pair(ngbr->ncells,simparms->iperd,hmat,(int)ibegin,(int)iend,
		    &ncell_now,&ncell_nowe,rmax,rmaxe);
    ngbr->icell = (int *)cmalloc(sizeof(int)*ncell_now);
    ngbr->jcell = (int *)cmalloc(sizeof(int)*ncell_now);
    ngbr->icelle = (int *)cmalloc(sizeof(int)*ncell_nowe);
    ngbr->jcelle = (int *)cmalloc(sizeof(int)*ncell_nowe);
    simparms->mem_bytes += 2*sizeof(int)*(ncell_now+ncell_nowe);
  } else {
    simparms->mem_bytes -= 2*sizeof(int)*(ncell_now+ncell_nowe);
    count_link_pair(ngbr->ncells,simparms->iperd,hmat,(int)ibegin,(int)iend,
		    &ncell_now,&ncell_nowe,rmax,rmaxe);
    ngbr->icell = (int *)realloc(ngbr->icell,sizeof(int)*ncell_now);
    ngbr->jcell = (int *)realloc(ngbr->jcell,sizeof(int)*ncell_now);
    ngbr->icelle = (int *)realloc(ngbr->icelle,sizeof(int)*ncell_nowe);
    ngbr->jcelle = (int *)realloc(ngbr->jcelle,sizeof(int)*ncell_nowe);
    simparms->mem_bytes += 2*sizeof(int)*(ncell_now+ncell_nowe);
  }

  ddot = hmat[0]*hmat[3] + hmat[1]*hmat[4] + hmat[2]*hmat[5];
  ddot /= (sqrt(hmat[0]*hmat[0]+hmat[1]*hmat[1]+hmat[2]*hmat[2])*
	   sqrt(hmat[3]*hmat[3]+hmat[4]*hmat[4]+hmat[5]*hmat[5]));

  sinangle = sin(acos(ddot));

  ddot = hmat[0]*hmat[6] + hmat[1]*hmat[7] + hmat[2]*hmat[8];
  ddot /= (sqrt(hmat[0]*hmat[0]+hmat[1]*hmat[1]+hmat[2]*hmat[2])*
	   sqrt(hmat[6]*hmat[6]+hmat[7]*hmat[7]+hmat[8]*hmat[8]));

  ddot = sin(acos(ddot));
  if(ddot<sinangle) sinangle = ddot;

  ddot = hmat[3]*hmat[6] + hmat[4]*hmat[7] + hmat[5]*hmat[8];
  ddot /= (sqrt(hmat[3]*hmat[3]+hmat[4]*hmat[4]+hmat[5]*hmat[5])*
	   sqrt(hmat[6]*hmat[6]+hmat[7]*hmat[7]+hmat[8]*hmat[8]));

  ddot = sin(acos(ddot));
  if(ddot<sinangle) sinangle = ddot;

  for(icell=ibegin;icell<iend;icell++){ /* loop over all cells */
    ic = icell/ncell2;               /* get cell i,j,k index */
    jc = (icell - ic*ncell2)/ncell;
    kc = (icell - ic*ncell2 - jc*ncell);
    for(i=0;i<=nmax;i++){            /* get jcell from i,j,k */
      for(j = (i==0)?0:-nmax;j<=nmax;j++){
	for(k= (i==0 && j==0)?0:-nmax;k<=nmax;k++){
	  in = (ic+i)%ncell;
	  jn = (jc+j)%ncell;
	  kn = (kc+k)%ncell;

	  in = abs(ic-in);
	  jn = abs(jc-jn);
	  kn = abs(kc-kn);

	  /* apply boundry conditions */
	  
	  if(iperd>0 && in>ncell/2) in = ncell-in;
	  if(iperd>1 && jn>ncell/2) jn = ncell-jn;
	  if(iperd>2 && kn>ncell/2) kn = ncell-kn;

	  /* subtract one cell length so that particles on edges interact */
	  if(in>0)in--;
	  if(jn>0)jn--;
	  if(kn>0)kn--;

	  /* get the real distance */
	  dx = hmat[0]*in*rcell+hmat[1]*jn*rcell+hmat[2]*kn*rcell;
	  dy = hmat[3]*in*rcell+hmat[4]*jn*rcell+hmat[5]*kn*rcell;
	  dz = hmat[6]*in*rcell+hmat[7]*jn*rcell+hmat[8]*kn*rcell;

	  /* scale by the sine of the angle for shortest distance */
	  r = sqrt(dx*dx+dy*dy+dz*dz)*sinangle;

	  /* if cells interact add to interaction list */
	  if(r<rmax){
	    in = (ic+i);
	    jn = (jc+j);
	    kn = (kc+k);

	    if(in<0)in += ncell;
	    if(jn<0)jn += ncell;
	    if(kn<0)kn += ncell;

	    if(in>=ncell)in -= ncell;
	    if(jn>=ncell)jn -= ncell;
	    if(kn>=ncell)kn -= ncell;

	    ngbr->icell[ngbr->npaircell] = icell;
	    ngbr->jcell[ngbr->npaircell] = in*ncell2+jn*ncell+kn;

	    ngbr->npaircell++;
	  } /*endif*/
	  if(r<rmaxe){
	    in = (ic+i);
	    jn = (jc+j);
	    kn = (kc+k);

	    if(in<0)in += ncell;
	    if(jn<0)jn += ncell;
	    if(kn<0)kn += ncell;

	    if(in>=ncell)in -= ncell;
	    if(jn>=ncell)jn -= ncell;
	    if(kn>=ncell)kn -= ncell;

	    ngbr->icelle[ngbr->npaircelle] = icell;
	    ngbr->jcelle[ngbr->npaircelle] = in*ncell2+jn*ncell+kn;

	    ngbr->npaircelle++;
	  } /*endif*/
	}
      }
    }/*end i loop */
  }


  if(mem_now != simparms->mem_bytes){
    sprintf(line,"Number of interacting cells (procs #%d) is %d (out of %d)",
	    simparms->rank,ngbr->npaircell,(num_cells*(ncells+1)/2));
    md_stdout(line);
    sprintf(line,"Number of electrostatic cells (procs #%d) is %d (out of %d)",
	    simparms->rank,ngbr->npaircelle,(num_cells*(ncells+1)/2));
    md_stdout(line);
    
    if(simparms->rank==0){
      sprintf(line,"Bytes of memory allocated so far = %d",
	      simparms->mem_bytes);
      md_stdout(line);
    }
  }
#ifdef DEBUG
  printf("Cell neighbors\n");
  for(i=0;i<ngbr->npaircell;i++){
    printf("%d %d %d\n",i,ngbr->icell[i],ngbr->jcell[i]);
  }
#endif
}
/*--------------------------------------------------------*/
void get_lnklist(int natoms,int ncell,int *head,int *llist,int iperd,
		 double *hmat,double *x,double *y,double *z)
{
  int i,ic,jc,kc,ncell3,ncell2,ncellm1,icell;
  double sx,sy,sz,hmati[9];

  ncell2 = ncell*ncell;
  ncell3 = ncell*ncell2;
  ncellm1 = ncell-1;

  gethinv9(hmat,hmati);
  for(i=0;i<ncell3;i++){
    head[i] = -1;
  }

  for(i=natoms-1;i>=0;i--){
    sx = x[i]*hmati[0]+y[i]*hmati[1]+z[i]*hmati[2];
    sy = x[i]*hmati[3]+y[i]*hmati[4]+z[i]*hmati[5];
    sz = x[i]*hmati[6]+y[i]*hmati[7]+z[i]*hmati[8];
    if(iperd>0) sx -= anint(sx);
    if(iperd>1) sy -= anint(sy);
    if(iperd>2) sz -= anint(sz);
    ic = (int)((sx+.5)*ncell);
    jc = (int)((sy+.5)*ncell);
    kc = (int)((sz+.5)*ncell);
    ic = MAX(ic,0);
    jc = MAX(jc,0);
    kc = MAX(kc,0);
    ic = MIN(ic,ncellm1);
    jc = MIN(jc,ncellm1);
    kc = MIN(kc,ncellm1);
    icell = ic*ncell2+jc*ncell+kc;
    llist[i] = head[icell];
    head[icell] = i;
  }
#ifdef DEBUG
  printf("Link head\n");
  for(i=0;i<ncell3;i++){
    printf("%d %d\n",i,head[i]);
  }
  printf("Link List\n");

  for(i=0;i<natoms;i++){
    printf("%d %d\n",i,llist[i]);
  }
#endif
}
/*--------------------------------------------------------*/
void get_lnklist_e(int natoms,int ncell,int *head,int *llist,int iperd,
		   double *hmat,double *x,double *y,double *z,double *q)
{
  int i,ic,jc,kc,ncell3,ncell2,ncellm1,icell;
  double sx,sy,sz,hmati[9];

  ncell2 = ncell*ncell;
  ncell3 = ncell*ncell2;
  ncellm1 = ncell-1;

  gethinv9(hmat,hmati);
  for(i=0;i<ncell3;i++){
    head[i] = -1;
  }

  for(i=natoms-1;i>=0;i--){
    if(q[i]!=0.0){
      sx = x[i]*hmati[0]+y[i]*hmati[1]+z[i]*hmati[2];
      sy = x[i]*hmati[3]+y[i]*hmati[4]+z[i]*hmati[5];
      sz = x[i]*hmati[6]+y[i]*hmati[7]+z[i]*hmati[8];

      if(iperd>0) sx -= anint(sx);
      if(iperd>1) sy -= anint(sy);
      if(iperd>2) sz -= anint(sz);
      ic = (int)((sx+.5)*ncell);
      jc = (int)((sy+.5)*ncell);
      kc = (int)((sz+.5)*ncell);
      ic = MAX(ic,0);
      jc = MAX(jc,0);
      kc = MAX(kc,0);
      ic = MIN(ic,ncellm1);
      jc = MIN(jc,ncellm1);
      kc = MIN(kc,ncellm1);
      icell = ic*ncell2+jc*ncell+kc;
      llist[i] = head[icell];
      head[icell] = i;
    }
  }
#ifdef DEBUG
  printf("Link head Elec\n");
  for(i=0;i<ncell3;i++){
    printf("%d %d\n",i,head[i]);
  }
  printf("Link List Elec\n");

  for(i=0;i<natoms;i++){
    printf("%d %d\n",i,llist[i]);
  }
#endif
}
/*-----------------------------------------------------------------------*/
void count_link_pair(int ncell,int iperd,double hmat[],int ibegin,int iend,
		     int *nvdw,int *nelec,double rmax,double rmaxe)
{
  int icell,nmax,ncell2,ic,jc,kc,i,j,k,in,jn,kn;
  double r,dx,dy,dz,rcell,ddot,sinangle;

  *nvdw = *nelec = 0;
  nmax = ncell/2;
  rcell = 1./(double)ncell;
  ncell2 = ncell*ncell;

  ddot = hmat[0]*hmat[3] + hmat[1]*hmat[4] + hmat[2]*hmat[5];
  ddot /= (sqrt(hmat[0]*hmat[0]+hmat[1]*hmat[1]+hmat[2]*hmat[2])*
	   sqrt(hmat[3]*hmat[3]+hmat[4]*hmat[4]+hmat[5]*hmat[5]));

  sinangle = sin(acos(ddot));

  ddot = hmat[0]*hmat[6] + hmat[1]*hmat[7] + hmat[2]*hmat[8];
  ddot /= (sqrt(hmat[0]*hmat[0]+hmat[1]*hmat[1]+hmat[2]*hmat[2])*
	   sqrt(hmat[6]*hmat[6]+hmat[7]*hmat[7]+hmat[8]*hmat[8]));

  ddot = sin(acos(ddot));
  if(ddot<sinangle) sinangle = ddot;

  ddot = hmat[3]*hmat[6] + hmat[4]*hmat[7] + hmat[5]*hmat[8];
  ddot /= (sqrt(hmat[3]*hmat[3]+hmat[4]*hmat[4]+hmat[5]*hmat[5])*
	   sqrt(hmat[6]*hmat[6]+hmat[7]*hmat[7]+hmat[8]*hmat[8]));

  ddot = sin(acos(ddot));
  if(ddot<sinangle) sinangle = ddot;

  for(icell=ibegin;icell<iend;icell++){ /* loop over all cells */
    ic = icell/ncell2;               /* get cell i,j,k index */
    jc = (icell - ic*ncell2)/ncell;
    kc = (icell - ic*ncell2 - jc*ncell);
    for(i=0;i<=nmax;i++){            /* get jcell from i,j,k */
      for(j = (i==0)?0:-nmax;j<=nmax;j++){
	for(k= (i==0 && j==0)?0:-nmax;k<=nmax;k++){
	  in = (ic+i)%ncell;
	  jn = (jc+j)%ncell;
	  kn = (kc+k)%ncell;

	  in = abs(ic-in);
	  jn = abs(jc-jn);
	  kn = abs(kc-kn);

	  /* apply boundry conditions */
	  
	  if(iperd>0 && in>ncell/2) in = ncell-in;
	  if(iperd>1 && jn>ncell/2) jn = ncell-jn;
	  if(iperd>2 && kn>ncell/2) kn = ncell-kn;

	  /* subtract one cell length so that particles 
	     in the corner are included */
	  if(in>0)in--;
	  if(jn>0)jn--;
	  if(kn>0)kn--;

	  /* get the real distance */
	  dx = hmat[0]*in*rcell+hmat[1]*jn*rcell+hmat[2]*kn*rcell;
	  dy = hmat[3]*in*rcell+hmat[4]*jn*rcell+hmat[5]*kn*rcell;
	  dz = hmat[6]*in*rcell+hmat[7]*jn*rcell+hmat[8]*kn*rcell;

	  /* scale by the sine of the angle for shortest distance */
	  r = sqrt(dx*dx+dy*dy+dz*dz)*sinangle;

	  /* if cells interact add to interaction list */
	  if(r<rmax){
	    (*nvdw)++;
	  } /*endif*/
	  if(r<rmaxe){
	    (*nelec)++;
	  } /*endif*/
	}
      }
    }/*end i loop */
  }
}
/*-------------------------------------------------------------------*/
