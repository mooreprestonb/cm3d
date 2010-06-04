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

/* subroutines dealing with one_four interactions
   search_14_base - searches the data base and sets up initial vectors
   fonefour - gets the force of the 1-4 interactions
   getvonefour - gets the potential energy of the 1-4 interactions
   onefour   - check to see if to atoms are a special 1-4 interactions
*/

/* #define DEBUG*/

/* #define PRINT_ONEFOUR */
#include "md.h"

static int iatom,jatom,i14type;
static double qiatom,qjatom;
static double **v14tabs,*dx214tabs,*r14max2s,*r14min2s;

/*-----------------------------------------------------------------------*/
int onfo_pair(int iatom,int jatom,ONFO *onfo)
{
  int i,j,itemp,nb;

  if(iatom>jatom){itemp = iatom;iatom = jatom; jatom = itemp;}

  for(nb=0;nb<onfo->n14s;nb++){
    i = onfo->i14[nb]; j = onfo->j14[nb];
    if(i>j){itemp = i; i = j; j = itemp;}
    if(iatom == i && jatom == j) return 1;
  }
  return 0;
}
    
/*-----------------------------------------------------------------------*/
void exclude_onfo(int *nexclude,int **exclude,ONFO *onfo)
{
  int i,j,nb;

  for(nb=0;nb<onfo->n14s;nb++){
    if(onfo->i14[nb]<onfo->j14[nb]){i=onfo->i14[nb];j=onfo->j14[nb];
    } else {j=onfo->i14[nb];i=onfo->j14[nb];}
    insert(exclude[i],&nexclude[i],j);
  }
}
/*-----------------------------------------------------------------------*/
void fonfo(SIMPARMS *simparms,COORDS *coords)
{
  int k14,i,j,k,ii,ierror,jer=0;
  long ibegin,iend;
  double qij;
  double f,p,fm1,f0,f1;
  double xdis,ydis,zdis,dis2,fcc;
  double *px,*py,*pz,*fx,*fy,*fz;
#ifdef CUBIC
  double f2;
#endif
  int n14s,*i14,*j14,*i14dx;
  double **dv14tab,*r14min2,*r14max2,*dx214tab;
  double scale_onefour,scale_onefour_e;
  LINE line;

  n14s = coords->onfo.n14s;
  i14 = coords->onfo.i14;
  j14 = coords->onfo.j14;
  i14dx = coords->onfo.i14dx;
  dv14tab = coords->onfo.dv14tab;
  r14min2 = coords->onfo.r14min2;
  r14max2 = coords->onfo.r14max2;
  dx214tab = coords->onfo.dx214tab;
  scale_onefour = coords->onfo.scale_onefour;
  scale_onefour_e = coords->onfo.scale_onefour_e;

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
  ierror = 0;

  decomp1d((long) n14s,simparms->size,simparms->rank,&ibegin,&iend);
  for(k14=ibegin;k14<iend;k14++){
    i = i14[k14];
    j = j14[k14];
    ii = i14dx[k14];

    xdis = px[i] - px[j];
    ydis = py[i] - py[j];
    zdis = pz[i] - pz[j];

    dis2=xdis*xdis+ydis*ydis+zdis*zdis;

    if(dis2 < r14max2[ii]){
      p = (dis2-r14min2[ii])*dx214tab[ii];
      k = (int)p;
      if(k<1){ierror=1;jer = k14;k=1;p=1.;}
      p = p-(double)k;
      
      fm1 = dv14tab[ii][k-1];
      f0  = dv14tab[ii][k];
      f1  = dv14tab[ii][k+1];
      
#ifdef CUBIC
      f2  = dv14tab[ii][k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      fcc = f*scale_onefour;

      if((qij = coords->qch[i]*coords->qch[j]) != 0){
        fcc -= qij*scale_onefour_e/(sqrt(dis2)*dis2);
      }
      
      xdis  *= fcc;  ydis  *= fcc;  zdis  *= fcc;
      fx[i] -= xdis; fy[i] -= ydis; fz[i] -= zdis;
      fx[j] += xdis; fy[j] += ydis; fz[j] += zdis;
    }
  }

  if(ierror){
    sprintf(line,"distance out of range in onefour force routine\n");
    sprintf(line,"%s between atoms %d and %d",line,i14[jer],j14[jer]);
    md_warning(line);
  }
}
/*------------------------------------------------------------------*/
void getvonfo(SIMPARMS *simparms,COORDS *coords,
	      double *vonfo,double *vonfo_e,double *wonfo,
	      double *wonfo_e,double wonfotensor[9],double wonfo_etensor[9])
{
  int k14,i,j,k,ii,ierror,jer;
  long ibegin,iend;
  double qij;
  double f,p,fm1,f0,f1;
  double xdis,ydis,zdis,dis2;
  double spot,sprs,swelec,spelec;
  double sprstensor[9],swelectensor[9];
#ifdef PARA
  double slocal[22],rlocal[22];
#endif
#ifdef CUBIC
  double f2;
#endif
  int n14s,*i14,*j14,*i14dx;
  double **v14tab,**dv14tab;
  double *r14min2,*r14max2,*dx214tab;
  double scale_onefour,scale_onefour_e;
  LINE line;

  n14s = coords->onfo.n14s;
  i14 = coords->onfo.i14;
  j14 = coords->onfo.j14;
  i14dx = coords->onfo.i14dx;
  v14tab = coords->onfo.v14tab;
  dv14tab = coords->onfo.dv14tab;
  r14min2 = coords->onfo.r14min2;
  r14max2 = coords->onfo.r14max2;
  dx214tab = coords->onfo.dx214tab;
  scale_onefour = coords->onfo.scale_onefour;
  scale_onefour_e = coords->onfo.scale_onefour_e;
  
  spot = sprs = swelec = spelec = 0.;
  for(i=0;i<9;i++) sprstensor[i] = swelectensor[i]=0.0;
  jer = ierror = 0;

  decomp1d((long)n14s,simparms->size,simparms->rank,&ibegin,&iend);
  for(k14=ibegin;k14<iend;k14++){
    i = i14[k14];    j = j14[k14];    ii = i14dx[k14];
    
    xdis = coords->px[i] - coords->px[j];
    ydis = coords->py[i] - coords->py[j];
    zdis = coords->pz[i] - coords->pz[j];

    dis2=xdis*xdis+ydis*ydis+zdis*zdis;

    if(dis2 < r14max2[ii]){
      p = (dis2-r14min2[ii])*dx214tab[ii];
      k = (int)p;
      if(k<1) {ierror= 1; jer = k14;k = 1; p=1.;}

      p = p-(double)k;
      
      fm1 = v14tab[ii][k-1];
      f0  = v14tab[ii][k];
      f1  = v14tab[ii][k+1];
      
#ifdef CUBIC
      f2  = v14tab[ii][k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      spot += f;
      
#ifdef PRINT_ONEFOUR
      printf("onefour %d-%d distance = %g, energy = %g, electrostatic = ",
         i,j,sqrt(dis2),f/KCAL);
#endif

      fm1 = dv14tab[ii][k-1];
      f0  = dv14tab[ii][k];
      f1  = dv14tab[ii][k+1];

#ifdef CUBIC
      f2  = dv14tab[ii][k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      sprs -= dis2*f;
      sprstensor[0] -= xdis*xdis*f;
      sprstensor[1] -= ydis*xdis*f;
      sprstensor[2] -= zdis*xdis*f;

      sprstensor[3] -= xdis*ydis*f;
      sprstensor[4] -= ydis*ydis*f;
      sprstensor[5] -= zdis*ydis*f;

      sprstensor[6] -= xdis*zdis*f;
      sprstensor[7] -= ydis*zdis*f;
      sprstensor[8] -= zdis*zdis*f;

      if((qij = coords->qch[i]*coords->qch[j]) != 0){
        f = qij/sqrt(dis2);
#ifdef PRINT_ONEFOUR
        printf("%g\n",f/KCAL);
#endif

        spelec += f;
        swelec += f;
        
        swelectensor[0] += f*xdis*xdis/dis2;
        swelectensor[1] += f*ydis*xdis/dis2;
        swelectensor[2] += f*zdis*xdis/dis2;

        swelectensor[3] += f*xdis*ydis/dis2;
        swelectensor[4] += f*ydis*ydis/dis2;
        swelectensor[5] += f*zdis*ydis/dis2;

        swelectensor[6] += f*xdis*zdis/dis2;
        swelectensor[7] += f*ydis*zdis/dis2;
        swelectensor[8] += f*zdis*zdis/dis2;

      }
    }
  }

  if(ierror){
    sprintf(line,"distance out of range in onefour energy routine\n");
    sprintf(line,"%s between atoms %d and %d",line,i14[jer],j14[jer]);
    md_warning(line);
  }

#ifdef PARA
  slocal[0] = spot;
  slocal[1] = sprs;
  slocal[2] = spelec;
  slocal[3] = swelec;
  for(i=0;i<9;i++) slocal[4+i] = sprstensor[i];
  for(i=0;i<9;i++) slocal[13+i] = swelectensor[i];
  MPI_Allreduce(slocal,rlocal,22,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  spot = rlocal[0];
  sprs = rlocal[1];
  spelec = rlocal[2];
  swelec = rlocal[3];
  for(i=0;i<9;i++) sprstensor[i] = rlocal[4+i];
  for(i=0;i<9;i++) swelectensor[i] = rlocal[13+i];
#endif

  *vonfo   += spot*scale_onefour;
  *wonfo   += sprs*scale_onefour;
  *vonfo_e += spelec*scale_onefour_e;
  *wonfo_e += swelec*scale_onefour_e;

  for(i=0;i<9;i++) {
   wonfotensor[i] += sprstensor[i]*scale_onefour;
   wonfo_etensor[i] += swelectensor[i]*scale_onefour_e;
  }
}
/*------------------------------------------------------------------*/

void fconefour(COORDS *coords,double **fcmat)
{
  int i,j,it,ii,jj,k,k14;
  double ax,ay,az,r,r2,r3;
  double dxx,dyy,dzz,dxy,dxz,dyz,dudr,du2dr2;
  double f,f0,fm1,f1,p;
#ifdef CUBIC
  double f2;
#endif
  int n14s,*i14,*j14,*i14dx;
  double **dv14tab,**d2v14tab;
  double *r14min2,*r14max2,*dx214tab;
  double scale_onefour,scale_onefour_e;

  n14s = coords->onfo.n14s;
  i14 = coords->onfo.i14;
  j14 = coords->onfo.j14;
  i14dx = coords->onfo.i14dx;
  dv14tab = coords->onfo.dv14tab;
  d2v14tab = coords->onfo.d2v14tab;
  r14min2 = coords->onfo.r14min2;
  r14max2 = coords->onfo.r14max2;
  dx214tab = coords->onfo.dx214tab;
  scale_onefour = coords->onfo.scale_onefour;
  scale_onefour_e = coords->onfo.scale_onefour_e;
  
  for(k14=0;k14<n14s;k14++){
    i = i14[k14];    j = j14[k14];    it = i14dx[k14];
    
    ax = coords->px[i] - coords->px[j];
    ay = coords->py[i] - coords->py[j];
    az = coords->pz[i] - coords->pz[j];

    r2 = (ax*ax+ay*ay+az*az);
    
#ifdef DEBUG
    printf("%d %d %d %g\n",i,j,it,r2); 
#endif
    if(r2<r14min2[it]){
      fprintf(stderr,
	      "ERROR: distance is shorter then short range cutoff(in 14)\n");
      fprintf(stderr,"r2 = %g < %g = r14min2\n",r2,r14min2[it]);
      fprintf(stderr,"between atoms %d and %d (start count from 0)\n",i,j);
      fprintf(stderr,"%d %g %g %g\n",i,
	      coords->px[i],coords->py[i],coords->pz[i]);
      fprintf(stderr,"%d %g %g %g\n",j,
	      coords->px[j],coords->py[j],coords->pz[j]);
      exit(1);
    }
    
    if(r2<r14max2[it]){
      p = (r2-r14min2[it])*dx214tab[it];
      k = (int)p;
#ifdef DEBUG
      if(k<0 || k > coords->onfo.n14table-4){
	fprintf(stderr,"ERROR distance out of range %g %d\n",r2,k);
	exit(1);
      }
#endif
      p = p-(double)k;
      
      fm1 = dv14tab[it][k-1];
      f0  = dv14tab[it][k];
      f1  = dv14tab[it][k+1];
      
#ifdef CUBIC
      f2  = dv14tab[it][k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      r2 = 1./r2;
      r = sqrt(r2);
      r3 = r2*r;

      /* r*f because dv14tab is NOT dv/dr but (1/r)dv/dr */
      dudr = (scale_onefour*r*f-
	      (coords->qch[i]*coords->qch[j]*r2)*scale_onefour_e);
      
      fm1 = d2v14tab[it][k-1];
      f0  = d2v14tab[it][k];
      f1  = d2v14tab[it][k+1];
    
#ifdef CUBIC
      f2  = d2v14tab[it][k+2];
      f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		  p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
      f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
      du2dr2 = (scale_onefour*f + 
		scale_onefour_e*2.*coords->qch[i]*coords->qch[j]*r3);
      
      dxx = du2dr2*ax*ax*r2+dudr*(-ax*ax*r3 + r);
      dyy = du2dr2*ay*ay*r2+dudr*(-ay*ay*r3 + r);
      dzz = du2dr2*az*az*r2+dudr*(-az*az*r3 + r);
      dxy = du2dr2*ax*ay*r2+dudr*(-ax*ay*r3);
      dxz = du2dr2*ax*az*r2+dudr*(-ax*az*r3);
      dyz = du2dr2*ay*az*r2+dudr*(-ay*az*r3);
      
      ii=3*i-1;
      jj=3*j-1;
    
      /* fill in force matrix  */
      fcmat[ii+1][jj+1] += -dxx;
      fcmat[ii+1][jj+2] += -dxy;
      fcmat[ii+1][jj+3] += -dxz;
      fcmat[ii+2][jj+1] += -dxy;
      fcmat[ii+2][jj+2] += -dyy;
      fcmat[ii+2][jj+3] += -dyz;
      fcmat[ii+3][jj+1] += -dxz;
      fcmat[ii+3][jj+2] += -dyz;
      fcmat[ii+3][jj+3] += -dzz;
      
      /* fill in diagonal elements of matrix  */
      
      fcmat[ii+1][ii+1] += dxx;
      fcmat[ii+1][ii+2] += dxy;
      fcmat[ii+1][ii+3] += dxz;
      fcmat[ii+2][ii+1] += dxy;
      fcmat[ii+2][ii+2] += dyy;
      fcmat[ii+2][ii+3] += dyz;
      fcmat[ii+3][ii+1] += dxz;
      fcmat[ii+3][ii+2] += dyz;
      fcmat[ii+3][ii+3] += dzz;
      
      fcmat[jj+1][jj+1] += dxx;
      fcmat[jj+1][jj+2] += dxy;
      fcmat[jj+1][jj+3] += dxz;
      fcmat[jj+2][jj+1] += dxy;
      fcmat[jj+2][jj+2] += dyy;
      fcmat[jj+2][jj+3] += dyz;
      fcmat[jj+3][jj+1] += dxz;
      fcmat[jj+3][jj+2] += dyz;
      fcmat[jj+3][jj+3] += dzz;
    }
  }
}
/*-------------------------------------------------------------*/
void fconefourn(COORDS *coords,double **fcmat)
{
  int k14;

  r14max2s = coords->onfo.r14max2;
  r14min2s = coords->onfo.r14min2;
  v14tabs  = coords->onfo.v14tab;
  dx214tabs  = coords->onfo.dx214tab;
  for(k14=0;k14<coords->onfo.n14s;k14++){
    iatom = coords->onfo.i14[k14];
    jatom = coords->onfo.j14[k14];
    i14type = coords->onfo.i14dx[k14];
    qiatom = coords->qch[iatom];
    qjatom = coords->qch[jatom];
    ngradv2(coords->px,coords->py,coords->pz,iatom,iatom,pot_onefour,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,iatom,jatom,pot_onefour,fcmat);
    ngradv2(coords->px,coords->py,coords->pz,jatom,jatom,pot_onefour,fcmat);
  }
}
/*-------------------------------------------------------------*/
double pot_onefour(double *px,double *py,double *pz)
{
  int k;
  double ax,ay,az,r2,spot,p;
  double f,f0,fm1,f1;
#ifdef CUBIC
  double f2;
#endif

  spot = 0.;
  
  ax = px[iatom] - px[jatom];
  ay = py[iatom] - py[jatom];
  az = pz[iatom] - pz[jatom];

  r2 = (ax*ax+ay*ay+az*az);
  
  if(r2<r14min2s[i14type]){
    fprintf(stderr,
	    "ERROR: distance is shorter then short range cutoff (in 14)\n");
    fprintf(stderr,"r2 = %g < %g = r14min2\n",r2,r14min2s[i14type]);
    exit(1);
  }
  
  if(r2<r14max2s[i14type]){
    p = (r2-r14min2s[i14type])*dx214tabs[i14type];
    k = (int)p;
#ifdef DEBUG
    if(k<1){
      fprintf(stderr,"ERROR distance out of range %g %d\n",r2,k);
      exit(1);
    }
#endif
    p = p-(double)k;
    
    fm1 = v14tabs[i14type][k-1];
    f0  = v14tabs[i14type][k];
    f1  = v14tabs[i14type][k+1];
    
#ifdef CUBIC
    f2  = v14tabs[i14type][k+2];
    f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
    f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
    spot = f+qiatom*qjatom/sqrt(r2);
  }
  
  return spot;
}
/*-------------------------------------------------------------*/
