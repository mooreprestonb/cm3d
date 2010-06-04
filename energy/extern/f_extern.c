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

/*  External electrical potential.... */

#include "md.h"

/* #define VERBOSE  */
/*-----------------------------------------------------------------------*/

void force_extern(SIMPARMS *simparms,COORDS *coords,int ntable)
{
  int i,j,natoms,ierror,jer,k,*ityp,is;
  long ibegin,iend;
  double *xdis,*ydis,*zdis,*dv,*dve,*px,*py,*pz,*qch,**dtab;
  double f,p,fm1,f0,f1,fcc,xp,rdx,zmin,*fz;
  double hmati[9];
  LINE line;

  natoms = simparms->natoms;

  decomp1d((long)natoms,simparms->size,simparms->rank,&ibegin,&iend);
  is = iend-ibegin;
  xdis = (double *)cmalloc(3*is*sizeof(double));
  ydis = xdis +   is;
  zdis = xdis + 2*is;
  gethinv9(coords->hmat,hmati);

  px = coords->px;  py = coords->py;  pz = coords->pz;
  qch = coords->qch;
  dve = coords->dvtab_extern_e;
  rdx = 1./coords->dxtab_extern;
  dtab = coords->dvtab_extern;
  zmin = coords->zmin_extern;
  ityp = simparms->itype;
  fz = coords->fzr;

  for(i=ibegin;i<iend;i++){
    j = i-ibegin;
    xdis[j] = px[i];        ydis[j] = py[i];        zdis[j] = pz[i];
  } 
  period(is,xdis,ydis,zdis,coords->hmat,hmati,simparms->iperd,simparms->ivol);

  jer = ierror = 0;
#include "vectorize.h"
  for(i=ibegin;i<iend;i++){
    dv = dtab[ityp[i]];
    xp = zdis[i-ibegin];
    p = rdx*(xp-zmin);
    k = (int)p;
    if(k<1){ierror=1;jer = i;k=1;p=1.;}
    if(k>ntable-2){ierror=1;jer = i;k=ntable-2;p=(double)k;}
    p = p-(double)k;

    /* external force  */
    fm1 = dv[k-1];    f0  = dv[k];    f1  = dv[k+1];
#ifdef CUBIC
    f2  = dv[k+2];
    f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
    f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
    fcc = -f;

    /* external electrostatic force */
    fm1 = dve[k-1];    f0  = dve[k];    f1  = dve[k+1];
#ifdef CUBIC
    f2  = dve[k+2];
    f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
    f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif

#ifdef VERBOSE
    printf("ext_force %d %g %g %g %d %d\n",i,xp,fcc,f,ityp[i],k);
#endif

    fcc += -qch[i]*f;
    /* add force */
    fz[i] += fcc;
  }

  if(ierror){
    sprintf(line,"distance out of range in extern force routine - ");
    sprintf(line,"%s atom # %d ",line,jer);
    md_warning(line);
  }

#ifdef VERBOSE 
  if(simparms->rank==0){
    printf("External force\n");
    for(i=0;i<natoms;i++){
      printf("%d total force = %g \n",i,fz[i]);
    }
  }
  exit(1);
#endif

  free(xdis);
}

/* #define VERBOSE */
/*------------------------------------------------------------------------*/
void getvextern(SIMPARMS *simparms,COORDS *coords,double *vext,double *wext,
		double *vext_e,double *wext_e,int ntable)
{
  int i,j,natoms,ierror,jer,k,*itype,icount;
  long ibegin,iend;
  double *xdis,*ydis,*zdis,*px,*py,*pz,*v,*ve,*dv,*dve,*qch;
  double xp,f,fe,p,fm1,f0,f1,zmin,rdx;
  double **vtab,**dvtab,vex,vex_e,wex,wex_e;
  double hmati[9];
  LINE line;
#ifdef CUBIC
  double f2;
#endif

  natoms = simparms->natoms;
  itype  = simparms->itype;
  px  = coords->px;  py  = coords->py;  pz  = coords->pz;
  qch = coords->qch;
  ve  = coords->vtab_extern_e;  vtab = coords->vtab_extern;
  dvtab = coords->dvtab_extern;  dve = coords->dvtab_extern_e;
  zmin= coords->zmin_extern;  rdx = 1./coords->dxtab_extern;

  decomp1d((long)natoms,simparms->size,simparms->rank,&ibegin,&iend);
  icount = iend-ibegin;
  xdis = (double *)cmalloc(3*icount*sizeof(double));
  ydis = xdis +   icount;
  zdis = xdis + 2*icount;
  gethinv9(coords->hmat,hmati);

  for(i=ibegin;i<iend;i++){
    j = i-ibegin;
    xdis[j] = px[i];    
    ydis[j] = py[i];    
    zdis[j] = pz[i];
  } 
  period(icount,xdis,ydis,zdis,coords->hmat,hmati,
	 simparms->iperd,simparms->ivol);

  jer = ierror = 0;
  vex = vex_e = wex = wex_e = 0.;
  
  for(i=ibegin;i<iend;i++){
    v  = vtab[itype[i]];
    xp = zdis[i-ibegin];
    p = (xp-zmin)*rdx;
    k = (int)p;
    if(k<1){ierror=1;jer = i;k=1;p=1.0;}
    if(k>ntable-2){ierror=1;jer = i;k=ntable-2;p=(double)k;}
    p = p-(double)k;
    
    /* external potential */
    fm1 = v[k-1];    f0  = v[k];    f1  = v[k+1];
#ifdef CUBIC
    f2  = v[k+2];
    f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
    f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
    vex += f;

    /* external electrstatic potential */
    fm1 = ve[k-1];    f0  = ve[k];    f1  = ve[k+1];
#ifdef CUBIC
    f2  = ve[k+2];
    fe=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
    fe = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
    vex_e += qch[i]*fe;

#ifdef VERBOSE
    printf("ext_pot %d %g %g %g %g %g %d %d\n",i,xp,f,fe,vex,vex_e,itype[i],k);
#endif
    /* pressure calculation external potential */
    dv  = dvtab[itype[i]];
    fm1 = dv[k-1];    f0  = dv[k];    f1  = dv[k+1];
#ifdef CUBIC
    f2  = dv[k+2];
    f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
    f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
    wex -= xp*f;

    /* pressure calculation electrostatic external potential */
    fm1 = dve[k-1];    f0  = dve[k];    f1  = dve[k+1];
#ifdef CUBIC
    f2  = dve[k+2];
    f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
		p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
    f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
    wex_e -= qch[i]*xp*f;
  }
  free(xdis);

#ifdef PARA
  {
    double hmato[4];
    for(i=0;i<4;i++) hmato[i] = 0.0;
    hmati[0] = vex;  hmati[1] = vex_e;  hmati[2] = wex;  hmati[3] = wex_e;
    MPI_Allreduce(hmati,hmato,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    vex = hmato[0];  vex_e = hmato[1];  wex = hmato[2];  wex_e = hmato[3];
  }
#endif

  *vext   += vex;
  *vext_e += vex_e;
  *wext   += wex;
  *wext_e += wex_e;

  if(ierror){
    sprintf(line,"distance out of range in extern force routine - ");
    sprintf(line,"%s atom # %d ",line,jer);
    md_warning(line);
  }
}

