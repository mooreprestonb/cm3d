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


/* routine to calculate force constant matrix of 
   nonbonded interactions */


#include "md.h"

/* #define DEBUG */
#define CUBIC

/* global variables for numerical force constant matrix evaluation */
int iatom,jatom,iperdn,iensembln;
double *hmatn,*hmatin,qiatom,qjatom;
INTER *intern;

/*----------------------------------------------------------------*/
void fcinter(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **fcmat)
{
  int i,j,it,jt,ii,jj,intr,*excl,nexcl,k;
  int *map,*itype;
  double ax,ay,az,r,r2,r3;
  double dxx,dyy,dzz,dxy,dxz,dyz,dudr,du2dr2;
  double f,f0,fm1,f1,p;
  double qi,qj,*q,**dv,**d2v,*dx2;
  double *x,*y,*z,*rmax,*rmin;

#ifdef CUBIC
  double f2;
#endif

  /* derivatives of vdw + coulomb cross terms */
  x = coords->px;
  y = coords->py;
  z = coords->pz;
  q = coords->qch;
  itype = inter->itype;
  rmax = inter->rmaxl_cut2;
  rmin = inter->rmins_cut2;
  dx2  = inter->dx2tab;
  dv   = inter->dvtab;
  d2v  = inter->d2vtab;

  for(i=0;i<simparms->natoms-1;i++){
    it = itype[i];
    excl = &(inter->exclude[i][0]);
    nexcl = *excl;
    qi = q[i];
    map = inter->map[it];
    for(j=i+1;j<simparms->natoms;j++){
      if (nexcl == j){ 
        nexcl = *(++excl);
      } else if (map[itype[j]] == -1){
        /* skip null interactions (check if also bonded exclusion)*/
      } else {
        qj = q[j];
        jt = itype[j];
        intr = map[jt];
        
        ax = x[i] - x[j];
        ay = y[i] - y[j];
        az = z[i] - z[j];
        
        period(1,&ax,&ay,&az,coords->hmat,coords->hmati,
           simparms->iperd,simparms->ivol);
        
        r2 = (ax*ax+ay*ay+az*az);

        if(r2<rmax[intr]){
          p = (r2-rmin[intr])*dx2[intr];
          k = (int)p;
          
          r2 = 1./r2;
          r = sqrt(r2);
          r3 = r2*r;
          
#ifdef DEBUG
          if(k<1 || k > inter->ntable-4){
            fprintf(stderr,"ERROR distance out of range %g %d\n",r2,k);
            exit(1);
          }
#endif
          p = p-(double)k;
          
          fm1 = dv[intr][k-1];  f0  = dv[intr][k];  f1  = dv[intr][k+1];
#ifdef CUBIC
          f2  = dv[intr][k+2];
          f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
             p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
          f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
          /* r*f because dvtab is NOT dv/dr but (1/r)dv/dr */
          dudr = f/r-qi*qj*r2;
          
          fm1 = d2v[intr][k-1]; f0 = d2v[intr][k]; f1 = d2v[intr][k+1];
#ifdef CUBIC
          f2  = inter->d2vtab[intr][k+2];
          f=f0+(p/2)*((2*f1-2*fm1/3-f0-f2/3)+
             p*((fm1-2*f0+f1)+p*(f0-fm1/3-f1+f2/3)));
#else
          f = f0 + .5*p*(f1-fm1 + p*(fm1-2.*f0+f1));
#endif
          du2dr2 = f+2.*qi*qj*r3;
          
          dxx = du2dr2*ax*ax*r2+dudr*(-ax*ax*r3 + r);
          dyy = du2dr2*ay*ay*r2+dudr*(-ay*ay*r3 + r);
          dzz = du2dr2*az*az*r2+dudr*(-az*az*r3 + r);
          dxy = du2dr2*ax*ay*r2+dudr*(-ax*ay*r3);
          dxz = du2dr2*ax*az*r2+dudr*(-ax*az*r3);
          dyz = du2dr2*ay*az*r2+dudr*(-ay*az*r3);
          
          ii=3*i;
          jj=3*j;

          /* fill in force matrix  */
          fcmat[ii  ][jj  ] += -dxx;
          fcmat[ii  ][jj+1] += -dxy;
          fcmat[ii  ][jj+2] += -dxz;
          fcmat[ii+1][jj  ] += -dxy;
          fcmat[ii+1][jj+1] += -dyy;
          fcmat[ii+1][jj+2] += -dyz;
          fcmat[ii+2][jj  ] += -dxz;
          fcmat[ii+2][jj+1] += -dyz;
          fcmat[ii+2][jj+2] += -dzz;
	  
          /* fill in diagonal elements of matrix  */
          
          fcmat[ii  ][ii  ] += dxx;
          fcmat[ii  ][ii+1] += dxy;
          fcmat[ii  ][ii+2] += dxz;
          fcmat[ii+1][ii  ] += dxy;
          fcmat[ii+1][ii+1] += dyy;
          fcmat[ii+1][ii+2] += dyz;
          fcmat[ii+2][ii  ] += dxz;
          fcmat[ii+2][ii+1] += dyz;
          fcmat[ii+2][ii+2] += dzz;
	  
          fcmat[jj  ][jj  ] += dxx;
          fcmat[jj  ][jj+1] += dxy;
          fcmat[jj  ][jj+2] += dxz;
          fcmat[jj+1][jj  ] += dxy;
          fcmat[jj+1][jj+1] += dyy;
          fcmat[jj+1][jj+2] += dyz;
          fcmat[jj+2][jj  ] += dxz;
          fcmat[jj+2][jj+1] += dyz;
          fcmat[jj+2][jj+2] += dzz;
        }
      }
    }
  }
}
/*-------------------------------------------------------------*/
void fcintern(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **fcmat)
{
  int *excl,nexcl;

  intern = inter;
  iperdn = simparms->iperd;
  iensembln  = simparms->ivol;
  hmatn = coords->hmat;
  hmatin= coords->hmati;

  for(iatom=0;iatom<simparms->natoms-1;iatom++){
    qiatom = coords->qch[iatom];
    excl = &(inter->exclude[iatom][0]);
    nexcl = *(excl);
    for(jatom=iatom+1;jatom<simparms->natoms;jatom++){
      if (nexcl == jatom){ 
        nexcl = *(++excl);
      } else if (inter->map[iatom][inter->itype[jatom]] == -1){
        /* skip null interactions (check if also bonded exclusion)*/
      } else {
        qjatom = coords->qch[jatom];
        ngradv2(coords->px,coords->py,coords->pz,iatom,iatom,pot_inter,fcmat);
        ngradv2(coords->px,coords->py,coords->pz,iatom,jatom,pot_inter,fcmat);
        ngradv2(coords->px,coords->py,coords->pz,jatom,jatom,pot_inter,fcmat);
      }
    }
  }
}
/*-------------------------------------------------------------*/
double pot_inter(double *px,double *py,double *pz)
{
  int it,jt,k,intr;
  double ax,ay,az,r2,spot,p;
  double f,f0,fm1,f1;
#ifdef CUBIC
  double f2;
#endif

  spot = 0.;
  it = intern->itype[iatom];
  jt = intern->itype[jatom];
  intr = intern->map[it][jt];

  ax = px[iatom] - px[jatom];
  ay = py[iatom] - py[jatom];
  az = pz[iatom] - pz[jatom];

  period(1,&ax,&ay,&az,hmatn,hmatin,iperdn,iensembln);

  r2 = (ax*ax+ay*ay+az*az);

  if(r2<intern->rmaxl_cut2[intr]){
    if(r2<intern->rmins_cut2[intr]){
      fprintf(stderr,
	      "ERROR: distance is shorter then short range cutoff\n");
      fprintf(stderr,"r2 = %g < %g = rmins2\n",r2,intern->rmins_cut2[intr]);
      exit(1);
    }

    p = (r2-intern->rmins_cut2[intr])*intern->dx2tab[intr];
    k = (int)p;
#ifdef DEBUG
    if(k<1 || k > intern->ntable-4){
      fprintf(stderr,"ERROR distance out of range %g %d\n",r2,k);
      exit(1);
    }
#endif
    p = p-(double)k;

    fm1 = intern->vtab[intr][k-1];
    f0  = intern->vtab[intr][k];
    f1  = intern->vtab[intr][k+1];
    
#ifdef CUBIC
    f2  = intern->vtab[intr][k+2];
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


