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


/* #define DEBUG */
#define H (1.e-7)

#ifdef DEBUG
#include <stdio.h>
#endif

/* calculate second derivative numericaly */

void ngradv2(double *px,double *py,double *pz,int i,int j,
	     double (*pot)(double *,double *,double *),double **fcmat)
{
  int ii,jj;
  double vpp,vpm,vmp,vmm,oldi,oldj,gradv;

  ii = i*3;
  jj = j*3;

  /* diagonal elements (also = -sum off diag elements for i==j*/ 
#ifdef DEBUG
  printf("potential = %.15g\n",(*pot)(px,py,pz));
#endif

  if(i == j){
    vmp = (*pot)(px,py,pz);

    oldi = px[i];
    px[i] = oldi + H;
    vpp = (*pot)(px,py,pz);
    px[i] = oldi - H;
    vmm = (*pot)(px,py,pz);
    px[i] = oldi;
    gradv = (vpp+vmm-2.*vmp)/(H*H);
    fcmat[ii][ii] += gradv;

    oldi = py[i];
    py[i] = oldi + H;
    vpp = (*pot)(px,py,pz);
    py[i] = oldi - H;
    vmm = (*pot)(px,py,pz);
    py[i] = oldi;
    gradv = (vpp+vmm-2.*vmp)/(H*H);
    fcmat[ii+1][ii+1] += gradv;

    oldi = pz[i];
    pz[i] = oldi + H;
    vpp = (*pot)(px,py,pz);
    pz[i] = oldi - H;
    vmm = (*pot)(px,py,pz);
    pz[i] = oldi;
    gradv = (vpp+vmm-2.*vmp)/(H*H);
    fcmat[ii+2][ii+2] += gradv;

  } else {

    /* X-X */
    oldi = px[i];
    oldj = px[j];
    px[i] = oldi + H;
    px[j] = oldj + H;
    vpp = (*pot)(px,py,pz);
    px[i] = oldi - H;
    px[j] = oldj + H;
    vmp = (*pot)(px,py,pz);
    px[i] = oldi + H;
    px[j] = oldj - H;
    vpm = (*pot)(px,py,pz);
    px[i] = oldi - H;
    px[j] = oldj - H;
    vmm = (*pot)(px,py,pz);
    px[i] = oldi;
    px[j] = oldj;
    gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
    fcmat[ii][jj] += gradv;

    /* Y-Y */
    oldi = py[i];
    oldj = py[j];
    py[i] = oldi + H;
    py[j] = oldj + H;
    vpp = (*pot)(px,py,pz);
    py[i] = oldi - H;
    py[j] = oldj + H;
    vmp = (*pot)(px,py,pz);
    py[i] = oldi + H;
    py[j] = oldj - H;
    vpm = (*pot)(px,py,pz);
    py[i] = oldi - H;
    py[j] = oldj - H;
    vmm = (*pot)(px,py,pz);
    py[i] = oldi;
    py[j] = oldj;
    gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
    fcmat[ii+1][jj+1] += gradv;
    
    /* Z-Z */
    oldi = pz[i];
    oldj = pz[j];
    pz[i] = oldi + H;
    pz[j] = oldj + H;
    vpp = (*pot)(px,py,pz);
    pz[i] = oldi - H;
    pz[j] = oldj + H;
    vmp = (*pot)(px,py,pz);
    pz[i] = oldi + H;
    pz[j] = oldj - H;
    vpm = (*pot)(px,py,pz);
    pz[i] = oldi - H;
    pz[j] = oldj - H;
    vmm = (*pot)(px,py,pz);
    pz[i] = oldi;
    pz[j] = oldj;
    gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
    fcmat[ii+2][jj+2] += gradv;
  }
  
  /* X-Y */
  oldi = px[i];
  oldj = py[j];
  px[i] = oldi + H;
  py[j] = oldj + H;
  vpp = (*pot)(px,py,pz);
  px[i] = oldi - H;
  py[j] = oldj + H;
  vmp = (*pot)(px,py,pz);
  px[i] = oldi + H;
  py[j] = oldj - H;
  vpm = (*pot)(px,py,pz);
  px[i] = oldi - H;
  py[j] = oldj - H;
  vmm = (*pot)(px,py,pz);
  px[i] = oldi;
  py[j] = oldj;
  gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
  fcmat[ii][jj+1] += gradv;
  
  /* X-Z */
  oldi = px[i];
  oldj = pz[j];
  px[i] = oldi + H;
  pz[j] = oldj + H;
  vpp = (*pot)(px,py,pz);
  px[i] = oldi - H;
  pz[j] = oldj + H;
  vmp = (*pot)(px,py,pz);
  px[i] = oldi + H;
  pz[j] = oldj - H;
  vpm = (*pot)(px,py,pz);
  px[i] = oldi - H;
  pz[j] = oldj - H;
  vmm = (*pot)(px,py,pz);
  px[i] = oldi;
  pz[j] = oldj;
  gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
  fcmat[ii][jj+2] += gradv;
  
  /* Y-X */
  oldi = py[i];
  oldj = px[j];
  py[i] = oldi + H;
  px[j] = oldj + H;
  vpp = (*pot)(px,py,pz);
  py[i] = oldi - H;
  px[j] = oldj + H;
  vmp = (*pot)(px,py,pz);
  py[i] = oldi + H;
  px[j] = oldj - H;
  vpm = (*pot)(px,py,pz);
  py[i] = oldi - H;
  px[j] = oldj - H;
  vmm = (*pot)(px,py,pz);
  py[i] = oldi;
  px[j] = oldj;
  gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
  fcmat[ii+1][jj] += gradv;
  
  /* Y-Z */
  oldi = py[i];
  oldj = pz[j];
  py[i] = oldi + H;
  pz[j] = oldj + H;
  vpp = (*pot)(px,py,pz);
  py[i] = oldi - H;
  pz[j] = oldj + H;
  vmp = (*pot)(px,py,pz);
  py[i] = oldi + H;
  pz[j] = oldj - H;
  vpm = (*pot)(px,py,pz);
  py[i] = oldi - H;
  pz[j] = oldj - H;
  vmm = (*pot)(px,py,pz);
  py[i] = oldi;
  pz[j] = oldj;
  gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
  fcmat[ii+1][jj+2] += gradv;
    
  /* Z-X */
  oldi = pz[i];
  oldj = px[j];
  pz[i] = oldi + H;
  px[j] = oldj + H;
  vpp = (*pot)(px,py,pz);
  pz[i] = oldi - H;
  px[j] = oldj + H;
  vmp = (*pot)(px,py,pz);
  pz[i] = oldi + H;
  px[j] = oldj - H;
  vpm = (*pot)(px,py,pz);
  pz[i] = oldi - H;
  px[j] = oldj - H;
  vmm = (*pot)(px,py,pz);
  pz[i] = oldi;
  px[j] = oldj;
  gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
  fcmat[ii+2][jj] += gradv;
  
  /* Z-Y */
  oldi = pz[i];
  oldj = py[j];
  pz[i] = oldi + H;
  py[j] = oldj + H;
  vpp = (*pot)(px,py,pz);
  pz[i] = oldi - H;
  py[j] = oldj + H;
  vmp = (*pot)(px,py,pz);
  pz[i] = oldi + H;
  py[j] = oldj - H;
  vpm = (*pot)(px,py,pz);
  pz[i] = oldi - H;
  py[j] = oldj - H;
  vmm = (*pot)(px,py,pz);
  pz[i] = oldi;
  py[j] = oldj;
  gradv = (vpp+vmm-vpm-vmp)/(4.*H*H);
  fcmat[ii+2][jj+1] += gradv;
}
