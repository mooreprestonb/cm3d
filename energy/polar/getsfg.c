/*   Copyright (C) 1997 1998 1999 2000 Dr. Preston B. Moore

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
#define WATER
/* #define MUTEST */
/* routines to sum polarizability tensor and dipole along interfaces */

/*-------------------------------------------------------------------*/
void getalpha_interface(SIMPARMS *simparms,COORDS *coords,INTER *inter,
   NMODES *nmodes,double **akeep,int imode)
{
  int i,j,k,ii,ihalf,kcount1,kcount2,natoms;
  double **alpha1,**alpha2;

  alpha1 = dmatrix(0,2,0,2);
  alpha2 = dmatrix(0,2,0,2);

  natoms = simparms->natoms;
  
  printf("WARNING this needs to be changed for anything but water (3 atom molecules) \n");
  printf("WARNING from getalpha_interface \n");

  ihalf = 0;
  /* zero variables */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      alpha1[i][j] = alpha2[i][j] = 0.0;
  kcount1 = kcount2 = 0;
  for(i=0;i<natoms;i++){
    ii = i*3;
    /* first atom in molecule determines if in top or bottom half */
    /* only coded for 3 atom molecules !!!!! */
    
    if (!(i%3)) {
      if(coords->pz[i] > coords->cmz){
	ihalf = 1;
      }else{
        ihalf = 2;
      }
    }
    /* extract 3x3 block = mol. polar. tensor &eqn 2.13 N&P */
    if(ihalf == 1){
      kcount1 += 1;
      for(k=0;k<3;k++)
        for(j=0;j<3;j++)
          alpha1[k][j] += akeep[ii+k][j];
    }else if(ihalf == 2) {
      kcount2 += 1;
      for(k=0;k<3;k++)
        for(j=0;j<3;j++)
          alpha2[k][j] += akeep[ii+k][j];
    } else {
      printf("Warning: atom %d not assigned to top or bottom half\n",i);
      exit(0);
    }
  }

  
  /* saving derivatives for inm analysis */
  if(kcount1 !=0)nmodes->A_xx1[imode] = alpha1[0][0]/(double)kcount1;
  if(kcount1 !=0)nmodes->A_zz1[imode] = alpha1[2][2]/(double)kcount1;
  if(kcount1 !=0)nmodes->A_xz1[imode] = alpha1[0][2]/(double)kcount1;
  if(kcount1 !=0)nmodes->A_zx1[imode] = alpha1[2][0]/(double)kcount1;
  if(kcount2 !=0)nmodes->A_xx2[imode] = alpha2[0][0]/(double)kcount2;
  if(kcount2 !=0)nmodes->A_zz2[imode] = alpha2[2][2]/(double)kcount2;
  if(kcount2 !=0)nmodes->A_xz2[imode] = alpha2[0][2]/(double)kcount2;
  if(kcount2 !=0)nmodes->A_zx2[imode] = alpha2[2][0]/(double)kcount2;

  free_dmatrix(alpha1,0,2,0,2);
  free_dmatrix(alpha2,0,2,0,2);
}
/*-------------------------------------------------------------------*/
void getmu_interface(SIMPARMS *simparms,COORDS *coords,INTER *inter,
   NMODES *nmodes,double **dipder)
{
  int imode,i,j,k,ii;
  int ihalf,kcount1,kcount2,natoms;
  double qm,dx1,dy1,dz1,dx2,dy2,dz2;
  double **mat;
#ifdef WATER
  int ll,iii;
  double qmp,sum;
  double dqmpo,dqmph1,dqmph2;
  double masfac;
  double *dqo,*dqh1,*dqh2;
  
  ihalf = 0;
  qmp = sum=dqmpo=dqmph1=dqmph2=masfac = 0.0;
  dqo= malloc(9*sizeof(double));
  dqh1= malloc(9*sizeof(double));
  dqh2= malloc(9*sizeof(double));
#endif

  natoms = simparms->natoms;
  printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  printf("WARNING getmu_interface using Raw coordinates (not in simulation box) in dipole calculation \n");
  printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  
    /* MUTEST sets the fcmat to identity and masses to 1 to get dmu/dx_i instead dmu/dQ_i */
#ifdef MUTEST
  for(j=0;j<natoms;j++){
          coords->amass[j]=1. ;
  }
  
  for(j=0;j<natoms*3;j++){
    for(k=0;k<natoms*3;k++){
        nmodes->fcmat[j][k]=0.  ;
    } 
  }
        
  for(j=0;j<natoms*3;j++){
    nmodes->fcmat[j][j]=1.  ;
  } 
  
  printf("$$$$$$$Temp fmat identity rdamp$$$$$$$$$$$$$$$\n");
  printmatrix(stdout,natoms*3,nmodes->fcmat);
  printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

#endif

  /* getcm(natoms,coords->px,coords->py,coords->pz,coords->amass,
     &coords->cmx,&coords->cmy,&coords->cmz); */
  /*printf("com= %lg %lg %lg\n",coords->cmx,coords->cmy,coords->cmz);*/
  /* allocate temporary vector */
  mat = dmatrix(0,3*natoms-1,0,3*natoms-1);
  
  for(i=0;i<3*natoms;i++){
    for(j=0;j<3*natoms;j++){
      qm = 0.;
      for(k=0;k<3*natoms;k++)
        qm += dipder[j][k]*nmodes->fcmat[k][i]/sqrt(coords->amass[k/3]);
      mat[j][i] = qm;
    }
  }

  printf("WARNING this needs to be changed for anything but water (3 atom molecules) \n");
  printf("WARNING from getmu_interface, interface MUST be perpendicular to z direction \n");
  
/* BEGIN loop over 3N INM modes */
  for(imode=0;imode<natoms*3;imode++){
    /* zero variables */
    kcount1 = kcount2 = 0;
    dx1 = dy1 = dz1 = 0.0;
    dx2 = dy2 = dz2 = 0.0;
    for(i=0;i<natoms;i++){
      ii = i*3;
      /* first atom in molecule determines if in top or bottom half */
      /* only coded for 3 atom molecules !!!!! */
      if(!(i%3)) {
        if(coords->pz[i] > coords->cmz)
          ihalf = 1;
        else
          ihalf = 2;
      }
#ifdef MUTEST
      ihalf =1;
#endif      
      
      /* add induced & permanent charge contributions */
      if(ihalf == 1){
        kcount1 += 1;
        qm = coords->qch[i]/sqrt(coords->amass[i]);
#ifndef MUTEST
        dx1 += mat[ii  ][imode];
        dy1 += mat[ii+1][imode];
        dz1 += mat[ii+2][imode];
#endif
        dx1 += qm*nmodes->fcmat[ii+0][imode];
        dy1 += qm*nmodes->fcmat[ii+1][imode];
        dz1 += qm*nmodes->fcmat[ii+2][imode];
      }
      else if(ihalf == 2) {
        kcount2 += 1;
        qm = coords->qch[i]/sqrt(coords->amass[i]);
#ifndef MUTEST
        dx2 += mat[ii  ][imode];
        dy2 += mat[ii+1][imode];
        dz2 += mat[ii+2][imode];
#endif        
        dx2 += qm*nmodes->fcmat[ii+0][imode];
        dy2 += qm*nmodes->fcmat[ii+1][imode];
        dz2 += qm*nmodes->fcmat[ii+2][imode];
      } else {
        printf("Warning: atom %d not assigned to top or bottom half\n",i);
        exit(0);
      }
      
    }
    
#ifdef WATER  

      /* APBS add dipole derivative depending on geometry via Hynes water dipole mode */

    
    for(j=0,iii=0;iii<natoms;iii++,j+=3){
      if(iii%3==0&&imode==0){
/* APBS now get d(dq)/dx gerivatives */ 
#ifdef MUTEST
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("CALLING getdipwater j = %d \n",j);
	printf("iii=%d \n",iii);
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
#endif
/*        getdipWater(j,x,y,z,dqo,dqh1,dqh2,&dqmpo,&dqmph1,&dqmph2); */
/*        getdipWater(iii,coords->px,coords->py,coords->pz,dqo,dqh1,dqh2,&dqmpo,&dqmph1,&dqmph2);
 */
#ifdef MUTEST
	/*	getnumdipWater(iii,x,y,z,dqo,dqh1,dqh2,&dqmpo,&dqmph1,&dqmph2); */
#endif
        if(coords->pz[iii] > coords->cmz)
          ihalf = 1;
        else
          ihalf = 2;
      }
      
#ifdef MUTEST
      ihalf =1;
#endif      

      if((iii%3)==0)qmp = dqmpo/sqrt(coords->amass[iii]);
      if((iii%3)==1)qmp = dqmph1/sqrt(coords->amass[iii]);
      if((iii%3)==2)qmp = dqmph2/sqrt(coords->amass[iii]);

      if(ihalf==1){
      dx1 += qmp*nmodes->fcmat[j+0][imode];
      dy1 += qmp*nmodes->fcmat[j+1][imode];
      dz1 += qmp*nmodes->fcmat[j+2][imode];    
      }
      if(ihalf==2){
      dx2 += qmp*nmodes->fcmat[j+0][imode];
      dy2 += qmp*nmodes->fcmat[j+1][imode];
      dz2 += qmp*nmodes->fcmat[j+2][imode];    
      }
      

      if((iii%3)==0){
        sum=0.;
        for(ll=0;ll<9;ll++){
          if(ll==0)masfac=coords->amass[iii];
          if(ll==1)masfac=coords->amass[iii];
          if(ll==2)masfac=coords->amass[iii];
          if(ll==3)masfac=coords->amass[iii+1];
          if(ll==4)masfac=coords->amass[iii+1];
          if(ll==5)masfac=coords->amass[iii+1];
          if(ll==6)masfac=coords->amass[iii+2];
          if(ll==7)masfac=coords->amass[iii+2];
          if(ll==8)masfac=coords->amass[iii+2];
          sum += 1./sqrt(masfac)*nmodes->fcmat[j+ll][imode]*dqo[ll];
        }
      }
      
      if((iii%3)==1){
        sum=0.;
        for(ll=0;ll<9;ll++){
          if(ll==0)masfac=coords->amass[iii-1];
          if(ll==1)masfac=coords->amass[iii-1];
          if(ll==2)masfac=coords->amass[iii-1];
          if(ll==3)masfac=coords->amass[iii];
          if(ll==4)masfac=coords->amass[iii];
          if(ll==5)masfac=coords->amass[iii];
          if(ll==6)masfac=coords->amass[iii+1];
          if(ll==7)masfac=coords->amass[iii+1];
          if(ll==8)masfac=coords->amass[iii+1];
          sum += 1./sqrt(masfac)*nmodes->fcmat[j+ll-3][imode]*dqh1[ll];
        }
      }
      
      if((iii%3)==2){
        sum=0.;
        for(ll=0;ll<9;ll++){
          if(ll==0)masfac=coords->amass[iii-2];
          if(ll==1)masfac=coords->amass[iii-2];
          if(ll==2)masfac=coords->amass[iii-2];
          if(ll==3)masfac=coords->amass[iii-1];
          if(ll==4)masfac=coords->amass[iii-1];
          if(ll==5)masfac=coords->amass[iii-1];
          if(ll==6)masfac=coords->amass[iii];
          if(ll==7)masfac=coords->amass[iii];
          if(ll==8)masfac=coords->amass[iii];
          sum += 1./sqrt(masfac)*nmodes->fcmat[j+ll-6][imode]*dqh2[ll];
        }
      }

      if(ihalf==1){
        dx1 += coords->px[iii] *sum;
        dy1 += coords->py[iii] *sum;
        dz1 += coords->pz[iii] *sum;
      }
      
      if(ihalf==2){
        dx2 += coords->px[iii] *sum;
        dy2 += coords->py[iii] *sum;
        dz2 += coords->pz[iii] *sum;
      }
    }
    

#ifdef MUTEST
    kcount1=1;
    kcount2=1;
#endif
/*endif for WATER */
#endif
    
   /* saving derivatives for inm analysis */
    nmodes->mux1[imode]  = dx1/(double)kcount1;
    nmodes->muy1[imode]  = dy1/(double)kcount1;
    nmodes->muz1[imode]  = dz1/(double)kcount1;
    nmodes->mux2[imode]  = dx2/(double)kcount2;
    nmodes->muy2[imode]  = dy2/(double)kcount2;
    nmodes->muz2[imode]  = dz2/(double)kcount2;
#ifdef MUTEST
    kcount1=1;
    kcount2=1;
    
    printf("==================in e- ===================\n");
    printf("mode imode = %d dmu_i/dx = %lg \n",imode,dx1/sqrt(CCONV));
    printf("mode imode = %d dmu_i/dy = %lg \n",imode,dy1/sqrt(CCONV));
    printf("mode imode = %d dmu_i/dz = %lg \n",imode,dz1/sqrt(CCONV));
    printf("mode imode = %d dmu_i/dx = %lg \n",imode,dx2/sqrt(CCONV));
    printf("mode imode = %d dmu_i/dy = %lg \n",imode,dy2/sqrt(CCONV));
    printf("mode imode = %d dmu_i/dz = %lg \n",imode,dz2/sqrt(CCONV));
    printf("=======================================\n");
#endif
/* close loop over imode */
  }
  
}

  




/*-------------------------------------------------------------------*/





