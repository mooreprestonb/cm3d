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

#define LOAD_PVEC
#define STRATT_METH

/* #define DEBUG */

#ifdef DEBUG
#define PRINT_MOMENTS
#define PRINT_AXIS
#define PRINT_EIGVEC
#define PRINT_SUM
#define PRINT_MODI 
#endif


#ifndef DBL_MAX
#define DBL_MAX 1e+37
#endif

/* #define CHAR_MODE */
/*------------------------------------------------------------------------*/
void get_axis(int natoms,int nspec,int *nmol,int *napm,COORDS *coords,
	      double **molmat)
{
  int namol,namol3,ispec,ioff,*idx;
  int i,j,k,l,imol,imode;
  double cx,cy,cz,mass,mmass,sum,xnew,ynew,znew;
  static double *dx,*dy,*dz,*r2,**axis;
  static double *d,*e,**pmolmat,*amass;
  
  ioff = 0;
  for(ispec=0;ispec<nspec;ispec++){
    namol = napm[ispec];
    namol3 = namol*3;

#ifdef DEBUG
    printf("spec %d %d napm = %d\n",ispec,nmol[ispec],namol);
#endif
    
    idx = (int *)cmalloc(IMODES*sizeof(int));
    axis = dmatrix(0,(int)DIM-1,0,(int)DIM-1);
    dx = (double *)cmalloc(namol*sizeof(double));
    dy = (double *)cmalloc(namol*sizeof(double));
    dz = (double *)cmalloc(namol*sizeof(double));
    r2 = (double *)cmalloc(namol*sizeof(double));
    d = (double *)cmalloc(namol3*sizeof(double));
    e = (double *)cmalloc(namol3*sizeof(double));
    pmolmat = dmatrix(0,namol3-1,0,namol3-1);
    
    for(imol=0;imol<nmol[ispec];imol++){
      /* calculate center of mass */      
      cx = cy = cz = 0.;
      mmass = 0.;
      for(j=ioff,i=0;i<namol;i++,j++){
	mmass += mass = coords->amass[j];
	cx += coords->px[j]*mass; 
	cy += coords->py[j]*mass; 
	cz += coords->pz[j]*mass;
      }
      cx /= mmass; cy /= mmass; cz /= mmass;
      
      for(i=0,j=ioff;i<namol;i++,j++){
	dx[i] = coords->px[j] - cx; 
	dy[i] = coords->py[j] - cy; 
	dz[i] = coords->pz[j] - cz;
	r2[i] = dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i];
      }

      for(i=0;i<(int)DIM;i++) for(j=0;j<(int)DIM;j++) axis[i][j] = 0.;
      
      for(i=0;i<namol;i++){
	mass = coords->amass[i+ioff];
	axis[0][0] += mass*(r2[i] - dx[i]*dx[i]);
	axis[1][1] += mass*(r2[i] - dy[i]*dy[i]);
	axis[2][2] += mass*(r2[i] - dz[i]*dz[i]);
	
	axis[0][1] -= mass*dx[i]*dy[i];
	axis[0][2] -= mass*dx[i]*dz[i];
	axis[1][2] -= mass*dy[i]*dz[i];
      }
      
      axis[1][0] = axis[0][1];
      axis[2][0] = axis[0][2];
      axis[2][1] = axis[1][2];
      
      /* solve for eigenvectors which are the principle axis*/
      
#ifdef PRINT_AXIS
      printf("{");
      for(i=0;i<(int)DIM-1;i++){
	printf("{");
	for(j=0;j<(int)DIM-1;j++)
	  printf("%10g,",axis[i][j]);
	printf("%10g},\n",axis[i][j]);
      }
      printf("{");
      for(j=0;j<(int)DIM-1;j++)
	printf("%10g,",axis[i][j]);
      printf("%10g}}\n",axis[i][j]);
#endif
      
      if(namol >1){
	rs_me((int)DIM,d,axis,1);
	ceigsrtv(d,axis,(int)DIM);

#ifdef PRINT_MOMENTS
	printf("Moments of intertia:principle frame!!!\n");
	printf("%g %g %g",d[0], d[1], d[2]);
#endif
	
#ifdef PRINT_AXIS
	printf("\n{");
	for(i=0;i<namol3-1;i++){
	  printf("{");
	  for(j=0;j<namol3-1;j++)
	    printf("%10g,",molmat[i+ioff*(int)DIM][j]);
	  printf("%10g},\n",molmat[i+ioff*(int)DIM][j]);
	}
	printf("{");
	for(j=0;j<namol3-1;j++)
	  printf("%10g,",molmat[i+ioff*(int)DIM][j]);
	printf("%10g}}\n",molmat[i+ioff*(int)DIM][j]);
#endif
	
	/* get eigenvectors of isolated molecule */

	/* load molmat into pmolmat so that memory is contiguous for
	   library matrix call */
	for(l=i=0;i<namol;i++,l+=3){
	  k = ioff*(int)DIM+l;
	  for(j=0;j<namol3;j++){
	    pmolmat[l  ][j] = molmat[k  ][j];
	    pmolmat[l+1][j] = molmat[k+1][j];
	    pmolmat[l+2][j] = molmat[k+2][j];
	  }
	} 


	/* get and sort the eigenvectors */
	rs_me(namol3,d,pmolmat,1);
	vsrtasnd(d,pmolmat,namol3,namol3);

	/* load the solution back into molmat */
	for(l=i=0;i<namol;i++,l+=3){
	  k = ioff*(int)DIM+l;
	  for(j=0;j<namol3;j++){
	    molmat[k  ][j] = pmolmat[l  ][j];
	    molmat[k+1][j] = pmolmat[l+1][j];
	    molmat[k+2][j] = pmolmat[l+2][j];
	  }
	} 

#ifdef PRINT_EIGVEC
	for(i=0;i<namol3;i++) printf("%g ",d[i]);
	printf("\n{");
	for(i=0;i<namol3-1;i++){
	  printf("{");
	  for(j=0;j<namol3-1;j++)
	    printf("%8.5g,",molmat[i+ioff*(int)DIM][j]);
	  printf("%8.5g},\n",molmat[i+ioff*(int)DIM][j]);
	}
	printf("{");
	for(j=0;j<namol3-1;j++)
	  printf("%8.5g,",molmat[i+ioff*(int)DIM][j]);
	printf("%8.5g}}\n",molmat[i+ioff*(int)DIM][j]);
#endif
	/* find three lowest magnitude eigenvalues */
	idx[2] = idx[1] = idx[0] = -1;
	cx = cy = cz = DBL_MAX;
	for(i=0;i<namol3;i++){
	  if(fabs(d[i])<cx){
	    cz = cy;cy = cx;cx = fabs(d[i]);
	    idx[2] = idx[1];idx[1] = idx[0];idx[0] = i;
	  } else if(fabs(d[i])<cy){
	    cz = cy;cy = fabs(d[i]);
	    idx[2] = idx[1];idx[1] = i;
	  } else if(fabs(d[i])<cz){
	    cz = fabs(d[i]);
	    idx[2] = i;
	  }
	}
#ifdef DEBUG
	printf("lowest eigenvectors are %d %d %d with values %g %g %g\n",
	       idx[0],idx[1],idx[2],d[idx[0]],d[idx[1]],d[idx[2]]);
#endif
      } else {
	d[0] = d[1] = d[2] = 0.;
	idx[0] = 0;idx[1] = 1;idx[2] = 2;
	axis[0][0] = axis[1][1] = axis[2][2] = 1.;
	axis[0][1] = axis[0][2] = axis[1][0] = 0.;
	axis[1][2] = axis[2][0] = axis[2][1] = 0.;
      }
#ifdef LOAD_PVEC
      /* load principle axis into vectors */
      for(i=0;i<namol;i++){
	k = (i+ioff)*(int)DIM;
	mass = sqrt(coords->amass[i+ioff]);
	for(j=0;j<(int)DIM;j++){
	  molmat[k+j][idx[0]] = axis[j][2]*mass;
	  molmat[k+j][idx[1]] = axis[j][1]*mass;
	  molmat[k+j][idx[2]] = axis[j][0]*mass;
	}
      }
      xnew = ynew = znew = 0.;
      for(i=0,j=ioff*(int)DIM;i<namol3;i++,j++){
	xnew += molmat[j][idx[0]]*molmat[j][idx[0]];
	ynew += molmat[j][idx[1]]*molmat[j][idx[1]];
	znew += molmat[j][idx[2]]*molmat[j][idx[2]];
      }
      xnew = 1./sqrt(xnew);
      ynew = 1./sqrt(ynew);
      znew = 1./sqrt(znew);
      for(i=0,j=ioff*(int)DIM;i<namol3;i++,j++){
	molmat[j][idx[0]] *= xnew;
	molmat[j][idx[1]] *= ynew;
	molmat[j][idx[2]] *= znew;
      }
#endif

      /* load molmat into pmolmat and unmass weight them */
      for(l=i=0;i<namol;i++,l+=3){
	k = ioff*(int)DIM+l;
	mass = 1./sqrt(coords->amass[i+ioff]);
	for(j=0;j<namol3;j++){
	  pmolmat[l  ][j] = molmat[k  ][j]*mass;
	  pmolmat[l+1][j] = molmat[k+1][j]*mass;
	  pmolmat[l+2][j] = molmat[k+2][j]*mass;
	}
      } 
      /* renormalize to SUM a^2 = 1 */
      for(i=0;i<namol3;i++){
	sum = 0.;
	for(j=0;j<namol3;j++){
	  sum += pmolmat[j][i]*pmolmat[j][i];
	}
	sum = 1./sqrt(sum);
	for(j=0;j<namol3;j++){
	  pmolmat[j][i] *= sum;
	}
      }

      /* rotate into principle frame */
      for(i=0;i<namol;i++){
	k = i*3;
	for(j=0;j<namol3;j++){
	  xnew = ynew = znew = 0.;
	  for(l=0;l<(int)DIM;l++){
	    xnew += pmolmat[k+l][j]*axis[l][2];
	    ynew += pmolmat[k+l][j]*axis[l][1];
	    znew += pmolmat[k+l][j]*axis[l][0];
	  }
	  pmolmat[k  ][j] = xnew;
	  pmolmat[k+1][j] = ynew;
	  pmolmat[k+2][j] = znew;
	}
	xnew = (dx[i]*axis[0][2]+dy[i]*axis[1][2]+dz[i]*axis[2][2]);
	ynew = (dx[i]*axis[0][1]+dy[i]*axis[1][1]+dz[i]*axis[2][1]);
	znew = (dx[i]*axis[0][0]+dy[i]*axis[1][0]+dz[i]*axis[2][0]);
	dx[i] = xnew; dy[i] = ynew;  dz[i] = znew;
      }

#ifdef PRINT_EIGVEC
      for(i=0;i<namol3;i++) printf("%9.5f ",d[i]);
      printf("\n{\n");
      for(i=0;i<namol3-1;i++){
	printf("{");
	for(j=0;j<namol3-1;j++)
	  printf("%9.5f,",pmolmat[i][j]);
	printf("%9.5f},\n",pmolmat[i][j]);
      }
      printf("{");
      for(j=0;j<namol3-1;j++)
	printf("%9.5f,",pmolmat[i][j]);
      printf("%9.5f}}\n",pmolmat[i][j]);
#endif

      for(imode=0;imode<IMODES;imode++) idx[imode] = -1;
      for(i=0;i<namol3;i++) e[i]=0.;
      /* transition */
      for(j=0;j<namol3;j++){
	cx = cy = cz = 0.; 
	for(i=0;i<namol3;i+=3){
	  cx += pmolmat[i  ][j];
	  cy += pmolmat[i+1][j]; 
	  cz += pmolmat[i+2][j];
	}
	cx = fabs(cx);	cy = fabs(cy);	cz = fabs(cz);
#ifdef CHAR_MODE
	if(e[0] < cx && cx>cy && cx>cz){e[0] = cx; idx[0] = j;} 
	if(e[1] < cy && cy>cx && cy>cz){e[1] = cy; idx[1] = j;} 
	if(e[2] < cz && cz>cx && cz>cy){e[2] = cz; idx[2] = j;}
#endif
#ifdef STRATT_METH
	sum = sqrt(cx*cx+cy*cy+cz*cz);
	if(e[0] < sum){
	  e[2] = e[1];e[1] = e[0];e[0] = sum;
	  idx[2] = idx[1];idx[1] = idx[0];idx[0] = j;
	} else if(e[1] < sum){
	  e[2] = e[1];e[1] = sum;
	  idx[2] = idx[1];idx[1] = j;
	} else if(e[2]<sum){
	  e[2] = sum;
	  idx[2] = j;
	}
#endif /* STRATT_METH */
      }

      /* rotation calculate "torque" of the mode */
      if(namol > 1){
	for(j=0;j<namol3;j++){
	  cx = cy = cz = 0; /* rotation */
	  for(i=k=0;i<namol;i++,k+=3){
	    amass = coords->amass;
	    cx += amass[i]*(dy[i]*pmolmat[k+2][j]-dz[i]*pmolmat[k+1][j]);
	    cy += amass[i]*(dz[i]*pmolmat[k  ][j]-dx[i]*pmolmat[k+2][j]);
	    cz += amass[i]*(dx[i]*pmolmat[k+1][j]-dy[i]*pmolmat[k  ][j]);
	    /* 
	       cx += (dy[i]*pmolmat[k+2][j]-dz[i]*pmolmat[k+1][j]);
	       cy += (dz[i]*pmolmat[k  ][j]-dx[i]*pmolmat[k+2][j]);
	       cz += (dx[i]*pmolmat[k+1][j]-dy[i]*pmolmat[k  ][j]);
	       */
	  }
	  
	  cx = fabs(cx);cy = fabs(cy);cz = fabs(cz);
#ifdef CHAR_MODE
	  if(e[3] < cx && cx>cy && cx>cz){e[3] = cx; idx[3] = j;} 
	  if(e[4] < cy && cy>cx && cy>cz){e[4] = cy; idx[4] = j;} 
	  if(e[5] < cz && cz>cx && cz>cy){e[5] = cz; idx[5] = j;}
#endif
#ifdef STRATT_METH
	  sum = sqrt(cx*cx+cy*cy+cz*cz);
	  if(e[3]<sum){
	    e[5] = e[4];e[4] = e[3];e[3] = sum;
	    idx[5] = idx[4];idx[4] = idx[3];idx[3] = j;
	  } else if(e[4]<sum){
	    e[5] = e[4];e[4] = sum;
	    idx[5] = idx[4];idx[4] = j;
	  } else if(e[5]<sum){
	    e[5] = sum;
	    idx[5] = j;
	  }
#endif
#ifdef PRINT_SUM
	  printf("rotation %d %g %d %g %d %g %d %g\n",j,sum,
		 idx[3],e[3],idx[4],e[4],idx[5],e[5]);
#endif
	}
      } /* endif (namol>1) */
#ifdef PRINT_MODI
      if(namol==1){
	/* printf("atom (no principle axis calculated)\n"); */
      } else {
	for(i=0;i<6;i++){
	  printf("%d %d %10g %9.5f\n",i,idx[i],e[i],d[idx[i]]);
	}
      }
#endif
      if(namol==1) idx[3] = idx[4] = idx[5] = -1;
#ifdef STRATT_METH
      if(namol==2) idx[5] = -1;
#endif
#ifdef CHAR_MODE
#ifdef LOAD_PVEC
      if(namol==2) idx[3] = -1;
#else
      if(namol==2) idx[5] = -1;
#endif
#endif
      /* sort gas phase modes */
      for(i=0;i<6;i++){
	if(idx[i] != -1){
	  sum = d[i];
	  d[i] = d[idx[i]];
	  d[idx[i]] = sum;
	  for(j=0,k=ioff*(int)DIM;j<namol3;j++,k++){
	    sum = molmat[k][i];
	    molmat[k][i] = molmat[k][idx[i]]; 
	    molmat[k][idx[i]] = sum;
	  }
	  for(j=i;j<6;j++){
	    if(i == idx[j]) {
	      idx[j] = idx[i];
	    }
	  }
	  idx[i] = i;
	}
      }
      ioff += namol; /* update the atom offset */
    }	/* loop over next molecule */
    /* free up temporary vectors */
    free(idx);
    free_dmatrix(axis,0,(int)DIM-1,0,(int)DIM-1);
    free(dx); free(dy); free(dz); free(r2); free(d); free(e);
    free_dmatrix(pmolmat,0,namol3-1,0,namol3-1);

  } /* loop over next species */
}

/*---------------------------------------------------------------------*/
void overper(int natoms,int nspec,int *nmol,int *napm,
	     double **mmat,double **nm,double **permod)
{
  int i,j,k,l,m,ispec,imol,natoms3,napm3,ioff,mod_off;
  double *sumi,sum,tsum;
  
/* mmat is the force constant matrix fcmat, nm is molmat */

  natoms3 = 3*natoms;
  sumi = (double *)cmalloc(natoms3*sizeof(double));
  for(i=0;i<natoms3;i++) for(j=0;j<nspec*IMODES;j++) permod[j][i] = 0.0;
  
  for(l=0;l<natoms3;l++){ 
    tsum = 0.;
    ioff = 0;
    for(ispec=0;ispec<nspec;ispec++){ /* loop over species */
      napm3 = napm[ispec]*3;
      for(imol=0;imol<nmol[ispec];imol++){ /* loop over molecules */
	for(k=0;k<napm3;k++){ /* loop over atoms*3 */
	  sum = 0.;
	  j = ioff*3+napm3;
	  for(m=ioff*3;m<j;m++){sum += nm[m][k]*mmat[m][l];}
	  tsum += fabs(sum);
	  sumi[ioff*3+k] = fabs(sum);
	}
	ioff += napm[ispec];
      }
    }
    if(tsum != 0.0){
      tsum = 1./tsum;
      ioff = 0;
      for(ispec=0;ispec<nspec;ispec++){
	mod_off = IMODES*ispec;
	napm3 = napm[ispec]*3;
	for(imol=0;imol<nmol[ispec];imol++){
	  for(k=0;k<napm3;k++){
	    sum = sumi[ioff*3+k]*tsum;
	    if(k<(IMODES-1)){
	      permod[k+mod_off][l] += sum;
	    }else {
	      permod[(IMODES-1)+mod_off][l] += sum;
	    }
	  }
	  ioff += napm[ispec];
	}
      }
    }
  }
  free(sumi);
}
/*--------------------------------------------------------------------*/
