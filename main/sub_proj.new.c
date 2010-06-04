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


/* Calculate the summed projection of the instantaeous modes
   onto the fixed subspace basis 
   */

#include "md.h"

static double govl_new(int il,int ih,int jl,int jh,double **rov);

#define PROG_CO2
/* #define LIFE */

#define FCUTOFF (1.1)    /* 1.1 shuts it off cutoff for writing to over.dat */
#define NTIME 50        /* number of configurations held in memory */
#define NCONF 500        /* number of configurations to corrilate */

#ifdef JUNK
#define NCOMP 10

static int numm_l[NCOMP] = { 0, 0, 50,100,100,200,250,350,400, 0};
static int numm_h[NCOMP] = {50,100,100,200,250,250,350,400,450,250};
static int numn_l[NCOMP] = { 0, 0, 50,100,100,200,250,350,400, 0};
static int numn_h[NCOMP] = {50,100,100,200,250,250,350,400,450,250};
#endif

#define NCOMP 1

static int numm_l[NCOMP] = {100};
static int numm_h[NCOMP] = {200};
static int numn_l[NCOMP] = {100};
static int numn_h[NCOMP] = {200};


/* static variables that can be seen only functions in this file */
static int nsave;
static int *time_store,*ilife;
static double **rov,***cormat,**rov_dcomp;
static double *comp[NCOMP],*compt[NCOMP],*compr[NCOMP],*compb[NCOMP];
static double *comps[NCOMP],*compa[NCOMP];

#ifdef LIFE
static double **clife2,**afreq;
#endif

void init_subproj_new(SIMPARMS *simparms)
{
  int i,j;
  LINE line;

  i = simparms->natoms*3-1;
  rov = dmatrix(0,i,0,i);
  rov_dcomp = dmatrix(0,i,0,8);
  cormat = d3tensor(0,NTIME-1,0,i,0,i);
  time_store = malloc(NTIME*sizeof(int));

  if( rov==NULL || cormat ==NULL || time_store==NULL){
    fprintf(stderr,"SORRY, can't allocate memory for storage of configs\n");
    exit(1);
  }

#ifdef LIFE
  afreq = dmatrix(0,i,0,NCONF-1);
  clife2 = dmatrix(0,i,0,NCONF-1);
  if(afreq==NULL || clife2==NULL || ilife==NULL ){
    fprintf(stderr,"SORRY, can't allocate memory for storage to corr_func\n");
    exit(1);
  }
  simparms->mem_bytes += 3*((i+1)*NCONF*sizeof(double));
#endif

  ilife = malloc(NCONF*sizeof(int));

  for(i=0;i<NCOMP;i++){
    comp[i] = malloc(NCONF*sizeof(double));
    compt[i] = malloc(NCONF*sizeof(double));
    compr[i] = malloc(NCONF*sizeof(double));
    compb[i] = malloc(NCONF*sizeof(double));
    comps[i] = malloc(NCONF*sizeof(double));
    compa[i] = malloc(NCONF*sizeof(double));
    
    if(comp[i]==NULL || compt[i]==NULL || compr[i]==NULL || compb[i]==NULL || 
       comps[i]==NULL || compa[i]==NULL ){
      fprintf(stderr,"SORRY, can't allocate memory for storage to comps.%d\n",
	      i);
      exit(1);
    }
    simparms->mem_bytes += 9*NCONF*sizeof(double);
  }

  simparms->mem_bytes += (NTIME+1)*(i+1)*(i+1)*sizeof(double);
  simparms->mem_bytes += (i+1)*8*sizeof(double);
  simparms->mem_bytes += (NTIME+NCONF)*sizeof(int);

  sprintf(line,"Bytes of memory allocated so far is %d",simparms->mem_bytes);
  md_stdout(line);

#ifdef DEBUG
  if(NTIME>NCONF) {
    fprintf(stderr,
	    "GOOFBALL! you can't save more configs than you want to avg\n");
    exit(1);
  }
#endif
  nsave = NCONF/NTIME; /* how often to save a new configuration */

  /* ZERO corr functions */
  for(i=0;i<NCONF;i++) ilife[i] = 0;

  for(j=0;j<NCOMP;j++){
    for(i=0;i<NCONF;i++){
      comp[j][i]=compt[j][i]=compr[j][i]= 0.0;
      compb[j][i]=comps[j][i]=compa[j][i]=0.0;
    }
  }
#ifdef LIFE
  for(i=0;i<3*simparms->natoms;i++){
    for(j=0;j<NCONF;j++) afreq[i][j] = clife2[i][j] = 0.;
  }
#endif
}

/*************************************************************************/

void sub_proj_new(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *nmodes)
{
  int i,j,k,l,itime,ntime,natoms3,nconf,icomp,nl,nh;
  double brn,rke;
  double **tmat,**pmat;
  char name[50];
  FILE *fphs;
#ifdef LIFE
  double arn;
#endif

  natoms3 = simparms->natoms*3;
  nconf = nmodes->nconf-1;
  pmat = nmodes->fcmat;

  /* sort eigenvalues in ascending order */
  vsrtasnd(nmodes->dn,pmat,natoms3,natoms3);

#ifdef PROG_CO2

  if((nconf%nsave)==0){
    printf("**** storing new eigenvectors **** \n");
    /* move matrix pointers forward one */
    tmat = cormat[NTIME-1]; 
    k = time_store[NTIME-1];
    for(i=NTIME-1;i>=1;i--){
      cormat[i] = cormat[i-1];
      time_store[i] = time_store[i-1];
    }
    cormat[0] = tmat;
    time_store[0] = k;

    /* store new subspace */
    time_store[0] = nconf;
    for(i=0;i<natoms3;i++){
      for(j=0;j<natoms3;j++){
        cormat[0][i][j] = pmat[i][j];
      }
    }
  }

  printf("@@ Calculating projection @@ \n");
  
  ntime = MIN(nconf/nsave+1,NTIME);
  for(itime=0;itime<ntime;itime++){

    tmat = cormat[itime];
    pmat = nmodes->fcmat;
    for(i=0;i<natoms3;i++){
      for(j=0;j<natoms3;j++){
	rke=0.;
	for(k=0;k<natoms3;k++) rke += pmat[k][i]*tmat[k][j];   
	rov[i][j] = rke*rke; /* store the squares */
      }
    }

    pmat = nmodes->molmat;
    tmat = cormat[itime];
    for(i=0;i<natoms3;i++){  /* all eigenvectors */
      for(j=0;j<9;j++){      /* all trans rots bend sym asym */
	rov_dcomp[i][j] = 0.;
	for(k=0;k<natoms3;k+=9){  /* each gas phase modes */
	  rke=0.;
	  for(l=k;l<k+9;l++) rke += pmat[l][j]*tmat[l][i];  
	  rov_dcomp[i][j] += rke*rke; /* sum the squares */
	}
	if(rov_dcomp[i][j]>1.)printf("dcomp = %d %d %g\n",i,j,rov_dcomp[i][j]);
      }
    }  

    k = nconf - time_store[itime];
    /* printf(" k = %d %d %d %d\n",k,nconf,time_store[itime],itime); */

    if(k>=NCONF){
      fprintf(stdout,"ERROR you twit nconf > time_store\n");
      exit(1);
    }
    ilife[k]++;

#ifdef LIFE
    for(i=0;i<natoms3;i++){     
      arn=FCONV*sqrt(fabs(nmodes->dn[i]));
      if(fabs(nmodes->dn[i])!=nmodes->dn[i]) arn=-arn;
      
      clife2[i][k] += brn;
      afreq[i][k] += arn;
    }
#endif

    /*---------------- COMPs ------------------------*/

    for(icomp=0;icomp<NCOMP;icomp++){
      nl = numn_l[icomp];
      nh = MIN(numn_h[icomp],natoms3);

      comp[icomp][k] += govl(nl,nh,numm_l[icomp],numm_h[icomp],rov);
      
      /* translations */
      compt[icomp][k] += govl(nl,nh,0,3,rov_dcomp);
      
      /* rotations */
      compr[icomp][k] += govl(nl,nh,3,5,rov_dcomp);
      
      /* bending */
      compb[icomp][k] += govl(nl,nh,5,7,rov_dcomp);
      
      /* s-stretch */
      comps[icomp][k] += govl(nl,nh,0,9,rov_dcomp);

      /* a-stretch */
      compa[icomp][k] += govl(nl,nh,0,5,rov_dcomp);
    }
  }
  
  /* AVERAGE DONE WRITE OUT CORRELATION FUNCTIONS !!! */
  
  for(icomp=0;icomp<NCOMP;icomp++){
    sprintf(name,"sub%d.%d_%d.%d.dat",
	    numm_l[icomp],numm_h[icomp],numn_l[icomp],numn_h[icomp]);
    if((fphs=fopen(name,"w"))==NULL){
      fprintf(stderr,"ERROR:can't open \"%s\"\n",name);
    }
    fprintf(fphs,"# totals all states  step nstate sum_yke sum_xke\n");
    ntime = MIN(nconf+1,NCONF);
    for(k=0;k<ntime;k++){
      brn =  1./(double)(ilife[k]);
      fprintf(fphs,"%9g %9g %9g %9g %9g %9g %9g\n",k*simparms->dt,
	      comp[icomp][k]*brn,compt[icomp][k]*brn,compr[icomp][k]*brn,
	      compb[icomp][k]*brn,comps[icomp][k]*brn,compa[icomp][k]*brn);
    }
    fclose(fphs);
  }
  
#ifdef LIFE
  for(i=0;i<natoms3;i++){
    sprintf(name,"life/%d",i);
    if((fphs=fopen(name,"w"))==NULL){
      fprintf(stderr,"can't open \"%s\"\n",name);
      exit(1);
    }
    fprintf(fphs,"# time  abs(rov[i][i]) rov[i][i]**2 <freq>\n");

    ntime = MIN(nconf+1,NCONF);
    for(k=0;k<ntime;k++){
      brn = 1./(double)(ilife[k]);
      fprintf(fphs,"%g %g %g\n",k*simparms->dt,
	      clife2[i][k]*brn,afreq[i][k]*brn);
    }
    fclose(fphs);
  }
#endif

#endif /* PROG_CO2 */

}
/*-------------------------------------------------------------*/
static double govl_new(int il,int ih,int jl,int jh,double **rov)
{
  int i,j;
  double arn;
  
  arn = 0.;
  for(i=il;i<ih;i++){
    for(j=jl;j<jh;j++){
      arn += rov[i][j];
    }
  }
  return arn;
}
/*-------------------------------------------------------------*/
