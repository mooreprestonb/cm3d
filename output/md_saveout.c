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
/*---------------------------------------------------------*/
/* configuration */
void save_conf(FILE *fconf,int natoms,COORDS *coords)
{
  int ii;

  if(fconf==NULL){
    md_warning("No file name for configuration file, specify in input file\n");
    return;
  }

  for(ii=0;ii<natoms;ii++){
    fprintf(fconf,"%g %g %g\n",coords->px[ii],coords->py[ii],coords->pz[ii]);
  }
  fprintf(fconf,"%g %g %g %g %g %g %g %g %g\n",
	  coords->hmat[0],coords->hmat[1],coords->hmat[2],
	  coords->hmat[3],coords->hmat[4],coords->hmat[5],
	  coords->hmat[6],coords->hmat[7],coords->hmat[8]);

  fflush(fconf);
}
/*--------------------------------------------------------*/
void save_vel(FILE *fvel,int natoms,COORDS *coords)
{
  int ii;

  if(fvel==NULL){
    md_warning("No file name for velocity file, specify in input file\n");
    return;
  }

  for(ii=0;ii<natoms;ii++){
    fprintf(fvel,"%g %g %g\n",coords->vx[ii],coords->vy[ii],coords->vz[ii]);
  }
  fflush(fvel);
}
/*--------------------------------------------------------*/
void save_force(FILE *fforce,int natoms,COORDS *coords)
{
  int ii;
  double am;
  double *fx,*fy,*fz;
  
  if(fforce==NULL){
    md_warning("No file name for force file, specify in input file\n");
    return;
  }

  fx = coords->fxt;  fy = coords->fyt;  fz = coords->fzt;
  for(ii=0;ii<natoms;ii++) fx[ii] = fy[ii] = fz[ii] = 0.;

  for(ii=0;ii<natoms;ii++){
    fx[ii] += coords->fxl[ii]+coords->fxa[ii]+coords->fxr[ii];
    fy[ii] += coords->fyl[ii]+coords->fya[ii]+coords->fyr[ii];
    fz[ii] += coords->fzl[ii]+coords->fza[ii]+coords->fzr[ii];
  }

  for(ii=0;ii<natoms;ii++){
    am = coords->amass[ii];
    fprintf(fforce,"%g %g %g\n",fx[ii]*am,fy[ii]*am,fz[ii]*am);
  }
  fflush(fforce);
}
/*--------------------------------------------------------*/
void save_freezeforce(int nfreeze,COORDS *coords)
{
  int ii,j;
  double am;

  double *fx = coords->ffreeze;
  static FILE *fp=NULL;

  if(fp==NULL){
    fp = fopen("freeze_force","w");
  }

  for(ii=0;ii<nfreeze;++ii){
    j = coords->ifreeze[ii];
    am = coords->amass[j];
    fprintf(fp,"%g %g %g %d\n",fx[ii*DIM]*am,fx[ii*DIM+1]*am,fx[ii*DIM+2]*am,j);
  }
}

/*--------------------------------------------------------*/
void save_eigvec(char *filename,int m,int n,double *dn,double **d2v)
{
  char line[MAXLINELEN];
  int i,j;
  FILE *fp;
  
  if((fp=fopen(filename,"w"))==NULL){
    sprintf(line,"can't open \"%s\" for saving eigenvalues",filename);
    md_error(line);
  }
  /* write header */
  fprintf(fp,"%d %d\n",m,n);

  for(i=0;i<m;i++){
    fprintf(fp,"%g\n",dn[i]);
  }
  for(i=0;i<m;i++){
    fprintf(fp,"\n");
    for(j=0;j<n;j++){
      fprintf(fp,"%g\n",d2v[j][i]);
    }
  }
  fclose(fp);
}
/*---------------------------------------------------------*/
#ifndef BUFSIZ
#define BUFSIZ 1024
#endif

void keep(char *name)
{
  char save[MAXLINELEN];
  int n,buf[BUFSIZ];
  FILE *fpi,*fpo;
  LINE line;

  sprintf(save,"%s.O",name);

  sprintf(line,"Saving file initial files as %s",save);
  md_stdout(line);

  if((fpi= fopen(name,"r"))==NULL){
    fprintf(stderr,"ERROR: can't open file %s for reading\n",name);
    return;
  }
  if((fpo= fopen(save,"w"))==NULL){
    fprintf(stderr,"ERROR: can't open file %s for writing\n",save);
    return;
  }
  while( (n = fread(buf,sizeof(int),BUFSIZ,fpi)) > 0 ){
    fwrite(buf,sizeof(int),n,fpo);
  }
  fwrite(buf,sizeof(int),n,fpo);

  fclose(fpi);
  fclose(fpo);
}
/*-----------------------------------------------------------*/
