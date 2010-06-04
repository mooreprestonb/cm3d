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

/* subroutines dealing with input for normal mode analysis */

#include "md.h"

static FILE *fp_conf;

void open_nma_in(char *confile)
{
  char line[MAXLINELEN];
  if((fp_conf = fopen(confile,"r"))==NULL){
    sprintf(line,"can't open configuration file %s",confile);
    md_error(line);
  }
}
/*------------------------------------------------------------------*/
void read_header_conf(int natoms,int *nconf,double *dt)
{
  int npart;
  char line[MAXLINELEN];

  if(fgets(line,MAXLINELEN,fp_conf) == NULL){
    md_error("can't read header of the configuration file");
  }
  if(sscanf(line,"%*c%d %d %lg",&npart,nconf,dt)!=3){
    md_error("something is wrong with the header in configuration file");
  }
  *dt *= (double) *nconf;
  if(npart != natoms){
    md_error("Number of atoms don't agree between configuration header and the input files!");
  }
}
/*------------------------------------------------------------------*/
int read_config(int natoms,double *px,double *py,double *pz,double *hmat)
{
  char line[MAXLINELEN];
  int i;
  static double cell;
  static double cell_old = 0.0;

  if(fgets(line,MAXLINELEN,fp_conf) == NULL) {
    return 0;
  }
  if(sscanf(line,"%lg %lg %lg",&px[0],&py[0],&pz[0]) != 3){
    md_error("while reading in configurations");
  }
  for(i=1;i<natoms;i++){
    if(fgets(line,MAXLINELEN,fp_conf) == NULL){
      md_error("while reading in configurations");
    }
    if(sscanf(line,"%lg %lg %lg",&px[i],&py[i],&pz[i]) != 3){
      md_error("while reading in configurations");
    }
  }
  if(fgets(line,MAXLINELEN,fp_conf) == NULL){
    md_error("while reading in configurations");
  }
  if(sscanf(line,"%lg %lg %lg %lg %lg %lg %lg %lg %lg",
	    &hmat[0],&hmat[1],&hmat[2],
	    &hmat[3],&hmat[4],&hmat[5],
	    &hmat[6],&hmat[7],&hmat[8]) != 9){
    md_error("while reading in configurations");
  }
  
  cell = get_deth(hmat);
  if(cell_old == 0.0){
    cell_old = cell;
  }
  
  /* check that volume hasn't moved by much */
  if(fabs(1.- cell/cell_old)>.5){
    sprintf(line,"volume has changed more the 1/2 %g %g",cell,cell_old);
    md_error(line);
  }
  cell_old = cell;

  return 1;
}
/*-----------------------------------------------------------------------*/
void read_cont(char *outfile,int natoms,double *frmin,double *frmax,
	       int npt,NMODES *nmodes)
{
  char line[MAXLINELEN];
  int i,j,nptd,natomd;
  double fr,*x;
  FILE *fp_nma;
  
  /* read in data so far */
  if((fp_nma = fopen(outfile,"r"))==NULL){
    sprintf(line,"Can't open %s Thus, one can't continue",outfile);
    md_error(line);
  }
  if(fgets(line,MAXLINELEN,fp_nma) == NULL){
    sprintf(line,"Can't read header line in file %s thus, one can't continue",
	    outfile);
    md_error(line);
  }
  if(sscanf(line,"%*s %d %d",&natomd,&nmodes->nconf_off) != 2){
    sprintf(line,"in header line in file %s thus, one can't continue",
	    outfile);
    md_error(line);
  }
  if(natoms != natomd){
    sprintf(line,"File to be continued does not have the same # of atoms");
    sprintf(line,"%s NOT continuing.\n\t%s %d != %d",
	    line,outfile,natomd,natoms);
    md_error(line);
  }
  /* count lines */
  nptd = 0;
  while(fgets(line,MAXLINELEN,fp_nma) != NULL) nptd++;
  rewind(fp_nma);
  if(nptd != npt){
    sprintf(line,"number of points in %s %d != %d input number",
	    outfile,nptd,npt);
    md_error(line);
  }

  /* reread header */
  fgets(line,MAXLINELEN,fp_nma);

  fscanf(fp_nma,"%lg %lg",frmin,&nmodes->freq[0]);    
  for(j=0;j<nmodes->imodes;j++) fscanf(fp_nma,"%lg",&nmodes->decmod[j][0]);
  for(i=1;i<npt-1;i++){
    fscanf(fp_nma,"%*g %lg",&nmodes->freq[i]);
    for(j=0;j<nmodes->imodes;j++) fscanf(fp_nma,"%lg",&nmodes->decmod[j][i]);
  }
  fscanf(fp_nma,"%lg %lg",frmax,&nmodes->freq[npt-1]);
  for(j=0;j<nmodes->imodes;j++) fscanf(fp_nma,"%lg",&nmodes->decmod[j][i]);
  fclose(fp_nma);
  
  nmodes->df = (*frmax - *frmin)/(double)(npt-1);
  fr = 3.*natoms*(nmodes->nconf_off)*(nmodes->df);
  for(i=0;i<npt;i++){
    nmodes->freq[i] *= fr;
    for(j=0;j<nmodes->imodes;j++) nmodes->decmod[j][i] *= fr;
  }
  
  /* allocate dummy positions */
  x = (double *)cmalloc((3*natoms+9)*sizeof(double));
  
  for(j=0;j<nmodes->nconf_off;j++){
    if(read_config(natoms,&x[0],&x[natoms],&x[2*natoms],&x[3*natoms])){
      md_error("EOF reached in config file, not continuing");
    }
  }
  /* free temporary memory */
  free(x);
}
/*-----------------------------------------------------------------------*/
