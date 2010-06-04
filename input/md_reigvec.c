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
void read_eigvec(char *filename,int m,int n,double *dn,double **d2v)
{
  char line[MAXLINELEN];
  int i,j,mt,nt;
  FILE *fp;
  
  if((fp=fopen(filename,"r"))==NULL){
    sprintf(line,"can't open \"%s\" for reading eigenvalues",filename);
    md_error(line);
  }
  /* write header */
  fscanf(fp,"%d %d\n",&mt,&nt);
  if(m!=mt || n != nt ){
    fprintf(stderr,"ERROR: header in eigenvector file \"%s\" does not match\n",
	    filename);
    fprintf(stderr,"\twith input values\n");
    fprintf(stderr,"nstate   = %d == %d\n",m,mt);
    fprintf(stderr,"nmol*3   = %d == %d\n",n,nt);
    exit(1);
  }
  for(i=0;i<m;i++){ fscanf(fp,"%lg",&dn[i]);}
  for(i=0;i<m;i++) for(j=0;j<n;j++){
    fscanf(fp,"%lg",&d2v[j][i]);
  }
  
  fclose(fp);
}
