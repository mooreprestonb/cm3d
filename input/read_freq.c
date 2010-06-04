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


/* file to read in frequency data */

#include "md.h"

void read_freqfile(char *filename,int ndim,double *svec)
{
  int n,i;
  double sv;
  FILE *fp;
  

  printf("%s %d\n",filename,ndim);
  if((fp = fopen(filename,"r"))==NULL){
    fprintf(stderr,"ERROR: can't open %s\n",filename);
    exit(1);
  }
  /* read header */

  fscanf(fp,"%d",&n);

  if(n!=ndim){
    fprintf(stderr,"ERROR: %d != %d (in file %s)\n",ndim,n,filename);
    exit(1);
  }

  for(i=0;i<n;i++){
    fscanf(fp,"%lg",&(svec[i]));
  }

  for(i=0;i<n;i++){
    sv = svec[i];
    sv /= FCONV;     /* convert to picoseconds */
    sv *= sv;      /* square them up for eigenvalues */
    svec[i] = sv;
  }
}

