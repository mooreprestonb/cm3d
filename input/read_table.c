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

/*
This file changes "read_table" to read an external "wall" potential with
four columns in it instead of two. The four columns are the distance from
the wall (same as before), the potential energy (same as before), the
derivative of the potential energy with respect to the distance from
the wall (which was calculated before by taking a numerical derivative
of the potential), and the derivative of the potential energy with respect
to the sphere radius (which is infinity), which is the quantity I need
for my free energy method for colloids.
*/

/* subroutine to read in a table of numbers and spline then on a new table */

#include <stdlib.h>
#include <stdio.h>
#include "md.h"

#define MAXLINE 256 
#define NORDER 4 /* degree N-1 */
  /* #define SPLINE */

void read_table(char *file,int ntable,double zmin,double zmax,
		double *vtab,double *dvtab,double *dv2tab)
{
  char line[MAXLINE];
  int i,j,num;
  double *x,*y,dx,r,r2,pt;
  FILE *fp;
#ifdef SPLINE
  double *z;
#endif
  /* open data base file */
  if((fp = fopen(file,"r")) == NULL){
    fprintf(stderr,"ERROR: can't open external file (%s)\n",file);
    exit(1);
  }
  /* read header */
  if(fgets(line,MAXLINE,fp)== NULL ){
    fprintf(stderr,"ERROR: can't read first line of (%s)\n",file);
    exit(1);
  }

  /* get number of lines */
  if(sscanf(line,"%d",&num)!=1){
    fprintf(stderr,"ERROR:can't read in header of potential file (%s)",file);
    exit(1);
  }

  /* allocate temp array */
  if((x=malloc(3*num*sizeof(double)))==NULL){
    fprintf(stderr,"ERROR:can't allocate memory (%ld) for potential (%s)\n",
	    (long)2*num*sizeof(double),file);
    exit(1);
  }
  y = x +   num; /* point y to correct place */
#ifdef SPLINE
  z = x + 2*num; /* point z to correct place */
#endif
  
  fprintf(stdout,"Reading external table file \"%s\" with %d numbers\n",
	  file,num);

  for(i=0;i<num;i++){
    if(fgets(line,MAXLINE,fp)==NULL){   /* read file */
      fprintf(stderr,"ERROR: while reading in %s at line %d\n",file,i);
      exit(1);
    } else {
      if((sscanf(line,"%lg %lg",&x[i],&y[i]))!= 2){
	sprintf(line,"something went wrong while reading in potential");
	sprintf(line,"%s at line %d",line,i+1);
	md_error(line);
      }
    }
  }
  fclose(fp);
  
  /* find max/min */
  
  if(zmin != x[0] || zmax != x[num-1]){
    fprintf(stderr,"ERROR, Zmin or Zmax givin %g != %g or %g != %g found\n",
	    zmin,x[0],zmax,x[num-1]);
    exit(1);
  }
  dx = (zmax- zmin)/(double)(ntable-1);

  /* tabulate points */

#ifdef SPLINE
  spline(x,y,num, 2e30, 2e30,z);

  for(i=0;i<ntable;i++){
    r = i*dx+ zmin;
    splint(x,y,z,num,r,&pt); 
    vtab[i] = pt;
  }
#else
  j = 0;
  for(i=0;i<ntable;i++){
    r = i*dx+ zmin;
    while(r>x[j+NORDER/2] && j<(num-NORDER)) j++;
    polint(&(x[j])-1,&(y[j])-1,NORDER,r,&vtab[i], &pt);
  }
#endif
  free(x);   /* free memory */

  /* calculate derivatives numerically */
  r = 1./(2.*dx);
  r2 = 1./(dx*dx);
  dvtab[0] =  (vtab[1]-vtab[0])/dx;
  dv2tab[0] = (vtab[2]-2.*vtab[1]+vtab[0])/(dx*dx);
  for(i=1;i<ntable-1;i++) {
    dvtab[i] =  r*(vtab[i+1]-vtab[i-1]);
    dv2tab[i] = r2*(vtab[i+1]-2.*vtab[i]+vtab[i-1]);
  }
  dvtab[ntable-1] = -(vtab[ntable-1]-vtab[ntable-2])/dx;
  dv2tab[ntable-1] = r2*(vtab[ntable-3]-2.*vtab[ntable-2]+vtab[ntable-1]);
}
