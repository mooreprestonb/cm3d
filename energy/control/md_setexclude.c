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


/* routine to set up the exclusion lists */

/*-----------------------------------------------------------------------*/
#include "md.h"

/* #define DEBUG */
/* #define WRITE_EXCLUSION */
/* #define ALL_POS */

/*-----------------------------------------------------------------------*/
void set_exclude(SIMPARMS *simparms,COORDS *coords,INTER *inter)
{
  int i,j,num_excl,mem_exclude;
  LINE line;

  /* set up exclusion list */
  if(simparms->rank==0){ md_stdout("Allocating exclusion list memory");}

  inter->exclude  = (int**)cmalloc(simparms->natoms*sizeof(int *));
  inter->nexclude = (int *)cmalloc(simparms->natoms*sizeof(int));
  num_excl = 0;
  simparms->mem_bytes += simparms->natoms*sizeof(int);
  simparms->mem_bytes += simparms->natoms*sizeof(int *);

  mem_exclude = simparms->max_exclude;

  for(i=0;i<simparms->natoms;i++){
    inter->nexclude[i] = 0;
    inter->exclude[i] = (int *)cmalloc((mem_exclude+1)*sizeof(int));
#include "vectorize.h"
    for(j=0;j<=mem_exclude;j++) inter->exclude[i][j] = -1;

#ifdef ALL_POS
    inter->exclude[i] = (int *)cmalloc((simparms->natoms+1-i)*sizeof(int));
#include "vectorize.h"
    for(j=0;j<=simparms->natoms-i;j++) inter->exclude[i][j] = -1;
#endif
  }

  if(simparms->rank==0) md_stdout("Excluding bonded interactions");
  exclude_bonds(inter->nexclude,inter->exclude,&coords->bonds);
  exclude_bends(inter->nexclude,inter->exclude,&coords->bends);
  exclude_tors(inter->nexclude,inter->exclude,&coords->tors);
  exclude_onfo(inter->nexclude,inter->exclude,&coords->onfo);
  
  for(i=0;i<simparms->natoms;i++){
    if(inter->nexclude[i]>=mem_exclude){
      sprintf(line,"memory of exclusion is not enough increase mem_exclude");
      sprintf(line,"from %d to above %d (atom number %d has %d exclusions)",
	      mem_exclude,inter->nexclude[i],i,inter->nexclude[i]);
      md_error(line);
    }
    num_excl += inter->nexclude[i];
  }
  if(simparms->rank==0){
    sprintf(line,"There were %d interaction exclusions found",num_excl);
    md_stdout(line);
  }

  set_ecorr(simparms,coords,inter->nexclude,inter->exclude);

  /* exclude NULL interactions */

#ifdef EXCLUDE_NULL
  num_null = 0;
  /* need more memory */
  for(i=0;i<simparms->natoms-1;i++){
    for(j=i+1;j<simparms->natoms;j++){
      if(inter->map[inter->itype[i]][inter->itype[j]] == -1){
	num_null += insert(inter->exclude[i],&inter->nexclude[i],j);
      }
    }
  }
  if(simparms->rank==0){
    sprintf(line,"There were %d null interactions added to exclusions",
	    num_null);
    md_stdout(line);
  }
#endif

  for(i=0;i<simparms->natoms;i++){
    inter->exclude[i] = (int *)realloc(inter->exclude[i],
				       (inter->nexclude[i]+1)*sizeof(int));
    simparms->mem_bytes += (inter->nexclude[i]+1)*sizeof(int);
  }
#ifdef WRITE_EXCLUSION
  printf("exclusion list %d\n",num_excl);
  num_excl = 0;
  for(i=0;i<simparms->natoms;i++){
    for(j=0;j<inter->nexclude[i];j++){
      printf("%d = %d %d %d\n",num_excl++,i,inter->exclude[i][j],j);
    }
  }
#endif
  if(simparms->rank==0){
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
}
/*-----------------------------------------------------------------------*/
