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

void save_colv(FILE *fcolv,int istep,double dt,COLVAR *colvar)
{
  int i;
  fprintf(fcolv," %g ",istep*dt);
  for(i=0;i<colvar->ncolvar;i++){
    fprintf(fcolv," %g %g ",colvar->pcolvar[i] , colvar->pistcolvar[i]);
  }
  fprintf(fcolv," \n ");
  fflush(fcolv);


}

void save_hills(char *hillfile,int istep,double dt,COLVAR *colvar)
{
  FILE *fhill;
  int ihill,icol;
  char line[MAXLINELEN];

  fhill = fopen(hillfile,"w");
  if(fhill==NULL){
    sprintf(line,"Can't open hill file \"%s\"",hillfile);
    md_error(line);
  }

  fprintf(fhill,"# hillfile ( %d , %d ) @ %g ps: # depth(K) (pos width)*ncol\n",
	  colvar->nhills,colvar->ncolvar,dt*istep);
  for(ihill=0;ihill<colvar->nhills;++ihill){
    fprintf(fhill,"%d %g ",ihill, colvar->t_hilldepth[ihill]);
    for(icol=0;icol<colvar->ncolvar;++icol){
      fprintf(fhill," %g %g ",colvar->t_pcolvar[icol][ihill],
	      colvar->t_hillwidth[icol][ihill]);
    }
    fprintf(fhill,"\n");
  }
  /* write current posistion */
  fprintf(fhill,"# %d %g",istep,istep*dt);
  for(icol=0;icol<colvar->ncolvar;++icol){
    fprintf(fhill," %.30lg %.30lg",colvar->pcolvar[icol],colvar->vcolvar[icol]);
  }
  fprintf(fhill,"\n");
  
  fclose(fhill);
}

void read_hills(char *hillfile,int istep,double dt,COLVAR *colvar)
{
  FILE *fhill;
  int ihill,icol,iread,nhills,ncolvar;
  char line[MAXLINELEN];
  double depth,*thwidth;
  double tdepth;

  fhill = fopen(hillfile,"r");
  if(fhill==NULL){
    sprintf(line,"Can't open hill file \"%s\"",hillfile);
    md_error(line);
  }

  fgets(line,MAXLINELEN,fhill);
  iread = sscanf(line,"# hillfile ( %d , %d )",&nhills,&ncolvar);
  if((ncolvar != colvar->ncolvar) || (iread!=2)){
    sprintf(line,"Header for hill restart is inconsistant, expecting\n# hillfile ( %d , %d ) ...",0,colvar->ncolvar);
    md_error(line);
  }
  thwidth = (double *)malloc(ncolvar*sizeof(double));

  /* store current values */
  tdepth = colvar->hilldepth;
  for(icol=0;icol<colvar->ncolvar;++icol){
    thwidth[icol] = colvar->hillwidth[icol];
  }

  for(ihill=0;ihill<nhills;++ihill){
    fscanf(fhill,"%d %lg ",&iread, &depth);

    colvar->hilldepth = depth;
    for(icol=0;icol<colvar->ncolvar;++icol){
      fscanf(fhill," %lg %lg",&colvar->pcolvar[icol],
	     &colvar->hillwidth[icol]);
    }
    add_hills(colvar);
  }
  
  /* read current position */
  iread = fscanf(fhill,"%*s %d %lg",&istep,&dt);
  for(icol=0;icol<colvar->ncolvar;++icol){
    fscanf(fhill,"%lg %lg",&colvar->pcolvar[icol],&colvar->vcolvar[icol]);
  }


  fclose(fhill);
  free(thwidth);

  /* reset values */
  colvar->hilldepth = tdepth;
  for(icol=0;icol<colvar->ncolvar;++icol){
    colvar->hillwidth[icol] = thwidth[icol];
  }
}

