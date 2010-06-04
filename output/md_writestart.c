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

/*----------------------------------------------------------------------*/
void writestart(FILENAMES *filenames,int istep,SIMPARMS *simparms,
		ENERGIES *energies,COORDS *coords)
{
  int i,j;
  FILE *frest;

  frest  = cfopenw(filenames->restart);
  fprintf(frest,"%d %d %d %d %d %12g %12g %12g\n",
	  istep,simparms->natoms,simparms->ntherm,simparms->nchain,
	  simparms->nbar,simparms->temp,simparms->pext,simparms->dt);
  fprintf(frest,"%12g %12g %12g %12g %12g %12g %12g %12g\n",
	  energies->avham,energies->avpotra,energies->avpoter,energies->avzke,
	  energies->avpotn,energies->avzken,
	  energies->avpotv,energies->avzkev);
  fprintf(frest,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	  energies->avpot_inter,energies->avpot_bond,energies->avpot_bend,
	  energies->avpot_tors,energies->avpot_onfo,energies->avpot_onfo_e,
	  energies->avpot_elec,energies->avpot_recip,
	  energies->avpot_extern,energies->avpot_extern_e);
  fprintf(frest,"%12g %12g %12g %12g %12g\n",
	  energies->avtemp,energies->avprs,energies->avvol,
	  energies->econv,energies->econ2);
  for(i=0;i<9;i++) fprintf(frest,"%12g ",energies->avprs_tensor[i]);
  fprintf(frest,"\n");
  
  for(i=0;i<simparms->natoms;i++){
    fprintf(frest,"%.18g %.18g %.18g %.18g %.18g %.18g\n",
	    coords->px[i],coords->py[i],coords->pz[i],
	    coords->vx[i],coords->vy[i],coords->vz[i]);
  }
  fprintf(frest,"%16.10g %16.10g %16.10g\n",
	  coords->hmat[0],coords->hmat[1],coords->hmat[2]);
  fprintf(frest,"%16.10g %16.10g %16.10g\n", 
	  coords->hmat[3],coords->hmat[4],coords->hmat[5]);
  fprintf(frest,"%16.10g %16.10g %16.10g\n",
	  coords->hmat[6],coords->hmat[7],coords->hmat[8]);
  
  fprintf(frest,"%16.10g %16.10g %16.10g\n",
	  coords->vvol9[0],coords->vvol9[1],coords->vvol9[2]);
  fprintf(frest,"%16.10g %16.10g %16.10g\n", 
	  coords->vvol9[3],coords->vvol9[4],coords->vvol9[5]);
  fprintf(frest,"%16.10g %16.10g %16.10g\n",
	  coords->vvol9[6],coords->vvol9[7],coords->vvol9[8]);
  
  for(j=0;j<simparms->ntherm;j++){
    for(i=0;i<simparms->nchain;i++){
      fprintf(frest," %.12g %.12g",
	      coords->eta[j][i],coords->veta[j][i]);
    }
    fprintf(frest,"\n");
  }
#ifdef ISO
  fprintf(frest,"%.12g %.12g\n",coords->pvol,coords->vvol);
#endif
  for(i=0;i<simparms->nbar;i++){
    fprintf(frest,"%.12g %.12g\n",coords->bar[i],coords->vbar[i]);
  }
  fprintf(frest,"\n");
  fflush(frest);
  fclose(frest);
}
/*-------------------------------------------------------------------*/

void finish(FILENAMES *filenames,SIMPARMS *simparms,ENERGIES *energies,
	    COORDS *coords)
{
  md_stdout("Closing files");

  writestart(filenames,simparms->istep,simparms,energies,coords);

  write_coord(filenames->initfile,simparms->natoms,coords);

  if(filenames->fconf != NULL) fclose(filenames->fconf);  
  if(filenames->fvel != NULL) fclose(filenames->fvel);  
  if(filenames->fham != NULL) fclose(filenames->fham);  
  if(filenames->feng != NULL) fclose(filenames->feng);
  if(filenames->fext != NULL) fclose(filenames->fext);
  if(filenames->fcolv != NULL) fclose(filenames->fcolv);
  if(simparms->icalc_type==4) {
    if(filenames->hillfile != NULL) {
      save_hills(filenames->hillfile,simparms->istep,simparms->dt,&coords->colvar);
    }
  }
}
/*-----------------------------------------------------------------*/
void write_coord(char *name,int natoms,COORDS *coords)
{
  int i;
  FILE *fp;
  WORD filename;
  LINE line;

  sprintf(filename,"%s_N",name);
  if((fp = fopen(filename,"w")) == NULL){
    fprintf(stderr,"ERROR: can't open file (%s) to write coords\n",filename);
    exit(1);
  }
  /* print header */
  sprintf(line,"Saving new initial file (%s) with final coordinates",filename);
  md_stdout(line);

  fprintf(fp,"%d %d\n",natoms,1);

  for(i=0;i<natoms;i++){
    fprintf(fp,"%10g %10g %10g\n",coords->px[i],coords->py[i],coords->pz[i]);
  }
  for(i=0;i<3;i++){
    fprintf(fp,"%g %g %g\n",coords->hmat[i*3],
	    coords->hmat[i*3+1],coords->hmat[i*3+2]);
  }
  fclose(fp);
}
/*-----------------------------------------------------------------*/
void write_all_coords(SIMPARMS *simparms,COORDS *coords)
{
  int i,j;

  for(i=0;i<simparms->natoms;i++){
    printf("%d %9g %9g %9g\n%9g %9g %9g\n%9g %9g %9g\n",i,
	   coords->px[i],coords->py[i],coords->pz[i],
	   coords->vx[i],coords->vy[i],coords->vz[i],
	   coords->fxt[i],coords->fyt[i],coords->fzt[i]);
  }
  for(i=0;i<simparms->ntherm;i++){
    for(j=0;j<simparms->nchain;j++){
      printf("%d %d %g %g\n",i,j,coords->eta[i][j],coords->veta[i][j]);
    }
  }
  printf("%g %g\n",coords->pvol,coords->vvol);
  for(i=0;i<simparms->nbar;i++){
    printf("%g %g\n",coords->bar[i],coords->vbar[i]);
  }
}
/*-----------------------------------------------------------------*/
