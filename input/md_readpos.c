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

/* routine to read in the positions for initfile */

#include "md.h"

/*------------------------------------------------------------------*/
void readpos(char *command,char *initfile,SIMPARMS *simparms,
	     int *estep,ENERGIES *energies,COORDS *coords)
{
  char line[MAXLINELEN];
  int i,j,n,nt,nc,nb,ist;
  double ttemp,tdt,tpext;
  double *px,*py,*pz,*vx,*vy,*vz,**eta,**veta,*bar,*vbar;
  FILE *fp;

  *estep   = 0;
  px = coords -> px;  py = coords -> py;  pz = coords -> pz;
  vx = coords -> vx;  vy = coords -> vy;  vz = coords -> vz;
  eta = coords -> eta;  veta = coords -> veta;
  bar = coords -> bar;  vbar = coords -> vbar;
  energies->deltae2 = 0.;

#ifdef DEBUG
  sprintf(line,"istart = %d",simparms->istart);
  md_stdout(line);
#endif
  if(simparms->istart == 1){
    sprintf(line,"Doing cold start from initial file \"%s\"",initfile);
    md_stdout(line);
    
    if((fp = fopen(initfile,"r")) == NULL){
      sprintf(line,"%s: can't open control file (%s)",command,initfile);
      md_error(line);
    }
    /* read header  */
    if(fgets(line,MAXLINELEN,fp) == NULL){
      sprintf(line,"can't read header from infile (%s)",initfile);
      md_error(line);
    }
    if((sscanf(line,"%d %d",&n,&ist)) != 2){
      sprintf(line,"Something wrong with header from infile (%s)",initfile);
      md_error(line);
    }
    if(simparms->natoms != n){
      sprintf(line,"Different number of atoms natoms\n\tN in set file and parameter file = %d\n\tN in %s = %d",
	      simparms->natoms,initfile,n);
      md_error(line);
    }
    if(ist != simparms->istart){
      sprintf(line,"istart in initial file %d != %d in input file",
	      ist,simparms->istart);
      md_error(line);
    }
    
    for(i=0;i<simparms->natoms;i++){
      if(fgets(line,MAXLINELEN,fp) == NULL){
        sprintf(line,"Something is wrong with infile (%s) at line %d.",
		initfile,i+1);
        md_error(line);
      }
      
      if((sscanf(line,"%lg %lg %lg",&px[i],&py[i],&pz[i])) != 3){
        sprintf(line,"Something wrong with infile (%s) at line %d.",
           initfile,i+1);
        md_error(line);
      }
      vx[i] = vy[i] = vz[i] = 0.0;
    }
#ifdef ORIG_HMAT
    for(i=0;i<3;i++){
      if(fgets(line,MAXLINELEN,fp) == NULL){
        sprintf(line,"Something is wrong with infile (%s) at line %d. (hmat)",
           initfile,i+simparms->natoms+1);
        md_error(line);
      }
      if((sscanf(line,"%lg %lg %lg",
         &coords->hmat[i*3],&coords->hmat[i*3+1],&coords->hmat[i*3+2]))
         != 3){
        sprintf(line,"Something is wrong with infile (%s) at line %d. (hmat)",
           initfile,i+simparms->natoms+1);
        md_error(line);
      }
    }
#else
    for(i=0;i<9;i++){
      if((fscanf(fp,"%lg",&coords->hmat[i]))!=1){
        sprintf(line,"Something is wrong with infile (%s) at line %d. (hmat)",
           initfile,i+simparms->natoms+1);
        md_error(line);
      }
    }
    
#endif
    coords->pvol = log(get_deth(coords->hmat))/DIM;
    coords->vvol = 0.;
    gethinv9(coords->hmat,coords->hmati);
    for(i=0;i<simparms->nbar;i++) coords->bar[i] = coords->vbar[i] = 0.;
    fclose(fp);
  } else {
    sprintf(line,"Doing warm start from restart file \"%s\"",initfile);
    md_stdout(line);
    /* read dump file */
    if((fp = fopen(initfile,"r"))==NULL){
      sprintf(line,"can't open restart file %s",initfile);
      md_error(line);
    }
    /* read header  */
    if(fgets(line,MAXLINELEN,fp) == NULL){
      sprintf(line,"can't read 1st line of header from infile (%s)",initfile);
      md_error(line);
    }
    if(sscanf(line,"%d %d %d %d %d %lg %lg %lg",&simparms->istep,
       &n,&nt,&nc,&nb,&ttemp,&tpext,&tdt) != 8){
      sprintf(line,"can't read in header file of %s (line 1)",initfile);
      md_error(line);
    }
    if(simparms->natoms != n){
      sprintf(line,"Different number of atoms %d (initfile)",n);
      sprintf(line,"%s != %d (setfile+parmfile",line,simparms->natoms);
      md_error(line);
    }
    fgets(line,MAXLINELEN,fp);
    if(sscanf(line,"%lg %lg %lg %lg %lg %lg %lg %lg",
	      &(energies->avham),&(energies->avpotra),&(energies->avpoter),
	      &(energies->avzke),&(energies->avpotn),&(energies->avzken),
	      &(energies->avpotv),&(energies->avzkev)) != 8){
      sprintf(line,"while reading averages (line 2) (%s)",initfile);
      md_error(line);
    }
    fgets(line,MAXLINELEN,fp);
    if(sscanf(line,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
	      &(energies->avpot_inter),&(energies->avpot_bond),
	      &(energies->avpot_bend),&(energies->avpot_tors),
	      &(energies->avpot_onfo),&(energies->avpot_onfo_e),
	      &(energies->avpot_elec),&(energies->avpot_recip),
	      &(energies->avpot_extern),&(energies->avpot_extern_e))!= 10){
      sprintf(line,"while reading averages (line 3) (missing extern?) (%s)",
	      initfile);
      md_error(line);
    }

    fgets(line,MAXLINELEN,fp);
    if(sscanf(line,"%lg %lg %lg %lg %lg",&(energies->avtemp),
	      &(energies->avprs),&(energies->avvol),
	      &(energies->econv),&(energies->econ2)) != 5){
      sprintf(line,"while reading averages (line 4) (%s)",initfile);
      md_error(line);
    }
    fgets(line,MAXLINELEN,fp);
    if(sscanf(line,"%lg %lg %lg %lg %lg %lg %lg %lg %lg",
	      &(energies->avprs_tensor[0]),&(energies->avprs_tensor[1]),
	      &(energies->avprs_tensor[2]),&(energies->avprs_tensor[3]),
	      &(energies->avprs_tensor[4]),&(energies->avprs_tensor[5]),
	      &(energies->avprs_tensor[6]),&(energies->avprs_tensor[7]),
	      &(energies->avprs_tensor[8])) != 9){
      sprintf(line,"while reading averages (line 5) (%s)\n",initfile);
      sprintf(line,"possible old restart file (add 9 elements after line 4)");
      md_error(line);
    }
    for(i=0;i<n;i++){
      fgets(line,MAXLINELEN,fp);
      if(sscanf(line,"%lg %lg %lg %lg %lg %lg"
         ,&px[i],&py[i],&pz[i],&vx[i],&vy[i],&vz[i]) != 6){
        sprintf(line,"while reading in positions and velocites");
        sprintf(line,"%s at line %d (atom # %d)",line,i,i+5);
        md_error(line);
      }
    }
    /* read in hmatrix */
    for(i=0;i<3;i++){
      if(fgets(line,MAXLINELEN,fp) == NULL){
        sprintf(line,"something is wrong with infile \"%s\" at line %d (hmat)",
           initfile,i+simparms->natoms+5);
        md_error(line);
      }
      if(sscanf(line,"%lg %lg %lg",
         &coords->hmat[i*3],&coords->hmat[i*3+1],&coords->hmat[i*3+2])
         !=3){
        sprintf(line,"something is wrong with infile (%s) at line %d (hmat)",
           initfile,i+simparms->natoms+5);
        md_error(line);
      }
    }
    gethinv9(coords->hmat,coords->hmati);
    
    /* read in vol of hmat */
    for(i=0;i<3;i++){
      if(fgets(line,MAXLINELEN,fp) == NULL){
        sprintf(line,"something is wrong in file \"%s\" at line %d(vvol9)",
           initfile,i+simparms->natoms+5+3);
        md_error(line);
      }
      if(sscanf(line,"%lg %lg %lg",&coords->vvol9[i*3],&coords->vvol9[i*3+1],
         &coords->vvol9[i*3+2]) !=3){
        sprintf(line,"something is wrong in file \"%s\" at line %d (vvol9)\n",
           initfile,i+simparms->natoms+5+3);
        sprintf(line,
           "%s possible old restart file (add nine elements after hmat)",
           line);
        md_error(line);
      }
    }

    if(fabs((ttemp-simparms->temp)/TZ(ttemp))>ERRMAX){
      sprintf(line,"input temp %g != %g stored temp",simparms->temp,ttemp);
      md_warning(line);
      if(simparms->istart==2) samvel(simparms,coords,0);
    }
    if(fabs((tpext-simparms->pext)/TZ(tpext))>(ERRMAX*10000)){
      sprintf(line,"input pressure %g != %g stored pressure",
         simparms->pext,tpext);
      md_warning(line);
    }
    if(fabs((tdt-simparms->dt)/TZ(tdt))>ERRMAX){
      sprintf(line,"old time step %g != %g new time step",tdt,simparms->dt);
      md_warning(line);
    }
    if(simparms->nchain != nc || nt != simparms->ntherm){
      sprintf(line,"Different number of thermostats (will zero all)\n\t ");
      sprintf(line,"%s nchain = %d != %d (being used)\n\t",line,nc,
	      simparms->nchain);
      sprintf(line,"%s ntherm = %d != %d (being used)\n\t",line,nt,
	      simparms->ntherm);
      md_warning(line);
      for(j=0;j<simparms->ntherm;j++)
        for(i=0;i<simparms->nchain;i++) eta[j][i] = veta[j][i] = 0.;
      
      for(j=0;j<nt;j++) for(i=0;i<nc;i++) fscanf(fp,"%*g %*g");
    } else { 
      for(j=0;j<nt;j++){
        for(i=0;i<nc;i++){
          fscanf(fp,"%lg %lg",&eta[j][i],&veta[j][i]);
        }
      }
      if(simparms->istart==2){
	sprintf(line,"Type 2 start -> zeroing thermostats\n\t ");	
	for(j=0;j<simparms->ntherm;j++)
	  for(i=0;i<simparms->nchain;i++) eta[j][i] = veta[j][i] = 0.;
      }
    }

#ifdef ISO
    fscanf(fp,"%lg %lg",&coords->pvol,&coords->vvol);
    
    if(fabs(coords->pvol-log(get_deth(coords->hmat))/DIM)
       /get_deth(coords->hmat)>ERRMAX){
      sprintf(line,"Cell size don't match %12g (from cell coords) != %12g",
	      log(get_deth(coords->hmat))/DIM,coords->pvol);
      sprintf(line,"%s\n\t Resetting Volume to %g",
	      line,get_deth(coords->hmat));
      md_warning(line);
      coords->pvol = log(get_deth(coords->hmat))/DIM;
      coords->vvol = 0.;
    }
#endif
    if(simparms->nbar != nb){
      sprintf(line,"Different number of barostats\n\t");
      sprintf(line,"%s nbar = %d != %d (being used)",line,nb,simparms->nbar);
      md_warning(line);
    }

    /* zero new barostats */
    for(i=nb;i<simparms->nbar;i++) bar[i] = vbar[i] = 0.;
    /* read all available barostats */
    nb = MIN(simparms->nbar,nb);
    for(i=0;i<nb;i++) fscanf(fp,"%lg %lg",&bar[i],&vbar[i]);
    fclose(fp);
    keep(initfile);
  }

#ifdef DEBUG
  sprintf(line,"%d %p %p %p",simparms->natoms,px,py,pz);
  md_stdout(line);
  for(i=0;i<n;i++){
    sprintf(line,"%i %g %g %g %g %g %g",i,px[i],py[i],pz[i],vx[i],vy[i],vz[i]);
    md_stdout(line);
  }
  sprintf(line,"cell coords %g %g %g\t %g %g %g\t %g %g %g",
	  coords->hmat[0],coords->hmat[1],coords->hmat[2],
	  coords->hmat[3],coords->hmat[4],coords->hmat[5],
	  coords->hmat[6],coords->hmat[7],coords->hmat[8]);
  md_stdout(line);
  sprintf(line,"cellinv coords %g %g %g\t %g %g %g\t %g %g %g",
	  coords->hmati[0],coords->hmati[1],coords->hmati[2],
	  coords->hmati[3],coords->hmati[4],coords->hmati[5],
	  coords->hmati[6],coords->hmati[7],coords->hmati[8]);
  md_stdout(line);

  sprintf(line,"thermostats = %d by %d",simparms->ntherm,simparms->nchain);
  md_stdout(line);
  for(i=0;i<simparms->ntherm;i++){
    for(j=0;j<simparms->nchain;j++){
      sprintf(line,"%d %d %g %g",i,j,eta[i][j],veta[i][j]);
      md_stdout(line);
    }
  }
  sprintf(line,"barostats = %d",simparms->nbar); md_stdout(line);
  for(i=0;i<simparms->nbar;i++){
    sprintf(line,"%d %g %g",i,bar[i],vbar[i]); md_stdout(line);
  }
#endif

  sprintf(line,"%d Atoms read in with a Cell Volume of %g",
     simparms->natoms,get_deth(coords->hmat));
  md_stdout(line);
}
/*------------------------------------------------------------------*/
