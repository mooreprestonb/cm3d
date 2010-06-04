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
   C program to calculate the normal mode spectrum for configurations
   dumped from an md run.  It reads in the parameter file and configuration
   files 
  */

#include "md.h"

/* #define DEBUG */
/* #define FRIC */
/* #define SUBPROJ */
/* #define SPEC_ONLY */

#define SPEC

void nma_control_new(char *command,FILENAMES *filenames,
		 SIMPARMS *simparms,COORDS *coords,INTER *inter,
		 WRITE_STEP *write_step,SUBSPACE *subspace)
{
  int i,j,n3;
  NMODES nmodes;
  LINE line;
  double **dipder,**temp;

  printf("%s running INM analysis \n",command);
  nmodes.maxfreq = subspace->freq_min;  
  nmodes.minfreq = subspace->freq_max;
  nmodes.nconf_off = 0;
  nmodes.imodes = IMODES;

#ifdef SPEC_ONLY
  nmodes.flo = 0;
#else
  nmodes.flo = 1;
#endif

  /* allocate local memory */
  n3 = 3*simparms->natoms;
  coords->ux = cmalloc(n3*sizeof(double));
  coords->uy = coords->ux+  simparms->natoms;
  coords->uz = coords->ux+2*simparms->natoms;
  coords->scr_buf=(double *)malloc(6*simparms->natoms*sizeof(double));
  coords->scr_rec=coords->scr_buf + 3*simparms->natoms;

  dipder = dmatrix(0,n3-1,0,n3-1);
  temp = dmatrix(0,n3-1,0,n3-1);
  nma_allocate(simparms,&nmodes);
  sprintf(line,"Bytes of memory allocated so far is %d",simparms->mem_bytes);
  md_stdout(line);

  nmodes.df = (subspace->freq_max-subspace->freq_min)/
    (double)(simparms->npoints-1);

  read_header_conf(simparms->natoms,&nmodes.nconf,&simparms->dt);

  if(simparms->istart==2){
    read_cont(filenames->nmfile,simparms->natoms,
	      &subspace->freq_min,&subspace->freq_max,
	      simparms->npoints,&nmodes);
  }
  write_nma_start(simparms->natoms,simparms->npoints,
		  simparms->dt,subspace->freq_max,
		  subspace->freq_min,simparms->iunits,
		  write_step->ndump,&nmodes);

#ifdef SUBPROJ
  init_subproj(simparms);
#endif

  nmodes.nconf = 0;
  nmodes.acpu = cputime();


  sprintf(line,"\n***** Starting nma loop ***** (dt = %g)\n",simparms->dt);
  md_stdout(line);

  /* loop over all stored configurations */
  while(read_config(simparms->natoms,coords->px,coords->py,coords->pz,
     coords->hmat)){
    nmodes.nconf++;
    /* Get force constant matrix */
    fmat(simparms,coords,inter,&nmodes);

    for(i=0;i<n3;i++) for(j=0;j<n3;j++) temp[i][j] = nmodes.fcmat[i][j];
    
    /* diagonalize system get eigenvectors and eigenvalues */
#ifdef SPEC_ONLY
    rs_me(simparms->natoms*3,nmodes.dn,nmodes.fcmat,0); /* no eigenvectors */
#else
    rs_me(simparms->natoms*3,nmodes.dn,nmodes.fcmat,1);
#endif
    
    /* sort eigenvalues in ascending order */
    vsrtasnd(nmodes.dn,nmodes.fcmat,simparms->natoms*3,simparms->natoms*3);
    
#ifdef DEBUG
    printf("nconf = %d\n",nmodes.nconf);
    printf("FMAT:\n"); 
    printmatrix(simparms->natoms*3,nmodes.fcmat);

    for(i=0;i<n3;i++) for(j=0;j<n3;j++) dipder[i][j] = 0.;
    for(i=0;i<n3;i++)
      for(j=0;j<n3;j++){
	int k;
        for(k=0;k<n3;k++) dipder[i][j] += temp[i][k]*nmodes.fcmat[k][j];
      }
    printf("F V= L V\n");
    for(i=0;i<n3;i++)
      for(j=0;j<n3;j++)
	printf("%lg\t%lg\t%lg\n",dipder[j][i],nmodes.fcmat[j][i],nmodes.dn[i]);
#endif

    partratio(nmodes.fcmat,nmodes.dn,nmodes.part,simparms->natoms);

    if(write_step->psysvec){
      /* print eigenvectors */
      psysvec(filenames,nmodes.nconf,simparms->natoms,simparms->iunits,
	      simparms->npoints,coords->px,coords->py,coords->pz,coords->amass,
	      coords->hmat,nmodes.fcmat,nmodes.dn,nmodes.molmat,
	      nmodes.nspec,nmodes.nmol,nmodes.napm);
      
      /* write_coord(filenames->initfile,simparms->natoms,coords); */
      
      sprintf(line,"Wrote out eigenvectors system of config %d",nmodes.nconf);
      md_stdout(line);
    }
    /* get principle axis and angles then rotate known eigvec to mol axis*/
    if(nmodes.flo){
      get_axis(simparms->natoms,nmodes.nspec,nmodes.nmol,nmodes.napm,
         coords,nmodes.molmat);
      overper(simparms->natoms,nmodes.nspec,nmodes.nmol,nmodes.napm,
         nmodes.fcmat,nmodes.molmat,nmodes.permod); 
#ifdef FRIC
      friction(simparms,coords,inter,&nmodes);
#endif
#ifdef SUB_PROJ
      sub_proj(simparms,coords,inter,&nmodes);
#endif
      
#ifdef SPEC

/* calculate induced dipoles for IR spectrum */
      if(simparms->nma_type == 0 || simparms->nma_type == 2){  
        for(i=0;i<n3;i++) for(j=0;j<n3;j++) dipder[i][j] = 0.;
        getdid(simparms,coords,inter,dipder);
        /* for(i=0;i<n3;i++) dipder[i][i] = coords->qch[i/3]; */
        /* get dipole induced dipole contribution */
        get_spectrum(simparms->natoms,coords->amass,coords->qch,nmodes.fcmat,
           nmodes.fricmod,dipder);
        /* get charge induced dipole contribution */
        get_specq(simparms->natoms,coords->amass,coords->qch,nmodes.fcmat,
           nmodes.fricmodq);
        /* get both dipole induced dipole and charge induced dipole contribution */
        get_spect_full(simparms->natoms,coords->amass,coords->qch,nmodes.fcmat,
           nmodes.fricmodcq,dipder);
        avg_spec(simparms->natoms,&nmodes);
      }

/* calculate polarizability derivatives for RAMAN spectrum */
      if(simparms->nma_type == 1 || simparms->nma_type == 2){
        getpolder(simparms,coords,inter,&nmodes);
      }
#endif

      if(write_step->peigval) /* store frequencies */
        store_freq(simparms,n3,nmodes.dn,nmodes.nconf*simparms->dt);

/* bin frequencies, create density of states and calculate spectra */
      freqbin(simparms,simparms->npoints,subspace->freq_min,&nmodes);
      
      nmodes.acpu += cputime();
      /* saving along the way... */
      if(write_step->ndump && !(nmodes.nconf%write_step->ndump)){
        save_nma(simparms,filenames,simparms->npoints,subspace->freq_min,
           &nmodes);
      }
    
    }/* endif write out vectors */
  }
  /* Print out final Data */ 
  save_nma(simparms,filenames,simparms->npoints,subspace->freq_min,
   &nmodes);
  free_dmatrix(dipder,0,n3-1,0,n3-1);
  free_dmatrix(temp,0,n3-1,0,n3-1);
  printf("***Finished INM Analysis***\n");
}
/*------------------------------------------------------------------------*/
void nma_allocate_new(SIMPARMS *simparms,NMODES *nmodes)
{
  int i,ispec,imol,ioff;
  
  /* count species */
  nmodes->nspec = simparms->ispecies[simparms->natoms-1]+1;
  fprintf(stdout,"Allocating nma arrays, nspec = %d\n",nmodes->nspec);

  /* count number of atom/molecule */
  nmodes->napm = (int *)cmalloc(nmodes->nspec*sizeof(int));
  nmodes->nmol = (int *)cmalloc(nmodes->nspec*sizeof(int));
  
  simparms->mem_bytes += simparms->nspec*2*sizeof(int);

  for(i=0;i<nmodes->nspec;i++)  nmodes->napm[i] = nmodes->nmol[i] = 0.;

  i=ispec=ioff = 0;
  while(i<simparms->natoms){
    while(i+nmodes->napm[ispec] < simparms->natoms && 
	  ispec==simparms->ispecies[i+nmodes->napm[ispec]]){
      if(simparms->imolecule[i+nmodes->napm[ispec]]!=0) break;
      nmodes->napm[ispec]++;
    }
    while(i<simparms->natoms){
      if(ispec!=simparms->ispecies[i]) break;
      i++;
    }
    nmodes->nmol[ispec] = (i-ioff)/nmodes->napm[ispec];
    ioff = i;
    ispec++;
  }
  
  nmodes->suspx = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  nmodes->suspy = nmodes->suspx+  simparms->natoms;
  nmodes->suspz = nmodes->suspx+2*simparms->natoms;
  simparms->mem_bytes += simparms->natoms*3*sizeof(double);

  fprintf(stdout,"Allocated susceptability %d\n",simparms->mem_bytes);

  nmodes->fcmat = dmatrix(0,3*simparms->natoms-1,0,3*simparms->natoms-1);
  simparms->mem_bytes += 9*simparms->natoms*simparms->natoms*sizeof(double);

  fprintf(stdout,"Allocated Force Matrix %d\n",simparms->mem_bytes);

  nmodes->dn  = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  nmodes->part  = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  nmodes->en  = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  nmodes->en2 = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  nmodes->fpart = (double *)cmalloc(simparms->npoints*sizeof(double));
  nmodes->freq = (double *)cmalloc(simparms->npoints*sizeof(double));

  simparms->mem_bytes += 12*simparms->natoms*sizeof(double);
  simparms->mem_bytes += 2*simparms->npoints*sizeof(double);
  
/* IR */
  if(simparms->nma_type == 0 || simparms->nma_type == 2){  
    nmodes->fric = (double *)cmalloc(simparms->npoints*sizeof(double));
    nmodes->fpos = (double *)cmalloc(simparms->npoints*sizeof(double));
    nmodes->ftot = (double *)cmalloc(simparms->npoints*sizeof(double));
    for(i=0;i<simparms->npoints;i++) {
      nmodes->freq[i] = nmodes->fric[i] = nmodes->fpos[i] = nmodes->ftot[i] = 0.;
      simparms->mem_bytes += 3*simparms->npoints*sizeof(double);
    }
  }
  /* Raman */
  if(simparms->nma_type == 1 || simparms->nma_type == 2){  
    nmodes->redram = (double *)cmalloc(simparms->npoints*sizeof(double));
    nmodes->raman = (double *)cmalloc(simparms->npoints*sizeof(double));
    nmodes->oke = (double *)cmalloc(simparms->npoints*sizeof(double));
    nmodes->iso = (double *)cmalloc(simparms->npoints*sizeof(double));
    for(i=0;i<simparms->npoints;i++) {
      nmodes->redram[i] = nmodes->raman[i] = nmodes->oke[i] = nmodes->oke[i] = 0.;
      simparms->mem_bytes += 4*simparms->npoints*sizeof(double);
    }
  }
  
  fprintf(stdout,"Allocated points for nma plots %d\n",simparms->mem_bytes);
  
  ioff = nmodes->nspec*nmodes->imodes;
  nmodes->permod = dmatrix(0,ioff-1,0,simparms->natoms*3-1);
  nmodes->decmod = dmatrix(0,ioff-1,0,simparms->npoints-1);
  nmodes->molmat = (double **)cmalloc(simparms->natoms*3*sizeof(double *));
  nmodes->avg_freq = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  nmodes->avg_strn = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  for(i=0;i<simparms->natoms*3;i++){
    nmodes->avg_freq[i] =  nmodes->avg_strn[i] = 0.;
  }
  /* tally memory usage */
  simparms->mem_bytes += ioff*simparms->natoms*3*sizeof(double);
  simparms->mem_bytes += ioff*simparms->npoints*sizeof(double);
  simparms->mem_bytes += 3*simparms->natoms*3*sizeof(double);
  /* IR */
  if(simparms->nma_type == 0 || simparms->nma_type == 2){  
    nmodes->decmodd = dmatrix(0,ioff-1,0,simparms->npoints-1);
    nmodes->decmodq = dmatrix(0,ioff-1,0,simparms->npoints-1);
    nmodes->decmoddq = dmatrix(0,ioff-1,0,simparms->npoints-1);
    nmodes->fricmod = (double *)cmalloc(simparms->natoms*9*sizeof(double));
    nmodes->fricmodq = nmodes->fricmod + 3*simparms->natoms;
    nmodes->fricmodcq = nmodes->fricmod + 6*simparms->natoms;
    /* tally memory usage */
    simparms->mem_bytes += 3*ioff*simparms->npoints*sizeof(double);
    simparms->mem_bytes += 3*simparms->natoms*3*sizeof(double);
  }
  /* Raman */
  if(simparms->nma_type == 1 || simparms->nma_type == 2){  
    nmodes->decmodrr = dmatrix(0,ioff-1,0,simparms->npoints-1);
    nmodes->aniso_raman = (double *)cmalloc(simparms->natoms*3*sizeof(double));
    nmodes->decmod_isorr = dmatrix(0,ioff-1,0,simparms->npoints-1);
    nmodes->iso_raman = (double *)cmalloc(simparms->natoms*3*sizeof(double));
    /* tally memory usage */
    simparms->mem_bytes += 2*ioff*simparms->npoints*sizeof(double);
    simparms->mem_bytes += simparms->natoms*6*sizeof(double);
  }
  
  fprintf(stdout,"Allocated arrays for IR & Raman spectra %d\n",simparms->mem_bytes);

  ioff = 0;
  for(ispec=0;ispec< nmodes->nspec;ispec++){
    for(imol=0;imol< nmodes->nmol[ispec]; imol++){
      for(i=0;i<nmodes->napm[ispec];i++){
        nmodes->molmat[ioff++]=cmalloc((nmodes->napm)[ispec]*3*sizeof(double));
        nmodes->molmat[ioff++]=cmalloc((nmodes->napm)[ispec]*3*sizeof(double));
        nmodes->molmat[ioff++]=cmalloc((nmodes->napm)[ispec]*3*sizeof(double));
        simparms->mem_bytes += (nmodes->napm)[ispec]*9*sizeof(double);
      }
    }
  }
}





