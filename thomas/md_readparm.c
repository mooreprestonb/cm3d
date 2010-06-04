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

/* read in the parameter files */

#include "md.h"

/* #define PRINT_ATOMS */

/*--------------------------------------------------------------*/

void read_parmfile(SIMPARMS *simparms,COORDS *coords,STARTUP *startup,
		   INTER *inter,NGBR *ngbr)
{
  int i,j,k,ispec,ntypes,ioff;
  int iatom,ngroup,natoms,ngpm;
  int **ibo,**jbo; WORD **tbo;
  int **ibe,**jbe,**kbe; WORD **tbe;
  int **ito,**jto,**kto,**lto; WORD **tto;
  int **iof,**jof;
  int **ibonx,**jbonx,**kbonx; WORD **tbonx;
  WORD tgroup;
  SPECIES *spec_now;
  ATOM_TOPOL **atom_topol;
  LINE line;
  
  /* allocate temparary memory */
  atom_topol = (ATOM_TOPOL **)cmalloc(startup->nspec*sizeof(ATOM_TOPOL *));

  /* bonds */
  ibo = (int **)cmalloc(startup->nspec*sizeof(int*));
  jbo = (int **)cmalloc(startup->nspec*sizeof(int*));
  tbo = (WORD **)cmalloc(startup->nspec*sizeof(WORD*));

  /* bends */
  ibe = (int **)cmalloc(startup->nspec*sizeof(int*));
  jbe = (int **)cmalloc(startup->nspec*sizeof(int*));
  kbe = (int **)cmalloc(startup->nspec*sizeof(int*));
  tbe = (WORD **)cmalloc(startup->nspec*sizeof(WORD*));

  /* torsions */
  ito = (int **)cmalloc(startup->nspec*sizeof(int*));
  jto = (int **)cmalloc(startup->nspec*sizeof(int*));
  kto = (int **)cmalloc(startup->nspec*sizeof(int*));
  lto = (int **)cmalloc(startup->nspec*sizeof(int*));
  tto = (WORD **)cmalloc(startup->nspec*sizeof(WORD*));

  /* onefour */
  iof = (int **)cmalloc(startup->nspec*sizeof(int*));
  jof = (int **)cmalloc(startup->nspec*sizeof(int*));

  /* cross bonds */
  ibonx = (int **)cmalloc(startup->nspec*sizeof(int*));
  jbonx = (int **)cmalloc(startup->nspec*sizeof(int*));
  kbonx = (int **)cmalloc(startup->nspec*sizeof(int*));
  tbonx = (WORD **)cmalloc(startup->nspec*sizeof(WORD*));

  /* read in all connectivity data from files */
  read_top_file(simparms,startup,atom_topol,
		ibo,jbo,tbo,ibe,jbe,kbe,tbe,
		ito,jto,kto,lto,tto,iof,jof,ibonx,jbonx,kbonx,tbonx);
  
#ifdef DEBUG
  printf("Number of total atoms = %d\n",simparms->natoms);
#endif
  natoms = simparms->natoms;

  allocate_coords(simparms,coords);

  /* allocate memeory for maps and charges and masses */
  simparms->atom_map = (ATOM_MAP *)cmalloc(natoms*
					   sizeof((simparms->atom_map)[0]));
  
  simparms->atom = (WORD *)cmalloc(natoms*sizeof(WORD));
  simparms->group = (WORD *)cmalloc(natoms*sizeof(WORD));
  simparms->mole = (WORD *)cmalloc(startup->nspec*sizeof(WORD));
  simparms->nspec = startup->nspec;

  /* get all the atom types */
  ngroup = ntypes = iatom = 0;

  /* fill in the maps */
  for(ispec=0,spec_now=&startup->spec_root;ispec < startup->nspec;
      ispec++,spec_now=spec_now->next){
    
    /* save molecule type */
    strcpy((simparms->mole)[ispec],spec_now->name);

    /* set up to check for group type */
    strcpy((simparms->group)[ngroup],atom_topol[ispec][0].group);
    strcpy(tgroup,atom_topol[ispec][0].group);
    ngpm = 0;

    for(j=0;j<spec_now->napm;j++){

      /* check for new atom type */
      for(i=0;i<ntypes;i++){
	if(strcasecmp(atom_topol[ispec][j].type,(simparms->atom)[i]) == 0) {
	  break;
	}
      }
      if(i==ntypes){
	strcpy((simparms->atom)[ntypes],atom_topol[ispec][j].type);
	ntypes++;
      }

      /* check for new group type */
      if(strcasecmp(atom_topol[ispec][j].group,tgroup) != 0) {
	strcpy(tgroup,atom_topol[ispec][j].group);
	ngpm++;
      }

      for(k=0;k<spec_now->nmol;k++){
	ioff = iatom + j + k*spec_now->napm;
	(simparms->atom_map)[ioff].type = i;
	(simparms->atom_map)[ioff].atom = j;
	(simparms->atom_map)[ioff].molecule = k;
	(simparms->atom_map)[ioff].species = ispec;
	(simparms->atom_map)[ioff].group = ngpm+ngroup+k;
	(coords->amass)[ioff] = atom_topol[ispec][j].mass*MCONV;
	(coords->qch)[ioff]   = atom_topol[ispec][j].charge*sqrt(CCONV)*simparms->scalecharge;
	(coords->alpha)[ioff]   = atom_topol[ispec][j].alpha;
	strcpy((simparms->group)[ngpm+ngroup+k],atom_topol[ispec][j].group);
      }
    }
    iatom += spec_now->napm*spec_now->nmol;
    ngroup += ngpm*spec_now->nmol;
  }
  /* reallocate atom_types to just the number needed */
  simparms->atom = (WORD *)realloc(simparms->atom,ntypes*sizeof(WORD));
  
#ifdef PRINT_ATOMS
  for(i=0;i<natoms;i++){
    printf("%d %d %d %d %d %d %s %s %s %lg %lg\n",
	   i,(simparms->atom_map)[i].species,
	   (simparms->atom_map)[i].molecule,(simparms->atom_map)[i].group,
	   (simparms->atom_map)[i].atom,(simparms->atom_map)[i].type,
	   (label->mole)[(simparms->atom_map)[i].species],
	   (label->atom)[(simparms->atom_map)[i].type],
	   (label->group)[(simparms->atom_map)[i].group],
	   (coords->amass)[i],(coords->qch)[i]);
  }
#endif

  inter->ntype = ntypes;
  if(simparms->rank==0){
    sprintf(line,"There are %d groups and %d types of atoms",ngroup,ntypes);
    md_stdout(line);
  }
  /* let the intermolecular force routine know the atom types */
  inter->itype = (int *)cmalloc(natoms*sizeof(int));
  for(i=0;i<natoms;i++){inter->itype[i] = (simparms->atom_map)[i].type;}

  /* search the data base for the interactions */
  search_ter_base(simparms,ntypes,simparms->atom,startup,inter,ngbr);
  
  /* get the types of intra interactions */

  search_bond_base(simparms,startup->bond_file,ibo,jbo,tbo,
		   atom_topol,startup->spec_root,startup->nspec,
		   &coords->bonds);
  search_bend_base(simparms,startup->bend_file,ibe,jbe,kbe,tbe,atom_topol,
		   startup->spec_root,startup->nspec,&coords->bends);
  search_tors_base(simparms,startup->tors_file,ito,jto,kto,lto,tto,atom_topol,
		   startup->spec_root,startup->nspec,&coords->tors);
  search_onfo_base(simparms,startup->onfo_file,simparms->icalc_type,
		   iof,jof,atom_topol,
		   startup->spec_root,startup->nspec,startup->scale_onfo,
		   startup->scale_onfo_e,inter->ntable,&coords->onfo);

  /* cross potential interactions */
  search_bocross_base(simparms,startup->bond_file,ibonx,jbonx,kbonx,tbonx,
		      atom_topol,startup->spec_root,startup->nspec,
		      &coords->bondxs);

  /* setup the charge charge interaction  */

  set_charges(simparms,coords,inter,ngbr,startup);

  /* setup exclusions */
  set_exclude(simparms,coords,inter);

  /* free temporary storage */
  for(i=0;i<startup->nspec;i++){
    free(atom_topol[i]);
    if(ibo[i]!=NULL){free(ibo[i]); free(jbo[i]); free(tbo[i]);}
    if(ibe[i]!=NULL){free(ibe[i]); free(jbe[i]); free(kbe[i]); free(tbe[i]);}
    if(ito[i]!=NULL){
      free(ito[i]); free(jto[i]); free(kto[i]); free(lto[i]); free(tto[i]);
    }
    if(iof[i]!=NULL){free(iof[i]); free(jof[i]);}
    if(ibonx[i]!=NULL){
      free(ibonx[i]);free(jbonx[i]);free(kbonx[i]);free(tbonx[i]); 
    }
  }
  free(atom_topol);
  free(ibo); free(jbo); free(tbo);
  free(ibe); free(jbe); free(kbe); free(tbe);
  free(ito); free(jto); free(kto); free(lto); free(tto);
  free(iof); free(jof);
  free(ibonx);free(jbonx);free(kbonx);free(tbonx);
}

/*------------------------------------------------------------------------*/
void read_top_file(SIMPARMS *simparms,STARTUP *startup,
		   ATOM_TOPOL **atom_topol,
		   int **ibo,int **jbo, WORD **tbo,
		   int **ibe,int **jbe,int **kbe, WORD **tbe,
		   int **ito,int **jto,int **kto,int **lto, WORD **tto,
		   int **iof,int **jof,
		   int **ibonx,int **jbonx,int **kbonx, WORD **tbonx)
{
  int i,j,ispec;
  
  FILE *fp;
  SPECIES *spec_now;
  
  simparms->natoms = 0;
  spec_now = &(startup->spec_root);

  for(ispec=0;ispec<startup->nspec;ispec++,spec_now=spec_now->next){
    
    /* open parameter file */
    fp = fopen(spec_now->filename,"r");
    if (fp == NULL){
      fprintf(stderr,"ERROR: can't open set file (%s)\n",spec_now->filename);
      exit(1);
    } else {
      if(simparms->rank==0){
	fprintf(stdout,"Reading in parameter from file %s\n",
		spec_now->filename);
      }
    }

    get_species(fp,spec_now);

    if(spec_now->napm <= 0)
      fprintf(stderr,"ERROR, %s has %d atoms!\n",
	      spec_now->filename,spec_now->napm);

    atom_topol[ispec] = (ATOM_TOPOL *)cmalloc(spec_now->napm*
					     sizeof(atom_topol[0][0]));
    get_atom(fp,spec_now->filename,spec_now->napm,atom_topol[ispec]);

    /* bonding */
    ibo[ispec] = jbo[ispec] = NULL;
    tbo[ispec] = NULL;
    if(spec_now->nbond>0) {
      ibo[ispec] = (int *)cmalloc(spec_now->nbond*sizeof(int));
      jbo[ispec] = (int *)cmalloc(spec_now->nbond*sizeof(int));
      tbo[ispec] = (WORD *)cmalloc(spec_now->nbond*sizeof(WORD));
      
      get_bond(fp,spec_now->filename,spec_now->nbond,
	       ibo[ispec],jbo[ispec],tbo[ispec]);
      for(i=0;i<spec_now->nbond;i++){
	for(j=0;j<spec_now->napm;j++){
	  if(atom_topol[ispec][j].atm_idx == ibo[ispec][i]) ibo[ispec][i] = j;
	  if(atom_topol[ispec][j].atm_idx == jbo[ispec][i]) jbo[ispec][i] = j;
	}
      }
    }

    /* bending */
    ibe[ispec] = jbe[ispec] = kbe[ispec] = NULL;
    tbe[ispec] = NULL;
    if(spec_now->nbend>0){
      ibe[ispec] = (int *)cmalloc(spec_now->nbend*sizeof(int));
      jbe[ispec] = (int *)cmalloc(spec_now->nbend*sizeof(int));
      kbe[ispec] = (int *)cmalloc(spec_now->nbend*sizeof(int));
      tbe[ispec] = (WORD *)cmalloc(spec_now->nbend*sizeof(WORD));
      get_bend(fp,spec_now->filename,spec_now->nbend,
	       ibe[ispec],jbe[ispec],kbe[ispec],tbe[ispec]);
      for(i=0;i<spec_now->nbend;i++){
	for(j=0;j<spec_now->napm;j++){
	  if(atom_topol[ispec][j].atm_idx == ibe[ispec][i]) ibe[ispec][i] = j;
	  if(atom_topol[ispec][j].atm_idx == jbe[ispec][i]) jbe[ispec][i] = j;
	  if(atom_topol[ispec][j].atm_idx == kbe[ispec][i]) kbe[ispec][i] = j;
	}
      }
    }

    /* torsions */
    ito[ispec] = jto[ispec] = kto[ispec] = lto[ispec] = NULL;
    tto[ispec] = NULL;
    if(spec_now->ntors>0) {
      ito[ispec] = (int *)cmalloc(spec_now->ntors*sizeof(int));
      jto[ispec] = (int *)cmalloc(spec_now->ntors*sizeof(int));
      kto[ispec] = (int *)cmalloc(spec_now->ntors*sizeof(int));
      lto[ispec] = (int *)cmalloc(spec_now->ntors*sizeof(int));
      tto[ispec] = (WORD *)cmalloc(spec_now->ntors*sizeof(WORD));
      get_tors(fp,spec_now->filename,spec_now->ntors,
	       ito[ispec],jto[ispec],kto[ispec],lto[ispec],tto[ispec]);
      for(i=0;i<spec_now->ntors;i++){
	for(j=0;j<spec_now->napm;j++){
	  if(atom_topol[ispec][j].atm_idx == ito[ispec][i]) ito[ispec][i] = j;
	  if(atom_topol[ispec][j].atm_idx == jto[ispec][i]) jto[ispec][i] = j;
	  if(atom_topol[ispec][j].atm_idx == kto[ispec][i]) kto[ispec][i] = j;
	  if(atom_topol[ispec][j].atm_idx == lto[ispec][i]) lto[ispec][i] = j;
	}
      }
    }
    /* onefour */
    iof[ispec] = jof[ispec] = NULL;
    if(spec_now->n14>0){
      iof[ispec] = (int *)cmalloc(spec_now->n14*sizeof(int));
      jof[ispec] = (int *)cmalloc(spec_now->n14*sizeof(int));
      get_onefour(fp,spec_now->filename,spec_now->n14,iof[ispec],jof[ispec]);
      for(i=0;i<spec_now->n14;i++){
	for(j=0;j<spec_now->napm;j++){
	  if(atom_topol[ispec][j].atm_idx == iof[ispec][i]) iof[ispec][i] = j;
	  if(atom_topol[ispec][j].atm_idx == jof[ispec][i]) jof[ispec][i] = j;
	}
      }
    }
    /* cross bonding */
    ibonx[ispec] = jbonx[ispec] = kbonx[ispec] = NULL;
    tbonx[ispec] = NULL;
    if(spec_now->nbondx>0){
      ibonx[ispec] = (int *)cmalloc(spec_now->nbondx*sizeof(int));
      jbonx[ispec] = (int *)cmalloc(spec_now->nbondx*sizeof(int));
      kbonx[ispec] = (int *)cmalloc(spec_now->nbondx*sizeof(int));
      tbonx[ispec] = (WORD *)cmalloc(spec_now->nbondx*sizeof(WORD));
      get_bondx(fp,spec_now->filename,spec_now->nbondx,
		ibonx[ispec],jbonx[ispec],kbonx[ispec],tbonx[ispec]);
      for(i=0;i<spec_now->nbondx;i++){
	for(j=0;j<spec_now->napm;j++){
	  if(atom_topol[ispec][j].atm_idx==ibonx[ispec][i]) ibonx[ispec][i]=j;
	  if(atom_topol[ispec][j].atm_idx==jbonx[ispec][i]) jbonx[ispec][i]=j;
	  if(atom_topol[ispec][j].atm_idx==kbonx[ispec][i]) kbonx[ispec][i]=j;
	}
      }
    }
    fclose(fp);
    
    simparms->natoms += spec_now->napm*spec_now->nmol;
  }
}
