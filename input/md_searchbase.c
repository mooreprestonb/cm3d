/* 
   Copyright (C) Dr. Preston B. Moore

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

/* #define WRITE_TABLE  */
/* #define DEBUG */

/*---------------------------------------------------------------*/
void  search_ter_base(SIMPARMS *simparms,STARTUP *startup,
   INTER *inter,NGBR *ngbr)
{
  char line[MAXLINELEN];
  int nmeta_key;
  int i,j,ni,nj,k,l,ninter,inull,ntable,ntypes,isize;
  double vshift,skin,*vlrc,vlrc_tot,*wlrc,wlrc_tot;
  double rmaxs,rmins,rmaxl,rminl;
  WORD *atypes;
  META_KEY *meta_key;
  FILE *fp;
#ifdef WRITE_TABLE
  double xt;
  WORD table;
#endif

  ntypes = simparms->ntypes;
  atypes = simparms->atom;
  skin = inter->skin_ter;

  meta_key = (META_KEY *)cmalloc(sizeof(*meta_key));
  /* open data base file */
  if((fp = fopen(startup->inter_file,"r")) == NULL){
    sprintf(line,"can't open data base (%s)",startup->inter_file);
    md_error(line);
  } else {
    if(simparms->rank==0){
      sprintf(line,"Searching intermolecular data base %s for interactions",
	      startup->inter_file);
      md_stdout(line);
    }
  }
  fclose(fp);
  parse_file(startup->inter_file,&nmeta_key,meta_key);

  /* allocate interaction map */

  inter->map = (int **)cmalloc(ntypes*sizeof(int*));
  inter->map[0] = (int *)cmalloc(ntypes*ntypes*sizeof(int));
  for(i=0;i<ntypes*ntypes;i++) inter->map[0][i] = 0;
  for(i=1;i<ntypes;i++) inter->map[i] = inter->map[i-1]+ntypes;
  ninter = ntypes*(ntypes+1)/2+1;
  
  inter->rmaxs_cut2  = (double *)cmalloc(ninter*sizeof(double));
  inter->rmins_cut2  = (double *)cmalloc(ninter*sizeof(double));
  inter->rmaxl_cut2  = (double *)cmalloc(ninter*sizeof(double));
  inter->rminl_cut2  = (double *)cmalloc(ninter*sizeof(double));
  inter->rmaxs_skin2 = (double *)cmalloc(ninter*sizeof(double));
  inter->rmaxl_skin2 = (double *)cmalloc(ninter*sizeof(double));
  inter->rminl_skin2 = (double *)cmalloc(ninter*sizeof(double));
  inter->dx2tab      = (double *)cmalloc(ninter*sizeof(double));
  vlrc               = (double *)cmalloc(ninter*sizeof(double));
  wlrc               = (double *)cmalloc(ninter*sizeof(double));

  inter->vtab    = (double **)cmalloc(ninter*sizeof(double *));
  inter->dvtab   = (double **)cmalloc(ninter*sizeof(double *));
  inter->d2vtab  = (double **)cmalloc(ninter*sizeof(double *));
  inter->switchs = (double **)cmalloc(ninter*sizeof(double *));

  /* allocate table interactions */
  ntable = inter->ntable*ninter;
  inter->vtab[0]    = (double *)cmalloc(ntable*sizeof(double));
  inter->dvtab[0]   = (double *)cmalloc(ntable*sizeof(double));
  inter->d2vtab[0]  = (double *)cmalloc(ntable*sizeof(double));
  inter->switchs[0] = (double *)cmalloc(ntable*sizeof(double));

  ntable = inter->ntable;
  for(i=1;i<ninter;i++){
    inter->vtab[i] = inter->vtab[0] + i*ntable;
    inter->dvtab[i] = inter->dvtab[0] + i*ntable;
    inter->d2vtab[i] = inter->d2vtab[0] + i*ntable;
    inter->switchs[i] = inter->switchs[0] + i*ntable;
  }

  k = 0;
  for(i=0;i<ntypes;i++){
    for(j=i;j<ntypes;j++){
      get_ter_parm(nmeta_key,meta_key,startup->inter_file,atypes[i],atypes[j],
		   &rminl,&rmaxl,&rmins,&rmaxs,&vlrc[k],&wlrc[k],
		   inter,k,&inull,simparms->scaleeps);
      if(inull){
        inter->map[i][j] = inter->map[j][i] = -1;	
      } else {
        inter->rmaxs_cut2[k]  = rmaxs*rmaxs;
        inter->rmins_cut2[k]  = rmins*rmins;
        inter->rmaxl_cut2[k]  = rmaxl*rmaxl;
        inter->rminl_cut2[k]  = rminl*rminl;
        inter->rmaxs_skin2[k] =(rmaxs+skin)*(rmaxs+skin);
        inter->rmaxl_skin2[k] =(rmaxl+skin)*(rmaxl+skin);
        inter->rminl_skin2[k] =(rminl-skin)*(rminl-skin);
        inter->map[i][j] = inter->map[j][i] = k;
        k++;
      }
    }
  }
  free_meta_keys(nmeta_key,&meta_key);

  if(k==0) {
    md_warning("NO interactions????");
    k=1;
  }
    
  ntable = inter->ntable*k*sizeof(double);
  /* reallocate memory (so we don't waste space) */
  simparms->mem_bytes += 5*ntable;

  /* reduce to the correct number of interaction excluding nulls */
  inter->ntype = ninter = k;
  isize = ninter*sizeof(double);
  inter->rmaxs_cut2 = (double *)crealloc(inter->rmaxs_cut2,isize);
  inter->rmins_cut2 = (double *)crealloc(inter->rmins_cut2,isize);
  inter->rmaxl_cut2 = (double *)crealloc(inter->rmaxl_cut2,isize);
  inter->rminl_cut2 = (double *)crealloc(inter->rminl_cut2,isize);
  inter->rmaxs_skin2= (double *)crealloc(inter->rmaxs_skin2,isize);
  inter->rmaxl_skin2= (double *)crealloc(inter->rmaxl_skin2,isize);
  inter->rminl_skin2= (double *)crealloc(inter->rminl_skin2,isize);
  inter->dx2tab     = (double *)crealloc(inter->dx2tab,isize);
  vlrc              = (double *)crealloc(vlrc,isize);
  wlrc              = (double *)crealloc(wlrc,isize);

  simparms->mem_bytes += 8*isize;

  isize = ninter*sizeof(double *);
  inter->vtab    = (double **)crealloc(inter->vtab,isize);
  inter->dvtab   = (double **)crealloc(inter->dvtab,isize);
  inter->d2vtab  = (double **)crealloc(inter->d2vtab,isize);
  inter->switchs = (double **)crealloc(inter->switchs,isize);

  simparms->mem_bytes += 4*isize;

  isize = inter->ntable*k*sizeof(double);
  /* reallocate memory (so we don't waste space) */
  inter->vtab[0]    = (double *)crealloc(inter->vtab[0],isize);
  inter->dvtab[0]   = (double *)crealloc(inter->dvtab[0],isize);
  inter->d2vtab[0]  = (double *)crealloc(inter->d2vtab[0],isize);
  inter->switchs[0] = (double *)crealloc(inter->switchs[0],isize);

  simparms->mem_bytes += 4*isize;

  ntable = inter->ntable;
  for(i=1;i<ninter;i++){
    inter->vtab[i] = inter->vtab[0] + i*ntable;
    inter->dvtab[i] = inter->dvtab[0] + i*ntable;
    inter->d2vtab[i] = inter->d2vtab[0] + i*ntable;
    inter->switchs[i] = inter->switchs[0] + i*ntable;
  }

  if(startup->ishift){
    for(k=0;k<ninter;k++){
      vshift=inter->vtab[k][inter->ntable-3];
      for(i=0;i<inter->ntable;i++){
        inter->vtab[k][i]-=vshift;
      }
    }
  }

#ifdef WRITE_TABLE
  for(i=0;i<ninter;i++){
    sprintf(table,"table%d.dat",i); 
    sprintf(line,"Saving interaction %d between types in file %s",i,table);
    md_stdout(line);
    fp = fopen(table,"w");
    for(k=0;k<inter->ntable;k++){
      xt = k/inter->dx2tab[i] + inter->rmins_cut2[i];
      fprintf(fp,"%d %g %g %g %g %g\n",k,sqrt(xt),inter->vtab[i][k],
	      inter->dvtab[i][k]*sqrt(xt),inter->d2vtab[i][k],
	      inter->switchs[i][k]);
    }
    fclose(fp);
  }
#endif

  inter->dxn = (double *)cmalloc(3*simparms->natoms*sizeof(double));
  inter->dyn = inter->dxn+  simparms->natoms;
  inter->dzn = inter->dxn+2*simparms->natoms;

  simparms->mem_bytes += 3*simparms->natoms*sizeof(double);

  wlrc_tot = vlrc_tot = 0.0;
  for(i=0;i<ntypes;i++){
    for(ni=0,j=0;j<simparms->natoms;j++) if(inter->itype[j] == i) ni++;
    for(j=0;j<ntypes;j++){
      for(nj=0,k=0;k<simparms->natoms;k++) if(inter->itype[k] == j) nj++;
      l = inter->map[i][j];
      if(l != -1)vlrc_tot += ni*nj*M_PI*2.*vlrc[l]; 
      if(l != -1)wlrc_tot += ni*nj*M_PI*2.*wlrc[l]/3.; 
    }
  }

#ifdef DEBUG
  sprintf(line,"ninter = %d\nvlrc_tot = %g",ninter,vlrc_tot);
  md_stdout(line);
#endif

  if(simparms->vlrc == -1){
    simparms->vlrc=vlrc_tot=0.;
    simparms->wlrc=wlrc_tot=0.;
  } else {
    if(simparms->iperd==0) vlrc_tot = wlrc_tot = 0.;
    simparms->vlrc=vlrc_tot;
    simparms->wlrc=wlrc_tot;
  }
  free(vlrc);
  free(wlrc);

  /* free d2vtab if second derivative are not needed */
  if(simparms->icalc_type == 0 || 
     (simparms->icalc_type == 2 && simparms->imin_type!=3)){
    free(inter->d2vtab[0]);
    free(inter->d2vtab);
    simparms->mem_bytes -= inter->ntable*sizeof(double);
    simparms->mem_bytes -= ninter*sizeof(double *);
  }
  if(simparms->rank==0){
    sprintf(line,"There were %d atom types and %d interations (out of %d)",
	    ntypes,ninter,ntypes*(ntypes+1)/2);
    md_stdout(line);
    sprintf(line,"Long range correction factor V = %g K A^3, W = %g K A^3",
	    vlrc_tot,wlrc_tot);
    md_stdout(line);
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
}

/*---------------------------------------------------------------*/
void search_bond_base(SIMPARMS *simparms,char *filename,
		      int **ibo,int **jbo,WORD **tbo,
		      ATOM_TOPOL **atom_topol,SPECIES spec_root,int nspec,
		      BONDS *bonds)
{
  char line[MAXLINELEN];
  int i,j,ispec,ib,ifix,nmeta_key;
  int ib_now,ib_off,ioff,nb_types,iflag,it_now,iatom;
  double eq_now,fk_now,bm_now;
  SPECIES *spec_now;
  WORD *atm1,*atm2,typ1,typ2,*bond_type,btyp,tword;
  META_KEY *meta_key;
  FILE *fp;

  meta_key = (META_KEY *)cmalloc(sizeof(*meta_key));
  spec_now = &spec_root;
  /* allocate some memory for the static arrays */
  bonds->nbonds = 0;
  bonds->idxbond = (int *)cmalloc(sizeof(int));
  bonds->ibond = (int *)cmalloc(sizeof(int));
  bonds->jbond = (int *)cmalloc(sizeof(int));
  atm1 = (WORD *)cmalloc(sizeof(WORD));
  atm2 = (WORD *)cmalloc(sizeof(WORD));
  bond_type = (WORD *)cmalloc(sizeof(WORD));

  nmeta_key = nb_types = iatom =  0;
  ib_now = ib_off = 0;
  
  /* loop over the species and determine bonded types */
  for(ispec=0,spec_now=&spec_root;ispec<nspec;
      ispec++,spec_now=spec_now->next){
    
    ib = (ib_now+spec_now->nbond*spec_now->nmol);

    if(ib){
      bonds->idxbond =(int *)crealloc(bonds->idxbond,ib*sizeof(int));
      bonds->ibond   =(int *)crealloc(bonds->ibond,ib*sizeof(int));
      bonds->jbond   =(int *)crealloc(bonds->jbond,ib*sizeof(int));
    }
    for(ib=0;ib<spec_now->nbond;ib++){
      /* get the actual types of atoms into strings */
      strcpy(typ1,atom_topol[ispec][ibo[ispec][ib]].type);
      strcpy(typ2,atom_topol[ispec][jbo[ispec][ib]].type);
      strcpy(btyp,tbo[ispec][ib]);

      /* put in some order for ease of comparison */
      if(strcasecmp(typ1,typ2)>0){
        strcpy(tword,typ1);
        strcpy(typ1,typ2);
        strcpy(typ2,tword);
      }
      
      /* check to see if we have a new bonding interaction */
      iflag = 1;
      for(i=0;i<nb_types;i++){
        if(!strcasecmp(typ1,atm1[i]) &&  !strcasecmp(typ2,atm2[i]) &&
           !strcasecmp(btyp,bond_type[i])){
          iflag = 0;
          break;
        }
      }

      /* fill up bond lists */
      for(j=0;j<spec_now->nmol;j++){
        ioff = iatom + j*spec_now->napm;
        bonds->idxbond[ib_now] = i;
        bonds->ibond[ib_now] = ioff + ibo[ispec][ib];
        bonds->jbond[ib_now] = ioff + jbo[ispec][ib];
        ib_now++;
      }
      /* if we found a new type -- include it */
      if(iflag){
        nb_types++;
        atm1 = (WORD *)crealloc(atm1,nb_types*sizeof(WORD));
        atm2 = (WORD *)crealloc(atm2,nb_types*sizeof(WORD));
        bond_type = (WORD *)crealloc(bond_type,nb_types*sizeof(WORD));
        strcpy(atm1[nb_types-1],typ1);
        strcpy(atm2[nb_types-1],typ2);
        strcpy(bond_type[nb_types-1],btyp);
      }
    }
    if(ib_now != (ib_off += spec_now->nbond*spec_now->nmol)){
      sprintf(line,"in set up bonding routine %d %d",ib_off,ib_now);
      md_error(line);
    }
    iatom += spec_now->napm*spec_now->nmol;
  }
  bonds->nbonds = ib_off;
  bonds->itypes = nb_types;

  if(bonds->nbonds>0){
    bonds->eqbond = (double *)cmalloc(bonds->nbonds*sizeof(double));
    bonds->fkbond = (double *)cmalloc(bonds->nbonds*sizeof(double));
    bonds->bmorse = (double *)cmalloc(bonds->nbonds*sizeof(double));
    bonds->itypbond = (int *)cmalloc(bonds->nbonds*sizeof(int));
    bonds->ifix = (int *)cmalloc(bonds->nbonds*sizeof(int));
    bonds->aq  = (double *)cmalloc(NP_BOND*bonds->nbonds*sizeof(double));
    simparms->mem_bytes += ((3+NP_BOND)*bonds->nbonds*sizeof(double) + 
			    2*bonds->nbonds*sizeof(int));
    
    
    /* open data base file */
    if ( (fp = fopen(filename,"r")) == NULL){
      sprintf(line,"can't open bond data base (%s)",filename);
      md_error(line);
    } else {
      if(simparms->rank==0){
        sprintf(line,"Searching data base %s for bonded interactions",
           filename);
        md_stdout(line);
      }
    }
    fclose(fp);
    parse_file(filename,&nmeta_key,meta_key);

    for(i=0;i<nb_types;i++){
      get_bondparm(nmeta_key,meta_key,filename,atm1[i],atm2[i],bond_type[i],
		   &it_now,&eq_now,&fk_now,&bm_now,&ifix,
		   &(bonds->aq[i*NP_BOND]));
      bonds->itypbond[i] = it_now;  
      bonds->eqbond[i] = eq_now;  
      bonds->fkbond[i] = fk_now;
      bonds->bmorse[i] = bm_now;
      bonds->ifix[i] = ifix;
      /* aq already set */
    }
  }
  if(nmeta_key>0) free_meta_keys(nmeta_key,&meta_key);
  if(simparms->rank==0){
    sprintf(line,"There were %d bonds types and %d bonds",
	    nb_types,bonds->nbonds);
    md_stdout(line);
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
#ifdef DEBUG
  for(i=0;i<nb_types;i++){
    int ii = i*NP_BOND;
    sprintf(line,"%d %s %s %s %g %g %d aqs = %g %g %g %g %g",i,atm1[i],atm2[i],
	    bond_type[i],bonds->eqbond[i],bonds->fkbond[i],bonds->ifix[i],
	    bonds->aq[ii],bonds->aq[ii+1],bond->aq[ii+2],bond->aq[ii+3],
	    bond->aq[ii+4]));
    md_stdout(line);
  }
#endif
  free(atm1);  free(atm2);  free(bond_type);
}
/*---------------------------------------------------------------*/
/* \brief search data base frile for bend interactions
 *
 */
void search_bend_base(SIMPARMS *simparms,char *filename,
		      int **ibe,int **jbe,int **kbe,int **tybe,WORD **tbe,
		      ATOM_TOPOL **atom_topol,SPECIES spec_root,int nspec,
		      BENDS *bends)
{
  char line[MAXLINELEN];
  int i,j,ispec,ib,ifix,it_now;
  int ib_now,ib_off,ioff,nb_types,iflag,iatom;
  int *iexamp,*jexamp,*kexamp,*tyexamp;
  SPECIES *spec_now;
  WORD *atm1,*atm2,*atm3,typ1,typ2,typ3,tword;
  WORD *bend_type,btyp;
  FILE *fp;
  double eq_now,fk_now;

  spec_now = &spec_root;
  /* allocate some memory for the static arrays */
  bends->nbends = 0;
  bends->idxbend = (int *)cmalloc(sizeof(int));
  bends->ibend = (int *)cmalloc(sizeof(int));
  bends->jbend = (int *)cmalloc(sizeof(int));
  bends->kbend = (int *)cmalloc(sizeof(int));
  bends->itypbend = (int *)cmalloc(sizeof(int));
  iexamp = (int *)cmalloc(sizeof(int));
  jexamp = (int *)cmalloc(sizeof(int));
  kexamp = (int *)cmalloc(sizeof(int));
  tyexamp = (int *)cmalloc(sizeof(int));
  atm1 = (WORD *)cmalloc(sizeof(WORD));
  atm2 = (WORD *)cmalloc(sizeof(WORD));
  atm3 = (WORD *)cmalloc(sizeof(WORD));
  bend_type = (WORD *)cmalloc(sizeof(WORD));
  
  nb_types = iatom =  0;
  ib_now = ib_off = 0;
  
  /* loop over the species and determine bended types */
  for(ispec=0,spec_now=&spec_root;ispec<nspec;
      ispec++,spec_now=spec_now->next){
    
    ib = (ib_now+spec_now->nbend*spec_now->nmol);
    if(ib!=0){
      bends->idxbend =(int *)crealloc(bends->idxbend,ib*sizeof(int));
      bends->ibend   =(int *)crealloc(bends->ibend,ib*sizeof(int));
      bends->jbend   =(int *)crealloc(bends->jbend,ib*sizeof(int));
      bends->kbend   =(int *)crealloc(bends->kbend,ib*sizeof(int));
    }

    for(ib=0;ib<spec_now->nbend;ib++){
      /* get the actual types of atoms into strings */
      strcpy(typ1,atom_topol[ispec][ibe[ispec][ib]].type);
      strcpy(typ2,atom_topol[ispec][jbe[ispec][ib]].type);
      strcpy(typ3,atom_topol[ispec][kbe[ispec][ib]].type);
      strcpy(btyp,tbe[ispec][ib]);
      
      /* put in some order for ease of comparison */
      if(strcasecmp(typ1,typ3)>0){
        strcpy(tword,typ1);strcpy(typ1,typ3);strcpy(typ3,tword);
      }
      
      /* check to see if we have a new bending interaction */
      iflag = 1;
      for(i=0;i<nb_types;i++){
        if(!strcasecmp(typ1,atm1[i]) && !strcasecmp(typ2,atm2[i]) &&
           !strcasecmp(typ3,atm3[i]) && !strcasecmp(btyp,bend_type[i])){
          iflag = 0;
          break;
        }
      }

      /* fill up bend lists */
      for(j=0;j<spec_now->nmol;j++){
        ioff = iatom + j*spec_now->napm;
        bends->idxbend[ib_now] = i;
        bends->ibend[ib_now] = ioff + ibe[ispec][ib];
        bends->jbend[ib_now] = ioff + jbe[ispec][ib];
        bends->kbend[ib_now] = ioff + kbe[ispec][ib];
        // bends->ibendtype[ib_now] = tybe[ispec][ib];
        ib_now++;	  
      }
      /* if we found a new type -- include it */
      if(iflag){
        nb_types++;
        atm1 = (WORD *)crealloc(atm1,nb_types*sizeof(WORD));
        atm2 = (WORD *)crealloc(atm2,nb_types*sizeof(WORD));
        atm3 = (WORD *)crealloc(atm3,nb_types*sizeof(WORD));
        bend_type = (WORD *)crealloc(bend_type,nb_types*sizeof(WORD));
        iexamp = (int *)crealloc(iexamp,nb_types*sizeof(int));
        jexamp = (int *)crealloc(jexamp,nb_types*sizeof(int));
        kexamp = (int *)crealloc(kexamp,nb_types*sizeof(int));
        tyexamp = (int *)crealloc(tyexamp,nb_types*sizeof(int));
        strcpy(atm1[nb_types-1],typ1);
        strcpy(atm2[nb_types-1],typ2);
        strcpy(atm3[nb_types-1],typ3);
        strcpy(bend_type[nb_types-1],btyp);
        iexamp[nb_types-1]=ibe[ispec][ib];
        jexamp[nb_types-1]=jbe[ispec][ib];
        kexamp[nb_types-1]=kbe[ispec][ib];
        tyexamp[nb_types-1]=tybe[ispec][ib];
      }
    }
    if(ib_now != (ib_off += spec_now->nbend*spec_now->nmol)){
      sprintf(line,"in set up bending routine %d %d",ib_off,ib_now);
      md_error(line);
    }
    iatom += spec_now->napm*spec_now->nmol;
  }
  bends->nbends = ib_off;
  bends->itypes = nb_types;

  if(bends->nbends>0){
    bends->eqbend = (double *)cmalloc(nb_types*sizeof(double));
    bends->fkbend = (double *)cmalloc(nb_types*sizeof(double));
    bends->itypbend = (int *)cmalloc(nb_types*sizeof(int));
    bends->ifix = (int *)cmalloc(nb_types*sizeof(int));
    bends->aq = (double *)cmalloc(NP_BEND*nb_types*sizeof(double));
    simparms->mem_bytes += ((2+NP_BEND)*nb_types*sizeof(double) + 
			    2*nb_types*sizeof(int));

    /* open data base file */
    if ( (fp = fopen(filename,"r")) == NULL){
      sprintf(line,"can't open bend data base (%s)",filename);
      md_error(line);
    } else {
      if(simparms->rank==0){
        sprintf(line,"Searching data base %s for bended interactions",
           filename);
        md_stdout(line);
      }
    }
    
    for(i=0;i<nb_types;i++){
      get_bendparm(fp,filename,atm1[i],atm2[i],atm3[i],&it_now,bend_type[i],
		   &eq_now,&fk_now,&ifix,iexamp[i],jexamp[i],kexamp[i],
		   &(bends->aq[i*NP_BEND]));
      bends->eqbend[i] = eq_now;   
      bends->fkbend[i] = fk_now;
      bends->itypbend[i] = it_now;
      bends->ifix[i]   = ifix;
      /* bends->aq[i] already set */
    }
    fclose(fp);
  }
  if(simparms->rank==0){
    sprintf(line,"There were %d types and %d bends",nb_types,bends->nbends);
    md_stdout(line);
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
#ifdef DEBUG
  for(i=0;i<nb_types;i++){
    double *aqt = bends-aq[i];
    sprintf(line,"%d %s %s %s label:%s eq:%g %g typ:%d %d a0:%g %g %g %g %g",
	    i,atm1[i],atm2[i],atm3[i],
	    bend_type[i],bends->eqbend[i],bends->fkbend[i],
	    bends->itypbend[i], bends->ifix[i],
	    aqt[0],aqt[1],aqt[2],aqt[3],aqt[4]);
    md_stdout(line);
  }
#endif
  /* switch to radians */
  for(i=0;i<nb_types;i++){bends->eqbend[i] *= M_PI/180.0;}
  /* free temp arrays */
  free(atm1); free(atm2); free(atm3); free(bend_type);
  free(iexamp);  free(jexamp);  free(kexamp);  free(tyexamp);
}
/*---------------------------------------------------------------*/
void search_tors_base(SIMPARMS *simparms,char *filename,
		      int **ito,int **jto,int **kto,int **lto,
		      WORD **tto,ATOM_TOPOL **atom_topol,
		      SPECIES spec_root,int nspec,TORS *tors)
{
  char line[MAXLINELEN];
  int i,j,ispec,it,ifix;
  int it_now,it_off,ioff,nt_types,iflag,iatom;
  int num_harm,num_pow,int_pot;
  SPECIES *spec_now;
  WORD *atm1,*atm2,*atm3,*atm4,*ttors;
  WORD ttyp,typ1,typ2,typ3,typ4;
  FILE *fp;
  double eq_now,fk_now;
  VPHI vphi_now;

  spec_now = &spec_root;
  /* allocate some memory for the static arrays */
  tors->ntorss = 0;
  tors->idxtors = (int *)cmalloc(sizeof(int));
  tors->itors = (int *)cmalloc(sizeof(int));
  tors->jtors = (int *)cmalloc(sizeof(int));
  tors->ktors = (int *)cmalloc(sizeof(int));
  tors->ltors = (int *)cmalloc(sizeof(int));
  atm1 = (WORD *)cmalloc(sizeof(WORD));
  atm2 = (WORD *)cmalloc(sizeof(WORD));
  atm3 = (WORD *)cmalloc(sizeof(WORD));
  atm4 = (WORD *)cmalloc(sizeof(WORD));
  ttors = (WORD *)cmalloc(sizeof(WORD));
  
  nt_types = iatom =  0;
  it_now = it_off = 0;
  
  /* loop over the species and determine torsed types */
  for(ispec=0,spec_now=&spec_root;ispec<nspec;
      ispec++,spec_now=spec_now->next){
    
    it = (it_now+spec_now->ntors*spec_now->nmol);
    if(it!=0){
       tors->idxtors =(int *)crealloc(tors->idxtors,it*sizeof(int));
       tors->itors   =(int *)crealloc(tors->itors,it*sizeof(int));
       tors->jtors   =(int *)crealloc(tors->jtors,it*sizeof(int));
       tors->ktors   =(int *)crealloc(tors->ktors,it*sizeof(int));
       tors->ltors   =(int *)crealloc(tors->ltors,it*sizeof(int));
     }
    
    for(it=0;it<spec_now->ntors;it++){
      /* get the actual types of atoms into strings */
      strcpy(typ1,atom_topol[ispec][ito[ispec][it]].type);
      strcpy(typ2,atom_topol[ispec][jto[ispec][it]].type);
      strcpy(typ3,atom_topol[ispec][kto[ispec][it]].type);
      strcpy(typ4,atom_topol[ispec][lto[ispec][it]].type);
      strcpy(ttyp,tto[ispec][it]);
      
      /* check to see if we have a new torsing interaction */
      iflag = 1;
      for(i=0;i<nt_types;i++){
	if(!strcasecmp(typ1,atm1[i]) && !strcasecmp(typ2,atm2[i]) &&
	   !strcasecmp(typ3,atm3[i]) && !strcasecmp(typ4,atm4[i]) &&
	   !strcasecmp(ttyp,ttors[i])) {
	  iflag = 0;
	  break;
	}
	if(!strcasecmp(typ4,atm1[i]) && !strcasecmp(typ3,atm2[i]) &&
	   !strcasecmp(typ2,atm3[i]) && !strcasecmp(typ1,atm4[i]) &&
	   !strcasecmp(ttyp,ttors[i])) {
	  iflag = 0;
	  break;
	}
      }

      /* fill up tors lists */
      for(j=0;j<spec_now->nmol;j++){
	ioff = iatom + j*spec_now->napm;
	tors->idxtors[it_now] = i;
	tors->itors[it_now] = ioff + ito[ispec][it];
	tors->jtors[it_now] = ioff + jto[ispec][it];
	tors->ktors[it_now] = ioff + kto[ispec][it];
	tors->ltors[it_now] = ioff + lto[ispec][it];
	it_now++;	  
      }
	/* if we found a new type -- include it */
      if(iflag){
	nt_types++;
	atm1 = (WORD *)crealloc(atm1,nt_types*sizeof(WORD));
	atm2 = (WORD *)crealloc(atm2,nt_types*sizeof(WORD));
	atm3 = (WORD *)crealloc(atm3,nt_types*sizeof(WORD));
	atm4 = (WORD *)crealloc(atm4,nt_types*sizeof(WORD));
	ttors= (WORD *)crealloc(ttors,nt_types*sizeof(WORD));
	strcpy(atm1[nt_types-1],typ1);
	strcpy(atm2[nt_types-1],typ2);
	strcpy(atm3[nt_types-1],typ3);
	strcpy(atm4[nt_types-1],typ4);
	strcpy(ttors[nt_types-1],ttyp);
      }
    }
    if(it_now != (it_off += spec_now->ntors*spec_now->nmol)){
      sprintf(line,"in set up torsing routine %d %d",it_off,it_now);
      md_error(line);
    }
    iatom += spec_now->napm*spec_now->nmol;
  }
#ifdef DEBUG
  sprintf(line,"number of torsions types = %d",nt_types);
  md_stdout(line);
  for(i=0;i<nt_types;i++){
    sprintf(line,"%d %s %s %s %s {%s}",i,atm1[i],atm2[i],atm3[i],atm4[i],
	    ttors[i]);
    md_stdout(line);
  }
#endif

  tors->ntorss = it_off;
  tors->itypes = nt_types;
  num_harm = 0;
  num_pow  = 0;

  if(tors->ntorss>0){
    tors->ind_pot = (int *)cmalloc(tors->ntorss*sizeof(int));
    tors->ifix = (int *)cmalloc(tors->ntorss*sizeof(int));
    tors->ind_pot_num = (int *)cmalloc(tors->ntorss*sizeof(int));
    
    simparms->mem_bytes += 3*tors->ntorss*sizeof(int);
    
    tors->eqtors = (double *)cmalloc(sizeof(double));
    tors->fktors = (double *)cmalloc(sizeof(double));
    tors->vphi   = (VPHI *)cmalloc(sizeof(VPHI));
    
    /* open data base file */
    if ( (fp = fopen(filename,"r")) == NULL){
      sprintf(line,"can't open torsion data base (%s)",filename);
      md_error(line);
    } else {
      if(simparms->rank==0){
	sprintf(line,"Searching data base %s for torsions interactions",
		filename);
	md_stdout(line);
      }
    }
    
    for(i=0;i<nt_types;i++){
      get_torsparm(fp,filename,atm1[i],atm2[i],atm3[i],atm4[i],ttors[i],
		   &int_pot,&eq_now,&fk_now,vphi_now,&ifix);
      tors->ind_pot[i] = int_pot;
      tors->ifix[i] = ifix;
      switch(int_pot){
      case 0: /* harmonic */
	num_harm++;
	tors->eqtors = (double *)crealloc(tors->eqtors,num_harm*sizeof(double));
	tors->fktors = (double *)crealloc(tors->fktors,num_harm*sizeof(double));
	tors->eqtors[num_harm-1] = eq_now;
	tors->fktors[num_harm-1] = fk_now;
	tors->ind_pot_num[i] = num_harm-1;
	break;
      case 1: /* cosine series (both power and w/ phase */
      case 2:
	num_pow++;
	tors->vphi = (VPHI *)crealloc(tors->vphi,num_pow*sizeof(VPHI));
	for(j=0;j<MAXPOWER_TORS;j++) tors->vphi[num_pow-1][j] = vphi_now[j];
	tors->ind_pot_num[i] = num_pow-1;
	break;
      default:
	sprintf(line,"while setting up torsion-pot_type %d incorrect",int_pot);
	md_error(line);
	break;
      }
    }
    fclose(fp);
  }
  simparms->mem_bytes += 2*num_harm*sizeof(double);
  simparms->mem_bytes += num_pow*sizeof(VPHI);
  
  if(simparms->rank==0){
    sprintf(line,"Number of torss types = %d, Number of total torss = %d",
	    nt_types,tors->ntorss);
    md_stdout(line);
    sprintf(line,"Number of harmonic = %d, Number of power series = %d",
	    num_harm,num_pow);
    md_stdout(line);
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
#ifdef DEBUG
  for(i=0;i<nt_types;i++){
    sprintf(line,"%d %s %s %s %s %d",i,
	    atm1[i],atm2[i],atm3[i],atm4[i],tors->ind_pot[i]);
    md_stdout(line);
    switch(tors->ind_pot[i]){
    case 0:
      j = tors->ind_pot_num[i];
      sprintf(line,"%d %g %g",j,tors->eqtors[j],tors->fktors[j]);
      md_stdout(line);
      break;
    case 1:
      int_pot = tors->ind_pot_num[i];
      sprintf(line,"%d ",int_pot);
      for(j=0;j<MAXPOWER_TORS;j++) sprintf(line,"%s %g ",
					   line,tors->vphi[int_pot][j]);
      md_stdout(line);
      break;
    }
  }
#endif

  /* switch to radians */
  for(i=0;i<num_harm;i++){tors->eqtors[i] *= M_PI/180.0;}

  free(atm1); free(atm2); free(atm3); free(atm4); free(ttors);
}

/*------------------------------------------------------------------------*/
void search_onfo_base(SIMPARMS *simparms,char *filename,int icalc_type,
		      int **iof,int **jof,
		      ATOM_TOPOL **atom_topol,SPECIES spec_root,int nspec,
		      double scale_onfo,double scale_onfo_e,int ntab,
		      ONFO *onfo)
{
  char line[MAXLINELEN];
  int i,j,ispec,it;
  int it_now,it_off,ioff,iatom,n14_types,iflag;
  WORD *atm1,*atm2,typ1,typ2,tword;
  SPECIES *spec_now;
  FILE *fp;

  onfo->scale_onefour = scale_onfo;
  onfo->scale_onefour_e = scale_onfo_e;
  onfo->n14table = ntab;

  iatom = it_now = it_off = it = n14_types = 0;
  onfo->i14 =(int *)cmalloc(sizeof(int));
  onfo->j14 =(int *)cmalloc(sizeof(int));
  onfo->i14dx = (int *)cmalloc(sizeof(int));
  atm1 = (WORD*)cmalloc(sizeof(WORD));
  atm2 = (WORD*)cmalloc(sizeof(WORD));

  /* loop over the species and determine all one-four interactions */
  for(ispec=0,spec_now=&spec_root;ispec<nspec;
      ispec++,spec_now=spec_now->next){
    
    it = (it_now+spec_now->n14*spec_now->nmol);
    
    if(it!=0){
      onfo->i14   =(int *)crealloc(onfo->i14,it*sizeof(int));
      onfo->j14   =(int *)crealloc(onfo->j14,it*sizeof(int));
      onfo->i14dx = (int *)crealloc(onfo->i14dx,it*sizeof(int));
    }

    /* fill up 1-4 lists */
    for(it=0;it<spec_now->n14;it++){

      /* get the actual types of atoms into strings */
      strcpy(typ1,atom_topol[ispec][iof[ispec][it]].type);
      strcpy(typ2,atom_topol[ispec][jof[ispec][it]].type);
      /* put in some order for ease of comparison */
      if(strcasecmp(typ1,typ2)>0){
	strcpy(tword,typ1);strcpy(typ1,typ2);strcpy(typ2,tword);
      }
      
      /* check to see if we have a new 14 interaction */
      iflag = 1;
      for(i=0;i<n14_types;i++){
	if(!strcasecmp(typ1,atm1[i]) &&  !strcasecmp(typ2,atm2[i])){
	  iflag = 0;
	  break;
	}
      }
      /* fill up the interaction list */
      for(j=0;j<spec_now->nmol;j++){
	ioff = iatom + j*spec_now->napm;
	onfo->i14dx[it_now] = i;
	onfo->i14[it_now] = ioff + iof[ispec][it];
	onfo->j14[it_now] = ioff + jof[ispec][it];
	it_now++;
      }
      /* if we found a new type -- include it */
      if(iflag){
	n14_types++;
	atm1 = (WORD *)crealloc(atm1,n14_types*sizeof(WORD));
	atm2 = (WORD *)crealloc(atm2,n14_types*sizeof(WORD));
	strcpy(atm1[n14_types-1],typ1);
	strcpy(atm2[n14_types-1],typ2);
      }
    }
    if(it_now != (it_off += spec_now->n14*spec_now->nmol)){
      sprintf(line,"in set up one-four routine %d %d",it_off,it_now);
      md_error(line);
    }
    iatom += spec_now->napm*spec_now->nmol;
  }
  onfo->n14s = it_off;
  onfo->itypes = n14_types;
  
  if(simparms->rank==0){
    sprintf(line,"# of one-fours = %d, and n14_types = %d",it_now,n14_types);
    md_stdout(line);
    sprintf(line,"scale_onefour = %g and scale_onefour_e = %g",
	    onfo->scale_onefour,onfo->scale_onefour_e);
    md_stdout(line);
  }
#ifdef DEBUG
  for(i=0;i<it_now;i++){
    sprintf(line,"%d %d %d %d",i,onfo->i14[i],onfo->j14[i],onfo->i14dx[i]);
    md_stdout(line);
  }
#endif

  if(onfo->n14s>0){
    /* open data base file */
    if ( (fp = fopen(filename,"r")) == NULL){
      sprintf(line,"can't open data base (%s)",filename);
      md_error(line);
    } else {
      if(simparms->rank==0){
	sprintf(line,"Searching data base %s for one-four interactions",
		filename);
	md_stdout(line);
      }
    }
    
    onfo->r14max2  = (double *)cmalloc(n14_types*sizeof(double));
    onfo->r14min2  = (double *)cmalloc(n14_types*sizeof(double));
    onfo->dx214tab = (double *)cmalloc(n14_types*sizeof(double));
    
    onfo->v14tab = dmatrix(0,n14_types,0,onfo->n14table-1);
    onfo->dv14tab = dmatrix(0,n14_types,0,onfo->n14table-1);
    onfo->d2v14tab = dmatrix(0,n14_types,0,onfo->n14table-1);
    
    for(i=0;i<n14_types;i++){
      get_onfo_parm(fp,filename,atm1[i],atm2[i],
		    &onfo->r14min2[i],&onfo->r14max2[i],&onfo->dx214tab[i],
		    onfo->v14tab[i],onfo->dv14tab[i],
		    onfo->d2v14tab[i],onfo->n14table);
      onfo->r14min2[i] = onfo->r14min2[i]*onfo->r14min2[i];
      onfo->r14max2[i] = onfo->r14max2[i]*onfo->r14max2[i];
    }
    fclose(fp);

    /* free second derivative if not needed during calc type */
    if(icalc_type == 0 || icalc_type == 2){
      free_dmatrix(onfo->d2v14tab,0,n14_types,0,onfo->n14table-1);
    }
  }
  simparms->mem_bytes += 3*n14_types*sizeof(double *);
  if(simparms->rank==0){
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }

  free(atm1);free(atm2);
}
/*------------------------------------------------------------------------*/
void search_bocross_base(SIMPARMS *simparms,char *filename,
			 int **ibonx,int **jbonx,int **kbonx,
			 WORD **tbonx,ATOM_TOPOL **atom_topol,
			 SPECIES spec_root,int nspec,BONDXS *bondxs)
{
  char line[MAXLINELEN];
  int i,j,ispec,ib,itemp,ifix;
  int ib_now,ib_off,ioff,nb_types,iflag,iatom;
  SPECIES *spec_now;
  WORD *atm1,*atm2,*atm3,typ1,typ2,typ3,*box_type,btyp,tword;
  FILE *fp;
  double eq1_now,eq2_now,fk_now;

  spec_now = &spec_root;
  /* allocate some memory for the static arrays */
  bondxs->nbondxs = 0;
  bondxs->idxbox = (int *)cmalloc(sizeof(int));
  bondxs->ibondx = (int *)cmalloc(sizeof(int));
  bondxs->jbondx = (int *)cmalloc(sizeof(int));
  bondxs->kbondx = (int *)cmalloc(sizeof(int));
  atm1 = (WORD *)cmalloc(sizeof(WORD));
  atm2 = (WORD *)cmalloc(sizeof(WORD));
  atm3 = (WORD *)cmalloc(sizeof(WORD));
  box_type = (WORD *)cmalloc(sizeof(WORD));
  
  nb_types = iatom =  0;
  ib_now = ib_off = 0;
  
  /* loop over the species and determine bonded types */
  for(ispec=0,spec_now=&spec_root;ispec<nspec;
      ispec++,spec_now=spec_now->next){
    
    ib = (ib_now+spec_now->nbondx*spec_now->nmol);
    if(ib){
      bondxs->idxbox =(int *)crealloc(bondxs->idxbox,ib*sizeof(int));
      bondxs->ibondx =(int *)crealloc(bondxs->ibondx,ib*sizeof(int));
      bondxs->jbondx =(int *)crealloc(bondxs->jbondx,ib*sizeof(int));
      bondxs->kbondx =(int *)crealloc(bondxs->kbondx,ib*sizeof(int));
    }
    for(ib=0;ib<spec_now->nbondx;ib++){
      /* get the actual types of atoms into strings */
      strcpy(typ1,atom_topol[ispec][ibonx[ispec][ib]].type);
      strcpy(typ2,atom_topol[ispec][jbonx[ispec][ib]].type);
      strcpy(typ3,atom_topol[ispec][kbonx[ispec][ib]].type);
      strcpy(btyp,tbonx[ispec][ib]);
      
      /* put in some order for ease of comparison */
      if(strcasecmp(typ1,typ3)>0){
        strcpy(tword,typ1);strcpy(typ1,typ3);strcpy(typ3,tword);
        itemp = ibonx[ispec][ib];
        ibonx[ispec][ib] = kbonx[ispec][ib];
        kbonx[ispec][ib] = itemp;
      }
      
      /* check to see if we have a new bonding cross term interaction */
      iflag = 1;
      for(i=0;i<nb_types;i++){
	if(!strcasecmp(typ1,atm1[i]) &&  !strcasecmp(typ2,atm2[i]) &&
	   !strcasecmp(typ3,atm3[i]) && !strcasecmp(btyp,box_type[i])){
	  iflag = 0;
	  break;
	}
      }

      /* fill up bond lists */
      for(j=0;j<spec_now->nmol;j++){
	ioff = iatom + j*spec_now->napm;
	bondxs->idxbox[ib_now] = i;
	bondxs->ibondx[ib_now] = ioff + ibonx[ispec][ib];
	bondxs->jbondx[ib_now] = ioff + jbonx[ispec][ib];
	bondxs->kbondx[ib_now] = ioff + kbonx[ispec][ib];
	ib_now++;
      }
      /* if we found a new type -- include it */
      if(iflag){
	nb_types++;
	atm1 = (WORD *)crealloc(atm1,nb_types*sizeof(WORD));
	atm2 = (WORD *)crealloc(atm2,nb_types*sizeof(WORD));
	atm3 = (WORD *)crealloc(atm3,nb_types*sizeof(WORD));
	box_type = (WORD *)crealloc(box_type,nb_types*sizeof(WORD));
	strcpy(atm1[nb_types-1],typ1);
	strcpy(atm2[nb_types-1],typ2);
	strcpy(atm3[nb_types-1],typ3);
	strcpy(box_type[nb_types-1],btyp);
      }
    }
    if(ib_now != (ib_off += spec_now->nbondx*spec_now->nmol)){
      sprintf(line,"in set up cross bonding routine %d %d",ib_off,ib_now);
      md_error(line);
    }
    iatom += spec_now->napm*spec_now->nmol;
  }
  bondxs->nbondxs = ib_off;
  bondxs->itypes = nb_types;

  if(bondxs->nbondxs>0){
    bondxs->eq1bondx = (double *)cmalloc(bondxs->nbondxs*sizeof(double));
    bondxs->eq2bondx = (double *)cmalloc(bondxs->nbondxs*sizeof(double));
    bondxs->fkbondx = (double *)cmalloc(bondxs->nbondxs*sizeof(double));
    bondxs->ifix = (int *)cmalloc(bondxs->nbondxs*sizeof(int));
    
    /* open data base file */
    if((fp = fopen(filename,"r")) == NULL){
      sprintf(line,"can't open data base (%s)",filename);
      md_error(line);
    } else {
      if(simparms->rank==0){
	sprintf(line,"Searching data base %s for cross bonded interactions",
		filename);
	md_stdout(line);
      }
    }
    
    for(i=0;i<nb_types;i++){
      get_bondxparm(fp,filename,atm1[i],atm2[i],atm3[i],box_type[i],
		    &eq1_now,&eq2_now,&fk_now,&ifix);
      bondxs->eq1bondx[i] = eq1_now; 
      bondxs->eq2bondx[i] = eq2_now; 
      bondxs->fkbondx[i] = fk_now;
      bondxs->ifix[i] = ifix;
    }
    fclose(fp);
  }
  simparms->mem_bytes += 3*bondxs->nbondxs*sizeof(double);
  if(simparms->rank==0){
    sprintf(line,"number of bondx types = %d number of total cross bonds = %d"
	    ,nb_types,bondxs->nbondxs);
    md_stdout(line);
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }

#ifdef DEBUG
  for(i=0;i<bondxs->nbondxs;i++){
    j = bondxs->idxbox[i];
    sprintf(line,"%d %d %s %s %s %s %g %g %g",i,j,
	    atm1[j],atm2[j],atm3[j],box_type[j],
	    bondxs->eq1bondx[j],bondxs->eq2bondx[j],bondxs->fkbondx[j]);
    md_stdout(line);
  }
#endif
  free(atm1); free(atm2); free(atm3); free(box_type);
}
/*------------------------------------------------------------------------*/
