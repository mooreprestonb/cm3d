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

/* subroutine to read in free format simuluation input file
   and assign values or use default vaules */

#include "md.h"

#define DEBUG
#define TCONV       (1./1000.) /* fs -> ps */

/*---------------------------------------------------------*/
void read_sim_inputs(SIMPARMS *simparms,SUBSPACE *subspace,
		     WRITE_STEP *write_step,FILENAMES *filenames,
		     INTER *inter,NGBR *ngbr,COORDS *coords,STARTUP *startup)
{
  char line[MAXLINELEN];
  int num_dict,*num_found;
  WORD *dict;
  int i,num_keys=0,icase,iseed;
  KEY *key_root,*key_now;

#ifdef DMALLOC
  dmalloc_verify(NULL);
#endif

/* set defaults */
  set_sim_default(simparms,write_step,subspace,filenames,
		  inter,ngbr,coords,startup);

  /* set up key words */
  set_sim_keyword(&num_dict,&dict);
  num_found = (int *)cmalloc(num_dict*sizeof(int));
  for(i=0;i<num_dict;i++) num_found[i] = 0;
  
#ifdef DMALLOC
  dmalloc_verify(NULL);
#endif
  key_root = (KEY *)cmalloc(sizeof(KEY));
#ifdef DMALLOC
  dmalloc_verify(NULL);
#endif
  get_sim_keys(filenames->command,filenames->infile,&num_keys,key_root);

  for(i=0,key_now=key_root;i<num_keys;i++,key_now=key_now->next){
    icase = get_dict_num(num_dict,dict,key_now->keyword);
    if(icase == -1){
      print_dict(key_now->keyword,num_dict,dict);
    }
    num_found[icase]++;

    switch (icase) {
    case 0: /* calculation */
      if(strcasecmp("md",key_now->keyarg) == 0) {
	simparms->icalc_type = 0; 
      } else if (strcasecmp("subspace",key_now->keyarg) == 0){
        simparms->icalc_type = 1;
      } else if (strcasecmp("minimize",key_now->keyarg) == 0){
        simparms->icalc_type = 2;
      } else if (strcasecmp("nma",key_now->keyarg) == 0){
        simparms->icalc_type = 3;
      } else if (strcasecmp("colvar",key_now->keyarg) == 0){
        simparms->icalc_type = 4;
      } else {
        strcpy(line,
           "Calculation type must be md, subspace, minimize, NMA or colvar not ");
        strcat(line,key_now->keyarg);
        md_error(line);
      }
      break; 
    case 1: /* ensemble */
      if(strcasecmp("nve",key_now->keyarg) == 0){
        simparms->iensemble = 0; 
	simparms->nchain = simparms->nbar = simparms->ivol = 0;
      } else if (strcasecmp("nvt",key_now->keyarg) == 0){
        simparms->iensemble = 1; 
	simparms->nbar = simparms->ivol = 0;
      } else if (strcasecmp("npt_i",key_now->keyarg) == 0){
        simparms->iensemble = 2; 
	simparms->ivol = 1;
      } else if (strcasecmp("npt_ixy",key_now->keyarg) == 0){
        simparms->iensemble = 2; 
	simparms->ivol = 4;
      } else if (strcasecmp("npt_l",key_now->keyarg) == 0){
        simparms->iensemble = 2; 
	simparms->ivol = 2;
      } else if (strcasecmp("npt_f",key_now->keyarg) == 0){
        simparms->iensemble = 2; 
	simparms->ivol = 3;
      } else if (strcasecmp("min",key_now->keyarg) == 0){
        simparms->iensemble = 3; 
	simparms->nchain = simparms->nbar = simparms->ivol = 0;
      } else {
	sprintf(line,"ensemble must be nve,nvt,npt_{i,ixy,l,f} or min != %s",
	       key_now->keyarg);
	md_error(line);
      }
      break;
    case 2: /* ntime */
      sscanf(key_now->keyarg,"%d",&(simparms->nstep)); break; 
    case 3: /*istart */
      sscanf(key_now->keyarg,"%d",&(simparms->istart));break;
    case 4:  /* iperd */
      sscanf(key_now->keyarg,"%d",&(simparms->iperd)); break;
    case 5: /* dt */
      sscanf(key_now->keyarg,"%lg",&(simparms->dt));  break;
    case 6: /* t_ext */
      sscanf(key_now->keyarg,"%lg",&(simparms->temp)); break; 
    case 7: /*p_ext */
      sscanf(key_now->keyarg,"%lg",&(simparms->pext)); break; 
    case 8: /*nlen */
      sscanf(key_now->keyarg,"%d",&(simparms->nlen)); break; 
    case 9: /* neighbor_list */
      if(strcasecmp("nolist",key_now->keyarg) == 0){
        ngbr->ilist = 0; 
      } else if (strcasecmp("verlist",key_now->keyarg) == 0){
        ngbr->ilist = 1; 
      } else if (strcasecmp("linklist",key_now->keyarg) == 0){
        ngbr->ilist = 2; 
      } else {
	sprintf(line,"neighbor_list must be nolist,verlist or linklist != %s ",
		key_now->keyarg);
	md_error(line);
      }
      break;

    case 10: /*iwr_screen */
      sscanf(key_now->keyarg,"%d",&(write_step->nscrn));break; 
    case 11: /*iwr_dump */
      sscanf(key_now->keyarg,"%d",&(write_step->ndump));break; 
    case 12: /*iwr_inst */
      sscanf(key_now->keyarg,"%d",&(write_step->ninst));break; 
    case 13: /*iwr_confp */
      sscanf(key_now->keyarg,"%d",&(write_step->nconf));break; 
    case 14: /*iwr_confv */
      sscanf(key_now->keyarg,"%d",&(write_step->nvel)); break; 
    case 93: /*iwr_colv */
      sscanf(key_now->keyarg,"%d",&(write_step->ncolv)); break; 
    case 15: /*resamvel */
      sscanf(key_now->keyarg,"%d",&(write_step->resamvel)); break; 
      
    case 16: /* get initial subspace vecs from config(0) or file(1) */
      sscanf(key_now->keyarg,"%d",&(subspace->igetvec));break;/* ivec_init */
    case 17: /* subspce update frequency in steps (-1 is def and no update) */
      sscanf(key_now->keyarg,"%d",&(write_step->nupdss)); break;

    case 18: /* whether to print out eigen vectors for system*/
      if(strcasecmp("on",key_now->keyarg) == 0) {
	write_step->psysvec = 1;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	write_step->psysvec = 0;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;

    case 19: /* number of points in a histigram */
      sscanf(key_now->keyarg,"%d",&(simparms->npoints)); break; 

    case 20:                                           /* sim_name */
      md_warning("Keyword sim_name does nothing");
      break;
    case 21: strcpy(filenames->restart,key_now->keyarg); break;/* dpo_name */
    case 22: strcpy(filenames->initfile,key_now->keyarg); break;/* dpi_name */
    case 23: strcpy(filenames->instham,key_now->keyarg); break;/* ins_name */
    case 24: strcpy(filenames->configs,key_now->keyarg); break;/* pos_file */
    case 25: strcpy(filenames->velfile,key_now->keyarg); break;/* vel_file */
    case 94: strcpy(filenames->colvfile,key_now->keyarg); break;/* colv_file */
    case 97: strcpy(filenames->hillfile,key_now->keyarg); break;/* hill_file */
    case 26: strcpy(filenames->insteng,key_now->keyarg); break;/* eng_file */
    case 27: strcpy(filenames->setfile,key_now->keyarg); break;/*mol_set_file*/
      
    case 28: strcpy(startup->inter_file,key_now->keyarg); break;/*inter_file */
    case 29: strcpy(startup->bond_file,key_now->keyarg);break;/*bond_file */
    case 30: strcpy(startup->bend_file,key_now->keyarg);break;/*bend_file */
    case 31: strcpy(startup->tors_file,key_now->keyarg);break;/*torsion_file */
    case 32: strcpy(startup->onfo_file,key_now->keyarg);break;/*onefour_file */

    case 33: /* file to write out the normal mode spectra */
      sscanf(key_now->keyarg,"%s",filenames->nmfile);
      break;
    case 34: /* file to read in initial subspace vectors  */
      sscanf(key_now->keyarg,"%s",subspace->vecfile);
      break;
    case 35: /* a history file of subspace vectors used */
      sscanf(key_now->keyarg,"%s",subspace->vecconf);
      break;
    case 36: /* system vector files */
      sscanf(key_now->keyarg,"%s",filenames->sysvec);
      break;
    case 37: /* system vector files */
      sscanf(key_now->keyarg,"%s",filenames->molvec);
      break;
    case 38: /* file to write out the normal mode spectra */
      sscanf(key_now->keyarg,"%s",filenames->specfile);
      break;
    case 39: strcpy(filenames->forcefile,key_now->keyarg); break;/*force_file*/
   
    case 40: /*nres_lrf */
      sscanf(key_now->keyarg,"%d",&(simparms->ninter));
      if(simparms->ninter<=0) md_error("nres_lrf must be > 0");
      break;
    case 41: /*nres_intra*/
      sscanf(key_now->keyarg,"%d",&(simparms->ninnra));
      if(simparms->ninnra<=0) md_error("nres_intra must be > 0");
      break;
    case 42: /*len_nhc */
      sscanf(key_now->keyarg,"%d",&(simparms->nchain));   
      if(simparms->nchain<0) md_error("nchain must be >= 0");
      break;
    case 43: /* nbar */
      sscanf(key_now->keyarg,"%d",&(simparms->nbar));  
      if(simparms->nbar<0) md_error("nbar must be >= 0");
      break;
    case 44: /* rheal_res */
      inter->rheal = atof(key_now->keyarg); 
      if(inter->rheal<0) md_error("rheal_res must be >= 0");
      break;
    case 45: /* skin */
      sscanf(key_now->keyarg,"%lg",&inter->skin_ter);      
      if(inter->skin_ter<0) md_error("Skin must be >= 0");
      break;
    case 46: /* tau_nhc */
      sscanf(key_now->keyarg,"%lg",&startup->tau_nhc);   break;
    case 47: sscanf(key_now->keyarg,"%lg",&startup->tau_vol);  /* tau_vol */
      break;
    case 48: /* alp_ewald */
      sscanf(key_now->keyarg,"%lg",&startup->alp_ewald);  break;
    case 49: sscanf(key_now->keyarg,"%d",&startup->kmax); break;/* kmax */
    case 50: /*iseed*/
      sscanf(key_now->keyarg,"%d",&iseed);
      srandme(iseed);break;
    case 51: sscanf(key_now->keyarg,"%d",&startup->ishift); break;/* shift */
    case 52:                               /*iunits */
      if(strcasecmp("kelvin",key_now->keyarg) == 0) {
	simparms->iunits = 0; 
      } else if (strcasecmp("atomic",key_now->keyarg) == 0){
        simparms->iunits = 1;
      } else if (strcasecmp("kcal",key_now->keyarg) == 0){
        simparms->iunits = 2;
      } else if (strcasecmp("wavenumber",key_now->keyarg) == 0){
        simparms->iunits = 3;
      } else {
	strcpy(line,"Units must be kelvin, atomic, Kcal or wavenumber not ");
	strcat(line,key_now->keyarg);
	md_error(line);
      }
      break; 
    case 53: /* ntable */
      sscanf(key_now->keyarg,"%d",&inter->ntable);      break;
    case 54: /*scale_onfo */
      sscanf(key_now->keyarg,"%lg",&startup->scale_onfo); break;
    case 55: /*scale_onfo_e */
      sscanf(key_now->keyarg,"%lg",&startup->scale_onfo_e);break;
    case 56: /*rcute_max */
      sscanf(key_now->keyarg,"%lg",&startup->rcute_max);  break;
    case 57: /*rcute_min*/
      sscanf(key_now->keyarg,"%lg",&startup->rcute_min);  break;
    case 58: /*rcute_resp */
      sscanf(key_now->keyarg,"%lg",&startup->rcute_resp); break;
    case 59: sscanf(key_now->keyarg,"%d",&(subspace->nstate_max));break;
    case 60: strcpy(startup->nvecreal,key_now->keyarg);  break;
    case 61: strcpy(startup->nvecimag,key_now->keyarg);  break;
    case 62: sscanf(key_now->keyarg,"%lg",&(subspace->freq_min)); break;
    case 63: sscanf(key_now->keyarg,"%lg",&(subspace->freq_max)); break;
    case 64: sscanf(key_now->keyarg,"%d",&(subspace->num_subst)); break;
    case 65: strcpy(startup->nvecstore,key_now->keyarg); 
      break;/* vecstore temp array for numsubvec */
      
    case 66: /*resamvel */
      sscanf(key_now->keyarg,"%d",&(write_step->rescalevel)); break; 
    case 67: /* scale_subspace */
      if(strcasecmp("on",key_now->keyarg) == 0) {
	subspace->stemp_update = 1;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	subspace->stemp_update = 0;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;
    case 68: /*ncell_div */
      ngbr->ncells = atoi(key_now->keyarg);
      if(ngbr->ncells%2!=1){
	sprintf(line,"Cell divisions must be odd not %d\n",ngbr->ncells);
	md_error(line);
      }
      break;
    case 69: /*long_rc */
      if(strcasecmp("on",key_now->keyarg) == 0) {
	simparms->vlrc = 0;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	simparms->vlrc = -1;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;
    case 70: /*iwr_conff */
      sscanf(key_now->keyarg,"%d",&(write_step->nforce));break; 
    case 71: /* min_type */
      if(strcasecmp("powell",key_now->keyarg) == 0){
        simparms->imin_type = 0; 
      } else if (strcasecmp("cg",key_now->keyarg) == 0){
        simparms->imin_type = 1; 
      } else if (strcasecmp("simplex",key_now->keyarg) == 0){
        simparms->imin_type = 2; 
      } else if (strcasecmp("ir",key_now->keyarg) == 0){
        simparms->imin_type = 3; 
      } else {
	sprintf(line,"%s must be powell,cg,simplex or ir != %s",
		key_now->keyword,key_now->keyarg);
	md_error(line);
      }
      break;
    case 72: /* use cell */
      if(strcasecmp("on",key_now->keyarg) == 0) {
	simparms->nocell_all = 1;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	simparms->nocell_all = 0;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;
    case 73: /* max_exclude */
      sscanf(key_now->keyarg,"%d",&(simparms->max_exclude)); break; 
    case 74: /* iextern */
      if(strcasecmp("on",key_now->keyarg) == 0) {
	simparms->iextern = 1;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	simparms->iextern = 0;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;
    case 75: strcpy(filenames->extfile,key_now->keyarg);break;/*extern_file */
    case 76: /* min_tolerance */
      sscanf(key_now->keyarg,"%lg",&(simparms->min_tol)); break;
    case 77: /* wall_clock */
      sscanf(key_now->keyarg,"%lg",&(simparms->wallclock)); break; 
    case 78: strcpy(filenames->instext,key_now->keyarg); break;/* eval_file */
    case 79: /*nres_tors*/
      sscanf(key_now->keyarg,"%d",&(simparms->ninrra));
      if(simparms->ninrra<=0) md_error("nres_tors must be > 0");
      break;
    case 80: /*zmin_extern*/
      sscanf(key_now->keyarg,"%lg",&(coords->zmin_extern));
      break;
    case 81: /*zmax_extern*/
      sscanf(key_now->keyarg,"%lg",&(coords->zmax_extern));
      break;
    case 82:/* nres_nhc */ 
      sscanf(key_now->keyarg,"%d",&(simparms->nhc));
      if(simparms->nhc<=0) md_error("nres_nhc must be > 0");
      break;
    case 83:/* yoshidas */ 
      if(strcasecmp("on",key_now->keyarg) == 0) {
	simparms->nyoshida = 3;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	simparms->nyoshida = 1;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;
    case 84: /* scaleeps */
      sscanf(key_now->keyarg,"%lg",&(simparms->scaleeps)); break; 
    case 85: /* scalecharge */
      sscanf(key_now->keyarg,"%lg",&(simparms->scalecharge)); break;
    case 86: /* iwr_eigen_values */
      if(strcasecmp("on",key_now->keyarg) == 0) {
	write_step->peigval = 1;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	write_step->peigval = 0;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;
    case 87: /* file to write out participation ration */
      sscanf(key_now->keyarg,"%s",filenames->partrat);
      break;
    case 88: /* specify what type of NMA calculation */
      if(strcasecmp("ir",key_now->keyarg) == 0) {
        simparms->nma_type = 0; 
      } else if (strcasecmp("raman",key_now->keyarg) == 0){
        simparms->nma_type = 1;
      } else if (strcasecmp("all",key_now->keyarg) == 0){
        simparms->nma_type = 2;
      } else {
        strcpy(line,
	       "NMA calculation must be IR, Raman, or both, not ");
        strcat(line,key_now->keyarg);
        md_error(line);
      }
      printf("Running NMA calculation: nma_type=%d\n",simparms->nma_type);
      break;    
    case 89:/* polar forces */
      if(strcasecmp("on",key_now->keyarg) == 0) {
	simparms->ipolar = 1;
      } else if (strcasecmp("off",key_now->keyarg) == 0) {
	simparms->ipolar = 0;
      } else {
	sprintf(line,"Allowed arguments for %s are \"on\" or \"off\"",
		key_now->keyword);
	md_error(line);
      }
      break;
    case 90:/* ncolvar */
      sscanf(key_now->keyarg,"%d",&(coords->colvar.ncolvar));
      break;
    case 91:/* hilldepth */
      sscanf(key_now->keyarg,"%lg",&(coords->colvar.hilldepth));
      break;
    case 92:/* colvarcntrfile */
      strcpy(filenames->colvarcntrfile,key_now->keyarg); 
      break;
    case 95:
      sscanf(key_now->keyarg,"%d",&(coords->colvar.minstephill)); /* minstephill */
      break;
    case 96:
      sscanf(key_now->keyarg,"%d",&(coords->colvar.maxstephill)); /* maxstephill */
      break;
    case 98:
      sscanf(key_now->keyarg,"%d",&(coords->colvar.ltunehills)); /* tunehills */
      break;
    case 99: /* scaletemp */
      sscanf(key_now->keyarg,"%lg",&(simparms->scaletemp)); break; 
    default:
      print_dict(key_now->keyword,num_dict,dict);break;
    }
#ifdef DEBUG
    write_key(key_now);
#endif
  }

  /******** all variables must be read in by this point! ********/
  free(dict);

  /* convert from femtoseconds to picoseconds */
  simparms->dt *= TCONV;
  startup->tau_nhc *= TCONV;  
  startup->tau_vol *= TCONV; 

  simparms->pext /= PCONV;   /* convert form atm to K/A^3 */
  
  if(simparms->icalc_type==1){
    subspace->num_vecsub = (int *)cmalloc((subspace->num_subst)*sizeof(int));
    subspace->num_real = (int *)cmalloc((subspace->num_subst)*sizeof(int));
    subspace->num_imag = (int *)cmalloc((subspace->num_subst)*sizeof(int));
    subspace->num_vecsub_max=(int *)cmalloc((subspace->num_subst)*sizeof(int));
    subspace->num_real_max = (int *)cmalloc((subspace->num_subst)*sizeof(int));
    subspace->num_imag_max = (int *)cmalloc((subspace->num_subst)*sizeof(int));
    
    /* must be seperate for getnxtint to work correctly */
    for (i=0;i<subspace->num_subst;i++){
      subspace->num_vecsub_max[i]=getnxtint(startup->nvecstore);
    }
    for (i=0;i<subspace->num_subst;i++){
      subspace->num_real_max[i]=getnxtint(startup->nvecreal);
    }
    for (i=0;i<subspace->num_subst;i++){
      subspace->num_imag_max[i]=getnxtint(startup->nvecimag);
    }
  }
  /* ---------------------------- free KEYS memory --------------- */
  while(num_keys-->0){
    key_now = key_root->next; free(key_root);key_root = key_now;
  }
  free(key_root);

  if(num_found[50]==0) srandme(0);

  /* consistancy checks */

  /* MD calculation */
  if(simparms->icalc_type==0){
    switch(simparms->iensemble){
    case 0:   /* NVE */
      if(simparms->nchain != 0 ){
	md_error("You can't have a chain length with an NVE ensemble");
      }
      if(simparms->nbar != 0){
	md_error("You can't have a vol. barostat length with an NVE ensemble");
      }
      break;
    case 1:   /* NVT */
      if(simparms->nchain <= 0){
	md_error("You must have a chain length > 0  in the NVT ensemble");
      }
      if(simparms->nbar != 0){
	md_error("You can't have a vol. barostat length with an NVT ensemble");
      }
      break;
    case 2:   /* NPT */
      if(simparms->nchain <= 0){
	md_error("You must have a chain length > 0 in the NPT ensemble");
      }
      if(simparms->nbar < 0){
	md_error("You must have a volume baro length with an NPT ensemble");
      }
      if(simparms->nbar == 0){
	md_warning("No barostat - using only volume coordinate");
      }
      if(simparms->iperd != 3){
	md_error("You can't barostat without periodic boundies");
      }
      if(!simparms->ivol){
	if(simparms->nbar > 1){
	  sprintf(line,"You can't barostat a constant volume!\n\t");
	  sprintf(line,"%s Change either ensemble or nbar in file \"%s\"",
		  line,filenames->initfile);
	  md_error(line);
	}
      }
      break;
    case 3:   /* MIN */
      break;
    default:
      md_error("in the ensemble value?");
      break;
    }
    if( simparms->ninter <=0 || simparms->ninnra <=0 ||	simparms->ninrra <= 0){
      md_error("You must have a positive number of respa steps");
    }
    if( coords->colvar.ncolvar !=0 ) {
      md_error("You can't set ncolvar with an md simulation");
    }
  }
  /* subspace dynamics */
  if(simparms->icalc_type==1){
    if(simparms->iensemble != 0){
      md_error("Only nve ensemble is implimented in subspace dynamics :-(");
    }
    if( simparms->ninter > 1 || simparms->ninnra > 1  || simparms->ninrra > 1){
      md_error("Respa steps not implemented in subspace dynamics");
    }
    if( simparms->ninter <=0 || simparms->ninnra <=0 || simparms->ninrra <=0 ){
      md_error("You must have a positive number of respa steps");
    }
  }

  /* Minimize calculation */
  if(simparms->icalc_type==2){
    /* use nolist */
    simparms->nocell_all = 0;    
    ngbr->ilist = 0; 
  }
  /* normal mode calculation */
  if(simparms->icalc_type==3){
    simparms->nocell_all = 0;    
    ngbr->ilist = 0; 
    simparms->nchain = 0.; 
    simparms->nbar = 1.;
    if(simparms->npoints==0){
      md_error("You must have some points for the plots");
    }
  }

  if(simparms->icalc_type==4){ /* colvar calculation */
    if(coords->colvar.ncolvar < 1){
      md_error("You must have at least 1 colvar");
    }
    if(coords->colvar.hilldepth < 0){
      md_error("Hill depth must be >= 0");
    }
  }

  if(subspace->freq_max < subspace->freq_min){
    fprintf(stderr,"%g %g\n",subspace->freq_max,subspace->freq_min);
    md_error("You can't have the max frequency larger then the min frequency");
  }
  /* neighbor list check */
  if(ngbr->ilist==2 && simparms->nocell_all==0){
    md_warning("Link list used with nocell! -- changing nocell option");
    simparms->nocell_all = 1;
  }
  if(ngbr->ilist==0 && simparms->nocell_all==1){
    md_warning("Nolist used with link cell! (this is a waste of memory)");
    md_warning("\t -- changing nocell option");
    simparms->nocell_all = 0;
  }
  free(num_found);
}
/*---------------------------------------------------------------------*/ 
