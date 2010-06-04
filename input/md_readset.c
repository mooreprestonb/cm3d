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

/* subroutines to read in potential parameters */

#include "md.h"
#define MIN_MEM 100
#define FREEZE_DEF "mol_freeze"
/* #define DEBUG */
/*--------------------------------------------------------------*/
void read_setfile(SIMPARMS *simparms,char *setfile,
		  int *nspec,SPECIES *spec_root,COORDS *coords)
{
  int i,ich,nline,nch,num_keys,icase,ispec;
  SPECIES *spec_now=spec_root;

  char line[MAXLINELEN];
  int  num_dict=5,num_found[5],mfz;
  WORD metakey,dict[5];
  KEY key_root,*key_now;  
  FILE *fp;

  strcpy(dict[0],"mol_parm_file"); strcpy(dict[1],"mol_therm_opt");
  strcpy(dict[2],"nmol");          strcpy(dict[3],"mol_index");
  strcpy(dict[4],"atom_freeze");

  spec_root->nmol  = spec_root->molidx = spec_root->napm  = 0;
  spec_root->nbond = spec_root->nbend  = spec_root->ntors = 0;
  spec_root->n14   = spec_root->nbondx = 0;
  spec_root->next  = NULL;

  simparms->nfreeze=0;
  mfz = MIN_MEM;
  coords->ifreeze = (int *)malloc(mfz*sizeof(int));

  ispec = 0;
 
  /* open setfile */
  if ((fp=fopen(setfile,"r"))==NULL){
    sprintf(line,"can't open set file \"%s\"",setfile);
    md_error(line);
  } else {
    if(simparms->rank==0){
      sprintf(line,"Reading in set from file \"%s\"",setfile);
      md_stdout(line);
    }
  }
  /* loop until end of file */
  nline = nch = ispec = 0;
  while( (ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      get_meta_key(fp,setfile,line,&nch,&nline,metakey,&key_root,&num_keys);

#ifdef DEBUG
      md_stdout(metakey);
      for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
	sprintf(line,"%d %s %s",i,key_now->keyword,key_now->keyarg);
	md_stdout(line);
      }
#endif
      if(!strcasecmp(metakey,MOL_DEF)){
	for(i=0;i<num_dict;i++) num_found[i] = 0;
	for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
	  icase = get_dict_num(num_dict,dict,key_now->keyword);
	  if(icase == -1){
	    print_dict(key_now->keyword,num_dict,dict);
	  }
	  num_found[icase]++;
	  switch (icase) {
	  case 0:strcpy(spec_now->filename,key_now->keyarg); /*mol_parm_file*/
	    break;
	  case 1:strcpy(spec_now->thermopt,key_now->keyarg); /*mol_therm_opt*/
	    break;
	  case 2:sscanf(key_now->keyarg,"%d",&(spec_now->nmol)); /*nmol*/
	    break;
	  case 3:sscanf(key_now->keyarg,"%d",&(spec_now->molidx));/*mol_index*/
	    break;
	  default:
	    print_dict(key_now->keyword,num_dict,dict); break;
	  }
	}
	for(i=0;i<3;i++){
	  if(num_found[i] != 1){
	    err_found(dict[i],num_found[i],setfile,nline);
	  }
	}
	++ispec;
	spec_now->next = (SPECIES *)cmalloc(sizeof(SPECIES));
	spec_now = spec_now->next;
	spec_now->nmol  = spec_now->molidx = spec_now->napm = 0;
	spec_now->nbond = spec_now->nbend  = spec_now->ntors = 0;
	spec_now->n14   = spec_now->nbondx = 0;
	spec_now->next = NULL;
      }
      if(!strcasecmp(metakey,FREEZE_DEF)){
	for(i=0;i<num_dict;i++) num_found[i] = 0;
	for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
	  icase = get_dict_num(num_dict,dict,key_now->keyword);
	  if(icase == -1){
	    print_dict(key_now->keyword,num_dict,dict);
	  }
	  num_found[icase]++;
	  switch (icase) {
	  case 4:
	    if(mfz-2<simparms->nfreeze){
	      mfz += MIN_MEM;
	      coords->ifreeze=(int *)realloc(coords->ifreeze,mfz*sizeof(int));
	    }
	    coords->ifreeze[simparms->nfreeze] = (int)atof(key_now->keyarg);
	    simparms->nfreeze++;
	    break;
	  default:
	    print_dict(key_now->keyword,num_dict,dict); break;
	  }
	}
	for(i=4;i<5;i++){
	  if(num_found[i] != 1){
	    err_found(dict[i],num_found[i],setfile,nline);
	  }
	}
      }
      free_key(num_keys,&key_root);
    }
    if(ich == NEWLINE){nch = 0; nline++;}
  }
  fclose(fp);

  if(simparms->nfreeze==0){
    free(coords->ifreeze);
  } else {
    coords->ifreeze = (int*)realloc(coords->ifreeze,
				    simparms->nfreeze*sizeof(int));
    coords->ffreeze = (double *)malloc(simparms->nfreeze*sizeof(double)*DIM);
    for(i=0;i<DIM*(simparms->nfreeze);++i) coords->ffreeze[i] = 0.;
  }
#ifdef DEBUG
  sprintf(line,"species = %d",ispec);
  md_stdout(line);
  for(i=0,spec_now=spec_root;i<ispec;i++,spec_now = spec_now->next){
    sprintf(line,"%s %s %d %d",spec_now->filename,spec_now->thermopt,
	   spec_now->nmol,spec_now->molidx);
    md_stdout(line);
  }

  printf("frozen atoms = %d\n",simparms->nfreeze);
  for(i=0;i<simparms->nfreeze;i++){
    printf("%d %d\n",i,coords->ifreeze[i]);
  }
#endif
  *nspec = ispec;
  if(*nspec <= 0 ) {
    md_error("you must specify at least one species in the setfile");
  }
}

/*--------------------------------------------------------------*/
