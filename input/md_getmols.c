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

/* routines to read in the different arguments for the metakeys in a
   parameter file */

#include "md.h"

#define NUM_METAKEYS 8
WORD dictmetakey[NUM_METAKEYS] = {
  "mol_name_def","atom_def","bond_def","bend_def",
  "torsion_def","onefour_def","bondx_def","path_def"
};

/*---------------------------------------------------------------*/
void get_species(FILE *fp,SPECIES *spec_now)
{
  char line[MAXLINELEN];
  int i,ich,nline,nch,num_keys,icase;
  int num_dict=8,num_found[8],ifound;
  int napm,nbond,nbend,ntors,n14,nbondx,npath;
  WORD dict[8],metakey;
  KEY key_root,*key_now;

  strcpy(dict[0],"mol_name");  strcpy(dict[1],"natom");
  strcpy(dict[2],"nbond");     strcpy(dict[3],"nbend");
  strcpy(dict[4],"ntors");     strcpy(dict[5],"n14");
  strcpy(dict[6],"nbondx");    strcpy(dict[7],"npath");

  /* loop until end of file */
  rewind(fp);
  napm=nbond=nbend=ntors=n14=nbondx=npath=0;
  nline = nch = 0;
  while((ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      ifound = 0;
      get_meta_key(fp,spec_now->filename,
		   line,&nch,&nline,metakey,&key_root,&num_keys);
      
      if(!strcasecmp(dictmetakey[0],metakey)){
        ifound = 1;
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1){
            print_dict(key_now->keyword,num_dict,dict); exit(1);
          }
          num_found[icase]++;
          switch (icase) {
             case 0:/*mol_name_def*/
               strcpy(spec_now->name,key_now->keyarg);break;
             case 1:/*nmol*/
               sscanf(key_now->keyarg,"%d",&(spec_now->napm));break;
             case 2:/*bond*/
               sscanf(key_now->keyarg,"%d",&(spec_now->nbond));break;
             case 3:/*bend*/
               sscanf(key_now->keyarg,"%d",&(spec_now->nbend));break;
             case 4:/*tors*/
               sscanf(key_now->keyarg,"%d",&(spec_now->ntors));break;
             case 5:/*n14*/
               sscanf(key_now->keyarg,"%d",&(spec_now->n14));break; 
             case 6:/*bondx*/
               sscanf(key_now->keyarg,"%d",&(spec_now->nbondx));break;
             case 7:/*npath - do nothing for now....*/
               md_warning("Key word npath being ignored");
               break; 
             default:
               print_dict(key_now->keyword,num_dict,dict);  exit(1);  break;
          }
        }
      }
      if(!strcasecmp(dictmetakey[1],metakey)) {napm++;  ifound=1;}
      if(!strcasecmp(dictmetakey[2],metakey)) {nbond++; ifound=1;}
      if(!strcasecmp(dictmetakey[3],metakey)) {nbend++; ifound=1;}
      if(!strcasecmp(dictmetakey[4],metakey)) {ntors++; ifound=1;}
      if(!strcasecmp(dictmetakey[5],metakey)) {n14++;   ifound=1;}
      if(!strcasecmp(dictmetakey[6],metakey)) {nbondx++;ifound=1;}
      if(!strcasecmp(dictmetakey[7],metakey)) {npath++; ifound=1;}

      if(!ifound) print_dict(metakey,NUM_METAKEYS,dictmetakey);
      
      /* free link list of keys */
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch = 0;nline++;}
  }

  if(spec_now->napm != napm){
    sprintf(line,"# of atoms specified %d != %d # of atoms found\n\tsetting # of atoms to %d",spec_now->napm,napm,napm);
    md_warning(line);
    spec_now->napm = napm;
  }

  if(spec_now->nbond != nbond){
    sprintf(line,"# of bonds specified %d != %d # of bonds found\n\tsetting # of bond to %d",spec_now->nbond,nbond,nbond);
    md_warning(line);
    spec_now->nbond = nbond;
  }

  if(spec_now->nbend != nbend){
    sprintf(line,"# of bends specified %d != %d # of bends found\n\tsetting # of bend to %d",spec_now->nbend,nbend,nbend);
    md_warning(line);
    spec_now->nbend = nbend;
  }
  
  if(spec_now->ntors != ntors){
    sprintf(line,"# of torsions specified %d != %d # of torsons found\n\tsetting # of torsions to %d",spec_now->ntors,ntors,ntors);
    md_warning(line);
    spec_now->ntors = ntors;
  }

  if(spec_now->n14 != n14){
    sprintf(line,"# of 1-4 specified %d != %d # of 1-4 found\n\tsetting # of 1-4 to %d",spec_now->n14,n14,n14);
    md_warning(line);
    spec_now->n14 = n14;
  }

  if(spec_now->nbondx != nbondx){
    sprintf(line,"# of cross bonds specified %d != %d # of cross bonds found\n\tsetting # of cross bonds to %d",spec_now->nbondx,nbondx,nbondx);
    md_warning(line);
    spec_now->nbondx = nbondx;
  }

  rewind(fp);
}
/*-------------------------------------------------------------------------*/

void get_atom(FILE *fp,char *filename,int napm,ATOM_TOPOL *atom_topol)
{
  int i,ich,nline,nch,num_keys,icase;
  int natoms;
  char line[MAXLINELEN];
  int  num_dict = 8,num_found[8];
  WORD dict[8];
  WORD metakey;
  KEY key_root,*key_now;

  strcpy(dict[0],"atom_typ");  strcpy(dict[1],"atom_ind");
  strcpy(dict[2],"mass");      strcpy(dict[3],"charge");
  strcpy(dict[4],"group");     strcpy(dict[5],"valence");
  strcpy(dict[6],"alpha");     strcpy(dict[7],"res");


  /* loop until end of file */
  rewind(fp);
  nline = nch = natoms = 0;
  while((ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      
      if(!strncasecmp(dictmetakey[1],metakey,strlen(dictmetakey[1]))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1){
            print_dict(key_now->keyword,num_dict,dict);
            exit(1);
          }
          num_found[icase]++;
          switch (icase) {
             case 0: /* atom_typ */
               strcpy(atom_topol[natoms].type,key_now->keyarg);
               if(natoms+1 > napm) {
                 sprintf(line,"Number of atoms defined exceeds %d",napm);
                 md_error(line);
               }
               break;
               
             case 1:  /* atom index */
               atom_topol[natoms].atm_idx = (int) atof(key_now->keyarg);
               break;
             case 2: /* mass */
               atom_topol[natoms].mass = atof(key_now->keyarg);
               break;
             case 3: /* charge */
               atom_topol[natoms].charge = atof(key_now->keyarg);
               break;
             case 4: /* group name */
               strcpy(atom_topol[natoms].group,key_now->keyarg);
               break;
             case 5: /* valence */
               break;
             case 6: /* alpha polarlizablitity */
               atom_topol[natoms].alpha = atof(key_now->keyarg);
               break;
             case 7: /* residue */
               break;
             default:
               break; /* do nothing now */
               print_dict(key_now->keyword,num_dict,dict);
               exit(1); break;
          }
        }
        /* make sure all keys were found */
        for(i=0;i<3;i++){
          if(num_found[i] != 1){
            err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(atom_topol[natoms].mass <= 0.0){
          sprintf(line,"You have a mass of %g on atom %s (%d) in file %s",
             atom_topol[natoms].mass,atom_topol[natoms].type,
             atom_topol[natoms].atm_idx,filename);
          md_error(line);
        }
        if(num_found[3] == 0) atom_topol[natoms].charge = 0.;
        if(num_found[4] == 0) atom_topol[natoms].group[0] = '\0';
        if(num_found[6] == 0) atom_topol[natoms].alpha = 0.;
        
        natoms++;
        if(natoms>napm){
          sprintf(line,"number of atoms found %d > %d specified",natoms,napm);
          md_error(line);
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch = 0;nline++;}
  }
  rewind(fp);
  if(natoms!=napm){
    sprintf(line,"number of atoms found %d != %d specified\n",natoms,napm);
    md_error(line);
  }
}
/*-------------------------------------------------------------------------*/
void get_bond(FILE *fp,char *filename,int nbond,int *ibo,int *jbo,WORD *tbo)
{
  int i,ich,nline,nch,num_keys,icase;
  int nbonds;
  char line[MAXLINELEN];
  int  num_dict = 4,num_found[4];
  WORD dict[4],metakey;
  KEY key_root,*key_now;

  strcpy(dict[0],"atom1");  strcpy(dict[1],"atom2");
  strcpy(dict[2],"label");  strcpy(dict[3],"cons");

  /* loop until end of file */
  rewind(fp);
  nline = nch = nbonds = 0;
  while( (ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      if(!strncasecmp(dictmetakey[2],metakey,strlen(dictmetakey[2]))){
	for(i=0;i<num_dict;i++) num_found[i] = 0;
	for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
	  icase = get_dict_num(num_dict,dict,key_now->keyword);
	  if(icase == -1){
	    print_dict(key_now->keyword,num_dict,dict);
	    exit(1);
	  }
	  num_found[icase]++;
	  switch (icase) {
        case 0:sscanf(key_now->keyarg,"%d",&ibo[nbonds]); break; /*atom1*/
        case 1:sscanf(key_now->keyarg,"%d",&jbo[nbonds]); break; /*atom2*/
        case 2:strcpy(tbo[nbonds],key_now->keyarg);  break;
        case 3:
          sprintf(line,"constraints not implimented yet");
          md_error(line); break;
        default:
          print_dict(key_now->keyword,num_dict,dict);  exit(1); break;
	  }
	}
	for(i=0;i<2;i++){
	  if(num_found[i] != 1){
	    err_found(dict[i],num_found[i],filename,nline); 
	  }
	}
	if(num_found[2]==0) strcpy(tbo[nbonds],"");
	nbonds++;
	if(nbonds>nbond){
	  sprintf(line,"number of bonds found %d > %d specified",nbonds,nbond);
	  md_error(line);
	}
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch = 0;nline++;}
  }
  rewind(fp);
  if(nbonds!=nbond){
    sprintf(line,"number of bonds found %d != %d specified",nbonds,nbond);
    md_error(line);
  }
}
/*-------------------------------------------------------------------------*/
void get_bend(FILE *fp,char *filename,int nbend,int *ibe,int *jbe,int *kbe,
	      int *tybe,WORD *tbe)
{
  int i,ich,nline,nch,num_keys,icase,nbends;
  char line[MAXLINELEN];
  int  num_dict=5,num_found[5];
  WORD dict[5],metakey;
  KEY key_root,*key_now;

  strcpy(dict[0],"atom1");  strcpy(dict[1],"atom2");
  strcpy(dict[2],"atom3");  strcpy(dict[3],"label");
  strcpy(dict[4],"bendtype");
  /* loop until end of file */
  rewind(fp);
  nline = nch = nbends = 0;
  while((ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      if(!strncasecmp(dictmetakey[3],metakey,strlen(dictmetakey[3]))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1){
            print_dict(key_now->keyword,num_dict,dict);
            exit(1);
          }
          num_found[icase]++;
          switch (icase) {
	  case 0:sscanf(key_now->keyarg,"%d",&ibe[nbends]);break; /*atom1*/
	  case 1:sscanf(key_now->keyarg,"%d",&jbe[nbends]);break; /*atom2*/
	  case 2:sscanf(key_now->keyarg,"%d",&kbe[nbends]);break; /*atom3*/
	  case 3:strcpy(tbe[nbends],key_now->keyarg); break;
	  case 4:sscanf(key_now->keyarg,"%d",&tybe[nbends]);break; /*bendtype*/
	  default:
	    print_dict(key_now->keyword,num_dict,dict);exit(1); break;
          }
        }
        for(i=0;i<3;i++){
          if(num_found[i] != 1){
            err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(num_found[3]==0) strcpy(tbe[nbends],"");
        if(num_found[4]==0) tybe[nbends] =0;
        nbends++;
        if(nbends>nbend){
          sprintf(line,"number of bends found %d > %d specified",nbends,nbend);
          md_error(line);
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch=0;nline++;}
  }
  rewind(fp);
  if(nbends!=nbend){
    sprintf(line,"number of bends found %d != %d specified",nbends,nbend);
    md_error(line);
  }
}
/*-------------------------------------------------------------------------*/
void get_tors(FILE *fp,char *filename,int ntors,int *ito,int *jto,
   int *kto,int *lto,WORD *tto)
{
  int i,ich,nline,nch,num_keys,icase,ntorsion;
  char line[MAXLINELEN];
  int  num_dict = 5,num_found[5];
  WORD dict[5],metakey;
  KEY key_root,*key_now;

  strcpy(dict[0],"atom1");    strcpy(dict[1],"atom2");
  strcpy(dict[2],"atom3");    strcpy(dict[3],"atom4");
  strcpy(dict[4],"label");

  /* loop until end of file */
  rewind(fp);
  nline = nch = ntorsion = 0;
  while( (ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      
      if(!strncasecmp(dictmetakey[4],metakey,strlen(dictmetakey[4]))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1){
            print_dict(key_now->keyword,num_dict,dict);
            exit(1);
          }
          num_found[icase]++;
          switch (icase) {
             case 0:sscanf(key_now->keyarg,"%d",&ito[ntorsion]); break; /*atom1*/
             case 1:sscanf(key_now->keyarg,"%d",&jto[ntorsion]); break; /*atom2*/
             case 2:sscanf(key_now->keyarg,"%d",&kto[ntorsion]); break; /*atom3*/
             case 3:sscanf(key_now->keyarg,"%d",&lto[ntorsion]); break; /*atom4*/
             case 4:strcpy(tto[ntorsion],key_now->keyarg); break;  /*label*/
             default:
               print_dict(key_now->keyword,num_dict,dict);  exit(1); break;
          }
        }
        for(i=0;i<4;i++){
          if(num_found[i] != 1){
            err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(num_found[4] == 0) strcpy(tto[ntorsion],"");
        ntorsion++;
        if(ntorsion>ntors){
          sprintf(line,"number of torsions found %d > %d specified",
             ntorsion,ntors);
          md_error(line);
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch=0;nline++;}
  }
  rewind(fp);
  if(ntorsion!=ntors){
    sprintf(line,"number of torsions found %d != %d specified",ntorsion,ntors);
    md_error(line);
  }
}
/*-------------------------------------------------------------------------*/
void get_onefour(FILE *fp,char *filename,int nof,int *iof,int *jof)
{
  int i,ich,nline,nch,num_keys,icase,nofs;
  char line[MAXLINELEN];
  int  num_dict=2,num_found[2];
  WORD dict[2],metakey;
  KEY key_root,*key_now;
  
  strcpy(dict[0],"atom1");  strcpy(dict[1],"atom2");

  /* loop until end of file */
  rewind(fp);
  nline = nch = nofs = 0;
  while( (ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      
      if(!strncasecmp(dictmetakey[5],metakey,strlen(dictmetakey[5]))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1){
            print_dict(key_now->keyword,num_dict,dict);
            exit(1);
          }
          num_found[icase]++;
          switch (icase) {
             case 0:sscanf(key_now->keyarg,"%d",&iof[nofs]);break; /*atom1*/
             case 1:sscanf(key_now->keyarg,"%d",&jof[nofs]);break; /*atom2*/
             default:
               print_dict(key_now->keyword,num_dict,dict); exit(1); break;
          }
        }
        for(i=0;i<2;i++){
          if(num_found[i] != 1){
            err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        nofs++;
        if(nofs>nof){
          sprintf(line,"number of 1-4s found %d > %d specified",nofs,nof);
          md_error(line);
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch = 0;nline++;}
  }
  rewind(fp);
  if(nofs!=nof){
    sprintf(line,"number of 1-4s found %d != %d specified",nofs,nof);
    md_error(line);
  }
}
/*---------------------------------------------------------------------*/
void get_bondx(FILE *fp,char *filename,int nbondx,
   int *ibondx,int *jbondx,int *kbondx,WORD *tbondx)
{
  char line[MAXLINELEN];
  int i,ich,nline,nch,num_keys,icase,nbondxs;
  int  num_dict=4,num_found[4];
  WORD dict[4],metakey;
  KEY key_root,*key_now;

  strcpy(dict[0],"atom1");  strcpy(dict[1],"atom2");
  strcpy(dict[2],"atom3");  strcpy(dict[3],"label");
  /* loop until end of file */
  rewind(fp);
  nline = nch = nbondxs = 0;
  while((ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == META_KEY_CHAR){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      if(!strncasecmp(dictmetakey[6],metakey,strlen(dictmetakey[6]))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1){
            print_dict(key_now->keyword,num_dict,dict); exit(1);
          }
          num_found[icase]++;
          switch (icase) {
             case 0:sscanf(key_now->keyarg,"%d",&ibondx[nbondxs]);break;/*atom1*/
             case 1:sscanf(key_now->keyarg,"%d",&jbondx[nbondxs]);break;/*atom2*/
             case 2:sscanf(key_now->keyarg,"%d",&kbondx[nbondxs]);break;/*atom3*/
             case 3:strcpy(tbondx[nbondxs],key_now->keyarg); break;
             default:
               print_dict(key_now->keyword,num_dict,dict);exit(1); break;
          }
        }
        for(i=0;i<3;i++){
          if(num_found[i] != 1){
            err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(num_found[3]==0) strcpy(tbondx[nbondxs],"");
        nbondxs++;
        if(nbondxs>nbondx){
          sprintf(line,"number of cross bonds found %d>%d specified",
             nbondxs,nbondx);
          md_error(line);
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch=0;nline++;}
  }
  rewind(fp);
  if(nbondxs!=nbondx){
    sprintf(line,"numbor of cross bonds found %d != %d specified",
	    nbondxs,nbondx);
    md_error(line);
  }
}
/*------------------------------------------------------------------------*/
