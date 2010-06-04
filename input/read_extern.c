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

#define POT_TYPE_HAUTMAN "hautman"
#define POT_TYPE_TABLE   "table"
#define META_KEY_EXTN    "external"
#define META_KEY_EXTNE   "elec_external"

/* #define DEBUG */
/* #define WRITE_TABLE */

void read_extern_pot(char *filename,int ntypes,WORD types[],
		     int ntable,COORDS *coords)
{
  int i,indx,ifound,ich,nline,nch,num_keys,icase;
  int num_dict=6,num_found[6];
  double c12,c3,xt,f,df,d2f;
  double zmin,zmax,dx,z0,conv;
  double **dvtab,**vtab,*dvtab_e,*vtab_e,*dv2tab;
  WORD dict[6],metakey,atom1,pot_type,tfile;
  KEY key_root,*key_now;
  LINE line;
  FILE *fp;

  /* open data base file */
  if((fp = fopen(filename,"r")) == NULL){
    sprintf(line,"can't open external potential file (%s)",filename);
    md_error(line);
  } else {
    sprintf(line,"Reading %s for the external potential parameters",filename);
    sprintf(line,"%s\nElectrostatic Potential must be in Volts!!",line);
    md_stdout(line);
  }

  zmin = coords->zmin_extern;
  zmax = coords->zmax_extern;
  dx = (zmax- zmin)/(double)(ntable-1);
  coords->dxtab_extern= dx;
    
  strcpy(dict[0],"atom");         strcpy(dict[1],"pot_type");     
  strcpy(dict[2],"z0");           strcpy(dict[3],"c12");          
  strcpy(dict[4],"c3");           strcpy(dict[5],"tab_file");

  vtab_e  = coords->vtab_extern_e;
  dvtab_e = coords->dvtab_extern_e;
  vtab    = coords->vtab_extern;
  dvtab   = coords->dvtab_extern;

  dv2tab = cmalloc(ntable*sizeof(double));
  /* loop until end of file */
  for(indx=0;indx<ntypes;indx++){
    rewind(fp);
    nline = nch = ifound = 0;
    while( (ich = fgetc(fp)) != EOF && !(ifound)) {
      line[nch++] = (char)ich;
      if(ich == '~'){
	get_meta_key(fp,filename,
		     line,&nch,&nline,metakey,&key_root,&num_keys);
	if(!strncasecmp(META_KEY_EXTN,metakey,strlen(META_KEY_INTER))){
	  for(i=0;i<num_dict;i++) num_found[i] = 0;
	  for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
	    icase = get_dict_num(num_dict,dict,key_now->keyword);
	    if(icase == -1){
	      print_dict(key_now->keyword,num_dict,dict); exit(1);
	    }
	    num_found[icase]++;
	    switch (icase) {
	    case 0:strcpy(atom1,key_now->keyarg); break;	/* atom1 */
	    case 1:strcpy(pot_type,key_now->keyarg); break;/* pot_type */
	    case 2:sscanf(key_now->keyarg,"%lg",&z0);break;/*min_dist */
	    case 3:sscanf(key_now->keyarg,"%lg",&c12);      break;/* c12 */
	    case 4:sscanf(key_now->keyarg,"%lg",&c3);       break;/* c3 */
	    case 5:strcpy(tfile,key_now->keyarg); break;   /* tab_file */
	    default:
	      print_dict(key_now->keyword,num_dict,dict); exit(1); break;
	    }
	  }
	  /* found and got keys see if match the ones we are looking for*/
	  if(!strcasecmp(atom1,types[indx])){
	    if(!(!strcasecmp(POT_TYPE_NULL,pot_type) ||
		 !strcasecmp(POT_TYPE_HAUTMAN,pot_type) ||
		 !strcasecmp(POT_TYPE_TABLE,pot_type))){
	      fprintf(stderr,"ERROR: pot_type %s not implimented\n",pot_type);
	      fprintf(stderr,"in file %s at line %d\n",filename,nline);
	      exit(1);
	    }
#ifdef DEBUG
	    fprintf(stdout,"pot_type for atom type %s is %s\n",atom1,pot_type);
#endif
	    /* check for required key words */
	    for(i=0;i<2;i++){
	      if(num_found[i] != 1)
		err_found(dict[i],num_found[i],filename,nline); 
	    }
	    if(!strcasecmp(POT_TYPE_HAUTMAN,pot_type)){
	      for(i=2;i<5;i++){
		if(num_found[i] != 1) 
		  err_found(dict[i],num_found[i],filename,nline); 
	      }
	    }
	    ifound = 1;
	    if(!strcasecmp(POT_TYPE_TABLE,pot_type)){
	      if(num_found[5] != 1) 
		err_found(dict[i],num_found[i],filename,nline); 
	      read_table(tfile,ntable,zmin,zmax,vtab[indx],dvtab[indx],
			 dv2tab);
	    } else {
	      /* tabulate the energy and forces */
	      for(i=0;i<ntable;i++){
		xt = i*dx + zmin;
		if(!strcasecmp(POT_TYPE_HAUTMAN,pot_type))
		  func_hautman(xt,&f,&df,&d2f,z0,c12,c3);
		if(!strcasecmp(POT_TYPE_NULL,pot_type)) f = df = d2f = 0.;
		vtab[indx][i] = f; 
		dvtab[indx][i] = df; 
		/* d2vtab[indx][i] = d2f; */
	      }
	    }
	  }
	}
	free_key(num_keys,&key_root);
      }
      if(ich == '\n'){nch = 0;nline++;}
    }
    if(!ifound){
      fprintf(stderr,
	      "ERROR: External interaction for type %s not found in \"%s\"\n",
	      types[indx],filename);
      exit(1);
    }
  }
#ifdef DEBUG
  md_stdout("Assigned external potentials");
#endif  

  conv = VCONV/sqrt(CCONV); /*the sqrt is to compensate for the charge unit */

  /* assign electrostatic potential */
  rewind(fp);
  nline = nch = ifound = 0;
  while( (ich = fgetc(fp)) != EOF && !(ifound)) {
    line[nch++] = (char)ich;
    if(ich == '~'){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      if(!strncasecmp(META_KEY_EXTNE,metakey,strlen(META_KEY_INTER))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1){
            print_dict(key_now->keyword,num_dict,dict); exit(1);
          }
          num_found[icase]++;
          switch (icase) {
             case 0:strcpy(atom1,key_now->keyarg); break;	/* atom1 */
             case 1:strcpy(pot_type,key_now->keyarg); break;/* pot_type */
             case 2:sscanf(key_now->keyarg,"%lg",&z0);break;/*min_dist */
             case 3:sscanf(key_now->keyarg,"%lg",&c12); break;/* c12 */
             case 4:sscanf(key_now->keyarg,"%lg",&c3);  break;/* c3 */
             case 5:strcpy(tfile,key_now->keyarg); break;   /* tab_file */
             default:
               print_dict(key_now->keyword,num_dict,dict); exit(1); break;
          }
        }
        /* found and got keys see if match the ones we are looking for*/
        if(!(!strcasecmp(POT_TYPE_NULL,pot_type) ||
           !strcasecmp(POT_TYPE_HAUTMAN,pot_type) ||
           !strcasecmp(POT_TYPE_TABLE,pot_type))){
          fprintf(stderr,"ERROR: pot_type %s not implimented\n",pot_type);
          fprintf(stderr,"in file %s at line %d\n",filename,nline);
          exit(1);
        }
#ifdef DEBUG
        fprintf(stdout,"Electro pot type is %s\n",pot_type);
#endif
        ifound = 1;
        /* check for required key words we don't need atoms def */
        for(i=1;i<2;i++){
          if(num_found[i]!=1)err_found(dict[i],num_found[i],filename,nline); 
        }
        if(!strcasecmp(POT_TYPE_HAUTMAN,pot_type)){
          for(i=2;i<4;i++){
            if(num_found[i] != 1)
              err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(!strcasecmp(POT_TYPE_TABLE,pot_type)){
          if(num_found[5] != 1) 
            err_found(dict[i],num_found[i],filename,nline); 
          read_table(tfile,ntable,zmin,zmax,vtab_e,dvtab_e,dv2tab);
        } else {
          /* tabulate the energy and forces */
          for(i=0;i<ntable;i++){
            xt = i*dx + zmin;
            if(!strcasecmp(POT_TYPE_HAUTMAN,pot_type))
              func_hautman(xt,&f,&df,&d2f,z0,c12,c3);
            if(!strcasecmp(POT_TYPE_NULL,pot_type)) f = df = d2f = 0.;
            vtab_e[i] = f; 
            dvtab_e[i] = df; 
            /* d2vtab_e[i] = d2f; */
          }
        }
        for(i=0;i<ntable;i++){
          vtab_e[i] *= conv; 
          dvtab_e[i] *= conv; 
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch = 0;nline++;}
  }
  if(!ifound){
    fprintf(stderr,"ERROR: External electric field was not found in \"%s\"\n",
	    filename);
    fprintf(stderr,"\t keyword ~%s should be in the file\n",META_KEY_EXTNE);
    exit(1);
  }
  free(dv2tab);

#ifdef WRITE_TABLE
  for(i=0;i<ntypes;i++){
    sprintf(tfile,"table%d_e.dat",i); 
    sprintf(line,"Saving external %d in file %s",i,tfile);
    md_stdout(line);
    fp = fopen(tfile,"w");
    for(indx=0;indx<ntable;indx++){
      xt = indx*coords->dxtab_extern + coords->zmin_extern;
      fprintf(fp,"%g %g %g\n",xt,coords->vtab_extern[i][indx],
	      coords->dvtab_extern[i][indx]);
    }
    fclose(fp);
  }
  sprintf(tfile,"tablee_e.dat"); 
  sprintf(line,"Saving external %d in file %s",i,tfile);
  md_stdout(line);
  fp = fopen(tfile,"w");
  for(indx=0;indx<ntable;indx++){
    xt = indx*coords->dxtab_extern + coords->zmin_extern;
    fprintf(fp,"%g %g %g\n",xt,coords->vtab_extern_e[indx],
	    coords->dvtab_extern_e[indx]);
  }
  fclose(fp);
#endif
}

