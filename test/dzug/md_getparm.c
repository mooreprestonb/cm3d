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
   subroutines that search a file for parameter key words 
   and sets up returns potential paramets
   */

#define RMAX_DIST (10.)
#define RMIN_DIST (1.)
#define RRES_DIST (7.)

#include "md.h"

/* #define DEBUG */

/*-----------------------------------------------------------------------*/
void get_ter_parm(int nmeta_key,META_KEY *meta_key,
   char *filename,WORD atm1,WORD atm2,
   double *rminl,double *rmaxl,double *rmins,double *rmaxs,
   double *vlrc,double *wlrc,INTER *inter,int indx,int *inull,
   double scaleeps)
{
  char line[MAXLINELEN];
  int imkey,i,ifound,nline,num_keys,icase;
  int num_dict=20,num_found[20];
  WORD dict[20];
  WORD atcm1,atcm2,atom1,atom2,pot_type,tword;
  META_KEY *meta_now;
  KEYS *key_now;
  double eps_now,sig_now,rmin_dist,rmax_dist,rres_dist;
  double awill,bwill,c6,c8,c9,c10,yamma,gamma,drmexp,c12,d10;
  double xt,rch,rhl,f,df,d2f,br,sm;
  char fname[MAXLINELEN];
  FILE *ftab;

  if(strcasecmp(atm1,atm2)<=0){ 
    strcpy(atcm1,atm1); strcpy(atcm2,atm2);
  } else {
    strcpy(atcm2,atm1); strcpy(atcm1,atm2);
  }
  
  /* set up DEFAULTS */
  
  rmax_dist = RMAX_DIST;
  rmin_dist = RMIN_DIST;
  rres_dist = RRES_DIST;
  *wlrc = *vlrc = 0.;
  
  strcpy(dict[0],"atom1");        strcpy(dict[1],"atom2"); 
  strcpy(dict[2],"pot_type");     strcpy(dict[3],"min_dist");
  strcpy(dict[4],"max_dist");     strcpy(dict[5],"res_dist");
  strcpy(dict[6],"eps");          strcpy(dict[7],"sig");
  strcpy(dict[8],"awill");        strcpy(dict[9],"bwill");
  strcpy(dict[10],"c6");          strcpy(dict[11],"c8");
  strcpy(dict[12],"c9");          strcpy(dict[13],"c10");
  strcpy(dict[14],"yamma");       strcpy(dict[15],"drmexp");
  strcpy(dict[16],"gamma");       strcpy(dict[17],"c12");         
  strcpy(dict[18],"d10");         strcpy(dict[19],"tabfile");
  
  /* loop until end of linklist keywords */
  nline = ifound = 0;
  meta_now=meta_key;
  
  for(imkey=0;imkey<nmeta_key && ifound==0;imkey++,meta_now=meta_now->next){
    if(!strncasecmp(META_KEY_INTER,meta_now->meta,strlen(META_KEY_INTER))){
      for(i=0;i<num_dict;i++) num_found[i] = 0;
      num_keys = meta_now->nkeys;
      key_now=meta_now->keys;
      for(i=0;i<num_keys;i++,key_now=key_now->next){
        icase = get_dict_num(num_dict,dict,key_now->name);
        if(icase == -1){
          print_dict(key_now->name,num_dict,dict); exit(1);
        }
        num_found[icase]++;
        switch (icase) {
           case 0:strcpy(atom1,key_now->arg); break;	/* atom1 */
           case 1:strcpy(atom2,key_now->arg); break;	/* atom2 */
           default:
             print_dict(key_now->arg,num_dict,dict); exit(1); break;
        }
        if(num_found[0]==1 && num_found[1]==1) break;
      }
      /* found inter_parm and got keys
         see if match the ones we are looking for*/
      if(strcasecmp(atom1,atom2)>0){
        strcpy(tword,atom1); strcpy(atom1,atom2); strcpy(atom2,tword);
      }
      
      if(!strcasecmp(atom1,atcm1) && !strcasecmp(atom2,atcm2)){
        /* we have match so we assign the rest of the keys */
        num_keys = meta_now->nkeys;
        key_now=meta_now->keys;
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        for(i=0;i<num_keys;i++,key_now=key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->name);
          if(icase == -1){
            print_dict(key_now->name,num_dict,dict); exit(1);
          }
          num_found[icase]++;
          switch (icase) {
	  case 0:	case 1: break;	/* atom1 and atom2 already done */
	  case 2:strcpy(pot_type,key_now->arg); break;/* pot_type */
	  case 3:sscanf(key_now->arg,"%lg",&rmin_dist);break;/*min*/
	  case 4:sscanf(key_now->arg,"%lg",&rmax_dist);break;/*max*/
	  case 5:sscanf(key_now->arg,"%lg",&rres_dist);break;/*res*/
	  case 6:sscanf(key_now->arg,"%lg",&eps_now);  break;/* eps */
	  case 7:sscanf(key_now->arg,"%lg",&sig_now);  break;/* sig */
	  case 8:sscanf(key_now->arg,"%lg",&awill);    break;/* awill*/
	  case 9:sscanf(key_now->arg,"%lg",&bwill);    break;/* bwill*/
	  case 10:sscanf(key_now->arg,"%lg",&c6);      break;/* c6 */
	  case 11:sscanf(key_now->arg,"%lg",&c8);      break;/* c8 */
	  case 12:sscanf(key_now->arg,"%lg",&c9);      break;/* c9 */
	  case 13:sscanf(key_now->arg,"%lg",&c10);     break;/* c10 */
	  case 14:sscanf(key_now->arg,"%lg",&yamma);   break;/*yamma*/
	  case 15:sscanf(key_now->arg,"%lg",&drmexp);  break;/*drmexp*/
	  case 16:sscanf(key_now->arg,"%lg",&gamma);   break;/*gamma*/
	  case 17:sscanf(key_now->arg,"%lg",&c12);     break;/*c12*/
	  case 18:sscanf(key_now->arg,"%lg",&d10);     break;/*d10*/
	  case 19:sscanf(key_now->arg,"%s",fname);     break;/*tabname*/
	  default:
	    print_dict(key_now->arg,num_dict,dict); exit(1); break;
          }
        }
        if(!(!strcasecmp(POT_TYPE_NULL,pot_type) ||
	     !strcasecmp(POT_TYPE_LJ,pot_type) || 
	     !strcasecmp(POT_TYPE_LJ64,pot_type) || 
	     !strcasecmp(POT_TYPE_LJ96,pot_type) || 
	     !strcasecmp(POT_TYPE_WILL,pot_type) || 
	     !strcasecmp(POT_TYPE_TAB,pot_type) || 
	     !strcasecmp(POT_TYPE_HYDBND,pot_type) || 
	     !strcasecmp(POT_TYPE_AZIZ,pot_type) || 
	     !strcasecmp(POT_TYPE_DZUG,pot_type))){
          fprintf(stderr,"ERROR: pot_type %s not implimented\n",pot_type);
          fprintf(stderr,"in file %s ",filename);
          fprintf(stderr,"between atoms %s and %s\n",atom1,atom2);
          exit(1);
        }
#ifdef DEBUG
        fprintf(stdout,"pot_type between %s and %s is %s\n",
           atom1,atom2,pot_type);
#endif
        /* check for required key words */
        for(i=0;i<3;i++){
          if(num_found[i] != 1)
            err_found(dict[i],num_found[i],filename,nline); 
        }
        if(!strcasecmp(POT_TYPE_NULL,pot_type)){
          *inull = 1;
          return;
        } else {
          *inull = 0;
        }
        if(!strcasecmp(POT_TYPE_LJ,pot_type) || 
	   !strcasecmp(POT_TYPE_LJ64,pot_type) ||
	   !strcasecmp(POT_TYPE_LJ96,pot_type)){
          for(i=6;i<8;i++){
            if(num_found[i] != 1) 
              err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(!strcasecmp(POT_TYPE_WILL,pot_type)){
          for(i=8;i<13;i++){
            if(num_found[i] != 1) 
              err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(!strcasecmp(POT_TYPE_AZIZ,pot_type)){
          for(i=8;i<=16;i++){
            if(num_found[i] != 1) 
              err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(!strcasecmp(POT_TYPE_HYDBND,pot_type)){
          for(i=17;i<19;i++){
            if(num_found[i] != 1) 
              err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        if(!strcasecmp(POT_TYPE_TAB,pot_type)){
          for(i=19;i<20;i++){
            if(num_found[i] != 1)
              err_found(dict[i],num_found[i],filename,nline);
          }
        }
        if(!strcasecmp(POT_TYPE_DZUG,pot_type)){
          for(i=6;i<8;i++){ // right now we are just using sig and eps, like the LJ pot.
            if(num_found[i] != 1)
              err_found(dict[i],num_found[i],filename,nline);
          }
        }
        ifound = 1;
        if(rres_dist<rmax_dist){
          *rmins = rmin_dist;  	    
          *rmaxs = rres_dist;
          *rminl = rres_dist-inter->rheal;
          *rmaxl = rmax_dist;
        } else {
          *rmins = rmin_dist;  	    
          *rmaxs = rmax_dist;
          *rmaxl = rmax_dist;  	    
          *rminl = rmax_dist+2.*inter->skin_ter;
          fprintf(stderr,"WARNING: respa cutoff not being used!");
          fprintf(stderr," between atom types %s and %s\n",atm1,atm2);
          rres_dist = rmax_dist+2.*inter->rheal;
        }
        inter->dx2tab[indx]=((double)inter->ntable-3.)/
          (rmax_dist*rmax_dist-rmin_dist*rmin_dist);
	 
        /* cutoff */
        if(*rminl-inter->skin_ter < rmin_dist) {
          fprintf(stderr,"ERROR:Respa distance, healing length, skin and");
          fprintf(stderr," minimum distance are inconsistant\n");
          fprintf(stderr,"on atoms %s %s with values %g %g %g %g\n",
             atm1,atm2,rres_dist,inter->rheal,
             inter->skin_ter,rmin_dist);
          exit(1);
        } else {
          rch = (rres_dist-inter->rheal);
        }
        if (inter->rheal == 0.) {
          rhl = 1.;
        } else{
          rhl = 1./inter->rheal;
        }
        if(rmin_dist > rmax_dist){
          md_error("short cutoff can't be longer then long cutoff");
        }
        /* set long range correction */
        if(!strcasecmp(POT_TYPE_LJ,pot_type)){
          *vlrc = (4.*eps_now*pow(sig_now,6.)*
             ((1./9.)*pow(sig_now,6.)*pow(1./(rmax_dist),9.)-
                (1./3.)*pow(1./(rmax_dist),3.)));
          *wlrc = (8.*eps_now*pow(sig_now,6.)*
             (pow(1./(rmax_dist),3.) -
                (2./3.)*pow(sig_now,6.)*pow(1./(rmax_dist),9.)));
        }
        /* set long range correction *CHECK* fudge */
        if(!strcasecmp(POT_TYPE_LJ96,pot_type)){
          *vlrc = ((27./4.)*eps_now*pow(sig_now,6.)*
             ((1./6.)*pow(sig_now,3.)*pow(1./(rmax_dist),6.)-
                (1./3.)*pow(1./(rmax_dist),3.)));
          *wlrc = (2.*(27./4.)*eps_now*pow(sig_now,6.)*
             (2.*pow(1./(rmax_dist),3.) -
                (3./2.)*pow(sig_now,3.)*pow(1./(rmax_dist),6.)));
        }
	/* no correction applied at present
	   if(!strcasecmp(POT_TYPE_LJ64,pot_type)){
	   }
	*/
        if(!strcasecmp(POT_TYPE_WILL,pot_type)){
          sm = bwill*rmax_dist;
          *vlrc=((awill*(sm*sm+2*sm+2.)/pow(bwill,3.))*exp(-bwill*rmax_dist)-
             c10*pow(rmax_dist,-7.)/7. - c9*pow(rmax_dist,-6.)/6. -
             c8*pow(rmax_dist,-5.)/5. - c6*pow(rmax_dist,-3.)/3.);   
          
          *wlrc =(-(awill*exp(-bwill*rmax_dist)*
             (sm*sm*sm+3.*sm*sm+6.*sm+6.)/pow(bwill,3.)) +
             c10*pow(rmax_dist,-7.)*10./7.+c9*pow(rmax_dist,-6.)*9./6.+
             c8*pow(rmax_dist,-5.)*8./5.+c6*pow(rmax_dist,-3.)*6./3.);
          
        }
        if(!strcasecmp(POT_TYPE_AZIZ,pot_type)){
          /* using only dispersion part of potential without the switch */
          *vlrc = -(c10/(7.*pow(rmax_dist,7.))+c8/(5.*pow(rmax_dist,5.))
             +c6/(3.*pow(rmax_dist,3.)));
          *wlrc =(c10*pow(rmax_dist,-7.)*10./7.+c9*pow(rmax_dist,-6.)*9./6.+
             c8*pow(rmax_dist,-5.)*8./5.+c6*pow(rmax_dist,-3.)*6./3.);
        }
        if(!strcasecmp(POT_TYPE_HYDBND,pot_type)){
          *vlrc = c12/9.*pow(1./(rmax_dist),9.)-
            d10/7.*pow(1./(rmax_dist),7.);
          *wlrc = -c12*pow(1./(rmax_dist),9.)*12./9.+
            d10*pow(1./(rmax_dist),7.)*10./7.;
        }
        /* no correction applied for vlrc and wlrc at present
	      if(!strcasecmp(POT_TYPE_DZUG,pot_type))
	      {
	      }
        */
        *vlrc *= scaleeps;
        *wlrc *= scaleeps;
        /* tabulate the energy and forces */
        if(!strcasecmp(POT_TYPE_TAB,pot_type)) {
          printf("reading POT_TYPE_TAB file name %s\n",fname);
          if((ftab = fopen(fname,"r")) == NULL){
            sprintf(line,"can't open tabpot file (%s)",fname);
            md_error(line);
          }
        }
        for(i=0;i<inter->ntable;i++){
          xt = i/(inter->dx2tab[indx]) + (*rmins)*(*rmins);
          if(xt<0){
            sprintf(line,"Interaction table values are negative!\n");
            sprintf(line,"%s xt = %g (%s %s %s)",line,
               xt,atom1,atom2,pot_type);
            md_error(line);
          }
          if(!strcasecmp(POT_TYPE_TAB,pot_type)) {
            func_tab_pot(ftab,sqrt(xt),&f,&df,&d2f);
	    /*
	      printf("%e %e %e %e\n",sqrt(xt),f,df,d2f);
	    */
          }
          if(!strcasecmp(POT_TYPE_LJ,pot_type))
            func_lj(sqrt(xt),&f,&df,&d2f,eps_now,sig_now);
          if(!strcasecmp(POT_TYPE_LJ96,pot_type))
            func_lj96(sqrt(xt),&f,&df,&d2f,eps_now,sig_now);
          if(!strcasecmp(POT_TYPE_LJ64,pot_type))
            func_lj64(sqrt(xt),&f,&df,&d2f,eps_now,sig_now);
          if(!strcasecmp(POT_TYPE_WILL,pot_type))
            func_will(sqrt(xt),&f,&df,&d2f,awill,bwill,c6,c8,c10);
          if(!strcasecmp(POT_TYPE_AZIZ,pot_type))
            func_aziz(sqrt(xt),&f,&df,&d2f,awill,bwill,c6,c8,c9,c10,
               yamma,drmexp,gamma);
          if(!strcasecmp(POT_TYPE_HYDBND,pot_type))
            func_hydbnd(sqrt(xt),&f,&df,&d2f,c12,d10);
          if(!strcasecmp(POT_TYPE_DZUG,pot_type))
            func_dzug(sqrt(xt),&f,&df,&d2f,eps_now,sig_now);
          
          inter->vtab[indx][i] = f*scaleeps;
          inter->dvtab[indx][i] = df*scaleeps; 
          inter->d2vtab[indx][i] = d2f*scaleeps;
          if(xt<=rch*rch){
            inter->switchs[indx][i] = 1.0; 
          }
          if(xt>rch*rch && xt < rres_dist*rres_dist){
            br = (sqrt(xt)-rch)*rhl; sm = (1.+ br*br*(2.*br-3.));
            inter->switchs[indx][i] = sm;         
          }
          if(xt >= rres_dist*rres_dist){
            inter->switchs[indx][i] = 0.0;       
          }
        }
        if(!strcasecmp(POT_TYPE_TAB,pot_type)) fclose(ftab);
      }
    }
  }
  if(!ifound){
    fprintf(stderr,
	    "ERROR: Inter interaction between %s and %s not found is %s\n"
	    ,atm1,atm2,filename);
    exit(1);
  }
}
/*-----------------------------------------------------------------------*/
void get_bondparm(int nmeta_key,META_KEY *meta_key,
   char *filename,WORD atm1,WORD atm2,WORD bond_label,
   int *it_now,double *eq_now,double *fk_now,double *bm_now,int *ifix)
{
  int i,ifound,nline,num_keys,icase,ifx,imkey;
  int num_dict=8,num_found[8],ptype;
  WORD dict[8];
  WORD atom1,atom2,pot_type,tword,tbo;
  META_KEY *meta_now;
  KEYS *key_now;
  double eq,fk,bm;

  strcpy(dict[0],"atom1");     strcpy(dict[1],"atom2");
  strcpy(dict[2],"pot_type");  strcpy(dict[3],"eq");
  strcpy(dict[4],"fk");        strcpy(dict[5],"label");
  strcpy(dict[6],"ifix");      strcpy(dict[7],"bmorse");

  /* loop until end of linklist keywords */
  nline = ifound = 0;
  meta_now=meta_key;
  
  for(imkey=0;imkey<nmeta_key && ifound==0;imkey++,meta_now=meta_now->next){
    if(!strncasecmp(META_KEY_BOND,meta_now->meta,strlen(META_KEY_BOND))){
      for(i=0;i<num_dict;i++) num_found[i] = 0;
      ifx=0;
      num_keys = meta_now->nkeys;
      key_now=meta_now->keys;
      for(i=0;i<num_keys;i++,key_now=key_now->next){
        icase = get_dict_num(num_dict,dict,key_now->name);
        if(icase == -1){
          print_dict(key_now->name,num_dict,dict); exit(1);
        }
        num_found[icase]++;
        switch (icase) {
          
           case 0: strcpy(atom1,key_now->arg); break;/* atom1 */
           case 1: strcpy(atom2,key_now->arg); break;/* atom2 */
           case 2: strcpy(pot_type,key_now->arg); break;/* pot_type */
           case 3: sscanf(key_now->arg,"%lg",&eq); break;/* eq */
           case 4: sscanf(key_now->arg,"%lg",&fk); break;/* fk */
           case 5: sscanf(key_now->arg,"%s",tbo);  break;/* label */
           case 6: sscanf(key_now->arg,"%d",&ifx); break;/* ifix */
           case 7: sscanf(key_now->arg,"%lg",&bm); break;/* bmorse */
           default:
             print_dict(key_now->name,num_dict,dict);
             exit(1); break;
        }
      }
      /* check to see if we found required keys */
      for(i=0;i<=2;i++){
        if(num_found[i] != 1){
          err_found(dict[i],num_found[i],filename,nline); 
        }
      }
      /* if match the ones we are looking for*/
      if(!strcasecmp(POT_TYPE_HARM,pot_type)) ptype = 0;
      else if (!strcasecmp(POT_TYPE_MORSE,pot_type)) ptype = 1;
      else {
        fprintf(stderr,"ERROR: pot_type %s (in file %s not implimented\n"
           ,pot_type,filename);
        exit(1);
      }
      for(i=3;i<=4;i++){
        if(num_found[i] != 1){
          err_found(dict[i],num_found[i],filename,nline); 
        }
      }
      if(ptype==0) bm = 0.;
      if(ptype==1 && num_found[7]==0){
        err_found(dict[i],num_found[i],filename,nline);
      }
      if(num_found[5]==0) tbo[0] = NILL;
      if(strcasecmp(atom1,atom2)>0){
        strcpy(tword,atom1);  strcpy(atom1,atom2);  strcpy(atom2,tword);
      }
      if(!strcasecmp(atom1,atm1) && !strcasecmp(atom2,atm2) &&
         !strcasecmp(bond_label,tbo)){
        ifound = 1;  *eq_now = eq;  *fk_now = fk; *ifix = ifx;
        *it_now = ptype; *bm_now = bm;
      }
    }
  }
  if(!ifound){
    fprintf(stderr,"ERROR: Bonded interaction between %s and %s",atm1,atm2);
    fprintf(stderr," with type \"%s\" not found in %s\n",bond_label,filename);
    exit(1);
  }
}
/*-----------------------------------------------------------------------*/
void get_bondxparm(FILE *fp,char *filename,WORD atm1,WORD atm2,WORD atm3,
		   WORD box_type,double *eq1_now,double *eq2_now,
		   double *fk_now,int *ifix)
{
  char line[MAXLINELEN];
  int i,ifound,ich,nline,nch,num_keys,icase,ifx;
  int num_dict=9,num_found[9];
  WORD dict[9];
  WORD metakey,atom1,atom2,atom3,pot_type,tword,label;
  KEY key_root,*key_now;
  double eq1,eq2,fk,temp;

  strcpy(dict[0],"atom1");    strcpy(dict[1],"atom2");
  strcpy(dict[2],"atom3");    strcpy(dict[3],"pot_type");
  strcpy(dict[4],"eq1");      strcpy(dict[5],"eq2");
  strcpy(dict[6],"fk");       strcpy(dict[7],"label");
  strcpy(dict[8],"ifix");

  /* loop until end of file */
  rewind(fp);
  nline = nch = ifound = 0;
  while((ich = fgetc(fp)) != EOF && !(ifound)) {
    line[nch++] = (char)ich;
    if(ich == '~'){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      if(!strncasecmp(META_KEY_BONDX,metakey,strlen(META_KEY_BONDX))){
	for(i=0;i<num_dict;i++) num_found[i] = 0;
	ifx = 0;
	for(i=0,key_now=&key_root;i<num_keys;i++,key_now=key_now->next){
	  icase = get_dict_num(num_dict,dict,key_now->keyword);
	  if(icase==-1){
	    print_dict(key_now->keyword,num_dict,dict); exit(1);
	  }
	  num_found[icase]++;
	  switch (icase) {
	  case 0: strcpy(atom1,key_now->keyarg); break;	/* atom1 */
	  case 1: strcpy(atom2,key_now->keyarg); break;	/* atom2 */
	  case 2: strcpy(atom3,key_now->keyarg); break;	/* atom3 */
	  case 3: strcpy(pot_type,key_now->keyarg); break;/* pot_type */
	  case 4: sscanf(key_now->keyarg,"%lg",&eq1); break; /* eq1 */
	  case 5: sscanf(key_now->keyarg,"%lg",&eq2); break; /* eq2 */
	  case 6: sscanf(key_now->keyarg,"%lg",&fk); break; /* fk */
	  case 7: strcpy(label,key_now->keyarg); break; /* label */
	  case 8: sscanf(key_now->keyarg,"%d",&ifx); break; /* ifix */
	  default:
	    print_dict(key_now->keyword,num_dict,dict); exit(1); break;
	  }
	}
	/* check to see if we found required keys */
	for(i=0;i<=6;i++){
	  if(num_found[i] != 1){
	    err_found(dict[i],num_found[i],filename,nline); 
	  }
	}
	if(num_found[7]==0)label[0] = NILL;
	/* found inter_parm and got keys
	   see if match the ones we are looking for*/
	if(!strcasecmp(POT_TYPE_BONDX,pot_type)){
	  if(strcasecmp(atom1,atom3)>0){
	    strcpy(tword,atom1); strcpy(atom1,atom3); strcpy(atom3,tword);
	    temp=eq1;eq1=eq2;eq2=temp;
	  }
	  if(!strcasecmp(atom1,atm1) && !strcasecmp(atom2,atm2) &&
	     !strcasecmp(atom3,atm3) && !strcasecmp(label,box_type)){
	    ifound = 1; *eq1_now=eq1; *eq2_now=eq2; *fk_now = fk; *ifix = ifx;
	  }
	} else {
	  fprintf(stderr,"ERROR: pot_type %s (in file %s) not implimented\n"
		  ,pot_type,filename);
	  exit(1);
	}
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch=0;nline++;}
  }
  if(!ifound){
    fprintf(stderr,
	    "ERROR: Cross bonded interaction between %s %s %s not found (%s)\n"
	    ,atm1,atm2,atm3,filename);
    exit(1);
  }
}
/*-----------------------------------------------------------------------*/
void get_bendparm(FILE *fp,char *filename,WORD atm1,WORD atm2,WORD atm3,
		  WORD bend_label,double *eq_now,double *fk_now,int *ifix,
		  int ibend,int jbend,int kbend)
{
  char line[MAXLINELEN];
  int i,ifound,ich,nline,nch,num_keys,icase,ifx;
  int num_dict=8,num_found[8];
  WORD dict[8];
  WORD metakey,atom1,atom2,atom3,pot_type,tword,label;
  KEY key_root,*key_now;
  double eq,fk;

  strcpy(dict[0],"atom1");   strcpy(dict[1],"atom2");
  strcpy(dict[2],"atom3");   strcpy(dict[3],"pot_type");
  strcpy(dict[4],"eq");      strcpy(dict[5],"fk");
  strcpy(dict[6],"label");   strcpy(dict[7],"ifix");

  /* loop until end of file */
  rewind(fp);
  nline = nch = ifound = 0;
  while( (ich = fgetc(fp)) != EOF && !(ifound)) {
    line[nch++] = (char)ich;
    if(ich == '~'){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      if(!strncasecmp(META_KEY_BEND,metakey,strlen(META_KEY_BEND))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        ifx = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase==-1){
            print_dict(key_now->keyword,num_dict,dict); exit(1);
          }
          num_found[icase]++;
          switch (icase) {
             case 0:strcpy(atom1,key_now->keyarg);break;		/* atom1 */
             case 1:strcpy(atom2,key_now->keyarg);break;		/* atom2 */
             case 2:strcpy(atom3,key_now->keyarg);break;		/* atom3 */
             case 3:strcpy(pot_type,key_now->keyarg);break;	/* pot_type */
             case 4:sscanf(key_now->keyarg,"%lg",&eq);break;	/* eq */
             case 5:sscanf(key_now->keyarg,"%lg",&fk);break;	/* fk */
             case 6:strcpy(label,key_now->keyarg);break;	/* label */
             case 7:sscanf(key_now->keyarg,"%d",&ifx);break;	/* ifix */
             default:
               print_dict(key_now->keyword,num_dict,dict); exit(1);break;
          }
        }
        if(num_found[6]==0) label[0] = NILL;
        /* check to see if we found required keys */
        for(i=0;i<=3;i++){
          if(num_found[i] != 1){
            err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        /* found inter_parm and got keys
           see if match the ones we are looking for*/
        if(!strcasecmp(POT_TYPE_HARM,pot_type)){
          for(i=4;i<=5;i++){
            if(num_found[i] != 1){
              err_found(dict[i],num_found[i],filename,nline); 
            }
          }
          if(strcasecmp(atom1,atom3)>0){
            strcpy(tword,atom1); strcpy(atom1,atom3);  strcpy(atom3,tword);
          }
          if(!strcasecmp(atom1,atm1) && !strcasecmp(atom2,atm2) &&
             !strcasecmp(atom3,atm3) && !strcasecmp(bend_label,label)){
            ifound = 1; *eq_now = eq; *fk_now = fk; *ifix = ifx;
          }
        } else {
          fprintf(stderr,"ERROR: pot_type %s (in file %s) not implimented\n"
             ,pot_type,filename);
          exit(1);
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch=0;nline++;}
  }
  if(!ifound){
    fprintf(stderr,
	    "ERROR: Bended interaction between %s (%d) %s (%d) %s (%d) (%s) not found in %s\n"
	    ,atm1,ibend,atm2,jbend,atm3,kbend,bend_label,filename);
      exit(1);   
  }
}
/*-----------------------------------------------------------------------*/
void get_torsparm(FILE *fp,char *filename,WORD atm1,WORD atm2,WORD atm3,
		  WORD atm4,WORD ttors,int *int_pot,
		  double *eq_now,double *fk_now,VPHI vphi_now,int *ifix)
{
  char line[MAXLINELEN];
  int i,ifound,ich,nline,nch,num_keys,icase,ifx;
  int num_dict=16,num_found[16];
  WORD dict[16];
  WORD metakey,atom1,atom2,atom3,atom4,pot_type,tto;
  KEY key_root,*key_now;
  double eq,fk;
  VPHI vphi_temp;
  
  strcpy(dict[0],"atom1");    strcpy(dict[1],"atom2");
  strcpy(dict[2],"atom3");    strcpy(dict[3],"atom4");
  strcpy(dict[4],"pot_type"); strcpy(dict[5],"eq");
  strcpy(dict[6],"fk");       strcpy(dict[7],"p0");
  strcpy(dict[8],"p1");       strcpy(dict[9],"p2");
  strcpy(dict[10],"p3");      strcpy(dict[11],"p4");
  strcpy(dict[12],"p5");      strcpy(dict[13],"p6");
  strcpy(dict[14],"label");   strcpy(dict[15],"ifix");
  
  /* loop until end of file */
  rewind(fp);
  nline = nch = ifound = 0;
  while( (ich = fgetc(fp)) != EOF && !(ifound)) {
    line[nch++] = (char)ich;
    if(ich == '~'){
      get_meta_key(fp,filename,line,&nch,&nline,metakey,&key_root,&num_keys);
      if(!strncasecmp(META_KEY_TORS,metakey,strlen(META_KEY_TORS))){
        for(i=0;i<num_dict;i++) num_found[i] = 0;
        ifx = 0;
        for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
          icase = get_dict_num(num_dict,dict,key_now->keyword);
          if(icase == -1) print_dict(key_now->keyword,num_dict,dict); 
          num_found[icase]++;
          switch (icase) {
             case 0:strcpy(atom1,key_now->keyarg); break;		/* atom1 */
             case 1:strcpy(atom2,key_now->keyarg); break;		/* atom2 */
             case 2:strcpy(atom3,key_now->keyarg); break;		/* atom3 */
             case 3:strcpy(atom4,key_now->keyarg); break;		/* atom4 */
             case 4:strcpy(pot_type,key_now->keyarg);break;	/* pot_type */
             case 5:sscanf(key_now->keyarg,"%lg",&eq);break;	/* eq */
             case 6:sscanf(key_now->keyarg,"%lg",&fk);break;	/* fk */
             /*p0-p6*/
             case 7: case 8: case 9: case 10: case 11: case 12: case 13:
               sscanf(key_now->keyarg,"%lg",&(vphi_temp[icase-7]));
               break;
             case 14: sscanf(key_now->keyarg,"%s",tto);            /* label */
               break;
             case 15:sscanf(key_now->keyarg,"%d",&ifx);break;	/* ifix */
             default:
               print_dict(key_now->keyword,num_dict,dict); exit(1); break;
          }
        }
        for(i=0;i<=4;i++){
          if(num_found[i] != 1){
            err_found(dict[i],num_found[i],filename,nline); 
          }
        }
        
        if(num_found[14] == 0) tto[0] = NILL;
        
        if((!strcasecmp(atom1,atm1) || !strcasecmp(atom1,WILDCARD)) &&
           (!strcasecmp(atom2,atm2) || !strcasecmp(atom2,WILDCARD)) &&
           (!strcasecmp(atom3,atm3) || !strcasecmp(atom3,WILDCARD)) &&
           (!strcasecmp(atom4,atm4) || !strcasecmp(atom4,WILDCARD)) &&
           (!strcasecmp(ttors,tto))){
          ifound=1;
        } else if((!strcasecmp(atom1,atm4) || !strcasecmp(atom1,WILDCARD)) &&
           (!strcasecmp(atom2,atm3) || !strcasecmp(atom2,WILDCARD)) &&
           (!strcasecmp(atom3,atm2) || !strcasecmp(atom3,WILDCARD)) &&
           (!strcasecmp(atom4,atm1) || !strcasecmp(atom4,WILDCARD)) &&
           (!strcasecmp(ttors,tto))){
          ifound=1;
        }
        
#ifdef DEBUG
        printf("%s %s %s %s %s vs %s %s %s %s %s\n",atm1,atm2,atm3,atm4,ttors,
           atom1,atom2,atom3,atom4,tto);
#endif
        if(ifound){
          if(strcasecmp(POT_TYPE_HARM,pot_type) && 
             strcasecmp(POT_TYPE_POWER,pot_type)){
            fprintf(stderr,"ERROR: pot_type %s (in file %s) not implimented\n"
               ,pot_type,filename);
            exit(1);
          }
          if(!strcasecmp(POT_TYPE_HARM,pot_type)){
            for(i=5;i<=6;i++){
              if(num_found[i] != 1){
                err_found(dict[i],num_found[i],filename,nline); 
              }
            }
            *int_pot = 0; *eq_now = eq; *fk_now = fk; *ifix = ifx;
          }
          if(!strcasecmp(POT_TYPE_POWER,pot_type)){
            for(i=7;i<=13;i++){
              if(num_found[i] != 1){
#ifdef ALLPOWER
                err_found(dict[i],num_found[i],filename,nline); 
#else
                vphi_temp[i-7] = 0.;
#endif
              }
            }
            *int_pot = 1;
            for(i=0;i<MAXPOWER_TORS;i++)vphi_now[i] = vphi_temp[i];
            *ifix = ifx;
          }
        }
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch=0;nline++;}
  }

  if(!ifound){
    fprintf(stderr,"ERROR: Torsion interaction between %s %s %s %s ",
	    atm1,atm2,atm3,atm4);
    fprintf(stderr,"with label \"%s\" not found in %s\n",
	    ttors,filename);
    exit(1);
  }
}
/*-----------------------------------------------------------------------*/
void get_onfo_parm(FILE *fp,char *filename,WORD atm1,WORD atm2,double *rmin,
		   double *rmax,double *dx2tab,double *vtab,double *dvtab,
		   double *d2vtab,int ntable)
{
  char line[MAXLINELEN];
  int i,ifound,ich,nline,nch,num_keys,icase;
  int num_dict=20,num_found[20];
  WORD dict[20];
  WORD metakey,atcm1,atcm2,atom1,atom2,pot_type,tword;
  KEY key_root,*key_now;
  double eps_now,sig_now,rmin_dist,rmax_dist;
  double awill,bwill,c6,c8,c9,c10,gamma,yamma,drmexp,c12,d10;
  double xt,f,df,d2f;
  char fname[MAXLINELEN];
  FILE *ftab;

  if(strcasecmp(atm1,atm2)<=0){
    strcpy(atcm1,atm1);  strcpy(atcm2,atm2);
  } else {
    strcpy(atcm2,atm1);  strcpy(atcm1,atm2);
  }

  strcpy(dict[0],"atom1");     strcpy(dict[1],"atom2");
  strcpy(dict[2],"pot_type");  strcpy(dict[3],"min_dist");
  strcpy(dict[4],"max_dist");  strcpy(dict[5],"res_dist");
  strcpy(dict[6],"eps");       strcpy(dict[7],"sig");
  strcpy(dict[8],"awill");     strcpy(dict[9],"bwill");
  strcpy(dict[10],"c6");       strcpy(dict[11],"c8");
  strcpy(dict[12],"c9");       strcpy(dict[13],"c10");
  strcpy(dict[14],"yamma");    strcpy(dict[15],"drmexp");   
  strcpy(dict[16],"gamma");    strcpy(dict[17],"c12");      
  strcpy(dict[18],"d10");      strcpy(dict[19],"tabfile");

  rewind(fp);
  nline = nch = ifound = 0;
  while( (ich = fgetc(fp)) != EOF && !(ifound)) {
    line[nch++] = (char)ich;
    if(ich == '~'){
      get_meta_key(fp,filename,
		   line,&nch,&nline,metakey,&key_root,&num_keys);
      
      if(!strncasecmp(META_KEY_ONFO,metakey,strlen(META_KEY_ONFO))){
	for(i=0;i<num_dict;i++) num_found[i] = 0;
	for(i=0,key_now=&key_root;i<num_keys;i++,key_now = key_now->next){
	  icase = get_dict_num(num_dict,dict,key_now->keyword);
	  if(icase == -1){
	    print_dict(key_now->keyword,num_dict,dict); exit(1);
	  }
	  num_found[icase]++;
	  switch (icase) {
	  case 0:strcpy(atom1,key_now->keyarg); break;	/* atom1 */
     	  case 1:strcpy(atom2,key_now->keyarg); break; /* atom2 */
	  case 2:strcpy(pot_type,key_now->keyarg); break;/* pot_type */
	  case 3:sscanf(key_now->keyarg,"%lg",&rmin_dist); break;/*min_dist*/
	  case 4:sscanf(key_now->keyarg,"%lg",&rmax_dist);  break;/* max_dist*/
	  case 5:break; /* res_dist ignored in 1-4 interactions*/
	  case 6:sscanf(key_now->keyarg,"%lg",&eps_now);break;	/* eps */
	  case 7:sscanf(key_now->keyarg,"%lg",&sig_now);break;	/* sig */
	  case 8:sscanf(key_now->keyarg,"%lg",&awill);  break;	/* awill */
	  case 9:sscanf(key_now->keyarg,"%lg",&bwill);  break;	/* bwill */
	  case 10:sscanf(key_now->keyarg,"%lg",&c6);    break;	/* c6 */
	  case 11:sscanf(key_now->keyarg,"%lg",&c8);    break;	/* c8 */
	  case 12:sscanf(key_now->keyarg,"%lg",&c9);    break;	/* c9 */
	  case 13:sscanf(key_now->keyarg,"%lg",&c10);   break;	/* c10 */
	  case 14:sscanf(key_now->keyarg,"%lg",&yamma); break;	/* yamma */
	  case 15:sscanf(key_now->keyarg,"%lg",&drmexp);break;	/* drmexp */
	  case 16:sscanf(key_now->keyarg,"%lg",&gamma); break;	/* gamma */
	  case 17:sscanf(key_now->keyarg,"%lg",&c12);   break;	/* c12 */
	  case 18:sscanf(key_now->keyarg,"%lg",&d10);   break;	/* d10 */
    case 19:sscanf(key_now->keyarg,"%s",fname);   break;/*tabname*/
	  default:
	    print_dict(key_now->keyword,num_dict,dict);  exit(1);   break;
	  }
	}
	/* found inter_parm and got keys
	   see if match the ones we are looking for*/
	if(strcasecmp(atom1,atom2)>0){
	  strcpy(tword,atom1); strcpy(atom1,atom2); strcpy(atom2,tword);
	}

	if(!strcasecmp(atom1,atcm1) && !strcasecmp(atom2,atcm2)){
	  if(!(!strcasecmp(POT_TYPE_LJ,pot_type) ||
	       !strcasecmp(POT_TYPE_LJ64,pot_type) || 
	       !strcasecmp(POT_TYPE_LJ96,pot_type) || 
	       !strcasecmp(POT_TYPE_WILL,pot_type) || 
	       !strcasecmp(POT_TYPE_HYDBND,pot_type) || 
	       !strcasecmp(POT_TYPE_AZIZ,pot_type))){
	    fprintf(stderr,"ERROR: pot_type %s not implimented\n",pot_type);
	    fprintf(stderr,"in file %s at line %d\n",filename,nline);
	    fprintf(stderr,"between atoms %s and %s\n",atom1,atom2);
	    exit(1);
	  }
#ifdef DEBUG
	  fprintf(stdout,"pot_type (1-4) between %s and %s is %s\n",
		  atom1,atom2,pot_type);
#endif
	  /* check for required key words */
	  for(i=0;i<5;i++){
	    if(num_found[i] != 1)
	      err_found(dict[i],num_found[i],filename,nline); 
	  }
	  if(!strcasecmp(POT_TYPE_LJ,pot_type) || 
	     !strcasecmp(POT_TYPE_LJ64,pot_type) ||
	     !strcasecmp(POT_TYPE_LJ96,pot_type)){
	    for(i=6;i<8;i++){
	      if(num_found[i] != 1) 
		err_found(dict[i],num_found[i],filename,nline); 
	    }
	  }
	  if(!strcasecmp(POT_TYPE_WILL,pot_type)){
	    for(i=8;i<13;i++){
	      if(num_found[i] != 1) 
		err_found(dict[i],num_found[i],filename,nline); 
	    }
	  }
	  if(!strcasecmp(POT_TYPE_AZIZ,pot_type)){
	    for(i=8;i<=16;i++){
	      if(num_found[i] != 1) 
		err_found(dict[i],num_found[i],filename,nline); 
	    }
	  }
	  if(!strcasecmp(POT_TYPE_HYDBND,pot_type)){
	    for(i=17;i<19;i++){
	      if(num_found[i] != 1) 
		err_found(dict[i],num_found[i],filename,nline); 
	    }
	  }

	  ifound = 1;	  *rmin = rmin_dist;  *rmax = rmax_dist;
	  *dx2tab=((double)ntable-3.)/((*rmax)*(*rmax)- (*rmin)*(*rmin));
	  
	  if(rmin_dist > rmax_dist){
	    md_error("short cutoff can't be longer then long cutoff");
	  }
          if(!strcasecmp(POT_TYPE_TAB,pot_type)) {
            if((fp = fopen(fname,"r")) == NULL){
              sprintf(line,"can't open tabpot file (%s)",fname);
              md_error(line);
            }
          }

	  /* tabulate the energy and forces */
	  for(i=0;i<ntable;i++){
	    xt = i/(*dx2tab) + (*rmin)*(*rmin);
            if(!strcasecmp(POT_TYPE_TAB,pot_type))
              func_tab_pot(ftab,sqrt(xt),&f,&df,&d2f);
	    if(!strcasecmp(POT_TYPE_LJ,pot_type))
	      func_lj(sqrt(xt),&f,&df,&d2f,eps_now,sig_now);
	    if(!strcasecmp(POT_TYPE_LJ64,pot_type))
	      func_lj64(sqrt(xt),&f,&df,&d2f,eps_now,sig_now);
	    if(!strcasecmp(POT_TYPE_LJ96,pot_type))
	      func_lj96(sqrt(xt),&f,&df,&d2f,eps_now,sig_now);
	    if(!strcasecmp(POT_TYPE_WILL,pot_type))
	      func_will(sqrt(xt),&f,&df,&d2f,awill,bwill,c6,c8,c10);
	    if(!strcasecmp(POT_TYPE_AZIZ,pot_type))
	      func_aziz(sqrt(xt),&f,&df,&d2f,awill,bwill,c6,c8,c9,c10,
			yamma,drmexp,gamma);
	    if(!strcasecmp(POT_TYPE_HYDBND,pot_type))
	      func_hydbnd(sqrt(xt),&f,&df,&d2f,c12,d10);
	    
	    vtab[i] = f; dvtab[i] = df; d2vtab[i] = d2f;
	  }
          if(!strcasecmp(POT_TYPE_TAB,pot_type)) fclose(ftab);
	}
      }
      free_key(num_keys,&key_root);
    }
    if(ich == '\n'){nch=0;nline++;}
  }
  if(!ifound){
    fprintf(stderr,
	    "ERROR: Onefour interaction between %s and %s not found is %s\n"
	    ,atm1,atm2,filename);
    exit(1);
  }
}
/*-----------------------------------------------------------------------*/
