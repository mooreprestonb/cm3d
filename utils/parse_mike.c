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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "system.h"

static int nmeta;

KEY *key_now;
META_KEY *meta_now;

/* #define DEBUG */

/*--------------------------------------------*/
void addmeta(char *name)
{
  nmeta++;
  meta_now->meta = strdup(name);
  meta_now->next = malloc(sizeof(META_KEY));

  meta_now = meta_now->next;

  meta_now->nkeys = 0;
  meta_now->next = NULL;  
  key_now = meta_now->keys = malloc(sizeof(KEY));
}
/*--------------------------------------------*/
void addkey(char *name)
{
  int i;

  meta_now->nkeys++;
  key_now->name = strdup(name);
  key_now->next = malloc(sizeof(KEY));
  key_now = key_now->next;
}
/*--------------------------------------------*/
void addarg(char *s,int itype)
{
  int i;

  key_now->itype = itype;
  
  itype = 2; /* return all strings */
#ifdef DEBUG
  printf("%d %s\n",key_now->itype,s); 
#endif
  switch(itype){
  case 0:
    key_now->arg.ival = (int)atof(s);
    break;
  case 1:
    key_now->arg.dval = atof(s);
    break;
  case 2:
    key_now->arg.sval = strdup(s);
    break;
  }
}
/*--------------------------------------------*/
void free_meta_keys(int nmeta,META_KEY *meta_key)
{
  int i,j;
  META_KEY *mp;
  KEY *kp,*kpt;
  
  for(i=0;i<nmeta;i++){
    kp = meta_key->keys;
    for(j=0;j<meta_key->nkeys;j++){
      kpt = kp->next;
      free(kp);
      kp = kpt;
    }
    mp = meta_key->next;
    if(i!=0) free(meta_key);
    meta_key = mp;
  }
}
/*------------------------------------------------------------------------*/
void parse_file(char *filename,int *num_meta,META_KEY *meta_root)
{
  int nch,nline,lenmeta,lenkey,ich;
  char line[MAXWORDLEN];
  WORD metakey,keyword,keyarg;
  FILE *fp;
  
  meta_now = meta_root;
  meta_root->nkeys = 0;
  meta_root->next = NULL;
  key_now = meta_root->keys = malloc(sizeof(KEY));

  nmeta = 0;

  /* open file */
  if((fp = fopen(filename,"r")) == NULL){
    fprintf(stderr,"ERROR: can't open file (%s)\n",filename);
    exit(1);
  }

  /* loop until end of file */
  nline = nch = 0;
  while( (ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if(ich == '~'){
      lenmeta=0;
      /* tilde found loop until [ is found and assign metakey */
      while( (ich = fgetc(fp)) != '[') {
        line[nch++] = (char)ich;
        /* ERROR HANDLING */
        if(ich == '\0' || ich == '~' || ich == '\n' || ich == EOF )
          syntax_error(filename,fp,nline,line,nch);
        
        if(ich != ' ' && ich != '\t'){
          metakey[lenmeta++] = (char)ich;
        }
      }
      if(lenmeta == 0) {
        printf("No meta key!\n");
        syntax_error(filename,fp,nline,line,nch);
      }      
      metakey[lenmeta++] = '\0';
      
      /* loop until we finnish this metakey */
      while( (ich = fgetc(fp)) != ']') {
        line[nch++] = (char)ich;
        
        if(nch >= MAXWORDLEN || ich==EOF){
          printf("Maxword reached or EOF reached\n");
          syntax_error(filename,fp,nline,line,nch);
        }
        
        if( ich == '\\' ){ 
          lenkey = 0;           /* get new key word */
          while( (ich = fgetc(fp)) != '{') {
            line[nch++] = (char)ich;
            if( ich == '\0' || ich == '\\' || ich == '\n'){
              printf("Illegal character in key word\n");
              syntax_error(filename,fp,nline,line,nch);
            }
            if( ich != ' ' && ich != '\t') keyword[lenkey++] = (char)ich;
          }
          keyword[lenkey] = '\0';
          line[nch++] = (char)ich;
          
          lenkey = 0;           /* get new key argument */
          while((ich= fgetc(fp)) != '}' ){
            line[nch++] = (char)ich;
            if( ich == '\0' || ich == '\\' || ich == '\n' || ich == '{'){
              printf("Illegal character in key argument\n");
              syntax_error(filename,fp,nline,line,nch);
            }
            if( ich != ' ' && ich != '\t') keyarg[lenkey++] = (char)(ich);	
          }
          keyarg[lenkey] = '\0';
          line[nch++] = (char)ich;
          addarg(keyarg,2);
          addkey(keyword);
        }
        /* if new line update nline and reset nch */
        if(ich == '\n'){nch = 0; nline++; }
      }
      /* we have complete metakey add it */
      addmeta(metakey);
    }
    if(ich == '\n'){nch = 0; nline++; }
  }
  printf("num_meta = %d\n",nmeta);
  *num_meta = nmeta;
#ifdef DEBUG
  {
    int i,j;
    
    meta_now = meta_root;
    for(i=0;i<nmeta;i++,meta_now=meta_now->next){
      printf("%d %s\n",i,meta_now->meta);
      key_now = meta_now->keys;
      for(j=0;j<meta_now->nkeys;j++,key_now=key_now->next){
        printf("%d \\%s{%s}\n",j,key_now->name,key_now->arg.sval);      
      }
    }
  }
  
#endif
}
/*------------------------------------------------------------------------*/

