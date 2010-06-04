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

/* from fp get the meta key and all it key arguments */

#include "md.h"

/*--------------------------------------------------------------*/
void get_meta_key(FILE *fp,char *filename,char *line,int *nch,int *nline,
	     char *metakey,KEY *key_root,int *num_keys)
{
  int lenmeta=0,ich;

  /* tilde found loop until [ is found and assign metakey */
  while( (ich = fgetc(fp)) != '[') {
    line[(*nch)++] = (char)ich;
    /* ERROR HANDLING */
    if(ich == '\0' || ich == '~' || ich == '\n' || ich == EOF )
      syntax_error(filename,fp,*nline,line,*nch);

    if(ich != ' ' && ich != '\t'){
      metakey[lenmeta++] = (char)ich;
    }
  }
  if(lenmeta == 0) syntax_error(filename,fp,*nline,line,*nch);
  metakey[lenmeta++] = '\0';

  get_keys(fp,filename,line,nline,nch,num_keys,key_root);
}
/*---------------------------------------------------------*/ 

void get_keys(FILE *fp,char *filename,char *line,int *nline,int *nch,
	 int *num_keys,KEY *key_root)
{
  int ich;
  int lenkey,lenarg;
  KEY *key_now;
  
  *num_keys = 0;
  lenkey = lenarg = 0;
  key_now = key_root;
  key_now->next = NULL;
  strcpy(key_now->keyword,"");
  strcpy(key_now->keyarg,"");

  while( (ich = fgetc(fp)) != ']') {
    line[(*nch)++] = (char)ich;
    
    if((*nch) >= MAXLINELEN || ich==EOF)
      syntax_error(filename,fp,*nline,line,*nch);
    
    if( ich == '\\' ){
      /* get new key word */
      (*num_keys)++;
      while( (ich = fgetc(fp)) != '{') {
	line[(*nch)++] = (char)ich;
	if( ich == '\0' || ich == '[' || ich == '\\' || 
	   ich == '\n'  || ich == '~')
	  syntax_error(filename,fp,*nline,line,*nch);
	
	if( ich != ' ' && ich != '\t')
	  key_now->keyword[lenkey++] = (char)ich;
      }
      key_now->keyword[lenkey] = '\000';
      line[(*nch)++] = (char)ich;
      /* get new key argument */
      while((ich= fgetc(fp)) != '}' ){
	line[(*nch)++] = (char)ich;
	if(ich == '\0' || ich == '[' || ich == '\\' || 
	   ich == '\n' || ich == '{' || ich == '~')
	  syntax_error(filename,fp,*nline,line,*nch);
	
	if( ich != ' ' && ich != '\t')
	  key_now->keyarg[lenarg++] = (char)(ich);	
      }
      if(lenarg==0) syntax_error(filename,fp,*nline,line,*nch);
      key_now->keyarg[lenarg] = '\000';
      line[(*nch)++] = (char)ich;

      key_now->next = (KEY *)malloc(sizeof(KEY));
      key_now = key_now->next;  lenkey = 0;
      lenarg = 0;      key_now->next = NULL;

    }
    /* if new line update nline and reset nch */
    if(ich == '\n'){*nch = 0;(*nline)++;}
  }
}
/*------------------------------------------------------------------------*/
void free_key(int num_keys,KEY *key_now)
{
  KEY *key_next;

  /* set over key_root */
  key_now = key_now->next;

  /* loop over allocated keys */
  while(num_keys-->0){
    key_next = key_now->next; free(key_now); key_now = key_next;
  }
}
/*------------------------------------------------------------------------*/

void syntax_error(const char *infile,FILE *fp,int nline,char *line,int nch)
{
  int ich;

  fprintf(stderr,"ERROR: Syntax error in file %s at line %d\n",infile,nline);
  for(ich=0;ich<nch;ich++) fprintf(stderr,"%c",line[ich]);
  if(line[--nch] != '\n')  fprintf(stderr,"%c",'\n');
  for(ich=0;ich<nch;ich++) fprintf(stderr,"%c",'_');
  fprintf(stderr,"%c\n\n",'^');
  fclose(fp);
  exit(1);
}
/*---------------------------------------------------------------------*/
int get_dict_num(int num_dict,WORD dict[],char *name)
{
  int i;
  
  for(i=0;i<num_dict;i++){
    if(!(strcasecmp(dict[i],name))){
      return i;
    }
  }
  return -1;
}
/*-----------------------------------------------------------------------*/ 
void print_dict(WORD keyword,int num,WORD dict[])
{
  int i,j;
  fprintf(stderr,"ERROR: keyword \"%s\" not in the dictionary\n",keyword);
  for(j=0;j<num/2;j++){
    i = j*2;
    if(dict[i][0] != '\0') fprintf(stdout,"dict[%d] = %s\t",i,dict[i]);
    if(dict[i+1][0] != '\0')fprintf(stdout,"\tdict[%d] = %s",i+1,dict[i+1]);
    fprintf(stdout,"\n");
  }
  if(num%2 !=0){
    if(dict[num-1][0] != '\0') 
      fprintf(stdout,"dict[%d] = %s\t",num-1,dict[num-1]);
    fprintf(stdout,"\n");
  }
  exit(1);
}
/*-----------------------------------------------------------------------*/ 
void err_found(WORD dict,int num_found,char *filename,int nline)
{
  fprintf(stderr,"ERROR: keyword %s found %d times in file %s at line %d\n",
	  dict,num_found,filename,nline);
  exit(1);
}  
/*-----------------------------------------------------------------------*/ 
void write_key(KEY *key){fprintf(stdout,"%s %s\n",key->keyword,key->keyarg);}
/*-----------------------------------------------------------------------*/ 
void get_sim_keys(const char *command,const char *infile,
   int *num_keys,KEY *key_root)
{
  char backslash,nil,newline,lbrace,rbrace,space,tab;
  char line[MAXLINELEN],string[MAXLINELEN];
  int nline,nch,ich;
  int lenkey,lenarg;
  FILE *fp;
  KEY *key_now;
  
#ifdef DMALLOC
  dmalloc_verify(NULL);
#endif
  key_now = key_root; lenkey = 0; lenarg = 0;
  
  backslash = '\\'; nil = '\0'; newline = '\n'; lbrace = '{';
  rbrace = '}'; space = ' '; tab = '\t';
  
  nline =  1; nch = 0;

#ifdef DMALLOC
  dmalloc_verify(NULL);
#endif
  
  fp = fopen(infile,"r");
  if (fp == NULL){
    sprintf(string,"%s can't open file \"%s\"",command,infile);
    md_error(string);
  }
  
  while( (ich = fgetc(fp)) != EOF) {
    line[nch++] = (char)ich;
    if( ich == backslash){
      (*num_keys)++;      
      while( (ich = fgetc(fp)) != lbrace){
        line[nch++] = (char)ich;
        if(ich == EOF || ich == backslash || ich == newline)
          syntax_error(infile,fp,nline,line,nch);
        if(ich != space && ich != tab)
          key_now->keyword[lenkey++] = (char)ich;
      }
      key_now->keyword[lenkey] = nil;
      line[nch++] = (char)ich;
      while( (ich = fgetc(fp)) != rbrace){
        line[nch++] = (char)ich;
        if(ich == EOF || ich == backslash || ich == newline || ich == lbrace)
          syntax_error(infile,fp,nline,line,nch);
        if(ich != space && ich != tab)
          key_now->keyarg[lenarg++] = (char)ich; 
      }
      key_now->keyarg[lenarg] = nil;
      line[nch++] = (char)ich;
      key_now->next = (KEY *)malloc(sizeof(KEY));
      key_now = key_now->next;
      lenkey = 0;
      lenarg = 0;
    }
    if(ich == '\n'){nch = 0;nline++;}
  }
  fclose(fp);
}
/*---------------------------------------------------------------------*/ 
int getnxtint (char *charinfo)
{
  static char *charinfoold=NULL;
  static int iloc;
  static int digit;
  int n;

  if(charinfoold != charinfo){
    charinfoold=charinfo;
    iloc = 0;
    digit = 0;
  }

  if(charinfo[iloc] == '\0') return digit;
  if (strcmp(charinfoold,charinfo)==0){
    while(!isdigit(charinfo[iloc]) && (charinfo[iloc] != '-')) {
      if(charinfo[iloc] == '\0') return digit;
      iloc++;
    }
  } else {
    iloc = 0;
    charinfoold = charinfo;
  }
  sscanf(&charinfo[iloc],"%d%n",&digit,&n);

  iloc += n;
  return digit;
}




