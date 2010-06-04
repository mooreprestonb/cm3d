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

/* routines to output strings to the screen */

#include "md.h"

FILE *errnow=NULL;
FILE *outnow=NULL;

/*---------------------------------------------------------*/
void set_stderrnow(char *s)
{
  FILE *errtemp=errnow;
  if(errnow==NULL) errnow = stderr;
  
  fprintf(errtemp,"Setting stderr to %s\n",s);
  if((errnow = fopen(s,"a+"))==NULL){
    fprintf(errtemp,"ERROR: can't open %s to write our stderr\n",s);
    fprintf(errtemp,"\texiting system\n");
    exit(1);
  }
  if(errtemp!=errnow && errtemp!=stderr) fclose(errtemp);
}
/*---------------------------------------------------------*/
void set_stdoutnow(char *s)
{
  LINE line;  
  if(errnow==NULL) errnow = stderr;
  if(outnow==NULL) outnow = stdout;
  fprintf(errnow,"Setting stdout to %s\n",s);
  if(outnow != stdout) fclose(outnow);
  if((outnow = fopen(s,"a+"))==NULL){
    sprintf(line,"can't open %s to write our stdout",s);
    md_error(line);
    exit(1);
  }  
}
/*---------------------------------------------------------*/
void usage(char *s)
{
  if(errnow==NULL) errnow = stderr;
  fprintf(errnow,"Usage: %s [flags] <infile> [outfile]\n",s);
  exit(1);
}
/*--------------------------------------------------------------*/
void md_warning(char *s)
{
  if(errnow==NULL) errnow = stderr;
  fprintf(errnow,"WARNING: %s\n",s);
}
/*-------------------------------------------------------------*/
void md_error(char *s)
{
  if(errnow==NULL) errnow = stderr;
  fprintf(errnow,"ERROR: %s\n",s);
  fprintf(errnow,"\texiting system\n");
#ifdef PARA
  MPI_Abort(MPI_COMM_WORLD,1);
#endif
  exit(1);
}
/*-------------------------------------------------------------*/
void md_stdout(char *s)
{
  if(outnow==NULL) outnow = stdout;
  fprintf(outnow,"%s\n",s); fflush(stdout); fflush(stderr);
}
/*-------------------------------------------------------------*/
