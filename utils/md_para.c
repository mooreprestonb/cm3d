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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#ifndef LONG_MAX
#define LONG_MAX        2147483647      /* max value of a "long int" */
#endif

#ifndef INT_MAX
#define INT_MAX        INT_MAX     /* max value of a "long int" */
#endif

/*-------------------------------------------------------------*/

void decomp1d(int n,int size,int rank,long *s,long *e)
{
  long nlocal, deficit;

  nlocal= n/size;     /*nice local size */
  *s=rank * nlocal;
  deficit=n%size;     /* remainder */
  *s += ((rank < deficit) ? rank : deficit); /* add deficit to s if nec*/
  if (rank < deficit) nlocal++;
  *e = *s + nlocal;            /* calculate end point */
  if (*e > n || rank == size-1) *e = n;
}

/*-------------------------------------------------------------*/
void decomp1d_old(int n,int size,int rank,int *s,int *e)
{
  int nlocal, deficit;

  nlocal= n/size;     /*nice local size */
  *s=rank * nlocal + 1;
  deficit=n%size;     /* remainder */
  *s= *s + ((rank < deficit) ? rank : deficit); /* add deficit to s if nec*/
  if (rank < deficit) nlocal++;
  *e= *s + nlocal - 1;
  if (*e > n || rank == size-1) *e = n;
  --*s;  /* fix s so that it starts at zero */
}
/*-------------------------------------------------------------*/
/* decompose top i loop of i=0,n-1 and j=i+1,n-1 about evenly 
   (some j factors off) */

/*-------------------------------------------------------------*/
void decomp_trig(int n,int size,int rank,long *ibegin,long *iend)
{ 
  unsigned int nlocal,nl;
  int s,e;

  if(n==0){
    *ibegin = 0;
    *iend = 0;
  }else if(n>(INT_MAX/n)){
    /* npairs to big for integer so we spread i evenly */
#ifdef WARNTRIG
    fprintf(stderr,"Warning: npairs (%d) is too big for integer 2\n",n);
    fprintf(stderr,"\t try using different neighbor lists\n");
    fprintf(stderr,"\t %d %d\n",n,INT_MAX);
#endif
    decomp1d(n,size,rank,ibegin,iend);
  } else {
    nl = (n%2) ? n*((n-1)/2) : (n/2)*(n-1);
    nlocal = nl/size;
    rank = size-1-rank;
    s = rank * nlocal;
    e = s + nlocal;
    if(rank==size-1) *ibegin=0;
    else *ibegin = n-(int)(sqrt((double)(2*e+.25))+.5);
    if(rank==0) *iend = n;
    else  *iend = n-(int)(sqrt((double)(2*s+.25))+.5);
  }
}

void decomp_trig_old(int n,int size,int rank,long *ibegin,long *iend)
{ 
  unsigned long i;
  long s,e;

  if(n==0){
    *ibegin = 0;
    *iend = 0;
  } else if(2.>((((LONG_MAX/n)/n)/n)/n)){
    /* npairs to big for integer so we spread i evenly */
#ifdef WARNTRIG
    fprintf(stderr,"Warning: npairs (%d) is too big for integer\n",n);
    fprintf(stderr,"\t try using different neighbor lists\n");
    exit(1); 
#endif
    decomp1d(n,size,rank,ibegin,iend);
  } else {    
    i = n*(n-1)/2;
    decomp1d(i,size,rank,&s,&e);
    i = 1+4*n*(n-1);
    *ibegin = (int)((-1.+2.*n-sqrt((double)(i-8*s)))/2.);
    *iend = (int)((-1.+2.*n-sqrt((double)(i-8*e)))/2.);
  }
}


#ifdef MPI_NONEXISTANT
void MPI_Init(int *argc,char ***argc){ printf("MPI_NONEXISTANT\n");}
void MPI_Comm_rank(int comm, int *rank) {*rank = 0;}
void MPI_Comm_size(int comm, int *size) {*size = 1;}
void MPI_Finalize(void){printf("MPI_Finalize\n");}
void MPI_Bcast(void *junk,int n,int,int,int comm){}
void MPI_Allreduce(void *send,void *rec,int n,int type,int op,int comm){
  int i;
  for(i=0;i<n;i++){
    rec[i]=send[i];
  }
}
#endif

