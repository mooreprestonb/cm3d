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

/* subroutines to return memory of matrix or tensor3 */

#include <stdlib.h>
#include <stdio.h>

#ifndef MAX
#define MAX(A,B) ((A>B)?(A):(B))
#endif

#ifndef MIN
#define MIN(A,B) ((A<B)?(A):(B))
#endif
double *dvector(int nl,int nh);
void free_dvector(double *v,int nl,int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
double **drematrix(double **,int ,int ,int ,int ,int ,int ,int ,int );
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
double ***d3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh);
void free_d3tensor(double ***,int,int,int,int,int,int);

/*---------------------------------------------------------------------*/
/* allocate an double vector with subscipt range v[nl...nh] */
double *dvector(int nl,int nh)
{
  double *v; 
  v = (double *)malloc((nh-nl+1)*sizeof(double));
  return v-nl;
}
/*---------------------------------------------------------------------*/
void free_dvector(double *v,int nl,int nh){ free(v+nl);}

/*---------------------------------------------------------------------*/
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow,ncol;
  double **m;

  nrow = nrh-nrl+1;  ncol = nch-ncl+1;

  if(nrow*ncol==0) return (double **)NULL;
  /* allocate pointers to rows */

  m = (double **)malloc(nrow*sizeof(double*));
  if(!m){
    fprintf(stderr,"ERROR: can't allocate memory for matrix (%ld bytes)\n",
	    (long) nrow*sizeof(double *));
    exit(1);
  }
  m -= nrl;

  /* allocate rows */
  m[nrl] = (double *)calloc(nrow*ncol,sizeof(double));
  if(!m[nrl]){
    fprintf(stderr,
	    "ERROR: can't allocate memory for arrays for matrix (%ld bytes)\n",
	    (long)nrow*ncol*sizeof(double));
    exit(1);
  }
  m[nrl] -= ncl;

  /* set pointers to rows */
  for(i=nrl+1;i<=nrh;i++) m[i] = m[i-1]+ncol;

  return m;
}

/*----------------------------------------------------------------------*/
/* 
   reallocate a double matrix with subscript range m[nrlo..nrho][nclo..ncho] 
   to a double matrix with subscript range m[nrln..nrhn][ncln..nlhn]*/

double **drematrix(double **mat,int nrlo,int nrho,int nclo,int ncho,
		   int nrln,int nrhn,int ncln,int nchn)
{
  int i,j;
  double **m;

  if(nrlo != nrln || nclo != ncln){
    fprintf(stderr,"ERROR: in realocating matrix ... indicies don't match\n");
    exit(1);
  }
  /* allocate new space */
  m = dmatrix(nrln,nrhn,ncln,nchn);

  /* copy memory */
  for(i=nrlo;i<=MIN(nrho,nrhn);i++){
    for(j=nclo;j<=MIN(ncho,nchn);j++){
      m[i][j] = mat[i][j];
    }
  }
  free_dmatrix(mat,nrlo,ncho,nclo,ncho);
  return m;
}

/*----------------------------------------------------------------------*/
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
  free((char *)(m[nrl]+ncl)); free((char *)(m+nrl));
}

/*----------------------------------------------------------------------*/

double ***d3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
  int i,j,nrow,ncol,ndep;
  double ***t;

  nrow = nrh-nrl+1; ncol = nch-ncl+1; ndep = ndh-ndl+1;

  /* allocate top pointer */

  t = (double ***)malloc(nrow*sizeof(double**));
  if(!t){
    fprintf(stderr,"ERROR: (1) can't allocate memory for tensor\n");
    exit(1);
  }
  t -= nrl;

  /* allocate matrix pointers */
  t[nrl] = (double **)malloc(nrow*ncol*sizeof(double*));
  if(!t[nrl]){
    fprintf(stderr,"ERROR: (2) can't allocate memory for arrays in tensor\n");
    exit(1);
  }
  t[nrl] -= ncl;

  /* allocate final pointers and set all pointers */
  t[nrl][ncl] = (double *)calloc(nrow*ncol*ndep,sizeof(double));
  if(!t[nrl][ncl]){
    fprintf(stderr,"ERROR: (3) can't allocate memory for rows in tensor\n");
    exit(1);
  }
  t[nrl][ncl] -= ndl;


  /* set pointers to rows */
  for(j=ncl+1;j<=nch;j++) t[nrl][j] = t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++){
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j] = t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */

  return t;
}
/*----------------------------------------------------------------------*/

void free_d3tensor(double ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
  free((char *)(t[nrl][ncl]+ndl)); 
  free((char *)(t[nrl]+ncl));
  free((char *)(t+nrl));
}
/*----------------------------------------------------------------------*/
