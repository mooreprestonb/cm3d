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


/* math routines that might be useful */

double gethinv9(double hmat[], double hmati[])
{
  double rvol;

  /* gets inverse, hmati, of the iperd dimensional matrix hmat */
  rvol = 1./(hmat[0]*(hmat[4]*hmat[8] - hmat[5]*hmat[7]) + 
	     hmat[1]*(hmat[5]*hmat[6] - hmat[3]*hmat[8]) + 
	     hmat[2]*(hmat[3]*hmat[7] - hmat[4]*hmat[6]));
  
  hmati[0] = (hmat[4]*hmat[8] - hmat[5]*hmat[7])*rvol;
  hmati[1] = (hmat[2]*hmat[7] - hmat[1]*hmat[8])*rvol;
  hmati[2] = (hmat[1]*hmat[5] - hmat[2]*hmat[4])*rvol;
  hmati[3] = (hmat[5]*hmat[6] - hmat[3]*hmat[8])*rvol;
  hmati[4] = (hmat[0]*hmat[8] - hmat[2]*hmat[6])*rvol;
  hmati[5] = (hmat[2]*hmat[3] - hmat[0]*hmat[5])*rvol;
  hmati[6] = (hmat[3]*hmat[7] - hmat[4]*hmat[6])*rvol;
  hmati[7] = (hmat[1]*hmat[6] - hmat[0]*hmat[7])*rvol;
  hmati[8] = (hmat[0]*hmat[4] - hmat[1]*hmat[3])*rvol;

  return (1./rvol);
}

/*-------------------------------------------------------------------------*/
double get_deth(double *hmat)
{
  /* gets det of hmat */
  
  return (hmat[0]*(hmat[4]*hmat[8] - hmat[5]*hmat[7]) + 
	  hmat[1]*(hmat[5]*hmat[6] - hmat[3]*hmat[8]) + 
	  hmat[2]*(hmat[3]*hmat[7] - hmat[4]*hmat[6]));
} 
/*-------------------------------------------------------------------------*/
void sym33(double mat[9])
{
  double tmp;

  tmp = .5*(mat[1] + mat[3]);
  mat[1] = mat[3] = tmp;
  tmp = .5*(mat[2] + mat[6]);
  mat[2] = mat[6] = tmp;
  tmp = .5*(mat[5] + mat[7]); 
  mat[5] = mat[7] = tmp;
}
/*-------------------------------------------------------------------------*/
#ifndef DIM
#define DIM 3
#endif

void v1fv3(int n,double *all, double *x, double *y, double *z)
{
  int i;
  for(i=0;i<n;++i){
    all[i*DIM] = x[i]; all[i*DIM+1] = y[i]; all[i*DIM+2] = z[i];
  }
}

/*-------------------------------------------------------------------------*/
void v3fv1(int n,double *all, double *x, double *y, double *z)
{
  int i;
  for(i=0;i<n;++i){
    x[i] = all[i*DIM];
    y[i] = all[i*DIM+1];
    z[i] = all[i*DIM+2];
  }
}
/*-------------------------------------------------------------------------*/

double ddot1a(int n,double *a,int astep,double *b,int bstep)
{
  int i,j,k;
  double ddot1;

  ddot1 = 0.;
  for(i=0,j=0,k=0; i<n ;i++,j+=astep,k+=bstep){
    ddot1 += a[j]*b[k];
  }

  return ddot1;
}

/*-------------------------------------------------------------------------*/
double ddot1(int n,double *a,double *b)
{
  int i;
  double ddot1;

  ddot1 = 0.;
  for(i=0; i<n ;i++) ddot1 += a[i]*b[i];

  return ddot1;
}

/*-------------------------------------------------------------------------*/

double dsum1a(int n,double *a,int astep)
{
  int i,j;
  double dsum1;

  dsum1 = 0.;
  for(i=0,j=0;i<n;i++,j+=astep){
    dsum1 += a[j];
  }

  return dsum1;
}
/*-------------------------------------------------------------------------*/

double dsum1(int n,double *a)
{
  int i;
  double dsum1;

  dsum1 = 0.;
  for(i=0; i<n ;i++) dsum1 += a[i];

  return dsum1;
}
/*-------------------------------------------------------------------------*/
void matmulv(int n,double a[],double b[],double c[]) /* mat a*b=c */
{
  int n1,m2,k;
  double sum;
  
  for(n1=0;n1<n;n1++){
    for(m2=0;m2<n;m2++){
      sum = 0.;
      for(k=0;k<n;k++){
	sum += a[n1*n+k]*b[k*n+m2];
      }
      c[n1+m2*n] = sum;
    }
  }
}
/*-------------------------------------------------------------------------*/
void transmatv(int n,double a[],double b[]) /* transpost a and put into b */
{
  int i,j;

  for(i=0;i<n;i++){
    b[i*n+i] = a[i*n+i];
    for(j=i+1;j<n;j++){
      b[j*n+i] = a[i*n+j];
      b[i*n+j] = a[j*n+i];
    }
  }
}
