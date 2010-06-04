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
#include <math.h>
#include <mpi.h>
#include <assert.h>

#include "parautil.h"
#include "mpi3dfft.h"

#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265359
#endif


typedef struct 
{
   double re;
   double im;
} dcomplex;


int nProc,myRank;

static void init(dcomplex *u,int nx,int ny,int nz,
   int rx,int ry,int rz)
/****************************************************************/
/* Initializes u[] (nx*ny*nz distributed on 1st index)
   to (1+i)*delta(r-(rx,ry,rz))
   */
{
   int i,dnx,i0nx;
   decomp1d(nx,nProc,myRank,&i0nx,&dnx);
   for ( i=0; i<dnx*ny*nz; i++ ) {
      u[i].re=u[i].im=0;
   }
   if ( i0nx<=rx && rx<i0nx+dnx ) {
      int ir=((rx-i0nx)*ny+ry)*nz+rz;
      u[ir].re=u[ir].im=0.5;
   }
}


static void check(dcomplex *ut,int nx,int ny,int nz,
   int rx,int ry,int rz,int isign)
/*****************************************************************/
/* Assumes u[] was loaded with (1+i)*delta(R-(rx,ry,rz)) and checks
   that ut[] contains transform of this, which is
   (1+i)*exp(isign*i*(gx*rx+gy*ry+gz*rz))
   */
{
   double eps=1e-12;
   int i,j,k,dnx,i0nx;

   decomp1d(nx,nProc,myRank,&i0nx,&dnx);
   for ( i=0; i<dnx; i++ ) {
      int ii=i0nx+i;
      for ( j=0; j<ny; j++ ) {
         for ( k=0; k<nz; k++ ) {
            double a,b,gr = 2*PI*(((double)ii*rx)/nx+((double)j*ry)/ny+((double)k*rz)/nz);
            int ig=(i*ny+j)*nz+k;
            a = cos(gr)-isign*sin(gr);
            b = cos(gr)+isign*sin(gr);
            a -= ut[ig].re;
            b -= ut[ig].im;
            assert(fabs(a)<=eps && fabs(b)<=eps);
         }
      }
   }
}



int main(int argc,char *argv[])
/*****************************************************************/
{
   dcomplex *uR,*uG;

   int n1 = 4;
   int n2 = 8;
   int n3 = 16;
   int nR,distNr;
   int nG,distNg;

   int dn1,dn2,i10,i20;
   int i,j,k;
   
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nProc);
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   decomp1d(n1,nProc,myRank,&i10,&dn1);
   nR = n1*n2*n3;
   distNr = dn1*n2*n3;
   
   decomp1d(n2,nProc,myRank,&i20,&dn2);
   nG = n1*n2*n3;
   distNg = dn2*n1*n3;

   uR = (dcomplex *) malloc(distNr*sizeof(*uR));
   assert(uR);
   uG = (dcomplex *) malloc(distNg*sizeof(*uG));
   assert(uG);

/* check transform R to G */
   for ( i=0; i<n1; i++ ) {
      for ( j=0; j<n2; j++ ) {
         for ( k=0; k<n3; k++ ) {
            init(uR,n1,n2,n3,i,j,k);
            mpiZ3dFFT(uR,n1,n2,n3,uG,-1,2.0,MPI_COMM_WORLD);
            check(uG,n2,n1,n3,j,i,k,-1);
         }
      }
   }
   
/* check transform G to R */
   for ( j=0; j<n2; j++ ) {
      for ( i=0; i<n1; i++ ) {
         for ( k=0; k<n3; k++ ) {
            init(uG,n2,n1,n3,j,i,k);
            mpiZ3dFFT(uG,n2,n1,n3,uR, 1,2.0,MPI_COMM_WORLD);
            check(uR,n1,n2,n3,i,j,k, 1);
         }
      }
   }

   free(uR);
   free(uG);
   
   MPI_Finalize();

   return 0;
}



