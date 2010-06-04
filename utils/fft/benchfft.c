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
#include <mpi.h>
#include <assert.h>

#include "parautil.h"
#include "mpi3dfft.h"
#include "prf.h"

typedef struct {double re;double im;} dcomplex;

int nProc,myRank;

static void loadZMRandom(dcomplex *u,int n1,int n2,int n3);
static void moments(dcomplex *u,int n1,int n2,int n3,
   dcomplex *um,double *u2);

static void usage(void)
/*****************************************************************/
{
   printf("Usage:\n\n");
   printf("benchfft\n");
   printf("    does a single pair of 64*64*64 transforms\n");
   printf("benchfft n1 n2 n3\n");
   printf("    does a single pair of n1*n2*n3 transforms\n");
   printf("benchfft n1 n2 n3 niter\n");
   printf("    does niter pairs of n1*n2*n3 transforms\n");
}
   

int main(int argc,char *argv[])
/*****************************************************************/
{
   static Profiler prfMain;
   static Profiler prfMainOther;
   static Profiler prfFFTTotal;
   
   dcomplex *uR,*uG;
   dcomplex urm,ugm;
   double ssqr1,ssqg,ssqr2;
   int n1,n2,n3,nR,distNr,nG,distNg;
   int dn1,dn2,i10,i20;
   int i,ip,nIter;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nProc);
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   switch ( argc ) {
      case 1:
         n1=n2=n3=64;
         nIter=1;
         break;
      case 4:
         n1 = atoi(argv[1]);
         n2 = atoi(argv[2]);
         n3 = atoi(argv[3]);
         nIter=1;
         break;
      case 5:
         n1 = atoi(argv[1]);
         n2 = atoi(argv[2]);
         n3 = atoi(argv[3]);
         nIter = atoi(argv[4]);
         break;
      default:
         usage();
         return(0);
   }

   n1 = selectFFTGrid(n1);
   n2 = selectFFTGrid(n2);
   n3 = selectFFTGrid(n3);

   if ( myRank==0 ) {
      printf("\n");
      printf("Selected real space FFT Grid is %d*%d*%d\n",n1,n2,n3);
      printf(" transposed g space FFT Grid is %d*%d*%d\n\n",n2,n1,n3);

      printf("Distribution:\n");
      printf("--------------------------------------\n");
      printf("Rank   Real Space             Transform Space\n");
      for ( ip=0; ip<nProc; ip++ ) {
         int off1,len1,off2,len2;
         decomp1d(n1,nProc,ip,&off1,&len1);
         decomp1d(n1,nProc,ip,&off2,&len2);
         printf("%2d     (%3d:%3d,%3d,%3d)      (%3d:%3d,%3d,%3d)\n",ip,
            off1,off1+len1-1,n2,n3,off2,off2+len2-1,n1,n3);
      }
      printf("\n");
   }
   
   NewProfiler("Main (total)",&prfMain);
   NewProfiler("Main (other)",&prfMainOther);
   NewProfiler("FFT (total)",&prfFFTTotal);
   ProfileStart(prfMain);
   ProfileStart(prfMainOther);

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

   /* fill uR with zero mean random numbers */
   loadZMRandom(uR,n1,n2,n3);
   moments(uR,n1,n2,n3,&urm,&ssqr1);
   
   /* transform to G-space then back to R-space */
   ProfileStart(prfFFTTotal);
   for ( i=0; i<nIter; i++ ) {
      mpiZ3dFFT(uR,n1,n2,n3,uG,-1,1.0/nR,MPI_COMM_WORLD);
      mpiZ3dFFT(uG,n2,n1,n3,uR, 1,1.0,   MPI_COMM_WORLD);
   }
   ProfileStopI(prfFFTTotal);

   /* check moments */
   moments(uG,n2,n1,n3,&ugm,&ssqg);
   ssqg *= n1*n2*n3;
   moments(uR,n1,n2,n3,&urm,&ssqr2);
   
   
   ProfileStop(prfMainOther);
   ProfileStopI(prfMain);

   if ( myRank==0 ) {
      printf("Before transform,     ssqr1 = %f\n",ssqr1);
      printf("After transform,       ssqg = %f\n",ssqg);
      printf("Back again transform, ssqr1 = %f\n",ssqr2);
   }

   if ( myRank==0 ) {
      DumpProfileResults(stdout);
   }
   
   free(uR);
   free(uG);
   
   MPI_Finalize();

   return 0;
}


static void loadZMRandom(dcomplex *u,int n1,int n2,int n3)
/***********************************************************/
/* fill u with zero mean random numbers */
{
   dcomplex sum,mySum;
   int i,dn1,i0n1;

   decomp1d(n1,nProc,myRank,&i0n1,&dn1);

   mySum.re=0;
   mySum.im=0;
   for ( i=0; i<dn1*n2*n3; i++ ) {
      dcomplex u0;
      u0.re = rand();
      u0.im = rand();
      mySum.re += (u[i].re = (u0.re/RAND_MAX-0.5));
      mySum.im += (u[i].im = (u0.im/RAND_MAX-0.5));
   }
   MPI_Allreduce(&mySum,&sum,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   sum.re/=n1*n2*n3;
   sum.im/=n1*n2*n3;

   for ( i=0; i<dn1*n2*n3; i++ ) {
      u[i].re-=sum.re;
      u[i].im-=sum.im;
   }
}

static void moments(dcomplex *u,int n1,int n2,int n3,
   dcomplex *um,double *u2)
/***********************************************************/
{
   dcomplex sum,mySum;
   double s2,mys2;
   int i,dn1,i0n1;

   decomp1d(n1,nProc,myRank,&i0n1,&dn1);

   mySum.re=0;
   mySum.im=0;
   
   for ( i=0; i<dn1*n2*n3; i++ ) {
      mySum.re += u[i].re;
      mySum.im += u[i].im;
   }
   MPI_Allreduce(&mySum,&sum,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   sum.re /= n1*n2*n3;
   sum.im /= n1*n2*n3;

   mys2=0;
   for ( i=0; i<dn1*n2*n3; i++ ) {
      double a,b;
      a = u[i].re;
      b = u[i].im;
      mys2 += a*a+b*b;
   }
   
   MPI_Allreduce(&mys2,&s2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   *um = sum;
   *u2 = s2;
}

