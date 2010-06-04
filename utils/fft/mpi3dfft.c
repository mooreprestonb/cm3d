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
#include <mpi.h>
#include <assert.h>
#include "prf.h"

/* define this here so that we don't force a particular name
   for this type on anyone. *I* may like "dcomplex" but maybe
   their are some really strange people who like "zomplex"
   or some nonsense like that
   */
typedef struct {
   double re;
   double im;
} dcomplex;


static void zTranspose12(dcomplex *ut,dcomplex *u,int n1,int n2,int n3,
   MPI_Datatype mpitype,MPI_Comm comm);



void mpiZ3dFFT(void *u,int n1,int n2,int n3,void *v,int isign,
   double scale,MPI_Comm comm)
/*********************************************************************/
/* This routine takes an n1*n2*n3 3-d double complex array that is
   distributed across processors on its first (n1) index (with MPI
   communicator comm), and performs a 3d-FFT.

   The calculated transform is n(r) = scale * sum(n(g)*exp(isign*I*g*r),g)
   (where r and g have no real meaning here -- choose your own convention!),
   so isign sets the sign of I in the exponential, and scale is the
   multiplicative factor.

   The transformed data is placed in v (which may be the same as u), but it
   is an n2*n1*n3 array, also distributed on the first (n2, now) index.

   ADDITIONAL ASSUMPTIONS:

   - Distribution is according to my decomp1d routine.
   - Each node has exactly dist(n1)*n2*n3 elements, with no extra "wings"
     or holes.
   - (n1,n2,n3) satisfy the length requirements of the underlying
     1d-FFT routines
   - The data are arranged as array of struct { double re,im; } in
     standard C order.
   */
{
   static Profiler prfFFT1;
   static Profiler prfFFT2;
   static Profiler prfFFT3;
   static Profiler prfFFTAll;
   static MPI_Datatype MPIDComplex;
   static int isInit=0;
   dcomplex *pu=(dcomplex *)u;
   dcomplex *pv=(dcomplex *)v;
   dcomplex *pvbuf=0;
   MPI_Comm myComm;
   int nProc;
   int myRank;
   int i0n1,dn1,i0n2,dn2,distNtot;
   int i,islice;
   
   NewProfiler("FFT on 1st",&prfFFT1);
   NewProfiler("FFT on 2nd",&prfFFT2);
   NewProfiler("FFT on 3rd",&prfFFT3);
   NewProfiler("FFT rest",  &prfFFTAll);

   ProfileStart(prfFFTAll);
   
   /* Get our own unique communication space to avoid our messages
      colliding with any belonging to the caller */
   MPI_Comm_dup(comm,&myComm);
   MPI_Comm_size(myComm,&nProc);
   MPI_Comm_rank(myComm,&myRank);

   /* Make a type for complex numbers */
   if ( !isInit ) {
      MPI_Type_contiguous(2,MPI_DOUBLE,&MPIDComplex);
      MPI_Type_commit(&MPIDComplex);
   }

   decomp1d(n1,nProc,myRank,&i0n1,&dn1);
   decomp1d(n2,nProc,myRank,&i0n2,&dn2);
   distNtot = dn1*n2*n3;

   pvbuf = (dcomplex *) malloc((size_t)(distNtot*sizeof(*pvbuf)));
   assert(pvbuf);

   for ( i=0; i<distNtot; i++ ) {
      pvbuf[i]=pu[i];
   }

   /* Do all transforms along 3rd index */
   ProfileStart(prfFFT3);
   for ( islice=0; islice<dn1; islice++ ) {
      m1dFFT(isign,&pvbuf[islice*n2*n3],n2,n3,n3,1,1.0);
   }
   ProfileStop(prfFFT3);

   /* Do all transforms along 2nd index */
   ProfileStart(prfFFT2);
   for ( islice=0; islice<dn1; islice++ ) {
      m1dFFT(isign,&pvbuf[islice*n2*n3],n3,1,n2,n3,1.0);
   }
   ProfileStop(prfFFT2);

   /* Global transpose first two dimensions */
   zTranspose12(pv,pvbuf,n1,n2,n3,MPIDComplex,myComm);

   /* Do all transforms along (new) second index */
   ProfileStart(prfFFT1);
   for ( islice=0; islice<dn2; islice++ ) {
      m1dFFT(isign,&pv[islice*n1*n3],n3,1,n1,n3,scale);
   }
   ProfileStop(prfFFT1);
   
   /* Free everything up (except for the MPIDComplex type) */
   free(pvbuf);
   MPI_Comm_free(&myComm);

   ProfileStop(prfFFTAll);
}


static void zTranspose12(dcomplex *ut,dcomplex *u,int n1,int n2,int n3,
   MPI_Datatype mpitype,MPI_Comm comm)
/******************************************************************/
/* Semigeneric function that transposes the first two indicies
   of an n1*n2*n3 array that is distributed on the 1st index.
   The transposed array is n2*n1*n3, also distributed on (its) first
   index.
   */
{
   static Profiler prfTransGlob;
   static Profiler prfTransLoc;
   int nProc,myRank;
   int i,j,k,phase;
   int dn1,i0n1,dn2,i0n2;

   NewProfiler("local transpose",&prfTransLoc);
   NewProfiler("global transpose",&prfTransGlob);
   ProfileStart(prfTransGlob);
/*
 * get local dimensions of input and output arrays
 */
   MPI_Comm_size(comm,&nProc);
   MPI_Comm_rank(comm,&myRank);
   decomp1d(n1,nProc,myRank,&i0n1,&dn1);
   decomp1d(n2,nProc,myRank,&i0n2,&dn2);
/*
 * transpose local
 */
   ProfileStart(prfTransLoc);
   for ( i=0; i<dn1; i++ ) {
      for ( j=0; j<dn2; j++ ) {
         dcomplex *utji=&ut[(j*n1+i+i0n1)*n3];
         dcomplex *uij =&u[(i*n2+j+i0n2)*n3];
#pragma ivdep
         for ( k=0; k<n3; k++ ) {
#if 0
            ut[(j*n1+i+i0n1)*n3+k] = u[(i*n2+j+i0n2)*n3+k];
#endif
            utji[k].re = uij[k].re;
            utji[k].im = uij[k].im;
         }
      }
   }
   ProfileStop(prfTransLoc);
/*
 * Exchange pieces treating the processors as a ring.
 * In each phase, a box is sent phase processors up the ring
 * while one is received from the same distance down the ring. 
 */
   for ( phase=1; phase<nProc; phase++ ) {
      MPI_Datatype send_box,recv_slice,recv_box;
      MPI_Status stat;
      int d_dn1,d_i0n1,d_dn2,d_i0n2;
      int s_dn1,s_i0n1,s_dn2,s_i0n2;
      int dest,srce;

      dest = (myRank+phase) % nProc;
      srce = (nProc+myRank-phase) % nProc;
/*
 *    make receive box type from srce
 */
      decomp1d(n1,nProc,srce,&s_i0n1,&s_dn1);
      decomp1d(n2,nProc,srce,&s_i0n2,&s_dn2);
      MPI_Type_vector(dn2,n3,n1*n3,mpitype,&recv_slice);
      MPI_Type_hvector(s_dn1,1,n3*sizeof(*ut),recv_slice,&recv_box);
      MPI_Type_commit(&recv_box);
/*
 *    make send box type to dest
 */
      decomp1d(n1,nProc,dest,&d_i0n1,&d_dn1);
      decomp1d(n2,nProc,dest,&d_i0n2,&d_dn2);
      MPI_Type_vector(dn1,d_dn2*n3,n2*n3,mpitype,&send_box);
      MPI_Type_commit(&send_box);
/*
 *    send box up the ring, recv box from down the ring
 */
      MPI_Sendrecv(
         &u[d_i0n2*n3], 1,send_box,dest,0,
         &ut[s_i0n1*n3],1,recv_box,srce,0,comm,&stat);
/*
 *    free types
 */
      MPI_Type_free(&send_box);
      MPI_Type_free(&recv_box);
      MPI_Type_free(&recv_slice);
   }

   ProfileStop(prfTransGlob);
}








