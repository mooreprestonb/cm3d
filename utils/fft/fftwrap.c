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
 * It provides a common interface via the m1dFFT() function, which
 * performs simultaneous FFTs of multiple sequences. Valid values for
 * the length of the transform can be obtained from selectFFTGrid(),
 * which gives the nearest allowed length that is greater than or equal
 * to the requested length. By default, an internal radix-2 FFT routine
 * will be used if none of SGI_COMPLIB, IBM_ESSL, and HP_VECLIB is
 * defined. It will be slow and nasty, but not the worst it could be!
 * Beware that this internal routine is derived from Numerical Recipes
 * four1() so copyrights are involved.
 */
/********************************************************************/

#include "fftwrap.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* define this here so that we don't foist
   a particular naming convention on the rest of the
   world
   */
typedef struct 
{
   double re;
   double im;
} dcomplex;

/*
  Here are declarations for the vendor FFT and other library routines and
  some ugly magic inspired by IBM's essl.h that gives typical C calling
  semantics to some FORTRAN interfaces.
  */

static int __fwi1,__fwi2,__fwi3,__fwi4,__fwi5,__fwi6,__fwi7,__fwi8,__fwi9,__fwi10;
static double __fwd1;

#if defined(SGI_COMPLIB)
   dcomplex *zfftm1di(int n_elem,dcomplex *work);
   int zfftm1d(int isign,int n_elem,int n_seq,dcomplex *z,int inc_elem,int inc_seq,
	       dcomplex *work);
   int zscalm1d(int n_elem,int n_seq,double alpha,dcomplex *z,int inc_elem,int inc_seq);
   dcomplex *zfft3di(int,int,int,dcomplex *);
   int zfft3d(int,int,int,int,dcomplex *,int,int,dcomplex *);

#elif defined(IBM_ESSL)
   void dcft(int *init,dcomplex *x,int *incx,int *incsx,dcomplex *y,int *incy,int *incsy,
             int *n,int *m,int *isign,double *scale,double *work1,int *nw1,double *work2,
             int *nw2);
   void dcft3(dcomplex *x,int *i2x,int *i3x,dcomplex *y,int *i2y,int *i3y,int *n1,
              int *n2,int *n3,int *isign,double *scale,double *aux,int *naux);
#  define dcft(in,x,ix,sx,y,iy,sy,n,m,isign,scal,w1,nw1,w2,nw2) \
   ( __fwi1=in,__fwi2=ix,__fwi3=sx,__fwi4=iy,__fwi5=sy,__fwi6=n,__fwi7=m, \
     __fwi8=isign,__fwd1=scale,__fwi9=nw1,__fwi10=nw2, \
     dcft(&__fwi1,x,&__fwi2,&__fwi3,y,&__fwi4,&__fwi5,&__fwi6,&__fwi7,&__fwi8, \
          &__fwd1,w1,&__fwi9,w2,&__fwi10) )
#  define dcft3(x,i2x,i3x,y,i2y,i3y,n1,n2,n3,isign,scale,aux,naux) \
   ( __fwi1=i2x,__fwi2=i3x,__fwi3=i2y,__fwi4=i3y,__fwi5=n1,__fwi6=n2,__fwi7=n3, \
     __fwi8=isign,__fwd1=scale,__fwi9=naux, \
     dcft3(x,&__fwi1,&__fwi2,y,&__fwi3,&__fwi4,&__fwi5,&__fwi6,&__fwi7,&__fwi8, \
           &__fwd1,aux,&__fwi9) )
#elif defined(HP_VECLIB)
   extern void zffts(dcomplex *z,int *nelem,int *incelem,int *nseq,
                     int *incseq,int *iopt,int *ier);
   extern z3dfft(dcomplex *z,int *n1,int *n2,int *n3,int *ldz,int *mdz,int *iopt,
                 int *ier);
#  define zffts(z,nel,incel,nseq,incseq,iopt,ier) \
      ( __fwi1=nel,__fwi2=incel,__fwi3=nseq,__fwi4=incseq,__fwi5=iopt, \
        zffts(z,&__fwi1,&__fwi2,&__fwi3,&__fwi4,&__fwi5,ier) )
#  define z3dfft(z,n1,n2,n3,ldz,mdz,iopt,ier) \
      ( __fwi1=n1,__fwi2=n2,__fwi3=n3,__fwi4=ldz,__fwi5=mdz,__fwi6=iopt, \
        z3dfft(z,&__fwi1,&__fwi2,&__fwi3,&__fwi4,&__fwi5,&__fwi6,ier) )
#else
   static void NRm1dFFT(dcomplex *u,int nseq,int iseq,int nelem,int ielem,
      int isign);
   
#endif


/*
  Workspace for some of the FFT routines (allocated on first
  call). It should be big enough for most applications.
  */
static int nwork=20000;
static double *work[2];
static int work_space_allocated=0;
static void allocate_workspace(void);

static void multiple_zdscal(int n_seq,int inc_seq,int n_elem,int inc_elem,
                            double s,dcomplex *z);



int selectFFTGrid(int n)
/******************************************************************/
/* This function returns the smallest FFT size available that is
 * >= n.
 * Only sizes up to 128 are put in here, but larger ones could be
 * added. That would be unusually large for many 3d-FFT applications
 * but hardly unthinkable.
 */
{
#if defined(SGI_COMPLIB) || defined(HP_VECLIB)
   static const int size[] = { 0, 4,8,12,16,20,24,   30,32,36,40,   48,
                               60,64,   72,80,   90,96,100,108,120,128, -1 };
#elif defined(IBM_ESSL)
   static const int size[] = { 0, 4,8,12,16,20,24,28,30,32,36,40,42,48,
                               56,60,64,70,72,80,84,90,96,112,120,126,128, -1 };
#else
   static const int size[] = { 0, 4,8,16,32,64,128, -1 };
#endif
   int i;
   assert(n>0);
   assert(n<128);
   
   for ( i=0; size[i]>=0; i++ ) {
      if ( size[i]>=n ) break;
   }
   
   return size[i];
}


void m1dFFT(int isign,void *vz,int n_seq,int inc_seq,
   int n_elem,int inc_elem,double scale)
/******************************************************************/
/* This function does multiple 1-d FFTs on sequences of numbers in x.
 *    isign  = +/- 1 = sign of exponent in exp(isign*i*w*t)
 *    x              pointer to start of input data
 *    inc_elem       stride between elements of each sequence
 *    inc_seq        stride between the first elements of each sequence
 *    n_elem         number of elements within each sequence
 *    n_seq          number of sequences to transform
 */
{
   dcomplex *z = (dcomplex *) vz;
   int i,j,ier;

   allocate_workspace();

#if defined(IBM_ESSL)
   dcft(1,z,inc_elem,inc_seq,z,inc_elem,inc_seq,n_elem,n_seq,-isign,scale,work[0],
        nwork,work[1],nwork);
   dcft(0,z,inc_elem,inc_seq,z,inc_elem,inc_seq,n_elem,n_seq,-isign,scale,work[0],
        nwork,work[1],nwork);

#elif defined(SGI_COMPLIB)
   zfftm1di(n_elem,(dcomplex *)work[0]);
   zfftm1d(isign,n_elem,n_seq,z,inc_elem,inc_seq,(dcomplex *)work[0]);
   if ( scale!=1 ) {
      zscalm1d(n_elem,n_seq,scale,z,inc_elem,inc_seq);
   }

#elif defined(HP_VECLIB)
   if ( isign==-1 ) {
      zffts(z,n_elem,inc_elem,n_seq,inc_seq,1,&ier);
      if ( scale!=1 ) {
         multiple_zdscal(n_seq,inc_seq,n_elem,inc_elem,scale,z);
      }
   } else if ( isign==1 ) { /* note that zffts always scales here! */
      dcomplex *zseq=z;
      scale*=n_elem;
      zffts(z,n_elem,inc_elem,n_seq,inc_seq,-1,&ier);
      multiple_zdscal(n_seq,inc_seq,n_elem,inc_elem,scale,z);
   }

#else
   NRm1dFFT(z,n_seq,inc_seq,n_elem,inc_elem,isign);
   if ( scale!=1 ) {
      multiple_zdscal(n_seq,inc_seq,n_elem,inc_elem,scale,z);
   }
   
#endif   
}


extern void Three_d_FFT(int isign,dcomplex *z,int n1,int n2,int n3,
   int td2,int td3,double scale)
/********************************************************************/
{
   int i,j,ier;

   allocate_workspace();

#if defined(IBM_ESSL)
   dcft3(z,td3,td3*td2,z,td3,td3*td2,n3,n2,n1,-isign,scale,work[0],nwork);
#elif defined(SGI_COMPLIB)
   zfft3di(n3,n2,n1,(dcomplex *)work[0]);
   zfft3d(isign,n3,n2,n1,z,td3,td2,(dcomplex *)work[0]);
#elif defined(HP_VECLIB)
   if ( isign==-1 ) {
      z3dfft(z,n3,n2,n1,td3,td2,1,&ier);
      if ( scale!=1 ) {
         dcomplex *z1=z;
         int is1=td2*td3;
         for ( i=0; i<n1; i++,z1+=is1 ) {
            multiple_zdscal(n2,td3,n3,1,scale,z1);
         }
      }
   } else if ( isign==1 ) {
      z3dfft(z,n3,n2,n1,td3,td2,-1,&ier);
      if ( scale!=1/(double)(n1*n2*n3) ) {
         dcomplex *z1;
         int is1=td2*td3;
         scale*=n1*n2*n3;
         for ( i=0; i<n1; i++,z1+=is1 ) {
            multiple_zdscal(n2,td3,n3,1,scale,z1);
         }
      }         
   }
#else
#endif
}


static void allocate_workspace(void)
/*******************************************************************/
{
   if ( !work_space_allocated ) {
      double *pd = malloc(2*nwork*sizeof(*pd));
      assert(pd);
      work[0] = pd;
      work[1] = pd+nwork;
      work_space_allocated = 1;
   }
}

static void multiple_zdscal(int n_seq,int inc_seq,int n_elem,int inc_elem,
                            double s,dcomplex *z)
/************************************************************************/
/* Portable multiple scale dcomplex by double routine. It makes an
   attempt at good cache usage by looping over the smaller stride in
   the inner loop
   */
{
   if ( inc_seq<=inc_elem ) {
      dcomplex *zj=z;
      int j;
      for ( j=0; j<n_elem; j++,zj+=inc_elem ) {
         dcomplex *zi=zj;
         int i;
         for ( i=0; i<n_seq; i++,zi+=inc_seq ) {
            zi->re *= s;
            zi->im *= s;
         }
      }
   } else {
      dcomplex *zi=z;
      int i;
      for ( i=0; i<n_seq; i++,zi+=inc_seq ) {
         dcomplex *zj=zi;
         int j;
         for ( j=0; j<n_elem; j++,zj+=inc_elem ) {
            zj->re *= s;
            zj->im *= s;
         }
      }
   }
}


#if !(defined(IBM_ESSL) || defined(SGI_COMPLIB) || defined(HP_VECLIB))

#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265359
#endif

static void NRm1dFFT(dcomplex *u,int nseq,int iseq,int nelem,
   int ielem,int isign)
/********************************************************************/
/* A somewhat hacked-up Numerical Recipes routine (c.f. four1()) that
   does nseq FFT's of nelem elements each. The stride between the
   starting points of each sequence is iseq. The stride between
   elements of a sequence is ielem.

   The inner loop is over the sequences, so this generally gives
   better cache behavior if iseq<ielem (but not always).
   */
{
   int n,i,istep,j,mmax,n2;

   n = nelem;
   n2 = n/2;
   j  = n2;

   /* put data into bit-reverse order */
   for ( i=1; i<=n-2; i++ ) {
      int m;
      if ( j>i ) {
         int ii=ielem*i;
         int jj=ielem*j;
         int is;
         for ( is=0; is<nseq; is++ ) {
            dcomplex temp;
            int iss=iseq*is;
            temp = u[ii+iss];
            u[ii+iss]=u[jj+iss];
            u[jj+iss]=temp;
         }
      }
      m = n2;
      for (;;) {
         if ( m<2 || j<m ) break;
         j -= m;
         m /= 2;
      }
      j += m;
   }

   /* do multiple FFT's with inner loop
      over sequences */
   mmax=1;
   while ( mmax<n ) {
      dcomplex w,wp,wt;
      double theta,s;
      int m;
      istep = 2*mmax;
      theta = PI/(isign*mmax);
      s     = sin(0.5*theta);
      wp.re = -2*s*s;
      wp.im = sin(theta);
      w.re  = 1;
      w.im  = 0;
      for ( m=1; m<=mmax; m++ ) {
         dcomplex ws = w;
         for ( i=m; i<=n; i+=istep ) {
            int is,ii,jj;
            j = i+mmax;
            ii = ielem*(i-1);
            jj = ielem*(j-1);
            for ( is=0; is<nseq; is++ ) {
               dcomplex temp;
               int iss=iseq*is;
               temp.re = ws.re*u[jj+iss].re-ws.im*u[jj+iss].im;
               temp.im = ws.re*u[jj+iss].im+ws.im*u[jj+iss].re;
               u[jj+iss].re = u[ii+iss].re-temp.re;
               u[jj+iss].im = u[ii+iss].im-temp.im;
               u[ii+iss].re += temp.re;
               u[ii+iss].im += temp.im;
            }
         }
         wt.re = (w.re*wp.re-w.im*wp.im)+w.re;
         wt.im = (w.re*wp.im+w.im*wp.re)+w.im;
         w     = wt;
      }
      mmax = istep;
   }

}


#endif




