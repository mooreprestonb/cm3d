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
/* Normal C Path   */
 
#ifndef __essl
 
  #define __essl 1
 
/**********************************************************************
 *  Licensed Materials - Property of IBM                              *
*    LICENSED MATERIALS - PROPERTY OF IBM                             *
*    THIS MODULE IS "RESTRICTED MATERIALS OF IBM"                     *
*                                                                     *
*    COPYRIGHT = 5765-042 (C) COPYRIGHT IBM CORP. 1991, 1994.         *
*    ALL RIGHTS RESERVED.                                             *
*                                                                     *
*    U.S. GOVERNMENT USERS RESTRICTED RIGHTS - USE, DUPLICATION       *
*    OR DISCLOSURE RESTRICTED BY GSA ADP SCHEDULE CONTRACT WITH       *
*    IBM CORP.                                                        *
*    SEE COPYRIGHT INSTRUCTIONS.                                      *
*                                                                     *
*    THE SOURCE CODE FOR THIS PROGRAM IS NOT PUBLISHED OR OTHERWISE   *
*    DIVESTED OF ITS TRADE SECRETS, IRRESPECTIVE OF WHAT HAS BEEN     *
*    DEPOSITED WITH THE U.S. COPYRIGHT OFFICE.                        *
 *                                                                    *
 *  Program name - <essl.h> header file                               *
 *  Descriptive name - ESSL release 4 C language header file          *
 *                                                                    *
 *  Function : This file must be included in any C file               *
 *             containing ESSL calls in order for the ESSL calls      *
 *             to work as documented in the ESSL Guide and            *
 *             Reference (SC23-0526).                                 *
 *  Change activity -                                                 *
 *                    Added new routines for V2.2                     *
 *                    Total = 425+16=441                (5/93)        *
 *                    Added iessl()                      3/94         *
 *                    Added extern for C++ path          6/94         *
 *                    and fixed esvi4=i3 problem         6/94         *
 *                    If CMPLX and DCMPLX not defined, defined  5/95  *
 *********************************************************************/
 
/*  Definition of complex data types */
 
#ifndef  _CMPLX
#define  _CMPLX 1
#ifndef  _REIM
#define _REIM   1
#endif
typedef union { struct { float  __re, __im; } __data; double __align; }  cmplx;
#endif
 
#ifndef  _DCMPLX
#define  _DCMPLX 1
#ifndef  _REIM
#define _REIM   1
#endif
typedef union { struct { double __re, __im; } __data; double __align; } dcmplx;
#endif
 
#ifdef  _REIM
#define RE(x)   ((x).__data.__re)
#define IM(x)   ((x).__data.__im)
#endif
 
 
/* The following temporary variables are used to create fortran linkage */
 
int  __esvi1,__esvi2,__esvi3,__esvi4,__esvi5;
int  __esvi6,__esvi7,__esvi8,__esvi9,__esvi10,__esvi11;
int  __esvi12,__esvi13;
float __esvf1,__esvf2,__esvf3,__esvf4;
double __esvd1,__esvd2,__esvd3,__esvd4;
cmplx __esvc1,__esvc2,__esvc3,__esvc4,__esvctmp;
dcmplx __esvdc1,__esvdc2,__esvdc3,__esvdc4,__esvdtmp;
 
void  (*__ESVFP)();
void  (*__ESVFP2)();
void  (*__ESVFP3)();
 
#ifdef __ESVERR
#define __ESVI    int
#else
#define __ESVI    void
#endif
 
/*  Linear Algebra Subprograms  */
 
/*  Vector-Scalar Subprograms  */
int isamax(int*, float *, int*);
int idamax(int*, double *, int*);
int icamax(int*,  cmplx *, int*);
int izamax(int*, dcmplx *, int*);
 
int isamin(int*, float *, int*);
int idamin(int*, double *, int*);
 
int ismax(int*, float *, int*);
int idmax(int*, double *, int*);
 
int ismin(int*, float *, int*);
int idmin(int*, double *, int*);
 
int iessl(void);
 
float   sasum(int*,  float *, int*);
double  dasum(int*, double *, int*);
float  scasum(int*,  cmplx *, int*);
double dzasum(int*, dcmplx *, int*);
 
void   saxpy(int*,  float*,  float *, int*,  float *, int*);
void   daxpy(int*, double*, double *, int*, double *, int*);
void   caxpy(int*,  cmplx*,  cmplx *, int*,  cmplx *, int*);
void   zaxpy(int*, dcmplx*, dcmplx *, int*, dcmplx *, int*);
 
void   scopy(int*,  float *, int*,  float *, int*);
void   dcopy(int*, double *, int*, double *, int*);
void   ccopy(int*,  cmplx *, int*,  cmplx *, int*);
void   zcopy(int*, dcmplx *, int*, dcmplx *, int*);
 
float  sdot(int*, float *, int*, float *, int*);
double ddot(int*, double *, int*, double *, int*);
 
void cdotu( int*,  cmplx *, int*,  cmplx *, int*, cmplx *);
void zdotu( int*, dcmplx *, int*, dcmplx *, int*, dcmplx *);
void cdotc( int*, cmplx  *, int*, cmplx  *, int*, cmplx *);
void zdotc( int*, dcmplx *, int*, dcmplx *, int*, dcmplx *);
 
 
void   snaxpy(int*,int*,float*, int*, void *, int*, int*,
              void *, int*, int*);
void   dnaxpy(int*,int*,double*, int*,void *, int*, int*,
              void *, int*, int*);
 
void   sndot(int *,int *,float *,int *,int *,void *,
             int *,int *,void *,int *,int *);
void   dndot(int *,int *,double *,int *,int *,void *,int *,
             int *,void *,int *,int *);
 
float snrm2(int*, float *, int*);
double dnrm2(int*, double *, int*);
float  scnrm2(int*,  cmplx *, int*);
double dznrm2(int*, dcmplx *, int*);
 
float  snorm2(int*,  float *, int*);
double dnorm2(int*, double *, int*);
float  cnorm2(int*,  cmplx *, int*);
double znorm2(int*, dcmplx *, int*);
 
void   srotg( float *,  float *,  float *,  float *);
void   drotg(double *, double *, double *, double *);
void   crotg( cmplx *,  cmplx *,  float *,  cmplx *);
void   zrotg(dcmplx *, dcmplx *, double *, dcmplx *);
 
void    srot(int*,  float *, int*,  float *, int*,  float*,  float*);
void    drot(int*, double *, int*, double *, int*, double*, double*);
void    crot(int*,  cmplx *, int*,  cmplx *, int*,  float*,  cmplx*);
void    zrot(int*, dcmplx *, int*, dcmplx *, int*, double*, dcmplx*);
void   csrot(int*,  cmplx *, int*,  cmplx *, int*,  float*,  float*);
void   zdrot(int*, dcmplx *, int*, dcmplx *, int*, double*, double*);
 
void    sscal(int*,  float*,  float *, int*);
void    dscal(int*, double*, double *, int*);
void    cscal(int*,  cmplx*,  cmplx *, int*);
void    zscal(int*, dcmplx*, dcmplx *, int*);
void   csscal(int*,  float*,  cmplx *, int*);
void   zdscal(int*, double*, dcmplx *, int*);
 
void   sswap(int*,  float *, int*,  float *, int*);
void   dswap(int*, double *, int*, double *, int*);
void   cswap(int*,  cmplx *, int*,  cmplx *, int*);
void   zswap(int*, dcmplx *, int*, dcmplx *, int*);
 
void    syax(int*,  float*,  float *, int*,  float *, int*);
void    dyax(int*, double*, double *, int*, double *, int*);
void    cyax(int*,  cmplx*,  cmplx *, int*,  cmplx *, int*);
void    zyax(int*, dcmplx*, dcmplx *, int*, dcmplx *, int*);
void   csyax(int*,  float*,  cmplx *, int*,  cmplx *, int*);
void   zdyax(int*, double*, dcmplx *, int*, dcmplx *, int*);
 
void   szaxpy(int*,  float *,  float *, int*,  float *, int*,
                     float *, int*);
void   dzaxpy(int*, double *, double *, int*, double *, int*,
                     double *, int*);
void   czaxpy(int*,  cmplx *,  cmplx *, int*,  cmplx *, int*,
                     cmplx *, int*);
void   zzaxpy(int*, dcmplx *, dcmplx *, int*, dcmplx *, int*,
                    dcmplx *, int*);
 
void svea(int *,  float *, int *,  float *, int *,  float *, int *);
void dvea(int *, double *, int *, double *, int *, double *, int *);
void cvea(int *,  cmplx *, int *,  cmplx *, int *,  cmplx *, int *);
void zvea(int *, dcmplx *, int *, dcmplx *, int *, dcmplx *, int *);
 
void sves(int *,  float *, int *,  float *, int *,  float *, int *);
void dves(int *, double *, int *, double *, int *, double *, int *);
void cves(int *,  cmplx *, int *,  cmplx *, int *,  cmplx *, int *);
void zves(int *, dcmplx *, int *, dcmplx *, int *, dcmplx *, int *);
 
void svem(int *,  float *, int *,  float *, int *,  float *, int *);
void dvem(int *, double *, int *, double *, int *, double *, int *);
void cvem(int *,  cmplx *, int *,  cmplx *, int *,  cmplx *, int *);
void zvem(int *, dcmplx *, int *, dcmplx *, int *, dcmplx *, int *);
 
/*  Sparse Vector-Scalar Subroutines  */
 
void   ssctr(int*,  float *, int *,  float *);
void   dsctr(int*, double *, int *, double *);
void   csctr(int *,  cmplx *, int *,  cmplx *);
void   zsctr(int *, dcmplx *, int *, dcmplx *);
 
void   sgthr(int*,  float *,  float *, int *);
void   dgthr(int*, double *, double *, int *);
void   cgthr(int *,  cmplx *,  cmplx *, int *);
void   zgthr(int *, dcmplx *, dcmplx *, int *);
 
void   sgthrz(int*,  float *,  float *, int *);
void   dgthrz(int*, double *, double *, int *);
void   cgthrz(int *,  cmplx *,  cmplx *, int *);
void   zgthrz(int *, dcmplx *, dcmplx *, int *);
 
void   saxpyi(int*,  float*,  float *, int *,  float *);
void   daxpyi(int*, double*, double *, int *, double *);
void   caxpyi(int *,  cmplx *,  cmplx *, int *,  cmplx *);
void   zaxpyi(int *, dcmplx *, dcmplx *, int *, dcmplx *);
 
float  sdoti(int*,  float *, int *,  float *);
double ddoti(int*, double *, int *, double *);
 
cmplx cdotui(int *,  cmplx *, int *,  cmplx *);
dcmplx zdotui(int *, dcmplx *, int *, dcmplx *);
 
cmplx cdotci(int *,  cmplx *, int *,  cmplx *);
dcmplx zdotci(int *, dcmplx *, int *, dcmplx *);
 
/*  Dense Matrix-Vector Subroutines  */
 
void   sgemv(char *,int*, int*,  float*, void *, int*,  float *, int*,  float*,
             float *, int*);
void   dgemv(char *,int*, int*, double*, void *, int*, double *, int*, double*,
             double *, int*);
void   cgemv(char *,int*, int*,  cmplx*, void *, int*,  cmplx *, int*,  cmplx*,
             cmplx *, int*);
void   zgemv(char *,int*, int*, dcmplx*, void *, int*, dcmplx *, int*, dcmplx*,
             dcmplx *, int*);
 
void   sgemx(int *,int *,float *,void *,int *,float *,int *,
             float *,int *);
void   dgemx(int*,int *,double *,void *,int *, double *,int *,
             double *,int *);
 
void   sgemtx(int *,int *,float *,void *,int *,float *,int *,
              float *,int *);
void   dgemtx(int *,int *,double *,void *,int *,double *,int *,
              double *,int *);
 
#define sger  sger1
#define dger  dger1
void   sger1(int *,int *,float *,float *,int *,float *,
             int *, void *,int *);
void   dger1(int *,int *,double *,double *,int *,double *,
             int *,void *,int *);
void   cgeru(int *,int *,cmplx *,cmplx *,int *,cmplx *,int *,
             void *,int*);
void   zgeru(int*,int *,dcmplx *,dcmplx *,int *,dcmplx *,int *,
             void *,int*);
void   cgerc(int *,int *,cmplx *,cmplx *,int *,cmplx *,int *,
             void *,int *);
void   zgerc(int *,int *,dcmplx *,dcmplx *,int *,dcmplx *,int *,
             void *,int *);
 
void   sslmx(int *,float *,float *,float *,int *,float *,int *);
void   dslmx(int *,double *,double *,double *,int *,double *,int *);
 
void   sslr1(int *,  float *,  float *, int *,  float *);
void   dslr1(int *, double *, double *, int *, double *);
 
void   sslr2(int *, float *, float *, int *, float *, int *, float *);
void   dslr2(int *, double *, double *, int *, double *, int*, double *);
void strmv(char *, char *, char *, int *, void *, int *, float *, int *);
void dtrmv(char *, char *, char *, int *, void *, int *, void *, int *);
void ctrmv(char *, char *, char *, int *, void *, int *, void *, int *);
void ztrmv(char *, char *, char *, int *, void *, int *, void *, int *);
 
void sspmv(char *, int *,  float *,  float *,  float *, int *,  float *,
            float *, int *);
void dspmv(char *, int *, double *, double *, double *, int *, double *,
           double *, int *);
void chpmv(char *, int *,  cmplx *,  cmplx *,  cmplx *, int *,  cmplx *,
            cmplx *, int *);
void zhpmv(char *, int *, dcmplx *, dcmplx *, dcmplx *, int *, dcmplx *,
           dcmplx *, int *);
 
void ssymv(char *, int *,  float *, void *, int *,  float *, int *,
            float *,  float *, int *);
void dsymv(char *, int *, double *, void *, int *, double *, int *,
           double *, double *, int *);
void chemv(char *, int *,  cmplx *, void *, int *,  cmplx *, int *,
            cmplx *,  cmplx *, int *);
void zhemv(char *, int *, dcmplx *, void *, int *, dcmplx *, int *,
           dcmplx *, dcmplx *, int *);
 
void sspr(char *, int *,  float *,  float *, int *,   float *);
void dspr(char *, int *, double *, double *, int *,  double *);
void chpr(char *, int *,  float *,  cmplx *, int *,   cmplx *);
void zhpr(char *, int *, double *, dcmplx *, int *,  dcmplx *);
 
void ssyr(char *, int *,  float *,  float *, int *,  void *, int *);
void dsyr(char *, int *, double *, double *, int *,  void *, int *);
void cher(char *, int *,  float *,  cmplx *, int *,  void *, int *);
void zher(char *, int *, double *, dcmplx *, int *,  void *, int *);
 
void sspr2(char *, int *,  float *,  float *, int *, float *, int *,
            float *);
void dspr2(char *, int *, double *, double *, int *,double *, int *,
           double *);
void chpr2(char *, int *,  cmplx *,  cmplx *, int *, cmplx *, int *,
            cmplx *);
void zhpr2(char *, int *, dcmplx *, dcmplx *, int *,dcmplx *, int *,
           dcmplx *);
 
void ssyr2(char *, int *,  float *,  float *, int *,  float *, int *,
           void *, int *);
void dsyr2(char *, int *, double *, double *, int *, double *, int *,
           void *, int *);
void cher2(char *, int *,  cmplx *,  cmplx *, int *,  cmplx *, int *,
           void *, int *);
void zher2(char *, int *, dcmplx *, dcmplx *, int *, dcmplx *, int *,
           void *, int *);
 
void sgbmv(char *, int *, int *, int *, int *,  float *, void *, int *,
            float *, int *,  float *,  float *, int *);
void dgbmv(char *, int *, int *, int *, int *, double *, void *, int *,
           double *, int *, double *, double *, int *);
void cgbmv(char *, int *, int *, int *, int *,  cmplx *, void *, int *,
            cmplx *, int *,  cmplx *,  cmplx *, int *);
void zgbmv(char *, int *, int *, int *, int *, dcmplx *, void *, int *,
           dcmplx *, int *, dcmplx *, dcmplx *, int *);
 
void ssbmv(char *, int *, int *,  float *, void *, int *,  float *,
           int *,  float *,  float *, int *);
void dsbmv(char *, int *, int *, double *, void *, int *, double *,
           int *, double *, double *, int *);
void chbmv(char *, int *, int *,  cmplx *, void *, int *,  cmplx *,
           int *,  cmplx *,  cmplx *, int *);
void zhbmv(char *, int *, int *, dcmplx *, void *, int *, dcmplx *,
           int *, dcmplx *, dcmplx *, int *);
 
void stbmv(char *, char *, char *, int *, int *, void *, int *, float *,
           int *);
void dtbmv(char *, char *, char *, int *, int *, void *, int *,double *,
           int *);
void ctbmv(char *, char *, char *, int *, int *, void *, int *, cmplx *,
           int *);
void ztbmv(char *, char *, char *, int *, int *, void *, int *,dcmplx *,
           int *);
 
void stpmv(char *, char *, char *, int *,  float *,  float *, int *);
void dtpmv(char *, char *, char *, int *, double *, double *, int *);
void ctpmv(char *, char *, char *, int *,  cmplx *,  cmplx *, int *);
void ztpmv(char *, char *, char *, int *, dcmplx *, dcmplx *, int *);
 
/*  Sparse Matrix-Vector Subroutines  */
 
void   dsmmx(int *, int *, void *, void *, int *, double *, double *);
 
__ESVI dsmtm(int *, int *, void *, void *, int *, int *, int *, void *, void *,
             int *, float *, int *);
 
void   dsdmx(int *, int *, int *, void *, int *, char *, int *,
             double *, double *);
 
/*  Matrix Operation Subroutines  */
 
void   sgeadd(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *);
void   dgeadd(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *);
void   cgeadd(void *, int *, char *, void *, int *, char *, void *, int *,
              int *, int *);
void   zgeadd(void *, int *, char *, void *, int *, char *, void *, int *,
              int *, int *);
 
void   sgesub(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *);
void   dgesub(void *, int *, char *, void *, int *, char *, void *, int *,
              int *, int *);
void   cgesub(void *, int *, char *, void *, int *, char *, void *, int *,
              int *, int *);
void   zgesub(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *);
 
void   sgemul(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *, int *);
void   dgemul(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *, int *);
void   cgemul(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int*, int*);
void   zgemul(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *, int *);
 
__ESVI sgemms(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *, int *,  float *, int *);
__ESVI dgemms(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *, int*, double *, int *);
__ESVI cgemms(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *, int *,  float *, int *);
__ESVI zgemms(void *, int *, char *, void *, int *, char *, void *, int *,
               int *, int *, int *, double *, int *);
 
void   sgemm(char *, char *, int *, int *, int *,  float *, void *, int *,
             void *, int *,  float *, void *, int *);
void   dgemm(char *, char *, int *, int *, int *, double *, void *, int *,
             void *, int *, double *, void *, int *);
void   cgemm(char *, char *, int *, int *, int *,  cmplx *, void *, int *,
             void *, int *,  cmplx *, void *, int *);
void   zgemm(char *, char *, int *, int *, int *, dcmplx *, void *, int *,
             void *, int *, dcmplx *, void *, int *);
 
void   ssyrk(char *, char *, int *, int *, float *, void *, int *, float *,
             void *, int *);
void   dsyrk(char *, char *, int *, int *, double *, void *, int *, double *,
             void *, int *);
 
void   sgetmi(void *, int *, int *);
void   dgetmi(void *, int *, int *);
void   cgetmi(void *, int *, int *);
void   zgetmi(void *, int *, int *);
 
void   sgetmo(void *, int *, int *, int *, void *, int *);
void   dgetmo(void *, int *, int *, int *, void *, int *);
void   cgetmo(void *, int *, int *, int *, void *, int *);
void   zgetmo(void *, int *, int *, int *, void *, int *);
void strmm(char *, char *, char *, char *, int *, int *, float *, void*,
                                                  int *, void *, int *);
void dtrmm(char *, char *, char *, char *, int *, int *, double *, void*,
                                                  int *, void *, int *);
void ctrmm(char *, char *, char *, char *, int *, int *, cmplx *, void*,
                                                  int *, void *, int *);
void ztrmm(char *, char *, char *, char *, int *, int *, dcmplx *, void*,
                                                  int *, void *, int *);
void ssymm(char *, char *, int *, int *, float *, void *, int *, void *,
                                           int *, float *, void *, int *);
void dsymm(char *, char *, int *, int *, double *, void *, int *, void *,
                                           int *, double *, void *, int *);
void dsyr2k(char *, char *, int *, int *, double *, void *, int *, void *,
            int *, double *, void *, int *);
void ssyr2k(char *, char *, int *, int *, float *, void *, int *, void *,
            int *, float *, void *, int *);
 
void csyrk(char *, char *, int *, int *,  cmplx *, void *, int *,
            cmplx *, void *, int *);
void zsyrk(char *, char *, int *, int *, dcmplx *, void *, int *,
           dcmplx *, void *, int *);
 
void cherk(char *, char *, int *, int *,  float *, void *, int *,
            float *, void *, int *);
void zherk(char *, char *, int *, int *, double *, void *, int *,
           double *, void *, int *);
 
void csyr2k(char *, char *, int *, int *,  cmplx *, void *, int *,
            void *, int *, cmplx *, void *, int *);
void zsyr2k(char *, char *, int *, int *, dcmplx *, void *, int *,
            void *, int *, dcmplx *, void *, int *);
void cher2k(char *, char *, int *, int *,  cmplx *, void *, int *,
            void *, int *,  float *, void *, int *);
void zher2k(char *, char *, int *, int *, dcmplx *, void *, int *,
            void *, int *, double *, void *, int *);
 
 
void csymm(char *, char *, int *, int *,  cmplx *, void *, int *, void *,
           int *,   cmplx *, void *, int *);
void zsymm(char *, char *, int *, int *, dcmplx *, void *, int *, void *,
           int *,  dcmplx *, void *, int *);
void chemm(char *, char *, int *, int *,  cmplx *, void *, int *, void *,
           int *,   cmplx *, void *, int *);
void zhemm(char *, char *, int *, int *, dcmplx *, void *, int *, void *,
           int *,  dcmplx *, void *, int *);
 
 
/*  Dense Linear Algebraic Equation Subroutines  */
 
__ESVI sgef(void *, int *, int *, int *);
__ESVI dgef(void *, int *, int *, int *);
__ESVI cgef(void *, int *, int *, int *);
__ESVI zgef(void *, int *, int *, int *);
 
void   sges(void *, int *, int *, int *,  float *, int *);
void   dges(void *, int *, int *, int *, double *, int *);
void   cges(void *, int *, int *, int *,  cmplx *, int *);
void   zges(void *, int *, int *, int *, dcmplx *, int *);
 
void sgesm(char *, void *, int *, int *, int *, void *, int *, int *);
void dgesm(char *, void *, int *, int *, int *, void *, int *, int *);
void cgesm(char *, void *, int *, int *, int *, void *, int *, int *);
void zgesm(char *, void *, int *, int *, int *, void *, int *, int *);
 
__ESVI sgefcd(void *, int *, int *, int *, int *, float *, float *, float *,
               int *);
__ESVI dgefcd(void *, int *, int *, int *, int *, double *, double *, double *,
               int *);
 
__ESVI sppf( float *, int *, int *);
__ESVI dppf(double *, int *, int *);
 
void   spps( float *, int *,  float *, int *);
void   dpps(double *, int *, double *, int *);
 
__ESVI sppfcd( float *, int *, int *, float *, float *, float *, int *);
__ESVI dppfcd(double *,int *, int *, double *, double *, double *, int *);
 
__ESVI sgeicd(void *, int *, int *, int *,  float *,  float *,  float *,
               int *);
__ESVI dgeicd(void *, int *, int *, int *, double *, double *, double *,
               int *);
 
__ESVI sppicd( float *, int *, int *,  float *,  float *,  float *,
               int *);
__ESVI dppicd(double *, int *, int *, double *, double *, double *,
               int *);
 
void   strsm(char *, char *, char *, char *, int *, int *,  float *, void *,
              int *, void *, int *);
void   dtrsm(char *, char *, char *, char *, int *, int *, double *, void *,
              int *, void *, int *);
void ctrsm(char *, char *, char *, char *, int *, int *, cmplx *, void *,
            int *, void *, int *);
void ztrsm(char *, char *, char *, char *, int *, int *, dcmplx *, void *,
            int *, void *, int *);
void strsv(char *, char *, char *, int *, void *, int *, void *, int *);
void dtrsv(char *, char *, char *, int *, void *, int *, void *, int *);
void ctrsv(char *, char *, char *, int *, void *, int *, void *, int *);
void ztrsv(char *, char *, char *, int *, void *, int *, void *, int *);
 
 
void spof(char *, void *, int *, int *);
void dpof(char *, void *, int *, int *);
void cpof(char *, void *, int *, int *);
void zpof(char *, void *, int *, int *);
 
void sposm(char *, void *, int *, int *, void *, int *, int *);
void dposm(char *, void *, int *, int *, void *, int *, int *);
void cposm(char *, void *, int *, int *, void *, int *, int *);
void zposm(char *, void *, int *, int *, void *, int *, int *);
 
__ESVI spofcd(char *, void *, int *, int *, int *,  float *,  float *,
             float *, int *);
__ESVI dpofcd(char *, void *, int *, int *, int *, double *, double *,
            double *, int *);
 
__ESVI spoicd(char *, void *, int *, int *, int *,  float *,  float *,
             float *, int *);
__ESVI dpoicd(char *, void *, int *, int *, int *, double *, double *,
            double *, int *);
 
__ESVI stri(char *, char *, void *, int *, int *);
__ESVI dtri(char *, char *, void *, int *, int *);
 
__ESVI stpi(char *, char *,  float *, int *);
__ESVI dtpi(char *, char *, double *, int *);
 
 
void stpsv(char *, char *, char *, int *,  float *,  float *, int *);
void dtpsv(char *, char *, char *, int *, double *, double *, int *);
void ctpsv(char *, char *, char *, int *,  cmplx *,  cmplx *, int *);
void ztpsv(char *, char *, char *, int *, dcmplx *, dcmplx *, int *);
 
void stbsv(char *, char *, char *, int *, int *, void *, int *, float *,
           int *);
void dtbsv(char *, char *, char *, int *, int *, void *, int *,double *,
           int *);
void ctbsv(char *, char *, char *, int *, int *, void *, int *, cmplx *,
           int *);
void ztbsv(char *, char *, char *, int *, int *, void *, int *,dcmplx *,
           int *);
 
void cgtnpf(int *,  cmplx *,  cmplx *,  cmplx *, int *);
void zgtnpf(int *, dcmplx *, dcmplx *, dcmplx *, int *);
 
void cgtnps(int *,  cmplx *,  cmplx *,  cmplx *,  cmplx *);
void zgtnps(int *, dcmplx *, dcmplx *, dcmplx *, dcmplx *);
 
__ESVI dsris(char *, char *, int *, double *, int *, int *, double *,
           double *, int *, double *, double *, int *, double *, int *);
 
/*  Banded Linear Algebraic Equation Subroutines  */
 
__ESVI sgbf(void *, int *, int *, int *, int *, int *);
__ESVI dgbf(void *, int *, int *, int *, int *, int *);
 
void   sgbs(void *, int *, int *, int *, int *, int *,  float *);
void   dgbs(void *, int *, int *, int *, int *, int *, double *);
 
__ESVI spbf(void *, int *, int *, int *);
__ESVI dpbf(void *, int *, int *, int *);
 
void   spbs(void *, int *, int *, int *,  float *);
void   dpbs(void *, int *, int *, int *, double *);
 
__ESVI spbchf(void *, int *, int *, int *);
__ESVI dpbchf(void *, int *, int *, int *);
 
void   spbchs(void *, int *, int *, int *,  float *);
void   dpbchs(void *, int *, int *, int *, double *);
 
__ESVI sgtf(int *,  float *,  float *,  float *,  float *, int *);
__ESVI dgtf(int *, double *, double *, double *, double *, int *);
 
void   sgts(int *,  float *,  float *,  float *,  float *, int *,
            float *);
void   dgts(int *, double *, double *, double *, double *, int *,
            double *);
 
void   sgtnpf(int *,  float *,  float *,  float *, int *);
void   dgtnpf(int *, double *, double *, double *, int *);
 
void   sgtnps(int *,  float *,  float *,  float *,  float *);
void   dgtnps(int *, double *, double *, double *, double *);
 
void   sgtnp(int *,  float *,  float *,  float *,  float *);
void   dgtnp(int *, double *, double *, double *, double *);
void   cgtnp(int *,  cmplx *,  cmplx *,  cmplx *,  cmplx *);
void   zgtnp(int *, dcmplx *, dcmplx *, dcmplx *, dcmplx *);
 
void   sptf(int *,  float *,  float *, int*);
void   dptf(int *, double *, double *, int*);
 
void   spts(int *,  float *,  float *,  float *);
void   dpts(int *, double *, double *, double *);
 
/*  Sparse Linear Algebraic Equation Subroutines  */
 
__ESVI dgsf(int *,int *, int *, double *, int *, int *, int *, int *, double *,
            double *, double *, int *);
 
__ESVI dgss(int *, int *, double *, int *, int *, int *, double *, double *,
            int *);
 
__ESVI dgkfs(int *, double *, int *, int *, double *, int *, int *, int *,
             double *, double *, int *, void *, int *, int *);
 
__ESVI dskfs(int *, double *, int *, int *, int *, double *, double *,
             int *, void *, int *, int *);
 
__ESVI dsmcg(int *, int *, void *, void *, int *, double *, double *, int *,
           double *, double *, int *, double *, int *);
 
__ESVI dsdcg(int *, int *, int *, void *, int *, int *, double *, double *,
             int *, double *, double *, int *, double *, int *);
 
__ESVI dsmgcg(int *, int *, void *, void *, int *, double *, double *, int *,
              double *, double *, int *, double *, int *);
 
__ESVI dsdgcg(int *, int *, void *, int *, int *, double *, double *, int *,
              double *, double *, int *, double *, int *);
 
/*  Linear Least Squares Subroutines  */
 
__ESVI sgesvf(int *, void *, int *, void *, int *, int *,  float *, int *,
              int *, float *, int *);
__ESVI dgesvf(int *, void *, int *, void *, int *, int *, double *, int *,
              int *, double *, int *);
 
void   sgesvs(void *, int *, void *, int *, int *,  float *, void *, int *,
               int *, int *,  float *);
void   dgesvs(void *, int *, void *, int *, int *, double *, void *, int *,
               int *, int *, double *);
 
__ESVI sgells(int *, void *, int *, void *, int *, void *, int *,  float *,
              float *, int *, int *, int *, int *,  float *, int *);
__ESVI dgells(int *, void *, int *, void *, int *, void *, int *, double *,
              double *, int *, int *, int *, int *, double *, int *);
 
 
/*  Eigensystem Analysis Subroutines  */
 
 
__ESVI sgeev(int *, void *, int *,  cmplx *, void *, int *, int *, int *,
             float *, int *);
__ESVI dgeev(int *, void *, int *, dcmplx *, void *, int *, int *, int *,
             double *, int *);
__ESVI cgeev(int *, void *, int *,  cmplx *, void *, int *, int *, int *,
             float *, int *);
__ESVI zgeev(int *, void *, int *, dcmplx *, void *, int *, int *, int *,
             double *, int *);
 
#define sslev  sspev
#define dslev  dspev
#define chlev  chpev
#define zhlev  zhpev
__ESVI sspev(int *,  float *,  float *, void *, int *, int *,  float *,
             int *);
__ESVI dspev(int *, double *, double *, void *, int *, int *, double *,
             int *);
__ESVI chpev(int *,  cmplx *,  float *, void *, int *, int *,  float *,
             int *);
__ESVI zhpev(int *, dcmplx *, double *, void *, int *, int *, double *,
             int *);
 
__ESVI sspsv(int *,  float *,  float *, void *, int *, int *, int *,  float *,
             int *);
__ESVI dspsv(int *, double *, double *, void *, int *, int *, int *, double *,
             int *);
__ESVI chpsv(int *,  cmplx *,  float *, void *, int *, int *, int *,  float *,
             int *);
__ESVI zhpsv(int *, dcmplx *, double *, void *, int *, int *, int *, double *,
             int *);
 
__ESVI sgegv(int *, void *, int *, void *, int *,  cmplx *,  float *, void *,
             int *, int *,  float *, int *);
__ESVI dgegv(int *, void *, int *, void *, int *, dcmplx *, double *, void *,
             int *, int *, double *, int *);
 
__ESVI ssygv(int *, void *, int *, void *, int *,  float *, void *, int *,
             int *, float *, int *);
__ESVI dsygv(int *, void *, int *, void *, int *, double *, void *, int *,
             int *, double *, int *);
 
/*  Fourier Transform Subroutines  */
 
__ESVI scft(int *, void *, int *, int *, void *, int *, int *, int *, int *,
            int *, float *, double *, int *, double *, int *);
 
__ESVI dcft(int *, void *, int *, int *, void *, int *, int *, int *, int *,
            int *, double *, double *, int *, double *, int *);
 
__ESVI srcft(int *, void *, int *, void *, int *, int *, int *, int *, float *,
             double *, int *, double *, int *, double *, int *);
 
__ESVI drcft(int *,void *, int *, void *, int *, int *, int *, int *, double *,
             double *, int *, double *, int *);
 
__ESVI scrft(int *,void *, int *, void *, int *, int *, int *, int *, float *,
             double *, int *, double *, int *, double *, int *);
 
__ESVI dcrft(int *,void *, int *, void *, int *, int *, int *, int *, double *,
             double *, int *, double *, int *);
 
__ESVI dcft2(int *, void *, int *, int *, void *, int *, int *, int *,
              int *, int *, double *, double *, int *, double *, int *);
 
__ESVI dcft3(void *, int *, int *, void *, int *, int *, int *, int *,
              int *, int *, double *, double *, int *);
 
__ESVI drcft2(int *, void *, int *, void *, int *, int *, int *, int *,
            double *, double *, int *, double *, int *);
 
__ESVI dcrft2(int *, void *, int *, void *, int *, int *, int *, int *,
            double *, double *, int *, double *, int *);
 
__ESVI drcft3(void *, int *, int *, void *, int *, int *, int *, int *,
              int *, int *, double *, double *, int *);
 
__ESVI dcrft3(void *, int *, int *, void *, int *, int *, int *, int *,
               int *, int *, double *, double *, int *);
 
__ESVI scosft(int *, void *, int *, int *, void *, int *, int *, int *, int *,
              float *, double *, int *, double *, int *);
 
__ESVI scosf(int *, void *, int *, int *, void *, int *, int *, int *, int *,
              float *, double *, int *, double *, int *);
 
__ESVI dcosf(int *, void *, int *, int *, void *, int *, int *, int *, int *,
             double *, double *, int *, double *, int *);
 
__ESVI ssinf(int *, void *, int *, int *, void *, int *, int *, int *, int *,
              float *, double *, int *, double *, int *);
 
__ESVI dsinf(int *, void *, int *, int *, void *, int *, int *, int *, int *,
             double *, double *, int *, double *, int *);
 
__ESVI scft2(int *, void *, int *, int *, void *, int *, int *, int *, int *,
             int *, float *, double *, int *, double *, int *);
 
__ESVI srcft2(int *,void *, int *, void *, int *, int *, int *, int *, float *,
              double *, int *, double *, int *, double *, int *);
 
__ESVI scrft2(int *,void *, int *, void *, int *, int *, int *, int *, float *,
              double *, int *, double *, int *, double *, int *);
 
__ESVI scft3(void *, int *, int *, void *, int *, int *, int *, int *, int *,
             int *, float *, double *, int *);
 
__ESVI srcft3(void *, int *, int *, void *, int *, int *, int *, int *, int *,
              int *, float *, double *, int *);
 
__ESVI scrft3(void *, int *, int *, void *, int *, int *, int *, int *, int *,
              int *, float *, double *, int *);
 
/*  Convolutions/Correlations Subroutines  */
 
__ESVI scon(int *, float *, int *, void *, int *, int *, void *, int *, int *,
            int *, int *, int *, int *, int *, double *, int *, double *,
            int *);
__ESVI scor(int *, float *, int *, void *, int *, int *, void *, int *, int *,
            int *, int *, int *, int *, int *, double *, int *, double *,
            int *);
 
void   scond(float *, int *, void *, int *, void *, int *, int *, int *,
               int *, int *);
void   scord(float *, int *, void *, int *, void *, int *, int *, int *,
               int *, int *);
 
__ESVI sconf(int *, float *, int *, void *, int *, int *, void *, int *,
             int *, int *, int *, int *, int *, int *, double *, int *,
             double *, int *);
__ESVI scorf(int *, float *, int *, void *, int *, int *, void *, int *,
             int *, int *, int *, int *, int *, int *, double *, int *,
             double *, int *);
 
void   sdcon( float *, int *, void *, int *, void *, int *, int *, int *,
              int *, int *, int *);
void   ddcon(double *, int *, void *, int *, void *, int *, int *, int *,
              int *, int *, int *);
void   sdcor( float *, int *, void *, int *, void *, int *, int *, int *,
              int *, int *, int *);
void   ddcor(double *, int *, void *, int *, void *, int *, int *, int *,
              int *, int *, int *);
 
__ESVI sacor(int *, void *, int *, int *, void *, int *, int *, int *,
             int *, int *, double *, int *, double *, int *);
 
__ESVI sacorf(int *, void *, int *, int *, void *, int *, int *, int *,
              int *, int *, double *, int *, double *, int *);
 
/*  Related Computations Subroutines  */
 
void   spoly( float *, int *, int *,  float *, int *,  float *, int *,
              int *);
void   dpoly(double *, int *, int *, double *, int *, double *, int *,
              int *);
 
void   sizc( float *, int *, int *, int *, int *);
void   dizc(double *, int *, int *, int *, int *);
 
void   strec( float *,  float *, int *,  float *, int *,  float *,
              int *, int *, int *);
void   dtrec(double *, double *, int *, double *, int *, double *, int *,
             int *, int *);
 
__ESVI sqint( float *, float *, float *, float *, int *, int *, float *,
              int *, float *, int *, int *);
__ESVI dqint(double *, double *, double *, double *, int *, int *, double *,
             int *, double *, int *, int *);
 
__ESVI swlev( float *, int *,  float *, int *,  float *, int *, int *,
              double *, int *);
__ESVI dwlev(double *, int *, double *, int *, double *, int *, int *,
             double *, int *);
 
/*  Sorting and Searching Subroutines  */
 
void   isort(   int *, int *, int *);
void   ssort( float *, int *, int *);
void   dsort(double *, int *, int *);
 
void   isortx(   int *, int *, int *, int *);
void   ssortx( float *, int *, int *, int *);
void   dsortx(double *, int *, int *, int *);
 
void   ibsrch(   int *, int *, int *,    int *, int *, int *, int *,
                 int *, int *);
void   sbsrch( float *, int *, int *,  float *, int *, int *, int *,
                 int *, int *);
void   dbsrch(double *, int *, int *, double *, int *, int *, int *,
                 int *, int *);
 
void   issrch(   int *, int *, int *,    int *, int *, int *, int *,
                 int *);
void   sssrch( float *, int *, int *,  float *, int *, int *, int *,
                 int *);
void   dssrch(double *, int *, int *, double *, int *, int *, int *,
                 int *);
 
void   isorts(   int *, int *, int *, int *,    int *, int *);
void   ssorts( float *, int *, int *, int *,  float *, int *);
void   dsorts(double *, int *, int *, int *, double *, int *);
 
/*  Interpolation Subroutines  */
 
void   spint( float *,  float *, int *,  float *, int *,  float *,  float *,
             int *);
void   dpint(double *, double *, int *, double *, int *, double *, double *,
             int *);
 
__ESVI stpint( float *,  float *, int *, int *,  float *,  float *, int *,
               float *, int *);
__ESVI dtpint(double *, double *, int *, int *, double *, double *, int *,
              double *, int *);
 
void   scsint( float *,  float *, void *, int *, int *,  float *,
               float *, int *);
void   dcsint(double *, double *, void *, int *, int *, double *,
              double *, int *);
 
__ESVI scsin2( float *,  float *, void *, int *, int *, int *,  float *,
               float *, int *, int *, void *, int *,  float *, int *);
__ESVI dcsin2(double *, double *, void *, int *, int *, int *, double *,
              double *, int *, int *, void *, int *, double *, int *);
 
 
/*  Numerical Quadrature Subroutines  */
 
void   sptnq( float *,  float *, int *,  float *,  float *);
void   dptnq(double *, double *, int *, double *, double *);
 
float  sglnq(void *(),  float *,  float *, int *);
double dglnq(void *(), double *, double *, int *);
 
float  sglnq2(void *(), float *, float *, int *, float *, float *, int *,
              void *, int *);
double dglnq2(void *(), double *, double *, int *, double *, double *, int *,
              void *, int *);
 
float  sglgq(void *(),  float *,  float *, int *);
double dglgq(void *(), double *, double *, int *);
 
float  sgraq(void *(),  float *,  float *, int *);
double dgraq(void *(), double *, double *, int *);
 
float  sghmq(void *(),  float *,  float *, int *);
double dghmq(void *(), double *, double *, int *);
 
/*  Random Number Generation Subroutines  */
 
void   surand(double *, int *,  float *);
void   durand(double *, int *, double *);
 
__ESVI snrand(double *, int *,  float *,  float *, int *);
__ESVI dnrand(double *, int *, double *, double *, int *);
 
void surxor(int *, int *,  float *,  float *);
void durxor(int *, int *, double *, double *);
 
/*  Parallel Processing Subroutines  */
 
void   dgemlp(void *,int *, char *, void *, int *, char *, void *, int *, int*,
              int *, int *);
 
__ESVI dgefp(void *, int *, int *, int *, double *, int *);
 
__ESVI dppfp(double *, int *, double *, int *);
 
__ESVI dgkfsp(int *, double *, int *, int *, double *, int *, int *, int *,
              double *, double *, int *, void *, int *, int *);
 
__ESVI dskfsp(int *, double *, int *, int *, int *, double *, double *,
              int *, void *, int *, int *);
 
__ESVI scftp(int *, void *, int *, int *, void *, int *, int *, int *,
             int *, int *, float *, double *, int *, double *, int *);
 
__ESVI scft2p(int *, void *, int *, int *, void *, int *, int *, int *, int *,
              int *, float *, double *, int *, double *, int *);
 
__ESVI scft3p(void *, int *, int *, void *, int *, int *, int *, int *, int *,
              int *, float *, double *, int *);
 
 
/*  Utility Subroutines  */
 
void   einfo(int *, int *, int *);
void   dsrsm(int *, double *, int *, int *, int *, int *, void *, void *,
             int *);
void   stride(int *, int *, int *, char *, int *);
void   ivsset(int *);
 
__ESVI dgktrn(int *, double *, int *, int *, double *, int *, int *, int *,
            double *, int *);
__ESVI dsktrn(int *, double *, int *, int *, int *, double *, int *);
void ievops(int *);
 
 
 
#define isamax(i,x,j)    (__esvi1=i,__esvi2=j,isamax(&__esvi1,x,&__esvi2))
#define idamax(i,x,j)    (__esvi1=i,__esvi2=j,idamax(&__esvi1,x,&__esvi2))
#define icamax(i,x,j)    (__esvi1=i,__esvi2=j,icamax(&__esvi1,x,&__esvi2))
#define izamax(i,x,j)    (__esvi1=i,__esvi2=j,izamax(&__esvi1,x,&__esvi2))
 
#define isamin(i,x,j)    (__esvi1=i,__esvi2=j,isamin(&__esvi1,x,&__esvi2))
#define idamin(i,x,j)    (__esvi1=i,__esvi2=j,idamin(&__esvi1,x,&__esvi2))
 
#define ismax(i,x,j)     (__esvi1=i,__esvi2=j,ismax(&__esvi1,x,&__esvi2))
#define idmax(i,x,j)     (__esvi1=i,__esvi2=j,idmax(&__esvi1,x,&__esvi2))
 
#define ismin(i,x,j)     (__esvi1=i,__esvi2=j,ismin(&__esvi1,x,&__esvi2))
#define idmin(i,x,j)     (__esvi1=i,__esvi2=j,idmin(&__esvi1,x,&__esvi2))
 
#define sasum(i,x,j)     (__esvi1=i,__esvi2=j,sasum(&__esvi1,x,&__esvi2))
#define dasum(i,x,j)     (__esvi1=i,__esvi2=j,dasum(&__esvi1,x,&__esvi2))
#define scasum(i,x,j)    (__esvi1=i,__esvi2=j,scasum(&__esvi1,x,&__esvi2))
#define dzasum(i,x,j)    (__esvi1=i,__esvi2=j,dzasum(&__esvi1,x,&__esvi2))
 
#define saxpy(i1,f1,x1,i2,x2,i3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                         __esvf1=f1,saxpy(&__esvi1,&__esvf1, \
                                                  x1,&__esvi2,x2,&__esvi3))
 
#define daxpy(i1,d1,x1,i2,x2,i3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                         __esvd1=d1,daxpy(&__esvi1,&__esvd1, \
                                                  x1,&__esvi2,x2,&__esvi3))
 
#define caxpy(i1,c1,x1,i2,x2,i3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                         __esvc1=c1,caxpy(&__esvi1,&__esvc1, \
                                                  x1,&__esvi2,x2,&__esvi3))
 
#define zaxpy(i1,dc1,x1,i2,x2,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                      __esvdc1=dc1,zaxpy(&__esvi1,&__esvdc1, \
                                                   x1,&__esvi2,x2,&__esvi3))
 
 
#define scopy(i1,x1,i2,x2,i3)              (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   scopy(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define dcopy(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   dcopy(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define ccopy(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   ccopy(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define zcopy(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   zcopy(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define sdot(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   sdot(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define ddot(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   ddot(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define cdotu(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   esvcdtu(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                        &__esvctmp),__esvctmp)
 
#define zdotu(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   esvzdtu(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                        &__esvdtmp),__esvdtmp)
 
#define cdotc(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   esvcdtc(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                        &__esvctmp),__esvctmp)
 
#define zdotc(i1,x1,i2,x2,i3)             (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   esvzdtc(&__esvi1,x1,&__esvi2,x2,&__esvi3,  \
                                        &__esvdtmp),__esvdtmp)
 
#define snaxpy(i1,i2,x1,i3,x2,i4,i5,x3,i6,i7)        (__esvi1=i1,__esvi2=i2, \
                      __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvi6=i6,__esvi7=i7, \
                           snaxpy(&__esvi1,&__esvi2,x1,&__esvi3,x2,&__esvi4, \
                                             &__esvi5,x3,&__esvi6,&__esvi7))
 
#define dnaxpy(i1,i2,x1,i3,x2,i4,i5,x3,i6,i7)        (__esvi1=i1,__esvi2=i2, \
                      __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvi6=i6,__esvi7=i7, \
                           dnaxpy(&__esvi1,&__esvi2,x1,&__esvi3,x2,&__esvi4, \
                                             &__esvi5,x3,&__esvi6,&__esvi7))
 
 
#define sndot(i1,i2,x1,i3,i4,x2,i5,i6,x3,i7,i8)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                            __esvi6=i6,__esvi7=i7,__esvi8=i8, \
                            sndot(&__esvi1,&__esvi2,x1,&__esvi3,&__esvi4,x2, \
                                    &__esvi5,&__esvi6,x3,&__esvi7,&__esvi8))
 
#define dndot(i1,i2,x1,i3,i4,x2,i5,i6,x3,i7,i8)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                           __esvi6=i6,__esvi7=i7,__esvi8=i8,  \
                            dndot(&__esvi1,&__esvi2,x1,&__esvi3,&__esvi4,x2, \
                                    &__esvi5,&__esvi6,x3,&__esvi7,&__esvi8))
 
#define snrm2(i,x,j)       (__esvi1=i,__esvi2=j,snrm2(&__esvi1,x,&__esvi2))
#define dnrm2(i,x,j)       (__esvi1=i,__esvi2=j,dnrm2(&__esvi1,x,&__esvi2))
#define scnrm2(i,x,j)     (__esvi1=i,__esvi2=j,scnrm2(&__esvi1,x,&__esvi2))
#define dznrm2(i,x,j)     (__esvi1=i,__esvi2=j,dznrm2(&__esvi1,x,&__esvi2))
 
 
#define snorm2(i,x,j)  (__esvi1=i,__esvi2=j,snorm2(&__esvi1,x,&__esvi2))
#define dnorm2(i,x,j)  (__esvi1=i,__esvi2=j,dnorm2(&__esvi1,x,&__esvi2))
#define cnorm2(i,x,j)  (__esvi1=i,__esvi2=j,cnorm2(&__esvi1,x,&__esvi2))
#define znorm2(i,x,j)  (__esvi1=i,__esvi2=j,znorm2(&__esvi1,x,&__esvi2))
 
#define srotg(x1,x2,x3,x4)            srotg(x1,x2,x3,x4)
#define drotg(x1,x2,x3,x4)            drotg(x1,x2,x3,x4)
#define crotg(x1,x2,x3,x4)            crotg(x1,x2,x3,x4)
#define zrotg(x1,x2,x3,x4)            zrotg(x1,x2,x3,x4)
 
#define srot(i1,x1,i2,x2,i3,f1,f2)      (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvf1=f1,__esvf2=f2,srot(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvf1,&__esvf2))
 
#define drot(i1,x1,i2,x2,i3,d1,d2)      (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvd1=d1,__esvd2=d2,drot(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvd1,&__esvd2))
 
 
#define crot(i1,x1,i2,x2,i3,f1,c2)      (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvf1=f1,__esvc2=c2,crot(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvf1,&__esvc2))
 
#define zrot(i1,x1,i2,x2,i3,d1,dc2)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                  __esvd1=d1,__esvdc2=dc2,zrot(&__esvi1,x1, \
                                 &__esvi2,x2,&__esvi3,&__esvd1,&__esvdc2))
 
#define csrot(i1,x1,i2,x2,i3,f1,f2)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   __esvf1=f1,__esvf2=f2,csrot(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvf1,&__esvf2))
 
#define zdrot(i1,x1,i2,x2,i3,d1,d2)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   __esvd1=d1,__esvd2=d2,zdrot(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvd1,&__esvd2))
 
#define sscal(i1,f1,x1,i2)              (__esvi1=i1,__esvi2=i2,__esvf1=f1, \
                                     sscal(&__esvi1,&__esvf1,x1,&__esvi2))
 
#define dscal(i1,d1,x1,i2)             (__esvi1=i1,__esvi2=i2,__esvd1=d1, \
                                    dscal(&__esvi1,&__esvd1,x1,&__esvi2))
 
#define cscal(i1,c1,x1,i2)             (__esvi1=i1,__esvi2=i2,__esvc1=c1, \
                                    cscal(&__esvi1,&__esvc1,x1,&__esvi2))
 
#define zscal(i1,dc1,x1,i2)          (__esvi1=i1,__esvi2=i2,__esvdc1=dc1, \
                                   zscal(&__esvi1,&__esvdc1,x1,&__esvi2))
 
#define csscal(i1,f1,x1,i2)            (__esvi1=i1,__esvi2=i2,__esvf1=f1, \
                                   csscal(&__esvi1,&__esvf1,x1,&__esvi2))
 
#define zdscal(i1,d1,x1,i2)            (__esvi1=i1,__esvi2=i2,__esvd1=d1, \
                                   zdscal(&__esvi1,&__esvd1,x1,&__esvi2))
 
#define sswap(i1,x1,i2,x2,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                 sswap(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define dswap(i1,x1,i2,x2,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                 dswap(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define cswap(i1,x1,i2,x2,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                 cswap(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define zswap(i1,x1,i2,x2,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                 zswap(&__esvi1,x1,&__esvi2,x2,&__esvi3))
 
 
 
 
#define syax(i1,f1,x1,i2,x2,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                     __esvf1=f1,syax(&__esvi1,&__esvf1, \
                                             x1,&__esvi2,x2,&__esvi3))
 
#define dyax(i1,d1,x1,i2,x2,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                     __esvd1=d1,dyax(&__esvi1,&__esvd1, \
                                             x1,&__esvi2,x2,&__esvi3))
 
#define cyax(i1,c1,x1,i2,x2,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                     __esvc1=c1,cyax(&__esvi1,&__esvc1, \
                                             x1,&__esvi2,x2,&__esvi3))
 
#define zyax(i1,dc1,x1,i2,x2,i3)    (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                  __esvdc1=dc1,zyax(&__esvi1,&__esvdc1, \
                                             x1,&__esvi2,x2,&__esvi3))
 
#define csyax(i1,f1,x1,i2,x2,i3)    (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvf1=f1,csyax(&__esvi1,&__esvf1, \
                                             x1,&__esvi2,x2,&__esvi3))
 
#define zdyax(i1,d1,x1,i2,x2,i3)    (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvd1=d1,zdyax(&__esvi1,&__esvd1, \
                                             x1,&__esvi2,x2,&__esvi3))
 
 
 
#define szaxpy(i1,f1,x1,i2,x2,i3,x3,i4) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                            __esvi4=i4,__esvf1=f1,szaxpy(&__esvi1,&__esvf1, \
                                     x1,&__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define dzaxpy(i1,d1,x1,i2,x2,i3,x3,i4) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                            __esvi4=i4,__esvd1=d1,dzaxpy(&__esvi1,&__esvd1, \
                                     x1,&__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define czaxpy(i1,c1,x1,i2,x2,i3,x3,i4) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                            __esvi4=i4,__esvc1=c1,czaxpy(&__esvi1,&__esvc1, \
                                     x1,&__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define zzaxpy(i1,dc1,x1,i2,x2,i3,x3,i4) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                                    __esvi4=i4,__esvdc1=dc1, \
                                              zzaxpy(&__esvi1,&__esvdc1,x1, \
                                         &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
 
#define svea(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,svea(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define dvea(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,dvea(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define cvea(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,cvea(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define zvea(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,zvea(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define sves(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,sves(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define dves(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,dves(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define cves(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,cves(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define zves(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,zves(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define svem(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,svem(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define dvem(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,dvem(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define cvem(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,cvem(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
#define zvem(i1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                                      __esvi4=i4,zvem(&__esvi1,x1,    \
                                      &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
 
 
 
 
/*  Sparse Vector-Scalar Subroutines  */
 
#define ssctr(i1,x1,x2,x3)        (__esvi1=i1,ssctr(&__esvi1,x1,x2,x3))
#define dsctr(i1,x1,x2,x3)        (__esvi1=i1,dsctr(&__esvi1,x1,x2,x3))
#define csctr(i1,x1,x2,x3)   (__esvi1=i1,csctr(&__esvi1,x1,x2,x3))
#define zsctr(i1,x1,x2,x3)   (__esvi1=i1,zsctr(&__esvi1,x1,x2,x3))
 
#define sgthr(i1,x1,x2,x3)        (__esvi1=i1,sgthr(&__esvi1,x1,x2,x3))
#define dgthr(i1,x1,x2,x3)        (__esvi1=i1,dgthr(&__esvi1,x1,x2,x3))
#define cgthr(i1,x1,x2,x3)   (__esvi1=i1,cgthr(&__esvi1,x1,x2,x3))
#define zgthr(i1,x1,x2,x3)   (__esvi1=i1,zgthr(&__esvi1,x1,x2,x3))
 
#define sgthrz(i1,x1,x2,x3)      (__esvi1=i1,sgthrz(&__esvi1,x1,x2,x3))
#define dgthrz(i1,x1,x2,x3)      (__esvi1=i1,dgthrz(&__esvi1,x1,x2,x3))
#define cgthrz(i1,x1,x2,x3)  (__esvi1=i1,cgthrz(&__esvi1,x1,x2,x3))
#define zgthrz(i1,x1,x2,x3)  (__esvi1=i1,zgthrz(&__esvi1,x1,x2,x3))
 
#define saxpyi(i1,f1,x1,x2,x3)  (__esvi1=i1,__esvf1=f1,saxpyi(&__esvi1,  \
                                                     &__esvf1,x1,x2,x3))
#define daxpyi(i1,d1,x1,x2,x3)  (__esvi1=i1,__esvd1=d1,daxpyi(&__esvi1,  \
                                                     &__esvd1,x1,x2,x3))
#define caxpyi(i1,c1,x1,x2,x3)  (__esvi1=i1,__esvc1=c1,caxpyi(&__esvi1, \
                                &__esvc1,x1,x2,x3))
#define zaxpyi(i1,dc1,x1,x2,x3)  (__esvi1=i1,__esvdc1=dc1,zaxpyi(&__esvi1, \
                                &__esvdc1,x1,x2,x3))
 
#define sdoti(i1,x1,x2,x3)        (__esvi1=i1,sdoti(&__esvi1,x1,x2,x3))
#define ddoti(i1,x1,x2,x3)        (__esvi1=i1,ddoti(&__esvi1,x1,x2,x3))
 
#define cdotui(i1,x1,x2,x3)  (__esvi1=i1,esvcdtui(&__esvi1,x1,x2,x3,  \
                              &__esvctmp),__esvctmp)
#define zdotui(i1,x1,x2,x3)  (__esvi1=i1,esvzdtui(&__esvi1,x1,x2,x3,  \
                              &__esvdtmp),__esvdtmp)
 
#define cdotci(i1,x1,x2,x3)  (__esvi1=i1,esvcdtci(&__esvi1,x1,x2,x3,  \
                              &__esvctmp),__esvctmp)
#define zdotci(i1,x1,x2,x3)  (__esvi1=i1,esvzdtci(&__esvi1,x1,x2,x3,  \
                              &__esvdtmp),__esvdtmp)
 
 
/* Matrix-Vector linear algebra subroutines  */
 
#define sgemv(ch1,i1,i2,f1,x1,i3,x2,i4,f2,x3,i5)  (__esvi1=i1,__esvi2=i2, \
                   __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvf1=f1,__esvf2=f2, \
                     sgemv(ch1,&__esvi1,&__esvi2,&__esvf1,x1,&__esvi3,x2, \
                                          &__esvi4,&__esvf2,x3,&__esvi5))
 
#define dgemv(ch1,i1,i2,d1,x1,i3,x2,i4,d2,x3,i5)  (__esvi1=i1,__esvi2=i2, \
                   __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvd1=d1,__esvd2=d2, \
                     dgemv(ch1,&__esvi1,&__esvi2,&__esvd1,x1,&__esvi3,x2, \
                                          &__esvi4,&__esvd2,x3,&__esvi5))
 
#define cgemv(ch1,i1,i2,c1,x1,i3,x2,i4,c2,x3,i5)  (__esvi1=i1,__esvi2=i2, \
                   __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvc1=c1,__esvc2=c2, \
                     cgemv(ch1,&__esvi1,&__esvi2,&__esvc1,x1,&__esvi3,x2, \
                                          &__esvi4,&__esvc2,x3,&__esvi5))
 
#define zgemv(ch1,i1,i2,dc1,x1,i3,x2,i4,dc2,x3,i5)  (__esvi1=i1,__esvi2=i2, \
                 __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvdc1=dc1,__esvdc2=dc2, \
                      zgemv(ch1,&__esvi1,&__esvi2,&__esvdc1,x1,&__esvi3,x2, \
                                           &__esvi4,&__esvdc2,x3,&__esvi5))
 
 
#define sgemx(i1,i2,f1,x1,i3,x2,i4,x3,i5)          (__esvi1=i1,__esvi2=i2, \
                              __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvf1=f1,  \
                                      sgemx(&__esvi1,&__esvi2,&__esvf1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
 
#define dgemx(i1,i2,d1,x1,i3,x2,i4,x3,i5)          (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvd1=d1, \
                                      dgemx(&__esvi1,&__esvi2,&__esvd1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
 
#define sgemtx(i1,i2,f1,x1,i3,x2,i4,x3,i5)  (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvf1=f1, \
                                     sgemtx(&__esvi1,&__esvi2,&__esvf1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
#define dgemtx(i1,i2,d1,x1,i3,x2,i4,x3,i5)         (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvd1=d1, \
                                     dgemtx(&__esvi1,&__esvi2,&__esvd1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
 
#define sger1(i1,i2,f1,x1,i3,x2,i4,x3,i5)          (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvf1=f1, \
                                      sger1(&__esvi1,&__esvi2,&__esvf1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
#define dger1(i1,i2,d1,x1,i3,x2,i4,x3,i5)          (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvd1=d1, \
                                      dger1(&__esvi1,&__esvi2,&__esvd1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
#define cgeru(i1,i2,c1,x1,i3,x2,i4,x3,i5)          (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvc1=c1, \
                                      cgeru(&__esvi1,&__esvi2,&__esvc1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
#define zgeru(i1,i2,dc1,x1,i3,x2,i4,x3,i5)         (__esvi1=i1,__esvi2=i2, \
                             __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvdc1=dc1, \
                                     zgeru(&__esvi1,&__esvi2,&__esvdc1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
#define cgerc(i1,i2,c1,x1,i3,x2,i4,x3,i5)          (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvc1=c1, \
                                      cgerc(&__esvi1,&__esvi2,&__esvc1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
 
#define zgerc(i1,i2,dc1,x1,i3,x2,i4,x3,i5)         (__esvi1=i1,__esvi2=i2, \
                             __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvdc1=dc1, \
                                     zgerc(&__esvi1,&__esvi2,&__esvdc1,x1, \
                                        &__esvi3,x2,&__esvi4,x3,&__esvi5))
 
 
#define sslmx(i1,f1,x1,x2,i2,x3,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvf1=f1,sslmx(&__esvi1,&__esvf1,x1, \
                                                 x2,&__esvi2,x3,&__esvi3))
 
#define dslmx(i1,d1,x1,x2,i2,x3,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvd1=d1,dslmx(&__esvi1,&__esvd1,x1, \
                                                 x2,&__esvi2,x3,&__esvi3))
 
#define sslr1(i1,f1,x1,i2,x2)           (__esvi1=i1,__esvi2=i2,__esvf1=f1, \
                                 sslr1(&__esvi1,&__esvf1,x1,&__esvi2,x2))
 
#define dslr1(i1,d1,x1,i2,x2)           (__esvi1=i1,__esvi2=i2,__esvd1=d1, \
                                 dslr1(&__esvi1,&__esvd1,x1,&__esvi2,x2))
 
 
#define sslr2(i1,f1,x1,i2,x2,i3,x3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                       __esvf1=f1,sslr2(&__esvi1,&__esvf1, \
                                              x1,&__esvi2,x2,&__esvi3,x3))
 
#define dslr2(i1,d1,x1,i2,x2,i3,x3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                       __esvd1=d1,dslr2(&__esvi1,&__esvd1, \
                                              x1,&__esvi2,x2,&__esvi3,x3))
#define strmv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1, __esvi2=i2,\
                                          __esvi3=i3, strmv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
#define dtrmv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1, __esvi2=i2,\
                                          __esvi3=i3, dtrmv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
#define ctrmv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1, __esvi2=i2,\
                                          __esvi3=i3, ctrmv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
#define ztrmv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1, __esvi2=i2,\
                                          __esvi3=i3, ztrmv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define sspmv(ch1,i1,f1,x1,x2,i2,f2,x3,i3) (__esvi1=i1,__esvi2=i2,     \
                                    __esvi3=i3,__esvf1=f1,__esvf2=f2,  \
                                    sspmv(ch1,&__esvi1,&__esvf1,x1,x2,\
                                    &__esvi2,&__esvf2,x3,&__esvi3))
 
#define dspmv(ch1,i1,d1,x1,x2,i2,d2,x3,i3) (__esvi1=i1,__esvi2=i2,     \
                                    __esvi3=i3,__esvd1=d1,__esvd2=d2,  \
                                    dspmv(ch1,&__esvi1,&__esvd1,x1,x2,\
                                    &__esvi2,&__esvd2,x3,&__esvi3))
 
#define chpmv(ch1,i1,c1,x1,x2,i2,c2,x3,i3) (__esvi1=i1,__esvi2=i2,     \
                                    __esvi3=i3,__esvc1=c1,__esvc2=c2,  \
                                    chpmv(ch1,&__esvi1,&__esvc1,x1,x2,\
                                    &__esvi2,&__esvc2,x3,&__esvi3))
 
#define zhpmv(ch1,i1,dc1,x1,x2,i2,dc2,x3,i3)  (__esvi1=i1,__esvi2=i2, \
                                __esvi3=i3,__esvdc1=dc1,__esvdc2=dc2, \
                                  zhpmv(ch1,&__esvi1,&__esvdc1,x1,x2,\
                                    &__esvi2,&__esvdc2,x3,&__esvi3))
 
#define ssymv(ch1,i1,f1,x1,i2,x2,i3,f2,x3,i4) (__esvi1=i1,__esvi2=i2,  \
                                    __esvi3=i3,__esvi4=i4,__esvf1=f1,  \
                                    __esvf2=f2,ssymv(ch1,&__esvi1,    \
                                    &__esvf1,x1,&__esvi2,x2,           \
                                    &__esvi3,&__esvf2,x3,&__esvi4))
 
#define dsymv(ch1,i1,d1,x1,i2,x2,i3,d2,x3,i4) (__esvi1=i1,__esvi2=i2,  \
                                    __esvi3=i3,__esvi4=i4,__esvd1=d1,  \
                                    __esvd2=d2,dsymv(ch1,&__esvi1,    \
                                    &__esvd1,x1,&__esvi2,x2,           \
                                    &__esvi3,&__esvd2,x3,&__esvi4))
 
#define chemv(ch1,i1,c1,x1,i2,x2,i3,c2,x3,i4) (__esvi1=i1,__esvi2=i2,  \
                                    __esvi3=i3,__esvi4=i4,__esvc1=c1,  \
                                    __esvc2=c2,chemv(ch1,&__esvi1,    \
                                    &__esvc1,x1,&__esvi2,x2,           \
                                    &__esvi3,&__esvc2,x3,&__esvi4))
 
 
#define zhemv(ch1,i1,dc1,x1,i2,x2,i3,dc2,x3,i4) (__esvi1=i1,__esvi2=i2,\
                                    __esvi3=i3,__esvi4=i4,__esvdc1=dc1,\
                                    __esvdc2=dc2,zhemv(ch1,&__esvi1,  \
                                    &__esvdc1,x1,&__esvi2,x2,          \
                                    &__esvi3,&__esvdc2,x3,&__esvi4))
 
 
#define sspr(ch1,i1,f1,x1,i2,x2)  (__esvi1=i1,__esvi2=i2,__esvf1=f1,  \
                                   sspr(ch1,&__esvi1,&__esvf1,x1,    \
                                   &__esvi2,x2))
 
#define dspr(ch1,i1,d1,x1,i2,x2)  (__esvi1=i1,__esvi2=i2,__esvd1=d1,  \
                                   dspr(ch1,&__esvi1,&__esvd1,x1,    \
                                   &__esvi2,x2))
 
 
#define chpr(ch1,i1,f1,x1,i2,x2)  (__esvi1=i1,__esvi2=i2,__esvf1=f1,  \
                                   chpr(ch1,&__esvi1,&__esvf1,x1,    \
                                   &__esvi2,x2))
 
#define zhpr(ch1,i1,d1,x1,i2,x2)  (__esvi1=i1,__esvi2=i2,__esvd1=d1,\
                                   zhpr(ch1,&__esvi1,&__esvd1,x1,    \
                                   &__esvi2,x2))
 
#define ssyr(ch1,i1,f1,x1,i2,x2,i3) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                             __esvf1=f1,ssyr(ch1,&__esvi1,&__esvf1,x1,\
                                   &__esvi2,x2,&__esvi3))
 
#define dsyr(ch1,i1,d1,x1,i2,x2,i3) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                             __esvd1=d1,dsyr(ch1,&__esvi1,&__esvd1,x1,\
                                   &__esvi2,x2,&__esvi3))
 
#define cher(ch1,i1,f1,x1,i2,x2,i3) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                             __esvf1=f1,cher(ch1,&__esvi1,&__esvf1,x1,\
                                   &__esvi2,x2,&__esvi3))
 
#define zher(ch1,i1,d1,x1,i2,x2,i3) (__esvi1=i1,__esvi2=i2,__esvi3=i3,\
                          __esvd1=d1,zher(ch1,&__esvi1,&__esvd1,x1,\
                                   &__esvi2,x2,&__esvi3))
 
#define sspr2(ch1,i1,f1,x1,i2,x2,i3,x3)       (__esvi1=i1,__esvi2=i2, \
                                               __esvi3=i3,__esvf1=f1, \
                                      sspr2(ch1,&__esvi1,&__esvf1,x1,\
                                       &__esvi2,x2,&__esvi3,x3))
 
#define dspr2(ch1,i1,d1,x1,i2,x2,i3,x3)       (__esvi1=i1,__esvi2=i2, \
                                               __esvi3=i3,__esvd1=d1, \
                                      dspr2(ch1,&__esvi1,&__esvd1,x1,\
                                       &__esvi2,x2,&__esvi3,x3))
 
#define chpr2(ch1,i1,c1,x1,i2,x2,i3,x3)       (__esvi1=i1,__esvi2=i2, \
                                               __esvi3=i3,__esvc1=c1, \
                                      chpr2(ch1,&__esvi1,&__esvc1,x1,\
                                       &__esvi2,x2,&__esvi3,x3))
 
#define zhpr2(ch1,i1,dc1,x1,i2,x2,i3,x3)      (__esvi1=i1,__esvi2=i2, \
                                             __esvi3=i3,__esvdc1=dc1, \
                                     zhpr2(ch1,&__esvi1,&__esvdc1,x1,\
                                      &__esvi2,x2,&__esvi3,x3))
 
#define ssyr2(ch1,i1,f1,x1,i2,x2,i3,x3,i4)    (__esvi1=i1,__esvi2=i2, \
                                    __esvi3=i3,__esvi4=i4,__esvf1=f1, \
                                    ssyr2(ch1,&__esvi1,&__esvf1,x1,  \
                                  &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
 
#define dsyr2(ch1,i1,d1,x1,i2,x2,i3,x3,i4)    (__esvi1=i1,__esvi2=i2, \
                                    __esvi3=i3,__esvi4=i4,__esvd1=d1, \
                                    dsyr2(ch1,&__esvi1,&__esvd1,x1,  \
                                  &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define cher2(ch1,i1,c1,x1,i2,x2,i3,x3,i4)    (__esvi1=i1,__esvi2=i2, \
                                    __esvi3=i3,__esvi4=i4,__esvc1=c1, \
                                    cher2(ch1,&__esvi1,&__esvc1,x1,  \
                                  &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define zher2(ch1,i1,dc1,x1,i2,x2,i3,x3,i4)   (__esvi1=i1,__esvi2=i2, \
                                  __esvi3=i3,__esvi4=i4,__esvdc1=dc1, \
                                    zher2(ch1,&__esvi1,&__esvdc1,x1, \
                                  &__esvi2,x2,&__esvi3,x3,&__esvi4))
 
#define sgbmv(ch1,i1,i2,i3,i4,f1,x1,i5,x2,i6,f2,x3,i7)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                    __esvf1=f1,__esvf2=f2,sgbmv(ch1, \
                                  &__esvi1,&__esvi2,&__esvi3,&__esvi4,\
                                  &__esvf1,x1,&__esvi5,x2,&__esvi6,   \
                                  &__esvf2,x3,&__esvi7))
 
#define dgbmv(ch1,i1,i2,i3,i4,d1,x1,i5,x2,i6,d2,x3,i7)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                    __esvd1=d1,__esvd2=d2,dgbmv(ch1, \
                                  &__esvi1,&__esvi2,&__esvi3,&__esvi4,\
                                  &__esvd1,x1,&__esvi5,x2,&__esvi6,   \
                                  &__esvd2,x3,&__esvi7))
 
#define cgbmv(ch1,i1,i2,i3,i4,c1,x1,i5,x2,i6,c2,x3,i7)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                    __esvc1=c1,__esvc2=c2,cgbmv(ch1, \
                                  &__esvi1,&__esvi2,&__esvi3,&__esvi4,\
                                  &__esvc1,x1,&__esvi5,x2,&__esvi6,   \
                                  &__esvc2,x3,&__esvi7))
 
#define zgbmv(ch1,i1,i2,i3,i4,dc1,x1,i5,x2,i6,dc2,x3,i7) (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                 __esvdc1=dc1,__esvdc2=dc2,zgbmv(ch1,\
                                  &__esvi1,&__esvi2,&__esvi3,&__esvi4,\
                                  &__esvdc1,x1,&__esvi5,x2,&__esvi6,  \
                                  &__esvdc2,x3,&__esvi7))
 
#define ssbmv(ch1,i1,i2,f1,x1,i3,x2,i4,f2,x3,i5)         (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvf1=f1,__esvf2=f2, \
                                        ssbmv(ch1,&__esvi1,&__esvi2, \
                                  &__esvf1,x1,&__esvi3,x2,&__esvi4,   \
                                  &__esvf2,x3,&__esvi5))
 
#define dsbmv(ch1,i1,i2,d1,x1,i3,x2,i4,d2,x3,i5)         (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvd1=d1,__esvd2=d2, \
                                        dsbmv(ch1,&__esvi1,&__esvi2, \
                                  &__esvd1,x1,&__esvi3,x2,&__esvi4,   \
                                  &__esvd2,x3,&__esvi5))
 
#define chbmv(ch1,i1,i2,c1,x1,i3,x2,i4,c2,x3,i5)         (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvc1=c1,__esvc2=c2, \
                                        chbmv(ch1,&__esvi1,&__esvi2, \
                                  &__esvc1,x1,&__esvi3,x2,&__esvi4,   \
                                  &__esvc2,x3,&__esvi5))
 
#define zhbmv(ch1,i1,i2,dc1,x1,i3,x2,i4,dc2,x3,i5)       (__esvi1=i1,  \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4,  \
                                 __esvi5=i5,__esvdc1=dc1,__esvdc2=dc2, \
                                        zhbmv(ch1,&__esvi1,&__esvi2,  \
                                  &__esvdc1,x1,&__esvi3,x2,&__esvi4,   \
                                  &__esvdc2,x3,&__esvi5))
 
#define stbmv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        stbmv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
#define dtbmv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        dtbmv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
#define ctbmv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        ctbmv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
#define ztbmv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        ztbmv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
 
#define stpmv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               stpmv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
#define dtpmv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               dtpmv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
#define ctpmv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               ctpmv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
#define ztpmv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               ztpmv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
 
/*  Sparse Matrix-Vector linear algebra Subroutines  */
 
 
#define dsmmx(i1,i2,x1,x2,i3,x3,x4)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                            dsmmx(&__esvi1,&__esvi2,x1,x2, \
                                                          &__esvi3,x3,x4))
 
#ifdef  __ESVERR
#define dsmtm(i1,i2,x1,x2,i3,x6,x7,x3,x4,i6,x5,x8)  (__esvi1=i1,__esvi2=i2, \
                                                      __esvi3=i3,__esvi6=i6, \
                                             dsmtm(&__esvi1,&__esvi2,x1,x2, \
                                                       &__esvi3,x6,x7,x3,x4, \
                                                           &__esvi6,x5,x8))
#else
#define dsmtm(i1,i2,x1,x2,i3,x6,x7,x3,x4,i6,x5,i7)  (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi6=i6,__esvi7=i7, \
                                            dsmtm(&__esvi1,&__esvi2,x1,x2,  \
                                                       &__esvi3,x6,x7,x3,x4, \
                                                     &__esvi6,x5,&__esvi7))
#endif
 
#define dsdmx(i1,i2,i3,x1,i4,ch1,x2,x3,x4)          (__esvi1=i1,__esvi2=i2, \
                                      __esvi3=i3,__esvi4=i4,dsdmx(&__esvi1, \
                               &__esvi2,&__esvi3,x1,&__esvi4,ch1,x2,x3,x4))
 
 
/*  Matrix Operation Subroutines  */
 
#define sgeadd(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         sgeadd(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
#define dgeadd(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         dgeadd(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
#define cgeadd(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         cgeadd(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
#define zgeadd(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         zgeadd(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
 
 
#define sgesub(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         sgesub(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
#define dgesub(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         dgesub(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
#define cgesub(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         cgesub(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
#define zgesub(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5)      (__esvi1=i1,__esvi2=i2, \
                                            __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         zgesub(x1,&__esvi1,ch1,x2,&__esvi2, \
                                         ch2,x3,&__esvi3,&__esvi4,&__esvi5))
 
#define sgemul(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6)   (__esvi1=i1,__esvi2=i2, \
                                                       __esvi3=i3,__esvi4=i4, \
                                                       __esvi5=i5,__esvi6=i6, \
                                                  sgemul(x1,&__esvi1,ch1,x2, \
                                                    &__esvi2,ch2,x3,&__esvi3, \
                                                &__esvi4,&__esvi5,&__esvi6))
 
#define dgemul(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6)  (__esvi1=i1,__esvi2=i2,  \
                                                       __esvi3=i3,__esvi4=i4, \
                                                       __esvi5=i5,__esvi6=i6, \
                                                  dgemul(x1,&__esvi1,ch1,x2, \
                                                    &__esvi2,ch2,x3,&__esvi3, \
                                                &__esvi4,&__esvi5,&__esvi6))
 
#define cgemul(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6)  (__esvi1=i1,__esvi2=i2,  \
                                                       __esvi3=i3,__esvi4=i4, \
                                                       __esvi5=i5,__esvi6=i6, \
                                                  cgemul(x1,&__esvi1,ch1,x2, \
                                                    &__esvi2,ch2,x3,&__esvi3, \
                                                &__esvi4,&__esvi5,&__esvi6))
 
#define zgemul(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6)  (__esvi1=i1,__esvi2=i2,  \
                                                       __esvi3=i3,__esvi4=i4, \
                                                       __esvi5=i5,__esvi6=i6, \
                                                  zgemul(x1,&__esvi1,ch1,x2, \
                                                    &__esvi2,ch2,x3,&__esvi3, \
                                                &__esvi4,&__esvi5,&__esvi6))
 
#ifdef __ESVERR
#define sgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,x5)       (__esvi1=i1, \
                                           __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 sgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                         &__esvi4,&__esvi5,&__esvi6,x4,x5))
#else
#define sgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,i7)       (__esvi1=i1, \
                                                     __esvi2=i2,__esvi7=i7,  \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 sgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                   &__esvi4,&__esvi5,&__esvi6,x4,&__esvi7))
#endif
 
#ifdef __ESVERR
#define dgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,x5)       (__esvi1=i1, \
                                           __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 dgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                         &__esvi4,&__esvi5,&__esvi6,x4,x5))
#else
#define dgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,i7)       (__esvi1=i1, \
                                                     __esvi2=i2,__esvi7=i7,  \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 dgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                   &__esvi4,&__esvi5,&__esvi6,x4,&__esvi7))
#endif
 
#ifdef __ESVERR
#define cgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,x5)       (__esvi1=i1, \
                                           __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 cgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                         &__esvi4,&__esvi5,&__esvi6,x4,x5))
#else
#define cgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,i7)       (__esvi1=i1, \
                                                     __esvi2=i2,__esvi7=i7,  \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 cgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                   &__esvi4,&__esvi5,&__esvi6,x4,&__esvi7))
#endif
 
#ifdef __ESVERR
#define zgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,x5)       (__esvi1=i1, \
                                           __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 zgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                         &__esvi4,&__esvi5,&__esvi6,x4,x5))
#else
#define zgemms(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6,x4,i7)       (__esvi1=i1, \
                                                     __esvi2=i2,__esvi7=i7,  \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                 zgemms(x1,&__esvi1,ch1,x2, \
                                                   &__esvi2,ch2,x3,&__esvi3, \
                                   &__esvi4,&__esvi5,&__esvi6,x4,&__esvi7))
#endif
 
#define sgemm(ch1,ch2,i1,i2,i3,f1,x1,i4,x2,i5,f2,x3,i6)        (__esvi1=i1, \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                            __esvi6=i6,__esvf1=f1,__esvf2=f2,sgemm(ch1,ch2, \
                                     &__esvi1,&__esvi2,&__esvi3,&__esvf1,x1, \
                                &__esvi4,x2,&__esvi5,&__esvf2,x3,&__esvi6))
 
#define dgemm(ch1,ch2,i1,i2,i3,d1,x1,i4,x2,i5,d2,x3,i6)        (__esvi1=i1, \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                            __esvi6=i6,__esvd1=d1,__esvd2=d2,dgemm(ch1,ch2, \
                                     &__esvi1,&__esvi2,&__esvi3,&__esvd1,x1, \
                                &__esvi4,x2,&__esvi5,&__esvd2,x3,&__esvi6))
 
#define cgemm(ch1,ch2,i1,i2,i3,c1,x1,i4,x2,i5,c2,x3,i6)        (__esvi1=i1, \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                            __esvi6=i6,__esvc1=c1,__esvc2=c2,cgemm(ch1,ch2, \
                                     &__esvi1,&__esvi2,&__esvi3,&__esvc1,x1, \
                                &__esvi4,x2,&__esvi5,&__esvc2,x3,&__esvi6))
 
#define zgemm(ch1,ch2,i1,i2,i3,dc1,x1,i4,x2,i5,dc2,x3,i6)      (__esvi1=i1, \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                        __esvi6=i6,__esvdc1=dc1,__esvdc2=dc2,zgemm(ch1,ch2, \
                                    &__esvi1,&__esvi2,&__esvi3,&__esvdc1,x1, \
                               &__esvi4,x2,&__esvi5,&__esvdc2,x3,&__esvi6))
 
 
#define ssyrk(ch1,ch2,i1,i2,f1,x1,i3,f2,x2,i4)      (__esvi1=i1,__esvi2=i2, \
                                __esvi3=i3,__esvi4=i4,__esvf1=f1,__esvf2=f2, \
                               ssyrk(ch1,ch2,&__esvi1,&__esvi2,&__esvf1,x1, \
                                            &__esvi3,&__esvf2,x2,&__esvi4))
 
#define dsyrk(ch1,ch2,i1,i2,d1,x1,i3,d2,x2,i4)      (__esvi1=i1,__esvi2=i2, \
                                __esvi3=i3,__esvi4=i4,__esvd1=d1,__esvd2=d2, \
                               dsyrk(ch1,ch2,&__esvi1,&__esvi2,&__esvd1,x1, \
                                            &__esvi3,&__esvd2,x2,&__esvi4))
 
#define sgetmi(x1,i1,i2)                  (__esvi1=i1,__esvi2=i2,sgetmi(x1, \
                                                        &__esvi1,&__esvi2))
#define dgetmi(x1,i1,i2)                  (__esvi1=i1,__esvi2=i2,dgetmi(x1, \
                                                        &__esvi1,&__esvi2))
#define cgetmi(x1,i1,i2)                  (__esvi1=i1,__esvi2=i2,cgetmi(x1, \
                                                        &__esvi1,&__esvi2))
#define zgetmi(x1,i1,i2)                  (__esvi1=i1,__esvi2=i2,zgetmi(x1, \
                                                        &__esvi1,&__esvi2))
 
 
#define sgetmo(x1,i1,i2,i3,x2,i4)        (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                             __esvi4=i4,sgetmo(x1,&__esvi1, \
                                            &__esvi2,&__esvi3,x2,&__esvi4))
 
#define dgetmo(x1,i1,i2,i3,x2,i4)        (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                             __esvi4=i4,dgetmo(x1,&__esvi1, \
                                            &__esvi2,&__esvi3,x2,&__esvi4))
 
#define cgetmo(x1,i1,i2,i3,x2,i4)        (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                             __esvi4=i4,cgetmo(x1,&__esvi1, \
                                            &__esvi2,&__esvi3,x2,&__esvi4))
 
#define zgetmo(x1,i1,i2,i3,x2,i4)        (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                             __esvi4=i4,zgetmo(x1,&__esvi1, \
                                            &__esvi2,&__esvi3,x2,&__esvi4))
 
#define strmm(ch1,ch2,ch3,ch4,i1,i2,f1,x1,i3,x2,i4)   (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvf1=f1,\
                                          strmm(ch1,ch2,ch3,ch4,&__esvi1,\
                                          &__esvi2,&__esvf1,x1,&__esvi3,\
                                          x2,&__esvi4))
#define dtrmm(ch1,ch2,ch3,ch4,i1,i2,d1,x1,i3,x2,i4)   (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvd1=d1,\
                                          dtrmm(ch1,ch2,ch3,ch4,&__esvi1,\
                                          &__esvi2,&__esvd1,x1,&__esvi3,\
                                          x2,&__esvi4))
#define ctrmm(ch1,ch2,ch3,ch4,i1,i2,c1,x1,i3,x2,i4)   (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvc1=c1,\
                                          ctrmm(ch1,ch2,ch3,ch4,&__esvi1,\
                                          &__esvi2,&__esvc1,x1,&__esvi3,\
                                          x2,&__esvi4))
#define ztrmm(ch1,ch2,ch3,ch4,i1,i2,dc1,x1,i3,x2,i4)   (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvdc1=dc1,\
                                          ztrmm(ch1,ch2,ch3,ch4,&__esvi1,\
                                          &__esvi2,&__esvdc1,x1,&__esvi3,\
                                          x2,&__esvi4))
#define ssymm(ch1,ch2,i1,i2,f1,x1,i3,x2,i4,f2,x3,i5)  (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5,\
                                          __esvf1=f1,__esvf2=f2, ssymm(ch1,\
                                          ch2,&__esvi1,&__esvi2,&__esvf1,x1,\
                                          &__esvi3,x2,&__esvi4,&__esvf2,x3,\
                                          &__esvi5))
#define dsymm(ch1,ch2,i1,i2,d1,x1,i3,x2,i4,d2,x3,i5)  (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5,\
                                          __esvd1=d1,__esvd2=d2, dsymm(ch1,\
                                          ch2,&__esvi1,&__esvi2,&__esvd1,x1,\
                                          &__esvi3,x2,&__esvi4,&__esvd2,x3,\
                                          &__esvi5))
#define dsyr2k(ch1,ch2,i1,i2,d1,x1,i3,x2,i4,d2,x3,i5)  (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5,\
                                          __esvd2=d2, __esvd1=d1, dsyr2k(ch1,\
                                          ch2,&__esvi1,&__esvi2,&__esvd1,x1,\
                                          &__esvi3,x2,&__esvi4,&__esvd2,x3,\
                                          &__esvi5))
#define ssyr2k(ch1,ch2,i1,i2,f1,x1,i3,x2,i4,f2,x3,i5)  (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5,\
                                          __esvf1=f1, __esvf2=f2, ssyr2k(ch1,\
                                          ch2,&__esvi1,&__esvi2,&__esvf1,x1,\
                                          &__esvi3,x2,&__esvi4,&__esvf2,x3,\
                                          &__esvi5))
 
#define csyrk(ch1,ch2,i1,i2,c1,x1,i3,c2,x2,i4) (__esvi1=i1,__esvi2=i2,\
                                               __esvi3=i3,__esvi4=i4, \
                                               __esvc1=c1,__esvc2=c2, \
                                     csyrk(ch1,ch2,&__esvi1,&__esvi2,\
                         &__esvc1,x1,&__esvi3,&__esvc2,x2,&__esvi4))
 
#define zsyrk(ch1,ch2,i1,i2,dc1,x1,i3,dc2,x2,i4)        (__esvi1=i1,  \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                          __esvdc1=dc1,__esvdc2=dc2,  \
                                     zsyrk(ch1,ch2,&__esvi1,&__esvi2,\
                       &__esvdc1,x1,&__esvi3,&__esvdc2,x2,&__esvi4))
 
#define cherk(ch1,ch2,i1,i2,f1,x1,i3,f2,x2,i4)          (__esvi1=i1,  \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                              __esvf1=f1,__esvf2=f2,  \
                                     cherk(ch1,ch2,&__esvi1,&__esvi2,\
                         &__esvf1,x1,&__esvi3,&__esvf2,x2,&__esvi4))
 
#define zherk(ch1,ch2,i1,i2,d1,x1,i3,d2,x2,i4)          (__esvi1=i1,  \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                              __esvd1=d1,__esvd2=d2,  \
                                     zherk(ch1,ch2,&__esvi1,&__esvi2,\
                         &__esvd1,x1,&__esvi3,&__esvd2,x2,&__esvi4))
 
 
#define csyr2k(ch1,ch2,i1,i2,c1,x1,i3,x2,i4,c2,x3,i5)    (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvc1=c1,__esvc2=c2, \
                                    csyr2k(ch1,ch2,&__esvi1,&__esvi2,\
                         &__esvc1,x1,&__esvi3,x2,&__esvi4,&__esvc2,   \
                                                       x3,&__esvi5))
 
#define zsyr2k(ch1,ch2,i1,i2,dc1,x1,i3,x2,i4,dc2,x3,i5)  (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                __esvi5=i5,__esvdc1=dc1,__esvdc2=dc2, \
                                    zsyr2k(ch1,ch2,&__esvi1,&__esvi2,\
                         &__esvdc1,x1,&__esvi3,x2,&__esvi4,&__esvdc2, \
                                                       x3,&__esvi5))
 
#define cher2k(ch1,ch2,i1,i2,c1,x1,i3,x2,i4,f1,x3,i5)    (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvf1=f1,__esvc1=c1, \
                                    cher2k(ch1,ch2,&__esvi1,&__esvi2,\
                         &__esvc1,x1,&__esvi3,x2,&__esvi4,&__esvf1,   \
                                                       x3,&__esvi5))
 
#define zher2k(ch1,ch2,i1,i2,dc1,x1,i3,x2,i4,d1,x3,i5)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                  __esvi5=i5,__esvd1=d1,__esvdc1=dc1, \
                                    zher2k(ch1,ch2,&__esvi1,&__esvi2,\
                         &__esvdc1,x1,&__esvi3,x2,&__esvi4,&__esvd1,  \
                                                       x3,&__esvi5))
 
 
#define csymm(ch1,ch2,i1,i2,c1,x1,i3,x2,i4,c2,x3,i5)     (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvc1=c1,__esvc2=c2, \
                                    csymm(ch1,ch2,&__esvi1,&__esvi2, \
                         &__esvc1,x1,&__esvi3,x2,&__esvi4,&__esvc2,   \
                                                       x3,&__esvi5))
 
#define zsymm(ch1,ch2,i1,i2,dc1,x1,i3,x2,i4,dc2,x3,i5)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                __esvi5=i5,__esvdc1=dc1,__esvdc2=dc2, \
                                    zsymm(ch1,ch2,&__esvi1,&__esvi2, \
                         &__esvdc1,x1,&__esvi3,x2,&__esvi4,&__esvdc2, \
                                                       x3,&__esvi5))
 
#define chemm(ch1,ch2,i1,i2,c1,x1,i3,x2,i4,c2,x3,i5)     (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvc1=c1,__esvc2=c2, \
                                    chemm(ch1,ch2,&__esvi1,&__esvi2, \
                         &__esvc1,x1,&__esvi3,x2,&__esvi4,&__esvc2,   \
                                                       x3,&__esvi5))
 
#define zhemm(ch1,ch2,i1,i2,dc1,x1,i3,x2,i4,dc2,x3,i5)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                __esvi5=i5,__esvdc1=dc1,__esvdc2=dc2, \
                                    zhemm(ch1,ch2,&__esvi1,&__esvi2, \
                         &__esvdc1,x1,&__esvi3,x2,&__esvi4,&__esvdc2, \
                                                       x3,&__esvi5))
 
/*  Dense Linear Algebraic Equation Subroutines  */
 
 
#define sgef(x1,i1,i2,x2)                           (__esvi1=i1,__esvi2=i2, \
                                            sgef(x1,&__esvi1,&__esvi2,x2))
 
#define dgef(x1,i1,i2,x2)                           (__esvi1=i1,__esvi2=i2, \
                                            dgef(x1,&__esvi1,&__esvi2,x2))
 
#define cgef(x1,i1,i2,x2)                           (__esvi1=i1,__esvi2=i2, \
                                            cgef(x1,&__esvi1,&__esvi2,x2))
 
#define zgef(x1,i1,i2,x2)                           (__esvi1=i1,__esvi2=i2, \
                                            zgef(x1,&__esvi1,&__esvi2,x2))
 
 
#define sges(x1,i1,i2,x2,x3,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                sges(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
 
#define dges(x1,i1,i2,x2,x3,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                dges(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
 
#define cges(x1,i1,i2,x2,x3,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                cges(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
 
#define zges(x1,i1,i2,x2,x3,i3)          (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                zges(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
 
#define sgesm(ch1,x1,i1,i2,x2,x3,i3,i4)                  (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                               sgesm(ch1,x1,&__esvi1,&__esvi2,x2,x3, \
                                                 &__esvi3,&__esvi4))
 
#define dgesm(ch1,x1,i1,i2,x2,x3,i3,i4)                  (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                               dgesm(ch1,x1,&__esvi1,&__esvi2,x2,x3, \
                                                 &__esvi3,&__esvi4))
 
#define cgesm(ch1,x1,i1,i2,x2,x3,i3,i4)                  (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                               cgesm(ch1,x1,&__esvi1,&__esvi2,x2,x3, \
                                                 &__esvi3,&__esvi4))
 
#define zgesm(ch1,x1,i1,i2,x2,x3,i3,i4)                  (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                               zgesm(ch1,x1,&__esvi1,&__esvi2,x2,x3, \
                                                 &__esvi3,&__esvi4))
 
#ifdef __ESVERR
#define sgefcd(x1,i1,i2,x2,i3,x3,x4,x5,x6)          (__esvi1=i1,__esvi2=i2, \
                                             __esvi3=i3,sgefcd(x1,&__esvi1, \
                                         &__esvi2,x2,&__esvi3,x3,x4,x5,x6))
#else
#define sgefcd(x1,i1,i2,x2,i3,x3,x4,x5,i4)          (__esvi1=i1,__esvi2=i2, \
                                  __esvi3=i3,__esvi4=i4,sgefcd(x1,&__esvi1, \
                                   &__esvi2,x2,&__esvi3,x3,x4,x5,&__esvi4))
#endif
 
#ifdef __ESVERR
#define dgefcd(x1,i1,i2,x2,i3,x3,x4,x5,x6)          (__esvi1=i1,__esvi2=i2, \
                                             __esvi3=i3,dgefcd(x1,&__esvi1, \
                                         &__esvi2,x2,&__esvi3,x3,x4,x5,x6))
#else
#define dgefcd(x1,i1,i2,x2,i3,x3,x4,x5,i4)          (__esvi1=i1,__esvi2=i2, \
                                  __esvi3=i3,__esvi4=i4,dgefcd(x1,&__esvi1, \
                                   &__esvi2,x2,&__esvi3,x3,x4,x5,&__esvi4))
#endif
 
#define sppf(x1,i1,i2)  (__esvi1=i1,__esvi2=i2,sppf(x1,&__esvi1,&__esvi2))
 
#define dppf(x1,i1,i2)  (__esvi1=i1,__esvi2=i2,dppf(x1,&__esvi1,&__esvi2))
 
#define spps(x1,i1,x2,i2)                   (__esvi1=i1,__esvi2=i2,spps(x1, \
                                                     &__esvi1,x2,&__esvi2))
 
#define dpps(x1,i1,x2,i2)                   (__esvi1=i1,__esvi2=i2,dpps(x1, \
                                                     &__esvi1,x2,&__esvi2))
 
 
#ifdef __ESVERR
#define sppfcd(x1,i1,i2,x2,x3,x4,x5)      (__esvi1=i1,__esvi2=i2,sppfcd(x1, \
                                            &__esvi1,&__esvi2,x2,x3,x4,x5))
#else
#define sppfcd(x1,i1,i2,x2,x3,x4,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                           sppfcd(x1,&__esvi1,&__esvi2,x2,  \
                                                           x3,x4,&__esvi3))
#endif
 
#ifdef __ESVERR
#define dppfcd(x1,i1,i2,x2,x3,x4,x5)      (__esvi1=i1,__esvi2=i2,dppfcd(x1, \
                                            &__esvi1,&__esvi2,x2,x3,x4,x5))
#else
#define dppfcd(x1,i1,i2,x2,x3,x4,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                            dppfcd(x1,&__esvi1,&__esvi2,x2, \
                                                           x3,x4,&__esvi3))
#endif
 
#ifdef __ESVERR
#define sgeicd(x1,i1,i2,i3,x2,x3,x4,x5)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                      sgeicd(x1,&__esvi1,&__esvi2,&__esvi3, \
                                                              x2,x3,x4,x5))
#else
#define sgeicd(x1,i1,i2,i3,x2,x3,x4,i4)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvi4=i4,sgeicd(x1,&__esvi1,&__esvi2, \
                                               &__esvi3,x2,x3,x4,&__esvi4))
#endif
 
#ifdef __ESVERR
#define dgeicd(x1,i1,i2,i3,x2,x3,x4,x5)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                      dgeicd(x1,&__esvi1,&__esvi2,&__esvi3, \
                                                              x2,x3,x4,x5))
#else
#define dgeicd(x1,i1,i2,i3,x2,x3,x4,i4)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvi4=i4,dgeicd(x1,&__esvi1,&__esvi2, \
                                               &__esvi3,x2,x3,x4,&__esvi4))
#endif
 
#ifdef __ESVERR
#define sppicd(x1,i1,i2,x2,x3,x4,x5)     (__esvi1=i1,__esvi2=i2,sppicd(x1, \
                                           &__esvi1,&__esvi2,x2,x3,x4,x5))
#else
#define sppicd(x1,i1,i2,x2,x3,x4,i3)    (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                           sppicd(x1,&__esvi1,&__esvi2,x2, \
                                                          x3,x4,&__esvi3))
#endif
 
#ifdef __ESVERR
#define dppicd(x1,i1,i2,x2,x3,x4,x5)     (__esvi1=i1,__esvi2=i2,dppicd(x1, \
                                           &__esvi1,&__esvi2,x2,x3,x4,x5))
#else
#define dppicd(x1,i1,i2,x2,x3,x4,i3)    (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                          dppicd(x1,&__esvi1,&__esvi2,x2, \
                                                          x3,x4,&__esvi3))
#endif
 
#define strsm(ch1,ch2,ch3,ch4,i1,i2,f1,x1,i3,x2,i4) (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvf1=f1, \
                                            strsm(ch1,ch2,ch3,ch4,&__esvi1, \
                                              &__esvi2,&__esvf1,x1,&__esvi3, \
                                                              x2,&__esvi4))
 
#define dtrsm(ch1,ch2,ch3,ch4,i1,i2,d1,x1,i3,x2,i4) (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvd1=d1, \
                                            dtrsm(ch1,ch2,ch3,ch4,&__esvi1, \
                                              &__esvi2,&__esvd1,x1,&__esvi3, \
                                                              x2,&__esvi4))
#define ctrsm(ch1,ch2,ch3,ch4,i1,i2,c1,x1,i3,x2,i4)    (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvc1=c1,\
                                          ctrsm(ch1,ch2,ch3,ch4,&__esvi1,\
                                          &__esvi2,&__esvc1,x1,&__esvi3,x2,\
                                          &__esvi4))
#define ztrsm(ch1,ch2,ch3,ch4,i1,i2,dc1,x1,i3,x2,i4)    (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3,__esvi4=i4,__esvdc1=dc1,\
                                          ztrsm(ch1,ch2,ch3,ch4,&__esvi1,\
                                          &__esvi2,&__esvdc1,x1,&__esvi3,x2,\
                                          &__esvi4))
 
#define strsv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3, strsv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
#define dtrsv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3, dtrsv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
#define ctrsv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3, ctrsv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
#define ztrsv(ch1,ch2,ch3,i1,x1,i2,x2,i3)    (__esvi1=i1,__esvi2=i2,\
                                          __esvi3=i3, ztrsv(ch1,ch2,ch3,\
                                          &__esvi1,x1,&__esvi2,x2,&__esvi3))
 
#define spof(ch1,x1,i1,i2)       (__esvi1=i1,__esvi2=i2,spof(ch1,x1, \
                                           &__esvi1,&__esvi2))
 
#define dpof(ch1,x1,i1,i2)       (__esvi1=i1,__esvi2=i2,dpof(ch1,x1, \
                                           &__esvi1,&__esvi2))
 
#define cpof(ch1,x1,i1,i2)       (__esvi1=i1,__esvi2=i2,cpof(ch1,x1, \
                                           &__esvi1,&__esvi2))
 
#define zpof(ch1,x1,i1,i2)       (__esvi1=i1,__esvi2=i2,zpof(ch1,x1, \
                                           &__esvi1,&__esvi2))
 
#define sposm(ch1,x1,i1,i2,x2,i3,i4)                     (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                  sposm(ch1,x1,&__esvi1,&__esvi2,x2, \
                                                 &__esvi3,&__esvi4))
 
#define dposm(ch1,x1,i1,i2,x2,i3,i4)                     (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                  dposm(ch1,x1,&__esvi1,&__esvi2,x2, \
                                                 &__esvi3,&__esvi4))
 
#define cposm(ch1,x1,i1,i2,x2,i3,i4)                     (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                  cposm(ch1,x1,&__esvi1,&__esvi2,x2, \
                                                 &__esvi3,&__esvi4))
 
#define zposm(ch1,x1,i1,i2,x2,i3,i4)                     (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                  zposm(ch1,x1,&__esvi1,&__esvi2,x2, \
                                                 &__esvi3,&__esvi4))
 
#ifdef __ESVERR
#define spofcd(ch1,x1,i1,i2,i3,x2,x3,x4,x5)  (__esvi1=i1,__esvi2=i2,  \
                                  __esvi3=i3,spofcd(ch1,x1,&__esvi1, \
                                  &__esvi2,&__esvi3,x2,x3,x4,x5))
#else
#define spofcd(ch1,x1,i1,i2,i3,x2,x3,x4,i4)  (__esvi1=i1,__esvi2=i2,  \
                                __esvi3=i3,__esvi4=i4,spofcd(ch1,x1, \
                                &__esvi1,&__esvi2,&__esvi3,x2,x3,x4,  \
                                                          &__esvi4))
#endif
 
#ifdef __ESVERR
#define dpofcd(ch1,x1,i1,i2,i3,x2,x3,x4,x5)  (__esvi1=i1,__esvi2=i2,  \
                                  __esvi3=i3,dpofcd(ch1,x1,&__esvi1, \
                                  &__esvi2,&__esvi3,x2,x3,x4,x5))
#else
#define dpofcd(ch1,x1,i1,i2,i3,x2,x3,x4,i4)  (__esvi1=i1,__esvi2=i2,  \
                                __esvi3=i3,__esvi4=i4,dpofcd(ch1,x1, \
                                &__esvi1,&__esvi2,&__esvi3,x2,x3,x4,  \
                                                          &__esvi4))
#endif
 
#ifdef __ESVERR
#define spoicd(ch1,x1,i1,i2,i3,x2,x3,x4,x5)  (__esvi1=i1,__esvi2=i2,  \
                                  __esvi3=i3,spoicd(ch1,x1,&__esvi1, \
                                  &__esvi2,&__esvi3,x2,x3,x4,x5))
#else
#define spoicd(ch1,x1,i1,i2,i3,x2,x3,x4,i4)  (__esvi1=i1,__esvi2=i2,  \
                                __esvi3=i3,__esvi4=i4,spoicd(ch1,x1, \
                                &__esvi1,&__esvi2,&__esvi3,x2,x3,x4,  \
                                                          &__esvi4))
#endif
 
#ifdef __ESVERR
#define dpoicd(ch1,x1,i1,i2,i3,x2,x3,x4,x5)  (__esvi1=i1,__esvi2=i2,  \
                                  __esvi3=i3,dpoicd(ch1,x1,&__esvi1, \
                                  &__esvi2,&__esvi3,x2,x3,x4,x5))
#else
#define dpoicd(ch1,x1,i1,i2,i3,x2,x3,x4,i4)  (__esvi1=i1,__esvi2=i2,  \
                                __esvi3=i3,__esvi4=i4,dpoicd(ch1,x1, \
                                &__esvi1,&__esvi2,&__esvi3,x2,x3,x4,  \
                                                          &__esvi4))
#endif
 
#define stri(ch1,ch2,x1,i1,i2)              (__esvi1=i1,__esvi2=i2,  \
                               stri(ch1,ch2,x1,&__esvi1,&__esvi2))
 
#define dtri(ch1,ch2,x1,i1,i2)              (__esvi1=i1,__esvi2=i2,  \
                               dtri(ch1,ch2,x1,&__esvi1,&__esvi2))
 
#define stpi(ch1,ch2,x1,i1)     (__esvi1=i1,stpi(ch1,ch2,x1,&__esvi1))
#define dtpi(ch1,ch2,x1,i1)     (__esvi1=i1,dtpi(ch1,ch2,x1,&__esvi1))
 
 
#define stpsv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               stpsv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
#define dtpsv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               dtpsv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
#define ctpsv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               ctpsv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
#define ztpsv(ch1,ch2,ch3,i1,x1,x2,i2)    (__esvi1=i1,__esvi2=i2,     \
                                               ztpsv(ch1,ch2,ch3,    \
                                           &__esvi1,x1,x2,&__esvi2))
 
 
#define stbsv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        stbsv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
#define dtbsv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        dtbsv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
#define ctbsv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        ctbsv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
#define ztbsv(ch1,ch2,ch3,i1,i2,x1,i3,x2,i4)             (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                        ztbsv(ch1,ch2,ch3,&__esvi1,  \
                                  &__esvi2,x1,&__esvi3,x2,&__esvi4))
 
#define cgtnpf(i1,x1,x2,x3,i2) (__esvi1=i1,__esvi2=i2,cgtnpf(&__esvi1,\
                                                  x1,x2,x3,&__esvi2))
 
#define zgtnpf(i1,x1,x2,x3,i2) (__esvi1=i1,__esvi2=i2,zgtnpf(&__esvi1,\
                                                  x1,x2,x3,&__esvi2))
 
#define cgtnps(i1,x1,x2,x3,x4)          (__esvi1=i1,cgtnps(&__esvi1, \
                                                       x1,x2,x3,x4))
 
#define zgtnps(i1,x1,x2,x3,x4)          (__esvi1=i1,zgtnps(&__esvi1, \
                                                       x1,x2,x3,x4))
 
#ifdef __ESVERR
#define dsris(ch1,ch2,i1,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)          \
        (__esvi1=i1,dsris(ch1,ch2,&__esvi1,x1,x2,x3,x4,x5,x6,x7,x8,  \
         x9,x10,x11))
#else
#define dsris(ch1,ch2,i1,x1,x2,x3,x4,x5,x6,x7,x8,i2,x9,i3)            \
        (__esvi1=i1,__esvi2=i2,__esvi3=i3,dsris(ch1,ch2,&__esvi1,x1, \
                         x2,x3,x4,x5,x6,x7,x8,&__esvi2,x9,&__esvi3))
#endif
 
 
/*  Banded Linear Algebraic Equation Subroutines  */
 
 
#define sgbf(x1,i1,i2,i3,i4,x2) (__esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                           sgbf(x1,&__esvi1,&__esvi2,&__esvi3,&__esvi4,x2))
 
#define dgbf(x1,i1,i2,i3,i4,x2) (__esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                           dgbf(x1,&__esvi1,&__esvi2,&__esvi3,&__esvi4,x2))
 
#define sgbs(x1,i1,i2,i3,i4,x2,x3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                       __esvi4=i4,sgbs(x1,&__esvi1,&__esvi2, \
                                                    &__esvi3,&__esvi4,x2,x3))
 
#define dgbs(x1,i1,i2,i3,i4,x2,x3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                       __esvi4=i4,dgbs(x1,&__esvi1,&__esvi2, \
                                                    &__esvi3,&__esvi4,x2,x3))
 
 
#define spbf(x1,i1,i2,i3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   spbf(x1,&__esvi1,&__esvi2,&__esvi3))
 
#define dpbf(x1,i1,i2,i3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   dpbf(x1,&__esvi1,&__esvi2,&__esvi3))
 
#define spbs(x1,i1,i2,i3,x2)      (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   spbs(x1,&__esvi1,&__esvi2,&__esvi3,x2))
 
#define dpbs(x1,i1,i2,i3,x2)      (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   dpbs(x1,&__esvi1,&__esvi2,&__esvi3,x2))
 
#define spbchf(x1,i1,i2,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                 spbchf(x1,&__esvi1,&__esvi2,&__esvi3))
 
#define dpbchf(x1,i1,i2,i3)     (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                               dpbchf(x1,&__esvi1,&__esvi2,&__esvi3))
 
#define spbchs(x1,i1,i2,i3,x2)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                               spbchs(x1,&__esvi1,&__esvi2,&__esvi3,x2))
 
#define dpbchs(x1,i1,i2,i3,x2)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                               dpbchs(x1,&__esvi1,&__esvi2,&__esvi3,x2))
 
 
#define sgtf(i1,x1,x2,x3,x4,x5)  (__esvi1=i1,sgtf(&__esvi1,x1,x2,x3,x4,x5))
#define dgtf(i1,x1,x2,x3,x4,x5)  (__esvi1=i1,dgtf(&__esvi1,x1,x2,x3,x4,x5))
 
#define sgts(i1,x1,x2,x3,x4,x5,x6)             (__esvi1=i1,sgts(&__esvi1,x1, \
                                                            x2,x3,x4,x5,x6))
 
#define dgts(i1,x1,x2,x3,x4,x5,x6)             (__esvi1=i1,dgts(&__esvi1,x1, \
                                                            x2,x3,x4,x5,x6))
 
#define sgtnpf(i1,x1,x2,x3,i2)       (__esvi1=i1,__esvi2=i2,sgtnpf(&__esvi1, \
                                                         x1,x2,x3,&__esvi2))
 
#define dgtnpf(i1,x1,x2,x3,i2)       (__esvi1=i1,__esvi2=i2,dgtnpf(&__esvi1, \
                                                         x1,x2,x3,&__esvi2))
 
#define sgtnps(i1,x1,x2,x3,x4)    (__esvi1=i1,sgtnps(&__esvi1,x1,x2,x3,x4))
 
#define dgtnps(i1,x1,x2,x3,x4)    (__esvi1=i1,dgtnps(&__esvi1,x1,x2,x3,x4))
 
#define sgtnp(i1,x1,x2,x3,x4)     (__esvi1=i1,sgtnp(&__esvi1,x1,x2,x3,x4))
#define dgtnp(i1,x1,x2,x3,x4)     (__esvi1=i1,dgtnp(&__esvi1,x1,x2,x3,x4))
#define cgtnp(i1,x1,x2,x3,x4)     (__esvi1=i1,cgtnp(&__esvi1,x1,x2,x3,x4))
#define zgtnp(i1,x1,x2,x3,x4)     (__esvi1=i1,zgtnp(&__esvi1,x1,x2,x3,x4))
#define sptf(i1,x1,x2,i2)              (__esvi1=i1,__esvi2=i2,sptf(&__esvi1, \
                                                            x1,x2,&__esvi2))
 
#define dptf(i1,x1,x2,i2)              (__esvi1=i1,__esvi2=i2,dptf(&__esvi1, \
                                                            x1,x2,&__esvi2))
 
#define spts(i1,x1,x2,x3)              (__esvi1=i1,spts(&__esvi1,x1,x2,x3))
 
#define dpts(i1,x1,x2,x3)              (__esvi1=i1,dpts(&__esvi1,x1,x2,x3))
 
 
 
/* Sparse Linear Algebraic Equations Subroutines   */
 
 
#ifdef  __ESVERR
#define dgsf(i1,i2,i3,x1,x2,x3,i4,x4,x5,x6,x7,x8)   (__esvi1=i1,__esvi2=i2, \
                                                      __esvi3=i3,__esvi4=i4, \
                                           dgsf(&__esvi1,&__esvi2,&__esvi3, \
                                         x1,x2,x3,&__esvi4,x4,x5,x6,x7,x8))
#else
#define dgsf(i1,i2,i3,x1,x2,x3,i4,x4,x5,x6,x7,i5)   (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                           dgsf(&__esvi1,&__esvi2,&__esvi3, \
                                              x1,x2,x3,&__esvi4,x4,x5,x6,x7, \
                                                                 &__esvi5))
#endif
 
#ifdef  __ESVERR
#define dgss(i1,i2,x1,x2,x3,i3,x4,x5,x6)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                                     dgss(&__esvi1,&__esvi2, \
                                                x1,x2,x3,&__esvi3,x4,x5,x6))
#else
#define dgss(i1,i2,x1,x2,x3,i3,x4,x5,i4)  (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                          __esvi4=i4,dgss(&__esvi1,&__esvi2, \
                                          x1,x2,x3,&__esvi3,x4,x5,&__esvi4))
#endif
 
#ifdef __ESVERR
#define dgkfs(i1,x1,i2,x2,x3,i3,x4,x5,x6,x7,x9,x8,i5,i6)       (__esvi1=i1, \
                                           __esvi5=i5,__esvi2=i2,__esvi6=i6, \
                                     __esvi3=i3,dgkfs(&__esvi1,x1,&__esvi2, \
                       x2,x3,&__esvi3,x4,x5,x6,x7,x9,x8,&__esvi5,&__esvi6))
#else
#define dgkfs(i1,x1,i2,x2,x3,i3,x4,x5,x6,x7,i4,x8,i5,i6)       (__esvi1=i1, \
                                           __esvi5=i5,__esvi2=i2,__esvi6=i6, \
                          __esvi3=i3,__esvi4=i4,dgkfs(&__esvi1,x1,&__esvi2, \
                 x2,x3,&__esvi3,x4,x5,x6,x7,&__esvi4,x8,&__esvi5,&__esvi6))
#endif
 
#ifdef  __ESVERR
#define dskfs(i1,x1,i2,x2,x3,x4,x5,x7,x6,i4,i5)     (__esvi1=i1,__esvi2=i2, \
                          __esvi4=i4,__esvi5=i5,dskfs(&__esvi1,x1,&__esvi2, \
                                      x2,x3,x4,x5,x7,x6,&__esvi4,&__esvi5))
#else
#define dskfs(i1,x1,i2,x2,x3,x4,x5,i3,x6,i4,i5)     (__esvi1=i1,__esvi2=i2, \
               __esvi3=i3,__esvi4=i4,__esvi5=i5,dskfs(&__esvi1,x1,&__esvi2, \
                                x2,x3,x4,x5,&__esvi3,x6,&__esvi4,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define dsmcg(i1,i2,x1,x2,i3,x3,x4,x5,x6,x7,x9,x8,x10)        (__esvi1=i1,  \
                                                      __esvi2=i2,__esvi3=i3, \
                                 dsmcg(&__esvi1,&__esvi2,x1,x2,&__esvi3,x3, \
                                                    x4,x5,x6,x7,x9,x8,x10))
#else
#define dsmcg(i1,i2,x1,x2,i3,x3,x4,x5,x6,x7,i4,x8,i5)         (__esvi1=i1,  \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                 dsmcg(&__esvi1,&__esvi2,x1,x2,&__esvi3,x3, \
                                         x4,x5,x6,x7,&__esvi4,x8,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define dsdcg(i1,i2,i3,x1,i4,x2,x3,x4,x5,x6,x7,x9,x8,x10)     (__esvi1=i1,  \
                                           __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                          dsdcg(&__esvi1,&__esvi2,&__esvi3, \
                                              x1,&__esvi4,x2,x3,x4,x5,x6,x7, \
                                                                x9,x8,x10))
#else
#define dsdcg(i1,i2,i3,x1,i4,x2,x3,x4,x5,x6,x7,i5,x8,i6)      (__esvi1=i1,  \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                               __esvi6=i6,dsdcg(&__esvi1,&__esvi2,&__esvi3, \
                                              x1,&__esvi4,x2,x3,x4,x5,x6,x7, \
                                                     &__esvi5,x8,&__esvi6))
#endif
#ifdef  __ESVERR
#define dsmgcg(i1,i2,x1,x2,i3,x3,x4,x5,x6,x7,x9,x8,x10)        (__esvi1=i1, \
                                                      __esvi2=i2,__esvi3=i3, \
                                dsmgcg(&__esvi1,&__esvi2,x1,x2,&__esvi3,x3, \
                                                    x4,x5,x6,x7,x9,x8,x10))
#else
#define dsmgcg(i1,i2,x1,x2,i3,x3,x4,x5,x6,x7,i4,x8,i5)         (__esvi1=i1, \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                dsmgcg(&__esvi1,&__esvi2,x1,x2,&__esvi3,x3, \
                                         x4,x5,x6,x7,&__esvi4,x8,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define dsdgcg(i1,i2,x1,i3,x2,x3,x4,x5,x6,x7,x9,x8,x10)        (__esvi1=i1, \
                                                      __esvi2=i2,__esvi3=i3, \
                                dsdgcg(&__esvi1,&__esvi2,x1,&__esvi3,x2,x3, \
                                                    x4,x5,x6,x7,x9,x8,x10))
#else
#define dsdgcg(i1,i2,x1,i3,x2,x3,x4,x5,x6,x7,i4,x8,i5)         (__esvi1=i1, \
                                __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                dsdgcg(&__esvi1,&__esvi2,x1,&__esvi3,x2,x3, \
                                         x4,x5,x6,x7,&__esvi4,x8,&__esvi5))
#endif
 
 
/* Linear Least Square Subroutines  */
 
 
#ifdef  __ESVERR
#define sgesvf(i1,x1,i2,x2,i3,i4,x3,i5,i6,x4,x5)    (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                              __esvi6=i6,sgesvf(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvi4,x3,&__esvi5, \
                                                           &__esvi6,x4,x5))
#else
#define sgesvf(i1,x1,i2,x2,i3,i4,x3,i5,i6,x4,i7)    (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                  __esvi6=i6,__esvi7=i7,sgesvf(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvi4,x3,&__esvi5, \
                                                     &__esvi6,x4,&__esvi7))
#endif
 
#ifdef  __ESVERR
#define dgesvf(i1,x1,i2,x2,i3,i4,x3,i5,i6,x4,x5)    (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                             __esvi6=i6,dgesvf(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvi4,x3,&__esvi5, \
                                                           &__esvi6,x4,x5))
#else
#define dgesvf(i1,x1,i2,x2,i3,i4,x3,i5,i6,x4,i7)    (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                   __esvi6=i6,__esvi7=i7,dgesvf(&__esvi1,x1, \
                                  &__esvi2,x2,&__esvi3,&__esvi4,x3,&__esvi5, \
                                                     &__esvi6,x4,&__esvi7))
#endif
 
#define sgesvs(x1,i1,x2,i2,i3,x3,x4,i4,i5,i6,f1)    (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                  __esvi6=i6,__esvf1=f1,sgesvs(x1,&__esvi1, \
                                        x2,&__esvi2,&__esvi3,x3,x4,&__esvi4, \
                                               &__esvi5,&__esvi6,&__esvf1))
 
 
#define dgesvs(x1,i1,x2,i2,i3,x3,x4,i4,i5,i6,d1)    (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                  __esvi6=i6,__esvd1=d1,dgesvs(x1,&__esvi1, \
                                        x2,&__esvi2,&__esvi3,x3,x4,&__esvi4, \
                                               &__esvi5,&__esvi6,&__esvd1))
 
 
#ifdef  __ESVERR
#define sgells(i1,x1,i2,x2,i3,x3,i4,x4,f1,i5,i6,i7,x6,x5,x7)  (__esvf1=f1,  \
                                __esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                           __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                sgells(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                     &__esvi4,x4,&__esvf1,&__esvi5,&__esvi6, \
                                                        &__esvi7,x6,x5,x7))
#else
#define sgells(i1,x1,i2,x2,i3,x3,i4,x4,f1,i5,i6,i7,x6,x5,i9)  (__esvf1=f1,  \
                                __esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                __esvi5=i5,__esvi6=i6,__esvi7=i7,__esvi9=i9, \
                                sgells(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                     &__esvi4,x4,&__esvf1,&__esvi5,&__esvi6, \
                                                  &__esvi7,x6,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define dgells(i1,x1,i2,x2,i3,x3,i4,x4,d1,i5,i6,i7,x6,x5,x7)   (__esvd1=d1, \
                                __esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                           __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                dgells(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                     &__esvi4,x4,&__esvd1,&__esvi5,&__esvi6, \
                                                        &__esvi7,x6,x5,x7))
#else
#define dgells(i1,x1,i2,x2,i3,x3,i4,x4,d1,i5,i6,i7,x6,x5,i9)   (__esvd1=d1, \
                                __esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                __esvi5=i5,__esvi6=i6,__esvi7=i7,__esvi9=i9, \
                                dgells(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                     &__esvi4,x4,&__esvd1,&__esvi5,&__esvi6, \
                                                  &__esvi7,x6,x5,&__esvi9))
#endif
 
 
/* Eigensystem Analysis   Subroutines   */
 
 
#ifdef  __ESVERR
#define sgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,x6)       (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                         sgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                              &__esvi3,x4,&__esvi4,x5,x6))
#else
#define sgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,i5)       (__esvi1=i1,__esvi2=i2, \
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         sgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                        &__esvi3,x4,&__esvi4,x5,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define dgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,x6)       (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                         dgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                              &__esvi3,x4,&__esvi4,x5,x6))
#else
#define dgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,i5)       (__esvi1=i1,__esvi2=i2, \
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         dgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                        &__esvi3,x4,&__esvi4,x5,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define cgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,x6)       (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                         cgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                              &__esvi3,x4,&__esvi4,x5,x6))
#else
#define cgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,i5)       (__esvi1=i1,__esvi2=i2, \
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         cgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                        &__esvi3,x4,&__esvi4,x5,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define zgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,x6)       (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                         zgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                              &__esvi3,x4,&__esvi4,x5,x6))
#else
#define zgeev(i1,x1,i2,x2,x3,i3,x4,i4,x5,i5)       (__esvi1=i1,__esvi2=i2, \
                                          __esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                         zgeev(&__esvi1,x1,&__esvi2,x2,x3, \
                                        &__esvi3,x4,&__esvi4,x5,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define sspev(i1,x1,x2,x3,i2,i3,x4,x5)             (__esvi1=i1,__esvi2=i2, \
                                                __esvi3=i3,sspev(&__esvi1, \
                                        x1,x2,x3,&__esvi2,&__esvi3,x4,x5))
#else
#define sspev(i1,x1,x2,x3,i2,i3,x4,i4)             (__esvi1=i1,__esvi2=i2, \
                                     __esvi3=i3,__esvi4=i4,sspev(&__esvi1, \
                                                x1,x2,x3,&__esvi2,&__esvi3, \
                                                             x4,&__esvi4))
#endif
 
#ifdef  __ESVERR
#define dspev(i1,x1,x2,x3,i2,i3,x4,x5)             (__esvi1=i1,__esvi2=i2, \
                                                __esvi3=i3,dspev(&__esvi1, \
                                       x1, x2,x3,&__esvi2,&__esvi3,x4,x5))
#else
#define dspev(i1,x1,x2,x3,i2,i3,x4,i4)             (__esvi1=i1,__esvi2=i2, \
                                     __esvi3=i3,__esvi4=i4,dspev(&__esvi1, \
                                               x1, x2,x3,&__esvi2,&__esvi3, \
                                                             x4,&__esvi4))
#endif
 
#ifdef  __ESVERR
#define chpev(i1,x1,x2,x3,i2,i3,x4,x5)             (__esvi1=i1,__esvi2=i2, \
                                                __esvi3=i3,chpev(&__esvi1, \
                                        x1,x2,x3,&__esvi2,&__esvi3,x4,x5))
#else
#define chpev(i1,x1,x2,x3,i2,i3,x4,i4)             (__esvi1=i1,__esvi2=i2, \
                                     __esvi3=i3,__esvi4=i4,chpev(&__esvi1, \
                                                x1,x2,x3,&__esvi2,&__esvi3, \
                                                             x4,&__esvi4))
#endif
 
#ifdef  __ESVERR
#define zhpev(i1,x1,x2,x3,i2,i3,x4,x5)             (__esvi1=i1,__esvi2=i2, \
                                                __esvi3=i3,zhpev(&__esvi1, \
                                        x1,x2,x3,&__esvi2,&__esvi3,x4,x5))
#else
#define zhpev(i1,x1,x2,x3,i2,i3,x4,i4)             (__esvi1=i1,__esvi2=i2, \
                                     __esvi3=i3,__esvi4=i4,zhpev(&__esvi1, \
                                                x1,x2,x3,&__esvi2,&__esvi3, \
                                                             x4,&__esvi4))
#endif
 
#ifdef  __ESVERR
#define sspsv(i1,x1,x2,x3,i2,i3,i4,x4,x5)          (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                sspsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                          &__esvi4,x4,x5))
#else
#define sspsv(i1,x1,x2,x3,i2,i3,i4,x4,i5)          (__esvi1=i1,__esvi2=i2, \
                                          __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                sspsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                    &__esvi4,x4,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define dspsv(i1,x1,x2,x3,i2,i3,i4,x4,x5)          (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                dspsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                          &__esvi4,x4,x5))
#else
#define dspsv(i1,x1,x2,x3,i2,i3,i4,x4,i5)          (__esvi1=i1,__esvi2=i2, \
                                          __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                dspsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                    &__esvi4,x4,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define chpsv(i1,x1,x2,x3,i2,i3,i4,x4,x5)          (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                chpsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                          &__esvi4,x4,x5))
#else
#define chpsv(i1,x1,x2,x3,i2,i3,i4,x4,i5)          (__esvi1=i1,__esvi2=i2, \
                                          __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                chpsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                    &__esvi4,x4,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define zhpsv(i1,x1,x2,x3,i2,i3,i4,x4,x5)          (__esvi1=i1,__esvi2=i2, \
                                                     __esvi3=i3,__esvi4=i4, \
                                zhpsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                          &__esvi4,x4,x5))
#else
#define zhpsv(i1,x1,x2,x3,i2,i3,i4,x4,i5)          (__esvi1=i1,__esvi2=i2, \
                                          __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                zhpsv(&__esvi1,x1,x2,x3,&__esvi2,&__esvi3, \
                                                    &__esvi4,x4,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define sgegv(i1,x1,i2,x2,i3,x3,x4,x5,i4,i5,x6,x7)  (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                 sgegv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                            x4,x5,&__esvi4,&__esvi5,x6,x7))
#else
#define sgegv(i1,x1,i2,x2,i3,x3,x4,x5,i4,i5,x6,i6)  (__esvi1=i1,__esvi2=i2, \
                                __esvi5=i5,__esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                 sgegv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                      x4,x5,&__esvi4,&__esvi5,x6,&__esvi6))
#endif
 
#ifdef  __ESVERR
#define dgegv(i1,x1,i2,x2,i3,x3,x4,x5,i4,i5,x6,x7)  (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                 dgegv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                            x4,x5,&__esvi4,&__esvi5,x6,x7))
#else
#define dgegv(i1,x1,i2,x2,i3,x3,x4,x5,i4,i5,x6,i6)  (__esvi1=i1,__esvi2=i2, \
                                __esvi5=i5,__esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                 dgegv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                      x4,x5,&__esvi4,&__esvi5,x6,&__esvi6))
#endif
 
#ifdef  __ESVERR
#define ssygv(i1,x1,i2,x2,i3,x3,x4,i4,i5,x5,x6)     (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                 ssygv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                               x4,&__esvi4,&__esvi5,x5,x6))
#else
#define ssygv(i1,x1,i2,x2,i3,x3,x4,i4,i5,x5,i6)     (__esvi1=i1,__esvi2=i2, \
                                __esvi5=i5,__esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                 ssygv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                         x4,&__esvi4,&__esvi5,x5,&__esvi6))
#endif
 
#ifdef  __ESVERR
#define dsygv(i1,x1,i2,x2,i3,x3,x4,i4,i5,x5,x6)     (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                 dsygv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                               x4,&__esvi4,&__esvi5,x5,x6))
#else
#define dsygv(i1,x1,i2,x2,i3,x3,x4,i4,i5,x5,i6)     (__esvi1=i1,__esvi2=i2, \
                                __esvi5=i5,__esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                 dsygv(&__esvi1,x1,&__esvi2,x2,&__esvi3,x3, \
                                         x4,&__esvi4,&__esvi5,x5,&__esvi6))
#endif
 
 
/*  Fourier Transforms  Subroutines   */
 
 
#ifdef  __ESVERR
#define scft(i1,x1,i2,i3,x2,i4,i5,x7,i7,i8,f1,x3,x5,x4,x6)     (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                __esvi3=i3,__esvi4=i4,__esvi8=i8,__esvf1=f1, \
                                     scft(&__esvi1,x1,&__esvi2,&__esvi3,x2, \
                                              &__esvi4,&__esvi5,x7,&__esvi7, \
                                            &__esvi8,&__esvf1,x3,x5,x4,x6))
#else
#define scft(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,f1,x3,i9,x4,i10)    (__esvi1=i1, \
                                                      __esvi2=i2,__esvi7=i7, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi8=i8,__esvi9=i9, \
                                                    __esvi10=i10,__esvf1=f1, \
                                                          scft(&__esvi1,x1, \
                                                       &__esvi2,&__esvi3,x2, \
                                                          &__esvi4,&__esvi5, \
                                                          &__esvi6,&__esvi7, \
                                                       &__esvi8,&__esvf1,x3, \
                                                    &__esvi9,x4,&__esvi10))
#endif
 
#ifdef  __ESVERR
#define dcft(i1,x1,i2,i3,x2,i4,i5,x7,i7,i8,d1,x3,x5,x4,x6)     (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                __esvi3=i3,__esvi4=i4,__esvi8=i8,__esvd1=d1, \
                                     dcft(&__esvi1,x1,&__esvi2,&__esvi3,x2, \
                                              &__esvi4,&__esvi5,x7,&__esvi7, \
                                            &__esvi8,&__esvd1,x3,x5,x4,x6))
#else
#define dcft(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,d1,x3,i9,x4,i10)    (__esvi1=i1, \
                                                      __esvi2=i2,__esvi7=i7, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi8=i8,__esvi9=i9, \
                                                    __esvi10=i10,__esvd1=d1, \
                                                          dcft(&__esvi1,x1, \
                                                       &__esvi2,&__esvi3,x2, \
                                                          &__esvi4,&__esvi5, \
                                                          &__esvi6,&__esvi7, \
                                                       &__esvi8,&__esvd1,x3, \
                                                    &__esvi9,x4,&__esvi10))
#endif
 
#ifdef  __ESVERR
#define srcft(i1,x1,i2,x2,i3,x9,i5,i6,f1,x3,x6,x4,x7,x5,x8)    (__esvi1=i1, \
                                           __esvi2=i2,__esvi5=i5,__esvi6=i6, \
                                                      __esvi3=i3,__esvf1=f1, \
                                    srcft(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                           x9,&__esvi5,&__esvi6,&__esvf1,x3, \
                                                           x6,x4,x7,x5,x8))
#else
#define srcft(i1,x1,i2,x2,i3,i4,i5,i6,f1,x3,i7,x4,i8,x5,i9)    (__esvi1=i1, \
                                                      __esvi2=i2,__esvi7=i7, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi8=i8,__esvi9=i9, \
                                                                 __esvf1=f1, \
                                                         srcft(&__esvi1,x1, \
                                                       &__esvi2,x2,&__esvi3, \
                                                          &__esvi4,&__esvi5, \
                                                       &__esvi6,&__esvf1,x3, \
                                                       &__esvi7,x4,&__esvi8, \
                                                              x5,&__esvi9))
#endif
 
#ifdef __ESVERR
#define drcft(i1,x1,i2,x2,i3,x5,i5,i6,d1,x3,x6,x4,x7)  (__esvi1=i1, \
                         __esvi2=i2,__esvi3=i3,__esvi5=i5,__esvi6=i6, \
                                       __esvd1=d1,drcft(&__esvi1,x1, \
                           &__esvi2,x2,&__esvi3,x5,&__esvi5,&__esvi6, \
                                              &__esvd1,x3,x6,x4,x7))
#else
#define drcft(i1,x1,i2,x2,i3,i4,i5,i6,d1,x3,i7,x4,i8)  (__esvi1=i1, \
                     __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5,__esvi6=i6, \
                     __esvi7=i7,__esvi8=i8,__esvd1=d1,drcft(&__esvi1,x1, \
                     &__esvi2,x2,&__esvi3,&__esvi4,&__esvi5,&__esvi6, \
                     &__esvd1,x3,&__esvi7,x4,&__esvi8))
#endif
 
#ifdef  __ESVERR
#define scrft(i1,x1,i2,x2,i3,x9,i5,i6,f1,x3,x6,x4,x7,x5,x8)    (__esvi1=i1, \
                                           __esvi2=i2,__esvi5=i5,__esvi6=i6, \
                                                      __esvi3=i3,__esvf1=f1, \
                                    scrft(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                           x9,&__esvi5,&__esvi6,&__esvf1,x3, \
                                                           x6,x4,x7,x5,x8))
#else
#define scrft(i1,x1,i2,x2,i3,i4,i5,i6,f1,x3,i7,x4,i8,x5,i9)    (__esvi1=i1, \
                                                      __esvi2=i2,__esvi7=i7, \
                                                      __esvi5=i5,__esvi6=i6, \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi8=i8,__esvi9=i9, \
                                                                 __esvf1=f1, \
                                                         scrft(&__esvi1,x1, \
                                                       &__esvi2,x2,&__esvi3, \
                                                          &__esvi4,&__esvi5, \
                                                       &__esvi6,&__esvf1,x3, \
                                                       &__esvi7,x4,&__esvi8, \
                                                              x5,&__esvi9))
#endif
 
#ifdef __ESVERR
#define dcrft(i1,x1,i2,x2,i3,x5,i5,i6,d1,x3,x6,x4,x7)      (__esvi1=i1, \
                            __esvi2=i2,__esvi3=i3,__esvi5=i5,__esvi6=i6, \
                                           __esvd1=d1,dcrft(&__esvi1,x1, \
                              &__esvi2,x2,&__esvi3,x5,&__esvi5,&__esvi6, \
                                                 &__esvd1,x3,x6,x4,x7))
#else
#define dcrft(i1,x1,i2,x2,i3,i4,i5,i6,d1,x3,i7,x4,i8)  (__esvi1=i1, \
                     __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5,__esvi6=i6, \
                     __esvi7=i7,__esvi8=i8,__esvd1=d1,dcrft(&__esvi1,x1, \
                     &__esvi2,x2,&__esvi3,&__esvi4,&__esvi5,&__esvi6, \
                     &__esvd1,x3,&__esvi7,x4,&__esvi8))
#endif
 
 
#ifdef  __ESVERR
#define dcft2(i1,x1,i2,i3,x2,i4,i5,x3,x4,i6,d1,x5,x6,x7,x8)           \
                                   (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                    __esvi4=i4,__esvi5=i5,__esvi6=i6, \
                                       __esvd1=d1,dcft2(&__esvi1,x1, \
                                       &__esvi2,&__esvi3,x2,&__esvi4, \
                                                      &__esvi5,x3,x4, \
                                     &__esvi6,&__esvd1,x5,x6,x7,x8))
#else
#define dcft2(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,d1,x3,i9,x4,i10)          \
                        (__esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                  __esvi8=i8,__esvi9=i9,__esvi10=i10, \
                                      __esvd1=d1,dcft2(&__esvi1,x1,  \
                                       &__esvi2,&__esvi3,x2,&__esvi4, \
                                          &__esvi5,&__esvi6,&__esvi7, \
                                    &__esvi8,&__esvd1,x3,&__esvi9,x4, \
                                                         &__esvi10))
#endif
 
#ifdef  __ESVERR
#define drcft2(i1,x1,i2,x2,i3,x3,x4,i4,d1,x5,x6,x7,x8)         \
                        (__esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                      __esvd1=d1,drcft2(&__esvi1,x1, \
                                         &__esvi2,x2,&__esvi3,x3,x4,  \
                              &__esvi4,&__esvd1,x5,x6,x7,x8))
#else
#define drcft2(i1,x1,i2,x2,i3,i4,i5,i6,d1,x3,i7,x4,i8)          \
                        (__esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                    __esvi8=i8,__esvd1=d1, \
                            drcft2(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                              &__esvi4,&__esvi5,&__esvi6,&__esvd1,x3, \
                                  &__esvi7,x4,&__esvi8))
#endif
 
#ifdef  __ESVERR
#define dcrft2(i1,x1,i2,x2,i3,x3,x4,i4,d1,x5,x6,x7,x8)         \
                                   (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                               __esvi4=i4,__esvd1=d1, \
                            dcrft2(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                          x3,x4,&__esvi4,&__esvd1,x5, \
                                                   x6,x7,x8))
#else
#define dcrft2(i1,x1,i2,x2,i3,i4,i5,i6,d1,x3,i7,x4,i8)          \
                        (__esvi1=i1,__esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                    __esvi8=i8,__esvd1=d1, \
                            dcrft2(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                              &__esvi4,&__esvi5,&__esvi6,&__esvd1,x3, \
                                  &__esvi7,x4,&__esvi8))
#endif
 
#ifdef  __ESVERR
#define dcft3(x1,i1,i2,x2,i3,i4,x3,x4,x5,i5,d1,x6,x7)               \
                                           (__esvi1=i1, __esvi2=i2, \
                                             __esvi3=i3,__esvi4=i4, \
                                             __esvi5=i5,__esvd1=d1, \
                           dcft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                                 &__esvi4,x3,x4,x5, \
                                         &__esvi5,&__esvd1,x6,x7))
#else
#define dcft3(x1,i1,i2,x2,i3,i4,i5,i6,i7,i8,d1,x3,i9)  (__esvi1=i1,  \
                                   __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                   __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                   __esvi8=i8,__esvi9=i9,__esvd1=d1, \
                            dcft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                    &__esvi8,&__esvd1,x3,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define drcft3(x1,i1,i2,x2,i3,i4,x3,x4,x5,i5,d1,x6,x7)   (__esvi1=i1, \
                                                          __esvi2=i2, \
                                               __esvi3=i3,__esvi4=i4, \
                                               __esvi5=i5,__esvd1=d1, \
                            drcft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                                   &__esvi4,x3,x4,x5, \
                                           &__esvi5,&__esvd1,x6,x7))
#else
#define drcft3(x1,i1,i2,x2,i3,i4,i5,i6,i7,i8,d1,x3,i9)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                    __esvi8=i8,__esvi9=i9,__esvd1=d1, \
                            drcft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                 &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                     &__esvi8,&__esvd1,x3,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define dcrft3(x1,i1,i2,x2,i3,i4,x3,x4,x5,i5,d1,x6,x7)   (__esvi1=i1, \
                                                          __esvi2=i2, \
                                               __esvi3=i3,__esvi4=i4, \
                                               __esvi5=i5,__esvd1=d1, \
                            dcrft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                                   &__esvi4,x3,x4,x5, \
                                           &__esvi5,&__esvd1,x6,x7))
#else
#define dcrft3(x1,i1,i2,x2,i3,i4,i5,i6,i7,i8,d1,x3,i9)   (__esvi1=i1, \
                                    __esvi2=i2,__esvi3=i3,__esvi4=i4, \
                                    __esvi5=i5,__esvi6=i6,__esvi7=i7, \
                                    __esvi8=i8,__esvi9=i9,__esvd1=d1, \
                            dcrft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                 &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                     &__esvi8,&__esvd1,x3,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define scosft(i1,x1,i2,i3,x2,i4,i5,x8,i7,f1,x4,x6,x5,x7)      (__esvi1=i1, \
                                           __esvf1=f1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                      scosft(&__esvi1,x1,&__esvi2,&__esvi3, \
                                           x2,&__esvi4,&__esvi5,x8,&__esvi7, \
                                                     &__esvf1,x4,x6,x5,x7))
#else
#define scosft(i1,x1,i2,i3,x2,i4,i5,i6,i7,f1,x4,i8,x5,i9)      (__esvi1=i1, \
                                           __esvf1=f1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                      scosft(&__esvi1,x1,&__esvi2,&__esvi3, \
                                     x2,&__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                         &__esvf1,x4,&__esvi8,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define scosf(i1,x1,i2,i3,x2,i4,i5,x8,i7,f1,x4,x6,x5,x7)      (__esvi1=i1, \
                                           __esvf1=f1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                      scosf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                           x2,&__esvi4,&__esvi5,x8,&__esvi7, \
                                                     &__esvf1,x4,x6,x5,x7))
#else
#define scosf(i1,x1,i2,i3,x2,i4,i5,i6,i7,f1,x4,i8,x5,i9)      (__esvi1=i1, \
                                           __esvf1=f1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                      scosf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                     x2,&__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                         &__esvf1,x4,&__esvi8,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define dcosf(i1,x1,i2,i3,x2,i4,i5,x8,i7,d1,x4,x6,x5,x7)      (__esvi1=i1, \
                                           __esvd1=d1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                      dcosf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                           x2,&__esvi4,&__esvi5,x8,&__esvi7, \
                                                     &__esvd1,x4,x6,x5,x7))
#else
#define dcosf(i1,x1,i2,i3,x2,i4,i5,i6,i7,d1,x4,i8,x5,i9)      (__esvi1=i1, \
                                           __esvd1=d1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                      dcosf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                     x2,&__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                         &__esvd1,x4,&__esvi8,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define ssinf(i1,x1,i2,i3,x2,i4,i5,x8,i7,f1,x4,x6,x5,x7)      (__esvi1=i1, \
                                           __esvf1=f1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                      ssinf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                           x2,&__esvi4,&__esvi5,x8,&__esvi7, \
                                                     &__esvf1,x4,x6,x5,x7))
#else
#define ssinf(i1,x1,i2,i3,x2,i4,i5,i6,i7,f1,x4,i8,x5,i9)      (__esvi1=i1, \
                                           __esvf1=f1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                      ssinf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                     x2,&__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                         &__esvf1,x4,&__esvi8,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define dsinf(i1,x1,i2,i3,x2,i4,i5,x8,i7,d1,x4,x6,x5,x7)      (__esvi1=i1, \
                                           __esvd1=d1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                      dsinf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                           x2,&__esvi4,&__esvi5,x8,&__esvi7, \
                                                     &__esvd1,x4,x6,x5,x7))
#else
#define dsinf(i1,x1,i2,i3,x2,i4,i5,i6,i7,d1,x4,i8,x5,i9)      (__esvi1=i1, \
                                           __esvd1=d1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                      dsinf(&__esvi1,x1,&__esvi2,&__esvi3, \
                                     x2,&__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                         &__esvd1,x4,&__esvi8,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define scft2(i1,x1,i2,i3,x2,i4,i5,x7,x8,i8,f1,x3,x5,x4,x6)    (__esvi1=i1, \
                                                      __esvi2=i2,__esvi5=i5, \
                                           __esvi3=i3,__esvi4=i4,__esvi8=i8, \
                                              __esvf1=f1,scft2(&__esvi1,x1, \
                                              &__esvi2,&__esvi3,x2,&__esvi4, \
                                                             &__esvi5,x7,x8, \
                                            &__esvi8,&__esvf1,x3,x5,x4,x6))
#else
#define scft2(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,f1,x3,i9,x4,i10)   (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                         __esvi8=i8,__esvi9=i9,__esvi10=i10, \
                                             __esvf1=f1,scft2(&__esvi1,x1,  \
                                              &__esvi2,&__esvi3,x2,&__esvi4, \
                                                 &__esvi5,&__esvi6,&__esvi7, \
                                           &__esvi8,&__esvf1,x3,&__esvi9,x4, \
                                                                &__esvi10))
#endif
 
#ifdef  __ESVERR
#define srcft2(i1,x1,i2,x2,i3,x9,x10,i6,f1,x3,x6,x4,x7,x5,x8)  (__esvi1=i1, \
                                           __esvi2=i2,__esvi6=i6,__esvi3=i3, \
                                             __esvf1=f1,srcft2(&__esvi1,x1, \
                                                &__esvi2,x2,&__esvi3,x9,x10, \
                                      &__esvi6,&__esvf1,x3,x6,x4,x7,x5,x8))
#else
#define srcft2(i1,x1,i2,x2,i3,i4,i5,i6,f1,x3,i7,x4,i8,x5,i9)   (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                           __esvi8=i8,__esvi9=i9,__esvf1=f1, \
                                   srcft2(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                     &__esvi4,&__esvi5,&__esvi6,&__esvf1,x3, \
                                         &__esvi7,x4,&__esvi8,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define scrft2(i1,x1,i2,x2,i3,x9,x10,i6,f1,x3,x6,x4,x7,x5,x8)  (__esvi1=i1, \
                                                      __esvi2=i2,__esvi6=i6, \
                                                      __esvi3=i3,__esvf1=f1, \
                                   scrft2(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                                x9,x10,&__esvi6,&__esvf1,x3, \
                                                           x6,x4,x7,x5,x8))
#else
#define scrft2(i1,x1,i2,x2,i3,i4,i5,i6,f1,x3,i7,x4,i8,x5,i9)   (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                           __esvi8=i8,__esvi9=i9,__esvf1=f1, \
                                   scrft2(&__esvi1,x1,&__esvi2,x2,&__esvi3, \
                                     &__esvi4,&__esvi5,&__esvi6,&__esvf1,x3, \
                                         &__esvi7,x4,&__esvi8,x5,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define scft3(x1,i1,i2,x2,i3,i4,x5,x6,x7,i8,f1,x3,x4)          (__esvi1=i1, \
                                                                 __esvi2=i2, \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi8=i8,__esvf1=f1, \
                                    scft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                                          &__esvi4,x5,x6,x7, \
                                                  &__esvi8,&__esvf1,x3,x4))
#else
#define scft3(x1,i1,i2,x2,i3,i4,i5,i6,i7,i8,f1,x3,i9)         (__esvi1=i1,  \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                           __esvi8=i8,__esvi9=i9,__esvf1=f1, \
                                    scft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                        &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                            &__esvi8,&__esvf1,x3,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define srcft3(x1,i1,i2,x2,i3,i4,x5,x6,x7,i8,f1,x3,x4)         (__esvi1=i1, \
                                                                 __esvi2=i2, \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi8=i8,__esvf1=f1, \
                                   srcft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                                          &__esvi4,x5,x6,x7, \
                                                  &__esvi8,&__esvf1,x3,x4))
#else
#define srcft3(x1,i1,i2,x2,i3,i4,i5,i6,i7,i8,f1,x3,i9)         (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                           __esvi8=i8,__esvi9=i9,__esvf1=f1, \
                                   srcft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                        &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                            &__esvi8,&__esvf1,x3,&__esvi9))
#endif
 
#ifdef  __ESVERR
#define scrft3(x1,i1,i2,x2,i3,i4,x5,x6,x7,i8,f1,x3,x4)         (__esvi1=i1, \
                                                                 __esvi2=i2, \
                                                      __esvi3=i3,__esvi4=i4, \
                                                      __esvi8=i8,__esvf1=f1, \
                                   scrft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                                          &__esvi4,x5,x6,x7, \
                                                  &__esvi8,&__esvf1,x3,x4))
#else
#define scrft3(x1,i1,i2,x2,i3,i4,i5,i6,i7,i8,f1,x3,i9)         (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                           __esvi8=i8,__esvi9=i9,__esvf1=f1, \
                                   scrft3(x1,&__esvi1,&__esvi2,x2,&__esvi3, \
                                        &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                            &__esvi8,&__esvf1,x3,&__esvi9))
#endif
 
 
/*  Convulutions/Correlation Subroutines   */
 
 
#ifdef  __ESVERR
#define scon(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,x6,x5,x7)       \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                                  __esvi10=i10,__esvi11=i11, \
                                                 scon(&__esvi1,x1,&__esvi2, \
                                  x2,&__esvi3,&__esvi4,x3,&__esvi5,&__esvi6, \
                                       &__esvi7,&__esvi8,&__esvi9,&__esvi10, \
                                                    &__esvi11,x4,x6,x5,x7))
#else
#define scon(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,i12,x5,i13)     \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                     __esvi10=i10,__esvi11=i11,__esvi12=i12, \
                                    __esvi13=i13,scon(&__esvi1,x1,&__esvi2, \
                                  x2,&__esvi3,&__esvi4,x3,&__esvi5,&__esvi6, \
                                       &__esvi7,&__esvi8,&__esvi9,&__esvi10, \
                                      &__esvi11,x4,&__esvi12,x5,&__esvi13))
#endif
 
#ifdef  __ESVERR
#define scor(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,x6,x5,x7)       \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                                  __esvi10=i10,__esvi11=i11, \
                                                 scor(&__esvi1,x1,&__esvi2, \
                                  x2,&__esvi3,&__esvi4,x3,&__esvi5,&__esvi6, \
                                       &__esvi7,&__esvi8,&__esvi9,&__esvi10, \
                                                    &__esvi11,x4,x6,x5,x7))
#else
#define scor(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,i12,x5,i13)     \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                     __esvi10=i10,__esvi11=i11,__esvi12=i12, \
                                    __esvi13=i13,scor(&__esvi1,x1,&__esvi2, \
                                  x2,&__esvi3,&__esvi4,x3,&__esvi5,&__esvi6, \
                                       &__esvi7,&__esvi8,&__esvi9,&__esvi10, \
                                      &__esvi11,x4,&__esvi12,x5,&__esvi13))
#endif
 
 
 
#define scond(x1,i1,x2,i2,x3,i3,i4,i5,i6,i7)        (__esvi1=i1,__esvi2=i2, \
                                           __esvi7=i7,__esvi5=i5,__esvi6=i6, \
                                            __esvi3=i3,__esvi4=i4,scond(x1, \
                                           &__esvi1,x2,&__esvi2,x3,&__esvi3, \
                                                          &__esvi4,&__esvi5, \
                                                        &__esvi6,&__esvi7))
 
 
#define scord(x1,i1,x2,i2,x3,i3,i4,i5,i6,i7)        (__esvi1=i1,__esvi2=i2, \
                                           __esvi7=i7,__esvi5=i5,__esvi6=i6, \
                                            __esvi3=i3,__esvi4=i4,scord(x1, \
                                           &__esvi1,x2,&__esvi2,x3,&__esvi3, \
                                      &__esvi4,&__esvi5,&__esvi6,&__esvi7))
 
#ifdef  __ESVERR
#define sconf(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,x6,x5,x7)      \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                                  __esvi10=i10,__esvi11=i11, \
                                                sconf(&__esvi1,x1,&__esvi2, \
                                           x2,&__esvi3,&__esvi4,x3,&__esvi5, \
                                        &__esvi6,&__esvi7,&__esvi8,&__esvi9, \
                                          &__esvi10,&__esvi11,x4,x6,x5,x7))
#else
#define sconf(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,i12,x5,i13)    \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                     __esvi10=i10,__esvi11=i11,__esvi12=i12, \
                                   __esvi13=i13,sconf(&__esvi1,x1,&__esvi2, \
                                           x2,&__esvi3,&__esvi4,x3,&__esvi5, \
                                        &__esvi6,&__esvi7,&__esvi8,&__esvi9, \
                                        &__esvi10,&__esvi11,x4,&__esvi12,x5, \
                                                                &__esvi13))
#endif
 
#ifdef   __ESVERR
#define scorf(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,x6,x5,x7)      \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                                  __esvi10=i10,__esvi11=i11, \
                                                scorf(&__esvi1,x1,&__esvi2, \
                                  x2,&__esvi3,&__esvi4,x3,&__esvi5,&__esvi6, \
                                       &__esvi7,&__esvi8,&__esvi9,&__esvi10, \
                                                    &__esvi11,x4,x6,x5,x7))
#else
#define scorf(i1,x1,i2,x2,i3,i4,x3,i5,i6,i7,i8,i9,i10,i11,x4,i12,x5,i13)    \
                                          (__esvi1=i1,__esvi2=i2,__esvi7=i7, \
                                           __esvi5=i5,__esvi6=i6,__esvi3=i3, \
                                           __esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                     __esvi10=i10,__esvi11=i11,__esvi12=i12, \
                                   __esvi13=i13,scorf(&__esvi1,x1,&__esvi2, \
                                  x2,&__esvi3,&__esvi4,x3,&__esvi5,&__esvi6, \
                                       &__esvi7,&__esvi8,&__esvi9,&__esvi10, \
                                      &__esvi11,x4,&__esvi12,x5,&__esvi13))
#endif
 
 
#define sdcon(x1,i1,x2,i2,x3,i3,i4,i5,i6,i7,i8)     (__esvi1=i1,__esvi8=i8, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                          sdcon(x1,&__esvi1,x2,&__esvi2,x3, \
                                        &__esvi3,&__esvi4,&__esvi5,&__esvi6, \
                                                        &__esvi7,&__esvi8))
 
 
#define ddcon(x1,i1,x2,i2,x3,i3,i4,i5,i6,i7,i8)     (__esvi1=i1,__esvi8=i8, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                          ddcon(x1,&__esvi1,x2,&__esvi2,x3, \
                                        &__esvi3,&__esvi4,&__esvi5,&__esvi6, \
                                                        &__esvi7,&__esvi8))
 
#define sdcor(x1,i1,x2,i2,x3,i3,i4,i5,i6,i7,i8)     (__esvi1=i1,__esvi8=i8, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                          sdcor(x1,&__esvi1,x2,&__esvi2,x3, \
                                        &__esvi3,&__esvi4,&__esvi5,&__esvi6, \
                                                        &__esvi7,&__esvi8))
 
#define ddcor(x1,i1,x2,i2,x3,i3,i4,i5,i6,i7,i8)     (__esvi1=i1,__esvi8=i8, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                          ddcor(x1,&__esvi1,x2,&__esvi2,x3, \
                                        &__esvi3,&__esvi4,&__esvi5,&__esvi6, \
                                                        &__esvi7,&__esvi8))
 
#ifdef  __ESVERR
#define sacor(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,x3,x5,x4,x6)       (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                __esvi6=i6,__esvi3=i3,__esvi4=i4,__esvi8=i8, \
                                    sacor(&__esvi1,x1,&__esvi2,&__esvi3,x2, \
                                        &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                                     &__esvi8,x3,x5,x4,x6))
#else
#define sacor(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,x3,i9,x4,i10)     (__esvi1=i1,  \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi6=i6,__esvi3=i3,__esvi4=i4, \
                                         __esvi8=i8,__esvi9=i9,__esvi10=i10, \
                                    sacor(&__esvi1,x1,&__esvi2,&__esvi3,x2, \
                                        &__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                        &__esvi8,x3,&__esvi9,x4,&__esvi10))
#endif
 
#ifdef  __ESVERR
#define sacorf(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,x3,x5,x4,x6)      (__esvi1=i1, \
                                __esvi2=i2,__esvi7=i7,__esvi5=i5,__esvi6=i6, \
                                           __esvi3=i3,__esvi4=i4,__esvi8=i8, \
                                               sacorf(&__esvi1,x1,&__esvi2, \
                                     &__esvi3,x2,&__esvi4,&__esvi5,&__esvi6, \
                                            &__esvi7,&__esvi8,x3,x5,x4,x6))
#else
#define sacorf(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,x3,i9,x4,i10)    (__esvi1=i1,  \
                                __esvi2=i2,__esvi7=i7,__esvi5=i5,__esvi6=i6, \
                                __esvi3=i3,__esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                  __esvi10=i10,sacorf(&__esvi1,x1,&__esvi2, \
                                     &__esvi3,x2,&__esvi4,&__esvi5,&__esvi6, \
                               &__esvi7,&__esvi8,x3,&__esvi9,x4,&__esvi10))
#endif
 
 
 
/*  Related Computations  Subroutines   */
 
 
 
#define spoly(x1,i1,i2,x2,i3,x3,i4,i5)   (__esvi1=i1,__esvi2=i2,__esvi5=i5, \
                                            __esvi3=i3,__esvi4=i4,spoly(x1, \
                                           &__esvi1,&__esvi2,x2,&__esvi3,x3, \
                                                        &__esvi4,&__esvi5))
 
#define dpoly(x1,i1,i2,x2,i3,x3,i4,i5)   (__esvi1=i1,__esvi2=i2,__esvi5=i5, \
                                            __esvi3=i3,__esvi4=i4,dpoly(x1, \
                                           &__esvi1,&__esvi2,x2,&__esvi3,x3, \
                                                        &__esvi4,&__esvi5))
 
 
#define sizc(x1,i1,i2,i3,x2)     (__esvi1=i1,__esvi2=i2,__esvi3=i3,sizc(x1, \
                                            &__esvi1,&__esvi2,&__esvi3,x2))
 
 
#define dizc(x1,i1,i2,i3,x2)     (__esvi1=i1,__esvi2=i2,__esvi3=i3,dizc(x1, \
                                            &__esvi1,&__esvi2,&__esvi3,x2))
 
 
#define strec(f1,x1,i1,x2,i2,x3,i3,i4,i5) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                        __esvi4=i4,__esvi5=i5,__esvf1=f1,strec(&__esvf1,x1, \
                                   &__esvi1,x2,&__esvi2,x3,&__esvi3,&__esvi4, \
                                                                  &__esvi5))
 
 
#define dtrec(d1,x1,i1,x2,i2,x3,i3,i4,i5) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                        __esvi4=i4,__esvi5=i5,__esvd1=d1,dtrec(&__esvd1,x1, \
                                   &__esvi1,x2,&__esvi2,x3,&__esvi3,&__esvi4, \
                                                                  &__esvi5))
 
 
#define sqint(f1,f2,f3,x1,i1,i2,x2,i3,x3,i4,i5)     (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                           __esvf1=f1,__esvf2=f2,__esvf3=f3, \
                              sqint(&__esvf1,&__esvf2,&__esvf3,x1,&__esvi1, \
                                &__esvi2,x2,&__esvi3,x3,&__esvi4,&__esvi5))
 
#define dqint(d1,d2,d3,x1,i1,i2,x2,i3,x3,i4,i5)     (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                           __esvd1=d1,__esvd2=d2,__esvd3=d3, \
                              dqint(&__esvd1,&__esvd2,&__esvd3,x1,&__esvi1, \
                                &__esvi2,x2,&__esvi3,x3,&__esvi4,&__esvi5))
 
#ifdef  __ESVERR
#define swlev(x1,i1,x2,i2,x3,i3,i4,x4,x5)           (__esvi1=i1,__esvi2=i2, \
                                                      __esvi3=i3,__esvi4=i4, \
                                 swlev(x1,&__esvi1,x2,&__esvi2,x3,&__esvi3, \
                                                           &__esvi4,x4,x5))
#else
#define swlev(x1,i1,x2,i2,x3,i3,i4,x4,i5)           (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                 swlev(x1,&__esvi1,x2,&__esvi2,x3,&__esvi3, \
                                                     &__esvi4,x4,&__esvi5))
#endif
 
#ifdef  __ESVERR
#define dwlev(x1,i1,x2,i2,x3,i3,i4,x4,x5)           (__esvi1=i1,__esvi2=i2, \
                                                      __esvi3=i3,__esvi4=i4, \
                                          dwlev(x1,&__esvi1,x2,&__esvi2,x3, \
                                                  &__esvi3,&__esvi4,x4,x5))
#else
#define dwlev(x1,i1,x2,i2,x3,i3,i4,x4,i5)           (__esvi1=i1,__esvi2=i2, \
                                           __esvi5=i5,__esvi3=i3,__esvi4=i4, \
                                          dwlev(x1,&__esvi1,x2,&__esvi2,x3, \
                                            &__esvi3,&__esvi4,x4,&__esvi5))
#endif
 
 
/*  Sorting and Searching  Subroutines   */
 
 
 
#define isort(x1,i1,i2) (__esvi1=i1,__esvi2=i2,isort(x1,&__esvi1,&__esvi2))
#define ssort(x1,i1,i2) (__esvi1=i1,__esvi2=i2,ssort(x1,&__esvi1,&__esvi2))
#define dsort(x1,i1,i2) (__esvi1=i1,__esvi2=i2,dsort(x1,&__esvi1,&__esvi2))
 
#define isortx(x1,i1,i2,x2)                (__esvi1=i1,__esvi2=i2,isortx(x1, \
                                                      &__esvi1,&__esvi2,x2))
#define ssortx(x1,i1,i2,x2)                (__esvi1=i1,__esvi2=i2,ssortx(x1, \
                                                      &__esvi1,&__esvi2,x2))
#define dsortx(x1,i1,i2,x2)                (__esvi1=i1,__esvi2=i2,dsortx(x1, \
                                                      &__esvi1,&__esvi2,x2))
 
 
 
#define ibsrch(x1,i1,i2,x2,i3,i4,x3,x4,i5) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                  __esvi4=i4,__esvi5=i5,ibsrch(x1,&__esvi1, \
                                 &__esvi2,x2,&__esvi3,&__esvi4,x3,x4,&__esvi5))
 
#define sbsrch(x1,i1,i2,x2,i3,i4,x3,x4,i5) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                   __esvi4=i4,__esvi5=i5,sbsrch(x1,&__esvi1, \
                                  &__esvi2,x2,&__esvi3,&__esvi4,x3,x4,&__esvi5))
 
#define dbsrch(x1,i1,i2,x2,i3,i4,x3,x4,i5) (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                  __esvi4=i4,__esvi5=i5,dbsrch(x1,&__esvi1, \
                                  &__esvi2,x2,&__esvi3,&__esvi4,x3,x4,&__esvi5))
 
#define issrch(x1,i1,i2,x2,i3,i4,i5,x3) (__esvi1=i1,__esvi2=i2,__esvi5=i5, \
                                  __esvi3=i3,__esvi4=i4,issrch(x1,&__esvi1, \
                         &__esvi2,x2,&__esvi3,&__esvi4,&__esvi5,x3))
 
#define sssrch(x1,i1,i2,x2,i3,i4,i5,x3) (__esvi1=i1,__esvi2=i2,__esvi5=i5, \
                                   __esvi3=i3,__esvi4=i4,sssrch(x1,&__esvi1, \
                          &__esvi2,x2,&__esvi3,&__esvi4,&__esvi5,x3))
 
#define dssrch(x1,i1,i2,x2,i3,i4,i5,x3) (__esvi1=i1,__esvi2=i2,__esvi5=i5, \
                                   __esvi3=i3,__esvi4=i4,dssrch(x1,&__esvi1, \
                           &__esvi2,x2,&__esvi3,&__esvi4,&__esvi5,x3))
 
 
#define isorts(x1,i1,i2,x2,x3,i3)   (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                  isorts(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
#define ssorts(x1,i1,i2,x2,x3,i3)   (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                  ssorts(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
#define dsorts(x1,i1,i2,x2,x3,i3)   (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                                  dsorts(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
 
/*  INTERPOLATION  Subroutines    */
 
 
 
#define spint(x1,x2,i1,x3,x7,x5,x6,i3)      (__esvi1=i1,__esvi3=i3,spint(x1, \
                                          x2,&__esvi1,x3,x7,x5,x6,&__esvi3))
 
 
#define dpint(x1,x2,i1,x3,x7,x5,x6,i3)      (__esvi1=i1,__esvi3=i3,dpint(x1, \
                                          x2,&__esvi1,x3,x7,x5,x6,&__esvi3))
 
#ifdef  __ESVERR
#define stpint(x1,x2,i1,i2,x3,x4,i3,x5,x6)           (__esvi1=i1,__esvi2=i2, \
                                           __esvi3=i3,stpint(x1,x2,&__esvi1, \
                                             &__esvi2,x3,x4,&__esvi3,x5,x6))
#else
#define stpint(x1,x2,i1,i2,x3,x4,i3,x5,i4)           (__esvi1=i1,__esvi2=i2, \
                                __esvi3=i3,__esvi4=i4,stpint(x1,x2,&__esvi1, \
                                       &__esvi2,x3,x4,&__esvi3,x5,&__esvi4))
#endif
 
#ifdef  __ESVERR
#define dtpint(x1,x2,i1,i2,x3,x4,i3,x5,x6)           (__esvi1=i1,__esvi2=i2, \
                                          __esvi3=i3,dtpint(x1,x2,&__esvi1, \
                                            &__esvi2,x3,x4,&__esvi3,x5,x6))
#else
#define dtpint(x1,x2,i1,i2,x3,x4,i3,x5,i4)           (__esvi1=i1,__esvi2=i2, \
                               __esvi3=i3,__esvi4=i4,dtpint(x1,x2,&__esvi1, \
                                      &__esvi2,x3,x4,&__esvi3,x5,&__esvi4))
#endif
 
 
#define scsint(x1,x2,x3,i1,x7,x5,x6,i3)    (__esvi1=i1,__esvi3=i3,scsint(x1, \
                                          x2,x3,&__esvi1,x7,x5,x6,&__esvi3))
 
 
#define dcsint(x1,x2,x3,i1,x7,x5,x6,i3)    (__esvi1=i1,__esvi3=i3,dcsint(x1, \
                                          x2,x3,&__esvi1,x7,x5,x6,&__esvi3))
 
#ifdef  __ESVERR
#define scsin2(x1,x2,x3,i1,i2,i3,x4,x5,i4,i5,x6,i6,x7,x8)       (__esvi1=i1, \
                                 __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                        __esvi6=i6,scsin2(x1,x2,x3,&__esvi1, \
                                   &__esvi2,&__esvi3,x4,x5,&__esvi4,&__esvi5, \
                                                         x6,&__esvi6,x7,x8))
#else
#define scsin2(x1,x2,x3,i1,i2,i3,x4,x5,i4,i5,x6,i6,x7,i7)       (__esvi1=i1, \
                                 __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                             __esvi6=i6,__esvi7=i7,scsin2(x1,x2,x3,&__esvi1, \
                                   &__esvi2,&__esvi3,x4,x5,&__esvi4,&__esvi5, \
                                                   x6,&__esvi6,x7,&__esvi7))
#endif
 
#ifdef  __ESVERR
#define dcsin2(x1,x2,x3,i1,i2,i3,x4,x5,i4,i5,x6,i6,x7,x8)       (__esvi1=i1, \
                                 __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                                        __esvi6=i6,dcsin2(x1,x2,x3,&__esvi1, \
                                   &__esvi2,&__esvi3,x4,x5,&__esvi4,&__esvi5, \
                                                         x6,&__esvi6,x7,x8))
#else
#define dcsin2(x1,x2,x3,i1,i2,i3,x4,x5,i4,i5,x6,i6,x7,i7)       (__esvi1=i1, \
                                 __esvi2=i2,__esvi3=i3,__esvi4=i4,__esvi5=i5, \
                             __esvi6=i6,__esvi7=i7,dcsin2(x1,x2,x3,&__esvi1, \
                                   &__esvi2,&__esvi3,x4,x5,&__esvi4,&__esvi5, \
                                                   x6,&__esvi6,x7,&__esvi7))
#endif
 
 
 
/*  NUMERICAL QUADRATURE  Subroutines    */
 
 
 
#define sptnq(x1,x2,i1,x3,x4)     (__esvi1=i1,sptnq(x1,x2,&__esvi1,x3,x4))
 
#define dptnq(x1,x2,i1,x3,x4)     (__esvi1=i1,dptnq(x1,x2,&__esvi1,x3,x4))
 
 
#define sglnq(__ESVFP,f1,f2,i1)    (__esvi1=i1,__esvf1=f1,__esvf2=f2, \
                              sglnq(__ESVFP,&__esvf1,&__esvf2,&__esvi1))
 
 
#define dglnq(__ESVFP,d1,d2,i1)    (__esvi1=i1,__esvd1=d1,__esvd2=d2, \
                              dglnq(__ESVFP,&__esvd1,&__esvd2,&__esvi1))
 
 
#define sglgq(__ESVFP,f1,f2,i1)    (__esvi1=i1,__esvf1=f1,__esvf2=f2, \
                              sglgq(__ESVFP,&__esvf1,&__esvf2,&__esvi1))
 
#define dglgq(__ESVFP,d1,d2,i1)    (__esvi1=i1,__esvd1=d1,__esvd2=d2, \
                              dglgq(__ESVFP,&__esvd1,&__esvd2,&__esvi1))
 
#define sgraq(__ESVFP,f1,f2,i1)    (__esvi1=i1,__esvf1=f1,__esvf2=f2, \
                              sgraq(__ESVFP,&__esvf1,&__esvf2,&__esvi1))
 
#define dgraq(__ESVFP,d1,d2,i1)    (__esvi1=i1,__esvd1=d1,__esvd2=d2, \
                              dgraq(__ESVFP,&__esvd1,&__esvd2,&__esvi1))
 
#define sghmq(__ESVFP,f1,f2,i1)    (__esvi1=i1,__esvf1=f1,__esvf2=f2, \
                              sghmq(__ESVFP,&__esvf1,&__esvf2,&__esvi1))
 
#define dghmq(__ESVFP,d1,d2,i1)    (__esvi1=i1,__esvd1=d1,__esvd2=d2, \
                              dghmq(__ESVFP,&__esvd1,&__esvd2,&__esvi1))
 
#define sglnq2(__ESVFP,f1,f2,i1,f3,f4,i2,x1,i3)      (__esvi1=i1,__esvi2=i2, \
                                __esvi3=i3,__esvf1=f1,__esvf2=f2,__esvf3=f3, \
                                    __esvf4=f4,sglnq2(__ESVFP, \
                                        &__esvf1,&__esvf2,&__esvi1,&__esvf3, \
                                            &__esvf4,&__esvi2,x1,&__esvi3))
 
#define dglnq2(__ESVFP,d1,d2,i1,d3,d4,i2,x1,i3)      (__esvi1=i1,__esvi2=i2, \
                                __esvi3=i3,__esvd1=d1,__esvd2=d2,__esvd3=d3, \
                                     __esvd4=d4,dglnq2(__ESVFP, \
                                        &__esvd1,&__esvd2,&__esvi1,&__esvd3, \
                                            &__esvd4,&__esvi2,x1,&__esvi3))
 
 
/*  Random Number Generator  Subroutines    */
 
 
#define surand(x1,i1,x2)               (__esvi1=i1,surand(x1,&__esvi1,x2))
 
#define durand(x1,i1,x2)               (__esvi1=i1,durand(x1,&__esvi1,x2))
 
#ifdef  __ESVERR
#define snrand(x1,i1,x2,x3,x4)                       (__esvi1=i1,snrand(x1, \
                                                        &__esvi1,x2,x3,x4))
#else
#define snrand(x1,i1,x2,x3,i2)            (__esvi1=i1,__esvi2=i2,snrand(x1, \
                                                  &__esvi1,x2,x3,&__esvi2))
#endif
 
#ifdef  __ESVERR
#define dnrand(x1,i1,x2,x3,x4)                       (__esvi1=i1,dnrand(x1, \
                                                        &__esvi1,x2,x3,x4))
#else
#define dnrand(x1,i1,x2,x3,i2)            (__esvi1=i1,__esvi2=i2,dnrand(x1, \
                                                  &__esvi1,x2,x3,&__esvi2))
#endif
 
#define surxor(x1,i1,x2,x3)  (__esvi1=i1,surxor(x1,&__esvi1,x2,x3))
#define durxor(x1,i1,x2,x3)  (__esvi1=i1,durxor(x1,&__esvi1,x2,x3))
 
 
/*  Utility Subroutines     */
 
 
#define dsrsm(i1,x1,x2,x3,i2,x7,x5,x6,i4)           (__esvi1=i1,__esvi2=i2, \
                                       __esvi4=i4,dsrsm(&__esvi1,x1,x2,x3, \
                                               &__esvi2,x7,x5,x6,&__esvi4))
 
#define stride(i1,i2,x,ch1,i4)    (__esvi1=i1,__esvi2=i2,__esvi4=i4,stride( \
                                         &__esvi1,&__esvi2,x,ch1,&__esvi4))
 
#define ivsset(i1)                           (__esvi1=i1,ivsset(&__esvi1))
 
#define einfo(i1,x1,x2)                 (__esvi1=i1,einfo(&__esvi1,x1,x2))
 
#ifdef __ESVERR
#define dgktrn(i1,x1,i2,x2,x3,i3,x4,x5,x6,x7)           (__esvi1=i1, \
                             __esvi2=i2,__esvi3=i3,dgktrn(&__esvi1, \
                          x1,&__esvi2,x2,x3,&__esvi3,x4,x5,x6,x7))
#else
#define dgktrn(i1,x1,i2,x2,x3,i3,x4,x5,x6,i4) (__esvi1=i1,__esvi2=i2, \
                             __esvi3=i3,__esvi4=i4,dgktrn(&__esvi1,  \
                      x1,&__esvi2,x2,x3,&__esvi3,x4,x5,x6,&__esvi4))
#endif
 
#ifdef __ESVERR
#define dsktrn(i1,x1,i2,x2,x3,x4,x5)         (__esvi1=i1,__esvi2=i2, \
                         dsktrn(&__esvi1,x1,&__esvi2,x2,x3,x4,x5))
#else
#define dsktrn(i1,x1,i2,x2,x3,x4,i3)         (__esvi1=i1,__esvi2=i2, \
                                       __esvi3=i3,dsktrn(&__esvi1,  \
                                    x1,&__esvi2,x2,x3,x4,&__esvi3))
#endif
 
#define ievops(i1)  (__esvi1=i1,ievops(&__esvi1))
 
 
 
 
 
/*  Parallel Processing Subroutines    */
 
 
#define dgemlp(x1,i1,ch1,x2,i2,ch2,x3,i3,i4,i5,i6)  (__esvi1=i1,__esvi2=i2, \
                                __esvi3=i3,__esvi4=i4,__esvi5=i5,__esvi6=i6, \
                                 dgemlp(x1,&__esvi1,ch1,x2,&__esvi2,ch2,x3, \
                                      &__esvi3,&__esvi4,&__esvi5,&__esvi6))
 
#ifdef  __ESVERR
#define dgefp(x1,i1,i2,x2,x3,x4)                    (__esvi1=i1,__esvi2=i2, \
                                    dgefp(x1,&__esvi1,&__esvi2,x2,x3,x4))
#else
#define dgefp(x1,i1,i2,x2,x3,i3)         (__esvi1=i1,__esvi2=i2,__esvi3=i3, \
                              dgefp(x1,&__esvi1,&__esvi2,x2,x3,&__esvi3))
#endif
 
#ifdef  __ESVERR
#define dppfp(x1,i1,x2,x3)           (__esvi1=i1,dppfp(x1,&__esvi1,x2,x3))
#else
#define dppfp(x1,i1,x2,i2)                 (__esvi1=i1,__esvi2=i2,dppfp(x1, \
                                                     &__esvi1,x2,&__esvi2))
#endif
 
#ifdef  __ESVERR
#define dskfsp(i1,x1,i2,x2,x3,x4,x5,x7,x6,i4,i5)    (__esvi1=i1,__esvi2=i2, \
                          __esvi4=i4,__esvi5=i5,dskfsp(&__esvi1,x1,&__esvi2, \
                                      x2,x3,x4,x5,x7,x6,&__esvi4,&__esvi5))
#else
#define dskfsp(i1,x1,i2,x2,x3,x4,x5,i3,x6,i4,i5)    (__esvi1=i1,__esvi2=i2, \
               __esvi3=i3,__esvi4=i4,__esvi5=i5,dskfsp(&__esvi1,x1,&__esvi2, \
                                x2,x3,x4,x5,&__esvi3,x6,&__esvi4,&__esvi5))
#endif
 
#ifdef __ESVERR
#define dgkfsp(i1,x1,i2,x2,x3,i3,x4,x5,x6,x7,x9,x8,i5,i6)      (__esvi1=i1, \
                                           __esvi5=i5,__esvi2=i2,__esvi6=i6, \
                                     __esvi3=i3,dgkfsp(&__esvi1,x1,&__esvi2, \
                       x2,x3,&__esvi3,x4,x5,x6,x7,x9,x8,&__esvi5,&__esvi6))
#else
#define dgkfsp(i1,x1,i2,x2,x3,i3,x4,x5,x6,x7,i4,x8,i5,i6)      (__esvi1=i1, \
                                           __esvi5=i5,__esvi2=i2,__esvi6=i6, \
                         __esvi3=i3,__esvi4=i4,dgkfsp(&__esvi1,x1,&__esvi2, \
                 x2,x3,&__esvi3,x4,x5,x6,x7,&__esvi4,x8,&__esvi5,&__esvi6))
#endif
 
#ifdef  __ESVERR
#define scftp(i1,x1,i2,i3,x2,i4,i5,x7,i7,i8,f1,x3,x5,x4,x6)    (__esvi1=i1, \
                                           __esvi2=i2,__esvi7=i7,__esvi5=i5, \
                                           __esvi3=i3,__esvi4=i4,__esvi8=i8, \
                                              __esvf1=f1,scftp(&__esvi1,x1, \
                                     &__esvi2,&__esvi3,x2,&__esvi4,&__esvi5, \
                                           x7,&__esvi7,&__esvi8,&__esvf1,x3, \
                                                                 x5,x4,x6))
#else
#define scftp(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,f1,x3,i9,x4,i10)  (__esvi1=i1,  \
                                __esvi2=i2,__esvi7=i7,__esvi5=i5,__esvi6=i6, \
                                __esvi3=i3,__esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                 __esvi10=i10,__esvf1=f1,scftp(&__esvi1,x1, \
                                     &__esvi2,&__esvi3,x2,&__esvi4,&__esvi5, \
                                     &__esvi6,&__esvi7,&__esvi8,&__esvf1,x3, \
                                                    &__esvi9,x4,&__esvi10))
#endif
 
 
#ifdef  __ESVERR
#define scft2p(i1,x1,i2,i3,x2,i4,i5,x7,x8,i8,f1,x3,x5,x4,x6)   (__esvi1=i1, \
                                                      __esvi2=i2,__esvi5=i5, \
                                           __esvi3=i3,__esvi4=i4,__esvi8=i8, \
                                             __esvf1=f1,scft2p(&__esvi1,x1, \
                                     &__esvi2,&__esvi3,x2,&__esvi4,&__esvi5, \
                                                 x7,x8,&__esvi8,&__esvf1,x3, \
                                                                 x5,x4,x6))
#else
#define scft2p(i1,x1,i2,i3,x2,i4,i5,i6,i7,i8,f1,x3,i9,x4,i10)  (__esvi1=i1, \
                                __esvi2=i2,__esvi7=i7,__esvi5=i5,__esvi6=i6, \
                                __esvi3=i3,__esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                 __esvi10=i10,__esvf1=f1,scft2p(&__esvi1,x1, \
                                     &__esvi2,&__esvi3,x2,&__esvi4,&__esvi5, \
                                     &__esvi6,&__esvi7,&__esvi8,&__esvf1,x3, \
                                                    &__esvi9,x4,&__esvi10))
#endif
 
#ifdef  __ESVERR
#define scft3p(x1,i1,i2,x2,i3,i4,x5,x6,x7,i8,f1,x3,x4)         (__esvi1=i1, \
                                                                 __esvi2=i2, \
                                           __esvi3=i3,__esvi4=i4,__esvi8=i8, \
                                 __esvf1=f1,scft3p(x1,&__esvi1,&__esvi2,x2, \
                                                 &__esvi3,&__esvi4,x5,x6,x7, \
                                                  &__esvi8,&__esvf1,x3,x4))
#else
#define scft3p(x1,i1,i2,x2,i3,i4,i5,i6,i7,i8,f1,x3,i9)         (__esvi1=i1, \
                                __esvi2=i2,__esvi7=i7,__esvi5=i5,__esvi6=i6, \
                                __esvi3=i3,__esvi4=i4,__esvi8=i8,__esvi9=i9, \
                                 __esvf1=f1,scft3p(x1,&__esvi1,&__esvi2,x2, \
                               &__esvi3,&__esvi4,&__esvi5,&__esvi6,&__esvi7, \
                                            &__esvi8,&__esvf1,x3,&__esvi9))
#endif
 
 
/* The following three definitions are added for error recovery migartion
   among different platforms     */
 
#define errset_ errset
#define errsav_ errsav
#define errstr_ errstr
#endif
 
#else
 
/* C++ Path   */
 
#ifndef _essl
 
  #define _essl 1
extern "C++" {
/**********************************************************************
*    LICENSED MATERIALS - PROPERTY OF IBM                             *
*    THIS MODULE IS "RESTRICTED MATERIALS OF IBM"                     *
*                                                                     *
*    COPYRIGHT = 5765-042 (C) COPYRIGHT IBM CORP. 1991, 1994.         *
*    ALL RIGHTS RESERVED.                                             *
*                                                                     *
*    U.S. GOVERNMENT USERS RESTRICTED RIGHTS - USE, DUPLICATION       *
*    OR DISCLOSURE RESTRICTED BY GSA ADP SCHEDULE CONTRACT WITH       *
*    IBM CORP.                                                        *
*    SEE COPYRIGHT INSTRUCTIONS.                                      *
*                                                                     *
*    THE SOURCE CODE FOR THIS PROGRAM IS NOT PUBLISHED OR OTHERWISE   *
*    DIVESTED OF ITS TRADE SECRETS, IRRESPECTIVE OF WHAT HAS BEEN     *
*    DEPOSITED WITH THE U.S. COPYRIGHT OFFICE.                        *
 *                                                                    *
 *  Program name - <essl.h> header file                               *
 *  Descriptive name - ESSL V2.2 C++ language header file             *
 *                                                                    *
 *  Function : This file must be included in any C++ file             *
 *             containing ESSL calls in order for the ESSL calls      *
 *             to work as documented in the ESSL Guide and            *
 *             Reference (SC23-0526).                                 *
 *                                                                    *
 *  Change activity - Added V2.2 routines   (AM-2/93)                 *
 *                    Total = 425+16=441                              *
 *                    Added iessl()           3/94                    *
 *********************************************************************/
 
#include <iostream.h>
#include <complex.h>
 
#ifdef _ESVERR
#define _ESVI    int
#else
#define _ESVI    void
#endif
 
/*  Definition of complex data types */
 
#ifndef  _CMPLX
#define  _CMPLX 1
class cmplx
  {
   private:
     union { struct { float  _re,_im; } _data; double _esvalign;};
   public:
      cmplx() { _data._re = 0.0; _data._im = 0.0; }
      cmplx(float r, float i = 0.0) { _data._re = r; _data._im = i; }
      cmplx(cmplx &c) { _data._re = c._data._re;   \
                        _data._im = c._data._im; }  //copy constructor
      friend inline float sreal(const cmplx& a) { return a._data._re; }
      friend inline float simag(const cmplx& a) { return a._data._im; }
  };
#endif
 
cmplx esvctmp;
complex esvdtmp;
 
extern "FORTRAN" {
 
/*  Linear Algebra Subprograms  */
 
/*  Vector-Scalar Subprograms  */
int isamax(const int &,  float *, const int &);
int idamax(const int &, double *, const int &);
int icamax(const int &,  cmplx *, const int &);
int izamax(const int &, complex *, const int &);
 
int isamin(const int &, float *, const int &);
int idamin(const int &, double *, const int &);
 
int ismax(const int &, float *, const int &);
int idmax(const int &, double *, const int &);
 
int ismin(const int &, float *, const int &);
int idmin(const int &, double *, const int &);
 
int iessl(void);
 
float   sasum(const int &,  float *, const int &);
double  dasum(const int &, double *, const int &);
float  scasum(const int &,  cmplx *, const int &);
double dzasum(const int &, complex *, const int &);
 
void  saxpy(const int &, const  float &,  float *, const int &,  float *,
            const int &);
void  daxpy(const int &, const double &, double *, const int &, double *,
            const int &);
void  caxpy(const int &, const  cmplx &,  cmplx *, const int &,  cmplx *,
            const int &);
void  zaxpy(const int &, const complex &, complex *, const int &, complex *,
            const int &);
 
void   scopy(const int &,  float *, const int &,  float *, const int &);
void   dcopy(const int &, double *, const int &, double *, const int &);
void   ccopy(const int &,  cmplx *, const int &,  cmplx *, const int &);
void   zcopy(const int &, complex *, const int &, complex *, const int &);
 
float  sdot(const int &, float *, const int &, float *, const int &);
double ddot(const int &, double *, const int &, double *, const int &);
 
void esvcdtu( const int &,   cmplx *, const int &,   cmplx *, const int &,
              cmplx &);
void esvzdtu( const int &, complex *, const int &, complex *, const int &,
            complex &);
void esvcdtc( const int &,  cmplx  *, const int &,  cmplx  *, const int &,
              cmplx &);
void esvzdtc( const int &, complex *, const int &, complex *, const int &,
            complex &);
#define cdotu(i1,x1,i2,x2,i3) (esvcdtu(i1,x1,i2,x2,i3,esvctmp),esvctmp)
#define zdotu(i1,x1,i2,x2,i3) (esvzdtu(i1,x1,i2,x2,i3,esvdtmp),esvdtmp)
#define cdotc(i1,x1,i2,x2,i3) (esvcdtc(i1,x1,i2,x2,i3,esvctmp),esvctmp)
#define zdotc(i1,x1,i2,x2,i3) (esvzdtc(i1,x1,i2,x2,i3,esvdtmp),esvdtmp)
 
void   snaxpy(const int &, const int &,  float*, const int &, void *,
              const int &, const int &,  void *, const int &, const int &);
void   dnaxpy(const int &, const int &, double*, const int &, void *,
              const int &, const int &,  void *, const int &, const int &);
 
void   sndot(const int &, const int &, float *, const int &, const int &,void *,
             const int &, const int &, void *, const int &, const int &);
void   dndot(const int &, const int &, double *,const int &, const int &,void *,
             const int &, const int &, void *, const int &, const int &);
 
float snrm2(const int &, float *, const int &);
double dnrm2(const int &, double *, const int &);
float  scnrm2(const int &,  cmplx *, const int &);
double dznrm2(const int &, complex *, const int &);
 
float  snorm2(const int &,  float *, const int &);
double dnorm2(const int &, double *, const int &);
float  cnorm2(const int &,  cmplx *, const int &);
double znorm2(const int &, complex *, const int &);
 
void   srotg( float *,  float *,  float *,  float *);
void   drotg(double *, double *, double *, double *);
void   crotg( cmplx *,  cmplx *,  float *,  cmplx *);
void   zrotg(complex *, complex *, double *, complex *);
 
void    srot(const int &,  float *, const int &,  float *, const int &,
             const float &, const float &);
void    drot(const int &, double *, const int &, double *, const int &,
             const double &, const double &);
void    crot(const int &,  cmplx *, const int &,  cmplx *, const int &,
             const  float &, const cmplx &);
void    zrot(const int &, complex *, const int &, complex *, const int &,
             const double &, const complex &);
void   csrot(const int &,  cmplx *, const int &,  cmplx *, const int &,
             const float &, const float &);
void   zdrot(const int &, complex *, const int &, complex *, const int &,
             const double &, const double &);
 
void    sscal(const int &, const   float &,   float *, const int &);
void    dscal(const int &, const  double &,  double *, const int &);
void    cscal(const int &, const   cmplx &,   cmplx *, const int &);
void    zscal(const int &, const complex &, complex *, const int &);
void   csscal(const int &, const   float &,   cmplx *, const int &);
void   zdscal(const int &, const  double &, complex *, const int &);
 
void   sswap(const int &,   float *, const int &,   float *, const int &);
void   dswap(const int &,  double *, const int &,  double *, const int &);
void   cswap(const int &,   cmplx *, const int &,   cmplx *, const int &);
void   zswap(const int &, complex *, const int &, complex *, const int &);
 
void    syax(const int &, const   float &,   float *, const int &,   float *,
             const int &);
void    dyax(const int &, const  double &,  double *, const int &,  double *,
             const int &);
void    cyax(const int &, const   cmplx &,   cmplx *, const int &,   cmplx *,
             const int &);
void    zyax(const int &, const complex &, complex *, const int &, complex *,
             const int &);
void   csyax(const int &, const   float &,   cmplx *, const int &,   cmplx *,
             const int &);
void   zdyax(const int &, const  double &, complex *, const int &, complex *,
             const int &);
 
void   szaxpy(const int &, const   float &,   float *, const int &,   float *,
              const int &, float *, const int &);
void   dzaxpy(const int &, const  double &,  double *, const int &,  double *,
              const int &, double *, const int &);
void   czaxpy(const int &, const   cmplx &,   cmplx *, const int &,   cmplx *,
              const int &, cmplx *, const int &);
void   zzaxpy(const int &, const complex &, complex *, const int &, complex *,
              const int &, complex *, const int &);
 
void svea(const int &,  float *, const int &,  float *, const int &,  float *,
          const int &);
void dvea(const int &, double *, const int &, double *, const int &, double *,
          const int &);
void cvea(const int &,  cmplx *, const int &,  cmplx *, const int &,  cmplx *,
          const int &);
void zvea(const int &, complex *, const int &, complex *, const int &,
          complex *, const int &);
 
void sves(const int &,  float *, const int &,  float *, const int &,  float *,
          const int &);
void dves(const int &, double *, const int &, double *, const int &, double *,
          const int &);
void cves(const int &,  cmplx *, const int &,  cmplx *, const int &,  cmplx *,
          const int &);
void zves(const int &, complex *, const int &, complex *, const int &,
          complex *, const int &);
 
void svem(const int &,  float *, const int &,  float *, const int &,  float *,
          const int &);
void dvem(const int &, double *, const int &, double *, const int &, double *,
          const int &);
void cvem(const int &,  cmplx *, const int &,  cmplx *, const int &,  cmplx *,
          const int &);
void zvem(const int &, complex *, const int &, complex *, const int &,
          complex *, const int &);
 
 
/*  Sparse Vector-Scalar Subroutines  */
 
void   ssctr(const int &,  float *,  int *,  float *);
void   dsctr(const int &, double *,  int *, double *);
void csctr(const int &,  cmplx *,  int *,  cmplx *);
void zsctr(const int &, complex *, int *, complex *);
 
void   sgthr(const int &,  float *,  float *,  int *);
void   dgthr(const int &, double *, double *,  int *);
void cgthr(const int &,  cmplx *,   cmplx *,  int *);
void zgthr(const int &, complex *, complex *, int *);
 
void   sgthrz(const int &,  float *,  float *,  int *);
void   dgthrz(const int &, double *, double *,  int *);
void cgthrz(const int &,  cmplx *,  cmplx *,  int *);
void zgthrz(const int &, complex *, complex *,  int *);
 
void   saxpyi(const int &, const  float &,  float *, int *,  float *);
void   daxpyi(const int &, const double &, double *, int *, double *);
void caxpyi(const int &, const   cmplx &,   cmplx *,  int *,   cmplx *);
void zaxpyi(const int &, const complex &, complex *,  int *, complex *);
 
float  sdoti(const int &,  float *, int *,  float *);
double ddoti(const int &, double *, int *, double *);
 
void  esvcdtui(const int &,  cmplx *,  int *,  cmplx *, cmplx &);
void  esvzdtui(const int &, complex *, int *, complex *, complex &);
void  esvcdtci(const int &,  cmplx *,  int *,  cmplx *, cmplx &);
void  esvzdtci(const int &, complex *, int *, complex *, complex &);
 
#define cdotui(i1,x1,x2,x3) (esvcdtui(i1,x1,x2,x3,esvctmp),esvctmp)
#define zdotui(i1,x1,x2,x3) (esvzdtui(i1,x1,x2,x3,esvdtmp),esvdtmp)
#define cdotci(i1,x1,x2,x3) (esvcdtci(i1,x1,x2,x3,esvctmp),esvctmp)
#define zdotci(i1,x1,x2,x3) (esvzdtci(i1,x1,x2,x3,esvdtmp),esvdtmp)
/*  Dense Matrix-Vector Subroutines  */
 
void   sgemv(char *, const int &, const int &, const  float &, void *,
             const int &,  float *, const int &, const  float &,  float *,
             const int &);
void   dgemv(char *, const int &, const int &, const double &, void *,
             const int &, double *, const int &, const double &, double *,
             const int &);
void   cgemv(char *, const int &, const int &, const  cmplx &, void *,
             const int &,  cmplx *, const int &, const  cmplx &, cmplx *,
             const int &);
void   zgemv(char *, const int &, const int &, const complex &, void *,
             const int &, complex *, const int &, const complex &, complex *,
             const int &);
 
void   sgemx(const int &, const int &, const float &, void *, const int &,
              float *, const int &,  float *, const int &);
void   dgemx(const int &, const int &,const double &,void *,const int &,
             double *, const int &, double *, const int &);
 
void   sgemtx(const int &,const int &,const  float &,void *,const int &,
               float *, const int &, float *, const int &);
void   dgemtx(const int &,const int &,const double &,void *,const int &,
             double *, const int &, double *, const int &);
 
#define sger  sger1
#define dger  dger1
void   sger1(const int &, const int &, const  float &,  float *, const int &,
              float *, const int &, void *, const int &);
void   dger1(const int &, const int &, const double &, double *, const int &,
             double *, const int &, void *, const int &);
void   cgeru(const int &, const int &, const cmplx &,  cmplx *, const int &,
               cmplx *, const int &, void *, const int &);
void   zgeru(const int &, const int &, const complex &, complex *, const int &,
             complex *, const int &, void *, const int &);
void   cgerc(const int &, const int &, const cmplx &,  cmplx *, const int &,
               cmplx *, const int &, void *, const int &);
void   zgerc(const int &, const int &, const complex &, complex *, const int &,
             complex *, const int &, void *, const int &);
 
void   sslmx(const int &, const float &, float *, float *, const int &,
              float *, const int &);
void   dslmx(const int &, const double &, double *, double *, const int &,
             double *,const int &);
 
void   sslr1(const int &, const  float &,  float *, const int &,  float *);
void   dslr1(const int &, const double &, double *, const int &, double *);
 
void   sslr2(const int &, const  float &,  float *, const int &,  float *,
             const int &, float *);
void   dslr2(const int &, const double &, double *, const int &, double *,
             const int &, double *);
void strmv(char *, char *, char *, const int &, void *, const int &,
           float *, const int &);
void dtrmv(char *, char *, char *, const int &, void *, const int &,
           void *, const int &);
void ctrmv(char *, char *, char *, const int &, void *, const int &,
           void *, const int &);
void ztrmv(char *, char *, char *, const int &, void *, const int &,
           void *, const int &);
 
void sspmv(char *, const int &,  const float &, float *,  float *, const int &,
           const float &, float *, const int &);
void dspmv(char *, const int &, const double &, double *, double *,const int &,
           const double &, double *, const int &);
void chpmv(char *, const int &,  const cmplx &,  cmplx *, cmplx *, const int &,
           const cmplx &, cmplx *, const int &);
void zhpmv(char *, const int &, const complex &, complex *, complex *,
           const int &, const complex &, complex *, const int &);
 
void ssymv(char *, const int &,  const float &, void *, const int &,  float *,
           const int &, const float &,  float *, const int &);
void dsymv(char *, const int &, const double &, void *, const int &, double *,
           const int &, const double &, double *, const int &);
void chemv(char *, const int &, const  cmplx &, void *, const int &,  cmplx *,
           const int &, const cmplx &,  cmplx *, const int &);
void zhemv(char *, const int &, const complex &, void *, const int &,complex *,
           const int &, const complex &, complex *, const int &);
 
void sspr(char *, const int &, const float &,  float *, const int &, float *);
void dspr(char *, const int &, const double &, double *, const int &,double *);
void chpr(char *, const int &, const  float &, cmplx *, const int &, cmplx *);
void zhpr(char *, const int &, const double &,complex *,const int &,complex *);
 
void ssyr(char *, const int &, const float &,  float *, const int &,  void *,
          const int &);
void dsyr(char *, const int &, const double &, double *, const int &,  void *,
          const int &);
void cher(char *, const int &, const float &,  cmplx *, const int &,  void *,
          const int &);
void zher(char *, const int &, const double &, complex *, const int &,  void *,
          const int &);
 
void sspr2(char *, const int &, const float &,  float *, const int &, float *,
           const int &, float *);
void dspr2(char *, const int &, const double &, double *, const int &,double *,
           const int &, double *);
void chpr2(char *, const int &, const cmplx &,  cmplx *, const int &, cmplx *,
           const int &, cmplx *);
void zhpr2(char *, const int &, const complex &, complex *, const int &,
           complex *, const int &, complex *);
 
void ssyr2(char *, const int &, const float &,  float *, const int &,  float *,
           const int &, void *, const int &);
void dsyr2(char *, const int &, const double &, double *, const int &,double *,
           const int &, void *, const int &);
void cher2(char *, const int &, const cmplx &,  cmplx *, const int &,  cmplx *,
           const int &, void *, const int &);
void zher2(char *, const int &, const complex &, complex *, const int &,
           complex *, const int &, void *, const int &);
 
void sgbmv(char *, const int &, const int &, const int &, const int &,
           const float &, void *, const int &, float *, const int &,
           const float &,  float *, const int &);
void dgbmv(char *, const int &, const int &, const int &, const int &,
           const double &, void *, const int &, double *, const int &,
           const double &, double *, const int &);
void cgbmv(char *, const int &, const int &, const int &, const int &,
           const cmplx &, void *, const int &, cmplx *, const int &,
           const cmplx &,  cmplx *, const int &);
void zgbmv(char *, const int &, const int &, const int &, const int &,
           const complex &, void *, const int &, complex *, const int &,
           const complex &, complex *, const int &);
 
void ssbmv(char *, const int &, const int &, const float &, void *,const int &,
           float *, const int &, const float &,  float *, const int &);
void dsbmv(char *, const int &, const int &, const double &,void *,const int &,
           double *, const int &, const double &, double *, const int &);
void chbmv(char *, const int &, const int &, const cmplx &, void *,const int &,
           cmplx *, const int &, const  cmplx &,  cmplx *, const int &);
void zhbmv(char *, const int &, const int &, const complex &, void *,
           const int &, complex *, const int &, const complex &, complex *,
           const int &);
 
void stbmv(char *, char *, char *, const int &, const int &, void *,
           const int &, float *, const int &);
void dtbmv(char *, char *, char *, const int &, const int &, void *,
           const int &,double *, const int &);
void ctbmv(char *, char *, char *, const int &, const int &, void *,
           const int &, cmplx *, const int &);
void ztbmv(char *, char *, char *, const int &, const int &, void *,
           const int &,complex *, const int &);
 
 
void stpmv(char *, char *, char *, const int &,  float *,  float *,
           const int &);
void dtpmv(char *, char *, char *, const int &, double *, double *,
           const int &);
void ctpmv(char *, char *, char *, const int &,  cmplx *,  cmplx *,
           const int &);
void ztpmv(char *, char *, char *, const int &, complex *, complex *,
           const int &);
 
/*  Sparse Matrix-Vector Subroutines  */
 
void   dsmmx(const int &, const int &, void *, void *, const int &, double *,
             double *);
 
_ESVI dsmtm(const int &, const int &, void *, void *, const int &, int &,
             int &, void *, void *, const int &, float *, const int &);
 
void   dsdmx(const int &, const int &, const int &, void *, const int &,
             char *, int *, double *, double *);
 
/*  Matrix Operation Subroutines  */
 
void   sgeadd(void *, const int &, char *, void *, const int &, char *,
              void *, const int &, const int &, const int &);
void   dgeadd(void *, const int &, char *, void *, const int &, char *,
              void *, const int &, const int &, const int &);
void   cgeadd(void *, const int &, char *, void *, const int &, char *,
              void *, const int &, const int &, const int &);
void   zgeadd(void *, const int &, char *, void *, const int &, char *,
              void *, const int &, const int &, const int &);
 
void   sgesub(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &);
void   dgesub(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &);
void   cgesub(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &);
void   zgesub(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &);
 
void   sgemul(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &, const int &);
void   dgemul(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &, const int &);
void   cgemul(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &, const int &);
void   zgemul(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &, const int &);
 
_ESVI sgemms(void *, const int &, char *, void *, const int &, char *, void *,
             const int &, const int &, const int &, const int &,  float *,
             const int &);
_ESVI dgemms(void *, const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &, const int &, double *,
              const int &);
_ESVI cgemms(void *, const int &, char *, void *, const int &, char *, void *,
             const int &, const int &, const int &, const int &,  float *,
             const int &);
_ESVI zgemms(void *, const int &, char *, void *, const int &, char *, void *,
             const int &, const int &, const int &, const int &, double *,
             const int &);
 
void   sgemm(char *, char *, const int &, const int &, const int &,
             const float &, void *, const int &, void *, const int &,
             const float &, void *, const int &);
void   dgemm(char *, char *, const int &, const int &, const int &,
             const double &, void *, const int &, void *, const int &,
             const double &, void *, const int &);
void   cgemm(char *, char *, const int &, const int &, const int &,
             const cmplx &, void *, const int &, void *, const int &,
             const cmplx &, void *, const int &);
void   zgemm(char *, char *, const int &, const int &, const int &,
             const complex &, void *, const int &, void *, const int &,
             const  complex &, void *, const int &);
 
void   ssyrk(char *, char *, const int &, const int &, const float &, void *,
             const int &, const float &, void *, const int &);
void   dsyrk(char *, char *, const int &, const int &, const double &, void *,
             const int &, const double &, void *, const int &);
 
void   sgetmi(void *, const int &, const int &);
void   dgetmi(void *, const int &, const int &);
void   cgetmi(void *, const int &, const int &);
void   zgetmi(void *, const int &, const int &);
 
void   sgetmo(void *, const int &, const int &, const int &, void *,
              const int &);
void   dgetmo(void *, const int &, const int &, const int &, void *,
              const int &);
void   cgetmo(void *, const int &, const int &, const int &, void *,
              const int &);
void   zgetmo(void *, const int &, const int &, const int &, void *,
              const int &);
void strmm(char *, char *, char *, char *, const int &, const int &,
           const  float &, void*, const int &, void *, const int &);
void dtrmm(char *, char *, char *, char *, const int &, const int &,
           const  double &, void*, const int &, void *, const int &);
void ctrmm(char *, char *, char *, char *, const int &, const int &,
           const cmplx &, void*, const int &, void *, const int &);
void ztrmm(char *, char *, char *, char *, const int &, const int &,
           const complex &, void*, const int &, void *, const int &);
void ssymm(char *, char *, const int &, const int &, const float &, void *,
           const int &, void *, const int &, const float &, void *,
           const int &);
void dsymm(char *, char *, const int &, const int &, const double &, void *,
           const int &, void *, const int &, const double &, void *,
           const int &);
void dsyr2k(char *, char *, const int &, const int &, const double &, void *,
            const int &, void *, const int &, const double &, void *,
            const int &);
void ssyr2k(char *, char *, const int &, const int &, const float &, void *,
            const int &, void *, const int &, const float &, void *,
            const int &);
 
 
void csyrk(char *, char *, const int &, const int &, const cmplx &, void *,
           const int &, const cmplx &, void *, const int &);
void zsyrk(char *, char *, const int &, const int &, const complex &, void *,
           const int &, const complex &, void *, const int &);
 
void cherk(char *, char *, const int &, const int &, const float &, void *,
           const int &, const float &, void *, const int &);
void zherk(char *, char *, const int &, const int &, const double &, void *,
           const int &, const double &, void *, const int &);
 
void csyr2k(char *, char *, const int &, const int &, const cmplx &, void *,
            const int &, void *, const int &, const cmplx &, void *,
            const int &);
void zsyr2k(char *, char *, const int &, const int &, const complex &, void *,
            const int &, void *, const int &, const complex &, void *,
            const int &);
void cher2k(char *, char *, const int &, const int &, const cmplx &, void *,
            const int &, void *, const int &, const float &, void *,
            const int &);
void zher2k(char *, char *, const int &, const int &, const complex &, void *,
            const int &, void *, const int &, const double &, void *,
            const int &);
 
 
void csymm(char *, char *, const int &, const int &, const cmplx &, void *,
           const int &, void *, const int &, const  cmplx &, void *,
           const int &);
void zsymm(char *, char *, const int &, const int &, const complex &, void *,
           const int &, void *, const int &, const complex &, void *,
           const int &);
void chemm(char *, char *, const int &, const int &, const cmplx &, void *,
           const int &, void *, const int &, const  cmplx &, void *,
           const int &);
void zhemm(char *, char *, const int &, const int &, const complex &, void *,
           const int &, void *, const int &, const complex &, void *,
           const int &);
 
/*  Dense Linear Algebraic Equation Subroutines  */
 
_ESVI sgef(void *, const int &, const int &, int *);
_ESVI dgef(void *, const int &, const int &, int *);
_ESVI cgef(void *, const int &, const int &, int *);
_ESVI zgef(void *, const int &, const int &, int *);
 
void   sges(void *, const int &, const int &, int *,  float *,
            const int &);
void   dges(void *, const int &, const int &, int *, double *,
            const int &);
void   cges(void *, const int &, const int &, int *,  cmplx *,
            const int &);
void   zges(void *, const int &, const int &, int *, complex *,
            const int &);
 
void sgesm(char *, void *, const int &, const int &, int *, void *,const int &,
           const int &);
void dgesm(char *, void *, const int &, const int &, int *, void *,const int &,
           const int &);
void cgesm(char *, void *, const int &, const int &, int *, void *,const int &,
           const int &);
void zgesm(char *, void *, const int &, const int &, int *, void *,const int &,
           const int &);
 
_ESVI sgefcd(void *, const int &, const int &, int *, const int &,
             float *, float *, float *, const int &);
_ESVI dgefcd(void *, const int &, const int &, int *, const int &,
             double *, double *, double *, const int &);
 
_ESVI sppf( float *, const int &, const int &);
_ESVI dppf(double *, const int &, const int &);
 
void   spps( float *, const int &,  float *, const int &);
void   dpps(double *, const int &, double *, const int &);
 
_ESVI sppfcd( float *, const int &, const int &, float *, float *, float *,
              const int &);
_ESVI dppfcd(double *,const int &, const int &, double *, double *, double *,
              const int &);
 
_ESVI sgeicd(void *, const int &, const int &, const int &, float &,
             float *,  float *, const int &);
_ESVI dgeicd(void *, const int &, const int &, const int &, double &,
              double *, double *, const int &);
 
_ESVI sppicd( float *, const int &, const int &,  float *,  float *,  float *,
               const int &);
_ESVI dppicd(double *, const int &, const int &, double *, double *, double *,
               const int &);
 
void   strsm(char *, char *, char *, char *, const int &, const int &,
             const float &, void *, const int &, void *, const int &);
void   dtrsm(char *, char *, char *, char *, const int &, const int &,
             const double &, void *, const int &, void *, const int &);
void ctrsm(char *, char *, char *, char *, const int &, const int &,
           const cmplx &, void *, const int &, void *, const int &);
void ztrsm(char *, char *, char *, char *, const int &, const int &,
           const complex &, void *, const int &, void *, const int &);
void strsv(char *, char *, char *, const int &, void *, const int &,
           void *, const int &);
void dtrsv(char *, char *, char *, const int &, void *, const int &,
           void *, const int &);
void ctrsv(char *, char *, char *, const int &, void *, const int &,
           void *, const int &);
void ztrsv(char *, char *, char *, const int &, void *, const int &,
           void *, const int &);
 
void spof(char *, void *, const int &, const int &);
void dpof(char *, void *, const int &, const int &);
void cpof(char *, void *, const int &, const int &);
void zpof(char *, void *, const int &, const int &);
 
void sposm(char *, void *, const int &, const int &, void *, const int &,
           const int &);
void dposm(char *, void *, const int &, const int &, void *, const int &,
           const int &);
void cposm(char *, void *, const int &, const int &, void *, const int &,
           const int &);
void zposm(char *, void *, const int &, const int &, void *, const int &,
           const int &);
 
_ESVI spofcd(char *, void *, const int &, const int &, const int &,  float *,
             float *, float *, const int &);
_ESVI dpofcd(char *, void *, const int &, const int &, const int &, double *,
             double *, double *, const int &);
 
_ESVI spoicd(char *, void *, const int &, const int &, const int &,  float *,
             float *, float *, const int &);
_ESVI dpoicd(char *, void *, const int &, const int &, const int &, double *,
             double *, double *, const int &);
 
_ESVI stri(char *, char *, void *, const int &, const int &);
_ESVI dtri(char *, char *, void *, const int &, const int &);
 
_ESVI  stpi(char *, char *,  float *, const int &);
_ESVI  dtpi(char *, char *, double *, const int &);
 
void stpsv(char *, char *, char *, const int &,  float *,  float *,
           const int &);
void dtpsv(char *, char *, char *, const int &, double *, double *,
           const int &);
void ctpsv(char *, char *, char *, const int &,  cmplx *,  cmplx *,
           const int &);
void ztpsv(char *, char *, char *, const int &, complex *, complex *,
           const int &);
 
void stbsv(char *, char *, char *, const int &, const int &, void *,
           const int &, float *, const int &);
void dtbsv(char *, char *, char *, const int &, const int &, void *,
           const int &,double *, const int &);
void ctbsv(char *, char *, char *, const int &, const int &, void *,
           const int &, cmplx *, const int &);
void ztbsv(char *, char *, char *, const int &, const int &, void *,
           const int &, complex *, const int &);
 
void cgtnpf(const int &,  cmplx *,  cmplx *,  cmplx *, const int &);
void zgtnpf(const int &, complex *, complex *, complex *, const int &);
 
void cgtnps(const int &,  cmplx *,  cmplx *,  cmplx *,  cmplx *);
void zgtnps(const int &, complex *, complex *, complex *, complex *);
 
_ESVI dsris(char *, char *, const int &, double *, int *, int *, double *,
            double *, int *, double *, double *, const int &, double *,
            const int &);
 
/*  Banded Linear Algebraic Equation Subroutines  */
 
_ESVI sgbf(void *, const int &, const int &, const int &, const int &, int *);
_ESVI dgbf(void *, const int &, const int &, const int &, const int &, int *);
 
void   sgbs(void *, const int &, const int &, const int &, const int &,
            int *,  float *);
void   dgbs(void *, const int &, const int &, const int &, const int &,
            int *, double *);
 
_ESVI spbf(void *, const int &, const int &, const int &);
_ESVI dpbf(void *, const int &, const int &, const int &);
 
void   spbs(void *, const int &, const int &, const int &,  float *);
void   dpbs(void *, const int &, const int &, const int &, double *);
 
_ESVI spbchf(void *, const int &, const int &, const int &);
_ESVI dpbchf(void *, const int &, const int &, const int &);
 
void   spbchs(void *, const int &, const int &, const int &,  float *);
void   dpbchs(void *, const int &, const int &, const int &, double *);
 
_ESVI sgtf(const int &,  float *,  float *,  float *,  float *, int *);
_ESVI dgtf(const int &, double *, double *, double *, double *, int *);
 
void   sgts(const int &,  float *,  float *,  float *,  float *, int *,
            float *);
void   dgts(const int &, double *, double *, double *, double *, int *,
            double *);
 
void   sgtnpf(const int &,  float *,  float *,  float *, const int &);
void   dgtnpf(const int &, double *, double *, double *, const int &);
 
void   sgtnps(const int &,  float *,  float *,  float *,  float *);
void   dgtnps(const int &, double *, double *, double *, double *);
 
void   sgtnp(const int &,  float *,  float *,  float *,  float *);
void   dgtnp(const int &, double *, double *, double *, double *);
void   cgtnp(const int &,  cmplx *,  cmplx *,  cmplx *,  cmplx *);
void   zgtnp(const int &, complex *, complex *, complex *, complex *);
 
void   sptf(const int &,  float *,  float *, const int &);
void   dptf(const int &, double *, double *, const int &);
 
void   spts(const int &,  float *,  float *,  float *);
void   dpts(const int &, double *, double *, double *);
 
/*  Sparse Linear Algebraic Equation Subroutines  */
 
_ESVI dgsf(const int &,const int &, const int &, double *, int *,
           int *, const int &, int *, double *, double *,
           double *, const int &);
 
_ESVI dgss(const int &, const int &, double *, int *, int *,
           const int &, double *, double *, const int &);
 
_ESVI dgkfs(const int &, double *, const int &, int *, double *,
            const int &, int *, int *, double *, double *,
            const int &, void *, const int &, const int &);
 
_ESVI dskfs(const int &, double *, const int &, int *, int *, double *,
            double *, const int &, void *, const int &, const int &);
 
_ESVI dsmcg(const int &, const int &, void *, void *, const int &, double *,
            double *, int *, double *, double *, const int &, double *,
            const int &);
 
_ESVI dsdcg(const int &, const int &, const int &, void *, const int &,
            int *, double *, double *, int *, double *, double *,
            const int &, double *, const int &);
 
_ESVI dsmgcg(const int &, const int &, void *, void *, const int &, double *,
             double *, int *, double *, double *, const int &, double *,
             const int &);
 
_ESVI dsdgcg(const int &, const int &, void *, const int &, int *, double *,
             double *, int *, double *, double *, const int &, double *,
             const int &);
 
/*  Linear Least Squares Subroutines  */
 
_ESVI sgesvf(const int &, void *, const int &, void *, const int &,
             const int &,  float *, const int &, const int &, float *,
             const int &);
_ESVI dgesvf(const int &, void *, const int &, void *, const int &,
             const int &, double *, const int &, const int &, double *,
             const int &);
 
void   sgesvs(void *, const int &, void *, const int &, const int &,  float *,
              void *, const int &, const int &, const int &, const  float &);
void   dgesvs(void *, const int &, void *, const int &, const int &, double *,
              void *, const int &, const int &, const int &, const double &);
 
_ESVI sgells(const int &, void *, const int &, void *, const int &, void *,
             const int &,  float *, const float &, const int &, const int &,
             const int &, int *,  float *, const int &);
_ESVI dgells(const int &, void *, const int &, void *, const int &, void *,
             const int &, double *, const  double &, const int &, const int &,
             const int &, int *, double *, const int &);
 
 
/*  Eigensystem Analysis Subroutines  */
 
 
_ESVI sgeev(const int &, void *, const int &,  cmplx *, void *, const int &,
            int *, const int &, float *, const int &);
_ESVI dgeev(const int &, void *, const int &, complex *, void *, const int &,
            int *, const int &, double *, const int &);
_ESVI cgeev(const int &, void *, const int &,  cmplx *, void *, const int &,
            int *, const int &, float *, const int &);
_ESVI zgeev(const int &, void *, const int &, complex *, void *, const int &,
            int *, const int &, double *, const int &);
 
#define sslev  sspev
#define dslev  dspev
#define chlev  chpev
#define zhlev  zhpev
_ESVI sspev(const int &,  float *,  float *, void *, const int &, const int &,
            float *, const int &);
_ESVI dspev(const int &, double *, double *, void *, const int &, const int &,
            double *, const int &);
_ESVI chpev(const int &,  cmplx *,  float *, void *, const int &, const int &,
            float *, const int &);
_ESVI zhpev(const int &, complex *, double *, void *, const int &, const int &,
            double *, const int &);
 
_ESVI sspsv(const int &,  float *,  float *, void *, const int &, const int &,
            const int &,  float *, const int &);
_ESVI dspsv(const int &, double *, double *, void *, const int &, const int &,
            const int &, double *, const int &);
_ESVI chpsv(const int &,  cmplx *,  float *, void *, const int &, const int &,
            const int &,  float *, const int &);
_ESVI zhpsv(const int &, complex *, double *, void *, const int &, const int &,
            const int &, double *, const int &);
 
_ESVI sgegv(const int &, void *, const int &, void *, const int &,  cmplx *,
            float *, void *, const int &, const int &,  float *, const int &);
_ESVI dgegv(const int &, void *, const int &, void *, const int &, complex *,
            double *, void *, const int &, const int &, double *, const int &);
 
_ESVI ssygv(const int &, void *, const int &, void *, const int &,  float *,
            void *, const int &, const int &, float *, const int &);
_ESVI dsygv(const int &, void *, const int &, void *, const int &, double *,
            void *, const int &, const int &, double *, const int &);
 
/*  Fourier Transform Subroutines  */
 
_ESVI scft(const int &, void *, const int &, const int &, void *, const int &,
           const int &, const int &, const int &, const int &, const float &,
           double *, const int &, double *, const int &);
 
_ESVI dcft(const int &, void *, const int &, const int &, void *, const int &,
           const int &, const int &, const int &, const int &, const double &,
           double *, const int &, double *, const int &);
 
_ESVI srcft(const int &, void *, const int &, void *, const int &, const int &,
            const int &, const int &, const float &, double *, const int &,
            double *, const int &, double *, const int &);
 
_ESVI drcft(const int &,void *, const int &, void *, const int &, const int &,
            const int &, const int &, const double &, double *, const int &,
            double *, const int &);
 
_ESVI scrft(const int &,void *, const int &, void *, const int &, const int &,
            const int &, const int &, const float &, double *, const int &,
            double *, const int &, double *, const int &);
 
_ESVI dcrft(const int &,void *, const int &, void *, const int &, const int &,
            const int &, const int &, const double &, double *, const int &,
            double *, const int &);
 
 
_ESVI dcft2(const int &, void *, const int &, const int &, void *, const int &,
            const int &, const int &, const int &, const int &, const double &,
            double *, const int &, double *, const int &);
 
_ESVI dcft3(void *, const int &, const int &, void *, const int &, const int &,
            const int &, const int &, const int &, const int &, const double &,
            double *, const int &);
 
_ESVI drcft2(const int &, void *, const int &, void *, const int &,
            const int &, const int &, const int &, const double &,
            double *, const int &, double *, const int &);
 
_ESVI dcrft2(const int &, void *, const int &, void *, const int &,
             const int &, const int &, const int &, const double &,
             double *, const int &, double *, const int &);
 
_ESVI drcft3(void *, const int &, const int &, void *, const int &,
             const int &, const int &, const int &, const int &, const int &,
             const double &, double *, const int &);
 
_ESVI dcrft3(void *, const int &, const int &, void *, const int &,
             const int &, const int &, const int &, const int &, const int &,
             const double &, double *, const int &);
 
_ESVI scosft(const int &, void *, const int &, const int &, void *, const int &,
             const int &, const int &, const int &, const float &, double *,
             const int &, double *, const int &);
 
_ESVI scosf(const int &, void *, const int &, const int &, void *, const int &,
             const int &, const int &, const int &, const float &, double *,
             const int &, double *, const int &);
 
_ESVI dcosf(const int &, void *, const int &, const int &, void *, const int &,
             const int &, const int &, const int &, const double &, double *,
             const int &, double *, const int &);
 
_ESVI ssinf(const int &, void *, const int &, const int &, void *, const int &,
             const int &, const int &, const int &, const float &, double *,
             const int &, double *, const int &);
 
_ESVI dsinf(const int &, void *, const int &, const int &, void *, const int &,
             const int &, const int &, const int &, const double &, double *,
             const int &, double *, const int &);
 
_ESVI scft2(const int &, void *, const int &, const int &, void *, const int &,
            const int &, const int &, const int &, const int &, const float &,
            double *, const int &, double *, const int &);
 
_ESVI srcft2(const int &,void *, const int &, void *, const int &, const int &,
             const int &, const int &, const float &, double *, const int &,
             double *, const int &, double *, const int &);
 
_ESVI scrft2(const int &,void *, const int &, void *, const int &, const int &,
             const int &, const int &, const float &, double *, const int &,
             double *, const int &, double *, const int &);
 
_ESVI scft3(void *, const int &, const int &, void *, const int &, const int &,
            const int &, const int &, const int &, const int &, const float &,
            double *, const int &);
 
_ESVI srcft3(void *, const int &, const int &, void *, const int &, const int &,
             const int &, const int &, const int &, const int &, const float &,
             double *, const int &);
 
_ESVI scrft3(void *, const int &, const int &, void *, const int &, const int &,
             const int &, const int &, const int &, const int &, const float &,
             double *, const int &);
 
/*  Convolutions/Correlations Subroutines  */
 
_ESVI scon(const int &, float *, const int &, void *, const int &, const int &,
           void *, const int &, const int &, const int &, const int &,
           const int &, const int &, const int &, double *, const int &,
           double *, const int &);
_ESVI scor(const int &, float *, const int &, void *, const int &, const int &,
           void *, const int &, const int &, const int &, const int &,
           const int &, const int &, const int &, double *, const int &,
           double *, const int &);
 
void   scond(float *, const int &, void *, const int &, void *, const int &,
             const int &, const int &, const int &, const int &);
void   scord(float *, const int &, void *, const int &, void *, const int &,
             const int &, const int &, const int &, const int &);
 
_ESVI sconf(const int &, float *, const int &, void *, const int &,const int &,
            void *, const int &, const int &, const int &, const int &,
            const int &, const int &, const int &, double *, const int &,
            double *, const int &);
_ESVI scorf(const int &, float *, const int &, void *, const int &,const int &,
            void *, const int &, const int &, const int &, const int &,
            const int &, const int &, const int &, double *, const int &,
            double *, const int &);
 
void   sdcon( float *, const int &, void *, const int &, void *, const int &,
              const int &, const int &, const int &, const int &, const int &);
void   ddcon(double *, const int &, void *, const int &, void *, const int &,
             const int &, const int &, const int &, const int &, const int &);
void   sdcor( float *, const int &, void *, const int &, void *, const int &,
             const int &, const int &, const int &, const int &, const int &);
void   ddcor(double *, const int &, void *, const int &, void *, const int &,
             const int &, const int &, const int &, const int &, const int &);
 
_ESVI sacor(const int &, void *, const int &, const int &, void *, const int &,
            const int &, const int &, const int &, const int &, double *,
            const int &, double *, const int &);
 
_ESVI sacorf(const int &, void *, const int &, const int &, void *,const int &,
             const int &, const int &, const int &, const int &, double *,
             const int &, double *, const int &);
 
/*  Related Computations Subroutines  */
 
void   spoly( float *, const int &, const int &,  float *, const int &,
             float *, const int &, const int &);
void   dpoly(double *, const int &, const int &, double *, const int &,
             double *, const int &, const int &);
 
void   sizc( float *, const int &, const int &, const int &, int *);
void   dizc(double *, const int &, const int &, const int &, int *);
 
void   strec(const float &,  float *, const int &,  float *, const int &,
             float *, const int &, const int &, const int &);
void   dtrec(const double &, double *, const int &, double *, const int &,
             double *, const int &, const int &, const int &);
 
_ESVI sqint(const float &, const float &, const float &, float *, const int &,
            const int &, float *, const int &, float *, const int &,
            const int &);
_ESVI dqint(const double &, const double &, const double &, double *,
            const int &, const int &, double *, const int &, double *,
            const int &, const int &);
 
_ESVI swlev( float *, const int &,  float *, const int &, float *, const int &,
            const int &, double *, const int &);
_ESVI dwlev(double *, const int &, double *, const int &, double *, const int &,
            const int &, double *, const int &);
 
/*  Sorting and Searching Subroutines  */
 
void   isort(   int *, const int &, const int &);
void   ssort( float *, const int &, const int &);
void   dsort(double *, const int &, const int &);
 
void   isortx(   int *, const int &, const int &, int *);
void   ssortx( float *, const int &, const int &, int *);
void   dsortx(double *, const int &, const int &, int *);
 
void   ibsrch( int *, const int &, const int &, int *, const int &,
              const int &, int *, int *, const int &);
void   sbsrch( float *, const int &, const int &,  float *, const int &,
               const int &, int *, int *, const int &);
void   dbsrch(double *, const int &, const int &, double *, const int &,
              const int &, int *, int *, const int &);
 
void   issrch( int *, const int &, const int &, int *, const int &,
              const int &, const int &, int *);
void   sssrch( float *, const int &, const int &,  float *, const int &,
              const int &, const int &, int *);
void   dssrch(double *, const int &, const int &, double *, const int &,
              const int &, const int &, int *);
 
void   isorts( int *, const int &, const int &,  int *, int *, const int &);
void   ssorts( float *, const int &, const int &, int *, float *, const int &);
void   dsorts(double *, const int &, const int &, int *, double *, const int &);
 
/*  Interpolation Subroutines  */
 
void   spint( float *,  float *, const int &,  float *, int *, float *,
              float *, const int &);
void   dpint(double *, double *, const int &, double *, int *, double *,
             double *, const int &);
 
_ESVI stpint( float *,  float *, const int &, const int &,  float *,  float *,
              const int &, float *, const int &);
_ESVI dtpint(double *, double *, const int &, const int &, double *, double *,
             const int &, double *, const int &);
 
void   scsint( float *,  float *, void *, const int &, int *,  float *,
               float *, const int &);
void   dcsint(double *, double *, void *, const int &, int *, double *,
              double *, const int &);
 
_ESVI scsin2( float *,  float *, void *, const int &, const int &, const int &,
              float *, float *, const int &, const int &, void *, const int &,
              float *, const int &);
_ESVI dcsin2(double *, double *, void *, const int &, const int &, const int &,
             double *, double *, const int &, const int &, void *, const int &,
             double *, const int &);
 
 
/*  Numerical Quadrature Subroutines  */
 
void   sptnq( float *,  float *, const int &,  float *,  float *);
void   dptnq(double *, double *, const int &, double *, double *);
 
float  sglnq(void *, const float &, const float &, const int &);
double dglnq(void *, const double &, const double &, const int &);
 
float  sglnq2(void *, const float &, const float &, const int &, const float &,
              const float &, const int &, void *, const int &);
double dglnq2(void *, const double &,const double &,const int &,const double &,
              const double &, const int &, void *, const int &);
 
float  sglgq(void *,  const float &,  const float &, const int &);
double dglgq(void *, const double &, const double &, const int &);
 
float  sgraq(void *,  const float &, const  float &, const int &);
double dgraq(void *, const double &, const double &, const int &);
 
float  sghmq(void *,  const float &,  const float &, const int &);
double dghmq(void *, const double &, const double &, const int &);
 
/*  Random Number Generation Subroutines  */
 
void   surand(double *, const int &,  float *);
void   durand(double *, const int &, double *);
 
_ESVI snrand(double *, const int &,  float *,  float *, const int &);
_ESVI dnrand(double *, const int &, double *, double *, const int &);
 
void surxor(int *, const int &,  float *,  float *);
void durxor(int *, const int &, double *, double *);
 
/*  Parallel Processing Subroutines  */
 
void   dgemlp(void *,const int &, char *, void *, const int &, char *, void *,
              const int &, const int &, const int &, const int &);
 
_ESVI dgefp(void *, const int &, const int &, int *, double *,
            const int &);
 
_ESVI dppfp(double *, const int &, double *, const int &);
 
_ESVI dgkfsp(const int &, double *, const int &, int *, double *,
            const int &, int *, int *, double *, double *,
            const int &, void *, const int &, const int &);
 
_ESVI dskfsp(const int &, double *, const int &, int *, int *, double *,
            double *, const int &, void *, const int &, const int &);
 
_ESVI scftp(const int &, void *, const int &, const int &, void *, const int &,
           const int &, const int &, const int &, const int &, const float &,
           double *, const int &, double *, const int &);
 
_ESVI scft2p(const int &, void *, const int &, const int &, void *, const int &,
            const int &, const int &, const int &, const int &, const float &,
            double *, const int &, double *, const int &);
 
_ESVI scft3p(void *, const int &, const int &, void *, const int &, const int &,
            const int &, const int &, const int &, const int &, const float &,
            double *, const int &);
 
 
/*  Utility Subroutines  */
 
void   einfo(int &, int &, int &);
void   dsrsm(const int &, double *, int *, int *, const int &, int *, void *,
             void *, const int &);
void   stride(const int &, const int &, int *, char *, const int &);
void   ivsset(const int &);
 
_ESVI dgktrn(const int &, double *, const int &, int *, double *, const int &,
             int *, int *, double *, const int &);
_ESVI dsktrn(const int &, double *, const int &, int *, int *, double *,
             const int &);
void ievops(const int &);
 
/* The following three definitions are added for error recovery migartion
   among different platforms     */
 
void errset(const int &, const int &, const int &, const int &,
            void *, const int &);
void errsav(const int &, void *);
void errstr(const int &, void *);
void einfo(int &, int &, int &);
 
#define errset_  errset
#define errsav_  errsav
#define errstr_  errstr
}
}
#endif
#endif
