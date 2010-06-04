C  
C   
C   Copyright 1997 Preston Moore and University of Pennsylvania
C
C

C#define BLAS

c**************************************************
      subroutine lwdsyev(ijz,iul,n,a,lda,w,work,
     $     lwork,info)
c     
      integer ijz,iul,n,lda,lwork,info
      real*8 a(lda,n),work(lwork),w(n)
      character*1 cjz(0:1),cul(0:1)
      data cjz/'N','V'/,cul/'U','L'/
c
#ifdef BLAS
      call dsyev(cjz(ijz),cul(iul),n,a,lda,w,work,
     $     lwork,info)
#endif

#ifdef ESSL
      ijz = 1
      call dspev(ijz,a,w,lda,n,work,lwork)
#endif

#ifndef BLAS
#ifndef ESSL
      print *,"lwdsyev needs blas or ESSL at compile time"
      stop
#endif
#endif

      return
      end
c**************************************************
      subroutine bwdgemm(ita,itb,l,n,m,zalpha,za,lda,zb,ldb,
     $     zbeta,zc,ldc)
c
      integer ita,itb,l,m,n,lda,ldb,ldc
      real*8 zalpha,za(*),zb(*),zbeta,zc(*)
      character*1 op(0:2)
      data op/'n','t','c'/
c
#ifdef BLAS
      call dgemm(op(ita),op(itb),l,n,m,zalpha,za,lda,zb,ldb,
     $     zbeta,zc,ldc)
#else
      print *,"bwdgemm needs blas at compile time"
      stop
#endif
      return
      end
C***************************************************
      subroutine bwdgemv(ita,m,n,alpha,a,lda,x,incx,beta,y,incy)
c
      integer ita,m,n,lda,incx,incy
      real*8 alpha,a(*),x(*),beta,y(*)
      character*1 op(0:2)
      data op/'n','t','c'/
c
#ifdef BLAS
      call dgemv(op(ita),m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
      print *,"bwdgemv needs blas at compile time"
      stop
#endif
      return
      end
C***************************************************
      subroutine bwdsytrf(uplo,n,a,lda,ipiv,work,lwork,info)

      character*1 uplo
      integer n,lda,ipiv(*),lwork,info
      real*8 a(*),work(*)

#ifdef BLAS
      call dsytrf(uplo,n,a,lda,ipiv,work,lwork,info) 
#endif
      return
      end

C***************************************************
      subroutine bwdsytri(uplo,n,a,lda,ipiv,work,info)

      character*1 uplo
      integer n,lda,ipiv(*),info
      real*8 a(*),work(*)

#ifdef BLAS
      call dsytri(uplo,n,a,lda,ipiv,work,info) 
#else
      print *,"bwdsytri needs blas at compile time"
      stop
#endif
      return
      end

C***************************************************
      subroutine bwdsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
      character*1 uplo
      integer info,lda,ldb,n,nrhs,ipiv(*)
      real*8 a(*),b(*)

#ifdef BLAS
      call dsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
#else
      print *,"bwdsytrs needs blas at compile time"
      stop
#endif
      return
      end
