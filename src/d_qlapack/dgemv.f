cTL off
*
* The 3-Clause BSD License
*
* Copyright (c) 2022 TurboRVB group based on the following works
*
* This fortran implementation is derived in parts on the reference
* implementation of LAPACK and BLAS on www.netlib.org.
* LAPACK is subject to the following license:
*
* Copyright (c) 1992-2010 The University of Tennessee and The University
*                         of Tennessee Research Foundation.
*                         All rights reserved.
* Copyright (c) 2000-2010 The University of California Berkeley.
*                         All rights reserved.
* Copyright (c) 2006-2010 The University of Colorado Denver.
*                         All rights reserved.
*
      SUBROUTINE DGEMV_(trans,M,N,alpha,A,lda,X,incx,beta
     +,Y,incy,yes_ontarget)

      IMPLICIT NONE

      REAL*8 :: alpha, beta, A(lda, *), X(*), Y(*)
      INTEGER*4 :: incx, incy, lda, M, N
      CHARACTER :: trans
      LOGICAL :: yes_ontarget

      REAL*8 :: one, zero
      PARAMETER (one=1.0D+0,zero=0.0D+0)
      INTEGER*4 :: I, J, info, lenx, leny

      INTRINSIC MAX

#ifndef _OFFLOAD
      CALL DGEMV(trans,M,N,alpha,A,lda,X,incx,beta,Y,incy)
#else
#ifdef _CUBLAS
      IF(.NOT.yes_ontarget) THEN
      CALL DGEMV(trans,M,N,alpha,A,lda,X,incx,beta,Y,incy)
      ELSE
#ifdef RISC
           CALL cublas_dgemv_offload_(trans,M,N,alpha,A,lda,X
     +                               ,incx,beta,Y,incy)
           CALL cudasync_
#else
           CALL cublas_dgemv_offload(trans,M,N,alpha,A,lda,X
     +                              ,incx,beta,Y,incy)
           CALL cudasync
#endif
      END IF
#else
      IF(.NOT.yes_ontarget) THEN
          CALL DGEMV(trans,M,N,alpha,A,lda,X,incx,beta,Y,incy)
          RETURN
      END IF

      INFO = 0
      IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (lda.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (incx.EQ.0) THEN
          INFO = 8
      ELSE IF (incy.EQ.0) THEN
          INFO = 11
      END IF

      IF (INFO.NE.0 .OR.
     +    M.EQ.0 .OR.
     +    N.EQ.0 .OR.
     +    ((alpha.EQ.ZERO) .AND. (beta.EQ.ONE))) RETURN

      IF (trans.EQ.'N'.OR.trans.EQ.'n') THEN
          CALL dgemvn_offload(M,N,alpha,A,lda,X,incx,beta,Y,incy)
      ELSE
          CALL dgemvt_offload(M,N,alpha,A,lda,X,incx,beta,Y,incy)
      END IF
#endif
#endif

      RETURN

      END

      SUBROUTINE DGEMV__(trans,M,N,alpha,A,lda,X,incx,beta
     +,Y,incy,yes_ontarget)

      IMPLICIT NONE

      REAL*8 :: alpha, beta
      REAL*8 :: A(lda,*),X(*),Y(*)
      INTEGER*4 :: incx, incy, lda, M, N
      CHARACTER :: trans
      LOGICAL :: yes_ontarget

      INTEGER*4 :: XLEN, YLEN

      IF(.NOT.yes_ontarget) THEN
          CALL DGEMV(trans,M,N,alpha,A,lda,X,incx,beta
     +    ,Y,incy)
      ELSE
*         THE MATRIX IS ASSUMED ON THE GPU VECTORS IN CPU
          XLEN = -1
          YLEN = -1
          IF(trans.EQ.'N'.OR.trans.EQ.'n') THEN
              XLEN = incx*N
              YLEN = incy*M
          END IF
          IF(trans.EQ.'T'.OR.trans.EQ.'t') THEN
              XLEN = incx*M
              YLEN = incy*N
          END IF
!TL off
#ifdef _OFFLOAD
!$omp target update to (x(1:xlen))
!$omp target update to (y(1:ylen)) if(beta.ne.0)
#endif
          CALL DGEMV_(trans,M,N,alpha,A,lda,X,incx,beta
     +                ,Y,incy,.true.)
#ifdef _OFFLOAD
!$omp target update from  (y(1:ylen))
#endif
      END IF
      RETURN
      END
