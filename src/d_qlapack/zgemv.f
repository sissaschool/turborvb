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
      SUBROUTINE ZGEMV_(trans,M,N,alpha,A,lda,X,incx,beta
     +,Y,incy,YES_ONTARGET)

      IMPLICIT NONE

      COMPLEX*16 :: alpha, beta, A(lda, *), X(*), Y(*)
      INTEGER*4 :: incx, incy, lda, M, N
      CHARACTER :: trans
      LOGICAL :: yes_ontarget

      COMPLEX*16 :: one, zero
      PARAMETER (one=(1.0D+0,0.0D+0), zero=(0.0D+0,0.0D+0))
      INTEGER*4 :: I, J, info, lenx, leny

      INTRINSIC MAX

#ifndef _OFFLOAD
      CALL ZGEMV(trans,M,N,alpha,A,lda,X,incx,beta,Y,incy)
#else
#ifdef _CUBLAS
      IF(.NOT.yes_ontarget) THEN
          CALL ZGEMV(trans,M,N,alpha,A,lda,X,incx,beta,Y,incy)
      ELSE
#ifdef RISC
           CALL cublas_zgemv_offload_(trans,M,N,alpha,A,lda,x
     +,incx,beta,y,incy)
           CALL cudasync_
#else
           CALL cublas_zgemv_offload(trans,M,N,alpha,A,lda,x
     +,incx,beta,y,incy)
           CALL cudasync
#endif
      END IF
#else
      IF (.NOT.yes_ontarget) THEN
          CALL ZGEMV(trans,M,N,alpha,A,lda,X,incx,beta,Y,incy)
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
     +    ((alpha.EQ.zero) .AND. (beta.EQ.one))) RETURN

      IF (trans.EQ.'N'.OR.trans.EQ.'n') THEN
          CALL zgemvn_offload(M,N,alpha,A,lda,X,incx,beta,Y,incy)
      END IF
      IF (trans.EQ.'T'.OR.trans.EQ.'t') THEN
          CALL zgemvt_offload(M,N,alpha,A,lda,X,incx,beta,Y,incy)
      END IF
      IF (trans.EQ.'C'.OR.trans.EQ.'c') THEN
          CALL zgemvc_offload(M,N,alpha,A,lda,X,incx,beta,Y,incy)
      END IF
#endif
#endif

      RETURN

      END
