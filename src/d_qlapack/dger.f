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
      SUBROUTINE DGER_(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE allio, ONLY: yes_ontarget
c
      IMPLICIT NONE
c
      REAL*8 :: ALPHA
      REAL*8 :: A(LDA*(N-1)+M), X(INCX*(M-1)+1), Y(INCY*(N-1)+1)
      INTEGER*4 :: INCX, INCY, LDA, M, N
      INTEGER I, INFO, J
c
      IF(.NOT.yes_ontarget) THEN
          CALL DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
          RETURN
      END IF
c
#ifdef _OFFLOAD
#ifdef _CUBLAS
#ifdef RISC
      CALL cublas_dger_offload_(M,N,alpha,X,incx,Y
     +                          ,incy,A,LDA)
      CALL cudasync_
#else
      CALL cublas_dger_offload(M,N,alpha,X,incx,Y
     +                        ,incy,A,LDA)
      CALL cudasync
#endif
#else
c
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
c
      IF ((INFO.NE.0)
     +    .OR.(M.EQ.0)
     +    .OR.(N.EQ.0)
     +    .OR.(ALPHA.EQ.0))
     +  RETURN
c
!$omp target teams distribute parallel do collapse(2)
      DO J = 1,N
          DO I = 1,M
              A(LDA*(J-1)+I) = A(LDA*(J-1)+I) +
     1                         ALPHA*X(INCX*(I-1)+1)*Y(INCY*(J-1)+1)
          END DO
      END DO
#endif
#endif
      RETURN
      END
