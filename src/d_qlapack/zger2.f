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
      SUBROUTINE ZGER2(M,N,ALPHA,X,LDX,Y,LDY,A,LDA)
      USE constants,  ONLY: yes_ontarget

      IMPLICIT NONE
c
      COMPLEX*16 ALPHA
      COMPLEX*16 :: A(LDA*(N-1)+M),X(LDX,2),Y(LDY,2)
      INTEGER*4 :: LDX, LDY, LDA, M, N
c
      COMPLEX*16 :: ONE, ZERO
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      PARAMETER (ZERO= (0.0D+0,0.0D+0))

#ifdef _OFFLOAD
      INTEGER*4 :: I, J, INFO
#endif
c
      IF(.NOT.yes_ontarget) THEN
          CALL ZGEMM ('N','T',M, N, 2, ALPHA
     +               ,X, LDX, Y, LDY, ONE, A, LDA)
          RETURN
      END IF
c
#ifdef _OFFLOAD
#ifdef _CUBLAS
#ifdef RISC
      CALL cublas_zgemm_offload_('N','T',M,N,2,alpha,X,LDX,Y
     +                          ,LDY,ONE,A,LDA)
      CALL cudasync_
#else
      CALL cublas_zgemm_offload('N','T',M,N,2,alpha,X,LDX,Y
     +                         ,LDY,ONE,A,LDA)
      CALL cudasync
#endif
#else
c
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
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
     +                         ALPHA*(X(I,1)*Y(J,1)+X(I,2)*Y(J,2))
          END DO
      END DO
#endif
#endif
      RETURN
      END
