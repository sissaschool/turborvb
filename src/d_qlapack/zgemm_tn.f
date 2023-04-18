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
      SUBROUTINE ZGEMM_TN(M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
c
      USE constants, ONLY: yes_ontarget
c
      IMPLICIT NONE
c
      REAL*8 :: A(LDA,*), B(LDB,*), C(LDC,*)
      REAL*8 :: ALPHA, BETA
      INTEGER*4 :: M, N, K, LDA, LDB, LDC
c
      IF(n.EQ.0.OR.m.EQ.0.OR.k.EQ.0) RETURN
c
      IF(.NOT.yes_ontarget) then
          CALL ZGEMM('T','N',M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
          RETURN
      END IF
c
#ifndef _OFFLOAD
      CALL ZGEMM('T','N',M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#else
#ifdef _CUBLAS
#ifdef RISC
      CALL cublas_zgemm_offload_('T','N',M,N,K,alpha,A,LDA
     +                          ,b,LDB,beta,c,LDC)
      CALL cudasync_
#else
      CALL cublas_zgemm_offload('T','N',M,N,K,alpha,A,LDA
     +                         ,b,LDB,beta,c,LDC)
      CALL cudasync
#endif
#else
!$omp target update from(A(1:lda,1:n))
!$omp target update from(B(1:lda,1:n))
      CALL ZGEMM('T','N',M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!$omp target update to(A(1:lda,1:n))
!$omp target update to(B(1:lda,1:n))
#endif
#endif
c
      RETURN
      END

      SUBROUTINE ZGEMM_(TRANA,TRANB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
c
      IMPLICIT NONE
c
      COMPLEX*16 :: A(LDA,*), B(LDB,*), C(LDC,*)
      COMPLEX*16 :: ALPHA, BETA
      INTEGER*4 :: M, N, K, LDA, LDB, LDC
      CHARACTER(len=1) :: TRANA, TRANB
c
      IF(n.EQ.0.OR.m.EQ.0.OR.k.EQ.0) RETURN
c
#ifndef _OFFLOAD
      CALL ZGEMM(TRANA,TRANB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#else
#ifdef _CUBLAS
#ifdef RISC
      CALL cublas_zgemm_offload_(TRANA,TRANB,M,N,K,alpha,A,LDA
     +                          ,b,LDB,beta,c,LDC)
      CALL cudasync_
#else
      CALL cublas_zgemm_offload(TRANA,TRANB,M,N,K,alpha,A,LDA
     +                         ,b,LDB,beta,c,LDC)
      CALL cudasync
#endif
#endif
#endif
      RETURN

      END

      SUBROUTINE ZTRSM_(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
c
      IMPLICIT NONE
c
      COMPLEX*16 :: ALPHA
      COMPLEX*16 :: A(LDA,*),B(LDB,*)
      INTEGER*4 :: lda,ldb,M,N
      CHARACTER(len=1) :: SIDE, UPLO, TRANSA, DIAG
c
#ifndef _OFFLOAD
      CALL  ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
#else
#ifdef _CUBLAS
#ifdef RISC
      CALL cublas_ztrsm_offload_(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA
     +                          ,A,LDA,B,LDB)
      CALL cudasync_
#else
      CALL cublas_ztrsm_offload(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA
     +                         ,A,LDA,B,LDB)
      CALL cudasync
#endif
#else
!$omp target update from(A(1:lda,1:n))
!$omp target update from(B(1:lda,1:n))
      CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!$omp target update to(A(1:lda,1:n))
!$omp target update to(B(1:lda,1:n))
#endif
#endif
c
      RETURN
c
      END

