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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     dgetrf_
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE dgetrf_(M, N, A, lda, ipiv, info)
#ifdef _CUSOLVER
       use allio, only: handle, dev_dgetrf_workspace, dev_Info
#endif
       IMPLICIT NONE
       INTEGER*4   :: M,stat, N, lda, info, lipiv,i,j
       REAL*8, DIMENSION(lda,*) :: A
       INTEGER*4, DIMENSION(*)  :: ipiv
#ifdef _CUSOLVER
#ifdef _OFFLOAD
#ifdef RISC
       CALL cusolver_dgetrf_(handle, stat, M, N, A, lda
     1 , dev_dgetrf_workspace, ipiv, dev_Info)
       CALL cudasync_
#else
       CALL cusolver_dgetrf(handle, stat, M, N, A, lda
     1 , dev_dgetrf_workspace, ipiv, dev_Info)
       CALL cudasync
#endif
!$omp target update from(dev_Info)
       info = dev_Info(1)
#endif
#else
       CALL dgetrf(M, N, A, lda, ipiv, info)
#endif
      END SUBROUTINE dgetrf_
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     zgetrf_
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE zgetrf_(M, N, A, lda, ipiv, info)
#ifdef _CUSOLVER
       use allio, only: handle, dev_zgetrf_workspace, dev_Info
#endif
       IMPLICIT NONE
       INTEGER*4   :: M,stat, N, lda, info, lipiv,i,j
       COMPLEX*16, DIMENSION(lda,*) :: A
       INTEGER*4, DIMENSION(*)  :: ipiv
#ifdef _CUSOLVER
#ifdef _OFFLOAD
#ifdef RISC
       CALL cusolver_zgetrf_(handle, stat, M, N, A, lda
     1 , dev_zgetrf_workspace, ipiv, dev_Info)
       CALL cudasync_
#else
       CALL cusolver_zgetrf(handle, stat, M, N, A, lda
     1 , dev_zgetrf_workspace, ipiv, dev_Info)
       CALL cudasync
#endif
!$omp target update from(dev_Info)
       info = dev_Info(1)
#endif
#else
       CALL zgetrf(M, N, A, lda, ipiv, info)
#endif
      END SUBROUTINE zgetrf_
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     dgetri_
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE dgetri_(N, A, lda, ipiv, work, lwork, info)
#ifdef _CUSOLVER
       use allio, only: handle, dev_dgetri_workspace, dev_Info
#endif
       INTEGER*4 :: N, lda, lwork, info
       INTEGER*4, DIMENSION(*) :: ipiv
       REAL*8, DIMENSION(*) :: work
       REAL*8, DIMENSION(lda, N) :: A
c
       INTEGER*4 :: stat, i, j
       stat = 0
c
#ifdef _CUSOLVER
       IF (lwork.LT.0) THEN
         ! Do nothing if it is lwork query
         work(1) = 1
         info = 0
         RETURN
       END IF
!$omp target teams distribute parallel do collapse(2)
       DO i = 1, N
         DO j = 1, N
           dev_dgetri_workspace(i,j) = 0.0
         END DO
       END DO
!$omp end target teams distribute parallel do
!$omp target teams distribute parallel do
       DO i = 1, N
         dev_dgetri_workspace(i,i) = 1.0
       END DO
!$omp end target teams distribute parallel do
#ifdef RISC
       CALL cusolver_dgetrs_(handle, stat, "N", N, N, A, lda, ipiv
     1 , dev_dgetri_workspace, N, dev_Info)
       CALL cudasync_
#else
       CALL cusolver_dgetrs(handle, stat, "N", N, N, A, lda, ipiv
     1 , dev_dgetri_workspace, N, dev_Info)
       CALL cudasync
#endif
!$omp target teams distribute parallel do collapse(2)
       DO i = 1, N
         DO j = 1, N
           A(i,j) = dev_dgetri_workspace(i,j)
         END DO
       END DO
!$omp end target teams distribute parallel do
!$omp target update from(dev_Info)
       info = dev_Info(1)
#else
       CALL dgetri(N, A, lda, ipiv, work, lwork, info)
#endif
      END SUBROUTINE dgetri_
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     zgetri_
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE zgetri_(N, A, lda, ipiv, work, lwork, info)
#ifdef _CUSOLVER
       use allio, only: handle, dev_zgetri_workspace, dev_Info
#endif
       INTEGER*4 :: N, lda, lwork
       INTEGER*4, DIMENSION(*) :: ipiv
       COMPLEX*16, DIMENSION(*) :: work
       COMPLEX*16, DIMENSION(lda, N) :: A
c
       INTEGER*4 :: stat, i, j
       stat = 0
c
#ifdef _CUSOLVER
       IF (lwork.LT.0) THEN
         ! Do nothing if it is lwork query
         work(1) = 1
         info = 0
         RETURN
       END IF
!$omp target teams distribute parallel do collapse(2)
       DO i = 1, N
         DO j = 1, N
           dev_zgetri_workspace(i,j) = 0.0
         END DO
       END DO
!$omp end target teams distribute parallel do
!$omp target teams distribute parallel do
       DO i = 1, N
         dev_zgetri_workspace(i,i) = 1.0
       END DO
!$omp end target teams distribute parallel do
#ifdef RISC
       CALL cusolver_zgetrs_(handle, stat, "N", N, N, A, lda, ipiv
     1 , dev_zgetri_workspace, N, dev_Info)
       CALL cudasync_
#else
       CALL cusolver_zgetrs(handle, stat, "N", N, N, A, lda, ipiv
     1 , dev_zgetri_workspace, N, dev_Info)
       CALL cudasync
#endif
!$omp target teams distribute parallel do collapse(2)
       DO i = 1, N
         DO j = 1, N
           A(i,j) = dev_zgetri_workspace(i,j)
         END DO
       END DO
!$omp end target teams distribute parallel do
!$omp target update from(dev_Info)
       info = dev_Info(1)
#else
       CALL zgetri(N, A, lda, ipiv, work, lwork, info)
#endif
      END SUBROUTINE zgetri_
