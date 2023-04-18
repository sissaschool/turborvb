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
* Copyright (c) 1992-2010 The University of Tennessee and
*                         The University of of Tennessee
*                         Research Foundation.
*                         All rights reserved.
* Copyright (c) 2000-2010 The University of California Berkeley.
*                         All rights reserved.
* Copyright (c) 2006-2010 The University of Colorado Denver.
*                         All rights reserved.
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     zgemvt
c
c     Subroutine performing lapack's zgemv with operation 'T'.
c     It has the same arguments except the first one is missing
c     and it is assumed to be 'T'.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE zgemvt_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
c
      IMPLICIT NONE
c
      COMPLEX*16 :: alpha, beta
      INTEGER*4 :: incx, incy, lda, M, N
      COMPLEX*16 :: A(lda*(N-1)+M), X(incx*(N-1)+1), Y(incy*(M-1)+1)
c
      REAL*8 :: one, zero
      PARAMETER (one=1.0D+0,zero=0.0D+0)
c
      INTEGER*4 ::I, J
      COMPLEX*16 :: temp
c
#ifdef _OFFLOAD
      IF (beta.NE.one) THEN
          IF (beta.EQ.zero) THEN
!$omp target teams distribute parallel do
              DO  I = 1,N
                  Y(incy*(I-1)+1) = zero
              END DO
          ELSE
!$omp target teams distribute parallel do
              DO  I = 1,N
                  Y(incy*(I-1)+1) = beta*Y(incy*(I-1)+1)
              END DO
          END IF
      END IF
c
      IF (alpha.NE.zero) THEN
!$omp target teams distribute private(temp)
          DO  J = 1,N
              temp = zero
!$omp parallel do reduction(+:temp)
              DO  I = 1,M
                  temp = temp + A(lda*(J-1)+I)*X(incx*(I-1)+1)
              END DO
              Y(incy*(J-1)+1) = Y(incy*(J-1)+1) + alpha*temp
          END DO
      END IF
#endif

      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     zgemvc
c
c     Subroutine performing lapack's zgemv with operation 'C'.
c     It has the same arguments except the first one is missing
c     and it is assumed to be 'C'.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE zgemvc_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
c
      IMPLICIT NONE
c
      COMPLEX*16 :: alpha, beta
      INTEGER*4 :: incx, incy, lda, M, N
      COMPLEX*16 :: A(lda*(N-1)+M), X(incx*(N-1)+1), Y(incy*(M-1)+1)
c
      REAL*8 :: one, zero
      PARAMETER (one=1.0D+0,zero=0.0D+0)
c
      INTEGER*4 ::I, J
      COMPLEX*16 :: temp
c
#ifdef _OFFLOAD
      IF (beta.NE.one) THEN
          IF (beta.EQ.zero) THEN
!$omp target teams distribute parallel do
              DO  I = 1,N
                  Y(incy*(I-1)+1) = zero
              END DO
          ELSE
!$omp target teams distribute parallel do
              DO  I = 1,N
                  Y(incy*(I-1)+1) = beta*Y(incy*(I-1)+1)
              END DO
          END IF
      END IF
c
      IF (alpha.NE.zero) THEN
!$omp target teams distribute private(temp)
          DO  J = 1,N
              temp = zero
!$omp parallel do reduction(+:temp)
              DO  I = 1,M
                  temp = temp + DCONJG(A(lda*(J-1)+I)) * X(incx*(I-1)+1)
              END DO
              Y(incy*(J-1)+1) = Y(incy*(J-1)+1) + alpha * temp
          END DO
      END IF
#endif

      RETURN

      END
