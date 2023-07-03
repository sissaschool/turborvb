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
c     dgemvt
c
c     Subroutine performing lapack's dgemv with operation 'T'.
c     It has the same arguments except the first one is missing
c     and it is assumed to be 'T'.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE dgemvt_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
c
      IMPLICIT NONE
      REAL*8 :: alpha, beta, a(lda*(n-1)+m)
     +        , x((n-1)*incx+1), y((m-1)*incy+1)
      INTEGER*4 :: incx, incy, lda, m, n
c
      DOUBLE PRECISION one,zero
      PARAMETER (one=1.0D+0,zero=0.0D+0)
c
      INTEGER i, j
      REAL*8 :: temp
c
#ifdef _OFFLOAD
      IF (beta.NE.one) THEN
              IF (beta.eq.zero) THEN
!$omp target teams distribute parallel do
                  DO I = 1, M
                      Y(INCY*(I-1)+1) = ZERO
                  END DO
              ELSE
!$omp target teams distribute parallel do
                  DO I = 1,M
                      Y(INCY*(I-1)+1) = BETA*Y(INCY*(I-1)+1)
                  END DO
              END IF
      END IF
c
      IF (alpha.NE.zero) THEN
!$omp target teams distribute private(temp)
              DO J = 1,N
                  temp=0.d0
!$omp parallel do reduction(+:temp)
                  DO I = 1,M
                    temp = temp + A((J-1)*lda+I)*X(incx*(I-1)+1)
                  END DO
                  Y(incy*(J-1)+1) = Y(incy*(J-1)+1) + alpha*temp
            END DO
      END IF
#endif
      RETURN
      END
