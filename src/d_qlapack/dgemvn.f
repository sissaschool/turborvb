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
c     dgemvn
c
c     Subroutine performing lapack's dgemv with operation 'N'.
c     It has the same arguments except the first one is missing
c     and it is assumed to be 'N'.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE dgemvn_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
c
      IMPLICIT NONE
      REAL*8 :: alpha, beta, A(lda*(n-1)+m)
     +        , X((n-1)*incx+1), Y((m-1)*incy+1)
      INTEGER*4 incx, incy, lda, M, N
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
                      Y(incy*(I-1)+1) = zero
                  END DO
              ELSE
!$omp target teams distribute parallel do
                  DO I = 1, M
                      Y(incy*(I-1)+1) = BETA*Y(incy*(I-1)+1)
                  END DO
              END IF
      END IF
c
      IF (alpha.NE.zero) THEN
!$omp target teams distribute private(temp)
          DO I = 1,M
              temp=0.d0
!$omp parallel do reduction(+:temp)
              DO J = 1,N
                  temp = temp+A(lda*(J-1)+I)*X(incx*(J-1)+1)
              END DO
              Y(incy*(I-1)+1) = Y(incy*(I-1)+1)+alpha*temp
          END DO
      END IF
#endif
      RETURN
      END
