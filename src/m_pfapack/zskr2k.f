cTL off
*
* The 3-Clause BSD License
*
* Copyright (c) 2022 TurboRVB group based on the following works
*
* Pfapack: https://michaelwimmer.org/downloads.html
* Copyright (c) 2010 Michael Wimmer, Universiteit Leiden.
*
* Any use of this software is allowed as far as the author is concerned.
* The author only asks users that use this software and obtain
* scientific results, to honor using this software by citing the paper
* accompanying this library: http://arxiv.org/abs/1102.3440 and
* ACM Trans. Math. Software 38 30 (2012).
*
* The FORTRAN implementation is derived in parts on the reference
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
*
      SUBROUTINE ZSKR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA
      DOUBLE COMPLEX BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
*  MODIFIED BY SANDRO SORELLA ON 20/6/2018. NOWADAYS DOES NOT MAKE ANY SENSE TO WRITE DOWN 
*  HOMEMADE BLAS-3 ROUTINES TO SAVE A FACTOR 2 FLOPS AND/OR MEMORY. I JUST USED STANDARD BLAS3
*  BY DOING TWICE THE COMPUTATIONAL EFFORT BUT OBTAINING AN EFFICIENT CODE ALSO OPENMP THREADED

*     ..
*
*  Purpose
*  =======
*
*  ZSKR2K  performs one of the skew-symmetric rank 2k operations
*
*     C := alpha*A*B^T - alpha*B*A^T + beta*C,
*
*  or
*
*     C := alpha*A^T*B - alpha*B^T*A + beta*C,
*
*  where  alpha and beta  are scalars,  C is an  n by n
*  skew-symmetric matrix and  A and B  are  n by k matrices in the first case
*  and  k by n  matrices in the second case.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'    C := alpha*A*B^T   -
*                                         alpha*B*A^T   +
*                                         beta*C.
*
*              TRANS = 'T' or 't'    C := alpha*A^T*B   -
*                                         alpha*B^T*A   +
*                                         beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns  of the  matrices  A and B,  and on  entry  with
*           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
*           matrices  A and B.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE COMPLEX  array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE COMPLEX   array of DIMENSION ( LDB, kb ), where kb is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  k by n  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE COMPLEX          .
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - DOUBLE COMPLEX  array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  skew-symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  skew-symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*           Note that the diagonal elements need
*           not be set,  they are assumed to be zero,  and on exit they
*           are set to zero.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 10/22/2010
*     Michael Wimmer, Universiteit Leiden
*     Based on ZHER2K from BLAS (www.netlib.org)
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))

*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZSKR2K',INFO)
          RETURN
      END IF
*
      IF (LSAME(TRANS,'N')) THEN
*     C := alpha*A*B^T - alpha*B*A^T + beta*C,
!     allocate(ab(N,4*k))
!     ab(1:N,1:k)=A(1:N,1:k)
!     ab(1:N,k+1:2*k)=-B(1:N,1:k)
!     ab(1:N,2*k+1:3*k)=B(1:N,1:k)
!     ab(1:N,3*k+1:4*k)=A(1:N,1:k)
!     call zgemm('N','T',N,N,2*k,alpha,ab,N,ab(1,2*k+1),N,beta,C,LDC)
      call zgemm('N','T',N,N,k,alpha,a,lda,b,ldb,beta,C,LDC)
      call zgemm('N','T',N,N,k,-alpha,b,ldb,a,lda,one,C,LDC)
      else
*     C := alpha*A^T*B - alpha*B^T*A + beta*C,
!     allocate(ab(4*k,N))
!     ab(1:k,1:N)=A(1:k,1:N)
!     ab(k+1:2*k,1:N)=-B(1:k,1:N)
!     ab(2*k+1:3*k,1:N)=B(1:k,1:N)
!     ab(3*k+1:4*k,1:N)=A(1:k,1:N)
!     call zgemm('T','N',N,N,2*k,alpha,ab,4*k,ab(2*k+1,1),4*k
!    *,beta,C,LDC)
      call zgemm('T','N',N,N,k,alpha,a,lda,b,ldb,beta,C,LDC)
      call zgemm('T','N',N,N,k,-alpha,b,ldb,a,lda,one,C,LDC)
      endif

!     deallocate(ab)

*
      RETURN
*
*     End of ZSKR2K.
*
      END
