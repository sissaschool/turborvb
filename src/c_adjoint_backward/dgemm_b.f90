!TL off
subroutine dgemm_b(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, AB, LDAB, B, LDB, BB&
        &, LDBB, BETA, CB, LDCB)
    implicit none
    CHARACTER*1        TRANSA, TRANSB
    CHARACTER*1        TRANSAB, TRANSBB
    INTEGER            I,J,M, N, K, LDA, LDAB, LDB, LDBB, LDCB
    DOUBLE PRECISION   ALPHA, BETA
    !     .. Array Arguments ..
    DOUBLE PRECISION   A(LDA, *), B(LDB, *)
    DOUBLE PRECISION   AB(LDAB, *), BB(LDBB, *), CB(LDCB, *)

    !         true dimensions:
    !                        C(M,N)

    !         CASE      'N'             'T'
    !                 A(M,K)           A(K,M)
    !                    'N'            'T'
    !                 B(K,N)           B(N,K)
    !         C is not used and is not passed.
    !    CB is instead used as an input and is modified in output.
    !    A and B are input not modified in output.
    !    Output  AB and BB and CB

    if(transa.eq.'n'.or.transa.eq.'N') then
        transab = 'N'
    else
        transab = 'T'
    endif

    if(transb.eq.'n'.or.transb.eq.'N') then
        transbb = 'N'
    else
        transbb = 'T'
    endif

    if(transab.eq.'N'.and.transbb.eq.'N') then

        call dgemm_('N', 'T', m, k, n, alpha, cb, ldcb, b, ldb, 1.d0, ab, ldab)
        call dgemm_('T', 'N', k, n, m, alpha, a, lda, cb, ldcb, 1.d0, bb, ldbb)

    elseif(transab.ne.'N'.and.transbb.eq.'N') then

        call dgemm_('N', 'T', k, m, n, alpha, b, ldb, cb, ldcb, 1.d0, ab, ldab)
        call dgemm_('N', 'N', k, n, m, alpha, a, lda, cb, ldcb, 1.d0, bb, ldbb)

    elseif(transab.eq.'N'.and.transbb.ne.'N') then

        call dgemm_('N', 'N', m, k, n, alpha, cb, ldcb, b, ldb, 1.d0, ab, ldab)
        call dgemm_('T', 'N', n, k, m, alpha, cb, ldcb, a, lda, 1.d0, bb, ldbb)

    elseif(transab.ne.'N'.and.transbb.ne.'N') then

        call dgemm_('T', 'T', k, m, n, alpha, b, ldb, cb, ldcb, 1.d0, ab, ldab)
        call dgemm_('T', 'T', n, k, m, alpha, cb, ldcb, a, lda, 1.d0, bb, ldbb)

    endif
    call dscalmatrix(m,n,beta,cb,ldcb)
    return
end
subroutine dscalmatrix(m,n,alpha,a,lda)
implicit none
integer i,j,n,m,lda
real*8 alpha,a(lda*(n-1)+m)
if(alpha.ne.0.d0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
do j=1,n
 do i=1,m
 a(lda*(j-1)+i)=a(lda*(j-1)+i)*alpha
 enddo
enddo
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do 
#endif
else
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
do j=1,n
 do i=1,m
 a(lda*(j-1)+i)=0.d0
 enddo
enddo
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do 
#endif
endif
return
end

