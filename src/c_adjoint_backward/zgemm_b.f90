!TL off
subroutine zgemm_b(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, AB, LDAB, B, LDB, BB&
        &, LDBB, BETA, CB, LDCB)
    use constants, only : zzero, zone
    implicit none
    CHARACTER*1        TRANSA, TRANSB
    CHARACTER*1        TRANSAB, TRANSBB
    INTEGER            I,J,M, N, K, LDA, LDAB, LDB, LDBB, LDCB
    COMPLEX*16   ALPHA, BETA, alphac,betac
    !     .. Array Arguments ..
    COMPLEX*16  A(LDA, *), B(LDB, *)
    COMPLEX*16   AB(LDAB, *), BB(LDBB, *), CB(LDCB, *)

    !         true dimensions:

    !                        C(M,N)

    !         CASE      'N'             'T'
    !                 A(M,K)           A(K,M)
    !                    'N'            'T'
    !                 B(K,N)           B(N,K)
    !         C is not used and is not passed.
    !    CB is instead used as an input and is modified in output
    !    A and B are input not modified output  with some trick (see conja conjb ..)
    !    Output  AB and BB and CB
    alphac = conjg(alpha)
    betac = conjg(beta)

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
        call zgemm_('N', 'C', m, k, n, alphac, cb, ldcb, b, ldb, zone, ab, ldab)
        call zgemm_('C', 'N', k, n, m, alphac, a, lda, cb, ldcb, zone, bb, ldbb)
    elseif(transab.eq.'T'.and.transbb.eq.'N') then
        call conjar
        call conjbr
        call zgemm_('N', 'T', k, m, n, alphac, b, ldb, cb, ldcb, zone, ab, ldab)
        call zgemm_('N', 'N', k, n, m, alphac, a, lda, cb, ldcb, zone, bb, ldbb)
        call conjar
        call conjbr
    elseif(transab.eq.'N'.and.transbb.eq.'T') then
        call conja
        call conjb
        call zgemm_('N', 'N', m, k, n, alphac, cb, ldcb, b, ldb, zone, ab, ldab)
        call zgemm_('T', 'N', n, k, m, alphac, cb, ldcb, a, lda, zone, bb, ldbb)
        call conja
        call conjb
    elseif(transab.eq.'T'.and.transbb.eq.'T') then
        call zgemm_('C', 'T', k, m, n, alphac, b, ldb, cb, ldcb, zone, ab, ldab)
        call zgemm_('T', 'C', n, k, m, alphac, cb, ldcb, a, lda, zone, bb, ldbb)
    endif

    call zscalmatrix(m,n,betac,cb,ldcb)
    return

contains

    subroutine conja
        implicit none
        integer i, j
!       do j = 1, k
!           do i = 1, m
!               a(i, j) = dconjg(a(i, j))
!           enddo
!       enddo
    call conjmat_(m,k,a,lda)
    end subroutine conja

    subroutine conjar
        implicit none
!       integer i, j
!       do j = 1, m
!           do i = 1, k
!               a(i, j) = dconjg(a(i, j))
!           enddo
!       enddo
    call conjmat_(k,m,a,lda)
    end subroutine conjar

    subroutine conjb
        implicit none
!       integer i, j
!       do j = 1, k
!           do i = 1, n
!               b(i, j) = dconjg(b(i, j))
!           enddo
!       enddo
    call conjmat_(n,k,b,ldb)
    end subroutine conjb

    subroutine conjbr
        implicit none
!       integer i, j
!       do j = 1, n
!           do i = 1, k
!               b(i, j) = dconjg(b(i, j))
!           enddo
!       enddo
     call conjmat_(k,n,b,ldb)
    end subroutine conjbr
end subroutine zgemm_b
subroutine zscalmatrix(m,n,alpha,a,lda)
implicit none
integer i,j,n,m,lda
complex*16  alpha,a(lda*(n-1)+m)
complex*16, parameter:: zero=dcmplx(0.d0,0.d0)
if(alpha.ne.zero) then
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
 a(lda*(j-1)+i)=zero
 enddo
enddo
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do 
#endif
endif
return
end
subroutine conjmat_(m,n,a,lda)
implicit none
integer i,j,n,m,lda
complex*16  a(lda*(n-1)+m)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
do j=1,n
 do i=1,m
 a(lda*(j-1)+i)=dconjg(a(lda*(j-1)+i))
 enddo
enddo
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do 
#endif
return
end


