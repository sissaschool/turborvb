!TL off
subroutine dgemv_b(TRANS, M, N, ALPHA, A, LDA, AB, LDAB, X, LDX, XB &
        &, LDXB, BETA, YB, LDYB)
    implicit none
    character*1 trans
    integer M, N, LDA, LDAB, LDX, LDXB, LDYB
    real*8 alpha, beta
    real*8 A(LDA, *), AB(LDAB, *), X(*), XB(*), YB(*)
    if(trans.eq.'n'.or.trans.eq.'N') then
        call dgemv('T', M, N, alpha, A, LDA, YB, LDYB, 1.d0, XB, LDXB)
        call dger(M, N, alpha, YB, LDYB, X, LDX, AB, LDAB)
        call dscal(M, beta, yb, ldyb)
    else
        call dgemv('N', M, N, alpha, A, LDA, YB, LDYB, 1.d0, XB, LDXB)
        call dger(M, N, alpha, X, LDX, YB, LDYB, AB, LDAB)
        call dscal(N, beta, yb, ldyb)
    endif
    return
end

subroutine zgemv_b(TRANS, M, N, ALPHA, A, LDA, AB, LDAB, X, LDX, XB &
        &, LDXB, BETA, YB, LDYB)
    implicit none
    character*1 trans
    integer M, N, LDA, LDAB, LDX, LDXB, LDYB
    complex*16, parameter :: one=(1.d0,0.d0)
    complex*16 alpha, beta,alphac,betac
    complex*16 A(LDA, *), AB(LDAB, *), X(*), XB(*), YB(*)
    alphac=conjg(alpha)
    betac=conjg(beta)
    if(trans.eq.'n'.or.trans.eq.'N') then
        call zgemv('C', M, N, alphac, A, LDA, YB, LDYB, one, XB, LDXB)
        call zgerc(M, N, alphac, YB, LDYB, X, LDX, AB, LDAB)
        call zscal(M,betac, yb, ldyb)
    else
!       to avoid conjugate of matrix
        call conjmat(1,M,XB,LDXB)
        call conjmat(1,N,YB,LDYB)
        call zgemv('N', M, N, alpha, A, LDA, YB, LDYB, one, XB, LDXB)
        call conjmat(1,M,XB,LDXB)
        call conjmat(1,N,YB,LDYB)

        call conjmat(1,M,X,LDX)
        call zgeru(M,N,alphac, X, LDX, YB, LDYB, AB, LDAB)
        call conjmat(1,M,X,LDX)
        call zscal(N,betac,yb,ldyb)
    endif
    return
end

 
