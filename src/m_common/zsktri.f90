! Copyright (C) 2022 TurboRVB group
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

!> This subroutine manipulates a skewsymmetric matrix A.
subroutine zsktri(uplo, n, a, lda, ainv, ldinv, ipiv, work, info)
    implicit none

    ! argument parameters
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda, ldinv, ipiv(*)
    complex*16, intent(in out) :: a(lda, *), work(*)
    complex*16, intent(out) :: ainv(ldinv, *)
    integer, intent(out) :: info

    ! local variables
    integer i, j
    real*8 rcond
    logical yeslap

    !     First factorization assumed
    !     lwork=3*n
    !     CALL ZSKTRF( UPLO, 'N', N, A, LDA, IPIV, WORK, LWORK, INFO)
    !     Identity matrix
    !     ipiv dimension required 3n
    !     work dimension required n^2+12*n-2
    info = 0
    yeslap = .false. ! if false the homemade algorithm is done.
    do i = 1, n
        do j = 1, i - 1
            ainv(j, i) = dcmplx(0.d0, 0.d0)
        end do
        ainv(i, i) = dcmplx(1.d0, 0.d0)
        do j = i + 1, n
            ainv(j, i) = dcmplx(0.d0, 0.d0)
        end do
    end do
    ! fisrt permutation of ainv
    do i = n, 1, -1
        work(1:n) = ainv(i, 1:n)
        ainv(i, 1:n) = ainv(ipiv(i), 1:n)
        ainv(ipiv(i), 1:n) = work(1:n)
    end do
    ! second
    if (UPLO .eq. 'u' .or. UPLO .eq. 'U') then
        do i = 1, N - 1
            work(i) = a(i, i + 1)
            work(n + i - 1) = -a(i, i + 1)
        end do
        do j = 2, n - 1
            a(1:j - 1, j) = a(1:j - 1, j + 1)
        end do
        a(1:n - 1, n) = dcmplx(0.d0, 0.d0)
    else
        do i = 1, N - 1
            work(i) = -a(i + 1, i)
            work(n + i - 1) = a(i + 1, i)
        end do
        do j = n - 1, 2, -1
            a(j + 1:n, j) = a(j + 1:n, j - 1)
        end do
        a(2:n, 1) = dcmplx(0.d0, 0.d0)
    end if
    work(2*N - 1:3*N - 2) = dcmplx(0.d0, 0.d0) ! diagonal elements of skew matrix , obviously set to zero.
    call ZTRTRS(UPLO, 'N', 'U', N, N, A, LDA, ainv, LDINV, INFO) ! ainv = A^-1
    ! USE STANDARD LAPACK DGTSVX or the simples DGTSV?
    if (yeslap) then
        call ZGTSVX('N', 'N', N, N, work(N), work(2*N - 1), work, work(3*N)&
                &, work(4*N), work(5*N), work(6*N), IPIV(N + 1), AINV, LDINV, work(7*N), N&
                &, RCOND, work(7*N + N*N - 1), work(8*N + N*N - 1)&
                &, work(9*N + N*N - 1), work(11*N + N*N - 1), INFO)
        do i = 1, N
            ainv(1:N, i) = work(7*N + (i - 1)*N:7*N + i*N - 1)
        end do
    else
        call ZSKTRS(UPLO, N, N, WORK, AINV, LDINV, WORK(N), INFO)
    end if
    call ZTRTRS(UPLO, 'T', 'U', N, N, A, LDA, ainv, LDINV, INFO)
    !     last permutation of ainv
    do i = 1, n
        work(1:n) = ainv(i, 1:n)
        ainv(i, 1:n) = ainv(ipiv(i), 1:n)
        ainv(ipiv(i), 1:n) = work(1:n)
    end do
    return
end

