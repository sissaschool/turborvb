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

subroutine dsktri(uplo, n, a, lda, ainv, ldinv, ipiv, work, info)
    implicit none
    character uplo
    integer n, i, j, lda, ldinv, info, lwork, ipiv(*)
    real*8 a(lda, *), ainv(ldinv, *), work(*)
    real*8 rcond
    logical yeslap
    !     First factorization assumed
    !     CALL DSKTRF( UPLO, 'N', N, A, LDA, IPIV, WORK, LWORK, INFO)
    !     Identity matrix
    !     ipiv dimension required 3n
    !     work dimension required n^2+12*n-2
    yeslap = .false. ! if .false. the homemade algorithm is done.
    do i = 1, n
        do j = 1, i - 1
            ainv(j, i) = 0.d0
        end do
        ainv(i, i) = 1.d0
        do j = i + 1, n
            ainv(j, i) = 0.d0
        end do
    end do
    !     fisrt permutation of ainv
    do i = n, 1, -1
        work(1:n) = ainv(i, 1:n)
        ainv(i, 1:n) = ainv(ipiv(i), 1:n)
        ainv(ipiv(i), 1:n) = work(1:n)
    end do
    !      SECOND
    if (UPLO .eq. 'u' .or. UPLO .eq. 'U') then
        do i = 1, N - 1
            work(i) = a(i, i + 1)
            work(n + i - 1) = -a(i, i + 1)
        end do
        do j = 2, n - 1
            a(1:j - 1, j) = a(1:j - 1, j + 1)
        end do
        a(1:n - 1, n) = 0.d0
    else
        do i = 1, N - 1
            work(i) = -a(i + 1, i)
            work(n + i - 1) = a(i + 1, i)
        end do
        do j = n - 1, 2, -1
            a(j + 1:n, j) = a(j + 1:n, j - 1)
        end do
        a(2:n, 1) = 0.d0
    end if
    work(2*N - 1:3*N - 2) = 0.d0 ! diagonal elements of skew matrix , obviously set to zero.
    call DTRTRS(UPLO, 'N', 'U', N, N, A, LDA, ainv, LDINV, INFO)
    ! USE STANDARD LAPACK DGTSVX or the simpler DGTSV?
    if (yeslap) then
        call DGTSVX('N', 'N', N, N, work(N), work(2*N - 1), work, work(3*N)&
                &, work(4*N), work(5*N), work(6*N), IPIV(N + 1), AINV, LDA, work(7*N), N&
                &, RCOND, work(7*N + N*N - 1), work(8*N + N*N - 1)&
                &, work(9*N + N*N - 1), IPIV(2*N + 1), INFO)
        do i = 1, N
            ainv(1:N, i) = work(7*N + (i - 1)*N:7*N + i*N - 1)
        end do
    else
        call DSKTRS(UPLO, N, N, WORK, AINV, LDINV, WORK(N), INFO)
    end if
    call DTRTRS(UPLO, 'T', 'U', N, N, A, LDA, ainv, LDINV, INFO)
    !     last permutation of ainv
    do i = 1, n
        work(1:n) = ainv(i, 1:n)
        ainv(i, 1:n) = ainv(ipiv(i), 1:n)
        ainv(ipiv(i), 1:n) = work(1:n)
    end do

    return
end

