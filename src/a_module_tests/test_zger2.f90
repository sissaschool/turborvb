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

program test_zger2

#if defined(_OFFLOAD)
    use constants, only: yes_ontarget
#endif

    implicit none

    complex*16, allocatable, dimension(:, :) :: A, A_orig
    complex*16, allocatable, dimension(:, :) :: X, Y
    real*8, allocatable, dimension(:, :) :: helper_r, helper_c
    complex*16 :: one = 1.d0, zero = 0.d0
    integer :: s, ii, jj

#if defined(_OFFLOAD)
    yes_ontarget = .true.
#endif

    read (*, *) s

    allocate (A(s, s))
    allocate (A_orig(s, s))
    allocate (X(s, 2))
    allocate (Y(s, 2))
    allocate (helper_r(s, s))
    allocate (helper_c(s, s))

    call random_number(helper_r)
    call random_number(helper_c)
    A = cmplx(helper_r, helper_c)
    call random_number(helper_r(:, 1:2))
    call random_number(helper_c(:, 1:2))
    X = cmplx(helper_r(:, 1:2), helper_c(:, 1:2))
    call random_number(helper_r(:, 1:2))
    call random_number(helper_c(:, 1:2))
    y = cmplx(helper_r(:, 1:2), helper_c(:, 1:2))

    A_orig = A

#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp target data map(tofrom:A, X, Y)
#endif
    call zger2(s, s&
              &, one&
              &, X, s&
              &, Y, s&
              &, A&
              &, s)
#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp end target data
#endif

    ! Permutate rows on the original matrix
    do ii = 1, s
        do jj = 1, s
            A(jj, ii) = A(jj, ii) - one*(X(jj, 1)*Y(ii, 1) + X(jj, 2)*Y(ii, 2))
        end do
    end do

    ! Compare product with original matrix
    A = A_orig - A
    if (maxval(abs(A)) > 1.0d-10) then
        print *, "ERROR", maxval(abs(A))
    else
        print *, "OK"
    end if

    deallocate (A, A_orig, X, Y, helper_r, helper_c)

end program test_zger2
