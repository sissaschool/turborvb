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

program test_dger

#if defined(_OFFLOAD)
    use constants, only: yes_ontarget
#endif

    implicit none

    real*8, allocatable, dimension(:, :) :: A, A_orig
    real*8, allocatable, dimension(:) :: X, Y
    real*8, allocatable, dimension(:, :) :: helper_r, helper_c
    real*8 :: one = 1.d0, zero = 0.d0
    integer :: s, ii, jj

#if defined(_OFFLOAD)
    yes_ontarget = .true.
#endif

    read (*, *) s

    allocate (A(s, s))
    allocate (A_orig(s, s))
    allocate (X(s))
    allocate (Y(s))
    allocate (helper_r(s, s))
    allocate (helper_c(s, s))

    call random_number(A)
    call random_number(X)
    call random_number(Y)

    A_orig = A

#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp target data map(tofrom:A, X, Y)
#endif
    call dger_(s, s&
              &, one&
              &, X, 1&
              &, Y, 1&
              &, A&
              &, s)
#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp end target data
#endif

    ! Permutate rows on the original matrix
    do ii = 1, s
        do jj = 1, s
            A(jj, ii) = A(jj, ii) - one*X(jj)*Y(ii)
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

end program test_dger
