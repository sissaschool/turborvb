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

program test_zgemv
    implicit none

    complex*16, allocatable, dimension(:, :) :: B
    complex*16, allocatable, dimension(:) :: A, C, C_orig
    real*8, allocatable, dimension(:, :) :: helper_r, helper_c
    complex*16 :: one = 1.d0, zero = 0.d0
    integer :: s, gen
    character(len=1) :: operation
    logical :: yes_on_target

    yes_on_target = .false.

#if defined(_OFFLOAD) && defined(_CUBLAS)
    yes_on_target = .true.
#endif

    ! gen = 1 : Generate matrices
    ! gen = 0 : Compare matrices

    read (*, *) s, operation, gen

    allocate (A(s))
    allocate (B(s, s))
    allocate (C(s))
    allocate (C_orig(s))
    allocate (helper_r(s, s))
    allocate (helper_c(s, s))

    if (gen .eq. 1) then
        call random_number(helper_r(:, 1))
        call random_number(helper_c(:, 1))
        A = cmplx(helper_r(:, 1), helper_c(:, 1))
        call random_number(helper_r)
        call random_number(helper_c)
        B = cmplx(helper_r, helper_c)
    else
        open (unit=10, form="unformatted", file="A", action="read")
        read (10) A
        close (10)
        open (unit=10, form="unformatted", file="B", action="read")
        read (10) B
        close (10)
        open (unit=10, form="unformatted", file="C", action="read")
        read (10) C_orig
        close (10)
    end if

    C = 0

#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp target data map(to:A,B)
!$omp target data map(from:C)
#endif
    call zgemv_(operation&
               &, s, s&
               &, one&
               &, B, s&
               &, A, 1&
               &, zero&
               &, C, 1&
               &, yes_on_target)
#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp end target data
!$omp end target data
#endif

    if (gen .eq. 1) then
        open (unit=10, form="unformatted", file="A", action="write")
        write (10) A
        close (10)
        open (unit=10, form="unformatted", file="B", action="write")
        write (10) B
        close (10)
        open (unit=10, form="unformatted", file="C", action="write")
        write (10) C
        close (10)
    else
        C = C - C_orig
        if (maxval(abs(C)) > 1.0d-10) then
            print *, "ERROR"
        else
            print *, "OK"
        end if
    end if

    deallocate (A, B, C, C_orig, helper_r, helper_c)

end program test_zgemv
