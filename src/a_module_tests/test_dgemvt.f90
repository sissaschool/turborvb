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

program test_dgemvt
    implicit none

    real*8, allocatable, dimension(:, :) :: B
    real*8, allocatable, dimension(:) :: A, C, C_orig
    real*8 :: one = 1.d0, zero = 0.d0
    integer :: s, gen

    ! gen = 1 : Generate matrices
    ! gen = 0 : Compare matrices

    read (*, *) s, gen

    allocate (A(s))
    allocate (B(s, s))
    allocate (C(s))
    allocate (C_orig(s))

    if (gen .eq. 1) then
        call random_number(A)
        call random_number(B)
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

#if defined(_OFFLOAD)
!$omp target data map(to:A,B)
!$omp target data map(from:C)
#endif
    call dgemvt_offload(s, s&
                       &, one&
                       &, B, s&
                       &, A, 1&
                       &, zero&
                       &, C, 1)
#if defined(_OFFLOAD)
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
        if (maxval(C) > 1.0d-10) then
            print *, "ERROR"
        else
            print *, "OK"
        end if
    end if

    deallocate (A, B, C, C_orig)

end program test_dgemvt
