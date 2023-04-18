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

program test_dtrsm
    implicit none

    real*8, allocatable, dimension(:, :) :: A, B, C, C_orig
    real*8 :: one = 1.d0, zero = 0.d0
    integer :: s, gen, ii, jj
    character(len=1) :: uplo, side, trans

    ! gen = 1 : Generate matrices
    ! gen = 0 : Compare matrices

    read (*, *) s, uplo, side, trans, gen

    allocate (A(s, s))
    allocate (B(s, s))
    allocate (C(s, s))
    allocate (C_orig(s, s))

    if (gen .eq. 1) then
        call random_number(A)
        call random_number(B)

        ! A is upper and lower unitriangular matrix
        do ii = 1, s
            do jj = 1, ii - 1
                A(jj, ii) = A(ii, jj)
            end do
        end do
        do ii = 1, s
            A(ii, ii) = 1
        end do
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

    C = B

#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp target data map(to:A)
!$omp target data map(tofrom:C)
#endif
    call dtrsm_(side, uplo, trans, "U"&
               &, s, s&
               &, one&
               &, A, s&
               &, C, s)
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
        if (maxval(abs(C)) > 1.0d-7) then
            print *, "ERROR", maxval(abs(C))
        else
            print *, "OK"
        end if
    end if

    deallocate (A, B, C, C_orig)

end program test_dtrsm
