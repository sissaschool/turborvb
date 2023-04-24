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

program test_dgemm
    implicit none

    real*8, allocatable, dimension(:, :) :: A, B, C, C_orig
    real*8, parameter :: one = 1.d0, zero = 0.d0
    real*8 :: time_sec = 0
    integer*4 :: s, gen, ii, repeats
    integer*8 :: T_start, T_end, rate, N_ops
    character(len=1) :: first, second
    character(len=100) :: message, str_tmp1, str_tmp2

    ! gen = 1 : Generate matrices
    ! gen = 0 : Compare matrices

    read (*, *) s, first, second, gen

    allocate (A(s, s))
    allocate (B(s, s))
    allocate (C(s, s))
    allocate (C_orig(s, s))

    if (gen .gt. 0) then
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

    rate = 1
    if (gen .lt. 0 .or. gen .gt. 1) then
        repeats = abs(gen)
        call system_clock(count_rate=rate)
    else
        repeats = 1
    end if

    C = 0

#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp target data map(to:A,B)
!$omp target data map(from:C)
#endif
    if (gen .lt. 0 .or. gen .gt. 1) then
        ! If speed test is on, do 10% extra warm up steps
        do ii = 1, repeats/10
            call dgemm_(first, second&
                       &, s, s, s&
                       &, one&
                       &, A, s&
                       &, B, s&
                       &, zero&
                       &, C, s)
        end do
    end if
    call system_clock(T_start)
    do ii = 1, repeats
        call dgemm_(first, second&
                   &, s, s, s&
                   &, one&
                   &, A, s&
                   &, B, s&
                   &, zero&
                   &, C, s)
    end do
    call system_clock(T_end)
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
    else if (gen .eq. 0) then
        C = C - C_orig
        if (maxval(C) > 1.0d-10) then
            print *, "ERROR"
        else
            print *, "OK"
        end if
    else
        time_sec = (T_end - T_start)/real(rate)
        N_ops = s
        N_ops = N_ops*N_ops*(2*N_ops + 2)
        N_ops = N_ops*repeats

        write (str_tmp1, "(I5)") repeats
        str_tmp1 = adjustl(str_tmp1)
        write (str_tmp2, *) "dgemm_ ", trim(str_tmp1), "x", first, second
        write (str_tmp1, "(I5)") s
        str_tmp1 = adjustl(str_tmp1)
        write (message, *) trim(str_tmp2), "(", trim(str_tmp1), "x", trim(str_tmp1), "):"
        str_tmp2 = message
        write (str_tmp1, "(10F5.2)") time_sec
        str_tmp1 = adjustl(str_tmp1)
        write (message, *) trim(str_tmp2), "  ", trim(str_tmp1), " sec"
        str_tmp2 = message
        write (str_tmp1, "(10F6.2)") N_ops/(1000**3*time_sec)
        str_tmp1 = adjustl(str_tmp1)
        write (message, *) trim(str_tmp2), "  (", trim(str_tmp1), " GFlops)"
        write (6, *) trim(message)

    end if

    deallocate (A, B, C, C_orig)

end program test_dgemm
