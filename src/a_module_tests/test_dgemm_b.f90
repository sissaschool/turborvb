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

program test_dgemm_b
    implicit none

    real*8, allocatable, dimension(:, :) :: A, B, AB, BB, CB, CB_i
    real*8, allocatable, dimension(:, :) :: AB_orig, BB_orig, CB_orig
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
    allocate (AB(s, s))
    allocate (BB(s, s))
    allocate (CB(s, s))
    allocate (CB_i(s, s))
    allocate (AB_orig(s, s))
    allocate (BB_orig(s, s))
    allocate (CB_orig(s, s))

    if (gen .gt. 0) then
        call random_number(A)
        call random_number(B)
        call random_number(CB)
    else
        open (unit=10, form="unformatted", file="A", action="read")
        read (10) A
        close (10)
        open (unit=10, form="unformatted", file="B", action="read")
        read (10) B
        close (10)
        open (unit=10, form="unformatted", file="CB", action="read")
        read (10) CB
        close (10)
        open (unit=10, form="unformatted", file="AB", action="read")
        read (10) AB_orig
        close (10)
        open (unit=10, form="unformatted", file="BB", action="read")
        read (10) BB_orig
        close (10)
        open (unit=10, form="unformatted", file="CB_out", action="read")
        read (10) CB_orig
        close (10)
    end if

    A(1, 1) = 0

    rate = 1
    if (gen .lt. 0 .or. gen .gt. 1) then
        repeats = abs(gen)
        call system_clock(count_rate=rate)
    else
        repeats = 1
    end if

    AB = 0
    BB = 0
    CB_i = CB

#if defined(_OFFLOAD) && defined(_CUBLAS)
!$omp target data map(to:A,B)
!$omp target data map(tofrom:AB,BB,CB)
#endif
    if (gen .lt. 0 .or. gen .gt. 1) then
        ! If speed test is on, do 10% extra warm up steps
        do ii = 1, repeats/10
            call dgemm_b(first, second&
                        &, s, s, s&
                        &, one&
                        &, A, s&
                        &, AB, s&
                        &, B, s&
                        &, BB, s&
                        &, zero&
                        &, CB, s)
        end do
    end if
    call system_clock(T_start)
    do ii = 1, repeats
        call dgemm_b(first, second&
                    &, s, s, s&
                    &, one&
                    &, A, s&
                    &, AB, s&
                    &, B, s&
                    &, BB, s&
                    &, zero&
                    &, CB, s)
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
        open (unit=10, form="unformatted", file="AB", action="write")
        write (10) AB
        close (10)
        open (unit=10, form="unformatted", file="BB", action="write")
        write (10) BB
        close (10)
        open (unit=10, form="unformatted", file="CB", action="write")
        write (10) CB_i
        close (10)
        open (unit=10, form="unformatted", file="CB_out", action="write")
        write (10) CB
        close (10)
    else if (gen .eq. 0) then
        AB = AB - AB_orig
        BB = BB - BB_orig
        CB = CB - CB_orig
        if (maxval(abs(AB)) > 1.0d-10) then
            print *, "ERROR AB", maxval(abs(AB))
        end if
        if (maxval(abs(BB)) > 1.0d-10) then
            print *, "ERROR BB", maxval(abs(BB))
        end if
        if (maxval(abs(CB)) > 1.0d-10) then
            print *, "ERROR CB", maxval(abs(CB))
        end if
        if (maxval(abs(AB)) <= 10d-10&
           &.and. maxval(abs(BB)) <= 10d-10&
           &.and. maxval(abs(CB)) <= 10d-10) then
            print *, "OK"
        end if
    else
        time_sec = (T_end - T_start)/real(rate)
        N_ops = s
        N_ops = N_ops*N_ops*(2*N_ops + 2) ! First dgemm
        N_ops = N_ops*2 ! There are actually two dgemms
        N_ops = N_ops + s*s ! Plus one scale
        N_ops = N_ops*repeats ! Times number of repeats

        write (str_tmp1, "(I5)") repeats
        str_tmp1 = adjustl(str_tmp1)
        write (str_tmp2, *) "dgemm_b ", trim(str_tmp1), "x", first, second
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

    deallocate (A, B, AB, BB, CB, CB_i, AB_orig, BB_orig, CB_orig)

end program test_dgemm_b
