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

program test_dgetrf

#if defined(_OFFLOAD) && defined(_CUSOLVER)
    use allio, only: handle, dev_dgetrf_workspace, dev_Info
#endif

    implicit none

    real*8, allocatable, dimension(:, :) :: A, C, upper, lower, dummy
    real*8, allocatable, dimension(:) :: temp
    integer*4, allocatable, dimension(:) :: ipiv
    real*8 :: one = 1.d0, zero = 0.d0
    integer :: s, gen, info, ii, jj, kk

#if defined(_OFFLOAD) && defined(_CUSOLVER)
    integer :: lworkspace, stat
#endif

    ! gen = 1 : Generate matrices
    ! gen = 0 : Compare matrices

    read (*, *) s, gen

#if defined(_OFFLOAD) && defined(_CUSOLVER)
    lworkspace = 1
#ifdef RISC
    call cusolver_handle_init_(handle)
#else
    call cusolver_handle_init(handle)
#endif

    allocate (dummy(s, s))
!$omp target data map(alloc:dummy)

#ifdef RISC
    call cusolver_dgetrf_buffersize_(handle, stat, s, s&
                                      &, dummy, s, lworkspace)
#else
    call cusolver_dgetrf_buffersize(handle, stat, s, s&
                                      &, dummy, s, lworkspace)
#endif

!$omp end target data
    deallocate (dummy)

    allocate (dev_dgetrf_workspace(lworkspace))

!$omp target data map(alloc:dev_dgetrf_workspace, dev_Info)
#endif

    allocate (A(s, s))
    allocate (ipiv(s))
    allocate (temp(s))
    allocate (C(s, s))
    allocate (upper(s, s))
    allocate (lower(s, s))

    if (gen .eq. 1) then
        call random_number(A)
        open (unit=10, form="unformatted", file="A", action="write")
        write (10) A
        close (10)
        print *, "GENERATED"
        stop 0
    else
        open (unit=10, form="unformatted", file="A", action="read")
        read (10) A
        close (10)
    end if

    C = A

#if defined(_OFFLOAD) && defined(_CUSOLVER)
!$omp target data map(tofrom:C,ipiv)
#endif
    call dgetrf_(s, s&
                &, C, s&
                &, ipiv&
                &, info)
#if defined(_OFFLOAD) && defined(_CUSOLVER)
!$omp end target data
#endif

    if (info .ne. 0) then
        print *, "ERROR non-zero exit status", info
        stop 0
    end if

    ! Permutate rows on the original matrix
    do ii = 1, s
        temp = A(ii, :)
        A(ii, :) = A(ipiv(ii), :)
        A(ipiv(ii), :) = temp
    end do

    ! Create upper and lower matrices
    upper = 0
    lower = 0
    do ii = 1, s
        do jj = 1, s
            if (ii .le. jj) then
                upper(ii, jj) = C(ii, jj)
            end if
            if (ii .gt. jj) then
                lower(ii, jj) = C(ii, jj)
            end if
            if (ii .eq. jj) then
                lower(ii, jj) = 1.0
            end if
        end do
    end do

    ! Multiply lower and upper matrix
    C = 0
    do ii = 1, s
        do jj = 1, s
            do kk = 1, s
                C(ii, jj) = C(ii, jj) + lower(ii, kk)*upper(kk, jj)
            end do
        end do
    end do

    ! Compare product with original matrix
    C = C - A
    if (maxval(C) > 1.0d-6) then
        print *, "ERROR"
    else
        print *, "OK"
    end if

#if defined(_OFFLOAD) && defined(_CUSOLVER)
!$omp end target data
#ifdef RISC
    call cusolver_handle_destroy_(handle)
#else
    call cusolver_handle_destroy(handle)
#endif
#endif

    deallocate (A, ipiv, temp, upper, lower, C)

#if defined(_OFFLOAD) && defined(_CUSOLVER)
    deallocate (dev_dgetrf_workspace)
#endif

end program test_dgetrf
