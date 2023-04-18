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

program test_zgetrfi

#if defined(_OFFLOAD) && defined(_CUSOLVER)
    use allio, only: handle, dev_zgetrf_workspace&
                  &, dev_zgetri_workspace, dev_Info
#endif

    implicit none

    complex*16, allocatable, dimension(:, :) :: A, C, C_ans, dummy
    complex*16, allocatable, dimension(:) :: temp, work
    real*8, allocatable, dimension(:, :) :: helper_r, helper_c
    integer*4, allocatable, dimension(:) :: ipiv
    complex*16 :: one = 1.d0, zero = 0.d0
    integer :: s, gen, info, ii, jj, kk, lwork

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
    call cusolver_zgetrf_buffersize_(handle, stat, s, s&
                                      &, dummy, s, lworkspace)
#else
    call cusolver_zgetrf_buffersize(handle, stat, s, s&
                                      &, dummy, s, lworkspace)
#endif

!$omp end target data
    deallocate (dummy)

    allocate (dev_zgetri_workspace(s, s))
    allocate (dev_zgetrf_workspace(lworkspace))

!$omp target data map(alloc:dev_zgetrf_workspace, dev_zgetri_workspace&
!$omp                    &, dev_Info)
#endif

    allocate (A(s, s))
    allocate (ipiv(s))
    allocate (temp(s))
    allocate (C(s, s))
    allocate (C_ans(s, s))
    allocate (helper_r(s, s))
    allocate (helper_c(s, s))

    if (gen .eq. 1) then
        call random_number(helper_r)
        call random_number(helper_c)
        A = cmplx(helper_r, helper_c)
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

    call zgetrf(s, s&
               &, C, s&
               &, ipiv&
               &, info)

    if (info .ne. 0) then
        print *, "ERROR non-zero exit status from zgetrf", info
        stop 0
    end if

#if defined(_OFFLOAD) && defined(_CUSOLVER) && defined(_X_)
!$omp target data map(tofrom:C, ipiv)
#endif
    lwork = -1
    allocate (work(1))
    call zgetri(s&
                &, C, s&
                &, ipiv&
                &, work, lwork&
                &, info)
    if (info .ne. 0) then
        print *, "ERROR non-zero exit status from zgetri work size assesment"&
              &, info
        stop 0
    end if
    lwork = int(work(1))
    deallocate (work)
    allocate (work(lwork))
    call zgetri(s&
                &, C, s&
                &, ipiv&
                &, work, lwork&
                &, info)
    deallocate (work)
#if defined(_OFFLOAD) && defined(_CUSOLVER) && defined(_X_)
!$omp end target data
#endif

    if (info .ne. 0) then
        print *, "ERROR non-zero exit status from zgetri_ wrapper", info
        stop 0
    end if

    ! Multiply lower and upper matrix
    !C_ans = 0
    !do ii = 1, s
    !    do jj = 1, s
    !        do kk = 1, s
    !            C_ans(ii, jj) = C_ans(ii, jj) + A(ii, kk) * C(kk, jj)
    !        end do
    !    end do
    !end do
    call zgemm("N", "N", s, s, s, one, A, s, C, s, zero, C_ans, s)

    print *, C_ans
    ! Remove one form the diagonal
    C = C_ans
    do ii = 1, s
        C(ii, ii) = C_ans(ii, ii) - 1
    end do

    if (maxval(abs(C)) > 1.0d-6) then
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

    deallocate (A, ipiv, temp, C, C_ans, helper_r, helper_c)

#if defined(_OFFLOAD) && defined(_CUSOLVER)
    deallocate (dev_zgetrf_workspace, dev_zgetri_workspace)
#endif

end program test_zgetrfi
