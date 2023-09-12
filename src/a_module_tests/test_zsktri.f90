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

program test_zsktri
    implicit none

    complex*16, allocatable, dimension(:, :) :: A, A_inv, A_inv_orig
    complex*16, allocatable, dimension(:) :: W
    integer, allocatable, dimension(:) :: ipiv
    real*8, allocatable, dimension(:, :) :: helper_r, helper_c
    complex*16 :: one = 1.d0, zero = 0.d0
    integer :: s, gen, ii, jj, info
    character(len=1) :: uplo
    
    ! s dimension of the test matrix, s x s.
    ! uplo: U->upper triangular, L->lower triangular
    ! gen = 0 : Compare matrices, gen = 1 : Generate matrices
    write (*, *) 's, uplo, gen'
    read (*, *) s, uplo, gen

    allocate (ipiv(s))
    allocate (A(s, s))
    allocate (A_inv(s, s))
    allocate (A_inv_orig(s, s))
    allocate (W(s**2+12*s-2))
    allocate (helper_r(s, s))
    allocate (helper_c(s, s))
    
    if (gen .eq. 1) then

        ! generate a skew_symmetric matrix A
        call random_number(helper_r)
        call random_number(helper_c)
    
        do ii = 1, s
            do jj = 1, ii - 1
                helper_c(jj, ii) = -helper_c(ii, jj)
                helper_r(jj, ii) = -helper_r(ii, jj)
            end do
        end do
        
        do ii = 1, s
            helper_r(ii, ii) = 0
            helper_c(ii, ii) = 0
        end do
    
        A = cmplx(helper_r, helper_c)

    else
        open (unit=10, form="unformatted", file="A", action="read")
        read (10) A
        close (10)
        open (unit=10, form="unformatted", file="A_inv_orig", action="read")
        read (10) A_inv_orig
        close (10)
    end if

    W = 0

    do ii = 1, s
        ipiv(ii)=ii
    end do

    call zsktri(uplo, s, A, s, A_inv, s, ipiv, W, info)

    if (gen .eq. 1) then
        open (unit=10, form="unformatted", file="A", action="write")
        write (10) A
        close (10)
        open (unit=10, form="unformatted", file="A_inv_orig", action="write")
        write (10) A_inv_orig
        close (10)
    else
        A_inv = A_inv - A_inv_orig
        if (maxval(abs(A_inv)) > 1.0d-10) then
            print *, "ERROR"
        else
            print *, "OK"
        end if
    end if
    
    deallocate (A, A_inv, A_inv_orig, W, helper_r, helper_c)

end program test_zsktri
