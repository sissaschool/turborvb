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

program test_upwinvp_complex
    use constants, only: yes_ontarget
    implicit none
    integer, parameter :: nel = 3, indt = 4
    complex*16 :: psi(indt, nel, 2), ainv(nel), ainvn(nel), winv(nel, indt), expected(nel, indt)
    integer :: i, j, k
    real*8, parameter :: eps = 1.0e-10 ! A small value for numerical comparison

    yes_ontarget = .true.

    ! sleep to avoid unexpected corruptions.
    call sleep(1)

    ! Initialize test data
    winv = cmplx(0.0, 0.0)
    do i = 1, indt
        do j = 1, nel
            do k = 1, 2
                psi(i, j, k) = cmplx(i*j*k, i*j*k)
            end do
        end do
    end do
    ainv = (/cmplx(1.0, 1.0), cmplx(2.0, 2.0), cmplx(3.0, 3.0)/)
    ainvn = (/cmplx(4.0, 4.0), cmplx(5.0, 5.0), cmplx(6.0, 6.0)/)

    call upwinvp_complex(nel, indt, winv, ainv, ainvn, psi)

    ! Calculate expected result
    do i = 1, indt
        do j = 1, nel
            expected(j, i) = ainv(j)*psi(i, j, 1) + ainvn(j)*psi(i, j, 2)
        end do
    end do

    ! Check the result
    do i = 1, indt
        do j = 1, nel
            if (abs(winv(j, i) - expected(j, i)) > eps) then
                print *, "Test failed at i=", i, ", j=", j
                print *, "OK"
                stop 1
            end if
        end do
    end do

    print *, "OK"
    print *, "All tests passed."
end program test_upwinvp_complex
