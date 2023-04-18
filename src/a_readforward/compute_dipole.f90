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

module Dipole_module

    implicit none
    logical :: ifdipole ! flag for dipole moment calculation
    integer :: ndipole ! the number of dipole components
    integer :: nquad ! the dimension of the quadrupole matrix
    real(8), dimension(:), allocatable :: nuclear_dipole, dipole
    real(8), dimension(:), allocatable :: cc_chg, quad_diag
    real(8), dimension(:, :), allocatable :: quad

contains

    subroutine calc_dipole(r_el, vdim, ell)

        use allio, only: rion, zetar, nion, nel
        implicit none

        real(8), intent(in) :: r_el(3, *)
        real(8), intent(in) :: ell(*)
        integer, intent(in) :: vdim(*)
        real(8) :: r_tmp(3)
        integer :: n, ll, l
        logical :: check_inside

        dipole(1:ndipole) = nuclear_dipole(1:ndipole)

        do n = 1, nel

            dipole(1:ndipole) = dipole(1:ndipole) - r_el(1:ndipole, n)

        end do

    end subroutine calc_dipole

    subroutine calc_quadrupole(r_el, vdim, ell)

        ! to be modified!!!

        use allio, only: rion, zetar, nion, nel
        implicit none

        real(8), intent(in) :: r_el(3, *)
        integer, intent(in) :: vdim(*)
        real(8), intent(in) :: ell(*)
        integer :: i, j, n
        real(8) :: qtot, dist, cion(3, nion), cel(3, nel), sumw, nuc_charge

        !!Calculation of the center of the nuclear charge
        !!(use this until the problem of the asymmetry of the electronic charge is solved)
        !!Warning: this way to calculate the quadrupole is correct only if the dipole of the system is zero
        cc_chg = 0.d0
        do i = 1, nion
            cc_chg(:) = cc_chg(:) + zetar(i)*rion(:, i)
        end do
        cc_chg(:) = cc_chg(:)/sum(zetar(:))

        ! Reference frame of center of charge
        do i = 1, nion
            cion(:, i) = rion(:, i) - cc_chg(:)
        end do
        do i = 1, nel
            cel(:, i) = r_el(:, i) - cc_chg(:)
        end do

        !Formula of quadrupole moment is as following by E.Coccia
        !Q_ij = \sum_n^nel 0.5*(q_k*3*x_in*x_jn - r_n^2d_ij)
        !Nuclear contribution to the quadrupole moment
        quad = 0.d0
        do n = 1, nion
            dist = cion(1, n)**2 + cion(2, n)**2 + cion(3, n)**2
            do i = 1, 3
                do j = i, 3
                    if (i .eq. j) then
                        quad(i, j) = quad(i, j) + zetar(n)*(3.d0*cion(i, n)*cion(j, n) - dist)
                    else
                        quad(i, j) = quad(i, j) + zetar(n)*3.d0*cion(i, n)*cion(j, n)
                        quad(j, i) = quad(i, j)
                    end if
                end do
            end do
        end do

        !Electronic contribution to the quadrupole moment
        do n = 1, nel
            dist = cel(1, n)**2 + cel(2, n)**2 + cel(3, n)**2
            do i = 1, 3
                do j = i, 3
                    if (i .eq. j) then
                        quad(i, j) = quad(i, j) - (3.d0*cel(i, n)*cel(j, n) - dist)
                    else
                        quad(i, j) = quad(i, j) - 3.d0*cel(i, n)*cel(j, n)
                        quad(j, i) = quad(i, j)
                    end if
                end do
            end do
        end do

        quad = 0.5d0*quad
        return
    end subroutine

    subroutine diag_quad

        integer :: i, j, n, INFO
        real*8 :: E(2), TAU(2), WORK(10)

        call DSYTRD('U', nquad, quad, nquad, quad_diag, E, TAU, WORK, 10, INFO)
        call DSTERF(nquad, quad_diag, E, INFO)

    end subroutine diag_quad

end module Dipole_module
