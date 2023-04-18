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

module symm_data
    implicit none

    !***** Symmetries variables
    integer :: isymm(3, 3, 48), nrot, isymm_full(3, 3, 48), isymm_norot(3, 3, 48)
    !the symmetry matrices
    !the number of symmetry matrice

    character :: isname(48)*45
    ! full name of the rotational part of each selected symmetry operation

    integer :: nsym
    ! used to find symmetries between parameters

    integer :: nsym_full
    ! output of sgama.f90 : real number of symmetry operations

    ! Indipendent integer translation vectors between the unitary cells
    type translation
        double precision :: vec(3)
        double precision :: tsign
        integer :: z(2)
    end type translation

    real*8 at_qe(3, 3), car2cry(3, 3)

    ! Coordinated of all unitary cells in reduced units
    double precision, allocatable :: pos(:, :)

    ! Map the distance between two cells in another one usint a symmetry
    integer, pointer :: cellmap_jas(:, :, :), cellmap_det(:, :, :)
    ! cellmap(1,i1,i2) = i3     index referes to a translation tral(:,i3)
    ! cellmap(2,i1,i2) = isymm  index referes to a symmetry (with APBC sing)
    ! cellmap(3,i1,i2) = index referes to a translation without any symm
    ! cellmap(4,i1,i2) = record the sign in APBC case

    double precision :: dist(3), test(3), mtest(3)

    ! Flag to not use symmetries
    logical :: nosym, notra, nosym_contr, forces_sym, nosym_forces, notra_forces

    ! Inversion symmetry flag
    !logical :: invsym

    ! Symmetry rotation flag
    !logical :: norotation

    ! Not use time-reversal
    !logical :: time_reversal

    double precision, parameter :: eps = 1e-5

    ! if .false. not consider all atoms inside the basic cell different
    logical :: eqatoms, eq_intatoms
    !******* NOT USED VARIABLES ********

    !********** Log file unit ************
    integer, parameter :: lunit = 25
    logical :: write_log
contains

    subroutine init_cell
        implicit none
        integer ipiv(3), i, info
        real*8 matscra(9)
        car2cry = at_qe
        call dgetrf(3, 3, car2cry, 3, ipiv, info)
        if (info .ne. 0) then
            write (6, *) ' ERROR in initialization cell (dgetrf) !!! ', info
            do i = 1, 3
                write (6, *) i, ipiv(i), car2cry(i, i)
            end do
        else
            call dgetri(3, car2cry, 3, ipiv, matscra, 9, info)
            if (info .ne. 0) write (6, *) ' ERROR in initialization cell (dgetri) !!! '
        end if
    end subroutine init_cell
    subroutine CartesianToCrystal(r, howmany)
        integer, intent(in) :: howmany
        double precision, dimension(3, howmany), intent(inout) :: r
        double precision s(3)
        integer i
        do i = 1, howmany
            s(:) = r(:, i)
            call dgemv('N', 3, 3, 1.d0, car2cry, 3, s, 1, 0.d0, r(1, i), 1)
        end do
    end subroutine CartesianToCrystal
end module symm_data
