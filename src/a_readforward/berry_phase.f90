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

module berry_phase

    real(8), allocatable :: berry_comp_vec(:)
    logical ifberry
    real*8 berry_exp(6)

contains

    ! This subroutine allocate and initialise
    ! the constant vector with reciprocal cell dimensions
    ! used in z_N caluculations
    ! nel: # of electrons
    subroutine init_berry_phase(nel)

        use cell, only: cellscale
        use constants, only: TWO_PI
        use allio, only: rank
        implicit none

        integer, intent(in) :: nel
        integer :: i

        ! This vector contains nel copies of the 3D vector
        ! made by the orthogonal dimensions of the reciprocal cell
        ! BUGFIX: true for orthorombic cells only!!!

        if (allocated(berry_comp_vec)) deallocate (berry_comp_vec)
        allocate (berry_comp_vec(3*nel))
        do i = 1, nel
            berry_comp_vec(3*(i - 1) + 1) = TWO_PI/cellscale(1)
            berry_comp_vec(3*(i - 1) + 2) = TWO_PI/cellscale(2)
            berry_comp_vec(3*(i - 1) + 3) = TWO_PI/cellscale(3)
        end do

        !if( rank.eq.0 ) then
        !  write(6,*) '[Berry] ---> initialise berry_comp_vec'
        !  write(6,*) berry_comp_vec(:)
        !  write(6,*) '---> end'
        !endif

    end subroutine

    ! This subroutine compute z_N (here N = nel)
    ! according to Resta & Sorella PRL 82, 370 (1999)
    ! kelw:       electronic positions for a single walker
    ! nel:        # electrons
    ! nw:         # walkers
    ! wconf:      walker's weight
    ! berry_exp:  z_N

    subroutine berry_phase_update(kelw, nel)

        implicit none

        real*8, intent(in) :: kelw(3*nel)
        integer, intent(in) :: nel
        real*8 berry_exp_dummy
        real*8 DDOT
        integer i

        ! write(6,*) '---> berry_phase arguments'
        ! write(6,*) 'kelw      ', kelw
        ! write(6,*) 'nel       ', nel
        ! write(6,*) 'wconf     ', wconf
        ! write(6,*) 'berry_exp ', berry_exp
        ! write(6,*) '---> end'

        ! compute k.x, the argumet of the exponential
        ! berry_comp_vec defined in cell.f90
        do i = 1, 3
            berry_exp_dummy = DDOT(nel, berry_comp_vec(i), 3, kelw(i), 3)

            ! compute e^{i*k.x}
            berry_exp(2*i - 1) = dcos(berry_exp_dummy)
            berry_exp(2*i) = dsin(berry_exp_dummy)

        end do

        ! write(6,*) 'berry_exp [updated]', berry_exp

    end subroutine

end module
