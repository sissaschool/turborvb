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

module atomsorb
    implicit none

    !!**  Max number of orbitals for atom
    integer, parameter :: maxorb = 20
    !
    !!**  Max total number of parameters for each orbitals
    integer, parameter :: maxpar = 10
    !
    !!*** Max number of shell for atoms
    integer, parameter :: max_shell = 5

    integer, parameter :: lambda_dim = maxorb*max_shell

    type orbital
        integer ioptorb
        integer npar
        integer mult
        double precision parm(maxpar)
        ! *** logical flag for parameter that has no to be optimized
        logical fixparm(maxpar)
    end type orbital

    type single_atom
        integer kion
        !*** number of configuration for an atom
        integer nconf
        !*** number of pairing orbitals for a given configuration and atom
        integer :: norb
        !*** number of threebody orbitals for a given configuration and atom
        integer norbj

        !***** Position of the orbital in the shell
        integer shell_pos(maxorb)
        integer jshell_pos(maxorb)
        integer nlambda, njlambda

        integer map_lambda(lambda_dim), map_jlambda(lambda_dim)

        double precision :: jonebody(lambda_dim)

        double precision :: onsite_lambda(lambda_dim, lambda_dim)
        double precision :: onsite_jlambda(lambda_dim, lambda_dim)
        !****** Pairing orbitals for each configuration
        type(orbital) :: orb_list(maxorb)
        !****** Three-body orbitals for each configuration
        type(orbital) :: jorb_list(maxorb)
    end type single_atom

    type(single_atom), allocatable :: atom_list(:)
end module atomsorb
