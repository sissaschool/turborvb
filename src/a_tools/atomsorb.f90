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

!> @brief Module for handling atomic orbital parameters and configurations
!>
!> This module defines data structures and parameters for representing atomic orbitals,
!> their configurations, and related properties in quantum chemistry calculations.
module atomsorb
    implicit none

    !> Maximum number of orbitals per atom
    integer, parameter :: maxorb = 20
    
    !> Maximum number of parameters for each orbital
    integer, parameter :: maxpar = 10
    
    !> Maximum number of shells per atom
    integer, parameter :: max_shell = 5

    !> Dimension for lambda arrays (maxorb * max_shell)
    integer, parameter :: lambda_dim = maxorb*max_shell

    !> @brief Type definition for orbital parameters and properties
    !>
    !> Contains optimization flags, number of parameters, multiplicity and parameter values
    type orbital
        integer ioptorb                    !< Orbital optimization flag
        integer npar                       !< Number of parameters
        integer mult                       !< Multiplicity
        double precision parm(maxpar)      !< Array of orbital parameters
        logical fixparm(maxpar)            !< Flags for fixed (non-optimizable) parameters
    end type orbital

    !> @brief Type definition for atomic configuration and orbitals
    !>
    !> Contains atomic properties, orbital configurations, and interaction parameters
    type single_atom
        integer kion                       !< Ion index
        integer nconf                      !< Number of configurations
        integer norb                       !< Number of pairing orbitals
        integer norbj                      !< Number of three-body orbitals

        integer shell_pos(maxorb)          !< Orbital positions in shell
        integer jshell_pos(maxorb)         !< J-orbital positions in shell
        integer nlambda, njlambda          !< Lambda and J-lambda dimensions

        integer map_lambda(lambda_dim)      !< Lambda mapping array
        integer map_jlambda(lambda_dim)     !< J-lambda mapping array

        double precision jonebody(lambda_dim) !< One-body J parameters

        !> Onsite interaction parameters for lambda orbitals
        double precision onsite_lambda(lambda_dim, lambda_dim)
        !> Onsite interaction parameters for j-lambda orbitals
        double precision onsite_jlambda(lambda_dim, lambda_dim)
        
        type(orbital) orb_list(maxorb)     !< List of pairing orbitals
        type(orbital) jorb_list(maxorb)    !< List of three-body orbitals
    end type single_atom

    !> Array of atomic configurations
    type(single_atom), allocatable :: atom_list(:)
end module atomsorb
