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

!> @brief Module defining derived data types for TurboRVB
!>
!> This module contains the definition of various derived data types used
!> throughout the TurboRVB codebase for managing arrays, ion components,
!> k-point grids, and wave function factors.
!>
!> Types
!> -----
!> ARRAY_INT : Dynamic integer array
!>     Contains an allocatable integer array for flexible storage.
!>
!> ION_COMP : Ion component information
!>     Stores multiplicity and arrays of ion indices and components.
!>
!> NKGRID : k-point grid structure
!>     Manages k-point grid dimensions and processing status.
!>
!> WF_FACTOR : Wave function factor object
!>     Contains parameters and matrices for Jastrow/determinant factors.
!>
!> Notes
!> -----
!> - All types use allocatable arrays for dynamic memory management.
!> - The WF_FACTOR type is used for both Jastrow and determinant factors.
!> - Provides a structured approach to data organization in TurboRVB.
!>
!> Example
!> -------
!> Used throughout TurboRVB for structured data management and wave function handling.
module types
    implicit none

!> @brief Dynamic integer array type
!>
!> A simple wrapper for an allocatable integer array, providing
!> a structured way to handle dynamic integer arrays.
    type array_int
        integer, allocatable :: col(:)  !< Allocatable integer array
    end type array_int

!> @brief Ion component information type
!>
!> Stores information about ion components including multiplicity
!> and arrays of ion indices and their corresponding components.
    type ion_comp
        integer :: mult                 !< Multiplicity of the component
        integer, allocatable :: ion(:)  !< Array of ion indices
        integer, allocatable :: comp(:) !< Array of component values
    end type ion_comp

!> @brief k-point grid structure type
!>
!> Manages k-point grid information including grid dimensions,
!> k-point indices, and processing status flags.
    type nkgrid
        integer :: dimshell             !< Dimension of the shell
        integer, allocatable :: kpip(:, :) !< k-point indices (2D array)
        logical, allocatable :: tobedone(:) !< Processing status flags
    end type nkgrid

!> @brief Wave function factor object type
!>
!> Comprehensive data structure for storing wave function factor
!> parameters and matrices. Used for both Jastrow and determinant factors.
!> Contains two-body parameters, exponents, and various matrices for
!> different spin configurations.
    type wf_factor
        real(8), allocatable :: twobody_par(:)    !< Two-body parameters
        real(8), allocatable :: exps(:)           !< Exponent parameters
        real(8), allocatable :: bas_mat(:)        !< Basis matrix
        real(8), allocatable :: exp_mat(:)        !< Exponent matrix
        real(8), allocatable :: exp_mat_sz(:)     !< Exponent matrix for Sz
        real(8), allocatable :: exp_mat_c(:)      !< Complex exponent matrix
        real(8), allocatable :: exp_mat_sz_c(:)   !< Complex exponent matrix for Sz
    end type

end module types
