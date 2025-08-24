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

!> @brief External potential module for quantum Monte Carlo calculations
!>
!> This module provides functionality for handling external potentials in
!> quantum Monte Carlo calculations, including cube file potentials, QM/MM
!> interfaces, and various molecular mechanics interactions.
!>
!> @details
!> The module contains variables for:
!> - Cube file potential data (grid, atoms, charges)
!> - Statistical averages and error estimates
!> - MPI parallelization support
!> - Link atom definitions for QM/MM interfaces
!> - Molecular mechanics restraints and forces
!>
!> @note Used in VMC and DMC calculations with external potentials
!> @note Supports QM/MM hybrid calculations and molecular mechanics
module extpot
    character(len=80) :: filename_cube, title_cube(2)
    integer :: n_x, n_y, n_z
    integer :: n_atoms
    real*8 :: x0(3), delta(3)
    integer, dimension(:), allocatable :: id_atom
    real*8, dimension(:, :), allocatable :: x_atom
    real*8, dimension(:), allocatable :: chrg_atom
    real*8, dimension(:, :, :), allocatable :: pot
    real*8, dimension(:), allocatable :: xdata, ydata, zdata
    integer :: nout = 0, ncount = 0
    real*8 :: ave = 0, ave2 = 0
    real*8 :: err_el
    ! for MPI runs
    integer :: t_nout = 0, t_ncount = 0
    real*8 :: t_ave = 0, t_ave2 = 0
    logical :: ext_pot
    ! for nuclei
    integer :: ncount_ion = 0, t_ncount_ion = 0
    real*8 :: ave_ion = 0, ave2_ion = 0
    real*8 :: err_ion
    real*8 :: t_ave_ion = 0, t_ave2_ion = 0
    real*8 :: total_ave, total_err
    ! link atoms
    logical :: link_atom
    integer :: latoms
    ! MM restraints
    logical :: mm_restr
    ! write random walk
    logical :: write_rwalk
end module extpot

!> @brief Link atoms module for QM/MM interface
!>
!> This module defines data structures and variables for handling link atoms
!> in quantum mechanics/molecular mechanics (QM/MM) hybrid calculations.
!> Link atoms are used to cap the QM region at the boundary with the MM region.
!>
!> @details
!> The module contains:
!> - Link atom definitions (cap, QM, MM atoms)
!> - Capping atom parameters and flags
!> - Arrays for managing link atom connectivity
!>
!> @note Used in QM/MM hybrid calculations
!> @note Supports various capping schemes for different force fields
module link_atoms
    ! link atoms
    integer :: maxcap
    real*8 calpha(5)
    ! cap   -> capping atom
    ! qm    -> qm-link atom
    ! mm    -> mm-link atom
    type itriplet
        integer :: cap, qm, mm
    end type itriplet
    type(itriplet), dimension(:), allocatable :: capping
    logical, dimension(:), allocatable :: log_cap
    integer, dimension(:), allocatable :: qm, prt, cap
end module link_atoms

!> @brief Link angle module for molecular mechanics interactions
!>
!> This module provides data structures and variables for handling molecular
!> mechanics interactions in QM/MM calculations, including bonds, angles,
!> dihedrals, and improper dihedrals.
!>
!> @details
!> The module contains:
!> - Two-body radial (bond) contributions
!> - Three-body angular contributions
!> - Four-body dihedral contributions
!> - Improper dihedral contributions
!> - Restraint parameters for MM region
!>
!> @note Used in QM/MM calculations with molecular mechanics force fields
!> @note Supports various force field types (GROMOS, AMBER, etc.)
module link_angle
    use constants, only: pi
    character(20) :: filename_link
    !two-body radial contribution
    type ifour
        integer :: i, j
        real*8 :: kbond
        real*8 :: req
    end type ifour
    ! Three-body angular contribution
    type ifive
        integer :: i, j, k
        real*8 :: ktheta
        real*8 :: thetaeq
    end type ifive
    type(ifive), dimension(:, :), allocatable :: linangle
    real*8 :: mm_pot_theta, mm_pot_dihed, mm_pot_impr
    integer :: maxth = 8, maxphi = 28, maxqhi = 3
    integer, dimension(:), allocatable :: ntheta
    ! Proper dihedral contribution
    type iseven
        integer :: i, j, k, l, mult
        real*8 :: kphi
        real*8 :: pcos
    end type iseven
    type(iseven), dimension(:, :), allocatable :: lindhd
    integer, dimension(:), allocatable :: nphi
    ! Improper dihedral contribution
    type isix
        integer :: i, j, k, l
        real*8 :: kqhi
        real*8 :: qcos
    end type isix
    type(isix), dimension(:, :), allocatable :: linimp
    integer, dimension(:), allocatable :: nimp
    ! only if mm_restr=.true.
    type(ifour), dimension(:), allocatable :: cl_bond
    type(ifive), dimension(:), allocatable :: cl_angle
    type(isix), dimension(:), allocatable :: cl_dimp
    type(isix), dimension(:), allocatable :: cl_dihe
end module link_angle

!> @brief Classical restraints module for QM/MM calculations
!>
!> This module provides variables for handling classical restraints in
!> QM/MM calculations, including bond, angle, dihedral, and improper
!> dihedral restraints.
!>
!> @details
!> The module contains:
!> - Number of different types of restraints
!> - Restraint force constants and equilibrium values
!> - Restraint forces for each type of interaction
!> - Scaling factor for MM contributions
!>
!> @note Used in QM/MM calculations with restraint potentials
!> @note Supports harmonic and other restraint functional forms
module cl_restr
    integer :: nbonds, nth, ndihe, ndimp
    real*8 :: mm_fact
    real*8 :: restr_bond, restr_angle, restr_dihe, restr_dimp
    real*8, dimension(:, :), allocatable :: restr_f_bond, restr_f_angle, restr_f_dihe, restr_f_dimp
end module cl_restr

!> @brief Total angle module for molecular mechanics energy averages
!>
!> This module provides variables for accumulating and averaging molecular
!> mechanics energy contributions in QM/MM calculations.
!>
!> @details
!> The module contains:
!> - Counters for different types of interactions
!> - Running averages and squared averages
!> - Error estimates for each interaction type
!> - MPI parallelization support for statistics
!>
!> @note Used for statistical analysis of MM energy contributions
!> @note Supports both serial and parallel calculations
module tot_angle
    ! Potential energy average
    integer :: ncount_bond, ncount_angle, ncount_dihed, ncount_impr
    real*8 :: ave_bond, ave2_bond, ave_angle, ave2_angle, ave_dihed, ave2_dihed, ave_impr, ave2_impr
    integer :: t_ncount_bond, t_ncount_angle, t_ncount_dihed, t_ncount_impr
    real*8 :: t_ave_bond, t_ave2_bond, t_ave_angle, t_ave2_angle, t_ave_dihed, t_ave2_dihed, t_ave_impr, t_ave2_impr
    real*8 :: err_bond, err_angle, err_dihed, err_impr
end module tot_angle

!> @brief Splines module for interpolation of external potentials
!>
!> This module provides data structures and variables for handling spline
!> interpolation of external potentials on 3D grids.
!>
!> @details
!> The module contains:
!> - Spline orders and knot numbers for each dimension
!> - Spline coefficients for 3D interpolation
!> - Knot positions and grid spacing
!> - Current position for interpolation
!>
!> @note Used for efficient evaluation of external potentials
!> @note Supports cubic spline interpolation in 3D
module splines
    integer :: kxord, kyord, kzord, nxknot, nyknot, nzknot
    integer :: nxcoef, nycoef, nzcoef
    real*8 :: x, y, z
    real*8 :: deltavec(3)
    real*8, dimension(:, :, :), allocatable :: bscoef
    real*8, dimension(:), allocatable :: &
            &        xknot, yknot, zknot
end module splines

!> @brief Vector module for 3D vector field interpolation
!>
!> This module provides data structures and variables for handling 3D
!> vector field interpolation and evaluation.
!>
!> @details
!> The module contains:
!> - Grid dimensions for vector field
!> - Vector field values on 3D grid
!> - Coordinate vectors for interpolation
!> - First call flag for initialization
!>
!> @note Used for vector field potentials and forces
!> @note Supports 3D vector interpolation
module vector
    integer :: nxvec, nyvec, nzvecz
    logical :: first_call_vec
    real*8, dimension(:), allocatable :: xvec, yvec, zvec
    real*8, dimension(:, :, :), allocatable :: value
end module vector

!> @brief External forces module for QM/MM calculations
!>
!> This module provides data structures and variables for handling external
!> forces in QM/MM calculations, including spline-interpolated forces and
!> van der Waals interactions.
!>
!> @details
!> The module contains:
!> - Spline parameters for force interpolation
!> - External force coefficients and knots
!> - Force arrays for different interaction types
!> - QM/MM region force contributions
!>
!> @note Used in QM/MM calculations with external force fields
!> @note Supports various force field types and interpolation schemes
module ext_forces
    integer :: f_kx, f_ky, f_kz, f_nx, f_ny, f_nz
    real*8, dimension(:, :, :), allocatable :: f_coef
    real*8, dimension(:), allocatable :: f_xknot, f_yknot, f_zknot
    real*8, dimension(:, :), allocatable :: forcext, forcext_el
    real*8, dimension(:, :), allocatable :: force_vdw
    ! QMC/MM region
    real*8, dimension(:, :), allocatable :: mm_f_theta, mm_f_dihed, mm_f_impr
end module ext_forces

!> @brief Van der Waals module for intermolecular interactions
!>
!> This module provides data structures and variables for handling van der
!> Waals interactions in QM/MM calculations, including Lennard-Jones
!> parameters and statistical averages.
!>
!> @details
!> The module contains:
!> - Lennard-Jones parameters (C12, C6)
!> - Van der Waals coordinates and connectivity
!> - Statistical averages and error estimates
!> - Support for GROMOS and CPMD force fields
!>
!> @note Used in QM/MM calculations with van der Waals interactions
!> @note Supports various force field parameter sets
module van_der_waals
    logical :: vdw
    real*8, dimension(:, :), allocatable :: c12, c6, cs12, cs6
    character(20) :: filename_vdw
    integer :: nratt, nat_nn, nat_tot
    real*8, dimension(:, :), allocatable :: coord_nn
    integer, dimension(:), allocatable :: qmc_vdw, nn_vdw
    ! Averages
    integer :: ncount_vdw = 0, t_ncount_vdw = 0
    real*8 :: ave_vdw = 0, ave2_vdw = 0
    real*8 :: err_vdw
    real*8 :: t_ave_vdw = 0, t_ave2_vdw = 0
    real*8 :: sum_pot, sum_err
    ! Gromos and CPMD
    type cp
        integer :: ind, it
    end type cp
    type(cp), dimension(:), allocatable :: cpmd
    type gromos
        integer :: ind, it
        logical :: qm
    end type gromos
    type(gromos), dimension(:), allocatable :: grom
end module van_der_waals

!> @brief Exclusion list module for molecular mechanics
!>
!> This module provides data structures for handling exclusion lists in
!> molecular mechanics calculations, including 1-2, 1-3, and 1-4 exclusions.
!>
!> @details
!> The module contains:
!> - Exclusion list data structures
!> - Maximum numbers for different exclusion types
!> - Arrays for storing excluded atom pairs
!> - Support for 1-4 interactions with scaling
!>
!> @note Used in molecular mechanics force field calculations
!> @note Supports standard exclusion rules for bonded interactions
module exc_list
    integer, parameter :: nmax = 140, nmax14 = 60
    type list
        integer :: ref, n
    end type list
    type(list), dimension(:), allocatable :: exc, exc14
    integer, dimension(:, :), allocatable :: exc_at, exc_14
end module exc_list
