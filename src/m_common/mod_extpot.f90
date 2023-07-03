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

module cl_restr
    integer :: nbonds, nth, ndihe, ndimp
    real*8 :: mm_fact
    real*8 :: restr_bond, restr_angle, restr_dihe, restr_dimp
    real*8, dimension(:, :), allocatable :: restr_f_bond, restr_f_angle, restr_f_dihe, restr_f_dimp
end module cl_restr

module tot_angle
    ! Potential energy average
    integer :: ncount_bond, ncount_angle, ncount_dihed, ncount_impr
    real*8 :: ave_bond, ave2_bond, ave_angle, ave2_angle, ave_dihed, ave2_dihed, ave_impr, ave2_impr
    integer :: t_ncount_bond, t_ncount_angle, t_ncount_dihed, t_ncount_impr
    real*8 :: t_ave_bond, t_ave2_bond, t_ave_angle, t_ave2_angle, t_ave_dihed, t_ave2_dihed, t_ave_impr, t_ave2_impr
    real*8 :: err_bond, err_angle, err_dihed, err_impr
end module tot_angle

module splines
    integer :: kxord, kyord, kzord, nxknot, nyknot, nzknot
    integer :: nxcoef, nycoef, nzcoef
    real*8 :: x, y, z
    real*8 :: deltavec(3)
    real*8, dimension(:, :, :), allocatable :: bscoef
    real*8, dimension(:), allocatable :: &
            &        xknot, yknot, zknot
end module splines

module vector
    integer :: nxvec, nyvec, nzvecz
    logical :: first_call_vec
    real*8, dimension(:), allocatable :: xvec, yvec, zvec
    real*8, dimension(:, :, :), allocatable :: value
end module vector

module ext_forces
    integer :: f_kx, f_ky, f_kz, f_nx, f_ny, f_nz
    real*8, dimension(:, :, :), allocatable :: f_coef
    real*8, dimension(:), allocatable :: f_xknot, f_yknot, f_zknot
    real*8, dimension(:, :), allocatable :: forcext, forcext_el
    real*8, dimension(:, :), allocatable :: force_vdw
    ! QMC/MM region
    real*8, dimension(:, :), allocatable :: mm_f_theta, mm_f_dihed, mm_f_impr
end module ext_forces

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

module exc_list
    integer, parameter :: nmax = 140, nmax14 = 60
    type list
        integer :: ref, n
    end type list
    type(list), dimension(:), allocatable :: exc, exc14
    integer, dimension(:, :), allocatable :: exc_at, exc_14
end module exc_list
