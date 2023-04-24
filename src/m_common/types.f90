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

module types
    implicit none

    type array_int
        integer, allocatable :: col(:)
    end type array_int

    type ion_comp
        integer :: mult
        integer, allocatable :: ion(:), comp(:)
    end type ion_comp

    type nkgrid
        integer :: dimshell
        integer, allocatable :: kpip(:, :)
        logical, allocatable :: tobedone(:)
    end type nkgrid

    ! jastrow/determinant object
    ! ex. for Jastrow I have:
    ! vj,vju,winvj,jasmat,jasmat_sz,jasmat_c,jasmatsz_c
    type wf_factor
        real(8), allocatable :: twobody_par(:)
        real(8), allocatable :: exps(:)
        real(8), allocatable :: bas_mat(:)
        real(8), allocatable :: exp_mat(:)
        real(8), allocatable :: exp_mat_sz(:)
        real(8), allocatable :: exp_mat_c(:)
        real(8), allocatable :: exp_mat_sz_c(:)
    end type

end module types
