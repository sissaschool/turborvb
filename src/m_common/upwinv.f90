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

!> @brief Update wave function inverse matrix for real wave functions
!>
!> This subroutine updates the wave function inverse matrix WINV when
!> a single electron (JEL) is moved. The update follows the Sherman-Morrison
!> formula for rank-1 updates of matrix inverses.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> JEL : integer, in
!>     Index of the electron being moved.
!> INDT : integer, in
!>     Number of basis functions (dimension of orbital space).
!> NELORB : integer, in
!>     Number of orbitals (unused in current implementation).
!> WINV : real*8 array, inout
!>     Wave function inverse matrix (NEL × INDT).
!>     On entry: current inverse matrix.
!>     On exit: updated inverse matrix.
!> V : real*8 array, in
!>     Vector of coefficients for the update (NEL).
!> PSI : real*8 array, in
!>     Wave function matrix (INDT × NEL).
!>
!> Notes
!> -----
!> - Uses OpenMP offloading for parallel computation on accelerators.
!> - The update follows the formula: WINV_new = WINV_old + v * psi^T.
!> - Special handling for the moved electron (JEL) with factor (1 + v(JEL)).
!> - Optimized for quantum Monte Carlo applications.
!>
!> Algorithm
!> ---------
!> 1. Update all electrons except JEL: WINV(j,i) += v(j) * psi(i,j)
!> 2. Update electron JEL: WINV(jel,i) = psi(i,jel) * (1 + v(jel))
!>
!> Example
!> -------
!> Used in variational Monte Carlo for efficient wave function updates.
subroutine upwinv(nel, jel, indt, nelorb, winv, v, psi)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nel, indt, jel, nelorb
    real*8, intent(in) :: psi(indt, nel), v(nel)
    real*8, intent(inout) :: winv(nel, indt)

    ! local variables
    integer i, j

#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
!> @brief Update wave function inverse for all electrons except JEL
    do j = 1, nel
        !              if(j.ne.jel) then
        do i = 1, indt
            winv(j, i) = winv(j, i) + v(j)*psi(i, j)
        end do
        !              else
        !              do i=1,indt
        !              winv(j,i)=psi(i,j)*(1.d0+v(jel))
        !              enddo
        !              endif
    end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
!> @brief Special update for the moved electron JEL
    do i = 1, indt
        winv(jel, i) = psi(i, jel)*(1.d0 + v(jel))
    end do
    return
end

!> @brief Update wave function inverse matrix for complex wave functions
!>
!> This subroutine updates the complex wave function inverse matrix WINV when
!> a single electron (JEL) is moved. The update follows the Sherman-Morrison
!> formula for rank-1 updates of complex matrix inverses.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> JEL : integer, in
!>     Index of the electron being moved.
!> INDT : integer, in
!>     Number of basis functions (dimension of orbital space).
!> NELORB : integer, in
!>     Number of orbitals (unused in current implementation).
!> WINV : complex*16 array, inout
!>     Complex wave function inverse matrix (NEL × INDT).
!>     On entry: current inverse matrix.
!>     On exit: updated inverse matrix.
!> V : complex*16 array, in
!>     Vector of complex coefficients for the update (NEL).
!> PSI : complex*16 array, in
!>     Complex wave function matrix (INDT × NEL).
!>
!> Notes
!> -----
!> - Uses OpenMP offloading for parallel computation on accelerators.
!> - The update follows the formula: WINV_new = WINV_old + v * psi^T.
!> - Special handling for the moved electron (JEL) with factor (1 + v(JEL)).
!> - Complex arithmetic is used throughout the computation.
!> - Optimized for quantum Monte Carlo applications with complex wave functions.
!>
!> Algorithm
!> ---------
!> 1. Update all electrons except JEL: WINV(j,i) += v(j) * psi(i,j)
!> 2. Update electron JEL: WINV(jel,i) = psi(i,jel) * (1 + v(jel))
!>
!> Example
!> -------
!> Used in variational Monte Carlo for efficient complex wave function updates.
subroutine upwinv_complex(nel, jel, indt, nelorb, winv, v, psi)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nel, indt, jel, nelorb
    complex*16, intent(in) :: psi(indt, nel), v(nel)
    complex*16, intent(inout) :: winv(nel, indt)

    ! local variables
    integer i, j
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
!> @brief Update complex wave function inverse for all electrons except JEL
    do j = 1, nel
        !              if(j.ne.jel) then
        do i = 1, indt
            winv(j, i) = winv(j, i) + v(j)*psi(i, j)
        end do
        !              else
        !              do i=1,indt
        !              winv(j,i)=psi(i,j)*(1.d0+v(jel))
        !              enddo
        !              endif
    end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
!> @brief Special update for the moved electron JEL
    do i = 1, indt
        winv(jel, i) = psi(i, jel)*(1.d0 + v(jel))
    end do
    return
end
