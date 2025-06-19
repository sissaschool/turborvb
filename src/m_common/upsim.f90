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

!> @brief Matrix update utilities for quantum Monte Carlo calculations
!>
!> This module provides subroutines for updating matrix elements in quantum
!> Monte Carlo calculations, supporting both real and complex matrices.
!> It handles different matrix types including standard matrices and
!> Pfaffian matrices for paired and unpaired electrons.
!>
!> The module includes:
!> - Real matrix update with symmetry handling
!> - Complex matrix update with Hermitian/transpose symmetry
!> - Support for Pfaffian matrices with up/down electron blocks
!> - Proper handling of paired and unpaired electron configurations
!>
!> @author TurboRVB group
!> @version 1.0
!> @date 2022

!> @brief Update real matrix elements with symmetry handling
!>
!> This subroutine updates matrix elements in a real matrix, handling
!> different matrix types and symmetry conditions. It supports both
!> standard matrices and Pfaffian matrices for quantum Monte Carlo
!> calculations.
!>
!> @param[in,out] amat Matrix to be updated (ndim, *)
!> @param[in] ndim Leading dimension of the matrix
!> @param[in] ind Linear index of the matrix element to update
!> @param[in] value New value for the matrix element
!> @param[in] symmagp Flag for symmetric matrix updates
!> @param[in] ipf Flag indicating matrix type:
!>                - ipf = 1: Standard matrix
!>                - ipf ≠ 1: Pfaffian matrix
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Converts linear index to (i,j) matrix coordinates
!> 2. For standard matrices (ipf = 1):
!>    - Sets amat(i,j) = value
!>    - If symmagp is true, sets amat(j,i) = value
!> 3. For Pfaffian matrices (ipf ≠ 1):
!>    - Sets amat(i,j) = value and amat(j,i) = -value
!>    - Handles up/down electron blocks with proper symmetry
!>    - Manages paired and unpaired electron configurations
!>    - Applies symmetry conditions based on ion types
!>
!> @note The subroutine handles antisymmetry of Pfaffian matrices
!> @note Up/down electron blocks are handled separately for paired electrons
!> @note The subroutine respects ion type constraints (kiontot)
!> @note Used in wave function optimization and matrix updates
subroutine upsim(amat, ndim, ind, value, symmagp, ipf)
    use allio, only: nelorb_at, pfaffup, kiontot
    implicit none

    integer ndim, ndimh, ind, i, j, ipf
    real(8) amat(ndim, *), value
    logical symmagp
    if (ipf .eq. 1) then
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value
        if (symmagp .and. j .le. ndim) amat(j, i) = value
    else
        ! if symmagp
        !   Pfaff(up,down)=A with A= A^T  and Pfaff(up,up)=Pfaff(do,do)
        !  In the complex case two cases are possible
        !  if yeshermite      --> A=A^+  and Pfaff(up,up)=Pfaff(do,do)^*
        !  if .not.yeshermite --> A=A^T  and Pfaff(up,up)=Pfaff(do,do)
        !  On top of that
        !  if pfaffup=.true. the down-down part of the pfaffian is absent

        ndimh = nelorb_at/2
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value
        if (j .le. ndim) amat(j, i) = -value ! With the unpaired j>ndim
        if (symmagp .and. kiontot(i) .ne. 0 .and. kiontot(j) .ne. 0) then
            if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                amat(i + ndimh, j + ndimh) = value
                amat(j + ndimh, i + ndimh) = -value
            elseif (j .gt. ndimh .and. j .le. nelorb_at .and. i .le. ndimh) then
                amat(j - ndimh, i + ndimh) = value
                amat(i + ndimh, j - ndimh) = -value
                !       elseif(j.le.ndimh.and.i.gt.ndimh) then
                !          amat(j+ndimh, i-ndimh)=value
                !          amat(i-ndimh, j+ndimh)=-value
            end if
        end if
    end if
    return
end subroutine upsim

!> @brief Update complex matrix elements with symmetry handling
!>
!> This subroutine updates matrix elements in a complex matrix, handling
!> different matrix types and symmetry conditions. It supports both
!> standard matrices and Pfaffian matrices with proper handling of
!> Hermitian and transpose symmetries.
!>
!> @param[in,out] amat Complex matrix to be updated (ndim, *)
!> @param[in] ndim Leading dimension of the matrix
!> @param[in] ind Linear index of the matrix element to update
!> @param[in] value Array containing new complex value for the matrix element
!> @param[in] symmagp Flag for symmetric matrix updates
!> @param[in] ipf Flag indicating matrix type:
!>                - ipf = 1: Standard matrix
!>                - ipf ≠ 1: Pfaffian matrix
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Converts linear index to (i,j) matrix coordinates
!> 2. For standard matrices (ipf = 1):
!>    - Sets amat(i,j) = value(1)
!>    - If symmagp is true, applies Hermitian or transpose symmetry
!> 3. For Pfaffian matrices (ipf ≠ 1):
!>    - Sets amat(i,j) = value(1) and amat(j,i) = -value(1)
!>    - Handles up/down electron blocks with proper symmetry
!>    - Applies Hermitian or transpose symmetry based on yes_hermite flag
!>    - Manages paired and unpaired electron configurations
!>
!> @note For Hermitian matrices (yes_hermite = .true.), uses conjugate symmetry
!> @note For non-Hermitian matrices, uses transpose symmetry
!> @note The subroutine handles antisymmetry of Pfaffian matrices
!> @note Up/down electron blocks are handled separately for paired electrons
!> @note Used in complex wave function optimization and matrix updates
subroutine upsim_complex(amat, ndim, ind, value, symmagp, ipf)
    use allio, only: yes_hermite, nelorb_at, pfaffup, kiontot
    implicit none
    integer ndim, ndimh, ind, i, j, ipf
    complex*16 amat(ndim, *), value(*)
    logical symmagp
    !      write(6,*) value
    if (ipf .eq. 1) then
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value(1)
        if (symmagp .and. j .le. nelorb_at .and. i .le. nelorb_at) then
            if (yes_hermite) then
                amat(j, i) = conjg(value(1))
            else
                amat(j, i) = value(1)
            end if
        end if
    else
        ndimh = nelorb_at/2
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value(1)
        if (j .le. ndim) amat(j, i) = -value(1) !   if it is unpaired element j,i is not present

        if (symmagp .and. kiontot(i) .ne. 0 .and. kiontot(j) .ne. 0) then
            if (yes_hermite) then
                if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                    amat(i + ndimh, j + ndimh) = conjg(value(1))
                    amat(j + ndimh, i + ndimh) = -conjg(value(1))
                elseif (j .gt. ndimh .and. i .le. ndimh) then
                    amat(j - ndimh, i + ndimh) = conjg(value(1))
                    amat(i + ndimh, j - ndimh) = -conjg(value(1))
                    !          elseif(j.le.ndimh.and.i.gt.ndimh) then
                    !             amat(j+ndimh, i-ndimh)=conjg(value(1))
                    !             amat(i-ndimh, j+ndimh)=-conjg(value(1))
                end if
            else
                if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                    amat(i + ndimh, j + ndimh) = value(1)
                    amat(j + ndimh, i + ndimh) = -value(1)
                elseif (j .gt. ndimh .and. i .le. ndimh) then
                    amat(j - ndimh, i + ndimh) = value(1)
                    amat(i + ndimh, j - ndimh) = -value(1)
                    !          elseif(j.le.ndimh.and.i.gt.ndimh) then
                    !             amat(j+ndimh, i-ndimh)=value(1)
                    !             amat(i-ndimh, j+ndimh)=-value(1)
                end if
            end if
        end if
    end if
    return
end subroutine upsim_complex
