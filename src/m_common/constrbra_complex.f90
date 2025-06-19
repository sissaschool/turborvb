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

!> @brief Complex wave function constraint builder with sparse matrix indexing
!>
!> This subroutine constructs constraints for complex wave functions using sparse
!> matrix indexing. It handles Jastrow factor derivatives and maps them to
!> wave function parameters with support for different symmetry configurations
!> and parameter types.
!>
!> @param[in] iesfreer Number of free parameters in the wave function
!> @param[in] n3body Number of three-body terms (Jastrow factors)
!> @param[in] jbraj Array mapping Jastrow terms to parameter indices
!> @param[in] nozeroj Array of sparse matrix indices for non-zero elements
!> @param[in] derjas Array of Jastrow derivatives
!> @param[in,out] econf Energy configuration array to be updated
!> @param[in] nw Number of walkers
!> @param[in] iopt Operation mode flag
!>
!> @details
!> The subroutine handles different configurations based on the following parameters:
!>
!> **Parameter Types (ipc):**
!> - ipc = 1: Real parameters (single array element per parameter)
!> - ipc = 2: Complex parameters (two array elements per parameter)
!>
!> **Symmetry Configurations (symmagp):**
!> - symmagp = .true.: AGP symmetry (4 array elements per parameter)
!> - symmagp = .false.: Standard complex (2 array elements per parameter)
!>
!> **Operation Modes (iopt):**
!> - iopt = 0: Accumulate derivatives (add to existing values)
!> - iopt ≠ 0: Replace derivatives (set new values)
!>
!> **Parameter Mapping Schemes:**
!> - Real parameters (ipc = 1): Direct mapping to single array element
!> - Complex AGP (ipc = 2, symmagp = .true.): 4-element mapping (4*j-3, 4*j-1)
!> - Complex standard (ipc = 2, symmagp = .false.): 2-element mapping (2*j-1, 2*j)
!>
!> @note The subroutine always initializes the energy configuration array to zero.
!> @note Complex parameters require special handling for real and imaginary parts.
!> @note AGP symmetry introduces additional parameter constraints.
subroutine constrbra_complex(iesfreer, n3body, jbraj, nozeroj&
        &, derjas, econf, nw, iopt)
    use Constants, only: ipc
    use allio, only: symmagp, nelorb_c
    implicit none
    integer iesfreer, jbraj(*), nozeroj(*), nw, i, j, n3body, iopt
    real*8 econf(nw, *), derjas(*)

    !> @brief Early return if no free parameters
    if (iesfreer .eq. 0) return

    !> @brief Initialize energy configuration array to zero
    do j = 1, iesfreer
        econf(1, j) = 0.d0
    end do

    !> @brief Handle complex parameters (ipc = 2)
    if (ipc .eq. 2) then
        !> @brief AGP symmetry configuration
        if (symmagp) then
            if (iopt .eq. 0) then
                !> @brief Accumulate derivatives for AGP symmetry
                do i = 1, n3body
                    j = jbraj(i)
                    if (j .gt. 0) then
                        !> @brief Add real and imaginary parts for positive parameter
                        econf(1, 4*j - 3) = econf(1, 4*j - 3) + derjas(2*nozeroj(i) - 1)
                        econf(1, 4*j - 1) = econf(1, 4*j - 1) + derjas(2*nozeroj(i))
                    elseif (j .lt. 0) then
                        !> @brief Subtract real and imaginary parts for negative parameter
                        econf(1, -4*j - 3) = econf(1, -4*j - 3) - derjas(2*nozeroj(i) - 1)
                        econf(1, -4*j - 1) = econf(1, -4*j - 1) - derjas(2*nozeroj(i))
                    end if
                end do
            else
                !> @brief Replace derivatives for AGP symmetry
                do i = 1, n3body
                    j = jbraj(i)
                    if (j .gt. 0) then
                        !> @brief Set real and imaginary parts for positive parameter
                        econf(1, 4*j - 3) = derjas(2*nozeroj(i) - 1)
                        econf(1, 4*j - 1) = derjas(2*nozeroj(i))
                    elseif (j .lt. 0) then
                        !> @brief Set negative real and imaginary parts for negative parameter
                        econf(1, -4*j - 3) = -derjas(2*nozeroj(i) - 1)
                        econf(1, -4*j - 1) = -derjas(2*nozeroj(i))
                    end if
                end do
            end if
        else ! if symmagp
            !> @brief Standard complex configuration
            if (iopt .eq. 0) then
                !> @brief Accumulate derivatives for standard complex
                do i = 1, n3body
                    j = jbraj(i)
                    if (j .gt. 0) then
                        !> @brief Add real and imaginary parts for positive parameter
                        econf(1, 2*j - 1) = econf(1, 2*j - 1) + derjas(2*nozeroj(i) - 1)
                        econf(1, 2*j) = econf(1, 2*j) + derjas(2*nozeroj(i))
                    elseif (j .lt. 0) then
                        !> @brief Subtract real and imaginary parts for negative parameter
                        econf(1, -2*j - 1) = econf(1, -2*j - 1) - derjas(2*nozeroj(i) - 1)
                        econf(1, -2*j) = econf(1, -2*j) - derjas(2*nozeroj(i))
                    end if
                end do
            else
                !> @brief Replace derivatives for standard complex
                do i = 1, n3body
                    j = jbraj(i)
                    if (j .gt. 0) then
                        !> @brief Set real and imaginary parts for positive parameter
                        econf(1, 2*j - 1) = derjas(2*nozeroj(i) - 1)
                        econf(1, 2*j) = derjas(2*nozeroj(i))
                    elseif (j .lt. 0) then
                        !> @brief Set negative real and imaginary parts for negative parameter
                        econf(1, -2*j - 1) = -derjas(2*nozeroj(i) - 1)
                        econf(1, -2*j) = -derjas(2*nozeroj(i))
                    end if
                end do
            end if
        end if

    else ! ipc=2
        !> @brief Handle real parameters (ipc = 1)

        if (iopt .eq. 0) then
            !> @brief Accumulate derivatives for real parameters

            do i = 1, n3body
                j = jbraj(i)

                if (j .gt. 0) then
                    !> @brief Add derivative for positive parameter
                    econf(1, j) = econf(1, j) + derjas(nozeroj(i))
                elseif (j .lt. 0) then
                    !> @brief Subtract derivative for negative parameter
                    econf(1, -j) = econf(1, -j) - derjas(nozeroj(i))
                end if
            end do

        else
            !> @brief Replace derivatives for real parameters

            do i = 1, n3body
                j = jbraj(i)
                if (j .gt. 0) then
                    !> @brief Set derivative for positive parameter
                    econf(1, j) = derjas(nozeroj(i))
                elseif (j .lt. 0) then
                    !> @brief Set negative derivative for negative parameter
                    econf(1, -j) = -derjas(nozeroj(i))
                end if
            end do
        end if

    end if
    return
end
