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

!> @brief Complex wave function constraint builder for quantum Monte Carlo
!>
!> This subroutine constructs constraints for complex wave functions in quantum
!> Monte Carlo calculations. It handles both real and complex wave function
!> representations, applying proper sign conventions and Hermitian constraints
!> when necessary.
!>
!> @param[in] iesfreer Number of free parameters in the wave function
!> @param[in] n3body Number of three-body terms (Jastrow factors)
!> @param[in] jbraj Array mapping Jastrow terms to parameter indices
!> @param[in] derjas Array of Jastrow derivatives
!> @param[in,out] econf Energy configuration array to be updated
!> @param[in] nw Number of walkers
!> @param[in] iopt Operation mode flag
!>
!> @details
!> The subroutine performs the following operations based on the input parameters:
!>
!> **Operation Modes (iopt):**
!> - iopt = 0 or 3: Accumulate derivatives (add to existing values)
!> - iopt ≠ 0,3: Replace derivatives (set new values)
!>
!> **Wave Function Types (ipc):**
!> - ipc = 1: Real wave function
!> - ipc = 2: Complex wave function
!>
!> **Parameter Mapping (jbraj):**
!> - jbraj(i) > 0: Add derivative to parameter jbraj(i)
!> - jbraj(i) < 0: Subtract derivative from parameter |jbraj(i)|
!>
!> **Complex Wave Function Handling:**
!> For complex wave functions (ipc = 2), the subroutine handles real and
!> imaginary components separately:
!> - Real component: stored at index 2*j-1
!> - Imaginary component: stored at index 2*j
!>
!> **Hermitian Constraints:**
!> When yes_hermite is true, the subroutine applies Hermitian constraints:
!> - For negative jbraj: Real part is added, imaginary part is subtracted
!> - This ensures proper symmetry in the complex wave function
!>
!> @note The subroutine returns early if iesfreer = 0 (no free parameters).
!> @note For complex wave functions, each parameter requires two array elements
!>       (real and imaginary parts).
!> @note The Hermitian constraint ensures that the wave function satisfies
!>       proper quantum mechanical symmetries.
subroutine constrbr_complex(iesfreer, n3body, jbraj, derjas               &
        &, econf, nw, iopt)
    use Constants, only: ipc
    use allio, only: yes_hermite
    implicit none
    integer iesfreer, jbraj(*), nw, i, j, n3body, iopt
    real*8 econf(nw, *), derjas(*)

    !> @brief Early return if no free parameters
    if (iesfreer .eq. 0) return

    !> @brief Initialize energy configuration array for non-accumulation modes
    if (iopt .ne. 3) then
        do j = 1, iesfreer
            econf(1, j) = 0.d0
        end do
    end if

    !> @brief Handle real wave function case (ipc = 1)
    if (ipc .eq. 1) then

        !> @brief Accumulate derivatives (iopt = 0 or 3)
        if (iopt .eq. 0 .or. iopt .eq. 3) then

            do i = 1, n3body
                j = jbraj(i)
                !            write(*,*) ' j = ',j,derjas(i)
                if (j .gt. 0) then
                    !> @brief Add derivative to positive parameter index
                    econf(1, j) = econf(1, j) + derjas(i)
                elseif (j .lt. 0) then
                    !> @brief Subtract derivative from absolute parameter index
                    econf(1, -j) = econf(1, -j) - derjas(i)
                end if
            end do

        else
            !> @brief Replace derivatives (iopt ≠ 0,3)

            do i = 1, n3body
                j = jbraj(i)
                if (j .gt. 0) then
                    !> @brief Set derivative for positive parameter index
                    econf(1, j) = derjas(i)
                elseif (j .lt. 0) then
                    !> @brief Set negative derivative for absolute parameter index
                    econf(1, -j) = -derjas(i)
                end if
            end do

        end if

    else
        !> @brief Handle complex wave function case (ipc = 2)

        !> @brief Accumulate derivatives for complex wave function
        if (iopt .eq. 0 .or. iopt .eq. 3) then

            do i = 1, n3body
                j = jbraj(i)
                !            write(*,*) ' j = ',j,derjas(i)
                if (j .gt. 0) then
                    !> @brief Add real and imaginary derivatives to positive parameter
                    econf(1, 2*j - 1) = econf(1, 2*j - 1) + derjas(2*i - 1)
                    econf(1, 2*j) = econf(1, 2*j) + derjas(2*i)
                elseif (j .lt. 0) then
                    if (yes_hermite) then
                        !> @brief Apply Hermitian constraint: real part added, imaginary part subtracted
                        econf(1, -2*j - 1) = econf(1, -2*j - 1) + derjas(2*i - 1)
                        econf(1, -2*j) = econf(1, -2*j) - derjas(2*i)
                    else
                        !> @brief Standard complex constraint: both parts subtracted
                        econf(1, -2*j - 1) = econf(1, -2*j - 1) - derjas(2*i - 1)
                        econf(1, -2*j) = econf(1, -2*j) - derjas(2*i)
                    end if
                end if
            end do

        else
            !> @brief Replace derivatives for complex wave function

            do i = 1, n3body
                j = jbraj(i)
                if (j .gt. 0) then
                    !> @brief Set real and imaginary derivatives for positive parameter
                    econf(1, 2*j - 1) = derjas(2*i - 1)
                    econf(1, 2*j) = derjas(2*i)
                elseif (j .lt. 0) then
                    if (yes_hermite) then
                        !> @brief Apply Hermitian constraint: real part set, imaginary part negated
                        econf(1, -2*j - 1) = derjas(2*i - 1)
                        econf(1, -2*j) = -derjas(2*i)
                    else
                        !> @brief Standard complex constraint: both parts negated
                        econf(1, -2*j - 1) = -derjas(2*i - 1)
                        econf(1, -2*j) = -derjas(2*i)
                    end if
                end if
            end do

        end if
    end if ! endif ipc

    return
end
