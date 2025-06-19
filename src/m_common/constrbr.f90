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

!> @brief Real wave function constraint builder for quantum Monte Carlo
!>
!> This subroutine constructs constraints for real wave functions in quantum
!> Monte Carlo calculations. It handles Jastrow factor derivatives and maps
!> them to the appropriate wave function parameters using sign conventions.
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
!> **Parameter Mapping (jbraj):**
!> - jbraj(i) > 0: Add derivative to parameter jbraj(i)
!> - jbraj(i) < 0: Subtract derivative from parameter |jbraj(i)|
!>
!> **Algorithm:**
!> 1. Early return if no free parameters (iesfreer = 0)
!> 2. Initialize energy configuration array to zero for non-accumulation modes
!> 3. Loop through all three-body terms (Jastrow factors)
!> 4. Apply parameter mapping with appropriate sign conventions
!> 5. Update energy configuration array accordingly
!>
!> @note This subroutine is the real wave function version of constrbr_complex.
!> @note The subroutine returns early if iesfreer = 0 (no free parameters).
!> @note For real wave functions, each parameter requires only one array element.
subroutine constrbr(iesfreer, n3body, jbraj, derjas               &
        &, econf, nw, iopt)
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

    return
end
