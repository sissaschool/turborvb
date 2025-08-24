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

!> @file bconstrbr.f90
!> @brief Module containing constraint handling for real Jastrow parameters
!> @author TurboRVB group
!> @date 2022
!> @version 1.0

!> @brief Translates constrained Jastrow parameter changes to unconstrained representation
!> @details This subroutine handles the transformation of constrained Jastrow parameter 
!>          changes back to their original unconstrained representation. It is specifically 
!>          designed for handling real wave functions and 3-body Jastrow terms.
!> 
!>          The subroutine performs the following operations:
!>          1. Maps constrained Jastrow parameters to their original positions using jbraj array
!>          2. Handles real wave functions only (simplified version of bconstrbr_complex)
!>          3. Supports sign-flipped parameter mapping for constraint handling
!>          4. Specifically designed for 3-body Jastrow correlation functions
!> 
!>          The algorithm uses a mapping array jbraj to translate between constrained
!>          and unconstrained parameter spaces, where:
!>          - Positive values: direct mapping from ddw to derjas
!>          - Negative values: sign-flipped mapping from ddw to derjas
!> 
!>          This is the real-valued counterpart to bconstrbr_complex, optimized for
!>          cases where complex arithmetic is not required.
!> 
!> @param[in] iesfreer Switch flag for constraint handling (0: no constraints, >0: apply constraints)
!> @param[in] n3body Number of 3-body Jastrow parameters
!> @param[in] jbraj Mapping array from constrained to unconstrained Jastrow parameters
!> @param[out] derjas Output array for unconstrained Jastrow derivatives
!> @param[in] ddw Input array of constrained parameter changes
!> 
!> @note If iesfreer = 0, the subroutine returns immediately without any operations
!> 
!> @warning The subroutine assumes that ddw array contains valid constrained parameter values
!> 
!> @par Algorithm Details:
!> The subroutine implements simple parameter mapping for real wave functions:
!> 1. <b>Direct Mapping:</b> For positive jbraj values, copy ddw(j) to derjas(i)
!> 2. <b>Sign-Flipped Mapping:</b> For negative jbraj values, copy -ddw(-j) to derjas(i)
!> 
!> @par Parameter Mapping:
!> - <b>Positive jbraj values:</b> derjas(i) = ddw(j)
!> - <b>Negative jbraj values:</b> derjas(i) = -ddw(-j)
!> 
!> @par 3-Body Jastrow Terms:
!> This subroutine is specifically designed for 3-body Jastrow correlation functions,
!> which describe electron-electron-electron correlations in quantum Monte Carlo
!> calculations. The n3body parameter specifies the number of such terms.
!> 
!> @par Performance Considerations:
!> This real-valued version is more efficient than bconstrbr_complex for real wave
!> functions as it avoids unnecessary complex arithmetic operations.
!> 
!> @see bconstrbr_complex subroutine for complex wave function handling
!> @see bconstraint subroutine for similar constraint handling in determinant matrices
subroutine bconstrbr(iesfreer, n3body, jbraj, derjas              &
        &, ddw)
    implicit none
    
    !> @param[in] iesfreer Switch flag for constraint handling (0: no constraints, >0: apply constraints)
    integer, intent(in) :: iesfreer
    
    !> @param[in] n3body Number of 3-body Jastrow parameters
    integer, intent(in) :: n3body
    
    !> @param[in] jbraj Mapping array from constrained to unconstrained Jastrow parameters
    integer, intent(in) :: jbraj(*)
    
    !> @param[out] derjas Output array for unconstrained Jastrow derivatives
    real*8, intent(out) :: derjas(*)
    
    !> @param[in] ddw Input array of constrained parameter changes
    real*8, intent(in) :: ddw(*)
    
    ! Local variables
    integer :: i, j
    
    ! Early return if no constraints are applied
    if (iesfreer .eq. 0) return
    
    ! Map real Jastrow parameters from constrained to unconstrained space
    do i = 1, n3body
        j = jbraj(i)
        if (j .gt. 0) then
            ! Direct mapping for positive indices
            derjas(i) = ddw(j)
        elseif (j .lt. 0) then
            ! Sign-flipped mapping for negative indices
            derjas(i) = -ddw(-j)
        end if
    end do
    
    return
end
