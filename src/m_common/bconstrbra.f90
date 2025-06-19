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

!> @file bconstrbra.f90
!> @brief Module containing constraint handling for Jastrow matrix operations
!> @author TurboRVB group
!> @date 2022
!> @version 1.0

!> @brief Translates constrained Jastrow parameter changes to matrix representation
!> @details This subroutine handles the transformation of constrained Jastrow parameter 
!>          changes back to their original unconstrained representation in matrix form.
!>          It uses the upsim subroutine to update the Jastrow matrix elements directly.
!> 
!>          The subroutine performs the following operations:
!>          1. Maps constrained Jastrow parameters to their original positions using jbraj array
!>          2. Updates Jastrow matrix elements using upsim subroutine
!>          3. Supports sign-flipped parameter mapping for constraint handling
!>          4. Specifically designed for 3-body Jastrow correlation functions
!>          5. Uses matrix update operations instead of direct array assignment
!> 
!>          The algorithm uses a mapping array jbraj to translate between constrained
!>          and unconstrained parameter spaces, where:
!>          - Positive values: direct mapping from ddw to matrix elements
!>          - Negative values: sign-flipped mapping from ddw to matrix elements
!> 
!>          This version uses matrix update operations (upsim) which may include
!>          additional matrix structure maintenance compared to direct assignment.
!> 
!> @param[in] iesfreer Switch flag for constraint handling (0: no constraints, >0: apply constraints)
!> @param[in] n3body Number of 3-body Jastrow parameters
!> @param[in] jbraj Mapping array from constrained to unconstrained Jastrow parameters
!> @param[in] nozeroj Array of non-zero element indices for Jastrow matrix
!> @param[in,out] derjas Jastrow matrix to be updated
!> @param[in] nelorbj Number of Jastrow orbitals
!> @param[in] ddw Input array of constrained parameter changes
!> 
!> @note If iesfreer = 0, the subroutine returns immediately without any operations
!> 
!> @warning The subroutine modifies the derjas matrix in-place using upsim operations
!> 
!> @par Algorithm Details:
!> The subroutine implements matrix-based parameter mapping:
!> 1. <b>Direct Mapping:</b> For positive jbraj values, call upsim with ddw(j)
!> 2. <b>Sign-Flipped Mapping:</b> For negative jbraj values, call upsim with -ddw(-j)
!> 
!> @par Matrix Update Operations:
!> Uses upsim subroutine with parameters:
!> - derjas: Target Jastrow matrix
!> - nelorbj: Number of Jastrow orbitals
!> - nozeroj(i): Non-zero element index
!> - ddw(j) or -ddw(-j): Parameter value
!> - .true.: Symmetry flag
!> - 1: Matrix type flag
!> 
!> @par 3-Body Jastrow Terms:
!> This subroutine is specifically designed for 3-body Jastrow correlation functions,
!> which describe electron-electron-electron correlations in quantum Monte Carlo
!> calculations. The n3body parameter specifies the number of such terms.
!> 
!> @see bconstrbr subroutine for direct array assignment version
!> @see bconstrbra_sparse subroutine for sparse matrix version
!> @see upsim subroutine for matrix update operations
subroutine bconstrbra(iesfreer, n3body, jbraj, nozeroj, derjas     &
        &, nelorbj, ddw)
    implicit none
    
    !> @param[in] iesfreer Switch flag for constraint handling (0: no constraints, >0: apply constraints)
    integer, intent(in) :: iesfreer
    
    !> @param[in] n3body Number of 3-body Jastrow parameters
    integer, intent(in) :: n3body
    
    !> @param[in] jbraj Mapping array from constrained to unconstrained Jastrow parameters
    integer, intent(in) :: jbraj(*)
    
    !> @param[in] nozeroj Array of non-zero element indices for Jastrow matrix
    integer, intent(in) :: nozeroj(*)
    
    !> @param[in] nelorbj Number of Jastrow orbitals
    integer, intent(in) :: nelorbj
    
    !> @param[in] ddw Input array of constrained parameter changes
    real*8, intent(in) :: ddw(*)
    
    !> @param[in,out] derjas Jastrow matrix to be updated
    real*8, intent(inout) :: derjas(*)
    
    ! Local variables
    integer :: i, j
    
    ! Early return if no constraints are applied
    if (iesfreer .eq. 0) return
    
    ! Map constrained Jastrow parameters to matrix elements using upsim
    do i = 1, n3body
        j = jbraj(i)
        if (j .gt. 0) then
            ! Direct mapping for positive indices using matrix update
            call upsim(derjas, nelorbj, nozeroj(i), ddw(j), .true., 1)
        elseif (j .lt. 0) then
            ! Sign-flipped mapping for negative indices using matrix update
            call upsim(derjas, nelorbj, nozeroj(i), -ddw(-j), .true., 1)
        end if
    end do
    
    return
end

!> @brief Translates constrained Jastrow parameter changes using sparse matrix operations
!> @details This subroutine handles the transformation of constrained Jastrow parameter 
!>          changes back to their original unconstrained representation using sparse
!>          matrix operations. It is an optimized version that avoids matrix update
!>          operations for better performance with sparse matrices.
!> 
!>          The subroutine performs the following operations:
!>          1. Maps constrained Jastrow parameters using sparse indexing (nozeroj)
!>          2. Direct assignment to sparse matrix elements
!>          3. Supports sign-flipped parameter mapping for constraint handling
!>          4. Specifically designed for 3-body Jastrow correlation functions
!>          5. Optimized for sparse matrix operations
!> 
!>          The algorithm uses a mapping array jbraj with sparse indexing through nozeroj
!>          to translate between constrained and unconstrained parameter spaces, where:
!>          - Positive values: direct mapping from ddw to sparse matrix elements
!>          - Negative values: sign-flipped mapping from ddw to sparse matrix elements
!> 
!>          This version is optimized for sparse matrices by avoiding matrix update
!>          operations and using direct assignment instead.
!> 
!> @param[in] iesfreer Switch flag for constraint handling (0: no constraints, >0: apply constraints)
!> @param[in] n3body Number of 3-body Jastrow parameters
!> @param[in] jbraj Mapping array from constrained to unconstrained Jastrow parameters
!> @param[in] nozeroj Array of non-zero element indices for sparse Jastrow matrix
!> @param[out] derjas Sparse Jastrow matrix to be updated
!> @param[in] nelorbj Number of Jastrow orbitals
!> @param[in] ddw Input array of constrained parameter changes
!> 
!> @note If iesfreer = 0, the subroutine returns immediately without any operations
!> 
!> @warning The subroutine uses direct assignment instead of matrix update operations
!> 
!> @par Algorithm Details:
!> The subroutine implements sparse matrix parameter mapping:
!> 1. <b>Direct Mapping:</b> For positive jbraj values, direct assignment derjas(i) = ddw(j)
!> 2. <b>Sign-Flipped Mapping:</b> For negative jbraj values, direct assignment derjas(i) = -ddw(-j)
!> 
!> @par Sparse Matrix Optimization:
!> Uses direct assignment operations instead of upsim calls:
!> - derjas(i) = ddw(j) for positive jbraj values
!> - derjas(i) = -ddw(-j) for negative jbraj values
!> - Avoids matrix structure maintenance overhead
!> 
!> @par Performance Benefits:
!> - Faster execution for sparse matrices
!> - Reduced memory access patterns
!> - Avoids unnecessary matrix update operations
!> 
!> @par 3-Body Jastrow Terms:
!> This subroutine is specifically designed for 3-body Jastrow correlation functions,
!> which describe electron-electron-electron correlations in quantum Monte Carlo
!> calculations. The n3body parameter specifies the number of such terms.
!> 
!> @see bconstrbra subroutine for matrix update version
!> @see bconstrbr subroutine for direct array assignment version
subroutine bconstrbra_sparse(iesfreer, n3body, jbraj, nozeroj, derjas&
        &, nelorbj, ddw)
    implicit none
    
    !> @param[in] iesfreer Switch flag for constraint handling (0: no constraints, >0: apply constraints)
    integer, intent(in) :: iesfreer
    
    !> @param[in] n3body Number of 3-body Jastrow parameters
    integer, intent(in) :: n3body
    
    !> @param[in] jbraj Mapping array from constrained to unconstrained Jastrow parameters
    integer, intent(in) :: jbraj(*)
    
    !> @param[in] nozeroj Array of non-zero element indices for sparse Jastrow matrix
    integer, intent(in) :: nozeroj(*)
    
    !> @param[in] nelorbj Number of Jastrow orbitals
    integer, intent(in) :: nelorbj
    
    !> @param[in] ddw Input array of constrained parameter changes
    real*8, intent(in) :: ddw(*)
    
    !> @param[out] derjas Sparse Jastrow matrix to be updated
    real*8, intent(out) :: derjas(*)
    
    ! Local variables
    integer :: i, j
    
    ! Early return if no constraints are applied
    if (iesfreer .eq. 0) return
    
    ! Map constrained Jastrow parameters to sparse matrix elements using direct assignment
    do i = 1, n3body
        j = jbraj(nozeroj(i))
        if (j .gt. 0) then
            ! Direct mapping for positive indices using sparse indexing
            derjas(i) = ddw(j)
!           call  upsim(derjas, nelorbj, nozeroj(i), ddw(j), .true., 1)
        elseif (j .lt. 0) then
            ! Sign-flipped mapping for negative indices using sparse indexing
            derjas(i) = -ddw(-j)
!           call  upsim(derjas, nelorbj, nozeroj(i), -ddw(-j), .true., 1)
        end if
    end do
    
    return
end
