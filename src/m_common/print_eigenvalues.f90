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

!> @brief Print eigenvalues and determine matrix rank for quantum Monte Carlo
!>
!> This subroutine prints eigenvalues to output and determines the effective
!> rank of a matrix by counting non-zero eigenvalues. It is used in quantum
!> Monte Carlo calculations for analyzing matrix properties and numerical
!> stability.
!>
!> @param[in] rank MPI rank of the calling process
!> @param[in] eig Array of eigenvalues to be printed
!> @param[in] maxdimeig Maximum dimension of eigenvalue array
!> @param[out] irankdet Effective rank of the matrix (number of non-zero eigenvalues)
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Initializes irankdet to 0
!> 2. Iterates through all eigenvalues from 1 to maxdimeig:
!>    - Checks if eigenvalue magnitude is greater than 1.0e-11
!>    - Increments irankdet for each non-zero eigenvalue
!>    - Prints eigenvalue index and value (only from rank 0)
!> 3. Returns the effective rank of the matrix
!>
!> @note Only rank 0 prints eigenvalues to avoid duplicate output in parallel
!> @note Uses a threshold of 1.0e-11 to determine non-zero eigenvalues
!> @note The effective rank represents the number of linearly independent
!>       components in the matrix
!> @note Used for analyzing overlap matrices, Hamiltonian matrices, and
!>       other matrices in quantum Monte Carlo calculations
!> @note The subroutine is useful for debugging and numerical stability analysis
subroutine print_eigenvalues(rank, eig, maxdimeig, irankdet)
    implicit none
    integer, intent(in) :: rank, maxdimeig
    integer i, irankdet
    real*8 eig(*)
    irankdet = 0
    do i = 1, maxdimeig
        if (abs(eig(i)) .gt. 1.d-11) then
            irankdet = irankdet + 1
        end if
        if (rank .eq. 0) write (6, *) i, eig(i)
    end do
    return
end subroutine print_eigenvalues
