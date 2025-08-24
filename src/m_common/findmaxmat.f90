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

!> @brief Find maximum matrix indices for sparse matrix operations
!>
!> This subroutine determines the maximum row and column indices used in
!> a sparse matrix representation, which is useful for allocating memory
!> and optimizing matrix operations in quantum Monte Carlo calculations.
!>
!> @param[in] iessw Number of active orbitals to consider (0 = skip calculation)
!> @param[in] nnozero Number of non-zero elements in the sparse matrix
!> @param[in] nozero Array containing linear indices of non-zero matrix elements
!> @param[in] nelorbh Number of orbitals (matrix dimension)
!> @param[in] jbradet Array mapping sparse matrix indices to orbital indices
!> @param[in,out] lastmol Maximum row/column index found (updated on output)
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Checks if calculation should be performed (iessw > 0)
!> 2. Iterates through all non-zero matrix elements
!> 3. For each element with valid jbradet mapping (1 ≤ jbradet(j) ≤ iessw):
!>    - Converts linear index to (row, column) coordinates
!>    - Updates lastmol with the maximum of current lastmol, row, and column
!> 4. Only considers elements where both row and column are ≤ nelorbh
!>
!> @note The subroutine uses linear indexing where element (i,j) corresponds
!>       to index (j-1)*nelorbh + i
!> @note If iessw = 0, the subroutine returns immediately without calculation
!> @note The subroutine updates lastmol in-place with the maximum value found
!> @note This is typically used for memory allocation and optimization purposes
subroutine findmaxmat(iessw, nnozero, nozero, nelorbh, jbradet, lastmol)
    implicit none
    integer iessw, nelorbh, nozero(*), i, j, nnozero, jbradet(*), lastmol, ix, iy

    if (iessw .eq. 0) return

    do j = 1, nnozero
        i = jbradet(j)
        if (i .ne. 0 .and. i .le. iessw) then
            iy = (nozero(j) - 1)/nelorbh + 1
            ix = nozero(j) - (iy - 1)*nelorbh
            if (iy .le. nelorbh) lastmol = max(lastmol, ix, iy)
        end if
    end do

    return
end
