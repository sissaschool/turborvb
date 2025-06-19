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

!> @brief Real wave function constraint builder with sparse matrix indexing
!>
!> This subroutine constructs constraints for real wave functions using sparse
!> matrix indexing. It handles Jastrow factor derivatives and maps them to
!> wave function parameters using a sparse matrix representation.
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
!> The subroutine performs the following operations:
!>
!> **Operation Modes (iopt):**
!> - iopt = 0: Accumulate derivatives (add to existing values)
!> - iopt ≠ 0: Replace derivatives (set new values)
!>
!> **Sparse Matrix Handling:**
!> The nozeroj array provides indices into the sparse derivative array,
!> allowing efficient handling of only non-zero matrix elements.
!>
!> **Parameter Mapping:**
!> - jbraj(i) > 0: Add derivative to parameter jbraj(i)
!> - jbraj(i) < 0: Subtract derivative from parameter |jbraj(i)|
!>
!> @note This subroutine always initializes the energy configuration array to zero.
!> @note The sparse indexing allows efficient handling of large, sparse matrices.
subroutine constrbra(iesfreer, n3body, jbraj, nozeroj, derjas      &
        &, econf, nw, iopt)
    implicit none
    integer iesfreer, jbraj(*), nozeroj(*), nw, i, j, n3body, iopt
    real*8 econf(nw, *), derjas(*)

    !> @brief Early return if no free parameters
    if (iesfreer .eq. 0) return

    !> @brief Initialize energy configuration array to zero
    do j = 1, iesfreer
        econf(1, j) = 0.d0
    end do

    !> @brief Accumulate derivatives (iopt = 0)
    if (iopt .eq. 0) then

        do i = 1, n3body
            j = jbraj(i)

            if (j .gt. 0) then
                !> @brief Add derivative using sparse matrix indexing
                econf(1, j) = econf(1, j) + derjas(nozeroj(i))
            elseif (j .lt. 0) then
                !> @brief Subtract derivative using sparse matrix indexing
                econf(1, -j) = econf(1, -j) - derjas(nozeroj(i))
            end if
        end do

    else
        !> @brief Replace derivatives (iopt ≠ 0)

        do i = 1, n3body
            j = jbraj(i)
            if (j .gt. 0) then
                !> @brief Set derivative using sparse matrix indexing
                econf(1, j) = derjas(nozeroj(i))
            elseif (j .lt. 0) then
                !> @brief Set negative derivative using sparse matrix indexing
                econf(1, -j) = -derjas(nozeroj(i))
            end if
        end do

    end if

    return
end

!> @brief Sparse matrix constraint builder with alternative indexing
!>
!> This subroutine constructs constraints using an alternative sparse matrix
!> indexing scheme. It differs from constrbra in the order of array indexing
!> for the sparse matrix elements.
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
!> The subroutine performs the following operations:
!>
!> **Operation Modes (iopt):**
!> - iopt = 0: Accumulate derivatives (add to existing values)
!> - iopt ≠ 0: Replace derivatives (set new values)
!>
!> **Alternative Sparse Indexing:**
!> This version uses jbraj(nozeroj(i)) instead of jbraj(i), providing
!> a different mapping between sparse matrix elements and parameter indices.
!>
!> **Parameter Mapping:**
!> - jbraj(nozeroj(i)) > 0: Add derivative to parameter jbraj(nozeroj(i))
!> - jbraj(nozeroj(i)) < 0: Subtract derivative from parameter |jbraj(nozeroj(i))|
!>
!> @note This subroutine always initializes the energy configuration array to zero.
!> @note The alternative indexing scheme may be more efficient for certain
!>       sparse matrix structures.
subroutine constrbra_sparse(iesfreer, n3body, jbraj, nozeroj, derjas&
        &, econf, nw, iopt)
    implicit none
    integer iesfreer, jbraj(*), nozeroj(*), nw, i, j, n3body, iopt
    real*8 econf(nw, *), derjas(*)

    !> @brief Early return if no free parameters
    if (iesfreer .eq. 0) return

    !> @brief Initialize energy configuration array to zero
    do j = 1, iesfreer
        econf(1, j) = 0.d0
    end do

    !> @brief Accumulate derivatives (iopt = 0)
    if (iopt .eq. 0) then

        do i = 1, n3body
            j = jbraj(nozeroj(i))

            if (j .gt. 0) then
                !> @brief Add derivative using alternative sparse indexing
                econf(1, j) = econf(1, j) + derjas(i)
            elseif (j .lt. 0) then
                !> @brief Subtract derivative using alternative sparse indexing
                econf(1, -j) = econf(1, -j) - derjas(i)
            end if
        end do

    else
        !> @brief Replace derivatives (iopt ≠ 0)

        do i = 1, n3body
            j = jbraj(nozeroj(i))
            if (j .gt. 0) then
                !> @brief Set derivative using alternative sparse indexing
                econf(1, j) = derjas(i)
            elseif (j .lt. 0) then
                !> @brief Set negative derivative using alternative sparse indexing
                econf(1, -j) = -derjas(i)
            end if
        end do

    end if

    return
end
