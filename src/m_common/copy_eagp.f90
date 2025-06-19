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

!> @brief Copy and transform AGP (Antisymmetrized Geminal Power) matrices
!>
!> This subroutine handles copying and transformation between different
!> representations of AGP matrices. It can convert between a small matrix
!> representation and a full matrix representation, including Pfaffian
!> matrix elements.
!>
!> @param[in] small2big Direction flag for matrix transformation
!> @param[in] ipc Parameter type flag (1=real, 2=complex)
!> @param[in] nelorb_c Number of orbitals in the core
!> @param[in] nelcol_c Number of columns in the full matrix
!> @param[in,out] detmat_small Small matrix representation (nelorb_c × nelcol_c)
!> @param[in,out] eagp_pfaff Pfaffian matrix elements ((nelcol_c-nelorb_c) × (nelcol_c-nelorb_c))
!> @param[in,out] detmat_big Full matrix representation (nelcol_c × nelcol_c)
!>
!> @details
!> The subroutine performs different operations based on the small2big flag:
!>
!> **small2big = .true. (Small to Big transformation):**
!> 1. Copy the small matrix to the upper-left block of the big matrix
!> 2. Apply antisymmetry: detmat_big(i,j) = -detmat_small(j,i) for i > nelorb_c, j ≤ nelorb_c
!> 3. Copy Pfaffian elements to the lower-right block of the big matrix
!>
!> **small2big = .false. (Big to Small transformation):**
!> 1. Copy the upper-left block of the big matrix to the small matrix
!> 2. Extract Pfaffian elements from the lower-right block of the big matrix
!>
!> **Matrix Structure:**
!> The full matrix (detmat_big) has the structure:
!> ```
!> [ detmat_small        ]
!> [ -detmat_small^T     eagp_pfaff ]
!> ```
!>
!> **Parameter Handling:**
!> - ipc = 1: Real matrices (single precision elements)
!> - ipc = 2: Complex matrices (double precision elements, real/imaginary parts)
!>
!> @note The subroutine assumes proper matrix dimensions and memory allocation.
!> @note For complex matrices (ipc=2), the transformation handles real and imaginary parts.
!> @note The antisymmetry property is enforced during the transformation.
!> @note This subroutine is essential for AGP wave function calculations.
subroutine copy_eagp(small2big, ipc, nelorb_c, nelcol_c, detmat_small, eagp_pfaff, detmat_big)
    implicit none
    logical small2big
    integer ipc, nelorb_c, nelcol_c, i, j
    real*8 detmat_small(ipc*nelorb_c, nelcol_c), detmat_big(ipc*nelcol_c, nelcol_c)&
            &, eagp_pfaff(ipc*(nelcol_c - nelorb_c), nelcol_c - nelorb_c)
    
    !> @brief Transform from small to big matrix representation
    if (small2big) then
        !> @brief Copy small matrix to upper-left block of big matrix
        do j = 1, nelcol_c
            do i = 1, ipc*nelorb_c
                detmat_big(i, j) = detmat_small(i, j)
            end do
        end do
        
        !> @brief Apply antisymmetry and copy Pfaffian elements
        do i = nelorb_c + 1, nelcol_c
            !> @brief Apply antisymmetry for upper-right block
            do j = 1, nelorb_c
                detmat_big(ipc*(i - 1) + 1:ipc*i, j) = -detmat_small(ipc*(j - 1) + 1:ipc*j, i)
            end do
            !> @brief Copy Pfaffian elements to lower-right block
            do j = nelorb_c + 1, nelcol_c
                detmat_big(ipc*(i - 1) + 1:ipc*i, j) = &
                        &eagp_pfaff(ipc*(i - 1 - nelorb_c) + 1:ipc*(i - nelorb_c), j - nelorb_c)
            end do
        end do
    else
        !> @brief Transform from big to small matrix representation
        !> @brief Copy upper-left block of big matrix to small matrix
        do j = 1, nelcol_c
            do i = 1, ipc*nelorb_c
                detmat_small(i, j) = detmat_big(i, j)
            end do
        end do
        
        !> @brief Extract Pfaffian elements from lower-right block
        do i = nelorb_c + 1, nelcol_c
            do j = nelorb_c + 1, nelcol_c
                eagp_pfaff(ipc*(i - nelorb_c - 1) + 1:ipc*(i - nelorb_c), j - nelorb_c) = &
                        &detmat_big(ipc*(i - 1) + 1:ipc*i, j)
            end do
        end do
    end if
end subroutine copy_eagp
