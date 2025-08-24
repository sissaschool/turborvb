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

!> @brief Vector scaling and mapping for periodic boundary conditions
!>
!> This subroutine applies periodic mapping to a set of 3D vectors using
!> the mapping function from the cell module. It is designed for efficient
!> processing of multiple vectors in quantum Monte Carlo calculations
!> with periodic boundary conditions.
!>
!> @param[in] n Number of 3D vectors to process
!> @param[in] cellfat Cell dimensions for periodic mapping (3)
!> @param[in,out] vect Array of 3D vectors to be mapped (3*n)
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Iterates through n 3D vectors stored in the vect array
!> 2. For each vector (i, i+1, i+2), applies the map function:
!>    - vect(i) = map(vect(i), cellfat(1)) for x-component
!>    - vect(i+1) = map(vect(i+1), cellfat(2)) for y-component
!>    - vect(i+2) = map(vect(i+2), cellfat(3)) for z-component
!> 3. Uses OpenMP parallelization for efficient processing of large arrays
!>
!> @note The vect array contains n consecutive 3D vectors stored as:
!>       [x1, y1, z1, x2, y2, z2, ..., xn, yn, zn]
!> @note The map function applies periodic boundary conditions to each component
!> @note OpenMP parallelization is used for performance on multi-core systems
!> @note The subroutine modifies the input vect array in-place
!> @note Used for coordinate transformations in periodic quantum Monte Carlo
subroutine scalevect(n, cellfat, vect)
    use cell, only: map
    implicit none
    integer n, i
    real*8 cellfat(3), vect(3*n)
!$omp parallel do default(shared) private(i)
    do i = 1, 3*n, 3
        vect(i) = map(vect(i), cellfat(1))
        vect(i + 1) = map(vect(i + 1), cellfat(2))
        vect(i + 2) = map(vect(i + 2), cellfat(3))
    end do
!$omp end parallel do
    return
end
