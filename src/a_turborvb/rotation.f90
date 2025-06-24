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

!> @brief      Apply 3x3 rotation matrix to a 3D vector
!> @details    This subroutine applies a 3x3 rotation matrix to a 3D vector,
!>             performing the matrix-vector multiplication rotvec = rotmatrix * vec.
!> @param[in]  rotmatrix  3x3 rotation matrix
!> @param[in]  vec        3D input vector
!> @param[out] rotvec     3D rotated output vector
subroutine rotation(rotmatrix, vec, rotvec)

    real*8 rotmatrix(3, 3), vec(3), rotvec(3)
    integer i, j

    do i = 1, 3
        rotvec(i) = 0.d0
        do j = 1, 3
            rotvec(i) = rotvec(i) + rotmatrix(i, j)*vec(j)
        end do
    end do

    return
end
