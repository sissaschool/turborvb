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

subroutine make_u(alpha, x, y, z, u)
    implicit none
    real*8 alpha, x, y, z, uc
    real*8 u(3, 3)

    uc = 1.d0 - dcos(alpha)

    u(1, 1) = dcos(alpha) + uc*x**2
    u(1, 2) = uc*x*y - dsin(alpha)*z
    u(1, 3) = uc*x*z + dsin(alpha)*y

    u(2, 1) = uc*y*x + dsin(alpha)*z
    u(2, 2) = dcos(alpha) + uc*y**2
    u(2, 3) = uc*y*z - dsin(alpha)*x

    u(3, 1) = uc*z*x - dsin(alpha)*y
    u(3, 2) = uc*z*y + dsin(alpha)*x
    u(3, 3) = dcos(alpha) + uc*z**2

    return
end subroutine make_u
