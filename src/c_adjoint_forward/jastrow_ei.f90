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

function jastrow_ei(r, vj, iesd)
    implicit none
    real*8 r, rz, jastrow_ei, vj(*)
    integer iesd
    select case (iesd)

        ! exponential form
    case (4)
        !     r=dsqrt(rc(1)**2+rc(2)**2+rc(3)**2)
        jastrow_ei = 0.5d0/vj(1)*(1.d0 - dexp(-vj(1)*r))
        ! 1 parameter 0.5*r/(1+b*r)
    case (1)
        !     r=dsqrt(rc(1)**2+rc(2)**2+rc(3)**2)
        jastrow_ei = 0.5d0*r/(1.d0 + vj(1)*r)
        ! no jastrow
    case (8)
        ! jastrow for pseudo soft
        jastrow_ei = 1.d0/vj(1)*(1.d0 - dexp(-vj(1)*r**3))
    case (0)
        jastrow_ei = 1.d0
    end select

    return
end
