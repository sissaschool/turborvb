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

function damping(r, n)
    implicit none
    real*8 damping, r
    integer n
    !     damping=1.d0
    !     return
    !     From G. Scoles et al. J. Chem. Phys. 76, 3057 (1982);
    if (n .ne. 0) then
        damping = (1.d0 - dexp(-r*2.1d0/n - 0.109d0*r**2/dsqrt(dble(n))))**n
    else
        damping = 1.d0
    end if
    return
end
