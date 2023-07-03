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

function mapping(nion, dist, zeta, b, c, imin, imax)
    !
    implicit none
    integer nion, i, imin, imax
    real*8 dist(nion), zeta(*), mapping, b, c, mindist, cost, maxdist

    imin = 1
    imax = 1
    mindist = 0.d0
    maxdist = 0.d0
    do i = 1, nion
        if (zeta(i) .ne. 0.d0) then ! the ghosts are not included in the mapping
            if (dist(i) .lt. mindist .or. mindist .eq. 0.d0) then
                imin = i
                mindist = dist(i)
            end if
            if (dist(i) .gt. maxdist .or. maxdist .eq. 0.d0) then
                imax = i
                maxdist = dist(i)
            end if
        end if
    end do

    cost = max(zeta(imin)**2, 1.d0)

    !
    mapping = (1.d0 + c*cost*mindist)/(1.d0 + b*mindist)

    !! at large distance the mapping has to be the same.
    mapping = mapping/cost

    return
end
