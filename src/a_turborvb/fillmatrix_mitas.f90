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

subroutine fillmatrix(rotmatrix)
    implicit none

    real*8 x1, x2, x3, xsum, xsum2, theta, sthet                            &
            &, cfi, sfi, yy1, yy2, yy3, zz1, zz2, zz3, u1, u2, usum, uu, y1, y2, y3, z1, z2, z3  &
            &, drand1, rotmatrix(3, 3)

2   x1 = 1.d0 - 2.d0*drand1()
    x2 = 1.d0 - 2.d0*drand1()
    xsum = x1*x1 + x2*x2
    if (xsum .ge. 1.d0) goto 2
    xsum2 = 2.d0*dsqrt(dabs(1.d0 - xsum))
    x1 = x1*xsum2
    x2 = x2*xsum2
    x3 = 1.d0 - 2.d0*xsum

    theta = dacos(x3)
    sthet = dsin(theta)
    if (sthet .lt. 1.d-07) then
        cfi = 1.d0
        sfi = 0.d0
        x1 = 0.d0
        x2 = 0.d0
        x3 = 1.d0
        sthet = 0.d0
    else
        cfi = x1/sthet
        sfi = x2/sthet
    end if

    yy1 = x3*cfi
    yy2 = x3*sfi
    yy3 = -sthet
    zz1 = -sfi
    zz2 = cfi
    zz3 = 0.d0

3   u1 = drand1()*2.d0 - 1.d0
    u2 = drand1()*2.d0 - 1.d0
    usum = u1*u1 + u2*u2
    if (usum .ge. 1.d0 .or. usum .le. 0.d0) goto 3
    uu = dsqrt(usum)

    u1 = u1/uu
    u2 = u2/uu
    y1 = yy1*u1 + zz1*u2
    y2 = yy2*u1 + zz2*u2
    y3 = yy3*u1 + zz3*u2
    z1 = yy1*u2 - zz1*u1
    z2 = yy2*u2 - zz2*u1
    z3 = yy3*u2 - zz3*u1

    rotmatrix(1, 1) = x1
    rotmatrix(2, 1) = x2
    rotmatrix(3, 1) = x3
    rotmatrix(1, 2) = y1
    rotmatrix(2, 2) = y2
    rotmatrix(3, 2) = y3
    rotmatrix(1, 3) = z1
    rotmatrix(2, 3) = z2
    rotmatrix(3, 3) = z3

    return
end
