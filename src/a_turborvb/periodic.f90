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

subroutine Apply_periodic(rion, nion, cellscale, dt4, mass, sumx2)
    implicit none
    real*8 rion(3, nion), mass(3, nion), cellscale(3), dt4, sumx2&
            &, px, px2, dt4u
    integer nion, i, j
    sumx2 = 0.d0
    do i = 1, nion
        do j = 1, 3
            dt4u = dt4/mass(j, i)
            call periodic(rion(j, i), cellscale(j), dt4u, px, px2)
            rion(j, i) = px
            sumx2 = sumx2 + px2*mass(j, i)
        end do
    end do
    return
end
subroutine periodic(x, L, dt4, periodicx, periodicx2)
    implicit none
    real*8 periodicx, x, L, dt4, arg, xval, maxerror, num, den, num2&
            &, periodicx2, xmin
    real*8, external :: dlamch
    integer n
    maxerror = dlamch('E')
    !          Finds the representative in [-L/2,L/2]
    xmin = x - anint(x/L)*L
    num = xmin
    den = 1.d0
    num2 = xmin**2
    n = 0
    arg = 1.d0
    !          Positive n
    do while (arg .gt. maxerror)
        n = n + 1
        xval = xmin + L*n
        arg = dexp(-(xval**2 - xmin**2)/dt4)
        num = num + xval*arg
        num2 = num2 + xval**2*arg
        den = den + arg
    end do
    !          Negative n
    n = 0
    arg = 1.d0
    do while (arg .gt. maxerror)
        n = n - 1
        xval = xmin + L*n
        arg = dexp(-(xval**2 - xmin**2)/dt4)
        num = num + xval*arg
        num2 = num2 + xval**2*arg
        den = den + arg
    end do
    periodicx = num/den
    periodicx2 = num2/den
    return
end
