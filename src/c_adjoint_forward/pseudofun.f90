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

function pseudofun(nmax, r, param, psip)
    !
    implicit none
    real*8 pseudofun, param(3, *), r, logr, r2, psip(*)
    integer nmax, i
    !
    if (r .lt. 1d-9) r = 1d-9
    r2 = r**2
    logr = dlog(r)

    do i = 1, nmax
        psip(i) = dexp(-param(3, i)*r2 + logr*param(2, i))
    end do
    !        call my_dexp(nmax,psip,psip)
    !

    !        do i=1,nmax
    !        psip(i)=param(1,i)*psip(i)
    !        enddo
    pseudofun = 0.d0
    do i = 1, nmax
        !        pseudofun=pseudofun+param(1,i)*r**param(2,i)
        !    &               *dexp(-param(3,i)*r2)
        pseudofun = pseudofun + psip(i)*param(1, i)
    end do
    pseudofun = pseudofun/r2
    !
    return
end
