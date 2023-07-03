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

subroutine findzero(eig_min, eig_max, npm, ndim, psi, eig         &
        &, emin, emax, green, e)
    integer npm, ndim, k, i, maxit
    real*8 psi(npm, *), eig(*), green(*)                            &
            &, eig_min, eig_max, emin(*), emax(*), e(*)
    maxit = 50
    call dscalzero(ndim, eig_min, emin, 1)
    call dscalzero(ndim, eig_max, emax, 1)
    do i = 1, maxit
        do k = 1, ndim
            e(k) = (emin(k) + emax(k))*0.5d0
        end do
        call evalgreen(npm, ndim, psi, eig, e, green)
        do k = 1, ndim
            if (green(k) .gt. 0.d0) then
                emax(k) = e(k)
            else
                emin(k) = e(k)
            end if
        end do
    end do

    return
end

subroutine evalgreen(npm, ndim, psi, eig, e, green)
    implicit none
    integer npm, ndim, k_no, i, j, k, l
    real*8 psi(npm, *), eig(*), green(*), e(*)
    do k = 1, ndim
        green(k) = psi(k, 1)**2/(eig(1) - e(k))
    end do
    do l = 2, ndim
        do k = 1, ndim
            green(k) = green(k) + psi(k, l)**2/(eig(l) - e(k))
        end do
    end do
    return
end
