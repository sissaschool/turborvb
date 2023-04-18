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

subroutine upvpot_ei(rkel, zeta, rion, vpot, nion, LBox, epsvpot)
    implicit none
    integer nel, k, nion
    real(8) rkel(3)
    real(8) zeta(nion), vpot, ngivej, LBox, rion(3, *), cost, epsvpot
    if (epsvpot .eq. 0) then
        vpot = 1.d0
    else

        vpot = 0.d0

        do k = 1, nion
            !      Yukawa potential
            cost = ngivej(rion(1, k), rkel, LBox) + epsvpot
            vpot = vpot + 2.d0*zeta(k)/cost
        end do
        vpot = dsqrt(vpot)
    end if
    return
end
