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

function enercutoff(etry, srpar, costcut, ener)
    implicit none
    real*8 etry, costcut, ener, srpar, enercutoff
    if (costcut .ne. 0.d0) then
        enercutoff = costcut*tanh((ener - etry)/costcut) + etry
    else
        if (1.d0 + srpar*ener .ge. 0) then
            enercutoff = 0.d0
        else
            enercutoff = -2.d0/srpar
        end if
    end if
    return
end
