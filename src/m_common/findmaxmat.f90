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

subroutine findmaxmat(iessw, nnozero, nozero, nelorbh, jbradet, lastmol)
    implicit none
    integer iessw, nelorbh, nozero(*), i, j, nnozero, jbradet(*), lastmol, ix, iy

    if (iessw .eq. 0) return

    do j = 1, nnozero
        i = jbradet(j)
        if (i .ne. 0 .and. i .le. iessw) then
            iy = (nozero(j) - 1)/nelorbh + 1
            ix = nozero(j) - (iy - 1)*nelorbh
            if (iy .le. nelorbh) lastmol = max(lastmol, ix, iy)
        end if
    end do

    return
end
