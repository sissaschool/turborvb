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

subroutine upjbra(nw, jbra, ipip, indz)
    implicit none
    integer jbra(*), ipip(*), indz(*), j, i, nw, ind
    !        before I make another permutation to optimize the reshuff
    !        in1=# walkers in each processor nw/in1=integer>=1

    do j = 1, nw
        ipip(j) = 0
    end do
    do j = 1, nw
        ind = jbra(j)
        ipip(ind) = ipip(ind) + 1
    end do
    ind = 0
    !        store the position of the killed walkers
    do j = 1, nw
        if (ipip(j) .eq. 0) then
            ind = ind + 1
            indz(ind) = j
        end if
    end do
    !        now update jbra
    ind = 0
    do j = 1, nw
        if (ipip(j) .ne. 0) then
            jbra(j) = j
            !              the replicated walkers take the position of the killed ones
            do i = 2, ipip(j)
                ind = ind + 1
                jbra(indz(ind)) = j
            end do
        end if
    end do
    return
end subroutine upjbra
