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

subroutine random(L, Lz, ztry, table, diag, iout, indvic, sign           &
        &, npow, gamma, istart, indtm, ipc)
    implicit none
    integer ipc, i, L, Lz, indvic, iout, index, ipart, istart, indtm(*)
    real*8 diag, table(ipc, *), try, ztry, sign, gamma, cost, npow

    !     choose a random configuration given the table,wtot
    !     it is assumed that  diag < ztry <= wtot ,i.e. non diagonal move

    try = diag
    i = 0
    do while (ztry .ge. try .and. i .lt. Lz)
        i = i + 1
        index = (i - 1)/L + 1
        ipart = i - (index - 1)*L
        if (npow .ne. 0.d0 .and. table(1, i) .gt. 0.d0 .and. &
            (index .ge. istart .and. index .le. indtm(ipart))) then
            try = try + table(1, i)*(1.d0 - npow*(1.d0 + gamma))
        elseif (table(1, i) .gt. 0.d0) then
            try = try + table(1, i)
        else
            cost = -gamma*table(1, i)
            try = try + cost
        end if
    end do

    !     It is possible by chance and roundoff that the final i=Lz and table(Lz)=0
    !     if for roundoff i=Lz set i the maximum with table(i)=/0
    if (i .eq. Lz) then
        do while ((table(1, i) .eq. 0.d0 .or. (table(1, i) .lt. 0 .and. gamma .eq. 0)&
                &.or. (npow*(1.d0 + gamma) .eq. 1.d0 .and. index .ge. istart)) .and. i .gt. 1)
            i = i - 1
            index = (i - 1)/L + 1
        end do
    end if

    !     if(i.eq.0) i=1
    if (table(1, i) .lt. 0.) then
        sign = -1.d0
    else
        sign = 1.d0
    end if

    indvic = index
    iout = ipart
    !      indvic=(i-1)/L+1
    !      iout=i-(indvic-1)*L

    if (table(1, i) .eq. 0.d0) indvic = -indvic !  Unrecoverable error

    return
end

