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

character(6) function intc(i)
    !
    character(6) temp
    integer, intent(in) :: i
    !
    if (i .lt. 10 .and. i .ge. 0) then
        write (temp, '(i1)') i
    else if (i .lt. 100 .and. i .gt. -10) then
        write (temp, '(i2)') i
    else if (i .lt. 1000 .and. i .gt. -100) then
        write (temp, '(i3)') i
    else if (i .lt. 10000 .and. i .gt. -1000) then
        write (temp, '(i4)') i
    else if (i .lt. 100000 .and. i .gt. -10000) then
        write (temp, '(i5)') i
    else if (i .lt. 1000000 .and. i .gt. -100000) then
        write (temp, '(i6)') i
    else
        write (temp, '(a6)') "******"
    end if
    intc = temp
    !
end function intc
