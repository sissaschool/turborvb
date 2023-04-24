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

subroutine convertdec(imax, idigit)
    implicit none
    integer imax, irest, i, idigit(6)
    !       conversion of the digit
    idigit(1) = imax/100000
    irest = imax - idigit(1)*100000
    idigit(2) = irest/10000
    irest = irest - idigit(2)*10000
    idigit(3) = irest/1000
    irest = irest - idigit(3)*1000
    idigit(4) = irest/100
    irest = irest - idigit(4)*100
    idigit(5) = irest/10
    idigit(6) = irest - 10*idigit(5)
    do i = 1, 6
        idigit(i) = idigit(i) + 48
    end do
    return
end
