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

subroutine PUSHINTEGER4(adi4ibuf, adi4buf, x)
    integer x, adi4ibuf, adi4buf(*)
    if (adi4ibuf .ge. 1024) then
        write (6, *) ' Error too small buffer !!! '
    else
        adi4buf(adi4ibuf) = x
        adi4ibuf = adi4ibuf + 1
    end if
end subroutine PUSHINTEGER4

subroutine POPINTEGER4(adi4ibuf, adi4buf, x)
    integer x, adi4ibuf, adi4buf(*)
    if (adi4ibuf .le. 1) then
        write (6, *) ' Error check PUSH/POP !!! '
    else
        adi4ibuf = adi4ibuf - 1
        x = adi4buf(adi4ibuf)
    end if
end subroutine POPINTEGER4

subroutine PUSHREAL8(adr8ibuf, adr8buf, x)
    integer adr8ibuf
    real*8 x, adr8buf(*)
    if (adr8ibuf .ge. 1024) then
        write (6, *) ' Error too small buffer !!! '
    else
        adr8buf(adr8ibuf) = x
        adr8ibuf = adr8ibuf + 1
    end if
end subroutine PUSHREAL8

subroutine POPREAL8(adr8ibuf, adr8buf, x)
    integer adr8ibuf
    real*8 x, adr8buf(*)
    if (adr8ibuf .le. 1) then
        write (6, *) ' Error check PUSH/POP !!! '
    else
        adr8ibuf = adr8ibuf - 1
        x = adr8buf(adr8ibuf)
    end if
end subroutine POPREAL8

