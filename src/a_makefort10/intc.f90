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

character(6) function intc(a)
    !
    character(6) temp
    real*8, intent(in) :: a
    integer :: i
    real*8 :: thr = 1.0d-9
    !
    i = int(a) ! not nint!! because we want to distinguish 1.00, 1.01, .... 1.99.

    if (i .ne. a) then
        ! K.N. 12 May 2022. Here, not int but nint is needed!! for avoiding precision underflow.
        ! I noticed that int does not work for a=2 (He atom).
        !if(a*10.eq.int(a*10)) then
        if (abs((nint(a*10) - a*10)) .le. thr) then
            if (i .ge. 0) then
                if (i .lt. 10) then
                    write (temp, '(f3.1)') a
                elseif (i .lt. 100) then
                    write (temp, '(f4.1)') a
                elseif (i .lt. 1000) then
                    write (temp, '(f5.1)') a
                elseif (i .lt. 10000) then
                    write (temp, '(f6.1)') a
                end if
            else
                if (-i .lt. 10) then
                    write (temp, '(f4.1)') a
                elseif (-i .lt. 100) then
                    write (temp, '(f5.1)') a
                elseif (-i .lt. 1000) then
                    write (temp, '(f6.1)') a
                end if
            end if
            ! K.N. 12 May 2022. Here, not int but nint is needed!! for avoiding precision underflow.
            ! I noticed that int does not work for a=2 (He atom).
            !elseif(a*100.eq.int(a*100)) then
        elseif (abs((nint(a*100) - a*100)) .le. thr) then
            if (i .ge. 0) then
            if (i .lt. 10) then
                write (temp, '(f4.2)') a
            elseif (i .lt. 100) then
                write (temp, '(f5.2)') a
            elseif (i .lt. 1000) then
                write (temp, '(f6.2)') a
            end if
            else
            if (-i .lt. 10) then
                write (temp, '(f5.2)') a
            elseif (-i .lt. 100) then
                write (temp, '(f6.2)') a
            end if
            end if
        else
        end if

    else

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

    end if

    intc = temp
    !
end function intc
