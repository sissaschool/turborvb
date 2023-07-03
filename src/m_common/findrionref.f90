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

subroutine findrionfref(nion, a, rion, ref, mindist)
    implicit none
    integer i, nion
    real*8 rion(3, nion), ref, a, setint, shiftref, mindist
    real*8, dimension(:), allocatable :: gap, psip
    integer, dimension(:), allocatable :: ipsip

    allocate (gap(nion), ipsip(nion), psip(nion))

    !         write(6,*) ' Find modulo 0,a ',0,a
    do i = 1, nion
        psip(i) = setint(rion(1, i) + ref, a)
        !         write(6,*) i,psip(i)
    end do

    if (nion .eq. 1) then

        gap(1) = a
        ipsip(1) = 1

    else

        call dsortx(psip, 1, nion, ipsip)

        !         write(6,*) ' gap '
        do i = 1, nion - 1
            gap(ipsip(i)) = psip(i + 1) - psip(i)
            !         write(6,*) i,ipsip(i),gap(ipsip(i))
        end do
        gap(ipsip(nion)) = a - psip(nion) + psip(1)
        !         write(6,*) nion,ipsip(nion),gap(ipsip(nion))
        !         find the maximum gap
        !         recompute psip
        do i = 1, nion
            psip(i) = setint(rion(1, i) + ref, a)
        end do

        call dsortx(gap, 1, nion, ipsip)

        !         write(6,*) ' Ordered coordinate/gap '
        !         do i=1,nion
        !         write(6,*) i,ipsip(i),psip(ipsip(i)),gap(i)
        !         enddo

    end if

    shiftref = a - psip(ipsip(nion)) - gap(nion)/2.d0
    shiftref = setint(shiftref, a)
    if (abs(shiftref - a) .lt. shiftref) shiftref = shiftref - a

    mindist = gap(nion)/2.d0

    !         write(6,*) ' shiftref/mindist  inside =',shiftref,mindist

    ref = ref + shiftref

    deallocate (gap, psip, ipsip)

    return
end

function setint(x, a)
    real*8 setint, x, a
    integer n
    if (x .ge. 0.d0) then
        n = x/a
        setint = x - n*a
    else
        n = -x/a
        setint = x + n*a + a
    end if
    return
end

