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

subroutine hopping(nion, jn, dstep, rcart                         &
        &, rion, imin, itry, tjasder, LBox)
    use Cell
    implicit none
    integer nion, jn, i, itry, imin
    real(8) rcart(3), LBox, scart(3)
    real(8) rion(3, *), dstep, drand1, tjasder, zeta

    !       Make an hopping move with probability 50%
    !       In the remaining case make the usual move around the atom

    zeta = drand1()

    if (zeta .lt. abs(tjasder) .and. nion .gt. 1) then

        !       Generate  a random itry ne imin

        !       imin=nion*drand1()+1

        itry = (nion - 1)*drand1() + 1
        if (itry .ge. imin) itry = itry + 1

        !     In order to satisfy detailed balance both the
        !     hopping directions  have to be taken into account.
        if (drand1() .gt. 0.5d0) then
            do i = 1, 3
                rcart(i) = rcart(i) + rion(i, itry) - rion(i, imin)
            end do
        else
            do i = 1, 3
                rcart(i) = rcart(i) - rion(i, itry) + rion(i, imin)
            end do
        end if
        dstep = 0.d0
    else
        rcart(jn) = rcart(jn) + dstep
    end if

    !********** Periodic Boundary Conditions **********
    !       if(LBox.gt.0) call ApplyPBC(rcart,1)
    !********** Periodic Boundary Conditions **********

    return
end
