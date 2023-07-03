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

subroutine eval_iond(iond, rion, nion, LBox, psip, iond_cart)
    use Cell
    implicit none
    integer nion, i, j, d, counter, inds, kkk
    real(8) iond(nion, nion), rion(3, *), LBox, psip(3, *), &
            &iond_cart(3, nion, nion)

    !    cdrion(3,(nion*(nion-1))/2)
    !    1,srion(3,(nion*(nion-1))/2)

    if (LBox .gt. 0.d0) then
        !pbc case
        !vectors between ions
        counter = 0
        do i = 1, nion
            do j = i + 1, nion
                counter = counter + 1
                psip(:, counter) = rion(:, i) - rion(:, j)
            end do
        end do

        call ApplyPBC(psip, counter)

        !distances between ions
        counter = 0
        do i = 1, nion
            do j = i + 1, nion
                counter = counter + 1
                iond(i, j) = dsqrt(sum(psip(:, counter)**2))
                iond_cart(:, i, j) = psip(:, counter)
            end do
        end do
    else
        !open system case
        do i = 1, nion
            do j = i + 1, nion
                !distances between ions
                iond(i, j) = dsqrt(sum((rion(:, i) - rion(:, j))**2))
                iond_cart(:, i, j) = rion(:, i) - rion(:, j)
            end do
        end do
    end if

    ! put the correspondent (j,i) from (i,j)
    do i = 1, nion
        iond(i, i) = 0.d0
        iond_cart(:, i, i) = 0.d0
        !distances between ions
        do j = i + 1, nion
            iond(j, i) = iond(i, j)
            iond_cart(:, j, i) = -iond_cart(:, i, j)
        end do
    end do

    return
end
