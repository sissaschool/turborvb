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

subroutine branchingo(nw, wconfn, weight, zeta, icdiff, ipip, jbra)

    implicit none

    integer :: nw, i, j, ipip(*), ind, icdiff, ni
    integer :: jbra(*)
    real(8) :: zeta(*), weight, try, tryp, wm, dstep, cost, wconfn(*)
    !
    ! in case of decoupled VMC/DMC calculations, perform
    ! the branching within each pool which contains
    ! nw/nk walkers.
    !
    weight = dabs(wconfn(1))
    do i = 2, nw
        weight = weight + dabs(wconfn(i))
    end do

    !       optimized branching with or without sorttf
    dstep = zeta(nw + 1)/dble(nw)
    !       Optimized branching
    do i = 1, nw
        zeta(i) = weight*(dstep + (i - 1)/dble(nw))
    end do
    zeta(nw + 1) = weight + 1.d0

    ind = 1
    try = 0.d0
    icdiff = 0
    do i = 1, nw
        tryp = try + dabs(wconfn(i))
        ni = 0
        do while (zeta(ind) .lt. tryp .and. zeta(ind) .ge. try)
            jbra(ind) = i
            ind = ind + 1
            ni = ni + 1
        end do
        try = tryp
        if (ni .ne. 0) icdiff = icdiff + 1
    end do

    !       now I have to permute the indices according to jbra
    !       vectorized reshuffling
    !       before I make a permutation to optimize the reshuffling
    call upjbra(nw, jbra, ipip, ipip(nw + 1))

    do i = 1, nw
        wconfn(i) = 1.d0
    end do
    weight = weight/nw

    return
end
