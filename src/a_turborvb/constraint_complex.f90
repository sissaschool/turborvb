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

subroutine constraint_complex(iessw, derl, n, nozero&
        &, psip, econf, nw, jbradet)
    use constants, only: ipc, ipf, deps
    use cell, only: cellscale, phase2pi
    use allio, only: rank, kiontot, rion, yes_crystal, ndim_detmat&
            &, symmagp, yes_hermite, yes_correct, nelorb_at, pfaffup, nnozero_eagp&
            &, eagp_pfaffb, ndiff
    implicit none
    integer ind_transp, ind, ndim, iessw, nozero(*), i, ii, j, ix, iy, n, nw, jbradet(*)
    real*8 econf(nw, *), psip(*), derl(*)

    if (iessw .eq. 0) return

    do i = 1, iessw
        econf(1, i) = 0.d0
    end do
    if (ipc .eq. 1) then
        do i = 1, n
            psip(i) = derl(nozero(i))
        end do
        do i = 1, nnozero_eagp
            ind = n + i
            iy = (nozero(ind) - 1)/ndiff + 1
            ix = nozero(ind) - (iy - 1)*ndiff
            psip(ind) = eagp_pfaffb(ix, iy)
        end do
        do j = 1, n + nnozero_eagp
            i = jbradet(j)
            if (i .gt. 0) then
                econf(1, i) = econf(1, i) + psip(j)
            elseif (i .lt. 0) then
                econf(1, -i) = econf(1, -i) - psip(j)
            end if
        end do
    else

        if (symmagp .and. yes_correct) then
            do i = 1, n
                iy = (nozero(i) - 1)/ndim_detmat + 1
                ix = nozero(i) - (iy - 1)*ndim_detmat
                if (ipf .eq. 2) then
                    ndim = nelorb_at/2
                    if (ix .le. ndim .and. iy .le. ndim .and. .not. pfaffup) then
                        ind_transp = (iy + ndim - 1)*ndim_detmat + ix + ndim
                    elseif (iy .gt. ndim .and. ix .le. ndim .and. iy .le. ndim_detmat) then
                        ind_transp = (ix + ndim - 1)*ndim_detmat + iy - ndim
                    else
                        ind_transp = -nozero(i) ! this case should be treated as normal
                    end if
                else
                    ind_transp = (ix - 1)*ndim_detmat + iy
                end if
                if (iy .le. ndim_detmat .and. ind_transp .ne. -nozero(i)) then
                    psip(4*i - 3) = derl(2*nozero(i) - 1)
                    psip(4*i - 2) = derl(2*nozero(i))
                    psip(4*i - 1) = derl(2*ind_transp - 1)
                    psip(4*i) = derl(2*ind_transp)
                else
                    psip(4*i - 3) = derl(2*nozero(i) - 1)
                    psip(4*i - 2) = derl(2*nozero(i))
                    psip(4*i - 1) = 0.d0
                    psip(4*i) = 0.d0
                end if
                if (ipf .eq. 2) then
                    if (iy .gt. ndim .and. ix .le. ndim .and. iy - ndim .eq. ix) psip(4*i - 3:4*i) = psip(4*i - 3:4*i)/2.d0
                else
                    if (ix .eq. iy) psip(4*i - 3:4*i) = psip(4*i - 3:4*i)/2.d0 ! to avoid double counting
                end if
            end do
            do i = 1, nnozero_eagp
                ind = n + i
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                ! The ghost are always considered normal variables with no Hermitian relation.
                psip(4*ind - 3) = eagp_pfaffb(2*ix - 1, iy)
                psip(4*ind - 2) = eagp_pfaffb(2*ix, iy)
                psip(4*ind - 1) = 0.d0
                psip(4*ind) = 0.d0
            end do
        else
            do i = 1, n
                psip(2*i - 1) = derl(2*nozero(i) - 1)
                psip(2*i) = derl(2*nozero(i))
            end do
            do i = 1, nnozero_eagp
                ind = n + i
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                psip(2*ind - 1) = eagp_pfaffb(2*ix - 1, iy)
                psip(2*ind) = eagp_pfaffb(2*ix, iy)
            end do
        end if

        do j = 1, n
            i = jbradet(j)
            iy = (nozero(j) - 1)/ndim_detmat + 1
            ix = nozero(j) - (iy - 1)*ndim_detmat
            if (i .gt. 0) then
                if (symmagp .and. yes_correct) then
                    econf(1, 4*i - 3) = econf(1, 4*i - 3) + psip(4*j - 3) + psip(4*j - 1)
                    econf(1, 4*i - 2) = econf(1, 4*i - 2) + psip(4*j - 2) + psip(4*j)
                    if (yes_hermite) then
                        econf(1, 4*i - 1) = econf(1, 4*i - 1) + psip(4*j - 2) - psip(4*j)
                        econf(1, 4*i) = econf(1, 4*i) - psip(4*j - 3) + psip(4*j - 1)
                    else
                        econf(1, 4*i - 1) = econf(1, 4*i - 1) + psip(4*j - 2) + psip(4*j)
                        econf(1, 4*i) = econf(1, 4*i) - psip(4*j - 3) - psip(4*j - 1)
                    end if
                else
                    econf(1, 2*i - 1) = econf(1, 2*i - 1) + psip(2*j - 1)
                    econf(1, 2*i) = econf(1, 2*i) + psip(2*j)
                end if
            elseif (i .lt. 0) then
                if (symmagp .and. yes_correct) then
                    if (yes_hermite) then
                        econf(1, -4*i - 3) = econf(1, -4*i - 3) - psip(4*j - 3) - psip(4*j - 1)
                        econf(1, -4*i - 2) = econf(1, -4*i - 2) - psip(4*j - 2) - psip(4*j)
                        econf(1, -4*i - 1) = econf(1, -4*i - 1) - psip(4*j - 2) + psip(4*j)
                        econf(1, -4*i) = econf(1, -4*i) + psip(4*j - 3) - psip(4*j - 1)
                    else
                        econf(1, -4*i - 3) = econf(1, -4*i - 3) - psip(4*j - 3) - psip(4*j - 1)
                        econf(1, -4*i - 2) = econf(1, -4*i - 2) - psip(4*j - 2) - psip(4*j)
                        econf(1, -4*i - 1) = econf(1, -4*i - 1) - psip(4*j - 2) - psip(4*j)
                        econf(1, -4*i) = econf(1, -4*i) + psip(4*j - 3) + psip(4*j - 1)
                    end if
                else
                    econf(1, -2*i - 1) = econf(1, -2*i - 1) - psip(2*j - 1)
                    econf(1, -2*i) = econf(1, -2*i) - psip(2*j)
                end if
            end if
        end do
        do ii = 1, nnozero_eagp
            j = n + ii
            i = jbradet(j)
            if (i .gt. 0) then
                if (symmagp .and. yes_correct) then
                    econf(1, 4*i - 3) = econf(1, 4*i - 3) + psip(4*j - 3)
                    econf(1, 4*i - 2) = econf(1, 4*i - 2) + psip(4*j - 2)
                    econf(1, 4*i - 1) = 0.d0
                    econf(1, 4*i) = 0.d0
                else
                    econf(1, 2*i - 1) = econf(1, 2*i - 1) + psip(2*j - 1)
                    econf(1, 2*i) = econf(1, 2*i) + psip(2*j)
                end if
            elseif (i .lt. 0) then
                if (symmagp .and. yes_correct) then
                    econf(1, -4*i - 3) = econf(1, -4*i - 3) - psip(4*j - 3)
                    econf(1, -4*i - 2) = econf(1, -4*i - 2) - psip(4*j - 2)
                    econf(1, -4*i - 1) = 0.d0
                    econf(1, -4*i) = 0.d0
                else
                    econf(1, -2*i - 1) = econf(1, -2*i - 1) - psip(2*j - 1)
                    econf(1, -2*i) = econf(1, -2*i) - psip(2*j)
                end if
            end if
        end do
    end if ! endif ipc
    return
end
