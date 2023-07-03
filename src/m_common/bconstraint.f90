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

subroutine bconstraint(iessw, derl, nelorb, n, nozero              &
        &, psip, econf, nw, jbradet, symmagp, yes_update)
    use constants, only: ipc, ipf, deps
    use cell, only: cellscale, phase2pi
    use allio, only: rank, yes_correct, nnozero_eagp, eagp_pfaff, ndiff
    implicit none
    !        This subroutine translate back the change in the
    !        original unconstrained representation
    integer iessw, nozero(*), i, ix, iy, j, n, nw, nelorb, ind&
            &, jbradet(*), iesswread
    real*8 econf(nw, *), psip(*), derl(*)
    logical symmagp, yes_update
    if (iessw .eq. 0) return

    if (ipc .eq. 1) then

        call dscalzero(n + nnozero_eagp, 0.d0, psip, 1)
        do i = 1, n + nnozero_eagp
            j = jbradet(i)
            if (j .gt. 0) then
                psip(i) = econf(1, j)
            elseif (j .lt. 0) then
                psip(i) = -econf(1, -j)
            end if
        end do

        if (yes_update) then
            do i = 1, n
                if (psip(i) .ne. 0.d0)                                         &
                        &  call upsimp(derl, nelorb, nozero(i), psip(i), symmagp, ipf)
            end do
            do i = 1, nnozero_eagp
                ind = i + n
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    eagp_pfaff(ix, iy) = eagp_pfaff(ix, iy) + psip(ind)
                    eagp_pfaff(iy, ix) = -eagp_pfaff(ix, iy)
                elseif (j .lt. 0) then
                    eagp_pfaff(ix, iy) = eagp_pfaff(ix, iy) - psip(ind)
                    eagp_pfaff(iy, ix) = -eagp_pfaff(ix, iy)
                end if
            end do

        else
            do i = 1, n
                if (abs(jbradet(i)) .ne. 0)                                    &
                        &  call upsim(derl, nelorb, nozero(i), psip(i), symmagp, ipf)
            end do
            do i = 1, nnozero_eagp
                ind = i + n
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    eagp_pfaff(ix, iy) = psip(ind)
                    eagp_pfaff(iy, ix) = -psip(ind)
                elseif (j .lt. 0) then
                    eagp_pfaff(ix, iy) = -psip(ind)
                    eagp_pfaff(iy, ix) = psip(ind)
                end if
            end do
        end if

    else

        call dscalzero(2*n + 2*nnozero_eagp, 0.d0, psip, 1)
        if (symmagp .and. yes_correct) then
            do i = 1, n + nnozero_eagp
                j = jbradet(i)
                if (j .gt. 0) then
                    psip(2*i - 1) = econf(1, 4*j - 3)
                    psip(2*i) = econf(1, 4*j - 1)
                elseif (j .lt. 0) then
                    psip(2*i - 1) = -econf(1, -4*j - 3)
                    psip(2*i) = -econf(1, -4*j - 1)
                end if
            end do

        else

            do i = 1, n + nnozero_eagp
                j = jbradet(i)
                if (j .gt. 0) then
                    psip(2*i - 1) = econf(1, 2*j - 1)
                    psip(2*i) = econf(1, 2*j)
                elseif (j .lt. 0) then
                    psip(2*i - 1) = -econf(1, -2*j - 1)
                    psip(2*i) = -econf(1, -2*j)
                end if
            end do

        end if

        if (yes_update) then
            do i = 1, n
                if (sum(abs(psip(2*i - 1:2*i))) .ne. 0.d0)                  &
                        & call upsimp_complex(derl, nelorb, nozero(i), psip(2*i - 1), symmagp, ipf)
            end do
            do i = 1, nnozero_eagp
                ind = i + n
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    eagp_pfaff(2*ix - 1, iy) = eagp_pfaff(2*ix - 1, iy) + psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = eagp_pfaff(2*ix, iy) + psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)
                elseif (j .lt. 0) then
                    eagp_pfaff(2*ix - 1, iy) = eagp_pfaff(2*ix - 1, iy) - psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = eagp_pfaff(2*ix, iy) - psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)
                end if
            end do

        else
            do i = 1, n
                if (abs(jbradet(i)) .ne. 0)&
                        & call upsim_complex(derl, nelorb, nozero(i), psip(2*i - 1), symmagp, ipf)
            end do

            do i = 1, nnozero_eagp
                ind = i + n
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    eagp_pfaff(2*ix - 1, iy) = psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)
                elseif (j .lt. 0) then
                    eagp_pfaff(2*ix - 1, iy) = -psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = -psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)
                end if
            end do

        end if ! endif yesupdate

    end if ! endif ipc=1

    return
end
