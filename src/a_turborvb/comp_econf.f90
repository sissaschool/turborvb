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

!-------------------------------------------------------------------------
! subroutines comp_econf/comp_econf_complex compute the energy and other
! correlation functions after one loop over the walkers
!-------------------------------------------------------------------------
! Main variables:
! econf = array containing the values of energy and correlations functions
! iese =  determines which parts of the total energy are stored
!         iese=1 compute the total energy only (default option)
! isfix = store the variance <H^2> if different from 0

subroutine comp_econf(j, js, econf, enert, diffkin, vpot, vcut, voffpseudo, table)

    use constants, only: ipc
    use allio, only: itest, iese, in1, yes_complex, nel, nelup, neldo, indt, &
                     psip, indtj, istart, gamma, ivic, indksj, isfix, rank

    implicit none

    integer j, js, kk
    real(8) econf(*), diffkin(3, *), enert(ipc, *), vpot(*), &
        voffpseudo, table(*), vcut(*)

    if (iese .ne. 0) then
        ! energy
        econf(j) = enert(1, j)
        if (ipc .eq. 2) econf(j + in1) = enert(2, j)
        if (iese .ge. 1 + ipc) then
            econf(j + in1*ipc) = diffkin(2, j)*0.5d0 ! unit Ha
        end if
        if (iese .ge. 2 + ipc) then
            call diffus(nel, indt, gamma, ivic(1, 1, indksj), table(indtj), psip, istart)
            econf(j + (1 + ipc)*in1) = psip(1)
        end if
        if (iese .eq. 3 + ipc) then
            econf(j + (2 + ipc)*in1) = voffpseudo*0.5d0
        end if
    end if
    ! correlation functions
    if (iese .ne. 3 + ipc) then
        if (iese .ge. 3 + ipc) econf(j + (2 + ipc)*in1) = vcut(j)
        if (iese .ge. 4 + ipc) econf(j + (3 + ipc)*in1) = diffkin(1, j)*0.5d0 !unit Ha
        if (iese .ge. 5 + ipc) econf(j + (4 + ipc)*in1) = diffkin(2, j)*0.5d0 !unit Ha
    end if

    return
end subroutine comp_econf
