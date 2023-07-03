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

!------------------------------------------------------------------------
! subroutines updiag/updiag_complex compute the local energy both
! diagonal and off-diagonal part. Regularization of the Coulomb potential
! is also computed.
!------------------------------------------------------------------------

subroutine updiag(table, tmu, diag, enerc, winvup, winvdo, tabpip&
        &, nelup, neldo, nel, indt, indtupt, indteff, istart, vpot, vpotregr, parcutg&
        &, vcut, kin, novar, flagsign, enertrue, vpotoff, kince, vpotge)
    use allio, only: itestrfn, fncont, safelrdmc
    implicit none
    integer i, nelup, neldo, indt, indtupt, jj, ncore                                &
            &, nel, istart, istartmin, indteff, novar, parcutg
    real*8 table(nel, *), tmu(nel, *), enerc, etry, enereff                          &
            &, winvup(nelup, *), winvdo(neldo, *), tabpip(*), kinc                      &
            &, diag, vpot, psiln, costtmu, enertrue_sav, enertrue, vpotoff, vpotc, vpotreg &
            &, diagr, vpotge(nel), vpotchosen, vcut, kin(3), vpotregr(2, nel)&
            &, vpotcreg, costtmue, vpotg, cost, kince(nel)
    logical flagsign, flagsignt, flagparg

    call upkinc(indt, nelup, neldo, winvup, winvdo, tabpip, kinc, kince)
    vpotreg = sum(vpotregr(2, :))
    vpotoff = 0.d0
    do jj = istart, indtupt ! to compute  the off diagonal  pseudo part put in vpot
        do i = 1, nel
            vpotoff = vpotoff - table(i, jj)
        end do
    end do
    enertrue = kinc + vpot + vpotoff
    enerc = enertrue

    if (indteff .eq. istart - 1) then ! locality approximation ill defined
        ! left for tests
        vpotc = vpot + vpotoff
        vpotcreg = vpotreg + vpotoff
    else
        vpotc = vpot
        vpotcreg = vpotreg
    end if
    costtmu = 0.d0
    do jj = 1, istart - 1
        do i = 1, nel
            costtmu = costtmu + tmu(i, jj)
        end do
    end do

    diag = -kinc - vpotc

    !     for dmclrdmc
    do jj = indtupt + 1, indt
        do i = 1, nel
            diag = diag - table(i, jj)
        end do
    end do

    do i = 1, nel
        vpotge(i) = kince(i) + vpotregr(1, i)
    end do
    do i = 1, nel
        cost = 0.d0
        costtmue = 0.d0
        do jj = 1, istart - 1 ! do for all indt
            costtmue = costtmue + tmu(i, jj)
            cost = cost - table(i, jj)
        end do
        diag = diag + cost
        vpotge(i) = vpotge(i) - cost - costtmue
    end do
    diagr = diag

    vpotg = sum(vpotge(:))

    !     The effective potential is ill defined if it can be infinite.

    !     vpotcreg = bar Coulomb regularized this is finite when a node is approached
    !     vpotg = -diag-costtmu    ! this is finite when a nuclei is approached
    !                                if psi_T satisfies the e-i/e-e cusp conditions
    !                                as is  always the case in TurboRVB

    !     vpotchosen is bounded from below and H^a has a well defined ground state
    !     H^a--> H +O(a^2) almost everywhere apart in the nodal surface or
    !     when kincut = parcutg.
    !     vcut: the right quantity that is sensitive to the proximity to the node

    if (abs(parcutg) .eq. 2) then
        vcut = vpotg - vpotcreg
        kin(1) = diagr + costtmu + vpotcreg + kinc ! k^a
    else
        vcut = vpotg - vpotc
        kin(1) = diagr + costtmu + vpotc + kinc ! k^a
    end if
    kin(2) = kinc ! k

    ! parcutg=0  no cutoff
    ! parcutg=ne 0  choose max(vpotc,vpota) only when sign in the local part

    flagsignt = .false.

    if (parcutg .ne. 0) then

        vpotchosen = 0.d0
        istartmin = min(6, istart - 1)
        do i = 1, nel
            !       check if one electron cross the nodal surface
            flagsign = .false.
            do jj = 1, istartmin
                if (table(i, jj) .lt. 0.00) then
                    flagsign = .true.
                    flagsignt = .true.
                end if
            end do
            if (safelrdmc) flagsign = .true. ! apply bound for any conf to be sure H_a is bounded from below

            if (flagsign) then
                if (abs(parcutg) .eq. 1) then
                    if (vpotge(i) .lt. vpotregr(1, i)) then
                        vpotchosen = vpotchosen + vpotregr(1, i)
                    else
                        vpotchosen = vpotchosen + vpotge(i)
                    end if
                elseif (abs(parcutg) .eq. 2) then
                    if (vpotge(i) .lt. vpotregr(2, i)) then
                        vpotchosen = vpotchosen + vpotregr(2, i)
                    else
                        vpotchosen = vpotchosen + vpotge(i)
                    end if
                else
                    vpotchosen = vpotchosen + vpotge(i)
                end if
            else
                vpotchosen = vpotchosen + vpotge(i)
            end if
        end do

        diag = -vpotchosen - costtmu
        enereff = enerc - vpotg + vpotchosen ! always the local energy of H^a
        if (novar .eq. 1) then
            enertrue = enertrue - vpotg + vpotchosen
            enerc = enerc - vpotg + vpotchosen
        end if
        if (novar .eq. 2) then
            enertrue = enertrue + 0.5d0*(-vpotg + vpotchosen)
            enerc = enerc - vpotg + vpotchosen
        end if
    else
        enereff = enerc
    end if

    flagsign = flagsignt

    if (abs(parcutg) .le. 3) then
        kin(3) = enereff ! The local energy of H^a
    else
        kin(3) = enertrue ! The local energy of H^a
    end if

    if (fncont .and. itestrfn .ne. -3) then
        do jj = istart, indtupt ! to compute  the off diagonal  pseudo part put in vpot
            do i = 1, nel
                vpotge(i) = vpotge(i) - table(i, jj)
            end do
        end do
    end if

    !stop

    return
end subroutine updiag

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!

subroutine updiag_complex(table, tmu, diag, enerc, winvup, winvdo, tabpip &
                          , nelup, neldo, nel, indt, indtupt, indteff, istart, vpot, vpotregr, parcutg &
                          , vcut, kin, novar, flagsign, enertrue, vpotoff, kince, vpotge)
    use allio, only: itestrfn, fncont, safelrdmc
    use Constants, only: zone, zzero, ipc

    implicit none
    integer i, nelup, neldo, indt, indtupt, jj, ncore &
        , nel, istart, istartmin, indteff, novar, parcutg, ip4
    real(8) tmu(nel, *), etry, tabpip(*), vpot, psiln, costtmu, vpotreg &
        , vpotchosen, vcut, vpotregr(2, nel), costtmue
    logical flagsign, flagsignt, flagparg
    ! complex variables
    complex(8) table(nel, *), winvup(nelup, *), winvdo(neldo, *), kinc, kince(nel), enertrue, vpotoffc, enerc
    real(8) enereff, diag, diagr, kin(3), &
        vpotoff, vpotc, vpotcreg, vpotge(nel), vpotg, cost

    call upkinc_complex(indt, nelup, neldo, winvup, winvdo, tabpip, kinc, kince)

    vpotreg = sum(vpotregr(2, :))
    vpotoffc = zzero
    do jj = istart, indtupt ! to compute the off diagonal pseudo part put in vpot
        do i = 1, nel
            vpotoffc = vpotoffc - table(i, jj)
        end do
    end do
    enertrue = kinc + vpot + vpotoffc ! total energy
    enerc = enertrue
    vpotoff = vpotoffc

    if (indteff .eq. istart - 1) then ! locality approximation ill defined
        ! left for tests
        vpotc = vpot + vpotoffc
        vpotcreg = vpotreg + vpotoffc
    else
        vpotc = vpot
        vpotcreg = vpotreg
    end if
    costtmu = 0.d0
    do jj = 1, istart - 1
        do i = 1, nel
            costtmu = costtmu + tmu(i, jj)
        end do
    end do

    diag = -kinc - vpotc ! local (diagonal) part of the Hamiltonian

    ! for dmclrdmc
    do jj = indtupt + 1, indt
        do i = 1, nel
            diag = diag - dreal(table(i, jj))
        end do
    end do

    do i = 1, nel
        vpotge(i) = dreal(kince(i)) + vpotregr(1, i)
    end do
    do i = 1, nel
        cost = 0.d0
        costtmue = 0.d0
        do jj = 1, istart - 1
            costtmue = costtmue + tmu(i, jj)
            cost = cost - table(i, jj)
        end do
        diag = diag + cost
        vpotge(i) = vpotge(i) - cost - costtmue
    end do
    diagr = diag
    vpotg = sum(vpotge(:))

    !     The effective potential is ill defined if it can be infinite.

    !     vpotcreg = bar Coulomb regularized this is finite when a node is approached
    !     vpotg = -diag-costtmu    ! this is finite when a nuclei is approached
    !                                if psi_T satisfies the e-i/e-e cusp conditions
    !                                as is  always the case in TurboRVB

    !     vpotchosen is bounded from below and H^a has a well defined ground state
    !     H^a--> H +O(a^2) almost everywhere apart in the nodal surface or
    !     when kincut = parcutg.
    !     vcut: the right quantity that is sensitive to the proximity to the node

    if (abs(parcutg) .eq. 2) then
        vcut = vpotg - vpotcreg
        kin(1) = diagr + costtmu + vpotcreg + kinc ! k^a
    else
        vcut = vpotg - vpotc
        kin(1) = diagr + costtmu + vpotc + kinc ! k^a
    end if

    kin(2) = kinc ! k
    flagsignt = .false.

    ! for dmclrdmc
    if (parcutg .ne. 0) then

        vpotchosen = 0.d0
        istartmin = min(6, istart - 1)
        do i = 1, nel
            ! check if one electron cross the nodal surface
            flagsign = .false.
            do jj = 1, istartmin
                if (dreal(table(i, jj)) .lt. 0.d0) then
                    flagsign = .true.
                    flagsignt = .true.
                end if
            end do
            if (safelrdmc) flagsign = .true.
            if (flagsign) then
                if (abs(parcutg) .eq. 1) then
                    if (vpotge(i) .lt. vpotregr(1, i)) then
                        vpotchosen = vpotchosen + vpotregr(1, i)
                    else
                        vpotchosen = vpotchosen + vpotge(i)
                    end if
                elseif (abs(parcutg) .eq. 2) then
                    if (vpotge(i) .lt. vpotregr(2, i)) then
                        vpotchosen = vpotchosen + vpotregr(2, i)
                    else
                        vpotchosen = vpotchosen + vpotge(i)
                    end if
                else
                    vpotchosen = vpotchosen + vpotge(i)
                end if
            else
                vpotchosen = vpotchosen + vpotge(i)
            end if
        end do

        diag = -vpotchosen - costtmu
        enereff = enerc - vpotg + vpotchosen ! always the local energy of H^a
        if (novar .eq. 1) then
            enertrue = enertrue - vpotg + vpotchosen
            enerc = enerc - vpotg + vpotchosen
        end if
        if (novar .eq. 2) then
            enertrue = enertrue + 0.5d0*(-vpotg + vpotchosen)
            enerc = enerc + 0.5d0*(-vpotg + vpotchosen)
        end if
    else
        enereff = enerc
    end if

    flagsign = flagsignt

    if (abs(parcutg) .le. 3) then
        kin(3) = enereff ! The local energy of H^a
    else
        kin(3) = enertrue ! The local energy of H^a
    end if
    if (fncont .and. itestrfn .ne. -3) then
        do jj = istart, indtupt ! to compute  the off diagonal  pseudo part put in vpot
            do i = 1, nel
                vpotge(i) = vpotge(i) - table(i, jj)
            end do
        end do
    end if
    return

end subroutine updiag_complex
