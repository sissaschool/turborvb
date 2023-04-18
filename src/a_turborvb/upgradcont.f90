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

subroutine upgradcont(gradpsibar, gradpsi, indt, nelup, neldo         &
        &, winvup, winvdo, tabpip, ttry, gradtotbar, gradtot, rcart, dist          &
        &, rion, nion, zeta, LBox)
    use Cell
    use allio, only: zetamin, distmin
    implicit none
    integer i, jj, jn, iout, indt, nelup, neldo, ion, nion, nel
    real*8 winvup(nelup, *), winvdo(neldo, *), tabpip(nelup + neldo, *)      &
            &, gradpsibar(3, *), gradpsi(3, *), ttry, dist(nion, nelup + neldo)            &
            &, gradtot, gradtotbar, gradpsimod, rion(3, *), rdiff(3)         &
            &, rdiffmod, cutoff, zeta(*), rdiffmod2, gradpsimod2, scale              &
            &, gradpsibarmod2
    real*8 rcart(3, *), LBox, sdiff(3)

    nel = nelup + neldo
    gradtotbar = 0.d0
    gradtot = 0.d0

    do iout = 1, nel

        ! compute gradient for the particle iout
        if (iout .le. nelup) then
            do jn = 1, 3
                jj = indt + jn
                gradpsi(jn, iout) = tabpip(iout, jj) + winvup(iout, jj)
            end do
        else
            do jn = 1, 3
                jj = indt + jn
                gradpsi(jn, iout) = tabpip(iout, jj) + winvdo(iout - nelup, jj)
            end do
        end if

        gradpsimod2 = 0.d0
        do jn = 1, 3
            gradpsimod2 = gradpsimod2 + gradpsi(jn, iout)**2
        end do
        ! updating the modulus square of total gradient
        gradtot = gradtot + gradpsimod2
        gradpsimod = sqrt(gradpsimod2)

        ! find the closest nucleus ion to the iout electron
        ion = 1
        distmin(iout) = dist(ion, iout)
        do jn = 2, nion
            if (dist(jn, iout) .lt. distmin(iout)) then
                ion = jn
                distmin(iout) = dist(jn, iout)
            end if
        end do
        zetamin(iout) = zeta(ion)

        ! rdiff vector connecting ion nucleus to iout electron
        rdiffmod2 = 0.d0
        rdiff(:) = rcart(:, iout) - rion(:, ion)
        if (LBox .gt. 0.d0) then
            call ApplyPBC(rdiff, 1)
        end if
        rdiffmod2 = sum(rdiff(:)**2)
        rdiffmod = sqrt(rdiffmod2)

        ! compute gradpsibar (with cutoff)
        cutoff = 0.d0
        do jn = 1, 3
            cutoff = cutoff + gradpsi(jn, iout)*rdiff(jn)
        end do
        !   Eq.36 Umrigar J.Chem.Phys. vol.99 p.2866
        cutoff = cutoff/rdiffmod/gradpsimod
        cutoff = 0.5d0*(1.d0 + cutoff) + zeta(ion)**2*rdiffmod2/10.d0        &
                & /(4.d0 + zeta(ion)**2*rdiffmod2)

        scale = (sqrt(1.d0 + 4.d0*cutoff*gradpsimod2*ttry) - 1.d0)           &
                & /cutoff/2.d0/gradpsimod2/ttry

        gradpsibarmod2 = 0.d0
        do jn = 1, 3
            gradpsibar(jn, iout) = scale*gradpsi(jn, iout)
            gradpsibarmod2 = gradpsibarmod2 + gradpsibar(jn, iout)**2
        end do

        ! updating the modulus square of total cutoff gradient
        gradtotbar = gradtotbar + gradpsibarmod2

    end do

    return
end
