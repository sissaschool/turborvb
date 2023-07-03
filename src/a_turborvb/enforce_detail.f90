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

subroutine enforce_detail(indt, nel, identity, gamma, tabler, norm_tab)
    implicit none
    integer i, j, nel, indt
    real*8 gamma, veff, norm_tab, identity, allid, newid, tabler(nel, indt)
    allid = 0.d0
    norm_tab = 0.d0
    veff = 0.d0
    do j = 1, indt - 1
        do i = 1, nel
            if (tabler(i, j) .gt. 0) then
                norm_tab = norm_tab + tabler(i, j)
            else
                norm_tab = norm_tab - gamma*tabler(i, j)
                veff = veff + tabler(i, j)
            end if
            allid = allid + tabler(i, j)
        end do
    end do
    allid = allid/nel
    veff = veff/nel
    newid = max(identity - allid + (1.d0 + gamma)*veff, 0.d0)
    norm_tab = norm_tab + newid*nel
    tabler(:, indt) = newid
    return
end

subroutine cutwstodmc(yesalfe, nel, wsto, lambda, diffkin, vpotge, cutreg)
    use allio, only: zetamin, distmin
    implicit none
    integer nel, j
    logical yesalfe
    real*8 wsto, wstosav, diffkin, vpotge(nel), cutreg, cutregu, costz&
            &, lambda
    if (cutreg .lt. 0.d0) return ! no cutoff
    if (yesalfe) then
        wstosav = wsto
        if (wsto - lambda .gt. cutreg) then
            wsto = cutreg + lambda
        elseif (wsto - lambda .lt. -cutreg) then
            wsto = -cutreg + lambda
        end if
        diffkin = diffkin - wsto + wstosav
    else

        !     if(abs(wsto+sum(vpotge(1:nel))).gt.1d-8) &
        !    &write(6,*) ' ERROR wsto ! ',wsto,sum(vpotge(1:nel))

        do j = 1, nel
            cutregu = -cutreg*zetamin(j)**2*(1 + distmin(j)**2)&
                    & /(1 + (zetamin(j)*distmin(j))**2)
            if (vpotge(j) .lt. cutregu) then
                wsto = wsto + vpotge(j) - cutregu
                diffkin = diffkin - vpotge(j) + cutregu
                !      elseif(vpotge(j).gt.-cutregu) then
                !      wsto=wsto+vpotge(j)+cutregu
                !      diffkin=diffkin-vpotge(j)-cutregu
            end if
        end do
    end if
    return
end

