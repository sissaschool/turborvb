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

subroutine ratiofn_psi(parcut, psiold, ratio, tabler, typereg)
    implicit none
    real*8 psiold, ratio, parcut, psinew, tabler, lpsinew, lpsiold
    integer typereg
    !     Input psiold = log(|wfold|),ratio=wfnew/wfold
    !     Output (wfnew/wfnewr)/(wfold/wfoldr)
    !     wf*r means the regularized wf
    !     The sign does not play any role since Sign wf*r= Sign wf*.
    if (ratio .ne. 0.d0) then ! We do not need to introduce any cutoff
        psinew = dlog(abs(ratio)) + psiold

        if (typereg .eq. 2) then

            if (psinew .lt. parcut) then
                tabler = 0.d0
            else
                if (typereg .le. 3) then
                    tabler = dexp(psinew - psiold) ! Here wfnewr = dexp(parcut)
                else
                    tabler = 1.d0
                end if
            end if
        elseif (typereg .eq. 4) then
            !   wf/wfr= wf/sqrt(wf^2 + epsilon^2), epsilon=exp(parcut)
            lpsiold = psiold - parcut - 0.5d0*dlog(dexp(2.d0*(psiold - parcut)) + 1.d0)
            lpsinew = psinew - parcut - 0.5d0*dlog(dexp(2.d0*(psinew - parcut)) + 1.d0)
            tabler = dexp(lpsinew - lpsiold)
        elseif (typereg .eq. 8) then
            !   wf/wfr= |wf|. wfr=sgn (wf)
            tabler = abs(ratio)
        else
            if (psiold .lt. parcut) then ! Here wfoldr= dexp(parcut)
                if (psinew .lt. parcut) then
                    tabler = dexp(psinew - psiold) ! Here wfnewr = dexp(parcut)
                else
                    tabler = exp(parcut - psiold) ! Here wfnewr=wfnew
                end if
            else ! Here wfoldr=wfold
                if (psinew .lt. parcut) then
                    tabler = dexp(psinew - parcut) ! Here wfnewr = dexp(parcut)
                else
                    tabler = 1.d0 ! Here wfnewr=wfnew
                end if
            end if
        end if
    else
        tabler = 0.d0
    end if
    return
end

