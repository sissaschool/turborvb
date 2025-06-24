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

!> @brief      Determine the sign of a wave function value
!> @details    This function determines the sign of a wave function value based on
!>             the precision mode (ipc). For real wave functions (ipc=1), it returns
!>             the value directly. For complex wave functions (ipc=2), it determines
!>             the sign based on the phase modulo pi.
!> @param[in]  psisn    Wave function value to determine sign for
!> @return     wf_sign  Sign of the wave function value (1 or -1 for complex, actual value for real)
integer function wf_sign(psisn)
    use constants, only: ipc, pi
    implicit none
    real*8 psisn
    if (ipc .eq. 1) then
        wf_sign = psisn
    else
        if (mod(nint(mod(abs(psisn), pi)), 2) .eq. 0) then
            wf_sign = 1
        else
            wf_sign = -1
        end if
    end if
    return
end

