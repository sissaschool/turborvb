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

function jastrow(rc, vj, iesd, ispin)
    implicit none
    real*8 r, rc(*), rz, jastrow, vj(*)
    integer iesd, ispin
    r = dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2)
    select case (iesd)
        ! exponential form
    case (4)
        if (ispin .lt. 0) then
            jastrow = 0.5d0/vj(1)*(1.d0 - dexp(-vj(1)*r))
        else
            jastrow = 0.25d0/vj(1)*(1.d0 - dexp(-vj(1)*r))
        end if

        jastrow = dexp(jastrow)

        ! 1 parameter 0.5*r/(1+b*r)
    case (1)
        jastrow = 0.5d0*r/(1.d0 + vj(1)*r)
        jastrow = dexp(jastrow)

        ! 2 parameters rescaled 0.5*r/(1+b*r)  r=(1-exp(-k*r))/k
    case (6)
        if (vj(2) .gt. 1d-9) r = (1.d0 - dexp(-vj(2)*r))/vj(2)
        jastrow = 0.5d0*r/(1.d0 + vj(1)*r)
        jastrow = dexp(jastrow)

    case (7)
        if (vj(2) .gt. 1d-9) r = (1.d0 - dexp(-vj(2)*r))/vj(2)
        if (ispin .gt. 0) then
            jastrow = 0.25d0*r/(1.d0 + vj(1)*r)
        else
            jastrow = 0.5d0*r/(1.d0 + vj(1)*r)
        end if
        jastrow = dexp(jastrow)

    case (5)
        jastrow = 0.5d0/(vj(1) + vj(3)*vj(2))* &
                &   (1.d0 + vj(3) - dexp(-vj(1)*r) - vj(3)*dexp(-vj(2)*r))
        jastrow = dexp(jastrow)

        ! spin contaminated: up-up and up-down cusps
        ! 1 parameter
    case (-1)
        ! parallel  spins
        if (ispin .gt. 0) then
            jastrow = 0.25d0*r/(1.d0 + vj(1)*r)
        else
            ! opposite  spins
            jastrow = 0.5d0*r/(1.d0 + vj(1)*r)
        end if
        jastrow = dexp(jastrow)
        ! asymmetric jastrow (xy ne z component)
    case (-2)
        rz = dsqrt(vj(1)**2*(rc(1)**2 + rc(2)**2) + (vj(2)*rc(3))**2)
        jastrow = 0.5d0*r/(1.d0 + rz)
        jastrow = dexp(jastrow)

    case (3)
        jastrow = 0.5d0*r*(1.d0/(1.d0 + vj(1)*r)                           &
                & + vj(3)*r/(1.d0 + vj(2)*r)**2)
        jastrow = dexp(jastrow)

        ! spin contaminated: up-up and up-down cusps
        ! 2 parameters
    case (2)
        ! parallel  spins
        if (ispin .gt. 0) then
            jastrow = 0.25d0*r/(1.d0 + vj(2)*r)
            ! opposite spins
        else
            jastrow = 0.5d0*r/(1.d0 + vj(1)*r)
        end if
        jastrow = dexp(jastrow)

        ! no jastrow
    case (0)
        jastrow = 1.d0
    end select

    return
end
