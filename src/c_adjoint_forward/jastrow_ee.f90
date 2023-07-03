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

function jastrow_ee(rc, vj, iesd, ispin)
    use cell, only: metric
    use allio, only: iespbc, norm_metric
    implicit none
    real*8 r, rs, rc(3), rz, jastrow_ee, vj(*)
    integer iesd, ispin, j
    if (iespbc) then
        r = norm_metric(rc, metric)
    else
        r = dsqrt(sum(rc(1:3)**2))
    end if
    select case (iesd)
        ! exponential form
    case (4)
        if (ispin .lt. 0) then
            jastrow_ee = 0.5d0/vj(1)*(1.d0 - dexp(-vj(1)*r))
        else
            jastrow_ee = 0.25d0/vj(1)*(1.d0 - dexp(-vj(1)*r))
        end if
    case (-4)
        !     No spin contamination
        jastrow_ee = 0.5d0/vj(1)*(1.d0 - dexp(-vj(1)*r))
    case (1)
        ! 1 parameter 0.5*r/(1+b*r)
        jastrow_ee = 0.5d0*r/(1.d0 + vj(1)*r)
    case (6)
        ! 2 parameters rescaled 0.5*r/(1+b*r)  r=(1-exp(-k*r))/k
        if (vj(2) .gt. 1d-9) then
            rs = (1.d0 - dexp(-vj(2)*r))/vj(2)
        else
            rs = r
        end if
        jastrow_ee = 0.5d0*rs/(1.d0 + vj(1)*rs)
    case (8)
        jastrow_ee = 1.d0/vj(1)*(1.d0 - dexp(-vj(1)*r**3))
    case (9)
        !       RVB log jastrow r-->0  Const. +r/2 ,r --> infty --> -Const. >0 log r
        jastrow_ee = -log(1 + vj(1)*(1 - r/vj(2))**2)*vj(2)&
                    & *(1.d0 + vj(1))/(4.d0*vj(1))
    case (10) ! The same with cusp condition for equal spin electrons r-> r/2
        if (ispin .gt. 0) then
            jastrow_ee = -log(1 + vj(1)*(1 - r/(2.d0*vj(2)))**2)*vj(2)&
                        & *(1.d0 + vj(1))/(4.d0*vj(1))
        else
            jastrow_ee = -log(1 + vj(1)*(1 - r/vj(2))**2)*vj(2)&
                        & *(1.d0 + vj(1))/(4.d0*vj(1))
        end if
    case (7)
        if (ispin .gt. 0) then
            jastrow_ee = 0.25d0*r/(1.d0 + vj(1)*r)
        else
            jastrow_ee = 0.5d0*r/(1.d0 + vj(1)*r)
        end if
    case (5)
        jastrow_ee = 0.5d0/(vj(1) + vj(3)*vj(2))* &
                &   (1.d0 + vj(3) - dexp(-vj(1)*r) - vj(3)*dexp(-vj(2)*r))

        ! spin contaminated: up-up and up-down cusps
        ! 1 parameter equal to -7
    case (-1)
        ! parallel  spins
        if (ispin .gt. 0) then
            jastrow_ee = 0.25d0*r/(1.d0 + vj(1)*r)
        else
            ! opposite  spins
            jastrow_ee = 0.5d0*r/(1.d0 + vj(1)*r)
        end if
    case (-2)
        ! asymmetric jastrow_ee (xy ne z component)
        rz = dsqrt(vj(1)**2*(rc(1)**2 + rc(2)**2) + (vj(2)*rc(3))**2)
        jastrow_ee = 0.5d0*r/(1.d0 + rz)

    case (3)
        jastrow_ee = 0.5d0*r*(1.d0/(1.d0 + vj(1)*r)                        &
                & + vj(3)*r/(1.d0 + vj(2)*r)**2)

        ! spin contaminated: up-up and up-down cusps
        ! 2 parameters
    case (2)
        ! parallel  spins
        if (ispin .gt. 0) then
            jastrow_ee = 0.25d0*r/(1.d0 + vj(2)*r)
            ! opposite spins
        else
            jastrow_ee = 0.5d0*r/(1.d0 + vj(1)*r)
        end if

        ! no jastrow_ee
    case (0)
        jastrow_ee = 0.d0
    end select
    return
end
