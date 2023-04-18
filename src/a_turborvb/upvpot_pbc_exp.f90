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

subroutine upvpot(npar, iesdr, nw, nel, nelup, rcart, econf, vj          &
        &, rion, nion, costz, costz3, LBox)
    use Constants
    use Cell
    implicit none
    integer npar, nw, L, nel, keli, kelj, iesdr2iesd    &
            &, k, i, ii, jj, ind, indo, iopt, j, jold, iesd, iesdr, nelup, nion
    real*8 xi, yi, zi, xj, yj, zj, nr, nrz, dx, dy, dz
    real*8 econf(nw, *), vj(*), cost, fcost                               &
            &, fcost1, fcost2, fcost3, fat1, fat2, fat3, fat4, costk                   &
            &, costz(*), costz3(*), rion(3, *), LBox
    real*8 rcart(3, *), PI_L, INV_PI_L, rc(3), rt(3), src(3)
    logical iesspin

    !     iesd=iesdr
    !     if(iesdr.eq.-5) iesd=1
    !     if(iesdr.eq.-6) iesd=4
    !     if(iesdr.eq.-7) iesd=4

    iesd = iesdr2iesd(iesdr)

    if (iesd .lt. 0 .or. iesd .eq. 2 .or. iesdr .eq. -7) then
        iesspin = .true.
    else
        iesspin = .false.
    end if

    if (LBox .gt. 0.d0) then
        PI_L = PI
        INV_PI_L = 1.d0/PI
    end if

    costk = 0.5d0
    do k = 1, npar
        econf(1, k) = 0.d0
    end do

    !          This is to avoid negative values
    if (iesd .ne. 6 .and. iesd .ne. 7) then
        if (npar .eq. 1) then
            vj(1) = abs(vj(1))
        elseif (npar .gt. 1) then
            vj(1) = abs(vj(1))
            vj(2) = abs(vj(2))
        end if
    end if

    do j = 1, nel
        zi = rcart(3, j)
        yi = rcart(2, j)
        xi = rcart(1, j)

        if (iesd .eq. -2) then
            if (LBox .gt. 0.d0) then
                write (*, *) ' Jastrow -2 not supported in Periodic Systems!'
                stop
            end if
            do i = j + 1, nel
                zj = rcart(3, i)
                yj = rcart(2, i)
                xj = rcart(1, i)

                dx = (xi - xj)**2
                dy = (yi - yj)**2
                dz = (zi - zj)**2
                nr = dsqrt(dx + dy + dz)
                nrz = dsqrt(vj(1)**2*(dx + dy) + vj(2)**2*dz)

                fcost = 1.d0/(2.d0*nrz*(1.d0 + nrz)**2)

                econf(1, 1) = econf(1, 1) - fcost*vj(1)*(dx + dy)*nr
                econf(1, 2) = econf(1, 2) - fcost*vj(2)*dz*nr

            end do
        else

            do i = j + 1, nel

                if (iesspin) then
                    if ((i .le. nelup .and. j .le. nelup) .or. (i .gt. nelup .and. j .gt. nelup)) then
                        ! parallel spins
                        costk = 0.25d0
                    else
                        costk = 0.5d0
                    end if
                end if

                zj = rcart(3, i)
                yj = rcart(2, i)
                xj = rcart(1, i)

                !*********** PERIODIC WORLD ***************
                if (LBox .le. 0.d0) then
                    nr = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2
                else
                    rc(1) = xi - xj
                    rc(2) = yi - yj
                    rc(3) = zi - zj
                    if (yes_tilted) then
                        call CartesianToCrystal(rc, 1)
                    end if
                    nr = INV_PI_L**2*sum((cellscale(1:3)                               &
                                         & *dsin(PI_L*rc(1:3)/cellscale(1:3)))**2)
                    !                 nr=0.d0
                    !                 nr=nr+dsin(PI_L*(xi-xj))**2
                    !                 nr=nr+dsin(PI_L*(yi-yj))**2
                    !                 nr=nr+dsin(PI_L*(zi-zj))**2
                    !                 nr=INV_PI_L**2*nr
                end if
                !*********** PERIODIC WORLD ***************

                if (npar .eq. 1) then

                    if (abs(iesd) .eq. 1) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            fcost = -costk*(cost/(1.d0 + vj(1)*cost))**2
                            econf(1, npar) = econf(1, npar) + fcost
                        end if

                    elseif (iesd .eq. 4) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            !                 Unstable for overflaw
                            !                  fcost=dexp(vj(1)*cost)
                            !                  fcost=costk*(1.d0-fcost+vj(1)*cost)/(vj(1)**2*fcost)
                            fcost = dexp(-vj(1)*cost)
                            fcost = costk*(fcost - 1.d0 + fcost*vj(1)*cost)/vj(1)**2
                            econf(1, npar) = econf(1, npar) + fcost
                        end if

                    elseif (iesd .eq. 6) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            if (vj(2) .gt. 1d-9) cost = (1.d0 - dexp(-vj(2)*cost))/vj(2)
                            fcost = -0.5d0*(cost/(1.d0 + vj(1)*cost))**2
                            econf(1, npar) = econf(1, npar) + fcost
                        end if

                    elseif (iesd .eq. 7) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            if (vj(2) .gt. 1d-9) cost = (1.d0 - dexp(-vj(2)*cost))/vj(2)
                            fcost = -costk*(cost/(1.d0 + vj(1)*cost))**2
                            econf(1, npar) = econf(1, npar) + fcost
                        end if

                    end if

                elseif (npar .eq. 2) then

                    if (iesd .eq. 2) then

                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            if (costk .eq. 0.5d0) then
                                fcost = -costk*(cost/(1.d0 + vj(1)*cost))**2
                                econf(1, 1) = econf(1, 1) + fcost
                            elseif (costk .eq. 0.25d0) then
                                fcost = -costk*(cost/(1.d0 + vj(2)*cost))**2
                                econf(1, 2) = econf(1, 2) + fcost
                            end if
                        end if

                    elseif (iesd .eq. 6) then

                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            if (vj(2) .gt. 1d-9) then
                                fcost2 = ((1 + vj(2)*cost)*dexp(-vj(2)*cost) - 1.d0)/vj(2)**2
                                cost = (1.d0 - dexp(-vj(2)*cost))/vj(2)
                            else
                                fcost2 = -0.5d0*cost**2
                            end if
                            fcost1 = -0.5d0*(cost/(1.d0 + vj(1)*cost))**2
                            fcost2 = fcost2*0.5d0/(1.d0 + vj(1)*cost)**2
                            econf(1, 1) = econf(1, 1) + fcost1
                            econf(1, 2) = econf(1, 2) + fcost2
                        end if
                    elseif (abs(iesd) .eq. 1) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            fcost = -costk*(cost/(1.d0 + vj(1)*cost))**2
                            econf(1, 1) = econf(1, 1) + fcost
                        end if
                    elseif (iesd .eq. 4) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            !                 Unstable for overflaw
                            !                  fcost=dexp(vj(1)*cost)
                            !                  fcost=costk*(1.d0-fcost+vj(1)*cost)/(vj(1)**2*fcost)
                            fcost = dexp(-vj(1)*cost)
                            fcost = costk*(fcost - 1.d0 + fcost*vj(1)*cost)/vj(1)**2
                            econf(1, 1) = econf(1, 1) + fcost
                        end if

                    elseif (iesd .eq. 7) then

                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)
                            if (vj(2) .gt. 1d-9) then
                                fcost2 = ((1 + vj(2)*cost)*dexp(-vj(2)*cost) - 1.d0)/vj(2)**2
                                cost = (1.d0 - dexp(-vj(2)*cost))/vj(2)
                            else
                                fcost2 = -0.5d0*cost**2
                            end if
                            fcost1 = -costk*(cost/(1.d0 + vj(1)*cost))**2
                            fcost2 = fcost2*costk/(1.d0 + vj(1)*cost)**2
                            econf(1, 1) = econf(1, 1) + fcost1
                            econf(1, 2) = econf(1, 2) + fcost2
                        end if

                    end if

                elseif (npar .eq. 3) then

                    if (nr .ne. 0.d0) then
                        cost = dsqrt(nr)

                        if (iesd .eq. 3) then
                            fcost1 = cost/(1.d0 + vj(1)*cost)
                            fcost2 = nr/(1.d0 + vj(2)*cost)**2
                            econf(1, 1) = econf(1, 1) - costk*fcost1**2
                            econf(1, 2) = econf(1, 2) - 2.d0*costk*vj(3)*fcost2**1.5d0
                            econf(1, 3) = econf(1, 3) + costk*fcost2

                        elseif (iesd .eq. 5) then
                            fat1 = dexp(-vj(1)*cost)
                            fat2 = dexp(-vj(2)*cost)
                            fat3 = vj(1) + vj(2)*vj(3)
                            fat4 = 0.5d0/fat3**2
                            fcost1 = -1.d0 - vj(3) + fat1*(1.d0 + fat3*cost) + vj(3)*fat2
                            fcost1 = fcost1*fat4
                            fcost2 = -1.d0 - vj(3) + fat1 + (vj(3) + fat3*cost)*fat2
                            fcost2 = fcost2*fat4*vj(3)
                            fcost3 = vj(1)*(1.d0 - fat2) + vj(2)*(-1.d0 + fat1)
                            fcost3 = fcost3*fat4

                            econf(1, 1) = econf(1, 1) + fcost1
                            econf(1, 2) = econf(1, 2) + fcost2
                            econf(1, 3) = econf(1, 3) + fcost3

                        end if

                    end if
                end if
            end do
        end if

        if (iesdr .le. -5) then
            do i = 1, nion
                zj = rion(3, i)
                yj = rion(2, i)
                xj = rion(1, i)

                !************** PERIODIC WORLD **************
                if (LBox .le. 0.d0) then
                    nr = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2
                else
                    rc(1) = xi - xj
                    rc(2) = yi - yj
                    rc(3) = zi - zj
                    if (yes_tilted) then
                        call CartesianToCrystal(rc, 1)
                    end if
                    nr = INV_PI_L**2*sum((cellscale(1:3)                       &
                                         & *dsin(PI_L*rc(1:3)/cellscale(1:3)))**2)
                    !                 nr=0.d0
                    !                 nr=nr+dsin(PI_L*(xi-xj))**2
                    !                 nr=nr+dsin(PI_L*(yi-yj))**2
                    !                 nr=nr+dsin(PI_L*(zi-zj))**2
                    !                 nr=nr*INV_PI_L**2
                end if
                !************** PERIODIC WORLD **************

                !              costz=(2.d0*zeta(i))**0.25d0
                !              costz3=costz**3

                if (npar .eq. 1) then

                    if (iesdr .eq. -5 .or. iesdr .eq. -17) then
                        !                    costk=0.5d0
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)*costz(i)
                            fcost = -0.5d0*(cost/(1.d0 + vj(1)*cost))**2
                            econf(1, npar) = econf(1, npar) - costz3(i)*fcost
                        end if
                    end if
                    if (iesdr .eq. -6 .or. iesdr .eq. -7 .or. iesdr .eq. -15) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)*costz(i)
                            !                 Unstable for overflaw
                            !                 fcost=dexp(vj(1)*cost)
                            !                 fcost=(1.d0-fcost+vj(1)*cost)/(2.d0*vj(1)**2*fcost)
                            fcost = dexp(-vj(1)*cost)
                            fcost = (fcost - 1.d0 + fcost*vj(1)*cost)/(2.d0*vj(1)**2)
                            econf(1, npar) = econf(1, npar) - costz3(i)*fcost
                        end if
                    end if
                elseif (npar .eq. 2) then
                    if (iesdr .eq. -5 .or. iesdr .eq. -17) then
                        !                    costk=0.5d0
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)*costz(i)
                            fcost = -0.5d0*(cost/(1.d0 + vj(2)*cost))**2
                            econf(1, 2) = econf(1, 2) - costz3(i)*fcost
                        end if
                    end if
                    if (iesdr .eq. -6 .or. iesdr .eq. -7 .or. iesdr .eq. -15) then
                        if (nr .ne. 0.d0) then
                            cost = dsqrt(nr)*costz(i)
                            !                 Unstable for overflaw
                            !                 fcost=dexp(vj(1)*cost)
                            !                 fcost=(1.d0-fcost+vj(1)*cost)/(2.d0*vj(1)**2*fcost)
                            fcost = dexp(-vj(2)*cost)
                            fcost = (fcost - 1.d0 + fcost*vj(2)*cost)/(2.d0*vj(2)**2)
                            econf(1, 2) = econf(1, 2) - costz3(i)*fcost
                        end if
                    end if
                end if
            end do
        end if
    end do

    return
end
