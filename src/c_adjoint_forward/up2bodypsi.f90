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

subroutine up2bodypsi(nel, nelup, psiln, kel, vj, iesdr          &
        &, rion, nion, costz, costz3, LBox, niesd)
    use Constants
    use Cell
    use allio, only: pointvj, n_body_on, norm_metric
    implicit none

    integer nel, nelup, iesd, iesdr, j1, j2, ispin, nion, iesdr2iesd    &
            &, niesd, k
    real*8 kel(3, nel), psiln
    real*8 vj(max(1, niesd)), rc(3), rion(3, nion), costz(nion), costz3(nion)&
            &, costz0, jastrow_ei, jastrow_ee, r0
    logical iesspin

    real(8) LBox

    iesd = iesdr2iesd(iesdr)

    if (iesdr .eq. -7 .or. iesd .eq. 2 .or. iesd .lt. 0) then
        iesspin = .true.
    else
        iesspin = .false.
    end if

!$omp parallel do default(shared) reduction(+: psiln) private(costz0,j1,j2,ispin,rc,r0)
    do j1 = 1, nel
        costz0 = 0.d0
        !do for the second particle
        do j2 = j1 + 1, nel
            if (iesspin) then
                if ((j1 .le. nelup .and. j2 .le. nelup) .or. (j1 .gt. nelup .and. j2 .gt. nelup)) then
                    ! parallel spins
                    ispin = 1
                else
                    ispin = -1
                end if
            else
                ispin = -1
            end if
            rc(:) = kel(:, j1) - kel(:, j2)
            !************ PERIODIC WORLD ***********
            if (LBox .gt. 0.d0) then
                !       call CartesianToCrystal(rc,1)
                rc(:) = car2cry(:, 1)*rc(1) + car2cry(:, 2)*rc(2) + car2cry(:, 3)*rc(3)
                !       rc(1:3)=cellpi(1:3)*dsin(rc(1:3)/cellpi(1:3))
                do k = 1, 3
                    rc(k) = map(rc(k), cellscale(k))
                end do
            end if
            !************ PERIODIC WORLD ***********

            costz0 = costz0 + jastrow_ee(rc, vj, iesd, ispin)
        end do

        !        psiln=costz0+psiln

        if (n_body_on .ne. 0) then

            !        costz0=0.d0

            do j2 = 1, nion

                !           costz=(2.d0*zeta(j2))**0.25d0
                !           costz3=costz**3

                if (LBox .gt. 0.d0) then
                    rc(:) = kel(:, j1) - rion(:, j2)
                    !         call CartesianToCrystal(rc,1)
                    rc(:) = car2cry(:, 1)*rc(1) + car2cry(:, 2)*rc(2) + car2cry(:, 3)*rc(3)
                    !         rc(1:3)=costz(j2)*cellpi(1:3)*dsin(rc(1:3)/cellpi(1:3))
                    do k = 1, 3
                        rc(k) = costz(j2)*map(rc(k), cellscale(k))
                    end do
                    r0 = norm_metric(rc, metric)
                else
                    rc(:) = costz(j2)*(kel(:, j1) - rion(:, j2))
                    r0 = dsqrt(sum(rc(:)**2))
                end if

                costz0 = costz0 - jastrow_ei(r0, vj(pointvj(1, j2)), pointvj(2, j2))*costz3(j2)

            end do
        end if
        psiln = costz0 + psiln

    end do

    return
end
