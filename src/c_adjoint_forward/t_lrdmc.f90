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

function t_lrdmc(kel, nion, rion, indvic, alat, ivic, plat, cellscale, iespbc)
    use cell, only: yes_tilted, s2r, car2cry
    use allio, only: zmin, zetar
    implicit none
    integer nion, nel, i
    integer indvic, imu, jmu
    real*8 t_lrdmc, alat, plat(*), ivic(3, *), distmu, distmin, rion(3, *), fun     &
            &, rmu(3), srmu(3), cellscale(3)
    real*8 kel(3)
    logical iespbc
    ! the kinetic part on the lattice is required (first 12 ivic)
    !      Compute the minimum distance between the electron
    !      and the nuclei's with maximum Z
    if (plat(1) .ne. 0.d0) then
        rmu(:) = kel(:) + 0.5d0*ivic(:, indvic) - rion(:, 1)
        if (iespbc) then
            if (yes_tilted) then
                !        call CartesianToCrystal(srmu,1)
                srmu(:) = car2cry(:, 1)*rmu(1) + car2cry(:, 2)*rmu(2) + car2cry(:, 3)*rmu(3)
                srmu(1) = anint(srmu(1)/cellscale(1))
                srmu(2) = anint(srmu(2)/cellscale(2))
                srmu(3) = anint(srmu(3)/cellscale(3))
                rmu(:) = rmu(:) - s2r(:, 1)*srmu(1) - s2r(:, 2)*srmu(2) - s2r(:, 3)*srmu(3)
                !        call dgemv('N',3,3,-1.d0,s2r,3,srmu,1,1.d0,rmu,1)
            else
                rmu(:) = rmu(:) - cellscale(:)*anint(rmu(:)/cellscale(:))
            end if
        end if
        distmin = 0.d0
        if (zetar(1) .ge. zmin) distmin = sum(rmu(:)**2)
        !     Calculate min el-ion  distance --> distmu
        do i = 2, nion
            if (zetar(i) .ge. zmin) then
                rmu(:) = kel(:) + 0.5d0*ivic(:, indvic) - rion(:, i)
                if (iespbc) then
                    if (yes_tilted) then
                        srmu(:) = car2cry(:, 1)*rmu(1) + car2cry(:, 2)*rmu(2) + car2cry(:, 3)*rmu(3)
                        srmu(1) = anint(srmu(1)/cellscale(1))
                        srmu(2) = anint(srmu(2)/cellscale(2))
                        srmu(3) = anint(srmu(3)/cellscale(3))
                        rmu(:) = rmu(:) - s2r(:, 1)*srmu(1) - s2r(:, 2)*srmu(2) - s2r(:, 3)*srmu(3)
                        !          srmu(:)=rmu(:)
                        !          call CartesianToCrystal(srmu,1)
                        !          srmu(1)=anint(srmu(1)/cellscale(1))
                        !          srmu(2)=anint(srmu(2)/cellscale(2))
                        !          srmu(3)=anint(srmu(3)/cellscale(3))
                        !          call dgemv('N',3,3,-1.d0,s2r,3,srmu,1,1.d0,rmu,1)
                    else
                        rmu(:) = rmu(:) - cellscale(:)*anint(rmu(:)/cellscale(:))
                    end if
                end if
                distmu = sum(rmu(:)**2)
                if (distmu .lt. distmin .or. distmin .eq. 0.d0) distmin = distmu
            end if
        end do
        if (plat(1) .gt. 0.d0) then
            fun = 1.d0/(1.d0 + plat(1)*distmin)
        else
            fun = exp(distmin/plat(1))
        end if
        if (indvic .le. 6) then
            t_lrdmc = plat(2)*fun/sum(ivic(:, indvic)**2)
        else
            t_lrdmc = plat(3)*(1.d0 - fun)/sum(ivic(:, indvic)**2)
        end if
    else
        if (indvic .le. 6) then
            t_lrdmc = plat(2)/sum(ivic(:, indvic)**2)
        else
            t_lrdmc = plat(3)/sum(ivic(:, indvic)**2)
        end if
    end if

    return
end
