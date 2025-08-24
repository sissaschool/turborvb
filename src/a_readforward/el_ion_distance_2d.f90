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

!> @brief Computes electron-ion distances in 2D systems
!> @details This subroutine calculates the distances between electrons
!>          and ions in 2D systems. It projects the distances onto the
!>          z=0 plane and handles out-of-plane filtering. For periodic
!>          boundary conditions, it applies PBC to the distance vectors.
!>          For non-periodic systems, it filters out electrons that are
!>          too close to the z-plane.
!> @param[out] eliond Distance matrix between ions and electrons (nion, nel)
!> @param[in] kel Electron positions (3, nel)
!> @param[in] nel Number of electrons
!> @param[in] rion Ion positions (3, nion)
!> @param[in] nion Number of ions
!> @param[in] iespbc Flag for periodic boundary conditions
!> @param[in] outofplane Out-of-plane threshold distance
subroutine el_ion_distance_2d(eliond, kel, nel, rion, nion, iespbc, outofplane)

    use Cell
    implicit none

    integer nion, nel, k, j
    real(8) eliond(nion, *), rion(3, *), kel(3, *)
    real(8) dist_kel(3)
    real(8) outofplane
    logical iespbc

    if (iespbc) then !not implemented for PCB

        do j = 1, nel
            do k = 1, nion
                dist_kel(:) = kel(:, j) - rion(:, k)
                call ApplyPBC(dist_kel, 1)
                eliond(k, j) = dsqrt(max(sum(dist_kel(1:2)**2), 1d-18))
                !          if(eliond(k,j).lt.1d-18) eliond(k,j)=1d-18
            end do
        end do

    else

        do j = 1, nel
            do k = 1, nion
                if (abs(kel(3, j) - rion(3, k)) .ge. outofplane) then
                    !distance between k and j projected on the z=0 plane
                    eliond(k, j) = dsqrt(max((kel(1, j) - rion(1, k))**2 + (kel(2, j) - rion(2, k))**2, 1d-18))
                    !             if(eliond(k,j).lt.1d-18) eliond(k,j)=1d-18
                else
                    eliond(k, j) = 1.d8 ! put it out-of-range
                end if
            end do
        end do

    end if

end subroutine el_ion_distance_2d
