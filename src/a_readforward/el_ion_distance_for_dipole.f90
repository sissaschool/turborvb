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

subroutine el_ion_distance_for_dipole(eliond, kel, nel, rion, nion, iespbc)

    use Cell
    implicit none

    integer nion, nel, k, j
    real(8) eliond(nion, *), rion(3, *), kel(3, *)
    real(8) dist_kel(3)
    logical iespbc

    if (iespbc) then !Calcola la distanza lungo Z in pbc taking in to account also the sign!

        do j = 1, nel
            do k = 1, nion
                dist_kel(:) = kel(:, j) - rion(:, k)
                call ApplyPBC(dist_kel, 1)
                eliond(k, j) = dist_kel(3) !Modified here from el_ion_distance_2D
                if (abs(eliond(k, j)) .lt. 1d-9) eliond(k, j) = 1d-9
            end do
        end do

    else

        do j = 1, nel
            do k = 1, nion
                eliond(k, j) = kel(3, j) - rion(3, k) !distance along z
                if (abs(eliond(k, j)) .lt. 1d-9) eliond(k, j) = 1d-9
            end do
        end do

    end if

end subroutine el_ion_distance_for_dipole
