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

!> @brief      Scale electron and ion positions by cell scale factors
!> @details    Multiplies electron positions (kel) and ion positions (rion) by
!>             the cell scale factors (scalecell). This is used for cell scaling
!>             operations in periodic boundary conditions or unit cell transformations.
!> @param[in]  nion        Number of ions
!> @param[in]  nel         Number of electrons
!> @param[in]  indt        Index offset for derivatives
!> @param[in]  scalecell   Cell scale factors (3 components)
!> @param[in,out] kel      Electron positions array
!> @param[in,out] rion     Ion positions array
subroutine upcellkel(nion, nel, indt, scalecell, kel, rion)
    implicit none
    integer nion, nel, i, indt
    real*8 scalecell(3), kel(3, *), rion(3, nion)
    do i = 1, nel*(indt + 1)
        kel(:, i) = scalecell(:)*kel(:, i)
    end do
    do i = 1, nion
        rion(:, i) = rion(:, i)*scalecell(:)
    end do
    return
end
