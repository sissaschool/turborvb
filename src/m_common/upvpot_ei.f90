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

!> @brief Calculate electron-ion potential energy with Yukawa potential
!>
!> This subroutine computes the electron-ion potential energy for a given
!> electron position using a Yukawa potential with periodic boundary conditions.
!> The potential includes a screening parameter and can handle both isolated
!> and periodic systems.
!>
!> Parameters
!> ----------
!> RKEL : real*8 array, in
!>     Position of the electron (3 coordinates).
!> ZETA : real*8 array, in
!>     Array of ion charges (NION).
!> RION : real*8 array, in
!>     Positions of ions (3 × NION).
!> VPOT : real*8, out
!>     Computed electron-ion potential energy.
!> NION : integer, in
!>     Number of ions in the system.
!> LBOX : real*8, in
!>     Box length for periodic boundary conditions.
!> EPSVPOT : real*8, in
!>     Screening parameter for Yukawa potential.
!>     If EPSVPOT = 0: returns VPOT = 1.0 (no interaction).
!>     If EPSVPOT > 0: computes Yukawa potential.
!>
!> Notes
!> -----
!> - Uses Yukawa potential: V(r) = Z*exp(-κr)/r for each ion.
!> - Handles periodic boundary conditions through the NGIVEJ function.
!> - The total potential is the square root of the sum of individual contributions.
!> - For EPSVPOT = 0, returns a constant value (no electron-ion interaction).
!> - Optimized for quantum Monte Carlo applications.
!>
!> Algorithm
!> ---------
!> If EPSVPOT = 0:
!>   VPOT = 1.0
!> Else:
!>   VPOT = √(Σᵢ 2*zeta(i)/(rᵢ + epsvpot))
!>   where rᵢ is the distance between electron and ion i
!>
!> Example
!> -------
!> Used in quantum Monte Carlo for electron-ion interaction calculations.
subroutine upvpot_ei(rkel, zeta, rion, vpot, nion, LBox, epsvpot)
    implicit none
    integer nel, k, nion
    real(8) rkel(3)
    real(8) zeta(nion), vpot, ngivej, LBox, rion(3, *), cost, epsvpot
!> @brief Handle case with no electron-ion interaction
    if (epsvpot .eq. 0) then
        vpot = 1.d0
    else
!> @brief Initialize potential energy
        vpot = 0.d0
!> @brief Sum contributions from all ions using Yukawa potential
        do k = 1, nion
            cost = ngivej(rion(1, k), rkel, LBox) + epsvpot
            vpot = vpot + 2.d0*zeta(k)/cost
        end do
!> @brief Take square root of total potential
        vpot = dsqrt(vpot)
    end if
    return
end
