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

!> @brief Calculate ion-ion potential energy with Ewald summation
!>
!> This function computes the ion-ion electrostatic potential energy
!> using either direct Coulomb interaction (for isolated systems) or
!> Ewald summation (for periodic boundary conditions). The calculation
!> includes both the real-space and self-interaction contributions.
!>
!> Parameters
!> ----------
!> ZETA : real*8 array, in
!>     Array of ion charges (NION).
!> IOND : real*8 array, in
!>     Matrix of ion-ion distances (NION × NION).
!> NION : integer, in
!>     Number of ions in the system.
!> KAPPA : real*8, in
!>     Ewald parameter for splitting the Coulomb interaction.
!> LBOX : real*8, in
!>     Box length for periodic boundary conditions.
!>     If LBOX ≤ 0: isolated system (direct Coulomb).
!>     If LBOX > 0: periodic system (Ewald summation).
!>
!> Returns
!> -------
!> real*8
!>     Total ion-ion potential energy including self-interaction correction.
!>
!> Notes
!> -----
!> - For isolated systems (LBOX ≤ 0): uses direct Coulomb interaction.
!> - For periodic systems (LBOX > 0): uses Ewald summation with:
!>   * Real-space contribution using complementary error function.
!>   * Self-interaction correction to avoid double-counting.
!> - Includes nearest-neighbor interactions for periodic systems.
!> - Uses SIMD vectorization for performance optimization.
!> - The self-interaction term is: -2*κ*ε₀/√π * Σ(zeta²).
!>
!> Algorithm
!> ---------
!> For isolated systems:
!>   E = Σᵢⱼ 2*zeta(i)*zeta(j)*veps(iond(i,j))
!>
!> For periodic systems:
!>   E_real = Σᵢⱼ 2*zeta(i)*zeta(j)*erfc(κ*rᵢⱼ)/rᵢⱼ
!>   E_self = -2*κ*ε₀/√π * Σ(zeta²)
!>   E_total = E_real + E_self
!>
!> Example
!> -------
!> Used in quantum Monte Carlo calculations for ion-ion interactions.
function upvpotaa(zeta, iond, nion, kappa, LBox)
    use dielectric
    use allio, only: iond_cart, x_neigh, neigh, dist_shift, rank
    implicit none
    integer nel, i, j, ii, jj, nion
    real(8) :: derfc, x_shift(3), cost_z
    real(8) zeta(nion), iond(nion, nion), kappa, LBox, pot_aa, upvpotaa, eself1b
    double precision, parameter :: PI = 3.14159265358979323846d0

!> @brief Initialize potential energy
    pot_aa = 0.d0

    if (LBox .le. 0.d0) then
!> @brief Isolated system: direct Coulomb interaction
        do i = 1, nion
            do j = i + 1, nion
                if (zeta(i)*zeta(j) .ne. 0.d0) &
       &pot_aa = pot_aa + 2.d0*zeta(i)*zeta(j)*veps(iond(i, j))
            end do
        end do
!> @brief No self-interaction correction for isolated systems
        eself1b = 0.d0
    else
!> @brief Periodic system: Ewald summation
        do i = 1, nion
            do j = i + 1, nion
                cost_z = 2.d0*zeta(i)*zeta(j)
                jj = nion*(j - 1) + i
                if (cost_z .ne. 0.d0) then
!> @brief Calculate real-space contribution with nearest neighbors
#ifdef _SIMD
!$omp simd
#endif
                    do ii = 1, neigh
                        dist_shift(ii) = dsqrt((iond_cart(1, jj)&
                           & + x_neigh(ii, 1))**2 + (iond_cart(2, jj)&
                           & + x_neigh(ii, 2))**2 + (iond_cart(3, jj)&
                           & + x_neigh(ii, 3))**2)
                        dist_shift(ii) = cost_z*rep_erfc(dist_shift(ii), kappa)
                    end do
                    pot_aa = pot_aa + sum(dist_shift(1:neigh))
                end if
            end do
        end do
!> @brief Self-interaction correction for periodic systems
        eself1b = -2*kappa*epsilon0/dsqrt(pi)*sum(zeta(1:nion)**2)

    end if
!> @brief Return total potential energy including self-interaction
    upvpotaa = pot_aa + eself1b

    return
end
