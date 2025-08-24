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

!> @brief Calculate damping function for quantum chemistry interactions
!>
!> This function computes a damping factor used in quantum chemistry
!> calculations to modify the strength of interactions based on distance
!> and quantum number. The damping function is based on the work by
!> G. Scoles et al. in J. Chem. Phys. 76, 3057 (1982).
!>
!> @param[in] r Distance parameter (typically interatomic distance)
!> @param[in] n Quantum number parameter (0 for no damping)
!> @return damping Damping factor (0 ≤ damping ≤ 1)
!>
!> @details
!> The damping function implements the following formula:
!>
!> **For n ≠ 0:**
!> ```
!> damping = (1 - exp(-r*2.1/n - 0.109*r²/√n))^n
!> ```
!>
!> **For n = 0:**
!> ```
!> damping = 1.0
!> ```
!>
!> **Physical Interpretation:**
!> - The damping function reduces the strength of interactions at short distances
!> - The reduction is more pronounced for larger quantum numbers (n)
!> - The function approaches 1.0 as distance increases
!> - For n = 0, no damping is applied (damping = 1.0)
!>
!> **Parameters:**
!> - r: Distance parameter (typically in atomic units)
!> - n: Quantum number that controls the damping strength
!> - 2.1: Linear damping coefficient
!> - 0.109: Quadratic damping coefficient
!>
!> @note The function is commonly used in van der Waals and dispersion interactions.
!> @note The damping factor is always in the range [0, 1].
!> @note For n = 0, the function returns 1.0 regardless of the distance r.
!> @note The exponential form ensures smooth behavior across all distances.
function damping(r, n)
    implicit none
    real*8 damping, r
    integer n
    !     damping=1.d0
    !     return
    !     From G. Scoles et al. J. Chem. Phys. 76, 3057 (1982);
    !> @brief Return unity damping for n = 0 (no damping)
    if (n .ne. 0) then
        !> @brief Calculate damping factor using Scoles formula
        !> @note Based on G. Scoles et al. J. Chem. Phys. 76, 3057 (1982)
        damping = (1.d0 - dexp(-r*2.1d0/n - 0.109d0*r**2/dsqrt(dble(n))))**n
    else
        !> @brief No damping applied for n = 0
        damping = 1.d0
    end if
    return
end
