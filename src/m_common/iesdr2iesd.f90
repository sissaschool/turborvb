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

!> @brief Jastrow factor type mapping function for two-body terms
!>
!> This function maps Jastrow factor type codes (iesdr) to specific two-body
!> Jastrow factor implementations used in quantum Monte Carlo calculations.
!> It provides a standardized way to determine the appropriate two-body
!> Jastrow factor type based on the input parameter.
!>
!> @param[in] iesdr Jastrow factor type code
!> @return Integer code representing the specific two-body Jastrow factor type
!>
!> @details
!> The function performs the following mapping:
!> - Default: Returns the input iesdr value unchanged
!> - iesdr = -15: Maps to Type 1 two-body Jastrow factor
!> - iesdr = -5: Maps to Type 1 two-body Jastrow factor
!> - iesdr = -6: Maps to Type 4 two-body Jastrow factor
!> - iesdr = -7: Maps to Type 4 two-body Jastrow factor
!> - iesdr = -17: Maps to Type -1 two-body Jastrow factor
!> - iesdr = -20: Maps to Type -1 two-body Jastrow factor
!> - iesdr = -21: Maps to Type 2 two-body Jastrow factor
!>
!> @note The function provides backward compatibility for different
!>       Jastrow factor implementations in TurboRVB
!> @note Type 1 corresponds to standard two-body Jastrow factors
!> @note Type 4 corresponds to extended two-body Jastrow factors
!> @note Type -1 corresponds to alternative two-body Jastrow factor forms
!> @note Type 2 corresponds to specialized two-body Jastrow factor forms
!> @note The function is used in wave function setup and optimization
function iesdr2iesd(iesdr)
    implicit none
    integer iesdr, iesdr2iesd
    iesdr2iesd = iesdr
    if (iesdr .eq. -15) iesdr2iesd = 1
    if (iesdr .eq. -5) iesdr2iesd = 1
    if (iesdr .eq. -6) iesdr2iesd = 4
    if (iesdr .eq. -7) iesdr2iesd = 4
    if (iesdr .eq. -17) iesdr2iesd = -1
    if (iesdr .eq. -20) iesdr2iesd = -1
    if (iesdr .eq. -21) iesdr2iesd = 2
    return
end
