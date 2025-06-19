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

!> @brief Jastrow factor type mapping function for one-body terms
!>
!> This function maps Jastrow factor type codes (iesdr) to specific one-body
!> Jastrow factor implementations used in quantum Monte Carlo calculations.
!> It provides a standardized way to determine the appropriate one-body
!> Jastrow factor type based on the input parameter.
!>
!> @param[in] iesdr Jastrow factor type code
!> @return Integer code representing the specific one-body Jastrow factor type
!>
!> @details
!> The function performs the following mapping:
!> - iesdr = 0: No one-body Jastrow factor (return 0)
!> - iesdr = -1: Type 1 one-body Jastrow factor
!> - iesdr = -5: Type 1 one-body Jastrow factor
!> - iesdr = -15: Type 4 one-body Jastrow factor (default short range)
!> - iesdr = -6: Type 4 one-body Jastrow factor
!> - iesdr = -7: Type 4 one-body Jastrow factor
!> - iesdr = -17: Type 1 one-body Jastrow factor
!> - iesdr = -20: Type 4 one-body Jastrow factor
!> - iesdr = -21: Type 4 one-body Jastrow factor
!> - Default: Type 4 one-body Jastrow factor (short range exponential)
!>
!> @note The function provides backward compatibility for different
!>       Jastrow factor implementations in TurboRVB
!> @note Type 4 corresponds to the default short range one-body exponential
!>       Jastrow factor
!> @note Type 1 corresponds to alternative one-body Jastrow factor forms
!> @note The function is used in wave function setup and optimization
function iesdr1iesd(iesdr)
    implicit none
    integer iesdr, iesdr1iesd
    iesdr1iesd = 4 ! default short range one-body
    !     one body exp.
    if (iesdr .eq. 0) iesdr1iesd = 0 ! No one body Jastrow
    if (iesdr .eq. -1) iesdr1iesd = 1
    if (iesdr .eq. -15) iesdr1iesd = 4
    if (iesdr .eq. -5) iesdr1iesd = 1
    if (iesdr .eq. -6) iesdr1iesd = 4
    if (iesdr .eq. -7) iesdr1iesd = 4
    if (iesdr .eq. -17) iesdr1iesd = 1
    if (iesdr .eq. -20) iesdr1iesd = 4
    if (iesdr .eq. -21) iesdr1iesd = 4
    return
end
