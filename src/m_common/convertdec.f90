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

!> @brief Convert decimal integer to ASCII digit array
!>
!> This subroutine converts a decimal integer to an array of ASCII digit
!> characters. It extracts each digit from the integer and converts it
!> to its corresponding ASCII representation.
!>
!> @param[in] imax Input integer to be converted (0 ≤ imax ≤ 999999)
!> @param[out] idigit Array of 6 ASCII digit characters (ASCII codes 48-57)
!>
!> @details
!> The subroutine performs the following operations:
!>
!> **Algorithm:**
!> 1. Extract the most significant digit (100000s place)
!> 2. Calculate remainder and extract next digit (10000s place)
!> 3. Continue for all 6 digits (1000s, 100s, 10s, 1s places)
!> 4. Convert each digit to ASCII by adding 48 (ASCII '0' = 48)
!>
!> **Digit Extraction:**
!> - idigit(1): 100000s place (imax/100000)
!> - idigit(2): 10000s place (remainder/10000)
!> - idigit(3): 1000s place (remainder/1000)
!> - idigit(4): 100s place (remainder/100)
!> - idigit(5): 10s place (remainder/10)
!> - idigit(6): 1s place (final remainder)
!>
!> **ASCII Conversion:**
!> Each digit is converted to ASCII by adding 48:
!> - 0 → 48 (ASCII '0')
!> - 1 → 49 (ASCII '1')
!> - ... → ...
!> - 9 → 57 (ASCII '9')
!>
!> @note The input integer should be in the range [0, 999999].
!> @note The output array contains 6 elements regardless of input value.
!> @note Leading zeros are preserved in the output array.
!> @note This subroutine is commonly used for file naming and string formatting.
subroutine convertdec(imax, idigit)
    implicit none
    integer imax, irest, i, idigit(6)
    
    !> @brief Extract 100000s digit
    idigit(1) = imax/100000
    irest = imax - idigit(1)*100000
    
    !> @brief Extract 10000s digit
    idigit(2) = irest/10000
    irest = irest - idigit(2)*10000
    
    !> @brief Extract 1000s digit
    idigit(3) = irest/1000
    irest = irest - idigit(3)*1000
    
    !> @brief Extract 100s digit
    idigit(4) = irest/100
    irest = irest - idigit(4)*100
    
    !> @brief Extract 10s digit
    idigit(5) = irest/10
    
    !> @brief Extract 1s digit
    idigit(6) = irest - 10*idigit(5)
    
    !> @brief Convert all digits to ASCII representation
    do i = 1, 6
        idigit(i) = idigit(i) + 48
    end do
    return
end
