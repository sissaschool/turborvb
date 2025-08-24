! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2002-2004 quantum-ESPRESSO group
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

!
!------------------------------------------------------------------------------!
!> @file kind.f90
!> @brief Data type definitions and precision information module
!> @details This module defines portable data types for TurboRVB, ensuring
!>          consistent precision across different platforms and compilers.
!>          Based on quantum-ESPRESSO's kind module, it provides high-precision
!>          real types and utility functions for data type information.
!> @author TurboRVB group (based on quantum-ESPRESSO group code)
!> @date 2022
!> @see print_kind_info
module kinds
    !------------------------------------------------------------------------------!

    implicit none
    save
    ! ... kind definitions
    !> @brief Double precision real type (14 digits, range up to 10^200)
    integer, parameter :: DP = selected_real_kind(14, 200)
    !> @brief Single precision real type (6 digits, range up to 10^30)
    integer, parameter :: sgl = selected_real_kind(6, 30)
    !> @brief 32-bit integer type (range up to 10^9)
    integer, parameter :: i4b = selected_int_kind(9)
    private
    public :: i4b, sgl, DP, print_kind_info
    !
    !------------------------------------------------------------------------------!
    !
contains
    !
    !------------------------------------------------------------------------------!
    !
    !> @brief Print detailed information about the used data types
    !> @details This subroutine provides comprehensive information about the
    !>          precision, range, and characteristics of the defined data types.
    !>          Useful for debugging and ensuring numerical accuracy across
    !>          different platforms and compilers.
    !> @param[in] stdout Integer unit number for output (typically 6 for stdout)
    !> @note The subroutine prints information for:
    !>       - DP (double precision): precision, epsilon, tiny, huge values
    !>       - sgl (single precision): precision, epsilon, tiny, huge values
    !>       - INTEGER: kind value, bit size, largest representable number
    !>       - LOGICAL: kind value
    !>       - CHARACTER: kind value
    !> @see DP, sgl, i4b
    subroutine print_kind_info(stdout)
        !
        !------------------------------------------------------------------------------!
        !
        implicit none
        integer, intent(IN) :: stdout
        !
        write (stdout, '(/,T2,A)') 'DATA TYPE INFORMATION:'
        !
        write (stdout, '(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E16.8))') &
            'REAL: Data type name:', 'DP', '      Kind value:', kind(0.0_dp), &
            '      Precision:', precision(0.0_dp), &
            '      Smallest nonnegligible quantity relative to 1:', &
            epsilon(0.0_dp), '      Smallest positive number:', tiny(0.0_dp), &
            '      Largest representable number:', huge(0.0_dp)
        write (stdout, '(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E16.8))') &
            '      Data type name:', 'sgl', '      Kind value:', kind(0.0_sgl), &
            '      Precision:', precision(0.0_sgl), &
            '      Smallest nonnegligible quantity relative to 1:', &
            epsilon(0.0_sgl), '      Smallest positive number:', tiny(0.0_sgl), &
            '      Largest representable number:', huge(0.0_sgl)
        write (stdout, '(/,T2,A,T72,A,4(/,T2,A,T61,I20))') &
            'INTEGER: Data type name:', '(default)', '         Kind value:', &
            kind(0), '         Bit size:', bit_size(0), &
            '         Largest representable number:', huge(0)
        write (stdout, '(/,T2,A,T72,A,/,T2,A,T75,I6,/)') 'LOGICAL: Data type name:', &
            '(default)', '         Kind value:', kind(.true.)
        write (stdout, '(/,T2,A,T72,A,/,T2,A,T75,I6,/)') &
            'CHARACTER: Data type name:', '(default)', '           Kind value:', &
            kind('C')
        !
    end subroutine print_kind_info
    !
    !------------------------------------------------------------------------------!
end module kinds
!------------------------------------------------------------------------------!
