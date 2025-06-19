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

!> @brief Portable random number generator for quantum Monte Carlo
!>
!> This function implements a portable random number generator based on
!> the Park-Miller minimal standard multiplicative linear congruential
!> generator combined with a Marsaglia shift sequence. It provides
!> uniform random numbers in the interval (0,1) for use in quantum
!> Monte Carlo calculations.
!>
!> @param[in,out] idum Seed for random number generation:
!>                      - idum ≤ 0: Initialize generator with |idum|
!>                      - idum > 0: Generate next random number
!> @return Real random number in the interval (0,1)
!>
!> @details
!> The function implements a combined generator with the following features:
!> 1. Park-Miller multiplicative linear congruential generator:
!>    - Parameters: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
!>    - Recurrence: iy = IA * (iy - k*IQ) - IR*k (mod IM)
!> 2. Marsaglia shift sequence for ix:
!>    - Uses XOR and bit shifts for additional randomness
!> 3. Combined output: ran = am * (ix XOR iy OR 1)
!>
!> @note The generator has a period of approximately 2^31 - 1
!> @note The function uses saved variables to maintain state between calls
!> @note Initialization occurs when idum ≤ 0 or iy < 0
!> @note The function is portable across different platforms and compilers
!> @note Used extensively in Monte Carlo sampling and random walk algorithms
function ran(idum)
    implicit none
    integer, parameter :: K4B = selected_int_kind(9)
    integer(K4B), intent(INOUT) :: idum
    real :: ran
    integer(K4B), parameter :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2
    real, save :: am
    integer(K4B), save :: ix = -1, iy = -1, k
    if (idum <= 0 .or. iy < 0) then
        am = nearest(1.0, -1.0)/IM
        iy = ior(ieor(888889999, abs(idum)), 1)
        ix = ieor(777755555, abs(idum))
        idum = abs(idum) + 1
    end if
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    k = iy/IQ
    iy = IA*(iy - k*IQ) - IR*k
    if (iy < 0) iy = iy + IM
    ran = am*ior(iand(IM, ieor(ix, iy)), 1)
end
