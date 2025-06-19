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

!> @brief Ion reference position optimization for periodic boundary conditions
!>
!> This module provides functionality for finding optimal reference positions
!> for ions in periodic systems, particularly useful for quantum Monte Carlo
!> calculations with periodic boundary conditions.
!>
!> The module includes:
!> - Subroutine for finding optimal reference positions that maximize
!>   minimum distances between ions
!> - Function for computing periodic modulo operations
!> - Algorithms for gap analysis in periodic coordinate systems
!>
!> @author TurboRVB group
!> @version 1.0
!> @date 2022

!> @brief Find optimal reference position for ions in periodic system
!>
!> This subroutine determines an optimal reference position that maximizes
!> the minimum distance between ions in a periodic system. It analyzes the
!> gaps between ion positions along a given direction and shifts the reference
!> to the center of the largest gap.
!>
!> @param[in] nion Number of ions in the system
!> @param[in] a Periodicity length along the direction of interest
!> @param[in] rion Array of ion positions (3, nion)
!> @param[in,out] ref Reference position (updated to optimal position)
!> @param[out] mindist Minimum distance between ions after optimization
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Computes periodic coordinates for all ions using setint function
!> 2. Sorts ions by their periodic coordinates
!> 3. Calculates gaps between consecutive ions in periodic space
!> 4. Handles the wrap-around gap between the last and first ion
!> 5. Finds the largest gap and its center position
!> 6. Shifts the reference position to the center of the largest gap
!> 7. Updates mindist to half the largest gap size
!>
!> @note The subroutine works in 1D along the x-direction of ion positions
!> @note For single ion systems, the gap is set to the full periodicity length
!> @note The optimal reference maximizes the minimum distance between any two ions
!> @note The subroutine modifies the input ref parameter to the optimal position
!> @note mindist represents the radius of the largest possible sphere around
!>       any ion that doesn't contain other ions
subroutine findrionfref(nion, a, rion, ref, mindist)
    implicit none
    integer i, nion
    real*8 rion(3, nion), ref, a, setint, shiftref, mindist
    real*8, dimension(:), allocatable :: gap, psip
    integer, dimension(:), allocatable :: ipsip

    allocate (gap(nion), ipsip(nion), psip(nion))

    !         write(6,*) ' Find modulo 0,a ',0,a
    do i = 1, nion
        psip(i) = setint(rion(1, i) + ref, a)
        !         write(6,*) i,psip(i)
    end do

    if (nion .eq. 1) then

        gap(1) = a
        ipsip(1) = 1

    else

        call dsortx(psip, 1, nion, ipsip)

        !         write(6,*) ' gap '
        do i = 1, nion - 1
            gap(ipsip(i)) = psip(i + 1) - psip(i)
            !         write(6,*) i,ipsip(i),gap(ipsip(i))
        end do
        gap(ipsip(nion)) = a - psip(nion) + psip(1)
        !         write(6,*) nion,ipsip(nion),gap(ipsip(nion))
        !         find the maximum gap
        !         recompute psip
        do i = 1, nion
            psip(i) = setint(rion(1, i) + ref, a)
        end do

        call dsortx(gap, 1, nion, ipsip)

        !         write(6,*) ' Ordered coordinate/gap '
        !         do i=1,nion
        !         write(6,*) i,ipsip(i),psip(ipsip(i)),gap(i)
        !         enddo

    end if

    shiftref = a - psip(ipsip(nion)) - gap(nion)/2.d0
    shiftref = setint(shiftref, a)
    if (abs(shiftref - a) .lt. shiftref) shiftref = shiftref - a

    mindist = gap(nion)/2.d0

    !         write(6,*) ' shiftref/mindist  inside =',shiftref,mindist

    ref = ref + shiftref

    deallocate (gap, psip, ipsip)

    return
end

!> @brief Compute periodic modulo operation
!>
!> This function computes the periodic modulo operation, mapping any real
!> number to the interval [0, a) where a is the periodicity length.
!> It handles both positive and negative input values correctly.
!>
!> @param[in] x Input value to be mapped to periodic interval
!> @param[in] a Periodicity length
!> @return Real value in the interval [0, a)
!>
!> @details
!> The function performs the following operations:
!> 1. For positive x: computes x mod a using integer division
!> 2. For negative x: computes (x + n*a + a) mod a where n = floor(-x/a)
!> 3. Ensures the result is always in the interval [0, a)
!>
!> @note The function is equivalent to fmod(x, a) but handles negative
!>       values correctly for periodic boundary conditions
!> @note For x ≥ 0: result = x - floor(x/a) * a
!> @note For x < 0: result = x + floor(-x/a) * a + a
!> @note The result is always non-negative and less than a
function setint(x, a)
    real*8 setint, x, a
    integer n
    if (x .ge. 0.d0) then
        n = x/a
        setint = x - n*a
    else
        n = -x/a
        setint = x + n*a + a
    end if
    return
end

