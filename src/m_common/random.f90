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

!> @brief Random configuration selection for quantum Monte Carlo moves
!>
!> This subroutine implements a random configuration selection algorithm
!> for quantum Monte Carlo calculations, particularly useful for
!> non-diagonal moves in variational Monte Carlo and diffusion Monte Carlo.
!> It uses a cumulative probability table to select configurations
!> based on their transition probabilities.
!>
!> @param[in] L Number of particles or components per configuration
!> @param[in] Lz Total number of possible configurations
!> @param[in] ztry Random number in [diag, wtot] for configuration selection
!> @param[in] table Probability table containing transition weights (ipc, *)
!> @param[in] diag Diagonal element weight (lower bound for ztry)
!> @param[out] iout Selected particle index within the configuration
!> @param[out] indvic Selected configuration index
!> @param[out] sign Sign of the selected transition (+1 or -1)
!> @param[in] npow Power parameter for transition weight modification
!> @param[in] gamma Gamma parameter for transition weight modification
!> @param[in] istart Starting index for configuration selection
!> @param[in] indtm Array of maximum indices for each particle
!> @param[in] ipc Leading dimension of table array
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Initializes cumulative probability with diagonal element
!> 2. Iterates through configurations until cumulative probability ≥ ztry:
!>    - Calculates configuration and particle indices
!>    - Applies power and gamma modifications to transition weights
!>    - Accumulates modified weights in cumulative probability
!> 3. Handles edge cases and roundoff errors:
!>    - Adjusts index if final configuration has zero weight
!>    - Ensures valid configuration selection
!> 4. Determines sign based on transition weight sign
!> 5. Returns selected configuration and particle indices
!>
!> @note The algorithm assumes diag < ztry ≤ wtot (non-diagonal move)
!> @note Configuration index is calculated as: index = (i-1)/L + 1
!> @note Particle index is calculated as: ipart = i - (index-1)*L
!> @note Negative indvic indicates unrecoverable error (zero weight)
!> @note The subroutine handles both positive and negative transition weights
!> @note Power and gamma parameters modify transition weights for optimization
subroutine random(L, Lz, ztry, table, diag, iout, indvic, sign           &
        &, npow, gamma, istart, indtm, ipc)
    implicit none
    integer ipc, i, L, Lz, indvic, iout, index, ipart, istart, indtm(*)
    real*8 diag, table(ipc, *), try, ztry, sign, gamma, cost, npow

    !     choose a random configuration given the table,wtot
    !     it is assumed that  diag < ztry <= wtot ,i.e. non diagonal move

    try = diag
    i = 0
    do while (ztry .ge. try .and. i .lt. Lz)
        i = i + 1
        index = (i - 1)/L + 1
        ipart = i - (index - 1)*L
        if (npow .ne. 0.d0 .and. table(1, i) .gt. 0.d0 .and. &
            (index .ge. istart .and. index .le. indtm(ipart))) then
            try = try + table(1, i)*(1.d0 - npow*(1.d0 + gamma))
        elseif (table(1, i) .gt. 0.d0) then
            try = try + table(1, i)
        else
            cost = -gamma*table(1, i)
            try = try + cost
        end if
    end do

    !     It is possible by chance and roundoff that the final i=Lz and table(Lz)=0
    !     if for roundoff i=Lz set i the maximum with table(i)=/0
    if (i .eq. Lz) then
        do while ((table(1, i) .eq. 0.d0 .or. (table(1, i) .lt. 0 .and. gamma .eq. 0)&
                &.or. (npow*(1.d0 + gamma) .eq. 1.d0 .and. index .ge. istart)) .and. i .gt. 1)
            i = i - 1
            index = (i - 1)/L + 1
        end do
    end if

    !     if(i.eq.0) i=1
    if (table(1, i) .lt. 0.) then
        sign = -1.d0
    else
        sign = 1.d0
    end if

    indvic = index
    iout = ipart
    !      indvic=(i-1)/L+1
    !      iout=i-(indvic-1)*L

    if (table(1, i) .eq. 0.d0) indvic = -indvic !  Unrecoverable error

    return
end

