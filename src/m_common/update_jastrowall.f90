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

!> @brief Update Jastrow factor gradient and Laplacian components
!>
!> This subroutine updates the gradient and Laplacian components of the
!> Jastrow factor matrix by adding contributions from the PSIP array.
!> Used for unrestricted wave functions (no spin separation).
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> JASTROWALL : real*8 array, inout
!>     Jastrow factor matrix (NEL × NEL).
!>     On entry: current Jastrow matrix.
!>     On exit: updated Jastrow matrix.
!> PSIP : real*8 array, in
!>     Update contributions (NEL × NEL).
!>
!> Notes
!> -----
!> - Updates off-diagonal elements only (k ≠ j).
!> - Used for gradient and Laplacian calculations in quantum Monte Carlo.
!> - Applies to unrestricted wave functions.
!>
!> Algorithm
!> ---------
!> For each electron pair (j,k) with j ≠ k:
!>   JASTROWALL(k,j) += PSIP(j,k)
!>
!> Example
!> -------
!> Used in variational Monte Carlo for Jastrow factor updates.
subroutine upjastrowall(nel, jastrowall, psip)
    implicit none
    integer j, k, nel
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nel
        do k = 1, nel
            if (k .ne. j) jastrowall(k, j) = jastrowall(k, j) + psip(j, k)
        end do
    end do
    return
end subroutine upjastrowall

!> @brief Update Jastrow factor for Sz-conserving wave functions
!>
!> This subroutine updates the Jastrow factor matrix for Sz-conserving
!> wave functions, handling spin-up and spin-down electrons separately
!> with appropriate sign changes for opposite spins.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Total number of electrons.
!> NELUP : integer, in
!>     Number of spin-up electrons.
!> JASTROWALL : real*8 array, inout
!>     Jastrow factor matrix (NEL × NEL).
!>     On entry: current Jastrow matrix.
!>     On exit: updated Jastrow matrix.
!> PSIP : real*8 array, in
!>     Update contributions (NEL × NEL).
!>
!> Notes
!> -----
!> - Handles spin-up and spin-down electrons separately.
!> - Same-spin interactions: add contribution.
!> - Opposite-spin interactions: subtract contribution.
!> - Used for Sz-conserving wave functions.
!>
!> Algorithm
!> ---------
!> For spin-up electrons (j ≤ NELUP):
!>   Same spin (k ≤ NELUP): JASTROWALL(k,j) += PSIP(j,k)
!>   Opposite spin (k > NELUP): JASTROWALL(k,j) -= PSIP(j,k)
!> For spin-down electrons (j > NELUP):
!>   Opposite spin (k ≤ NELUP): JASTROWALL(k,j) -= PSIP(j,k)
!>   Same spin (k > NELUP): JASTROWALL(k,j) += PSIP(j,k)
!>
!> Example
!> -------
!> Used in variational Monte Carlo for Sz-conserving Jastrow updates.
subroutine upjastrowall_sz(nel, nelup, jastrowall, psip)
    implicit none
    integer j, k, nel, nelup
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nelup
        do k = 1, nelup
            if (k .ne. j) jastrowall(k, j) = jastrowall(k, j) + psip(j, k)
        end do
        do k = nelup + 1, nel
            jastrowall(k, j) = jastrowall(k, j) - psip(j, k)
        end do
    end do
    do j = nelup + 1, nel
        do k = 1, nelup
            jastrowall(k, j) = jastrowall(k, j) - psip(j, k)
        end do
        do k = nelup + 1, nel
            if (k .ne. j) jastrowall(k, j) = jastrowall(k, j) + psip(j, k)
        end do
    end do
    return
end subroutine upjastrowall_sz

!> @brief Update Jastrow factor for pseudo and lattice positions in LRDMC
!>
!> This subroutine updates the (1:INDT) part of the Jastrow factor matrix
!> for pseudo and lattice positions in LRDMC calculations. Used for
!> unrestricted wave functions.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> JASTROWALL : real*8 array, inout
!>     Jastrow factor matrix (NEL × NEL).
!>     On entry: current Jastrow matrix.
!>     On exit: updated Jastrow matrix.
!> PSIPMU : real*8 array, in
!>     Update contributions for pseudo/lattice positions (NEL × NEL).
!>
!> Notes
!> -----
!> - Updates off-diagonal elements only (j ≠ i).
!> - Specifically for pseudo and lattice position updates in LRDMC.
!> - Used for unrestricted wave functions.
!>
!> Algorithm
!> ---------
!> For each electron pair (i,j) with j ≠ i:
!>   JASTROWALL(j,i) += PSIPMU(i,j)
!>
!> Example
!> -------
!> Used in LRDMC for pseudo and lattice position updates.
subroutine upjastrowallfat(nel, jastrowall, psipmu)
    implicit none
    integer nel, i, j
    real*8 jastrowall(nel, *), psipmu(nel, *)

    do i = 1, nel
        do j = 1, nel
            if (j .ne. i) jastrowall(j, i) = jastrowall(j, i) + psipmu(i, j)
        end do
    end do

    return
end subroutine upjastrowallfat

!> @brief Update Jastrow factor for Sz-conserving pseudo/lattice positions
!>
!> This subroutine updates the (1:INDT) part of the Jastrow factor matrix
!> for pseudo and lattice positions in LRDMC calculations, handling
!> Sz-conserving wave functions with proper spin separation.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Total number of electrons.
!> NELUP : integer, in
!>     Number of spin-up electrons.
!> JASTROWALL : real*8 array, inout
!>     Jastrow factor matrix (NEL × NEL).
!>     On entry: current Jastrow matrix.
!>     On exit: updated Jastrow matrix.
!> PSIPMU : real*8 array, in
!>     Update contributions for pseudo/lattice positions (NEL × NEL).
!>
!> Notes
!> -----
!> - Handles spin-up and spin-down electrons separately.
!> - Same-spin interactions: add contribution.
!> - Opposite-spin interactions: subtract contribution.
!> - Specifically for pseudo and lattice position updates in LRDMC.
!>
!> Algorithm
!> ---------
!> For spin-up electrons (i ≤ NELUP):
!>   Same spin (j ≤ NELUP): JASTROWALL(j,i) += PSIPMU(i,j)
!>   Opposite spin (j > NELUP): JASTROWALL(j,i) -= PSIPMU(i,j)
!> For spin-down electrons (i > NELUP):
!>   Opposite spin (j ≤ NELUP): JASTROWALL(j,i) -= PSIPMU(i,j)
!>   Same spin (j > NELUP): JASTROWALL(j,i) += PSIPMU(i,j)
!>
!> Example
!> -------
!> Used in LRDMC for Sz-conserving pseudo and lattice position updates.
subroutine upjastrowallfat_sz(nel, nelup, jastrowall, psipmu)
    implicit none
    integer nel, nelup, i, j
    real*8 jastrowall(nel, *), psipmu(nel, *)

    do i = 1, nelup
        do j = 1, nelup
            if (j .ne. i) jastrowall(j, i) = jastrowall(j, i) + psipmu(i, j)
        end do
        do j = nelup + 1, nel
            jastrowall(j, i) = jastrowall(j, i) - psipmu(i, j)
        end do
    end do

    do i = nelup + 1, nel
        do j = 1, nelup
            jastrowall(j, i) = jastrowall(j, i) - psipmu(i, j)
        end do
        do j = nelup + 1, nel
            if (j .ne. i) jastrowall(j, i) = jastrowall(j, i) + psipmu(i, j)
        end do
    end do

    return
end subroutine upjastrowallfat_sz

!> @brief Update Jastrow factor wave function components
!>
!> This subroutine updates the wave function-related components of the
!> Jastrow factor matrix by adding contributions from the PSIP array.
!> Used for unrestricted wave functions.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> JASTROWALL : real*8 array, inout
!>     Jastrow factor matrix (NEL × NEL).
!>     On entry: current Jastrow matrix.
!>     On exit: updated Jastrow matrix.
!> PSIP : real*8 array, in
!>     Update contributions for wave function components (NEL × NEL).
!>
!> Notes
!> -----
!> - Updates off-diagonal elements only (j ≠ k).
!> - Specifically for wave function component updates.
!> - Used for unrestricted wave functions.
!>
!> Algorithm
!> ---------
!> For each electron pair (j,k) with j ≠ k:
!>   JASTROWALL(j,k) += PSIP(k,j)
!>
!> Example
!> -------
!> Used in variational Monte Carlo for wave function component updates.
subroutine upjastrowallpsi(nel, jastrowall, psip)
    implicit none
    integer j, k, nel
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nel
        do k = 1, nel
            if (j .ne. k) jastrowall(j, k) = jastrowall(j, k) + psip(k, j)
        end do
    end do

    return
end subroutine upjastrowallpsi

!> @brief Update Jastrow factor wave function components for Sz-conserving case
!>
!> This subroutine updates the wave function-related components of the
!> Jastrow factor matrix for Sz-conserving wave functions, handling
!> spin-up and spin-down electrons with proper sign changes and symmetry.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Total number of electrons.
!> NELUP : integer, in
!>     Number of spin-up electrons.
!> JASTROWALL : real*8 array, inout
!>     Jastrow factor matrix (NEL × NEL).
!>     On entry: current Jastrow matrix.
!>     On exit: updated Jastrow matrix.
!> PSIP : real*8 array, in
!>     Update contributions for wave function components (NEL × NEL).
!>
!> Notes
!> -----
!> - Handles spin-up and spin-down electrons separately.
!> - Same-spin interactions: add contribution.
!> - Opposite-spin interactions: subtract contribution.
!> - Ensures matrix symmetry by copying upper triangle to lower triangle.
!> - Used for Sz-conserving wave functions.
!>
!> Algorithm
!> ---------
!> 1. Update spin-up electrons (j ≤ NELUP):
!>    Same spin (k > j): JASTROWALL(k,j) += PSIP(k,j)
!>    Opposite spin (k > NELUP): JASTROWALL(k,j) -= PSIP(k,j)
!> 2. Update spin-down electrons (j > NELUP):
!>    Same spin (k > j): JASTROWALL(k,j) += PSIP(k,j)
!> 3. Ensure symmetry: JASTROWALL(j,k) = JASTROWALL(k,j)
!>
!> Example
!> -------
!> Used in variational Monte Carlo for Sz-conserving wave function updates.
subroutine upjastrowallpsi_sz(nel, nelup, jastrowall, psip)
    implicit none
    integer j, k, nel, nelup
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nelup
        do k = j + 1, nelup
            jastrowall(k, j) = jastrowall(k, j) + psip(k, j)
        end do
        do k = nelup + 1, nel
            jastrowall(k, j) = jastrowall(k, j) - psip(k, j)
        end do
    end do
    do j = nelup + 1, nel
        do k = j + 1, nel
            jastrowall(k, j) = jastrowall(k, j) + psip(k, j)
        end do
    end do

    do j = 1, nel
        do k = j + 1, nel
            jastrowall(j, k) = jastrowall(k, j)
        end do
    end do

    return
end subroutine upjastrowallpsi_sz
