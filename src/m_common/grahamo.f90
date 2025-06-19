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

!> @brief Gram-Schmidt orthogonalization with standard Euclidean metric
!>
!> This module provides Gram-Schmidt orthogonalization subroutines for
!> quantum Monte Carlo calculations using the standard Euclidean metric
!> (identity overlap matrix). It includes both real and complex versions
!> with optimized memory usage and numerical stability considerations.
!>
!> The module implements:
!> - Real Gram-Schmidt orthogonalization with Euclidean metric
!> - Complex Gram-Schmidt orthogonalization with Euclidean metric
!> - Numerical stability checks using machine precision
!> - Efficient memory usage with workspace arrays
!>
!> @author TurboRVB group
!> @version 1.0
!> @date 2022

!> @brief Real Gram-Schmidt orthogonalization with Euclidean metric
!>
!> This subroutine performs Gram-Schmidt orthogonalization of real vectors
!> using the standard Euclidean metric (identity overlap matrix). The
!> algorithm includes numerical stability checks and proper handling of
!> linearly dependent vectors.
!>
!> @param[in,out] psi Matrix of vectors to be orthogonalized (ndime, *)
!> @param[in,out] sc Workspace array for projections (mh + ndime)
!> @param[out] rn Array of normalization factors for each vector
!> @param[in] ndime Number of components in each vector
!> @param[in] mh Number of vectors to orthogonalize
!> @param[out] info Information about orthogonalization success:
!>                   - info = 0: Successful orthogonalization
!>                   - info > 0: Number of linearly dependent vectors found
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Initializes workspace array sc with -1 values
!> 2. Normalizes all input vectors to unit length
!> 3. Checks if the first vector is linearly independent
!> 4. For each subsequent vector i = 2, ..., mh:
!>    - Projects vector i onto all previous orthogonalized vectors
!>    - Subtracts the projections to make it orthogonal
!>    - Normalizes the resulting orthogonal vector
!>    - Checks for linear dependence using machine precision
!> 5. Returns normalization factors and information about linear dependencies
!>
!> @note Uses standard Euclidean scalar product: a · b = sum_i a_i * b_i
!> @note Numerical stability is ensured using machine precision thresholds
!> @note Linearly dependent vectors are set to zero
!> @note The workspace array sc is used efficiently to avoid memory allocation
!> @note Uses BLAS operations (DNRM2, DGEMV, DSCAL, DCOPY) for efficiency
subroutine GRAHAMO(PSI, SC, RN, NDIME, MH, INFO)
    implicit none
    integer mh, ndime, mx, i, mhh, info
    real*8 PSI(NDIME, *), SC(MH + NDIME), RN(MH), dnrm2, cost, eps, epsmin
    real*8, external :: dlamch
    eps = 10000.d0*dlamch('E')
    epsmin = dlamch('S')
    MHH = MH + 1
    do I = 1, MH
        SC(I) = -1.d0
    end do

    !          First normalize the input orbitals
    do i = 1, MH
        COST = DNRM2(NDIME, PSI(1, i), 1)
        if (cost .gt. epsmin) then
            cost = 1.d0/cost
            call DSCAL(NDIME, COST, PSI(1, i), 1)
        else
            psi(:, i) = 0.d0
        end if
    end do
    RN(1) = DNRM2(NDIME, PSI, 1)
    if (RN(1) .gt. eps) then
        info = 0
    else
        info = 1
    end if
    do I = 2, MH
        call DGEMV('T', ndime, I - 1, 1.d0, PSI, ndime, PSI(1, I), 1, 0.d0, SC, 1)
        call DGEMV('N', NDIME, I, 1.d0, PSI, ndime, SC, 1, 0.d0, SC(MHH), 1)
        RN(I) = DNRM2(NDIME, SC(MHH), 1)
        if (RN(I) .gt. eps) then
            COST = -1.d0/RN(I)
            call DSCAL(NDIME, COST, SC(MHH), 1)
            call DCOPY(NDIME, SC(MHH), 1, PSI(1, I), 1)
        else
            psi(:, i) = 0.d0
            info = info + 1
        end if

    end do
    return
end subroutine GRAHAMO

!> @brief Complex Gram-Schmidt orthogonalization with Euclidean metric
!>
!> This subroutine performs Gram-Schmidt orthogonalization of complex vectors
!> using the standard Euclidean metric (identity overlap matrix). The
!> algorithm includes numerical stability checks and proper handling of
!> linearly dependent vectors using complex conjugate scalar products.
!>
!> @param[in,out] psi Matrix of complex vectors to be orthogonalized (ndime, *)
!> @param[in,out] sc Complex workspace array for projections (mh + ndime)
!> @param[out] rn Array of normalization factors for each vector
!> @param[in] ndime Number of components in each vector
!> @param[in] mh Number of vectors to orthogonalize
!> @param[out] info Information about orthogonalization success:
!>                   - info = 0: Successful orthogonalization
!>                   - info > 0: Number of linearly dependent vectors found
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Initializes workspace array sc with -1 values
!> 2. Normalizes all input vectors to unit length using complex norm
!> 3. Checks if the first vector is linearly independent
!> 4. For each subsequent vector i = 2, ..., mh:
!>    - Projects vector i onto all previous orthogonalized vectors using
!>      complex conjugate scalar product
!>    - Subtracts the projections to make it orthogonal
!>    - Normalizes the resulting orthogonal vector
!>    - Checks for linear dependence using machine precision
!> 5. Returns normalization factors and information about linear dependencies
!>
!> @note Uses complex conjugate scalar product: a · b = sum_i a_i * conjg(b_i)
!> @note Numerical stability is ensured using machine precision thresholds
!> @note Linearly dependent vectors are set to zero
!> @note The workspace array sc is used efficiently to avoid memory allocation
!> @note Uses BLAS operations (DZNRM2, ZGEMV, ZSCAL, ZCOPY) for efficiency
subroutine GRAHAMO_COMPLEX(PSI, SC, RN, NDIME, MH, INFO)
    use constants, only: zone, zmone, zzero
    implicit none
    integer mh, ndime, mx, i, mhh, info
    complex*16 PSI(NDIME, *), SC(MH + NDIME), cost
    real*8 RN(MH), eps, epsmin
    real*8, external :: dlamch, dznrm2
    eps = 10000.d0*dlamch('E') ! Relative machine precision
    epsmin = dlamch('S') ! Safe minimum that allows the inverse.

    MHH = MH + 1
    do I = 1, MH
        SC(I) = zmone
    end do
    !          First normalize the input orbitals
    do i = 1, MH
        COST = DZNRM2(NDIME, PSI(1, i), 1)
        if (dreal(cost) .gt. epsmin) then
            cost = zone/cost
            call ZSCAL(NDIME, COST, PSI(1, i), 1)
        else
            psi(:, i) = zzero
        end if
    end do

    RN(1) = DZNRM2(NDIME, PSI, 1)
    if (RN(1) .gt. eps) then
        info = 0
    else
        info = 1
    end if
    do I = 2, MH
        call ZGEMV('C', ndime, I - 1, zone, PSI, ndime, PSI(1, I), 1, zzero, SC, 1)
        call ZGEMV('N', NDIME, I, zone, PSI, ndime, SC, 1, zzero, SC(MHH), 1)
        RN(I) = DZNRM2(NDIME, SC(MHH), 1)
        if (RN(I) .gt. eps) then
            COST = zmone/RN(I)
            call ZSCAL(NDIME, COST, SC(MHH), 1)
            call ZCOPY(NDIME, SC(MHH), 1, PSI(1, I), 1)
        else
            psi(:, i) = zzero
            info = info + 1
        end if

    end do
    return
end
