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

!> @brief Invert a complex skew-symmetric matrix using LU decomposition
!>
!> This subroutine computes the inverse of a complex skew-symmetric matrix A
!> using the LU decomposition A = L*U*P, where P is a permutation matrix.
!> The inverse is computed as A^(-1) = P^T * U^(-1) * L^(-1).
!>
!> Parameters
!> ----------
!> UPLO : character, in
!>     Specifies which part of the matrix A is stored:
!>     'U' or 'u': Upper triangular part is stored.
!>     'L' or 'l': Lower triangular part is stored.
!> N : integer, in
!>     Order of the matrix A.
!> A : complex*16 array, in
!>     The complex skew-symmetric matrix A, stored in packed format.
!>     On entry, contains the LU factors from ZSKTRF.
!> LDA : integer, in
!>     Leading dimension of array A.
!> AINV : complex*16 array, out
!>     The inverse of matrix A.
!> LDINV : integer, in
!>     Leading dimension of array AINV.
!> IPIV : integer array, in
!>     Pivot indices from ZSKTRF; dimension at least N.
!> WORK : complex*16 array, work
!>     Workspace array; dimension at least N^2 + 12*N - 2.
!> INFO : integer, out
!>     = 0: successful exit.
!>     < 0: if INFO = -i, the i-th argument had an illegal value.
!>     > 0: if INFO = i, U(i,i) is exactly zero; the matrix is singular.
!>
!> Notes
!> -----
!> - Assumes that ZSKTRF has been called previously to factorize A.
!> - Uses a three-step process: permutation, triangular solve, permutation.
!> - Supports both upper and lower triangular storage formats.
!> - The algorithm handles the special structure of complex skew-symmetric matrices.
!> - Complex arithmetic is used throughout the computation.
!>
!> Algorithm
!> ---------
!> 1. Initialize AINV as identity matrix
!> 2. Apply first permutation (reverse order)
!> 3. Extract skew-symmetric elements and prepare for triangular solve
!> 4. Solve triangular system using ZTRTRS
!> 5. Solve tridiagonal system using ZSKTRS or ZGTSVX
!> 6. Solve triangular system with transpose
!> 7. Apply final permutation
!>
!> Example
!> -------
!> Used in quantum Monte Carlo calculations for complex matrix operations.
subroutine zsktri(uplo, n, a, lda, ainv, ldinv, ipiv, work, info)
    implicit none

    ! argument parameters
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda, ldinv, ipiv(*)
    complex*16, intent(in out) :: a(lda, *), work(*)
    complex*16, intent(out) :: ainv(ldinv, *)
    integer, intent(out) :: info

    ! local variables
    integer i, j
    real*8 rcond
    logical yeslap

!> @brief Initialize INFO and choose algorithm
    info = 0
    yeslap = .false. ! if false the homemade algorithm is done.
!> @brief Initialize AINV as identity matrix
    do i = 1, n
        do j = 1, i - 1
            ainv(j, i) = dcmplx(0.d0, 0.d0)
        end do
        ainv(i, i) = dcmplx(1.d0, 0.d0)
        do j = i + 1, n
            ainv(j, i) = dcmplx(0.d0, 0.d0)
        end do
    end do
!> @brief Apply first permutation in reverse order
    do i = n, 1, -1
        work(1:n) = ainv(i, 1:n)
        ainv(i, 1:n) = ainv(ipiv(i), 1:n)
        ainv(ipiv(i), 1:n) = work(1:n)
    end do
!> @brief Extract skew-symmetric elements and prepare triangular solve
    if (UPLO .eq. 'u' .or. UPLO .eq. 'U') then
!> @brief Upper triangular case: extract superdiagonal elements
        do i = 1, N - 1
            work(i) = a(i, i + 1)
            work(n + i - 1) = -a(i, i + 1)
        end do
        do j = 2, n - 1
            a(1:j - 1, j) = a(1:j - 1, j + 1)
        end do
        a(1:n - 1, n) = dcmplx(0.d0, 0.d0)
    else
!> @brief Lower triangular case: extract subdiagonal elements
        do i = 1, N - 1
            work(i) = -a(i + 1, i)
            work(n + i - 1) = a(i + 1, i)
        end do
        do j = n - 1, 2, -1
            a(j + 1:n, j) = a(j + 1:n, j - 1)
        end do
        a(2:n, 1) = dcmplx(0.d0, 0.d0)
    end if
!> @brief Set diagonal elements to zero (skew-symmetric property)
    work(2*N - 1:3*N - 2) = dcmplx(0.d0, 0.d0) ! diagonal elements of skew matrix , obviously set to zero.
!> @brief First triangular solve: U * X = B
    call ZTRTRS(UPLO, 'N', 'U', N, N, A, LDA, ainv, LDINV, INFO) ! ainv = A^-1
!> @brief Solve tridiagonal system using LAPACK or homemade algorithm
    if (yeslap) then
!> @brief Use LAPACK ZGTSVX for complex tridiagonal solve
        call ZGTSVX('N', 'N', N, N, work(N), work(2*N - 1), work, work(3*N)&
                &, work(4*N), work(5*N), work(6*N), IPIV(N + 1), AINV, LDINV, work(7*N), N&
                &, RCOND, work(7*N + N*N - 1), work(8*N + N*N - 1)&
                &, work(9*N + N*N - 1), work(11*N + N*N - 1), INFO)
        do i = 1, N
            ainv(1:N, i) = work(7*N + (i - 1)*N:7*N + i*N - 1)
        end do
    else
!> @brief Use homemade ZSKTRS for complex skew-symmetric tridiagonal solve
        call ZSKTRS(UPLO, N, N, WORK, AINV, LDINV, WORK(N), INFO)
    end if
!> @brief Second triangular solve: L^T * X = B
    call ZTRTRS(UPLO, 'T', 'U', N, N, A, LDA, ainv, LDINV, INFO)
!> @brief Apply final permutation to complete A^(-1) = P^T * U^(-1) * L^(-1)
    do i = 1, n
        work(1:n) = ainv(i, 1:n)
        ainv(i, 1:n) = ainv(ipiv(i), 1:n)
        ainv(ipiv(i), 1:n) = work(1:n)
    end do
    return
end

