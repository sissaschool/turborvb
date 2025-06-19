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

!> @brief Solve complex skew-symmetric tridiagonal system A*X = B with quad precision
!>
!> This subroutine solves the linear system A*X = B where A is a
!> non-singular complex skew-symmetric tridiagonal matrix of even dimension n.
!> The matrix A is stored in packed format with A(i+1,i) = a(i) and
!> A(i,i+1) = -a(i) for i = 1, ..., n-1. Uses quad-precision arithmetic
!> for improved numerical stability.
!>
!> Parameters
!> ----------
!> UPLO : character, in
!>     Specifies which part of the matrix A is stored:
!>     'U' or 'u': Upper triangular part is stored (A(i,i+1) = a(i)).
!>     'L' or 'l': Lower triangular part is stored (A(i+1,i) = a(i)).
!> N : integer, in
!>     Order of the matrix A (must be even).
!> NHRS : integer, in
!>     Number of right-hand sides.
!> A : complex*16 array, in
!>     The complex skew-symmetric tridiagonal matrix stored in packed format.
!>     For UPLO = 'U': A(i,i+1) = a(i), i = 1, ..., n-1.
!>     For UPLO = 'L': A(i+1,i) = a(i), i = 1, ..., n-1.
!> B : complex*16 array, inout
!>     Right-hand side matrix (N × NHRS). On exit, contains the solution X.
!> LDB : integer, in
!>     Leading dimension of array B.
!> X : complex*16 array, work
!>     Workspace array for storing intermediate results.
!> INFO : integer, out
!>     = 0: successful exit.
!>     < 0: if INFO = -i, the i-th argument had an illegal value.
!>     > 0: if INFO = i, A(i,i+1) is exactly zero; the matrix is singular.
!>
!> Notes
!> -----
!> - The matrix A must be of even dimension n.
!> - The solution is computed using forward and backward substitution.
!> - Even-indexed variables are solved forward, odd-indexed backward.
!> - The algorithm exploits the special structure of complex skew-symmetric matrices.
!> - For UPLO = 'L', the signs of the matrix elements are temporarily flipped.
!> - Uses quad-precision arithmetic for improved numerical accuracy.
!> - Processes each right-hand side separately for better cache performance.
!>
!> Algorithm
!> ---------
!> 1. For UPLO = 'L', flip signs of matrix elements
!> 2. For each right-hand side:
!>    a. Solve even-indexed variables forward: x(2), x(4), ..., x(n)
!>    b. Solve odd-indexed variables backward: x(n-1), x(n-3), ..., x(1)
!> 3. Copy solution back to B
!> 4. Restore original matrix signs if needed
!>
!> Example
!> -------
!> Used in high-precision complex skew-symmetric matrix operations.
subroutine zsktrs(uplo, n, nhrs, a, b, ldb, x, info)
    implicit none

    ! argument variables
    integer, intent(in) :: n, nhrs, ldb
    integer, intent(out) :: info
    complex*16, intent(inout) :: a(*), b(ldb, *)
    complex*16, intent(out) :: x(*)
    character, intent(in) :: uplo

    ! local variables
    integer i, j
#ifdef __PORT
    complex*16 aq, a2q, xo
#else
    complex*32 aq, a2q, xo
#endif

!> @brief Initialize INFO and flip signs for lower triangular storage
    info = 0
    if (uplo .eq. 'l' .or. uplo .eq. 'L') then
        a(1:n - 1) = -a(1:n - 1)
    end if
!> @brief Solve A*X = B for complex skew-symmetric tridiagonal matrix A
!> A is stored as A(i+1,i) = a(i), A(i,i+1) = -a(i), i = 1, ..., n-1
!> Matrix dimension n must be even
!> @brief Process each right-hand side separately
    do j = 1, nhrs
!> @brief Solve even-indexed variables forward with quad precision
        aq = a(1)
        xo = b(1, j)/aq
        x(2) = xo
        do i = 4, n, 2
            aq = a(i - 1)
            a2q = a(i - 2)
            xo = (a2q*xo + b(i - 1, j))/aq
            x(i) = xo
        end do
!> @brief Solve odd-indexed variables backward with quad precision
        aq = a(n - 1)
        xo = -b(n, j)/aq
        x(n - 1) = xo
        do i = n - 3, 1, -2
            aq = a(i)
            a2q = a(i + 1)
            xo = -(-a2q*xo + b(i + 1, j))/aq
            x(i) = xo
        end do
!> @brief Copy solution back to B
        do i = 1, n
            b(i, j) = x(i)
        end do
    end do
!> @brief Restore original matrix signs
    if (uplo .eq. 'l' .or. uplo .eq. 'L') then
        a(1:n - 1) = -a(1:n - 1)
    end if
    return
end
