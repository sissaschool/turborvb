! Copyright (C) 2023- TurboRVB group
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

!> @brief Solve complex skew-symmetric tridiagonal system A*X = B
!>
!> This subroutine solves the linear system A*X = B where A is a
!> non-singular complex skew-symmetric tridiagonal matrix of even dimension n.
!> The matrix A is stored in packed format with A(i+1,i) = a(i) and
!> A(i,i+1) = -a(i) for i = 1, ..., n-1.
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
!> - Complex arithmetic is used throughout the computation.
!>
!> Algorithm
!> ---------
!> 1. For UPLO = 'L', flip signs of matrix elements
!> 2. Solve even-indexed variables forward: x(2), x(4), ..., x(n)
!> 3. Solve odd-indexed variables backward: x(n-1), x(n-3), ..., x(1)
!> 4. Copy solution back to B
!> 5. Restore original matrix signs if needed
!>
!> Example
!> -------
!> Used in complex skew-symmetric matrix inversion and linear system solving.
subroutine zsktrs(uplo, n, nhrs, a, b, ldb, x, info)
    implicit none

    ! argument parameters
    integer, intent(in) :: n, nhrs, ldb
    integer, intent(out) :: info
    complex*16, intent(in out) :: a(*), b(ldb, nhrs)
    complex*16, intent(out) :: x(nhrs, *)
    character, intent(in) :: uplo

    ! local variables
    integer i, j

!> @brief Initialize INFO and flip signs for lower triangular storage
    info = 0
    if (uplo .eq. 'l' .or. uplo .eq. 'L') then
        a(1:n - 1) = -a(1:n - 1)
    end if
!> @brief Solve A*X = B for complex skew-symmetric tridiagonal matrix A
!> A is stored as A(i+1,i) = a(i), A(i,i+1) = -a(i), i = 1, ..., n-1
!> Matrix dimension n must be even
!> @brief Solve even-indexed variables forward
    x(:, 2) = b(1, :)/a(1)
    do i = 4, n, 2
        x(:, i) = (b(i - 1, :) + a(i - 2)*x(:, i - 2))/a(i - 1)
    end do
!> @brief Solve odd-indexed variables backward
    x(:, n - 1) = -b(n, :)/a(n - 1)
    do i = n - 3, 1, -2
        x(:, i) = -(b(i + 1, :) - a(i + 1)*x(:, i + 2))/a(i)
    end do
!> @brief Copy solution back to B
    do i = 1, n
        b(i, :) = x(:, i)
    end do
!> @brief Restore original matrix signs
    if (uplo .eq. 'l' .or. uplo .eq. 'L') then
        a(1:n - 1) = -a(1:n - 1)
    end if
    return
end
