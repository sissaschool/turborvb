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

!> @brief Parallel matrix-matrix multiplication with MPI distribution
!>
!> This subroutine performs matrix-matrix multiplication C = α*A*B + β*C
!> in parallel using MPI. The K dimension is distributed across processors,
!> with each processor computing a partial result that is then reduced
!> to form the final result.
!>
!> Parameters
!> ----------
!> TRANSA : character*1, in
!>     Transpose flag for matrix A ('N'=no transpose, 'T'=transpose).
!> TRANSB : character*1, in
!>     Transpose flag for matrix B ('N'=no transpose, 'T'=transpose).
!> M : integer, in
!>     Number of rows in matrices A and C.
!> N : integer, in
!>     Number of columns in matrices B and C.
!> K : integer, in
!>     Number of columns in A and rows in B (distributed dimension).
!> ALPHA : double precision, in
!>     Scalar multiplier for A*B.
!> A : double precision array, in
!>     Input matrix A (size depends on TRANSA).
!> LDA : integer, in
!>     Leading dimension of matrix A.
!> B : double precision array, in
!>     Input matrix B (size depends on TRANSB).
!> LDB : integer, in
!>     Leading dimension of matrix B.
!> BETA : double precision, in
!>     Scalar multiplier for C.
!> C : double precision array, inout
!>     Input/output matrix C (M × N).
!> LDC : integer, in
!>     Leading dimension of matrix C.
!> nproc : integer, in
!>     Number of MPI processes.
!> rank : integer, in
!>     MPI rank of current process.
!> comm_mpi : integer, in
!>     MPI communicator.
!>
!> Notes
!> -----
!> - The K dimension is distributed across processors for parallel computation.
!> - Each processor computes a partial matrix multiplication with its local K slice.
!> - Results are reduced using MPI to form the final matrix C.
!> - Only rank 0 applies the BETA scaling; other ranks start with zero.
!> - Supports all four combinations of transpose operations (NN, NT, TN, TT).
!>
!> Example
!> -------
!> Used in distributed linear algebra operations throughout TurboRVB.
subroutine DGEMM_MY(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA&
        &, B, LDB, BETA, C, LDC, nproc, rank, comm_mpi)
    implicit none
#ifdef PARALLEL
    include 'mpif.h'
    integer ndim2, i, nu, nm, indr
#endif
    character*1 TRANSA, TRANSB
    integer M, N, K, LDA, LDB, LDC, nproc, rank, comm_mpi
    double precision ALPHA, BETA
    double precision A(LDA, *), B(LDB, *), C(LDC, *)
#ifdef PARALLEL
#ifdef __TEST
    integer dima, dimb, dimc
!> @brief Broadcast input matrices for testing consistency
    if (nproc .gt. 1) then
        if (transa .eq. 'N' .or. transa .eq. 'n') then
            dima = LDA*(K - 1) + M
        else
            dima = LDA*(M - 1) + K
        end if
        if (transb .eq. 'N' .or. transb .eq. 'n') then
            dimb = LDB*(N - 1) + K
        else
            dimb = LDB*(K - 1) + N
        end if
        dimc = LDC*(N - 1) + M
        call bcast_real(a, dima, 0, comm_mpi)
        call bcast_real(b, dimb, 0, comm_mpi)
        call bcast_real(c, dimc, 0, comm_mpi)
    end if
#endif
    ndim2 = LDC*(N - 1) + M
    if (nproc .gt. 1) then
!> @brief Calculate local K dimension for this processor
        nm = k/nproc
        if (nm*nproc .ne. k) nm = nm + 1
        indr = rank*nm + 1
        nu = nm
        if (indr + nm - 1 .gt. k) nu = k - indr + 1

        if (nu .gt. 0) then
!> @brief Handle all four transpose combinations with parallel computation
            if ((transa .eq. 'N' .or. transa .eq. 'n') .and. (transb .eq. 'N' .or. transb .eq. 'n')) then
!> @brief A*B (no transpose for either matrix)
                if (rank .ne. 0) then
                    call dgemm('N', 'N', m, n, nu, alpha, a(1, indr), lda&
         &, b(indr, 1), ldb, 0.d0, c, ldc)
                else
                    call dgemm('N', 'N', m, n, nu, alpha, a(1, indr), lda&
         &, b(indr, 1), ldb, beta, c, ldc)
                end if
            elseif ((transa .eq. 'N' .or. transa .eq. 'n') .and. (transb .eq. 'T' .or. transb .eq. 't')) then
!> @brief A*B^T (transpose B only)
                if (rank .ne. 0) then
                    call dgemm('N', 'T', m, n, nu, alpha, a(1, indr), lda&
         &, b(1, indr), ldb, 0.d0, c, ldc)
                else
                    call dgemm('N', 'T', m, n, nu, alpha, a(1, indr), lda&
         &, b(1, indr), ldb, beta, c, ldc)
                end if
            elseif ((transa .eq. 'T' .or. transa .eq. 't') .and. (transb .eq. 'N' .or. transb .eq. 'n')) then
!> @brief A^T*B (transpose A only)
                if (rank .ne. 0) then
                    call dgemm('T', 'N', m, n, nu, alpha, a(indr, 1), lda&
         &, b(indr, 1), ldb, 0.d0, c, ldc)
                else
                    call dgemm('T', 'N', m, n, nu, alpha, a(indr, 1), lda&
         &, b(indr, 1), ldb, beta, c, ldc)
                end if
            elseif ((transa .eq. 'T' .or. transa .eq. 't') .and. (transb .eq. 'T' .or. transb .eq. 't')) then
!> @brief A^T*B^T (transpose both matrices)
                if (rank .ne. 0) then
                    call dgemm('T', 'T', m, n, nu, alpha, a(indr, 1), lda&
         &, b(1, indr), ldb, 0.d0, c, ldc)
                else
                    call dgemm('T', 'T', m, n, nu, alpha, a(indr, 1), lda&
         &, b(1, indr), ldb, beta, c, ldc)
                end if
            end if
        elseif (rank .ne. 0) then
!> @brief Zero out result matrix for processors with no work
            do i = 1, n
                c(1:m, i) = 0.d0
            end do
        end if
        if (rank .ne. 0 .and. ldc .gt. m) then
!> @brief Zero out unused elements in leading dimension
            do i = 1, n - 1
                c(m + 1:ldc, i) = 0.d0
            end do
        end if
!> @brief Reduce partial results from all processors
        call reduce_base_real(ndim2, c, comm_mpi, -1)
    else
#endif
!> @brief Serial computation using standard BLAS DGEMM
        call DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA&
                &, B, LDB, BETA, C, LDC)
#ifdef PARALLEL
    end if
#endif
    return
end
