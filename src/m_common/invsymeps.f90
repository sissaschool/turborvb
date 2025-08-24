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

!=======================================================================
!> @file invsymeps.f90
!> @brief Symmetric (Hermitian) matrix inversion and trace utilities
!> @details This module provides routines for inverting symmetric or Hermitian matrices
!>          using eigenvalue decomposition, as well as utility functions for computing
!>          traces and trace-like quantities for real and complex matrices.
!> @author TurboRVB group
!> @date 2022
!> @section Features
!> - Robust inversion of symmetric/Hermitian matrices with eigenvalue thresholding
!> - Handles both real and complex matrices (ipc=1 or 2)
!> - Computes traces and trace-like quantities for matrix products
!> - Parallelization support (MPI)
!=======================================================================

!-----------------------------------------------------------------------
!> @brief Invert a symmetric (Hermitian) matrix using eigenvalue decomposition
!> @details Computes the inverse of a symmetric (real or complex Hermitian) matrix
!>          by diagonalizing and thresholding small eigenvalues. Handles both real
!>          and complex cases, and supports parallel execution with MPI.
!> @param[in] ipc 1 for real, 2 for complex
!> @param[in] n Matrix dimension
!> @param[in,out] a Matrix to invert (overwritten with inverse)
!> @param[in] lda Leading dimension of a
!> @param[out] info Error flag (0=success)
!> @param[in] eps Eigenvalue threshold for inversion
!> @param[out] mine Number of excluded eigenvalues
!> @param[out] umat Eigenvector matrix
!> @param[out] eigmat Eigenvalues (and their inverses)
!> @param[in] nproc Number of MPI processes
!> @param[in] rank MPI rank
!> @param[in] comm_mpi MPI communicator
!> @note Uses DSYEV_MY or ZSYEV_MY for diagonalization
!> @note Eigenvalues below eps are excluded from inversion
!> @note Writes warnings for small or negative eigenvalues
!> @note Returns inverse in-place in a
subroutine invsymeps(ipc, n, a, lda, info, eps, mine, umat, eigmat, nproc, rank, comm_mpi)
    use constants, only: zone, zzero
    implicit none
    integer n, lda, info, i, j, mine, neig, rank, ierr, countexc, dimorb&
            &, dime, minen, nproc, comm_mpi, ipc
    real*8 a(ipc*lda, *), eps, umat(ipc*n, n), eigmat(n), abstol, dlamch, condnum
#ifdef PARALLEL
    include 'mpif.h'
#endif
    real*8, dimension(:, :), allocatable :: b
    real*8, dimension(:), allocatable :: eig
    !#ifndef __CASO
    !#endif
    allocate (b(ipc*n, n), eig(n))
    b = 0.d0
    eig = 0.d0
    mine = 1

    !       Now the inverse of this matrix using SDV diagonalization
    !       by eliminating the irrelevant directions
    ! maximum accuracy
    abstol = 0.d0

    do j = 1, n
        b(1:ipc*n, j) = a(1:ipc*n, j) ! save original matrix
    end do

    do j = 1, n
        umat(1:ipc*n, j) = a(1:ipc*n, j)
    end do
    if (ipc .eq. 1) then
        call DSYEV_MY('V', 'L', n, umat, n, eig, info, nproc, rank, comm_mpi)
    else
        call ZSYEV_MY('V', 'L', n, umat, n, eig, info, nproc, rank, comm_mpi)
    end if
    countexc = 1
    do i = 1, n
        if (eig(i) .le. -1d12) countexc = countexc + 1
    end do
    if (mine .ne. countexc) then
        if (rank .eq. 0) write (6, *) 'ERROR in dsyevx, some eigenv<= -1d12 '
        info = -1
    end if
    condnum = eig(1)/eig(n)
    i = 1
    do while (condnum .le. 0.d0)
        i = i + 1
        condnum = eig(i)/eig(n)
    end do
    if (rank .eq. 0) write (6, *) ' first non zero/Inverse condition number basis = ', i, condnum
    mine = 1
    do i = 1, n
        if (eig(i)/eig(n) .gt. eps) then ! the condition number criterium
            !        if(rank.eq.0) write(6,*) i,eig(i)
            eigmat(i) = dsqrt(1.d0/eig(i))
            eig(i) = dsqrt(eig(i))
        else
            mine = mine + 1
            if (rank .eq. 0 .and. eig(i) .le. eps .and. eig(i) .ne. -1d12)&
                    & write (6, *) ' warning small  eigenvalue !!! ', eig(i)
            !        do j=1,n
            !        write(6,*) j,umat(j,i)
            !        enddo
            eig(i) = 0.d0
            eigmat(i) = 0.d0
        end if
    end do
    if (info .gt. 0 .and. rank .eq. 0)                                    &
            & write (6, *) ' info > 0 in dsyevx !!! ', info

    !        define  b

    do j = 1, n
        do i = 1, ipc*n
            b(i, j) = umat(i, j)*eigmat(j)
        end do
    end do

    !        it should be authomatically symmetric/hermitian
    if (ipc .eq. 1) then
        call dgemm_my('N', 'T', n, n, n, 1.d0, b, n, b, n, 0.d0, a, lda, nproc, rank, comm_mpi)
    else
        call zgemm_my('N', 'C', n, n, n, zone, b, n, b, n, zzero, a, lda, nproc, rank, comm_mpi)
    end if

    if (mine .ne. 1 .and. rank .eq. 0) then
        write (6, *) ' Warning neglecting', mine - 1, 'over', n, 'directions in SVD'
    end if

    deallocate (b, eig)

    return
end subroutine invsymeps

!-----------------------------------------------------------------------
!> @brief Compute the sum of squared elements (trace-like) for a real matrix
!> @details Returns sum_{i,j} a(i,j)*a(j,i) for a real matrix a.
!> @param[in] n Matrix dimension
!> @param[in] a Input matrix (lda, *)
!> @param[in] lda Leading dimension of a
!> @return Sum of squared elements (trace-like quantity)
function tracemat(n, a, lda)
    implicit none
    integer i, j, n, lda
    real*8 tracemat, a(lda, *)
    tracemat = 0.d0
    do i = 1, n
        do j = 1, n
            tracemat = tracemat + a(i, j)*a(j, i)
        end do
    end do
    return
end

!-----------------------------------------------------------------------
!> @brief Compute the sum of squared elements (trace-like) for a real or complex matrix
!> @details Returns sum_{i,j} a(i,j)*a(j,i) for a real or complex matrix a.
!> @param[in] n Matrix dimension
!> @param[in] a Input matrix (ipc*lda, *)
!> @param[in] lda Leading dimension of a
!> @return Sum of squared elements (trace-like quantity)
function tracematc(n, a, lda)
    use constants, only: ipc
    implicit none
    integer i, j, n, lda
    real*8 tracematc, a(ipc*lda, *)
    tracematc = 0.d0
    if (ipc .eq. 1) then
        do i = 1, n
            do j = 1, n
                tracematc = tracematc + a(i, j)*a(j, i)
            end do
        end do
    else
        do i = 1, n
            do j = 1, n
                tracematc = tracematc + a(2*i - 1, j)*a(2*j - 1, i) + a(2*i, j)*a(2*j, i)
            end do
        end do
    end if
    return
end

!-----------------------------------------------------------------------
!> @brief Compute the sum of products of two matrices (trace-like) for real or complex matrices
!> @details Returns sum_{i,j} a(i,j)*b(j,i) for real or complex matrices a, b.
!> @param[in] n Number of rows
!> @param[in] m Number of columns
!> @param[in] a Input matrix (ipc*lda, *)
!> @param[in] lda Leading dimension of a
!> @param[in] b Input matrix (ipc*ldb, *)
!> @param[in] ldb Leading dimension of b
!> @return Sum of products (trace-like quantity)
function tracemat2(n, m, a, lda, b, ldb)
    use constants, only: ipc
    implicit none
    integer i, j, n, m, lda, ldb
    real*8 tracemat2, tracemat2c, a(ipc*lda, *), b(ipc*ldb, *)
    tracemat2 = 0.d0
    if (ipc .eq. 1) then
        do i = 1, n
            do j = 1, m
                tracemat2 = tracemat2 + a(i, j)*b(j, i)
            end do
        end do
    else
        tracemat2 = tracemat2c(n, m, a, lda, b, ldb)
    end if
    return
end

!-----------------------------------------------------------------------
!> @brief Compute the sum of products of two complex matrices (trace-like)
!> @details Returns sum_{i,j} a(i,j)*b(j,i) for complex matrices a, b.
!> @param[in] n Number of rows
!> @param[in] m Number of columns
!> @param[in] a Input complex matrix (lda, *)
!> @param[in] lda Leading dimension of a
!> @param[in] b Input complex matrix (ldb, *)
!> @param[in] ldb Leading dimension of b
!> @return Sum of products (trace-like quantity)
function tracemat2c(n, m, a, lda, b, ldb)
    implicit none
    integer i, j, n, m, lda, ldb
    real*8 tracemat2c
    complex*16 a(lda, *), b(ldb, *)
    tracemat2c = 0.d0
    do i = 1, n
        do j = 1, m
            tracemat2c = tracemat2c + a(i, j)*b(j, i)
        end do
    end do
    return
end

!-----------------------------------------------------------------------
!> @brief Compute the trace of a real matrix
!> @details Returns sum_{i} a(i,i) for a real matrix a.
!> @param[in] n Matrix dimension
!> @param[in] a Input matrix (lda, *)
!> @param[in] lda Leading dimension of a
!> @return Trace of the matrix
function tracetrue(n, a, lda)
    implicit none
    integer i, j, n, lda
    real*8 tracetrue, a(lda, *)
    tracetrue = 0.d0
    do i = 1, n
        tracetrue = tracetrue + a(i, i)
    end do
    return
end
