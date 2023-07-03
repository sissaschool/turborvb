! Copyright (C) 2022 TurboRVB based on a LAPACK routine by
! Copyright (C) 1994 Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
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

subroutine DSYEV_MY(JOBZ, UPLO, N, A, LDA, W, INFO, nproc, rank, comm_mpi)
    use allio, only: sub_comm_diag, nproc_diag
    implicit none
    logical yesbcast
    integer nproc, rank, comm_mpi
    double precision, dimension(:, :), allocatable :: z
    double precision, dimension(:), allocatable :: workp
    integer, dimension(:), allocatable :: iwork, ifail
    real*8 vl, vu, abstol
    integer il, iu, lwork
    integer, external :: ilaenv
    real*8, external :: dlamch
    !
    !     .. Scalar Arguments ..
    character JOBZ, UPLO
    integer INFO, LDA, N, M, I, J
    !     ..
    !     .. Array Arguments ..
    double precision A(LDA, *), W(*)

    info = 0

    if (UPLO .eq. 'L' .or. UPLO .eq. 'l') then
        do i = 1, N
            do j = i + 1, N
                A(i, j) = A(j, i)
            end do
        end do
    else
        do i = 1, N
            do j = i + 1, N
                A(j, i) = A(i, j)
            end do
        end do
    end if

    if (nproc .ne. 1 .and. nproc_diag .ne. 1) then

        if (nproc_diag .eq. 0) then
            !write(6,*) "enter caso with all processor"
            yesbcast = .false.
            call DSYEV_CASO(JOBZ, UPLO, N, A, LDA, W, INFO, nproc, rank, comm_mpi)
        else
            if (sub_comm_diag%parent .ne. comm_mpi) then
                !         write(6,*) "Warning! dsyev communicator mismatch",sub_comm_diag%parent,comm_mpi
                !         INFO=1000
                !         return
                yesbcast = .true.
            else
                yesbcast = .false.
            end if
            !       else
            !          write(6,*) "using sub_comm", sub_comm_diag%nproc, sub_comm_diag%rank
            if (sub_comm_diag%yesin) then
                !            write(6,*) "using A"
                call DSYEV_CASO(JOBZ, UPLO, N, A, LDA, W, INFO, sub_comm_diag%nproc, sub_comm_diag%rank, sub_comm_diag%comm)
                !            write(6,*) "using A done"
            else
                !            write(6,*) "using B"
                do i = 1, N - 1
                    A(1:LDA, i) = 0.d0
                end do
                A(1:N, N) = 0.d0
                !            write(6,*) "using B done"
            end if
            !       endif
        end if
#ifdef PARALLEL
        call bcast_real(W, N, 0, comm_mpi)
        if (yesbcast) then
            if (sub_comm_diag%yesin) then
                call reduce_base_real(LDA*(N - 1) + N, A, sub_comm_diag%comm, -1)
            end if
            call bcast_real(A, LDA*(N - 1) + N, 0, comm_mpi)
        else
            call reduce_base_real(LDA*(N - 1) + N, A, comm_mpi, -1)
        end if
#endif

        !     if(rank.eq.0) write(6,*) A(1:LDA,1:N)

    else

        !      DOUBLE PRECISION, dimension(:,:), allocatable:: z
        !      double PRECISION, dimension(:), allocatable:: workp
        !      integer, dimension(:), allocatable:: iwork,ifail
        if (UPLO .eq. 'L' .or. UPLO .eq. 'l') then
            lwork = ILAENV(1, 'DSYTRD', 'AL', N, N, N, N)
        else
            lwork = ILAENV(1, 'DSYTRD', 'AU', N, N, N, N)
        end if
        lwork = max((lwork + 2)*N, 3*N)
        !     allocate(z(LDA,N),iwork(5*N),ifail(N),workp(lwork))
        allocate (workp(lwork))
        !     vl=0.d0
        !     vu=0.d0
        !     il=1
        !     iu=1
        !     abstol=2*dlamch('s') ! maximum accuracy
        !       CALL DSYEVX( JOBZ, 'A', UPLO, N, A, LDA, VL, VU, il, iu,&
        !    &                   abstol, M, W, Z, LDA, WORKP, lwork, IWORK,&
        !    &                   IFAIL, INFO )
        call DSYEV(JOBZ, UPLO, N, A, LDA, W, WORKP, lwork, INFO)

        !       do i=1,N
        !       A(1:N,i)=Z(1:N,i)
        !       enddo

        !     deallocate(z,iwork,ifail,workp)
        deallocate (workp)
#ifdef PARALLEL
        call bcast_real(W, N, 0, comm_mpi)
        call bcast_real(A, LDA*(N - 1) + N, 0, comm_mpi)
#endif

    end if

    return
end subroutine DSYEV_MY

subroutine DSYEV_CASO(JOBZ, UPLO, N, A, LDA, W, INFO, nproc, rank, comm_mpi)
    use dspev_module

    implicit none
    integer nproc, rank, comm_mpi
    !
    !  -- LAPACK driver routine (version 2.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     September 30, 1994
    !
    !     .. Scalar Arguments ..
    character JOBZ, UPLO
    integer INFO, LDA, N, M, I
    !     ..
    !     .. Array Arguments ..
    double precision A(LDA, *), W(*)

    double precision, dimension(:, :), allocatable :: z, ap

    real*8, external :: dlamch
    real*8 machp
    integer nrl, firsti, j, ind, ierr, dima

#ifdef PARALLEL
    include 'mpif.h'
#endif

    if (N .le. 0) return
    !      NB the parallel input is passed to the wrapper through the
    !      WORK space ad follows:
    machp = DLAMCH('E')
    !      Now standard 1D distribution.
    nrl = N/nproc
    if (nrl*nproc .ne. N) nrl = nrl + 1
    allocate (ap(nrl, N), z(nrl, N))
    firsti = rank + 1
    i = 0
    ap = 0.d0
    z = 0.d0
    ind = firsti
    do while (ind .le. N)
        do j = 1, N
            ap(i + 1, j) = a(ind, j)
        end do
        i = i + 1
        ind = nproc*i + firsti
    end do

    call pdspev_drv_ss(jobz, ap, nrl, w, z, nrl, &
            &             nrl, n, nproc, rank, comm_mpi, machp, 0)

    info = 0

    !      now gather a
    !      remove all the part not used in the processor
    if (rank .ne. 0) then
        do i = 1, N - 1
            a(1:LDA, i) = 0.d0
        end do
        a(1:N, N) = 0.d0
    else
        do i = 1, N
            a(1:N, i) = 0.d0
        end do
    end if
    i = 0
    ind = firsti
    do while (ind .le. N)
        do j = 1, N
            a(ind, j) = z(i + 1, j)
        end do
        i = i + 1
        ind = nproc*i + firsti
    end do

    !#ifdef PARALLEL
    !       dima=LDA*(N-1)+N
    !       CALL reduce_base_real(dima,a,comm_mpi,-1)
    !#endif
    deallocate (ap, z)

    return
end subroutine DSYEV_CASO
