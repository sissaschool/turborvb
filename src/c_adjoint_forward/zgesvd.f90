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

subroutine ZGESVD_MY(JOBZ, N, M, A, LDA, ZR, W, nproc, rank, comm_mpi)
    use dspev_module
    use allio, only: sub_comm_diag
    use constants, only: zzero, zone
    implicit none
    real*8, external :: dlamch, dznrm2
    real*8 machp, cnorm
    integer nrl, firsti, j, ind, ierr, dima, comm_mpi
    logical yesbcast
#ifdef PARALLEL
    include 'mpif.h'
#endif
    integer nproc, rank
    !     INPUT A common to all processors
    !     OUTPUT A contains the left eigenvectors in SVD stored culumnwise
    !            ZR contains the right eigenvectors in SVD stored rowwise
    !            (as in standard dgesvd in lapack).
    !
    !  -- LAPACK driver routine (version 2.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     September 30, 1994
    !   NB At variance of lapack the singular values W are sorted in ascending order W(i+1)>=W(i).
    ! Also the right eigenvectors of A are not conjugated as in lapack
    !   A= (Output A) lambda ZR^* in this notations.
    !
    !     .. Scalar Arguments ..
    character JOBZ
    integer LDA, N, M, I
    !     ..
    !     .. Array Arguments ..
    complex*16 A(LDA, N), ZR(LDA, *)
    real*8 W(*)

    complex*16, dimension(:, :), allocatable :: z, b, ap
    complex*16, dimension(:), allocatable :: workp
    integer, dimension(:), allocatable :: iwork, ifail
    if (N .le. 0) return

    if (sub_comm_diag%parent .ne. comm_mpi) then
        !         write(6,*) "Warning! dsyev communicator mismatch",sub_comm_diag%parent,comm_mpi
        !         INFO=1000
        !         return
        yesbcast = .true.
    else
        yesbcast = .false.
    end if
    if (sub_comm_diag%yesin) then

        !      NB the parallel input is passed to the wrapper through the
        !      WORK space ad follows:
        allocate (b(lda, N))
        b = zzero
        call zgemm_my('N', 'C', n, n, m, zone, a, lda, a, lda, zzero, b, lda&
                &, sub_comm_diag%nproc, sub_comm_diag%rank, sub_comm_diag%comm)
        !&,nproc,rank,comm_mpi)

        machp = DLAMCH('E')
        !      Now standard 1D distribution.
        !      compute the left eigenvectors
        nrl = N/sub_comm_diag%nproc
        if (nrl*sub_comm_diag%nproc .ne. N) nrl = nrl + 1
        allocate (ap(nrl, N), z(nrl, N))
        firsti = sub_comm_diag%rank + 1
        i = 0
        ap = zzero
        z = zzero
        ind = firsti
        do while (ind .le. N)
            do j = 1, N
                ap(i + 1, j) = b(ind, j)
            end do
            i = i + 1
            ind = sub_comm_diag%nproc*i + firsti
        end do

        call zdspev_drv_ss(jobz, ap, nrl, w, z, nrl, &
                &nrl, n, sub_comm_diag%nproc, sub_comm_diag%rank, sub_comm_diag%comm, machp, 0)
        !    &             nrl, n, nproc,rank, comm_mpi, machp , 0 )

        !      save conjugate original a in b

        b = dconjg(a)

        !      now gather a
        !      remove all the part not used in the processor
        if (sub_comm_diag%rank .ne. 0) then
            do i = 1, N - 1
                a(1:LDA, i) = zzero
            end do
            a(1:N, N) = zzero
        else
            do i = 1, N
                a(1:N, i) = zzero
            end do
        end if
        i = 0
        ind = firsti
        do while (ind .le. N)
            do j = 1, N
                a(ind, j) = z(i + 1, j)
            end do
            i = i + 1
            ind = sub_comm_diag%nproc*i + firsti
        end do

        deallocate (ap, z)
    else

        do i = 1, N - 1
            A(1:LDA, i) = zzero
        end do
        A(1:N, N) = zzero

    end if

#ifdef PARALLEL
    dima = 2*(LDA*(N - 1) + N)
    if (yesbcast) then
        if (sub_comm_diag%yesin) then
            call reduce_base_real(dima, A, sub_comm_diag%comm, -1)
        end if
        call bcast_real(A, dima, 0, comm_mpi)
    else
        call reduce_base_real(dima, A, comm_mpi, -1)
    end if
#endif

    if (sub_comm_diag%yesin) then
        !      Now compute the right eigenvectors and put in zr stored rowwise
        call zgemm_my('T', 'N', n, m, n, zone, a, lda, b, lda, zzero, zr, lda&
                &, sub_comm_diag%nproc, sub_comm_diag%rank, sub_comm_diag%comm)
        !      To make the same output of zgesvd
        !      call conjmat(n,m,zr,lda)
        !&,nproc,rank,comm_mpi)
        deallocate (b)
    end if

#ifdef PARALLEL
    dima = 2*(lda*(M - 1) + N)
    call bcast_real(zr, dima, 0, comm_mpi)
    call bcast_real(W, N, 0, comm_mpi)
#endif
    do i = 1, n
        if (w(i) .gt. 0.d0) then
            w(i) = dsqrt(w(i))
        else
            w(i) = dznrm2(m, zr(i, 1), lda)
        end if
        if (w(i) .gt. 0.d0) then
            do j = 1, m
                zr(i, j) = zr(i, j)/w(i)
            end do
        end if
    end do
    return
end
