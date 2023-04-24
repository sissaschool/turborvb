! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2002 FPMD group
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

module descriptors
    !
    implicit none
    save

    integer ldim_cyclic, ldim_block_sca
    integer lind_block_sca
    integer gind_block_sca
    external ldim_cyclic, ldim_block_sca
    external lind_block_sca
    external gind_block_sca

    !  Descriptor for Cannon's algorithm
    !
    !  Parameters to define and manage the Descriptor
    !  of square matricxes block distributed on a square grid of processors
    !  to be used with Cannon's algorithm for matrix multiplication
    !
    integer, parameter :: descla_siz_ = 16
    integer, parameter :: ilar_ = 1
    integer, parameter :: nlar_ = 2
    integer, parameter :: ilac_ = 3
    integer, parameter :: nlac_ = 4
    integer, parameter :: nlax_ = 5
    integer, parameter :: lambda_node_ = 6
    integer, parameter :: la_n_ = 7
    integer, parameter :: la_nx_ = 8
    integer, parameter :: la_npr_ = 9
    integer, parameter :: la_npc_ = 10
    integer, parameter :: la_myr_ = 11
    integer, parameter :: la_myc_ = 12
    integer, parameter :: la_comm_ = 13
    integer, parameter :: la_me_ = 14
    integer, parameter :: la_nrl_ = 15
    integer, parameter :: la_nrlx_ = 16
    !
    ! desc( ilar_ )  global index of the first row in the local block of lambda
    ! desc( nlar_ )  number of row in the local block of lambda ( the "2" accounts for spin)
    ! desc( ilac_ )  global index of the first column in the local block of lambda
    ! desc( nlac_ )  number of column in the local block of lambda
    ! desc( nlax_ )  leading dimension of the distribute lambda matrix
    ! desc( lambda_node_ )  if > 0 the proc holds a block of the lambda matrix
    ! desc( la_n_ )     global dimension of the matrix
    ! desc( la_nx_ )    global leading dimension
    ! desc( la_npr_ )   number of row processors
    ! desc( la_npc_ )   number of column processors
    ! desc( la_myr_ )   processor row index
    ! desc( la_myc_ )   processor column index
    ! desc( la_comm_ )  communicator
    ! desc( la_me_ ) processor index ( from 0 to desc( la_npr_ ) * desc( la_npc_ ) - 1 )
    ! desc( la_nrl_ ) number of local row, when the matrix is cyclically distributed across proc
    ! desc( la_nrlx_ ) leading dimension, when the matrix is distributed by row

    integer :: descla(descla_siz_)

contains

    !------------------------------------------------------------------------
    !
    subroutine descla_local_dims(i2g, nl, n, nx, np, me)
        implicit none
        integer, intent(OUT) :: i2g !  global index of the first local element
        integer, intent(OUT) :: nl !  local number of elements
        integer, intent(IN) :: n !  number of actual element in the global array
        integer, intent(IN) :: nx !  dimension of the global array (nx>=n) to be distributed
        integer, intent(IN) :: np !  number of processors
        integer, intent(IN) :: me !  taskid for which i2g and nl are computed
        !
        !  note that we can distribute a global array larger than the
        !  number of actual elements. This could be required for performance
        !  reasons, and to have an equal partition of matrix having different size
        !  like matrixes of spin-up and spin-down
        !
        nl = ldim_block_sca(nx, np, me)
        i2g = gind_block_sca(1, nx, np, me)
        !
        ! This is to try to keep a matrix N * N into the same
        ! distribution of a matrix NX * NX, useful to have
        ! the matrix of spin-up distributed in the same way
        ! of the matrix of spin-down
        !
        if (i2g + nl - 1 > n) nl = n - i2g + 1 ! shifting the index to maintain right dimension
        if (nl < 0) nl = 0
        return
        !
    end subroutine descla_local_dims
    !
    !
    subroutine descla_init(desc, n, nx, np, me, comm, includeme)
        !
        implicit none
        integer, intent(OUT) :: desc(:)
        integer, intent(IN) :: n !  the size of this matrix
        integer, intent(IN) :: nx !  the max among different matrices sharing
        !  this descriptor or the same data distribution
        integer, intent(IN) :: np(2), me(2), comm
        integer, intent(IN) :: includeme
        integer :: ir, nr, ic, nc, lnode, nlax, nrl, nrlx
        integer :: ip, npp

        if (np(1) /= np(2)) &
            call errore(' descla_init ', ' only square grid of proc are allowed ', 2)
        if (n < 0) &
            call errore(' descla_init ', ' dummy argument n less than 1 ', 3)
        if (nx < n) &
            call errore(' descla_init ', ' dummy argument nx less than n ', 4)
        if (np(1) < 1) &
            call errore(' descla_init ', ' dummy argument np less than 1 ', 5)

        ! find the block maximum dimensions

        nlax = ldim_block_sca(nx, np(1), 0)
        !
        ! find local block dimensions, if appropriate
        ! only for processes involved in the distribution
        !
        if (includeme == 1) then
            !
            call descla_local_dims(ir, nr, n, nx, np(1), me(1))
            call descla_local_dims(ic, nc, n, nx, np(2), me(2))
            !
            lnode = 1
            !
        else
            !
            nr = 0
            nc = 0
            !
            ir = 0
            ic = 0
            !
            lnode = -1
            !
        end if

        desc(ilar_) = ir
        desc(nlar_) = nr
        desc(ilac_) = ic
        desc(nlac_) = nc
        desc(nlax_) = nlax
        desc(lambda_node_) = lnode
        desc(la_n_) = n
        desc(la_nx_) = nx
        desc(la_npr_) = np(1)
        desc(la_npc_) = np(2)
        desc(la_myr_) = me(1)
        desc(la_myc_) = me(2)
        desc(la_comm_) = comm
        desc(la_me_) = desc(la_myc_) + desc(la_myr_)*desc(la_npr_)

        npp = np(1)*np(2)

        !  Compute local dimension of the cyclically distributed matrix
        !
        if (includeme == 1) then
            nrl = ldim_cyclic(n, npp, desc(la_me_))
        else
            nrl = 0
        end if
        nrlx = n/npp + 1

        desc(la_nrl_) = nrl
        desc(la_nrlx_) = nrlx

        if (nr < 0 .or. nc < 0) &
            call errore(' descla_init ', ' wrong valune for computed nr and nc ', 1)
        if (nlax < 1) &
            call errore(' descla_init ', ' wrong value for computed nlax ', 2)
        if (nlax < nr) &
            call errore(' descla_init ', ' nlax < nr ', (nr - nlax))
        if (nlax < nc) &
            call errore(' descla_init ', ' nlax < nc ', (nc - nlax))
        if (nrlx < nrl) &
            call errore(' descla_init ', ' nrlx < nrl ', (nrl - nrlx))
        if (nrl < 0) &
            call errore(' descla_init ', ' nrl < 0 ', abs(nrl))

        ! WRITE(*,*) 'me1,me2,nr,nc,ir,ic= ', me(1), me(2), nr, nc, ir, ic

        return
    end subroutine descla_init

#ifdef __SCALAPACK

    subroutine dsqmsym(n, a, lda, desc)
        !
        ! Double precision SQuare Matrix SYMmetrization
        !
        implicit none
        !
        integer, intent(IN) :: n
        integer, intent(IN) :: lda
        real*8 :: a(lda, *)
        integer, intent(IN) :: desc(descla_siz_)
#if defined PARALLEL
        include 'mpif.h'
        integer :: istatus(MPI_STATUS_SIZE)
#endif
        integer :: i, j
        integer :: comm
        integer :: nr, nc, dest, sreq, ierr, sour
        real*8 :: atmp

#if defined PARALLEL

        if (desc(lambda_node_) <= 0) then
            return
        end if

        if (n /= desc(la_n_)) &
            call errore(" dsqmsym ", " wrong global dim n ", n)
        if (lda /= desc(nlax_)) &
            call errore(" dsqmsym ", " wrong leading dim lda ", lda)

        comm = desc(la_comm_)

        nr = desc(nlar_)
        nc = desc(nlac_)
        if (desc(la_myc_) == desc(la_myr_)) then
            !
            !  diagonal block, procs work locally
            !
            do j = 1, nc
                do i = j + 1, nr
                    a(i, j) = a(j, i)
                end do
            end do
            !
        else if (desc(la_myc_) > desc(la_myr_)) then
            !
            !  super diagonal block, procs send the block to sub diag.
            !
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), dest)
            call mpi_isend(a, lda*lda, MPI_DOUBLE_PRECISION, dest, 1, comm, sreq, ierr)
            !
            if (ierr /= 0) &
                call errore(" dsqmsym ", " in isend ", abs(ierr))
            !
        else if (desc(la_myc_) < desc(la_myr_)) then
            !
            !  sub diagonal block, procs receive the block from super diag,
            !  then transpose locally
            !
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), sour)
            call mpi_recv(a, lda*lda, MPI_DOUBLE_PRECISION, sour, 1, comm, istatus, ierr)
            !
            if (ierr /= 0) &
                call errore(" dsqmsym ", " in recv ", abs(ierr))
            !
            do j = 1, lda
                do i = j + 1, lda
                    atmp = a(i, j)
                    a(i, j) = a(j, i)
                    a(j, i) = atmp
                end do
            end do
            !
        end if

        if (desc(la_myc_) > desc(la_myr_)) then
            !
            call MPI_Wait(sreq, istatus, ierr)
            !
            if (ierr /= 0) &
                call errore(" dsqmsym ", " in wait ", abs(ierr))
            !
        end if

#else

        do j = 1, n
            !
            do i = j + 1, n
                !
                a(i, j) = a(j, i)
                !
            end do
            !
        end do

#endif

        return
    end subroutine dsqmsym

    subroutine zsqmher(n, a, lda, desc)
        !
        ! double complex (Z) SQuare Matrix HERmitianize
        !
        implicit none
        !
        integer, intent(IN) :: n
        integer, intent(IN) :: lda
        complex(8) :: a(lda, lda)
        integer, intent(IN) :: desc(descla_siz_)
#if defined PARALLEL
        include 'mpif.h'
        integer :: istatus(MPI_STATUS_SIZE)
#endif
        integer :: i, j
        integer :: comm, myid
        integer :: nr, nc, dest, sreq, ierr, sour
        complex(8) :: atmp

#if defined PARALLEL

        if (desc(lambda_node_) <= 0) then
            return
        end if

        if (n /= desc(la_n_)) &
            call errore(" zsqmher ", " wrong global dim n ", n)
        if (lda /= desc(nlax_)) &
            call errore(" zsqmher ", " wrong leading dim lda ", lda)

        comm = desc(la_comm_)

        nr = desc(nlar_)
        nc = desc(nlac_)
        if (desc(la_myc_) == desc(la_myr_)) then
            !
            !  diagonal block, procs work locally
            !
            do j = 1, nc
                a(j, j) = DCMPLX(dble(a(j, j)))
                do i = j + 1, nr
                    a(i, j) = conjg(a(j, i))
                end do
            end do
            !
        else if (desc(la_myc_) > desc(la_myr_)) then
            !
            !  super diagonal block, procs send the block to sub diag.
            !
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), dest)
            call mpi_isend(a, lda*lda, MPI_DOUBLE_COMPLEX, dest, 1, comm, sreq, ierr)
            !
            if (ierr /= 0) &
                call errore(" zsqmher ", " in mpi_isend ", abs(ierr))
            !
        else if (desc(la_myc_) > desc(la_myr_)) then
            !
            !  sub diagonal block, procs receive the block from super diag,
            !  then transpose locally
            !
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), sour)
            call mpi_recv(a, lda*lda, MPI_DOUBLE_COMPLEX, sour, 1, comm, istatus, ierr)
            !
            if (ierr /= 0) &
                call errore(" zsqmher ", " in mpi_recv ", abs(ierr))
            !
            do j = 1, lda
                do i = j + 1, lda
                    atmp = a(i, j)
                    a(i, j) = a(j, i)
                    a(j, i) = atmp
                end do
            end do
            do j = 1, nc
                do i = 1, nr
                    a(i, j) = conjg(a(i, j))
                end do
            end do
            !
        end if

        if (desc(la_myc_) > desc(la_myr_)) then
            !
            call MPI_Wait(sreq, istatus, ierr)
            !
            if (ierr /= 0) &
                call errore(" zsqmher ", " in MPI_Wait ", abs(ierr))
            !
        end if

#else

        do j = 1, n
            !
            a(j, j) = DCMPLX(dble(a(j, j)))
            !
            do i = j + 1, n
                !
                a(i, j) = conjg(a(j, i))
                !
            end do
            !
        end do

#endif

        return
    end subroutine zsqmher

#endif
    !

    subroutine cyc2blk_redist(n, a, lda, b, ldb, desc)
        !
        !  Parallel square matrix redistribution.
        !  A (input) is cyclically distributed by rows across processors
        !  B (output) is distributed by block across 2D processors grid
        !
        implicit none
        !
        integer, intent(IN) :: n
        integer, intent(IN) :: lda, ldb
        real*8 :: a(lda, *), b(ldb, *)
        integer :: desc(*)
        !
#if defined (PARALLEL)
        !
        include 'mpif.h'
        !
#endif
        !
        integer :: ierr, itag
        integer :: np, ip, me, nproc, comm_a
        integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
        integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
        !
        real*8, allocatable :: rcvbuf(:, :, :)
        real*8, allocatable :: sndbuf(:, :)
        integer, allocatable :: ip_desc(:, :)
        !
        character(len=256) :: msg
        !
#if defined (PARALLEL)

        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        if (np /= desc(la_npc_)) &
            call errore(' cyc2blk_redist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' cyc2blk_redist ', ' n less or equal zero ', 1)
        !  IF( desc( la_n_ ) < nproc ) &
        !  & CALL errore( ' cyc2blk_redist ', ' Dimension matrix less than the number of proc ', 1 )

        allocate (ip_desc(descla_siz_, nproc))
        ip_desc = 0

        call mpi_barrier(comm_a, ierr)

        call mpi_allgather(desc, descla_siz_, mpi_integer, ip_desc, descla_siz_, mpi_integer, comm_a, ierr)
        if (ierr /= 0) &
            call errore(" cyc2blk_redist ", " in mpi_allgather ", abs(ierr))
        !
        nbuf = (nb/nproc + 2)*nb
        !
        allocate (sndbuf(nb/nproc + 2, nb))
        allocate (rcvbuf(nb/nproc + 2, nb, nproc))

        sndbuf = 0.d0
        rcvbuf = 0.d0

        do ip = 0, nproc - 1
            !
            if (ip_desc(nlax_, ip + 1) /= nb) &
                call errore(' cyc2blk_redist ', ' inconsistent block dim nb ', 1)
            !
            if (ip_desc(lambda_node_, ip + 1) > 0) then

                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
                !
                do j = 1, ip_nc
                    jj = j + ip_ic - 1
                    il = 1
                    do i = 1, ip_nr
                        ii = i + ip_ir - 1
                        if (mod(ii - 1, nproc) == me) then
                            call check_sndbuf_index()
                            sndbuf(il, j) = a((ii - 1)/nproc + 1, jj)
                            il = il + 1
                        end if
                    end do
                end do

            end if

            call mpi_barrier(comm_a, ierr)

            call mpi_gather(sndbuf, nbuf, mpi_double_precision, &
                            rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr)
            if (ierr /= 0) &
                call errore(" cyc2blk_redist ", " in mpi_gather ", abs(ierr))

        end do

        !
        nr = desc(nlar_)
        nc = desc(nlac_)
        ir = desc(ilar_)
        ic = desc(ilac_)
        !

        do ip = 0, nproc - 1
            do j = 1, nc
                il = 1
                do i = 1, nr
                    ii = i + ir - 1
                    if (mod(ii - 1, nproc) == ip) then
                        call check_rcvbuf_index()
                        b(i, j) = rcvbuf(il, j, ip + 1)
                        il = il + 1
                    end if
                end do
            end do
        end do
        !
        !
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else

        b(1:n, 1:n) = a(1:n, 1:n)

#endif

        return

    contains

        subroutine check_sndbuf_index()
            character(LEN=38), save :: msg = ' check_sndbuf_index in cyc2blk_redist '
            if (j > size(sndbuf, 2)) call errore(msg, ' j > SIZE(sndbuf,2) ', ip + 1)
            if (il > size(sndbuf, 1)) call errore(msg, ' il > SIZE(sndbuf,1) ', ip + 1)
            if ((ii - 1)/nproc + 1 < 1) call errore(msg, ' ( ii - 1 )/nproc + 1 < 1 ', ip + 1)
            if ((ii - 1)/nproc + 1 > size(a, 1)) call errore(msg, ' ( ii - 1 )/nproc + 1 > SIZE(a,1) ', ip + 1)
            if (jj < 1) call errore(msg, ' jj < 1 ', ip + 1)
            if (jj > n) call errore(msg, ' jj > n ', ip + 1)
            return
        end subroutine check_sndbuf_index

        subroutine check_rcvbuf_index()
            character(LEN=38), save :: msg = ' check_rcvbuf_index in cyc2blk_redist '
            if (i > ldb) call errore(msg, ' i > ldb ', ip + 1)
            if (j > ldb) call errore(msg, ' j > ldb ', ip + 1)
            if (j > nb) call errore(msg, ' j > nb  ', ip + 1)
            if (il > size(rcvbuf, 1)) call errore(msg, ' il too large ', ip + 1)
            return
        end subroutine check_rcvbuf_index

    end subroutine cyc2blk_redist

    subroutine cyc2blk_zredist(n, a, lda, b, ldb, desc)
        !
        !  Parallel square matrix redistribution.
        !  A (input) is cyclically distributed by rows across processors
        !  B (output) is distributed by block across 2D processors grid
        !
        implicit none
        !
        integer, intent(IN) :: n
        integer, intent(IN) :: lda, ldb
        complex(8) :: a(lda, *), b(ldb, *)
        integer :: desc(*)
        !
#if defined (PARALLEL)
        !
        include 'mpif.h'
        !
#endif
        !
        integer :: ierr, itag
        integer :: np, ip, me, nproc, comm_a
        integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
        integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
        !
        complex(8), allocatable :: rcvbuf(:, :, :)
        complex(8), allocatable :: sndbuf(:, :)
        integer, allocatable :: ip_desc(:, :)
        !
        character(len=256) :: msg
        !
#if defined (PARALLEL)

        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        if (np /= desc(la_npc_)) &
            call errore(' cyc2blk_zredist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' cyc2blk_zredist ', ' n less or equal zero ', 1)
        !  IF( desc( la_n_ ) < nproc ) &
        !  & CALL errore( ' cyc2blk_redist ', ' Dimension matrix less than the number of proc ', 1 )

        allocate (ip_desc(descla_siz_, nproc))
        ip_desc = 0

        call mpi_barrier(comm_a, ierr)

        call mpi_allgather(desc, descla_siz_, mpi_integer, ip_desc, descla_siz_, mpi_integer, comm_a, ierr)
        if (ierr /= 0) &
            call errore(" cyc2blk_zredist ", " in mpi_allgather ", abs(ierr))
        !
        nbuf = (nb/nproc + 2)*nb
        !
        allocate (sndbuf(nb/nproc + 2, nb))
        allocate (rcvbuf(nb/nproc + 2, nb, nproc))

        sndbuf = (0.d0, 0.d0)
        rcvbuf = (0.d0, 0.d0)

        do ip = 0, nproc - 1
            !
            if (ip_desc(nlax_, ip + 1) /= nb) &
                call errore(' cyc2blk_zredist ', ' inconsistent block dim nb ', 1)
            !
            if (ip_desc(lambda_node_, ip + 1) > 0) then

                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
                !
                do j = 1, ip_nc
                    jj = j + ip_ic - 1
                    il = 1
                    do i = 1, ip_nr
                        ii = i + ip_ir - 1
                        if (mod(ii - 1, nproc) == me) then
                            call check_sndbuf_index()
                            sndbuf(il, j) = a((ii - 1)/nproc + 1, jj)
                            il = il + 1
                        end if
                    end do
                end do

            end if

            call mpi_barrier(comm_a, ierr)

            call mpi_gather(sndbuf, nbuf, mpi_double_complex, &
                            rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr)
            if (ierr /= 0) &
                call errore(" cyc2blk_zredist ", " in mpi_gather ", abs(ierr))

        end do

        !
        nr = desc(nlar_)
        nc = desc(nlac_)
        ir = desc(ilar_)
        ic = desc(ilac_)
        !

        do ip = 0, nproc - 1
            do j = 1, nc
                il = 1
                do i = 1, nr
                    ii = i + ir - 1
                    if (mod(ii - 1, nproc) == ip) then
                        call check_rcvbuf_index()
                        b(i, j) = rcvbuf(il, j, ip + 1)
                        il = il + 1
                    end if
                end do
            end do
        end do
        !
        !
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else

        b(1:n, 1:n) = a(1:n, 1:n)

#endif

        return

    contains

        subroutine check_sndbuf_index()
            character(LEN=38), save :: msg = ' check_sndbuf_index in cyc2blk_redist '
            if (j > size(sndbuf, 2)) call errore(msg, ' j > SIZE(sndbuf,2) ', ip + 1)
            if (il > size(sndbuf, 1)) call errore(msg, ' il > SIZE(sndbuf,1) ', ip + 1)
            if ((ii - 1)/nproc + 1 < 1) call errore(msg, ' ( ii - 1 )/nproc + 1 < 1 ', ip + 1)
            if ((ii - 1)/nproc + 1 > size(a, 1)) call errore(msg, ' ( ii - 1 )/nproc + 1 > SIZE(a,1) ', ip + 1)
            if (jj < 1) call errore(msg, ' jj < 1 ', ip + 1)
            if (jj > n) call errore(msg, ' jj > n ', ip + 1)
            return
        end subroutine check_sndbuf_index

        subroutine check_rcvbuf_index()
            character(LEN=38), save :: msg = ' check_rcvbuf_index in cyc2blk_redist '
            if (i > ldb) call errore(msg, ' i > ldb ', ip + 1)
            if (j > ldb) call errore(msg, ' j > ldb ', ip + 1)
            if (j > nb) call errore(msg, ' j > nb  ', ip + 1)
            if (il > size(rcvbuf, 1)) call errore(msg, ' il too large ', ip + 1)
            return
        end subroutine check_rcvbuf_index

    end subroutine cyc2blk_zredist

    subroutine blk2cyc_redist(n, a, lda, b, ldb, desc)
        !
        !  Parallel square matrix redistribution.
        !  A (output) is cyclically distributed by rows across processors
        !  B (input) is distributed by block across 2D processors grid
        !
        implicit none
        !
        integer, intent(IN) :: n
        integer, intent(IN) :: lda, ldb
        real*8 :: a(lda, *), b(ldb, *)
        integer :: desc(*)
        !
#if defined (PARALLEL)
        !
        include 'mpif.h'
        !
#endif
        !
        integer :: ierr, itag
        integer :: np, ip, me, comm_a, nproc
        integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
        integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
        !
        real*8, allocatable :: rcvbuf(:, :, :)
        real*8, allocatable :: sndbuf(:, :)
        integer, allocatable :: ip_desc(:, :)
        !
        character(len=256) :: msg
        !
#if defined (PARALLEL)

        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        if (np /= desc(la_npc_)) &
            call errore(' blk2cyc_redist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' blk2cyc_redist ', ' n less or equal zero ', 1)
        !  IF( desc( la_n_ ) < nproc ) &
        ! & CALL errore( ' blk2cyc_redist ', ' Dimension matrix less than the number of proc ', 1 )

        allocate (ip_desc(descla_siz_, nproc))

        ip_desc = 0

        call mpi_barrier(comm_a, ierr)

        call mpi_allgather(desc, descla_siz_, mpi_integer, ip_desc, descla_siz_, mpi_integer, comm_a, ierr)
        if (ierr /= 0) &
            call errore(" blk2cyc_redist ", " in mpi_allgather ", abs(ierr))
        !
        nbuf = (nb/nproc + 2)*nb
        !
        allocate (sndbuf(nb/nproc + 2, nb))
        allocate (rcvbuf(nb/nproc + 2, nb, nproc))

        sndbuf = 0.d0
        rcvbuf = 0.d0
        !
        nr = desc(nlar_)
        nc = desc(nlac_)
        ir = desc(ilar_)
        ic = desc(ilac_)
        !
        do ip = 0, nproc - 1
            do j = 1, nc
                il = 1
                do i = 1, nr
                    ii = i + ir - 1
                    if (mod(ii - 1, nproc) == ip) then
                        sndbuf(il, j) = b(i, j)
                        il = il + 1
                    end if
                end do
            end do
            call mpi_barrier(comm_a, ierr)
            call mpi_gather(sndbuf, nbuf, mpi_double_precision, &
                            rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr)
            if (ierr /= 0) &
                call errore(" blk2cyc_redist ", " in mpi_gather ", abs(ierr))
        end do
        !

        do ip = 0, nproc - 1
            !
            if (ip_desc(lambda_node_, ip + 1) > 0) then

                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
                !
                do j = 1, ip_nc
                    jj = j + ip_ic - 1
                    il = 1
                    do i = 1, ip_nr
                        ii = i + ip_ir - 1
                        if (mod(ii - 1, nproc) == me) then
                            a((ii - 1)/nproc + 1, jj) = rcvbuf(il, j, ip + 1)
                            il = il + 1
                        end if
                    end do
                end do

            end if

        end do
        !
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else

        a(1:n, 1:n) = b(1:n, 1:n)

#endif

        return

    end subroutine blk2cyc_redist

    subroutine blk2cyc_zredist(n, a, lda, b, ldb, desc)
        !
        !  Parallel square matrix redistribution.
        !  A (output) is cyclically distributed by rows across processors
        !  B (input) is distributed by block across 2D processors grid
        !
        implicit none
        !
        integer, intent(IN) :: n
        integer, intent(IN) :: lda, ldb
        complex(8) :: a(lda, *), b(ldb, *)
        integer :: desc(*)
        !
#if defined (PARALLEL)
        !
        include 'mpif.h'
        !
#endif
        !
        integer :: ierr, itag
        integer :: np, ip, me, comm_a, nproc
        integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
        integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
        !
        complex(8), allocatable :: rcvbuf(:, :, :)
        complex(8), allocatable :: sndbuf(:, :)
        integer, allocatable :: ip_desc(:, :)
        !
        character(len=256) :: msg
        !
#if defined (PARALLEL)

        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        if (np /= desc(la_npc_)) &
            call errore(' blk2cyc_zredist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' blk2cyc_zredist ', ' n less or equal zero ', 1)
        !  IF( desc( la_n_ ) < nproc ) &
        ! & CALL errore( ' blk2cyc_zredist ', ' Dimension matrix less than the number of proc ', 1 )

        allocate (ip_desc(descla_siz_, nproc))

        ip_desc = 0

        call mpi_barrier(comm_a, ierr)

        call mpi_allgather(desc, descla_siz_, mpi_integer, ip_desc, descla_siz_, mpi_integer, comm_a, ierr)
        if (ierr /= 0) &
            call errore(" blk2cyc_zredist ", " in mpi_allgather ", abs(ierr))
        !
        nbuf = (nb/nproc + 2)*nb
        !
        allocate (sndbuf(nb/nproc + 2, nb))
        allocate (rcvbuf(nb/nproc + 2, nb, nproc))

        sndbuf = (0.d0, 0.d0)
        rcvbuf = (0.d0, 0.d0)
        !
        nr = desc(nlar_)
        nc = desc(nlac_)
        ir = desc(ilar_)
        ic = desc(ilac_)
        !
        do ip = 0, nproc - 1
            do j = 1, nc
                il = 1
                do i = 1, nr
                    ii = i + ir - 1
                    if (mod(ii - 1, nproc) == ip) then
                        sndbuf(il, j) = b(i, j)
                        il = il + 1
                    end if
                end do
            end do
            call mpi_barrier(comm_a, ierr)
            call mpi_gather(sndbuf, nbuf, mpi_double_complex, &
                            rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr)
            if (ierr /= 0) &
                call errore(" blk2cyc_zredist ", " in mpi_gather ", abs(ierr))
        end do
        !

        do ip = 0, nproc - 1
            !
            if (ip_desc(lambda_node_, ip + 1) > 0) then

                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
                !
                do j = 1, ip_nc
                    jj = j + ip_ic - 1
                    il = 1
                    do i = 1, ip_nr
                        ii = i + ip_ir - 1
                        if (mod(ii - 1, nproc) == me) then
                            a((ii - 1)/nproc + 1, jj) = rcvbuf(il, j, ip + 1)
                            il = il + 1
                        end if
                    end do
                end do

            end if

        end do
        !
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else

        a(1:n, 1:n) = b(1:n, 1:n)

#endif

        return

    end subroutine blk2cyc_zredist

end module descriptors
