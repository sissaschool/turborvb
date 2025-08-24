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

!> @brief Matrix descriptor module for parallel block-cyclic distribution
!>
!> This module provides data structures and utility routines for describing
!> and managing block-cyclic distributed matrices on a 2D processor grid.
!> It is designed for use with Cannon's algorithm and other parallel matrix
!> operations in scientific computing. The module supports both real and
!> complex matrices, and provides routines for descriptor initialization,
!> local/global index mapping, and parallel redistribution.
!>
!> @details
!> The descriptor array encodes the local/global structure of a distributed
!> matrix, including block sizes, processor grid coordinates, and MPI
!> communicator information. This enables efficient parallel matrix
!> multiplication and redistribution, especially in quantum chemistry and
!> electronic structure calculations.
!>
!> @note This module is used throughout TurboRVB for distributed linear algebra.
module descriptors
    !
    implicit none
    save

    !> @brief External index and block size mapping functions
    integer ldim_cyclic, ldim_block_sca
    integer lind_block_sca
    integer gind_block_sca
    external ldim_cyclic, ldim_block_sca
    external lind_block_sca
    external gind_block_sca

    !> @brief Descriptor array parameter indices for block-cyclic matrices
    integer, parameter :: descla_siz_ = 16  !< Size of descriptor array
    integer, parameter :: ilar_ = 1        !< Global index of first local row
    integer, parameter :: nlar_ = 2        !< Number of local rows
    integer, parameter :: ilac_ = 3        !< Global index of first local column
    integer, parameter :: nlac_ = 4        !< Number of local columns
    integer, parameter :: nlax_ = 5        !< Leading dimension of distributed matrix
    integer, parameter :: lambda_node_ = 6 !< >0 if processor holds a block
    integer, parameter :: la_n_ = 7        !< Global matrix dimension
    integer, parameter :: la_nx_ = 8       !< Global leading dimension
    integer, parameter :: la_npr_ = 9      !< Number of processor rows
    integer, parameter :: la_npc_ = 10     !< Number of processor columns
    integer, parameter :: la_myr_ = 11     !< Processor row index
    integer, parameter :: la_myc_ = 12     !< Processor column index
    integer, parameter :: la_comm_ = 13    !< MPI communicator
    integer, parameter :: la_me_ = 14      !< Processor linear index
    integer, parameter :: la_nrl_ = 15     !< Local rows for cyclic distribution
    integer, parameter :: la_nrlx_ = 16    !< Leading dimension for row distribution
    !
    !> @brief Descriptor array for a distributed matrix
    integer :: descla(descla_siz_)

contains

    !------------------------------------------------------------------------
    !> @brief Compute local block indices and sizes for block-cyclic distribution
    !>
    !> Determines the global index of the first local element and the number of
    !> local elements for a given processor in a block-cyclic distributed array.
    !>
    !> Parameters
    !> ----------
    !> i2g : integer, out
    !>     Global index of the first local element.
    !> nl : integer, out
    !>     Number of local elements.
    !> n : integer, in
    !>     Number of actual elements in the global array.
    !> nx : integer, in
    !>     Dimension of the global array (nx >= n) to be distributed.
    !> np : integer, in
    !>     Number of processors.
    !> me : integer, in
    !>     Task ID for which i2g and nl are computed.
    !>
    !> Notes
    !> -----
    !> - Allows distributing a global array larger than the number of actual elements.
    !> - Ensures equal partitioning for matrices of different sizes (e.g., spin-up/down).
    !>
    !> Example
    !> -------
    !> Used internally by descla_init to set up block-cyclic descriptors.
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
    !------------------------------------------------------------------------
    !> @brief Initialize a block-cyclic matrix descriptor for parallel distribution
    !>
    !> Sets up the descriptor array for a block-cyclic distributed matrix on a
    !> square processor grid. Handles both local and global matrix properties,
    !> including block sizes, processor coordinates, and MPI communicator info.
    !>
    !> Parameters
    !> ----------
    !> desc : integer array, out
    !>     Descriptor array to be initialized.
    !> n : integer, in
    !>     Size of the matrix (number of rows/columns).
    !> nx : integer, in
    !>     Maximum size among matrices sharing this descriptor.
    !> np : integer array(2), in
    !>     Number of processors in each grid dimension (must be square).
    !> me : integer array(2), in
    !>     Processor coordinates in the grid.
    !> comm : integer, in
    !>     MPI communicator.
    !> includeme : integer, in
    !>     If 1, include this processor in the distribution; else, set as inactive.
    !>
    !> Notes
    !> -----
    !> - Only square processor grids are supported.
    !> - Handles both block and cyclic distributions for advanced parallelism.
    !> - Performs error checking on all input parameters and computed values.
    !>
    !> Example
    !> -------
    !> Used to initialize matrix descriptors for distributed matrix-matrix multiplication.
    subroutine descla_init(desc, n, nx, np, me, comm, includeme)
        implicit none
        integer, intent(OUT) :: desc(:)
        integer, intent(IN) :: n !  the size of this matrix
        integer, intent(IN) :: nx !  the max among different matrices sharing
        !  this descriptor or the same data distribution
        integer, intent(IN) :: np(2), me(2), comm
        integer, intent(IN) :: includeme
        integer :: ir, nr, ic, nc, lnode, nlax, nrl, nrlx
        integer :: ip, npp

        !> @brief Error checking for processor grid and matrix sizes
        if (np(1) /= np(2)) &
            call errore(' descla_init ', ' only square grid of proc are allowed ', 2)
        if (n < 0) &
            call errore(' descla_init ', ' dummy argument n less than 1 ', 3)
        if (nx < n) &
            call errore(' descla_init ', ' dummy argument nx less than n ', 4)
        if (np(1) < 1) &
            call errore(' descla_init ', ' dummy argument np less than 1 ', 5)

        !> @brief Find the block maximum dimensions
        nlax = ldim_block_sca(nx, np(1), 0)
        !
        !> @brief Find local block dimensions, if appropriate
        !> Only for processes involved in the distribution
        if (includeme == 1) then
            call descla_local_dims(ir, nr, n, nx, np(1), me(1))
            call descla_local_dims(ic, nc, n, nx, np(2), me(2))
            lnode = 1
        else
            nr = 0
            nc = 0
            ir = 0
            ic = 0
            lnode = -1
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

        !> @brief Compute local dimension of the cyclically distributed matrix
        if (includeme == 1) then
            nrl = ldim_cyclic(n, npp, desc(la_me_))
        else
            nrl = 0
        end if
        nrlx = n/npp + 1

        desc(la_nrl_) = nrl
        desc(la_nrlx_) = nrlx

        !> @brief Error checking for computed values
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

        return
    end subroutine descla_init

#ifdef __SCALAPACK

    !------------------------------------------------------------------------
    !> @brief Symmetrize a distributed real square matrix in parallel
    !>
    !> This subroutine enforces symmetry (A = A^T) on a real square matrix
    !> distributed in block-cyclic fashion across a 2D processor grid. It uses
    !> MPI communication to exchange off-diagonal blocks and ensures that the
    !> resulting matrix is symmetric on all processors.
    !>
    !> Parameters
    !> ----------
    !> n : integer, in
    !>     Global matrix dimension (number of rows/columns).
    !> a : real*8 array, inout
    !>     Local block of the distributed matrix (size lda × *).
    !> lda : integer, in
    !>     Leading dimension of the local block.
    !> desc : integer array, in
    !>     Descriptor array describing the matrix distribution.
    !>
    !> Notes
    !> -----
    !> - For diagonal blocks, symmetry is enforced locally.
    !> - For off-diagonal blocks, MPI is used to exchange and transpose blocks.
    !> - Only active processors (desc(lambda_node_) > 0) participate.
    !> - Error checking is performed for matrix dimensions and leading dimensions.
    !>
    !> Example
    !> -------
    !> Used in distributed matrix-matrix multiplication to ensure symmetric results.
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

        !> @brief Only active processors participate
        if (desc(lambda_node_) <= 0) then
            return
        end if

        !> @brief Error checking for matrix dimensions
        if (n /= desc(la_n_)) &
            call errore(" dsqmsym ", " wrong global dim n ", n)
        if (lda /= desc(nlax_)) &
            call errore(" dsqmsym ", " wrong leading dim lda ", lda)

        comm = desc(la_comm_)

        nr = desc(nlar_)
        nc = desc(nlac_)
        !> @brief Diagonal block: enforce symmetry locally
        if (desc(la_myc_) == desc(la_myr_)) then
            do j = 1, nc
                do i = j + 1, nr
                    a(i, j) = a(j, i)
                end do
            end do
        !> @brief Super-diagonal block: send block to sub-diagonal processor
        else if (desc(la_myc_) > desc(la_myr_)) then
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), dest)
            call mpi_isend(a, lda*lda, MPI_DOUBLE_PRECISION, dest, 1, comm, sreq, ierr)
            if (ierr /= 0) &
                call errore(" dsqmsym ", " in isend ", abs(ierr))
        !> @brief Sub-diagonal block: receive block, transpose locally
        else if (desc(la_myc_) < desc(la_myr_)) then
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), sour)
            call mpi_recv(a, lda*lda, MPI_DOUBLE_PRECISION, sour, 1, comm, istatus, ierr)
            if (ierr /= 0) &
                call errore(" dsqmsym ", " in recv ", abs(ierr))
            do j = 1, lda
                do i = j + 1, lda
                    atmp = a(i, j)
                    a(i, j) = a(j, i)
                    a(j, i) = atmp
                end do
            end do
        end if

        if (desc(la_myc_) > desc(la_myr_)) then
            call MPI_Wait(sreq, istatus, ierr)
            if (ierr /= 0) &
                call errore(" dsqmsym ", " in wait ", abs(ierr))
        end if

#else
        !> @brief Serial: enforce symmetry locally for full matrix
        do j = 1, n
            do i = j + 1, n
                a(i, j) = a(j, i)
            end do
        end do
#endif

        return
    end subroutine dsqmsym

    !------------------------------------------------------------------------
    !> @brief Hermitianize a distributed complex square matrix in parallel
    !>
    !> This subroutine enforces Hermitian symmetry (A = A^†) on a complex
    !> square matrix distributed in block-cyclic fashion across a 2D processor grid.
    !> It uses MPI communication to exchange off-diagonal blocks and ensures that
    !> the resulting matrix is Hermitian on all processors.
    !>
    !> Parameters
    !> ----------
    !> n : integer, in
    !>     Global matrix dimension (number of rows/columns).
    !> a : complex*16 array, inout
    !>     Local block of the distributed matrix (size lda × lda).
    !> lda : integer, in
    !>     Leading dimension of the local block.
    !> desc : integer array, in
    !>     Descriptor array describing the matrix distribution.
    !>
    !> Notes
    !> -----
    !> - For diagonal blocks, Hermitian symmetry is enforced locally.
    !> - For off-diagonal blocks, MPI is used to exchange and conjugate-transpose blocks.
    !> - Only active processors (desc(lambda_node_) > 0) participate.
    !> - Error checking is performed for matrix dimensions and leading dimensions.
    !>
    !> Example
    !> -------
    !> Used in distributed matrix-matrix multiplication to ensure Hermitian results.
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

        !> @brief Only active processors participate
        if (desc(lambda_node_) <= 0) then
            return
        end if

        !> @brief Error checking for matrix dimensions
        if (n /= desc(la_n_)) &
            call errore(" zsqmher ", " wrong global dim n ", n)
        if (lda /= desc(nlax_)) &
            call errore(" zsqmher ", " wrong leading dim lda ", lda)

        comm = desc(la_comm_)

        nr = desc(nlar_)
        nc = desc(nlac_)
        !> @brief Diagonal block: enforce Hermitian symmetry locally
        if (desc(la_myc_) == desc(la_myr_)) then
            do j = 1, nc
                a(j, j) = DCMPLX(dble(a(j, j)))
                do i = j + 1, nr
                    a(i, j) = conjg(a(j, i))
                end do
            end do
        !> @brief Super-diagonal block: send block to sub-diagonal processor
        else if (desc(la_myc_) > desc(la_myr_)) then
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), dest)
            call mpi_isend(a, lda*lda, MPI_DOUBLE_COMPLEX, dest, 1, comm, sreq, ierr)
            if (ierr /= 0) &
                call errore(" zsqmher ", " in mpi_isend ", abs(ierr))
        !> @brief Sub-diagonal block: receive block, transpose and conjugate locally
        else if (desc(la_myc_) > desc(la_myr_)) then
            call GRID2D_RANK('R', desc(la_npr_), desc(la_npc_), &
                             desc(la_myc_), desc(la_myr_), sour)
            call mpi_recv(a, lda*lda, MPI_DOUBLE_COMPLEX, sour, 1, comm, istatus, ierr)
            if (ierr /= 0) &
                call errore(" zsqmher ", " in mpi_recv ", abs(ierr))
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
        end if

        if (desc(la_myc_) > desc(la_myr_)) then
            call MPI_Wait(sreq, istatus, ierr)
            if (ierr /= 0) &
                call errore(" zsqmher ", " in MPI_Wait ", abs(ierr))
        end if

#else
        !> @brief Serial: enforce Hermitian symmetry locally for full matrix
        do j = 1, n
            a(j, j) = DCMPLX(dble(a(j, j)))
            do i = j + 1, n
                a(i, j) = conjg(a(j, i))
            end do
        end do
#endif

        return
    end subroutine zsqmher

#endif
    !

    !------------------------------------------------------------------------
    !> @brief Redistribute a real matrix from cyclic to block-cyclic distribution
    !>
    !> This subroutine redistributes a real matrix A, initially distributed
    !> cyclically by rows across processors, into a block-cyclic distributed
    !> matrix B on a 2D processor grid. It uses MPI communication to gather
    !> and scatter matrix blocks as needed.
    !>
    !> Parameters
    !> ----------
    !> n : integer, in
    !>     Global matrix dimension (number of rows/columns).
    !> a : real*8 array, in
    !>     Input matrix, cyclically distributed by rows (size lda × *).
    !> lda : integer, in
    !>     Leading dimension of the input matrix.
    !> b : real*8 array, out
    !>     Output matrix, block-cyclic distributed (size ldb × *).
    !> ldb : integer, in
    !>     Leading dimension of the output matrix.
    !> desc : integer array, in
    !>     Descriptor array describing the block-cyclic distribution.
    !>
    !> Notes
    !> -----
    !> - Only square processor grids are supported.
    !> - Uses MPI_Allgather and MPI_Gather for communication.
    !> - Performs error checking for block sizes and processor mesh.
    !> - Serial fallback simply copies the matrix.
    !>
    !> Example
    !> -------
    !> Used to convert between different distributed matrix layouts in parallel algorithms.
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

        !> @brief Only active processors participate
        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        !> @brief Error checking for processor mesh and matrix size
        if (np /= desc(la_npc_)) &
            call errore(' cyc2blk_redist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' cyc2blk_redist ', ' n less or equal zero ', 1)

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

        !> @brief Loop over all processors to gather and redistribute blocks
        do ip = 0, nproc - 1
            if (ip_desc(nlax_, ip + 1) /= nb) &
                call errore(' cyc2blk_redist ', ' inconsistent block dim nb ', 1)
            if (ip_desc(lambda_node_, ip + 1) > 0) then
                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
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

        nr = desc(nlar_)
        nc = desc(nlac_)
        ir = desc(ilar_)
        ic = desc(ilac_)

        !> @brief Unpack received blocks into the output matrix
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
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else
        !> @brief Serial: copy input matrix to output
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

    !------------------------------------------------------------------------
    !> @brief Redistribute a complex matrix from cyclic to block-cyclic distribution
    !>
    !> This subroutine redistributes a complex matrix A, initially distributed
    !> cyclically by rows across processors, into a block-cyclic distributed
    !> matrix B on a 2D processor grid. It uses MPI communication to gather
    !> and scatter matrix blocks as needed.
    !>
    !> Parameters
    !> ----------
    !> n : integer, in
    !>     Global matrix dimension (number of rows/columns).
    !> a : complex*16 array, in
    !>     Input matrix, cyclically distributed by rows (size lda × *).
    !> lda : integer, in
    !>     Leading dimension of the input matrix.
    !> b : complex*16 array, out
    !>     Output matrix, block-cyclic distributed (size ldb × *).
    !> ldb : integer, in
    !>     Leading dimension of the output matrix.
    !> desc : integer array, in
    !>     Descriptor array describing the block-cyclic distribution.
    !>
    !> Notes
    !> -----
    !> - Only square processor grids are supported.
    !> - Uses MPI_Allgather and MPI_Gather for communication.
    !> - Performs error checking for block sizes and processor mesh.
    !> - Serial fallback simply copies the matrix.
    !>
    !> Example
    !> -------
    !> Used to convert between different distributed matrix layouts in parallel algorithms.
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

        !> @brief Only active processors participate
        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        !> @brief Error checking for processor mesh and matrix size
        if (np /= desc(la_npc_)) &
            call errore(' cyc2blk_zredist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' cyc2blk_zredist ', ' n less or equal zero ', 1)

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

        !> @brief Loop over all processors to gather and redistribute blocks
        do ip = 0, nproc - 1
            if (ip_desc(nlax_, ip + 1) /= nb) &
                call errore(' cyc2blk_zredist ', ' inconsistent block dim nb ', 1)
            if (ip_desc(lambda_node_, ip + 1) > 0) then
                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
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

        nr = desc(nlar_)
        nc = desc(nlac_)
        ir = desc(ilar_)
        ic = desc(ilac_)

        !> @brief Unpack received blocks into the output matrix
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
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else
        !> @brief Serial: copy input matrix to output
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

    !------------------------------------------------------------------------
    !> @brief Redistribute a real matrix from block-cyclic to cyclic distribution
    !>
    !> This subroutine redistributes a real matrix B, initially distributed
    !> in block-cyclic fashion on a 2D processor grid, into a cyclically
    !> distributed matrix A by rows across processors. It uses MPI communication
    !> to gather and scatter matrix blocks as needed.
    !>
    !> Parameters
    !> ----------
    !> n : integer, in
    !>     Global matrix dimension (number of rows/columns).
    !> a : real*8 array, out
    !>     Output matrix, cyclically distributed by rows (size lda × *).
    !> lda : integer, in
    !>     Leading dimension of the output matrix.
    !> b : real*8 array, in
    !>     Input matrix, block-cyclic distributed (size ldb × *).
    !> ldb : integer, in
    !>     Leading dimension of the input matrix.
    !> desc : integer array, in
    !>     Descriptor array describing the block-cyclic distribution.
    !>
    !> Notes
    !> -----
    !> - Only square processor grids are supported.
    !> - Uses MPI_Allgather and MPI_Gather for communication.
    !> - Performs error checking for block sizes and processor mesh.
    !> - Serial fallback simply copies the matrix.
    !>
    !> Example
    !> -------
    !> Used to convert from block-cyclic to cyclic distribution for certain algorithms.
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

        !> @brief Only active processors participate
        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        !> @brief Error checking for processor mesh and matrix size
        if (np /= desc(la_npc_)) &
            call errore(' blk2cyc_redist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' blk2cyc_redist ', ' n less or equal zero ', 1)

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
        !> @brief Pack local blocks into send buffer for redistribution
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

        !> @brief Unpack received blocks into cyclic distribution
        do ip = 0, nproc - 1
            if (ip_desc(lambda_node_, ip + 1) > 0) then
                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
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
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else
        !> @brief Serial: copy input matrix to output
        a(1:n, 1:n) = b(1:n, 1:n)
#endif

        return
    end subroutine blk2cyc_redist

    !------------------------------------------------------------------------
    !> @brief Redistribute a complex matrix from block-cyclic to cyclic distribution
    !>
    !> This subroutine redistributes a complex matrix B, initially distributed
    !> in block-cyclic fashion on a 2D processor grid, into a cyclically
    !> distributed matrix A by rows across processors. It uses MPI communication
    !> to gather and scatter matrix blocks as needed.
    !>
    !> Parameters
    !> ----------
    !> n : integer, in
    !>     Global matrix dimension (number of rows/columns).
    !> a : complex*16 array, out
    !>     Output matrix, cyclically distributed by rows (size lda × *).
    !> lda : integer, in
    !>     Leading dimension of the output matrix.
    !> b : complex*16 array, in
    !>     Input matrix, block-cyclic distributed (size ldb × *).
    !> ldb : integer, in
    !>     Leading dimension of the input matrix.
    !> desc : integer array, in
    !>     Descriptor array describing the block-cyclic distribution.
    !>
    !> Notes
    !> -----
    !> - Only square processor grids are supported.
    !> - Uses MPI_Allgather and MPI_Gather for communication.
    !> - Performs error checking for block sizes and processor mesh.
    !> - Serial fallback simply copies the matrix.
    !>
    !> Example
    !> -------
    !> Used to convert from block-cyclic to cyclic distribution for certain algorithms.
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

        !> @brief Only active processors participate
        if (desc(lambda_node_) < 0) then
            return
        end if

        np = desc(la_npr_) !  dimension of the processor mesh
        nb = desc(nlax_) !  leading dimension of the local matrix block
        me = desc(la_me_) !  my processor id (starting from 0)
        comm_a = desc(la_comm_)
        nproc = desc(la_npr_)*desc(la_npc_)

        !> @brief Error checking for processor mesh and matrix size
        if (np /= desc(la_npc_)) &
            call errore(' blk2cyc_zredist ', ' works only with square processor mesh ', 1)
        if (n < 1) &
            call errore(' blk2cyc_zredist ', ' n less or equal zero ', 1)

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
        !> @brief Pack local blocks into send buffer for redistribution
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

        !> @brief Unpack received blocks into cyclic distribution
        do ip = 0, nproc - 1
            if (ip_desc(lambda_node_, ip + 1) > 0) then
                ip_nr = ip_desc(nlar_, ip + 1)
                ip_nc = ip_desc(nlac_, ip + 1)
                ip_ir = ip_desc(ilar_, ip + 1)
                ip_ic = ip_desc(ilac_, ip + 1)
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
        deallocate (ip_desc)
        deallocate (rcvbuf)
        deallocate (sndbuf)

#else
        !> @brief Serial: copy input matrix to output
        a(1:n, 1:n) = b(1:n, 1:n)
#endif

        return
    end subroutine blk2cyc_zredist

end module descriptors
