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

module parallel_module

    use allio, only: nproc, rank, nprocrep, rankrep, nproccolrep, rankcolrep, &
                     commrep_mpi, commcolrep_mpi, manyfort10, iflagerr, nbead, yeswrite10, min_block, writescratch
    use kpoints_mod, only: nk, wkp, kaverage, compute_bands
    use constants, only: old_threads

    implicit none

    type pool
        integer :: comm
        integer :: rank
        integer :: nproc
    end type pool

    integer, private :: ierr
    type(pool), public :: intra_pool, inter_pool

contains

    subroutine setup_para

        ! This subroutine setup the parallel environments involved
        ! in the DFT calculation. On top there is the splitting of
        ! MPI communicators in rows/columns to deal with k-points
        ! sampling:
        ! commrep_mpi --> communicator inside a processor row which deals with a
        !                 single k-point during the calculation.
        ! commcolrep_mpi --> communicator among different rows (columns communicator). It
        !                    has to be used in order to sum certain quantities among k-points
        !                    such as density on the grid of DFT energy.
        ! If specified at compilation with the option -D__SCALAPACK, inside a processor row the main
        ! matrices are divided in square blocks to speed up the diagonalization with SCALAPACK
        ! routines. Also SCALAPACK environment is initialized in this subroutine.

#if defined __SCALAPACK
        use allio, only: me_blacs, np_blacs, world_cntx, np_ortho1, np_ortho, ortho_cntx, &
                         leg_ortho, ortho_comm, me_ortho, me_ortho1, ortho_comm_id, nelorb
#endif

        implicit none

#if defined __SCALAPACK
        integer, allocatable :: blacsmap(:, :)
        integer, allocatable :: blacsmap_other(:, :)
        integer, allocatable :: blacstmp(:, :)
        integer :: nprow, npcol, myrow, mycol, color, key
#endif
        integer :: row_id, col_id, row_comm, col_comm, mcol, irow, &
                   jcol, nk_opt, mcol_rep
        logical :: split_comms
#ifdef PARALLEL
        include 'mpif.h'
#endif

        iflagerr = 0

        ! If kaverage=.true. then the calculation involves a k-points sampling (nk.gt.1 and SC calculation).
        ! When sampling is done we have:
        ! 1) initialization of parallel environment for splitting the k-points among pools
        ! 2) different printing of variables during the calculation
        ! 3) possibility to read/write multiple fort.10 for each k-point
        ! 4) Fermi energy and smeared occupations calculated within k-points framework
        ! If band structure calculation, communicators must not be splitted.
        split_comms = .false.
        if (manyfort10 .and. .not. compute_bands) split_comms = .true.
        !
        ! Default values of communicators/id
#if defined  PARALLEL
        ! row = comm for processors dealing with a single k-point
        row_comm = MPI_COMM_WORLD
        commrep_mpi = MPI_COMM_WORLD
        rankrep = rank
        nprocrep = nproc
        ! column = comm among k-points
        col_comm = MPI_COMM_WORLD
        commcolrep_mpi = MPI_COMM_WORLD
        rankcolrep = 0
        nproccolrep = 1
        mcol = nproc
#else
        row_comm = 0
        commrep_mpi = 0
        rankrep = rank
        nprocrep = nproc
        !
        col_comm = 0
        commcolrep_mpi = 0
        rankcolrep = rank
        nproccolrep = nproc
        mcol = nproc
#endif

        ! trivial check
#ifndef PARALLEL
        if (split_comms) &
            call error(' setup_parallel ', &
                       ' compile with PARALLEL option for more than one k-point ', 1, rank)
#endif

#if defined PARALLEL

        ! split communicators only if dealing with more than 1 k-point
        if (split_comms) then

            if (nk .gt. nproc) then
                if (rank .eq. 0) write (6, *) ' # MPI processes/ # of k-points ', nproc, nk
                call error(' setup_parallel ', ' choose an higher number of processors &
                        &or decrease the number of k-points! ', 1, rank)
            end if

            if (mod(nproc, nk) .ne. 0) then
                if (rank .eq. 0) write (6, *) ' # MPI processes/ # of k-points ', nproc, nk
                call error(' setup_parallel ', ' Choose a number of processors multiple of the &
                        &number of k-points! ', 1, rank)
            end if
            ! OLD version. Now stop the program if nk is not multiple of nproc.
            !           ! force number of k-points to be multiple of nproc if this is not the case
            !        else
            !           call find_nk_opt(nproc,nk,nk_opt)
            !           if(rank.eq.0) write(6,*) 'Warning: # of k-points rescaled from',nk,'to',nk_opt
            !           nk = nk_opt
            !        endif

            call mpi_barrier(MPI_COMM_WORLD, ierr)

#if defined DEBUG
            if (rank .eq. 0) write (*, *) 'Processors involved:', nproc, nproc/nk
            if (rank .eq. 0) write (*, *) '    Iam    irow    jcol  row-id  col-id'
#endif

            row_id = 0
            col_id = rank
            mcol = nproc/nk
            irow = rank/mcol
            jcol = mod(rank, mcol)
            ! split the communicators
            call MPI_comm_split(MPI_COMM_WORLD, irow, jcol, row_comm, ierr)
            call MPI_comm_split(MPI_COMM_WORLD, jcol, irow, col_comm, ierr)
            call MPI_comm_rank(row_comm, row_id, ierr)
            call MPI_comm_rank(col_comm, col_id, ierr)
            call MPI_barrier(MPI_COMM_WORLD, ierr) ! for unreliable networks
            ! intra-pools
            intra_pool%comm = row_comm
            intra_pool%rank = row_id
            intra_pool%nproc = mcol
            ! inter-pools
            inter_pool%comm = col_comm
            inter_pool%rank = col_id
            inter_pool%nproc = nk

            nproccolrep = nk
            mcol_rep = mcol
            commrep_mpi = row_comm
            nprocrep = mcol
            commcolrep_mpi = col_comm
            rankrep = row_id
            rankcolrep = col_id

            ! check for errors in splitting communicators
            if (ierr .ne. 0) &
                call error(' InitializeAll ', &
                           ' ERROR in splitting communicators for k-points sampling!', 1, rank)

            ! print processors grid for debugging
#ifdef DEBUG
            write (*, '(5I5)') rank, irow, jcol, rankrep, rankcolrep
#endif
        end if

#endif
        !
        ! setup SCALAPACK environment
        !
        ! in the case of k-points sampling the SCALAPACK grid is
        ! initialized within each processors pool.
        !
        ! SCALAPACK version needs the MPI environment.
#ifndef PARALLEL
#ifdef __SCALAPACK
        call error(' InitializeAll ', &
                   ' compile with -DPARALLEL option if you want to use SCALAPACK version', 1, rank)
#endif
#endif

#if defined __SCALAPACK && defined PARALLEL

        call BLACS_PINFO(me_blacs, np_blacs)
        ! me_blacs = BLACS process identifier
        call BLACS_GET(-1, 0, world_cntx)
        ! This subroutine factorizes the number of processors (NPROC)
        ! into NPROW and NPCOL according to the shape
        ! 'S' = square grid
        ! np_ortho(1) = # of row of the grid
        ! np_ortho(2) = # of column of the grid
        call grid2d_dims('S', nprocrep, np_ortho(1), np_ortho(2))

        !   propagate nelorb in any case common to all processors
        call bcast_integer(nelorb, 1, 0, mpi_comm_world)
        if (nelorb/np_ortho(1) .lt. min_block) then
            np_ortho(1) = nelorb/min_block
            if (np_ortho(1)*min_block .ne. nelorb) np_ortho(1) = np_ortho(1) + 1
            np_ortho(2) = np_ortho(1)
            if (rank .eq. 0) then
                write (6, *) ' Warning using less number of processor to distribute the matrix , block too small !!! ' &
                    , np_ortho(1), nelorb
            end if
        end if

        np_ortho1 = np_ortho(1)*np_ortho(2)

        !  here we choose the first processors omp threads can be set anyway
        !
        color = 0
        if (rankrep < np_ortho1) color = 1
        !
        leg_ortho = 1
        !
        !
        key = rankrep
        !
        !  initialize the communicator for the new group
        !  color = 1 : the group of processors performs matrix distribution
        !  color = 0 : processors not involved in the distribution
        !
#ifdef PARALLEL
        call MPI_COMM_SPLIT(commrep_mpi, color, key, ortho_comm, ierr)
#endif
        if (ierr /= 0) &
            call errore(" init_ortho_group ", " error splitting communicator ", ierr)
        !
        !  Computes coordinates of the processors, in row maior order
        !
#ifdef PARALLEL
        call mpi_comm_size(ortho_comm, np_ortho1, ierr)
        call mpi_comm_rank(ortho_comm, me_ortho1, ierr)
#endif
        ! np_ortho1 = size of the group associated with ortho_comm
        ! me_ortho1 = rank of the current MPI task
        if (color == 1 .and. np_ortho1 /= np_ortho(1)*np_ortho(2)) &
            call errore(" init_ortho_group ", " wrong number of proc in ortho_comm ", ierr)
        !
        if (rankrep == 0 .and. me_ortho1 /= 0) &
            call errore(" init_ortho_group ", " wrong root in ortho_comm ", ierr)
        !
        if (color == 1) then ! only for processes involved
            ortho_comm_id = 1
            call GRID2D_COORDS('R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2))
            ! me_ortho(1) = row of the MPI task in the grid
            ! me_ortho(2) = column of the MPI task
            call GRID2D_RANK('R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr)
            if (ierr /= me_ortho1) &
                call errore(" init_ortho_group ", " wrong coordinates in ortho_comm ", ierr)
            if (me_ortho1*leg_ortho /= rankrep) &
                call errore(" init_ortho_group ", " wrong rank assignment in ortho_comm ", ierr)
        else
            ortho_comm_id = 0
            me_ortho(1) = me_ortho1
            me_ortho(2) = me_ortho1
        end if

        if (ortho_comm_id > 0) then
            allocate (blacsmap(np_ortho(1), np_ortho(2)))
            allocate (blacstmp(np_ortho(1), np_ortho(2)))
            blacstmp = 0
            blacsmap = 0
            blacsmap(me_ortho(1) + 1, me_ortho(2) + 1) = me_blacs
            nprow = np_ortho(1)
            npcol = np_ortho(2)
        else ! processes not involved in distribution
            nprow = np_ortho1
            npcol = 1
            allocate (blacsmap(np_ortho1, 1))
            allocate (blacstmp(np_ortho1, 1))
            blacstmp = 0
            blacsmap = 0
            blacsmap(me_ortho1 + 1, 1) = me_blacs
        end if

#ifdef PARALLEL
        call mpi_allreduce(blacsmap, blacstmp, size(blacsmap), MPI_INTEGER, MPI_SUM, ortho_comm, ierr)
#endif
        blacsmap = blacstmp

        !  map processors to grid
        ortho_cntx = world_cntx
        call BLACS_GRIDMAP(ortho_cntx, blacsmap, nprow, nprow, npcol)

        !  get info from communicator ortho_cntx and check the values
        call BLACS_GRIDINFO(ortho_cntx, nprow, npcol, myrow, mycol)
        if (ortho_comm_id > 0) then
            if (np_ortho(1) /= nprow) call errore(' init ', ' problem with SCALAPACK, wrong nprow ', 1)
            if (np_ortho(2) /= npcol) call errore(' init ', ' problem with SCALAPACK, wrong npcol ', 1)
            if (me_ortho(1) /= myrow) call errore(' init ', ' problem with SCALAPACK, wrong myrow ', 1)
            if (me_ortho(2) /= mycol) call errore(' init ', ' problem with SCALAPACK, wrong mycol ', 1)
        end if

        deallocate (blacsmap)
        deallocate (blacstmp)

#endif

        ! final check for error flags
        call checkiflagerr(iflagerr, rank, "ERROR in initialization parallel environment")

        return

    end subroutine setup_para

    ! This subroutine check if the number of k-points in input
    ! is multiple of the processors requested and rescale this
    ! number otherwise. In order to avoid this rescaling, use
    ! the tool find_kpoints.x in order to know how many processor
    ! need to be allocated for the DFT run.

    subroutine find_nk_opt(nproc, nk, nk_opt)

        implicit none
        integer i, nproc, nk, nk_opt, diff

        nk_opt = nk
        diff = 0
        if (mod(nproc, nk) .eq. 0) return
        do i = 1, nproc
            if (mod(nproc, i) .eq. 0) then
                diff = nk - nproc/i
                if (diff .gt. 0) then
                    nk_opt = nproc/i
                    exit
                end if
            end if
        end do

        return

    end subroutine find_nk_opt

    ! -------------------------------------------------
    ! Set of routines to sum/average/collect quantities
    ! which are distributed across processor pools
    ! -------------------------------------------------

    ! This routine recovers eigenvalues and occupations of KS orbitals
    ! from the processor pools and collect then to the master in the variables
    ! occupations_sav and eigmol_sav.
    ! It is used in the evaluation of electron occupations and for variables
    ! printout at the end of SC cycle.

    subroutine collect_from_pools

        use constants, only: ipc
        use setup, only: yeslsda, bands, indk, eigmol, eigmoldo, occupations, occupationdo, eigmol_sav, &
                         eigmoldo_sav, occupations_sav, occupationsdo_sav, molecorb, molecorbdo, &
                         molecorb_sav, molecorbdo_sav
        implicit none
        ! local
        integer :: i
        logical :: double_occ

#if defined PARALLEL
        include "mpif.h"
#endif

        if (compute_bands) return

        double_occ = .false.
        if (yeslsda .or. ipc .eq. 2) double_occ = .true.

        eigmol_sav = 0.d0
        occupations_sav = 0.d0
        if (writescratch .eq. 0) molecorb_sav = 0.d0
        if (double_occ) then
            eigmoldo_sav = 0.d0
            occupationsdo_sav = 0.d0
            if (writescratch .eq. 0) molecorbdo_sav = 0.d0
        end if

        if (.not. kaverage .or. (kaverage .and. nk .eq. 1)) then ! therefore I have just one k-point
            if (rank .eq. 0) then
                eigmol_sav(1:bands, 1) = eigmol(1:bands)
                occupations_sav(1:bands, 1) = occupations(1:bands)
                if (writescratch .eq. 0) molecorb_sav(:, 1:bands, 1) = molecorb(:, 1:bands)
                if (double_occ) then
                    eigmoldo_sav(1:bands, 1) = eigmoldo(1:bands)
                    occupationsdo_sav(1:bands, 1) = occupationdo(1:bands)
                    if (writescratch .eq. 0) molecorbdo_sav(:, 1, 1) = molecorbdo(:, 1)
                end if
            end if

        else

#ifdef PARALLEL
            call mpi_gather(eigmol, bands, MPI_DOUBLE_PRECISION, eigmol_sav, &
                            bands, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
            call mpi_gather(occupations, bands, MPI_DOUBLE_PRECISION, occupations_sav, &
                            bands, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
            if (writescratch .eq. 0) &
        &   call mpi_gather(molecorb, size(molecorb), MPI_DOUBLE_PRECISION, molecorb_sav, &
                            size(molecorb), MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
            if (double_occ) then
                call mpi_gather(eigmoldo, bands, MPI_DOUBLE_PRECISION, eigmoldo_sav, &
                                bands, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
                call mpi_gather(occupationdo, bands, MPI_DOUBLE_PRECISION, occupationsdo_sav, &
                                bands, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
                if (writescratch .eq. 0) &
             &  call mpi_gather(molecorbdo, size(molecorbdo), MPI_DOUBLE_PRECISION, molecorbdo_sav, &
                                size(molecorbdo), MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
            end if
#endif

        end if
        !
        ! syncronize
#if defined PARALLEL && defined UNREL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
        return

    end subroutine collect_from_pools

    subroutine gather_from_grid(ps_in, buf, ps_out, meshproc, nx, ny, nz)
        implicit none
        ! input
        integer, intent(in) :: nx, ny, nz, meshproc
        real(8), intent(in) :: ps_in(meshproc)
        real(8), intent(inout) :: ps_out(nx, ny, nz), buf(meshproc, nprocrep)
        ! local
        integer :: indmesh, indproc, proc, i, j, k, ierr

#ifdef PARALLEL
        include "mpif.h"
        buf = 0.d0
        call mpi_gather(ps_in, meshproc, MPI_DOUBLE_PRECISION, buf, &
                        meshproc, MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
#endif
        if (rankrep .eq. 0) then
            indproc = 0
            indmesh = 1
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        indproc = indproc + 1
                        ps_out(i, j, k) = buf(indmesh, indproc)
                        if (indproc .eq. nprocrep) then
                            indproc = 0
                            indmesh = indmesh + 1
                        end if
                    end do
                end do
            end do
        end if

        return
    end subroutine gather_from_grid

    subroutine scatter_to_grid(ps_in, buf, ps_out, meshproc, nx, ny, nz)
        implicit none
        ! input
        integer, intent(in) :: nx, ny, nz, meshproc
        real(8), intent(in) :: ps_in(nx, ny, nz)
        real(8), intent(inout) :: ps_out(meshproc), buf(meshproc, nprocrep)
        ! local
        integer :: proc, i, j, k, indmesh, indproc, ierr
#ifdef PARALLEL
        include "mpif.h"
#endif
        buf = 0.d0
        if (rankrep .eq. 0) then
            indmesh = 1
            indproc = 0
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        indproc = indproc + 1
                        buf(indmesh, indproc) = ps_in(i, j, k)
                        if (indproc .eq. nprocrep) then
                            indmesh = indmesh + 1
                            indproc = 0
                        end if
                    end do
                end do
            end do
        end if
#ifdef PARALLEL
        call mpi_scatter(buf, meshproc, MPI_DOUBLE_PRECISION, ps_out, meshproc, &
                         MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
#endif
        return
    end subroutine scatter_to_grid

    subroutine collect_kpoints(ps, dim_ps, value, root)
        ! simple summation over pools
        implicit none
        integer, intent(in) :: dim_ps, root
        real(8), intent(in) :: value
        real(8), intent(inout) :: ps(dim_ps)
        integer :: ind_ps, ierr
#ifdef PARALLEL
        include "mpif.h"
#endif
        if (.not. manyfort10) return
        ind_ps = rankcolrep + 1
        ps(ind_ps) = value
#ifdef PARALLEL
        call reduce_base_real(dim_ps, ps, commcolrep_mpi, root)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
        return
    end subroutine collect_kpoints

end module parallel_module
