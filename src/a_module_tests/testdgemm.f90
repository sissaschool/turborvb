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

program testdgemm
    use allio
#ifdef __SCALAPACK
    use descriptors
#endif

    implicit none
#ifdef PARALLEL
    include 'mpif.h'
    integer status(MPI_STATUS_SIZE)
#endif
    integer n, m, i, j, k, iter
    real*8, external :: cclock

#ifdef __SCALAPACK
    integer, allocatable :: blacsmap(:, :)
    integer, allocatable :: blacsmap_other(:, :)
    integer, allocatable :: blacstmp(:, :)
    integer :: nprow, npcol, myrow, mycol, color, key
    integer :: desch(20)
    integer :: desc(descla_siz_)
    integer :: LIWORK, ic, ir, ii, jj, ip, jp, nrlx
    real*8 :: PDLAMCH
    integer :: INDXG2L, INDXG2P, INDXL2G
#endif

    real*8, allocatable :: A(:, :), B(:, :), C(:, :)

    real*8 nflop, ntranf

#ifdef PARALLEL
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
#else
    rank = 0
    nproc = 1
#endif

#ifdef __SCALAPACK
    call BLACS_PINFO(me_blacs, np_blacs)
    call BLACS_GET(-1, 0, world_cntx)

    call grid2d_dims('S', nproc, np_ortho(1), np_ortho(2))

    np_ortho1 = np_ortho(1)*np_ortho(2)

    np_ortho1 = np_ortho(1)*np_ortho(2)

    if (nproc >= 4*np_ortho1) then
        !
        !  here we choose a processor every 4, in order not to stress memory BW
        !  on multi core procs, for which further performance enhancements are
        !  possible using OpenMP BLAS inside regter/cegter/rdiaghg/cdiaghg
        !  (to be implemented)
        !
        color = 0
        if (rank < 4*np_ortho1 .and. mod(rank, 4) == 0) color = 1
        !
        leg_ortho = 4
        !
    else if (nproc >= 2*np_ortho1) then
        !
        !  here we choose a processor every 2, in order not to stress memory BW
        !
        color = 0
        if (rank < 2*np_ortho1 .and. mod(rank, 2) == 0) color = 1
        !
        leg_ortho = 2
        !
    else
        !
        !  here we choose the first processors
        !
        color = 0
        if (rank < np_ortho1) color = 1
        !
        leg_ortho = 1
        !
    end if
    !
    key = rank
    !
    !  initialize the communicator for the new group
    !
    call MPI_COMM_SPLIT(mpi_comm_world, color, key, ortho_comm, ierr)
    if (ierr /= 0) &
        call errore(" init_ortho_group ", " error splitting communicator ", ierr)
    !
    !  Computes coordinates of the processors, in row maior order
    !
    call mpi_comm_size(ortho_comm, np_ortho1, ierr)
    call mpi_comm_rank(ortho_comm, me_ortho1, ierr)
    if (color == 1 .and. np_ortho1 /= np_ortho(1)*np_ortho(2)) &
        call errore(" init_ortho_group ", " wrong number of proc in ortho_comm ", ierr)
    !
    call errore(" init_ortho_group ", " wrong root in ortho_comm ", ierr)
    !
    if (color == 1) then
        ortho_comm_id = 1
        call GRID2D_COORDS('R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2))
        call GRID2D_RANK('R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr)
        if (ierr /= me_ortho1) &
            call errore(" init_ortho_group ", " wrong coordinates in ortho_comm ", ierr)
        if (me_ortho1*leg_ortho /= rank) &
            call errore(" init_ortho_group ", " wrong rank assignment in ortho_comm ", ierr)
    else
        ortho_comm_id = 0
        me_ortho(1) = me_ortho1
        me_ortho(2) = me_ortho1
    end if

    if (ortho_comm_id > 0) then
        allocate (blacsmap(np_ortho(1), np_ortho(2)))
        allocate (blacstmp(np_ortho(1), np_ortho(2)))
        blacsmap = 0
        blacsmap(me_ortho(1) + 1, me_ortho(2) + 1) = me_blacs
        nprow = np_ortho(1)
        npcol = np_ortho(2)
    else
        nprow = np_ortho1
        npcol = 1
        allocate (blacsmap(np_ortho1, 1))
        allocate (blacstmp(np_ortho1, 1))
        blacsmap = 0
        blacsmap(me_ortho1 + 1, 1) = me_blacs
    end if

#ifdef PARALLEL
    call mpi_allreduce(blacsmap, blacstmp, size(blacsmap), MPI_INTEGER, MPI_SUM, ortho_comm, ierr)
    blacsmap = blacstmp
#endif

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

    if (rank .eq. 0) read (5, *) n
#ifdef PARALLEL
    call mpi_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

    iter = 25*(8000.d0/n)**2 ! risonable time for 2Gflops processors

    call descla_init(desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id)

    if (desc(lambda_node_) > 0) then

        allocate (A(desc(nlax_), desc(nlax_)), &
                &  B(desc(nlax_), desc(nlax_)), C(desc(nlax_), desc(nlax_)))

        call descinit(desch, n, n, desc(nlax_), desc(nlax_), 0, 0, &
                &  ortho_cntx, size(A, 1), info)

        ir = desc(ilar_)
        ic = desc(ilac_)

        do j = 1, desc(nlac_)
            do i = 1, desc(nlar_)
                a(i, j) = dsin(dble(i + ir - 1) + 0.33d0*dble(j + ic - 1)**2)
            end do
        end do
        do j = 1, desc(nlac_)
            do i = 1, desc(nlar_)
                b(i, j) = dsin(3.d0*dble(i + ir - 1)**2 - 0.13d0*dble(j + ic - 1))
            end do
        end do

        c(:, :) = 0.d0
        if (rank .eq. 0) write (6, *) ' Starting '

        time = cclock()
        do i = 1, iter

            call PDGEMM('N', 'N', n, n, n, 1.0d0, A, 1, 1, desch, B, 1, 1, &
                    & desch, 0.1d0, C, 1, 1, desch)

            !      if(rank.eq.0) write(6,*) ' Iteration = ',i

        end do

        time = cclock() - time

    end if

    nflop = 2*dble(n)**3*iter

    if (rank .eq. 0) then

        write (6, *) ' Wall time =', time
        write (6, *) ' Speed in Gflops =', nflop/time/1d9
        write (6, *) ' Speed in Gflops/proc  =', nflop/time/1d9/dble(nproc)

    end if

#ifdef PARALLEL
    call mpi_finalize(ierr)
#endif

    stop

end
