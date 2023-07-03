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

module scal_lins
    implicit none
    integer nproc_ortho

#ifdef PARALLEL

contains

    !Calculate the dimensions of the submatrices according to the convention used in
    !distribute_column
    subroutine calc_dime(m, siz, last, num_procs)
        implicit none
        integer :: m, siz, last, num_procs

        siz = m/num_procs
        last = siz
        if (mod(m, num_procs) .ne. 0) then
            siz = siz + 1
            last = mod(m, siz)
        end if
    end subroutine calc_dime

!  subroutine

    !This subroutine is the modified version of Cavazzoni's setup_para of the DFT that
    !is inside parallel_mod.f90
    subroutine set_env(commu)
        use allio, only: me_blacs, np_blacs, world_cntx, np_ortho1, np_ortho, ortho_cntx, &
                         leg_ortho, ortho_comm, me_ortho, me_ortho1, ortho_comm_id, iflagerr, max_ortho
        implicit none
        include 'mpif.h'

        integer :: commu, ierr, nproc, rank
        integer, allocatable :: blacsmap(:, :)
        integer, allocatable :: blacsmap_other(:, :)
        integer, allocatable :: blacstmp(:, :)
        integer :: nprow, npcol, myrow, mycol, color, key
        integer :: row_id, col_id, row_comm, col_comm, mcol, irow, &
                   jcol, nk_opt, mcol_rep

#ifdef __SCALAPACK
#ifdef PARALLEL

        iflagerr = 0

        call MPI_COMM_RANK(commu, rank, ierr)
        call MPI_COMM_SIZE(commu, nproc, ierr)
        !
        ! Default values of communicators/id
        ! row = comm for processors dealing with a single k-point
        ! column = comm among k-points

        !
        ! setup SCALAPACK environment
        !
        ! in the case of k-points sampling the SCALAPACK grid is
        ! initialized within each processors pool.
        !
        ! SCALAPACK version needs the MPI environment.

        call BLACS_PINFO(me_blacs, np_blacs)
        ! me_blacs = BLACS process identifier
        call BLACS_GET(-1, 0, world_cntx)
        ! This subroutine factorizes the number of processors (NPROC)
        ! into NPROW and NPCOL according to the shape
        ! 'S' = square grid
        ! np_ortho(1) = # of row of the grid
        ! np_ortho(2) = # of column of the grid
        if (nproc .le. max_ortho .or. max_ortho .le. 0) then
            call grid2d_dims('S', nproc, np_ortho(1), np_ortho(2))
        else
            call grid2d_dims('S', max_ortho, np_ortho(1), np_ortho(2))
        end if

!   propagate nelorb in any case common to all processors
!    if(nelorb/np_ortho(1).lt.min_block) then
!      np_ortho(1)=nelorb/min_block
!      if(np_ortho(1)*min_block.ne.nelorb) np_ortho(1)=np_ortho(1)+1
!      np_ortho(2)=np_ortho(1)
!!    if(rank.eq.0) write(6,*) &
!     ' Warning using less number of processor to distribute the matrix , block too small !!! ',np_ortho(1),nelorb
!    endif

        np_ortho1 = np_ortho(1)*np_ortho(2)
!   nproc_ortho is defined also for processors not belonging to the pool
!   NB np_ortho1 is modified later by mpi_comm_size.
        nproc_ortho = np_ortho1

        !  here we choose the first processors omp threads can be set anyway
        !
        color = 0
        if (rank < np_ortho1) color = 1
        !
        leg_ortho = 1
        !
        !
        key = rank
        !
        !  initialize the communicator for the new group
        !  color = 1 : the group of processors performs matrix distribution
        !  color = 0 : processors not involved in the distribution
        !

        call MPI_COMM_SPLIT(commu, color, key, ortho_comm, ierr)

        if (ierr /= 0) &
            call errore(" init_ortho_group ", " error splitting communicator ", ierr)
        !
        !  Computes coordinates of the processors, in row maior order

        !
        call mpi_comm_size(ortho_comm, np_ortho1, ierr)
        call mpi_comm_rank(ortho_comm, me_ortho1, ierr)

        ! np_ortho1 = size of the group associated with ortho_comm
        ! me_ortho1 = rank of the current MPI task
        if (color == 1 .and. np_ortho1 /= np_ortho(1)*np_ortho(2)) &
            call errore(" init_ortho_group ", " wrong number of proc in ortho_comm ", ierr)
        !
        if (rank == 0 .and. me_ortho1 /= 0) &
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
        else ! processes not involved in distribution
            nprow = np_ortho1
            npcol = 1
            allocate (blacsmap(np_ortho1, 1))
            allocate (blacstmp(np_ortho1, 1))
            blacsmap = 0
            blacsmap(me_ortho1 + 1, 1) = me_blacs
        end if

        call mpi_allreduce(blacsmap, blacstmp, size(blacsmap), MPI_INTEGER, MPI_SUM, ortho_comm, ierr)

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

        ! final check for error flags
        call checkiflagerr(iflagerr, rank, "ERROR in initialization parallel environment")

#endif
#endif
        return

    end subroutine set_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !It takes the block distributed matrix s (m,m), of block size
    !in_blocks and distribute it on the mesh of out_procs processors,
    !then it solves the linear system associated with the matrix s and
    !the known term solu (input) that is known to all processors, in
    !output the master process has the solution stored in solu.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine para_syst(m, lds, in_blocks, rank, out_procs, comm, commtot, ictxt, solu, s)
        use descriptors
        use kpoints_mod, only: kaverage
        implicit none
        include 'mpif.h'
        integer :: ictxt, myrow, mycol, comm, commtot, first_proc
        integer :: ierr, rank
        !Variabili ausiliarie
        integer :: i, j
        !side is the side of the procs matrix
        integer, dimension(2) :: total, me
        integer :: m, lds, in_blocks, out_procs, nc, side, includeme
        !Dimension of the submatrices
        integer :: sizef, lastf, sizesq, lastsq
        !s is the input matrix block distributed by column, solu is the vector
        !of the known term in input and the solution in output
        real(8), dimension(m) :: solu, sd
        real(8), dimension(lds, in_blocks) :: s
        real(8), allocatable :: a(:, :), b(:, :), buff(:, :), vec(:)
        !descriptor for the matrix and the solution
        integer, allocatable :: req(:)
        integer, dimension(9) :: desca, descb
        integer, dimension(16) :: cavazz
!   real(8):: cost
!   call MPI_COMM_RANK (comm, rank, ierr)
!   call MPI_COMM_SIZE (comm, num_procs, ierr)

!   write(6,*) ' comm inside =',rank,comm,out_procs
!   if(rank.lt.out_procs) then
#ifdef __SCALAPACK
#ifdef PARALLEL

        first_proc = 0
        includeme = 1
        side = nint(sqrt(dble(out_procs)))

        !Changing the distribution from columns block to rows cyclic
        call calc_dime(m, sizef, lastf, out_procs)
        allocate (buff(m, sizef), b(sizef, m))
        b = 0.d0
        buff = 0.d0
        call distribute_column(m, lds, sizef, in_blocks, rank, in_blocks, out_procs, s, buff, comm)
        call col2row(m, rank, out_procs, sizef, lastf, b, buff)

        !Blacs initialization for all processors
        call blacs_gridinfo(ictxt, side, side, myrow, mycol)

        if (rank .lt. out_procs) then

            !Distributing the matrix with Cavazzoni algorithm
            call calc_dime(m, sizesq, lastsq, side)
            total(:) = side
            me(1) = myrow
            me(2) = mycol
            call descla_init(cavazz, m, m, total, me, comm, includeme)
            allocate (a(sizesq, sizesq))
            a = 0.d0
            call cyc2blk_redist(m, b, sizef, a, sizesq, cavazz)

            !Distributing the vector of the known term (lapack name b, here solu and vec)
            allocate (vec(sizesq))
            vec = 0.d0
            if (mod(rank, side) .eq. 0) then
                j = sizesq
                i = rank/side
                if (i .eq. side - 1) j = lastsq
                vec(1:j) = solu(i*sizesq + 1:i*sizesq + j)
!      write(6,*) ' vec stored ',rank,j,i*sizesq+1,i*sizesq+j,sum(vec(1:j))
            end if

            !Descriptor initialization
            call descinit(desca, m, m, cavazz(nlax_), cavazz(nlax_), &
                          0, 0, ictxt, sizesq, ierr)
            call descinit(descb, m, 1, cavazz(nlax_), 1, &
                          0, 0, ictxt, sizesq, ierr)
!   cost=sum(a(:,:))
!   else
!   cost=0.d0
        end if

!   call reduce_base_real(1,cost,comm,-1)
!   if(rank.eq.0) write(6,*) ' sum matrix elements after reshuff =',cost

        if (rank .lt. out_procs) then

            !Cholesky
            call pdpotrf("U", m, a, 1, 1, desca, ierr)

            !Linear system
            call pdpotrs("U", m, 1, a, 1, 1, desca, vec, 1, 1, &
                         descb, ierr)

            !Rebuilding the solu vector
!!TEST
!    if ((rank-side).lt.0) then
!       i=mod(rank,side)
!       j=sizesq
!       if (i.eq.side-1) j=lastsq
!       solu(1+i*sizesq:j+i*sizesq)=vec(1:j)
!    end if
!    do i=0,side-1
!       j=sizesq
!       if (i.eq.side-1) j=lastsq
!       call mpi_bcast (solu(1+i*sizesq), j, MPI_DOUBLE, i, comm, ierr )
!    end do
!    write (6,*) rank, vec
            solu = 0.d0
            if (mod(rank, side) .eq. 0) then
                j = sizesq
                i = rank/side
                if (i .eq. side - 1) j = lastsq
                solu(i*sizesq + 1:i*sizesq + j) = vec(1:j)
!      write(6,*) ' vec stored ',rank,j,i*sizesq+1,i*sizesq+j,sum(vec(1:j))
            end if
            deallocate (vec, a)
        else
            solu = 0.d0
        end if
        call reduce_base_real(m, solu, commtot, -1)
        deallocate (b, buff)
#endif
#endif
    end subroutine para_syst

    !This subroutine reshuffle a matrix divided between the processors by column.
    !The initial configuration is over the first processes in block of in_blocks elements
    !and is stored in the elements S of dimension (m,sizei), the final configuration
    !is CYCLIC over num_procs processors and is stored in b(m,sizef). Last is the
    !dimension of the last block of columns of the input matrix.
    subroutine distribute_column(m, lds, dimm, dims, rank, in_blocks, out_procs, s, b, comm)
        implicit none
        include 'mpif.h'
        integer :: rank, dimm, dims, lds, out_procs, in_blocks, ierr, comm, in_procs
        !npone (n+1) is the number of row with one row more than m/out_procs
        integer :: m, sizei, sizef, last, npone
        real(8) :: s(lds, dims), b(m, dimm)
        !Adress to send/recv, pos is the original position in the mxm matrix
        integer :: adress, pos, last_el
        integer :: i, j, count
        integer, allocatable :: req(:), status(:, :)

#ifdef __SCALAPACK
#ifdef PARALLEL
        sizei = in_blocks
        last = mod(m, in_blocks)

        in_procs = m/sizei
        if (last .ne. 0) then
            in_procs = in_procs + 1
        else
            last = sizei
        end if
!    call calc_dime(m, sizei, last, in_procs)

        sizef = m/out_procs
        npone = mod(m, out_procs)
        last_el = sizef
        if (rank .lt. npone) last_el = last_el + 1
        count = 0
        allocate (req(sizei + sizef + 1), status(MPI_STATUS_SIZE, sizei + sizef + 1))
        !Receaving
        if (rank .lt. out_procs) then
            do i = 1, last_el
                pos = (i - 1)*out_procs + rank + 1
                adress = (pos - 1)/sizei
                count = count + 1
                call mpi_irecv(b(1, i), m, MPI_DOUBLE, adress, pos, comm, req(count), ierr)
            end do
        end if

        !Sending
        if (rank .lt. in_procs - 1) last = sizei
        if (rank .lt. in_procs) then
            do j = 1, last
                pos = j + sizei*rank
                adress = mod(pos - 1, out_procs)
                count = count + 1
                call mpi_isend(s(1, j), m, MPI_DOUBLE, adress, pos, comm, req(count), ierr)
            end do
        end if

        call MPI_WAITALL(count, req, status, ierr)
        deallocate (status, req)
        call mpi_barrier(comm, ierr)
#endif
#endif
    end subroutine distribute_column

    !Trasforming the coloumn in rows
    subroutine col2row(m, rank, procs, siz, last, b, buff)
        implicit none
        !nc is the number of column of the buffer
        integer :: m, rank, procs, siz, last, nc, i, j
        real(8), dimension(siz, m) :: b
        real(8), dimension(m, siz) :: buff

        if (rank .lt. procs) then
            nc = m/procs
            if (rank .lt. mod(m, procs)) nc = nc + 1
            do i = 1, m
                do j = 1, nc
                    b(j, i) = buff(i, j)
                end do
            end do
        end if

    end subroutine col2row

#endif

end module scal_lins

