!TL off
! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2002 FPMD group
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! Copyright (C) 2002-2008 Quantum ESPRESSO group
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

!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!     ==================================================================
subroutine DSORTX(COUNT, INUTILE, N, INDEX)
    !     ==--------------------------------------------------------------==
    !     == Sorting routine for the reciprocal space vectors (g)         ==
    !     == KB07AD HANDLES DOUBLE PRECISION VARIABLES                    ==
    !     == STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)            ==
    !     == THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2)   ==
    !     == NUMBERS TO BE SORTED. THIS IS MORE THAN THE IBM VIRTUAL      ==
    !     == MEMORY SPACE WILL HOLD .                                     ==
    !     ==--------------------------------------------------------------==
    implicit none
    integer N, I, M, IS, IA, LA, if, MLOOP, IFKA, IS1, J, INT, IY, INTEST, K, IFK, K1, IP, LNGTH
    integer index(N)
    real*8 count(N), AV, X
    integer INUTILE, MARK
    dimension MARK(50)
    !     ==--------------------------------------------------------------==
    !     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
    !     ==--------------------------------------------------------------==
    do I = 1, N
        index(I) = I
    end do
    !     ==--------------------------------------------------------------==
    !     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
    !     ==--------------------------------------------------------------==
    if (N .eq. 1) goto 200
    if (N .ge. 1) goto 30
    goto 200
    !     ==--------------------------------------------------------------==
    !     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
    !     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
    !     ==--------------------------------------------------------------==
30  M = 12
    !     ==--------------------------------------------------------------==
    !     == SET UP INITIAL VALUES.                                       ==
    !     ==--------------------------------------------------------------==
    LA = 2
    IS = 1
    if = N
    do 190 MLOOP = 1, N
        !     ==--------------------------------------------------------------==
        !     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
        !     ==--------------------------------------------------------------==
        IFKA = if - IS
        if ((IFKA + 1) .gt. M) goto 70
        !     ==--------------------------------------------------------------==
        !     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
        !     ==--------------------------------------------------------------==
        IS1 = IS + 1
        do 60 J = IS1, if
            I = J
40          if (count(I - 1) .lt. count(I)) goto 60
            if (count(I - 1) .gt. count(I)) goto 50
            if (index(I - 1) .lt. index(I)) goto 60
50          AV = count(I - 1)
            count(I - 1) = count(I)
            count(I) = AV
            INT = index(I - 1)
            index(I - 1) = index(I)
            index(I) = INT
            I = I - 1
            if (I .gt. IS) goto 40
60          continue
            LA = LA - 2
            goto 170
            !     ==--------------------------------------------------------------==
            !     ==                      *  QUICKSORT        **                  ==
            !     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
            !     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S==
            !     == HIGHEST ADDRESS.                                             ==
            !     ==--------------------------------------------------------------==
70          IY = (IS + if)/2
            X = count(IY)
            INTEST = index(IY)
            count(IY) = count(if)
            index(IY) = index(if)
            !     ==--------------------------------------------------------------==
            !     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
            !     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
            !     == OF X .                                                       ==
            !     ==--------------------------------------------------------------==
            K = 1
            IFK = if
            !     ==--------------------------------------------------------------==
            !     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
            !     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
            !     == NECESSARY, UNTIL THEY MEET .                                 ==
            !     ==--------------------------------------------------------------==
            do 110 I = IS, if
                if (X .gt. count(I)) goto 110
                if (X .lt. count(I)) goto 80
                if (INTEST .gt. index(I)) goto 110
80              if (I .ge. IFK) goto 120
                count(IFK) = count(I)
                index(IFK) = index(I)
                K1 = K
                do 100 K = K1, IFKA
                    IFK = if - K
                    if (count(IFK) .gt. X) goto 100
                    if (count(IFK) .lt. X) goto 90
                    if (INTEST .le. index(IFK)) goto 100
90                  if (I .ge. IFK) goto 130
                    count(I) = count(IFK)
                    index(I) = index(IFK)
                    GO TO 110
100                 continue
                    goto 120
110                 continue
                    !     ==--------------------------------------------------------------==
                    !     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
                    !     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
                    !     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
                    !     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
                    !     == INDEPENDENTLY .                                              ==
                    !     ==--------------------------------------------------------------==
120                 count(IFK) = X
                    index(IFK) = INTEST
                    IP = IFK
                    goto 140
130                 count(I) = X
                    index(I) = INTEST
                    IP = I
                    !     ==--------------------------------------------------------------==
                    !     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
                    !     ==--------------------------------------------------------------==
140                 if ((IP - IS) .gt. (if - IP)) goto 150
                    MARK(LA) = if
                    MARK(LA - 1) = IP + 1
                    if = IP - 1
                    goto 160
150                 MARK(LA) = IP - 1
                    MARK(LA - 1) = IS
                    IS = IP + 1
                    !     ==--------------------------------------------------------------==
                    !     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
                    !     ==--------------------------------------------------------------==
160                 LNGTH = if - IS
                    if (LNGTH .le. 0) goto 180
                    !     ==--------------------------------------------------------------==
                    !     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
                    !     ==--------------------------------------------------------------==
                    LA = LA + 2
                    goto 190
170                 if (LA .le. 0) goto 200
                    !     ==--------------------------------------------------------------==
                    !     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
                    !     ==--------------------------------------------------------------==
180                 if = MARK(LA)
                    IS = MARK(LA - 1)
190             end do
                !     ==--------------------------------------------------------------==
200             return
            end

            !
            ! Copyright (C) 2002 FPMD group
            ! This file is distributed under the terms of the
            ! GNU General Public License. See the file `License'
            ! in the root directory of the present distribution,
            ! or http://www.gnu.org/copyleft/gpl.txt .
            !

            subroutine GRID2D_DIMS(grid_shape, nproc, nprow, npcol)
                !
                ! This subroutine factorizes the number of processors (NPROC)
                ! into NPROW and NPCOL according to the shape
                !
                !    Written by Carlo Cavazzoni
                !
                implicit none
                character, intent(IN) :: grid_shape
                integer, intent(IN) :: nproc
                integer, intent(OUT) :: nprow, npcol
                integer :: sqrtnp, i
                !
                sqrtnp = int(sqrt(real(nproc) + 0.1))
                !
                if (grid_shape == 'S') then
                    ! Square grid
                    nprow = sqrtnp
                    npcol = sqrtnp
                else
                    ! Rectangular grid
                    do i = 1, sqrtnp + 1
                        if (mod(nproc, i) == 0) nprow = i
                    end do
                    npcol = nproc/nprow
                end if
                return
            end subroutine

            subroutine GRID2D_COORDS(order, rank, nprow, npcol, row, col)
                !
                !  this subroutine compute the cartesian coordinetes "row" and "col"
                !  of the processor whose MPI task id is "rank".
                !  Note that if the rank is larger that the grid size
                !  all processors whose MPI task id is greather or equal
                !  than nprow * npcol are placed on the diagonal extension of the grid itself
                !
                implicit none
                character, intent(IN) :: order
                integer, intent(IN) :: rank ! process index starting from 0
                integer, intent(IN) :: nprow, npcol ! dimensions of the processor grid
                integer, intent(OUT) :: row, col
                if (rank >= 0 .and. rank < nprow*npcol) then
                    if (order == 'C' .or. order == 'c') then
                        !  grid in COLUMN MAJOR ORDER
                        row = mod(rank, nprow)
                        col = rank/nprow
                    else
                        !  grid in ROW MAJOR ORDER
                        row = rank/npcol
                        col = mod(rank, npcol)
                    end if
                else
                    row = rank
                    col = rank
                end if
                return
            end subroutine

            subroutine GRID2D_RANK(order, nprow, npcol, row, col, rank)
                !
                !  this subroutine compute the processor MPI task id "rank" of the processor
                !  whose cartesian coordinate are "row" and "col".
                !  Note that the subroutine assume cyclic indexing ( row = nprow = 0 )
                !
                implicit none
                character, intent(IN) :: order
                integer, intent(OUT) :: rank ! process index starting from 0
                integer, intent(IN) :: nprow, npcol ! dimensions of the processor grid
                integer, intent(IN) :: row, col

                if (order == 'C' .or. order == 'c') then
                    !  grid in COLUMN MAJOR ORDER
                    rank = mod(row + nprow, nprow) + mod(col + npcol, npcol)*nprow
                else
                    !  grid in ROW MAJOR ORDER
                    rank = mod(col + npcol, npcol) + mod(row + nprow, nprow)*npcol
                end if
                !
                return
            end subroutine

            !
            ! Copyright (C) 2001-2013 Quantum ESPRESSO group
            ! This file is distributed under the terms of the
            ! GNU General Public License. See the file `License'
            ! in the root directory of the present distribution,
            ! or http://www.gnu.org/copyleft/gpl.txt .
            !-----------------------------------------------------------------------

            integer function ldim_block_sca(gdim, np, me)

                !   gdim = global dimension of distributed array
                !   np   = number of processors
                !   me   = index of the calling processor (starting from 0)
                !
                !   this function return the number of elements of the distributed array
                !   stored in the local memory of the processor "me" for equal block
                !   data distribution, all block have the same size but the last one.
                !   Example of the block distribution of 10 elements array a on 4 processors
                !   array elements  |  PEs
                !    a(1)           |   0
                !    a(2)           |   0
                !    a(3)           |   0
                !    a(4)           |   1
                !    a(5)           |   1
                !    a(6)           |   1
                !    a(7)           |   2
                !    a(8)           |   2
                !    a(9)           |   2
                !    a(10)          |   3

                implicit none
                integer :: gdim, np, me, nb

                if (me >= np .or. me < 0) then
                    write (6, *) ' ** ldim_block: arg no. 3 out of range '
                    stop
                end if

                nb = int(gdim/np)
                if (mod(gdim, np) /= 0) then
                    nb = nb + 1
                    ! ... last processor take the rest
                    if (me == (np - 1)) nb = gdim - (np - 1)*nb
                end if

                ldim_block_sca = nb

                return
            end function ldim_block_sca

            integer function lind_block_sca(ig, nx, np, me)
                !
                !   INPUT :
                !      ig  global index of the x dimension of array element
                !      nx  dimension of the global array
                !      np  number of processor in the x dimension of the processors grid
                !      me  index of the local processor in the processor grid
                !                (starting from zero)
                !
                !   OUTPUT :
                !
                !      lind_block_sca return the local index corresponding to the
                !      global index "ig" for an equal block distribution
                !

                implicit none

                integer :: ig, nx, np, me, nb

                nb = int(nx/np)
                if (mod(nx, np) /= 0) nb = nb + 1

                lind_block_sca = ig - me*nb

                return
            end function lind_block_sca

            integer function gind_block_sca(lind, n, np, me)

                !  This function computes the global index of a distributed array entry
                !  pointed to by the local index lind of the process indicated by me.
                !  lind      local index of the distributed matrix entry.
                !  N         is the size of the global array.
                !  me        The coordinate of the process whose local array row or
                !            column is to be determined.
                !  np        The total number processes over which the distributed
                !            matrix is distributed.

                integer, intent(IN) :: lind, n, me, np
                integer nb

                if (me >= np .or. me < 0) then
                    write (6, *) ' ** ldim_block: arg no. 3 out of range '
                    stop
                end if

                nb = int(n/np)
                if (mod(n, np) /= 0) nb = nb + 1

                gind_block_sca = lind + me*nb

                return

            end function gind_block_sca

            integer function ldim_cyclic(gdim, np, me)

                !   gdim = global dimension of distributed array
                !   np   = number of processors
                !   me   = index of the calling processor (starting from 0)
                !
                !   this function return the number of elements of the distributed array
                !   stored in the local memory of the processor "me" for a cyclic
                !   data distribution.
                !   Example of the cyclic distribution of a 10 elements array on 4 processors
                !   array elements  |  PEs
                !    a(1)           |   0
                !    a(2)           |   1
                !    a(3)           |   2
                !    a(4)           |   3
                !    a(5)           |   0
                !    a(6)           |   1
                !    a(7)           |   2
                !    a(8)           |   3
                !    a(9)           |   0
                !    a(10)          |   1

                implicit none
                integer :: gdim, np, me, r, q

                if (me >= np .or. me < 0) then
                    write (6, *) ' ** ldim_cyclic: arg no. 3 out of range '
                    stop
                end if

                q = int(gdim/np)
                r = mod(gdim, np)

                if (me .lt. r) then
                    ldim_cyclic = q + 1
                else
                    ldim_cyclic = q
                end if

                return
            end function ldim_cyclic

            !
            ! Copyright (C) 2002-2008 Quantum ESPRESSO group
            ! This file is distributed under the terms of the
            ! GNU General Public License. See the file `License'
            ! in the root directory of the present distribution,
            ! or http://www.gnu.org/copyleft/gpl.txt .

            !----------------------------------------------------------------------------
            subroutine reduce_base_real(dim, ps, comm, root)
                use buffer, only: bufdim, buff_reduce
                !----------------------------------------------------------------------------
                !
                ! ... sums a distributed variable ps(dim) over the processors.
                ! ... This version uses a fixed-length buffer of appropriate (?) dim
                ! ...              uses SHMEM if available, MPI otherwhise
                !
                implicit none
                !
                integer, intent(IN) :: dim
                real*8, intent(INOUT) :: ps(dim)
                integer, intent(IN) :: comm ! communicator
                integer, intent(IN) :: root ! if root <  0 perform a reduction to all procs
                ! if root >= 0 perform a reduce only to root proc.
                !
#if defined (PARALLEL)

                include 'mpif.h'
                !
                integer :: info, n, nbuf, nproc, myid, maxb
                !

                maxb = min(dim, bufdim)

                ! allocate(buff_reduce(maxb))
                !
                call mpi_comm_size(comm, nproc, info)
                if (info /= 0) call errore('reduce_base_real', 'error in mpi_comm_size', info)

                call mpi_comm_rank(comm, myid, info)
                if (info /= 0) call errore('reduce_base_real', 'error in mpi_comm_rank', info)
                !
                if (dim <= 0 .or. nproc <= 1) return
                !
                ! ... synchronize processes
                !
#ifdef UNREL
                call mpi_barrier(comm, info)
!$omp barrier
                if (info /= 0) call errore('reduce_base_real', 'error in mpi_barrier', info)
#endif
                !
                if (dim .eq. maxb) then
                    buff_reduce(1:maxb) = ps(1:maxb)

                    if (root >= 0) then
                        call MPI_REDUCE(buff_reduce, ps, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info)
                        if (info /= 0) call errore('reduce_base_real', 'error in mpi_reduce 1', info)
                    else
                        call MPI_ALLREDUCE(buff_reduce, ps, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info)
                        if (info /= 0) call errore('reduce_base_real', 'error in mpi_allreduce 1', info)
                    end if
                    !

                else

                    nbuf = dim/maxb
                    !

                    do n = 1, nbuf
                        !
                        if (root >= 0) then
                            call MPI_REDUCE(ps(1 + (n - 1)*maxb), buff_reduce, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info)
                            if (info /= 0) call errore('reduce_base_real', 'error in mpi_reduce 1', info)
                        else
                            call MPI_ALLREDUCE(ps(1 + (n - 1)*maxb), buff_reduce, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info)
                            if (info /= 0) call errore('reduce_base_real', 'error in mpi_allreduce 1', info)
                        end if
                        !
                        if (root < 0) then
                            ps((1 + (n - 1)*maxb):(n*maxb)) = buff_reduce(1:maxb)
                        else if (root == myid) then
                            ps((1 + (n - 1)*maxb):(n*maxb)) = buff_reduce(1:maxb)
                        end if
                        !
                    end do
                    !
                    ! ... possible remaining elements < maxb
                    !
                    if ((dim - nbuf*maxb) > 0) then
                        !
                        if (root >= 0) then
                            call MPI_REDUCE(ps(1 + nbuf*maxb), buff_reduce, (dim - nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info)
                            if (info /= 0) call errore('reduce_base_real', 'error in mpi_reduce 2', info)
                        else
                            call MPI_ALLREDUCE(ps(1 + nbuf*maxb), buff_reduce, (dim - nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info)
                            if (info /= 0) call errore('reduce_base_real', 'error in mpi_allreduce 2', info)
                        end if
                        !
                        if (root < 0) then
                            ps((1 + nbuf*maxb):dim) = buff_reduce(1:(dim - nbuf*maxb))
                        else if (root == myid) then
                            ps((1 + nbuf*maxb):dim) = buff_reduce(1:(dim - nbuf*maxb))
                        end if
                        !
                    end if
                    !

                end if

#endif
                !
                return
                !
            end subroutine reduce_base_real

            !----------------------------------------------------------------------------
            subroutine reduce_base_complex(dim, ps, comm, root)
                use buffer, only: bufdim, buff_reducec
                !----------------------------------------------------------------------------
                !
                ! ... sums a distributed variable ps(dim) over the processors.
                ! ... This version uses a fixed-length buffer of appropriate (?) dim
                ! ...              uses SHMEM if available, MPI otherwhise
                !
                implicit none
                !
                integer, intent(IN) :: dim
                complex(8), intent(INOUT) :: ps(dim)
                integer, intent(IN) :: comm ! communicator
                integer, intent(IN) :: root ! if root <  0 perform a reduction to all procs
                ! if root >= 0 perform a reduce only to root proc.
                !
#ifdef PARALLEL

                include 'mpif.h'
                !
                integer :: info, n, nbuf, nproc, myid, maxb
                !
                !
                maxb = min(dim, bufdim)
                !
                call mpi_comm_size(comm, nproc, info)
                if (info /= 0) call errore('reduce_base_complex', 'error in mpi_comm_size', info)

                call mpi_comm_rank(comm, myid, info)
                if (info /= 0) call errore('reduce_base_complex', 'error in mpi_comm_rank', info)
                !
                if (dim <= 0 .or. nproc <= 1) return ! trivial output
                !
                ! ... synchronize processes
                !
#ifdef UNREL
                call mpi_barrier(comm, info)
!$omp barrier
                if (info /= 0) call errore('reduce_base_complex', 'error in mpi_barrier', info)
#endif
                !
                if (dim .eq. maxb) then

                    buff_reducec(1:maxb) = ps(1:maxb)
                    if (root >= 0) then
                        call MPI_REDUCE(buff_reducec, ps, maxb, MPI_COMPLEX16, MPI_SUM, root, comm, info)
                        if (info /= 0) call errore('reduce_base_complex', 'error in mpi_reduce 1', info)
                    else
                        call MPI_ALLREDUCE(buff_reducec, ps, maxb, MPI_COMPLEX16, MPI_SUM, comm, info)
                        if (info /= 0) call errore('reduce_base_complex', 'error in mpi_allreduce 1', info)
                    end if
                    !
                else
                    !
                    nbuf = dim/maxb
                    !
                    do n = 1, nbuf
                        !
                        if (root >= 0) then
                            call MPI_REDUCE(ps(1 + (n - 1)*maxb), buff_reducec, maxb, MPI_COMPLEX16, MPI_SUM, root, comm, info)
                            if (info /= 0) call errore('reduce_base_complex', 'error in mpi_reduce 1', info)
                        else
                            call MPI_ALLREDUCE(ps(1 + (n - 1)*maxb), buff_reducec, maxb, MPI_COMPLEX16, MPI_SUM, comm, info)
                            if (info /= 0) call errore('reduce_base_complex', 'error in mpi_allreduce 1', info)
                        end if
                        !
                        if (root < 0) then
                            ps((1 + (n - 1)*maxb):(n*maxb)) = buff_reducec(1:maxb)
                        else if (root == myid) then
                            ps((1 + (n - 1)*maxb):(n*maxb)) = buff_reducec(1:maxb)
                        end if
                        !
                    end do
                    !
                    ! ... possible remaining elements < maxb
                    !
                    if ((dim - nbuf*maxb) > 0) then
                        !
                        if (root >= 0) then
                            call MPI_REDUCE(ps(1 + nbuf*maxb), buff_reducec, (dim - nbuf*maxb), MPI_COMPLEX16, MPI_SUM, root, comm, info)
                            if (info /= 0) call errore('reduce_base_complex', 'error in mpi_reduce 2', info)
                        else
                            call MPI_ALLREDUCE(ps(1 + nbuf*maxb), buff_reducec, (dim - nbuf*maxb), MPI_COMPLEX16, MPI_SUM, comm, info)
                            if (info /= 0) call errore('reduce_base_complex', 'error in mpi_allreduce 2', info)
                        end if
                        !
                        if (root < 0) then
                            ps((1 + nbuf*maxb):dim) = buff_reducec(1:(dim - nbuf*maxb))
                        else if (root == myid) then
                            ps((1 + nbuf*maxb):dim) = buff_reducec(1:(dim - nbuf*maxb))
                        end if
                        !
                    end if
                    !
                end if

#endif
                !
                return
                !
            end subroutine reduce_base_complex

!
! Copyright (C) 2002-2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
!  Wrapper for MPI implementations that have problems with large messages
!
!  In some Cray XD1 systems the communicaction subsystem
!  crashes when message exceed a given size, so we need
!  to break down MPI communications in smaller pieces
!

!=----------------------------------------------------------------------------=!
!
            subroutine BCAST_INTEGER(array, n, root, gid)
                !  In most architectures the mpi_allreduce is much more efficient than a bcast.
                !  Thus I have replaced the bcast with an mpi_allreduce, also safer because much
                !  more used.
                use buffer, only: bufdim, buffi_reduce
                implicit none
                integer :: n, root, gid, ierr
                integer :: array(n)
#if defined PARALLEL
                include 'mpif.h'
                integer :: nblk, blksiz, msgsiz, iblk, istart, i, myid, nproc
                !       common /mp_base_integer/ buffi
                call mpi_comm_rank(gid, myid, ierr)
                call mpi_comm_size(gid, nproc, ierr)
#ifdef UNREL
                call mpi_barrier(gid, ierr)
!$omp barrier
                if (ierr /= 0) call errore('BCAST_INTEGER', 'error in mpi_barrier', ierr)
#endif

#ifdef __KCOMP
                if (nproc .ge. 8) then
#else
                    if (nproc .ge. 128) then
#endif
                        if (n <= bufdim) then
                            if (myid .eq. root) then
                                buffi_reduce(1:n) = array(1:n)
                            else
                                buffi_reduce(1:n) = 0
                            end if
                            !          CALL MPI_BCAST( buffi_reduce, n, MPI_INTEGER, root, gid, ierr )
                            !          array(1:n)=buffi_reduce(1:n)
                            call MPI_ALLREDUCE(buffi_reduce, array, n, MPI_INTEGER, MPI_SUM, gid, ierr)
                            if (ierr /= 0) call errore(' bcast_integer ', ' error in mpi_allreduce 1 ', ierr)
                        else
                            nblk = n/bufdim
                            blksiz = bufdim
                            do iblk = 1, nblk
                                istart = (iblk - 1)*bufdim + 1
                                if (myid .eq. root) then
                                    buffi_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                else
                                    buffi_reduce(1:blksiz) = 0
                                end if
                                !           CALL MPI_BCAST( buffi_reduce, blksiz, MPI_INTEGER, root, gid, ierr )
                                !           array(istart:istart+blksiz-1)=buffi_reduce(1:blksiz)
                                call MPI_ALLREDUCE(buffi_reduce, array(istart), blksiz, MPI_INTEGER, MPI_SUM, gid, ierr)
                                if (ierr /= 0) call errore(' bcast_integer ', ' error in mpi_allreduce 2 ', ierr)
                            end do
                            blksiz = mod(n, bufdim)
                            if (blksiz > 0) then
                                istart = nblk*bufdim + 1
                                if (myid .eq. root) then
                                    buffi_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                else
                                    buffi_reduce(1:blksiz) = 0
                                end if
                                !       CALL MPI_BCAST( buffi_reduce, blksiz, MPI_INTEGER, root, gid, ierr )
                                !             array(istart:istart+blksiz-1)=buffi_reduce(1:blksiz)
                                call MPI_ALLREDUCE(buffi_reduce, array(istart), blksiz, MPI_INTEGER, MPI_SUM, gid, ierr)
                                if (ierr /= 0) call errore(' bcast_integer ', ' error in mpi_allreduce 3 ', ierr)
                            end if
                        end if

                    else

                        if (n <= bufdim) then
                            buffi_reduce(1:n) = array(1:n)
                            call MPI_BCAST(buffi_reduce, n, MPI_INTEGER, root, gid, ierr)
                            array(1:n) = buffi_reduce(1:n)
                            if (ierr /= 0) call errore(' bcast_integer ', ' error in mpi_bcast 1 ', ierr)
                        else
                            nblk = n/bufdim
                            blksiz = bufdim
                            do iblk = 1, nblk
                                istart = (iblk - 1)*bufdim + 1
                                buffi_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                call MPI_BCAST(buffi_reduce, blksiz, MPI_INTEGER, root, gid, ierr)
                                array(istart:istart + blksiz - 1) = buffi_reduce(1:blksiz)
                                if (ierr /= 0) call errore(' bcast_integer ', ' error in mpi_bcast 2 ', ierr)
                            end do
                            blksiz = mod(n, bufdim)
                            if (blksiz > 0) then
                                istart = nblk*bufdim + 1
                                buffi_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                call MPI_BCAST(buffi_reduce, blksiz, MPI_INTEGER, root, gid, ierr)
                                array(istart:istart + blksiz - 1) = buffi_reduce(1:blksiz)
                                if (ierr /= 0) call errore(' bcast_integer ', ' error in mpi_bcast 3 ', ierr)
                            end if
                        end if
                    end if
#endif
                    return
                    end subroutine BCAST_INTEGER

                    subroutine BCAST_REAL(array, n, root, gid)
                        !  In most architectures the mpi_allreduce is much more efficient than a bcast.
                        !  Thus I have replaced the bcast with an mpi_allreduce, also safer because much
                        !  more used.
                        use buffer, only: bufdim, buff_reduce
                        implicit none
                        integer :: n, root, gid, ierr
                        real*8 :: array(n)
#if defined PARALLEL
                        include 'mpif.h'
                        integer :: nblk, blksiz, msgsiz, iblk, istart, i, myid, nproc
                        !       common /mp_base_real/ buff
                        call mpi_comm_size(gid, nproc, ierr)
                        call mpi_comm_rank(gid, myid, ierr)
#ifdef UNREL
                        call mpi_barrier(gid, ierr)
!$omp barrier
                        if (ierr /= 0) call errore('BCAST_REAL', 'error in mpi_barrier', ierr)
#endif

#if __KCOMP
                        if (nproc .ge. 8) then
#else
                            if (nproc .ge. 128) then
#endif
                                if (n <= bufdim) then
                                    if (myid .eq. root) then
                                        buff_reduce(1:n) = array(1:n)
                                    else
                                        buff_reduce(1:n) = 0.d0
                                    end if
                                    !        CALL MPI_BCAST( buff_reduce, n, MPI_DOUBLE_PRECISION, root, gid, ierr )
                                    !        array(1:n)=buff_reduce(1:n)
                                    call MPI_ALLREDUCE(buff_reduce, array, n, MPI_DOUBLE_PRECISION, MPI_SUM, gid, ierr)
                                    if (ierr /= 0) call errore('bcast_real', 'error in mpi_allreduce 1', ierr)
                                else
                                    nblk = n/bufdim
                                    blksiz = bufdim
                                    do iblk = 1, nblk
                                        istart = (iblk - 1)*bufdim + 1
                                        if (myid .eq. root) then
                                            buff_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                        else
                                            buff_reduce(1:blksiz) = 0.d0
                                        end if
                                        call MPI_ALLREDUCE(buff_reduce, array(istart), blksiz, MPI_DOUBLE_PRECISION, MPI_SUM, gid, ierr)
                                        !   CALL MPI_BCAST( buff_reduce, blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
                                        !             array(istart:istart+blksiz-1)=buff_reduce(1:blksiz)
                                        if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_allreduce 2 ', ierr)
                                    end do
                                    blksiz = mod(n, bufdim)
                                    if (blksiz > 0) then
                                        istart = nblk*bufdim + 1
                                        if (myid .eq. root) then
                                            buff_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                        else
                                            buff_reduce(1:blksiz) = 0.d0
                                        end if
                                        call MPI_ALLREDUCE(buff_reduce, array(istart), blksiz, MPI_DOUBLE_PRECISION, MPI_SUM, gid, ierr)
                                        !   CALL MPI_BCAST( buff_reduce, blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
                                        !             array(istart:istart+blksiz-1)=buff_reduce(1:blksiz)
                                        if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_allreduce 3 ', ierr)
                                    end if
                                end if
                            else

                                if (n <= bufdim) then

                                    buff_reduce(1:n) = array(1:n)
                                    call MPI_BCAST(buff_reduce, n, MPI_DOUBLE_PRECISION, root, gid, ierr)
                                    array(1:n) = buff_reduce(1:n)
                                    if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_bcast 1 ', ierr)
                                else
                                    nblk = n/bufdim
                                    blksiz = bufdim
                                    do iblk = 1, nblk
                                        istart = (iblk - 1)*bufdim + 1
                                        buff_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                        call MPI_BCAST(buff_reduce, blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr)
                                        array(istart:istart + blksiz - 1) = buff_reduce(1:blksiz)

                                        if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_bcast 2 ', ierr)
                                    end do
                                    blksiz = mod(n, bufdim)
                                    if (blksiz > 0) then
                                        istart = nblk*bufdim + 1
                                        buff_reduce(1:blksiz) = array(istart:istart + blksiz - 1)
                                        call MPI_BCAST(buff_reduce, blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr)
                                        array(istart:istart + blksiz - 1) = buff_reduce(1:blksiz)

                                        if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_bcast 3 ', ierr)
                                    end if
                                end if

                            end if
#endif
                            return
                            end subroutine BCAST_REAL

                            subroutine BCAST_COMPLEX(array, n, root, gid)
                                !  In most architectures the mpi_allreduce is much more efficient than a bcast.
                                !  Thus I have replaced the bcast with an mpi_allreduce, also safer because much
                                !  more used.
                                use buffer, only: bufdim, buff_reducec

                                implicit none
                                integer :: n, root, gid, ierr
                                complex(8) :: array(n)

#ifdef PARALLEL

                                include 'mpif.h'
                                integer :: nblk, blksiz, msgsiz, iblk, istart, i, myid, nproc
!       common /mp_base_real/ buff
                                call mpi_comm_size(gid, nproc, ierr)
                                call mpi_comm_rank(gid, myid, ierr)

#ifdef UNREL
                                call mpi_barrier(gid, ierr)
!$omp barrier
                                if (ierr /= 0) call errore('BCAST_COMPLEX', 'error in mpi_barrier', ierr)
#endif

                                if (nproc .ge. 1024) then
                                    if (n <= bufdim) then
                                        if (myid .eq. root) then
                                            buff_reducec(1:n) = array(1:n)
                                        else
                                            buff_reducec(1:n) = (0.d0, 0.d0)
                                        end if
                                        call MPI_ALLREDUCE(buff_reducec, array, n, MPI_COMPLEX16, MPI_SUM, gid, ierr)
                                        if (ierr /= 0) call errore('bcast_complex', 'error in mpi_allreduce 1', ierr)
                                    else
                                        nblk = n/bufdim
                                        blksiz = bufdim
                                        do iblk = 1, nblk
                                            istart = (iblk - 1)*bufdim + 1
                                            if (myid .eq. root) then
                                                buff_reducec(1:blksiz) = array(istart:istart + blksiz - 1)
                                            else
                                                buff_reducec(1:blksiz) = (0.d0, 0.d0)
                                            end if
                                            call MPI_ALLREDUCE(buff_reducec, array(istart), blksiz, MPI_COMPLEX16, MPI_SUM, gid, ierr)
!   CALL MPI_BCAST( buff_reduce, blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
!             array(istart:istart+blksiz-1)=buff_reduce(1:blksiz)
                                            if (ierr /= 0) call errore(' bcast_complex ', ' error in mpi_allreduce 2 ', ierr)
                                        end do
                                        blksiz = mod(n, bufdim)
                                        if (blksiz > 0) then
                                            istart = nblk*bufdim + 1
                                            if (myid .eq. root) then
                                                buff_reducec(1:blksiz) = array(istart:istart + blksiz - 1)
                                            else
                                                buff_reducec(1:blksiz) = (0.d0, 0.d0)
                                            end if
                                            call MPI_ALLREDUCE(buff_reducec, array(istart), blksiz, MPI_COMPLEX16, MPI_SUM, gid, ierr)
!   CALL MPI_BCAST( buff_reduce, blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
!             array(istart:istart+blksiz-1)=buff_reduce(1:blksiz)
                                            if (ierr /= 0) call errore(' bcast_complex ', ' error in mpi_allreduce 3 ', ierr)
                                        end if
                                    end if
                                else

                                    if (n <= bufdim) then

                                        buff_reducec(1:n) = array(1:n)
                                        call MPI_BCAST(buff_reducec, n, MPI_COMPLEX16, root, gid, ierr)
                                        array(1:n) = buff_reducec(1:n)
                                        if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_bcast 1 ', ierr)
                                    else
                                        nblk = n/bufdim
                                        blksiz = bufdim
                                        do iblk = 1, nblk
                                            istart = (iblk - 1)*bufdim + 1
                                            buff_reducec(1:blksiz) = array(istart:istart + blksiz - 1)
                                            call MPI_BCAST(buff_reducec, blksiz, MPI_COMPLEX16, root, gid, ierr)
                                            array(istart:istart + blksiz - 1) = buff_reducec(1:blksiz)

                                            if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_bcast 2 ', ierr)
                                        end do
                                        blksiz = mod(n, bufdim)
                                        if (blksiz > 0) then
                                            istart = nblk*bufdim + 1
                                            buff_reducec(1:blksiz) = array(istart:istart + blksiz - 1)
                                            call MPI_BCAST(buff_reducec, blksiz, MPI_COMPLEX16, root, gid, ierr)
                                            array(istart:istart + blksiz - 1) = buff_reducec(1:blksiz)
                                            if (ierr /= 0) call errore(' bcast_real ', ' error in mpi_bcast 3 ', ierr)
                                        end if
                                    end if

                                end if
#endif

                                return
                            end subroutine BCAST_COMPLEX

!
! ... "reduce"-like subroutines
!
!----------------------------------------------------------------------------
                            subroutine reduce_base_real_to(dim, ps, psout, comm, root)
                                !----------------------------------------------------------------------------
                                !
                                ! ... sums a distributed variable ps(dim) over the processors.
                                ! ... This version uses a fixed-length buffer of appropriate (?) dim
                                ! ...              uses SHMEM if available, MPI otherwhise
                                !
                                use buffer, only: bufdim
                                implicit none
                                !
                                integer, intent(IN) :: dim
                                real*8, intent(IN) :: ps(dim)
                                real*8, intent(OUT) :: psout(dim)
                                integer, intent(IN) :: comm ! communecator
                                integer, intent(IN) :: root ! if root <  0 perform a reduction to all procs
                                ! if root >= 0 perform a reduce only to root proc.
                                !
#if defined (PARALLEL)
                                !
                                include 'mpif.h'
                                !
                                integer :: info, n, nbuf, nproc, myid
                                integer :: maxb
                                !
                                call mpi_comm_size(comm, nproc, info)
                                if (info /= 0) call errore('reduce_base_real_to', 'error in mpi_comm_size', info)

                                call mpi_comm_rank(comm, myid, info)
                                if (info /= 0) call errore('reduce_base_real_to', 'error in mpi_comm_rank', info)

                                maxb = min(bufdim, dim)

                                !
                                if (dim <= 0 .or. nproc < 1) return
                                !
                                if (nproc == 1) then
                                    psout = ps
                                    return
                                end if
                                !
                                ! ... synchronize processes
                                !
#ifdef UNREL
                                call mpi_barrier(comm, info)
!$omp barrier
                                if (info /= 0) call errore('reduce_base_real_to', 'error in mpi_barrier', info)
#endif
                                !
                                nbuf = dim/maxb
                                !
                                do n = 1, nbuf
                                    !
                                    !
                                    if (root >= 0) then
                                        call MPI_REDUCE(ps(1 + (n - 1)*maxb), psout(1 + (n - 1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info)
                                        if (info /= 0) call errore('reduce_base_real_to', 'error in mpi_reduce 1', info)
                                    else
                                        call MPI_ALLREDUCE(ps(1 + (n - 1)*maxb), psout(1 + (n - 1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info)
                                        if (info /= 0) call errore('reduce_base_real_to', 'error in mpi_allreduce 1', info)
                                    end if
                                    !
                                end do
                                !
                                ! ... possible remaining elements < maxb
                                !
                                if ((dim - nbuf*maxb) > 0) then
                                    !
                                    if (root >= 0) then
                                        call MPI_REDUCE(ps(1 + nbuf*maxb), psout(1 + nbuf*maxb), (dim - nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info)
                                        if (info /= 0) call errore('reduce_base_real_to', 'error in mpi_reduce 2', info)
                                    else
                                        call MPI_ALLREDUCE(ps(1 + nbuf*maxb), psout(1 + nbuf*maxb), (dim - nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info)
                                        if (info /= 0) call errore('reduce_base_real_to', 'error in mpi_allreduce 2', info)
                                    end if
                                    !
                                    !
                                end if
                                !
#else
                                psout = ps
                                !
#endif
                                return
                                !
                            end subroutine reduce_base_real_to

!----------------------------------------------------------------------------
                            subroutine reduce_base_complex_to(dim, ps, psout, comm, root)
                                !----------------------------------------------------------------------------
                                !
                                ! ... sums a distributed variable ps(dim) over the processors.
                                ! ... This version uses a fixed-length buffer of appropriate (?) dim
                                ! ...              uses SHMEM if available, MPI otherwhise
                                !
                                use buffer, only: bufdim
                                implicit none
                                !
                                integer, intent(IN) :: dim
                                complex(8), intent(IN) :: ps(dim)
                                complex(8), intent(OUT) :: psout(dim)
                                integer, intent(IN) :: comm ! communicator
                                integer, intent(IN) :: root ! if root <  0 perform a reduction to all procs
                                ! if root >= 0 perform a reduce only to root proc.
                                !
#if defined (PARALLEL)
                                !
                                include 'mpif.h'
                                !
                                integer :: info, n, nbuf, nproc, myid
                                integer :: maxb
                                !
                                call mpi_comm_size(comm, nproc, info)
                                if (info /= 0) call errore('reduce_base_complex_to', 'error in mpi_comm_size', info)

                                call mpi_comm_rank(comm, myid, info)
                                if (info /= 0) call errore('reduce_base_complex_to', 'error in mpi_comm_rank', info)

                                maxb = min(bufdim, dim)
                                !
                                if (dim <= 0 .or. nproc < 1) return
                                !
                                if (nproc == 1) then
                                    psout = ps
                                    return
                                end if
                                !
                                ! ... synchronize processes
                                !
#ifdef UNREL
                                call mpi_barrier(comm, info)
!$omp barrier
                                if (info /= 0) call errore('reduce_base_complex_to', 'error in mpi_barrier', info)
#endif
                                !
                                nbuf = dim/maxb ! >= 1
                                !
                                do n = 1, nbuf
                                    !
                                    !
                                    if (root >= 0) then
                                        call MPI_REDUCE(ps(1 + (n - 1)*maxb), psout(1 + (n - 1)*maxb), maxb, MPI_COMPLEX16, MPI_SUM, root, comm, info)
                                        if (info /= 0) call errore('reduce_base_complex_to', 'error in mpi_reduce 1', info)
                                    else
                                        call MPI_ALLREDUCE(ps(1 + (n - 1)*maxb), psout(1 + (n - 1)*maxb), maxb, MPI_COMPLEX16, MPI_SUM, comm, info)
                                        if (info /= 0) call errore('reduce_base_complex_to', 'error in mpi_allreduce 1', info)
                                    end if
                                    !
                                end do
                                !
                                ! ... possible remaining elements < maxb
                                !
                                if ((dim - nbuf*maxb) > 0) then
                                    !
                                    if (root >= 0) then
                                        call MPI_REDUCE(ps(1 + nbuf*maxb), psout(1 + nbuf*maxb), (dim - nbuf*maxb), MPI_COMPLEX16, MPI_SUM, root, comm, info)
                                        if (info /= 0) call errore('reduce_base_complex_to', 'error in mpi_reduce 2', info)
                                    else
                                        call MPI_ALLREDUCE(ps(1 + nbuf*maxb), psout(1 + nbuf*maxb), (dim - nbuf*maxb), MPI_COMPLEX16, MPI_SUM, comm, info)
                                        if (info /= 0) call errore('reduce_base_complex_to', 'error in mpi_allreduce 2', info)
                                    end if
                                    !
                                    !
                                end if
                                !
#else
                                psout = ps
                                !
#endif

                                return
                                !
                            end subroutine reduce_base_complex_to
!
!
