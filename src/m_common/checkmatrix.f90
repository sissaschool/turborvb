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

subroutine checkmatrix_sparse(nnozero, nnozeron, nelorb, nozero, nozeron&
        &, index_out, rank, nelorb_c)
    implicit none
    integer nnozero, nnozeron, nelorb, nozero(*), nozeron(*), rank, dimmax&
            &, i, ind, icheck, iy, ix, ii, jj, ifound, ierr, nelorb_c
    real*8, dimension(:), allocatable :: sortvect
    integer, dimension(:), allocatable :: indexo, indexn
    logical doloop
    integer index_out(nnozero)
    real*8 val_nozero, val_nozeron
    integer*8 maxindex
#ifdef PARALLEL
    include 'mpif.h'
#endif

    maxindex = int(nelorb_c, 8)*nelorb

    dimmax = max(nnozero, nnozeron)

    if (dimmax .eq. 0) return

    !       write(6,*) ' nnozero nnozeron inside ',nnozero,nnozeron

    if (rank .eq. 0) then

        allocate (sortvect(dimmax), indexo(dimmax), indexn(dimmax))

        indexo = 1
        indexn = 1

        if (nnozero .gt. 0) then
            do i = 1, nnozero
                sortvect(i) = nozero(i) + dble(nelorb)*(nozero(i + nnozero) - 1)
            end do
            call dsortx(sortvect, 1, nnozero, indexo)

            if (sortvect(1) .lt. 1 .or. sortvect(nnozero) .gt. maxindex) then
                if (rank .eq. 0) write (6, *) ' Matrix elements out of range !!! '
#ifdef  PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if

        end if
        if (nnozeron .gt. 0) then
            do i = 1, nnozeron
                sortvect(i) = abs(nozeron(i)) + dble(nelorb)*(nozeron(i + nnozero) - 1)
            end do
            call dsortx(sortvect, 1, nnozeron, indexn)
            if (sortvect(1) .lt. 1 .or. sortvect(nnozero) .gt. maxindex) then
                if (rank .eq. 0) write (6, *) ' Matrix elements out of range !!! '
#ifdef  PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if
        end if
        index_out = 0
        !       search in ascending order
        ifound = 1
        do i = 1, nnozeron
            ii = indexn(i)
            jj = indexo(ifound)
            val_nozero = nozero(jj) + dble(nelorb)*(nozero(jj + nnozero) - 1)
            val_nozeron = abs(nozeron(ii)) + dble(nelorb)*(nozeron(ii + nnozero) - 1)

            do while (ifound .lt. nnozero .and. val_nozero .lt. val_nozeron)
                ifound = ifound + 1
                jj = indexo(ifound)
                val_nozero = nozero(jj) + dble(nelorb)*(nozero(jj + nnozero) - 1)
            end do
!           jj = indexo(ifound)
!           val_nozero=nozero(jj)+dble(nelorb)*(nozero(jj+nnozero)-1)
            val_nozeron = abs(nozeron(ii)) + dble(nelorb)*(nozeron(ii + nnozero) - 1)
            if (val_nozero .eq. val_nozeron) then
                !         scaleout(ii)=scalein(jj)
                if (nozeron(ii) .lt. 0) then
                    index_out(jj) = -ii
                else
                    index_out(jj) = ii
                end if
            else
                iy = nozeron(ii + nnozero)
                ix = abs(nozeron(ii))
                if (rank .eq. 0) write (6, *) 'Not found in given lambda table: ix, iy', ix, iy
            end if
        end do

        deallocate (sortvect, indexo, indexn)

    end if
#ifdef PARALLEL
    call bcast_integer(index_out, nnozero, 0, MPI_COMM_WORLD)
#endif

    return
end

subroutine checkmatrix(nnozero, nnozeron, nelorb, nozero, nozeron&
        &, index_out, rank, nelorb_c)
    implicit none
    integer nnozero, nnozeron, nelorb, nozero(*), nozeron(*), rank, dimmax&
            &, i, ind, icheck, iy, ix, ii, jj, ifound, ierr, maxindex, nelorb_c
    real*8, dimension(:), allocatable :: sortvect
    integer, dimension(:), allocatable :: indexo, indexn
    logical doloop
    integer index_out(nnozero)
#ifdef PARALLEL
    include 'mpif.h'
#endif

    maxindex = nelorb_c*nelorb

    dimmax = max(nnozero, nnozeron)

    if (dimmax .eq. 0) return

    !       write(6,*) ' nnozero nnozeron inside ',nnozero,nnozeron

    if (rank .eq. 0) then

        allocate (sortvect(dimmax), indexo(dimmax), indexn(dimmax))

        indexo = 1
        indexn = 1

        if (nnozero .gt. 0) then
            sortvect(1:nnozero) = nozero(1:nnozero)
            call dsortx(sortvect, 1, nnozero, indexo)

            if (sortvect(1) .lt. 1 .or. sortvect(nnozero) .gt. maxindex) then
                if (rank .eq. 0) write (6, *) ' Matrix elements out of range !!! '
#ifdef  PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if

        end if

        if (nnozeron .gt. 0) then
            sortvect(1:nnozeron) = abs(nozeron(1:nnozeron))
            call dsortx(sortvect, 1, nnozeron, indexn)
            if (sortvect(1) .lt. 1 .or. sortvect(nnozero) .gt. maxindex) then
                if (rank .eq. 0) write (6, *) ' Matrix elements out of range !!! '
#ifdef  PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if
        end if
        index_out = 0
        !       search in ascending order
        ifound = 1
        do i = 1, nnozeron
            ii = indexn(i)
            do while (ifound .lt. nnozero .and. nozero(indexo(ifound)) .lt. abs(nozeron(ii)))
                ifound = ifound + 1
            end do
            jj = indexo(ifound)
            if (nozero(jj) .eq. abs(nozeron(ii))) then
                !         scaleout(ii)=scalein(jj)
                if (nozeron(ii) .lt. 0) then
                    index_out(jj) = -ii
                else
                    index_out(jj) = ii
                end if
            else
                iy = (abs(nozeron(ii)) - 1)/nelorb + 1
                ix = abs(nozeron(ii)) - (iy - 1)*nelorb
                if (rank .eq. 0) write (6, *) 'Not found in given lambda table: ix, iy', ix, iy
            end if
        end do

        deallocate (sortvect, indexo, indexn)

    end if
#ifdef PARALLEL
    call bcast_integer(index_out, nnozero, 0, MPI_COMM_WORLD)
#endif

    return
end

subroutine checkrepeat(nnozero, nelorb, nozero, rank)
    implicit none
    integer nnozero, nelorb, nozero(*), rank&
            &, i, ind, iy, ix, ii, jj, ifound, j, ierr
    real*8, dimension(:), allocatable :: sortvect
    integer, dimension(:), allocatable :: index
    logical iflagerr
#ifdef PARALLEL
    include 'mpif.h'
#endif
    if (nnozero .eq. 0) return

    if (rank .eq. 0) then

        allocate (sortvect(nnozero), index(nnozero))
        iflagerr = .false.

        sortvect(1:nnozero) = abs(nozero(1:nnozero))
        call dsortx(sortvect, 1, nnozero, index)

        do i = 2, nnozero
            if (abs(nozero(index(i))) .eq. abs(nozero(index(i - 1)))) then
                ifound = 2
                if (i .lt. nnozero) then
                    j = i + 1
                    do while (abs(nozero(index(j))) .eq. abs(nozero(index(j - 1))) .and. j .lt. nnozero)
                        ifound = ifound + 1
                        j = j + 1
                    end do
                    if (j .eq. nnozero) then
                        if (abs(nozero((index(j)))) .eq. abs(nozero(index(j - 1)))) ifound = ifound + 1
                    end if
                end if
                ii = index(i)
                iy = (abs(nozero(ii)) - 1)/nelorb + 1
                ix = abs(nozero(ii)) - (iy - 1)*nelorb
                write (6, *) ' Repeated matrix elem. ', ix, iy
                write (6, *) ' # times   ', ifound
                iflagerr = .true.
            end if
        end do

        deallocate (sortvect, index)

    end if
#ifdef PARALLEL
    call mpi_bcast(iflagerr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (iflagerr) call mpi_finalize(ierr)
#endif
    if (iflagerr) stop
    return
end

