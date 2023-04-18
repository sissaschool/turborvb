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

program communicate
    use mpi
    implicit none

    integer :: ierr, rank, num_procs, status(MPI_STATUS_SIZE)
    !Variabili ausiliarie
    integer :: i, j
    integer :: m, in_procs, req, ipc
    !Dimension of the submatrices (new and final)
    integer :: sizei, last, sizef, lastf
    real(8), allocatable :: a(:, :), s(:, :), b(:, :)

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

    ipc = 2
    m = 100
    in_procs = 3
    !  num_procs=9
    !Size of the blocks to be sent
    call calc_dime(m, sizei, last, rank, in_procs)
    if (rank .lt. in_procs) then
        allocate (S(ipc*m, sizei))
    else
        allocate (S(ipc*m, 0))
    end if

    !Filling the matrix
    if (rank .eq. 0) then
        allocate (a(ipc*m, m))
        do j = 1, m
            do i = 1, ipc*m
                call random_number(a(i, j))
            end do
        end do
        !Sending the matrix to the first in_procs processors
        do i = 0, in_procs - 1
            if (i .eq. 0) then
                s(:, 1:sizei) = a(:, 1:sizei)
            else if (i .lt. in_procs - 1) then
                call mpi_send(a(1, i*sizei + 1), ipc*m*sizei, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, ierr)
            else
                call mpi_send(a(1, i*sizei + 1), ipc*m*last, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, ierr)
            end if
        end do
    end if
    if (rank .gt. 0) then
        if (rank .lt. in_procs - 1) then
            call mpi_recv(S(1, 1), ipc*m*sizei, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, status, ierr)
        else if (rank .eq. in_procs - 1) then
            call mpi_recv(S(1, 1), ipc*m*last, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, status, ierr)
        end if
    end if
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) write (6, *) "Completed the distribution of the initial matrix"

    !Size of the final blocks
    call calc_dime(m, sizef, lastf, rank, num_procs)
    allocate (b(ipc*m, sizef))
    !  write (6,*) rank, "size, last, sizef, lastf", size, last, sizef,lastf

    call distribute_column(m, ipc, rank, in_procs, num_procs, s, b)

    deallocate (s)
    allocate (s(ipc*m, sizef))
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (rank .eq. 0) then
        do i = 0, num_procs - 1
            if (i .eq. 0) then
                s(:, 1:sizef) = a(:, 1:sizef)
            else if (i .lt. num_procs - 1) then
                call mpi_send(a(1, i*sizef + 1), ipc*m*sizef, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, ierr)
            else
                call mpi_send(a(1, i*sizef + 1), ipc*m*lastf, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, ierr)
            end if
        end do
    end if
    if (rank .gt. 0) then
        if (rank .lt. num_procs - 1) then
            call mpi_recv(S(1, 1), ipc*m*sizef, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, status, ierr)
        else if (rank .eq. num_procs - 1) then
            call mpi_recv(S(1, 1), ipc*m*lastf, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, status, ierr)
        end if
    end if
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (rank .eq. 0) write (6, *) "Completed the distribution of the checking matrix"

    if (rank .eq. num_procs - 1) then
        do i = 1, lastf
            do j = 1, ipc*m
                if (s(j, i) .ne. b(j, i)) write (6, *) "ERROR", rank, j, i
            end do
        end do
    else
        do i = 1, sizef
            do j = 1, ipc*m
                if (s(j, i) .ne. b(j, i)) write (6, *) "ERROR", rank, j, i
            end do
        end do
    end if

    call mpi_finalize(ierr)
    stop

contains

    !C.G.
    !This subroutine reshuffle a matrix divided between the processors by column.
    !The initial configuration is over the first in_procs processes and is stored in
    !the elements S of dimension (ipc*m,sizei), the final configuration is over num_procs
    !processors and is stored in b(ipc*m,sizef). Last and lastf are the dimension of the
    !last blocks of columns of the matrix
    subroutine distribute_column(m, ipc, rank, in_procs, out_procs, s, b)
        use mpi
        implicit none
        integer :: rank, out_procs, in_procs, ierr
        integer :: m, sizei, sizef, last, lastf, ipc
        real(8) :: s(m*ipc, *), b(m*ipc, *)
        !Starting/final process to send/recv
        integer :: start_s, end_s, start_r, end_r
        integer :: first_el, last_el
        integer :: i, j, count
        integer, allocatable :: req(:), status(:, :)

        call calc_dime(m, sizei, last, rank, in_procs)
        call calc_dime(m, sizef, lastf, rank, out_procs)

        allocate (req(sizei/sizef + sizef/sizei + 4), status(MPI_STATUS_SIZE, sizei/sizef + sizef/sizei + 4))

        count = 0
        if (rank .lt. out_procs) then
            start_r = sizef*rank/sizei
            first_el = mod(sizef*rank, sizei)
            if (rank .eq. out_procs - 1) then
                end_r = (sizef*rank + lastf - 1)/sizei
                last_el = lastf - sizei*(end_r - start_r) + first_el
            else
                end_r = (sizef*(rank + 1) - 1)/sizei
                last_el = sizef - sizei*(end_r - start_r) + first_el
            end if
            do i = start_r, end_r
                if (i .eq. start_r .and. i .eq. end_r) then
                    if (rank .eq. out_procs - 1) then
                        count = count + 1
                        call mpi_irecv(b(1, 1), ipc*m*lastf, MPI_DOUBLE, i, rank*out_procs + i, MPI_COMM_WORLD, req(count), ierr)
                    else
                        count = count + 1
                        call mpi_irecv(b(1, 1), ipc*m*sizef, MPI_DOUBLE, i, rank*out_procs + i, MPI_COMM_WORLD, req(count), ierr)
                    end if
                elseif (i .eq. start_r) then
                    count = count + 1
                    call mpi_irecv(b(1, 1), ipc*m*(sizei - first_el), MPI_DOUBLE, i, &
                                   rank*out_procs + i, MPI_COMM_WORLD, req(count), ierr)
                else if (i .eq. end_r) then
                    count = count + 1
                    call mpi_irecv(b(1, 1 + sizei*(i - start_r) - first_el), ipc*m*last_el, &
                                   MPI_DOUBLE, i, rank*out_procs + i, MPI_COMM_WORLD, req(count), ierr)
                else
                    count = count + 1
                    call mpi_irecv(b(1, 1 + sizei*(i - start_r) - first_el), ipc*m*sizei, &
                                   MPI_DOUBLE, i, rank*out_procs + i, MPI_COMM_WORLD, req(count), ierr)
                end if
            end do
        end if

        !Sending the blocks to the processors
        if (rank .lt. in_procs) then
            start_s = sizei*rank/sizef
            first_el = mod(sizei*rank, sizef)
            if (rank .eq. in_procs - 1) then
                end_s = (sizei*rank + last - 1)/sizef
                last_el = last - sizef*(end_s - start_s) + first_el
            else
                end_s = (sizei*(rank + 1) - 1)/sizef
                last_el = sizei - sizef*(end_s - start_s) + first_el
            end if

            do i = start_s, end_s
                if (i .eq. start_s .and. i .eq. end_s) then
                    if (rank .eq. in_procs - 1) then
                        count = count + 1
                        call mpi_isend(s(1, 1), ipc*m*last, MPI_DOUBLE, i, rank + out_procs*i, MPI_COMM_WORLD, req(count), ierr)
                    else
                        count = count + 1
                        call mpi_isend(s(1, 1), ipc*m*sizei, MPI_DOUBLE, i, rank + out_procs*i, MPI_COMM_WORLD, req(count), ierr)
                    end if
                elseif (i .eq. start_s) then
                    count = count + 1
                    call mpi_isend(s(1, 1), ipc*m*(sizef - first_el), MPI_DOUBLE, i, &
                                   rank + out_procs*i, MPI_COMM_WORLD, req(count), ierr)
                else if (i .eq. end_s) then
                    count = count + 1
                    call mpi_isend(s(1, 1 + sizef*(i - start_s) - first_el), ipc*m*last_el, &
                                   MPI_DOUBLE, i, rank + out_procs*i, MPI_COMM_WORLD, req(count), ierr)
                else
                    count = count + 1
                    call mpi_isend(s(1, 1 + sizef*(i - start_s) - first_el), ipc*m*sizef, &
                                   MPI_DOUBLE, i, rank + out_procs*i, MPI_COMM_WORLD, req(count), ierr)
                end if
            end do
        end if

        call MPI_WAITALL(count, req, status, ierr)
        deallocate (status, req)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
    end subroutine distribute_column

    !Calculate the dimensions of the submatrices according to the convention used in
    !distribute_column
    subroutine calc_dime(m, siz, last, rank, num_procs)
        implicit none
        integer :: m, siz, last, rank, num_procs

        siz = m/num_procs
        last = siz
        if (mod(m, num_procs) .ne. 0) then
            siz = siz + 1
            last = mod(m, siz)
        end if
    end subroutine calc_dime

end program COMMUNICATE

