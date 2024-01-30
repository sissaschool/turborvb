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

module sub_comm
    implicit none
    type :: mpi_sub_comm
        integer rank, nproc
        integer comm, parent
        logical yesin
    end type

contains

    subroutine mpi_sub_comm_create(parent_comm, new_size, child, ierror)
        implicit none
#ifdef PARALLEL
        include 'mpif.h'
#endif
        integer, intent(in) :: parent_comm
        integer, intent(inout) :: new_size
        integer, intent(inout) :: ierror ! intent(out)?
        type(mpi_sub_comm), intent(inout) :: child

        integer orig_group, sub_group, i
        integer, dimension(:), allocatable :: ranks
#ifdef PARALLEL
        allocate (ranks(new_size))
#else
        if (new_size .ne. 1) then
            write (6, *) ' Warning changing number of processors to 1 in serial '
            new_size = 1
        end if
        allocate (ranks(new_size))
#endif
        do i = 1, new_size
            ranks(i) = i - 1
        end do

        child%parent = parent_comm
#ifdef PARALLEL
        call mpi_comm_group(parent_comm, orig_group, ierror)
        call mpi_group_incl(orig_group, new_size, ranks, sub_group, ierror)
        call mpi_comm_create(parent_comm, sub_group, child%comm, ierror)
        call mpi_group_rank(sub_group, child%rank, ierror)
        child%nproc = new_size
        if (child%rank .eq. MPI_UNDEFINED) then
            child%yesin = .false.
        else
            child%yesin = .true.
        end if
#else
        child%nproc = new_size
        child%yesin = .true.
#endif

        !    write(6,*) child%rank,"xx",new_size

#ifdef PARALLEL
        call mpi_group_free(orig_group, ierror)
        call mpi_group_free(sub_group, ierror)
#endif
        deallocate (ranks)

        return
    end subroutine mpi_sub_comm_create

    subroutine mpi_sub_comm_free(child, ierror)
        implicit none
        integer, intent(inout) :: ierror ! intent(out)?
        type(mpi_sub_comm), intent(inout) :: child
#ifdef PARALLEL
        include 'mpif.h'

        if (child%yesin) call mpi_comm_free(child%comm, ierror)
#endif
        return
    end subroutine mpi_sub_comm_free
end module sub_comm
