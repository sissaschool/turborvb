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

!> @brief MPI sub-communicator management for parallel quantum Monte Carlo
!>
!> This module provides functionality for creating and managing MPI
!> sub-communicators used in parallel quantum Monte Carlo calculations.
!> It allows for dynamic creation of smaller communication groups
!> from larger parent communicators, enabling efficient parallel
!> processing with different group sizes.
!>
!> The module includes:
!> - Derived type for sub-communicator information
!> - Subroutine for creating sub-communicators
!> - Subroutine for freeing sub-communicators
!> - Support for both parallel and serial execution modes
!>
!> @author TurboRVB group
!> @version 1.0
!> @date 2022

module sub_comm
    implicit none
    
    !> @brief Derived type for MPI sub-communicator information
    !> @details Contains all necessary information for managing a sub-communicator:
    !> - rank: Rank of current process in the sub-communicator
    !> - nproc: Number of processes in the sub-communicator
    !> - comm: MPI communicator handle for the sub-communicator
    !> - parent: Parent communicator handle
    !> - yesin: Logical flag indicating if current process is in the sub-communicator
    type :: mpi_sub_comm
        integer rank, nproc
        integer comm, parent
        logical yesin
    end type

contains
    !> @brief Create a new MPI sub-communicator from a parent communicator
    !>
    !> This subroutine creates a new MPI sub-communicator containing a subset
    !> of processes from the parent communicator. It is useful for creating
    !> smaller communication groups for specific parallel operations.
    !>
    !> @param[in] parent_comm Parent MPI communicator handle
    !> @param[in] new_size Number of processes in the new sub-communicator
    !> @param[out] child Sub-communicator information structure
    !> @param[out] ierror MPI error code
    !>
    !> @details
    !> The subroutine performs the following operations:
    !> 1. Allocates rank array for the new communicator
    !> 2. Assigns ranks 0 to new_size-1 to the sub-communicator
    !> 3. Creates MPI group from the parent communicator
    !> 4. Creates new group with specified ranks
    !> 5. Creates new communicator from the group
    !> 6. Determines rank and membership in the new communicator
    !> 7. Frees temporary groups and deallocates arrays
    !>
    !> @note In serial mode, new_size is forced to 1
    !> @note Processes not included in the sub-communicator have yesin = .false.
    !> @note The subroutine handles both parallel and serial execution modes
    !> @note Used for creating specialized communication groups in parallel QMC
    subroutine mpi_sub_comm_create(parent_comm, new_size, child, ierror)
        implicit none
#ifdef PARALLEL
        include 'mpif.h'
#endif
        integer parent_comm, new_size, ierror
        type(mpi_sub_comm) :: child

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

    !> @brief Free an MPI sub-communicator
    !>
    !> This subroutine properly frees an MPI sub-communicator and its
    !> associated resources. It ensures clean deallocation of MPI
    !> communicator handles to prevent memory leaks.
    !>
    !> @param[in,out] child Sub-communicator structure to be freed
    !> @param[out] ierror MPI error code
    !>
    !> @details
    !> The subroutine performs the following operations:
    !> 1. Checks if the current process is a member of the sub-communicator
    !> 2. If yesin = .true., calls MPI_COMM_FREE to deallocate the communicator
    !> 3. Handles both parallel and serial execution modes
    !>
    !> @note Only processes that are members of the sub-communicator
!>       participate in the MPI_COMM_FREE call
!> @note The subroutine is safe to call even if the communicator
!>       was not properly created
!> @note Used for cleanup after parallel operations are complete
    subroutine mpi_sub_comm_free(child, ierror)
        implicit none
        integer ierror
        type(mpi_sub_comm) :: child
#ifdef PARALLEL
        include 'mpif.h'

        if (child%yesin) call mpi_comm_free(child%comm, ierror)
#endif
        return
    end subroutine mpi_sub_comm_free
end module sub_comm
