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

!> @brief Error handling and reporting utilities for TurboRVB
!>
!> This module provides comprehensive error handling and reporting functionality
!> for the TurboRVB quantum Monte Carlo code. It includes subroutines for
!> warning messages, fatal error termination, and parallel error synchronization.
!>
!> The module supports both serial and parallel execution modes, with proper
!> MPI coordination for error handling across multiple processes.
!>
!> @author TurboRVB group
!> @version 1.0
!> @date 2022

!> @brief Report errors and warnings with MPI-aware handling
!>
!> This subroutine provides centralized error reporting for TurboRVB calculations.
!> It handles both warnings (ierror < 0) and fatal errors (ierror > 0) with
!> appropriate MPI coordination for parallel execution.
!>
!> @param[in] routin Name of the calling routine for error identification
!> @param[in] messag Error or warning message to display
!> @param[in] ierror Error code:
!>                   - ierror < 0: Warning (program continues)
!>                   - ierror = 0: No error (subroutine returns immediately)
!>                   - ierror > 0: Fatal error (program terminates)
!> @param[in] rank MPI rank of the calling process
!>
!> @details
!> The subroutine performs the following operations:
!> 1. For warnings (ierror < 0):
!>    - Displays warning message with routine name and error code
!>    - Only rank 0 outputs the message
!>    - Program continues execution
!> 2. For fatal errors (ierror > 0):
!>    - Displays error message with routine name and error code
!>    - For ierror = 3, all ranks output their rank number
!>    - Otherwise, only rank 0 outputs the error
!>    - Terminates all MPI processes and stops the program
!>
!> @note The subroutine uses formatted output with separator lines for
!>       clear visual distinction of error messages
!> @note For parallel execution, MPI_ABORT is called to ensure all processes
!>       terminate when a fatal error occurs
!> @note Error code 3 is special and causes all ranks to output their rank number
subroutine error(routin, messag, ierror, rank)
    use allio, only: iflagerr
    implicit none
    ! the name of the calling routine
    character(*) :: routin
    ! the output message
    character(*) :: messag
    ! the error flag
    integer :: ierror, ierr, info
    integer rank, ierror_mpi

#ifdef PARALLEL
    include 'mpif.h'
#endif

    if (ierror == 0) return

    if (ierror < 0) then
        if (rank .eq. 0) then
            write (*, *) ' '
            write (*, '(1x,78(''%''))')
            write (*, '(5x,''from '',a,'' : Warning #'',i10)') routin, ierror
            write (*, '(5x,a)') 'Warning  '//messag
            write (*, '(1x,78(''%''))')
        end if
    else
        if (rank .eq. 0 .or. ierror .eq. 3) then
            if (ierror .eq. 3) then
                write (*, *) ' From proc # ', rank
            else
                write (*, *) ' From Master '
            end if
            !       write(*,'(1x,78(''%''))')
            write (*, '(5x,''from '',a,'' : ERROR #'',i10)') routin, ierror
            write (*, '(5x,a)') 'ERROR  '//messag
            !       write(*,'(1x,78(''%''))')
        end if
    end if

    if (ierror < 0) then
        ! simple warning, the program goes on.
        return

    else
        !      ! increment error counter
        !       iflagerr=iflagerr+1

        !   Fatal error, abort all MPI processes and stop the program.
#ifdef PARALLEL
        call mpi_abort(MPI_COMM_WORLD, 0, ierr)
#endif
        stop

    end if

end subroutine error

!> @brief Simple error reporting subroutine for serial execution
!>
!> This subroutine provides basic error reporting functionality for serial
!> execution mode. It displays error information and terminates the program
!> when an error occurs.
!>
!> @param[in] a First error message string
!> @param[in] b Second error message string  
!> @param[in] ierr Error code (program terminates if > 0)
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Checks if ierr > 0 (error condition)
!> 2. If error exists, outputs both message strings and error code
!> 3. Terminates program execution
!>
!> @note This subroutine is designed for serial execution and does not
!>       include MPI coordination
!> @note The subroutine returns immediately if ierr <= 0 (no error)
subroutine errore(a, b, ierr)
    implicit none
    character(LEN=*) :: A
    character(LEN=*) :: B
    integer :: IERR

    if (ierr <= 0) return

    write (6, *) A
    write (6, *) B
    write (6, *) IERR
    stop
end subroutine errore

!> @brief Global error checking and synchronization for parallel execution
!>
!> This subroutine provides global error checking across all MPI processes
!> and ensures proper synchronization for error handling in parallel
!> quantum Monte Carlo calculations.
!>
!> @param[in] iflagerr Local error flag for current process
!> @param[in] rank MPI rank of the calling process
!> @param[in] messag Error message to display if global error exists
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Performs global reduction of error flags across all MPI processes
!> 2. If any process has an error (iflagerrall != 0):
!>    - Rank 0 displays the error message
!>    - For unreliable networks, synchronizes with barrier
!>    - Frees MPI communicators if k-point averaging is active
!>    - Finalizes MPI and terminates all processes
!> 3. If no global errors, subroutine returns normally
!>
!> @note This subroutine must be called by ALL MPI processes, not just rank 0
!> @note Do not use this subroutine inside rank.eq.0 regions
!> @note The subroutine handles both reliable and unreliable network configurations
!> @note For k-point averaging calculations, it properly cleans up MPI communicators
subroutine checkiflagerr(iflagerr, rank, messag)
    use kpoints_mod, only: kaverage
    use allio, only: commcolrep_mpi, commrep_mpi
    implicit none
    integer ierr, iflagerrall, iflagerr, rank
    character(*) :: messag

#ifdef PARALLEL
    include 'mpif.h'
#endif

#ifdef PARALLEL
    call mpi_allreduce(iflagerr, iflagerrall, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    iflagerrall = iflagerr
#endif

    if (iflagerrall .ne. 0) then
        if (rank .eq. 0) write (6, *) messag
#ifdef PARALLEL
#ifdef UNREL
!   For unreliable  networks.
        call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
        if (kaverage) then
            call mpi_comm_free(commrep_mpi, ierr)
            call mpi_comm_free(commcolrep_mpi, ierr)
        end if
        call mpi_finalize(ierr)
#endif
        stop
    end if

    return

end subroutine checkiflagerr

