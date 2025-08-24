! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
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
!----------------------------------------------------------------------------
!> @brief Write and handle error messages for parallel and serial execution
!> @details This subroutine writes an error message to output. If ierr <= 0, it does nothing.
!>          If ierr > 0, it prints the error and stops execution. In parallel execution,
!>          it also writes to a crash file and aborts MPI if available.
!> @param[in] calling_routine Name of the calling routine
!> @param[in] message Error message to output
!> @param[in] ierr Error flag (if > 0, triggers stop)
subroutine errore(calling_routine, message, ierr)
    !----------------------------------------------------------------------------
    !
    ! ... This is a simple routine which writes an error message to output:
    ! ... if ierr <= 0 it does nothing,
    ! ... if ierr  > 0 it stops.
    !
    ! ...          **** Important note for parallel execution ***
    !
    ! ... in parallel execution unit 6 is written only by the first node;
    ! ... all other nodes have unit 6 redirected to nothing (/dev/null).
    ! ... As a consequence an error not occurring on the first node
    ! ... will be invisible. For T3E and ORIGIN machines, this problem
    ! ... is solved by writing an error message to unit * instead of 6.
    ! ... Whenever possible (IBM SP machines), we write to the standard
    ! ... error, unit 0 (the message will appear in the error files
    ! ... produced by loadleveler).
    !
    use io_global, only: stdout
    use io_files, only: crashunit, crash_file
    use parallel_include
    !
    implicit none
    !
    character(LEN=*), intent(IN) :: calling_routine, message
    ! the name of the calling calling_routinee
    ! the output messagee
    integer, intent(IN) :: ierr
    ! the error flag
    integer :: mpime, mpierr
    ! the task id
    !
    logical :: exists
    !
    !
    if (ierr <= 0) return
    !
    ! ... the error message is written un the "*" unit
    !
    write (UNIT=*, FMT='(/,1X,78("%"))')
    write (UNIT=*, &
           FMT='(5X,"from ",A," : error #",I10)') calling_routine, ierr
    write (UNIT=*, FMT='(5X,A)') message
    write (UNIT=*, FMT='(1X,78("%"),/)')
    !
#if defined (__PARA) && defined (__AIX)
    !
    ! ... in the case of ibm machines it is also written on the "0" unit
    ! ... which is automatically connected to stderr
    !
    write (UNIT=0, FMT='(/,1X,78("%"))')
    write (UNIT=0, &
           FMT='(5X,"from ",A," : error #",I10)') calling_routine, ierr
    write (UNIT=0, FMT='(5X,A)') message
    write (UNIT=0, FMT='(1X,78("%"),/)')
    !
#endif
    !
    write (*, '("     stopping ...")')
    !
    call flush_unit(stdout)
    !
#if defined (__PARA) && defined (__MPI)
    !
    mpime = 0
    !
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpime, mpierr)
    !
    !  .. write the message to a file and close it before exiting
    !  .. this will prevent loss of information on systems that
    !  .. do not flush the open streams
    !  .. added by C.C.
    !
    open (UNIT=crashunit, FILE=crash_file, &
          POSITION='APPEND', STATUS='UNKNOWN')
    !
    write (UNIT=crashunit, FMT='(/,1X,78("%"))')
    write (UNIT=crashunit, FMT='(5X,"task #",I10)') mpime
    write (UNIT=crashunit, &
           FMT='(5X,"from ",A," : error #",I10)') calling_routine, ierr
    write (UNIT=crashunit, FMT='(5X,A)') message
    write (UNIT=crashunit, FMT='(1X,78("%"),/)')
    !
    close (UNIT=crashunit)
    !
    ! ... try to exit in a smooth way
    !
    call MPI_ABORT(MPI_COMM_WORLD, mpierr)
    !
    call MPI_FINALIZE(mpierr)
    !
#endif
    !
    stop 2
    !
    return
    !
end subroutine errore
!
!----------------------------------------------------------------------
!> @brief Write an informational message from a given routine
!> @details This subroutine writes an info message to output from the specified routine.
!> @param[in] routine Name of the calling routine
!> @param[in] message Informational message to output
subroutine infomsg(routine, message)
    !----------------------------------------------------------------------
    !
    ! ... This is a simple routine which writes an info message
    ! ... from a given routine to output.
    !
    use io_global, only: stdout, ionode
    !
    implicit none
    !
    character(LEN=*) :: routine, message
    ! the name of the calling routine
    ! the output message
    !
    if (ionode) then
        !
        write (stdout, '(5X,"Message from routine ",A,":")') routine
        write (stdout, '(5X,A)') message
        !
    end if
    !
    return
    !
end subroutine infomsg
!
!> @brief Error handling utilities for routine call stack and memory errors
!> @details This module provides utilities for tracking the call stack, reporting
!>          memory errors, and issuing warnings. It supports traceback and
!>          hierarchical error reporting for easier debugging.
module error_handler
    implicit none
    private

    public :: init_error, add_name, chop_name, error_mem, warning

    type chain
        character(len=35) :: routine_name
        type(chain), pointer :: previous_link
    end type chain

    type(chain), pointer :: routine_chain

contains

    !> @brief Initialize the error handler with the first routine name
    !> @param[in] routine_name Name of the first routine in the call stack
    subroutine init_error(routine_name)
        implicit none
        character(len=*), intent(in) :: routine_name

        allocate (routine_chain)

        routine_chain%routine_name = routine_name
        nullify (routine_chain%previous_link)

        return
    end subroutine init_error

    !> @brief Add a routine name to the call stack
    !> @param[in] routine_name Name of the routine to add
    subroutine add_name(routine_name)
        implicit none
        character(len=*), intent(in) :: routine_name
        type(chain), pointer :: new_link

        allocate (new_link)
        new_link%routine_name = routine_name
        new_link%previous_link => routine_chain
        routine_chain => new_link

        return
    end subroutine add_name

    !> @brief Remove the most recent routine name from the call stack
    subroutine chop_name
        implicit none
        type(chain), pointer :: chopped_chain

        chopped_chain => routine_chain%previous_link
        deallocate (routine_chain)
        routine_chain => chopped_chain

        return
    end subroutine chop_name

    !> @brief Recursively print the call stack for error tracing
    !> @param[in] error_code Error code to determine stop or return
    recursive subroutine trace_back(error_code)

        implicit none
        integer :: error_code

        write (unit=*, fmt=*) "   Called by ", routine_chain%routine_name
        if (.not. associated(routine_chain%previous_link)) then
            write (unit=*, fmt=*) &
                " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
            write (unit=*, fmt=*) " "
            if (error_code > 0) then
                stop
            else
                return
            end if
        end if

        routine_chain => routine_chain%previous_link
        call trace_back(error_code)

    end subroutine trace_back

    !> @brief Report a memory error or warning with traceback
    !> @param[in] message Error or warning message
    !> @param[in] error_code (Optional) Error code (default 1: fatal, -1: warning)
    subroutine error_mem(message, error_code)
        character(len=*), intent(in) :: message
        integer, intent(in), optional :: error_code
        integer :: action_code
        type(chain), pointer :: save_chain

        if (present(error_code)) then
            action_code = error_code
        else
            action_code = 1
        end if

        if (action_code /= 0) then
            write (unit=*, fmt=*) " "
            write (unit=*, fmt=*) &
                " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"

            if (action_code > 0) then
                write (unit=*, fmt=*) "   Fatal error in routine `", &
                    trim(routine_chain%routine_name), "': ", message
            else
                write (unit=*, fmt=*) "   Warning from routine `", &
                    trim(routine_chain%routine_name), "': ", message
                save_chain => routine_chain
            end if
            write (unit=*, fmt=*) &
                " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
            routine_chain => routine_chain%previous_link
            call trace_back(action_code)
            routine_chain => save_chain
        end if

        return
    end subroutine error_mem

    !> @brief Issue a warning message with traceback
    !> @param[in] message Warning message
    subroutine warning(message)
        character(len=*), intent(in) :: message
        call error_mem(message, -1)
        return
    end subroutine warning

end module error_handler
