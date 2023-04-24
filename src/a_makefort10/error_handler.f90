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

    subroutine init_error(routine_name)
        implicit none
        character(len=*), intent(in) :: routine_name

        allocate (routine_chain)

        routine_chain%routine_name = routine_name
        nullify (routine_chain%previous_link)

        return
    end subroutine init_error

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

    subroutine chop_name
        implicit none
        type(chain), pointer :: chopped_chain

        chopped_chain => routine_chain%previous_link
        deallocate (routine_chain)
        routine_chain => chopped_chain

        return
    end subroutine chop_name

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

    subroutine warning(message)
        character(len=*), intent(in) :: message
        call error_mem(message, -1)
        return
    end subroutine warning

end module error_handler
