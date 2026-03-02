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

subroutine error(routin, messag, ierror, rank)
    use allio, only: iflagerr
    use logger_io, only: log_warning, log_error
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
        call log_warning(' ')
        call log_warning(' from ', trim(routin), ' : Warning #', ierror)
        call log_warning(' Warning ', trim(messag))
    else
        if (rank .eq. 0 .or. ierror .eq. 3) then
            if (ierror .eq. 3) then
                call log_error(' From proc # ', rank)
            else
                call log_error(' From Master ')
            end if
            call log_error(' from ', trim(routin), ' : ERROR #', ierror)
            call log_error(' ERROR ', trim(messag))
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

subroutine errore(a, b, ierr)
    use logger_io, only: log_info
    implicit none
    character(LEN=*) :: A
    character(LEN=*) :: B
    integer :: IERR

    if (ierr <= 0) return

    call log_info(A)
    call log_info(B)
    call log_info(IERR)
    stop
end subroutine errore

!------------------ checkiflagerr -----------------
! This is a subroutine very easy to use.
! It's been designed to stop the whole program if there's any error.
! It is both parallel and serial supported.
! In parallel, it should be called by all mpi processes not only master.
! REMINDER: Don't use it inside rank.eq.0 region.
!--------------------------------------------------
subroutine checkiflagerr(iflagerr, rank, messag)
    use kpoints_mod, only: kaverage
    use allio, only: commcolrep_mpi, commrep_mpi
    use logger_io, only: log_info
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
        call log_info(messag)
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

