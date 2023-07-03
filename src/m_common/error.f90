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

