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

module mpiio
    implicit none

#ifdef PARALLEL

    include 'mpif.h'

    type file_obj
        ! file name
        character(len=200) name
        ! file pointer
        integer fp
        ! communicator
        integer comm
        ! mpi setup
        integer nproc, rank
        ! file view
        integer view
        ! view setup
        integer array_of_blocklengths(3)
        integer(kind=MPI_ADDRESS_KIND) array_of_displacements(3)
        integer array_of_types(3), etype
        ! view upper lower boundary
        !integer(kind=MPI_ADDRESS_KIND) lb, ub
        ! first displacement
        integer(kind=MPI_OFFSET_KIND) disp
    end type file_obj

contains
    subroutine mpiio_file_open(input_comm, filename, amode, myfile)
        implicit none

        integer, intent(in) :: input_comm, amode
        character(len=*), intent(in) :: filename
        type(file_obj), intent(out) :: myfile
        integer :: ierr

        myfile%comm = input_comm
        myfile%name = trim(filename)
        call MPI_Comm_size(myfile%comm, myfile%nproc, ierr)
        call MPI_Comm_rank(myfile%comm, myfile%rank, ierr)

        call MPI_File_open(myfile%comm, trim(filename), amode, MPI_INFO_NULL, myfile%fp, ierr)
        if (ierr .ne. MPI_SUCCESS) then
            if (myfile%rank .eq. 0) write (6, *) "ERROR: Fail to open file ", trim(filename), " with MPI-IO!"
            call MPI_Abort(myfile%comm, ierr)
        end if
        if (myfile%rank .eq. 0) write (6, *) "Successfully open file ", trim(myfile%name), " with MPI-IO!"

        myfile%view = MPI_DATATYPE_NULL
    end subroutine mpiio_file_open

    subroutine mpiio_file_get_disp(myfile)
        ! This subroutine update the file displacement myfile%disp used by reset_file_view
        implicit none

        type(file_obj), intent(in) :: myfile

        integer :: ierr

        call MPI_File_get_position(myfile%fp, myfile%disp, ierr)
        if (ierr .ne. MPI_SUCCESS) then
            if (myfile%rank .eq. 0) write (6, *) "ERROR: Fail to get file pointer position ", trim(myfile%name), " with MPI-IO!"
            call MPI_Abort(myfile%comm, ierr)
        end if
        !write(6,*) "my disp", myfile%disp, ", rank ", myfile%rank
    end subroutine mpiio_file_get_disp

    subroutine mpiio_file_create_view(myfile, size_of_data, datatype)
        implicit none

        type(file_obj), intent(inout) :: myfile
        integer, intent(in) :: size_of_data, datatype
        integer :: ierr, unit_size

        integer(kind=MPI_ADDRESS_KIND) record_size
        integer(kind=MPI_OFFSET_KIND) disp_zero

        disp_zero = 0
        myfile%etype = datatype
        call MPI_Type_size(myfile%etype, unit_size, ierr)
        !write(6,*) "Ye here", size_of_data, unit_size
        record_size = unit_size*size_of_data
        !write(6,*) "Ye record_size", record_size, myfile%rank, myfile%nproc
        myfile%array_of_blocklengths = (/1, size_of_data, 1/)
        myfile%array_of_displacements = (/disp_zero, myfile%rank*record_size, record_size*myfile%nproc/)
        myfile%array_of_types = (/MPI_LB, myfile%etype, MPI_UB/)

        !write(6,*) "Ye array", myfile%array_of_displacements

        call MPI_Type_create_struct(3, myfile%array_of_blocklengths, myfile%array_of_displacements, &
                                    myfile%array_of_types, myfile%view, ierr)
        if (ierr .ne. MPI_SUCCESS) then
            write (*, *) 'MPI_Type_create_struct failed'
            call MPI_Abort(myfile%comm, ierr)
        end if

        call MPI_Type_commit(myfile%view, ierr)
        if (ierr .ne. MPI_SUCCESS) then
            write (*, *) 'MPI_Type_myfile%commit failed'
            call MPI_Abort(myfile%comm, ierr)
        end if

        !call MPI_Type_get_extent(myfile%view, myfile%lb, myfile%ub, ierr)
        !if (ierr.ne.MPI_SUCCESS) then
        !  write (*,*) 'MPI_Type_extent failed'
        !  call MPI_Abort(myfile%comm, ierr)
        !else
        !    write(*,*) 'extent = ', myfile%lb, myfile%ub
        !endif
    end subroutine mpiio_file_create_view

    subroutine mpiio_file_set_zero(myfile)
        implicit none

        type(file_obj), intent(inout) :: myfile
        integer(kind=MPI_OFFSET_KIND) zero
        integer :: ierr

        zero = 0
        call MPI_File_set_size(myfile%fp, zero, ierr)
        if (ierr .ne. MPI_SUCCESS) then
            write (*, *) 'MPI_FILE_SET_ZERO failed'
            call MPI_Abort(myfile%comm, ierr)
        end if
        myfile%disp = zero
    end subroutine mpiio_file_set_zero

    subroutine mpiio_file_reset_view(myfile)
        implicit none

        type(file_obj), intent(inout) :: myfile
        integer :: ierr

        call MPI_File_set_view(myfile%fp, myfile%disp, myfile%etype, myfile%view, &
                               'native', MPI_INFO_NULL, ierr)
        if (ierr .ne. MPI_SUCCESS) then
            write (*, *) 'MPI_FILE_SET_VIEW failed'
            call MPI_Abort(myfile%comm, ierr)
        end if
    end subroutine mpiio_file_reset_view

    subroutine mpiio_file_close(myfile)
        implicit none

        type(file_obj), intent(inout) :: myfile
        integer :: ierr

        call MPI_Type_free(myfile%view, ierr)
        call MPI_File_close(myfile%fp, ierr)
        if (ierr .ne. MPI_SUCCESS) then
            write (*, *) 'MPI_FILE_CLOSE failed'
            call MPI_Abort(myfile%comm, ierr)
        end if
        if (myfile%rank .eq. 0) write (*, *) "Successfully closed file ", trim(myfile%name), " with MPI-IO!"
    end subroutine mpiio_file_close

#else

    type file_obj
        ! dummy type for serial version
    end type file_obj

#endif

end module mpiio
