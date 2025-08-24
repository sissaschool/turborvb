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

/**
 * @brief MPI-IO module for parallel file operations
 *
 * This module provides a high-level interface for parallel file I/O operations
 * using MPI-IO. It includes a file object type and subroutines for opening,
 * reading, writing, and closing files in parallel environments.
 *
 * @details
 * The module provides the following functionality:
 * - File object type with MPI-IO metadata
 * - File opening with specified access modes
 * - File view creation for parallel data distribution
 * - File pointer management and positioning
 * - File closing and cleanup
 * - Error handling for MPI-IO operations
 *
 * Key features:
 * - Parallel file access with collective operations
 * - Custom file views for data distribution
 * - Automatic error handling and abort on failure
 * - Support for both serial and parallel builds
 * - File pointer tracking and management
 *
 * @note
 * - Only available when PARALLEL is defined
 * - Serial version provides dummy type for compatibility
 * - All operations include error checking
 * - File views are automatically committed and freed
 *
 * @author TurboRVB group
 * @date 2022
 */
module mpiio
    implicit none

#ifdef PARALLEL

    include 'mpif.h'

    /**
     * @brief File object type for MPI-IO operations
     *
     * This derived type encapsulates all the information needed for
     * MPI-IO file operations, including file handle, communicator,
     * and view information.
     *
     * @details
     * The type contains:
     * - name: File name string
     * - fp: MPI file pointer
     * - comm: MPI communicator for collective operations
     * - nproc, rank: Process count and rank information
     * - view: MPI datatype for file view
     * - array_of_blocklengths, array_of_displacements, array_of_types: View parameters
     * - etype: Elementary datatype
     * - disp: File displacement for positioning
     *
     * @note
     * - Used by all MPI-IO subroutines in this module
     * - Automatically manages MPI datatype lifecycle
     * - Supports custom file views for parallel access patterns
     */
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
    /**
     * @brief Open a file for parallel I/O operations
     *
     * This subroutine opens a file for parallel I/O using MPI-IO.
     * It initializes the file object with the necessary metadata
     * and opens the file with the specified access mode.
     *
     * @param[in] input_comm MPI communicator for collective operations
     * @param[in] filename Name of the file to open
     * @param[in] amode Access mode (e.g., MPI_MODE_RDONLY, MPI_MODE_WRONLY, MPI_MODE_CREATE)
     * @param[out] myfile File object to be initialized
     *
     * @details
     * The subroutine:
     * 1. Stores the communicator and file name
     * 2. Gets process count and rank information
     * 3. Opens the file with MPI_File_open
     * 4. Initializes the view to MPI_DATATYPE_NULL
     * 5. Reports success or failure
     *
     * @note
     * - Collective operation (all processes must call)
     * - Aborts the program on file open failure
     * - Prints status messages from rank 0
     * - File view must be created separately if needed
     *
     * @see mpiio_file_close(), mpiio_file_create_view()
     */
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

    /**
     * @brief Get current file pointer position
     *
     * This subroutine retrieves the current position of the file pointer
     * and stores it in the file object's disp field. This position is
     * used by reset_file_view to restore the file view.
     *
     * @param[in] myfile File object containing the file handle
     *
     * @details
     * The subroutine calls MPI_File_get_position to obtain the current
     * file pointer position and stores it in myfile%disp for later use
     * in reset_file_view.
     *
     * @note
     * - Individual operation (each process can call independently)
     * - Aborts the program on failure
     * - Position is stored in myfile%disp
     * - Used in conjunction with reset_file_view
     *
     * @see mpiio_file_reset_view()
     */
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

    /**
     * @brief Create a custom file view for parallel data distribution
     *
     * This subroutine creates a custom MPI datatype for the file view
     * that defines how data is distributed across processes. It sets up
     * a structured view where each process has a specific region of the file.
     *
     * @param[in,out] myfile File object to be modified
     * @param[in] size_of_data Size of data per process
     * @param[in] datatype MPI datatype for the data elements
     *
     * @details
     * The subroutine creates a structured file view with:
     * - One element of the specified datatype per process
     * - Displacement based on process rank and data size
     * - Upper bound to define the view extent
     * - Automatic commitment of the datatype
     *
     * The view structure:
     * - Lower bound: 0
     * - Data block: size_of_data elements starting at rank*record_size
     * - Upper bound: nproc*record_size
     *
     * @note
     * - Collective operation (all processes must call)
     * - Aborts the program on datatype creation failure
     * - Automatically commits the datatype
     * - Previous view is automatically freed
     *
     * @see mpiio_file_reset_view(), mpiio_file_close()
     */
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

    /**
     * @brief Set file size to zero (truncate file)
     *
     * This subroutine truncates the file to zero size, effectively
     * clearing its contents. It also resets the file displacement
     * to zero.
     *
     * @param[in,out] myfile File object containing the file handle
     *
     * @details
     * The subroutine calls MPI_File_set_size with size zero to
     * truncate the file and resets myfile%disp to zero.
     *
     * @note
     * - Collective operation (all processes must call)
     * - Aborts the program on failure
     * - Useful for clearing file contents before writing
     * - Resets file pointer position
     *
     * @see mpiio_file_open(), mpiio_file_reset_view()
     */
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

    /**
     * @brief Reset file view to the stored displacement
     *
     * This subroutine sets the file view using the stored displacement
     * in myfile%disp. It is typically called after repositioning the
     * file pointer to restore the custom view.
     *
     * @param[in,out] myfile File object containing file handle and view
     *
     * @details
     * The subroutine calls MPI_File_set_view with the stored displacement
     * and the custom view datatype to restore the file view for parallel
     * I/O operations.
     *
     * @note
     * - Individual operation (each process can call independently)
     * - Aborts the program on failure
     * - Uses stored displacement from get_disp
     * - Requires a valid view datatype
     *
     * @see mpiio_file_get_disp(), mpiio_file_create_view()
     */
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

    /**
     * @brief Close file and cleanup MPI-IO resources
     *
     * This subroutine closes the file and frees all associated MPI-IO
     * resources, including the custom view datatype.
     *
     * @param[in,out] myfile File object to be closed and cleaned up
     *
     * @details
     * The subroutine:
     * 1. Frees the custom view datatype
     * 2. Closes the file
     * 3. Reports success or failure
     * 4. Aborts the program on failure
     *
     * @note
     * - Collective operation (all processes must call)
     * - Automatically frees the view datatype
     * - Aborts the program on close failure
     * - Prints status messages from rank 0
     *
     * @see mpiio_file_open(), mpiio_file_create_view()
     */
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

    /**
     * @brief Dummy file object type for serial builds
     *
     * This is a dummy type provided for compatibility when the code
     * is compiled without MPI support. It allows the same interface
     * to be used in both serial and parallel builds.
     */
    type file_obj
        ! dummy type for serial version
    end type file_obj

#endif

end module mpiio
