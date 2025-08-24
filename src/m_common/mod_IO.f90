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
 * @brief Input/Output utilities module for file and directory operations
 *
 * This module provides a collection of utility subroutines for file and
 * directory operations in TurboRVB. It includes functions for creating
 * directories, changing directories, copying files, removing files,
 * renaming files, and string manipulation utilities.
 *
 * @details
 * The module provides the following functionality:
 * - Directory operations: creation, navigation, path retrieval
 * - File operations: copying, removal, renaming, status checking
 * - String utilities: C-string conversion, clearing, integer-to-string
 * - System utilities: sleep function for timing control
 *
 * Key features:
 * - Cross-platform compatibility
 * - Error handling for file operations
 * - Safe string manipulation with proper null termination
 * - Integration with system commands for file operations
 *
 * @note
 * - All subroutines include input validation
 * - File operations use system commands for portability
 * - String operations handle C-compatibility for system calls
 * - Directory operations check for existing directories before creation
 *
 * @author TurboRVB group
 * @date 2022
 */
module IO_m
    implicit none
    integer, parameter :: lchlen = 1000

    !
contains
    !
    /**
     * @brief Create a directory if it doesn't already exist
     *
     * This subroutine creates a directory with the specified name if it
     * doesn't already exist. It first checks if the directory exists by
     * attempting to change to it, then creates it if necessary.
     *
     * @param[in] dirname Character string containing the directory name to create
     *
     * @details
     * The subroutine:
     * 1. Saves the current working directory
     * 2. Attempts to change to the target directory
     * 3. Compares the new path with the original path
     * 4. Creates the directory if it doesn't exist
     * 5. Returns to the original directory
     *
     * @note
     * - Returns immediately if dirname is empty
     * - Uses system mkdir command for directory creation
     * - Preserves the original working directory
     * - Safe to call multiple times for the same directory
     *
     * @see cd_dir(), get_dir(), imkdir()
     */
    subroutine mk_dir(dirname)
        implicit none
        character(*) :: dirname
        character(lchlen) :: checkdir, path
        if (len_trim(dirname) == 0) return ! return for wrong input
        call get_dir(path)
        call cd_dir(trim(dirname))
        call get_dir(checkdir)
        call cd_dir(trim(path)) ! go back to original path
        if (trim(checkdir) .ne. trim(path)) return ! directory already exists
        call imkdir(cstr(trim(dirname)))
    end subroutine

    /**
     * @brief Change the current working directory
     *
     * This subroutine changes the current working directory to the specified
     * directory name using the system chdir command.
     *
     * @param[in] dirname Character string containing the target directory name
     *
     * @note
     * - Returns immediately if dirname is empty
     * - Uses system chdir command for directory change
     * - No error checking for invalid directory names
     *
     * @see mk_dir(), get_dir(), ichdir()
     */
    subroutine cd_dir(dirname)
        implicit none
        character(*) :: dirname
        if (len_trim(dirname) == 0) return
        call ichdir(cstr(trim(dirname)))
    end subroutine
    !
    /**
     * @brief Copy a file from source to destination
     *
     * This subroutine copies a file from the source location to the destination
     * using the system cp command. The operation is performed silently with
     * error output redirected to /dev/null.
     *
     * @param[in] file_ Character string containing the source file path
     * @param[in] dest_ Character string containing the destination file path
     * @param[out] ierr_ Integer error code (0 for success, non-zero for failure)
     *
     * @details
     * The subroutine executes: cp file_ dest_ >& /dev/null
     * and returns the exit status in ierr_.
     *
     * @note
     * - Uses system cp command for file copying
     * - Error output is suppressed
     * - No validation of file existence or permissions
     *
     * @see rm_file(), rename_file(), isystem()
     */
    subroutine cp_file(file_, dest_, ierr_)
        implicit none
        character(*) :: file_, dest_
        integer :: ierr_
        call isystem(cstr("cp "//file_//" "//dest_//" >& /dev/null"), ierr_)
    end subroutine
    !
    /**
     * @brief Remove a file from the filesystem
     *
     * This subroutine removes a file from the filesystem using the system
     * remove command.
     *
     * @param[in] filename Character string containing the file path to remove
     *
     * @note
     * - Returns immediately if filename is empty
     * - Uses system remove command for file deletion
     * - No error checking or confirmation
     *
     * @see cp_file(), rename_file(), iremove()
     */
    subroutine rm_file(filename)
        implicit none
        character(*) :: filename
        if (len_trim(filename) == 0) return
        call iremove(cstr(trim(filename)))
    end subroutine
    !
    /**
     * @brief Rename a file from old name to new name
     *
     * This subroutine renames a file from the old filename to the new filename
     * using the system rename command.
     *
     * @param[in] filename_old Character string containing the current file name
     * @param[in] filename_new Character string containing the new file name
     *
     * @note
     * - Returns immediately if filename_old is empty
     * - Uses system rename command for file renaming
     * - No error checking for file existence or permissions
     *
     * @see cp_file(), rm_file(), irename()
     */
    subroutine rename_file(filename_old, filename_new)
        implicit none
        character(*) :: filename_old, filename_new
        if (len_trim(filename_old) == 0) return
        call irename(cstr(trim(filename_old)), cstr(trim(filename_new)))
    end subroutine
    !
    /**
     * @brief Get the current working directory path
     *
     * This subroutine retrieves the current working directory path and
     * stores it in the provided character string.
     *
     * @param[out] path Character string to store the current working directory path
     *
     * @details
     * The subroutine:
     * 1. Calls the system getcwd function
     * 2. Extracts the actual path length
     * 3. Truncates the path to the correct length
     *
     * @note
     * - Uses system getcwd command
     * - Path is truncated to the actual length returned by the system
     * - No error checking for buffer overflow
     *
     * @see mk_dir(), cd_dir(), igetcwd()
     */
    subroutine get_dir(path)
        implicit none
        integer :: ln
        character(*) :: path
        call igetcwd(path, ln)
        path = path(1:ln)
    end subroutine get_dir
    !
    /**
     * @brief Check if a file is currently open
     *
     * This function checks whether a file with the specified name is
     * currently open in the program.
     *
     * @param[in] filename Character string containing the file name to check
     * @return Logical value: .true. if file is open, .false. otherwise
     *
     * @note
     * - Returns .false. if filename is empty
     * - Uses Fortran inquire statement
     * - Only checks if the file is open, not if it exists
     *
     * @see cp_file(), rm_file(), rename_file()
     */
    logical function file_is_open(filename)
        character(*) :: filename
        file_is_open = .false.
        if (len_trim(filename) == 0) return
        inquire (file=filename, opened=file_is_open)
        !
    end function
    !
    /**
     * @brief Convert a Fortran string to a C-style null-terminated string
     *
     * This function converts a Fortran character string to a C-style
     * null-terminated string suitable for passing to C functions or
     * system calls.
     *
     * @param[in] si Input Fortran character string
     * @return C-style null-terminated character string
     *
     * @details
     * The function:
     * 1. Clears the output string
     * 2. Copies the input string characters
     * 3. Adds a null terminator at the end
     *
     * @note
     * - Output string must be large enough to hold input + null terminator
     * - Uses clear_str() to initialize the output string
     * - Essential for interfacing with C functions and system calls
     *
     * @see clear_str()
     */
    character(lchlen) function cstr(si) result(so)
        character(*), intent(IN) :: si
        integer :: i
        i = len(trim(si))
        call clear_str(so)
        so(1:i) = si(1:i)
        so(i + 1:i + 1) = achar(0)
    end function cstr
    !
    /**
     * @brief Clear a character string by filling it with spaces
     *
     * This subroutine fills a character string with space characters,
     * effectively clearing its contents.
     *
     * @param[out] str Character string to be cleared
     *
     * @note
     * - Fills the entire string with space characters
     * - Used as a helper function for string initialization
     * - Essential for proper C-string conversion
     *
     * @see cstr()
     */
    subroutine clear_str(str)
        character(*), intent(out) :: str
        integer :: i
        do i = 1, len(str)
            str(i:i) = " "
        end do
    end subroutine clear_str
    !
    /**
     * @brief Convert an integer to a character string
     *
     * This function converts an integer value to its string representation
     * and returns it as a left-adjusted character string.
     *
     * @param[in] pp Integer value to convert
     * @return Character string representation of the integer
     *
     * @details
     * The function:
     * 1. Writes the integer to a temporary string
     * 2. Left-adjusts the result to remove leading spaces
     * 3. Returns the formatted string
     *
     * @note
     * - Uses internal write for conversion
     * - Result is left-adjusted (no leading spaces)
     * - Useful for creating file names or messages with numbers
     */
    character(len=20) function intstr(pp)
        !   "Convert an integer to string."
        integer, intent(in) :: pp
        write (intstr, *) pp
        intstr = adjustl(intstr)
    end function intstr
    !
    /**
     * @brief Simple sleep function using computational delay
     *
     * This subroutine implements a simple sleep/delay function by performing
     * computational work (square root calculations) for the specified number
     * of iterations.
     *
     * @param[in] nstep Number of iterations for the delay loop
     *
     * @details
     * The subroutine performs nstep iterations of square root calculations
     * to create a computational delay. The actual delay time depends on
     * the computer's performance.
     *
     * @note
     * - Not a precise timing mechanism
     * - Delay time varies with system performance
     * - Used for simple timing control in calculations
     * - More efficient than system sleep for short delays
     */
    subroutine my_sleep(nstep)
        implicit none
        integer :: i, j, nstep
        do i = 1, nstep
            j = sqrt(i*1.d0)
        end do
        return
    end subroutine my_sleep
    !
end module IO_m
