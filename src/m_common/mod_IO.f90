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

module IO_m
    implicit none
    integer, parameter :: lchlen = 1000

    !
contains
    !
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

    subroutine cd_dir(dirname)
        implicit none
        character(*) :: dirname
        if (len_trim(dirname) == 0) return
        call ichdir(cstr(trim(dirname)))
    end subroutine
    !
    subroutine cp_file(file_, dest_, ierr_)
        implicit none
        character(*) :: file_, dest_
        integer :: ierr_
        call isystem(cstr("cp "//file_//" "//dest_//" >& /dev/null"), ierr_)
    end subroutine
    !
    subroutine rm_file(filename)
        implicit none
        character(*) :: filename
        if (len_trim(filename) == 0) return
        call iremove(cstr(trim(filename)))
    end subroutine
    !
    subroutine rename_file(filename_old, filename_new)
        implicit none
        character(*) :: filename_old, filename_new
        if (len_trim(filename_old) == 0) return
        call irename(cstr(trim(filename_old)), cstr(trim(filename_new)))
    end subroutine
    !
    subroutine get_dir(path)
        implicit none
        integer :: ln
        character(*) :: path
        call igetcwd(path, ln)
        path = path(1:ln)
    end subroutine get_dir
    !
    logical function file_is_open(filename)
        character(*) :: filename
        file_is_open = .false.
        if (len_trim(filename) == 0) return
        inquire (file=filename, opened=file_is_open)
        !
    end function
    !
    character(lchlen) function cstr(si) result(so)
        character(*), intent(IN) :: si
        integer :: i
        i = len(trim(si))
        call clear_str(so)
        so(1:i) = si(1:i)
        so(i + 1:i + 1) = achar(0)
    end function cstr
    !
    subroutine clear_str(str)
        character(*), intent(out) :: str
        integer :: i
        do i = 1, len(str)
            str(i:i) = " "
        end do
    end subroutine clear_str
    !
    character(len=20) function intstr(pp)
        !   "Convert an integer to string."
        integer, intent(in) :: pp
        write (intstr, *) pp
        intstr = adjustl(intstr)
    end function intstr
    !
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
