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

!> @brief Find a specific section in a file
!> @details This subroutine searches for a specific section name in a file
!>          by reading through the file line by line. It rewinds the file
!>          to the beginning and searches for an exact match with the
!>          provided section name. If the section is found, the file pointer
!>          is positioned at the beginning of that section. If not found,
!>          an error is raised.
!> @param[in] funit File unit number to search in
!> @param[in] section_name Name of the section to find
subroutine findsection(funit, section_name)
    implicit none
    logical done
    integer, intent(in) :: funit
    character(len=*), intent(in) :: section_name
    character(20) :: line
    character(50) :: empty_string

    done = .false.

    rewind (funit)
    do while (.not. done)
        read (funit, '(a)', end=100) line
        !
        if (trim(line) == trim(section_name)) done = .true.
        !
    end do

    return

100 write (empty_string, *) " Section ", section_name, " not found "
    call errore("findsection", empty_string, 1)
end subroutine findsection

