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

!> This subroutine returns typeorb for jastrow orbitals
subroutine write_type_orb(nshellj, multij, ioccj, typeorb)
    implicit none

    ! argument parameters
    integer, intent(in) :: nshellj, multij(*), ioccj(*)
    integer, intent(out) :: typeorb(*)

    ! local variables
    integer i, ii, j, ind_type

    ii = 0
    ind_type = 1

    !        write(*,*)'xxx nshell',nshellj

    do i = 1, nshellj
        !         write(*,*)'xxx',i,multij(i)
        do j = 1, multij(i)
            ii = ii + 1
            if (ioccj(ii) .ne. 0) then
                select case (multij(i))
                case (1)
                    typeorb(ind_type) = 0
                    ind_type = ind_type + 1
                case (3)
                    typeorb(ind_type) = j
                    ind_type = ind_type + 1
                case (5)
                    typeorb(ind_type) = 3 + j
                    ind_type = ind_type + 1
                case (7)
                    typeorb(ind_type) = 8 + j
                    ind_type = ind_type + 1
                case (9)
                    typeorb(ind_type) = 15 + j
                    ind_type = ind_type + 1
                case default
                    write (6, *) 'ERROR non existing orbital in Jastrow '
                end select
            end if
        end do
    end do
    !        write(*,*)'XXXX',ind_type-1
    !        do i=1,ind_type-1
    !        write(*,*)i,typeorb(i)
    !        enddo
    return
end
