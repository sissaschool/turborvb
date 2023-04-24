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

subroutine print_eigenvalues(rank, eig, maxdimeig, irankdet)
    implicit none
    integer, intent(in) :: rank, maxdimeig
    integer i, irankdet
    real*8 eig(*)
    irankdet = 0
    do i = 1, maxdimeig
        if (abs(eig(i)) .gt. 1.d-11) then
            irankdet = irankdet + 1
        end if
        if (rank .eq. 0) write (6, *) i, eig(i)
    end do
    return
end subroutine print_eigenvalues
