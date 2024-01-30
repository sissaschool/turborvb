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

subroutine scalevect(n, cellfat, vect)
    use cell, only: map
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: cellfat(3)
    real*8, intent(inout) :: vect(3*n)

    integer i
!$omp parallel do default(shared) private(i)
    do i = 1, 3*n, 3
        vect(i) = map(vect(i), cellfat(1))
        vect(i + 1) = map(vect(i + 1), cellfat(2))
        vect(i + 2) = map(vect(i + 2), cellfat(3))
    end do
!$omp end parallel do
    return
end
