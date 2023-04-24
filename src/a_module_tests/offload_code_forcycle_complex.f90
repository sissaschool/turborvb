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

program offload_data_real
    implicit none
    integer, parameter :: N = 5000
    integer :: ii
    complex*16 :: A(N), B(N)

    do ii = 1, N
        A(ii) = ii + 1
        B(ii) = 0
    end do

    !$omp target data map(tofrom: A, B)
    !$omp target teams distribute parallel do
    do ii = 1, N
        B(ii) = A(ii) - 1
    end do
    !$omp end target teams distribute parallel do
    !$omp end target data

    B = B - A + 1

    do ii = 1, N
        if (abs(B(ii)) .gt. 1.0e-32) then
            print *, "Error"
            stop 1
        end if
    end do

    print *, "Everything is OK!"
    stop 0

end program offload_data_real

