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

program offload_data_complex
    implicit none
    complex*16 :: A(1), B(1)

    A = 2
    B = 2

    !$omp target data map(tofrom: A, B)
    !$omp target
    A(1) = 1
    B(1) = 1
    !$omp end target
    !$omp end target data

    if (A(1) .ne. 1) then
        print *, "Offload of array A went wrong"
        stop 1
    end if

    if (B(1) .ne. 1) then
        print *, "Offload of array A went wrong"
        stop 1
    end if

    print *, "Everything is OK!"
    stop 0

end program offload_data_complex

