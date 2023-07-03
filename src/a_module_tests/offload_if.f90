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

program offload_if
    implicit none
    real*8 :: A(1), B(1)
    logical :: flag

    A = 5
    B = 5
    flag = .true.

    !$omp target data map(tofrom: A, B) if(flag)
    A(1) = 1
    B(1) = 1
    !$omp end target data

    if (A(1) .eq. 1) then
        print *, "Offload 1 of array A went wrong"
        stop 1
    end if

    if (B(1) .eq. 1) then
        print *, "Offload 1 of array B went wrong"
        stop 1
    end if

    flag = .false.

    !$omp target data map(tofrom: A, B) if(flag)
    A(1) = 5
    B(1) = 5
    !$omp end target data

    if (A(1) .ne. 5) then
        print *, "Offload 2 of array A went wrong"
        stop 1
    end if

    if (B(1) .ne. 5) then
        print *, "Offload 2 of array B went wrong"
        stop 1
    end if

    print *, "Everything is OK!"
    stop 0

end program offload_if
