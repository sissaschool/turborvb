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
    integer, parameter :: N = 1000
    integer :: ii, jj
    real*8 :: A(N, N), B(N), tmp

    do jj = 1, N
        B(jj) = 0
        do ii = 1, N
            A(ii, jj) = ii + jj
        end do
    end do

#ifdef _OFFLOAD
    !$omp target data map(tofrom: A, B)
    !$omp target teams distribute parallel do default(shared) reduction(+:B) collapse(2)
#else
    !$omp parallel do reduction(+:B) collapse(2)
#endif
    do jj = 1, N
        do ii = 1, N
            B(jj) = B(jj) + A(ii, jj)
        end do
    end do
#ifdef _OFFLOAD
    !$omp end target teams distribute parallel do
    !$omp end target data
#else
    !$omp end parallel do
#endif

    do jj = 1, N
        if (jj*N + (N*(N + 1))/2 .ne. B(jj)) then
            print *, "Error"
            stop 1
        end if
    end do

    print *, "Everything is OK!"
    stop 0

end program offload_data_real
