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

program openmp_test

    use omp_lib

    integer :: num_threads

    call omp_set_num_threads(2)
    !$omp parallel
    !$omp master
    num_threads = omp_get_num_threads()
    !$omp end master
    !$omp end parallel

    write (*,*) "Number of OMP threads:", num_threads

    if (num_threads.ne.2) then
        stop 1
    end if

    stop 0

end program openmp_test
