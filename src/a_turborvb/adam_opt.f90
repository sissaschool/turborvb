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

! This module was implemented by K.Nakano in Sep. 2019
! to adust tpar automatically
!

subroutine adam_opt(i_main, nweight, ndimp, first_moment, second_moment, alphab)

    implicit none

    integer ndimp, ierr, rank, i, i_main, nweight, iteration
    real*8 square_force(ndimp), alphab(ndimp), first_moment(ndimp), &
        second_moment(ndimp), adam_beta1, adam_beta2, adam_eps

#ifdef PARALLEL
    include 'mpif.h'
    call mpi_comm_rank(mpi_comm_world, rank, ierr)
#else
    rank = 0
#endif

    adam_beta1 = 0.9d0
    adam_beta2 = 0.999d0
    adam_eps = 0.00000001d0

    first_moment = adam_beta1*first_moment + (1 - adam_beta1)*alphab

    do i = 1, ndimp
        square_force(i) = alphab(i)**2
    end do

    second_moment = adam_beta2*second_moment + (1 - adam_beta2)*square_force
    iteration = i_main/nweight
    do i = 1, ndimp
        alphab(i) = (first_moment(i)*sqrt(1 - adam_beta2**iteration))/ &
                    ((1 - adam_beta1**iteration)*(sqrt(second_moment(i)) + adam_eps))
    end do

end subroutine adam_opt
