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

module grid_module

    real*8, allocatable :: out_grid(:, :)
    logical, allocatable :: active_points(:)
    real*8 :: grid_start(3), center_ion(3)
    real*8 :: da(3)
    integer :: grid_points
    logical :: ext_grid
    integer :: center

contains

    subroutine compute_grid(ngrid_l, ell, vell)

        use allio, only: rank, nion, rion
        implicit none

#ifdef PARALLEL
        include 'mpif.h'
        integer ierr
#endif

        integer :: ngrid_l(3)
        double precision :: ell(3), vell(3)
        logical :: correct
        double precision :: e_k

        integer :: i, j, k, ii

        allocate (out_grid(3, grid_points))
        allocate (active_points(1:grid_points))

        correct = .true.

        if (ext_grid) then

            if (rank .eq. 0) then

                write (6, *) ".....reading grid from turbo_grid.dat file"

                inquire (file='turbo_grid.dat', exist=correct)

                if (correct) then
                    open (unit=28, file='turbo_grid.dat', form='formatted', status='unknown')
                else
                    goto 100
                end if

                do i = 1, grid_points
                    read (28, *) (out_grid(j, i), j=1, 3)
                end do

                write (6, *) "turbo_grid.dat read correctly"

                close (28)
            end if

            active_points(:) = .true.

        else
            ii = 0

            vell(:) = ell(:)/(ngrid_l(:))

            if (rank .eq. 0) then
                write (6, *) "Grid vectors = ", vell
                write (6, *) "Grid points per dimension =", ngrid_l
            end if

            if (da(1) .eq. 0 .or. da(2) .eq. 0 .or. da(3) .eq. 0) da(:) = vell(:)

            do k = 1, ngrid_l(3)
                do j = 1, ngrid_l(2)
                    do i = 1, ngrid_l(1)
                        ii = ii + 1
                        out_grid(1, ii) = (i - 0.5)*vell(1)
                        out_grid(2, ii) = (j - 0.5)*vell(2)
                        out_grid(3, ii) = (k - 0.5)*vell(3)
                    end do
                end do
            end do

            ! Center=1 the center of the grid is the center of the axes
            ! Center=2 the center of the grid is the barycenter of the atoms
            ! Center=3 the origin of the grid is a point in space given externally
            !          In this last case do nothing.
            if (center .eq. 1) then
                grid_start(:) = -ell(:)/2.0
            elseif (center .eq. 2) then
                do i = 1, nion
                    center_ion(:) = center_ion(:) + rion(:, i)
                end do
                center_ion = center_ion/nion
                grid_start(:) = -ell(:)/2.0 + center_ion(:)
            end if

            if (vell(1) .eq. 0) then
                e_k = (ell(2) + grid_start(2))**2/2.0
            elseif (vell(2) .eq. 0 .or. vell(3) .eq. 0) then
                e_k = (ell(1) + grid_start(1))**2/2.0
            end if

            if (rank .eq. 0) write (6, *) e_k

            !      allocate(active_points(1:grid_points))
            if (e_k .gt. 0) then
                active_points(:) = .false.
            else
                active_points(:) = .true.
            end if

            do i = 1, grid_points
                out_grid(:, i) = out_grid(:, i) + grid_start(:)
                !            if(rank.eq.0) write(6,*) i,e_k, out_grid(1,i),out_grid(2,i),out_grid(3,i),&
                !                                      out_grid(1,i)**2+out_grid(2,i)**2+out_grid(3,i)**2
                if (e_k .gt. 0) then
                    if (vell(1) .eq. 0) then
                        if ((out_grid(2, i)**2 + out_grid(3, i)**2) .le. 2*e_k) then
                            out_grid(1, i) = sqrt(2*e_k - out_grid(2, i)**2 - out_grid(3, i)**2)
                            active_points(i) = .true.
                        end if
                    elseif (vell(2) .eq. 0) then
                        if ((out_grid(1, i)**2 + out_grid(3, i)**2) .le. 2*e_k) then
                            out_grid(2, i) = sqrt(2*e_k - out_grid(1, i)**2 - out_grid(3, i)**2)
                            active_points(i) = .true.
                        end if
                    elseif (vell(3) .eq. 0) then
                        if ((out_grid(1, i)**2 + out_grid(2, i)**2) .le. 2*e_k) then
                            out_grid(3, i) = sqrt(2*e_k - out_grid(1, i)**2 - out_grid(2, i)**2)
                            active_points(i) = .true.
                        end if
                    end if
                end if
                !        if(rank.eq.0) write(6,*) active_points(i)
            end do

        end if

#ifdef PARALLEL
        call mpi_bcast(correct, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

100     if (.not. correct) then
            if (rank .eq. 0) write (6, *) "ERROR! Missing turbo_grid.dat file"
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        elseif (ext_grid) then
#ifdef PARALLEL
            call bcast_real(out_grid, 3*grid_points, 0, MPI_COMM_WORLD)
#endif
        end if

#ifdef PARALLEL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

    end subroutine compute_grid

end module grid_module
