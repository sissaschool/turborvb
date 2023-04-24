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

subroutine printden(rion_ref, rion_shift, nxl, dent)

    use allio

    integer :: mesh(3), ind, i, j, k, proc
    real(8) :: cell_dens(3), mesh_origin(3), rion_ref(*), rion_shift(*), dent(*)
    real(8), allocatable, dimension(:, :, :) :: datagrid
    real(8), allocatable, dimension(:, :) :: rion_plot

#ifdef PARALLEL
    real(8), allocatable, dimension(:, :) :: buffer_grid
    include 'mpif.h'
#endif

    allocate (rion_plot(3, nion))
    rion_plot = 0.d0

    if (.not. iespbc) then

        ! for non period systems translate both mesh_origin and nuclear geometry (rion) by L(i)/2
        ! in order to place the fragment at the center of the cell drawn by Xcrysden
        ! note: rion_ref(i) = (n(i)-1)/n(i) * L(i)/2 + mesh_offset(i) (the mesh_offset is rion_shift)
        rion_plot(1, :) = rion(1, :) + rion_ref(1) - rion_shift(1)
        rion_plot(2, :) = rion(2, :) + rion_ref(2) - rion_shift(2)
        rion_plot(3, :) = rion(3, :) + rion_ref(3) - rion_shift(3)

        cell_dens(1) = (nx - 1)*ax
        cell_dens(2) = (ny - 1)*ay
        cell_dens(3) = (nz - 1)*az

        mesh(1) = nx
        mesh(2) = ny
        mesh(3) = nz

        mesh_origin(1) = rion_ref(1)
        mesh_origin(2) = rion_ref(2)
        mesh_origin(3) = rion_ref(3)

    else

        rion_plot(1, :) = rion(1, :)
        rion_plot(2, :) = rion(2, :)
        rion_plot(3, :) = rion(3, :)

        cell_dens(1) = nx*ax
        cell_dens(2) = ny*ay
        cell_dens(3) = nz*az

        mesh(1) = nx + 1
        mesh(2) = ny + 1
        mesh(3) = nz + 1

        mesh_origin(1) = rion_shift(1)
        mesh_origin(2) = rion_shift(2)
        mesh_origin(3) = rion_shift(3)

    end if

    if (rank .eq. 0) then
        allocate (datagrid(mesh(1), mesh(2), mesh(3)))
        datagrid = 0.d0
    end if

#ifdef PARALLEL
    if (rank .eq. 0) then
        allocate (buffer_grid(nxl*ny*nz, nproc))
        buffer_grid = 0.d0
    end if

    call mpi_gather(dent, nxl*ny*nz, MPI_DOUBLE_PRECISION, buffer_grid&
     &, nxl*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (rank .eq. 0) then
        do proc = 0, nproc - 1
            ind = 0
            do k = 1, nz
                do j = 1, ny
                    do i = proc + 1, nx, nproc
                        ind = ind + 1
                        datagrid(i, j, k) = buffer_grid(ind, proc + 1)
                    end do
                end do
            end do
        end do
        deallocate (buffer_grid)
    end if

#else

    ind = 0
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                ind = ind + 1
                datagrid(i, j, k) = dent(ind)
            end do
        end do
    end do

#endif

    if (rank .eq. 0) then

        if (iespbc) then
            ! periodic condition on the grid
            do k = 1, nz
                do j = 1, ny
                    datagrid(mesh(1), j, k) = datagrid(1, j, k)
                end do
            end do

            do k = 1, nz
                do j = 1, mesh(1)
                    datagrid(j, mesh(2), k) = datagrid(j, 1, k)
                end do
            end do

            do k = 1, mesh(2)
                do j = 1, mesh(1)
                    datagrid(j, k, mesh(3)) = datagrid(j, k, 1)
                end do
            end do

        end if

        ! nion, rion, atom_number, iespbc read from read_fort10 called in the main driver

        write (6, *) 'writing dft density on xcrysden file'
        call plot_3d_data(1, cell_dens, cell_dens, nion, rion_plot &
                          , atom_number, iespbc, mesh, mesh_origin, datagrid, 1, 'dft_density')
        deallocate (datagrid)
    end if

    deallocate (rion_plot)

end subroutine printden
