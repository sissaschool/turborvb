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

!> @brief 3D data visualization output in XCrySDen format for tilted cells
!>
!> This subroutine generates 3D data visualization files in XCrySDen (.xsf) format
!> for quantum Monte Carlo calculations with tilted (non-orthogonal) unit cells.
!> It outputs crystal structure information and 3D data grids that can be
!> visualized using XCrySDen or similar software.
!>
!> The subroutine handles unit conversion from atomic units to Angstroms and
!> supports both periodic and non-periodic systems with proper file naming
!> conventions. Unlike plot_3d_data, this version supports general 3x3
!> transformation matrices for tilted cells.
!>
!> @param[in] ipc Number of components per grid point (1 for scalar, 3 for vector)
!> @param[in] s2r_loc Local transformation matrix for grid coordinates (3,3)
!> @param[in] s2r Global transformation matrix for crystal structure (3,3)
!> @param[in] natoms Number of atoms in the system
!> @param[in] pos Atomic positions in atomic units (3, natoms)
!> @param[in] zeta Atomic numbers or charges (natoms)
!> @param[in] periodic Flag indicating periodic boundary conditions
!> @param[in] mesh Grid dimensions (3)
!> @param[in] origin Origin of the grid in atomic units (3)
!> @param[in] datagrid 3D data grid to be plotted (ipc*mesh(1), mesh(2), mesh(3))
!> @param[in] indmol Molecular index for file naming
!> @param[in] word Identifier string for file naming
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Checks file unit availability and opens output file
!> 2. Converts all units from atomic units to Angstroms:
!>    - Transformation matrices (s2r and s2r_loc)
!>    - Origin and atomic positions
!>    - Data grid values (divided by length_unit^3)
!> 3. Writes XCrySDen format header:
!>    - CRYSTAL section with primitive vectors from s2r matrix
!>    - PRIMCOORD section with atomic positions and charges
!> 4. Writes 3D data grid section:
!>    - BEGIN_BLOCK_DATAGRID_3D header
!>    - Grid dimensions and origin
!>    - Primitive vectors for grid orientation from s2r_loc matrix
!>    - Data values in optimized format (3 values per line)
!> 5. Closes the output file
!>
!> @note Output filename format: output_<word><indmol>.xsf
!> @note The subroutine handles both scalar (ipc=1) and vector (ipc=3) data
!> @note Data is written in chunks of 3 values per line for efficiency
!> @note Used for visualizing electron density, wave functions, and other
!>       quantum mechanical properties in tilted cell geometries
!> @note Compatible with XCrySDen, VESTA, and other visualization software
!> @note The s2r matrix defines the crystal structure, while s2r_loc defines
!>       the local grid coordinate system
subroutine plot_3d_data_tilted(ipc, s2r_loc, s2r, natoms, pos, zeta, periodic, mesh, origin, datagrid, indmol, word)

    use constants, only: length_unit

    implicit none
    integer, parameter :: ofile = 55555555
    !
    integer, intent(in) :: mesh(3), ipc
    double precision :: datagrid(ipc*mesh(1), mesh(2), mesh(3))
    double precision :: s2r(3, 3), s2r_rescaled(3, 3), s2r_loc(3, 3), s2r_loc_rescaled(3, 3)
    double precision :: origin(3), origin_rescaled(3)
    double precision :: pos(3, natoms), pos_rescaled(3, natoms)
    double precision, parameter :: eps = 1.d-6
    real*8, intent(in) :: zeta(natoms)
    integer, intent(in) :: natoms
    logical, intent(in) :: periodic
    character(14), intent(in) :: word
    !
    integer :: i1, i2, indmol, idigit(6), ix, iy, iz, i3
    logical :: control
    logical :: UNITOK, UNITOP

    inquire (unit=ofile, exist=UNITOK, opened=UNITOP)

    if (UNITOK .and. .not. UNITOP) then

        call convertdec(indmol, idigit)
        open (unit=ofile, file="output_"//trim(word)//char(idigit(1))//char(idigit(2))// &
              char(idigit(3))//char(idigit(4))//char(idigit(5))//char(idigit(6))//".xsf", form="formatted", status="unknown")

    else

        write (6, *) 'error from open unit in plot_3d_data, unit already opened'
        return

    end if

    ! rescale units from au to Angstrom
    s2r_rescaled = s2r*length_unit
    s2r_loc_rescaled = s2r_loc*length_unit
    origin_rescaled = origin*length_unit
    pos_rescaled = pos*length_unit
    datagrid = datagrid/length_unit**3

    !

    write (ofile, *) 'CRYSTAL '
    write (ofile, *) 'PRIMVEC '
    write (ofile, *) s2r_rescaled(:, 1)
    write (ofile, *) s2r_rescaled(:, 2)
    write (ofile, *) s2r_rescaled(:, 3)
    write (ofile, *) 'PRIMCOORD '
    write (ofile, *) natoms, 1
    do i1 = 1, natoms
        write (ofile, '(i3,6(F16.8,4X))') nint(zeta(i1)), pos_rescaled(:, i1)
    end do

    if (any(mesh(:) .ne. 0)) then
        write (ofile, *)
        write (ofile, *) 'BEGIN_BLOCK_DATAGRID_3D'
        write (ofile, *) '  my_datagrid'
        write (ofile, *) '  BEGIN_DATAGRID_3D_datagrid'
        write (ofile, '(3i5)') mesh(:)
        write (ofile, '(3f14.8)') origin_rescaled(:)
        write (ofile, '(3f14.8)') s2r_loc_rescaled(:, 1)
        write (ofile, '(3f14.8)') s2r_loc_rescaled(:, 2)
        write (ofile, '(3f14.8)') s2r_loc_rescaled(:, 3)
        !
        !  do ix=1,mesh(1)
        !    do iy=1,mesh(2)
        !      write(ofile,*) (datagrid(ix,iy,iz),iz=1,mesh(3))
        !    enddo
        !  enddo
        !
        !write(ofile,*) (((datagrid(ix,iy,iz),iz=1,mesh(3)),iy=1,mesh(2)),ix=1,mesh(1))
        do i1 = 1, mesh(3)
            do i2 = 1, mesh(2)
                !       write(ofile,'(1000(e14.6,4x))') (datagrid(i3,i2,i1),i3=1,mesh(1))
                do i3 = 1, mesh(1), 3
                    if (i3 + 2 .le. mesh(1)) then
                        write (ofile, *) datagrid(ipc*(i3 - 1) + 1, i2, i1), datagrid(ipc*i3 + 1, i2, i1)&
                                &, datagrid(ipc*(i3 + 1) + 1, i2, i1)
                    elseif (i3 + 1 .le. mesh(1)) then
                        write (ofile, *) datagrid(ipc*(i3 - 1) + 1, i2, i1), datagrid(ipc*i3 + 1, i2, i1)
                    else
                        write (ofile, *) datagrid(ipc*(i3 - 1) + 1, i2, i1)
                    end if
                end do
            end do
        end do
        !
        write (ofile, *) 'END_DATAGRID_3D'
        write (ofile, *) 'END_BLOCK_DATAGRID_3D'

    end if

    close (ofile)

end subroutine plot_3d_data_tilted
