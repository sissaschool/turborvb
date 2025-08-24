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

!> @brief Assign orbital types for Jastrow orbitals based on shell structure
!>
!> This subroutine assigns orbital type indices to Jastrow orbitals based on
!> their shell multiplicity and occupation. The orbital types follow a specific
!> numbering scheme for different angular momentum shells (s, p, d, f, g).
!>
!> Parameters
!> ----------
!> NSHELLJ : integer, in
!>     Number of shells in the Jastrow basis.
!> MULTIJ : integer array, in
!>     Multiplicity of each shell (1=s, 3=p, 5=d, 7=f, 9=g).
!> IOCCJ : integer array, in
!>     Occupation numbers for each orbital (0=unoccupied, 1=occupied).
!> TYPEORB : integer array, out
!>     Assigned orbital type indices for occupied orbitals.
!>
!> Notes
!> -----
!> - Only occupied orbitals (IOCCJ ≠ 0) are assigned type indices.
!> - The orbital type numbering follows the scheme:
!>   * s orbitals (multiplicity 1): type = 0
!>   * p orbitals (multiplicity 3): type = 1, 2, 3
!>   * d orbitals (multiplicity 5): type = 4, 5, 6, 7, 8
!>   * f orbitals (multiplicity 7): type = 9, 10, 11, 12, 13, 14, 15
!>   * g orbitals (multiplicity 9): type = 16, 17, 18, 19, 20, 21, 22, 23, 24
!> - The subroutine processes shells sequentially and assigns types incrementally.
!> - Used in Jastrow factor construction for quantum Monte Carlo.
!>
!> Algorithm
!> ---------
!> For each shell i with multiplicity multij(i):
!>   For each orbital j in the shell:
!>     If ioccj(ii) ≠ 0 (occupied):
!>       Assign typeorb based on multiplicity and orbital index
!>
!> Example
!> -------
!> Used in Jastrow factor parameterization and orbital basis construction.
subroutine write_type_orb(nshellj, multij, ioccj, typeorb)
    implicit none

    ! argument parameters
    integer, intent(in) :: nshellj, multij(*), ioccj(*)
    integer, intent(out) :: typeorb(*)

    ! local variables
    integer i, ii, j, ind_type

!> @brief Initialize counters
    ii = 0
    ind_type = 1

!> @brief Process each shell and assign orbital types
    do i = 1, nshellj
        !         write(*,*)'xxx',i,multij(i)
        do j = 1, multij(i)
            ii = ii + 1
            if (ioccj(ii) .ne. 0) then
!> @brief Assign orbital type based on shell multiplicity
                select case (multij(i))
                case (1)
!> @brief s orbital: type = 0
                    typeorb(ind_type) = 0
                    ind_type = ind_type + 1
                case (3)
!> @brief p orbitals: type = 1, 2, 3
                    typeorb(ind_type) = j
                    ind_type = ind_type + 1
                case (5)
!> @brief d orbitals: type = 4, 5, 6, 7, 8
                    typeorb(ind_type) = 3 + j
                    ind_type = ind_type + 1
                case (7)
!> @brief f orbitals: type = 9, 10, 11, 12, 13, 14, 15
                    typeorb(ind_type) = 8 + j
                    ind_type = ind_type + 1
                case (9)
!> @brief g orbitals: type = 16, 17, 18, 19, 20, 21, 22, 23, 24
                    typeorb(ind_type) = 15 + j
                    ind_type = ind_type + 1
                case default
!> @brief Error for unsupported orbital multiplicity
                    write (6, *) 'ERROR non existing orbital in Jastrow '
                end select
            end if
        end do
    end do
    !        write(*,*)'XXXX',ind_type-1
    !        do i=1,ind_type-1
    !        write(*,*)i,typeorb(i)
    !        enddo
    return
end
