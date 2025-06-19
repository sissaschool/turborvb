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

!> @file evaluate_invariant.f90
!> @brief Evaluate orbital invariants for Jastrow or symmetry-adapted wave functions.
!>
!> This file provides a subroutine to compute invariants (such as s-s, s-p, p-p)
!> between two orbitals, which are used in the construction of Jastrow factors or
!> symmetry-adapted wave functions in quantum Monte Carlo calculations. The invariants
!> are determined based on the types of the two orbitals and their relative position vector.
!>
!> @author TurboRVB group
!> @date 2022
!>
!> @section Usage
!> Typically used in the evaluation of Jastrow factors or when constructing
!> symmetry-adapted basis functions for correlated wave functions.
!>
!> @section Examples
!> @code{.f90}
!>   call evaluate_invariant(vec_r, ix, iy, typeorb, jas_invariant, orbps)
!> @endcode

!-------------------------------------------------------------------------------
!> @brief Evaluate invariants between two orbitals for Jastrow or symmetry purposes.
!>
!> This subroutine computes a set of invariants (e.g., s-s, s-p, p-p) between two orbitals
!> based on their types and the relative position vector. These invariants are used in the
!> construction of Jastrow factors or symmetry-adapted wave functions.
!>
!> @param[in] vec_r Relative position vector between the two orbitals (length >= 3)
!> @param[in] ix Index of the first orbital
!> @param[in] iy Index of the second orbital
!> @param[in] typeorb Array specifying the type of each orbital (0=s, 1-3=p_x/p_y/p_z)
!> @param[out] jas_invariant Array of computed invariants (length 4)
!> @param[out] orbps Logical flag, true if the order of orbitals is s-p (for symmetry)
!>
!> The invariants are:
!>   - jas_invariant(1): s-s overlap (always 1 if both are s)
!>   - jas_invariant(2): s-p or p-s overlap (component of vec_r)
!>   - jas_invariant(3): p-p overlap (product of components)
!>   - jas_invariant(4): 1 if both p orbitals are the same, 0 otherwise
!>
!> @note The routine is designed for use in QMC codes where orbital symmetry and
!>       Jastrow factor construction require such invariants.
!>
!> @section Examples
!> @code{.f90}
!>   real(8) :: vec_r(3), jas_invariant(4)
!>   integer :: typeorb(2)
!>   logical :: orbps
!>   typeorb = [0, 1]  ! s and p_x
!>   call evaluate_invariant(vec_r, 1, 2, typeorb, jas_invariant, orbps)
!> @endcode
subroutine evaluate_invariant(vec_r, ix, iy, typeorb, jas_invariant, orbps)
    implicit none

    real(8) jas_invariant(4)
    integer i, ix, iy, typeorb(*)
    real(8) vec_r(*)
    logical orbps

    !       vec_phi1=0.d0
    !       vec_phi2=0.d0

    jas_invariant = 0.d0
    orbps = .false.

    !############################################################
    !       s s
    if (typeorb(ix) .eq. 0 .and. typeorb(iy) .eq. 0) then
        jas_invariant(1) = 1.d0
    end if

    !#############################################################
    !       s p case
    if ((typeorb(ix) .gt. 0 .and. typeorb(ix) .le. 3) .and. typeorb(iy) .eq. 0) &
            & then

        !        vec_phi1(typeorb(ix))=1
        !        phi1_dot_r=0.d0

        !        do i=1,3
        !         phi1_dot_r=phi1_dot_r+vec_phi1(i)*vec_r(i)
        !        enddo

        jas_invariant(2) = vec_r(typeorb(ix))

        !        if(abs(jas_invariant(2)).gt.2)&
        !       & write(6,*) ' ERROR in jas invariant '

        !        jas_invariant(2)=jas_invariant(2)+10.d0  ! I pass in this way the information
        !            that typeorb(ix)>typeorb(iy)
        orbps = .true.

    elseif                                                           &
            & ((typeorb(iy) .gt. 0 .and. typeorb(iy) .le. 3) .and. typeorb(ix) .eq. 0)   &
            &then

        !        vec_phi1(typeorb(iy))=1
        !        phi1_dot_r=0.d0

        !        do i=1,3
        !         phi1_dot_r=phi1_dot_r+vec_phi1(i)*vec_r(i)
        !        enddo
        jas_invariant(2) = vec_r(typeorb(iy))

    end if

    !###################################################################
    !       p p case
    if ((typeorb(ix) .gt. 0 .and. typeorb(ix) .le. 3) .and.                 &
            &     (typeorb(iy) .gt. 0 .and. typeorb(iy) .le. 3)) then

        !         vec_phi1(typeorb(ix))=1
        !         vec_phi2(typeorb(iy))=1

        !         phi1_dot_r=0.d0
        !         phi2_dot_r=0.d0
        !
        !         do i=1,3
        !          phi1_dot_r =phi1_dot_r+ vec_phi1(i)*vec_r(i)
        !          phi2_dot_r =phi2_dot_r+ vec_phi2(i)*vec_r(i)
        !         enddo

        jas_invariant(3) = vec_r(typeorb(ix))*vec_r(typeorb(iy))

        if (typeorb(ix) .eq. typeorb(iy)) jas_invariant(4) = 1.d0

    end if

    return
end
