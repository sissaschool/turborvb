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
