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

!> @brief Computes cosine and sine values for plane wave evaluation
!> @details This subroutine calculates nk pairs of cosine and sine values
!>          for a particle located at position x. It evaluates plane wave
!>          contributions exp(i*k*r) = cos(k*r) + i*sin(k*r) for correlation
!>          function calculations. Supports 1D, 2D, and 3D systems.
!> @param[in] x Particle position coordinates
!> @param[out] sp Array containing cosine and sine values (n2, nk)
!> @param[in] n2 Dimension of sp array (should be 2 for cos/sin pairs)
!> @param[in] nk Number of wave vectors
!> @param[in] ndim Dimension of coordinate space
!> @param[in] rkcomp Wave vector components (ndim, nk)
!> @param[in] ddim Physical dimension of the system (1, 2, or 3)
!> @param[in] vdim Vector dimensions for coordinate mapping
subroutine cossin(x, sp, n2, nk, ndim, rkcomp, ddim, vdim)

    implicit none
    integer n2, nk, ndim, i, ddim, vdim(*)
    real*8 x(*), sp(n2, *), rkcomp(ndim, *), dum
    ! calculates nk sin and cos for a particle located at x
    ! and puts them in sp.
    if (ddim .eq. 3) then
        do i = 1, nk
            dum = rkcomp(1, i)*x(vdim(1)) + rkcomp(2, i)*x(vdim(2)) + rkcomp(3, i)*x(vdim(3))
            sp(1, i) = cos(dum)
            sp(2, i) = sin(dum)
        end do

    else if (ddim .eq. 2) then
        do i = 1, nk
            dum = rkcomp(1, i)*x(vdim(1)) + rkcomp(2, i)*x(vdim(2))
            sp(1, i) = cos(dum)
            sp(2, i) = sin(dum)
        end do

    else if (ddim .eq. 1) then
        do i = 1, nk
            dum = rkcomp(1, i)*x(vdim(1))
            sp(1, i) = cos(dum)
            sp(2, i) = sin(dum)
        end do

    end if
end
