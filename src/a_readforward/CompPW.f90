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

subroutine CompPW(rcord, ndim, nparts, nvects, rkcomp, pwmat, ddim, vdim)
    !    compute pw matrix pwmat

    implicit none

    integer j, ndim, nparts, nvects, mnkv, ddim, vdim(*)
    real*8 rcord(ndim, *), rkcomp(ndim, *)
    real*8 pwmat(2*nvects, *)

    ! electron contribution
    do j = 1, nparts
        call cossin(rcord(1, j), pwmat(1, j), 2, nvects, ndim, rkcomp, ddim, vdim)
    end do

    return
end
