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

subroutine CompRhok(rhoknew, rhoknew2, nparts, nvects, pwmat, nppss)
    !    compute rhok
    implicit none

    integer j, k, nparts, nvects, nvects2, nppss(*)
    real*8 rhoknew(*), rhoknew2(2*nvects, *)
    real*8 pwmat(2*nvects, *)

    nvects2 = 2*nvects

    do k = 1, nvects2
        rhoknew(k) = 0.d0
        rhoknew2(k, 1) = 0.d0
        rhoknew2(k, 2) = 0.d0
    end do

    ! electron contribution
    do j = 1, nparts
        do k = 1, nvects2
            rhoknew(k) = rhoknew(k) + pwmat(k, j)
            if (j .le. nppss(1)) then
                rhoknew2(k, 1) = rhoknew2(k, 1) + pwmat(k, j)
            else
                rhoknew2(k, 2) = rhoknew2(k, 2) + pwmat(k, j)
            end if
        end do
    end do

    return
end
