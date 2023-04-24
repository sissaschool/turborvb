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

subroutine copy_eagp(small2big, ipc, nelorb_c, nelcol_c, detmat_small, eagp_pfaff, detmat_big)
    implicit none
    logical small2big
    integer ipc, nelorb_c, nelcol_c, i, j
    real*8 detmat_small(ipc*nelorb_c, nelcol_c), detmat_big(ipc*nelcol_c, nelcol_c)&
            &, eagp_pfaff(ipc*(nelcol_c - nelorb_c), nelcol_c - nelorb_c)
    if (small2big) then
        do j = 1, nelcol_c
            do i = 1, ipc*nelorb_c
                detmat_big(i, j) = detmat_small(i, j)
            end do
        end do
        do i = nelorb_c + 1, nelcol_c
            do j = 1, nelorb_c
                detmat_big(ipc*(i - 1) + 1:ipc*i, j) = -detmat_small(ipc*(j - 1) + 1:ipc*j, i)
            end do
            do j = nelorb_c + 1, nelcol_c
                detmat_big(ipc*(i - 1) + 1:ipc*i, j) = &
                        &eagp_pfaff(ipc*(i - 1 - nelorb_c) + 1:ipc*(i - nelorb_c), j - nelorb_c)
            end do
        end do
    else
        do j = 1, nelcol_c
            do i = 1, ipc*nelorb_c
                detmat_small(i, j) = detmat_big(i, j)
            end do
        end do
        do i = nelorb_c + 1, nelcol_c
            do j = nelorb_c + 1, nelcol_c
                eagp_pfaff(ipc*(i - nelorb_c - 1) + 1:ipc*(i - nelorb_c), j - nelorb_c) = &
                        &detmat_big(ipc*(i - 1) + 1:ipc*i, j)
            end do
        end do
    end if
end subroutine copy_eagp
