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

subroutine bconstrbra(iesfreer, n3body, jbraj, nozeroj, derjas     &
        &, nelorbj, ddw)
    implicit none
    integer iesfreer, jbraj(*), nozeroj(*), i, j, n3body, nelorbj
    real*8 ddw(*), derjas(*)
    if (iesfreer .eq. 0) return
    do i = 1, n3body
        j = jbraj(i)
        if (j .gt. 0) then
            call upsim(derjas, nelorbj, nozeroj(i), ddw(j), .true., 1)
        elseif (j .lt. 0) then
            call upsim(derjas, nelorbj, nozeroj(i), -ddw(-j), .true., 1)
        end if
    end do
    return
end
subroutine bconstrbra_sparse(iesfreer, n3body, jbraj, nozeroj, derjas&
        &, nelorbj, ddw)
    implicit none
    integer iesfreer, jbraj(*), nozeroj(*), i, j, n3body, nelorbj
    real*8 ddw(*), derjas(*)
    if (iesfreer .eq. 0) return
    do i = 1, n3body
        j = jbraj(nozeroj(i))
        if (j .gt. 0) then
            derjas(i) = ddw(j)
!           call  upsim(derjas, nelorbj, nozeroj(i), ddw(j), .true., 1)
        elseif (j .lt. 0) then
            derjas(i) = -ddw(-j)
!           call  upsim(derjas, nelorbj, nozeroj(i), -ddw(-j), .true., 1)
        end if
    end do
    return
end
