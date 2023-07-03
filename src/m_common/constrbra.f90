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

subroutine constrbra(iesfreer, n3body, jbraj, nozeroj, derjas      &
        &, econf, nw, iopt)
    implicit none
    integer iesfreer, jbraj(*), nozeroj(*), nw, i, j, n3body, iopt
    real*8 econf(nw, *), derjas(*)

    if (iesfreer .eq. 0) return

    do j = 1, iesfreer
        econf(1, j) = 0.d0
    end do

    if (iopt .eq. 0) then

        do i = 1, n3body
            j = jbraj(i)

            if (j .gt. 0) then
                econf(1, j) = econf(1, j) + derjas(nozeroj(i))
            elseif (j .lt. 0) then
                econf(1, -j) = econf(1, -j) - derjas(nozeroj(i))
            end if
        end do

    else

        do i = 1, n3body
            j = jbraj(i)
            if (j .gt. 0) then
                econf(1, j) = derjas(nozeroj(i))
            elseif (j .lt. 0) then
                econf(1, -j) = -derjas(nozeroj(i))
            end if
        end do

    end if

    return
end
subroutine constrbra_sparse(iesfreer, n3body, jbraj, nozeroj, derjas&
        &, econf, nw, iopt)
    implicit none
    integer iesfreer, jbraj(*), nozeroj(*), nw, i, j, n3body, iopt
    real*8 econf(nw, *), derjas(*)

    if (iesfreer .eq. 0) return

    do j = 1, iesfreer
        econf(1, j) = 0.d0
    end do

    if (iopt .eq. 0) then

        do i = 1, n3body
            j = jbraj(nozeroj(i))

            if (j .gt. 0) then
                econf(1, j) = econf(1, j) + derjas(i)
            elseif (j .lt. 0) then
                econf(1, -j) = econf(1, -j) - derjas(i)
            end if
        end do

    else

        do i = 1, n3body
            j = jbraj(nozeroj(i))
            if (j .gt. 0) then
                econf(1, j) = derjas(i)
            elseif (j .lt. 0) then
                econf(1, -j) = -derjas(i)
            end if
        end do

    end if

    return
end
