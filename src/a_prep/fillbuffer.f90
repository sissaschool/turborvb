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

subroutine fillbuff_complex(nelorb, indt, istart, indtms, tcost, vpots, winv_c, buffer)

    use constants, only: zhalf, zone

    implicit none

    ! input variables
    integer, intent(in) :: nelorb, indt, istart, indtms
    real(8), intent(in) :: tcost(*), vpots
    complex(8), intent(in) :: winv_c(nelorb, 0:indt + 4)
    complex(8), intent(out) :: buffer(nelorb)

    ! local variables
    integer :: dims, i

    dims = indtms - istart + 1
    buffer(1:nelorb) = zhalf*dcmplx(vpots)*winv_c(1:nelorb, 0) ! unit Hartree

    if (dims .gt. 0) then
        ! adding non local part if any
        call dgemv('N', 2*nelorb, dims, -0.5d0, winv_c(1, istart), 2*nelorb, tcost(istart)&
                &, 1, 1.d0, buffer, 1)
    end if
    !      SSTEST
    !      buffer(1:nelorb)=0.d0 ! disregard potential for tests
    !      adding kinetic energy unit Hartree
    buffer(1:nelorb) = buffer(1:nelorb) - 0.5d0*winv_c(1:nelorb, indt + 4)

    return
end subroutine fillbuff_complex

subroutine fillbuff(nelorb, indt, istart, indtms, tcost, vpots, winv, buffer)
    implicit none
    integer nelorb, indt, istart, indtms, dims, i
    real*8 tcost(*), vpots, winv(nelorb, 0:indt + 4), buffer(*)
    dims = indtms - istart + 1
    buffer(1:nelorb) = 0.5d0*vpots*winv(1:nelorb, 0) ! unit Hartree
    if (dims .gt. 0) then
        !      adding non local part if any
        call dgemv('N', nelorb, dims, -0.5d0, winv(1, istart), nelorb, tcost(istart)&
                &, 1, 1.d0, buffer, 1)
        !      write(6,*) ' buffer inside fillbuf ',sum(buffer(1:nelorb))
    end if
    !      SSTEST
    !      buffer(1:nelorb)=0.d0 ! disregard potential for tests
    !      adding kinetic energy unit Hartree
    buffer(1:nelorb) = buffer(1:nelorb) - 0.5d0*winv(1:nelorb, indt + 4)

end subroutine fillbuff

