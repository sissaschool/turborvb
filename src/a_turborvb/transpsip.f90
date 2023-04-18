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

subroutine transpsip(psip, iesupr, iesupc, iesuptrans                &
        &, multranspip, transpip, derdet_mu)
    use types
    implicit none
    integer iesupr, imax, iesupc, i, iesuptrans(*), multranspip(*), j
    type(array_int), intent(in) :: transpip(*)
    real*8 psip(*), derdet_mu(*)

    imax = max(iesupr, iesupc)
    call dcopy(iesupr, psip, 1, psip(imax + 1), 1)
    call dscalzero(iesupc, 0.d0, psip, 1)

    ! first zeta
    do i = 1, iesupr
        if (iesuptrans(i) .ne. 0) psip(iesuptrans(i)) = psip(imax + i)
    end do

    ! then mu_c (coefficients)
    do i = 1, iesupc
        if (multranspip(i) .ne. 0) then
            psip(i) = 0.d0
            do j = 1, multranspip(i)
                psip(i) = psip(i) + derdet_mu(transpip(j)%col(i))
            end do
        end if
    end do

    return
end
