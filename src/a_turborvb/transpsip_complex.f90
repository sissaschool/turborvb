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

subroutine transpsip_complex(psip, iesupr, iesupc, iesuptrans&
        &, multranspip, transpip, derdet_mu)
    use types
    use Constants, only: ipc
    implicit none
    integer iesupr, iesupc, i, imax, iesuptrans(*), multranspip(*), j
    type(array_int), intent(in) :: transpip(*)
    real*8 psip(*), derdet_mu(*)

    if (ipc .eq. 2) then
        imax = max(2*iesupc, iesupr)
        call dcopy(iesupr, psip, 1, psip(imax + 1), 1)
        call dscalzero(2*iesupc, 0.d0, psip, 1)

        ! first zeta
        do i = 1, iesupr
            if (iesuptrans(i) .ne. 0) then
                psip(2*iesuptrans(i) - 1) = psip(imax + i)
                psip(2*iesuptrans(i)) = 0.d0
            end if
        end do

        ! then mu_c (coefficients)
        do i = 1, iesupc
            if (multranspip(i) .ne. 0) then
                psip(2*i - 1:2*i) = 0.d0
                do j = 1, multranspip(i)
                    psip(2*i - 1) = psip(2*i - 1) + derdet_mu(2*transpip(j)%col(i) - 1)
                    psip(2*i) = psip(2*i) + derdet_mu(2*transpip(j)%col(i))
                end do
            end if
        end do

    else
        imax = max(iesupc, iesupr)
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
    end if

    return
end
