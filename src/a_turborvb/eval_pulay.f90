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

subroutine eval_pulay(winv, nelorb, nel, nelup, neldo, indt           &
        &, ainvup, ainvdo, winvj, nelorbj, winvjbar, orbderiv                  &
        &, winvjbarsz, iessz, add_pulay, psip, nelorbh, nelorbjh)
    use Constants, only: ip4
    implicit none
    integer nel, nelorb, nelup, indt, neldo, i, l, nelorbj, j, m, d            &
            &, add_pulay, nelorbj5, nelorbh, nelorbjh, indt5, nelorb5
    real(8) winv(nelorb, 0:indt + ip4, *), orbderiv(3, *)                  &
            &, ainvup(nelup, *), ainvdo(neldo, *), ddot, winvjbar(nelorbjh, *)        &
            &, winvjbarsz(*), winvj(nelorbj, 0:indt + ip4, *), tmp1, tmp2              &
            &, psip(nelorbjh, *)
    logical iessz

    call dscalzero(3*(nelorb + nelorbj), 0.d0, orbderiv, 1)

    indt5 = indt + ip4 + 1
    nelorbj5 = nelorbj*indt5
    nelorb5 = nelorb*indt5

    if (add_pulay .gt. 1) then
        do l = 1, nelorbh
            do d = 1, 3
                orbderiv(d, l) = orbderiv(d, l)                                     &
                        & - ddot(nelup, winv(l, indt + d, 1), nelorb5, ainvup(1, l), 1)
                orbderiv(d, l) = orbderiv(d, l)                                       &
                        & - ddot(neldo, winv(l, indt + d, 1 + nelup), nelorb5, ainvdo(1, l), 1)
            end do
        end do
    end if

    do i = 1, nelorbjh
        psip(i, 1) = 0.d0
        psip(i, 2) = 0.d0
        psip(i, 3) = 0.d0
        psip(i, 4) = 0.d0
        do j = 1, nel
            psip(i, 1) = psip(i, 1) + winvj(i, indt + 1, j)
            psip(i, 2) = psip(i, 2) + winvj(i, indt + 2, j)
            psip(i, 3) = psip(i, 3) + winvj(i, indt + 3, j)
            psip(i, 4) = psip(i, 4) + winvjbar(i, j)
        end do
    end do

    do l = 1, nelorbjh
        orbderiv(1, nelorb + l) = orbderiv(1, nelorb + l) - psip(l, 1)*psip(l, 4)  &
                & + ddot(nel, winvj(l, indt + 1, 1), nelorbj5, winvjbar(l, 1), nelorbjh)
        orbderiv(2, nelorb + l) = orbderiv(2, nelorb + l) - psip(l, 2)*psip(l, 4)  &
                & + ddot(nel, winvj(l, indt + 2, 1), nelorbj5, winvjbar(l, 1), nelorbjh)
        orbderiv(3, nelorb + l) = orbderiv(3, nelorb + l) - psip(l, 3)*psip(l, 4) &
                & + ddot(nel, winvj(l, indt + 3, 1), nelorbj5, winvjbar(l, 1), nelorbjh)
    end do

    if (iessz) then
        do i = 1, nelorbjh
            psip(i, 1) = 0.d0
            psip(i, 2) = 0.d0
            psip(i, 3) = 0.d0
            psip(i, 4) = 0.d0
            do j = 1, nelup
                psip(i, 1) = psip(i, 1) + winvj(i, indt + 1, j)
                psip(i, 2) = psip(i, 2) + winvj(i, indt + 2, j)
                psip(i, 3) = psip(i, 3) + winvj(i, indt + 3, j)
                psip(i, 4) = psip(i, 4) + winvjbarsz(i + (j - 1)*nelorbjh)
            end do
            do j = nelup + 1, nel
                psip(i, 1) = psip(i, 1) - winvj(i, indt + 1, j)
                psip(i, 2) = psip(i, 2) - winvj(i, indt + 2, j)
                psip(i, 3) = psip(i, 3) - winvj(i, indt + 3, j)
                psip(i, 4) = psip(i, 4) - winvjbarsz(i + (j - 1)*nelorbjh)
            end do
        end do
        do l = 1, nelorbjh
            !         tmp1=ddot(nelorbj,jasmatsz(l),nelorbj,winvj(1,0,j),1)
            !         tmp2=ddot(nelorbj,jasmatsz(l),nelorbj,winvj(1,0,i),1)
            orbderiv(1, nelorb + l) = orbderiv(1, nelorb + l) - psip(l, 1)*psip(l, 4) &
                    & + ddot(nel, winvj(l, indt + 1, 1), nelorbj5, winvjbarsz(l), nelorbjh)
            orbderiv(2, nelorb + l) = orbderiv(2, nelorb + l) - psip(l, 2)*psip(l, 4) &
                    & + ddot(nel, winvj(l, indt + 2, 1), nelorbj5, winvjbarsz(l), nelorbjh)
            orbderiv(3, nelorb + l) = orbderiv(3, nelorb + l) - psip(l, 3)*psip(l, 4) &
                    & + ddot(nel, winvj(l, indt + 3, 1), nelorbj5, winvjbarsz(l), nelorbjh)
        end do
    end if

    return
end
