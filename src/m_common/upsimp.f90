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

subroutine upsimp(amat, ndim, ind, value, symmagp, ipf)
    use allio, only: nelorb_at, pfaffup, kiontot
    implicit none
    integer ndim, ndimh, ind, i, j, ipf
    real*8 amat(ndim, *), value
    logical symmagp
    if (ipf .eq. 1) then
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = amat(i, j) + value
        if (j .ne. i .and. j .le. ndim .and. symmagp) amat(j, i) = amat(j, i) + value
    else
        ndimh = nelorb_at/2
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = amat(i, j) + value
        if (j .le. ndim) amat(j, i) = amat(j, i) - value
        if (symmagp .and. kiontot(i) .ne. 0 .and. kiontot(j) .ne. 0) then
            if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                amat(i + ndimh, j + ndimh) = amat(i + ndimh, j + ndimh) + value
                amat(j + ndimh, i + ndimh) = amat(j + ndimh, i + ndimh) - value
            elseif (j .gt. ndimh .and. i .le. ndimh) then
                amat(j - ndimh, i + ndimh) = amat(j - ndimh, i + ndimh) + value
                amat(i + ndimh, j - ndimh) = amat(i + ndimh, j - ndimh) - value
                !          elseif(j.le.ndimh.and.i.gt.ndimh) then
                !             amat(j+ndimh, i-ndimh)=amat(j+ndimh, i-ndimh)+value
                !             amat(i-ndimh, j+ndimh)=amat(i-ndimh, j+ndimh)-value
            end if
        end if
    end if
    return
end
subroutine upsimp_complex(amat, ndim, ind, value, symmagp, ipf)
    use allio, only: yes_hermite, nelorb_at, pfaffup, kiontot
    implicit none
    integer ndim, ndimh, ind, ipf, i, j
    complex*16 amat(ndim, *), value(*)
    logical symmagp
    if (ipf .eq. 1) then
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        if (i .eq. j .and. yes_hermite .and. symmagp .and. j .le. nelorb_at) then
            amat(i, i) = amat(i, i) + real(value(1))
        else
            amat(i, j) = amat(i, j) + value(1)
        end if
        if (j .ne. i .and. j .le. nelorb_at .and. i .le. nelorb_at .and. symmagp) then
            if (yes_hermite) then
                amat(j, i) = amat(j, i) + dconjg(value(1))
            else
                amat(j, i) = amat(j, i) + value(1)
            end if
        end if
    else
        ndimh = nelorb_at/2
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = amat(i, j) + value(1)
        !if(j.le.2*nelorb_at.and.i.le.2*nelorb_at)
        if (j .le. ndim) amat(j, i) = amat(j, i) - value(1) !  it has to be replicated even if it is unpaired

        if (symmagp .and. kiontot(j) .ne. 0 .and. kiontot(i) .ne. 0) then
            if (yes_hermite) then
                if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                    amat(i + ndimh, j + ndimh) = amat(i + ndimh, j + ndimh) + conjg(value(1))
                    amat(j + ndimh, i + ndimh) = amat(j + ndimh, i + ndimh) - conjg(value(1))
                elseif (j .gt. ndimh .and. i .le. ndimh) then
                    amat(j - ndimh, i + ndimh) = amat(j - ndimh, i + ndimh) + conjg(value(1))
                    amat(i + ndimh, j - ndimh) = amat(i + ndimh, j - ndimh) - conjg(value(1))
                    !          elseif(j.le.ndimh.and.i.gt.ndimh) then
                    !             amat(j+ndimh, i-ndimh)=conjg(value(1))
                    !             amat(i-ndimh, j+ndimh)=-conjg(value(1))
                end if
            else
                if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                    amat(i + ndimh, j + ndimh) = amat(i + ndimh, j + ndimh) + value(1)
                    amat(j + ndimh, i + ndimh) = amat(j + ndimh, i + ndimh) - value(1)
                elseif (j .gt. ndimh .and. i .le. ndimh) then
                    amat(j - ndimh, i + ndimh) = amat(j - ndimh, i + ndimh) + value(1)
                    amat(i + ndimh, j - ndimh) = amat(i + ndimh, j - ndimh) - value(1)
                    !          elseif(j.le.ndimh.and.i.gt.ndimh) then
                    !             amat(j+ndimh, i-ndimh)=value(1)
                    !             amat(i-ndimh, j+ndimh)=-value(1)
                end if
            end if
        end if

    end if
    return
end
