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

subroutine evaldet(a, nl, nd, ipiv, detln, detsn)
    use constants, only: ipf

    implicit none

    integer nl, nd, isum, i, detsn
    real*8 a(nl, nd), detln
    integer ipiv(nd)
    if (ipf .eq. 1) then
        detsn = 1
        detln = 0.d0
        isum = 0
#ifdef _CUSOLVER
!$omp target teams distribute parallel do reduction(+:isum,detln) reduction(*:detsn)
#endif
        do i = 1, nd
            detln = log(abs(a(i, i))) + detln
            detsn = int(sign(1.d0, a(i, i)))*detsn
            if (ipiv(i) .ne. i) isum = isum + 1
        end do
#ifdef _CUSOLVER
! Workaround for nvfortran compiler. Does not update automatically scalars.
!$omp target update from(isum,detln,detsn)
#endif
        if (mod(isum, 2) .ne. 0) detsn = -detsn
    else
        detsn = int(sign(1.d0, a(1, 2)))
        detln = log(abs(a(1, 2)))

        isum = 0
        if (ipiv(1) .ne. 1) isum = 1
        do i = 2, nd
            if (ipiv(i) .ne. i) isum = isum + 1
        end do
        do i = 3, nd - 1, 2
            detln = log(abs(a(i, i + 1))) + detln
            detsn = int(sign(1.d0, a(i, i + 1)))*detsn
        end do
        if (mod(isum, 2) .ne. 0) detsn = -detsn

    end if

    return

end subroutine evaldet
