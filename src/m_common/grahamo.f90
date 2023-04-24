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

subroutine GRAHAMO(PSI, SC, RN, NDIME, MH, INFO)
    implicit none
    integer mh, ndime, mx, i, mhh, info
    real*8 PSI(NDIME, *), SC(MH + NDIME), RN(MH), dnrm2, cost, eps, epsmin
    real*8, external :: dlamch
    eps = 10000.d0*dlamch('E')
    epsmin = dlamch('S')
    MHH = MH + 1
    do I = 1, MH
        SC(I) = -1.d0
    end do

    !          First normalize the input orbitals
    do i = 1, MH
        COST = DNRM2(NDIME, PSI(1, i), 1)
        if (cost .gt. epsmin) then
            cost = 1.d0/cost
            call DSCAL(NDIME, COST, PSI(1, i), 1)
        else
            psi(:, i) = 0.d0
        end if
    end do
    RN(1) = DNRM2(NDIME, PSI, 1)
    if (RN(1) .gt. eps) then
        info = 0
    else
        info = 1
    end if
    do I = 2, MH
        call DGEMV('T', ndime, I - 1, 1.d0, PSI, ndime, PSI(1, I), 1, 0.d0, SC, 1)
        call DGEMV('N', NDIME, I, 1.d0, PSI, ndime, SC, 1, 0.d0, SC(MHH), 1)
        RN(I) = DNRM2(NDIME, SC(MHH), 1)
        if (RN(I) .gt. eps) then
            COST = -1.d0/RN(I)
            call DSCAL(NDIME, COST, SC(MHH), 1)
            call DCOPY(NDIME, SC(MHH), 1, PSI(1, I), 1)
        else
            psi(:, i) = 0.d0
            info = info + 1
        end if

    end do
    return
end subroutine GRAHAMO

subroutine GRAHAMO_COMPLEX(PSI, SC, RN, NDIME, MH, INFO)
    use constants, only: zone, zmone, zzero
    implicit none
    integer mh, ndime, mx, i, mhh, info
    complex*16 PSI(NDIME, *), SC(MH + NDIME), cost
    real*8 RN(MH), eps, epsmin
    real*8, external :: dlamch, dznrm2
    eps = 10000.d0*dlamch('E') ! Relative machine precision
    epsmin = dlamch('S') ! Safe minimum that allows the inverse.

    MHH = MH + 1
    do I = 1, MH
        SC(I) = zmone
    end do
    !          First normalize the input orbitals
    do i = 1, MH
        COST = DZNRM2(NDIME, PSI(1, i), 1)
        if (dreal(cost) .gt. epsmin) then
            cost = zone/cost
            call ZSCAL(NDIME, COST, PSI(1, i), 1)
        else
            psi(:, i) = zzero
        end if
    end do

    RN(1) = DZNRM2(NDIME, PSI, 1)
    if (RN(1) .gt. eps) then
        info = 0
    else
        info = 1
    end if
    do I = 2, MH
        call ZGEMV('C', ndime, I - 1, zone, PSI, ndime, PSI(1, I), 1, zzero, SC, 1)
        call ZGEMV('N', NDIME, I, zone, PSI, ndime, SC, 1, zzero, SC(MHH), 1)
        RN(I) = DZNRM2(NDIME, SC(MHH), 1)
        if (RN(I) .gt. eps) then
            COST = zmone/RN(I)
            call ZSCAL(NDIME, COST, SC(MHH), 1)
            call ZCOPY(NDIME, SC(MHH), 1, PSI(1, I), 1)
        else
            psi(:, i) = zzero
            info = info + 1
        end if

    end do
    return
end
