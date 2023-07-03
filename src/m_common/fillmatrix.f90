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

subroutine fillmatrix_vec(rotmatrix, vecr)

    !      Written by Sorella, define a unitary matrix of a rotation (rotmatrix)
    !      leaving unchanged the axis defined by vecr, that is not changed.
    use Constants, only: Pi, zzero, zone, zimg
    implicit none

    real*8 x3, drand1, cost, rotmatrix(3, 3)
    real*8 vecr(3), vec(3), rwork(7), w(3)
    complex*16 amat(3, 3), work(5), ascra(3, 3), aexp(3, 3), wc(3)
    real*8 rotmat_vec(3, 3), rotmat_scra(3, 3)
    integer info, i, j

    cost = sum(vecr(:)**2)

    if (cost .gt. 0.d0) then
        vec(:) = vecr(:)/sqrt(cost)
    else
        vec = 0.d0
    end if

    amat = zzero
    ascra = zzero
    aexp = zzero

    amat(2, 3) = -vec(1)*zimg
    amat(1, 3) = vec(2)*zimg
    amat(1, 2) = -vec(3)*zimg

    amat(3, 2) = dconjg(amat(2, 3))
    amat(3, 1) = dconjg(amat(1, 3))
    amat(2, 1) = dconjg(amat(1, 2))

    call ZHEEV('V', 'U', 3, amat, 3, w, work, 5, rwork, info)
    x3 = (drand1() - 0.5d0)*2.d0*pi
    do i = 1, 3
        wc(i) = dcmplx(dcos(x3*w(i)), dsin(x3*w(i)))
    end do

    do j = 1, 3
        do i = 1, 3
            ascra(i, j) = amat(i, j)*wc(j)
        end do
    end do

    call zgemm('N', 'C', 3, 3, 3, zone, ascra, 3, amat, 3, zzero, aexp, 3)

    do j = 1, 3
        do i = 1, 3
            rotmat_vec(i, j) = aexp(i, j)
        end do
    end do
    rotmat_scra = rotmatrix
    call dgemm('N', 'N', 3, 3, 3, 1.d0, rotmat_vec, 3, rotmat_scra, 3, 0.d0, rotmatrix, 3)

    return
end
subroutine fillmatrix(rotmatrix)

    !      Written by Sorella, Mazzola and Y. Luo on 25/7/2013. Mitas
    !      routine was wrong unfortunately. All previous calculation are
    !      correct only in the limit a-->0
    use Constants, only: Pi
    implicit none

    real*8 x1, x2, x3, xsum, theta, phi, alpha                  &
            &, y1, y2, y3, z1, z2, z3, drand1, rotmatrix(3, 3), yval(3), yvalo(3)

    x3 = 1.d0 - 2.d0*drand1()

    phi = 2.d0*pi*drand1()

    theta = dacos(x3)

    x1 = cos(phi)*sin(theta)
    x2 = sin(phi)*sin(theta)
    x3 = cos(theta)
    !    Find an arbitrary solution y2 orthogonal to x in a stable way
    if (abs(x1) .lt. abs(x2)) then
        if (abs(x2) .lt. abs(x3)) then
            !         x3 is the maximum
            y1 = 1.d0
            y2 = 0.d0
            y3 = -x1/x3
        else
            !         x2 is the maximum
            y1 = 1.d0
            y3 = 0.d0
            y2 = -x1/x2
        end if
    else
        if (abs(x1) .lt. abs(x3)) then
            !         x3 is the maximum
            y1 = 1.d0
            y2 = 0.d0
            y3 = -x1/x3
        else
            !         x1 is the maximum
            y2 = 1.d0
            y3 = 0.d0
            y1 = -x2/x1
        end if
    end if
    !       Now rotate by an angle 0< phi1 < 2 Pi around x1
    alpha = 2.d0*pi*drand1()
    call make_u(alpha, x1, x2, x3, rotmatrix)
    yval(1) = y1
    yval(2) = y2
    yval(3) = y3
    call dgemv('N', 3, 3, 1.d0, rotmatrix, 3, yval, 1, 0.d0, yvalo, 1)
    xsum = dsqrt(sum(yvalo(:)**2))
    y1 = yvalo(1)/xsum
    y2 = yvalo(2)/xsum
    y3 = yvalo(3)/xsum

    !       vector product

    z1 = x2*y3 - y2*x3
    z2 = -(x1*y3 - x3*y1)
    z3 = x1*y2 - x2*y1

    xsum = z1**2 + z2**2 + z3**2
    if (abs(1 - xsum) .gt. 1d-6) write (6, *) ' ERROR in fillmatrix '

    rotmatrix(1, 1) = x1
    rotmatrix(2, 1) = x2
    rotmatrix(3, 1) = x3
    rotmatrix(1, 2) = y1
    rotmatrix(2, 2) = y2
    rotmatrix(3, 2) = y3
    rotmatrix(1, 3) = z1
    rotmatrix(2, 3) = z2
    rotmatrix(3, 3) = z3

    return
end
