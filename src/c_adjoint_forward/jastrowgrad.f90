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

subroutine jastrowgrad(rc, vj, iesd, jastrow, gradv, grad2, ispin)
    implicit none
    real*8 r, rc(3), jastrow, vj(*), grad, grad2, aux, aux2, aux3, peffm       &
            &, gradv(3), cost, acost, bcost, rz, dxrz, dxr, fat, rs, fat0
    integer iesd, ispin, j
    !        Output jastrow= log (two body Jastrow)
    !        grad2= 2/r u'  + u''
    !        grad= u'
    !        gradv= grad u

    !     input r output jastrow grad e grad2
    select case (iesd)

    case (4)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        fat = dexp(-vj(1)*r)
        if (ispin .lt. 0) then
            jastrow = 0.5d0/vj(1)*(1.d0 - fat)
            !        jastrow=dexp(jastrow)
            grad = 0.5d0*fat
            grad2 = 0.5d0*fat*(2.d0 - vj(1)*r)
        else
            jastrow = 0.25d0/vj(1)*(1.d0 - fat)
            !        jastrow=dexp(jastrow)
            grad = 0.25d0*fat
            grad2 = 0.25d0*fat*(2.d0 - vj(1)*r)
        end if

        grad2 = grad2/r
        cost = grad/r
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
    case (-4)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        fat = dexp(-vj(1)*r)
        jastrow = 0.5d0/vj(1)*(1.d0 - fat)
        !        jastrow=dexp(jastrow)
        grad = 0.5d0*fat
        grad2 = 0.5d0*fat*(2.d0 - vj(1)*r)
        grad2 = grad2/r
        cost = grad/r
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
    case (8)
        ! jastrow for pseudo soft
        r = dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2)
        fat = dexp(-vj(1)*r**3)
        jastrow = 1.d0/vj(1)*(1.d0 - fat)
        cost = fat*3*r
        grad2 = -3*r*fat*(-4 + 3*vj(1)*r**3)
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
    case (9)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        Acost = vj(2)*(1 + vj(1))/(4.d0*vj(1))
        fat0 = (1.d0 - r/vj(2))**2
        fat = 1.d0 + vj(1)*fat0
        jastrow = -Acost*log(fat)
        grad = 2.d0*Acost*vj(1)*(1.d0 - r/vj(2))/(vj(2)*fat)
        grad2 = -Acost/vj(2)**2*vj(1)*(-4.d0*vj(1)*fat0/fat**2 + 2.d0/fat)
        cost = grad/r
        grad2 = 2*cost + grad2
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
    case (10)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        Acost = vj(2)*(1 + vj(1))/(4.d0*vj(1))
        if (ispin .gt. 0) then
            fat0 = (1.d0 - 0.5d0*r/vj(2))**2
        else
            fat0 = (1.d0 - r/vj(2))**2
        end if
        fat = 1.d0 + vj(1)*fat0
        jastrow = -Acost*log(fat)
        if (ispin .gt. 0) then
            grad = Acost*vj(1)*(1.d0 - 0.5d0*r/vj(2))/(vj(2)*fat)
            grad2 = -Acost/vj(2)**2*vj(1)*(-vj(1)*fat0/fat**2 + 0.5d0/fat)
        else
            grad = 2.d0*Acost*vj(1)*(1.d0 - r/vj(2))/(vj(2)*fat)
            grad2 = -Acost/vj(2)**2*vj(1)*(-4.d0*vj(1)*fat0/fat**2 + 2.d0/fat)
        end if
        cost = grad/r
        grad2 = 2*cost + grad2
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
    case (5)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        aux = dexp(-vj(1)*r)
        aux2 = dexp(-vj(2)*r)
        aux3 = 0.5d0/(vj(1) + vj(3)*vj(2))
        jastrow = aux3*(1.d0 + vj(3) - aux - vj(3)*aux2)
        !        jastrow=dexp(jastrow)
        grad = aux3*(vj(1)*aux + vj(3)*vj(2)*aux2)
        grad2 = aux3*(2.d0*vj(2)*vj(3)*aux2 - vj(2)**2*vj(3)*aux2*r        &
                & - vj(1)*aux*(-2.d0 + vj(1)*r))
        grad2 = grad2/r
        cost = grad/r
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do

    case (1)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        aux = 1.d0/(1.d0 + vj(1)*r)
        aux2 = aux*aux
        aux3 = aux2*aux
        jastrow = 0.5d0*r*aux
        !           jastrow=dexp(jastrow)
        grad = 0.5d0*aux2
        grad2 = aux3/r
        cost = grad/r
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
    case (6)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        if (vj(2) .gt. 1d-9) then
            rs = (1.d0 - dexp(-vj(2)*r))/vj(2)
        else
            rs = r
        end if
        aux = 1.d0/(1.d0 + vj(1)*rs)
        aux2 = aux*aux
        aux3 = aux2*aux
        jastrow = 0.5d0*rs*aux
        !           jastrow=dexp(jastrow)
        grad = 0.5d0*aux2*dexp(-vj(2)*r)
        grad2 = (-vj(1)*aux3*dexp(-vj(2)*r)                           &
                & - 0.5d0*vj(2)*aux2)*dexp(-vj(2)*r)
        grad2 = 2.d0/r*grad + grad2
        cost = grad/r
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do

    case (7)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        if (vj(2) .gt. 1d-9) then
            rs = (1.d0 - dexp(-vj(2)*r))/vj(2)
        else
            rs = r
        end if
        aux = 1.d0/(1.d0 + vj(1)*rs)
        aux2 = aux*aux
        aux3 = aux2*aux

        ! parallel
        if (ispin .gt. 0) then
            jastrow = 0.25d0*rs*aux
            grad = 0.25d0*aux2*dexp(-vj(2)*r)
            grad2 = 0.5d0*(-vj(1)*aux3*dexp(-vj(2)*r)                  &
                    & - 0.5d0*vj(2)*aux2)*dexp(-vj(2)*r)

            ! antiparallel
        else
            jastrow = 0.5d0*rs*aux
            grad = 0.5d0*aux2*dexp(-vj(2)*r)
            grad2 = (-vj(1)*aux3*dexp(-vj(2)*r)                        &
                    & - 0.5d0*vj(2)*aux2)*dexp(-vj(2)*r)
        end if

        grad2 = 2.d0/r*grad + grad2
        cost = grad/r
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
        !           jastrow=dexp(jastrow)

    case (2)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)

        if (ispin .gt. 0) then
            aux = 1.d0/(1.d0 + vj(2)*r)
            aux2 = aux*aux
            aux3 = aux2*aux
            jastrow = 0.25d0*r*aux
            grad = 0.25d0*aux2
            grad2 = 0.5d0*(-vj(2)*aux3 + aux2/r)
            cost = grad/r
            do j = 1, 3
                gradv(j) = cost*rc(j)
            end do
        else
            aux = 1.d0/(1.d0 + vj(1)*r)
            aux2 = aux*aux
            aux3 = aux2*aux
            jastrow = 0.5d0*r*aux
            grad = 0.5d0*aux2
            grad2 = -vj(1)*aux3 + aux2/r
            cost = grad/r
            do j = 1, 3
                gradv(j) = cost*rc(j)
            end do
        end if
        !           jastrow=dexp(jastrow)

    case (-1)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)
        !        cusp condition for parallel electrons
        aux = 1.d0/(1.d0 + vj(1)*r)
        aux2 = aux*aux
        aux3 = aux2*aux
        if (ispin .gt. 0) then
            jastrow = 0.25d0*r*aux
            grad = 0.25d0*aux2
            !        grad2=-vj(1)*(aux+aux2+aux3)
            grad2 = 0.5d0*(-vj(1)*aux3 + aux2/r)
            cost = grad/r
            do j = 1, 3
                gradv(j) = cost*rc(j)
            end do
        else
            jastrow = 0.5d0*r*aux
            grad = 0.5d0*aux2
            !        grad2=-vj(1)*(aux+aux2+aux3)
            grad2 = -vj(1)*aux3 + aux2/r
            cost = grad/r
            do j = 1, 3
                gradv(j) = cost*rc(j)
            end do
        end if
        !           jastrow=dexp(jastrow)

    case (-2)

        rz = dsqrt(vj(1)**2*(rc(1)**2 + rc(2)**2) + (vj(2)*rc(3))**2)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)

        aux = 1.d0/(1.d0 + rz)
        jastrow = 0.5d0*r*aux
        !        jastrow=dexp(jastrow)
        aux = 0.5d0/(1.d0 + rz)**2/rz/r
        aux2 = vj(2)**2 - vj(1)**2

        acost = aux2*rc(3)**2
        bcost = -aux2*(rc(1)**2 + rc(2)**2)

        !        write(6,*) ' New laplace -2 '

        do j = 1, 2
            gradv(j) = rc(j)*(acost + rz)*aux
            !         write(6,*) ' Grad =',j,gradv(j)
        end do
        gradv(3) = rc(3)*(bcost + rz)*aux
        !         write(6,*) ' Grad =',3,gradv(3)

        ! ! the trivial laplacian term
        grad2 = (2.d0*acost + bcost + 3.d0*rz)*aux

        dxrz = -aux/rz/(1.d0 + rz)*(acost + 3.d0*acost*rz + 2.d0*rz**2)
        dxr = -aux/r*(acost + rz)

        grad2 = grad2 + (rc(1)**2 + rc(2)**2)*(dxr/r + dxrz/rz*vj(1)**2)

        dxrz = -aux/rz/(1.d0 + rz)*(bcost + 3.d0*bcost*rz + 2.d0*rz**2)
        dxr = -aux/r*(bcost + rz)

        grad2 = grad2 + rc(3)**2*(dxr/r + dxrz/rz*vj(2)**2)

    case (3)
        r = max(dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2), 1d-9)

        aux = 1.d0/(1.d0 + vj(1)*r)
        aux2 = aux*aux
        aux3 = aux2*aux
        jastrow = 0.5d0*r*aux
        grad = 0.5d0*aux2
        !        grad2=-vj(1)*(aux+aux2+aux3)
        aux = r/(1.d0 + vj(2)*r)**2
        aux2 = aux/(1.d0 + vj(2)*r)
        jastrow = jastrow + vj(3)*0.5d0*r*aux
        !     jastrow=dexp(jastrow)
        grad = grad + vj(3)*aux2
        aux3 = aux3 + 3.d0*vj(3)*aux/(1.d0 + vj(2)*r)**2
        grad2 = aux3/r
        cost = grad/r
        do j = 1, 3
            gradv(j) = cost*rc(j)
        end do
    case (0)
        jastrow = 0.d0
        grad = 0.d0
        !        call dscalzero(3,0.d0,gradv,1)
        gradv = 0.d0
        grad2 = 0.d0

    case default
#ifndef _OFFLOAD
        write (6, *) 'ERROR Jastrow: wrong jastrow number', iesd
        stop
#endif
    end select

    return
end
