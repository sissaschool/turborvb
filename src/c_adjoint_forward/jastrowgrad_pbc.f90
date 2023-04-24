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

subroutine jastrowgrad_pbc(rc, vj, iesd, jastrow, gradv, grad2, ispin, constz)

    use Cell, only: cellscale, metric, map, dmap, ddmap, car2cry
    use allio, only: niesd, norm_metric
    ! map: mapping function acting on the crystal coordinates
    ! dmap: first derivative of the mapping function
    ! ddmap: second derivative of the mapping function

    implicit none

    real(8), intent(in) :: vj(max(1, niesd)), constz
    integer, intent(in) :: iesd, ispin
    real(8) :: rc(3), jastrow, grad, grad2, grad2_before &
               , gradv(3), r, rc_map(3), drc_map(3), ddrc_map(3)
    !        Output jastrow= log (two body Jastrow)
    !        grad= u'
    !        grad2= u''

    !        gradv --> final gradient
    !        grad2 --> final laplacian

    integer k, niesdu
    niesdu = max(niesd, 1)
    niesdu = min(2, niesdu) ! at most 2 parameters

    if (constz .eq. 0.d0) then
        !     define the output
        gradv = 0.d0
        grad2 = 0.d0
        jastrow = 0.d0
        return
    end if

    ! working in the crystal reference
    !  call CartesianToCrystal(rc,1)
    rc(:) = car2cry(:, 1)*rc(1) + car2cry(:, 2)*rc(2) + car2cry(:, 3)*rc(3)
    do k = 1, 3
        ! mapping function and its derivatives acting on crystal coordinates here
        rc_map(k) = map(rc(k), cellscale(k))
        drc_map(k) = dmap(rc(k), cellscale(k))
        ddrc_map(k) = ddmap(rc(k), cellscale(k))
        rc_map(k) = constz*rc_map(k)
        drc_map(k) = constz*drc_map(k)
        ddrc_map(k) = constz*ddrc_map(k)
    end do

    ! use metric to perform the square
    r = norm_metric(rc_map, metric)
    r = max(r, 1d-9)

    call update_jgrad(iesd, ispin, r, vj, niesdu, jastrow, grad, grad2)

    !  input r,jastrow_before,grad,grad2_before,rc_map,drc_map,ddrc_map,metric,car2cry
    !  output gradv,grad2
    grad2_before = grad2

    call update_all(r, grad, grad2_before, rc_map, drc_map, ddrc_map, metric, car2cry, gradv, grad2, constz)

    return
end subroutine jastrowgrad_pbc

subroutine update_all(r, grad, grad2_before, rc_map, drc_map, ddrc_map, metric, car2cry, gradv, grad2, constz)
    implicit none
    real(8) :: grad, grad2, grad_before, grad2_before &
               , gradv(3), cost, r, hess(3), rc_map(3), drc_map(3), ddrc_map(3) &
               , vprime(3), vsecond(3), vprime_prod(3, 3) &
               , metric(3, 3), car2cry(3, 3), constz
    integer :: i, j
    ! input -->r,grad,grad2_before,rc_map,drc_map,ddrc_map,metric,car2cry
    ! output --> gradv,grad2
    ! constants--> constz

    cost = grad/r
    grad2 = (grad2_before - cost)/r**2

    ! once radial functions computed, compute now gradient and hessian
    ! vprime(:)=0.d0
    ! do i=1,3
    ! vprime(:)=vprime(:)+metric(:,i)*rc_map(i)
    ! enddo
    ! call dgemv('N',3,3,1.d0,metric,3,rc_map,1,0.d0,vprime,1)
    vprime(:) = metric(:, 1)*rc_map(1) + metric(:, 2)*rc_map(2) + metric(:, 3)*rc_map(3)
    vsecond(:) = ddrc_map(:)*vprime(:)
    vprime(:) = drc_map(:)*vprime(:)
    do i = 1, 3
        do j = i + 1, 3
            vprime_prod(i, j) = drc_map(i)*metric(i, j)*drc_map(j)
            vprime_prod(j, i) = vprime_prod(i, j)
        end do
        vprime_prod(i, i) = drc_map(i)**2*metric(i, i)
    end do

    ! gradv(:)=0.d0
    ! do i=1,3
    !    gradv(:)=gradv(:)+car2cry(i,:)*vprime(i)
    ! enddo
    ! call dgemv('T',3,3,1.d0,car2cry,3,vprime,1,0.d0,gradv,1)
    gradv(:) = car2cry(1, :)*vprime(1) + car2cry(2, :)*vprime(2) + car2cry(3, :)*vprime(3)
!   hess = gradv
!   gradv(:) = gradv(:) * cost
    hess(:) = grad2*gradv(:)**2

    do i = 1, 3
        do j = i + 1, 3
            hess(:) = hess(:) + 2.d0*cost*car2cry(i, :)*car2cry(j, :)*vprime_prod(i, j)
        end do
        hess(:) = hess(:) + cost*car2cry(i, :)**2*(vprime_prod(i, i) + vsecond(i))
    end do

    ! do i=1,3
    !    hess(:)=hess(:)+cost*car2cry(i,:)**2*vsecond(i)
    ! enddo

    gradv(:) = gradv(:)*cost/constz
!   hess(:) = hess(:) / constz**2
    grad2 = sum(hess(:)/constz**2)

end subroutine update_all
subroutine update_jgrad(iesd, ispin, r, vj, niesd, jastrow, grad, grad2)
    implicit none
    integer iesd, ispin, niesd
    real*8 r, jastrow, grad, grad2
    real*8 fat, aux, aux2, aux3, vj(niesd)
    ! input r output jastrow, grad, grad2

    select case (iesd)

        ! compute radial Jastrow function and its first and second derivatives at r

        ! done only for case(4) and case(-4), work in progress for the others
    case (4)
        fat = dexp(-vj(1)*r)
        if (ispin .lt. 0) then
            jastrow = 0.5d0/vj(1)*(1.d0 - fat)
            grad = 0.5d0*fat
            grad2 = -0.5d0*fat*vj(1)
        else
            jastrow = 0.25d0/vj(1)*(1.d0 - fat)
            grad = 0.25d0*fat
            grad2 = -0.25d0*fat*vj(1)
        end if

    case (-4)
        fat = dexp(-vj(1)*r)
        jastrow = 0.5d0/vj(1)*(1.d0 - fat)
        grad = 0.5d0*fat
        grad2 = -0.5d0*fat*vj(1)

    case (8)

        fat = dexp(-vj(1)*r**3)
        jastrow = 1.d0/vj(1)*(1.d0 - fat)
        grad = 3.d0*r**2*fat
        grad2 = (6.d0*r - 9.d0*vj(1)*r**4)*fat

    case (1)
        aux = 1.d0/(1.d0 + vj(1)*r)
        aux2 = aux*aux
        aux3 = aux2*aux
        jastrow = 0.5d0*r*aux
        grad = 0.5d0*aux2
        grad2 = -vj(1)*aux3

    case (2)
        if (ispin .gt. 0) then
            aux = 1.d0/(1.d0 + vj(2)*r)
            aux2 = aux*aux
            aux3 = 0.5d0*aux2*aux
            jastrow = 0.25d0*r*aux
            grad = 0.25d0*aux2
            grad2 = -vj(2)*aux3
        else
            aux = 1.d0/(1.d0 + vj(1)*r)
            aux2 = aux*aux
            aux3 = aux2*aux
            jastrow = 0.5d0*r*aux
            grad = 0.5d0*aux2
            grad2 = -vj(1)*aux3
        end if

    case (-1)
        !        cusp condition for parallel electrons
        aux = 1.d0/(1.d0 + vj(1)*r)
        aux2 = aux*aux
        aux3 = aux2*aux
        if (ispin .gt. 0) then
            jastrow = 0.25d0*r*aux
            grad = 0.25d0*aux2
            grad2 = -0.5d0*vj(1)*aux3
        else
            jastrow = 0.5d0*r*aux
            grad = 0.5d0*aux2
            grad2 = -vj(1)*aux3
        end if

    case (9)
        ! not coded yet

        aux2 = vj(2)*(1 + vj(1))/(4.d0*vj(1))
        aux = (1.d0 - r/vj(2))**2
        fat = 1.d0 + vj(1)*aux
        jastrow = -aux2*log(fat)
        grad = 2.d0*aux2*vj(1)*(1.d0 - r/vj(2))/(vj(2)*fat)
        grad2 = -aux2/vj(2)**2*vj(1)*(-4.d0*vj(1)*aux/fat**2 + 2.d0/fat)

    case (10)
        ! not coded yet

        aux2 = vj(2)*(1 + vj(1))/(4.d0*vj(1))
        if (ispin .gt. 0) then
            aux = (1.d0 - 0.5d0*r/vj(2))**2
        else
            aux = (1.d0 - r/vj(2))**2
        end if
        fat = 1.d0 + vj(1)*aux
        jastrow = -aux2*log(fat)
        if (ispin .gt. 0) then
            grad = aux2*vj(1)*(1.d0 - 0.5d0*r/vj(2))/(vj(2)*fat)
            grad2 = -aux2/vj(2)**2*vj(1)*(-vj(1)*aux/fat**2 + 0.5d0/fat)
        else
            grad = 2.d0*aux2*vj(1)*(1.d0 - r/vj(2))/(vj(2)*fat)
            grad2 = -aux2/vj(2)**2*vj(1)*(-4.d0*vj(1)*aux/fat**2 + 2.d0/fat)
        end if

    case (0)
        jastrow = 0.d0
        grad = 0.d0
        grad2 = 0.d0

    case default
#ifndef _OFFLOAD
        write (6, *) 'ERROR Jastrow: wrong jastrow number'
        stop
#endif
    end select

end subroutine update_jgrad
