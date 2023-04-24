!TL off
! written by hand by S. Sorella on 11/8/2019
SUBROUTINE JASTROWGRAD_PBC_B(rc, rcb, vj, vjb, iesd, jastrow, gradv&
        &, gradvb, grad2, grad2b, ispin, constz, cellscaleb, metricb, car2cryb)
    USE CONSTANTS, ONLY : PI
    USE CELL, ONLY : cellscale, car2cry, metric, map, dmap, ddmap, CartesianToCrystal, CartesianToCrystal_b
    use allio, only : niesd, norm_metric
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: vj(max(niesd, 1)), constz
    INTEGER, INTENT(IN) :: iesd, ispin
    real(8) :: rc(3), jastrow, grad, grad2, grad_before, grad2_before&
            , gradv(3), r, rc_map(3), drc_map(3), ddrc_map(3)
    real*8 rcb(3), gradb, grad2b, grad2_beforeb, rc_before(3)&
            , gradvb(3), rb, rc_mapb(3), drc_mapb(3), ddrc_mapb(3)
    real*8 vjb(*), cellscaleb(3), metricb(3, 3), car2cryb(3, 3)
    integer k, niesdu




    !  first step vanish all adjoint arrays  not passed
    grad2_beforeb = 0.d0
    gradb = 0.d0
    rc_mapb = 0.d0
    drc_mapb = 0.d0
    ddrc_mapb = 0.d0
    niesdu = max(niesd, 1)
    niesdu = min(2, niesdu)

    !     fast output

    if(constz.eq.0) return

    ! working in the crystal reference
    rc_before = rc
    call CartesianToCrystal(rc, 1)

    do k = 1, 3
        ! mapping function and its derivatives acting on crystal coordinates here
        rc_map(k) = map(rc(k), cellscale(k))
        drc_map(k) = dmap(rc(k), cellscale(k))
        ddrc_map(k) = ddmap(rc(k), cellscale(k))
        rc_map(k) = constz * rc_map(k)
        drc_map(k) = constz * drc_map(k)
        ddrc_map(k) = constz * ddrc_map(k)
    enddo

    ! use metric to perform the square
    r = norm_metric(rc_map, metric)
    if(r.le.1.d-9) return ! there is nothing to differentiate

    call update_jgrad(iesd, ispin, r, vj, niesdu, jastrow, grad, grad2)

    grad2_before = grad2

    !  BEGIN reverse algorithm
    call  update_all_b(r, rb, grad, gradb, grad2_before, grad2_beforeb, rc_map, rc_mapb, drc_map, drc_mapb, ddrc_map, ddrc_mapb, metric, metricb, car2cry, car2cryb, gradvb, grad2b, constz)

    !  reverse of
    !  grad2_before=grad2

    grad2b = grad2b + grad2_beforeb

    Call  UPDATE_JGRAD_B(iesd, ispin, r, rb, vj, vjb, jastrow&
            &, grad, gradb, grad2, grad2b)



    !  reverse of r=norm_metric(rc_map,metric)

    call norm_metric_b(rc_map, rc_mapb, metric, metricb, rb)

    ! reverse of

    !  do k=1,3
    !     ! mapping function and its derivatives acting on crystal coordinates here
    !     rc_map(k)=constz*rc_map(k)
    !     drc_map(k)=constz*drc_map(k)
    !     ddrc_map(k)=constz*ddrc_map(k)
    !  enddo

    do k = 1, 3
        rc_mapb(k) = rc_mapb(k) * constz
        drc_mapb(k) = drc_mapb(k) * constz
        ddrc_mapb(k) = ddrc_mapb(k) * constz
    enddo


    ! reverse of
    ! do k=1,3
    !    rc_map(k)=map(rc(k),cellscale(k))
    !    drc_map(k)=dmap(rc(k),cellscale(k))
    !    ddrc_map(k)=ddmap(rc(k),cellscale(k))
    ! enddo
    do k = 1, 3
        call ddmap_b(rc(k), rcb(k), cellscale(k), cellscaleb(k), ddrc_mapb(k))
        call dmap_b(rc(k), rcb(k), cellscale(k), cellscaleb(k), drc_mapb(k))
        call map_b(rc(k), rcb(k), cellscale(k), cellscaleb(k), rc_mapb(k))
    enddo

    call CartesianToCrystal_b(rc_before, rcb, car2cryb, 1)

END SUBROUTINE JASTROWGRAD_PBC_B


subroutine update_all_b(r, rb, grad, gradb, grad2_before, grad2_beforeb, rc_map, rc_mapb, drc_map, drc_mapb, ddrc_map, ddrc_mapb, metric, metricb, car2cry, car2cryb, gradvb, grad2b, constz)
    implicit none
    real(8) :: grad, grad2, grad_before, grad2_before&
            , gradv(3), cost, r, rc_map(3), drc_map(3), ddrc_map(3) &
            , vprime(3), vsecond(3), vprime_prod(3, 3)&
            , metric(3, 3), car2cry(3, 3), constz
    real(8) :: gradb, grad2b, grad_beforeb, grad2_beforeb&
            , gradvb(3), costb, rb, hessb(3), rc_mapb(3), drc_mapb(3), ddrc_mapb(3) &
            , vprimeb(3), vsecondb(3), vprime_prodb(3, 3)&
            , metricb(3, 3), car2cryb(3, 3)
    real(8) :: aux(3), vprime_before(3), aux0
    integer i, j

    ! input -->r,grad,grad2_before,rc_map,drc_map,ddrc_map,metric,car2cry
    ! output --> gradv,grad2
    ! constants--> constz

    cost = grad / r
    grad2 = (grad2_before - cost) / r**2

!   grad2_after1 = grad2
    ! once radial functions computed, compute now gradient and hessian
    ! vprime(:)=0.d0
    ! do i=1,3
    ! vprime(:)=vprime(:)+metric(:,i)*rc_map(i)
    ! enddo
    call dgemv('N', 3, 3, 1.d0, metric, 3, rc_map, 1, 0.d0, vprime, 1)
    vprime_before = vprime
    vsecond(:) = ddrc_map(:) * vprime(:)
    vprime(:) = drc_map(:) * vprime(:)
    do i = 1, 3
        do j = i + 1, 3
            vprime_prod(i, j) = drc_map(i) * metric(i, j) * drc_map(j)
            vprime_prod(j, i) = vprime_prod(i, j)
        enddo
        vprime_prod(i, i) = drc_map(i)**2 * metric(i, i)
    enddo


    ! gradv(:)=0.d0
    ! do i=1,3
    !    gradv(:)=gradv(:)+car2cry(i,:)*vprime(i)
    ! enddo
    call dgemv('T', 3, 3, 1.d0, car2cry, 3, vprime, 1, 0.d0, gradv, 1)
!   hess(:) = grad2 * gradv(:)**2

!   do i = 1, 3
!       do j = i + 1, 3
!           hess(:) = hess(:) + 2.d0 * cost * car2cry(i, :) * car2cry(j, :) * vprime_prod(i, j)
!       enddo
!       hess(:) = hess(:) + cost * car2cry(i, :)**2 * (vprime_prod(i, i) + vsecond(i))
!   enddo

    ! do i=1,3
    !    hess(:)=hess(:)+cost*car2cry(i,:)**2*vsecond(i)
    ! enddo

!   done before updating gradv
!   gradv(:) = gradv(:)*cost/ constz
!   hess(:) = hess(:) / constz**2
!   grad2 = sum(hess(:)/constz**2)

    ! BEGIN REVERSE ALGORITHM

!   reverse of  grad2 = sum(hess(:)/constz**2)
!   hessb(:) = grad2b
    hessb(:) = grad2b/ constz**2
    grad2b = 0.d0
!  reverse of gradv(:) = gradv(:)*cost/ constz
    costb=sum(gradv(:)*gradvb(:))/constz
    gradvb(:) = gradvb(:)*cost/ constz
    ! vanish all internal variables
    rc_mapb = 0.d0
    drc_mapb = 0.d0
    ddrc_mapb = 0.d0
    vprimeb = 0.d0
    vsecondb = 0.d0
    vprime_prodb = 0.d0

    do i = 1, 3
        !    hess(:)=hess(:)+cost*car2cry(i,:)**2*(vprime_prod(i,i)+vsecond(i))

        aux(:) = car2cry(i, :)**2 * (vprime_prod(i, i) + vsecond(i))
        costb = costb + sum(aux(:) * hessb(:))
        aux(:) = 2.d0 * cost * car2cry(i, :) * (vprime_prod(i, i) + vsecond(i))
        car2cryb(i, :) = car2cryb(i, :) + aux(:) * hessb(:)
        aux0 = cost * sum(car2cry(i, :)**2 * hessb(:))
        vprime_prodb(i, i) = vprime_prodb(i, i) + aux0
        vsecondb(i) = vsecondb(i) + aux0

        do j = i + 1, 3

            !       hess(:)=hess(:)+2.d0*cost*car2cry(i,:)*car2cry(j,:)*vprime_prod(i,j)

            aux(:) = 2.d0 * car2cry(i, :) * car2cry(j, :) * vprime_prod(i, j)
            costb = costb + sum(aux(:) * hessb(:))
            car2cryb(i, :) = car2cryb(i, :) + 2.d0 * cost * car2cry(j, :) * vprime_prod(i, j) * hessb(:)
            car2cryb(j, :) = car2cryb(j, :) + 2.d0 * cost * car2cry(i, :) * vprime_prod(i, j) * hessb(:)
            aux0 = 2.d0 * cost * sum(car2cry(i, :) * car2cry(j, :) * hessb(:))
            vprime_prodb(i, j) = vprime_prodb(i, j) + aux0
        enddo
    enddo
    ! reverse of
    ! hess(:)=grad2*gradv(:)**2
    grad2b = grad2b + sum(gradv(:)**2 * hessb(:))
    gradvb(:) =gradvb(:)+2.d0 * grad2 * gradv(:) * hessb(:)
    hessb(:) = 0.d0
    ! reverse of
    !  call dgemv('T',3,3,1.d0,car2cry,3,vprime,1,0.d0,gradv,1)
    call dgemv_b('T', 3, 3, 1.d0, car2cry, 3, car2cryb, 3, vprime, 1, vprimeb, 1, 0.d0, gradvb, 1)
    ! reverse of
    ! do i=1,3
    !    do j=i+1,3
    !       vprime_prod(i,j)=drc_map(i)*metric(i,j)*drc_map(j)
    !       vprime_prod(j,i)=vprime_prod(i,j)
    !    enddo
    !    vprime_prod(i,i)=drc_map(i)**2*metric(i,i)
    ! enddo
    do i = 1, 3
        !    vprime_prod(i,i)=drc_map(i)**2*metric(i,i)
        drc_mapb(i) = drc_mapb(i) + 2.d0 * drc_map(i) * metric(i, i) * vprime_prodb(i, i)
        metricb(i, i) = metricb(i, i) + drc_map(i)**2 * vprime_prodb(i, i)
        do j = i + 1, 3
            vprime_prodb(i, j) = vprime_prodb(j, i) + vprime_prodb(i, j)
            drc_mapb(j) = drc_mapb(j) + vprime_prodb(i, j) * drc_map(i) * metric(i, j)
            drc_mapb(i) = drc_mapb(i) + vprime_prodb(i, j) * drc_map(j) * metric(i, j)
            metricb(i, j) = metricb(i, j) + vprime_prodb(i, j) * drc_map(i) * drc_map(j)
        enddo
    enddo
    ! reverse of   vprime(:)=drc_map(:)*vprime(:)
    drc_mapb(:) = drc_mapb(:) + vprime_before(:) * vprimeb(:)
    vprimeb(:) = vprimeb(:) * drc_map(:)
    !  reverse of vsecond(:)=ddrc_map(:)*vprime(:)
    vprimeb(:) = vprimeb(:) + vsecondb(:) * ddrc_map(:)
    ddrc_mapb(:) = ddrc_mapb(:) + vsecondb(:) * vprime_before(:)
    ! reverse of call dgemv('N',3,3,1.d0,metric,3,rc_map,1,0.d0,vprime,1)
    call dgemv_b('N', 3, 3, 1.d0, metric, 3, metricb, 3, rc_map, 1, rc_mapb, 1, 0.d0, vprimeb, 1)
    ! reverse of
    ! grad2=(grad2_before-cost)/r**2
    costb = -grad2b / r**2 + costb
    grad2_beforeb = grad2b / r**2
    rb = rb - 2.d0 * (grad2_before - cost) / r**3 * grad2b

    ! reverse of
    !  cost=grad/r
    gradb = gradb + costb / r
    rb = rb - grad / r**2 * costb

    grad2b = 0.d0
    gradvb = 0.d0
    !
end subroutine update_all_b
!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.3 (r3163) - 09/25/2009 09:03
!
!  Differentiation of update_jgrad in reverse (adjoint) mode:
!   gradient, with respect to input variables: grad2 r grad vj
!   of linear combination of output variables: grad2 r grad vj
SUBROUTINE UPDATE_JGRAD_B(iesd, ispin, r, rb, vj, vjb, jastrow, grad, &
        &  gradb, grad2, grad2b)
    IMPLICIT NONE
    INTEGER :: iesd, ispin
    REAL*8 :: r, jastrow, grad, grad2
    REAL*8 :: rb, gradb, grad2b
    REAL*8 :: fat, aux, aux2, aux3, vj(*)
    REAL*8 :: fatb, auxb, aux2b, aux3b, vjb(*)
    INTEGER :: branch
    REAL*8 :: temp3
    REAL*8 :: temp2
    REAL*8 :: temp1
    REAL*8 :: temp0
    REAL*8 :: temp13b
    REAL*8 :: temp21b
    REAL*8 :: temp24
    INTRINSIC DEXP
    REAL*8 :: temp23
    REAL*8 :: temp22
    REAL*8 :: temp9b0
    REAL*8 :: temp21
    REAL*8 :: temp20
    REAL*8 :: temp19b
    REAL*8 :: tempb1
    REAL*8 :: tempb0
    REAL*8 :: temp3b
    REAL*8 :: temp19
    REAL*8 :: temp18
    REAL*8 :: temp17
    REAL*8 :: temp16
    REAL*8 :: temp12b
    REAL*8 :: temp6b
    REAL*8 :: temp15
    REAL*8 :: temp20b
    REAL*8 :: temp14
    REAL*8 :: temp13
    REAL*8 :: temp12
    REAL*8 :: temp11
    REAL*8 :: temp10
    REAL*8 :: temp15b
    REAL*8 :: temp9b
    REAL*8 :: temp18b
    REAL*8 :: tempb
    REAL*8 :: temp2b
    REAL*8 :: temp11b
    REAL*8 :: temp5b
    INTRINSIC LOG
    REAL*8 :: temp14b
    REAL*8 :: temp22b
    REAL*8 :: temp1b
    REAL*8 :: temp
    REAL*8 :: temp1b5
    REAL*8 :: temp1b4
    REAL*8 :: temp1b3
    REAL*8 :: temp1b2
    REAL*8 :: temp9
    REAL*8 :: temp1b1
    REAL*8 :: temp8
    REAL*8 :: temp1b0
    REAL*8 :: temp7
    REAL*8 :: temp10b
    REAL*8 :: temp6
    REAL*8 :: temp4b
    REAL*8 :: temp5
    REAL*8 :: temp4
    ! input r output jastrow, grad, grad2
    SELECT CASE  (iesd)
    CASE (4)
        ! compute radial Jastrow function and its first and second derivatives at r
        ! done only for case(4) and case(-4), work in progress for the others
        fat = DEXP(-(vj(1) * r))
        IF (ispin .LT. 0) THEN
            fatb = 0.5d0 * gradb - 0.5d0 * vj(1) * grad2b
            vjb(1) = vjb(1) - 0.5d0 * fat * grad2b
        ELSE
            fatb = 0.25d0 * gradb - 0.25d0 * vj(1) * grad2b
            vjb(1) = vjb(1) - 0.25d0 * fat * grad2b
        END IF
        tempb = fat * fatb
        vjb(1) = vjb(1) - r * tempb
        rb = rb - vj(1) * tempb
    CASE (-4)
        fat = DEXP(-(vj(1) * r))
        fatb = 0.5d0 * gradb - 0.5d0 * vj(1) * grad2b
        tempb0 = fat * fatb
        vjb(1) = vjb(1) - r * tempb0 - 0.5d0 * fat * grad2b
        rb = rb - vj(1) * tempb0
    CASE (8)
        fat = DEXP(-(vj(1) * r**3))
        temp1b = fat * grad2b
        temp0 = r**4
        fatb = 3.d0 * r**2 * gradb + (6.d0 * r - 9.d0 * (vj(1) * temp0)) * grad2b
        temp = r**3
        tempb1 = fat * fatb
        rb = rb + 3.d0 * fat * 2 * r * gradb - vj(1) * 3 * r**2 * tempb1 + (6.d0 - 9.d0 * vj(1&
                &) * 4 * r**3) * temp1b
        vjb(1) = vjb(1) - temp * tempb1 - 9.d0 * temp0 * temp1b
    CASE (1)
        aux = 1.d0 / (1.d0 + vj(1) * r)
        aux2 = aux * aux
        aux3 = aux2 * aux
        aux3b = -(vj(1) * grad2b)
        aux2b = aux * aux3b + 0.5d0 * gradb
        auxb = 2 * aux * aux2b + aux2 * aux3b
        temp1b0 = -(auxb / (vj(1) * r + 1.d0)**2)
        vjb(1) = vjb(1) + r * temp1b0 - aux3 * grad2b
        rb = rb + vj(1) * temp1b0
    CASE (2)
        IF (ispin .GT. 0) THEN
            aux = 1.d0 / (1.d0 + vj(2) * r)
            aux2 = aux * aux
            aux3 = 0.5d0 * aux2 * aux
            aux3b = -(vj(2) * grad2b)
            aux2b = 0.5d0 * aux * aux3b + 0.25d0 * gradb
            auxb = 2 * aux * aux2b + 0.5d0 * aux2 * aux3b
            temp1b1 = -(auxb / (vj(2) * r + 1.d0)**2)
            vjb(2) = vjb(2) + r * temp1b1 - aux3 * grad2b
            rb = rb + vj(2) * temp1b1
        ELSE
            aux = 1.d0 / (1.d0 + vj(1) * r)
            aux2 = aux * aux
            aux3 = aux2 * aux
            aux3b = -(vj(1) * grad2b)
            aux2b = aux * aux3b + 0.5d0 * gradb
            auxb = 2 * aux * aux2b + aux2 * aux3b
            temp1b2 = -(auxb / (vj(1) * r + 1.d0)**2)
            vjb(1) = vjb(1) + r * temp1b2 - aux3 * grad2b
            rb = rb + vj(1) * temp1b2
        END IF
    CASE (-1)
        !        cusp condition for parallel electrons
        aux = 1.d0 / (1.d0 + vj(1) * r)
        aux2 = aux * aux
        aux3 = aux2 * aux
        IF (ispin .GT. 0) THEN
            vjb(1) = vjb(1) - 0.5d0 * aux3 * grad2b
            aux3b = -(0.5d0 * vj(1) * grad2b)
            aux2b = 0.25d0 * gradb
        ELSE
            vjb(1) = vjb(1) - aux3 * grad2b
            aux3b = -(vj(1) * grad2b)
            aux2b = 0.5d0 * gradb
        END IF
        aux2b = aux2b + aux * aux3b
        auxb = 2 * aux * aux2b + aux2 * aux3b
        temp1b3 = -(auxb / (vj(1) * r + 1.d0)**2)
        vjb(1) = vjb(1) + r * temp1b3
        rb = rb + vj(1) * temp1b3
    CASE (9)
        ! not coded yet
        aux2 = vj(2) * (1 + vj(1)) / (4.d0 * vj(1))
        aux = (1.d0 - r / vj(2))**2
        fat = 1.d0 + vj(1) * aux
        temp8 = fat**2
        temp5 = vj(1) * aux / temp8
        temp7 = vj(2)**2
        temp6b = -((2.d0 / fat - 4.d0 * temp5) * grad2b / temp7)
        temp6 = aux2 * vj(1) / temp7
        temp5b = temp6 * 4.d0 * grad2b / temp8
        temp2 = vj(2) * fat
        temp4 = r / vj(2)
        temp3b = 2.d0 * (1.d0 - temp4) * gradb / temp2
        aux2b = vj(1) * temp3b + vj(1) * temp6b
        vjb(1) = vjb(1) + aux * temp5b + aux2 * temp6b
        temp3 = aux2 * vj(1) / temp2
        temp4b = -(2.d0 * temp3 * gradb / vj(2))
        temp2b = -(temp3 * temp3b)
        vjb(2) = vjb(2) + fat * temp2b - temp4 * temp4b - temp6 * 2 * vj(2) * temp6b
        fatb = vj(2) * temp2b - temp5 * 2 * fat * temp5b + temp6 * 2.d0 * grad2b / fat**2
        auxb = vj(1) * fatb + vj(1) * temp5b
        temp1 = r / vj(2)
        temp1b4 = -(2 * (1.d0 - temp1) * auxb / vj(2))
        rb = rb + temp1b4 + temp4b
        vjb(1) = vjb(1) + aux * fatb + aux2 * temp3b
        temp1b5 = aux2b / (4.d0 * vj(1))
        vjb(2) = vjb(2) + (vj(1) + 1) * temp1b5 - temp1 * temp1b4
        vjb(1) = vjb(1) + (vj(2) - vj(2) * (vj(1) + 1) / vj(1)) * temp1b5
    CASE (10)
        ! not coded yet
        aux2 = vj(2) * (1 + vj(1)) / (4.d0 * vj(1))
        IF (ispin .GT. 0) THEN
            aux = (1.d0 - 0.5d0 * r / vj(2))**2
        ELSE
            aux = (1.d0 - r / vj(2))**2
        END IF
        fat = 1.d0 + vj(1) * aux
        IF (ispin .GT. 0) THEN
            temp17 = fat**2
            temp14 = vj(1) * aux / temp17
            temp16 = vj(2)**2
            temp15b = -((0.5d0 / fat - temp14) * grad2b / temp16)
            temp15 = aux2 * vj(1) / temp16
            temp14b = temp15 * grad2b / temp17
            temp11 = vj(2) * fat
            temp13 = r / vj(2)
            temp12b = (1.d0 - 0.5d0 * temp13) * gradb / temp11
            aux2b = vj(1) * temp12b + vj(1) * temp15b
            vjb(1) = vjb(1) + aux * temp14b + aux2 * temp15b
            temp12 = aux2 * vj(1) / temp11
            temp13b = -(temp12 * 0.5d0 * gradb / vj(2))
            temp11b = -(temp12 * temp12b)
            vjb(2) = vjb(2) + fat * temp11b - temp13 * temp13b - temp15 * 2 * vj(2) * &
                    &        temp15b
            fatb = vj(2) * temp11b - temp14 * 2 * fat * temp14b + temp15 * 0.5d0 * grad2b / &
                    &        fat**2
            auxb = vj(1) * temp14b
            rb = rb + temp13b
            vjb(1) = vjb(1) + aux2 * temp12b
        ELSE
            temp24 = fat**2
            temp21 = vj(1) * aux / temp24
            temp23 = vj(2)**2
            temp22b = -((2.d0 / fat - 4.d0 * temp21) * grad2b / temp23)
            temp22 = aux2 * vj(1) / temp23
            temp21b = temp22 * 4.d0 * grad2b / temp24
            temp18 = vj(2) * fat
            temp20 = r / vj(2)
            temp19b = 2.d0 * (1.d0 - temp20) * gradb / temp18
            aux2b = vj(1) * temp19b + vj(1) * temp22b
            vjb(1) = vjb(1) + aux * temp21b + aux2 * temp22b
            temp19 = aux2 * vj(1) / temp18
            temp20b = -(2.d0 * temp19 * gradb / vj(2))
            temp18b = -(temp19 * temp19b)
            vjb(2) = vjb(2) + fat * temp18b - temp20 * temp20b - temp22 * 2 * vj(2) * &
                    &        temp22b
            fatb = vj(2) * temp18b - temp21 * 2 * fat * temp21b + temp22 * 2.d0 * grad2b / &
                    &        fat**2
            auxb = vj(1) * temp21b
            rb = rb + temp20b
            vjb(1) = vjb(1) + aux2 * temp19b
        END IF
        vjb(1) = vjb(1) + aux * fatb
        auxb = auxb + vj(1) * fatb
        IF (ispin.gt.0) THEN
            temp9 = r / vj(2)
            temp9b0 = -(2 * (1.d0 - 0.5d0 * temp9) * 0.5d0 * auxb / vj(2))
            rb = rb + temp9b0
            vjb(2) = vjb(2) - temp9 * temp9b0
        ELSE
            temp10 = r / vj(2)
            temp10b = -(2 * (1.d0 - temp10) * auxb / vj(2))
            rb = rb + temp10b
            vjb(2) = vjb(2) - temp10 * temp10b
        END IF
        temp9b = aux2b / (4.d0 * vj(1))
        vjb(2) = vjb(2) + (vj(1) + 1) * temp9b
        vjb(1) = vjb(1) + (vj(2) - vj(2) * (vj(1) + 1) / vj(1)) * temp9b
    CASE DEFAULT
        STOP
    END SELECT
    grad2b = 0.0_8
    gradb = 0.0_8
END SUBROUTINE UPDATE_JGRAD_B
