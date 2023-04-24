!TL off
!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.3 (r3163) - 09/25/2009 09:03
!
!  Differentiation of jastrow_ee in reverse (adjoint) mode:
!   gradient, with respect to input variables: rc vj
!   of linear combination of output variables: rc vj jastrow_ee
SUBROUTINE JASTROW_EE_B(rc, rcb, vj, vjb, iesd, ispin, jastrow_eeb, metricb)
    use cell, only : metric
    use allio, only : iespbc, norm_metric
    IMPLICIT NONE
    REAL*8 :: r, rc(*), rz, jastrow_ee, vj(*)
    REAL*8 :: rb, fatb, fat, rcb(*), rzb, jastrow_eeb, vjb(*)
    INTEGER :: iesd, ispin
    REAL*8 :: aux, aux2, aux3, temp4b2, temp4b7, temp4b3, auxb, rs, rsb
    REAL*8 :: temp3
    REAL*8 :: temp2
    REAL*8 :: temp1
    REAL*8 :: temp0
    REAL*8 :: temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14
    REAL*8 :: temp6b, temp7b, temp10b, temp11b, temp12b, temp13b, temp14b
    INTRINSIC DEXP
    REAL*8 :: tempb9
    REAL*8 :: tempb8
    REAL*8 :: tempb7
    REAL*8 :: tempb6
    REAL*8 :: tempb5
    REAL*8 :: tempb4
    REAL*8 :: tempb3
    REAL*8 :: tempb2
    REAL*8 :: tempb1
    REAL*8 :: tempb0
    REAL*8 :: temp0b
    REAL*8 :: tempb12
    REAL*8 :: tempb11
    REAL*8 :: tempb10
    REAL*8 :: temp3b
    REAL*8 :: temp2b2
    REAL*8 :: temp2b1
    REAL*8 :: temp2b0
    REAL*8 :: temp5b2
    REAL*8 :: temp5b1
    REAL*8 :: temp5b0
    REAL*8 :: tempb
    REAL*8 :: temp0b1
    REAL*8 :: temp0b0
    REAL*8 :: temp2b
    REAL*8 :: temp5b
    REAL*8 :: temp1b
    INTRINSIC DSQRT
    REAL*8 :: temp
    REAL*8 :: temp1b4
    REAL*8 :: temp1b3
    REAL*8 :: temp1b2
    REAL*8 :: temp1b1
    REAL*8 :: temp1b0
    REAL*8 :: temp4b
    REAL*8 :: temp4
    real*8 metricb(3, 3)
    integer i, j

    if(iespbc) then
        r = norm_metric(rc, metric)
    else
        r = dsqrt(sum(rc(1:3)**2))
    endif

    SELECT CASE  (iesd)
    CASE (4)
        ! exponential form
        IF (ispin .LT. 0) THEN
            tempb0 = 0.5d0 * jastrow_eeb / vj(1)
            tempb1 = -(DEXP(-(vj(1) * r)) * tempb0)
            vjb(1) = vjb(1) - (1.d0 - DEXP(-(vj(1) * r))) * tempb0 / vj(1) - r * tempb1
            rb = -(vj(1) * tempb1)
        ELSE
            tempb2 = 0.25d0 * jastrow_eeb / vj(1)
            tempb3 = -(DEXP(-(vj(1) * r)) * tempb2)
            vjb(1) = vjb(1) - (1.d0 - DEXP(-(vj(1) * r))) * tempb2 / vj(1) - r * tempb3
            rb = -(vj(1) * tempb3)
        END IF
    CASE (-4)
        ! exponential form
        tempb0 = 0.5d0 * jastrow_eeb / vj(1)
        tempb1 = -(DEXP(-(vj(1) * r)) * tempb0)
        vjb(1) = vjb(1) - (1.d0 - DEXP(-(vj(1) * r))) * tempb0 / vj(1) - r * tempb1
        rb = -(vj(1) * tempb1)
    CASE (1)
        tempb4 = 0.5d0 * jastrow_eeb / (vj(1) * r + 1.d0)
        tempb5 = -(r * tempb4 / (vj(1) * r + 1.d0))
        rb = vj(1) * tempb5 + tempb4
        vjb(1) = vjb(1) + r * tempb5
    CASE (6)
        IF (vj(2) .GT. 1d-9) THEN
            rs = (1.d0 - DEXP(-(vj(2) * r))) / vj(2)
        ELSE
            rs = r
        END IF

        aux = 1.d0 / (1.d0 + vj(1) * rs)
        aux2 = aux * aux
        aux3 = aux2 * aux

        temp4 = DEXP(-(vj(2) * r))
        auxb = 0.5d0 * rs * jastrow_eeb
        temp4b7 = -(auxb / (vj(1) * rs + 1.d0)**2)
        rsb = vj(1) * temp4b7 + 0.5d0 * aux * jastrow_eeb
        vjb(1) = vjb(1) + rs * temp4b7

        IF (vj(2).LT.1d-9) THEN
            rb = rsb
        ELSE
            temp4b2 = rsb / vj(2)
            temp4b3 = -(DEXP(-(vj(2) * r)) * temp4b2)
            vjb(2) = vjb(2) - (1.d0 - DEXP(-(vj(2) * r))) * temp4b2 / vj(2) - r * &
                    &        temp4b3
            rb = - vj(2) * temp4b3
        END IF

    CASE (7)
        IF (ispin .GT. 0) THEN
            tempb8 = 0.25d0 * jastrow_eeb / (vj(1) * r + 1.d0)
            tempb9 = -(r * tempb8 / (vj(1) * r + 1.d0))
            rb = vj(1) * tempb9 + tempb8
            vjb(1) = vjb(1) + r * tempb9
        ELSE
            tempb10 = 0.5d0 * jastrow_eeb / (vj(1) * r + 1.d0)
            tempb11 = -(r * tempb10 / (vj(1) * r + 1.d0))
            rb = vj(1) * tempb11 + tempb10
            vjb(1) = vjb(1) + r * tempb11
        END IF
    CASE (9)
        temp3 = 4.d0 * vj(1)
        temp2 = vj(2) * (vj(1) + 1.d0)
        temp = temp2 / temp3
        temp0 = r / vj(2)
        temp1 = vj(1) * (-temp0 + 1)**2 + 1
        temp1b = -(temp * jastrow_eeb / temp1)
        temp0b = -(vj(1) * 2 * (1 - temp0) * temp1b / vj(2))
        tempb0 = -(LOG(temp1) * jastrow_eeb / temp3)
        vjb(1) = vjb(1) + (vj(2) - temp * 4.d0) * tempb0 + (1 - temp0)**2 * temp1b
        rb = temp0b
        vjb(2) = vjb(2) + (vj(1) + 1.d0) * tempb0 - temp0 * temp0b
    CASE (10)
        ! The same with cusp condition for equal spin electrons r-> r/2
        IF (ispin .GT. 0) THEN
            temp9 = 4.d0 * vj(1)
            temp8 = vj(2) * (vj(1) + 1.d0)
            temp4 = temp8 / temp9
            temp7 = 2.d0 * vj(2)
            temp5 = r / temp7
            temp6 = vj(1) * (-temp5 + 1)**2 + 1
            temp6b = -(temp4 * jastrow_eeb / temp6)
            temp5b = -(vj(1) * 2 * (1 - temp5) * temp6b / temp7)
            temp4b = -(LOG(temp6) * jastrow_eeb / temp9)
            vjb(1) = vjb(1) + (vj(2) - temp4 * 4.d0) * temp4b + (1 - temp5)**2 * temp6b
            rb = temp5b
            vjb(2) = vjb(2) + (vj(1) + 1.d0) * temp4b - temp5 * 2.d0 * temp5b
        ELSE
            temp14 = 4.d0 * vj(1)
            temp13 = vj(2) * (vj(1) + 1.d0)
            temp10 = temp13 / temp14
            temp11 = r / vj(2)
            temp12 = vj(1) * (-temp11 + 1)**2 + 1
            temp12b = -(temp10 * jastrow_eeb / temp12)
            temp11b = -(vj(1) * 2 * (1 - temp11) * temp12b / vj(2))
            temp10b = -(LOG(temp12) * jastrow_eeb / temp14)
            vjb(1) = vjb(1) + (vj(2) - temp10 * 4.d0) * temp10b + (1 - temp11)**2 * &
                    &        temp12b
            rb = temp11b
            vjb(2) = vjb(2) + (vj(1) + 1.d0) * temp10b - temp11 * temp11b
        END IF
    CASE (5)
        temp = vj(1) + vj(3) * vj(2)
        temp0 = DEXP(-(vj(2) * r))
        temp0b = 0.5d0 * jastrow_eeb / temp
        temp0b0 = -(DEXP(-(vj(1) * r)) * temp0b)
        temp0b1 = -(vj(3) * DEXP(-(vj(2) * r)) * temp0b)
        tempb12 = -((vj(3) - DEXP(-(vj(1) * r)) - vj(3) * temp0 + 1.d0) * temp0b / temp)
        vjb(3) = vjb(3) + vj(2) * tempb12 + (1.0_8 - temp0) * temp0b
        vjb(1) = vjb(1) + tempb12 - r * temp0b0
        rb = -(vj(2) * temp0b1) - vj(1) * temp0b0
        vjb(2) = vjb(2) + vj(3) * tempb12 - r * temp0b1
    CASE (8)
        !     Jastrow for pseudo soft
        !        direct algorithm intermediate variables needs to be recomputed.
        !        r=dsqrt(rc(1)**2+rc(2)**2+rc(3)**2)
        fat = dexp(-vj(1) * r**3)
        !        jastrow=1.d0/vj(1)*(1.d0-fat)
        !         cost is local so costb needs to be initialized to zero.
        ! adjoint of jastrow=1.d0/vj(1)*(1.d0-fat)
        vjb(1) = vjb(1) - jastrow_eeb / vj(1)**2 * (1.d0 - fat)
        fatb = -jastrow_eeb / vj(1)

        ! adjoint of fat=dexp(-vj(1)*r**3)
        vjb(1) = vjb(1) - fatb * fat * r**3
        rb = -fatb * fat * vj(1) * 3 * r**2
        ! adjoint of dsqrt(rc(1)**2+rc(2)**2+rc(3)**2), check on division by zero.
        !        if(r.gt.1d-9) then
        !        rcb(1)=rcb(1)+rb*rc(1)/r
        !        rcb(2)=rcb(2)+rb*rc(2)/r
        !        rcb(3)=rcb(3)+rb*rc(3)/r
        !        endif

    CASE (-1)
        ! spin contaminated: up-up and up-down cusps
        ! 1 parameter equal to -7
        ! parallel  spins
        IF (ispin .GT. 0) THEN
            temp1b = 0.25d0 * jastrow_eeb / (vj(1) * r + 1.d0)
            temp1b0 = -(r * temp1b / (vj(1) * r + 1.d0))
            rb = vj(1) * temp1b0 + temp1b
            vjb(1) = vjb(1) + r * temp1b0
        ELSE
            temp1b1 = 0.5d0 * jastrow_eeb / (vj(1) * r + 1.d0)
            temp1b2 = -(r * temp1b1 / (vj(1) * r + 1.d0))
            rb = vj(1) * temp1b2 + temp1b1
            vjb(1) = vjb(1) + r * temp1b2
        END IF
    CASE (-2)
        ! asymmetric jastrow_ee (xy ne z component)
        rz = DSQRT(vj(1)**2 * (rc(1)**2 + rc(2)**2) + (vj(2) * rc(3))**2)
        temp2b = 0.5d0 * jastrow_eeb / (rz + 1.d0)
        rb = temp2b
        rzb = -(r * temp2b / (rz + 1.d0))
        temp1 = rc(1)**2 + rc(2)**2
        IF (vj(1)**2 * temp1 + (vj(2) * rc(3))**2 .EQ. 0.0) THEN
            temp2b0 = 0.0
        ELSE
            temp2b0 = rzb / (2.D0 * DSQRT(vj(1)**2 * temp1 + (vj(2) * rc(3))**2))
        END IF
        temp1b3 = vj(1)**2 * temp2b0
        temp1b4 = 2 * vj(2) * rc(3) * temp2b0
        vjb(1) = vjb(1) + temp1 * 2 * vj(1) * temp2b0
        rcb(1) = rcb(1) + 2 * rc(1) * temp1b3
        rcb(2) = rcb(2) + 2 * rc(2) * temp1b3
        vjb(2) = vjb(2) + rc(3) * temp1b4
        rcb(3) = rcb(3) + vj(2) * temp1b4
    CASE (3)
        temp4 = (vj(2) * r + 1.d0)**2
        temp2 = vj(3) * r / temp4
        temp3 = vj(1) * r + 1.d0
        temp4b = 0.5d0 * r * jastrow_eeb
        temp3b = -(temp4b / temp3**2)
        temp2b1 = temp4b / temp4
        temp2b2 = -(temp2 * 2 * (vj(2) * r + 1.d0) * temp2b1)
        rb = vj(2) * temp2b2 + vj(3) * temp2b1 + vj(1) * temp3b + 0.5d0 * (1.0 / temp3&
                & + temp2) * jastrow_eeb
        vjb(1) = vjb(1) + r * temp3b
        vjb(3) = vjb(3) + r * temp2b1
        vjb(2) = vjb(2) + r * temp2b2
    CASE (2)
        ! spin contaminated: up-up and up-down cusps
        ! 2 parameters
        ! parallel  spins
        IF (ispin .GT. 0) THEN
            temp5b = 0.25d0 * jastrow_eeb / (vj(2) * r + 1.d0)
            temp5b0 = -(r * temp5b / (vj(2) * r + 1.d0))
            rb = vj(2) * temp5b0 + temp5b
            vjb(2) = vjb(2) + r * temp5b0
        ELSE
            temp5b1 = 0.5d0 * jastrow_eeb / (vj(1) * r + 1.d0)
            temp5b2 = -(r * temp5b1 / (vj(1) * r + 1.d0))
            rb = vj(1) * temp5b2 + temp5b1
            vjb(1) = vjb(1) + r * temp5b2
        END IF
    CASE (0)
        rb = 0.0_8
    END SELECT
    IF(r.ne.0.d0) THEN
        if(iespbc) then
            call norm_metric_b(rc, rcb, metric, metricb, rb)
        else
            tempb = rb / r
            rcb(1) = rcb(1) + rc(1) * tempb
            rcb(2) = rcb(2) + rc(2) * tempb
            rcb(3) = rcb(3) + rc(3) * tempb
        endif
    endif
    jastrow_eeb = 0.d0
END SUBROUTINE JASTROW_EE_B
