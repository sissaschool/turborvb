!TL off
!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.2 (r3024) - 06/17/2009 13:03
!  
!  Differentiation of pseudofun in reverse (adjoint) mode:
!   gradient, with respect to input variables: r psip
!   of linear combination of output variables: r pseudofun psip
SUBROUTINE PSEUDOFUN_B(nmax, r, rb, param, psip, psipb, pseudofunb)
    IMPLICIT NONE
    REAL*8 :: pseudofun, param(3, *), r, logr, r2, psip(*)
    REAL*8 :: pseudofunb, rb, logrb, r2b, psipb(*), rsav
    INTEGER :: nmax, i
    INTEGER :: branch
    INTRINSIC DEXP
    REAL*8 :: tempb
    INTRINSIC DLOG
    ! psipb is just used here and therefore it is initialized to zero.
    !
    rsav = r

    IF (r .LT. 1d-9) THEN
        r = 1d-9
        !   CALL PUSHINTEGER4(1)
        ! ELSE
        !   CALL PUSHINTEGER4(0)
    END IF
    r2 = r**2
    logr = DLOG(r)
    DO i = 1, nmax
        psip(i) = DEXP(-(param(3, i) * r2) + logr * param(2, i))
    END DO
    pseudofun = 0.d0
    DO i = 1, nmax
        pseudofun = pseudofun + psip(i) * param(1, i)
    END DO
    r2b = -(pseudofun * pseudofunb / r2**2)
    pseudofunb = pseudofunb / r2
    DO i = nmax, 1, -1
        psipb(i) = param(1, i) * pseudofunb
    END DO
    logrb = 0.0_8
    DO i = nmax, 1, -1
        tempb = DEXP(param(2, i) * logr - param(3, i) * r2) * psipb(i)
        logrb = logrb + param(2, i) * tempb
        r2b = r2b - param(3, i) * tempb
        psipb(i) = 0.0_8
    END DO
    IF (rsav.LT. 1d-9) THEN
        r = 1d-9
        !   do nothing to rb
        !   rb = 0.0_8
    else
        rb = rb + 2 * r * r2b + logrb / r
    END IF
    pseudofunb = 0.d0

END SUBROUTINE PSEUDOFUN_B
