! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2002-2004 quantum-ESPRESSO group
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

subroutine slater(rs, ex, vx)
    !-----------------------------------------------------------------------
    !        Slater exchange with alpha=2/3
    !
    implicit none
    real*8 :: rs, ex, vx
    !  atom units:  Hartree and Bohr radius.
    real*8, parameter :: f = -0.687247939924714d0, alpha = 2.0d0/3.0d0
    ! f = -9/8*(3/2pi)^(2/3)
    !
    ex = f*alpha/rs
    vx = 4.d0/3.d0*f*alpha/rs
    !
    return
end subroutine slater

subroutine pz(rs, iflagr, ec, vc)
    !-----------------------------------------------------------------------
    !     LDA parameterization from Monte Carlo data
    !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
    !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
    !
    implicit none
    real*8 :: rs, ec, vc
    integer :: iflagr, iflag
    !
    real*8 :: a(2), b(2), c(2), d(2), gc(2), b1(2), b2(2), corr(2)
    real*8 :: lnrs, rs12, ox, dox
    !
    data a/0.0311d0, 0.031091d0/, b/-0.048d0, -0.046644d0/, &
        c/0.0020d0, 0.00419d0/, d/-0.0116d0, -0.00983d0/
    data gc/-0.1423d0, -0.103756d0/, b1/1.0529d0, 0.56371d0/, &
        b2/0.3334d0, 0.27358d0/, corr/1.d0, 1.d0/
    !

    iflag = abs(iflagr)
    ! corrected intepolation formula to allow continuity of ec at rs=1.

    if (iflagr .lt. 0) then
        corr(iflag) = gc(iflag)/(b(iflag) + d(iflag)) - b1(iflag) - b2(iflag)
    else
        corr(iflag) = 1.d0
    end if

    if (rs .lt. 1.0d0) then
        ! high density formula
        lnrs = dlog(rs)
        ec = a(iflag)*lnrs + b(iflag) + c(iflag)*rs*lnrs + d( &
             iflag)*rs
        vc = a(iflag)*lnrs + (b(iflag) - a(iflag)/3.d0) + 2.d0/ &
             3.d0*c(iflag)*rs*lnrs + (2.d0*d(iflag) - c(iflag)) &
             /3.d0*rs
    else
        ! interpolation formula
        rs12 = dsqrt(rs)
        ox = corr(iflag) + b1(iflag)*rs12 + b2(iflag)*rs
        dox = corr(iflag) + 7.d0/6.d0*b1(iflag)*rs12 + 4.d0/3.d0* &
              b2(iflag)*rs
        ec = gc(iflag)/ox
        vc = ec*dox/ox
    end if
    !
    return
end subroutine pz
!

!-----------------------------------------------------------------------
subroutine slaterKZK(rs, ex, vx, vol, iopt)
    !-----------------------------------------------------------------------
    !        Slater exchange with alpha=2/3, Kwee, Zhang and Krakauer KE
    !        correction
    !
    ! USE kinds
    implicit none
    real*8 :: rs, ex, vx, dL, vol, ga, a0, a3dl5_true, ex_t
    real*8, parameter :: a1 = -2.2037d0, &
                         !  value of a3 published a3=-0.015 e-mail Zhang 29/12/2008 more accurate.
                         !             a2 = 0.4710d0, a3 = -0.01504001153d0, ry2h = 0.5d0
                         a2 = 0.4710d0, a3 = -1.504002727495861d-2, ry2h = 0.5d0, a4 = -2.530982213157182d-2
    real*8, parameter :: f = -0.687247939924714d0, alpha = 2.0d0/3.0d0, &
            & pi = 3.14159265358979323846d0
    integer :: iopt
    ! f = -9/8*(3/2pi)^(2/3)
    !
    !  vx = ex- 1/3 d ex/dlog rs
    !  pi = 4.d0 * atan(1.d0)
    !  a4= (ex-vx/3.d0)*ga**6/dl**5  for rs-->ga^-
    a0 = f*alpha*2.d0

    dL = vol**(1.d0/3.d0)
    ga = 0.5d0*dL*(3.d0/pi)**(1.d0/3.d0)
    !
    if (rs .le. ga) then
        ex = a0/rs + a1*rs/dL**2.d0 + a2*rs**2.d0/dL**3.d0
        vx = (4.d0*a0/rs + 2.d0*a1*rs/dL**2.d0 + &
              a2*rs**2.d0/dL**3.d0)/3.d0

        if (iopt .lt. 0) ex = ex - a4*dl**5/ga**9*rs**3 ! also the exchange is continuous

    else

        if (iopt .gt. 0) then
            ex = a0/ga + a1*ga/dL**2.d0 + a2*ga**2.d0/dL**3.d0 ! solids
            vx = ex

        else

            ex = a3*dL**5.d0/rs**6.d0 ! molecules & solids
            vx = 3.d0*ex

        end if

    end if

    ex = ry2h*ex ! Ry to Hartree
    vx = ry2h*vx
    !
    return
end subroutine slaterKZK

!-----------------------------------------------------------------------
subroutine pzKZK(rs, ec, vc, vol)
    !-----------------------------------------------------------------------
    !     LDA parameterization form Monte Carlo data
    !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
    !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
    !
    ! USE kinds
    implicit none
    real*8 :: rs, ec, vc, ec0(2), vc0(2), ec0p
    integer :: iflag, kr
    !
    real*8 :: a(2), b(2), c(2), d(2), gc(2), b1(2), b2(2)
    real*8 :: lnrs, rs12, ox, dox, lnrsk, rsk
    real*8 :: a1, grs, g1, g2, g3, g4, dL, vol, gh, gl, grsp
    real*8 :: f3, f2, f1, f0, pi
    real*8 :: D1, D2, D3, P1, P2, ry2h
    !
    data a/0.0311d0, 0.031091d0/, b/-0.048d0, -0.046644d0/, &
        c/0.0020d0, 0.00419d0/, d/-0.0116d0, -0.00983d0/
    data gc/-0.1423d0, -0.103756d0/, b1/1.0529d0, 0.56371d0/, &
        b2/0.3334d0, 0.27358d0/
    data a1/-2.2037d0/, g1/0.1182d0/, g2/1.1656d0/, g3/-5.2884d0/, &
        g4/-1.1233d0/
    data ry2h/0.5d0/
    !
    iflag = 1
    pi = 4.d0*atan(1.d0)
    dL = vol**(1.d0/3.d0)
    gh = 0.5d0*dL/(2.d0*pi)**(1.d0/3.d0)
    gl = dL*(3.d0/2.d0/pi)**(1.d0/3.d0)

    rsk = gh
    do kr = 1, 2
        lnrsk = dlog(rsk)
        if (rsk .lt. 1.0d0) then
            ! high density formula
            ec0(kr) = a(iflag)*lnrsk + b(iflag) + c(iflag)*rsk*lnrsk + d( &
                      iflag)*rsk
            vc0(kr) = a(iflag)*lnrsk + (b(iflag) - a(iflag)/3.d0) + 2.d0/ &
                      3.d0*c(iflag)*rsk*lnrsk + (2.d0*d(iflag) - c(iflag)) &
                      /3.d0*rsk
        else
            ! interpolation formula
            rs12 = dsqrt(rsk)
            ox = 1.d0 + b1(iflag)*rs12 + b2(iflag)*rsk
            dox = 1.d0 + 7.d0/6.d0*b1(iflag)*rs12 + 4.d0/3.d0* &
                  b2(iflag)*rsk
            ec0(kr) = gc(iflag)/ox
            vc0(kr) = ec0(kr)*dox/ox
        end if
        !
        grs = g1*rsk*lnrsk + g2*rsk + g3*rsk**1.5d0 + g4*rsk**2.d0
        grsp = g1*lnrsk + g1 + g2 + 1.5d0*g3*rsk**0.5d0 + &
               2.d0*g4*rsk
        ec0(kr) = ec0(kr) + (-a1*rsk/dL**2.d0 + grs/dL**3.d0)*ry2h
        vc0(kr) = vc0(kr) + (-2.d0*a1*rsk/dL**2.d0/3.d0 + &
                             grs/dL**3.d0 - grsp*rsk/3.d0/dL**3.d0)*ry2h
        !
        rsk = rs
    end do

    lnrs = dlog(rs)
    if (rs .le. gh) then
        ec = ec0(2)
        vc = vc0(2)
    else
        if (rs .le. gl) then
            ec0p = 3.d0*(ec0(1) - vc0(1))/gh
            P1 = 3.d0*ec0(1) - gh*ec0p
            P2 = ec0p
            D1 = gl - gh
            D2 = gl**2.d0 - gh**2.d0
            D3 = gl**3.d0 - gh**3.d0
            f2 = 2.d0*gl**2.d0*P2*D1 + D2*P1
            f2 = f2/(-(2.d0*gl*D1)**2.d0 + 4.d0*gl*D1*D2 - D2**2.d0)
            f3 = -(P2 + 2.d0*D1*f2)/(3.d0*D2)
            f1 = -(P1 + D2*f2)/(2.d0*D1)
            f0 = -gl*(gl*f2 + 2.d0*f1)/3.d0
            !
            ec = f3*rs**3.d0 + f2*rs**2.d0 + f1*rs + f0
            vc = f2*rs**2.d0/3.d0 + f1*2.d0*rs/3.d0 + f0
        else
            ec = 0.d0
            vc = 0.d0
        end if
    end if
    !
    return
end subroutine pzKZK
!
!---------------------------------------------------------------
subroutine pbex(rho, grho, iflagr, sx, v1x, v2x)
    !---------------------------------------------------------------
    !
    ! PBE exchange (without Slater exchange):
    ! iflagr=4  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
    ! iflagr=5  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)
    !
    implicit none
    real*8 :: rho, grho, sx, v1x, v2x
    ! input: charge and squared gradient
    ! output: energy
    ! output: potential
    integer :: iflag, iflagr
    ! local variables
    real*8 :: kf, agrho, s1, s2, ds, dsg, exunif, fx
    ! (3*pi2*|rho|)^(1/3)
    ! |grho|
    ! |grho|/(2*kf*|rho|)
    ! s^2
    ! n*ds/dn
    ! n*ds/d(gn)
    ! exchange energy LDA part
    ! exchange energy gradient part
    real*8 :: dxunif, dfx, f1, f2, f3, dfx1
    ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
    real*8 :: third, c1, c2, c5, pi
    parameter(pi=3.14159265358979323846d0)
    parameter(third=1.d0/3.d0, c1=0.75d0/pi, &
              c2=3.093667726280136d0, c5=4.d0*third)
    ! parameters of the functional
    real*8 :: k(2), mu
    data k/0.804d0, 1.2450d0/, mu/0.21951d0/
    !
    iflag = iflagr - 9

    agrho = dsqrt(grho)
    kf = c2*rho**third
    dsg = 0.5d0/kf
    s1 = agrho*dsg/rho
    s2 = s1*s1
    ds = -c5*s1
    !
    !   Energy
    !
    f1 = s2*mu/k(iflag)
    f2 = 1.d0 + f1
    f3 = k(iflag)/f2
    fx = k(iflag) - f3
    exunif = -c1*kf
    sx = exunif*fx
    !
    !   Potential
    !
    dxunif = exunif*third
    dfx1 = f2*f2
    dfx = 2.d0*mu*s1/dfx1
    v1x = sx + dxunif*fx + exunif*dfx*ds
    v2x = exunif*dfx*dsg/agrho

    sx = sx*rho
    return
end subroutine pbex
!
!---------------------------------------------------------------
subroutine pbec(rho, grho, sc, v1c, v2c)
    !---------------------------------------------------------------
    !
    ! PBE correlation (without LDA part)
    ! J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
    !
    implicit none
    real*8 :: rho, grho, sc, v1c, v2c
    real*8 :: ga, be
    parameter(ga=0.031091d0, be=0.066725d0)
    real*8 :: third, pi34, xkf, xks
    parameter(third=1.d0/3.d0, pi34=0.6203504908994d0)
    parameter(xkf=1.919158292677513d0, xks=1.128379167095513d0)
    ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
    real*8 :: kf, ks, rs, ec, vc, t, expe, af, bf, y, xy, qy
    real*8 :: s1, h0, dh0, ddh0
    !
    rs = pi34/rho**third
    call pw(rs, 1, ec, vc)
    kf = xkf/rs
    ks = xks*dsqrt(kf)
    t = dsqrt(grho)/(2.d0*ks*rho)
    expe = dexp(-ec/ga)
    af = be/ga*(1.d0/(expe - 1.d0))
    bf = expe*(vc - ec)
    y = af*t*t
    xy = (1.d0 + y)/(1.d0 + y + y*y)
    qy = y*y*(2.d0 + y)/(1.d0 + y + y*y)**2
    s1 = 1.d0 + be/ga*t*t*xy
    h0 = ga*dlog(s1)
    dh0 = be*t*t/s1*(-7.d0/3.d0*xy - qy*(af*bf/ &
                                         be - 7.d0/3.d0))
    ddh0 = be/(2.d0*ks*ks*rho)*(xy - qy)/s1
    sc = rho*h0
    v1c = h0 + dh0
    v2c = ddh0
    !
    return
end subroutine pbec

!-----------------------------------------------------------------------
subroutine pw(rs, iflag, ec, vc)
    !-----------------------------------------------------------------------
    !     iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
    !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
    !
    implicit none
    real*8 :: rs, ec, vc
    integer :: iflag
    !
    real*8 :: a, b1, b2, c0, c1, c2, c3, d0, d1
    parameter(a=0.031091d0, b1=7.5957d0, b2=3.5876d0, c0=a, &
              c1=0.046644d0, c2=0.00664d0, c3=0.01043d0, d0=0.4335d0, &
              d1=1.4408d0)
    real*8 :: lnrs, rs12, rs32, rs2, om, dom, olog
    real*8 :: a1(2), b3(2), b4(2)
    data a1/0.21370d0, 0.026481d0/, b3/1.6382d0, -0.46647d0/, &
        b4/0.49294d0, 0.13354d0/
    !
    ! high- and low-density formulae implemented but not used in PW case
    ! (reason: inconsistencies in PBE/PW91 functionals)
    !

    if (rs .lt. 1d0 .and. iflag .eq. 2) then
        ! high density formula
        lnrs = dlog(rs)
        ec = c0*lnrs - c1 + c2*rs*lnrs - c3*rs
        vc = c0*lnrs - (c1 + c0/3.d0) + 2.d0/3.d0*c2*rs* &
             lnrs - (2.d0*c3 + c2)/3.d0*rs
    elseif (rs .gt. 100.d0 .and. iflag .eq. 2) then
        ! low density formula
        ec = -d0/rs + d1/rs**1.5d0
        vc = -4.d0/3.d0*d0/rs + 1.5d0*d1/rs**1.5d0
    else
        ! interpolation formula
        rs12 = dsqrt(rs)
        rs32 = rs*rs12
        rs2 = rs**2
        om = 2.d0*a*(b1*rs12 + b2*rs + b3(iflag)*rs32 + b4( &
                     iflag)*rs2)
        dom = 2.d0*a*(0.5d0*b1*rs12 + b2*rs + 1.5d0*b3( &
                      iflag)*rs32 + 2.d0*b4(iflag)*rs2)
        olog = dlog(1.d0 + 1.0d0/om)
        ec = -2.d0*a*(1.d0 + a1(iflag)*rs)*olog
        vc = -2.d0*a*(1.d0 + 2.d0/3.d0*a1(iflag)*rs) &
             *olog - 2.d0/3.d0*a*(1.d0 + a1(iflag)*rs)*dom/ &
             (om*(om + 1.d0))
    end if
    !
    return
end subroutine pw
!
