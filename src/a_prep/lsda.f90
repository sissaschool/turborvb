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

subroutine pz_polarized(rs, ec, vc)
    !-----------------------------------------------------------------------
    !     J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
    !     spin-polarized energy and potential
    !
    ! USE kinds
    implicit none
    real*8 :: rs, ec, vc
    real*8 :: a, b, c, d, gc, b1, b2
    parameter(a=0.01555d0, b=-0.0269d0, c=0.0007d0, d= &
              -0.0048d0, gc=-0.0843d0, b1=1.3981d0, b2=0.2611d0)
    real*8 :: lnrs, rs12, ox, dox
    real*8, parameter :: xcprefact = 0.022575584d0, pi34 = 0.6203504908994d0
    ! REAL(DP) :: betha, etha, csi, prefact
    !
    if (rs .lt. 1.0d0) then
        ! high density formula
        lnrs = dlog(rs)
        ec = a*lnrs + b + c*rs*lnrs + d*rs
        vc = a*lnrs + (b - a/3.d0) + 2.d0/3.d0*c*rs*lnrs + &
             (2.d0*d - c)/3.d0*rs
    else
        ! interpolation formula
        rs12 = dsqrt(rs)
        ox = 1.d0 + b1*rs12 + b2*rs
        dox = 1.d0 + 7.d0/6.d0*b1*rs12 + 4.d0/3.d0*b2*rs
        ec = gc/ox
        vc = ec*dox/ox
    end if
    !
    !  IF ( lxc_rel ) THEN
    !     betha = prefact * pi34 / rs
    !     etha = DSQRT( 1 + betha**2 )
    !     csi = betha + etha
    !     prefact = 1.0D0 - (3.0D0/2.0D0) * ( (betha*etha - log(csi))/betha**2 )**2
    !     ec = ec * prefact
    !     vc = vc * prefact
    !  ENDIF
    return
end subroutine pz_polarized
!
!-----------------------------------------------------------------------
subroutine pz_spin(rs, zeta, ec, vcup, vcdw)
    !-----------------------------------------------------------------------
    !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
    !
    ! USE kinds
    implicit none
    real*8 :: rs, zeta, ec, vcup, vcdw
    !
    real*8 :: ecu, vcu, ecp, vcp, fz, dfz
    real*8 :: p43, third
    parameter(p43=4.0d0/3.d0, third=1.d0/3.d0)
    !
    ! unpolarized part (Perdew-Zunger formula)
    call pz(rs, 1, ecu, vcu)
    ! polarization contribution
    call pz_polarized(rs, ecp, vcp)
    !
    if (abs(zeta) .lt. 1.d0) then
        fz = ((1.0d0 + zeta)**p43 + (1.d0 - zeta)**p43 - 2.d0)/ &
             (2.d0**p43 - 2.d0)
        dfz = p43*((1.0d0 + zeta)**third - (1.d0 - zeta)**third) &
              /(2.d0**p43 - 2.d0)

        ec = ecu + fz*(ecp - ecu)
        vcup = vcu + fz*(vcp - vcu) + (ecp - ecu)*dfz*(1.d0 - zeta)
        vcdw = vcu + fz*(vcp - vcu) + (ecp - ecu)*dfz*(-1.d0 - &
                                                       zeta)
        !
    elseif (zeta .ge. 1.d0) then

        fz = 1.d0
        dfz = p43*(2.d0**third)/(2.d0**p43 - 2.d0)

        ec = ecu + fz*(ecp - ecu)
        vcup = vcu + fz*(vcp - vcu)
        vcdw = vcu + fz*(vcp - vcu) + (ecp - ecu)*dfz*(-2.d0)

    elseif (zeta .le. -1.d0) then

        fz = 1.d0
        dfz = -p43*2.d0**third/(2.d0**p43 - 2.d0)

        ec = ecu + fz*(ecp - ecu)
        vcup = vcu + fz*(vcp - vcu) + (ecp - ecu)*dfz*2.d0
        vcdw = vcu + fz*(vcp - vcu)

    end if

    !
    return
end subroutine pz_spin

subroutine pz_spinKZK(rs, zeta, ec, vcup, vcdw, vol)
    !-----------------------------------------------------------------------
    !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
    !
    ! USE kinds
    implicit none
    real*8 :: rs, zeta, ec, vcup, vcdw, vol
    !
    real*8 :: ecu, vcu, ecp, vcp, fz, dfz, mz, mz43
    real*8 :: p43, third
    parameter(p43=4.0d0/3.d0, third=1.d0/3.d0)
    !
    ! unpolarized part (Perdew-Zunger formula)
    call pzKZK(rs, ecu, vcu, vol)
    ! polarization contribution
    call pzKZK_polarized(rs, ecp, vcp, vol)
    !
    ! workaroud: if (1.d0-zeta) is negative, some compilers returns
    ! NaN when computing (1.d0-zeta)**p43
    mz = 1.d0 - zeta
    if (mz .le. 0) then
        mz43 = mz**(4.d0)
        mz43 = mz**(1.d0/3.d0)
    else
        mz43 = mz**p43
    end if
    fz = ((1.0d0 + zeta)**p43 + mz43 - 2.d0)/ &
         (2.d0**p43 - 2.d0)
    dfz = p43*((1.0d0 + zeta)**third - (1.d0 - zeta)**third) &
          /(2.d0**p43 - 2.d0)
    !
    ec = ecu + fz*(ecp - ecu)
    vcup = vcu + fz*(vcp - vcu) + (ecp - ecu)*dfz*(1.d0 - zeta)
    vcdw = vcu + fz*(vcp - vcu) + (ecp - ecu)*dfz*(-1.d0 - &
                                                   zeta)
    !
    return
end subroutine pz_spinKZK

subroutine pzKZK_polarized(rs, ec, vc, vol)
    !-----------------------------------------------------------------------
    !     LDA parameterization form Monte Carlo data
    !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
    !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
    !
    ! It is assumed that the finite size correction of the correlation energy
    ! are the same for the polarized and unpolarized electron gas at the same
    ! density
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
    ! real*8 :: a, b, c, d, gc, b1, b2
    ! parameter (a = 0.01555d0, b = - 0.0269d0, c = 0.0007d0, d = &
    !      - 0.0048d0, gc = - 0.0843d0, b1 = 1.3981d0, b2 = 0.2611d0)
    ! real*8 :: lnrs, rs12, ox, dox
    ! REAL*8, PARAMETER :: xcprefact = 0.022575584d0, pi34 = 0.6203504908994d0

    data a/0.01555d0, 0.031091d0/, b/-0.0269d0, -0.046644d0/, &
        c/0.0007d0, 0.00419d0/, d/-0.0048d0, -0.00983d0/
    data gc/-0.0843d0, -0.103756d0/, b1/1.3981d0, 0.56371d0/, &
        b2/0.2611d0, 0.27358d0/
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
end subroutine pzKZK_polarized
