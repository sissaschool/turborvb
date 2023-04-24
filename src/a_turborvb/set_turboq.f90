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

subroutine set_turboq(gamma_eig, kdyn_eig, dt, alphaall, alphaqmc, mnoise&
        &, Gn, Tn, Gni, Gnh, Gnih)
    implicit none
    real*8, intent(in) :: gamma_eig, kdyn_eig, alphaall, alphaqmc, dt
    real*8, intent(out) :: mnoise(2, 2), Gn, Tn, Gni(2, 2), Gnh, Gnih
    real*8 cost, costm1, xbar, sinhx, coshx, sinx, cosx, gammag, gammagx, gammagxm
    real*8 dergamma, der2gamma, costh, sinhxh, coshxh, costhm1
    real*8, external :: mygamma
    complex*16 gammagxc
    complex*16, external :: mygammac
    logical xbarpos

    if (kdyn_eig .gt. 1d-6) then
        cost = exp(-dt*gamma_eig/2.d0)
        costh = exp(-dt*gamma_eig/4.d0)
        xbar = (gamma_eig/2.d0)**2 - kdyn_eig
        if (xbar .ge. 0) then
            xbar = sqrt(xbar)
            coshx = cosh(xbar*dt)
            sinhx = sinh(xbar*dt)
            coshxh = cosh(xbar*dt/2.d0)
            sinhxh = sinh(xbar*dt/2.d0)
            xbarpos = .true.
        else
            xbar = sqrt(-xbar)
            coshx = cos(xbar*dt)
            sinhx = sin(xbar*dt)
            coshxh = cos(xbar*dt/2.d0)
            sinhxh = sin(xbar*dt/2.d0)
            xbarpos = .false.
        end if

        if (xbar*dt .gt. 1d-6) then
            Gn = cost/xbar*sinhx
            Gnh = costh/xbar*sinhxh
            Tn = (1 - cost*coshx)/kdyn_eig - 0.5d0*cost*gamma_eig*sinhx/(xbar*kdyn_eig)
            Gni(1, 1) = cost*(coshx - 0.5d0*gamma_eig/xbar*sinhx)
            Gni(1, 2) = -cost*kdyn_eig*sinhx/xbar
            Gni(2, 1) = Gn
            Gni(2, 2) = cost*(coshx + 0.5d0*gamma_eig/xbar*sinhx) - 1.d0 ! -1 as we remove the identity
            Gnih = costh*(coshxh + 0.5d0*gamma_eig/xbar*sinhxh) - 1.d0
        else
            Gn = cost*dt
            Gnh = costh*dt/2.d0
            Tn = (1 - cost*coshx)/kdyn_eig - 0.5d0*cost*gamma_eig*dt/kdyn_eig
            Gni(1, 1) = cost*(coshx - 0.5d0*gamma_eig*dt)
            Gni(1, 2) = -cost*kdyn_eig*dt
            Gni(2, 1) = Gn
            Gni(2, 2) = cost*(coshx + 0.5d0*gamma_eig*dt) - 1.d0 ! -1 as we remove the identity
            Gnih = costh*(coshxh + 0.25d0*gamma_eig*dt) - 1.d0
        end if

        !          Protection from roundoff
        if (xbar*dt .le. 1d-6) then
            xbar = 1d-6/dt
        end if
        if (xbarpos) then
            gammag = mygamma(gamma_eig, dt)
            gammagx = mygamma(gamma_eig + 2*xbar, dt)
            gammagxm = mygamma(gamma_eig - 2*xbar, dt)
            dergamma = (gammagx - gammagxm)/xbar
            der2gamma = (gammagx + gammagxm - 2*gammag)/xbar**2
            mnoise(1, 1) = alphaall*0.25d0*(gammagx + gammagxm + 2*gammag) + &
                    & 0.0625d0*alphaall*gamma_eig**2*der2gamma + &
                    & 0.25d0*alphaall*gamma_eig*dergamma
            mnoise(2, 2) = 0.25d0*alphaall*der2gamma
            mnoise(1, 2) = -0.125d0*alphaall*gamma_eig*der2gamma - &
                    & 0.25d0*alphaall*dergamma
        else
            gammag = mygamma(gamma_eig, dt)
            gammagxc = mygammac(dcmplx(gamma_eig, 2*xbar), dt)
            dergamma = 2*aimag(gammagxc)/xbar
            der2gamma = (2*gammagxc - 2*gammag)/xbar**2
            mnoise(1, 1) = alphaall*0.25d0*(2*gammagxc + 2*gammag) - &
                    & 0.0625d0*alphaall*gamma_eig**2*der2gamma + &
                    & 0.25d0*alphaall*gamma_eig*dergamma
            mnoise(2, 2) = -0.25d0*alphaall*der2gamma
            mnoise(1, 2) = 0.125d0*alphaall*gamma_eig*der2gamma - &
                    & 0.25d0*alphaall*dergamma
        end if
        !          Now the bar noise to be added to the force with noise
        !          correction

        mnoise(1, 1) = mnoise(1, 1)/Gn**2 - alphaqmc
        mnoise(2, 2) = mnoise(2, 2)/Tn**2 - alphaqmc
        mnoise(1, 2) = mnoise(1, 2)/Tn/Gn - alphaqmc
        mnoise(2, 1) = mnoise(1, 2)

    else

        costh = exp(-dt*gamma_eig/2.d0)
        if (costh .lt. 0.99999d0) then
            costhm1 = costh - 1.d0
        else
            costhm1 = -dt/2.d0*gamma_eig + (dt*gamma_eig/2.d0)**2/2.d0
        end if
        cost = exp(-dt*gamma_eig)
        if (cost .lt. 0.99999d0) then
            costm1 = cost - 1.d0
        else
            costm1 = -dt*gamma_eig + (dt*gamma_eig)**2/2.d0
        end if

        Gn = -costm1/gamma_eig
        Gnh = -costhm1/gamma_eig
        Gnih = 0.d0
        Tn = 1.d0/gamma_eig*(dt - Gn)
        Gni(1, 1) = cost
        Gni(2, 1) = Gn
        Gni(2, 2) = 0.d0
        Gni(1, 2) = 0.d0
        !        Tn=1.d0/gamma_eig**2*(dt*gamma_eig+costm1)
        mnoise(1, 1) = -0.5d0*alphaall/gamma_eig*(1.d0 + cost)*costm1

        !         mnoise(1,1)=-0.5d0*alphaall*gamma_eig*(1.d0+cost)/costm1-alphaqmc

        mnoise(2, 2) = alphaall/gamma_eig**2* &
                &(dt + 1.d0/gamma_eig*costm1*(2.d0 - 0.5d0*(1.d0 + cost)))

        !         mnoise(1,2)=alphaall/gamma_eig*&
        !    & (-costm1/gamma_eig*(1.-0.5d0*(1+cost)))

        mnoise(1, 2) = 0.5d0*alphaall*(costm1/gamma_eig)**2

        mnoise(1, 1) = mnoise(1, 1)/Gn**2 - alphaqmc
        mnoise(2, 2) = mnoise(2, 2)/Tn**2 - alphaqmc
        mnoise(1, 2) = mnoise(1, 2)/Tn/Gn - alphaqmc
        mnoise(2, 1) = mnoise(1, 2)

    end if

    return
end
function mygamma(x, dt)
    implicit none
    real*8 x, mygamma, dt, xdt
    xdt = x*dt
    if (abs(xdt) .gt. 1d-6) then
        mygamma = (1.d0 - exp(-xdt))/x
    else
        mygamma = dt - 0.5d0*xdt*dt
    end if
    return
end
function mygammac(x, dt)
    implicit none
    complex*16 x, xdt, mygammac
    real*8 dt
    xdt = x*dt
    if (abs(xdt) .gt. 1d-6) then
        mygammac = (1.d0 - exp(-xdt))/x
    else
        mygammac = dt - 0.5d0*xdt*dt
    end if
    return
end

