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

module dielectric
    real*8 epsilon0, Cgauss, dielectric_ratio, dielectric_length, vq0_diel, vgauss_diel
    integer case_diel
contains
    subroutine init_dielectric
        implicit none
        if (dielectric_ratio .eq. 1.d0) case_diel = 0
        if (dielectric_ratio .ne. 1.d0 .and. case_diel .eq. 0) case_diel = 2
        select case (case_diel)
        case (0)
            epsilon0 = 1.d0
            Cgauss = 1.d0
        case (1)
            epsilon0 = 1.d0/dielectric_ratio
            Cgauss = 0.5d0/dielectric_length**2
        case (2)
            epsilon0 = 1.d0/dielectric_ratio
            Cgauss = 1.d0/dielectric_length
        case default
            case_diel = 0
        end select
        if (dielectric_ratio .eq. -1.d0) epsilon0 = 0.d0 ! Default input no long range
    end subroutine init_dielectric
    subroutine init_dielectric_dft(kappa)
        use constants, only: Pi
        implicit none
        real*8 kappa
        select case (case_diel)
        case (0)
!  (q=0)  Int dr^3  v_shortrange (including Ewald)
            vq0_diel = pi/kappa**2
!   Int dr^3 v_total(r) x Exp[-r^2/2] (4 pi for pure Coulomb)
            vgauss_diel = 4.d0*pi
        case (1)
!  (q=0)  Int dr^3  v_shortrange (including Ewald)
            vq0_diel = epsilon0*pi/kappa**2 + 2.d0*pi*(1.d0 - epsilon0)/Cgauss
!   Int dr^3 v_total(r) x Exp[-r^2/2] (4 pi for pure Coulomb)
            vgauss_diel = 4.d0*pi*epsilon0 + (1.d0 - epsilon0)*2.d0*pi/(Cgauss + 0.5d0)
        case (2)
!  (q=0)  Int dr^3  v_shortrange (including Ewald)
            vq0_diel = epsilon0*pi/kappa**2 + (1.d0 - epsilon0)*pi/Cgauss**2
!   Int dr^3 v_total(r) x Exp[-r^2/2] (4 pi for pure Coulomb)
            vgauss_diel = 4.d0*pi*epsilon0&
           &+ (1.d0 - epsilon0)*4*pi*(1.d0 - Cgauss/sqrt(0.5d0 + Cgauss**2))
        end select
    end subroutine init_dielectric_dft
    function veps(r)
        implicit none
        real*8 r, veps, derfc
        select case (case_diel)
        case (0)
            veps = 1.d0/r
        case (1)
            veps = (epsilon0 + (1.d0 - epsilon0)*exp(-Cgauss*r*r))/r
        case (2)
            veps = (epsilon0 + (1.d0 - epsilon0)*derfc(Cgauss*r))/r
        end select
    end function veps
    subroutine veps_b(r, rb, vepsb)
!   Here rb=rb + vepsb * d/dr  veps(r)
        use constants, only: M_2_SQRTPI
        real*8 r, rb, vepsb, r2, derfc
        r2 = r*r
        select case (case_diel)
        case (0)
            rb = rb - vepsb/r2
        case (1)
            rb = rb - vepsb/r2*(epsilon0 + (1.d0 - epsilon0)*(1.d0 + 2*Cgauss*r2)*exp(-Cgauss*r2))
        case (2)
            rb = rb - epsilon0*vepsb/r2 - (1.d0 - epsilon0)*vepsb*(derfc(Cgauss*r)/r2&
           &+ M_2_SQRTPI*Cgauss*exp(-(Cgauss*r)**2)/r)
        end select
        vepsb = 0.d0
    end subroutine veps_b
    function rep_erfc(r, kappa)
        implicit none
        real*8 rep_erfc, r, kappa
        real*8 derfc
        select case (case_diel)
        case (0)
            rep_erfc = derfc(r*kappa)/r
        case (1)
            rep_erfc = epsilon0*derfc(r*kappa)/r + (1.d0 - epsilon0)*exp(-r*r*Cgauss)/r
        case (2)
            rep_erfc = epsilon0*derfc(r*kappa)/r + (1.d0 - epsilon0)*derfc(r*Cgauss)/r
        end select
    end function rep_erfc
    subroutine rep_erfc_b(r, rb, kappa, rep_erfcb)
!   Here rb=rb + rep_erfcb * d/dr  rep_erfc(r)
        use constants, only: M_2_SQRTPI
        implicit none
        real*8 r, rb, kappa, rep_erfcb, derfc, r2
        r2 = r*r
        select case (case_diel)
        case (0)
            rb = -rep_erfcb*(derfc(kappa*r)/r2 + kappa*exp(-kappa*kappa*r2)/r*M_2_SQRTPI)
        case (1)
            rb = -rep_erfcb*epsilon0*(derfc(kappa*r)/r2 + kappa*exp(-kappa*kappa*r2)/r*M_2_SQRTPI)
            rb = rb - rep_erfcb/r2*(1.d0 - epsilon0)*(1.d0 + 2*Cgauss*r2)*exp(-Cgauss*r2)
        case (2)
            rb = -rep_erfcb*epsilon0*(derfc(kappa*r)/r2 + kappa*exp(-kappa*kappa*r2)/r*M_2_SQRTPI)
            rb = rb - rep_erfcb*(1 - epsilon0)*(derfc(Cgauss*r)/r2 + Cgauss*exp(-Cgauss*Cgauss*r2)/r*M_2_SQRTPI)
        end select
!   rep_erfcb=0.d0
    end subroutine rep_erfc_b
end module dielectric

