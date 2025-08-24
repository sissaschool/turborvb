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

!> @brief Dielectric screening module for quantum Monte Carlo calculations
!>
!> This module provides functionality for handling dielectric screening effects
!> in quantum Monte Carlo calculations. It implements various dielectric models
!> for modifying Coulomb interactions in condensed matter systems, particularly
!> useful for simulating systems with different dielectric environments.
!>
!> The module supports three different dielectric models:
!> - case_diel = 0: Pure Coulomb interaction (no screening)
!> - case_diel = 1: Gaussian screening model
!> - case_diel = 2: Error function screening model
!>
!> @details
!> The module contains global variables for dielectric parameters and provides
!> functions for calculating screened potentials and their derivatives.
!> It is used in quantum Monte Carlo calculations to account for dielectric
!> screening effects in bulk materials, surfaces, and interfaces.
!>
!> @note Used in VMC and DMC calculations for systems with dielectric screening
!> @note Supports different screening models for various physical systems
module dielectric
    real*8 epsilon0, Cgauss, dielectric_ratio, dielectric_length, vq0_diel, vgauss_diel
    integer case_diel
contains

    !> @brief Initialize dielectric parameters based on screening model
    !>
    !> This subroutine initializes the dielectric parameters based on the
    !> specified dielectric ratio and screening model. It sets up the
    !> appropriate constants for the chosen dielectric screening approach.
    !>
    !> @details
    !> The subroutine handles three cases:
    !> - case_diel = 0: Pure Coulomb (ε₀ = 1, C_gauss = 1)
    !> - case_diel = 1: Gaussian screening (ε₀ = 1/ratio, C_gauss = 0.5/length²)
    !> - case_diel = 2: Error function screening (ε₀ = 1/ratio, C_gauss = 1/length)
    !>
    !> @note If dielectric_ratio = -1, ε₀ is set to 0 (no long-range interaction)
    !> @note Automatically determines case_diel if not explicitly set
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

    !> @brief Initialize dielectric parameters for DFT calculations
    !>
    !> This subroutine initializes dielectric parameters specifically for
    !> density functional theory (DFT) calculations, computing the q=0
    !> components and Gaussian integrals needed for DFT implementations.
    !>
    !> @param[in] kappa Screening parameter for Ewald summation
    !>
    !> @details
    !> The subroutine calculates:
    !> - vq0_diel: q=0 component of the screened potential
    !> - vgauss_diel: Gaussian integral of the total potential
    !>
    !> Different formulas are used for each dielectric case:
    !> - case 0: Pure Coulomb interaction
    !> - case 1: Gaussian screening model
    !> - case 2: Error function screening model
    !>
    !> @note Used in DFT calculations with dielectric screening
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

    !> @brief Calculate screened Coulomb potential
    !>
    !> This function calculates the screened Coulomb potential at a given
    !> distance r, using the appropriate dielectric screening model.
    !>
    !> @param[in] r Distance between particles (atomic units)
    !> @return Screened Coulomb potential V(r)
    !>
    !> @details
    !> The function implements different screening models:
    !> - case 0: V(r) = 1/r (pure Coulomb)
    !> - case 1: V(r) = (ε₀ + (1-ε₀)exp(-C_gauss*r²))/r (Gaussian screening)
    !> - case 2: V(r) = (ε₀ + (1-ε₀)erfc(C_gauss*r))/r (Error function screening)
    !>
    !> @note Used in quantum Monte Carlo calculations for screened interactions
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

    !> @brief Calculate screened Coulomb potential with backward propagation
    !>
    !> This subroutine calculates the screened Coulomb potential and its
    !> derivative for backward propagation in automatic differentiation.
    !> It modifies the input variables rb and vepsb according to the
    !> chain rule for derivatives.
    !>
    !> @param[in] r Distance between particles (atomic units)
    !> @param[in,out] rb Backward propagated variable for r
    !> @param[in,out] vepsb Backward propagated variable for potential
    !>
    !> @details
    !> The subroutine implements the backward propagation rule:
    !> rb = rb + vepsb * dV(r)/dr
    !> where dV(r)/dr is the derivative of the screened potential.
    !>
    !> Different formulas are used for each dielectric case:
    !> - case 0: Pure Coulomb derivative
    !> - case 1: Gaussian screening derivative
    !> - case 2: Error function screening derivative
    !>
    !> @note Used in automatic differentiation for quantum Monte Carlo
    !> @note vepsb is set to zero after use (chain rule completion)
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

    !> @brief Calculate screened complementary error function potential
    !>
    !> This function calculates the screened complementary error function
    !> potential, which is used in Ewald summation for long-range interactions
    !> with dielectric screening.
    !>
    !> @param[in] r Distance between particles (atomic units)
    !> @param[in] kappa Screening parameter for Ewald summation
    !> @return Screened complementary error function potential
    !>
    !> @details
    !> The function implements different screening models:
    !> - case 0: erfc(κr)/r (pure Coulomb with Ewald screening)
    !> - case 1: ε₀*erfc(κr)/r + (1-ε₀)*exp(-r²*C_gauss)/r (Gaussian screening)
    !> - case 2: ε₀*erfc(κr)/r + (1-ε₀)*erfc(C_gauss*r)/r (Error function screening)
    !>
    !> @note Used in Ewald summation for periodic systems with dielectric screening
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

    !> @brief Calculate screened complementary error function potential with backward propagation
    !>
    !> This subroutine calculates the screened complementary error function
    !> potential and its derivative for backward propagation in automatic
    !> differentiation, specifically for Ewald summation with dielectric screening.
    !>
    !> @param[in] r Distance between particles (atomic units)
    !> @param[in,out] rb Backward propagated variable for r
    !> @param[in] kappa Screening parameter for Ewald summation
    !> @param[in,out] rep_erfcb Backward propagated variable for potential
    !>
    !> @details
    !> The subroutine implements the backward propagation rule:
    !> rb = rb + rep_erfcb * dV(r)/dr
    !> where dV(r)/dr is the derivative of the screened complementary error function potential.
    !>
    !> Different formulas are used for each dielectric case, incorporating
    !> both the Ewald screening (κ) and dielectric screening (C_gauss) parameters.
    !>
    !> @note Used in automatic differentiation for Ewald summation with dielectric screening
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

