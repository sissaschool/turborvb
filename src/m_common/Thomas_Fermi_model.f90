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

module Thomas_Fermi_model
    ! Author Kosuke Nakano 2019
    ! Affiliation SISSA
    ! e-mail:kousuke_1123@icloud.com
    ! created 13 Dec. 2019
    ! modified 13 Feb. 2020. The accuracy of the integration has been improved.

    implicit none

contains

    ! function Thomas_Fermi_core_electron_number(atomic_number, r_c)
    ! to calculate the number of electrons inside r_c according to the Thomas-Fermi model.
    ! According to the textbook written in Takada (it can be also obtained from the Landau textbook),
    ! The number of electron inside rc is represented as
    !     b = (9 * np.pi ** 2 / 128) ** (1.0 / 3.0)
    !     r_TF = b * Z ** (-1.0 / 3.0)
    !     x_TF = r / r_TF
    !     rho_r = (32 * Z ** 2) / (9 * np.pi ** 3) * (Gross_Dreizler(x_TF) / x_TF) ** (3.0 / 2.0)
    !     p_r = exp(-r** 2.0 / 2.0 / rc ** 2.0)          # function p(r) in the double-grid scheme
    !     num_ele_core = integral (4 * np.pi * r ** 2 * rho_r * p_r)     * dr
    !     num_ele_val  = integral (4 * np.pi * r ** 2 * rho_r * (1-p_r)) * dr

    function Thomas_Fermi_core_electron_number(atomic_number, r_c) result(num_core_electron)

        real(8), intent(in) :: atomic_number, r_c
        real(8) num_core_electron, num_val_electron
        real(8) r_TF, x_TF, b, pi, rho_r, p_r
        real(8) :: r = 0.0d0, r_min = 1.0d-3, r_max = 25.0d0, dr = 1.0d-3

        ! define constant
        pi = 3.1415926535d0
        b = (9*pi**2/128)**(1.0d0/3.0d0)
        r_TF = b*atomic_number**(-1.0/3.0)

        ! results
        num_core_electron = 0.0d0
        num_val_electron = 0.0d0

        ! integral
        r = r_min
        do while (r .lt. r_max)
            x_TF = r/r_TF
            rho_r = (32.0d0*atomic_number**2.0d0)/(9.0d0*pi**3.0d0)*(Gross_Dreizler(x_TF)/x_TF)**(3.0d0/2.0d0)
            p_r = dexp(-r**2.0d0/2.0d0/r_c**2.0d0)
            num_core_electron = num_core_electron + (4.0d0*pi*r**2.0d0*rho_r*p_r)*dr
            num_val_electron = num_val_electron + (4.0d0*pi*r**2.0d0*rho_r*(1.0d0 - p_r))*dr
            r = r + dr
            ! write(6,*) "r = ", r
        end do

        ! output result (test)
        ! write(6,*) "Z = ", atomic_number
        ! write(6,*) "r_c = ", r_c
        ! write(6,*) "num_core_electron = ", num_core_electron
        ! write(6,*) "num_val_electron = ", num_val_electron

        ! end

    end function Thomas_Fermi_core_electron_number

    ! Gross Dreizler developed a good approximation of chi(x)

    function Gross_Dreizler(x) result(f)

        real(8), intent(in) :: x
        real(8) f

        f = 1.0d0/(1.0d0 + 1.4712d0*x - 0.4973d0*x**(3.0d0/2.0d0) + 0.3875d0*x**(2.0d0) + 0.002102d0*x**(3.0d0))

    end function Gross_Dreizler

end module Thomas_Fermi_model
