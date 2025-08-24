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

!> @brief      Calculate diffusion coefficient for real wave function
!> @details    This subroutine computes the diffusion coefficient for the
!>             quantum Monte Carlo walkers. It evaluates the quantity that
!>             should be exactly one for perfect short-time diffusion.
!> @param[in]  nel      Number of electrons
!> @param[in]  indt     Number of time steps
!> @param[in]  gamma    Diffusion parameter
!> @param[in]  ivic     Velocity array (3, indt, nel)
!> @param[in]  table    Wave function table
!> @param[out] fun      Diffusion coefficient
!> @param[in]  istart   Starting index for calculation
subroutine diffus(nel, indt, gamma, ivic, table, fun, istart)

    implicit none

    integer nel, indt, i, j, istart
    real(8) ivic(3, indt, *), table(nel, *), fun, amu2, gamma, cost

    !        this quantity should be exactly one for perfect
    !        short time diffusion
    fun = 0.d0
    do i = 1, istart - 1
        cost = 0.d0
        do j = 1, nel
            amu2 = sum(ivic(:, i, j)**2)
            !        this regularize the diffusion in the true H^a
            cost = cost + table(j, i)*amu2
        end do
        fun = fun + cost
    end do
    fun = fun/(6.d0*nel) ! This has to be one for perfect diffusion

    return
end subroutine diffus

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!

!> @brief      Calculate diffusion coefficient for complex wave function
!> @details    This subroutine computes the diffusion coefficient for complex
!>             wave functions in quantum Monte Carlo. Similar to diffus but
!>             handles complex wave function tables.
!> @param[in]  nel      Number of electrons
!> @param[in]  indt     Number of time steps
!> @param[in]  gamma    Diffusion parameter
!> @param[in]  ivic     Velocity array (3, indt, nel)
!> @param[in]  table    Complex wave function table
!> @param[out] fun      Complex diffusion coefficient
!> @param[in]  istart   Starting index for calculation
subroutine diffus_complex(nel, indt, gamma, ivic, table, fun, istart)

    use Constants, only: zzero

    implicit none

    integer nel, indt, i, j, istart
    real(8) ivic(3, indt, *), amu2, gamma
    complex(8) table(nel, *), fun, cost

    !        this quantity should be exactly one for perfect
    !        short time diffusion
    fun = zzero
    do i = 1, istart - 1
        cost = zzero
        do j = 1, nel
            amu2 = sum(ivic(:, i, j)**2)
            !        this regularize the diffusion in the true H^a
            cost = cost + table(j, i)*amu2
        end do
        fun = fun + cost
    end do
    fun = fun/(6.d0*nel) ! This has to be one for perfect diffusion

    return
end subroutine diffus_complex
