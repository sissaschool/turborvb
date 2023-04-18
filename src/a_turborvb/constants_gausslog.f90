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

! This module contains the definition of some physical constants

! Note about measure units: the program works internally with atomic units:
!   lenght -> bohr radius   (= 0.52917 angstrom)
!   energy -> hartree       (= 27.211383 eV)
!   mass   -> electron mass (proton mass = 1836.1527 electron mass)
!   time   -> 2.4188843e-17 s
!
! Instead the units used for reading and writing the result are:
!   lenght -> angstrom
!   energy -> eV
!   mass   -> amu (= 1/12 Carbon mass)
!   time   -> fs  (1e-15 s)
!   press  -> kbar

module Constants
    implicit none

    double precision, parameter :: PI = 3.14159265358979323846d0
    double precision, parameter :: TWO_PI = 6.28318530717958647692d0
    double precision, parameter :: PI2 = 9.86960440108935861881d0 !pi^2
    double precision, parameter :: M_E = 2.7182818284590452354d0 !e
    double precision, parameter :: LOG2E = 1.4426950408889634074d0 !log_2 e
    double precision, parameter :: PI_2 = 1.57079632679489661923d0 !pi/2
    double precision, parameter :: PI_4 = 0.78539816339744830962d0 !pi/4
    double precision, parameter :: INV_PI = 0.31830988618379067154d0 !1/pi
    double precision, parameter :: INV_PI2 = 0.10132118364233777144d0 !1/pi^2
    double precision, parameter :: INV2_PI = 0.63661977236758134308d0 !2/pi
    double precision, parameter :: SQRT2 = 1.41421356237309504880d0 !sqrt(2)
    double precision, parameter :: INV_SQRT2 = 0.70710678118654752440d0 !1/sqrt(2)
    double precision, parameter :: M_2_SQRTPI = 1.1283791670955125738961589031215452d0 ! 2/sqrt(pi)

    double precision, parameter :: lenght_unit = 0.52917d0
    double precision, parameter :: energy_unit = 27.211383d0
    double precision, parameter :: mass_unit = 5.4461702d-4
    double precision, parameter :: time_unit = 2.4188843d-2
    double precision, parameter :: press_unit = 294210.1d0
    double precision, parameter :: kboltz = 1.3806503D-23 ! kg m^2/K s^2
    double precision, parameter :: machine_precision = 10d-12
    double precision, parameter :: RYDBERG = 2.17987190389E-18 !kg m^2/s^2

    !     3d without cusp condition
    !     cost1d=0.5 cost2d=dsqrt(3.d0)*cost1d cost3d=2.d0*cost2d
    double precision, parameter :: cost1d = 0.5d0
    double precision, parameter :: cost2d = 0.866025403784439d0
    double precision, parameter :: cost3d = 1.73205080756888d0

    !     4f norm coeff
    !     cost1f=0.5 cost2f=dsqrt(6.d0)/2.d0*cost1f
    !     cost3f=dsqrt(15.d0)*cost1f
    !     cost4f=dsqrt(10.d0)/2.d0*cost1f
    double precision, parameter :: cost1f = 0.5d0
    double precision, parameter :: cost2f = 0.612372435695794d0
    double precision, parameter :: cost3f = 1.93649167310371d0
    double precision, parameter :: cost4f = 0.790569415042095d0

    !     5g norm coeff
    !     cost1g=1.d0/8.d0
    !     cost2g=sqrt(5.d0/2.d0)/2.d0
    !     cost3g=sqrt(5.d0)/4.d0
    !     cost4g=sqrt(35.d0/2.d0)/2.d0
    !     cost5g=sqrt(35.d0)/8.d0
    double precision, parameter :: cost1g = 0.125d0
    double precision, parameter :: cost2g = 0.79056941504209d0
    double precision, parameter :: cost3g = 0.55901699437494d0
    double precision, parameter :: cost4g = 2.09165006633518d0
    double precision, parameter :: cost5g = 0.73950997288745d0

    double precision, parameter :: ratios = 1.01376055947381d0
    double precision, parameter :: ratiop = 1.03390156130503d0
    double precision, parameter :: ratiod = 1.06467059108311d0
    double precision, parameter :: ratiof = 1.10812695335949d0
    double precision, parameter :: ratiog = 1.16688699096108d0
    double precision, parameter :: ratioh = 1.24430440600566d0
    double precision, parameter :: ratioi = 1.34472087186923d0

    double precision, parameter :: ratiocs = 0.445865236040760d0
    double precision, parameter :: ratiocp = 0.197650419299989d0
    double precision, parameter :: ratiocd = 0.0753377187877891d0
    double precision, parameter :: ratiocf = 0.0254886822445654d0
    double precision, parameter :: ratiocg = 0.00781819740791014d0
    double precision, parameter :: ratioch = 0.00220715840582833d0
    double precision, parameter :: ratioci = 0.000579955493883660d0

    integer, parameter :: dp = 8
    integer, parameter :: sp = 4
    integer  iflagerr, ip4, nbdgetri

    !
    ! STRING SIZE
    !
    !      integer, parameter  :: lchlen=300

end module Constants
