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
!   mass   -> electron mass   9.10938291*10^-31 Kg
!   unit mass = 1/12 C_12     1.660538921*10^-27 Kg
!   time   -> 2.4188843e-17 s
!
! Instead the units used for reading and writing the result are:
!   length -> angstrom
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

    double precision, parameter :: length_unit = 0.5291772109d0
    double precision, parameter :: energy_unit = 27.211386245988d0
    double precision, parameter :: mass_unit = 5.4857990901d-4  !  m_e/unit mass (Dalton)

!   double precision, parameter :: time_unit = 2.4188843d-2  ! ???
    double precision, parameter :: press_unit = 29421.02648438959d0 ! --> Gpa
    double precision, parameter :: kboltz = 1.380649D-23 ! kg m^2/K s^2
    double precision, parameter :: machine_precision = 10d-12
    double precision, parameter :: RYDBERG = 2.179872361104D-18 !kg m^2/s^2
    double precision, PARAMETER :: eps8 = 1.0D-8
    double precision, PARAMETER :: deps = 1.0D-5
    !double precision, PARAMETER :: epsbas = -1.d0*log(1.0D-6)
    double precision, PARAMETER :: one_6 = 0.16666666666666666667d0
    double precision, PARAMETER :: one_9 = 0.11111111111111111111d0
    double precision, PARAMETER :: one_27 = 0.03703703703703703703d0
    double precision, PARAMETER :: one_54 = 0.018518518518518518516d0
    double precision, PARAMETER :: four_9 = 0.44444444444444444444d0

    !     1s orbital, STO with cusp condition
    !       double precision, parameter :: b1s=3.d0
    !       double precision, parameter :: cost1s = 0.587719
    double precision, parameter :: b1s = 1.d0
    double precision, parameter :: cost1s = 0.783091668430495d0
    !     1s orbitail, with cut P(x)=(x+a)^n/(1+(x+z)^n)
    !       (cost so that Z^(3/2) * cost * P(Z*r) * Exp(-Z*r) is normalized )
    !      Table: n, a, cost
    !       {1, 0.61803398874989484820, 0.85180469537657067484},
    !       {2, 1.00000000000000000000, 0.67579506012732275754},
    !       {3,  1.1640351402897696035, 0.61041792776122593477},
    !       {4,  1.2263393530877080588, 0.58542132302621750732},
    !       {5,  1.2466281572105588837, 0.57486261341980385597},
    !       {6,  1.2493358153441433000, 0.56997156415356984638},
    !       {7,  1.2445431960326412355, 0.56752311307928818389},
    !       {8,  1.2366315951775263509, 0.56621483934227620209}
    double precision, parameter :: costSTO1s_n = 4.0d0
    double precision, parameter :: costSTO1s_a = 1.2263393530877080588d0
    double precision, parameter :: costSTO1s_c = 0.58542132302621750732d0

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

    !     6h norm coeff
    !     cost1h=sqrt(11/pi)/16
    !     cost2h=sqrt(165/4pi)/16 x 2
    !     cost3h=sqrt(1155/4pi)/8 x 2
    !     cost4h=sqrt(385/2pi)/32.d0 x 2
    !     cost5h=3*sqrt(385.d0/4pi)/16.d0 x 2
    !     cost6h=3*sqrt(77.d0/2pi)/32.d0 x 2

    double precision, parameter :: cost1h = 0.116950322453424d0
    double precision, parameter :: cost2h = 0.452946651195697d0
    double precision, parameter :: cost3h = 2.39676839248666d0
    double precision, parameter :: cost4h = 0.489238299435250d0
    double precision, parameter :: cost5h = 2.07566231488104d0
    double precision, parameter :: cost6h = 0.656382056840170d0

!    cost1i= 1/32 sqrt[13/pi]
!    cost2i= 1/16 sqrt[273/pi]
!    cost3i= 1/64 sqrt[2x1365/pi]
!    cost4i = cost3i x 2
!    cost5i = 3/32 Sqrt[91/pi]
!    cost6i = 3/32 Sqrt[2x 1001/pi]
!    cost7i = 1/64 Sqrt[ 2 x 3003/pi]
    double precision, parameter :: cost1i = .06356920226762842593d0
    double precision, parameter :: cost2i = .58262136251873138884d0
    double precision, parameter :: cost3i = .46060262975746174957d0
    double precision, parameter :: cost4i = .92120525951492349914d0
    double precision, parameter :: cost5i = .50456490072872415925d0
    double precision, parameter :: cost6i = 2.3666191622317520320d0
    double precision, parameter :: cost7i = .68318410519191432197d0


    !       double precision, parameter :: ratios=1.01376055947381d0
    double precision, parameter :: ratios = 0.722512696367119d0
    !      double precision, parameter :: ratiop=1.03390156130503d0
    double precision, parameter :: ratiop = 1.4737345971013202d0
    !       double precision, parameter :: ratiod=1.06467059108311d0
    double precision, parameter :: ratiod = 1.752365579751248d0
    !      double precision, parameter :: ratiof=1.10812695335949d0
    double precision, parameter :: ratiof = 1.63133805835853d0
    !      double precision, parameter :: ratiog=1.16688699096108d0
    double precision, parameter :: ratiog = 1.29856658549666d0
    double precision, parameter :: ratioh = 1.24430440600566d0
    double precision, parameter :: ratioi = 1.34472087186923d0

    double precision, parameter :: ratiocs = 0.31777059276736856d0
    double precision, parameter :: ratiocp = 0.28173307010611937d0
    double precision, parameter :: ratiocd = 0.12400006759498713d0
    double precision, parameter :: ratiocf = 0.03752336975191113d0
    double precision, parameter :: ratiocg = 0.00870045684918202d0
    double precision, parameter :: ratioch = 0.00220715840582833d0
    double precision, parameter :: ratioci = 0.000579955493883660d0

#ifdef __PORT
    integer, parameter :: qp = 8
#else
    integer, parameter :: qp = 16
#endif
    integer old_threads
    integer, parameter :: dp = 8
    integer, parameter :: sp = 4
    integer  iflagerr, ip4, ipc, ipj, ipf, ip_reshuff, nbdgetri&
            &, nelorbh_ip, nmol_ip

    logical yes_ontarget
    ! steepness of the Fermi function for computing bump orbitals 8xx
    double precision, parameter :: steep_fd = 0.005d0
    !      To be defined in Initializeall
    double precision safemin, epsmach

    ! complex numbers
    complex(8), parameter :: zzero = dcmplx(0.d0, 0.d0)
    complex(8), parameter :: zone = dcmplx(1.d0, 0.d0)
    complex(8), parameter :: zonea = dcmplx(1.d0, 1.d0)
    complex(8), parameter :: zimg = dcmplx(0.d0, 1.d0)
    complex(8), parameter :: zmone = dcmplx(-1.d0, 0.d0)
    complex(8), parameter :: zhalf = dcmplx(0.5d0, 0.d0)

end module Constants
