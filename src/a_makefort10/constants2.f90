! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2002-2006 quantum-ESPRESSO group
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

!----------------------------------------------------------------------------
module constants2
    !----------------------------------------------------------------------------
    !
    use kinds, only: DP
    !
    ! ... The constants needed everywhere
    !
    implicit none
    !
    save
    !
    ! ... Mathematical constants
    !
    real(DP), parameter :: pi = 3.14159265358979323846_dp
    real(DP), parameter :: tpi = 2.0_dp*pi
    real(DP), parameter :: fpi = 4.0_dp*pi
    real(DP), parameter :: sqrtpi = 1.77245385090551602729_dp
    real(DP), parameter :: sqrtpm1 = 1.0_dp/sqrtpi
    real(DP), parameter :: sqrt2 = 1.41421356237309504880_dp
    !
    ! ... Physical constants, SI (NIST CODATA 2006), Web Version 5.1
    !     http://physics.nist.gov/constants
    real(DP), parameter :: H_PLANCK_SI = 6.62606896e-34_dp ! J s
    real(DP), parameter :: K_BOLTZMANN_SI = 1.3806504e-23_dp ! J K^-1
    real(DP), parameter :: ELECTRON_SI = 1.602176487e-19_dp ! C
    real(DP), parameter :: ELECTRONVOLT_SI = 1.602176487e-19_dp ! J
    real(DP), parameter :: ELECTRONMASS_SI = 9.10938215e-31_dp ! Kg
    real(DP), parameter :: HARTREE_SI = 4.35974394e-18_dp ! J
    real(DP), parameter :: RYDBERG_SI = HARTREE_SI/2.0_dp ! J
    real(DP), parameter :: BOHR_RADIUS_SI = 0.52917720859e-10_dp ! m
    real(DP), parameter :: AMU_SI = 1.660538782e-27_dp ! Kg
    real(DP), parameter :: C_SI = 2.99792458e+8_dp ! m sec^-1
    !
    ! ... Physical constants, atomic units:
    ! ... AU for "Hartree" atomic units (e = m = hbar = 1)
    ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
    !
    real(DP), parameter :: K_BOLTZMANN_AU = K_BOLTZMANN_SI/HARTREE_SI
    real(DP), parameter :: K_BOLTZMANN_RY = K_BOLTZMANN_SI/RYDBERG_SI
    !
    ! ... Unit conversion factors: energy and masses
    !
    real(DP), parameter :: AUTOEV = HARTREE_SI/ELECTRONVOLT_SI
    real(DP), parameter :: RYTOEV = AUTOEV/2.0_dp
    real(DP), parameter :: AMU_AU = AMU_SI/ELECTRONMASS_SI
    real(DP), parameter :: AMU_RY = AMU_AU/2.0_dp
    !
    ! ... Unit conversion factors: atomic unit of time, in s and ps
    !
    real(DP), parameter :: AU_SEC = H_PLANCK_SI/tpi/HARTREE_SI
    real(DP), parameter :: AU_PS = AU_SEC*1.0e+12_dp
    !
    ! ... Unit conversion factors: pressure (1 Pa = 1 J/m^3, 1GPa = 10 Kbar )
    !
    real(DP), parameter :: AU_GPA = HARTREE_SI/BOHR_RADIUS_SI**3 &
                           /1.0e+9_dp
    real(DP), parameter :: RY_KBAR = 10.0_dp*AU_GPA/2.0_dp
    !
    ! ... Unit conversion factors: 1 debye = 10^-18 esu*cm
    ! ...                                  = 3.3356409519*10^-30 C*m
    ! ...                                  = 0.208194346 e*A
    ! ... ( 1 esu = (0.1/c) Am, c=299792458 m/s)
    !
    real(DP), parameter :: DEBYE_SI = 3.3356409519_dp*1.0e-30_dp ! C*m
    real(DP), parameter :: AU_DEBYE = ELECTRON_SI*BOHR_RADIUS_SI/ &
                           DEBYE_SI
    !
    real(DP), parameter :: eV_to_kelvin = ELECTRONVOLT_SI/K_BOLTZMANN_SI
    real(DP), parameter :: ry_to_kelvin = RYDBERG_SI/K_BOLTZMANN_SI
    !
    !  Speed of light in atomic units
    !
    real(DP), parameter :: C_AU = C_SI/BOHR_RADIUS_SI*AU_SEC
    !
    ! ... zero up to a given accuracy
    !
    real(DP), parameter :: eps4 = 1.0e-4_dp
    real(DP), parameter :: eps6 = 1.0e-6_dp
    real(DP), parameter :: eps8 = 1.0e-8_dp
    real(DP), parameter :: eps12 = 1.0e-12_dp
    real(DP), parameter :: eps14 = 1.0e-14_dp
    real(DP), parameter :: eps16 = 1.0e-16_dp
    real(DP), parameter :: eps24 = 1.0e-24_dp
    real(DP), parameter :: eps32 = 1.0e-32_dp
    !
    real(DP), parameter :: gsmall = 1.0e-12_dp
    !
    real(DP), parameter :: e2 = 2.0_dp ! the square of the electron charge
    real(DP), parameter :: degspin = 2.0_dp ! the number of spins per level
    !
    !!!!!! COMPATIBIILITY
    !
    real(DP), parameter :: amconv = AMU_RY
    real(DP), parameter :: uakbar = RY_KBAR
    real(DP), parameter :: bohr_radius_cm = bohr_radius_si*100.0_dp
    real(DP), parameter :: BOHR_RADIUS_ANGS = bohr_radius_cm*1.0e8_dp
    real(DP), parameter :: ANGSTROM_AU = 1.0_dp/BOHR_RADIUS_ANGS
    real(DP), parameter :: DIP_DEBYE = AU_DEBYE
    real(DP), parameter :: AU_TERAHERTZ = AU_PS
    real(DP), parameter :: AU_TO_OHMCMM1 = 46000.0_dp ! (ohm cm)^-1
    !

end module constants2

! perl script to create a program to list the available constants:
! extract with: grep '^!XX!' constants.f90 | sed 's,!XX!,,' > mkconstlist.pl
! then run: perl mkconstlist.pl constants.f90 > testme.f90
! and compile and run: testme.f90
!XX!#!/usr/bin/perl -w
!XX!
!XX!use strict;
!XX!
!XX!print <<EOF
!XX!! list all available constants and derived values
!XX!
!XX!PROGRAM list_constants
!XX!
!XX!  USE kinds, ONLY : DP
!XX!  USE constants
!XX!
!XX!EOF
!XX!;
!XX!
!XX!while(<>) {
!XX!  if ( /REAL\s*\(DP\)\s*,\s*PARAMETER\s*::\s*([a-zA-Z_0-9]+)\s*=.*$/ )  {
!XX!    print "  WRITE (*,'(A18,G24.17)') '$1:',$1\n";
!XX!  }
!XX!}
!XX!
!XX!print <<EOF
!XX!
!XX!END PROGRAM list_constants
!XX!EOF
!XX!;
!XX!
