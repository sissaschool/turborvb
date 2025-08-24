! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2004-2007 QUANTUM-ESPRESSO group
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

!> @file atom_weight.f90
!> @brief Module containing atomic weight lookup functionality
!> @author TurboRVB group
!> @date 2022
!> @version 1.0

!> @brief Returns the atomic weight (relative atomic mass) of an element given its atomic number
!> @details This function provides access to a lookup table of atomic weights for elements with 
!>          atomic numbers 1-103 (Hydrogen to Lawrencium). The atomic weights are based on the 
!>          standard atomic weights recommended by IUPAC, taking into account the natural abundance 
!>          of isotopes.
!> 
!>          The function performs bounds checking and returns an error if an invalid atomic number 
!>          is provided.
!> 
!>          Atomic weights are dimensionless quantities relative to 12C = 12, stored with double 
!>          precision (DP) for accuracy. The data array is indexed from 1 to 103, matching atomic numbers.
!> 
!> @param[in] atomic_number The atomic number Z of the element (1-103)
!> @return The atomic weight in atomic mass units (u)
!> 
!> @note If atomic_number < 1 or atomic_number > 103, the function calls errore() with an error 
!>       code of 1000 + atomic_number.
!> 
!> @warning The error message in the original code mentions 'atom_name' but this function is 
!>          'atom_weight'. This might be a copy-paste error from another function.
!> 
!> @par Examples:
!> @code
!> atom_weight(1)  ! returns 1.00794 (Hydrogen)
!> atom_weight(6)  ! returns 12.0107 (Carbon)  
!> atom_weight(79) ! returns 196.96655 (Gold)
!> @endcode
!> 
!> @par Data Source:
!> IUPAC standard atomic weights (natural abundance)
!> 
!> @see kinds module for DP precision definition
!> @see errore() for error handling
function atom_weight(atomic_number)
    ! ------------------------------------------------------------------
    !
    use kinds, only: DP
    implicit none
    
    !> @param[in] atomic_number The atomic number Z of the element (1-103)
    integer, intent(in) :: atomic_number
    
    !> @return The atomic weight in atomic mass units
    real(DP) :: atom_weight
    
    !> @var weights Array storing atomic weights for elements 1-103
    !> @details Index corresponds to atomic number Z. Values are based on IUPAC standard 
    !>          atomic weights (natural abundance).
    real(DP) :: weights(103)
    
    ! Atomic weight data array (in atomic mass units)
    ! Index corresponds to atomic number Z
    ! Values are based on IUPAC standard atomic weights (natural abundance)
    data weights/1.00794_dp, 4.00260_dp, &                    !< H, He
        6.941_dp, 9.01218_dp, 10.811_dp, 12.0107_dp, 14.00674_dp, &  !< Li, Be, B, C, N
        15.9994_dp, 18.99840_dp, 20.1797_dp, &                !< O, F, Ne
        22.98977_dp, 24.3050_dp, 26.98154_dp, 28.0855_dp, 30.97376_dp, &  !< Na, Mg, Al, Si, P
        32.066_dp, 35.4527_dp, 39.948_dp, &                   !< S, Cl, Ar
        39.0983_dp, 40.078_dp, 44.95591_dp, 47.867_dp, 50.9415_dp, &  !< K, Ca, Sc, Ti, V
        51.9961_dp, 54.93805_dp, 55.845_dp, &                 !< Cr, Mn, Fe
        58.93320_dp, 58.6934_dp, 63.546_dp, 65.39_dp, &       !< Co, Ni, Cu, Zn
        69.723_dp, 72.61_dp, 74.92160_dp, 78.96_dp, 79.904_dp, 83.80_dp, &  !< Ga, Ge, As, Se, Br, Kr
        85.4678_dp, 87.62_dp, 88.90585_dp, 91.224_dp, 92.90638_dp, &  !< Rb, Sr, Y, Zr, Nb
        95.94_dp, 98._dp, &                                    !< Mo, Tc
        101.07_dp, 102.90550_dp, 106.42_dp, 107.8682_dp, 112.411_dp, &  !< Ru, Rh, Pd, Ag, Cd
        114.818_dp, 118.710_dp, 121.760_dp, 127.60_dp, 126.90447_dp, &  !< In, Sn, Sb, Te, I
        131.29_dp, &                                           !< Xe
        132.90545_dp, 137.327_dp, 138.9055_dp, 140.116_dp, 140.90765_dp, &  !< Cs, Ba, La, Ce, Pr
        144.24_dp, 145._dp, 150.36_dp, 151.964_dp, 157.25_dp, &  !< Nd, Pm, Sm, Eu, Gd
        158.92534_dp, 162.50_dp, 164.93032_dp, 167.26_dp, &    !< Tb, Dy, Ho, Er
        168.93421_dp, 173.04_dp, 174.967_dp, &                 !< Tm, Yb, Lu
        178.49_dp, 180.9479_dp, 183.84_dp, 186.207_dp, 190.23_dp, &  !< Hf, Ta, W, Re, Os
        192.217_dp, 195.078_dp, 196.96655_dp, 200.59_dp, &     !< Ir, Pt, Au, Hg
        204.3833_dp, 207.2_dp, 208.98038_dp, 209._dp, 210._dp, 222._dp, &  !< Tl, Pb, Bi, Po, At, Rn
        223._dp, 226._dp, 227._dp, 232.0381_dp, 231.03588_dp, &  !< Fr, Ra, Ac, Th, Pa
        238.0289_dp, 237._dp, 244._dp, &                       !< U, Np, Pu
        243._dp, 247._dp, 247._dp, 251._dp, 252._dp, 257._dp, &  !< Am, Cm, Bk, Cf, Es, Fm
        258._dp, 259._dp, 262._dp/                             !< Md, No, Lr

    ! Input validation: check if atomic number is within valid range
    if (atomic_number < 1 .or. atomic_number > 103) then
        ! Error: atomic number out of range (1-103)
        ! Note: The error message mentions 'atom_name' but this is 'atom_weight'
        ! This might be a copy-paste error from another function
        call errore('atom_weight', 'invalid atomic number', 1000 + atomic_number)
    else
        ! Return the atomic weight for the given atomic number
        ! Array indexing: weights(1) = Hydrogen, weights(6) = Carbon, etc.
        atom_weight = weights(atomic_number)
    end if
    
    return

end function atom_weight
!
