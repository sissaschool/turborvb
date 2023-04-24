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

module atom_names
    implicit none
    character(2), dimension(:), allocatable :: AtomsNames

contains

    subroutine load_names
        implicit none
        allocate (AtomsNames(103))
        AtomsNames(1) = "H"
        AtomsNames(2) = "He"
        AtomsNames(3) = "Li"
        AtomsNames(4) = "Be"
        AtomsNames(5) = "B"
        AtomsNames(6) = "C"
        AtomsNames(7) = "N"
        AtomsNames(8) = "O"
        AtomsNames(9) = "F"
        AtomsNames(10) = "Ne"
        AtomsNames(11) = "Na"
        AtomsNames(12) = "Mg"
        AtomsNames(13) = "Al"
        AtomsNames(14) = "Si"
        AtomsNames(15) = "P"
        AtomsNames(16) = "S"
        AtomsNames(17) = "Cl"
        AtomsNames(18) = "Ar"
        AtomsNames(19) = "K"
        AtomsNames(20) = "Ca"
        AtomsNames(21) = "Sc"
        AtomsNames(22) = "Ti"
        AtomsNames(23) = "V"
        AtomsNames(24) = "Cr"
        AtomsNames(25) = "Mn"
        AtomsNames(26) = "Fe"
        AtomsNames(27) = "Co"
        AtomsNames(28) = "Ni"
        AtomsNames(29) = "Cu"
        AtomsNames(30) = "Zn"
        AtomsNames(31) = "Ga"
        AtomsNames(32) = "Ge"
        AtomsNames(33) = "As"
        AtomsNames(34) = "Se"
        AtomsNames(35) = "Br"
        AtomsNames(36) = "Kr"
        AtomsNames(37) = "Rb"
        AtomsNames(38) = "Sr"
        AtomsNames(39) = "Y"
        AtomsNames(40) = "Zr"
        AtomsNames(41) = "Nb"
        AtomsNames(42) = "Mo"
        AtomsNames(43) = "Tc"
        AtomsNames(44) = "Ru"
        AtomsNames(45) = "Rh"
        AtomsNames(46) = "Pd"
        AtomsNames(47) = "Ag"
        AtomsNames(48) = "Cd"
        AtomsNames(49) = "In"
        AtomsNames(50) = "Sn"
        AtomsNames(51) = "Sb"
        AtomsNames(52) = "Te"
        AtomsNames(53) = "I"
        AtomsNames(54) = "Xe"
        AtomsNames(55) = "Cs"
        AtomsNames(56) = "Ba"
        AtomsNames(57) = "La"
        AtomsNames(58) = "Ce"
        AtomsNames(59) = "Pr"
        AtomsNames(60) = "Nd"
        AtomsNames(61) = "Pm"
        AtomsNames(62) = "Sm"
        AtomsNames(63) = "Eu"
        AtomsNames(64) = "Gd"
        AtomsNames(65) = "Tb"
        AtomsNames(66) = "Dy"
        AtomsNames(67) = "Ho"
        AtomsNames(68) = "Er"
        AtomsNames(69) = "Tm"
        AtomsNames(70) = "Yb"
        AtomsNames(71) = "Lu"
        AtomsNames(72) = "Hf"
        AtomsNames(73) = "Ta"
        AtomsNames(74) = "W"
        AtomsNames(75) = "Re"
        AtomsNames(76) = "Os"
        AtomsNames(77) = "Ir"
        AtomsNames(78) = "Pt"
        AtomsNames(79) = "Au"
        AtomsNames(80) = "Hg"
        AtomsNames(81) = "Ti"
        AtomsNames(82) = "Pb"
        AtomsNames(83) = "Bi"
        AtomsNames(84) = "Po"
        AtomsNames(85) = "At"
        AtomsNames(86) = "Rn"
        AtomsNames(87) = "Fr"
        AtomsNames(88) = "Ra"
        AtomsNames(89) = "Ac"
        AtomsNames(90) = "Th"
        AtomsNames(91) = "Pa"
        AtomsNames(92) = "U"
        AtomsNames(93) = "Np"
        AtomsNames(94) = "Pu"
        AtomsNames(95) = "Am"
        AtomsNames(96) = "Cm"
        AtomsNames(97) = "Bk"
        AtomsNames(98) = "Cf"
        AtomsNames(99) = "Es"
        AtomsNames(100) = "Fm"
        AtomsNames(101) = "101"
        AtomsNames(102) = "No"
        AtomsNames(103) = "Lr"
    end subroutine load_names

    function get_number(symbol)
        implicit none
        character(2), intent(in) :: symbol
        integer i, get_number

        get_number = -1
        do i = 1, 103
            if (trim(symbol) == AtomsNames(i)) get_number = i
        end do
        return
    end function get_number

end module atom_names
