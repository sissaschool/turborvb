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

module atom_label
    implicit none
    character(20), dimension(:), allocatable :: AtomsLabels

contains

    subroutine load_label
        implicit none
        allocate (AtomsLabels(103))
        AtomsLabels(1) = "hydrogen"
        AtomsLabels(2) = "helium"
        AtomsLabels(3) = "lithium"
        AtomsLabels(4) = "beryllium"
        AtomsLabels(5) = "boron"
        AtomsLabels(6) = "carbon"
        AtomsLabels(7) = "nitrogen"
        AtomsLabels(8) = "oxygen"
        AtomsLabels(9) = "fluorine"
        AtomsLabels(10) = "neon"
        AtomsLabels(11) = "sodium"
        AtomsLabels(12) = "magnesium"
        AtomsLabels(13) = "aluminium"
        AtomsLabels(14) = "silicon"
        AtomsLabels(15) = "phoshorus"
        AtomsLabels(16) = "sulfur"
        AtomsLabels(17) = "chlorine"
        AtomsLabels(18) = "argon"
        AtomsLabels(19) = "potassium"
        AtomsLabels(20) = "calcium"
        AtomsLabels(21) = "scandium"
        AtomsLabels(22) = "titanium"
        AtomsLabels(23) = "vanadium"
        AtomsLabels(24) = "chromium"
        AtomsLabels(25) = "manganese"
        AtomsLabels(26) = "iron"
        AtomsLabels(27) = "cobalt"
        AtomsLabels(28) = "nickel"
        AtomsLabels(29) = "copper"
        AtomsLabels(30) = "zinc"
        AtomsLabels(31) = "gallium"
        AtomsLabels(32) = "germanium"
        AtomsLabels(33) = "arsenic"
        AtomsLabels(34) = "selenium"
        AtomsLabels(35) = "bromine"
        AtomsLabels(36) = "krypton"
        AtomsLabels(37) = "rubidium"
        AtomsLabels(38) = "strontium"
        AtomsLabels(39) = "yttrium"
        AtomsLabels(40) = "zirconium"
        AtomsLabels(41) = "niobium"
        AtomsLabels(42) = "molybdenum"
        AtomsLabels(43) = "technetium"
        AtomsLabels(44) = "ruthenium"
        AtomsLabels(45) = "rhodium"
        AtomsLabels(46) = "palladium"
        AtomsLabels(47) = "silver"
        AtomsLabels(48) = "cadmium"
        AtomsLabels(49) = "indium"
        AtomsLabels(50) = "tin"
        AtomsLabels(51) = "antimony"
        AtomsLabels(52) = "tellurium"
        AtomsLabels(53) = "iodine"
        AtomsLabels(54) = "xenon"
        AtomsLabels(55) = "caesium"
        AtomsLabels(56) = "barium"
        AtomsLabels(57) = "lanthanum"
        AtomsLabels(58) = "Cerium"
        AtomsLabels(59) = "Praseodymium"
        AtomsLabels(60) = "Neodymium"
        AtomsLabels(61) = "Promethium"
        AtomsLabels(62) = "samarim"
        AtomsLabels(63) = "Europium"
        AtomsLabels(64) = "gadolinium"
        AtomsLabels(65) = "terbium"
        AtomsLabels(66) = "Dysprosium"
        AtomsLabels(67) = "Ho"
        AtomsLabels(68) = "Er"
        AtomsLabels(69) = "Tm"
        AtomsLabels(70) = "Yb"
        AtomsLabels(71) = "Lu"
        AtomsLabels(72) = "Hf"
        AtomsLabels(73) = "Ta"
        AtomsLabels(74) = "W"
        AtomsLabels(75) = "Re"
        AtomsLabels(76) = "Os"
        AtomsLabels(77) = "Ir"
        AtomsLabels(78) = "Pt"
        AtomsLabels(79) = "Au"
        AtomsLabels(80) = "Hg"
        AtomsLabels(81) = "Ti"
        AtomsLabels(82) = "Pb"
        AtomsLabels(83) = "Bi"
        AtomsLabels(84) = "Po"
        AtomsLabels(85) = "At"
        AtomsLabels(86) = "Rn"
        AtomsLabels(87) = "Fr"
        AtomsLabels(88) = "Ra"
        AtomsLabels(89) = "Ac"
        AtomsLabels(90) = "Th"
        AtomsLabels(91) = "Pa"
        AtomsLabels(92) = "U"
        AtomsLabels(93) = "Np"
        AtomsLabels(94) = "Pu"
        AtomsLabels(95) = "Am"
        AtomsLabels(96) = "Cm"
        AtomsLabels(97) = "Bk"
        AtomsLabels(98) = "Cf"
        AtomsLabels(99) = "Es"
        AtomsLabels(100) = "Fm"
        AtomsLabels(101) = "101"
        AtomsLabels(102) = "No"
        AtomsLabels(103) = "Lr"
    end subroutine load_label

end module atom_label
