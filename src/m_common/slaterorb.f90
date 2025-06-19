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

!> @brief Orbital type classification for Slater vs Gaussian orbitals
!>
!> This function determines whether a given orbital type code corresponds
!> to a Slater-type orbital (STO) or a Gaussian-type orbital (GTO) in
!> quantum Monte Carlo calculations. It provides a systematic way to
!> classify orbital types based on their numerical identifiers.
!>
!> @param[in] ioptorb Orbital type code to be classified
!> @return Logical value: .true. for Slater orbitals, .false. for Gaussian orbitals
!>
!> @details
!> The function performs the following classification:
!> 1. Slater-type orbitals (return .true.):
!>    - Ranges: 1:8, 10:15, 20:29, 30:35, 38:43, 50, 55:57, 66, 70, 71, 80:89
!>    - Extended ranges: 121:123, 125:130, 133:144
!> 2. Gaussian-type orbitals (return .false.):
!>    - Ranges: 16:19, 36, 37, 44:49, 51:54, 58:59, 60:65, 68:69, 72, 73
!>    - Extended ranges: 100:105, 108, 109, 131, 132, 145:155, 161
!>    - Large ranges: 1000:1099, 2000:2099, 1100:1199, 2100:2199, 1200:1299, 2200:2299
!>    - Special case: 200 (constant orbital considered Gaussian)
!> 3. Unknown orbital types: Sets error flag and reports the unknown code
!>
!> @note The function uses case statements for efficient classification
!> @note Slater orbitals include s, p, d, f, g orbitals and their variants
!> @note Gaussian orbitals include contracted Gaussians and specialized forms
!> @note The constant orbital (ioptorb=200) is classified as Gaussian
!> @note Unknown orbital types trigger an error condition (iflagerr=1)
!> @note Used in wave function setup and orbital basis set management
function slaterorb(ioptorb)
    use constants, only: iflagerr
    implicit none
    integer ioptorb
    logical slaterorb
    select case (ioptorb)
        !     case(34,10,12,28,57,80) !  Slater s
        !     slaterorb=.true.
        !     case(20,22,82)             ! Slater p
        !     slaterorb=.true.
        !     case(30,33,84)             ! Slater d
        !     slaterorb=.true.
        !     case(70,86)                ! Slater f
        !     slaterorb=.true.
        !     case(55,88)                ! Slater g
        !     slaterorb=.true.
    case (1:8, 10:15, 20:29, 30:35, 38:43, 50, 55:57, 66, 70, 71, 80:89, &
          121:123, 125:130, 133:144)
        slaterorb = .true.
    case (16:19, 36, 37, 44:49, 51:54, 58:59, 60:65, 68:69, 72, 73, 100:105, 108, 109, &
          131, 132, 145:155, 161, 1000:1099, 2000:2099, 1100:1199, &
          2100:2199, 1200:1299, 2200:2299, 200)
! the constant orbital is considered Gaussian
        slaterorb = .false.
    case default
        !     Orbital not found
        write (6, *) ' ERROR orbital not Slater neither Gaussian, &
                &              check slaterorb function, ioptorb= ', ioptorb
        iflagerr = 1
    end select
    return
end
