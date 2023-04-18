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

! ------------------------------------------------------------------
function atom_weight(atomic_number)
    ! ------------------------------------------------------------------
    !
    use kinds, only: DP
    implicit none
    integer :: atomic_number
    real(DP) :: atom_weight

    real(DP) :: weights(103)
    data weights/1.00794_dp, 4.00260_dp, &
        6.941_dp, 9.01218_dp, 10.811_dp, 12.0107_dp, 14.00674_dp, &
        15.9994_dp, 18.99840_dp, 20.1797_dp, &
        22.98977_dp, 24.3050_dp, 26.98154_dp, 28.0855_dp, 30.97376_dp, &
        32.066_dp, 35.4527_dp, 39.948_dp, &
        39.0983_dp, 40.078_dp, 44.95591_dp, 47.867_dp, 50.9415_dp, &
        51.9961_dp, 54.93805_dp, 55.845_dp, &
        58.93320_dp, 58.6934_dp, 63.546_dp, 65.39_dp, &
        69.723_dp, 72.61_dp, 74.92160_dp, 78.96_dp, 79.904_dp, 83.80_dp, &
        85.4678_dp, 87.62_dp, 88.90585_dp, 91.224_dp, 92.90638_dp, &
        95.94_dp, 98._dp, &
        101.07_dp, 102.90550_dp, 106.42_dp, 107.8682_dp, 112.411_dp, &
        114.818_dp, 118.710_dp, 121.760_dp, 127.60_dp, 126.90447_dp, &
        131.29_dp, &
        132.90545_dp, 137.327_dp, 138.9055_dp, 140.116_dp, 140.90765_dp, &
        144.24_dp, 145._dp, 150.36_dp, 151.964_dp, 157.25_dp, &
        158.92534_dp, 162.50_dp, 164.93032_dp, 167.26_dp, &
        168.93421_dp, 173.04_dp, 174.967_dp, &
        178.49_dp, 180.9479_dp, 183.84_dp, 186.207_dp, 190.23_dp, &
        192.217_dp, 195.078_dp, 196.96655_dp, 200.59_dp, &
        204.3833_dp, 207.2_dp, 208.98038_dp, 209._dp, 210._dp, 222._dp, &
        223._dp, 226._dp, 227._dp, 232.0381_dp, 231.03588_dp, &
        238.0289_dp, 237._dp, 244._dp, &
        243._dp, 247._dp, 247._dp, 251._dp, 252._dp, 257._dp, &
        258._dp, 259._dp, 262._dp/

    if (atomic_number < 1 .or. atomic_number > 103) then
        call errore('atom_name', 'invalid atomic number', 1000 + atomic_number)
    else
        atom_weight = weights(atomic_number)
    end if
    return

end function atom_weight
!
