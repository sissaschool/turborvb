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

function atom_weight(atomic_number)
    ! ------------------------------------------------------------------
    !
    implicit none
    integer :: atomic_number
    real*8 :: atom_weight

    real*8 :: weights(103)
    data weights/1.00794d0, 4.00260d0, &
        6.941d0, 9.01218d0, 10.811d0, 12.0107d0, 14.00674d0, &
        15.9994d0, 18.99840d0, 20.1797d0, &
        22.98977d0, 24.3050d0, 26.98154d0, 28.0855d0, 30.97376d0, &
        32.066d0, 35.4527d0, 39.948d0, &
        39.0983d0, 40.078d0, 44.95591d0, 47.867d0, 50.9415d0, &
        51.9961d0, 54.93805d0, 55.845d0, &
        58.93320d0, 58.6934d0, 63.546d0, 65.39d0, &
        69.723d0, 72.61d0, 74.92160d0, 78.96d0, 79.904d0, 83.80d0, &
        85.4678d0, 87.62d0, 88.90585d0, 91.224d0, 92.90638d0, &
        95.94d0, 98.d0, &
        101.07d0, 102.90550d0, 106.42d0, 107.8682d0, 112.411d0, &
        114.818d0, 118.710d0, 121.760d0, 127.60d0, 126.90447d0, &
        131.29d0, &
        132.90545d0, 137.327d0, 138.9055d0, 140.116d0, 140.90765d0, &
        144.24d0, 145.d0, 150.36d0, 151.964d0, 157.25d0, &
        158.92534d0, 162.50d0, 164.93032d0, 167.26d0, &
        168.93421d0, 173.04d0, 174.967d0, &
        178.49d0, 180.9479d0, 183.84d0, 186.207d0, 190.23d0, &
        192.217d0, 195.078d0, 196.96655d0, 200.59d0, &
        204.3833d0, 207.2d0, 208.98038d0, 209.d0, 210.d0, 222.d0, &
        223.d0, 226.d0, 227.d0, 232.0381d0, 231.03588d0, &
        238.0289d0, 237.d0, 244.d0, &
        243.d0, 247.d0, 247.d0, 251.d0, 252.d0, 257.d0, &
        258.d0, 259.d0, 262.d0/

    if (atomic_number < 1 .or. atomic_number > 103) then
        call errore('atom_name', 'invalid atomic number', 1000 + atomic_number)
    else
        atom_weight = weights(atomic_number)
    end if
    return

end function atom_weight
!
