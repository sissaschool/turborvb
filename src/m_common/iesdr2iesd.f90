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

function iesdr2iesd(iesdr)
    implicit none
    integer iesdr, iesdr2iesd
    iesdr2iesd = iesdr
    if (iesdr .eq. -15) iesdr2iesd = 1
    if (iesdr .eq. -5) iesdr2iesd = 1
    if (iesdr .eq. -6) iesdr2iesd = 4
    if (iesdr .eq. -7) iesdr2iesd = 4
    if (iesdr .eq. -17) iesdr2iesd = -1
    if (iesdr .eq. -20) iesdr2iesd = -1
    if (iesdr .eq. -21) iesdr2iesd = 2
    return
end
