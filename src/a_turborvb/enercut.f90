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

!> @brief      Apply energy cutoff using arctangent function
!> @details    This function applies a smooth energy cutoff using an arctangent
!>             function. If cut=0, it returns the original energy. Otherwise,
!>             it applies the transformation: etry + cut*atan((ener-etry)/cut).
!> @param[in]  ener     Input energy value
!> @param[in]  etry     Reference energy for cutoff
!> @param[in]  cut      Cutoff parameter (0 for no cutoff)
!> @return     enercut  Energy value after cutoff transformation
function enercut(ener, etry, cut)
    implicit none
    real*8 enercut, ener, etry, cut, datan, argtan
    if (cut .eq. 0.d0) then
        enercut = ener
        return
    end if
    argtan = (ener - etry)/cut
    enercut = etry + cut*datan(argtan)
    return
end
