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

function ngivej(rcart1, rcart2, LBox)

    use Cell, only: cellscale, s2r, car2cry

    implicit none
    real*8 dxyz(3), ngivej, dsqrt
    real*8 rcart1(3), rcart2(3), vecscra(3), LBox
    !     mapping on a box 0<= x < Lx , 0 <= y < Ly , 0 <= z < Lz
    dxyz(:) = rcart1(:) - rcart2(:)

    if (LBox .gt. 0.d0) then
        !      vecscra(:)=dxyz(:)
        !      call CartesianToCrystal(vecscra,1)
        vecscra(:) = car2cry(:, 1)*dxyz(1) + car2cry(:, 2)*dxyz(2) + car2cry(:, 3)*dxyz(3)
        vecscra(1) = anint(vecscra(1)/cellscale(1))
        vecscra(2) = anint(vecscra(2)/cellscale(2))
        vecscra(3) = anint(vecscra(3)/cellscale(3))
        dxyz(:) = dxyz(:) - s2r(:, 1)*vecscra(1) - s2r(:, 2)*vecscra(2) - s2r(:, 3)*vecscra(3)
        !      call dgemv('N',3,3,-1.d0,s2r,3,vecscra,1,1.d0,dxyz,1)
    end if

    ngivej = dsqrt(dxyz(1)**2 + dxyz(2)**2 + dxyz(3)**2)

    if (ngivej .lt. 1d-9) ngivej = 1d-9

    !      write(6,*) ' ngive j =',ngivej
    return
end
subroutine dgivej(rcart1, rcart2, LBox, dxyz)

    use Cell, only: cellscale, s2r, car2cry

    implicit none
    real*8 dxyz(3), dsqrt
    real*8 rcart1(3), rcart2(3), vecscra(3), LBox
    !     mapping on a box 0<= x < Lx , 0 <= y < Ly , 0 <= z < Lz
    dxyz(:) = rcart1(:) - rcart2(:)

    if (LBox .gt. 0.d0) then
        !      vecscra(:)=dxyz(:)
        !      call CartesianToCrystal(vecscra,1)
        vecscra(:) = car2cry(:, 1)*dxyz(1) + car2cry(:, 2)*dxyz(2) + car2cry(:, 3)*dxyz(3)
        vecscra(1) = anint(vecscra(1)/cellscale(1))
        vecscra(2) = anint(vecscra(2)/cellscale(2))
        vecscra(3) = anint(vecscra(3)/cellscale(3))
        dxyz(:) = dxyz(:) - s2r(:, 1)*vecscra(1) - s2r(:, 2)*vecscra(2) - s2r(:, 3)*vecscra(3)
        !      call dgemv('N',3,3,-1.d0,s2r,3,vecscra,1,1.d0,dxyz,1)
    end if

    return
end
