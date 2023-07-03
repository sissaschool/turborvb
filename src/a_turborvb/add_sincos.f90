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
!
! Sandro Sorella created on 27th Nov. 2007.

subroutine add_sincos(zeta, sumImgReal, kSq, Vk                      &
        &, kxbnd, kybnd, kzbnd, vsin, vcos)
    implicit none
    real(8) startsiny, startcosy, sine, cosn, bufsine                     &
            &, deltacosx, deltasinx                                              &
            &, deltacosy, deltasiny, deltacosz, deltasinz, startx, startsinx         &
            &, startcosx, k2, zeta, sumImgReal(*), DX, kSQ, Vk(*)                     &
            &, kcf, vsin(*), vcos(*)
    integer kx, ky, kz, kybnd, kxbnd, kzbnd
    parameter(DX=1.d-12)
    integer cnt

    cnt = 1

    deltacosx = vcos(1)
    deltasinx = vsin(1)
    deltacosy = vcos(2)
    deltasiny = vsin(2)
    deltacosz = vcos(3)
    deltasinz = vsin(3)
    startsinx = vsin(4)
    startcosx = vcos(4)

    do kx = 0, kxbnd

        startsiny = startsinx
        startcosy = startcosx
        startsinx = startsinx*deltacosx + startcosx*deltasinx
        startcosx = startcosx*deltacosx - startsiny*deltasinx

        do ky = -kybnd, kybnd

            sine = startsiny
            cosn = startcosy
            startsiny = startsiny*deltacosy + startcosy*deltasiny
            startcosy = startcosy*deltacosy - sine*deltasiny

            do kz = -kzbnd, kzbnd

                k2 = kx*kx + ky*ky + kz*kz
                if (k2 .lt. kSq .and. k2 .gt. DX) then
                    sumImgReal(2*cnt) = sumImgReal(2*cnt) + zeta*sine*Vk(cnt)
                    sumImgReal(2*cnt - 1) = sumImgReal(2*cnt - 1) + zeta*cosn*Vk(cnt)
                    cnt = cnt + 1

                end if
                bufsine = sine
                sine = sine*deltacosz + cosn*deltasinz
                cosn = cosn*deltacosz - bufsine*deltasinz
            end do
        end do
    end do

    return
end
