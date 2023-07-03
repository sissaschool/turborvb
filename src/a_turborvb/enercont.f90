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

function enercont(indt, nelup, neldo, winvup, winvdo, tabpip)
    implicit none
    integer i, j, jj, indt, nelup, neldo
    real*8 winvup(nelup, *), winvdo(neldo, *)                            &
            &, tabpip(nelup + neldo, *), enercont

    enercont = 0.d0

    !       write(6,*) ' Laplacian up '

    do i = 1, nelup
        enercont = enercont - winvup(i, indt + 4) - tabpip(i, indt + 4)
        !       write(6,*) i,tabpip(i,indt+4),winvup(i,indt+4)

    end do

    !       write(6,*) ' Gradient up  '
    do j = 1, 3
        jj = indt + j
        do i = 1, nelup
            enercont = enercont - 2.d0*tabpip(i, jj)*winvup(i, jj)                &
                    & - tabpip(i, jj)*tabpip(i, jj)
            !        write(6,*) i,jj,tabpip(i,jj),winvup(i,jj)
        end do
    end do

    !       write(6,*) ' Laplacian down '

    do i = 1, neldo
        enercont = enercont - winvdo(i, indt + 4) - tabpip(i + nelup, indt + 4)
        !       write(6,*) i,tabpip(i+nelup,indt+4),winvdo(i,indt+4)
    end do
    !       write(6,*) ' Gradient down   '

    do j = 1, 3
        jj = indt + j
        do i = 1, neldo
            enercont = enercont - 2.d0*tabpip(i + nelup, jj)*winvdo(i, jj)          &
                    & - tabpip(i + nelup, jj)*tabpip(i + nelup, jj)
            !       write(6,*) i,jj,tabpip(i+nelup,jj),winvdo(i,jj)
        end do
    end do

    return
end
