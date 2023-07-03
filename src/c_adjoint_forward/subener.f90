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

subroutine subener(indt, nelup, neldo, winvup, winvdo, tabpip, enercont)

    use constants, only: ipc

    implicit none
    integer i, j, jj, indt, nelup, neldo
    real(8) winvup(ipc*nelup, *), winvdo(max(ipc*neldo, 1), *), tabpip(nelup + neldo, *), enercont(ipc)
    integer indel

    enercont = 0.d0

    ! defining indices for real/complex winvup/winvdo

    if (ipc .eq. 1) then

        do i = 1, nelup
            enercont(1) = enercont(1) - winvup(i, indt + 4) - tabpip(i, indt + 4)
        end do

        do j = 1, 3
            jj = indt + j
            do i = 1, nelup
                enercont(1) = enercont(1) - 2.d0*tabpip(i, jj)*winvup(i, jj) &
                              - tabpip(i, jj)*tabpip(i, jj)
                !                write(6,*) i,jj,tabpip(i,jj),winvup(i,jj)
            end do
        end do

        do i = 1, neldo
            enercont(1) = enercont(1) - winvdo(i, indt + 4) - tabpip(i + nelup, indt + 4)
            !    write(6,*) i,tabpip(i+nelup,indt+4),winvdo(i,indt+4)
        end do

        do j = 1, 3
            jj = indt + j
            do i = 1, neldo
                enercont(1) = enercont(1) - 2.d0*tabpip(i + nelup, jj)*winvdo(i, jj)          &
                        & - tabpip(i + nelup, jj)*tabpip(i + nelup, jj)
                !     write(6,*) i,jj,tabpip(i+nelup,jj),winvdo(i,jj)
            end do
        end do

    else ! complex case

        do i = 1, nelup
            indel = 2*i - 1
            enercont(1) = enercont(1) - winvup(indel, indt + 4) - tabpip(i, indt + 4)
            enercont(2) = enercont(2) - winvup(indel + 1, indt + 4)
        end do

        do j = 1, 3
            jj = indt + j
            do i = 1, nelup
                indel = 2*i - 1
                enercont(1) = enercont(1) - 2.d0*tabpip(i, jj)*winvup(indel, jj) &
                              - tabpip(i, jj)*tabpip(i, jj)
                enercont(2) = enercont(2) - 2.d0*tabpip(i, jj)*winvup(indel + 1, jj)
                !                write(6,*) i,jj,tabpip(i,jj),winvup(indel,jj)
            end do
        end do

        do i = 1, neldo
            indel = 2*i - 1
            enercont(1) = enercont(1) - winvdo(indel, indt + 4) - tabpip(i + nelup, indt + 4)
            enercont(2) = enercont(2) - winvdo(indel + 1, indt + 4)
        end do

        do j = 1, 3
            jj = indt + j
            do i = 1, neldo
                indel = 2*i - 1
                enercont(1) = enercont(1) - 2.d0*tabpip(i + nelup, jj)*winvdo(indel, jj)          &
                        & - tabpip(i + nelup, jj)*tabpip(i + nelup, jj)
                enercont(2) = enercont(2) - 2.d0*tabpip(i + nelup, jj)*winvdo(indel + 1, jj)
            end do
        end do

    end if

    return

end subroutine subener
