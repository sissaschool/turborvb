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

!--------------------------------------------------------------------------!
! upkinc/upkinc_complex compute the kinetic (diagonal) part of the
! local energy for the global wave function.
! Main varibles: winvup(i,indt+1:indt+4) = gradients and laplacians for
!                                          up-spin electrons
!                winvdo(i+nelup,indt+1:indt+4) = gradients and laplacians
!                                                for down-spin electrons
!                enercont = total kinetic energy
!--------------------------------------------------------------------------!

!> @brief      Compute kinetic energy contribution for real wave function
!> @details    This subroutine computes the kinetic energy contribution from
!>             the Laplacian and gradient terms for both up and down electrons.
!>             It stores individual contributions in kinc array and sums them
!>             in enercont.
!> @param[in]  indt      Index offset for derivatives
!> @param[in]  nelup     Number of up electrons
!> @param[in]  neldo     Number of down electrons
!> @param[in]  winvup    Inverse matrix for up electrons
!> @param[in]  winvdo    Inverse matrix for down electrons
!> @param[in]  tabpip    Wave function derivatives table
!> @param[out] enercont  Total kinetic energy contribution
!> @param[out] kinc      Individual kinetic energy contributions per electron
subroutine upkinc(indt, nelup, neldo, winvup, winvdo, tabpip, enercont, kinc)
    implicit none
    integer i, j, jj, indt, nelup, neldo
    real*8 winvup(nelup, *), winvdo(neldo, *)                            &
            &, tabpip(nelup + neldo, *), enercont, kinc(nelup + neldo)

    ! write(6,*) ' Laplacian up '
    do i = 1, nelup
        kinc(i) = -winvup(i, indt + 4) - tabpip(i, indt + 4)
    end do
    !       write(6,*) ' Gradient up  '
    do j = 1, 3
        jj = indt + j
        do i = 1, nelup
            kinc(i) = kinc(i) - 2.d0*tabpip(i, jj)*winvup(i, jj)                &
                    & - tabpip(i, jj)*tabpip(i, jj)
        end do
    end do

    !       write(6,*) ' Laplacian down '

    do i = 1, neldo
        kinc(i + nelup) = -winvdo(i, indt + 4) - tabpip(i + nelup, indt + 4)
    end do
    !       write(6,*) ' Gradient down   '
    do j = 1, 3
        jj = indt + j
        do i = 1, neldo
            kinc(i + nelup) = kinc(i + nelup) - 2.d0*tabpip(i + nelup, jj)*winvdo(i, jj)&
                    & - tabpip(i + nelup, jj)*tabpip(i + nelup, jj)
        end do
    end do
    enercont = sum(kinc(:))
    return
end subroutine upkinc

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!

!> @brief      Compute kinetic energy contribution for complex wave function
!> @details    This subroutine computes the kinetic energy contribution from
!>             the Laplacian and gradient terms for complex wave functions.
!>             Similar to upkinc but handles complex matrices and returns
!>             complex energy contributions.
!> @param[in]  indt      Index offset for derivatives
!> @param[in]  nelup     Number of up electrons
!> @param[in]  neldo     Number of down electrons
!> @param[in]  winvup    Complex inverse matrix for up electrons
!> @param[in]  winvdo    Complex inverse matrix for down electrons
!> @param[in]  tabpip    Wave function derivatives table
!> @param[out] enercont  Total complex kinetic energy contribution
!> @param[out] kinc      Individual complex kinetic energy contributions per electron
subroutine upkinc_complex(indt, nelup, neldo, winvup, winvdo, tabpip, enercont, kinc)

    implicit none
    integer i, j, jj, indt, nelup, neldo
    real(8) tabpip(nelup + neldo, *)
    complex(8) winvup(nelup, *), winvdo(neldo, *), enercont, kinc(nelup + neldo)

    ! write(6,*) ' Laplacian up '
    do i = 1, nelup
        kinc(i) = -winvup(i, indt + 4) - tabpip(i, indt + 4)
    end do
    ! write(6,*) ' Gradient up  '
    do j = 1, 3
        jj = indt + j
        do i = 1, nelup
            kinc(i) = kinc(i) - 2.d0*tabpip(i, jj)*winvup(i, jj) &
                    & - tabpip(i, jj)*tabpip(i, jj)
        end do
    end do

    ! write(6,*) ' Laplacian down '
    do i = 1, neldo
        kinc(i + nelup) = -winvdo(i, indt + 4) - tabpip(i + nelup, indt + 4)
    end do
    ! write(6,*) ' Gradient down   '
    do j = 1, 3
        jj = indt + j
        do i = 1, neldo
            kinc(i + nelup) = kinc(i + nelup) - 2.d0*tabpip(i + nelup, jj)*winvdo(i, jj)&
                    & - tabpip(i + nelup, jj)*tabpip(i + nelup, jj)
        end do
    end do
    enercont = sum(kinc(:))
    return
end subroutine upkinc_complex
