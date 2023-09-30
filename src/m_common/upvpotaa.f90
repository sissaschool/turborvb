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

function upvpotaa(zeta, iond, nion, kappa, LBox)
    use dielectric, only: veps, rep_erfc, epsilon0
    use allio, only: iond_cart, x_neigh, neigh, dist_shift, rank

    implicit none

    ! argument parameters
    integer, intent(in) :: nion
    real*8, intent(in) :: zeta(nion), iond(nion, nion), kappa, LBox

    ! local variables
    integer nel, i, j, ii, jj
    real*8 :: derfc, x_shift(3), cost_z
    real*8 pot_aa, eself1b
    real*8 upvpotaa

    double precision, parameter :: PI = 3.14159265358979323846d0

    pot_aa = 0.d0

    if (LBox .le. 0.d0) then
        do i = 1, nion
            do j = i + 1, nion
                if (zeta(i)*zeta(j) .ne. 0.d0) &
                    pot_aa = pot_aa + 2.d0*zeta(i)*zeta(j)*veps(iond(i, j))
            end do
        end do
        ! ********  EWALD REAL PART *************
        eself1b = 0.d0
    else
        do i = 1, nion
            do j = i + 1, nion
                cost_z = 2.d0*zeta(i)*zeta(j)
                jj = nion*(j - 1) + i
                if (cost_z .ne. 0.d0) then
!                       & pot_aa = pot_aa + 2.d0 * zeta(i) * zeta(j) * derfc(kappa * iond(i, j)) / iond(i, j)
!                        pot_aa = pot_aa + cost_z*rep_erfc(iond(i, j),kappa)
#ifdef _SIMD
!$omp simd
#endif
                    do ii = 1, neigh
                        dist_shift(ii) = dsqrt((iond_cart(1, jj) + x_neigh(ii, 1))**2 + &
                                               (iond_cart(2, jj) + x_neigh(ii, 2))**2 + &
                                               (iond_cart(3, jj) + x_neigh(ii, 3))**2)
                        dist_shift(ii) = cost_z*rep_erfc(dist_shift(ii), kappa)
                    end do
!                   pot_aa = pot_aa + cost_z*rep_erfc(dist_shift(ii),kappa)
                    pot_aa = pot_aa + sum(dist_shift(1:neigh))
                end if
            end do
        end do
        eself1b = -2*kappa*epsilon0/dsqrt(pi)*sum(zeta(1:nion)**2)

    end if
    ! subtracting the self interaction contribution

    upvpotaa = pot_aa + eself1b
    !     write(6,*) ' kappa inside =',kappa
    !     write(6,*) ' E_sr Qbox conv =',pot_aa/2.d0
    !     write(6,*) ' Eself Qbox conv (H) =',-eself1b/2.d0

    return
end
