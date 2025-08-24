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

!> @brief Find zero crossings of Green's function using bisection method
!> @details This subroutine implements a bisection algorithm to find the zero crossings
!> of the Green's function. It iteratively refines the interval [eig_min, eig_max] by
!> evaluating the Green's function at the midpoint and updating the interval bounds
!> based on the sign of the function value. The algorithm uses a maximum of 50 iterations
!> to ensure convergence.
!> @param[in] eig_min Minimum eigenvalue bound
!> @param[in] eig_max Maximum eigenvalue bound
!> @param[in] npm Leading dimension of matrix psi
!> @param[in] ndim Dimension of the problem
!> @param[in] psi Eigenvector matrix
!> @param[in] eig Eigenvalue array
!> @param[in,out] emin Array of minimum bounds for each dimension
!> @param[in,out] emax Array of maximum bounds for each dimension
!> @param[in,out] green Green's function values
!> @param[in,out] e Energy values for Green's function evaluation
subroutine findzero(eig_min, eig_max, npm, ndim, psi, eig         &
        &, emin, emax, green, e)
    integer npm, ndim, k, i, maxit
    real*8 psi(npm, *), eig(*), green(*)                            &
            &, eig_min, eig_max, emin(*), emax(*), e(*)
    maxit = 50
    call dscalzero(ndim, eig_min, emin, 1)
    call dscalzero(ndim, eig_max, emax, 1)
    do i = 1, maxit
        do k = 1, ndim
            e(k) = (emin(k) + emax(k))*0.5d0
        end do
        call evalgreen(npm, ndim, psi, eig, e, green)
        do k = 1, ndim
            if (green(k) .gt. 0.d0) then
                emax(k) = e(k)
            else
                emin(k) = e(k)
            end if
        end do
    end do

    return
end

!> @brief Evaluate Green's function for given energy values
!> @details This subroutine computes the Green's function values for given energy points e.
!> The Green's function is calculated as a sum over eigenstates using the formula:
!> G(e) = sum_l |psi(k,l)|^2 / (eig(l) - e(k)), where psi(k,l) are the eigenvector
!> components and eig(l) are the eigenvalues.
!> @param[in] npm Leading dimension of matrix psi
!> @param[in] ndim Dimension of the problem
!> @param[in] psi Eigenvector matrix
!> @param[in] eig Eigenvalue array
!> @param[in] e Energy values for evaluation
!> @param[out] green Green's function values
subroutine evalgreen(npm, ndim, psi, eig, e, green)
    implicit none
    integer npm, ndim, k_no, i, j, k, l
    real*8 psi(npm, *), eig(*), green(*), e(*)
    do k = 1, ndim
        green(k) = psi(k, 1)**2/(eig(1) - e(k))
    end do
    do l = 2, ndim
        do k = 1, ndim
            green(k) = green(k) + psi(k, l)**2/(eig(l) - e(k))
        end do
    end do
    return
end
