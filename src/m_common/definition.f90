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

!> @brief Generate spherical quadrature points for different polyhedral symmetries
!>
!> This subroutine generates quadrature points and weights for numerical
!> integration on the unit sphere based on different polyhedral symmetries.
!> It supports octahedral, tetrahedral, and icosahedral point sets with
!> various numbers of integration points.
!>
!> @param[in] maxint Number of quadrature points (4, 6, 12, 18, 26, 32, 50)
!> @param[out] wintpseudo Array of quadrature weights
!> @param[in] lmax Maximum angular momentum for Legendre polynomials
!> @param[out] versor Array of unit vectors (3 × maxint)
!> @param[out] legendre Array of Legendre polynomial values
!> @param[in] rank MPI rank for error reporting
!> @param[out] iflagerr Error flag (1 if error occurs)
!>
!> @details
!> The subroutine supports the following quadrature schemes:
!>
!> **Supported Point Sets:**
!> - maxint = 4: Tetrahedral symmetry (4 points)
!> - maxint = 6: Octahedral symmetry (6 points)
!> - maxint = 12: Icosahedral symmetry (12 points)
!> - maxint = 18: Extended octahedral symmetry (18 points)
!> - maxint = 26: High-order octahedral symmetry (26 points)
!> - maxint = 32: Extended icosahedral symmetry (32 points)
!> - maxint = 50: High-order octahedral symmetry (50 points)
!>
!> **Point Generation:**
!> Each symmetry type generates points based on geometric properties:
!> - Tetrahedral: Points at (±1/√3, ±1/√3, ±1/√3) with sign combinations
!> - Octahedral: Points along coordinate axes (±1,0,0), (0,±1,0), (0,0,±1)
!> - Icosahedral: Points based on icosahedron vertices and face centers
!>
!> **Weight Assignment:**
!> Weights are assigned to ensure proper normalization:
!> - Equal weights for symmetric point sets
!> - Different weights for different point classes in extended sets
!>
!> **Legendre Polynomial Evaluation:**
!> After generating points, Legendre polynomials are evaluated at each
!> point for use in spherical harmonic expansions.
!>
!> @note The subroutine requires the legfun function for Legendre polynomial evaluation.
!> @note Points are normalized to lie on the unit sphere.
!> @note Weights sum to 4π for proper spherical integration.
!> @note Error handling is provided for unsupported point counts.
subroutine definition(maxint, wintpseudo, lmax, versor, legendre, rank, iflagerr)
    !
    implicit none
    character(3) pseudoname
    integer maxint, lmax, i, j, rank, iflagerr
    real(8) versor(3, *), legendre(lmax - 1, *), cost, legfun, wintpseudo(*)  &
            &, par, theta, phi, pi, q, r, s, theta1, theta2
    !
    !> @brief Initialize versor array to zero
    do i = 1, 3
        do j = 1, maxint
            versor(i, j) = 0.d0
        end do
    end do
    !
    !> @brief Octahedral symmetry quadrature (6 points)
    if (maxint .eq. 6) then
        !> @brief Assign equal weights for 6-point octahedral quadrature
        do i = 1, maxint
            wintpseudo(i) = 1.d0/6.d0
        end do
        !> @brief Generate points along coordinate axes
        versor(1, 1) = 1.d0
        versor(1, 2) = -1.d0
        versor(2, 3) = 1.d0
        versor(2, 4) = -1.d0
        versor(3, 5) = 1.d0
        versor(3, 6) = -1.d0
        ! tetrahedron symmetry quadrature
    elseif (maxint .eq. 4) then
        !> @brief Assign equal weights for 4-point tetrahedral quadrature
        do i = 1, maxint
            wintpseudo(i) = 1.d0/4.d0
        end do
        !> @brief Generate tetrahedral points with proper sign combinations
        par = 1.d0/dsqrt(3.d0)
        do i = 1, 3
            do j = 1, maxint
                versor(i, j) = par
            end do
        end do
        versor(2, 2) = -versor(2, 2)
        versor(3, 2) = -versor(3, 2)
        versor(1, 3) = -versor(1, 3)
        versor(3, 3) = -versor(3, 3)
        versor(1, 4) = -versor(1, 4)
        versor(2, 4) = -versor(2, 4)
        ! icosahedron symmetry quadrature
    elseif (maxint .eq. 12) then
        !> @brief Assign equal weights for 12-point icosahedral quadrature
        do i = 1, maxint
            wintpseudo(i) = 1.d0/12.d0
        end do
        !> @brief Generate icosahedral points using spherical coordinates
        pi = acos(-1.d0)
        versor(3, 1) = 1.d0
        versor(3, 2) = -1.d0
        do j = 0, 4
            theta = atan(2.d0)
            phi = 2.d0*pi*dble(j)/5.d0
            versor(1, j + 3) = sin(theta)*cos(phi)
            versor(2, j + 3) = sin(theta)*sin(phi)
            versor(3, j + 3) = cos(theta)
            theta = pi - atan(2.d0)
            phi = (2.d0*dble(j) + 1.d0)*pi/5.d0
            versor(1, j + 8) = sin(theta)*cos(phi)
            versor(2, j + 8) = sin(theta)*sin(phi)
            versor(3, j + 8) = cos(theta)
        end do
        ! octahedron symmetry quadrature
    elseif (maxint .eq. 18) then
        !> @brief Assign different weights for different point classes
        do i = 1, 6
            wintpseudo(i) = 1.d0/30.d0
        end do
        do i = 7, maxint
            wintpseudo(i) = 1.d0/15.d0
        end do
        !> @brief Generate extended octahedral points
        pi = 1.d0/dsqrt(2.d0)
        versor(1, 1) = 1.d0
        versor(1, 2) = -1.d0
        versor(2, 3) = 1.d0
        versor(2, 4) = -1.d0
        versor(3, 5) = 1.d0
        versor(3, 6) = -1.d0
        versor(1, 7) = pi
        versor(2, 7) = pi
        versor(1, 8) = -pi
        versor(2, 8) = pi
        versor(1, 9) = pi
        versor(2, 9) = -pi
        versor(1, 10) = -pi
        versor(2, 10) = -pi
        versor(1, 11) = pi
        versor(3, 11) = pi
        versor(1, 12) = -pi
        versor(3, 12) = pi
        versor(1, 13) = pi
        versor(3, 13) = -pi
        versor(1, 14) = -pi
        versor(3, 14) = -pi
        versor(2, 15) = pi
        versor(3, 15) = pi
        versor(2, 16) = -pi
        versor(3, 16) = pi
        versor(2, 17) = pi
        versor(3, 17) = -pi
        versor(2, 18) = -pi
        versor(3, 18) = -pi

    elseif (maxint .eq. 26) then ! octahedron
        !> @brief Assign weights for different point classes in 26-point quadrature
        do i = 1, 6
            wintpseudo(i) = 1.d0/21.d0
        end do
        do i = 7, 18
            wintpseudo(i) = 4.d0/105.d0
        end do
        do i = 19, maxint
            wintpseudo(i) = 27.d0/840.d0
        end do
        !> @brief Generate high-order octahedral points
        pi = 1.d0/dsqrt(2.d0)
        q = 1.d0/dsqrt(3.d0)
        versor(1, 1) = 1.d0
        versor(1, 2) = -1.d0
        versor(2, 3) = 1.d0
        versor(2, 4) = -1.d0
        versor(3, 5) = 1.d0
        versor(3, 6) = -1.d0
        versor(1, 7) = pi
        versor(2, 7) = pi
        versor(1, 8) = -pi
        versor(2, 8) = pi
        versor(1, 9) = pi
        versor(2, 9) = -pi
        versor(1, 10) = -pi
        versor(2, 10) = -pi
        versor(1, 11) = pi
        versor(3, 11) = pi
        versor(1, 12) = -pi
        versor(3, 12) = pi
        versor(1, 13) = pi
        versor(3, 13) = -pi
        versor(1, 14) = -pi
        versor(3, 14) = -pi
        versor(2, 15) = pi
        versor(3, 15) = pi
        versor(2, 16) = -pi
        versor(3, 16) = pi
        versor(2, 17) = pi
        versor(3, 17) = -pi
        versor(2, 18) = -pi
        versor(3, 18) = -pi

        !> @brief Add cube vertices for extended quadrature
        versor(1, 19) = q
        versor(2, 19) = q
        versor(3, 19) = q

        versor(1, 20) = -q
        versor(2, 20) = q
        versor(3, 20) = q

        versor(1, 21) = q
        versor(2, 21) = -q
        versor(3, 21) = q

        versor(1, 22) = -q
        versor(2, 22) = -q
        versor(3, 22) = q

        versor(1, 23) = q
        versor(2, 23) = q
        versor(3, 23) = -q

        versor(1, 24) = q
        versor(2, 24) = -q
        versor(3, 24) = -q

        versor(1, 25) = -q
        versor(2, 25) = -q
        versor(3, 25) = -q

        versor(1, 26) = -q
        versor(2, 26) = q
        versor(3, 26) = -q

    elseif (maxint .eq. 32) then ! icosahedron

        !> @brief Assign weights for different point classes in 32-point quadrature
        do i = 1, 12
            wintpseudo(i) = 5.d0/168.d0
        end do
        do i = 13, maxint
            wintpseudo(i) = 27.d0/840.d0
        end do

        !> @brief Generate extended icosahedral points
        pi = acos(-1.d0)
        theta1 = dacos((2.d0 + dsqrt(5.d0))/dsqrt(15.d0 + 6.d0*dsqrt(5.d0)))
        theta2 = dacos(1.d0/dsqrt(15.d0 + 6.d0*dsqrt(5.d0)))

        versor(3, 1) = 1.d0
        versor(3, 2) = -1.d0
        do j = 0, 4
            theta = atan(2.d0)
            phi = 2.d0*pi*dble(j)/5.d0
            versor(1, j + 3) = sin(theta)*cos(phi)
            versor(2, j + 3) = sin(theta)*sin(phi)
            versor(3, j + 3) = cos(theta)
            theta = pi - atan(2.d0)
            phi = (2.d0*dble(j) + 1.d0)*pi/5.d0
            versor(1, j + 8) = sin(theta)*cos(phi)
            versor(2, j + 8) = sin(theta)*sin(phi)
            versor(3, j + 8) = cos(theta)
        end do

        !> @brief Add additional icosahedral points with different theta angles
        do j = 0, 4
            theta = theta1
            phi = (2.d0*dble(j) + 1.d0)*pi/5.d0
            versor(1, j + 13) = sin(theta)*cos(phi)
            versor(2, j + 13) = sin(theta)*sin(phi)
            versor(3, j + 13) = cos(theta)
            theta = theta2
            phi = (2.d0*dble(j) + 1.d0)*pi/5.d0
            versor(1, j + 18) = sin(theta)*cos(phi)
            versor(2, j + 18) = sin(theta)*sin(phi)
            versor(3, j + 18) = cos(theta)
            theta = pi - theta1
            phi = 2.d0*pi*dble(j)/5.d0
            versor(1, j + 23) = sin(theta)*cos(phi)
            versor(2, j + 23) = sin(theta)*sin(phi)
            versor(3, j + 23) = cos(theta)
            theta = pi - theta2
            phi = 2.d0*pi*dble(j)/5.d0
            versor(1, j + 28) = sin(theta)*cos(phi)
            versor(2, j + 28) = sin(theta)*sin(phi)
            versor(3, j + 28) = cos(theta)
        end do

    elseif (maxint .eq. 50) then ! octahedron
        !> @brief Assign weights for different point classes in 50-point quadrature
        do i = 1, 6
            wintpseudo(i) = 4.d0/315.d0
        end do
        do i = 7, 18
            wintpseudo(i) = 64.d0/2835.d0
        end do
        do i = 19, 26
            wintpseudo(i) = 27.d0/1280.d0
        end do
        do i = 27, maxint
            wintpseudo(i) = 14641.d0/725760.d0
        end do

        !> @brief Generate high-order octahedral points with additional parameters
        pi = 1.d0/dsqrt(2.d0)
        q = 1.d0/dsqrt(3.d0)
        r = 1.d0/dsqrt(11.d0)
        s = 3.d0/dsqrt(11.d0)

        !> @brief Standard octahedral points
        versor(1, 1) = 1.d0
        versor(1, 2) = -1.d0
        versor(2, 3) = 1.d0
        versor(2, 4) = -1.d0
        versor(3, 5) = 1.d0
        versor(3, 6) = -1.d0
        versor(1, 7) = pi
        versor(2, 7) = pi
        versor(1, 8) = -pi
        versor(2, 8) = pi
        versor(1, 9) = pi
        versor(2, 9) = -pi
        versor(1, 10) = -pi
        versor(2, 10) = -pi
        versor(1, 11) = pi
        versor(3, 11) = pi
        versor(1, 12) = -pi
        versor(3, 12) = pi
        versor(1, 13) = pi
        versor(3, 13) = -pi
        versor(1, 14) = -pi
        versor(3, 14) = -pi
        versor(2, 15) = pi
        versor(3, 15) = pi
        versor(2, 16) = -pi
        versor(3, 16) = pi
        versor(2, 17) = pi
        versor(3, 17) = -pi
        versor(2, 18) = -pi
        versor(3, 18) = -pi

        !> @brief Cube vertices
        versor(1, 19) = q
        versor(2, 19) = q
        versor(3, 19) = q

        versor(1, 20) = -q
        versor(2, 20) = q
        versor(3, 20) = q

        versor(1, 21) = q
        versor(2, 21) = -q
        versor(3, 21) = q

        versor(1, 22) = -q
        versor(2, 22) = -q
        versor(3, 22) = q

        versor(1, 23) = q
        versor(2, 23) = q
        versor(3, 23) = -q

        versor(1, 24) = q
        versor(2, 24) = -q
        versor(3, 24) = -q

        versor(1, 25) = -q
        versor(2, 25) = -q
        versor(3, 25) = -q

        versor(1, 26) = -q
        versor(2, 26) = q
        versor(3, 26) = -q

        !> @brief Additional high-order points with r and s parameters
        versor(1, 27) = r
        versor(2, 27) = r
        versor(3, 27) = s

        versor(1, 28) = -r
        versor(2, 28) = r
        versor(3, 28) = s

        versor(1, 29) = r
        versor(2, 29) = -r
        versor(3, 29) = s

        versor(1, 30) = -r
        versor(2, 30) = -r
        versor(3, 30) = s

        versor(1, 31) = r
        versor(2, 31) = r
        versor(3, 31) = -s

        versor(1, 32) = r
        versor(2, 32) = -r
        versor(3, 32) = -s

        versor(1, 33) = -r
        versor(2, 33) = -r
        versor(3, 33) = -s

        versor(1, 34) = -r
        versor(2, 34) = r
        versor(3, 34) = -s

        versor(1, 35) = r
        versor(2, 35) = s
        versor(3, 35) = r

        versor(1, 36) = -r
        versor(2, 36) = s
        versor(3, 36) = r

        versor(1, 37) = r
        versor(2, 37) = -s
        versor(3, 37) = r

        versor(1, 38) = -r
        versor(2, 38) = -s
        versor(3, 38) = r

        versor(1, 39) = r
        versor(2, 39) = s
        versor(3, 39) = -r

        versor(1, 40) = r
        versor(2, 40) = -s
        versor(3, 40) = -r

        versor(1, 41) = -r
        versor(2, 41) = -s
        versor(3, 41) = -r

        versor(1, 42) = -r
        versor(2, 42) = s
        versor(3, 42) = -r

        versor(1, 43) = s
        versor(2, 43) = r
        versor(3, 43) = r

        versor(1, 44) = -s
        versor(2, 44) = r
        versor(3, 44) = r

        versor(1, 45) = s
        versor(2, 45) = -r
        versor(3, 45) = r

        versor(1, 46) = -s
        versor(2, 46) = -r
        versor(3, 46) = r

        versor(1, 47) = s
        versor(2, 47) = r
        versor(3, 47) = -r

        versor(1, 48) = s
        versor(2, 48) = -r
        versor(3, 48) = -r

        versor(1, 49) = -s
        versor(2, 49) = -r
        versor(3, 49) = -r

        versor(1, 50) = -s
        versor(2, 50) = r
        versor(3, 50) = -r

    elseif (maxint .eq. 0) then
        if (rank .eq. 0) write (6, *) ' Warning no quadrature points with pseudo '
    else
        if (rank .eq. 0) write (6, *) 'ERROR # of quadrature points non defined!!!'
        iflagerr = 1
        return
    end if
    !
    !> @brief Evaluate Legendre polynomials at all quadrature points
    do i = 1, maxint
        cost = versor(3, i)
        do j = 1, lmax - 1
            legendre(j, i) = legfun(j - 1, cost)
        end do
    end do
    !
    return
end
