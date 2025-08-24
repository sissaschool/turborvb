! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2000 Wolfgang Schadow
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

!> @file bspline90_22.f90
!> @brief B-spline interpolation library for one, two, and three dimensions
!> @author Wolfgang Schadow (original), TurboRVB group (modifications)
!> @date 2000 (original), 2022 (modifications)
!> @version 2.2
!> @details
!> This library contains routines for B-spline interpolation in one, two, and
!> three dimensions. Part of the routines are based on the book by Carl de Boor:
!> "A practical guide to Splines" (Springer, New-York 1978) and have the same
!> calling sequence and names as the corresponding routines from the IMSL library.
!> 
!> The library provides:
!> - 1D B-spline interpolation and evaluation
!> - 2D tensor-product B-spline interpolation and evaluation  
!> - 3D tensor-product B-spline interpolation and evaluation
!> - Derivative evaluation for all dimensions
!> - Grid-based evaluation capabilities
!>
!> @note Results may vary slightly on different architectures due to floating-point
!> precision differences.
!>
!> @see de Boor, C. (1978). A practical guide to Splines. Springer-Verlag.

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   VERSION 2.2
!
!   f90 VERSION
!
!   This library contains routines for B-spline interpolation in
!   one, two, and three dimensions. Part of the routines are based
!   on the book by Carl de Boor: A practical guide to Splines (Springer,
!   New-York 1978) and have the same calling sequence and names as
!   the corresponding routines from the IMSL library. For documen-
!   tation see the additional files. NOTE: The results in the demo
!   routines may vary slightly on different architectures.
!
!   by W. Schadow 12/04/99
!   last changed by W. Schadow 07/28/2000
!
!
!   Wolfgang Schadow
!   TRIUMF
!   4004 Wesbrook Mall
!   Vancouver, B.C. V6T 2A3
!   Canada
!
!   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
!
!   www  : http://www.triumf.ca/people/schadow
!
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   Copyright (C) 2000 Wolfgang Schadow
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Library General Public
!   License as published by the Free Software Foundation; either
!   version 2 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public
!   License along with this library; if not, write to the
!   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!   Boston, MA  02111-1307, USA.
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!> @brief Module defining numeric precision constants
!> @details Provides single and double precision kind parameters for consistent
!> precision handling across the B-spline library.
module numeric

    !> @var sgl Single precision kind parameter
    integer, parameter :: sgl = kind(1.0)
    !> @var dbl Double precision kind parameter  
    integer, parameter :: dbl = kind(1.0d0)

end module numeric

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!> @brief Main B-spline interpolation module
!> @details Provides comprehensive B-spline functionality for interpolation,
!> evaluation, and derivative computation in 1D, 2D, and 3D.
!>
!> The following routines are included:
!> - @ref dbsnak: Compute "not-a-knot" spline knot sequence
!> - @ref dbsint: Compute spline interpolant and B-spline coefficients
!> - @ref dbsval: Evaluate spline at a point
!> - @ref dbsder: Evaluate spline derivatives
!> - @ref dbs1gd: Evaluate spline derivatives on a grid
!> - @ref dbs2in: 2D tensor-product spline interpolation
!> - @ref dbs2dr: 2D spline derivative evaluation
!> - @ref dbs2vl: 2D spline evaluation
!> - @ref dbs2gd: 2D spline evaluation on a grid
!> - @ref dbs3in: 3D tensor-product spline interpolation
!> - @ref dbs3vl: 3D spline evaluation
!> - @ref dbs3dr: 3D spline derivative evaluation
!> - @ref dbs3gd: 3D spline evaluation on a grid
module bspline

    private

    public dbsnak
    public dbsint, dbsval, dbsder, dbs1gd
    public dbs2in, dbs2dr, dbs2vl, dbs2gd
    public dbs3in, dbs3vl, dbs3dr, dbs3gd

contains

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Compute the "not-a-knot" spline knot sequence
    !> @details Generates a knot sequence for B-spline interpolation using the
    !> "not-a-knot" end condition. This condition ensures smooth interpolation
    !> by setting the first and last k knots to be equal, where k is the spline order.
    !> 
    !> The algorithm follows de Boor's method (p. 167) for constructing
    !> appropriate knot sequences that avoid oscillations at the boundaries.
    !>
    !> @param[in] nx Number of data points
    !> @param[in] xvec Array of length nx containing the location of data points
    !> @param[in] kxord Order of the spline (must be <= nx)
    !> @param[out] xknot Array of length nx+kxord containing the knot sequence
    !>
    !> @note The knot sequence is non-decreasing and follows the "not-a-knot" condition
    !> @warning kxord must satisfy 0 <= kxord <= nx
    !> @see de Boor, C. (1978). A practical guide to Splines. Springer-Verlag.
    subroutine dbsnak(nx, xvec, kxord, xknot)

        use numeric

        implicit none

        integer, intent(in) :: nx, kxord

        real(kind=dbl), dimension(nx), intent(in) :: xvec
        real(kind=dbl), dimension(nx + kxord), intent(out) :: xknot

        real(kind=dbl) :: eps
        integer :: ix
        logical :: first = .true.

        save first, eps

        ! Initialize epsilon for numerical stability
        if (first) then
            first = .false.
            eps = epsilon(1.0_dbl)
            !write(6,*) "subroutine dbsnak: "
            !write(6,*) "eps = ",eps
        end if

        ! Validate input parameters
        if ((kxord .lt. 0) .or. (kxord .gt. nx)) then
            write (6, *) "subroutine dbsnak: error"
            write (6, *) "0 <= kxord <= nx is required."
            write (6, *) "kxord = ", kxord, " and nx = ", nx, " is given."
            stop
        end if

        ! Set first k knots to the leftmost data point
        do ix = 1, kxord
            xknot(ix) = xvec(1)
        end do

        ! Set interior knots based on spline order parity
        if (mod(kxord, 2) .eq. 0) then
            ! Even order: use data points directly
            do ix = kxord + 1, nx
                xknot(ix) = xvec(ix - kxord/2)
            end do
        else
            ! Odd order: use midpoints between data points
            do ix = kxord + 1, nx
                xknot(ix) = 0.5_dbl*(xvec(ix - kxord/2) + xvec(ix - kxord/2 - 1))
            end do
        end if

        ! Set last k knots to the rightmost data point with small offset
        do ix = nx + 1, nx + kxord
            xknot(ix) = xvec(nx)*(1.0_dbl + eps)
        end do

    end subroutine dbsnak

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Compute the spline interpolant, returning the B-spline coefficients
    !> @details Computes B-spline coefficients that interpolate the given data points.
    !> The algorithm constructs a linear system and solves it to find the coefficients
    !> that ensure the spline passes through all data points.
    !>
    !> The method uses the B-spline basis functions to construct a banded linear system
    !> which is then solved using LU decomposition for banded matrices.
    !>
    !> @param[in] nx Number of data points
    !> @param[in] xvec Array of length nx containing the data point abscissas
    !> @param[in] xdata Array of length nx containing the data point ordinates
    !> @param[in] kx Order of the spline (must be <= nx)
    !> @param[in] xknot Array of length nx+kx containing the knot sequence (non-decreasing)
    !> @param[out] bcoef Array of length nx containing the B-spline coefficients
    !>
    !> @note The knot sequence must be non-decreasing
    !> @warning The linear system must be solvable (iflag = 1)
    !> @see de Boor, C. (1978). A practical guide to Splines. Springer-Verlag.
    subroutine dbsint(nx, xvec, xdata, kx, xknot, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: nx, kx
        real(kind=dbl), dimension(nx), intent(in) :: xdata, xvec
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(nx), intent(out) :: bcoef

        integer :: nxp1, kxm1, kpkm2, leftx, lenq
        integer :: ix, ik, ilp1mx, jj, iflag
        real(kind=dbl) :: xveci
        real(kind=dbl), dimension((2*kx - 1)*nx) :: work

        ! Initialize parameters for banded matrix construction
        nxp1 = nx + 1
        kxm1 = kx - 1
        kpkm2 = 2*kxm1
        leftx = kx
        lenq = nx*(kx + kxm1)

        ! Initialize work array
        do ix = 1, lenq
            work(ix) = 0.0_dbl
        end do

        ! Construct the banded matrix for the linear system
        do ix = 1, nx
            xveci = xvec(ix)
            ilp1mx = min0(ix + kx, nxp1)
            leftx = max0(leftx, ix)
            if (xveci .lt. xknot(leftx)) goto 998
30          if (xveci .lt. xknot(leftx + 1)) go to 40
            leftx = leftx + 1
            if (leftx .lt. ilp1mx) go to 30
            leftx = leftx - 1
            if (xveci .gt. xknot(leftx + 1)) goto 998
40          call bsplvb(xknot, nx + kx, kx, 1, xveci, leftx, bcoef)
            jj = ix - leftx + 1 + (leftx - kx)*(kx + kxm1)
            do ik = 1, kx
                jj = jj + kpkm2
                work(jj) = bcoef(ik)
            end do
        end do

        ! Factor the banded matrix
        call banfac(work, kx + kxm1, nx, kxm1, kxm1, iflag)

        ! Check if factorization was successful
        if (iflag .ne. 1) then
            write (6, *) "subroutine dbsint: error"
            write (6, *) "no solution of linear equation system !!!"
            stop
        end if

        ! Copy data to coefficient array
        do ix = 1, nx
            bcoef(ix) = xdata(ix)
        end do

        ! Solve the linear system
        call banslv(work, kx + kxm1, nx, kxm1, kxm1, bcoef)

        return

998     write (6, *) "subroutine dbsint:"
        write (6, *) "xknot(ix) <= xknot(ix+1) required."
        write (6, *) ix, xknot(ix), xknot(ix + 1)

        stop

    end subroutine dbsint

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Evaluate a spline at a given point
    !> @details Evaluates a B-spline at a specified point using the de Boor algorithm.
    !> The function finds the appropriate knot interval and computes the spline value
    !> using the B-spline basis functions and coefficients.
    !>
    !> The evaluation uses the recursive de Boor algorithm which is numerically
    !> stable and efficient for B-spline evaluation.
    !>
    !> @param[in] x Point at which the spline is to be evaluated
    !> @param[in] kx Order of the spline
    !> @param[in] xknot Array of length nx+kx containing the knot sequence (non-decreasing)
    !> @param[in] nx Number of B-spline coefficients
    !> @param[in] bcoef Array of length nx containing the B-spline coefficients
    !> @return Value of the spline at x
    !>
    !> @note The knot sequence must be non-decreasing
    !> @warning x must be within the knot range [xknot(1), xknot(nx+kx)]
    !> @see de Boor, C. (1978). A practical guide to Splines. Springer-Verlag.
    function dbsval(x, kx, xknot, nx, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: nx, kx
        real(kind=dbl) :: dbsval
        real(kind=dbl) :: x
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(nx), intent(in) :: bcoef

        integer :: il, ik, ix, leftx
        real(kind=dbl) :: save1, save2
        real(kind=dbl), dimension(kx) :: work, dl, dr

        ! Check knot sequence and find appropriate interval
        leftx = 0

        do ix = 1, nx + kx - 1
            if (xknot(ix) .gt. xknot(ix + 1)) then
                write (6, *) "subroutine dbsval:"
                write (6, *) "xknot(ix) <= xknot(ix+1) required."
                write (6, *) ix, xknot(ix), xknot(ix + 1)
                stop
            end if
            if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix + 1))) leftx = ix
        end do

        if (leftx .eq. 0) then
            write (6, *) "subroutine dbsval:"
            write (6, *) "ix with xknot(ix) <= x < xknot(ix+1) required."
            write (6, *) "x = ", x
            stop
        end if

        ! Initialize work arrays for de Boor algorithm
        do ik = 1, kx - 1
            work(ik) = bcoef(leftx + ik - kx)
            dl(ik) = x - xknot(leftx + ik - kx)
            dr(ik) = xknot(leftx + ik) - x
        end do

        work(kx) = bcoef(leftx)
        dl(kx) = x - xknot(leftx)

        ! Apply de Boor algorithm recursively
        do ik = 1, kx - 1
            save2 = work(ik)
            do il = ik + 1, kx
                save1 = work(il)
                work(il) = (dl(il)*work(il) + dr(il - ik)*save2)                  &
                        & /(dl(il) + dr(il - ik))
                save2 = save1
            end do
        end do

        dbsval = work(kx)

    end function dbsval

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Evaluate the derivative of a spline at a given point
    !> @details Evaluates the derivative of a B-spline at a specified point.
    !> The function can compute derivatives of any order up to (spline order - 1).
    !> For iderx = 0, it returns the spline value itself.
    !>
    !> The derivative computation uses the relationship between B-spline derivatives
    !> and lower-order B-splines, combined with the de Boor algorithm.
    !>
    !> @param[in] iderx Order of the derivative to be evaluated (0 = function value)
    !> @param[in] x Point at which the spline is to be evaluated
    !> @param[in] kx Order of the spline
    !> @param[in] xknot Array of length nx+kx containing the knot sequence (non-decreasing)
    !> @param[in] nx Number of B-spline coefficients
    !> @param[in] bcoef Array of length nx containing the B-spline coefficients
    !> @return Value of the iderx-th derivative of the spline at x
    !>
    !> @note For iderx >= kx, the result is zero
    !> @warning x must be within the knot range [xknot(1), xknot(nx+kx)]
    !> @see de Boor, C. (1978). A practical guide to Splines. Springer-Verlag.
    function dbsder(iderx, x, kx, xknot, nx, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: iderx, kx, nx
        real(kind=dbl) :: dbsder
        real(kind=dbl), intent(in) :: x
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(nx), intent(in) :: bcoef

        integer :: ix, ik, il, leftx
        real(kind=dbl) :: save, save1, save2, y, sum, dik
        real(kind=dbl), dimension(kx) :: work, dl, dr, bsp

        ! Check knot sequence and find appropriate interval
        leftx = 0
        do ix = 1, nx + kx - 1
            if (xknot(ix) .gt. xknot(ix + 1)) then
                write (6, *) "subroutine dbsder:"
                write (6, *) "xknot(ix) <= xknot(ix+1) required."
                stop
            end if
            if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix + 1))) leftx = ix
        end do

        if (leftx .eq. 0) then
            write (6, *) "subroutine dbsder:"
            write (6, *) "ix with xknot(ix) <= x < xknot(ix+1) required."
            write (6, *) "xknot(1)     = ", xknot(1)
            write (6, *) "xknot(nx+kx) = ", xknot(nx + kx)
            write (6, *) "         x   = ", x
            stop
        end if

        ! Handle function value (iderx = 0)
        if (iderx .eq. 0) then

            do ik = 1, kx - 1
                work(ik) = bcoef(leftx + ik - kx)
                dl(ik) = x - xknot(leftx + ik - kx)
                dr(ik) = xknot(leftx + ik) - x
            end do

            work(kx) = bcoef(leftx)
            dl(kx) = x - xknot(leftx)

            do ik = 1, kx - 1
                save2 = work(ik)
                do il = ik + 1, kx
                    save1 = work(il)
                    work(il) = (dl(il)*work(il) + dr(il - ik)*save2)               &
                            & /(dl(il) + dr(il - ik))
                    save2 = save1
                end do
            end do

            dbsder = work(kx)

        ! Handle derivative computation (1 <= iderx < kx)
        elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

            ! Compute B-spline basis functions for derivative
            bsp(1) = 1.0_dbl
            do ik = 1, kx - iderx - 1
                dr(ik) = xknot(leftx + ik) - x
                dl(ik) = x - xknot(leftx + 1 - ik)
                save = bsp(1)
                bsp(1) = 0.0_dbl
                do il = 1, ik
                    y = save/(dr(il) + dl(ik + 1 - il))
                    bsp(il) = bsp(il) + dr(il)*y
                    save = bsp(il + 1)
                    bsp(il + 1) = dl(ik + 1 - il)*y
                end do
            end do

            ! Compute derivative coefficients
            do ik = 1, kx
                work(ik) = bcoef(leftx + ik - kx)
                dr(ik) = xknot(leftx + ik) - x
                dl(ik) = x - xknot(leftx + ik - kx)
            end do

            ! Apply derivative operator
            do ik = 1, iderx
                dik = dble(kx - ik)
                save2 = work(ik)
                do il = ik + 1, kx
                    save1 = work(il)
                    work(il) = dik*(work(il) - save2)/(dl(il) + dr(il - ik))
                    save2 = save1
                end do
            end do

            ! Compute final derivative value
            sum = 0.0_dbl
            do ix = 1, kx - iderx
                sum = sum + bsp(ix)*work(iderx + ix)
            end do

            dbsder = sum

        else
            ! Higher order derivatives are zero
            dbsder = 0.0_dbl
        end if

    end function dbsder

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Evaluate the derivative of a spline on a grid
    !> @details Evaluates the derivative of a B-spline at multiple points efficiently.
    !> This routine is optimized for evaluating splines on regular grids by
    !> reusing knot interval searches and basis function computations.
    !>
    !> The grid evaluation uses the hunt algorithm to efficiently locate knot
    !> intervals for consecutive grid points, significantly improving performance
    !> compared to individual point evaluation.
    !>
    !> @param[in] iderx Order of the derivative to be evaluated (0 = function value)
    !> @param[in] nxvec Length of vector xvec
    !> @param[in] xvec Array of length nxvec containing evaluation points (strictly increasing)
    !> @param[in] kx Order of the spline
    !> @param[in] xknot Array of length nx+kx containing the knot sequence (non-decreasing)
    !> @param[in] nx Number of B-spline coefficients
    !> @param[in] bcoef Array of length nx containing the B-spline coefficients
    !> @param[out] val Array of length nxvec containing the derivative values
    !>
    !> @note xvec should be strictly increasing for optimal performance
    !> @warning All points in xvec must be within the knot range
    !> @see huntn subroutine for efficient knot interval searching
    subroutine dbs1gd(iderx, nxvec, xvec, kx, xknot, nx, bcoef, val)

        use numeric

        implicit none

        integer, intent(in) :: iderx, nxvec, kx, nx
        real(kind=dbl), dimension(nxvec), intent(in) :: xvec
        real(kind=dbl), dimension(nx), intent(in) :: bcoef
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(nxvec), intent(out) :: val

        integer :: i, il, ik, ix
        integer, dimension(nxvec) :: leftx
        real(kind=dbl) :: dik
        real(kind=dbl), dimension(nxvec, kx) :: dl, dr, biatx, work
        real(kind=dbl), dimension(nxvec) :: save1, save2, term

        logical :: same, next

        ! Initialize knot interval search
        leftx(1) = 0
        call huntn(xknot, nx + kx, kx, xvec(1), leftx(1))

        ! Efficiently find knot intervals for all grid points
        do ix = 2, nxvec
            leftx(ix) = leftx(ix - 1)
            same = (xknot(leftx(ix)) .le. xvec(ix))                                &
                    &        .and. (xvec(ix) .le. xknot(leftx(ix) + 1))
            if (.not. same) then
                leftx(ix) = leftx(ix) + 1
                next = (xknot(leftx(ix)) .le. xvec(ix))                        &
                        &           .and. (xvec(ix) .le. xknot(leftx(ix) + 1))
                if (.not. next)                                                     &
                        &           call huntn(xknot, nx + kx, kx, xvec(ix), leftx(ix))
            end if
        end do

        ! Validate knot sequence
        do ix = 1, nx + kx - 1
            if (xknot(ix) .gt. xknot(ix + 1)) then
                write (6, *) "subroutine dbs1gd:"
                write (6, *) "xknot(ix) <= xknot(ix+1) required."
                write (6, *) ix, xknot(ix), xknot(ix + 1)
                write (6, *)
                write (6, *) xknot
                stop
            end if
        end do

        ! Validate evaluation points
        do ix = 1, nxvec
            if ((xvec(ix) .lt. xknot(1)) .or. (xvec(ix) .gt. xknot(nx + kx))) then
                write (6, *) "subroutine dbs1gd:"
                write (6, *) "ix with xknot(ix) <= x < xknot(ix+1) required."
                write (6, *) "x = ", xvec(ix)
                stop
            end if
        end do

        ! Handle function value evaluation (iderx = 0)
        if (iderx .eq. 0) then

            ! Initialize basis functions
            do ix = 1, nxvec
                biatx(ix, 1) = 1._dbl
                val(ix) = 0._dbl
            end do

            ! Compute B-spline basis functions for all grid points
            do ik = 1, kx - 1
                do ix = 1, nxvec
                    dr(ix, ik) = xknot(leftx(ix) + ik) - xvec(ix)
                    dl(ix, ik) = xvec(ix) - xknot(leftx(ix) + 1 - ik)
                    save1(ix) = 0._dbl
                end do

                do il = 1, ik
                    do ix = 1, nxvec
                        term(ix) = biatx(ix, il)                                   &
                                & /(dr(ix, il) + dl(ix, ik + 1 - il))
                        biatx(ix, il) = save1(ix) + dr(ix, il)*term(ix)
                        save1(ix) = dl(ix, ik + 1 - il)*term(ix)
                    end do
                end do

                do ix = 1, nxvec
                    biatx(ix, ik + 1) = save1(ix)
                end do
            end do

            ! Compute spline values using basis functions and coefficients
            do ik = 1, kx
                do ix = 1, nxvec
                    val(ix) = val(ix) + biatx(ix, ik)*bcoef(leftx(ix) - kx + ik)
                end do
            end do

        ! Handle derivative evaluation (1 <= iderx < kx)
        elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

            ! Initialize basis functions for derivative
            do ix = 1, nxvec
                biatx(ix, 1) = 1._dbl
                val(ix) = 0._dbl
            end do

            ! Compute B-spline basis functions for derivative
            do ik = 1, kx - iderx - 1
                do ix = 1, nxvec
                    dr(ix, ik) = xknot(leftx(ix) + ik) - xvec(ix)
                    dl(ix, ik) = xvec(ix) - xknot(leftx(ix) + 1 - ik)
                    save1(ix) = biatx(ix, 1)
                    biatx(ix, 1) = 0.0_dbl
                    do il = 1, ik
                        term(ix) = save1(ix)                                    &
                                & /(dr(ix, il) + dl(ix, ik + 1 - il))
                        biatx(ix, il) = biatx(ix, il) + dr(ix, il)*term(ix)
                        save1(ix) = biatx(ix, il + 1)
                        biatx(ix, il + 1) = dl(ix, ik + 1 - il)*term(ix)
                    end do
                end do
            end do

            ! Compute derivative coefficients
            do ik = 1, kx
                do ix = 1, nxvec
                    work(ix, ik) = bcoef(leftx(ix) + ik - kx)
                    dr(ix, ik) = xknot(leftx(ix) + ik) - xvec(ix)
                    dl(ix, ik) = xvec(ix) - xknot(leftx(ix) + ik - kx)
                end do
            end do

            ! Apply derivative operator
            do ik = 1, iderx
                dik = dble(kx - ik)
                do ix = 1, nxvec
                    save2(ix) = work(ix, ik)
                    do il = ik + 1, kx
                        save1(ix) = work(ix, il)
                        work(ix, il) = dik*(work(ix, il) - save2(ix))                 &
                                & /(dl(ix, il) + dr(ix, il - ik))
                        save2(ix) = save1(ix)
                    end do
                end do
            end do

            ! Compute final derivative values
            do i = 1, kx - iderx
                do ix = 1, nxvec
                    val(ix) = val(ix) + biatx(ix, i)*work(ix, iderx + i)
                end do
            end do

        else
            ! Higher order derivatives are zero
            do ix = 1, nxvec
                val(ix) = 0.0_dbl
            end do

        end if

    end subroutine dbs1gd

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dbsdca(iderx, x, kx, xknot, nx, bcoef, leftx)

        !
        ! This routine is equivalent to the routine dbsder, but it does not
        ! check the parameters!!!
        !
        ! Evaluates the derivative of a spline, given its B-spline representation.
        !
        !
        !   iderx  - order of the derivative to be evaluated.  (input)
        !            in particular, iderx = 0 returns the value of the
        !            spline.
        !   x      - point at which the spline is to be evaluated.  (input)
        !   kx     - order of the spline.  (input)
        !   xknot  - array of length nx+kx containing the knot
        !            sequence.  (input)
        !            xknot must be nondecreasing.
        !   nx     - number of B-spline coefficients.  (input)
        !   bcoef  - array of length nx containing the B-spline
        !            coefficients.  (input)
        !   leftx  - number of the intervall of xknot that includes x
        !   dbsdca - value of the ideriv-th derivative of the spline at x.
        !            (output)
        !

        use numeric

        implicit none

        integer, intent(in) :: iderx, kx, nx
        real(kind=dbl) :: dbsdca
        real(kind=dbl), intent(in) :: x
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(nx), intent(in) :: bcoef

        integer :: i, ik, il, leftx
        real(kind=dbl) :: save, save1, save2, y, sum, dik
        real(kind=dbl), dimension(kx) :: work, dl, dr, bsp

        if (iderx .eq. 0) then

            do ik = 1, kx - 1
                work(ik) = bcoef(leftx + ik - kx)
                dl(ik) = x - xknot(leftx + ik - kx)
                dr(ik) = xknot(leftx + ik) - x
            end do

            work(kx) = bcoef(leftx)
            dl(kx) = x - xknot(leftx)

            do ik = 1, kx - 1
                save2 = work(ik)
                do il = ik + 1, kx
                    save1 = work(il)
                    work(il) = (dl(il)*work(il) + dr(il - ik)*save2)               &
                            & /(dl(il) + dr(il - ik))
                    save2 = save1
                end do
            end do

            dbsdca = work(kx)

        elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
            bsp(1) = 1.0_dbl
            do ik = 1, kx - iderx - 1
                dr(ik) = xknot(leftx + ik) - x
                dl(ik) = x - xknot(leftx + 1 - ik)
                save = bsp(1)
                bsp(1) = 0.0_dbl
                do il = 1, ik
                    y = save/(dr(il) + dl(ik + 1 - il))
                    bsp(il) = bsp(il) + dr(il)*y
                    save = bsp(il + 1)
                    bsp(il + 1) = dl(ik + 1 - il)*y
                end do
            end do

            do ik = 1, kx
                work(ik) = bcoef(leftx + ik - kx)
                dr(ik) = xknot(leftx + ik) - x
                dl(ik) = x - xknot(leftx + ik - kx)
            end do

            do ik = 1, iderx
                dik = dble(kx - ik)
                save2 = work(ik)
                do il = ik + 1, kx
                    save1 = work(il)
                    work(il) = dik*(work(il) - save2)/(dl(il) + dr(il - ik))
                    save2 = save1
                end do
            end do

            sum = 0.0_dbl

            do i = 1, kx - iderx
                sum = sum + bsp(i)*work(iderx + i)
            end do

            dbsdca = sum

        else
            dbsdca = 0.0_dbl
        end if

    end function dbsdca

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Compute a two-dimensional tensor-product spline interpolant
    !> @details Computes B-spline coefficients for 2D tensor-product spline interpolation.
    !> The algorithm constructs a 2D spline that interpolates data on a rectangular grid
    !> by applying 1D spline interpolation sequentially in each direction.
    !>
    !> The tensor-product approach allows efficient construction of 2D splines
    !> by combining 1D B-splines in the x and y directions.
    !>
    !> @param[in] nx Number of data points in the x-direction
    !> @param[in] xvec Array of length nx containing x-direction data points (strictly increasing)
    !> @param[in] ny Number of data points in the y-direction
    !> @param[in] yvec Array of length ny containing y-direction data points (strictly increasing)
    !> @param[in] xydata Array of size nx by ny containing values to be interpolated
    !> @param[in] ldf Leading dimension of xydata as specified in calling program
    !> @param[in] kx Order of the spline in the x-direction (must be <= nx)
    !> @param[in] ky Order of the spline in the y-direction (must be <= ny)
    !> @param[in] xknot Array of length nx+kx containing x-direction knot sequence (non-decreasing)
    !> @param[in] yknot Array of length ny+ky containing y-direction knot sequence (non-decreasing)
    !> @param[out] bcoef Array of length nx*ny containing tensor-product B-spline coefficients
    !>
    !> @note Both xvec and yvec must be strictly increasing
    !> @warning All knot sequences must be non-decreasing
    !> @see dbsint for 1D spline interpolation details
    subroutine dbs2in(nx, xvec, ny, yvec, xydata, ldf, kx, ky, xknot, yknot, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: nx, ny, kx, ky, ldf

        real(kind=dbl), dimension(nx), intent(in) :: xvec
        real(kind=dbl), dimension(ny), intent(in) :: yvec
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(ny + ky), intent(in) :: yknot
        real(kind=dbl), dimension(ldf, *), intent(in) :: xydata
        real(kind=dbl), dimension(nx, ny), intent(out) :: bcoef

        real(kind=dbl), dimension(max(nx, ny), max(nx, ny)) :: work1
        real(kind=dbl), dimension(max(nx, ny)) :: work2
        real(kind=dbl), dimension(max((2*kx - 1)*nx, (2*ky - 1)*ny)) :: work3

        ! First interpolate in x-direction for each y value
        call spli2d(xvec, ldf, xydata, xknot, nx, kx, ny, work2, work3, work1)
        ! Then interpolate in y-direction using the x-interpolated results
        call spli2d(yvec, ny, work1, yknot, ny, ky, nx, work2, work3, bcoef)

    end subroutine dbs2in

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Helper subroutine for 2D spline interpolation
    !> @details Performs 1D spline interpolation along one direction for 2D data.
    !> This is an internal subroutine used by dbs2in to construct tensor-product splines.
    !>
    !> @param[in] xyvec Array containing data points in the interpolation direction
    !> @param[in] ld Leading dimension of the data array
    !> @param[in] xydata 2D data array to be interpolated
    !> @param[in] xyzknot Knot sequence for the interpolation direction
    !> @param[in] n Number of data points in the interpolation direction
    !> @param[in] k Order of the spline in the interpolation direction
    !> @param[in] m Number of data points in the other direction
    !> @param[out] work2 Work array for temporary storage
    !> @param[out] work3 Work array for banded matrix operations
    !> @param[out] bcoef Output B-spline coefficients
    subroutine spli2d(xyvec, ld, xydata, xyzknot, n, k, m, work2, work3, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: ld, n, k, m
        real(kind=dbl), dimension(n), intent(in) :: xyvec
        real(kind=dbl), dimension(n + k), intent(in) :: xyzknot
        real(kind=dbl), dimension(ld, m), intent(in) :: xydata
        real(kind=dbl), dimension(m, n), intent(out) :: bcoef

        real(kind=dbl), dimension(n), intent(out) :: work2
        real(kind=dbl), dimension((2*k - 1)*n), intent(out) :: work3

        integer :: np1, km1, kpkm2, left, lenq, i, iflag, ilp1mx, j, jj
        real(kind=dbl) :: xyveci

        ! Initialize parameters for banded matrix construction
        np1 = n + 1
        km1 = k - 1
        kpkm2 = 2*km1
        left = k
        lenq = n*(k + km1)

        ! Initialize work array
        do i = 1, lenq
            work3(i) = 0.0_dbl
        end do

        ! Construct the banded matrix for the linear system
        do i = 1, n
            xyveci = xyvec(i)
            ilp1mx = min0(i + k, np1)
            left = max0(left, i)
            if (xyveci .lt. xyzknot(left)) go to 998
30          if (xyveci .lt. xyzknot(left + 1)) go to 40
            left = left + 1
            if (left .lt. ilp1mx) go to 30
            left = left - 1
            if (xyveci .gt. xyzknot(left + 1)) go to 998
40          call bsplvb(xyzknot, n + k, k, 1, xyveci, left, work2)
            jj = i - left + 1 + (left - k)*(k + km1)
            do j = 1, k
                jj = jj + kpkm2
                work3(jj) = work2(j)
            end do
        end do

        ! Factor the banded matrix
        call banfac(work3, k + km1, n, km1, km1, iflag)

        ! Check if factorization was successful
        if (iflag .ne. 1) then
            write (6, *) "subroutine dbs2in: error"
            write (6, *) "no solution of linear equation system !!!"
            stop
        end if

        ! Solve the linear system for each column
        do j = 1, m
            do i = 1, n
                work2(i) = xydata(i, j)
            end do

            call banslv(work3, k + km1, n, km1, km1, work2)

            do i = 1, n
                bcoef(j, i) = work2(i)
            end do
        end do

        return

998     write (6, *) "subroutine db2in:"
        write (6, *) "i with knot(i) <= x/y < knot(i+1) required."
        write (6, *) "knot(1)   = ", xyzknot(1)
        write (6, *) "knot(n+k) = ", xyzknot(n + k)
        write (6, *) "      x/y = ", xyveci

        stop

    end subroutine spli2d

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Evaluate a two-dimensional tensor-product spline at a point
    !> @details Evaluates a 2D tensor-product B-spline at a specified point (x,y).
    !> The evaluation uses the tensor-product structure to efficiently compute
    !> the spline value by evaluating 1D B-splines in each direction.
    !>
    !> @param[in] x x-coordinate of the evaluation point
    !> @param[in] y y-coordinate of the evaluation point
    !> @param[in] kx Order of the spline in the x-direction
    !> @param[in] ky Order of the spline in the y-direction
    !> @param[in] xknot Array of length nx+kx containing x-direction knot sequence (non-decreasing)
    !> @param[in] yknot Array of length ny+ky containing y-direction knot sequence (non-decreasing)
    !> @param[in] nx Number of B-spline coefficients in the x-direction
    !> @param[in] ny Number of B-spline coefficients in the y-direction
    !> @param[in] bcoef Array of length nx*ny containing tensor-product B-spline coefficients
    !> @return Value of the spline at (x,y)
    !>
    !> @note Both knot sequences must be non-decreasing
    !> @warning (x,y) must be within the knot ranges
    !> @see dbsval for 1D spline evaluation details
    function dbs2vl(x, y, kx, ky, xknot, yknot, nx, ny, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: nx, ny, kx, ky
        real(kind=dbl), intent(in) :: x, y
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(ny + ky), intent(in) :: yknot
        real(kind=dbl), dimension(nx, ny), intent(in) :: bcoef
        real(kind=dbl) :: dbs2vl

        integer :: ix, iy, iky, leftx, lefty
        real(kind=dbl), dimension(ky) :: work

        ! Check x-direction knot sequence and find interval
        leftx = 0

        do ix = 1, nx + kx - 1
            if (xknot(ix) .gt. xknot(ix + 1)) then
                write (6, *) "subroutine dbs2vl:"
                write (6, *) "xknot(ix) <= xknot(ix+1) required."
                write (6, *) ix, xknot(ix), xknot(ix + 1)
                write (6, *)
                write (6, *) xknot
                stop
            end if
            if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix + 1))) leftx = ix
        end do

        if (leftx .eq. 0) then
            write (6, *) "subroutine dbs2vl:"
            write (6, *) "ix with xknot(ix) <= x < xknot(ix+1) required."
            write (6, *) "x = ", x
            write (6, *)
            write (6, *) xknot
            stop
        end if

        ! Check y-direction knot sequence and find interval
        lefty = 0

        do iy = 1, ny + ky - 1
            if (yknot(iy) .gt. yknot(iy + 1)) then
                write (6, *) "subroutine dbs2vl:"
                write (6, *) "yknot(iy) <= yknot(iy+1) required."
                write (6, *) iy, yknot(iy), yknot(iy + 1)
                stop
            end if
            if ((yknot(iy) .le. y) .and. (y .lt. yknot(iy + 1))) lefty = iy
        end do

        if (lefty .eq. 0) then
            write (6, *) "subroutine dbs2vl:"
            write (6, *) "iy with yknot(iy) <= y < yknot(iy+1) required."
            write (6, *) "yknot(iy)   = ", yknot(iy)
            write (6, *) "  y         = ", y
            write (6, *) "yknot(iy+1) = ", yknot(iy + 1)
            stop
        end if

        ! Evaluate 1D splines in x-direction for each y basis function
        do iky = 1, ky
            work(iky) = dbsdca(0, x, kx, xknot, nx, bcoef(1, lefty - ky + iky), leftx)
        end do

        ! Evaluate the final result using y-direction spline
        dbs2vl = dbsval(y, ky, yknot(lefty - ky + 1), ky, work)

    end function dbs2vl

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dbs2dr(iderx, idery, x, y, kx, ky, xknot, yknot, nx, ny, bcoef)

        !
        !  Evaluates the derivative of a two-dimensional tensor-product spline,
        !  given its tensor-product B-spline representation.
        !
        !   iderx  - order of the derivative in the x-direction.  (input)
        !   idery  - order of the derivative in the y-direction.  (input)
        !   x      - x-coordinate of the point at which the spline is to be
        !            evaluated.  (input)
        !   y      - y-coordinate of the point at which the spline is to be
        !            evaluated.  (input)
        !   kx     - order of the spline in the x-direction.  (input)
        !   ky     - order of the spline in the y-direction.  (input)
        !   xknot  - array of length nx+kx containing the knot
        !            sequence in the x-direction.  (input)
        !            xknot must be nondecreasing.
        !   yknot  - array of length ny+ky containing the knot
        !            sequence in the y-direction.  (input)
        !            yknot must be nondecreasing.
        !   nx     - number of B-spline coefficients in the x-direction.
        !            (input)
        !   ny     - number of B-spline coefficients in the y-direction.
        !            (input)
        !   bcoef  - array of length nx*ny containing the
        !            tensor-product B-spline coefficients.  (input)
        !            bscoef is treated internally as a matrix of size nx
        !            by ny.
        !   dbs2dr  - value of the (iderx,idery) derivative of the spline at
        !            (x,y).  (output)
        !

        use numeric

        implicit none

        integer, intent(in) :: iderx, idery
        integer, intent(in) :: kx, nx, ky, ny
        real(kind=dbl) :: dbs2dr
        real(kind=dbl), intent(in) :: x, y
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(ny + ky), intent(in) :: yknot
        real(kind=dbl), dimension(nx, ny), intent(in) :: bcoef

        integer :: ix, iy, iky, nintx, ninty
        real(kind=dbl), dimension(ky) :: work

        !
        !     check if knot(i) <= knot(i+1) and calculation of i so that
        !     knot(i) <= x < knot(i+1)
        !

        nintx = 0

        do ix = 1, nx + kx - 1
            if (xknot(ix) .gt. xknot(ix + 1)) then
                write (6, *) "subroutine dbs2dr:"
                write (6, *) "xknot(ix) <= xknot(ix+1) required."
                write (6, *) ix, xknot(ix), xknot(ix + 1)
                stop
            end if
            if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix + 1))) nintx = ix
        end do

        if (nintx .eq. 0) then
            write (6, *) "subroutine dbs2dr:"
            write (6, *) "ix with xknot(ix) <= x < xknot(ix+1) required."
            write (6, *) "x = ", x
            stop
        end if

        ninty = 0

        do iy = 1, ny + ky - 1
            if (yknot(iy) .gt. yknot(iy + 1)) then
                write (6, *) "subroutine dbs2dr:"
                write (6, *) "yknot(iy) <= yknot(iy+1) required."
                write (6, *) iy, yknot(iy), yknot(iy + 1)
                stop
            end if
            if ((yknot(iy) .le. y) .and. (y .lt. yknot(iy + 1))) ninty = iy
        end do

        if (ninty .eq. 0) then
            write (6, *) "subroutine dbs2dr:"
            write (6, *) "iy with yknot(iy) <= y < yknot(iy+1) required."
            write (6, *) "y = ", y
            stop
        end if

        do iky = 1, ky
            work(iky) = dbsdca(iderx, x, kx, xknot, nx, bcoef(1, ninty - ky + iky), nintx)
        end do

        dbs2dr = dbsder(idery, y, ky, yknot(ninty - ky + 1), ky, work)

    end function dbs2dr

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Evaluate a two-dimensional tensor-product spline on a grid
    !> @details Evaluates the derivative of a 2D tensor-product spline at multiple points.
    !> This routine is optimized for evaluating splines on regular grids by
    !> reusing knot interval searches and basis function computations.
    !>
    !> The grid evaluation uses the hunt algorithm to efficiently locate knot
    !> intervals for consecutive grid points, significantly improving performance
    !> compared to individual point evaluation.
    !>
    !> @param[in] iderx Order of the x-derivative
    !> @param[in] idery Order of the y-derivative
    !> @param[in] nxvec Length of vector xvec (strictly increasing)
    !> @param[in] xvec Array of x-coordinates for evaluation
    !> @param[in] nyvec Length of vector yvec (strictly increasing)
    !> @param[in] yvec Array of y-coordinates for evaluation
    !> @param[in] kx Order of the spline in the x-direction
    !> @param[in] ky Order of the spline in the y-direction
    !> @param[in] xknot x-direction knot sequence (non-decreasing)
    !> @param[in] yknot y-direction knot sequence (non-decreasing)
    !> @param[in] nx Number of B-spline coefficients in the x-direction
    !> @param[in] ny Number of B-spline coefficients in the y-direction
    !> @param[in] bcoef Tensor-product B-spline coefficients (nx by ny)
    !> @param[out] val Array of derivative values on the grid
    !> @param[in] ldf Leading dimension of val
    !>
    !> @note Both xvec and yvec must be strictly increasing
    !> @warning All points must be within the knot ranges
    subroutine dbs3gd(iderx, idery, iderz, nxvec, xvec, nyvec, yvec, nzvec, zvec, kx, ky, kz, xknot, yknot, zknot, &
            & nx, ny, nz, bcoef, val, ldf, mdf)

        !
        !  Evaluates the derivative of a three-dimensional tensor-product spline,
        !  given its tensor-product B-spline representation on a grid.
        !
        !   iderx  - order of the x-derivative.  (input)
        !   idery  - order of the y-derivative.  (input)
        !   iderz  - order of the z-derivative.  (input)
        !   nxvec  - number of grid points in the x-direction.  (input)
        !   xvec   - array of length nx containing the x-coordinates at
        !            which the spline is to be evaluated.  (input)
        !            the points in xvec should be strictly increasing.
        !   nyvec  - number of grid points in the y-direction.  (input)
        !   yvec   - array of length ny containing the y-coordinates at
        !            which the spline is to be evaluated.  (input)
        !            the points in yvec should be strictly increasing.
        !   nzvec  - number of grid points in the z-direction.  (input)
        !   zvec   - array of length nz containing the z-coordinates at
        !            which the spline is to be evaluated.  (input)
        !            the points in zvec should be strictly increasing.
        !   kx     - order of the spline in the x-direction.  (input)
        !   ky     - order of the spline in the y-direction.  (input)
        !   kz     - order of the spline in the z-direction.  (input)
        !   xknot  - array of length nx+kx containing the knot
        !            sequence in the x-direction.  (input)
        !            xknot must be nondecreasing.
        !   yknot  - array of length ny+ky containing the knot
        !            sequence in the y-direction.  (input)
        !            yknot must be nondecreasing.
        !   zknot  - array of length nz+kz containing the knot
        !            sequence in the z-direction.  (input)
        !            zknot must be nondecreasing.
        !   nx     - number of B-spline coefficients in the x-direction.
        !            (input)
        !   ny     - number of B-spline coefficients in the y-direction.
        !            (input)
        !   nz     - number of B-spline coefficients in the z-direction.
        !            (input)
        !   bcoef  - array of length nx*ny*nz containing the
        !            tensor-product B-spline coefficients.  (input)
        !            bscoef is treated internally as a matrix of size nx
        !            by ny by nz.
        !   val    - array of size nx by ny by nz containing the values of
        !            the (iderx,idery,iderz) derivative of the spline on the
        !            nx by ny by nz grid.  (output)
        !            value(i,j,k) contains the derivative of the spline at the
        !            point (xvec(i),yvec(j),zvec(k)).
        !   ldf    - leading dimension of value exactly as specified in the
        !            dimension statement of the calling program.  (input)
        !   mdf    - middle dimension of value exactly as specified in the
        !            dimension statement of the calling program.  (input)
        !

        use numeric

        implicit none

        integer, intent(in) :: iderx, idery, iderz
        integer, intent(in) :: nxvec, nyvec, nzvec
        integer, intent(in) :: kx, nx, ky, ny, kz, nz
        integer, intent(in) :: ldf, mdf

        real(kind=dbl), dimension(nxvec), intent(in) :: xvec
        real(kind=dbl), dimension(nyvec), intent(in) :: yvec
        real(kind=dbl), dimension(nzvec), intent(in) :: zvec
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(ny + ky), intent(in) :: yknot
        real(kind=dbl), dimension(nz + kz), intent(in) :: zknot
        real(kind=dbl), dimension(nx, ny, nz), intent(in) :: bcoef
        real(kind=dbl), dimension(ldf, mdf, *), intent(out) :: val

        integer :: i, ik, il, ix, iy, iz, ikx, iky, ikz
        integer, dimension(nxvec) :: leftx
        integer, dimension(nyvec) :: lefty
        integer, dimension(nzvec) :: leftz
        real(kind=dbl), dimension(nxvec, kx) :: biatx
        real(kind=dbl), dimension(nyvec, ky) :: biaty
        real(kind=dbl), dimension(nzvec, kz) :: biatz
        real(kind=dbl), dimension(max(nxvec, nyvec, nzvec)) :: term, save1

        real(kind=dbl), dimension(max(nxvec, nyvec, nzvec), max(kx, ky, kz)) :: dl, dr

        logical :: same, next

        do i = 1, nx + kx - 1
            if (xknot(i) .gt. xknot(i + 1)) then
                write (6, *) "subroutine dbs3gd:"
                write (6, *) "xknot(i) <= xknot(i+1) required."
                write (6, *) i, xknot(i), xknot(i + 1)
                write (6, *)
                write (6, *) xknot
                stop
            end if
        end do

        do i = 1, nxvec
            if ((xvec(i) .lt. xknot(1)) .or. (xvec(i) .gt. xknot(nx + kx))) then
                write (6, *) "subroutine dbs3gd:"
                write (6, *) "ix with xknot(ix) <= x < xknot(ix+1) required."
                write (6, *) "x = ", xvec(i)
                stop
            end if
        end do

        leftx(1) = 0

        call huntn(xknot, nx + kx, kx, xvec(1), leftx(1))

        do ix = 2, nxvec
            leftx(ix) = leftx(ix - 1)
            same = (xknot(leftx(ix)) .le. xvec(ix))                                &
                    &        .and. (xvec(ix) .le. xknot(leftx(ix) + 1))
            if (.not. same) then
                leftx(ix) = leftx(ix) + 1
                next = (xknot(leftx(ix)) .le. xvec(ix))                        &
                        &           .and. (xvec(ix) .le. xknot(leftx(ix) + 1))
                if (.not. next) call huntn(xknot, nx + kx, kx, xvec(ix), leftx(ix))
            end if
        end do

        do i = 1, ny + ky - 1
            if (yknot(i) .gt. yknot(i + 1)) then
                write (6, *) "subroutine dbs3gd:"
                write (6, *) "yknot(i) <= yknot(i+1) required."
                write (6, *) i, yknot(i), yknot(i + 1)
                write (6, *)
                write (6, *) yknot
                stop
            end if
        end do

        do i = 1, nyvec
            if ((yvec(i) .lt. yknot(1)) .or. (yvec(i) .gt. yknot(ny + ky))) then
                write (6, *) "subroutine dbs3gd:"
                write (6, *) "iy with yknot(iy) <= y < yknot(iy+1) required."
                write (6, *) "y = ", yvec(i)
                stop
            end if
        end do

        lefty(1) = 0

        call huntn(yknot, ny + ky, ky, yvec(1), lefty(1))

        do iy = 2, nyvec
            lefty(iy) = lefty(iy - 1)
            same = (yknot(lefty(iy)) .le. yvec(iy))                                &
                    &        .and. (yvec(iy) .le. yknot(lefty(iy) + 1))
            if (.not. same) then
                lefty(iy) = lefty(iy) + 1
                next = (yknot(lefty(iy)) .le. yvec(iy))                        &
                        &           .and. (yvec(iy) .le. yknot(lefty(iy) + 1))
                if (.not. next) call huntn(yknot, ny + ky, ky, yvec(iy), lefty(iy))
            end if
        end do

        do i = 1, nz + kz - 1
            if (zknot(i) .gt. zknot(i + 1)) then
                write (6, *) "subroutine dbs3gd:"
                write (6, *) "zknot(i) <= zknot(i+1) required."
                write (6, *) i, zknot(i), zknot(i + 1)
                write (6, *)
                write (6, *) zknot(:)
                stop
            end if
        end do

        do i = 1, nzvec
            if ((zvec(i) .lt. zknot(1)) .or. (zvec(i) .gt. zknot(nz + kz))) then
                write (6, *) "subroutine dbs3gd:"
                write (6, *) "iz with zknot(iz) <= z < zknot(iz+1) required."
                write (6, *) "z = ", zvec(i)
                stop
            end if
        end do

        leftz(1) = 0

        call huntn(zknot, nz + kz, kz, zvec(1), leftz(1))

        do iz = 2, nzvec
            leftz(iz) = leftz(iz - 1)
            same = (zknot(leftz(iz)) .le. zvec(iz))                                &
                    &        .and. (zvec(iz) .le. zknot(leftz(iz) + 1))
            if (.not. same) then
                leftz(iz) = leftz(iz) + 1
                next = (zknot(leftz(iz)) .le. zvec(iz))                        &
                        &           .and. (zvec(iz) .le. zknot(leftz(iz) + 1))
                if (.not. next) call huntn(zknot, nz + kz, kz, zvec(iz), leftz(iz))
            end if
        end do
        ! by E. Coccia (4/1/11): evaluate the function
        if ((iderx .eq. 0) .and. (idery .eq. 0) .and. (iderz .eq. 0)) then

            do ix = 1, nxvec
                biatx(ix, 1) = 1.0_dbl
            end do

            do ik = 1, kx - 1
                do ix = 1, nxvec
                    dr(ix, ik) = xknot(leftx(ix) + ik) - xvec(ix)
                    dl(ix, ik) = xvec(ix) - xknot(leftx(ix) + 1 - ik)
                    save1(ix) = 0._dbl
                end do

                do il = 1, ik
                    do ix = 1, nxvec
                        term(ix) = biatx(ix, il)/(dr(ix, il) + dl(ix, ik + 1 - il))
                        biatx(ix, il) = save1(ix) + dr(ix, il)*term(ix)
                        save1(ix) = dl(ix, ik + 1 - il)*term(ix)
                    end do
                end do

                do ix = 1, nxvec
                    biatx(ix, ik + 1) = save1(ix)
                end do
            end do

            do iy = 1, nyvec
                biaty(iy, 1) = 1.0_dbl
            end do

            do ik = 1, ky - 1
                do iy = 1, nyvec
                    dr(iy, ik) = yknot(lefty(iy) + ik) - yvec(iy)
                    dl(iy, ik) = yvec(iy) - yknot(lefty(iy) + 1 - ik)
                    save1(iy) = 0._dbl
                end do

                do il = 1, ik
                    do iy = 1, nyvec
                        term(iy) = biaty(iy, il)/(dr(iy, il) + dl(iy, ik + 1 - il))
                        biaty(iy, il) = save1(iy) + dr(iy, il)*term(iy)
                        save1(iy) = dl(iy, ik + 1 - il)*term(iy)
                    end do
                end do

                do iy = 1, nyvec
                    biaty(iy, ik + 1) = save1(iy)
                end do
            end do

            do iz = 1, nzvec
                biatz(iz, 1) = 1.0_dbl
            end do

            do ik = 1, kz - 1
                do iz = 1, nzvec
                    dr(iz, ik) = zknot(leftz(iz) + ik) - zvec(iz)
                    dl(iz, ik) = zvec(iz) - zknot(leftz(iz) + 1 - ik)
                    save1(iz) = 0._dbl
                end do

                do il = 1, ik
                    do iz = 1, nzvec
                        term(iz) = biatz(iz, il)/(dr(iz, il) + dl(iz, ik + 1 - il))
                        biatz(iz, il) = save1(iz) + dr(iz, il)*term(iz)
                        save1(iz) = dl(iz, ik + 1 - il)*term(iz)
                    end do
                end do

                do iz = 1, nzvec
                    biatz(iz, ik + 1) = save1(iz)
                end do
            end do

            do iz = 1, nzvec
                do iy = 1, nyvec
                    do ix = 1, nxvec
                        val(ix, iy, iz) = 0.0_dbl
                    end do
                end do
            end do

            do ikz = 1, kz
                do iky = 1, ky
                    do ikx = 1, kx
                        do iz = 1, nzvec
                            do iy = 1, nyvec
                                do ix = 1, nxvec
                                    val(ix, iy, iz) = val(ix, iy, iz)                        &
                                            & + biatx(ix, ikx)*biaty(iy, iky)              &
                                                    & *biatz(iz, ikz)                              &
                                                    & *bcoef(leftx(ix) - kx + ikx, &
                                                            &          lefty(iy) - ky + iky, leftz(iz) - kz + ikz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            ! by E. Coccia (4/1/11): evaluate the derivatives
        else

            do iz = 1, nzvec
                do iy = 1, nyvec
                    do ix = 1, nxvec
                        val(ix, iy, iz) = dbs3dr(iderx, idery, iderz, xvec(ix), &
                                &  yvec(iy), zvec(iz), kx, ky, kz, xknot, yknot, &
                                &  zknot, nx, ny, nz, bcoef)
                    end do
                end do
            end do

        end if

    end subroutine dbs3gd

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Compute a three-dimensional tensor-product spline interpolant
    !> @details Computes B-spline coefficients for 3D tensor-product spline interpolation.
    !> The algorithm constructs a 3D spline that interpolates data on a rectangular grid
    !> by applying 1D spline interpolation sequentially in each direction.
    !>
    !> The tensor-product approach extends 2D splines to 3D by combining 1D B-splines
    !> in the x, y, and z directions.
    !>
    !> @param[in] nx Number of data points in the x-direction
    !> @param[in] xvec Array of length nx containing x-direction data points (increasing)
    !> @param[in] ny Number of data points in the y-direction
    !> @param[in] yvec Array of length ny containing y-direction data points (increasing)
    !> @param[in] nz Number of data points in the z-direction
    !> @param[in] zvec Array of length nz containing z-direction data points (increasing)
    !> @param[in] xyzdata Array of size nx by ny by nz containing values to be interpolated
    !> @param[in] ldf Leading dimension of xyzdata as specified in calling program
    !> @param[in] mdf Middle dimension of xyzdata as specified in calling program
    !> @param[in] kx Order of the spline in the x-direction (must be <= nx)
    !> @param[in] ky Order of the spline in the y-direction (must be <= ny)
    !> @param[in] kz Order of the spline in the z-direction (must be <= nz)
    !> @param[in] xknot Array of length nx+kx containing x-direction knot sequence (non-decreasing)
    !> @param[in] yknot Array of length ny+ky containing y-direction knot sequence (non-decreasing)
    !> @param[in] zknot Array of length nz+kz containing z-direction knot sequence (non-decreasing)
    !> @param[out] bcoef Array of length nx*ny*nz containing tensor-product B-spline coefficients
    !>
    !> @note All data point arrays must be increasing
    !> @warning All knot sequences must be non-decreasing
    !> @see dbs2in for 2D spline interpolation details
    subroutine dbs3in(nx, xvec, ny, yvec, nz, zvec, xyzdata, ldf, mdf, kx, ky, kz, &
            & xknot, yknot, zknot, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: nx, ny, nz, kx, ky, kz
        integer, intent(in) :: ldf, mdf

        real(kind=dbl), dimension(nx), intent(in) :: xvec
        real(kind=dbl), dimension(ny), intent(in) :: yvec
        real(kind=dbl), dimension(nz), intent(in) :: zvec
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(ny + ky), intent(in) :: yknot
        real(kind=dbl), dimension(nz + kz), intent(in) :: zknot
        real(kind=dbl), dimension(ldf, mdf, nz), intent(in) :: xyzdata
        real(kind=dbl), dimension(nx, ny, nz), intent(out) :: bcoef

        integer :: iz
        real(kind=dbl), dimension(nx, ny, nz) :: work1
        real(kind=dbl), dimension(nz) :: work2
        real(kind=dbl), dimension((2*kz - 1)*nz) :: work3

        ! First interpolate in z-direction for each (x,y) pair
        call spli3d(zvec, ldf, mdf, xyzdata, zknot, nz, kz, nx, ny, work2, work3, work1, &
                &     nx, ny, nz)

        ! Then interpolate in x and y directions using 2D interpolation
        do iz = 1, nz
            call dbs2in(nx, xvec, ny, yvec, work1(1, 1, iz), nx, kx, ky, xknot, yknot, &
                    &        bcoef(1, 1, iz))
        end do

    end subroutine dbs3in

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Helper subroutine for 3D spline interpolation
    !> @details Performs 1D spline interpolation along the z-direction for 3D data.
    !> This is an internal subroutine used by dbs3in to construct tensor-product splines.
    !>
    !> @param[in] xyzvec Array containing z-direction data points
    !> @param[in] ldf Leading dimension of the data array
    !> @param[in] mdf Middle dimension of the data array
    !> @param[in] xyzdata 3D data array to be interpolated
    !> @param[in] xyzknot Knot sequence for the z-direction
    !> @param[in] n Number of data points in the z-direction
    !> @param[in] k Order of the spline in the z-direction
    !> @param[in] m Number of data points in the x-direction
    !> @param[in] l Number of data points in the y-direction
    !> @param[out] work2 Work array for temporary storage
    !> @param[out] work3 Work array for banded matrix operations
    !> @param[out] bcoef Output B-spline coefficients
    !> @param[in] nx Number of x-direction data points
    !> @param[in] ny Number of y-direction data points
    !> @param[in] nz Number of z-direction data points
    subroutine spli3d(xyzvec, ldf, mdf, xyzdata, xyzknot, n, k, m, l, work2, work3, &
            & bcoef, nx, ny, nz)

        use numeric

        implicit none

        integer, intent(in) :: ldf, mdf, n, k, m, l
        integer, intent(in) :: nx, ny, nz
        real(kind=dbl), dimension(n), intent(in) :: xyzvec
        real(kind=dbl), dimension(n + k), intent(in) :: xyzknot
        real(kind=dbl), dimension(ldf, mdf, *), intent(in) :: xyzdata
        real(kind=dbl), dimension(nx, ny, nz), intent(out) :: bcoef
        real(kind=dbl), dimension(n), intent(out) :: work2
        real(kind=dbl), dimension((2*k - 1)*n), intent(out) :: work3

        integer :: np1, km1, kpkm2, left, lenq, i, ilp1mx, j, jj, iflag, in
        real(kind=dbl) :: xyzveci

        ! Initialize parameters for banded matrix construction
        np1 = n + 1
        km1 = k - 1
        kpkm2 = 2*km1
        left = k
        lenq = n*(k + km1)

        ! Initialize work array
        do i = 1, lenq
            work3(i) = 0._dbl
        end do

        ! Construct the banded matrix for the linear system
        do i = 1, n
            xyzveci = xyzvec(i)
            ilp1mx = min0(i + k, np1)
            left = max0(left, i)
            if (xyzveci .lt. xyzknot(left)) go to 998
30          if (xyzveci .lt. xyzknot(left + 1)) go to 40
            left = left + 1
            if (left .lt. ilp1mx) go to 30
            left = left - 1
            if (xyzveci .gt. xyzknot(left + 1)) go to 998
40          call bsplvb(xyzknot, n + k, k, 1, xyzveci, left, work2)
            jj = i - left + 1 + (left - k)*(k + km1)
            do j = 1, k
                jj = jj + kpkm2
                work3(jj) = work2(j)
            end do
        end do

        ! Factor the banded matrix
        call banfac(work3, k + km1, n, km1, km1, iflag)

        ! Check if factorization was successful
        if (iflag .ne. 1) then
            write (6, *) "subroutine dbs3in: error"
            write (6, *) "no solution of linear equation system !!!"
            stop
        end if

        ! Solve the linear system for each (x,y) pair
        do j = 1, l
            do i = 1, m
                do in = 1, n
                    work2(in) = xyzdata(i, j, in)
                end do

                call banslv(work3, k + km1, n, km1, km1, work2)

                do in = 1, n
                    bcoef(i, j, in) = work2(in)
                end do

            end do
        end do

        return

998     write (6, *) "subroutine db3in:"
        write (6, *) "i with knot(i) <= x/y/z < knot(i+1) required."
        write (6, *) "knot(1)   = ", xyzknot(1)
        write (6, *) "knot(n+k) = ", xyzknot(n + k)
        write (6, *) "    x/y/z = ", xyzveci

        stop

    end subroutine spli3d

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Evaluate a three-dimensional tensor-product spline at a point
    !> @details Evaluates a 3D tensor-product B-spline at a specified point (x,y,z).
    !> The evaluation uses the tensor-product structure to efficiently compute
    !> the spline value by evaluating 1D B-splines in each direction.
    !>
    !> @param[in] x x-coordinate of the evaluation point
    !> @param[in] y y-coordinate of the evaluation point
    !> @param[in] z z-coordinate of the evaluation point
    !> @param[in] kx Order of the spline in the x-direction
    !> @param[in] ky Order of the spline in the y-direction
    !> @param[in] kz Order of the spline in the z-direction
    !> @param[in] xknot Array of length nx+kx containing x-direction knot sequence (non-decreasing)
    !> @param[in] yknot Array of length ny+ky containing y-direction knot sequence (non-decreasing)
    !> @param[in] zknot Array of length nz+kz containing z-direction knot sequence (non-decreasing)
    !> @param[in] nx Number of B-spline coefficients in the x-direction
    !> @param[in] ny Number of B-spline coefficients in the y-direction
    !> @param[in] nz Number of B-spline coefficients in the z-direction
    !> @param[in] bcoef Array of length nx*ny*nz containing tensor-product B-spline coefficients
    !> @return Value of the spline at (x,y,z)
    !>
    !> @note All knot sequences must be non-decreasing
    !> @warning (x,y,z) must be within the knot ranges
    !> @see dbs2vl for 2D spline evaluation details
    function dbs3vl(x, y, z, kx, ky, kz, xknot, yknot, zknot, nx, ny, nz, bcoef)

        use numeric

        implicit none

        integer, intent(in) :: nx, ny, nz, kx, ky, kz
        real(kind=dbl), intent(in) :: x, y, z
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(ny + ky), intent(in) :: yknot
        real(kind=dbl), dimension(nz + kz), intent(in) :: zknot
        real(kind=dbl), dimension(nx, ny, nz), intent(in) :: bcoef
        real(kind=dbl) :: dbs3vl

        integer :: iz, nintz
        real(kind=dbl), dimension(kz) :: work

        ! Check z-direction knot sequence and find interval
        nintz = 0

        do iz = 1, nz + kz - 1
            if (zknot(iz) .gt. zknot(iz + 1)) then
                write (6, *) "subroutine dbs3vl:"
                write (6, *) "zknot(iz) <= zknot(iz+1) required."
                write (6, *) iz, zknot(iz), zknot(iz + 1)
                stop
            end if
            if ((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
        end do

        if (nintz .eq. 0) then
            write (6, *) "subroutine dbs3vl:"
            write (6, *) "iz with zknot(iz) <= z < zknot(iz+1) required."
            write (6, *) "zknot(iz)   = ", zknot(iz)
            write (6, *) "  z         = ", z
            write (6, *) "zknot(iz+1) = ", zknot(iz + 1)
            stop
        end if

        ! Evaluate 2D splines in (x,y) for each z basis function
        do iz = 1, kz
            work(iz) = dbs2vl(x, y, kx, ky, xknot, yknot, nx, ny, bcoef(1, 1, nintz - kz + iz))
        end do

        ! Evaluate the final result using z-direction spline
        dbs3vl = dbsval(z, kz, zknot(nintz - kz + 1), kz, work)

    end function dbs3vl

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dbs3dr(iderx, idery, iderz, x, y, z, kx, ky, kz, xknot, yknot, zknot, &
            & nx, ny, nz, bcoef)

        !
        !  Evaluates the derivative of a three-dimensional tensor-product spline,
        !  given its tensor-product B-spline representation.
        !
        !   iderx  - order of the x-derivative.  (input)
        !   idery  - order of the y-derivative.  (input)
        !   iderz  - order of the z-derivative.  (input)
        !   x      - x-coordinate of the point at which the spline is to be
        !            evaluated.  (input)
        !   y      - y-coordinate of the point at which the spline is to be
        !            evaluated.  (input)
        !   z      - z-coordinate of the point at which the spline is to be
        !            evaluated.  (input)
        !   kx     - order of the spline in the x-direction.  (input)
        !   ky     - order of the spline in the y-direction.  (input)
        !   kz     - order of the spline in the z-direction.  (input)
        !   xknot  - array of length nx+kx containing the knot
        !            sequence in the x-direction.  (input)
        !            xknot must be nondecreasing.
        !   yknot  - array of length ny+ky containing the knot
        !            sequence in the y-direction.  (input)
        !            yknot must be nondecreasing.
        !   zknot  - array of length nz+kz containing the knot
        !            sequence in the z-direction.  (input)
        !            zknot must be nondecreasing.
        !   nx     - number of B-spline coefficients in the x-direction.
        !            (input)
        !   ny     - number of B-spline coefficients in the y-direction.
        !            (input)
        !   nz     - number of B-spline coefficients in the z-direction.
        !            (input)
        !   bcoef  - array of length nx*ny*nz containing the
        !            tensor-product B-spline coefficients.  (input)
        !            bscoef is treated internally as a matrix of size nx
        !            by ny by nz.
        !   dbs3dr - value of the (iderx,idery,iderz) derivative of the
        !            spline at (x,y,z).  (output)
        !

        use numeric

        implicit none

        integer, intent(in) :: iderx, idery, iderz
        integer, intent(in) :: nx, ny, nz, kx, ky, kz
        real(kind=dbl), intent(in) :: x, y, z
        real(kind=dbl), dimension(nx + kx), intent(in) :: xknot
        real(kind=dbl), dimension(ny + ky), intent(in) :: yknot
        real(kind=dbl), dimension(nz + kz), intent(in) :: zknot
        real(kind=dbl), dimension(nx, ny, nz), intent(in) :: bcoef
        real(kind=dbl) :: dbs3dr

        integer :: iz, nintz
        real(kind=dbl), dimension(kz) :: work

        !
        !     check if knot(i) <= knot(i+1) and calculation of i so that
        !     knot(i) <= x < knot(i+1)
        !

        nintz = 0

        do iz = 1, nz + kz - 1
            if (zknot(iz) .gt. zknot(iz + 1)) then
                write (6, *) "subroutine dbs3vl:"
                write (6, *) "zknot(iz) <= zknot(iz+1) required."
                write (6, *) iz, zknot(iz), zknot(iz + 1)
                stop
            end if
            if ((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
        end do

        if (nintz .eq. 0) then
            write (6, *) "subroutine dbs3dr:"
            write (6, *) "iz with zknot(iz) <= z < zknot(iz+1) required."
            write (6, *) "zknot(iz)   = ", zknot(iz)
            write (6, *) "  z         = ", z
            write (6, *) "zknot(iz+1) = ", zknot(iz + 1)
            stop
        end if

        do iz = 1, kz
            work(iz) = dbs2dr(iderx, idery, x, y, kx, ky, xknot, yknot, nx, ny, &
                    &        bcoef(1, 1, nintz - kz + iz))
        end do

        dbs3dr = dbsder(iderz, z, kz, zknot(nintz - kz + 1), kz, work)

    end function dbs3dr

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Compute B-spline basis function values
    !> @details Computes the values of B-spline basis functions at a given point.
    !> This is a core subroutine used by most B-spline evaluation routines.
    !> The algorithm uses the recursive de Boor formula to compute basis functions.
    !>
    !> @param[in] t Array containing the knot sequence
    !> @param[in] n Length of the knot sequence
    !> @param[in] jhigh Highest order of basis functions to compute
    !> @param[in] index Index for controlling computation (1 for initialization)
    !> @param[in] x Point at which to evaluate basis functions
    !> @param[in] left Index of the knot interval containing x
    !> @param[out] biatx Array containing the computed basis function values
    !>
    !> @note This subroutine is called internally by other B-spline routines
    !> @see de Boor, C. (1978). A practical guide to Splines. Springer-Verlag.
    subroutine bsplvb(t, n, jhigh, index, x, left, biatx)

        use numeric

        implicit none

        integer, intent(in) :: n, jhigh, index, left

        real(kind=dbl), intent(in) :: x
        real(kind=dbl), dimension(n), intent(in) :: t
        real(kind=dbl), dimension(jhigh), intent(out) :: biatx

        integer :: j = 1
        integer :: i, jp1
        real(kind=dbl) :: saved, term
        real(kind=dbl), dimension(jhigh) :: dl, dr

        ! Initialize basis functions if index = 1
        if (index .eq. 1) then
            j = 1
            biatx(1) = 1.0_dbl
            if (j .ge. jhigh) return
        end if

20      jp1 = j + 1

        ! Compute distances for de Boor algorithm
        dr(j) = t(left + j) - x
        dl(j) = x - t(left + 1 - j)
        saved = 0._dbl

        ! Apply de Boor recursion formula
        do i = 1, j
            term = biatx(i)/(dr(i) + dl(jp1 - i))
            biatx(i) = saved + dr(i)*term
            saved = dl(jp1 - i)*term
        end do

        biatx(jp1) = saved
        j = jp1

        ! Continue until all required basis functions are computed
        if (j .lt. jhigh) go to 20

    end subroutine bsplvb

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Factor a banded matrix for linear system solution
    !> @details Performs LU factorization of a banded matrix for efficient
    !> solution of linear systems. This is used in B-spline interpolation
    !> to solve the system of equations for computing spline coefficients.
    !>
    !> @param[in,out] w Banded matrix to be factored (overwritten with factorization)
    !> @param[in] nroww Number of rows in the work array w
    !> @param[in] nrow Number of rows in the actual matrix
    !> @param[in] nbandl Number of subdiagonals
    !> @param[in] nbandu Number of superdiagonals
    !> @param[out] iflag Status flag (1 = success, 2 = failure)
    !>
    !> @note The matrix is stored in banded format for efficiency
    !> @warning The matrix must be non-singular for successful factorization
    subroutine banfac(w, nroww, nrow, nbandl, nbandu, iflag)

        use numeric

        implicit none

        integer, intent(in) :: nroww, nrow
        integer, intent(in) :: nbandl, nbandu
        integer, intent(out) :: iflag
        real(kind=dbl), dimension(nroww, nrow), intent(inout) :: w

        real(kind=dbl) :: pivot, factor
        integer :: middle, nrowm1, jmax, kmax, ipk, midmk, i, j, k

        iflag = 1
        middle = nbandu + 1
        nrowm1 = nrow - 1

        ! Handle special cases
        if (nrowm1 .lt. 0) goto 999
        if (nrowm1 .eq. 0) goto 900
        if (nrowm1 .gt. 0) goto 10

10      if (nbandl .gt. 0) go to 30

        do i = 1, nrowm1
            if (w(middle, i) .eq. 0._dbl) go to 999
        end do

        go to 900

30      if (nbandu .gt. 0) go to 60

        do i = 1, nrowm1
            pivot = w(middle, i)
            if (pivot .eq. 0._dbl) go to 999
            jmax = min0(nbandl, nrow - i)
            do j = 1, jmax
                w(middle + j, i) = w(middle + j, i)/pivot
            end do
        end do

        return

60      do i = 1, nrowm1
            pivot = w(middle, i)
            if (pivot .eq. 0._dbl) go to 999
            jmax = min0(nbandl, nrow - i)
            do j = 1, jmax
                w(middle + j, i) = w(middle + j, i)/pivot
            end do

            kmax = min0(nbandu, nrow - i)

            do k = 1, kmax
                ipk = i + k
                midmk = middle - k
                factor = w(midmk, ipk)
                do j = 1, jmax
                    w(midmk + j, ipk) = w(midmk + j, ipk) - w(middle + j, i)                  &
                                       & *factor
                end do
            end do
        end do

900     if (w(middle, nrow) .ne. 0._dbl) return
999     iflag = 2

    end subroutine banfac

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Solve a banded linear system using LU factorization
    !> @details Solves a banded linear system Ax = b using the LU factorization
    !> computed by banfac. This is used in B-spline interpolation to compute
    !> spline coefficients from the interpolation conditions.
    !>
    !> @param[in] w Banded matrix in LU factored form (from banfac)
    !> @param[in] nroww Number of rows in the work array w
    !> @param[in] nrow Number of rows in the actual matrix
    !> @param[in] nbandl Number of subdiagonals
    !> @param[in] nbandu Number of superdiagonals
    !> @param[in,out] b Right-hand side vector (overwritten with solution)
    !>
    !> @note The matrix must be factored by banfac before calling this routine
    !> @warning The matrix must be non-singular
    subroutine banslv(w, nroww, nrow, nbandl, nbandu, b)

        use numeric

        implicit none

        integer, intent(in) :: nroww, nrow
        integer, intent(in) :: nbandl, nbandu
        real(kind=dbl), dimension(nroww, nrow), intent(in) :: w
        real(kind=dbl), dimension(nrow), intent(inout) :: b

        integer :: middle, nrowm1, jmax, i, j

        middle = nbandu + 1
        if (nrow .eq. 1) goto 99
        nrowm1 = nrow - 1
        if (nbandl .eq. 0) goto 30

        do i = 1, nrowm1
            jmax = min0(nbandl, nrow - i)
            do j = 1, jmax
                b(i + j) = b(i + j) - b(i)*w(middle + j, i)
            end do
        end do

30      if (nbandu .gt. 0) goto 50

        do i = 1, nrow
            b(i) = b(i)/w(1, i)
        end do

        return

50      do i = nrow, 2, -1
            b(i) = b(i)/w(middle, i)
            jmax = min0(nbandu, i - 1)
            do j = 1, jmax
                b(i - j) = b(i - j) - b(i)*w(middle - j, i)
            end do
        end do

99      b(1) = b(1)/w(middle, 1)

    end subroutine banslv

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> @brief Efficiently find knot interval for a given point
    !> @details Uses the hunt algorithm to efficiently locate the knot interval
    !> containing a given point. This is optimized for sequential evaluation
    !> where points are ordered, as it can reuse previous interval information.
    !>
    !> The hunt algorithm is particularly efficient for grid-based evaluation
    !> where consecutive points are likely to be in the same or adjacent intervals.
    !>
    !> @param[in] xx Array containing the knot sequence
    !> @param[in] n Length of the knot sequence
    !> @param[in] kord Order of the B-spline
    !> @param[in] x Point for which to find the interval
    !> @param[in,out] jlo Index of the interval containing x (updated on output)
    !>
    !> @note This routine is optimized for sequential point evaluation
    !> @warning The knot sequence must be non-decreasing
    subroutine huntn(xx, n, kord, x, jlo)

        use numeric

        implicit none

        integer, intent(in) :: n, kord
        real(kind=dbl), intent(in) :: x
        real(kind=dbl), dimension(n), intent(in) :: xx

        integer, intent(inout) :: jlo

        integer :: max, null, jhi, jm, inc

        ! Set bounds for valid intervals
        max = n - kord
        null = kord

        ! Initialize search if jlo is out of bounds
        if (jlo .le. null .or. jlo .gt. max) then
            jlo = null
            jhi = max + 1
            goto 30
        end if

        inc = 1

        ! Hunt forward if x is in or after current interval
        if (x .ge. xx(jlo)) then
10          jhi = jlo + inc
            if (jhi .gt. max) then
                jhi = max + 1
            else if (x .ge. xx(jhi)) then
                jlo = jhi
                inc = inc + inc
                goto 10
            end if
        else
            ! Hunt backward if x is before current interval
            jhi = jlo
20          jlo = jhi - inc
            if (jlo .le. null) then
                jlo = null
            else if (x .lt. xx(jlo)) then
                jhi = jlo
                inc = inc + inc
                goto 20
            end if
        end if

30      if (jhi - jlo .eq. 1) return

        jm = (jhi + jlo)/2
        if (x .gt. xx(jm)) then
            jlo = jm
        else
            jhi = jm
        end if

        goto 30

    end subroutine huntn

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module bspline
