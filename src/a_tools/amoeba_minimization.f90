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

module amoeba_minimization

    implicit none

    integer :: param_number, simplex_dim, points_number, itmax, iter, param_number_read, points_number_read
    real(8), dimension(:), allocatable :: param0
    real(8), dimension(:, :), allocatable :: points
    logical, dimension(:), allocatable :: sel
    real(8), dimension(:, :), allocatable :: p
    real(8), dimension(:), allocatable :: y, x, xscratch, vecscratch, vecsigma
    real(8) :: delta, tolerance
    real(8), external :: pseudofun

! INPUT ###################################################

! PARAMETERS
! param_number: number of parameters
! param0: initial set of parameters
! sel: constraint on the parameter freedom
! true --> free to move
! false --> fixed

! DATA SET
! points_number: number of data points to fit
! points: data to fit
! points(1,i): x_i
! points(2,i): y_i

! MAIN FUNCTION
! pseudofun defined elsewhere

! CONTROL PARAMETERS
! sigma: estimated error in the fit
! (vecsigma = vecsigma(i) --> it is renormalized based on the density of the x_i points such that
! \sum_i vecsigma(i) / points_number = sigma )
! delta: initial offset for the simplex
! tolerance: tolerance for amoeba convergence
! itmax: max number of amoeba iterations
! alpha, beta, gamma: internal amoeba parameters

! INTERNAL #################################################

! INTERNAL ARRAYS AND VARIABLES
! simplex_dim: number of parameters modified
! p(i,j): i-th simplex vertex vec_i(j)
! y(i): value of the chi squared at the i-th vertex of the simplex
! x(i): parameters modified during minimization

    private
    public :: param0, sel, points, simplex_dim, &
              delta, tolerance, itmax, p, y, x, vecscratch, chi, &
              alloc_datavec, dealloc_datavec, simplex_constructor, amoeba, &
              dealloc_simplex, update_param, iter, vecsigma

contains

    subroutine alloc_datavec(param_number_read, points_number_read)
        integer, intent(in) :: param_number_read, points_number_read
        param_number = param_number_read
        points_number = points_number_read
        allocate (param0(param_number), sel(param_number))
        allocate (points(2, points_number))
        allocate (vecsigma(points_number))
    end subroutine alloc_datavec

    subroutine dealloc_datavec
        if (allocated(param0)) deallocate (param0)
        if (allocated(sel)) deallocate (sel)
        if (allocated(points)) deallocate (points)
        if (allocated(vecsigma)) deallocate (vecsigma)
    end subroutine dealloc_datavec

    subroutine dealloc_simplex
        if (allocated(y)) deallocate (y)
        if (allocated(p)) deallocate (p)
        if (allocated(x)) deallocate (x)
        if (allocated(xscratch)) deallocate (xscratch)
        if (allocated(vecscratch)) deallocate (vecscratch)
    end subroutine dealloc_simplex

    subroutine update_param(p, sel, param0)

        real(8), intent(in) :: p(simplex_dim + 1, *)
        logical, intent(in) :: sel(*)
        real(8), intent(out) :: param0(*)
        integer :: i, j

        j = 0
        do i = 1, param_number
        if (sel(i)) then
            j = j + 1
            param0(i) = p(1, j)
        end if
        end do

    end subroutine update_param

    subroutine simplex_constructor(delta, vecsigma, weights)

! to be called at the beginning to initialize the simplex correctly!!!!!
        real(8) :: vecsigma(*), weights(*), delta, sigma, norm
        integer :: i, j

! define vecsigma
! pin the value at the origin
!    vecsigma(1)=1.d0/abs(points(1,1)-points(1,2))/dble(points_number)/1000.d0
!    do i=2,points_number-1
!       vecsigma(i)=1.d0/abs(points(1,i)-points(1,i+1))/dble(points_number)
!    enddo
!    vecsigma(points_number)=vecsigma(points_number-1)

!    norm=1.d0
!    do i=1,points_number
!       norm=norm+vecsigma(i)
!    enddo
!    norm=norm*abs(points(1,points_number)-points(1,1))/dble(points_number)**2
!    do i=1,points_number
!       vecsigma(i)=vecsigma(i)*sigma/norm
!    enddo

!use same matrics as in the diagonalization algorithm for consistency
!    dx(i)=1/vecsigma(i)**2
!    vecsigma(i)=1.d0/sqrt(dx(i))

! old integration
!    vecsigma(1)=(points(1,2)-points(1,1))/2.d0+points(1,1)
!    vecsigma(points_number)=(points(1,points_number)-points(1,points_number-1))/2.d0
!    do i=2,points_number-1
!       vecsigma(i)=0.5d0*(points(1,i+1)-points(1,i-1))
!    enddo

        do i = 1, points_number
            vecsigma(i) = points(1, i)
        end do
        do i = 1, points_number
! normalize vecsigma --> dx / Ltot
            vecsigma(i) = vecsigma(i)*weights(i)/abs(points(1, points_number) - points(1, 1))
            vecsigma(i) = 1.d0/sqrt(vecsigma(i))
        end do

! define the simplex
        j = 0
        do i = 1, param_number
            if (sel(i)) j = j + 1
        end do
        simplex_dim = j

        allocate (y(simplex_dim + 1), p(simplex_dim + 1, simplex_dim))
        allocate (x(simplex_dim))
        allocate (xscratch(param_number))
        allocate (vecscratch(param_number/3))

        j = 0
        do i = 1, param_number
        if (sel(i)) then
            j = j + 1
            p(1, j) = param0(i)
            x(j) = p(1, j)
        end if
        end do
        y(1) = chi(x, vecsigma)

        do i = 1, simplex_dim
        do j = 1, simplex_dim
        if (j .ne. i) then
            p(i + 1, j) = p(1, j)
        else
            p(i + 1, j) = p(1, j) + delta
        end if
        x(j) = p(i + 1, j)
        end do
        y(i + 1) = chi(x, vecsigma)
        end do

        write (6, *) 'initial chi square', y(1)

    end subroutine simplex_constructor

    function chi(x, vecsigma)

! chi: cost function for minimization
! chi "measures" the distance of the fitted form with respect to the datapoints
! vecsigma: control parameter to tune the distance
        real(8) :: chi
        real(8) :: x(*)
        real(8) :: vecsigma(*)
        integer :: i, j

        do i = 1, param_number
            xscratch(i) = param0(i)
        end do

        j = 0
        do i = 1, param_number
        if (sel(i)) then
            j = j + 1
            xscratch(i) = x(j)
        end if
        end do

        chi = 0.d0
        do i = 1, points_number
            chi = chi + ((pseudofun(param_number/3, points(1, i), xscratch, vecscratch) - points(2, i))/vecsigma(i))**2
        end do
!    chi=chi/(points_number-simplex_dim)

    end function chi

    subroutine amoeba(p, y, tolerance, itmax)
! kernel of the minimization module

        real(8) :: p(simplex_dim + 1, *), y(*)
        integer :: itmax
        real(8) :: tolerance, ftol

        real(8), parameter :: alpha = 1.0
        real(8), parameter :: beta = 0.5
        real(8), parameter :: gamma = 2.0

        integer :: mpts, ndim, mp, np
        integer :: ilo, ihi, inhi, i, j

        real(8), dimension(:), allocatable :: pr, prr, pbar
        real(8) :: ypr, yprr, rtol

        ndim = simplex_dim
        np = ndim
        mp = np + 1
        mpts = np + 1 !mpts logical dimension of the simplex

        allocate (pr(ndim), prr(ndim), pbar(ndim))

!    real(8) :: p(mp,np),y(mp)

        ftol = tolerance
        iter = 1

        do while (iter .le. itmax)
            ilo = 1
            if (y(1) .gt. y(2)) then
                ihi = 1
                inhi = 2
            else
                ihi = 2
                inhi = 1
            end if
            do i = 1, mpts
                if (y(i) .lt. y(ilo)) ilo = i
                if (y(i) .gt. y(ihi)) then
                    inhi = ihi
                    ihi = i
                elseif (y(i) .gt. y(inhi)) then
                    if (i .ne. ihi) inhi = i
                end if
            end do

            rtol = 2.d0*abs(y(ihi) - y(ilo))/(abs(y(ihi)) + abs(y(ilo)))
            if (rtol .le. ftol) go to 10
            do j = 1, ndim
                pbar(j) = 0.d0
            end do
            do i = 1, mpts
            if (i .ne. ihi) then
            do j = 1, ndim
                pbar(j) = pbar(j) + p(i, j)
            end do
            end if
            end do
            do j = 1, ndim
                pbar(j) = pbar(j)/dble(ndim)
                pr(j) = (1 + alpha)*pbar(j) - alpha*p(ihi, j)
            end do
            ypr = chi(pr, vecsigma)
            if (ypr .le. y(ilo)) then
            do j = 1, ndim
                prr(j) = gamma*pr(j) + (1.d0 - gamma)*pbar(j)
            end do
            yprr = chi(prr, vecsigma)
            if (yprr .lt. y(ilo)) then
            do j = 1, ndim
                p(ihi, j) = prr(j)
            end do
            y(ihi) = yprr
            else
            do j = 1, ndim
                p(ihi, j) = pr(j)
            end do
            y(ihi) = ypr
            end if
            elseif (ypr .ge. y(inhi)) then
            if (ypr .lt. y(ihi)) then
            do j = 1, ndim
                p(ihi, j) = pr(j)
            end do
            y(ihi) = ypr
            end if
            do j = 1, ndim
                prr(j) = beta*p(ihi, j) + (1.d0 - beta)*pbar(j)
            end do
            yprr = chi(prr, vecsigma)
            if (yprr .lt. y(ihi)) then
            do j = 1, ndim
                p(ihi, j) = prr(j)
            end do
            y(ihi) = yprr
            else
            do i = 1, mpts
                if (i .ne. ilo) then
                    do j = 1, ndim
                        pr(j) = 0.5*(p(i, j) + p(ilo, j))
                        p(i, j) = pr(j)
                    end do
                    y(i) = chi(pr, vecsigma)
                end if
            end do
            end if
            else
            do j = 1, ndim
                p(ihi, j) = pr(j)
            end do
            y(ihi) = ypr
            end if
!         write(*,*) iter
!         do j=1,ndim+1
!         write(*,*) y(j),(p(j,i),i=1,ndim)
!         enddo
            iter = iter + 1
        end do

        write (*, *) 'WARNING: amoeba exceeding maximum iterations'

10      deallocate (pr, prr, pbar)

        return

    end subroutine amoeba

end module amoeba_minimization

