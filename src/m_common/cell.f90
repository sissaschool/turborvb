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

module Cell

    use constants
    use symmetries

    implicit none

    ! cell(6) contains (a,b,c) and (alpha,beta,gamma) (in radians)
    ! notice that cell(2)=b/a and cell(3)=c/a
    ! by convention, a is placed parallel to x, and b lies on the xy plane
    real(8), dimension(6) :: celldm(6), celldm2(3)

    ! constant for pressure
    real(8) :: costpr, metric_min

    ! r2s is the matrix for cartesian -> crystal conversion
    ! s2r is the matrix for crystal -> cartesian conversion
    ! metric is the metric matrix of the crystal
    ! at are the direct lattice vectors
    ! recip is the reciprocal lattice
    real(8), dimension(3, 3) :: r2s, s2r, metric, at, recip, bg, car2cry

    ! Vector used to rescale derivatives
    real(8), dimension(12) :: cellscale
    real(8), dimension(3) :: cellpi
    real(8), allocatable :: cellscalep(:, :), x_neigh(:, :), rphasep(:, :)&
   &, s2rp(:, :, :), dist_shift(:), distreg_shift(:), disto_shift(:), distrego_shift(:)
    ! celldm volume
    real(8) :: omega, unit_volume

    ! L min shorter side of the box
    real(8) :: LMin

    ! flag to check if derivatives respect to cell are evaluated
    logical :: cellderiv, givens2r, yes_tilted, chosen_map

    ! t_rev contains a flag which check if time reversal can be applied to each point group symmetry
    integer :: case_map, t_rev(48), neigh

    ! various vectors which store the phase of the wave function
    ! phase/phase_down = phase of the wave function in crystal coordinates for up/down electrons
    ! phase2pi/phase2pi_down = 2*PI*phase/phase_down
    ! cell_phase = (2*PI*phase)/L
    ! cell_phase2 = sum(cell_phase(:)**2)
    ! rphase = (2*PI*phase)/L used only in calculation with old PBC (not Crystal) basis set.
    !          In the present version, only the Jastrow uses this old basis set.
    real(8) :: phase(3), phase2pi(3), phase_down(3), phase2pi_down(3)
    real(8) :: rphase(3), cell_phase(3), cell_phase2

    real(8), allocatable :: cosphase(:, :), sinphase(:, :, :)
    real(8), allocatable :: cosphaseb(:, :), sinphaseb(:, :, :)

    ! flag to check if we are doing a Gamma Point calculation
    logical :: gamma_point, yes2d, yes1d
    double precision :: amap, bmap
    double precision, parameter :: x_c = 0.25d0

    !  real(8) :: map, dmap, ddmap

contains
    !=====================================================================
    ! Given celldm(1:6), calculate the above matrices
    !=====================================================================
    subroutine InitCell(nion, nel, yes_complex)
        integer i, nion, nel, info, ipiv(3)
        double precision :: a, b, c, alpha, beta, gamma, vec_mod
        double precision matscra(3, 3), work(9), eigscra(3)
        logical :: yes_complex

        if (.not. givens2r) then
            !   This part holds only for orthorombic supercells.
            yes2d = .false.
            yes1d = .false.

            a = celldm(1)
            b = celldm(2)*a
            c = celldm(3)*a
            celldm(4:6) = 90.d0*PI/180.d0
            alpha = celldm(4); beta = celldm(5); gamma = celldm(6)
            if (c .eq. 0.d0 .and. b .eq. 0.d0) then
                omega = a
                yes1d = .true.
            elseif (c .eq. 0.d0) then
                omega = a*b*dsin(gamma)
                yes2d = .true.
            else
                omega = (a*b*c)*dsqrt(&
                        & 1.d0 - dcos(alpha)**2.d0 - dcos(beta)**2.d0 - dcos(gamma)**2.d0 + &
                                & 2.d0*dcos(alpha)*dcos(beta)*dcos(gamma))
            end if
            celldm2(1) = celldm(1)**2
            celldm2(2) = celldm(2)**2
            celldm2(3) = celldm(3)**2
            t_rev(:) = 0

            ! direct space cell
            s2r(1, 1) = a
            s2r(2, 1) = 0.d0
            s2r(3, 1) = 0.d0
            s2r(1, 2) = b*dcos(gamma)
            s2r(2, 2) = b*dsin(gamma)
            s2r(3, 2) = 0.d0
            s2r(1, 3) = c*dcos(beta)
            s2r(2, 3) = c*(dcos(alpha) - dcos(beta)*dcos(gamma))       &
                    & /dsin(gamma)
            s2r(3, 3) = omega/(a*b*dsin(gamma))
            ! direct lattice versors
            if (yes2d .or. yes1d) s2r(3, 3) = 1.d0
            if (yes1d) s2r(2, 2) = 1.d0

        else
            !    compute volume given s2r
            matscra = s2r
            call dgetrf(3, 3, matscra, 3, ipiv, info)
            omega = matscra(1, 1)
            do i = 2, 3
                omega = omega*matscra(i, i)
            end do
            omega = abs(omega)
        end if

        !   calculation cellscale given s2r
        do i = 1, 3
            cellscale(i) = sqrt(s2r(1, i)**2 + s2r(2, i)**2 + s2r(3, i)**2)
        end do
        if (yes2d .or. yes1d) cellscale(3) = 0.d0
        if (yes1d) cellscale(2) = 0.d0

        !   calculation  at  given s2r and cellscale
        do i = 1, 3
            if (cellscale(i) .ne. 0.d0) then
                at(:, i) = s2r(:, i)/cellscale(i)
            else
                at(:, i) = s2r(:, i)
            end if
        end do
        unit_volume = omega/(cellscale(1)*cellscale(2)*cellscale(3))

        ! reciprocal space cell
        r2s(1, :) = cross_product(s2r(1, 2), s2r(1, 3))/omega
        r2s(2, :) = cross_product(s2r(1, 3), s2r(1, 1))/omega
        r2s(3, :) = cross_product(s2r(1, 1), s2r(1, 2))/omega
        if (yes2d .or. yes1d) r2s(3, :) = 0.d0
        if (yes1d) r2s(2, :) = 0.d0
        recip(:, :) = r2s(:, :)*(2.d0*PI)
        ! reciprocal lattice versors
        bg = 0.d0
        do i = 1, 3
            vec_mod = sqrt(recip(1, i)**2 + recip(2, i)**2 + recip(3, i)**2)
            if (vec_mod .ne. 0.d0) then
                bg(1, i) = recip(1, i)/vec_mod
                bg(2, i) = recip(2, i)/vec_mod
                bg(3, i) = recip(3, i)/vec_mod
            end if
        end do

        ! metric matrix
        metric = matmul(transpose(at), at)
        !   compute the minimum eigenvalue of the metric < 1
        matscra = metric
        call dsyev('N', 'L', 3, matscra, 3, eigscra, work, 9, info)
        metric_min = dsqrt(eigscra(1))

        car2cry = at
        call dgetrf(3, 3, car2cry, 3, ipiv, info)
        if (info .ne. 0) then
            write (6, *) ' ERROR in initialization cell (dgetrf) !!! '
        else
            call dgetri(3, car2cry, 3, ipiv, matscra, 9, info)
            if (info .ne. 0) write (6, *) ' ERROR in initialization cell (dgetri) !!! '
        end if

        !   if(.not.givens2r)  then
        !   cellscale(1)=celldm(1)
        !   cellscale(2:3)=celldm(2:3)*celldm(1)
        !   endif
        !   if(.not.allocated(cellpi)) allocate(cellpi(3))
        !   s2r saved in cellscale
        call dcopy(9, s2r, 1, cellscale(4), 1)

        cellpi(1) = cellscale(1)/Pi
        cellpi(2) = cellscale(2)/Pi
        cellpi(3) = cellscale(3)/Pi

        if (cellscale(3) .ne. 0.d0 .and. cellscale(2) .ne. 0.d0) then
            lmin = min(cellscale(1), cellscale(2), cellscale(3))
        elseif (cellscale(2) .ne. 0.d0) then
            lmin = min(cellscale(1), cellscale(2))
        else
            lmin = cellscale(1)
        end if

        ! rphase used for real boundary conditions only
        ! rphase superceded by cell_phase in the case of complex wave function (yes_complex=.true.)
        rphase(:) = 0.d0
        if (.not. yes_complex) rphase(1:3) = phase(1:3)/cellscale(1:3)*TWO_PI
        cell_phase(1:3) = phase(1:3)/cellscale(1:3)*TWO_PI
        cell_phase2 = sum(cell_phase(:)**2)
        phase2pi(:) = phase(:)*TWO_PI
        phase2pi_down(:) = phase_down(:)*TWO_PI

        costpr = 1.d0/3.d0/omega

        if (allocated(cellscalep)) deallocate (cellscalep, rphasep, s2rp)
        allocate (cellscalep(3, nion), rphasep(3, nion), s2rp(3, 3, nion))

        bmap = 2.d0*(4.d0*x_c**2 - x_c)/(6.d0*x_c - 1.d0)
        amap = 0.25d0*(bmap - x_c)**2*(1.d0 - 2.d0*x_c)**3/(1.d0 + 4.d0*bmap - 6.d0*x_c)

        if (allocated(x_neigh)) deallocate (x_neigh, dist_shift, distreg_shift, disto_shift, distrego_shift)
        allocate (x_neigh(neigh, 3), dist_shift(neigh), distreg_shift(neigh), disto_shift(neigh), distrego_shift(neigh))
        dist_shift = 0.d0
        distreg_shift = 0.d0
        disto_shift = 0.d0
        distrego_shift = 0.d0
        x_neigh(1, :) = 0.d0
        if (neigh .ge. 7) then
            x_neigh(2, :) = s2r(:, 1)
            x_neigh(3, :) = -s2r(:, 1)
            x_neigh(4, :) = s2r(:, 2)
            x_neigh(5, :) = -s2r(:, 2)
            x_neigh(6, :) = s2r(:, 3)
            x_neigh(7, :) = -s2r(:, 3)
        end if
        if (neigh .ge. 19) then
            x_neigh(8, :) = s2r(:, 1) + s2r(:, 2)
            x_neigh(9, :) = -(s2r(:, 1) + s2r(:, 2))
            x_neigh(10, :) = (s2r(:, 1) + s2r(:, 3))
            x_neigh(11, :) = -(s2r(:, 1) + s2r(:, 3))
            x_neigh(12, :) = (s2r(:, 2) + s2r(:, 3))
            x_neigh(13, :) = -(s2r(:, 2) + s2r(:, 3))
            x_neigh(14, :) = s2r(:, 1) - s2r(:, 2)
            x_neigh(15, :) = -(s2r(:, 1) - s2r(:, 2))
            x_neigh(16, :) = (s2r(:, 1) - s2r(:, 3))
            x_neigh(17, :) = -(s2r(:, 1) - s2r(:, 3))
            x_neigh(18, :) = (s2r(:, 2) - s2r(:, 3))
            x_neigh(19, :) = -(s2r(:, 2) - s2r(:, 3))
        end if
        if (neigh .ge. 27) then
            x_neigh(20, :) = s2r(:, 1) + s2r(:, 2) + s2r(:, 3)
            x_neigh(21, :) = -(s2r(:, 1) + s2r(:, 2) + s2r(:, 3))
            x_neigh(22, :) = s2r(:, 1) - s2r(:, 2) + s2r(:, 3)
            x_neigh(23, :) = -(s2r(:, 1) - s2r(:, 2) + s2r(:, 3))
            x_neigh(24, :) = s2r(:, 1) + s2r(:, 2) - s2r(:, 3)
            x_neigh(25, :) = -(s2r(:, 1) + s2r(:, 2) - s2r(:, 3))
            x_neigh(26, :) = -s2r(:, 1) + s2r(:, 2) + s2r(:, 3)
            x_neigh(27, :) = -(-s2r(:, 1) + s2r(:, 2) + s2r(:, 3))
        end if
        if (neigh .ge. 33) then
            x_neigh(28, :) = 2.d0*s2r(:, 1)
            x_neigh(29, :) = -2.d0*s2r(:, 1)
            x_neigh(30, :) = 2.d0*s2r(:, 2)
            x_neigh(31, :) = -2.d0*s2r(:, 2)
            x_neigh(32, :) = 2.d0*s2r(:, 3)
            x_neigh(33, :) = -2.d0*s2r(:, 3)
        end if
    end subroutine InitCell

    ! Perform the cross product between two vectors
    !=====================================================================
    function cross_product(a, b)
        double precision, dimension(3), intent(in) :: a, b
        double precision, dimension(3) :: cross_product
        cross_product(1) = a(2)*b(3) - a(3)*b(2)
        cross_product(2) = a(3)*b(1) - a(1)*b(3)
        cross_product(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    !====================================================================
    ! Convert cartesian coordinates to crystal
    !====================================================================
    subroutine CartesianToCrystal(r, howmany)
        integer, intent(in) :: howmany
        double precision, dimension(3, howmany), intent(inout) :: r
        real*8 s(3)
        integer i
!$omp parallel do default(shared) private(i,s)
        do i = 1, howmany
            !     call dgemv('N',3,3,1.d0,car2cry,3,s,1,0.d0,r(1,i),1)
            s(1) = r(1, i)
            s(2) = r(2, i)
            s(3) = r(3, i)
            r(1, i) = car2cry(1, 1)*s(1) + car2cry(1, 2)*s(2) + car2cry(1, 3)*s(3)
            r(2, i) = car2cry(2, 1)*s(1) + car2cry(2, 2)*s(2) + car2cry(2, 3)*s(3)
            r(3, i) = car2cry(3, 1)*s(1) + car2cry(3, 2)*s(2) + car2cry(3, 3)*s(3)
        end do
!$omp end parallel do
    end subroutine CartesianToCrystal
    subroutine CartesianToCrystal_b(rbefore, rb, car2cryb, howmany)
        integer, intent(in) :: howmany
        double precision, dimension(3, howmany), intent(in) :: rbefore
        double precision, dimension(3, howmany), intent(inout) :: rb
        double precision, dimension(3, 3), intent(inout) :: car2cryb
        !   double precision s(3),sb(3)
        integer i
        do i = 1, howmany
            !     s(:)=rbefore(:,i)
            !     sb=0.d0
            !     call dgemv_b('N',3,3,1.d0,car2cry,3,car2cryb,3,s,1,sb,1,0.d0,rb(1,i),1)
            !     rb(:,i)=sb(:)
            car2cryb(:, 1) = car2cryb(:, 1) + rb(:, i)*rbefore(1, i)
            car2cryb(:, 2) = car2cryb(:, 2) + rb(:, i)*rbefore(2, i)
            car2cryb(:, 3) = car2cryb(:, 3) + rb(:, i)*rbefore(3, i)
            rb(:, i) = car2cry(1, :)*rb(1, i) + car2cry(2, :)*rb(2, i) + car2cry(3, :)*rb(3, i)
        end do
    end subroutine CartesianToCrystal_b

    ! adapted from QuantumESPRESSO
    subroutine cryst_to_cart(nvec, vec, trmat, iflag)
        !
        !     This routine transforms the atomic positions or the k-point
        !     components from crystallographic to cartesian coordinates
        !     ( iflag=1 ) and viceversa ( iflag=-1 ) for a set of vectors.
        !     Output cartesian coordinates are stored in the input ('vec') array
        !
        implicit none
        !
        integer, intent(in) :: nvec, iflag
        ! nvec:  number of vectors (atomic positions or k-points)
        !        to be transformed from crystal to cartesian and vice versa
        ! iflag: gives the direction of the transformation
        real(DP), intent(in) :: trmat(3, 3)
        ! trmat: transformation matrix
        ! if iflag=1:
        !    trmat = at,s2r ,    basis (or vectors) of the real-space lattice,       for atomic positions
        !          = bg,recip ,  basis (or vectors) of the reciprocal-space lattice, for k-points
        ! if iflag=-1: the opposite
        real(DP), intent(inout) :: vec(3, nvec)
        ! coordinates of the vector (atomic positions or k-points) to be
        ! transformed - overwritten on output
        !
        !    local variables
        !
        integer :: nv, kpol
        ! counter on vectors
        ! counter on polarizations
        real(DP) :: vau(3)
        ! workspace
        !
        !     Compute the cartesian coordinates of each vectors
        !     (atomic positions or k-points components)
        !
        do nv = 1, nvec
            if (iflag .eq. 1) then
                do kpol = 1, 3
                    vau(kpol) = trmat(kpol, 1)*vec(1, nv) + trmat(kpol, 2) &
                                *vec(2, nv) + trmat(kpol, 3)*vec(3, nv)
                end do
            else
                do kpol = 1, 3
                    vau(kpol) = trmat(1, kpol)*vec(1, nv) + trmat(2, kpol) &
                                *vec(2, nv) + trmat(3, kpol)*vec(3, nv)
                end do
            end if
            do kpol = 1, 3
                vec(kpol, nv) = vau(kpol)
            end do
        end do
        !
        return
    end subroutine cryst_to_cart

    subroutine ApplyPBC(s, howmany)
        implicit none
        integer, intent(in) :: howmany
        double precision, dimension(3, howmany) :: s
        double precision vecscra(3)
        integer i
        !   NB  is never in the GPU
        ! In principle one has to find the nearest image inside the Wigner-Seitz unit
        ! but there will be a difference as compared to the conventional distance below
        ! (the one valid for an ortho supercell) only at the boundary where all contributions
        ! vanish.
        do i = 1, howmany
            vecscra(:) = car2cry(:, 1)*s(1, i) + car2cry(:, 2)*s(2, i) + car2cry(:, 3)*s(3, i)
            !     vecscra(:)=s(:,i)
            !     call CartesianToCrystal(vecscra,1)
            vecscra(1) = anint(vecscra(1)/cellscale(1))
            vecscra(2) = anint(vecscra(2)/cellscale(2))
            vecscra(3) = anint(vecscra(3)/cellscale(3))
            s(:, i) = s(:, i) - s2r(:, 1)*vecscra(1) - s2r(:, 2)*vecscra(2) - s2r(:, 3)*vecscra(3)
            !     call dgemv('N',3,3,-1.d0,s2r,3,vecscra,1,1.d0,s(1,i),1)
        end do
    end subroutine ApplyPBC

    function map(x, cell_period)
        implicit none
        ! argument variables
        real(8), intent(in) :: x, cell_period
        ! local variables
        real(8) :: map
        if (cell_period .eq. 0.d0) then
            map = x
        else
            map = cell_period*map0(x/cell_period)
        end if
    end function map

    function dmap(x, cell_period)
        implicit none
        real(8) :: x, cell_period
        real(8) :: dmap
        !   dmap=dcos(x/cell_period)
        if (cell_period .eq. 0.d0) then
            dmap = 1.d0
        else
            dmap = dmap0(x/cell_period)
        end if
    end function dmap

    function ddmap(x, cell_period)
        implicit none
        real(8) :: x, cell_period
        real(8) :: ddmap
        if (cell_period .eq. 0.d0) then
            ddmap = 0.d0
        else
            ddmap = ddmap0(x/cell_period)/cell_period
        end if
    end function ddmap

    function map0(x)
        implicit none
        ! argument variables
        real(8), intent(in) :: x
        ! local variables
        real(8) xc, map0
        integer p
        ! this function depend only on x and is such that f'=1 and f(1/2)=0
        select case (case_map)
        case (0)
            map0 = sin(x*Pi)/Pi
        case (1)
            xc = x - anint(x) !   -1/2 < x < 1/2
            map0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. x_c) then
                map0 = xc
            elseif (xc .gt. x_c .and. xc .lt. 0.5d0) then
                map0 = amap/((0.5d0 - xc)**2*(bmap - xc))
            elseif (xc .gt. -0.5d0) then
                map0 = -amap/((0.5d0 + xc)**2*(bmap + xc))
            end if
        case (2)
            xc = x - anint(x) !   -1/2 < x < 1/2
            map0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. one_6) then
                map0 = xc
            elseif (xc .gt. one_6 .and. xc .lt. 0.5d0) then
                map0 = one_54/(0.5d0 - xc)**2
            elseif (xc .gt. -0.5d0) then
                map0 = -one_54/(0.5d0 + xc)**2
            end if
        case (3)
            map0 = sin(x*TWO_Pi)/TWO_Pi
        case (4)
            xc = x - anint(x) !   -1/2 < x < 1/2
            map0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                map0 = xc
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                map0 = 1.d0/(8.d0 - 16.d0*xc)
            elseif (xc .gt. -0.5d0) then
                map0 = -1.d0/(8.d0 + 16.d0*xc)
            end if
        case (5)
            xc = x - anint(x) !   -1/2 < x < 1/2
            map0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                map0 = xc
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                map0 = 1.d0/(-32.d0*(xc - 0.5d0) - 128d0*(xc - 0.5d0)**2 - 256d0*(xc - 0.5d0)**3)
            elseif (xc .gt. -0.5d0) then
                map0 = -1.d0/(-32.d0*(-xc - 0.5d0) - 128d0*(-xc - 0.5d0)**2 - 256d0*(-xc - 0.5d0)**3)
            end if
        case default
            p = case_map - 5
            xc = x - anint(x) !   -1/2 < x < 1/2
            map0 = 0.d0
            if (xc .ne. 0.d0 .and. abs(xc) .lt. 0.5d0) then
                map0 = (1.d0 - exp(-1.d0/xc**2 + 4.d0))**p ! protected from division by zero
                if (map0 .ne. 0.d0) map0 = xc/map0
            end if
        end select
    end function map0

    function dmap0(x)
        real(8) :: x, xc, dmap0, dummy, ddummy, dummy0
        integer p
        ! this function depend only on x and is such that f'=1 and f(1/2)=0
        select case (case_map)
        case (0)
            dmap0 = cos(x*Pi)
        case (1)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. x_c) then
                dmap0 = 1.d0
            elseif (xc .gt. x_c .and. xc .lt. 0.5d0) then
                dmap0 = amap/(0.5d0 - xc)**2/(bmap - xc)**2 + 2.d0*amap/(0.5d0 - xc)**3/(bmap - xc)
            elseif (xc .gt. -0.5d0) then
                dmap0 = amap/(0.5d0 + xc)**2/(bmap + xc)**2 + 2.d0*amap/(0.5d0 + xc)**3/(bmap + xc)
            end if
        case (2)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. one_6) then
                dmap0 = 1.d0
            elseif (xc .gt. one_6 .and. xc .lt. 0.5d0) then
                dmap0 = one_27/(0.5d0 - xc)**3
            elseif (xc .gt. -0.5d0) then
                dmap0 = one_27/(0.5d0 + xc)**3
            end if
        case (3)
            dmap0 = cos(x*TWO_Pi)
        case (4)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                dmap0 = 1.d0
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                dmap0 = 0.25d0/(1.d0 - 2.d0*xc)**2
            elseif (xc .gt. -0.5d0) then
                dmap0 = 0.25d0/(1.d0 + 2.d0*xc)**2
            end if
        case (5)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                dmap0 = 1.d0
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                dmap0 = -(-32.d0 - 256.d0*(xc - 0.5d0) - 768.d0*(xc - 0.5d0)**2)/&
                     &(-32.d0*(xc - 0.5d0) - 128d0*(xc - 0.5d0)**2 - 256d0*(xc - 0.5d0)**3)**2
            elseif (xc .gt. -0.5d0) then
                dmap0 = -(-32.d0 - 256.d0*(-xc - 0.5d0) - 768.d0*(-xc - 0.5d0)**2)/&
                     &(-32.d0*(-xc - 0.5d0) - 128d0*(-xc - 0.5d0)**2 - 256d0*(-xc - 0.5d0)**3)**2
            end if
        case default
            p = case_map - 5
            xc = x - anint(x) !   -1/2 < x < 1/2
            dmap0 = 0.d0
            if (xc .ne. 0.d0 .and. abs(xc) .lt. 0.5d0) then
                dummy0 = exp(-1.d0/xc**2 + 4.d0)
                dummy = (1.d0 - dummy0) ! protected from division by zero
                if (dummy .ne. 0.d0) then
                    ddummy = -2.d0*p/xc**3*dummy0/dummy
                    dmap0 = 1.0/dummy**p - xc/dummy**p*ddummy
                end if
            end if
        end select
    end function dmap0

    function ddmap0(x)
        real(8) :: x, xc, ddmap0, dummy, ddummy, d2dummy, dummy0
        integer p
        select case (case_map)
        case (0)
            ddmap0 = -Pi*sin(x*Pi)
        case (1)
            xc = x - anint(x) !   -1/2 < x < 1/2
            ddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. x_c) then
                ddmap0 = 0.d0
            elseif (xc .gt. x_c .and. xc .lt. 0.5d0) then
                ddmap0 = 2.d0*amap/(0.5d0 - xc)**2/(bmap - xc)**3&
                        & + 4.d0*amap/(0.5d0 - xc)**3/(bmap - xc)**2&
                        & + 6.d0*amap/(0.5d0 - xc)**4/(bmap - xc)
            elseif (xc .gt. -0.5d0) then
                ddmap0 = -2.d0*amap/(0.5d0 + xc)**2/(bmap + xc)**3&
                        & - 4.d0*amap/(0.5d0 + xc)**3/(bmap + xc)**2&
                        & - 6.d0*amap/(0.5d0 + xc)**4/(bmap + xc)
            end if
        case (2)
            xc = x - anint(x) !   -1/2 < x < 1/2
            ddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. one_6) then
                ddmap0 = 0.d0
            elseif (xc .gt. one_6 .and. xc .lt. 0.5d0) then
                ddmap0 = one_9/(0.5d0 - xc)**4
            elseif (xc .gt. -0.5d0) then
                ddmap0 = -one_9/(0.5d0 + xc)**4
            end if
        case (3)
            ddmap0 = -TWO_Pi*sin(x*TWO_Pi)
        case (4)
            xc = x - anint(x) !   -1/2 < x < 1/2
            ddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                ddmap0 = 0.d0
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                ddmap0 = 1.d0/(1.d0 - 2.d0*xc)**3
            elseif (xc .gt. -0.5d0) then
                ddmap0 = -1.d0/(1.d0 + 2.d0*xc)**3
            end if
        case (5)
            xc = x - anint(x) !   -1/2 < x < 1/2
            ddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                ddmap0 = 0.d0
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                ddmap0 = -0.5d0*(1.d0 - 4.d0*xc)**2*(5.d0 - 20.d0*xc + 24.d0*xc**2)/&
              &(-1.d0 + 6.d0*xc - 16.d0*xc**2 + 16.d0*xc**3)**3
            elseif (xc .gt. -0.5d0) then
                ddmap0 = 0.5d0*(1.d0 + 4.d0*xc)**2*(5.d0 + 20.d0*xc + 24.d0*xc**2)/&
              &(-1.d0 - 6.d0*xc - 16.d0*xc**2 - 16.d0*xc**3)**3
            elseif (xc .gt. -0.5d0) then
            end if
        case default
            p = case_map - 5
            xc = x - anint(x) !   -1/2 < x < 1/2
            ddmap0 = 0.d0
            if (xc .ne. 0.d0 .and. abs(xc) .lt. 0.5d0) then
                dummy0 = exp(-1.d0/xc**2 + 4.d0)
                dummy = (1.d0 - dummy0) ! protected from division by zero
                if (dummy .ne. 0.d0) then
                    ddummy = -2.d0*p/xc**3*dummy0/dummy
                    d2dummy = 2.d0*p*dummy0/dummy**2*(3.d0/xc**4*dummy + 2.d0/xc**6*(p - 1.d0 - p*dummy))
                    ddmap0 = -2.0/dummy**p*ddummy + 2.0*xc/dummy**p*ddummy**2 - xc/dummy**p*d2dummy
                end if
            end if
        end select
    end function ddmap0

    function dddmap0(x)
        real(8) :: x, xc, dddmap0, dummy0, dummy
        integer p
        ! this function depend only on x and is such that f'=1 and f(1/2)=0
        select case (case_map)
        case (0)
            dddmap0 = -Pi2*cos(x*Pi)
        case (1)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. x_c) then
                dddmap0 = 0.d0
            elseif (xc .gt. x_c .and. xc .lt. 0.5d0) then
                dddmap0 = 6.d0*amap/(0.5d0 - xc)**2/(bmap - xc)**4 &
                        & + 12.d0*amap/(0.5d0 - xc)**3/(bmap - xc)**3&
                        & + 18.d0*amap/(0.5d0 - xc)**4/(bmap - xc)**2&
                        & + 24.d0*amap/(0.5d0 - xc)**5/(bmap - xc)
            elseif (xc .gt. -0.5d0) then
                dddmap0 = 6.d0*amap/(0.5d0 + xc)**2/(bmap + xc)**4 &
                        & + 12.d0*amap/(0.5d0 + xc)**3/(bmap + xc)**3&
                        & + 18.d0*amap/(0.5d0 + xc)**4/(bmap + xc)**2&
                        & + 24.d0*amap/(0.5d0 + xc)**5/(bmap + xc)
            end if
        case (2)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. one_6) then
                dddmap0 = 0.d0
            elseif (xc .gt. one_6 .and. xc .lt. 0.5d0) then
                dddmap0 = four_9/(0.5d0 - xc)**5
            elseif (xc .gt. -0.5d0) then
                dddmap0 = four_9/(0.5d0 + xc)**5
            end if
        case (3)
            dddmap0 = -4.d0*Pi2*cos(x*TWO_Pi)
        case (4)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                dddmap0 = 0.d0
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                dddmap0 = 6.d0/(1.d0 - 2.d0*xc)**4
            elseif (xc .gt. -0.5d0) then
                dddmap0 = 6.d0/(1.d0 + 2.d0*xc)**4
            end if
        case (5)
            xc = x - anint(x) !   -1/2 < x < 1/2
            dddmap0 = 0.d0 ! if xc=+/- 1/2
            if (abs(xc) .le. 0.25d0) then
                dddmap0 = 0.d0
            elseif (xc .gt. 0.25d0 .and. xc .lt. 0.5d0) then
                dddmap0 = 3.d0*(5.d0 - 112.d0*xc + 928.d0*xc**2 - 3840.d0*xc**3 + 8640.d0*xc**4&
              &- 10240.d0*xc**5 + 5120.d0*xc**6)/&
              &(-1.d0 + 6.d0*xc - 16.d0*xc**2 + 16.d0*xc**3)**4
            elseif (xc .gt. -0.5d0) then
                dddmap0 = 3.d0*(5.d0 + 112.d0*xc + 928.d0*xc**2 + 3840.d0*xc**3 + 8640.d0*xc**4&
              &+ 10240.d0*xc**5 + 5120.d0*xc**6)/&
              &(-1.d0 - 6.d0*xc - 16.d0*xc**2 - 16.d0*xc**3)**4
            end if
        case default
            p = case_map - 5
            xc = x - anint(x) !   -1/2 < x < 1/2
            dddmap0 = 0.d0
            if (xc .ne. 0.d0 .and. abs(xc) .lt. 0.5d0) then
                dummy0 = exp(-1.d0/xc**2 + 4.d0)
                dummy = (1.d0 - dummy0) ! protected from division by zero
                if (dummy .ne. 0.d0) then
                    dddmap0 = (8.d0*dummy0*dummy**(-1 - p)*p)/xc**8 - &
                    &(24.d0*dummy0**2*dummy**(-2 - p)*(-1.d0 - p)*p)/xc**8 + &
                    &(8.d0*dummy0**3*dummy**(-3 - p)*(-2.d0 - p)*(-1.d0 - p)*p)/xc**8 - &
                    &(24.d0*dummy0*dummy**(-1 - p)*p)/xc**6 + &
                    &(24.d0*dummy0**2*dummy**(-2 - p)*(-1.d0 - p)*p)/xc**6 + &
                    &(6.d0*dummy0*dummy**(-1 - p)*p)/xc**4
                end if
            end if
        end select
    end function dddmap0

end module Cell

subroutine map0_b(x, xb, map0b)
    use Cell, only: dmap0
    implicit none
    real*8 x, xb, map0b
    xb = xb + dmap0(x)*map0b
    map0b = 0.d0
end subroutine map0_b

subroutine dmap0_b(x, xb, dmap0b)
    use Cell, only: ddmap0
    implicit none
    real*8 x, xb, dmap0b
    xb = xb + ddmap0(x)*dmap0b
    dmap0b = 0.d0
end subroutine dmap0_b

subroutine ddmap0_b(x, xb, ddmap0b)
    use Cell, only: dddmap0
    implicit none
    real*8 x, xb, ddmap0b
    xb = xb + dddmap0(x)*ddmap0b
    ddmap0b = 0.d0
end subroutine ddmap0_b

subroutine map_b(x, xb, cell_period, cell_periodb, mapb)
    use Cell, only: map0
    implicit none
    real*8 x, xb, y, yb, z, zb, cell_period, cell_periodb, mapb
    !  map=cell_period*map0(x/cell_period)
    !  y=x/cell_period
    !  z=map0(y)
    !  map=cell_period*z
    if (cell_period .eq. 0.d0) then
        xb = xb + mapb
    else
        y = x/cell_period
        z = map0(y)
        cell_periodb = cell_periodb + z*mapb
        zb = cell_period*mapb
        mapb = 0.d0
        yb = 0.d0
        call map0_b(y, yb, zb)
        xb = xb + yb/cell_period
        cell_periodb = cell_periodb - yb*x/cell_period**2
    end if
    mapb = 0.d0
end subroutine map_b

subroutine dmap_b(x, xb, cell_period, cell_periodb, dmapb)
    implicit none
    real*8 x, xb, y, yb, cell_period, cell_periodb, dmapb
    !  dmap=dmap0(x/cell_period)
    !  y=x/cell_period
    !  dmap=dmap0(y)
    if (cell_period .ne. 0.d0) then
        yb = 0.d0
        y = x/cell_period
        call dmap0_b(y, yb, dmapb)
        xb = xb + yb/cell_period
        cell_periodb = cell_periodb - yb*x/cell_period**2
    end if
    dmapb = 0.d0
end subroutine dmap_b

subroutine ddmap_b(x, xb, cell_period, cell_periodb, ddmapb)
    use Cell, only: ddmap0
    implicit none
    real*8 x, xb, y, yb, z, zb, cell_period, cell_periodb, ddmapb
    !  dmap=ddmap0(x/cell_period)/cell_period
    !  y=x/cell_period
    !  z=ddmap0(y)
    !  ddmap=z/cell_period
    if (cell_period .ne. 0.d0) then
        yb = 0.d0
        y = x/cell_period
        z = ddmap0(y)
        zb = ddmapb/cell_period
        cell_periodb = cell_periodb - ddmapb*z/cell_period**2
        ddmapb = 0.d0
        call ddmap0_b(y, yb, zb)
        xb = xb + yb/cell_period
        cell_periodb = cell_periodb - yb*x/cell_period**2
    end if
    ddmapb = 0.d0
end subroutine ddmap_b
