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

!#include "mathvec.h"

module kpoints_mod

    use cell

    public
    logical yes_kpoints, same_phase, opposite_phase ! flag to switch on multiple read of fort.10's
    integer nk, nk1, nk2, nk3, k1, k2, k3, nks ! number of k-points and corresponding offsets for MP grid
    ! nk = total number of k-points
    ! nks = number of inequivalent kpoints in MP algorithm

    ! xkp,wkp = final k-points/weights for up/down electrons
    real(8), dimension(:, :), allocatable :: xkp, xkp_down
    real(8), dimension(:), allocatable :: wkp, wkp_down
    real(8) tot_wt, tot_wt_down

    integer kp_type ! methods used to generate k-points:
    ! 0 = gamma point or single phase calculation --> the phase is read from fort.10
    ! 1 = Monkhorst-Pack algorithm.
    ! 2 = read k-points in crystal coordinates from datasmin.input
    ! 3 = read k-points path in 1BZ in crystal coordinates as in QuantumESPRESSO.
    !     For band structure calculations. Useful also when searching special k-points
    !     along a diagonal direction in the 1BZ.
    ! 4 = generate random k-points.
    ! 5 = generate a nk*nk matrix of k-points for up/down spin electrons. Application of
    !     flavor twisted boundary conditions ( PRB 80, 180402(R) ). The matrix is generated
    !     by combining a normal and an offsetted MP grids. Alternatively, use kp_type=2 to give
    !     the matrix in input.

    integer :: ikpoint ! index which identify the k-point in the current pool ( = rankcolrep + 1)

    logical :: kaverage ! if .true. the twist average algorithm is turned ON.

    logical :: skip_equivalence ! if .false. the number of k-points is reduced using Bravais lattice symmetries

    logical :: time_reversal ! if .true. use the time reversal symmetry to reduce the number of k-points.

    logical :: double_kpgrid ! if .true. impose opposite boundary condition for UP/DOWN spin electrons.

    logical :: decoupled_run ! if .true. the code performs decoupled (i.e. indipendent k-points) VMC/DMC k-points calculations

    logical :: compute_bands ! flag determining whether the calculation is a band structure
    ! calculation or not (working for DFT only so far)

    private :: init_random_seed, found_sec, kpoint_grid, generate_kp_path

    ! namelist for k-points grid definition
    namelist /kpoints/ kp_type, nk1, nk2, nk3, k1, k2, k3, time_reversal, &
        skip_equivalence, double_kpgrid, compute_bands

contains

    !---------------------------------------------------------
    subroutine get_kpoints(iflag, nel, nion, rion, atom_number, rs)
        !---------------------------------------------------------
        !
        ! This routine initialize k-points for all kp_type
        !!!!!!! NB !!!!!!!
        ! k-points are left in crystal coordinates. Every subroutine
        ! which needs k-points must transform them back to carthesian
        ! coordinates (multiply by 2*pi/cell_scale)
        !

        use constants, only: PI

        implicit none

        integer :: iflag, nel, nion
        real(8) :: rs, t_comp, alphap, betap, gammap
        integer :: i, j, ikp, nks_up, nks_do
        real(8), dimension(:), allocatable :: wkp_, wkp_down_
        real(8), dimension(:, :), allocatable :: xkp_, xkp_down_
        real(8), external :: ran
        real*8, intent(in) :: rion(3, *), atom_number(*)
        logical smerding

        ! initialize cell in any case
        celldm(1) = 1.d0
        !  celldm(4:6) = 90d0*PI/180.d0

        omega = celldm(2)*celldm(3)

        !    omega=celldm(2)*celldm(3)
        celldm(1) = (PI*nel*4.d0/3.d0/omega)**(1.d0/3.d0)*rs

        !  omega=celldm(2)*celldm(3)
        !  celldm(1)=(PI*nel*4.d0/3.d0/omega)**(1.d0/3.d0)*rs
        ! setting up direct and reciprocal lattices
        if (yes_tilted) then
            givens2r = .true.
        else
            givens2r = .false.
        end if
        call InitCell(nion, nel, .true.)
        ! compute Bravais lattice symmetry operations
        call set_sym_bl(s2r)
        call purge_isymm(.true., nrot, isymm, sname, t_rev, rion, nion, atom_number, cellscale)
        write (6, *) 'number of point symmetries found', nrot
        write (6, *) ' Symmetries types found:'
        do i = 1, nrot
            write (6, *) sname(i)
            write (6, *) (isymm(j, 1, i), j=1, 3)
            write (6, *) (isymm(j, 2, i), j=1, 3)
            write (6, *) (isymm(j, 3, i), j=1, 3)
        end do
        smerding = .false.
        if (kp_type .lt. 0) then
            smerding = .true.
            kp_type = -kp_type
        end if
        !
        if (kp_type .eq. 0) then
            xkp(:, 1) = phase(:)
            xkp_down(:, 1) = phase_down(:)
            wkp(:) = 1.d0
            wkp_down(:) = 1.d0

        elseif (kp_type .eq. 1) then

            call kpoint_grid(nrot, time_reversal, skip_equivalence, isymm, t_rev, &
                             nk, k1, k2, k3, nk1, nk2, nk3, nks, xkp, wkp)
            nk = nks ! define the right number of k-points
            !
            ! in the case of an automatic grid, I restrict the choice of the
            ! down spin electrons to the opposite of up spin ones. This
            ! should probably be generalized.

            !
            if (opposite_phase) then
                xkp_down = -xkp
            else
                xkp_down = xkp
            end if
            wkp_down = wkp

        elseif (kp_type .eq. 2 .or. kp_type .eq. 3) then
            if (found_sec(5, 'KPOINTS', iflag)) then
                iflag = 1
                do i = 1, nk1
                    read (5, *, err=121) xkp(1, i), xkp(2, i), xkp(3, i), wkp(i)
                end do
                iflag = 0
                if (double_kpgrid) then
                    read (5, *) ! one blank
                    do i = 1, nk1
                        read (5, *, err=121) xkp_down(1, i), xkp_down(2, i), xkp_down(3, i), wkp_down(i)
                    end do
                else
                    if (opposite_phase) then
                        xkp_down = -xkp
                    else
                        xkp_down = xkp
                    end if
                    wkp_down = wkp
                end if
            end if
            !
            if (kp_type .eq. 3) then
                call generate_kp_path(xkp, wkp, nk1, nk2, iflag)
                if (double_kpgrid) then
                    call generate_kp_path(xkp_down, wkp_down, nk1, nk2, iflag)
                else
                    if (opposite_phase) then
                        xkp_down = -xkp
                    else
                        xkp_down = xkp
                    end if
                end if
                wkp(:) = 1.d0/nk
                wkp_down(:) = wkp(:)
            end if

        elseif (kp_type .eq. 4) then
            call init_random_seed
            nk = nk1 ! # of k-points
            wkp(:) = 1.d0/nk
            wkp_down(:) = wkp(:)
            do i = 1, nk
                do j = 1, 3
                    call random_number(t_comp)
                    t_comp = t_comp - 0.5d0
                    xkp(j, i) = t_comp
                    if (opposite_phase) then
                        xkp_down(j, i) = -t_comp
                    else
                        xkp_down(j, i) = t_comp
                    end if
                end do
            end do

        elseif (kp_type .eq. 5) then

            nk = nint(sqrt(dble(nk)))

            allocate (xkp_(3, nk), wkp_(nk))
            allocate (xkp_down_(3, nk), wkp_down_(nk))
            xkp_ = 0.d0
            wkp_ = 0.d0
            xkp_down_ = 0.d0
            wkp_down_ = 0.d0

            ! combine MP grids with and without offset to create
            ! the k-points matrix needed by flavor twist

            ! no offset
            call kpoint_grid(nrot, time_reversal, skip_equivalence, isymm, t_rev, &
                             nk, 0, 0, 0, nk1, nk2, nk3, nks, xkp_, wkp_)
            nks_up = nks

            ! with offset
            call kpoint_grid(nrot, time_reversal, skip_equivalence, isymm, t_rev, &
                             nk, 1, 1, 1, nk1, nk2, nk3, nks, xkp_down_, wkp_down_)
            nks_do = nks

            nk = nks_up*nks_do ! now I have a k-point matrix in the case of flavor twist

            ikp = 1
            do i = 1, nks_up
                do j = 1, nks_do
                    xkp(:, ikp) = xkp_(:, i)
                    xkp_down(:, ikp) = xkp_down_(:, j)
                    wkp(ikp) = wkp_(i)*wkp_down_(j)
                    wkp_down(ikp) = wkp_(i)*wkp_down_(j)
                    ikp = ikp + 1
                end do
            end do

            deallocate (xkp_, wkp_, xkp_down_, wkp_down_)

        end if
        if (smerding) then
            write (6, *) ' Warning adding a small pseudo random twist &
                    & to remove degeneracies '
            !     adding a pseudo random small contribution to remove all possible
            !     degeneracy of the eigenvalues
            xkp(1, :) = xkp(1, :) + 1.d-8/dsqrt(5.d0)
            xkp(2, :) = xkp(2, :) - 1.d-8/dsqrt(7.d0)
            xkp(3, :) = xkp(3, :) - 1.d-8/dsqrt(11.d0)
            if (opposite_phase) then
                xkp_down(1, :) = xkp_down(1, :) - 1.d-8/dsqrt(5.d0)
                xkp_down(2, :) = xkp_down(2, :) + 1.d-8/dsqrt(7.d0)
                xkp_down(3, :) = xkp_down(3, :) + 1.d-8/dsqrt(11.d0)
            else
                xkp_down(1, :) = xkp_down(1, :) + 1.d-8/dsqrt(5.d0)
                xkp_down(2, :) = xkp_down(2, :) - 1.d-8/dsqrt(7.d0)
                xkp_down(3, :) = xkp_down(3, :) - 1.d-8/dsqrt(11.d0)
            end if
        end if
121     return

    end subroutine get_kpoints

    !----------------------------------------------------------------------------
    subroutine kpoint_grid(nrot, time_reversal, skip_equivalence, s, t_rev, &
                           nkp, k1, k2, k3, nk1, nk2, nk3, nks, xk, wk)
        !----------------------------------------------------------------------------
        !
        !  Automatic generation of a uniform grid of k-points with MP algorithm.
        !  Ported from QuantumESPRESSO.
        !
        implicit none
        !
        integer, intent(in) :: nrot, k1, k2, k3, nk1, nk2, nk3, &
                               t_rev(48), s(3, 3, 48)
        logical, intent(in) :: time_reversal, skip_equivalence
        !
        integer, intent(inout) :: nks, nkp
        real(8), intent(out) :: xk(3, nkp)
        real(8), intent(out) :: wk(nkp)
        ! LOCAL:
        real(8), parameter :: eps = 1.0d-5
        real(8) :: xkr(3), fact, xx, yy, zz
        real(8), allocatable :: xkg(:, :), wkk(:)
        integer :: nkr, i, j, k, ns, n, nk
        integer, allocatable :: equiv(:)
        logical :: in_the_list
        !
        nkr = nk1*nk2*nk3
        allocate (xkg(3, nkr), wkk(nkr))
        allocate (equiv(nkr))
        !
        do i = 1, nk1
            do j = 1, nk2
                do k = 1, nk3
                    !  this is nothing but consecutive ordering
                    n = (k - 1) + (j - 1)*nk3 + (i - 1)*nk2*nk3 + 1
                    !  xkg are the components of the complete grid in crystal axis
                    xkg(1, n) = dble(i - 1)/nk1 + dble(k1)/2/nk1
                    xkg(2, n) = dble(j - 1)/nk2 + dble(k2)/2/nk2
                    xkg(3, n) = dble(k - 1)/nk3 + dble(k3)/2/nk3
                end do
            end do
        end do

        !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
        !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

        do nk = 1, nkr
            equiv(nk) = nk
        end do

        if (skip_equivalence) then
            write (6, *) ' Warning k-points grid: skip check of k-points equivalence '
            wkk = 1.d0
        else
            do nk = 1, nkr
                !  check if this k-point has already been found equivalent to another
                if (equiv(nk) == nk) then
                    wkk(nk) = 1.0d0
                    !  check if there are equivalent k-point to this in the list
                    !  (excepted those previously found to be equivalent to another)
                    !  check both k and -k
                    do ns = 1, nrot
                        do i = 1, 3
                            xkr(i) = s(i, 1, ns)*xkg(1, nk) &
                                     + s(i, 2, ns)*xkg(2, nk) &
                                     + s(i, 3, ns)*xkg(3, nk)
                            xkr(i) = xkr(i) - nint(xkr(i))
                        end do
                        if (t_rev(ns) == 1) xkr = -xkr
                        xx = xkr(1)*nk1 - 0.5d0*k1
                        yy = xkr(2)*nk2 - 0.5d0*k2
                        zz = xkr(3)*nk3 - 0.5d0*k3
                        in_the_list = abs(xx - nint(xx)) <= eps .and. &
                                      abs(yy - nint(yy)) <= eps .and. &
                                      abs(zz - nint(zz)) <= eps
                        if (in_the_list) then
                            i = mod(nint(xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1) + 1
                            j = mod(nint(xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2) + 1
                            k = mod(nint(xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3) + 1
                            n = (k - 1) + (j - 1)*nk3 + (i - 1)*nk2*nk3 + 1
                            if (n > nk .and. equiv(n) == n) then
                                equiv(n) = nk
                                wkk(nk) = wkk(nk) + 1.0d0
                            end if
                        end if
                        if (time_reversal) then
                            xx = -xkr(1)*nk1 - 0.5d0*k1
                            yy = -xkr(2)*nk2 - 0.5d0*k2
                            zz = -xkr(3)*nk3 - 0.5d0*k3
                            in_the_list = abs(xx - nint(xx)) <= eps .and. abs(yy - nint(yy)) <= eps &
                                          .and. abs(zz - nint(zz)) <= eps
                            if (in_the_list) then
                                i = mod(nint(-xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1) + 1
                                j = mod(nint(-xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2) + 1
                                k = mod(nint(-xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3) + 1
                                n = (k - 1) + (j - 1)*nk3 + (i - 1)*nk2*nk3 + 1
                                if (n > nk .and. equiv(n) == n) then
                                    equiv(n) = nk
                                    wkk(nk) = wkk(nk) + 1.0d0
                                end if
                            end if
                        end if
                    end do
                end if
            end do
        end if

        !  count irreducible points and order them
        nks = 0
        fact = 0.0d0
        do nk = 1, nkr
            if (equiv(nk) == nk) then
                nks = nks + 1
                wk(nks) = wkk(nk)
                fact = fact + wk(nks)
                !  bring back into to the first BZ
                do i = 1, 3
                    xk(i, nks) = xkg(i, nk) - nint(xkg(i, nk))
                end do
            end if
        end do
        !  normalize weights to one
        do nk = 1, nks
            wk(nk) = wk(nk)/fact
        end do

        write (6, *) 'number of inequivalent/total k-points found', nks, nkp
        ! set total number of k-points as the number of inequivalent ones
        nkp = nks

        deallocate (equiv)
        deallocate (xkg, wkk)

        return

    end subroutine kpoint_grid

    ! ----------------------------------------
    subroutine generate_kp_path(xkp, wkp, nk1, nk2, iflag)
        ! ----------------------------------------
        !
        !  This subroutine generates a k-points path in the 1BZ according
        !  to the points read in the section KPOINTS.

        implicit none

        integer, intent(in) :: nk1, nk2
        integer, intent(inout) :: iflag
        real(8), intent(out) :: xkp(3, nk), wkp(nk)
        real(8), dimension(:, :), allocatable :: klines
        integer :: i, j, k
        real(8) :: dk, t, xpar(3), xold(3), wold, xmod

        ! read lines of the path
        allocate (klines(3, nk1))
        do i = 1, nk1
            klines(:, i) = xkp(:, i)
        end do
        xkp = 0.d0
        wkp = 0.d0
        dk = 1.d0/nk2
        iflag = 0
        ! generates k-points along the lines
        ! using parametrization
        do i = 1, nk1 - 1
            t = 1.d0
            do j = 1, nk2
                if (t .lt. 0) then
                    iflag = 1
                    go to 123
                end if

                if (j .eq. 1) then
                    xold(:) = klines(:, i)
                else
                    xold(:) = xkp(:, (i - 1)*nk2 + j - 1)
                end if

                wold = wkp((i - 1)*nk2 + j)
                xpar(:) = t*(klines(:, i) - klines(:, i + 1)) + klines(:, i + 1)
                xkp(:, (i - 1)*nk2 + j) = xpar(:)
                xmod = dsqrt(sum((xpar(:) - xold(:))**2))
                wkp((i - 1)*nk2 + j + 1) = wold + xmod
                t = t - dk

            end do
        end do
        ! add the last point
        xkp(:, (nk1 - 1)*nk2 + 1) = klines(:, nk1)
        !
        !
123     if (iflag .ne. 0) write (6, *) 'ERROR in generate k-points path'
        deallocate (klines)

        return

    end subroutine generate_kp_path

    subroutine check_kpoints(id, comm, rank)

        implicit none
        integer, intent(in) :: id, comm, rank
        integer :: i, indk
        real(8), dimension(:, :), allocatable :: all_phases, all_phases_down
        real(8) :: phase_diff, phase_diff_down
        logical :: is_same(3), is_same_do(3)

        allocate (all_phases(3, nk), all_phases_down(3, nk))
        all_phases = 0.d0
        all_phases_down = 0.d0
        indk = id + 1

        all_phases(:, indk) = phase(:)
        all_phases_down(:, indk) = phase_down(:)

        call reduce_base_real(size(all_phases), all_phases, comm, -1)
        call reduce_base_real(size(all_phases_down), all_phases_down, comm, -1)

        do i = 1, nk

            is_same(:) = ((abs(xkp(:, i) - all_phases(:, i))) < 1.0d-5)
            is_same_do(:) = ((abs(xkp_down(:, i) - all_phases_down(:, i))) < 1.0d-5)

            if ((.not. all(is_same)) .or. (.not. all(is_same_do))) then
                phase_diff = sum(abs(xkp(:, i) - all_phases(:, i)))
                phase_diff_down = sum(abs(xkp_down(:, i) - all_phases_down(:, i)))
                call error('kpoints', ' Some wavefunction phases do not correspond to the &
                        & input k-points! Check &kpoints card of fort.10_000xxx ', 1, rank)
            end if

        end do

        deallocate (all_phases, all_phases_down)

        return
    end subroutine check_kpoints
    !
    ! weighted sums of over pools for different variables types
    ! To call those functions use the interface above!
    !
    subroutine sum_kpoints_array_int(ps, dim_ps, comm, root)
        implicit none
        integer, intent(in) :: dim_ps, comm, root
        integer, intent(inout) :: ps(*)
        integer :: ierr
        integer, dimension(:), allocatable :: ps_tmp
#ifdef PARALLEL
        include "mpif.h"
#endif
        if ((.not. kaverage) .or. decoupled_run) return
        allocate (ps_tmp(dim_ps))
        ps_tmp(1:dim_ps) = ps(1:dim_ps)
        ps_tmp(1:dim_ps) = ps_tmp(1:dim_ps)*wkp(ikpoint)
#ifdef PARALLEL
        call mpi_allreduce(ps_tmp, ps, dim_ps, MPI_INTEGER, MPI_SUM, comm, ierr)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
        deallocate (ps_tmp)
        return
    end subroutine sum_kpoints_array_int

    subroutine sum_kpoints_scalar_int(ps, comm, root)
        implicit none
        integer, intent(in) :: comm, root
        integer, intent(inout) :: ps
        integer :: ps_tmp, ierr
#ifdef PARALLEL
        include "mpif.h"
#endif
        if ((.not. kaverage) .or. decoupled_run) return
        ps_tmp = ps
        ps_tmp = ps_tmp*wkp(ikpoint)
#ifdef PARALLEL
        call mpi_allreduce(ps_tmp, ps, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
        return
    end subroutine sum_kpoints_scalar_int

    subroutine sum_kpoints_scalar_real8(ps, comm, root)
        implicit none
        integer, intent(in) :: root, comm
        real(8), intent(inout) :: ps
        integer :: ierr
#ifdef PARALLEL
        include "mpif.h"
#endif
        if ((.not. kaverage) .or. decoupled_run) return
        ps = ps*wkp(ikpoint)
#ifdef PARALLEL
        call reduce_base_real(1, ps, comm, root)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
        return
    end subroutine sum_kpoints_scalar_real8

    subroutine sum_kpoints_array_real8(ps, dim_ps, comm, root)
        implicit none
        integer, intent(in) :: dim_ps, root, comm
        real(8), intent(inout) :: ps(*)
        integer :: ierr
#ifdef PARALLEL
        include "mpif.h"
#endif
        if ((.not. kaverage) .or. decoupled_run) return
        ps(1:dim_ps) = ps(1:dim_ps)*wkp(ikpoint)
#ifdef PARALLEL
        call reduce_base_real(dim_ps, ps, comm, root)
#endif
        return
    end subroutine sum_kpoints_array_real8

    subroutine sum_kpoints_array_real16(ps, dim_ps, comm, root)
        implicit none
        integer, intent(in) :: dim_ps, comm, root
        real(qp), intent(inout) :: ps(*)
        real*8, dimension(:), allocatable :: ps_tmp, ps_double
        integer :: ierr
#ifdef PARALLEL
        include "mpif.h"
#endif
        if ((.not. kaverage) .or. decoupled_run) return
        allocate (ps_tmp(dim_ps), ps_double(dim_ps))
        ps_tmp(1:dim_ps) = ps(1:dim_ps)*wkp(ikpoint)
#ifdef PARALLEL
        call mpi_allreduce(ps_tmp, ps_double, dim_ps, MPI_REAL8, MPI_SUM, comm, ierr)
        ps(1:dim_ps) = ps_double(1:dim_ps)
#else
        ps(1:dim_ps) = ps_tmp(1:dim_ps)
#endif
        deallocate (ps_tmp, ps_double)
        return
    end subroutine sum_kpoints_array_real16

    subroutine sum_kpoints_scalar_real16(ps, comm, root)
        implicit none
        integer, intent(in) :: comm, root
        real(qp), intent(inout) :: ps
        real*8 :: ps_tmp, ps_double
        integer :: ierr
#ifdef PARALLEL
        include "mpif.h"
#endif
        if ((.not. kaverage) .or. decoupled_run) return
        ps_tmp = ps*wkp(ikpoint)
#ifdef PARALLEL
        call mpi_allreduce(ps_tmp, ps_double, 1, MPI_REAL8, MPI_SUM, comm, ierr)
        ps = ps_double
#else
        ps = ps_tmp
#endif
        return
    end subroutine sum_kpoints_scalar_real16

    ! -----------------------------------------------
    logical function found_sec(funit, section_name, iflag)
        ! -----------------------------------------------
        !
        ! This function search in the input file for section "section_name".
        ! If the section is not found then iflag is set to one.

        implicit none

        integer funit, iflag
        character(len=*), intent(in) :: section_name
        character(20) :: line

        found_sec = .false.

        do while (.not. found_sec)
            read (funit, '(a)', end=122) line
            !
            if (trim(line) == trim(section_name)) found_sec = .true.
            !
        end do
        return
        ! error
122     write (6, *) 'section  ', section_name, '  not found'
        iflag = 1
        return

    end function found_sec
    !
    ! initialize seed for random k-points generation
    ! based on system time.
    !
    subroutine init_random_seed
        implicit none
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
        call random_seed(size=n)
        allocate (seed(n))
        call system_clock(count=clock)
        seed = clock + 37*(/(i - 1, i=1, n)/)
        call random_seed(put=seed)
        deallocate (seed)
        return
    end subroutine init_random_seed

end module kpoints_mod
