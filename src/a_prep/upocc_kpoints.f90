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

!----------------------------------------
subroutine upocc_kpoints(optocc, dt, iopt)
!----------------------------------------

    ! This subroutine computes Fermi level and KS orbitals occupations for both
    ! fixed occupations (optocc=0) and with smearing (optocc=1).
    ! Fermi energy is evaluated by the master and then broadcast to the
    ! processors pools. This is a safe choice to avoid too much
    ! communication during iterative bisection method.
    !
    ! iopt = 1 -> used for printing initial occupations
    ! iopt = 0 -> default option for SC iterations

    use constants, only: ipc, qp
    use allio, only: rank, nelup, neldo, commcolrep_mpi, rankrep, decoupled_run
    use kpoints_mod, only: nk, wkp, xkp, xkp_down, wkp_down, sum_kpoints_scalar_real16, compute_bands
    use setup, only: yeslsda, bands, indk, eigmol, eigmoldo, occupations, &
                     occupationdo, eigmol_sav, eigmoldo_sav, occupations_sav, &
                     occupationsdo_sav, nelocc, neloccdo, oldocc, oldoccdo, print_energies
    use freeelmod_complex, only: count2, costalln, countall, costall
    use parallel_module, only: collect_from_pools
    use compute_efermi

    implicit none

#ifdef PARALLEL
    include 'mpif.h'
#endif

    integer, intent(in) :: optocc, iopt
    real(qp), intent(in) :: dt
    !
    integer :: ierr, ibnd, ikp
    integer, dimension(:), allocatable :: number_particles_, number_particles_do_
    logical :: double_occ
    real*8 cost_double
    if (compute_bands) return
    !
    ! collect eigeinvalues and occupations from pools
    ! for Efermi computation.
    call collect_from_pools
    !
    ! compute Fermi energies
    !
    double_occ = .false.
    if (yeslsda .or. ipc .eq. 2) double_occ = .true.
    !
    efermi = 0.0_qp
    ef_up = 0.0_qp
    ef_do = 0.0_qp
    smear_type = optocc

    if (.not. decoupled_run) then
        if (rank .eq. 0) then
            if (optocc .eq. 0) then
                call efermiFix(eigmol_sav, occupations_sav, bands, ef_up)
                if (double_occ) call efermiFix(eigmoldo_sav, occupationsdo_sav, bands, ef_do)
            else
                ! spin up electrons
                call efermiMet(eigmol_sav, bands, nelup, dt, ef_up)
                ! spin down electrons
                if (double_occ) then
                    call efermiMet(eigmoldo_sav, bands, neldo, dt, ef_do)
                else
                    call efermiMet(eigmol_sav, bands, neldo, dt, ef_do)
                end if
            end if
            efermi = ef_up ! assuming nelup > neldo
        end if
#if defined PARALLEL
        cost_double = efermi
        call mpi_bcast(cost_double, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        efermi = cost_double
        cost_double = ef_up
        call mpi_bcast(cost_double, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        ef_up = cost_double
        cost_double = ef_do
        call mpi_bcast(cost_double, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        ef_do = cost_double
#endif
    else
        if (optocc .eq. 0) then
            call efermiFix(eigmol, occupations, bands, ef_up)
            if (double_occ) call efermiFix(eigmoldo, occupationdo, bands, ef_do)
        else
            ! spin up electrons
            call efermiMet(eigmol, bands, nelup, dt, ef_up)
            ! spin down electrons
            if (double_occ) then
                call efermiMet(eigmoldo, bands, neldo, dt, ef_do)
            else
                call efermiMet(eigmol, bands, neldo, dt, ef_do)
            end if
        end if
        ! here no broadcast is needed since the Fermi energy is kept
        ! independent for each k-point!
        efermi = ef_up
    end if
    !
    ! compute number of particles for each k-point
    !
    if (print_energies) then
        allocate (number_particles_(nk))
        number_particles_ = 0
        number_particles = 0
        if (double_occ) then
            allocate (number_particles_do_(nk))
            number_particles_do_ = 0
            number_particles_do = 0
        end if
    end if
    !
    ! update occupations with smearing
    !
    if (optocc .ne. 0) then

        metS = 0.0_qp
        metSdo = 0.0_qp
        occupations = 0.d0
        if (double_occ) occupationdo = 0.d0

        do ibnd = 1, bands
            occupations(ibnd) = smear(eigmol(ibnd), dt, ef_up, smear_type)
            metS = metS + dt*smearS(eigmol(ibnd), dt, ef_up, smear_type)
            if (double_occ) then
                occupationdo(ibnd) = smear(eigmoldo(ibnd), dt, ef_do, smear_type)
                metSdo = metSdo + dt*smearS(eigmoldo(ibnd), dt, ef_do, smear_type)
            else
                occupations(ibnd) = occupations(ibnd) + smear(eigmol(ibnd), dt, ef_do, smear_type)
            end if
        end do
        if (rankrep .eq. 0) then
            number_particles_(indk) = nint(sum(occupations(1:bands)))
            if (double_occ) number_particles_do_(indk) = nint(sum(occupationdo(1:bands)))
        end if

        call sum_kpoints_scalar_real16(metS, commcolrep_mpi, -1)
        if (double_occ) call sum_kpoints_scalar_real16(metSdo, commcolrep_mpi, -1)
        !    metS = metS + metSdo
        ! When a Fermi-Dirac smearing is used
        ! for a given dt, the variational energy is really the free energy F,
        ! and F = E - TS , with E = \sum_i eig_i, -TS = metS
        metS = -metS
        metSdo = -metSdo

#ifdef PARALLEL
        if (print_energies) then
            call mpi_allreduce(number_particles_, number_particles, &
                               nk, MPI_INTEGER, MPI_SUM, commcolrep_mpi, ierr)
            if (double_occ) &
                call mpi_allreduce(number_particles_do_, number_particles_do, &
                                   nk, MPI_INTEGER, MPI_SUM, commcolrep_mpi, ierr)
        end if
#endif

    else
        number_particles = nint(sum(occupations(1:bands)))
        if (double_occ) number_particles_do = nint(sum(occupationdo(1:bands)))
    end if

    if (print_energies) then
        deallocate (number_particles_)
        if (double_occ) deallocate (number_particles_do_)
    end if

    ! print output
    call collect_from_pools

    if (rank .eq. 0) then
        call check_occupation
        if (iopt .eq. 1) then
            write (6, *)
            write (6, '(a)') ' Eigenvalues/occupations before starting SC cycle '
            if (double_occ) then
                do ikp = 1, nk
                    write (6, *)
                    write (6, '(a,3F12.8)') ' # k-point spin up:  ', xkp(:, ikp)
                    do ibnd = 1, bands
                        write (6, 100) ibnd, eigmol_sav(ibnd, ikp), occupations_sav(ibnd, ikp)
                    end do
                    write (6, '(a,3F12.8)') ' # k-point spin down:  ', xkp_down(:, ikp)
                    do ibnd = 1, bands
                        write (6, 100) ibnd, eigmoldo_sav(ibnd, ikp), occupationsdo_sav(ibnd, ikp)
                    end do
                end do
            else
                do ikp = 1, nk
                    write (6, *)
                    write (6, '(a,3F12.8)') ' # k-point:  ', xkp(:, ikp)
                    do ibnd = 1, bands
                        write (6, 100) ibnd, eigmol_sav(ibnd, ikp), occupations_sav(ibnd, ikp)
                    end do
                end do
            end if
        end if
        write (6, *)
        if (double_occ) then
            write (6, 101) ' Fermi energy up/down:         ', ef_up, ef_do
            write (6, 101) ' Entropic contribution up/down:', metS, metSdo
        else
            write (6, 101) ' Fermi energy:         ', ef_up
            write (6, 101) ' Entropic contribution:', metS
        end if
    end if
100 format(3x, I10, X, 2f20.10)
101 format(3x, a, 2f16.8)

    if (iopt .eq. 1) then
        if (optocc .ge. 1) then
            countall = bands
            costall = costalln
            oldocc = occupations
            if (double_occ) oldoccdo = occupationdo
            nelocc = bands
            neloccdo = bands
            count2 = 0
            do ibnd = 1, bands
                if (eigmol(ibnd) .le. efermi) count2 = count2 + 1
            end do
        else
            oldocc = occupations
            if (double_occ) oldoccdo = occupationdo
        end if
    end if

    return

contains

    subroutine check_occupation
        ! routine to check if constraint \sum_i nocc_i = N_el is satisfied.
        implicit none
        integer :: ikp
        real(8) :: checkocc, checkoccdo, sum1
        checkocc = 0.d0
        checkoccdo = 0.d0
        do ikp = 1, nk
            sum1 = sum(occupations_sav(1:bands, ikp))
            checkocc = checkocc + wkp(ikp)*sum1
            if (double_occ) then
                sum1 = sum(occupationsdo_sav(1:bands, ikp))
                checkoccdo = checkoccdo + wkp_down(ikp)*sum1
            end if
        end do
        write (6, '(a,2F20.10)') ' Total occupations found up/down:', checkocc, checkoccdo
        return
    end subroutine check_occupation

    !----------------------------
end subroutine upocc_kpoints
!----------------------------

!-----------------------------------------------------------------------------
subroutine read_occupations(occupations, occupationdo, nelocc, neloccdo, occread)
    !-----------------------------------------------------------------------------

    use constants, only: ipc
    use allio, only: rank, symmagp, neldo, rankcolrep, compute_bands
    use setup, only: optocc, bands, yeslsda, occupations_sav, occupationsdo_sav
    use freeelmod_complex, only: count2, costall, lastpaired, nmoltot, countall

    implicit none
    logical, intent(in) :: occread
    integer, intent(inout) :: nelocc, neloccdo
    real(8), intent(inout) :: occupations(*), occupationdo(*)
    integer i, ierr, iflag, iflagall, indk
#ifdef PARALLEL
    include "mpif.h"
#endif

    iflag = 0
    iflagall = 0

    if (occread) then

        if (rank .eq. 0) then

            iflag = 1
            read (5, *, err=100, end=100) (occupations(i), i=1, nelocc)
            iflag = 0
100         if (iflag .eq. 1) then
                call error(' read_occupations ', ' Error reading spin up occupations ', -1, rank)
                iflagall = iflagall + 1
            end if

            ! real wavefunction
            if (ipc .eq. 1) then
                if (yeslsda) then
                    if (neloccdo .ne. 0) then
                        iflag = 1
                        read (5, *, err=101, end=101) (occupationdo(i), i=1, neloccdo)
                        iflag = 0
101                     if (iflag .eq. 1) then
                            call error(' read_occupations ', ' Error reading spin down occupations ', -1, rank)
                            iflagall = iflagall + 1
                        end if
                    end if
                end if
            else
                ! complex wavefunction
                if (.not. symmagp) then
                    if (neloccdo .ne. 0) then
                        iflag = 1
                        read (5, *, err=102, end=102) (occupationdo(i), i=1, neloccdo)
                        iflag = 0
102                     if (iflag .eq. 1) then
                            call error(' read_occupations ', ' Provide spin down occupations for &
                                    &non-symmetric complex WF!', -1, rank)
                            iflagall = iflagall + 1
                        end if
                    end if
                end if
            end if
        end if
        if (rank .eq. 0) write (6, *)
        call checkiflagerr(iflagall, rank, " ERROR reading occupations! ")
        !
        ! in the case of a symmetric complex wavefunction, also occupations
        ! for down-spin electrons must be defined since it behaves
        ! like if it was a LSDA calculation with non symmetric AGP.
        !
        if (ipc .eq. 2 .and. symmagp) then
            neloccdo = nelocc
            do i = 1, nelocc
                occupations(i) = occupations(i)/2.d0
                occupationdo(i) = occupations(i)
            end do
        end if

        if (nelocc .gt. nmoltot) then
            if (rank .eq. 0) then
                write (6, *)
                write (6, '(a,I6)') ' INCREASE mol. orb. in fort.10 by ', nelocc - nmoltot
            end if
            call error(' read_occupations ', ' Too few molecular orbitals < nelocc ', 1, rank)
        end if
        if (nelocc .gt. bands) then
            if (rank .eq. 0) then
                write (6, *)
                write (6, '(a,I6)') ' INCREASE bands in input by ', nelocc - bands
            end if
            call error(' read_occupations ', ' Too few molecular orbitals < nelocc ', 1, rank)
        end if
#ifdef PARALLEL
        call mpi_bcast(occupations, bands, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (yeslsda .or. ipc .eq. 2) &
            call mpi_bcast(occupationdo, bands, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

        if (rank .eq. 0 .and. .not. compute_bands) then
            write (6, *) ' Total occupation UP read from std input'
            !do i=max(neldo-10,1),bands
            do i = 1, bands
                write (6, *) i, occupations(i)
            end do
            if (yeslsda .or. ipc .eq. 2) then
                write (6, *) ' Total occupation DOWN read from std input'
                !   do i=max(neldo-10,1),bands
                do i = 1, bands
                    write (6, *) i, occupationdo(i)
                end do
            end if
        end if
        !
        ! find last occupied orbital (lastpaired), the total number
        ! of occupied orbitals (countall) and set the number of
        ! fully (fractionally) occupied orbitals in count2 (costall).
        !
        lastpaired = 0
        count2 = 0
        countall = 0
        costall = 0.d0
        do i = 1, nelocc
            if (occupations(i) .ge. 0.d0) then
                lastpaired = i
                if (.not. yeslsda .and. ipc .eq. 1) then
                    if (occupations(i) .eq. 2.d0 .or. occupations(i) .eq. 0.d0) count2 = count2 + 1
                    if (occupations(i) .lt. 2.d0 .and. occupations(i) .ge. 0.d0) costall = costall + occupations(i)
                elseif (yeslsda .or. ipc .eq. 2) then
                    if (occupations(i) .eq. 1.d0 .or. occupations(i) .eq. 0.d0) count2 = count2 + 1
                    if (occupations(i) .lt. 1.d0 .and. occupations(i) .ge. 0.d0) costall = costall + occupations(i)
                end if
                countall = countall + 1
            end if
            occupations(i) = abs(occupations(i))
        end do
        if (rank .eq. 0) write (6, *) ' Read last occupied paired orbital =', lastpaired
        ! calculation costall
        if (countall .gt. count2) then
            costall = nint(costall)/dble(countall - count2)
            if (rank .eq. 0) write (6, *) ' Warning costall found =', costall
        else
            costall = 0.d0
        end if
        if (optocc .eq. 1) then
            countall = bands
            nelocc = bands
            neloccdo = bands
        end if

    end if ! endif occread

    return

end subroutine read_occupations
