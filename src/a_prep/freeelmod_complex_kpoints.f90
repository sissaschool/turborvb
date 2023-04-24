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

module freeelmod_complex

    use setup
    use parallel_module
    use compute_efermi, only: efermi, number_particles, number_particles_do
    use fourier_module
    use kpoints_mod, only: sum_kpoints_scalar_real8

    implicit none

    public
    logical :: yest
    integer :: nelorbc_in, nelorbjc_in, niesd_in, iesdrr_in, nelcolc_in, nshellmol, iesupr_cmol, occ_cmol, &
               iesup_cmol, iesupindmol, nmoltot, i_max, nnozeroc_in, nnozerojc_in, iesswrc_in, indnn_in, &
               inddiff, inds, iesfreerc_in, contraction_in, contractionj_in, count2, countall, count2n, &
               countalln, iter, lastpaired, dimjac, dimjac2, indproc
    integer, allocatable, dimension(:) :: ind_unpair, nozeroc_in, jbradetc_in, jbrajc_in, &
                                          nozerojc_in, jbradetnc_in, jbrajnc_in
    real(8) :: costall, edft_av, ecorr_av, exchange_av, spingrid_av
    real(8), dimension(:), allocatable :: dft_energies_all, corr_energy_all, exchange_energy_all, spin_grid_all

    logical, private :: check_kp
    integer, private :: i, j, k, ii, jj, kk, proc, indkp, is, jjs, ks, ind
    real(8), private :: sumh, sumo, sumhup, sumhdo, sum_eigv, sum_eigvdo, sumden, sumspin, &
                        sumhodo, sumodo, epsder_dft, lworkr, derivR, derivI, max_errdft, cycle_timep

    complex(8), private :: deriv
    real(8), dimension(:, :), allocatable :: xkp_sav

    real(8), external :: cclock

#ifdef __SCALAPACK
    integer, private :: ic, ir
#endif

contains

    subroutine self_consistent_run

        implicit none

        real*8, external :: ddot, drand1
        complex*16, external :: zdotu, zdotc_
        character(len=30) :: dum1, dum2, dum3
        logical :: is_same_kp(3)
#ifdef PARALLEL
        include "mpif.h"
#endif
        !
        ! read occupations from file occupationlevels.dat
        ! We have a single file for all k-points in the calculation.
        !
        if (iopt .eq. 1) then

            molecorb = 0.d0
            if (yeslsda .or. ipc .eq. 2) molecorbdo = 0.d0
            dent = dens0
            if (yeslsda) spint = spin0

        else

            if (occopen .and. .not. occread) then
                if (manyfort10) then
                    allocate (xkp_sav(3, nk)) ! for check that k-points are correct
                    xkp_sav = 0.d0
                end if
                if (rank .eq. 0) then
                    occupations = 0.d0
                    if (yeslsda .or. ipc .eq. 2) occupationdo = 0.d0
                    rewind (11)
                    read (11, *)
                    read (11, *) bandso, nproco, meshproco, nko
                    read (11, *)
                    do kk = 1, nko
                        if (manyfort10) read (11, *) dum1, xkp_sav(:, kk)
                        read (11, *)
                        if (bands .gt. bandso) then
                            do ii = 1, bands
                                read (11, *) i, eigmol_sav(ii, kk), occupations_sav(ii, kk)
                            end do
                            eigmol_sav(bandso + 1:bands, kk) = eigmol_sav(bandso, kk) + 100.d0
                        else
                            do ii = 1, bands
                                read (11, *) i, eigmol_sav(ii, kk), occupations_sav(ii, kk)
                            end do
                            do ii = bands + 1, bandso
                                read (11, *) i, psip(ii), psip(ii + 1)
                            end do
                        end if
                        if (yeslsda .or. ipc .eq. 2) then
                            read (11, *)
                            if (bands .gt. bandso) then
                                do ii = 1, bands
                                    read (11, *) i, eigmoldo_sav(ii, kk), occupationsdo_sav(ii, kk)
                                end do
                                eigmoldo_sav(bandso + 1:bands, kk) = eigmoldo_sav(bandso, kk) + 100.d0
                            else
                                do ii = 1, bands
                                    read (11, *) i, eigmoldo_sav(ii, kk), occupationsdo_sav(ii, kk)
                                end do
                                do ii = bands + 1, bandso
                                    read (11, *) i, psip(ii), psip(ii + 1)
                                end do
                            end if
                        end if
                        read (11, *) ! blank after one k-point
                    end do
                    write (6, *) ' Occupations and eigenvalues correctly read from file '
                end if ! endif rank.eq.0
#ifdef PARALLEL
                call mpi_bcast(nko, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(bandso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(nproco, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(meshproco, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                ! broadcast total eigenvalues and occupations to every processor
                if (manyfort10) call mpi_bcast(xkp_sav, size(xkp_sav), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(eigmol_sav, size(eigmol_sav), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(occupations_sav, size(occupations_sav), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                if (yeslsda .or. ipc .eq. 2) then
                    call mpi_bcast(occupationsdo_sav, size(occupationsdo_sav), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                    call mpi_bcast(eigmoldo_sav, size(eigmoldo_sav), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                end if
#endif
                !
                ! now distribute to the pools occupations and k-points
                !
                if (manyfort10) then
                    ! check if k-points read are correct for each processor, otherwise stop.
                    do ii = 1, nk
                        is_same_kp(:) = (abs(xkp_sav(:, ii) - xkp(:, ii)) < 1d-7)
                        if (all(is_same_kp)) then ! only if the right k-point
                            occupations(1:bands) = occupations_sav(1:bands, ii)
                            eigmol(1:bands) = eigmol_sav(1:bands, ii)
                            if (yeslsda .or. ipc .eq. 2) then
                                occupationdo(1:bands) = occupationsdo_sav(1:bands, ii)
                                eigmoldo(1:bands) = eigmoldo_sav(1:bands, ii)
                            end if
                            check_kp = .true.
                        else

                            if (rank .eq. 0) write (6, *) ' read, computed K =', (xkp_sav(jj, ii), xkp(jj, ii), jj=1, 3)

                            call error(' freeelmod_complex ', ' one or more k-points read from &
                                    &file do not coincide with input fort.10 ones.Check your input !! ', 1, rank)

                        end if
                    end do
                    deallocate (xkp_sav)
                end if

            else ! in this case one can read occupations from input file
                ! and file occupations.dat is not empty.
                ! occupations read form standard input override the ones from occupations.dat
                if (rank .eq. 0 .and. occopen) then
                    rewind (11)
                    read (11, *)
                    read (11, *) bandso, nproco, meshproco, nko
                elseif (rank .eq. 0) then
                    bandso = 0
                    nproco = 0
                    meshproco = 0
                    nko = 0
                end if
#ifdef PARALLEL
                call mpi_bcast(bandso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(nproco, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(meshproco, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(nko, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
                if (occread) then
                    oldocc = occupations
                    if (yeslsda .or. ipc .eq. 2) oldoccdo = occupationdo
                end if
            end if ! endif occopen.and..not.occread

            if (.not. occread) then
                oldocc = occupations
                if (yeslsda .or. ipc .eq. 2) oldoccdo = occupationdo
            end if

        end if ! endif if(iopt.eq.1)

        if ((((nproco .ne. nproc .or. meshproc .ne. meshproco) .and. .not. noread) .or. iopt .eq. 2)) then
            if (nproc .ne. nproco) then
                call error(' self_consistent_run ', &
                           ' Impossible to continue fast with different # processors ', -1, rank)
                noread = .true. ! compute again matrix elements
                yesreadf = .false.
            end if
        end if
        !
        ! initialize vhartreeq in reciprocal space
        ! Ewald correction to Hartree potential.

        !
        call evalvhartreeq()
        !

        ! initialize overlap and hamiltonian matrices
        ! once for all in a SC run. In a nonSC run they need
        ! to be initialized for each k-points.
        call initialize_mats_new()
        if (rank .eq. 0) then
            if (iespbc) then
                write (6, '(a,2E21.8)') ' Initial Hartree potential: ', sum(vhartree(1:meshproc)), sum(vhartreeq(1:meshproc))
                write (6, '(a,F15.8)') ' Kappa found:         ', kappa
                write (6, '(a,F15.8)') ' Ewald constant found:', ewaldion1b + (1.d0 - weightvh)*ewaldel1b
                write (6, '(a,F15.8)') ' eselfion =           ', ewaldion1b
                write (6, '(a,F15.8)') ' eself1bel =          ', ewaldel1b
            end if
        end if

        ! Recompute the charge/spin density with the MOs read from fort.10
        if (writescratch .eq. 0 .and. iopt .eq. 0) then
            if (rank .eq. 0) write (6, *) 'Warning reading initial density from file'
            rewind (unit_scratch_distributed)
            read (unit_scratch_distributed)
            read (unit_scratch_distributed)
            if (ipc .eq. 2) read (unit_scratch_distributed)
            if (double_mesh) read (unit_scratch_distributed)
            if (yeslsda) then
                read (unit_scratch_distributed) dent, spint
            else
                read (unit_scratch_distributed) dent
            end if
        else
            dent(:) = 0.d0
            if (yeslsda) spint(:) = 0.d0
            call updenorb_new
        end if

        if (typeopt .eq. 2) then
            dent_before(:, 1) = dens0
            dent_after(:, 1) = dent(:)
            dent_before(:, 2) = dent(:)
            if (yeslsda) then
                spint_before(:, 1) = spin0
                spint_after(:, 1) = spint(:)
                spint_before(:, 2) = spint(:)
            end if
        end if

        if (h_charge .ne. 0.d0) then
            indmesh = 0
            ind = 0
            !
            ! initialize the same spin grid for each k-point
            indproc = 0
            !
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        if (indproc .eq. rankrep) then
                            indmesh = indmesh + 1
                            ind = ind + 1
                            ks = dble(k - 1)/(dble(nz)/dble(nzs)) + 1
                            jjs = dble(j - 1)/(dble(ny)/dble(nys)) + 1
                            is = dble(i - 1)/(dble(nx)/dble(nxs)) + 1
                            if (charge_input(is, jjs, ks) .gt. 0) then
                                gridcharge(indmesh) = .true.
                                gridnocharge(indmesh) = .false.
                            elseif (charge_input(is, jjs, ks) .lt. 0) then
                                gridcharge(indmesh) = .false.
                                gridnocharge(indmesh) = .false.
                            else
                                gridnocharge(indmesh) = .true.
                                gridcharge(indmesh) = .false.
                            end if
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                    end do
                end do
            end do
        end if

        if (yeslsda .and. (iopt .eq. 1 .or. randspin .lt. 0.d0 .or. h_field .ne. 0)) then

            indmesh = 0
            ind = 0
            !
            ! initialize the same spin grid for each k-point
            !
            indproc = 0
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        if (indproc .eq. rankrep) then
                            indmesh = indmesh + 1
                            ind = ind + 1
                            ks = dble(k - 1)/(dble(nz)/dble(nzs)) + 1
                            jjs = dble(j - 1)/(dble(ny)/dble(nys)) + 1
                            is = dble(i - 1)/(dble(nx)/dble(nxs)) + 1
                            if (spin_input(is, jjs, ks) .gt. 0) then
                                gridspin(indmesh) = .true.
                                gridnospin(indmesh) = .false.
                            elseif (spin_input(is, jjs, ks) .lt. 0) then
                                gridspin(indmesh) = .false.
                                gridnospin(indmesh) = .false.
                            else
                                gridnospin(indmesh) = .true.
                                gridspin(indmesh) = .false.
                            end if
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                    end do
                end do
            end do
        end if

        if (iopt .eq. 1) then ! this is just the case of SC-run, since for NonSC run
            ! iopt is automatically set to 0.

            if (rank .eq. 0) write (6, *) ' Eigenvalue one body Ham  '
            !
            lworkr = 0
            !
            ! evaluate sum of hamiltonian and overlap matrices
            !
            call checksum
            !
            ! diagonalize KS Hamiltonian in non-orthogonal basis set
            ! Cavazzoni-Sorella diagonalization algorithm
            !
            call diagonalize_hamiltonian(1, lworkr)
            !
            ! KS orbitals occupation with or without smearing
            !
            call upocc_kpoints(optocc, epsshell, 1)
            !
            ! charge and spin density on the mesh
            !
            call updenorb_new

            if (yeslsda .and. randspin .lt. 0.d0) then
                if (rank .eq. 0) write (6, '(A)') ' Initialization spin magnetization '
                indmesh = 0
                indproc = 0
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            if (indproc .eq. rankrep) then
                                indmesh = indmesh + 1
                                if (gridnospin(indmesh)) then
                                    spint(indmesh) = 0.d0
                                else
                                    if (gridspin(indmesh)) then
                                        spint(indmesh) = dent(indmesh)
                                    else
                                        spint(indmesh) = -dent(indmesh)
                                    end if
                                end if
                            end if
                            indproc = indproc + 1
                            if (indproc .eq. nprocrep) indproc = 0
                        end do
                    end do
                end do
                if (h_field .eq. 0.d0 .and. randspin .ge. 0.d0 .and. allocated(gridspin)) &
                        & deallocate (gridspin, gridnospin)
            end if

            if (yeslsda) deallocate (spin_input)

            if (mod(typeopt, 2) .eq. 0) then

                ! total charge and spin and spin fluctuations
                call cutdens

                !if(rank.eq.0) write(6,*) ' Initial density ',denstot
                if (rank .eq. 0) write (6, *) ' Spin variance before =', spin2tot
                if (rank .eq. 0) write (6, *) ' |spin|  =', spingrid/2
                if (rank .eq. 0 .and. yeslsda) write (6, *) ' Initial polarization ', spintot
                if (abs(denstot - nel) .gt. 1d-4) then
                    if (rank .eq. 0) write (6, *) ' Warning  initial density wrong, rescaling'
                    dent(:) = dent(:)*dble(nel)/denstot
                end if
                if (yeslsda) then
                    if (rank .eq. 0) write (6, *) ' Initial spin polarization ', spintot
                    if (rank .eq. 0) write (6, *) ' Initial variance spin  ', spin2tot
                    if (abs(spintot - nelup + neldo) .gt. 1d-4) then
                        if (rank .eq. 0) write (6, *) ' Warning  initial polarization  wrong, shifting'
                        cost = sum(spint(:) - spin0)/meshproc
                        spint(:) = spint(:) - cost
                    end if
                end if
            end if ! endif typeopt=even
            !
            ! DFT total energy
            if (.not. occread) then
                occupations = oldocc
            else
                oldocc = occupations
            end if
            edft = sum(occupations(1:nelocc)*eigmol(1:nelocc))

            if (yeslsda .or. ipc .eq. 2) then
                if (.not. occread) then
                    occupationdo = oldoccdo
                else
                    oldoccdo = occupationdo
                end if
                edft = edft + sum(occupationdo(1:neloccdo)*eigmoldo(1:neloccdo))
            end if
            ! print out non interacting energies for each k-point in the file fort.21
            if (print_energies) then
                dft_energies_all = 0.d0
                call collect_kpoints(dft_energies_all, nk, edft, -1)
            end if
            ! sum energies over the k-points
            call sum_kpoints_scalar_real8(edft, commcolrep_mpi, -1)
            if (rank .eq. 0) write (6, *) ' Efree ', edft

            ! print non-interacting energies in the "fort.21" formatted file
            if (print_energies) then

                if (decoupled_run) then
                    edft_av = 0.d0
                    do i = 1, nk
                        edft_av = edft_av + dft_energies_all(i)
                    end do
                    edft_av = edft_av/nk
                else
                    edft_av = edft
                end if

                if (rank .eq. 0) then
                    write (unit_print_energies, '(a)') ' # kx   ky    kz    w(k)    E(k)  N_up(k)    N_do(k) '
                    write (unit_print_energies, '(a)') ' # Non interacting DFT energy '
                    write (unit_print_energies, '(a,F15.8)') ' # Efree averaged: ', edft_av
                    if (yeslsda .or. ipc .eq. 2) then
                        do ii = 1, nk
                            write (unit_print_energies, '(4F8.4,F15.8,2I6)') xkp(:, ii), wkp(ii), &
                                dft_energies_all(ii), number_particles(ii), number_particles_do(ii)
                        end do
                    else
                        do ii = 1, nk
                            write (unit_print_energies, '(4F8.4,F15.8,I6)') xkp(:, ii), wkp(ii), &
                                dft_energies_all(ii), number_particles(ii)
                        end do
                    end if
                    write (unit_print_energies, *)
                end if

            end if
            !
            ! End initialization part iopt=1
            ! evaluation of the hamiltonian given the read/computed density
            !
        elseif (iopt .ne. 1) then ! read from input iopt.ne.1

            if (rank .eq. 0) write (6, *) ' Reading orbitals from input '
            if (rank .eq. 0) write (6, *) ' Eigenvalue DFT  Ham  '

            if (optocc .eq. 1) then
                nelocc = bands
                neloccdo = bands
                countall = bands
            end if

            if (meshproc .ne. meshproco) then
                if (rank .eq. 0) write (6, *) ' Warning orthogonalizing orbitals'
#ifdef __SCALAPACK
                if (ipc .eq. 1) then
                    call graham_scalapack(molecorb, oversl, psip, nelorbu, nelorbu, &
                                          nelocc, info, descla, desch, size(oversl, 1), rank)
                    if (yeslsda) then
                        call graham_scalapack(molecorbdo, oversl, psip, nelorbu, nelorbu, &
                                              neloccdo, info, descla, desch, size(oversl, 1), rank)
                    end if
                else
                    call graham_scalapack_complex(molecorb, oversl, psip, nelorbu, nelorbu, &
                                                  nelocc, info, descla, desch, size(oversl, 1)/2, rankrep)
                    call graham_scalapack_complex(molecorbdo, oversldo, psip, nelorbu, nelorbu, &
                                                  neloccdo, info, descla, desch, size(oversldo, 1)/2, rankrep)
                end if
#else
                if (ipc .eq. 1) then
                    call graham(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                    if (yeslsda) &
                        call graham(molecorbdo, overs, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
                else
                    call graham_complex(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                    call graham_complex(molecorbdo, oversdo, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
                end if

#endif

#ifdef PARALLEL
                call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier
                !
                ! check errors/warning from orthogonalization
                !
                if (info .ne. 0) call error(' freeelmod_complex ', &
                                            ' info different from 0 in orthogonalization ', -1, rank)
                if (.not. occopen .and. .not. occread) &
                    call error(' freeelmod_complex ', ' Occupations not defined ', 1, rank)
                !
                ! charge and spin density on the mesh
                call updenorb_new
                ! total charge and spin density
                call cutdens

            end if

            call uphamilt_new
            call checksum

            if (rank .eq. 0) write (6, *) 'Input matrix H before = ', sumhup, sumhdo, sumh, sumo
            if (rank .eq. 0) write (6, *) 'Input dens^2 spin^2   = ', sumden, sumspin
            if (rank .eq. 0) write (6, *) ' Initial tovpot        = ', totvpot

            lworkr = 0
            !
            !  In order to evaluate lworkr: the optimal scratch size for dsyevx
            !
#ifdef __SCALAPACK
            if (ipc .eq. 1) then
                call eval_hamilt(nelorbu, oversl, hamiltl, molecorbl&
                     &, eigmol, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 1, premat, eigmat, info, bands&
                     &, size(molecorbl, 1), mincond)
                if (yeslsda) then
                    call eval_hamilt(nelorbu, oversl, hamiltldo, molecorbldo&
                         &, eigmoldo, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 0, premat, eigmat, info, bands&
                         &, size(molecorbldo, 1), mincond)
                end if
            else
                call eval_hamilt_complex(nelorbu, oversl, hamiltl, molecorbl, umatl&
                     &, eigmol, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 1, premat, eigmat, info, bands&
                     &, size(molecorbl, 1)/2, mincond)
                call eval_hamilt_complex(nelorbu, oversldo, hamiltldo, molecorbldo, umatldo&
                     &, eigmoldo, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 1, premat, eigmat_down, info, bands&
                     &, size(molecorbldo, 1)/2, mincond)
            end if

#else
            if (ipc .eq. 1) then
                call eval_hamilt(nelorbu, overs, hamilt, molecorb_old&
                        &, eigmol, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 1, &
                        premat, eigmat, info, bands, 1, mincond)
                if (yeslsda) &
                        call eval_hamilt(nelorbu, overs, hamiltdo, molecorbdo_old&
                                &, eigmoldo, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 0, &
                                premat, eigmat, info, bands, 1, mincond)
            else
                call eval_hamilt_complex(nelorbu, overs, hamilt, molecorb_old, umatl&
                        &, eigmol, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 1, premat, eigmat, info&
                        &, bands, 1, mincond)
                call eval_hamilt_complex(nelorbu, oversdo, hamiltdo, molecorbdo_old, umatldo&
                        &, eigmoldo, nelorbu, epsover, eps_mach, rankrep, 1, lworkr, 1, premat, eigmat_down, &
                        info, bands, 1, mincond)
            end if
#endif

            call upocc_kpoints(optocc, epsshell, 0)

            if (typeopt .eq. 2) then

                dent_before(:, 1) = dens0
                dent_after(:, 1) = dent(:)
                dent_before(:, 2) = dent(:)

                if (yeslsda) then
                    spint_before(:, 1) = spin0
                    spint_after(:, 1) = spint(:)
                    spint_before(:, 2) = spint(:)
                end if

            end if

            if (mod(typeopt, 2) .eq. 0) then

                if (.not. yesreadf) then

                    !            if(rank.eq.0) write(6,*) 'Warning reading initial density from file'
                    !            rewind(unit_scratch_distributed)
                    !            read(unit_scratch_distributed)
                    !            if(ipc.eq.2) read(unit_scratch_distributed)
                    !            if(yeslsda) then
                    !               read(unit_scratch_distributed) dent,spint
                    !            else
                    !               read(unit_scratch_distributed) dent
                    !            endif
                    !         else
                    call updenorb_new
                end if

                call cutdens

                if (abs(denstot - nel) .gt. 1d-4) then
                    if (rank .eq. 0) write (6, *) ' Warning  initial density wrong, rescaling'
                    dent(:) = dent(:)*dble(nel)/denstot
                end if

                if (yeslsda .and. randspin .lt. 0) then
                    if (rank .eq. 0) write (6, *) ' Initialization spin magnetization  '
                    !       indr=-1
                    indmesh = 0
                    indproc = 0
                    do k = 1, nz
                        do j = 1, ny
                            do i = 1, nx
                                if (indproc .eq. rankrep) then
                                    !           indr=indr+1
                                    !           if(indr-(indr/nproc)*nproc.eq.rank) then
                                    indmesh = indmesh + 1
                                    if (gridnospin(indmesh)) then
                                        spint(indmesh) = 0.d0
                                    else
                                        if (gridspin(indmesh)) then
                                            spint(indmesh) = dent(indmesh)
                                        else
                                            spint(indmesh) = -dent(indmesh)
                                        end if
                                    end if
                                    !           endif
                                end if
                                indproc = indproc + 1
                                if (indproc .eq. nprocrep) indproc = 0
                            end do
                        end do
                    end do

                    call cutdens
                    if (h_field .eq. 0.d0 .and. allocated(gridspin)) deallocate (gridspin, gridnospin)
                end if
                if (yeslsda) deallocate (spin_input)

                if (yeslsda) then
                    if (rank .eq. 0) write (6, *) ' Initial spin polarization ', spintot
                    if (rank .eq. 0) write (6, *) ' Initial  |spin |  = ', spingrid/2
                    if (rank .eq. 0) write (6, *) ' Initial variance spin  ', spin2tot
                    if (abs(spintot - nelup + neldo) .gt. 1d-4) then
                        if (rank .eq. 0) write (6, *) ' Warning  initial polarization  wrong, shifting'
                        cost = sum(spint(:) - spin0)/meshproc
                        spint(:) = spint(:) - cost
                    end if
                end if

            end if ! endif typeopt=4/0

            edft = sum(occupations(1:bands)*eigmol(1:bands))
            if (yeslsda .or. ipc .eq. 2) then
                edft = edft + sum(occupationdo(1:bands)*eigmoldo(1:bands))
            end if

            ! print out non interacting energies for each k-point in the file fort.21
            if (print_energies) then
                dft_energies_all = 0.d0
                call collect_kpoints(dft_energies_all, nk, edft, -1)
                if (rank .eq. 0) then
                    if (yeslsda .or. ipc .eq. 2) then
                        write (unit_print_energies, '(a)') ' # kx   ky    kz    w(k)    E(k)    N_up(k)    N_do(k)'
                        write (unit_print_energies, '(a)') ' # Initial total DFT energy '
                        do ii = 1, nk
                            write (unit_print_energies, '(4F8.4,F15.8,2I6)') xkp(:, ii), wkp(ii), &
                                dft_energies_all(ii), number_particles(ii), number_particles_do(ii)
                        end do
                    else
                        write (unit_print_energies, '(a)') ' # kx   ky    kz    w(k)    E(k)    N_up(k) '
                        write (unit_print_energies, '(a)') ' # Initial total DFT energy '
                        do ii = 1, nk
                            write (unit_print_energies, '(4F8.4,F15.8,I6)') xkp(:, ii), wkp(ii), &
                                dft_energies_all(ii), number_particles(ii)
                        end do
                    end if
                    write (unit_print_energies, *)
                end if
            end if

            call sum_kpoints_scalar_real8(edft, commcolrep_mpi, -1)
            edft = edft + totvpot

            if (rank .eq. 0) write (6, *) ' Initial dft energy =', edft

            if (occopen) then
                occupations = oldocc
                if (yeslsda .or. ipc .eq. 2) occupationdo = oldoccdo
            end if

            edft = 0.d0

        end if ! endif iopt.ne.0

        if (rank .eq. 0) write (6, '(a)') ' DFT initialization OK '

        ! add a random component to the spin
        if (yeslsda .and. randspin .gt. 0.d0) then

            if (rank .eq. 0) write (6, *) ' Warning adding a random component to the spin '
            if (mod(typeopt, 2) .eq. 0) then
                cost = 0.d0
                do ii = 1, meshproc
                    cost1 = randspin*(0.5d0 - drand1())*dble(nel)/mesh
                    if (abs(cost1 + spint(ii)) .ge. dent(ii)) cost1 = 0.d0 ! rejected
                    spint(ii) = spint(ii) + cost1
                    cost = cost + cost1
                end do
                cost = cost/meshproc
                spint = spint - cost
            else
                do jj = 1, bands
                    do ii = 1, nelorbu
                        ! here I add the random component to both real and imaginary part of MOs
                        if (ipc .eq. 1) then
                            molecorb(ii, jj) = molecorb(ii, jj) + randspin/dsqrt(dble(mesh))*(drand1() - 0.5d0)
                            molecorbdo(ii, jj) = molecorbdo(ii, jj) + randspin/dsqrt(dble(mesh))*(drand1() - 0.5d0)
                        else
                            molecorb(2*ii - 1:2*ii, jj) = molecorb(2*ii - 1:2*ii, jj) &
                                                          + zone*randspin/dsqrt(dble(mesh))*(drand1() - 0.5d0)
                            molecorbdo(2*ii - 1:2*ii, jj) = molecorbdo(2*ii - 1:2*ii, jj) &
                                                            + zone*randspin/dsqrt(dble(mesh))*(drand1() - 0.5d0)
                        end if
                    end do
                end do
            end if
        end if

        ! orthogonalize orbitals
        if (orthodiag) call improvediag

        if (mod(typeopt, 2) .eq. 1) then
#ifdef __SCALAPACK
            if (ipc .eq. 1) then
                call graham_scalapack(molecorb, oversl, psip, nelorbu, nelorbu, nelocc &
                                      , info, descla, desch, size(oversl, 1), rank)
                if (yeslsda) &
                    call graham_scalapack(molecorbdo, oversl, psip, nelorbu, nelorbu, neloccdo &
                                          , info, descla, desch, size(oversl, 1), rank)
            else
                call graham_scalapack_complex(molecorb, oversl, psip, nelorbu, nelorbu, nelocc &
                                              , info, descla, desch, size(oversl, 1)/2, rank)
                call graham_scalapack_complex(molecorbdo, oversldo, psip, nelorbu, nelorbu, neloccdo &
                                              , info, descla, desch, size(oversl, 1)/2, rank)
            end if
#else
            if (ipc .eq. 1) then
                call graham(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                if (yeslsda) &
                    call graham(molecorbdo, overs, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
            else
                call graham_complex(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                call graham_complex(molecorbdo, oversdo, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
            end if
#endif

#ifdef PARALLEL
            call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier

            molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
            if (yeslsda .or. ipc .eq. 2) molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)
        end if

        ! For whatever case molecorb_old=molecorb
        molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
        if (yeslsda .or. ipc .eq. 2) molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)

#ifdef PARALLEL
        ! all processor in the pool must have consistent w.f.
        call mpi_bcast(molecorb, size(molecorb), MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
        if (yeslsda .or. ipc .eq. 2) &
            call mpi_bcast(molecorbdo, size(molecorbdo), MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
#endif

        ! ---------------------------- Initialize DFT self-consistent cycle ---------------------- !

        ! write sum of the elements of initial eigenvectors for checking
        sum_eigv = sum(molecorb(:, 1:nelocc))
        if (yeslsda .or. ipc .eq. 2) sum_eigvdo = sum(molecorbdo(:, 1:nelocc))
#ifdef PARALLEL
        if (manyfort10) then
            call reduce_base_real(1, sum_eigv, commcolrep_mpi, -1)
            sum_eigv = sum_eigv/nk
            if (yeslsda .or. ipc .eq. 2) then
                call reduce_base_real(1, sum_eigvdo, commcolrep_mpi, -1)
                sum_eigvdo = sum_eigvdo/nk
            end if
        end if
#endif
        if (rank .eq. 0) write (6, '(a,F20.6)') ' # Molecorb start:      ', sum_eigv
        if (rank .eq. 0 .and. (yeslsda .or. ipc .eq. 2)) &
            write (6, '(a,F20.6)') ' # Molecorb start down: ', sum_eigvdo

        if (iopt .eq. 1) then
            dentold = dens0
            if (yeslsda) spintold = spin0
        else
            dentold = dent
            if (yeslsda) spintold = spint
        end if
        if (typeopt .eq. 2) then
            dentnew = dentold
            if (yeslsda) spintnew = spintold
        end if

        mincost = 1.d0
        mixingstep = mixing
        if (typeopt .eq. 3) mixingstep = 0.d0

        loading_time = 0.d0
        zgemm_time = 0.d0
        dgemm_time = 0.d0
        diag_time = 0.d0
        symtime = 0.d0
        cycle_time = 0.d0
        time_fft = 0.d0

        iter = 0
        itersto = 0
        errdft = 2.d0*epsdft

        ! ----------for test------------
        if (typeopt .lt. 0) then

            if (rank .eq. 0) write (6, *) ' Check DFT  functional '
            typeopt = -typeopt
            ! choose a random direction
            epsder_dft = mixing
            mixingstep = 0.d0

            call cyclesteep_complex

            edftp = edft

            molecorbs(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
            if (yeslsda .or. ipc .eq. 2) molecorbdos(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)

            do i = 1, bands
                do j = 1, nelorbu
                    if (ipc .eq. 1) then
                        molecorb(j, i) = molecorb(j, i) + 2.d0*epsder_dft*dcos(dble(i + j**2))
                        if (yeslsda) molecorbdo(j, i) = molecorbdo(j, i) + 2.d0*epsder_dft*dcos(dble(i**2 + j**3))
                    else
                        molecorb(2*j - 1, i) = molecorb(2*j - 1, i) + 2.d0*epsder_dft*dcos(dble(i + j**2))
                        molecorbdo(2*j - 1, i) = molecorbdo(2*j - 1, i) + 2.d0*epsder_dft*dcos(dble(i**2 + j**3))
                        molecorb(2*j, i) = molecorb(2*j, i) + 2.d0*epsder_dft*dcos(dble(i + j**2))
                        molecorbdo(2*j, i) = molecorbdo(2*j, i) + 2.d0*epsder_dft*dcos(dble(i**2 + j**3))
                    end if
                end do
            end do

            ! Orthogonalize before (more accurate)
#ifdef __SCALAPACK
            if (typeopt .eq. 1) then
                if (ipc .eq. 1) then
                    call graham_scalapack(molecorb, oversl, psip, nelorbu, nelorbu, nelocc &
                                          , info, descla, desch, size(oversl, 1), rank)
                    if (yeslsda) &
                        call graham_scalapack(molecorbdo, oversl, psip, nelorbu, nelorbu, neloccdo &
                                              , info, descla, desch, size(oversl, 1), rank)
                else
                    call graham_scalapack_complex(molecorb, oversl, psip, nelorbu, nelorbu, nelocc, &
                                                  info, descla, desch, size(oversl, 1)/2, rank)
                    call graham_scalapack_complex(molecorbdo, oversldo, psip, nelorbu, nelorbu, neloccdo, &
                                                  info, descla, desch, size(oversl, 1)/2, rank)
                end if
            end if
#else
            if (typeopt .eq. 1) then
                if (ipc .eq. 1) then
                    call graham(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                    if (yeslsda) &
                        call graham(molecorbdo, overs, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
                else
                    call graham_complex(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                    call graham_complex(molecorbdo, oversdo, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
                end if
            end if
#endif

#ifdef PARALLEL
            call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier

            molecorbs(1:ipc*nelorbu, 1:bands) = molecorbs(1:ipc*nelorbu, 1:bands) &
                                                - molecorb(1:ipc*nelorbu, 1:bands)
            if (yeslsda .or. ipc .eq. 2) &
                molecorbdos(1:ipc*nelorbu, 1:bands) = molecorbdos(1:ipc*nelorbu, 1:bands) &
                                                      - molecorbdo(1:ipc*nelorbu, 1:bands)

            costder = -1d0/epsder_dft* &
                      sum(molecorbs(1:ipc*nelorbu, 1:nelocc)*molecorb_old(1:ipc*nelorbu, 1:nelocc))
            if (yeslsda .or. ipc .eq. 2) then
                costder = costder - 1.d0/epsder_dft* &
                          sum(molecorbdos(1:ipc*nelorbu, 1:neloccdo)*molecorbdo_old(1:ipc*nelorbu, 1:neloccdo))
            end if

            !        Orthogonalize after
#ifdef __SCALAPACK
            if (typeopt .eq. 3) then
                if (ipc .eq. 1) then
                    call graham_scalapack(molecorb, oversl, psip, nelorbu, nelorbu, nelocc&
                         &, info, descla, desch, size(oversl, 1), rank)
                    if (yeslsda) call graham_scalapack(molecorbdo, oversl, psip, nelorbu, nelorbu&
                         &, neloccdo, info, descla, desch, size(oversl, 1), rank)
                else
                    call graham_scalapack_complex(molecorb, oversl, psip, nelorbu, nelorbu, nelocc, &
                                                  info, descla, desch, size(oversl, 1)/2, rank)
                    call graham_scalapack_complex(molecorbdo, oversl, psip, nelorbu, nelorbu, &
                                                  neloccdo, info, descla, desch, size(oversl, 1)/2, rank)
                end if
            end if
#else
            if (typeopt .eq. 3) then
                if (ipc .eq. 1) then
                    call graham(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                    if (yeslsda) &
                        call graham(molecorbdo, overs, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
                else
                    call graham_complex(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                    call graham_complex(molecorbdo, oversdo, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
                end if
            end if
#endif

#ifdef PARALLEL
            call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier

            call cyclesteep_complex

            ! The more accurate symmetric expression for derivative is taken.

            if (rank .eq. 0) then
                if (yeslsda .or. ipc .eq. 2) then
                    if (ipc .eq. 1) then
                        deriv = costder - 1d0/epsder_dft*sum(molecorbs(1:nelorbu, 1:nelocc)*molecorb_old(1:nelorbu, 1:nelocc))&
                                & - 1d0/epsder_dft*sum(molecorbdos(1:nelorbu, 1:neloccdo)*molecorbdo_old(1:nelorbu, 1:neloccdo))
                        write (6, '(a,2F20.8)') ' Analitycal derivative DFT functional = ', deriv, costder
                    else
                        derivR = 0.d0
                        derivI = 0.d0
                        do ii = 1, nelorbu
                            ! real part
                            derivR = derivR + costder &
                                     - 1d0/epsder_dft*sum(molecorbs(2*ii - 1, 1:nelocc)*molecorb_old(2*ii - 1, 1:nelocc)) &
                                     - 1d0/epsder_dft*sum(molecorbdos(2*ii - 1, 1:neloccdo)*molecorbdo_old(2*ii - 1, 1:neloccdo))
                            ! imaginary part
                            derivI = derivI + costder - 1d0/epsder_dft* &
                                     sum(molecorbs(2*ii, 1:nelocc)*molecorb_old(2*ii, 1:nelocc)) &
                                     - 1d0/epsder_dft*sum(molecorbdos(2*ii, 1:neloccdo)*molecorbdo_old(2*ii, 1:neloccdo))
                        end do
                        deriv = derivR + zimg*derivI
                        write (6, '(a,3F20.8)') ' Analitycal derivative DFT functional = ', deriv, costder
                    end if
                else
                    write (6, '(a,2F20.8)') ' Analitycal derivative DFT functional = ', &
                            &costder - 1d0/epsder_dft*sum(molecorbs(1:nelorbu, 1:nelocc)* &
                            molecorb_old(1:nelorbu, 1:nelocc)), costder
                end if
                write (6, '(a,F20.8)') ' Numerical  derivative DFT functional = ', (edft - edftp)/epsder_dft, edft - edftp
            end if

#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop

        end if

#ifdef UNREL

#ifdef PARALLEL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier
#endif

        ! --------------------------------------------------
        !              main self-consistent cycle
        ! --------------------------------------------------

        if (rank .eq. 0) then
            write (6, *) ' '
            write (6, *) ' ------------------------------ '
            write (6, *) '        Starting SC cycle       '
            write (6, *) ' ------------------------------ '
            write (6, *) ' '
        end if

        cycle_timep = 0.d0
        cycle_timep = cclock()

        do while (errdft .gt. epsdft .and. iter .lt. maxit)

            iter = iter + 1

            if (mod(typeopt, 2) .eq. 0) then

                call cycle_complex ! typeopt = 2 -> linear mixing scheme
                ! typeopt = 4 -> Anderson mixing scheme

            elseif (typeopt .eq. 1) then

                call cyclesteep_complex ! typeopt = 1 -> steepest descent method

            end if

        end do

        cycle_time = cycle_time + (cclock() - cycle_timep)

        if (rank .eq. 0) then
            write (6, *) ' '
            write (6, *) ' ------------------------ '
            write (6, *) '     SC cycle completed   '
            write (6, *) ' ------------------------ '
            write (6, *) ' '
        end if

        if (mod(typeopt, 2) .eq. 1 .and. rank .eq. 0 .and. .not. manyfort10) then
            write (6, *) ' Estimated eigenvalues '
            do i = 1, bands
                if (occupations(i) .ne. 0.d0) then
                    write (6, *) i, eigmol(i)/occupations(i)
                end if
            end do
            if (yeslsda .or. ipc .eq. 2) then
                write (6, *) ' Estimated eigenvalues spin down  '
                do i = 1, bands
                    if (occupationdo(i) .ne. 0.d0) then
                        write (6, *) i, eigmoldo(i)/occupationdo(i)
                    end if
                end do
            end if
        end if

        if (mod(typeopt, 2) .eq. 1) then
            yespassed = .false.
        else
            yespassed = .true.
            edftp = edft
        end if

        if (mod(typeopt, 2) .eq. 1 .and. iter .lt. maxit) then ! only if converged
            yespassed = .true.
            edftp = edft
            mixing = 1.d0
            maxcg = 0
            mixingder = 0.d0
            optocc = 0
            optpar = .false. ! to avoid memory conflicts in parallel
            molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
            if (yeslsda .or. ipc .eq. 2) molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)
            if (.not. allocated(molecorbs)) then
                allocate (molecorbs(ipc*nelorbu, bands))
                if (yeslsda .or. ipc .eq. 2) allocate (molecorbdos(ipc*nelorbu, bands))
            end if

            molecorbs(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
            oldocc = occupations

            if (yeslsda .or. ipc .eq. 2) then
                molecorbdos(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)
                oldoccdo = occupationdo
            end if

            !recover matrix hamilt and compute eigenvalues eigenvectors
            call cycle_complex
        end if

        if (rank .eq. 0) then
            !
            ! write eigenvalues and occupations in the file
            ! "occupationlevels.dat".
            !
            rewind (11)
            write (11, *) ' bands / # of processors / mesh / # of k-points '
            write (11, *) bands, nproc, meshproc, nk ! writing total mesh dimension for NonSC calculations
            if (corr_hartree .and. scale_hartree .gt. 0.d0) then
                write (11, *) scale_hartreen, '<-- Final scale_hartree '
            else
                write (11, *)
            end if
            do jj = 1, nk
                if (manyfort10) write (11, 301) 'k-point:', xkp(1, jj), xkp(2, jj), xkp(3, jj)
                if (yeslsda .or. ipc .eq. 2) then
                    write (11, *) ' Eigenvalues up / Occupations up'
                    do ii = 1, bands
                        write (11, 300) ii, eigmol_sav(ii, jj), occupations_sav(ii, jj)
                    end do
                    write (11, '(a)') ' Eigenvalues down / Occupations down'
                    do ii = 1, bands
                        write (11, 300) ii, eigmoldo_sav(ii, jj), occupationsdo_sav(ii, jj)
                    end do
                else
                    write (11, '(a)') ' Eigenvalues / Occupations '
                    do ii = 1, bands
                        write (11, 300) ii, eigmol_sav(ii, jj), occupations_sav(ii, jj)
                    end do
                end if
                write (11, *)
            end do
            close (11)
300         format(3x, I3, X, 2f20.10)
301         format(2x, A8, X, 3f12.8)

        end if
        !
        ! sorting according to the occupation + a bit of energetic
        if (nelorbh .gt. 1 .and. (optocc .eq. 1 .or. occread) .and. yespassed) then
            ii = 0
            do i = 1, bands
                if (occupations(i) .ne. 0.d0) ii = ii + 1
                if (ii .le. count2) then ! no matter how the weight>0
                    ! the paired  orbitals come first
                    if (occupations(i) .ne. 0.d0) then
                        if (cost .ne. 2.d0) cost = 2.d0
                    else
                        cost = 0.d0
                    end if
                else
                    cost = occupations(i)
                    if (cost .ne. 0.d0) cost = 1.d0
                end if
                !      paired orbitals   cost=2
                !      others            cost=1

                psip(i) = cost - 0.5d0*(eigmol(i) - eigmol(1))/(eigmol(bands) - eigmol(1)) ! < 1
            end do

            call dsortx(psip, 1, bands, ipsip)

            if (yeslsda .or. ipc .eq. 2) then
                ii = 0
                do i = 1, bands
                    if (occupationdo(i) .ne. 0.d0) ii = ii + 1
                    if (ii .le. count2) then
                        if (occupationdo(i) .ne. 0.d0) then
                            cost = 2.d0
                        else
                            cost = 0.d0
                        end if
                    else
                        cost = 1.d0
                    end if
                    !      paired orbitals  cost=2
                    !      others           cost=1
                    psip(i + bands) = cost - 0.5d0*(eigmoldo(i) - eigmoldo(1))/(eigmoldo(bands) - eigmoldo(1)) ! < 1
                end do

                call dsortx(psip(bands + 1), 1, bands, ipsip(bands + 1))

            end if

            psip(1:bands) = eigmol(1:bands)
            do i = 1, bands
                eigmol(i) = psip(ipsip(bands - i + 1))
            end do
            if (rank .eq. 0 .and. .not. manyfort10) write (6, *) ' Sorted eigenvalues/occupations up '
            psip(1:bands) = occupations(1:bands)
            do i = 1, bands
                occupations(i) = psip(ipsip(bands - i + 1))
                if (rank .eq. 0 .and. .not. manyfort10) write (6, *) i, eigmol(i), occupations(i)
            end do
            if (yeslsda .or. ipc .eq. 2) then
                if (rank .eq. 0 .and. .not. manyfort10) then
                    write (6, *)
                    write (6, *) ' Sorted eigenvalues/occupations down '
                end if
                psip(bands + 1:2*bands) = eigmoldo(1:bands)
                do i = 1, bands
                    eigmoldo(i) = psip(bands + ipsip(2*bands - i + 1))
                end do
                psip(bands + 1:2*bands) = occupationdo(1:bands)
                do i = 1, bands
                    occupationdo(i) = psip(bands + ipsip(2*bands - i + 1))
                    if (rank .eq. 0 .and. .not. manyfort10) write (6, *) i, eigmoldo(i), occupationdo(i)
                end do
                if (rank .eq. 0 .and. .not. manyfort10) write (6, *)
            end if

            if (mod(typeopt, 2) .ne. 1) then

                molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
                do i = 1, bands
                    molecorb(1:ipc*nelorbu, i) = molecorb_old(1:ipc*nelorbu, ipsip(bands - i + 1))
                end do
                if (yeslsda .or. ipc .eq. 2) then
                    molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)
                    do i = 1, bands
                        molecorbdo(1:ipc*nelorbu, i) = molecorbdo_old(1:ipc*nelorbu, ipsip(2*bands - i + 1))
                    end do
                end if

            else
                ! do not introduce garbage in the occupied molecular orbitals
                occupations = oldocc
                molecorb(1:ipc*nelorbu, 1:bands) = molecorbs(1:ipc*nelorbu, 1:bands)
                molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
                if (yeslsda .or. ipc .eq. 2) then
                    occupationdo = oldoccdo
                    molecorbdo(1:ipc*nelorbu, 1:bands) = molecorbdos(1:ipc*nelorbu, 1:bands)
                    molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)
                end if
                !     fill also the unoccupied for other purposes
            end if
        end if ! endif yespassed=.true.

        if (mod(typeopt, 2) .eq. 0) then

            call updenorb_new
            molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
            if (yeslsda .or. ipc .eq. 2) molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)

            call uphamilt_new
            call evalevar
            !
            ! check orthogonality final eigenvectors
            !
            if (rankrep .eq. 0) then
                mincost = 0.d0
                do i = 1, nelocc
                    do j = i + 1, nelocc
                        if (ipc .eq. 1) then
                            cost = ddot(nelorbu, molecorb(1, i), 1, molecorb_old(1, j), 1)
                        else
                            cost = zdotc_(nelorbu, molecorb(1, i), 1, molecorb_old(1, j), 1)
                        end if
                        cost = dabs(cost)
                        if (cost .gt. mincost) mincost = cost
                    end do
                end do
                if (yeslsda .or. ipc .eq. 2) then
                    mincosts = 0.d0
                    do i = 1, neloccdo
                        do j = i + 1, neloccdo
                            if (ipc .eq. 1) then
                                cost = ddot(nelorbu, molecorbdo(1, i), 1, molecorbdo_old(1, j), 1)
                            else
                                cost = zdotc_(nelorbu, molecorbdo(1, i), 1, molecorbdo_old(1, j), 1)
                            end if
                            cost = dabs(cost)
                            if (cost .gt. mincosts) mincosts = cost
                        end do
                    end do
                end if
            end if

            if (rank .eq. 0) then
                write (6, *) ' Check orthogonality h psi with psi i =/ j ', mincost
                if (yeslsda .or. ipc .eq. 2) &
                    write (6, *) ' Check orthogonality h psi with psi i down =/ j ', mincosts
            end if
            !
            ! end orthogonality check
            !
            if (rank .eq. 0) write (6, *) ' Variational energy without orth. ', edftvar
            if (rank .eq. 0) write (6, *) ' Variational const with no orth. ', totvpot

            if (memlarge .and. orthodiag) then
#ifdef __SCALAPACK
                if (ipc .eq. 1) then
                    call graham_scalapack(molecorb_old, oversl, psip, nelorb, nelorbu, &
                                          nelocc, info, descla, desch, size(oversl, 1), rankrep)
                    if (yeslsda) &
                        call graham_scalapack(molecorbdo_old, oversl, psip, nelorb, nelorbu, &
                                              neloccdo, info, descla, desch, size(oversl, 1), rankrep)
                else
                    call graham_scalapack_complex(molecorb_old, oversl, psip, nelorb, nelorbu, &
                                                  nelocc, info, descla, desch, size(oversl, 1)/2, rankrep)
                    call graham_scalapack_complex(molecorbdo_old, oversldo, psip, nelorb, nelorbu, &
                                                  neloccdo, info, descla, desch, size(oversldo, 1)/2, rankrep)
                end if
#else
                if (ipc .eq. 1) then
                    call graham(molecorb_old, overs, nelorbu, psip, nelorb, nelorbu, nelocc, info)
                    if (yeslsda) &
                        call graham(molecorbdo_old, overs, nelorbu, psip, nelorb, nelorbu, neloccdo, info)
                else
                    call graham_complex(molecorb_old, overs, nelorbu, psip, nelorb, nelorbu, nelocc, info)
                    call graham_complex(molecorbdo_old, oversdo, nelorbu, psip, nelorb, nelorbu, neloccdo, info)
                end if
#endif

#ifdef PARALLEL
                call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier

                molecorb(1:ipc*nelorbu, 1:bands) = molecorb_old(1:ipc*nelorbu, 1:bands)
                if (yeslsda .or. ipc .eq. 2) molecorbdo(1:ipc*nelorbu, 1:bands) = molecorbdo_old(1:ipc*nelorbu, 1:bands)

                call updenorb_new
                call uphamilt_new
                call evalevar

            end if

            ! restore molecorb modified by evalevar
            molecorb(1:ipc*nelorbu, 1:bands) = molecorb_old(1:ipc*nelorbu, 1:bands)
            if (yeslsda .or. ipc .eq. 2) molecorbdo(1:ipc*nelorbu, 1:bands) = molecorbdo_old(1:ipc*nelorbu, 1:bands)

        elseif (yespassed) then
            molecorb(1:ipc*nelorbu, 1:nelocc) = molecorbs(1:ipc*nelorbu, 1:nelocc)
            if (yeslsda .or. ipc .eq. 2) molecorbdo(1:ipc*nelorbu, 1:neloccdo) = molecorbdos(1:ipc*nelorbu, 1:neloccdo)
        end if

        if (rank .eq. 0 .and. .not. manyfort10) then
            write (6, *) ' Final molecorb written ', sum(molecorb(:, 1:nelocc))
            if (yeslsda .or. ipc .eq. 2) &
                write (6, *) ' Final molecorb down written ', sum(molecorbdo(:, 1:neloccdo))
        end if

        if (print_energies) then

            if (decoupled_run) then
                edft_av = 0.d0
                do i = 1, nk
                    edft_av = edft_av + dft_energies_all(i)
                end do
                edft_av = edft_av/nk
            else
                edft_av = edft
            end if

            if (rank .eq. 0) then
                write (unit_print_energies, *)
                write (unit_print_energies, '(a)') '# Total DFT energies after SC cycle '
                write (unit_print_energies, '(a,F15.8)') '# E_LDA averaged ', edft_av
                if (yeslsda .or. ipc .eq. 2) then
                    do ii = 1, nk
                        write (unit_print_energies, '(4F8.4,F15.8,2I6)') xkp(:, ii), wkp(ii), &
                            dft_energies_all(ii), number_particles(ii), number_particles_do(ii)
                    end do
                else
                    do ii = 1, nk
                        write (unit_print_energies, '(4F8.4,F15.8,I6)') xkp(:, ii), wkp(ii), &
                            dft_energies_all(ii), number_particles(ii)
                    end do
                end if
                write (unit_print_energies, *)
            end if
            deallocate (dft_energies_all)
            deallocate (number_particles)
            if (yeslsda .or. ipc .eq. 2) deallocate (number_particles_do)
        end if

        return

    end subroutine self_consistent_run

    subroutine cycle_complex

        !
        ! Density mixing using simple mixing (typeopt == 2) or
        ! accelerated Anderson/Broyden scheme (typeopt == 4).
        !

        implicit none

        real(8) :: normcorr, normcorrb, normcorra, cost, &
                   mixingnew, scallast, normlast, tmp_norm
        integer :: maxdim, ii, nsvd
        real(8), external :: ddot, cclock

#ifdef PARALLEL
        include 'mpif.h'
#endif

        timep = cclock()

        ! save the values of the density in the previous iteration(s)
        ! in order to perform the mixing
        if (typeopt .eq. 4) then

            itersto = itersto + 1

            if (itersto .gt. maxold) then
                ! reshuffhling
                itersto = maxold
                do ii = 1, maxold - 1
                    dent_before(:, ii) = dent_before(:, ii + 1)
                    dent_after(:, ii) = dent_after(:, ii + 1)
                    dent_aftern(:, ii) = dent_aftern(:, ii + 1)
                    mixingtrue(ii) = mixingtrue(ii + 1)
                end do
                if (yeslsda) then
                    do ii = 1, maxold - 1
                        spint_before(:, ii) = spint_before(:, ii + 1)
                        spint_after(:, ii) = spint_after(:, ii + 1)
                        spint_aftern(:, ii) = spint_aftern(:, ii + 1)
                    end do
                end if
            end if

            dent_before(:, itersto) = dent(:)
            if (yeslsda) spint_before(:, itersto) = spint(:)

        else

            if (typeopt .eq. 2) then
                dent_before(:, 1) = dent_before(:, 2)
                if (yeslsda) spint_before(:, 1) = spint_before(:, 2)
                if (iter .gt. 2) then
                    dent = dentnew ! restoring the (charge/spin) density
                    if (yeslsda) spint = spintnew
                end if
                dent_before(:, 2) = dent(:)
                if (yeslsda) spint_before(:, 2) = spint(:)
            end if

        end if
        !
        ! main cycle steps:
        ! 1) update hamiltonian with the new density
        ! 2) diagonalize the new hamiltonian
        ! 3) update electron occupations
        ! 4) compute the new density
        !

        call uphamilt_new
        !
        time = cclock()
        loading_time = loading_time + time - timep
        !
        call diagonalize_hamiltonian(0, lworkr)
        !
#ifdef DEBUG
        if (rank .eq. 0) then
            write (6, *)
            write (6, *) 'Eigenvalues Hamiltonian (spin up). Iteration/k-point:', iter, xkp(:, indk)
            do i = 1, bands
                write (6, *) i, eigmol(i)
            end do
            if (ipc .eq. 2 .or. yeslsda) then
                write (6, *) 'Eigenvalues Hamiltonian (spin down), Iteration/k-point:', iter, xkp(:, indk)
                do i = 1, bands
                    write (6, *) i, eigmoldo(i)
                end do
            end if
        end if
#endif

        call upocc_kpoints(optocc, epsshell, 0)
        !
        timep = cclock()
        diag_time = diag_time + timep - time

        if (typeopt .eq. 2 .and. iter .ge. 2) then

            call updenorb_new

            dent_after(:, 1) = dent_after(:, 2)
            dent_after(:, 2) = dent(:)

            if (yeslsda) then
                spint_after(:, 1) = spint_after(:, 2)
                spint_after(:, 2) = spint(:)
            end if

            dent(:) = dent_after(:, 2) - dent_before(:, 2) - dent_after(:, 1) + dent_before(:, 1)
            betaden = sum(dent(1:meshproc)**2)
            if (yeslsda) then
                spint(:) = spint_after(:, 2) - spint_before(:, 2) - spint_after(:, 1) + spint_before(:, 1)
                betaspin = sum(spint(1:meshproc)**2)
            end if
            !
            ! evaluation mixing coefficient
            !
#ifdef PARALLEL
            call reduce_base_real(1, betaden, commrep_mpi, -1)
            if (yeslsda) &
                call reduce_base_real(1, betaspin, commrep_mpi, -1)
#endif
            call sum_kpoints_scalar_real8(betaden, commcolrep_mpi, -1)
            if (yeslsda) call sum_kpoints_scalar_real8(betaspin, commcolrep_mpi, -1)
            !
            beta = sum((dent_after(1:meshproc, 2) - dent_before(1:meshproc, 2))*dent(1:meshproc))/betaden
            if (yeslsda) then
                betas = sum((spint_after(1:meshproc, 2) - spint_before(1:meshproc, 2))*spint(1:meshproc))/betaspin
            end if

#ifdef PARALLEL
            call reduce_base_real(1, beta, commrep_mpi, -1)
            if (yeslsda) &
                call reduce_base_real(1, betas, commrep_mpi, -1)
#endif
            dentold(:) = (1.d0 - beta)*dent_before(:, 2) + beta*dent_before(:, 1)
            dent(:) = (1.d0 - beta)*dent_after(:, 2) + beta*dent_after(:, 1)
            dentnew = mixing*dent + (1.d0 - mixing)*dentold

            if (yeslsda) then
                spintold(:) = (1.d0 - betas)*spint_before(:, 2) + betas*spint_before(:, 1)
                spint(:) = (1.d0 - betas)*spint_after(:, 2) + betas*spint_after(:, 1)
                spintnew = mixing*spint + (1.d0 - mixing)*spintold
            end if

            normcorrb = sum((dent_after(1:meshproc, 2) - dent_before(1:meshproc, 2))**2)/mesh
            if (yeslsda) normcorrb = normcorrb + &
                    &sum((spint_after(1:meshproc, 2) - spint_before(1:meshproc, 2))**2)/mesh
#ifdef PARALLEL
            call reduce_base_real(1, normcorrb, commrep_mpi, -1)
            tmp_norm = normcorrb
            call mpi_allreduce(tmp_norm, normcorrb, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                               commcolrep_mpi, ierr)
#endif
            if (rank .eq. 0 .and. optocc .ge. 0) write (6, *) ' Full norm correction =', normcorrb

        end if ! end if typeopt.eq.2

        molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
        if (yeslsda .or. ipc .eq. 2) molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)
        !
        ! Computation DFT energy after iteration
        !
        errdft = edft
        edft = sum(occupations(1:nelocc)*eigmol(1:nelocc))
        if (yeslsda .or. ipc .eq. 2) edft = edft + sum(occupationdo(1:neloccdo)*eigmoldo(1:neloccdo))

        edft = edft + totvpot
        !
        ! save DFT energies for each k-point
        !
        if (print_energies) then
            dft_energies_all = 0.d0
            corr_energy_all = 0.d0
            exchange_energy_all = 0.d0
            call collect_kpoints(dft_energies_all, nk, edft, -1)
            call collect_kpoints(corr_energy_all, nk, ecorr, -1)
            call collect_kpoints(exchange_energy_all, nk, exchange, -1)
            if (yeslsda) then
                spin_grid_all = 0.d0
                call collect_kpoints(spin_grid_all, nk, spingrid, -1)
            end if
        end if

        call sum_kpoints_scalar_real8(edft, commcolrep_mpi, -1)
        !
        oldocc = occupations
        if (yeslsda .or. ipc .eq. 2) oldoccdo = occupationdo
        !
        errdft = dabs(errdft - edft)
        etry = edft
        !
#ifdef PARALLEL
        if (decoupled_run) then
            max_errdft = 0.d0
            call mpi_allreduce(errdft, max_errdft, 1, MPI_DOUBLE_PRECISION, MPI_MAX, commcolrep_mpi, ierr)
            errdft = max_errdft
        end if
#endif
        !
        ! in the case of decoupled k-points calculations, show the arithmetic average of the energy
        ! among all k-points.
        !
        if (decoupled_run) then
            edft_av = sum(dft_energies_all(:))/nk
            ecorr_av = sum(corr_energy_all(:))/nk
            exchange_av = sum(exchange_energy_all(:))/nk
            if (yeslsda) spingrid_av = sum(spin_grid_all(:))/nk
        else
            edft_av = edft
            ecorr_av = ecorr
            exchange_av = exchange
            if (yeslsda) spingrid_av = spingrid
        end if

        if (rank .eq. 0) write (6, *)
        if (yeslsda) then
            if (rank .eq. 0) then
                write (6, '(A27,I6,5f18.7)') &
                    ' Iter,E,xc,corr, |spin|, errdft = ' &
                    , iter, edft_av, exchange_av, ecorr_av, spingrid_av/2, errdft
            end if
        else
            if (rank .eq. 0) then
                write (6, '(A15,I6,4f18.7)') &
                    ' Iter,E,xc,corr,errdft= ' &
                    , iter, edft_av, exchange_av, ecorr_av, errdft
            end if
        end if
        if (rank .eq. 0) write (6, *)

        loading_time = loading_time + cclock() - timep

        if (typeopt .eq. 4) then

            mixingtrue(itersto) = mixingder

            call updenorb_new
            call cutdens

            normcorrb = sum((dent(1:meshproc) - dent_before(1:meshproc, itersto))**2)/mesh
            if (yeslsda) normcorrb = normcorrb + &
                    &sum((spint(1:meshproc) - spint_before(1:meshproc, itersto))**2)/mesh
#ifdef PARALLEL
            call reduce_base_real(1, normcorrb, commrep_mpi, -1)
#endif
            if (rank .eq. 0 .and. optocc .ge. 0) write (6, *) ' Full norm correction before =', normcorrb
            !        if(rankrep.eq.0.and.optocc.ge.0) write(6,*) ' Full norm correction before =',normcorrb,rankcolrep
            dent_after(:, itersto) = (1.d0 - mixingder)*dent_before(:, itersto) + mixingder*dent(:)

            if (yeslsda) then
                spint_after(:, itersto) = (1.d0 - mixingder)*spint_before(:, itersto)&
                        & + mixingder*spint(:)
            end if

            dent(:) = dent_after(:, itersto)
            if (yeslsda) spint(:) = spint_after(:, itersto)

            ! 1)
            call uphamilt_new

            time = cclock()

            loading_time = loading_time + time - timep

            eigmolo = eigmol
            if (yeslsda .or. ipc .eq. 2) eigmolodo = eigmoldo
            ! 2)
            call diagonalize_hamiltonian(0, lworkr)
            ! 3)
            call upocc_kpoints(optocc, epsshell, 0)
            ! 4)
            call updenorb_new

            normcorra = sum((dent(1:meshproc) - dent_after(1:meshproc, itersto))**2)/mesh
            if (yeslsda) normcorra = normcorra + &
                    &sum((spint(1:meshproc) - spint_after(1:meshproc, itersto))**2)/mesh
#ifdef PARALLEL
            call reduce_base_real(1, normcorra, commrep_mpi, -1)
#endif
            if (rank .eq. 0 .and. optocc .ge. 0) write (6, *) ' Full norm correction after  =', normcorra
            !        if(rankrep.eq.0.and.optocc.ge.0) write(6,*) ' Full norm correction after  =',normcorra,rankcolrep
            eigmol = eigmolo
            occupations = oldocc
            if (yeslsda .or. ipc .eq. 2) then
                eigmoldo = eigmolodo
                occupationdo = oldoccdo
            end if

            ! -------------------------------------------
            ! perform mixing between dent and dent_after
            ! -------------------------------------------

            dent_aftern(:, itersto) = mixingder*dent(:) + (1.d0 - mixingder)*dent_after(:, itersto)

            if (yeslsda) spint_aftern(:, itersto) = mixingder*spint(:)&
                    & + (1.d0 - mixingder)*spint_after(:, itersto)

            call dgemm('T', 'N', itersto, itersto, meshproc, volmesh, dent_after, size(dent_after, 1)&
                    &, dent_aftern, size(dent_aftern, 1), 0.d0, jac, dimjac)

            call dgemm('T', 'N', itersto, itersto, meshproc, -volmesh, dent_after, size(dent_after, 1)&
                    &, dent_after, size(dent_after, 1), 1.d0, jac, dimjac)

            if (yeslsda) then
                call dgemm('T', 'N', itersto, itersto, meshproc, volmesh, spint_after&
                        &, size(spint_after, 1), spint_aftern, size(spint_aftern, 1), 1.d0, jac, dimjac)
                call dgemm('T', 'N', itersto, itersto, meshproc, -volmesh, spint_after&
                        &, size(spint_after, 1), spint_after, size(spint_after, 1), 1.d0, jac, dimjac)
            end if

            do ii = 1, itersto
                jac(:, ii) = jac(:, ii)*mixingtrue(itersto)/mixingtrue(ii)
            end do

            call dgemm('T', 'N', itersto, itersto, meshproc, volmesh, dent_after&
                    &, size(dent_after, 1), dent_after, size(dent_after, 1), 1.d0, jac, dimjac)

            if (yeslsda) then
                call dgemm('T', 'N', itersto, itersto, meshproc, volmesh, spint_after&
                        &, size(spint_after, 1), spint_after, size(spint_after, 1), 1.d0, jac, dimjac)
            end if

            call dgemm('T', 'N', itersto, itersto, meshproc, volmesh, dent_after&
                    &, size(dent_after, 1), dent_after, size(dent_after, 1), 0.d0, overjac, dimjac)

            if (yeslsda) then
                call dgemm('T', 'N', itersto, itersto, meshproc, volmesh, spint_after&
                        &, size(spint_after, 1), spint_after, size(spint_after, 1), 1.d0, overjac, dimjac)
            end if

            !      shifting the trivial vectors
            call dgemv('T', meshproc, itersto, volmesh, dent_aftern, size(dent_aftern, 1)&
                    &, dent_before(1, itersto), 1, 0.d0, vetaux, 1)

            call dgemv('T', meshproc, itersto, -volmesh, dent_after, size(dent_after, 1)&
                    &, dent_before(1, itersto), 1, 1.d0, vetaux, 1)

            if (yeslsda) then
                call dgemv('T', meshproc, itersto, volmesh, spint_aftern, size(spint_aftern, 1)&
                        &, spint_before(1, itersto), 1, 1.d0, vetaux, 1)
                call dgemv('T', meshproc, itersto, -volmesh, spint_after, size(spint_after, 1)&
                        &, spint_before(1, itersto), 1, 1.d0, vetaux, 1)
            end if

            do ii = 1, itersto
                vetaux(ii, 1) = vetaux(ii, 1)*mixingtrue(itersto)/mixingtrue(ii)
            end do

            call dgemv('T', meshproc, itersto, volmesh, dent_after, size(dent_after, 1)&
                    &, dent_before(1, itersto), 1, 1.d0, vetaux, 1)

            if (yeslsda) then
                call dgemv('T', meshproc, itersto, volmesh, spint_after, size(spint_after, 1)&
                        &, spint_before(1, itersto), 1, 1.d0, vetaux, 1)
            end if

            vetaux(1:itersto, 1) = vetaux(1:itersto, 1)&
                    & - volmesh*ddot(meshproc, dent_before(1, itersto), 1, dent_after(1, itersto), 1)

            if (yeslsda) then
                vetaux(1:itersto, 1) = vetaux(1:itersto, 1)&
                        & - volmesh*ddot(meshproc, spint_before(1, itersto), 1, spint_after(1, itersto), 1)
            end if

            call dgemv('T', meshproc, itersto, volmesh, dent_after, size(dent_after, 1)&
                    &, dent_after(1, itersto), 1, 0.d0, vetaux(1, 2), 1)

            if (yeslsda) then
                call dgemv('T', meshproc, itersto, volmesh, spint_after, size(spint_after, 1)&
                        &, spint_after(1, itersto), 1, 1.d0, vetaux(1, 2), 1)
            end if

            do ii = 1, itersto
                jac(ii, 1:itersto) = jac(ii, 1:itersto) - vetaux(1:itersto, 1)
            end do

            do ii = 1, itersto
                jac(1:itersto, ii) = jac(1:itersto, ii) - vetaux(1:itersto, 2)
            end do

            ! The same for the overlap matrix
            call dgemv('T', meshproc, itersto, volmesh, dent_after, size(dent_after, 1)&
                    &, dent_before(1, itersto), 1, 0.d0, vetaux, 1)

            if (yeslsda) then
                call dgemv('T', meshproc, itersto, volmesh, spint_after, size(spint_after, 1)&
                        &, spint_before(1, itersto), 1, 1.d0, vetaux, 1)
            end if

            vetaux(1:itersto, 2) = vetaux(1:itersto, 1)

            vetaux(1:itersto, 1) = vetaux(1:itersto, 1)&
                    & - volmesh*ddot(meshproc, dent_before(1, itersto), 1, dent_before(1, itersto), 1)

            if (yeslsda) then
                vetaux(1:itersto, 1) = vetaux(1:itersto, 1)&
                        & - volmesh*ddot(meshproc, spint_before(1, itersto), 1, spint_before(1, itersto), 1)
            end if

            do ii = 1, itersto
                overjac(ii, 1:itersto) = overjac(ii, 1:itersto) - vetaux(1:itersto, 1)
            end do

            do ii = 1, itersto
                overjac(1:itersto, ii) = overjac(1:itersto, ii) - vetaux(1:itersto, 2)
            end do

#ifdef PARALLEL
            call reduce_base_real(dimjac2, jac, commrep_mpi, -1)
            call reduce_base_real(dimjac2, overjac, commrep_mpi, -1)
#endif

            ! Matrix to be inverted

            ! More consistent deceleration for mixing-->0 sstandard mixing small.
            if (mixing .ne. 1.d0) then
                cost = 1.d0 + (1.d0/mixing - 1.d0)*mixingder
                jac(1:itersto, 1:itersto) = jac(1:itersto, 1:itersto) - &
                        &cost*overjac(1:itersto, 1:itersto)
            else
                jac(1:itersto, 1:itersto) = jac(1:itersto, 1:itersto) - &
                        &overjac(1:itersto, 1:itersto)
            end if

            ! Known term

            !        maxold=1
            maxdim = itersto

            if (jaccond .gt. 0.d0) then
                ! first diagonalize the overlap matrix
                call dsyev('V', 'U', maxdim, overjac, dimjac, sojac, psip, lworkjac, info)

                if (rank .eq. 0) write (6, *) ' Overlap cond num. =', sojac(1)/sojac(maxdim)
                do ii = 1, maxdim - 1
                    if (sojac(ii)/sojac(maxdim) .lt. jaccond) sojac(ii) = 0.d0
                end do
                !        change basis for jac
                call dgemm('N', 'N', maxdim, maxdim, maxdim, 1.d0, jac, dimjac, overjac&
                        &, dimjac, 0.d0, mataux, dimjac)

                call dgemm('T', 'N', maxdim, maxdim, maxdim, 1.d0, overjac, dimjac, mataux&
                        &, dimjac, 0.d0, jac, dimjac)

                do ii = 1, maxdim
                    if (sojac(ii) .ne. 0.d0) then
                        sjac(ii) = 1.d0/dsqrt(sojac(ii))
                    else
                        sjac(ii) = 0.d0
                    end if
                end do

                do ii = 1, maxdim
                    do jj = 1, maxdim
                        jac(ii, jj) = jac(ii, jj)*sjac(ii)*sjac(jj)
                    end do
                    if (sojac(ii) .gt. 0.d0) then
                        alpha(ii) = -dsqrt(sojac(ii))*overjac(maxdim, ii)
                    else
                        alpha(ii) = 0.d0
                    end if
                end do

                ! SVD of matrix jac
                call dgesvd('A', 'A', maxdim, maxdim, jac, dimjac, sjac, ujac, dimjac&
                        &, vjac, dimjac, psip, lworkjac, info)
                !        now applying the inverse to alpha
                call dgemv('T', maxdim, maxdim, 1.d0, ujac, dimjac, alpha, 1, 0.d0, psip, 1)
                !        Removing singular directions
                if (rank .eq. 0) write (6, *) ' svd  jac =', (sjac(ii), ii=1, maxdim)
                if (rank .eq. 0) write (6, *) 'svd   overlap left/right=', &
                        &(sum(ujac(1:maxdim, ii)*vjac(1:maxdim, ii)), ii=1, maxdim)

                nsvd = 0
                do ii = 2, maxdim
                    if (sjac(ii)/sjac(1) .lt. jaccond) then
                        sjac(ii) = 0.d0
                        nsvd = nsvd + 1
                    else
                        sjac(ii) = 1.d0/sjac(ii)
                    end if
                end do
                sjac(1) = 1.d0/sjac(1)
                if (rank .eq. 0 .and. nsvd .ne. 0) write (6, *) ' Warning removing', nsvd&
                        &, 'singular direction jac '

                alpha(1:maxdim) = psip(1:maxdim)*sjac(1:maxdim)

                call dgemv('T', maxdim, maxdim, 1.d0, vjac, dimjac, alpha, 1, 0.d0, psip, 1)
                !        Now going back to the original basis
                do ii = 1, maxdim
                    if (sojac(ii) .ne. 0.d0) then
                        psip(ii) = psip(ii)/dsqrt(sojac(ii))
                    else
                        psip(ii) = 0.d0
                    end if
                end do

                call dgemv('N', maxdim, maxdim, 1.d0, overjac, dimjac, psip, 1, 0.d0, alpha, 1)

            else

                alpha(1:itersto) = -overjac(1:itersto, itersto)
                call dgetrf(maxdim, maxdim, jac, dimjac, ipsip, info)
                if (rank .eq. 0) write (6, *) ' jac diag=', (jac(ii, ii), ii=1, maxdim)
                call dgetrs('N', maxdim, 1, jac, dimjac, ipsip, alpha, maxdim, info)

            end if

            mixingnew = normcorrb/normcorra
            dent(:) = 0.d0

            do ii = 1, itersto
                dent(:) = dent(:) + alpha(ii)*(dent_after(:, ii) - dent_before(:, itersto))
            end do

            if (yeslsda) then
                spint(:) = 0.d0

                do ii = 1, itersto
                    spint(:) = spint(:) + alpha(ii)*(spint_after(:, ii) - spint_before(:, itersto))
                end do

            end if
            !
            normcorr = sum(dent(1:meshproc)**2)/mesh
            if (maxold .gt. 1) then
                normlast = sum((dent_after(1:meshproc, itersto) - dent_before(1:meshproc, itersto))**2)/mesh
                scallast = sum(dent(1:meshproc)*(dent_after(1:meshproc, itersto) - dent_before(1:meshproc, itersto)))/mesh
            end if
            if (yeslsda) then
                normcorr = normcorr + sum(spint(1:meshproc)**2)/mesh
                if (maxold .gt. 1) then
                    scallast = scallast + &
                            &sum(spint(1:meshproc)*(spint_after(1:meshproc, itersto) - spint_before(1:meshproc, itersto)))/mesh
                    normlast = normlast + &
                            &sum((spint_after(1:meshproc, itersto) - spint_before(1:meshproc, itersto))**2)/mesh
                end if
            end if

#ifdef PARALLEL
            call reduce_base_real(1, normcorr, commrep_mpi, -1)
            if (manyfort10) then
                tmp_norm = normcorr
                call mpi_allreduce(tmp_norm, normcorr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                                   commcolrep_mpi, ierr)
            end if
            if (maxold .gt. 1) then
                call reduce_base_real(1, scallast, commrep_mpi, -1)
                call reduce_base_real(1, normlast, commrep_mpi, -1)
            end if
#endif

            if (normcorr .gt. 0) then
                normcorr = dsqrt(normcorr)
            else
                normcorr = 0.d0
            end if

            if (changemix) then
                if (maxold .eq. 1) then
                    if (alpha(itersto) .gt. 0) then
                        mixingder = (alpha(itersto)*mixing*mixingder + mixingder)/2.d0
                        if (rank .eq. 0) write (6, *) ' Warning new mixingder =', mixingder
                    end if
                else
                    if (scallast .gt. 0.d0) then
                        mixingder = (scallast/normlast*mixingder + mixingder)/2.d0
                        if (rank .eq. 0) write (6, *) ' Warning new mixingder =', mixingder
                    end if
                end if
            end if

            if (rank .eq. 0) then
                write (6, '(A12,I6,1f20.10)') ' Norm corr.=', iter, normcorr
            end if

            if (rank .eq. 0) then
                if (itersto .eq. 1) then
                    write (6, *) ' Warning simple mixing =', alpha(itersto)*mixingder
                else
                    write (6, *) ' Warning last mixing =', alpha(itersto)
                end if
            end if

            dent(:) = dent_before(:, itersto) + dent(:)
            if (yeslsda) spint(:) = spint_before(:, itersto) + spint(:)

            ! restore molecorb anyway for final output molecorb,eigmol
            molecorb(1:ipc*nelorbu, 1:bands) = molecorb_old(1:ipc*nelorbu, 1:bands)
            if (yeslsda .or. ipc .eq. 2) molecorbdo(1:ipc*nelorbu, 1:bands) = molecorbdo_old(1:ipc*nelorbu, 1:bands)

        end if ! end if typeopt.eq.4 - DEFAULT mixing scheme

        ! --------------------------
        ! End default MIXING scheme
        ! --------------------------

        if (typeopt .eq. 0 .or. (typeopt .eq. 2 .and. iter .le. 1)) then
            dentold = dent
            if (yeslsda) spintold = spint
        end if

        if ((iter .le. 1 .and. typeopt .eq. 2) .or. typeopt .eq. 0) then

            call updenorb_new

            ! the first two iterations make the standard stuff: simple mixing between
            ! new density and uniform density, defined at the beginning.
            normcorrb = sum((dent(1:meshproc) - dentold(1:meshproc))**2)/mesh
            if (yeslsda) normcorrb = normcorrb + &
                    &sum((spint(1:meshproc) - spintold(1:meshproc))**2)/mesh
#ifdef PARALLEL
            call reduce_base_real(1, normcorrb, commrep_mpi, -1)
#endif
            if (rank .eq. 0 .and. optocc .ge. 0) write (6, *) ' Full norm correction =', normcorrb
            ! low mixing at first 2 iterations
            if (typeopt .eq. 2) then
                dent = 0.25d0*mixing*dent + (1.d0 - 0.25d0*mixing)*dentold
                if (yeslsda) spint = 0.25d0*mixing*spint + (1.d0 - 0.25d0*mixing)*spintold
            else
                dent = mixing*dent + (1.d0 - mixing)*dentold
                if (yeslsda) spint = mixing*spint + (1.d0 - mixing)*spintold
            end if

        end if

        call cutdens

        occupations = oldocc
        if (yeslsda .or. ipc .eq. 2) occupationdo = oldoccdo

        return

    end subroutine cycle_complex

    subroutine cyclesteep_complex

        implicit none

        !
        ! Density mixing using steepest-descent method
        !

        real*8, external :: ddot, dnrm2, cclock
        complex*16, external :: zdotu, zdotc_
#ifdef __SCALAPACK
        integer :: irow, icol
#endif
#ifdef PARALLEL
        include 'mpif.h'
#endif

        timep = cclock()

        call updenorb_new
        call cutdens
        call uphamilt_new

#ifdef __SCALAPACK
        molecorb_old = 0.d0

        if (descla(lambda_node_) > 0) then

            allocate (molecorb_part(ipc*descla(nlax_), descla(nlax_)))
            molecorb_part = 0.d0
            irow = descla(ilar_)
            icol = descla(ilac_)
            do j = 1, descla(nlac_)
                if ((j + icol - 1) <= bands) then
                    do i = 1, descla(nlar_)
                        if (ipc .eq. 1) then
                            molecorbl(i, j) = molecorb((i + irow - 1), (j + icol - 1))
                        else
                            molecorbl(2*i - 1:2*i, j) = molecorb(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1))
                        end if
                    end do
                end if
            end do

            if (ipc .eq. 1) then
                call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.d0, hamiltl, 1, 1,&
                     & desch, molecorbl, 1, 1, desch, 0.d0, molecorb_part, 1, 1, desch)
            else
                call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, (1.0d0, 0.d0), hamiltl, 1, 1,&
                     & desch, molecorbl, 1, 1, desch, (0.0d0, 0.d0), molecorb_part, 1, 1, desch)
            end if

            irow = descla(ilar_)
            icol = descla(ilac_)
            do j = 1, descla(nlac_)
                if ((j + icol - 1) <= bands) then
                    do i = 1, descla(nlar_)
                        if (ipc .eq. 1) then
                            molecorb_old((i + irow - 1), (j + icol - 1)) = molecorb_part(i, j)
                        else
                            molecorb_old(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1)) = molecorb_part(2*i - 1:2*i, j)
                        end if
                    end do
                end if
            end do
            deallocate (molecorb_part)
        end if
#ifdef PARALLEL
        call reduce_base_real(size(molecorb_old), molecorb_old, commrep_mpi, -1)
#endif

        if (yeslsda .or. ipc .eq. 2) then
            molecorbdo_old = 0.d0
            if (descla(lambda_node_) > 0) then
                allocate (molecorb_part(ipc*descla(nlax_), descla(nlax_)))
                molecorb_part = 0.d0

                irow = descla(ilar_)
                icol = descla(ilac_)

                do j = 1, descla(nlac_)
                    if ((j + icol - 1) <= bands) then
                        do i = 1, descla(nlar_)
                            if (ipc .eq. 1) then
                                molecorbldo(i, j) = molecorbdo((i + irow - 1), (j + icol - 1))
                            else
                                molecorbldo(2*i - 1:2*i, j) = molecorbdo(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1))
                            end if
                        end do
                    end if
                end do

                if (ipc .eq. 1) then
                    call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.d0, hamiltldo, 1, 1,&
                         & desch, molecorbldo, 1, 1, desch, 0.d0, molecorb_part, 1, 1, desch)
                else
                    call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, (1.0d0, 0.d0), hamiltldo, 1, 1,&
                         & desch, molecorbldo, 1, 1, desch, (0.0d0, 0.d0), molecorb_part, 1, 1, desch)
                end if
                irow = descla(ilar_)
                icol = descla(ilac_)
                do j = 1, descla(nlac_)
                    if ((j + icol - 1) <= bands) then
                        do i = 1, descla(nlar_)
                            if (ipc .eq. 1) then
                                molecorbdo_old((i + irow - 1), (j + icol - 1)) = molecorb_part(i, j)
                            else
                                molecorbdo_old(2*(i + irow - 1) &
                                               - 1:2*(i + irow - 1), (j + icol - 1)) = molecorb_part(2*i - 1:2*i, j)
                            end if
                        end do
                    end if
                end do
                deallocate (molecorb_part)
            end if
#ifdef PARALLEL
            call reduce_base_real(size(molecorbdo_old), molecorbdo_old, commrep_mpi, -1)
#endif
        end if

#else

        !       Steepest desceent
        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, hamilt, nelorbu&
                    &, molecorb, nelorbu, 0.d0, molecorb_old, nelorb)
            if (yeslsda) then
                call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, hamiltdo, nelorbu&
                        &, molecorbdo, nelorbu, 0.d0, molecorbdo_old, nelorb)
            end if
        else
            call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, hamilt, nelorbu&
                    &, molecorb, nelorbu, zzero, molecorb_old, nelorb)
            call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, hamiltdo, nelorbu&
                    &, molecorbdo, nelorbu, zzero, molecorbdo_old, nelorb)
        end if
#endif

        do i = 1, bands
            molecorb_old(:, i) = occupations(i)*molecorb_old(:, i)
        end do
        if (yeslsda .or. ipc .eq. 2) then
            do i = 1, bands
                molecorbdo_old(:, i) = occupationdo(i)*molecorbdo_old(:, i)
            end do
        end if

#ifdef __SCALAPACK
        allocate (molecorb_part(ipc*nelorbu, bands))
        molecorb_part = 0.d0
        if (descla(lambda_node_) > 0) then
            !        destroy overs mat_in
            if (ipc .eq. 1) then
                call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.0d0, &
                            oversl, 1, 1, desch, molecorbl, 1, 1, desch, 0.0d0, hamiltl, 1, 1, desch)
            else
                call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, (1.0d0, 0.d0), &
                            oversl, 1, 1, desch, molecorbl, 1, 1, desch, (0.0d0, 0.d0), hamiltl, 1, 1, desch)
            end if
            irow = descla(ilar_)
            icol = descla(ilac_)
            do j = 1, descla(nlac_)
                if ((j + icol - 1) <= bands) then
                    do i = 1, descla(nlar_)
                        if (ipc .eq. 1) then
                            molecorb_part((i + irow - 1), (j + icol - 1)) = hamiltl(i, j)
                        else
                            molecorb_part(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1)) = hamiltl(2*i - 1:2*i, j)
                        end if
                    end do
                end if
            end do
        end if
#ifdef PARALLEL
        call reduce_base_real(size(molecorb_part), molecorb_part, commrep_mpi, -1)
#endif

        errdft = edft
        edft = 0.d0

        do i = 1, nelocc
            if (ipc .eq. 1) then
                eigmol(i) = ddot(nelorbu, molecorb(1, i), 1, molecorb_old(1, i), 1)
            else
                eigmol(i) = zdotc_(nelorbu, molecorb(1, i), 1, molecorb_old(1, i), 1)
            end if
            molecorb_old(1:ipc*nelorbu, i) = molecorb_old(1:ipc*nelorbu, i) - eigmol(i)*molecorb_part(1:ipc*nelorbu, i)
            edft = edft + eigmol(i)
        end do

        if (yeslsda .or. ipc .eq. 2) then

            molecorb_part = 0.d0
            if (descla(lambda_node_) > 0) then

                if (ipc .eq. 1) then
                    call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.d0, &
                                oversl, 1, 1, desch, molecorbldo, 1, 1, desch, 0.d0, hamiltldo, 1, 1, desch)
                else
                    call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, (1.0d0, 0.d0), &
                                oversl, 1, 1, desch, molecorbldo, 1, 1, desch, (0.0d0, 0.d0), hamiltldo, 1, 1, desch)
                end if
                irow = descla(ilar_)
                icol = descla(ilac_)
                do j = 1, descla(nlac_)
                    if ((j + icol - 1) <= bands) then
                        do i = 1, descla(nlar_)
                            if (ipc .eq. 1) then
                                molecorb_part((i + irow - 1), (j + icol - 1)) = hamiltldo(i, j)
                            else
                                molecorb_part(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1)) = hamiltldo(2*i - 1:2*i, j)
                            end if
                        end do
                    end if
                end do
            end if
#ifdef PARALLEL
            call reduce_base_real(size(molecorb_part), molecorb_part, commrep_mpi, -1)
#endif

            do i = 1, neloccdo
                if (ipc .eq. 1) then
                    eigmoldo(i) = ddot(nelorbu, molecorbdo(1, i), 1, molecorbdo_old(1, i), 1)
                else
                    eigmoldo(i) = zdotc_(nelorbu, molecorbdo(1, i), 1, molecorbdo_old(1, i), 1)
                end if
                molecorbdo_old(1:ipc*nelorbu, i) = molecorbdo_old(1:ipc*nelorbu, i) - eigmoldo(i)*molecorb_part(1:ipc*nelorbu, i)
                edft = edft + eigmoldo(i)
            end do

        end if

        deallocate (molecorb_part)

#else
        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, overs, nelorbu, molecorb&
                    &, nelorbu, 0.d0, hamilt, nelorbu)
        else
            call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, overs, nelorbu, molecorb&
                &, nelorbu, zzero, hamilt, nelorbu)
        end if

        errdft = edft
        edft = 0.d0
        do i = 1, nelocc
        if (ipc .eq. 1) then
            eigmol(i) = ddot(nelorbu, molecorb(1, i), 1, molecorb_old(1, i), 1)
        else
            eigmol(i) = zdotc_(nelorbu, molecorb(1, i), 1, molecorb_old(1, i), 1)
        end if
        molecorb_old(1:ipc*nelorbu, i) = molecorb_old(1:ipc*nelorbu, i) - eigmol(i)*hamilt(1:ipc*nelorbu, i)
        edft = edft + eigmol(i)
        end do
        if (yeslsda .or. ipc .eq. 2) then
        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, overs, nelorbu, molecorbdo&
                &, nelorbu, 0.d0, hamiltdo, nelorbu)
        else
            call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, overs, nelorbu, molecorbdo&
                &, nelorbu, zzero, hamiltdo, nelorbu)
        end if

        do i = 1, neloccdo
        if (ipc .eq. 1) then
            eigmoldo(i) = ddot(nelorbu, molecorbdo(1, i), 1, molecorbdo_old(1, i), 1)
        else
            eigmoldo(i) = zdotc_(nelorbu, molecorbdo(1, i), 1, molecorbdo_old(1, i), 1)
        end if
        molecorbdo_old(1:ipc*nelorbu, i) = molecorbdo_old(1:ipc*nelorbu, i) - eigmoldo(i)*hamiltdo(1:ipc*nelorbu, i)
        edft = edft + eigmoldo(i)
        end do

        end if
#endif
        ! average the energy over the k-points
        edft = edft + totvpot

        if (print_energies) then
            dft_energies_all = 0.d0
            corr_energy_all = 0.d0
            exchange_energy_all = 0.d0
            call collect_kpoints(dft_energies_all, nk, edft, -1)
            call collect_kpoints(corr_energy_all, nk, ecorr, -1)
            call collect_kpoints(exchange_energy_all, nk, exchange, -1)
            if (yeslsda) then
                spin_grid_all = 0.d0
                call collect_kpoints(spin_grid_all, nk, spingrid, -1)
            end if
        end if

        call sum_kpoints_scalar_real8(edft, commcolrep_mpi, -1)
        errdft = dabs(errdft - edft)

#ifdef PARALLEL
        if (decoupled_run) then
            max_errdft = 0.d0
            call mpi_allreduce(errdft, max_errdft, 1, MPI_DOUBLE_PRECISION, MPI_MAX, commcolrep_mpi, ierr)
            errdft = max_errdft
        end if
#endif

        !
        ! in the case of decoupled k-points calculations, show the arithmetic average of the energy
        ! among all k-points.
        !
        if (decoupled_run) then
            edft_av = sum(dft_energies_all(:))/nk
            ecorr_av = sum(corr_energy_all(:))/nk
            exchange_av = sum(exchange_energy_all(:))/nk
            if (yeslsda) spingrid_av = sum(spin_grid_all(:))/nk
        else
            edft_av = edft
            ecorr_av = ecorr
            exchange_av = exchange
            if (yeslsda) spingrid_av = spingrid
        end if

        if (mixingstep .ne. 0.d0) then

            if (yeslsda) then
                if (rank .eq. 0) then
                    write (6, '(A27,I6,5f18.7)') &
                        ' Iter,E,xc,corr, |spin| = ' &
                        , iter, edft_av, exchange_av, ecorr_av, spingrid_av/2, errdft
                end if
            else
                if (rank .eq. 0) then
                    write (6, '(A15,I6,4f18.7)') &
                        ' Iter,E,xc,corr= ' &
                        , iter, edft_av, exchange_av, ecorr_av, errdft
                end if
            end if

            do i = 1, nelocc
            if (ipc .eq. 1) then
                call daxpy(nelorbu, -mixingstep, molecorb_old(1, i), 1, molecorb(1, i), 1)
            else
                call zaxpy(nelorbu, -mixingstep, molecorb_old(1, i), 1, molecorb(1, i), 1)
            end if
            end do
            if (yeslsda .or. ipc .eq. 2) then
            do i = 1, neloccdo
            if (ipc .eq. 1) then
                call daxpy(nelorbu, -mixingstep, molecorbdo_old(1, i), 1, molecorbdo(1, i), 1)
            else
                call zaxpy(nelorbu, -mixingstep, molecorbdo_old(1, i), 1, molecorbdo(1, i), 1)
            end if
            end do
            end if
            !
            ! orthogonalization molecular orbitals
            !
#ifdef __SCALAPACK
            if (ipc .eq. 1) then
                call graham_scalapack(molecorb, oversl, psip, nelorbu, nelorbu, &
                                      nelocc, info, descla, desch, size(oversl, 1), rankrep)
                if (yeslsda) &
                    call graham_scalapack(molecorbdo, oversl, psip, nelorbu, nelorbu, &
                                          neloccdo, info, descla, desch, size(oversl, 1), rankrep)
            else
                call graham_scalapack_complex(molecorb, oversl, psip, nelorbu, nelorbu, &
                                              nelocc, info, descla, desch, size(oversl, 1)/2, rankrep)
                call graham_scalapack_complex(molecorbdo, oversldo, psip, nelorbu, nelorbu, &
                                              neloccdo, info, descla, desch, size(oversldo, 1)/2, rankrep)
            end if
#else
            if (ipc .eq. 1) then
                call graham(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                if (yeslsda) &
                    call graham(molecorbdo, overs, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
            else
                call graham_complex(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
                call graham_complex(molecorbdo, oversdo, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
            end if
#endif

#ifdef PARALLEL
            call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier

        end if
        !       compute the new sum of eigenvalue
        !       check orthogonality
        if (info .ne. 0 .and. rank .eq. 0) write (6, *) ' Warning info>0 after graham ', info

#ifdef DEBUG
        !       check orthogonality
        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbh, nelup, nelorbh, 1.d0, overs, nelorb&
                 &, molecorb, nelorb, 0.d0, hamilt, nelorb)
        else
            call zgemm('N', 'N', nelorbh, nelup, nelorbh, zone, overs, nelorb&
                 &, molecorb, nelorb, zzero, hamilt, nelorb)
        end if
        if (rank .eq. 0) write (6, *) ' scalar products after graham '
        do i = 1, nelup
            do j = i, nelup
                if (ipc .eq. 1) then
                    if (rank .eq. 0) write (6, *) i, j, ddot(nelorbh, molecorb(1, i), 1, hamilt(1, j), 1)
                else
                    if (rank .eq. 0) write (6, *) i, j, zdotc_(nelorbh, molecorb(1, i), 1, hamilt(1, j), 1)
                end if
            end do
        end do
#endif
        loading_time = loading_time + cclock() - timep
        return

    end subroutine cyclesteep_complex

    subroutine checksum

        implicit none
        real(8) costs, costd, costh, costo, costodo, costhup, costhdo, costhodo
#ifdef PARALLEL
        include 'mpif.h'
#endif

        costh = 0.d0
        costo = 0.d0
        costhup = 0.d0
        costhdo = 0.d0
        costodo = 0.d0
        costhodo = 0.d0
        costs = 0.d0
        costd = 0.d0

        ! check charge/spin density sum on the grid
        if (yeslsda) costs = sum(spint(1:meshproc)**2)
        costd = sum(dent(1:meshproc)**2)

#ifdef __SCALAPACK
!         write(6,*) ' # Processor ',rank,descla( lambda_node_ ),ir,ic,nelorb
!         write(6,*) ' Small matrix limits ',rank,descla( nlac_ ),descla( nlar_ )
        if (descla(lambda_node_) > 0) then
            ir = descla(ilar_)
            ic = descla(ilac_)
            do j = 1, descla(nlac_)
                do i = 1, descla(nlar_)
                    if (i + ir - 1 .le. nelorb .and. i + ir - 1 .ge. 1 .and. j + ic - 1 .le. nelorb .and. j + ic - 1 .ge. 1) then
                        if (ipc .eq. 1) then
                            costh = costh + overhaml(i, j)
                            costo = costo + oversl(i, j)
                            costhup = costhup + hamiltl(i, j)
                            if (yeslsda) &
                                costhdo = costhdo + hamiltldo(i, j)
                        else
                            costh = costh + overhaml(2*i - 1, j) + overhaml(2*i, j)
                            costo = costo + oversl(2*i, j) + oversl(2*i - 1, j)
                            costhup = costhup + hamiltl(2*i - 1, j) + hamiltl(2*i, j)
                            costhdo = costhdo + hamiltldo(2*i - 1, j) + hamiltldo(2*i, j)
                            costodo = costodo + oversldo(2*i - 1, j) + oversldo(2*i, j)
                            costhodo = costhodo + overhamldo(2*i - 1, j) + overhamldo(2*i, j)
                        end if
                    end if
                end do
            end do
        end if

#ifdef PARALLEL
        call mpi_reduce(costh, sumh, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
        call mpi_reduce(costo, sumo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
        call mpi_reduce(costhup, sumhup, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
        call mpi_reduce(costhdo, sumhdo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
        call mpi_reduce(costhodo, sumhodo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
        call mpi_reduce(costodo, sumodo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
        call mpi_reduce(costs, sumspin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
        call mpi_reduce(costd, sumden, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
#else
        costh = sum(overham(:, :))
        costo = sum(overs(:, :))
        costhup = sum(hamilt(:, :))
        if (yeslsda) costhdo = sum(hamiltdo(:, :))
        if (ipc .eq. 2) then
            costodo = sum(oversdo(:, :))
            costhodo = sum(overhamdo(:, :))
        end if

        sumh = costh
        sumhup = costhup
        sumhdo = costhdo
        sumo = costo
        sumden = costd
        sumspin = costs
        sumodo = costodo
        sumhodo = costhodo
#endif

#else

        sumh = 0.d0
        sumhup = 0.d0
        sumhdo = 0.d0
        sumo = 0.d0
        sumhodo = 0.d0
        sumodo = 0.d0
        sumden = costd
        sumspin = costs
        do i = 1, nelorbu
            do j = 1, nelorbu
            if (ipc .eq. 1) then
                sumh = sumh + overham(i, j)
                sumo = sumo + overs(i, j)
                sumhup = sumhup + hamilt(i, j)
                if (yeslsda) sumhdo = sumhdo + hamiltdo(i, j)
            else
                sumh = sumh + overham(2*i - 1, j)
                sumo = sumo + overs(2*i - 1, j)
                sumhodo = sumhodo + overhamdo(2*i - 1, j)
                sumodo = sumodo + oversdo(2*i - 1, j)
                sumhup = sumhup + hamilt(2*i - 1, j)
                sumhdo = sumhdo + hamiltdo(2*i - 1, j)

                sumh = sumh + overham(2*i, j)
                sumo = sumo + overs(2*i, j)
                sumhodo = sumhodo + overhamdo(2*i, j)
                sumodo = sumodo + oversdo(2*i, j)
                sumhup = sumhup + hamilt(2*i, j)
                sumhdo = sumhdo + hamiltdo(2*i, j)
            end if
            end do
        end do
#endif
        !
        ! collect the sums for all k-points and normalize
#ifdef PARALLEL
        if (manyfort10) then
            call reduce_base_real(1, sumh, commcolrep_mpi, -1)
            sumh = sumh/nk
            call reduce_base_real(1, sumo, commcolrep_mpi, -1)
            sumo = sumo/nk
            call reduce_base_real(1, sumhup, commcolrep_mpi, -1)
            sumhup = sumhup/nk
            call reduce_base_real(1, sumden, commcolrep_mpi, -1)
            call reduce_base_real(1, sumhdo, commcolrep_mpi, -1)
            sumhdo = sumhdo/nk
            if (yeslsda) &
                call reduce_base_real(1, sumspin, commcolrep_mpi, -1)
        end if
#endif

        if (rank .eq. 0) then
            write (6, *)
            if (yeslsda .or. ipc .eq. 2) then
                write (6, '(a)') ' Check matrix elements sum:'
                write (6, *) ' overs    ', sumo
                write (6, *) ' overham  ', sumh
                write (6, *) ' oversdo  ', sumodo
                write (6, *) ' overhamdo', sumhodo
                write (6, *) ' hamilt   ', sumhup
                write (6, *) ' hamiltdo ', sumhdo
            else
                write (6, '(a)') ' Check matrix elements sum: '
                write (6, *) ' overs    ', sumo
                write (6, *) ' overham  ', sumh
                write (6, *) ' oversdo  ', sumodo
                write (6, *) ' overhamdo', sumhodo
                write (6, *) ' hamilt   ', sumhup
            end if
            write (6, *)
        end if

        return

    end subroutine checksum

    subroutine diagonalize_hamiltonian(iopt, lworkr)

        implicit none
        integer, intent(in) :: iopt
        real(8), intent(inout) :: lworkr
        integer :: optprint
        logical :: from_scratch
#ifdef __SCALAPACK
        integer :: irow, icol, bdim, bdimdo
#endif
#ifdef PARALLEL
        include "mpif.h"
#endif
        !
        ! if optprint=1 eigenvalues Ham and overlaps
        !
        optprint = 0
        if (iopt .eq. 1) optprint = 1
        !
        ! for SC: compute initial solution from potential
        ! without XC.
        !
        from_scratch = .false.
        if (iopt .eq. 1) from_scratch = .true.
        if (compute_bands) from_scratch = .false.
        !
        ! initialize eigenvectors
        !
        molecorb = 0.d0
        if (yeslsda .or. ipc .eq. 2) molecorbdo = 0.d0
        !
#ifdef __SCALAPACK
        !
        ! initialize correct block dimension
        !
        if (ipc .eq. 1) then
            bdim = size(molecorbl, 1)
            bdimdo = size(molecorbldo, 1)
        else
            bdim = size(molecorbl, 1)/2
            bdimdo = size(molecorbldo, 1)/2
        end if

        if (ipc .eq. 1) then
            if (from_scratch) hamiltl = overhaml
            call eval_hamilt(nelorbu, oversl, hamiltl, molecorbl &
                 &, eigmol, nelorbu, epsover, eps_mach, rank, optprint, lworkr, &
                 iopt, premat, eigmat, info, bands, bdim, mincond)
            ! save new eigenvectors
            call eqmat_scalapack(nelorbu, molecorbl, bdim, descla, molecorb, bands)
        else
            if (from_scratch) hamiltl = overhaml
            call eval_hamilt_complex(nelorbu, oversl, hamiltl, molecorbl, umatl, &
                                     eigmol, nelorbu, epsover, eps_mach, rankrep, optprint, lworkr, &
                                     iopt, premat, eigmat, info, bands, bdim, mincond)
            ! save new eigenvectors
            call eqmat_scalapack_complex(nelorbu, molecorbl, bdim, descla, molecorb, bands)
        end if

#ifdef PARALLEL
        call reduce_base_real(size(molecorb), molecorb, commrep_mpi, -1)
#endif

        if (ipc .eq. 1) then ! real case

            if (yeslsda) then
                if (from_scratch) hamiltldo = overhaml
                call eval_hamilt(nelorbu, oversl, hamiltldo, molecorbldo &
                     &, eigmoldo, nelorbu, epsover, eps_mach, rank, optprint, lworkr, &
                     iopt, premat, eigmat, info, bands, bdimdo, mincond)
                molecorbdo = 0.d0
                ! save new eigenvectors
                call eqmat_scalapack(nelorbu, molecorbldo, bdimdo, descla, molecorbdo, bands)
            end if

        else ! complex case

            if (from_scratch) hamiltldo = overhamldo
            call eval_hamilt_complex(nelorbu, oversldo, hamiltldo, molecorbldo, umatldo, &
                                     eigmoldo, nelorbu, epsover, eps_mach, rankrep, optprint, lworkr, &
                                     iopt, premat, eigmat_down, info, bands, bdimdo, mincond)
            ! save new eigenvectors
            call eqmat_scalapack_complex(nelorbu, molecorbldo, bdimdo, descla, molecorbdo, bands)

        end if

#ifdef PARALLEL
        if (yeslsda .or. ipc .eq. 2) call reduce_base_real(size(molecorbdo), molecorbdo, commrep_mpi, -1)
#endif

#else

        if (ipc .eq. 1) then

            if (from_scratch) hamilt = overham

            call eval_hamilt(nelorbu, overs, hamilt, molecorb &
            &, eigmol, nelorbu, epsover, eps_mach, rank, optprint, &
            lworkr, iopt, premat, eigmat, info, bands, 1, mincond)
            if (yeslsda) then
                if (from_scratch) hamiltdo = overham
                call eval_hamilt(nelorbu, overs, hamiltdo, molecorbdo &
                &, eigmoldo, nelorbu, epsover, eps_mach, rank, optprint, lworkr, &
                iopt, premat, eigmat, info, bands, 1, mincond)
            end if

        else

            if (from_scratch) hamilt = overham
            call eval_hamilt_complex(nelorbu, overs, hamilt, molecorb, umatl &
            &, eigmol, nelorbu, epsover, eps_mach, rankrep, optprint, &
            lworkr, iopt, premat, eigmat, info, bands, ipc, mincond)
            if (from_scratch) hamiltdo = overhamdo
            call eval_hamilt_complex(nelorbu, oversdo, hamiltdo, molecorbdo, umatldo &
            &, eigmoldo, nelorbu, epsover, eps_mach, rankrep, optprint, lworkr, &
            iopt, premat, eigmat_down, info, bands, ipc, mincond)

        end if

#endif
        !
        ! orthogonalize orbitals and evaluate improved eigenvalues
        !
        if (orthodiag .and. iopt .ne. 1) call improvediag

#ifdef PARALLEL
        call mpi_bcast(molecorb, size(molecorb), MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
        if (yeslsda .or. ipc .eq. 2) &
            call mpi_bcast(molecorbdo, size(molecorbdo), MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
#endif

        return
    end subroutine diagonalize_hamiltonian

end module freeelmod_complex

subroutine eqmat_scalapack_complex(n, a, lda, desc, b, lcol)
    !nelorbu,molecorbldo,size(molecorbldo,1)/2,descla,molecorbdo
    use descriptors, only: descla_siz_, lambda_node_, nlar_, nlac_, &
                           ilac_, ilar_, nlax_, la_n_
    implicit none
    integer, intent(in) :: n, lda, desc(descla_siz_), lcol
    complex(8), intent(inout) :: a(lda, *), b(n, *)
    integer :: nr, nc, ir, ic, i, j

#ifdef PARALLEL

    if (desc(lambda_node_) <= 0) then
        return
    end if
    if (n /= desc(la_n_)) &
        call errore(" eqmat_scalapack_complex ", " wrong global dim n ", n)
    if (lda /= desc(nlax_)) &
        call errore(" eqmat_scalapack_complex ", " wrong leading dim lda ", lda)

    nr = desc(nlar_)
    nc = desc(nlac_)
    ir = desc(ilar_)
    ic = desc(ilac_)

    do j = 1, nc
        do i = 1, nr
            if ((j + ic - 1) <= lcol) then
                b((i + ir - 1), (j + ic - 1)) = a(i, j)
            end if
        end do
    end do

#endif

    return
end subroutine eqmat_scalapack_complex

subroutine eqmat_scalapack(n, a, lda, desc, b, lcol)
    !nelorbu,molecorbldo,size(molecorbldo,1)/2,descla,molecorbdo
    use descriptors, only: descla_siz_, lambda_node_, nlar_, nlac_, &
                           ilac_, ilar_, nlax_, la_n_
    implicit none
    integer, intent(in) :: n, lda, desc(descla_siz_), lcol
    real(8), intent(inout) :: a(lda, *), b(n, *)
    integer :: nr, nc, ir, ic, i, j

#ifdef PARALLEL

    if (desc(lambda_node_) <= 0) then
        return
    end if
    if (n /= desc(la_n_)) &
        call errore(" eqmat_scalapack ", " wrong global dim n ", n)
    if (lda /= desc(nlax_)) &
        call errore(" eqmat_scalapack ", " wrong leading dim lda ", lda)

    nr = desc(nlar_)
    nc = desc(nlac_)
    ir = desc(ilar_)
    ic = desc(ilac_)

    do j = 1, nc
        do i = 1, nr
            if ((j + ic - 1) <= lcol) then
                b((i + ir - 1), (j + ic - 1)) = a(i, j)
            end if
        end do
    end do

#endif

    return
end subroutine eqmat_scalapack
