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

subroutine write_output_and_finalize

    use setup
    use freeelmod_complex
    use fourier_module
    use parallel_module, only: collect_from_pools

    implicit none
#ifdef PARALLEL
    include "mpif.h"
#endif
    integer :: i, j, ir, ic
    real(8), dimension(:, :), allocatable :: eigvU, eigvD
    real*8, dimension(:, :), allocatable :: molecorbtimeso, over_part, over_new
    logical is_same_up(3), is_same_down(3)
    !
    ! write scratch files with final quantities
    !
    ! writing the final distributed density/spin in the second record
    ! of scratch files for each processor (for continuation)
    if (writescratch .eq. 0) then
        if (typeopt .eq. 2) then
            dent = dentnew
            if (yeslsda) spint = spintnew
        end if
        rewind (unit_scratch_distributed)
        read (unit_scratch_distributed)
        read (unit_scratch_distributed)
        if (ipc .eq. 2) read (unit_scratch_distributed)
        if (double_mesh) read (unit_scratch_distributed)
        if (yeslsda) then
            write (unit_scratch_distributed) dent, spint
        else
            write (unit_scratch_distributed) dent
        end if
    end if
    !
    ! collects eigenvalues and occupations from pools
    !
    errsav = errdft
    if (.not. compute_bands) then
        call collect_from_pools
        if (rank .eq. 0) call print_calculation_outputs()
    end if
    !
    ! be sure all processors in the pool have consistent wavefunction
    !
#ifdef PARALLEL
    call mpi_bcast(molecorb, size(molecorb), MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
    if (yeslsda .or. ipc .eq. 2) &
        call mpi_bcast(molecorbdo, size(molecorbdo), MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
#endif
    !
    ! use molecorb_old/molecorbdo_old to save the converged eigenvectors
    !
    if (.not. symmagp .and. (opposite_phase .or. ipc .eq. 1) .and. try_translation) then
        if (rank .eq. 0) write (6, *) ' Warning sorting down molecular orbitals '
        !  Here molecorb --> over x molecorb
        allocate (molecorbtimeso(ipc*nelorbu, bands))
        molecorbtimeso = 0.d0
#ifdef __SCALAPACK
        ir = 0
        ic = 0
        if (descla(lambda_node_) > 0) then
            allocate (over_part(ipc*nlax, nlax))
            allocate (over_new(ipc*nlax, nlax))
            ir = descla(ilar_)
            ic = descla(ilac_)
            over_part = 0.d0
            over_new = 0.d0
            do j = 1, descla(nlac_)
                do i = 1, descla(nlar_)
                    if ((j + ic - 1) <= bands) then
                        if (ipc .eq. 2) then
                            over_part(2*i - 1, j) = molecorb(2*(i + ir - 1) - 1, (j + ic - 1))
                            over_part(2*i, j) = molecorb(2*(i + ir - 1), (j + ic - 1))
                        else
                            over_part(i, j) = molecorb(i + ir - 1, (j + ic - 1))
                        end if
                    end if
                end do
            end do
            if (ipc .eq. 2) then
                call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, zone, oversl, 1, 1,&
                & desch, over_part, 1, 1, desch, zzero, over_new, 1, 1, desch)
            else
                call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.d0, oversl, 1, 1,&
                & desch, over_part, 1, 1, desch, 0.d0, over_new, 1, 1, desch)
            end if
! MMM
            do j = 1, descla(nlac_)
                do i = 1, descla(nlar_)
                    if ((j + ic - 1) <= bands) then
                    if (ipc .eq. 2) then
                        molecorbtimeso(2*(i + ir - 1) - 1, (j + ic - 1)) = over_new(2*i - 1, j)
                        molecorbtimeso(2*(i + ir - 1), (j + ic - 1)) = over_new(2*i, j)
                    else
                        molecorbtimeso(i + ir - 1, (j + ic - 1)) = over_new(i, j)
                    end if
                    end if
                end do
            end do
            deallocate (over_part, over_new)
        end if

#ifdef PARALLEL
        call reduce_base_real(ipc*nelorbu*bands, molecorbtimeso, commrep_mpi, -1)
#endif

#else
        if (ipc .eq. 2) then
            call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, overs, nelorbu, molecorb, nelorbu, zzero, molecorbtimeso, nelorbu)
        else
            call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, overs, nelorbu, molecorb, nelorbu, 0.d0, molecorbtimeso, nelorbu)
        end if
#endif

        if (rank .eq. 0) then
            write (6, *) ' Molecular x over '
            do i = 1, bands
                write (6, *) i, sum(abs(molecorbtimeso(:, i)))
            end do
        end if
        if (ipc .eq. 2) then
            call sort_molecular(molecorbdo, molecorbtimeso&
                    &, eigmol, eigmoldo, nelorbu, bands, .false., rank)
        else
            call sortr_molecular(molecorbdo, molecorbtimeso&
                    &, eigmol, eigmoldo, nelorbu, bands, rank)
        end if
        deallocate (molecorbtimeso)
    end if

    if (contracted_on) then

        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbh, nelocc, nelorbu, 1.d0, mu_c, nelorbh, &
                       molecorb, nelorbu, 0.d0, molecorb_old, nelorbh)
            if (yeslsda) &
                call dgemm('N', 'N', nelorbh, neloccdo, nelorbu, 1.d0, mu_c, &
                           nelorbh, molecorbdo, nelorbu, 0.d0, molecorbdo_old, nelorbh)
        else

            call zgemm('N', 'N', nelorbh, nelocc, nelorbu, zone, mu_c, nelorbh, &
                       molecorb, nelorbu, zzero, molecorb_old, nelorbh)
            call zgemm('N', 'N', nelorbh, neloccdo, nelorbu, zone, mu_c, &
                       nelorbh, molecorbdo, nelorbu, zzero, molecorbdo_old, nelorbh)

        end if
    else
        molecorb_old(1:ipc*nelorbu, 1:bands) = molecorb(1:ipc*nelorbu, 1:bands)
        if (yeslsda .or. ipc .eq. 2) molecorbdo_old(1:ipc*nelorbu, 1:bands) = molecorbdo(1:ipc*nelorbu, 1:bands)
    end if

    !   Make some space for future allocations
    deallocate (molecorb)
    if (allocated(molecorbdo)) deallocate (molecorbdo)
    if (memlarge) then
        if (allocated(wf)) deallocate (wf)
    end if
#ifdef PARALLEL
    deallocate (fp, sndbuf, rcvbuf)
#else
    deallocate (fp)
#endif

    !
    ! only for SC update the dup_c vector with the new molecular orbitals
    !
    if (.not. compute_bands) then
        !
        ! NB: the QMC assumes an hermitian Lambda_ij in this case, therefore
        ! if this is not the case (as for molecorbdo != dconj(molecorb)) then
        ! the calculation of the phase is not working.
        !
        if (rank .eq. 0) write (6, *)
        ! this case corresponds to:
        ! symmagp = .true.
        ! yes_hermite = .true.
        if (opposite_phase) then

            !
            ! conjugate the down-spin orbital in the case of a yes_hermite=.true.
            ! for both symmetric and non symmetric AGP

            if (symmagp) then
                !
                if (ipc .eq. 2 .or. yeslsda) then
                    molecorbdo_old = molecorb_old
                    if (ipc .eq. 2) call conjmat(nelorb, bands, molecorbdo_old, nelorb)
                    if (rank .eq. 0) then
                        write (6, '(a)') &
                            ' Warning: eigenvectors of down spin electrons are complex conjugate of the up spin ones! '
                    end if
                end if
            elseif (try_translation) then
                if (rank .eq. 0) write (6, '(a)') ' Warning: sorting down eigenvectors &
                        & to attempt a translation invariant AGP (use random twists) '
            end if

            ! this case corresponds to:
            ! symmagp = .true.
            ! yes_hermite = .false.
            ! It is never used since Hamiltonian does not conserve this symmetry.
        elseif (ipc .eq. 2 .and. symmagp) then
            !
            ! in the case of a non-hermitian and symmetric AGP, the
            ! eigenvectors up are simply equal to the eigenvectors down
            !
            if (rank .eq. 0) write (6, '(a)') ' Warning: eigenvectors of down spin electrons are equal to the up spin ones! '
            molecorbdo_old = molecorb_old

#ifdef DEBUG
            ! this case is only for test
        elseif (ipc .eq. 2 .and. .not. symmagp) then

            if (opposite_phase) then

                if (rank .eq. 0) then
                    write (6, '(a)') ' Warning: eigenvectors of down spin electrons are complex conjugate of the up spin ones! '
                end if
                molecorbdo_old = molecorb_old
                call conjmat(nelorb, bands, molecorbdo_old, nelorb)

            elseif (same_phase) then

                molecorbdo_old = molecorb_old
                if (rank .eq. 0) write (6, '(a)') ' Warning: eigenvectors of down spin electrons are equal to the up spin ones! '

            end if
#endif

        elseif (.not. symmagp) then

            if (rank .eq. 0) write (6, '(a)') ' Warning: non-symmetric AGP. No constraints are imposed on the molecular orbitals. '

        end if
        if (rank .eq. 0) write (6, *)
        !
        ! fill fort.10 with new molecular orbitals and update all variables
        ! necessary to write_fort10()
        !

        call update_fort10
        !
        ! write the final wave function. If k-points sampling, write one
        ! wf for each k-point if the flag yeswrite10 is set to .true.

        !
        if (zero_jas) then
            if (contractionj .ne. 0) then
                jasmat_c = 0.d0
            else
                jasmat = 0.d0
            end if
        end if
        if (newj_onebody(1) .gt. 0 .and. abs(iesdr) .gt. 4) then
            if (abs(niesd) .le. 2) then
                vj(abs(niesd)) = newj_onebody(1)
            else
                do j = 2, niesd
                    vj(j) = newj_onebody(j - 1)
                end do
            end if
        end if

        if (manyfort10 .and. yeswrite10 .and. rankrep .eq. 0) then
            if (rank .eq. 0) write (6, '(a)') ' Writing final wave functions for all k-points!'
            rewind (unit_scratch_fort10)
            if (yeslsda .or. ipc .eq. 2) then
                nelup = nint(sum(occupations(1:bands)))
                neldo = nint(sum(occupationdo(1:bands)))
                nel = nelup + neldo
            else
                neldo = (nint(sum(occupations(1:bands))) - ndiff)/2
                nelup = neldo + ndiff
                nel = nelup + neldo
            end if
            call write_fort10(unit_scratch_fort10)
        end if
        !
        ! write the final wavefunction for the first k-point in fort.10_new
        !
        if (rank .eq. 0) then
            if (.not. manyfort10) then
                write (6, '(a)') ' Write the parameters of the final wavefunction!'
                close (10)
                open (unit=10, file='fort.10_new', form='formatted', status='unknown', position='rewind')
                call write_fort10(ufort10)
            end if
        end if

    end if ! endif .not. compute_bands
#ifdef PARALLEL
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier
    !
    ! deallocate all local variables
    !
#ifdef __SCALAPACK
    deallocate (molecorbl, oversl, overhaml, hamiltl)
    if (yeslsda .or. ipc .eq. 2) then
        deallocate (molecorbldo, hamiltldo)
        if (ipc .eq. 2) &
            deallocate (oversldo, umatldo, overhamldo, umatl)
    end if

    deallocate (molecorb_old, overs, overham, hamilt)
    if (yeslsda .or. ipc .eq. 2) then
        deallocate (molecorbdo_old, hamiltdo)
        if (ipc .eq. 2) deallocate (oversdo, overhamdo)
    end if
    if (allocated(overinvsl)) deallocate (overinvsl)
#else
    if (epssr .ne. 0) deallocate (overinvs, overinvsl)
    deallocate (molecorb_old, overs, overham, hamilt)
    if (yeslsda .or. ipc .eq. 2) then
        deallocate (molecorbdo_old, hamiltdo)
        if (ipc .eq. 2) &
            deallocate (oversdo, overhamdo, umatl, umatldo)
    end if
#endif

    if (typeopt .lt. 0 .or. typeopt .eq. 3) then
        if (yeslsda .or. ipc .eq. 2) deallocate (molecorbdos)
        deallocate (molecorbs)
    end if

    deallocate (nozeroc_in, jbradetc_in, jbradetnc_in)
    if (allocated(nozerojc_in)) deallocate (nozerojc_in)
    if (allocated(jbrajc_in)) deallocate (jbrajc_in)
    if (allocated(jbrajnc_in)) deallocate (jbrajnc_in)

    deallocate (eigmol, newocc, oldocc, occupations, eigmoldo, oldoccdo, occupationdo, &
                occupations_sav, eigmol_sav, occupationsdo_sav, eigmat_down)

    if (yeslsda .or. ipc .eq. 2) deallocate (eigmoldo_sav)
    deallocate (vhartree, dentold, premat, vhartreeq)
    if (yeslsda) deallocate (spintold)
    if (yesgrad) deallocate (gradt)

    if (typeopt .eq. 2) then
        deallocate (dent_after, dent_before, dentnew)
        if (yeslsda) deallocate (spint_after, spint_before, spintnew)
    elseif (typeopt .eq. 4) then
        deallocate (dent_after, dent_before, dent_aftern, mixingtrue)
        if (yeslsda) deallocate (spint_after, spint_before, spint_aftern)
        deallocate (overjac, jac, mataux, vetaux)
        if (jaccond .gt. 0.d0) deallocate (sjac, sojac, ujac, vjac)
    end if
    if (allocated(volmesh_proc)) deallocate (volmesh_proc)
    if (typeopt .eq. 4) then
        deallocate (eigmolo)
        if (yeslsda .or. ipc .eq. 2) deallocate (eigmolodo)
    end if
    if (allocated(spin_input)) deallocate (spin_input)
    if (allocated(charge_input)) deallocate (charge_input)
    if (allocated(gridspin)) deallocate (gridspin, gridnospin)
    if (allocated(gridcharge)) deallocate (gridcharge, gridnocharge)
    if (allocated(ipsip)) deallocate (ipsip)
    if (allocated(psip)) deallocate (psip)
    ! useful only for writing the w.f.
    ! useless in the NonSC case
    if (nmoltot .gt. molecular) then
        iscraipsip = max(iscraipsip, iesupr_c + iesupind)
        molecular = nmoltot ! redefine the number or molecular
    end if

    if (double_mesh) then
        deallocate (nx8, nxny8, ind_init, nx_proc, nxny_proc, ncub_min, ncub_max)
    end if

    allocate (ipsip(iscraipsip))
    allocate (psip(iscramax))
    ipsip = 0
    psip = 0.d0
    !
    ! write quantities to be saved to scratch file "TurboDFT.sav"
    !
    if (writescratch .eq. 0) then
        if (rank .eq. 0) write (6, '(a)') ' Writing scratch files '
        if (write_den) call write_total_density()
        if (rank .eq. 0) then
            allocate (eigvU(ipc*nelorbu, bands))
            eigvU = 0.d0
            if (yeslsda .or. ipc .eq. 2) then
                allocate (eigvD(ipc*nelorbu, bands))
                eigvD = 0.d0
            end if
            do i = 1, nk
                if (yeslsda .or. ipc .eq. 2) then
                    eigvU(1:ipc*nelorbu, 1:bands) = molecorb_sav(1:ipc*nelorbu, 1:bands, i)
                    eigvD(1:ipc*nelorbu, 1:bands) = molecorbdo_sav(1:ipc*nelorbu, 1:bands, i)
                    write (unit_scratch_eigenvects) eigvU, eigvD
                else
                    eigvU(1:ipc*nelorbu, 1:bands) = molecorb_sav(1:ipc*nelorbu, 1:bands, i)
                    write (unit_scratch_eigenvects) eigvU
                end if
            end do
            deallocate (eigvU)
            if (yeslsda .or. ipc .eq. 2) deallocate (eigvD)
            close (unit_scratch_densities)
        end if
    end if
    !
    ! print final density in xcrysden format
    !
    if (write_den) call printden(rion_ref, rion_shift, nxl, dent)

    deallocate (dent)
    if (yeslsda) deallocate (spint)

    if (allocated(buffer_grid)) deallocate (buffer_grid)
    if (allocated(dens_grid)) deallocate (dens_grid)
    if (allocated(spint_grid)) deallocate (spint_grid)
    if (allocated(molecorb_sav)) deallocate (molecorb_sav)
    if (allocated(molecorbdo_sav)) deallocate (molecorbdo_sav)
    !
    ! close all files
    !
    if (rank .eq. 0) then
        close (ufort10)
        close (11)
        if (npsa .ne. 0) close (8) ! pseudopotential
    end if
    !
    ! close scratch files
    !
    !   if(writescratch.eq.0) then
    if (writescratch .eq. 0) close (unit_scratch_distributed)
    if (yeswrite10 .and. rankrep .eq. 0 .and. manyfort10) close (unit_scratch_fort10)
    if (rank .eq. 0) then
        close (unit_scratch_densities)
        close (unit_scratch_eigenvects)
        if (print_energies) close (unit_print_energies)
    end if
    !   endif
    !
    ! deallocate all variables allocated by DFT/QMC routines
    !
    call deallocate_all

    return

contains

    subroutine print_calculation_outputs()

        implicit none

        if (iter .ge. maxit) then
            write (6, 100) ' Warning Turbo-DFT  terminates without convergence, Error= ', real(errsav)
        else
            write (6, 101) ' OK Turbo-DFT converged  with energy tollerance  ', real(errsav), '<', real(epsdft)
            write (6, '(A,X,I5)') ' # Iterations =', iter
        end if
        if (mod(typeopt, 2) .eq. 1 .and. yespassed) then
            write (6, 100) ' Final variational DFT  energy (Ha) = ', edftp
            edftvar = edftp
        elseif (mod(typeopt, 2) .eq. 0) then
            write (6, 100) ' Final variational DFT  energy (Ha) = ', edftvar
        end if
        write (6, 100) ' Final self consistent energy (Ha) =', edft
        if (double_mesh .and. corr_hartree .and. iespbc) then
            write (6, 100) ' Final estimated energy (corrected Hartree) ', edftvar + ehartree - eh_ew
        end if

        write (6, 100) ' Final exchange  energy            =', exchange
        write (6, 100) ' Final correlation  energy         =', ecorr
        write (6, 101) ' Final Fermi energy  =', efermi, 'Ha =', efermi*energy_unit, 'eV'
        if (yeslsda) then
            write (6, 100) ' Final total magnetization (a.u.) ', spingrid/2.d0
        end if

        if (iespbc) then
            write (6, 100) ' Final Hartree energy (QE conv. sum q=/0) (Ha) =', vh_test
            write (6, 100) ' Final Hartree energy on a mesh  =', vh_att
            if (.not. corr_hartree) write (6, 100) ' Final Hartree energy no Ewald corrected (Ha) =', eh_ew
            write (6, *) ' Ewald contribution (QE conv) = ', vpotaa/2.d0
            write (6, 100) ' Eself (Qbox conv) = ', kappanew/dsqrt(pi)*sum(zetar(1:nion)**2)
            write (6, 100) ' E_sr (Qbox conv) = ', &
                    &vpotaa/2.d0 + kappanew/dsqrt(pi)*sum(zetar(1:nion)**2)
        else
            write (6, 100) ' Final Hartree energy(H)  =', ehartree
            if (double_mesh) write (6, 100) ' Final Hartree larger mesh energy (Ha) =', vh_test
        end if
        if (corr_hartree .and. scale_hartree .gt. 0) write (6, 100) ' Final rescaled Hartree density factor =', scale_hartreen
        write (6, *)
        write (6, '(A)') ' Turbo-DFT TIMINGS '
        write (6, 100) ' Total time (sec.)  =', time_total
        write (6, 100) ' Total initialization time (sec.)  =', init_time
        write (6, 100) ' Total loading time matrices (sec.)  =', loading_time
        write (6, 100) ' Total diagonalization  time (sec.) =', diag_time
        write (6, 102) ' Total/Upload FFT  time (sec.) =', time_fft, time_uploadfft
        write (6, 100) ' Total self-consistent cycle time (sec.)  =', cycle_time
        write (6, 100) ' Total density symmetrization time (sec.) =', symtime
        if (ipc .eq. 1) then
            write (6, 100) ' Total time dgemm =', dgemm_time
        else
            write (6, 100) ' Total time zgemm =', zgemm_time
        end if
        if (diag_time .gt. loading_time) then
            write (6, *) ' Warning you should run with a smaller number of processors !!! '
#ifdef __DOOMP
#else
            write (6, *) ' Remind the diagonalization part is not OpenMP parallel so far '
#endif
        end if

        return

100     format(A, X, F28.15)
101     format(A, X, F28.15, X, A, F28.15, X, A)
102     format(A, X, 2f28.15)

    end subroutine print_calculation_outputs

end subroutine write_output_and_finalize

subroutine write_total_density

    use allio, only: rank, commrep_mpi, nprocrep, nx, ny, nz, writescratch
    use setup, only: spint, dent, buffer_grid, meshproc, yeslsda, spint_grid, dens_grid, &
                     unit_scratch_densities

    implicit none
    integer :: proc, ind, i, j, k, ierr, indproc
#ifdef PARALLEL
    include "mpif.h"
#endif
    !
    ! save the total density on a file for NonSC calculations
    !
    ! fill total charge density
    allocate (dens_grid(nx, ny, nz))
    dens_grid = 0.d0
    if (yeslsda) then
        allocate (spint_grid(nx, ny, nz))
        spint_grid = 0.d0
    end if
#ifdef PARALLEL
    buffer_grid = 0.0d0
    call mpi_gather(dent, meshproc, MPI_DOUBLE_PRECISION, buffer_grid&
         &, meshproc, MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
    if (rank .eq. 0) then
        indproc = 0
        do proc = 0, nprocrep - 1
            ind = 0
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        if (indproc .eq. proc) then
                            ind = ind + 1
                            dens_grid(i, j, k) = buffer_grid(ind, proc + 1)
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                    end do
                end do
            end do
        end do
    end if
    ! fill total spin density
    if (yeslsda) then
        buffer_grid = 0.d0
        call mpi_gather(spint, meshproc, MPI_DOUBLE_PRECISION, buffer_grid&
             &, meshproc, MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
        if (rank .eq. 0) then
            spint_grid = 0.d0
            indproc = 0
            do proc = 0, nprocrep - 1
                ind = 0
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            if (proc .eq. indproc) then
                                ind = ind + 1
                                spint_grid(i, j, k) = buffer_grid(ind, proc + 1)
                            end if
                            indproc = indproc + 1
                            if (indproc .eq. nprocrep) indproc = 0
                        end do
                    end do
                end do
            end do
        end if
    end if
#else
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                dens_grid(i, j, k) = dent(i)
                if (yeslsda) spint_grid(i, j, k) = spint(i)
            end do
        end do
    end do
#endif
    ! write total density on scratch files
    if (rank .eq. 0 .and. writescratch .eq. 0) then
        if (yeslsda) then
            write (unit_scratch_densities) dens_grid
            write (unit_scratch_densities) spint_grid
        else
            write (unit_scratch_densities) dens_grid
        end if
        close (unit_scratch_densities)
    end if

    deallocate (dens_grid)
    if (allocated(spint_grid)) deallocate (spint_grid)

    return
end subroutine write_total_density

!
! sort down spin molecular orbitals according to overlap (real version)
!
subroutine sortr_molecular(moleculardo, molecularup, eigup, eigdo, nelorbh, bands, rank)

    implicit none
    !
    ! input
    integer, intent(in) :: nelorbh, bands, rank
    real(8), intent(inout) :: moleculardo(nelorbh, bands), molecularup(nelorbh, bands), eigup(bands), eigdo(bands)
    !
    ! local
    integer :: i, j, j_min
    logical :: conjugate
    real(8) :: distmax, dist, cost, over
    logical, dimension(:), allocatable :: done
    real(8), dimension(:, :), allocatable :: molecular_sav
    real(8) :: phase_min, overc

    allocate (molecular_sav(nelorbh, bands))
    allocate (done(bands))
    done = .false.
    distmax = max(eigup(bands), eigdo(bands)) - min(eigup(1), eigdo(1)) + 1
    molecular_sav = moleculardo
    if (rank .eq. 0) write (6, *) ' Warning sorting molecular according to momentum conservation in AGP', distmax
    do i = 1, bands
        dist = distmax
        j_min = 0
        do j = 1, bands
            overc = sum(molecularup(:, i)*molecular_sav(:, j))&
                    & /sqrt(sum(molecularup(:, i)**2)*sum(molecular_sav(:, j)**2))
            over = abs(overc)

            !       cost=abs(eigup(i)-eigdo(j))-over
            cost = abs((eigup(i) - eigdo(j))*(1.d0 - over))
            if (.not. done(j) .and. cost .lt. dist) then
                dist = cost
                j_min = j
                phase_min = overc/over
            end if
        end do
        if (j_min .ne. 0) then
            done(j_min) = .true.
            moleculardo(:, i) = phase_min*molecular_sav(:, j_min)
            if (rank .eq. 0 .and. i .ne. j_min) write (6, *) ' Warning changed order ', i, j_min
        else
            if (rank .eq. 0) write (6, *) ' ERROR not found molecular orbital #', i
        end if
    end do
    deallocate (done, molecular_sav)
    return

end subroutine sortr_molecular

!
! sort down spin molecular orbitals according to overlap (complex version)
!
subroutine sort_molecular(moleculardo, molecularup, eigup, eigdo, nelorbh, bands, conjugate, rank)

    implicit none
    !
    ! input
    integer, intent(in) :: nelorbh, bands, rank
    logical, intent(in) :: conjugate
    complex(8), intent(in) :: molecularup(nelorbh, bands)
    complex(8), intent(inout) :: moleculardo(nelorbh, bands)
    real(8), intent(in) :: eigup(bands), eigdo(bands)
    !
    ! local
    integer :: i, j, j_min
    logical, dimension(:), allocatable :: done
    real(8) :: distmax, dist, cost, over
    complex(8), dimension(:, :), allocatable :: molecular_sav
    complex(8) :: overc, phase_min

    allocate (molecular_sav(nelorbh, bands))
    allocate (done(bands))

    done = .false.
    distmax = max(eigup(bands), eigdo(bands)) - min(eigup(1), eigdo(1)) + 1
    molecular_sav = moleculardo

    if (rank .eq. 0) write (6, *) ' Warning sorting molecular according to momentum conservation in AGP', distmax
    do i = 1, bands
        dist = distmax
        j_min = 0
        do j = 1, bands
            if (conjugate) then
                overc = sum(molecularup(:, i)*conjg(molecular_sav(:, j)))&
                        & /sqrt(sum(abs(molecularup(:, i))**2)*sum(abs(molecular_sav(:, j))**2))
            else
                overc = sum(molecularup(:, i)*molecular_sav(:, j))&
                        & /sqrt(sum(abs(molecularup(:, i))**2)*sum(abs(molecular_sav(:, j))**2))
            end if
            over = abs(overc)
            !          cost=abs(eigup(i)-eigdo(j))-over
            cost = abs((eigup(i) - eigdo(j))*(1.d0 - over))
            if (.not. done(j) .and. cost .lt. dist) then
                dist = cost
                if (conjugate) then
                    phase_min = overc/over
                else
                    phase_min = conjg(overc/over)
                end if
                j_min = j
            end if
        end do

        if (j_min .ne. 0) then
            done(j_min) = .true.

            moleculardo(:, i) = phase_min*molecular_sav(:, j_min)
            if (rank .eq. 0 .and. i .ne. j_min) write (6, *) ' Warning changed order ', i, j_min
        else
            if (rank .eq. 0) write (6, *) ' ERROR not found molecular orbital #', i
        end if
    end do
    deallocate (done, molecular_sav)
    return

end subroutine sort_molecular
