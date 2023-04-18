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

!-----------------------------
subroutine non_self_consistent_run()
    !-----------------------------

    use constants, only: TWO_PI, qp
    use kpoints_mod
    use setup
    use freeelmod_complex, only: diagonalize_hamiltonian, checksum
    use fourier_module
    implicit none

    integer :: indkp, mpi_err, i
    logical :: occs_file_e
    real(8) :: lworkr
#ifdef PARALLEL
    include "mpif.h"
#endif
    !
    ! initialize vhartreeq in reciprocal space
    ! Ewald correction to Hartree potential.
    !
    call evalvhartreeq()
    if (rank .eq. 0) then
        write (6, '(a,2F15.8)') ' Initial Hartree potential: ', sum(vhartree(:)), sum(vhartreeq(:))
        if (iespbc) then
            write (6, '(a,F15.8)') ' Kappa found:         ', kappa
            write (6, '(a,F15.8)') ' Ewald constant found:', ewaldion1b + (1.d0 - weightvh)*ewaldel1b
            write (6, '(a,F15.8)') ' eselfion =           ', ewaldion1b
            write (6, '(a,F15.8)') ' eself1bel =          ', ewaldel1b
        end if
    end if
    !
    ! read electronic density from self-consistent run
    !
    call read_density_from_file()

    ! just to be sure initialize all to zero (again)
    molecorb = 0.d0
    if (yeslsda .or. ipc .eq. 2) molecorbdo = 0.d0
    eigmol_sav = 0.d0
    molecorb_sav = 0.d0
    if (yeslsda .or. ipc .eq. 2) then
        eigmoldo_sav = 0.d0
        molecorbdo_sav = 0.d0
    end if

    k_loop: do indkp = 1, nk
        !
        ! initialize k-point
        !
        phase(:) = xkp(:, indkp)
        phase_down(:) = xkp_down(:, indkp)
        rphase(:) = 0.d0
        phase2pi(:) = phase(:)*TWO_PI
        phase2pi_down(:) = phase_down(:)*TWO_PI
        indk = indkp
        !
        ! compute basis and Hamiltonian overlaps:
        ! S_ij = < \phi_i | \phi_j >
        ! H_ij = < \phi_i | H | \phi_j >
        ! With the new basis set (yes_crystal=.true.) they are both dependent
        ! from the chosen k-point.
        !
        call initialize_mats_new()
        !
        ! update density-dependent part of the hamiltonian
        ! starting from the density read from file
        !
        call uphamilt_new()
        !
        ! diagonalize new hamiltonian
        !
        lworkr = 0
        call diagonalize_hamiltonian(1, lworkr)
        if (orthodiag) call improvediag()
        !
        ! check matrix elements sum in debug mode
        !
#ifdef DEBUG
        call checksum()
#endif

        if (rank .eq. 0) then
            write (6, *)
            write (6, '(A)') ' Eigenvalues up/down: '
            do i = 1, bands
                write (6, '(I2,2X,2F10.6)') i, eigmol(i), eigmoldo(i)
            end do
        end if

        ! save eigenvalues/eigenvectors
        ! for postprocessing tools
        eigmol_sav(1:bands, indkp) = eigmol(1:bands)
        molecorb_sav(:, 1:bands, indkp) = molecorb(:, 1:bands)
        occupations_sav(1:bands, indkp) = occupations(1:bands)
        if (yeslsda .or. ipc .eq. 2) then
            eigmoldo_sav(1:bands, indkp) = eigmoldo(1:bands)
            molecorbdo_sav(:, 1:bands, indkp) = molecorbdo(:, 1:bands)
            occupationsdo_sav(1:bands, indkp) = occupationdo(1:bands)
        end if
        if (rank .eq. 0) write (6, '(A,2X,I4,A,3F10.6)') ' k-point # ', indkp, ' completed: ', xkp(:, indkp)
        ! put to zero variables related to hamiltonian
        ! which will be updated next iteration
        eigmol = 0.d0
        molecorb = 0.d0
        if (yeslsda .or. ipc .eq. 2) then
            eigmoldo = 0.d0
            molecorbdo = 0.d0
        end if
#ifdef __SCALAPACK
        molecorbl = 0.d0
        if (yeslsda .or. ipc .eq. 2) molecorbldo = 0.d0
#endif

    end do k_loop
    !
    ! compute relevant quantities such as band structure
    ! or DOS. The type of computation is selected with the
    ! input variable "task":
    !                       task=1 --> band structure calculation
    !                       task=2 --> Density of states
    !                       task=3 --> Fermi surface calculation (to be implemented)
    !
    call postprocessing(task)

    return

end subroutine non_self_consistent_run

!----------------------------------
subroutine read_density_from_file()
    !----------------------------------

    use constants, only: ipc
    use allio, only: rank, nproc, nk, nprocrep, nx, ny, nz
    use setup, only: nproco, nko, yeslsda, spint, dent, unit_scratch_distributed, meshproc, unit_scratch_densities
    implicit none

    integer :: i, j, k, proc, indmesh, mpi_err, ierr
    real(8), dimension(:, :, :), allocatable :: dens_grid, spint_grid
    real(8), dimension(:, :), allocatable :: buffer_grid
#ifdef PARALLEL
    include "mpif.h"
#endif

    if (nproc/nk .ne. nproco/nko) then
        call error(' non_self_consistent_run ', ' Number of processor per k-point &
                &different with respect to self-consistent run ', -1, rank)
    end if

    if (nproc/nk .eq. nproco/nko) then
        ! if same number of processors (per k-point) as SC, then read
        ! usual scratch files "tmp0xxxxx"
        rewind (unit_scratch_distributed)
        ! first line(s) contains the overlap matrices
        read (unit_scratch_distributed)
        if (ipc .eq. 2) read (unit_scratch_distributed)
        if (yeslsda) then
            read (unit_scratch_distributed) dent, spint
        else
            read (unit_scratch_distributed) dent
        end if
        close (unit_scratch_distributed)
    else
        if (rank .eq. 0) write (6, *) ' ERROR run with the same number of proc.! '
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif
        stop

    end if

    !
    ! check if initial total charge read from file is OK
    call cutdens()

    return
end subroutine read_density_from_file

!----------------------------------
subroutine postprocessing(task)
    !----------------------------------
    use allio, only: rank
    implicit none
    integer, intent(in) :: task

    select case (task)

    case (1) ! band structure calculation
        if (rank .eq. 0) call plot_bands()

    case (2) ! DOS calculation
        if (rank .eq. 0) call evaluate_dos()

        !case(3) ! Fermi surface calculation (TO BE IMPLEMENTED)
        !if(rank.eq.0) call evaluate_fermisurf

    case (5) ! evaluate all
        if (rank .eq. 0) then
            call plot_bands()
            call evaluate_dos()
        end if

    case (0)
        call error(' non_sc_postprocessing ', ' No tasks required ', -1, rank)

    case default
        call error(' non_sc_postprocessing ', ' Task type not recognized ', -1, rank)

    end select

    return

end subroutine postprocessing

!---------------------
subroutine plot_bands
    !---------------------

    ! Routine for plotting band structures from NSC calculations. Adapted from the
    ! QuantumESPRESSO code.
    ! kp_type must be 2 or 3 in order to run this routine.
    ! In the case of a kpoints path (kp_type=3), one already knows the number of k-points and the
    ! number of lines in the path, which are "nk" and "nk1-1" respectively.

    use allio, only: rank
    use kpoints_mod, only: nk, xkp, kp_type
    use setup, only: eigmol_sav, eigmoldo_sav, yeslsda, bands
    use constants, only: ipc, energy_unit ! Ha to eV conversion

    implicit none
    logical, dimension(:), allocatable :: high_symmetry, in_range, in_rangedo
    integer i, j, k
    real(8), dimension(:), allocatable :: kp_path
    real(8) tmp_norm, k1(3), k2(3), dkl, dkl_sav, emin, emax
    character(len=80) :: filename

    if (kp_type .lt. 2 .or. kp_type .gt. 3) then
        call error(' postprocessing ', ' Use kp_type=2 or kp_type=3 for plotting &
                &              the DFT band structure. ', -1, rank)
        return
    end if

    ! initialization
    write (6, *) '---------------------------'
    write (6, *) ' Evaluating band structure '
    write (6, *) '---------------------------'

    filename = 'band_structure.dat'
    allocate (high_symmetry(nk), in_range(bands))
    in_range = .true.
    if (yeslsda .or. ipc .eq. 2) then
        allocate (in_rangedo(bands))
        in_rangedo = .true.
    end if
    allocate (kp_path(nk))
    high_symmetry = .false.
    kp_path(1) = 0.d0
    !
    ! find high symmetry points in the 2D plot
    !
    do k = 1, nk
        if (k .eq. 1 .or. k .eq. nk) then
            high_symmetry(k) = .true.
        else
            k1(:) = xkp(:, k) - xkp(:, k - 1)
            k2(:) = xkp(:, k + 1) - xkp(:, k)
            tmp_norm = (k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3))/ &
                       dsqrt(k1(1)*k1(1) + k1(2)*k1(2) + k1(3)*k1(3))/ &
                       dsqrt(k2(1)*k2(1) + k2(2)*k2(2) + k2(3)*k2(3))
            ! if changing direction in the 1BZ, then it's a high symmetry point
            high_symmetry(k) = (abs(tmp_norm - 1.d0) .gt. 1.0d-4)
        end if
        if (sum(xkp(:, k)**2) .lt. 1.0d-9) high_symmetry(k) = .true.
        ! distance between 2 kps.
        if (k .eq. 2) dkl_sav = dsqrt(k1(1)**2 + k1(2)**2 + k1(3)**2)
    end do
    !
    ! create 2D path
    !
    do k = 2, nk

        dkl = dsqrt((xkp(1, k) - xkp(1, k - 1))**2 + &
                    (xkp(2, k) - xkp(2, k - 1))**2 + &
                    (xkp(3, k) - xkp(3, k - 1))**2)

        if (dkl .gt. 5.d0*dkl_sav) then
            kp_path(k) = kp_path(k - 1)
        elseif (dkl .gt. 1.0d-5) then
            kp_path(k) = kp_path(k - 1) + dkl
            dkl_sav = dkl
        else
            kp_path(k) = kp_path(k - 1) + dkl
        end if

    end do
    ! find minimum and maximum values
    emin = eigmol_sav(1, 1)
    emax = eigmol_sav(1, 1)
    do k = 1, nk
        do j = 1, bands
            emin = min(emin, eigmol_sav(j, k))
            if (yeslsda .or. ipc .eq. 2) emin = min(emin, eigmoldo_sav(j, k))
            emax = max(emax, eigmol_sav(j, k))
            if (yeslsda .or. ipc .eq. 2) emax = max(emax, eigmoldo_sav(j, k))
        end do
    end do
    ! find bands values which are in the specified energy range
    in_range = .true.
    in_rangedo = .true.
    do i = 1, bands
        in_range(i) = any(eigmol_sav(i, 1:nk) >= emin .and. eigmol_sav(i, 1:nk) <= emax)
        in_rangedo(i) = any(eigmoldo_sav(i, 1:nk) >= emin .and. eigmoldo_sav(i, 1:nk) <= emax)
    end do

    write (6, *) ' Band structures INFORMATION: '
    write (6, *) 'Minimum/maximum bands values:', emin, emax
    write (6, *) 'High symmetry points: '
    do k = 1, nk
        if (high_symmetry(k)) then
            write (6, '(F10.6,3F10.6)') kp_path(k), xkp(:, k)
        end if
    end do
    !
    ! now draw bands
    open (unit=80, file=trim(filename), form='formatted', status='unknown')
    write (80, *) ' # Spin up DFT bands (eV) '
    do i = 1, bands
        if (in_range(i)) then
            if (mod(i, 2) /= 0) then
                write (80, 100) (kp_path(j), energy_unit*eigmol_sav(i, j), j=1, nk)
            else
                write (80, 100) (kp_path(j), energy_unit*eigmol_sav(i, j), j=nk, 1, -1)
            end if
        end if
    end do
    if (yeslsda .or. ipc .eq. 2) then
        ! add two carriage return for gnuplot
        write (80, *)
        write (80, *)
        write (80, *) ' # Spin down DFT bands (eV)'
        do i = 1, bands
            if (in_rangedo(i)) then
                if (mod(i, 2) /= 0) then
                    write (80, 100) (kp_path(j), energy_unit*eigmoldo_sav(i, j), j=1, nk)
                else
                    write (80, 100) (kp_path(j), energy_unit*eigmoldo_sav(i, j), j=nk, 1, -1)
                end if
            end if
        end do
    end if
    !
    close (80)
    deallocate (high_symmetry, kp_path, in_range)
    if (yeslsda .or. ipc .eq. 2) deallocate (in_rangedo)

100 format(2f10.4)

    return

end subroutine plot_bands

!-----------------------
subroutine evaluate_dos
    !-----------------------

    !
    ! Routine to calculate the density of states.
    ! Use optocc>0 and a non-zero value of the smearing (epsshell)
    ! in order to compute the density of states.
    !

    use constants, only: ipc, energy_unit, qp
    use kpoints_mod, only: nk, xkp, wkp
    use setup, only: eigmol_sav, eigmoldo_sav, epsshell, optocc, deltaE, &
                     emin, emax, yeslsda, bands
    use compute_efermi, only: smearD

    implicit none
    integer :: i, j, k, ibnd, indk, smear_type
    integer :: ndos ! # of points in the DOS
    real(8), dimension(:), allocatable :: tdos, dos_up, dos_down
    real(qp) :: e, dt
    character(len=80) :: filename

    ! initialization
    write (6, *) '------------------------------'
    write (6, *) ' Evaluating density of states '
    write (6, *) '------------------------------'
    filename = 'density_of_states.dat'
    if (optocc .eq. 0 .or. epsshell .eq. 0.d0) then
        write (6, *) 'Warning: DOS must be updated with a smearing function! Setting default values.'
        optocc = 1
        epsshell = 0.01d0
    end if
    smear_type = optocc
    dt = epsshell
    ! setting min/max values of the energies
    ! if not chosen from the user
    if (emin .eq. 0.d0 .and. emax .eq. 0.d0) then
        emin = eigmol_sav(1, 1)
        emax = eigmol_sav(1, 1)
        do k = 1, nk
            do j = 1, bands
                emin = min(emin, eigmol_sav(j, k))
                if (yeslsda .or. ipc .eq. 2) emin = min(emin, eigmoldo_sav(j, k))
                emax = max(emax, eigmol_sav(j, k))
                if (yeslsda .or. ipc .eq. 2) emax = max(emax, eigmoldo_sav(j, k))
            end do
        end do
    end if
    !
    ndos = int((emax - emin)/deltaE)
    allocate (dos_up(ndos + 1), dos_down(ndos + 1), tdos(ndos + 1))

    dos_up = 0.d0
    dos_down = 0.d0
    tdos = 0.d0
    e = 0.d0

    open (unit=80, file=trim(filename), form='formatted', status='unknown')
    write (80, *) '#  energy (eV)       Density_of_states (spin up)'
    do i = 0, ndos
        e = emin + deltaE*i
        do indk = 1, nk
            do ibnd = 1, bands
                dos_up(i + 1) = dos_up(i + 1) + wkp(indk)*smearD(eigmol_sav(ibnd, nk), dt, e, smear_type)
            end do
        end do
        write (80, 100) e*energy_unit, dos_up(i + 1)
    end do

    if (yeslsda .or. ipc .eq. 2) then
        write (80, *)
        write (80, *)
        write (80, *) '#  energy (eV)       Density_of_states (spin down)'
        do i = 0, ndos
            e = emin + deltaE*i
            do indk = 1, nk
                do ibnd = 1, bands
                    dos_down(i + 1) = dos_down(i + 1) + wkp(indk)*smearD(eigmoldo_sav(ibnd, indk), dt, e, smear_type)
                end do
            end do
            write (80, 100) e*energy_unit, dos_down(i + 1)
        end do
    end if

    deallocate (dos_up, dos_down, tdos)
    write (6, *)
    close (80)

100 format(2f10.4)

    return

end subroutine evaluate_dos
