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

subroutine initialize_environment()

    use setup
    use parallel_module
    use compute_efermi
    use fourier_module
    use freeelmod_complex
    use dielectric
    implicit none

    integer :: i, j, k
    real(8), external :: upvpotaa
    real*8 costkappa
#ifdef __SCALAPACK
    integer :: ic, ir
#endif

#ifdef PARALLEL
    include "mpif.h"
#endif
    !
    ! several DFT calculations are possible
    !
    ! The different options are determined by the following input flags:
    ! kaverage (&kpoints card)         --> if .true. it activates k-point sampling algorithm within DFT
    ! decoupled_run (&parameters card) --> if .true. and k-average is also .true. it performs decoupled calculations for each
    !                                      k-point, i.e. each k-point is computed separately from the others.
    ! In order to perform band structure calculations, perform at first a self-consistent calculation converged with k-points
    ! (kaverage=.true. and decoupled_run=.false.), then save the density and rerun a non self-consistent calculation
    ! (kaverage=.true. and decoupled_run=.true.) starting from the previously computed density on a much finer k-point path.
    !
    ! use the following function to printout the type of calculation performed in the current run
    call print_header(kaverage, decoupled_run, compute_bands)
    occopen = .false.
    bandso = 0
    nproco = 0
    meshproco = 0
    nko = 0
    if (rank .eq. 0) open (unit=11, file='occupationlevels.dat', form='formatted', status='unknown')

    if (rank .eq. 0 .and. (iopt .ne. 1 .or. compute_bands)) then
        occopen = .false.
        read (11, *, err=111, end=111)
        read (11, *, err=111, end=111) bandso, nproco, meshproco, nko
        if (corr_hartree .and. scale_hartree .gt. 0.d0) then
            read (11, *) scale_hartreen
            write (6, *) ' Warning read value of scale_hartree = ', scale_hartreen
        end if
        occopen = .true.
111     continue
    end if
    if (print_energies .and. rank .eq. 0) then
        open (unit=unit_print_energies, file="EDFT_vsk.dat", form='formatted', status='unknown')
    end if

#ifdef PARALLEL
    call mpi_bcast(occopen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(bandso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nproco, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(meshproco, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nko, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (corr_hartree .and. scale_hartree .gt. 0.d0 .and. iopt .ne. 1) then
        call mpi_bcast(scale_hartreen, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end if
#endif
    !
    ! initial checks
    !
    if (occopen) then
        if (rank .eq. 0) write (6, *) ' File occupationlevels.dat OK '
    elseif (iopt .ne. 1) then
        if (rank .eq. 0) write (6, *) ' Warning occupationlevels.dat empty '
    end if
    !
    ! initialize mesh for MOs integration
    !
    call initialize_mesh(nx, ny, nz, ax, ay, az, fx, fy, fz, mesh, i8cost, volmesh&
            &, meshproc, meshproc_tot)
    if (rank .eq. 0) write (6, *) ' after initialize_mesh '
    voltot = volmesh*nx
    voltot = voltot*ny
    voltot = voltot*nz

    if (double_mesh) then
        if (from_ions) then
            voltot_double = 0.d0
            do i = 1, nion
                if (zetar(i) .ge. minz_at) then
                    voltot_double = voltot_double + dble(nx0 - 1)*dble(ny0 - 1)*dble(nz0 - 1)*volmesh
                end if
            end do
        else
            voltot_double = dble(nx0 - 1)*dble(ny0 - 1)*dble(nz0 - 1)*volmesh
        end if
    end if

    ! uniform densities
    dens0 = nel/voltot
    spin0 = (nelup - neldo)/voltot
    !
    ! some checks for non self-consistent calculations
    !
    if (compute_bands) then
        if (.not. occopen) then ! need to start from a self-consistent calculation
            call error(' initialize_calculations ', ' occupationlevels.dat not found! Run a &
                    &self-consistent calculation before computing band        &
                    &structure ', -1, rank)
        end if
        if (meshproco*(nproco/nko) .ne. mesh) then ! need for same mesh
            call error(' initialize_calculations ', ' Use the same mesh &
                    &as SC calculation! ', 1, rank)
        end if
    end if
    !
    ! set indices of the main matrices and k-points:
    ! nelorbu,nelorbudo
    ! nelorb3
    ! indk
    ! double_overs/same_phase
    call set_addresses()
    !
    ! --------------------------------------
    ! General ALLOCATION and INITIALIZATION
    ! --------------------------------------
    !
    ! allocation main matrices (overlaps,hamiltonian,eigenvectors)
    !
#ifdef __SCALAPACK
    !
    ! distributed matrices: allocate the block with maximum local dimension "nlax"
    !
    if (descla(lambda_node_) > 0) then ! process holds a block
        if (.not. allocated(molecorbl)) then
            allocate (molecorbl(ipc*nlax, nlax), oversl(ipc*nlax, nlax), &
                      overhaml(ipc*nlax, nlax), hamiltl(ipc*nlax, nlax))
            molecorbl = 0.d0
            overhaml = 0.d0
            oversl = 0.d0
            hamiltl = 0.d0
            if (yeslsda .or. ipc .eq. 2) then
                allocate (molecorbldo(ipc*nlax, nlax), hamiltldo(ipc*nlax, nlax))
                molecorbldo = 0.d0
                hamiltldo = 0.d0
                if (ipc .eq. 2) then
                    allocate (oversldo(ipc*nlax, nlax), umatldo(ipc*nlax, nlax), &
                              overhamldo(ipc*nlax, nlax), umatl(ipc*nlax, nlax))
                    oversldo = 0.d0
                    umatldo = 0.d0
                    overhamldo = 0.d0
                    umatl = 0.d0
                end if
            end if
        end if
    else ! the process doesn't hold a block
        if (.not. allocated(molecorbl)) then
            allocate (molecorbl(ipc, 1), oversl(ipc, 1), overhaml(ipc, 1), hamiltl(ipc, 1))
            molecorbl = 0.d0
            overhaml = 0.d0
            oversl = 0.d0
            hamiltl = 0.d0

            if (yeslsda .or. ipc .eq. 2) then
                allocate (molecorbldo(ipc, 1), hamiltldo(ipc, 1))
                molecorbldo = 0.d0
                hamiltldo = 0.d0
                if (ipc .eq. 2) then
                    allocate (oversldo(ipc, 1), umatldo(ipc, 1), &
                              overhamldo(ipc, 1), umatl(ipc, 1))
                    oversldo = 0.d0
                    umatldo = 0.d0
                    overhamldo = 0.d0
                    umatl = 0.d0
                end if
            end if
        end if
    end if
    !
    ! non-distributed matrices needed for SCALAPACK
    !
    allocate (molecorb_old(ipc*nelorb, bands), molecorb(ipc*nelorbu, bands), &
              overs(ipc, 1), overham(ipc, 1), hamilt(ipc, 1))
    molecorb_old = 0.d0
    molecorb = 0.d0
    overham = 0.d0
    hamilt = 0.d0
    if (yeslsda .or. ipc .eq. 2) then
        allocate (molecorbdo(ipc*nelorbu, bands), molecorbdo_old(ipc*nelorb, bands), hamiltdo(ipc, 1))
        molecorbdo = 0.d0
        molecorbdo_old = 0.d0
        hamiltdo = 0.d0
        if (ipc .eq. 2) then
            allocate (oversdo(ipc, 1), overhamdo(ipc, 1))
            oversdo = 0.d0
            overhamdo = 0.d0
        end if
    end if
    ! defininig again local first index for
    ! rows and columns for the local block
    ir = descla(ilar_)
    ic = descla(ilac_)
#else
    !
    ! main matrices of the code for the normal version
    !
    allocate (molecorb_old(ipc*nelorb, nelorb), molecorb(ipc*nelorbu, nelorbu), overs(ipc*nelorbu, nelorbu), &
              overham(ipc*nelorbu, nelorbu), hamilt(ipc*nelorbu, nelorbu))
    molecorb_old = 0.d0
    molecorb = 0.d0
    overham = 0.d0
    hamilt = 0.d0
    if (yeslsda .or. ipc .eq. 2) then
        allocate (molecorbdo(ipc*nelorbu, nelorbu), molecorbdo_old(ipc*nelorb, nelorb), &
                  hamiltdo(ipc*nelorbu, nelorbu))
        molecorbdo = 0.d0
        molecorbdo_old = 0.d0
        hamiltdo = 0.d0
        if (ipc .eq. 2) then
            allocate (oversdo(ipc*nelorbu, nelorbu), overhamdo(ipc*nelorbu, nelorbu), &
                      umatl(ipc*nelorbu, nelorbu), umatldo(ipc*nelorbu, nelorbu))
            oversdo = 0.d0
            overhamdo = 0.d0
            umatl = 0.d0
            umatldo = 0.d0
        end if
    end if
#endif
    !
    ! eigenvalues and occupations
    !
    allocate (eigmol(nelorbu), newocc(bands + nelup - neldo), oldocc(bands), occupations(bands))
    eigmol = 0
    newocc = 0
    oldocc = 0
    occupations = 0
    allocate (eigmoldo(nelorbudo), oldoccdo(bandsdo), occupationdo(bandsdo))
    allocate (occupations_sav(bands, nk), eigmol_sav(bands, nk))

    if (writescratch .eq. 0) then
        allocate (molecorb_sav(ipc*nelorbu, bands, nk))
        molecorb_sav = 0.d0
        if (yeslsda .or. ipc .eq. 2) then
            allocate (molecorbdo_sav(ipc*nelorbudo, bands, nk))
            molecorbdo_sav = 0.d0
        end if
    end if
    allocate (occupationsdo_sav(bandsdo, nk), eigmoldo_sav(bandsdo, nk))
    allocate (eigmat(nelorbu), eigmat_down(nelorbu))

    eigmoldo = 0
    oldoccdo = 0
    occupationdo = 0
    occupations_sav = 0
    eigmol_sav = 0.d0
    occupationsdo_sav = 0.d0
    eigmoldo_sav = 0.d0
    eigmat = 0.d0
    eigmat_down = 0.d0
    !
    if (rank .eq. 0) write (6, *) ' after initialize molecorb ', write_den
    ! charge and spin densities
    !
    allocate (dent(meshproc), dentold(meshproc), premat(nelorbu))
    dent = 0.d0
    dentold = 0.d0
    premat = 0.d0
    if (double_mesh) then
        allocate (volmesh_proc(meshproc))
        volmesh_proc = volmesh
    end if
    if (yeslsda) then
        allocate (spint(meshproc), spintold(meshproc))
        spint = 0.d0
        spintold = 0.d0
    end if
    if (yeslsda .and. (iopt .eq. 1 .or. randspin .lt. 0.d0 .or. h_field .ne. 0)) then
        allocate (gridspin(meshproc), gridnospin(meshproc))
        gridspin = .true.
        gridnospin = .true.
    end if
    if (h_charge .ne. 0.d0) then
        allocate (gridcharge(meshproc), gridnocharge(meshproc))
        gridcharge = .true.
        gridnocharge = .true.
    end if
    if (write_den) then
        allocate (buffer_grid(meshproc, nprocrep))
        buffer_grid = 0.d0
    end if
    !
    ! variables required by mixing schemes
    !
    if (typeopt .eq. 2) then
        allocate (dent_after(meshproc, 2), dent_before(meshproc, 2), dentnew(meshproc))
        dent_after = 0.d0
        dent_before = 0.d0
        dentnew = 0.d0
        if (yeslsda) then
            allocate (spint_after(meshproc, 2), spint_before(meshproc, 2), spintnew(meshproc))
            spint_after = 0.d0
            spint_before = 0.d0
            spintnew = 0.d0
        end if
    end if
    if (typeopt .eq. 4) then
        allocate (dent_after(meshproc, maxold), dent_before(meshproc, maxold)&
                &, dent_aftern(meshproc, maxold), mixingtrue(maxold))
        dent_after = 0.d0
        dent_before = 0.d0
        dent_aftern = 0.d0
        mixingtrue = 0.d0
        if (yeslsda) then
            allocate (spint_after(meshproc, maxold), spint_before(meshproc, maxold), &
                      spint_aftern(meshproc, maxold))
            spint_after = 0.d0
            spint_before = 0.d0
            spint_aftern = 0.d0
        end if
        dimjac = maxold
        dimjac2 = dimjac**2
        allocate (overjac(dimjac, dimjac), jac(dimjac, dimjac), alpha(dimjac))
        allocate (mataux(dimjac, dimjac), vetaux(dimjac, 2))
        allocate (eigmolo(nelorbu))
        overjac = 0
        jac = 0
        alpha = 0
        mataux = 0
        vetaux = 0
        if (jaccond .gt. 0.d0) then
            allocate (sjac(dimjac), sojac(dimjac), ujac(dimjac, dimjac), vjac(dimjac, dimjac))
            sjac = 0
            sojac = 0
            ujac = 0
            vjac = 0
        end if
        if (yeslsda .or. ipc .eq. 2) then
            allocate (eigmolodo(nelorbu))
            eigmoldo = 0.d0
        end if
    end if
    ! save wave function values to speed up Hamiltonian update
    if (rank .eq. 0) write (6, *) ' after initialize dent ', memlarge
    if (memlarge) then
        if (ipc .eq. 1) then
            if (rank .eq. 0) write (6, '(A, E12.5, A)') " memlarge=.true. option allocates ", &
                8.d0*nelorbu*meshproc/1.d9, " Gbyte per MPI task for WF!"
            allocate (wf(nelorbu, meshproc))
            wf = 0.d0
            wf_dim = nelorbu
        else ! in this case I save the complex wave function for up/down spin electrons
            ! Allocation is halved for optimize_overs=.true.
            if (.not. double_overs) then
                if (rank .eq. 0) write (6, '(A, E12.5, A)') " memlarge=.true. option allocates ", &
                    16.d0*nelorbu*meshproc/1.d9, " Gbyte per MPI task for WF!"
                allocate (wf(2*nelorbu, meshproc))
                wf = 0.d0
                wf_dim = 2*nelorbu
            else
                if (rank .eq. 0) write (6, '(A, E12.5, A)') " memlarge=.true. option allocates ", &
                    32.d0*nelorbu*meshproc/1.d9, " Gbyte per MPI task for WF!"
                allocate (wf(4*nelorbu, meshproc))
                wf = 0.d0
                wf_dim = 4*nelorbu
            end if
        end if
    end if
    !
    ! print_energies vectors
    !
    if (print_energies) then
        allocate (dft_energies_all(nk))
        allocate (corr_energy_all(nk))
        allocate (exchange_energy_all(nk))
        allocate (number_particles(nk))
        dft_energies_all = 0.d0
        corr_energy_all = 0.d0
        exchange_energy_all = 0.d0
        number_particles = 0
        if (yeslsda .or. ipc .eq. 2) then
            if (yeslsda) then
                allocate (spin_grid_all(nk))
                spin_grid_all = 0.d0
            end if
            allocate (number_particles_do(nk))
            number_particles_do = 0
        end if
    end if
    !
    ! initialize indices for Fourier transform
    if (rank .eq. 0) write (6, *) ' Warning if you do not see after initialize fourier, try to use larger OMP_NUM_THREADS  '
    !
    call initialize_fourier()
    if (rank .eq. 0) write (6, *) ' after initialize fourier '

    ! initialize all variables
    spintot = 0.d0
    denstot = 0.d0
    mindist = 2.d-9 ! a but less accuracy as QMC calculations.
    eigmol = 0.d0
    newocc = 0.d0
    oldocc = 0.d0

    vhartree = 0.d0
    vhartreeq = 0.d0
    fp = zzero
#ifdef PARALLEL
    sndbuf = zzero
    rcvbuf = zzero
#endif

    vlocaltot = 0.d0
    vq0tots = 0.d0

    dent = 0.d0
    occupations = 0.d0
    molecorb = 0.d0
    molecorb_old = 0.d0
    premat = 0.d0
    eigmat = 0.d0
    eigmat_down = 0.d0
    if (memlarge) wf = 0.d0
    if (yeslsda .or. ipc .eq. 2) then
        eigmoldo = 0.d0
        if (yeslsda) spint = 0.d0
        oldoccdo = 0.d0
        occupationdo = 0.d0
        molecorbdo_old = 0.d0
    end if
    ! eigenvectors of the overlap matrix for
    ! up/down spin electrons
    if (ipc .eq. 2) then
        umatl = 0.d0
        umatldo = 0.d0
    end if

    occupations_sav = 0.d0
    eigmol_sav = 0.d0
    if (writescratch .eq. 0) molecorb_sav = 0.d0
    if (yeslsda .or. ipc .eq. 2) then
        occupationsdo_sav = 0.d0
        eigmoldo_sav = 0.d0
        if (writescratch .eq. 0) molecorbdo_sav = 0.d0
    end if

    hamilt = 0.d0
    if (yeslsda .or. ipc .eq. 2) hamiltdo = 0.d0
#ifdef __SCALAPACK
    hamiltl = 0.d0
    if (yeslsda .or. ipc .eq. 2) hamiltldo = 0.d0
#endif

    edft_av = 0.d0
    ecorr_av = 0.d0
    exchange_av = 0.d0
    spingrid_av = 0.d0

    ! initialization of the starting density for all
    ! mixing schemes

    if (typeopt .lt. 0) then
        molecorbs = 0.d0
        if (yeslsda .or. ipc .eq. 2) molecorbdos = 0.d0
    end if

    if (typeopt .eq. 2) then
        dent_after = dens0
        dent_before = dens0
        dentnew = dens0
        if (yeslsda) then
            spint_after = spin0
            spint_before = spin0
            spintnew = spin0
        end if

    elseif (typeopt .eq. 4) then ! default mixing scheme

        dent_after = dens0
        dent_aftern = dens0
        dent_before = dens0
        if (yeslsda) then
            spint_after = spin0
            spint_aftern = spin0
            spint_before = spin0
        end if
        jac = 0.d0
        overjac = 0.d0
        mataux = 0.d0
        vetaux = 0.d0

    end if

    if (typeopt .eq. 2) then
        dent_before = dens0
        dent_after = dens0
        if (yeslsda) then
            spint_before = spin0
            spint_after = spin0
        end if
    end if

    zgemm_time = 0.d0
    dgemm_time = 0.d0
    dgemm_timep = 0.d0
    zgemm_timep = 0.d0
    init_time = 0.d0

    call init_dielectric_dft(kappa)

    ! ---------------------------------------
    ! END ALLOCATION and INITIALIZATION
    ! ---------------------------------------

    nmoltot = nmol + ndiff
    if (.not. symmagp .or. ipc .eq. 2) nmoltot = nmoltot + nmol ! for complex wave function always use
    !
    ! reading occupations
    !
    occupations = 0.d0
    if (yeslsda .or. ipc .eq. 2) occupationdo = 0.d0
    occread = .false.
    if (nelocc .ne. 0) occread = .true.
    !
    ! read occupations from std input if required.
    !
    call read_occupations(occupations, occupationdo, nelocc, neloccdo, occread)
    !
    ! evaluate ion-ion distances
    !
    call eval_iond(iond, rion, nion, LBox, psip, iond_cart)
    kappanew = kappa
    ! updating ion-ion potential (density indipendent)
    vpotaa = upvpotaa(zetar, iond, nion, kappanew, LBox)
    if (rank .eq. 0) write (6, '(a,F20.10)') ' Asymptotic value of one-body dft pot =', vpotaa/2.d0
    !
    ! add  the q=0 contribution of the short range added  Ewald part.
    ! Correction to el-ion potential due to Ewald.
    !
    if (iespbc) then
        costkappa = vq0_diel

        vpotaa = vpotaa + weightvh*costkappa*nel**2/voltot
        if (rank .eq. 0) write (6, '(a,2F20.10)') ' Shift vpotaa =', vpotaa, weightvh*costkappa*nel**2/voltot
    end if
    !
    ! reading molecular orbitals from fort.10 put in molecorb_old
    !
    call read_molecular()
    !
    if (.not. occread) then

        countall = countalln
        count2 = count2n
        costall = costalln

        if (optocc .eq. 1) then
            countall = bands
            count2 = neldo
        end if

    end if

    if (.not. occread .and. .not. occopen .and. optocc .le. 0) then
        if (rank .eq. 0) write (6, *) ' Warning default occupation '
        occupations = 0.d0
        if (yeslsda .or. ipc .eq. 2) then
            occupationdo = 0.d0
            occupations(1:nelup) = 1.d0
            occupationdo(1:neldo) = 1.d0
        else
            occupations(1:neldo) = 2.d0
            if (nelup .gt. neldo) occupations(neldo + 1:nelup) = 1.d0
        end if
    end if

    oldocc = occupations
    if (yeslsda .or. ipc .eq. 2) oldoccdo = occupationdo
    !
    if (contracted_on) then

        if (ipc .eq. 1) then
            call dgemm('T', 'N', nelorbu, bands, nelorbh, 1.d0, mu_c, &
                       nelorbh, molecorb_old, nelorbh, 0.d0, molecorb, nelorbu)
            if (yeslsda) &
                call dgemm('T', 'N', nelorbu, bands, nelorbh, 1.d0, mu_c, &
                           nelorbh, molecorbdo_old, nelorbh, 0.d0, molecorbdo, nelorbu)
        else
            call zgemm('T', 'N', nelorbu, bands, nelorbh, zone, mu_c, &
                       nelorbh, molecorb_old, nelorbh, zzero, molecorb, nelorbu)
            call zgemm('T', 'N', nelorbu, bands, nelorbh, zone, mu_c, &
                       nelorbh, molecorbdo_old, nelorbh, zzero, molecorbdo, nelorbu)
        end if

        do i = 1, bands
            molecorb_old(1:ipc*nelorbu, i) = molecorb(1:ipc*nelorbu, i)
        end do

        if (yeslsda .or. ipc .eq. 2) then
            do i = 1, bands
                molecorbdo_old(1:ipc*nelorbu, i) = molecorbdo(1:ipc*nelorbu, i)
            end do
        end if

    else
        do i = 1, bands
            molecorb(1:ipc*nelorbu, i) = molecorb_old(1:ipc*nelorbu, i)
        end do
        if (yeslsda .or. ipc .eq. 2) then
            do i = 1, bands
                molecorbdo(1:ipc*nelorbu, i) = molecorbdo_old(1:ipc*nelorbu, i)
            end do
        end if
    end if

    if (rank .eq. 0) then
        write (6, '(a,1e18.8)') ' Molecorb read:      ', sum(abs(molecorb(:, 1:bands)))
        if (yeslsda .or. ipc .eq. 2) &
            write (6, '(a,1e18.8)') ' Molecorb read down: ', sum(abs(molecorbdo(:, 1:bands)))
    end if

    return

end subroutine initialize_environment

subroutine print_header(kaverage, decoupled_run, compute_bands)

    use allio, only: ipc, rank

    implicit none
    logical, intent(in) :: kaverage, decoupled_run, compute_bands

    if (rank .eq. 0) then
        if (kaverage) then
            write (6, *)
            write (6, *) ' --------------------------------------------------'
            write (6, *) '     DFT calculation - k-points sampling calculation '
            write (6, *) ' --------------------------------------------------'
            write (6, *)
        elseif (decoupled_run .and. kaverage) then
            write (6, *)
            write (6, *) ' -------------------------------------------------------'
            write (6, *) '     DFT calculation - indipendent k-points calculation '
            write (6, *) ' -------------------------------------------------------'
            write (6, *)
        elseif (.not. kaverage) then
            if (ipc .eq. 1) then
                write (6, *)
                write (6, *) ' --------------------------------------------'
                write (6, *) '     DFT calculation - Gamma point calculation '
                write (6, *) ' --------------------------------------------'
                write (6, *)
            else
                write (6, *)
                write (6, *) ' -----------------------------------------------'
                write (6, *) '     DFT calculation - single phase calculation '
                write (6, *) ' -----------------------------------------------'
                write (6, *)
            end if
        elseif (compute_bands) then
            write (6, *) ' '
            write (6, *) ' --------------------------------------------------------'
            write (6, *) '     DFT calculation - non self-consistent run           '
            write (6, *) ' --------------------------------------------------------'
            write (6, *) ' '
        end if
    end if

    return

end subroutine print_header

!-----------------------------------------------------------------
subroutine initialize_mesh(nx, ny, nz, ax, ay, az, fx, fy, fz, &
                           mesh, i8cost, volmesh, meshproc, meshproc_tot)
    !-----------------------------------------------------------------

    ! This subroutine initialize all variables related to the
    ! molecular orbitals integration mesh. The mesh is chosen to
    ! be as far as possible from ion positions in order to avoid
    ! divergences in the potential.
    ! For k-points sampling we use the same mesh for all k-points.

    use constants, only: pi
    use allio, only: rion, zmax, zetar, iespbc, nion, nprocrep, rank
    use setup, only: lxp, lyp, lzp, nbufd, rion_shift, shift_origin, &
                     rion_ref, totnbuf, bufbuf, mind, shiftx, shifty, shiftz, &
                     double_mesh, scale_z, l0_at, nx0, ny0, nz0, minz_at, nx_at, ny_at, nz_at, &
                     rion_from, from_ions
    use cell, only: yes_tilted, unit_volume, at, CartesianToCrystal

    implicit none

    integer, intent(inout) :: nx, ny, nz, meshproc, meshproc_tot
    integer(8), intent(out) :: mesh, i8cost
    real(8), intent(inout) :: ax, ay, az, volmesh, fx, fy, fz
    integer i, j, k, ind, ii, scalea
    real(8) axt, ayt, azt
    integer(8) mesh_try
    real*8, allocatable :: rion_sav(:, :)
#ifdef __DEBUG
    real(8) x(3)
#endif

    mesh = nx
    mesh = mesh*ny
    mesh = mesh*nz
    volmesh = ax*ay*az*unit_volume

    lxp = nx*ax
    lyp = ny*ay
    lzp = nz*az

    if (yes_tilted) then

        fx = 1.d0
        fy = 1.d0
        fz = 1.d0

    else

        fx = 2.0d0*pi/lxp
        fy = 2.0d0*pi/lyp
        fz = 2.0d0*pi/lzp
        !

    end if
    ! Compute the first point of the mesh i,j,k=1 in rion_ref
    !
    if (.not. iespbc) then
        ! the reference is the average ion position
        do j = 1, 3
            rion_ref(j) = sum(rion(j, :))/nion
        end do
        rion_ref(1) = rion_ref(1) - (nx - 1)/2.d0*ax
        rion_ref(2) = rion_ref(2) - (ny - 1)/2.d0*ay
        rion_ref(3) = rion_ref(3) - (nz - 1)/2.d0*az
    else
        rion_ref = 0.d0
    end if
    if (shiftx) rion_ref(:) = rion_ref(:) + ax/2.d0*at(:, 1)
    if (shifty) rion_ref(:) = rion_ref(:) + ay/2.d0*at(:, 2)
    if (shiftz) rion_ref(:) = rion_ref(:) + az/2.d0*at(:, 3)
    !
    ! Optimze the mesh to be as far as possible from the ions
    if (shift_origin) then

        if (double_mesh) then
            axt = ax/scale_z ! To avoid occasional collisions
            ayt = ay/scale_z ! To avoid occasional collisions
            azt = az/scale_z ! To avoid occasional collisions
        else
            axt = ax
            ayt = ay
            azt = az
        end if
        allocate (rion_sav(3, nion))
        rion_sav = rion

        if (iespbc) then
            call CartesianToCrystal(rion_sav, nion)
            call CartesianToCrystal(rion_ref, 1)
        end if

        call findrionfref(nion, axt, rion_sav, rion_ref, mind(1))
        call findrionfref(nion, ayt, rion_sav(2, 1), rion_ref(2), mind(2))
        call findrionfref(nion, azt, rion_sav(3, 1), rion_ref(3), mind(3))

        if (iespbc) rion_ref(:) = rion_ref(1)*at(:, 1) + rion_ref(2)*at(:, 2) + rion_ref(3)*at(:, 3)

        deallocate (rion_sav)

    else
        ! Since the error is quadratic it should give an error ~1e-8 on the DFT
        ! point symmetries. Remind the point symmetries refers to the origin in Turbo
        ! and a too large shift breaks the symmetries.
        axt = (1.d0 + sqrt(5.d0))/10.d0**6 ! To avoid occasional collisions
        ayt = (1.d0 + sqrt(5.d0))/10.d0**6 ! To avoid occasional collisions
        azt = (1.d0 + sqrt(5.d0))/10.d0**6 ! To avoid occasional collisions
        allocate (rion_sav(3, nion))
        rion_sav = rion

        if (iespbc) then
            call CartesianToCrystal(rion_sav, nion)
            call CartesianToCrystal(rion_ref, 1)
        end if

        call findrionfref(nion, axt, rion_sav, rion_ref, mind(1))
        call findrionfref(nion, ayt, rion_sav(2, 1), rion_ref(2), mind(2))
        call findrionfref(nion, azt, rion_sav(3, 1), rion_ref(3), mind(3))

        if (iespbc) rion_ref(:) = rion_ref(1)*at(:, 1) + rion_ref(2)*at(:, 2) + rion_ref(3)*at(:, 3)

        deallocate (rion_sav)

        !  The i=1 corresponds to the origin

    end if

    if (rank .eq. 0) write (6, *) ' Minimum ion-mesh distance =', dsqrt(sum(mind(:)**2))
    if (rank .eq. 0) write (6, *) ' Origin shift used =', rion_ref(:)
    rion_shift = rion_ref

    if (.not. from_ions .and. double_mesh .and. rion_from(1) .eq. 1.d23) then
        if (rank .eq. 0) write (6, *) ' Warning default center position as reference'
        do j = 1, 3
            rion_from(j) = sum(rion(j, :))/nion
        end do
    end if

    rion_ref(:) = rion_ref(:) + (nx - 1)/2.d0*ax*at(:, 1)
    rion_ref(:) = rion_ref(:) + (ny - 1)/2.d0*ay*at(:, 2)
    rion_ref(:) = rion_ref(:) + (nz - 1)/2.d0*az*at(:, 3)
    mesh_try = mesh
    i8cost = mesh_try/nprocrep
    meshproc = i8cost
    i8cost = i8cost*nprocrep
    i8cost = mesh_try - i8cost
    if (i8cost .ne. 0) meshproc = meshproc + 1

    if (rank .eq. 0) write (6, *) 'New center of mesh =', rion_ref(:)
    if (double_mesh) then

        if (2*l0_at .lt. ax) then
            l0_at = (ax + 1d-6)/2.d0
            if (rank .eq. 0) write (6, *) ' Warning l0_at changed to ', l0_at
        end if
        if (2*l0_at .lt. ay) then
            l0_at = (ay + 1d-6)/2.d0
            if (rank .eq. 0) write (6, *) ' Warning l0_at changed to ', l0_at
        end if
        if (2*l0_at .lt. az) then
            l0_at = (ay + 1d-6)/2.d0
            if (rank .eq. 0) write (6, *) ' Warning l0_at changed to ', l0_at
        end if

        if (nx_at .gt. 0) then
            nx0 = nx_at
            ny0 = ny_at
            nz0 = nz_at
        else
            nx0 = (2.d0*l0_at)/ax + 1
            ny0 = (2.d0*l0_at)/ay + 1
            nz0 = (2.d0*l0_at)/az + 1
        end if

        if (rank .eq. 0) then
            write (6, *) ' Small box for finer mesh (a.u.) =', (nx0 - 1)*ax, (ny0 - 1)*ay, (nz0 - 1)*az
            write (6, *) ' Small grid for finer mesh =', nx0, ny0, nz0
        end if

        scalea = scale_z
        if (.not. from_ions) then
            i8cost = scalea*(nx0 - 1) + 3
            i8cost = i8cost*(scalea*(ny0 - 1) + 3)
            i8cost = i8cost*(scalea*(nz0 - 1) + 3)
            mesh_try = i8cost/nprocrep
            i8cost = i8cost - nprocrep*mesh_try
            if (i8cost .ne. 0) mesh_try = mesh_try + 1
            meshproc = meshproc + mesh_try
            i8cost = nx0 + 2
            i8cost = i8cost*(ny0 + 2)
            i8cost = i8cost*(nz0 + 2)
            mesh_try = i8cost/nprocrep
            i8cost = i8cost - nprocrep*mesh_try
            if (i8cost .ne. 0) mesh_try = mesh_try + 1
            meshproc = meshproc + mesh_try
        else
            do ii = 1, nion
                if (zetar(ii) .ge. minz_at) then
                    i8cost = nx0 + 2
                    i8cost = i8cost*(ny0 + 2)
                    i8cost = i8cost*(nz0 + 2)
                    mesh_try = i8cost/nprocrep
                    i8cost = i8cost - nprocrep*mesh_try
                    if (i8cost .ne. 0) mesh_try = mesh_try + 1
                    meshproc = meshproc + mesh_try
                    i8cost = scalea*(nx0 - 1) + 3
                    i8cost = i8cost*(scalea*(ny0 - 1) + 3)
                    i8cost = i8cost*(scalea*(nz0 - 1) + 3)
                    mesh_try = i8cost/nprocrep
                    i8cost = i8cost - nprocrep*mesh_try
                    if (i8cost .ne. 0) mesh_try = mesh_try + 1
                    meshproc = meshproc + mesh_try
                end if
            end do
        end if
    end if
    meshproc_tot = meshproc

    if (meshproc .ge. nbufd) then
        bufbuf = nbufd
    else
        bufbuf = meshproc
        if (rank .eq. 0) write (6, *) ' Warning, your buffer is too large &
                &, working with the maximum one =', meshproc
    end if
    ! total number of buffers required
    totnbuf = ceiling(dble(meshproc)/dble(bufbuf))

    if (rank .eq. 0) write (6, *) ' Total number of buffers (lower bound) ', totnbuf

#ifdef __DEBUG
    if (rank .eq. 0) then
        write (6, *) ' Writing grid positions on a file '
        open (unit=91, file='turbogrid.dat', form='formatted', position='rewind', status='unknown')
        ind = 0
        do k = 1, nz
!     x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
!        x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
!           x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)
                    x(:) = (-(nz + 1)/2.d0 + k)*az + rion_ref(:) + (-(ny + 1)/2.d0 + j)*ay*at(:, 2) + (-(nz + 1)/2.d0 + k)*at(:, 3)
                    ind = ind + 1
                    write (91, '(I9,2X,3F12.6)') ind, x(:)
                end do
            end do
        end do
        close (91)

    end if
#endif

end subroutine initialize_mesh
