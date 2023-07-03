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

module setup

    use allio
    use sub_comm
    use IO_m
#ifdef __SCALAPACK
    use descriptors
#endif

    implicit none

    public
    !
    ! INPUT PARAMETERS
    !
    integer :: maxit ! maximum number of SC iterations allowed
    integer :: typedft ! type of DFT functional
    integer :: typeopt ! type of mixing scheme
    integer :: optocc ! electron occupation type (0 = fixed occupations, 1 = smearing)
    integer :: nelocc ! number of occupied electronic states
    integer :: neloccdo ! number of occupied down-spin electronic states (for option yeslsda=.true.)
    integer :: bands ! number of electronic levels/bands to be computed for each k-point
    integer :: bufbuf ! dimension of the buffer for real space integration
    integer :: bandso, nproco, meshproco, mesho, nproc_fft ! number of bands/processors/mesh points read from occupationlevels.dat
    integer :: nxs, nys, nzs ! sizes of input spin grid
    real(8) :: epsdft ! convergence thresold for SC cycle
    real(8) :: mixing ! value of mixing coefficient
    real(8) :: epsover ! precision on overlap preconditioning for diagonalization
    real(qp) :: epsshell ! value of electronic temperature (only Fermi-Dirac smearing implemented so far)
    real(8) :: eps_mach ! machine precision
    real(8) :: h_field, h_charge ! value of external magnetic field
    real(8) :: randspin ! if > 0 adds a random component to spin up/spin down starting eigenvectors
    real(8) :: deltaE ! energy bin for DOS calculation
    real(8) :: emin, emax ! maximum/minimum energy for bands/DOS calculations
    real(8) :: mixingder
    logical :: contracted_on ! if .true. enables calculation on contracted basis set
    logical :: memlarge ! if .true. keep teh wave function to speed up computations
    logical :: yeslsda ! flag for collinear spin calculations
    logical :: orthodiag ! orthogonalized eigenvectors after a SC cycle
    logical :: occopen ! flag for reading occupations from input (see optocc options)
    logical :: write_den ! if .true. print final density on a grid
    logical :: yesgrad ! evaluate the gradient of the electronic density (for PBE functional not implemented yet!)
    logical :: yespassed ! flag to determine if a SC is correctly converged
    logical :: noread ! if .true. all overlap matrix elements are computed from scratch (determined during the calculation)
    logical :: yesreadf ! if .true. charge/spin density are read from scratch file (determined during the calculation)
    logical :: changemix ! flag to enable adaptive mixing of the density
    logical :: optimize_overs ! if .true. choose to optimize overlaps computation (both time and memory)
    ! if possible, otherwise compute all from scratch
    logical :: double_overs ! if .true. the code computes overlaps for up/down spin electrons from scratch
    logical :: print_energies ! print out energies for each k-points: first the Efree(k) = non interacting energy
    ! then the Etot(k) = total DFT energy at the end of the self-consistent cycle
    !  logical  :: shift_origin,shiftx,shifty,shiftz     ! apply a shift to the origin of the real space grid, defined in allio.
    logical :: corr_hartree
    logical :: do_hartree
    logical :: try_translation ! sort molecular orbitals according to their overlap (in the case of k_down = - k_up)
    logical :: zero_jas ! put the one-body Jastrow to zero at the end of the calculation
    logical :: write_matrix ! write overlap matrices on a file for effective Hamiltonian calculations
    logical :: fix_density ! evolve k-points independently but with averaged density, to be used with decoupled_run=.true.
    logical :: linear ! Linear fit of the hartree potential
    logical :: from_ions !  if true reference is from ions
    logical :: commensurate_lattice
    integer :: maxcg, maxold, itersto, nmol_in, molecsw, mincond, i_new, j_new, lworkjac, nc, nx_at, ny_at, nz_at
    integer*2, allocatable :: spin_input(:, :, :) ! input spin grid for LSDA calculations
    integer*2, allocatable :: charge_input(:, :, :) ! input charge input grid for LSDA/LDA  calculations
    real(8) :: epssr, weightcorr, weightxc, costalln, spingrid, jaccond, &
               weightvh, vlocaltot, vq0tots, kappanew, eh_ew, newj_onebody(100) &
               , setj_onebody(100), l0_at, minz_at, vmax0, vmax0_in, scale_hartree, scale_hartreen

    !
    ! REAL SPACE INTEGRATION MESH
    !
    integer(8) :: mesh, i8cost ! mesh = total # of mesh points
    ! meshproc = # of mesh points per processor
    integer :: indmesh, nbuf, totnbuf, meshproc, meshproc_lsda, meshproc_wf, meshproc_tot
    real(8) :: volmesh, rion_ref(3), rion_from(3), rion_shift(3), mind(3), &
               lxp, lyp, lzp, rion_upload(3) ! reciprocal mesh size
    !
    ! K-POINTS
    !
    integer :: indk ! k-point index inside processors pool
    integer :: nko ! # of k-points read in occupationalevels.dat file
    integer :: nk_rest ! # of k-points which are disregarded due to wrong number of processors
    integer :: task ! defines the task to be performed after NonSC calculation.
    ! task=1 -> band structure calculation
    ! task=2 -> DOS calculation
    integer :: nx0, ny0, nz0 ! mesh grid for the atomic integration up to a cutoff l0_at
    !
    ! FILES UNITS
    !
    integer, parameter :: ufort10 = 10 ! fort.10/fort.10_new in input
    integer, parameter :: unit_scratch_distributed = 19 ! real space grid distributed quantities
    integer, parameter :: unit_scratch_fort10 = 20 ! fort.10 for each k-point
    integer, parameter :: unit_scratch_densities = 21 ! total charge/spin densities
    integer, parameter :: unit_scratch_eigenvects = 22 ! LDA eigenvectors
    integer, parameter :: unit_print_energies = 23 ! energies for each k-point (useful for searching special k_point)
    !
    ! DIMENSIONS
    !
    integer :: nelorbu, nelorbudo ! leading dimension of main matrices  (normal version)
    integer :: nlax, nlaxdo ! leading dimension of block matrices (__SCALAPACK)
    integer :: nelorb3 ! leading dimension of buffers
    integer :: wf_dim ! leading dimension of the wf array, memlarge option
    integer :: bandsdo
    integer :: nelorbdo
    integer :: scale_z
    ! 1body Jastrow parameter
    integer :: iesdvj
    type(mpi_sub_comm) :: sub_comm_fft
#ifdef __SCALAPACK
    integer, public :: desch(16)
#endif

    !
    ! MAIN MATRICES
    !
    real(8), allocatable, dimension(:, :) :: &
        overs, overinvs, overinvsl, oversdo, & ! matrix of non-orthogonal basis overlaps: < \Psi_i | \Psi_j >
        overham, overhamdo, & ! hamiltonian matrix independent from the density
        hamilt, hamiltdo, & ! hamiltonian matrix: < \Psi_i | H | \Psi_j >
        wf, & ! wave functions values saved for continuation with memlarge option
        molecorb, molecorbdo, molecorb_old, & ! eigenvectors of KS hamiltonian
        molecorbs, molecorb_part, molecorbdo_old, molecorbdos
    real(8), dimension(:, :), allocatable, public :: umatl, umatldo ! overlap matrix eigenvectors
#ifdef __SCALAPACK
    real(8), allocatable, dimension(:, :), public :: molecorbl, molecorbldo, &
                                                     hamiltl, overhaml, hamiltldo, oversldo, overhamldo
    real(8), allocatable, dimension(:, :), public :: oversl
#endif
    real(8), external :: ran
    !
    ! hamiltonian and density
    !
    logical :: occread, optpar, yescross, check, done
    logical, allocatable, dimension(:) :: occorb, gridspin, gridnospin&
            &, gridcharge, gridnocharge

    real(8), dimension(:), allocatable :: occupations, occupationdo, oldocc, oldoccdo, newocc, & ! occupations of KS eigenvalues
                                          eigmol, eigmoldo, eigmolo, eigmolodo, & ! eigenvalues (spin up/down) of KS Hamiltonian
                                          dent, volmesh_proc, dentnew, dentold, spint, spintold, spintnew, & ! loc charge/spin dens
                                          premat, eigmat, eigmat_down, zall, fh, fs
    real(8), dimension(:, :), allocatable :: occupations_sav, occupationsdo_sav, &
                                             eigmol_sav, eigmoldo_sav, buffer_grid ! to save important quantities
    real(8), dimension(:, :, :), allocatable :: dens_grid, spint_grid, molecorb_sav, molecorbdo_sav
    real(8) :: energy, vpotaa, vpotden, errdft, edft, eh, &
               mindist, dens, totdens, gradmod, &
               totvpot, voltot, voltot_double, dens0, spin0, rsdft, ecorr, exchange, lapden, gradden(3), totnorm, &
               dens_true, ehartree, betaden, edftp, edftb, costder, mixingstep, &
               der_after, der_before, der2, trystep, der1cg, der2cg, edft_ref, lambdacg, gammacg, &
               lambdastep, vh_shift, sumold, errsav, edftvar, mincost, mincosts, &
               denstot, vh_test, vh_att, spintot, spin2tot, betaspin, betas, &
               v_field, zgmax, fsmesh, fhmesh, ac, rmax2, zfirst, zfirsts, lprecdft, maxint, &
               loading_time, diag_time, dgemm_time, dgemm_timep, zgemm_time, zgemm_timep &
               , symtime, cycle_time, init_time, time_fft, time_uploadfft, time_total

    ! mixing schemes
    real(8), dimension(:, :), allocatable :: dent_before, dent_after, dent_aftern & ! needed by mixing typeopt=2
                                             , spint_before, spint_after, spint_aftern, & ! needed by mixing typeopt=2
                                             overjac, jac, mataux, vetaux, ujac, vjac, & ! needed by Anderson mixing typeopt=4
                                             gradt
    real(8), dimension(:), allocatable :: mixingtrue, alpha, sjac, sojac, weightx, weighty, weightz
    !
    ! definition of the namelist containing all DFT input parameters. This namelist
    ! is only used for DFT calculations and it does not appear in the main QMC code.
    !
    namelist /DFT/ mixing, typedft, maxit, epsdft, optocc, typeopt, memlarge, &
        epsover, nelocc, neloccdo, bands, mixingder, epssr, maxcg, weightcorr, &
        weightxc, maxold, epsshell, weightvh, orthodiag, mincond, &
        randspin, nzs, nys, nxs, jaccond, h_field, eps_mach, nc, rmax, &
        write_den, contracted_on, task, deltaE, emin, emax, optimize_overs, &
        shift_origin, try_translation, zero_jas, write_matrix, fix_density, &
        newj_onebody, setj_onebody, h_charge, shiftx, shifty, shiftz, &
        l0_at, scale_z, linear, corr_hartree, rion_from, minz_at, vmax0_in, &
        nx_at, ny_at, nz_at, from_ions, scale_hartree, do_hartree, nproc_fft

    namelist /band_structure/ task, deltaE, emin, emax

contains

    !-------------------------
    subroutine Initializeall
        !-------------------------

        implicit none

        integer :: i, j, k, kk, dimspin_input, nel_neutral, thread_active, rankrep_sav&
                &, commrep_mpi_sav, maxdiv, minn, jj
        character(len=80) :: charascratch
        character(lchlen) :: path, scratchpath
        real*8 rion_ref0(3), rc(3), r0, cost_double
        real*8, external :: jastrow_ei
#ifdef PARALLEL
        include "mpif.h"
#endif
#if defined(_OPENMP)
        integer, external :: omp_get_max_threads
        thread_active = omp_get_max_threads()
#else
        thread_active = 1
#endif
        !
        iflagerr = 0
        iflagerrall = 0
        iflagnorm = 3 ! normalization always computed
        !
        ! Open scratch files.
        ! They include one file "tmp00xxxx" on each processor. The first record
        ! of these files contains the partial density (or the partial density and spin density)
        ! associated to that specific processor. The second record
        ! contains the basis (variable overs) and hamiltonian (variable
        ! overham) overlap matrix elements.
        !
        if (writescratch .ne. 0) &
            call error(' prep ', ' Scratch files not written ', -1, rank)

        call check_scratch(rank, path, scratchpath)

        ! scratch files distributed on the real space grid (tmp000xxx)
        !
        if (writescratch .eq. 0) then
            open (unit=unit_scratch_distributed, file=trim(wherescratch)//'.'//chara, form='unformatted', status='unknown')
        end if
        !
        ! reading/writing wave functions for all k-points (fort10_000xxx)
        !

        if ((yeswrite10 .or. yesread10) .and. rankrep .eq. 0) then

#ifdef __KCOMP
            charascratch = trim(scratchpath)//'fort.10_'//trim(chara)
#else
            charascratch = trim(scratchpath)//'fort.10_'//trim(charaq)
#endif
            open (unit=unit_scratch_fort10, file=charascratch, form='formatted', position='rewind', status='unknown')
        end if
        !
        ! open the formatted files to save LDA quantities needed by post processing tools
        ! and nonSC run:
        ! total_densities.sav -> charge/spin densities
        ! wavefunction.sav    -> eigenvectors
        !
        if (rank .eq. 0 .and. writescratch .eq. 0) then
            ! save total charge/spin densities
            charascratch = trim(scratchpath)//'total_densities.sav'
            open (unit=unit_scratch_densities, file=charascratch, form='unformatted', position='rewind', status='unknown')
            charascratch = trim(scratchpath)//'wavefunction.sav'
            open (unit=unit_scratch_eigenvects, file=charascratch, form='unformatted', position='rewind', status='unknown')
        end if
        !
        ! reading pseudo potential if any
        !
        mindist = 1d-8 ! Initialization minimo e-i distance.

        call read_pseudo()
        pseudologic = pseudorandom
        !
        ! these indices are required by the "write_fort10" at the
        ! end of a SC calculation.
        !
        iseed = abs(iseedr)
        nmat = np + 1
        ndim = np + 1
        npm = np + 1
        npmn = nbinmax + 1
        nmats = nmat - iesm - iesd
        ndimp = np
        ndimpdim = max(ndimp, 1)
        nwm = nw + 1
        nws = in1 + 1
        nwm = nw
        nws = in1
        npp = np + 1
        npf = np
        nprest = np
        nindt = 0
        ! initialize random number generator anyway
        irstart = iseed
        do i = 1, rank
            irstart = ran(iseed)*2**29
        end do
        iseed = 2*irstart + 1
        call rand_init(iseed)
        !
        ! reading the wave function fort.10/fort.10_new
        ! Get all the parameters related to the system.
        !
        ! with the new convention for read_fort10, we have to open
        ! the file for each processor pool.
        if (rank .eq. 0) then
            if (iopt .eq. 1) then
                open (unit=ufort10, file='fort.10', form='formatted', status='unknown')
            else
                open (unit=ufort10, file='fort.10_new', form='formatted', status='unknown')
            end if
        end if

        if ((iopt .ne. 1 .and. kaverage .and. yesread10) .or. (manyfort10 .and. yesread10 .and. nbead .gt. 1)) then
            !
            ! In case of yeswrite10=.true. and continuation
            ! read molecular orbitals from each turborvb.scratch/fort.10_00xxxx and
            ! fill the matrix dup_c for all k-points separately.
            !
            indk = rankcolrep + 1
            call read_fort10(unit_scratch_fort10)
            ! if k-point read from input is not correct, substitute it with the w.f. phase
            if (all(abs(phase(:) - xkp(:, indk)) .gt. 1d-10) .or. &
                all(abs(phase_down(:) - xkp_down(:, indk)) .gt. 1d-10)) then
                xkp(:, indk) = phase(:)
                xkp_down(:, indk) = phase_down(:)
                call error(' prep ', ' Some input k-points are wrong, substituted with WF phase! ', -1, rank)
            end if
        else
            ! in the other cases (start from scratch, or wave functions for each k-point not
            ! saved to disk, then just read the same fort.10 for all.
            commrep_mpi_sav = commrep_mpi
            rankrep_sav = rankrep
#ifdef PARALLEL
            commrep_mpi = MPI_COMM_WORLD
#else
            commrep_mpi = 0
#endif
            call read_fort10(ufort10)
            commrep_mpi = commrep_mpi_sav
            rankrep = rankrep_sav
        end if

        !
        ! printout parameters related to Jastrow 1body
        if (rank .eq. 0) then
            write (6, *) ' symmagp =', symmagp
            write (6, *) 'Number 1body Jastrow parameters: ', niesd
            if (niesd .ge. 1) then
                do i = 1, niesd
                    write (6, *) i, vj(i)
                end do
            else
                write (6, *) '1body Jastrow not present'
            end if
            if (iesdr .le. -5) then
                write (6, *) ' initial costz, costz3, zeta_Q_Caffarel '
                do i = 1, nion
                    write (6, *) i, costz(i), costz3(i), zetaq(i)
                end do
            end if
        end if
        do j = 1, nw
            jbra(j) = j
        end do
        ! main check after reading the wave function
        if (contraction .eq. 0) &
            call error(' setup_calculation ', ' DFT works only with molecular orbitals ', 1, rank)
        !
        ! read card MOLECULAR from std input
        ! and set default values
        !
        ax = 0.d0
        ay = 0.d0
        az = 0.d0
        nx = 0
        ny = 0
        nz = 0
        nmolmax = 0
        nmolmin = 0
        nmol = 0
        yesmin = 0
        weight_loc = -1.d0
        power = 1.d0
        call read_datasmin_mol()
        !
        ! read DFT cards from std input
        ! and set default values
        !
        ! &DFT card
        maxit = 100
        mixing = -1.d0
        typedft = 1
        epsdft = 1.d-7
        optocc = -1 ! Do not optimize the occupation
        orthodiag = .false. ! no orthogonalization after diagonalization
        typeopt = 4 ! default Anderson mixing self-consistency
        memlarge = .true. ! default optimize CPU.
        epsover = 1d-13 ! Almost machine precision condition number works.
        mincond = 1 ! First orbital to consider regardless of cond number.
        epssr = -1.d0 ! no further regularization of overlap inverse
        eps_mach = -1.d0 ! default machine precision
        nelocc = 0
        neloccdo = 0 ! Used for LSDA
        bands = 2*neldo + nelup - neldo ! Assuming at most 8 fold degeneracy last shell
        mixingder = -1.d0
        maxcg = 40
        weightcorr = 1.d0 ! the amount of correlation, weightcorr=0--> fake HF.
        epsshell = 0.005d0 ! energy level difference to consider it degenerate
        weightxc = 1.d0 ! the amount of exchange
        weightvh = 1.d0 ! the weight for the Hartree potential
        if (iopt .eq. 0) then
            randspin = 0.d0 ! do nothing in continuation
        else
            randspin = -1.d0 ! if>0 add random component to magnetization on the defined grid.
            ! if < 0 initialize spin according to grid if defined
        end if
        h_field = 0.d0 ! external magnetic field
        h_charge = 0.d0 ! external charge  field
        maxold = 3 ! maximum number of previous iterations for which density is saved
        nzs = -1 ! nxs,nys,nzs are the sizes of the spin grid
        nys = -1
        nxs = -1
        jaccond = 0.d0
        nc = -1
        rmax = 0.d0
        write_den = .false.
        contracted_on = .false.
        optimize_overs = .true.
        print_energies = .true.
        fix_density = .false.
        newj_onebody = -1.d0
        setj_onebody = -1.d0
        ! &band_structure card
        task = 0
        emin = 0.d0
        emax = 0.d0
        deltaE = 0.d0
        linear = .false.
        from_ions = .true.
        nx_at = -1
        ny_at = -1
        nz_at = -1
        l0_at = 1.d0
        minz_at = -1.d0
        if (double_mesh) then
            scale_z = 2
        else
            scale_z = 1
        end if

        if (.not. skip_equivalence .and. yes_tilted .and. rank .eq. 0) then
            write (6, *) ' Warning Point group symmetries not implementes with yes_tilted'
        end if

        if (skip_equivalence .or. yes_tilted) then
            shift_origin = .true.
        else
            shift_origin = .false.
        end if
        if (rank .eq. 0) write (6, *) ' Default value of shift_origin', shift_origin
        shiftx = .false.
        shifty = .false.
        shiftz = .false.

        try_translation = .true.
        do_hartree = .false.
        if (double_mesh) then
            scale_hartree = 1.d0
            corr_hartree = .true.
        else
            scale_hartree = -1.d0
            corr_hartree = .false.
        end if

        zero_jas = .false.
        write_matrix = .false.
        rion_from(1) = 1.d23
        vmax0_in = -100.d0
        nproc_fft = -1
        !
        if (rank .eq. 0) then

            iflagerr = 1
            read (5, nml=DFT, err=200)
            iflagerr = 0
200         if (iflagerr .ne. 0) call error(' setup_calculation ', ' error in reading DFT card ', 1, rank)

            if (compute_bands) then
                iflagerr = 1
                read (5, nml=band_structure, err=201)
                iflagerr = 0
201             if (iflagerr .ne. 0) call error(' setup_calculation ', ' error in reading band_structure card ', 1, rank)
            end if
            if (nx_at .gt. 0 .and. ny_at .lt. 0) then
                ny_at = nx_at
                write (6, *) ' Default value of ny_at =', ny_at
            end if
            if (ny_at .gt. 0 .and. nz_at .lt. 0) then
                nz_at = ny_at
                write (6, *) ' Default value of nz_at =', nz_at
            end if
            if (double_mesh) then
                write (6, *) ' Input scale for finer mesh  =', scale_z
                if (from_ions) then
                    if (l0_at .gt. 0 .and. nx_at .lt. 0) then
                        write (6, *) ' Finer mesh on a cube of side ', 2*l0_at, &
                            'around EACH atom'
                    else
                        write (6, *) ' Finer mesh on a ortho cell of side ', nx_at, ny_at, nz_at &
                            , ' unit ax,ay,az original mesh around EACH atom '
                    end if
                else
                    if (l0_at .gt. 0 .and. nx_at .lt. 0) then
                        write (6, *) ' Finer mesh on a cube of side ', 2*l0_at&
                                &, ' referenced from the center'
                    else
                        write (6, *) ' Finer mesh on a ortho cell of side ', nx_at, ny_at, nz_at &
                            , ' unit ax,ay,az original mesh, referenced from the center  '
                    end if
                end if
            end if
            !    The integration mesh includes the points at the boarder
            if (nx_at .gt. 0) nx_at = nx_at + 1
            if (ny_at .gt. 0) ny_at = ny_at + 1
            if (nz_at .gt. 0) nz_at = nz_at + 1
            if (bands .le. nelup) then
                write (6, *) ' Warning # bands > neldo, restoring the default '
                bands = 2*neldo + nelup - neldo
            end if

        end if
        !
        ! check input and set default values
        !
        call check_dft_input()
        call checkiflagerr(iflagerr, rank, "ERROR in TurboDFT setup")

#ifdef PARALLEL
        ! DFT self-consistent calculation
        call mpi_bcast(epsdft, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(epsover, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(eps_mach, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(mincond, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(epssr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(weightcorr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(weightxc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(weightvh, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(randspin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(h_field, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(h_charge, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(jaccond, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(mixing, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        cost_double = epsshell
        call mpi_bcast(cost_double, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        epsshell = cost_double
        call mpi_bcast(mixingder, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(rmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(maxit, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nelocc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(neloccdo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(typedft, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(typeopt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(bands, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(maxcg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(maxold, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(optocc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(orthodiag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(memlarge, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nxs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nys, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nzs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(write_den, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(contracted_on, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(yeswrite10, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(yesread10, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(iopt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(decoupled_run, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(noread, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(yesreadf, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(changemix, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(yeslsda, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(kaverage, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(optimize_overs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(print_energies, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(shift_origin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(shiftx, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(shifty, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(shiftz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(try_translation, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(corr_hartree, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(do_hartree, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(scale_hartree, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(zero_jas, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(write_matrix, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(fix_density, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(linear, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(from_ions, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        ! band structure calculation
        call mpi_bcast(deltaE, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(emin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(emax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(newj_onebody, 100, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(setj_onebody, 100, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(l0_at, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(rion_from, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(minz_at, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(vmax0_in, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nx_at, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(ny_at, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nz_at, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(scale_z, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(task, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(nproc_fft, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
        if (setj_onebody(1) .gt. 0 .and. abs(iesdr) .gt. 4) then
            if (abs(niesd) .le. 2) then
                vj(abs(niesd)) = setj_onebody(1)
            else
                do j = 2, niesd
                    vj(j) = setj_onebody(j - 1)
                end do
            end if
            !     recompute scale_one_body
            !  The value of the exp for r=0
            if (n_body_on .ne. 0) then
                scale_one_body = 0.d0
                if (.not. iespbc) then
                    do j = 1, 3
                        rion_ref0(j) = sum(rion(j, :))/nion
                    end do
                else
                    rion_ref0(:) = 0.d0
                end if
                do jj = 1, nion
                    if (iespbc) then
                        rc(1:3) = rion_ref0(1:3) - rion(1:3, jj)
                        call CartesianToCrystal(rc, 1)
                        do kk = 1, 3
                            rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                        end do
                        r0 = norm_metric(rc, metric)
                    else
                        rc(:) = (rion_ref0(:) - rion(:, jj))*costz(jj)
                        r0 = dsqrt(sum(rc(:)**2))
                    end if
                    scale_one_body = scale_one_body - &
                            &jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                end do

            else
                scale_one_body = 0.d0
            end if
            if (rank .eq. 0) write (6, *) ' Reccomputed scale one body =', scale_one_body
        end if

        if (do_hartree .and. corr_hartree .or. .not. double_mesh) then
            if (rank .eq. 0 .and. do_hartree) write (6, *) ' Warning do_hartree incompatible with corr_hartree, set last to false'
            if (rank .eq. 0 .and. .not. double_mesh) write (6, *) ' Warning corr_hartree set to false'
            corr_hartree = .false.
        end if
        if (do_hartree .and. .not. double_mesh) then
            do_hartree = .false.
            if (rank .eq. 0) then
                write (6, *) ' Warning do_hartre=.true. possible only with double_mesh=.true., changed do_hartree=.false. '
            end if
            do_hartree = .false.
        end if

        if (nproc_fft .lt. 0) then
            if (do_hartree .and. scale_z .ne. 1) then
                minn = nx*scale_z
                if (nz*scale_z .lt. minn) minn = nz*scale_z
                if (nproc_fft .gt. 0 .and. nproc_fft .lt. nprocrep) then
                    if (minn .gt. nproc_fft) minn = nproc_fft
                else
                    if (minn .gt. nprocrep) minn = nprocrep
                end if
                maxdiv = 0
                do i = 1, minn
                    if (mod(nx*scale_z, i) .eq. 0 .and. mod(nz*scale_z, i) .eq. 0) maxdiv = i
                end do
            else
                minn = nx
                if (nz .lt. minn) minn = nz
                if (nproc_fft .gt. 0 .and. nproc_fft .lt. nprocrep) then
                    if (minn .gt. nproc_fft) minn = nproc_fft
                else
                    if (minn .gt. nprocrep) minn = nprocrep
                end if
                maxdiv = 0
                do i = 1, minn
                    if (mod(nx, i) .eq. 0 .and. mod(nz, i) .eq. 0) maxdiv = i
                end do
            end if

            !         if(iespbc) then
            nproc_fft = maxdiv
            !         else
            !           if(2*maxdiv.le.nprocrep) then
            ! Try to  double the number of processors
            !           nproc_fft=2*maxdiv
            !           else
            !           nproc_fft=maxdiv
            !           endif
            !         endif

        end if

        commensurate_lattice = .true.
        if (nproc_fft .lt. nprocrep .or. do_hartree) commensurate_lattice = .false.

#ifdef PARALLEL
        if (rank .eq. 0) write (6, *) "sub_comm_fft uses", nproc_fft, "processors"
        call mpi_sub_comm_create(commrep_mpi, nproc_fft, sub_comm_fft, ierr)
#endif

        if (ipc .eq. 2 .and. .not. double_overs .and. contracted_on) then
            do kk = 1, iesup_c - 2*nelorbh*molecular
                if (dup_c(2*kk) .ne. 0.d0) then
                    double_overs = .true.
                    if (rank .eq. 0) then
                        write (6, *) ' Warning you should have real contracted if &
                                & you want to use the fast algorithm! '
                    end if
                end if
            end do
#ifdef PARALLEL
! Just to avoid to be inconsistent.
            call mpi_bcast(double_overs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
        end if

        if (do_hartree) then
            if (rank .eq. 0) write (6, *) ' ERROR option do_hartree not implemented yet '
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        ! read input spin grid (lines right after the DFT card)
        if (yeslsda) then
            dimspin_input = nxs*nys*nzs
            allocate (spin_input(nxs, nys, nzs))
            spin_input = 0
            if (rank .eq. 0) then
                write (6, *) ' Read spin grid ', nzs, nys, nxs
                if (dimspin_input .gt. 2) then
                    do k = 1, nzs
                        do j = 1, nys
                            read (5, *) (spin_input(i, j, k), i=1, nxs)
                            write (6, *) (spin_input(i, j, k), i=1, nxs)
                        end do
                        if (k .ne. nzs) then
                            read (5, *)
                            write (6, *)
                        end if
                    end do
                else
                    spin_input = 1
                    if (randspin .lt. 0.d0 .and. dimspin_input .eq. 2) then
                        write (6, *) ' Warning initializing antiferromagnetic spin '
                        if (nzs .eq. 2) spin_input(1, 1, 2) = -1
                        if (nxs .eq. 2) spin_input(2, 1, 1) = -1
                        if (nys .eq. 2) spin_input(1, 2, 1) = -1
                    end if
                end if
            end if
#ifdef PARALLEL
            call mpi_bcast(spin_input, dimspin_input, MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr)
#endif
        end if
        if (h_charge .ne. 0.d0) then
            write (6, *) ' Charge grid input'
            dimspin_input = nxs*nys*nzs
            allocate (charge_input(nxs, nys, nzs))
            charge_input = 0
            if (rank .eq. 0) then
                write (6, *) ' Read charge grid ', nzs, nys, nxs
                if (dimspin_input .gt. 2) then
                    if (yeslsda) read (5, *)
                    do k = 1, nzs
                        do j = 1, nys
                            read (5, *) (charge_input(i, j, k), i=1, nxs)
                            write (6, *) (charge_input(i, j, k), i=1, nxs)
                        end do
                        if (k .ne. nzs) then
                            read (5, *)
                            write (6, *)
                        end if
                    end do
                else
                    charge_input = 1
                    if (randspin .lt. 0.d0 .and. dimspin_input .eq. 2) then
                        write (6, *) ' Warning initializing CDW commensurate '
                        if (nzs .eq. 2) charge_input(1, 1, 2) = -1
                        if (nxs .eq. 2) charge_input(2, 1, 1) = -1
                        if (nys .eq. 2) charge_input(1, 2, 1) = -1
                    end if
                end if
            end if
#ifdef PARALLEL
            call mpi_bcast(charge_input, dimspin_input, MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr)
#endif
        end if
        ! initialize Ewalds sums
        if (iespbc) then
            !      nel_neutral=sum(zetar(:)) !  for tests
            nel_neutral = nel
            call InitEwald(nion, zetar, nel_neutral, nws, thread_active)
            if (rank .eq. 0) then
                write (*, *) ' Number of G vectors ', n_gvec
                write (*, *) ' Ewald Self Energy ', eself
            end if
            kmax = n_gvec
            kmax2 = 2*kmax
            call EwaldSum1b(rion) ! Initialize one body ewald
        else
            ewaldion1b = 0.d0
            ewaldel1b = 0.d0
            n_gvec = 1
            kmax = 0
            kmax2 = 0
            at = 0.d0
            at(1, 1) = 1.d0
            at(2, 2) = 1.d0
            at(3, 3) = 1.d0
        end if
        !
        ! write information on the calculation
        !
        call write_info
        ! syncronize all processors
#ifdef  PARALLEL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier
        !
        ! define scratch vector psip/ipsip maximum dimension
        !
        iscramax = 3*nion*nel*npsamax*nintpseudo + 2*npseudopar + indt*nel
        iscramax = max(3*nion*nion + nel, 500, 10*nelorb, nelup*nelup, iscramax)
        iscramax = max(iscramax, nelorbh*(indt + 6) + 20*(indt + 1)) ! upnewwf
        if (jaccond .gt. 0.d0) then
            iscramax = max(iscramax, 5*maxold)
            lworkjac = 5*maxold
        end if
        if (yeslsda) then
            iscramax = max(iscramax, 2*bands)
        else
            iscramax = max(iscramax, bands)
        end if
        iscramax = max(iscramax, nelorb)
        !
        iscraipsip = max(nelorb, 2*maxit - 1)
        iscraipsip = max(iscraipsip, 2*nnozero_c + 1)
        iscraipsip = max(iscraipsip, 2*nnozeroj_c + 1, 9*nion + 1, iesupr + 1)
        if (yeslsda) then
            iscraipsip = max(iscraipsip, 2*bands)
        else
            iscraipsip = max(iscraipsip, bands)
        end if
        !
        ! initialize dimensions of distributed matrices (hamiltl,oversl)
        !
#ifdef __SCALAPACK
        if (contracted_on) then
            call descla_init(descla, nelorb_at, nelorb_at, np_ortho, me_ortho, ortho_comm, ortho_comm_id)
            if (descla(lambda_node_) > 0) then
                call descinit(desch, nelorb_at, nelorb_at, descla(nlax_), descla(nlax_), 0, 0, ortho_cntx, descla(nlax_), info)
                if (info /= 0) call errore(' prep  ', ' descinit ', abs(info))
            end if
        else
            call descla_init(descla, nelorb, nelorb, np_ortho, me_ortho, ortho_comm, ortho_comm_id)
            if (descla(lambda_node_) > 0) then
                call descinit(desch, nelorb, nelorb, descla(nlax_), descla(nlax_), 0, 0, ortho_cntx, descla(nlax_), info)
                if (info /= 0) call errore(' prep  ', ' descinit ', abs(info))
            end if
        end if
#endif
        vh_test = -1.d0 ! Initialization vh_test
        if (scale_hartree .gt. 0) then
            scale_hartreen = scale_hartree
        else
            scale_hartreen = 1.d0
        end if
        if (enforce_detailb) enforce_detailb = .false. ! option true is not considered in DFT

        return

        !-----------------------------
    end subroutine Initializeall
    !-----------------------------

    !---------------------------
    subroutine check_dft_input
        !---------------------------

        implicit none
        real(8), external :: dlamch
        !
        ! check all default values of input DFT
        ! parameters in the master and then broadcast to
        ! all processors.
        !
        if (nelup .lt. nel - nelup) &
            call error(' setup_calculation ', ' It is assumed that nelup > neldo ', 1, rank)

        if (rank .eq. 0) then
            ! initial dimension of the spin grid
            if (nxs .eq. -1) then
                nxs = 1
            end if
            if (nys .eq. -1) then
                nys = 1
            end if
            if (nzs .eq. -1) then
                nzs = 1
            end if
            ! number of bands
            if (nelocc .ne. 0 .and. bands .lt. nelocc) then
                bands = nelocc
                write (6, *) ' Default value of bands =', bands
            end if
            if (bands .gt. nelorbh .and. .not. contracted_on) then
                bands = nelorbh
                write (6, *) ' Warning # bands < =', bands
            end if
            if (bands .gt. nelorb_at .and. contracted_on) then
                bands = nelorb_at
                write (6, *) ' Warning # bands < =', bands
            end if
            ! machine precision
            if (eps_mach .eq. -1.d0) then
                eps_mach = DLAMCH('E')
                write (6, *) ' Default value for relative machine precision =', eps_mach
            end if
            ! mixing parameter
            if (mixing .eq. -1.d0) then
                if (typeopt .eq. 2) then
                    mixing = 0.25d0
                elseif (typeopt .eq. 0) then
                    mixing = 0.25d0
                elseif (typeopt .eq. 4) then
                    mixing = 1.d0
                else
                    mixing = 0.01d0
                end if
                write (6, *) ' Default value of mixing =', mixing
            else
                write (6, *) ' mixing used =', mixing
            end if
            ! exchange-correlation functional type
            if (typedft .eq. 1) then
                write (6, *) ' Default value of typedft  =', typedft
            else
                write (6, *) ' typedft used  =', typedft
            end if
            !
            if (mixingder .eq. -1.d0) then
                if (mod(typeopt, 2) .eq. 0) then
                    mixingder = 0.05d0
                else
                    mixingder = 100.d0
                end if
                write (6, *) ' Default value of mixingder =', mixingder
            else
                write (6, *) ' mixingder used =', mixingder
            end if
            write (6, *) ' Warning epsover used =', epsover

            if (mixingder .lt. 0) then
                changemix = .true.
                mixingder = -mixingder
                write (6, *) ' Warning adaptive changing mixingder !!! '
            else
                changemix = .false.
            end if
            !
            ! definition of the flag "yeslsda":
            ! exchange-correlation functionals for LSDA calculations.
            ! NB. If the variable typedft is negative, then the code
            ! perform a simple check of the XC fuctional.
            !
            if ((abs(typedft) .eq. 4 .or. abs(typedft) .eq. 5) .and. (.not. symmagp .or. ipc .eq. 2)) then
                yeslsda = .true.
            else
                yeslsda = .false.
                if ((abs(typedft) .eq. 4 .or. abs(typedft) .eq. 5))&
                        & call error(' setup_calculation ', &
                        & 'To use LSDA functional, start with a non-symmetric AGP. typedft&
                                & changed to corresponding LDA type ', -1, rank)
                if (typedft .eq. 4) typedft = 1
                if (typedft .eq. -4) typedft = -1
                if (typedft .eq. 5) typedft = 3
                if (typedft .eq. -5) typedft = -3
            end if

            write (6, *) ' maxold used =', maxold ! # of previous iteration densities saved
            if (mod(typeopt, 2) .eq. 0 .and. epssr .ne. 0.d0) then
                epssr = 0.d0
                write (6, *) ' Warning epssr set to zero in this case '
            elseif (epssr .eq. -1.d0) then
                epssr = -1.d0
                write (6, *) ' Default value of epssr =', epssr
            end if
            if (weightvh .ne. 1.d0) then
                write (6, *) ' Warning assuming one body Ewald sums'
            end if
            !
            if (optocc .ne. 0) then
                write (6, *) ' Smearing used epsshell=', epsshell
            else
                epsshell = 0.0d0 ! no smearing with fixed occupations
            end if
            !
            if (write_den) write (6, *) 'Warning write density in Xcrysden format'
            if (contracted_on) write (6, *) 'Warning using contracted basis'
            ! some check for k-points
            !       if(.not.kaverage.and.yeswrite10) then
            !          write(6,*) ' Warning: cannot write many fort.10 if no k-points sampling!! '
            !          yeswrite10=.false.
            !       endif
            if (manyfort10 .and. .not. yeswrite10) &
                write (6, *) ' Warning: writing only fort.10 related to first k-point '
            if (.not. manyfort10 .and. (yeswrite10 .or. yesread10)) then
                yeswrite10 = .false.
                yesread10 = .false.
                write (6, *) ' Warning yeswrite10/yesread10=.false. in this case, forced '
            end if

            ! some initialization for reading scratch files
            if (iopt .ne. 1 .and. writescratch .eq. 0) then
                noread = .false.
                yesreadf = .true.
            else
                noread = .true.
                yesreadf = .false.
            end if
            if (iopt .eq. 3) then
                yesreadf = .false.
                iopt = 0
            end if
            ! set default values for band structure calculations
            if (compute_bands) then
                iopt = 1 ! density is read from file in this case
                typeopt = 0 ! no mixing needed
                noread = .true. ! recompute overlaps form scratch since k-points are different
                print_energies = .false. ! do not print out energies for each k-point
                decoupled_run = .false. ! no meaning in band calculation
                kaverage = .true. ! need the k-points card
            end if
            if (.not. kaverage .and. yeswrite10 .and. nbead .gt. 1) decoupled_run = .true.

            ! some checks on decoupled k-points calculations
            if (decoupled_run .and. (.not. kaverage .and. .not. yeswrite10)) decoupled_run = .false.
            if (decoupled_run) print_energies = .true.
            if (fix_density .and. .not. decoupled_run) then
                fix_density = .false.
                write (6, *) ' Warning fix_density=.true. option can be used only with decoupled_run =.true., &
                        & forced fix_density= .false.'
            end if
            if (.not. kaverage .and. yeswrite10) fix_density = .false.

        end if ! endif rank.eq.0
        if (rank .eq. 0) write (6, *) ' decoupled run =', decoupled_run

        return

    end subroutine check_dft_input

    !----------------------
    subroutine write_info
        !----------------------
        !
        ! write information on the simulation before
        ! starting the SC cycle.
        !
        use allio, only: iespbc, yes_crystal
        implicit none
        character(len=5) :: bas_type
        character(len=30) :: calc_type

        if (iespbc) then
            if (yes_crystal) then
                bas_type = 'PBC_C'
            else
                bas_type = 'PBC'
            end if
        else
            bas_type = 'OPEN'
        end if

        calc_type = 'single phase'
        if (.not. decoupled_run) then
            if (kaverage) calc_type = 'kaverage'
        elseif (decoupled_run) then
            calc_type = 'independent kaverage'
        elseif (compute_bands) then
            calc_type = 'band structure'
        end if

        if (rank .eq. 0) then

            write (6, *)
            write (6, '(a)') '   CALCULATION INFORMATION '
            write (6, '(a)') ' # Hartree Atomic Units used'
            write (6, '(a,a)') ' # computation type  = ', calc_type
            write (6, '(a,I5)') ' # iopt              = ', iopt
            write (6, '(a,L2)') ' # memlarge          = ', memlarge
            write (6, '(a,a)') ' # basis type        = ', bas_type
            write (6, '(a,I10)') ' # basis dim         = ', nelorb
            write (6, '(a,I10)') ' # Number MOs        = ', nmol
            write (6, '(a,L2)') ' # contracted_on     = ', contracted_on
            write (6, '(a,I5)') ' # One-body Jastrow  = ', n_body_on
            if (yes_sparse) write (6, '(a)') ' Warning using SPARSE jastrow matrix in writing'
            if (yes_crystal) &
                write (6, '(a,F8.5)') ' # basis cutoff      = ', epsbas

            write (6, '(a,I5)') ' # typedft           = ', typedft
            write (6, '(a,I5)') ' # typeopt           = ', typeopt
            write (6, '(a,L2)') ' # spin calculation  = ', yeslsda

            if (yeslsda) &
                write (6, '(a,F8.5)') ' # external H field  = ', h_field
            write (6, '(a,F8.5)') ' # external H charge field  = ', h_charge

            write (6, '(a,I5)') ' # bands             = ', bands
            write (6, '(a,F9.5)') ' # mixing            = ', mixing
            write (6, '(a,I5)') ' # optocc            = ', optocc
            write (6, '(a,F12.8)') ' # smearing          = ', epsshell
            write (6, '(a,I5)') ' # k-points          = ', nk
            write (6, '(a,L2)') ' # Optimize overlaps   = ', optimize_overs
            write (6, *)

        end if

#ifdef __SCALAPACK
        call error(' prep ', ' Using SCALAPACK algorithm ', -1, rank)
#endif

        return

    end subroutine write_info

    subroutine set_addresses

        implicit none

        if (contracted_on) then
            nelorbu = nelorb_at
            if (yeslsda .or. ipc .eq. 2) then
                nelorb3 = 3*nelorb_at
            else
                nelorb3 = 2*nelorb_at
            end if
        else
            nelorbu = nelorb
            if (yeslsda .or. ipc .eq. 2) then
                nelorb3 = 3*nelorb
            else
                nelorb3 = 2*nelorb
            end if
        end if
        !
        if (yeslsda .or. ipc .eq. 2) then
            nelorbudo = nelorbu
            nelorbdo = nelorb
            bandsdo = bands
        else
            nelorbudo = 1
            nelorbdo = 1
            bandsdo = 1
        end if
        !
#ifdef __SCALAPACK
        nlax = descla(nlax_)
#else
        nlax = 1
#endif
        !
        ! initialize k-point of the current pool once for all
        !
        indk = rankcolrep + 1
        ikpoint = indk
        if (manyfort10) then
            phase(:) = xkp(:, indk)
            phase_down(:) = xkp_down(:, indk)
            phase2pi(:) = phase(:)*TWO_PI
            phase2pi_down(:) = phase_down(:)*TWO_PI
        end if
        !
        ! if phase up = +/- phase down not necessary to compute
        ! twice the overlap matrices. The following relations holds:
        ! same_phase = T -> oversdo = overs
        ! same_phase = F -> oversdo = dconjg(overs)
        ! The same for Ham overlaps and Ham matrix. This allows a decrease
        ! of the computational time of approximately a factor 2.
        !
        if (ipc .eq. 2) then

            !      is_same_up(:)   = ( ( abs( phase(:)-phase_down(:) ) ) < 1.0d-6 )
            !      is_same_down(:) = ( ( abs( phase(:)+phase_down(:) ) ) < 1.0d-6 )

            !      double_overs   = .true.
            !      same_phase     = .false.
            if (.not. same_phase .and. .not. opposite_phase) then
                double_overs = .true.
            else
                double_overs = .false.
            end if

        else

            double_overs = .false.
            !      same_phase = .true.

        end if

        if (optimize_overs .and. double_overs) then
            if (rank .eq. 0) write (6, *) ' Warning overlap optimization is not &
                    &possible with two arbitrarily different phases '
            optimize_overs = .false.
        end if

        return

    end subroutine set_addresses

end module setup
