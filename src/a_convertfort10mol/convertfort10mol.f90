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
!
! Sandro Sorella created on 25th Nov. 2009.

!> @brief Main program for converting fort.10 files and computing molecular orbitals
!> @details This program reads a fort.10 input file and performs molecular orbital
!>          conversion and optimization. It supports both serial and parallel execution
!>          using MPI, handles pseudopotential calculations, and generates optimized
!>          molecular orbitals for quantum Monte Carlo calculations.
!>          
!>          The program processes three main namelist inputs:
!>          - control: Optimization parameters and algorithm settings
!>          - mesh_info: 3D mesh grid parameters for orbital evaluation
!>          - molec_info: Molecular orbital specifications and constraints
!>          
!>          Key features:
!>          - Parallel processing with MPI support
!>          - Pseudopotential handling
!>          - Molecular orbital optimization
!>          - Periodic boundary condition support
!>          - Output generation for TurboRVB calculations
program convertfort10mol
    use convertmod
    use allio
    use sub_comm
    implicit none
    logical yesbig
    integer rankn, nprocn, ithread, comm_mpi, j, i
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    !> @brief Control parameters for molecular orbital optimization
    !> @details This namelist contains parameters that control the molecular orbital
    !>          optimization process, including convergence criteria, algorithm settings,
    !>          and output options.
    namelist /control/ epsdgm, molopt, weight_loc, power, orthoyes, yesbig, membig, gramyes &
        , epsbas, allowed_averagek, add_onebody2det, only_molecular, add_offmol
    
    !> @brief 3D mesh grid parameters for orbital evaluation
    !> @details This namelist defines the 3D mesh grid used for evaluating molecular
    !>          orbitals and computing overlap matrices. It includes grid dimensions,
    !>          spacing parameters, and origin shift options.
    namelist /mesh_info/ nbufd, nx, ny, nz, ax, ay, az, shift_origin, shiftx, shifty, shiftz
    
    !> @brief Molecular orbital specifications and constraints
    !> @details This namelist specifies the number and range of molecular orbitals
    !>          to be computed, including minimum and maximum orbital counts and
    !>          overlap printing options.
    namelist /molec_info/ nmol, nmolmin, nmolmax, printoverlap

#ifdef PARALLEL
    !> @brief Initialize MPI parallel processing environment
    !> @details This section initializes the MPI environment for parallel execution,
    !>          sets up communicators, and defines rank and size variables for
    !>          distributed computation across multiple processes.
    include 'mpif.h'
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nprocn, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rankn, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    comm_mpi = MPI_COMM_WORLD
    ! define also these values for the read_fort10
    rankrep = rank
    commrep_mpi = MPI_COMM_WORLD
    rankcolrep = 0
    commcolrep_mpi = MPI_COMM_WORLD
#else
    !> @brief Serial execution setup
    !> @details This section sets up variables for serial execution when MPI is not
    !>          available. It also handles command-line help options for the program.
    rankn = 0
    nprocn = 1
    comm_mpi = 0
    rankrep = 0
    commrep_mpi = 0
    rankcolrep = 0
    commcolrep_mpi = 0
    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        ! Input the name of the file exactly as it is in /doc
        name_tool = 'convertfort10mol'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added
#endif
    ! output version information
    if (rankn .eq. 0) call print_version

    nprocu = nprocn
    nproc_diag = nprocu
    call mpi_sub_comm_create(comm_mpi, nproc_diag, sub_comm_diag, ierr)

    if (rankn .eq. 0) then

        write (*, *) " * * * READ fort.10 and COMPUTE molecular orbitals * * * "

        open (unit=10, file='fort.10_in', form='formatted', status='unknown')

    end if

    rank = rankn
    nproc = nprocn
    call default_allocate ! default allocation contained in memOP.f90
    !  rankopt=rank
    !  nprocopt=nproc
    !#ifdef PARALLEL
    !  commpot_mpi=MPI_COMM_WORLD
    !#endif
    if (rank .eq. 0) then
        !> @brief Initialize default parameters for molecular orbital optimization
        !> @details This section sets up default values for optimization parameters
        !>          including convergence criteria, algorithm flags, and processing options.
        yesmin = 0 ! no optimization is employed here
        yesbig = .false.
        membig = .true.
        molopt = 0
        epsdgm = 1d-14 ! the sdv criterium  machine precision cond. number
        power = 1.d0
        weight_loc = -1.d0
        orthoyes = .true.
        gramyes = .false.
        iflagerr = 1
        add_onebody2det = .true.
        only_molecular = .false.
        add_offmol = .false.
        read (5, nml=control, err=115)
        iflagerr = 0
115     if (iflagerr .ne. 0) write (6, *) ' ERROR reading control !!! '

        !> @brief Set up mesh parameters based on system size and optimization settings
        !> @details This section determines appropriate mesh parameters for orbital
        !>          evaluation, including buffer sizes and grid dimensions based on
        !>          the number of orbitals and optimization criteria.
        if (epsdgm .ge. 0.d0) then
            if (nelorb .lt. 2000) then
                nbufd = 1000
            else
                nbufd = 100
            end if
        else
            nbufd = 1
        end if
        nx = 0
        ny = 0
        nz = 0
        ax = 0.1d0
        ay = 0.d0
        az = 0.d0
        shift_origin = .true.
        shiftx = .false.
        shifty = .false.
        shiftz = .false.

        if (iflagerr .eq. 0) then
            iflagerr = 1
            read (5, nml=mesh_info, err=116)
            iflagerr = 0
116         if (iflagerr .ne. 0) write (6, *) ' ERROR reading mesh_info !!! '
            if (nx .eq. 0) then
                write (6, *) ' The number of mesh point is zero !!! ', nx, ny, nz
                stop
            end if
            if (ny .eq. 0) ny = nx
            if (nz .eq. 0) nz = ny
        end if
    end if

    call checkiflagerr(iflagerr, rank, 'ERROR reading')
#ifdef PARALLEL
    call mpi_bcast(molopt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbufd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ay, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(az, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(power, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(weight_loc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsdgm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsbas, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesbig, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(membig, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(orthoyes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(gramyes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(add_onebody2det, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(only_molecular, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(add_offmol, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shift_origin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shiftx, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shifty, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shiftz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
    if (add_offmol .and. .not. only_molecular) then
        if (rank .eq. 0) write (6, *) ' Warning add_offmol works only with only_molecular on,  set to .true.  '
        only_molecular = .true.
    end if

    if (yesbig) yesfast = 0

    !> @brief Check for pseudopotential files and initialize pseudopotential handling
    !> @details This section performs a fast read of fort.10 to check for pseudopotential
    !>          data and sets up the pseudopotential processing environment if needed.
    ! read_fast to check if there are pseudo
    if (rank .eq. 0) call read_fort10_fast
#ifdef PARALLEL
    call mpi_bcast(npsar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
    npsa = npsar
    pseudofile = "pseudo.dat"
    iflagerrall = 0

    if (npsa .gt. 0) nintpsa = 6

    !> @brief Read pseudopotential data and wave function information
    !> @details This section reads pseudopotential files and the main fort.10 input
    !>          file containing wave function data, orbital information, and system
    !>          parameters for the molecular orbital conversion process.
    ! ----------------------
    ! reading pseudo and w.f.
    call read_pseudo
    call read_fort10(10)
    ! ----------------------
    if (allowed_averagek .and. rank .eq. 0) write (6, *) ' Warning attaching flux to det !!! '

    !> @brief Display and validate determinant matrix information
    !> @details This section prints the determinant matrix for small systems (nelorb_c < 100)
    !>          to verify the input data and check matrix properties for both contracted
    !>          and uncontracted basis sets.
    if (rank .eq. 0 .and. nelorb_c .lt. 100) then
        write (6, *) ' Read matrix detmat_c ', nelorb_c, nelorbh, sum(abs(detmat_c(:)))
        if (contraction .eq. 0) then
            do i = 1, nelorb_c
                do j = i, nelorb_c
                    if (ipc .eq. 1) then
                        write (6, *) i, j, detmat(nelorb_c*(j - 1) + i)
                    else
                        write (6, *) i, j, detmat(ipc*nelorb_c*(j - 1) + 2*i - 1), &
                            detmat(ipc*nelorb_c*(j - 1) + 2*i)
                    end if
                end do
            end do
        else
            do i = 1, nelorb_c
                do j = i, nelorb_c
                    if (ipc .eq. 1) then
                        write (6, *) i, j, detmat_c(nelorb_c*(j - 1) + i)
                    else
                        write (6, *) i, j, detmat_c(ipc*nelorb_c*(j - 1) + 2*i - 1), &
                            detmat_c(ipc*nelorb_c*(j - 1) + 2*i)
                    end if
                    !write(6,*) i,j,detmat_c(ipc*nelorb_c*(j-1)+2*i-1:ipc*nelorb_c*(j-1)+2*i)
                end do
            end do
        end if
    end if

    !> @brief Recompute determinant matrix if needed for uncontracted basis
    !> @details This section checks if the determinant matrix needs to be recomputed
    !>          for uncontracted basis sets and performs the matrix contraction
    !>          operation if necessary.
    if (nelorb_c .ge. nelorbh*ipf .and. yesfast .ne. 0) then ! nelorb_c=nelorbh in the case of uncontracted basis

        if (iscramax .le. nelorbh*ipf*nelcol_c) then
            deallocate (psip)
            iscramax = nelorbh*ipf*nelcol_c
            allocate (psip(iscramax))
        end if
        if (rank .eq. 0) write (6, *) ' Warning recomputing detmat '
        yesfast = 0
        if (allocated(detmat)) deallocate (detmat)
        allocate (detmat(nelorbh*ipf*nelcolh))
        yesdetmat = .true.
        detmat = 0.d0
        call scontract_mat_det(nelorbh, nelorbh, nelcolh&
                &, nelorb_c, nelcol_c, detmat, detmat_c, mu_c, psip)
    end if

    !> @brief Set up mesh parameters for periodic boundary conditions or free space
    !> @details This section determines the appropriate mesh spacing based on whether
    !>          periodic boundary conditions are used (cellscale) or free space
    !>          calculations (user-defined ax, ay, az parameters).
    if (iespbc) then
        ax = cellscale(1)/nx
        ay = cellscale(2)/ny
        az = cellscale(3)/nz
        if (rank .eq. 0) write (6, *) ' lattice mesh chosen ', ax, ay, az, Lbox
    else
        if (rank .eq. 0) write (6, *) ' lattice mesh ax,ay,az '
        if (ay .eq. 0.d0) ay = ax
        if (az .eq. 0.d0) az = ay
        if (rank .eq. 0) write (6, *) ax, ay, az
    end if
    !
    if (rank .eq. 0) then
        !> @brief Initialize molecular orbital parameters and constraints
        !> @details This section sets up the molecular orbital specifications including
        !>          the number of orbitals to compute, minimum and maximum orbital
        !>          counts, and overlap printing options.
        nmol = 0
        nmolmin = 0
        nmolmax = 0
        parr = epsdgm
        printoverlap = .false.
        iflagerr = 1
        read (5, nml=molec_info, err=117)
        iflagerr = 0
117     continue
        if (nmol .eq. 0 .or. epsdgm .eq. -1.d0) nmol = molecular - nelup + neldo
        if (nmolmax .eq. 0) nmolmax = nmol
        if (nmol .lt. nmolmax) then
            nmol = nmolmax
            write (6, *) ' Warning nmol>= nmolmax, changed nmol=nmolmax', nmol
        end if
        if (nmolmin .eq. 0 .and. nmol .gt. 0) nmolmin = 1
    end if
    call checkiflagerr(iflagerr, rank, 'ERROR reading molec_info')
#ifdef PARALLEL
    call mpi_bcast(nmol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nmolmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nmolmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(parr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(printoverlap, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    if (nmol .lt. neldo) then
        if (rank .eq. 0) write (6, *) ' Warning nmol>=neldo , changed to ', neldo
        nmol = neldo
        nmolmax = nmol
        nmolmin = nmol
    end if
    if (nmolmax .lt. neldo) then
        if (rank .eq. 0) write (6, *) ' Warning nmolmax>=neldo , changed to ', neldo
        nmolmax = neldo
        nmolmin = neldo
    end if
    if (nmolmin .eq. 0) then
        if (rank .eq. 0) write (6, *) ' Default value of nmolmin= ', neldo
        nmolmin = neldo
    end if
    nmolmaxw = nmolmax

    if (rank .eq. 0) write (6, *) ' Chosen nmolmin nmolmax =', nmolmin, nmolmax
    !> @brief Perform molecular orbital conversion and optimization
    !> @details This section executes the main molecular orbital conversion process
    !>          using the convertmol_fast subroutine. The computation is wrapped in
    !>          OpenMP target data directives for GPU acceleration if available.
#ifdef _OFFLOAD
!$omp target data map(to:mu_c)
#endif
    call convertmol_fast
#ifdef _OFFLOAD
!$omp end target data
#endif

    !> @brief Generate output files with optimized molecular orbitals
    !> @details This section creates the output fort.10_new file containing the
    !>          optimized molecular orbitals and updated wave function parameters
    !>          for use in subsequent TurboRVB calculations.
    if (rank .eq. 0) then
        close (10)
        open (unit=10, file='fort.10_new', form='formatted', status='unknown')
        contraction = 1 ! contraction assumed in any case
        call write_fort10(10)
        close (10)
    end if

#ifdef PARALLEL
    call mpi_finalize(ierr)
#endif

end program convertfort10mol
