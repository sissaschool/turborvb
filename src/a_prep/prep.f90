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

program main

    use setup
    use freeelmod_complex, only: self_consistent_run
    use parallel_module, only: old_threads, setup_para
    implicit none
    integer, external :: omp_get_max_threads
    character(lchlen) :: path
    real(8), external :: cclock
#ifdef PARALLEL
    include 'mpif.h'
#else
    integer status(1)
    character(20) str
    character(100) name_tool
    !
    ! open documentation if required
    !
    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then
        !          Input the name of the file exactly as it is in /doc
        name_tool = 'prep'
        call help_online(name_tool)
        stop
    end if
#endif
    !
    ! Initialize MPI
    !
#if defined PARALLEL
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) write (6, *) ' Number of mpi proc =', nproc
#else
    rank = 0
    nproc = 1
#endif
    !
    ! Initialize OpenMP
    !
#if defined _OPENMP
    old_threads = omp_get_max_threads()
#else
    old_threads = 1
#endif

    ! output version information
    if (rank .eq. 0) call print_version
    
    ! print threads and mpi info.
    if (rank .eq. 0) write (6, *) ' Number of threads /mpi proc =', old_threads

#if defined UNREL_SMP
#if defined PARALLEL
    call mpi_bcast(old_threads, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
#if defined _OPENMP
    call omp_set_num_threads(old_threads) ! force the same
#endif
#endif
    if (rank .eq. 0) write (6, *) ' Warning init. value of threads/mpi task', old_threads
    !
    call get_dir(path)
    if (rank .eq. 0) write (6, *) ' Initial path : ', path
    !
    ! Initialize error flags needed by read_datasmin
    yesdft = .true. ! important flag to avoid useless allocation and checks
    iflagerrall = 0
    iflagerr = 0
    !
    ! reading main input cards from std input, in common with QMC calculations.
    ! One needs to do it at first in order to gather necessary information
    ! on parallelization.
    !     Cards read: SIMULATIONS,PARAMETERS,PSEUDO,VMC,K-POINTS.
    call read_datasmin
    yesfast = 1
    if (itestr4 .eq. -8) then
        call error(' prep ', ' Optimization of Z is not possible with DFT ', -1, rank)
        yeszagp = .false.
        itestr4 = -4
    end if
    !
    ! initialize parallel environment at various levels:
    ! 1) processors pools in case of k-points sampling (flag kaverage)
    ! 2) SCALAPACK block matrices within the pools (compiling option __SCALAPACK)
    call setup_para
    !
    ! Read all necessary files to setup the DFT calculation:
    ! 1) read pseudopotential if any.
    ! 2) read wave function fort.10/fort.10_new
    ! 3) read DFT and MOLECULAR cards on std input.
    ! 4) generate scratch files for continuation.
    call Initializeall
    !
    ! deallocate useless variables allocated for QMC
    ! small reallocation just to be safe
    !
    deallocate (tabpip, table, tabler, ainv, winvdo, winvup, winvj, &
                winvbar, winvjbar, wint, diagfn, jastrowall_ee, &
                jastrowall_ei, jasnew_ee, jasnew_ei, winvjbarsz)
    allocate (tabpip(1), table(1), tabler(1), ainv(1), winvdo(1), winvup(1), &
              winvj(1), winvbar(1), winvjbar(1), wint(1), diagfn(1), &
              jastrowall_ee(1, 1, 1, 1), jastrowall_ei(1, 1, 1), jasnew_ee(1), &
              jasnew_ei(1), winvjbarsz(1), dup(max(ipc*iesup, 1)))
    tabpip = 0.d0
    table = 0.d0
    tabler = 0.d0
    ainv = 0.d0
    winvdo = 0.d0
    winvup = 0.d0
    winvj = 0.d0
    winvbar = 0.d0
    winvjbar = 0.d0
    wint = 0.d0
    diagfn = 0.d0
    jastrowall_ee = 0.d0
    jastrowall_ei = 0.d0
    jasnew_ee = 0.d0
    jasnew_ei = 0.d0
    winvjbarsz = 0.d0
    dup = 0.d0
    deallocate (detmat, nozero, nozerodet, scale)
    allocate (detmat(nelorbh), nozero(1), scale(1), nozerodet(1))
    detmat = 0.d0
    nozero = 0
    scale = 0
    nozerodet = 0
    if (contraction .gt. 0 .and. .not. contracted_on) then
        deallocate (mu_c)
        allocate (mu_c(1, 1))
        mu_c = 0.d0
    end if
    if (iespbc) then
        deallocate (sum_q_cos_gr, sum_q_sin_gr)
        allocate (sum_q_cos_gr(1, 1), sum_q_sin_gr(1, 1))
        sum_q_cos_gr = 0.d0
        sum_q_sin_gr = 0.d0
    end if
    deallocate (winv, psip, ipsip)
    !
    ! allocate global scratch vectors needed by DFT
    ! winv --> real scratch vector
    !
    if (.not. yes_complex) then
        allocate (winv(nelorb*(indt + 5)), psip(iscramax), ipsip(iscraipsip))
        winv = 0.d0
        psip = 0.d0
        ipsip = 0
    else
        call error(' prep ', ' Double allocation for complex algorithm ', -1, rank)
        allocate (psip(2*iscramax), ipsip(iscraipsip), winv(1))
        psip = 0.d0
        ipsip = 0
        winv = 0.d0
    end if
    !
    ! initialize all variables needed by self and non self-consistent
    ! calculations.
    time_total = cclock()
    !
    call initialize_environment()
    !
    ! update molecular orbitals using one-body Hamiltonian
    ! with self-consistent DFT cycle or perform a
    ! non self-consistent calculation to compute band structure
    ! and other electronic properties.
    !
    if (.not. compute_bands) then
        call self_consistent_run()
    else
        call non_self_consistent_run()
    end if
    !
    ! deallocate all variables and write the final
    ! outputs and wave functions.
    !
    time_total = cclock() - time_total

    call write_output_and_finalize()

#if defined __SCALAPACK
    call BLACS_GRIDEXIT(ortho_cntx)
    call mpi_comm_free(ortho_comm, ierr)
#endif
#ifdef PARALLEL
    if (manyfort10 .and. .not. compute_bands) then
        call mpi_comm_free(commrep_mpi, ierr)
        call mpi_comm_free(commcolrep_mpi, ierr)
    end if
    call mpi_finalize(ierr)
#endif
    stop

end program main
