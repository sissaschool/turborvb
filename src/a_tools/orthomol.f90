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

program orthomol
    use convertmod
    use allio
    implicit none
    logical yesbig
    integer rankn, nprocn, ithread, i
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    namelist /mesh_info/ nbufd, nx, ny, nz, ax, ay, az, shift_origin, shiftx, shifty, shiftz
    namelist /molec_info/ nummol

#ifdef PARALLEL
    include 'mpif.h'
    call mpi_init(ierr)
!     call mpi_init_thread(MPI_THREAD_FUNNELED,ithread,ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nprocn, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rankn, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    nprocu = nprocn
    ! define also these values for the read_fort10
    rankrep = rank
    commrep_mpi = MPI_COMM_WORLD
    rankcolrep = 0
    commcolrep_mpi = MPI_COMM_WORLD

#else
    ! define also these values for the read_fort10
    rankrep = 0
    commrep_mpi = 0
    rankcolrep = 0
    commcolrep_mpi = 0
    rankn = 0
    nprocn = 1
    nprocu = 1

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'orthomol'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added
#endif

    if (rankn .eq. 0) then
        write (*, *) " * * * READ fort.10 * * * "

        open (unit=10, file='fort.10', form='formatted', status='unknown')

    end if

    rank = rankn
    nproc = nprocn
    !     rankopt=rank
    !     nprocopt=nproc
    call default_allocate

    if (rank .eq. 0) then
        if (nelorb .lt. 2000) then
            nbufd = 1000
        else
            nbufd = 100
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
        read (5, nml=mesh_info)

        if (nx .eq. 0) then
            write (6, *) ' The number of mesh point is zero !!! ', nx, ny, nz
            stop
        end if

        if (ny .eq. 0) ny = nx
        if (nz .eq. 0) nz = ny

    end if

#ifdef PARALLEL
    call mpi_bcast(nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbufd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ay, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(az, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(add_onebody2det, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shift_origin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shiftx, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shifty, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shiftz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    !     read if there are pseudo
    call read_fort10_fast
    npsa = npsar

    pseudofile = "pseudo.dat"
    iflagerrall = 0

    if (npsa .gt. 0) nintpsa = 6

    call read_pseudo

    call read_fort10(10)

    if (iespbc) then
        ax = cellscale(1)/nx
        ay = cellscale(2)/ny
        az = cellscale(3)/nz
        if (rank .eq. 0) write (6, *) ' lattice mesh chosen ', ax, ay, az, Lbox
    else
        if (rank .eq. 0) write (6, *) ' lattice mesh ax,ay,az '
        if (ay .eq. 0.d0) ay = ax
        if (az .eq. 0.d0) az = ay
        write (6, *) ax, ay, az
    end if

    if (rank .eq. 0) then

        read (5, nml=molec_info)

        if (nummol .ne. 0) then
            allocate (wheremol(abs(nummol)))
            read (5, *) (wheremol(i), i=1, abs(nummol))
        else
            write (6, *) ' Nothing to be orthogonalized '
        end if

        if (nummol .lt. 0) then
            nummol = -nummol
            do i = 1, nummol
                wheremol(i) = wheremol(i) + nelorb_c - molecular
            end do
        end if

    end if

#ifdef PARALLEL
    call mpi_bcast(nummol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (nummol .eq. 0) then
        call mpi_finalize(ierr)
        stop
    end if
    if (rank .ne. 0) allocate (wheremol(nummol))
    call mpi_bcast(wheremol, nummol, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
    if (nummol .eq. 0) stop
#endif

    if (rank .eq. 0) write (6, *) ' Starting orthogonalization '

    call ortho_fast

    if (rank .eq. 0) then
        close (10)
        open (unit=10, file='fort.10_new', form='formatted', status='unknown')

        contraction = 1 ! contracted assumed in any case

        call write_fort10(10)

        close (10)
    end if

#ifdef PARALLEL
    call mpi_finalize(ierr)
#endif

    stop
end
