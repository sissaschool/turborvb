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

program readforward

    !by E. Coccia (22/11/11): read the external potential and vdW parameters
    use extpot, only: ext_pot
    use van_der_waals, only: vdw

    use Constants
    use IO_m
    use allio, only: wherescratch, ngen, nw, npsa &
            &, iflagerr, iflagerrall, yesdft, nel, nelup, nion, zetar_fast, rion_fast, atom_number_fast, zetar&
            &, rion, rionsav, atom_number, pseudorandom, rank, nproc, ieskin, yespress, nwm, nws, in1&
            &, iesbra, kmax2, LBox, iespbc, iend, iendr, yesmin&
            &, ifreqdump, iopt, symmagp, molecular, yesfast, itestr&
            &, nelorb_c, neldo, detmat_c, detmat, iessw, nnozero_c, nozero_c, jbradet, contraction, lastmol&
            &, iesup, nmolfn, ireadmin, nelorbh, nelcol_c, mu_c, firstmol, projm&
            &, celldm, rs, iseedr, ngenc, nbrar, npow, etryr, tbrar, npr&
            &, npbrar, nmatbr, nnozerojr, yesminr, membigr&
            &, ireadr, nmpr, nfatr, np3r, npsovr, nwr&
            &, wcort, wcorwt, nweightr, ngs, iendr, indcor, parbest                   &
            &, enermin, varmin, iesconv, fmax, jmax                        &
            &, parcutr, inext, pippoc, pippo, nprocr, nbinread, stepcg, ncgread         &
            &, epscutur, epstlur, tmes, winvbarfn, winvfn, nintpseudo, nintpsa, indt, epscut &
            &, dists_kel, epscuttype, cellscale, nshell, mult, ipip_adr, dists, yes_complex&
            &, io_level, details_SP, details_DP, membig, rankrep, rankcolrep, commrep_mpi&
            &, commcolrep_mpi, nk, mcol, chara, charaq, wkp, molyes, add_onebody2det&
            &, ndiff, vpotsav_ee
#ifdef _OFFLOAD
    use allio, only: handle, dev_dgetri_workspace, dev_zgetri_workspace, ipsip, psip, jasmat, muj_c &
    &, jasmat_c, eagp_pfaff, winv, winvj, agp, agpn, ainv, winvbar, winvjbar, ainvup&
    &, ainvdo, ldworkspace, lzworkspace, dev_info, dev_dgetrf_workspace, dev_zgetrf_workspace
#endif
    use cell
    use berry_phase

#ifdef PARALLEL
    ! for mpiio by Y. Luo (24/2/15)
    use mpiio, only: mpiio_file_reset_view, mpiio_file_set_zero, mpiio_file_get_disp&
                         &, mpiio_file_create_view, mpiio_file_open, mpiio_file_close
#endif

    use grid_module
    use Ewald, only: InitEwald, kmax, n_gvec, eself, sum_q_cos_gr, sum_q_sin_gr

    use Assar_module
    use Spin2, only: ifspin2, nspin2, spin2_local, inits2pfaff, inits2pfaff_c, ratiospin
    use qpwf_module, only: ifqpwf, ifqpwf_h, ifqpwf_k, ifqpwf_extr, &
                           n_extr_points, extr_point, qpwf_image, &
                           decouple_k, decouple_files
    use Rho_corr_module, only: ifrho_corr
    use Dipole_module, only: ifdipole, ndipole, dipole, quad, nquad, nuclear_dipole, cc_chg, quad_diag

    implicit none
#ifdef _CUSOLVER
    integer :: stat
#endif
    type table_spc
        integer, dimension(:), allocatable :: indd
        real(8), dimension(:), allocatable :: dists
        integer :: nsp
    end type

    integer maxk, k, kt, kkt, i, j, ng, l, nrest, nbin, maxf &
        , nwbuf, nbias, maxj, ibinit, icount, ibin, nmis, lbin, nbuf, ibin_start, ibin_av &
        , nbra, nrestr, ib, kk, iskip, ngenr, idim, iel &
        , ndim, ddim, nind, rflag, skip_read, nelsan &
        , ngind, ng6ind, nskind, ind, iw, ngen_local &
        , off, qpwftable, rhotable, spintable, indspin0, nrhoind, kkt_r, shift, maxf_r, k_index &
        , j_shift, size_per_proc, nw_per_file, write_start, ntot &
        , pairtable, spairtable, npairind, ii, i_init, i_length, ncorrsamp, nind_corrfun &
        , ioptread, ngendone, lbinr, idigit(6), jj, indcheck, ll, jcol, irow, indkmom, multcell, m
    ! by E. Coccia (24/6/11)
    integer :: itmp
    real*8 history_change_send, history_change_rec &
        , history_change_onsite, history_total &
        , history_change_send_s, history_change_rec_s &
        , history_change_onsite_s, history_total_s, xold, yold &
        , shiftlog, logsamp, logav, countlog, corr_norm, cellscale_sav(3), maxwsk
    ! ndim dimension
    double precision, external :: dlamch
    parameter(ndim=3)
    ! vectors of coordinates
    real(4), dimension(:, :, :), allocatable :: iconfm, iconfmm
    real(4), dimension(:, :, :), allocatable :: walk
    real(8), dimension(:, :), allocatable :: dwalk
    ! pointers for forward walking (population history)
    integer, dimension(:, :), allocatable :: jbra, jbram, jbra_temp
    integer, dimension(:), allocatable :: jbra_old, jbra_step
    ! weights of the walkers
    real(4), dimension(:, :), allocatable :: wconf, wconfm
    real(4), dimension(:), allocatable :: swconf
    real(8), dimension(:), allocatable :: seconfd, spsiln
    ! average weight of generation
    real(8), dimension(:), allocatable :: wbuf, wbufm
    ! vectors of correlation functions
    real(4), dimension(:, :), allocatable :: amisc
    real(4), dimension(:), allocatable :: opar
    real*8 weight
    ! vectors of averaged corr functions per time step propagation
    ! by E. Coccia (22/6/11): For each process!
    ! ek,wk --> average over the bin
    ! ebin,ebin2,wbin --> global average (with binning technique)
    real(8), dimension(:), allocatable :: wk, wbin
    real(8), dimension(:, :), allocatable :: ek, ebin, ebin2
    ! by E. Coccia (22/6/11): Total variables
    real(8), dimension(:), allocatable :: wk_tot
    real(8), dimension(:, :), allocatable :: ek_tot
    ! other vectors and variables
    ! rion(3,nion)
    !real(8), DIMENSION(:,:), ALLOCATABLE :: rion
    real*8 wsk, r_offset(3), r_offset_pair(3), rx(3), ry(3), nooffset(3), rmax, drmax
    real*4 dx(3)
    integer ncell(3), nshlls, max_shells, lptable, nvects, ngrid_l(3), ngrid_p &
        , ind_offset(3), ind_offset_pair(3), vdim(3), iix(3)

    integer, dimension(:), allocatable :: ipart

    real*8 ell(3), vell(3), elli(3), maxg6r, tpiell(3), el2(3), cutr, cutr2 &
        , dtable, dtablei, rho, dxil(3), dxil_p(3), signk, sphere_radius, vecscra(3)

    real(8), dimension(:), allocatable :: density, sdensity, corrsamp
    real(8), dimension(:), allocatable :: paircorr, spaircorr
    real(8), dimension(:), allocatable :: datagrid
    real(8), dimension(:), allocatable :: psip_for
    !real(8), DIMENSION(:), ALLOCATABLE :: atom_number,zetar

    real(8), dimension(:, :), allocatable :: rkcomp, rkcomp_crys
    real(8), dimension(:), allocatable :: rknorm, wtk
    integer, dimension(:), allocatable :: kmult
    real(8), dimension(:), allocatable :: rhok, rhoks, pwmat
    real(8), dimension(:, :), allocatable :: sofk

    integer, dimension(:), allocatable :: simap, simap_tmp
    real(8), dimension(:), allocatable :: kvec, imap
    real(4), dimension(:, :), allocatable :: sangle
    real(8), dimension(:, :), allocatable :: wcorr, ecorr

    real(8), dimension(:), allocatable :: ioniond
    integer, dimension(:), allocatable :: iioniond
    ! variables needed to distinguish atomic species sectors
    integer :: sp1, sp2, nspec, nspec2, tmp_adr, nshellsp
    integer, dimension(:), allocatable :: indsp, adr_spc, sec_spc
    real(8), dimension(:), allocatable :: atom_spec
    logical, dimension(:), allocatable :: check_spec
    type(table_spc), dimension(:), allocatable :: spc_table
    !
    real(8) :: dist
    integer :: ishell
    integer :: nel_read, nelup_read
    integer read_start, size_run, sizep, nstart, nend, nstart_shift, nend_shift &
        , ist, ien, id1, iw_r, idummy

    character(2) chars
    character(3) checkpbc
    character(20) answer
    character(lchlen) :: path, scratchpath, charascratch
    integer ianswer

    real(8) zero(3), rion_center(3), cutk, kspin(3)
    real(4), dimension(:), allocatable :: signpsi

    logical ifrho, fermi_flag, ifpair, ifsofk, ifspin, ifkspin, err_stop, ifcorrs, iffluct
    logical ife, err_read, oldscra, ifsan, longio, noeloc, if_compute, allshells, check_inside

    real(8), parameter :: eps = 1.d-6

    real(8) :: outofplane

    ! (Kosuke Nakano) an logical variable for computing Radial distribution function of atom
    logical rdf_for_atom

    ! optimal allocation:
    ! nbuf = maxk
    ! nm = #_corr_functions
    ! nw = nwr

#ifdef PARALLEL
    include 'mpif.h'
    integer n1, ierr, nrank, srank, status(MPI_STATUS_SIZE), ithread
    integer(kind=MPI_OFFSET_KIND) disp
    real(4), allocatable :: SP_buffer(:, :)
    real(8), allocatable :: DP_buffer(:, :)
    integer :: buffer_counter, buffer_depth, SP_block_size, DP_block_size
#else
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'readforward'
        call help_online(name_tool)
        stop
    end if
#endif

#ifdef PARALLEL
    call mpi_init(ierr)
!call mpi_init_thread(MPI_THREAD_FUNNELED,ithread,ierr)
    call mpi_comm_size(MPI_COMM_WORLD, sizep, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (rank .eq. 0) then
        write (6, *) 'PARALLEL CALCULATION'
        write (6, *) 'largest tag value for point-to-point communication', MPI_TAG_UB
!  if(rank.eq.0) write(6,*) ' Initial mpi value of threads',ithread
    end if
! define also these values for the read_fort10
    rankrep = rank
    commrep_mpi = MPI_COMM_WORLD
    rankcolrep = 0
    commcolrep_mpi = MPI_COMM_WORLD
#else
    rank = 0
    rankrep = 0
    rankcolrep = 0
    commrep_mpi = 0
    commcolrep_mpi = 0
    sizep = 1
#endif

    ! output version information
    if (rank .eq. 0) call print_version

    zero = 0.d0
    max_shells = 20000 ! max number of k-shells allowed

    ! Obviously we are not doing minimization
    yesmin = 0

    maxwsk = DLAMCH('O')/100.d0 ! Maximum threshold for overflow

    if (rank .eq. 0) call read_fort11_begin

    call init_variables(nel, nelup, nion, iespbc, celldm, rs)

    if (rank .eq. 0) then
        call read_corr_fun(nel, nelup, nion, iespbc, celldm, rs, cellscale &
                           , ngenr, ell, ncell, nbias, maxf, ibinit, lbin, iskip, ddim, cutk, vdim, ngrid_l, r_offset, ngrid_p &
                           , ifrho, ifspin, ifkspin, kspin, ifpair, iffluct, ifsofk, fermi_flag, err_stop, ifcorrs, shiftlog &
                           , ioptread, longio, noeloc, sphere_radius, allshells, outofplane, rdf_for_atom)
        size_run = nprocr
        nelup_read = nelup
        nel_read = nel
    end if

    allocate (winvbarfn(1), winvfn(1), vpotsav_ee(1, 1))

    add_onebody2det = .false.

#ifdef PARALLEL
    call mpi_bcast(nelup_read, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nel_read, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nwr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ioptread, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(size_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ncell, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ell, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(kspin, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(cellscale, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbias, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(maxf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(lbin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ibinit, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iskip, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(lptable, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifrho, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifcorrs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(noeloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifspin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifspin2, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifrho_corr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(decouple_k, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(decouple_files, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
! QPWF BCASTING variables ----------------------------------------------
    call mpi_bcast(ifqpwf, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifqpwf_h, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifqpwf_k, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifqpwf_extr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifberry, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_extr_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(extr_point, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!-----------------------------------------------------------------------
    call mpi_bcast(ifpair, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iffluct, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(fermi_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ddim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ngrid_l, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ngrid_p, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vdim, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(r_offset, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(shiftlog, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifsofk, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(cutk, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(err_stop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifkspin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(longio, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifrho_assar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nswitch, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(kswitch, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(assar_cut, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(assar_parr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ext_grid, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(grid_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(center, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(da, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(grid_start, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(sphere_radius, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(allshells, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifdipole, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(outofplane, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rdf_for_atom, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nel, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nelup, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nion, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iespbc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(molyes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    multcell = ncell(1)*ncell(2)*ncell(3)

    cellscale_sav(1:3) = cellscale(1:3)*ncell(1:3)

!test
!    write(6,*) 'cellscale test',cellscale_sav(1:3), cellscale(1:3), ncell(1:3)

    if (.not. decouple_k .and. decouple_files) then
        if (rank .eq. 0) write (6, *) &
            ' Warning decouple_files only with decouple_k = true!'
        decouple_files = .false.
    end if

#ifdef PARALLEL
    oldscra = .false.
#else
    oldscra = .true.
#endif

    if (err_stop) then
        if (rank .eq. 0) write (6, *) 'I am sorry: something went wrong in reading readforward.input file'
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif
        stop
    end if

    ! now read qmc simulation input
    ! Initialize flags error
    yesdft = .false.
    iflagerrall = 0
    iflagerr = 0
    nproc = sizep

    call get_dir(path)
    if (rank .eq. 0) write (*, *) ' Initial path : ', path
    epsmach = 1000.d0*dlamch('e')
    safemin = 10000.d0*dlamch('s')

    call read_datasmin

    ! Restore possible change of cellscale by read_datasmin
    cellscale(1:3) = cellscale_sav(1:3)/ncell(1:3)

    if (nk .le. 1) decouple_k = .false.
    if (.not. decouple_k .and. decouple_files) decouple_files = .false.

    if (nk .gt. 1 .and. decouple_k) then

        if (rank .eq. 0) then
            write (6, *) ' weight given '
            do i = 1, nk
                write (6, *) ' K / weight =', i, wkp(i)
            end do
        end if
        mcol = nproc/nk
        if (mcol*nk .ne. nproc .or. nproc .ne. size_run) then
            if (rank .eq. 0) write (6, *) ' # processors / # k-points :', nproc, nk
            if (mcol*nk .ne. nproc) then
                call error(' Initializeall ', ' # of processors must be multiple of &
                        &  the # of k-points !!', 1, rank)
            else
                call error(' Initializeall ', ' # of processors must be the same  &
                        &  of the vmc !!', 1, rank)
            end if
        end if
#ifdef PARALLEL
        ! split the communicators
        irow = rank/mcol
        jcol = mod(rank, mcol)
        call mpi_comm_split(mpi_comm_world, irow, jcol, commrep_mpi, ierr)
        call mpi_comm_split(mpi_comm_world, jcol, irow, commcolrep_mpi, ierr)
        call mpi_comm_rank(commrep_mpi, rankrep, ierr)
        call mpi_comm_rank(commcolrep_mpi, rankcolrep, ierr)
        call mpi_barrier(MPI_COMM_WORLD, ierr) ! for unreliable networks
#endif

    else

        nk = 1

    end if

    if (decouple_k) then
        indkmom = (rank + nproc/nk)/(nproc/nk)
        call convertdec(indkmom, idigit)

        if (indkmom .gt. 99999) then
            chara = char(idigit(1))//char(idigit(2))//char(idigit(3))// &
                    &char(idigit(4))//char(idigit(5))//char(idigit(6))
        elseif (indkmom .gt. 9999) then
            chara = char(idigit(2))//char(idigit(3))// &
                    &char(idigit(4))//char(idigit(5))//char(idigit(6))
        elseif (indkmom .gt. 999) then
            chara = char(idigit(3))//char(idigit(4))//char(idigit(5))//char(idigit(6))
        elseif (indkmom .gt. 99) then
            chara = char(idigit(4))//char(idigit(5))//char(idigit(6))
        elseif (indkmom .gt. 9) then
            chara = char(idigit(5))//char(idigit(6))
        else
            chara = char(idigit(6))
        end if

    end if

    if (longio .or. ifcorrs) then
        if (decouple_k) then
            open (unit=14, file='binsK'//trim(chara), form='formatted', status='unknown')
        else
            open (unit=14, file='bins.dat', form='formatted', status='unknown')
        end if
    end if

    ! In ANY event do not use regularization.
    epscut = 0.d0
    epscuttype = 0
    if (ifrho_assar) then
        nintpseudo = 0
        nintpsa = 0
    end if

    iopt = ioptread

    if (oldscra) then
        scratchpath = './'
        wherescratch = './scratch'
    else
#ifdef __KCOMP
        scratchpath = "./"
#else
        scratchpath = trim(path)//"/turborvb.scratch/"
#endif
        wherescratch = trim(scratchpath)//"tmp"
    end if

    if_compute = .false.

    if (ifpair) if_compute = .true.

    !if(ifcorrs.or.ifrho_assar.or.ifspin2.or.ifqpwf.or.ifrho_corr.or.ifdipole) then
    if (ifcorrs .or. ifrho_assar .or. ifspin2 .or. ifqpwf .or. ifrho_corr) then
        ! MC: when I compute dipole moments I do not need the wave function!
        ! need of pseudo and full wave function

        call read_pseudo
        if_compute = .true.

        if (rankrep .eq. 0) then
            !   LEGGI IL FORT.10 APPROPRIATO
#ifdef __KCOMP
            call convertdec(rank, idigit)
#else
            call convertdec(rankcolrep, idigit)
#endif
            charaq = char(idigit(1))//char(idigit(2))//char(idigit(3))// &
                     char(idigit(4))//char(idigit(5))//char(idigit(6))
            !      if(ifrho_assar.or.ifspin2.or.ifdipole) then
            if (ifrho_assar .or. ifspin2 .or. ifspin) then
                if (decouple_k .and. molyes) then
                    if (rank .eq. 0) write (6, *) ' Warning reading correlated wf from turborvb.scratch folder'
                    open (unit=10, file=trim(scratchpath)//'fort.10_'//trim(charaq), form='formatted', status='unknown')
                else
                    if (rank .eq. 0) write (6, *) ' Warning reading correlated wf from fort.10_corr '
                    open (unit=10, file='fort.10', form='formatted', status='unknown')
                end if
            else
                if (decouple_k .and. molyes) then
                    if (rank .eq. 0) write (6, *) ' Warning reading correlated wf from turborvb.scratch folder'
                    open (unit=10, file=trim(scratchpath)//'fort.10_'//trim(charaq), form='formatted', status='unknown')
                else
                    if (rank .eq. 0) write (6, *) ' Warning reading correlated wf from fort.10_corr '
                    open (unit=10, file=trim(scratchpath)//'fort.10_'//trim(charaq), form='formatted', status='unknown')
                    open (unit=10, file='fort.10_corr', form='formatted', status='unknown')
                end if
            end if
        end if
        nwm = nw
        nws = in1

        if (ifrho_assar) npsa = 0 !(Matteo) Excluding pseudo

        call read_fort10(10)
        off = 1
        do l = 1, ddim
            ind_offset(l) = off
            off = off*ngrid_l(l)
            dxil(l) = ngrid_l(l)/ell(vdim(l))
        end do
        if (.not. iespbc) then
            rion_center(:) = rion(:, 1)
            do i = 2, nion
                rion_center(:) = rion_center(:) + rion(:, i)
            end do
            rion_center(:) = rion_center(:)/nion
            if (rank .eq. 0) then
                write (6, *) 'geometrical center of the molecule', rion_center(:)
                write (6, *) 'offset shifted by geometrical center of the molecule'
            end if
            r_offset(:) = r_offset(:) + rion_center(:)
            r_offset_pair(1:3) = -cellscale(1:3)/2.d0 - 0.5d0/dxil(:)
        else
            r_offset_pair(:) = -0.5d0/dxil(:)
        end if

        ! (Matteo) Building the grid for Assaraf or STM calculations
        if (ifrho_assar .or. ifqpwf .or. ifrho_corr) then
            call compute_grid(ngrid_l, ell, vell)
        else
            allocate (out_grid(3, 1))
            out_grid = 0.d0
        end if

        call update_nmolfn

        if (rank .eq. 0) write (6, *) ' firstmol nmolfn =', firstmol, nmolfn

        if (LBox .gt. 0) then
            call InitEwald(nion, zetar, nel, nws)
            if (rank .eq. 0) then
                write (*, *) ' Number of G vectors ', n_gvec
                write (*, *) ' Ewald Self Energy ', eself
            end if
            kmax = n_gvec
            kmax2 = 2*kmax
        else
            n_gvec = 1
            kmax = 0
            kmax2 = 0
            if (iesbra) then
                allocate (sum_q_cos_gr(1, 1), sum_q_sin_gr(1, 1))
                sum_q_cos_gr = 0.d0
                sum_q_sin_gr = 0.d0
            end if
        end if

    else

        allocate (zetar(nion), atom_number(nion), rion(3, nion))

        if (rank .eq. 0) then
            open (unit=10, file='fort.10', form='formatted', status='unknown')
            call read_fort10_fast
            !     write(6,*) ' after read fort10 fast 1'
            zetar = zetar_fast
            rion = rion_fast
            atom_number = atom_number_fast
            close (10)
        end if

#ifdef PARALLEL
        call bcast_real(zetar, nion, 0, MPI_COMM_WORLD)
        call bcast_real(celldm, 6, 0, MPI_COMM_WORLD)
        call bcast_real(s2r, 9, 0, MPI_COMM_WORLD)
        call bcast_real(rs, 1, 0, MPI_COMM_WORLD)
        call bcast_real(atom_number, nion, 0, MPI_COMM_WORLD)
        call bcast_real(rion, 3*nion, 0, MPI_COMM_WORLD)
        call mpi_bcast(ipc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(yes_tilted, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    end if
#ifdef _CUSOLVER
#ifdef RISC
    call cusolver_handle_init_(handle)
#else
    call cusolver_handle_init(handle)
#endif
#endif
#ifdef _OFFLOAD
    yes_ontarget = .false. ! everything should be back to cpu no matter how
!$omp target data map(to: jasmat,muj_c,jasmat_c,detmat,detmat_c&
!$omp &,projm,mu_c,eagp_pfaff)&
!$omp &  map(to:psip&
!$omp &,winv,winvj&
!$omp &,agp,agpn,ainv,winvbar,winvjbar,winvfn,winvbarfn,ainvup,ainvdo)
#ifdef _CUSOLVER
    if (ipf .ne. 2) then
        if (rank .eq. 0) write (6, *) ' Warning, using cusolver routines '
        ldworkspace = 1
        lzworkspace = 1
        !
        if (ipc .eq. 1) then
#ifdef RISC
            call cusolver_dgetrf_buffersize_(handle, stat, nelup, nelup, psip, nelup, ldworkspace)
#else
            call cusolver_dgetrf_buffersize(handle, stat, nelup, nelup, psip, nelup, ldworkspace)
#endif
            allocate (dev_dgetri_workspace(nelup, nelup))
            allocate (dev_zgetri_workspace(1, 1))
        else
#ifdef RISC
            call cusolver_zgetrf_buffersize_(handle, stat, nelup, nelup, psip, nelup, lzworkspace)
#else
            call cusolver_zgetrf_buffersize(handle, stat, nelup, nelup, psip, nelup, lzworkspace)
#endif
            allocate (dev_dgetri_workspace(1, 1))
            allocate (dev_zgetri_workspace(nelup, nelup))
        end if
        call checkiflagerr(stat, rank, ' Something went wrong in calculating GPU buffer space')
        allocate (dev_dgetrf_workspace(ldworkspace))
        allocate (dev_zgetrf_workspace(lzworkspace))
        dev_Info = 0
        dev_dgetrf_workspace = 0
        dev_zgetrf_workspace = 0
        dev_dgetri_workspace = 0
        dev_zgetri_workspace = 0
    end if
    !$omp target data map(to:ipsip, dev_Info&
    !$omp &,dev_dgetri_workspace,dev_zgetri_workspace&
    !$omp &,dev_dgetrf_workspace,dev_zgetrf_workspace) if(ipf.ne.2)
    if (rank .eq. 0) write (6, *) ' GPU memory for cusolver allocated'

#endif
#else
    yes_ontarget = .false.
#endif

    ! CCC
    if (ifspin2) then
        allocate (ratiospin(nelup, neldo))
        !call read_fort10(10)
        !       if (rank.eq.0) then
        if (ipf .eq. 2) then
            if (contraction .eq. 0) call inits2pfaff(detmat)
            if (contraction .ne. 0) call inits2pfaff_c(detmat_c, mu_c)
        end if
        !       end if
        !       call mpi_finalize(ierr)
        !       stop
    end if

    if (ifberry) then
        if (ipc .eq. 2) then
            yes_complex = .true.
        else
            yes_complex = .false.
        end if

        if (rs .lt. 0.d0) then
            rs = (-3.d0/4.d0/Pi*rs**3*celldm(2)*celldm(3)/nel)**(1.d0/3.d0)
        end if
        celldm(1) = 1.d0
        !      celldm(4:6) = 90d0*PI/180.d0
        omega = celldm(2)*celldm(3)
        celldm(1) = (PI*nel*4.d0/3.d0/omega)**(1.d0/3.d0)*rs

        if (yes_tilted) then
            givens2r = .true.
        else
            givens2r = .false.
        end if
        call InitCell(nion, nel, yes_complex)
        call init_berry_phase(nel)
        if (rank .eq. 0) then
            write (*, "(a13,2x,i15)") '# [nel]      ', nel
            write (*, "(a13,3(2x,f15.8))") '# [rec. cell]', berry_comp_vec(1), berry_comp_vec(2), berry_comp_vec(3)
        end if
        !          write( 120,"(a13,2x,i15)" )     '# [nel]      ', nel
        !          write( 120,"(a13,3(2x,f15.8))" ) '# [rec. cell]',berry_comp_vec(1),berry_comp_vec(2),berry_comp_vec(3)
    end if

#ifdef PARALLEL
    if (mod(size_run, sizep) .ne. 0) then
        if (rank .eq. 0) then
            write (6, *) '# files must be the same for every process !!'
            write (6, *) '# processors', sizep
            write (6, *) '# files', size_run
        end if
        call mpi_finalize(ierr)
        stop
    else
        size_per_proc = size_run/sizep
        if (rank .eq. 0) write (6, *) '# files per process', size_per_proc
    end if
#else
    size_per_proc = size_run
#endif

    read_start = 30
    write_start = read_start + size_per_proc

    !!!!!STARTING ALLOCATION CORR FUN VECTORS

    allocate (datagrid(1))

    if (ifkspin) then
        kspin(1:3) = 2.d0*pi*kspin(1:3)/cellscale(1:3)
        if (rank .eq. 0) write (6, *) ' Chosen momentum structure factor =', kspin(1:3)
    end if

    !by E. Coccia (22/11/11): read the external electric field
    if (ext_pot) then
        call extpot_read()
    end if
    !by E. Coccia (22/11/11): read the vdW parameters
    if (vdw) then
        call vdw_read()
    end if

    ! RHO Assaraf: G(r),am,bm .......................................................
    ! Calling initialization for Assaraf module by O. Chernomor (Modified by M. Barborini)

    if (ifrho_assar) then
        call initialize_assaraf(vell, ngrid_l)
    end if
    !END rho Assaraf coef's.................................................................................

    nrhoind = 0
    npairind = 0
    if (ifrho .or. ifspin) then

        if (sphere_radius .eq. 0.d0) then
            if (.not. ifspin2) then
                off = 1
                do l = 1, ddim
                    ind_offset(l) = off
                    off = off*ngrid_l(l)
                    dxil(l) = ngrid_l(l)/ell(vdim(l))
                end do
            end if
            nrhoind = off

        else

#ifdef PARALLEL
            call mpi_bcast(nel, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nelup, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nion, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(iespbc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

            if (.not. allocated(zetar)) allocate (zetar(nion))
            if (.not. allocated(atom_number)) allocate (atom_number(nion))
            if (.not. allocated(rion)) allocate (rion(3, nion))
            if (.not. allocated(atom_spec)) allocate (atom_spec(nion))
            atom_spec = 0.d0
            nspec = 0

            if (rank .eq. 0) then
                open (unit=10, file='fort.10', form='formatted', status='unknown')
                !If ifspin2 it has already be done
                if (.not. ifspin2) call read_fort10_fast
                zetar = zetar_fast
                rion = rion_fast
                atom_number = atom_number_fast
                close (10)
                ! find and save atomic species
                allocate (check_spec(nion))
                check_spec = .false.
                do i = 1, nion
                    if (.not. check_spec(i)) then
                        nspec = nspec + 1
                        atom_spec(nspec) = atom_number(i)
                    end if
                    do j = 1, nion
                        if (atom_number(i) .eq. atom_number(j) .and. .not. check_spec(j)) &
                            check_spec(j) = .true.
                    end do
                end do
                nspec2 = nspec*nspec
                deallocate (check_spec)
            end if
            !
            if (rank .eq. 0) then
                write (6, *) ' Atomic species found:'
                do i = 1, nspec
                    write (6, *) i, atom_spec(i)
                end do
            end if

#ifdef PARALLEL
            call bcast_real(zetar, nion, 0, MPI_COMM_WORLD)
            call bcast_real(atom_number, nion, 0, MPI_COMM_WORLD)
            call bcast_real(rion, 3*nion, 0, MPI_COMM_WORLD)
            call bcast_real(atom_spec, nion, 0, MPI_COMM_WORLD)
            call bcast_real(nspec, 1, 0, MPI_COMM_WORLD)
            call bcast_real(nspec2, 1, 0, MPI_COMM_WORLD)
#endif

            nrhoind = nion
            if (allocated(dists_kel)) deallocate (dists_kel)
            allocate (dists_kel(nion, nel))

        end if

    end if

    if (.not. if_compute .and. allocated(rion)) then
        if (.not. allocated(rionsav)) allocate (rionsav(3, nion))
        ! shift ions by -0.5d0/dxil(:) for plotting with xcrysden
        ! (bin center is shifted by 0.5d0/dxil(:) with respect to bin index)
        if (.not. yes_tilted) then
        do i = 1, nion
            rionsav(:, i) = rion(:, i) - r_offset(:) - 0.5d0/dxil(:)
        end do
        else
        vecscra(:) = 0.d0
        do j = 1, 3
            vecscra(:) = vecscra(:) - s2r(:, j)*0.5d0/dxil(j)/cellscale(j)
        end do
        vecscra(:) = vecscra(:) - r_offset(:)
        do i = 1, nion
            rionsav(:, i) = rion(:, i) + vecscra(:)
        end do
        end if
        if (iespbc) call DMyApplyPBC(rionsav, nion, cellscale_sav) ! Put the ion in (0,L)
    end if

    if (ifpair) then

        if (sphere_radius .eq. 0.d0) then

            do l = 1, ddim
                ind_offset_pair(l) = ind_offset(l)
                dxil_p(l) = dxil(l)
            end do
            npairind = nrhoind

            rmax = ell(vdim(1))/2.d0
            do l = 2, ddim
                rmax = max(rmax, ell(vdim(l))/2.d0)
            end do
            drmax = dble(ngrid_p)/rmax

            if (rank .eq. 0) write (6, *) 'rmax and radial step found in g(r)', rmax, 1.d0/drmax

        else

            ! build a ion-ion mapping to distance and
            ! compute ion-ion distances
            allocate (ioniond(nion*nion), iioniond(nion*nion), mult(2*nion*nion), &
                      ipip_adr(2*nion*nion), dists(2*nion*nion))
            allocate (indsp(nion*nion))
            !
            mult = 0
            ipip_adr = 0
            nshell = 0
            !
            call el_ion_distance(ioniond, rion, nion, rion, nion, iespbc)

            call dsortx(ioniond, 1, nion*nion, iioniond)

            ! shell degeneracy
            dist = ioniond(1)
            dists(1) = dist
            ipip_adr(iioniond(1)) = 1
            nshell = 1
            mult(1) = 1
            ishell = 1
            i = 1
            do i = 2, nion*nion
                if (ioniond(i) .gt. dist + eps .or. allshells) then
                    ishell = ishell + 1
                    dist = ioniond(i)
                    dists(ishell) = dist
                    nshell = nshell + 1
                end if
                mult(ishell) = mult(ishell) + 1
                ipip_adr(iioniond(i)) = ishell
            end do
            ! shell degenerancy for species sector
            !
            ! find species sector for each atomic distance
            !
            do i = 1, nion
                do j = 1, nion
                    ind = (i - 1)*nion + j
                    do k = 1, nspec
                        if (atom_number(j) .eq. atom_spec(k)) sp1 = k
                        if (atom_number(i) .eq. atom_spec(k)) sp2 = k
                    end do
                    indsp(ind) = (sp1 - 1)*nspec + sp2
                end do
            end do

            !
            ! fill atomic species table to build species-species map
            !
            allocate (spc_table(nspec2), adr_spc(nion*nion + 1), sec_spc(nion*nion + 1))

            do i = 1, nspec2
                allocate (spc_table(i)%dists(nion*nion + 1), spc_table(i)%indd(nion*nion + 1))
            end do

            ind = 0
            spc_table(:)%nsp = 0
            do i = 1, nspec2
                tmp_adr = 0
                do j = 1, nion
                    do k = 1, nion
                        ind = k + (j - 1)*nion
                        dist = ioniond(ind)
                        if (indsp(iioniond(ind)) .eq. i) then
                            tmp_adr = tmp_adr + 1
                            spc_table(i)%dists(tmp_adr) = dist
                            spc_table(i)%indd(tmp_adr) = iioniond(ind)
                            spc_table(i)%nsp = spc_table(i)%nsp + 1
                        end if
                    end do
                end do
            end do

            !
            ! build the mapping given the species table. Use the same variables
            ! used for the general pair corr functions with suitable shift
            !
            ishell = nshell + 1 ! right after general pair corr
            adr_spc = 0
            sec_spc = 0
            do i = 1, nspec2
                dist = spc_table(i)%dists(1)
                dists(ishell) = dist
                sec_spc(ishell) = i
                do j = 1, spc_table(i)%nsp
                    if (spc_table(i)%dists(j) .gt. dist + eps) then
                        ishell = ishell + 1
                        dist = spc_table(i)%dists(j)
                        dists(ishell) = dist
                        sec_spc(ishell) = i
                    end if
                    mult(ishell) = mult(ishell) + 1
                    ! address for comp_corr_fun
                    tmp_adr = spc_table(i)%indd(j)
                    adr_spc(tmp_adr) = ishell
                end do
                ishell = ishell + 1
            end do

            nshellsp = ishell - nshell - 1
            if (rank .eq. 0) then
                write (6, *) 'number of non-inequivalent shells in site-site corr fun', nshell + nshellsp
                write (6, *) 'total number of inequivalent shells for atomic species pairs', nshellsp
                write (6, *) 'total number of site-site pairs', sum(mult(1:nshell))
            end if
            npairind = nshell + nshellsp
        end if

    end if

    if (.not. allocated(adr_spc)) allocate (adr_spc(1))
    if (.not. allocated(sec_spc)) allocate (sec_spc(1))
    if (.not. allocated(iioniond)) allocate (iioniond(1))

    if (ifsofk) then
        ! structure factor
        tpiell = 0.d0
        do l = 1, ndim
            tpiell(l) = (2.d0*Pi)/cellscale(vdim(l))
        end do

        allocate (rkcomp(ndim, max_shells), kmult(0:max_shells), rkcomp_crys(ndim, max_shells))
        allocate (rknorm(0:max_shells), wtk(0:max_shells))
        rkcomp = 0.d0
        rkcomp_crys = 0.d0

        ! build k-grid based on cutk in crystal space
        call shells(ddim, tpiell, cutk, nshlls, rkcomp, rknorm, kmult, nvects, max_shells, max_shells, ndim, rank)

        nvects = kmult(nshlls - 1)
        if (rank .eq. 0) write (6, *) '# k-points with |k| < k_cut ', nvects

        ! sorting of the kvectors
        ! determine simap
        allocate (kvec(2*nvects), simap(2*nvects), imap(2*nvects), simap_tmp(2*nvects))

        do i = 1, nvects
            kvec(i) = rkcomp(1, i)
            imap(i) = dble(i)
            simap(i) = i
            kvec(i + nvects) = -rkcomp(1, i)
            imap(i + nvects) = dble(i + nvects)
            simap(i + nvects) = -i
        end do
        simap_tmp = simap

        ! sorting x
        call dsort(kvec, imap, 2*nvects, 2)
        do i = 1, 2*nvects
            simap(i) = simap_tmp(int(imap(i)))
        end do
        simap_tmp = simap

        if (ddim .gt. 1) then

            do i = 1, 2*nvects
                imap(i) = dble(i)
                signk = abs(simap(i))/simap(i)
                kvec(i) = signk*rkcomp(2, abs(simap(i)))
            end do

            ! sorting y
            do i = 1, 2*nvects

                if (i .eq. 1) then

                    i_init = 1
                    i_length = 1
                    signk = abs(simap(i))/simap(i)
                    xold = signk*rkcomp(1, abs(simap(i)))

                elseif (abs(simap(i))/simap(i)*rkcomp(1, abs(simap(i))) .ne. xold .or. i .eq. 2*nvects) then

                    if (i .eq. 2*nvects) then
                        i_length = i_length + 1
                    end if

                    call dsort(kvec(i_init), imap(i_init), i_length, 2)
                    do ii = i_init, i_init + i_length - 1
                        simap_tmp(ii) = simap(int(imap(ii)))
                    end do

                    i_init = i
                    i_length = 1
                    signk = abs(simap(i))/simap(i)
                    xold = signk*rkcomp(1, abs(simap(i)))

                else

                    i_length = i_length + 1
                    signk = abs(simap(i))/simap(i)
                    xold = signk*rkcomp(1, abs(simap(i)))

                end if

            end do
            simap = simap_tmp

            if (ddim .eq. 3) then

                do i = 1, 2*nvects
                    imap(i) = dble(i)
                    signk = abs(simap(i))/simap(i)
                    kvec(i) = signk*rkcomp(3, abs(simap(i)))
                end do

                ! sorting z
                do i = 1, 2*nvects

                    if (i .eq. 1) then

                        i_init = 1
                        i_length = 1
                        signk = abs(simap(i))/simap(i)
                        xold = signk*rkcomp(1, abs(simap(i)))
                        yold = signk*rkcomp(2, abs(simap(i)))

                    elseif (abs(simap(i))/simap(i)*rkcomp(1, abs(simap(i))) .ne. xold .or. &
                            abs(simap(i))/simap(i)*rkcomp(2, abs(simap(i))) .ne. yold .or. &
                            i .eq. 2*nvects) then

                        if (i .eq. 2*nvects) then
                            i_length = i_length + 1
                        end if

                        call dsort(kvec(i_init), imap(i_init), i_length, 2)
                        do ii = i_init, i_init + i_length - 1
                            simap_tmp(ii) = simap(int(imap(ii)))
                        end do

                        i_init = i
                        i_length = 1
                        signk = abs(simap(i))/simap(i)
                        xold = signk*rkcomp(1, abs(simap(i)))
                        yold = signk*rkcomp(1, abs(simap(i)))

                    else

                        i_length = i_length + 1
                        signk = abs(simap(i))/simap(i)
                        xold = signk*rkcomp(1, abs(simap(i)))
                        yold = signk*rkcomp(2, abs(simap(i)))

                    end if

                end do
                simap = simap_tmp

            end if

        end if

        rkcomp_crys = rkcomp
        if (yes_tilted) then
            do i = 1, nvects
                vecscra(:) = 0.d0
                do j = 1, 3
                    vecscra(:) = vecscra(:) + recip(j, :)*rkcomp(j, i)/tpiell(j)
                end do
                rkcomp(:, i) = vecscra(:)
            end do

            !test scalar
!           if(rank.eq.0) then
!           write(6,*) 'k_i in crystal index'
!           do i=1,nvects
!              !              write(6,*) i,rkcomp(:,i) * cellscale(vdim(:)) / 2.0 / Pi
!              write(6,*) i,rkcomp_crys(:,i) / tpiell(:)
!           enddo
!           write(6,*) 'orthogonality TEST: reciprocal_lattice_vectors_i dot lattice_vector_j / 2 / Pi'
!           do i=1,3
!              do j=1,3
!                 write(6,*) i,j,sum(s2r(:,j)*recip(i,:)) / 2.0 / Pi
!              enddo
!           enddo
!           write(6,*) 'orthogonality TEST: k_i dot lattice_vector_j / 2 / Pi'
!           do i=1,nvects
!              do j=1,3
!                 write(6,*) i,j,sum(s2r(:,j)*rkcomp(:,i)) / 2.0 / Pi
!              enddo
!           enddo
!           endif

        end if

        deallocate (simap_tmp, kvec, imap)

    else

        allocate (rkcomp(ndim, 1), kmult(0:1), rkcomp_crys(ndim, 1))
        allocate (rknorm(0:1), wtk(0:1))
        allocate (simap(1))

    end if

    if (ifdipole) then

        allocate (nuclear_dipole(ddim))

        if (rank .eq. 0) write (6, *) 'offset applied', r_offset

        nuclear_dipole = 0.d0
        do i = 1, nion

            nuclear_dipole(1:ddim) = nuclear_dipole(1:ddim) + zetar(i)*rion(1:ddim, i)
            !         write(6,*) i,el2,zetar(i)

        end do

    end if

    countlog = 0.d0
    logav = 0.d0

    if (rank .eq. 0) write (6, *) '***************************************************'
    qpwftable = 0
    rhotable = 0
    spintable = 0
    pairtable = 0
    spairtable = 0
    nskind = 0
    nind = 0
    ncorrsamp = 0

    !(Ye) Adding Spin^2 calculation
    if (ifspin2) then
        if (rank .eq. 0) write (6, *) "COMPUTING Spin^2"
        nind = nind + 1
        nspin2 = 1
    else
        nspin2 = 0
    end if

    if (ifberry) then
        nind = nind + 6
    end if

    !(Matteo) Adding Rho calculation of the electron density with Assaraf Method
    if (ifrho_assar) then
        if (rank .eq. 0) then
            write (6, *) "COMPUTING Rho with Assaraf Method"
        end if
        nrhoind = grid_points
        rhotable = grid_points
        nind = nind + grid_points
        allocate (density(max(grid_points, 1)))
        pairtable = npairind
        nind = nind + npairind
        allocate (paircorr(max(npairind, 1)))
    end if

    ! Adding calculation of dipole moment
    if (ifdipole) then
        if (rank .eq. 0) write (6, *) "COMPUTING Dipole and Quadrupole Moment"
        ndipole = ddim
        nquad = ddim
        allocate (dipole(ndipole), quad(nquad, nquad), cc_chg(nquad), quad_diag(nquad))
        nind = nind + ndipole + nquad**2
    else
        ndipole = 0
        nquad = 0
    end if

    !(Matteo) Rho calculation with new method
    if (ifrho_corr) then
        if (rank .eq. 0) then
            write (6, *) "COMPUTING Rho with new method"
        end if
        nrhoind = grid_points
        rhotable = grid_points
        nind = nind + grid_points
        allocate (density(max(grid_points, 1)))
    end if

    !-----------------------------------------------------------
    if (ifrho .and. .not. allocated(density)) then
        if (rank .eq. 0) then
            write (6, *) 'COMPUTING rho'
            if (ddim .eq. 2) then
                write (6, *) 'density contour plot'
            end if
        end if
        !  nrhoind=nrhoind+3   ! adding order parameters
        rhotable = nrhoind
        nind = nind + nrhoind
        allocate (density(max(nrhoind, 1)))
        if (allocated(datagrid)) deallocate (datagrid)
        allocate (datagrid((ngrid_l(1) + 1)*(ngrid_l(2) + 1)*(ngrid_l(3) + 1)))

        if (ifpair) then
            if (rank .eq. 0) then
                write (6, *) 'COMPUTING charge-charge corr fun'
            end if
            pairtable = npairind + ngrid_p
        else
            pairtable = npairind
        end if

        nind = nind + pairtable
        allocate (paircorr(max(pairtable, 1)))
        allocate (ipart(nel))

    elseif (.not. allocated(density)) then
        allocate (density(1), paircorr(1))
        allocate (ipart(1))

    end if

    if (ifspin) then
        if (rank .eq. 0) then
            write (6, *) 'COMPUTING spin'
            if (ddim .eq. 2) then
                write (6, *) 'spin density contour plot'
            end if
        end if
        if (.not. ifkspin) then
            spintable = nrhoind + 1
        else
            spintable = nrhoind + 2
        end if
        nind = nind + spintable
        allocate (sdensity(spintable))
        if (allocated(datagrid)) deallocate (datagrid)
        allocate (datagrid((ngrid_l(1) + 1)*(ngrid_l(2) + 1)*(ngrid_l(3) + 1)))

        if (ifpair) then
            if (rank .eq. 0) then
                write (6, *) 'COMPUTING spin-spin corr fun'
            end if
            spairtable = npairind + ngrid_p
        else
            spairtable = npairind
        end if

        nind = nind + spairtable
        allocate (spaircorr(max(spairtable, 1)))

    else
        allocate (sdensity(1), spaircorr(1))
    end if

    if (ifsofk) then
        if (rank .eq. 0) write (6, *) 'COMPUTING S(k)'
        nskind = nvects
        if (fermi_flag) then
            nind = nind + 5*nvects
        else
            nind = nind + nvects
        end if
        allocate (rhok(2*nvects), rhoks(4*nvects), pwmat(2*nvects*nel), sofk(nvects, 5))
    else
        nvects = 1
        allocate (rhok(1), rhoks(1), pwmat(1), sofk(1, 1))
    end if

    !-----------------------------------------
    !(Matteo) QPWF calculations
    if (ifqpwf .and. .not. ifqpwf_extr .and. .not. ifqpwf_k) then
        if (rank .eq. 0) write (6, *) 'COMPUTING QPWF'
        nrhoind = grid_points
        qpwftable = ipc*grid_points
        nind = nind + ipc*grid_points
        allocate (qpwf_image(max(ipc*grid_points, 1)))
    elseif (ifqpwf_k) then
        if (rank .eq. 0) write (6, *) 'COMPUTING QPWF IN MOMENTUM SPACE'
        nrhoind = 2*grid_points
        qpwftable = 2*grid_points
        nind = nind + 2*grid_points
        allocate (qpwf_image(max(2*grid_points, 1)))
    elseif (ifqpwf_extr) then
        if (rank .eq. 0) write (6, *) 'COMPUTING QPWF'
        nrhoind = n_extr_points
        qpwftable = ipc*n_extr_points
        nind = nind + ipc*n_extr_points
        allocate (qpwf_image(max(ipc*n_extr_points, 1)))
    else
        allocate (qpwf_image(1))
    end if

    !-----------------------------------------
    if (ifcorrs) then
        if (rank .eq. 0) write (6, *) 'COMPUTING Correlated Sampling'
        ncorrsamp = 3 + ipc
        nind = nind + ncorrsamp
        allocate (corrsamp(ncorrsamp))
        if (rankrep .eq. 0) then
            call rand_init(abs(iseedr))
        end if
    else
        allocate (corrsamp(1))
    end if

    if (nwr .ne. 0) nw = nwr

    maxf_r = maxf/iskip
    maxk = maxf + nbias

    if (ifcorrs) then
        nbuf = lbin
    else
        nbuf = max(maxk, 10) ! store at least a buffer of 10 lines
    end if

    nwbuf = nw*nbuf

    if (ifcorrs) then
        nind_corrfun = nind + 2
    else
        nind_corrfun = nind
    end if

    if (rank .eq. 0) write (6, *) ' nind (n_correlation_functions)=', nind, nw, sizep

#ifdef PARALLEL
!cccccccccccccccccccccccccccccccccccccccccccc
! each process goes from ist to ien walker
    n1 = mod(nw, sizep)
    if (n1 .ne. 0) then
        if (rank .eq. 0) write (6, *) 'nw must be multiple of number of processors!', nw, sizep
        call mpi_finalize(ierr)
        stop
    end if
    ist = rank*(nw/sizep) + 1
    ien = (rank + 1)*(nw/sizep)
    id1 = rank*(nw/sizep)
    in1 = ien - ist + 1
    nw_per_file = in1/size_per_proc
    itmp = nw/sizep
!cccccccccccccccccccccccccccccccccccccccccccc
#else
    ist = 1
    ien = nw
    id1 = 0
    in1 = nw
    nw_per_file = in1/size_per_proc
#endif

    !global vectors: jbra, jbram, jbra_old
    !local vectors: walk, iconfm, iconfmm, amisc, wconf, wconfm, swconf
    allocate (iconfm(nind, nbuf, in1), iconfmm(nind, nbuf, in1))
    allocate (walk(ndim, nel_read, in1), dwalk(ndim, nel_read))
    allocate (jbra(nbuf, nw), jbram(nbuf, nw), jbra_old(nw), jbra_step(nw), jbra_temp(nbuf, in1))
    allocate (wconf(nbuf, in1), wconfm(nbuf, in1))
    allocate (swconf(ipc*in1), seconfd(in1), spsiln(in1), signpsi(in1))
    allocate (wbuf(nbuf), wbufm(nbuf))
    allocate (amisc(nind, in1), opar(nind))
    allocate (ek(0:maxf_r, nind), ebin(0:maxf_r, nind_corrfun), ebin2(0:maxf_r, nind_corrfun), wk(0:maxf_r) &
              , wbin(0:maxf_r))
    ! by E. Coccia (23/6/11)
    if (rankrep .eq. 0) allocate (ek_tot(0:maxf_r, nind), wk_tot(0:maxf_r))

    allocate (psip_for(2*nind_corrfun*(maxf_r + 1) + 4)) ! scratch vector for averaging and writing

    if (ifcorrs .and. pseudorandom) then
        ifsan = .true.
        allocate (sangle(18, nel*in1))
    else
        ifsan = .false.
        allocate (sangle(18, 1))
    end if

    ibin_start = 0
    iconfm = 0.d0
    iconfmm = 0.d0
    walk = 0.d0
    dwalk = 0.d0
    jbra = 0
    jbram = 0
    jbra_old = 0
    jbra_step = 0
    wconf = 0.d0
    wconfm = 0.d0
    wbuf = 0.d0
    wbufm = 0.d0
    amisc = 0.d0
    opar = 0.d0
    ek = 0.d0
    ebin = 0.d0
    ebin2 = 0.d0
    wk = 0.d0
    wbin = 0.d0
    psip_for = 0.d0
    if (rankrep .eq. 0) then
        ! by E. Coccia (23/6/11)
        ek_tot = 0.d0
        wk_tot = 0.d0
    end if

    history_change_send = 0.d0
    history_change_rec = 0.d0
    history_change_onsite = 0.d0
    history_total = 0.d0

    err_read = .false.
    ibin_start = 0

    ! here control the continuation
    if (rank .eq. 0) then
        ! iopt read in the call read_datasmin
        if (iopt .eq. 0) then
            inquire (file='fort.readforward', exist=ife)
            if (.not. ife) then
                write (6, *) 'file fort.1 does not exist'
                write (6, *) 'you cannot continue the readforward.x run'
                write (6, *) 'taking into account only the generations in the last QMC run'
                err_read = .true.
            else
                write (6, *) 'CONTINUING A PREVIOUS readforward.x RUN'
            end if
            ! check also the length here (to be done)
        end if

        if (iopt .eq. 0) then
            if (.not. err_read) then
                ! read from file previously computed values
                open (1, file='fort.readforward', form='unformatted', status='unknown')
                read (1) ibin_start, lbinr
                read (1) ((ebin(i, j), i=0, maxf_r), j=1, nind_corrfun)
                read (1) ((ebin2(i, j), i=0, maxf_r), j=1, nind_corrfun)
                read (1) (wbin(i), i=0, maxf_r)
                close (1)

                if (longio .or. ifcorrs) then
                    do i = 1, ibin_start
                        read (14, *)
                    end do
                end if

            end if

            ! last number of generations in datasmin (ngen)
            ! total number of generations in fort.11 file (ngenc)
            ! requested number of generations in readfowrad.input (ngenr)
#ifdef PARALLEL
            call mpi_bcast(ibin_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(lbinr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

            if (iopt .eq. 0 .and. lbinr .ne. lbin) then
                write (6, *) ' Warning  continuing  readforward with the same bin length', lbinr
                lbin = lbinr
            end if

            ngendone = ibin_start*lbinr*ifreqdump

            if (mod(lbin, ifreqdump) .ne. 0) then
                write (6, *) ' ERROR the code does not work for mod(lbin,ifreqdump) =/0 '
#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if

            write (6, *) 'last number of generations ', ngendone
            write (6, *) 'total number of generations (in fort.11)', ngenc
            write (6, *) 'requested number of generations (in readfoward.input)', ngenr

            if (ngenr .eq. 0 .or. ngenr + ngendone .gt. ngenc) then
                if (ngenr + ngendone .gt. ngenc) write (6, *) ' Warning: Max number of generations in last run=', ngen
                ngenr = ngenc - ngendone
            end if
            skip_read = ngendone/ifreqdump
            ngen = ngenr/ifreqdump
            write (6, *) 'previously dumped generations', skip_read
        else
            if (ngenr .eq. 0 .or. ngenr .gt. ngenc) then
                if (ngenr .gt. ngenc) write (6, *) ' Warning max number of generations =', ngenc
                ngenr = ngenc
            end if

            skip_read = 0
            ngen = ngenr/ifreqdump

        end if

    end if

#ifdef PARALLEL
    call mpi_bcast(ngen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(skip_read, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

    if (lbin .eq. 0) lbin = max(ngen/100, 1)

    j_shift = rank*size_per_proc

#ifdef PARALLEL
    if (io_level .eq. 1) then
        do i = 0, size_per_proc - 1
            j = read_start + i
            call convertdec(i + j_shift, idigit)
            charascratch = trim(scratchpath)//'details.'//&
                 &char(idigit(1))//char(idigit(2))//char(idigit(3))//&
                 &char(idigit(4))//char(idigit(5))//char(idigit(6))
!  !      write(6,*)charascratch
            open (unit=j, file=charascratch, form='unformatted')
            do k = 1, skip_read
                read (j)
            end do
        end do
    else if (io_level .eq. 2) then
        call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'details_SP.all', MPI_MODE_RDONLY, details_SP)
        call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'details_DP.all', MPI_MODE_RDONLY, details_DP)
        ! skipping the first lines during the first read
    end if
#else
    open (read_start, file='fort.12.new', form='unformatted')
    do k = 1, skip_read
        read (read_start)
    end do
#endif

    ng = (ngen - 1)/nbuf + 1
    nrest = ngen - (ng - 1)*nbuf

    if (rank .eq. 0) then
        write (6, *) '***************************************************'
        write (6, *) '# ions', nion
        write (6, *) '# electrons', nel
        write (6, *) '# generations', ngen
        write (6, *) '# walkers', nw
        write (6, *) '# datafiles', size_run
        write (6, *) '# processors', sizep
        write (6, *) 'population bias correction', nbias
        write (6, *) 'bin length (bin_length) =', lbin
        write (6, *) 'initial bin for averages =', ibinit
        write (6, *) 'forward walking propagation', maxf
        write (6, *) 'averages computed every', iskip, 'forwarded iterations'
        write (6, *) 'forwarded iterations written', maxf_r
        write (6, *) '***************************************************'
    end if

    ! end INITIALIZATION

    icount = 0
    ibin = 0

    ! loop over the number of buffers
    ! ng -> number of generations
    do kt = 1, ng
        do ib = 1, nbuf
            wbufm(ib) = wbuf(ib)
        end do
        do iw = ist, ien
            iw_r = iw - id1
            do ib = 1, nbuf
                do i = 1, nind
                    iconfmm(i, ib, iw_r) = iconfm(i, ib, iw_r)
                end do
            end do
        end do
        do iw = ist, ien
            iw_r = iw - id1
            do ib = 1, nbuf
                wconfm(ib, iw_r) = wconf(ib, iw_r)
            end do
        end do

        call copyi4(nwbuf, jbra, jbram)

        if (kt .eq. ng) then
            maxj = nrest
        else
            maxj = nbuf
        end if

        do ib = 1, maxj
            !wbuf is the same for all processes (files), does not need to be communicated
            !#ifdef PARALLEL
            nstart_shift = ist - id1
            nstart = ist
            if (io_level .eq. 1) then
                do j = read_start, read_start + size_per_proc - 1
                    nend_shift = nstart_shift + nw_per_file - 1
                    nend = nstart + nw_per_file - 1
                    if (.not. ifcorrs .and. .not. ifqpwf .and. .not. ifrho_corr) then
                        read (j) wbuf(ib), &
                            (((walk(idim, iel, k), idim=1, ndim), iel=1, nel_read), k=nstart_shift, nend_shift), &
                            (jbra(ib, k), k=nstart, nend), (swconf(k), k=ipc*(nstart_shift - 1) + 1, ipc*nend_shift)
                    elseif (ifcorrs .and. pseudorandom) then

                        read (j) wbuf(ib), &
                                (((walk(idim, iel, k), idim=1, ndim), iel=1, nel_read), k=nstart_shift, nend_shift), &
                                (jbra(ib, k), k=nstart, nend), (swconf(k), k=ipc*(nstart_shift - 1) + 1, ipc*nend_shift), &
                                &(seconfd(k), k=nstart_shift, nend_shift), (spsiln(k), k=nstart_shift, nend_shift), &
                                ((sangle(idim, k), idim=1, 18), k=nel_read*(nstart_shift - 1) + 1, nel_read*nend_shift)
                        if (ipc .eq. 1) then
                            do k = nstart_shift, nend_shift
                                signpsi(k) = 1.0
                                if (swconf(k) .lt. 0) signpsi(k) = -1.0
                            end do
                        else
                            do k = nstart_shift, nend_shift
                                signpsi(k) = swconf(2*k)
                            end do
                        end if
                    else
                        read (j) wbuf(ib), &
                                (((walk(idim, iel, k), idim=1, ndim), iel=1, nel_read), k=nstart_shift, nend_shift), &
                                (jbra(ib, k), k=nstart, nend), (swconf(k), k=ipc*(nstart_shift - 1) + 1, ipc*nend_shift), &
                                &(seconfd(k), k=nstart_shift, nend_shift), (spsiln(k), k=nstart_shift, nend_shift)
                        if (ipc .eq. 1) then
                            do k = nstart_shift, nend_shift
                                signpsi(k) = 1.d0
                                if (swconf(k) .lt. 0) signpsi(k) = -1.0
                            end do
                        else
                            do k = nstart_shift, nend_shift
                                signpsi(k) = swconf(2*k)
                            end do
                        end if
                    end if
                    ntot = nel_read*(nend_shift - nstart_shift + 1)

! test
! write(6,*) 'test before',(((walk(idim, iel, k), idim = 1, ndim), iel = 1,nel_read), k = nstart_shift, nend_shift)
! write(6,*) 'ntot, cellscale, r_offset', ntot, cellscale, r_offset

                    if (.not. if_compute)&
                            &     call MyApplyPBC(walk(1, 1, nstart_shift), ntot, cellscale, r_offset, iespbc)

! write(6,*) 'test after',(((walk(idim, iel, k), idim = 1, ndim), iel = 1,nel_read), k = nstart_shift, nend_shift)

                    nstart_shift = nstart_shift + nw_per_file
                    nstart = nstart + nw_per_file
                end do

! test
!                stop

#ifndef PARALLEL
            end if
#else
            else if (io_level .eq. 2) then
!         if(.not.ifcorrs) then
!            ! to be copied from the following "else" YYYYYY
!         else
            if (details_SP%view == MPI_DATATYPE_NULL) then
                jj = 21*nel_read + 1 + ipc
                ! using disp to skip records
                SP_block_size = jj*in1
                disp = jj*nw*skip_read*4
                details_SP%disp = disp
                !call mpiio_file_create_view(details_SP, jj*nw_per_file,
                !MPI_REAL)
                call mpiio_file_create_view(details_SP, SP_block_size, MPI_REAL)
                if (rank == 0) write (6, *) "details_SP data size per file per iteration", jj*nw_per_file
                if (rank == 0) write (6, *) "details_SP data size per proc per iteration", SP_block_size
                if (rank == 0) write (6, *) "details_SP total data size per iteration", jj*nw
                call mpiio_file_reset_view(details_SP)

                buffer_depth = 1048576/(SP_block_size*4)
                if (buffer_depth > ngen) buffer_depth = ngen
                if (buffer_depth < 1) buffer_depth = 1
                call mpi_bcast(buffer_depth, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                if (rank .eq. 0) write (6, *) "mpiio: buffer_depth", buffer_depth
                if (rank .eq. 0) write (6, *) "mpiio: buffer_size (Byte)", buffer_depth*SP_block_size*4

                allocate (SP_buffer(SP_block_size, buffer_depth))
                buffer_counter = buffer_depth
            end if
            if (details_DP%view == MPI_DATATYPE_NULL) then
                ! using disp to skip records
                DP_block_size = (1 + nw_per_file*2)*size_per_proc
                disp = (1 + nw_per_file*2)*size_run*skip_read*8
                details_DP%disp = disp
                !call mpiio_file_create_view(details_DP, 1+nw_per_file*4,MPI_DOUBLE_PRECISION)
                call mpiio_file_create_view(details_DP, DP_block_size, MPI_DOUBLE_PRECISION)
                if (rank == 0) write (6, *) "details_DP data size per file per iteration", 1 + nw_per_file*2
                if (rank == 0) write (6, *) "details_DP data size per proc per iteration", DP_block_size
                if (rank == 0) write (6, *) "details_DP total data size per iteration", (1 + nw_per_file*2)*size_run
                call mpiio_file_reset_view(details_DP)
                allocate (DP_buffer(DP_block_size, buffer_depth))
            end if

            if (buffer_counter == buffer_depth) then
                if (ngen - (kt - 1)*nbuf - ib + 1 .lt. buffer_depth) then
                    buffer_depth = ngen - (kt - 1)*nbuf - ib + 1
                    if (rank .eq. 0) write (6, *) "mpiio: last buffer depth", buffer_depth
                end if
                call MPI_File_read_all(details_SP%fp, SP_buffer, SP_block_size*buffer_depth, MPI_REAL, status, ierr)
                call MPI_File_read_all(details_DP%fp, DP_buffer, DP_block_size*buffer_depth, MPI_DOUBLE_PRECISION, status, ierr)
                buffer_counter = 0
            end if

            buffer_counter = buffer_counter + 1

            nstart_shift = ist - id1
            nstart = ist
            kk = 0
            jj = 0

            do j = 1, size_per_proc
                nend_shift = nstart_shift + nw_per_file - 1
                nend = nstart + nw_per_file - 1

                do k = nstart_shift, nend_shift
                    do iel = 1, nel_read
                        kk = kk + 3
                        walk(1:3, iel, k) = SP_buffer(kk - 2:kk, buffer_counter)
                    end do
                end do
!              if(pseudorandom) then
                do k = nel_read*(nstart_shift - 1) + 1, nel_read*nend_shift
                    kk = kk + 18
                    sangle(1:18, k) = SP_buffer(kk - 17:kk, buffer_counter)
                end do
!              endif
                jbra(ib, nstart:nend) = int(SP_buffer(kk + 1:kk + nw_per_file, buffer_counter))
                swconf(ipc*(nstart_shift - 1) + 1:ipc*nend_shift) &
                    = SP_buffer(kk + nw_per_file + 1:kk + nw_per_file*(1 + ipc), buffer_counter)
                kk = kk + nw_per_file*(1 + ipc)
                wbuf(ib) = DP_buffer(jj + 1, buffer_counter)
                seconfd(nstart_shift:nend_shift) = DP_buffer(jj + 2:jj + 1 + nw_per_file, buffer_counter)
                spsiln(nstart_shift:nend_shift) = DP_buffer(jj + 2 + nw_per_file:jj + 1 + nw_per_file*2, buffer_counter)
                jj = jj + 1 + nw_per_file*2

                if (ipc .eq. 1) then
                    do k = nstart_shift, nend_shift
                        signpsi(k) = 1.0
                        if (swconf(k) .lt. 0) signpsi(k) = -1.0
                    end do
                else
                    do k = nstart_shift, nend_shift
                        signpsi(k) = swconf(2*k)
                    end do
                end if

                ntot = nel_read*nw_per_file

!test
!      write(6,*) 'test before',(((walk(idim, iel, k), idim = 1, ndim), iel = 1,nel_read), k = nstart_shift, nend_shift)
!      write(6,*) 'ntot, cellscale, r_offset', ntot, cellscale, r_offset

                if (.not. if_compute) call MyApplyPBC(walk(1, 1, nstart_shift), ntot, cellscale, r_offset, iespbc)

!     write(6,*) 'test after',(((walk(idim, iel, k), idim = 1, ndim), iel = 1,nel_read), k = nstart_shift, nend_shift)

                nstart_shift = nstart_shift + nw_per_file
                nstart = nstart + nw_per_file
            end do

! test
!                stop

!         endif
            end if
#endif

            !    dwalk and wconf are in double precision !

            do iw = ist, ien
                iw_r = iw - id1

                do iel = 1, nel_read
                    do idim = 1, ndim
                        dwalk(idim, iel) = walk(idim, iel, iw_r)
                    end do
                end do
                wconf(ib, iw_r) = abs(swconf(ipc*(iw_r - 1) + 1))
                ! here compute the correlation function (local operator)

                if (ifsan) then
                    nelsan = nel_read*(iw_r - 1) + 1
                else
                    nelsan = 1
                end if

                call compute_corr_fun(dwalk, ndim, nel, ell, nel_read, nelup_read&
                        &, ifrho, ifspin, ifkspin, kspin, nrhoind, ddim, dxil&
                        &, ind_offset, density, sdensity, ngrid_l, vell, ifpair, npairind, ngrid_p, drmax&
                        &, paircorr, spaircorr, r_offset_pair, r_offset, vdim, ifsofk&
                        &, nvects, rhok, rhoks, pwmat, kmult, rkcomp, sofk, ifcorrs, corrsamp&
                        &, seconfd(iw_r), spsiln(iw_r), signpsi(iw_r), sangle(1, nelsan)&
                        &, shiftlog, logsamp, noeloc, ipart, sphere_radius, adr_spc, multcell, outofplane, rdf_for_atom)

                logav = logav + logsamp
                countlog = countlog + 1.d0
                ind = 0
                do i = 1, rhotable
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = density(i)
                end do
                do i = 1, pairtable
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = paircorr(i)
                end do
                !if(rank.eq.0) write(6,*) 'check rho index',ind
                indspin0 = ind + 1
                if (ifspin) then
                    do i = 1, spintable
                        ind = ind + 1
                        iconfm(ind, ib, iw_r) = sdensity(i)
                    end do
                    !        ind=ind+1
                    !        iconfm(ind,ib,iw_r)=0.d0
                end if
                do i = 1, spairtable
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = spaircorr(i)
                end do
                !if(rank.eq.0) write(6,*) 'check spin index',ind
                do i = 1, nskind
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = sofk(i, 1)
                end do
                if (fermi_flag) then
                    do kk = 1, 4
                        do i = 1, nskind
                            ind = ind + 1
                            iconfm(ind, ib, iw_r) = sofk(i, kk + 1)
                        end do
                    end do
                end if
                !if(rank.eq.0) write(6,*) 'check sofk index',ind
                do i = 1, ncorrsamp
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = corrsamp(i)
                end do
                !if(rank.eq.0) write(6,*) 'check corrs index',ind
                !(Ye) collect the local value of spin^2
                do i = 1, nspin2
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = spin2_local
                end do
                !if(rank.eq.0) write(6,*) 'check spin2 index',ind
                !(Matteo) collect the STM value
                do i = 1, qpwftable
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = qpwf_image(i)
                end do
                !         if(rank.eq.0) write(6,*) 'check stm index',ind

                ! collect the local value of dipole moment
                do i = 1, ndipole
                    ind = ind + 1
                    iconfm(ind, ib, iw_r) = dipole(i)
                end do
                ! asign the nquad*nquad quadrupole moment components into psip
                do i = 1, nquad
                    do j = 1, nquad
                        ind = ind + 1
                        iconfm(ind, ib, iw_r) = quad(i, j)
                    end do
                end do
                if (ifberry) then
                    do i = 1, 6
                        ind = ind + 1
                        iconfm(ind, ib, iw_r) = berry_exp(i)
                    end do
                end if
                !if(rank.eq.0) write(6,*) 'check dipole index',ind
            end do ! iw (number of walkers)

        end do !ib (block length)

#ifdef PARALLEL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        jbra_temp(:, 1:in1) = jbra(:, ist:ien)
        call mpi_gather(jbra_temp, nbuf*in1, MPI_INTEGER &
                        , jbra, nbuf*in1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        call bcast_integer(jbra, nwbuf, 0, MPI_COMM_WORLD)
#endif
        if (kt .eq. 1) then

            icount = maxk
            do j = maxk + 1, maxj ! loop over the generations in the buffer
                icount = icount + 1
                do k = maxk, 0, -1
                    if (k .eq. maxk) then
                        wsk = 1.d0
                    else
                        if (wsk .le. maxwsk .and. wbuf(j - k - 1) .gt. 1) wsk = wsk*wbuf(j - k - 1)
                    end if
                    !           !SSS added to check
                    !            wsk=1.d0

                    if (k .eq. maxk - nbias) then

                        do iw = 1, nw
                            jbra_old(iw) = iw
                            jbra_step(iw) = iw
                        end do
                        do iw = ist, ien
                            iw_r = iw - id1
                            do i = 1, nind
                                amisc(i, iw_r) = iconfm(i, j - k, iw_r)
                            end do
                        end do
                        k_index = maxf_r

                    end if

                    if (k .le. maxk - nbias) then
                        !           forward walking propagation corr fun

                        if (mod(maxf - k, iskip) .eq. 0) then
                            ! skip calculation averages

                            !########################################################################
                            ! forward walking communication
                            history_total = history_total + dble(in1)

                            do iw = 1, nw
                                kkt = jbra_old(iw)
                                if (jbra_step(iw) .ne. iw) then
                                    iw_r = iw - id1
                                    kkt_r = kkt - id1

#ifdef PARALLEL
                                    ! find the process (nrank) which iw belongs to
                                    ! and find the process (srank) which kkt belongs to
                                    nrank = (iw - 1)/in1
                                    srank = (kkt - 1)/in1

                                    ! belonging to the same process
                                    if (nrank .eq. srank .and. rank .eq. nrank) then
#endif
                                        do kk = 1, nind
                                            amisc(kk, iw_r) = iconfm(kk, j - maxf, kkt_r)
                                        end do
                                        history_change_onsite = history_change_onsite + 1.d0

#ifdef PARALLEL
! belonging to another process: sender
                                    elseif (nrank .ne. srank .and. rank .eq. srank) then

                                        call mpi_send(iconfm(1, j - maxf, kkt_r), nind, MPI_REAL4 &
                                                      , nrank, 1, MPI_COMM_WORLD, ierr)
                                        history_change_send = history_change_send + 1.d0

! belonging to another process: receiver
                                    elseif (nrank .ne. srank .and. rank .eq. nrank) then

                                        call mpi_recv(amisc(1, iw_r), nind, MPI_REAL4 &
                                                      , srank, 1, MPI_COMM_WORLD, status, ierr)
                                        history_change_rec = history_change_rec + 1.d0

                                    end if
#endif

                                end if !jbra(j-k,iw).ne.iw
                            end do !iw=1,nw

                            !#########################################################################

                            !               calculation of opar
                            do kk = 1, nind
                                opar(kk) = 0.d0
                            end do
                            weight = 0.d0
                            do iw = ist, ien
                                iw_r = iw - id1
                                weight = weight + wconf(j - k, iw_r)
                            end do
                            call sgemv('N', nind, in1, 1.0, amisc, nind, wconf(j - k, 1), nbuf, 0.0, opar, 1)

                            do iw = 1, nw
                                jbra_step(iw) = iw
                            end do
                            !                average
                            do kk = 1, nind
                                ek(k_index, kk) = ek(k_index, kk) + opar(kk)*wsk
                            end do
                            wk(k_index) = wk(k_index) + wsk*weight

                            k_index = k_index - 1

                        end if !mod(maxf-k,iskip).eq.0

                        !            update history (jbra_old)
                        do iw = 1, nw
                            jbra_old(iw) = jbra_old(jbra(j - k, iw))
                            jbra_step(iw) = jbra_step(jbra(j - k, iw))
                        end do

                    end if !k.le.maxk-nbias

                end do !end k

                if (mod(icount, lbin) .eq. 0) then
                    ibin = ibin + 1

#ifdef PARALLEL

#ifdef UNREL
                    call reduce_base_real_to(maxf_r + 1, wk(0), psip_for, commrep_mpi, -1)
                    call reduce_base_real_to((maxf_r + 1)*nind, ek(0, 1), psip_for(maxf_r + 2), commrep_mpi, -1)

                    shift = (maxf_r + 1)*(nind + 1)
                    call mpi_allreduce(history_total, psip_for(shift + 1), 1, MPI_DOUBLE_PRECISION&
                         &, MPI_SUM, commrep_mpi, ierr)
                    call mpi_allreduce(history_change_onsite, psip_for(shift + 2), 1&
                         &, MPI_DOUBLE_PRECISION, MPI_SUM, commrep_mpi, ierr)
                    call mpi_allreduce(history_change_send, psip_for(shift + 3), 1&
                         &, MPI_DOUBLE_PRECISION, MPI_SUM, commrep_mpi, ierr)
                    call mpi_allreduce(history_change_rec, psip_for(shift + 4), 1&
                         &, MPI_DOUBLE_PRECISION, MPI_SUM, commrep_mpi, ierr)

#else

                    call reduce_base_real_to(maxf_r + 1, wk(0), psip_for, commrep_mpi, 0)
                    call reduce_base_real_to((maxf_r + 1)*nind, ek(0, 1), psip_for(maxf_r + 2), commrep_mpi, 0)

                    shift = (maxf_r + 1)*(nind + 1)
                    call mpi_reduce(history_total, psip_for(shift + 1), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)
                    call mpi_reduce(history_change_onsite, psip_for(shift + 2), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)
                    call mpi_reduce(history_change_send, psip_for(shift + 3), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)
                    call mpi_reduce(history_change_rec, psip_for(shift + 4), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)
#endif

                    if (rankrep .eq. 0) then
                        do k = 0, maxf_r
                            wk_tot(k) = psip_for(k + 1)
                        end do
                        ind = maxf_r + 1
                        do kk = 1, nind
                            do k = 0, maxf_r
                                ind = ind + 1
                                ek_tot(k, kk) = psip_for(ind)
                            end do
                        end do

                        history_total_s = psip_for(ind + 1)
                        history_change_onsite_s = psip_for(ind + 2)
                        history_change_send_s = psip_for(ind + 3)
                        history_change_rec_s = psip_for(ind + 4)
! by E. Coccia (23/6/11)
                        wk_tot = wk_tot/dble(nw/nk)
                        ek_tot = ek_tot/dble(nw/nk)
                        !ek_tot = ek_tot/real(itmp)
                    end if

#else

                    history_total_s = history_total
                    history_change_onsite_s = history_change_onsite
                    history_change_send_s = history_change_send
                    history_change_rec_s = history_change_rec

#endif
                    ! by E. Coccia (24/6/11)
#ifndef PARALLEL
                    ek_tot = ek/dble(nw/nk)
                    wk_tot = wk/dble(nw/nk)
#endif

                    if (ibin .eq. ibinit - 1 .and. shiftlog .eq. 0.d0) then

                        !     write(6,*) ' ibin here =',ibin,logav,countlog

#ifdef PARALLEL
                        psip_for(1) = logav
                        psip_for(2) = countlog

                        call mpi_allreduce(psip_for, psip_for(3), 2, MPI_DOUBLE_PRECISION &
                                           , MPI_SUM, commrep_mpi, ierr)
                        shiftlog = psip_for(3)/psip_for(4)
#else
                        shiftlog = logav/countlog
#endif

                        if (rank .eq. 0) write (6, *) ' Default value of shiftlog =', shiftlog

                    end if
                    !        write(6,*) ' ibin here after init =',ibin,logav,countlog

                    !if(rank.eq.0) write(*,*) "before computing averages"
                    !only the master computes and writes the final averages
                    if (rankrep .eq. 0) then

                        if (ibin .ge. ibinit) then
                            do k = 0, maxf_r
                                if (wk_tot(k) .ne. 0.d0) then
                                    !              wk_tot(k)=wk_tot(k)/nw
                                    !DEBUG questo wk_tot sicuro va bene?
                                    if (ifspin .and. .not. ifkspin) ek_tot(k, indspin0 + spintable - 1) = &
                                            &sum(abs(ek_tot(k, indspin0:indspin0 + spintable - 2)))

                                    do kk = 1, nind
                                        ebin(k, kk) = ebin(k, kk) + ek_tot(k, kk)
                                        ebin2(k, kk) = ebin2(k, kk) + ek_tot(k, kk)**2/wk_tot(k)

                                    end do

                                    if (.not. ifcorrs .and. k .eq. 0 .and. longio) then
#ifdef __KCOMP
                                        write (14, '(32767(e20.10,1x))') wk_tot(k), (ek_tot(k, kk)/wk_tot(k), kk=1, nind)
#else
                                        write (14, '(1000000(e20.10,1x))') wk_tot(k), (ek_tot(k, kk)/wk_tot(k), kk=1, nind)
#endif
                                    end if

                                    if (ifcorrs .and. k .eq. 0) then
                                        if (ipc .eq. 1) then
                                            write (14, '(5(e20.10,1x))') ek_tot(k, nind - 2)/wk_tot(k), wk_tot(k)&
                                                    &, ek_tot(k, nind - 1)/ek_tot(k, nind), ek_tot(k, nind)&
                                                    &, ek_tot(k, nind - 3)/wk_tot(k)
                                            if (ek_tot(k, nind) .ne. 0.d0) then
                                                ! correlated sampling with reweighting of a bosonic wave function
                                                ! (given by weights w(i))
                                                ! notations
                                                ! < O > =  \sum_i O(i) w(i) / \sum_i w(i)
                                                ! ek(nind-2) = \sum_i E(i) w(i)
                                                ! ek(nind-1) = \sum_i E_new(i) * |\Psi_new(i) / \Psi(i)|^2  w(i)
                                                ! corr_norm = < |\Psi_new / \Psi|^2 >
                                                ! correlated energy  ebin(k,nind+2)
                                                ! < E_new * |\Psi_new / \Psi|^2 > / corr_norm
                                                ! correlated energy difference ebin(k,nind+1)
                                                ! < E - E_new * |\Psi_new / \Psi|^2 / corr_norm > = ek(k,nind-2)
                                                ! - ek(k,nind-1) / corr_norm
                                                corr_norm = ek_tot(k, nind)/wk_tot(k)
                                                ebin(k, nind + 1) = ebin(k, nind + 1) &
                                                                    + ek_tot(k, nind - 2) - ek_tot(k, nind - 1)/corr_norm
                                                ebin2(k, nind + 1) = ebin2(k, nind + 1) &
                                                                     + (ek_tot(k, nind - 2) &
                                                                        - ek_tot(k, nind - 1)/corr_norm)**2/wk_tot(k)
                                                ebin(k, nind + 2) = ebin(k, nind + 2) &
                                                                    + ek_tot(k, nind - 1)/corr_norm
                                                ebin2(k, nind + 2) = ebin2(k, nind + 2) &
                                                                     + (ek_tot(k, nind - 1)/corr_norm)**2/wk_tot(k)
                                            end if
                                        else

                                            write (14, '(6(e20.10,1x))') ek_tot(k, nind - 3)/wk_tot(k), wk_tot(k)&
                                                    &, ek_tot(k, nind - 2)/ek_tot(k, nind - 1), ek_tot(k, nind - 1)&
                                                    &, ek_tot(k, nind - 4)/wk_tot(k), ek_tot(k, nind)/wk_tot(k)
                                            if (ek_tot(k, nind - 1) .ne. 0.d0) then
                                                ! correlated sampling with reweighting of a bosonic wave function
                                                ! (given by weights w(i))
                                                ! notations
                                                ! < O > =  \sum_i O(i) w(i) / \sum_i w(i)
                                                ! ek(nind-2) = \sum_i E(i) w(i)
                                                ! ek(nind-1) = \sum_i E_new(i) * |\Psi_new(i) / \Psi(i)|^2  w(i)
                                                ! corr_norm = < |\Psi_new / \Psi|^2 >
                                                ! correlated energy  ebin(k,nind+2)
                                                ! < E_new * |\Psi_new / \Psi|^2 > / corr_norm
                                                ! correlated energy difference ebin(k,nind+1)
                                                ! < E - E_new * |\Psi_new / \Psi|^2 / corr_norm > = ek(k,nind-2)
                                                ! - ek(k,nind-1) / corr_norm
                                                corr_norm = ek_tot(k, nind - 1)/wk_tot(k)
                                                ebin(k, nind + 1) = ebin(k, nind + 1) + ek_tot(k, nind - 3) &
                                                                    - ek_tot(k, nind - 2)/corr_norm
                                                ebin2(k, nind + 1) = ebin2(k, nind + 1) &
                                                                     + (ek_tot(k, nind - 3) &
                                                                        - ek_tot(k, nind - 2)/corr_norm)**2/wk_tot(k)
                                                ebin(k, nind + 2) = ebin(k, nind + 2) &
                                                                    + ek_tot(k, nind - 2)/corr_norm
                                                ebin2(k, nind + 2) = ebin2(k, nind + 2) &
                                                                     + (ek_tot(k, nind - 2)/corr_norm)**2/wk_tot(k)

                                            end if
                                        end if
                                    end if

                                    wbin(k) = wbin(k) + wk_tot(k)
                                end if
                            end do

                            ibin_av = ibin + ibin_start
                            if (rank .eq. 0) then
                                write (6, *) 'total bin for averaging', ibin_av - ibinit + 1
                                write (6, *) 'frequency of onsite change', history_change_onsite_s/history_total_s
#ifdef PARALLEL
                                write (6, *) 'frequency of send change', history_change_send_s/history_total_s
                                write (6, *) 'frequency of receive change', history_change_rec_s/history_total_s
#endif
                            end if
                            if (longio) then

                                call write_corr_fun(ebin, ebin2, wbin, ibin_av, ibinit, nind_corrfun, maxf_r, ddim &
                                                    , ell, nel, nelup, nrhoind, dxil, ind_offset, psip_for, write_start &
                                                    , ncell, ifrho, ifspin, ifkspin, iespbc, atom_number, datagrid &
                                                    , npairind, dxil_p, ngrid_p, ind_offset_pair, ifpair, iffluct &
                                                    , drmax, r_offset, vdim, nskind, ifsofk, fermi_flag, rkcomp_crys &
                                                    , ndim, simap, ifcorrs, sphere_radius, sec_spc, nshellsp, iioniond &
                                                    , allshells)
                            end if

                            open (1, file='fort.readforward', form='unformatted', status='unknown')
                            write (1) ibin_av - ibinit + 1, lbin
                            write (1) ((ebin(k, kk), k=0, maxf_r), kk=1, nind_corrfun)
                            write (1) ((ebin2(k, kk), k=0, maxf_r), kk=1, nind_corrfun)
                            write (1) (wbin(k), k=0, maxf_r)
                            close (1)

                        end if

                    end if !rank.eq.0

                    ! set to zero for the next bin
                    ! by E. Coccia (22/6/11): initializing total arrays ek_tot and wk_tot
                    do k = 0, maxf_r
                        do kk = 1, nind
                            ek(k, kk) = 0.d0
                        end do
                        wk(k) = 0.d0
                    end do
                    if (rankrep .eq. 0) then
                        do k = 0, maxf_r
                            do kk = 1, nind
                                ek_tot(k, kk) = 0.d0
                            end do
                            wk_tot(k) = 0.d0
                        end do
                    end if

                end if !if(mod(icount,lbin).eq.0)

            end do ! end j

        else

            do j = 1, maxj
                icount = icount + 1
                do k = maxk, 0, -1

                    if (k .eq. maxk) then
                        wsk = 1.d0
                    else
                        if (j - k .gt. 1) then
                            if (wsk .le. maxwsk .and. wbuf(j - k - 1) .gt. 1) wsk = wsk*wbuf(j - k - 1)
                        else
                            if (wsk .le. maxwsk .and. wbuf(j - k - 1 + nbuf) .gt. 1) wsk = wsk*wbufm(j - k - 1 + nbuf)
                        end if
                    end if
                    !        !SSS per controllo
                    !         wsk=1.d0

                    if (k .eq. maxk - nbias) then
                        !           reset jbra_old and k_index

                        do iw = 1, nw
                            jbra_old(iw) = iw
                            jbra_step(iw) = iw
                        end do
                        if (j - k .ge. 1) then
                            do iw = ist, ien
                                iw_r = iw - id1
                                do i = 1, nind
                                    amisc(i, iw_r) = iconfm(i, j - k, iw_r)
                                end do
                            end do
                        else
                            do iw = ist, ien
                                iw_r = iw - id1
                                do i = 1, nind
                                    amisc(i, iw_r) = iconfmm(i, j - k + nbuf, iw_r)
                                end do
                            end do
                        end if
                        k_index = maxf_r

                    end if

                    if (k .le. maxk - nbias) then

                        if (mod(maxf - k, iskip) .eq. 0) then
                            !               forward walking propagation corr fun

                            if (j - k .ge. 1) then
                                ! skip calculation averages

                                !########################################################################
                                ! forward walking communication
                                history_total = history_total + dble(in1)

                                do iw = 1, nw
                                    kkt = jbra_old(iw)
                                    if (jbra_step(iw) .ne. iw) then

                                        iw_r = iw - id1
                                        kkt_r = kkt - id1

#ifdef PARALLEL
                                        ! find the process (nrank) which iw belongs to
! and find the process (srank) which kkt belongs to
                                        nrank = (iw - 1)/in1
                                        srank = (kkt - 1)/in1

! belonging to the same process
                                        if (nrank .eq. srank .and. rank .eq. nrank) then
#endif
                                            if (j - maxf .ge. 1) then
                                                do kk = 1, nind
                                                    amisc(kk, iw_r) = iconfm(kk, j - maxf, kkt_r)
                                                end do
                                            else
                                                do kk = 1, nind
                                                    amisc(kk, iw_r) = iconfmm(kk, j - maxf + nbuf, kkt_r)
                                                end do
                                            end if
                                            history_change_onsite = history_change_onsite + 1.d0

#ifdef PARALLEL
! belonging to another process: sender
                                        elseif (nrank .ne. srank .and. rank .eq. srank) then

                                            if (j - maxf .ge. 1) then
                                                call mpi_send(iconfm(1, j - maxf, kkt_r), nind, MPI_REAL4 &
                                                              , nrank, 1, MPI_COMM_WORLD, ierr)
                                            else
                                                call mpi_send(iconfmm(1, j - maxf + nbuf, kkt_r), nind, MPI_REAL4 &
                                                              , nrank, 1, MPI_COMM_WORLD, ierr)
                                            end if
                                            history_change_send = history_change_send + 1.d0
! belonging to another process: receiver
                                        elseif (nrank .ne. srank .and. rank .eq. nrank) then

                                            call mpi_recv(amisc(1, iw_r), nind, MPI_REAL4 &
                                                          , srank, 1, MPI_COMM_WORLD, status, ierr)
                                            history_change_rec = history_change_rec + 1.d0

                                        end if
#endif

                                    end if !if(kkt.ne.iw)
                                end do
                                !#########################################################################

                                !               calculation of opar
                                do kk = 1, nind
                                    opar(kk) = 0.d0
                                end do
                                weight = 0.d0
                                do iw = ist, ien
                                    iw_r = iw - id1
                                    weight = weight + wconf(j - k, iw_r)
                                end do
                                call sgemv('N', nind, in1, 1.0, amisc, nind, wconf(j - k, 1), nbuf, 0.0, opar, 1)
                                !            if(ifspin)&
                                !     opar(indspin0+spintable)=sum(abs(opar(indspin0:indspin0+spintable-1)))

                            else

                                !########################################################################
                                ! forward walking communication

                                history_total = history_total + dble(in1)

                                do iw = 1, nw
                                    kkt = jbra_old(iw)
                                    if (jbra_step(iw) .ne. iw) then

                                        iw_r = iw - id1
                                        kkt_r = kkt - id1

#ifdef PARALLEL
                                        ! find the process (nrank) which iw belongs to
                                        ! and find the process (srank) which kkt belongs to
                                        nrank = (iw - 1)/in1
                                        srank = (kkt - 1)/in1

                                        ! belonging to the same process
                                        if (nrank .eq. srank .and. rank .eq. nrank) then
#endif

                                            do kk = 1, nind
                                                amisc(kk, iw_r) = iconfmm(kk, j - maxf + nbuf, kkt_r)
                                            end do
                                            history_change_onsite = history_change_onsite + 1.d0

#ifdef PARALLEL
                                            ! belonging to another process: sender
                                        elseif (nrank .ne. srank .and. rank .eq. srank) then

                                            call mpi_send(iconfmm(1, j - maxf + nbuf, kkt_r), nind, MPI_REAL4 &
                                                          , nrank, 1, MPI_COMM_WORLD, ierr)
                                            history_change_send = history_change_send + 1.d0

                                            ! belonging to another process: receiver
                                        elseif (nrank .ne. srank .and. rank .eq. nrank) then

                                            call mpi_recv(amisc(1, iw_r), nind, MPI_REAL4 &
                                                          , srank, 1, MPI_COMM_WORLD, status, ierr)
                                            history_change_rec = history_change_rec + 1.d0

                                        end if
#endif

                                    end if !if(kkt.ne.iw)
                                end do
                                !#########################################################################

                                !               calculation of opar
                                do kk = 1, nind
                                    opar(kk) = 0.d0
                                end do
                                weight = 0.d0
                                do iw = ist, ien
                                    iw_r = iw - id1
                                    weight = weight + wconfm(j - k + nbuf, iw_r)
                                end do
                                call sgemv('N', nind, in1, 1.0, amisc, nind, wconfm(j - k + nbuf, 1), nbuf, 0.0, opar, 1)
                                !            if(ifspin)&
                                !     opar(indspin0+spintable)=sum(abs(opar(indspin0:indspin0+spintable-1)))

                            end if

                            do iw = 1, nw
                                jbra_step(iw) = iw
                            end do

                            do kk = 1, nind
                                ek(k_index, kk) = ek(k_index, kk) + opar(kk)*wsk
                            end do
                            wk(k_index) = wk(k_index) + wsk*weight

                            k_index = k_index - 1

                        end if !mod(maxf-k,iskip).eq.0

                        ! update history (jbra_old)
                        if (j - k .ge. 1) then
                            do iw = 1, nw
                                jbra_old(iw) = jbra_old(jbra(j - k, iw))
                                jbra_step(iw) = jbra_step(jbra(j - k, iw))
                            end do
                        else
                            do iw = 1, nw
                                jbra_old(iw) = jbra_old(jbram(j - k + nbuf, iw))
                                jbra_step(iw) = jbra_step(jbram(j - k + nbuf, iw))
                            end do
                        end if

                    end if !k.le.maxk-nbias
                end do ! end do k

                if (mod(icount, lbin) .eq. 0) then
                    ibin = ibin + 1

#ifdef PARALLEL

#ifdef UNREL
                    call reduce_base_real_to(maxf_r + 1, wk(0), psip_for, commrep_mpi, -1)
                    call reduce_base_real_to((maxf_r + 1)*nind, ek(0, 1), psip_for(maxf_r + 2), commrep_mpi, -1)
                    shift = (maxf_r + 1)*(nind + 1)
                    call mpi_allreduce(history_total, psip_for(shift + 1), 1, MPI_DOUBLE_PRECISION&
                         &, MPI_SUM, commrep_mpi, ierr)
                    call mpi_allreduce(history_change_onsite, psip_for(shift + 2), 1&
                         &, MPI_DOUBLE_PRECISION, MPI_SUM, commrep_mpi, ierr)
                    call mpi_allreduce(history_change_send, psip_for(shift + 3), 1&
                         &, MPI_DOUBLE_PRECISION, MPI_SUM, commrep_mpi, ierr)
                    call mpi_allreduce(history_change_rec, psip_for(shift + 4), 1&
                         &, MPI_DOUBLE_PRECISION, MPI_SUM, commrep_mpi, ierr)
#else

                    call reduce_base_real_to(maxf_r + 1, wk(0), psip_for, commrep_mpi, 0)
                    call reduce_base_real_to((maxf_r + 1)*nind, ek(0, 1), psip_for(maxf_r + 2), commrep_mpi, 0)
                    shift = (maxf_r + 1)*(nind + 1)
                    call mpi_reduce(history_total, psip_for(shift + 1), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)
                    call mpi_reduce(history_change_onsite, psip_for(shift + 2), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)
                    call mpi_reduce(history_change_send, psip_for(shift + 3), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)
                    call mpi_reduce(history_change_rec, psip_for(shift + 4), 1, MPI_DOUBLE_PRECISION &
                                    , MPI_SUM, 0, commrep_mpi, ierr)

#endif

                    if (rankrep .eq. 0) then
                        do k = 0, maxf_r
                            wk_tot(k) = psip_for(k + 1)
                        end do
                        ind = maxf_r + 1
                        do kk = 1, nind
                            do k = 0, maxf_r
                                ind = ind + 1
                                ek_tot(k, kk) = psip_for(ind)
                                !DEBUG
                            end do
                        end do

                        history_total_s = psip_for(ind + 1)
                        history_change_onsite_s = psip_for(ind + 2)
                        history_change_send_s = psip_for(ind + 3)
                        history_change_rec_s = psip_for(ind + 4)
                        ! by E. Coccia (23/6/11)
                        wk_tot = wk_tot/dble(nw/nk)
                        ek_tot = ek_tot/dble(nw/nk)
                        !DEBUG
                        !ek_tot = ek_tot/real(itmp)
!DEBUG come sopra
                    end if

#else

                    history_total_s = history_total
                    history_change_onsite_s = history_change_onsite
                    history_change_send_s = history_change_send
                    history_change_rec_s = history_change_rec

#endif

                    ! by E. Coccia (24/6/11)
#ifndef PARALLEL
                    ek_tot = ek/dble(nw/nk)
                    wk_tot = wk/dble(nw/nk)
                    !DEBUG
#endif
                    !     if(rank.eq.0) write(6,*) ibin, ibinit-1,shiftlog
                    if (ibin .eq. ibinit - 1 .and. shiftlog .eq. 0.d0) then
#ifdef PARALLEL
                        psip_for(1) = logav
                        psip_for(2) = countlog
                        call mpi_allreduce(psip_for, psip_for(3), 2, MPI_DOUBLE_PRECISION &
                                           , MPI_SUM, commrep_mpi, ierr)
                        shiftlog = psip_for(3)/psip_for(4)
#else
                        shiftlog = logav/countlog
#endif
                        if (rank .eq. 0) write (6, *) ' Default value of shiftlog =', shiftlog

                    end if

                    !only the master computes and writes the final averages
                    if (rankrep .eq. 0) then

                        if (ibin .ge. ibinit) then
                            do k = 0, maxf_r
                                if (wk_tot(k) .ne. 0.d0) then
                                    !              wk_tot(k)=wk_tot(k)/nw
                                    if (ifspin .and. .not. ifkspin) ek_tot(k, indspin0 + spintable - 1) = &
                                            &sum(abs(ek_tot(k, indspin0:indspin0 + spintable - 2)))
                                    !DEBUG

                                    do kk = 1, nind
                                        !                ek_tot(k,kk)=ek_tot(k,kk)/nw
                                        ebin(k, kk) = ebin(k, kk) + ek_tot(k, kk)
                                        ebin2(k, kk) = ebin2(k, kk) + ek_tot(k, kk)**2/wk_tot(k)
                                        !DEBUG
                                        !                       write(22,*) ek(k,kk)/wk(k),wk(k)
                                    end do
                                    if (.not. ifcorrs .and. k .eq. 0 .and. longio) then
#ifdef __KCOMP
                                        write (14, '(32767(e20.10,1x))') wk_tot(k), (ek_tot(k, kk)/wk_tot(k), kk=1, nind)
#else
                                        write (14, '(1000000(e20.10,1x))') wk_tot(k), (ek_tot(k, kk)/wk_tot(k), kk=1, nind)
#endif
                                    end if
                                    if (ifcorrs .and. k .eq. 0) then
                                        if (ipc .eq. 1) then
                                            write (14, '(5(e20.10,1x))') ek_tot(k, nind - 2)/wk_tot(k), wk_tot(k)&
                                                    &, ek_tot(k, nind - 1)/ek_tot(k, nind), ek_tot(k, nind)&
                                                    &, ek_tot(k, nind - 3)/wk_tot(k)

                                            if (ek_tot(k, nind) .ne. 0.d0) then
                                                ! correlated sampling with reweighting of a bosonic wave function
                                                ! (given by weights w(i))
                                                ! notations
                                                ! < O > =  \sum_i O(i) w(i) / \sum_i w(i)
                                                ! ek(nind-2) = \sum_i E(i) w(i)
                                                ! ek(nind-1) = \sum_i E_new(i) * |\Psi_new(i) / \Psi(i)|^2  w(i)
                                                ! corr_norm = < |\Psi_new / \Psi|^2 >
                                                ! correlated energy  ebin(k,nind+2)
                                                ! < E_new * |\Psi_new / \Psi|^2 > / corr_norm
                                                ! correlated energy difference ebin(k,nind+1)
                                                ! < E - E_new * |\Psi_new / \Psi|^2 / corr_norm > = ek(k,nind-2)
                                                ! - ek(k,nind-1) / corr_norm
                                                corr_norm = ek_tot(k, nind)/wk_tot(k)
                                                ebin(k, nind + 1) = ebin(k, nind + 1) + ek_tot(k, nind - 2) &
                                                                    - ek_tot(k, nind - 1)/corr_norm
                                                ebin2(k, nind + 1) = ebin2(k, nind + 1) &
                                                                     + (ek_tot(k, nind - 2) &
                                                                        - ek_tot(k, nind - 1)/corr_norm)**2/wk_tot(k)
                                                ebin(k, nind + 2) = ebin(k, nind + 2) &
                                                                    + ek_tot(k, nind - 1)/corr_norm
                                                ebin2(k, nind + 2) = ebin2(k, nind + 2) &
                                                                     + (ek_tot(k, nind - 1)/corr_norm)**2/wk_tot(k)

                                            end if
                                        else
                                            write (14, '(6(e20.10,1x))') ek_tot(k, nind - 3)/wk_tot(k), wk_tot(k)&
                                                    &, ek_tot(k, nind - 2)/ek_tot(k, nind - 1), ek_tot(k, nind - 1)&
                                                    &, ek_tot(k, nind - 4)/wk_tot(k), ek_tot(k, nind)/wk_tot(k)
                                            if (ek_tot(k, nind - 1) .ne. 0.d0) then
                                                ! correlated sampling with reweighting of a bosonic wave function
                                                ! (given by weights w(i))
                                                ! notations
                                                ! < O > =  \sum_i O(i) w(i) / \sum_i w(i)
                                                ! ek(nind-2) = \sum_i E(i) w(i)
                                                ! ek(nind-1) = \sum_i E_new(i) * |\Psi_new(i) / \Psi(i)|^2  w(i)
                                                ! corr_norm = < |\Psi_new / \Psi|^2 >
                                                ! correlated energy  ebin(k,nind+2)
                                                ! < E_new * |\Psi_new / \Psi|^2 > / corr_norm
                                                ! correlated energy difference ebin(k,nind+1)
                                                ! < E - E_new * |\Psi_new / \Psi|^2 / corr_norm > = ek(k,nind-2)
                                                ! - ek(k,nind-1) / corr_norm

                                                corr_norm = ek_tot(k, nind - 1)/wk_tot(k)
                                                ebin(k, nind + 1) = ebin(k, nind + 1) &
                                                                    + ek_tot(k, nind - 3) - ek_tot(k, nind - 2)/corr_norm
                                                ebin2(k, nind + 1) = ebin2(k, nind + 1) &
                                                                     + (ek_tot(k, nind - 3) &
                                                                        - ek_tot(k, nind - 2)/corr_norm)**2/wk_tot(k)
                                                ebin(k, nind + 2) = ebin(k, nind + 2) &
                                                                    + ek_tot(k, nind - 2)/corr_norm
                                                ebin2(k, nind + 2) = ebin2(k, nind + 2) &
                                                                     + (ek_tot(k, nind - 2)/corr_norm)**2/wk_tot(k)
                                            end if
                                        end if ! endif ipc
                                    end if
                                    wbin(k) = wbin(k) + wk_tot(k)
                                end if
                            end do

                            !         write(6,*) 'bin',ibin
                            ibin_av = ibin + ibin_start
                            if (rank .eq. 0) then
                                write (6, *) 'total bin for averaging ', ibin_av - ibinit + 1

                                if (history_change_onsite_s .ne. 0) then
                                    write (6, *) 'frequency of onsite change', history_change_onsite_s/history_total_s
#ifdef PARALLEL
                                    write (6, *) 'frequency of send change', history_change_send_s/history_total_s
                                    write (6, *) 'frequency of receive change', history_change_rec_s/history_total_s

#endif
                                end if
                            end if

                            if (longio) then
                                !DEBUG
                                if (rankrep .eq. 0) write (*, *) "ebin", ebin(1:3, 1)
                                call write_corr_fun(ebin, ebin2, wbin, ibin_av, ibinit &
                                                    , nind_corrfun, maxf_r, ddim, ell, nel, nelup &
                                                    , nrhoind, dxil, ind_offset, psip_for, write_start, ncell, ifrho &
                                                    , ifspin, ifkspin, iespbc &
                                                    , atom_number, datagrid, npairind, dxil_p, ngrid_p, ind_offset_pair &
                                                    , ifpair, iffluct, drmax, r_offset, vdim &
                                                    , nskind, ifsofk, fermi_flag, rkcomp_crys &
                                                    , ndim, simap, ifcorrs, sphere_radius, sec_spc &
                                                    , nshellsp, iioniond, allshells)
                            end if

                            if (decouple_k) then
                                open (1, file='fort.readforwardK'//trim(chara), form='unformatted', status='unknown')
                            else
                                open (1, file='fort.readforward', form='unformatted', status='unknown')
                            end if
                            write (1) ibin_av - ibinit + 1, lbin
                            write (1) ((ebin(k, kk), k=0, maxf_r), kk=1, nind_corrfun)
                            write (1) ((ebin2(k, kk), k=0, maxf_r), kk=1, nind_corrfun)
                            write (1) (wbin(k), k=0, maxf_r)
                            close (1)

                        end if

                    end if !rankrep.eq.0

                    ! set to zero for the next bin
                    ! by E. Coccia (22/6/11): initializing total arrays ek_tot and wk_tot
                    do k = 0, maxf_r
                        do kk = 1, nind
                            ek(k, kk) = 0.d0
                            !ek_tot(k,kk)=0.d0
                        end do
                        wk(k) = 0.d0
                        !wk_tot(k)=0.d0
                    end do
                    ! by E. Coccia (22/6/11): initializing total arrays ek_tot and wk_tot
                    if (rankrep .eq. 0) then
                        do k = 0, maxf_r
                            do kk = 1, nind
                                ek_tot(k, kk) = 0.d0
                            end do
                            wk_tot(k) = 0.d0
                        end do
                    end if
                end if !(mod(icount,lbin).eq.0)
            end do !j=1,maxj
        end if !kt.eq.1
        end do

#ifdef PARALLEL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

        !            calculation error bars
        if (rankrep .eq. 0) then
            nmis = ibin_av - ibinit + 1
            if (rank .eq. 0) then
                write (6, *) 'writing final averages'
                write (6, *) 'bin length', lbin
                write (6, *) 'total # bin considered', nmis
            end if
            !DEBUG
            call write_corr_fun(ebin, ebin2, wbin, ibin_av, ibinit, nind_corrfun, maxf_r, ddim &
                                , ell, nel, nelup, nrhoind, dxil, ind_offset, psip_for, write_start, ncell, ifrho &
                                , ifspin, ifkspin, iespbc, atom_number, datagrid, npairind, dxil_p, ngrid_p &
                                , ind_offset_pair, ifpair, iffluct, drmax, r_offset, vdim, nskind, ifsofk, fermi_flag &
                                , rkcomp_crys, ndim, simap, ifcorrs, sphere_radius, sec_spc, nshellsp, iioniond, allshells)

        end if

        if (ifcorrs .or. ifqpwf .or. ifrho_corr) then

#ifdef PARALLEL
            call mpi_reduce(logav, psip_for(1), 1, MPI_DOUBLE_PRECISION &
                            , MPI_SUM, 0, commrep_mpi, ierr)
            call mpi_reduce(countlog, psip_for(2), 1, MPI_DOUBLE_PRECISION &
                            , MPI_SUM, 0, commrep_mpi, ierr)
#else
            psip_for(1) = logav
            psip_for(2) = countlog
#endif

            if (rank .eq. 0) then
                write (6, *) ' Warning use the number below to determine shiftlog input '
                write (6, *) ' Average log exponent =', (psip_for(1)/psip_for(2))
            end if

        end if

#ifdef PARALLEL
        if (io_level .eq. 1) then
            do j = read_start, read_start + size_per_proc - 1
                close (j)
            end do
        else if (io_level .eq. 2) then
            if (allocated(DP_buffer)) deallocate (DP_buffer)
            if (allocated(SP_buffer)) deallocate (SP_buffer)
            call mpiio_file_close(details_SP)
            call mpiio_file_close(details_DP)
        end if
#else
        close (read_start)
#endif

        if (rankrep .eq. 0 .and. (longio .or. ifcorrs)) close (14)

        if (ifcorrs) then
            call deallocate_all
            if (rankrep .eq. 0) then
                close (10)
            end if
        end if

        if (ifspin2) deallocate (ratiospin)

        !by E. Coccia (22/6/11): deallocate total arrays for MPI
        if (rankrep .eq. 0) deallocate (ek_tot, wk_tot)

        !by E Coccia (22/11/11: deallocate external fields
        if (ext_pot) then
            call deallocate_extpot()
        end if
        if (vdw) then
            call deallocate_vdw()
        end if

#ifdef _OFFLOAD
!$omp end target data
#ifdef _CUSOLVER
        if (ipf .ne. 2) then
            deallocate (dev_dgetrf_workspace)
            deallocate (dev_zgetrf_workspace)
            deallocate (dev_dgetri_workspace)
            deallocate (dev_zgetri_workspace)
        end if
        !$omp end target data
#endif
#endif
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif
        stop

        contains

        subroutine update_nmolfn
            implicit none
            logical forceyes
            integer firstmolt
            real*8 cost

            firstmolt = 1
            if (molecular .ne. 0) then
                firstmolt = nelorb_c - molecular + 1
                !       check if only the molecular diagonal are present
                firstmol = firstmolt
                cost = 0.d0
                !   check odd raws even columns

                !   compute the last relevant molecular orbital
                lastmol = firstmolt - 1 ! It could be that all molec. orb. are not used.
                if (ipc .eq. 2) then
                    do i = firstmolt, nelorb_c - ndiff
                        do j = firstmolt, nelorb_c - ndiff
                            if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0&
                                    &.and. max(i, j) .gt. lastmol) lastmol = max(i, j)
                        end do
                    end do
                else
                    do i = firstmolt, nelorb_c - ndiff
                        do j = firstmolt, nelorb_c - ndiff
                            if (detmat_c(nelorb_c*(i - 1) + j) .ne. 0 .and. max(i, j) .gt. lastmol) lastmol = max(i, j)
                        end do
                    end do
                end if

            elseif (contraction .ne. 0) then
                lastmol = 0
                if (ipc .eq. 2) then
                    do j = 1, nelorb_c
                        do i = 1, nelorb_c
                            if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0 .and. &
                                max(i, j) .gt. lastmol) lastmol = max(i, j)
                        end do
                        do i = nelorb_c + 1, nelcol_c
                            if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0 .and. &
                                j .gt. lastmol) lastmol = j
                        end do
                    end do
                else
                    do j = 1, nelorb_c
                        do i = 1, nelorb_c
                            if (detmat_c(nelorb_c*(i - 1) + j) .ne. 0 .and. max(i, j) .gt. lastmol) lastmol = max(i, j)
                        end do
                        do i = nelorb_c + 1, nelcol_c
                            if (detmat_c(nelorb_c*(i - 1) + j) .ne. 0 .and. j .gt. lastmol) lastmol = j
                        end do
                    end do
                end if
            end if

            if (rank .eq. 0 .and. contraction .ne. 0 .and. lastmol .lt. nelorb_c - ndiff .and. molecular .ne. 0)&
                        & write (6, *) ' Warning estimated last relevant molecular orbital =', lastmol

            if (yesfast .eq. -1) then

                forceyes = .true.

                if ((itestr .ne. -5 .or. (iesup .eq. 0 .and. iessw .eq. 0)) .and. molecular .ne. 0) then
                    yesfast = 1
                    firstmol = firstmolt
                    nmolfn = lastmol - firstmolt + 1

                    cost = 0.d0
                    do i = 1, firstmol - 1
                        do j = 1, firstmol - 1
                            cost = cost + abs(detmat_c((ipc*nelorb_c)*(i - 1) + ipc*(j - 1) + 1))
                            if (ipc .gt. 1) cost = cost + abs(detmat_c((ipc*nelorb_c)*(i - 1) + ipc*j))
                        end do
                    end do

                    if (cost .ne. 0.d0) then
                        yesfast = 2
                        firstmol = 1
                        nmolfn = lastmol
                    end if

                elseif (contraction .ne. 0) then

                    if (ireadmin .eq. 1 .and. membig) then
                        yesfast = 0

                    else
                        yesfast = 2
                        firstmol = 1
                        nmolfn = lastmol
                    end if

                else

                    yesfast = 0

                end if

                if (nelorbh .le. 2*nmolfn .and. symmagp .and. ireadmin .eq. 0 .and. membig) yesfast = 0
                ! the basis is too small to be conv.

            else

                forceyes = .false.

                ! here yesfast in changed to zero if there are no MO.
                ! DO NOT USE yesfast=1 if no MO are present!!

                if (yesfast .eq. 1 .and. molecular .ne. 0) then
                    firstmol = firstmolt
                    nmolfn = lastmol - firstmolt + 1
                elseif (yesfast .eq. 2 .and. contraction .ne. 0) then
                    firstmol = 1
                    nmolfn = lastmol
                else
                    yesfast = 0
                    firstmol = 1
                    nmolfn = lastmol
                end if

            end if

            if (contraction .ne. 0) then
                if (ipc .eq. 1) then
                    call dgemm('N', 'N', nelorbh*ipf, nelcol_c, nelorb_c, 1.d0, mu_c, nelorbh*ipf&
                            &, detmat_c, nelorb_c, 0.d0, projm, nelorbh*ipf)
                else
                    call zgemm('N', 'N', nelorbh*ipf, nelcol_c, nelorb_c, zone, mu_c, nelorbh*ipf&
                            &, detmat_c, nelorb_c, zzero, projm, nelorbh*ipf)
                end if
            end if

            if (rank .eq. 0) write (6, *) ' Chosen yesfast =', yesfast

        end subroutine update_nmolfn

    end program readforward

    subroutine init_variables(nelr, nelupr, nionr, iespbcr, celldmr, rsrr)
        use allio, only: nel, nelup, nion, iespbc, celldm, rs, rion_fast
        implicit none

        integer, intent(out) :: nelr, nelupr, nionr
        real(8), intent(out) :: celldmr(6), rsrr
        logical, intent(out) :: iespbcr

        open (unit=10, file='fort.10', form='formatted', status='unknown')

        call read_fort10(10)

        nelr = nel
        nelupr = nelup
        nionr = nion
        iespbcr = iespbc
        celldmr = celldm
        rsrr = rs

        call deallocate_all
        close (10)

    end subroutine init_variables

    subroutine MyApplyPBC(s, howmany, cellscale, r_offset, iespbc)
        use cell, only: car2cry, s2r
        implicit none
        integer :: i, k, j
        integer, intent(in) :: howmany
        real*8, intent(in) :: cellscale(3), r_offset(3)
        real*8 vecscra(3), scale_fat(3)
        real(4), dimension(3, howmany) :: s
        logical iespbc
        ! Below input x periodic mod_L=cellscale -->
        ! maps the    offset  < x < L+offset mod(L) to
        ! output x' : 0 <x'< L
        s(1, :) = s(1, :) - r_offset(1)
        s(2, :) = s(2, :) - r_offset(2)
        s(3, :) = s(3, :) - r_offset(3)

        if (iespbc) then
            do i = 1, 3
                scale_fat(i) = dsqrt(sum(s2r(:, i)**2))/cellscale(i)
            end do
            do i = 1, howmany
                vecscra(:) = car2cry(:, 1)*s(1, i) + car2cry(:, 2)*s(2, i) + car2cry(:, 3)*s(3, i)
                !     vecscra(:)=s(:,i)
                !     call CartesianToCrystal(vecscra,1)
                vecscra(1) = anint(vecscra(1)/cellscale(1) - 0.5d0)/scale_fat(1)
                vecscra(2) = anint(vecscra(2)/cellscale(2) - 0.5d0)/scale_fat(2)
                vecscra(3) = anint(vecscra(3)/cellscale(3) - 0.5d0)/scale_fat(3)
                s(:, i) = s(:, i) - s2r(:, 1)*vecscra(1) - s2r(:, 2)*vecscra(2) - s2r(:, 3)*vecscra(3)
            end do
            !   s(1,:)=s(1,:)-cellscale(1)*nint(dble(s(1,:)-cellscale(1)/2.d0-r_offset(1))/cellscale(1))-r_offset(1)
            !   s(2,:)=s(2,:)-cellscale(2)*nint(dble(s(2,:)-cellscale(2)/2.d0-r_offset(2))/cellscale(2))-r_offset(2)
            !   s(3,:)=s(3,:)-cellscale(3)*nint(dble(s(3,:)-cellscale(3)/2.d0-r_offset(3))/cellscale(3))-r_offset(3)
        end if
    end subroutine MyApplyPBC

    subroutine DDMyApplyPBC(s, howmany, cellscale, r_offset, iespbc)
        use cell, only: car2cry, s2r
        implicit none
        integer :: i, k, j
        integer, intent(in) :: howmany
        real*8, intent(in) :: cellscale(3), r_offset(3)
        real(8), dimension(3, howmany) :: s
        real*8 vecscra(3), scale_fat(3)
        logical iespbc
        ! Below input x periodic mod_L=cellscale -->
        ! maps the    offset  < x < L+offset mod(L) to
        ! output x' : 0 <x'< L
        !  rion_center mapped to L/2,L/2,L/2
        s(1, :) = s(1, :) - r_offset(1)
        s(2, :) = s(2, :) - r_offset(2)
        s(3, :) = s(3, :) - r_offset(3)

        if (iespbc) then
            do i = 1, 3
                scale_fat(i) = dsqrt(sum(s2r(:, i)**2))/cellscale(i)
            end do
            do i = 1, howmany
                vecscra(:) = car2cry(:, 1)*s(1, i) + car2cry(:, 2)*s(2, i) + car2cry(:, 3)*s(3, i)
                !     vecscra(:)=s(:,i)
                !     call CartesianToCrystal(vecscra,1)
                vecscra(1) = anint(vecscra(1)/cellscale(1) - 0.5d0)/scale_fat(1)
                vecscra(2) = anint(vecscra(2)/cellscale(2) - 0.5d0)/scale_fat(2)
                vecscra(3) = anint(vecscra(3)/cellscale(3) - 0.5d0)/scale_fat(3)
                s(:, i) = s(:, i) - s2r(:, 1)*vecscra(1) - s2r(:, 2)*vecscra(2) - s2r(:, 3)*vecscra(3)
            end do
        end if
    end subroutine DDMyApplyPBC

    subroutine DMyApplyPBC(s, howmany, cellscale)
        use cell, only: car2cry, s2r
        implicit none
        integer :: i, k, j
        integer, intent(in) :: howmany
        real*8, intent(in) :: cellscale(3)
        real*8, dimension(3, howmany) :: s
        real*8 vecscra(3), scale_fat(3)
        do i = 1, 3
            scale_fat(i) = dsqrt(sum(s2r(:, i)**2))/cellscale(i)
        end do
        ! Below input x periodic mod_L=cellscale -->
        ! maps the    0  < x < L mod(L) to
        ! output x' : 0 <x'< L
        do i = 1, howmany
            vecscra(:) = car2cry(:, 1)*s(1, i) + car2cry(:, 2)*s(2, i) + car2cry(:, 3)*s(3, i)
            !     vecscra(:)=s(:,i)
            !     call CartesianToCrystal(vecscra,1)
            vecscra(1) = anint(vecscra(1)/cellscale(1) - 0.5d0)/scale_fat(1)
            vecscra(2) = anint(vecscra(2)/cellscale(2) - 0.5d0)/scale_fat(2)
            vecscra(3) = anint(vecscra(3)/cellscale(3) - 0.5d0)/scale_fat(3)
            s(:, i) = s(:, i) - s2r(:, 1)*vecscra(1) - s2r(:, 2)*vecscra(2) - s2r(:, 3)*vecscra(3)
        end do
    end subroutine DMyApplyPBC

    subroutine copyi4(nw, jbra, jbram)
        implicit none
        integer i, nw
        integer jbra(*), jbram(*)
        do i = 1, nw
            jbram(i) = jbra(i)
        end do
        return
    end subroutine copyi4

