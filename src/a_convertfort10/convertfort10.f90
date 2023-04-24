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

program convertfort10
    use allio
    use convertmod, only: shift_originref, rion_ref
    use constants
    use IO_m
    use sub_comm
    implicit none
    integer*8 indr, mesh, proc8
    integer i, j, k, ind, nbuf, nleft, jion, max_iter&
            &, nelorbc_in, nelorbcu_in, nelorbjc_in, niesd_in, iesdrr_in, imin, imax, dimo     &
            &, nelcolc_in, maximp, nelorbtot, kk, ii, jj, mine, imain, indparmain       &
            &, i_max, iyt, nelorb2, rankn, nprocn, type_lambda        &
            &, norb_min, norb_max, dim_mat, rang, adrcost_in, adrcost, contractionj_in&
            &, molecularc, molecularjc, nion_in, indorb, indpar, indproc, ipf_in, ipc_in, contraction_in&
            &, adrcostc, adrcostc_in, comm_mpi, nbufp, bufdim, inddetmat, dimdistp, nbufrep, nelorbpf, nelorb_cu, ndiff_in
    real*8 x(3), volmesh, overlapsquare, overo, overo_norm, overon_norm &
            &, overodot_norm, overo_norm_inout   &
            &, overlapsquarej, overlapsquarejsz, overoj, overojsz, maxnorm &
            &, normgem, overunpaired, overunpairedc, overold, overoldsz, dnrm2  &
            &, overmax, epsvpot, overojloc, overojszloc, tracematloc&
            &, overojnloc, overojsznloc, tracematnloc, normunpaired, maxpsi, overlap_ref&
            &, overlapnew, over_new, prec
    complex*16 volmeshc, costc
    double precision, allocatable, dimension(:, :) :: buffer, overs, rion_in&
            &, buffer_c, buffero, oversn, overon, mat, oversj, oversnj, overonj        &
            &, over_sav, molecorb, oversnl, overonl, oversnjl, overonjl, psi_out      &
            &, inv_sav, gammamat, umat, buffer_w, mujc_in            &
            &, dmrg_prod
    double precision, allocatable, dimension(:, :, :) :: buff, buffj
    real*8, allocatable, dimension(:) :: detmatc_in, wdist&
            &, atom_number_in, jasmatc_in, jasmatszc_in, vj_in, distp, eigmat&
            &, weightbuf, tik, inf_box, sup_box
    integer, allocatable, dimension(:) :: multpointer, wpsip
    integer, allocatable, dimension(:, :) :: mupointer
    integer, allocatable, dimension(:) :: address
    logical, allocatable, dimension(:) :: occorb, orbcostn_in, occion
    logical iessz_in, eqcontr, nomolecular, nomolecularj, force_real&
            &, checkall, flag_open, yesopen, eqion, overlap, change_contr, change_jas&
            &, add_onebody2det_read

    character(lchlen) :: path, scratchpath
    real*8 :: psilnn, r0, rc(3), scale_unp
    integer :: meshp
    real*8, external :: jastrow_ei
    logical bigram, yespardiag
    real*8, external :: tracemat2, tracemat, tracematc
    integer, dimension(:, :), allocatable :: where_basis
    integer, dimension(:), allocatable :: nbasis_local
    integer :: max_loc_basis
    logical, allocatable, dimension(:) :: optimize

    !########################################################
    namelist /option/ wherescratch, eqion, bigram, symiesup

    namelist /control/ epsdgel, epsvpot, overlap, change_contr, change_jas &
        , yespardiag, epsbas, double_kpgrid, scale_unp, force_real &
        , real_agp, rmax, rmaxj, rmaxinv, max_iter, prec

    namelist /mesh_info/ nbufd, nx, ny, nz, ax, ay, az, add_onebody2det, shift_origin&
            &, shiftx, shifty, shiftz

    !#######################################################

#ifdef PARALLEL
    include 'mpif.h'
    integer skip, status(MPI_STATUS_SIZE), nrank, srank, ithread
#else
    integer status(1), mpi_comm_world
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'convertfort10'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added
#endif
#ifdef PARALLEL
    call mpi_init(ierr)
!     call mpi_init_thread(MPI_THREAD_FUNNELED,ithread,ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nprocn, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rankn, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
!     if(rank.eq.0) write(6,*) ' Initial mpi value of threads',ithread
    comm_mpi = MPI_COMM_WORLD
    ! define also these values for the read_fort10
    rankrep = rank
    commrep_mpi = MPI_COMM_WORLD
    rankcolrep = 0
    commcolrep_mpi = MPI_COMM_WORLD
    commopt_mpi = MPI_COMM_WORLD
    rankopt = rank
#else
    comm_mpi = 0
    rankn = 0
    rankrep = 0
    commrep_mpi = 0
    rankcolrep = 0
    commcolrep_mpi = 0
    nprocn = 1
    commopt_mpi = 0
#endif
    proc8 = nprocn
    !#ifdef __CASO
    !#else
    !      nprocu=1
    !#endif

    !#ifdef PARALLEL
    !
    !       iesscra=.false.
    !       if(rankn.eq.0) then
    !       open(unit=9,file='wherescratch.dat',form='formatted',status='unknown')
    !       read(9,'(A61)',end=1153) charaadd
    !       write(6,*) ' Scratch dir defined in ',charaadd
    !       iesscra=.true.
    ! 1153  continue
    !       endif
    !       if(.not.iesscra) charaadd='scratch.'
    !      call mpi_bcast(charaadd,61,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    !      call mpi_bcast(iesscra,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    !#endif
    eqion = .true.
    bigram = .true.
    symiesup = .true.

    if (rankn .eq. 0) then
        wherescratch = './' ! Default value for output scratch extension
        oldscra = .false.
        iflagerr = 1
        read (5, nml=option, err=115)
        iflagerr = 0
115     continue
    end if
    call checkiflagerr(iflagerr, rankn, 'ERROR reading option')
#ifdef PARALLEL
    call mpi_bcast(wherescratch, 60 + lchlen, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(eqion, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(bigram, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(symiesup, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    call check_scratch(rankn, './', scratchpath)

    if (rankn .eq. 0) then
        open (unit=10, file='fort.10_in', form='formatted', status='unknown')
        open (unit=11, file='outmesh.bin', form='unformatted'               &
                &, status='unknown')
        open (unit=12, file='outmeshj.bin', form='unformatted'               &
                &, status='unknown')
    end if

    epsdgel = 1d-13 ! condition number criterium
    epsbas = 1d-7 ! default value for epsbas
    epsvpot = 0.d0
    scale_unp = 1.d0 ! weight of unpaired.
    force_real = .false.
    flag_open = .false.
    overlap = .false.
    change_contr = .true.
    change_jas = .false.
    yespardiag = .true.
    double_kpgrid = .true.
    rmax = -1.d0
    rmaxj = -1.d0
    rmaxinv = -1.d0
    max_iter = 25000
    prec = 1.d-6
    real_agp = .false.
    if (rankn .eq. 0) then
        iflagerr = 1
        read (5, nml=control, err=116)
        iflagerr = 0
116     continue
    end if
    call checkiflagerr(iflagerr, rankn, 'ERROR reading control')
#ifdef PARALLEL
    call mpi_bcast(epsdgel, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsbas, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsvpot, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(scale_unp, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rmax, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rmaxj, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rmaxinv, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(prec, 1, MPI_DOUBLE_PRECISION, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(max_iter, 1, MPI_INTEGER, 0                   &
   &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(overlap, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(change_contr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(change_jas, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yespardiag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(double_kpgrid, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(force_real, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(real_agp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    if (yespardiag) then
        nprocu = nprocn
    else
        nprocu = 1
    end if
    nproc_diag = nprocu
    call mpi_sub_comm_create(comm_mpi, nproc_diag, sub_comm_diag, ierr)

    !     FIRST PART END OF READ fort.10_in

    call load_fort10in

    if (.not. iespbc) then
        yesopen = .true.
    else
        yesopen = .false.
    end if

    if (rankn .eq. 0) then
        if (overlap .and. .not. symmagp .and. ipc .ne. 1) then
            write (6, *) ' Overlap with spin dependent or complex AGP not implemented yet '
            iflagerr = 1
        end if
    end if
    call checkiflagerr(iflagerr, rankn, 'ERROR reading opt_zeta')

    !     stored info overs  oversj write outmesh outmeshj (if any)
    if (rankn .eq. 0)                                                    &
            &open (unit=10, file='fort.10_out', form='formatted', status='unknown')

    nion_in = nion

    if (eqion) then
        allocate (rion_in(3, nion), atom_number_in(nion), occion(nion))
        atom_number_in = atom_number
        rion_in = rion
    end if
    call deallocate_all

    if (rankn .eq. 0) write (6, *) ' # processor read ', nprocn

    rank = rankn
    nw = nprocn
    nproc = nprocn
    call default_allocate
    add_onebody2det = add_onebody2det_read
    in1 = 1

    if (rank .eq. 0) then
        call read_fort10_fast
        npsa = npsar
    end if
#ifdef PARALLEL
    call mpi_bcast(npsa, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    npsar = npsa
#endif

    if (.not. eqion) then
        pseudofile = 'pseudo_out'
    else
        pseudofile = 'pseudo.dat'
    end if
    rewind (8)
    call read_pseudo
    call read_fort10(10)

    nelorbpf = nelorbh*ipf
    nelorb_cu = nelorb_c
    if (npar_eagp .gt. 0) then
        nelorbpf = nelorbpf + ndiff
        nelorb_cu = nelorb_c + ndiff
    end if

    if (ipj .eq. 2 .and. change_jas) then
        if (rank .eq. 0) &
                &write (6, *) ' Warning change_jas is set to .false., as it does not &
                &work for generic Jastrow  -12/-22'
        change_jas = .false.
    end if

    !     allocate in any event big matrix detmat
    if (allocated(detmat)) deallocate (detmat)
    if (npar_eagp .gt. 0) then
        allocate (detmat(ipc*nelorbpf*nelorbpf))
    else
        allocate (detmat(ipc*nelorbpf*(nelorbpf + ndiff)))
    end if
    detmat = 0.d0

    adrcost = 0
    do i = 1, nelorbjh
        if (orbcostl(i)) adrcost = i
    end do
    if (contractionj .ne. 0) then
        adrcostc = 0
        do i = 1, nelorbj_c
            if (orbcostn(i)) adrcostc = i
        end do
    else
        adrcostc = adrcost
    end if

    !     to avoid memory conflict in case of non contracted
    if (.not. allocated(iesuptransb)) allocate (iesuptransb(1))
    if (.not. allocated(iesuptransbj)) allocate (iesuptransbj(1))

    if (eqion) then

        if (nion_in .ne. nion) then
            if (rank .eq. 0) write (6, *) ' The two wf should have the same number of ions', nion, nion_in

#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop

        end if

        occion = .true.
        !  They have to be in the same order otherwise the basis is generally changed
        do j = 1, nion
            !      do i=1,nion
            if (atom_number(j) .eq. atom_number_in(j) .and. occion(j)) then
                rion(1:, j) = rion_in(:, j)
                occion(j) = .false.
            end if
            !       enddo
        end do
        do i = 1, nion

            if (occion(i)) then

                if (rank .eq. 0) write (6, *) 'The input wf do not match with fort.10_out,&
                        & try eqion=.false. '

#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif
                stop

            end if

        end do

        deallocate (rion_in, atom_number_in, occion)
    end if

    nomolecular = .true.

    if (allocated(distp)) deallocate (distp)
    dimdistp = ipc*(indt + 5)*nelorbh + nelorbh*(indt + 6) + 27*(indt + 1)*nshell + nelorbh
    dimdistp = max(dimdistp, (indt + 5)*nelorbjh + nelorbjh*(indt + 5) + 27*(indt + 1)*nshellj + nelorbjh)
    allocate (distp(dimdistp))
    distp = 0.d0

    if (ipf_in .ne. ipf) then
        if (rankn .eq. 0) write (6, *) ' Pfaffian in/out =/  Normal in/out CASE not implemented:&
                & please transform both in Pfaffians with convertpfaff.x tool'
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif
        stop
    end if
    if (ipc_in .ne. ipc) then
        if (rankn .eq. 0) write (6, *) ' Complex  in/out =/  Real in/out CASE not implemented:&
                & please transform both in complex  with real_to_complex.x  tool'
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif
        stop
    end if

    if (contraction .ne. 0) then

        !     all about  determinant
        if (allocated(buffer)) deallocate (buffer)
        if (allocated(buffero)) deallocate (buffero)
        if (allocated(weightbuf)) deallocate (weightbuf)

        if (mesh .gt. nbufd) then
            allocate (buffer(ipc*nelorb, max(ipc, ipf)*nbufd), buffero(ipc*nelorbc_in, max(ipc, ipf)*nbufd)          &
                    &, weightbuf(nbufd))
        else
            allocate (buffer(ipc*nelorb, max(ipc, ipf)*mesh), buffero(ipc*nelorbc_in, max(ipc, ipf)*mesh)            &
                    &, weightbuf(mesh))
        end if
        buffer = 0.d0
        buffero = 0.d0
        weightbuf = 0.d0

        allocate (oversnl(ipc*nelorbpf, ipc*nelorbpf), overonl(ipc*nelorbcu_in, ipc*nelorbpf))
        oversnl = 0.d0
        overonl = 0.d0

        !         evaluation matrix M
        deallocate (psip)
        allocate (mat(ipc*nelorbpf, max(nelorbpf, nelorbcu_in)))
        allocate (over_sav(ipc*nelorbpf, ipc*nelorbpf)&
                &, umat(ipc*nelorbpf, nelorbpf), eigmat(nelorbpf))
        allocate (inv_sav(ipc*nelorbpf, nelorbpf))
        if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
            allocate (psip(max(nelorbh, nelorbc_in)*max(nelorb, nelorbc_in)))
        else
            allocate (psip(2*ipc*max(nelorbpf, nelorbcu_in)*max(nelorbpf, nelorbcu_in)))
        end if

        psip = 0.d0
        mat = 0.d0
        over_sav = 0.d0
        umat = 0.d0
        eigmat = 0.d0
        inv_sav = 0.d0

        call load_gemz

        overlap_ref = overlapsquare

        !       Now the unpaired orbitals....
        if (nelcol_c .gt. nelorb_c .and. npar_eagp .eq. 0) then
            allocate (oversn(ipc*nelorb*ipf, nelorb*ipf))
            allocate (gammamat(ipc*ndiff, ndiff))
            gammamat = 0.d0
            oversn = 0.d0

            if (ndiff .ne. nelcolc_in - nelorbc_in) then
                if (rank .eq. 0) write (6, *) ' The two wf should  have             &
                        &the same number of unpaired orb.'
                !!!'
#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if

            !       calculation overlap matrix Gamma

            if (ipc .eq. 1) then
                call dgemm_my('N', 'N', nelorbc_in, ndiff, nelorbc_in, 1.d0, overs       &
                        &, nelorbc_in, detmatc_in(nelorbc_in*nelorbc_in + 1), nelorbc_in, 0.d0   &
                        &, psip, nelorbc_in, nprocu, rank, comm_mpi)

                call dgemm_my('T', 'N', ndiff, ndiff, nelorbc_in, 1.d0                 &
                        &, detmatc_in(nelorbc_in*nelorbc_in + 1), nelorbc_in, psip, nelorbc_in, 0.d0&
                        &, gammamat, ndiff, nprocu, rank, comm_mpi)

            else

                call zgemm_my('N', 'N', nelorbc_in, ndiff, nelorbc_in, zone, overs       &
                        &, nelorbc_in, detmatc_in(ipc*nelorbc_in*nelorbc_in + 1), nelorbc_in, zzero   &
                        &, psip, nelorbc_in, nprocu, rank, comm_mpi)

                call zgemm_my('C', 'N', ndiff, ndiff, nelorbc_in, zone                 &
                        &, detmatc_in(ipc*nelorbc_in*nelorbc_in + 1), nelorbc_in, psip, nelorbc_in, zzero&
                        &, gammamat, ndiff, nprocu, rank, comm_mpi)

            end if

            call invsym(ndiff, gammamat, ndiff, info)

            if (ipc .eq. 1) then

                call dgemm_my('T', 'N', ipf*nelorbh, ndiff, nelorbc_in, 1.d0                &
                        &, overonl, nelorbc_in, detmatc_in(nelorbc_in*nelorbc_in + 1)           &
                        &, nelorbc_in, 0.d0, mat, ipf*nelorbh, nprocu, rank, comm_mpi)

                !      now building the new matrix --> oversn

                call dgemm_my('N', 'T', ndiff, ipf*nelorbh, ndiff, 1.d0, gammamat, ndiff       &
                        &, mat, ipf*nelorbh, 0.d0, psip, ndiff, nprocu, rank, comm_mpi)

                call dgemm_my('N', 'N', ipf*nelorbh, ipf*nelorbh, ndiff, 1.d0, mat, ipf*nelorbh   &
                        &, psip, ndiff, 0.d0, oversn, ipf*nelorb, nprocu, rank, comm_mpi)

            else

                call zgemm_my('C', 'N', ipf*nelorbh, ndiff, nelorbc_in, zone                &
                        &, overonl, nelorbc_in, detmatc_in(2*nelorbc_in*nelorbc_in + 1)           &
                        &, nelorbc_in, zzero, mat, ipf*nelorbh, nprocu, rank, comm_mpi)

                !      now building the new matrix --> oversn

                call zgemm_my('N', 'C', ndiff, ipf*nelorbh, ndiff, zone, gammamat, ndiff       &
                        &, mat, ipf*nelorbh, zzero, psip, ndiff, nprocu, rank, comm_mpi)
                call zgemm_my('N', 'N', ipf*nelorbh, ipf*nelorbh, ndiff, zone, mat, ipf*nelorbh&
                        &, psip, ndiff, zzero, oversn, ipf*nelorb, nprocu, rank, comm_mpi)
            end if

            deallocate (psip)
            lwork = ipc*ipf*ipf*nelorb*nelorb + ipf*nelorb
            allocate (psip(lwork))

            psip = 0.d0

            !  NB here we use the umat computed in load_gemz for the inverse of over_sav
            !  with doover=.false. over_sav is not changed
            call dsygv_my(ipc, ipf*nelorbh, oversn, ipf*nelorb, over_sav, ipf*nelorb&
                    &, .false., umat, eigmat, psip, mine, psip(ipf*nelorb + 1), epsdgel, .true., info&
                    &, nprocu, rank, comm_mpi)

            overunpaired = 1.d0
            do i = 1, ndiff
                overunpaired = overunpaired*psip(ipf*nelorbh - i + 1)
            end do

            if (rank .eq. 0)                                                   &
                    &  write (6, *) ' Overlap square unpaired uncontr. =', overunpaired
            if (.not. overlap) then
                do i = 1, ndiff
                    call dcopy(ipc*ipf*nelorbh, oversn(1, ipf*nelorbh - ndiff + i), 1&
                            &, detmat(ipc*nelorbh*ipf*(ipf*nelorbh + i - 1) + 1), 1)
                end do
            end if

            !       stop

            deallocate (oversn)
            ! endif unpaired orbitals
        end if

        if (contraction .gt. 0 .and. .not. overlap .and. change_contr) then
            ! do not do DMRG when calculating overlap !!!!
            !         evaluation of the density matrix put it in oversln
            !     Now the real DMRG Case D= lambda lamda^dag

            if (nelup .gt. ndiff .or. npar_eagp .gt. 0) then

                if (ndiff .ne. 0 .and. npar_eagp .eq. 0) then

                    deallocate (psip)

                    if (ipc .eq. 1) then
                        allocate (psip(ipf*ipf*nelorbh*nelorbh))
                        psip = 0.d0
                        call dgemm_my('N', 'N', ipf*nelorbh, ipf*nelorbh, ipf*nelorbh, 1.d0&
                                &, over_sav, ipf*nelorb, detmat, ipf*nelorbh, 0.d0, psip, ipf*nelorbh, nprocu&
                                &, rank, comm_mpi)
                        normgem = tracemat(ipf*nelorbh, psip, ipf*nelorbh)
                    else
                        allocate (psip(4*ipf*ipf*nelorbh*nelorbh))
                        !  S_L^T \lambda^*=psip

                        call conjmat(ipf*nelorbh, ipf*nelorbh, detmat, ipf*nelorbh)
                        call zgemm_my('T', 'N', ipf*nelorbh, ipf*nelorbh, ipf*nelorbh, zone&
                                &, over_sav, ipf*nelorb, detmat, ipf*nelorbh, zzero, psip, ipf*nelorbh, nprocu&
                                &, rank, comm_mpi)
                        call conjmat(ipf*nelorbh, ipf*nelorbh, detmat, ipf*nelorbh)
                        ! S_R \lambda^T = psip(2nelorbh^2+:)
                        call zgemm_my('N', 'T', ipf*nelorbh, ipf*nelorbh, ipf*nelorbh, zone&
                                &, over_sav(1, ipf*nelorb + 1), ipf*nelorb, detmat, ipf*nelorbh, zzero&
                                &, psip(2*ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh, nprocu, rank, comm_mpi)
                        normgem = tracemat2(ipf*nelorbh, ipf*nelorbh, psip, ipf*nelorbh&
                                &, psip(2*ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh)
                    end if
                    normgem = abs(normgem) ! it is in any case positive
                end if ! endif ndiff

                !  lambda S_R^*
                if (ipc .eq. 1) then
                    call dgemm_my('N', 'N', nelorbpf, nelorbpf, nelorbpf, 1.d0, detmat&
                            &, nelorbpf, over_sav, nelorbpf, 0.d0, mat, nelorbpf, nprocu&
                            &, rank, comm_mpi)
                else
                    call conjmat(nelorbpf, nelorbpf, over_sav(1, nelorbpf + 1), nelorbpf)
                    call zgemm_my('N', 'N', nelorbpf, nelorbpf, nelorbpf, zone, detmat&
                            &, nelorbpf, over_sav(1, nelorbpf + 1), nelorbpf, zzero, mat, nelorbpf&
                            &, nprocu, rank, comm_mpi)
                    call conjmat(nelorbpf, nelorbpf, over_sav(1, nelorbpf + 1), nelorbpf)
                end if

                !  lambda S^*_R lambda^dag
                if (ipc .eq. 1) then
                    call dgemm_my('N', 'T', nelorbpf, nelorbpf, nelorbpf, 1.d0&
                            &, mat, nelorbpf, detmat, nelorbpf, 0.d0, oversnl, nelorbpf, nprocu&
                            &, rank, comm_mpi)
                else
                    call zgemm_my('N', 'C', nelorbpf, nelorbpf, nelorbpf, zone, mat&
                            &, nelorbpf, detmat, nelorbpf, zzero, oversnl, nelorbpf, nprocu&
                            &, rank, comm_mpi)
                end if

                !       Now we sum the density matrix corresponding to the unpaired
                if (ndiff .ne. 0 .and. npar_eagp .eq. 0) then
                    if (ipc .eq. 1) then

                        call dgemm_my('N', 'T', ipf*nelorbh, ipf*nelorbh, ndiff, 1.d0&
                                &, detmat(ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh&
                                &, detmat(ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh, 0.0d0, psip, ipf*nelorbh&
                                &, nprocu, rank, comm_mpi)

                        call dgemm_my('N', 'N', ipf*nelorbh, ipf*nelorbh, ipf*nelorbh, 1.d0&
                                &, over_sav, ipf*nelorb, psip, ipf*nelorbh, 0.d0, mat, ipf*nelorbh, nprocu&
                                &, rank, comm_mpi)

                    else

                        call zgemm_my('N', 'C', ipf*nelorbh, ipf*nelorbh, ndiff, zone&
                                &, detmat(2*ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh&
                                &, detmat(2*ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh, zzero, psip, ipf*nelorbh&
                                &, nprocu, rank, comm_mpi)
                        !
                        call conjmat(ipf*nelorbh, ipf*nelorbh, psip, ipf*nelorbh)
                        call zgemm_my('T', 'N', ipf*nelorbh, ipf*nelorbh, ipf*nelorbh, zone&
                                &, over_sav, ipf*nelorb, psip, ipf*nelorbh, zzero, mat, ipf*nelorbh, nprocu&
                                &, rank, comm_mpi)
                        call conjmat(ipf*nelorbh, ipf*nelorbh, psip, ipf*nelorbh)
                    end if

                    normunpaired = tracematc(ipf*nelorbh, mat, ipf*nelorbh)

                    cost = scale_unp*normgem/normunpaired*dble(ndiff)/dble(neldo)
                    call daxpy(ipc*ipf*ipf*nelorbh*nelorbh, cost, psip, 1, oversnl, 1)
                end if ! endif ndiff

                !      end sum density matrix

            else ! nelup>ndiff (ndiff=nelup and npar_eagp=0 below )
                !            only unpaired orbital contribute to the density matrix

                if (ipc .eq. 1) then
                    call dgemm_my('N', 'T', ipf*nelorbh, ipf*nelorbh, ndiff, 1.d0&
                            &, detmat(ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh&
                            &, detmat(ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh, 0.d0, oversnl, ipf*nelorbh&
                            &, nprocu, rank, comm_mpi)
                else

                    call zgemm_my('N', 'C', ipf*nelorbh, ipf*nelorbh, ndiff, zone&
                            &, detmat(2*ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh&
                            &, detmat(2*ipf*ipf*nelorbh*nelorbh + 1), ipf*nelorbh, zzero, oversnl&
                            &, ipf*nelorbh, nprocu, rank, comm_mpi)

                end if

            end if

            if (force_real .and. ipc .eq. 2) then
                call dscalzero(nelorbpf*nelorbpf, 0.d0, oversnl(2, 1), 2)
                ! NB nelorb=nelorbh always
                call dscalzero(nelorbpf*nelorbpf, 0.d0, over_sav(2, 1), 2)
            end if

            if (allocated(psip)) deallocate (psip)
            allocate (psip(2*ipc*nelorbh*(nelorbh + 1)))
            psip = 0.d0

            allocate (multpointer(nelorb_c), mupointer(nelorbh, nelorb_c))
            allocate (address(nelorb_c), occorb(nelorb_c))
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                allocate (psi_out(nelorbh, nelorbh))
            else
                allocate (psi_out(ipc*nelorbh, 2*nelorbh))
            end if
            psi_out = 0.d0
            if (allocated(mat)) deallocate (mat)
            allocate (mat(ipc*nelorbh, nelorbh))
            mat = 0.d0

            mupointer = 0
            multpointer = 0

            do kk = 1, nelorb_c
                do ii = 1, iesup_c
                    do jj = 1, multranspip(ii)
                        iy = (transpip(jj)%col(ii) - 1)/(ipf*nelorbh) + 1
                        ix = transpip(jj)%col(ii) - (iy - 1)*nelorbh*ipf
                        if (iy .eq. kk) then
                            multpointer(kk) = multpointer(kk) + 1
                            mupointer(multpointer(kk), kk) = ix
                        end if
                    end do
                end do
            end do

            do ii = 1, nelorb_c
                if (multpointer(ii) .eq. nelorbh) nomolecular = .false.
            end do

            occorb = .true.
            do ii = 1, nelorb_c

                if (occorb(ii) .and. ((multpointer(ii) .gt. 1 .and. nomolecular)&
                        &.or. multpointer(ii) .eq. nelorbh)) then
                    ! only for contracted any ti

                    ind = 0
                    do kk = 1, nelorb_c
                        eqcontr = .true.
                        if (multpointer(ii) .ne. multpointer(kk)) then
                            eqcontr = .false.
                        else
                            do jj = 1, multpointer(ii)
                                if (mupointer(jj, ii) .ne. mupointer(jj, kk)) eqcontr = .false.
                            end do
                        end if
                        if (eqcontr) then
                            ind = ind + 1
                            address(ind) = kk
                            occorb(kk) = .false.
                        end if

                    end do ! enddo kk
                    !             now load   the small matrices
                    !            if(multpointer(ii).le.nelorbh) then

                    do i = 1, multpointer(ii)
                        do j = 1, multpointer(ii)
                            psi_out(ipc*(i - 1) + 1:ipc*i, j) = &
                                    &oversnl(ipc*(mupointer(i, ii) - 1) + 1:ipc*mupointer(i, ii), mupointer(j, ii))
                        end do
                    end do

                    psip = 0.d0
                    do i = 1, multpointer(ii)
                        do j = 1, multpointer(ii)
                            psip(nelorbh + ipc*nelorbh*(j - 1) + ipc*(i - 1) + 1:nelorbh + ipc*nelorbh*(j - 1) + ipc*i) = &
                                    &over_sav(ipc*(mupointer(i, ii) - 1) + 1:ipc*mupointer(i, ii), mupointer(j, ii))
                        end do
                    end do

                    call dsygv_my(ipc, multpointer(ii), psi_out, nelorbh, psip(nelorbh + 1), nelorbh&
                            &, .true., umat, eigmat, psip, mine, mat, epsdgel, .false., info, nprocu, rank, comm_mpi)
                    if (rank .eq. 0) then
                        write (6, *) ' DMRG AGP eigs atom=', kiontot(ii)
                        do j = 1, multpointer(ii)
                            write (6, *) j, psip(j)
                        end do
                    end if
                    if (psip(1) .lt. -epsdgel .and. rank .eq. 0) then
                        write (6, *) ' Warning possible   garbage !!!&
                                & Negative DMRG eigenv increase epsdgel '
                        !      write(6,*) 'eigenvalues dmrg =',ii,(psip(j),j=1,multpointer(ii))
                    end if

                    if (ind .le. multpointer(ii)) then
                        do i = 1, ind

                            !             choose the gauge
                            maxpsi = 0.d0
                            imax = 0
                            if (rank .eq. 0) write (6, *) ' dimension =', multpointer(ii), ind
                            do kk = 1, ipc*multpointer(ii)
                                if (abs(psi_out(kk, multpointer(ii) - i + 1)) .gt. maxpsi) then
                                    maxpsi = abs(psi_out(kk, multpointer(ii) - i + 1))
                                    imax = kk
                                end if
                            end do
                            if (psi_out(imax, multpointer(ii) - i + 1) .gt. 0) then
                                do j = 1, multpointer(ii)
                                    mu_c(ipc*(mupointer(j, address(i)) - 1) + 1:ipc*mupointer(j, address(i)), address(i))&
                                            & = psi_out(ipc*(j - 1) + 1:ipc*j, multpointer(ii) - i + 1)
                                end do
                            else
                                do j = 1, multpointer(ii)
                                    mu_c(ipc*(mupointer(j, address(i)) - 1) + 1:ipc*mupointer(j, address(i)), address(i))&
                                            & = -psi_out(ipc*(j - 1) + 1:ipc*j, multpointer(ii) - i + 1)
                                end do
                            end if
                        end do
                    else
                        if (rank .eq. 0) write (6, *) ' Too many contracted orbitals    &
                                &  they are dependent !!! ', ind, multpointer(ii)
#ifdef PARALLEL
                        call mpi_finalize(ierr)
#endif
                        stop
                    end if

                    ! main if to consider the orbital
                end if

                !! enddo for all contracted  ii
            end do

            !      now definition dup_c and dupr
            do ii = 1, iesup_c
                if (multranspip(ii) .ne. 0) then
                    dup_c(ipc*(ii - 1) + 1:ipc*ii) = 0.d0
                    do jj = 1, multranspip(ii)
                        iy = (transpip(jj)%col(ii) - 1)/(ipf*nelorbh) + 1
                        ix = transpip(jj)%col(ii) - (iy - 1)*(ipf*nelorbh)
                        dup_c(ipc*(ii - 1) + 1:ipc*ii) = dup_c(ipc*(ii - 1) + 1:ipc*ii) + mu_c(ipc*(ix - 1) + 1:ipc*ix, iy)
                    end do
                    dup_c(ipc*(ii - 1) + 1:ipc*ii) = dup_c(ipc*(ii - 1) + 1:ipc*ii)/multranspip(ii)
                end if
            end do
            !       Symmetrization
            if (iesupind .gt. 0 .and. symiesup) then
                if (allocated(jbraiesup)) deallocate (jbraiesup)
                if (allocated(dup)) deallocate (dup)
                allocate (dup(ipc*iesupind))
                allocate (jbraiesup(iesup_c))
                jbraiesup = 0
                dup = 0.d0
                ind = 0
                do i = 1, iesupind
                    ind = ind + 1
                    ii = jbraiesup_sav(ind)
                    do j = 1, abs(ii)
                        if (jbraiesup_sav(ind + j) .gt. 0) then
                            jbraiesup(jbraiesup_sav(ind + j)) = i
                        else
                            jbraiesup(-jbraiesup_sav(ind + j)) = -i
                        end if
                    end do
                    ind = ind + abs(ii)
                end do
                call constrbr_complex(iesupind, iesup_c, jbraiesup, dup_c&
                        &, dup, 1, 1)
                call bconstrbr_complex(iesupind, iesup_c, jbraiesup, dup_c, dup)
            end if
            !      now definition dup_c and dupr
            do ii = 1, iesup_c
                if (multranspip(ii) .ne. 0) then
                    do jj = 1, multranspip(ii)
                        iy = (transpip(jj)%col(ii) - 1)/(ipf*nelorbh) + 1
                        ix = transpip(jj)%col(ii) - (iy - 1)*ipf*nelorbh
                        mu_c(ipc*(ix - 1) + 1:ipc*ix, iy) = dup_c(ipc*(ii - 1) + 1:ipc*ii)
                    end do
                end if
            end do

            deallocate (mupointer, multpointer, occorb, psi_out, address)

        end if ! .not.overlap

        deallocate (over_sav, inv_sav, umat, eigmat)
        deallocate (mat)
        deallocate (oversnl)
        deallocate (overonl)
        deallocate (buffer, buffero, weightbuf)

        !  endif contraction ne 0
    end if

    nomolecularj = .true.

    if (change_jas .and. nelorbj .ne. 0 .and. nelorbjc_in .ne. 0&
            &.and. sum(abs(jasmatc_in(:))) .ne. 0) then

        !     all about jastrow

        if (mesh .gt. nbufd) then
            allocate (buffer(nelorbj, nbufd), buffero(nelorbjc_in, nbufd)        &
                    &, weightbuf(nbufd))
            bufdim = nbufd
        else
            allocate (buffer(nelorbj, mesh), buffero(nelorbjc_in, mesh)          &
                    &, weightbuf(mesh))
            bufdim = mesh
        end if
        allocate (oversnjl(nelorbj, nelorbj)                                &
                &, overonjl(nelorbjc_in, nelorbj))
        allocate (over_sav(nelorbj, nelorbj)                                &
                &, umat(nelorbjh, nelorbjh), eigmat(nelorbjh))
        allocate (inv_sav(nelorbj, nelorbj))
        deallocate (psip)
        if (allocated(mat)) deallocate (mat)
        allocate (psip(max(nelorbjh, nelorbjc_in)* &
                &    max(nelorbj, nelorbjc_in)), mat(nelorbjh, max(nelorbjh, nelorbjc_in)))

        call load_jasz
        overlap_ref = overlapsquarej
        if (iessz .and. iessz_in) overlap_ref = overlapsquarej + overlapsquarejsz
        !            hcg=derz
        !            gcg=derzsz

        if (allocated(inv_sav)) deallocate (inv_sav)

        !           ABOVE  the optimization of Z if any
        !           the following does not change input jasmat only

        if (contractionj .gt. numcost) then

            !         evaluation of the density matrix put it in oversnjl

            call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, 1.d0, over_sav&
                    &, nelorbj, jasmat, nelorbjh, 0.d0, mat, nelorbjh, nprocu, rank, comm_mpi)

            call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, 1.d0, jasmat&
                    &, nelorbjh, mat, nelorbjh, 0.d0, oversnjl, nelorbj, nprocu, rank, comm_mpi)

            ! does not make a big effect
            if (iessz .and. iessz_in) then

                cost = abs(overoj/overojsz)

                if (rank .eq. 0) write (6, *) ' Warning changing density matrix for Sz Jastrow weight= ', cost

                !        add to the density matrix the Jastrow sz with the same weight
                call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, 1.d0, over_sav&
                        &, nelorbj, jasmatsz, nelorbjh, 0.d0, mat, nelorbjh, nprocu, rank, comm_mpi)

                call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, cost, jasmatsz&
                        &, nelorbjh, mat, nelorbjh, 1.d0, oversnjl, nelorbj, nprocu, rank, comm_mpi)

                !      call dgemm('N','N',nelorbjh,nelorbjh,nelorbjh,1.d0,jasmatsz&
                !    &,nelorbjh,mat,nelorbjh,0.d0,oversnjl,nelorbj)

            end if

            if (allocated(psip)) deallocate (psip)
            allocate (psip(nelorbjh + nelorbjh*nelorbjh))
            psip = 0.d0

            allocate (multpointer(nelorbj_c), mupointer(nelorbjh, nelorbj_c))
            allocate (address(nelorbj_c), occorb(nelorbj_c))
            allocate (psi_out(nelorbjh, nelorbjh))
            if (allocated(mat)) deallocate (mat)
            allocate (mat(nelorbjh, nelorbjh))
            mat = 0.d0
            psi_out = 0.d0

            !       allocate(occorb(nelorbjh,nelorbjh))
            !       occorb=.true.
            mupointer = 0
            multpointer = 0

            do kk = 1, nelorbj_c
                do ii = 1, npar3body_c
                    do jj = 1, multranspipj(ii)
                        iy = (transpipj(jj)%col(ii) - 1)/nelorbjh + 1
                        ix = transpipj(jj)%col(ii) - (iy - 1)*nelorbjh
                        if (iy .eq. kk) then
                            multpointer(kk) = multpointer(kk) + 1
                            mupointer(multpointer(kk), kk) = ix
                        end if
                    end do
                end do
            end do

            !      among all the orbitals take the one with maximum eigenvalue of
            !      the density matrix (see DMRG)
            molecularjc = 0
            nomolecularj = .true.
            do ii = 1, nelorbj_c
                if (multpointer(ii) .eq. nelorbjh) then
                    nomolecularj = .false.
                    molecularjc = molecularjc + 1
                end if
            end do

            if (iessz .and. molecularj .gt. 0) then
                norb_max = nelorbj_c - molecularj/2
                norb_min = nelorbj_c - molecularj + 1
            elseif (molecularj .gt. 0) then
                norb_max = nelorbj_c
                norb_min = nelorbj_c - molecularj + 1
            else
                norb_max = nelorbj_c
                norb_min = 1
            end if

            occorb = .true.

            do ii = norb_min, norb_max
                if (occorb(ii) .and. ((multpointer(ii) .gt. 1 .and. nomolecularj)  &
                        &.or. multpointer(ii) .eq. nelorbjh)) then
                    ! only for contracted any t
                    ind = 0
                    do kk = norb_min, norb_max
                        !             determine the contracted orbital kk of the same type of ii
                        eqcontr = .true.
                        if (multpointer(ii) .ne. multpointer(kk)) then
                            eqcontr = .false.
                        else
                            do jj = 1, multpointer(ii)
                                if (mupointer(jj, ii) .ne. mupointer(jj, kk)) eqcontr = .false.
                            end do
                        end if
                        if (eqcontr) then
                            ind = ind + 1
                            address(ind) = kk
                            occorb(kk) = .false.
                        end if
                        ! enddo kk
                    end do

                    !           now load   the small matrices

                    do i = 1, multpointer(ii)
                        do j = 1, multpointer(ii)
                            psi_out(i, j) = oversnjl(mupointer(i, ii), mupointer(j, ii))
                            psip(nelorbjh*j + i)                                &
                                    & = over_sav(mupointer(i, ii), mupointer(j, ii))
                        end do
                    end do
                    if (rank .eq. 0) write (6, *) ' DMRG Jas eigs', ii, multpointer(ii)
                    call dsygv_my(1, multpointer(ii), psi_out, nelorbjh, psip(nelorbjh + 1), nelorbjh&
                            &, .true., umat, eigmat, psip, mine, mat, epsdgel, .false., info, nprocu, rank, comm_mpi)
                    if (rank .eq. 0) then
                        write (6, *) ' DMRG Jas eigs atom=', kiontotj(ii)
                        do j = 1, multpointer(ii)
                            write (6, *) j, psip(j)
                        end do
                    end if
                    !          call  eval_molec_epsdgel(multpointer(ii),psip(nelorbjh+1),mat&
                    !    &,psi_out,psip,nelorbjh,epsdgel,1,rank,rank,comm_mpi,0,.true.)
                    if (psip(1) .lt. -epsdgel .and. rank .eq. 0) then
                        write (6, *) ' Warning possible   garbage !!!                 &
                                &    Negative DMRG Jastow eigenv                                   &
                                &       increase epsdgel '
                        write (6, *) 'eigenvalues dmrg =', ii, (psip(j), j=1, multpointer(ii))
                    end if

                    if (ind .le. multpointer(ii)) then

                        do i = 1, ind
                            !             choose the gauge
                            maxpsi = 0.d0
                            imax = 0
                            do kk = 1, multpointer(ii)
                                if (abs(psi_out(kk, multpointer(ii) - i + 1)) .gt. maxpsi) then
                                    maxpsi = abs(psi_out(kk, multpointer(ii) - i + 1))
                                    imax = kk
                                end if
                            end do
                            if (psi_out(imax, multpointer(ii) - i + 1) .gt. 0) then
                                do j = 1, multpointer(ii)
                                    muj_c(mupointer(j, address(i)), address(i))              &
                                            & = psi_out(j, multpointer(ii) - i + 1)
                                end do
                            else
                                do j = 1, multpointer(ii)
                                    muj_c(mupointer(j, address(i)), address(i))              &
                                            & = -psi_out(j, multpointer(ii) - i + 1)
                                end do
                            end if
                        end do

                    else
                        if (rank .eq. 0) write (6, *) ' Too many contracted Jastrow     &
                                & orbitals  they are dependent !!! ', ind, multpointer(ii)
#ifdef PARALLEL
                        call mpi_finalize(ierr)
#endif
                        stop
                    end if

                    ! main if to consider the orbital
                end if

                !! enddo for all contracted  ii
            end do

            if (iessz .and. molecularj .gt. 0) then
                norb_max = nelorbj_c
                norb_min = nelorbj_c - molecularj/2 + 1

                occorb = .true.

                do ii = norb_min, nelorbj_c
                    if (occorb(ii) .and. ((multpointer(ii) .gt. 1 .and. nomolecularj)  &
                            &.or. multpointer(ii) .eq. nelorbjh)) then
                        ! only for contracted any t
                        ind = 0
                        do kk = norb_min, norb_max
                            !             determine the contracted orbital kk of the same type of ii
                            eqcontr = .true.
                            if (multpointer(ii) .ne. multpointer(kk)) then
                                eqcontr = .false.
                            else
                                do jj = 1, multpointer(ii)
                                    if (mupointer(jj, ii) .ne. mupointer(jj, kk)) eqcontr = .false.
                                end do
                            end if
                            if (eqcontr) then
                                ind = ind + 1
                                address(ind) = kk
                                occorb(kk) = .false.
                            end if
                            ! enddo kk
                        end do

                        !           now load   the small matrices

                        !           if(multpointer(ii).lt.nelorbjh) then

                        do i = 1, multpointer(ii)
                            do j = 1, multpointer(ii)
                                psi_out(i, j) = oversnjl(mupointer(i, ii), mupointer(j, ii))
                                psip(nelorbjh*j + i) = over_sav(mupointer(i, ii), mupointer(j, ii))
                            end do
                        end do
                        if (rank .eq. 0) write (6, *) ' DMRG Jas Sz  eigs'
                        call dsygv_my(1, multpointer(ii), psi_out, nelorbjh, psip(nelorbjh + 1), nelorbjh&
                                &, .true., umat, eigmat, psip, mine, mat, epsdgel, .false., info, nprocu, rank, comm_mpi)
                        !          call  eval_molec_epsdgel(multpointer(ii),psip(nelorbjh+1),mat&
                        !    &,psi_out,psip,nelorbjh,epsdgel,1,rank,rank,comm_mpi,0,.true.)
                        if (psip(1) .lt. -epsdgel .and. rank .eq. 0) then
                            write (6, *) ' Warning possible   garbage !!!                 &
                                    &    Negative DMRG Jastow eigenv                                   &
                                    &       increase epsdgel '
                            write (6, *) 'eigenvalues dmrg =', ii, (psip(j), j=1, multpointer(ii))
                        end if

                        if (ind .le. multpointer(ii)) then

                            do i = 1, ind
                                !             choose the gauge
                                if (psi_out(1, multpointer(ii) - i + 1) .gt. 0) then
                                    do j = 1, multpointer(ii)
                                        muj_c(mupointer(j, address(i)), address(i))              &
                                                & = psi_out(j, multpointer(ii) - i + 1)
                                    end do
                                else
                                    do j = 1, multpointer(ii)
                                        muj_c(mupointer(j, address(i)), address(i))              &
                                                & = -psi_out(j, multpointer(ii) - i + 1)
                                    end do
                                end if
                            end do

                        else
                            if (rank .eq. 0) write (6, *) ' Too many contracted Jastrow     &
                                    & orbitals  they are dependent !!! ', ind, multpointer(ii)
#ifdef PARALLEL
                            call mpi_finalize(ierr)
#endif
                            stop
                        end if

                        ! main if to consider the orbital
                    end if

                    !! enddo for all contracted  ii
                end do

            end if ! molecularj and iessz

            !      now definition dup_c and dupr  NB no change of Z exponents in dup
            do ii = 1, npar3body_c
                if (multranspipj(ii) .ne. 0) then
                    vju_c(ii) = 0.d0
                    do jj = 1, multranspipj(ii)
                        iy = (transpipj(jj)%col(ii) - 1)/nelorbjh + 1
                        ix = transpipj(jj)%col(ii) - (iy - 1)*nelorbjh
                        vju_c(ii) = vju_c(ii) + muj_c(ix, iy)
                    end do
                    vju_c(ii) = vju_c(ii)/multranspipj(ii)
                    do jj = 1, multranspipj(ii)
                        iy = (transpipj(jj)%col(ii) - 1)/nelorbjh + 1
                        ix = transpipj(jj)%col(ii) - (iy - 1)*nelorbjh
                        muj_c(ix, iy) = vju_c(ii)
                    end do
                end if
            end do

            deallocate (mupointer, multpointer, occorb, psi_out, address)
            ! endif contractionj
        end if

        deallocate (oversnjl)
        deallocate (overonjl)
        deallocate (buffer, buffero, weightbuf)
        deallocate (over_sav, umat, eigmat)

        !   endif contraction update oversnjl overonjl
    end if

    !       END  EVALUATION BEST NEW MATRICES detmat jasmat jasmatsz
    !       in the uncontracted basis
    !     all about  determinant
    if (allocated(oversn)) deallocate (oversn)
    if (allocated(overon)) deallocate (overon)
    if (allocated(over_sav)) deallocate (over_sav)
    if (allocated(inv_sav)) deallocate (inv_sav)
    if (allocated(umat)) deallocate (umat, eigmat)
    if (allocated(psip)) deallocate (psip)
    allocate (psip(2*ipc*max(nelorb_cu, nelorbcu_in)*max(nelorb_cu, nelorbcu_in)))
    psip = 0.d0

    if (mesh .gt. nbufd) then
        allocate (buffer(ipc*nelorbh, max(ipc, ipf)*nbufd)&
                &, buffero(ipc*nelorbc_in, max(ipc, ipf)*nbufd), weightbuf(nbufd))
        if (contraction .ne. 0) allocate (buffer_c(ipc*nelorb_c, max(ipc, ipf)*nbufd))
        bufdim = nbufd
    else
        allocate (buffer(ipc*nelorbh, max(ipc, ipf)*mesh), buffero(ipc*nelorbc_in, max(ipc, ipf)*mesh)           &
                &, weightbuf(mesh))
        if (contraction .ne. 0) allocate (buffer_c(ipc*nelorb_c, ipc*mesh))
        bufdim = mesh
    end if
    nbufp = bufdim + 1
    buffer = 0.d0
    buffero = 0.d0
    weightbuf = 0.d0

    if (contraction .ne. 0) buffer_c = 0.d0

    allocate (oversn(ipc*nelorb_cu, ipc*nelorb_cu), overon(ipc*nelorbcu_in, ipc*nelorb_cu))
    allocate (over_sav(ipc*nelorb_cu, ipc*nelorb_cu)&
            &, umat(ipc*nelorb_cu, nelorb_cu), eigmat(nelorb_cu))

    if (ipc .eq. 1) then
        allocate (inv_sav(nelorb_cu, nelorb_cu))
        inv_sav = 0.d0
    end if
    oversn = 0.d0
    overon = 0.d0
    if (npar_eagp .gt. 0) then
        do k = 1, max(ndiff, ndiff_in)
            if (k .le. ndiff) oversn(ipc*(k - 1) + 1 + ipc*nelorb_c, k + nelorb_c) = 1.d0
            if (k .le. min(ndiff, ndiff_in)) &
                    &overon(ipc*(k - 1) + 1 + ipc*nelorbc_in, k + nelorb_c) = 1.d0
        end do
    end if

    over_sav = 0.d0
    umat = 0.d0
    eigmat = 0.d0

    if (.not. bigram) then
#ifdef PARALLEL
        rewind (100)
!     read(100)
#else
        rewind (11)
        ! read first useless record
        read (11)
#endif
    end if
    oversn = 0.d0
    overon = 0.d0
    iflagnorm = 3
    ind = 0
    indtot = 0
    indproc = 0
    !       indr=0
    nbuf = 0
#ifdef _OFFLOAD
!$omp target data map(oversn,overon) map(to:mu_c) map(alloc:buffer,buffero,buffer_c)
#endif
    do k = 1, nz
        !       x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
        do j = 1, ny
            !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
            do i = 1, nx
                indtot = indtot + 1
                !           indr=indr+1
                !           if(indr-(indr/proc8)*proc8.eq.rank) then
                if (indproc .eq. rankn) then
                    ind = ind + 1
                    !           x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)

                    x(:) = (-(nz + 1)/2.d0 + k)*az*at(:, 3) + (-(ny + 1)/2.d0 + j)*ay*at(:, 2) &
                           + (-(nx + 1)/2.d0 + i)*ax*at(:, 1) + rion_ref(:)
                    call upnewwf(1, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                            &, dupr, zetar, rion, distp, buffer(1, ind), nelorbh, nion, kion            &
                            &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                            &, indpar_tab, indorb_tab, indshell_tab, .true.)
                    if (ipc .eq. 2) then
                        call upnewwf(1, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                &, dupr, zetar, rion, distp, buffer(1, bufdim + ind), nelorbh, nion, kion            &
                                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                &, indpar_tab, indorb_tab, indshell_tab, .false.)
                    end if

                    if (add_onebody2det) then
                        psilnn = -scale_one_body
                        do jj = 1, nion
                            if (iespbc) then
                                rc(:) = x(:) - rion(:, jj)
                                call CartesianToCrystal(rc, 1)
                                do kk = 1, 3
                                    rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                                end do
                                r0 = norm_metric(rc, metric)
                            else
                                rc(:) = (x(:) - rion(:, jj))*costz(jj)
                                r0 = dsqrt(sum(rc(:)**2))
                            end if
                            psilnn = psilnn - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                        end do
                        buffer(1:ipc*nelorbh, ind) = buffer(1:ipc*nelorbh, ind)*dexp(psilnn)
                        if (ipc .eq. 2) buffer(1:ipc*nelorbh, bufdim + ind) = buffer(1:ipc*nelorbh, bufdim + ind)*dexp(psilnn)
                    end if

                    call upvpot_ei(x, zetar, rion, weightbuf(ind), nion        &
                            &, LBox, epsvpot)
                    call dscal(ipc*nelorbh, weightbuf(ind), buffer(1, ind), 1)
                    if (ipc .eq. 2) call dscal(2*nelorbh, weightbuf(ind), buffer(1, bufdim + ind), 1)

                end if ! endif indproc

                if (indtot .eq. nbufrep .or. (i .eq. nx .and. j .eq. ny .and. k .eq. nz)) then
#ifdef _OFFLOAD
!$omp target update to (buffer)
#endif
                    if (contraction .eq. 0) then
                        if (ipc .eq. 1) then
                            call dgemm_('N', 'T', nelorbh, nelorbh, ind, volmesh, buffer  &
                                       &, nelorbh, buffer, nelorbh, 1.d0, oversn, nelorb_cu)
                        else
                            call zgemm_('N', 'C', nelorbh, nelorbh, ind, volmeshc, buffer  &
                                    &, nelorbh, buffer, nelorbh, zone, oversn, nelorb_cu)
                            if (ipf .eq. 2) then
                                call zgemm_('N', 'C', nelorbh, nelorbh, ind, volmeshc, buffer(1, nbufp)  &
                                        &, nelorbh, buffer(1, nbufp), nelorbh, zone, oversn(2*nelorbh + 1, nelorbh + 1), nelorb_cu)
                            else
                                call zgemm_('N', 'C', nelorb_c, nelorb_c, ind, volmeshc, buffer(1, nbufp)  &
                                        &, nelorbh, buffer(1, nbufp), nelorbh, zone, oversn(1, nelorb_c + 1), nelorb_cu)
                            end if
                        end if
                    else
                        !CCC
                        if (ipc .eq. 1) then
                            call dgemm_('T', 'N', nelorb_c, ind, nelorbh, 1.d0, mu_c&
                                    &, ipf*nelorbh, buffer, nelorbh, 0.d0, buffer_c, nelorb_c)
                            call dgemm_('N', 'T', nelorb_c, nelorb_c, ind, volmesh&
                                    &, buffer_c, nelorb_c, buffer_c, nelorb_c, 1.d0, oversn, nelorb_cu)
                            if (ipf .eq. 2) then
                                call dgemm_('T', 'N', nelorb_c, ind, nelorbh, 1.d0, mu_c(nelorbh + 1, 1)&
                                        &, 2*nelorbh, buffer, nelorbh, 0.d0, buffer_c(1, nbufp), nelorb_c)
                                call dgemm_('N', 'T', nelorb_c, nelorb_c, ind, volmesh&
                                        &, buffer_c(1, nbufp), nelorb_c, buffer_c(1, nbufp), nelorb_c, 1.d0&
                                        &, oversn, nelorb_cu)
                            end if
                        else
                            call zgemm_('T', 'N', nelorb_c, ind, nelorbh, zone, mu_c&
                                    &, ipf*nelorbh, buffer, nelorbh, zzero, buffer_c, nelorb_c)
                            call zgemm_('N', 'C', nelorb_c, nelorb_c, ind, volmeshc&
                                    &, buffer_c, nelorb_c, buffer_c, nelorb_c, zone, oversn, nelorb_cu)
                            if (ipf .eq. 2) then
                                call zgemm_('T', 'N', nelorb_c, ind, nelorbh, zone, mu_c(2*nelorbh + 1, 1)&
                                        &, 2*nelorbh, buffer(1, nbufp), nelorbh, zzero, buffer_c(1, nbufp), nelorb_c)
                                call zgemm_('N', 'C', nelorb_c, nelorb_c, ind, volmeshc&
                                        &, buffer_c(1, nbufp), nelorb_c, buffer_c(1, nbufp), nelorb_c, zone&
                                        &, oversn, nelorb_cu)
                            else
                                call zgemm_('T', 'N', nelorb_c, ind, nelorbh, zone, mu_c&
                                            &, nelorbh, buffer(1, nbufp), nelorbh, zzero, buffer_c(1, nbufp), nelorb_c)
                                call zgemm_('N', 'C', nelorb_c, nelorb_c, ind, volmeshc, buffer_c(1, nbufp) &
                                            &, nelorb_c, buffer_c(1, nbufp), nelorb_c, zone, oversn(1, nelorb_c + 1), nelorb_c)
                            end if
                        end if

                    end if

                    nbuf = nbuf + 1
                    if (bigram) then
                        buffero(:, :) = buff(:, :, nbuf)
                    else

#ifdef PARALLEL
                        read (100) buffero
#else
                        read (11) buffero
#endif
                    end if

                    do ii = 1, ind
                        call dscal(ipc*nelorbc_in, weightbuf(ii), buffero(1, ii), 1)
                        if (ipc .eq. 2 .or. ipf .eq. 2) call dscal(ipc*nelorbc_in, weightbuf(ii), buffero(1, ii + bufdim), 1)
                    end do
#ifdef _OFFLOAD
!$omp target update to (buffero)
#endif
                    if (contraction .eq. 0) then
                        !CCC
                        if (ipc .eq. 1) then
                            call dgemm_('N', 'T', nelorbc_in, nelorbh, ind, volmesh       &
                                    &, buffero, nelorbc_in, buffer, nelorbh, 1.d0, overon, nelorbcu_in)
                            if (ipf .eq. 2) then
                                call dgemm_('N', 'T', nelorbc_in, nelorbh, ind, volmesh, buffero(1, nbufp) &
                                            &, nelorbc_in, buffer, nelorbh, 1.d0, overon(1, nelorbh + 1), nelorbcu_in)
                            end if
                        else
                            call zgemm_('N', 'C', nelorbc_in, nelorbh, ind, volmeshc       &
                                    &, buffero, nelorbc_in, buffer, nelorbh, zone, overon, nelorbcu_in)
                            if (ipf .eq. 2) then
                                call zgemm_('N', 'C', nelorbc_in, nelorbh, ind, volmeshc, buffero(1, nbufp)&
                                        &, nelorbc_in, buffer(1, nbufp), nelorbh, zone, overon(1, nelorbh + 1), nelorbcu_in)
                            else
                                call zgemm_('N', 'C', nelorbc_in, nelorb_c, ind, volmeshc, buffero(1, nbufp) &
                                            &, nelorbc_in, buffer(1, nbufp), nelorbh, zone, overon(1, nelorb_c + 1), nelorbc_in)
                            end if
                        end if
                    else
                        if (ipc .eq. 1) then
                            call dgemm_('N', 'T', nelorbc_in, nelorb_c, ind, volmesh&
                                    &, buffero, nelorbc_in, buffer_c, nelorb_c, 1.d0, overon, nelorbcu_in)
                            if (ipf .eq. 2) then
                                call dgemm_('N', 'T', nelorbc_in, nelorb_c, ind, volmesh       &
                                        &, buffero(1, nbufp), nelorbc_in, buffer_c(1, nbufp), nelorb_c, 1.d0&
                                        &, overon, nelorbcu_in)
                            end if
                        else
                            call zgemm_('N', 'C', nelorbc_in, nelorb_c, ind, volmeshc      &
                                    &, buffero, nelorbc_in, buffer_c, nelorb_c, zone, overon, nelorbcu_in)
                            if (ipf .eq. 2) then
                                call zgemm_('N', 'C', nelorbc_in, nelorb_c, ind, volmeshc, buffero(1, nbufp) &
                                            &, nelorbc_in, buffer_c(1, nbufp), nelorb_c, zone, overon, nelorbcu_in)
                            else
                                call zgemm_('N', 'C', nelorbc_in, nelorb_c, ind, volmeshc, buffero(1, nbufp) &
                                            &, nelorbc_in, buffer_c(1, nbufp), nelorb_c, zone, overon(1, nelorb_c + 1), nelorbc_in)
                            end if
                        end if
                    end if
                    ind = 0
                    indtot = 0
                end if ! endif load buf
                indproc = indproc + 1
                if (indproc .eq. nprocn) indproc = 0
            end do
        end do
    end do
#ifdef _OFFLOAD
!$omp end target data
#endif
    !      collect overon oversn
    if (ipf .eq. 2) then
        if (ipc .eq. 1 .and. contraction .eq. 0) then
            !CCC
            !        overon(nelorbc_in/2+1:nelorbc_in,nelorb_c/2+1:nelorb_c)=overon(1:nelorbc_in/2,1:nelorb_c/2)
            oversn(nelorb_c/2 + 1:nelorb_c, nelorb_c/2 + 1:nelorb_c) = oversn(1:nelorb_c/2, 1:nelorb_c/2)
        elseif (ipc .eq. 2) then
            do k = nelorb_c + 1, nelorb_cu
                oversn(2*k - 1, k) = 1.d0
                overon(2*(k - nelorb_c + nelorbc_in) - 1, k) = 1.d0
            end do
            overon(1:2*nelorbcu_in, nelorb_cu + 1:2*nelorb_cu) = overon(1:2*nelorbcu_in, 1:nelorb_cu)
            oversn(1:2*nelorb_cu, nelorb_cu + 1:2*nelorb_cu) = oversn(1:2*nelorb_cu, 1:nelorb_cu)
        end if
        if (ipc .eq. 1) then
            do k = nelorb_c + 1, nelorb_cu
                oversn(k, k) = 1.d0
                overon(k - nelorb_c + nelorbc_in, k) = 1.d0
            end do
        end if

    end if
#ifdef PARALLEL
    nelorb2 = ipc*ipc*nelorb_cu*nelorb_cu
    call reduce_base_real(nelorb2, oversn, MPI_COMM_WORLD, -1)
    nelorb2 = ipc*ipc*nelorbcu_in*nelorb_cu
    call reduce_base_real(nelorb2, overon, MPI_COMM_WORLD, -1)
#endif
    if (ipc .eq. 2) then
        call conjmat(nelorb_cu, 2*nelorb_cu, oversn, nelorb_cu)
        call conjmat(nelorbcu_in, 2*nelorb_cu, overon, nelorbcu_in)
    end if

    over_sav = oversn
    !        write(6,*) ' diagonal overlaps '
    !          do i=1,2*nelorbh
    !          write(6,*) i,oversn(2*i-1,i),oversn(2*i,i)
    !          enddo
    !        write(6,*) ' overon ',nelorbc_in
    !          do i=1,nelorbc_in
    !             do j=1,2*nelorbh
    !             write(6,*) i,j,overon(2*i-1,j),overon(2*i,j)
    !             enddo
    !          enddo
    !          stop

    if (real_agp .or. rmax .gt. 0) then
        if (ipf .eq. 2) then
            write (6, *) ' ERROR not implemented option real_agp/rmax with  pfaffian !!! '
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        !CG
        !    definire optimize
        allocate (optimize(ipc*nelorb_c*nelorb_c))
        optimize = .false.
        call eval_iond(iond, rion, nion, LBox, psip, iond_cart)
        do ix = 1, nelorb_c
            do iy = 1, nelorb_c
                if (iond(kiontot(ix) + (kiontot(iy) - 1)*nion) .le. rmax .or. rmax .lt. 0) then
                    ind = ipc*(ix - 1) + 1 + ipc*nelorb_c*(iy - 1)
                    optimize(ind) = .true.
                    if (ipc .eq. 2 .and. .not. real_agp) optimize(ind + 1) = .true.
                end if
            end do
        end do
        if (real_agp .and. ipc .eq. 1) then
            do i = 1, nnozero_c
                if (sjbradet(i)) optimize(nozero_c(i)) = .false.
            end do
        end if
        if (rank .eq. 0) then
            !write(6,*) ' Optimized vector '
            !do i=1,ipc*nelorb_c*nelorb_c
            !   write(6,*) i,optimize(i)
            !enddo
        end if
        if (allocated(mat)) deallocate (mat)
        allocate (mat(ipc*nelorb_c, nelorb_c))
        mat = 0.d0
        do i = 1, nelorb_c
            mat(ipc*(i - 1) + 1, i) = 1.d0
        end do
        !       if(rank.eq.0) then
        !          write(6,*) ' overo =',overo
        !          write(6,*) ' detmatc_in=',sum(abs(detmatc_in(:)))
        !          write(6,*) ' over left =',sum(abs(oversn(:,1:nelorb_c)))
        !          if(ipc.eq.2) then
        !             write(6,*) ' over right =',sum(abs(oversn(:,nelorb_c+1:2*nelorb_c)))
        !          endif
        !          write(6,*) ' overon left =',sum(abs(overon(:,1:nelorb_c)))
        !          if(ipc.eq.2) then
        !             write(6,*) ' overon right =',sum(abs(overon(:,nelorb_c+1:2*nelorb_c)))
        !          endif

        !         endif
        type_lambda = 1
        if ((ipc .eq. 2) .and. (symmagp .eqv. .true.) .and. (yes_hermite .eqv. .false.)) type_lambda = 2
        if (symmagp .eqv. .false.) type_lambda = 3

        if ((yespardiag .eqv. .false.) .and. (nprocu .eq. 1)) then
            if (ipc .eq. 2) then
                call max_ovlp(nelorbc_in, nelorb_c, type_lambda, nnozero_c, nozero_c, jbradet, ipc, max_iter &
                              , optimize, detmatc_in, overo, oversn, oversn(1, nelorb_c + 1), overon, overon(1, nelorb_c + 1) &
                              , mat, prec, overlapsquare, rank, 1, 1)
            else
                call max_ovlp(nelorbc_in, nelorb_c, type_lambda, nnozero_c, nozero_c, jbradet, ipc, max_iter &
                              , optimize, detmatc_in, overo, oversn, oversn, overon, overon &
                              , mat, prec, overlapsquare, rank, 1, 1)
            end if
        else
            if (ipc .eq. 2) then
                call max_ovlp(nelorbc_in, nelorb_c, type_lambda, nnozero_c, nozero_c, jbradet, ipc, max_iter &
                              , optimize, detmatc_in, overo, oversn, oversn(1, nelorb_c + 1), overon, overon(1, nelorb_c + 1) &
                              , mat, prec, overlapsquare, rank, nprocu, mpi_comm_world)
            else
                call max_ovlp(nelorbc_in, nelorb_c, type_lambda, nnozero_c, nozero_c, jbradet, ipc, max_iter &
                              , optimize, detmatc_in, overo, oversn, oversn, overon, overon &
                              , mat, prec, overlapsquare, rank, nprocu, mpi_comm_world)
            end if
        end if

        deallocate (optimize)

        if (rank .eq. 0) write (6, *) ' Overlap square =', overlapsquare

    else

        if (overlap) then
            if (allocated(mat)) deallocate (mat)
            allocate (mat(nelorb_cu, 2*nelorb_cu))
            mat = 0.d0
            if (npar_eagp .gt. 0) then
                call copy_eagp(.true., ipc, nelorb_c, nelcol_c, detmat_c, eagp_pfaff, mat)
                deallocate (detmat_c)
                allocate (detmat_c(ipc*nelcol_c*nelcol_c))
                do k = 1, nelorb_cu
                    detmat_c(nelorb_cu*(k - 1) + 1:nelorb_cu*k) = mat(1:nelorb_cu, k)
                end do
            end if

            if (contraction .ne. 0) then
                ! compute S_out*lambda_out
                call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0           &
                        &, oversn, nelorb_cu, detmat_c, nelorb_cu, 0.d0, mat, nelorb_cu, nprocu, rank, comm_mpi)
                ! compute S_out*lambda_out*S_out
                call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0           &
                        &, mat, nelorb_cu, oversn, nelorb_cu, 0.d0, mat(1, nelorb_cu + 1), nelorb_cu, nprocu, rank, comm_mpi)
                ! compute S_out*lambda_out*S_out*lambda_out
                call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0           &
                        &, mat(1, nelorb_cu + 1), nelorb_cu, detmat_c, nelorb_cu, 0.d0, mat, nelorb_cu, nprocu, rank, comm_mpi)
            else
                ! compute S_out*lambda_out
                call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0           &
                        &, oversn, nelorb_cu, detmat, nelorb_cu, 0.d0, mat, nelorb_cu, nprocu, rank, comm_mpi)
                ! compute S_out*lambda_out*S_out
                call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0           &
                        &, mat, nelorb_cu, oversn, nelorb_cu, 0.d0, mat(1, nelorb_cu + 1), nelorb_cu, nprocu, rank, comm_mpi)
                ! compute S_out*lambda_out*S_out*lambda_out
                call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0           &
                        &, mat(1, nelorb_cu + 1), nelorb_cu, detmat, nelorb_cu, 0.d0, mat, nelorb_cu, nprocu, rank, comm_mpi)
            end if

            ! compute Tr [S_out*lambda_out*S_out*lambda_out]
            overon_norm = tracematc(nelorb_cu, mat, nelorb_cu)

        end if

        !         evaluation matrix M
        if (allocated(mat)) deallocate (mat)
        allocate (mat(ipc*max(nelorb_cu, nelorbcu_in), nelorb_cu))
        dim_mat = max(nelorb_cu, nelorbcu_in)
        mat = 0.d0

        ! compute lambda_in*O
        if (ipc .eq. 1) then
            call dgemm_my('N', 'N', nelorbcu_in, nelorb_cu, nelorbcu_in, 1.d0           &
                    &, detmatc_in, nelorbcu_in, overon, nelorbcu_in, 0.d0, mat, dim_mat, nprocu, rank, comm_mpi)
        end if

        !       calculation inverse of overlap matrix
        lwork = nelorb_cu*nelorb_cu

        if (npar_eagp .gt. 0) then
            imin = 1
            imax = nelorb_cu
            dimo = nelorb_cu
        else
            call findminmax(nnozero_c, nozero_c, nelorb_c, imin, imax, dimo)
            if (molecular .gt. 0) then
                imin = imax - molecular + 1
                dimo = molecular - nelup + neldo
                if (rank .eq. 0) write (6, *) ' Changed imin =', imin, imax, dimo
            end if
        end if
        over_sav(1:ipc*imax, 1:imax) = oversn(1:ipc*imax, 1:imax)
        call invsymeps(ipc, dimo, oversn(ipc*(imin - 1) + 1, imin)&
                &, nelorb_cu, info, epsdgel, mine, umat, eigmat, nprocu, rank, comm_mpi)

        if (ipc .eq. 2) then
            over_sav(1:2*imax, nelorb_cu + 1:nelorb_cu + imax) = oversn(1:2*imax, nelorb_cu + 1:nelorb_cu + imax)
            call invsymeps(ipc, dimo, oversn(2*imin - 1, nelorb_cu + imin)&
                    &, nelorb_cu, info, epsdgel, mine, umat, eigmat, nprocu, rank, comm_mpi)
        else
            call dgemm_my('N', 'N', nelorbcu_in, dimo, dimo, 1.d0, overon(1, imin)&
                    &, nelorbcu_in, oversn(imin, imin), nelorb_cu, 0.d0, psip, nelorbcu_in, nprocu, rank, comm_mpi)
        end if

        !
        if (ipc .eq. 1) then
            inv_sav = 0.d0
            inv_sav(imin:imin + dimo - 1, imin:imin + dimo - 1) = &
                    &oversn(imin:imin + dimo - 1, imin:imin + dimo - 1)

            call dgemm_my('T', 'N', dimo, dimo, nelorbcu_in, 1.d0, psip           &
                    &, nelorbcu_in, mat(1, imin), dim_mat, 0.d0, oversn, dimo, nprocu, rank, comm_mpi)
        end if

        if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
            overlapsquare = tracemat(dimo, oversn, dimo)/overo
        elseif (ipc .eq. 1) then

            call dgemm_my('T', 'N', nelorbcu_in, nelorb_cu, nelorbcu_in, 1.d0           &
                    &, detmatc_in, nelorbcu_in, overon, nelorbcu_in, 0.d0, mat, dim_mat, nprocu, rank, comm_mpi)

            call dgemm_my('T', 'N', dimo, dimo, nelorbcu_in, 1.d0, psip           &
                    &, nelorbcu_in, mat(1, imin), dim_mat, 0.d0, psip(nelorbcu_in*dimo + 1), dimo, nprocu, rank, comm_mpi)

            overlapsquare = tracemat2(dimo, dimo, oversn, dimo, psip(nelorbcu_in*dimo + 1), dimo)/overo

        else

            !  dimo is the new basis restricted to the non zero possible values according to findiminimax
            !  nelorbc_in is the old basis
            !  New algorithm using mat and psip as auxiliary matrices

            ! \bar S_L (S^\prime_L)^-1
            call zgemm_my('N', 'N', nelorbcu_in, dimo, dimo, zone, overon(1, imin)&
                    &, nelorbcu_in, oversn(2*imin - 1, imin), nelorb_cu, zzero, psip, nelorbcu_in, nprocu, rank, comm_mpi)
            !
            !   lambda \bar S_R^*
            call conjmat(nelorbcu_in, dimo, overon(1, imin + nelorb_cu), nelorbcu_in)
            call zgemm_my('N', 'N', nelorbcu_in, dimo, nelorbcu_in, zone, detmatc_in, nelorbcu_in&
                    &, overon(1, imin + nelorb_cu), nelorbcu_in, zzero, mat, nelorbcu_in, nprocu, rank, comm_mpi)
            call conjmat(nelorbcu_in, dimo, overon(1, imin + nelorb_cu), nelorbcu_in)
            !  umat=A = psip^dag mat= (S^\prime_L)^-1 \bar S_L^dag lambda \bar S_R^*
            call zgemm_my('C', 'N', dimo, dimo, nelorbcu_in, zone, psip&
                    &, nelorbcu_in, mat, nelorbcu_in, zzero, umat, dimo, nprocu, rank, comm_mpi)
            !! S_R \lambda^T
            !    call zgemm_my('N','T',dimo,dimo,dimo,zone,overs(2*nelorb_at+1,nelorb_diag+nelorb_at+1),nelorb_diag&
            !             &,detmat_c(inddetc),nelorb_diag,zzero,mat,dimo,nprocu,rankopt,commopt_mpi)
            !            compute inverse right. Here we cannot assume R and L are the same.
            !   Definition of mat=\bar S_R^T \lambda^dag
            call zgemm_my('T', 'C', dimo, nelorbcu_in, nelorbcu_in, zone&
                    &, overon(1, nelorb_cu + imin), nelorbcu_in&
                    &, detmatc_in, nelorbcu_in, zzero, mat, dimo, nprocu, rank, comm_mpi)
            !    B= psip=S_R^*-1 \bar S_R^T \lambda^\dag  \bar S_L
            call conjmat(dimo, dimo, oversn(2*imin - 1, nelorb_cu + imin), nelorb_cu)
            call zgemm_my('N', 'N', dimo, nelorbcu_in, dimo, zone, oversn(2*imin - 1, nelorb_cu + imin), nelorb_cu&
                    &, mat, dimo, zzero, psip, dimo, nprocu, rank, comm_mpi)
            call conjmat(dimo, dimo, oversn(2*imin - 1, nelorb_cu + imin), nelorb_cu)
            !  Use mat
            call zgemm_my('N', 'N', dimo, dimo, nelorbcu_in, zone, psip&
                    &, dimo, overon(1, imin), nelorbcu_in, zzero, mat, dimo, nprocu&
                    &, rankopt, commopt_mpi)
            !
            overlapsquare = tracemat2(dimo, dimo, umat, dimo, mat, dimo)/overo
        end if

        if (rank .eq. 0 .and. overo .ne. 0.d0) write (6, *)&
                &  ' Overlap square Geminal  found =', overlapsquare
        !       normalization
        !OK

        !       to have the same normalization of the previous one
        !       cost=1.d0/dsqrt(overlapsquare)
        !       to be closest in L2 norm
        cost = 1.d0
        deallocate (mat)
        allocate (mat(ipc*nelorb_cu, nelorb_cu))
        mat = 0.d0 ! vanish all matrix elements untouched

        if (overo .ne. 0.d0) then

            ! compute lambda_proj (lambda_in projected onto the basis set out)
            if (ipc .eq. 1) then
                call dgemm_my('N', 'N', dimo, dimo, dimo, cost, oversn&
                        &, dimo, inv_sav(imin, imin), nelorb_cu, 0.d0, mat(imin, imin), nelorb_cu, nprocu, rank, comm_mpi)
                call checkmat(nelorb_c, mat, nelorb_cu, nozero_c, nnozero_c, checkall, rank, 0.d0, symmagp, ipf)
            else
                call zgemm_my('N', 'T', dimo, dimo, dimo, zone, umat&
                        &, dimo, oversn(2*imin - 1, nelorb_cu + imin), nelorb_cu, zzero, mat(2*imin - 1, imin), nelorb_cu&
                        &, nprocu, rank, comm_mpi)
                call checkmat_complex(nelorb_c, mat, nelorb_cu, nozero_c, nnozero_c, checkall, rank, 0.d0, symmagp, ipf)
            end if

            if (.not. checkall .or. overlap) then

                ! compute S_out*lambda_proj
                if (ipc .eq. 1 .and. symmagp .and. ipf .eq. 1) then
                    call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, over_sav&
                            &, nelorb_cu, mat, nelorb_cu, 0.d0, psip, nelorb_cu, nprocu, rank, comm_mpi)
                    over_new = tracemat(nelorb_cu, psip, nelorb_cu)/overo
                elseif (ipc .eq. 1) then
                    call dgemm_my('T', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, over_sav, nelorb_cu&
                            &, mat, nelorb_cu, 0.d0, psip, nelorb_cu, nprocu, rank, comm_mpi)
                    call dgemm_my('N', 'T', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, over_sav, nelorb_cu&
                            &, mat, nelorb_cu, 0.d0, psip(nelorb_cu*nelorb_cu + 1), nelorb_cu, nprocu, rank, comm_mpi)
                    over_new = tracemat2(nelorb_cu, nelorb_cu, psip, nelorb_cu, psip(nelorb_cu*nelorb_cu + 1), nelorb_cu)/overo
                else
                    call conjmat(nelorb_cu, nelorb_cu, mat, nelorb_cu)
                    !  S_L^T \lambda^*
                    call zgemm_my('T', 'N', nelorb_cu, nelorb_cu, nelorb_cu, zone, over_sav, nelorb_cu&
                            &, mat, nelorb_cu, zzero, psip, nelorb_cu, nprocu, rank, comm_mpi)
                    call conjmat(nelorb_cu, nelorb_cu, mat, nelorb_cu)
                    ! S_R \lambda^T ! destroy the left overlape no longer used
                    call zgemm_my('N', 'T', nelorb_cu, nelorb_cu, nelorb_cu, zone, over_sav(1, nelorb_cu + 1), nelorb_cu&
                            &, mat, nelorb_cu, zzero, over_sav, nelorb_cu, nprocu, rank, comm_mpi)

                    over_new = tracemat2(nelorb_cu, nelorb_cu, psip, nelorb_cu, over_sav, nelorb_cu)/overo

                end if

                if (.not. checkall .and. rank .eq. 0) then
                    write (6, *) ' Warning the matrix  AGP has some extra non zero element '
                    write (6, *) 'Overlap square with no zero', over_new
                end if

                if (overlap) then
                    ! compute S_out*lambda_proj*S_out
                    call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, psip &
                                  , nelorb_cu, over_sav, nelorb_cu, 0.d0, psip(nelorb_cu*nelorb_cu + 1) &
                                  , nelorb_cu, nprocu, rank, comm_mpi)

                    ! compute S_out*lambda_proj*S_out*lambda_out
                    if (contraction .ne. 0) then
                        call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, psip(nelorb_cu*nelorb_cu + 1)&
                                &, nelorb_cu, detmat_c, nelorb_cu, 0.d0, psip, nelorb_cu, nprocu, rank, comm_mpi)
                    else
                        call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, psip(nelorb_cu*nelorb_cu + 1)&
                                &, nelorb_cu, detmat, nelorb_cu, 0.d0, psip, nelorb_cu, nprocu, rank, comm_mpi)
                    end if

                    ! compute Tr [S_out*lambda_proj*S_out*lambda_out]
                    overodot_norm = tracemat(nelorb_cu, psip, nelorb_cu)

                    ! compute S_out*lambda_proj*S_out*lambda_proj
                    call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, psip(nelorb_cu*nelorb_cu + 1)&
                            &, nelorb_cu, mat, nelorb_cu, 0.d0, psip, nelorb_cu, nprocu, rank, comm_mpi)
                    ! compute Tr [S_out*lambda_proj*S_out*lambda_proj]
                    overo_norm_inout = tracemat(nelorb_cu, psip, nelorb_cu)

                    if (rank .eq. 0) then
                        write (6, *) 'Overlap between _in and _out geminals'
                        write (6, *) 'absolute normalization', &
                                &overodot_norm/sqrt(overo_norm*overon_norm)
                        write (6, *) 'normalization in _out basis set', &
                                &overodot_norm/sqrt(overo_norm_inout*overon_norm)
                    end if
                end if

            end if

        end if

    end if ! end if optmat

    if (contraction .ne. 0) then
        if (npar_eagp .gt. 0) then
            deallocate (detmat_c)
            allocate (detmat_c(ipc*nelorb_cu*nelorb_cu))
        end if
        call dcopy(ipc*nelorb_cu*nelorb_cu, mat, 1, detmat_c, 1)
    else
        call dcopy(ipc*nelorb_cu*nelorb_cu, mat, 1, detmat, 1)
    end if

    !      Now the unpaired orbitals....
    ! number of unpaired orbitals
    !       ndiff=nelcol_c-nelorb_c
    if (nelcol_c .gt. nelorb_c .and. npar_eagp .eq. 0) then
        if (molecular .gt. 0) dimo = dimo + nelup - neldo
        if (ndiff .ne. nelcolc_in - nelorbc_in) then
            if (rank .eq. 0) write (6, *) ' The two wf should  have the same&
                    &number of unpaired orb.  !!!'
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if
        !       now the new matrix  <O| N unpaired > --> mat
        !       restoring the right dimension if molecular>0

        if (ipc .eq. 1) then
            call dgemm_my('T', 'N', nelorb_c, ndiff, nelorbc_in, 1.d0               &
                    &, overon, nelorbc_in, detmatc_in(nelorbc_in*nelorbc_in + 1)            &
                    &, nelorbc_in, 0.d0, mat, nelorb_c, nprocu, rank, comm_mpi)
        else
            call conjmat(nelorbc_in, ndiff, detmatc_in(2*nelorbc_in*nelorbc_in + 1), nelorbc_in)
            call zgemm_my('C', 'N', nelorb_c, ndiff, nelorbc_in, zone               &
                    &, overon, nelorbc_in, detmatc_in(2*nelorbc_in*nelorbc_in + 1)            &
                    &, nelorbc_in, zzero, mat, nelorb_c, nprocu, rank, comm_mpi)
            call conjmat(nelorbc_in, ndiff, detmatc_in(2*nelorbc_in*nelorbc_in + 1)&
                    &, nelorbc_in)
        end if

        if (contraction .eq. 0) then
            if (allocated(gammamat)) deallocate (gammamat)
            allocate (gammamat(ipc*ndiff, ndiff))
            gammamat = 0.d0
            if (ipc .eq. 1) then
                !       calculation overlap matrix Gamma
                call dgemm_my('N', 'N', nelorbc_in, ndiff, nelorbc_in, 1.d0, overs       &
                        &, nelorbc_in, detmatc_in(nelorbc_in*nelorbc_in + 1), nelorbc_in, 0.d0   &
                        &, psip, nelorbc_in, nprocu, rank, comm_mpi)

                call dgemm_my('T', 'N', ndiff, ndiff, nelorbc_in, 1.d0                  &
                        &, detmatc_in(nelorbc_in*nelorbc_in + 1), nelorbc_in, psip, nelorbc_in, 0.d0&
                        &, gammamat, ndiff, nprocu, rank, comm_mpi)

            else

                !       calculation overlap matrix Gamma
                call zgemm_my('N', 'N', nelorbc_in, ndiff, nelorbc_in, zone, overs       &
                        &, nelorbc_in, detmatc_in(2*nelorbc_in*nelorbc_in + 1), nelorbc_in, zzero   &
                        &, psip, nelorbc_in, nprocu, rank, comm_mpi)

                call zgemm_my('C', 'N', ndiff, ndiff, nelorbc_in, zone                  &
                        &, detmatc_in(2*nelorbc_in*nelorbc_in + 1), nelorbc_in, psip, nelorbc_in, zzero&
                        &, gammamat, ndiff, nprocu, rank, comm_mpi)

            end if

            !       write(6,*) ' Overlap matrix gamma '
            !       do i=1,ndiff
            !         do j=1,ndiff
            !         write(6,*) i,j,gammamat(i,j)
            !         enddo
            !       enddo

            !       Now the inverse of this matrix
            call invsym(ndiff, gammamat, ndiff, info)
        end if

        !      now building the new matrix --> oversn

        if (ipc .eq. 1) then

            call dgemm_my('N', 'T', ndiff, nelorb_c, ndiff, 1.d0, gammamat, ndiff      &
                    &, mat, nelorb_c, 0.d0, psip, ndiff, nprocu, rank, comm_mpi)

            call dgemm_my('N', 'N', nelorb_c, nelorb_c, ndiff, 1.d0, mat, nelorb_c     &
                    &, psip, ndiff, 0.d0, oversn, nelorb_c, nprocu, rank, comm_mpi)

        else

            call zgemm_my('N', 'C', ndiff, nelorb_c, ndiff, zone, gammamat, ndiff      &
                    &, mat, nelorb_c, zzero, psip, ndiff, nprocu, rank, comm_mpi)

            call zgemm_my('N', 'N', nelorb_c, nelorb_c, ndiff, zone, mat, nelorb_c     &
                    &, psip, ndiff, zzero, oversn, nelorb_c, nprocu, rank, comm_mpi)

        end if

        deallocate (psip)
        lwork = ipc*nelorb_c*nelorb_c + nelorb_c
        allocate (psip(lwork))
        psip = 0.d0

        if (rank .eq. 0) then
            write (6, *) ' input martices ', nelorb_c
            write (6, *) ' minimum index  ', imin
            write (6, *) ' dimension uncontr =  ', dimo
        end if

        if (molecular .gt. 0 .or. ipc .eq. 2) then
            call invsymeps(ipc, dimo, over_sav(imin, imin)&
                    &, nelorb_c, info, epsdgel, mine, umat, eigmat, nprocu, rank, comm_mpi)
        end if

        call dsygv_my(ipc, dimo, oversn(ipc*(imin - 1) + 1, imin), nelorb_c, over_sav(ipc*(imin - 1) + 1, imin)&
                &, nelorb_c, .false., umat, eigmat, psip, mine, psip(nelorb_c + 1), epsdgel, .true., info, nprocu&
                &, rank, comm_mpi)

        if (imin .gt. 1) then
            oversn(1:ipc*(imin - 1), imin:imin + dimo - 1) = 0.d0
        end if

        overunpairedc = 1.d0
        do i = 1, ndiff
            overunpairedc = overunpairedc*psip(dimo - i + 1)
        end do

        if (rank .eq. 0) write (6, *) ' Overlap square unpaired contr. ='    &
                &, overunpairedc
        !        if(rank.eq.0) then
        !        write(6,*) ' Eigenvalues found ',ndiff
        !        do i=1,dimo
        !        write(6,*) i,psip(i)
        !        enddo
        !        endif

        !       stop

        !       write(6,*) ' Eigenvectors written ',imin
        !       do i=imin,imin+dimo-1
        !        write(6,*) ' eigenvector # ',i
        !        do j=imin,imin+dimo-1
        !        write(6,*) j,oversn(j,i)
        !        enddo
        !       enddo

        !       Here we have to check if orbitals are allowed

        if (allocated(occorb)) deallocate (occorb)
        allocate (occorb(ndiff))

        occorb = .true.
        do i = 1, ndiff
            !       among the last dimo orbitals take the one with the maximum
            !       overlap with iyt (hoping that has the some symmetry)
            iyt = nelorb_c + i
            overmax = 0.d0
            i_max = 0
            do ii = 1, ndiff
                if (occorb(ii)) then
                    cost = 0.d0
                    do j = 1, nnozero_c
                        iy = (nozero_c(j) - 1)/nelorb_c + 1
                        ix = nozero_c(j) - (iy - 1)*nelorb_c
                        if (iy .eq. iyt) cost = cost + sum(oversn(ipc*(ix - 1) + 1:ipc*ix, imin - 1 + dimo - ndiff + ii)**2)
                    end do
                    !            normalization in l2
                    cost = cost/sum(oversn(:, imin - 1 + dimo - ndiff + ii)**2)
                    if (cost .gt. overmax) then
                        overmax = cost
                        i_max = ii
                    end if
                end if
            end do
            occorb(i_max) = .false.
            if (i_max .eq. 0) then
                if (rank .eq. 0) write (6, *) ' There should be some error          &
                        &   in unpaired orbitals !!'
#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if

            if (contraction .ne. 0) then
                call dcopy(ipc*nelorb_c, oversn(1, imin - 1 + dimo - ndiff + i_max), 1       &
                        &, detmat_c(ipc*nelorb_c*(iyt - 1) + 1), 1)
            else
                call dcopy(ipc*nelorb_c, oversn(1, dimo + imin - 1 - ndiff + i_max), 1       &
                        &, detmat(ipc*nelorb_c*(iyt - 1) + 1), 1)
            end if
        end do

        deallocate (gammamat, occorb)

        ! endif unpaired orbitals
    end if

    if (allocated(overs)) deallocate (overs)
    if (allocated(oversn)) deallocate (oversn)
    if (allocated(overon)) deallocate (overon)

    !         end all about determinant
    if (nelorbjc_in .gt. 0 .and. nelorbj .gt. 0 .and. change_jas) then

        !     all about jastrow
        if (allocated(buffer)) deallocate (buffer)
        if (allocated(weightbuf)) deallocate (weightbuf)
        if (allocated(buffer_c)) deallocate (buffer_c)
        if (allocated(buffero)) deallocate (buffero)
        if (allocated(over_sav)) deallocate (over_sav)
        if (allocated(inv_sav)) deallocate (inv_sav)

        if (mesh .gt. nbufd) then
            allocate (buffer(nelorbjh, nbufd), buffero(nelorbjc_in, nbufd)       &
                    &, weightbuf(nbufd))
            if (contractionj .ne. 0) allocate (buffer_c(nelorbj_c, nbufd))
        else
            allocate (buffer(nelorbjh, mesh), buffero(nelorbjc_in, mesh)         &
                    &, weightbuf(mesh))
            if (contractionj .ne. 0) allocate (buffer_c(nelorbj_c, mesh))
        end if
        buffer = 0.d0
        buffero = 0.d0
        weightbuf = 0.d0
        if (contractionj .ne. 0) buffer_c = 0.d0

        allocate (oversnj(nelorbj_c, nelorbj_c), inv_sav(nelorbj_c, nelorbj_c)&
                &, overonj(nelorbjc_in, nelorbj_c), over_sav(nelorbj_c, nelorbj_c))
        oversnj = 0.d0
        inv_sav = 0.d0
        overonj = 0.d0
        over_sav = 0.d0

        if (allocated(umat)) deallocate (umat, eigmat)
        allocate (umat(nelorbj_c, nelorbj_c), eigmat(nelorbj_c))
        umat = 0.d0
        eigmat = 0.d0

        if (allocated(psip)) deallocate (psip)
        allocate (psip(max(nelorbj_c, nelorbjc_in)* &
                &    max(nelorbj_c, nelorbjc_in)))
        psip = 0.d0

        if (.not. bigram) then
#ifdef PARALLEL
            rewind (101)
#else
            rewind (12)
            ! read first useless record
            read (12)
#endif
        end if
        if (contractionj .ne. 0) then
            !      allocate(muj_sav(nelorbj_c))
            do i = 1, nelorbj_c
                if (orbcostn(i) .and. i .ne. adrcostc) then
                    !      muj_c(adrcost,i)=0.d0
                    muj_c(:, i) = 0.d0
                end if
            end do
        end if

        oversnj = 0.d0
        overonj = 0.d0

        iflagnorm = 3
        ind = 0

        !       indr=0
        indtot = 0
        indproc = 0
        nbuf = 0
#ifdef _OFFLOAD
!$omp target data map(oversnj,overonj) map(to:muj_c) map(alloc:buffer,buffero,buffer_c)
#endif
        do k = 1, nz
            !       x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
                !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
                    indtot = indtot + 1
                    !           indr=indr+1
                    !           if(indr-(indr/proc8)*proc8.eq.rank) then
                    if (indproc .eq. rankn) then
                        ind = ind + 1

                        x(:) = (-(nz + 1)/2.d0 + k)*az*at(:, 3) + (-(ny + 1)/2.d0 + j)*ay*at(:, 2) &
                               + (-(nx + 1)/2.d0 + i)*ax*at(:, 1) + rion_ref(:)
                        !           x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)

                        call upnewwf(1, 0, 0, 1, nshelljr, ioptorbj, ioccj, x, 1, r, rmu       &
                                &, vjur, zetar, rion, distp, buffer(1, ind), nelorbjh, nion, kionj          &
                                &, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                                &, indparj_tab, indorbj_tab, indshellj_tab, .true.)

                        call upvpot_ei(x, zetar, rion, weightbuf(ind), nion        &
                                &, LBox, epsvpot)

                        call dscal(nelorbjh, weightbuf(ind), buffer(1, ind), 1)

                    end if

                    if (indtot .eq. nbufrep .or. (i .eq. nx .and. j .eq. ny .and. k .eq. nz)) then
#ifdef _OFFLOAD
!$omp target update to (buffer)
#endif
                        if (contractionj .eq. 0) then
                            call dgemm_('N', 'T', nelorbj_c, nelorbj_c, ind, volmesh, buffer&
                                    &, nelorbjh, buffer, nelorbjh, 1.d0, oversnj, nelorbj_c)
                        else
                            call dgemm_('T', 'N', nelorbj_c, ind, nelorbjh, 1.d0, muj_c     &
                                    &, nelorbjh, buffer, nelorbjh, 0.d0, buffer_c, nelorbj_c)
                            call dgemm_('N', 'T', nelorbj_c, nelorbj_c, ind, volmesh       &
                                    &, buffer_c, nelorbj_c, buffer_c, nelorbj_c, 1.d0, oversnj, nelorbj_c)
                        end if

                        nbuf = nbuf + 1
                        if (bigram) then
                            buffero(:, :) = buffj(:, :, nbuf)
                        else
#ifdef PARALLEL
                            read (101) buffero
                            flush (101)
#else
                            read (12) buffero
#endif
                        end if
                        do ii = 1, ind
                            call dscal(nelorbjc_in, weightbuf(ii), buffero(1, ii), 1)
                        end do
#ifdef _OFFLOAD
!$omp target update to (buffero)
#endif
                        if (contractionj .eq. 0) then
                            call dgemm_('N', 'T', nelorbjc_in, nelorbj_c, ind, volmesh     &
                                    &, buffero, nelorbjc_in, buffer, nelorbjh, 1.d0, overonj, nelorbjc_in)
                        else
                            call dgemm_('N', 'T', nelorbjc_in, nelorbj_c, ind, volmesh     &
                                    &, buffero, nelorbjc_in, buffer_c, nelorbj_c, 1.d0, overonj, nelorbjc_in)
                        end if
                        ind = 0
                        indtot = 0
                    end if ! endif load buf
                    indproc = indproc + 1
                    if (indproc .eq. nprocn) indproc = 0
                end do
            end do
        end do
#ifdef _OFFLOAD
!$omp end target data
#endif
#ifdef PARALLEL
!      collect overonj,oversnj
        nelorb2 = nelorbj_c*nelorbj_c
        call reduce_base_real(nelorb2, oversnj, MPI_COMM_WORLD, -1)
        nelorb2 = nelorbjc_in*nelorbj_c
        call reduce_base_real(nelorb2, overonj, MPI_COMM_WORLD, -1)
#endif

        over_sav = oversnj

        !         evaluation matrix M
        if (allocated(mat)) deallocate (mat)
        dim_mat = max(nelorbjc_in, nelorbj_c)
        allocate (mat(dim_mat, nelorbj_c))
        mat = 0.d0

        call dgemm_my('N', 'N', nelorbjc_in, nelorbj_c, nelorbjc_in, 1.d0        &
                &, jasmatc_in, nelorbjc_in, overonj, nelorbjc_in, 0.d0, mat, dim_mat, nprocu, rank, comm_mpi)

        call findminmax(nnozeroj_c, nozeroj_c, nelorbj_c, imin, imax, dimo)

        !       if(rank.eq.0)                                                   &
        !    &  write(6,*) ' Minimum and Maximum orbital allowed in new Jas ='  &
        !    &,imin,imax,dimo,nelorbj_c

        if (molecularj .gt. 0) then
            imin = imax - molecularj + 1
            dimo = molecularj
            if (iessz) then
                imax = imax - molecularj/2
                dimo = molecularj/2
            end if
            if (rank .eq. 0) write (6, *) ' Changed imin =', imin, imax, dimo
        end if

        !       calculation inverse of overlap matrix
        !       lwork=nelorbj_c*nelorbj_c
        !       call dgetrf(dimo,dimo,oversnj(imin,imin),nelorbj_c,ipsip,info)
        !       call dgetri(dimo,oversnj(imin,imin),nelorbj_c,ipsip,psip
        !    1,lwork,info)

        !       write(6,*) ' before invsymeps ',dimo,imin

        !       write(6,*) ' size oversn =',size(oversnj)
        !       write(6,*) ' size over_sav =',size(over_sav)
        !       write(6,*) ' size umat =',size(umat)
        !       write(6,*) ' size eigmat =',size(eigmat)
        !      write(6,*) ' input oversnj  overonj ',imin,imax,dimo,nelorbj_c

        !      do i=1,nelorbjc_in
        !        do j=imin,imax
        !        write(6,*) i,j,overonj(i,j),oversnj(i,j)
        !        enddo
        !      enddo

        over_sav(1:imax, 1:imax) = oversnj(1:imax, 1:imax)
        call invsymeps(1, dimo, oversnj(imin, imin)&
                &, nelorbj_c, info, epsdgel, mine, umat, eigmat, nprocu, rank, comm_mpi)

        call dgemm_my('N', 'N', nelorbjc_in, dimo, dimo, 1.d0, overonj(1, imin)&
                &, nelorbjc_in, oversnj(imin, imin), nelorbj_c, 0.d0, psip, nelorbjc_in, nprocu, rank, comm_mpi)

        inv_sav = oversnj

        call dgemm_my('T', 'N', dimo, dimo, nelorbjc_in, 1.d0, psip           &
                &, nelorbjc_in, mat(1, imin), dim_mat, 0.d0, oversnj, dimo, nprocu, rank, comm_mpi)

        overlapsquarej = tracemat(dimo, oversnj, dimo)/overoj

        if (rank .eq. 0 .and. overoj .ne. 0.d0)                                                   &
                &  write (6, *) ' Overlap square Jas  found =', overlapsquarej
        !       cost=1.d0/dsqrt(overlapsquarej)
        !       to be closest in L2 norm
        cost = 1.d0
        deallocate (mat)
        allocate (mat(nelorbj_c, nelorbj_c))
        mat = 0.d0

        call dgemm_my('N', 'N', dimo, dimo, dimo, cost, oversnj&
                &, dimo, inv_sav(imin, imin), nelorbj_c, 0.d0, mat(imin, imin), nelorbj_c, nprocu, rank, comm_mpi)

        call setorbcost(nelorbj_c, mat, orbcostn, nozeroj_c, nnozeroj_c)
        call purify(mat, nelorbj_c, over_sav, nelorbj_c, orbcostn)

        call checkmat(nelorbj_c, mat, nelorbj_c, nozeroj_c, nnozeroj_c, checkall, rank, 0.d0, .true., 1)
        if (.not. checkall) then
            if (rank .eq. 0) write (6, *) ' Warning the matrix Jas has some extra&
                    & non zero element '
            call dgemm_my('N', 'N', nelorbj_c, nelorbj_c, nelorbj_c, 1.d0, over_sav&
                    &, nelorbj_c, mat, nelorbj_c, 0.d0, psip, nelorbj_c, nprocu, rank, comm_mpi)
            if (rank .eq. 0)&
                    &write (6, *) 'Overlap square with no zero', tracemat(nelorbj_c, psip, nelorbj_c)/overoj
        end if

        if (contractionj .ne. 0) then
            call dcopy(nelorbj_c*nelorbj_c, mat, 1, jasmat_c, 1)
        else
            call dcopy(nelorbj_c*nelorbj_c, mat, 1, jasmat, 1)
        end if
        if (iessz .and. iessz_in) then

            deallocate (mat)
            allocate (mat(dim_mat, nelorbj_c))
            mat = 0.d0

            call dgemm_my('N', 'N', nelorbjc_in, nelorbj_c, nelorbjc_in, 1.d0        &
                    &, jasmatszc_in, nelorbjc_in, overonj, nelorbjc_in, 0.d0                &
                    &, mat, dim_mat, nprocu, rank, comm_mpi)

            !       recover old def of psip
            call dgemm_my('N', 'N', nelorbjc_in, dimo, dimo, 1.d0, overonj(1, imin)&
                    &, nelorbjc_in, inv_sav(imin, imin), nelorbj_c, 0.d0, psip, nelorbjc_in, nprocu, rank, comm_mpi)

            !       calculation inverse of overlap matrix
            !       lwork=nelorbj_c*nelorbj_c
            !       call dgetrf(nelorbj_c,nelorbj_c,oversnj,nelorbj_c,ipsip,info)
            !       call dgetri(nelorbj_c,oversnj,nelorbj_c,ipsip,psip,lwork,info)

            if (molecularj .gt. 0) then
                imin = imax + 1
                imax = imax + molecularj/2
                dimo = molecularj/2
                if (rank .eq. 0) write (6, *) ' Changed imin Sz =', imin, imax, dimo
                over_sav(1:imax, 1:imax) = oversnj(1:imax, 1:imax)
                call invsymeps(1, dimo, oversnj(imin, imin)&
                        &, nelorbj_c, info, epsdgel, mine, umat, eigmat, nprocu, rank, comm_mpi)
                inv_sav = oversnj
                !       recover old def of psip
                call dgemm_my('N', 'N', nelorbjc_in, dimo, dimo, 1.d0, overonj(1, imin)&
                        &, nelorbjc_in, inv_sav(imin, imin), nelorbj_c, 0.d0, psip, nelorbjc_in, nprocu, rank, comm_mpi)
            end if

            call dgemm_my('T', 'N', dimo, dimo, nelorbjc_in, 1.d0, psip           &
                    &, nelorbjc_in, mat(1, imin), dim_mat, 0.d0, oversnj, dimo, nprocu, rank, comm_mpi)

            overlapsquarejsz = tracemat(dimo, oversnj, dimo)/overojsz

            if (rank .eq. 0 .and. overojsz .ne. 0.d0)                             &
                    &  write (6, *) ' Overlap square Jas Sz found =', overlapsquarejsz

            !       cost=1.d0/dsqrt(overlapsquarejsz)
            !       to be closest in L2 norm
            cost = 1.d0
            deallocate (mat)
            allocate (mat(nelorbj_c, nelorbj_c))
            mat = 0.d0
            call dgemm_my('N', 'N', dimo, dimo, dimo, cost, oversnj, dimo&
                    &, inv_sav(imin, imin), nelorbj_c, 0.d0, mat(imin, imin), nelorbj_c, nprocu, rank, comm_mpi)

            call setorbcost(nelorbj_c, mat, orbcostn, nozeroj_c, nnozeroj_c)
            call purify(mat, nelorbj_c, over_sav, nelorbj_c, orbcostn)

            call checkmat(nelorbj_c, mat, nelorbj_c, nozeroj_c, nnozeroj_c, checkall, rank, 0.d0, .true., 1)
            if (.not. checkall) then
                if (rank .eq. 0) write (6, *) ' Warning the matrix Jas-Sz  has some extra non zero element '
                call dgemm_my('N', 'N', nelorbj_c, nelorbj_c, nelorbj_c, 1.d0, over_sav&
                        &, nelorbj_c, mat, nelorbj_c, 0.d0, psip, nelorbj_c, nprocu, rank, comm_mpi)
                if (rank .eq. 0) &
                    write (6, *) 'Overlap square with no zero', tracemat(nelorbj_c, psip, nelorbj_c)/overojsz
            end if

            if (contractionj .ne. 0) then
                call dcopy(nelorbj_c*nelorbj_c, mat, 1, jasmatsz_c, 1)
            else
                call dcopy(nelorbj_c*nelorbj_c, mat, 1, jasmatsz, 1)
            end if

            !      write(6,*) ' input jasmat Sz '
            !      do ix=1,nelorbjc_in
            !         do iy=ix,nelorbjc_in
            !         write(6,*) ix,iy,jasmatszc_in(nelorbjc_in*(iy-1)+ix)
            !         enddo
            !      enddo
            !      write(6,*) ' output jasmatsz_c ',nelorbj_c
            !      do ix=1,nelorbj_c
            !         do iy=ix,nelorbj_c
            !      write(6,*) ix,iy,jasmatsz_c(nelorbj_c*(iy-1)+ix)
            !      enddo
            !      enddo

            !      stop

        else
            jasmatsz_c = 0.d0
            jasmatsz = 0.d0
        end if

    else

        if (contractionj .ne. 0) then
            jasmat_c = 0.d0
            jasmatsz_c = 0.d0
        else
            jasmat = 0.d0
            jasmatsz = 0.d0
        end if

    end if

    !         end all about determinant
    if (npar_eagp .gt. 0) then
        if (allocated(mat)) deallocate (mat)
        allocate (mat(ipc*nelcol_c, nelcol_c))
        if (contraction .ne. 0) then
            call dcopy(ipc*nelcol_c*nelcol_c, detmat_c, 1, mat, 1)
            call copy_eagp(.false., ipc, nelorb_c, nelcol_c, detmat_c, eagp_pfaff, mat)
        else
            call dcopy(ipc*nelcol_c*nelcol_c, detmat, 1, mat, 1)
            call copy_eagp(.false., ipc, nelorb_c, nelcol_c, detmat, eagp_pfaff, mat)
        end if
        deallocate (mat)
    end if

    if (rank .eq. 0) then
        close (10)
        open (unit=10, file='fort.10_new', form='formatted', status='unknown')
    end if
    if (change_jas .and. niesd .eq. niesd_in) then
        vj(1:niesd) = vj_in(1:niesd)
    end if
    if (rank .eq. 0) then
        !     allowed_averagek=.false.
        call write_fort10(10)
        close (10)
    end if

#ifdef PARALLEL
    if (rank .eq. 0) then
        close (11, status='DELETE')
        close (12, status='DELETE')
    end if
    if (.not. bigram) then
        close (100, status='DELETE')
        close (101, status='DELETE')
    end if
#else
    close (11, status='DELETE')
    close (12, status='DELETE')
#endif
    if (allocated(buff)) deallocate (buff)
    if (allocated(buffj)) deallocate (buffj)

#ifdef PARALLEL
    call mpi_finalize(ierr)
#endif
    stop

contains

    subroutine load_fort10in
        use allio, only: norm_metric
        implicit none
        real*8 :: mind(3), psilnn, r0, rc(3)
        real*8, external :: jastrow_ei, tracematnloc, tracematloc
        integer indmax

        !     begin load_fort10in

        if (rankn .eq. 0)                                                    &
                & open (unit=10, file='fort.10_in', form='formatted', status='unknown')

        rank = rankn
        nw = nprocn
        nproc = nprocn
        in1 = 1
        call default_allocate
        yesfast = 0 ! it works only with allocation of detmat
        if (rank .eq. 0) then
            call read_fort10_fast
            npsa = npsar
        end if
#ifdef PARALLEL
        call mpi_bcast(npsa, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        npsar = npsa
#endif
        pseudofile = 'pseudo.dat'
        rewind (8)
        call read_pseudo
        call read_fort10(10)

        nelorbpf = nelorbh*ipf
        nelorb_cu = nelorb_c
        if (npar_eagp .gt. 0) then
            nelorbpf = nelorbpf + ndiff
            nelorb_cu = nelorb_c + ndiff
        end if

        allocate (orbcostn_in(max(nelorbj_c, 1)))
        if (npar_eagp .gt. 0) then
            allocate (detmatc_in(ipc*nelcol_c*nelcol_c))
        else
            allocate (detmatc_in(ipc*nelorb_c*nelcol_c))
        end if
        allocate (jasmatc_in(max(nelorbj_c*nelorbj_c, 1)))
        jasmatc_in = 0.d0
        if (iessz) then
            allocate (jasmatszc_in(max(nelorbj_c*nelorbj_c, 1)))
            jasmatszc_in = 0.d0
        end if
        dimdistp = ipc*(indt + 5)*nelorbh + nelorbh*(indt + 6) + 27*(indt + 1)*nshell + nelorbh
        dimdistp = max(dimdistp, (indt + 5)*nelorbjh + nelorbjh*(indt + 6) + 27*(indt + 1)*nshellj + nelorbjh)
        allocate (distp(dimdistp))
        distp = 0.d0

        contractionj_in = contractionj

        ipf_in = ipf
        ipc_in = ipc
        contraction_in = contraction
        if (contractionj .gt. 0) then
            allocate (mujc_in(nelorbjh, nelorbj_c))
            mujc_in = muj_c
            jasmatc_in = jasmat_c
            if (iessz) jasmatszc_in = jasmatsz_c
        else
            jasmatc_in = jasmat
            if (iessz) jasmatszc_in = jasmatsz
        end if
        orbcostn_in = orbcostn
        adrcost = 0
        do i = 1, nelorbjh
            if (orbcostl(i)) adrcost = i
        end do
        if (contractionj .ne. 0) then
            adrcostc = 0
            do i = 1, nelorbj_c
                if (orbcostn(i)) adrcostc = i
            end do
        else
            adrcostc = adrcost
        end if
        adrcost_in = adrcost
        adrcostc_in = adrcostc

        iessz_in = iessz
        niesd_in = niesd
        iesdrr_in = iesdrr
        allocate (vj_in(niesd_in))
        vj_in(1:niesd) = vj(1:niesd)

        nelcolc_in = nelcol_c
        nelorbc_in = nelorb_c
        ndiff_in = ndiff
        nelorbcu_in = nelorb_c
        if (npar_eagp .gt. 0) nelorbcu_in = nelorb_c + ndiff
        nelorb_cu = nelorbcu_in
        nelorbjc_in = nelorbj_c
        if (ipj .eq. 2 .and. change_jas) then
            if (rank .eq. 0) write (6, *) ' Warning change_jas is set to .false., as it does not &
                    &work for generic Jastrow  -12/-22'
            change_jas = .false.
        end if
        if (rank .eq. 0) then
            !      write(6,*) ' Input buffer dim '
            !      read(5,*) nbufd
            nx = 0
            ny = 0
            nz = 0
            ax = 0.
            ay = 0.
            az = 0
            nbufd = -1
            !       Default values
            if (n_body_on .ne. 0) then
                add_onebody2det = .true.
            else
                add_onebody2det = .false.
            end if
            shift_origin = .true.
            shiftx = .false.
            shifty = .false.
            shiftz = .false.

            iflagerr = 1
            read (5, nml=mesh_info, err=118)
            iflagerr = 0
118         continue
            if (ny .eq. 0) then
                ny = nx
                write (6, *) ' Default value for ny=', ny
            end if
            if (nz .eq. 0) then
                nz = ny
                write (6, *) ' Default value for nz=', nz
            end if
            if (ay .eq. 0. .and. .not. iespbc) then
                ay = ax
                write (6, *) ' Default value for ay=', ay
            end if
            if (az .eq. 0. .and. .not. iespbc) then
                az = ay
                write (6, *) ' Default value for az=', az
            end if
            if (nbufd .eq. -1) then
                nbufd = 1024
                !        if(nelorb.lt.2000) then
                !        nbufd=1000
                !        else
                !        nbufd=100  ! there may be memory problems
                !        endif
                write (6, *) ' Default value for buffer dimension=', nbufd
            end if
        end if
        call checkiflagerr(iflagerr, rankn, 'ERROR reading mesh_info')
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
        add_onebody2det_read = add_onebody2det

        if (iespbc) then
            ax = cellscale(1)/nx
            ay = cellscale(2)/ny
            az = cellscale(3)/nz
            if (rank .eq. 0) write (6, *) ' lattice mesh chosen ', ax, ay, az, Lbox
        end if

        volmesh = ax*ay*az*unit_volume
        volmeshc = dcmplx(volmesh)
        if ((nbufd == 0) .or. (nx == 0) .or. (ny == 0) .or. (nz == 0) .or. volmesh .eq. 0.d0) then
            if (rank .eq. 0) then
                if ((nx == 0) .or. (ny == 0) .or. (nz == 0)) then
                    write (*, *) 'ERROR in the mesh input nx>0,ny>0,nz>0', nx, ny, nz
                end if
                if (nbufd == 0) write (*, *) 'ERROR  buffer dimension nbufd>0', nbufd
                if (volmesh .eq. 0.d0) write (6, *) ' ERROR Volmesh =0 , you should define ax>0,ay>0,az>0'&
                        &, ax, ay, az
                write (*, *) 'Program ends'
            end if
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        mesh = nx
        mesh = mesh*ny
        mesh = mesh*nz

        if (nbufd .lt. mesh) then
            bufdim = nbufd
        else
            bufdim = mesh
        end if
        nbufrep = nprocn*bufdim

        x = 0.d0
        call shift_originref

        if (rank .eq. 0) write (6, *) 'New center of mesh =', rion_ref(:)

#ifdef PARALLEL
        if (.not. bigram) open (100, file=trim(wherescratch)//'.'//chara, form='unformatted')
#else
        rewind (11)
        if (contraction .ne. 0) then
            write (11) nx, ny, nz, ax, ay, az, nbufd, nelorb_c, rion_ref
        else
            write (11) nx, ny, nz, ax, ay, az, nbufd, nelorbh, rion_ref
        end if
#endif
        !       first part calculation overlap determinant

        if (mesh .gt. nbufd) then
            bufdim = nbufd
            allocate (buffer(ipc*ipf*nelorbh, max(ipc, ipf)*nbufd), weightbuf(nbufd))
            if (contraction .ne. 0) allocate (buffer_c(ipc*nelorb_c, max(ipf, ipc)*nbufd))
        else
            bufdim = mesh
            !     Here for ipf=2 buffer is overdimensioned to be consistent with contracted case
            allocate (buffer(ipc*ipf*nelorbh, max(ipc, ipf)*mesh), weightbuf(mesh))
            if (contraction .ne. 0) allocate (buffer_c(ipc*nelorb_c, max(ipf, ipc)*mesh))
        end if
        nbufp = bufdim + 1

        buffer = 0.d0
        weightbuf = 0.d0
        if (contraction .ne. 0) buffer_c = 0.d0
        meshp = (mesh - 1)/proc8 + 1
        nbuf = (meshp - 1)/bufdim + 1
        indr = nbuf
        indr = indr*bufdim
        if (contraction .gt. 0) then
            indr = indr*nelorb_c
        else
            indr = indr*nelorbh
        end if

        if (indr .gt. 2147483647 .and. bigram) then
            if (rank .eq. 0) write (6, *) ' Warning not enough address  per mpi processor in AGP'&
                    &, indr, '<=', 2147483647
            bigram = .false.
        end if

        if (bigram) then
            if (contraction .ne. 0) then
                allocate (buff(ipc*nelorb_c, max(ipc, ipf)*bufdim, nbuf))
            else
                !  For ipf=2 buff is overdimensioned to be consistent with contracted case
                !  In such case buff(ipc*nelorbh+1:2*ipc*nelorbh,:) is identically equal to buff(1:ipc*nelorbh,:)
                !  that in the contracted case may be different.
                allocate (buff(ipc*ipf*nelorbh, max(ipc, ipf)*bufdim, nbuf))
            end if
            buff = 0.d0
        end if

        allocate (overs(ipc*nelorb_cu, ipc*nelorb_cu))
        overs = 0.d0

        iflagnorm = 3
        ind = 0
        !       indr=0
        nbuf = 0
        indproc = 0
        indtot = 0
        !       write(6,*) ' nbufrep before =',rank,nbufrep
#ifdef _OFFLOAD
!$omp target data map(overs) map(to:mu_c) map(alloc:buffer,buffer_c)
#endif
        do k = 1, nz
            !       x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
                !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
                    indtot = indtot + 1
                    !           indr=indr+1
                    !           if(indr-(indr/proc8)*proc8.eq.rank) then

                    if (indproc .eq. rankn) then
                        ind = ind + 1

                        !           x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)

                        x(:) = (-(nz + 1)/2.d0 + k)*az*at(:, 3) + (-(ny + 1)/2.d0 + j)*ay*at(:, 2) &
                               + (-(nx + 1)/2.d0 + i)*ax*at(:, 1) + rion_ref(:)

                        call upnewwf(1, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                &, dupr, zetar, rion, distp, buffer(1, ind), nelorbh, nion, kion            &
                                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                &, indpar_tab, indorb_tab, indshell_tab, .true.)
                        if (ipc .eq. 2) then
                            call upnewwf(1, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                    &, dupr, zetar, rion, distp, buffer(1, bufdim + ind), nelorbh, nion, kion            &
                                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                    &, indpar_tab, indorb_tab, indshell_tab, .false.)
                        end if
                        if (add_onebody2det_read) then
                            psilnn = -scale_one_body
                            do jj = 1, nion
                                if (iespbc) then
                                    rc(:) = x(:) - rion(:, jj)
                                    call CartesianToCrystal(rc, 1)
                                    do kk = 1, 3
                                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                                    end do
                                    r0 = norm_metric(rc, metric)
                                else
                                    rc(:) = (x(:) - rion(:, jj))*costz(jj)
                                    r0 = dsqrt(sum(rc(:)**2))
                                end if
                                psilnn = psilnn - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                            end do
                            buffer(1:ipc*nelorbh, ind) = buffer(1:ipc*nelorbh, ind)*dexp(psilnn)
                            if (ipc .eq. 2) buffer(1:2*nelorbh, bufdim + ind) = buffer(1:2*nelorbh, bufdim + ind)*dexp(psilnn)

                        end if
                        call upvpot_ei(x, zetar, rion, weightbuf(ind), nion        &
                                &, LBox, epsvpot)
                    end if ! endif indproc

                    if (indtot .eq. nbufrep .or. (i .eq. nx .and. j .eq. ny .and. k .eq. nz)) then
                        if (ipf .eq. 2) then
                            if (ipc .eq. 1) then
                                buffer(nelorbh + 1:2*nelorbh, nbufp:2*bufdim) = buffer(1:nelorbh, 1:bufdim)
                                buffer(nelorbh + 1:2*nelorbh, 1:bufdim) = 0.d0
                                buffer(1:nelorbh, nbufp:2*bufdim) = 0.d0
                            else
                                buffer(2*nelorbh + 1:4*nelorbh, nbufp:2*bufdim) = buffer(1:2*nelorbh, nbufp:2*bufdim)
                                buffer(1:2*nelorbh, nbufp:2*bufdim) = 0.d0
                                buffer(2*nelorbh + 1:4*nelorbh, 1:bufdim) = 0.d0
                            end if
                        end if
                        if (contraction .eq. 0) then
                            if (bigram) then
                                buff(:, :, nbuf + 1) = buffer(:, :)
                            else
#ifdef PARALLEL
                                write (100) buffer
                                flush (100)
#else
                                write (11) buffer
#endif
                            end if
                            !   scale after writing as in the contracted case
                            do ii = 1, ind
                                call dscal(ipc*nelorbh*ipf, weightbuf(ii), buffer(1, ii), 1)
                                if (ipc .eq. 2) call dscal(2*nelorbh*ipf, weightbuf(ii), buffer(1, ii + bufdim), 1)
                            end do
#ifdef _OFFLOAD
!$omp target update to (buffer)
#endif
                            if (ipc .eq. 1) then
                                call dgemm_('N', 'T', nelorbh, nelorbh, ind, volmesh, buffer  &
                              &, ipf*nelorbh, buffer, ipf*nelorbh, 1.d0, overs, nelorb_cu)
                            else
                                call zgemm_('N', 'C', nelorbh, nelorbh, ind, volmeshc, buffer  &
                              &, ipf*nelorbh, buffer, ipf*nelorbh, zone, overs, nelorb_cu)
                                if (ipf .eq. 2) then
                                    call zgemm_('N', 'C', nelorbh, nelorbh, ind, volmeshc, buffer(2*nelorbh + 1, nbufp) &
                                                , ipf*nelorbh, buffer(2*nelorbh + 1, nbufp), ipf*nelorbh, zone &
                                                , overs(2*nelorbh + 1, nelorbh + 1), nelorb_cu)
                                else
                                    call zgemm_('N', 'C', nelorbh, nelorbh, ind, volmeshc, buffer(1, nbufp)  &
                                  &, nelorbh, buffer(1, nbufp), nelorbh, zone, overs(1, nelorb_cu + 1), nelorb_cu)
                                end if
                            end if
                        else ! contraction > 0 below
#ifdef _OFFLOAD
!$omp target update to (buffer)
#endif

                            !               if(ipf.eq.2) buffer(ipc*nelorbh+1:2*ipc*nelorbh,:)=buffer(1:ipc*nelorbh,:)
                            !!!CCC guarda bk.f90
                            if (ipc .eq. 1) then
                                call dgemm_('T', 'N', nelorb_c, ind, nelorbh, 1.d0, mu_c&
                              &, ipf*nelorbh, buffer, ipf*nelorbh, 0.d0, buffer_c, nelorb_c)
                                if (ipf .eq. 2) then
                                    call dgemm_('T', 'N', nelorb_c, ind, nelorbh, 1.d0, mu_c(nelorbh + 1, 1)&
                                  &, ipf*nelorbh, buffer(nelorbh + 1, nbufp), ipf*nelorbh, 0.d0, buffer_c(1, nbufp), nelorb_c)
                                end if
                            else

                                if (ipf .eq. 2) then
                                    call zgemm_('T', 'N', nelorb_c, ind, nelorbh, zone&
                                   &, mu_c, 2*nelorbh, buffer, ipf*nelorbh, zzero&
                                   &, buffer_c, nelorb_c)
                                    call zgemm_('T', 'N', nelorb_c, ind, nelorbh, zone&
                                   &, mu_c(2*nelorbh + 1, 1), 2*nelorbh, buffer(2*nelorbh + 1, nbufp)&
                                   &, ipf*nelorbh, zzero, buffer_c(1, nbufp), nelorb_c)
                                    ! CCC Controlla la definizione di questo buffer
                                    ! buffer_c(1:ipc*nelorb_c,nbufp:nbufp+ind-1)=buffer_c(1:ipc*nelorb_c,1:ind)
                                else

                                    call zgemm_('T', 'N', nelorb_c/ipf, ind, nelorbh, zone, mu_c&
                                   &, ipf*nelorbh, buffer, ipf*nelorbh, zzero, buffer_c, nelorb_c)

                                    call zgemm_('T', 'N', nelorb_c, ind, nelorbh, zone, mu_c        &
                                   &, nelorbh, buffer(1, nbufp), ipf*nelorbh, zzero, buffer_c(1, nbufp), nelorb_c)
                                end if

                            end if
#ifdef _OFFLOAD
!$omp target update from (buffer_c)
#endif
                            if (bigram) then
                                buff(:, :, nbuf + 1) = buffer_c(:, :)
                            else
#ifdef PARALLEL
                                write (100) buffer_c
                                flush (100)
#else
                                write (11) buffer_c
#endif
                            end if
                            do ii = 1, ind
                                call dscal(ipc*nelorb_c, weightbuf(ii), buffer_c(1, ii), 1)
                                if (ipc .eq. 2 .or. ipf .eq. 2) then
                                    call dscal(ipc*nelorb_c, weightbuf(ii), buffer_c(1, ii + bufdim), 1)
                                end if
                            end do
#ifdef _OFFLOAD
!$omp target update to (buffer_c)
#endif
                            if (ipc .eq. 1) then
                                call dgemm_('N', 'T', nelorb_c, nelorb_c, ind, volmesh         &
                             &, buffer_c, nelorb_c, buffer_c, nelorb_c, 1.d0, overs, nelorb_cu)
                                if (ipf .eq. 2) then
                                    call dgemm_('N', 'T', nelorb_c, nelorb_c, ind, volmesh         &
                                 &, buffer_c(1, nbufp), nelorb_c, buffer_c(1, nbufp), nelorb_c&
                                 &, 1.d0, overs, nelorb_cu)
                                end if
                            else
                                call zgemm_('N', 'C', nelorb_c, nelorb_c, ind, volmeshc, buffer_c &
                                            , nelorb_c, buffer_c, nelorb_c, zone, overs, nelorb_cu)
                                if (ipf .eq. 2) then
                                    call zgemm_('N', 'C', nelorb_c, nelorb_c, ind, volmeshc, buffer_c(1, nbufp) &
                                                , nelorb_c, buffer_c(1, nbufp), nelorb_c, zone, overs, nelorb_cu)
                                else
                                    call zgemm_('N', 'C', nelorb_c, nelorb_c, ind, volmeshc, buffer_c(1, nbufp) &
                                                , nelorb_c, buffer_c(1, nbufp), nelorb_c, zone, overs(1, nelorb_cu + 1) &
                                                , nelorb_cu)
                                end if
                            end if

                        end if

                        nbuf = nbuf + 1
                        ind = 0
                        indtot = 0
                    end if ! endif buffer load
                    indproc = indproc + 1
                    if (indproc .eq. nprocn) indproc = 0
                end do
            end do
        end do
#ifdef _OFFLOAD
!$omp end target data
#endif
#ifdef PARALLEL
        call mpi_allreduce(ind, indmax, 1                     &
       &, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
        indmax = ind
#endif
        if (indmax .ne. 0) then
            write (6, *) ' ERROR load_fort10in AGP check input nbufd and/or code '
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if
        !CCC
        if (ipf .eq. 2) then
            if (contraction .eq. 0 .and. ipc .eq. 1) then
                overs(nelorb_c/2 + 1:nelorb_c, nelorb_c/2 + 1:nelorb_c) = &
                    overs(1:nelorb_c/2, 1:nelorb_c/2)
            elseif (ipc .eq. 2) then
                do k = nelorb_c + 1, nelorb_cu
                    overs(2*k - 1, k) = 1.d0
                end do
                overs(1:2*nelorb_cu, nelorb_cu + 1:2*nelorb_cu) = overs(1:2*nelorb_cu, 1:nelorb_cu)
            end if
            if (ipc .eq. 1) then
                do k = nelorb_c + 1, nelorb_cu
                    overs(k, k) = 1.d0
                end do
            end if
        end if
#ifdef PARALLEL
        nelorb2 = ipc*ipc*nelorb_cu*nelorb_cu
!           collect overs
        call reduce_base_real(nelorb2, overs, MPI_COMM_WORLD, -1)
#endif

        !     Standard definition scalar product
        if (ipc .eq. 2) call conjmat(nelorb_cu, 2*nelorb_cu, overs, nelorb_cu)

        !     save original contracted matrices  and dimensions
        !    nelorb_c_sav=nelorb_c
        if (contraction .gt. 0) then
            if (npar_eagp .gt. 0) then
                call copy_eagp(.true., ipc, nelorb_c, nelcol_c, detmat_c, eagp_pfaff, detmatc_in)
                !       nelorb_c=nelcol_c
            else
                detmatc_in = detmat_c
            end if
        else
            if (npar_eagp .gt. 0) then
                call copy_eagp(.true., ipc, nelorb_c, nelcol_c, detmat, eagp_pfaff, detmatc_in)
                !       nelorb_c=nelcol_c
            else
                detmatc_in = detmat
            end if
        end if

        deallocate (psip)
        if ((symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) .and. .not. overlap) then
            allocate (psip(ipc*nelorb_c*nelorb_c))
        else
            allocate (psip(2*ipc*nelorb_cu*nelorb_cu))
        end if
        psip = 0.d0

        if (contraction .gt. 0 .or. npar_eagp .gt. 0) then
            ! compute S_in*lambda_in
            if (ipc .eq. 1) then
                call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, overs       &
                        &, nelorb_cu, detmatc_in, nelorb_cu, 0.d0, psip, nelorb_cu, nprocu, rank, comm_mpi)
            else
                ! S_L^T lambda^*
                call conjmat(nelorb_cu, nelorb_cu, detmatc_in, nelorb_cu)
                call zgemm_my('T', 'N', nelorb_cu, nelorb_cu, nelorb_cu, zone, overs&
                        &, nelorb_cu, detmatc_in, nelorb_cu, zzero, psip, nelorb_cu, nprocu, rank, comm_mpi)
                call conjmat(nelorb_cu, nelorb_cu, detmatc_in, nelorb_cu)
            end if
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                overo = tracemat(nelorb_cu, psip, nelorb_cu)
                if (overlap) then
                    ! compute S_in*lambda_in*S_in
                    call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, psip       &
                            &, nelorb_cu, overs, nelorb_cu, 0.d0, psip(nelorb_cu*nelorb_cu + 1), nelorb_cu, nprocu, rank, comm_mpi)
                    ! compute S_in*lambda_in*S_in*lambda_in
                    call dgemm_my('N', 'N', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, psip(nelorb_cu*nelorb_cu + 1)       &
                            &, nelorb_cu, detmatc_in, nelorb_cu, 0.d0, psip, nelorb_cu, nprocu, rank, comm_mpi)
                    ! compute Tr [S_in*lambda_in*S_in*lambda_in]
                    overo_norm = tracemat(nelorb_cu, psip, nelorb_cu)
                end if

            else
                if (ipc .eq. 1) then
                    call dgemm_my('N', 'T', nelorb_cu, nelorb_cu, nelorb_cu, 1.d0, overs&
                            &, nelorb_cu, detmatc_in, nelorb_cu, 0.d0, psip(nelorb_cu*nelorb_cu + 1), nelorb_cu&
                            &, nprocu, rank, comm_mpi)
                else
                    ! S_R \lambda^T
                    call zgemm_my('N', 'T', nelorb_cu, nelorb_cu, nelorb_cu, zone&
                            &, overs(1, nelorb_cu + 1), nelorb_cu, detmatc_in, nelorb_cu, zzero&
                            &, psip(2*nelorb_cu*nelorb_cu + 1), nelorb_cu, nprocu, rank, comm_mpi)
                end if
                overo = tracemat2(nelorb_cu, nelorb_cu, psip, nelorb_cu&
                        &, psip(ipc*nelorb_cu*nelorb_cu + 1), nelorb_cu)
                ! opzione overlap NON implementata in questo caso
            end if ! symmagp

        else ! contraction =0 below
            ! compute S_in*lambda_in
            if (ipc .eq. 1) then
                call dgemm_my('N', 'N', nelorbpf, nelorbpf, nelorbpf, 1.d0, overs          &
                        &, nelorbpf, detmat, nelorbpf, 0.d0, psip, nelorbpf, nprocu, rank, comm_mpi)
            else
                !  S_L^T \lambda^*
                call conjmat(nelorbpf, nelorbpf, detmat, nelorbpf)
                call zgemm_my('T', 'N', nelorbpf, nelorbpf, nelorbpf, zone, overs          &
                        &, nelorbpf, detmat, nelorbpf, zzero, psip, nelorbpf, nprocu, rank, comm_mpi)
                call conjmat(nelorbpf, nelorbpf, detmat, nelorbpf)
            end if

            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                overo = tracemat(nelorbh, psip, nelorbh)

                if (overlap) then
                    ! compute S_in*lambda*S_in
                    call dgemm_my('N', 'N', nelorbh, nelorbh, nelorbh, 1.d0, psip          &
                            &, nelorbh, overs, nelorbh, 0.d0, psip(nelorbh*nelorbh + 1), nelorbh, nprocu, rank, comm_mpi)
                    ! compute S_in*lambda_in*S_in*lambda_in
                    call dgemm_my('N', 'N', nelorbh, nelorbh, nelorbh, 1.d0, psip(nelorbh*nelorbh + 1)      &
                            &, nelorbh, detmat, nelorbh, 0.d0, psip, nelorbh, nprocu, rank, comm_mpi)
                    ! compute Tr [S_in*lambda_in*S_in*lambda_in]
                    overo_norm = tracemat(nelorb_c, psip, nelorb_c)
                end if

            else

                if (ipc .eq. 1) then
                    call dgemm_my('N', 'T', nelorbpf, nelorbpf, nelorbpf, 1.d0, overs          &
                            &, nelorbpf, detmat, nelorbpf, 0.d0, psip(nelorbpf*nelorbpf + 1), nelorbpf, nprocu, rank, comm_mpi)
                else
                    ! S_R \lambda^T
                    call zgemm_my('N', 'T', nelorbpf, nelorbpf, nelorbpf, zone&
                            &, overs(1, nelorbpf + 1), nelorbpf, detmat, nelorbpf, zzero&
                            &, psip(2*nelorbpf*nelorbpf + 1), nelorbpf, nprocu, rank, comm_mpi)
                end if
                overo = tracemat2(nelorbpf, nelorbpf, psip, nelorbpf, psip(ipc*nelorbpf*nelorbpf + 1), nelorbpf)

            end if

            ! opzione overlap NON implementata in questo caso

        end if

        !         end about determinant

        if (nelorbj_c .gt. 0 .and. change_jas) then
            !       second  part calculation overlap jastrow
            if (allocated(buffer)) deallocate (buffer)
            if (allocated(weightbuf)) deallocate (weightbuf)
            if (allocated(buffer_c)) deallocate (buffer_c)

            if (mesh .gt. nbufd) then
                bufdim = nbufd
                allocate (buffer(nelorbjh, nbufd), weightbuf(nbufd))
                if (contractionj .ne. 0) allocate (buffer_c(nelorbj_c, nbufd))
            else
                bufdim = mesh
                allocate (buffer(nelorbjh, mesh), weightbuf(mesh))
                if (contractionj .ne. 0) allocate (buffer_c(nelorbj_c, mesh))
            end if

            nbuf = (meshp - 1)/bufdim + 1
            if (contractionj .gt. 0) then
                indr = nbuf
                indr = indr*bufdim
                indr = indr*nelorbj_c
            else
                indr = nbuf
                indr = indr*bufdim
                indr = indr*nelorbjh
            end if

            if (indr .gt. 2147483647 .and. bigram) then
                if (rank .eq. 0) write (6, *) ' Warning not enough address  per mpi processor in Jas '&
                        &, indr, '<=', 2147483647
                bigram = .false.
            end if
            if (bigram) then
                if (contractionj .ne. 0) then
                    allocate (buffj(nelorbj_c, bufdim, nbuf))
                else
                    allocate (buffj(nelorbjh, bufdim, nbuf))
                end if
            end if
            allocate (oversj(nelorbj_c, nelorbj_c))

            buffer = 0.d0
            weightbuf = 0.d0
            if (contractionj .ne. 0) buffer_c = 0.d0
            oversj = 0.d0

#ifdef PARALLEL
            if (.not. bigram) open (101, file=trim(wherescratch)//'buf.'//chara, form='unformatted')
#else
            rewind (12)
            if (contractionj .ne. 0) then
                ! first rec
                write (12) nx, ny, nz, ax, ay, az, nbufd, nelorbj_c, rion_ref
            else
                ! first reco
                write (12) nx, ny, nz, ax, ay, az, nbufd, nelorbjh, rion_ref
            end if
#endif

            if (contractionj .ne. 0) then
                do i = 1, nelorbj_c
                    if (orbcostn(i) .and. i .ne. adrcostc) then
                        muj_c(:, i) = 0.d0
                    end if
                end do
                mujc_in = muj_c
            end if

            oversj = 0.d0
            iflagnorm = 3
            ind = 0
            !       indr=0
            indtot = 0
            nbuf = 0
            indproc = 0
#ifdef _OFFLOAD
!$omp target data map(oversj) map(to:muj_c) map(alloc:buffer,buffer_c)
#endif
            do k = 1, nz
                !       x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
                do j = 1, ny
                    !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                    do i = 1, nx
                        !           indr=indr+1
                        !           if(indr-(indr/proc8)*proc8.eq.rank) then
                        indtot = indtot + 1
                        if (rankn .eq. indproc) then
                            ind = ind + 1
                            !           x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)
                            x(:) = (-(nz + 1)/2.d0 + k)*az*at(:, 3) + (-(ny + 1)/2.d0 + j)*ay*at(:, 2) &
                                   + (-(nx + 1)/2.d0 + i)*ax*at(:, 1) + rion_ref(:)

                            call upnewwf(1, 0, 0, 1, nshelljr, ioptorbj, ioccj, x, 1, r, rmu       &
                                    &, vjur, zetar, rion, distp, buffer(1, ind), nelorbjh, nion, kionj          &
                                    &, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                                    &, indparj_tab, indorbj_tab, indshellj_tab, .true.)

                            call upvpot_ei(x, zetar, rion, weightbuf(ind), nion        &
                                    &, LBox, epsvpot)

                        end if ! endif indproc

                        if (indtot .eq. nbufrep .or. (i .eq. nx .and. j .eq. ny .and. k .eq. nz)) then
                            if (contractionj .eq. 0) then
                                if (bigram) then
                                    buffj(:, :, nbuf + 1) = buffer(:, :)
                                else
#ifdef PARALLEL
                                    write (101) buffer
                                    flush (101)
#else
                                    write (12) buffer
#endif
                                end if
                                do ii = 1, ind
                                    call dscal(nelorbjh, weightbuf(ii), buffer(1, ii), 1)
                                end do
#ifdef  _OFFLOAD
!$omp target update to (buffer)
#endif
                                call dgemm_('N', 'T', nelorbj_c, nelorbj_c, ind, volmesh, buffer&
                                        &, nelorbjh, buffer, nelorbjh, 1.d0, oversj, nelorbj_c)
                            else
#ifdef  _OFFLOAD
!$omp target update to (buffer)
#endif
                                call dgemm_('T', 'N', nelorbj_c, ind, nelorbjh, 1.d0, muj_c     &
                                            &, nelorbjh, buffer, nelorbjh, 0.d0, buffer_c, nelorbj_c)
#ifdef  _OFFLOAD
!$omp target update from (buffer_c)
#endif
                                if (bigram) then
                                    buffj(:, :, nbuf + 1) = buffer_c(:, :)
                                else
#ifdef PARALLEL
                                    write (101) buffer_c
                                    flush (101)
#else
                                    write (12) buffer_c
#endif
                                end if
                                do ii = 1, ind
                                    call dscal(nelorbj_c, weightbuf(ii), buffer_c(1, ii), 1)
                                end do
#ifdef  _OFFLOAD
!$omp target update to (buffer_c)
#endif
                                call dgemm_('N', 'T', nelorbj_c, nelorbj_c, ind, volmesh       &
                                                   &, buffer_c, nelorbj_c, buffer_c, nelorbj_c, 1.d0, oversj, nelorbj_c)
                            end if
                            nbuf = nbuf + 1
                            ind = 0
                            indtot = 0
                        end if ! endif load buf
                        indproc = indproc + 1
                        if (indproc .eq. nprocn) indproc = 0
                    end do
                end do
            end do
#ifdef _OFFLOAD
!$omp end target data
#endif
#ifdef PARALLEL
            call mpi_allreduce(ind, indmax, 1                     &
           &, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
            indmax = ind
#endif
            if (indmax .ne. 0) then
                write (6, *) ' ERROR load_fort10in Jas check input nbufd and/or code '
#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if
            nelorb2 = nelorbj_c*nelorbj_c
#ifdef  PARALLEL
!        collect oversj
            call reduce_base_real(nelorb2, oversj, MPI_COMM_WORLD, -1)
#endif

            !        calculation overlap old
            deallocate (psip)
            allocate (psip(nelorbj_c*nelorbj_c))
            psip = 0.d0
            if (contractionj .gt. 0) then
                call setorbcost(nelorbj_c, jasmat_c, orbcostn, nozeroj_c, nnozeroj_c)
                call purify(jasmat_c, nelorbj_c, oversj, nelorbj_c, orbcostn)
                jasmatc_in = jasmat_c
                call dgemm_my('N', 'N', nelorbj_c, nelorbj_c, nelorbj_c, 1.d0, oversj   &
                        &, nelorbj_c, jasmat_c, nelorbj_c, 0.d0, psip, nelorbj_c, nprocu, rank, comm_mpi)
                overoj = tracemat(nelorbj_c, psip, nelorbj_c)

                overojnloc = tracematnloc(nelorbj_c, psip, nelorbj_c, orbcostn)
                overojloc = tracematloc(nelorbj_c, oversj, nelorbj_c, jasmat_c, nelorbj_c, orbcostn)

                if (iessz) then
                    call setorbcost(nelorbj_c, jasmatsz_c, orbcostn, nozeroj_c, nnozeroj_c)
                    call purify(jasmatsz_c, nelorbj_c, oversj, nelorbj_c, orbcostn)
                    jasmatszc_in = jasmatsz_c
                    call dgemm_my('N', 'N', nelorbj_c, nelorbj_c, nelorbj_c, 1.d0, oversj   &
                            &, nelorbj_c, jasmatsz_c, nelorbj_c, 0.d0, psip, nelorbj_c, nprocu, rank, comm_mpi)
                    overojsznloc = tracematnloc(nelorbj_c, psip, nelorbj_c, orbcostn)
                    overojsz = tracemat(nelorbj_c, psip, nelorbj_c)
                    overojszloc = tracematloc(nelorbj_c, oversj, nelorbj_c, jasmatsz_c, nelorbj_c&
                            &, orbcostn)
                end if

            else

                call setorbcost(nelorbjh, jasmat, orbcostl, nozeroj, nnozeroj)
                call purify(jasmat, nelorbjh, oversj, nelorbj, orbcostl)
                jasmatc_in = jasmat

                call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, 1.d0, oversj&
                        &, nelorbjh, jasmat, nelorbjh, 0.d0, psip, nelorbjh, nprocu, rank, comm_mpi)
                overoj = tracemat(nelorbjh, psip, nelorbjh)
                overojnloc = tracematnloc(nelorbjh, psip, nelorbjh, orbcostn)
                overojloc = tracematloc(nelorbjh, oversj, nelorbjh, jasmat, nelorbjh, orbcostn)
                if (iessz) then
                    call setorbcost(nelorbjh, jasmatsz, orbcostl, nozeroj, nnozeroj)
                    call purify(jasmatsz, nelorbjh, oversj, nelorbj, orbcostl)
                    jasmatszc_in = jasmatsz
                    call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, 1.d0, oversj      &
                            &, nelorbjh, jasmatsz, nelorbjh, 0.d0, psip, nelorbjh, nprocu, rank, comm_mpi)
                    overojsznloc = tracematnloc(nelorbjh, psip, nelorbjh, orbcostn)
                    overojsz = tracemat(nelorbjh, psip, nelorbjh)
                    overojszloc = tracematloc(nelorbjh, oversj, nelorbjh, jasmatsz, nelorbjh, orbcostn)
                end if
            end if
        end if
        !         end about jastrow

        if (allocated(buffer)) deallocate (buffer)
        if (allocated(weightbuf)) deallocate (weightbuf)
        if (allocated(buffer_c)) deallocate (buffer_c)

        close (10)
    end subroutine load_fort10in

    subroutine load_gemz
        use allio, only: norm_metric
        implicit none
        real*8 :: psilnn, r0, rc(3)
        integer :: indmax
        real*8, external :: jastrow_ei

        !        input dupr  ! remember to change the derivatives
        if (.not. bigram) then
#ifdef PARALLEL
            rewind (100)
#else
            rewind (11)
            ! read first useless record
            read (11)
#endif
        end if

        oversnl = 0.d0
        overonl = 0.d0
        iflagnorm = 3
        ind = 0
        indtot = 0
        indproc = 0
        !       indr=0
        nbuf = 0
#ifdef _OFFLOAD
!$omp target data map(oversnl,overonl) map(alloc:buffer,buffero)
#endif
        do k = 1, nz
            !       x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
                !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
                    indtot = indtot + 1
                    !           indr=indr+1
                    !           if(indr-(indr/proc8)*proc8.eq.rank) then
                    if (indproc .eq. rankn) then
                        ind = ind + 1
                        !           x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)
                        x(:) = (-(nz + 1)/2.d0 + k)*az*at(:, 3) + (-(ny + 1)/2.d0 + j)*ay*at(:, 2) &
                               + (-(nx + 1)/2.d0 + i)*ax*at(:, 1) + rion_ref(:)

                        call upnewwf(1, 0, 0, 1, nshell, ioptorb, ioccup, x, 1, r, rmu         &
                                &, dupr, zetar, rion, distp, buffer(1, ind), nelorb, nion, kion             &
                                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                &, indpar_tab, indorb_tab, indshell_tab, .true.)
                        if (ipc .eq. 2) then
                            call upnewwf(1, 0, 0, 1, nshell, ioptorb, ioccup, x, 1, r, rmu         &
                                    &, dupr, zetar, rion, distp, buffer(1, ind + bufdim), nelorb, nion, kion             &
                                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                    &, indpar_tab, indorb_tab, indshell_tab, .false.)
                        end if

                        if (add_onebody2det_read) then
                            psilnn = -scale_one_body
                            do jj = 1, nion
                                if (iespbc) then
                                    rc(:) = x(:) - rion(:, jj)
                                    call CartesianToCrystal(rc, 1)
                                    do kk = 1, 3
                                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                                    end do
                                    r0 = norm_metric(rc, metric)
                                else
                                    rc(:) = (x(:) - rion(:, jj))*costz(jj)
                                    r0 = dsqrt(sum(rc(:)**2))
                                end if
                                psilnn = psilnn - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                            end do
                            buffer(1:ipc*nelorbh, ind) = buffer(1:ipc*nelorbh, ind)*dexp(psilnn)
                            if (ipc .eq. 2) buffer(1:2*nelorbh, bufdim + ind) = &
                                    &buffer(1:2*nelorbh, bufdim + ind)*dexp(psilnn)
                        end if

                        call upvpot_ei(x, zetar, rion, weightbuf(ind), nion, LBox, epsvpot)
                        call dscal(ipc*nelorb, weightbuf(ind), buffer(1, ind), 1)
                        if (ipc .eq. 2) call dscal(2*nelorb, weightbuf(ind), buffer(1, bufdim + ind), 1)

                    end if ! endif indproc

                    if (indtot .eq. nbufrep .or. (i .eq. nx .and. j .eq. ny .and. k .eq. nz)) then
#ifdef _OFFLOAD
!$omp target  update to (buffer)
#endif
                        if (ipc .eq. 1) then
                            call dgemm_('N', 'T', nelorb, nelorb, ind, volmesh, buffer      &
                                    &, nelorb, buffer, nelorb, 1.d0, oversnl, nelorbpf)
                        else
                            call zgemm_('N', 'C', nelorb, nelorb, ind, volmeshc, buffer      &
                                    &, nelorb, buffer, nelorb, zone, oversnl, nelorbpf)
                            if (ipf .eq. 2) then
                                call zgemm_('N', 'C', nelorb, nelorb, ind, volmeshc, buffer(1, nbufp)  &
                                        &, nelorb, buffer(1, nbufp), nelorb, zone, oversnl(2*nelorbh + 1, nelorbh + 1), nelorbpf)
                            else
                                call zgemm_('N', 'C', nelorb, nelorb, ind, volmeshc, buffer(1, nbufp)  &
                                        &, nelorb, buffer(1, nbufp), nelorb, zone, oversnl(1, nelorb + 1), nelorb)
                            end if
                        end if
                        nbuf = nbuf + 1

                        if (bigram) then
                            buffero(:, :) = buff(:, :, nbuf)
                        else
#ifdef PARALLEL
                            read (100) buffero

#else
                            read (11) buffero
#endif
                        end if
                        do ii = 1, ind
                            call dscal(ipc*nelorbc_in, weightbuf(ii), buffero(1, ii), 1)
                            if (ipc .eq. 2 .or. ipf .eq. 2) then
                                call dscal(ipc*nelorbc_in, weightbuf(ii), buffero(1, ii + bufdim), 1)
                            end if
                        end do
#ifdef _OFFLOAD
!$omp target update to (buffero)
#endif
                        if (ipc .eq. 1) then
                            call dgemm_('N', 'T', nelorbc_in, nelorb, ind, volmesh         &
                                    &, buffero, nelorbc_in, buffer, nelorb, 1.d0, overonl, nelorbcu_in)
                            if (ipf .eq. 2) then
                                call dgemm_('N', 'T', nelorbc_in, nelorb, ind, volmesh, buffero(1, nbufp) &
                                            , nelorbc_in, buffer, nelorb, 1.d0, overonl(1, nelorbh + 1), nelorbcu_in)
                            end if
                        else
                            call zgemm_('N', 'C', nelorbc_in, nelorb, ind, volmeshc         &
                                    &, buffero, nelorbc_in, buffer, nelorb, zone, overonl, nelorbcu_in)
                            if (ipf .eq. 2) then
                                call zgemm_('N', 'C', nelorbc_in, nelorb, ind, volmeshc         &
                                        &, buffero(1, nbufp), nelorbc_in, buffer(1, nbufp), nelorb&
                                        &, zone, overonl(1, nelorb + 1), nelorbcu_in)
                            else
                                call zgemm_('N', 'C', nelorbc_in, nelorb, ind, volmeshc, buffero(1, nbufp) &
                                            , nelorbc_in, buffer(1, nbufp), nelorb, zone, overonl(1, nelorb + 1), nelorbcu_in)
                            end if
                        end if

                        ind = 0
                        indtot = 0
                    end if ! endif load buf
                    indproc = indproc + 1
                    if (indproc .eq. nprocn) indproc = 0
                end do
            end do
        end do
#ifdef _OFFLOAD
!$omp end target data
#endif
#ifdef PARALLEL
        call mpi_allreduce(ind, indmax, 1&
      &, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
        indmax = ind
#endif
        if (indmax .ne. 0) then
            write (6, *) ' ERROR load_gemz check input nbufd and/or code '
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        !       collect overonl oversnl
        if (ipf .eq. 2) then
            if (ipc .eq. 1) then
                oversnl(nelorbh + 1:2*nelorbh, nelorbh + 1:2*nelorbh) = &
                    oversnl(1:nelorbh, 1:nelorbh)
                do k = ipf*nelorbh + 1, nelorbpf
                    oversnl(k, k) = 1.d0
                end do
                do k = ipf*nelorbh + 1, nelorbpf
                    overonl(k - ipf*nelorbh + nelorbc_in, k) = 1.d0
                end do
            else
                do k = ipf*nelorbh + 1, nelorbpf
                    oversnl(2*k - 1, k) = 1.d0
                    overonl(2*(k - ipf*nelorbh) - 1 + 2*nelorbc_in, k) = 1.d0
                end do
                oversnl(1:2*nelorbpf, nelorbpf + 1:2*nelorbpf) = &
                    oversnl(1:2*nelorbpf, 1:nelorbpf)
                overonl(1:2*nelorbcu_in, nelorbpf + 1:2*nelorbpf) = &
                    overonl(1:2*nelorbcu_in, 1:nelorbpf)
            end if
        end if
#ifdef PARALLEL
        nelorb2 = ipc*ipc*nelorbpf*nelorbpf
        call reduce_base_real(nelorb2, oversnl, MPI_COMM_WORLD, -1)
        nelorb2 = ipc*ipc*nelorbcu_in*nelorbpf
        call reduce_base_real(nelorb2, overonl, MPI_COMM_WORLD, -1)
#endif
        if (ipc .eq. 2) then
            call conjmat(nelorbpf, 2*nelorbpf, oversnl, nelorbpf)
            call conjmat(nelorbcu_in, 2*nelorbpf, overonl, nelorbcu_in)
        end if

        !   determination best matrix approximating fort.10_in in the new basis

        over_sav = oversnl

        if (ipc .eq. 1) then
            call dgemm_my('N', 'N', nelorbcu_in, nelorbpf, nelorbcu_in, 1.d0, detmatc_in &
                    &, nelorbcu_in, overonl, nelorbcu_in, 0.d0, mat, nelorbcu_in, nprocu, rank, comm_mpi)
        else
            call invsymeps(ipc, nelorbpf, oversnl(1, nelorbpf + 1), nelorbpf&
                    &, info, epsdgel, mine, umat, eigmat, nprocu, rank, comm_mpi)
        end if
        call invsymeps(ipc, nelorbpf, oversnl, nelorbpf, info, epsdgel&
                &, mine, umat, eigmat, nprocu, rank, comm_mpi)
        if (info .ne. 0) then
            if (rank .eq. 0) write (6, *) ' SDV  failed  !!! ', info
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        if (ipc .eq. 1) then
            call dgemm_my('N', 'N', nelorbcu_in, nelorbpf, nelorbpf, 1.d0, overonl&
                    &, nelorbcu_in, oversnl, nelorbpf, 0.d0, psip, nelorbcu_in, nprocu, rank, comm_mpi)
        else
            !      different algorithm

            ! \bar S_L (S^\prime_L)^-1
            call zgemm_my('N', 'N', nelorbcu_in, nelorbpf, nelorbpf, zone, overonl&
                    &, nelorbcu_in, oversnl, nelorbpf, zzero, psip, nelorbcu_in, nprocu, rank, comm_mpi)
            !
            !   lambda \bar S_R^*
            call conjmat(nelorbcu_in, nelorbpf, overonl(1, nelorbpf + 1), nelorbcu_in)
            call zgemm_my('N', 'N', nelorbcu_in, nelorbpf, nelorbcu_in, zone, detmatc_in&
                    &, nelorbcu_in, overonl(1, nelorbpf + 1), nelorbcu_in, zzero, mat&
                    &, nelorbcu_in, nprocu, rank, comm_mpi)
            call conjmat(nelorbcu_in, nelorbpf, overonl(1, nelorbpf + 1), nelorbcu_in)
            !  umat=A = psip^dag mat= (S^\prime_L)^-1 \bar S_L^dag lambda \bar S_R^*
            call zgemm_my('C', 'N', nelorbpf, nelorbpf, nelorbcu_in, zone, psip&
                    &, nelorbcu_in, mat, nelorbcu_in, zzero&
                    &, psip(2*nelorbcu_in*nelorbpf + 1), nelorbpf, nprocu, rank, comm_mpi)
            ! umat = psip(2*nelorbc_in*nelorbh+1) umat here refers to notes parbcs.
        end if

        if (ipc .eq. 1) inv_sav = oversnl

        if (ipc .eq. 1) then
            call dgemm_my('T', 'N', nelorbpf, nelorbpf, nelorbcu_in, 1.d0, mat&
                    &, nelorbcu_in, psip, nelorbcu_in, 0.d0, oversnl, nelorbpf, nprocu, rank, comm_mpi)
        end if

        if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
            overlapsquare = tracemat(nelorbh, oversnl, nelorbh)/overo
        else
            if (ipc .eq. 1) then
                call dgemm_my('T', 'N', nelorbcu_in, nelorbpf, nelorbcu_in, 1.d0, detmatc_in &
                        &, nelorbcu_in, overonl, nelorbcu_in, 0.d0, mat, nelorbcu_in, nprocu, rank, comm_mpi)

                call dgemm_my('T', 'N', nelorbpf, nelorbpf, nelorbcu_in, 1.d0, mat&
                        &, nelorbcu_in, psip, nelorbcu_in, 0.d0, psip(nelorbcu_in*nelorbpf + 1)&
                        &, nelorbpf, nprocu, rank, comm_mpi)
            else

                !            compute inverse right. Here we cannot assume R and L are the same.
                !   Definition of mat=\bar S_R^T \lambda^dag
                call zgemm_my('T', 'C', nelorbpf, nelorbcu_in, nelorbcu_in, zone&
                        &, overonl(1, nelorbpf + 1), nelorbcu_in, detmatc_in, nelorbcu_in, zzero, mat&
                        &, nelorbpf, nprocu, rank, comm_mpi)
                !    B= oversnl=S_R^*-1 \bar S_R^T \lambda^\dag  \bar S_L
                call conjmat(nelorbpf, nelorbpf, oversnl(1, nelorbpf + 1), nelorbpf)
                call zgemm_my('N', 'N', nelorbpf, nelorbcu_in, nelorbpf, zone&
                        &, oversnl(1, nelorbpf + 1), nelorbpf, mat, nelorbpf, zzero, psip, nelorbpf&
                        &, nprocu, rank, comm_mpi)
                call conjmat(nelorbpf, nelorbpf, oversnl(1, nelorbpf + 1), nelorbpf)
                !  Use oversnl The left inverse is no longer required
                call zgemm_my('N', 'N', nelorbpf, nelorbpf, nelorbcu_in, zone, psip&
                        &, nelorbpf, overonl, nelorbcu_in, zzero, oversnl, nelorbpf, nprocu&
                        &, rankopt, commopt_mpi)

            end if

            if (ipc .eq. 1) then
                overlapsquare = tracemat2(nelorbpf, nelorbpf, oversnl, nelorbpf, psip(nelorbcu_in*nelorbpf + 1), nelorbpf)/overo
            else
                overlapsquare = tracemat2(nelorbpf, nelorbpf, psip(2*nelorbcu_in*nelorbpf + 1), nelorbpf, oversnl, nelorbpf)/overo
            end if

        end if

        if (rank .eq. 0 .and. overo .ne. 0) write (6, *)&
                & ' Overlap square Geminal uncontracted found =', overlapsquare
        !OK

        !       to be closest in L2 norm
        if (ipc .eq. 1) then
            cost = 1.d0
            call dgemm_my('N', 'N', nelorbpf, nelorbpf, nelorbpf, cost, inv_sav&
                    &, nelorbpf, oversnl, nelorbpf, 0.d0, mat, nelorbpf, nprocu, rank, comm_mpi)
        else
            !  The optimal detmat is umat (of notes)^ S^R^T^{-1}
            call zgemm_my('N', 'T', nelorbpf, nelorbpf, nelorbpf, zone&
                    &, psip(2*nelorbcu_in*nelorbpf + 1), nelorbpf, oversnl(1, nelorbpf + 1)&
                    &, nelorbpf, zzero, mat, nelorbpf, nprocu, rank, comm_mpi)
        end if

        !         write(6,*) ' Output matrix ='
        !           do i=1,ipf*nelorbh
        !              do j=i,ipf*nelorbh
        !              write(6,*) i,j,detmat(ipf*ipc*nelorbh*(i-1)+j),detmat(ipf*ipc*nelorbh*(j-1)+i)
        !              enddo
        !           enddo

        if (.not. overlap) then
            detmat = 0.d0
            do i = 1, nelorbpf
                call dcopy(ipc*nelorbpf, mat(1, i), 1, detmat(ipc*(i - 1)*nelorbpf + 1), 1)
            end do
        end if

    end subroutine load_gemz
    subroutine load_jasz
        implicit none
        integer indmax
        if (.not. bigram) then
#ifdef PARALLEL
            rewind (101)
#else
            rewind (12)
            ! read first useless record
            read (12)
#endif
        end if

        oversnjl = 0.d0
        overonjl = 0.d0

        iflagnorm = 3
        ind = 0
        indtot = 0
        indproc = 0
        !       indr=0
        nbuf = 0
#ifdef _OFFLOAD
!$omp target data map(oversnjl,overonjl) map(alloc:buffer,buffero)
#endif
        do k = 1, nz
            !       x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
                !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
                    !           indr=indr+1
                    !         if(indr-(indr/proc8)*proc8.eq.rank) then
                    indtot = indtot + 1
                    if (indproc .eq. rankn) then
                        ind = ind + 1
                        !           x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)
                        x(:) = (-(nz + 1)/2.d0 + k)*az*at(:, 3) + (-(ny + 1)/2.d0 + j)*ay*at(:, 2) &
                               + (-(nx + 1)/2.d0 + i)*ax*at(:, 1) + rion_ref(:)
                        call upnewwf(1, 0, 0, 1, nshellj, ioptorbj, ioccj, x, 1, r, rmu        &
                                &, vjur, zetar, rion, distp, buffer(1, ind), nelorbj, nion, kionj           &
                                &, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                                &, indparj_tab, indorbj_tab, indshellj_tab, .true.)

                        call upvpot_ei(x, zetar, rion, weightbuf(ind), nion        &
                                &, LBox, epsvpot)
                        call dscal(nelorbj, weightbuf(ind), buffer(1, ind), 1)
                    end if ! endif indproc

                    if (indtot .eq. nbufrep .or. (i .eq. nx .and. j .eq. ny .and. k .eq. nz)) then
#ifdef _OFFLOAD
!$omp target update to (buffer)
#endif

                        call dgemm_('N', 'T', nelorbj, nelorbj, ind, volmesh, buffer    &
                                &, nelorbj, buffer, nelorbj, 1.d0, oversnjl, nelorbj)
                        nbuf = nbuf + 1
                        if (bigram) then
                            buffero(:, :) = buffj(:, :, nbuf)
                        else
#ifdef PARALLEL
                            read (101) buffero
#else
                            read (12) buffero
#endif
                        end if
                        do ii = 1, ind
                            call dscal(nelorbjc_in, weightbuf(ii), buffero(1, ii), 1)
                        end do

#ifdef _OFFLOAD
!$omp target  update to (buffero)
#endif
                        call dgemm_('N', 'T', nelorbjc_in, nelorbj, ind, volmesh       &
                                &, buffero, nelorbjc_in, buffer, nelorbj, 1.d0, overonjl, nelorbjc_in)
                        ind = 0
                        indtot = 0
                    end if ! endif load buf
                    indproc = indproc + 1
                    if (indproc .eq. nprocn) indproc = 0
                end do
            end do
        end do
#ifdef _OFFLOAD
!$omp end target data
#endif
#ifdef PARALLEL
        call mpi_allreduce(ind, indmax, 1&
      &, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
        indmax = ind
#endif
        if (indmax .ne. 0) then
            write (6, *) ' ERROR load_jas check input nbufd and/or code '
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

#ifdef PARALLEL
!       collect overonl oversnl
        nelorb2 = nelorbj*nelorbj
        call reduce_base_real(nelorb2, oversnjl, MPI_COMM_WORLD, -1)
        nelorb2 = nelorbjc_in*nelorbj
        call reduce_base_real(nelorb2, overonjl, MPI_COMM_WORLD, -1)
#endif

        over_sav = oversnjl

        call dgemm_my('N', 'N', nelorbjc_in, nelorbjh, nelorbjc_in, 1.d0         &
                &, jasmatc_in, nelorbjc_in, overonjl, nelorbjc_in, 0.d0                 &
                &, mat, nelorbjc_in, nprocu, rank, comm_mpi)

        over_sav(1:nelorbj, 1:nelorbj) = oversnjl(1:nelorbj, 1:nelorbj)
        call invsymeps(1, nelorbjh, oversnjl, nelorbj, info, epsdgel&
                &, mine, umat, eigmat, nprocu, rank, comm_mpi)

        call dgemm_my('N', 'N', nelorbjc_in, nelorbjh, nelorbjh, 1.d0, overonjl   &
                &, nelorbjc_in, oversnjl, nelorbj, 0.d0, psip, nelorbjc_in, nprocu, rank, comm_mpi)

        if (info .ne. 0) then
            if (rank .eq. 0) write (6, *) ' SDV  failed  !!! '
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        inv_sav = oversnjl

        call dgemm_my('T', 'N', nelorbjh, nelorbjh, nelorbjc_in, 1.d0, mat&
                &, nelorbjc_in, psip, nelorbjc_in, 0.d0, oversnjl, nelorbjh, nprocu, rank, comm_mpi)

        overlapsquarej = tracemat(nelorbjh, oversnjl, nelorbjh)/overoj

        if (rank .eq. 0 .and. overoj .ne. 0.d0) write (6, *)&
                & ' Overlap square Jas uncontracted found =', overlapsquarej

        cost = 1.d0
        call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, cost, inv_sav&
                &, nelorbj, oversnjl, nelorbjh, 0.d0, mat, nelorbjh, nprocu, rank, comm_mpi)

        if (.not. overlap) then
            jasmat = 0.d0
            do i = 1, nelorbjh
                call dcopy(nelorbjh, mat(1, i), 1, jasmat(nelorbjh*(i - 1) + 1), 1)
            end do
        end if

        call dcopy(nelorbjh*nelorbjh, oversnjl, 1, psip, 1)

        if (iessz .and. iessz_in) then

            !       recover old def of psip
            call dgemm_my('N', 'N', nelorbjc_in, nelorbjh, nelorbjh, 1.d0, overonjl   &
                    &, nelorbjc_in, inv_sav, nelorbj, 0.d0, psip, nelorbjc_in, nprocu, rank, comm_mpi)

            call dgemm_my('N', 'N', nelorbjc_in, nelorbjh, nelorbjc_in, 1.d0         &
                    &, jasmatszc_in, nelorbjc_in, overonjl, nelorbjc_in, 0.d0               &
                    &, mat, nelorbjc_in, nprocu, rank, comm_mpi)

            call dgemm_my('T', 'N', nelorbjh, nelorbjh, nelorbjc_in, 1.d0, mat           &
                    &, nelorbjc_in, psip, nelorbjc_in, 0.d0, oversnjl, nelorbjh, nprocu, rank, comm_mpi)

            overlapsquarejsz = tracemat(nelorbjh, oversnjl, nelorbjh)/overojsz

            if (rank .eq. 0 .and. overojsz .ne. 0.d0) write (6, *)&
                    & ' Overlap square uncontracted Jas Sz found =', overlapsquarejsz

            !       cost=1.d0/dsqrt(overlapsquarejsz)
            !       to be closest in L2 norm

            cost = 1.d0
            call dgemm_my('N', 'N', nelorbjh, nelorbjh, nelorbjh, cost, inv_sav&
                    &, nelorbj, oversnjl, nelorbjh, 0.d0, mat, nelorbjh, nprocu, rank, comm_mpi)

            if (.not. overlap) then
                jasmatsz = 0.d0
                do i = 1, nelorbjh
                    call dcopy(nelorbjh, mat(1, i), 1, jasmatsz(nelorbjh*(i - 1) + 1), 1)
                end do
            end if

            call dcopy(nelorbjh*nelorbjh, oversnjl, 1, psip, 1)

            !      write(6,*) ' analitycal der ',nelorbjh,nshelljr
            !      do i=1,nshelljr-1
            !      write(6,*) i,derzsz(i),derz(i)
            !      enddo

        else
            if (.not. overlap) then
                jasmatsz_c = 0.d0
                jasmatsz = 0.d0
            end if
        end if

    end subroutine load_jasz

end program convertfort10

function tracematloc(n, over, ldo, jasmat, ldj, orbcost)
    implicit none
    integer i, j, k, n, ldo, ldj
    real*8 tracematloc, over(ldo, *), jasmat(ldj, *), vj, vk
    logical orbcost(*)
    tracematloc = 0.d0
    do j = 1, n
        if (.not. orbcost(j)) then
            do k = 1, n
                if (.not. orbcost(k)) then
                    vj = 0.d0
                    vk = 0.d0
                    do i = 1, n
                        if (orbcost(i)) then
                            vj = vj + jasmat(i, j)
                            vk = vk + jasmat(i, k)
                        end if
                    end do
                    tracematloc = tracematloc + vj*vk*over(j, k)
                end if
            end do
        end if
    end do
    return
end

function tracematnloc(n, over, ldo, orbcost)
    implicit none
    integer i, j, k, n, ldo
    real*8 tracematnloc, over(ldo, *)
    logical orbcost(*)
    tracematnloc = 0.d0
    do j = 1, n
        if (.not. orbcost(j)) then
            do k = 1, n
                if (.not. orbcost(k)) tracematnloc = tracematnloc + over(j, k)*over(k, j)
            end do
        end if
    end do
    return
end

subroutine findminmax(nnozero, nozero, nelorb, imin, imax, len)
    implicit none
    integer nnozero, nozero(*), nelorb, imin, imax, i, ix, iy, len
    imin = nelorb
    imax = 0
    do i = 1, nnozero
        iy = (nozero(i) - 1)/nelorb + 1
        ix = nozero(i) - (iy - 1)*nelorb
        if (ix .gt. imax .and. imax .le. nelorb) imax = ix
        if (iy .gt. imax .and. iy .le. nelorb) imax = iy
        if (ix .lt. imin .and. imin .ge. 1) imin = ix
        if (iy .lt. imin .and. imin .ge. 1) imin = iy
    end do
    len = imax - imin + 1
    return
end subroutine findminmax

subroutine dsygv_my(ipc, n, a, lda, over, ldo, doover, umat, eigmat&
        &, eig, mine, work, epsdgel, effham, info, nprocr, rank, comm_mpi)
    use constants, only: zone, zzero
    implicit none
    integer n, ipc, lda, ldo, mine, lwork, info, i, j, ndim, rank, nprocr, nproc&
            &, dimorb, ierr, comm_mpi
    real*8 a(ipc*lda, *), umat(ipc*n, n), eigmat(n), work(*), eig(n)
    real*8 over(ipc*ldo, *)
    real*8 epsdgel
    logical doover, effham
    !     This is standard generalized eigenvalue equation for an eignevector v with eigenvalue e:
    !      a v = e s v  effham=.true.
    !      a s v = e v     effham=.false.
    !      In input over is given already in diagonal form over = umat eig umat^dag
    !      where eig = eigmat^{-2} if |eig(i)/eig_max| > epsdgel  eig(i)=0 otherwise (and also eigmat(i)=0)
    !      so that the corresponding direction is disregarded from the calculation
    !     change basis for the
#ifdef PARALLEL
    include 'mpif.h'
#endif
    if (n .le. 1000) then
        nproc = 1 ! scalar code
    else
        nproc = nprocr
    end if
    if (doover) then
        call invsymeps(ipc, n, over, ldo, info, epsdgel, mine, umat, eigmat, nproc, rank, comm_mpi)
    end if
    !     if(rank.eq.0.or.nproc.gt.1) then
    if (ipc .eq. 1) then
        call dgemm_my('N', 'N', n, n, n, 1.d0, a, lda, umat, n, 0.d0, work, n, nproc, rank, comm_mpi)
        call dgemm_my('T', 'N', n, n, n, 1.d0, umat, n, work, n, 0.d0, a, lda, nproc, rank, comm_mpi)
    else
        call zgemm_my('N', 'N', n, n, n, zone, a, lda, umat, n, zzero, work, n, nproc, rank, comm_mpi)
        call zgemm_my('C', 'N', n, n, n, zone, umat, n, work, n, zzero, a, lda, nproc, rank, comm_mpi)
    end if
    if (effham) then
        do i = 1, n
            do j = 1, n
                a(ipc*(i - 1) + 1:ipc*i, j) = eigmat(i)*eigmat(j)*a(ipc*(i - 1) + 1:ipc*i, j)
            end do
        end do
    else
        do i = 1, n
            do j = 1, n
                if (eigmat(i) .ne. 0.d0 .and. eigmat(j) .ne. 0.d0) then
                    a(ipc*(i - 1) + 1:ipc*i, j) = a(ipc*(i - 1) + 1:ipc*i, j)/eigmat(i)/eigmat(j)
                else
                    a(ipc*(i - 1) + 1:ipc*i, j) = 0.d0
                end if
            end do
        end do
    end if

    do i = 1, mine - 1
        eig(i) = 0.d0
    end do

    ndim = n - mine + 1
    if (ipc .eq. 1) then
        call dsyev_my('V', 'U', ndim, a(mine, mine), lda, eig(mine), info, nproc, rank, comm_mpi)
    else
        call zsyev_my('V', 'U', ndim, a(2*mine - 1, mine), lda, eig(mine), info, nproc, rank, comm_mpi)
    end if

    !     write(6,*) ' Output eigenvalues '
    !     do i=1,n
    !     write(6,*) i,eig(i)
    !     enddo

    !     build the new eigenvector in the original basis
    do j = 1, ndim
        do i = 1, ndim
            work(ipc*ndim*(j - 1) + ipc*(i - 1) + 1:ipc*ndim*(j - 1) + ipc*i) = &
                    &a(ipc*(i + mine - 2) + 1:ipc*(i + mine - 1), j + mine - 1)*eigmat(i + mine - 1)
        end do
    end do

    if (ipc .eq. 1) then
        call dgemm_my('N', 'N', n, ndim, ndim, 1.d0, umat(1, mine), n                &
                &, work, ndim, 0.d0, a(1, mine), lda, nproc, rank, comm_mpi)
    else
        call zgemm_my('N', 'N', n, ndim, ndim, zone, umat(1, mine), n                &
                &, work, ndim, zzero, a(1, mine), lda, nproc, rank, comm_mpi)
    end if

    !     endif
#ifdef  PARALLEL
!        bcast just the relevant info to the nodes
    dimorb = ipc*(lda*(n - 1) + n)
    call bcast_real(a, dimorb, 0, COMM_MPI)
    call bcast_real(eig, n, 0, COMM_MPI)
#endif
    return
end subroutine dsygv_my

subroutine setorbcost(nelorbj_c, jasmat_c, orbcostn, nozeroj_c, nnozeroj_c)
    implicit none
    integer ix, iy, ixt, iyt, k, nelorbj_c, nnozeroj_c, indsto, indstos, indadd&
            &, indadds, nozeroj_c(*), icost
    real*8 jasmat_c(*)
    logical orbcostn(*)
    !      This subroutine try to fix the possibility that there exist more
    !      than one constant orbital in fort.10. In that case only the last
    !      constant orbital is used (address icost) and all the matrix
    !      rewritten.
    icost = 0
    do ixt = 1, nelorbj_c
        if (orbcostn(ixt)) icost = ixt
    end do
    if (icost .eq. 0) return
    ix = icost
    do ixt = 1, nelorbj_c
        if (orbcostn(ixt) .and. ixt .ne. icost) then
            do iyt = 1, nelorbj_c
                if (jasmat_c(nelorbj_c*(iyt - 1) + ixt) .ne. 0.d0) then
                    if (orbcostn(iyt)) then
                        iy = icost
                    else
                        iy = iyt
                    end if
                    indadd = nelorbj_c*(iy - 1) + ix
                    indadds = nelorbj_c*(ix - 1) + iy
                    indsto = nelorbj_c*(iyt - 1) + ixt
                    indstos = nelorbj_c*(ixt - 1) + iyt
                    jasmat_c(indadd) = jasmat_c(indadd) + jasmat_c(indsto)
                    if (ix .ne. iy) then
                        jasmat_c(indadds) = jasmat_c(indadd)
                    elseif (indsto .ne. indstos) then
                        jasmat_c(indadd) = jasmat_c(indadd) + jasmat_c(indstos)
                    end if
                    jasmat_c(indsto) = 0.d0
                    jasmat_c(indstos) = 0.d0
                end if
            end do
        elseif (ixt .eq. icost) then
            do iyt = 1, nelorbj_c
                if (jasmat_c(nelorbj_c*(iyt - 1) + ixt) .ne. 0.d0 .and. orbcostn(iyt) .and. iyt .ne. icost) then
                    iy = icost
                    indadd = nelorbj_c*(iy - 1) + ix
                    indsto = nelorbj_c*(iyt - 1) + ixt
                    indstos = nelorbj_c*(ixt - 1) + iyt
                    jasmat_c(indadd) = jasmat_c(indadd) + jasmat_c(indsto)
                    if (indsto .ne. indstos) jasmat_c(indadd) = jasmat_c(indadd) + jasmat_c(indstos)
                    jasmat_c(indsto) = 0.d0
                    jasmat_c(indstos) = 0.d0
                end if
            end do
        end if
    end do
end subroutine setorbcost
subroutine checkmat_complex(nelorbj_c, jasmat_c, lead, nozeroj_c, nnozeroj_c, checkall, rank, epsr, symmagp, ipf)
    use allio, only: yes_hermite, pfaffup, kiontot
    implicit none
    integer ix, iy, ixt, iyt, k, ndim, nelorbj_c, lead, nnozeroj_c, rank&
            &, indadds, imax, jmax, ipf, ndimh, nozeroj_c(*)
    complex*16 jasmat_c(*)
    real*8 eps, epsr, error, errormax, cost, scalem
    integer count
    logical check, checkall, symmagp
    real*8, external :: dlamch
    complex*16, dimension(:), allocatable :: jasmat_sav
    complex*16 value
    checkall = .true.
    allocate (jasmat_sav(lead*nelorbj_c))
    ndim = nelorbj_c*lead
    jasmat_sav = 0.d0
    jasmat_sav(1:ndim) = jasmat_c(1:ndim)
    jasmat_c(1:ndim) = 0.d0
    scalem = 0.d0
    do k = 1, nnozeroj_c
        iy = (nozeroj_c(k) - 1)/nelorbj_c + 1
        ix = nozeroj_c(k) - (iy - 1)*nelorbj_c
        if (iy .le. nelorbj_c) then
            jasmat_c(lead*(iy - 1) + ix) = jasmat_sav(lead*(iy - 1) + ix)
            if (ipf .eq. 2) then
                jasmat_c(lead*(ix - 1) + iy) = -jasmat_sav(lead*(iy - 1) + ix)
                ndimh = nelorbj_c/2
                if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                    value = jasmat_sav(lead*(iy - 1) + ix)
                    if (yes_hermite) then
                        if (iy .le. ndimh .and. ix .le. ndimh .and. .not. pfaffup) then
                            jasmat_c(lead*(iy + ndimh - 1) + ix + ndimh) = conjg(value)
                            jasmat_c(lead*(ix + ndimh - 1) + iy + ndimh) = -conjg(value)
                        elseif (iy .gt. ndimh .and. ix .le. ndimh) then
                            jasmat_c(lead*(ix + ndimh - 1) + iy - ndimh) = conjg(value)
                            jasmat_c(lead*(iy - ndimh - 1) + ix + ndimh) = -conjg(value)
                        end if
                    else
                        if (iy .le. ndimh .and. ix .le. ndimh .and. .not. pfaffup) then
                            jasmat_c(lead*(iy + ndimh - 1) + ix + ndimh) = value
                            jasmat_c(lead*(ix + ndimh - 1) + iy + ndimh) = -value
                        elseif (iy .gt. ndimh .and. iy .le. nelorbj_c .and. ix .le. ndimh) then
                            jasmat_c(lead*(ix + ndimh - 1) + iy - ndimh) = value
                            jasmat_c(lead*(iy - ndimh - 1) + ix + ndimh) = -value
                        end if
                    end if
                end if
            else
                if (symmagp) then
                    if (yes_hermite) then
                        jasmat_c(lead*(ix - 1) + iy) = conjg(jasmat_sav(lead*(iy - 1) + ix))
                    else
                        jasmat_c(lead*(ix - 1) + iy) = jasmat_sav(lead*(iy - 1) + ix)
                    end if
                end if
            end if
            cost = abs(jasmat_c(lead*(iy - 1) + ix))
            if (cost .gt. scalem) scalem = cost
        end if
    end do
    eps = epsr
    if (eps .eq. 0.d0) eps = dlamch('E')

    error = 0.d0
    count = 0
    imax = 0
    jmax = 0
    errormax = 0
    do ixt = 1, nelorbj_c
        do iyt = 1, nelorbj_c
            cost = abs(jasmat_sav(lead*(iyt - 1) + ixt))
            if (cost .gt. eps*scalem) then
                check = .false.
                if (abs(jasmat_c(lead*(iyt - 1) + ixt)) .ne. 0.d0 .or.&
                        &((symmagp .or. ipf .eq. 2) .and. abs(jasmat_c(lead*(ixt - 1) + iyt)) .ne. 0.d0)) check = .true.
                if (.not. check) then
                    checkall = .false.
                    count = count + 1
                    error = error + cost
                    if (cost .gt. errormax) then
                        errormax = cost
                        imax = ixt
                        jmax = iyt
                    end if
                    if (rank .eq. 0 .and. abs(jasmat_sav(lead*(iyt - 1) + ixt)) .gt. 1d-7) &
                            &write (6, *) ' Not found element ', ixt, iyt&
                            &, jasmat_sav(lead*(iyt - 1) + ixt)
                end if
            end if
        end do
    end do
    if (rank .eq. 0 .and. .not. checkall) then
        write (6, *) ' Number of wrong matrix elements x~2 =', count
        write (6, *) ' Average error =', error/count/scalem
        write (6, *) ' Max  error =', errormax/scalem, imax, jmax
    end if

    deallocate (jasmat_sav)
end subroutine checkmat_complex
subroutine checkmat(nelorbj_c, jasmat_c, lead, nozeroj_c, nnozeroj_c, checkall, rank, epsr, symmagp, ipf)
    use allio, only: pfaffup, kiontot
    implicit none
    integer ix, iy, ixt, iyt, k, ipf, ndim, ndimh, nelorbj_c, lead, nnozeroj_c, rank&
            &, indadds, imax, jmax, nozeroj_c(*)
    real*8 jasmat_c(*), eps, epsr, error, errormax, cost, scalem, value
    integer count
    logical check, checkall, symmagp
    real*8, external :: dlamch
    real*8, dimension(:), allocatable :: jasmat_sav
    checkall = .true.
    allocate (jasmat_sav(lead*nelorbj_c))
    ndim = nelorbj_c*lead
    jasmat_sav = 0.d0
    jasmat_sav(1:ndim) = jasmat_c(1:ndim)
    jasmat_c(1:ndim) = 0.d0
    scalem = 0.d0
    do k = 1, nnozeroj_c
        iy = (nozeroj_c(k) - 1)/nelorbj_c + 1
        ix = nozeroj_c(k) - (iy - 1)*nelorbj_c
        if (iy .le. nelorbj_c) then
            jasmat_c(lead*(iy - 1) + ix) = jasmat_sav(lead*(iy - 1) + ix)
            if (symmagp .and. ipf .eq. 1) jasmat_c(lead*(ix - 1) + iy) = jasmat_sav(lead*(iy - 1) + ix)
            if (ipf .eq. 2) then
                jasmat_c(lead*(ix - 1) + iy) = -jasmat_sav(lead*(iy - 1) + ix)
                ndimh = nelorbj_c/2
                if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                    value = jasmat_sav(lead*(iy - 1) + ix)
                    if (iy .le. ndimh .and. ix .le. ndimh .and. .not. pfaffup) then
                        jasmat_c(lead*(iy + ndimh - 1) + ix + ndimh) = value
                        jasmat_c(lead*(ix + ndimh - 1) + iy + ndimh) = -value
                    elseif (iy .gt. ndimh .and. ix .le. ndimh) then
                        jasmat_c(lead*(ix + ndimh - 1) + iy - ndimh) = value
                        jasmat_c(lead*(iy - ndimh - 1) + ix + ndimh) = -value
                    end if
                end if
            end if
            cost = abs(jasmat_c(lead*(iy - 1) + ix))
            if (cost .gt. scalem) scalem = cost
        end if
    end do
    eps = epsr
    if (eps .eq. 0.d0) eps = dlamch('E')

    error = 0.d0
    count = 0
    imax = 0
    jmax = 0
    errormax = 0
    do ixt = 1, nelorbj_c
        do iyt = 1, nelorbj_c
            cost = dabs(jasmat_sav(lead*(iyt - 1) + ixt))
            if (cost .gt. eps*scalem) then
                check = .false.
                if (jasmat_c(lead*(iyt - 1) + ixt) .ne. 0.d0 .or.&
                        &((symmagp .or. ipf .eq. 2) .and. jasmat_c(lead*(ixt - 1) + iyt) .ne. 0.d0)) check = .true.
                if (.not. check) then
                    checkall = .false.
                    count = count + 1
                    error = error + cost
                    if (cost .gt. errormax) then
                        errormax = cost
                        imax = ixt
                        jmax = iyt
                    end if
                    if (rank .eq. 0 .and. abs(jasmat_sav(lead*(iyt - 1) + ixt)) .gt. 1.d-7) &
                            &write (6, *) ' Not found element ', ixt, iyt, jasmat_sav(lead*(iyt - 1) + ixt)
                end if
            end if
        end do
    end do
    if (rank .eq. 0 .and. .not. checkall) then
        write (6, *) ' Number of wrong matrix elements x~2 =', count
        write (6, *) ' Average error =', error/count/scalem
        write (6, *) ' Max  error =', errormax/scalem, imax, jmax
    end if

    deallocate (jasmat_sav)
end subroutine checkmat

subroutine purify(jasmat, nelorbjh, overj, nelorbj, orbcostn)
    implicit none
    integer nelorbj, nelorbjh, icost, i, j
    real*8 jasmat(nelorbjh, *), overj(nelorbj, *), gam, alf
    logical orbcostn(nelorbjh)
    !      orthogonalization with respect to the constant orbital.
    icost = 0
    do i = 1, nelorbjh
        if (orbcostn(i)) icost = i
    end do
    if (icost .eq. 0) return
    gam = 0.d0
    do i = 1, nelorbjh
        do j = 1, nelorbjh
            gam = gam + overj(i, icost)*overj(j, icost)*jasmat(i, j)
        end do
    end do
    alf = gam/overj(icost, icost)**2
    jasmat(icost, icost) = jasmat(icost, icost) - alf
end subroutine purify
