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

!#include "mathvec.h"
program main
    use allio
    use convertmod
    use IO_m
    use trexio
    ! by E. Coccia (8/11/10)
    use extpot
    ! by E. Coccia (28/12/10)
    use van_der_waals
    ! by E. Coccia (27/5/11)
    use tot_angle
    use scal_lins
#ifdef PARALLEL
    ! for mpiio by Y. Luo (24/2/15)
    use mpiio, only: mpiio_file_reset_view, mpiio_file_set_zero          &
           &, mpiio_file_get_disp, mpiio_file_create_view
#endif

    real*8 timeppp, drand1, inittime, enercont, mapping, enthalpy        &
           &, weight_vir, costwnn, costexpn, coeff_nw, enerold           &
           &, costfact, ttrysav, ratiowconfn, wconfnsav, ttryt           &
           &, enercuto, psioldinp, voffpseudo, srforce, wforce           &
           &, srforcew, wforcew, costmpi(3), costsr, identity            &
           &, timemcp, timeoptp, time1p, vr, vl, celldmsav, costprsav    &
           &, celllogb(3), cellelb(3), dist_min, eloc_old, eloc_new      &
           &, psi_old, psi_new, energyq, ekinq, dt4, ekinqp, rs_write    &
           &, celldm_write(3), scale_spsi, countav_all, countt_all       &
           &, psirav_all, avratio, nmovet, minratio, maxratio            &
           &, target_ratio, avenernum, avenerden, avener, wbraw          &
           &, cellscale_write(3), cclock_zero, reweight_dmc, enerdiff    &
           &, sr2elb(9), sr2logb(9), gpu_size, gpu_sizej
    integer nleft, i, j, k, ind, ii, jj, kk, kkf, kkfg, imin, imax, imaxn&
           &, stodim, iese_eff, np3p3, maxrank, maxdimeig, nmollimit     &
           &, firstmolt, iflagpip, ind_min, ieskinion, new_threads       &
           &, nprogress, dimfk, perbin, jcol, irow, comm_mpi, nweightu   &
           &, np3p4, np3p5, np3p6, np3p7, ifreqchange, indtju            &
           &, ifreqchanger, np3p8, np3p9, np3p10, np3p11, nwnk, nwnkp    &
           &, commrep_mpi_sav, rankrep_sav, nmol_count

    !       integer rank_av,nw_av,comm_av

    logical noproj, flagmu, flagcont, forceyes                           &
           &, flag, iesmeas, dir_exist, yesprojm                         &
           &, yeswrite, yeswrite12, acc_dyn, yescomm, someparameter      &
           &, yesforce, yesdodet, yesdodet_nof, someparameterdet         &
           &, yesdo_imag, yesupdate_ion, wdone, yes_deallocate

    integer, external              :: wf_sign
    real(8), external              :: cclock, dlamch
    integer, external              :: omp_get_max_threads

    real(8), dimension(:, :), allocatable :: kelsav, rion_write, rion_sav
    real(8), dimension(:, :, :), allocatable :: rionall
    real(8), dimension(:), allocatable :: detmat_sav, work, zagp_imag
    real(8), dimension(:), allocatable :: wbra_t, ener_t, wconfn_kps
    real(8), dimension(:, :), allocatable :: etot_t, wtotf_t
    integer, dimension(:), allocatable :: jbrasymiesup
    character(lchlen)              :: path, scratchpath
    character(lchlen + 20)         :: ranseedfilename

    !     memory required for Algorithmic Differentiation (AAD)

    integer                        :: nelorbjmax, indtmax, nshelljmax    &
           &, iscramaxold, isdistp, iscramaxb, indshell, prep_try
    real*8 totforce, norm_tab, totpulay, elocb(2), logpsib(2), eloc(2)   &
           &, logpsi(2)
    real*8, dimension(:, :), allocatable :: kelind
    real*8, allocatable            :: kelb(:, :), kelindb(:, :)          &
           &, rionb(:, :), tabpipb(:), winvupb(:), winvdob(:), ainvb(:)  &
           &, ainvupbb(:), ainvdobb(:), psipb(:), distb(:), rb(:)        &
           &, rmub(:), iond_cartb(:), winvb(:), winvjb(:), winvbarb(:)   &
           &, winvjbarszb(:), winvjbarb(:), prefactorb(:), wpseudob(:)   &
           &, legendreb(:), rmucosb(:), rmusinb(:), tmub(:), winvfnn(:)  &
           &, winvjfn(:), ivicb(:), tabpipsav(:), kelindlb(:, :), vjb(:) &
           &, duprb(:), vjurb(:), rionlb(:, :), vjlb(:), duprlb(:)       &
           &, vjurlb(:), detmatb(:), mu_cb(:), jasmatb(:), jasmatszb(:)  &
           &, detmatlb(:), mu_clb(:), jasmatlb(:), jasmatszlb(:)         &
           &, forcedw(:, :), pulaydw(:, :), fbead(:, :), mass_ion(:, :)  &
           &, detmat_cb(:), muj_cb(:), jasmat_cb(:), jasmatsz_cb(:)

#ifdef _OFFLOAD
    integer omp_get_num_devices
#endif
#ifdef _CUSOLVER
    integer                        :: stat
#endif
#ifdef PARALLEL
    include 'mpif.h'
    integer skip, status(MPI_STATUS_SIZE), ithread
    real(8) naccmpi, nmovetmpi, naccpseudompi, acclargempi
    real(4), allocatable           :: SP_buffer(:, :)
    real(8), allocatable           :: DP_buffer(:, :)
    integer                        :: buffer_counter, buffer_depth       &
           &, kelcont_size, SP_block_size, DP_block_size
#else
    integer status(1), nn
    integer                        :: iargc
    character(20) str
    character(100) name_tool
    call getarg(1, str)
    nn = iargc()

    if (nn .gt. 0) then
        select case (str)
            case ('help', '-help', '--help')
                write (*, *) 'Choose among the following cases'
                write (*, *) 'vmc      Input file for VMC'
                write (*, *) 'dmc      Input file for DMC'
                write (*, *) 'lrdmc    Input file for lattice reguralized DMC'
                write (*, *) 'opt      Input file for the optimization'
                write (*, *) 'optmol   Input file for the optimization with'// &
                   & ' molecular orbitals'
                write (*, *) 'dyn      Input file for the dynamics'
                write (*, *) 'quantum  Input file for the quantum dynamics'
                write (*, *) 'test     Input file for testing TurboRVB'
            case ('vmc')
                write (*, *) ' Sample file for the Variational Monte Carlo '
                name_tool = 'datasvmc'
                call help_online(name_tool)
            case ('dmc')
                write (*, *) ' Sample file for the DMC '
                name_tool = 'datasdmc'
                call help_online(name_tool)
            case ('lrdmc')
                write (*, *) ' Sample file for the LRDMC '
                name_tool = 'datasfn'
                call help_online(name_tool)
            case ('opt')
                write (*, *) ' Sample file for the optimization'
                name_tool = 'datasmin'
                call help_online(name_tool)
            case ('optmol')
                write (*, *) ' Sample file for the optimization with'// &
                   & ' molecular orbitals'
                name_tool = 'datasminmol'
                call help_online(name_tool)
            case ('dyn')
                write (*, *) ' Sample file for the molecular dynamics'
                name_tool = 'datasdyn'
                call help_online(name_tool)
            case ('quantum')
                write (*, *) ' Sample file for the quantum molecular dynamics'
                name_tool = 'datasquantum'
                call help_online(name_tool)
            case ('test')
                write (*, *) ' Sample file for testing TurboRVB'
                name_tool = 'datastest'
                call help_online(name_tool)
                case default
                write (*, *) ' help not available'
        end select
        stop

    end if
#endif

#ifdef PARALLEL
    call mpi_init(ierr)
    ! call mpi_init_thread(MPI_THREAD_FUNNELED,ithread,ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#else
    rank = 0
    nproc = 1
#endif

#ifdef _OFFLOAD
    if (rank .eq. 0) write (6, *) ' #  GPU used/1mpi process=:'&
        &, omp_get_num_devices()
#endif

    inittime = cclock()
#ifdef PARALLEL
    ! all have to be syncronized with the master
    call mpi_bcast(inittime, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
    !     output version information
    if (rank .eq. 0) call print_version

    if (nproc .gt. 1 .and. rank .eq. 0) write (6, *)&
       & ' Number of mpi proc =', nproc
#ifdef _OPENMP
    old_threads = omp_get_max_threads()
#else
    old_threads = 1
#endif
    if (rank .eq. 0) write (6, *) ' Number of threads/mpi proc =', old_threads
    new_threads = old_threads
#ifdef UNREL_SMP
#ifdef PARALLEL
    call mpi_bcast(old_threads, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
    call omp_set_num_threads(old_threads) ! force the same.
    if (new_threads .ne. old_threads) write (6, *) &
         &' Warning error in number of threads (recovered) in proc # '&
         &, rank, old_threads, new_threads
    if (rank .eq. 0) write (6, *) ' Warning init. value of threads/mpi task'&
         &, old_threads
#endif
    timepp = cclock()
    time_main = 0.d0
    ! by E. Coccia (8/11/10): call read_cube into Initializeall and
    ! call spline interpolation if ext_pot=.true.
    ! by E. Coccia (4/2/11): read vdw.dat if the van der Waals term
    ! in the external potential is present
    call Initializeall
    yes_ontarget = .false.
  if(rank.eq.0) write(6,*) ' Size arrays ', size(jasmat), size(muj_c)&
      &, size(jasmat_c), size(detmat), size(detmat_c), size(projm), size(mu_c)&
      &, size(eagp_pfaff), size(psip), size(ainvs), size(winvs), size(winvsj)&
      &, size(winv), size(winvj), size(psinew)+size(agp), size(agpn)&
      &, size(ainv), size(winvbar), size(winvjbar), size(winvfn)&
      &, size(winvbarfn),size(ainvdo),size(ainvup)

#ifdef _OFFLOAD
    if (itest .eq. 2) then
        cost = 1.d0
    else
        !       In the dmc it is more efficient
        cost = indt
    end if
    if (ipc*nelorbh*nelup_mat*cost .gt. max_target) then
        yes_ontarget = .true.
        !         if(.not.yes_fastbranch) then
        !         yes_fastbranch=.true.
        !       if(rank.eq.0.and.itest.ne.2) &
        !      & write(6,*) ' Warning Blas2 in target works only with &
        !      & yes_fastbranch option on --> forced to .true. '
        !         endif
    else
        yes_ontarget = .false.
    end if
    if (iessz) then
        write (6, *) ' Warning  GPU is not optimized  with Jsz, use the more general e.g. -27 '
    end if
    gpu_size = size(jasmat) + size(muj_c) + size(jasmat_c) + size(detmat) + size(detmat_c)
    gpu_size = gpu_size + size(projm) + size(mu_c) + size(eagp_pfaff) + size(psip)
    gpu_size = gpu_size + size(ainvs) + size(winvs) + size(winvsj) + size(winv) + size(winvj)
    gpu_size = gpu_size + size(psinew) + size(agp) + size(agpn) + size(ainv) + size(winvbar)
    gpu_size = gpu_size + size(winvjbar) + size(winvfn) + size(winvbarfn) + size(ainvdo)
    gpu_size = gpu_size + size(ainvup)
    gpu_sizej = size(jasmat) + size(muj_c) + size(jasmat_c)
    gpu_sizej = gpu_sizej + size(winvsj) + size(winvj) + size(winvjbar)

    if (rank .eq. 0) write (6, *)&
       & ' Memory allocated in  the GPU (Gb) ALL/Jastwow =  '&
       &, gpu_size*8.d-9, gpu_sizej*8d-9

    gpu_sizej = size(jasmat) + size(muj_c) + size(jasmat_c)
    gpu_size = gpu_sizej + size(detmat) + size(detmat_c) + size(projm) + size(mu_c)
    if (rank .eq. 0) write (6, *)&
       & ' Memory common in the GPU All/Jastrow='&
       &, gpu_size*8d-9, gpu_sizej*8d-9

! matrix common to all walkers and unchanged during the Markov chain
! May be changed in optimization of dynamics:
! jasmat
! muj_c
! jasmat_c
! detmat
! detmat_c
! projm
! mu_c
! eagp_pfaff
!$omp target data map(to: jasmat,muj_c,jasmat_c,detmat,detmat_c&
!$omp &,projm,mu_c,eagp_pfaff)&
!$omp &  map(to:psip&
!$omp &,ainvs,winvs,winvsj,winv,winvj,psinew&
!$omp &,agp,agpn,ainv,winvbar,winvjbar,winvfn,winvbarfn,ainvup,ainvdo)
    if (rank .eq. 0) then
        if (yes_ontarget) then
            write (6, *) ' Warning Blas2 in target  '
        else
            if (rank .eq. 0) write (6, *) ' Warning Blas2 in cpu '
        end if
    end if

#ifdef _CUSOLVER
    if (ipf .ne. 2) then
        if (rank .eq. 0) write (6, *) ' Warning, using cusolver routines '
        ldworkspace = 1
        lzworkspace = 1
        !
        if (ipc .eq. 1) then
#ifdef RISC
            call cusolver_dgetrf_buffersize_&
                &(handle, stat, nelup, nelup, psip, nelup, ldworkspace)
#else
            call cusolver_dgetrf_buffersize&
                &(handle, stat, nelup, nelup, psip, nelup, ldworkspace)
#endif
            allocate (dev_dgetri_workspace(nelup, nelup))
            allocate (dev_zgetri_workspace(1, 1))
        else
#ifdef RISC
            call cusolver_zgetrf_buffersize_&
                &(handle, stat, nelup, nelup, psip, nelup, lzworkspace)
#else
            call cusolver_zgetrf_buffersize&
                &(handle, stat, nelup, nelup, psip, nelup, lzworkspace)
#endif
            allocate (dev_dgetri_workspace(1, 1))
            allocate (dev_zgetri_workspace(nelup, nelup))
        end if
        call checkiflagerr&
            &(stat, rank, ' Something went wrong in calculating GPU buffer space')
        !
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
#endif

    !======================================================================
    !======================================================================
    !  This is the main loop. IMPORTANT the main index of the code is I
    !  Any subroutine after ''contains'' has not to use the index I inside any loop.
    !  The second important index is J, it runs over the walkers of each processor.
    !  JS (=J for serial code) instead is the walker index for all processors.
    !  The other indices ii,jj,..., are dummy indices.

    ! added test for maxtime if no optimization is performed (as in this case one
    ! has to close the cycle) itestr.eq.-5
    cclock_zero = cclock() - inittime
    ! all have to be syncronized with the master
#ifdef PARALLEL
    call mpi_bcast(cclock_zero, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

    !PRINT *, "Continue 1"

    do while (i_main .lt. ngen + iend .and. (cclock_zero .lt. maxtime .or. itestr .eq. -5))

        i_main = i_main + 1 ! MC generations

        if (i_main .eq. iend + 1) then
            ! compute by scratch the correcting factor
            if (nfat .ne. 0) then
                wcort = wcorw(1)
                do jj = 2, nfat
                    wcort = wcort*wcorw(jj)
                end do
            end if

            if (iopt .eq. 1) then
                pseudologic = pseudorandom
                iesrandoml = iesrandoma
            else
                pseudologic = .false.
                iesrandoml = .false.
            end if
            flagcont = .true.

            timepp = cclock()

            !
            call makeallmeas

            !
            time_meas = time_meas + cclock() - timepp
            ! restore logic
            pseudologic = pseudorandom
            iesrandoml = iesrandoma

            ! write(6,*) ' after makeallmeas ',(indpar_tab(ii),ii=1,nshell)

            call checkiflagerr(iflagerr, rank, 'Something went wrong in QMC initialization!')

#ifdef  _TIME
            timings = 0.d0
            timingsb = 0.d0
#endif
            if (rank .eq. 0) write (6, *) ' Initialization OK '

            !       if(rank.eq.0) write(6,*) ' Initial projm =',sum(abs(projm(1:nelorbh*nelorb_c)))
            !       if(rank.eq.0) write(6,*) ' Initial mu_c =',sum(abs(mu_c(:,1:nelorb_c)))

            if (itestr .eq. -5) then

                call preprpar(rpar, kp_ion, iond, nion, kiontotj         &
                       &, nozeroj_c, nnozeroj_c, nelorbj_c, jbraj        &
                       &, iesfree, indfree, jbrajsz, iesinv, indinv      &
                       &, orbcostn, kiontot, nozero_c, nnozero_c         &
                       &, nelorb_c, jbradet, iessw, indsw, adrlambda     &
                       &, whereiesm, iesm, whereiesup, iesup             &
                       &, iond_cart, typeorb, nshellj_c, multj_c         &
                       &, ioccj_c, occj_c, jas_invariant, orbps)

                ! Untouched parameters below a certain treshold
                do ii = 1, kp_ion
                    if ((rpar(ii) .lt. rmaxj .or. rmaxj .eq. 0)                           &
                         &.and. (ii .gt. endinv .and. rpar(ii) .gt. 0)) rpar(ii) = 0.d0
                    if ((rpar(ii) .lt. rmaxinv .or. rmaxinv .eq. 0)                       &
                         &.and. ii .le. endinv) rpar(ii) = 0.d0
                    if ((abs(rpar(ii)) .lt. rmax .or. rmax .eq. 0) .and. rpar(ii) .lt. 0.d0)    &
                         &   rpar(ii) = 0.d0
                end do

                if (rank .eq. 0 .and. ncg_adr .gt. 0) then
                    if (npar .gt. 0) write (6, *)&
                       & ' Parametrization charge-Jastrow ', npar
                    if (nparinv .gt. 0) write (6, *)&
                       & ' Parametrization spin-Jastrow ', nparinv
                    if (nparsw .gt. 0) write (6, *)&
                       & ' Parametrization AGP matrix ', nparsw
                end if
#ifdef PARALLEL
                if (yesquantum .and. commrep_mpi .ne. commsr_mpi) then
                    do ii = 1, kp_ion
                        psip(kp_ion + ii) = abs(rpar(ii))
                    end do
                    ! If there is some accepted parameter in some bead accept the same for all
                    call mpi_allreduce&
                        &(psip(kp_ion + 1), psip, kp_ion, MPI_DOUBLE_PRECISION, MPI_MIN, commsr_mpi, ierr)
                    do ii = 1, kp_ion
                        if (psip(ii) .eq. 0.d0) rpar(ii) = 0.d0
                    end do
                end if
#ifdef UNREL_DIAG
                ! set consistency among processors
                call bcast_real(rpar, kp_ion, 0, commsr_mpi)
#endif
#endif

                if (ncg_adr .gt. 0) then
                    call pareff(npar, initpar, nparsw, initparsw, nparinv, initparinv&
                         &, endinv, ncg_adr, kp0, rpar, reducel, jas_invariant, adrlambda, nmax_ion&
                         &, type_atom, allfit, orbps)
                    call project_v(.true.)
                    call project_alphavar
                    !     elseif(rmaxj.ne.0.or.rmaxinv.ne.0.or.rmax.ne.0.and.ireadmin.eq.1) then
                    !    no parametrization but locality
                    !     if(rank.eq.0) write(6,*) ' Warning setting to zero unoptimazible > rmax '
                    !     call project_rmax
                    !     if(contraction.ne.0) &
                    !    & call dgemm_my('N','N',nelorbh,nelcol_c,nelorb_c,1.d0,mu_c,nelorbh&
                    !    &,detmat_c,nelorb_c,0.d0,projm,nelorbh,nprocu,rank)
                end if

            end if ! if itestr.eq.-5

            ! timescra=timescra+cclock()-timep excluded initialization time
            ! computing by scratch the factors involved

            timeinit = cclock() - timepp
            time = cclock() ! count the time without inizialization
            ! progress indicator
            nprogress = max(1, int(ngen*0.01))
            if (rank .eq. 0) write (6, *) 0, "% progress starts!"
            time1p = time

        end if ! end if for scratch i.eq.iend+1

        !**************** MAIN OPERATIONS **************

        timemcp = cclock()
        ener = 0.d0

        ! make the single electron move for all the walkers
        call makeqmcsingleel
        !--------------------------------------------------
        timeppp = cclock()
        timemc = timemc + timeppp - timemcp

#ifdef PARALLEL
        if (freqcheck .ne. 0) then
            if (mod(i_main, freqcheck) .eq. 0) then ! do the check each few steps
                call checkiflagerr(iflagerr, rank, 'Something wrong in QMC, sing Det or similar!')
            end if
        end if
#else
        call checkiflagerr(iflagerr, rank, 'Something wrong in QMC, sing Det or similar!')
#endif

        if ((itest .eq. 2 .or. (itest .eq. 1 .and. iesmeas))) then

            if (itest .eq. 2) then
                pseudologic = pseudorandom
                iesrandoml = iesrandoma
            else
                ! do not use extra random number
                iesrandoml = .false.
                pseudologic = .false.
                !          if(itest.eq.1.and..not.pseudorandom) pseudologic=.true.
            end if
            flagcont = .false.
            timepp = cclock()

            ! measures all quantities related to the simulation for all walkers separately
            time_main = time_main + timepp - timeppp
            call makeallmeas
            !--------------------------------------------------

            ! restore logic (for uptabtot)
            pseudologic = pseudorandom
            iesrandoml = iesrandoma
            timeppp = cclock()
            time_meas = time_meas + timeppp - timepp
        end if ! end if itest.eq.2

        ! calculation average quantities before branching
        if (itest .ne. 2) then
            do j = 1, in1
                counttot = counttot + 1.d0
                if (yescut(j)) then
                    countreg = countreg + 1.d0
                end if
                if (mod(typereg, 2) .eq. 0) then
                    psiav = psiav + psidetln(j)
                    psisav = psisav + psidetln(j)**2
                else
                    psiav = psiav + psiln(j)
                    psisav = psisav + psiln(j)**2
                end if
            end do
        end if
        !
        ! collect energy and weights after one generation
        !
        ener = 0.d0
        wtotf(1) = 0.d0
        do j = 1, in1
! Changed otherwise Pulay forces could be inconsistent when enert=/enertrue
!       ener = ener + enertrue(j) * factorsr(j)
            ener = ener + enert(1, j)*factorsr(j)
            wtotf(1) = wtotf(1) + factorsr(j)
        end do

        wtotf(2) = 0.d0
        do j = 1, in1
            wtotf(2) = wtotf(2) + wint(j) ! integrated weight for processor
        end do
        !
        ! computation of the total energy after MC generation
        ! it is assumed that the integrated weights acts only on iese <= 6 variables
        ! Perform the weighted integration for all in1 walkers per processor.
        !
        ! correlation functions to be averaged
        ! NB in the case of complex wave function the imaginary part is disregarded!
        call dgemv('T', in1, iese_eff, 1.d0, econf, in1, wint, 1, 0.d0, etot, 1)
        ! writing also the variance + writing O^k operators
        if (np3 - iese_eff .gt. 0) then
            call dgemv('T', in1, np3 - iese_eff, 1.d0, econf(iese_eff*in1 + 1), in1, &
                       factorsr, 1, 0.d0, etot(iese_eff + 1), 1)
        end if
        if (i_main .eq. ibinit .and. iopt .eq. 1) then
            avenernum = 0.d0
            avenerden = 0.d0
            tave_cyrus = 0.d0
            tcount_cyrus = 0.d0
        end if
        tave_cyrus = tave_cyrus + sum(t_cyrus(1:in1))
        tcount_cyrus = tcount_cyrus + in1
        avenernum = avenernum + ener
        avenerden = avenerden + wtotf(1)

#ifdef PARALLEL
        wtot(1) = ener
        wtot(2:3) = wtotf(1:2)
        wtot(4:np3p3) = etot(1:np3)
        wtot(np3p4) = nacc
        wtot(np3p5) = nmovet
        wtot(np3p6) = tave_cyrus
        wtot(np3p7) = tcount_cyrus
        wtot(np3p8) = countav
        wtot(np3p9) = countt
        wtot(np3p10) = avenernum
        wtot(np3p11) = avenerden
        if ((commrep_mpi .eq. mpi_comm_world) .or. decoupled_run) then
            call reduce_base_real_to(np3p11, wtot, psip, commrep_mpi, -1)
        else
            call reduce_base_real_to(np3p7, wtot, psip, MPI_COMM_WORLD, -1)
            call reduce_base_real_to(4, wtot(np3p8), psip(np3p8), commrep_mpi, -1)
        end if
        call dcopy(np3, psip(4), 1, etot, 1)
        ener = psip(1)
        wtotf(1) = psip(2)
        wtotf(2) = psip(3)
        avreweight = psip(np3p8)/psip(np3p9)
        avratio = psip(np3p4)/psip(np3p5)
        avener = psip(np3p10)/psip(np3p11)
#else
        avreweight = countav/countt
        avratio = nacc/nmovet
        avener = avenernum/avenerden
#endif
        if (changelambda) then
            if (yesnleft .and. (i_main .ge. ibinit .or. cclock_zero .ge. maxtime/20.d0&
                 &.or. iopt .eq. 0 .or. iopt .eq. 3)) lambda = -avener
            if (yesnleft .and. (i_main .ge. ibinit .or. cclock_zero .ge. maxtime/20.d0)&
                 &.and. (iopt .eq. 1 .or. iopt .eq. 2) .and. rank .eq. 0 .and. wdone) then
                write (6, *) ' Warning beginning to change trial energy=', i_main, avener*ris(2)
                wdone = .false.
            end if
        end if
        if (change_tstep&
           & .and. (avratio .lt. minratio .or. avratio .gt. maxratio)&
           & .and. tstep .gt. 1.d0&
           & .and. tstep .lt. 10.d0) then
            if ((i_main .gt. 1 .or. (iopt .eq. 0 .or. iopt .eq. 3))&
                & .and. mod(i_main, ifreqchanger) .eq. 0&
                & .and. i_main .ne. ngen) then
                tstep = tstep*avratio/target_ratio
                if (rank .eq. 0) write (6, *) ' Warning changing tstep/Acceptance ='&
                     &, i_main, real(tstep), real(avratio)
                nacc = 0.d0
                nmovet = 0.d0
            end if
        end if

        if ((avreweight .lt. 0.65d0 .or. avreweight .gt. .95d0)&
             &.and. epscuttype .ne. 0 .and. change_epscut .and.&
             &(itestr .ne. -5 .or. nweight .gt. 20)) then
            if ((i_main .gt. 1 .or. (iopt .eq. 0 .or. iopt .eq. 3)) .and. mod(i_main, ifreqchange) .eq. 0) then
                !   if(i_main.gt.1.or.iopt.ne.1) then
                if (avreweight .gt. 0.95d0 .and. avreweight .lt. .99d0) then
                    cost = .2d0/(1.d0 - avreweight)
                elseif (avreweight .ge. .99d0) then
                    cost = 20.d0
                else
                    cost = avreweight/0.8d0
                end if
                !     There may be convergence problems and therefore I put a cutoff below.
                if (cost .gt. 3d0) cost = 3.d0
                if (cost .lt. 0.3333333333d0) cost = 0.3333333333d0
                epscutu = epscutu*cost
                epscut = epscutu
                epstlu = epstlrat*epscutu
                epstl = epstlrat*epscut
                if (rank .eq. 0) write (6, *) ' Warning changing epscut/reweight ', i_main&
                     &, real(epscutu), real(avreweight)
                countt = 0.d0
                countav = 0.d0
                psirav = 0.d0
            end if
        end if

        !     Inizialization counters
        if (i_main .eq. 1 .and. iopt .eq. 1) then
            if (itest .eq. 2) then
                countt = 0.d0
                countav = 0.d0
                psirav = 0.d0
            end if
            nacc = 0.d0
            nmovet = 0.d0
            avenernum = 0.d0
            avenerden = 0.d0
            tave_cyrus = 0.d0
            tcount_cyrus = 0.d0
        end if

        ! Rescaling the unit energy.
        ! Put everything in Hartree and Atomic units
        if (rankrep .eq. 0) then

            cost = 1.d0/wtotf(1)
            ener = ener*cost
            call dscal(np3 - iese_eff, cost, etot(iese_eff + 1), 1)
            cost = 1.d0/wtotf(2)
            call dscal(iese_eff, cost, etot, 1)
            if (decoupled_run) then
                wtotf = wtotf/(nw/nk)
            else
                wtotf = wtotf/nw
            end if
            if (iese .ge. 1) etot(1:ipc) = etot(1:ipc)*ris(2)
            if (isfix .ge. 1) etot(indfix) = etot(indfix)*ris(4)
            if (isfix .ge. 2) etot(indfix + 1) = etot(indfix + 1)*ris(2)
            if (nrep .ne. 0) call dscal(nrep, ris(2), etot(repf), 1)
            if (ieskin .ne. 0) then
                call dscal(ieskin, ris(1), etot(kinf), skipforce)
                if (skipforce .eq. 3) then
                    call dscal(ieskin, ris(3), etot(kinf + 1), skipforce)
                    call dscal(ieskin, ris(3), etot(kinf + 2), skipforce)
                end if
            end if

        end if

        if (itestr .eq. -5) then ! for optimization only
            timeoptp = cclock()
            call updatewfopt
            timeopt = timeopt + cclock() - timeoptp
        end if

        call writeandbranch ! branching and writing fort.12/fort.12.fn/fort.12.new

        ! progress indicator
        if ((i_main - iend) .ge. nprogress) then
            if (rank .eq. 0)&
               & write (6, '(I15," steps, ",I3," % done in ", f10.2," sec")')&
               & i_main - iend, int(nprogress*100.d0/ngen), cclock() - time1p
            time1p = cclock()
            do while ((i_main - iend) .ge. nprogress)
                nprogress = nprogress + max(1, int(ngen*0.01))
            end do
        end if
        !

        if (itestr .ne. -5) then
            cclock_zero = cclock() - inittime
            ! all have to be syncronized with the master
#ifdef PARALLEL
            call mpi_bcast(cclock_zero, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
        end if
    end do ! end of the main loop over i_main
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
    timepp = cclock()
    time = timepp - time
    time_main = time_main + timepp - timeppp
    ! Total time without the end writing.

    !======== finalizing QM/MM scheme ==========!

    ! by. E. Coccia (29/11/10): write the final additional energy
    if (ext_pot .and. itestr .ne. -5) then
#ifdef PARALLEL
        ! electrons
        call mpi_reduce(ave, t_ave, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ave2, t_ave2, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ncount, t_ncount, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(nout, t_nout, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        ! nuclei
        call mpi_reduce(ave_ion, t_ave_ion, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ave2_ion, t_ave2_ion, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ncount_ion, t_ncount_ion, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
        ! electrons
        t_ave = ave
        t_ave2 = ave2
        t_ncount = ncount
        t_nout = nout
        ! nuclei
        t_ave_ion = ave_ion
        t_ave2_ion = ave2_ion
        t_ncount_ion = ncount_ion
#endif

        ! write the external potentials
        if (rank .eq. 0) then
            ! electronic potential
            call extpot_final(nel)
            ! by E. Coccia (13/1/11): write the nuclear potential
            call ion_final(nion)
        end if

        ! by E. Coccia (4/2/11): write the vdw energy
        if (vdw) then
#ifdef PARALLEL
            call mpi_reduce(ave_vdw, t_ave_vdw, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
            call mpi_reduce(ave2_vdw, t_ave2_vdw, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
            call mpi_reduce(ncount_vdw, t_ncount_vdw, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
            t_ave_vdw = ave_vdw
            t_ave2_vdw = ave2_vdw
            t_ncount_vdw = ncount_vdw
#endif
            ! by E. Coccia (27/5/11): write the classical potential energies
            if (link_atom) then
                ! Bond angle
#ifdef PARALLEL
                call mpi_reduce(ave_angle, t_ave_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_angle, t_ave2_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_angle, t_ncount_angle, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                t_ave_angle = ave_angle
                t_ave2_angle = ave2_angle
                t_ncount_angle = ncount_angle
#endif
                ! Proper dihedral
#ifdef PARALLEL
                call mpi_reduce(ave_dihed, t_ave_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_dihed, t_ave2_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_dihed, t_ncount_dihed, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                t_ave_dihed = ave_dihed
                t_ave2_dihed = ave2_dihed
                t_ncount_dihed = ncount_dihed
#endif
                ! Improper dihedral
#ifdef PARALLEL
                call mpi_reduce(ave_impr, t_ave_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_impr, t_ave2_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_impr, t_ncount_impr, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                t_ave_impr = ave_impr
                t_ave2_impr = ave2_impr
                t_ncount_impr = ncount_impr
#endif
            end if
            if (rank .eq. 0) then
                call vdw_final()
            end if
        end if
    end if
    ! by E. Coccia (10/12/11): MM restraints
    if (mm_restr .and. itestr .ne. -5) then
        ! Bond distance
#ifdef PARALLEL
        call mpi_reduce(ave_bond, t_ave_bond, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ave2_bond, t_ave2_bond, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ncount_bond, t_ncount_bond, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
        t_ave_bond = ave_bond
        t_ave2_bond = ave2_bond
        t_ncount_bond = ncount_bond
#endif
        ! Bond angle
#ifdef PARALLEL
        call mpi_reduce(ave_angle, t_ave_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ave2_angle, t_ave2_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ncount_angle, t_ncount_angle, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
        t_ave_angle = ave_angle
        t_ave2_angle = ave2_angle
        t_ncount_angle = ncount_angle
#endif
        ! Proper dihedral
#ifdef PARALLEL
        call mpi_reduce(ave_dihed, t_ave_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ave2_dihed, t_ave2_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ncount_dihed, t_ncount_dihed, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
        t_ave_dihed = ave_dihed
        t_ave2_dihed = ave2_dihed
        t_ncount_dihed = ncount_dihed
#endif
        ! Improper dihedral
#ifdef PARALLEL
        call mpi_reduce(ave_impr, t_ave_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ave2_impr, t_ave2_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ncount_impr, t_ncount_impr, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
        t_ave_impr = ave_impr
        t_ave2_impr = ave2_impr
        t_ncount_impr = ncount_impr
#endif
        if (rank .eq. 0) then
            call vdw_final()
            write (6, *) ' New Energy (no MM) = ', ener_true(1)*ris(2) - sum_pot, sigma_true(1)*ris(2)
            write (*, *) '|******************************************|'
            write (*, *) '|         EXTERNAL QMC/MM POTENTIAL        |'
            write (*, *) '|******************************************|'
            write (*, *) ''
        end if
    end if

    !       count the final time inside the routine
    if (rank .eq. 0) write (6, *) ' Before finalizeall '

    call Finalizeall
    stop

contains
    subroutine update_projm
        implicit none
#ifdef _OFFLOAD
!$omp target update to(mu_c,detmat_c)
        if (yes_complex) then
            call zgemm_('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, zone, mu_c, nelorbh*ipf&
                 &, detmat_c, nelorb_c, zzero, projm, ipf*nelorbh)
        else
            call dgemm_('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf*nelorbh&
                 &, detmat_c, nelorb_c, 0.d0, projm, ipf*nelorbh)
        end if
!$omp target update from (projm)
#else
        if (yes_complex) then
            call zgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, zone, mu_c, nelorbh*ipf&
                 &, detmat_c, nelorb_c, zzero, projm, ipf*nelorbh, nprocu, rankrep, commrep_mpi)
        else
            call dgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf*nelorbh&
                 &, detmat_c, nelorb_c, 0.d0, projm, ipf*nelorbh, nprocu, rankrep, commrep_mpi)
        end if
#endif
    end subroutine update_projm
    subroutine update_projm_
        implicit none
        !  Outside the target region
        if (yes_complex) then
            call zgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, zone, mu_c, nelorbh*ipf&
                 &, detmat_c, nelorb_c, zzero, projm, ipf*nelorbh, nprocu, rankrep, commrep_mpi)
        else
            call dgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf*nelorbh&
                 &, detmat_c, nelorb_c, 0.d0, projm, ipf*nelorbh, nprocu, rankrep, commrep_mpi)
        end if
    end subroutine update_projm_
    subroutine makeqmcsingleel

        implicit none
        real*8, external :: dnrm2
        real*8 drand1, enercont, mapping, zeta_random, arg1, arg2
        logical :: heat_bath_dmc
        real(8) :: diffuse_norm, costexpn_rej
#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(1) ! scalar code
#endif

        !       This subroutine update the electron coordinates of all walkers
        !       with single electron moves using all possible options (VMC,DMC,...).
        !       During this loop all tables to make fast updates are computed.
        !       No comunication between walkers (and processors) is done inside.
        !       This routine is the core of the code.

        do js = ist, ien

            j = js - istm
            jtype2 = j
            !        calculation indices tables
            indtj = Lztab*(j - 1) + 1
            indtjr = Lztabr*(j - 1) + 1
            indtabj = Ltab*(j - 1) + 1
            indtabbj = Ltabb*(j - 1) + 1
            indkj = nel*(indt + 1)*(j - 1) + 1
            indksj = nel*(j - 1) + 1
            indksij = nelnion*(j - 1) + 1
            indwwj = nel2wt*(j - 1) + 1
            indwwjfn = nel2wtfn*(j - 1) + 1
            indwwjj = nel2wtj*(j - 1) + 1
            indwupj = nel2upt*(j - 1) + 1
            indwdoj = nel2dot*(j - 1) + 1
            indaupj = nel2up*(j - 1) + 1
            indbar = nel2bar*(j - 1) + 1
            indbarfn = nel2barfn*(j - 1) + 1
            indjbar = nel2jbar*(j - 1) + 1
            if (iessz) then
                indjbarsz = indjbar
            else
                indjbarsz = 1
            end if

            nleft = nbra
            tleft = tbra
            vpotint = 0.d0
            diffint = 0.d0

            wint(j) = 0.d0
            enerint(j) = 0.d0
            !      enerintw(j)=0.d0
            irej = 0
            icount = 0
            diffuse_norm = 1.d0
            heat_bath_dmc = .false.
            ttry = tbra

            if (tleft .eq. 0.d0) tleft = 1.d0
            if (itest .eq. 1) nmovet = nmovet + 1.d0 ! count the number of steps

            t_cyrus(j) = 0.d0

            do while ((.not. yesnleft .and. tleft .gt. 0.d0) .or. (yesnleft .and. nleft .gt. 0))

                if (itest .eq. 1) then

                    costwn = wsto(j) - lambda
                    if (fncont) then
                        !  Add the following line to recover the same output of Turbo version #<686.
                        !                 zeta(1)=drand1()
                        if (itestrfn .eq. 1) then
                            if (better_dmc .and. cutreg .ne. -2.d0) then
                                costwn = costwn/nel
                            else
                                ! Eq.39 Umrigar JCP '93 neglecting the irrelevant constant.
                                costwn = costwn*sqrt(gradtotbar(j)/gradtot(j))/nel
                            end if
                            costwnn = costwn
                        elseif (itestrfn .eq. 6) then
                            ! New regularized with the single electron gradient only
                            ! since we have only one electron iout that moves
                            iout = mod(icount, nel) + 1
                            costwn = costwn*dnrm2(3, gradpsibar(1, indksj + iout - 1), 1) &
                                     /dnrm2(3, gradpsi(1, indksj + iout - 1), 1)/nel
                            costwnn = costwn
                        elseif (itestrfn .eq. -2 .and. heat_bath_dmc) then
                            ! update the weights only before the heat bath move
                            ! diffuse_norm corrects for the proper diffusion in the presence of rejection (local variable)
                            costwn = diffuse_norm*costwn*sqrt(gradtotbar(j)/gradtot(j))
                            costwnn = costwn
                            diffuse_norm = 1.d0
                        elseif (itestrfn .eq. -3 .and. heat_bath_dmc) then
                            ! update the weights only before the heat bath move
                            ! in this case the heat bath is performed after every single particle move
                            if (irej .eq. 1 .and. rejweight) then
                                costwn = 0.d0
                                ! I can apply the standard renormalization
                                ! of the weights in case of rejection
                                ! namely if irej.eq.1 (previous diffusion
                                ! move rejected) and rejweight option is true
                                ! do not weight the walker
                            else
                                if (better_dmc .and. cutreg .ne. -2.d0) then
                                    costwn = costwn/nel
                                else
                                    costwn = costwn*sqrt(gradtotbar(j)/gradtot(j))/nel
                                end if
                                ! Removed non size consistent term
                                !                     costwn=costwn/nel
                            end if
                            costwnn = costwn
                        end if
                    else
                        costwnn = costwn ! given by the local energy of the fermionic wf Sign Det x |psi_G| in any case
                        costwn = costwn + diag(j) - diagfn(j)
                    end if

                    if (fncont) then
                        if ((itestrfn .eq. -2 .or. itestrfn .eq. -3) .and. icount .eq. nbra) then
                            if (heat_bath_dmc) tleft = ttry
                            ! the last iteration is a non local pseudo move
                        elseif ((itestrfn .ne. -2 .and. itestrfn .ne. -3) .and. icount .eq. nbram) then
                            tleft = ttry
                            ! the last iteration is a diffusion move
                        else
                            tleft = ttry*nbra
                            ! I don't care what is tleft, I don't need to exit the do while
                        end if
                    else
                        pdiag = diagfn(j) - wsto(j)
                        zeta(1) = drand1()
                        if (pdiag .lt. 0.d0) then
                            ttry = dlog(1.d0 - zeta(1))/pdiag
                            if (ttry - tleft .ge. 0.d0 .and. .not. yesnleft) ttry = tleft
                        else
                            ttry = tleft
                            ! kill the walker
                            wconfn(js) = 0.d0
                        end if
                        ttrysav = ttry
                        wconfnsav = wconfn(js)
                    end if

                else !itest.eq.1

                    nmovet = nmovet + 1.d0
                    iout = drand1()*nel + 1
                    indold = indkj + iout - 1
                    jn = drand1()*3 + 1

                    call dcopy(3, kel(1, indold), 1, rcart, 1)

                    f = mapping(nion, dist(indksij + (iout - 1)*nion), zetar, bcost, ccost   &
                         &, imin, imax)

                    ! gaussian move
                    arg1 = 1.d0 - drand1()
                    arg2 = TWO_PI*drand1()

                    dstep(jn) = f*tstep*dsqrt(-2.d0*dlog(arg1))*dcos(arg2)

                    call hopping(nion, jn, dstep(jn), rcart, rion, imin, itry, hopfraction &
                         &, LBox)

                    if (LBox .le. 0.d0) then
                        do kk = 1, nion
                            dists(kk) = (rcart(1) - rion(1, kk))**2 + &
                                 &      (rcart(2) - rion(2, kk))**2 + (rcart(3) - rion(3, kk))**2
                            dists(kk) = dsqrt(max(dists(kk), 1d-18))
                            !                     if(dists(kk).le.1d-18) dists(kk)=1d-18
                        end do
                    else
                        !                  call dscalzero(nion,0.d0,dists,1)
                        !                  call dscalzero(3*nion,0.d0,dists_kel,1)
                        do kk = 1, nion
                            dists_kel(:, kk) = rcart(:) - rion(:, kk)
                        end do

                        call ApplyPBC(dists_kel, nion)

                        dists(:) = dists_kel(1, :)**2                                  &
                             & + dists_kel(2, :)**2 + dists_kel(3, :)**2
                        dists(:) = dsqrt(max(dists(:), 1d-18))
                        !                  do kk=1,nion
                        !                     if(dists(kk).le.1d-18) dists(kk)=1d-18
                        !                  enddo

                    end if
                    !               call my_dsqrt(nion,dists,dists)

                    fb = mapping(nion, dists, zetar, bcost, ccost, imin, imaxn)

                    if (dstep(jn) .ne. 0.d0 .or. imin .eq. itry) then ! accepted move

                        ratiodet = psidetln(j)

#if defined (_OPENMP) && defined (__NOOMP)
                        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

                        timep = cclock()

                        call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                             &, ainv(indaupj), kel(1, indkj), rcart, iesdr, vj, dupr                   &
                             &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr            &
                             &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur                                 &
                             &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar(indbar)        &
                             &, winvjbar(indjbar), winvjbarsz(indjbarsz), iflagnorm, cnorm, iflagerr &
                             &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne       &
                             &, jastrowall_ee(1, iout, 0, j), jasnew_ee, jastrowall_ei(1, iout, j)      &
                             &, jasnew_ei, n_body_on, niesd, nshell, epscuttype&
                             &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                        time_ratiovar = time_ratiovar + cclock() - timep

                        !write(6,*) 'ratio after ratiovar =',ratior,ratiodet

#if defined (_OPENMP) && defined (__NOOMP)
                        call omp_set_num_threads(1) ! restore the previous threads
#endif

                        if (epscuttype .gt. 0) then
                            if (ratiodet .ne. 0.d0) then
                                call ratio_psi(psidetln(j), ratiodet, epscutu, ratioreg        &
                                     &, spsi, epstlu)
                                ratiorn(1:ipc) = (ratioreg/(ratiodet/psidetln(j)))*ratior(1:ipc)
                            else
                                ratiorn = 0.d0
                            end if
                        else
                            ratiorn = ratior
                        end if

                        !          write(6,*) ' conf =',i,kel(jn,iout)
                        !          write(*,*) 'ratio,ratiodet',ratior,ratiodet
                        !          write(6,*) ' drand1() =',drand1()

                        !    Correct the acceptance for the non trivial move
                        if (dstep(jn) .ne. 0.d0 .and. epscuttype .ne. -100) then
                            cost = 0.5d0/tstep**2
                            ratio = sum(ratiorn(1:ipc)**2)* &
                                 &  f/fb*dexp(-dstep(jn)**2                                         &
                                 & *(cost/fb**2 - cost/f**2))
                        else
                            ratio = sum(ratiorn(1:ipc)**2)
                        end if

                        !write(6,*) ' Final ratio =',ratiorn,ratio
                        !stop

                        ! Thermal bath
                        if (tstep .lt. 0) ratio = ratio/(1.d0 + ratio)

                    else ! if dstep(jn).ne.0.d0.or.imin.eq.itry

                        ! reject proposed hopping
                        ratio = 0.d0

                    end if

                    zeta_random = drand1()

                    if (ratio .le. zeta_random .and. epscuttype .ne. -100) then
                        ntry = 1 ! rejected move
                    else
                        ntry = 0 ! accepted move
                        if (itry .eq. imax .and. dstep(jn) .eq. 0) acclarge = acclarge + 1.d0
                    end if

                end if ! end if for itest.eq.1

                !write(6,*) 'INFO gaussian move',itest,ttry,tleft,fncont,yesnleft

                if ((itest .eq. 1 .and. ttry .ge. tleft .and.                           &
                     &.not. fncont .and. .not. yesnleft) .or. (itest .eq. 2 .and. ntry .eq. 1)) then

                    if (itest .eq. 1) then

                        if (ttry .ne. tleft) then
                            write (errmsg, *) ' error impossible condition ttry > tleft ', ttry, tleft
                            call checkiflagerr(1, rank, errmsg)
                        end if

                        if (wconfn(js) .ne. 0.d0) then

                            if (parcutg .ge. 0) then
                                !                     costexp=dexp(costwn*tleft)
                                if (yesnleft) then
                                    if (lambda .gt. diagfn(j)) then
                                        costexpn = -(diagfn(j) - wsto(j))/(lambda - diagfn(j))
                                    else
                                        costexpn = 0.d0
                                    end if
                                else
                                    costexpn = dexp(costwnn*tleft)
                                end if
                            else ! linear approximation more stable but biased
                                !                     costexp=max(1.d0+costwn*tleft,0.d0)
                                costexpn = max(1.d0 + costwnn*tleft, 0.d0)
                            end if

                            if (.not. fncont) then

                                if (parcutg .ge. 0) then
                                    if (yesnleft) then
                                        if (lambda - diagfn(j) .gt. 0.d0) then
                                            if (.not. yes_cutweight) then
                                                cost = wconfn(js)/(lambda - diagfn(j))
                                            else
                                                cost = wconfn(js)*(1.d0/(lambda - diagfn(j)))**(1.d0 + cutweight*alat**2)
                                            end if
                                        else
                                            cost = 0.d0
                                        end if
                                    else
                                        if (abs(costwnn*tleft) .gt. 1d-7) then
                                            cost = (costexpn - 1.d0)/costwnn*wconfn(js)
                                        else
                                            cost = tleft*wconfn(js)
                                        end if
                                    end if
                                elseif (parcutg .lt. 0) then
                                    cost = tleft*wconfn(js)
                                end if

                                enerint(j) = enerint(j) + enertrue(j)*cost
                                ! In case epscutdmc>0 use the bosonic local energy which is smooth.
                                vpotint = vpotint + cost*diffkin(2, j)
                                diffint = diffint + cost*diffuse(j)
                                wint(j) = wint(j) + cost

                                wconfn(js) = wconfn(js)*costexpn

                            else

                                enerint(j) = enerint(j) + enert(1, j)*wconfn(js)
                                !                     enerintw(j)=enerintw(j)+enert(j)*wconf(js)
                                vpotint = vpotint + wconfn(js)*diffkin(2, j)
                                diffint = diffint + wconfn(js)*diffuse(j)
                                wint(j) = wint(j) + wconfn(js)
                                !                     wintw(j)=wintw(j)+wconf(js)

                                !                     wconf(js)=wconf(js)*costexp
                                if (rejweight) costexpn_rej = wconfn(js) ! save the old weight
                                wconfn(js) = wconfn(js)*costexpn

                                if (.not. rejweight) costexpn_rej = wconfn(js) ! save the old weight

                                !                     costexp_rej=costexp ! it is used later for rejection

                                ! end if itestr.ne.1
                            end if

                            ! end if wconf.ne.0
                        end if

                        tleft = 0.

                        ! end if itest.eq.1
                    end if

                    if (yesnleft) nleft = nleft - 1

                    t_cyrus(j) = t_cyrus(j) + ttry

                else ! accept a new move

                    if (itest .eq. 1) then

                        if (parcutg .ge. 0) then
                            !                  costexp=dexp(costwn*ttry)
                            if (yesnleft) then
                                if (lambda .gt. diagfn(j)) then
                                    costexpn = -(diagfn(j) - wsto(j))/(lambda - diagfn(j))
                                else
                                    costexpn = 0.d0
                                end if
                            else
                                costexpn = dexp(costwnn*ttry)
                            end if
                        else
                            !                  costexp=max(1.d0+costwn*ttry,0.d0)
                            costexpn = max(1.d0 + costwnn*ttry, 0.d0)
                        end if

                        if (.not. fncont) then

                            if (parcutg .ge. 0) then

                                if (yesnleft) then
                                    if (lambda - diagfn(j) .gt. 0.d0) then
                                        if (.not. yes_cutweight) then
                                            cost = wconfn(js)/(lambda - diagfn(j))
                                        else
                                            cost = wconfn(js)*(1.d0/(lambda - diagfn(j)))**(1.d0 + cutweight*alat**2)
                                        end if
                                    else
                                        cost = 0.d0
                                    end if
                                else
                                    if (abs(costwnn*ttry) .gt. 1d-7) then
                                        cost = (costexpn - 1.d0)/costwnn*wconfn(js)
                                    else
                                        cost = ttry*wconfn(js)
                                    end if
                                end if
                            elseif (parcutg .lt. 0) then

                                cost = ttry*wconfn(js)

                            end if

                            enerint(j) = enerint(j) + enertrue(j)*cost
                            vpotint = vpotint + cost*diffkin(2, j)
                            diffint = diffint + cost*diffuse(j)
                            wint(j) = wint(j) + cost

                            !             if(yesnleft.and.better_lrdmc) then
                            !              wconfn(js)=wconfn(js)*(wsto(j)-diagfn(j))/(lambda-diagfn(j))
                            !              else
                            wconfn(js) = wconfn(js)*costexpn
                            !             endif

                            t_cyrus(j) = t_cyrus(j) + ttry

                        else ! if(not.fncont)

                            if ((itestrfn .ne. -2 .and. itestrfn .ne. -3) .or. heat_bath_dmc) then
                                ! in standard dmc weight the walkers after each single particle move
                                ! in non local dmc weight the walkers before a heat bath move
                                enerint(j) = enerint(j) + enert(1, j)*wconfn(js)
                                !                     enerintw(j)=enerintw(j)+enert(j)*wconf(js)
                                vpotint = vpotint + wconfn(js)*diffkin(2, j)
                                diffint = diffint + wconfn(js)*diffuse(j)
                                wint(j) = wint(j) + wconfn(js)
                                !                     wintw(j)=wintw(j)+wconf(js)
                                !                     wconf(js)=wconf(js)*costexp
                                if (rejweight) costexpn_rej = wconfn(js)
                                wconfn(js) = wconfn(js)*costexpn
                                if (.not. rejweight) costexpn_rej = wconfn(js)
                                !                        costexp_rej=costexp  ! it is used later for rejection in standard dmc
                            end if

                        end if ! if(not.fncont)

                        tleft = tleft - ttry

                    end if ! if(itest.eq.1)

                    if (yesnleft) nleft = nleft - 1
                    if (itest .ne. 1) then
                        psidetln(j) = ratiodet
                        if (spsi .ne. 1.d0) nontr = nontr + 1.d0
                        !               No longer needed this dirty trick below
                        !               if(kaverage) then
                        !                  wconfn(js)=spsi*scale_spsi*wkp(ikpoint)
                        !               else
                        wconfn(js) = spsi
                        !               endif
                    end if

                    if (.not. heat_bath_dmc) nacc = nacc + 1.d0

                    if (itest .eq. 1) then

                        if (indt .ne. 0 .and. (.not. fncont .or. heat_bath_dmc)) then
                            ! if DMC with non local moves, do non local move when heat_bath_dmc is .true.
                            ! if LRDMC pass here always

                            if (better_dmc) then
                                call enforce_detail(indt, nel, identity, gamma, tabler(indtjr), norm_tab)
                                zeta(1) = norm_tab*(1.d0 - drand1())
                            else
                                zeta(1) = (1.d0 - drand1())*(wsto(j) - diagfn(j))
                            end if

                            if (typereg .ge. 6 .or. epscutdmc .ne. 0) then
                                call random(nel, Lzeff, zeta(1), table(indtj), 0.d0, iout, indvic, &
                                     & psign, npow, gamma, istart, indtm(1, j), ipc)
                            else
                                call random(nel, Lzeff, zeta(1), tabler(indtjr), 0.d0, iout, indvic, &
                                     & psign, npow, gamma, istart, indtm(1, j), 1)
                            end if

                            if (epscutdmc .ne. 0.d0 .and. novar .eq. 0) then
                                ii = indtjr + (indvic - 1)*nel + iout - 1
                                jj = indtj + (indvic - 1)*nel*ipc + ipc*(iout - 1)
                                reweight_dmc = table(jj)/tabler(ii)
                            end if

                            if (indvic .lt. 0) then
                                iflagerr = 1
                                write (6, *) ' Error in random in processor #,walker #', rank, j
                                indvic = -indvic
                            end if
                            if (.not. fncont) then
                                enerold = diffkin(3, j)
                                psioldinp = psidetln(j)
                                if (npsa .gt. 0 .and. indvic .ge. istart) naccpseudo = naccpseudo + 1
                            else
                                ! if(heat_bath_dmc) then
                                ! non local move for pseudopotentials with DMC
                                if (indvic .le. indtm(iout, j)) then
                                    ! perform the non local move only for non trivial rotations (do nothing for the identity)
                                    naccpseudo = naccpseudo + 1
                                    irej = 0

                                    indp = indksj + iout - 1
                                    ind = indkj + iout - 1
                                    do jn = 1, 3
                                        rcart(jn) = kel(jn, ind) + ivic(jn, indvic, indp)
                                    end do
                                else
                                    ! reject the move (identity or zero step move)
                                    irej = 1
                                end if
                            end if !if(.not.fncont)
                        else !if(indt.ne.0)
                            ! generate a new coordinate with the diffusion algorithm
                            iout = mod(icount, nel) + 1
                            indold = indkj + iout - 1
                            icount = icount + 1
                            indvic = 1
                            irej = 0
                            call dcopy(3, kel(1, indold), 1, rcarto, 1)
                            ! gaussian move+drift

                            do jn = 1, 3
                                arg1 = 1.d0 - drand1()
                                arg2 = TWO_PI*drand1()
                                dstep(jn) = dsqrt(-4.*ttry*dlog(arg1))*dcos(arg2)
                                rcart(jn) = kel(jn, indold) + dstep(jn) &
                                            + 2.d0*gradpsibar(jn, indksj + iout - 1)*ttry
                            end do

                        end if ! endif for indt.ne.0
                    end if !endif itest.eq.1

                    if (fncont .and. irej .eq. 0) then
                        timep = cclock()
#if defined (_OPENMP) && defined (__NOOMP)
                        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
                        call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior       &
                             &, ainv(indaupj), kel(1, indkj), rcart, iesdr, vj, dupr      &
                             &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr            &
                             &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur                                 &
                             &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar(indbar)         &
                             &, winvjbar(indjbar), winvjbarsz(indjbarsz), iflagnorm, cnorm, iflagerr &
                             &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne       &
                             &, jastrowall_ee(1, iout, 0, j), jasnew_ee, jastrowall_ei(1, iout, j)      &
                             &, jasnew_ei, n_body_on, niesd, nshell, epscuttype&
                             &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                        time_ratiovar = time_ratiovar + cclock() - timep

#if defined (_OPENMP) && defined (__NOOMP)
                        call omp_set_num_threads(1) ! restore the scalar
#endif
                        if (ratiodet .le. 0.d0 .and. .not. heat_bath_dmc) then
                            ! cross the node reject conf
                            ! note the condition .not.heat_bath_dmc
                            ! if heat_bath_dmc the accepted heat bath move can cross the node
                            ! provided the sign of table is positive and so we do not need to reject the move
                            irej = 1
                            if (itestrfn .ne. -2 .and. itestrfn .ne. -3) then
                                !                     wconf(js)=wconf(js)/costexp_rej ! Reject and rescale weight in standard dmc
                                wconfn(js) = costexpn_rej ! Reject and rescale weight in standard dmc
                            elseif (itestrfn .eq. -2 .and. rejweight) then
                                diffuse_norm = diffuse_norm - 1.d0/dble(nel)
                                ! renormalization of the time step in weighting
                                ! used for non local dmc itestrfn=-2 (non local move after all electron sweep)
                                ! rejweight false -> do not renormalize the weight after rejection
                            end if
                        end if
                    end if ! if(fncont.and.irej.eq.0)

                    if (irej .eq. 0) then
                        naccm(j) = mod(naccm(j) + 1, nscra)
                        timep = cclock()
#if defined (_OPENMP) && defined (__NOOMP)
                        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
                        ! fast updating of the w.f. (VMC case)
                        call uptabtot_global(js, pseudologic)
#if defined (_OPENMP) && defined (__NOOMP)
                        call omp_set_num_threads(1) ! restore the scalar
#endif
                        time_uptabtot = time_uptabtot + cclock() - timep
                        ! Here is the rejection scheme
                        if (fncont .and. .not. heat_bath_dmc) then
                            ! condition .not.heat_bath_dmc : only the diffusion move can be rejected

                            gradbarold = gradtotbar(j)
                            gradold = gradtot(j)
                            call dcopy(3*nel, gradpsibar(1, indksj), 1, gradpsibarold, 1)
                            call dcopy(3*nel, gradpsi(1, indksj), 1, gradpsiold, 1)

                            ! evaluation gradient  psi
                            call upgradcont(gradpsibar(1, indksj), gradpsi(1, indksj)      &
                                 &, indt, nelup, neldo, winvup(indwupj), winvdo(indwdoj)                 &
                                 &, tabpip(indtabj), ttry, gradtotbar(j), gradtot(j), kel(1, indkj)       &
                                 &, dist(indksij), rion, nion, zetar, LBox)

                            do jn = 1, 3
                                dx_old(jn) = dstep(jn)
                                dx_new(jn) = -(dstep(jn)                         &
                                     & + 2.d0*gradpsibarold(jn, iout)*ttry)      &
                                     & - 2.d0*gradpsibar(jn, indksj + iout - 1)*ttry
                            end do

                            pdiag = sum(ratior(1:ipc)**2)*dexp(-0.25d0*(dx_new(1)**2 + dx_new(2)**2 + &
                                 & dx_new(3)**2 - dx_old(1)**2 - dx_old(2)**2 - dx_old(3)**2)/ttry)

#if defined (_OPENMP) && defined (__NOOMP)
                            call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

                            if (drand1() .ge. pdiag) then
                                ! reject the configuration
                                ! we compute again the wf on the previous position, which did not change
                                ! we assume that the rejection is very unlikely (true in the small time step limit)
                                irej = 1
                                call dcopy(3, rcarto, 1, rcart, 1)
                                gradtotbar(j) = gradbarold
                                gradtot(j) = gradold
                                call dcopy(3*nel, gradpsibarold, 1, gradpsibar(1, indksj), 1)
                                call dcopy(3*nel, gradpsiold, 1, gradpsi(1, indksj), 1)

                                timep = cclock()

                                call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior       &
                                     &, ainv(indaupj), kel(1, indkj), rcart, iesdr, vj, dupr      &
                                     &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr            &
                                     &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur                                 &
                                     &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar(indbar)         &
                                     &, winvjbar(indjbar), winvjbarsz(indjbarsz), iflagnorm, cnorm, iflagerr &
                                     &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne       &
                                     &, jastrowall_ee(1, iout, 0, j), jasnew_ee, jastrowall_ei(1, iout, j)      &
                                     &, jasnew_ei, n_body_on, niesd, nshell, epscuttype&
                                     &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                                time_ratiovar = time_ratiovar + cclock() - timep

                                timep = cclock()

                                call uptabtot_global(js, pseudologic)

                                time_uptabtot = time_uptabtot + cclock() - timep

                                if (itestrfn .ne. -2 .and. itestrfn .ne. -3) then
                                    ! wconf(js)=wconf(js)/costexp_rej
                                    ! Reject and rescale the weight in standard dmc
                                    wconfn(js) = costexpn_rej ! Reject and rescale the weight in standard dmc
                                elseif (itestrfn .eq. -2 .and. rejweight) then
                                    diffuse_norm = diffuse_norm - 1.d0/dble(nel)
                                    ! renormalization of the time step in weighting
                                    ! used for non local dmc itestrfn=-2 (non local move after all electron sweep)
                                    ! rejweight false -> do not renormalize the weight after rejection
                                end if

                                naccm(j) = mod(naccm(j) + 1, nscra)

                            end if !if(drand1().ge.pdiag)
                        end if ! if(fncont.and..not.heat_bath_dmc)

                        ! endif irej.eq.0
                    end if

#if defined (_OPENMP) && defined (__NOOMP)
                    call omp_set_num_threads(old_threads) ! restore the threads
#endif

                    if (irej .eq. 1 .and. .not. heat_bath_dmc) then
                        nacc = nacc - 1.d0
                    end if

                    !if(psisn(j).eq.0.d0)  then
                    if (singdet(j)) then
                        ifz = 1
                    else
                        ifz = 0
                    end if

                    if (itest .eq. 2 .and. naccm(j) .eq. 0) then

                        call scratchdet(winv(indwwj), winvbar(indbar), ainv(indaupj), psidetln(j)&
                             &, psip, psip(ipc*nelup_mat*nelup_mat + 1), ipsip, iflagerr)
                        if (itest .ne. 1 .and. iflagerr .eq. 0) then

                            if (epscuttype .eq. 2) then
                                call psireg(ainv(indaupj), agpn, psip, nelup_mat, neldo, psidetln(j), yes_ontarget)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
                                do jj = 1, nelup_mat
                                    do ii = 1, ipc*nelup_mat
                                        agp(ii, jj, jtype2) = agpn(ii, jj)
                                    end do
                                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif

                            end if
                            ! Otherwise one may  repeat scratchdet even with no acceptance
                            naccm(j) = 1
                        end if

                        !write(6,*) 'detlogA',psidetln(j),ainv(indaupj)

                        if (iflagerr .ne. 0) then
                            if (nw .gt. 1) then
                                write (6, *) ' Warning walker singular !!! ', js
                                !                     wconf(js)=0.d0
                                wconfn(js) = 0.d0
                                if (iesbra) iflagerr = 0 ! trying to continue
                            end if
                        end if
                        !
                    end if

                    if (ifz .eq. 1 .or. (naccm(j) .eq. 0 .and. itest .ne. 2)) then
                        timep = cclock()

                        call upscratch_global(js, .false., .false.)

                        timescra = timescra + cclock() - timep

                        !if(psisn(j).eq.0) then
                        if (singdet(j)) then
                            write (6, *) ' Warning psi singular !!! # walker= ', js
                            if (nw .gt. 1) then
                                !                     wconf(js)=0.d0
                                wconfn(js) = 0.d0
                                if (iesbra) iflagerr = 0 ! trying to continue
                            end if
                        end if

                        ! fine if ifz
                    end if

#if defined (_OPENMP) && defined (__NOOMP)
                    call omp_set_num_threads(1) ! restore the scalar
#endif

                    if (itestrfn .eq. -2 .or. itestrfn .eq. -3) then
                        if (((itestrfn .eq. -2 .and. mod(icount, nel) .eq. 0) .or. itestrfn .eq. -3) &
                            .and. .not. heat_bath_dmc) then
                            ! itestrfn.eq.-2: after an all electron diffusion sweep,
                            ! perform (next) a heat bath move based on non local potentials
                            ! itestrfn.eq.-3: after every single particle diffusion move, perform (next) a heat bath move
                            heat_bath_dmc = .true.

                        elseif (heat_bath_dmc) then
                            ! after a heat bath move, update gradients and do (next) a standard diffusion move
                            call upgradcont(gradpsibar(1, indksj), gradpsi(1, indksj)      &
                                 &, indt, nelup, neldo, winvup(indwupj), winvdo(indwdoj)                 &
                                 &, tabpip(indtabj), ttry, gradtotbar(j), gradtot(j), kel(1, indkj)       &
                                 &, dist(indksij), rion, nion, zetar, LBox)
                            heat_bath_dmc = .false.

                        end if
                    end if

                    if (itest .ne. 2) then

                        !if(psisn(j).ne.0) then
                        if (.not. singdet(j)) then
                            if (yes_complex) then
                                call uptable_complex(nelup, neldo, indtupt, table(indtj), tabler(indtjr)&
                                     &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), tmu(indtabbj)&
                                     &, epscutdmc, psiln(j), psidetln(j), costa, parcutg, istart, typereg)

                                call updiag_complex(table(indtj), tmu(indtabbj), diag(j), enert(1, j)&
                                     &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), nelup, neldo, nel&
                                     &, indt, indtupt, indteff, istart, vpot(j), vpotreg(1, indksj), parcutg, vcut(j) &
                                     &, diffkin(1, j), novar, yescut(j), ener_c, voffpseudo, psip, psip(2*nel + 1))
                                enertrue(j) = dreal(ener_c)
                            else

                                call uptable(nelup, neldo, indtupt, table(indtj), tabler(indtjr)&
                                     &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), tmu(indtabbj)&
                                     &, epscutdmc, psiln(j), psidetln(j), costa, parcutg, istart, typereg)

                                call updiag(table(indtj), tmu(indtabbj), diag(j), enert(1, j)&
                                     &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), nelup, neldo, nel&
                                     &, indt, indtupt, indteff, istart, vpot(j), vpotreg(1, indksj), parcutg, vcut(j) &
                                     &, diffkin(1, j), novar, yescut(j), enertrue(j), voffpseudo, psip, psip(nel + 1))

                            end if

                            if (iese .ge. 3) then
                                call diffus(nel, indt, gamma, ivic(1, 1, indksj)&
                                     &, table(indtj), diffuse(j), istart)
                            end if

                            if (.not. yes_complex) then
                                call energy(Lzeff, table(indtj), tabler(indtjr), diag(j), wsto(j), veff&
                                     &, veffright, enerdiff, itest, npow, gamma, nel, istart, indtm(1, j), epscutdmc)
                            else
                                call energy_complex(Lzeff, table(indtj), tabler(indtjr), diag(j), wsto(j), veff&
                                     &, veffright, enerdiff, itest, npow, gamma, nel, istart, indtm(1, j), epscutdmc)
                            end if

                            if (add_diff) then
                                !                  change the definition of the local energy
                                enertrue(j) = enertrue(j) + enerdiff
                                enert(1, j) = enert(1, j) + enerdiff
                                if (ipc .eq. 2) ener_c = ener_c + enerdiff
                            end if
                            diffkin(3, j) = diffkin(3, j) + enerdiff

                            if (novar .eq. 3) enertrue(j) = -wsto(j)

                            if (isfix .ne. 0) then
                                if (ipc .eq. 2) then
                                    econf(j + nwfix) = ener_c*dconjg(ener_c)
                                else
                                    econf(j + nwfix) = enertrue(j)**2
                                end if
                            end if

                            diagfn(j) = diag(j) + (1.d0 + gamma)*veff + npow*(1.d0 + gamma)*veffright
                            if (better_dmc) call cutwstodmc(yesalfe, nel&
                                                         &, wsto(j), lambda&
                                                         &, diffkin(3, j)&
                                                         &, psip(ipc*nel + 1)&
                                                         &, cutreg)

                            if (rsignr .ne. 0.d0) then
                                diagfn(j) = diagfn(j) + rsignr*(diffkin(3, j) + lambda)
                                wsto(j) = wsto(j) + rsignr*(diffkin(3, j) + lambda)
                            end if

                            ! fine if psi>0
                        end if
                        if (epstldmc .ne. 0.d0 .and. psidetln(j) .lt. epstldmc) then
                            ! kill the walker if the Determinant is dirty
                            countcut = countcut + 1.d0
                            !                  wconf(js)=0.d0
                            wconfn(js) = 0.d0
                        end if

                        enercuto = enerold
                        enerold = 1.d0

                    end if

                end if ! endif   accepted meas

#if defined (_OPENMP) && defined (__NOOMP)
                call omp_set_num_threads(1) ! restore the scalar  anyway
#endif

                ! end do for the tleft or nleft
            end do

            if (itest .ne. 2) then
                if (wint(j) .ne. 0.d0) then
                    enerint(j) = enerint(j)/wint(j)

                    if (iese .ge. 1) econf(j) = enerint(j)
                    if (iese .ge. 2) econf(j + in1) = vpotint/wint(j)*ris(2)
                    if (iese .ge. 3) econf(j + 2*in1) = diffint/wint(j)

                end if
                !         if(wintw(j).ne.0.d0) enerintw(j)=enerintw(j)/wintw(j)
                if (iese .ge. 4) econf(j + 3*in1) = vcut(j)
                if (iese .ge. 5) econf(j + 4*in1) = diffkin(1, j)
                if (iese .ge. 6) econf(j + 5*in1) = diffkin(2, j)
            else
                wint(j) = wconfn(js)
            end if
            ! Definition factorsr
            if (itest .eq. 2 .or. fncont .or. .not. yesnleft) then
                !         VMC or DMC
                factorsr(j) = wconfn(js)
            else
                !   LRDMC integrated weight in the future average  random time step interval
                if (lambda - diagfn(j) .gt. 0.d0) then
                    !            1/(V(x)-E_0)
                    if (.not. yes_cutweight) then
                        ttrysav = 1.d0/(lambda - diagfn(j))
                    else
                        ttrysav = (1.d0/(lambda - diagfn(j)))**(1.d0 + cutweight*alat**2) ! The new time on average
                    end if
                    factorsr(j) = ttrysav*wconfn(js)
                else
                    ttrysav = 0.d0 ! Kill this weight.
                end if
                factorsr(j) = ttrysav*wconfn(js)
                ! Reweighting factor when Psi_T =/ Psi_G
                factorsr(j) = factorsr(j)*reweight_dmc
            end if
        end do

        ! end parallel evolution
        ! --------------fine DO sui walkers--------------------------
#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

    end subroutine makeqmcsingleel

    subroutine makeallmeas

        !----------------------------------------------------
        use extpot, only: write_rwalk

        implicit none
        !        This subroutine compute by scratch all correlation functions
        !        used for the minimization of the energy or for the dynamic
        !        if flagcont=.false. it is a bit faster because does not
        !        recompute by scratch the normalization of the coefficients
        !        in the orbitals defining the AGP.
        !      This object  should be divided in three:
        !  a)  calculation of all necessary for energy and few other correlations
        !  b)  calculation of all necessary for energy minimization
        !  c)  calculation of ionic forces
        !  But this inside the main loop over js
        real*8 drand1, enercont, mapping
        integer firstmolc, nmolfnc
        logical passicomp

        ! parallel scratch

        do js = ist, ien

            j = js - istm
            indtj = Lztab*(j - 1) + 1
            indtjr = Lztabr*(j - 1) + 1
            indtabj = Ltab*(j - 1) + 1
            indtabbj = Ltabb*(j - 1) + 1
            indkj = nel*(indt + 1)*(j - 1) + 1
            indksj = nel*(j - 1) + 1
            indksij = nelnion*(j - 1) + 1
            indwwj = nel2wt*(j - 1) + 1
            indwwjfn = nel2wtfn*(j - 1) + 1
            indwwjj = nel2wtj*(j - 1) + 1
            indwupj = nel2upt*(j - 1) + 1
            indwdoj = nel2dot*(j - 1) + 1
            indaupj = nel2up*(j - 1) + 1
            indbar = nel2bar*(j - 1) + 1
            indbarfn = nel2barfn*(j - 1) + 1
            indjbar = nel2jbar*(j - 1) + 1
            if (iessz) then
                indjbarsz = indjbar
            else
                indjbarsz = 1
            end if
            timep = cclock()

            if (flagcont) iflagnorm = 3 ! initialize by scratch normalizations
            passicomp = .false.

            ! Restore winvbar when not updated

            if (developer .eq. -1 .and. flagcont .and. rank .eq. 0) then
                if (i_main .eq. 1) then
                    kelsav(1:3, 1:nel) = kel(1:3, indkj:indkj + nel - 1)
                    write (6, *) ' Initial configuration ='
                    do jj = 1, nel
                        write (6, *) jj, kelsav(1:3, indkj + jj - 1)
                    end do
                else
                    write (6, *) ' Reset configuration , iteration =', i
                    kel(1:3, indkj:indkj + nel - 1) = kelsav(1:3, 1:nel)
                end if
            end if
            timep = cclock()
            ! upscratch when doing standard VMC
            call upscratch_global(js, pseudologic, iesrandoml)
            if (flagcont .and. developer .eq. -1 .and. rank .eq. 0) then
                write (6, *) ' Jastrowall-ee =', sum(jastrowall_ee(:, :, 0, j))
                write (6, *) ' Jastrowall-ei =', sum(jastrowall_ei(:, :, j))
            end if
            timescra = timescra + cclock() - timep

            naccm(j) = 1 ! Since the upscratch has been already done

            if (itest .eq. 2) then
                ! variance of the w.f.
                psirav = psirav + psidetln(j)*wconfn(js) ! inverse of the w.f.
                countav = countav + wconfn(js)
                if (kaverage) then
                    countt = countt + wkp(ikpoint)*scale_spsi
                else
                    countt = countt + 1.d0
                end if
                counttot = counttot + wconfn(js)
                psiav = psiav + psidetln(j)*wconfn(js)
                psisav = psisav + psidetln(j)**2*wconfn(js)
            end if

            if (.not. singdet(j)) then
                ! local energy computation
                if (yes_complex) then

                    call uptable_complex(nelup, neldo, indtupt, table(indtj), tabler(indtjr)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), tmu(indtabbj)&
                         &, epscutdmc, psiln(j), psidetln(j), costa, parcutg, istart, typereg)

                    call updiag_complex(table(indtj), tmu(indtabbj), diag(j), enert(1, j)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), nelup, neldo, nel&
                         &, indt, indtupt, indteff, istart, vpot(j), vpotreg(1, indksj), parcutg, vcut(j) &
                         &, diffkin(1, j), novar, yescut(j), ener_c, voffpseudo, psip, psip(2*nel + 1))
                    ! take only the real part of the energy/pseudopotential
                    enertrue(j) = real(ener_c)

                    call energy_complex(Lzeff, table(indtj), tabler(indtjr), diag(j), wsto(j), veff&
                         &, veffright, enerdiff, itest, npow, gamma, nel, istart, indtm(1, j), epscutdmc)

                else

                    call uptable(nelup, neldo, indtupt, table(indtj), tabler(indtjr)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), tmu(indtabbj)&
                         &, epscutdmc, psiln(j), psidetln(j), costa, parcutg, istart, typereg)

                    call updiag(table(indtj), tmu(indtabbj), diag(j), enert(1, j)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), nelup, neldo, nel&
                         &, indt, indtupt, indteff, istart, vpot(j), vpotreg(1, indksj), parcutg, vcut(j) &
                         &, diffkin(1, j), novar, yescut(j), enertrue(j), voffpseudo, psip, psip(nel + 1))

                    call energy(Lzeff, table(indtj), tabler(indtjr), diag(j), wsto(j), veff&
                         &, veffright, enerdiff, itest, npow, gamma, nel, istart, indtm(1, j), epscutdmc)

                end if

                if (add_diff) then
                    enertrue(j) = enertrue(j) + enerdiff
                    enert(1, j) = enert(1, j) + enerdiff
                    if (ipc .eq. 2) ener_c = ener_c + enerdiff
                end if
                diffkin(3, j) = diffkin(3, j) + enerdiff

                if (flagcont .and. developer .eq. -1 .and. rank .eq. 0) then
                    write (6, *) ' Energy same conf =', j, enertrue(j), diag(j), tmu(indtabbj), veff, veffright
                end if

                !           write(6,'(A,2X,I3,2X,3F10.6,2X,2I3,2X,3F10.6)') ' Energy same conf =',j,enertrue(j),diag(j),&
                !           tmu(indtabbj),rankrep,rankcolrep,xkp(:,rankcolrep+1)
                !           call mpi_finalize(ierr)
                !           stop

                !------------------------ parameters variation -------------------------!

                if (novar .eq. 3) enertrue(j) = -wsto(j)

                if (itest .eq. 1) then
                    diagfn(j) = diag(j) + (1.d0 + gamma)*veff + npow*(1.d0 + gamma)*veffright
                    if (better_dmc) call cutwstodmc(yesalfe, nel, wsto(j), lambda, diffkin(3, j), psip(ipc*nel + 1), cutreg)
                    if (rsignr .ne. 0.d0) then
                        diagfn(j) = diagfn(j) + rsignr*(diffkin(3, j) + lambda)
                        wsto(j) = wsto(j) + rsignr*(diffkin(3, j) + lambda)
                    end if
                end if

                if (ieskint .gt. 0 .and. (mod((i_main - 1)/nweight, iskipdyn) .ge. iskipdyn - nmore_force)) then
                    yeswrite = .true.
                else
                    if (ieskint .eq. 0) then
                        !             when no dynamics write forces.dat
                        yeswrite = .true.
                    else
                        if (signalnoise .and. idyn .eq. 0) then
                            yeswrite = .true.
                        else
                            yeswrite = .false.
                        end if
                    end if
                end if
                !  The wf has to be written when it is most accurate after reweight0.
                if (iskipdyn .gt. 1) then
                    if (ieskint .gt. 0 .and. (mod(i_main/nweight, iskipdyn) .eq. iskipdyn - nmore_force)) then
                        yeswrite12 = .true.
                    else
                        if (ieskint .eq. 0) then
                            !             when no dynamics write forces.dat
                            yeswrite12 = .true.
                        else
                            if (signalnoise .and. idyn .eq. 0) then
                                yeswrite12 = .true.
                            else
                                yeswrite12 = .false.
                            end if
                        end if
                    end if
                else
                    yeswrite12 = .true. ! no other possibility than writing each step.
                end if

                if (iesking .gt. 0 .or. (yeswrite .and. ieskint .gt. 0)) then
                    yesforce = .true.
                    yesdodet = .true.
                else
                    yesforce = .false.
                    yesdodet = yesdodet_nof
                end if

                if (detc_proj) then
                    firstmolc = 1
                    nmolfnc = nelorb_at
                else
                    firstmolc = firstmol
                    nmolfnc = nmolfn
                end if

                if (yeszagp .and. (.not. yesforce .or. itestr .ne. -5) .and. ipc .eq. 2 .and. yes_correct) then
                    yesdo_imag = .true.
                else
                    yesdo_imag = .false.
                end if

                if (yesdo_imag) then
                    elocb = 0.d0
                    logpsib = 0.d0
                    logpsib(2) = -1.d0
                    if (detc_proj) then
                        call computeb_global(js, contraction, contractionj&
                             &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc, nmolfnc, detmat_proj)
                    else
                        call computeb_global(js, contraction, contractionj&
                             &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc&
                             &, nmolfnc, detmat_c)
                    end if
                    if (contraction .ne. 0) then
                        call dcopy(iesupr_2, duprb, 1, psip, 1)
                        call transpsip_complex(psip, iesupr_2, iesup_c, iesuptrans&
                             &, multranspip, transpip, mu_cb)
                        call constrbr_complex(iesup, iesup_c, jbraiesup, psip                     &
                             &, econf(nwup + j), in1, 0)
                    else
                        psip(1:ipc*iesup_c) = 0.d0
                        call dcopy(iesup_c, duprb, 1, psip, ipc)
                        call constrbr_complex(iesup, iesup_c, jbraiesup, psip, econf(nwup + j), in1, 0)
                    end if
                    !          Save  the derivative of the imaginary part with respect to all parameters
                    do kk = 1, iesup
                        zagp_imag(kk) = econf(nwup + in1*(kk - 1) + j)
                    end do
                    call deallocate_computeb
                end if

                if (someparameter .or. yesforce) then
                    elocb = 0.d0
                    logpsib = 0.d0
                    logpsib(1) = 1.d0
                    if (detc_proj) then
                        call computeb_global(js, contraction, contractionj&
                             &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc, nmolfnc, detmat_proj)
                    else
                        call computeb_global(js, contraction, contractionj&
                             &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc, nmolfnc, detmat_c)
                    end if

                    !      Load all necessary to minimization+ forces if necessary
                    !      kelindlb,rionlb, celllogb + econf
                    if (iesd .ne. 0) then
                        call dcopy(iesd, vjb, 1, econf(nwdim + j), in1)
                    end if
                    if (iesfree .ne. 0) then
                        if (contractionj .eq. 0) then
                            if (yes_sparse) then
                                call constrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder, jasmatb&
                                      &, econf(nwfree + j), in1, 0)
                            else
                                call constrbra(iesfree, nnozeroj, jbraj, nozerojder, jasmatb&
                                     &, econf(nwfree + j), in1, 0)
                            end if
                        else
                            call constrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c          &
                                 &, jasmat_cb, econf(nwfree + j), in1, 0)
                        end if
                    end if
                    if (iessz) then
                        if (contractionj .eq. 0) then
                            call constrbra(iesinv, nnozeroj, jbrajsz, nozerojder          &
                                 &, jasmatszb, econf(nwinv + j), in1, 0)
                        else
                            call constrbra(iesinv, nnozeroj_c, jbrajsz, nozeroj_c         &
                                 &, jasmatsz_cb, econf(nwinv + j), in1, 0)
                        end if
                    end if
                    if (yeszj) then
                        call dcopy(size(vjurb), vjurb, 1, psip, 1)
                        !               if(rank.eq.0) then
                        !                write(6,*) ' Z Jastrow der wf ',size(vjurb)
                        !                do jj=1,size(vjurb)
                        !                write(6,*) jj,vjurb(jj)
                        !                enddo
                        !               endif
                        if (contractionj .ne. 0) then
                            ! npar3body_2 = npar3body/2 <= npar3body_c
                            call transpsip(psip, npar3body_2, npar3body_c, iesuptransj&
                                 &, multranspipj, transpipj, muj_cb)
                        end if
                        call constrbr(iesm, npar3body_c, jbraiesm, psip        &
                             &, econf(nwm2 + j), in1, 0)

                    elseif (iesm .ne. 0 .and. contractionj .ne. 0) then

                        ! in this case you optimize the contraction coefficients
                        do ii = 1, npar3body_c
                            psip(ii) = 0.d0
                            do jj = 1, multranspipj(ii)
                                psip(ii) = psip(ii) + muj_cb(transpipj(jj)%col(ii))
                            end do
                        end do

                        call constrbr(iesm, npar3body_c, jbraiesm, psip                   &
                             &, econf(nwm2 + j), in1, 0)
                    end if

                    if (yesdodet) then

                        ! derivative of the log wave function with respect to det lambda

                        !       alldet(iopt,indt,nelorb,nelup,neldo,ainv,winv        &
                        !    &,ainvupb,derl,der,contraction)
                        if (iessw .ne. 0) then
                            if (contraction .eq. 0) then
                                !       from normal to effective
                                !   if(rank.eq.0) write(6,*) ' Passi qui XV real-eff der '
                                if (allowed_averagek) call attach_phase2det(.false., detmatb)
                                call constraint_complex(iessw, detmatb, nnozero&
                                     &, nozerodet, psip, econf(nwsw + j), in1, jbradet)
                            else

                                if (yesmin .ne. 0 .and. detc_proj) then
                                    call projectder(nelorb_at, nelorb_c, nmolmatdo, nmolmat, projmat_c&
                                         &, nmolmat, detmat_cb, psip, yesmin, symmagp)
                                end if

                                if (allowed_averagek) call attach_phase2det(.false., detmat_cb)

                                call constraint_complex(iessw, detmat_cb, nnozero_c  &
                                     &, nozero_c, psip, econf(nwsw + j), in1, jbradet)
                            end if
                        end if

                        if (yeszagp) then
                            if (contraction .ne. 0) then
                                call dcopy(iesupr_2, duprb, 1, psip, 1)

                                call transpsip_complex(psip, iesupr_2, iesup_c, iesuptrans&
                                     &, multranspip, transpip, mu_cb)
                                call constrbr_complex(iesup, iesup_c, jbraiesup, psip                     &
                                     &, econf(nwup + j), in1, 0)
                            else
                                psip(1:ipc*iesup_c) = 0.d0
                                call dcopy(iesup_c, duprb, 1, psip, ipc)
                                call constrbr_complex(iesup, iesup_c, jbraiesup, psip, econf(nwup + j), in1, 0)
                            end if
                            if (yesdo_imag) then
                                !          Put  in the wf imaginary part der. with respect to Z  calculated previously. NB
                                !          the other  are consistent as far
                                !          as the contracted coefficient are concerned. So no need to distinguish them
                                !          in the loop below
                                do kk = 1, iesup, 2
                                    !                if(rank.eq.0) then
                                    !                write(6,*) ' before after =',econf(nwup+in1*kk+j),zagp_imag(kk)
                                    !                write(6,*) ' before after imag =',econf(nwup+in1*(kk-1)+j),-zagp_imag(kk+1)
                                    econf(nwup + in1*kk + j) = zagp_imag(kk)
                                    !                endif
                                end do
                            end if
                            !            call mpi_finalize(ierr)
                            !            stop

                        elseif (iesup .ne. 0 .and. contraction .ne. 0) then
                            ! in this case you optimize the contraction coefficients
                            if (ipc .eq. 1) then
                                do ii = 1, iesup_c
                                    psip(ii) = 0.d0
                                    do jj = 1, multranspip(ii)
                                        psip(ii) = psip(ii) + mu_cb(transpip(jj)%col(ii))
                                    end do
                                end do
                            else
                                do ii = 1, iesup_c
                                    psip(2*ii - 1:2*ii) = 0.d0
                                    do jj = 1, multranspip(ii)
                                        psip(2*ii - 1) = psip(2*ii - 1) + mu_cb(2*transpip(jj)%col(ii) - 1)
                                        psip(2*ii) = psip(2*ii) + mu_cb(2*transpip(jj)%col(ii))
                                    end do
                                end do
                            end if

                            call constrbr_complex(iesup, iesup_c, jbraiesup, psip               &
                                 &, econf(nwup + j), in1, 0)

                        end if

                    end if ! endif yesdodet
                    if (yesforce) then
                        celllogb = cellelb
                        sr2logb = sr2elb
                        allocate (rionlb(3, nion), kelindlb(3, nel))
                        rionlb = rionb
                        kelindlb = kelindb
                    end if
                    call deallocate_computeb
                end if ! yesparameter or yesforce

                if ((itestrr .eq. -4 .and. someparameter)&
                    & .or. yesforce&
                    & .or. lrdmc_der&
                    & .or. (someparameterdet .and. (yes_real .or. itest .eq. 1))) then
                    if (yesdo_imag) then
                        elocb = 0.d0
                        logpsib = 0.d0
                        elocb(2) = -1.d0
                        if (detc_proj) then
                            call computeb_global(js, contraction, contractionj&
                                 &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc, nmolfnc, detmat_proj)
                        else
                            call computeb_global(js, contraction, contractionj&
                                 &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc, nmolfnc, detmat_c)
                        end if
                        if (contraction .ne. 0) then
                            call dcopy(iesupr_2, duprb, 1, psip, 1)
                            call transpsip_complex(psip, iesupr_2, iesup_c, iesuptrans&
                                 &, multranspip, transpip, mu_cb)
                            call constrbr_complex(iesup, iesup_c, jbraiesup, psip                     &
                                 &, econfh(nwup + j), in1, 0)
                        else
                            psip(1:ipc*iesup_c) = 0.d0
                            call dcopy(iesup_c, duprb, 1, psip, ipc)
                            call constrbr_complex(iesup, iesup_c, jbraiesup, psip, econfh(nwup + j), in1, 0)
                        end if
                        !          Put  in the imaginary part
                        do kk = 1, iesup
                            zagp_imag(kk) = econfh(nwup + in1*(kk - 1) + j)
                        end do
                        call deallocate_computeb
                    end if

                    passicomp = .true.
                    elocb = 0.d0
                    logpsib = 0.d0
                    elocb(1) = 1.d0
                    if (detc_proj) then
                        call computeb_global(js, contraction, contractionj&
                             &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc, nmolfnc, detmat_proj)
                    else
                        call computeb_global(js, contraction, contractionj&
                             &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmolc, nmolfnc, detmat_c)
                    end if
                    !      Load all necessary to minimization+ forces if necessary
                    !      kelindb,rionb, cellelb + econf
                    if (iesd .ne. 0) then
                        call dcopy(iesd, vjb, 1, econfh(nwdim + j), in1)
                    end if
                    if (iesfree .ne. 0) then
                        if (contractionj .eq. 0) then
                            if (yes_sparse) then
                                call constrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder, jasmatb&
                                              &, econfh(nwfree + j), in1, 0)
                            else
                                call constrbra(iesfree, nnozeroj, jbraj, nozerojder, jasmatb&
                                     &, econfh(nwfree + j), in1, 0)
                            end if
                        else
                            call constrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c          &
                                 &, jasmat_cb, econfh(nwfree + j), in1, 0)
                        end if
                    end if
                    if (iessz) then
                        if (contractionj .eq. 0) then
                            call constrbra(iesinv, nnozeroj, jbrajsz, nozerojder          &
                                 &, jasmatszb, econfh(nwinv + j), in1, 0)
                        else
                            call constrbra(iesinv, nnozeroj_c, jbrajsz, nozeroj_c         &
                                 &, jasmatsz_cb, econfh(nwinv + j), in1, 0)
                        end if
                    end if
                    if (yeszj) then
                        call dcopy(size(vjurb), vjurb, 1, psip, 1)
                        !               if(rank.eq.0) then
                        !                write(6,*) ' Z Jastrow der en ',size(vjurb)
                        !                do jj=1,size(vjurb)
                        !                write(6,*) jj,vjurb(jj)
                        !                enddo
                        !               endif
                        if (contractionj .ne. 0) then
                            ! npar3body_2 = npar3body/2 <= npar3body_c
                            call transpsip(psip, npar3body_2, npar3body_c, iesuptransj&
                                 &, multranspipj, transpipj, muj_cb)
                        end if

                        call constrbr(iesm, npar3body_c, jbraiesm, psip        &
                             &, econfh(nwm2 + j), in1, 0)

                    elseif (iesm .ne. 0 .and. contractionj .ne. 0) then

                        ! in this case you optimize the contraction coefficients
                        do ii = 1, npar3body_c
                            psip(ii) = 0.d0
                            do jj = 1, multranspipj(ii)
                                psip(ii) = psip(ii) + muj_cb(transpipj(jj)%col(ii))
                            end do
                        end do

                        call constrbr(iesm, npar3body_c, jbraiesm, psip                   &
                             &, econfh(nwm2 + j), in1, 0)
                    end if

                    if (yesdodet) then

                        ! derivative of the log wave function with respect to det lambda

                        !       alldet(iopt,indt,nelorb,nelup,neldo,ainv,winv        &
                        !    &,ainvupb,derl,der,contraction)
                        if (iessw .ne. 0) then
                            if (contraction .eq. 0) then
                                !       from normal to effective
                                !   if(rank.eq.0) write(6,*) ' Passi qui XVII real-eff der'
                                if (allowed_averagek) call attach_phase2det(.false., detmatb)
                                call constraint_complex(iessw, detmatb, nnozero&
                                     &, nozerodet, psip, econfh(nwsw + j), in1, jbradet)
                            else

                                if (yesmin .ne. 0 .and. detc_proj) then
                                    call projectder(nelorb_at, nelorb_c, nmolmatdo, nmolmat, projmat_c&
                                         &, nmolmat, detmat_cb, psip, yesmin, symmagp)
                                end if
                                !       from normal to effective
                                !   if(rank.eq.0) write(6,*) ' Passi qui XVIII real-eff der'
                                if (allowed_averagek) call attach_phase2det(.false., detmat_cb)
                                call constraint_complex(iessw, detmat_cb, nnozero_c  &
                                     &, nozero_c, psip, econfh(nwsw + j), in1, jbradet)

                            end if
                        end if

                        if (yeszagp) then
                            if (contraction .ne. 0) then
                                call dcopy(iesupr_2, duprb, 1, psip, 1)
                                call transpsip_complex(psip, iesupr_2, iesup_c, iesuptrans&
                                     &, multranspip, transpip, mu_cb)
                                call constrbr_complex(iesup, iesup_c, jbraiesup, psip                     &
                                     &, econfh(nwup + j), in1, 0)
                            else
                                psip(1:ipc*iesup_c) = 0.d0
                                call dcopy(iesup_c, duprb, 1, psip, ipc)
                                call constrbr_complex(iesup, iesup_c, jbraiesup, psip, econfh(nwup + j), in1, 0)
                            end if
                            if (yesdo_imag) then
                                !          Put  in the imaginary part calculated previously
                                do kk = 1, iesup, 2
                                    econfh(nwup + in1*kk + j) = zagp_imag(kk)
                                end do
                            end if
                        elseif (iesup .ne. 0 .and. contraction .ne. 0) then
                            ! in this case you optimize the contraction coefficients

                            if (ipc .eq. 1) then
                                do ii = 1, iesup_c
                                    psip(ii) = 0.d0
                                    do jj = 1, multranspip(ii)
                                        psip(ii) = psip(ii) + mu_cb(transpip(jj)%col(ii))
                                    end do
                                end do
                            else
                                do ii = 1, iesup_c
                                    psip(2*ii - 1:2*ii) = 0.d0
                                    do jj = 1, multranspip(ii)
                                        psip(2*ii - 1) = psip(2*ii - 1) + mu_cb(2*transpip(jj)%col(ii) - 1)
                                        psip(2*ii) = psip(2*ii) + mu_cb(2*transpip(jj)%col(ii))
                                    end do
                                end do
                            end if

                            call constrbr_complex(iesup, iesup_c, jbraiesup, psip               &
                                 &, econfh(nwup + j), in1, 0)

                        end if

                    end if ! endif yesdodet

                end if

                if (yesforce) then
                    allocate (forcedw(3, nion), pulaydw(3, nion))
                    forcedw = 0.d0
                    pulaydw = 0.d0
                    call updatedwarp(nel, nion, kelind, rion, rmu, r, kelindb, kelindlb&
                         &, rionb, rionlb, cellelb, celllogb, forcedw, pulaydw&
                         &, iespbc, warp, powerwarp, atom_number, warpmat, niong)

                    if (nbra_cyrus .gt. 0) then
                        ! Update the sum registers by definition at ncyrus+1,
                        ! ncyrus+2 contains the previous gradients of local energy
                        if (lrdmc_der) then

                            !     store the present gradients for the next step (we need the previous  gradients )
                            queue_cyrus(:, 1:nion, nbra_cyrus + 2, j) = forcedw(:, 1:nion)
                            queue_cyrus(:, nion + 1:2*nion, nbra_cyrus + 2, j) = pulaydw(:, 1:nion)
                            queue_cyrus(:, 2*nion + 1, nbra_cyrus + 2, j) = cellelb(:)
                            queue_cyrus(:, 2*nion + 2, nbra_cyrus + 2, j) = celllogb(:)

                            forcedw(:, 1:nion) = queue_cyrus(:, 1:nion, first_cyrus(j), j)
                            cellelb(:) = queue_cyrus(:, 2*nion + 1, first_cyrus(j), j)
                            pulaydw(:, 1:nion) = (1.d0 - weight_moroni)&
                                &*queue_cyrus(:, nion + 1:2*nion, first_cyrus(j), j)&
                                & + weight_moroni*pulaydw(:, 1:nion)
                            celllogb(:) = (1.d0 - weight_moroni)&
                                &*queue_cyrus(:, 2*nion + 2, first_cyrus(j), j)&
                               & + weight_moroni*celllogb(:)

                            queue_cyrus(:, 1:2*nion + 2, first_cyrus(j), j) = queue_cyrus(:, 1:2*nion + 2, nbra_cyrus + 2, j)

                            !               Update the pointer of the queue
                            first_cyrus(j) = mod(first_cyrus(j), nbra_cyrus) + 1

                        else
                            !     store the present gradients for the next step (we need the previous  gradients )
                            queue_cyrus(:, nion + 1:2*nion, nbra_cyrus + 2, j) = pulaydw(:, 1:nion)
                            queue_cyrus(:, 2*nion + 2, nbra_cyrus + 2, j) = celllogb(:)

                            !pulaydw(:, 1:nion) = queue_cyrus(:, nion+1:2*nion, first_cyrus(j), j)
                            !celllogb(:) = queue_cyrus(:, 2*nion + 2, first_cyrus(j), j)

                            !   end our Pulay correction.

                            queue_cyrus(:, 1:nion, nbra_cyrus + 1, j) = queue_cyrus(:, 1:nion, nbra_cyrus + 1, j)&
    & - queue_cyrus(:, 1:nion, first_cyrus(j), j) - (forcedw(:, 1:nion)&
    & + queue_cyrus(:, 1:nion, nbra_cyrus + 2, j))/2.d0*t_cyrus(j)
                            !        the same for pressures
                            queue_cyrus(:, 2*nion + 1, nbra_cyrus + 1, j) = queue_cyrus(:, 2*nion + 1, nbra_cyrus + 1, j)&
    & - queue_cyrus(:, 2*nion + 1, first_cyrus(j), j) - (cellelb(:)&
    & + queue_cyrus(:, 2*nion + 1, nbra_cyrus + 2, j))/2.d0*t_cyrus(j)

                            !   Replace the new gradients in the position of the oldest in the queue
                            queue_cyrus(:, 1:nion, first_cyrus(j), j) = -(forcedw(:, 1:nion)&
                               & + queue_cyrus(:, 1:nion, nbra_cyrus + 2, j))/2.d0 * t_cyrus(j)

                            queue_cyrus(:, 2*nion + 1, first_cyrus(j), j) = -(cellelb(:)&
                               & + queue_cyrus(:, 2*nion + 1, nbra_cyrus + 2, j))/2.d0*t_cyrus(j)
                            !               Update the pointer
                            first_cyrus(j) = mod(first_cyrus(j), nbra_cyrus) + 1
                            !     store the present gradients for the next step (we need the previous  gradients )
                            queue_cyrus(:, 1:nion, nbra_cyrus + 2, j) = forcedw(:, 1:nion)
                            queue_cyrus(:, 2*nion + 1, nbra_cyrus + 2, j) = cellelb(:)

                            !     Update the pulay

                            pulaydw(:, 1:nion) = (1.d0 - weight_moroni)*queue_cyrus(:, nion + 1:2*nion, first_cyrus(j), j)&
                               & + weight_moroni*pulaydw(:, 1:nion) + queue_cyrus(:, 1:nion, nbra_cyrus + 1, j)
                            celllogb(:) = (1.d0 - weight_moroni)*queue_cyrus(:, 2*nion + 2, first_cyrus(j), j)&
                               & + weight_moroni*celllogb(:) + queue_cyrus(:, nion + 1, nbra_cyrus + 1, j)

                            queue_cyrus(:, nion + 1:2*nion, first_cyrus(j), j) =&
                                & queue_cyrus(:, nion + 1:2*nion, nbra_cyrus + 2, j)
                            queue_cyrus(:, 2*nion + 2, first_cyrus(j), j) =&
                                & queue_cyrus(:, 2*nion + 2, nbra_cyrus + 2, j)
                        end if
                    end if

                    ! by E. Coccia (3/7/11): updating forces for capping atoms
                    if (link_atom) then
                        call force_capping(forcedw)
                        call force_capping(pulaydw)
                    end if

                    kkf = 0

                    do kk = 1, ieskinr
                        ind = abs(ion_table(kk)%mult)
                        if (ion_table(kk)%mult .gt. 0) then
                            totpulay = 0.d0
                            totforce = 0.d0
                            k_ion = abs(ion_table(kk)%ion(1))

                            if (atom_number(k_ion) .gt. 0) then
                                do jj = 1, ind
                                    i_ion = ion_table(kk)%comp(jj)
                                    k_ion = abs(ion_table(kk)%ion(jj))
                                    if (ion_table(kk)%ion(jj) .gt. 0) then
                                        totpulay = totpulay + pulaydw(i_ion, k_ion)
                                        totforce = totforce + forcedw(i_ion, k_ion)
                                    else
                                        totpulay = totpulay - pulaydw(i_ion, k_ion)
                                        totforce = totforce - forcedw(i_ion, k_ion)
                                    end if
                                end do
                                kkf = kkf + 1

                                econf(nwkin + j + (kkf - 1)*in1*skipforce) = totpulay*scalepulay
                                econfion(in1*(kkf - 1) + j) = totforce

                                if (itestr .ne. -5) then
                                    econf(nwkin + j + in1*((kkf - 1)*skipforce + 1)) = econfion(in1*(kkf - 1) + j)
                                    econf(nwkin + j + in1*((kkf - 1)*skipforce + 2)) = &
                                         &econf(nwkin + j + in1*(kkf - 1)*skipforce)*enert(1, j)
                                end if
                                !
                            end if ! endif kkf

                        end if ! if it is allowed

                        ! kk=1,ieskinr
                    end do

                    if (iespbc) then

                        if (fixc) then
                            !  Surface pressure
                            derEV(j) = (cellelb(1)/cellscale(2) + cellelb(2)/cellscale(1)&
                                 &)/(2.d0*unit_volume)

                            p_pulay(j) = (celllogb(1)/cellscale(2) + celllogb(2)/cellscale(1)&
                                 &)/(2.d0*unit_volume)*scalepulay

                        else

                            derEV(j) = (cellelb(1)/cellscale(2)/cellscale(3) + cellelb(2)/cellscale(1)&
                                 & /cellscale(3) + cellelb(3)/cellscale(1)/cellscale(2))/(3.d0*unit_volume)

                            p_pulay(j) = (celllogb(1)/cellscale(2)/cellscale(3) + celllogb(2)/cellscale(1)&
                                 & /cellscale(3) + celllogb(3)/cellscale(1)/cellscale(2))/(3.d0*unit_volume)*scalepulay

                        end if

                    end if
                    !        write(6,*) ' Volume Derivative local energy DWarp ',i,derEV(j)
                    !        write(6,*) ' Volume Derivative logpsi  DWarp ',i,p_pulay(j)

                    if (ieskint .gt. ieskinr_pos) then

                        if (ieskint .eq. ieskinrp .and. yespress) then

                            kk = ieskinr_pos + 1 - iesking
                            econf(nwkin + j + in1*((kk - 1)*skipforce)) = p_pulay(j)
                            econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = derEV(j)
                            econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = p_pulay(j)*enert(1, j)

                        elseif (ieskint .eq. ieskinrp + 1) then

                            !     derivative with respect to Lx at fixed Ly/Lx Lz/Lx
                            kk = ieskinr_pos + 1 - iesking
                            deriv1 = (derEV(j) + pressfixed)*3.d0*cellscale(2)*cellscale(3)
                            pulay1 = p_pulay(j)*3.d0*cellscale(2)*cellscale(3)
                            econfion(in1*(kk - 1) + j) = deriv1
                            econf(nwkin + j + (kk - 1)*in1*skipforce) = pulay1

                            if (itestr .ne. -5) then
                                !     save the pressure
                                kk = ieskinr_pos + 1 - iesking

                                econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = deriv1
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = pulay1*enert(1, j)

                                if (yespress) then
                                    kk = ieskinr_pos + 2 - iesking
                                    econf(nwkin + j + in1*((kk - 1)*skipforce)) = p_pulay(j)
                                    econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = derEV(j)
                                    econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = p_pulay(j)*enert(1, j)
                                end if

                            end if

                        elseif (ieskint .eq. ieskinrp + 2) then

                            deriv1 = cellelb(2) - cellscale(1)/cellscale(2)*cellelb(1)

                            pulay1 = celllogb(2) - cellscale(1)/cellscale(2)*celllogb(1)

                            deriv2 = cellelb(3) - cellscale(1)/cellscale(3)*cellelb(1)

                            pulay2 = celllogb(3) - cellscale(1)/cellscale(3)*celllogb(1)

                            pulay1 = pulay1*scalepulay
                            pulay2 = pulay2*scalepulay

                            !      load econfion/econf
                            kk = ieskinr_pos + 1 - iesking
                            econfion(in1*(kk - 1) + j) = deriv1
                            econf(nwkin + j + (kk - 1)*in1*skipforce) = pulay1

                            kk = ieskinr_pos + 2 - iesking
                            econfion(in1*(kk - 1) + j) = deriv2

                            econf(nwkin + j + (kk - 1)*in1*skipforce) = pulay2

                            if (itestr .ne. -5) then ! forces for VMC or DMC, GFMC

                                !     derivative energy with respect to b with fixed volume=axbxc

                                kk = ieskinr_pos + 1 - iesking

                                econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = deriv1
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = pulay1*enert(1, j)

                                kk = ieskinr_pos + 2 - iesking

                                !     derivative energy with respect to c constant axbxc
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = deriv2
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = pulay2*enert(1, j)

                                if (yespress) then
                                    kk = ieskinr_pos - iesking + 3
                                    econf(nwkin + j + in1*((kk - 1)*skipforce)) = p_pulay(j)
                                    econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = derEV(j)
                                    econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = p_pulay(j)*enert(1, j)
                                end if

                            end if

                        elseif (ieskint .eq. ieskinrp + 3) then

                            deriv1 = cellelb(1) + pressfixed*cellscale(2)*cellscale(3)

                            pulay1 = celllogb(1)

                            deriv2 = cellelb(2) + pressfixed*cellscale(1)*cellscale(3)

                            pulay2 = celllogb(2)

                            deriv3 = cellelb(3) + pressfixed*cellscale(1)*cellscale(2)

                            pulay3 = celllogb(3)

                            pulay1 = pulay1*scalepulay
                            pulay2 = pulay2*scalepulay
                            pulay3 = pulay3*scalepulay

                            !      load econfion/econf
                            kk = ieskinr_pos + 1 - iesking
                            econfion(in1*(kk - 1) + j) = deriv1
                            econf(nwkin + j + (kk - 1)*in1*skipforce) = pulay1

                            kk = ieskinr_pos + 2 - iesking
                            econfion(in1*(kk - 1) + j) = deriv2

                            econf(nwkin + j + (kk - 1)*in1*skipforce) = pulay2

                            kk = ieskinr_pos + 3 - iesking
                            econfion(in1*(kk - 1) + j) = deriv3

                            econf(nwkin + j + (kk - 1)*in1*skipforce) = pulay3

                            ! forces for VMC or DMC, GF
                            if (itestr .ne. -5) then

                                kk = ieskinr_pos + 1 - iesking

                                econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = deriv1
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = pulay1*enert(1, j)

                                kk = ieskinr_pos + 2 - iesking
                                !     derivative energy with respect to c constant axbxc
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = deriv2
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = pulay2*enert(1, j)

                                kk = ieskinr_pos + 3 - iesking
                                !     derivative energy with respect to c constant axbxc
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = deriv3
                                econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = pulay3*enert(1, j)

                                if (yespress) then
                                    kk = ieskinr_pos - iesking + 4
                                    econf(nwkin + j + in1*((kk - 1)*skipforce)) = p_pulay(j)
                                    econf(nwkin + j + in1*((kk - 1)*skipforce + 1)) = derEV(j)
                                    econf(nwkin + j + in1*((kk - 1)*skipforce + 2)) = p_pulay(j)*enert(1, j)
                                end if

                            end if

                        end if ! endif ieskin eq ieskinr+

                    end if ! endif ieskin> ieskinr

                end if ! ieskin.gt.0

                !     END CALCULATIONS OF IONIC FORCES
                if (passicomp) call deallocate_computeb

                ! energy and correlations functions after makeallmeas
                if (itest .eq. 2) then
                    if (.not. yes_complex) then
                        call comp_econf(j, js, econf, enert, diffkin, vpot, vcut, voffpseudo, table)
                        if (isfix .ne. 0) econf(j + nwfix) = enertrue(j)**2
                    else
                        call comp_econf(j, js, econf, enert, diffkin, vpot, vcut, voffpseudo, tabler)
                        if (isfix .ne. 0) econf(j + nwfix) = ener_c*dconjg(ener_c)
                    end if

                elseif (fncont) then
                    ttry = tbra
                    call upgradcont(gradpsibar(1, indksj), gradpsi(1, indksj)                  &
                         &, indt, nelup, neldo, winvup(indwupj), winvdo(indwdoj)                 &
                         &, tabpip(indtabj), ttry, gradtotbar(j), gradtot(j), kel(1, indkj)       &
                         &, dist(indksij), rion, nion, zetar, LBox)
                else
                    call diffus(nel, indt, gamma, ivic(1, 1, indksj), table(indtj)&
                         &, diffuse(j), istart)
                end if

                !-----------------------Berry phase---------------------------!

                !   calculation derivatives with respect to iesking only case
                if (iesking .gt. 0) then
                    kkf = 0
                    do kk = 1, ieskinr
                        ind = abs(ion_table(kk)%mult)
                        if (ion_table(kk)%mult .gt. 0) then
                            totpulay = 0.d0
                            totforce = 0.d0
                            k_ion = abs(ion_table(kk)%ion(1))

                            if (atom_number(k_ion) .lt. 0) then
                                do jj = 1, ind
                                    k_ion = abs(ion_table(kk)%ion(jj))
                                    i_ion = ion_table(kk)%comp(jj)
                                    if (ion_table(kk)%ion(jj) .gt. 0) then
                                        totpulay = totpulay + pulaydw(i_ion, k_ion)
                                        totforce = totforce + forcedw(i_ion, k_ion)
                                    else
                                        totpulay = totpulay - pulaydw(i_ion, k_ion)
                                        totforce = totforce - forcedw(i_ion, k_ion)
                                    end if
                                end do
                                kkf = kkf + 1

                                econf(nwking + j + (kkf - 1)*in1) = totpulay
                                econfh(nwking + in1*(kkf - 1) + j) = totforce

                            end if ! endif kkf

                        end if ! endif allowed component

                        ! kk=1,ieskinr
                    end do
                end if

                ! by E. Coccia (20/12/11): write electronic random walk
                if (rank .eq. 0 .and. write_rwalk) then
                    call write_traj(nel, nelup, kel(1, indkj), enert(1, j))
                end if

                if (allocated(pulaydw)) deallocate (pulaydw, forcedw, kelindlb, rionlb)

                if (nrep .ne. 0) then
                    !  Correlation functions times local energy
                    if (yes_correct) then

                        if (lrdmc_der .and. .not. lrdmc_nonodes) then
                            do jj = 1, ndimj
                                econf(nwrep + in1*(jj - 1) + j) = econfh(j + in1*(jj - 1 + iese))
                                econf(j + in1*(jj - 1 + iese)) = 0.d0
                            end do
                            do jj = ndimjp, nrep, 2
                                econf(nwrep + in1*(jj - 1) + j) = econfh(j + in1*(jj - 1 + iese))
                                econf(j + in1*(jj - 1 + iese)) = 0.d0
                                econf(nwrep + in1*jj + j) = econfh(j + in1*(jj + iese))
                                econf(j + in1*(jj + iese)) = 0.d0
                            end do
                        else
                            do jj = 1, ndimj
                                econf(nwrep + in1*(jj - 1) + j) = econf(j + in1*(jj - 1 + iese))*enert(1, j)
                            end do
                            do jj = ndimjp, nrep, 2
                                ! e_L (O_R-i O_I) | NB O_I = - econf
                                econf(nwrep + in1*(jj - 1) + j) = econf(j + in1*(jj - 1 + iese))*enert(1, j) - &
                                     &             econf(j + in1*(jj + iese))*enert(2, j)
                                ! By Cauchy relations O_Rnew=-OI O_Inew=O_R
                                ! e_L (O_Rnew-i O_Inew)=e_L (-O_I-i O_R)  | NB O_I = - econf
                                econf(nwrep + in1*jj + j) = econf(j + in1*(jj - 1 + iese))*enert(2, j) + &
                                     &             econf(j + in1*(jj + iese))*enert(1, j)
                            end do
                        end if
                    elseif (yes_real) then
                        !            do jj=1,ndimj
                        !               econf(nwrep+in1*(jj-1)+j)=econf(j+in1*(jj-1+iese))*enert(1,j)
                        !            enddo
                        !            do jj=ndimjp,nrep
                        do jj = 1, nrep
                            econf(nwrep + in1*(jj - 1) + j) = econf(j + in1*(jj - 1 + iese))*enert(1, j) + &
                                 &scaleeloc*econfh(j + in1*(jj - 1 + iese))
                        end do
                    else
                        if (lrdmc_der .and. .not. lrdmc_nonodes) then
                            do jj = 1, nrep
                                econf(nwrep + in1*(jj - 1) + j) = econfh(j + in1*(jj - 1 + iese))
                                econf(j + in1*(jj - 1 + iese)) = 0.d0
                            end do
                        else
                            do jj = 1, nrep
                                econf(nwrep + in1*(jj - 1) + j) = econf(j + in1*(jj - 1 + iese))*enert(1, j)
                            end do
                        end if
                    end if
                end if
            end if ! fine if psi(j)>0
            !      call mpi_finalize(ierr)
            !      stop

        end do ! end do for the walkers
        ! end parallel scratch
        return
    end subroutine makeallmeas

    subroutine updatewfopt

        !By E. Coccia (13/12/11)
        use van_der_waals, only: sum_pot

        implicit none
        real*8 drand1, enercont, jacobian, mapping, time0
        !       This subroutine updates all necessary information concerning
        !       the wavefunction when doing optimization or dynamics.
        timep = cclock()

        ! optimization #########################################################
        if (i_main .le. inext + ibinit - nweight*perbin) wtotf = 0.d0

        if (i_main .eq. inext) then
            !      save the alpha before changing it

            if (pippo .eq. ibinit + 1) then
                pippoc = 1
                call dscalzero((1 + ipc)*nbindim, 0.d0, efenergy, 1)
                call dscalzero(nbindim*ndimpdim*2, 0.d0, efp, 1)
                call dscalzero(nbindim*ieskindim*3, 0.d0, ef, 1)
                if (LBox .gt. 0.d0) call dscalzero(3*nbindim, 0.d0, efpress, 1)
                srforce = 0.d0
                srforcew = 0.d0
                wforce = 0.d0
                wforcew = 0.d0
            end if

            if (LBox .gt. 0.d0) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    efpress(1, indpippo) = derEV(cpippo - istm)*factorsr(cpippo - istm)*wcort   &
                         & + efpress(1, indpippo)
                    efpress(2, indpippo) = p_pulay(cpippo - istm)*factorsr(cpippo - istm)  &
                                          & *wcort + efpress(2, indpippo)
                    efpress(3, indpippo) = p_pulay(cpippo - istm)*factorsr(cpippo - istm)  &
                                          & *wcort*enert(1, cpippo - istm) + efpress(3, indpippo)
                end do
            end if

            if (ipc .eq. 2) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    efenergy(1:2, indpippo) = enert(1:2, cpippo - istm)                             &
                                             & *factorsr(cpippo - istm)*wcort + efenergy(1:2, indpippo)
                    efenergy(3, indpippo) = factorsr(cpippo - istm)*wcort + efenergy(3, indpippo)
                end do
            else
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    efenergy(1, indpippo) = enert(1, cpippo - istm)                             &
                                           & *factorsr(cpippo - istm)*wcort + efenergy(1, indpippo)
                    efenergy(2, indpippo) = factorsr(cpippo - istm)*wcort + efenergy(2, indpippo)
                end do
            end if

            !************* PULAY PRESSURE ****************
            if (LBox .gt. 0.d0) then
                call bootpress(dimfk, efenergy, efpress, derEVp, errEVp             &
                     &, derEVnopulay, errEVnopulay, LBox, omega, ris, rank, nproc)
            end if

            if (ndimp .ge. 1) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    if (yes_correct) then

                        if (lrdmc_der .and. .not. lrdmc_nonodes) then
                            do kk = 1, ndimj
                                efp(kk, 1, indpippo) = econfh((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do
                            do kk = ndimjp, ndimp, 2
                                efp(kk, 1, indpippo) = econfh((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk + 1, 1, indpippo) = econfh(kk*in1 + cpippo - istm)              &
                                                          & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 1, indpippo)
                                efp(kk, 2, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                                efp(kk + 1, 2, indpippo) = econf(kk*in1 + cpippo - istm)              &
                                                          & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 2, indpippo)
                            end do

                        else
                            do kk = 1, ndimj
                                efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = enert(1, cpippo - istm)                                  &
                                                      & *econf(cpippo - istm + (kk - 1)*in1)                                    &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do
                            do kk = ndimjp, ndimp, 2
                                efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk + 1, 1, indpippo) = econf(kk*in1 + cpippo - istm)              &
                                                          & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 1, indpippo)
                                efp(kk, 2, indpippo) = (enert(1, cpippo - istm)*econf(cpippo - istm + (kk - 1)*in1)&
                                     & - enert(2, cpippo - istm)*econf(cpippo - istm + kk*in1))&
                                     & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                                efp(kk + 1, 2, indpippo) = (enert(2, cpippo - istm)*econf(cpippo - istm + (kk - 1)*in1)&
                                     & + enert(1, cpippo - istm)*econf(cpippo - istm + kk*in1))&
                                     & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 2, indpippo)
                            end do
                        end if
                    elseif (yes_real) then
                        !                  do kk=1,ndimj
                        !                     efp(kk,1,indpippo)=econf((kk-1)*in1+cpippo-istm)              &
                        !                          &*factorsr(cpippo-istm)*wcort+efp(kk,1,indpippo)
                        !                     efp(kk,2,indpippo)=enert(1,cpippo-istm)                                  &
                        !                          &*econf(cpippo-istm+(kk-1)*in1)                                    &
                        !                          &*factorsr(cpippo-istm)*wcort+efp(kk,2,indpippo)
                        !                  enddo

                        !                  do kk=ndimjp,ndimp
                        do kk = 1, ndimp

                            efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                  & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                            efp(kk, 2, indpippo) = (enert(1, cpippo - istm)&
                                 & *econf(cpippo - istm + (kk - 1)*in1) + &
                                 & scaleeloc*econfh(cpippo - istm + (kk - 1)*in1))&
                                 & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                        end do
                    else

                        if (lrdmc_der .and. .not. lrdmc_nonodes) then
                            do kk = 1, ndimp
                                efp(kk, 1, indpippo) = econfh((kk - 1)*in1 + cpippo - istm)  &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = econf((kk - 1)*in1 + cpippo - istm)  &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do
                        else
                            do kk = 1, ndimp
                                efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = enert(1, cpippo - istm)                                  &
                                                      & *econf(cpippo - istm + (kk - 1)*in1)                                    &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do
                        end if
                    end if
                end do
            end if ! ndimp.ge.1

            if (ieskin .ge. 1) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    do kk = 1, ieskin
                        ! der energy
                        ef(kk, 1, indpippo) = econfion((kk - 1)*in1 + cpippo - istm)             &
                                             & *factorsr(cpippo - istm)*wcort + ef(kk, 1, indpippo)

                        ! lor Or
                        ef(kk, 2, indpippo) = econf(nwkin + (kk - 1)*in1 + cpippo - istm)          &
                                             & *factorsr(cpippo - istm)*wcort + ef(kk, 2, indpippo)
                        ! log Or*energy
                        ef(kk, 3, indpippo) = enert(1, cpippo - istm)                                &
                                             & *econf(nwkin + cpippo - istm + (kk - 1)*in1)                           &
                                             & *factorsr(cpippo - istm)*wcort + ef(kk, 3, indpippo)
                        !                      ef(kk,4,indpippo)=econf(nwkin+cpippo-istm+(kk-1)*in1)**2       &
                        !                           &   *factorsr(cpippo-istm)*wcort+ef(kk,4,indpippo)
                    end do
                end do ! end do for kk=1,ieskin
            end if ! endif for ieskin.ge.1

            pippo = 1
            pippoc = 1

            if (mod(i_main/nweight, iskipdyn) .eq. 0) then
                idynu = idyn
            else
                idynu = -idyn
            end if

            if (idynu .gt. 0 .or. (idyn .eq. 0 .and. signalnoise .and. ieskin .gt. 0)) then
                yesupdate_ion = .true.
            else
                yesupdate_ion = .false.
            end if

            if (ndimp .gt. 0 .or. ieskint .gt. 0) then
#ifdef UNREL
#ifdef PARALLEL
                call mpi_barrier(MPI_COMM_WORLD, ierr) ! for unreliable networks
#endif
!$omp barrier
#endif

                if (yesupdate_ion) then

                    call bootforcecov(nbin, efenergy, ndimp, ef, ieskin, force&
                         &, err, psip, fk, dimfk, okav, rank, commrep_mpi, nprocrep, commcov_mpi, nproccov)

#ifdef PARALLEL
                    if (yesquantum) then
                        !            Calculation virial before correcting forces.
                        !   Add quantum corrections to energy
                        !  Calculation centroid coordinates
                        call mpi_allreduce(rion, fbead, nion*3                      &
                             &, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                        fbead = fbead/nproc
                        fbead = fbead - rion
                        if (iespbc) call Apply_periodic(fbead, nion, cellscale, dt4, mass_ion, cost)
                        cost = 0.d0
                        !          sumbeta=0.d0
                        do jj = 1, nion
                            do kk = 1, 3
                                ind = ndimp + kk + (jj - 1)*3
                                cost = cost + 0.5d0*fbead(kk, jj)*force(ind)
                                !           sumbeta=sumbeta-fbead(kk,jj)*force(ind)
                            end do
                        end do

                        call mpi_allreduce(cost, ekinq, 1                      &
                             &  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                        ekinq = ekinq/nproc + 1.5d0*nion*temp/nbead
                        !                  cost=sumbeta
                        !                  call mpi_allreduce(cost,sumbeta,1                      &
                        !                       &  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                        !                  sumbeta=sumbeta/3.d0/dble(nion)/dble(nproc)
                        !                  tmes=sumbeta
                    end if
#endif

                end if

                !   The averages of forces have to be done within each bead and only after
                !   can be averaged again
                call bootparameter(nbin, efenergy, efp, ndimj, ndimp, force, err, psip, sigma_true&
                     &, ener_true, fk, dimfk, okav, weightall, 1, ndimp, rank, commrep_mpi, nprocrep&
                     &, commrep_mpi, nprocrep, nproc)

#ifdef PARALLEL

                if (yesquantum) then

                    if (nprocsr .gt. nprocrep) then
                        ! In this way we sum the SR matrix after MC average each independent of the other
                        weightall = weightall*dble(nprocsr)/dble(nprocrep)
                        call reduce_base_real(ndimp, force, commsr_mpi, -1)
                        !         The average force.
                        force(1:ndimp) = force(1:ndimp)/nprocsr
                        call reduce_base_real(ndimp, err, commsr_mpi, -1)
                        !         The average error
                        err(1:ndimp) = err(1:ndimp)/nprocsr/dsqrt(dble(nprocsr/nprocrep))
                        !
                        fk(:, 1:ndimp) = fk(:, 1:ndimp)/dble(nprocsr/nprocrep)
                    end if

                end if

                if (kaverage .and. .not. decoupled_run) then

                    !           Here each k has its own weight proportinal to wkp.
                    !           NB: k-average presently works only with yesavsr=.true. option
                    if (nprocsr .gt. nprocrep) then
                        !            Total weight (not averaged over k-points)
                        weightall = weightall/wkp(ikpoint)
                        ener_true = ener_true*wkp(ikpoint)
                        call reduce_base_real(1, ener_true, commsr_mpi, -1)
                        ener_true = ener_true/nprocrep
                        sigma_true = (sigma_true*wkp(ikpoint))**2
                        call reduce_base_real(1, sigma_true, commsr_mpi, -1)
                        sigma_true = dsqrt(sigma_true/dble(nprocrep))
                        !            The average force
                        force(:) = force(:)*wkp(ikpoint)
                        call reduce_base_real(size(force), force, commsr_mpi, -1)
                        force(:) = force(:)/dble(nprocrep)
                        !            The average error
                        err(:) = (err(:)*wkp(ikpoint))**2
                        call reduce_base_real(size(err), err, commsr_mpi, -1)
                        err(:) = dsqrt(err(:)/dble(nprocrep))
                        fk(:, :) = fk(:, :)*wkp(ikpoint)
                    end if

                end if
#endif

#ifdef PARALLEL
                if (yesquantum .and. nbead .gt. 1 .and. idynu .gt. 0) then

                    call mpi_allgather(rion, 3*nion, MPI_DOUBLE_PRECISION, rionall&
                         &, 3*nion, MPI_DOUBLE_PRECISION, commcolrep_mpi, ierr)

                    !    modify dek depending on rank value

                    jj = rankcolrep + 1
                    ekinqp = 0.d0
                    if (iespbc) then
                        if (jj .gt. 1 .and. jj .lt. nbead) then
                            rionsav(:, :) = rionall(:, :, jj + 1) - rionall(:, :, jj)
                            if (yesperiodize) then
                                call Apply_periodic(rionsav, nion, cellscale, dt4, mass_ion, cost)
                            else
                                call ApplyPBC(rionsav, nion)
                            end if
                            fbead(:, :) = rionall(:, :, jj - 1) - rionall(:, :, jj)
                            if (yesperiodize) then
                                call Apply_periodic(fbead, nion, cellscale, dt4, mass_ion, cost)
                                ekinqp = ekinqp + cost
                            else
                                call ApplyPBC(fbead, nion)
                                ekinqp = ekinqp + sum(mass_ion(:, :)*fbead(:, :)**2)
                            end if
                            fbead(:, :) = fbead(:, :) + rionsav(:, :)
                            !                  fbead(:,:)=rionall(:,:,jj+1) - 2*rionall(:,:,jj) + rionall(:,:,jj-1)
                        elseif (jj .eq. 1) then
                            rionsav(:, :) = rionall(:, :, jj + 1) - rionall(:, :, jj)
                            if (yesperiodize) then
                                call Apply_periodic(rionsav, nion, cellscale, dt4, mass_ion, cost)
                            else
                                call ApplyPBC(rionsav, nion)
                            end if
                            fbead(:, :) = rionall(:, :, nbead) - rionall(:, :, jj)
                            if (yesperiodize) then
                                call Apply_periodic(fbead, nion, cellscale, dt4, mass_ion, cost)
                                ekinqp = ekinqp + cost
                            else
                                call ApplyPBC(fbead, nion)
                                ekinqp = ekinqp + sum(mass_ion(:, :)*fbead(:, :)**2)
                            end if
                            fbead(:, :) = fbead(:, :) + rionsav(:, :)
                            !                  fbead(:,:) =rionall(:,:,jj+1) - 2*rionall(:,:,jj) + rionall(:,:,nbead)
                        else
                            rionsav(:, :) = rionall(:, :, 1) - rionall(:, :, jj)
                            if (yesperiodize) then
                                call Apply_periodic(rionsav, nion, cellscale, dt4, mass_ion, cost)
                            else
                                call ApplyPBC(rionsav, nion)
                            end if
                            fbead(:, :) = rionall(:, :, jj - 1) - rionall(:, :, jj)
                            if (yesperiodize) then
                                call Apply_periodic(fbead, nion, cellscale, dt4, mass_ion, cost)
                                ekinqp = ekinqp + cost
                            else
                                call ApplyPBC(fbead, nion)
                                ekinqp = ekinqp + sum(mass_ion(:, :)*fbead(:, :)**2)
                            end if
                            fbead(:, :) = fbead(:, :) + rionsav(:, :)
                            !                  fbead(:,:) =rionall(:,:,1) - 2*rionall(:,:,jj) + rionall(:,:,jj-1)
                        end if
                        fbead(:, :) = fbead(:, :)*mass_ion(:, :)

                    else ! if iespbc

                        if (jj .gt. 1 .and. jj .lt. nbead) then
                            fbead(:, :) = mass_ion(:, :)*(rionall(:, :, jj + 1) - 2*rionall(:, :, jj) + rionall(:, :, jj - 1))
                            ekinqp = ekinqp + sum(mass_ion(:, :)*(rionall(:, :, jj) - rionall(:, :, jj - 1))**2)
                        elseif (jj .eq. 1) then
                            ekinqp = ekinqp + sum(mass_ion(:, :)*(rionall(:, :, jj) - rionall(:, :, nbead))**2)
                            ! fix ring boundary conditions
                            fbead(:, :) = mass_ion(:, :)*(rionall(:, :, jj + 1) - 2*rionall(:, :, jj) + rionall(:, :, nbead))
                        else
                            fbead(:, :) = mass_ion(:, :)*(rionall(:, :, 1) - 2*rionall(:, :, jj) + rionall(:, :, jj - 1))
                            ekinqp = ekinqp + sum(mass_ion(:, :)*(rionall(:, :, jj) - rionall(:, :, jj - 1))**2)
                        end if

                    end if ! endif iespbc

                    cost = ekinqp
                    call mpi_allreduce(cost, ekinqp, 1                      &
                         &  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                    ekinqp = -ekinqp/nproc*temp**2*0.25d0 + 1.5d0*nion*temp

                    if (.not. yesturboq) call daxpy(ieskin, 0.5d0*temp**2, fbead, 1, force(ndimp + 1), 1)

                end if !yesquantum
#endif

            end if ! fine if ndimp>0

            !      DEALLOCATE A LOT OF USELESS STAFF TO OPTIMIZE MEMORY
#ifdef _OFFLOAD
            deallocate (tabpip, table, tabler, winvdo, winvup&
                 &, wint, dist, tmu, diagfn &
                 &, jastrowall_ee, jastrowall_ei, jasnew_ee, jasnew_ei, winvjbarsz)
            allocate (psip_reweight(iscrapip))
            ! initialize psip
            psip_reweight = 0.d0
#else
            deallocate (tabpip, table, tabler, winv, ainv, winvdo, winvup, winvj    &
                 &, winvbar, winvjbar, wint, dist, tmu, diagfn&
                 &, jastrowall_ee, jastrowall_ei, jasnew_ee, jasnew_ei, winvjbarsz)

            !deallocate (psip)
            allocate (psip_reweight(iscrapip))
            ! initialize psip
            psip_reweight = 0.d0
#endif
#ifdef UNREL
#ifdef PARALLEL
            call mpi_barrier(MPI_COMM_WORLD, ierr) ! for unreliable  networks.
#endif
!$omp barrier
#endif
            time0 = cclock()
            if (rank .eq. 0) write (6, *) ' Time boots =', time0 - timep
            timep = time0
            nweightu = perbin*nweight
            call reweight0(Nw, in1, np, npmn, factorsr                      &
                 &, ipsip, psip_reweight, alphab, sov, econf, econfh, econfion                &
                 &, ieskin, ncg, epsdgel, epstion, 1, wcort, itestrr, enert           &
                 &, scalpar, parcut, iflagerr, force, err, iesconv, itouch                 &
                 &, parcutmin, parcutpar, npbra, kl, etry                                &
                 &, epsi, tpar, beta, ist, ien, rank, ierr, parr, parcute            &
                 &, nbin, fk, dimfk, fkav, okav, skdiag, weightall, reduce, nbinmax, nweightu, ibinit &
                 &, 19, stepcg, lwork, idynu, temp, weight_vir, friction, scalecov, delta0, delta0q, delta0k&
                 &, dt, velion, ris, tmes, cov, rpar, npar, initpar, nparsw, initparsw &
                 &, nparinv, initparinv, endinv, iond, nion, adrlambda, rmax, 20&
                 &, writescratch, countscra, bufscra, jas_invariant, tolcg, ieskinion     &
                 &, rion, iespbc, atom_number, eps_dyn5, maxdev_dyn, acc_dyn, normcorr&
                 &, row_comm, row_id, yescomm)

            time0 = cclock()
            if (rank .eq. 0) write (6, *) ' Total time reweight =', time0 - timep
            timep = time0

            if (change_epscut) then
                wtot(1) = countt
                wtot(2) = countav
                wtot(3) = psirav
                call reduce_base_real_to(3, wtot, costmpi, commrep_mpi, -1)
                if (costmpi(1) .ne. 0.d0) then
                    avreweight = costmpi(2)/costmpi(1)
                else
                    avreweight = 0.8d0 ! do not pass below
                end if
                if (costmpi(2) .ne. 0.d0) then
                    psirav_all = costmpi(3)/costmpi(2)
                else
                    psirav_all = 0.d0
                end if
! renormalizing  the weight of the previous steps assuming that all the steps
! before matters as only 10 steps.
! This is to avoid that the epscut is not changed for long simulations when
! the wf can be instead changed a lot.
                if (countt .gt. tpar_buffer_len*nweight) then
                    countav = countav/countt*nweight*tpar_buffer_len
                    psirav = psirav/countt*nweight*tpar_buffer_len
                    countt = nweight*tpar_buffer_len
                end if
                if ((avreweight .lt. 0.65d0 .or. avreweight .gt. .95d0)&
                     &.and. epscuttype .ne. 0) then
                    if (avreweight .gt. 0.95d0 .and. avreweight .lt. .99d0) then
                        cost = .2d0/(1.d0 - avreweight)
                    elseif (avreweight .ge. .99d0) then
                        cost = 20.d0
                    else
                        cost = avreweight/0.8d0
                    end if
                    if (cost .gt. 3d0) cost = 3.d0
                    if (cost .lt. 0.3333333333d0) cost = 0.3333333333d0
                    epscutu = epscutu*cost
                    epscut = epscutu
                    epstlu = epstlrat*epscutu
                    epstl = epstlrat*epscut
                    if (rank .eq. 0) write (6, *) ' Warning changing epscut in reweight', i_main, epscutu, avreweight
                    countt = 0.d0
                    countav = 0.d0
                    psirav = 0.d0
                end if
            end if

#ifdef PARALLEL

            if (yesquantum .and. nbead .gt. 1 .and. yeswrite .and. acc_dyn) then

                if (idyn .ne. 5) then
                    cost = tmes/nbead**2 ! The true target temperature
                    call mpi_reduce(cost, tmes, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, commcolrep_mpi, ierr)
                end if

                energyq = ener_true(1) + ekinq

                if (rank .eq. 0) write (6, *) 'New quantum energy/kin =', energyq*ris(2), ekinq*ris(2)

                !            Correction to the pressure  + 2 the quantum Kinetic energy.
                if (iespbc) then
                    ! Unit Hartree set in bootpress
                    derEVp = derEVp + 2.d0*costpr*ekinq*ris(2)
                    derEVnopulay = derEVnopulay + 2.d0*costpr*ekinq*ris(2)
                end if

            end if
#endif

            if (rank .eq. 0 .and. yeswrite .and. acc_dyn) then
                if (iespbc) then
                    write (23, 123) derEVp, errEVp, derEVnopulay, errEVnopulay
                    if (flush_write) flush (23)
                end if
                if (yesquantum) then
                    write (16, 123) energyq*ris(2), sigma_true(1)*ris(2), tmes&
                         &, ekinq*ris(2), ekinqp*ris(2)& ! added quantum kin energy,virial/primitive
                         &, (force(cpippo)*ris(3), err(cpippo)*ris(3), cpippo=ndimp + 1, ndimp + ieskin)
                else
                    if (pressfixed .ne. 0.d0) then
                        enthalpy = (ener_true(1) + pressfixed*cellscale(1)*cellscale(2)*cellscale(3))*ris(2)
                        write (16, 123) enthalpy, sigma_true(1)*ris(2), tmes, weight_vir, ener_true(1)*ris(2)&
                             &, (force(cpippo)*ris(3), err(cpippo)*ris(3), cpippo=ndimp + 1, ndimp + ieskin)
                    else
                        write (16, 123) ener_true(1)*ris(2), sigma_true(1)*ris(2), tmes, weight_vir, 0.d0&
                             &, (force(cpippo)*ris(3), err(cpippo)*ris(3), cpippo=ndimp + 1, ndimp + ieskin)
                    end if
                end if
                if (flush_write) flush (16)
                if (write_cov .and. flush_write) flush (18)
                if (idyn .gt. 1) then
                    write (22, 123) (velion(3, jj)/ris(1), jj=1, ieskindim)
                    if (flush_write) flush (22)
                end if
                if (ieskint .ne. 0) then
                    if (yesquantum .and. nbead .gt. 1 .and. yeswritebead) then
                        do ccpippo = 1, nbead
                            if (LBox .gt. 0) then
                                write (17, 123) ((rionall(jj, cpippo, ccpippo)/ris(1), jj=1, 3), cpippo=1, nion)  &
                                     &, cellscale(1:3)
                            else
                                write (17, 123) ((rionall(jj, cpippo, ccpippo)/ris(1), jj=1, 3), cpippo=1, nion)
                            end if
                        end do
                        if (flush_write) flush (17)
                    else
                        if (LBox .gt. 0) then
                            write (17, 123) ((rion(jj, cpippo)/ris(1), jj=1, 3), cpippo=1, nion)  &
                                 &, cellscale(1:3)
                        else
                            write (17, 123) ((rion(jj, cpippo)/ris(1), jj=1, 3), cpippo=1, nion)
                        end if
                        if (flush_write) flush (17)
                    end if

                    !                              if(yesquantum.and.nbead.gt.1) then

                    !             write(27,123) (((rionall(jj,cpippo,ccpippo)/ris(1),jj=1,3),cpippo=1,nion),ccpippo=1,nbead)

                    !                    if(flush_write)  flush(27)

                    !                              endif

                end if
            end if

            !     close and repone the file by scratch
            if (writescratch .eq. 0) then
                rewind (19)
            else
                countscra = 0
            end if

            indc = iesinv + iesm
            if (change_tpar) then
                if (rank .eq. 0) write (6, *) 'tpar before adjust', tpar
                call adjust_tpar(i_main, nweight, ener_true(1)*ris(2), sigma_true(1)*ris(2), tpar, ngentry, itestr4)
            end if

            if (iesd .ne. 0) then
                do k = indc + 1, indc + iesd
                    vj(k - indc) = vj(k - indc) + alphab(k)*tpar

                    if (vj(k - indc) .le. minjonetwobody) then
                        if (rank .eq. 0) write (6, *) ' Warning one/two body too small', vj(k - indc), minjonetwobody
                        vj(k - indc) = minjonetwobody

                    end if
                    alphavar(k) = vj(k - indc)
                end do
            end if

            !       Rescale alphab according to the interval of Z chosen
            indc = iesinv + iesm + iesd + iesfree + iessw
            if (yeszagp)&
                 & call cutminzmaxz(ipc, contraction, iesup, iesupr_c, iesuptransb, jbraiesup&
                 &, minz, maxz, indc, alphavar, alphab, scalpar, fixpar, tpar, nmat, rank)
            indc = iesinv
            if (yeszj)&
                 & call cutminzmaxz(1, contractionj, iesm, npar3bodyr_c, iesuptransbj, jbraiesm&
                 &, minzj, maxzj, indc, alphavar, alphab, scalpar, fixpar, tpar, nmat, rank)

            if (nmolmax .gt. 0) then
                !      changing detmat_c and dup_c according to alphab (correction)
                !      Changing wf before convertmol
                !      One has to define the new Z and the new detmat,jasmat,jasmatsz
                !      used by convertmol
                noproj = .false.
                if (rank .eq. 0) write (6, *) ' Projecting '
                if (iessw .ne. 0) then
                    indc = iesinv + iesm + iesd + iesfree
                    do k = indc + 1, indc + iessw
                        dsw(k - indc) = alphavar(k) + alphab(k)*tpar
                        alphavar(k) = dsw(k - indc)
                        psip(k - indc) = alphab(k)*tpar
                    end do
                    !                 from real to effective
                    if (rank .eq. 0) write (6, *) ' Passi qui XIX real-eff'
                    if (allowed_averagek) call attach_phase2det(.false., detmat_c)
                    call bconstraint(iessw, detmat_c, nelorb_c, nnozero_c&
                         &, nozero_c, psip(iessw + 1), psip, 1, jbradet, symmagp, .true.)

                    !  delayed later
                    !                 from effective to real
                    !                   if(rank.eq.0) write(6,*) ' Passi qui XX eff-real'
                    !                   call attach_phase2det(.true.,detmat_c)
                    !                  endif
                    !                 Here one should check whether is symmetric...

                end if
                if (iesinv .ne. 0) then
                    indc = 0
                    do k = indc + 1, indc + iesfreesz
                        ddwsz(k - indc) = alphavar(k) + alphab(k)*tpar
                        alphavar(k) = ddwsz(k - indc)
                    end do
                    ! update jasmatsz_c
                    if (contractionj .ne. 0) then
                        call bconstrbra(iesfreesz, nnozeroj_c, jbrajsz, nozeroj_c       &
                             &, jasmatsz_c, nelorbj_c, ddwsz)
                    else
                        call bconstrbra(iesfreesz, nnozeroj, jbrajsz, nozeroj           &
                             &, jasmatsz, nelorbjh, ddwsz)
                    end if
                end if

                if (iesfree .ne. 0) then
                    indc = iesinv + iesm + iesd
                    do k = indc + 1, indc + iesfree
                        ddw(k - indc) = alphavar(k) + alphab(k)*tpar
                        alphavar(k) = ddw(k - indc)
                    end do
                    ! update jasmat_c
                    if (contractionj .ne. 0) then
                        call bconstrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c, jasmat_c  &
                             &, ipj*nelorbj_c, ddw)
                    else
                        if (yes_sparse) then
                            call bconstrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder, jasmat&
                          &, ipj*nelorbjh, ddw)
                        else
                            call bconstrbra(iesfree, nnozeroj, jbraj, nozeroj, jasmat&
                                 &, ipj*nelorbjh, ddw)
                        end if
                    end if
                end if
                if (iesup .ne. 0) then
                    indc = iesinv + iesm + iesd + iesfree + iessw
                    do k = indc + 1, indc + iesup
                        dup(k - indc) = alphavar(k) + tpar*alphab(k)
                        alphavar(k) = dup(k - indc)
                    end do
                    call bconstrbr_complex(iesup, iesupr_c, jbraiesup, dup_c, dup)
                end if

                if (contraction .eq. 0) then
                    if (yes_complex) then
                        do jj = 1, iesupr_c
                            dupr(jj) = dup_c(2*jj - 1)
                        end do
                    else
                        do jj = 1, iesupr_c
                            dupr(jj) = dup_c(jj)
                        end do
                    end if
                else
                    ! restore the possible uncontracted coefficients
                    !                        if(allocated(muc_np)) mu_c=muc_np
                    if (yes_complex) then
                        do jj = 1, iesup_c
                            do kk = 1, multranspip(jj)
                                iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                ix = transpip(kk)%col(jj) - (iy - 1)*ipf*nelorbh
                                mu_c(2*ix - 1, iy) = dup_c(2*jj - 1)
                                mu_c(2*ix, iy) = dup_c(2*jj)
                            end do
                        end do
                        do jj = 1, iesupr_c
                            if (iesuptransb(jj) .ne. 0) then
                                dupr(iesuptransb(jj)) = dup_c(2*jj - 1)
                            end if
                        end do
                    else
                        do jj = 1, iesup_c
                            do kk = 1, multranspip(jj)
                                iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                ix = transpip(kk)%col(jj) - (iy - 1)*(ipf*nelorbh)
                                mu_c(ix, iy) = dup_c(jj)
                            end do
                        end do
                        do jj = 1, iesupr_c
                            if (iesuptransb(jj) .ne. 0) then
                                dupr(iesuptransb(jj)) = dup_c(jj)
                            end if
                        end do
                    end if
                    !                       if(allocated(muc_np)) muc_np=mu_c
                end if ! endif contraction
                if (iesm .ne. 0) then
                    indc = iesinv
                    do k = indc + 1, indc + iesm
                        vju(k - indc) = vju(k - indc) + alphab(k)*tpar
                        alphab(k) = 0.d0
                        alphavar(k) = vju(k - indc)
                    end do
                    call bconstrbr(iesm, npar3bodyr_c, jbraiesm, vju_c, vju)
                end if

                if (contractionj .eq. 0) then

                    do jj = 1, npar3bodyr_c
                        vjur(jj) = vju_c(jj)
                    end do

                else

                    ! change muj and vjur
                    do jj = 1, npar3body_c
                        do kk = 1, multranspipj(jj)
                            iy = (transpipj(kk)%col(jj) - 1)/nelorbjh + 1
                            ix = transpipj(kk)%col(jj) - (iy - 1)*nelorbjh
                            muj_c(ix, iy) = vju_c(jj)
                        end do
                    end do

                    do jj = 1, npar3bodyr_c
                        if (iesuptransbj(jj) .ne. 0) then
                            vjur(iesuptransbj(jj)) = vju_c(jj)
                        end if
                    end do

                end if
                !                     Store the old ion positions before writing it, namely
                !                     consistent with the last step optimization.

                call update_ionpos
                if (yeszagp .or. cellderiv) call update_kgrid
                if (nmolmax .gt. 0 .and. iessw .gt. 0) then
                    if (rank .eq. 0) write (6, *) ' Passi qui XX eff-real'
                    if (allowed_averagek) call attach_phase2det(.true., detmat_c)
                    !                 Here one should check whether is symmetric...
                end if

                if (contractionj .ne. 0) then
                    !                            call scontract_mat_jas(nelorbjh,nelorbjh,nelorbjh&
                    !                                &,nelorbj_c,nelorbj_c,jasmat,jasmat_c,muj_c,psip) Ichanged

!            if(ipj.eq.2) then
!               call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!            else
!               call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh&
!                    &, nelorbj_c, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!            endif
                    if (iessz) call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
                         &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)
                end if

                if (contraction .gt. 0) then

                    timepp = cclock()
                    !  if(yesfast.ne.0.and.(2*nelorb_c-molecular.lt.nelorbh&
                    ! &.or.projnmin)) then
                    !      if(yesfast.ne.0) then
                    !                         if(dofast) then
                    !                           if(mod(i/nweight,max(ifreqdmrg,1)).ne.0.or.ifreqdmrg.le.0) then

                    if (detc_proj) then
                        call convertmol_c
                    else
                        call convertmol_fast
                    end if
                    if (yesdetmatc) then
                        call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                             &, nelcol_c, detmat, detmat_c, mu_c, psip)
                    end if
                    if (rank .eq. 0) write (6, *) ' Time convertmol =', cclock() - timepp
                    if (iessw .gt. 0) then
                        !       from real to effective load dsw effective parameters if allowed_averagek
                        if (rank .eq. 0) write (6, *) ' Passi qui XXI real-eff'
                        if (allowed_averagek) call attach_phase2det(.false., detmat_c)
                        call constrbra_complex(iessw, nnozero_c, jbradet, nozero_c, detmat_c      &
                             &, dsw, 1, 1)
                        !       from effective to real
                        if (allowed_averagek) call attach_phase2det(.true., detmat_c)
                        if (rank .eq. 0) write (6, *) ' Passi qui XXII eff-real'
                        indc = iesinv + iesm + iesd + iesfree
                        do k = indc + 1, indc + iessw
                            alphavar(k) = dsw(k - indc)
                        end do
                    end if ! endif contraction>0

                    !       Symmetrization
                    if (symiesup) then
                        call constrbr_complex(iesupind, iesupr_c, jbrasymiesup, dup_c, psip, 1, 1)
                        call bconstrbr_complex(iesupind, iesupr_c, jbrasymiesup, dup_c, psip)
                        if (yes_complex) then
                            do jj = 1, iesupr_c
                                if (iesuptransb(jj) .ne. 0) then
                                    dupr(iesuptransb(jj)) = dup_c(2*jj - 1)
                                end if
                            end do
                        else
                            do jj = 1, iesupr_c
                                if (iesuptransb(jj) .ne. 0) then
                                    dupr(iesuptransb(jj)) = dup_c(jj)
                                end if
                            end do
                        end if
                    end if
                    if (contraction .ne. 0) then
                        !                           if(allocated(muc_np)) mu_c=muc_np
                        ! restore the possible uncontracted coefficients
                        if (yes_complex) then
                            do jj = 1, iesup_c
                                do kk = 1, multranspip(jj)
                                    iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                    ix = transpip(kk)%col(jj) - (iy - 1)*(ipf*nelorbh)
                                    mu_c(2*ix - 1, iy) = dup_c(2*jj - 1)
                                    mu_c(2*ix, iy) = dup_c(2*jj)
                                end do
                            end do
                        else
                            do jj = 1, iesup_c
                                do kk = 1, multranspip(jj)
                                    iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                    ix = transpip(kk)%col(jj) - (iy - 1)*(ipf*nelorbh)
                                    mu_c(ix, iy) = dup_c(jj)
                                end do
                            end do
                        end if
                        !                           if(allocated(muc_np)) muc_np=mu_c

                        call update_projm

                    end if

                    !       Put in any event the consistency between alphavar and dup dup_c

                    if (iesup .ne. 0) then
                        indc = iesinv + iesm + iesd + iesfree + iessw
                        call constrbr_complex(iesup, iesup_c, jbraiesup, dup_c, dup, 1, 1)
                        do k = indc + 1, indc + iesup
                            alphavar(k) = dup(k - indc)
                        end do
                    end if

                end if

                if (contractionj .ne. 0) then
                    ! change muj and vjur
                    do jj = 1, npar3body_c
                        do kk = 1, multranspipj(jj)
                            iy = (transpipj(kk)%col(jj) - 1)/nelorbjh + 1
                            ix = transpipj(kk)%col(jj) - (iy - 1)*nelorbjh
                            muj_c(ix, iy) = vju_c(jj)
                        end do
                    end do
                    if (iesm .ne. 0) then
                        call constrbr(iesm, npar3body_c, jbraiesm, vju_c, vju, 1, 1)
                        indc = iesinv
                        do k = indc + 1, indc + iesm
                            alphavar(k) = vju(k - indc)
                        end do
                    end if
                end if

                !

                !
                !      alphab is the correction put to zero and alphavar
                !      hould be restored
                !      if parallel send the info dup_c detmat_c to all proc
            else
                noproj = .true.
            end if
            !               endif
            time0 = cclock()
            if (rank .eq. 0) write (6, *) ' Time around convertmol =', time0 - timep

            timep = time0

#ifdef  _OFFLOAD
            allocate (tabpip((indt + ip4)*nel*nws), table(max(ipc*nel*indt, ipc)*nws)&
                 &, tabler(max(nel*indt, 1)*nws)&
                 &, winvdo(ipc*(indt + ip4)*max(neldo, 1)*nws)             &
                 &, winvup(ipc*(indt + ip4)*nelup*nws) &
                 &, wint(nws), dist(nel*nws*nion), tmu(nel*max(indt, 1)*nws)&
                 &, diagfn(nws) &
                                     !    &,derdet_mu(nelorb*nelorb_c),derdet_c(nelorb_c*max(nelcol_c,nmaxder))&
                                     !    &,derjas_mu(nelorbj*nelorbj_c),derjas_c(nelorbj_c*nelorbj_c)       &
                 &, jastrowall_ee(nel, nel, 0:indt4j, nws), jasnew_ee(nel)                    &
                 &, jastrowall_ei(nion, nel, nws), jasnew_ei(nion))
            deallocate (psip_reweight)
#else
            !      REALLOCATE A LOT OF USEFUL STAFF
            deallocate (psip_reweight)
            !allocate (psip(iscramax))
            allocate (tabpip((indt + ip4)*nel*nws), table(max(ipc*nel*indt, ipc)*nws)&
                 &, tabler(max(nel*indt, 1)*nws), winv(ipc*nelorb*nel*(indt4 + 1)*nws)&
                 &, ainv(ipc*nelup_mat*nelup_mat*nws), winvdo(ipc*(indt + ip4)*max(neldo, 1)*nws)             &
                 &, winvup(ipc*(indt + ip4)*nelup*nws), winvj(max(nelorbj, 1)*nel*(indt4j + 1)*nws) &
                 &, winvbar(ipc*ipf*nelorbh*nel_mat*nws), winvjbar(ipj*nelorbjh*nel*nws + nws)           & !Ichanged
                 &, wint(nws), dist(nel*nws*nion), tmu(nel*max(indt, 1)*nws)&
                 &, diagfn(nws) &
                                     !    &,derdet_mu(nelorb*nelorb_c),derdet_c(nelorb_c*max(nelcol_c,nmaxder))&
                                     !    &,derjas_mu(nelorbj*nelorbj_c),derjas_c(nelorbj_c*nelorbj_c)       &
                 &, jastrowall_ee(nel, nel, 0:indt4j, nws), jasnew_ee(nel)                    &
                 &, jastrowall_ei(nion, nel, nws), jasnew_ei(nion))
#endif

            if (iessz) then
                allocate (winvjbarsz(max(nelorbjh, 1)*nel*nws))
            else
                allocate (winvjbarsz(1))
            end if

            ! reinitialize all  from zero
            psip = 0.d0
            winvbar = 0.d0
            winvjbar = 0.d0
            tabpip = 0.d0
            table = 0.d0
            tabler = 0.d0
            winv = 0.d0
            ainv = 0.d0
            winvdo = 0.d0
            winvup = 0.d0
            winvj = 0.d0
            wint = 0.d0
            !               wintw=0.d0
            dist = 0.d0
            tmu = 0.d0
            diagfn = 0.d0
            jastrowall_ee = 0.d0
            jastrowall_ei = 0.d0
            jasnew_ee = 0.d0
            jasnew_ei = 0.d0
            winvjbarsz = 0.d0

            !======================================================
            ! start automatic adjustments of tpar and parr
            ! added by K.Nakano on 23. Sep. 2019
            !======================================================
            !write(6,*) 'i_main', i_main
            !write(6,*) 'ngen', ngen
            !write(6,*) 'equil_steps', equil_steps
            !write(6,*) 'ngen-equil_steps', ngen-equil_steps
            !if (change_tpar .and. i_main .gt. (ngen-equil_steps)) then
            !    if(rank.eq.0) write(6,*) ' Warning: tpar, entering an equilibrium step.'
            !    if(rank.eq.0) write(6,*) ' Warning: tpar, no longer changes automatically.'
            !    change_tpar=.false.
            !end if

            !======================================================
            ! end automatic adjustments of tpar and parr
            !======================================================

            if (rank .eq. 0) then
                if (inext + nweight - iend .gt. ngen) then
                    write (6, *) 'Warning  stopping the program to terminate VMC bin'
                    ngen = inext - iend
                end if
                open (unit=7, file='stop.dat', form='formatted', status='unknown')
                if (cclock() - inittime .gt. maxtime .or. ngentry .eq. -5) then
                    ngentry = 0
                else
                    read (7, *, end=1155) ngentry
                end if
                if (ngentry .eq. 0) then
                    ngen = inext - iend
                    write (6, *) ' The program will stop at iteration', ngen
                elseif (ngentry .gt. 0) then
                    ngen = ngentry
                    write (6, *) ' The program will stop at iteration', ngen + iend
                elseif (ngentry .eq. -1) then
                    !        read also parr and epsi
                    if (ncg .eq. 0) then
                        read (7, *, end=1155) epsdgel, epsi
                        parr = epsdgel
                    else
                        read (7, *, end=1155) parr, epsi
                    end if
                    if (abs(parr)/10.d0 .lt. tolcg) tolcg = abs(parr)/10.d0
!            if(default_epsdgel.and.ncg.ne.0) epsdgel = abs(parr) / 10.d0

                    write (6, *) ' Warning changing parr and epsi on fly  ', parr, epsi
                elseif (ngentry .eq. -2) then
                    if (ncg .eq. 0) then
                        read (7, *, end=1155) epsdgel, epsi, tpar
                        parr = epsdgel
                    else
                        read (7, *, end=1155) parr, epsi, tpar
                    end if
                    write (6, *) ' Warning changing parr,epsi,tpar on fly  ', parr, epsi, tpar
                    if (abs(parr)/10.d0 .lt. tolcg) tolcg = abs(parr)/10.d0
!            if(default_epsdgel.and.ncg.ne.0) epsdgel = abs(parr) / 10.d0
                elseif (ngentry .eq. -3) then
                    if (ncg .eq. 0) then
                        read (7, *, end=1155) epsdgel, epsi, tpar, rweight
                        parr = epsdgel
                    else
                        read (7, *, end=1155) parr, epsi, tpar, rweight
                    end if
                    if (abs(parr)/10.d0 .lt. tolcg) tolcg = abs(parr)/10.d0
!            if(default_epsdgel.and.ncg.ne.0) epsdgel = abs(parr) / 10.d0
                    write (6, *) ' Warning changing parr,epsi,tpar,nweight on fly '&
                         &, parr, epsi, tpar, rweight
                    if (iskipdyn .gt. 1 .and. rweight .ne. nweight) then
                        write (6, *) ' ERROR nweight cannot be changed during dynamic with iskipdyn>1'
                        iflagerr = 1
                    end if
                elseif (ngentry .eq. -4) then
                    read (7, *, end=1155) cost, epstion
                    cost = cost/(2.d0*dt)
                    dt = dt*cost
                    call dscal(ieskin, cost, scalpar(np - ieskin + 1), 1)
                    write (6, *) ' Warning changing tion on fly New tion (H)=', dt*2.d0
                end if
1155            continue
                if (ngentry .ne. -3 .and. ngentry .ne. -4) ngentry = 0
                close (7)
            end if
#ifdef  PARALLEL
            call mpi_bcast(ngen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(parr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(tolcg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(epsdgel, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(epsi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(tpar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(tion, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(rweight, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ngentry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            if (ngentry .eq. -4) then
                call mpi_bcast(scalpar(np - ieskin + 1), ieskin&
                     &, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(epstion, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            end if
#endif
            if (ngentry .eq. -4) ngentry = 0

            if (itest .eq. 2 .and. itestr .ne. -5) then

#ifdef PARALLEL
                comm_mpi = commrep_mpi ! always the smallest communicator.

                call mpi_allreduce(psirav, psirav_all, 1                                &
                     &  , MPI_DOUBLE_PRECISION, MPI_SUM, COMM_mpi, ierr)
                call mpi_allreduce(countav, countav_all, 1&
                     &  , MPI_DOUBLE_PRECISION, MPI_SUM, COMM_mpi, ierr)
                call mpi_allreduce(countt, countt_all, 1&
                     &  , MPI_DOUBLE_PRECISION, MPI_SUM, COMM_mpi, ierr)
#else
                countav_all = countav
                psirav_all = psirav
                countt_all = countt
#endif
                if (countav_all .ne. 0.d0) then
                    psirav_all = psirav_all/countav_all
                else
                    psirav_all = 0.d0
                end if

            end if

            ! by. E. Coccia (7/12/10): write the final additional energy if optimization is running
            if (ext_pot) then
#ifdef PARALLEL
                ! electrons
                call mpi_reduce(ave, t_ave, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2, t_ave2, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount, t_ncount, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(nout, t_nout, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                ! nuclei
                call mpi_reduce(ave_ion, t_ave_ion, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_ion, t_ave2_ion, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_ion, t_ncount_ion, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                ! electrons
                t_ave = ave
                t_ave2 = ave2
                t_ncount = ncount
                t_nout = nout
                ! nuclei
                t_ave_ion = ave_ion
                t_ave2_ion = ave2_ion
                t_ncount_ion = ncount_ion
#endif
                ! by E. Coccia (4/2/11):  write the vdw energy if optimization is running
                if (vdw) then
#ifdef PARALLEL
                    call mpi_reduce(ave_vdw, t_ave_vdw, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                    call mpi_reduce(ave2_vdw, t_ave2_vdw, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                    call mpi_reduce(ncount_vdw, t_ncount_vdw, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                    t_ave_vdw = ave_vdw
                    t_ave2_vdw = ave2_vdw
                    t_ncount_vdw = ncount_vdw
#endif
                    ! by E. Coccia (27/5/11): write the classical potential energies
                    if (link_atom) then
                        ! Bond angle
#ifdef  PARALLEL
                        call mpi_reduce(ave_angle, t_ave_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                        call mpi_reduce(ave2_angle, t_ave2_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                        call mpi_reduce(ncount_angle, t_ncount_angle, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                        t_ave_angle = ave_angle
                        t_ave2_angle = ave2_angle
                        t_ncount_angle = ncount_angle
#endif
                        ! Proper dihedral
#ifdef PARALLEL
                        call mpi_reduce(ave_dihed, t_ave_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                        call mpi_reduce(ave2_dihed, t_ave2_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                        call mpi_reduce(ncount_dihed, t_ncount_dihed, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                        t_ave_dihed = ave_dihed
                        t_ave2_dihed = ave2_dihed
                        t_ncount_dihed = ncount_dihed
#endif
                        ! Improper dihedral
#ifdef PARALLEL
                        call mpi_reduce(ave_impr, t_ave_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                        call mpi_reduce(ave2_impr, t_ave2_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                        call mpi_reduce(ncount_impr, t_ncount_impr, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                        t_ave_impr = ave_impr
                        t_ave2_impr = ave2_impr
                        t_ncount_impr = ncount_impr
#endif
                    end if
                end if
            end if
            ! by E. Coccia (10/12/11): MM restraints
            if (mm_restr) then
                ! Bond distance
#ifdef  PARALLEL
                call mpi_reduce(ave_bond, t_ave_bond, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_bond, t_ave2_bond, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_bond, t_ncount_bond, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                t_ave_bond = ave_bond
                t_ave2_bond = ave2_bond
                t_ncount_bond = ncount_bond
#endif

                ! Bond angle
#ifdef  PARALLEL
                call mpi_reduce(ave_angle, t_ave_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_angle, t_ave2_angle, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_angle, t_ncount_angle, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                t_ave_angle = ave_angle
                t_ave2_angle = ave2_angle
                t_ncount_angle = ncount_angle
#endif
                ! Proper dihedral
#ifdef PARALLEL
                call mpi_reduce(ave_dihed, t_ave_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_dihed, t_ave2_dihed, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_dihed, t_ncount_dihed, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                t_ave_dihed = ave_dihed
                t_ave2_dihed = ave2_dihed
                t_ncount_dihed = ncount_dihed
#endif
                ! Improper dihedral
#ifdef PARALLEL
                call mpi_reduce(ave_impr, t_ave_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ave2_impr, t_ave2_impr, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(ncount_impr, t_ncount_impr, 1, MPI_INTEGER, mpi_sum, 0, MPI_COMM_WORLD, ierr)
#else
                t_ave_impr = ave_impr
                t_ave2_impr = ave2_impr
                t_ncount_impr = ncount_impr
#endif
            end if

            if (rank .eq. 0) then

                !by E. Coccia (9/4/14)

                !

                ! print the energy of the last iteration
                if (pressfixed .ne. 0.d0) then
                    enthalpy = (ener_true(1) + pressfixed*cellscale(1)*cellscale(2)*cellscale(3))*ris(2)
                    write (6, *) ' New Enthalpy/Energy = ', enthalpy, ener_true(1)*ris(2), sigma_true(1)*ris(2)
                else
                    write (6, *) ' New Energy = ', ener_true(1)*ris(2), sigma_true(1)*ris(2)
                end if
                if (itest .eq. 2 .and. psirav_all .ne. 0.d0) write (6, *) ' Average inverse A wf =', psirav_all
                ! by E. Coccia (7/12/10): write the external potential
                if (ext_pot) then
                    call extpot_final(nel)
                    ! by E. Coccia (13/1/11): nuclear potential
                    call ion_final(nion)
                    ! by E. Coccia (4/2/11): vdw energy
                    if (vdw) then
                        call vdw_final()
                    end if
                end if
                ! by E. Coccia (10/12/11): MM restraints
                if (mm_restr) then
                    call vdw_final()
                    write (6, *) ' New Energy (no MM) = ', ener_true(1)*ris(2) - sum_pot, sigma_true(1)*ris(2)
                    write (*, *) '|******************************************|'
                    write (*, *) '|         EXTERNAL QMC/MM POTENTIAL        |'
                    write (*, *) '|******************************************|'
                    write (*, *) ''
                end if

                !         check the energy is lower and the variance is not too large
                if (ngentry .eq. -3) then
                    if (ngentry .eq. -3) then
                        write (6, *) ' Initializing again '
                        ndone = rweight - ibinit
                        lbin = ndone/nbinr
                        ndone = nbinr*lbin
                        rweight = ndone + ibinit
                        write (6, *) ' Changing nweight = ', rweight
                        iesconv = 0
                    end if
                end if
            end if ! rank.eq.0

            ! by E. Coccia (2/3/11)
            if (ext_pot) then
                ! initialize variables for the next optimization step
                ! electrons
                ave = 0.d0; ave2 = 0.d0; ncount = 0; nout = 0
                t_ave = 0.d0; t_ave2 = 0.d0; t_ncount = 0; t_nout = 0
                !nuclei
                ave_ion = 0.d0; ave2_ion = 0.d0; ncount_ion = 0
                t_ave_ion = 0.d0; t_ave2_ion = 0.d0; t_ncount_ion = 0
                if (vdw) then
                    ave_vdw = 0.d0; ave2_vdw = 0.d0; ncount_vdw = 0
                    t_ave_vdw = 0.d0; t_ave2_vdw = 0.d0; t_ncount_vdw = 0
                    if (link_atom) then
                        ave_angle = 0.d0; ave2_angle = 0.d0; ncount_angle = 0
                        t_ave_angle = 0.d0; t_ave2_angle = 0.d0; t_ncount_angle = 0
                        ave_dihed = 0.d0; ave2_dihed = 0.d0; ncount_dihed = 0
                        t_ave_dihed = 0.d0; t_ave2_dihed = 0.d0; t_ncount_dihed = 0
                        ave_impr = 0.d0; ave2_impr = 0.d0; ncount_impr = 0
                        t_ave_impr = 0.d0; t_ave2_impr = 0.d0; t_ncount_impr = 0
                    end if
                end if
            end if
            ! by E. Coccia (10/12/11): MM restraints
            if (mm_restr) then
                ave_bond = 0.d0; ave2_bond = 0.d0; ncount_bond = 0
                t_ave_bond = 0.d0; t_ave2_bond = 0.d0; t_ncount_bond = 0
                ave_angle = 0.d0; ave2_angle = 0.d0; ncount_angle = 0
                t_ave_angle = 0.d0; t_ave2_angle = 0.d0; t_ncount_angle = 0
                ave_dihed = 0.d0; ave2_dihed = 0.d0; ncount_dihed = 0
                t_ave_dihed = 0.d0; t_ave2_dihed = 0.d0; t_ncount_dihed = 0
                ave_impr = 0.d0; ave2_impr = 0.d0; ncount_impr = 0
                t_ave_impr = 0.d0; t_ave2_impr = 0.d0; t_ncount_impr = 0
            end if

            if (rank .eq. 0) write (6, *) ' Used epscut,epstl =', epscutu, epstlu

#ifdef PARALLEL
            call mpi_bcast(tpar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(rweight, 1, MPI_DOUBLE_PRECISION                     &
                 &, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ndone, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(lbin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ngen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
            !   For unreliable  networks.
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
#endif

            time0 = cclock()

            if (rank .eq. 0) write (6, *) ' Time before changing guiding ', time0 - timep

            timep = time0

            !by Andrea Tirelli: ADAM optimizer
            if (yes_adams) then
                if (rank .eq. 0) write (6, *) ' Perfoming ADAM optimization'
                call adam_opt(i_main, nweight, ndimp, first_moment, second_moment, alphab)
            end if

            !        CHANGING THE GUIDING FUNCTION

            if (np .gt. 0) then

                indc = 0
                if (iesinv .ne. 0 .and. noproj) then
                    do k = 1, iesfreesz
                        alphavar(k) = tpar*alphab(k) + alphavar(k)
                    end do
                    do k = indc + 1, indc + iesfreesz
                        ddwsz(k - indc) = alphavar(k)
                    end do
                    ! update jasmatsz_c
                    if (contractionj .ne. 0) then
                        call bconstrbra(iesfreesz, nnozeroj_c, jbrajsz, nozeroj_c       &
                             &, jasmatsz_c, nelorbj_c, ddwsz)
                    else
                        call bconstrbra(iesfreesz, nnozeroj, jbrajsz, nozeroj           &
                             &, jasmatsz, nelorbjh, ddwsz)
                    end if
                end if
                indc = indc + iesinv
                if (iesm .ne. 0 .and. noproj) then
                    do k = indc + 1, indc + iesm
                        vju(k - indc) = vju(k - indc) + alphab(k)*tpar
                        alphavar(k) = vju(k - indc)
                    end do
                    call bconstrbr(iesm, npar3bodyr_c, jbraiesm, vju_c, vju)
                    if (contractionj .eq. 0) then
                        do jj = 1, npar3bodyr_c
                            vjur(jj) = vju_c(jj)
                        end do
                    else
                        ! change muj and vjur
                        do jj = 1, npar3body_c
                            do kk = 1, multranspipj(jj)
                                iy = (transpipj(kk)%col(jj) - 1)/nelorbjh + 1
                                ix = transpipj(kk)%col(jj) - (iy - 1)*nelorbjh
                                muj_c(ix, iy) = vju_c(jj)
                            end do
                        end do
                        do jj = 1, npar3bodyr_c
                            if (iesuptransbj(jj) .ne. 0) then
                                vjur(iesuptransbj(jj)) = vju_c(jj)
                            end if
                        end do
                    end if
                end if ! iesm.ne.0
                indc = indc + iesm
                if (iesd .ne. 0 .and. noproj) then
                    do k = indc + 1, indc + iesd
                        vj(k - indc) = vj(k - indc) + alphab(k)*tpar
                        if (vj(k - indc) .le. minjonetwobody) then
                            if (rank .eq. 0) write (6, *) ' Warning one/two body too small', vj(k - indc), minjonetwobody
                            vj(k - indc) = minjonetwobody

                        end if
                        alphavar(k) = vj(k - indc)
                    end do
                end if
                indc = indc + iesd
                if (iesfree .ne. 0 .and. noproj) then
                    do k = indc + 1, indc + iesfree
                        ddw(k - indc) = alphavar(k) + tpar*alphab(k)
                        alphavar(k) = ddw(k - indc)
                    end do
                    ! update jasmat_c
                    if (contractionj .ne. 0) then
                        call bconstrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c, jasmat_c  &
                             &, ipj*nelorbj_c, ddw)
                    else
                        if (yes_sparse) then
                            call bconstrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder, jasmat&
                          &, ipj*nelorbjh, ddw)
                        else
                            call bconstrbra(iesfree, nnozeroj, jbraj, nozeroj, jasmat&
                          &, ipj*nelorbjh, ddw)
                        end if
                    end if
                end if
                indc = indc + iesfree

                if ((iesinv .ne. 0 .or. iesfree .ne. 0 .or. iesm .ne. 0)&
                     &.and. contractionj .ne. 0 .and. noproj) then
                    ! update jasmat, jasmatsz

                    !                      call scontract_mat_jas(nelorbjh,nelorbjh,nelorbjh,nelorbj_c&
                    !                           &,nelorbj_c,jasmat,jasmat_c,muj_c,psip)                            Ichanged

!            if(ipj.eq.2) then
!               call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!            else
!               call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh&
!                    &, nelorbj_c, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!            endif
                    if (iessz) then
                        call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
                             &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)
                    end if
                end if

                if (iessw .ne. 0 .and. noproj) then
                    do k = indc + 1, indc + iessw
                        dsw(k - indc) = alphavar(k) + tpar*alphab(k)
                        psip(k - indc) = tpar*alphab(k) ! correction in this case
                        alphavar(k) = dsw(k - indc)
                    end do

                    if (contraction .eq. 0) then
                        !    From real to effective
                        if (rank .eq. 0) write (6, *) ' Passi qui XXIII real-eff'
                        if (allowed_averagek) call attach_phase2det(.false., detmat)

                        call bconstraint(iessw, detmat, ipf*nelorbh, nnozero&
                             &, nozero, psip(iessw + 1), psip, 1, jbradet, symmagp, .true.)
                        if (symmetrize_agp) then
                            call symmetrizeagp(nnozero_c, nozero_c, jbradet&
                                 &, jbradetn, dsw&
                                 &, iessw0, psip, ipsip, detmat, nelorb_c, nelorb_at, nelcol_c&
                                 &, symmagp, yes_hermite)
                            !                  Restore dsw used.
                            dsw(1:iessw) = alphavar(indc + 1:indc + iessw)
                        end if
                        !    back to real

                        !    if(rank.eq.0) write(6,*) ' Passi qui XXIV eff-real'
                        !               if(allowed_averagek) call attach_phase2det(.true.,detmat)

                    else
                        !    From real to effective
                        !         if(rank.eq.0) then
                        !      write(6,*) ' Before symmetrize correction=',rank,sum(abs(psip(1:iessw)))
                        !           write(6,*) ' Fixed phase elements '
                        !           do ix=1,nelorb_at
                        !            do iy=ix,nelorb_at
                        !            if(kiontot(ix).ne.kiontot(iy)) write(6,*)&
                        !                & ix,iy,kiontot(ix),kiontot(iy)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix-1)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix)
                        !            enddo
                        !           enddo
                        !           write(6,*) ' Unfixed phase elements '
                        !           do ix=1,nelorb_at
                        !            do iy=ix,nelorb_at
                        !            if(kiontot(ix).eq.kiontot(iy)) write(6,*)&
                        !                & ix,iy,kiontot(ix),kiontot(iy)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix-1)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix)
                        !            enddo
                        !           enddo
                        !         endif
                        if (rank .eq. 0) write (6, *) ' Passi qui XXV real-eff'
                        !      from real to effective
                        if (allowed_averagek) call attach_phase2det(.false., detmat_c)
                        call bconstraint(iessw, detmat_c, nelorb_c, nnozero_c&
                             &, nozero_c, psip(iessw + 1), psip, 1, jbradet, symmagp, .true.)
                        !               if(allowed_averagek) call attach_phase2det(.true.,detmat_c)
                        !                cost=sum(abs(detmat_c(1:ipc*nelorb_c*nelorb_c)))
                        !               if(allowed_averagek) call attach_phase2det(.false.,detmat_c)

                        if (symmetrize_agp) then

                            call symmetrizeagp(nnozero_c, nozero_c, jbradet&
                                 &, jbradetn, dsw&
                                 &, iessw0, psip, ipsip, detmat_c, nelorb_c, nelorb_at, nelcol_c&
                                 &, symmagp, yes_hermite)
                            !                  Restore dsw used.
                            dsw(1:iessw) = alphavar(indc + 1:indc + iessw)
                        end if

                        !                if(rank.eq.0) write(6,*) ' Before attach after symm '
                        !    back to real
                        !    if(rank.eq.0) write(6,*) ' Passi qui XXVI eff-real'
                        !               if(allowed_averagek) call attach_phase2det(.true.,detmat_c)
                        !        if(rank.eq.0) then
                        !              write(6,*) ' After symmetrize ',rank,abs(cost-sum(abs(detmat_c(1:ipc*nelorb_c*nelorb_c))))
                        !           write(6,*) ' Fixed phase elements '
                        !           do ix=1,nelorb_at
                        !            do iy=ix,nelorb_at
                        !            if(kiontot(ix).ne.kiontot(iy)) write(6,*)&
                        !                & ix,iy,kiontot(ix),kiontot(iy)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix-1)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix)&
                        !                &,detmat_c(ipc*nelorb_c*(ix-1)+2*iy-1)&
                        !                &,detmat_c(ipc*nelorb_c*(ix-1)+2*iy)
                        !            enddo
                        !           enddo
                        !           write(6,*) ' Unfixed phase elements '
                        !           do ix=1,nelorb_at
                        !            do iy=ix,nelorb_at
                        !            if(kiontot(ix).eq.kiontot(iy)) write(6,*)&
                        !                & ix,iy,kiontot(ix),kiontot(iy)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix-1)&
                        !                &,detmat_c(ipc*nelorb_c*(iy-1)+2*ix)
                        !            enddo
                        !           enddo
                        !         endif
                    end if
                    !                    call mpi_finalize(ierr)
                    !                    stop
                    ! if iessw.ne.0
                end if
                indc = indc + iessw

                if (iesup .ne. 0 .and. noproj) then

                    !                     if(rank.eq.0) write(6,*) ' dup done '
                    do k = indc + 1, indc + iesup
                        dup(k - indc) = alphavar(k) + tpar*alphab(k)
                        alphavar(k) = dup(k - indc)
                        !                     if(rank.eq.0) write(6,*) k-indc,dup(k-indc)
                    end do

                    call bconstrbr_complex(iesup, iesupr_c, jbraiesup, dup_c, dup)

                    if (contraction .eq. 0) then
                        if (yes_complex) then
                            do jj = 1, iesupr_c
                                dupr(jj) = dup_c(2*jj - 1)
                            end do
                        else
                            do jj = 1, iesupr_c
                                dupr(jj) = dup_c(jj)
                            end do
                        end if
                    else
                        !                        if(allocated(muc_np)) mu_c=muc_np
                        ! restore the possible uncontracted coefficients

                        if (yes_complex) then
                            ! change mu and dupr
                            do jj = 1, iesup_c
                                do kk = 1, multranspip(jj)
                                    iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                    ix = transpip(kk)%col(jj) - (iy - 1)*(ipf*nelorbh)
                                    mu_c(2*ix - 1, iy) = dup_c(2*jj - 1)
                                    mu_c(2*ix, iy) = dup_c(2*jj)
                                end do
                            end do
                            do jj = 1, iesupr_c
                                if (iesuptransb(jj) .ne. 0) then
                                    dupr(iesuptransb(jj)) = dup_c(2*jj - 1)
                                end if
                            end do
                        else
                            do jj = 1, iesup_c
                                do kk = 1, multranspip(jj)
                                    iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                    ix = transpip(kk)%col(jj) - (iy - 1)*(ipf*nelorbh)
                                    mu_c(ix, iy) = dup_c(jj)
                                end do
                            end do
                            do jj = 1, iesupr_c
                                if (iesuptransb(jj) .ne. 0) then
                                    dupr(iesuptransb(jj)) = dup_c(jj)
                                end if
                            end do
                        end if
                        !                       if(allocated(muc_np)) muc_np=mu_c
                    end if

                end if ! iesup.ne.0

                !     indc=indc+iesup
                if (noproj) call update_ionpos

                if (noproj .and. iessw .ne. 0) then
                    if (contraction .ne. 0) then
                        if (rank .eq. 0) write (6, *) ' Passi qui XXVI eff-real'
                        if (allowed_averagek) call attach_phase2det(.true., detmat_c)
                    else
                        if (rank .eq. 0) write (6, *) ' Passi qui XXIV eff-real'
                        if (allowed_averagek) call attach_phase2det(.true., detmat)
                    end if
                end if

                if (noproj .and. (yeszagp .or. cellderiv)) call update_kgrid

                if (np .ne. 0) then ! begin change guiding
                    !         calcolo wavefunction di ogni walker before changing it
                    !       update all matrix changed by the v,vsz
                    do k = 1, nmat - 1
                        alphab(k) = 0.d0
                    end do
                    if ((yesmin .eq. 0 .or. .not. noproj) .and. yesupdate_ion) then
                        if (rank .eq. 0) write (6, *) ' Recomputing rpar '
                        call eval_iond(iond, rion, nion, LBox, psip, iond_cart)
                        psip(1:kp_ion) = rpar(1:kp_ion) ! store old rpar
                        call preprpar(rpar, kp_ion, iond, nion, kiontotj                &
                             &, nozeroj_c, nnozeroj_c, nelorbj_c, jbraj, iesfree, indfree             &
                             &, jbrajsz, iesinv, indinv, orbcostn, kiontot, nozero_c, nnozero_c        &
                             &, nelorb_c, jbradet, iessw, indsw, adrlambda, whereiesm, iesm, whereiesup &
                             &, iesup, iond_cart, typeorb, nshellj_c, multj_c, ioccj_c, occj_c         &
                             &, jas_invariant, orbps)
                        !        Untouched parameters below a certain treshold
                        do ii = 1, kp_ion
                            if ((rpar(ii) .lt. rmaxj .or. rmaxj .eq. 0)                           &
                                 &.and. (ii .gt. endinv .and. rpar(ii) .gt. 0)) rpar(ii) = 0.d0
                            if ((rpar(ii) .lt. rmaxinv .or. rmaxinv .eq. 0)                       &
                                 &.and. ii .le. endinv) rpar(ii) = 0.d0
                            if ((abs(rpar(ii)) .lt. rmax .or. rmax .eq. 0) .and. rpar(ii) .lt. 0.d0)    &
                                 &   rpar(ii) = 0.d0
                        end do

                        !       Do not cancel  variational parameters that where previously
                        !       optimized if yescut* is false.
                        if (.not. yescutjas) then
                            do ii = 1, kp_ion
                                if (psip(ii) .eq. 0.d0 .and. rpar(ii) .gt. 0) rpar(ii) = 0.d0
                            end do
                        end if

                        if (.not. yescutdet) then
                            do ii = 1, kp_ion
                                if (psip(ii) .eq. 0.d0 .and. rpar(ii) .lt. 0) rpar(ii) = 0.d0
                            end do
                        end if
#ifdef PARALLEL
                        if (yesquantum .and. commrep_mpi .ne. commsr_mpi) then
                            do ii = 1, kp_ion
                                psip(kp_ion + ii) = abs(rpar(ii))
                            end do
                            ! If there is some accepted parameter in some bead accept the same for all
                            call mpi_allreduce(psip(kp_ion + 1), psip, kp_ion, MPI_DOUBLE_PRECISION&
                                 &, MPI_MIN, commsr_mpi, ierr)
                            do ii = 1, kp_ion
                                if (psip(ii) .eq. 0.d0) rpar(ii) = 0.d0
                            end do
                        end if
#ifdef UNREL_DIAG
                        ! set consistency among processors
                        call bcast_real(rpar, kp_ion, 0, commsr_mpi)
#endif
#endif

                        if (ncg_adr .gt. 0) then
                            !    The parameters have to be computed before changing reducel
                            call project_v(.true.)
                            call pareff(npar, initpar, nparsw, initparsw, nparinv, initparinv&
                                 &, endinv, ncg_adr, kp0, rpar, reducel, jas_invariant, adrlambda, nmax_ion&
                                 &, type_atom, allfit, orbps)
                            call project_alphavar
                        elseif (rmaxj .ne. 0 .or. rmaxinv .ne. 0 .or. rmax .ne. 0) then
                            !    no parametrization but locality
                            if (rank .eq. 0) write (6, *) ' Warning setting to zero unoptimazible > rmax '
                            call project_rmax
                        end if

                    else
                        !        Compute and write variational parameters (do not recompute mat)
                        if (ncg_adr .gt. 0) call project_v(.false.)
                        if (smoothcut .ne. 0.d0 .and. killcut) then
                            if (rank .eq. 0) write (6, *) &
                                 &' Warning damping  to zero unoptimazible > rmax by ', smoothcut
                            call project_rmax
                        end if
                    end if ! endif idyn>0

                    if ((yesmin .ne. 0 .or. iesup .ne. 0 .or. iessw .ne. 0) .and. contraction .ne. 0) then
                        ! update detmat
                        if (yesdetmatc) call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                             &, nelcol_c, detmat, detmat_c, mu_c, psip)

                        call update_projm
                    end if ! yesmin and contraction
                    !          make all meas and reinitialize normalization coefficients
                    pseudologic = pseudorandom
                    iesrandoml = iesrandoma
                    flagcont = .true.
                    singdet(1:in1) = .true.
                    time0 = cclock()
                    if (rank .eq. 0) write (6, *) ' Time around preprpar ', time0 - timep
                    timep = time0
                    ! Updating matrix before recomputing by scratch.
#ifdef _OFFLOAD
!$omp target update to (jasmat,muj_c,jasmat_c,detmat,detmat_c,projm,mu_c,eagp_pfaff)
#endif
                    call makeallmeas_fast
                    time0 = cclock()
                    if (rank .eq. 0) write (6, *) ' Time makeallmeas_fast =', time0 - timep
                    time_meas = time_meas + time0 - timep
                end if ! end changhing GUIDING (np.ne.0)
            end if ! fine if (np.gt.0)

            !       END CHANGING THE GUIDING

            !          writing data only if necessary

            if (rank .eq. 0 .and. (ieskint .eq. 0 .or. (yeswrite12 .and. acc_dyn)) .and.&
                 &.not. nowrite12) then

                write (6, *) ' Warning writing WF '

                psip(nmat) = alphab(nmat)
                do k = 1, nmat - 1
                    psip(k) = alphavar(k)
                end do

                if (iessw .ne. 0 .or. ieskint .ne. 0) then

                    if (contraction .ne. 0) then

                        if (detc_proj) then
                            !               From real to effective
                            ! From real to effective
                            if (rank .eq. 0) write (6, *) ' Passi qui XXVI real-eff'
                            if (allowed_averagek) call attach_phase2det(.false., detmat_proj)

                            if (yes_complex) then
                                do ii = 1, nnozero_c
                                    psip(2*ii - 1 + inddsw - 1) = detmat_proj(2*nozero_c(ii) - 1)
                                    psip(2*ii + inddsw - 1) = detmat_proj(2*nozero_c(ii))
                                end do
                            else
                                do ii = 1, nnozero_c
                                    psip(ii + inddsw - 1) = detmat_proj(nozero_c(ii))
                                end do
                            end if
                   !!               Back to real
                            if (rank .eq. 0) write (6, *) ' Passi qui XXVII eff-real'
                            if (allowed_averagek) call attach_phase2det(.true., detmat_proj)
                        else
                   !!               From real to effective
                            if (rank .eq. 0) write (6, *) ' Passi qui XXIX real-eff'
                            if (allowed_averagek) call attach_phase2det(.false., detmat_c)
                            if (yes_complex) then
                                do ii = 1, nnozero_c
                                    psip(2*ii - 1 + inddsw - 1) = detmat_c(2*nozero_c(ii) - 1)
                                    psip(2*ii + inddsw - 1) = detmat_c(2*nozero_c(ii))
                                end do
                            else
                                do ii = 1, nnozero_c
                                    psip(ii + inddsw - 1) = detmat_c(nozero_c(ii))
                                end do
                            end if
                   !!               Back to real
                            if (rank .eq. 0) write (6, *) ' Passi qui XXVIII eff-real'
                            if (allowed_averagek) call attach_phase2det(.true., detmat_c)
                        end if ! endif detc_proj
                        ind = ipc*nnozero_c + inddsw - 1
                    else ! if contraction
                !!               From real to effective
                        if (rank .eq. 0) write (6, *) ' Passi qui XXX real-eff'
                        if (allowed_averagek) call attach_phase2det(.false., detmat)
                        if (yes_complex) then
                            do ii = 1, nnozero
                                psip(2*ii - 1 + inddsw - 1) = detmat(2*nozero(ii) - 1)
                                psip(2*ii + inddsw - 1) = detmat(2*nozero(ii))
                            end do
                        else
                            do ii = 1, nnozero
                                psip(ii + inddsw - 1) = detmat(nozero(ii))
                            end do
                        end if
                !!               Back to real
                        if (rank .eq. 0) write (6, *) ' Passi qui XXXI eff-real'
                        if (allowed_averagek) call attach_phase2det(.true., detmat)
                        ind = ipc*nnozero + inddsw - 1
                    end if ! endif contraction
                    if (npar_eagp .gt. 0) then
                        do ii = nnozero_c + 1, nnozero_c + nnozero_eagp
                            !                      do iy=1,ndiff
                            !                        do ix=iy+1,ndiff
                            iy = (nozero_c(ii) - 1)/ndiff + 1
                            ix = nozero_c(ii) - (iy - 1)*ndiff
                            if (ipc .eq. 1) then
                                ind = ind + 1
                                psip(ind) = eagp_pfaff(ix, iy)
                            else
                                ind = ind + 1
                                psip(ind) = eagp_pfaff(2*ix - 1, iy)
                                ind = ind + 1
                                psip(ind) = eagp_pfaff(2*ix, iy)
                            end if
                            !                        enddo
                        end do
                        !                     ind=ind+ipc*npar_eagp
                    end if

                else ! iessw>0

                    ind = inddsw - 1

                end if ! endif iessw

                !                  endif ! endif molopt

                if (iesup .ne. 0 .or. ieskint .ne. 0) then
                    !                                       write(6,*) ' dup_c written '

                    do ii = 1, iesup_read
                        psip(ii + ind) = dup_c(ii)
                        !                                       write(6,*) ii,dup_c(ii)
                    end do
                    ind = ind + iesup_read
                end if

                ! updating ionic positions
                if (ieskint .ne. 0) then
                    if (idyn .ne. 0) then
                        do ii = 1, nion
                            do kk = 1, 3
                                ind = ind + 1
                                psip(ind) = rion_write(kk, ii)
                            end do
                        end do

                        if (cellderiv) then
                            ind = ind + 1
                            psip(ind) = rs_write
                            ind = ind + 1
                            psip(ind) = celldm_write(2)
                            ind = ind + 1
                            psip(ind) = celldm_write(3)
                        end if
                    else
                        do ii = 1, nion
                            do kk = 1, 3
                                ind = ind + 1
                                psip(ind) = rion(kk, ii)
                            end do
                        end do

                        if (cellderiv) then
                            ind = ind + 1
                            psip(ind) = rs
                            ind = ind + 1
                            psip(ind) = celldm(2)
                            ind = ind + 1
                            psip(ind) = celldm(3)
                        end if
                    end if
                end if

                nmatb = ind
                !         alphasto eliminated
                ! short output with no wf opt
                if (iread .eq. 1) then
                    write (12) 1.d0, 1.d0
                else
                    write (12) 1.d0, 1.d0, (psip(k), k=1, nmatb)
                end if
                if (flush_write) flush (12)
                srforce = 0.d0
                srforcew = 0.d0
                wforce = 0.d0
                wforcew = 0.d0
            end if ! rank.eq.0

            ngg = ngg + 1
            if (ngentry .eq. -3) then
                nweight = rweight
                ngentry = 0
            end if

            if (mod(i_main/nweight, iskipdyn) .eq. iskipdyn - nmore_force) then
                inext = inext + nmore_force*nweight
                perbin = nmore_force
            else
                inext = inext + nweight
                perbin = 1
            end if

        else ! i.eq.inext

            if (pippo .eq. ibinit + 1) then
                avenernum = 0.d0
                avenerden = 0.d0
                tave_cyrus = 0.d0
                tcount_cyrus = 0.d0
                pippoc = 1
                call dscalzero((1 + ipc)*nbindim, 0.d0, efenergy, 1)
                call dscalzero(nbindim*ndimpdim*2, 0.d0, efp, 1)
                call dscalzero(nbindim*ieskindim*3, 0.d0, ef, 1)
                if (LBox .gt. 0.d0) call dscalzero(3*nbindim, 0.d0, efpress, 1)
            end if

            if (LBox .gt. 0.d0) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    efpress(1, indpippo) = derEV(cpippo - istm)*factorsr(cpippo - istm)*wcort   &
                         & + efpress(1, indpippo)
                    efpress(2, indpippo) = p_pulay(cpippo - istm)*factorsr(cpippo - istm)  &
                                          & *wcort + efpress(2, indpippo)
                    efpress(3, indpippo) = p_pulay(cpippo - istm)*factorsr(cpippo - istm)  &
                                          & *wcort*enert(1, cpippo - istm) + efpress(3, indpippo)
                end do
            end if

            if (ipc .eq. 2) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    efenergy(1:2, indpippo) = enert(1:2, cpippo - istm)                             &
                                             & *factorsr(cpippo - istm)*wcort + efenergy(1:2, indpippo)
                    efenergy(3, indpippo) = factorsr(cpippo - istm)*wcort + efenergy(3, indpippo)
                end do
            else
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    efenergy(1, indpippo) = enert(1, cpippo - istm)                             &
                                           & *factorsr(cpippo - istm)*wcort + efenergy(1, indpippo)
                    efenergy(2, indpippo) = factorsr(cpippo - istm)*wcort + efenergy(2, indpippo)
                end do
            end if

            if (ndimp .ge. 1) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    if (yes_correct) then

                        if (lrdmc_der .and. .not. lrdmc_nonodes) then
                            do kk = 1, ndimj
                                efp(kk, 1, indpippo) = econfh((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do
                            do kk = ndimjp, ndimp, 2
                                efp(kk, 1, indpippo) = econfh((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk + 1, 1, indpippo) = econfh(kk*in1 + cpippo - istm)              &
                                                          & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 1, indpippo)
                                efp(kk, 2, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                                efp(kk + 1, 2, indpippo) = econf(kk*in1 + cpippo - istm)              &
                                                          & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 2, indpippo)
                            end do

                        else
                            do kk = 1, ndimj
                                efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = enert(1, cpippo - istm)                              &
                                                      & *econf(cpippo - istm + (kk - 1)*in1)                                  &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do
                            do kk = ndimjp, ndimp, 2
                                efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk + 1, 1, indpippo) = econf(kk*in1 + cpippo - istm)              &
                                                          & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 1, indpippo)
                                efp(kk, 2, indpippo) = (enert(1, cpippo - istm)*econf(cpippo - istm + (kk - 1)*in1)&
                                     & - enert(2, cpippo - istm)*econf(cpippo - istm + kk*in1))&
                                     & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                                efp(kk + 1, 2, indpippo) = (enert(2, cpippo - istm)*econf(cpippo - istm + (kk - 1)*in1)&
                                     & + enert(1, cpippo - istm)*econf(cpippo - istm + kk*in1))&
                                     & *factorsr(cpippo - istm)*wcort + efp(kk + 1, 2, indpippo)
                            end do
                        end if
                    elseif (yes_real) then
                        !                  do kk=1,ndimj
                        !                     efp(kk,1,indpippo)=econf((kk-1)*in1+cpippo-istm)              &
                        !                          &*factorsr(cpippo-istm)*wcort+efp(kk,1,indpippo)
                        !                     efp(kk,2,indpippo)=enert(1,cpippo-istm)                              &
                        !                          &  *econf(cpippo-istm+(kk-1)*in1)                                  &
                        !                          &  *factorsr(cpippo-istm)*wcort+efp(kk,2,indpippo)
                        !                  enddo
                        !                  do kk=ndimjp,ndimp
                        do kk = 1, ndimp
                            efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                  & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                            efp(kk, 2, indpippo) = (enert(1, cpippo - istm)&
                                 & *econf(cpippo - istm + (kk - 1)*in1) + &
                                 &scaleeloc*econfh(cpippo - istm + (kk - 1)*in1))&
                                 & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                        end do
                    else

                        if (lrdmc_der .and. .not. lrdmc_nonodes) then

                            do kk = 1, ndimp
                                efp(kk, 1, indpippo) = econfh((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do

                        else
                            do kk = 1, ndimp
                                efp(kk, 1, indpippo) = econf((kk - 1)*in1 + cpippo - istm)              &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 1, indpippo)
                                efp(kk, 2, indpippo) = enert(1, cpippo - istm)                              &
                                                      & *econf(cpippo - istm + (kk - 1)*in1)                                  &
                                                      & *factorsr(cpippo - istm)*wcort + efp(kk, 2, indpippo)
                            end do
                        end if
                    end if
                end do
            end if

            if (ieskin .ge. 1) then
                do cpippo = ist, ien
                    indpippo = pippoc + nbinr*(cpippo - ist)
                    do kk = 1, ieskin
                        ! der energy
                        ef(kk, 1, indpippo) = econfion((kk - 1)*in1 + cpippo - istm)             &
                                             & *factorsr(cpippo - istm)*wcort + ef(kk, 1, indpippo)
                        ! lor Or
                        ef(kk, 2, indpippo) = econf(nwkin + (kk - 1)*in1 + cpippo - istm)          &
                                             & *factorsr(cpippo - istm)*wcort + ef(kk, 2, indpippo)
                        ! log Or*energy
                        ef(kk, 3, indpippo) = enert(1, cpippo - istm)                                &
                                             & *econf(nwkin + cpippo - istm + (kk - 1)*in1)                           &
                                             & *factorsr(cpippo - istm)*wcort + ef(kk, 3, indpippo)
                        !                      ef(kk,4,indpippo)=econf(nwkin+cpippo-istm+(kk-1)*in1)**2       &
                        !                           &   *factorsr(cpippo-istm)*wcort+ef(kk,4,indpippo)
                    end do
                    !        ef(kk,1,pippo)=ef(kk,1,pippo)/nw !der energy
                    !        ef(kk,2,pippo)=ef(kk,2,pippo)/nw !lor Or
                    !        ef(kk,3,pippo)=ef(kk,3,pippo)/nw !log Or*energy
                    ! end do for kk=1,ieskin
                end do
                ! endif ieskin
            end if

            pippo = pippo + 1
            !     Update the bin counter for the next measure
            if (mod(pippo - ibinit - 1 + ndone*perbin, lbin*perbin) .eq. 0) then
                if (pippo .ne. ibinit + 1) pippoc = pippoc + 1
                !               if(perbin.gt.1.and.rank.eq.0) write(6,*) 'updated pippoc',i,pippo,pippoc
            end if
            if (i_main .eq. inext + ibinit + 1 - nweight) then
                srforce = 0.d0
                srforcew = 0.d0
                wforce = 0.d0
                wforcew = 0.d0
            end if ! i.eq.inext+ibinit-nweight

            call reweight0(Nw, in1, np, npmn, factorsr                       &
                 &, ipsip, psip, alphab, sov, econf, econfh, econfion                &
                 &, ieskin, ncg, epsdgel, epstion, 0, wcort, itestrr, enert           &
                 &, scalpar, parcut, iflagerr, force, err, iesconv, itouch                 &
                 &, parcutmin, parcutpar, npbra, kl, etry                                &
                 &, epsi, tpar, beta, ist, ien, rank, ierr, parr, parcute            &
                 &, nbin, fk, dimfk, fkav, okav, skdiag, weightall, reduce, nbinmax, nweight, ibinit &
                 &, 19, stepcg, lwork, idyn, temp, weight_vir, friction, scalecov, delta0, delta0q, delta0k  &
                 &, dt, velion, ris, tmes, cov, rpar, npar, initpar, nparsw, initparsw &
                 &, nparinv, initparinv, endinv, iond, nion, adrlambda, rmax, 20&
                 &, writescratch, countscra, bufscra, jas_invariant, tolcg, ieskinion     &
                 &, rion, iespbc, atom_number, eps_dyn5, maxdev_dyn, acc_dyn, normcorr&
                 &, row_comm, row_id, yescomm)

#ifdef _OPENMP
            new_threads = omp_get_max_threads()
#else
            new_threads = 1
#endif

            if (old_threads .ne. new_threads) then
#ifdef UNREL_SMP
                write (6, *) 'Warning number of threads not conserved', new_threads, old_threads
#else
                write (6, *) 'ERROR in number of threads !!! ', new_threads, old_threads
                iflagerr = 1
#endif
            end if

        end if ! fine if i.eq.inext

        do j = 1, in1
            econfw(j) = econf(j)
        end do

        if (nfat .ne. 0) then
            indcor = mod(indcor, nfat) + 1
            wcort = wcort*wbra/wcorw(indcor)
            wcorw(indcor) = wbra
        end if

#ifdef __KCOMP
123     format(32767e15.7)
#else
123     format(1000000e15.7)
#endif
    end subroutine updatewfopt

    subroutine writeandbranch

        implicit none
        real(8) drand1, enercont, jacobian, mapping, wbra_av
        ! variables needed to collect quantitis from pools
        integer :: ikp, ist_kp, ien_kp, icdiff_tot
        logical pseudologic_fast, iesrandoml_fast

#if defined (_OPENMP) && defined (__NOOMP)
        integer, external :: omp_get_max_threads
        call omp_set_num_threads(1) ! scalar code
#endif

        if (iesbra) then

            timepp = cclock()
            !
            ! communication of the weights to the master
            ! In the case of decoupled VMC/DMC runs, communicate the
            ! weights to the master but within each processor pool. Steps:
            ! - create local variable wconfn_kps with dimension nw/nk
            ! - collect weights from the global wconfn within each pool
            ! - perform the branching internally to each pool
            ! - redistribute the local (to pool) weights to the global array wconfn
            !

#ifdef PARALLEL

            if (decoupled_run) then
                call gather_wconfn(ist, ien, wconfn, wconfn_kps, nw, nk)
            else
                psip(1:in1) = wconfn(ist:ien)
                call mpi_gather(psip, in1, MPI_DOUBLE_PRECISION, wconfn, &
                                in1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            end if

#endif

            ! only the master performs the branching
            ! for standard DMC calculation.
            if (.not. decoupled_run .and. rank .eq. 0) then
                if (nw_max .le. 0) then
                    zeta(nw + 1) = 1.d0 - drand1()
                    call branchingo(nw, wconfn, wbra, zeta, icdiff, ipsip, jbra)
                else
                    wbra_av = 0.d0
                    icdiff_tot = 0
                    do j = 1, nw, nw_max
                        zeta(nw_max + 1) = 1.d0 - drand1()
                        call branchingo(nw_max, wconfn(j), wbra, zeta, icdiff, ipsip, jbra(j))
                        wbra_av = wbra_av + wbra
                        icdiff_tot = icdiff_tot + icdiff
                        !  jbra output is from 1-nw_max but  should  point to j-->j+nw_max-1
                        do k = 0, nw_max - 1
                            jbra(j + k) = jbra(j + k) + j - 1
                        end do
                    end do
                    icdiff = icdiff_tot
                    wbra = wbra_av/dble(nw/nw_max)
                end if
                ! only the master within each pool perform the branching in
                ! the case of a decoupled k-points DMC calculation.
            elseif (decoupled_run .and. rankrep .eq. 0) then
                zeta(nwnkp) = 1.d0 - drand1()
                call branchingo(nwnk, wconfn_kps, wbra, zeta, icdiff, ipsip, jbra)
            end if

#ifdef PARALLEL
            if (decoupled_run) then
                call scatter_wconfn(ist, ien, wconfn, wconfn_kps, nw, nk)
                call mpi_bcast(wbra, 1, MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)

                !         deallocate(wconfn_kps)
            else
                psip(1:in1) = wconfn(ist:ien)
                call mpi_scatter(wconfn(ist), in1, MPI_DOUBLE_PRECISION, &
                                 psip, in1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                wconfn(ist:ien) = psip(1:in1)
                call mpi_bcast(wbra, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            end if
            if (decoupled_run) then
                jbra(nwnkp) = ist
                call bcast_integer(jbra, nwnkp, 0, commrep_mpi)
                !         setting the global address
                ipsip(1:nwnk) = jbra(1:nwnk) + jbra(nwnkp) - 1
                call mpi_allgather(ipsip, nwnk, MPI_INTEGER, jbra&
                     &, nwnk, MPI_INTEGER, commcolrep_mpi, ierr)

            else
                call bcast_integer(jbra, nw, 0, MPI_COMM_WORLD)
            end if

#ifdef UNREL
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif

#endif

            ! NB: no dumping in the case of decoupled VMC/DMC calculation
            if (iread .eq. 2 .or. iread .eq. 3) then

                if (mod(i_main, ifreqdump) .eq. 0) then
                    ! dump configurations
                    if (ipc .eq. 1) then
                        do kk = 1, in1
                            !wconfw(kk)=wconfn(kk+istm)*psisn(kk) ! Take the sign into account
                            !              wconfw(kk)=wconfn(kk+istm)*wf_sign(psisn(kk)) ! Take the sign into account
                            wconfw(kk) = factorsr(kk)*wf_sign(psisn(kk)) ! Take the sign into account
                            psilnw(kk) = psiln(kk)
                        end do
                    else
                        do kk = 1, in1
                            !              wconfw(kk)=wconfn(kk+istm) ! Do not take the sign into account
                            wconfw(kk) = factorsr(kk) ! Do not take the sign into account
                            psilnw(kk) = psiln(kk)
                            wconfw(kk + in1) = psisn(kk) ! Take the phase into account
                        end do
                    end if
                    do kk = 1, nw
                        ipsip(kk) = jbraw(jbra(kk))
                    end do
                    do kk = 1, nw
                        jbraw(kk) = ipsip(kk)
                    end do
                    wbraw = wbraw*wbra

                    if (iread .eq. 2) then
                        if (io_level .eq. 1) then
                            if (ipc .eq. 2) then
                                write (15) wbraw, (((real(kel(k, (j - 1)*nrnel + jj)), k=1, 3), jj=1, nel), j=1, in1)&
                                     &, (jbraw(k), k=ist, ien), (wconfw(k), wconfw(k + in1), k=1, in1)
                            else
                                write (15) wbraw, (((real(kel(k, (j - 1)*nrnel + jj)), k=1, 3), jj=1, nel), j=1, in1)&
                                     &, (jbraw(k), k=ist, ien), (wconfw(k), k=1, in1)
                            end if
                        end if
                    elseif (iread .eq. 3) then
                        ! correlated sampling (dump configurations, weight, energy, logpsi)
                        if (io_level .eq. 1) then
                            !               if(pseudorandom) then
                            ! dump also the random angles
                            if (ipc .eq. 2) then

                                write (15) wbraw, (((real(kel(k, (j - 1)*nrnel + jj)), k=1, 3), jj=1, nel), j=1, in1)        &
                                     &, (jbraw(k), k=ist, ien), (wconfw(k), wconfw(k + in1), k=1, in1), (econf(k), k=1, in1)     &
                                     &, (psilnw(k), k=1, in1), ((real(angle(jj, k)), jj=1, 18), k=1, nel*in1)
                            else
                                write (15) wbraw, (((real(kel(k, (j - 1)*nrnel + jj)), k=1, 3), jj=1, nel), j=1, in1)        &
                                     &, (jbraw(k), k=ist, ien), (wconfw(k), k=1, in1), (econf(k), k=1, in1)      &
                                     &, (psilnw(k), k=1, in1), ((real(angle(jj, k)), jj=1, 18), k=1, nel*in1)
                            end if

                            !               else
                            !                  if(ipc.eq.2) then
                            !                   write(15) wbraw,(((real(kel(k,(j-1)*nrnel+jj)),k=1,3),jj=1,nel),j=1,in1)        &
                            !                        &,(jbraw(k),k=ist,ien),(wconfw(k),wconfw(k+in1),k=1,in1),(econf(k),k=1,in1)&
                            !                        &,(psilnw(k),k=1,in1)
                            !                   else
                            !                   write(15) wbraw,(((real(kel(k,(j-1)*nrnel+jj)),k=1,3),jj=1,nel),j=1,in1)        &
                            !                        &,(jbraw(k),k=ist,ien),(wconfw(k),k=1,in1),(econf(k),k=1,in1)     &
                            !                        &,(psilnw(k),k=1,in1)
                            !                   endif
                            !                endif
#ifndef PARALLEL
                        end if
#else
                    else if (io_level .eq. 2) then
                        ! check the first write and set the blocksize
                        if (details_SP%view == MPI_DATATYPE_NULL) then
                            call mpiio_file_get_disp(details_SP)
                            !                  if(pseudorandom) then
                            SP_block_size = (1 + ipc + 21*nel)*in1
                            !                  else
                            !                    SP_block_size=(1+ipc+3*nel)*in1
                            !                  endif
                            call mpiio_file_create_view(details_SP, SP_block_size, MPI_REAL)
                            if (rank == 0) write (6, *) "mpiio: details SP part size", SP_block_size
                            call mpiio_file_reset_view(details_SP)
                            buffer_depth = 10.0/(cclock() - time1p) ! dump data every 10s
                            if (buffer_depth*SP_block_size*4 > 1048576) then
                                if (rank .eq. 0) write (6, *) "mpiio: buffer_size is limited below 1MB!"
                                buffer_depth = 1048576/(SP_block_size*4)
                            end if
                            if (buffer_depth < 1) buffer_depth = 1
                            call mpi_bcast(buffer_depth, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                            if (rank .eq. 0) write (6, *) "mpiio: buffer_depth", buffer_depth
                            if (rank .eq. 0) write (6, *) "mpiio: buffer_size (Byte)", buffer_depth*SP_block_size*4
                            buffer_counter = 0
                            allocate (SP_buffer(SP_block_size, buffer_depth))
                        end if
                        if (details_DP%view == MPI_DATATYPE_NULL) then
                            DP_block_size = 1 + in1*2
                            if (rank == 0) write (6, *) "mpiio: details DP part size", DP_block_size
                            call mpiio_file_get_disp(details_DP)
                            call mpiio_file_create_view(details_DP, DP_block_size, MPI_DOUBLE_PRECISION)
                            call mpiio_file_reset_view(details_DP)
                            allocate (DP_buffer(DP_block_size, buffer_depth))
                        end if

                        buffer_counter = buffer_counter + 1

                        kk = 0
                        do j = 1, in1
                            do jj = 1, nel
                                kk = kk + 3
                                SP_buffer(kk - 2:kk, buffer_counter) = real(kel(1:3, (j - 1)*nrnel + jj))
                            end do
                        end do
                        !               if(pseudorandom) then
                        do jj = 1, nel*in1
                            kk = kk + 18
                            SP_buffer(kk - 17:kk, buffer_counter) = real(angle(1:18, jj))
                        end do
                        !               endif

                        SP_buffer(kk + 1:kk + in1, buffer_counter) = real(jbraw(ist:ien))
                        if (ipc .eq. 2) then
                            !   The real part and the imaginary part are consecutive, as they belong to a walker.
                            call scopy(in1, wconfw, 1, SP_buffer(kk + in1 + 1, buffer_counter), 2)
                            call scopy(in1, wconfw(in1 + 1), 1, SP_buffer(kk + in1 + 2, buffer_counter), 2)
                        else
                            SP_buffer(kk + in1 + 1:kk + 2*in1, buffer_counter) = wconfw(1:in1)
                        end if

                        DP_buffer(1, buffer_counter) = wbraw
                        DP_buffer(2:in1 + 1, buffer_counter) = econf(1:in1)
                        DP_buffer(in1 + 2:in1*2 + 1, buffer_counter) = psilnw(1:in1)

                        if (buffer_counter .eq. buffer_depth) then
                            call MPI_File_write_all(details_SP%fp, SP_buffer&
                                &, SP_block_size*buffer_counter, MPI_REAL, status, ierr)
                            call MPI_File_write_all(details_DP%fp, DP_buffer&
                                &, DP_block_size*buffer_counter&
                                &, MPI_DOUBLE_PRECISION, status, ierr)
                            buffer_counter = 0
!                   if(rank.eq.0) write(6,*) "Writing the partial details"
                        end if
                        if (i_main .eq. (ngen + iend) .and. buffer_counter .ne. 0) then
                            call MPI_File_write_all(details_SP%fp, SP_buffer&
                                &, SP_block_size*buffer_counter, MPI_REAL, status, ierr)
                            call MPI_File_write_all(details_DP%fp, DP_buffer&
                                &, DP_block_size*buffer_counter&
                                &, MPI_DOUBLE_PRECISION, status, ierr)
                            deallocate (SP_buffer)
                            deallocate (DP_buffer)
                            buffer_counter = 0
!                   if(rank.eq.0) write(6,*) "Writing the last part of details"
                        end if
                    end if ! endif io_level
#endif

                end if ! endif iread=2

                ! reset cumulative wconf and jbra after dumping
                do kk = 1, nw
                    wconfw(kk) = 1.d0
                    jbraw(kk) = kk
                end do
                wbraw = 1.d0

            else ! if (mod(ifreqdump...)

                ! update cumulative wconf and jbra
                do j = 1, nw
                    wconfw(j) = wconfw(j)*wbra
                    ipsip(j) = jbraw(jbra(j))
                end do
                do j = 1, nw
                    jbraw(j) = ipsip(j)
                end do
                wbraw = wbraw*wbra
            end if !  if(mod(i,ifreqdump).eq.0)

        end if ! iread.eq.2  or iread=3

        ! parallel reshuffling cccccccccccccccccccccccccccccccccccccccccc
#ifdef _OFFLOAD
!$omp target update from (winv,winvj,ainv,winvbar,winvjbar,winvfn,winvbarfn)   if(.not.yes_fastbranch.and.yes_ontarget)
#endif

        if (noblocking) then

            call reshuffhub_noblock(Lztab, Lztabr, Ltab, Ltabb, nelnion, nw, npf, jbra, kel         &
                 &, dist, econf, table, tabler, tabpip, winv, nel2wt, winvj, nel2wtj, winvup  &
                 &, nel2upt, winvdo, nel2dot, ainv, nel2up, wsto, diagfn                   &
                 &, diag, psiln, psisn                                                 &
                 &, enert, vpot, vpotreg, enertrue, diffuse, nelkel, tmu, naccm, winvbar, nel2bar&
                 &, winvjbar, winvjbarsz, nel2jbar, nel2jbarsz, ivic, pseudolocal, nel     &
                 &, indt, npsa, gradtot, gradtotbar, angle, gradpsi, gradpsibar            &
                 &, n_gvec, sum_q_cos_gr, sum_q_sin_gr, rank, nproc, ierr, status, psip      &
                 &, skipreshuff, iessz, Lbox, psidetln, jastrowall_ee, dimee              &
                 &, jastrowall_ei, dimei, indtm, yesivic, vcut, diffkin                   &
                 &, winvfn, nel2wtfn, winvbarfn, nel2barfn, vpotsav_ee, nelsquare)

        else

            call reshuffhub(Lztab, Lztabr, Ltab, Ltabb, nelnion, nw, npf, jbra, kel         &
                 &, dist, econf, table, tabler, tabpip, winv, nel2wt, winvj, nel2wtj, winvup  &
                 &, nel2upt, winvdo, nel2dot, ainv, nel2up, wsto, diagfn                   &
                 &, diag, psiln, psisn                                                 &
                 &, enert, vpot, vpotreg, enertrue, diffuse, nelkel, tmu, naccm, winvbar, nel2bar&
                 &, winvjbar, winvjbarsz, nel2jbar, nel2jbarsz, ivic, pseudolocal, nel     &
                 &, indt, npsa, gradtot, gradtotbar, angle, gradpsi, gradpsibar            &
                 &, n_gvec, sum_q_cos_gr, sum_q_sin_gr, rank, nproc, ierr, status, psip      &
                 &, skipreshuff, iessz, Lbox, psidetln, jastrowall_ee, dimee              &
                 &, jastrowall_ei, dimei, indtm, yesivic, vcut, diffkin                   &
                 &, winvfn, nel2wtfn, winvbarfn, nel2barfn, vpotsav_ee, nelsquare)

            if (yes_fastbranch) then
                flagcont = .true.
                pseudologic_fast = pseudologic
                iesrandoml_fast = iesrandoml
                pseudologic = .false.
                iesrandoml = .false.
                call makeallmeas_fast
                pseudologic = pseudologic_fast
                iesrandoml = iesrandoml_fast
            end if

        end if
#ifdef _OFFLOAD
!$omp target update to (winv,winvj,ainv,winvbar,winvjbar,winvfn,winvbarfn)   if(.not.yes_fastbranch.and.yes_ontarget)
#endif

        ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        iflagerr = 0 ! by sure all walkers are with non zero det

        time_branch = time_branch + cclock() - timepp

        elseif (iread .eq. 2 .or. iread .eq. 3) then ! no branching

        if (mod(i_main, ifreqdump) .eq. 0) then
            ! dump configurations
            if (ipc .eq. 1) then
                do kk = 1, in1
                    !wconfw(kk)=wconfn(kk+istm)*psisn(kk)
                    wconfw(kk) = wconfn(kk + istm)*wf_sign(psisn(kk))
                    psilnw(kk) = psiln(kk)
                    econfw(kk) = econf(kk)
                end do
            else
                do kk = 1, in1
                    !wconfw(kk)=wconfn(kk+istm)*psisn(kk)
                    wconfw(kk) = wconfn(kk + istm)
                    wconfw(kk + in1) = psisn(kk)
                    psilnw(kk) = psiln(kk)
                    econfw(kk) = econf(kk)
                end do
            end if

            !         if(iesbra) then
            !         do kk=1,nw
            !         ipsip(kk)=jbraw(jbra(kk))
            !         enddo
            !         do kk=1,nw
            !         jbraw(kk)=ipsip(kk)
            !         enddo
            !         endif
            if (iread .eq. 2) then
                if (io_level .eq. 1) then
                    write (15) wbra, (((real(kel(k, (j - 1)*nrnel + jj)), k=1, 3), jj=1, nel), j=1, in1)       &
                         &, (jbraw(k), k=ist, ien), (wconfw(k), k=1, ipc*in1)
                end if
            elseif (iread .eq. 3) then
                ! correlated sampling (dump configurations, weight, energy, logpsi)
                if (io_level .eq. 1) then
                    !               if(pseudorandom) then
                    ! dump also the random angles
                    if (ipc .eq. 1) then
                        write (15) wbra, (((real(kel(k, (j - 1)*nrnel + jj)), k=1, 3), jj=1, nel), j=1, in1)        &
                             &, (jbraw(k), k=ist, ien), (wconfw(k), k=1, in1), (econf(k), k=1, in1)      &
                             &, (psilnw(k), k=1, in1), ((real(angle(jj, k)), jj=1, 18), k=1, nel*in1)
                    else
                        write (15) wbra, (((real(kel(k, (j - 1)*nrnel + jj)), k=1, 3), jj=1, nel), j=1, in1)        &
                             &, (jbraw(k), k=ist, ien), (wconfw(k), wconfw(k + in1), k=1, in1), (econf(k), k=1, in1)      &
                             &, (psilnw(k), k=1, in1), ((real(angle(jj, k)), jj=1, 18), k=1, nel*in1)
                    end if
                    if (flush_write) flush (15)
                    !               else
                    !                  write(15) wbra,(((real(kel(k,(j-1)*nrnel+jj)),k=1,3),jj=1,nel),j=1,in1)        &
                    !                       &,(jbraw(k),k=ist,ien),(wconfw(k),k=1,in1),(econf(k),k=1,in1)      &
                    !                       &,(psilnw(k),k=1,in1)
                    !               endif
#ifndef PARALLEL
                end if
#else
            else if (io_level .eq. 2) then
                ! check the first write and set the blocksize
                if (details_SP%view == MPI_DATATYPE_NULL) then
                    call mpiio_file_get_disp(details_SP)
                    !                  if(pseudorandom) then
                    SP_block_size = (1 + ipc + 21*nel)*in1
                    !                  else
                    !                    SP_block_size=(1+ipc+3*nel)*in1
                    !                  endif
                    call mpiio_file_create_view(details_SP, SP_block_size, MPI_REAL)
                    if (rank == 0) write (6, *) "mpiio: details SP part size", SP_block_size
                    call mpiio_file_reset_view(details_SP)
                    buffer_depth = 10.0/(cclock() - time1p) ! dump data every 10s
                    if (buffer_depth*SP_block_size*4 > 1048576) then
                        if (rank .eq. 0) write (6, *) "mpiio: buffer_size is limited below 1MB!"
                        buffer_depth = 1048576/(SP_block_size*4)
                    end if
                    if (buffer_depth < 1) buffer_depth = 1
                    call mpi_bcast(buffer_depth, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                    if (rank .eq. 0) write (6, *) "mpiio: buffer_depth", buffer_depth
                    if (rank .eq. 0) write (6, *) "mpiio: buffer_size (Byte)", buffer_depth*SP_block_size*4
                    buffer_counter = 0
                    allocate (SP_buffer(SP_block_size, buffer_depth))
                end if
                if (details_DP%view == MPI_DATATYPE_NULL) then
                    DP_block_size = 1 + in1*2
                    if (rank == 0) write (6, *) "mpiio: details DP part size", DP_block_size
                    call mpiio_file_get_disp(details_DP)
                    call mpiio_file_create_view(details_DP, DP_block_size, MPI_DOUBLE_PRECISION)
                    call mpiio_file_reset_view(details_DP)
                    allocate (DP_buffer(DP_block_size, buffer_depth))
                end if

                buffer_counter = buffer_counter + 1

                kk = 0
                do j = 1, in1
                    do jj = 1, nel
                        kk = kk + 3
                        SP_buffer(kk - 2:kk, buffer_counter) = real(kel(1:3, (j - 1)*nrnel + jj))
                    end do
                end do
                !               if(pseudorandom) then
                do jj = 1, nel*in1
                    kk = kk + 18
                    SP_buffer(kk - 17:kk, buffer_counter) = real(angle(1:18, jj))
                end do
                !               endif

                SP_buffer(kk + 1:kk + in1, buffer_counter) = real(jbraw(ist:ien))
                if (ipc .eq. 2) then
                    !   The real part and the imaginary part are consecutive, as they belong to a walker.
                    call scopy(in1, wconfw, 1, SP_buffer(kk + in1 + 1, buffer_counter), 2)
                    call scopy(in1, wconfw(in1 + 1), 1, SP_buffer(kk + in1 + 2, buffer_counter), 2)
                else
                    SP_buffer(kk + in1 + 1:kk + 2*in1, buffer_counter) = wconfw(1:in1)
                end if

                DP_buffer(1, buffer_counter) = wbra
                DP_buffer(2:in1 + 1, buffer_counter) = econf(1:in1)
                DP_buffer(in1 + 2:in1*2 + 1, buffer_counter) = psilnw(1:in1)

                if (buffer_counter .eq. buffer_depth) then
                    call MPI_File_write_all(details_SP%fp, SP_buffer&
                        &, SP_block_size*buffer_counter, MPI_REAL, status, ierr)
                    call MPI_File_write_all(details_DP%fp, DP_buffer&
                        &, DP_block_size*buffer_counter, MPI_DOUBLE_PRECISION&
                        &, status, ierr)
                    buffer_counter = 0
!             if(rank.eq.0) write(6,*) "Writing the partial details"
                end if

                if (i_main .eq. (ngen + iend) .and. buffer_counter .ne. 0) then
                    call MPI_File_write_all(details_SP%fp, SP_buffer&
                        &, SP_block_size*buffer_counter, MPI_REAL, status, ierr)
                    call MPI_File_write_all(details_DP%fp, DP_buffer&
                        &, DP_block_size*buffer_counter, MPI_DOUBLE_PRECISION&
                        &, status, ierr)
                    deallocate (SP_buffer)
                    deallocate (DP_buffer)
                    buffer_counter = 0
!             if(rank.eq.0) write(6,*) "Writing the last part of details"
                end if
            end if ! endif io_level
#endif
        end if ! elseif iread.eq.3

        ! reset cumulative wconf and jbra after dumping
        do kk = 1, nw
            wconfw(kk) = 1.d0
            jbraw(kk) = kk
        end do
        !      elseif(iesbra) then
        !         ! update cumulative wconf and jbra
        !         do j=1,nw
        !            wconfw(j)=wconfw(j)*wbra
        !            ipsip(j)=jbraw(jbra(j))
        !         enddo
        !         do j=1,nw
        !         jbraw(j)=ipsip(j)
        !         enddo
        end if !  if(mod(i,ifreqdump).eq.0)

        !   elseif(iread.ge.6.and.rank.eq.0) then
        !      write(15) wbra,(jbra(k),k=1,nw),(wconfw(k),k=1,4*nw)&
        !           &,(real(ris(2)*enertrue(k)),k=1,nw),(real(ris(2)*enertrue(k)),k=1,nw)

        end if !  if iesbra

        if (rank .eq. 0) sumdiff = sumdiff + icdiff

        if (itestr .ne. -5) then

            ngg = ngg + 1

            !----------- writing fort.12 and fort.12.fn -----------!

            if (decoupled_run) then
                call gather_avgs(wbra, ener, etot, wtotf, wbra_t, ener_t, etot_t, wtotf_t, np3, nk)
            end if

            ! in the case of decoupled VMC/LRDMC runs write a line for each k-point.
            if (rank .eq. 0) then
                if (decoupled_run) then
                    ! in this case write quantities for each k-point!
                    do ikp = 1, nk
                        if (itest .ne. 2) then
                            ! dmclrdmc and optimization
                            write (12) wbra_t(ikp), wtotf_t(2, ikp)&
                                &, wtotf_t(1, ikp), etot_t(1, ikp)&
                                &, ener_t(ikp)*ris(2)&
                                &, (etot_t(k, ikp), k=2, np3)
                        else
                            ! vmc
                            write (12) wbra_t(ikp), wtotf_t(1, ikp), ener_t(ikp)*ris(2), (etot_t(k, ikp), k=1, np3)
                        end if
                    end do
                else
                    if (itest .ne. 2) then
                        ! dmclrdmc and optimization
                        write (12) wbra, wtotf(2), wtotf(1), etot(1), ener*ris(2), (etot(k), k=2, np3)
                    else
                        ! vmc
                        write (12) wbra, wtotf(1), ener*ris(2), (etot(k), k=1, np3)
                    end if
                end if
                if (flush_write) flush (12)
            end if
            !      if(decoupled_run) deallocate(wbra_t,etot_t,ener_t,wtotf_t)

        else
            if (iread .ne. 2 .and. iread .ne. 3) ngn = ngn + 1
            if (rank .eq. 0) write (13) wbra, wtotf(1), ener*ris(2), (etot(k), k=1, np3), tmes, derEVp
        end if

        if (iread .eq. 2 .or. iread .eq. 3) ngn = ngn + 1

#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
    end subroutine writeandbranch

    subroutine Initializeall
        use allio
        use qmckl
        ! by E. Coccia (9/11/10)
        use extpot, only: mm_restr, n_x, n_y, n_z, delta, x0, ext_pot, link_atom, write_rwalk
        use splines, only: bscoef
! by E. Coccia (23/12/10)
! by E. Coccia (4/2/11)
        use van_der_waals, only: vdw

        implicit none
        real*8, external :: dnrm2
        real(8) drand1, enercont, jacobian, mapping, sumdet, arg1, arg2
        real(8) rion_ref(3), mind(3)
        real, external :: ran
        real(8), external :: atom_weight

        integer ind, ii, jj, nmoltry, ref_comp, read_seeds, read_seeds_mpiio, ntpar, inddet, min_pointvj

        logical yes_ipreshuff, yeschange_nweight, read_ok
#if defined (_OPENMP) && defined (__NOOMP)
        integer, external :: omp_get_max_threads
        call omp_set_num_threads(1) ! scalar code
#endif

#ifdef _QMCKL
        qmckl_ctx = qmckl_context_create()
        if (qmckl_ctx.eq.0_8) then
            write (0,*) "QMCKL context is a null pointer, but it should never happen"
            stop 1
        else
            write (6,*) "QMCKL context created"
        end if
#endif

!  Definition once for all machine precision and safemin as in lapack (more strict)
        epsmach = dlamch('e')
        safemin = 10000.d0*dlamch('s')
!          Initialize flags error
        yesdft = .false.
        iflagerr = 0
!          iflagerr is the local processor flag for error.
!          Each subroutine does not make anything if iflagerr is non zero.
!          Inside any non parallel subroutine iflag is set to 1 if the
!          subroutine went in error. In the output, the subroutine "checkiflagerr"
!          will make an mpi_allreduce of iflagerr summing to iflagerrall.
!          All processor will stop if iflagerrall>0, namely the run
!          will continue if and only if all processor have not found any error.
!
        call get_dir(path)
        if (rank == 0) write (6, *) ' Initial path : ', trim(path)

! Read all input cards excpet &molecul from standard input.
! Set default values of main quantities.
        call read_datasmin
        yesmin_read = .false.
        if (molopt .ne. 0) yesmin_read = .true.

#ifdef _CUSOLVER
#ifdef RISC
        call cusolver_handle_init_(handle)
#else
        call cusolver_handle_init(handle)
#endif
#endif

#ifdef PARALLEL
        row_comm = MPI_COMM_WORLD
        col_comm = MPI_COMM_WORLD
        if ((yesquantum .or. kaverage) .and. in1 .ne. 1) then
            call error(' Initializeall ', ' Multiple walkers per processor &
                 &are not allowed in the kaverage/quantum cases. Set nw = nproc !!! ', 1, rank)
        end if
#endif
        row_id = 0
        col_id = rank
        yescomm = .true.

!       Default values of communicators/id
#ifdef  PARALLEL
        commopt_mpi = MPI_COMM_WORLD
        rankopt = rank
        nprocopt = nproc

        commrep_mpi = MPI_COMM_WORLD
        rankcolrep = rank
        rankrep = rank
        nprocrep = nproc

        row_comm = MPI_COMM_WORLD
        mcol = nproc

        commsr_mpi = MPI_COMM_WORLD
        commcolsr_mpi = 0
        nprocsr = nproc
        ranksr = rank

        commcov_mpi = MPI_COMM_WORLD
        nproccov = nproc

!       commjas_mpi=MPI_COMM_WORLD
!       nprocjas=nproc
#else
        commopt_mpi = 0
        rankopt = rank
        nprocopt = nproc

        commrep_mpi = 0
        rankcolrep = rank
        rankrep = rank
        nprocrep = nproc

        row_comm = 0
        mcol = nproc

        commsr_mpi = 0
        commcolsr_mpi = 0
        nprocsr = nproc
        ranksr = rank

        commcov_mpi = 0
        nproccov = nproc

        comm_col = 0
        comm_raw = 0
        rankraw = 0
        rankcol = 0
!       commjas_mpi=0
!       nprocjas=nproc
#endif

        if (yesquantum .and. nbead .eq. nproc .and. nrep_bead .eq. 1) yescomm = .false.
        if (kaverage) yescomm = .false.

#ifdef PARALLEL

        if (yesquantum .and. kaverage) then
            call error(' Initializeall ', ' Quantum molecular dynamics &
                 & and k-average together is not implemented !!', 1, rank)
        end if
        if (kaverage .and. nrep_bead .gt. 1) nrep_bead = 1 ! to be safe
!
! split MPI_COMM_WORLD in rows/columns groups if:
! 1) quantum calculations with multiple beads
! 2) k-points calculation with multiple k-points
!                  1)                            2)
        if ((yesquantum .and. nbead .gt. 1) .or. kaverage) then
            !
            ! define mcol = # of rows within the 2D communicators grid
            !
            if (yesquantum) then ! quantum MD
                mcol = nproc/(nbead/nrep_bead)
                if (mcol*(nbead/nrep_bead) .ne. nproc .or. (nbead/nrep_bead)*nrep_bead .ne. nbead) then
                    call error(' Initializeall ', ' # of processors must be multiple of &
                         &  the # of beads !!', 1, rank)
                end if
            else ! k-points sampling
                mcol = nproc/nk
                if (mcol*nk .ne. nproc) then
                    if (rank .eq. 0) write (6, *) ' # processors / # k-points :', nproc, nk
                    call error(' Initializeall ', ' # of processors must be multiple of &
                         &  the # of k-points !!', 1, rank)
                end if
            end if

            ! split the communicators
            irow = rank/mcol
            jcol = mod(rank, mcol)

            call mpi_comm_split(mpi_comm_world, irow, jcol, row_comm, ierr)
            call mpi_comm_split(mpi_comm_world, jcol, irow, col_comm, ierr)
            call mpi_comm_rank(row_comm, row_id, ierr)
            call mpi_comm_rank(col_comm, col_id, ierr)
            call mpi_barrier(MPI_COMM_WORLD, ierr) ! for unreliable networks

            ! number of processor per column
            if (yesquantum) then
                nproccolrep = nbead
            elseif (kaverage) then
                nproccolrep = nk
            end if

            ! NB: in the case of k-points nrep_bead is always .eq. 1
            if (nrep_bead .gt. 1) then
                mcol_rep = mcol/nrep_bead ! i.e. = nproc/nbead per bead
                irow = rank/mcol_rep
                nprocrep = mcol_rep
                jcol = mod(rank, mcol_rep)
                call mpi_comm_split(mpi_comm_world, irow, jcol, commrep_mpi, ierr)
                call mpi_comm_split(mpi_comm_world, jcol, irow, commcolrep_mpi, ierr)
                call mpi_comm_rank(commrep_mpi, rankrep, ierr)
                call mpi_comm_rank(commcolrep_mpi, rankcolrep, ierr)
            else
                mcol_rep = mcol ! i.e. = nproc/(nbead/nrep_bead) per bead group
                commrep_mpi = row_comm
                nprocrep = mcol
                commcolrep_mpi = col_comm
                rankrep = row_id
                rankcolrep = col_id
            end if

            !         row_id = column number
            !         col_id = row number
            if (.not. yesavopt) then
                nprocopt = mcol
                commopt_mpi = row_comm
                commcolopt_mpi = col_comm
                rankopt = row_id
            end if
            !         if(.not.yesavcov) then
            !         Always the covariance is averaged only in the consistent pool
            commcov_mpi = commrep_mpi
            nproccov = mcol_rep
            !         endif
            if (.not. yesavsr) then
                commsr_mpi = row_comm
                commcolsr_mpi = col_comm
                nprocsr = mcol
                ranksr = row_id
            end if

        end if ! END if(yesquantum.and.nbead.gt.1 .or. kaverage
#else
        nprocu = 1
        nproc_diag = 1
#endif

#ifdef  PARALLEL
        if (prep .gt. 0 .and. mod(nprocsr, prep) .ne. 0) then
            !       Choose the closest divisor of nprocsr
            prep_try = 1
            do k = 2, 2*prep
                if (mod(nprocsr, k) .eq. 0) prep_try = k
            end do
            if (prep_try .gt. 1) then
                prep = prep_try
            else
                prep = -1
            end if
            if (rank .eq. 0) write (6, *) ' Warning prep should be a divisor of', nprocsr, ' Chosen =', prep
        end if
        if (prep .gt. 0) then
            ! another split
            irow = ranksr/prep
            jcol = mod(ranksr, prep)
            call mpi_comm_split(commsr_mpi, irow, jcol, comm_raw, ierr)
            call mpi_comm_split(commsr_mpi, jcol, irow, comm_col, ierr)
            call mpi_comm_rank(comm_raw, rankraw, ierr)
            call mpi_comm_rank(comm_col, rankcol, ierr)
        end if
#endif

!     ! define new rank/nw/mpi_comm for computing the energy within each pool
!     ! in the case of a decoupled k-point calculation.
!     rank_av = rank
!     nw_av   = nw
!     comm_av = commsr_mpi ! = MPI_COMM_WORLD in the case of yesavsr=.true. as for kaverage.
!     if(kaverage.and.decoupled_run) then
!        rank_av = rankrep
!        nw_av   = nw/nk ! nw is always multiple of nk
!        comm_av = commrep_mpi
!     endif
!
        call check_scratch(rank, path, scratchpath)
        ranseedfilename = trim(scratchpath)//'randseed.'//trim(chara)
!
        iese_eff = min(iese, 3) ! no more than 3 averaged on time corr fun so far
        if (rank .eq. 0) write (6, *) ' iese_eff=', iese_eff
!
! reading pseudo potential if any (npsa>0 read by datasmin)
        call read_pseudo
!
!   if(alat.eq.0.d0) then
!      alat=1.d0
!      if(rank.eq.0) write(6,*) ' warning alat set to one ',alat
!   endif
        if (rank .eq. 0) then
            if (alat2 .ne. 0.d0) then
                write (6, *) ' lattice spacing a1 a2 =', abs(alat)              &
                     &, abs(alat*alat2)
            else
                write (6, *) ' Single mesh  =', abs(alat)
            end if
            write (6, *) ' scratch of determinant each ', nscra
        end if
! rsignr=0    fixed node
! rsignr=1     VMC ref.
        if (tstepfn .eq. 0.d0) then
            rsignr = -rsignr*(1.d0 + gamma)
        elseif (itestr .ne. -5 .or. itestr4 .lt. -10) then
            rsignr = tstepfn
        end if

        test_aad = .false. ! flag to include FN and jastrowall inside compute_eloc_logpsi

        if (abs(itestr) .eq. 1 .or. abs(itestr) .eq. 6 .or. itestr .eq. -2 .or. itestr .eq. -3 .or. itestr4 .lt. -10) then

            if (rsignr .ne. 0 .and. rank .eq. 0)                              &
                 &  write (6, *) ' Warning non standard Fixed node !!!!', rsignr
            if (rsignr .eq. 0 .and. rank .eq. 0)                              &
                 &  write (6, *) ' Standard Fixed node '

        end if

        if ((itestr .eq. -5 .and. itestr4 .ge. -10) .and. rsignr .lt. 0) then
            call checkiflagerr(1, rank, ' Error rsign >= 0 in this case !!!')
        end if

        iseed = abs(iseedr)

        ibinit = iboot

! setting main indices for
! energy and correlation functions

        nmat = np + 1
        ndim = np + 1
        npm = np + 1
        npmn = max(nbinmax + 1 - ieskin, 1)
        iesconv = 0
        iesdelay = 0

        if (npmn .gt. npm) npm = npmn

        if (npm .lt. 6) npm = 6

!   Define nprocu
        if (nproc_diag .eq. 0 .and. npmn .ge. 1024 .and. itestr .eq. -5) then
            nprocu = min(nprocopt, nprocrep)
            if (nprocu .le. 0) nprocu = 1
        else
            nprocu = 1
        end if
        if (nproc_diag .eq. 0) then
            nproc_diag = nprocu
            if (npmn/nproc_diag .lt. 32) nproc_diag = min(npmn/32, nproc_diag) ! the maximum we can parallelize efficiently.
            if (nproc_diag .le. 0) nproc_diag = 1
#ifdef __KCOMP
            if (nproc_diag .gt. 48) nproc_diag = (nproc_diag/48)*48
#else
            if (nproc_diag .gt. 32) nproc_diag = (nproc_diag/32)*32
#endif
        end if
        if (rank .eq. 0) write (6, *) "sub_comm_diag uses", nproc_diag, "processors"
        call mpi_sub_comm_create(commrep_mpi, nproc_diag, sub_comm_diag, ierr)
! create sub communicator for diagonalization
! protect the size of nproc_diag
        nmats = nmat - iesm - iesd
! Number of VMC parameters
        ndimp = np - ieskin
        ndimpdim = max(ndimp, 1)

        nwm = nw
        nws = in1
!  The GLOBAL variables are:
!  wconfn,jbra,zeta

! ############### ALLOCATION ############################
        np3p3 = np3 + 3
        np3p4 = np3 + 4
        np3p5 = np3 + 5
        np3p6 = np3 + 6
        np3p7 = np3 + 7
        np3p8 = np3 + 8
        np3p9 = np3 + 9
        np3p10 = np3 + 10
        np3p11 = np3 + 11
        allocate (alphavar(npm), wcorw(nfat + 1), alphab(npm))
        alphavar = 0.d0
        wcorw = 0.d0
        alphab = 0.d0
! array containing all correlation functions and O^k operators
! needed to perform the average.
        allocate (etot(npm), wtot(np3p11))
        etot = 0.d0
        wtot = 0.d0
        if (itestrr .eq. -4) then
            if (idyn .lt. 2) then
                allocate (cov(1))
                allocate (sov(npmn, npmn, 5))
            else
                allocate (sov(npmn, npmn, 5))
                allocate (cov(max(ieskin*ieskin, 1)))
                if (idyn .eq. 5 .and. addrognoso) then
                    allocate (cov_old(max(ieskin*ieskin, 1)))
                    cov_old = 0.d0
                end if
            end if
        else
            if (idyn .lt. 2) then
                allocate (cov(1))
                allocate (sov(npmn, npmn, 3))
            else
                allocate (sov(npmn, npmn, 4))
                allocate (cov(max(ieskin*ieskin, 1)))
                if (idyn .eq. 5 .and. addrognoso) then
                    allocate (cov_old(max(ieskin*ieskin, 1)))
                    cov_old = 0.d0
                end if
            end if
        end if
        cov = 0.d0
        sov = 0.d0
! ########################################################

        do i = 1, nmat
            alphab(i) = 0.d0
        end do
        alphab(nmat) = 1.d0
        cost = 1.d0
        np3m = (np3 - 1)/nmat + 1

        if (rank .eq. 0) then
            write (6, *) ' Number of corr functions written =', np3m*nmat
            write (6, *) ' iopt =', iopt
        end if

        ngn = 0
        ngg = 0

        npp = np + 1
        npf = np

        if (npbra .gt. np) then
            call checkiflagerr(1, rank, ' npbra <= np !!!')
        end if

        nwm2 = 0
        nwdim = 0
        nwfix = 0
        nwinv = 0
        nwnel = 0
        nwrep = 0
        nwdw = 0
        nwfree = 0
        nwsw = 0
        nwkin = 0
        nwup = 0
        nwking = 0
        nprest = np
        nindt = 0
#ifdef PARALLEL
        if (nproc .gt. 1) then
            call rand_init(iseed)
            !    define  random number with certainly different initial iseed
            allocate (ipsip(2*nproc), psip(2*nproc))
            ipsip = 0
            do i = 1, 2*nproc
                irstart = (1.d0 - drand1())*2**29
                psip(i) = 2*irstart + 1
            end do
            call dsortx(psip, 1, 2*nproc, ipsip)
            !   To be consistent with the master
            call bcast_real(psip, 2*nproc, 0, MPI_COMM_WORLD)

            !     Disregard the lowest iseed (e.g. drand1=1)
            if (psip(1) .eq. 1.d0) then
                j = 2
                do while (psip(j) .eq. psip(1))
                    j = j + 1
                end do
            else
                j = 1
            end if
            iseed = psip(j)
            i = 1
            do while (i .le. rank .and. j .lt. 2*nproc)
                j = j + 1
                iseed = psip(j)
                if (psip(j) .ne. psip(j - 1)) i = i + 1
            end do
            iflagerr = 0
            if (i .ne. rank + 1) iflagerr = 1
            call mpi_allreduce(iflagerr, iflagerrall, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (iflagerrall .ne. 0) then
                if (rank .eq. 0) write (6, *) ' ERROR in random seed generation, try with different iseed !!! '
                deallocate (ipsip, psip)
                call mpi_finalize(ierr)
                stop
            end if
            deallocate (ipsip, psip)
        else
            irstart = iseed
            iseed = 2*irstart + 1
        end if
#else
        irstart = iseed
        iseed = 2*irstart + 1
#endif
        if (rank .eq. 0) write (6, *) ' initial iseed =', rank, iseed
!     write(6,*) ' Chosen iseed =',rank,iseed
        call rand_init(iseed)
        if (nproc .eq. 1) then ! for tests
            if (rank .eq. 0) write (6, *) ' initial random number =', drand1()
            call rand_init(iseed)
        end if

        if (itestr .eq. -5) then
            skipforce = 1
        else
            skipforce = 3
        end if

! open all files
        call open_files(rank, scratchpath)
        if (rankrep .eq. 0) open (unit=10, file='fort.10', form='formatted', position='REWIND')
!
! reading the wave function! Distinguish from single k-point or
! multiple k-points case.
! molyes = flag defined in read_fort10_fast subroutine which is .true. only if
! MOs are present and therefore the code assumes you are starting from DFT
! wavefunctions.
! NB this flag will be defined again within the subroutine read_fort10()

!
        if (kaverage .and. molyes .or. (yesquantum .and. yesread10)) then
            call read_fort10(27)
            ! check if input k-points/weights are correct
            if (.not. yesquantum) call check_kpoints(rankcolrep, commcolrep_mpi, rank)
        else
            !      In case K-average with AGP we have to read all geminal
            !      with the same wf. only change the phase

            if (kaverage) then
                call read_fort10(10)
            else
                commrep_mpi_sav = commrep_mpi
                rankrep_sav = rankrep
#ifdef PARALLEL
                commrep_mpi = mpi_comm_world
#else
                commrep_mpi = 0
#endif
                rankrep = rank
                call read_fort10(10)
                commrep_mpi = commrep_mpi_sav
                rankrep = rankrep_sav
            end if
        end if

        maxoutput = (dble(in1)*(2*dble(npdim) + ipc + 1) + 1)*dble(nweight)*8.d0/1d9
        if (rank .eq. 0) then
            if (writescratch .eq. 0) then
                write (6, *) ' Warning TurboRVB needs ', maxoutput, &
                     &' Gygabyte disc space per processor '
            else
                write (6, *) ' Warning TurboRVB needs ', maxoutput, &
                     &' Gygabyte RAM per processor '
            end if
        end if

        if (writescratch .ne. 0 .and. itestr .eq. -5) then
            allocate (bufscra(max((in1*(2*npdim + 1 + ipc) + 1)*nweight*nmore_force, 1)))
        else
            allocate (bufscra(1))
        end if
        bufscra = 0.d0

! set the index for choosing the k-point within a pool
        if (kaverage) then
            ikpoint = rankcolrep + 1
        else
            ikpoint = 1
        end if

        if (defparcutg .and. .not. fncont) then
            parcutg = 1
            if (rank .eq. 0) write (6, *) ' Default value of parcutg= ', parcutg
        end if
        if (rank .eq. 0) then
            if (.not. defparcutg .and. itest .eq. 1 .and. .not. fncont .and. n_body_on .ne. 0 .and. parcutg .le. 1) then
                write (6, *) ' Warning parcutg should be equal to 2 in this case !!! You are doing a test? '
            elseif (.not. defparcutg .and. itest .eq. 1 .and. .not. fncont .and. n_body_on .eq. 0 .and. parcutg .ne. 1) then
                write (6, *) ' Warning parcutg should be equal to 1 in this case !!! You are doing a test? '
            end if
        end if
        if (allfit) then
            ntpar = nmax_ion*nmax_ion*4
        else
            ntpar = (nmax_ion*(nmax_ion + 1))/2*4
        end if
        if (npower .gt. 0) then
            initpar = -2 - powermin
            npar = ntpar*npower
            if (rank .eq. 0) write (6, *) ' Default #parameters in the long range Jastrow', npar
        elseif (initpar .lt. -1) then
            if (mod(npar, ntpar) .ne. 0) then
                write (errmsg, *) 'ERROR Jastrow #parameters npar multiple of', ntpar
                call checkiflagerr(1, rank, errmsg)
            end if
        end if

        if (npowersz .gt. 0) then
            initparinv = -2 - powerminsz
            nparinv = ntpar*npowersz
            if (rank .eq. 0) write (6, *) ' Default #parameters in the long range spin Jastrow', nparinv
        elseif (initparinv .lt. -1) then
            if (mod(nparinv, ntpar) .ne. 0) then
                write (errmsg, *) 'ERROR spin Jastrow #parameters nparinv multiple of', ntpar
                call checkiflagerr(1, rank, errmsg)
            end if
        end if
!      Recomputing ncg_adr as it can be changed above
        ncg_adr = npar + nparsw + nparinv
        if (ncg_adr .gt. 0) then
            allocate (reducel(ncg_adr, kp0), v_adr(ncg_adr), mat_adr(ncg_adr, ncg_adr)&
                 &, ipip_adr(ncg_adr))
            reducel = 0.d0
            v_adr = 0.d0
            mat_adr = 0.d0
            ipip_adr = 0
            if (rank .eq. 0 .and. iopt .eq. 1) write (25, '(a1,10I7)') '#', ntpar, nmax_ion, npower&
                 &, npowersz, initpar, initparinv, initparsw, npar, nparinv, nparsw
        end if

        if (vdw) then
            call vdw_read()
        end if
! by E. Coccia (9/5/11): define QM-link and capping atoms
! in the Turbo list, and MM-link atoms in the CPMD list
! link_read reads the link.dat from CPMD containing
! the classical force field for angle and dihedrals
! involving the QM-link atom
        if (link_atom) then
            call link_read()
            ! by E. Coccia (31/5/11): exclusion list for vdW interactions
            call exclusion_list()
            ! by E. Coccia (10/5/11): the capping atom
            ! is not an independent object
            call r_capping()
        end if

        firstmol = 1
        if (contraction .gt. 0) then
            nmolfn = nelorb_c
        else
            nmolfn = ipf*nelorbh
        end if

        if (contraction .gt. 0) then
            allocate (allowcontr(nelorb_c, 2))
            allowcontr(:, :) = .true.
        end if
        if (contractionj .gt. 0) then
            allocate (allowcontrj(nelorbj_c, 2))
            allowcontrj(:, :) = .true.
        end if

        if (developer .eq. -1) then
            allocate (kelsav(3, nel))
            kelsav = 0.d0
        end if
!       stop for trivial input
        if (nelup .le. 0 .or. nelorb .le. 0 .or. nelup .lt. neldo .or. neldo .lt. 0) then

            if (nelup .eq. 0) write (6, *) ' No electrons up  ', nelup
            if (nelup .lt. neldo) write (6, *) ' Please put  #electrons up > #neldo ', nelup, neldo

            if (nelorb .le. 0) write (6, *) ' The electron basis is empty ', nelorb

            if (nelup .lt. 0) write (6, *) ' Negative # electrons up? ', nelup
            if (neldo .lt. 0) write (6, *) ' Negative # electrons down? ', neldo

            call checkiflagerr(1, rank, "Stop for trivial input!")

        end if

!  for FN  ip_reshuff can be put to 2 only for ndiff=0 unpolarized case.
        yes_ipreshuff = .false.
        if ((.not. symmagp .or. ipc .eq. 2)&
           & .and. ipf .eq. 1&
           & .and. (ndiff .eq. 0 .or. itest .eq. 2)&
           & .and. molecular .gt. 0&
           & .and. yesfast .ne. 0) yes_ipreshuff = .true.

        firstmolt = 1
        if (molecular .ne. 0) then
            firstmolt = nelorb_c - molecular + 1
            !       check if only the molecular diagonal are present
            firstmol = firstmolt
            cost = 0.d0
            !   check odd raws even columns
            if (yes_ipreshuff) then
                do i = firstmol + 1, nelorb_c - ndiff, 2
                    do j = firstmol, nelorb_c - ndiff, 2
                        cost = cost + abs(detmat_c(ipc*nelorb_c*(i - 1) + ipc*(j - 1) + 1))
                        if (ipc .ne. 1) cost = cost + abs(detmat_c(ipc*nelorb_c*(i - 1) + ipc*j))
                    end do
                end do
                !          check the atomic are zero
                do i = 1, nelorb_at
                    do j = 1, nelorb_at
                        cost = cost + abs(detmat_c(ipc*nelorb_c*(i - 1) + ipc*(j - 1) + 1))
                        if (ipc .ne. 1) cost = cost + abs(detmat_c(ipc*nelorb_c*(i - 1) + ipc*j))
                    end do
                end do
                if (cost .ne. 0) yes_ipreshuff = .false.
                if ((.not. symmagp .or. ipc .eq. 2) .and. molopt .eq. 0 .and. itestr .eq. -5 .and. iessw .gt. 0) then
                    yes_ipreshuff = .false.
                end if
            end if

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

            !       Update the last molecular orbitals according to the one to be optimized
            if (itestr .eq. -5) then
                call findmaxmat(iessw, nnozero_c, nozero_c, nelorb_c, jbradet, lastmol)
            end if
        elseif (contraction .ne. 0) then
            lastmol = 0
            if (ipc .eq. 2) then
                do j = 1, nelorb_c
                    do i = 1, nelorb_c
                        if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j&
                           & - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0&
                           & .and. max(i, j) .gt. lastmol) lastmol = max(i, j)
                    end do
                    do i = nelorb_c + 1, nelcol_c
                        if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0&
                           & .and. j .gt. lastmol) lastmol = j
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
            !       Update the last molecular orbitals according to the one to be optimized
            if (itestr .eq. -5) then
                call findmaxmat(iessw, nnozero_c, nozero_c, nelorb_c, jbradet, lastmol)
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

                if (ireadmin .eq. 1 .and. membig) yesfast = 0

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

            ! the basis is too small to be convenient
            if (nelorbh .le. 2*nmolfn&
               & .and. symmagp&
               & .and. ipf .eq. 1&
               & .and. ipc .eq. 1&
               & .and. ireadmin .eq. 0&
               & .and. membig) yesfast = 0

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
            end if

        end if

        if (contraction .ne. 0) then
#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
            call update_projm_
#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(1) ! restore the scalar
#endif
        end if

        ax = 0.d0
        ay = 0.d0
        az = 0.d0
        nx = 0
        ny = 0
        nz = 0
        nmolmax = 0
        nmolmin = 0
        noproj = .true.
        yesmin = 0
        weight_loc = -1.d0

        if (yesfast .eq. 1 .and. molecular .eq. 0) then
            call checkiflagerr(1, rank, ' ERROR yesfast=1 only with molecular orbitals!, Do not &
                 & define yesfast if you have no idea what this variable means !!! ')
        end if

!    in ANY event put to false this option.

! if molopt=0  optimize  in the usual way th
        if (molopt .ne. 0 .or. read_molecul) then

            call read_datasmin_mol

            !        Replacing value of nmol not assumed to change from input
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                if (nmol .ne. molecular - ndiff) then
                    nmol = molecular - ndiff
                    if (rank .eq. 0) write (6, *) ' Warning replacing nmol =', nmol
                end if
            else
                if (nmol .ne. (molecular - ndiff)/2) then
                    nmol = (molecular - ndiff)/2
                    if (rank .eq. 0) write (6, *) ' Warning replacing nmol =', nmol
                end if
            end if

            if (.not. molyes) then
                call checkiflagerr(1, rank, ' ERROR you cannot run |molopt|>0  without &
                     &molecular orbitals in your fort.10, use convertfort10mol.x for &
                     & converting it !!! ')
            end if

            if (.not. symmagp .or. ipc .eq. 2 .or. ipf .eq. 2) then
                lastmol = firstmolt + 2*nmolmatw - 1
            else
                lastmol = firstmolt + nmolmatw - 1
            end if

            if (forceyes) then

                if (yesmin .ne. 0 .and. ireadmin .ne. 1) then
                    yesfast = 1 ! the fastest code
                    firstmol = nelorb_c - molecular + 1
                    nmolfn = lastmol - firstmol + 1
                    if (2*nmolfn .gt. nelorbh) yesfast = 0
                end if

            else ! forceyes=.false.
                !        Define firstmol and nmolfn also in this case.
                if (yesfast .ne. 0) then
                    if (yesfast .eq. 1) then
                        firstmol = nelorb_c - molecular + 1
                    elseif (yesfast .eq. 2) then
                        if (rank .eq. 0) write (6, *) ' Warning it should be faster with yesfast= 1 !!! '
                        firstmol = 1
                    end if
                    nmolfn = lastmol - firstmol + 1
                end if
            end if ! endif forceyes

            if (yesmin .eq. 1) then
                nmoltry = 2*nelorb_c - molecular
            else
                nmoltry = nelorb_c
            end if

        end if ! endif molopt>0

        if (yesmin .ne. 0 .or. read_molecul) then

            if ((yesmin .ne. 0. .or. read_molecul) .and. rank .eq. 0) write (6, *) ' Projection scheme !!! '

            if (detc_proj .and. itestr .eq. -5 .or. read_molecul) then
                if (molecular .eq. 0) then
                    write (errmsg, *) ' ERROR molopt=-1,5 works only&
                         &   with molecular orbitals !!!'
                    call checkiflagerr(1, rank, errmsg)
                end if
                allocate (detmat_proj(ipc*nelorb_c*max(nelcol_c, nel)))
                if (rank .eq. 0) write (6, *) ' Warning molecular orbitals optimization &
                     &  with contracted coefficients '
                detmat_proj = 0.d0
            end if

            if (iopt .eq. 1 .or. iopt .eq. 3 .or. read_molecul) then
                !       Input detmat  output NEW detmat rank mol, mu_c detmat_c same Z

                ! initialize projmat or orthogonalize the orbit

                if (ndiff .ne. 0) then
                    !        check consistency molecular orbitals
                    flag = .true.
                    do i = 1, ndiff
                        do j = 1, ndiff
                            inddet = ipc*nelorb_c*(nelorb_c + i - 1) + ipc*(nelorb_c - ndiff + j) - ipc + 1
                            if (j .ne. i) then
                                if (ipc .eq. 2) then
                                    if (detmat_c(inddet) .ne. 0.d0 .or. detmat_c(inddet + 1) .ne. 0.d0) then
                                        if (ireadmin .ne. 1) then
                                            flag = .false.
                                        else
                                            detmat_c(inddet) = 0.d0
                                            detmat_c(inddet + 1) = 0.d0
                                        end if
                                    end if
                                else
                                    if (detmat_c(inddet) .ne. 0.d0) then
                                        if (ireadmin .ne. 1) then
                                            flag = .false.
                                        else
                                            detmat_c(inddet) = 0.d0
                                        end if
                                    end if
                                end if
                            else
                                if (ipc .eq. 2) then
                                    if (detmat_c(inddet) .eq. 0.d0 .and. detmat_c(inddet + 1) .eq. 0.d0) then
                                        if (ireadmin .ne. 1) then
                                            flag = .false.
                                        else
                                            detmat_c(inddet) = 1.d0
                                            detmat_c(inddet + 1) = 0.d0
                                        end if
                                    end if
                                else
                                    if (detmat_c(inddet) .eq. 0.d0) then
                                        if (ireadmin .ne. 1) then
                                            flag = .false.
                                        else
                                            detmat_c(inddet) = 1.d0
                                        end if
                                    end if
                                end if
                            end if
                        end do
                    end do

                    if (.not. flag) then
                        call checkiflagerr(1, rank, ' The unpaired molecular orbitals should be &
                             & appropriately ordered in input ')
                    end if

                end if ! endif nelup.ne.neldo

                if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                    nmollimit = molecular - ndiff
                else
                    nmollimit = (molecular - ndiff)/2
                end if

                if (nmolmatw .gt. nmollimit) then
                    write (errmsg, *) ' ERROR not enough molecular orbital in &
                         & input fort.10, nmol/nmolmax  should be less than  ', nmollimit

                    call checkiflagerr(1, rank, errmsg)
                end if

                if (ireadmin .gt. 0 .or. yesfast .eq. 0 .or. yesmin .eq. 1) then ! first big if

                    if (ireadmin .eq. 0 .or. membig .or. .not. detc_proj) then !  III big if

                        call convertmol_fast

                        !      Preparing projection matrix project_c
                        !          Test detmat_proj prepared well
                        !                       allocate(projmat_c(ipc*nelorb_at,4*nmolmat))
                        !                       projmat_c=0.d0
                        !                       detmat_c=detmat_proj

                        if (detc_proj) then
                            detmat_c = 0.d0
                            call convertmol_c
                        end if ! endif detc_proj
                    else ! referred to the III big if
                        yesmin = 0
                        !       detc_proj=.false.
                        if (detc_proj) then
                            !   Here ireadmin =/0 , meaning that I am reading the average wf.
                            detc_proj = .false.
                            call convertmol_fast
                            !  below defined slowest more general output for avoiding errors
                            nmolfn = nelorb_c
                            firstmol = 1
                            yes_ipreshuff = .false.
                        end if
                    end if ! endif III big if ireadmin
                    !               endif ! endif second big if ireadmin
                elseif (ireadmin .gt. 0) then ! referred to the first big if

                    if (rank .eq. 0) write (6, *) ' Warning this part is not tested !!! '

                    call convertmol_fast

                    !     Now is no longer needed the projection

                    yesmin = 0
                    molopt = 0

                end if ! first big if ireadmin

#if defined (_OPENMP) && defined (__NOOMP)
                call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

                !       In this way readalles makes the average of detmat_c and not mu_c
                !       mu_c is calculated below
                if (contraction .ne. 0) then

                    !            if(allocated(muc_np)) mu_c=muc_np
                    if (yes_complex) then
                        do jj = 1, iesup_c
                            do kk = 1, multranspip(jj)
                                iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                ix = transpip(kk)%col(jj) - (iy - 1)*ipf*nelorbh
                                mu_c(2*ix - 1, iy) = dup_c(2*jj - 1)
                                mu_c(2*ix, iy) = dup_c(2*jj)
                            end do
                        end do
                    else
                        do jj = 1, iesup_c
                            do kk = 1, multranspip(jj)
                                iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                                ix = transpip(kk)%col(jj) - (iy - 1)*ipf*nelorbh
                                mu_c(ix, iy) = dup_c(jj)
                            end do
                        end do
                    end if
                    !            if(allocated(muc_np)) muc_np=mu_c
                    call update_projm_

                    if (yesdetmatc) call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                         &, nelcol_c, detmat, detmat_c, mu_c, psip)

                end if

                if (contractionj .ne. 0) then

                    do jj = 1, npar3body_c
                        do kk = 1, multranspipj(jj)
                            iy = (transpipj(kk)%col(jj) - 1)/nelorbjh + 1
                            ix = transpipj(kk)%col(jj) - (iy - 1)*nelorbjh
                            muj_c(ix, iy) = vju_c(jj)
                        end do
                    end do

                    !             call scontract_mat_jas(nelorbjh,nelorbjh,nelorbjh,nelorbj_c&
                    !                  &,nelorbj_c,jasmat,jasmat_c,muj_c,psip) Ichanged
!      if(ipj.eq.2) then
!         call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!      else
!         call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh&
!              &, nelorbj_c, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!      endif

                    if (iessz) call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
                         &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)
                end if

#if defined (_OPENMP) && defined (__NOOMP)
                call omp_set_num_threads(1) ! restore the scalar
#endif

                !       updating scale used later on
                if (yesdetmat) then
                    if (yes_complex) then
                        do ii = 1, nnozero
                            scale(2*ii - 1) = detmat(2*nozero(ii) - 1)
                            scale(2*ii) = detmat(2*nozero(ii))
                        end do
                    else
                        do ii = 1, nnozero
                            scale(ii) = detmat(nozero(ii))
                        end do
                    end if
                end if
                if (contractionj .ne. 0) then
                    !       updating scalej used later on
                    do ii = 1, nnozeroj
                        scalej(ii) = jasmat(nozeroj(ii))
                        !       write(6,*) ii,scalej(ii)
                    end do
                    !       updating scalejsz used later on
                    if (iessz) then
                        do ii = 1, nnozeroj
                            scalejsz(ii) = jasmatsz(nozeroj(ii))
                        end do
                    end if
                end if

            elseif (yesmin .ne. 0) then
                !         if(.not.allocated(projmat)) then
                if (symmagp .and. ipc .eq. 1) then
                    if (detc_proj) then
                        allocate (projmat_c(nelorb_at, 2*nmolmat))
                        projmat_c = 0.d0
                    end if
                else
                    if (detc_proj) then
                        allocate (projmat_c(ipc*nelorb_at, 4*nmolmat))
                        projmat_c = 0.d0
                    end if
                end if
                !         endif ! allocated projmat
            end if
        end if ! endif molopt ne 0

        if (itest .ne. 2) then
            if (ndiff .ne. 0 .and. molecular .ne. 0) then
                nmolfn = nelorb_c - firstmol + 1
                lastmol = nelorb_c
                if (rank .eq. 0) write (6, *) ' Warning changing nmolfn for LRDMC/DMC =', nmolfn
            end if
            if (contraction .ne. 0 .and. firstmol + nmolfn - 1 .gt. nelorb_c) then
                nmolfn = nelorb_c - firstmol + 1
            end if

        end if

        if ((nelorbh .le. nmolmax .or. iessw .eq. 0) .and. rank .eq. 0 .and. yesmin .ne. 0) then
            write (6, *) ' Warning using fast algorithm, but not necessary !!! '
        end if

        if (molecular .ne. 0) then
            nmol_count = 0
            do i = firstmol, nelorb_c
                do j = firstmol, nelorb_c
                    if (detmat_c((ipc*nelorb_c)*(i - 1) + ipc*(j - 1) + 1) .ne. 0.d0) &
                         &nmol_count = nmol_count + 1
                end do
            end do
        else
            nmol_count = nelorb_c
        end if

        if (nmol_count .gt. neldo .or. ipf .eq. 2) then
            nosingledet = .true.
        else
            nosingledet = .false.
        end if

        if (rank .eq. 0) then
            write (6, *) ' Default chosen yesfast ', yesfast
            if (yesfast .ne. 0) then
                write (6, *) ' Warning nmolfn, Speeding factor =  ', nmolfn&
                     &, nelorbh/dble(2*nmolfn)
            end if
            if (itest .ne. 2 .and. ipf .ne. 2) then
                if (nosingledet) then
                    write (6, *) ' Warning AGP algorithm update '
                else
                    write (6, *) ' Warning Single determinant  algorithm update (faster) '
                end if
            end if
        end if

        if (yesfast .eq. 0 .and. contraction .ne. 0) then
            yesdetmat = .false.
            yesdetmatc = .true.
            if (allocated(detmat)) deallocate (detmat)
            allocate (detmat(ipc*ipf*nelorbh*nelcol))
            detmat = 0.d0
            if (iscramax .le. ipf*ipc*nelorbh*nelcol_c) then
                iscramax = ipf*ipc*nelorbh*nelcol_c
                deallocate (psip)
                allocate (psip(iscramax))
                psip = 0.d0
            end if

            call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                 &, nelcol_c, detmat, detmat_c, mu_c, psip)
        end if

!       For the time being detmat and nozero have to be allocated in order
!       to be compatible with readalles when yesmin>0 even when yesfast>0

        if (.not. yesdetmat .and. .not. yesdetmatc) then

            !         in some cases it is also not necessary detmat and nozero

            deallocate (detmat, nozero)
            if (yes_complex) then
                allocate (detmat(2*nelorbh), nozero(1))
            else
                allocate (detmat(nelorbh), nozero(1))
            end if
            detmat = 0.d0
            nozero = 0
            if (itestr .ne. -5) then
                deallocate (nozerodet)
                allocate (nozerodet(1))
                nozerodet = 0
            end if

        end if

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

        if (rank .eq. 0) then
            write (6, *) '*********************************************'
            write (6, *) '*********************************************'
            if (itestr .eq. -5) then

                if (itestr4 .ge. -10) then
                    write (6, *) 'VMC with ENERGY MINIMIZATION'
                else
                    write (6, *) 'FN   with ENERGY MINIMIZATION'
                end if

                if (itestrr .eq. -4) then
                    write (6, *) 'STOCHASTIC RECONFIGURATION WITH HESSIAN'
                    write (6, *) ' Beta used =', beta
                elseif (itestrr .eq. -5) then
                    write (6, *) 'SIMPLE STOCHASTIC RECONFIGURATION'
                end if
                write (6, *) ' Stochastic acceleration =', tpar
                if (itestr3 .eq. -9) write (6, *)                                   &
                     &' No optimization Z with contracted '
            elseif (itestr .eq. 2) then
                write (6, *) 'VARIATIONAL MONTE CARLO'
            elseif (itestr .eq. 1 .or. itestr .eq. -2 .or. itestr .eq. -3) then
                write (6, *) 'DIFFUSION MONTE CARLO'
                if (rejweight) then
                    write (6, *) 'weights updated according to acceptance/rejection in diffusion'
                else
                    write (6, *) 'weights always updated'
                end if
            elseif (itestr .eq. -6) then
                write (6, *) 'LATTICE REGULARIZED DIFFUSION MONTE CARLO'
            end if
            write (6, *) '*********************************************'
            write (6, *) '*********************************************'
        end if

        if (itestrr .eq. -4 .or. itestrr .eq. -5) then

            if (contraction .ne. 0 .or. contractionj .ne. 0) then
                if (iesup .ne. 0 .and. .not. yeszagp) then
                    if (rank .eq. 0) write (6, *) 'zeta exponents in determinant excluded'
                end if
                if (iesm .ne. 0 .and. .not. yeszj) then
                    if (rank .eq. 0) write (6, *) 'zeta exponents in Jastrow  excluded'
                end if
                if ((iesm .ne. 0 .and. contractionj .eq. 0)&
                     &.and. (.not. yeszj .or. itestr3 .eq. -4 .or. itestr3 .eq. -9)) then
                    if (rank .eq. 0) write (6, *) ' Warning optimization Z Jastrow on '
                    yeszj = .true.
                end if
                if ((iesup .ne. 0 .and. contraction .eq. 0)                           &
                     &.and. (.not. yeszagp .or. itestr3 .eq. -4 .or. itestr3 .eq. -9)) then
                    write (6, *) ' Warning optimization Z AGP  on '
                    yeszagp = .true.
                end if
            end if

        end if ! itestrr.eq.-4.or.itestrr.eq.-5

        inddsw = iesfree + iesinv + iesm + iesd + 1
        stodim = iesm + iesd + iesup_read + ipc*nnozero_c + iesinv + iesfree + 3*nion + 3

        stodim = max(stodim, nmat)

! allocation cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        allocate (dup(max(iesup, 1)))
        dup = 0.d0
        if (itestr .eq. -5) then ! only for the optimization this memory is required
            if (stodim .gt. iscramax) then
                if (rank .eq. 0) then
                    write (6, *) ' Warning changing iscramax for the master ', stodim
                    deallocate (psip)
                    iscramax = stodim
                    allocate (psip(iscramax))
                    psip = 0.d0
                end if
            end if
            !      ALLOCATE(alphasto(stodim))
            allocate (ef(ieskindim, 3, nbindim), efenergy(1 + ipc, nbindim)             &
                 &, err(npm), force(npm), efp(ndimpdim, 2, nbindim))
            if (yes_adams) then
                allocate (first_moment(npm), second_moment(npm))
            else
                allocate (first_moment(1), second_moment(1))
            end if
            first_moment = 0.d0
            second_moment = 0.d0
            ef = 0.d0
            efenergy = 0.d0
            err = 0.d0
            force = 0.d0
            efp = 0.d0
            allocate (fk(nbindim, npdim), fkav(npm), reduce(ncgdim, npdim)&
                 &, okav(npm), skdiag(npm))
            allocate (velion(3, ieskindim), efpress(3, nbindim))
            efpress = 0.d0
            velion = 0.d0
            !     Initialize randomly according to the given temperature
            if (temp .ne. 0.d0) then
                if (yesquantum) then
                    if (rank .eq. 0) write (6, *) ' Initializing velocities at temp (a.u.)', nbead*temp*ris(2)
                    do ii = 1, ieskindim
                        arg1 = 1.d0 - drand1()
                        arg2 = TWO_PI*drand1()
                        velion(3, ii) = dsqrt(-2.d0*dlog(arg1)*nbead*temp)*dcos(arg2)
                    end do
                else
                    if (rank .eq. 0) write (6, *) ' Initializing velocities at temp (a.u.)', temp*ris(2)
                    do ii = 1, ieskindim
                        arg1 = 1.d0 - drand1()
                        arg2 = TWO_PI*drand1()
                        velion(3, ii) = dsqrt(-2.d0*dlog(arg1)*temp)*dcos(arg2)
                    end do
                end if
            end if
            fk = 0.d0
            fkav = 0.d0
            okav = 0.d0
            skdiag = 0.d0
            reduce = 0.d0
        end if !endif itestr=-5

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        stepcg = 0

        deallocate (scale)
        allocate (scale(1))

#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

        if (allocated(work)) deallocate (work)

        if (yesfast .eq. 0 .and. rank .eq. 0) then ! main test diagonalization
            deallocate (psip) ! to save memory
            ! calculation eigenvalues
            if (contraction .eq. 0) then
                if (.not. yes_complex) then
                    allocate (work(ipf*ipf*nelorbh*nelorbh + 3*ipf*nelorbh))
                    call dcopy(ipf*ipf*nelorbh*nelorbh, detmat, 1, work, 1)
                    allocate (eig(2*ipf*nelorbh))
                    eig = 0.d0
                    if (symmagp .and. ipf .eq. 1) then
                        call dsyev('N', 'U', nelorbh, work, nelorbh, eig&
                             &, work(nelorbh*nelorbh + 1), 3*nelorbh, info)
                    else
                        call dgeev('N', 'N', ipf*nelorbh, work, ipf*nelorbh&
                             &, eig, eig(ipf*nelorbh + 1), vl, 1, vr, 1, work(ipf*ipf*nelorbh*nelorbh + 1)&
                             &, 3*ipf*nelorbh, info)
                        do i = 1, ipf*nelorbh
                            eig(i) = sqrt(eig(i)**2 + eig(ipf*nelorbh + i)**2)
                        end do
                    end if
                else ! complex eigenvalues
                    allocate (work(2*ipf*ipf*nelorbh*nelorbh + 8*ipf*nelorbh))
                    call dcopy(2*ipf*ipf*nelorbh*nelorbh, detmat, 1, work, 1)
                    allocate (eig(2*ipf*nelorbh))
                    call zgeev('N', 'N', ipf*nelorbh, work, ipf*nelorbh, eig, vl, 1, vr, 1&
                         &, work(2*ipf*ipf*nelorbh*nelorbh + 1), 3*ipf*nelorbh&
                         &, work(2*ipf*ipf*nelorbh*nelorbh + 6*ipf*nelorbh + 1), info)
                    do i = 1, 2*ipf*nelorbh, 2
                        eig(i) = sqrt(eig(i)**2 + eig(i + 1)**2)
                    end do
                end if
            else ! contracted orbitals
                if (.not. yes_complex) then
                    allocate (work(nelorb_c*nelorb_c + 5*nelorb_c))
                    call dcopy(nelorb_c*nelorb_c, detmat_c, 1, work, 1)
                    if (symmagp) then
                        allocate (eig(nelorb_c))
                        eig = 0.d0
                        call dsyev('N', 'U', nelorb_c, work, nelorb_c, eig &
                             &, work(nelorb_c*nelorb_c + 1), 3*nelorb_c, info)
                    else
                        allocate (eig(nelorb_c))
                        eig = 0.d0
                        call dgesvd('N', 'N', nelorb_c, nelorb_c, work, nelorb_c, eig, vl, 1, vr, 1&
                             &, work(nelorb_c*nelorb_c + 1), 5*nelorb_c, info)
                    end if
                else ! complex eigenvalues
                    allocate (work(2*nelorb_c*nelorb_c + 11*nelorb_c))
                    call dcopy(2*nelorb_c*nelorb_c, detmat_c, 1, work, 1)
                    allocate (eig(nelorb_c))
                    eig = 0.d0
                    call zgesvd('N', 'N', nelorb_c, nelorb_c, work, nelorb_c, eig, vl, 1, vr, 1&
                         &, work(2*nelorb_c*nelorb_c + 1), 3*nelorb_c, work(2*nelorb_c*nelorb_c + 6*nelorb_c + 1), info)
                end if
            end if

#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(1) ! restore the scalar code
#endif

            write (6, *) '%%%%%%%%%%%%%%% DETERMINANTAL GEMINAL %%%%%%%%%%%%%%'
            write (6, *) ' Eigenvalues matrix lambda ', info, rank

            if (contraction .eq. 0) then
                maxdimeig = nelorbh*ipf
            else
                maxdimeig = nelorb_c
            end if

            call print_eigenvalues(rank, eig, maxdimeig, irankdet)

            if (allocated(work)) deallocate (work)
            deallocate (eig)
            allocate (psip(iscramax))
            psip = 0.d0

        end if ! endif test diagonalization

#ifdef PARALLEL
        call mpi_bcast(irankdet, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

        if (yesfast .eq. 0) then ! main if test diagonalization

            if (irankdet .lt. neldo) then
                write (errmsg, *) ' Singular matrix !!! ', irankdet
                call checkiflagerr(1, rank, errmsg)
            elseif (irankdet .eq. neldo) then
                if (rank .eq. 0) write (6, *) ' No resonance simple det !!! ', irankdet
            else
                if (rank .eq. 0) write (6, *) ' RVB is starting !!! ', irankdet
            end if

        end if ! if yesfast=0

!     copy scalej in jasmat

        if (contraction .eq. 0 .and. .not. yes_sparse) then
            call dscalzero(ipj*ipj*nelorbjh*nelorbjh, 0.d0, jasmat, 1)
            do i = 1, nnozeroj
                call upsim(jasmat, ipj*nelorbjh, nozeroj(i), scalej(i), .true., 1)
            end do
        end if
        if (developer .eq. 0) then
            deallocate (scalej)
            allocate (scalej(1))
            scalej = 0.d0
        end if

        if (iessz) then
            call dscalzero(nelorbjh*nelorbjh, 0.d0, jasmatsz, 1)
            do i = 1, nnozeroj
                call upsim(jasmatsz, nelorbjh, nozeroj(i), scalejsz(i), .true., 1)
            end do
            if (developer .eq. 0) then
                deallocate (scalejsz)
                allocate (scalejsz(1))
                scalejsz = 0.d0
            end if
        end if

        if (contractionj .ne. 0) then
            call constrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c, jasmat_c   &
                 &, ddw, 1, 1)
            if (iessz)                                                     &
                 &    call constrbra(iesinv, nnozeroj_c, jbrajsz, nozeroj_c            &
                 &, jasmatsz_c, ddwsz, 1, 1)
        else
            if (yes_sparse) then
                call constrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder, jasmat&
                     &, ddw, 1, 1)
            else
                call constrbra(iesfree, nnozeroj, jbraj, nozeroj, jasmat         &
                     &, ddw, 1, 1)
            end if
            if (iessz)                                                     &
                 &    call constrbra(iesinv, nnozeroj, jbrajsz, nozeroj                &
                 &, jasmatsz, ddwsz, 1, 1)
        end if

        call constrbr_complex(iesup, iesup_c, jbraiesup, dup_c                  &
           &, dup, 1, 1)

        call constrbr(iesm, npar3body_c, jbraiesm, vju_c                &
           &, vju, 1, 1)

        if (rank .eq. 0 .and. allowed_averagek) then
            write (6, *) ' Warning allowed translation via phase flux attaching'
        end if

        if (contraction .ne. 0) then
            !  From real to effective load dsw
            if (rank .eq. 0) write (6, *) ' Passi qui I real-eff '
            if (allowed_averagek) call attach_phase2det(.false., detmat_c)
            call constrbra_complex(iessw, nnozero_c, jbradet, nozero_c, detmat_c&
                 &, dsw, 1, 1)
            !  back to real
            if (allowed_averagek) call attach_phase2det(.true., detmat_c)
            if (rank .eq. 0) write (6, *) ' Passi qui III eff-real '
        else
            !  From real to effective load dsw
            if (rank .eq. 0) write (6, *) ' Passi qui II real-eff '
            if (allowed_averagek) call attach_phase2det(.false., detmat)
            call constrbra_complex(iessw, nnozero, jbradet, nozero, detmat            &
                 &, dsw, 1, 1)
            !  back to real
            if (rank .eq. 0) write (6, *) ' Passi qui IV eff-real '
            if (allowed_averagek) call attach_phase2det(.true., detmat)
        end if

        if (rank .eq. 0) then
            write (6, *) '%%%%%%%%%%%%%%%%%% 2 BODY JASTROW %%%%%%%%%%%%%%%%'
            write (6, *) 'Number 2 body Jastrow parameters ', niesd
            if (niesd .ge. 1) then
                do i = 1, niesd
                    write (6, *) i, vj(i)
                end do
            else
                write (6, *) '2 body Jastrow not present'
            end if
            if (iesdr .le. -5) then
                write (6, *) ' initial costz, costz3, zeta_Q_Caffarel '
                do i = 1, nion
                    write (6, *) i, costz(i), costz3(i), zetaq(i)
                end do
            end if
            if (iesfree .ne. 0 .or. iesm .ne. 0 .or. iesfreesz .ne. 0) then
                write (6, *) '%%%%%%%%%%%%%%%%%% 3 BODY JASTROW %%%%%%%%%%%%%%%%'
            end if
            write (6, *) 'Total independent parameters in Jastrow sector', &
                 &              inddsw - 1
            write (6, *)

        end if ! endif rank=0

        iesupskip = iesupr + npar3body
        iskip = iesupskip + 1

!----------------------------------------------
! Definition of main indices related to energy,
! corr functions and variational paramers
!----------------------------------------------

! first the energy index
        if (iese .ne. 0) then
            !if(.not.yes_complex) then
            nwinv = nwinv + in1*iese
            nwm2 = nwm2 + in1*iese
            nwdim = nwdim + in1*iese
            nwfix = nwfix + in1*iese
            nwnel = nwnel + in1*iese
            nwrep = nwrep + in1*iese
            nwdw = nwdw + in1*iese
            nwfree = nwfree + iese*in1
            nwsw = nwsw + iese*in1
            nwkin = nwkin + iese*in1
            nwup = nwup + iese*in1
            nwking = nwking + iese*in1
            nprest = nprest - iese
            nindt = nindt + iese
        end if

! then indices related to the w.f. parameters
        if (iesinv .ne. 0) then
            nwm2 = nwm2 + in1*iesinv
            nwdim = nwdim + in1*iesinv
            nwrep = nwrep + in1*iesinv
            nwfix = nwfix + in1*iesinv
            nwnel = nwnel + in1*iesinv
            nwdw = nwdw + in1*iesinv
            nwfree = nwfree + iesinv*in1
            nwsw = nwsw + iesinv*in1
            nwkin = nwkin + iesinv*in1
            nwup = nwup + iesinv*in1
            nwking = nwking + iesinv*in1
            nprest = nprest - iesinv
            nindt = nindt + iesinv
        end if

        if (iesm .ne. 0) then
            nwdim = nwdim + abs(iesm)*in1
            nwrep = nwrep + in1*abs(iesm)
            nwfix = nwfix + abs(iesm)*in1
            nwnel = nwnel + in1*abs(iesm)
            nwdw = nwdw + in1*abs(iesm)
            nwfree = nwfree + abs(iesm)*in1
            nwsw = nwsw + abs(iesm)*in1
            nwkin = nwkin + abs(iesm)*in1
            nwup = nwup + abs(iesm)*in1
            nwking = nwking + abs(iesm)*in1
            nprest = nprest - abs(iesm)
            nindt = nindt + abs(iesm)
        end if

        if (iesd .ne. 0) then
            nwrep = nwrep + in1*abs(iesd)
            nwfix = nwfix + abs(iesd)*in1
            nwnel = nwnel + in1*abs(iesd)
            nwdw = nwdw + in1*abs(iesd)
            nwfree = nwfree + abs(iesd)*in1
            nwsw = nwsw + abs(iesd)*in1
            nwkin = nwkin + abs(iesd)*in1
            nwup = nwup + abs(iesd)*in1
            nwking = nwking + abs(iesd)*in1
            nprest = nprest - abs(iesd)
            nindt = nindt + abs(iesd)
        end if

        if (iesfree .ne. 0) then
            nwrep = nwrep + in1*iesfree
            nwnel = nwnel + in1*iesfree
            nwfix = nwfix + iesfree*in1
            nwkin = nwkin + iesfree*in1
            nwup = nwup + iesfree*in1
            nwking = nwking + iesfree*in1
            nwsw = nwsw + iesfree*in1
            nwdw = nwdw + in1*iesfree
            nprest = nprest - iesfree
            nindt = nindt + iesfree
        end if

        if (iessw .ne. 0) then
            nwrep = nwrep + in1*iessw
            nwnel = nwnel + in1*iessw
            nwkin = nwkin + iessw*in1
            nwfix = nwfix + iessw*in1
            nwup = nwup + iessw*in1
            nwking = nwking + iessw*in1
            nwdw = nwdw + in1*iessw
            nprest = nprest - iessw
            nindt = nindt + iessw
        end if

        if (iesup .ne. 0) then
            nwkin = nwkin + iesup*in1
            nwking = nwking + iesup*in1
            nwrep = nwrep + in1*iesup
            nwfix = nwfix + iesup*in1
            nwnel = nwnel + in1*iesup
            nwdw = nwdw + in1*iesup
            nprest = nprest - iesup
            nindt = nindt + iesup
        end if

        if (iesking .ne. 0) then
            nwkin = nwkin + iesking*in1
            nwrep = nwrep + in1*iesking
            nwfix = nwfix + iesking*in1
            nwnel = nwnel + in1*iesking
            nwdw = nwdw + in1*iesking
            nprest = nprest - iesking
            nindt = nindt + iesking
        end if

        if (ieskin .ne. 0) then
            nwrep = nwrep + in1*ieskin*skipforce
            nwfix = nwfix + ieskin*in1*skipforce
            nwnel = nwnel + in1*ieskin*skipforce
            nwdw = nwdw + in1*ieskin*skipforce
            nprest = nprest - ieskin*skipforce
            nindt = nindt + ieskin*skipforce
        end if

        indfix = nindt + 1
        if (isfix .ne. 0) then
            nwrep = nwrep + in1*isfix
            nprest = nprest - isfix
            nindt = nindt + isfix
        end if
        indberry = nindt
        repf = nwrep/in1 + 1
        kinf = nwkin/in1 + 1

        if (itestr .ne. -5) then ! if the minimization is not assumed
            nrep = iesm + iesd + iesfree + iessw + iesinv + iesup + iesking
            nprest = nprest - nrep
            nindt = nindt + nrep
        else
            nrep = 0
        end if

!---------------------------
! End of indices definitions
!---------------------------

        if (rank .eq. 0) then
            write (6, *) ' nwfix given =', nwfix, nwrep
            write (6, *) ' Number of parameters in SR =', nindt
            write (6, *) ' Hartree Atomic Units '
        end if

        if (nprest .lt. 0) then
            write (errmsg, *) ' ERROR INCREASE NP !!! by ', -nprest, '=', np - nprest
            call checkiflagerr(1, rank, errmsg)
        elseif (nprest .gt. 0) then
            write (errmsg, *) ' ERROR DECREASE NP !!! by ', nprest, '=', np - nprest
            call checkiflagerr(1, rank, errmsg)
        end if

!---------------------------------------------
! Definition of dimension andaddresses of the
! main matrices
!---------------------------------------------

        nelnw = nel*nw*(indt + 1)
        nelnw0 = nel*nw
        nrnel = (indt + 1)*nel
        nelkel = 3*nrnel
        dimee = nel*nel*(indt4j + 1)
        dimei = nel*nion
        nrest = 0

        Lzeff = nel*indteff
        Lztab = nel*indt
        if (Lztab .eq. 0) Lztab = 1
        Lztabr = Lztab
        if (yes_complex) Lztab = Lztab*2
        Ltab = nel*(indt + ip4)
        Ltabb = nel*max(indt, 1)

        nwnp = in1*nmat

! save the starting ionic configurations
        call dscalzero(ieskindim, 0.d0, dek, 1)
        call dscalzero(ieskingdim, 0.d0, dekg, 1)

        if (iespbc .and. itestr .eq. -5) then
            if (ieskint .eq. ieskinr_pos + 2) then
                dek(ieskinr_pos + 1 - iesking) = cellscale(2)
                dek(ieskinr_pos + 2 - iesking) = cellscale(3)
                if (rank .eq. 0) write (6, *) ' PBC Updating cellscale b,c=', ieskinr_pos &
                     &, dek(ieskinr_pos + 1 - iesking), dek(ieskinr_pos + 2 - iesking)
            elseif (ieskint .eq. ieskinr_pos + 1) then
                dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                if (rank .eq. 0) write (6, *) ' PBC Updating cellscale a=', ieskinr_pos   &
                     &, dek(ieskinr_pos + 1 - iesking)
            elseif (ieskint .eq. ieskinr_pos + 3) then
                dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                dek(ieskinr_pos + 2 - iesking) = cellscale(2)
                dek(ieskinr_pos + 3 - iesking) = cellscale(3)
                if (rank .eq. 0) write (6, *) ' PBC Updating cellscale a,b,c=', ieskinr_pos  &
                     &, dek(ieskinr_pos + 1 - iesking), dek(ieskinr_pos + 2 - iesking), dek(ieskinr_pos + 3 - iesking)
            end if
        end if

        do j = 1, nw
            jbra(j) = j
        end do

! doubled indices for complex allocation
! only determinantal part
        if (yes_complex) then
            nel2wt = 2*nelorb*nel*(indt4 + 1) ! address for winv
            nel2up = 2*nelup_mat*nelup_mat ! address for ainv
            nel2upt = 2*nelup*(indt + ip4) ! address for winvup
            nel2dot = 2*neldo*(indt + ip4) ! address for winvdo
            nel2bar = 2*ipf*nel_mat*nelorbh ! address for winvbar
        else
            nel2wt = nelorb*nel*(indt4 + 1)
            nel2up = nelup_mat*nelup_mat
            nel2upt = nelup*(indt + ip4)
            nel2dot = neldo*(indt + ip4)
            nel2bar = ipf*nel_mat*nelorbh
        end if

! indices for the jastrow
        nel2wtj = nelorbj*nel*(indt4j + 1)
        nel2jbar = (nelup + neldo)*nelorbjh*ipj
        if (iessz) then
            nel2jbarsz = nel2jbar
        else
            nel2jbarsz = 1
        end if

        nelnion = (nelup + neldo)*nion

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccc parallelization: each process goes from
!cccccc ist to ien walkers
        ist = rank*(nw/nproc) + 1
        istm = ist - 1
        ien = (rank + 1)*(nw/nproc)
        id1 = rank*(nw/nproc)
#ifdef PARALLEL
        if (yes_fastbranch) then
            skipreshuff = 21*nel
        else
            skipreshuff = 14 + ipc + Ltab + Ltabb + Lztab + Lztabr + nel2upt + nel2dot + 8*nel + nelkel + npf
            skipreshuff = skipreshuff + 18*nel

            if (yesivic) skipreshuff = skipreshuff + 3*indt*nel
            if (npsa .gt. 0) skipreshuff = skipreshuff + nel
        end if
        skip = skipreshuff*in1
        if (iscramax .lt. skip .and. iesbra) then
            iscramax = skip
            deallocate (psip)
            allocate (psip(skip))
            psip = 0.d0
            if (rank .eq. 0) write (6, *) 'iscramax changed!', skip
        end if
#endif
        if (iscramax .lt. nparshellmax + nelorbpp) then
            iscramax = nparshellmax + nelorbpp
            deallocate (psip)
            allocate (psip(iscramax))
            psip = 0.d0
        end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Initialization
        table = 0.d0
        tabler = 0.d0

        if (itestr .eq. -2) then
            ! heat bath after all electron swap
            ! identity splitted into nel particle contributions
            identity = 1.d0/dble(nel)/tbra
        elseif (itestr .eq. -3) then
            ! heat bath after single particle diffusion
            ! non local operator splitted into a product of nel factors < e^{-tbra/nel V_nl} >
            ! tbra --> tbra/nel
            identity = 1.d0/tbra
        end if

        if (itestr .eq. -2 .or. itestr .eq. -3) then
            ! fill table and tabler with the identity in the last entry per particle
            ! this value will never be touched during the simulation
            do kk = 1, nws
                indtju = nel*indt*(kk - 1)
                if (yes_complex) then
                    do k = 1, 2*nel
                        table(k + 2*(indt - 1)*nel + 2*indtju) = dcmplx(identity)
                    end do
                    do k = 1, nel
                        tabler(k + (indt - 1)*nel + indtju) = identity
                    end do
                else
                    do k = 1, nel
                        table(k + (indt - 1)*nel + indtju) = identity
                        tabler(k + (indt - 1)*nel + indtju) = identity
                    end do
                end if
            end do
        end if

        bcost = 1.d0
        ccost = 1.d0
!     Initialize to zero the global variables
        kel = 0.d0
        wconfn = 0.d0

        if (.not. iesrandoma .and. itest .ne. 2 .and. alat .ne. 0.d0) then

            if (.not. allocated(rion_sav)) then
                allocate (rion_sav(3, nion))
                yes_deallocate = .true.
            else
                yes_deallocate = .false.
            end if
            rion_sav = rion
            if (iespbc) then
                call CartesianToCrystal(rion_sav, nion)
                call CartesianToCrystal(rion_ref, 1)
            end if

            call findrionfref(nion, alat, rion_sav, rion_ref, mind(1))
            call findrionfref(nion, alat, rion_sav(2, 1), rion_ref(2), mind(2))
            call findrionfref(nion, alat, rion_sav(3, 1), rion_ref(3), mind(3))
            if (rank .eq. 0) write (6, *) ' Minimum ion-mesh distance =', dsqrt(sum(mind(:)**2))
            if (iespbc) rion_ref(:) = rion_ref(1)*at(:, 1) + rion_ref(2)*at(:, 2) + rion_ref(3)*at(:, 3)
            if (yes_deallocate) then
                deallocate (rion_sav)
            else
                rion_sav = rion
            end if

        end if

        call initconf(nel, nelup, psiln, psisn, kel, wconfn, rion, dist&
           &, indt, nion, zetar, in1, rank, ierr, LBox, alat, rion_ref, iesrandoma, itest)

        imin = 0
        itry = 0

        psisn = 0.d0
        singdet = .true.
        psiln = 0.d0
        countscra = 0
        spsi = 1.d0
        vcut = 0.d0
        yescut = .false.
        flagmu = .false.
        wconfsav(:) = 1.d0

        ntry = 0
        nmatb = 0
        ttry = 0.d0
        ngentry = 0
        rata = 0.d0
        signflip = 0.d0
        nontr = 0.d0
        naccpseudo = 0.d0
        tmes = temp*ris(2)
        epscutu = epscut
        epstlu = epstl
        ngs = 0
        nmp = 0
        ng = 0
        icdiff = 0
        wbra = 1.d0
        varmin = 0.d0
        enermin = 0.d0
        acclarge = 0.d0
        indvic = 1
        parbest = 0.d0
        iflagnorm = 3
! iflagnorm=3 compute the normalization of orbitals and distances
! iflagnorm=2 compute only distances (rmu = el-nu , r=|el-nu|)
        coeff = 0.d0
        jmax = 0
        fmax = 0.d0

        if (mod(nbra, nel) .ne. 0 .and. fncont) then
            nbra = ((nbra - 1)/nel + 1)*nel
            if (rank .eq. 0)                                                &
                 &     write (6, *) ' Warning changing nbra multiple of nel =', nbra
        end if
        nbram = nbra - 1

        call dscalzero(nelorbh*(indt + 1 + ip4), 0.d0, psinew, 1)
        psip = 0.d0

        do j = 1, nwnp
            econf(j) = 1.d0
            econfh(j) = 0.d0
        end do
        do j = 1, in1
            econfh(in1*(npp - 1) + j) = 1.d0
        end do

        iend = 0
        costpassed = 1.d0

        do i = 1, nfat
            wcorw(i) = 1.d0
        end do
        indcor = nfat
        wcort = 1.d0

        psiav = 0.d0
        psisav = 0.d0
        psirav = 0.d0
        psirav_all = 0.d0
        countav = 0.d0
        countt = 0.d0
        avenernum = 0.d0
        avenerden = 0.d0
        tave_cyrus = 0.d0
        tcount_cyrus = 0.d0
        counttot = 0.d0
        countcut = 0.d0
        countreg = 0.d0

        if (idyn .ge. 2 .and. idyn .ne. 5) then
            ldynsecond = .true.
        else
            ldynsecond = .false.
        end if

        if (ldynsecond) then
            call dscal(ieskin, ris(5), scalpar(np - ieskin + 1), 1)
        elseif (idyn .ne. 5) then
            call dscal(ieskin, ris(2), scalpar(np - ieskin + 1), 1)
        elseif (.not. yesquantum) then
            call dscal(ieskin, 2.d0, scalpar(np - ieskin + 1), 1)
        end if

        dt = 0.d0
        ii = 1
!           The first one non zero
        do while (dt .eq. 0.d0 .and. ii .le. ieskin)
            dt = scalpar(np - ieskin + ii)
            ii = ii + 1
        end do
        ref_comp = ii - 1
        ind = 0
        ref_atom = 0
        do ii = 1, ieskinr
            if (ion_table(ii)%mult .gt. 0 .and. atom_number(abs(ion_table(ii)%ion(1))) .gt. 0) then
                ind = ind + 1
                if (ind .eq. ref_comp) ref_atom = ion_table(ii)%ion(1)
            end if
        end do

        if (rank .eq. 0 .and. ldynsecond .and. idyn .gt. 0) write (6, *)&
            & ' Warning Chosen Dt /component (a.u.) =', dt/ris(5), ref_comp
        if (rank .eq. 0 .and. .not. ldynsecond .and. idyn .gt. 0) write (6, *)&
           & ' Warning Chosen Dt /component (a.u.) =', dt/ris(2), ref_comp
        if (dt .eq. 0 .and. idyn .gt. 0) then
            if (rank .eq. 0) write (6, *) ' ERROR no component to do dynamics  Dt=0!!!!', dt
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

!   Inizialize and change if reading
        nacc = 0.d0
        nmovet = 0.d0
        yeschange_nweight = .false.
        if (nbra_cyrus .gt. 0) then
            dim_cyrus = 3*(2*nion + 2)*(nbra_cyrus + 2)
            allocate (first_cyrus(in1), queue_cyrus(3, 2*nion + 2, nbra_cyrus + 2, in1))
            queue_cyrus = 0.d0
            first_cyrus = 1
        else
            dim_cyrus = 0
        end if
        nbra_cyrus_read = nbra_cyrus

! iopt=0   continue the run starting from a given conf
! iopt=3   contine  the run but restart projection matrices with molopt=-1
! iopt=1   begin the calculation from scratch
! iopt=2   begin with old conf but restart averages
        if (iopt .eq. 0 .or. iopt .eq. 2 .or. iopt .eq. 3) then

            if (rank .eq. 0) call read_fort11_begin
#ifdef PARALLEL
            call mpi_bcast(ngenc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nwr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nbrar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npow, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(etryr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(tstep, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nrest, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nwr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(tbrar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npbrar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ireadr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nmpr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nfatr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(np3r, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npsovr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(wcort, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nweightr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ngs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(iendr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(indcor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(parbest, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                 &, ierr)
            call mpi_bcast(iflagerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(enermin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                 &, ierr)
            call mpi_bcast(varmin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
                 &, ierr)
            call mpi_bcast(fmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(jmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(inext, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(pippoc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(pippo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nprocr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nbinread, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(stepcg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ncgread, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(epscutur, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD   &
                 &, ierr)
            call mpi_bcast(epstlur, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                 &, ierr)
            call mpi_bcast(tmes, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD       &
                 &, ierr)
            if (scalermax) then
                call mpi_bcast(rmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                     &, ierr)
                call mpi_bcast(rmaxj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                     &, ierr)
            end if
            call mpi_bcast(iskipdynr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(epscuttyper, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(parr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD   &
                 &, ierr)
            call mpi_bcast(iesconv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(iesdelay, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(delay_changeparr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(tpar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD       &
                 &, ierr)
            call mpi_bcast(nbra_cyrus_read, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
#endif

            if (nbra_cyrus .ne. nbra_cyrus_read) then
                if (rank .eq. 0) write (6, *) ' ERROR nbra_cyrus cannot be changed! use the same nbra_cyrus also in VMC!!!'
#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            end if

            if ((iopt .eq. 0 .or. iopt .eq. 3) .and. (nprocr .eq. nproc)) then
                if (io_level .eq. 1) iflagerr = read_seeds(cstr(trim(ranseedfilename)))
#ifdef PARALLEL
                if (io_level .eq. 2) iflagerr = read_seeds_mpiio(cstr(trim(scratchpath)//'randseed.all'), rank)
#endif
                call checkiflagerr(iflagerr, rank, "Fail in reading random seeds!")
            elseif (iopt .eq. 0 .or. iopt .eq. 3) then
                if (rank .eq. 0) then
                    write (6, *) ' Warning continuing with different number of processors '
                    write (6, *) ' Warning the random number are initialized again '
                end if
            end if

#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(1) ! restore the scalar code
#endif

            if (itestr .eq. -5 .and. writescratch .eq. 0) rewind (19)
            nmatr = npr + 1
            if (iopt .eq. 0 .or. iopt .eq. 3) then
                !         nfat=nfatr
                !         npsov=npsovr
                if (nweight .eq. nweightr .and. iskipdyn .eq. iskipdynr) then
                    iend = iendr
                else
                    if (rank .eq. 0) write (6, *) ' Warning starting again the record counter '
                    iend = 0
                    yeschange_nweight = .true.
                end if

                !         np3r=np3
                !         np=npr
                !         npbra=npbrar
                ! fine iopt.eq.0
            end if
            !---------------------------------------------------

            if (rank .eq. 0) call read_fort11_end
            call checkiflagerr(iflagerr, rank, 'Error in reading fort.11')

            if (nwr .eq. nw) then
#ifdef PARALLEL
                if (yesquantum) then
                    if (rankrep .eq. 0) then
                        if (io_level .eq. 1) then
                            rewind (8)
                            if (.not. yesread10) then
                                if (contraction .ne. 0) then
                                    read (8) rion, velion, alphavar, detmat_c, projm, mu_c&
                                         &, jasmat, jasmatsz, jasmat_c, jasmatsz_c, vjur, dupr
                                else
                                    read (8) rion, velion, alphavar, detmat, projm, mu_c&
                                         &, jasmat, jasmatsz, jasmat_c, jasmatsz_c, vjur, dupr
                                end if
                                if (allocated(muj_c)) read (8) muj_c
                                if (npar_eagp .gt. 0) read (8) eagp_pfaff
                            else
                                read (8) rion, velion
                            end if
                            if (iopt .ne. 3) then
                                if (detc_proj) then
                                    read (8) detmat_proj, projmat_c
                                end if
                                if (itestr .eq. -5) read (8) reduce
                            end if
                        elseif (io_level .eq. 2) then
                            if (.not. yesread10) then
                                call MPI_File_read_all(quantcont%fp, rion, size(rion), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, velion, size(velion), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, alphavar, size(alphavar), MPI_DOUBLE_PRECISION, status, ierr)
                                if (contraction .ne. 0) then
                                    call MPI_File_read_all(quantcont%fp&
                                         &, detmat_c, size(detmat_c)&
                                         &, MPI_DOUBLE_PRECISION, status, ierr)
                                else
                                    call MPI_File_read_all(quantcont%fp, detmat, size(detmat), MPI_DOUBLE_PRECISION, status, ierr)
                                end if
                                call MPI_File_read_all(quantcont%fp, projm, size(projm), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, mu_c, size(mu_c), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, jasmat, size(jasmat), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, jasmat_c, size(jasmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, jasmatsz, size(jasmatsz), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, jasmatsz_c&
                                    &, size(jasmatsz_c), MPI_DOUBLE_PRECISION&
                                    &, status, ierr)
                                call MPI_File_read_all(quantcont%fp, vjur, size(vjur), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, dupr, size(dupr), MPI_DOUBLE_PRECISION, status, ierr)
                                if (allocated(muj_c))& !  read(8) muj_c
                                     &call MPI_File_read_all(quantcont%fp, muj_c, size(muj_c), MPI_DOUBLE_PRECISION, status, ierr)
                                if (npar_eagp .gt. 0)& !  read(8) muj_c
                             &call MPI_File_read_all(quantcont%fp, eagp_pfaff&
                             &, size(eagp_pfaff), MPI_DOUBLE_PRECISION, status, ierr)

                            else
                                call MPI_File_read_all(quantcont%fp, rion, size(rion), MPI_DOUBLE_PRECISION, status, ierr)
                                call MPI_File_read_all(quantcont%fp, velion, size(velion), MPI_DOUBLE_PRECISION, status, ierr)
                            end if ! endif .not.yeswrite10
                            if (iopt .ne. 3) then
                                if (detc_proj) then ! read(8) detmat_proj,projmat_c
                                    call MPI_File_read_all(quantcont%fp&
                                        &, detmat_proj, size(detmat_proj)&
                                        &, MPI_DOUBLE_PRECISION, status, ierr)
                                    call MPI_File_read_all(quantcont%fp&
                                        &, projmat_c, size(projmat_c)&
                                        &, MPI_DOUBLE_PRECISION, status, ierr)
                                end if
                                if (itestr .eq. -5) & !  read(8) reduce
                                     &call MPI_File_read_all(quantcont%fp&
                                     &, reduce, size(reduce), MPI_DOUBLE_PRECISION, status, ierr)
                                detmat_c = detmat_proj
                            end if
                        end if ! endif io_level
                    end if ! endif rankrep

                    if (yescomm) then
                        call bcast_real(rion, size(rion), 0, commrep_mpi)
                        call bcast_real(velion, size(velion), 0, commrep_mpi)
                        if (iopt .ne. 3) then
                            call bcast_real(alphavar, size(alphavar), 0, commrep_mpi)
                            if (contraction .ne. 0) then
                                call bcast_real(detmat_c, size(detmat_c), 0, commrep_mpi)
                            else
                                call bcast_real(detmat, size(detmat), 0, commrep_mpi)
                            end if
                            if (npar_eagp .gt. 0) call bcast_real(eagp_pfaff, size(eagp_pfaff), 0, commrep_mpi)
                            call bcast_real(projm, size(projm), 0, commrep_mpi)
                            call bcast_real(mu_c, size(mu_c), 0, commrep_mpi)
                            call bcast_real(jasmat, size(jasmat), 0, commrep_mpi)
                            call bcast_real(jasmatsz, size(jasmatsz), 0, commrep_mpi)
                            call bcast_real(jasmat_c, size(jasmat_c), 0, commrep_mpi)
                            call bcast_real(jasmatsz_c, size(jasmatsz_c), 0, commrep_mpi)
                            if (allocated(muj_c)) call bcast_real(muj_c, size(muj_c), 0, commrep_mpi)
                            call bcast_real(vjur, size(vjur), 0, commrep_mpi)
                            call bcast_real(dupr, size(dupr), 0, commrep_mpi)
                            if (detc_proj) then
                                call bcast_real(detmat_proj, size(detmat_proj), 0, commrep_mpi)
                                call bcast_real(projmat_c, size(projmat_c), 0, commrep_mpi)
                            end if
                            if (itestr .eq. -5) call bcast_real(reduce, size(reduce), 0, commrep_mpi)
                        end if
                    end if
                    if (yesdetmatc) then
                        call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                             &, nelcol_c, detmat, detmat_c, mu_c, psip)
                    end if
                end if

                if (io_level .eq. 1) then
                    rewind (9)
                    read (9) (((kel(kk, jj), kk=1, 3), jj=(ii - 1)*nrnel + 1, (ii - 1)*nrnel + nel), ii=1, in1)&
                         &, angle, epscutur, epstlur, countav, countt, nacc, nmovet, avenernum, avenerden
                    if (kaverage .and. allocated(detmat_proj) .and. iopt .ne. 3) then
                        read (9) detmat_proj, projmat_c
                    end if
                    if (nbra_cyrus .gt. 0) then
                        read (9) queue_cyrus, first_cyrus
                    end if
                    !            read(9) (kel(1:3,(ii-1)*nrnel+1:(ii-1)*nrnel+nel),ii=1,in1),angle
                else if (io_level .eq. 2) then
                    kelcont%disp = 0
                    ! couting the kelcont size per MPI task to create the right view.
                    kelcont_size = nel3*in1 + size(angle) + 8
                    if (kaverage .and. allocated(detmat_proj) .and. iopt .ne. 3) then
                        kelcont_size = kelcont_size + size(detmat_proj) + size(projmat_c)
                    end if
                    if (nbra_cyrus .gt. 0) then
                        kelcont_size = kelcont_size + in1*(dim_cyrus + 1)
                    end if
                    if (rank .eq. 0) write (6, '(A, A, /, A, I10, A)')&
                        & " Reading electronic configurations from ", trim(kelcont%name), &
                        & " Each MPI task reads ", kelcont_size, " double precision data."
                    call mpiio_file_create_view(kelcont, kelcont_size, MPI_DOUBLE_PRECISION)
                    call mpiio_file_reset_view(kelcont)
                    do ii = 1, in1
                        call MPI_File_read_all(kelcont%fp, kel(1, (ii - 1)*nrnel + 1), nel3,&
                             & MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                    call MPI_File_read_all(kelcont%fp, angle, size(angle), MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, epscutur, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, epstlur, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, countav, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, countt, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, nacc, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, nmovet, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, avenernum, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_read_all(kelcont%fp, avenerden, 1, MPI_DOUBLE_PRECISION, status, ierr)
                    if (kaverage .and. allocated(detmat_proj) .and. iopt .ne. 3) then
                        call MPI_File_read_all(kelcont%fp, detmat_proj, size(detmat_proj), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_read_all(kelcont%fp, projmat_c, size(projmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                    end if
                    if (nbra_cyrus .gt. 0) then
                        call MPI_File_read_all(kelcont%fp, queue_cyrus, in1*dim_cyrus, MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_read_all(kelcont%fp, psip, in1, MPI_DOUBLE_PRECISION, status, ierr)
                        first_cyrus(1:in1) = psip(1:in1)
                    end if
                end if
                !write(6,*) "check Read kel 1", rank, kel(1:3,1)
#else

                read (11) (((kel(kk, jj), kk=1, 3), jj=(ii - 1)*nrnel + 1, (ii - 1)*nrnel + nel), ii=1, in1)&
                    &, angle, epscutur, epstlur, countav, countt&
                                 &, nacc, nmovet, avenernum, avenerden
                if (nbra_cyrus .gt. 0) read (11) queue_cyrus, first_cyrus
#endif
            end if

            call checkiflagerr(iflagerr, rank, 'Error in main')
            if (rank .eq. 0) write (6, *) ' All files are read correctly ...'
#ifdef PARALLEL
            if (contraction .gt. 0) &
                 &call mpi_bcast(allowcontr, 2*nelorb_c, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call bcast_real(wconfn, nw, 0, MPI_COMM_WORLD)
            call bcast_real(wcorw, nfatr, 0, MPI_COMM_WORLD)
            if (.not. yesquantum) then
                call bcast_real(dek, ieskindim, 0, MPI_COMM_WORLD)
                call bcast_real(dekg, ieskingdim, 0, MPI_COMM_WORLD)
            end if
            if (detc_proj .and. .not. yesquantum .and. .not. kaverage) then
                call bcast_real(detmat_proj, size(detmat_proj), 0, MPI_COMM_WORLD)
                call bcast_real(projmat_c, size(projmat_c), 0, MPI_COMM_WORLD)
            end if
            if (itestr .eq. -5) then
                if (.not. yesquantum) then
                    call bcast_real(velion, ieskindim*3, 0, MPI_COMM_WORLD)
                    call bcast_real(reduce, ncgdim*npdim, 0, MPI_COMM_WORLD)
                end if
                call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
                if (rank .eq. 0 .and. idyn .gt. 2 .and. idyn .ne. 5) write (6, *) &
                     &' Temp velocities read =', dnrm2(ieskin, velion(3, 1), 3)**2/ieskin*ris(2)
            end if
#endif

#ifdef PARALLEL
            !   The old covariance matrix, if allocated, is always common to all processors
            !   even in the quantum case.
            if (allocated(cov_old)) then
                if (rank .eq. 0) write (6, *) ' Warning read old covariance matrix '
                call bcast_real(cov_old, size(cov_old), 0, MPI_COMM_WORLD)
            end if
#endif

            alphab = 0.d0

            if (contraction .ne. 0) then
                call update_projm_

                if (yesdetmatc) call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                     &, nelcol_c, detmat, detmat_c, mu_c, psip)
            end if

            if (nfat .ne. nfatr) then
                do i = 1, nfat
                    wcorw(i) = 1.d0
                end do
            end if

            ng = ngenc
            if (iopt .eq. 2) then
                ng = 0
                ngs = 0
            end if
            if (iopt .ne. 2 .and. (epscuttyper .eq. epscuttype)) then
                !      read epscut from file
                epscutu = epscutur
                epstlu = epstlur
            else
                epscutu = epscut
                epstlu = epstl
            end if
            if (rank .eq. 0) write (6, *) ' Used epscut,epstl,tstep =', epscutu, epstlu, tstep

            if (rank .eq. 0) then
                if (itestr .eq. -5) write (6, *) ' Cut off used (parr) for sr when kl=6,7 =', parr
                write (6, *) ' Number of QMC iterations  in this run', ngen
                if (itestr .eq. -5) write (6, *) ' Number of previous optimization steps', ng
                write (6, *) ' Number of previous QMC iterations ', ngs
            end if

            ! ngs number of optimization steps, ng number of QMC iterations.
            ! rewind the old file to continue to write if iopt.eq.0

            ! fine if iopt=0 or iopt=2
        end if

        call read_alphavar

        if (rank .eq. 0 .and. ieskin .ne. 0) then
            write (6, *) ' Warning rescaling acc forces ion  in Hartree'
        end if

        if (idyn .eq. 0) delta0 = 0.d0

        if (idyn .eq. 4) then
            if (delta0 .lt. 4.d0/3.d0*dt*normcorr) then
                !          if(delta0.lt.2.d0*dt) then
                !          delta0=2.d0*dt
                delta0 = 4.d0/3.d0*dt*normcorr
                if (rank .eq. 0) write (6, *) ' Warning  delta_0 >= 4/3 dt !!! Changed to ' &
                     &, delta0
            end if
        end if

        if (idyn .eq. 6 .and. .not. yesrootc) then
            if (delta0 .lt. dt*normcorr) then
                !          if(delta0.lt.2.d0*dt) then
                !          delta0=2.d0*dt
                delta0 = dt*normcorr
                if (rank .eq. 0) write (6, *) ' Warning  delta_0 >= dt !!! Changed to ' &
                     &, delta0
            end if
        end if

        if (idyn .eq. 3 .or. idyn .eq. 7 .or. idyn .eq. 8) then
            if (yesrootc) then
                if (temp .ne. 0.d0 .and. delta0 .lt. 2.d0*dsqrt(max(normcorr, 0.d0)/temp)) then
                    !        delta0=2.d0*dsqrt(normcorr/temp)
                    if (rank .eq. 0) write (6, *) ' Warning  delta_0 >  2/sqrt(T) !!! Change to &
                         & if you see many Warning, alias reduce tion ', 2.d0*dsqrt(normcorr/temp)
                end if
            else
                if (delta0 .lt. dt*normcorr) then
                    delta0 = dt*normcorr
                    if (yessecond) delta0 = delta0/2.d0
                    if (rank .eq. 0) write (6, *) ' Warning  delta_0 >  dt !!! Changed to ', delta0/dt, 'dt'
                end if
            end if
        end if

        if (rank .eq. 0 .and. itestr .eq. -5) then
            write (6, *) ' Initial scaling '
            do i = 1, np
                if (scalpar(i) .ge. 0.d0) write (6, *) i, scalpar(i)
            end do
        end if
        if (typedyncell .eq. 2 .and. fixa) scalpar(np - 2) = 0.d0
        if (typedyncell .eq. 2 .and. fixb) scalpar(np - 1) = 0.d0
        if (typedyncell .eq. 2 .and. fixc) scalpar(np) = 0.d0

        srforce = 0.d0
        srforcew = 0.d0
        wforce = 0.d0
        wforcew = 0.d0

!          walkers with weight wconf=1
!   if etry=0 set etry to the classical energy per site
        if (iopt .eq. 2 .or. nbra .ne. nbrar) then
            nacc = 0.d0
            nmovet = 0.d0
        end if
        do j = 1, nws
            naccm(j) = 1
        end do
        sumdiff = 0.d0

        time_ratiovar = 0.d0
        time_uptabtot = 0.d0
        time_meas = 0.d0
        timemc = 0.d0
        timeopt = 0.d0
        timescra = 0.d0
        timewf = 0.d0
        timepip = 0.d0
        if (itest .eq. 2) then
            lambda = dabs(-etry)
        else
            lambda = -etry
            if (rank .eq. 0) write (6, *) ' lambda chosen =', lambda
        end if
        ngg = ng
        ngn = ngs
        do j = ist, ien
            wconfw(j - istm) = wconfn(j)
        end do
        rweight = nweight
        if (rank .eq. 0 .and. itestr .eq. -5) then
            write (6, *) ' nweight before starting ', nweight
        end if
!       write(6,*) 'before  nweight nweightr # proc  ='
!    1,nweight,nweightr,rank+1

        ndone = nweight - ibinit
        lbin = ndone/nbinr

        if (lbin .eq. 0) then
            write (errmsg, *) ' ERROR # bins larger than nweight,&
                 & reduce nbinr in input !!! <=', ndone
            call checkiflagerr(1, rank, errmsg)
        end if
        flag = .false.
        if (nbinr*lbin .ne. ndone .and. itestr .eq. -5) flag = .true.
        ndone = nbinr*lbin
        nweight = ndone + ibinit
        if (rank .eq. 0 .and. flag)                                                 &
           &write (6, *) ' Warning nweight changed  !!!', nweight
        i_main = iend
        if ((iopt .ne. 0 .and. iopt .ne. 3) .or. yeschange_nweight) then
            inext = iend + nweight
            pippoc = 1
            pippo = 1
        end if

#ifdef  PARALLEL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

        if (iesm + abs(iesfree) + abs(iesinv) + abs(iessw) + iesup + ieskin + iesking .gt. 0) then
            iesmeas = .true.
        else
            iesmeas = .false.
        end if

!       allocating the minimum memory for ipsip

        if (itestr .ne. -5) then
            iscraipsip = nelup_mat + 9*nion
            if (iesbra) iscraipsip = max(iscraipsip, 2*nw)
            deallocate (ipsip)
            allocate (ipsip(iscraipsip))
            ipsip = 0
        end if

#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
        isdistp = 20*indt + 20
        nelorbjmax = max(nelorbj, 1)
        indtmax = max(indt, 1)
        nshelljmax = max(nshellj, 1)

        iscramaxb = isdistp
        iscramaxb = max(iscramaxb, (nion*(nion + 1))/2)

!  ALLOCATE all required for reverse mode.
        if (iessz) then
            iscramaxb = max(iscramaxb, 2*max(nel, indt + 5) + 2*max(nelorbjh, 1) + 20*(indt + 1)&
                 & + max(nelorbjh, 1)*(2*indt + 11)) + nel + nion
        else
            iscramaxb = max(iscramaxb, 2*max(nel, indt + 5) + max(ipj*nelorbjh, 1) + 20*(indt + 1)&
                 & + max(nelorbjh, 1)*(2*indt + 11)) + nel + nion
        end if
!  task9
        if (iessz) then
            isdistp = 27*(indt + 1)*nshelljmax + nelorbjmax*(indt + 6) + 1 + (7 + indt)*nelorbjmax + 2*max(nel, indt + 5)
        else
            isdistp = 27*(indt + 1)*nshelljmax + nelorbjmax*(indt + 6) + 1 + (ipj + 5 + indt)*nelorbjmax + 2*max(nel, indt + 5)
        end if
        iscramaxb = max(iscramaxb, isdistp)

!   NB in TASK8_B reverse mode psip is not used

        iscramaxb = max(iscramaxb, nparshellmax) ! task2 nparshellmax defined  after read_pseudo -->

        iscramaxb = max(iscramaxb, 2*ipc*nelup_mat*nelorb_c) ! task5 upper bound

        iscramaxb = max(iscramaxb, 2*nelorbj_c*nel) ! task5 Jastrow
        iscramaxb = max(iscramaxb, (nion*(nion - 1))/2, ipc*nelup*nelup)
        iscramaxb = max(iscramaxb, ipc*(indt + 5)*nelorbh + nelorbh*(indt + 6) + 27*(indt + 1)*nshell + nelorbh) ! task10

!   iscramaxold=iscramax
!   if(yesfast.ne.0) then
!      iscramax=max(iscramax,nmolfn*(nelup+nelorbh))
!   else
!      iscramax=max(iscramax,nelorbh*nelorbh)
!   endif
!   iscramax=max(iscramax,nelorbjh*nelorbjh)

!   if(iscramaxold.lt.iscramax) then
!      deallocate(psip)
!      allocate(psip(iscramax))
!      if(rank.eq.0) write(6,*) ' Warning increasing allocation scratch for&
!           & psip', iscramax
!   endif

        if (ieskint .ne. 0) then
            if (rank .eq. 0) write (6, *) ' DIFFERENTIAL WARP ALGORITHM FOR FORCES  '
        end if

!        write(6,*) ' Initial tabs AGP '
!        do ii=1,nshell
!        write(6,*) ii,indpar_tab(ii),indshell_tab(ii),indorb_tab(ii)
!        enddo

!        write(6,*) ' Initial tabs Jastrow  '
!        do ii=1,nshellj
!        write(6,*) ii,indparj_tab(ii),indshellj_tab(ii),indorbj_tab(ii)
!        enddo

!         stop

        if (rank .eq. 0 .and. .not. yespulay) &
           &write (6, *) ' Warning calculation of pulay suppressed '

        if (ieskin .ne. 0) then
            ieskinion = ieskin - (ieskint - ieskinr_pos)
            if (rank .eq. 0) write (6, *) ' Number of components ions =', ieskinion
        end if
        nelsquare = nel*nel
        if (itest .ne. 2) then

            allocate (vpotsav_ee(nelsquare, nws))
            if (yes_complex) then
                if (molecular .gt. 0 .and. ipf .eq. 2) then
                    allocate (winvfn(2*nmolfn*(indt + ip4)*nel*nws), winvbarfn(nmolfn*2*nel_mat*nws))
                    nel2wtfn = 2*nmolfn*nel*(indt + ip4)
                    nel2barfn = 2*nel_mat*nmolfn

                else
                    allocate (winvfn(2*nmolfn/ipf*(indt + ip4)*nel*nws), winvbarfn(nmolfn*2*nel_mat*nws))
                    nel2wtfn = 2*nmolfn/ipf*nel*(indt + ip4)
                    nel2barfn = 2*nel_mat*nmolfn
                end if
            else
                if (molecular .gt. 0 .and. ipf .eq. 2) then
                    allocate (winvfn(nmolfn*(indt + ip4)*nel*nws), winvbarfn(nmolfn*nel_mat*nws))
                    nel2wtfn = nmolfn*nel*(indt + ip4)
                    nel2barfn = nel_mat*nmolfn
                else
                    allocate (winvfn(nmolfn/ipf*(indt + ip4)*nel*nws), winvbarfn(nmolfn*nel_mat*nws))
                    nel2wtfn = nmolfn/ipf*nel*(indt + ip4)
                    nel2barfn = nel_mat*nmolfn
                end if
            end if
            !        indbarfn=nel2barfn*(j-1)+1

        else
            nel2wtfn = 0
            nel2barfn = 0
            allocate (vpotsav_ee(1, 1))
            if (yes_complex) then
                allocate (winvfn(2), winvbarfn(2))
            else
                allocate (winvfn(1), winvbarfn(1))
            end if
        end if
        winvbarfn = 0.d0
        winvfn = 0.d0
        vpotsav_ee = 0.d0

! by E. Coccia (8/11/10): read_cube and spline interpolation
        if (ext_pot) then
            call extpot_read
            ! by E. Coccia (4/1/11): spline interpolation for the forces
            ! geometry optimization (idyn.ne.0)
            ! only for AD
            if (idyn .ne. 0) then
                call forces_interpolate()
            end if
        end if
! by E. Coccia (10/12/11): read restr.dat
        if (mm_restr) then
            call restr_read()
        end if

!by E. Coccia (20/12/11): writing electronic random walk
        if (rank .eq. 0 .and. write_rwalk) open (671, file='rwalk.xyz')

        yeswrite = .false.
        yeswrite12 = .false.
!       definition jbrasymiesup

        if (symiesup) then
            allocate (jbrasymiesup(iesupr_c))
            jbrasymiesup = 0
            ind = 0
            do jj = 1, iesupind
                ind = ind + 1
                ii = jbraiesup_sav(ind)
                do j = 1, abs(ii)
                    if (jbraiesup_sav(ind + j) .gt. 0) then
                        jbrasymiesup(jbraiesup_sav(ind + j)) = jj
                        if (yeszagp) jbrasymiesup(jbraiesup_sav(ind + j) + iesup_c) = jj
                    else
                        jbrasymiesup(-jbraiesup_sav(ind + j)) = -jj
                        if (yeszagp) jbrasymiesup(-jbraiesup_sav(ind + j) + iesup_c) = -jj
                    end if
                end do
                ind = ind + abs(ii)
            end do
        end if

        if (scalepulay .ne. 1.d0 .and. rank .eq. 0) &
           &write (6, *) ' Warning biased Pulay scheme !!!  ', scalepulay
        dt4 = 1.d0
        if (rank .eq. 0 .and. idyn .ne. 8 .and. idyn .ne. 0) then
            if (ref_atom .ne. 0) then
                write (6, *) 'Reference unit mass atom (in your fort.10 in ascending order) ', ref_atom
                if (ldynsecond) then
                    write (6, *) ' Then your real time step used in dynamics is (a.u.)=', dt&
                       & *sqrt(atom_weight(nint(atom_number(ref_atom)))/mass_unit*scale_mass*2.d0)
                    write (6, *) 'Remind 1a.u.= 0.0242fs'
                elseif (idyn .ne. 5) then
                    write (6, *) ' Then your real time step used in dynamics is (a.u.)=', dt&
                       & *atom_weight(nint(atom_number(ref_atom)))/mass_unit*scale_mass*2 ! the factor 2 is for Rydberg unit
                    write (6, *) 'Remind 1a.u.= 0.0242fs'
                else
                    write (6, *) ' Then your real time step used in dynamics is (a.u.)=', dt&
                       & *atom_weight(nint(atom_number(ref_atom)))/mass_unit*scale_mass/2 ! the factor 2 is for Rydberg unit
                    write (6, *) 'Remind 1a.u.= 1H'
                end if
            end if
        elseif (rank .eq. 0 .and. idyn .eq. 8) then
            write (6, *) 'Your time units are a.u.!!!'
            write (6, *) 'Your time step is dt=', dt*2.d0*0.0241888, 'fs'
        end if

        if (yesquantum) then
            temp = temp*nbead
            allocate (rionall(3, nion, nbead))
            allocate (fbead(3, nion))
            allocate (mass_ion(3, nion))
            rionall = 0.d0
            fbead = 0.d0
            mass_ion = 0.d0
            dt4 = 4.d0/temp
            if (rank .eq. 0) write (6, *) ' Mass particles (unit m_e) '

            do ii = 1, 3
                do jj = 1, nion
                    mass_ion(ii, jj) = atom_weight(nint(atom_number(jj)))/mass_unit*scale_mass
                    if (rank .eq. 0 .and. ii .eq. 1) write (6, *) jj, mass_ion(ii, jj)
                    !         Use the same mass in the classical dynamics
                end do
            end do

            if (yesturboq) then
                allocate (kdyn(nbead, nbead), kdyn_eig(nbead))
                kdyn = 0.d0
                kdyn_eig = 0.d0
                do ii = 1, nbead
                    kdyn(ii, mod(ii, nbead) + 1) = -1.d0
                    kdyn(mod(ii, nbead) + 1, ii) = -1.d0
                    kdyn(ii, ii) = 2.d0
                end do

                if (idyn .eq. 8) then
                    kdyn = kdyn*temp**2
                    do ii = 1, 3
                        do jj = 1, nion
                            ind = ii + (jj - 1)*3
                            scalpar(np - ieskin + ind) = 2.d0*scalpar(np - ieskin + ind)/mass_ion(ii, jj)
                            ! The factor 2 comes from the Rydberg units
                        end do
                    end do
                else
                    kdyn = kdyn*temp**2*0.5d0*mass_ion(1, 1) ! The factor 1/2 comes from the Rydberg units
                    !          setting all the quantum masses equal to the first one
                    do ii = 1, 3
                        do jj = 1, nion
                            ind = ii + (jj - 1)*3
                            if (ind .ne. 1) scalpar(np - ieskin + ind) = scalpar(np - ieskin + 1)*mass_ion(1, 1)/mass_ion(ii, jj)
                        end do
                    end do
                end if

                if (size(psip) .lt. 3*nbead) then
                    iscramax = 3*nbead
                    deallocate (psip)
                    allocate (psip(iscramax))
                    psip = 0.d0
                end if

                call dsyev('V', 'L', nbead, kdyn, nbead, kdyn_eig&
                     &, psip, 3*nbead, info)
    !!         Identity test
                !          kdyn=0.d0
                !          do ii=1,nbead
                !          kdyn(ii,ii)=1.d0
                !          kdyn_eig(ii)=0.d0
                !          enddo
#ifdef     PARALLEL
                !          All should have exactly the same matrix to avoid roundoff
                call bcast_real(kdyn_eig, nbead, 0, MPI_COMM_WORLD)
                call bcast_real(kdyn, size(kdyn), 0, MPI_COMM_WORLD)
#endif
                if (rank .eq. 0) then
                    write (6, *) ' Eigenvalues elastic quantum term '
                    do ii = 1, nbead
                        write (6, *) ii, kdyn_eig(ii)
                    end do
                end if
            end if

        elseif (idyn .eq. 8) then

            if (rank .eq. 0) write (6, *) ' Classical dynamics idyn 8'
            if (rank .eq. 0) write (6, *) ' Mass particles (unit m_e) '
            ! classical dynamics in a.u. (no mass reference!)
            allocate (mass_ion(3, nion))

            do ii = 1, 3
                do jj = 1, nion
                    mass_ion(ii, jj) = atom_weight(nint(atom_number(jj)))/mass_unit*scale_mass
                    if (rank .eq. 0 .and. ii .eq. 1) write (6, *) jj, mass_ion(ii, jj)
                    !         Use the same mass rescaling as in the quantum dynamics for idyn.eq.8
                end do
            end do

            do ii = 1, 3
                do jj = 1, nion
                    ind = ii + (jj - 1)*3
                    scalpar(np - ieskin + ind) = 2.d0*scalpar(np - ieskin + ind)/mass_ion(ii, jj)
                    ! The factor 2 comes from the Rydberg units
                end do
            end do

        end if

        dimfk = nbin
        perbin = 1
        ndimj = iesinv + iesm + iesd + iesfree
        ndimjp = ndimj + 1
        ndims = iesinv + iesm + iesd + iesfree + iessw
        ndimsp = ndims + 1
        ndimiesup = iesinv + iesm + iesd + iesfree + iessw + iesup

        if (ieskint .ne. 0) then
            allocate (rion_write(3, nion))
            rion_write = rion
        end if
        celldm_write(1:3) = celldm(1:3)
        rs_write = rs
        if (iesup + iessw + iesking + ieskint .gt. 0) then
            if (rank .eq. 0) write (6, *) ' Warning AAD  determinant or forces !!!'
            yesdodet = .true.
        else
            yesdodet = .false.
        end if
        if (iesup + iessw .gt. 0) then
            if (rank .eq. 0) write (6, *) ' Warning AAD  determinant !!!'
            yesdodet_nof = .true.
        else
            yesdodet_nof = .false.
        end if
        yes_hessc = .false.
        if (ipc .eq. 2 .and. yesdodet_nof .and. itestr .eq. -5) then
            if (rank .eq. 0) write (6, *) ' Warning complex Hessian ! '
            yes_hessc = .true.
        end if
        if (iesup .gt. 0 .and. yeszagp) then
            if (rank .eq. 0) write (6, *) ' Warning AAD Z  determinant !!!'
        end if
        if (iesm .gt. 0 .and. yeszj) then
            if (rank .eq. 0) write (6, *) ' Warning AAD Z  Jastrow !!!'
        end if
        if (iesm + iesd + iessw + iesup + iesking + iesinv + iesfree .gt. 0) then
            someparameter = .true.
            if (rank .eq. 0) write (6, *) ' Warning AAD for some parameter more than  energy !!!'
        else
            someparameter = .false.
        end if
        if (iessw + iesup + iesking .gt. 0) then
            someparameterdet = .true.
        else
            someparameterdet = .false.
        end if
        yesprojm = .true.

        if (ip_reshuff .eq. -1) then
            if (yes_ipreshuff) then
                ip_reshuff = 2
                nelorbh_ip = nelorbh*2
                nmol_ip = nmolfn/2
                if (rank .eq. 0) write (6, *) ' Warning ip_reshuff set to 2 '
            else
                ip_reshuff = 1
                nmol_ip = nmolfn
                nelorbh_ip = nelorbh*ipf
            end if
        else
            if (ip_reshuff .gt. 2) ip_reshuff = 2
            if (ip_reshuff .lt. 1) ip_reshuff = 1
            if (rank .eq. 0) write (6, *) ' Warning forced ip_reshuff =', ip_reshuff
            if (ip_reshuff .eq. 1) then
                nmol_ip = nmolfn
                nelorbh_ip = nelorbh
            else
                nelorbh_ip = nelorbh*2
                nmol_ip = nmolfn/2
            end if
        end if
        if (rank .eq. 0) write (6, *) ' nmol_ip nelorbh_ip used =', nmol_ip, nelorbh_ip
        if (someparameterdet .and. yes_complex) then
            if (rank .eq. 0) then
                if (yes_correct) then
                    write (6, *) ' Warning corrected derivatives for complex case '
                else
                    write (6, *) ' Warning corrected complex derivatives with REAL algorithm '
                end if
            end if
        elseif (.not. yes_complex) then
            yes_correct = .false.
        end if
        if (yeszagp .and. ipc .eq. 2) then
            allocate (zagp_imag(iesup))
            zagp_imag = 0.d0
        end if

        allocate (allowed_par(ndimiesup + ieskin))
        allowed_par(:) = .true.
        do i = 1, ndimiesup + ieskin
            if (scalpar(i) .eq. 0.d0) then
                if (rank .eq. 0) write (6, *) ' Not allowed par =', i, scalpar(i)
                allowed_par(i) = .false.
            end if
        end do

        if (iesd .ne. 0 .and. iesdtwobodyoff) then
            !       the two body may have more than one parameter
            min_pointvj = -1
            do jj = 1, nion
                if (pointvj(1, jj) - 1 .lt. min_pointvj .or. min_pointvj .eq. -1) &
                     & min_pointvj = pointvj(1, jj) - 1
            end do
            if (min_pointvj .gt. 1) then
                allowed_par(1 + iesm + iesinv:iesm + iesinv + min_pointvj) = .false.
            else
                allowed_par(1 + iesm + iesinv) = .false.
            end if
        end if
        if (iesd .ne. 0 .and. iesdonebodyoff) then
            !       the two body may have more than one parameter
            allowed_par(iesm + iesinv + pointvj(1, 1):iesm + iesinv + iesd) = .false.
        end if

        if (iesd .ge. 2 .and. noopt_onebody) then
            allowed_par(2 + iesm + iesinv:iesd + iesm + iesinv) = .false.
        end if

        if (ipc .eq. 2) then
            if (symmagp) then
                if (yes_correct) then
                    do kk = ndimjp, ndims, 4
                        allowed_par(kk + 1) = .false.
                        if (real_agp) allowed_par(kk + 2) = .false.
                        allowed_par(kk + 3) = .false.
                    end do
                elseif (real_agp) then
                    do kk = ndimjp, ndims, 2
                        allowed_par(kk + 1) = .false.
                    end do
                end if
                !    Vanishing imaginary part with srcomplex=.true. also if yes_hermitian=.false.

                if (srcomplex .or. yes_hermite .or. contraction .eq. 0 .or. real_agp) then
                    !   Imaginary part of exponent and contracted orbitals set to zero
                    if (.not. real_contracted) then
                        if (rank .eq. 0) write (6, *) &
                             &' Warning real_contracted forced to .true. in this case'
                        real_contracted = .true.
                    end if
                    do kk = 1, iesup_c - 2*ipf*nelorbh*molecular
                        if (dup_c(2*kk) .ne. 0.d0) then
                            if (rank .eq. 0.d0) write (6, *) ' NON zero dup_c =', dup_c(2*kk)
                            call error(' Initializeall ', ' Contracted orbitals &
                                 & should be real in this case. ', 1, rank)
                        end if
                    end do

                    do kk = ndimsp, ndimiesup, 2
                        allowed_par(kk + 1) = .false.
                        if (contraction .eq. 0 .and. .not. yeszagp) allowed_par(kk) = .false.
                    end do
                end if
                if (contraction .ne. 0 .and. iesup .ne. 0) then
                    !         Vanish Imaginary part of Z exponents only
                    do kk = 1, iesupr_2
                        if (iesuptrans(kk) .ne. 0) then
                            jj = abs(jbraiesup(iesuptrans(kk)))
                            if (jj .ne. 0) then
                                allowed_par(2*jj + ndims) = .false.
                                if (.not. yeszagp) allowed_par(2*jj - 1 + ndims) = .false.
                            end if
                        end if
                    end do
                end if
            else ! symmagp
                !         Vanish Imaginary part of Z forces only

                if (real_agp) then
                    do kk = ndimjp, ndims, 2
                        allowed_par(kk + 1) = .false.
                    end do
                end if

                if (iesup .ne. 0) then
                    if (contraction .ne. 0 .and. .not. real_contracted) then
                        do kk = 1, iesupr_2
                            if (iesuptrans(kk) .ne. 0) then
                                jj = abs(jbraiesup(iesuptrans(kk)))
                                if (jj .ne. 0) then
                                    allowed_par(2*jj + ndims) = .false.
                                    if (.not. yeszagp) allowed_par(2*jj - 1 + ndims) = .false.
                                end if
                            end if
                        end do
                    else
                        do kk = 1, iesup_c - 2*ipf*nelorbh*molecular
                            if (dup_c(2*kk) .ne. 0.d0) then
                                if (rank .eq. 0.d0) write (6, *) ' NON zero dup_c =', dup_c(2*kk)
                                call error(' Initializeall ', ' Contracted orbitals &
                                     & should be real in this case. Use real_contracted=.false. ! ', 1, rank)
                            end if
                        end do
                        do kk = ndimsp, ndimiesup, 2
                            allowed_par(kk + 1) = .false.
                        end do
                    end if
                end if
            end if ! symmagp
            !        end Z and contracted orbitals begin lambda matrix
            if (yes_hermite .and. symmagp .and. ndims .gt. ndimj) then
                !  Vanishing imaginary  part of diagonal terms forces
                do kk = 1, nnozero_c
                    jj = abs(jbradet(kk))
                    iy = (nozero_c(kk) - 1)/nelorb_c + 1
                    ix = nozero_c(kk) - (iy - 1)*nelorb_c
                    if (ix .eq. iy .and. jj .ne. 0) then
                        if (yes_correct) then
                            allowed_par(4*jj - 2 + ndimj:4*jj + ndimj) = .false.
                        else
                            allowed_par(2*jj + ndimj) = .false.
                        end if
                    end if
                end do
            end if
            if (yes_crystal .and. ndims .gt. ndimj) then
                do kk = 1, nnozero_c
                    jj = abs(jbradet(kk))
                    ! There is only one effective parameter in this case
                    if (sjbradet(kk) .and. jj .ne. 0) then
                        if (symmagp .and. yes_correct) then
                            if (yes_hermite .or. no_sjbra) allowed_par(4*jj - 2 + ndimj:4*jj + ndimj) = .false.
                            if (no_sjbra) allowed_par(4*jj - 3) = .false.
                        else
                            if (yes_hermite .or. no_sjbra) allowed_par(2*jj + ndimj) = .false.
                            if (no_sjbra) allowed_par(2*jj - 1) = .false.
                        end if
                    end if
                end do
            end if
        else ! if ipc=2
            if (.not. yeszagp .and. iesup .ne. 0) then
                if (contraction .ne. 0) then
                    do kk = 1, iesupr_2
                        if (iesuptrans(kk) .ne. 0) then
                            jj = abs(jbraiesup(iesuptrans(kk)))
                            if (jj .ne. 0) then
                                allowed_par(jj + ndims) = .false.
                            end if
                        end if
                    end do
                else
                    do kk = ndimsp, ndimiesup
                        allowed_par(kk) = .false.
                    end do
                end if
            end if
        end if ! ipc=2

        if (rank .eq. 0 .and. itestr .eq. -5) then
            write (6, *) ' Not Allowed parameters '
            do kk = 1, ndimiesup + ieskin
                if (.not. allowed_par(kk)) write (6, *) kk, allowed_par(kk)
            end do
        end if

        if (ipc .eq. 2 .and. .not. yes_correct .and. (yesdodet .or. itest .eq. 1)) then
            yes_real = .true.
            if (scaleeloc .eq. -1.d0) scaleeloc = 0.5d0 ! Umrigar choice.
            if (srcomplex) then
                if (rank .eq. 0) write (6, *) ' Warning srcomplex turned to .false. (real var) '
                srcomplex = .false.
            end if
        else
            yes_real = .false.
        end if
        if (((itestr .eq. -5 .and. itest .eq. 1 .and. scaleeloc .ne. 0.d0)&
           &.or. scaleeloc .gt. 0) .and. .not. yes_real .and. ipc .eq. 1) then
            if (rank .eq. 0) write (6, *) ' Warning computing der local energy '
            if (scaleeloc .eq. -1.d0) scaleeloc = 0.5d0 ! Umrigar choice.
            yes_real = .true.
        end if

        if (kaverage) then
            cost = wkp(ikpoint)
#ifdef PARALLEL
            call mpi_allreduce(cost, scale_spsi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            !    The average value of the weight
            scale_spsi = scale_spsi/nproc
            scale_spsi = 1.d0/scale_spsi
#endif
        else
            scale_spsi = 1.d0
        end if

        target_ratio = 0.5d0*(1.d0 - hopfraction)
        minratio = 0.9d0*target_ratio
        maxratio = 1.1d0*target_ratio

        if (in1*nprocrep .gt. 10) then
            ifreqchange = 1
        else
            ifreqchange = nint(10/dble(in1*nprocrep))
        end if
        if (in1*nproc .gt. 10) then
            ifreqchanger = 1
        else
            ifreqchanger = nint(10/dble(in1*nproc))
        end if
        if (.not. fncont) then
            better_dmc = .false.
        else
            allocate (zetamin(nel), distmin(nel))
            zetamin = 0.d0
            distmin = 0.d0
        end if
        if (itestr .ne. -5 .and. iopt .ge. 1 .and. itest .ne. 2 .and. yesnleft .and. changelambda) then
            ibinit = min(nint(0.04999999d0*ngen + 1), 100)
            if (rank .eq. 0) write (6, *) ' Warning starting averaging energies from ', ibinit
        end if
        if ((iopt .eq. 0 .or. iopt .eq. 3) .and. yesnleft .and. changelambda) then
#ifdef PARALLEL
            wtot(1) = avenernum
            wtot(2) = avenerden
            if ((commrep_mpi .ne. mpi_comm_world) .and. decoupled_run) then
                call reduce_base_real_to(2, wtot, psip, commrep_mpi, -1)
            else
                call reduce_base_real_to(2, wtot, psip, MPI_COMM_WORLD, -1)
            end if
            lambda = -psip(1)/psip(2)
#else
            lambda = -avenernum/avenerden
#endif
            if (rank .eq. 0) write (6, *) ' Warning restoring value of trial energy='&
                 &, -lambda*ris(2)
        end if
        nwnk = nw/nk
        nwnkp = nwnk + 1
        wbraw = 1.d0
        if (allocated(jbraw)) then
            !   Forgot to inizialize jbraw
            do kk = 1, nw
                jbraw(kk) = kk
            end do
        end if
        if (yes_hessc) then
            !    Check iscrapip
            iscrapip = iscrapip + 8*(ncg + npbra) ! just to be sure
            if (iscrapip .le. 8*(ncg + npbra) + maxall) then
                iscrapip = 8*(ncg + npbra) + maxall
                write (6, *) ' Warning increasing iscrapip for Hessian complex '
            end if
        end if
        if (rank .eq. 0) write (6, *) ' firstmol nmolfn after all ', firstmol, nmolfn
!  Define always acc_dyn
        acc_dyn = .true.
        weight_vir = 1.d0
        reweight_dmc = 1.d0
#if defined PARALLEL && defined __SCALAPACK
        if (itestr .eq. -5 .and. k6gen .and. kl .eq. -6) then
            call set_env(commsr_mpi)
        end if
#endif
        wdone = .true.
        allocate (t_cyrus(in1))
        t_cyrus = 0.d0
        if (decoupled_run) then
            if (itestr .ne. -5) then
                allocate (wbra_t(nk), ener_t(nk), etot_t(np3, nk), wtotf_t(2, nk))
                etot_t = 0.d0
                wtotf_t = 0.d0
                wbra_t = 0.d0
                ener_t = 0.d0
            end if
            if (iesbra) then
                allocate (wconfn_kps(nwnk))
                wconfn_kps = 0.d0
            end if
        end if
        if (npsa .gt. 0 .and. enforce_detailb) then
            if (nintpsa .eq. 4) then
                call error(' Initializeall ', ' Detailed balance is not possible &
                     & with such a quadrature pseudo mesh  (change nintpsa>4)  ! ', 1, rank)
            end if
        end if
        only_molecular = .false.

        if (cutweight .eq. 0. .or. true_wagner .gt. 0) then
            yes_cutweight = .false.
        else
            yes_cutweight = .true.
        end if
        if (iespbc) then
            cost = 1.d0 - metric(1, 2)**2 - metric(1, 3)**2 - metric(2, 3)**2 + 2.d0*metric(1, 2)*metric(1, 3)*metric(2, 3)
            if (cost .lt. 0.d0) then
                call error(' Initializeall ', ' Metric read is wrong because not positive definite !!! ', 1, rank)
            end if
        end if
        if (nbra_cyrus .gt. 0) then
            read_ok = .true.
            do ii = 1, in1
                if (first_cyrus(ii) .lt. 1 .or. first_cyrus(ii) .gt. nbra_cyrus) read_ok = .false.
            end do
            if (.not. read_ok) call error(' Initializeall ', ' first_cyrus wrong in read or input ', 1, rank)
        end if
        count_zerowf = 0.d0
        count_allwf = 0.d0
        if (yes_sparse .and. rank .eq. 0 .and. .not. iessz .and. contractionj .eq. 0 .and. nelorbjh .ne. 0) then
            write (6, *) ' Warning using SPARSE matrix algorithm for Jastrow '
        end if
        enerdiff = 0.d0 ! just to be sure it is initialized.
    end subroutine Initializeall

    subroutine Finalizeall

        ! by E. Coccia (30/12/10): deallocate arrays for ext_pot
        use extpot, only: ext_pot, mm_restr, write_rwalk
! by E. Coccia (4/2/11): deallocate arrays for van_der_waals
        use van_der_waals, only: vdw
        use qmckl

        implicit none
        real*8, external :: dnrm2
        real*8 drand1, enercont, jacobian, mapping
        integer iend_sav
        integer(kind=qmckl_exit_code) :: rc
#if defined (_OPENMP) && defined (__NOOMP)
        integer, external :: omp_get_max_threads
        call omp_set_num_threads(1) ! scalar code
#endif

#ifdef _QMCKL
        rc = qmckl_context_destroy(qmckl_ctx)
        if (rc.ne.QMCKL_SUCCESS) then
            write (0, *) "Unable to destroy QMCkl context"
        end if
#endif
        

! by E. Coccia (20/12/11): writing electronic random walk
        if (rank .eq. 0 .and. write_rwalk) close (671)

!     at the end write the final  configuration file

        iend_sav = iend

        iend = i_main

#ifdef _CUSOLVER
#ifdef RISC
        call cusolver_handle_destroy_(handle)
#else
        call cusolver_handle_destroy(handle)
#endif
#endif

#ifdef PARALLEL
!  sum the acceptances for each process
        call mpi_reduce(nacc, naccmpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0 &
           &, MPI_COMM_WORLD, ierr)
        call mpi_reduce(nmovet, nmovetmpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0 &
           &, MPI_COMM_WORLD, ierr)
        call mpi_reduce(naccpseudo, naccpseudompi, 1                    &
           &, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(nontr, cost, 1                                  &
           &, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(acclarge, acclargempi, 1                        &
           &, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

! communicate to the master the configurations needed for the
!   final write
!      call mpi_gather(kel(1,nel*(indt+1)*(ist-1)+1),3*nel*(indt+1)*in1 &
!    &,MPI_DOUBLE_PRECISION,kel(1,nel*(indt+1)*(ist-1)+1)               &
!    &,3*nel*(indt+1)*in1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!   psip(1:in1)=wconf(ist:ien)
!   call mpi_gather(psip,in1,MPI_DOUBLE_PRECISION,wconf              &
!        &,in1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        psip(1:in1) = wconfn(ist:ien)
        call mpi_gather(psip, in1, MPI_DOUBLE_PRECISION, wconfn, in1, &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if (rank .eq. 0) then
            !      nacc=naccmpi   ! saved processor accepted moves on a file
            naccpseudo = naccpseudompi
            nontr = cost
            acclarge = acclargempi
        end if

        if (kaverage .and. .not. decoupled_run) then

            ! sum over k-points the values of the inverse wavefunction
            call mpi_reduce(psiav, cost, 1                                    &
                 &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
            call sum_kpoints_scalar_real8(cost, commcolrep_mpi, -1)
            if (rank .eq. 0) psiav = cost

            call mpi_reduce(psisav, cost, 1                                   &
                 &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commrep_mpi, ierr)
            call sum_kpoints_scalar_real8(cost, commcolrep_mpi, -1)
            if (rank .eq. 0) psisav = cost

        else

            call mpi_reduce(psiav, cost, 1                                    &
                 &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMMREP_MPI, ierr)
            if (rank .eq. 0) psiav = cost

            call mpi_reduce(psisav, cost, 1                                   &
                 &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMMREP_MPI, ierr)
            if (rank .eq. 0) psisav = cost

        end if

        call mpi_reduce(counttot, cost, 1                                 &
           &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMMREP_MPI, ierr)
        if (rank .eq. 0) counttot = cost

        call mpi_reduce(countcut, cost, 1                                 &
           &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMMREP_MPI, ierr)
        if (rank .eq. 0) countcut = cost

        call mpi_reduce(countreg, cost, 1                                 &
           &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMMREP_MPI, ierr)
        if (rank .eq. 0) countreg = cost

#endif

        if (yesquantum .and. rankrep .eq. 0 .and. nbead .gt. 1) then
            if (yeswrite10) call write_fort10(27)
        elseif ((kaverage .and. molyes) .and. rankrep .eq. 0 .and. itestr .eq. -5) then
            call write_fort10(27)
        end if

        if (rank .eq. 0) then
            ! only with minimization
            if (itestr .eq. -5) then
                write (6, *) ' Write the parameters of the final wavefunction '
                ! to speed up the writing
                close (10, status='DELETE')
                open (unit=10, file='fort.10', form='formatted', status='unknown')
                call write_fort10(10)

                open (unit=7, file='stop.dat', form='formatted', status='unknown')
                close (unit=7, status='DELETE')
                !          delete the final stop.dat to avoid errors in continuation
            end if
            ! to speed up the writing
            close (11, status='DELETE')
            open (unit=11, file='fort.11', form='unformatted', status='unknown')
            call write_fort11_begin

            if (ireadmin .ne. 1) call write_fort11_end
            if (rank .eq. 0 .and. idyn .ge. 2 .and. idyn .ne. 5) write (6, *) ' Temp velocities write ='&
                 &, dnrm2(ieskin, velion(3, 1), 3)**2/ieskin*ris(2)

        end if ! endif rank.eq.0

#ifdef __KCOMP
123     format(32767e15.7)
#else
123     format(1000000e15.7)
#endif

        if (io_level .ne. 0) then
#ifdef PARALLEL
            if (io_level .eq. 1) then
                rewind (9)
                !        write(9) (kel(1:3,(ii-1)*nrnel+1:(ii-1)*nrnel+nel),ii=1,in1),angle
                write (9) (((kel(kk, jj), kk=1, 3), jj=(ii - 1)*nrnel + 1, (ii - 1)*nrnel + nel), ii=1, in1)&
                     &, angle, epscutu, epstlu, countav, countt, nacc, nmovet, avenernum, avenerden
                if (rank .eq. 0) write (6, *) ' write  epscut =', epscutu, epstlu
                if (kaverage .and. allocated(detmat_proj)) then
                    write (9) detmat_proj, projmat_c
                end if
                if (nbra_cyrus .gt. 0) write (9) queue_cyrus, first_cyrus
            else if (io_level .eq. 2) then

                !        if(nbra_cyrus.gt.0) then
                !           call error( ' Initializeall ', ' nbra_cyrus not allowed with mpiio  &
                !        & ', 1,rank)
                !        endif
                call mpiio_file_set_zero(kelcont)
                kelcont%disp = 0
                ! couting the kelcont size per MPI task to create the right view.
                kelcont_size = nel3*in1 + size(angle) + 8
                if (kaverage .and. allocated(detmat_proj)) then
                    kelcont_size = kelcont_size + size(detmat_proj) + size(projmat_c)
                end if
                if (nbra_cyrus .gt. 0) then
                    kelcont_size = kelcont_size + in1*(dim_cyrus + 1)
                end if

                if (rank .eq. 0) write (6, '(A, A, /, A, I10, A)') " Writing electronic configurations to ", trim(kelcont%name), &
                     & " Each MPI task writes ", kelcont_size, " double precision data."
                call mpiio_file_create_view(kelcont, kelcont_size, MPI_DOUBLE_PRECISION)
                call mpiio_file_reset_view(kelcont)
                do ii = 1, in1
                    call MPI_File_write_all(kelcont%fp, kel(1, (ii - 1)*nrnel + 1), nel3, MPI_DOUBLE_PRECISION, status, ierr)
                end do
                call MPI_File_write_all(kelcont%fp, angle, size(angle), MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, epscutu, 1, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, epstlu, 1, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, countav, 1, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, countt, 1, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, nacc, 1, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, nmovet, 1, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, avenernum, 1, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_File_write_all(kelcont%fp, avenerden, 1, MPI_DOUBLE_PRECISION, status, ierr)
                if (kaverage .and. allocated(detmat_proj)) then
                    call MPI_File_write_all(kelcont%fp, detmat_proj, size(detmat_proj), MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_File_write_all(kelcont%fp, projmat_c, size(projmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                end if
                if (nbra_cyrus .gt. 0) then
                    call MPI_File_write_all(kelcont%fp, queue_cyrus, dim_cyrus*in1, MPI_DOUBLE_PRECISION, status, ierr)
                    psip(1:in1) = first_cyrus(1:in1)
                    call MPI_File_write_all(kelcont%fp, psip, in1, MPI_DOUBLE_PRECISION, status, ierr)
                end if
            end if
            !write(6,*) "check Write kel 1", rank, kel(1:3,1)
            if (allocated(kelind)) deallocate (kelind) ! needed for new compute_eloc_logpsi
            if (yesquantum .and. rankrep .eq. 0) then
                if (io_level .eq. 1) then
                    rewind (8)
                    if (.not. yeswrite10) then
                        if (contraction .ne. 0) then
                            write (8) rion, velion, alphavar, detmat_c, projm, mu_c&
                                 &, jasmat, jasmatsz, jasmat_c, jasmatsz_c, vjur, dupr
                        else
                            write (8) rion, velion, alphavar, detmat, projm, mu_c&
                                 &, jasmat, jasmatsz, jasmat_c, jasmatsz_c, vjur, dupr
                        end if
                        if (allocated(muj_c)) write (8) muj_c
                        if (npar_eagp .gt. 0) write (8) eagp_pfaff
                    else
                        write (8) rion, velion
                    end if
                    if (detc_proj) write (8) detmat_proj, projmat_c
                    if (itestr .eq. -5) write (8) reduce
                elseif (io_level .eq. 2) then
                    call mpiio_file_set_zero(quantcont)
                    call mpiio_file_reset_view(quantcont)
                    if (.not. yeswrite10) then
                        call MPI_File_write_all(quantcont%fp, rion, size(rion), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, velion, size(velion), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, alphavar, size(alphavar), MPI_DOUBLE_PRECISION, status, ierr)
                        if (contraction .ne. 0) then
                            call MPI_File_write_all(quantcont%fp, detmat_c, size(detmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                        else
                            call MPI_File_write_all(quantcont%fp, detmat, size(detmat), MPI_DOUBLE_PRECISION, status, ierr)
                        end if
                        call MPI_File_write_all(quantcont%fp, projm, size(projm), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, mu_c, size(mu_c), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, jasmat, size(jasmat), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, jasmat_c, size(jasmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, jasmatsz, size(jasmatsz), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, jasmatsz_c, size(jasmatsz_c), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, vjur, size(vjur), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, dupr, size(dupr), MPI_DOUBLE_PRECISION, status, ierr)
                        if (allocated(muj_c))& !  read(8) muj_c
                             &call MPI_File_write_all(quantcont%fp, muj_c, size(muj_c), MPI_DOUBLE_PRECISION, status, ierr)
                        if (npar_eagp .gt. 0)& !  read(8) muj_c
                            &call MPI_File_write_all(quantcont%fp&
                            &, eagp_pfaff, size(eagp_pfaff), MPI_DOUBLE_PRECISION, status, ierr)
                    else
                        call MPI_File_write_all(quantcont%fp, rion, size(rion), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, velion, size(velion), MPI_DOUBLE_PRECISION, status, ierr)
                    end if
                    if (detc_proj) then ! read(8) detmat_proj,projmat_c
                        call MPI_File_write_all(quantcont%fp, detmat_proj, size(detmat_proj), MPI_DOUBLE_PRECISION, status, ierr)
                        call MPI_File_write_all(quantcont%fp, projmat_c, size(projmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                    end if
                    if (itestr .eq. -5) & !  read(8) reduce
                         &call MPI_File_write_all(quantcont%fp, reduce, size(reduce), MPI_DOUBLE_PRECISION, status, ierr)
                end if

            end if
            if (rank .eq. 0) then
                nacc = naccmpi ! for write_output
                rata = nmovetmpi ! for write_output
            end if
#else
            write (11) (((kel(kk, jj), kk=1, 3), jj=(ii - 1)*nrnel + 1, (ii - 1)*nrnel + nel), ii=1, in1)&
                &, angle, epscutu, epstlu, countav, countt&
                            &, nacc, nmovet, avenernum, avenerden
            if (nbra_cyrus .gt. 0) write (11) queue_cyrus, first_cyrus
            rata = nmovet
#endif

            if (io_level .eq. 1) call write_seeds(cstr(trim(ranseedfilename)))
#ifdef PARALLEL
            if (io_level .eq. 2) call write_seeds_mpiio(cstr(trim(scratchpath)//'randseed.all'), rank)
#endif
        end if
        iend = iend_sav ! for writing
#ifdef PARALLEL
#ifdef _TIME
        call reduce_base_real(11, timings, mpi_comm_world, 0)
        call reduce_base_real(11, timingsb, mpi_comm_world, 0)
        call reduce_base_real(1, count_zerowf, mpi_comm_world, 0)
        call reduce_base_real(1, count_allwf, mpi_comm_world, 0)
        timings = timings/nproc ! the average time among proc.
        timingsb = timingsb/nproc ! the average time among proc.
#endif
        call reduce_base_real(1, timewf, mpi_comm_world, 0)
        timewf = timewf/nproc
#endif

        if (rank .eq. 0) call write_output

        call close_files(rank)

        if (unreliable .eq. 0) then
            call deallocate_all ! sometimes does not work on fucking infiniband network
            if (allocated(kelsav)) deallocate (kelsav)
            if (allocated(zagp_imag)) deallocate (zagp_imag)
            if (change_tpar .and. rank .eq. 0) &
                 &  deallocate (energy_list, error_energy_list, tpar_buffer_filled)
        end if
        if (yesquantum) then
            deallocate (rionall, fbead, mass_ion)
            if (yesturboq) deallocate (kdyn, kdyn_eig)
        end if
        if (allocated(rion_write)) deallocate (rion_write)
!    if(allocated(rion_sav)) deallocate(rion_sav)

! by. E. Coccia (30/12/10): deallocate arrays for ext_pot
        if (ext_pot) then
            call deallocate_extpot()
            ! by E. Coccia (4/1/11): deallocate arrays for the forces
            ! geometry optimization (idyn.ne.0)
            ! only for AD
            if (idyn .ne. 0) then
                call deallocate_forces()
            end if
            ! by E. Coccia (4/2/11): deallocate arrays for the vdw term
            if (vdw) then
                call deallocate_vdw()
            end if
            ! by E. Coccia (31/5/11): deallocate arrays for link_atom
            if (link_atom) then
                call deallocate_link()
            end if
        end if
        if (mm_restr) then
            call deallocate_restr()
        end if
        if (symiesup) deallocate (jbrasymiesup)
        deallocate (t_cyrus)
        if (nbra_cyrus .gt. 0) then
            deallocate (first_cyrus, queue_cyrus)
        end if
        if (decoupled_run) then
            if (itestr .ne. -5) deallocate (wbra_t, ener_t, etot_t, wtotf_t)
            if (iesbra) deallocate (wconfn_kps)
        end if
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif

    end subroutine Finalizeall

    subroutine makeallmeas_fast
        implicit none
!        This subroutine compute by scratch all correlation functions
!        used for the minimization of the energy or for the dynamic
!        if flagcont=.false. it is a bit faster because does not
!        recompute by scratch the normalization of the coefficients
!        in the orbitals defining the AGP.
!      This object  should be divided in three:
!  a)  calculation of all necessary for energy and few other correlations
!  b)  calculation of all necessary for energy minimization
!  c)  calculation of ionic forces
!  But this inside the main loop over js

! parallel scratch

        do js = ist, ien
            j = js - istm
            indtabj = Ltab*(j - 1) + 1
            indtabbj = Ltabb*(j - 1) + 1
            indkj = nel*(indt + 1)*(j - 1) + 1
            indksj = nel*(j - 1) + 1
            indksij = nelnion*(j - 1) + 1
            indwwjj = nel2wtj*(j - 1) + 1

            indtj = Lztab*(j - 1) + 1
            indtjr = Lztabr*(j - 1) + 1
            indwwj = nel2wt*(j - 1) + 1
            indwwjfn = nel2wtfn*(j - 1) + 1
            indwupj = nel2upt*(j - 1) + 1
            indwdoj = nel2dot*(j - 1) + 1
            indaupj = nel2up*(j - 1) + 1
            indbar = nel2bar*(j - 1) + 1
            indbarfn = nel2barfn*(j - 1) + 1

            indjbar = nel2jbar*(j - 1) + 1
            if (iessz) then
                indjbarsz = indjbar
            else
                indjbarsz = 1
            end if

            timep = cclock()

            if (flagcont) iflagnorm = 3 ! initialize by scratch normalizations

            !        Restore winvbar when not updated

            singdet(j) = .true.

            call upscratch_global(js, pseudologic, iesrandoml)

            if (flagcont .and. developer .eq. -1 .and. rank .eq. 0) then
                write (6, *) ' Jastrowall-ee =', sum(jastrowall_ee(:, :, 0, j))
                write (6, *) ' Jastrowall-ei =', sum(jastrowall_ei(:, :, j))
            end if

            timescra = timescra + cclock() - timep

            naccm(j) = 1 ! Since the upscratch has been already done

            if (itest .eq. 2) then
                psirav = psirav + psidetln(j)*wconfn(js)
                countav = countav + wconfn(js)
                if (kaverage) then
                    countt = countt + wkp(ikpoint)*scale_spsi
                else
                    countt = countt + 1.d0
                end if
                counttot = counttot + wconfn(js)
                psiav = psiav + psidetln(j)*wconfn(js)
                psisav = psisav + psidetln(j)**2*wconfn(js)

            end if

            if (.not. singdet(j)) then
                !if(psisn(j).ne.0) then

                if (yes_complex) then
                    call uptable_complex(nelup, neldo, indtupt, table(indtj), tabler(indtjr)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), tmu(indtabbj)&
                         &, epscutdmc, psiln(j), psidetln(j), costa, parcutg, istart, typereg)

                    call updiag_complex(table(indtj), tmu(indtabbj), diag(j), enert(1, j)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), nelup, neldo, nel&
                         &, indt, indtupt, indteff, istart, vpot(j), vpotreg(1, indksj), parcutg, vcut(j) &
                         &, diffkin(1, j), novar, yescut(j), ener_c, voffpseudo, psip, psip(2*nel + 1))
                    enertrue(j) = real(ener_c)

                    call energy_complex(Lzeff, table(indtj), tabler(indtjr), diag(j), wsto(j), veff&
                         &, veffright, enerdiff, itest, npow, gamma, nel, istart, indtm(1, j), epscutdmc)

                else

                    call uptable(nelup, neldo, indtupt, table(indtj), tabler(indtjr)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), tmu(indtabbj)&
                         &, epscutdmc, psiln(j), psidetln(j), costa, parcutg, istart, typereg)

                    call updiag(table(indtj), tmu(indtabbj), diag(j), enert(1, j)&
                         &, winvup(indwupj), winvdo(indwdoj), tabpip(indtabj), nelup, neldo, nel&
                         &, indt, indtupt, indteff, istart, vpot(j), vpotreg(1, indksj), parcutg, vcut(j) &
                         &, diffkin(1, j), novar, yescut(j), enertrue(j), voffpseudo, psip, psip(nel + 1))

                    call energy(Lzeff, table(indtj), tabler(indtjr), diag(j), wsto(j), veff&
                         &, veffright, enerdiff, itest, npow, gamma, nel, istart, indtm(1, j), epscutdmc)

                end if

                !-------------------------------------------------------------------
                if (add_diff) then
                    enertrue(j) = enertrue(j) + enerdiff
                    enert(1, j) = enert(1, j) + enerdiff
                    if (ipc .eq. 2) ener_c = ener_c + enerdiff
                end if
                diffkin(3, j) = diffkin(3, j) + enerdiff

                if (novar .eq. 3) enertrue(j) = -wsto(j)

                if (itest .eq. 1) then
                    diagfn(j) = diag(j) + (1.d0 + gamma)*veff + npow*(1.d0 + gamma)*veffright
                    if (better_dmc) call cutwstodmc(yesalfe, nel, wsto(j), lambda, diffkin(3, j), psip(ipc*nel + 1), cutreg)
                    if (rsignr .ne. 0.d0) then
                        diagfn(j) = diagfn(j) + rsignr*(diffkin(3, j) + lambda)
                        wsto(j) = wsto(j) + rsignr*(diffkin(3, j) + lambda)
                    end if
                end if

                if (itest .eq. 2) then
                    if (.not. yes_complex) then
                        call comp_econf(j, js, econf, enert, diffkin, vpot, vcut, voffpseudo, table)
                        if (isfix .ne. 0) econf(j + nwfix) = enertrue(j)**2
                    else
                        call comp_econf(j, js, econf, enert, diffkin, vpot, vcut, voffpseudo, tabler)
                        if (isfix .ne. 0) econf(j + nwfix) = ener_c*dconjg(ener_c)
                    end if

                elseif (fncont) then
                    ttry = tbra
                    call upgradcont(gradpsibar(1, indksj), gradpsi(1, indksj)                  &
                         &, indt, nelup, neldo, winvup(indwupj), winvdo(indwdoj)                 &
                         &, tabpip(indtabj), ttry, gradtotbar(j), gradtot(j), kel(1, indkj)       &
                         &, dist(indksij), rion, nion, zetar, LBox)
                else
                    call diffus(nel, indt, gamma, ivic(1, 1, indksj), table(indtj), diffuse(j), istart)
                end if

            end if ! fine if psi(j)>0

        end do ! end do for the walkers

! end parallel scratch

    end subroutine makeallmeas_fast

    subroutine project_v(doinv)
        implicit none
        logical doinv
        integer i, j
        real*8 check, check0, temp
!       compute effective parameters and change alphavar according to reducel
!       scriviamo su file anche i parametri variazionali.

        if (doinv) then
            call dgemm_my('N', 'T', ncg_adr, ncg_adr, kp0, 1.d0, reducel, ncg_adr, reducel&
                 &, ncg_adr, 0.d0, mat_adr, ncg_adr, nprocu, rankrep, commrep_mpi)
            do i = 1, ncg_adr
                if (sum(abs(reducel(i, 1:kp0))) .eq. 0.d0) then
                    mat_adr(i, i) = 1.d0
                end if
            end do
            call dgetrf(ncg_adr, ncg_adr, mat_adr, ncg_adr, ipip_adr, info)
            if (info .ne. 0) then
                write (errmsg, *) ' ERROR singular parametrization &
                     &(check your fort.10 and stdinput)  ', info
                call checkiflagerr(1, rank, errmsg)
            end if
        end if
        call dgemv('N', ncg_adr, kp0, 1.d0, reducel, ncg_adr, alphavar, 1, 0.d0, v_adr, 1)
        call dgetrs('N', ncg_adr, 1, mat_adr, ncg_adr, ipip_adr, v_adr, ncg_adr, info)
        if (rank .eq. 0) write (25, 123) (v_adr(i), i=1, ncg_adr)

#ifdef DEBUG
! check consistency here
        check = 0.d0
        do j = 1, kp0
            if (rpar(j) .ne. 0.d0) then
                temp = 0.d0
                do i = 1, ncg_adr
                    temp = temp + v_adr(i)*reducel(i, j)
                end do
                check0 = abs(temp - alphavar(j))
                if (check0 .gt. 1d-6 .and. rank .eq. 0) write (6, *) ' ERROR in parametr =', j, check0&
                     &, alphavar(j), temp
                check = check + check0
            end if
        end do
        if (rank .eq. 0) write (6, *) ' Error parametrization =', check
#endif
#ifdef __KCOMP
123     format(32767e15.7)
#else
123     format(1000000e15.7)
#endif
    end subroutine project_v

    subroutine project_rmax
        implicit none
        integer j, k, ind, iy, ix
        do j = 1, kp0
            if ((rpar(j) .lt. 0.d0 .and. yescutdet) .or. (rpar(j) .gt. 0.d0 .and. yescutjas)) then
                alphavar(j) = alphavar(j)*smoothcut
                psip(j) = 0.d0
                do k = 1, ncg
                    reduce(k, j) = 0.d0
                end do
            end if
        end do

        if (rmaxj .ne. 0.d0) then
            indc = iesinv + iesm + iesd
            do k = indc + 1, indc + iesfree
                ddw(k - indc) = alphavar(k)
            end do
            ! update jasmat_c
            if (contractionj .ne. 0) then
                call bconstrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c, jasmat_c  &
                     &, ipj*nelorbj_c, ddw)
!   if(ipj.eq.2) then
!      call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!   else
!      call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh&
!           &, nelorbj_c, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!   endif
            else
                if (yes_sparse) then
                    call bconstrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder, jasmat&
                         &, ipj*nelorbjh, ddw)
                else
                    call bconstrbra(iesfree, nnozeroj, jbraj, nozeroj, jasmat        &
                         &, ipj*nelorbjh, ddw)
                end if
            end if
        end if
        if (iessz .and. rmaxinv .ne. 0.d0) then
            do k = 1, iesfreesz
                ddwsz(k) = alphavar(k)
            end do
            ! update jasmatsz_c
            if (contractionj .ne. 0) then
                call bconstrbra(iesfreesz, nnozeroj_c, jbrajsz, nozeroj_c       &
                     &, jasmatsz_c, nelorbj_c, ddwsz)
                call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
                     &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)
            else
                call bconstrbra(iesfreesz, nnozeroj, jbrajsz, nozeroj           &
                     &, jasmatsz, nelorbjh, ddwsz)
            end if
        end if
        if (iessw .gt. 0 .and. rmax .ne. 0.d0 .and. .not. detc_proj) then
            indc = iesinv + iesm + iesd + iesfree
            do k = indc + 1, indc + iessw
                dsw(k - indc) = alphavar(k)
            end do
            if (contraction .eq. 0) then
                ! From real to effective
                if (rank .eq. 0) write (6, *) ' Passi qui VII real-eff '
                if (allowed_averagek) call attach_phase2det(.false., detmat)
                call bconstraint(iessw, detmat, ipf*nelorbh, nnozero&
                     &, nozero, psip, dsw, 1, jbradet, symmagp, .false.)
                !  Back to real
                if (rank .eq. 0) write (6, *) ' Passi qui V eff-real '
                if (allowed_averagek) call attach_phase2det(.true., detmat)
            else
                ! From real to effective
                if (rank .eq. 0) write (6, *) ' Passi qui VIII real-eff '
                if (allowed_averagek) call attach_phase2det(.false., detmat_c)
                call bconstraint(iessw, detmat_c, nelorb_c, nnozero_c&
                     &, nozero_c, psip, dsw, 1, jbradet, symmagp, .false.)
                !  Back to real
                if (rank .eq. 0) write (6, *) ' Passi qui VI eff-real '
                if (allowed_averagek) call attach_phase2det(.true., detmat_c)
            end if
        end if
    end subroutine project_rmax

    subroutine project_alphavar
        implicit none
        integer i, j, k, indc
!$omp parallel do shared(kp0,rpar,alphavar,ncg_adr,v_adr,reducel) private(i,j)
        do j = 1, kp0
            if (rpar(j) .ne. 0.d0) then
                alphavar(j) = 0.d0
                do i = 1, ncg_adr
                    alphavar(j) = alphavar(j) + v_adr(i)*reducel(i, j)
                end do
            end if
        end do
!$omp end parallel do
#ifdef PARALLEL
#ifdef UNREL_DIAG
        call bcast_real(alphavar, kp0, 0, commopt_mpi)
#endif
#endif
        if (npar .gt. 0 .or. rmaxj .ne. 0.d0) then
            indc = iesinv + iesm + iesd
            do k = indc + 1, indc + iesfree
                ddw(k - indc) = alphavar(k)
            end do
            ! update jasmat_c
            if (contractionj .ne. 0) then
                call bconstrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c, jasmat_c  &
                     &, ipj*nelorbj_c, ddw)
!   if(ipj.eq.2) then
!      call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!   else
!      call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh&
!           &, nelorbj_c, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!   endif
            else
                if (yes_sparse) then
                    call bconstrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder, jasmat&
                         &, ipj*nelorbjh, ddw)
                else
                    call bconstrbra(iesfree, nnozeroj, jbraj, nozeroj, jasmat&
                         &, ipj*nelorbjh, ddw)
                end if
            end if
        end if
        if ((nparinv .gt. 0 .or. rmaxinv .ne. 0.d0) .and. iessz) then
            do k = 1, iesfreesz
                ddwsz(k) = alphavar(k)
            end do
            ! update jasmatsz_c
            if (contractionj .ne. 0) then
                call bconstrbra(iesfreesz, nnozeroj_c, jbrajsz, nozeroj_c       &
                     &, jasmatsz_c, nelorbj_c, ddwsz)
                call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
                     &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)
            else
                call bconstrbra(iesfreesz, nnozeroj, jbrajsz, nozeroj           &
                     &, jasmatsz, nelorbjh, ddwsz)
            end if
        end if
        if (iessw .gt. 0 .and. (nparsw .gt. 0 .or. rmax .ne. 0.d0) .and. .not. detc_proj) then
            indc = iesinv + iesm + iesd + iesfree
            do k = indc + 1, indc + iessw
                !          dsw(k-indc)=alphavar(k)+tpar*alphab(k)
                !          psip(k-indc)=tpar*alphab(k)    ! correction in this case
                dsw(k - indc) = alphavar(k)
            end do
            if (contraction .eq. 0) then
                !    From real to effective
                if (rank .eq. 0) write (6, *) ' Passi qui IX real-eff '
                if (allowed_averagek) call attach_phase2det(.false., detmat)
                call bconstraint(iessw, detmat, ipf*nelorbh, nnozero&
                     &, nozero, psip, dsw, 1, jbradet, symmagp, .false.)
                !    Back to real
                if (rank .eq. 0) write (6, *) ' Passi qui XIII eff-real '
                if (allowed_averagek) call attach_phase2det(.true., detmat)
            else
                !    From real to effective
                if (rank .eq. 0) write (6, *) ' Passi qui XI real-eff '
                if (allowed_averagek) call attach_phase2det(.false., detmat_c)
                call bconstraint(iessw, detmat_c, nelorb_c, nnozero_c&
                     &, nozero_c, psip, dsw, 1, jbradet, symmagp, .false.)
                !    Back to real
                if (rank .eq. 0) write (6, *) ' Passi qui XII eff-real '
                if (allowed_averagek) call attach_phase2det(.true., detmat_c)
            end if
        end if
    end subroutine project_alphavar

    subroutine update_ionpos
        implicit none
        real*8 rc(3), r0
        real*8, external :: jastrow_ei, cond_find
        if (yesupdate_ion) then
            if (ieskint .ne. 0) rion_write = rion
            if (cellderiv) then
                rs_write = rs
                celldm_write(1:3) = celldm(1:3)
                cellscale_write(1:3) = cellscale(1:3)
            end if
            indc = iesfreesz + iesm + iesd + iesfree + iessw + iesup
            if (iesking .ne. 0) then
                do k = indc + 1, indc + iesking
                    dekg(k - indc) = tpar*alphab(k)
                    !          alphavar(k)=dekg(k-indc)
                end do
                indc = indc + iesking
            end if

            if (ieskin .ne. 0) then
                do k = indc + 1, indc + ieskin
                    dek(k - indc) = tparf*alphab(k)
                    !         alphavar(k)=dek(k-indc)
                end do
                indc = indc + ieskin
            end if

            if (ieskint .gt. ieskinr_pos .and. cellderiv) then
                if (ieskint .eq. ieskinr_pos + 2) then
                    cellscale(2) = cellscale(2) + dek(ieskinr_pos + 1 - iesking)
                    cellscale(3) = cellscale(3) + dek(ieskinr_pos + 2 - iesking)
                    cellscale(1) = omega/cellscale(2)/cellscale(3)
                    if (rank .eq. 0) write (6, *) ' New PBC a,b (a=V/cb) =', cellscale(1:3)
                    scalecell(1) = cellscale(1)/celldm(1)
                    scalecell(2) = cellscale(2)/celldm(1)/celldm(2)
                    scalecell(3) = cellscale(3)/celldm(1)/celldm(3)
                elseif (ieskint .eq. ieskinr_pos + 3) then
                    cellscale(1) = cellscale(1) + dek(ieskinr_pos + 1 - iesking)
                    cellscale(2) = cellscale(2) + dek(ieskinr_pos + 2 - iesking)
                    cellscale(3) = cellscale(3) + dek(ieskinr_pos + 3 - iesking)
                    if (eqcellab) then
                        cellscale(1) = (cellscale(1) + cellscale(2))/2.d0
                        cellscale(2) = cellscale(1)
                    elseif (eqcellbc) then
                        cellscale(2) = (cellscale(2) + cellscale(3))/2.d0
                        cellscale(3) = cellscale(2)
                    elseif (eqcellac) then
                        cellscale(1) = (cellscale(1) + cellscale(3))/2.d0
                        cellscale(3) = cellscale(1)
                    end if
                    if (rank .eq. 0) write (6, *) ' New PBC a,b,c =', cellscale(1:3)
                    scalecell(1) = cellscale(1)/celldm(1)
                    scalecell(2) = cellscale(2)/celldm(1)/celldm(2)
                    scalecell(3) = cellscale(3)/celldm(1)/celldm(3)
                elseif (ieskint .eq. ieskinr_pos + 1) then
                    cellscale(1) = cellscale(1) + dek(ieskinr_pos + 1 - iesking)
                    cellscale(2) = cellscale(1)*celldm(2)
                    cellscale(3) = cellscale(1)*celldm(3)
                    if (rank .eq. 0) write (6, *) ' New PBC a,b,c =', cellscale(1:3)
                    scalecell(1) = cellscale(1)/celldm(1)
                    scalecell(2) = cellscale(2)/celldm(1)/celldm(2)
                    scalecell(3) = cellscale(3)/celldm(1)/celldm(3)
                end if

                call upcellkel(nion, nel, indt, scalecell, kel(1, indkj), rion)

                celldm(1) = cellscale(1)
                celldm(2:3) = cellscale(2:3)/cellscale(1)
                givens2r = .false.
                call InitCell(nion, nel, yes_complex)

                if (scalermax) then
                    if (fixa .and. fixb) then
                        cost = scalecell(3)
                    elseif (fixa) then
                        cost = scalecell(2)
                    else
                        cost = scalecell(1)
                    end if
                    if (cost .ne. 1.d0) then
                        if (rmax .gt. 1d-10) then
                            rmax = rmax*cost
                            if (rank .eq. 0) write (6, *) ' Warning scaled rmax =', rmax
                        end if
                        if (rmaxj .gt. 1d-10) then
                            rmaxj = rmaxj*cost
                            if (rank .eq. 0) write (6, *) ' Warning scaled rmaxj =', rmaxj
                        end if
                    end if
                end if

                kappa = kappar/lmin
                if (yes_tilted .and. ksq .ne. 0.5d0) kappa = kappa/cond_find(metric(1, 2), metric(1, 3), metric(2, 3))

                call InitEwald(nion, zetar, nel, nws)

                kmax = n_gvec
                kmax2 = 2*kmax

                rs = (omega/nel*3.d0/4.d0/pi)**(1.d0/3.d0)
                if (rank .eq. 0) write (6, *) ' Old/New rs =', rs_write, rs
                !         Update ax,ay,az
                ax = cellscale(1)/nx
                ay = cellscale(2)/ny
                az = cellscale(3)/nz

            end if ! endif  ieskin > ieskinr

            if (ieskint .gt. 0) then
                ! updating ionic positions
                kkf = 0
                kkfg = 0
                do kk = 1, ieskinr
                    ind = ion_table(kk)%mult
                    if (ind .gt. 0) then
                        ! ionic move
                        k_ion = abs(ion_table(kk)%ion(1))
                        if (atom_number(k_ion) .gt. 0) then
                            kkf = kkf + 1
                            do jj = 1, ind
                                i_ion = ion_table(kk)%comp(jj)
                                k_ion = abs(ion_table(kk)%ion(jj))
                                rion(i_ion, k_ion) = rion(i_ion, k_ion)                      &
                                     & + dek(kkf)*sign(1.d0, dble(ion_table(kk)%ion(jj)))
                            end do
                        else
                            kkfg = kkfg + 1
                            ! updating ghosts
                            do jj = 1, ind
                                i_ion = ion_table(kk)%comp(jj)
                                k_ion = abs(ion_table(kk)%ion(jj))
                                rion(i_ion, k_ion) = rion(i_ion, k_ion)                      &
                                     & + dekg(kkfg)*sign(1.d0, dble(ion_table(kk)%ion(jj)))
                            end do
                        end if
                    end if ! if it is allowed
                end do

                ! by E. Coccia (10/5/11): capping atom definition
                if (link_atom) then
                    call r_capping()
                end if

            end if ! end if ieskint.ne.0

#ifdef PARALLEL
#ifdef UNREL_DIAG
            if (yesquantum .and. nproc .gt. nbead) then
                ! in each row we have the same rion coordinates
                call bcast_real(rion, 3*nion, 0, commrep_mpi)
            end if
#endif
#endif
        end if ! yesupdate_ion
!  scale_one_body always updated to be sure and avoid overflows/underflows (minor cpu)
        if (add_onebody2det) then
            scale_one_body = 0.d0
            do jj = 1, nion
                if (iespbc) then
                    rc(:) = -rion(:, jj)
                    call CartesianToCrystal(rc, 1)
                    do kk = 1, 3
                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                    end do
                    r0 = norm_metric(rc, metric)
                else
                    rc(:) = (-rion(:, jj))*costz(jj)
                    r0 = dsqrt(sum(rc(:)**2))
                end if
                scale_one_body = scale_one_body - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
            end do
        end if

    end subroutine update_ionpos

    subroutine copywinv(nelorb, nel, indt4, winv, winvfn)
        implicit none
        integer nelorb, indt4, nel

        real*8 winv(nelorb, 0:indt4, nel), winvfn(nelorb, nel)

        winvfn(:, :) = winv(:, 0, :)

    end subroutine copywinv

!===================================
! wrappers for the main subroutines
!===================================
    subroutine computeb_global(js, yesfast, yesfastj&
     &, yeszagp, yesdodet, yeszj, yesforce, elocb, logpsib, membig, firstmol, nmolfn, detmat_c)
        implicit none
        logical yeszagp, yesdodet, yeszj, yesforce, membig
        integer j, js, indtj, indtjr, indtabj, indtabbj, indkj, indksj, indksij, indwwj&
           &, indwwjfn, indwwjj, indwupj, indwdoj, indaupj, indbar, indbarfn, indjbar, indjbarsz&
           &, yesfast, yesfastj, firstmol, nmolfn
        real*8 detmat_c(*), elocb(2), logpsib(2)

! calculation of the addresses depending on the current walker
        j = js - istm
        indtj = Lztab*(j - 1) + 1
        indtjr = Lztabr*(j - 1) + 1
        indtabj = Ltab*(j - 1) + 1
        indtabbj = Ltabb*(j - 1) + 1
        indkj = nel*(indt + 1)*(j - 1) + 1
        indksj = nel*(j - 1) + 1
        indksij = nelnion*(j - 1) + 1
        indwwj = nel2wt*(j - 1) + 1
        indwwjfn = nel2wtfn*(j - 1) + 1
        indwwjj = nel2wtj*(j - 1) + 1
        indwupj = nel2upt*(j - 1) + 1
        indwdoj = nel2dot*(j - 1) + 1
        indaupj = nel2up*(j - 1) + 1
        indbar = nel2bar*(j - 1) + 1
        indbarfn = nel2barfn*(j - 1) + 1
        indjbar = nel2jbar*(j - 1) + 1
        if (iessz) then
            indjbarsz = indjbar
        else
            indjbarsz = 1
        end if

        allocate (kelind(3, nel))
        call dcopy(nel3, kel(1, indkj), 1, kelind, 1)
!   ALlocate the minimum according to the options
        if (iessz) then
            allocate (winvjbarszb(nelorbjh*nel))
        else
            allocate (winvjbarszb(1))
        end if
        winvjbarszb = 0.d0
        allocate (kelb(3, nel*(indt + 1)), kelindb(3, nel), rionb(3, nion)&
           &, tabpipb(nel*(indt + 4)), winvupb(ipc*nelup*(indt + 4))&
           &, winvdob(neldomax*(indt + 4))&
           &, ainvb(ipc*nelup_mat*nelup_mat), ainvupbb(ipc*nelup*nelorbh), ainvdobb(neldomax*nelorbh)&
           &, psipb(iscramaxb), distb(nel*nion), rb(nion*(indt + 1)), rmub(3*nion*(indt + 1))&
           &, iond_cartb(3*nion*nion))
        allocate (vjb(size(vj)), vjurb(size(vjur)), duprb(size(dupr)))
        allocate (winvb(ipc*nelorb*nel), winvjb(nelorbjmax*nel))

        kelb = 0.d0
        kelindb = 0.d0
        rionb = 0.d0
        tabpipb = 0.d0
        winvupb = 0.d0
        winvdob = 0.d0
        ainvb = 0.d0
        ainvupbb = 0.d0
        ainvdobb = 0.d0
        psipb = 0.d0
        distb = 0.d0
        rb = 0.d0
        rmub = 0.d0
        iond_cartb = 0.d0
        vjb = 0.d0
        vjurb = 0.d0
        duprb = 0.d0
        winvb = 0.d0
        winvjb = 0.d0

        if (yesfast .eq. 0) then
            allocate (detmatb(ipc*ipf*nelorbh*nelcol))
            allocate (detmat_cb(1), mu_cb(1))
        else
            allocate (detmatb(ipc*ipf*nelorbh))
            allocate (detmat_cb(ipc*nelorb_c*nelcol_c), mu_cb(ipc*ipf*nelorbh*nelorb_c))
        end if
        if (npar_eagp .gt. 0) then
            allocate (eagp_pfaffb(ndiff*ipc, ndiff))
        else
            allocate (eagp_pfaffb(1, 1))
        end if
        eagp_pfaffb = 0.d0
        detmatb = 0.d0
        detmat_cb = 0.d0
        mu_cb = 0.d0

        if (nelorbjh .gt. 0) then
            if (yesfastj .eq. 0) then
                allocate (jasmatb(size(jasmat)))
                if (iessz) then
                    allocate (jasmatszb(nelorbjh*nelorbjh))
                else
                    allocate (jasmatszb(1))
                end if
                allocate (muj_cb(1), jasmat_cb(1), jasmatsz_cb(1))
            else
                allocate (jasmatb(1), jasmatszb(1))
                allocate (jasmat_cb(ipj*ipj*nelorbj_c*nelorbj_c), muj_cb(nelorbjh*nelorbj_c))
                if (iessz) then
                    allocate (jasmatsz_cb(nelorbj_c*nelorbj_c))
                else
                    allocate (jasmatsz_cb(1))
                end if
            end if
        else
            allocate (muj_cb(1), jasmat_cb(1), jasmatsz_cb(1))
            allocate (jasmatszb(1), jasmatb(1))
        end if

        jasmatb = 0.d0
        jasmat_cb = 0.d0
        jasmatszb = 0.d0
        jasmatsz_cb = 0.d0
        muj_cb = 0.d0

        allocate (winvbarb(ipf*ipc*nelorbh*nel_mat), winvjbarb(ipj*nelorbjh*nel + 1)&
           &, prefactorb((indt - istart + 1)*nel), wpseudob(2*lmax)&
           &, legendreb((lmax - 1)*nintpseudo), rmucosb(3*nion*(indt + 1))&
           &, rmusinb(3*nion*(indt + 1)), tmub(max(nel*indt, 1)), ivicb(3*indtmax*nel)&
           &, tabpipsav(nel*(indt + 4)))

        winvbarb = 0.d0
        winvjbarb = 0.d0
        prefactorb = 0.d0
        wpseudob = 0.d0
        legendreb = 0.d0
        rmucosb = 0.d0
        rmusinb = 0.d0
        tmub = 0.d0
        ivicb = 0.d0
        tabpipsav = 0.d0
! eloc and logpsi have been already defined by the direct algorithm
! and in any case are not  used by the reverse algorithm
!eloc(1) = enert(1, j)
!logpsi(1) = psiln(j)
!if(ipc.eq.2) then
! eloc(2) = enert(2, j)
! logpsi(2) = psisn(j)
!endif
#ifdef _OFFLOAD
!$omp target data map(from:jasmatb,muj_cb,jasmat_cb,detmatb,detmat_cb,mu_cb&
!$omp& ,eagp_pfaffb) map(alloc:psipb,winvb,winvjb,ainvb,winvbarb,winvjbarb&
!$omp& ,ainvupbb,ainvdobb)
#endif
        if (membig .or. (indt4 .eq. 0 .and. indt4j .eq. 0)) then

            call COMPUTE_ELOC_LOGPSI_B(indt, nelorb, nelup, neldo&
                 &, tabpip(indtabj), tabpipb, tabpipsav, kelind, kelindb, kel(1, indkj), kelb, &
                 & winv(indwwj), winvb, indt4, winvup(indwupj), winvupb, winvdo(indwdoj)&
                 &, winvdob, ainv(indaupj), ainvb, ainvup, ainvupbb, ainvdo, ainvdobb&
                 &, psip, psipb, ipsip, psisn(j), iesdr, vj, vjb, size(vj), dupr, duprb&
                 &, size(dupr), zetar, rion, rionb, dist(indksij), distb, ioccup, ioccdo, ioptorb&
                 &, nshell, nshelldo, ivic(1, 1, indksj), ivicb, vpot(j), tmu(indtabbj), tmub, nion, r, rb&
                 &, rmu, rmub, kion, iond, iond_cartb, winvj(indwwjj), winvjb, indt4j, ioccj, kionj, vjur, vjurb, size(vjur)&
                 &, nelorbj, ioptorbj, nshellj, winvbar(indbar), winvbarb, detmat, detmatb, eagp_pfaffb&
                 &, winvjbar(indjbar), winvjbarb, winvjbarsz(indjbarsz), winvjbarszb, jasmat, jasmatb&
                 &, muj_c, muj_cb, jasmat_c, jasmat_cb, jasmatsz, jasmatszb, jasmatsz_c, jasmatsz_cb, nelorbj_c, yesfastj&
                 &, iessz, cnorm, iflagerr, npsa, lmax, nintpseudo, prefactor, prefactorb, rcutoff, parshell&
                 &, nparpshell, kindion, pshell, wpseudo, wpseudob&
                 &, legendre, legendreb, versor, wintpseudo, jpseudo, pseudolocal(indksj), istart&
                 &, costz, costz3, angle(1, indksj), indtm(1, j), lbox, rmucos, rmucosb, rmusin, rmusinb&
                 &, kappa, vpotreg(1, indksj), cutreg, psidetln(j), j, nelorbh, nelorbjh, niesd, iond_cart&
                 &, mu_c, mu_cb, detmat_c, detmat_cb, nelorb_c, firstmol, nmolfn, yesfast, yeszagp, yesdodet&
                 &, yeszj, yesforce, eloc, elocb, logpsi&
                 &, logpsib, nelorbjmax, neldomax, indtmax, nshelljmax, cellelb, sr2elb, iflagnorm&
                 &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab&
                 &, adr_nion, ind_nion, adrj_nion, indj_nion, .true.)

        else

#ifdef  _OFFLOAD
            ! The update below is likely not necessary because after upscratch_global
            ! winv and winvj are both  in CPU & GPU
!$omp target update from(winv(indwwj:indwwj+nel2wt-1)) if(yes_ontarget)
!$omp target update from(winvj(indwwjj:indwwjj+nel2wtj-1)) if(yes_ontarget.and.nel2wtj.gt.0)
#endif
            allocate (winvfnn(ipc*nelorb*nel))
            winvfnn = 0.d0
            call copywinv(ipc*nelorb, nel, indt4, winv(indwwj), winvfnn)
            allocate (winvjfn(nelorbjmax*nel))
            winvjfn = 0.d0
            if (nelorbj .ne. 0) then
                call copywinv(nelorbj, nel, indt4j, winvj(indwwjj), winvjfn)
            end if
#ifdef _OFFLOAD
!$omp target data map(to:winvfnn,winvjfn)
#endif
            call COMPUTE_ELOC_LOGPSI_B(indt, nelorb, nelup, neldo&
                 &, tabpip(indtabj), tabpipb, tabpipsav, kelind, kelindb, kel(1, indkj), kelb, &
                 & winvfnn, winvb, 0, winvup(indwupj), winvupb, winvdo(indwdoj)&
                 &, winvdob, ainv(indaupj), ainvb, ainvup, ainvupbb, ainvdo, ainvdobb&
                 &, psip, psipb, ipsip, psisn(j), iesdr, vj, vjb, size(vj), dupr, duprb&
                 &, size(dupr), zetar, rion, rionb, dist(indksij), distb, ioccup, ioccdo, ioptorb&
                 &, nshell, nshelldo, ivic(1, 1, indksj), ivicb, vpot(j), tmu(indtabbj), tmub, nion, r, rb&
                 &, rmu, rmub, kion, iond, iond_cartb, winvjfn, winvjb, 0, ioccj, kionj, vjur, vjurb, size(vjur)&
                 &, nelorbj, ioptorbj, nshellj, winvbar(indbar), winvbarb, detmat, detmatb, eagp_pfaffb&
                 &, winvjbar(indjbar), winvjbarb, winvjbarsz(indjbarsz), winvjbarszb, jasmat, jasmatb&
                 &, muj_c, muj_cb, jasmat_c, jasmat_cb, jasmatsz, jasmatszb, jasmatsz_c, jasmatsz_cb, nelorbj_c, yesfastj&
                 &, iessz, cnorm, iflagerr, npsa, lmax, nintpseudo, prefactor, prefactorb, rcutoff, parshell&
                 &, nparpshell, kindion, pshell, wpseudo, wpseudob&
                 &, legendre, legendreb, versor, wintpseudo, jpseudo, pseudolocal(indksj), istart&
                 &, costz, costz3, angle(1, indksj), indtm(1, j), lbox, rmucos, rmucosb, rmusin, rmusinb&
                 &, kappa, vpotreg(1, indksj), cutreg, psidetln(j), j, nelorbh, nelorbjh, niesd, iond_cart&
                 &, mu_c, mu_cb, detmat_c, detmat_cb, nelorb_c, firstmol, nmolfn, yesfast, yeszagp, yesdodet, yeszj&
                 &, yesforce, eloc, elocb, logpsi&
                 &, logpsib, nelorbjmax, neldomax, indtmax, nshelljmax, cellelb, sr2elb, iflagnorm&
                 &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab&
                 &, adr_nion, ind_nion, adrj_nion, indj_nion, .true.)
#ifdef _OFFLOAD
!$omp end target data
#endif
            deallocate (winvfnn, winvjfn)
        end if

#ifdef _OFFLOAD
!$omp end target data
#endif
! Not possible to vanish otherwise the derivatives will not work.
!!       if(yes_hermite.and.symmagp.and.ipc.eq.2.and.yesfast.ne.0) then
!!      constrained to be real to allow an hermitian matrix
!!      Only the atomic orbitals. The molecular ones are free to be complex as detmat_c.
!       do j=2,2*nelorbh*nelorb_at,2
!       mu_cb(j)=0.d0
!       enddo
!       endif
        return

    end subroutine computeb_global
    subroutine deallocate_computeb
        implicit none
        deallocate (winvjbarszb)
        deallocate (kelb, kelind, kelindb, rionb&
           &, tabpipb, winvupb, winvdob, ainvb, ainvupbb, ainvdobb&
           &, psipb, distb, rb, rmub, iond_cartb, winvb, winvjb)
        deallocate (vjb, vjurb, duprb)
        deallocate (detmatb, detmat_cb, mu_cb)

        deallocate (eagp_pfaffb)
        deallocate (jasmatb, jasmatszb)
        deallocate (muj_cb, jasmat_cb, jasmatsz_cb)

        deallocate (winvbarb, winvjbarb&
           &, prefactorb, wpseudob, legendreb, rmucosb&
           &, rmusinb, tmub, ivicb, tabpipsav)
    end subroutine deallocate_computeb

    subroutine upscratch_global(js, pseudologicr, iesrandomlr)

        use allio, only: pseudologic, iesrandoml
        implicit none
        logical pseudologicr, iesrandomlr, pseudologic_sav, iesrandoml_sav
        integer j, js, k, l, indtj, indtjr, indtabj, indtabbj, indkj, indksj, indksij, indwwj&
           &, indwwjj, indwupj, indwdoj, indaupj, indbar, indjbar, indjbarsz&
           &, indbarfn, indwwjfn
!   Not destroy the info in the module, restored at the end
        pseudologic_sav = pseudologic
        iesrandoml_sav = iesrandoml
        pseudologic = pseudologicr
        iesrandoml = iesrandomlr

! calculation of the addresses depending on the current walker
        j = js - istm
        indtj = Lztab*(j - 1) + 1
        indtjr = Lztabr*(j - 1) + 1
        indtabj = Ltab*(j - 1) + 1
        indtabbj = Ltabb*(j - 1) + 1
        indkj = nrnel*(j - 1) + 1
        indksj = nel*(j - 1) + 1
        indksij = nelnion*(j - 1) + 1
        indwwj = nel2wt*(j - 1) + 1
        indwwjfn = nel2wtfn*(j - 1) + 1
        indwwjj = nel2wtj*(j - 1) + 1
        indwupj = nel2upt*(j - 1) + 1
        indwdoj = nel2dot*(j - 1) + 1
        indaupj = nel2up*(j - 1) + 1
        indbar = nel2bar*(j - 1) + 1
        indbarfn = nel2barfn*(j - 1) + 1
        indjbar = nel2jbar*(j - 1) + 1
        if (iessz) then
            indjbarsz = indjbar
        else
            indjbarsz = 1
        end if
!   if(.not.allocated(kelind)) allocate(kelind(3,nel))
!   call dcopy(nel3,kel(1,indkj),1,kelind,1)
        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
           &, tabpip(indtabj), kel(1, indkj), kel(1, indkj), winv(indwwj), winvup(indwupj), winvdo(indwdoj)&
           &, ainv(indaupj), ainvup, ainvdo, psip(nel3 + 1), ipsip, wconfn(js), psisn(j), iesdr&
           &, vj, dupr, zetar, rion, dist(indksij), ioccup, ioccdo, ioptorb&
           &, nshell, nshelldo, ivic(1, 1, indksj), alat, plat, vpot(j), tmu(indtabbj)&
           &, nion, r, rmu, kion, iond, winvj(indwwjj), ioccj, kionj, vjur, nelorbj   &
           &, ioptorbj, nshellj, winvbar(indbar), detmat, winvjbar(indjbar), winvjbarsz(indjbarsz), jasmat&
           &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax&
           &, nintpseudo, prefactor, rcutoff, parshell                            &
           &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
           &, jpseudo, pseudolocal(indksj), istart, costz, costz3         &
           &, angle(1, indksj), indtm(1, j), LBox, rmucos, rmusin, kappa, vpotreg(1, indksj), cutreg           &
           &, psidetln(j), j, nelorbh, nelorbjh                         &
           &, niesd, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast, eloc, logpsi&
           &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
           &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

!   deallocate(kelind)

        pseudologic = pseudologic_sav
        iesrandoml = iesrandoml_sav

        return

    end subroutine upscratch_global

    subroutine uptabtot_global(js, pseudologic)
        implicit none
        logical pseudologic
        integer j, js, indtj, indtjr, indtabj, indtabbj, indkj, indksj, indksij, indwwj&
           &, indwwjfn, indwwjj, indwupj, indwdoj, indaupj, indbar, indbarfn, indjbar, indjbarsz, jpot

! calculation of the addresses depending on the current walker
        j = js - istm
        indtj = Lztab*(j - 1) + 1
        indtjr = Lztabr*(j - 1) + 1
        indtabj = Ltab*(j - 1) + 1
        indtabbj = Ltabb*(j - 1) + 1
        indkj = nel*(indt + 1)*(j - 1) + 1
        indksj = nel*(j - 1) + 1
        indksij = nelnion*(j - 1) + 1
        indwwj = nel2wt*(j - 1) + 1
        indwwjfn = nel2wtfn*(j - 1) + 1
        indwwjj = nel2wtj*(j - 1) + 1
        indwupj = nel2upt*(j - 1) + 1
        indwdoj = nel2dot*(j - 1) + 1
        indaupj = nel2up*(j - 1) + 1
        indbar = nel2bar*(j - 1) + 1
        indbarfn = nel2barfn*(j - 1) + 1
        indjbar = nel2jbar*(j - 1) + 1
        if (iessz) then
            indjbarsz = indjbar
        else
            indjbarsz = 1
        end if
!
        if (itest .ne. 2) then
            jpot = j
        else
            jpot = 1
        end if

        call uptabtot(nelup, neldo, nelorb, nelorbh, iout &
                      , kel(1, indkj), winv(indwwj), winvup(indwupj) &
                      , winvdo(indwdoj), ainv(indaupj), psiln(j), psisn(j) &
                      , epst, psip(nelorbpp), rcart, dist(indksij), dists, ainvs, winvs, psinew &
                      , indt, indt4, indt4j, indvic, tabpip(indtabj), alat, ivic(1, 1, indksj) &
                      , plat, rion, vpot(j) &
                      , iesdr, vj, zetar, dupr, nshellr, nshellr, ioptorb, ioccup, ioccdo &
                      , itestrfn, tmu(indtabbj), nion, r, rmu, kion, winvj(indwwjj), ioccj &
                      , kionj, vjur, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj &
                      , winvbar(indbar), detmat, winvjbar(indjbar), psip &
                      , winvjbarsz(indjbarsz), psip(nelorbp), jasmat, jasmatsz, iflagnorm &
                      , cnorm, iflagerr, ratior, npsa, lmax, nintpseudo, prefactor, rcutoff &
                      , parshell, nparpshell, kindion, pshell, wpseudo, legendre, versor &
                      , wintpseudo, jpseudo, pseudolocal(indksj), tcost, istart &
                      , costz, costz3, pseudologic, angle(1, indksj), iessz, keln, kappa &
                      , LBox, rmucos, rmusin, j, n_body_on &
                      , jastrowall_ee(1, 1, 0, j), jasnew_ee, jastrowall_ei(1, 1, j) &
                      , jasnew_ei, timepip, timewf, niesd, versoralat, iesrandoma &
                      , nshell, indtm(1, j), psidetln(j), mu_c, projm, nelorb_c, firstmol, nmolfn &
                      , yesfast, vpotreg(1, indksj), cutreg &
                      , indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab &
                      , detmat_c, winvbarfn(indbarfn), winvfn(indwwjfn), contraction, vpotsav_ee(1, jpot))

    end subroutine uptabtot_global

    subroutine read_alphavar
        implicit none
!       put the initial parameter wavefunction in alphavar

        if (np .gt. 0) then

            if (iopt .eq. 1 .or. .not. yesquantum .or. (yesquantum .and. yesread10)) then

                dek(1:ieskinr_pos) = 0.d0

                indc = 0
                if (iesinv .ne. 0) then
                    do k = indc + 1, indc + iesinv
                        alphavar(k) = ddwsz(k - indc)
                    end do
                end if
                indc = indc + iesinv
                if (iesm .ne. 0) then
                    do k = indc + 1, indc + iesm
                        alphavar(k) = vju(k - indc)
                    end do
                end if
                indc = indc + iesm

                if (iesd .ne. 0) then
                    do k = indc + 1, indc + iesd
                        alphavar(k) = vj(k - indc)
                    end do
                end if
                indc = indc + iesd

                if (iesfree .ne. 0) then
                    do k = indc + 1, indc + iesfree
                        alphavar(k) = ddw(k - indc)
                    end do
                end if
                indc = indc + iesfree

                if (iessw .ne. 0) then
                    do k = indc + 1, indc + iessw
                        alphavar(k) = dsw(k - indc)
                    end do
                end if
                indc = indc + iessw

                if (iesup .ne. 0) then
                    do k = indc + 1, indc + iesup
                        alphavar(k) = dup(k - indc)
                    end do
                end if

                indc = indc + iesup

                if (iesking .ne. 0) then
                    do k = indc + 1, indc + iesking
                        alphavar(k) = dekg(k - indc)
                    end do
                end if
                indc = indc + iesking

                if (ieskin .ne. 0) then
                    ! The ion mapping coordinates
                    do k = indc + 1, indc + ieskin
                        alphavar(k) = dek(k - indc)
                    end do
                end if

                if (iespbc .and. itestr .eq. -5) then
                    !       Restoring consistency with fort.10
                    if (ieskint .eq. ieskinr_pos + 2) then
                        dek(ieskinr_pos + 1 - iesking) = cellscale(2)
                        dek(ieskinr_pos + 2 - iesking) = cellscale(3)
                    elseif (ieskint .eq. ieskinr_pos + 1) then
                        dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                    elseif (ieskint .eq. ieskinr_pos + 3) then
                        dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                        dek(ieskinr_pos + 2 - iesking) = cellscale(2)
                        dek(ieskinr_pos + 3 - iesking) = cellscale(3)
                    end if
                end if

            else
                !    Read parameters from alphavar
                indc = 0
                if (iesinv .ne. 0) then
                    do k = indc + 1, indc + iesinv
                        ddwsz(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesinv
                if (iesm .ne. 0) then
                    do k = indc + 1, indc + iesm
                        vju(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesm

                if (iesd .ne. 0) then
                    do k = indc + 1, indc + iesd
                        vj(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesd

                if (iesfree .ne. 0) then
                    do k = indc + 1, indc + iesfree
                        ddw(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesfree

                if (iessw .ne. 0) then
                    do k = indc + 1, indc + iessw
                        dsw(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iessw

                if (iesup .ne. 0) then
                    do k = indc + 1, indc + iesup
                        dup(k - indc) = alphavar(k)
                    end do
                end if

                indc = indc + iesup

                if (iesking .ne. 0) then
                    do k = indc + 1, indc + iesking
                        dekg(k - indc) = 0.d0
                        alphavar(k) = 0.d0
                    end do
                end if
                indc = indc + iesking

                if (iespbc .and. itestr .eq. -5) then
                    !       Restoring consistency with fort.10
                    if (ieskint .eq. ieskinr_pos + 2) then
                        cellscale(2) = dek(ieskinr_pos + 1 - iesking)
                        cellscale(3) = dek(ieskinr_pos + 2 - iesking)
                        cellscale(1) = omega/cellscale(2)/cellscale(3)
                    elseif (ieskint .eq. ieskinr_pos + 1) then
                        cellscale(1) = dek(ieskinr_pos + 1 - iesking)
                        cellscale(2:3) = 1.d0
                    elseif (ieskint .eq. ieskinr_pos + 3) then
                        cellscale(1) = dek(ieskinr_pos + 1 - iesking)
                        cellscale(2) = dek(ieskinr_pos + 2 - iesking)
                        cellscale(3) = dek(ieskinr_pos + 3 - iesking)
                    end if
                end if

                if (ieskin .ne. 0) then
                    ! The ion mapping coordinates
                    do k = indc + 1, indc + ieskin
                        if (k - indc .le. ieskinr_pos - iesking) then
                            dek(k - indc) = 0.d0
                            alphavar(k) = 0.d0
                        else
                            dek(k - indc) = alphavar(k)
                        end if
                    end do
                end if

            end if ! not yesquantum or iopt=1

        end if ! np.gt.0

    end subroutine read_alphavar

end program
