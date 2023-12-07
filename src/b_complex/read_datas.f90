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

subroutine read_datasmin
    use allio
    use Thomas_Fermi_model !!!! new !!!! added by K.Nakano 11/09/2019
    use dielectric
    use qmckl
    use cartesian2spherical
    implicit none
    integer :: i, ii, j, ind
    real(8) :: num_ele_core_r_c, r_c, beta_for_r_c, kappa_for_r_c !!!! new !!!! added by K.Nakano 11/09/2019
    real(8), external :: dlamch
    logical noreadnpbra, def_nscra
    real*8 derfc, err_derfc
    integer(kind=qmckl_exit_code) :: rc
    integer(kind=8) :: qmckl_shell_num, qmckl_prim_num, qmckl_index(10000), npoints_qmckl, ao_num
    integer(kind=4) :: qmckl_index4(10000)
    real*8 :: electrons(3, 2)
    real*8, allocatable :: ao_value_qmckl(:)
#ifdef _OFFLOAD
    integer*8, parameter :: max_sparse = 40
#else
    integer*8, parameter :: max_sparse = 10
#endif
#ifdef PARALLEL
    include 'mpif.h'
#endif

    !by E. Coccia (9/5/11): namelist for link atoms in case of QMC/MM
    !default values
    ! dimension for calpha (5 capping atoms as maximum)
    maxcap = 5
    calpha(:) = 0.d0
    noreadnpbra = .true.

    if (rank .eq. 0) then
        !          read(5,*) etry,rsignr,nintpsa,npsa,npsamax
        !          read(5,*) ngen,nscra,nbra,npow,ncore,nweight,nfat,iboot
        !          read(5,*) iseedr,iopt,iread,nw,beta,ncg
        !          read(5,*) itestr4,tbra,gamma,nbinr,npbra,eta,epsi,epsdgel
        !    1,tstep,tpar,hopfraction,epscut,epstl,testderiv,kl
        !          read(5,*) ieser,iesinv,iesm,iesd,isfix,np3
        !    1,iesfree,iessw,iesup,ieskin
        !          read(5,*) parr,parcutr,parcutmin,minz,maxz,plat(1),plat(2)
        !    1,plat(3),alat2,alat,parcutpar

        !          defining default values for simulation
        ! do not read the unused section
        membig = .true. ! if .true. assuming large memory is given. cpu optimization
        ! when membig=.false. GPU memory is minimized but CPU
        ! memory is assumed large if membigcpu=.true.
        membigcpu = .true. ! Always faster algorithm
        freqcheck = -1
        developer = 0
        kSq = 0.d0
        kappar = 0.d0
        neigh = -1
        nscra = 0
        nbra = 0
        iseedr = 236413883
        nw = nproc
        yesfast = -1
        ip_reshuff = -1
        ! vmc
        itestr4 = 2
        ! startup
        ngen = 100
        ! begin the simulation
        iopt = 1
        iflagerrall = iflagerr + iflagerrall
        compute_bands = .false.
        ! Default minimum block size for computational efficiency at expence of memory. A smaller value will require less memory.
        min_block = 100

        ! This section should always exist
        iflagerr = 1
        maxtime = -1.d0
        ! Parallel diag is typically useless when optimization molecular  orbitals  is done (contracted basis is small)
        nproc_diag = 0

        ! control io behavior
        io_level = 1
        disk_io = 'default'
        double_mesh = .false.
#ifdef _OFFLOAD
        max_target = 200000 ! it is usually worth
#else
        max_target = 0
#endif
        max_targetsr = -1
        dielectric_ratio = 1.d0
        dielectric_length = 1.d0
        case_diel = 0
        novec_loop1 = .true.
        yes_sparse = .false.
        if (yesdft) then
            yes_sparse_choose = .false.
        else
            yes_sparse_choose = .true.
        end if
#ifdef _OFFLOAD
        max_sparse_choice = 100
#else
        max_sparse_choice = 20
#endif

        read (5, nml=simulation, err=111)

        write (6, *) ' After  reading simulation  '
        iflagerr = 0
111     if (iflagerr .ne. 0) write (6, *) ' ERROR reading simulation '
        iflagerrall = iflagerr + iflagerrall
        !          calculation itestr
        if (maxtime .lt. 0.d0) then
            maxtime = 86000
            write (6, *) ' Default maxtime (s)= ', maxtime
            write (6, *) ' Warning the program will stop after <~24h, &
                 & otherwise change maxtime in simulation section '
        end if

#ifndef _OFFLOAD
        if (.not. membig .and. membigcpu) then
            write (6, *) ' Warning this option is only  for test'//&
                        &', otherwise put also mambigcpu=.false.'//&
                        &' to minimize RAM memory allocation !!! '
        end if
#endif

        select case (trim(disk_io))
#ifdef PARALLEL
        case ('mpiio')
            io_level = 2
#endif
        case ('default')
            io_level = 1
        case ('nocont')
            io_level = 0
        case DEFAULT
            write (6, *) " Warning: unrecognized disk_io, set as 'default' !"
            io_level = 1
        end select

        if (io_level .eq. 0) write (6, *) " Warning: you cannot continue this run !!! Avoid disk_io='nocont' otherwise!"
        if (io_level .eq. 0 .and. iopt .ne. 1) then
            write (6, *) " Warning: disk_io changed to 'default' in continuation !!!"
            io_level = 1
        end if

        if (itestr4 .ge. -10) then
            !      Normal minimization
            !         definition simulation variable
            !         itestr=-5 minimization/dynamic
            !         itest=1  FN/LRDMC
            !         itest=2  VMC
            itestr3 = itestr4
        elseif (itestr4 .ge. -20) then
            itestr3 = itestr4 + 10
        elseif (itestr4 .ge. -30) then
            itestr3 = itestr4 + 20
        end if
        if (freqcheck .eq. -1 .and. nproc .gt. 1) then
            if (ngen .gt. 1000) then
                freqcheck = 100
            elseif (ngen .gt. 100) then
                freqcheck = 10
            else
                freqcheck = 0
            end if

            write (6, *) ' Default value of check flag error ', freqcheck
        end if

        itestrr = itestr3 ! a number < |10|
        if (itestr3 .eq. -8) itestrr = -4
        if (itestr3 .eq. -9) itestrr = -5
        itestr = itestrr
        if (itestr .eq. -4) itestr = -5
        itest = abs(itestr)
        if (itestr .eq. -5) then
            if (itestr4 .ge. -10) then
                itest = 2
            else
                itest = 1
            end if
        end if
        if (abs(itestr) .eq. 6 .or. itestr .eq. -2 .or. itestr .eq. -3) itest = 1
        if (itestr .eq. -7) itest = 2

        ! defining default values for simulation
        nintpsa = 0
        npsamax = 0
        npsa = -1
        if (yesdft .and. double_mesh) then
            pseudorandom = .false.
        else
            pseudorandom = .true.
        end if
        write (6, *) ' Default value of pseudorandom =', pseudorandom
        iflagerr = 1
        pseudofile = "pseudo.dat"
        read (5, nml=pseudo, err=112)
        write (6, *) ' After reading pseudo '
        write (6, *) ' Pseudopotential file name : ', trim(pseudofile)
        iflagerr = 0
112     if (iflagerr .ne. 0) write (6, *) ' ERROR reading pseudo '
        iflagerrall = iflagerr + iflagerrall

        ! , you have to know if you are using pseudo or not (empty section).
        !          write(6,*) ' after pseudo  ',npsa,nintpsa,npsamax
        tstep = -1.d0
        ! no large hopping move
        hopfraction = 0.d0
        ! no regularization
        epscut = 0.00001d0 ! default regularization
        !    epstlrat=dsqrt(epsmach) ! It looks important even if epsvar>0, machine precision error under control.
        epstlrat = 0.d0 ! with theta_reg<1/2 it should not be important.
        theta_reg = 0.375d0
        epscuttype = -1000
        change_epscut = .true.
        change_tstep = .true.
        epstl = 0.d0
        shift = 0.d0
        alat2v = 0.d0
        epsvar = 0.0001d0
        true_wagner = -1
        cutweight = -1.d6
        typereg = 0
        npow = 0.d0
        if (itest .eq. 2) then
            iflagerr = 1
            read (5, nml=vmc, err=113)
            if (hopfraction .ne. 0.d0) write (6, *) ' Read value of hopfraction ', hopfraction
            if (epscut .ne. 0.d0) write (6, *) ' Read value of epscut ', epscut
            if (epscuttype .ne. -1000) write (6, *) ' Read value of epscuttype ', epscuttype
            write (6, *) ' After reading vmc '
            iflagerr = 0
113         if (iflagerr .ne. 0) write (6, *) ' ERROR reading vmc '
            iflagerrall = iflagerr + iflagerrall
        end if

        if (theta_reg .eq. 0.5d0 .and. epstlrat .eq. 0.d0) then
            epstlrat = dsqrt(epsmach) ! It looks important even if epsvar>0, machine precision error under control.
            write (6, *) ' Warning epstlrat should be > 0 in this case, chosen one=', epstlrat
        end if

        !          defining default values for lrdmc
        ! no idea
        etry = 0.d0
        ! no locality
        ! no idea
        tbra = 0.1d0
        gamma = 0.d0
        plat = 0.d0
        alat2 = 0.d0
        alat = 0.d0
        tstepfn = 0.d0
        Klrdmc = 0.d0
        optbra = 0
        parcutg = 0
        defparcutg = .false.
        novar = 0
        epscutdmc = 0.d0
        epstldmc = 0.d0
        changelambda = .false.
        lrdmc_der = .false.
        lrdmc_nonodes = .false.
        enforce_detailb = .false.
        iesrandoma = .true.
        weight_moroni = 1.0d0
        add_diff = .true.
        if (itest .eq. 1 .or. itestr4 .eq. -7) then
            plat = -1000.d0
            alat2 = -1.d0
            Klrdmc = -1000.d0
            novar = -1
            parcutg = -1000
            epstldmc = -1.d0
            rejweight = .true.
            cutreg = 0.d0
            epscuttype = 0
            epscut = 0.d0
            better_dmc = .true.
            gauge_fixing = .false.
            scalermax = .false.
            yesalfe = .false.
            changelambda = .true.
            noblocking = .false.
            nbra_cyrus = 0
            if (itestr .ne. -3) yesalfe = .true.
            safelrdmc = .false.
            if (itestr .eq. -3) rejweight = .false.
            zmin = 0.d0
            ! default value for rejweight
            ! (rejecting or rescaling the weights according to the acceptance/rejection step)
            ! in standard dmc or non local dmc with heat bath after all electron diffusion
            ! rejecting or rescaling the weights according to acceptance should be the best choice.
            ! in non local dmc with heat bath after single particle diffusion the best choice
            ! is not to reject the weights.
            yes_fastbranch = .false.
            l0_kousuke = 2.d0
            nw_max = -1
            iflagerr = 1
            read (5, nml=dmclrdmc, err=114)
            write (6, *) ' After reading dmclrdmc '
            iflagerr = 0
114         if (iflagerr .ne. 0) write (6, *) ' ERROR reading dmclrdmc '
            iflagerrall = iflagerr + iflagerrall
            if (nw_max .gt. 1) then
                if (mod(nw, nw_max) .ne. 0) then
                    write (6, *) ' ERROR nw_max has to be divisor of nw  !!! '
                    iflagerr = 1
                    iflagerrall = iflagerr + iflagerrall
                end if
            end if

            if (parcutg .eq. -1000) then
                if ((itestr .eq. -6 .or. itestr4 .eq. -7 .or. itestr4 .lt. -20)&
                     &.and. .not. (lrdmc_der .and. .not. lrdmc_nonodes)) then
                    ! To avoid errors, also equivalent to parcutg=1 for soft pseudo,
                    ! parcutg=1 is only slightly faster in this case.
                    parcutg = 1
                else
                    parcutg = 0
                end if
                defparcutg = .true.
                write (6, *) ' Default value of parcutg=', parcutg
            end if

            if (epscutdmc .ne. 0.d0 .and. epstldmc .eq. -1.d0) then
                epstldmc = 0.d0
                write (6, *) ' Default value for epstldmc (no cutoff) =', epstldmc
            elseif (epscutdmc .eq. 0.d0) then
                epstldmc = 0.d0
            end if

            if (epscutdmc .ne. 0.d0 .and. gamma .eq. 0.d0) then
                write (6, *) ' Warning change gamma>0 with epscutdmc>0 !!!'
            end if

            write (6, *) ' after dmc  '
        end if
        !          defining default values for optimization

        ! typical acceleration with Hessian
        tpar = 0.35d0
        ! nweight does not change
        ! no correcting factore
        nfat = 0
        iboot = 0
        nweight = 1
        ! no cutoff on norm
        nmore_force = -1
        epsi = 10000d0
        epstion = 0.d0
        ! typical value to control conj grad and ill con
        epsdgel = -1.d0
        ! it works with parametrizations and molecular
        tolcg = -1.d0 ! default = 1d-6
        kl = -7
        ! no dynamic
        idyn = 0
        ! zero is not allowed
        nbinr = 1
        ! no normal parameters optimization
        npbra = 0
        ! Only one hessian direction
        if (itestr3 .eq. -8 .or. itestr3 .eq. -5) then
            yeszj = .true.
            yeszagp = .true.
        else
            yeszj = .false.
            yeszagp = .false.
        end if
        ncg = 1
        ! no cutoff on Z
        !minz=-1000.d0
        minz = 1.d-5
        maxz = -1000.d0
        !minzj=-1000.d0
        minzj = 1.d-5
        maxzj = -1000.d0
        minjonetwobody = -1.d0
        srcomplex = .true.
        power = 1.d0
        parr = 0.0d0 ! cutoff on inverse for S x = f (SR step)
        parcute = 0.d0 ! no cutoff on Hessian for Nightingale linear method

        parcut = epsmach*100.d0 ! safe relative machine precision  at least 12 digts but usually 8
        write (6, *) ' Default value of parcut = ', parcut

        parcutmin = 0.d0 ! no cutoff on devmax collective parameters

        parcutpar = 0.d0 !  ==                 normal     parameters

        tion = -1.d0 ! acceleration ion (unit mass)

        eps_dyn5 = 0.d0 ! no cut off for maximum change in norm  in idyn=5

        tcell = 0.d0 ! acceleration cell

        oldscaling = .false.

        molopt = 0 ! no optimization molecular with projection if an
        symiesup = .false.
        onebodysz = .false. ! optimize also the two-body Sz
        twobodyoff = .false. ! optimize two-three-four body Jastrow
        iesdtwobodyoff = .false. ! do not optimize two- body Jastrow (first parameter in iesd)
        iesdonebodyoff = .false. ! do not optimize one- body Jastrow (last parameters in iesd)
        yescutdet = .true.
        yescutjas = .true.
        fixpar = .false.
        symmetrize_agp = .true.
        nbead = -1
        yesquantum = .false.
        yeswrite10 = .false.
        yesread10 = .false.
        yeswritebead = .false.
        yesperiodize = .false.
        signalnoise = .false.
        noopt_onebody = .false.
        change_parr = .false.
        !    Option for automatic optimization with change_parr=.true.
        parr_min = 1d-4 ! The minimum parr set for convergence
        parr_max = 1d-1 ! Do not increase again parr> parr_max
        delay_changeparr = 20 ! Each time  parr is changed remain there for at  least 20 iterations
        maxiter_changeparr = -1 ! If positive the run stops after maxiter  iterations with devmax < parr_min

        !    Option for equilibration
        equil_steps = -1 ! the default value is set later (after reading nweight)

        ! added by Andrea Tirelli: additional parameters for adjust_tpar subroutine

        change_tpar = .false.
        use_stable_tpar = .true.
        ! added by Andrea Tirelli: additional parameters for adjust_tpar subroutine
        tpar_increased = .false.
        stop_increasing_tpar = .false.
        inc_tpar_frequency = 100
        divide_tpar = 4.0d0
        multiply_tpar = -1.0d0
        tpar_unstble_stop = 3
        tpar_buffer_len = 40
        len_shorter_buffer = 10
        times_tpar_decreased = 0
        cut_sigma = 3.5d0
        n_sigmas_tpar = 5.0d0
        len_tpar_stable_list = 10

        k6gen = .false. ! Default old reliable
        max_ortho = -1 ! Default all possible processors are involved in the k6gen linear system
        prep = -1 ! Default do not use prep
        beta_learning = -1.d0
        eps_umrigar = -1.d0
        yes_adams = .false.
        tpar_max = 1.d0 ! maximum tpar in linear method with change_tpar=.true.
        yes_dgelscut = .true.
        if (itestr .eq. -5) then
            iflagerr = 1
            !          if(molopt.ne.0) then
            !          membig=.false.
            !          write(6,*) ' Warning default membig set to .false., if something &
            !         &strange put membig=.true. '
            !          endif
            read (5, nml=optimization, err=115)
            write (6, *) ' After reading optimization   '
            iflagerr = 0
115         if (iflagerr .ne. 0) write (6, *) ' ERROR reading optimization '
        end if
        if (oldscaling) then
            write (6, *) ' Warning old definition wrong scaling of time sqrt(2) smaller than standard'
        end if

        if (use_stable_tpar) inc_tpar_frequency = 50

        if (change_tpar) then

            allocate (energy_list(tpar_buffer_len))
            allocate (error_energy_list(tpar_buffer_len))
            allocate (tpar_buffer_filled(tpar_buffer_len))
            allocate (tpar_stable_list(len_tpar_stable_list))

            tpar_stable_list(:) = -1

            if (iopt .ne. 0) then
                tpar_buffer_filled(:) = .false.
            else
                tpar_buffer_filled(:) = .true.
            end if

        end if

        if (multiply_tpar .eq. -1.0d0) then
            multiply_tpar = exp(log(divide_tpar)/dble(inc_tpar_frequency)*4)
        end if

#if defined PARALLEL && defined __SCALAPACK
        if (k6gen .and. min_block .eq. 100) then
            min_block = 1000
            write (6, *) ' Default value of min_block =', min_block
        end if
#else
        if (k6gen) then
            k6gen = .false.
            write (6, *) ' Warning to use k6gen recompile with __SCALAPACK precompiler &
                 &  flag !!! '
        end if
#endif

        if (parr .lt. 0.d0 .and. eps_umrigar .eq. -1.d0) then
            eps_umrigar = -parr/10.d0
            write (6, *) ' Default value of eps_umrigar ', eps_umrigar
        end if
        if (parr .gt. 0.d0 .and. eps_umrigar .ne. 0.d0) then
            eps_umrigar = 0.d0
            write (6, *) ' Warning   eps_umrigar=0 in this case , reset to  0'
        end if
        if (beta_learning .eq. -1.d0) then
            !     if(change_tpar) then
            !      beta_learning=0.99d0
            !      write(6,*) ' Default value of beta_learning=',beta_learning
            !     else
            beta_learning = 0.d0
            !     endif
        end if

        if (signalnoise .and. tion .eq. -1.d0) then
            write (6, *) ' Default value of tion =', tpar
            tion = tpar
        end if
        default_epsdgel = .false.
        if (epsdgel .eq. -1.d0) then
            default_epsdgel = .true.
            epsdgel = 0.001d0
            write (6, *) ' Default value of epsdgel =', epsdgel
        end if

        if (npbra .ne. 0) noreadnpbra = .false.
        iflagerrall = iflagerr + iflagerrall
        if (nmore_force .eq. -1) then
            write (6, *) ' Default value of nweight for forces nweight X ', nmore_force
            nmore_force = 1
        end if
        if (nmore_force .gt. 1 .and. iboot .ne. 0) then
            write (6, *) ' Warning iboot=/0 is not possible with nmore_force>0, forced to 0'
            iboot = 0
        end if
        if (ncg .gt. 1 .and. itestrr .ne. -4) then
            write (6, *) ' Warning with SR optimization ncg>1 not possible, ncg changed to 1 '
            ncg = 1
        end if

        if (iboot .lt. 0) then
            write (6, *) "Warning: negative iboot is not allowed! is reset to 0 !"
            iboot = 0
        end if

        if (nweight .le. iboot) then
            write (6, *) "ERROR: nweight too small !!! > iboot >=0 "
            iflagerrall = iflagerr + 1
        end if

        if (tolcg .eq. -1.d0) then
            tolcg = 1d-6 ! does not make any sense a too small accuracy, otherwise the energy
            !                    may increase.
            write (6, *) ' Default value for tolcg   =', tolcg
        end if

        if (minzj .eq. -1000.d0 .and. minz .ne. -1000.d0) then
            minzj = minz
            write (6, *) ' Default value for minzj =', minzj
        end if
        if (maxzj .eq. -1000.d0 .and. minz .ne. -1000.d0) then
            maxzj = maxz
            write (6, *) ' Default value for minzj =', maxzj
        end if
        if (twobodyoff .and. iessz) then
            if (.not. onebodysz) then
                write (6, *) ' Setting one body only also for spin Jastrow '
                onebodysz = .true.
            end if
        end if
        if (minjonetwobody .eq. -1.d0) then
            minjonetwobody = 0.05d0
            write (6, *) ' Default value of minimum one-two body Jastrow =', minjonetwobody
        end if

        !added by K. Nakano for automatic gutta cavat lapidem (to adjust tpar)
        if (change_tpar) then
            write (6, *) ' Warning: tpar is automatically changed on the fly.'
            write (6, *) ' multiply_tpar = ', multiply_tpar
        end if
        !added by K.Nakano

        !          defining default values for readio
        writescratch = 1 ! use Ram memory  to perform optimization
        ! no detailed acceleration
#ifdef   PARALLEL
        wherescratch = './' ! Default path for output scratch
#else
        wherescratch = 'old'
#endif
        unreliable = 0 ! assume that we can deallocate all at the end

        ncore = 0
        ! minimum output for vmc with calculation of variance
        np3 = 0
        ! number of correlation function measured
        np = -1
        ! minimum output
        iread = -1
        ! mandatory section you have to know what
        ifreqdump = 1 ! frequency of dumping

        ! iread=0 standard (no dumping on the file)
        ! iread=2 dumping configurations and weights
        ! iread=3 dumping configurations, weights, local energy and logpsi
        ! (if randomnumber is true dump also the pseudomesh angles)

        nowrite12 = .false.
        iflagerr = 1
#ifdef __KCOMP
        flush_write = .true.
#else
        flush_write = .false.
#endif
        read (5, nml=readio, err=116)
        iflagerr = 0
        write (6, *) ' After reading readio '
116     if (iflagerr .ne. 0) write (6, *) ' ERROR reading readio '
        iflagerrall = iflagerr + iflagerrall

        if (iread .eq. -1 .and. optbra .gt. 3) then
            iread = 6
            write (6, *) ' Default value of  iread ', iread
        elseif (iread .eq. -1) then
            iread = 0
            write (6, *) ' Default value of  iread ', iread
        end if

        if (io_level .eq. 2 .and. iread .eq. 2) then
            write (6, *) ' Warning iread promoted to 3 (you can do everything you do with iread=2) in mpiio case !!!'
            iread = 3
        end if

        if (iread .eq. 3) write (6, *) 'CORRELATED SAMPLING'
        if (iread .eq. 2 .or. iread .eq. 3) then
            write (6, *) 'dumping configurations each', ifreqdump, 'generations'
        end if

        if (trexiofile.ne.'') then
#ifdef _QMCKL
            use_qmckl = (QMCKL_SUCCESS.eq.qmckl_trexio_read(qmckl_ctx, trim(trexiofile), 1_8*len(trim(trexiofile))))
            if (use_qmckl) then
                write (6, *) "Loading TREXIO file:", trexiofile
            else
                write (6, *) "Failed to load TREXIO file:", trexiofile
            end if

            rc = qmckl_get_ao_basis_shell_num(qmckl_ctx, qmckl_shell_num)

            if (rc.ne.QMCKL_SUCCESS) then
                write (0, *) "Failed to get number of shells from TREXIO file"
                use_qmckl = .false.
            end if

            allocate(shell_ang_moms(qmckl_shell_num))

            rc = qmckl_get_ao_basis_shell_ang_mom(qmckl_ctx&
                                               &, shell_ang_moms&
                                               &, 1_8*qmckl_shell_num)

            if (rc.ne.QMCKL_SUCCESS) then
                write (0, *) "Failed to get angular momentum of shells from TREXIO file"
                use_qmckl = .false.
            end if

            allocate(cart_shell_inds(qmckl_shell_num))
            allocate(spher_shell_inds(qmckl_shell_num))

            ! Calculate indeces of cartesian and spherical shells
            cart_shell_inds(1) = 1
            spher_shell_inds(1) = 1
            do ii = 2, qmckl_shell_num
                cart_shell_inds(ii) = cart_shell_inds(ii-1) + cart_multiplicity(shell_ang_moms(ii-1))
                spher_shell_inds(ii) = spher_shell_inds(ii-1) + spher_multiplicity(shell_ang_moms(ii-1))
            end do

#else
            write (6, *) "Ignoring TREXIO file, TurboRVB was not compiled with QMCkl support"
            use_qmckl = .false.
#endif
        end if

        !          default values parameters vmc with calculation of variance

        ieser = -1
        isfix = -1
        iesinv = 0
        iesm = 0
        iesd = 0
        iesfree = 0
        iessw = 0
        iesup = 0
        ieskin = 0
        yes_correct = .true.

        ! mandatory section you have to know
        iflagerr = 1
        yespress = .false. ! no computation of pressure
        warp = .true.
        yespulay = .true.
        add_pulay = -1
        powerwarp = 4.d0
        typedyncell = -2 !  no cell derivatives
        if (lrdmc_der .and. .not. lrdmc_nonodes) then
            scalepulay = 0.d0
        else
            scalepulay = 1.d0
        end if
        scaleeloc = -1.d0 ! standard calculation of energy derivatives without local                     energy der.
        ! by E. Coccia (1/7/11); default values for the external field
        ext_pot = .false.
        vdw = .false.
        link_atom = .false.
        ! by E. Coccia (7/12/11): default value for the MM restraints
        mm_restr = .false.
        ! by E. Coccia (20/12/11): number of initial QMC steps neglected
        ! for averaing dipole, quadrupole and density
        ! default value -> all QMC steps for the average
        ! default value for writing random walk
        write_rwalk = .false.
        ! default value for the reference point for the dipole
        yesavsr = .true.
        yesavcov = .false.

        if (itestr .ne. -5) then
            yesavopt = .false.
        else
            yesavopt = .true.
        end if
        nrep_bead = 1
        yes_kpoints = .false.
        ! summing energies/weigths over k-points.
        epsbas = -1.d0

        cutoff_p = 3.d0 ! user-defined value for the cutoff of bump orbitals in a.u. (numbers 8xx)
        cutoff_p = cutoff_p**2
        if (.not. yesdft) then
            ! with k-points, perform nk decoupled VMC/DMC independent runs without averaging
            if (itestr .ne. -5) then
                decoupled_run = .true.
            else
                decoupled_run = .false.
            end if
        else
            decoupled_run = .false.
        end if
        fixa = .false.
        fixb = .false.
        fixc = .false.
        real_contracted = .true.
        real_agp = .false.
        no_sjbra = .false.
        pressfixed = 0.d0
        read_molecul = .false.
        epsder = 1d-5
        yes_scemama = .true.
        yes_scemama_open = .false. ! Better do not  do approximation by  default.

        read (5, nml=parameters, err=117)
        if (.not. ext_pot) then
            vdw = .false.
            link_atom = .false.
        end if

        if (epsbas .eq. -1.d0) then
            if (yesdft) then
                epsbas = 1.0d-8 ! default value for cutoff in lattice vector summation
                ! periodic basis set with k-points
            else
                epsbas = 1.0d-7
            end if
            if (yes_scemama) epsbas = epsbas/10.d0
            write (6, *) ' Default value of epsbas =', epsbas
        end if
        if (link_atom) vdw = .true.
        if (link_atom) mm_restr = .false.
        if (mm_restr) link_atom = .false.

        write (6, *) ' After reading parameters '
        iflagerr = 0
117     if (iflagerr .ne. 0) write (6, *) ' ERROR reading parameters '
        iflagerrall = iflagerr + iflagerrall
        if (yes_kpoints) then
            manyfort10 = .true.
        else
            manyfort10 = .false.
        end if

        if (nbead .gt. 1 .and. yesdft .or. (.not. yesdft .and. (yesread10 .or. yeswrite10))) then
            manyfort10 = .true.
        end if

        if (nrep_bead .gt. 1 .and. yesquantum) then
            if (yesavsr) then
                write (6, *) ' Warning yesavsr forced to .false. with nrep_bead>1'
                yesavsr = .false.
            end if
            if (yesavopt) then
                write (6, *) ' Warning yesavopt forced to .false. with nrep_bead>1'
                yesavopt = .false.
            end if
        elseif (.not. yesavsr .and. yesavopt) then
            write (6, *) ' Warning yesavsr forced to .true.  with yesavopt=.true.'
            yesavsr = .true.
        end if

        !---------Berry phase----------------------------------------

        ! Berry phase should not be computed
        ! during a minimization run
        if (add_pulay .eq. -1) then
            add_pulay = 2
            write (6, *) ' Default value of add_pulay =', add_pulay
        end if
        !------------------------------------------------------------

        if (itestr .ne. -5) then

            if (ieser .eq. -1) then
                if (itestr .eq. -6 .or. itestr4 .eq. -7) then
                    ieser = 6
                else
                    ieser = 1
                end if
                write (6, *) ' Default value for ieser=', ieser
            end if
            if (isfix .eq. -1) then
                if (parcutg .eq. 0) then
                    isfix = 1
                    write (6, *) ' Default value for isfix=', isfix
                else
                    if (gamma .ne. 0.d0 .and. optbra .ne. 3) then
                        isfix = 3
                    else
                        isfix = 1
                    end if

                    write (6, *) ' Default value for isfix=', isfix
                end if
            end if
        else
            if (ieser .eq. -1) then
                ieser = 0
                write (6, *) ' Default value for ieser=', ieser
            end if
            isfix = 0 ! no variance during minimization
        end if

        write (6, *) ' Basis set cutoff chosen:', epsbas
        write (6, *) ' after parameters   '
        !          default values for unused
        ! no special type of DMC/LRDMC
        rsignr = 0.d0
        ! no special way to deal with the Hessian
        beta = 0.d0
        ! no testing derivative
        testderiv = 0
        if (developer .gt. 0) then
            iflagerr = 1
            read (5, nml=unused, err=118)
            iflagerr = 0
            write (6, *) ' After reading unused    '
118         if (iflagerr .ne. 0) write (6, *) ' ERROR reading unused '
            iflagerrall = iflagerr + iflagerrall
        end if
        !          if(molopt.eq.0.and.testderiv.ge.0) testderiv=-1-testderiv

        ! by E. Coccia (22/11/10): read the flags
        ! default value
        !ext_pot = .false.
        !vdw = .false.
        !link_atom=.false.
        !iflagerr=1
        !read(5,nml=pot_ext,err=900)
        !iflagerr=0
        !if (.not.ext_pot) then
        !   vdw=.false.
        !   link_atom=.false.
        !endif
        !if (link_atom) vdw=.true.
        !write(6,*) 'After reading the external potential'
        if (ext_pot) then
            write (*, *) ''
            write (6, *) '|**************************************|'
            write (6, *) '| ADDING AN EXTERNAL QMC/MM POTENTIAL  |'
            write (6, *) '|**************************************|'
            write (*, *) ''
        end if !900        if (iflagerr.ne.0) write(6,*) 'ERROR reading the external potential'

        if (link_atom) then
            iflagerr = 1
            read (5, nml=link, err=898)
            iflagerr = 0
898         if (iflagerr .ne. 0) write (6, *) 'ERROR reading the link namelist'
        end if

        write (6, *) ' Parameters: iesinv,iesm,iesd,iesfree,iessw,iesup,ieskin '
        write (6, *) ' Parameters before read ', iesinv, iesm, iesd, iesfree, iessw, iesup, ieskin

        !
        ! correct input parameters if necessary
        ! by reading the first part of the wave function
        !
        call read_fort10_fast

        if (parcutg .lt. 2) then
            if (itestr .ne. -6 .and. itestr4 .ne. -7 .and. itestr4 .ge. -20) then
                if (cutreg .eq. 0) then
                    if (yesalfe) then
                        cutreg = 0.4d0*dsqrt(dble(nel)/tbra)
                    else
                        cutreg = 1.d0/dsqrt(tbra)
                    end if
                elseif (cutreg .gt. 0.d0) then
                    cutreg = 2.d0*cutreg ! One can put by hand whatever he wants, changed unit from H to Ry.
                end if
                if (.not. better_dmc) cutreg = 0.d0
            elseif (cutreg .ne. 0.d0) then
                cutreg = 0.d0
                write (6, *) ' Warning  cutreg set to zero', cutreg
            end if
        end if

        if (better_dmc) then
            if (cutreg .gt. 0) then
                if (.not. yesalfe) then
                    write (6, *) ' Warning DMC cutoff on local energy (Ry)=', cutreg*zmax**2
                else
                    write (6, *) ' Warning DMC cutoff on local energy (Ry)=', cutreg
                end if
            else
                if (cutreg .ne. -2.d0) then
                    write (6, *) ' Warning NO DMC cutoff on the local energy'
                else
                    write (6, *) ' Warning Umrigar 1993 cutoff on the local energy'
                end if
            end if
        else
            if (cutreg .ne. -2.d0) then
                write (6, *) ' Warning NO DMC cutoff on the local energy'
            else
                write (6, *) ' Warning Umrigar 1993 cutoff on the local energy'
            end if
        end if

        if (parcutg .eq. 0) then
            novar = 0
        else
            !  With parcutg =/0  H^a within LRDMC is regularized with the Guiding function
            !  so that we can sample the region where psi_T ~0.
            if (itestr .ne. -6 .and. itestr4 .ne. -7 .and. itestr4 .ge. -20) then
                write (6, *) ' Warning you cannot use the guiding =/ trial in DMC !!!'
                parcutg = 0
                novar = 0
                write (6, *) ' parcutg set to zero', parcutg
            else
                if (alat .eq. 0.d0) then
                    alat = -1.d0/zmax
                    write (6, *) ' Default value for alat =', alat
                end if
                if (novar .eq. -1) then
                    novar = 0
                    write (6, *) ' Default value of novar =', novar
                end if
            end if
        end if

        if (parcutg .eq. 2 .and. cutreg .eq. 0.d0) then
            if (zmax .gt. 0.d0) then
                cutreg = abs(alat)/sqrt(zmax)
                write (6, *) ' Default value of cutreg =', 1.d0/sqrt(zmax)
            end if
        elseif (parcutg .eq. 2) then
            write (6, *) ' Warning input cutreg =', cutreg
            cutreg = cutreg*abs(alat)
        end if

        ! defining default values for vmc
        if (epscuttype .eq. -1000) then
            if (epscut .eq. 0.d0) then
                epscuttype = 0
            elseif (epscut .gt. 0.d0) then
                epscuttype = 2
            end if
            write (6, *) ' Default value of epscuttype= ', epscuttype
        else
            write (6, *) ' Read value of epscuttype= ', epscuttype
        end if
        if (epscuttype .gt. 2 .or. epscuttype .lt. 0) then
            write (6, *) ' ERROR this epscuttype does not exist !!! '
            iflagerr = 1
            iflagerrall = iflagerr + iflagerrall
            !            check if the option exists otherwise set iflagerr=1
        end if

        if (epscut .eq. 0 .and. epscuttype .ne. 0) then
            if (epscuttype .le. 0) epscut = epscuttype
            if (epscuttype .gt. 0) then
                write (6, *) ' ERROR you should define epscut in this case !!! '
                iflagerr = 1
                iflagerrall = iflagerrall + 1
            else
                write (6, *) ' Default value of epscut=', epscut
            end if
        else
            write (6, *) ' Read value of epscut=', epscut
        end if

        if (tstep .eq. -1.d0) then
            if (itest .eq. 2) then
                tstep = 2.d0 ! reasonable universal
            else
                tstep = 0.d0
            end if
            write (6, *) ' Default value of tstep=', tstep
        else
            write (6, *) ' Read value of tstep=', tstep
        end if

        epstl = epstlrat*epscut

        if (iflagerr .eq. 0) then ! if the read was OK
            if (ireadminr .gt. 0 .and. yesfast .gt. 0) then
                write (6, *) ' Warning in reading parameters yesfast changed to Default '
                yesfast = -1
            end if

            if (ieskinold .ne. 0) then
                if (warp) then
                    write (6, *) ' Default used warp for forces '
                else
                    write (6, *) ' Warning no warp used for forces '
                end if
            end if

            if (yespress .and. iespbc .and. typedyncell .eq. -2 .and. itestr .eq. -5) then
                typedyncell = 3
                write (6, *) ' Warning  changed typedyncell with yespress=.true.=', 3
            end if

            if (.not. iespbc .and. (yespress .or. typedyncell .gt. 0)) then
                if (yespress) &
                     &write (6, *) ' Warning no calculation of pressure with open systems '
                yespress = .false.
                if (typedyncell .gt. 0) &
                     &write (6, *) ' Warning no cell derivatives with open systems'
                typedyncell = 0
            end if

            if ((typedyncell .ge. 0 .or. yespress) .and. ieskinold .eq. 0) then
                if (yespress) &
                    write (6, *) ' Warning no calculation of pressure without ieskin=1 '
                yespress = .false.
                if (typedyncell .ge. 0) &
                     &write (6, *) ' Warning no calculation of cell derivatives/forces  &
                     &without ieskin=1&
                     & in the parameters section '
                typedyncell = -1
            end if

            ieskinrp = ieskinr_pos
            if (yespress .and. iespbc .and. itestr .ne. -5) then
                ieskinrp = ieskinr_pos + 1
                ieskin = ieskin + 1
            end if

            if (ireadminr .ne. 0 .and. ngen .ge. nweight) then
                write (6, *) ' Warning ngen replaced to the maximum allowed ', nweight - 1
                ngen = nweight - 1
            elseif (mod(ngen, nweight) .ne. 0) then
                call error('main', 'you cannot continue the run!!!', -1, rank)
            end if

            if (npsa .eq. -1) npsa = npsar
            if (npsa .ne. 0 .and. nintpsa .eq. 0) then
                npsa = npsar
                if (yesdft) then
                    nintpsa = 50
                else
                    if (zmax .gt. 18) then ! from N on. up to third raw  included
                        nintpsa = 18
                    elseif (zmax .gt. 6) then
                        nintpsa = 12
                    else
                        nintpsa = 6
                    end if
                end if
                write (6, *) ' Default value for nintpsa ', nintpsa
            else
                write (6, *) ' Read value for nintpsa ', nintpsa
            end if
            if (npsa .ne. 0 .and. npsamax .eq. 0) then
                npsamax = 2
                if (npsa .eq. 1) npsamax = 1
                write (6, *) ' Default value for npsamax ', npsamax
            end if

            if (itest .eq. 1 .or. itestr4 .eq. -7) then ! begin if DMC/LRDMC

                if (itestr .eq. -6 .or. itestr4 .eq. -7 .or. itestr4 .lt. -20) then ! test if LRDMC
                    if (plat(1) .eq. -1000.d0) then
                        ! default value
                        if (alat .gt. 0.d0 .or. alat2 .gt. 0) then

!!!! new !!!! added by K.Nakano 11/09/2019
                            !                 if(npsa.gt.0) then
                            !                  plat(1)=-4.d0/25.d0
                            !                 else
                            if (npsa .ne. 0) then
                                r_c = 1.d0
                                plat(1) = -2.d0*r_c**2 ! default length of valence region r> r_c
                                if (rank .eq. 0) write (6, *)&
                                   & ' Default length of valence region r> rc,'//&
                                   & ' change it by  plat(1)=-2xrc^2 in input  =', r_c
                            else

                                beta_for_r_c = 0.75d0
                                kappa_for_r_c = 2.50d0
                                r_c = beta_for_r_c*zmax**(-5.0d0/7.0d0)&
                                    &*(kappa_for_r_c*(zmax*alat)**2.0d0 + 1.0d0)&
                                    &/((zmax*alat)**2.0d0 + 1.0d0)
                                plat(1) = -2.0d0*r_c**2.0d0
                                write (6, *) ' r_c and alat2 are given by the Thomas-Fermi theory'
                                !                 endif
!!!! new !!!! added by K.Nakano 11/09/2019

                            end if

                        else
                            plat(1) = 0.d0
                        end if
                        write (6, *) ' Default value for plat(1)=', plat(1)
                    elseif (alat2 .eq. -1.d0) then
                        if (plat(1) .gt. 0) then
                            r_c = plat(1)
                            plat(1) = -2.0d0*r_c**2.0d0
                        end if
                        write (6, *) ' Default value for plat(1)=', plat(1)
                    end if
                    if (alat2 .eq. -1.d0) then ! default value two lattice spaces
                        if (alat .gt. 0) then
                            !                  if(npsa.gt.0)  then
                            !                      alat2=sqrt(7.d0)
                            !                  else
!!!! new !!!! added by K.Nakano 11/09/2019
                            num_ele_core_r_c = Thomas_Fermi_core_electron_number(zmax, r_c) - core_pseudo
                            if (npsa .eq. 0) then
                                write (6, *) ' num_ele_core=', num_ele_core_r_c
                                write (6, *) ' num_ele_core/atomic_number=', num_ele_core_r_c/zmax

                                alat2 = dsqrt(l0_kousuke*(zmax - core_pseudo - num_ele_core_r_c)/num_ele_core_r_c)
!!!! new !!!! added by K.Nakano 11/09/2019
                            else
                                if (zmax .le. 18) then
                                    alat2 = sqrt(5.d0)
                                    write (6, *)&
                                       & ' Warning default alatp/alat for pseudo sqrt(5.)  !!!'//&
                                       & ', change it by alat2=better ratio in input'
                                else
                                    alat2 = sqrt(7.d0)
                                    write (6, *)&
                                       & ' Warning default alatp/alat for pseudo sqrt(7.)  !!!'//&
                                       & ', change it by alat2=better ratio in input'
                                end if
                            end if

                            !alat2=sqrt(7.d0)*(zmax/8.d0)**(1.D0/6.D0)

                            !                  endif
                            if (abs(nint(alat2) - alat2) .lt. epsmach) then
                                alat2 = sqrt(alat2**2 + 0.1d0)
                            end if
                        else
                            alat2 = 0.d0
                        end if
                        if (alat2 .eq. 0) then
                            write (6, *) ' Default single mesh '
                        else
                            write (6, *) ' Default double mesh alat2= ', alat2
                        end if
                    end if
                    if (alat2 .ne. 0.d0 .and. alat .gt. 0) then

                        !      check if alat2 is irrational
                        cost = 0.d0
                        i = 1
                        do while (cost .eq. 0.d0 .and. i .le. 1000000)
                            if (abs(nint(i*alat2) - i*alat2) .lt. 1d-6) cost = i*alat2
                            i = i + 1
                        end do

                        if (cost .ne. 0.d0) then
                            write (6, *) ' Warning alat2 rational number ', nint(cost), '/', nint(cost/alat2)
                            write (6, *) ' Minimum lattice space =', alat/(cost/alat2)
                        end if
                    end if
                    if (klrdmc .eq. -1000.d0) then
                        if (npsa .ne. 0) then
                            klrdmc = 1.d0
                        else
                            klrdmc = 0.d0
                        end if
                        write (6, *) ' Default Value for Klrdmc=', klrdmc
                    end if
                    if (plat(2) .eq. -1000.d0) then
                        plat(2) = 1.d0 + Klrdmc*alat**2
                        write (6, *) ' Default Value for eta=plat(2)=', plat(2)
                    end if

                    if (plat(3) .eq. -1000.d0 .and. alat2 .ne. 0.d0) then
                        plat(3) = plat(2)
                        write (6, *) ' Default value for plat(3)=', plat(3)
                    end if
                    if (plat(3) .eq. -1000.d0 .and. alat2 .eq. 0.d0) then
                        plat(3) = 0.d0
                        write (6, *) ' Default value for plat(3)=', plat(3)
                    end if
                    if (tbra .eq. 0.d0) then
                        tbra = 10.d0/zmax**2
                        write (6, *) ' Default value for tbra =', tbra
                    end if
                else
                    if (tbra .eq. 0.d0) then
                        tbra = 1.d0/zmax**2
                        write (6, *) ' Default value for tbra =', tbra
                    end if
                end if ! endif for lrdmc

            end if ! endif for itest
            if (ieskinold .ne. 0) then
                if ((idyn .ne. 0 .or. signalnoise) .and. typedyncell .eq. -2) then
                    typedyncell = 0
                    write (6, *) ' Default value of typedyncell (NVT) ', typedyncell
                end if

                if (typedyncell .eq. 0) then
                    ! default value
                    if (idyn .ne. 0) npbra = npbra + ieskin
                    ! default value for ncore
                    if (ncore .eq. 0 .and. (idyn .ne. 0 .or. (signalnoise .and. ieskin .ne. 0))) ncore = -1
                end if

                if (typedyncell .eq. -1) then
                    npbra = 0
                    ncore = 0
                end if

                if (typedyncell .eq. 1) then
                    !   Volume constant
                    ieskin = ieskin + 2
                    ! default value
                    if (idyn .ne. 0) npbra = npbra + ieskin
                    ! default value
                    if (ncore .eq. 0 .and. (idyn .ne. 0 .or. (signalnoise .and. ieskin .ne. 0))) ncore = -2
                end if

                if (typedyncell .eq. 2) then
                    !   Pressure constant
                    ieskin = ieskin + 3
                    ! default value
                    if (ncore .eq. 0 .and. (idyn .ne. 0 .or. (signalnoise .and. ieskin .ne. 0))) ncore = -2
                    ! default value
                    if (idyn .ne. 0) npbra = npbra + ieskin
                end if

                if (typedyncell .eq. 3) then
                    !   Isotropic Pressure constant
                    ieskin = ieskin + 1
                    ! default value
                    if (ncore .eq. 0 .and. (idyn .ne. 0 .or. (signalnoise .and. ieskin .ne. 0))) ncore = -2
                    ! default value
                    if (idyn .ne. 0) npbra = npbra + ieskin
                end if

            end if ! endif ieskin ne 0

            if (typedyncell .gt. 1 .and. yesquantum) then
                write (6, *) ' Warning quantum presssure correction not included in dynamics !!! &
                     & but only in output writing pressure.dat'
            end if

            if (ipc .eq. 1 .or. itest .ne. 2) then
                iese = abs(ieser)
            elseif (ieser .ne. 0) then
                iese = abs(ieser) + 1 ! only in vmc we compute the imaginary part also
                ! for all observables required
            else
                iese = 0
            end if

            write (6, *) ' Parameters after  read '                            &
                 &, iese, iesinv, iesm, iesd, iesfree, iessw, iesup, ieskin

            ieskint = ieskin + iesking

            ! default value for np np3
            if (np .eq. -1) then
                if (itestr .ne. -5) then
                    np = iese + abs(isfix) +&
                         & 2*(abs(iesinv) + iesm + iesd + abs(iesfree) + abs(iessw) + iesup + iesking)&
                         &+ 3*ieskin
                    np3 = np
                else
                    np = iese + abs(iesinv) + abs(isfix) + iesm + iesd + abs(iesfree)      &
                         &+ abs(iessw) + iesup + iesking + ieskin

                    ! no output in the optimization
                    if (ncg .eq. 0 .and. noreadnpbra .and. itestrr .eq. -4) then
                        npbra = npbra + np - ieskin
                        write (6, *) ' Default full Hessian npbra=', npbra
                    end if
                    np3 = 0
                end if
            end if

            if (abs(kl) .eq. 7 .and. itestr .eq. -5 .and. prep .eq. -1) then
                prep = nint(sqrt(dble(np)/dble(nweight)))
                write (6, *) ' Default value of prep', prep
                if (prep .le. 1) prep = -1
            end if

            def_nscra = .false.
            if (nscra .eq. 0) then
                def_nscra = .true.
                nscra = 2*nel + 2 ! +2 , Otherwise sometimes it calls an unnecessary scratchdet before compute_fast
                if (ipf .eq. 2) nscra = nel/2
                write (6, *) ' Default value for nscra ', nscra
            else
                write (6, *) ' Read value for nscra ', nscra
            end if
            if (nbra .eq. 0) then
                if (itest .eq. 2) then
                    if (itestr .eq. -5) then
                        nbra = 4*nel
                    else
                        nbra = 2*nel
                        if (ieskin .ne. 0) nbra = nbra*2
                    end if
                elseif (itest .eq. 1) then
                    if (itestr .ne. -6 .and. itestr .ne. -7 .and. itestr4 .ge. -20) then
                        nbra = nel
                    else
                        nbra = 1 ! default for lrdmc used only for output acceptance moves
                    end if
                end if
                write (6, *) ' Default value for nbra ', nbra
            else
                write (6, *) ' Read value for nbra ', nbra
            end if

            if (ncg .gt. np) then
                ncg = np
                write (6, *) ' Warning ncg<= #parameters!, changed to ', ncg
            end if

            if (npbra .gt. np) then
                npbra = max(np - ncg, 0)
                write (6, *) ' Warning npbra<= #parameters!, changed to ', npbra
            end if

            if (.not. yesdft) then
                open (unit=7, file='parminimized.d', form='formatted', status='unknown')
                j = abs(iesinv) + iesm + iesd + abs(iesfree) + abs(iessw) + iesup + iesking
                if (ipc .eq. 2 .and. .not. yes_correct) then
                    write (7, '(12I13)') iesinv, iesm, iesd, iesfree, iessw, iesup, ieskin, j, isfix, iese, -ipc, iesking
                else
                    write (7, '(12I13)') iesinv, iesm, iesd, iesfree, iessw, iesup, ieskin, j, isfix, iese, ipc, iesking
                end if

                if (yesquantum) then
                    if (itest .eq. 1 .and. itestr .ne. -6 .and. itestr .ne. -7 .and. itestr4 .ge. -20) then
                        !      for standard FN
                        write (7, '(2I13,6e18.9)') - nw, optbra, gamma, tbra*dble(nbra/nel), etry, plat(2), alat, tstepfn
                    else
                        write (7, '(2I13,6e18.9)') - nw, optbra, gamma, tbra, etry, plat(2), alat, tstepfn
                    end if
                else
                    if (itest .eq. 1 .and. itestr .ne. -6 .and. itestr .ne. -7 .and. itestr4 .ge. -20) then
                        !      for standard FN
                        write (7, '(2I13,6e18.9)') nw, optbra, gamma, tbra*dble(nbra/nel), etry, plat(2), alat, tstepfn
                    else
                        write (7, '(2I13,6e18.9)') nw, optbra, gamma, tbra, etry, plat(2), alat, tstepfn
                    end if
                end if
                if (idyn .eq. 0) close (7)
            end if

            if (itestr4 .ge. -10) then
                !      Normal minimization
                itestr3 = itestr4
                !      fncont=.false.
            elseif (itestr4 .ge. -20) then
                itestr3 = itestr4 + 10
                !      Minimization with continuous DMC, working without pseudo, deprecated
                !      fncont=.true.
            elseif (itestr4 .ge. -30) then
                !      Minimization with LRDMC, working also  with pseudo
                itestr3 = itestr4 + 20
                !      fncont=.false.
            end if
            powermin = 0
            npower = 0
            powerminsz = 0
            npowersz = 0

            !      Default values for section fitpar
            ! no fit of the long range part > rmaxinv
            nparinv = 0
            ! piecewise constant parametrization
            initparinv = -1
            ! only the local part
            rmaxinv = 1d-10
            ! no fit of the long range part > rmaxj
            npar = 0
            ! piecewise constant parametrization
            initpar = -1
            ! only the local part
            rmaxj = 1d-10
            ! no fit of the long range part > rmax
            nparsw = 0
            ! piecewise constant parametrization
            initparsw = -1
            ! only the local part
            rmax = 1d-10

            allfit = .true.

            if (iesfree .lt. 0 .or. iessw .lt. 0 .or. iesinv .lt. 0) then
                iflagerr = 1

                read (5, nml=fitpar, err=119)
                iflagerr = 0
                write (6, *) ' After reading fitpar '
119             if (iflagerr .ne. 0) write (6, *) ' ERROR reading fitpar '
                iflagerrall = iflagerr + iflagerrall
            end if
            if (npower .gt. 0 .and. powermin .lt. 0) then
                write (6, *) ' Warning  divergent Jastrow powermin<0 '
            end if
            if (npowersz .gt. 0 .and. powerminsz .lt. 0) then
                write (6, *) ' Warning  divergent spin Jastrow powermin<0 '
            end if

            if (iesfree .lt. 0 .and. iessw .lt. 0 .and. iesinv .lt. 0) then
                write (6, *) ' Parametrization lambda matrices Jas Sz Jas and Det '
                !      read(5,*)   nparinv,initparinv,rmaxinv,npar,initpar,rmaxj
                !    1,nparsw,initparsw,rmax
                iesinv = abs(iesinv)
                iesfree = abs(iesfree)
                iessw = abs(iessw)
                if (rmax .le. 0) rmax = 1d-10
                if (rmaxj .le. 0) rmaxj = 1d-10
                if (rmaxinv .le. 0) rmaxinv = 1d-10
            elseif (iesfree .lt. 0 .and. iessw .lt. 0) then
                write (6, *) ' Parametrization lambda matrices Jas and Det '
                !      read(5,*) npar,initpar,rmaxj,nparsw,initparsw,rmax
                iesfree = abs(iesfree)
                iessw = abs(iessw)
                nparinv = 0
                initparinv = 0
                rmaxinv = 0.d0
                if (rmax .le. 0) rmax = 1d-10
                if (rmaxj .le. 0) rmaxj = 1d-10
            elseif (iesinv .lt. 0 .and. iessw .lt. 0) then
                write (6, *) ' Parametrization lambda matrices Jas Sz and Det '
                !      read(5,*) nparinv,initparinv,rmaxinv,nparsw,initparsw,rmax
                iesinv = abs(iesinv)
                iessw = abs(iessw)
                npar = 0
                initpar = 0
                rmaxj = 0.d0
                if (rmax .le. 0) rmax = 1d-10
                if (rmaxinv .le. 0) rmaxinv = 1d-10
            elseif (iesinv .lt. 0 .and. iesfree .lt. 0) then
                write (6, *) ' Parametrization lambda matrices Jas Sz and Jas  '
                !      read(5,*) nparinv,initparinv,rmaxinv,npar,initpar,rmaxj
                iesinv = abs(iesinv)
                iesfree = abs(iesfree)
                nparsw = 0
                initparsw = 0
                rmax = 0.d0
                if (rmaxj .le. 0) rmaxj = 1d-10
                if (rmaxinv .le. 0) rmaxinv = 1d-10
            elseif (iesfree .lt. 0) then
                write (6, *) ' Parametrization lambda matrix  Jastrow   '
                !      read(5,*) npar,initpar,rmaxj
                iesfree = abs(iesfree)
                nparsw = 0
                initparsw = 0
                rmax = 0.d0
                nparinv = 0
                rmaxinv = 0.d0
                initparinv = 0
                if (rmaxj .le. 0) rmaxj = 1d-10
            elseif (iessw .lt. 0) then
                write (6, *) ' Parametrization lambda matrix  Det  '
                !      read(5,*) nparsw,initparsw,rmax
                iessw = abs(iessw)
                npar = 0
                initpar = 0
                rmaxj = 0.d0
                nparinv = 0
                initparinv = 0
                rmaxinv = 0.d0
                if (rmax .le. 0) rmax = 1d-10
            elseif (iesinv .lt. 0) then
                write (6, *) ' Parametrization lambda matrix  Jas-Sz  '
                !      read(5,*) nparinv,initparinv,rmaxinv
                iesinv = abs(iesinv)
                npar = 0
                initpar = 0
                rmaxj = 0.d0
                nparsw = 0
                initparsw = 0
                rmax = 0.d0
                if (rmaxinv .le. 0) rmaxinv = 1d-10
            else
                npar = 0
                nparsw = 0
                nparinv = 0
                initpar = 0
                initparsw = 0
                initparinv = 0
                rmax = 0.d0
                rmaxj = 0.d0
                rmaxinv = 0.d0
            end if

            if (initpar .lt. -1) then
                write (6, *) ' VdW parametrization of the charge Jastrow factor '
                if (mod(npar, 4) .ne. 0) then
                    write (6, *) ' ERROR you should have npar multiple of 4 !!! '
                    iflagerr = 1
                    iflagerrall = iflagerr + iflagerrall
                end if
            end if
            if ((npar .ne. 0 .or. nparsw .ne. 0 .or. nparinv .ne. 0 .or. npower .ne. 0&
                 &.or. npowersz .ne. 0) .and. abs(kl) .ne. 7) then
                write (6, *) ' Warning |kl| changed to 7 with parametrization '
                if (kl .gt. 0) then
                    kl = 7
                else
                    kl = -7
                end if
            end if
            write (6, *) ' kl read =', kl

            iflagk = 0
            if (iespbc) then
                !        read(5,*,end=1056) kSq,kappar
                !        iflagk=1
                !1056   continue
                if (kSq .eq. 0.d0) then
                    kSq = 1.d-8
                    write (6, *) ' Default values for kSq (precision Ewald) ', kSq
                end if
                if (kappar .eq. 0.d0) then
                    !             if(yes_tilted)  then
                    !             kappar=16
                    !             else
                    if (ksq .eq. 0.5d0) then
                        if (yes_tilted) then
                            kappar = 16.d0
                        else
                            kappar = 8.d0
                        end if
                    else
                        kappar = 1.d0
                        err_derfc = derfc(kappar)
                        do while (err_derfc .gt. ksq)
                            kappar = kappar + 0.5d0
                            err_derfc = derfc(kappar)/kappar
                        end do
                        kappar = kappar*2
                    end if
                    !             endif
                    write (6, *) ' Default value for kappar=', kappar
                end if
                if (neigh .eq. -1) then
                    neigh = 1
                    write (6, *) ' Default value of neighbors in Ewald ', neigh
                end if
            else
                kSq = 0.d0
                kappar = 0.d0
                rs = -1.d0
            end if
            write (6, *) ' before dynamic '

            !       default values
            temp = 0.d0
            if (yesquantum) then
                yesturboq = .true.
                yessecond = .true.
            else
                yesturboq = .false.
                yessecond = .false.
            end if
            yesrootc = .false.
            cleanrognoso = .false.
            if (idyn .eq. 5) then
                addrognoso = .true.
            else
                addrognoso = .false.
            end if
            normcorr = 1.d0
            maxdev_dyn = 0.d0
            friction = 0.d0
            stepcg_recount = .false.
            write_cov = .false.
            iskipdyn = 1
            scalecov = 1.d0

            delta0 = 0.d0
            delta0q = 0.d0
            delta0k = 0.d0
            smoothcut = 0.d0
            killcut = .false.
            scale_mass = 1.d0

            eqcellab = .false.
            eqcellac = .false.
            eqcellbc = .false.

            if (idyn .ne. 0) then
                iflagerr = 1
                read (5, nml=dynamic, err=120)
                write (6, *) ' After reading dynamic '
                iflagerr = 0
120             if (iflagerr .ne. 0) write (6, *) ' ERROR reading dynamic '
                iflagerrall = iflagerr + iflagerrall
                if (scalecov .le. 0.d0) scalecov = 1.d0
                ! if read value of temp is negative, its absolute value is interpreted
                if ((idyn .eq. 7 .or. idyn .eq. 8) .and. delta0q .eq. 0) then
                    delta0q = friction
                    write (6, *) ' Default value of Ceriotti 1/tau0 =', delta0q
                end if
                if ((idyn .eq. 7 .or. idyn .eq. 8) .and. delta0k .eq. 0) then
                    delta0k = 1.d0
                    write (6, *) ' Default value of delta0k (Ceriotti choice) =', delta0k
                end if
                ! as the temperature in K
                if (temp .le. 0.d0) temp = -0.5*temp*kboltz/rydberg

                write (*, *) ' Temperature (Kelvin) = ', temp*2*rydberg/kboltz

                !       if(idyn.eq.5.and.normcorr.ne.0.d0.and.tion.gt.2.d0/normcorr) then
                !      write(6,*) ' Warning you should use input tion <=2 Temp/normcorr&
                !       to apply the noise correction! , i.e. tion<=',2.d0*temp/normcorr
                !           normcorr=0.d0
                !          endif

                if (idyn .eq. 8 .and. nbead .gt. 1) then
                    yesquantum = .true.
                    yesturboq = .true.
                    yessecond = .true.
                    oldscaling = .true. ! idyn 8 does not care about mass reference!
                elseif (idyn .eq. 8) then
                    yessecond = .true.
                    oldscaling = .true. ! idyn 8 does not care about mass reference !
                end if
                if (idyn .eq. 5 .and. yesquantum .and. .not. yesturboq) then
                    write (6, *) ' Warning yesturboq = true with idyn=5 and quantum, changed '
                    yesturboq = .true.
                end if

                if (yesturboq .and. .not. yesquantum) then
                    yesturboq = .false.
                    write (6, *) ' Warning yesturboq = .false. in the classical case !!! '
                end if
                if (nbead .eq. -1 .and. yesquantum) then
                    nbead = min(nproc, 16)
                    write (6, *) ' Default value of beads =', nbead
                end if
                if (idyn .ne. 0) then
                    ! writing information about the quantum simulation
                    write (7, *) temp, nbead, nion, yeswritebead, yeswrite10
                    close (7)
                end if

                !       if(yesquantum.and.yesturboq.and..not.yesavcov) then
                !       write(6,*) ' Warning yesavcov forced to true in this case'
                !       yesavcov=.true.
                !       endif

                if (yesquantum .and. yesturboq .and. yesperiodize) then
                    write (6, *) ' Warning yesturboq force to false as periodization &
                         & does not work with yesturboq true'
                    yesperiodize = .false.
                end if

                if (yesquantum&
                    & .and. yesturboq&
                    & .and. idyn&
                    & .ne. 0&
                    & .and. (idyn .ne. 6 .and. idyn .ne. 7 .and. idyn .ne. 8 .and. idyn .ne. 5)) then
                    write (6, *) ' Warning forced to idyn=7 in this case'
                    idyn = 7
                end if

                ! correct obvious misprints
                if (iskipdyn .le. 0) iskipdyn = 1
                if (nmore_force .gt. 1) iskipdyn = iskipdyn + nmore_force - 1

                if (npbra .lt. ieskin) then
                    write (6, *) ' npbra > = ieskin !!! ', npbra, ieskin
                    iflagerr = 1
                    iflagerrall = iflagerr + iflagerrall
                end if

                if (stepcg_recount .and. ncg .gt. iskipdyn) ncg = iskipdyn

            end if

            if (idyn .eq. 2) then
                write (6, *) ' 2nd order Langevin dynamic with friction =', friction
                write (6, *) ' Move ions each ', iskipdyn, 'steps'
            elseif (idyn .eq. 3 .or. idyn .eq. 4 .or. idyn .eq. 6) then
                write (6, *) ' 2nd order Langevin dynamic with noise correction,    &
                     &  friction =', friction
                if (idyn .eq. 4) write (6, *) ' Delta t < 3/4 Delta_0 (unit Ry^-1) =', delta0
                if (idyn .eq. 3) write (6, *) ' Delta t < Delta_0= (unit Ry^-1) ', delta0
                write (6, *) ' Move ions each ', iskipdyn, 'steps'
                write (6, *) ' Matrix covariance scaled by ', scalecov
            elseif (idyn .eq. 8) then
                if (yesquantum) then
                    write (6, *) ' IDYN=8: Langevin based PIMD with Trotter breakup'
                    write (6, *) ' exact harmonic propagation of quantum part with optimal damping'
                    write (6, *) ' centroid friction =', delta0q
                    write (6, *) ' Born-Oppenheimer friction =', friction
                else
                    write (6, *) ' IDYN=8: classical 2nd order Langevin dynamic with noise correction'
                    write (6, *) ' friction =', friction
                end if
                write (6, *) ' Move ions each ', iskipdyn, 'steps'
                write (6, *) ' Matrix covariance scaled by ', scalecov
            elseif (idyn .eq. 5) then
                write (6, *) ' Structural optimization or dynamics  with covariance matrix '
                write (6, *) ' Move ions each ', iskipdyn, 'steps'
            end if

            if (eqcellab) then
                eqcellbc = .false.
                eqcellac = .false.
            elseif (eqcellbc) then
                eqcellac = .false.
                eqcellab = .false.
            elseif (eqcellac) then
                eqcellab = .false.
                eqcellbc = .false.
            end if

            if (ieskint .ne. 0) then

                write (*, *) ' Differential-Warp  nuclear Forces!'

                if (add_pulay .eq. 0) then
                    write (*, *) ' No pulay for ionic forces '
                elseif (add_pulay .eq. 1) then
                    write (*, *) ' Only three-body and Jastro pulay               &
                         & for ionic forces '
                else
                    write (*, *) ' Full Pulay for ionic forces '
                end if

            end if

            if (itestr .eq. -5) then ! only for minimization

                if (kl .ge. 0) then
                    write (6, *) ' Fixing the parameters FOREVER '
                else
                    write (6, *) ' Fixing the parameters EACH TIME '
                end if

                !          if((abs(kl).eq.6.or.abs(kl).eq.7).and.ncg.eq.0) then
                !             write(6,*) ' Option kl=6,7 must be used with ncg>0 !!! '
                !             iflagerr=1
                !             iflagerrall=iflagerr+iflagerrall
                !          endif

            end if

            if (npsamax .le. 1) npsamax = 1

            if (itestr3 .eq. -8 .or. itestr3 .eq. -4 .or. itestr3 .eq. -5 .or. itestr3 .eq. -9) then
                nbin = nbinr*nw/nproc
                write (6, *) ' Total number of bins per processor =', nbin
            end if

            if (itestr3 .eq. -8 .or. itestr3 .eq. -4 .or. itestr3 .eq. -5 .or. itestr3 .eq. -9) then
                np = iese + abs(iesinv) + iesm + iesd + isfix + abs(iesfree) + abs(iessw) + iesup&
                     &+ ieskint + isfix
                write (6, *) ' Collective + normal  parameters '
                write (6, *) ' np read =', np
            else
                ! np=nbinr  ! ok what is read in readio
                ! no use of binning is done
                nbin = 0
            end if
        end if ! endif iflagerr > 0
        iflagerr = max(iflagerrall, iflagerr)

        ! reading k-points card
        kp_type = 0
        k1 = 0
        k2 = 0
        k3 = 0
        nk = 1
        nk1 = -1
        nk2 = -1
        nk3 = -1
        time_reversal = .false.
        skip_equivalence = .true.
        double_kpgrid = .false.
        ! trivial checks
        ! complex w.f. read with read_fort10_fast
        !    if(manyfort10.and.nbead.le.1..and..not.yes_complex) then
        !       write(6,*) ' ERROR: complex w.f. needed for k-points sampling'
        !       manyfort10=.false.
        !       iflagerr=1
        !    endif
        ! periodic w.f.
        !    if(manyfort10.and.nbead.le.1.and..not.iespbc) then
        !       write(6,*) ' Warning: k-points card ignored for open systems '
        !       manyfort10=.false.
        !       kp_type=0
        !    endif

        !
        if (manyfort10) then

            if (nbead .gt. 1) then
                nk = nbead

            else

                iflagerr = 1
                write (6, *) ' before kpoints '
                read (5, nml=kpoints, err=121)
                iflagerr = 0
                if (nk2 .eq. -1) nk2 = nk1
                if (nk3 .eq. -1) nk3 = nk1
                !
                ! initialization main quantities related to k-points
                select case (kp_type)

                case (0) ! gamma point calculation.
                    nk = 1
                    if (.not. gamma_point) gamma_point = .true.

                case (1, -1) ! MP grid of k-points. No KPOINTS section is needed.
                    nk = nk1*nk2*nk3
                    if (nk .le. 0) then
                        write (6, *) 'Check k-points grid !!'
                        iflagerr = 1
                    end if

                    ! For kp_type=2 and kp_type=3 one needs to specify KPOINTS section which is organized as follows:
                    !
                    ! KPOINTS
                    ! xkp(1,1) xkp(2,1) xkp(3,1) wkp(1)
                    ! xkp(1,2) xkp(2,2) xkp(3,2) wkp(2)
                    ! ......
                    ! ......
                    ! xkp(1,nk) xkp(2,nk) xkp(3,nk) wkp(nk)
                    !
                    ! In the case of kp_type=3, the points represent initial and final points
                    ! of the lines which determines the k-points path for band structure calculations.
                    !
                    !
                case (2, -2) ! k-points given in crystal coordinates in section KPOINTS.
                    ! nk1 = number of k-points
                    nk = nk1
                    if (nk .le. 0) then
                        write (6, *) 'Specify nk1 greater than 0 if using kp_type=2 !!'
                        iflagerr = 1
                    end if

                case (3, -3) ! k-points path given in crystal coordinates.
                    ! nk1 = number of points that select lines in reciprocal space (nlines-1).
                    !       Points are given in the section KPOINTS. Weigths are disregarded.
                    ! nk2 = number of k-points per line
                    nk = (nk1 - 1)*nk2 + 1
                    if (nk .le. 0) then
                        write (6, *) 'Specify nk1 and nk2 greater than 0 if using kp_type=3 !!'
                        iflagerr = 1
                    end if

                case (4, -4) ! random k-points generation
                    nk = nk1
                    if (nk .le. 0) then
                        write (6, *) 'Specify nk1 greater than 0 if using kp_type=4 !!'
                        iflagerr = 1
                    end if

                case (5, -5) ! flavor-twisted boundary conditions
                    nk = (nk1*nk2*nk3)**2
                    if (.not. double_kpgrid) double_kpgrid = .true.
                    if (nk .le. 0) then
                        write (6, *) 'Check k-points grid !!'
                        iflagerr = 1
                    end if

                case default
                    write (6, *) 'kp_type not recognized !!'
                    iflagerr = 1

                end select
                !
                write (6, *) ' after kpoints '
121             if (iflagerr .ne. 0) write (6, *) 'ERROR reading k-points'
                if (iflagerr .eq. 0) write (6, *) 'after reading k-points, kp_type =', kp_type
                iflagerrall = iflagerr + iflagerrall

            end if

        end if

    end if ! endif rank=0.

    ! final error check
    call checkiflagerr(iflagerr, rank, "ERROR in read_datasmin")

    if (itest .ne. 2) then
        change_tstep = .false.
        change_epscut = .false.
    end if

#ifdef PARALLEL
    !
    ! Broadcasting the values read by the master
    !
    ! by E. Coccia (20/12/11)
    call mpi_bcast(write_rwalk, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    ! by E. Coccia (7/12/11)
    call mpi_bcast(mm_restr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    ! by E. Coccia (22/11/10) and (3/2/11)
    call mpi_bcast(ext_pot, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vdw, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    ! by E. Coccia (9/5/11)
    call mpi_bcast(link_atom, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(calpha, maxcap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(wherescratch, 60 + lchlen, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(writescratch, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(freqcheck, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ifreqdump, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(unreliable, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(etry, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rsignr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesfast, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ip_reshuff, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idyn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(molopt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iskipdyn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbead, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nrep_bead, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(stepcg_recount, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(write_cov, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(temp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD       &
         &, ierr)
    call mpi_bcast(max_target, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(max_targetsr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(friction, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD   &
         &, ierr)
    call mpi_bcast(cutreg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD       &
         &, ierr)
    call mpi_bcast(pressfixed, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD &
         &, ierr)
    call mpi_bcast(maxdev_dyn, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD &
         &, ierr)
    call mpi_bcast(scalecov, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD   &
         &, ierr)
    call mpi_bcast(delta0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(delta0q, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(delta0k, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(scalepulay, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(scaleeloc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(powerwarp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(beta_learning, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD&
         &, ierr)
    call mpi_bcast(epsder, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD&
         &, ierr)
    call mpi_bcast(nintpsa, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(typedyncell, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npsa, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npsar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ngen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(min_block, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nscra, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(def_nscra, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(testderiv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbra, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npow, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nweight, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nmore_force, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nfat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iboot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iseedr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iopt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nproc_diag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(typereg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iread, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(beta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ncg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(itestr3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(itestr4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(developer, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(tbra, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epscutdmc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD   &
         &, ierr)
    call mpi_bcast(epstldmc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD   &
         &, ierr)
    call mpi_bcast(gamma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD      &
         &, ierr)
    call mpi_bcast(zmin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(np, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbinr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npbra, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(optbra, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nbra_cyrus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsdgel, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
         &, ierr)
    call mpi_bcast(epstion, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
         &, ierr)
    call mpi_bcast(tolcg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
         &, ierr)
    call mpi_bcast(tstep, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(tstepfn, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(tpar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD       &
         &, ierr)
    call mpi_bcast(hopfraction, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD&
         &, ierr)
    call mpi_bcast(epscut, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(change_epscut, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(change_tstep, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(shift, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(alat2v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(epstl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epstlrat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(kl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epscuttype, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(theta_reg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iese, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ieser, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npsamax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(isfix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(np3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesfree, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesinv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iessw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesup, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ieskin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesking, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ieskint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ieskinr_pos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ieskinrp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(io_level, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ipc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(parr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(power, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(parcut, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
         &, ierr)
    call mpi_bcast(parcutmin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD  &
         &, ierr)
    call mpi_bcast(parcutpar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD  &
         &, ierr)
    call mpi_bcast(parcutg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD  &
         &, ierr)
    call mpi_bcast(novar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(minz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(maxz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(minzj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(maxzj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(minjonetwobody, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(plat, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(alat2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(add_pulay, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ncore, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rs, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(kSq, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(kappar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
         &, ierr)
    call mpi_bcast(neigh, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(allfit, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(normcorr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(initpar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(powermin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npower, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nparsw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(initparsw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nparinv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(initparinv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(powerminsz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npowersz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rmaxj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rmaxinv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(smoothcut, 1, MPI_DOUBLE_PRECISION, 0                   &
         &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(killcut, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(scale_mass, 1, MPI_DOUBLE_PRECISION, 0                   &
         &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(parr_min, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(parr_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(delay_changeparr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(maxiter_changeparr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(equil_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! added by K. Nakano
    call mpi_bcast(yes_adams, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(eps_umrigar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    ! added by K. Nakano
    call mpi_bcast(max_ortho, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(prep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(rejweight, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(better_dmc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(gauge_fixing, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(safelrdmc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(noblocking, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(add_diff, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesalfe, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(pseudorandom, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yeszj, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yeszagp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(srcomplex, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(membig, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(membigcpu, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(symiesup, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(onebodysz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(twobodyoff, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesdtwobodyoff, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesdonebodyoff, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yespress, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(warp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yespulay, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_correct, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(k6gen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsbas, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(eps_dyn5, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(maxtime, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(epsvar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yescutdet, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yescutjas, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(oldscaling, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(defparcutg, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(fixpar, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(symmetrize_agp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesquantum, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesavopt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesavcov, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesavsr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yeswrite10, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesread10, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yeswritebead, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesperiodize, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesturboq, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yessecond, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yesrootc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(addrognoso, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(cleanrognoso, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(eqcellab, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(eqcellac, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(eqcellbc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(fixa, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(fixb, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(fixc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(real_contracted, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(real_agp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(no_sjbra, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_complex, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(changelambda, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    ! --------------------------- k-points -------------------------------
    call mpi_bcast(skip_equivalence, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(time_reversal, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(double_kpgrid, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_kpoints, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(manyfort10, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(kp_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(k1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(k2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(k3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nk1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nk2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nk3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(decoupled_run, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(compute_bands, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    ! --------------------------------------------------------------------
    ! broadcast iespbc flag needed by various subroutines (ex. open_files)
    call mpi_bcast(iespbc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    !  molyes defined in read_fort10_fast
    call mpi_bcast(molyes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(signalnoise, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(noopt_onebody, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(scalermax, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ireadminr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(double_mesh, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(change_parr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(default_epsdgel, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(read_molecul, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(lrdmc_der, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(lrdmc_nonodes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(enforce_detailb, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iesrandoma, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nowrite12, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(flush_write, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_fastbranch, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(true_wagner, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! For module dielectric
    call mpi_bcast(dielectric_length, 1, MPI_DOUBLE_PRECISION, 0&
         &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(dielectric_ratio, 1, MPI_DOUBLE_PRECISION, 0&
         &, MPI_COMM_WORLD, ierr)
    call mpi_bcast(case_diel, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(weight_moroni, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_scemama, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_scemama_open, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(novec_loop1, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_sparse, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_sparse_choose, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(max_sparse_choice, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(yes_dgelscut, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
    ! For unreliable  networks.
    call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif

    min_running_ave = 1d20
    min_running_std = 1d20
    min_running_ave_energy = 1.d20
    min_running_std_energy = 1.d20

    if (max_targetsr .eq. -1) then
        max_targetsr = max_target
#ifdef _OFFLOAD
        if (rank .eq. 0) write (6, *) ' Default value of max_targetsr= ', max_targetsr
#endif
    end if
    ! # walkers must be a multiple of # processors
    if (mod(nw, nproc) .ne. 0 .or. nproc .gt. 1000000) then
        if (rank .eq. 0) then
            write (*, *) 'change number of processors!!!'
            write (*, *) 'walkers', nw, ' processes', nproc
            if (nproc .gt. 1000000) then
                write (6, *) ' too many processors, change convertdec'
                write (6, *) ' with more than 6 digits  '
            end if
        end if
        call mpi_finalize(ierr)
        stop
    end if

#endif

    call init_dielectric

    if (enforce_detailb .or. alat .lt. 0.d0) then
        if (.not. iesrandoma .and. rank .eq. 0) write (6, *) ' Warning random lattice forced (iesrandoma=.true.) '
        iesrandoma = .true.
    end if
    !
    if (ipc .eq. 1 .and. srcomplex) then
        if (rank .eq. 0) write (6, *) ' Warning srcomplex turned to false (real case) '
        srcomplex = .false.
    end if
    !

    ! initialize k-point grid for all processors
    !
    allocate (xkp(3, nk), wkp(nk))
    allocate (xkp_down(3, nk), wkp_down(nk))
    xkp = 0.d0
    wkp = 1.d0
    xkp_down = 0.d0
    wkp_down = 1.d0
    kaverage = .false.
    tot_wt = 0.d0
    tot_wt_down = 0.d0

    if (rank .eq. 0) then
        ! k-points calculation
        if (manyfort10 .and. nbead .le. 1) then

            call get_kpoints(iflagerr, nel, nion, rion_fast, atom_number_fast, rs)

            ! check if k-points are normalized, otherwise rescale to obtain
            ! a correct normalization
            do i = 1, nk
                tot_wt = tot_wt + wkp(i)
                tot_wt_down = tot_wt_down + wkp_down(i)
            end do
            if (abs(tot_wt - 1.d0) .gt. 1d-6) then
                write (6, *) ' Warning: wrong k-points up weights,rescaling! '
                wkp(:) = 1.d0/nk
                tot_wt = 1.d0
            end if
            if (abs(tot_wt_down - 1.d0) .gt. 1d-6) then
                write (6, *) ' Warning: wrong k-points down weights,rescaling! '
                wkp_down(:) = 1.d0/nk
                tot_wt_down = 1.d0
            end if

            if (double_kpgrid .and. rank .eq. 0) &
                write (6, *) 'Warning: up spin k-points might be different from down spin ones!'
            write (6, *) ' '
            write (6, *) ' type k-points/# k-points/total weight ', kp_type, nk
            write (6, *) ' k-points up/weights:'
            do i = 1, nk
                write (6, 300) i, xkp(1, i), xkp(2, i), xkp(3, i), wkp(i)
            end do
            write (6, *) ' '
            write (6, *) ' k-points down/weights:'
            do i = 1, nk
                write (6, 300) i, xkp_down(1, i), xkp_down(2, i), xkp_down(3, i), wkp_down(i)
            end do
            write (6, *) ' '
300         format(3x, I6, 4x, F10.7, 3x, F10.7, 3x, F10.7, 3x, F10.7)
            ! single phase calculation, xkp not relevant
        elseif (iespbc) then
            xkp(:, 1) = phase(:)
            xkp_down(:, 1) = phase_down(:)
            wkp(1) = 1.d0
            wkp_down(1) = 1.d0
        end if

        ! turning on k-points calculation.
        if (manyfort10 .and. nbead .le. 1) kaverage = .true.
        ! if no k-sampling, no sense to perform decoupled runs.
        if (.not. kaverage .and. decoupled_run) decoupled_run = .false.

        if (kaverage .and. nbead .le. 1) then
            ! write file "kp_weights.dat" needed to perform averages
            open (unit=37, file='kp_info.dat', form='formatted', status='unknown', position='rewind')
            write (6, *) ' Writing k-points information on file '
            write (37, '(I6)') nk
            write (37, *) '# up spin electrons '
            do i = 1, nk
                write (37, 301) i, xkp(1, i), xkp(2, i), xkp(3, i), wkp(i)
            end do
            write (37, *) '# down spin electrons '
            do i = 1, nk
                write (37, 301) i, xkp_down(1, i), xkp_down(2, i), xkp_down(3, i), wkp_down(i)
            end do
            close (37)
301         format(I6, 3f13.8, 1f19.14)
        end if

    end if

    call checkiflagerr(iflagerr, rank, "ERROR in k-points generation")

    ! broadcast k-points generated from the master
#ifdef PARALLEL
    call mpi_bcast(xkp, size(xkp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(wkp, size(wkp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(xkp_down, size(xkp_down), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(wkp_down, size(wkp_down), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(tot_wt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(tot_wt_down, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(kaverage, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(decoupled_run, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    ! broadcast again nk, it can change if # of k-points is reduced by BL symmetries.
    call mpi_bcast(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! broadcast again iread, it can change if decoupled_run is active.
    call mpi_bcast(iread, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

    if (kaverage .and. .not. decoupled_run .and. .not. yesavsr) then
        if (rank .eq. 0) write (6, *) ' Warning yesavsr true in this case '
        yesavsr = .true.
    end if

    if (kaverage .and. .not. decoupled_run .and. yesavcov) then
        if (rank .eq. 0) write (6, *) ' Warning yesavcov false in this case '
        yesavcov = .false.
    end if

    if (kaverage .and. .not. decoupled_run .and. yesavopt .and. molyes .and. yesavsr) then
        if (rank .eq. 0) write (6, *) ' Warning yesavopt false in this case '
        yesavopt = .false.
    end if

    if (rank .eq. 0) then
        write (6, *)
        if (kaverage .and. .not. decoupled_run) then
            write (6, *) ' Warning: starting a twist-averaged calculation!'
        elseif (kaverage .and. decoupled_run) then
            write (6, *) ' Warning: starting a decoupled k-points calculation!'
        end if
    end if
    !
    !
    ncg_adr = npar + nparsw + nparinv
    kp0 = np - ieskin

    if (idyn .eq. 0) then
        tparf = tpar
    else
        tparf = 1.d0
    end if
    if (rank .eq. 0) write (6, *) ' tparf =', tparf

    if (trim(wherescratch) .eq. 'old') then
        oldscra = .true.
    else
        oldscra = .false.
    end if

    if (itestr4 .ge. -10) then
        fncont = .false.
    elseif (itestr4 .ge. -20) then
        fncont = .true.
    elseif (itestr4 .ge. -30) then
        fncont = .false.
    end if

    nintpseudo = nintpsa

    ! # of walkers per processor
    in1 = nw/nproc

    ieskindim = max(ieskin, 1)
    ieskingdim = max(iesking, 1)
    nbindim = max(nbin, 1)
    ncgdim = max(ncg, 1)
    npdim = max(np, 1)

    nbinmax = ncg + npbra

    indinv = 0
    endinv = iesinv
    indfree = iesm + iesinv + iesd + iese
    indsw = indfree + iesfree

    !          itestr4 == actual input containing info about FN/LRDMC/VMC optimization
    !          according to its range modulo 10
    !          itestr3  ==  but without this info
    !          itestrr == determines the type of optimization
    !          as input in reweight0 itestrr=-5 --> SR
    !                                itestrr=-4 --> Hessian
    !          itestr=-5  means we are doing an optimization
    !          itest=2  uses VMC in all methods, itest=1 uses FN/LRDMC

    itestrr = itestr3
    if (itestr3 .eq. -8) itestrr = -4
    if (itestr3 .eq. -9) itestrr = -5
    itestr = itestrr
    if (itestr .eq. -4) itestr = -5
    itest = abs(itestr)
    if (itestr .eq. -5) then
        if (itestr4 .ge. -10) then
            itest = 2
        else
            itest = 1
        end if
    end if

    if (abs(itestr) .eq. 6 .or. itestr .eq. -2 .or. itestr .eq. -3) itest = 1
    if (itestr .eq. -7) itest = 2

    if (abs(itestr) .eq. 1 .or. itestr .eq. 6 .or. itestr .eq. -2 .or. itestr .eq. -3) fncont = .true.

    !          definition yesnleft

    if (fncont) then ! standard dmc
        yesnleft = .false.
    elseif (itest .ne. 1) then !standard vmc
        yesnleft = .true.
    elseif (itest .eq. 1 .and. nbra .ne. 1) then ! lrdmc with fixed # of acc.
        yesnleft = .true.
        if (yes_fastbranch .and. def_nscra) then
            nscra = 2*nbra
            if (rank .eq. 0) write (6, *) ' Warning upscratch done only after branching !!! '
        end if
    else
        yesnleft = .false. ! standard lrdmc
    end if

    if (true_wagner .eq. -1 .and. lrdmc_der .and. .not. lrdmc_nonodes) then
        true_wagner = 2
        if (rank .eq. 0) write (6, *) ' Default value of true_wagner =', true_wagner
    end if
    if (yesnleft .and. cutweight .eq. -1.d6 .and. rank .eq. 0 .and. true_wagner .le. 0) then
        if (idyn .ne. 0) then
            cutweight = 4.d0*zmax
        else
            cutweight = 0.d0
        end if
        write (6, *) 'Default cutoff on the weight =', cutweight
    end if
    if (true_wagner .gt. 0) then
        if (cutweight .gt. 0.d0) then
            if (lrdmc_der .and. .not. lrdmc_nonodes) then
                cutweight = cutweight/abs(alat)
            else
                cutweight = cutweight*abs(alat)**(2.d0/3.d0)
            end if
            if (rank .eq. 0) write (6, *) ' Warning Wagner regularization eps (scaled by alat^2/3)  = ', cutweight
        elseif (cutweight .ne. -1.d6) then
            cutweight = -cutweight
            if (rank .eq. 0) write (6, *) ' Warning Wagner regularization eps = ', cutweight
        else
            cutweight = 0.d0
        end if
    end if
#ifdef PARALLEL
    call mpi_bcast(cutweight, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD       &
         &, ierr)
#endif

    if (itestr4 .ge. -10) then
        itestrfn = itestr
    elseif (itestr4 .ge. -20) then
        itestrfn = 1
    elseif (itestr4 .ge. -30) then
        itestrfn = -6
    end if

    if (itest .ne. 2) then
        !       No branching with VMC
        iesbra = .true.
    else
        iesbra = .false.
    end if

    if (rank .eq. 0) write (6, *) 'itestrfn =', itestrfn

    ! Sandro, please update this comment. I think there are opnions no longer used
    ! see itestr=-1 (I guess..)

    !          itestr=2   Standard VMC
    !          itestr=-5  SR full minimizator with metropolis
    !          itestr=-1  FN lattice Reg.  with upscratch
    !          itestr=1   FN  standard   with Dt=tbra (rsignr=0, tstep=0, ga
    !          itestr=1   VMC with Dt=tbra (rsignr=0, tstep=1, gamma=0)
    !          itestr=-6  FN  lattice Reg. fast
    !          itestr=6   FN LR with locality approximation
    !          itestr=-7  VMC with kinetic energy discretized
    !          itestr=-4  SR full minimizator with metropolis  + Hessian
    !          itestr=-8  SRH full minimizator + zeta minimization
    !          itestr=-9  SR minimization with no Z
    !          itestr=-2  FN DMC with non local moves (heat bath after an all electron diffusion)
    !          itestr=-3  FN DMC with non local moves (heat bath after single electron diffusion)

    ! Fixed at machine precision
    epst = 1d-11
    if (kaverage .and. decoupled_run) then
        if (iread .gt. 0 .and. (epscut .ne. 0 .or. itest .ne. 2)) then
            if (rank .eq. 0) write (6, *) ' Warning iread=0 in this case !!! '
        end if
    end if

    if (ncore .lt. 0) then

        ncore = abs(ncore)
        allocate (icore(np), tcore(np))
        tcore = 0.d0
        icore = 0

        if (rank .eq. 0) then
            write (6, *) '----------------------------------------------'
            write (6, *) ' Used ', ncore + 1, 'different energy scales'

            ind = 0

            if (tion .eq. -1.d0) then

                do i = 1, ncore

                    read (5, *) nmp, tcore(ind + 1), (icore(j + ind), j=1, max(nmp, 1))
                    if (nmp .gt. 0) then
                        do j = 2, nmp
                            tcore(ind + j) = tcore(ind + 1)
                        end do
                    else
                        nmp = -nmp
                        do j = 2, nmp
                            tcore(ind + j) = tcore(ind + 1)
                            icore(ind + j) = icore(ind + 1) + j - 1
                        end do
                    end if

                    ind = ind + nmp

                end do

                ncore = ind

                write (6, *) ' Core parameters =', (icore(i), i=1, ncore)

            else
                ncore = ieskin
                do i = 1, ieskin
                    icore(i) = np - ieskin + i
                    tcore(i) = tion
                end do

                if (tcell .eq. 0.d0 .and. typedyncell .ne. 0) then
                    write (6, *) ' Warning default tcell=tion ', tion
                    tcell = tion
                end if

                if (typedyncell .eq. 1) then
                    do i = ieskin - 1, ieskin
                        tcore(i) = tcell
                    end do
                    if (fixc) tcore(ieskin) = 0.d0 ! Fixed surface changing b a=V/bxc
                elseif (typedyncell .eq. 2) then
                    do i = ieskin - 2, ieskin
                        tcore(i) = tcell
                    end do
                    if (fixa) tcore(ieskin - 2) = 0.d0
                    if (fixb) tcore(ieskin - 1) = 0.d0
                    if (fixc) tcore(ieskin) = 0.d0
                elseif (typedyncell .eq. 3) then
                    tcore(ieskin) = tcell
                end if
            end if
            ! rank.eq.0
        end if

#ifdef PARALLEL
        call mpi_bcast(ncore, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(tcell, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(tion, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
        !   For unreliable  networks.
        call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
        call bcast_real(tcore, ncore, 0, MPI_COMM_WORLD)
        call bcast_integer(icore, ncore, 0, MPI_COMM_WORLD)
#ifdef UNREL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
#endif
    else
        allocate (icore(1))
        allocate (tcore(1))
        tcore(1) = 0.d0
        icore = 0
    end if
end subroutine read_datasmin

subroutine read_datasmin_mol
    use allio
    use convertmod, only: nmolmatdo
    implicit none

#ifdef PARALLEL
    include 'mpif.h'
#endif
    if (itestr .eq. -5 .or. read_molecul) then
!

        if (rank .eq. 0) then
!       default values
            epsdgm = 1d-14 !  No more than 12 digits in diagonalization.
            smearing = 1d-5
            nbufd = -1
!       nmol=molecular-(nelup-neldo)
            nmolmax = 0
            nmolmaxw = 0
            nmolmin = 0
            nx = 0
            ny = 0
            nz = 0
            ax = 0.d0
            ay = 0.d0
            az = 0.d0
            weight_loc = -1.d0
            orthoyes = .true.
            gramyes = .true.
            iflagerr = 1
            if (npsar .gt. 0) then
                add_onebody2det = .false.
            else
                add_onebody2det = .true.
            end if
            epsrem_contr = epsdgel
            shift_origin = .true.
            shiftx = .false.
            shifty = .false.
            shiftz = .false.
            read (5, nml=molecul, err=121)
            write (6, *) ' After reading molecul '
            if (yesavopt) then
                if (rank .eq. 0) write (6, *) ' Warning yesavopt forced to false with mol optimiz.!!! '
                yesavopt = .false.
            end if
            iflagerr = 0
            if (molecular .eq. 0) then
                write (6, *) ' Warning   fort.10 should have molecular orbitals,&
              & please run again with the output fort.10  !'
                if (nmol .eq. -1 .or. nmol .lt. neldo) then
                    iflagerr = 1
                    write (6, *) ' ERROR you should have molecular orbitals in fort.10 '
                    write (6, *) ' ERROR please use convertfort10mol or rerun with nmol>neldo '
                end if
            end if
            if (contraction .eq. 0) then
                iflagerr = 1
                write (6, *) ' ERROR you should have contracted orbitals in fort.10 '
                write (6, *) ' ERROR please introduce contraction (even fake) in your AGP '
            end if

121         if (iflagerr .ne. 0) then
                write (6, *) ' ERROR reading molecul '
                iflagerrall = iflagerr + iflagerrall
            else
                if (ny .eq. 0) then
                    ny = nx
                    write (6, *) ' Default value for ny=', ny
                end if
                if (nz .eq. 0) then
                    nz = ny
                    write (6, *) ' Default value for nz=', nz
                end if
                if (.not. iespbc) then
                    if (ay .eq. 0.d0) then
                        ay = ax
                        write (6, *) ' Default value for ay=', ay
                    end if
                    if (az .eq. 0.d0) then
                        az = ay
                        write (6, *) ' Default value for az=', az
                    end if
                end if
                if (symmagp .and. ipc .eq. 1) then
                    nmol = molecular - ndiff
                else
                    nmol = (molecular - ndiff)/2
                end if
                write (6, *) ' Default value of nmol ', nmol
                if (nmolmin .eq. 0) then
                    nmolmin = neldo
                    write (6, *) ' Default value of nmolmin ', nmolmin
                end if
                if (nmolmax .eq. 0) then
                    nmolmax = neldo
                    write (6, *) ' Default value of nmolmax ', nmolmax
                end if

                write (6, *) ' after  read molec '

                if (weight_loc .eq. 0.d0) then
                    if (epsdgm .ne. 0.d0) then
                        weight_loc = epsdgm**2
                    else
                        weight_loc = 1d-8
                    end if
                    write (6, *) ' Default value for weight_loc =', weight_loc
                end if

                if (nmol .ne. 0 .and. nmolmax .eq. 0) then
                    nmolmax = nmol
                    write (6, *) ' Default value for nmolmax =', nmolmax
                end if

                if (nmolmaxw .eq. 0) then
                    nmolmaxw = nmolmax
                    write (6, *) ' Default value of nmolmaxw=', nmolmaxw
                end if

                if (nmolmax .lt. neldo) then
                    iflagerrall = iflagerrall + 1
                    write (6, *) ' Too small  nmolmax> =', neldo
                end if

!       read(5,*) epsdgm
                write (6, *) ' error converter  ', epsdgm
!       read(5,*) nx,ny,nz
                write (6, *) ' # mesh read ', nx, ny, nz

                if (nbufd .eq. -1) then
#ifdef __SCALAPACK
                    if (nelorb .gt. 2000 .and. .not. yesdft) then ! with SCALAPACK no memory problem
#else
                        if (nelorb .gt. 2000) then
#endif
                            nbufd = 100 ! there may be memory problems
                        else
                            nbufd = 1000 ! almost maximum efficiency dgemm
                        end if
                        write (6, *) 'Default value for buffer =', nbufd
                    end if

                end if

            end if ! endif rank.eq.0

!       never print the overlap in this case
            printoverlap = .false.

#ifdef PARALLEL
            call mpi_bcast(epsdgm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(smearing, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(epsrem_contr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(weight_loc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nbufd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(orthoyes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(gramyes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(add_onebody2det, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(shift_origin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(shiftx, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(shifty, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(shiftz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(iflagerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nw_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
!   For unreliable  networks.
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
#endif

            if (iessw .eq. 0 .and. iesup .eq. 0 .and. .not. yesdft .and. .not. read_molecul) then
                call checkiflagerr(1, rank,&
            &' ERROR you cannot use molopt>0 without optimizing AGP !!!')
            end if

            if (nx .eq. 0 .or. ny .eq. 0 .or. nz .eq. 0) then
                write (errmsg, *) ' Mesh should be finite !!!', nx, ny, nz
                call checkiflagerr(1, rank, errmsg)
            end if

            if (gramyes .and. .not. orthoyes) then
                if (rank .eq. 0) write (6, *) &
         &' Warning with Gram-Schmidt ortho, changed orthoyes= true'
                orthoyes = .true.
            end if
            if (iespbc) then
                ax = cellscale(1)/nx
                ay = cellscale(2)/ny
                az = cellscale(3)/nz
                if (rank .eq. 0)                                                  &
            &write (6, *) ' lattice mesh chosen ', ax, ay, az, cellscale(1)
            else
                if (rank .eq. 0) write (6, *) ' lattice mesh read ax,ay,az ', ax, ay, az
!       if(rank.eq.0) read(5,*) ax,ay,az
#ifdef PARALLEL
                call mpi_bcast(ax, 1, MPI_DOUBLE_PRECISION                        &
             &, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(ay, 1, MPI_DOUBLE_PRECISION                        &
             &, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(az, 1, MPI_DOUBLE_PRECISION                        &
             &, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
!   For unreliable  networks.
                call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
#endif
            end if
            if (ax .eq. 0 .or. ay .eq. 0 .or. az .eq. 0) then
                write (errmsg, *) ' Lattice constants should be finite !!!', ax, ay, az
                call checkiflagerr(1, rank, errmsg)
            end if

            if (rank .eq. 0) then
                write (6, *) '# molecular orbital Det considered/projected'
                write (6, *) nmol, nmolmin, nmolmax
            end if

#ifdef PARALLEL
            call mpi_bcast(nmol, 1, MPI_INTEGER                               &
         &, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nmolmin, 1, MPI_INTEGER                            &
         &, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nmolmax, 1, MPI_INTEGER                            &
         &, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nmolmaxw, 1, MPI_INTEGER                            &
         &, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
!   For unreliable  networks.
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
#endif
            detc_proj = .false.
            yesmin = 0
            if (molopt .ne. 0) then
                yesmin = 1
                if (rank .eq. 0 .and. .not. yesdft) write (6, *) ' Warning molecular orbitals&
                & are constraint to be written in terms of contracted orbitals'
                detc_proj = .true.
                molopt = 1
            end if
        else
            yesmin = 0
            detc_proj = .false.
        end if ! closed main if itestr.eq.-5

!  Here we should interchange nmol_min with nmolmin because they have the
!  opposite meaning in the code. In the code nmolmin is used for projection
!, whereas nmol_min set to one the eigenvalues <= nmol_min. In input the
!  meaning is opposite.
        if (yesmin .eq. 1) then
            nmolmat = ipf*nmolmax
            nmolmatw = nmolmaxw
            if (symmagp .or. ipf .eq. 2) then
                nmolmatdo = nmolmat
            else
                nmolmatdo = nmolmax
            end if
        else
            nmolmatdo = 0
            nmolmat = 0
            nmolmatw = 0
        end if
        end subroutine read_datasmin_mol
