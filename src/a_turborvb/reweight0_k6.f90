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

subroutine reweight0(nw, in1, npr, npm, factorsr                 &
        &, ipip, psip, alpha, sov, econf, econfh, econfion                  &
        &, ieskin, ncg, epst, epstion, iweight, wcort, itest, enert          &
        &, scalpar, parcutr, iflagerr, forza, err, iesconv                       &
        &, itouch, parcutmin, parcutpar, numparmin, klr                         &
        &, etry, epsi, tjas, beta, ist, ien, rank, ierr, parr, parcute       &
        &, nbin, fk, dimfk, fkav, okav, skdiag, weightall, reduce, nbinmax, nweight        &
        &, ibinit, indopen3, stepcg, lwork, idyn, temp, weight_vir, friction, scalecov  &
        &, delta0, delta0q, delta0k, dtr, velion, ris, tmes, cov, rpar, npar, initpar, nparsw   &
        &, initparsw, nparinv, initpowerinv, endinv                            &
        &, iond, nion, adrlambda, rmax, indopen4, writescratch, countscra, bufscra &
        &, jas_invariant, tolcg, ieskinion, rion, iespbc, atom_number     &
        &, eps_dyn5, maxdev_dyn, acc_dyn, normcorr, row_comm, row_id, yescomm)
    use constants, only: ipc, TWO_PI
    use allio, only: reducel, ncg_adr, kp0, stepcg_recount, write_cov, yesquantum, nbead&
            &, commopt_mpi, nprocopt, rankopt, nproc, yesavcov, commcov_mpi, nproccov, commsr_mpi&
            &, nprocsr, ranksr, ndimj, ndimjp, yesavopt, col_id, yesavsr, col_comm, mcol, commcolrep_mpi&
            &, mcol_rep, commrep_mpi, commcolsr_mpi, nrep_bead, rankrep, rankcolrep, nprocrep&
            &, kdyn, kdyn_eig, yesturboq, yessecond, ndims, ndimsp, symmagp, yes_hermite, allowed_par&
            &, yes_real, srcomplex, yes_correct, kaverage, yesrootc, cellpi&
            &, addrognoso, cov_old, ion_table, ieskinr, cleanrognoso, signalnoise&
            &, tion, tcell, ieskinr_pos, yes_hessc, wkp, ikpoint, commcolrep_mpi, npbra&
            &, real_agp, cellscale, change_parr, default_epsdgel, parr_min, parr_max&
            &, delay_changeparr, iesdelay, min_block, ortho_cntx, ortho_comm, k6gen, ngentry&
            &, maxiter_changeparr, prep, comm_raw, comm_col, rankcol, rankraw, beta_learning&
            &, eps_umrigar, max_targetsr, scale_grad, norm_corr, change_tpar, yes_dgelscut, nproc_diag
    use scal_lins
    implicit none
    integer Nw, nion, nwr, np, npp, nps, npps, npsr, i, j, ipip(*), k, info, k_par       &
            &, ncg, npm, iweight, iflagerr, iesconv, kk, kkk                        &
            &, ind, itest, kp_ion, kp_complex, ieskin, numpar, indopen3, indopen4     &
            &, itouch(*), numparmin, jmax, npr, nweight, ibinit                    &
            &, nptot, n2, n3, n4, n5, maxnm, jj, ii, nbinmax, dimmat                   &
            &, klr, indi, indj, nbin, nbinp, nbinr, lwork, lworkd, indr      &
            &, nbintot, numparp, indscali, indscalj, indmax, indmin, stepcg           &
            &, kkmin, in1, irep, iterk, irest, ndsr, maxit, imin, indin, indjn    &
            &, idyn, ieskin2, nprocm, blockdim, iblock, jblock, block_r           &
            &, blockread, irepr, dimbuf, npar, initpar, nparsw, initparsw, ncoll       &
            &, initpowerinv, nparinv, endinv, kn_ion, kn_ionp, ku_ion, ku_ionp        &
            &, adrlambda(2, *), indscra, ieskinion, mineig, niont, npk, dimfk, row_comm, row_id&
            &, nprepm, mini, maxi
    real*8 cclock, timep, beta_used, time_load
    real*8 wcort, weightall, normcorr, maxsr, condnum, Picost, weight_vir&
            &, psip(*), alpha(*), alph, cost1, dtstep, dtnoise, eigS&
            &, factorsr(*), sov(npm, npm, *), cost, costw, costwe, epst, parr, parcute&
            &, econf(in1, *), econfh(in1, *), econfion(in1, *), error        &
            &, parcutmin, parcutpar, devmax, enert(ipc, *), devmaxp, devinv, devmaxc       &
            &, scalpar(*), parcut, forza(*), err(*), tjas, etry, epsi, eps_dyn5 &
            &, costmy, corrn2, beta, betap, dnrm2, epsin, cut, cut_umrigar, argtan         &
            &, fk(dimfk, *), cost0, reduce(max(ncg, 1), *), scali, scalj                &
            &, fkav(*), okav(*), skdiag(*), costtot, parcutr, time                   &
            &, parcut2, dummy, eigmin, costn, cnorm, ddot, velion(3, *)                &
            &, temp, zeta, drand1, pi, dt, friction, scalecov, ris(4), Tmes, dtr   &
            &, tolcg, summpi, cov(*), rpar(*), force_coll, delta0, delta0q, delta0k, epstion&
            &, iond(nion, *), rmax, nzsr, nzsro, jas_invariant(*), rioncm(3)          &
            &, rion(3, nion), atom_number(*), over2, maxdev_dyn, cov_sav(ieskin*ieskin)
    real*8 mnoise(2, 2), zetan(2), costh, lnoise(2), eta_v, eta_r        &
            &, Tn, Gn, Gni(2, 2), alphaqmc, alphaall, scalef, maxsn&
            &, dth, Gnh, Gnih, normsr, costf, costc
    real(8) bufscra(*)
    integer errnoise, writescratch, countscra, nbufscra, n3n
    real(8), dimension(:, :), allocatable :: mat, mat_prep, eb, econf_coll&
            &, riondcm, vecrottr, projrottr, covpurif, psipbead, mat_buf
    real(8), dimension(:), allocatable :: v1, v2, v3, v4, sr, srg, srb, work&
            &, okav_coll, eig, alfa, betal, forza_sav, workGS, vel0, buff1, buff2&
            &, factorsr_sav
    integer(4), dimension(:), allocatable :: index
    logical, dimension(:), allocatable :: truecomplex
    integer*1, dimension(:), allocatable :: type_complex
    real(8) sinut, sinut_over_omega, sinut_times_omega, cosit, omega_harm
    parameter(pi=3.14159265358979323846d0)
    integer rank, ist, ien, istm, ierr, block, blockl, blockr, maxdim_inv&
   &, indk, indi_save
    logical iespbc, acc_dyn, yescomm, yes_umrigar, yes_targetprep
#ifdef PARALLEL
    include 'mpif.h'
#endif

    if (iflagerr .ne. 0) return

    errnoise = 0
    np = npr
    istm = ist - 1
    npp = np + 1
    nbinp = nbinmax + 1
    nbufscra = (2*np + ipc + 1)*in1 + 1
    maxnm = max(npp, nbinp)
    if (yes_hessc) maxnm = max(maxnm, (npbra + ncg)*ipc)
    n2 = maxnm + 1
    n3 = maxnm + n2
    n4 = maxnm + n3
    n5 = maxnm + n4
    npk = np - ieskin
    kp_ion = npk
    ieskin2 = ieskin*ieskin
    nptot = npm*npm
    nbintot = max(ncg, 1)*np
    parcut = abs(parcutr)
    !       Machine precision on input

    parcut2 = parcut**2

    if (parr .lt. 0) then
        yes_umrigar = .true.
        cut = -parr
        cut_umrigar = eps_umrigar*(1.d0 + cut)
    else
        cut_umrigar = 0.d0
        cut = parr
        yes_umrigar = .false.
    end if
    !       first compute the total weights weight and wm

    if (np .ne. 0) then
        !   if(rank.eq.0)  write(6,*) ' iweight inside =',iweight
        !       now calculation average correlation functions without sign
        !       update the matrix sov with rank one date

        !       write(6,*) 'wcort inside =',wcort

        countscra = countscra + 1
        !       if(rank.eq.0) write(6,*) ' Updated countscra =',countscra
        cost = wcort
        if (writescratch .ne. 0) then
            indscra = (countscra - 1)*nbufscra + 1
            bufscra(indscra) = cost
            do i = ist, ien
                do j = 1, np
                    bufscra(indscra + 1) = econf(i - istm, j)
                    bufscra(indscra + 2) = econfh(i - istm, j)
                    indscra = indscra + 2
                end do
!               bufscra(indscra + 1) = wconfn(i)
                bufscra(indscra + 1) = factorsr(i - istm) ! to be consistent with fn optimization.
                bufscra(indscra + 2) = enert(1, i - istm)
                if (ipc .eq. 2) then
                    bufscra(indscra + 3) = enert(2, i - istm)
                    indscra = indscra + 3
                else
                    indscra = indscra + 2
                end if
            end do
        else
            write (indopen3) cost, ((econf(i - istm, j), econfh(i - istm, j), j=1, np) &
                    &, factorsr(i - istm), enert(1:ipc, i - istm), i=ist, ien)
        end if

        if (iweight .eq. 1) then

            allocate (forza_sav(npr))
            call dcopy(npr, forza, 1, forza_sav, 1)

            if (idyn .ge. 2 .and. ieskin .gt. 0) then
                allocate (vel0(ieskin))
                vel0(1:ieskin) = velion(3, 1:ieskin)
            end if
            if (writescratch .eq. 0) rewind (indopen3)
            do i = 1, np
                ipip(i + np) = 1
            end do

            do i = ncg, 2, -1
                call dcopy(np, reduce(i - 1, 1), ncg, reduce(i, 1), ncg)
            end do
            call dscalzero(ncg, 1.d0, psip, 1)
            !          the hessian is not included in the collective

            if (ncg .gt. 0 .and. beta_learning .eq. 0.d0) call dscalzero(np, 0.d0, reduce, ncg)

            !          Use the other random directions
            if (stepcg .lt. ncg) then
                stepcg = stepcg + 1
            elseif (stepcg .gt. ncg) then
                stepcg = ncg
                if (rank .eq. 0) write (6, *) ' Warning setting ncg steps = ncg', ncg
            end if

            if (rank .eq. 0) write (6, *) ' Number of cg steps=', stepcg

            ind = ncg
            indmin = ind
            indmax = nbinmax
            ipip(n3:n3 + nbinmax - 1) = 0

            call dcopy(np, forza, 1, psip(n3), 1)

            !        write(6,*) ' Ratio err/skdiag '
            !        do i=1,kp_ion
            !        if(err(i).gt.0) write(6,*) i,skdiag(i)/err(i)**2
            !        enddo

            if (ncg .gt. 0) call dgemv('N', ncg, np, 1.d0, reduce, ncg, psip(n3), 1, 0.d0, psip(n2), 1)

            !         calculation devinv  PARALLEL
            devinv = 0.d0
            if (ncg .gt. 0) then
                !         calculation devmax in the direction of the last force
                kkmin = 1
                if (abs(klr) .eq. 6 .or. abs(klr) .eq. 7 .or. abs(klr) .eq. 2) kkmin = 2
                !!      if(idyn.ge.2) then
                do kk = kkmin, ncg
                    call dgemv('N', nbin, np, 1.d0, fk, nbin, reduce(kk, 1), ncg          &
                            &, 0.d0, psip(n4), 1)
                    cost = psip(n4)**2
                    do jj = 1, nbin - 1
                        cost = cost + psip(n4 + jj)**2
                    end do

#ifdef PARALLEL
                    call mpi_allreduce(cost, summpi, 1, MPI_DOUBLE_PRECISION           &
                         &  , MPI_SUM, commsr_mpi, ierr)
                    cost = dsqrt(summpi)
#else
                    cost = dsqrt(cost)
#endif

                    !                if(rank.eq.0) then

                    if (cost .ne. 0.d0) cost = abs(psip(n2 + kk - 1))/cost
                    if (cost .gt. devinv) devinv = cost
                    !

                    if (cost .lt. parcutmin) then
                        if (rank .eq. 0) write (6, *) ' disregarded coll =', kk, cost
                        ! disregard this collective parameter
                        ipip(n3 + kk - 1) = -1
                    else
                        if (cost .ne. 0.d0 .and. rank .eq. 0) write (6, *) ' collective =', kk, cost
                    end if
                    ! endif rank=0
                    !                endif
                    !         devinv=psip(n2)/dnrm2(nbin,psip(n4),1)
                end do
                !!      if(idyn.ge.2.and.ieskin2.gt.0) then

                if (rank .eq. 0) write (6, *) ' devmax force = ', devinv
                ! endif  devinv
            end if

            !   Calculation covariance matrix

            if (idyn .ge. 2 .and. ieskin2 .gt. 0) then
                call dgemm('T', 'N', ieskin, ieskin, nbin, 1.d0, fk(1, kp_ion + 1)     &
                        &, dimfk, fk(1, kp_ion + 1), dimfk, 0.d0, cov, ieskin)
#ifdef PARALLEL
                if ((yesquantum .and. .not. yesavcov) .or. (yesturboq .and. yesavcov)) then
!     the default for quantum
                    call reduce_base_real(ieskin2, cov, mpi_comm_world, -1)
!     All beads will have the same average covariance
                    cov(1:ieskin2) = cov(1:ieskin2)/nbead
                elseif (nrep_bead .gt. 1 .and. yesquantum) then
                    call reduce_base_real(ieskin2, cov, commopt_mpi, -1)
                    cov(1:ieskin2) = cov(1:ieskin2)/nrep_bead
                elseif (kaverage) then
                    call reduce_base_real(ieskin2, cov, mpi_comm_world, -1)
                else
                    call reduce_base_real(ieskin2, cov, commcov_mpi, -1)
                end if
#endif
                call dscal(ieskin2, scalecov, cov, 1)
            end if

            if (ncg .gt. 0) then
                do i = 1, npr
                    !          write(6,*) ' Variance parameter ',i,err(i)
                    if (scalpar(i) .eq. 0.d0 .or. err(i) .eq. 0.d0) then
                        !           if(rank.eq.0) write(6,*) ' Setting to zero parameter =',i
                        call dscalzero(ncg, 0.d0, reduce(1, i), 1)
                    end if
                end do
            end if

            !          if(rank.eq.0) then

            do i = 1, nbinmax
                ipip(i) = i
            end do

            !        now replace the less relevant directions with
            !        standard parameters directions

            !        write(6,*) ' nbinmax read =',nbinmax,indmax,ind
            cost = 0.d0
            do i = 1, npk
                if (err(i) .ne. 0.d0) then
                    psip(n4 + i - 1) = abs(forza(i)/err(i))
                else
                    psip(n4 + i - 1) = 0.d0
                end if
            end do

            call dsortx(psip(n4), 1, npk, ipip(n4))

            do ii = 1, npk
                i = ipip(n4 + npk - ii)
                ! at least the zero variance is left
                if (ind + ieskin .lt. indmax) then
                    if (err(i) .ne. 0.d0) then
                        cost = psip(n4 + npk - ii)
                        !   if parcutr ge 0 take in any case the forces on ions
                        if (cost .gt. parcutpar .and.&
                                &(rpar(i) .eq. 0 .or. ncg_adr .gt. 0) .and. allowed_par(i)) then
                            ind = ind + 1
                            ipip(n3 + ind - 1) = i
                            psip(ind) = err(i)**2
                            psip(n2 + ind - 1) = psip(n3 + i - 1)
                        end if
                    end if
                end if
            end do
            !             Take in any event the forces
            do i = npk + 1, np
                ind = ind + 1
                ipip(n3 + ind - 1) = i
                psip(ind) = err(i)**2
                psip(n2 + ind - 1) = psip(n3 + i - 1)
            end do

            if (ind - indmin .eq. 0 .and. indmin .eq. 0) then
                !        take at least one parameter
                if (rank .eq. 0) write (6, *) ' Warning I take at least one parameter '
                i = ipip(n4 + npk - 1)
                if (err(i) .ne. 0.d0) then
                    cost = abs(forza(i)/err(i))
                    ind = ind + 1
                    ipip(n3 + ind - 1) = i
                    psip(ind) = err(i)**2
                    psip(n2 + ind - 1) = psip(n3 + i - 1)
                else
                    if (rank .eq. 0) write (6, *) ' There should be some error !!! '
                    errnoise = 10
                end if
            end if

            numpar = ind
            if (rank .eq. 0) write (6, *) ' Normal parameters considered ', ind - indmin
            !        nbinmax-ncg=npbra= Max  # Normal parameters  input unchanged

            !        if(ind.eq.indmax) write(6,*)
            !    1' Warning increase nbin or parcutmin !!! '
            do i = indmin + 1, ind
                if (rank .eq. 0) write (6, *) i - indmin, ipip(n3 + i - 1)
            end do

            devmaxp = 0.d0
            if (ieskin .ne. 0 .and. idyn .ge. 0) then
                do i = np - ieskin + 1, np
                    if (err(i) .ne. 0.d0 .and. allowed_par(i)) then
                        cost = abs(forza(i)/err(i))
                        if (cost .gt. devmaxp) then
                            jmax = i
                            devmaxp = cost
                        end if
                    end if
                end do
                if (rank .eq. 0) write (6, *) ' devmax par ions  =', devmaxp, jmax
            end if

            devmaxc = devmaxp
            jmax = 0

            !               write(6,*) 'PARAMETERS per KPOINT:',forza(1:np-ieskin),err(1:np-ieskin),rank,rankcolrep

            devmaxp = 0.d0
            do i = 1, np - ieskin
                if (err(i) .ne. 0.d0 .and. (rpar(i) .eq. 0.d0 .or. ncg_adr .gt. 0)&
                        &.and. allowed_par(i)) then
                    cost = abs(forza(i)/err(i))
                    if (cost .gt. devmaxp) then
                        jmax = i
                        devmaxp = cost
                    end if
                end if
            end do
#ifdef PARALLEL
            cost = devmaxp
            call mpi_allreduce(cost, devmaxp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif

#ifdef __DEBUG
            if (rank .eq. 0) write (6, *) ' i_max / force / error / deviation: ', &
                jmax, forza(jmax), err(jmax), forza(jmax)/cost
#endif
            if (rank .eq. 0) write (6, *) ' devmax par Normal   =', devmaxp, jmax, iesconv
            !             devmaxp=max(devmaxp,devmaxc)

            !         Replace the vector forza with the sorted ones in the order chosen.
            call dcopy(nbinmax, psip(n2), 1, forza, 1)
            call dcopy(nbinmax, psip, 1, err, 1)

            !        end reduction
            !        err(i)=0.d0 takes into account the "dirty" parameters
            !        NOT the zero variance ones.
            if (np .ge. ncg) then
                do i = 1, nbinmax
                    if (err(i) .gt. 0.d0) then
                        err(i) = dsqrt(abs(err(i)))
                    else
                        err(i) = 0.d0
                    end if
                end do

            else

                do i = 1, npr
                    if (err(i) .gt. 0.d0) then
                        err(i) = dsqrt(abs(err(i)))
                    else
                        err(i) = 0.d0
                    end if
                end do
                do i = npr + 1, ncg
                    err(i) = 0.d0
                end do

            end if

            do i = 1, nbinmax
                itouch(i) = 0
            end do
            itouch(nbinp) = 1

            !         replace psip with the criteria of devmax

            !        Take in any event the parameter with largest fluctuation

            do i = 1, numpar
                !         numpar=numpar+1
                itouch(i) = 1
            end do

            numparp = numpar + 1

            devmaxc = devmaxp
            iesdelay = iesdelay + 1
            if (devmaxc .le. dble(parcutpar + 1.5) .and. devmaxc .ne. 0.d0) then
                if (iesconv .ge. 0) then
                    iesconv = iesconv + 1
                else
                    iesconv = 0
                    if (rank .eq. 0) write (6, *) ' Warning reinitializing iesconv'
                end if
                if (iesconv .ge. maxiter_changeparr .and. maxiter_changeparr .gt. 0&
                        &.and. change_parr .and. cut .lt. parr_min*1.0000001d0) then
                    if (rank .eq. 0) write (6, *) ' Wanderful Warning:  Optimization converged !!!'
                    ngentry = -5
                end if
                if (iesconv .ge. 3 .and. change_parr .and. iesdelay .ge. delay_changeparr) then
                    iesdelay = 0
                    cost = abs(parr)/10.d0
                    if (cost .ge. parr_min*0.9999999d0) then
                        parr = cost
                        if (default_epsdgel) then
                            epst = parr/10.d0
                        end if
                        if (yes_umrigar) parr = -parr
                    end if
                end if
                if (rank .eq. 0 .and. abs(parr) .lt. cut) then
                    write (6, *) ' New decreased parr =', abs(parr)
                    iesconv = 0
                end if
            elseif (devmaxc .ne. 0) then
                if (iesconv .gt. 0) then
                    iesconv = 0
                    if (rank .eq. 0) write (6, *) ' Warning reinitializing iesconv'
                elseif (devmaxc .gt. dble(parcutpar + 3)) then
                    iesconv = iesconv - 1
                end if
                if (iesconv .le. -15 .and. change_parr .and. iesdelay .ge. delay_changeparr) then
                    iesdelay = 0
                    cost = abs(parr)*10.d0
                    if (cost .le. parr_max*1.0000001d0) then
                        parr = cost
                        if (default_epsdgel) then
                            epst = parr/10.d0
                        end if
                        if (yes_umrigar) parr = -parr
                    end if
                end if
                if (rank .eq. 0 .and. abs(parr) .gt. cut) then
                    write (6, *) ' New increased parr =', abs(parr)
                    iesconv = 0
                end if
            end if

            do i = 1, numpar
                itouch(i) = 1
                ipip(i) = 1
            end do
            !      kill the component which are set to zero by scalpar

            !  fine if rank.eq.0
            !          endif

            !   if(rank.eq.0) write(6,*) ' force before kl ',np,sum(abs(psip(n3:n3+np-1)))
            k_par = 0

            if (abs(klr) .eq. 6 .or. abs(klr) .eq. 2) then
                if (rank .eq. 0) time = cclock()

                !   ridotto  il numero di parametri se rpar=/0
                call update_index

                if (ncg .gt. 0) then

#ifdef    PARALLEL

!         nprocsr=nw/in1

                    if (abs(klr) .ne. 2) then
                        if (k6gen) then
                        if (abs(parr) .ge. 1d-3) then
                            maxdim_inv = 1000
                        else
                            maxdim_inv = 2000
                        end if
                        else
                        if (abs(parr) .ge. 1d-3) then
                            maxdim_inv = 5000
                        else
                            maxdim_inv = 10000
                        end if
                        end if

                        nprocm = nprocsr - 1
                        maxit = 20*ku_ionp
                        if (nprocm .ne. 0) then
                            block = max(ku_ionp/nint(sqrt(dble(nprocm + 1))), 1) ! Block has to be at least 1
                            if (block .lt. min_block .and. ku_ionp .gt. maxdim_inv) block = min_block
                            nprocm = ku_ionp/block
                        else
                            block = ku_ionp
                        end if
                        if (block*(nprocm + 1) .lt. ku_ionp .or. ku_ionp .le. maxdim_inv) then
                            if (ku_ionp .le. maxdim_inv) then
                                block = ku_ionp
                            else
                                block = block + 1
                            end if
!           In this case some processor do not store the matrix
                            nprocm = ku_ionp/block
                        end if
                        if (ipc .eq. 2 .and. srcomplex .and. mod(block, 2) .ne. 0 .and. block .ne. ku_ionp) then
                            if (rank .eq. 0) write (6, *) ' Warning block even in this case !!! '
                            block = block + 1 !even blocking
                            nprocm = ku_ionp/block
                        end if
                        if (block*nprocm .eq. ku_ionp) then
!           In this case it is convenient not to store the matrix
!           in the last processor
                            nprocm = nprocm - 1
                        end if
                        if (k6gen .and. nprocm .gt. nproc_ortho) then
                            if (rank .eq. 0) write (6, *) ' Warning # processors changed for k6gen on pool! '
                            block = max(ku_ionp/nproc_ortho, 1)
                            if (block*nproc_ortho .ne. ku_ionp) block = block + 1
                            if (ipc .eq. 2 .and. mod(block, 2) .eq. 1) block = block + 1 ! in the complex case always even
                            nprocm = ku_ionp/block - 1
                            if ((nprocm + 1)*block .lt. ku_ionp) nprocm = nprocm + 1
                        end if

!         The following choice is for maximize speed up
!         blockread << the choice below to optimize memory.
!         However at most a factor two of memory is saved.

                        ! Has to be multiple of in1
                        blockread = in1*(max(block, in1)/in1)

                        ! make a block large enough
                        if (blockread .le. 1024 .and. in1 .lt. 1024) then
                            blockread = min((1024/in1)*in1, in1*(nweight - ibinit))
                        end if

                        npsr = ku_ionp*block
                        blockdim = block*block

!         EACH processor has its own block Columns

                        allocate (sr(npsr), srg(blockdim), srb(blockdim)                 &
                   &, eb(ku_ionp, blockread))

#ifdef UNREL
                        allocate (work(blockdim))
#endif

                        allocate (v1(ku_ion), v2(ku_ion), v3(ku_ion)                     &
                   &, v4(2*block*(nprocsr)))

#else
                        ! has to be a multiple of in1
                        blockread = in1*(max(ku_ionp, in1)/in1)

                        npsr = ku_ionp*ku_ionp
                        allocate (sr(npsr), eb(ku_ionp, blockread))
#endif

                        call dscalzero(npsr, 0.d0, sr, 1)

                        if (writescratch .ne. 0) then
                            countscra = ibinit
                        else
                            do kk = 1, ibinit
                                read (indopen3)
                            end do
                        end if
                        allocate (factorsr_sav(in1))

                        irep = blockread/in1
                        iterk = (nweight - ibinit)/irep
                        irest = nweight - ibinit - iterk*irep
                        if (irest .ne. 0) iterk = iterk + 1

                        do kk = 1, iterk
                            if (kk .eq. iterk .and. irest .ne. 0) then
                                irepr = irest
                            else
                                irepr = irep
                            end if

                            do jj = 1, irepr
                                if (writescratch .ne. 0) then
                                    countscra = countscra + 1
                                    indscra = (countscra - 1)*nbufscra + 1
                                    cost = bufscra(indscra)
                                    do i = 1, in1
                                        do j = 1, np
                                            econf(i, j) = bufscra(indscra + 1)
                                            econfh(i, j) = bufscra(indscra + 2)
                                            indscra = indscra + 2
                                        end do
!                               factorsr(i) = bufscra(indscra + 1)
                                        factorsr_sav(i) = bufscra(indscra + 1)
                                        if (ipc .eq. 2) then
                                            enert(1, i) = bufscra(indscra + 2)
                                            enert(2, i) = bufscra(indscra + 3)
                                            indscra = indscra + 3
                                        else
                                            enert(1, i) = bufscra(indscra + 2)
                                            indscra = indscra + 2
                                        end if
                                    end do
                                else
                                    read (indopen3) cost, ((econf(i, j), econfh(i, j), j=1, np)            &
                                            &, factorsr_sav(i), enert(1:ipc, i), i=1, in1)
                                end if

                                !          Output a given Processor i all block columns
                                do i = 1, in1
                                    costw = dsqrt(cost*factorsr_sav(i)/weightall)
                                    ind = in1*(jj - 1) + i
                                    do ii = 1, ku_ion
                                        eb(ii, ind) = (econf(i, index(ii)) - okav(index(ii)))*costw
                                    end do
                                    eb(ku_ionp, ind) = costw
                                    !         call dscal(npp,costw,econf(i,1),in1)
                                end do
                                !enddo irepr
                            end do
                            dimbuf = irepr*in1

                            !         STANDARD (SCALAPACK?)  MATRIX_MATRIX
                            !         SR = SR+EB (EB)^T
                            !         BLOCK MATRIX DIM = block
                            !         OUTPUT SR   #Column=block  for each Processor
                            deallocate (factorsr_sav)
#ifdef PARALLEL

!      Initialize block matrices
                            call dscalzero(blockdim, 0.d0, srg, 1)
                            call dscalzero(blockdim, 0.d0, srb, 1)

                            do j = 0, nprocm

                                jblock = j*block + 1

                                blockl = block
                                if (j .eq. nprocm) blockl = ku_ionp - jblock + 1

                                do i = 0, nprocm
                                    iblock = i*block + 1
                                    blockr = block
                                    if (i .eq. nprocm) blockr = ku_ionp - iblock + 1

! IN THE LOOP BELOW ALL PROCESSORS WORK FOR THE i-th PROCESSOR sr-matrix

                                    call dgemm('N', 'T', blockl, blockr, dimbuf, 1.d0                    &
                                 &, eb(jblock, 1), ku_ionp, eb(iblock, 1), ku_ionp, 0.d0, srb, block)

!  I asume that srg is updated only for the ROOT processor=i
#ifdef UNREL
                                    call reduce_base_real_to(blockdim, srb, work, commsr_mpi, -1)
                                    if (ranksr .eq. i) srg = work
#else
                                    call reduce_base_real_to(blockdim, srb, srg, commsr_mpi, i)
#endif
                                end do

!       Here we sum alltoghether
!       block_r=block
!       if(rank.eq.nprocm) block_r=npp-nprocm*block
                                do ii = 1, block
                                    do jj = jblock, jblock + blockl - 1
                                        sr(ku_ionp*(ii - 1) + jj) = sr(ku_ionp*(ii - 1) + jj)  &
                                    &   + srg(block*(ii - 1) + jj - jblock + 1)
                                    end do
                                end do
                            end do
!     END STANDARD PARALLEL (SCALAPACK?)  MATRIX-MATRIX

#else

                            !   STANDARD SCALAR MATRIX-NATRIX
                            call dgemm('N', 'T', ku_ionp, ku_ionp, dimbuf, 1.d0         &
                                    &, eb, ku_ionp, eb, ku_ionp, 1.d0, sr, ku_ionp)

                            !            do i=ist,ien
                            !           costw=cost*wconfn(i)
                            !           call dsyr('U',npp,costw,econf(i-istm,1),in1,sr,npp)
                            !           enddo
#endif
                            ! enddo iterk
                        end do

#ifdef PARALLEL
                        if (ranksr .eq. nprocm) then
                            ind = (ku_ionp - block*nprocm - 1)*ku_ionp + 1
                            call dcopy(ku_ionp, sr(ind), 1, psip, 1)
                        end if

                        cost = 1.d0/psip(ku_ionp)
!        Only nprocm has computed psip as above
                        call bcast_real(cost, 1, nprocm, commsr_mpi)

                        call dscal(npsr, cost, sr, 1)
                    else
                        nprocm = 0 ! go to the non parallel part as the matrix is not necessary
                    end if ! endif abs(klr)==2

                    skdiag(1:ku_ion) = 1.d0
                    do i = 1, ku_ion
                        ii = index(i)
                        if (allowed_par(ii)) then
                            psip(i) = forza_sav(ii)
                            if (scalpar(ii) .eq. 0.d0) skdiag(i) = 0.d0
                        else
                            psip(i) = 0.d0
                            skdiag(i) = 0.d0
                        end if
                    end do
                    if (rank .eq. 0) then
                        write (6, *) ' loading time  ', cclock() - time
                        time = cclock()
                    end if
!        Save known term.
                    psip(n4:n4 + ku_ion - 1) = psip(1:ku_ion)
                    if (nprocm .gt. 0) then
                        if (srcomplex) then
                            do i = 1, ku_ion
                                do jj = 1, block
                                    j = ranksr*block + jj
                                    if (j .le. ku_ion) then
!          do i=1,ku_ion
!           do j=1,ku_ion
                                        if (truecomplex(i) .and. truecomplex(j)) then
                                            if (symmagp) then
                                                cost = (sr(ku_ionp*(jj - 1) + i) + sr(ku_ionp*jj + i + 1)) ! s^rr+s^ii
                                                sr(ku_ionp*(jj - 1) + i) = cost
                                                sr(ku_ionp*jj + i + 1) = 0.d0
                                                sr(ku_ionp*jj + i) = 0.d0
                                                sr(ku_ionp*(jj - 1) + i + 1) = 0.d0
                                            else
                                                cost = (sr(ku_ionp*(jj - 1) + i) + sr(ku_ionp*jj + i + 1)) ! s^rr+s^ii
                                                sr(ku_ionp*(jj - 1) + i) = cost
                                                sr(ku_ionp*jj + i + 1) = cost
                                                cost = (sr(ku_ionp*jj + i) - sr(ku_ionp*(jj - 1) + i + 1)) ! s^ri-s^ir
                                                sr(ku_ionp*jj + i) = cost
                                                sr(ku_ionp*(jj - 1) + i + 1) = -cost
                                            end if
                                        elseif (i .gt. kp_complex .and. truecomplex(j) .and. symmagp) then
                                            sr(ku_ionp*jj + i) = 0.d0
                                        elseif (j .gt. kp_complex .and. truecomplex(i) .and. symmagp) then
                                            sr(ku_ionp*(jj - 1) + i + 1) = 0.d0
                                        end if
                                    end if
                                end do
                            end do
                        end if
                        cost = cut
                        if (k6gen) then
                            nbinr = ku_ion - ranksr*block
                            if (nbinr .gt. block) nbinr = block
                            if (nbinr .gt. 0) then
                                if (ranksr .ne. 0) skdiag(1:ranksr*block) = 0.d0
                                if (rank*block + nbinr .lt. ku_ion) skdiag(ranksr*block + nbinr + 1:ku_ion) = 0.d0
                                do j = 1, nbinr
                                    indj = ranksr*block + j
                                    if (skdiag(indj) .ne. 0.d0) skdiag(indj) = sr(ku_ionp*(j - 1) + indj)
                                end do
                            else
                                skdiag(1:ku_ion) = 0.d0
                            end if
#ifdef PARALLEL
                            call reduce_base_real(ku_ion, skdiag, commsr_mpi, -1)
#endif
                            maxsr = 0.d0
                            do j = 1, ku_ion
                                if (skdiag(j) .gt. maxsr) maxsr = skdiag(j)
                            end do
                            do j = 1, ku_ion
                                if (skdiag(j) .gt. maxsr*parcut2 .or. yes_umrigar) then
                                    skdiag(j) = dsqrt(skdiag(j) + eps_umrigar)
                                else
                                    if(rank.eq.0.and.allowed_par(index(j))&
                                        & .and..not.yes_umrigar) write(6,*)&
                                        & ' Warning sr parameter too small !!!'&
                                        & ,index(j),skdiag(j)
                                    if (yes_umrigar .and. allowed_par(index(j))) then
                                        skdiag(j) = dsqrt(maxsr*parcut2)
                                    else
                                        skdiag(j) = 0.d0
                                    end if
                                end if
                            end do
                            ind = 0

!        write(6,*) ' INFO =',rank,ku_ion,block,nproc_ortho

                            do j = 1, ku_ion
                                if (skdiag(j) .ne. 0.d0) then
!         ind=ind+1
!          IF(rank.eq.0) write(6,*) ' Known k6gen =',ind,j,ku_ion,psip(j),skdiag(j)
                                    psip(j) = psip(j)/skdiag(j)
                                    do k = 1, nbinr
                                        indk = ranksr*block + k
                                        if (skdiag(indk) .ne. 0.d0) then
                                            if (indk .eq. j .and. yes_umrigar) &
                                           &sr(ku_ionp*(k - 1) + j) = sr(ku_ionp*(k - 1) + j)*(1.d0 + cut) + cut_umrigar
                                            sr(ku_ionp*(k - 1) + j) = sr(ku_ionp*(k - 1) + j)/skdiag(j)/skdiag(indk)
                                            if (indk .eq. j .and. .not. yes_umrigar) &
                                           &sr(ku_ionp*(k - 1) + j) = sr(ku_ionp*(k - 1) + j) + cut
                                        else
                                            if (indk .eq. j) then
                                                sr(ku_ionp*(k - 1) + j) = 1.d0
                                            else
                                                sr(ku_ionp*(k - 1) + j) = 0.d0
                                            end if
                                        end if
                                    end do
                                else
                                    psip(j) = 0.d0
                                    do k = 1, nbinr
                                        indk = ranksr*block + k
                                        if (indk .eq. j) then
                                            sr(ku_ionp*(k - 1) + j) = 1.d0
                                        else
                                            sr(ku_ionp*(k - 1) + j) = 0.d0
                                        end if
                                    end do
                                end if
                            end do
!         if(nbinr.gt.0) then
!         cost=sum(sr(:))
!         else
!         cost=0.d0
!         endif
!     call reduce_base_real(1,cost,commsr_mpi,-1)
!         if(rank.eq.0) write(6,*) ' Sum matrix elements  before input k6gen ',cost,sum(psip(1:ku_ion))
!       write(6,*) ' Solution before ',rank,sum(psip(1:ku_ion)),sum(psip(1:601)),sum(psip(602:1202))

                            if (rank .eq. 0) write (6, *) ' Parallel linear system Blacs k6gen=true '

!      call mpi_finalize(ierr)
!      stop

!      write(6,*) ' Final comm =',ranksr,ortho_comm

                            call para_syst(ku_ion, ku_ionp, block, ranksr, nproc_ortho, ortho_comm, commsr_mpi, ortho_cntx&
                           &, psip, sr)
!        if(rank.eq.0) write(6,*) ' Output vec '
                            do i = 1, ku_ion
                            if (skdiag(i) .ne. 0.d0) then
                                v3(i) = psip(i)/skdiag(i)
                            else
                                v3(i) = 0.d0
                            end if
!        if(rank.eq.0) write(6,*) i,skdiag(i),v3(i)
                            end do

!        call mpi_finalize(ierr)
!        stop

                        else
                            call conjginvs(ku_ion, ku_ionp, block, rank, ranksr, commsr_mpi, sr, psip&
                        &, v1, v2, v4, maxit, tolcg, cost, skdiag, v3, parcut2, nprocsr, eps_umrigar)
                        end if
                        normsr = sum(psip(n4:n4 + ku_ion - 1)*v3(1:ku_ion))

                    else
!#ifdef PARALLEL
!         call bcast_real(sr,size(sr),0,row_comm)
!#endif
!       if(rank.eq.0) write(6,*) ' Warning scalar inversion '
!       add regularization

! The maximum diagonal of the SR matrix set the scale of the accuracy
!       if(rank.eq.0) then
!       write(6,*) ' Known term ',ku_ionp
!       do i=1,ku_ion
!       write(6,*) i,index(i),psip(i),okav(index(i)),sr(ku_ionp*(i-1)+i)
!       enddo
!       endif
!    D psi = sum_k  O_k dalpha_k
!  SR= < Ok^* O_k'>
!  sum_k'  SR_kk' dalpha_k' = fk
!  In the code Ok^*= is stored =X^r+i X^i
!  SR = s^rr + s^i,i -i (s^ri-s^ir)
!  dalpha = dalpha^r +i dalpha^i
!  is equivalent to a 2 times larger real system of equation:
!    |  s^rr+s^ii, s^ri-s^ir | | dalpha^r |= | f^r|
!    |  s^ir-s^ri, s^rr+s^ii | | dalpha^i |= | f^i|
!  When a parameter (O_k) is real it can be shown nothing change.

                        if (abs(klr) .ne. 2) then
                            if (srcomplex) then
                                do i = 1, ku_ion
                                    do j = 1, ku_ion
                                        if (truecomplex(i) .and. truecomplex(j)) then
                                            if (symmagp) then
                                                cost = sr(ku_ionp*(j - 1) + i) + sr(ku_ionp*j + i + 1)
                                                sr(ku_ionp*(j - 1) + i) = cost
                                                sr(ku_ionp*j + i + 1) = 0.d0
                                                sr(ku_ionp*j + i) = 0.d0
                                                sr(ku_ionp*(j - 1) + i + 1) = 0.d0
                                            else
                                                cost = sr(ku_ionp*(j - 1) + i) + sr(ku_ionp*j + i + 1)
                                                sr(ku_ionp*(j - 1) + i) = cost
                                                sr(ku_ionp*j + i + 1) = cost
                                                cost = sr(ku_ionp*j + i) - sr(ku_ionp*(j - 1) + i + 1)
                                                sr(ku_ionp*j + i) = cost
                                                sr(ku_ionp*(j - 1) + i + 1) = -cost
                                            end if
                                        elseif (i .gt. kp_complex .and. truecomplex(j) .and. symmagp) then
                                            sr(ku_ionp*j + i) = 0.d0 ! cancel imaginary part
                                        elseif (j .gt. kp_complex .and. truecomplex(i) .and. symmagp) then
                                            sr(ku_ionp*(j - 1) + i + 1) = 0.d0 ! cancel imaginary part
                                        end if
                                    end do
                                end do
                            end if

#ifdef __DEBUG
                            if (rank .eq. 0) then
!       check symmetry
                                cost = 0.d0
                                do i = 1, ku_ion
                                    do j = i + 1, ku_ion
                                        cost = cost + abs(sr(ku_ionp*(j - 1) + i) - sr(ku_ionp*(i - 1) + j))
                                    end do
                                end do

                                write (6, *) ' Check symmetry =', cost
                                call dsyev('V', 'L', ku_ion, sr, ku_ionp, psip, psip(n2), 3*ku_ion, info)
                                write (6, *) ' diagonal part '
                                do i = 1, ku_ion
                                if (psip(i) .le. 0.) then
                                    write (6, *) ' ERROR negative eigenvalue '
                                end if
                                write (6, *) i, psip(i)
                                end do
                            end if
                            if (rank .eq. 0) write (6, *) ' Everything works ! '
#ifdef PARALLEL
                            call mpi_finalize(ierr)
#endif
                            stop
#endif
                            cost = cut
                            maxsr = 0.d0
                            do i = 1, ku_ion
                                if (sr(ku_ionp*(i - 1) + i) .gt. maxsr) maxsr = sr(ku_ionp*(i - 1) + i)
                            end do
                            do i = 1, ku_ion
                                ii = index(i)
                                if (sr(ku_ionp*(i - 1) + i) .gt. parcut2*maxsr .and. allowed_par(ii)&
                             &.or. (yes_umrigar .and. allowed_par(ii))) then
                                    if (yes_umrigar) then
                                        sr(ku_ionp*(i - 1) + i) = sr(ku_ionp*(i - 1) + i)*(1.d0 + cost) + eps_umrigar
                                    else
                                        sr(ku_ionp*(i - 1) + i) = sr(ku_ionp*(i - 1) + i)*(1.d0 + cost)
                                    end if
                                else
                                    if (rank .eq. 0 .and. allowed_par(ii)) write (6, *) ' Warning sr parameter too small !!!'&
                                  &, index(i), forza_sav(index(i)), sr(ku_ionp*(i - 1) + i)
                                    psip(i) = 0.d0
                                    do j = 1, ku_ion
                                        sr(ku_ionp*(i - 1) + j) = 0.d0
                                        sr(ku_ionp*(j - 1) + i) = 0.d0
                                    end do
                                    sr(ku_ionp*(i - 1) + i) = 1.d0
                                end if
                            end do
!       now calculation inverse
                        end if ! abs(klr).eq.2
                        if (.not. allocated(v3)) allocate (v3(ku_ion))
                        if (abs(klr) .eq. 2) then
                        do ii = 1, ku_ion
                            if (allowed_par(index(ii))) then
                                v3(ii) = psip(ii)
                            else
                                v3(ii) = 0.d0
                            end if
                        end do
                        else
                        v3(1:ku_ion) = psip(1:ku_ion)
                        end if

!        if(rank.eq.0) then
!         do i=1,npsr
!         if(sr(i).eq.1.d0) write(6,*) ' sr =1 ', i-((i-1)/ku_ionp)*ku_ionp,(i-1)/ku_ionp+1,sr(i)
!         enddo
!        write(6,*) ' Sum matrix elements right =',sum(sr(:)),sum(v3(1:ku_ion))
!        do i=1,ku_ion
!        write(6,*) ' Known right =',i,v3(i)
!        enddo
!        endif
!       Choleski decomposition
                        if (abs(klr) .eq. 2) then
                            if (rank .eq. 0) write (6, *) ' Warning steepest descend '
                        else
                            call dpotrf('L', ku_ion, sr, ku_ionp, info)
                            call dpotrs('L', ku_ion, 1, sr, ku_ionp, psip, ku_ion, info)
                        end if
!       LU slower and less stable
!        call dgetrf(ku_ion,ku_ion,sr,ku_ionp,ipip(n4),info)
!        call dgetrs('N',ku_ion,1,sr,ku_ionp,ipip(n4),psip,ku_ion,info)
                        normsr = sum(v3(1:ku_ion)*psip(1:ku_ion))
                        v3(1:ku_ion) = psip(1:ku_ion)
!         if(rank.eq.0) then
!         do i=1,ku_ion
!         write(6,*) ' Sol right =',i,v3(i)
!         enddo
!         endif
!         call mpi_finalize(ierr)
!         stop

                    end if

                    if (rank .eq. 0) then
                        time = cclock() - time
                        write (6, *) ' Time inverse Conjugate  Gradients =', time
                    end if
!         if(rank.eq.0) write(6,*) ' Output reduce '
                    if (sum(abs(reduce(1, 1:npr))) .eq. 0.d0) then
                        beta_used = 0.d0
                    else
                        beta_used = beta_learning
                        do i = 1, npr
                            reduce(1, i) = beta_used*reduce(1, i)
                        end do
                    end if
                    do i = 1, ku_ion
                        ii = index(i)
                        if (allowed_par(ii)) reduce(1, ii) = reduce(1, ii) + (1.d0 - beta_used)*v3(i)
!         if(rank.eq.0)     write(6,*) i,index(i),reduce(1,index(i))
                    end do

!         call mpi_finalize(ierr)
!         stop

                    if (allocated(sr)) deallocate (sr, srg, srb, eb, v1, v2, v4)
                    deallocate (v3)

#ifdef UNREL
                    deallocate (work)
#endif

#else
                    write (6, *) ' loading time  ', cclock() - time
                    !        normalization matrix sr
                    !        write(6,*) ' Normalization const =',sr(npsr)
                    cost = 1.d0/sr(npsr)
                    call dscal(npsr, cost, sr, 1)
                    ! the hessian is known only for the ion
                    ndsr = ku_ion
                    ! default value
                    if (srcomplex) then
                        !          if(rank.eq.0) then
                        !          write(6,*) ' Input truecomplex '
                        !          do i=1,ku_ion
                        !          write(6,*) i,truecomplex(i)
                        !          enddo
                        !          endif
                        do i = 1, ndsr
                            do j = 1, ndsr
                                if (truecomplex(i) .and. truecomplex(j)) then
                                    if (symmagp) then
                                        cost = sr(ku_ionp*(j - 1) + i) + sr(ku_ionp*j + i + 1)
                                        sr(ku_ionp*(j - 1) + i) = cost
                                        sr(ku_ionp*j + i + 1) = 0.d0
                                        sr(ku_ionp*j + i) = 0.d0
                                        sr(ku_ionp*(j - 1) + i + 1) = 0.d0
                                    else
                                        cost = sr(ku_ionp*(j - 1) + i) + sr(ku_ionp*j + i + 1)
                                        sr(ku_ionp*(j - 1) + i) = cost
                                        sr(ku_ionp*j + i + 1) = cost
                                        cost = sr(ku_ionp*j + i) - sr(ku_ionp*(j - 1) + i + 1)
                                        sr(ku_ionp*j + i) = cost
                                        sr(ku_ionp*(j - 1) + i + 1) = -cost
                                    end if
                                elseif (i .gt. kp_complex .and. truecomplex(j) .and. symmagp) then
                                    sr(ku_ionp*j + i) = 0.d0
                                elseif (j .gt. kp_complex .and. truecomplex(i) .and. symmagp) then
                                    sr(ku_ionp*(j - 1) + i + 1) = 0.d0
                                end if
                            end do
                        end do
                    end if
                    ! The maximum diagonal of the SR matrix set the scale of the accuracy
                    maxsr = 0.d0
                    do i = 1, ku_ion
                        if (sr(ku_ionp*(i - 1) + i) .gt. maxsr) maxsr = sr(ku_ionp*(i - 1) + i)
                    end do
                    do i = 1, ndsr
                        ii = index(i)
                        fkav(ii) = sr((i - 1)*ku_ionp + i)
                        if (yes_umrigar .and. fkav(ii) .le. parcut2*maxsr) fkav(ii) = parcut2*maxsr
                        if ((fkav(ii) .le. parcut2*maxsr .and. .not. yes_umrigar)&
                                &.or. scalpar(ii) .eq. 0.d0 .or. .not. allowed_par(ii)) then
                            fkav(ii) = 0.d0
                            do j = 1, ndsr
                                if (j .ne. i) then
                                    ind = (j - 1)*ku_ionp + i
                                    indr = (i - 1)*ku_ionp + j
                                    sr(ind) = 0.d0
                                    sr(indr) = 0.d0
                                else
                                    ind = (i - 1)*ku_ionp + i
                                    sr(ind) = 1.d0
                                end if
                            end do
                        end if
                    end do
                    cost = cut
                    do i = 1, ndsr
                        ind = (i - 1)*ku_ionp + i
                        if (yes_umrigar) then
                            sr(ind) = sr(ind)*(1.d0 + cost) + eps_umrigar
                        else
                            sr(ind) = sr(ind)*(1.d0 + cost)
                        end if
                    end do

                    time = cclock()

                    !        now calculation inverse
                    !        Choleski decomposition
                    call dpotrf('L', ndsr, sr, ku_ionp, info)
                    !        LU slower and less stable
                    !        call dgetrf(ndsr,ndsr,sr,ku_ionp,ipip(n4),info)

                    do i = 1, ku_ion
                        ii = index(i)
                        if (allowed_par(ii)) then
                            psip(i) = forza_sav(ii)
                        else
                            psip(i) = 0.d0
                        end if
                    end do

                    if (ndsr .gt. 0) call dpotrs('L', ndsr, 1, sr, ku_ionp, psip, ndsr, info)
                    !     if(ndsr.gt.0) call dgetrs('N',ndsr,1,sr,ku_ionp,ipip(n4),psip,ndsr,info)

                    time = cclock() - time

                    write (6, *) ' Inverse time kl=6 ', time
                    if (sum(abs(reduce(1, 1:npr))) .eq. 0.d0) then
                        beta_used = 0.d0
                    else
                        beta_used = beta_learning
                        do i = 1, npr
                            reduce(1, i) = beta_used*reduce(1, i)
                        end do
                    end if

                    !         write(6,*) ' Output reduce ',ndsr,npr
                    do i = 1, ku_ion
                        ii = index(i)
                        if (allowed_par(ii)) then
                            reduce(1, ii) = reduce(1, ii) + (1.d0 - beta_used)*psip(i)
                            !         write(6,*) i,ii,psip(i)
                        end if
                    end do
                    do i = 1, npr
                        if (fkav(i) .eq. 0.d0) reduce(1, i) = 0.d0 ! avoid the irrelevant par.
                    end do

                    deallocate (sr, eb)

#endif
                    if (writescratch .eq. 0) rewind (indopen3)
                    !        if(allocated(index)) deallocate(index)
                end if ! if ncg>0

            end if

            if (rank .eq. 0) time = cclock()

            if (abs(klr) .eq. 7) then
                in1 = ien - ist + 1
                dimmat = in1*(nweight - ibinit)
                irep = nw/in1
                call update_index

                if (ncg .gt. 0) then

                    if (ncg_adr .gt. 0) then
                        !         definition long range collective variables
                        allocate (okav_coll(ncg_adr), econf_coll(ncg_adr, in1))
                        kn_ion = ku_ion
                        ku_ion = ku_ion + ncg_adr
                    else
                        kn_ion = ku_ion
                        allocate (okav_coll(1), econf_coll(1, 1))
                    end if

                    kn_ionp = kn_ion + 1

                    allocate (mat(dimmat, ku_ion))

                    if (prep .gt. 0) then
                        nprepm = ku_ion/prep
                        if (nprepm*prep .ne. ku_ion) nprepm = nprepm + 1
                        if (mod(nprepm, 2) .ne. 0 .and. ipc .eq. 2 .and. srcomplex) then
                            nprepm = nprepm + 1
                        end if
                    else
                        nprepm = 0
                    end if

                    mat = 0.d0

                    maxit = min(ku_ion, dimmat*irep)*20
                    allocate (v1(max(kp_ion, ku_ion)), v2(ku_ion), v3(ku_ion)         &
                            &, v4(max(ipc*ku_ion, nprepm*(2*prep + 1)) + (dimmat + ipc - 1)*max(1, prep)))

                    if (ncg_adr .gt. 0 .and. .not. signalnoise)&
                            &call dgemv('N', ncg_adr, kp0, 1.d0, reducel, ncg_adr, okav, 1, 0.d0, okav_coll, 1)

                    if (.not. signalnoise) then
                        allocate (factorsr_sav(in1))
                        if (writescratch .eq. 0) then
                            do kk = 1, ibinit
                                read (indopen3)
                            end do
                        end if

                        do kk = 1, nweight - ibinit
                            if (writescratch .ne. 0) then
                                countscra = kk + ibinit
                                indscra = (countscra - 1)*nbufscra + 1
                                cost = bufscra(indscra)
                                do i = 1, in1
                                    do j = 1, np
                                        econf(i, j) = bufscra(indscra + 1)
                                        econfh(i, j) = bufscra(indscra + 2)
                                        indscra = indscra + 2
                                    end do
                                    factorsr_sav(i) = bufscra(indscra + 1)
                                    if (ipc .eq. 2) then
                                        enert(1, i) = bufscra(indscra + 2)
                                        enert(2, i) = bufscra(indscra + 3)
                                        indscra = indscra + 3
                                    else
                                        enert(1, i) = bufscra(indscra + 2)
                                        indscra = indscra + 2
                                    end if
                                end do
                            else
                                read (indopen3) cost, ((econf(i, j), econfh(i, j), j=1, np)            &
                                        &, factorsr_sav(i), enert(1:ipc, i), i=1, in1)
                            end if
                            if (ncg_adr .gt. 0) call dgemm('N', 'T', ncg_adr, in1, kp0, 1.d0     &
                                    &, reducel, ncg_adr, econf, in1, 0.d0, econf_coll, ncg_adr)
                            do i = 1, in1
                                !           do i=ist,ien
                                costw = dsqrt(cost*factorsr_sav(i)/weightall)
                                jj = in1*(kk - 1) + i
                                do ii = 1, kn_ion
                                    kkk = index(ii)
                                    mat(jj, ii) = costw*(econf(i, kkk) - okav(kkk))
                                end do
                                if (kn_ionp .le. ku_ion) then
                                    do ii = kn_ionp, ku_ion
                                        mat(jj, ii) = costw*(econf_coll(ii - kn_ion, i) - okav_coll(ii - kn_ion))
                                    end do
                                end if
                            end do
                        end do
                        deallocate (factorsr_sav)

                    else ! signalnoise main if, below if signalnoise=.true.
                        cost = dsqrt(dble(nbin)*nprocsr)
                        do jj = 1, nbin
                            do ii = 1, ku_ion
                                kkk = index(ii)
                                mat(jj, ii) = fk(jj, kkk)*cost
                            end do
                        end do
                    end if

                    do ii = 1, kn_ion
                        kk = index(ii)
                        if (allowed_par(kk)) then
                            psip(ii) = forza_sav(kk) ! input forza
                        else
                            psip(ii) = 0.d0
                        end if
                    end do

                    if (signalnoise .and. ieskin .ne. 0 .and. idyn .eq. 0) then
                        if (tion .ne. tjas) then
                            costf = dsqrt(tion/tjas)
                            do ii = kn_ion - ieskin + 1, kn_ion - ieskin + ieskinr_pos
                                psip(ii) = psip(ii)*costf
                            end do
                        end if
                        if (tcell .ne. tjas) then
                            costc = dsqrt(tcell/tjas)
                            do ii = kn_ion - ieskin + ieskinr_pos + 1, kn_ion
                                psip(ii) = psip(ii)*costc
                            end do
                        end if
                    end if

                    if (kn_ionp .le. ku_ion .and. .not. signalnoise) then
                        do ii = kn_ionp, ku_ion
                            force_coll = 0.d0
                            do jj = 1, kp0
                                force_coll = force_coll + forza_sav(jj)*reducel(ii - kn_ion, jj)
                            end do
                            psip(ii) = force_coll
                        end do
                    end if
                    !          Now computation  skdiag
                    skdiag(1:ku_ion) = 1.d0
                    do ii = 1, kn_ion
                        kk = index(ii)
                        if (.not. allowed_par(kk)) skdiag(ii) = 0.d0
                    end do
                    cost = cut
                    psip(n4:n4 + ku_ion - 1) = psip(1:ku_ion)
                    if (rank .eq. 0) then
                        write (6, *) ' loading time  ', cclock() - time
                        time = cclock()
                    end if
                    yes_targetprep = .false.
#ifdef   PARALLEL
                    if (prep .gt. 0) then
                        allocate (mat_prep(nprepm, dimmat*prep), mat_buf(nprepm, dimmat))
#ifdef _OFFLOAD
                        if (nprepm*dimmat*prep .gt. max_targetsr) yes_targetprep = .true.
                        if (rank .eq. 0 .and. yes_targetprep) write (6, *) 'Warning SR inversion with GPU'
#endif
                        mat_prep = 0.d0
                        mat_buf = 0.d0
!         now gathering the raw
                        timep = cclock()

                        do i = 0, prep - 1
                            mini = i*nprepm + 1
                            maxi = min(mini + nprepm - 1, ku_ion)
                            mat_buf = 0.d0
                            do j = mini, maxi
                                mat_buf(j - mini + 1, 1:dimmat) = mat(1:dimmat, j)
                            end do

                            call mpi_gather(mat_buf, nprepm*dimmat, MPI_DOUBLE_PRECISION, mat_prep&
                           &, nprepm*dimmat, MPI_DOUBLE_PRECISION, i, comm_raw, ierr)
                        end do

                        if (rank .eq. 0) write (6, *) ' Time reshuff matrix =', cclock() - timep

                        timep = cclock()
                        call conjginv_prep(ku_ion, prep, nprepm, kp_complex, symmagp, dimmat*prep, rank&
                         &, ranksr, commsr_mpi, comm_raw, comm_col, mat_prep, psip, v1, v2, v4, maxit, tolcg&
                         &, cost, skdiag, v3, parcut, eps_umrigar, yes_targetprep)

                        if (rank .eq. 0) write (6, *) ' Time conjugate grad new =', cclock() - timep
                        deallocate (mat_prep, mat_buf)

                    else
#endif

                        call conjginv(ku_ion, kp_complex, symmagp, dimmat, rank, ranksr, commsr_mpi, mat&
                                &, psip, v1, v2, v4, maxit, tolcg, cost, skdiag, v3, parcut, eps_umrigar)
#ifdef   PARALLEL
                    end if
#endif
                    normsr = sum(psip(n4:n4 + ku_ion - 1)*v3(1:ku_ion))

                    if (signalnoise .and. ieskin .ne. 0 .and. idyn .eq. 0) then
                        if (tion .ne. tjas) then
                            do ii = kn_ion - ieskin + 1, kn_ion - ieskin + ieskinr_pos
                                v3(ii) = v3(ii)*costf
                            end do
                        end if
                        if (tcell .ne. tjas) then
                            do ii = kn_ion - ieskin + ieskinr_pos + 1, kn_ion
                                v3(ii) = v3(ii)*costc
                            end do
                        end if
                    end if

                    if (rank .eq. 0) then
                        time = cclock() - time
                        write (6, *) ' Time inverse Conjugate  Gradients =', time
                    end if

                    if (sum(abs(reduce(1, 1:npr))) .eq. 0.d0) then
                        beta_used = 0.d0
                    else
                        beta_used = beta_learning
                        do i = 1, npr
                            reduce(1, i) = reduce(1, i)*beta_used
                        end do
                    end if

                    !         call dscalzero(npr,0.d0,reduce,ncg)
                    !         if(rank.eq.0) write(6,*) ' Output reduce '
                    do i = 1, kn_ion
                        ii = index(i)
                        if (allowed_par(ii)) reduce(1, ii) = reduce(1, ii) + (1.d0 - beta_used)*v3(i)
                        !         if(rank.eq.0) write(6,*) i,ii,allowed_par(ii),reduce(1,ii)
                    end do
                    !         call mpi_finalize(ierr)
                    !         stop
                    !         Updating collective
                    !         if(rank.eq.0) then
                    !         write(6,*) ' Collective directions ='
                    !         do i=kn_ionp,ku_ion
                    !         write(6,*) i-kn_ion,v3(i)
                    !         enddo
                    !         endif
                    if (ncg_adr .gt. 0) call dgemv('T', ncg_adr, kp0, 1.d0 - beta_used, reducel       &
                            &, ncg_adr, v3(kn_ionp), 1, 1.d0, reduce, ncg)
                    deallocate (mat, v1, v2, v3, v4, okav_coll, econf_coll)
                    if (writescratch .eq. 0) rewind (indopen3)
                end if ! endif ncg>0
            end if ! endif abs(klr)=7

            !  compute the devmax last needed fk  psip(n3) (force components)
            !   if(rank.eq.0) write(6,*) ' force after  kl ',np,sum(abs(psip(n3:n3+np-1)))

#ifdef PARALLEL
            call mpi_bcast(numpar, 1, MPI_INTEGER, 0, commopt_mpi, ierr)
            numparp = numpar + 1
!   For unreliable  networks.
            call bcast_integer(ipip(n3), nbinmax, 0, commopt_mpi)
            call bcast_real(reduce, nbintot, 0, commsr_mpi)
#ifdef UNREL
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
#endif

            if (ncg .gt. 0) then

                !       call bcast_real(reduce,nbintot,0,commsr_mpi)

                call dgemv('N', nbin, npr, 1.d0, fk, dimfk, reduce, ncg               &
                        &, 0.d0, psip(n4), 1)

                cost = psip(n4)**2
                do jj = 1, nbin - 1
                    cost = cost + psip(n4 + jj)**2
                end do

#ifdef PARALLEL
                summpi = cost
#ifdef UNREL
                call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
                call mpi_reduce(summpi, cost, 1, MPI_DOUBLE_PRECISION              &
             &  , MPI_SUM, 0, commsr_mpi, ierr)
#endif
                if (rank .eq. 0 .and. ncg .gt. 0 .and. cost .gt. 0.) then
                    cost0 = ddot(npr, reduce, ncg, forza_sav, 1)
                    write (6, *) ' devmax SR step =', abs(cost0/dsqrt(cost))
                end if

            end if

            !         load matrix sr and update reduce

            if (itest .ne. -5) call dscalzero(nptot, 0.d0, sov(1, 1, 5), 1)
            nps = numpar - ieskin
            npps = nps + 1

            if (itest .ne. -5) then
                time_load = cclock()

                if (allocated(index) .and. ncg .gt. 0 .and. k_par .lt. npk .or. yes_hessc) then
                    allocate (v1(k_par), mat(ncg, k_par))
                    do j = 1, k_par
                        mat(:, j) = reduce(:, index(j))
                    end do
                end if

                if (writescratch .eq. 0) then
                    do kk = 1, ibinit
                        read (indopen3)
                    end do
                end if

                call dscalzero(nptot, 0.d0, sov, 1)
                call dscalzero(nptot, 0.d0, sov(1, 1, 2), 1)
                allocate (factorsr_sav(in1))

                do kk = 1, nweight - ibinit
                    if (writescratch .ne. 0) then
                        countscra = kk + ibinit
                        indscra = (countscra - 1)*nbufscra + 1
                        cost = bufscra(indscra)
                        do i = 1, in1
                            do j = 1, np
                                econf(i, j) = bufscra(indscra + 1)
                                econfh(i, j) = bufscra(indscra + 2)
                                indscra = indscra + 2
                            end do
                            factorsr_sav(i) = bufscra(indscra + 1)
                            if (ipc .eq. 2) then
                                enert(1, i) = bufscra(indscra + 2)
                                enert(2, i) = bufscra(indscra + 3)
                                indscra = indscra + 3
                            else
                                enert(1, i) = bufscra(indscra + 2)
                                indscra = indscra + 2
                            end if
                        end do
                    else
                        read (indopen3) cost, ((econf(i, j), econfh(i, j), j=1, np)            &
                                &, factorsr_sav(i), enert(1:ipc, i), i=1, in1)
                    end if
                    !         reduce econf, econfh
                    !    update reduced sov(1,1,1),sov(1,1,2),sov(1,1,5)

                    do i = 1, in1
                        if (ncg .gt. 0) then
                            if (allocated(index) .and. k_par .lt. npk .or. yes_hessc) then
                                do jj = 1, k_par
                                    v1(jj) = econf(i, index(jj))
                                end do
                                if (yes_hessc) then
                                    call dgemv('N', ncg, k_par, 1.d0, mat, ncg, v1, 1, 0.d0, psip(n3), 2)
                                    do jj = 1, k_par
                                        ind = type_complex(index(jj))
                                        if (ind .eq. 0) then
                                            v1(jj) = 0.d0
                                        elseif (ind .eq. 1) then
                                            v1(jj) = econf(i, index(jj) + 1)
                                        elseif (ind .eq. 2) then
                                            v1(jj) = -econf(i, index(jj) - 1)
                                        end if
                                    end do
                                    call dgemv('N', ncg, k_par, 1.d0, mat, ncg, v1, 1, 0.d0, psip(n3 + 1), 2)
                                else
                                    call dgemv('N', ncg, k_par, 1.d0, mat, ncg, v1, 1, 0.d0, psip(n3), 1)
                                end if
                            else
                                call dgemv('N', ncg, np, 1.d0                         &
                                        &, reduce, ncg, econf(i, 1), in1, 0.d0, psip(n3), 1)
                            end if
                        end if
                        if (yes_hessc) then
                            do jj = ncg + 1, nps
                                psip(n3 + 2*jj - 2) = econf(i, ipip(n3 + jj - 1))
                            end do
                            do jj = ncg + 1, nps
                                ind = type_complex(ipip(n3 + jj - 1))
                                if (ind .eq. 0) then
                                    psip(n3 + 2*jj - 1) = 0.d0
                                elseif (ind .eq. 1) then
                                    psip(n3 + 2*jj - 1) = econf(i, ipip(n3 + jj - 1) + 1)
                                elseif (ind .eq. 2) then
                                    psip(n3 + 2*jj - 1) = -econf(i, ipip(n3 + jj - 1) - 1)
                                end if
                            end do
                            psip(n3 + 2*nps) = 1.d0
                            psip(n3 + 2*nps + 1) = 0.d0
                        else
                            do jj = ncg + 1, nps
                                psip(n3 + jj - 1) = econf(i, ipip(n3 + jj - 1))
                            end do
                            psip(nps + n3) = 1.d0
                        end if
                        if (ncg .gt. 0) then
                            if (allocated(index) .and. k_par .lt. nps .or. yes_hessc) then
                                do jj = 1, k_par
                                    v1(jj) = econfh(i, index(jj))
                                end do
                                if (yes_hessc) then
                                    call dgemv('N', ncg, k_par, 1.d0, mat, ncg, v1, 1, 0.d0, psip(n4), 2)
                                    do jj = 1, k_par
                                        ind = type_complex(index(jj))
                                        if (ind .eq. 0) then
                                            v1(jj) = 0.d0
                                        elseif (ind .eq. 1) then
                                            v1(jj) = econfh(i, index(jj) + 1)
                                        elseif (ind .eq. 2) then
                                            v1(jj) = -econfh(i, index(jj) - 1)
                                        end if
                                    end do
                                    call dgemv('N', ncg, k_par, 1.d0, mat, ncg, v1, 1, 0.d0, psip(n4 + 1), 2)
                                else
                                    call dgemv('N', ncg, k_par, 1.d0, mat, ncg, v1, 1, 0.d0, psip(n4), 1)
                                end if
                            else
                                call dgemv('N', ncg, np, 1.d0                         &
                                        &, reduce, ncg, econfh(i, 1), in1, 0.d0, psip(n4), 1)
                            end if
                        end if
                        if (yes_hessc) then
                            do jj = ncg + 1, nps
                                psip(n4 + 2*jj - 2) = econfh(i, ipip(n3 + jj - 1))
                            end do
                            do jj = ncg + 1, nps
                                ind = type_complex(ipip(n3 + jj - 1))
                                if (ind .eq. 0) then
                                    psip(n4 + 2*jj - 1) = 0.d0
                                elseif (ind .eq. 1) then
                                    psip(n4 + 2*jj - 1) = econfh(i, ipip(n3 + jj - 1) + 1)
                                elseif (ind .eq. 2) then
                                    psip(n4 + 2*jj - 1) = -econfh(i, ipip(n3 + jj - 1) - 1)
                                end if
                            end do
                            psip(2*nps + n4) = enert(1, i)
                            psip(2*nps + n4 + 1) = -enert(2, i) ! We should mantain the conventions adopted.
                        else
                            do jj = ncg + 1, nps
                                psip(n4 + jj - 1) = econfh(i, ipip(n3 + jj - 1))
                            end do
                            psip(nps + n4) = enert(1, i)
                        end if

                        !         now  update matrices
                        costw = cost*factorsr_sav(i)
                        costwe = costw*enert(1, i)
                        if (yes_hessc) then
                            call dsyrk('U', 'T', npps, 2, costw, psip(n3), 2, 1.d0, sov(1, 1, 2), npm)
                        else
                            call dsyr('U', npps, costw, psip(n3), 1, sov(1, 1, 2), npm)
                        end if
                        if (itest .ne. -5) then
                            if (yes_hessc) then
                                !   NB notice that we upload the transponce of the Hamiltonian matrix
                                !   elements here and then we make the transponce again in the line
                                !   identified by the symbol CRUCIAL
                                call dgemm('T', 'N', npps, npps, 2, costw, psip(n4), 2, psip(n3), 2, 1.d0, sov, npm)
                                !    call dsyrk('U','T',npps,2,costwe,psip(n3),2,1.d0,sov(1,1,5),npm)
                                call dgemm('T', 'N', npps, npps, 2, costwe, psip(n3), 2, psip(n3), 2&
                                        &, 1.d0, sov(1, 1, 5), npm)
                                !        Contribution imaginary part
                                costwe = costw*enert(2, i)
                                do j = 0, 2*npps - 1, 2
                                    psip(n4 + j) = -psip(n3 + j + 1)
                                    psip(n4 + j + 1) = psip(n3 + j)
                                end do
                                call dgemm('T', 'N', npps, npps, 2, costwe, psip(n4), 2, psip(n3), 2&
                                        &, 1.d0, sov(1, 1, 5), npm)
                            else
                                call dger(npps, npps, costw, psip(n4), 1, psip(n3), 1, sov, npm)
                                call dsyr('U', npps, costwe, psip(n3), 1, sov(1, 1, 5), npm)
                            end if
                        end if
                    end do
                end do
                deallocate (factorsr_sav)
                if (writescratch .eq. 0) rewind (indopen3)
                if (allocated(index)) then
                    deallocate (index)
                    if (allocated(truecomplex)) deallocate (truecomplex)
                    if (allocated(type_complex)) deallocate (type_complex)
                    if (ncg .gt. 0 .and. k_par .lt. nps .or. yes_hessc) deallocate (v1, mat)
                end if
#ifdef PARALLEL
#ifdef UNREL
                call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#endif
                call reduce_base_real(nptot, sov(1, 1, 2), commrep_mpi, -1)
                if (itest .ne. -5) then
                    call reduce_base_real(nptot, sov(1, 1, 5), commrep_mpi, -1)
                    call reduce_base_real(nptot, sov, commrep_mpi, -1)
                end if
#endif

                cost = 1.d0/sov(npps, npps, 2)

            else
                if (allocated(index)) deallocate (index)
                if (allocated(truecomplex)) deallocate (truecomplex)
                if (allocated(type_complex)) deallocate (type_complex)
            end if

            np = numpar
            npp = numparp

            if (itest .ne. -5) then
                !        matrix-5
                !        call  reduce_mat(nbin,np,npm,sov(1,1,5),reduce,sov(1,1,3))
                do i = 1, npps
                    sov(i, i, 2) = sov(i, i, 2)*cost
                    if (.not. yes_hessc) sov(i, i, 5) = sov(i, i, 5)*cost
                    do j = 1, npps
                        sov(j, i, 1) = sov(j, i, 1)*cost
                    end do
                    if (yes_hessc) then
                        do j = 1, npps
                            sov(j, i, 5) = sov(j, i, 5)*cost
                        end do
                    end if
                    do j = i + 1, npps
                        sov(i, j, 2) = sov(i, j, 2)*cost
                        sov(j, i, 2) = sov(i, j, 2)
                        if (.not. yes_hessc) then
                            sov(i, j, 5) = sov(i, j, 5)*cost
                            sov(j, i, 5) = sov(i, j, 5)
                        end if
                    end do
                end do

            end if ! endif itest.ne.-5

            if (itest .eq. -4) then

                !          now replace the reduced parts of s(*,*,1) and s(*,*,2)
                ! Strictly speaking in the complex case I am not using the semiorthogonal basis
                ! but the original basis O_k |Psi> and O_0|Psi> , O_0=I ,the Real part is taken
                ! after in a way that the strong zero variance property is preserved.
                ! This implies that <O_k> = Re < O_k O_0> --> Re <O_k>.

                do i = 1, nps
                    do j = 1, nps
                        sov(i, j, 1) = sov(i, j, 1) - sov(i, npps, 1)*sov(j, npps, 2) ! only scorrelat
                        sov(i, j, 2) = sov(i, j, 2) - sov(i, npps, 2)*sov(j, npps, 2)
                    end do
                end do

                !         Evaluation of the matrix G

                do i = 1, nps
                    do j = 1, nps
                        sov(i, j, 5) = sov(i, j, 5)                                       &
                                & - sov(npps, i, 2)*sov(npps, j, 5) - sov(npps, j, 2)*sov(i, npps, 5)              &
                                & + sov(npps, npps, 5)*(sov(npps, i, 2)*sov(npps, j, 2) - sov(i, j, 2))
                    end do
                end do

                call dcopy(nptot, sov(1, 1, 2), 1, sov(1, 1, 3), 1)

#ifdef PARALLEL
                if (kaverage) then
                    sov(1:npps, 1:npps, 3) = sov(1:npps, 1:npps, 3)*wkp(ikpoint)
                    call reduce_base_real(nptot, sov(1, 1, 3), commcolrep_mpi, -1)
                elseif (nprocopt .gt. nprocrep) then
!       Take the average according to optimization choice
                    call reduce_base_real(nptot, sov(1, 1, 3), commopt_mpi, -1)
                    sov(1:npps, 1:npps, 3) = sov(1:npps, 1:npps, 3)/nprocopt
                end if
#endif
                time_load = cclock() - time_load
                if (rank .eq. 0) write (6, *) ' Loading matrices linear method =', time_load
                do i = 1, nps
                    indscali = ipip(n3 + i - 1)
                    if (indscali .eq. 0) then
                        psip(i) = -1.d0
                    elseif (indscali .gt. 0) then
                        psip(i) = scalpar(indscali)
                    else
                        psip(i) = 0.d0
                    end if
                    !        write(6,*) i,ipip(i)
                end do

                !        modify the matrix sov for  dynamic (trivial sov)
                if (yes_dgelscut) then

                    call dgelscut(sov(1, 1, 3), nps, npm, psip(n2), epst, ipip, nps          &
                            &, info, psip, lwork, parcut2, ipip(n3), rank, nprocopt, rankopt, commopt_mpi)

                    !       write(6,*) ' ipip after dgelscut '
                    do i = 1, nps
                        itouch(i) = ipip(i + nps)
                        ipip(i) = ipip(i + nps)
                        !       write(6,*) ipip(i+nps)
                    end do

                else

                    do i = 1, nps
                        itouch(i) = 1
                        ipip(i) = 1
                        !       write(6,*) ipip(i+nps)
                    end do

                end if

                betap = beta + 1.d0

                !          no symm with zero variance property for the right eigenvector

                do i = 1, nps
                    do j = 1, nps
                        !   CRUCIAL  Notice we interchange the indices of sov(i,j,1) see CRUCIAL above
                        sov(i, j, 3) = sov(j, i, 1) + sov(i, j, 5)
                    end do
                end do

                do i = 1, nps
                    do j = i + 1, nps
                        !  NB the scalar product is instead symmetric within numerical accuracy.
                        sov(i, j, 2) = 0.5d0*(sov(i, j, 2) + sov(j, i, 2))
                        sov(j, i, 2) = sov(i, j, 2)
                    end do
                end do

                if (yes_real) then

                    do i = 1, nps
                        sov(npps, i, 5) = 0.5d0*sov(i, npps, 1) + &
                                &    sov(npps, i, 1) - sov(npps, npps, 1)*sov(npps, i, 2)
                        sov(i, npps, 5) = sov(npps, i, 5)
                    end do

                    ! symmetrization also other part, as there is no strong zero variance property.
                    do i = 1, nps
                        do j = i + 1, nps
                            sov(i, j, 3) = 0.5d0*(sov(i, j, 3) + sov(j, i, 3))
                            sov(j, i, 3) = sov(i, j, 3)
                            sov(j, i, 1) = sov(i, j, 1)
                        end do
                    end do

                else

                    !   Below is the calculation of matrix elements of H_kk' connecting k or k'=0
                    !   namely the calculation of energy derivatives.

                    do i = 1, nps
                        sov(npps, i, 5) = sov(i, npps, 1) + &
                                &    sov(npps, i, 1) - sov(npps, npps, 1)*sov(npps, i, 2)
                        sov(i, npps, 5) = sov(npps, i, 1) - sov(npps, npps, 1)*sov(npps, i, 2)
                    end do

                end if

                do i = 1, nps
                    sov(i, npps, 2) = 0.d0
                    sov(npps, i, 2) = 0.d0
                end do

                do i = 1, nps
                    do j = 1, nps
                        sov(i, j, 5) = sov(i, j, 3)
                    end do
                end do

                sov(npps, npps, 5) = 0.d0
#ifdef PARALLEL
                if (kaverage) then
!NB when kaverage the average of matrices H_kk' in this way
!(after the use of the semiorthogonal basis) implies correctly the calculation of the average
! energy derivative that should vanish at the minimum average energy.
                    sov(1:npps, 1:npps, 1:5) = sov(1:npps, 1:npps, 1:5)*wkp(ikpoint)
                    call reduce_base_real(5*nptot, sov, commcolrep_mpi, -1)
                elseif (nprocopt .gt. nprocrep) then
                    call reduce_base_real(5*nptot, sov, commopt_mpi, -1)
                    sov(:, :, 1:5) = sov(:, :, 1:5)/nprocopt
                end if
#endif

                ipip(npps) = 1
                indi = 0
                do i = 1, npps
                    if (ipip(i) .eq. 1) then
                        indi = indi + 1
                        indj = 0
                        do j = 1, npps
                            if (ipip(j) .eq. 1) then
                                indj = indj + 1
                                sov(indi, indj, 3) = sov(i, j, 5)
                                sov(indi, indj, 4) = sov(i, j, 2)
                            end if
                        end do
                    end if
                end do
                if (parcute .gt. 0) then
                    !     add a mass to the Hamiltonian in the direction orthogonal to Psi_VMC
                    indi = 0
                    do i = 1, nps
                        if (ipip(i) .eq. 1) then
                            indi = indi + 1
                            indj = 0
                            do j = 1, nps
                                if (ipip(j) .eq. 1) then
                                    indj = indj + 1
                                    sov(indi, indj, 3) = sov(indi, indj, 3) + parcute*sov(i, j, 2)
                                end if
                            end do
                        end if
                    end do
                    indi = indi + 1 ! restore previous value of indi
                end if

                !     diagonalize the overlap matrix

                if (nproc_diag .gt. 1) then
                    call dsyev_my('V', 'L', indi, sov(1, 1, 4), npm, psip, info&
          &, nprocopt, rankopt, commopt_mpi)
                else
                    call dsyev('V', 'L', indi, sov(1, 1, 4), npm, psip, psip(n2), lwork, info)
                end if
                if (info .ne. 0) then
                    if (rank .eq. 0) write (6, *) ' Error in lapack dsyev !!! ', rankopt

                    !       stop
                    !#ifdef PARALLEL
                    !       call mpi_finalize(ierr)
                    !#endif
                    errnoise = -1

                else

                    !      now calculate normalized directions

                    if (yes_dgelscut) then

                        do i = 1, indi
                            if (psip(i) .gt. 0.d0) then
                                cost = 1.d0/dsqrt(psip(i))
                                call dscal(indi, cost, sov(1, i, 4), 1)
                            else
                                call dscalzero(indi, 0.d0, sov(1, i, 4), 1)
                            end if
                        end do
                        indj = indi
                    else
                        if (rank .eq. 0) write (6, *) ' Lowest/largest eig SR mat =', psip(1), psip(indi)
                        indj = indi
                        indi = 0
                        do i = 1, indj
                            if (psip(i) .gt. epst) then
                                indi = indi + 1
                                cost = 1.d0/dsqrt(psip(i))
                                if (indi .ne. i) call dcopy(indj, sov(1, i, 4), 1, sov(1, indi, 4), 1)
                                call dscal(indj, cost, sov(1, indi, 4), 1)
                            end if
                        end do
                        if (rank .eq. 0 .and. indj - indi .gt. 0) write (6, *) ' Warning eliminated small eig up to', indj - indi
                    end if

                    !      now change basis for sov(1,1,3)

                    call dgemm('N', 'N', indj, indi, indj, 1.d0, sov(1, 1, 3), npm            &
                            &, sov(1, 1, 4), npm, 0.d0, sov(1, 1, 5), npm)

                    call dgemm('T', 'N', indi, indi, indj, 1.d0, sov(1, 1, 4), npm            &
                            &, sov(1, 1, 5), npm, 0.d0, sov(1, 1, 3), npm)

                    !      Now compute right eigenvectors

                    call dgeev('N', 'V', indi, sov(1, 1, 3), npm, psip, psip(n2)              &
                            &, dummy, npm, sov(1, 1, 5), npm, psip(n3), lwork, info)

                    ncoll = min(indi - 1, ncg)

                    if (info .ne. 0) then
                        write (6, *) ' Error in lapack dgeev !!! ', rankopt
                        errnoise = -2
                    end if

                end if ! endif info ne 0

                if (errnoise .eq. 0) then

                    !     now go back into the previous basis

                    call dcopy(nptot, sov(1, 1, 5), 1, sov(1, 1, 3), 1)

                    call dgemm('N', 'N', indj, indi, indi, 1.d0, sov(1, 1, 4), npm            &
                            &, sov(1, 1, 3), npm, 0.d0, sov(1, 1, 5), npm)

                    do i = indi + 1, indj
                        sov(:, i, 5) = 0.d0
                        psip(i) = 0.d0
                        psip(n2 - 1 + i) = 0.d0
                    end do

                    indi_save = indi
                    !       read again sov
                    indi = 0
                    do i = 1, npps
                        if (ipip(i) .eq. 1) then
                            indi = indi + 1
                            indj = 0
                            do j = 1, npps
                                if (ipip(j) .eq. 1) then
                                    indj = indj + 1
                                    sov(indi, indj, 4) = sov(i, j, 2)
                                end if
                            end do
                        end if
                    end do

                    !       rescaling eigenvector
!                   if(rank.eq.0) write(6, *) ' Eigenvalues, Real , Imag ,  weight old  '
                    do i = 1, indi
                        if (sov(indi, i, 5) .ne. 0.d0) then
                            cnorm = 0.d0
                            do j = 1, indi - 1
                                do k = 1, indi - 1
                                    cnorm = cnorm + sov(j, k, 4)*sov(j, i, 5)*sov(k, i, 5)
                                end do
                            end do
                            if (cnorm .ne. 0.d0) then
                                cnorm = 1.d0/dsqrt(cnorm)
                                call dscal(indi, cnorm, sov(1, i, 5), 1)
                                cost = sov(indi, i, 5)**2
                                cost = cost/(cost + 1.d0)
                            else
                                cost = 1.d0
                            end if
                        else
                            cost = 0.d0
                        end if
!                       if(rank.eq.0) write(6, 124) i, psip(i), psip(n2 - 1 + i), cost
                    end do
124                 format(I12, 3f15.8)

                    !       Choose the closest real eigenvalue < 0= reference energy
                    imin = 0
                    eigmin = 0.d0
                    do i = 1, indi
                        !       the gain in energy / norm orth component is maximum
                        cost = -psip(i)*abs(sov(indi, i, 5))
                        if ((cost .gt. eigmin .or. imin .eq. 0)                               &
                                &.and. psip(i) .lt. 0.d0) then
                            !        if((cost.gt.eigmin.or.imin.eq.0).and.psip(n2-1+i).eq.0.d0
                            imin = i
                            eigmin = cost
                        end if
                    end do

                    !       stop

                    if (imin .eq. 0 .and. rank .eq. 0) write (6, *) ' Warning all positive eigenvalues '
                    !       Then no constraint on energy gain, choose the closest real

                    over2 = sov(indi, max(imin, 1), 5)**2
                    over2 = over2/(1.d0 + over2)

!                       if(rank.eq.0) write(6, *) ' Warning  Chosen the closest '
                    eigmin = 0.d0
                    do i = 1, indi
                        cost = sov(indi, i, 5)**2/(1.d0 + sov(indi, i, 5)**2)
                        if (cost .gt. eigmin .or. imin .eq. 0) then
                            imin = i
                            eigmin = cost
                        end if
                    end do
                    over2 = eigmin

                    eigmin = psip(imin)

                    if(rank.eq.0) write(6, *) 'Chosen  Eigenvalue (new-previous)'//&
                       &' #, Real/Ima3g, dimension  H', imin, psip(imin) / 2.d0, psip(n2-1+imin)/2.d0,indi_save

                    do i = 1, indi - 1
                        psip(i + n2 - 1) = sov(i, imin, 5)/sov(indi, imin, 5)
                    end do

!                    if((1.d0 - over2).gt.-epsi.and.epsi.lt.0.d0) then
!                        if(rank.eq.0) write(6, *) ' Warning decelerated move too much large !!! '
!                        cost = -epsi / (1.d0 - over2)
!                        do i = 1, indi - 1
!                            psip(i + n2 - 1) = cost * psip(i + n2 - 1)
!                            !        write(6,*) i,psip(i+n2-1)
!                        enddo
!
!                    endif

                    !      now simple steepest descent  for the parameters without hessian
                    indi = indi - 1
                    indin = indi

                    if (idyn .gt. 0 .and. rank .eq. 0) write (6, *) ' Ion forces '
                    do i = npps, np
                        !        if(ipip(i).eq.2) then
                        indin = indin + 1
                        ! a factor two is needed
                        !        psip(indin)=-sov(i,npp,5)*2.d0
                        psip(indin) = forza(i)
                        if (idyn .gt. 0 .and. rank .eq. 0) write (6, *) i, forza(i), err(i)
                        !        endif
                    end do
                    ndsr = indin - indi
                    indin = indi + 1

                    !        write(6,*) ' ndsr FOUND =',ndsr

                    !        save cov for output before it is destroyed by ion_dyn
                    if (write_cov) then
                        cov_sav(1:ieskin*ieskin) = cov(1:ieskin*ieskin)
                    end if
#ifdef PARALLEL
                    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
                    if (rank .eq. 0) write (6, *) ' Before dynamic ', info, idyn
                    if (ndsr .ne. 0) call ion_dynamics
                    if (rank .eq. 0) write (6, *) ' After dynamic ', info

                    if (write_cov .and. acc_dyn .and. rank .eq. 0) then
                        write (18, 123) ((ris(4)*cov_sav(ii + (jj - 1)*ieskin), ii=jj, ieskin), jj=1, ieskin)
                    end if

                    !         write(6,*) ' ipip, Solution found '
                    !         write(6,*) ' Normal parameters '

                    indi = 0
                    do i = 1, nps
                        if (ipip(i) .eq. 1) then
                            indi = indi + 1
                            psip(i) = psip(n2 + indi - 1)
                        else
                            psip(i) = 0.d0
                        end if
                        !         write(6,*) i,ipip(i),psip(i)
                    end do

                    indin = indi
                    if (idyn .gt. 0 .and. rank .eq. 0) write (6, *) ' Forces '
                    do i = npps, np
                        !         if(ipip(i).eq.2) then
                        indi = indi + 1
                        psip(i) = psip(n2 + indi - 1)
                        !         endif
                        if (idyn .gt. 0 .and. rank .eq. 0) write (6, *) i, psip(i)
                    end do

                    !         evaluation norm correction

                    cost = 0.d0
                    do i = 1, indin
                        do j = 1, indin
                            cost = cost + psip(n2 + i - 1)*psip(n2 + j - 1)*sov(j, i, 4)
                        end do
                    end do

                    cost = dsqrt(cost)
                    costn = 1.d0

                    if (cost .gt. abs(epsi) .and. epsi .ne. 0) then
                        costn = abs(epsi)/cost

                        if (rank .eq. 0) write (6, *) ' Warning decelerated ', costn
                        call dscal(nps, costn, psip, 1)
                    end if

                    if (epsi .lt. 0 .and. cost*1.5 .lt. -epsi) then
                        epsi = -1.5*cost
                        if (rank .eq. 0) write (6, *) 'Warning changing epsi on the fly =', epsi
                    end if
                    if (rank .eq. 0) write (6, *) ' Norm correction = ', costn*cost
                    norm_corr = costn*cost

                    if (np .gt. ncg) then

                        if ((abs(klr) .eq. 6 .or. abs(klr) .eq. 7 .or. abs(klr) .eq. 2) .and. ncoll .eq. 1 .and. rank .eq. 0)&
                                & write (6, *) ' Dt found =', psip(1)*tjas
                        scale_grad = psip(1)

                        if (ncoll .eq. 1 .and. psip(1) .lt. 0 .and. npbra .eq. 0) then
                            if (rank .eq. 0) write (6, *) ' Move rejected negative Dt !!! '
                            if (.not. change_tpar) psip(1) = 0.d0
                        end if

                        if (ncg .gt. 0) then
                            call dgemv('T', ncg, npr, 1.d0, reduce, ncg, psip, 1, 0.d0, psip(n2), 1)
                        else
                            call dscalzero(npr, 0.d0, psip(n2), 1)
                        end if
                        !        write(6,*) ' Updated old corr ',np
                        do jj = ncg + 1, np
                            indscali = ipip(n3 + jj - 1)
                            !       write(6,*) ' before =',jj,indscali,psip(n2+indscali-1),psip(jj)
                            psip(n2 + indscali - 1) = psip(n2 + indscali - 1) + psip(jj)
                        end do

                    else

                        if ((abs(klr) .eq. 6 .or. abs(klr) .eq. 7 .or. abs(klr) .eq. 2) .and. ncoll .eq. 1 .and. rank .eq. 0)&
                                &  write (6, *) ' Dt found =', psip(1)*tjas
                        scale_grad = psip(1)

                        if (ncoll .eq. 1 .and. psip(1) .lt. 0) then
                            if (rank .eq. 0) write (6, *) ' Move rejected negative Dt !!! '
                            if (.not. change_tpar) psip(1) = 0.d0
                        end if
                        call dgemv('T', np, npr, 1.d0, reduce, ncg, psip, 1, 0.d0, psip(n2), 1)
                    end if

                    call dcopy(npr, psip(n2), 1, psip, 1)

                    np = npr
                    npp = npr + 1

                end if

            elseif (itest .eq. -5) then

                info = 0

                numpar = 1
                if (signalnoise .and. idyn .eq. 0) np = 1 ! ionic forces are not treated independently
                ipip(1) = 1
                psip(1) = 1.d0

                if (ieskin .ne. 0 .and. idyn .gt. 0) then
                    indi = numpar
                    indin = indi
                    !         store the solution found in psip(n2)
                    call dcopy(np, psip, 1, psip(n2), 1)
                    if (rank .eq. 0) write (6, *) ' Ion forces '
                    do i = npps, np
                        ipip(i) = 2
                        indin = indin + 1
                        ! the io
                        psip(indin) = forza(i)
                        if (rank .eq. 0) write (6, *) i, psip(indin), err(i)
                    end do
                    ndsr = indin - indi
                    if (ndsr .ne. 0) call ion_dynamics

                    !         write(6,*) ' Solution found '

                    indi = 0
                    do i = 1, nps
                        if (ipip(i) .ne. 0) then
                            indi = indi + 1
                            psip(i) = psip(n2 + indi - 1)
                        else
                            psip(i) = 0.d0
                        end if
                        !         write(6,*) i,ipip(i),psip(i)
                    end do
                    !         write(6,*) ' after loop ',indi

                    indin = indi
                    !         write(6,*) ' Normal parameters '
                    do i = npps, np
                        !         write(6,*) ' after ipip,psip',i,ipip(i),psip(n2+i-1)
                        !         if(ipip(i).eq.2) then
                        indi = indi + 1
                        psip(i) = psip(n2 + indi - 1)
                        !         endif
                        !         write(6,*) i,psip(i)
                    end do
                end if

                !       if(ncg.gt.1) then

                !        cost=0.d0
                !         do i=1,nps
                !           do j=1,nps
                !           cost=cost+psip(i)*psip(j)*sov(j,i,1)
                !           enddo
                !         enddo
                !
                !
                !       cost=dsqrt(cost)

                cost = tjas*dsqrt(normsr)
                if (.not. signalnoise) then
                    if (rank .eq. 0) write (6, *) ' Norm correction = ', cost
                    norm_corr = cost

                    if (cost .gt. abs(epsi) .and. epsi .ne. 0) then
                        if (rank .eq. 0) write (6, *) ' Warning decelerated by ', epsi/cost
                        cost = abs(epsi)/cost
                        psip(1) = cost
                    end if
                    if (epsi .lt. 0 .and. norm_corr*3 .lt. -epsi) then
                        epsi = -3*norm_corr
                        if (rank .eq. 0) write (6, *) 'Warning changing epsi on the fly =', epsi
                    end if
                end if
                !      endif

                !       if(parcut.gt.0) then
                if (rank .eq. 0) write (6, *) ' I consider # ', numpar, 'independent param. '

                if (np .gt. ncg) then
                    if (ncg .gt. 0) then
                        call dgemv('T', ncg, npr, 1.d0, reduce, ncg, psip, 1, 0.d0, psip(n2), 1)
                    else
                        call dscalzero(npr, 0.d0, psip(n2), 1)
                    end if

                    do jj = ncg + 1, np
                        indscali = ipip(n3 + jj - 1)
                        !        write(6,*) ' before =',jj,indscali,psip(n2+indscali-1),psip(jj)
                        psip(n2 + indscali - 1) = psip(n2 + indscali - 1) + psip(jj)
                    end do
                else
                    call dgemv('T', np, npr, 1.d0, reduce, ncg, psip, 1, 0.d0, psip(n2), 1)
                end if
                call dcopy(npr, psip(n2), 1, psip, 1)
                if (idyn .eq. 0 .and. ieskin .ne. 0 .and. signalnoise) then
                    cost = tjas*dnrm2(ieskin, psip(npr - ieskin + 1), 1)
                    if (rank .eq. 0) write (6, *) ' Norm change ions =', cost
                    if (cost .gt. abs(epsi) .and. epsi .ne. 0) then
                        if (rank .eq. 0) write (6, *) ' Warning decelerated ions by ', epsi/cost
                        cost = abs(epsi)/cost
                        psip(npr - ieskin + 1:npr) = cost*psip(npr - ieskin + 1:npr)
                    end if
                end if
                !        write(6,*) ' Direction chosen '
                !        do i=1,npr
                !        write(6,*) i,psip(i)
                !        enddo
                np = npr
                npp = npr + 1
                !       endif
            end if

            if (info .ne. 0) then
                if (rank .eq. 0) write (6, *) ' Warning info =', info
                if (rank .eq. 0) write (6, *) 'Continuing with previous param. !!!'
                call dscalzero(np, 0.d0, psip, 1)
            end if

            !       Now rescaling of the solution found
            !if(rank.eq.0)           write(6,*) ' Correction found '
            do j = 1, np
                alpha(j) = psip(j)
                !if(rank.eq.0)           write(6,*) j,psip(j)
            end do
            !          stop

            alpha(npp) = 1.d0

            !       endif ! rank.eq.0

            !  This change has to be put after the if rank==0, otherwise only the
            !  master makes the change.
            if (idyn .gt. 0 .and. stepcg_recount) then
                stepcg = 0
                reduce(1:ncg, 1:np) = 0.d0
            end if
#ifdef PARALLEL
            call mpi_bcast(errnoise, 1, MPI_INTEGER, 0, commopt_mpi, ierr)
#endif
            if (errnoise .ne. 0) then

                if (rank .eq. 0) write (6, *) ' Error in reweight0  =', errnoise

#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif

                stop

            end if

#ifdef PARALLEL
            if (.not. yesquantum) then
                call bcast_real(alpha, npr, 0, commsr_mpi) ! the maximum communication as commsr_mpi>=commopt_mpi
            else
                call bcast_real(alpha, npk, 0, commopt_mpi)
                if (yescomm) then
                    call bcast_real(alpha(npk + 1), ieskin, 0, commrep_mpi)
                end if
            end if ! endif yesquantum
!$omp barrier
#endif

            if (idyn .gt. 2) deallocate (vel0)
            if (allocated(forza_sav)) then
                call dcopy(npr, forza_sav, 1, forza, 1)
                deallocate (forza_sav)
            end if

        end if ! fine if iweight

    end if ! fine if np > 0

    return

#ifdef __KCOMP
123 format(32767e15.7)
#else
123 format(1000000e15.7)
#endif

contains

    subroutine update_index
        implicit none
        kp_complex = 0
        if (srcomplex .and. .not. signalnoise) then
            allocate (index(kp_ion), truecomplex(kp_ion), type_complex(npk))
            truecomplex = .false. ! Only the first real element of the complex
            type_complex = 0
            !           type_complex=0   real parameter
            !           type_complex=1   real or imaginary parameter of the normal type
            !         namely  the associated complex operator O_R +i O_I is determined by
            !    type_complex(ind) -->     econf(ind)-i econf(ind+1)
            !           type_complex=2   imaginary parameter constructed with Cauchy
            !    type_complex(ind) -->     econf(ind)+i econf(ind-1)
            index = 0
            ku_ion = 0
            if (symmagp .and. yes_correct) then
                !           First the complex parameters
                do i = ndimjp, ndims, 4
                    type_complex(i) = 1
                    type_complex(i + 2) = 1
                    if (rpar(i) .eq. 0 .and. allowed_par(i) .and. allowed_par(i + 2)) then
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i
                        !            if(rank.eq.0) write(6,*) ' touched parameter =',i
                        truecomplex(ku_ion) = .true.
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i + 1 ! imaginary part
                        ku_ion = ku_ion + 1
                        truecomplex(ku_ion) = .true.
                        index(ku_ion) = i + 2
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i + 3 ! imaginary part
                        !            if(rank.eq.0) write(6,*) ' touched parameter =',i
                    end if
                end do
                do i = ndimsp, kp0, 2
                    type_complex(i) = 1
                    if (rpar(i) .eq. 0 .and. allowed_par(i)) then
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i
                        truecomplex(ku_ion) = .true.
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i + 1 ! Imaginary part
                        !            if(rank.eq.0) write(6,*) ' touched parameter =',i
                    end if
                end do
                !            Now the modulus  one determined by sjbradet, again with real and imag part for the forces.
                do i = ndimjp, ndims, 4
                    if (rpar(i) .eq. 0 .and. allowed_par(i) .and. .not. allowed_par(i + 2)) then
                        ku_ion = ku_ion + 1
                        truecomplex(ku_ion) = .true.
                        index(ku_ion) = i
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i + 1 ! the imaginary part
                        !            if(rank.eq.0) write(6,*) ' touched parameter =',i
                    end if
                end do
            else
                do i = ndimjp, kp0, 2
                    type_complex(i) = 1
                    type_complex(i + 1) = 2
                    if (rpar(i) .eq. 0 .and. allowed_par(i)) then
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i
                        !            if(rank.eq.0) write(6,*) ' touched parameter =',i
                        truecomplex(ku_ion) = .true.
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i + 1 ! the imaginary part is always needed
                        !            if(rank.eq.0) write(6,*) ' touched parameter =',i
                    end if
                end do
            end if
            kp_complex = ku_ion ! The last complex parameter
            do i = 1, ndimj
                if (rpar(i) .eq. 0 .and. allowed_par(i)) then
                    ku_ion = ku_ion + 1
                    index(ku_ion) = i
                    !           if(rank.eq.0) write(6,*) ' touched parameter =',i
                end if
            end do
        else
            if (signalnoise .and. idyn .eq. 0) kp_ion = npr
            allocate (index(kp_ion))
            allocate (type_complex(npk))
            type_complex = 0 ! simple real
            !          only the parameters allowed by the cutoff are moved
            index = 0
            ku_ion = 0
            if (ipc .eq. 2) then
                if (symmagp .and. yes_correct) then

                    do i = ndimjp, ndims, 4
                        if (.not. signalnoise) then
                            type_complex(i) = 1
                            type_complex(i + 2) = 1
                        end if
                        if (rpar(i) .eq. 0) then
                            if (allowed_par(i)) then
                                ku_ion = ku_ion + 1
                                index(ku_ion) = i
                            end if
                            if (allowed_par(i + 2)) then
                                ku_ion = ku_ion + 1
                                index(ku_ion) = i + 2
                            end if
                        end if
                    end do
                    do i = ndimsp, kp0, 2
                        if (.not. signalnoise) then
                            type_complex(i) = 1
                        end if
                        if (rpar(i) .eq. 0 .and. allowed_par(i)) then
                            ku_ion = ku_ion + 1
                            index(ku_ion) = i
                        end if
                    end do
                else
                    do i = ndimjp, kp0, 2
                        if (.not. signalnoise) then
                            type_complex(i) = 1
                            type_complex(i + 1) = 2
                        end if
                        if (rpar(i) .eq. 0) then
                            if (allowed_par(i)) then
                                ku_ion = ku_ion + 1
                                index(ku_ion) = i
                            end if
                            if (allowed_par(i + 1)) then
                                ku_ion = ku_ion + 1
                                index(ku_ion) = i + 1
                            end if
                        end if
                    end do
                end if
                do i = 1, ndimj
                    if (rpar(i) .eq. 0 .and. allowed_par(i)) then
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i
                        !           if(rank.eq.0) write(6,*) ' touched parameter =',i
                    end if
                end do

            else
                do i = 1, kp0
                    if (rpar(i) .eq. 0 .and. allowed_par(i)) then
                        ku_ion = ku_ion + 1
                        index(ku_ion) = i
                    end if
                end do
            end if
        end if
        if (signalnoise .and. idyn .eq. 0) then
            do i = kp0 + 1, kp0 + ieskin
                if (allowed_par(i)) then
                    ku_ion = ku_ion + 1
                    index(ku_ion) = i
                else
                    if (rank .eq. 0) write (6, *) ' Warning untouched ion parameter =', i
                end if
            end do
        end if

        k_par = ku_ion
        ku_ionp = ku_ion + 1
        if (rank .eq. 0) then
            if (symmagp .or. kp_complex .eq. 0) then
                write (6, *) ' Leading dimension matrix =', ku_ion - kp_complex/2
            else
                if (real_agp) then
                    write (6, *) ' Leading dimension matrix =', ku_ion - kp_complex/2
                else
                    write (6, *) ' Leading dimension matrix =', ku_ion
                end if
            end if
        end if

#ifdef PARALLEL
#ifdef UNREL_DIAG
        if (srcomplex .and. .not. signalnoise) then
            call mpi_bcast(truecomplex, kp_ion, MPI_LOGICAL, 0, commsr_mpi, ierr)
        end if
        call mpi_bcast(index, kp_ion, MPI_INTEGER, 0, commsr_mpi, ierr)
        call mpi_bcast(kp_complex, 1, MPI_INTEGER, 0, commsr_mpi, ierr)
        call mpi_bcast(type_complex, npk, MPI_INTEGER1, 0, commsr_mpi, ierr)
#endif
#endif
    end subroutine update_index

    subroutine ion_dynamics

        ! by E. Coccia (18/1/11)
        use extpot, only: ext_pot

        implicit none
        real*8 drand1, dnrm2, costm1, costhm1, trace, ratio_dyn, all_dyn
        integer indm, i, j, ind, i_ion, k_ion
        real*8, dimension(:, :), allocatable :: sov4, sov5
        !       BEGIN ION DYNAMICS

        !        write(6,*) ' sov =',ndsr,indin,sov(1,1,4),psip(indin)

        !         call dgetrf(ndsr,ndsr,sov(indin,indin,4),npm,ipip(n4),info)

        !         call dgetrs('N',ndsr,1,sov(indin,indin,4),npp,ipip(n4)
        !     1,psip(indin),ndsr,info)

        !         rescaling the solution
        !         write(6,*) ' sr inside '
        ratio_dyn = 0.d0
        all_dyn = 0.d0

        !         write(6,*) ' DTR/Temp  inside =',dtr,temp
        if (idyn .ge. 2 .and. idyn .le. 4 .or. idyn .eq. 6 .or. idyn .eq. 7 .or. idyn .eq. 8) then
            allocate (sov4(npp, 2))
            sov4 = 0.d0
        end if
        if (yesturboq .or. idyn .eq. 7) then
            allocate (sov5(ieskin, 7))
            sov5 = 0.d0
        elseif (idyn .eq. 5 .and. yesquantum) then
            allocate (sov5(ieskin, 5))
            sov5 = 0.d0
        end if
        n3n = 3*nion
        !    Reinitialize acc_dyn only when dynamic step is done (not when idyn<0).
        if (idyn .ge. 0) acc_dyn = .true.
        ! move the ions
        if (idyn .ge. 0 .and. dtr .ne. 0.d0 .and. (maxdev_dyn .eq. 0 .or. devmaxp .lt. maxdev_dyn)) then
            dt = dtr
            if (idyn .le. 1) then
                !       do k=1,3
                !         cost=0.d0
                !         indin=indi
                !         do i=1,npp
                !         if(ipip(i).eq.2) then
                !         indin=indin+1
                !         ind=ipip(n3+i-1)-kp_ion
                !         write(6,*) ' index force =',i,ind
                !         if(mod(ind-1,3)+1.eq.k) cost=cost+psip(indin)
                !         endif
                !         enddo
                !       write(6,*) ' sum rule forces ',cost
                !       enddo

                indin = indi
                do j = npps, np
                    !         if(ipip(j).eq.2) then
                    indin = indin + 1
                    psip(n2 + indin - 1) = psip(indin)*scalpar(ipip(n3 + j - 1))
                    !        write(6,*) indin,psip(indin),scalpar(ipip(n3+j-1)),ipip(n3+j-1)
                    if (temp .ne. 0.d0) then
                        zeta = dsqrt(-4.d0*scalpar(ipip(n3 + j - 1))*dlog(1.d0 - drand1())*temp)*dcos(2.d0*pi*drand1())
                        psip(n2 + indin - 1) = psip(n2 + indin - 1) + zeta
                    end if
                    !         write(6,*) ' new displacement ',indin,psip(n2+indin-1)
                    !         endif
                end do
                cost = dnrm2(ieskin, psip(n2 + indi), 1)
                if (rank .eq. 0) write (6, *) ' Norm change ions =', cost

                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                ! zero temperature dynamics along direction of maximum signal/noise direction
            elseif (idyn .eq. 5) then

                lworkd = 3*ieskin
                allocate (work(lworkd), eig(ieskin))

                ! by E. Coccia (18/1/11): allow rototranslations
                ! if an external potential is added
                niont = 0
                do i = 1, nion
                    if (atom_number(i) .gt. 0) niont = niont + 1
                end do
                rioncm = 0.d0
                do i = 1, nion
                    if (atom_number(i) .gt. 0) then
                        do j = 1, 3
                            rioncm(j) = rioncm(j) + rion(j, i)/niont
                        end do
                    end if
                end do
                if (.not. ext_pot .and. .not. iespbc) then
                    call cleancov
                elseif (ieskin .ge. 3*niont) then
                    mineig = 4
                else
                    mineig = 1
                end if

                if (ieskin .eq. 3*niont .and. .not. yesquantum) then
                    !         Estimation effective temperature  with virial theorem
                    ind = indi
                    cost = 0.d0
                    if (iespbc) then
                        Picost = TWO_PI/(1.d0 - 1.d0/niont)
                        weight_vir = 0.d0
                        do i = 1, niont
                            do j = 1, 3
                                ind = ind + 1
                                cost = cost - cellscale(j)/Picost*sin(Picost*(rion(j, i) - rioncm(j))/cellscale(j))*psip(ind)
                                weight_vir = weight_vir + cos(Picost*(rion(j, i) - rioncm(j))/cellscale(j))
                            end do
                        end do
                        cost = cost/weight_vir
                    else
                        do i = 1, niont
                            do j = 1, 3
                                ind = ind + 1
                                cost = cost - (rion(j, i) - rioncm(j))*psip(ind)
                            end do
                        end do
                        cost = cost/(ieskin - 3)
                    end if
                    Tmes = cost*ris(2) ! From Ry to H
                    if (rank .eq. 0) write (6, *) ' Temperature (H)/weight= ', Tmes, weight_vir
                elseif (.not. yesquantum) then
                    ! assumed that some atom is fixed
                    cost = 0.d0
                    ind = indi
                    do kk = 1, ieskinr
                        if (ion_table(kk)%mult .eq. 1) then
                            k_ion = abs(ion_table(kk)%ion(1))
                            if (atom_number(k_ion) .gt. 0) then
                                ind = ind + 1
                                i_ion = ion_table(kk)%comp(1)
                                cost = cost - rion(i_ion, k_ion)*psip(ind)
                            end if
                        end if
                    end do
                    if (ind .gt. indi) then
                        cost = cost/(ind - indi)
                        Tmes = cost*ris(2) ! From Ry to H
                        if (rank .eq. 0) write (6, *) ' Temperature (H)= ', Tmes
                    end if
                end if

                indin = 0
                do j = npps, np
                    indin = indin + 1
                    psip(n3 + indin - 1) = dsqrt(scalpar(ipip(n3 + j - 1))/dt)
                end do
                do j = 1, ieskin
                    do i = 1, ieskin
                        cov(ieskin*(j - 1) + i) = cov(ieskin*(j - 1) + i)*psip(n3 + i - 1)*psip(n3 + j - 1)
                    end do
                    if (psip(n3 + j - 1) .eq. 0.d0) cov(ieskin*(j - 1) + j) = 1.d0
                    psip(indi + j) = psip(indi + j)*psip(n3 + j - 1)
                end do
                if (yesquantum) then
                    !         Transform the coordinates also
                    do j = 1, nion
                        do i = 1, 3
                            ind = i + (j - 1)*3
                            sov5(ind, 1) = rion(i, j)/psip(n3 + ind - 1)
                        end do
                    end do
                end if
                !        Regularization
                if (epstion .gt. 0.d0) then
                    do i = 1, ieskin
                        cov((i - 1)*ieskin + i) = cov((i - 1)*ieskin + i)*(1.d0 + epstion)
                    end do
                elseif (epstion .lt. 0.d0) then
                    trace = 0.d0
                    do i = 1, ieskin
                        trace = trace + cov((i - 1)*ieskin + i)
                    end do
                    trace = trace/ieskin
                    do i = 1, ieskin
                        cov((i - 1)*ieskin + i) = cov((i - 1)*ieskin + i) - epstion*trace
                    end do
                end if

                if (delta0 .eq. 0.d0) then
                    dtstep = dtr
                    dtnoise = dtr
                else
                    !  Correction time step, that is exact for an harmonic potential.
                    !  Consistent in the limit dtr-->0
                    if (yesquantum) then
                        if (dtr .gt. 1d-6) then
                            dtstep = (1.d0 - exp(-dtr))
                        else
                            dtstep = dtr - dtr**2/2.d0
                        end if
                        if (2.d0*dtr .gt. 1d-6) then
                            dtnoise = (1.d0 - exp(-2.d0*dtr))/2.d0
                        else
                            dtnoise = dtr - dtr**2
                        end if
                    else
                        !   Here delta0 is the best approximation of the hessian H via H=delta0 x cov
                        if (dtr*delta0 .gt. 1d-6) then
                            dtstep = (1.d0 - exp(-delta0*dtr))/delta0
                        else
                            dtstep = dtr - dtr**2*delta0/2.d0
                        end if
                        if (2.d0*dtr*delta0 .gt. 1d-6) then
                            dtnoise = (1.d0 - exp(-2.d0*delta0*dtr))/(2.d0*delta0)
                        else
                            dtnoise = dtr - dtr**2*delta0
                        end if
                    end if
                end if

                if (addrognoso) then

                    if (sum(abs(cov_old(1:ieskin2))) .ne. 0.d0) then
                        ! At the first iteration cov_old is zero
                        call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin, velion, 3&
                                &, 0.d0, velion(2, 1), 3)
                        call dgemv('N', ieskin, ieskin, -1.d0, cov_old, ieskin, velion, 3&
                                &, 1.d0, velion(2, 1), 3)
                        if (yesquantum) velion(2, 1:ieskin) = delta0*velion(2, 1:ieskin)
                    else
                        velion(2, 1:ieskin) = 0.d0
                    end if
                    cov_old(1:ieskin2) = cov(1:ieskin2)
#ifdef PARALLEL
                    if (yesquantum) call bcast_real(cov_old, ieskin2, 0, mpi_comm_world)
#endif

                end if

                call dsyev_my('V', 'L', ieskin, cov, ieskin, eig&
                        &, info, nprocrep, rankrep, commrep_mpi)
                if (rank .eq. 0) then
                    write (6, *) ' Eigenvalues covariance ion forces '
                    condnum = 0.d0
                    do i = 1, ieskin
                        write (6, *) i, eig(i)
                        if (condnum .eq. 0.d0 .and. i .ge. mineig) condnum = eig(i)
                    end do
                    condnum = eig(ieskin)/condnum
                    write (6, *) ' Condition number covariance matrix =', condnum
                end if

                if (yesquantum) then
#ifdef PARALLEL
!      To avoid roundoff all processors have the same cov,psip
                    call bcast_real(cov, ieskin2, 0, mpi_comm_world)
                    call bcast_real(eig, ieskin, 0, mpi_comm_world)
#endif
                end if
                if (info .ne. 0) write (6, *) &
                        &' ERROR in diagonalization of covariance matrix ', rank
                call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin, psip(indi + 1), 1, 0.d0, work, 1)
                if (addrognoso) then
                    cost = -0.5d0/dtstep
                    call dgemv('T', ieskin, ieskin, cost, cov, ieskin, velion(2, 1), 3, 1.d0, work, 1)
                end if

                if (yesquantum) then
                    !      Transform also coordinates
                    call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin, sov5, 1, 0.d0, sov5(1, 5), 1)
#ifdef  PARALLEL
! In the quantum case we have also to apply the transformation kdyn that
! diagonalizes the interaction between the beads. Since each bead is in a
! different processor, this trasformation is done in parallel only
! kdyn matrix elements are common to all processors.
                    do i = 0, nbead - 1
                        sov5(1:ieskin, 1) = kdyn(rankcolrep + 1, i + 1)*work(1:ieskin) !  forces
                        sov5(1:ieskin, 2) = kdyn(rankcolrep + 1, i + 1)*sov5(1:ieskin, 5) !  coordinates
                        call reduce_base_real_to(2*ieskin, sov5, sov5(1, 3), commcolrep_mpi, i)
                    end do
                    work(1:ieskin) = sov5(1:ieskin, 3) - kdyn_eig(rankcolrep + 1)*sov5(1:ieskin, 4) ! the total force
!       sov5(1:ieskin,5)=sov5(1:ieskin,4)
#endif
                    Tmes = 0.d0

                    do i = 1, ieskin
                        eigS = friction + delta0*eig(i) + kdyn_eig(rankcolrep + 1)
                        !     Only finite frequency modes are bounded and virial applies
                        if (rankcolrep .ne. 0) Tmes = Tmes - work(i)*sov5(i, 4)
                        if (temp .eq. 0.d0) then
                            work(i) = dtr*work(i)/eigS
                            if (cleanrognoso) psip(i + indi) = 0.d0
                        else

                            if (abs(eigS) .gt. 1d-8) then
                                work(i) = work(i)/eigS*dtstep
                                if (normcorr .gt. 0.d0) then
                                    cost = 2.d0*temp*(dtnoise/eigS - normcorr*dtnoise**2*eig(i)/eigS**2)
                                    if (cost .lt. 0.d0) then
                                        cost = 0.d0 ! The Umrigar choice
                                        if (rankrep .eq. 0) write (6, *) 'Warning noise correction not possible!!! '
                                    end if
                                else
                                    cost = 2.d0*temp/eigS*dtnoise
                                end if
                                psip(i + indi) = dsqrt(-2.d0*dlog(1.d0 - drand1())*cost)*dcos(2.d0*pi*drand1())
                                work(i) = work(i) + psip(i + indi)

                            else
                                work(i) = 0.d0
                                psip(indi + i) = 0.d0
                            end if
                        end if
                    end do
                    !       Now go back in the original basis
#ifdef  PARALLEL
                    if (cleanrognoso) then
                    do i = 0, nbead - 1
                        sov5(1:ieskin, 1) = kdyn(i + 1, rankcolrep + 1)*work(1:ieskin) !coordinate dstep
                        sov5(1:ieskin, 2) = kdyn(i + 1, rankcolrep + 1)*psip(indi + 1:indi + ieskin) !coordinate dstep only noise
                        call reduce_base_real_to(2*ieskin, sov5, sov5(1, 3), commcolrep_mpi, i)
                    end do
                    call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin&
                &, sov5(1, 3), 1, 0.d0, work, 1)
                    call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin&
                &, sov5(1, 4), 1, 0.d0, psip(indi + 1), 1)
                    else
                    do i = 0, nbead - 1
                        sov5(1:ieskin, 1) = kdyn(i + 1, rankcolrep + 1)*work(1:ieskin) !coordinate dstep
                        call reduce_base_real_to(ieskin, sov5, sov5(1, 2), commcolrep_mpi, i)
                    end do
                    call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin&
                &, sov5(1, 2), 1, 0.d0, work, 1)
                    end if
!       Tmes=dnrm2(ieskin,velion,3)**2/ieskin*ris(2)
                    call reduce_base_real(1, Tmes, commcolrep_mpi, -1)
                    Tmes = Tmes/dble(ieskin)/dble(nbead)/dble(nbead - 1)*ris(2)
                    if (rank .eq. 0) write (6, *) ' Temperature Virial', Tmes
#endif
                    do i = 1, ieskin
                        psip(n2 + indi + i - 1) = work(i)*psip(n3 + i - 1)
                    end do
                    if (cleanrognoso) then
                        !       Only the noisy change due to T finite
                        do i = 1, ieskin
                            psip(indi + i) = psip(i + indi)*psip(n3 + i - 1)
                        end do
                    end if
                else ! following Classical case
                    maxsn = 0.d0
                    !   Hessian is assumed to be delta0 x S
                    if (epstion .gt. 0) then
                        do i = mineig, ieskin
                            maxsn = maxsn + work(i)**2/eig(i)
                        end do
                    else
                        do i = mineig, ieskin
                            if (eig(i) + epstion*trace .gt. 0.d0)&
                                    &maxsn = maxsn + work(i)**2/(eig(i) + epstion*trace)
                        end do
                    end if
                    do i = 1, ieskin
                        if (abs(eig(i)) .gt. 1d-8) then
                            work(i) = work(i)/eig(i)
                            work(i) = work(i)*dtstep
                        else
                            work(i) = 0.d0
                        end if
                    end do
                    if (maxsn .gt. 0.d0) maxsn = dsqrt(maxsn)
                    if (rank .eq. 0) write (6, *) ' Maximum devmax signal/noise ratio =', maxsn
                    call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin, work, 1, 0.d0, psip(indi + 1), 1)
                    indin = indi
                    do j = npps, np
                        !            if(ipip(j).eq.2) then
                        indin = indin + 1
                        psip(n2 + indin - 1) = psip(indin)*psip(n3 + j - npps)
                        if (rank .eq. 0) write (6, *) 'Change ion comp =', indin, psip(indin)
                        !            endif
                    end do

                    if (temp .ne. 0.d0) then

                        if (normcorr .gt. 0.d0) then
                            cost = 2.d0*temp*dtnoise - dtnoise**2
                        else
                            cost = 2.d0*temp*dtnoise
                        end if

                        do j = 1, ieskin
                            psip(indi + j) = dsqrt(-2.d0*dlog(1.d0 - drand1())*cost)*dcos(2.d0*pi*drand1())
                        end do
                        call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin, psip(indi + 1), 1, 0.d0, work, 1)
                        do i = 1, ieskin
                            if (abs(eig(i)) .gt. 1d-8) then
                                work(i) = work(i)/dsqrt(eig(i))
                            else
                                work(i) = 0.d0
                            end if
                        end do
                        call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin, work, 1, 0.d0, psip(indi + 1), 1)
                        indin = indi
                        !  Adding the random noise correlated according to the covariance matrix
                        do j = npps, np
                            indin = indin + 1
                            psip(n2 + indin - 1) = psip(n2 + indin - 1) + psip(indin)*psip(n3 - npps + j)
                            psip(indin) = psip(indin)*psip(n3 - npps + j)
                        end do
                    elseif (cleanrognoso) then
                        do j = 1, ieskin
                            psip(indi + j) = 0.d0
                        end do
                    end if
                end if ! endif yesquantum
                cost = dnrm2(ieskin, psip(n2 + indi), 1)
                if (cost .gt. eps_dyn5 .and. eps_dyn5 .gt. 0.d0) then
                    costn = eps_dyn5/cost
                    if (rank .eq. 0) write (6, *) ' Warning decelerated ions by  ', costn
                    call dscal(ieskin, costn, psip(n2 + indi), 1)
                    if (cleanrognoso) call dscal(ieskin, costn, psip(indi + 1), 1)
                    cost = eps_dyn5
                end if

                !         save coordinate change and store in velion
                if (cleanrognoso) then
                    !  Only the temperature dependent part vanishing for T-->0
                    do j = 1, ieskin
                        velion(1, j) = psip(indi + j)
                    end do
                else
                    do j = 1, ieskin
                        velion(1, j) = psip(n2 + indi + j - 1)
                    end do
                end if
                if (rank .eq. 0) write (6, *) ' Ratio dyn =', cost
                if (rank .eq. 0) write (6, *) ' Norm change ions =', cost

                deallocate (work, eig)
                !+++++ end dynamics 5

            elseif (idyn .eq. 4) then

                if (temp .ne. 0.d0) then
                    cost = delta0*0.5d0/temp
                else
                    cost = 0.d0
                end if

                call dscal(ieskin2, cost, cov, 1)

                indin = 0
                do j = npps, np
                    indin = indin + 1
                    psip(n3 + indin - 1) = dsqrt(scalpar(ipip(n3 + j - 1))/dt)
                    !       write(6,*) ' Inverse Mass =',indin,ipip(n3+j-1),psip(n3+indin-1)**2
                end do

                if (indin .ne. ieskin) then
                    if (rank .eq. 0) write (6, *) ' Error in dynamic ieskin ne indin !!!  ', ieskin, indin
                    errnoise = 6
                end if
                !       stop
                !
                do j = 1, ieskin
                    do i = 1, ieskin
                        cov(ieskin*(j - 1) + i) = cov(ieskin*(j - 1) + i)*psip(n3 + i - 1)*psip(n3 + j - 1)
                    end do
                end do

                !        save psip the QMC forces
                call dcopy(np, psip, 1, sov4(1, 2), 1)
                !        scale the force by the Mass
                do i = 1, ieskin
                    sov4(indi + i, 2) = sov4(indi + i, 2)*psip(n3 + i - 1)
                end do

                !         write(6,*) ' velion force before '
                !         do i=1,ieskin
                !         write(6,*) i,velion(3,i),sov(npp,indi+i,4)
                !         enddo

                do i = 1, ieskin
                    cov(ieskin*(i - 1) + i) = cov(ieskin*(i - 1) + i) + friction
                end do

                lworkd = 3*ieskin
                call dsyev_my('V', 'L', ieskin, cov, ieskin, psip, info, nprocrep, rankrep, commrep_mpi)
                !     call dsyev('V','L',ieskin,cov,ieskin,psip,psip(n5),lworkd,info)
                if (rank .eq. 0) then
                    write (6, *) ' Eigenvalues covariance '
                    do i = 1, ieskin
                        write (6, *) i, psip(i)
                    end do
                end if
#ifdef DEBUG
!  Z2   Gauge fixing
                do i = 1, ieskin
                    cost = abs(cov(ieskin*(i - 1) + 1))
                    indm = 1
                    do j = 2, ieskin
                    if (abs(cov(ieskin*(i - 1) + j)) .gt. cost) then
                        cost = abs(cov(ieskin*(i - 1) + j))
                        indm = j
                    end if
                    end do
                    if (cov(indm + (i - 1)*ieskin) .lt. 0.d0) cov(1 + ieskin*(i - 1):ieskin*i) = -cov(1 + ieskin*(i - 1):ieskin*i)
                end do
#endif

                if (info .ne. 0) then
                    if (rank .eq. 0) write (6, *) ' Error in lapack dsyev  dynamic !!! '
                    errnoise = 4
                end if

                if (rank .eq. 0) write (6, *) ' Ratio dyn =', (psip(ieskin) - friction)/sqrt(temp*157.8873306d0)

                !         write(6,*) ' Spectrum gamma=',info,(psip(i)*dt,i=1,ieskin)

                !       Changing basis, using the one that diagonalizes cov

                call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                        &, velion(3, 1), 3, 0.d0, sov4, 1)
                call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                        &, sov4(indi + 1, 2), 1, 0.d0, velion(2, 1), 3)
                !     Here sov4(*,1) contains the tranformed velocities
                !     velion(2,*) the transformed forces

#ifdef DEBUG
                if (rank .eq. 0) write (6, *) ' changed velocities ',&
                  &sum(velion(2, 4:ieskin)), sum(sov4(4:ieskin, 1))
                if (rank .eq. 0) write (6, *) ' changed velocities squares ',&
                &sum(velion(2, 4:ieskin)**2), sum(sov4(4:ieskin, 1)**2)
#endif
                alphaqmc = 0.d0
                !!!! INIZIO NEW DYN
#ifdef DEBUG
                do i = 4, ieskin
#else
                    do i = 1, ieskin
#endif
                        zetan(1) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                        zetan(2) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())

                        if (psip(i) .gt. 0.d0) then
                            costh = exp(-dt*psip(i)/2.d0)
                            cost = exp(-dt*psip(i))
                            alphaall = 2.d0*temp*psip(i)
                            if (delta0 .ne. 0.d0 .and. normcorr .gt. 0.d0) then
                                alphaqmc = 2.d0*normcorr*temp*(psip(i) - friction)/delta0
                            else
                                alphaqmc = 0.d0
                            end if
                            if (cost .lt. 0.99999d0) then
                                costm1 = cost - 1.d0
                                costhm1 = costh - 1.d0
                            else
                                !         for numerical stability
                                costm1 = -dt*psip(i) + (dt*psip(i))**2/2.d0
                                costhm1 = -dt*psip(i)/2 + (dt*psip(i))**2/8.d0
                            end if
                            mnoise(1, 1) = -psip(i)*temp*psip(i)* &
                                    &  (1.d0 + cost)/costm1 - alphaqmc
                            Gn = -costm1/psip(i)

                            Tn = 1.d0/psip(i)*(dt + costh/psip(i)*costm1)

                            mnoise(2, 2) = alphaall/Tn**2* &
                                    &(dt/psip(i)**2 + 0.5d0/psip(i)**3*(2.d0 - cost)*costm1*(1.d0 + cost))&
                                    & - alphaqmc

                            mnoise(1, 2) = alphaall/(Tn*2.d0*psip(i)*(1.d0 + costh))* &
                                    & (2.d0 - cost*costh*costm1/costhm1) - alphaqmc

                            mnoise(2, 1) = mnoise(1, 2)

                            !      make a square root of a positive definite matrix to define the noise.
                            all_dyn = all_dyn + mnoise(1, 1) + mnoise(2, 2)
                            ratio_dyn = ratio_dyn + 2*alphaqmc

                            call root2mat(mnoise, errnoise)

                        else
                            if (rank .eq. 0) write (6, *) ' ERROR negative eigenv. ', psip(i) &
                                    &, alphaqmc
                            errnoise = 1
                        end if !if psip(i).gt.0

                        if (rank .eq. 0 .and. errnoise .ne. 0) write (6, *) ' Error negative noise matrix !!!'

                        eta_v = mnoise(1, 1)*zetan(1) + mnoise(1, 2)*zetan(2)
                        eta_r = mnoise(2, 1)*zetan(1) + mnoise(2, 2)*zetan(2)

                        !        Now the consistent change of coordinates put in velion(1,*)
                        velion(1, i) = sov4(i, 1)*costh*Gn + Tn*(velion(2, i) + eta_r)

                        velion(2, i) = sov4(i, 1)*cost + Gn*(velion(2, i) + eta_v)

                    end do
#ifdef DEBUG
                    velion(1:2, 1:3) = 0.d0
                    if (rank .eq. 0) write (6, *) ' changed velocities 2 ',&
             &sum(velion(2, 1:ieskin)), sum(velion(1, 1:ieskin))
                    if (rank .eq. 0) write (6, *) ' changed velocities squares 2  ',&
                    &sum(velion(2, 1:ieskin)**2), sum(velion(1, 1:ieskin)**2)
#endif

                    !         Now go back to the original basis
                    call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin               &
                            &, velion(2, 1), 3, 0.d0, velion(3, 1), 3)
                    call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin               &
                            &, velion(1, 1), 3, 0.d0, psip(n2 + indi), 1)

                    !       back to the scale of coordinates
                    !       write(6,*) '  Final correction  inside '
                    cost = 0.d0
                    do i = 1, ieskin
                        psip(n2 + indi + i - 1) = psip(n2 + indi + i - 1)*psip(n3 + i - 1)
                        if (psip(n3 + i - 1) .ne. 0.d0) then
                            cost = cost + 1.d0
                        else
                            velion(3, i) = 0.d0
                        end if
                        !       write(6,*) i,psip(n2+indi+i-1),psip(n3+i-1),velion(3,i)
                    end do
                    Tmes = dnrm2(ieskin, velion(3, 1), 3)**2/cost*ris(2)
                    if (rank .eq. 0) write (6, *) ' Temperature (H)= ', Tmes
                    cost = dnrm2(ieskin, psip(n2 + indi), 1)
                    if (rank .eq. 0) write (6, *) ' Norm change ions =', cost

                !!!! FINE  NEW DYN
                    elseif (idyn .eq. 7) then

                    !   In this routine we apply the Trotter in the Master equation according to
                    !   Ceriotti's paper. Thus we first apply the idyn=3 scheme by changing only
                    !   the velocities and without moving the ions, with the covariance matrix
                    !   defined for each bead independently.
                    !   Then we apply the exact damped dynamic for the elastic quantum case
                    !   without any force.
                    if (yesrootc) then
                        cost = (delta0/2.d0)**2
                    else
                        if (temp .ne. 0.d0) then
                            cost = delta0*0.5d0/temp
                        else
                            cost = 0.d0
                        end if
                    end if

                    call dscal(ieskin2, cost, cov, 1)

                    if (yessecond) then
                        dth = dt/2.d0
                    else
                        dth = dt
                    end if

                    indin = 0
                    do j = npps, np
                        indin = indin + 1
                        psip(n3 + indin - 1) = dsqrt(scalpar(ipip(n3 + j - 1))/dt)
                    end do

                    if (indin .ne. ieskin) then
                        if (rank .eq. 0) write (6, *) ' Error in dynamic ieskin ne indin !!!  ', ieskin, indin
                        errnoise = 6
                    end if

                    do j = 1, ieskin
                        do i = 1, ieskin
                            cov(ieskin*(j - 1) + i) = cov(ieskin*(j - 1) + i)*psip(n3 + i - 1)*psip(n3 + j - 1)
                        end do
                    end do

                    !        first diagonalize sov(1,1,4)
                    !        save psip the QMC forces
                    call dcopy(np, psip, 1, sov4(1, 2), 1)
                    !        scale the force by the Mass
                    do i = 1, ieskin
                        sov4(indi + i, 2) = sov4(indi + i, 2)*psip(n3 + i - 1)
                    end do

                    !        do i=1,ieskin
                    !        cov(ieskin*(i-1)+i)=cov(ieskin*(i-1)+i)+friction
                    !        enddo

                    lworkd = 3*ieskin

                    !     if((yesquantum.and.nrep_bead.gt.1)) then
                    call dsyev_my('V', 'L', ieskin, cov, ieskin, psip, info, nprocrep, rankrep, commrep_mpi)
                    !     else
                    !     call dsyev('V','L',ieskin,cov,ieskin,psip,psip(n5),lworkd,info)
                    !     endif
                    if (yesrootc) then
                        psip(i) = dsqrt(max(psip(i), 0.d0)) + friction
                    else
                        do i = 1, ieskin
                            psip(i) = psip(i) + friction
                        end do
                    end if

                    if (rank .eq. 0) then
                        write (6, *) ' Eigenvalues covariance '
                        do i = 1, ieskin
                            write (6, *) i, psip(i)
                        end do
                    end if
#ifdef DEBUG
!  Z2   Gauge fixing
                    do i = 1, ieskin
                        cost = abs(cov(ieskin*(i - 1) + 1))
                        indm = 1
                        do j = 2, ieskin
                        if (abs(cov(ieskin*(i - 1) + j)) .gt. cost) then
                            cost = abs(cov(ieskin*(i - 1) + j))
                            indm = j
                        end if
                        end do
                      if (cov(indm + (i - 1)*ieskin) .lt. 0.d0) cov(1 + ieskin*(i - 1):ieskin*i)&
                         & = -cov(1 + ieskin*(i - 1):ieskin*i)
                    end do
#endif

                    if (info .ne. 0) then
                        if (rank .eq. 0) write (6, *) ' Error in lapack dsyev  dynamic !!! '
                        errnoise = 4
                    end if

                    if (rank .eq. 0) write (6, *) ' Ratio dyn =', (psip(ieskin) - friction)/sqrt(temp*157.8873306d0)

                    do i = 1, ieskin
                        if (psip(i) .gt. 0.d0) then
                            psip(n5 + i - 1) = (1.d0 - dexp(-psip(i)*dth))/psip(i)
                        else
                            if (rank .eq. 0) write (6, *) ' Refused eigenvalue ', i
                            psip(n5 + i - 1) = dth
                        end if
                    end do
                    !         compute the noise correction
                    alphaqmc = 0.d0

                    do i = 1, ieskin
                        if (psip(i) .gt. 0.d0) then
                            !         psip(n4+i-1)=(psip(i)*2.d0*temp)*psip(i)*
                            !    1  sinh(dt*psip(i))/(4.d0*sinh(psip(i)*dt*.5d0)**2)
                            !         the following is protected for overflow
                            cost = exp(-dth*psip(i))
                            if (delta0 .ne. 0.d0 .and. normcorr .gt. 0.d0) then
                                if (yesrootc) then
                                    alphaqmc = normcorr*(2.d0*(psip(i) - friction)/delta0)**2
                                else
                                    alphaqmc = normcorr*2.d0*temp*(psip(i) - friction)/delta0
                                end if
                            else
                                alphaqmc = 0.d0
                            end if
                            psip(n4 + i - 1) = psip(i)*temp*psip(i)* &
                                    &  (1.d0 + cost)/(1.d0 - cost) - alphaqmc
                            !          subtracting the noise already present in the forces
                        else
                            psip(n4 + i - 1) = 0.d0
                        end if
                        all_dyn = all_dyn + psip(n4 + i - 1) + alphaqmc
                        ratio_dyn = ratio_dyn + alphaqmc

                        if (psip(n4 + i - 1) .gt. 0.d0) then
                            psip(n4 + i - 1) = dsqrt(psip(n4 + i - 1))
                        else

                            psip(n4 + i - 1) = 0.d0
                            if (yesrootc) write (6, *) ' Warning noise correction not possible ', i, psip(i), psip(n4 + i - 1)
                            if (friction .gt. 0.d0 .and. .not. yesrootc) then
                                if (rank .eq. 0) write (6, *) ' There should be some error in reweight0         &
                                        &  ', i, psip(n4 + i - 1)
                                errnoise = 5
                            end if
                        end if
                    end do

                    !         now we are able to update the velocity
                    !         first change basis actual velocity and force
                    call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                            &, velion(3, 1), 3, 0.d0, sov4, 1)
                    call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                            &, sov4(indi + 1, 2), 1, 0.d0, velion(2, 1), 3)

                    if (yessecond) then
                        ! Compute the temperature at half time intervals.
                        do i = 1, ieskin
                            zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                            sov4(i, 1) = sov4(i, 1)*dexp(-psip(i)*dth)                    &
                                    & + psip(n5 + i - 1)*(velion(2, i) + psip(n4 + i - 1)*zeta)
                        end do
                        Tmes = dnrm2(ieskin, sov4, 1)**2/ieskin*ris(2)
                        ! Second  half time interval.
                        do i = 1, ieskin
                            zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                            velion(2, i) = sov4(i, 1)*dexp(-psip(i)*dth)                    &
                                    & + psip(n5 + i - 1)*(velion(2, i) + psip(n4 + i - 1)*zeta)
                        end do
                    else
                        do i = 1, ieskin
                            zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                            velion(2, i) = sov4(i, 1)*dexp(-psip(i)*dth)                    &
                                    & + psip(n5 + i - 1)*(velion(2, i) + psip(n4 + i - 1)*zeta)
                        end do
                    end if

                    !         go back in the original basis

                    call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin              &
                            &, velion(2, 1), 3, 0.d0, velion(3, 1), 3)

                    if (.not. yessecond) Tmes = dnrm2(ieskin, velion(3, 1), 3)**2/ieskin*ris(2)
                    if (rank .eq. 0) write (6, *) ' Temperature (H)= ', Tmes

                    !          call dgemv('T',ieskin,ieskin,1.d0,cov,ieskin              &
                    !     &,velion(3,1),3,0.d0,sov4,1)
                    !          call dgemv('T',ieskin,ieskin,1.d0,cov,ieskin              &
                    !     &,sov4(indi+1,2),1,0.d0,velion(2,1),3)

                    !      cov is the identity for idyn==7

                    !      call dcopy(ieskin,sov4(indi+1,2),1,velion(2,1),3)
                    velion(2, 1:ieskin) = 0.d0 ! no forces here
                    call dcopy(ieskin, velion(3, 1), 3, sov4, 1)
                    !      Here sov4(*,1) contains the transfomed velocities
                    !           velion(2,*) the transformed forces
                    !    In the following velion(1,*) is unchanged

                    !    Transformation coordinates in the basis of eigenvectors
                    !    velion(1,*) represents run time the difference between the atomic position
                    !    when the forces are computed and the coordinates when we start to apply
                    !    the new LD integrator. velion(1,*) is non zero only for second order
                    !    dynamics (yessecond), more accurate with the time step integration.
                    !    Notice that, after that,  we have to scale by the sqrt of the mass (psip(n3)~1/sqrt(m))
                    do j = 1, nion
                        do i = 1, 3
                            ind = i + (j - 1)*3
                            sov5(ind, 4) = rion(i, j)/psip(n3 + ind - 1)
                        end do
                    end do
#ifdef  PARALLEL
! In the quantum case we have also to apply the transformation kdyn that
! diagonalizes the interaction between the beads. Since each bead is in a
! different processor, this trasformation is done in parallel only
! kdyn matrix elements are common to all processors.
                    do i = 0, nbead - 1
                        sov5(1:ieskin, 1) = kdyn(rankcolrep + 1, i + 1)*sov4(1:ieskin, 1) !  velocities
!       sov5(1:ieskin,2)=kdyn(rankcolrep+1,i+1)*velion(2,1:ieskin)    !  forces
                        sov5(1:ieskin, 2) = kdyn(rankcolrep + 1, i + 1)*sov5(1:ieskin, 4) !  coordinates
                        call reduce_base_real_to(2*ieskin, sov5, sov5(1, 6), commcolrep_mpi, i)
                    end do
                    sov4(1:ieskin, 1) = sov5(1:ieskin, 6)
!       velion(2,*)  is actually zero below
#endif
#ifdef DEBUG
                    if (rank .eq. 0) write (6, *) ' changed velocities ',&
                      &sum(velion(2, 4:ieskin)), sum(sov4(4:ieskin, 1))
                    if (rank .eq. 0) write (6, *) ' changed velocities squares ',&
                    &sum(velion(2, 4:ieskin)**2), sum(sov4(4:ieskin, 1)**2)
#endif
                    alphaqmc = 0.d0
                !!!! INIZIO NEW DYN : Now everything is in the diagonal basis and the integration
                    !can be done exactly for each of the ieskin (x nbead , index labelled by
                    ! rankcolrep+1 in the quantum case) eigenvalues

#ifdef DEBUG
                    do i = 4, ieskin
#else
                        do i = 1, ieskin
#endif
                            zetan(1) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                            zetan(2) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())

                            alphaqmc = 0.d0
                            alphaall = 0.d0

                            ! gammaall
                            psip(i) = 2.d0*delta0k*sqrt(abs(kdyn_eig(rankcolrep + 1)))
                            if (psip(i) .lt. delta0q) psip(i) = delta0q ! cutoff on the lowest eigenvalue. Ceriotti's choice.

                            alphaall = 2*temp*psip(i)

                            call set_turboq(psip(i), kdyn_eig(rankcolrep + 1), dt&
                                    &, alphaall, alphaqmc, mnoise, Gn, Tn, Gni, Gnh, Gnih)

                            Gnh = 0.d0
                            Gnih = 0.d0

                            call root2mat(mnoise, errnoise)

                            if (errnoise .ne. 0 .and. rank .eq. 0) write (6, *) ' Error negative definite matrix noise '

                            eta_v = mnoise(1, 1)*zetan(1) + mnoise(1, 2)*zetan(2)
                            eta_r = mnoise(2, 1)*zetan(1) + mnoise(2, 2)*zetan(2)

                            !        Now the consistent change of coordinates put in velion(1,*)
                            !       velion(3,*) represents an extra displacement of the positions by dt/2 without the noise
                            !       and the force (Gn=0), i.e. consistent up to order dt^3/2.

                            sov4(i, 2) = sov4(i, 1)*Gni(2, 1) + Gni(2, 2)*sov5(i, 7) + Tn*eta_r
                            velion(2, i) = sov4(i, 1)*Gni(1, 1) + Gni(1, 2)*sov5(i, 7) + Gn*eta_v

                        end do
#ifdef DEBUG
                        velion(1:2, 1:3) = 0.d0
                        if (rank .eq. 0) write (6, *) ' changed velocities 2 ',&
                 &sum(velion(2, 1:ieskin)), sum(velion(1, 1:ieskin))
                        if (rank .eq. 0) write (6, *) ' changed velocities squares 2  ',&
                        &sum(velion(2, 1:ieskin)**2), sum(velion(1, 1:ieskin)**2)
#endif
                        !         Now go back to the original basis
                        !       Tmes=dnrm2(ieskin,velion,3)**2/ieskin*ris(2)
                        !       call reduce_base_real(1,Tmes,commcolrep_mpi,-1)
                        !       if(rank.eq.0)  write(6,*) ' sum after0',Tmes/nbead
                        !    Transformation back in the original basis
#ifdef  PARALLEL
                        do i = 0, nbead - 1
                            sov5(1:ieskin, 1) = kdyn(i + 1, rankcolrep + 1)*sov4(1:ieskin, 2) ! coordinates
                            sov5(1:ieskin, 2) = kdyn(i + 1, rankcolrep + 1)*velion(2, 1:ieskin) !  velocities
!        sov5(1:ieskin,3)=kdyn(i+1,rankcolrep+1)*velion(3,1:ieskin)   ! extra prop coordinates
                            call reduce_base_real_to(2*ieskin, sov5, sov5(1, 4), commcolrep_mpi, i)
                        end do
                        sov4(1:ieskin, 2) = sov5(1:ieskin, 4)
                        velion(2, 1:ieskin) = sov5(1:ieskin, 5)
!       velion(3,1:ieskin)=sov5(1:ieskin,6)
!       Tmes=dnrm2(ieskin,velion,3)**2/ieskin*ris(2)
!       call reduce_base_real(1,Tmes,commcolrep_mpi,-1)
!       if(rank.eq.0)  write(6,*) ' sum after1',Tmes/nbead
#endif

                        !        call dgemv('N',ieskin,ieskin,1.d0,cov,ieskin               &
                        !    &,velion(3,1),3,0.d0,sov4,1)
                        !        call dcopy(ieskin,velion(3,1),3,sov4,1)
                        !        call dgemv('N',ieskin,ieskin,1.d0,cov,ieskin               &
                        !    &,velion(2,1),3,0.d0,velion(3,1),3)
                        call dcopy(ieskin, velion(2, 1), 3, velion(3, 1), 3)
                        !        call dgemv('N',ieskin,ieskin,1.d0,cov,ieskin               &
                        !    &,sov4(1,2),1,0.d0,psip(n2+indi),1)
                        call dcopy(ieskin, sov4(1, 2), 1, psip(n2 + indi), 1)

                        !       back to the scale of coordinates
                        !       write(6,*) '  Final correction  inside '
                        do i = 1, ieskin
                            psip(n2 + indi + i - 1) = psip(n2 + indi + i - 1)*psip(n3 + i - 1)
                            if (psip(n3 + i - 1) .eq. 0.d0) then
                                velion(3, i) = 0.d0
                            end if
                        end do
                        cost = dnrm2(ieskin, psip(n2 + indi), 1)
                        if (rank .eq. 0) write (6, *) ' Norm change ions =', cost

                !!!! FINE  NEW IDYN=7  !!!!!

                        elseif (idyn .eq. 8) then !!! Bare Ceriotti algorithm

                        if (temp .ne. 0.d0) then
                            cost = delta0*0.5d0/temp
                        else
                            cost = 0.d0
                        end if

                        dth = dt/2.d0

                        call dscal(ieskin2, cost, cov, 1)

                        indin = 0
                        do j = npps, np
                            indin = indin + 1
                            ! scaling as 1/sqrt(m)
                            psip(n3 + indin - 1) = dsqrt(scalpar(ipip(n3 + j - 1))/dt)
                        end do

                        ! rescale the covariance matrix by masses
                        do j = 1, ieskin
                            do i = 1, ieskin
                                cov(ieskin*(j - 1) + i) = cov(ieskin*(j - 1) + i)*psip(n3 + i - 1)*psip(n3 + j - 1)
                            end do
                        end do

                        ! scale the force by the Mass
                        call dcopy(np, psip, 1, sov4(1, 2), 1) ! forces stored in sov4(1:ieskin,2)
                        do i = 1, ieskin
                            sov4(indi + i, 2) = sov4(indi + i, 2)*psip(n3 + i - 1)
                        end do

                        ! add diagonal term to covariance
                        do i = 1, ieskin
                            cov(ieskin*(i - 1) + i) = cov(ieskin*(i - 1) + i) + friction
                        end do

                        ! Save covariance matrix before it is destroyed by diagonalization
                        cov_sav(1:ieskin*ieskin) = cov(1:ieskin*ieskin)

                        ! diagonalize gamma matrix in the physical bead space stored in cov
                        !       lworkd=3*ieskin

                        ! previous version
                        !       call dsyev_my('V','L',ieskin,cov,ieskin,psip,psip(n5) &
                        !            ,lworkd,info,nprocrep,rankrep,commrep_mpi)

                        ! WARNING: mind that in the latest version dsyev_my has changed!!!!!!
                        ! Use the following one instead
                        call dsyev_my('V', 'L', ieskin, cov, ieskin, psip, info, nprocrep, rankrep, commrep_mpi)

                        if (rank .eq. 0) then
                            write (6, *) ' Eigenvalues covariance '
                            do i = 1, ieskin
                                write (6, *) i, psip(i)
                            end do
                        end if
#ifdef DEBUG
                        !  Z2   Gauge fixing
                        do i = 1, ieskin
                            cost = abs(cov(ieskin*(i - 1) + 1))
                            indm = 1
                            do j = 2, ieskin
                                if (abs(cov(ieskin*(i - 1) + j)) .gt. cost) then
                                    cost = abs(cov(ieskin*(i - 1) + j))
                                    indm = j
                                end if
                            end do
                      if (cov(indm + (i - 1)*ieskin) .lt. 0.d0) cov(1 + ieskin*(i - 1):ieskin*i)&
                         & = -cov(1 + ieskin*(i - 1):ieskin*i)
                        end do
#endif

                        if (info .ne. 0) then
                            if (rank .eq. 0) write (6, *) ' Error in lapack dsyev  dynamic !!! ', info
                            errnoise = 4
                        end if

                        ! The scale invariant criterium is gamma/sqrt(T) , so I have divided by the
                        ! sqrt(T/T_0) where T_0=1000K
                        if (rank .eq. 0) write (6, *) ' Ratio dyn =', (psip(ieskin) - friction)/sqrt(temp*157.8873306d0)

                        ! now compute coefficients as function of gamma eigenvalues
                        do i = 1, ieskin
                            if (psip(i) .gt. 0.d0) then
                                psip(n5 + i - 1) = (1.d0 - dexp(-psip(i)*dth))/psip(i)
                            else
                                if (rank .eq. 0) write (6, *) ' Refused eigenvalue ', i
                                psip(n5 + i - 1) = dth
                            end if
                        end do

                        ! compute the noise correction
                        alphaqmc = 0.d0

                        do i = 1, ieskin
                            if (psip(i) .gt. 0.d0) then
                                cost = exp(-dth*psip(i))
                                if (delta0 .ne. 0.d0 .and. normcorr .gt. 0.d0) then
                                    alphaqmc = normcorr*2.d0*temp*(psip(i) - friction)/delta0
                                else
                                    alphaqmc = 0.d0
                                end if
                                psip(n4 + i - 1) = temp*psip(i)**2*(1.d0 + cost)/(1.d0 - cost) - alphaqmc
                                !          subtracting the noise already present in the forces
                            else
                                psip(n4 + i - 1) = 0.d0
                            end if

                            if (psip(n4 + i - 1) .gt. 0.d0) then
                                psip(n4 + i - 1) = dsqrt(psip(n4 + i - 1))
                            else

                                psip(n4 + i - 1) = 0.d0
                                if (friction .gt. 0.d0) then
                                    if (rank .eq. 0) write (6, *) ' There should be some error in reweight0', i, psip(n4 + i - 1)
                                    errnoise = 5
                                end if
                            end if
                        end do

                        ! now we are able to update the velocity
                        ! first change velocity and force into space which diagonalizes gamma
                        call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin &
                                   , velion(3, 1), 3, 0.d0, sov4, 1)
                        call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin &
                                   , sov4(indi + 1, 2), 1, 0.d0, velion(2, 1), 3)

                        ! first half iteration for Born-Oppenheimer damping
                        do i = 1, ieskin
                            zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                            sov4(i, 1) = sov4(i, 1)*dexp(-psip(i)*dth) &
                                         + psip(n5 + i - 1)*(velion(2, i) + psip(n4 + i - 1)*zeta)
                        end do

                        ! go back in the original (physical) basis

                        call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin &
                                   , sov4, 1, 0.d0, velion(3, 1), 3)

                        ! Here sov4(*,1) contains the transformed velocities
                        call dcopy(ieskin, velion(3, 1), 3, sov4, 1)

                        if (yesquantum) then
                            ! Quantum harmonic damping in normal modes as suggested by Ceriotti

#ifdef  PARALLEL
! In the quantum case we have also to apply the transformation kdyn that
! diagonalizes the interaction between the beads. Since each bead is in a
! different processor, this trasformation is done in parallel only
! kdyn matrix elements are common to all processors.
! Transformation coordinates in the basis of eigenvectors

                            do i = 0, nbead - 1
                                sov5(1:ieskin, 1) = kdyn(rankcolrep + 1, i + 1)*sov4(1:ieskin, 1) !  velocities
                                call reduce_base_real_to(ieskin, sov5, sov5(1, 6), commcolrep_mpi, i)
                            end do
! sov5(i,6) ---> velocities

#endif
#ifdef DEBUG
                            if (rank .eq. 0) write (6, *) ' changed velocities ', sum(sov5(4:ieskin, 6))
                            if (rank .eq. 0) write (6, *) ' changed velocities squares ', sum(sov5(4:ieskin, 6)**2)
#endif

#ifdef DEBUG
                            do i = 4, ieskin
#else
                                do i = 1, ieskin
#endif
                                    zetan(1) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())

                                    ! gammaall
                                    psip(i) = 2.d0*delta0k*sqrt(abs(kdyn_eig(rankcolrep + 1)))
                                    ! cutoff on the lowest eigenvalue. Ceriotti's choice.
                                    if (psip(i) .lt. delta0q) psip(i) = delta0q

                                    ! perform the Langevin step in the normal
                                    ! modes as before, except that here there
                                    ! is no intrisinc noise
                                    ! namely alphaqmc=0
                                    ! First half time iteration

                                    cost = exp(-dth*psip(i))
                                    psip(n4 + i - 1) = temp*psip(i)**2*(1.d0 + cost)/(1.d0 - cost)
                                    psip(n4 + i - 1) = sqrt(psip(n4 + i - 1))
                                    psip(n5 + i - 1) = (1.d0 - dexp(-psip(i)*dth))/psip(i)

                                    velion(2, i) = sov5(i, 6)*dexp(-psip(i)*dth) + psip(n5 + i - 1)*psip(n4 + i - 1)*zetan(1)

                                end do

                                ! Temperature estimation
                                Tmes = dnrm2(ieskin, velion(2, 1), 3)**2/ieskin*ris(2)
                                if (rank .eq. 0) write (6, *) ' Temperature (H)= ', Tmes

                                ! perform the Langevin step in the normal modes
                                ! as before, except that here there is no
                                ! intrisinc noise
                                ! namely alphaqmc=0
                                ! Second half time iteration

#ifdef DEBUG
                                do i = 4, ieskin
#else
                                    do i = 1, ieskin
#endif
                                        zetan(1) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                                        velion(2, i) = velion(2, i)&
                                            &*dexp(-psip(i)*dth) + psip(n5 + i - 1)*psip(n4 + i - 1)*zetan(1)
                                    end do

#ifdef DEBUG
                                    velion(1:2, 1:3) = 0.d0
                                    if (rank .eq. 0) write (6, *) ' changed velocities 2 ', &
                                        sum(velion(2, 1:ieskin))
                                    if (rank .eq. 0) write (6, *) ' changed velocities squares 2  ', &
                                        sum(velion(2, 1:ieskin)**2)
#endif
                                    !       Now go back to the original basis
                                    !       Tmes=dnrm2(ieskin,velion,3)**2/ieskin*ris(2)
                                    !       call reduce_base_real(1,Tmes,commcolrep_mpi,-1)
                                    !       if(rank.eq.0)  write(6,*) ' sum after0',Tmes/nbead
                                    !    Transformation back in the original basis from the normal modes of harmonic matrix
#ifdef PARALLEL
                                    do i = 0, nbead - 1
                                        sov5(1:ieskin, 1) = kdyn(i + 1, rankcolrep + 1)*velion(2, 1:ieskin) !  velocities
                                        call reduce_base_real_to(ieskin, sov5, sov5(1, 4), commcolrep_mpi, i)
                                    end do
                                    velion(2, 1:ieskin) = sov5(1:ieskin, 4)
#endif

                                    else ! classical case

#ifdef DEBUG
                                    do i = 4, ieskin
#else
                                        do i = 1, ieskin
#endif
                                            zetan(1) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())

                                            ! gammaall
                                            psip(i) = delta0q ! classical friction

                                            ! perform the classical Langevin
                                            ! step as before, except that here
                                            ! there is no intrisinc noise
                                            ! namely alphaqmc=0
                                            ! First half time iteration

                                            cost = exp(-dth*psip(i))
                                            psip(n4 + i - 1) = temp*psip(i)**2*(1.d0 + cost)/(1.d0 - cost)
                                            psip(n4 + i - 1) = sqrt(psip(n4 + i - 1))
                                            psip(n5 + i - 1) = (1.d0 - dexp(-psip(i)*dth))/psip(i)

                                            velion(2, i) = sov4(i, 1)&
                                                &*dexp(-psip(i)*dth) + psip(n5 + i - 1)*psip(n4 + i - 1)*zetan(1)

                                        end do

                                        ! Temperature estimation
                                        Tmes = dnrm2(ieskin, velion(2, 1), 3)**2/ieskin*ris(2)
                                        if (rank .eq. 0) write (6, *) ' Temperature (H)= ', Tmes

                                        ! perform the Langevin step in the
                                        ! normal modes as before, except that here there is no intrisinc noise
                                        ! namely alphaqmc=0
                                        ! Second half time iteration

#ifdef DEBUG
                                        do i = 4, ieskin
#else
                                            do i = 1, ieskin
#endif
                                                zetan(1) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                                         velion(2, i) = velion(2, i)&
                                             &*dexp(-psip(i)*dth) + psip(n5 + i - 1)*psip(n4 + i - 1)*zetan(1)

                                            end do

#ifdef DEBUG
                                            velion(1:2, 1:3) = 0.d0
                                            if (rank .eq. 0) write (6, *) ' changed velocities 2 ', &
                                                sum(velion(2, 1:ieskin))
                                            if (rank .eq. 0) write (6, *) ' changed velocities squares 2  ', &
                                                sum(velion(2, 1:ieskin)**2)
#endif

                                            end if

                                            ! Old velocities stored in velion(2,1)
                                            velion(3, 1:ieskin) = velion(2, 1:ieskin)

                                            ! Recover the former covariance matrix to apply second BO damping iteration later
                                            cov(1:ieskin*ieskin) = cov_sav(1:ieskin*ieskin)

                                            ! diagonalize gamma matrix (useless step in practice but we need to recover psip)
                                            call dsyev_my('V', 'L', ieskin, cov&
                                                &, ieskin, psip, info, nprocrep&
                                                &, rankrep, commrep_mpi)

                                            if (rank .eq. 0) then
                                                write (6, *) ' Eigenvalues covariance 2'
                                                do i = 1, ieskin
                                                    write (6, *) i, psip(i)
                                                end do
                                            end if

#ifdef DEBUG
                                            !  Z2   Gauge fixing
                                            do i = 1, ieskin
                                                cost = abs(cov(ieskin*(i - 1) + 1))
                                                indm = 1
                                                do j = 2, ieskin
                                                    if (abs(cov(ieskin*(i - 1) + j)) .gt. cost) then
                                                        cost = abs(cov(ieskin*(i - 1) + j))
                                                        indm = j
                                                    end if
                                                end do
                                                if (cov(indm + (i - 1)*ieskin) .lt. 0.d0) cov(1 + ieskin*(i - 1):ieskin*i)&
                                                    & = -cov(1 + ieskin*(i - 1):ieskin*i)
                                            end do
#endif

                                            if (info .ne. 0) then
                                                if (rank .eq. 0) write (6, *) ' Error in lapack dsyev  dynamic !!! ', info
                                                errnoise = 4
                                            end if

                                            ! The scale invariant criterium is gamma/sqrt(T) , so I have divided by the
                                            ! sqrt(T/T_0) where T_0=1000K
                                    if (rank .eq. 0) write (6, *)&
                                       & ' Ratio dyn =', (psip(ieskin) - friction)/sqrt(temp*157.8873306d0)

                                            ! now compute coefficients as function of gamma eigenvalues
                                            do i = 1, ieskin
                                                if (psip(i) .gt. 0.d0) then
                                                    psip(n5 + i - 1) = (1.d0 - dexp(-psip(i)*dth))/psip(i)
                                                else
                                                    if (rank .eq. 0) write (6, *) ' Refused eigenvalue ', i
                                                    psip(n5 + i - 1) = dth
                                                end if
                                            end do

                                            ! compute the noise correction
                                            alphaqmc = 0.d0

                                            do i = 1, ieskin
                                                if (psip(i) .gt. 0.d0) then
                                                    cost = exp(-dth*psip(i))
                                                    if (delta0 .ne. 0.d0 .and. normcorr .gt. 0.d0) then
                                                        alphaqmc = normcorr*2.d0*temp*(psip(i) - friction)/delta0
                                                    else
                                                        alphaqmc = 0.d0
                                                    end if
                                                    psip(n4 + i - 1) = temp*psip(i)**2*(1.d0 + cost)/(1.d0 - cost) - alphaqmc
                                                    !          subtracting the noise already present in the forces
                                                else
                                                    psip(n4 + i - 1) = 0.d0
                                                end if

                                                if (psip(n4 + i - 1) .gt. 0.d0) then
                                                    psip(n4 + i - 1) = dsqrt(psip(n4 + i - 1))
                                                else

                                                    psip(n4 + i - 1) = 0.d0
                                                    if (friction .gt. 0.d0) then
                                                        if (rank .eq. 0) write (6, *)&
                                                           & ' There should be some error in reweight0', i, psip(n4 + i - 1)
                                                        errnoise = 5
                                                    end if
                                                end if
                                            end do

                                            ! now we are able to update the velocity
                                            ! first change velocity and force into space which diagonalizes gamma
                                            call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin &
                                                       , velion(3, 1), 3, 0.d0, sov4, 1)
                                            call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin &
                                                       , sov4(indi + 1, 2), 1, 0.d0, velion(2, 1), 3)

                                            ! Second half iteration for Born-Oppenheimer damping
                                            do i = 1, ieskin
                                                zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                                                sov4(i, 1) = sov4(i, 1)*dexp(-psip(i)*dth) &
                                                             + psip(n5 + i - 1)*(velion(2, i) + psip(n4 + i - 1)*zeta)
                                            end do

                                            ! go back in the original (physical) basis

                                            call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin &
                                                       , sov4, 1, 0.d0, velion(3, 1), 3)

                                            ! Here sov4(*,1) contains the transformed velocities
                                            call dcopy(ieskin, velion(3, 1), 3, sov4, 1)

                                            if (yesquantum) then
                                                ! quantum harmonic propagation

                                                ! scale positions by masses
                                                do j = 1, nion
                                                    do i = 1, 3
                                                        ind = i + (j - 1)*3
                                                        sov4(ind, 2) = rion(i, j)/psip(n3 + ind - 1)
                                                    end do
                                                end do

#ifdef  PARALLEL
! In the quantum case we have also to apply the transformation kdyn that
! diagonalizes the interaction between the beads. Since each bead is in a
! different processor, this trasformation is done in parallel only
! kdyn matrix elements are common to all processors.
! Transformation coordinates in the basis of eigenvectors

                                                do i = 0, nbead - 1
                                                    sov5(1:ieskin, 1) = kdyn(rankcolrep + 1, i + 1)&
                                                        &*sov4(1:ieskin, 1) !  velocities
                                                    sov5(1:ieskin, 2) = kdyn(rankcolrep + 1, i + 1)&
                                                        &*sov4(1:ieskin, 2) !  coordinates
                                                    call reduce_base_real_to(2*ieskin, sov5, sov5(1, 6), commcolrep_mpi, i)
                                                end do
! sov5(i,6) ---> velocities
! sov5(i,7) ---> coordinates

#endif
#ifdef DEBUG
                                if (rank .eq. 0) write (6, *) ' changed velocities '&
                                    &, sum(sov5(4:ieskin, 6)), sum(sov5(4:ieskin, 7))
                  if (rank .eq. 0) write (6, *) ' changed velocities squares '&
                      &, sum(sov5(4:ieskin, 6)**2), sum(sov5(4:ieskin, 7)**2)
#endif

                                                omega_harm = sqrt(abs(kdyn_eig(rankcolrep + 1)))
                                                cosit = cos(omega_harm*dt)
                                                sinut = sin(omega_harm*dt)
                                                sinut_times_omega = omega_harm*sinut
                                                if (abs(omega_harm) .lt. 1d-6) then
                                                    sinut_over_omega = dt
                                                else
                                                    sinut_over_omega = sin(omega_harm*dt)/omega_harm
                                                end if

#ifdef DEBUG
                                                do i = 4, ieskin
#else
                                                    do i = 1, ieskin
#endif

                                                        ! exact integration of harmonic equations of motion for massless particles
                                                        sov4(i, 2) = cosit*sov5(i, 7) + sinut_over_omega*sov5(i, 6) & ! positions
                                                                     - sov5(i, 7) ! mind this is the variation!
                                                        velion(2, i) = cosit*sov5(i, 6)&
                                                           & - sinut_times_omega*sov5(i, 7) ! velocities

                                                    end do

#ifdef DEBUG
                                                    velion(1:2, 1:3) = 0.d0
                                                    if (rank .eq. 0) write (6, *) ' changed velocities 2 ', &
                                                        sum(velion(2, 1:ieskin))
                                                    if (rank .eq. 0) write (6, *) ' changed velocities squares 2  ', &
                                                        sum(velion(2, 1:ieskin)**2)
#endif
                                                    !       Now go back to the original basis
                                                    !       Tmes=dnrm2(ieskin,velion,3)**2/ieskin*ris(2)
                                                    !       call reduce_base_real(1,Tmes,commcolrep_mpi,-1)
                                                    !       if(rank.eq.0)  write(6,*) ' sum after0',Tmes/nbead
                                                    !    Transformation back in
                                                    !        the original basis from the normal modes of harmonic matrix
#ifdef PARALLEL
                                                    do i = 0, nbead - 1
                                                        sov5(1:ieskin, 1) = kdyn(i + 1, rankcolrep + 1)&
                                                            &*sov4(1:ieskin, 2) ! coordinates
                                                        sov5(1:ieskin, 2) = kdyn(i + 1, rankcolrep + 1)&
                                                            &*velion(2, 1:ieskin) !  velocities
                                                        call reduce_base_real_to(2*ieskin, sov5, sov5(1, 4), commcolrep_mpi, i)
                                                    end do
                                                    sov4(1:ieskin, 2) = sov5(1:ieskin, 4)
                                                    velion(2, 1:ieskin) = sov5(1:ieskin, 5)
#endif

                                                    else
                                                    ! classical propagation

                                                    do i = 1, ieskin
                                                        sov4(i, 2) = dt*sov4(i, 1)
                                                        velion(2, i) = sov4(i, 1)
                                                    end do

                                                    end if

                                                    call dcopy(ieskin, velion(2, 1), 3, velion(3, 1), 3)
                                                    call dcopy(ieskin, sov4(1, 2), 1, psip(n2 + indi), 1)

                                                    !       back to the physical coordinates (scaling back by masses)
                                                    do i = 1, ieskin
                                                        psip(n2 + indi + i - 1) = psip(n2 + indi + i - 1)*psip(n3 + i - 1)
                                                        if (psip(n3 + i - 1) .eq. 0.d0) then
                                                            velion(3, i) = 0.d0
                                                        end if
                                                    end do
                                                    cost = dnrm2(ieskin, psip(n2 + indi), 1)
                                                    if (rank .eq. 0) write (6, *) ' Norm change ions =', cost

                                                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    ! end idyn.eq.8

                                                elseif (idyn .eq. 6) then
                                                    if (yesrootc) then
                                                        cost = (delta0/2.d0)**2
                                                    else
                                                        if (temp .ne. 0.d0) then
                                                            cost = delta0*0.5d0/temp
                                                        else
                                                            cost = 0.d0
                                                        end if
                                                    end if

                                                    call dscal(ieskin2, cost, cov, 1)

                                                    if (yessecond) then
                                                        dth = dt/2.d0
                                                    else
                                                        dth = 0.d0
                                                    end if

                                                    indin = 0
                                                    do j = npps, np
                                                        indin = indin + 1
                                                        psip(n3 + indin - 1) = dsqrt(scalpar(ipip(n3 + j - 1))/dt)
                                                    end do

                                                    if (indin .ne. ieskin) then
                                              if (rank .eq. 0) write (6, *) ' Error in dynamic ieskin ne indin !!!  '&
                                                  &, ieskin, indin
                                                        errnoise = 6
                                                    end if

                                                    do j = 1, ieskin
                                                        do i = 1, ieskin
                                                 cov(ieskin*(j - 1) + i)&
                                                    & = cov(ieskin*(j - 1) + i)*psip(n3 + i - 1)*psip(n3 + j - 1)
                                                        end do
                                                    end do

                                                    !        first diagonalize sov(1,1,4)
                                                    !        save psip the QMC forces
                                                    call dcopy(np, psip, 1, sov4(1, 2), 1)
                                                    !        scale the force by the Mass
                                                    do i = 1, ieskin
                                                        sov4(indi + i, 2) = sov4(indi + i, 2)*psip(n3 + i - 1)
                                                    end do

                                                    !        do i=1,ieskin
                                                    !        cov(ieskin*(i-1)+i)=cov(ieskin*(i-1)+i)+friction
                                                    !        enddo

                                                    lworkd = 3*ieskin

                                                    !     if((yesquantum.and.nrep_bead.gt.1)) then
                                            call dsyev_my('V', 'L', ieskin, cov&
                                                &, ieskin, psip, info, nprocrep&
                                                &, rankrep, commrep_mpi)

                                                    if (yesrootc) then
                                                        do i = 1, ieskin
                                                            psip(i) = dsqrt(max(psip(i), 0.d0)) + friction
                                                        end do
                                                    else
                                                        do i = 1, ieskin
                                                            psip(i) = psip(i) + friction
                                                        end do
                                                    end if

                                                    !     else
                                                    !     call dsyev('V','L',ieskin,cov,ieskin,psip,psip(n5),lworkd,info)
                                                    !     endif
#ifdef PARALLEL
!      To avoid roundoff all processors have the same cov,psip
                                                    if (yesturboq) then
                                                        call bcast_real(cov, ieskin2, 0, mpi_comm_world)
                                                        call bcast_real(psip, ieskin, 0, mpi_comm_world)
                                                    end if
#endif

                                                    if (rank .eq. 0) then
                                                        write (6, *) ' Eigenvalues covariance '
                                                        do i = 1, ieskin
                                                            write (6, *) i, psip(i)
                                                        end do
                                                    end if
#ifdef DEBUG
!  Z2   Gauge fixing
                                                    do i = 1, ieskin
                                                        cost = abs(cov(ieskin*(i - 1) + 1))
                                                        indm = 1
                                                        do j = 2, ieskin
                                                        if (abs(cov(ieskin*(i - 1) + j)) .gt. cost) then
                                                            cost = abs(cov(ieskin*(i - 1) + j))
                                                            indm = j
                                                        end if
                                                        end do
                                                        if (cov(indm + (i - 1)*ieskin) .lt. 0.d0)
                                                            cov(1 + ieskin*(i - 1):ieskin*i) = &
                                                                &-cov(1 + ieskin*(i - 1):ieskin*i)
                                                    end do
#endif

                                                    if (info .ne. 0) then
                                                        if (rank .eq. 0) write (6, *) ' Error in lapack dsyev  dynamic !!! '
                                                        errnoise = 4
                                                    end if

                                    if (rank .eq. 0) write (6, *) ' Ratio dyn =',&
                                       &(psip(ieskin) - friction)/sqrt(temp*157.8873306d0)

                                                    ! Now we go in the basis that diagonalizes cov
                                                    ! NB for quantum dynamic (yesturboq) cov is assumed the same for all beads.
                                                    call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                                                            &, velion(3, 1), 3, 0.d0, sov4, 1)
                                                    call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                                                            &, sov4(indi + 1, 2), 1, 0.d0, velion(2, 1), 3)
                                                    !      Here sov4(*,1) contains the transfomed velocities
                                                    !           velion(2,*) the transformed forces
                                                    !    In the following velion(1,*) is unchanged

                                                    if (yesturboq) then

                                                        if (.not. yessecond .and. sum(abs(velion(1, 1:ieskin))) .ne. 0.d0) then
                                                   write (6, *) ' ERROR velion 1 not conserved !!!  '&
                                                   &, sum(abs(velion(1, 1:ieskin)))
                                                            errnoise = 23
                                                        end if
                                                        !  Transformation coordinates in the basis of eigenvectors
                                                        !  velion(1,*) represents run time
                                                        !  the difference between the atomic position
                                                        !  when the forces are computed
                                                        !  and the coordinates when we start to apply
                                                        !    the new LD integrator. velion(1,*) is non zero only for second order
                                                        !    dynamics (yessecond), more accurate with the time step integration.
                                                        !    Notice that, after that,
                                                        !    we have to scale by the sqrt
                                                        !    of the mass (psip(n3)~1/sqrt(m))
                                                        do j = 1, nion
                                                            do i = 1, 3
                                                                ind = i + (j - 1)*3
                                                                sov5(ind, 1) = (rion(i, j) - velion(1, ind))/psip(n3 + ind - 1)
                                                            end do
                                                        end do
                                                        call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                                                                &, sov5, 1, 0.d0, sov5(1, 4), 1)
#ifdef  PARALLEL
! In the quantum case we have also to apply the transformation kdyn that
! diagonalizes the interaction between the beads. Since each bead is in a
! different processor, this trasformation is done in parallel only
! kdyn matrix elements are common to all processors.
                                                        do i = 0, nbead - 1
                                                            !  velocities
                                                            sov5(1:ieskin, 1) = kdyn(rankcolrep + 1, i + 1)*sov4(1:ieskin, 1)
                                                            !  forces
                                                            sov5(1:ieskin, 2) = kdyn(rankcolrep + 1, i + 1)*velion(2, 1:ieskin)
                                                            !  coordinates
                                                            sov5(1:ieskin, 3) = kdyn(rankcolrep + 1, i + 1)*sov5(1:ieskin, 4)
                                                            call reduce_base_real_to(3*ieskin, sov5, sov5(1, 5)&
                                                                &, commcolrep_mpi, i)
                                                        end do
                                                        sov4(1:ieskin, 1) = sov5(1:ieskin, 5)
                                                        velion(2, 1:ieskin) = sov5(1:ieskin, 6)
#endif
                                                    end if
#ifdef DEBUG
                                                    if (rank .eq. 0) write (6, *) ' changed velocities ',&
                                                      &sum(velion(2, 4:ieskin)), sum(sov4(4:ieskin, 1))
                                                    if (rank .eq. 0) write (6, *) ' changed velocities squares ',&
                                                    &sum(velion(2, 4:ieskin)**2), sum(sov4(4:ieskin, 1)**2)
#endif
                                                    alphaqmc = 0.d0
                !!!! INIZIO NEW DYN : Now everything is in the diagonal basis and the integration
                                                    !can be done exactly for each of the ieskin (x nbead , index labelled by
                                                    ! rankcolrep+1 in the quantum case) eigenvalues

#ifdef DEBUG
                                                    do i = 4, ieskin
#else
                                                        do i = 1, ieskin
#endif
                                                            zetan(1) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())
                                                            zetan(2) = dsqrt(-2.d0*dlog(1.d0 - drand1()))*dcos(2.d0*pi*drand1())

                                                            if (psip(i) .gt. 0.d0) then
                                                                if (delta0 .ne. 0.d0 .and. normcorr .gt. 0.d0) then
                                                                    if (yesrootc) then
                                                                        alphaqmc = normcorr*(2.d0*(psip(i) - friction)/delta0)**2
                                                                    else
                                                                        alphaqmc = 2.d0*normcorr*temp*(psip(i) - friction)/delta0
                                                                    end if
                                                                else
                                                                    alphaqmc = 0.d0
                                                                end if
                                                                alphaall = 2.d0*temp*psip(i)

                                                                if (yesturboq) then

                                                                    ! gammaall
                                                                    psip(i) = psip(i)&
                                                                       & + (psip(i) - friction)&
                                                                       & * delta0q*kdyn_eig(rankcolrep + 1) + &
                                                                            &2.d0*delta0k*sqrt(abs(kdyn_eig(rankcolrep + 1)))

                                                                    alphaall = 2*temp*psip(i)

                                                                    call set_turboq(psip(i), kdyn_eig(rankcolrep + 1), dt&
                                                                            &, alphaall, alphaqmc, mnoise, Gn, Tn, Gni, Gnh, Gnih)

                                                                    if (.not. yessecond) then
                                                                        Gnh = 0.d0
                                                                        Gnih = 0.d0
                                                                    end if

                                                                else

                                                                    cost = exp(-dt*psip(i))
                                                                    if (cost .lt. 0.99999d0) then
                                                                        costm1 = cost - 1.d0
                                                                    else
                                                                        costm1 = -dt*psip(i) + (dt*psip(i))**2/2.d0
                                                                    end if
                                                                    costh = exp(-dth*psip(i))
                                                                    if (costh .lt. 0.99999d0) then
                                                                        costhm1 = costh - 1.d0
                                                                    else
                                                                        costhm1 = -dth*psip(i) + (dth*psip(i))**2/2.d0
                                                                    end if

                                                                    Gni(1, 1) = cost
                                                                    Gn = -costm1/psip(i)
                                                                    if (yessecond) then
                                                                        Gnh = -costhm1/psip(i)
                                                                    else
                                                                        Gnh = 0.d0
                                                                    end if
                                                                    Tn = 1.d0/psip(i)*(dt - Gn)

                                                                    mnoise(1, 1) = -alphaall/(2.d0*psip(i))*costm1*(1 + cost)

                                                                    mnoise(2, 2) = alphaall/psip(i)**2* &
                                                                            &(dt + 1.d0/psip(i)&
                                                                            &*costm1*(2.d0 - 0.5d0*(1.d0 + cost)))

                                                                    mnoise(1, 2) = 0.5d0*alphaall*(costm1/psip(i))**2

                                                                    mnoise(1, 1) = mnoise(1, 1)/Gn**2 - alphaqmc
                                                                    mnoise(2, 2) = mnoise(2, 2)/Tn**2 - alphaqmc
                                                                    mnoise(1, 2) = mnoise(1, 2)/Tn/Gn - alphaqmc

                                                                    mnoise(2, 1) = mnoise(1, 2)

                                                                end if

                                                                all_dyn = all_dyn + mnoise(1, 1) + mnoise(2, 2)
                                                                ratio_dyn = ratio_dyn + 2*alphaqmc

                                                                call root2mat(mnoise, errnoise)

                                                            else
                                                                if (rank .eq. 0) write (6, *) ' ERROR negative eigenv. ', psip(i)&
                                                                        &, alphaqmc
                                                                errnoise = 1
                                                            end if !if psip(i).gt.0

                                        if (errnoise .ne. 0 .and. rank .eq. 0) write (6, *)&
                                           & ' Error negative definite matrix noise '

                                                            eta_v = mnoise(1, 1)*zetan(1) + mnoise(1, 2)*zetan(2)
                                                            eta_r = mnoise(2, 1)*zetan(1) + mnoise(2, 2)*zetan(2)

                                                            !        Now the consistent change of coordinates put in velion(1,*)
                                                            !       velion(3,*) represents an extra
                                                            !       displacement of the positions by dt/2 without the noise
                                                            !       and the force (Gn=0), i.e. consistent up to order dt^3/2.

                                                            if (yesturboq) then
                                                sov4(i, 2) = sov4(i, 1)*Gni(2, 1) + Gni(2, 2)*sov5(i, 7)&
                                                   & + Tn*(velion(2, i) + eta_r)
                                              velion(2, i) = sov4(i, 1)*Gni(1, 1) + Gni(1, 2)*sov5(i, 7)&
                                                 & + Gn*(velion(2, i) + eta_v)
                                                                !      Propagate further by half dt
                                                                !      without noise and force, if necessary
                                                                !      (yessecond=.true.) in terms
                                                                !      of the new velocities and coordinates.
                                                                velion(3, i) = Gnh*velion(2, i) + Gnih*(sov5(i, 7) + sov4(i, 2))
                                                            else
                                                                sov4(i, 2) = sov4(i, 1)*Gn + Tn*(velion(2, i) + eta_r)
                                                                velion(2, i) = sov4(i, 1)*Gni(1, 1) + Gn*(velion(2, i) + eta_v)
                                                                !      Propagate further by half dt without noise and forces
                                                                velion(3, i) = Gnh*velion(2, i)
                                                            end if

                                                        end do
#ifdef DEBUG
                                                        velion(1:2, 1:3) = 0.d0
                                                        if (rank .eq. 0) write (6, *) ' changed velocities 2 ',&
                                                 &sum(velion(2, 1:ieskin)), sum(velion(1, 1:ieskin))
                                                        if (rank .eq. 0) write (6, *) ' changed velocities squares 2  ',&
                                                        &sum(velion(2, 1:ieskin)**2), sum(velion(1, 1:ieskin)**2)
#endif

                                                        !         Now go back to the original basis
                                                        if (yesturboq) then
                                                            !       Tmes=dnrm2(ieskin,velion,3)**2/ieskin*ris(2)
                                                            !       call reduce_base_real(1,Tmes,commcolrep_mpi,-1)
                                                            !       if(rank.eq.0)  write(6,*) ' sum after0',Tmes/nbead

                                                            !    Transformation back in the original basis
#ifdef  PARALLEL
                                                            do i = 0, nbead - 1
                                                                sov5(1:ieskin, 1) = kdyn(i + 1, rankcolrep + 1)&
                                                                    &*sov4(1:ieskin, 2) ! coordinates
                                                                sov5(1:ieskin, 2) = kdyn(i + 1, rankcolrep + 1)&
                                                                    &*velion(2, 1:ieskin) !  velocities
                                                                sov5(1:ieskin, 3) = kdyn(i + 1, rankcolrep + 1)&
                                                                    &*velion(3, 1:ieskin) ! extra prop coordinates
                                                             call reduce_base_real_to(3*ieskin, sov5, sov5(1, 4)&
                                                                 &, commcolrep_mpi, i)
                                                            end do
                                                            sov4(1:ieskin, 2) = sov5(1:ieskin, 4)
                                                            velion(2, 1:ieskin) = sov5(1:ieskin, 5)
                                                            velion(3, 1:ieskin) = sov5(1:ieskin, 6)

!       Tmes=dnrm2(ieskin,velion,3)**2/ieskin*ris(2)
!       call reduce_base_real(1,Tmes,commcolrep_mpi,-1)
!       if(rank.eq.0)  write(6,*) ' sum after1',Tmes/nbead
#endif
                                                        end if

                                                        call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin               &
                                                                &, velion(3, 1), 3, 0.d0, sov4, 1)
                                                        call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin               &
                                                                &, velion(2, 1), 3, 0.d0, velion(3, 1), 3)
                                                        call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin               &
                                                                &, sov4(1, 2), 1, 0.d0, psip(n2 + indi), 1)

                                                        !       back to the scale of coordinates
                                                        !       write(6,*) '  Final correction  inside '
                                                        cost = 0.d0
                                                        do i = 1, ieskin
                                                            psip(n2 + indi + i - 1)&
                                                              & = (psip(n2 + indi + i - 1)&
                                                              & + sov4(i, 1))*psip(n3 + i - 1) - velion(1, i)
                                                            velion(1, i) = sov4(i, 1)*psip(n3 + i - 1)
                                                            if (psip(n3 + i - 1) .ne. 0.d0) then
                                                                cost = cost + 1.d0
                                                            else
                                                                velion(3, i) = 0.d0
                                                            end if
                                                        end do
                                                        Tmes = dnrm2(ieskin, velion(3, 1), 3)**2/cost*ris(2)
                                                        if (rank .eq. 0) write (6, *) ' Temperature (H)= ', Tmes
                                                        cost = dnrm2(ieskin, psip(n2 + indi), 1)
                                                        if (rank .eq. 0) write (6, *) ' Norm change ions =', cost

                !!!! FINE  NEW IDYN=6

                                                        elseif (idyn .ge. 2) then
                                                        !         simple second order Langevin dynamic assumed correct up
                                                        !         to second order
                                                        !         first compute covariance matrix and put it in sov(*,*,4)
                                                        !

                                                        ! by E. Coccia (18/1/11): allow rototranslations
                                                        ! if an external potential is added
                                                        !            if (.not.ext_pot.and..not.iespbc) then
                                                        !            call cleancov
                                                        !            endif
                                                        !
                                                        !     Restoring the indexing  of velocities consistent with the one of
                                                        !     forces
                                                        !      indin=0
                                                        !      do i=1,npp
                                                        !      if(ipip(i).eq.2) then
                                                        !      indin=indin+1
                                                        !      j=ipip(n3+i-1)-kp_ion
                                                        !      velion(2,indin)=velion(3,j)
                                                        !      endif
                                                        !      enddo
                                                        !      velion(3,1:ieskin)=velion(2,1:ieskin)
                                                        !      if(indin.ne.ieskin) then
                                                        !      write(6,*) ' ERROR idyn>1 does not work inconsitent # forces '&
                                                        !     &,indin,ieskin
                                                        !      errnoise=13
                                                        !      endif

                                                        if (yesrootc) then
                                                            cost = (delta0/2.d0)**2
                                                        else
                                                            if (temp .ne. 0.d0) then
                                                                cost = delta0*0.5d0/temp
                                                            else
                                                                cost = 0.d0
                                                            end if
                                                        end if

                                                        !         do i=1,ieskin
                                                        !         call dcopy(ieskin,cov(ieskin*(i-1)+1),1,sov4(1,i),1)
                                                        !         call dscal(ieskin,cost,sov4(1,i),1)
                                                        !         enddo
                                                        call dscal(ieskin2, cost, cov, 1)

                                                        !         call dgemm('T','N',ieskin,ieskin,nbin,cost,fk(1,kp_ion+1)
                                                        !    1,nbin,fk(1,kp_ion+1),nbin,0.d0,sov(1,1,4),npm)

                                                        !        Mass correction
                                                        !        write(6,*) ' Mass correction ',dt
                                                        indin = 0
                                                        do j = npps, np
                                                            indin = indin + 1

                                                            psip(n3 + indin - 1) = dsqrt(scalpar(ipip(n3 + j - 1))/dt)
                                                            !        write(6,*) ' Inverse Mass =',indin,ipip(n3+j-1)
                                                            !                                    ,psip(n3+indin-1)**2
                                                        end do

                                                        do j = 1, ieskin
                                                            do i = 1, ieskin
                                                                cov(ieskin*(j - 1) + i)&
                                                                   & = cov(ieskin*(j - 1) + i)&
                                                                   &*psip(n3 + i - 1)*psip(n3 + j - 1)
                                                            end do
                                                        end do

                                                        ! MC
                                                        !        scale the force by the Mass
                                                        call dcopy(np, psip, 1, sov4(1, 2), 1)
                                                        do i = 1, ieskin
                                                            !         sov4(indi+i-1,2)=sov4(indi+i-1,2)*psip(n3+i-1)
                                                            sov4(indi + i, 2) = sov4(indi + i, 2)*psip(n3 + i - 1)
                                                        end do
                                                        ! MC

                                                        !         write(6,*) ' Old/now  velocity =',velion(2,1),velion(3,1)

                                                        !         calculation extra noise

                                                        !         do i=1,ieskin
                                                        !         sov(i,i,4)=sov(i,i,4)+friction
                                                        !         enddo

                                                        !        save previouus velocity
                                                        !        call dcopy(ieskin,velion(3,1),3,velion(1,1),3)

                                                        !         calculation -gamma  v_n  velion(3,i),put it in sov(1,npp,4)
                                                        if (idyn .eq. 2) then

                                                            call dgemv('N', ieskin, ieskin, -1.d0, cov, ieskin             &
                                                                    &, velion(3, 1), 3, 0.d0, sov4, 1)

                                                            do i = 1, ieskin
                                                                sov4(i, 1) = sov4(i, 1) - friction*velion(3, i)
                                                            end do

                                                            !        update velion(3,1)

                                                            cost = 2.d0*friction*temp*dt

                                                            do i = 1, ieskin
                                                                !        write(6,*) ' force inside ',psip(indi+i),psip(n3+i-1)
                                                               zeta = dsqrt(-2.d0*dlog(1.d0 - drand1())*cost)&
                                                                   &*dcos(2.d0*pi*drand1())
                                                                !        write(6,*) ' zeta corrected ',zeta
                                                                !        II order unstable
                                                                !        velion(3,i)=velion(2,i)+2.d0*dt*
                                                                !    1  (psip(indi+i)*psip(n3+i-1)+sov(i,npp,4))+zeta
                                                                velion(3, i) = velion(3, i) + dt* &
                                                                        &  (psip(indi + i)*psip(n3 + i - 1) + sov4(i, 1)) + zeta
                                                            end do

                                                            !        Now the variation of coordinates up to order dt**2

                                                            !        calculation temperature
                                                            !        (unit mass assumed for all velocities)
                                                            Tmes = dnrm2(ieskin, velion(3, 1), 3)**2/ieskin*ris(2)

                                                        elseif (idyn .eq. 3) then
                                                            !        first diagonalize sov(1,1,4)
                                                            if (yessecond) then
                                                                dth = dt/2.d0
                                                            else
                                                                dth = dt
                                                            end if

                                                            !        do i=1,ieskin
                                                            !        cov(ieskin*(i-1)+i)=cov(ieskin*(i-1)+i)+friction
                                                            !        enddo
                                                            lworkd = 3*ieskin

                                                            call dsyev_my('V', 'L', ieskin, cov, ieskin, psip, info, nprocrep&
                                                                    &, rankrep, commrep_mpi)

                                                            !      call dsyev('V','L',ieskin,cov,ieskin,psip,psip(n5),lworkd,info)
                                                            !
                                                            if (yesrootc) then
                                                                do i = 1, ieskin
                                                                    psip(i) = dsqrt(max(psip(i), 0.d0)) + friction
                                                                end do
                                                            else
                                                                do i = 1, ieskin
                                                                    psip(i) = psip(i) + friction
                                                                end do
                                                            end if

                                                            if (rank .eq. 0) then
                                                                write (6, *) ' Eigenvalues gamma chosen '
                                                                do i = 1, ieskin
                                                                    write (6, *) i, psip(i)
                                                                end do
                                                            end if

                                                            if (info .ne. 0) then
                                                                if (rank .eq. 0) write (6, *)&
                                                                   & ' Error in lapack dsyev  dynamic !!! '
                                                                errnoise = 4
                                                            end if

                                                            ! The scale invariant criterium
                                                            ! is gamma/sqrt(T) , so I have divided by the
                                                            ! sqrt(T/T_0) where T_0=1000K

                                                            if (rank .eq. 0) write (6, *)&
                                                               & ' Ratio dyn ='&
                                                               &, (psip(ieskin) - friction)/sqrt(temp*157.8873306d0)

                                                            !         write(6,*) ' Spectrum gamma=',info,(psip(i)*dt,i=1,ieskin)

                                                            do i = 1, ieskin
                                                                if (psip(i) .gt. 0.d0) then
                                                                    psip(n5 + i - 1) = (1.d0 - dexp(-psip(i)*dth))/psip(i)
                                                                else
                                                                    if (rank .eq. 0) write (6, *) ' Refused eigenvalue ', i
                                                                    psip(n5 + i - 1) = dth
                                                                end if
                                                            end do
                                                            !         compute the noise correction
                                                            alphaqmc = 0.d0

                                                            do i = 1, ieskin
                                                                if (psip(i) .gt. 0.d0) then
                                                                    !         psip(n4+i-1)=(psip(i)*2.d0*temp)*psip(i)*
                                                                    !    1  sinh(dt*psip(i))/(4.d0*sinh(psip(i)*dt*.5d0)**2)
                                                                    !         the following is protected for overflow
                                                                    cost = exp(-dth*psip(i))
                                                                    if (delta0 .ne. 0.d0 .and. normcorr .gt. 0.d0) then
                                                                        if (yesrootc) then
                                                                           alphaqmc = normcorr&
                                                                               &*(2.d0*(psip(i) - friction)/delta0)**2
                                                                        else
                                                                           alphaqmc =&
                                                                           &normcorr*2.d0*temp*(psip(i) - friction)/delta0
                                                                        end if
                                                                    else
                                                                        alphaqmc = 0.d0
                                                                    end if
                                                                    psip(n4 + i - 1) = psip(i)*temp*psip(i)* &
                                                                            &  (1.d0 + cost)/(1.d0 - cost) - alphaqmc
                                                                    !          subtracting the noise already present in the forces
                                                                else
                                                                    psip(n4 + i - 1) = 0.d0
                                                                end if

                                                                all_dyn = all_dyn + psip(n4 + i - 1)
                                                                ratio_dyn = ratio_dyn + alphaqmc

                                                                if (psip(n4 + i - 1) .gt. 0.d0) then
                                                                    psip(n4 + i - 1) = dsqrt(psip(n4 + i - 1))
                                                                else
                                                                    if (yesrootc) write (6, *)&
                                                                       & ' Warning noise correction not possible '&
                                                                       &, i, psip(i), psip(n4 + i - 1)
                                                                    psip(n4 + i - 1) = 0.d0
                                                                    if (friction .gt. 0.d0 .and. .not. yesrootc) then
                                                                        if (rank .eq. 0) write (6, *)&
                                                                           & ' There should be some error in reweight0         &
                                                                                &  ', i, psip(n4 + i - 1)
                                                                        errnoise = 5
                                                                    end if
                                                                end if
                                                            end do

                                                            !         now we are able to update the velocity
                                                            !         first change basis actual velocity and force
                                                            call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                                                                    &, velion(3, 1), 3, 0.d0, sov4, 1)
                                                            call dgemv('T', ieskin, ieskin, 1.d0, cov, ieskin              &
                                                                    &, sov4(indi + 1, 2), 1, 0.d0, velion(2, 1), 3)

                                                            if (yessecond) then
                                                                ! Compute the temperature at half time intervals.
                                                                do i = 1, ieskin
                                                                    zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))&
                                                                        &*dcos(2.d0*pi*drand1())
                                                                    sov4(i, 1) = sov4(i, 1)*dexp(-psip(i)*dth)&
                                                                         & + psip(n5 + i - 1)*(velion(2, i)&
                                                                         & + psip(n4 + i - 1)*zeta)
                                                                end do
                                                                Tmes = dnrm2(ieskin, sov4, 1)**2/ieskin*ris(2)
                                                                ! Second  half time interval.
                                                                do i = 1, ieskin
                                                                    zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))&
                                                                                     &*dcos(2.d0*pi*drand1())
                                                                    velion(2, i) = sov4(i, 1)*dexp(-psip(i)*dth) &
                                                                         & + psip(n5 + i - 1)*(velion(2, i)&
                                                                         & + psip(n4 + i - 1)*zeta)
                                                                end do
                                                            else
                                                                do i = 1, ieskin
                                                                    zeta = dsqrt(-2.d0*dlog(1.d0 - drand1()))&
                                                                        &*dcos(2.d0*pi*drand1())
                                                                    velion(2, i) = sov4(i, 1)*dexp(-psip(i)*dth)&
                                                                          & + psip(n5 + i - 1)*(velion(2, i)&
                                                                          & + psip(n4 + i - 1)*zeta)
                                                                end do
                                                            end if

                                                            !         go back in the original basis

                                                            call dgemv('N', ieskin, ieskin, 1.d0, cov, ieskin              &
                                                                    &, velion(2, 1), 3, 0.d0, velion(3, 1), 3)

                                                                if (.not. yessecond) Tmes = dnrm2(ieskin, velion(3, 1), 3)**2&
                                                                    &/ieskin*ris(2)
                                                            if (rank .eq. 0) write (6, *) ' Temperature (H)= ', Tmes

                                                            ! endif idyn=3
                                                        end if

                                                        indin = 0
                                                        do i = npps, npp
                                                            !        II order half integer times for velocities
                                                            indin = indin + 1
                                                            if (yessecond .and. idyn .eq. 3) then
                                                                !        Old half time step
                                                                psip(indi + n2 + indin - 1) = velion(1, indin)
                                                                !        New half time step
                                                                velion(1, indin) = dth*psip(n3 + indin - 1)*velion(3, indin)
                                                                !        Full time step
                                                        psip(indi + n2 + indin - 1) = psip(indi + n2 + indin - 1)&
                                                        &+ velion(1, indin)
                                                            else
                                                              psip(indi + n2 + indin - 1) = dt*psip(n3 + indin - 1)&
                                                                   &*velion(3, indin)
                                                            end if
                                                            !         write(6,*) 'TEST',indin,velion(3,indin)&
                                                            !                    ,psip(n3+indin-1),psip(indi+n2+indin-1)
                                                        end do

                                                        cost = dnrm2(ieskin, psip(n2 + indi), 1)
                                                        if (rank .eq. 0) write (6, *) ' Norm change ions =', cost
                                                        end if ! endif idyn.ge.2
                                                        if (rank .eq. 0 .and. all_dyn .gt. 0) write (6, *) &
                                                                &  ' Ratio QMC noise/ALL noise', ratio_dyn/all_dyn
                                                        ! do not move the ions and do not update velocities
                                                    else

                                                        if (maxdev_dyn .ne. 0&
                                                           & .and. devmaxp .gt. maxdev_dyn &
                                                           & .and. idyn .gt. 0) then
                            if (rank .eq. 0) write (6, *) ' Warning devmax too large Ions not moved this step '&
                                &, devmaxp, maxdev_dyn
                                                            acc_dyn = .false.
                                                        end if
                                                        do i = 1, ieskin
                                                            psip(indi + n2 + i - 1) = 0.d0
                                                        end do
                                                    end if ! endif idyn.ge.0
                                                    !        END ION DYNAMICS
                                                    if (allocated(sov4)) deallocate (sov4)
                                                    if (allocated(sov5)) deallocate (sov5)
                                                    end subroutine ion_dynamics

                                                    subroutine cleancov
                                                        implicit none
                                                        integer indx, indy, indr
                                                        real*8 scalprod
                                                        real(8), dimension(:), allocatable :: vecprod
                                                        real*8, external :: dlamch
                                                        ! if not in the box and all the forces are calculated:
                                                        ! purification of matrix
                                                        ! cov from translations and rotations of the center of mass

                                                        if (ieskin .ge. 3*niont) then
                                                            mineig = 4
                                                        else
                                                            mineig = 1
                                                        end if

                                                        if (.not. iespbc .and. ieskin .eq. 3*niont) then
                                                            allocate (riondcm(3, niont), vecrottr(3*niont, 6), &
                                                                   &projrottr(ieskin, ieskin), workGS(ieskin + 12)&
                                                                   &, vecprod(ieskin))
                                                            if (epstion .lt. 0.d0) allocate (covpurif(ieskin, ieskin))

                                                            !             write(6,*) ' force before ',psip(indi+1:indi+ieskin)
                                                            if (rank .eq. 0) write (6, *) 'Center of mass of the system: '
                                                            if (rank .eq. 0) write (6, *) rioncm

                                             if (rank .eq. 0) write (6, *) 'Ion coordinates from the Center of mass'//&
                                                &' of the system: '
                                                            riondcm = 0.d0
                                                            indr = 0
                                                            do i = 1, nion
                                                                if (atom_number(i) .gt. 0) then
                                                                    indr = indr + 1
                                                                    do j = 1, 3
                                                                        riondcm(j, indr) = rion(j, i) - rioncm(j)
                                                                    end do
                                                                    !                write(6,*) 'Ion ',i
                                                                    !                write(6,*) riondcm(:,i)
                                                                end if
                                                            end do

                                                            vecrottr = 0.d0
                                                            do j = 1, niont
                                                                !   3 rotation vectors of the ions
                                                                !                vecrottr(3*(j-1)+1,1)=  0.d0
                                                                vecrottr(3*(j - 1) + 2, 1) = -riondcm(3, j)
                                                                vecrottr(3*(j - 1) + 3, 1) = riondcm(2, j)
                                                                vecrottr(3*(j - 1) + 1, 2) = riondcm(3, j)
                                                                !                vecrottr(3*(j-1)+2,2)=  0.d0
                                                                vecrottr(3*(j - 1) + 3, 2) = -riondcm(1, j)
                                                                vecrottr(3*(j - 1) + 1, 3) = -riondcm(2, j)
                                                                vecrottr(3*(j - 1) + 2, 3) = riondcm(1, j)
                                                                !                vecrottr(3*(j-1)+3,3)=  0.d0
                                                                !    3 translation vectors of the ions
                                                                vecrottr(3*(j - 1) + 1, 4) = 1.d0
                                                                vecrottr(3*(j - 1) + 2, 5) = 1.d0
                                                                vecrottr(3*(j - 1) + 3, 6) = 1.d0
                                                            end do

                                                            !     do k=4,6
                                                            !        vecprod=0.d0
                                                            !        do i=1,ieskin
                                                            !           do j=1,ieskin
                                                            !              vecprod(i)=vecprod(i)+cov(i+ieskin*(j-1))*vecrottr(j,k)
                                                            !           enddo
                                                            !        enddo
            !!       write(6,*) 'cov*vecrottr',k
            !!       write(6,*) vecprod
                                                            !     enddo

                                                            ! check GM orthonorm
                                                            !        do j=1,6
                                                            !           write(6,*) ' vecrottr ',j,' before GS orthonormalization'
                                                            !           write(6,*) vecrottr(:,j)
                                                            !        enddo

                                                            ! check scalar product between rotations and translations
                                                            !     do i=1,3
                                                            !        do j=4,6
                                                            !           scalprod=0.d0
                                                            !           do k=1,ieskin
                                                            !              scalprod=scalprod+vecrottr(k,i)*vecrottr(k,j)
                                                            !           enddo
            !!          write(6,*) ' vecrot',i,' vectr',j-3,' = ',scalprod
                                                            !        enddo
                                                            !     enddo
                                                            !         if(rank.eq.0) then
                                                            !          write(6,*) ' Before ortho '
                                                            !            do j=1,6
                                                            !               write(6,*) ' vecrottr ',j&
                                                            !                        &,' after GS orthonormalization'
                                                            !               write(6,*) vecrottr(:,j)
                                                            !            enddo
                                                            !         endif
                                                            call grahamo(vecrottr, workGS(7), workGS, ieskin, 6, info)
                                                            !      info=0
                                                            !      do i=1,6
                                                            !       if(workGS(i).lt.dlamch('E')) then
                                                            !       info=info+1
                                                            !       vecrottr(:,i)=0.d0
                                                            !       endif
                                                            !      enddo

                                                            if (info .ne. 0) then

                                               if (rank .eq. 0) write (6, *) ' Warning found 6-', info&
                                                   &, ' independent rotations !! '
                                                  if (info .eq. 1 .and. rank .eq. 0) write (6, *) &
                                                      &' Warning collinear molecule !!! '
                                                        if (info .gt. 1) then
                                                            if (rank .eq. 0) write (6, *)&
                                                               & ' ERROR dependent roto-translations !!! '
                                                                   errnoise = 11
                                                            end if
                                                            end if
                                                            ! check GS orthonorm
                                                            !         if(rank.eq.0) then
                                                            !            do j=1,6
                                                            !               write(6,*) ' vecrottr ',j&
                                                            !                        &,' after GS orthonormalization'
                                                            !               write(6,*) vecrottr(:,j)
                                                            !            enddo
                                                            !            do j=1,6
                                                            !             do k=j,6
                                                            !             write(6,*) ' scal j k =',j,k,
                                                            !                        sum(vecrottr(:,j)*vecrottr(:,k))
                                                            !             enddo
                                                            !            enddo
                                                            !         endif

                                                            !   projector orthogonal space of the rotations
                                                            projrottr = 0.d0
                                                            do i = 1, ieskin
                                                                projrottr(i, i) = 1.d0
                                                                !               do j=1,ieskin
                                                                !                 do k=1,6
                                                                !                 projrottr(i,j) = projrottr(i,j)&
                                                                !                                &-vecrottr(i,k)*vecrottr(j,k)
                                                                !                 enddo
                                                                !               enddo
                                                            end do
                                                            !    Optimization
                                                  call dgemm('N', 'T', ieskin&
                                                      &, ieskin, 6, -1.d0&
                                                      &, vecrottr, ieskin&
                                                      &, vecrottr, ieskin&
                                                      &, 1.d0, projrottr, ieskin)
                                                            !
                                                            !        Purification of forces
                                                            vecprod(1:ieskin) = psip(indi + 1:indi + ieskin)

                                       call dgemv('N', ieskin, ieskin&
                                           &, 1.d0, projrottr, ieskin&
                                           &, vecprod, 1, 0.d0&
                                           &, psip(indi + 1), 1)

                                                            ! be biased.

                                                            !   purification of the cov matrix
                                                            if (epstion .lt. 0.d0) then
                                                 call dgemm('N', 'N', ieskin&
                                                      &, ieskin, ieskin, 1.d0&
                                                      &, cov, ieskin, projrottr&
                                                      &, ieskin, 0.d0, covpurif&
                                                      &, ieskin)
                                                call dgemm('N', 'N', ieskin&
                                                    &, ieskin, ieskin, 1.d0&
                                                    &, projrottr, ieskin&
                                                    &, covpurif, ieskin, 0.d0&
                                                    &, cov, ieskin)
                                                            end if
                                                            !     write(6,*) 'Purified cov matrix'
                                                            !     do i=1,ieskin
                                                            !        do j=1,ieskin
                                                            !           write(6,*) i,j,cov(i+ieskin*(j-1))
                                                            !        enddo
                                                            !     enddo
            !! check purificated cov orthogonal to vecrottr
                                                            !     do i=1,3
                                                            !        do j=4,6
                                                            !           scalprod=0.d0
                                                            !           do k=1,ieskin
                                                            !              scalprod=scalprod+vecrottr(k,i)*vecrottr(k,j)
                                                            !           enddo
                                                            !           if(rank.eq.0) write(6,*) ' vecrot'
                                                            !               ,i,' vectr',j-3,' = ',scalprod
                                                            !        enddo
                                                            !     enddo

                                                            mineig = 7 - info
                                                            mineig = min(mineig, ieskin) ! Cannot be larger than ieskin

                                                            deallocate (riondcm, vecrottr, projrottr, workGS, vecprod)
                                                            if (epstion .lt. 0.d0) deallocate (covpurif)
                                                        end if

                                                    end subroutine cleancov

                                                    end subroutine reweight0
