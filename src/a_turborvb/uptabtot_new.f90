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

subroutine uptabtot(nelup, neldo, nelorb, nelorbh, jel            &
        &, kel, winv, winvup, winvdo, ainv, psiln, psisn, epst                     &
        &, psip, rcart, dist, dists, ainvs, winvs, psinew                         &
        &, indt, indt4, indt4j, indtc, tabpip, alat, ivic, plat, rion, vpot, iesd, vj  &
        &, zeta, dup, nshell, nshelldo, ioptorb, ioccup, ioccdo, itest             &
        &, tmu, nion, r, rmu, kion, winvj, ioccj, kionj, vju, nelorbj, nelorbjh       &
        &, ioptorbj, nshellj, winvsj, winvbar, detmat                           &
        &, winvjbar, winvjbarn, winvjbarsz, winvjbarszn, jasmat, jasmatsz        &
        &, iflagnorm, cnorm, iflagtab                                         &
        &, ratio, ncore, lmax, nintpseudo, prefactor, rcutoff, parshell           &
        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
        &, jpseudo, pseudolocal, tcost, istart, costz, costz3, ioptpseudo, angle   &
        &, iessz, keln, oldkappa, LBox, rmucos, rmusin, walker, n_body_on          &
        &, jastrowall_ee, jasnew_ee, jastrowall_ei, jasnew_ei, timepip          &
        &, timewf, niesd, versoralat, iesrandoma, nshelltot, indtm, psidetdmc&
        &, mu_c, projm, nelorb_c, firstmol, nmol, yesfast, vpotreg, cutreg&
        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab       &
        &, indshellj_tab, detmat_c, winvbarfn, winvfn, contraction, vpotsav_ee)

    use Cell
    use Ewald
    use Constants, only: ip4, ipj, ipc, zzero, zone, nelorbh_ip, nmol_ip
    use allio, only: epscuttype, agp, agpn, jtype2, nelup_mat, iespbc&
     &, molecular, enforce_detailb, dim_uptabtot, yes_ontarget, pot_aasunel&
     &, yes_sparse, nnozeroj, nozeroj, nelorbjh2, jasmat_c, muj_c, nelorbj_c, contractionj
    use dielectric
    implicit none

    integer nelup, neldo, nel, jel, itest, nion, indtm(*), indjel        &
           &, nelorb, nelorbh, i, indt, indtc, kk, jnew, iesd, jn, ioptorb(*), ioccup(*)&
           &, ioccdo(*), nshell, nshelldo, j, nelorb5, k, kion(*)        &
           &, ioccj(*), kionj(*), nelorbj, nelorbjh, nshellj, ioptorbj(*), iflagnorm&
           &, iflagtab, ncore, lmax, istart, nintpseudo, nparpshell(lmax, *)&
           &, kindion(*), pshell(*), jpseudo(lmax, *)                    &
           &, d, walker, niesd, nshelltot, n_body_on, nmol, firstmol     &
           &, nelorb_c, yesfast, indt4, indt4j, nelorbpsi, contraction, indt4n&
           &, leadpsi, nmolipf, nmolshift, ii, jj, ind_store, ix, iy
    integer :: indpar_tab(*), indorb_tab(*), indshell_tab(*), indparj_tab(*)&
            &, indorbj_tab(*), indshellj_tab(*)

    logical ifdmc, ifvmc, no_dgemv

    complex*16 phs

    real*8 cclock

    real*8 ainv(*)                                                &
            &, winv(ipc*nelorb, 0:indt4, nelup + neldo), winvup(ipc*nelup, *), winvdo(max(ipc*neldo, 1), *) &
            &, epst, psip(*), rion(3, nion), winvbar(ipf*ipc*nelorbh, *)                      &
            &, winvjbar(max(ipj*nelorbjh, 1), nelup + neldo), detmat(ipc*ipf*nelorbh, *), jasmat(*)&
            &, dist(nion, nelup + neldo), winvs(ipc*(indt4 + 1)*nelorbh)&
            &, ainvs(ipc*nelup_mat*(2 + min(epscuttype, 1)*(nelup_mat - 2)))&
            &, winvsj(max(nelorbj, 1), 0:indt + ip4)    &
            &, tabpip(*), alat, ivic(3, indt, *), winvjbarn(max(ipj*nelorbjh, 1))&
            &, winvjbarszn(*)&
            &, vpot, vpotreg(2, nelup + neldo), vj(*), zeta(nion), disto, dup(*), ngivej    &
            &, psinew(ipc*nelorbh*(indt + ip4 + 1)), r(0:indt, nion), rmu(3, 0:indt, nion)  &
            &, plat(*), tmu(nelup + neldo, *), dists(nion), ratio(ipc)                   &
            &, winvj(max(nelorbj, 1), 0:indt4j, nelup + neldo), vju(*), cnorm(*)          &
            &, tcost(*), costz(nion), costz3(nion)                                   &
            &, prefactor(indt - istart + 1, *), pseudolocal(*), rcutoff(*)                &
            &, parshell(3, *), wpseudo(*), legendre(lmax - 1, *), versor(3, *)             &
            &, wintpseudo(*), pseudolocalnew, angle(18, *)                            &
            &, winvjbarsz(*), jasmatsz(*), deltaewald, timepip, timep, timewf &
            &, versoralat(3, *), Zreg2, cutreg, cost, distb, distbo     &
            &, winvfn(*), winvbarfn(*), detmat_c(*), psidetsn
    real*8 rcart(3, *), kel(3, nelup + neldo, 0:indt), keln(3, 0:indt)         &
            &, rrion(3), rpip(3, 0:1), psiln, psidetln, LBox, dist_kel(3)&
            &, dist_kelo(3), distn, rmucos(3, 0:indt, nion), rmusin(3, 0:indt, nion)&
            &, x_shift(3), oldkappa, tmp1, tmp2, derfc, psidetdmc&
            &, jastrowall_ee(nelup + neldo, nelup + neldo, 0:indt4j) &
            &, jastrowall_ei(nion, nelup + neldo)&
            &, jasnew_ee(nelup + neldo), jasnew_ei&
            &, mu_c(ipc*ipf*nelorbh, *), projm(ipc*ipf*nelorbh, *), psisn, cost0&
            &, deltaewald_single, vec(3), costreg(2)
    real*8, external :: t_lrdmc
    logical ioptpseudo, iessz, iesrandoma
    real*8 vpotsav_ee(nelup + neldo, nelup + neldo)

    !#ifdef __NOOMP
    !         integer, external ::  omp_get_max_threads
    !         integer old_threads
    !         old_threads=omp_get_max_threads()
    !#endif

    !         jel is the elctron that moves
    !         indtc   is the direction of the move

    if (iflagtab .ne. 0) return

    ifdmc = .true.
    if (itest .ne. 1 .and. itest .ne. 6 .and. itest .ne. -2 .and. itest .ne. -3) ifdmc = .false.
    ifvmc = .true.
    if (itest .ne. 2 .and. itest .ne. -5 .and. itest .ne. -7) ifvmc = .false.

    indtm(jel) = istart - 1
    nelorb5 = nelorb*(indt4 + 1)
    nel = nelup + neldo
    indt4n = indt + ip4 + 1
    nelorbpsi = nelorbh*indt4n

    if (iessz) then
        indjel = (jel - 1)*nelorbjh + 1
    else
        indjel = 1
    end if

    if (ipf .eq. 2 .and. molecular .gt. 0) then
        leadpsi = nmol*ipc
        nmolipf = nmol
        nmolshift = 0
    else
        leadpsi = nmol*ipc/ipf
        nmolipf = nmol/ipf
        nmolshift = (ipf - 1)*nmol/2
    end if

    if (.not. ifvmc) then

        !#ifdef __NOOMP
        !         call omp_set_num_threads(1)  ! serial OpenMP code
        !#endif

        !         compute the new local moves
        if (ifdmc) then
            !            call dcopy(3,rcart,1,keln(1,0),1)
            keln(1:3, 0) = rcart(1:3, 1)
            !           if(LBox.gt.0.d0) call ApplyPBC(keln(1,0),1)
        else
            do kk = 1, 3
                keln(kk, 0) = kel(kk, jel, 0) + ivic(kk, indtc, jel)
            end do

            if (iesrandoma) then
                !         randomize direction mesh
                if (enforce_detailb) then
                    !           define a random matrix that leaves unchanged the direction ivic
                    !           in principle  necessary only for indtc<  istart
                    !           for enforcing detailes balance. We adopt here this redundancy
                    !           as also for the pseudo move randomization for  indtc>=istart
                    call fillmatrix_vec(angle(1, jel), ivic(1, indtc, jel))
                else
                    call fillmatrix(angle(1, jel))
                end if
                do i = 1, istart - 1
                    call dgemv('N', 3, 3, 1.d0, angle(1, jel), 3, versoralat(1, i)              &
                            &, 1, 0.d0, ivic(1, i, jel), 1)
                end do
            end if

            !          if(LBox.gt.0.d0) call ApplyPBC(keln(1,0),1)
            ! non standard fixed node --> kinetic non local part
            do i = 1, istart - 1
                do kk = 1, 3
                    keln(kk, i) = keln(kk, 0) + ivic(kk, i, jel)
                end do
            end do
            !         if(LBox.gt.0.d0) call ApplyPBC(keln(1,1),12)
        end if

        costreg(1:2) = pot_aasunel

!       vpot=0.d0 ! recomputed by scratch updating vpotreg(1,*)
        if (LBox .le. 0.d0) then
!$omp  parallel do default(shared) reduction(+:costreg) &
!$omp  private(k,cost,Zreg2,distb,distbo)
            do k = 1, nion
!               disto = dist(k, jel)
                dist(k, jel) = dsqrt((keln(1, 0) - rion(1, k))**2 + &
                        & (keln(2, 0) - rion(2, k))**2 + (keln(3, 0) - rion(3, k))**2)
                if (dist(k, jel) .le. 1d-9) dist(k, jel) = 1d-9
!               cost = 2.d0 * zeta(k) *veps(disto) - 2.d0 * zeta(k) *veps(dist(k, jel))
                distbo = veps(dist(k, jel))
                cost = -2.d0*zeta(k)*distbo
!               vpot = vpot + cost
                costreg(1) = costreg(1) + cost
                costreg(2) = costreg(2) + cost
                Zreg2 = costz(k)*costz3(k)
                if (Zreg2*cutreg .gt. 0.d0) then
                    distb = max(dist(k, jel), cutreg)
!                   distbo = max(disto, cutreg)
!                   costreg(2) = costreg(2) + Zreg2 / dist(k, jel) - Zreg2 / distb&
                    costreg(2) = costreg(2) + Zreg2*distbo - Zreg2*veps(distb)
!                           & - Zreg2 *veps(disto) + Zreg2 *veps(distbo)
                end if
            end do
!$omp end parallel do
        else
!$omp  parallel do default(shared) reduction(+:costreg) &
!$omp  private(ii,k,cost,distreg_shift,Zreg2,dist_kel,vec,dist_shift)
            do k = 1, nion
                dist_kel(:) = keln(:, 0) - rion(:, k)
                vec(:) = car2cry(:, 1)*dist_kel(1) + car2cry(:, 2)*dist_kel(2) + car2cry(:, 3)*dist_kel(3)
                vec(1:3) = anint(vec(1:3)/cellscale(1:3))
                dist_kel(:) = dist_kel(:) - s2r(:, 1)*vec(1) - s2r(:, 2)*vec(2) - s2r(:, 3)*vec(3)
!               dist_kelo(:)=kel(:,jel,0)-rion(:,k)
!               vec(:) = car2cry(:, 1) * dist_kelo(1) + car2cry(:, 2) * dist_kelo(2) + car2cry(:, 3) * dist_kelo(3)
!               vec(1:3) = anint(vec(1:3) / cellscale(1:3))
!               dist_kelo(:) = dist_kelo(:) - s2r(:, 1) * vec(1) - s2r(:, 2) * vec(2) - s2r(:, 3) * vec(3)
#ifdef _SIMD
!$omp simd
#endif
                do ii = 2, neigh
        dist_shift(ii)=max(dsqrt((dist_kel(1)+x_neigh(ii,1))**2+(dist_kel(2)&
            &+x_neigh(ii,2))**2+(dist_kel(3)+x_neigh(ii,3))**2),1d-9)
                    distreg_shift(ii) = rep_erfc(dist_shift(ii), kappa)
                end do

                dist_shift(1) = max(dsqrt(dist_kel(1)**2 + dist_kel(2)**2 + dist_kel(3)**2), 1d-9)
                distreg_shift(1) = rep_erfc(dist_shift(1), kappa)
                dist(k, jel) = dist_shift(1) ! the identity
!              cost = 2 * zeta(k) * (sum(distrego_shift(1:neigh)-distreg_shift(1:neigh)))
                cost = -2*zeta(k)*sum(distreg_shift(1:neigh))
!               vpot = vpot + cost

                costreg(1) = costreg(1) + cost
                costreg(2) = costreg(2) + cost

                Zreg2 = costz(k)*costz3(k)
                if (Zreg2*cutreg .gt. 0.d0) then

#ifdef _SIMD
!$omp simd
#endif
                    do ii = 2, neigh
                        dist_shift(ii) = max(dist_shift(ii), cutreg)
                        distreg_shift(ii) = Zreg2*(distreg_shift(ii) - rep_erfc(dist_shift(ii), kappa))
                    end do

                    dist_shift(1) = max(dist_shift(1), cutreg)
                    distreg_shift(1) = Zreg2*(distreg_shift(1) - rep_erfc(dist_shift(1), kappa))

                    costreg(2) = costreg(2) + sum(distreg_shift(1:neigh))
                end if
            end do
!$omp end parallel do
        end if
        vpotreg(1:2, jel) = costreg(1:2)

        no_dgemv = .false.
        if (ncore .gt. 0) then
            ! non local pseudo

            vec(:) = ivic(:, indtc, jel)
            call pseudoset(jel, keln(1, 0), ivic(1, istart, jel), prefactor &
                    &, pseudolocalnew, nel, dist(1, jel), rion, wpseudo, ncore &
                    &, indt - istart + 1, indtm(jel), ioptpseudo, angle(10, jel), psip, Lbox &
                    &, 1d-9, iflagtab, vec)

            if (iflagtab .ne. 0) return
            do i = istart, indtm(jel)
                do kk = 1, 3
                    keln(kk, i) = keln(kk, 0) + ivic(kk, i, jel)
                end do
            end do

!           cost = 2.d0 * (pseudolocalnew - pseudolocal(jel))
            cost = 2.d0*pseudolocalnew
!           vpot = vpot + cost
            vpotreg(1, jel) = vpotreg(1, jel) + cost
            vpotreg(2, jel) = vpotreg(2, jel) + cost
            pseudolocal(jel) = pseudolocalnew

! no_dgemv=.false. leads only to some marginal improvement in the most  convenient
! case when indtm(i)=indt for all i. Thus the condition below.
! When using the GPU the condition is opposite. It is enough that only one is
! equal to indt and it is convenient to use dgemv because the number of independent
! threads is huge, and the elapsed time will be exhausted by the few that perform
! the calculation with indtm=indt.

#ifdef _OFFLOAD
            if (yes_ontarget .and. indt .gt. 0) then
                no_dgemv = .true.
                i = 0
                do while (no_dgemv .and. i .lt. nel)
                    i = i + 1
                    if (indtm(i) .eq. indt) no_dgemv = .false.
                end do
            elseif (indt .gt. 0) then
                i = 0
                do while (.not. no_dgemv .and. i .lt. nel)
                    i = i + 1
                    if (indtm(i) .lt. indt) no_dgemv = .true.
                end do
            end if
#else
            if (indt .gt. 0) then
                i = 0
                do while (.not. no_dgemv .and. i .lt. nel)
                    i = i + 1
                    if (indtm(i) .lt. indt) no_dgemv = .true.
                end do
            end if
#endif
        end if ! ncore.ne.0

        if (LBox .le. 0.d0) then
            cost = 0.d0
!$omp parallel do default(shared) reduction(+:cost) private(i,distbo,cost0)
            do i = 1, nel
                if (i .ne. jel) then
!                   cost0 = 2.d0 / ngivej(keln(1, 0), kel(1, i, 0), LBox)&
                    distbo = 2.d0*veps(ngivej(keln(1, 0), kel(1, i, 0), LBox))
                    cost0 = distbo - vpotsav_ee(jel, i)
!                           - 2.d0 / ngivej(kel(1, jel, 0), kel(1, i, 0), LBox)
!                  &       - 2.d0 *veps(ngivej(kel(1, jel, 0), kel(1, i, 0), LBox))
                    vpotreg(1, i) = vpotreg(1, i) + cost0/2.d0
                    vpotreg(2, i) = vpotreg(2, i) + cost0/2.d0
                    vpotsav_ee(jel, i) = distbo
                    vpotsav_ee(i, jel) = distbo
                    cost = cost + distbo
                end if
            end do
!           vpot = vpot + cost
            vpotreg(1, jel) = vpotreg(1, jel) + cost/2.d0
            vpotreg(2, jel) = vpotreg(2, jel) + cost/2.d0
            vpot = sum(vpotreg(1, 1:nel))

        else
            cost = 0.d0
!$omp parallel do default(shared) reduction(+:cost) private(i,ii,dist_kel,vec,cost0,distbo,dist_shift)
            do i = 1, nel
                if (i .ne. jel) then
                    dist_kel(:) = keln(:, 0) - kel(:, i, 0)
                    vec(:) = car2cry(:, 1)*dist_kel(1) + car2cry(:, 2)*dist_kel(2) + car2cry(:, 3)*dist_kel(3)
                    vec(1:3) = anint(vec(1:3)/cellscale(1:3))
                    dist_kel(:) = dist_kel(:) - s2r(:, 1)*vec(1) - s2r(:, 2)*vec(2) - s2r(:, 3)*vec(3)
#ifdef _SIMD
!$omp simd
#endif
                    do ii = 2, neigh
                        dist_shift(ii)=max(dsqrt((dist_kel(1)+x_neigh(ii,1))**2&
                            &+(dist_kel(2)+x_neigh(ii,2))**2+(dist_kel(3)+x_neigh(ii,3))**2),1d-9)
                        dist_shift(ii) = 2.d0*rep_erfc(dist_shift(ii), kappa)
                    end do

                    dist_shift(1) = max(dsqrt(dist_kel(1)**2 + dist_kel(2)**2 + dist_kel(3)**2), 1d-9)
                    dist_shift(1) = 2.d0*rep_erfc(dist_shift(1), kappa)

                    distbo = sum(dist_shift(1:neigh))
                    cost0 = distbo - vpotsav_ee(jel, i)
                    cost = cost + distbo
                    vpotsav_ee(jel, i) = distbo
                    vpotsav_ee(i, jel) = distbo
                else
                    cost0 = 0.d0
                end if
                vpotreg(1, i) = vpotreg(1, i) + cost0/2.d0
                vpotreg(2, i) = vpotreg(2, i) + cost0/2.d0
            end do
            vpotreg(1, jel) = vpotreg(1, jel) + cost/2.d0
            vpotreg(2, jel) = vpotreg(2, jel) + cost/2.d0

            ! ********** Updating Ewald Potential *************
            call EwaldUpdate(kel(1, jel, 0), keln, jel, deltaewald, deltaewald_single, walker)
            !       The Long range Ewald contribution to the energy is non singular and is averaged
            !       in equal footing over all the electros.
            deltaewald_single = deltaewald_single/nel
            deltaewald = (deltaewald + eself)/nel
            do i = 1, nel
                if (i .ne. jel) then
                    vpotreg(1, i) = vpotreg(1, i) + deltaewald_single
                    vpotreg(2, i) = vpotreg(2, i) + deltaewald_single
                else
                    vpotreg(1, i) = vpotreg(1, i) + deltaewald
                    vpotreg(2, i) = vpotreg(2, i) + deltaewald
                end if
            end do
            ! ********** Updating Ewald Potential *************
            vpot = sum(vpotreg(1, 1:nel))
            !         timepip=timepip+cclock()-timep
        end if

        !         definition hopping costants from the LRDMC
        do kk = 1, istart - 1
            tcost(kk) = 1.d0
            tmu(jel, kk) = t_lrdmc(keln(1, 0), nion, rion, kk, alat, ivic(1, 1, jel), plat, cellscale, iespbc)
        end do
        do kk = istart, indtm(jel)
            tcost(kk) = 1.d0
            tmu(jel, kk) = -prefactor(kk - istart + 1, jel)
        end do
        do kk = indtm(jel) + 1, indt
            tcost(kk) = 0.d0
            tmu(jel, kk) = 0.d0
        end do
        !         definition winvs e ainvs

        if (jel .le. nelup) then

            timep = cclock()

            if (indt .gt. indtm(jel)) psinew = 0.d0

            call upnewwf(indt, 0, indtm(jel), 0, nshell, ioptorb, ioccup, keln, 1, r   &
                    &, rmu, dup, zeta, rion, psip, psinew, nelorbh, nion, kion, iflagnorm           &
                    &, cnorm, LBox, rmucos, rmusin, 1d-9&
                    &, indpar_tab, indorb_tab, indshell_tab, .true.)

#ifdef _OFFLOAD
!$omp target  update to (psinew) if(yes_ontarget.or.contraction.ne.0)
#endif

            timewf = timewf + cclock() - timep

#ifdef  _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
            do i = 1, ipc*(indt4 + 1)*nelorbh
                winvs(i) = psinew(i)
            end do

            if (ifdmc) then
                if (ipc .eq. 1) then
                    if (yesfast .ne. 0) then
                        if (ipf .eq. 2) then
                            call dgemv_('T', nelorbh, nmol_ip, 1.d0, mu_c(1, firstmol), nelorbh_ip&
                           &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', 2*nelorbh, nmol_ip, 1.d0, projm(1, firstmol), nelorbh_ip&
                                &, psip, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                        else
                            call dgemv_('T', nelorbh, nmol_ip, 1.d0, projm(1, firstmol), nelorbh_ip&
                                &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', nelorbh, nmol_ip, 1.d0, mu_c(1, firstmol), nelorbh_ip&
                                &, psip, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                        end if
                    else
                        if (ipf .eq. 2) then
                            call dgemv_('N', 2*nelorbh, nelorbh, 1.d0, detmat, 2*nelorbh, &
                                &winvs, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                        else
                            call dgemv_('T', nelorbh, nelorbh, 1.d0, detmat, nelorbh, &
                                &winvs, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                        end if
                    end if
                else
                    if (yesfast .ne. 0) then
                        if (ipf .eq. 2) then
                            call zgemv_('T', nelorbh, nmol_ip, zone, mu_c(1, firstmol), nelorbh_ip&
                                &, winvs, 1, zzero, psip, 1, yes_ontarget)
                            call zgemv_('N', 2*nelorbh, nmol_ip, zone, projm(1, firstmol), nelorbh_ip&
                                &, psip, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                        else
                            call zgemv_('T', nelorbh, nmol_ip, zone, projm(1, firstmol), nelorbh_ip&
                                &, winvs, 1, zzero, psip, 1, yes_ontarget)
                            call zgemv_('N', nelorbh, nmol_ip, zone, mu_c(1, firstmol), nelorbh_ip&
                                &, psip, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                        end if
                    else
                        if (ipf .eq. 2) then
                            call zgemv_('N', 2*nelorbh, nelorbh, zone, detmat, 2*nelorbh, &
                                &winvs, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                        else
                            call zgemv_('T', nelorbh, nelorbh, zone, detmat, nelorbh, &
                                &winvs, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                        end if
                    end if
                end if
            end if

            if (contraction .ne. 0) then
                call dcopy_vec_(ipc*nelorbpsi, psinew, psip)
                if (ipc .eq. 1) then
!           if(yes_ontarget) then
                    call dgemm_('T', 'N', nmolipf, indt4n, nelorbh, 1.d0, mu_c(1, firstmol)&
                 &, ipf*nelorbh, psip, nelorbh, 0.d0, psinew, nmolipf)
!           else
!          call dgemm('T','N',nmolipf, indt4n, nelorbh, 1.d0, mu_c(1, firstmol)&
!       &, ipf * nelorbh, psip, nelorbh, 0.d0, psinew, nmolipf)
!           endif
#ifdef _OFFLOAD
!$omp target update from (psinew(1:nelorbpsi))  if(.not.yes_ontarget)
#endif
!          call dgemm_tn(nmolipf, indt4n, nelorbh, 1.d0, mu_c(1, firstmol)&
!       &, ipf * nelorbh, psip, nelorbh, 0.d0, psinew, nmolipf)
                else
!           if(yes_ontarget) then
                    call zgemm_('T', 'N', nmolipf, indt4n, nelorbh, zone, mu_c(1, firstmol)&
                 &, ipf*nelorbh, psip, nelorbh, zzero, psinew, nmolipf)
!           else
!          call zgemm('T','N',nmolipf, indt4n, nelorbh, zone, mu_c(1, firstmol)&
!       &, ipf * nelorbh, psip, nelorbh, zzero, psinew, nmolipf)
!           endif
#ifdef _OFFLOAD
!$omp target update from (psinew(1:2*nelorbpsi))  if(.not.yes_ontarget)
#endif

!          call zgemm_tn(nmolipf, indt4n, nelorbh, zone, mu_c(1, firstmol)&
!       &, ipf * nelorbh, psip, nelorbh, zzero, psinew, nmolipf)
                end if
            end if

            if (ipf .eq. 1) then

                if (ipc .eq. 1) then
                    call dgemv_('T', nmol, nelup, 1.d0, winvbarfn(nmol*nelup + 1)    &
                         &, nmol, psinew, 1, 0.d0, ainvs, 1, yes_ontarget)
                else
                    call zgemv_('T', nmol, nelup, zone, winvbarfn(ipc*nmol*nelup + 1)    &
                           &, nmol, psinew, 1, zzero, ainvs, 1, yes_ontarget)
                end if
            else
                if (ipc .eq. 1) then
                    call dgemv_('T', nmolipf, nelup_mat, 1.d0, winvbarfn    &
                                          &, nmol, psinew, 1, 0.d0, ainvs, 1, yes_ontarget)

                    ainvs(jel) = 0.d0
#ifdef _OFFLOAD
!$omp target update to(ainvs(jel))  if(yes_ontarget)
#endif
                else
                    call zgemv_('T', nmolipf, nelup_mat, zone, winvbarfn    &
                            &, nmol, psinew, 1, zzero, ainvs, 1, yes_ontarget)

                    ainvs(2*jel - 1:2*jel) = 0.d0
#ifdef _OFFLOAD
!$omp target update to (ainvs(2 * jel - 1:2 * jel)) if(yes_ontarget)
#endif
                end if
            end if

            if (nelorbj .ne. 0) then
                !          iflagnorm=-iflagnorm
                timep = cclock()
                if (indt .gt. indtm(jel)) winvsj = 0.d0

                call upnewwf(indt, 0, indtm(jel), 0, nshellj, ioptorbj, ioccj, keln, 1  &
                        &, r, rmu, vju, zeta, rion, psip, winvsj, nelorbj, nion, kionj          &
                        &, iflagnorm, cnorm(nshelltot + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
#ifdef _OFFLOAD
!$omp target  update to (winvsj)  if(yes_ontarget)
#endif
                timewf = timewf + cclock() - timep
                if (yes_sparse) then
                    winvjbarn = 0.d0
                    if (ipj .eq. 2) then
!$omp parallel do default(shared) private(i,ix,iy) reduction(+:winvjbarn)
                        do i = 1, nnozeroj
                            iy = nozeroj(i + nnozeroj)
                            ix = nozeroj(i)
                            if (iy .le. nelorbjh) winvjbarn(ix) = winvjbarn(ix) + jasmat(i)*winvsj(iy, 0)
                            if (ix .ne. iy .and. ix .le. nelorbjh) winvjbarn(iy) = winvjbarn(iy) + jasmat(i)*winvsj(ix, 0)
                        end do
                    else
!$omp parallel do default(shared) private(i,ix,iy) reduction(+:winvjbarn)
                        do i = 1, nnozeroj
                            iy = nozeroj(i + nnozeroj)
                            ix = nozeroj(i)
                            winvjbarn(ix) = winvjbarn(ix) + jasmat(i)*winvsj(iy, 0)
                            if (ix .ne. iy) winvjbarn(iy) = winvjbarn(iy) + jasmat(i)*winvsj(ix, 0)
                        end do
                    end if
#ifdef _OFFLOAD
!$omp target update to (winvjbarn) if(yes_ontarget)
#endif
                else

                    if (contractionj .ne. 0) then
                        call dgemv_('T', nelorbjh, nelorbj_c, 1.d0, muj_c &
               &, nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                        call dgemv_('N', ipj*nelorbj_c, nelorbj_c, 1.d0, jasmat_c &
               &, ipj*nelorbj_c, winvjbarn, 1, 0.d0, psip, 1, yes_ontarget)
                        call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
               &, nelorbjh, psip, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                        if (ipj .eq. 2)  &
             &        call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
             &, nelorbjh, psip(nelorbj_c + 1), 1, 0.d0, winvjbarn(nelorbjh + 1), 1, yes_ontarget)
                    else
                        call dgemv_('T', nelorbjh, ipj*nelorbjh, 1.d0, jasmat &
               &, ipj*nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                    end if
                end if

                if (iessz) call dgemv('N', nelorbjh, nelorbjh, 1.d0, jasmatsz        &
                        &, nelorbjh, winvsj, 1, 0.d0, winvjbarszn, 1)
            end if

            psiln = psiln - psidetdmc
            call upinvhop_fnf(nelorb, jel, indt, ainv, winv           &
                    &, winvup, winvdo, psidetdmc, psisn, epst, psip, nelup, neldo, ainvs, winvs  &
                    &, psinew, leadpsi, detmat, indtm, nelorbh, mu_c, nelorb_c&
                    &, firstmol, nmol, yesfast, winvbarfn, winvfn, detmat_c, contraction, indt4, no_dgemv)

            psiln = psiln + psidetdmc

            ! jel.le.nelup
        else
            timep = cclock()
            if (indt .gt. indtm(jel)) psinew = 0.d0

            call upnewwf(indt, 0, indtm(jel), 0, nshelldo, ioptorb, ioccdo, keln, 1   &
                    &, r, rmu, dup, zeta, rion, psip, psinew, nelorbh, nion, kion, iflagnorm       &
                    &, cnorm, LBox, rmucos, rmusin, 1d-9&
                    &, indpar_tab, indorb_tab, indshell_tab, .false.)

            timewf = timewf + cclock() - timep
!           if(indt.gt.indtm(jel)) then
!      call dscalzero__(ipc*nelorbh*(indt-indtm(jel)),0.d0&
!                   &, psinew(ipc * nelorbh * (indtm(jel) + 1) + 1), 1)
!           endif
#ifdef _OFFLOAD
!$omp target update to (psinew) if(yes_ontarget.or.contraction.ne.0)
#endif

#ifdef  _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
            do i = 1, ipc*(indt4 + 1)*nelorbh
                winvs(i) = psinew(i)
            end do

            !       winvs(1:ipc*(indt4+1)*nelorbh)=psinew(1:ipc*(indt4+1)*nelorbh)

            if (ifdmc) then
                ! CCC to be done
                if (ipc .eq. 1) then
                    if (yesfast .ne. 0) then
                        if (ipf .eq. 2) then
                            call dgemv_('T', nelorbh, nmol_ip, 1.d0, mu_c(1 + nelorbh, firstmol), nelorbh_ip&
                                    &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', 2*nelorbh, nmol_ip, 1.d0, projm(1, firstmol), nelorbh_ip&
                                    &, psip, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                        else
                            call dgemv_('T', nelorbh, nmol_ip, 1.d0, mu_c(1, firstmol), nelorbh_ip&
                                    &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', nelorbh, nmol_ip, 1.d0, projm(1, firstmol), nelorbh_ip&
                                    &, psip, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                        end if
                    else
                        call dgemv_('N', ipf*nelorbh, nelorbh, 1.d0, detmat(1, 1 + (ipf - 1)*nelorbh)&
                                &, ipf*nelorbh, winvs, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                    end if
                else
                    if (yesfast .ne. 0) then
                        if (ipf .eq. 2) then
                            call zgemv_('T', nelorbh, nmol_ip, zone, mu_c(1 + 2*nelorbh, firstmol), nelorbh_ip&
                                    &, winvs, 1, zzero, psip, 1, yes_ontarget)
                            call zgemv_('N', 2*nelorbh, nmol_ip, zone, projm(1, firstmol), nelorbh_ip&
                                    &, psip, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                        else
                            call zgemv_('T', nelorbh, nmol_ip, zone, mu_c(1, firstmol), nelorbh_ip&
                                    &, winvs, 1, zzero, psip, 1, yes_ontarget)
                            call zgemv_('N', nelorbh, nmol_ip, zone, projm(1, firstmol), nelorbh_ip&
                                    &, psip, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                        end if
                    else
                        !                if(ipf.eq.2) then
                        call zgemv_('N', ipf*nelorbh, nelorbh, zone, detmat(1, (ipf - 1) + nelorbh + 1), ipf*nelorbh, &
                                &winvs, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                        !                else
                        !                   call zgemv('N',nelorbh,nelorbh,zone,detmat,nelorbh,&
                        !                        &winvs,1,zzero,winvbar(1,jel),1)
                        !                end if
                    end if
                end if
            end if

            if (contraction .ne. 0) then
                !       contract psinew
                !          psip(1:ipc*nelorbpsi)=psinew(1:ipc*nelorbpsi)
                call dcopy_vec_(ipc*nelorbpsi, psinew, psip)
                if (ipc .eq. 1) then
!                 if(yes_ontarget) then
                    call dgemm_('T', 'N', nmolipf, indt4n, nelorbh, 1.d0, mu_c(1 + (ipf - 1)*nelorbh, firstmol + nmolshift)&
                                         &, ipf*nelorbh, psip, nelorbh, 0.d0, psinew, nmolipf)
!     else
!      call dgemm('T','N',nmolipf, indt4n, nelorbh, 1.d0, mu_c(1 + (ipf - 1) * nelorbh, firstmol + nmolshift)&
!                           &, ipf * nelorbh, psip, nelorbh, 0.d0, psinew, nmolipf)
!     endif
#ifdef _OFFLOAD
!$omp target update from (psinew(1:nelorbpsi))  if(.not.yes_ontarget)
#endif

!                    call dgemm_tn(nmolipf, indt4n, nelorbh, 1.d0, mu_c(1 + (ipf - 1) * nelorbh, firstmol + nmolshift)&
!                            &, ipf * nelorbh, psip, nelorbh, 0.d0, psinew, nmolipf)

                else
!         if(yes_ontarget) then
                    call zgemm_('T', 'N', nmolipf, indt4n, nelorbh, zone, mu_c(1 + (ipf - 1)*2*nelorbh, firstmol + nmolshift)&
                      &, ipf*nelorbh, psip, nelorbh, zzero, psinew, nmolipf)
!         else
!          call zgemm('T','N',nmolipf, indt4n, nelorbh, zone, mu_c(1 + (ipf - 1) * 2 * nelorbh, firstmol + nmolshift)&
!            &, ipf * nelorbh, psip, nelorbh, zzero, psinew, nmolipf)
!         endif
#ifdef _OFFLOAD
!$omp target update from (psinew(1:2*nelorbpsi))  if(.not.yes_ontarget)
#endif

!                   call zgemm_tn(nmolipf, indt4n, nelorbh, zone, mu_c(1 + (ipf - 1) * 2 * nelorbh, firstmol + nmolshift)&
!                           &, ipf * nelorbh, psip, nelorbh, zzero, psinew, nmolipf)
                end if
            end if
            if (ipf .eq. 2) then
                if (ipc .eq. 1) then
                    call dgemv_('T', nmolipf, nelup_mat, 1.d0, winvbarfn(nmolshift + 1)&
                         &, nmol, psinew, 1, 0.d0, ainvs, 1, yes_ontarget)

                    ainvs(jel) = 0.d0
#ifdef _OFFLOAD
!$omp target update to (ainvs(jel)) if(yes_ontarget)
#endif
                else
                    call zgemv_('T', nmolipf, nelup_mat, zone, winvbarfn(nmolshift*2 + 1)&
                            &, nmol, psinew, 1, zzero, ainvs, 1, yes_ontarget)

                    ainvs(2*jel - 1:2*jel) = 0.d0
#ifdef _OFFLOAD
!$omp target update to(ainvs(2 * jel - 1:2 * jel))  if(yes_ontarget)
#endif
                end if
            else
                if (ipc .eq. 1) then
                    call dgemv_('T', nmol, nelup_mat, 1.d0, winvbarfn&
                            &, nmol, psinew, 1, 0.d0, ainvs, 1, yes_ontarget)
                else
                    call zgemv_('T', nmol, nelup_mat, zone, winvbarfn&
                            &, nmol, psinew, 1, zzero, ainvs, 1, yes_ontarget)
                end if
            end if

            if (nelorbj .ne. 0) then
                !          iflagnorm=-iflagnorm
                timep = cclock()
                if (indt .gt. indtm(jel)) winvsj = 0.d0
                call upnewwf(indt, 0, indtm(jel), 0, nshellj, ioptorbj, ioccj, keln, 1    &
                        &, r, rmu, vju, zeta, rion, psip, winvsj, nelorbj, nion, kionj          &
                        &, iflagnorm, cnorm(nshelltot + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .false.)
                timewf = timewf + cclock() - timep
#ifdef _OFFLOAD
!$omp target  update to (winvsj) if(yes_ontarget)
#endif

!               if(indt.gt.indtm(jel)) then
!    call dscalzero__(nelorbj*(indt-indtm(jel)),0.d0,winvsj(1,indtm(jel)+1),1)
!               endif
                !         update winvjbar
                if (yes_sparse) then
                    winvjbarn = 0.d0
                    if (ipj .eq. 2) then
!$omp parallel do default(shared) private(i,ix,iy) reduction(+:winvjbarn)
                        do i = 1, nnozeroj
                            iy = nozeroj(i + nnozeroj)
                            ix = nozeroj(i)
                            if (iy .gt. nelorbjh) winvjbarn(ix) = winvjbarn(ix) + jasmat(i)*winvsj(iy - nelorbjh, 0)
                            if (ix .ne. iy .and. ix .gt. nelorbjh) winvjbarn(iy)&
                               & = winvjbarn(iy) + jasmat(i)*winvsj(ix - nelorbjh, 0)
                        end do
                    else
!$omp parallel do default(shared) private(i,ix,iy) reduction(+:winvjbarn)
                        do i = 1, nnozeroj
                            iy = nozeroj(i + nnozeroj)
                            ix = nozeroj(i)
                            winvjbarn(ix) = winvjbarn(ix) + jasmat(i)*winvsj(iy, 0)
                            if (ix .ne. iy) winvjbarn(iy) = winvjbarn(iy) + jasmat(i)*winvsj(ix, 0)
                        end do
                    end if
#ifdef _OFFLOAD
!$omp target update to (winvjbarn) if(yes_ontarget)
#endif
                else
                    if (contractionj .ne. 0) then
                        call dgemv_('T', nelorbjh, nelorbj_c, 1.d0, muj_c &
             &, nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                        if (ipj .eq. 1) then
                            call dgemv_('N', nelorbj_c, nelorbj_c, 1.d0, jasmat_c &
                 &, nelorbj_c, winvjbarn, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
                 &, nelorbjh, psip, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                        else
                            call dgemv_('T', nelorbj_c, 2*nelorbj_c, 1.d0, jasmat_c(1 + nelorbj_c) &
                        &, 2*nelorbj_c, winvjbarn, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
                 &, nelorbjh, psip, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                            call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
                 &, nelorbjh, psip(1 + nelorbj_c), 1, 0.d0, winvjbarn(1 + nelorbjh), 1, yes_ontarget)
                        end if
                    else
                        if (ipj == 1) then
                            call dgemv_('T', nelorbjh, nelorbjh, 1.d0, jasmat&
                                &, nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                        else
                            call dgemv_('T', nelorbjh, 2*nelorbjh, 1.d0&
                                &, jasmat(1 + nelorbjh), 2*nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                        end if
                    end if
                end if
                if (iessz) call dgemv('N', nelorbjh, nelorbjh, 1.d0, jasmatsz&
                        &, nelorbjh, winvsj, 1, 0.d0, winvjbarszn, 1)
            end if

            psiln = psiln - psidetdmc

            call upinvhop_fnf(nelorb, jel, indt, ainv                &
                    &, winv, winvup, winvdo, psidetdmc, psisn, epst, psip, nelup, neldo, ainvs&
                    &, winvs, psinew, leadpsi, detmat, indtm, nelorbh, mu_c, nelorb_c&
                    &, firstmol, nmol, yesfast, winvbarfn, winvfn, detmat_c, contraction, indt4, no_dgemv)

            psiln = psiln + psidetdmc
            ! jel.le.nelup

        end if

        timep = cclock()
        ind_store = (3*indt + 9)*nel + 1
        call uptabpip(indt, nel, nelup, jel, kel, keln, tabpip, tcost, iesd, vj, winvsj&
       &, winvj, nelorbj, psip, psip(ind_store), winvjbar, winvjbarsz, winvjbarn&
       &, winvjbarszn, rion, nion, costz, costz3, iessz, LBox, n_body_on, jastrowall_ee&
       &, jastrowall_ei(1, jel), rmu, niesd, indtm, nelorbjh, psiln, no_dgemv)

        if (nelorbj .ne. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
            do i = 1, nelorbjh
                winvj(i, 0:indt4j, jel) = winvsj(i, 0:indt4j)
            end do
        end if
        do i = 0, indt
            kel(1:3, jel, i) = keln(1:3, i)
        end do
        timepip = timepip + cclock() - timep

    else
        !        VMC case
        !        Update in any case winvbar and winvjbar if nelorbj>0
        ! Update agp
        if (epscuttype .eq. 2) then
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
        if (ipc .eq. 1) then
            if (jel .gt. nelup) then
                if (yesfast .ne. 0) then
                    if (ipf .eq. 2) then
                        call dgemv_('T', nelorbh, nmol_ip, 1.d0, mu_c(nelorbh + 1, firstmol), nelorbh_ip&
                                &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                    else
                        call dgemv_('T', nelorbh, nmol_ip, 1.d0, mu_c(1, firstmol), nelorbh_ip&
                                &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                    end if

                    call dgemv_('N', ipf*nelorbh, nmol_ip, 1.d0, projm(1, firstmol), nelorbh_ip&
                            &, psip, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)

                else
                    if (ipf .eq. 2) then
                        call dgemv_('N', 2*nelorbh, nelorbh, 1.d0, detmat(1, nelorbh + 1)&
                                &, 2*nelorbh, winvs, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                    else
                        call dgemv_('N', nelorbh, nelorbh, 1.d0, detmat, nelorbh, &
                                &winvs, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                    end if
                end if
            else ! jel > nelup
                if (yesfast .ne. 0) then
                    if (ipf .eq. 2) then
                        call dgemv_('T', nelorbh, nmol_ip, 1.d0, mu_c(1, firstmol), nelorbh_ip&
                                &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                        call dgemv_('N', 2*nelorbh, nmol_ip, 1.d0, projm(1, firstmol), nelorbh_ip&
                                &, psip, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                    else
                        call dgemv_('T', nelorbh, nmol_ip, 1.d0, projm(1, firstmol), nelorbh_ip&
                                &, winvs, 1, 0.d0, psip, 1, yes_ontarget)
                        call dgemv_('N', nelorbh, nmol_ip, 1.d0, mu_c(1, firstmol), nelorbh_ip&
                                &, psip, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                    end if
                else
                    if (ipf .eq. 2) then
                        call dgemv_('N', 2*nelorbh, nelorbh, 1.d0, detmat, 2*nelorbh, &
                                &winvs, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                    else
                        call dgemv_('T', nelorbh, nelorbh, 1.d0, detmat, nelorbh, &
                                &winvs, 1, 0.d0, winvbar(1, jel), 1, yes_ontarget)
                    end if
                end if
            end if

        else ! if ipc=1 , below the complex case

            if (jel .gt. nelup) then
                if (yesfast .ne. 0) then
                    if (ipf .eq. 2) then
                        call zgemv_('T', nelorbh, nmol_ip, zone, mu_c(1 + 2*nelorbh, firstmol), nelorbh_ip&
                                &, winvs, 1, zzero, psip, 1, yes_ontarget)
                    else
                        call zgemv_('T', nelorbh, nmol_ip, zone, mu_c(1, firstmol), nelorbh_ip&
                                &, winvs, 1, zzero, psip, 1, yes_ontarget)
                    end if
                    call zgemv_('N', ipf*nelorbh, nmol_ip, zone, projm(1, firstmol), nelorbh_ip&
                            &, psip, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                else
                    if (ipf .eq. 2) then
                        call zgemv_('N', 2*nelorbh, nelorbh, zone, detmat(1, nelorbh + 1), 2*nelorbh, &
                                &winvs, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                    else
                        call zgemv_('N', nelorbh, nelorbh, zone, detmat, nelorbh, &
                                &winvs, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                    end if
                end if
            else
                if (yesfast .ne. 0) then
                    if (ipf .eq. 2) then
                        call zgemv_('T', nelorbh, nmol_ip, zone, mu_c(1, firstmol), nelorbh_ip&
                                &, winvs, 1, zzero, psip, 1, yes_ontarget)
                        call zgemv_('N', 2*nelorbh, nmol_ip, zone, projm(1, firstmol), nelorbh_ip&
                                &, psip, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                    else
                        call zgemv_('T', nelorbh, nmol_ip, zone, projm(1, firstmol), nelorbh_ip&
                                &, winvs, 1, zzero, psip, 1, yes_ontarget)
                        call zgemv_('N', nelorbh, nmol_ip, zone, mu_c(1, firstmol), nelorbh_ip&
                                &, psip, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                    end if
                else
                    if (ipf .eq. 2) then
                        call zgemv_('N', 2*nelorbh, nelorbh, zone, detmat, 2*nelorbh, &
                                &winvs, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                    else
                        call zgemv_('T', nelorbh, nelorbh, zone, detmat, nelorbh, &
                                &winvs, 1, zzero, winvbar(1, jel), 1, yes_ontarget)
                    end if
                end if
            end if

        end if

        !         update winvjbar
        if (nelorbj .ne. 0) then
            if (yes_sparse) then
                winvjbarn = 0.d0
                if (ipj .eq. 2) then
                    if (jel .le. nelup) then
!$omp parallel do default(shared) private(i,ix,iy) reduction(+:winvjbarn)
                        do i = 1, nnozeroj
                            iy = nozeroj(i + nnozeroj)
                            ix = nozeroj(i)
                            if (iy .le. nelorbjh) winvjbarn(ix) = winvjbarn(ix) + jasmat(i)*winvsj(iy, 0)
                            if (ix .ne. iy .and. ix .le. nelorbjh) winvjbarn(iy) = winvjbarn(iy) + jasmat(i)*winvsj(ix, 0)
                        end do
                    else
!$omp parallel do default(shared) private(i,ix,iy) reduction(+:winvjbarn)
                        do i = 1, nnozeroj
                            iy = nozeroj(i + nnozeroj)
                            ix = nozeroj(i)
                            if (iy .gt. nelorbjh) winvjbarn(ix) = winvjbarn(ix) + jasmat(i)*winvsj(iy - nelorbjh, 0)
                            if (ix .ne. iy .and. ix .gt. nelorbjh) winvjbarn(iy)&
                               & = winvjbarn(iy) + jasmat(i)*winvsj(ix - nelorbjh, 0)
                        end do
                    end if
                else
!$omp parallel do default(shared) private(i,ix,iy) reduction(+:winvjbarn)
                    do i = 1, nnozeroj
                        iy = nozeroj(i + nnozeroj)
                        ix = nozeroj(i)
                        winvjbarn(ix) = winvjbarn(ix) + jasmat(i)*winvsj(iy, 0)
                        if (ix .ne. iy) winvjbarn(iy) = winvjbarn(iy) + jasmat(i)*winvsj(ix, 0)
                    end do
                end if
                call dcopy(ipj*nelorbjh, winvjbarn, 1, winvjbar(1, jel), 1)
#ifdef _OFFLOAD
!$omp target update to (winvjbar(1:ipj*nelorbjh,jel:jel)) if(yes_ontarget)
#endif
            else
                if (ipj .eq. 2) then
                    if (jel .le. nelup) then
                        if (contractionj .ne. 0) then
                            call dgemv_('T', nelorbjh, nelorbj_c, 1.d0, muj_c &
                   &, nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                            call dgemv_('N', 2*nelorbj_c, nelorbj_c, 1.d0, jasmat_c &
                   &, 2*nelorbj_c, winvjbarn, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
                   &, nelorbjh, psip, 1, 0.d0, winvjbar(1, jel), 1, yes_ontarget)
                            call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
                    &, nelorbjh, psip(nelorbj_c + 1), 1, 0.d0, winvjbar(nelorbjh + 1, jel), 1, yes_ontarget)
                        else
                            call dgemv_('T', nelorbjh, 2*nelorbjh, 1.d0, jasmat&
                       &, 2*nelorbjh, winvsj, 1, 0.d0, winvjbar(1, jel), 1, yes_ontarget)
                        end if
                    else
                        if (contractionj .ne. 0) then
                            call dgemv_('T', nelorbjh, nelorbj_c, 1.d0, muj_c &
                 &, nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                            call dgemv_('T', nelorbj_c, 2*nelorbj_c, 1.d0, jasmat_c(1 + nelorbj_c) &
                        &, 2*nelorbj_c, winvjbarn, 1, 0.d0, psip, 1, yes_ontarget)
                            call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
                 &, nelorbjh, psip, 1, 0.d0, winvjbar(1, jel), 1, yes_ontarget)
                            call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
                &, nelorbjh, psip(1 + nelorbj_c), 1, 0.d0, winvjbar(1 + nelorbjh, jel), 1, yes_ontarget)
                        else
                            call dgemv_('T', nelorbjh, 2*nelorbjh, 1.d0, jasmat(1 + nelorbjh)&
                           &, 2*nelorbjh, winvsj, 1, 0.d0, winvjbar(1, jel), 1, yes_ontarget)
                        end if
                    end if
                else
                    if (contractionj .ne. 0) then
                        call dgemv_('T', nelorbjh, nelorbj_c, 1.d0, muj_c &
             &, nelorbjh, winvsj, 1, 0.d0, winvjbarn, 1, yes_ontarget)
                        call dgemv_('N', nelorbj_c, nelorbj_c, 1.d0, jasmat_c &
             &, nelorbj_c, winvjbarn, 1, 0.d0, psip, 1, yes_ontarget)
                        call dgemv_('N', nelorbjh, nelorbj_c, 1.d0, muj_c &
             &, nelorbjh, psip, 1, 0.d0, winvjbar(1, jel), 1, yes_ontarget)
                    else
                        call dgemv_('T', nelorbjh, nelorbjh, 1.d0, jasmat&
                       &, nelorbjh, winvsj, 1, 0.d0, winvjbar(1, jel), 1, yes_ontarget)
                    end if
                end if
            end if
            if (iessz) call dgemv('N', nelorbjh, nelorbjh, 1.d0, jasmatsz &
                    &, nelorbjh, winvsj, 1, 0.d0, winvjbarsz(indjel), 1)
        end if
        if (ipc .eq. 1) then
            psidetsn = 1.d0
        else
            psidetsn = 0.d0
        end if
        psidetln = 0.d0

        if (epscuttype .eq. 0) then
            call upinvhop(nelorbh, nelorb5, jel, ainv, winv, psidetln, &
                          psidetsn, epst, psip, nelup, neldo, ainvs, winvs, .true.)
        else
            !    do i=1,ipc*nelorb
            !    winv(i,0,jel)=winvs(i)
            !    enddo
            call dcopy_vec(ipc*nelorb, winvs, winv(1, 0, jel))
            call dcopy_vec(ipc*nelup_mat*nelup_mat, ainvs, ainv)
        end if

        if (ipc .eq. 1) then
            psisn = psisn*sign(1.d0, ratio(1))
            psiln = psiln + log(abs(ratio(1)))
        else
            phs = log(dcmplx(ratio(1), ratio(2)))
            psiln = psiln + real(phs)
            psisn = psisn + aimag(phs)
        end if

        do i = 1, nion
            dist(i, jel) = dists(i)
        end do

        jastrowall_ei(1, jel) = jasnew_ei

        do i = 1, nel
            jastrowall_ee(i, jel, 0) = jasnew_ee(i)
            !    Use the symmetry
            jastrowall_ee(jel, i, 0) = jasnew_ee(i)
        end do
        ! Use the symmetry
        if (nelorbj .ne. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
            do j = 1, nelorbjh
                winvj(j, 0, jel) = winvsj(j, 0)
            end do
        end if
        kel(1:3, jel, 0) = rcart(1:3, 1)

    end if
    !
#ifdef DEBUG
    write (6, *)
    write (6, *) ' jastrowall_ei after uptabtot in DEBUG=', sum(abs(jastrowall_ei(1, :)))
#endif
    !stop
    return
end subroutine uptabtot
