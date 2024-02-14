!TL off
subroutine compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
        &, ipsip, wconfn, psisn, iesdr, vj, dd                       &
        &, zeta, rion, dist, ioccup, ioccdo, ioptorb                             &
        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vju, nelorbj   &
        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflag, ncore, lmax        &
        &, nintpseudo, prefactor, rcutoff, parshell                            &
        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
        &, jpseudo, pseudolocal, istart, costz, costz3&
        &, angle, indtm, LBox, rmucos, rmusin, oldkappa, vpotreg, cutreg           &
        &, psidetln, walker, nelorbh, nelorbjh, niesd&
        &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmol, yesfast, eloc, logpsi&
        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscalen&
        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
    !  NB here yesfast and yesfastj DO NOT necessarily coincide with global
    !  variables contraction and contractionj in the main, respectively
    !  (>0 if contracted orbitals are used).


    !      nelorbjmax=max(nelorbj,1)
    !      neldomax=max(neldo,1)
    !      indtmax=max(indt,1)
    !      nshelljmax=max(nshellj,1)
#ifdef _NVTX
    use nvtx
#endif
    use Ewald
    use Cell
    use Constants, only : yes_ontarget,ip4, zzero, zone, nbdgetri, ipc, ipj, ipf, pi, two_pi
    use allio, only : yes_complex, jastrowall_ei, jastrowall_ee, &
            &winvfn, nel2wtfn, winvbarfn, nel2barfn, itestrfn, psiln, contraction, &
            &test_aad, iespbc, pseudologic, iesrandoml, rank, versoralat, singdet, &
            &epscuttype, iflagnorm, agp, agpn, itest, pointvj, nelup_mat, nel_mat, ndiff&
            &, n_body_on, gamma, lrdmc_der, lrdmc_nonodes, molecular, npar_eagp, eagp_pfaff,timings,cutweight,true_wagner,npow,membig,membigcpu,count_zerowf,count_allwf,yes_crystalj,nelsquare,vpotsav_ee,yes_sparse,nnozeroj,nelorbjh2,nozeroj, norm_metric
    implicit none

    integer nelup, nelused, neldo, nel, i, j, k, ierr, j1, j2, info, ispin, indt4, indt4j&
            &, iesd, iesdr, iesdr2iesd, nelorb, nelorbj, nion         &
            &, kk, jj, nshell, nshelldo, indt, nelorb5, nelorbj5    &
            &, iflag, ncore, lmax, istart, nintpseudo&
            &, psidetsn, walker, firstmmu, firstmdet&
            &, nelorbh, nelorbjh, niesd, signold, signnew, nelorb_c, firstmol, nmol, yesfast&
            &, ind1, ind2, ind3, ind4, ind5, ind6, isdistp, ip5, sizework1&
            &, sizewsz, nelorbjmax, neldomax, indtmax, nshellj, nshelljmax, nelorbj_c&
            &, nshellh, nshelljh, nshelldoh, nmolipf, nmolshift, nelup1, zerop1&
            &, contractionj, nelorbjc, nelorbjcp, wf_sign, nelorbju, indexi, indexj, istartu&
            &, indpfaff,ix,iy

    integer ioccup(*), ioccdo(*), ioptorb(*), ioptorbj(*)&
            &, kion(*), kionj(*)&
            &, ioccj(*), nparpshell(lmax, ncore), kindion(ncore + 1)&
            &, indtm(nelup + neldo), pshell(ncore), jpseudo(lmax, ncore), ipsip(*)
    real*8 tabpip(nelup + neldo, indt + 4), tmu(nelup + neldo, max(indt, 1))        &
            &, psip(*), winv(ipc * nelorb, 0:indt4, nelup + neldo)          &
            &, winvj(nelorbjmax, 0:indt4j, nelup + neldo)                       &
            &, winvup(ipc * nelup, indt + 4), winvdo(max(ipc * neldo, 1), indt + 4)&
            &, ainv(ipc * nelup_mat, nelup_mat) &
            &, ainvup(ipc * nelup, nelorbh), ainvdo(max(ipc * neldo, 1), nelorbh), ukwald            &
            &, wconfn, vj(*), vju(*), dist_kel(3)                            &
            &, jastrow1, dd(*), zeta(nion), grad1(3), grad2, rion(3, nion)         &
            &, iond(nion, nion)                                                  &
            &, dist(nion, nelup + neldo), r(0:indt,nion), rmu(3,0:indt, nion)       &
            &, ivic(3, indtmax, nelup + neldo), alat, plat(3), vold         &
            &, vpot, vpotreg(2, *)                                         &
            &, rc(3), winvbar(ipf * ipc * nelorbh, nel_mat), detmat(ipc * ipf * nelorbh, *)&
            &, winvjbar(max(int(ipj * nelorbjh,8)*(nelup+neldo),1)), jasmat(*), cnorm(*)   &
            &, prefactor(indt - istart + 1, nelup + neldo), pseudolocal(nelup + neldo)&
            &, rcutoff(ncore)   &
            &, parshell(3, *), wpseudo(2 * lmax), legendre(lmax - 1, nintpseudo)&
            &, versor(3, nintpseudo)&
            &, wintpseudo(nintpseudo), vpseudolocal, costz(*), costz3(*), costz0 &
            &, angle(18, *), logpsi(ipc), eloc(ipc), psidetln, psidetlnt(ipc), jasmatsz(*)                  &
            &, winvjbarsz(*)                           &
            &, iond_cart(3, nion, nion), cutreg&
            &, mu_c(*), cellscalen(12), celldmo(6), detmat_c(*), projm(*)&
            &, muj_c(*), jasmat_c(*), jasmatsz_c(*), psisn, cost, psidetln_old&
            &, jastrowee_old, jastrowei_old, winvjbar_old, winvjbarsz_old, ainv_old&
            &, winvfn_old, winvbar_old, winvbarfn_old, winvup_old, winvdo_old&
            &, winvuplap_old, winvdolap_old, logpsi_old, agp_old, vpotreg_old, tmu_old, kel_old, tabpip_old(2),timep

    real*8 kel(3, nelup + neldo, 0:indt), kelind(3, nelup + neldo)
    logical iessz, iesiesd, yesprojm, lrdmc_deru
    !************* PERIODIC WORLD ************************+

    real(8) LBox, rmucos(3,0:indt, nion), rmusin(3,0:indt, nion), oldkappa
    integer indpar_tab(*), indorb_tab(*), indshell_tab(*), indparj_tab(*)&
            &, indorbj_tab(*), indshellj_tab(*)

    integer :: ind_winvfn, ind_winvbarfn,ind_vpotsav,firstmolr
    logical :: yesfn, save_jastrowall, yesupwf, yesupwfj
    real(8) cosalpha, cosbeta, cosgamma, singamma
    real*8,  external:: cclock
    real*8,dimension(:,:,:), allocatable:: winv_big,winvj_big

    ! real*8 winvup_sav(ipc*nelup,indt+4),winvdo_sav(ipc*neldo,indt+4)

    !  real*8, allocatable:: ainv_sav(:,:)
    if(lrdmc_der.and..not.lrdmc_nonodes) then
        lrdmc_deru = .true.
    else
        lrdmc_deru = .false.
    endif

    if(neldo.eq.0) then
        nelup1 = nelup
        zerop1 = 0
    else
        zerop1 = 1
        nelup1 = nelup + 1
    endif
    if(ipf.eq.2.and.molecular.gt.0) then
        nmolipf = nmol
        nmolshift = ipc * nelorbh
    else
        nmolipf = nmol / ipf
        nmolshift = (ipf - 1) * (nmol + 1) * ipc * nelorbh
    endif

    if(lrdmc_deru) then
        istartu = 1
    else
        istartu = istart
    endif

    nshelldoh = nshelldo
    nshellh = nshell
    nshelljh = nshellj
    ! indices for FN tables
    ind_winvfn = (walker - 1) * nel2wtfn + 1
    ind_winvbarfn = (walker - 1) * nel2barfn + 1
    ! compute everything from scratch

    nelorbju = ipj * nelorbjh

    indpfaff = ipc * nelup_mat * nelup_mat + 1

    ! if(rank.eq.0)  write(6,*) ' Input kel =',sum(abs(kel(:,:,0))),sum(abs(angle(:,1:nelup+neldo)))


    ! keep track of sign/phase
    if(ipc.eq.2) then
        if(abs(psisn - TWO_PI * nint(psisn / TWO_PI)).lt.Pi / 2.d0) then
            signold = 1
        else
            signold = -1
        endif
    else
        signold = nint(psisn)
        !    if(ipf.eq.2) signold=1 !  determinant is the square of the pfaffian
    endif
    if(singdet(walker)) signold = 0

    !
    if(iflag.ne.0) return
    nel = nelup + neldo
    ip5 = 5
    nelorbjc = nel * nelorbj_c
    nelorbjcp = nelorbjc + 1
    !
    if(iessz) then
        sizewsz = nelorbjh * nel
    else
        sizewsz = 1
    endif
    if(ipj.eq.2) then
        nelused = nelup
    else
        nelused = nel
    endif

#if defined (_OFFLOAD) && defined (_DEBUG)
!$omp target update from (winvbar,ainv,winvjbar,agp(:,:,walker:walker)&
!$omp &,agpn,winv,winvj,winvfn,winvbarfn) if(yes_ontarget)
#endif

#ifdef _DEBUG
  if(epscuttype.eq.2) agp_old=sum(abs(agp(:,:,walker)))
  vpotreg_old=sum(abs(vpotreg(:,1:nel)))
  kel_old=sum(kel(:,:,:))
  tabpip_old(1)=sum(abs(tabpip(:,1:indt)))
  tabpip_old(2)=sum(abs(tabpip(:,indt+1:indt+4)))
!allocate(ainv_sav(ipc*nel_mat,nel_mat))
!ainv_sav=ainv

  if(indt.gt.0) tmu_old=sum(abs(tmu(:,:)))

  ainv_old=sum(abs(ainv(:,:)))
  winvup_old=sum(abs(winvup(:,:)))
!  winvup_sav=winvup
  winvuplap_old=sum(abs(winvup(:,indt+1:indt+4)))
  winvdo_old=sum(abs(winvdo(:,:)))
!  winvdo_sav=winvdo
  winvdolap_old=sum(abs(winvdo(:,indt+1:indt+4)))
  jastrowee_old=sum(abs(jastrowall_ee(:,:,0,walker)))
  jastrowei_old=sum(abs(jastrowall_ei(1,:,walker)))
  winvjbar_old=sum(abs(winvjbar(1:int(ipj*nelorbjh,8)*nel)))
  if(iessz) winvjbarsz_old=sum(abs(winvjbarsz(1:nelorbjh*nel)))
  winvfn_old=sum(abs(winvfn(ind_winvfn:ind_winvfn+nel2wtfn-1)))
  winvbarfn_old=sum(abs(winvbarfn(ind_winvbarfn:ind_winvbarfn+nel2barfn-1)))
  winvbar_old=sum(abs(winvbar(:,:)))
  logpsi_old=psiln(walker)
#endif
    !

    ! independent input kelind rion
    kel(:, :, 0) = kelind(:, :)
    if(iespbc) then
        celldmo(1:6) = celldm(1:6)
        !   if(yes_tilted) then
        ! input  cellscalen
        call dcopy(9, cellscalen(4), 1, s2r, 1)
        givens2r = .true.
        !   else
        !   givens2r=.false.
        !   celldm(1)=cellscalen(1)
        !   celldm(2)=cellscalen(2)/cellscalen(1)
        !   celldm(3)=cellscalen(3)/cellscalen(1)
        !   endif
        call InitCell(nion, nel, yes_complex)
        if(test_AAD) then
            call Ggen_samegrid
        else
            call Ggen
        endif
    endif
    iesd = iesdr2iesd(iesdr)

    if(iesdr.eq.-7.or.iesd.lt.0.or.iesd.eq.2) then
        iesiesd = .true.
    else
        iesiesd = .false.
    endif

    nelorb5 = nelorb * (indt4 + 1)
    nelorbj5 = nelorbj * (indt4j + 1)

    ! define in any case indtm
    do j = 1, nel
        indtm(j) = istart - 1
    enddo

    if(itest.eq.2) then  ! the true label of DMC/LRDMC (1) or VMC (2)
        yesfn = .false.
    else
        yesfn = .true.
    endif
#ifdef _DEBUG 
! Otherwise I cannot check averything
  movedions=1
  yesupwf=.true.
  yesupwfj=.true.
#else
    if(iflagnorm.eq.3.or.test_AAD) then
        movedions = 1
        yesupwf = .true.
        yesupwfj = .true.
    else
        if(indt4.gt.0.and.yesfn) then
            yesupwf = .false.
        else
            yesupwf = .true.
        endif
        if(indt4j.gt.0.and.yesfn) then
            yesupwfj = .false.
        else
            yesupwfj = .true.
        endif
    endif
#endif

    if(membigcpu.and..not.membig) then
     if(indt4.eq.0.and.yesupwf) then
     allocate(winv_big(ipc * nelorb, 0:indt+4, nel))
     winv_big=0.d0
     endif
     if(nelorbjh.gt.0.and.indt4j.eq.0.and.yesupwfj) then
     allocate(winvj_big(nelorbjmax, 0:indt+4,nel))
     winvj_big=0.d0
     endif
    endif

! membigcpu=.false. -> old algorithm anyway


    ! testing AAD
    if(test_aad) then
        save_jastrowall = .false.
        yesfn = .false.
        ind_winvbarfn = 1
        ind_winvfn = 1
    else
        save_jastrowall = .true.
    endif
    !      checking update
    if(yesfn) then
    ind_vpotsav=walker
    else
    ind_vpotsav=1
    endif

    psidetln_old = psidetln


    ! task #1  input kel rion --> output dist el-ion distance to be used later.
    ! write(6,*)  ' input dist = ',sum(dist(:,:))
#ifdef _TIME
    timep=cclock()
#endif
    call electron_ion_dist_TASK1(kelind, rion, dist)
#ifdef _TIME
    timings(1)=timings(1)+cclock()-timep
    timep=cclock()
#endif
    ! write(6,*)  ' output  dist = ',sum(dist(:,:))


    !   end task #1

    ! SS  begin task #2 --> Update all necessary to pseudo
    !  Output vpseudolocal,prefactor

    call task2(kelind, rion, dist, kel, vpseudolocal, prefactor, ivic, tmu)
#ifdef _TIME
    timings(2)=timings(2)+cclock()-timep
    timep=cclock()
#endif


    ! end task #2


    ! begin task # 3

    !     input rion(3,nion)  output iond(nion,nion)  ion-ion distances
    call task3(rion, iond)
    ! end task #3
#ifdef _TIME
    timings(3)=timings(3)+cclock()-timep
    timep=cclock()
#endif


    !     begin task #4
    !     evaluation gradients and laplacians of orbitals in the electron
    !     positions kel rion (with pseudo also in the mesh).
    !     in principle depend only on kelind for this case.
    !     output winv,winvj

    call task4(kel, rion, dd, vju, winv, winvj)

#ifdef _OFFLOAD
!if(yesfn) then
!$omp target update to (winv) if(yesupwf.or..not.yes_ontarget) 
!$omp target update to (winvj) if(yesupwfj.or..not.yes_ontarget)
!$omp target update from (winv) if(.not.yesupwf.and.yes_ontarget) 
!$omp target update from (winvj) if(.not.yesupwfj.and.yes_ontarget)
!else
! Too much latency in  the loop below.  It  is  better to transfer everything  once.
!!  In order  to save some comm GPU-CPU. NB update target works for consecutive elements in a vector so loop is necessary.
!do j=1,nel
!!$omp target update to (winv(1:ipc*nelorbh,0:0,j)) if(yesupwf.or..not.yes_ontarget) 
!!$omp target update to (winvj(1:nelorbjmax,0:0,j)) if(yesupwfj.or..not.yes_ontarget)
!!$omp target update from (winv(1:ipc*nelorbh,0:0,j)) if(.not.yesupwf.and.yes_ontarget) 
!!$omp target update from (winvj(1:nelorbjmax,0:0,j)) if(.not.yesupwfj.and.yes_ontarget)
!enddo
!endif
#endif

#ifdef _TIME
    count_allwf=count_allwf+dble(nel)*dble(nelorb+nelorbj)
    do i=1,nel
       do j=1,nelorb
       if(winv(j,0,i).eq.0.d0) count_zerowf=count_zerowf+1.d0 
       enddo
    enddo
    if(yes_crystalj) then
    do i=1,nel
       do j=1,nelorbj
       if(winvj(j,0,i).eq.0.d0) count_zerowf=count_zerowf+1.d0 
       enddo
    enddo
    endif
    timings(4)=timings(4)+cclock()-timep
    timep=cclock()
#endif

    !     end  task #4

    ! begin task #5  --> Update tables for speeding up the code
    !    Output winvbar,winvjbar,winvjbarsz

    call task5(winv, winvj, jasmat, jasmatsz, muj_c, jasmat_c, jasmatsz_c, detmat, mu_c&
            &, detmat_c, winvbar, winvjbar, winvjbarsz, winvbarfn(ind_winvbarfn))
#ifdef _TIME
    timings(5)=timings(5)+cclock()-timep
    timep=cclock()
#endif
    ! end task #5
    !  begin task #6
    !   evaluate classical potential energy
    !   input kelind(3,nel), (NB the mesh is not used)
    !   input rion(3,nion)
    !   input iond(nion,nion)
    !   output vpot (classical Coulomb energy).
    call task6(kelind, rion, iond, vpot)
#ifdef _TIME
    timings(6)=timings(6)+cclock()-timep
    timep=cclock()
#endif
    !  end task #6
    !  task #7  Update vpot
    !  input rion(3,nion), kelind(3,nel)
    !  output vpot, the rest not used.

    call task7(kelind, rion, vpseudolocal, vpot)

#ifdef _TIME
    timings(7)=timings(7)+cclock()-timep
    timep=cclock()
#endif

    !        one needs here nelorbj*nel*ip5+2*nel**2  scratch dimension psip
    !        +nel**2 if iessz on

    isdistp = nelorbjh * (indt + 5) + 27 * (indt + 1) * max(nshellj,nion)  ! exact max memory required for distp
    if(iessz) then
        ind1 = isdistp + 1
        ind2 = ind1 + 2 * nelorbjmax
        ind3 = ind2 + 2 * max(nel, indt + 5)
        ind4=ind3+nel+nion
        sizework1 = 2
    else
        ind1 = isdistp + 1
        ind2 = ind1 + ipj * nelorbjmax
        ind3 = ind2 + 2 * max(nel, indt + 5)
        ind4=ind3+nel+nion
        sizework1 = 1
    endif

    !  write(6,*) ' Memory required for psip here =',ipj,indt4j,ind3+nelorbjmax*(indt+5)-1+3*nel

    !  begin task #8
    !  given winv winvbar previously computed compute cofactor matrix and
    !  determinant wavefunction (without Jastrow)

    !  output psidetlnt,ainv

    call task8(winv, winvbar, psidetlnt, ainv, psip, psip(indpfaff))

    !   end task # 8
#ifdef _TIME
    timings(8)=timings(8)+cclock()-timep
    timep=cclock()
#endif



    !  task #9 see inside subroutine
    !   output tabpip,logpsi,tmu
#ifdef  _OFFLOAD
!$omp target update from (winvjbar) if(yesupwfj.or.yes_ontarget)
!$omp target update from (winvj) if(yes_ontarget.and..not.yesupwfj)
#endif
!   task9 is all inside the CPU 
    call  task9(kel, rion, winvj, winvjbar, winvjbarsz, ivic, tabpip&
            &, vj, vju, psidetlnt, logpsi, jastrowall_ee(1, 1, 0, walker), jastrowall_ei(1, 1, walker)&
            &, psip(ind1), psip(ind4), psip, psip(ind2),ipsip,psip(ind3))
    !  end task #9
#ifdef _TIME
    timings(9)=timings(9)+cclock()-timep
    timep=cclock()
#endif


    !   begin task #10
    !   Computing winvup winvdo to be used later for local eloc
    !   output ainvup ainvdo winvup,winvdo
    call  task10(ainv, winvbar, kel, rion, ainvup, ainvdo, winvup, winvdo, winvfn(ind_winvfn), psip)
    !  end task #10
#ifdef _TIME
    timings(10)=timings(10)+cclock()-timep
    timep=cclock()
#endif


    !     ind2=(20*indt+40)/nelorb+1
    !     ind1=20*indt+41-(ind2-1)*nelorb
    ! begin task # 11
    !     compute local energy given winvup winvdo tabpip
    !     output Local energy eloc
    call task11(winvup, winvdo, tabpip, tmu, vpot, eloc)
    ! end task #12
#ifdef _TIME
    timings(11)=timings(11)+cclock()-timep
    timep=cclock()
#endif

    !     restore the previous Cell and ewald vectors
    if(iespbc) then
        celldm(1:6) = celldmo(1:6)
        call InitCell(nion, nel, yes_complex)
        if(test_AAD) then
            call Ggen_samegrid
        else
            call Ggen
        endif
    endif
    ! check sign of the w.f.
    if(ipc.eq.2) then
        psisn = psidetlnt(2)

        if(abs(psisn - TWO_PI * nint(psisn / TWO_PI)).lt.Pi / 2.d0) then
            signnew = 1
        else
            signnew = -1
        endif

    else
        psisn = psidetsn
        signnew = nint(psisn)
    endif
    wf_sign = signnew

    psiln(walker) = logpsi(1)

    if(yesfn) psidetln = psidetlnt(1)
    if(wf_sign.ne.signold.and.signold.ne.0) then
        write(6, *) ' Warning updating sign wrong !!!,    &
                &              try to reduce nscra to be safe !!! '
    endif
#ifdef DEBUG
#ifdef _OFFLOAD
!$omp target  update from (winvbar)
!$omp target update from (winvbarfn(ind_winvbarfn:ind_winvbarfn+nel2barfn-1),&
!$omp& winvfn(ind_winvfn:ind_winvfn+nel2wtfn-1)) if(yesfn)
#endif
    if(rank.eq.0) then
       write(6,*)
       write(6,*) ' Main determinant matrices after compute_fast:',walker
       write(6,*) ' winvbar after    =',sum(abs(winvbar(1:ipf*ipc*nelorbh,1:nel_mat)))
       write(6,*) ' ainv after       =',sum(abs(ainv(1:ipc*nelup_mat,1:nelup_mat)))
       write(6,*) ' ainvup after     =',sum(abs(ainvup(1:ipc*nelup,1:nelorbh)))
       write(6,*) ' ainvdo after     =',sum(abs(ainvdo(1:ipc*max(neldo,1),1:nelorbh)))
       write(6,*) ' winvup after     =',sum(abs(winvup(1:ipc*nelup,1:indt+4)))
       write(6,*) ' winvdo after     =',sum(abs(winvdo(1:ipc*neldo,1:indt+4)))
       write(6,*) ' winvjbar after  =',sum(abs(winvjbar(1:int(ipj*nelorbjh,8)*nel)))
       write(6,*) ' jastrowall_ei = ',sum(abs(jastrowall_ei(1,:,walker)))
       write(6,*) ' jastrowall_ee =',sum(abs(jastrowall_ee(:,:,0,walker)))
       write(6,*) ' tabpip  = ',sum(abs(tabpip(:,:)))
       write(6,*) ' psidetln = ',psiln(walker),psidetln
       if(yesfn) then
       write(6,*) ' winvfn after     =',sum(abs(winvfn(1:ipc*nmolipf*(indt+4)*nel)))
       write(6,*) ' winvbarfn after  =',sum(abs(winvbarfn(1:ipc*nmol*nel)))
       write(6,*) ' vpotreg = ',sum(abs(vpotreg(:,1:nel)))

        endif
       write(6,*)
     endif
#endif
#ifdef _DEBUG
#ifdef _OFFLOAD
!$omp target  update from (winvbar)
!$omp target update from (winvbarfn(ind_winvbarfn:ind_winvbarfn+nel2barfn-1)) if(yesfn)
#endif
     if(signold.ne.0.and.rank.eq.0) then
     cost=abs(ainv_old-sum(abs(ainv(:,:))))/ainv_old
     if(cost.gt.1d-6) write(6,*) ' ERROR updating ainv',cost,ainv_old,sum(abs(ainv(:,:)))
!     if(cost.gt.1d-6) then
!     write(6,*) ' ERROR updating ainv',cost,ainv_old,sum(abs(ainv(:,:)))
!      write(6,*) ' Output ainv '
!     cost=0.d0
!      do i=1,nelup_mat
!         do j=i,nelup_mat
!         write(6,*) i,j,ainv(i,j),ainv(j,i),ainv_sav(i,j),ainv_sav(j,i)
!         cost=cost+abs(ainv(i,j)-ainv_sav(i,j))+abs(ainv(j,i)-ainv_sav(j,i))
!         write(6,*) ' ERROR=',cost
!         enddo
!      enddo
!     else
!       write(6,*) ' Updating ainv OK ',cost
!     endif
!     deallocate(ainv_sav)
     cost=abs(winvbar_old-sum(abs(winvbar(:,:))))/winvbar_old
     if(cost.gt.1d-6.and.itestrfn.ne.-6) write(6,*) ' ERROR updating winvbar ',cost,winvbar_old,sum(abs(winvbar(:,:)))
     cost=abs(jastrowee_old-sum(abs(jastrowall_ee(:,:,0,walker))))/jastrowee_old
     if(cost.gt.1d-6) write(6,*) ' ERROR updating jastrow_ee',cost,jastrowee_old,sum(abs(jastrowall_ee(:,:,0,walker)))
     if(jastrowei_old.ne.0.d0) then
     cost=abs(jastrowei_old-sum(abs(jastrowall_ei(1,:,walker))))/jastrowei_old
     else
     cost=0.d0
     endif
     if(cost.gt.1d-6) write(6,*) ' ERROR updating jastrow_ei',cost,jastrowei_old,sum(abs(jastrowall_ei(1,:,walker)))
     cost=abs(winvjbar_old-sum(abs(winvjbar(1:int(ipj*nelorbjh,8)*nel))))/winvjbar_old
     if(cost.gt.1d-6) write(6,*) ' ERROR updating winvjbar ',cost,winvjbar_old,sum(abs(winvjbar(1:ipj*nelorbjh*nel)))
     if(iessz) then
     cost=abs(winvjbarsz_old-sum(abs(winvjbarsz(1:nelorbjh*nel))))/winvjbarsz_old
     if(cost.gt.1d-6) write(6,*) ' ERROR updating winvjbarsz ',cost,winvjbarsz_old,sum(abs(winvjbarsz(1:nelorbjh*nel)))
     endif
     if(yesfn) then
     cost=abs(sum(abs(winvfn(ind_winvfn:ind_winvfn+nel2wtfn-1)))-winvfn_old)/winvfn_old
     if(cost.gt.1d-6) write(6,*) ' ERROR updating winvfn ',cost,sum(abs(winvfn(ind_winvfn:ind_winvfn+nel2wtfn-1))),winvfn_old
     cost=abs(sum(abs(winvbarfn(ind_winvbarfn:ind_winvbarfn+nel2barfn-1)))-winvbarfn_old)/winvbarfn_old
     if(cost.gt.1d-6) write(6,*) ' ERROR updating winvbarfn ',cost,sum(abs(winvbarfn(ind_winvbarfn:ind_winvbarfn+nel2barfn-1))),winvbarfn_old
     cost=abs(winvup_old-sum(abs(winvup(:,:))))/winvup_old
     if(cost.gt.1d-6) then
     write(6,*) ' ERROR updating winvup ',cost,winvup_old,sum(abs(winvup(:,:)))
!    do i=1,nelup*ipc
!    write(6,*) i,(winvup(i,j)-winvup_sav(i,j),j=1,indt+4)
!    enddo
     endif
     cost=abs(winvdo_old-sum(abs(winvdo(:,:))))/winvdo_old
     if(cost.gt.1d-6) then
     write(6,*) ' ERROR updating winvdo ',cost,winvdo_old,sum(abs(winvdo(:,:)))
     endif

     cost=abs(tabpip_old(1)-sum(abs(tabpip(:,1:indt))))
     if(cost.gt.1d-6) write(6,*) ' ERROR updating tabpip table = ',tabpip_old(1),sum(abs(tabpip(:,1:indt)))
     cost=abs(tabpip_old(2)-sum(abs(tabpip(:,indt+1:indt+4))))
     if(cost.gt.1d-6) write(6,*) ' ERROR updating tabpip lap = ',tabpip_old(2),sum(abs(tabpip(:,indt+1:indt+4)))
!    do i=1,neldo*ipc
!    write(6,*) i,(winvdo(i,j)-winvdo_sav(i,j),j=1,indt+4)
!    enddo
     cost=abs(winvuplap_old-sum(abs(winvup(:,indt+1:indt+4))))/winvuplap_old
    if(cost.gt.1d-6) write(6,*) ' ERROR updating winvup Grad/Lap ',cost,winvuplap_old,sum(abs(winvup(:,indt+1:indt+4)))
     cost=abs(winvdolap_old-sum(abs(winvdo(:,indt+1:indt+4))))/winvdolap_old
    if(cost.gt.1d-6) write(6,*) ' ERROR updating winvdo Grad/Lap ',cost,winvdolap_old,sum(abs(winvdo(:,indt+1:indt+4)))
    endif
     cost=abs(logpsi_old-psiln(walker))
     if(cost.gt.1d-6) write(6,*) ' ERROR updating logpsi ',cost,logpsi_old,psiln(walker)
     if(epscuttype.eq.2)  then
      cost=abs(agp_old-sum(abs(agp(:,:,walker))))/agp_old
  if(cost.gt.1.d-6) write(6,*) ' ERROR updating agp ',cost,agp_old,sum(abs(agp(:,:,walker)))
     endif
     if(indt.gt.0) then
     cost=abs(sum(abs(tmu(:,:)))-tmu_old)/tmu_old
   if((yesfn.or.lrdmc_deru).and.cost.gt.1d-6) write(6,*) ' ERROR in updating tmu =',cost,tmu_old,sum(abs(tmu(:,:)))
     endif
    if(yesfn) then
     cost=abs(vpotreg_old-sum(abs(vpotreg(:,1:nel))))/vpotreg_old
  if(cost.gt.1.d-6) write(6,*) ' ERROR updating vpotreg ',cost,vpotreg_old,sum(abs(vpotreg(:,1:nel)))
     cost=abs(kel_old-sum(kel(:,:,:)))
     if(cost.gt.1d-6) write(6,*)  ' ERROR kel ',walker,cost,kel_old,sum(kel(:,:,:))
    endif
     if(rank.eq.0) write(6,*) ' After updating '
     endif
#endif

    !write(6,*) 'nmol,nelup',nmol,nelup
    !do i = 1,2*nmol*nelup
    !   write(6,*) 'ups',i,winvbarfn(i)
    !enddo
    !stop

    iflagnorm = 1  ! everything is computed here


#ifdef _OFFLOAD
  if(yes_ontarget) then
#ifdef _CUSOLVER
!$omp target update to (agp(:,:,walker:walker)) if(ipf.eq.2)
!$omp target update to (agpn) if(epscuttype.eq.2.and.ipf.eq.2)
#else
!$omp target update to (agp(:,:,walker:walker)) 
!$omp target update to (agpn) if(epscuttype.eq.2)
#endif
!$omp target update to (winv) if(yesupwf.and..not.membigcpu)  ! no longer necessary (done after task 4)
!$omp target update to (winvj) if(yesupwfj.and..not.membigcpu) ! no longer necessary  (done after task 4)
  else
#ifdef _CUSOLVER
!$omp target update from (agp(:,:,walker:walker)) if(ipf.eq.1)  
!$omp target update from (agpn) if(epscuttype.eq.2.and.ipf.eq.1)
!$omp target update from (ainv) if(ipf.eq.1) 
#endif
!$omp target update from (winvbar) 
!$omp target update from (winvbarfn(ind_winvbarfn:ind_winvbarfn+nel2barfn-1)) if(yesfn) 
!$omp target update from (winvfn(ind_winvfn:ind_winvfn+nel2wtfn-1)) if(yesfn.and.yesupwf)
  endif
#endif
     if(allocated(winv_big)) deallocate(winv_big)
     if(allocated(winvj_big)) deallocate(winvj_big)

    return

contains

    subroutine task9(kel, rion, winvj, winvjbar, winvjbarsz&
            &, ivic, tabpip, vj, vju, psidetlnt, logpsi, jastrowall_ee, jastrowall_ei, work1, work2, work3, vecu,ispin_n,vold_n)
        implicit none
        !      input kel rion winvj,winvjbar,winvjbarsz (if spin Jastrow on)
        !      output tabpip(nel,indt+4),logpsi
        !      logpsi is the logarithm |psi| without the Determinantal part
        !      all the rest is just used to evaluate tabpip
        real(8), intent(in) :: kel(3, nel), rion(3, nion)&
                &, winvjbar(max(int(ipj * nelorbjh,8) * (nelup + neldo), 1)), winvjbarsz(sizewsz), psidetlnt(*)&
                &, ivic(3, indtmax, nelup + neldo), vj(*), vju(*)
        real(8), intent(inout) :: winvj(nelorbjmax, 0:indt4j, nel)
        real(8), intent(inout) :: tabpip(nel, indt + 4), logpsi(ipc)
        real(8) work1(ipj * nelorbjmax, sizework1), work2(nelorbjmax, 0:*), vecu(max(nel, indt + 5), 2)
        real(8) jastrowall_ee(nelup + neldo, nelup + neldo, 0:indt4j), jastrowall_ei(nion, nelup + neldo)
        real(8) rco(3), r0o, r0, wj, w0, costzj, csum
        real(8) work3(isdistp)
        integer nelorbjw, indtaj, indtmp,l
        integer*8 maxii,ii,indj
        integer ispin_n(nel)
        real(8) vold_n(nel+nion)
       real(8), external :: t_lrdmc, jastrow_ei, jastrow_ee
        !      input kel rion, vj
        !      output tabpip(nel,indt+4),logpsi,tmu
        !      logpsi is the logarithm |psi| without the Determinantal part
        !      all the rest is just used to evaluate tabpip

        !      Jastrow begin

        !      call dscalzero(4*nel,0.d0,tabpip(1,indt+1),1)
        !      tabpip(1:nel,indt+1:indt+4)=0.d0
#ifdef _NVTX
        call nvtxStartRange("Jastrow 3 body")
#endif
        tabpip = 0.d0
        jastrowall_ei = 0.d0
        jastrowall_ee = 0.d0
        maxii=nel
        maxii=(maxii*(nel-1))/2
!$omp parallel do default(shared) private(ii,j1,j2,kk,ispin,rc,jastrow1,grad1,grad2) reduction(+:tabpip) 
!       do j1 = 1, nel
!           do j2 = 1, nel
         do ii=1,maxii
                call find_j1j2(nel,ii,j1,j2) ! to collapse the loop by hand
                if(iesiesd) then

                    if((j1.le.nelup.and.j2.le.nelup).or.(j1.gt.nelup.and.j2.gt.nelup)) then
                        ! parallel spins
                        ispin = 1
                    else
                        ispin = -1
                    endif

                else
                    ispin = -1
                endif

!               if(j2.ne.j1) then
                    rc(:) = kel(:, j1) - kel(:, j2)
                    ! 2body Jastrow PBC/open
                    if(LBox.gt.0.d0) then
                        call jastrowgrad_pbc(rc, vj, iesd, jastrow1, grad1, grad2, ispin, 1.d0)
                    else
                        call jastrowgrad(rc, vj, iesd, jastrow1, grad1, grad2, ispin)
                    endif
                    ! jastrowall --> speed up Jastrow updating in ratiovar/uptabtot subroutines
                    if(save_jastrowall) then
                        if(yesfn) then
                            jastrowall_ee(j2, j1, 0) = jastrow1
                            jastrowall_ee(j2, j1, indt + 1:indt + 3) = grad1(:)
                            jastrowall_ee(j2, j1, indt + 4) = grad2
                            jastrowall_ee(j1, j2, 0) = jastrow1
                            jastrowall_ee(j1, j2, indt + 1:indt + 3) = -grad1(:)
                            jastrowall_ee(j1, j2, indt + 4) = grad2
                        else
                            jastrowall_ee(j1, j2, 0) = jastrow1
                            jastrowall_ee(j2, j1, 0) = jastrow1
                        endif
                    endif
                    ! update Jastrow table (for energy computation)
                    do kk = 1, 3
                        tabpip(j1, indt + kk) = tabpip(j1, indt + kk) + grad1(kk)
                        tabpip(j2, indt + kk) = tabpip(j2, indt + kk) - grad1(kk)
                    enddo
                    tabpip(j1, indt + 4) = tabpip(j1, indt + 4) + grad2
                    tabpip(j2, indt + 4) = tabpip(j2, indt + 4) + grad2
!           enddo
           enddo
!$omp end parallel do
            if(n_body_on.ne.0) then
!$omp parallel do default(shared) private(j1,j2,kk,rc,jastrow1,grad1,grad2) reduction(+:tabpip) 
              do j1 = 1, nel
                do j2 = 1, nion
                    if(LBox.gt.0.d0) then
                        rc(:) = kel(:, j1) - rion(:, j2)
                        call jastrowgrad_pbc(rc, vj(pointvj(1, j2))&
                                &, pointvj(2, j2), jastrow1, grad1, grad2, -1, costz(j2))
                        if(save_jastrowall) then
                            if(yesfn) then
                                jastrowall_ei(j2, j1) = jastrow1
                            else
                                jastrowall_ei(1, j1) = jastrowall_ei(1, j1) + costz3(j2) * jastrow1
                            endif
                        endif
                    else
                        rc(:) = (kel(:, j1) - rion(:, j2)) * costz(j2)

                        call jastrowgrad(rc, vj(pointvj(1, j2)), pointvj(2, j2), jastrow1, grad1, grad2, -1)

                        if(save_jastrowall) then
                            if(yesfn) then
                                jastrowall_ei(j2, j1) = jastrow1
                            else
                                jastrowall_ei(1, j1) = jastrowall_ei(1, j1) + costz3(j2) * jastrow1
                            endif
                        endif

                    endif

                    do kk = 1, 3
                        tabpip(j1, indt + kk) = tabpip(j1, indt + kk)                        &
                                & - grad1(kk) * costz(j2) * costz3(j2)
                    enddo
                    tabpip(j1, indt + 4) = tabpip(j1, indt + 4)                             &
                            & - grad2 * costz(j2)*costz(j2)*costz3(j2)
                enddo
        enddo
!$omp end parallel do
        endif

        if((ncore.gt.0.or.lrdmc_deru).or.yesfn) then


            do i = 1, indt
              do j = 1, nel
              tabpip(j, i) = 0.d0
              enddo
            enddo

            ! only if you are not doing VMC or SR or you have a non local potential
!$omp parallel do default(shared) private(j1,i,j2,ispin_n,vold_n,rc,rco,kk,r0,r0o,costzj)
            !do for the first particle
            do j1 = 1, nel
                !do for the nearest neighbours
!  computation  rco and vold ispin once for all
                do j2=1,nel
                        if(iesiesd) then
                            if((j1.le.nelup.and.j2.le.nelup).or.(j1.gt.nelup.and.j2.gt.nelup)) then
                                ! parallel spins
                                ispin_n(j2) = 1
                            else
                                ispin_n(j2) = -1
                            endif
                        else
                            ispin_n(j2) = -1
                        endif

                        if(j2.ne.j1) then
                            rco(:) = kel(:, j1) - kel(:, j2)

                            if(LBox.gt.0.d0) then
!                               call CartesianToCrystal(rco, 1)
         rco(:)=car2cry(:,1)*rco(1)+car2cry(:,2)*rco(2)+car2cry(:,3)*rco(3)
                                do kk = 1, 3
                                rco(kk) = map(rco(kk), cellscale(kk))
                                enddo
                            endif

                            vold_n(j2) = jastrow_ee(rco, vj, iesd, ispin_n(j2)) ! recompute vold
                        endif

                enddo
                if(n_body_on.ne.0) then
                do j2 = 1, nion !do for the second particle
                            if(LBox.gt.0.d0) then
                                rco(:) = kel(:, j1) - rion(:, j2)
!                               call CartesianToCrystal(rco, 1)
         rco(:)=car2cry(:,1)*rco(1)+car2cry(:,2)*rco(2)+car2cry(:,3)*rco(3)
                                do kk = 1, 3
                                    rco(kk) = costz(j2) * map(rco(kk), cellscale(kk))
                                enddo
                                r0o = norm_metric(rco, metric)
                            else
                                rco(:) = costz(j2) * (kel(:, j1) - rion(:, j2))
                                r0o = dsqrt(sum(rco(:)*rco(:)))
                            endif
                vold_n(nel+j2) = jastrow_ei(r0o, vj(pointvj(1, j2)), pointvj(2, j2)) ! recompute vold
                enddo
                endif
                do j2 = 1, nel
                  if(j2.ne.j1) then
                   do i = 1, indtm(j1)
                    !do for the second particle

                            rc(:) = kel(:, j1) + ivic(:, i, j1) - kel(:, j2)

                            if(LBox.gt.0.d0) then
!                               call CartesianToCrystal(rc, 1)
          rc(:)=car2cry(:,1)*rc(1)+car2cry(:,2)*rc(2)+car2cry(:,3)*rc(3)
                                do kk = 1, 3
                                    rc(kk) = map(rc(kk), cellscale(kk))
                                enddo
                            endif
                            costzj = jastrow_ee(rc, vj, iesd, ispin_n(j2))
                            if(yesfn) then
                                jastrowall_ee(j2, j1, i) = costzj
                                tabpip(j1,i) = tabpip(j1,i) + costzj - vold_n(j2)
                            else
                                tabpip(j1,i) = tabpip(j1,i) + costzj - vold_n(j2)
                            endif
                    enddo
                   endif
                enddo
!                   tabpip(j1, i) = tabpip(j1, i) * dexp(costz0)

                    if(n_body_on.ne.0) then
!                       costz0 = 0.d0
                     do j2 = 1, nion !do for the second particle
                      do i = 1, indtm(j1)
                            if(LBox.gt.0.d0) then
                                rc(:) = kel(:, j1) + ivic(:, i, j1) - rion(:, j2)
!                               call CartesianToCrystal(rc, 1)
       rc(:)=car2cry(:,1)*rc(1)+car2cry(:,2)*rc(2)+car2cry(:,3)*rc(3)
                                do kk = 1, 3
                                    rc(kk) = costz(j2) * map(rc(kk), cellscale(kk))
                                enddo
                                r0 = norm_metric(rc, metric)
                            else
                                rc(:) = costz(j2) * (kel(:, j1) + ivic(:, i, j1) - rion(:, j2))
                                r0 = dsqrt(sum(rc(:)*rc(:)))
                            endif
                            tabpip(j1,i) = tabpip(j1,i) - costz3(j2) * (jastrow_ei(r0, vj(pointvj(1, j2)), pointvj(2, j2)) - vold_n(nel+j2))
                      enddo
!                       tabpip(j1, i) = tabpip(j1, i) * dexp(costz0)
                     enddo
                    endif
                do i=1,indtm(j1)
                tabpip(j1,i)=dexp(tabpip(j1,i)) 
                enddo
            enddo
!$omp end parallel do
        endif  ! ncore.gt.0 or yesfn
#ifdef _NVTX
        call nvtxEndRange
#endif


        !!!! Jastrow two body end !!!!

#ifdef _NVTX
        call nvtxStartRange("Jastrow 2 body")
#endif
        psisn = 1
        logpsi = 0.d0

        call  up2bodypsi(nel, nelup, logpsi, kel, vj, iesdr  &
                &, rion, nion, costz, costz3, LBox, niesd)


        !!!!! 3body JASTROW !!!!

        if(nelorbj.ne.0) then
            ! evaluate work1
            !      work1(:,1)=winvjbar(1:nelorbju)
            !      do j=2,nel
            !         work1(:,1)=work1(:,1)+winvjbar((j-1)*nelorbju+1:j*nelorbju)
            !      enddo

            do k = 1, nelorbju
                work1(k, 1) = 0.d0
            enddo
!$omp parallel do default(shared) private(k,j)
            do k = 1, nelorbju
                do j = 1, nel
                    work1(k, 1) = work1(k, 1) + winvjbar(int(j - 1,8) * nelorbju + k)
                enddo
            enddo
            if(iessz) then
!$omp parallel do default(shared) private(k,j)
                do k = 1, nelorbjh
                    work1(k, 2) = winvjbarsz(k)
                    do j = 2, nelup
                        work1(k, 2) = work1(k, 2) + winvjbarsz((j - 1) * nelorbjh + k)
                    enddo
                    do j = nelup + 1, nel
                        work1(k, 2) = work1(k, 2) - winvjbarsz((j - 1) * nelorbjh + k)
                    enddo
                enddo
            endif


            !     Calculation wavefunction
            !
            ! Charge Jastrow
            !
            call dgemv('T', nelorbjh, nelused, 1.d0, winvj, nelorbj5, work1, 1, 0.d0, vecu, 1)

            do j = 1, nel
                vecu(j, 2) = 0.d0
            enddo
            if (ipj==2) then
                call dgemv('T', nelorbjh, neldo, 1.d0, winvj(1, 0, 1 + nelup), nelorbj5, work1(1 + nelorbjh, 1), 1, 0.d0, vecu(1 + nelup, 1), 1)

!$omp parallel do default(shared) private(k,j)
                do j = 1, nelup
                    do k = 1, nelorbjh
                        vecu(j, 2) = vecu(j, 2) + winvj(k, 0, j) * winvjbar(int(j - 1,8) * nelorbju + k)
                    enddo
                enddo
!$omp end parallel do
!$omp parallel do default(shared) private(k,j)
                do j = nelup + 1, nel
                    !     vecu(j,2)=sum(winvj(:,0,j)*winvjbar((j-1)*nelorbju+nelorbjh+1:j*nelorbju))
                    do k = 1, nelorbju - nelorbjh
                        vecu(j, 2) = vecu(j, 2) + winvj(k, 0, j) * winvjbar(int(j - 1,8) * nelorbju + nelorbjh + k)
                    enddo
                enddo
!$omp end parallel do
            else
!$omp parallel do default(shared) private(j,k)
                do j = 1, nel
                    !        vecu(j,2)=sum(winvj(:,0,j)*winvjbar((j-1)*nelorbjh+1:j*nelorbjh))
                    do k = 1, nelorbjh
                        vecu(j, 2) = vecu(j, 2) + winvj(k, 0, j) * winvjbar(int(j - 1,8) * nelorbjh + k)
                    enddo
                enddo
!$omp end parallel do
            endif

            do j = 1, nel
                !         wj=sum(winvj(:,0,j)*winvjbar((j-1)*nelorbjh+1:j*nelorbjh))
                !         costzj=sum(work1(:,1)*winvj(:,0,j))-wj
                !         costzj=vecu(j,1)-vecu(j,2)
                !         logpsi(1)=logpsi(1)+0.5d0*costzj
                !         if(save_jastrowall.and..not.yesfn) jastrowall_ee(1,j,0)=jastrowall_ee(1,j,0)+costzj
                logpsi(1) = logpsi(1) + 0.5d0 * (vecu(j, 1) - vecu(j, 2))
            enddo

            !
            ! Charge+Spin Jastrow
            !
            if(iessz) then
                call dgemv('T', nelorbjh, nel, 1.d0, winvj, nelorbj5, work1(1, 2), 1, 0.d0, vecu, 1)

!$omp parallel do default(shared) private(j)
                do j = 1, nel
                    vecu(j, 2) = sum(winvj(:, 0, j) * winvjbarsz((j - 1) * nelorbjh + 1:j * nelorbjh))
                enddo
!$omp end parallel do
                do j = 1, nelup
                    costzj = vecu(j, 1) - vecu(j, 2)
                    !            if(save_jastrowall.and..not.yesfn) jastrowall_ee(1,j,0)=jastrowall_ee(1,j,0)+costzj
                    logpsi(1) = logpsi(1) + 0.5d0 * costzj
                enddo
                do j = nelup + 1, nel
                    costzj = -(vecu(j, 1) + vecu(j, 2))
                    !            if(save_jastrowall.and..not.yesfn) jastrowall_ee(1,j,0)=jastrowall_ee(1,j,0)+costzj
                    logpsi(1) = logpsi(1) + 0.5d0 * costzj
                enddo
            endif


            do j = 1, nel
                if(indt4j.eq.0) then
                    if(membigcpu) then
                    work2(:, 0:indtm(j)) = winvj_big(:, 0:indtm(j), j)
                    work2(:, indt + 1:indt + 4) = winvj_big(:, indt + 1:indt + 4, j)
                    else
                    call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j)  &
                            &, nel, r, rmu, vju, zeta, rion, work3, work2, nelorbjh, nion, kionj     &
                            &, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                            &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
                    endif
                elseif(yesupwfj.and..not.membigcpu) then
                    !         else
                    call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j)  &
                            &, nel, r, rmu, vju, zeta, rion, work3, work2, nelorbjh, nion, kionj     &
                            &, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                            &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
                    winvj(:, 1:indtm(j), j) = work2(:, 1:indtm(j))
                    winvj(:, indt + 1:indt + 4, j) = work2(:, indt + 1:indt + 4)
                    if(indt.gt.indtm(j)) winvj(:, indtm(j) + 1:indt, j) = 0.d0
                else
                    work2(:, 0:indtm(j)) = winvj(:, 0:indtm(j), j)
                    work2(:, indt + 1:indt + 4) = winvj(:, indt + 1:indt + 4, j)
                endif



                !   Put consecutive information
                indtaj = indtm(j) + 5
                indj = int(j - 1,8) * nelorbju + 1
                indtmp = indtm(j) + 1
                if(indt.ne.indtm(j)) then
                    do k = 1, 4
                        work2(:, indtm(j) + k) = work2(:, indt + k)
                    enddo
                endif

                !          call dgemv('T',nelorbjh,indtaj,1.d0,work2,nelorbjh,winvjbar(indj),1,0.d0,vecu,1)
                !          call dgemv('T',nelorbjh,indtaj,1.d0,work2,nelorbjh,work1,1,0.d0,vecu(1,2),1)
                if (ipj==2) then
                    if (j<=nelup) then
                        call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, winvjbar(indj), 1, 0.d0, vecu, 1)
                        call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, work1, 1, 0.d0, vecu(1, 2), 1)

                    else
                        call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, winvjbar(indj + nelorbjh), 1, 0.d0, vecu, 1)
                        call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, work1(1 + nelorbjh, 1), 1, 0.d0, vecu(1, 2), 1)

                    endif
                else
                    call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, winvjbar(indj), 1, 0.d0, vecu, 1)
                    call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, work1, 1, 0.d0, vecu(1, 2), 1)
                endif


                do k = 1, indtm(j)
                    if(tabpip(j, k).ne.0.d0) then
                        tabpip(j, k) = tabpip(j, k) * exp(vecu(k + 1, 2) - vecu(1, 2) - vecu(k + 1, 1) + vecu(1, 1))
                        !        tabpip(j,k)=tabpip(j,k)*exp(sum((work2(:,k)-winvj(:,0,j))*work1(:,1))-wj+w0)
                    endif
                enddo
                do k = 1, 4
                    tabpip(j, k + indt) = tabpip(j, k + indt) + vecu(k + indtmp, 2) - vecu(k + indtmp, 1)
                    !         tabpip(j,k)=tabpip(j,k)+sum(work2(:,k)*work1(:,1))-wj
                enddo
             
               

                if(iessz) then
                    call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, winvjbarsz(indj), 1, 0.d0, vecu, 1)
                    call dgemv('T', nelorbjh, indtaj, 1.d0, work2, nelorbjh, work1(1, 2), 1, 0.d0, vecu(1, 2), 1)
                    if(j.le.nelup) then
                        do k = 1, indtm(j)
                            if(tabpip(j, k).ne.0.d0) then
                                tabpip(j, k) = tabpip(j, k) * exp(vecu(k + 1, 2) - vecu(1, 2) - vecu(k + 1, 1) + vecu(1, 1))
                                !       tabpip(j,k)=tabpip(j,k)*exp(sum((work2(:,k)-winvj(:,0,j))*work1(:,2))-wj+w0)
                            endif
                        enddo
                        do k = 1, 4
                            !         do k=indt+1,indt+4
                            tabpip(j, k + indt) = tabpip(j, k + indt) + vecu(k + indtmp, 2) - vecu(k + indtmp, 1)
                            !         tabpip(j,k+indt)=tabpip(j,k+indt)+sum(work2(:,k)*work1(:,2))-wj
                        enddo
                    else
                        do k = 1, indtm(j)
                            if(tabpip(j, k).ne.0.d0) then
                                tabpip(j, k) = tabpip(j, k) * exp(-vecu(k + 1, 2) + vecu(1, 2) - vecu(k + 1, 1) + vecu(1, 1))
                                !       tabpip(j,k)=tabpip(j,k)*exp(-sum((work2(:,k)-winvj(:,0,j))*work1(:,2))-wj+w0)
                            endif
                        enddo
                        !         do k=indt+1,indt+4
                        do k = 1, 4
                            tabpip(j, k + indt) = tabpip(j, k + indt) - vecu(k + indtmp, 2) - vecu(k + indtmp, 1)
                            !         tabpip(j,k)=tabpip(j,k)-sum(work2(:,k)*work1(:,2))-wj
                        enddo
                    endif
                endif ! endif iessz
            enddo ! end do j
#ifdef _NVTX
            call nvtxEndRange
#endif
        endif ! endif nelorbj




        ! update log(det(\Psi))
        logpsi(1:ipc) = logpsi(1:ipc) + psidetlnt(1:ipc)


        !write(6,*) ' jastrowall before=',sum(abs(jastrowall_ee)),sum(abs(jastrowall_ei))
        ! update jastrowall_ee in case of FN-DMC
        if(.not.test_aad) then 
#ifdef _OFFLOAD
!$omp target update to (winvjbar) if(.not.yes_ontarget.and..not.yesupwfj)
!$omp target update to (winvj) if(.not.membigcpu)
#endif
        call save_jall(yesfn, jastrowall_ee, winvjbar, winvjbarsz, winvj, psip)
        endif
#ifdef DEBUG
    write(6,*)
    write(6,*) ' Main Jastrow matrices after compute_fast'
    write(6,*) ' tabpip after     =',sum(abs(tabpip(1:nelup+neldo,1:indt)))
    write(6,*) ' jastrowall after =',sum(abs(jastrowall_ee)),sum(abs(jastrowall_ei))
#endif
        return
    end subroutine task9

    subroutine electron_ion_dist_TASK1(kel, rion, dist)
        implicit none
        real*8, intent(in) :: kel(3, nel), rion(3, nion)
        real*8, intent(inout) :: dist(nion, nel)
        ! SS  task #1  input kel rion --> output dist el-ion distance to be used later.
        !      reupdating distances
        if(LBox.gt.0.d0) then
            do j = 1, nel
                do k = 1, nion
                    dist_kel(:) = kel(:, j) - rion(:, k)
                    call ApplyPBC(dist_kel, 1)
                    dist(k, j) = sum(dist_kel(:)*dist_kel(:))
                    if(dist(k, j).lt.1d-18) dist(k, j) = 1d-18
                    dist(k, j) = dsqrt(dist(k, j))
                enddo
            enddo
        else
            do j = 1, nel
                do k = 1, nion
                    dist(k, j) = (kel(1, j) - rion(1, k))*(kel(1, j) - rion(1, k))+ &
                            & (kel(2, j) - rion(2, k))*(kel(2, j)-rion(2, k))+(kel(3, j) - rion(3, k))*(kel(3, j) - rion(3, k))
                    if(dist(k, j).lt.1d-18) dist(k, j) = 1d-18
                    dist(k, j) = dsqrt(dist(k, j))
                enddo
            enddo
        endif


        !  end task 1
    end subroutine electron_ion_dist_TASK1

    subroutine task2(kelind, rion, dist, kel, vpseudolocal, prefactor, ivic, tmu)
        implicit none
        real*8, intent(in) :: kelind(3, nel), rion(3, nion), dist(nion, nel)
        real*8, intent(inout) :: kel(3, nel, 0:indt), vpseudolocal&
                &, prefactor(indt - istart + 1, nel), ivic(3, indtmax, nelup + neldo), tmu(nel, max(indt, 1))
        real(8), external :: t_lrdmc
        real*8 sc(6), rn(3)

        ! choose a random direction (fillmatrix) and compute the integration mesh
        ! randomness can be avoided if
        ! ioptpseudo&iesrandoma=false
        !  input kelind(3,nel), rion,  dist(nel,nion) --->
        !  output:
        !  kel(3,j,1:indtm(j)) j=1,nel  integration mesh pseudo
        !  ivic (same information as above) used later.
        !  vpseudolocal  (scalar)            local part potential due to pseudo
        !  prefactor(indt-istart+1,nel)  information about non local part to be
        !  computed later


        ! angle is always an input
        ! In  the LRDDMC or DMC the randomization  is done only in initialization
        !  by uptabtot
        ! for consistency with continuation and branching
        if((ncore.gt.0.and.pseudologic).and..not.test_aad) then
            do j = 1, nel
                call fillmatrix(angle(10, j))
            enddo
        endif
        if(yesfn) then
            !     orthogonalize orthogonal matrix to avoid loosing accuracy
            do j = 1, nel
                call grahamo(angle(1, j), sc, rn, 3, 3, info)
                if(ncore.gt.0) call grahamo(angle(10, j), sc, rn, 3, 3, info)
            enddo
        endif

        if((ncore.gt.0.or.lrdmc_deru).or.yesfn) then
            ! if you are not doing VMC or SR or you have a non local potential
            if(itestrfn.eq.-6.or.itestrfn.eq.-7.or.itestrfn.eq.-1.or.lrdmc_deru) then

                !         if(iesrandoml) then
                !  randomize direction mesh if angle is not  the identity and define
                !  in ANY event ivic
                do j = 1, nel
                    do i = 1, istart - 1
                        call dgemv('N', 3, 3, 1.d0, angle(1, j), 3, versoralat(1, i) &
                                &, 1, 0.d0, ivic(1, i, j), 1)
                    enddo
                enddo
                !         endif
                do j = 1, nel
                    do i = 1, istart - 1
                        do kk = 1, 3
                            kel(kk, j, i) = kelind(kk, j) + ivic(kk, i, j)
                        enddo
                    enddo
                enddo
            endif

            if(ncore.gt.0) then
                ! if you are not doing VMC or SR or you have a non local potential

                vpseudolocal = 0.d0
                do j = 1, nel
                    !            do k=1,nion
                    !            r(k,0)=dist(j,k)
                    !            enddo
                    !  rn below is not used with iopt=.false. in pseudo.
                    call pseudoset(j, kelind(1, j), ivic(1, istart, j), prefactor   &
                            &, pseudolocal(j), nel, dist(1, j), rion, wpseudo, ncore  &
                            &, indt - istart + 1, indtm(j), .false., angle(10, j) &
                            &, psip, Lbox, 1d-9, iflag, rn)
                    if(iflag.ne.0) return

                    vpseudolocal = vpseudolocal + pseudolocal(j)
                    do i = istart, indtm(j)
                        do kk = 1, 3
                            kel(kk, j, i) = kelind(kk, j) + ivic(kk, i, j)
                        enddo
                    enddo
                enddo
            endif
        endif

        if(lrdmc_deru.or.yesfn) then
            do j = 1, nel
                do i = 1, istart - 1
                    tmu(j, i) = t_lrdmc(kelind(1, j), nion, rion, i, alat, ivic(1, 1, j), plat, cellscale, iespbc)
                enddo
            enddo
        endif

        do j = 1, nel
            do i = istart, indtm(j)
                tmu(j, i) = -prefactor(i - istart + 1, j)
            enddo
            do i = indtm(j) + 1, indt
                tmu(j, i) = 0.d0
            enddo
        enddo

        return

    end subroutine task2

    subroutine task4(kel, rion, dd, vju, winv, winvj)
        implicit none
        real*8, intent(in) :: kel(3, nel, 0:indt), rion(3, nion)
        real*8, intent(inout) :: winv(ipc * nelorbh, 0:indt4, nel)&
                &, winvj(nelorbjmax, 0:indt4j, nel), dd(*), vju(*)

        !     input  kel(3,nel,1:indt)
        !     output winv(nelorb,nel,indt+4),winvj(nelorbj,nel,indtj+4)
        !     all the rest is not used later.
        if(.not.yesupwf.and..not.yesupwfj) return

        if(yesupwf) winv=0.d0
        if(yesupwfj.and.nelorbj.ne.0) winvj=0.d0
        

        if(membigcpu.and.indt4.ne.0.and.indt4j.ne.0) then
        do j = 1, nelup
            if(yesupwf) then
                call upnewwf(indt, 0, indtm(j), 0, nshellh, ioptorb, ioccup, kel(1, j, 0)   &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .true.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
            endif
            if(nelorbj.ne.0.and.yesupwfj) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)          &
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo

        do j = nelup + 1, nel
            if(yesupwf) then
            call upnewwf(indt, 0, indtm(j), 0, nshelldoh, ioptorb, ioccdo, kel(1, j, 0) &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .false.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
            endif
            if(nelorbj.ne.0.and.yesupwfj) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)&
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .false.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo
        elseif(membigcpu.and.indt4.eq.0.and.indt4j.eq.0) then

! Allocated winv_big && winvj_big
        do j = 1, nelup
!           if(yesupwf) then
                call upnewwf(indt, 0, indtm(j), 0, nshellh, ioptorb, ioccup, kel(1, j, 0)   &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv_big(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .true.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
!           endif
            if(nelorbj.ne.0) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)          &
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj_big(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo

        do j = nelup + 1, nel
!           if(yesupwf) then
            call upnewwf(indt, 0, indtm(j), 0, nshelldoh, ioptorb, ioccdo, kel(1, j, 0) &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv_big(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .false.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
!           endif
            if(nelorbj.ne.0) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)&
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj_big(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .false.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo
        winv(:,0,1:nel)=winv_big(:,0,1:nel)
        if(nelorbj.ne.0) winvj(:,0,1:nel)=winvj_big(:,0,1:nel)
   

        elseif(membigcpu.and.indt4.eq.0) then

! allocated winv_big only winvj is already big
        do j = 1, nelup
!           if(yesupwf) then
                call upnewwf(indt, 0, indtm(j), 0, nshellh, ioptorb, ioccup, kel(1, j, 0)   &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv_big(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .true.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
!           endif
            if(nelorbj.ne.0.and.yesupwfj) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)          &
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo

        do j = nelup + 1, nel
!           if(yesupwf) then
            call upnewwf(indt, 0, indtm(j), 0, nshelldoh, ioptorb, ioccdo, kel(1, j, 0) &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv_big(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .false.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
!           endif
            if(nelorbj.ne.0.and.yesupwfj) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)&
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .false.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo
        winv(:,0,1:nel)=winv_big(:,0,1:nel)

        elseif(membigcpu.and.indt4j.eq.0) then

!       allocated winvj_big
        do j = 1, nelup
            if(yesupwf) then
                call upnewwf(indt, 0, indtm(j), 0, nshellh, ioptorb, ioccup, kel(1, j, 0)   &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .true.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
            endif
            if(nelorbj.ne.0) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)          &
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj_big(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo

        do j = nelup + 1, nel
            if(yesupwf) then
            call upnewwf(indt, 0, indtm(j), 0, nshelldoh, ioptorb, ioccdo, kel(1, j, 0) &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .false.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
            endif
            if(nelorbj.ne.0) then
                call upnewwf(indt, 0, indtm(j), 0, nshelljh, ioptorbj, ioccj, kel(1, j, 0)&
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj_big(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .false.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo
        if(nelorbj.ne.0) winvj(:,0,1:nel)=winvj_big(:,0,1:nel)

        else
        do j = 1, nelup
            if(yesupwf) then
                call upnewwf(0, 0, 0, 1, nshellh, ioptorb, ioccup, kel(1, j, 0)   &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .true.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
            endif
            if(nelorbj.ne.0.and.yesupwfj) then
                call upnewwf(0, 0, 0, 1, nshelljh, ioptorbj, ioccj, kel(1, j, 0)          &
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo

        do j = nelup + 1, nel
            if(yesupwf) then
                call upnewwf(0, 0, 0, 1, nshelldoh, ioptorb, ioccdo, kel(1, j, 0) &
                        &, nel, r, rmu, dd, zeta, rion, psip, winv(1, 0, j), nelorb, nion, kion         &
                        &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                        &, indpar_tab, indorb_tab, indshell_tab, .false.)
                !       if(yesupwfj) iflagnorm=-iflagnorm
            endif
            if(nelorbj.ne.0.and.yesupwfj) then
                call upnewwf(0, 0, 0, 1, nshelljh, ioptorbj, ioccj, kel(1, j, 0)&
                        &, nel, r, rmu, vju, zeta, rion, psip, winvj(1, 0, j), nelorbj   &
                        &, nion, kionj, iflagnorm, cnorm(nshell + 1), LBox, rmucos, rmusin, 1d-9&
                        &, indparj_tab, indorbj_tab, indshellj_tab, .false.)
                !      elseif(yesupwf) then
                !         iflagnorm=-iflagnorm
                !         if(iflagnorm.ne.1) iflagnorm=2
            endif
        enddo
        endif
        return

    end subroutine task4

    subroutine  task5(winv, winvj, jasmat, jasmatsz, muj_c, jasmat_c, jasmatsz_c, detmat, mu_c&
            &, detmat_c, winvbar, winvjbar, winvjbarsz, winvbarfn)

        implicit none
        real*8, intent(in) :: winv(ipc * nelorbh, 0:indt4, nel)&
                &, winvj(nelorbjmax, 0:indt4j, nel), jasmat(*), jasmatsz(*), detmat(ipc * ipf * nelorbh, *)&
                &, mu_c(*), detmat_c(*), muj_c(*), jasmat_c(*), jasmatsz_c(*)
        real*8, intent(inout) :: winvbar(ipf * ipc * nelorbh, nel_mat)&
                &, winvjbar(max(int(nelorbjh,8)* nel * ipj, 1)), winvjbarsz(sizewsz), winvbarfn(*)
        integer :: ind_diff
        integer*8 indw
! in lrdmc winvbar is not updated so here should  be computed by scratch

        ! to speed up the calculation introduced winvbar and winjbar
        ! input winv, winvj
        ! output winvbar and winvjbar (in terms of detmat and jasmat given
        ! input not changing) winvjbarsz (if spin Jastrow on)

        ! Evaluation of winvbar
        firstmmu = (firstmol - 1) * nelorbh * ipf * ipc + 1
        firstmdet = (firstmol - 1) * nelorb_c * ipc + ipc * (firstmol - 1) + 1
        firstmolr = ipc * nelorb_c * nelorb_c + ipc * (firstmol - 1) + 1

!       for  the unpaired
        call dscalzero_(size(winvbar),0.d0,winvbar,1)

        !   write(6,*) ' Input detmat_c '
        !   do i=1,nmol
        !   write(6,*) i,detmat_c(firstmdet+nelorb_c*(i-1)+i-1)&
        !   &,sum(abs(detmat_c(firstmdet+nelorb_c*(i-1):firstmdet+nelorb_c*i-1)))
        !   enddo

        if(yesfast.eq.0) then

            if(ipf.eq.2) then
                if(ipc.eq.2) then
                    call zgemm_('N', 'N', 2 * nelorbh, neldo, nelorbh, zone, detmat(1, nelorbh + 1), 2 * nelorbh&
                            &, winv(1, 0, nelup1), nelorb5, zzero, winvbar(1, nelup1), 2 * nelorbh)
                    call zgemm_('N', 'N', 2 * nelorbh, nelup, nelorbh, zone, detmat, 2 * nelorbh&
                            &, winv, nelorb5, zzero, winvbar, 2 * nelorbh)
                else
                    call dgemm_('N', 'N', 2 * nelorbh, neldo, nelorbh, 1.d0, detmat(1, nelorbh + 1), 2 * nelorbh&
                            &, winv(1, 0, nelup1), nelorb5, 0.d0, winvbar(1, nelup1), 2 * nelorbh)
                    call dgemm_('N', 'N', 2 * nelorbh, nelup, nelorbh, 1.d0, detmat, 2 * nelorbh&
                            &, winv, nelorb5, 0.d0, winvbar, 2 * nelorbh)
                endif
            else

                if(ipc.eq.2) then
                    call zgemm_('N', 'N', nelorbh, neldo, nelorbh, zone, detmat, nelorbh&
                            &, winv(1, 0, nelup1), nelorb5, zzero, winvbar(1, nelup1), nelorbh)
                    call zgemm_('T', 'N', nelorbh, nelup, nelorbh, zone, detmat, nelorbh, winv &
                            &, nelorb5, zzero, winvbar, nelorbh)
                else
                    call dgemm_('N', 'N', nelorbh, neldo, nelorbh, 1.d0, detmat, nelorbh&
                            &, winv(1, 0, nelup1), nelorb5, 0.d0, winvbar(1, nelup1), nelorbh)
                    call dgemm_('T', 'N', nelorbh, nelup, nelorbh, 1.d0, detmat, nelorbh, winv &
                            &, nelorb5, 0.d0, winvbar, nelorbh)
                endif
            endif ! endif ipf

            if(yesfn) then
                if(contraction.ne.0) then
                    if(ipc.eq.2) then
                        !  Computing winvbarfn using that detmat= m_c detmat_c mu_c ^T
                        if (ipf.eq.2) then
                            call zgemm_('T', 'N', nmol, nelup, nelorbh, zone, mu_c(firstmmu), 2 * nelorbh&
                                    &, winv, nelorb5, zzero, psip, nmol)
                            call zgemm_('T', 'N', nmol, neldo, nelorbh, zone, mu_c(firstmmu + 2 * nelorbh), 2 * nelorbh&
                                    &, winv(1, 0, nelup1), nelorb5, zzero, psip(2 * nmol * nelup + zerop1), nmol)
                            call zgemm_('N', 'N', nmol, nel, nmol, zone, detmat_c(firstmdet), nelorb_c&
                                    &, psip, nmol, zzero, winvbarfn, nmol)
                        else
                            call zgemm_('T', 'N', nmol, neldo, nelorbh, zone, mu_c(firstmmu), nelorbh&
                                    &, winv(1, 0, nelup1), nelorb5, zzero, psip, nmol)
                            call zgemm_('N', 'N', nmol, neldo, nmol, zone, detmat_c(firstmdet)&
                                    &, nelorb_c, psip, nmol, zzero, winvbarfn(ipc * nmol * nelup + zerop1), nmol)
                            call zgemm_('T', 'N', nmol, nelup, nelorbh, zone, mu_c(firstmmu), nelorbh&
                                    &, winv, nelorb5, zzero, psip, nmol)
                            call zgemm_('T', 'N', nmol, nelup, nmol, zone, detmat_c(firstmdet)&
                                    &, nelorb_c, psip, nmol, zzero, winvbarfn, nmol)
                        endif
                    else
                        if (ipf.eq.2) then
                            call dgemm_('T', 'N', nmol, nelup, nelorbh, 1.d0, mu_c(firstmmu), 2 * nelorbh&
                                    &, winv, nelorb5, 0.d0, psip, nmol)
                            call dgemm_('T', 'N', nmol, neldo, nelorbh, 1.d0, mu_c(firstmmu + nelorbh), 2 * nelorbh&
                                    &, winv(1, 0, nelup1), nelorb5, 0.d0, psip(nmol * nelup + zerop1), nmol)
                            call dgemm_('N', 'N', nmol, nel, nmol, 1.d0, detmat_c(firstmdet), nelorb_c&
                                    &, psip, nmol, 0.d0, winvbarfn, nmol)
                        else
                            call dgemm_('T', 'N', nmol, neldo, nelorbh, 1.d0, mu_c(firstmmu), nelorbh&
                                    &, winv(1, 0, nelup1), nelorb5, 0.d0, psip, nmol)
                            call dgemm_('N', 'N', nmol, neldo, nmol, 1.d0, detmat_c(firstmdet)&
                                    &, nelorb_c, psip, nmol, 0.d0, winvbarfn(nmol * nelup + zerop1), nmol)
                            call dgemm_('T', 'N', nmol, nelup, nelorbh, 1.d0, mu_c(firstmmu), nelorbh&
                                    &, winv, nelorb5, 0.d0, psip, nmol)
                            call dgemm_('T', 'N', nmol, nelup, nmol, 1.d0, detmat_c(firstmdet)&
                                    &, nelorb_c, psip, nmol, 0.d0, winvbarfn, nmol)
                        end if !ipf.eq.2
                    end if ! endif ipc=2
                else
             call dcopy_vec_(ipc * nelorbh * ipf * nel_mat, winvbar, winvbarfn)

                endif
            endif

        else ! yesfast=0

            if(ipc.eq.2) then
                if(yesprojm.and..not.yesfn) then
                    if (ipf.eq.2) then
                        !CCC prima contraggo a destra mu_c con winv (up e down separatamente)
                        ! Poi moltiplico con il projm per avere il winvbar
                        call zgemm_('T', 'N', nmol, nelup, nelorbh, zone, mu_c(firstmmu), 2 * nelorbh&
                                &, winv, nelorb5, zzero, psip, nmol)
                        call zgemm_('T', 'N', nmol, neldo, nelorbh, zone, mu_c(firstmmu + 2 * nelorbh), 2 * nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, zzero, psip(2 * nmol * nelup + zerop1), nmol)
                        call zgemm_('N', 'N', 2 * nelorbh, nel, nmol, zone, projm(firstmmu), 2 * nelorbh&
                                &, psip, nmol, zzero, winvbar, 2 * nelorbh)
                    else
                        call zgemm_('T', 'N', nmol, neldo, nelorbh, zone, mu_c(firstmmu), nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, zzero, psip, nmol)
                        call zgemm_('N', 'N', nelorbh, neldo, nmol, zone, projm(firstmmu), nelorbh&
                                &, psip, nmol, zzero, winvbar(1, nelup1), nelorbh)
                    end if
                else
                    if (ipf.eq.2) then
                        !CCC lo psip fino a 2*nmol*nelup_mat è il winvbarfn, poi è il mu_c^T winvbar
                        call zgemm_('T', 'N', nmol, nelup, nelorbh, zone, mu_c(firstmmu), 2 * nelorbh&
                                &, winv, nelorb5, zzero, psip(2 * nmol * nelup_mat + 1), nmol)
                        call zgemm_('T', 'N', nmol, neldo, nelorbh, zone, mu_c(firstmmu + 2 * nelorbh), 2 * nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, zzero, psip(2 * nmol * (nelup + nelup_mat) + zerop1), nmol)
                        if (yesprojm) then
                            call zgemm_('N', 'N', 2 * nelorbh, nel, nmol, zone, projm(firstmmu), 2 * nelorbh&
                                    &, psip(2 * nmol * nelup_mat + 1), nmol, zzero, winvbar, 2 * nelorbh)
                            call zgemm_('N', 'N', nmol, nel, nmol, zone, detmat_c(firstmdet), nelorb_c&
                                    &, psip(2 * nmol * nelup_mat + 1), nmol, zzero, psip, nmol)
                        else
                            call zgemm_('N', 'N', nmol, nel, nmol, zone, detmat_c(firstmdet), nelorb_c&
                                    &, psip(2 * nmol * nelup_mat + 1), nmol, zzero, psip, nmol)
                            call zgemm_('N', 'N', 2 * nelorbh, nel, nmol, zone, mu_c(firstmmu), 2 * nelorbh&
                                    &, psip, nmol, zzero, winvbar, 2 * nelorbh)
                        end if
                    else
                        call zgemm_('T', 'N', nmol, neldo, nelorbh, zone, mu_c(firstmmu), nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, zzero, psip(ipc * nmol * neldo + 1), nmol)
                        if(yesprojm) then
                            call zgemm_('N', 'N', nelorbh, neldo, nmol, zone, projm(firstmmu), nelorbh&
                                    &, psip(ipc * nmol * neldo + 1), nmol, zzero, winvbar(1, nelup1), nelorbh)
                            call zgemm_('N', 'N', nmol, neldo, nmol, zone, detmat_c(firstmdet), nelorb_c&
                                    &, psip(2 * nmol * neldo + 1), nmol, zzero, psip, nmol)
                        else
                            call zgemm_('N', 'N', nmol, neldo, nmol, zone, detmat_c(firstmdet), nelorb_c&
                                    &, psip(2 * nmol * neldo + 1), nmol, zzero, psip, nmol)
                            call zgemm_('N', 'N', nelorbh, neldo, nmol, zone, mu_c(firstmmu), nelorbh&
                                    &, psip, nmol, zzero, winvbar(1, nelup1), nelorbh)
                        endif
                    end if !ipf.eq.2
                endif

            else

                if(yesprojm.and..not.yesfn) then
                    if (ipf.eq.2) then
                        call dgemm_('T', 'N', nmol, nelup, nelorbh, 1.d0, mu_c(firstmmu), 2 * nelorbh&
                                &, winv, nelorb5, 0.d0, psip, nmol)
                        call dgemm_('T', 'N', nmol, neldo, nelorbh, 1.d0, mu_c((firstmmu + nelorbh)), 2 * nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, 0.d0, psip(nmol * nelup + zerop1), nmol)
                        call dgemm_('N', 'N', 2 * nelorbh, nel, nmol, 1.d0, projm(firstmmu), 2 * nelorbh&
                                &, psip, nmol, 0.d0, winvbar, 2 * nelorbh)
                    else
                        call dgemm_('T', 'N', nmol, neldo, nelorbh, 1.d0, mu_c(firstmmu), nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, 0.d0, psip, nmol)

                        !            write(6,*) ' psip after projm ',sum(abs(psip(1:nmol*neldo)))
                        call dgemm_('N', 'N', nelorbh, neldo, nmol, 1.d0, projm(firstmmu), nelorbh&
                                &, psip, nmol, 0.d0, winvbar(1, nelup1), nelorbh)
                    end if
                else
                    if (ipf.eq.2) then
                        !CCC lo psip fino a 2*nmol*nelup_mat è il winvbarfn, poi è il mu_c^T winvbar
                        call dgemm_('T', 'N', nmol, nelup, nelorbh, 1.d0, mu_c(firstmmu), 2 * nelorbh&
                                &, winv, nelorb5, 0.d0, psip(nmol * nelup_mat + 1), nmol)
                        call dgemm_('T', 'N', nmol, neldo, nelorbh, 1.d0, mu_c(firstmmu + nelorbh), 2 * nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, 0.d0, psip(nmol * (nelup + nelup_mat) + zerop1), nmol)
                        if (yesprojm) then
                            call dgemm_('N', 'N', 2 * nelorbh, nel, nmol, 1.d0, projm(firstmmu), 2 * nelorbh&
                                    &, psip(nmol * nelup_mat + 1), nmol, 0.d0, winvbar, 2 * nelorbh)
                            call dgemm_('N', 'N', nmol, nel, nmol, 1.d0, detmat_c(firstmdet), nelorb_c&
                                    &, psip(nmol * nelup_mat + 1), nmol, 0.d0, psip, nmol)
                        else
                            call dgemm_('N', 'N', nmol, nel, nmol, 1.d0, detmat_c(firstmdet), nelorb_c&
                                    &, psip(nmol * nelup_mat + 1), nmol, 0.d0, psip, nmol)
                            call dgemm_('N', 'N', 2 * nelorbh, nel, nmol, 1.d0, mu_c(firstmmu), 2 * nelorbh&
                                    &, psip, nmol, 0.d0, winvbar, 2 * nelorbh)
                        end if
                    else
                        call dgemm_('T', 'N', nmol, neldo, nelorbh, 1.d0, mu_c(firstmmu), nelorbh&
                                &, winv(1, 0, nelup1), nelorb5, 0.d0, psip(nmol * neldo + 1), nmol)
                        if(yesprojm) then
                            call dgemm_('N', 'N', nelorbh, neldo, nmol, 1.d0, projm(firstmmu), nelorbh&
                                    &, psip(nmol * neldo + 1), nmol, 0.d0, winvbar(1, nelup1), nelorbh)
                            call dgemm_('N', 'N', nmol, neldo, nmol, 1.d0, detmat_c(firstmdet), nelorb_c&
                                    &, psip(nmol * neldo + 1), nmol, 0.d0, psip, nmol)
                        else
                            call dgemm_('N', 'N', nmol, neldo, nmol, 1.d0, detmat_c(firstmdet), nelorb_c&
                                    &, psip(nmol * neldo + 1), nmol, 0.d0, psip, nmol)
                            call dgemm_('N', 'N', nelorbh, neldo, nmol, 1.d0, mu_c(firstmmu), nelorbh&
                                    &, psip, nmol, 0.d0, winvbar(1, nelup1), nelorbh)
                        endif
                    end if !ipf.eq.2
                endif
            endif

            if (yesfn) then
                if(ipf.eq.2) then
                    call dcopy_vec_(ipc * nmol * nel, psip,  winvbarfn)
                else
                 call dcopy_vec_(ipc * nmol * neldo, psip,  winvbarfn(ipc * nmol * nelup + zerop1))
                end if
            end if

            if (ipf.eq.1) then
                if(ipc.eq.2) then
                    if(yesprojm.or.yesfn) then
                        call zgemm_('T', 'N', nmol, nelup, nelorbh, zone, projm(firstmmu), nelorbh&
                                &, winv, nelorb5, zzero, psip, nmol)
                    else
                        call zgemm_('T', 'N', nmol, nelup, nelorbh, zone, mu_c(firstmmu), nelorbh&
                                &, winv, nelorb5, zzero, psip(2 * nmol * nelup + 1), nmol)

                        call zgemm_('T', 'N', nmol, nelup, nmol, zone, detmat_c(firstmdet), nelorb_c&
                                &, psip(2 * nmol * nelup + 1), nmol, zzero, psip, nmol)
                    endif
                    call zgemm_('N', 'N', nelorbh, nelup, nmol, zone, mu_c(firstmmu), nelorbh&
                            &, psip, nmol, zzero, winvbar, nelorbh)

                else

                    if(yesprojm.or.yesfn) then
                        call dgemm_('T', 'N', nmol, nelup, nelorbh, 1.d0, projm(firstmmu), nelorbh&
                                &, winv, nelorb5, 0.d0, psip, nmol)
                        !write(6,*) ' psip  projm =',sum(abs(psip(1:nmol*nelup)))
                    else
                        call dgemm_('T', 'N', nmol, nelup, nelorbh, 1.d0, mu_c(firstmmu), nelorbh&
                                &, winv, nelorb5, 0.d0, psip(nmol * nelup + 1), nmol)
                        call dgemm_('T', 'N', nmol, nelup, nmol, 1.d0, detmat_c(firstmdet), nelorb_c&
                                &, psip(nmol * nelup + 1), nmol, 0.d0, psip, nmol)
                    endif
                    !         write(6,*) ' input mu_c, psip =',firstmmu,nelorbh,nmol&
                    !  &,sum(abs(mu_c(firstmmu:firstmmu+nelorbh*nelup-1))),sum(abs(psip(1:nelup*nmol)))
                    call dgemm_('N', 'N', nelorbh, nelup, nmol, 1.d0, mu_c(firstmmu), nelorbh&
                            &, psip, nmol, 0.d0, winvbar, nelorbh)
                    !         write(6,*) ' output winvbar =',sum(abs(winvbar(:,1:nelup)))

                endif
                if(yesfn) call dcopy_vec_(ipc * nmol * nelup, psip, winvbarfn)
                !write(6,*) ' winbarfn after 2 =',sum(abs(winvbarfn(1:nmol*nelup))),sum(abs(winvbarfn(1:nmol*nel)))

                !write(6,*) ' winbar after 2 =',sum(abs(winvbar(:,1:nelup))),sum(abs(winvbar(:,1:nel)))
            endif
        end if


        !       Evaluation of winvjbar
        if(nelorbj.ne.0.and.yesupwfj) then
!           winvjbar = 0.d0
            if(yes_sparse) then
            winvjbar=0.d0
            else
            call dscalzero_(size(winvjbar),0.d0,winvjbar,1)
            endif
            if(contractionj.eq.0) then

                !up-up + down-up in the case ipj=2
                if(yes_sparse) then
                  if(ipj.eq.2) then
!$omp parallel do default(shared) private(i,j,ix,iy,indw) 
                  do j=1,nelup
                    do i=1,nnozeroj
                    iy=nozeroj(i+nnozeroj)
                    ix=nozeroj(i)
                     if(iy.le.nelorbjh) then 
                     indw=int(j-1,8)*nelorbjh2+ix
                     winvjbar(indw)=winvjbar(indw)+jasmat(i)*winvj(iy,0,j)
                     endif
                     if(ix.ne.iy.and.ix.le.nelorbjh) then
                     indw=int(j-1,8)*nelorbjh2+iy
                     winvjbar(indw)=winvjbar(indw)+jasmat(i)*winvj(ix,0,j)
                     endif
                    enddo
                   enddo
!$omp parallel do default(shared) private(i,j,ix,iy,indw) 
                   do j=nelup+1,nel
                    do i=1,nnozeroj
                    iy=nozeroj(i+nnozeroj)
                    ix=nozeroj(i)
                     if(iy.gt.nelorbjh) then 
                     indw=int(j-1,8)*nelorbjh2+ix
          winvjbar(indw)=winvjbar(indw)+jasmat(i)*winvj(iy-nelorbjh,0,j)
                     endif
                     if(ix.ne.iy.and.ix.gt.nelorbjh) then
                     indw=int(j-1,8)*nelorbjh2+iy
          winvjbar(indw)=winvjbar(indw)+jasmat(i)*winvj(ix-nelorbjh,0,j)
                     endif
                    enddo
                   enddo
                 else
!$omp parallel do default(shared) private(i,j,ix,iy,indw) 
                  do j=1,nel
                    do i=1,nnozeroj
                    iy=nozeroj(i)
                    ix=nozeroj(i+nnozeroj)
                    indw=int(j-1,8)*nelorbjh+ix
                     winvjbar(indw)=winvjbar(indw)+jasmat(i)*winvj(iy,0,j)
                     if(ix.ne.iy) then
                     indw=int(j-1,8)*nelorbjh+iy
                     winvjbar(indw)=winvjbar(indw)+jasmat(i)*winvj(ix,0,j)
                     endif
                    enddo
                  enddo
                 endif
#ifdef _OFFLOAD
!$omp target update to (winvjbar)
#endif
                else
                call dgemm_('N', 'N', ipj * nelorbjh, nelused, nelorbjh, 1.d0, jasmat, ipj * nelorbjh   &
                        &, winvj(1, 0, 1), nelorbj5, 0.d0, winvjbar, ipj * nelorbjh)

                !up-down + down-down

                 if (ipj==2) then
                    call dgemm_('N', 'N', ipj * nelorbjh, neldo, nelorbjh, 1.d0&
                            &, jasmat(nelorbju * nelorbjh + 1), ipj * nelorbjh, winvj(1, 0, nelup1), nelorbj5&
                            &, 0.d0, winvjbar(zerop1 + nelorbju * nelup), ipj * nelorbjh)
                 endif
                endif

                if(iessz) then
#ifdef  _OFFLOAD
!$omp target data map(from:winvjbarsz(1:nelorbjh*nel))
#endif
                    call dgemm_('N', 'N', nelorbjh, nel, nelorbjh, 1.d0, jasmatsz, nelorbjh &
                            &, winvj, nelorbj5, 0.d0, winvjbarsz, nelorbjh)
#ifdef  _OFFLOAD
!$omp end target data
#endif
                endif
            else

                !up-up

                call dgemm_('T', 'N', nelorbj_c, nelused, nelorbjh, 1.d0, muj_c, nelorbjh&
                        &, winvj, nelorbj5, 0.d0, psip(nelorbjcp), nelorbj_c)
                call dgemm_('N', 'N', nelorbj_c, nelused, nelorbj_c, 1.d0, jasmat_c, ipj * nelorbj_c&
                        &, psip(nelorbjcp), nelorbj_c, 0.d0, psip, nelorbj_c)
                call dgemm_('N', 'N', nelorbjh, nelused, nelorbj_c, 1.d0, muj_c, nelorbjh&
                        &, psip, nelorbj_c, 0.d0, winvjbar, ipj * nelorbjh)

                if (ipj==2) then

                    !down-up

                    !         call dgemm_('T','N',nelorbj_c,nelup,nelorbjh,1.d0,muj_c,nelorbjh&
                    !              &,winvj,nelorbj5,0.d0,psip(nelorbjcp),nelorbj_c)
                    call dgemm_('N', 'N', nelorbj_c, nelup, nelorbj_c, 1.d0, jasmat_c(1 + nelorbj_c), 2 * nelorbj_c& !jasmat (1+nelorbj_c,1)
                            &, psip(nelorbjcp), nelorbj_c, 0.d0, psip, nelorbj_c)
                    call dgemm_('N', 'N', nelorbjh, nelup, nelorbj_c, 1.d0, muj_c, nelorbjh&
                            &, psip, nelorbj_c, 0.d0, winvjbar(1 + nelorbjh), 2 * nelorbjh)     !winvjbar(1+nelorbjh,1)

                    !up-down
                    call dgemm_('T', 'N', nelorbj_c, neldo, nelorbjh, 1.d0, muj_c, nelorbjh&
                            &, winvj(1, 0, nelup1), nelorbj5, 0.d0, psip(nelorbjcp), nelorbj_c)
                    call dgemm_('N', 'N', nelorbj_c, neldo, nelorbj_c, 1.d0, jasmat_c(1 + 2 * nelorbj_c * nelorbj_c), 2 * nelorbj_c& !jasmat(1,1+nelorbj_c)
                            &, psip(nelorbjcp), nelorbj_c, 0.d0, psip, nelorbj_c)
                    call dgemm_('N', 'N', nelorbjh, neldo, nelorbj_c, 1.d0, muj_c, nelorbjh&
                            &, psip, nelorbj_c, 0.d0, winvjbar(zerop1 + nelorbju * nelup), 2 * nelorbjh)             !winvjbar(1,1+nelup)

                    !down-down

                    !         call dgemm_('T','N',nelorbj_c,neldo,nelorbjh,1.d0,muj_c,nelorbjh&
                    !              &,winvj(1,0,1+nelup),nelorbj5,0.d0,psip(nelorbjcp),nelorbj_c)
                    call dgemm_('N', 'N', nelorbj_c, neldo, nelorbj_c, 1.d0, jasmat_c(1 + 2 * nelorbj_c * nelorbj_c + nelorbj_c), 2 * nelorbj_c& !jasmat(1+nelorbj_c,1+nelorbj_c)
                            &, psip(nelorbjcp), nelorbj_c, 0.d0, psip, nelorbj_c)
                    call dgemm_('N', 'N', nelorbjh, neldo, nelorbj_c, 1.d0, muj_c, nelorbjh&
                            &, psip, nelorbj_c, 0.d0, winvjbar(zerop1 + nelup * nelorbju + nelorbjh), 2 * nelorbjh)  !winvjbar(1+nelup,1+nelorbjh)

                endif

                if(iessz) then
#ifdef  _OFFLOAD
!$omp target data map(from:winvjbarsz(1:nelorbjh*nel))
#endif
                    call dgemm_('N', 'N', nelorbj_c, nel, nelorbj_c, 1.d0, jasmatsz_c, nelorbj_c&
                            &, psip(nelorbjc), nelorbj_c, 0.d0, psip, nelorbj_c)
                    call dgemm_('N', 'N', nelorbjh, nel, nelorbj_c, 1.d0, muj_c, nelorbjh&
                            &, psip, nelorbj_c, 0.d0, winvjbarsz, nelorbjh)
#ifdef  _OFFLOAD
!$omp end target data
#endif
                endif
            endif
        endif

        if(ndiff.ne.0) then
            if(yesfast.eq.0) then
                !          do i=1,ndiff
                !             winvbar(1:ipc*nelorbh,nel+i)=detmat(1:ipc*nelorbh,nelorbh+i)
                !          enddo
                do i = 1, ndiff
!                   winvbar(1:ipc * ipf * nelorbh, nel + i) = detmat(1:ipc * ipf * nelorbh, ipf * nelorbh + i)
      call dcopy_vec_(ipc*ipf*nelorbh,detmat(1,ipf*nelorbh+i),winvbar(1,nel+i))
                    if(yesfn) then
                        if(contraction.ne.0) then
  ind_diff = (firstmol - 1) * (nelorb_c + i) * ipc + ipc * (firstmol - 1) + 1
  call dcopy_vec_(ipc*nmol, detmat_c(ind_diff), winvbarfn(ipc*nmol * (nel + i - 1) + 1))
                        else
             call dcopy_vec_(ipc*ipf*nelorbh, detmat(1, ipf * nelorbh + i),  winvbarfn(ipf*ipc*nelorbh * (nel + i - 1) + 1))
                        endif
                    endif
                enddo
            else
                if(yesprojm) then
                    do i = 1, ndiff
!                        winvbar(1:ipc * ipf * nelorbh, nel + i) = &
!                                &projm(ipc * ipf * nelorbh * (nelorb_c + i - 1) + 1:ipc * ipf * nelorbh * (nelorb_c + i))
             call dcopy_vec_(ipc * ipf * nelorbh,projm(ipc * ipf * nelorbh *&
         & (nelorb_c + i - 1) + 1),winvbar(1,nel+i))

                        !        call dcopy(nelorbh,projm(nelorbh*(nelorb_c+i-1)+1),1&
                        !    &,winvbar(1,nel+i),1)
                    enddo
                else
                    if(ipc.eq.2) then
                        call zgemm_('N', 'N', ipf * nelorbh, ndiff, nelorb_c, zone, mu_c, ipf * nelorbh&
                                &, detmat_c(2 * nelorb_c * nelorb_c + 1), nelorb_c, zzero, winvbar(1, nel + 1), ipf * nelorbh)
                    else
                        call dgemm_('N', 'N', ipf * nelorbh, ndiff, nelorb_c, 1.d0, mu_c, ipf * nelorbh&
                                &, detmat_c(nelorb_c * nelorb_c + 1), nelorb_c, 0.d0, winvbar(1, nel + 1), ipf * nelorbh)
                    endif
                endif
                ! fill winvbarfn in any case for FN-DMC
                if(yesfn) then
                    do i = 1, ndiff
                        if(contraction.ne.0) then
!         winvbarfn(ipc * nmol * (nel + i - 1) + 1:ipc * nmol * (nel + i)) = &
!                                    detmat_c(firstmolr + ipc * nelorb_c * (i - 1):firstmolr + ipc * nelorb_c * (i - 1) - 1 + ipc * nmol)
      call dcopy_vec_(ipc*nmol,detmat_c(firstmolr + ipc * nelorb_c * (i - 1))&
     &,winvbarfn(ipc * nmol * (nel + i - 1) + 1))
                        else
!                            winvbarfn(ipc * ipf * nelorbh * (nel + i - 1) + 1:ipc * ipf * nelorbh * (nel + i)) = &
!                                    detmat(1:ipc * ipf * nelorbh, ipf * nelorbh + i)
      call dcopy_vec_(ipc * ipf * nelorbh,detmat(1,ipf*nelorbh+i)&
     &,winvbarfn(ipc * ipf * nelorbh * (nel + i - 1) + 1))
                        endif
                    enddo
                endif
            endif
        endif

        return

    end subroutine task5

    subroutine task8(winv, winvbar, psidetlnt, ainv, psip, scratch_getri)
        implicit none
        real*8, intent(in) :: winv(ipc * nelorb, 0:indt4, nel), winvbar(ipc * ipf * nelorbh, nel_mat)
        real*8, intent(inout) :: psidetlnt(*), ainv(ipc * nelup_mat, nelup_mat)
        real*8, external :: dnrmsq,dnrmsq_
        real*8 :: psip(ipc * nelup_mat, nelup_mat)&
                &, scratch_getri(ipc * nbdgetri)
        complex*16 pphase
        integer i
        !  given winv winvbar previously computed compute cofactor matrix and
        !  determinant wavefunction (without Jastrow)

        !  input  winv, winvbar
        !  output psidetlnt,ainv(nelup,nelup) (cofactor matrix)



        !        AGP part wf.
        !        now compute the pairing function
        !
        !        call dscalzero(nelorb*nelup,0.d0,psip,1)
         call dscalzero_(size(psip),0.d0,psip,1)
        if(ipc.eq.2) then
            if (ipf.eq.2) then
                call zgemm_('T', 'N', nelup, nelup_mat, nelorbh, zone, winv               &
                        &, nelorb5, winvbar, 2 * nelorbh, zzero, psip, nelup_mat)
                call zgemm_('T', 'N', neldo, nelup_mat, nelorbh, zone, winv(1, 0, nelup1)&
                        &, nelorb5, winvbar(2 * nelorbh + 1, 1), 2 * nelorbh, zzero, psip(2 * nelup + zerop1, 1), nelup_mat)
#ifdef _OFFLOAD
!$omp target update from (psip)
#endif
                if(ndiff.ne.0) then
                    do k = 1, ndiff
                        do j = 1, nelup_mat - ndiff
                            psip(2 * nelup_mat - 2 * k + 1, j) = -psip(2 * j - 1, nelup_mat - k + 1)
                            psip(2 * nelup_mat - 2 * k + 2, j) = -psip(2 * j, nelup_mat - k + 1)
                        enddo
                    enddo

                    if(npar_eagp.gt.0) then
                        do k = 1, ndiff
                            do j = 1, ndiff
                                psip(2 * nelup_mat - 2 * ndiff + 2 * k - 1, nelup_mat - ndiff + j) = eagp_pfaff(2 * k - 1, j)
                                psip(2 * nelup_mat - 2 * ndiff + 2 * k, nelup_mat - ndiff + j) = eagp_pfaff(2 * k, j)
                            enddo
                        enddo
                    endif
                endif
                !       project antisymmetric
                do i = 1, nelup_mat
                    do j = 1, nelup_mat
                        if(j.gt.i) then
                            psip(2 * i - 1, j) = 0.5d0 * (psip(2 * i - 1, j) - psip(2 * j - 1, i))
                            psip(2 * i, j) = 0.5d0 * (psip(2 * i, j) - psip(2 * j, i))
                            psip(2 * j - 1, i) = -psip(2 * i - 1, j)
                            psip(2 * j, i) = -psip(2 * i, j)
                        endif
                    enddo
                    psip(2 * i - 1, i) = 0.d0
                    psip(2 * i, i) = 0.d0
                enddo

                if(epscuttype.eq.2) then
                    !      agp(1:2*nelup_mat,1:nelup_mat,walker)=psip(1:2*nelup_mat,1:nelup_mat)
                    do i = 1, 2 * nelup_mat
                        do j = 1, nelup_mat
                            agp(i, j, walker) = psip(i, j)
                        enddo
                    enddo
                endif
                call zsktrf('U', 'N', nelup_mat, psip, nelup_mat, ipsip, scratch_getri, nbdgetri, info)
                !         call zgetrf(nelup_mat,nelup_mat,psip,nelup_mat,ipsip,info)
            else
                call zgemm_('T', 'N', nelup, nelup, nelorbh, zone, winv&
                        &, nelorb5, winvbar(1, nelup + 1), nelorbh, zzero, psip, nelup)
#ifdef _OFFLOAD
#ifndef _CUSOLVER
!$omp target update from (psip)
#endif
#endif
                if(epscuttype.eq.2) then
#ifdef _CUSOLVER
!$omp target teams distribute parallel do collapse(2)
#endif
                    do j = 1, nelup
                        do i = 1, 2 * nelup
                            !      agp(1:2*nelup,1:nelup,walker)=psip(1:2*nelup,1:nelup)
                            agp(i, j, walker) = psip(i, j)
                        enddo
                    enddo
                endif
             call zgetrf_(nelup,nelup,psip,nelup,ipsip,info)

            end if
        else
            if(ipf.eq.2) then

                call dgemm_('T', 'N', nelup, nelup_mat, nelorbh, 1.d0, winv               &
                        &, nelorb5, winvbar, 2 * nelorbh, 0.d0, psip, nelup_mat)
                call dgemm_('T', 'N', neldo, nelup_mat, nelorbh, 1.d0, winv(1, 0, nelup1)&
                        &, nelorb5, winvbar(nelorbh + 1, 1), 2 * nelorbh, 0.d0, psip(nelup1, 1), nelup_mat)
#ifdef _OFFLOAD
!$omp target update from (psip)
#endif

                if(ndiff.ne.0) then
                    !       Equivalent to
                    !       call dgemv('T',nelorbh,nelup,-1.d0,winv,nelorb5,winvbar(1,nelup_mat),1,psip(nelup_mat,1),nelup_mat)
                    !       call dgemv('T',nelorbh,neldo,-1.d0,winv(1,0,nelup+1),nelorb5,winvbar(nelorbh+1,nelup_mat),1,psip(nelup_mat,nelup+1),nelup_mat)
                    do k = 1, ndiff
                        do i = 1, nelup_mat - ndiff
                            psip(nelup_mat - k + 1, i) = -psip(i, nelup_mat - k + 1)
                        enddo
                    enddo
                    if(npar_eagp.gt.0) then
                        do k = 1, ndiff
                            do j = 1, ndiff
                                psip(nelup_mat - ndiff + k, nelup_mat - ndiff + j) = eagp_pfaff(k, j)
                            enddo
                        enddo
                    endif
                endif
                !        write(6,*) ' Input pfaffian matrix ',ndiff,sum(abs(psip(1:nelup_mat,1:nelup_mat)))
                do i = 1, nelup_mat
                    do j = 1, nelup_mat
                        if(j.gt.i) then
                            psip(i, j) = 0.5d0 * (psip(i, j) - psip(j, i))
                            psip(j, i) = -psip(i, j)
                        endif
                    enddo
                    psip(i, i) = 0.d0
                enddo


                if(epscuttype.eq.2) then
                    do j = 1, nelup_mat
                        do i = 1, nelup_mat
                            !      agp(1:nelup,1:nelup,walker)=psip(1:nelup,1:nelup)
                            agp(i, j, walker) = psip(i, j)
                        enddo
                    enddo

                endif


                !       if(epscuttype.eq.2) agp(1:nelup_mat,1:nelup_mat,walker)=psip(1:nelup_mat,1:nelup_mat)

                call dsktrf('U', 'N', nelup_mat, psip, nelup_mat, ipsip, scratch_getri, nbdgetri, info)


            else
                call dgemm_('T', 'N', nelup, nelup, nelorbh, 1.d0, winv&
                        &, nelorb5, winvbar(1, nelup + 1), nelorbh, 0.d0, psip, nelup)
#ifdef _OFFLOAD
#ifndef _CUSOLVER
!$omp target update from (psip)
#endif
#endif
                if(epscuttype.eq.2) then
#ifdef _CUSOLVER
!$omp target teams distribute parallel do collapse(2)
#endif
                    do j = 1, nelup
                        do i = 1, nelup
                            !      agp(1:nelup,1:nelup,walker)=psip(1:nelup,1:nelup)
                            agp(i, j, walker) = psip(i, j)
                        enddo
                     enddo

                endif
                call dgetrf_(nelup, nelup, psip, nelup, ipsip, info)
            endif
        endif
        if(.not.test_AAD) singdet(walker) = .false.

        if(info.ne.0.and.ipf.eq.1) then

            do i = 1, nel
                write(6, *) i, ' Electron coordinates =', (kel(jj, i, 0), jj = 1, 3)
            enddo
            write(6, *) 'ERROR zero det info =/0 IN Z/DGETRF '
            !        call dscalzero(nelorb*nelup,0.d0,psip,1)
            call dscalzero_(size(psip),0.d0,psip,1)

            if(ipf.eq.2) then
                if(ipc.eq.2) then
                    call zgemm_('T', 'N', nelup, nelup_mat, nelorbh, zone, winv               &
                            &, nelorb5, winvbar, 2 * nelorbh, zzero, psip, nelup_mat)
                    call zgemm_('T', 'N', neldo, nelup_mat, nelorbh, zone, winv(1, 0, nelup1)&
                            &, nelorb5, winvbar(2 * nelorbh + 1, 1), 2 * nelorbh, zzero, psip(2 * nelup + zerop1, 1), nelup_mat)
#ifdef _OFFLOAD
!$omp target update from (psip)
#endif
                    if(ndiff.ne.0) then
                        do k = 1, ndiff
                            do j = 1, nelup_mat - ndiff
                                psip(2 * nelup_mat - 2 * k + 1, j) = -psip(2 * j - 1, nelup_mat - k + 1)
                                psip(2 * nelup_mat - 2 * k + 2, j) = -psip(2 * j, nelup_mat - k + 1)
                            enddo
                        enddo
                        if(npar_eagp.gt.0)  then
                            do k = 1, ndiff
                                do j = 1, ndiff
                                    psip(2 * nelup_mat - 2 * ndiff + 2 * k - 1, nelup_mat - ndiff + j) = eagp_pfaff(2 * k - 1, j)
                                    psip(2 * nelup_mat - 2 * ndiff + 2 * k, nelup_mat - ndiff + j) = eagp_pfaff(2 * k, j)
                                enddo
                            enddo
                        endif
                    endif
                else
                    call dgemm_('T', 'N', nelup, nelup_mat, nelorbh, 1.d0, winv               &
                            &, nelorb5, winvbar, 2 * nelorbh, 0.d0, psip, nelup_mat)
                    call dgemm_('T', 'N', neldo, nelup_mat, nelorbh, 1.d0, winv(1, 0, nelup1)&
                            &, nelorb5, winvbar(nelorbh + 1, 1), 2 * nelorbh, 0.d0, psip(nelup1, 1), nelup_mat)
#ifdef _OFFLOAD
!$omp target update from (psip)
#endif
                    if(ndiff.ne.0) then
                        do k = 1, ndiff
                            do  j = 1, nelup_mat - ndiff
                                psip(nelup_mat - k + 1, j) = -psip(j, nelup_mat - k + 1)
                            enddo
                        enddo
                        if(npar_eagp.gt.0) then
                            do k = 1, ndiff
                                do j = 1, ndiff
                                    psip(nelup_mat - ndiff + k, nelup_mat - ndiff + j) = eagp_pfaff(k, j)
                                enddo
                            enddo
                        endif
                    endif
                endif
            else
                if(ipc.eq.2) then
                    call zgemm_('T', 'N', nelup, nelup, nelorbh, zone, winv               &
                            &, nelorb5, winvbar(1, nelup1), nelorbh, zzero, psip, nelup)
                else
                    call dgemm_('T', 'N', nelup, nelup, nelorbh, 1.d0, winv               &
                            &, nelorb5, winvbar(1, nelup1), nelorbh, 0.d0, psip, nelup)
                endif
#ifdef _OFFLOAD
!$omp target update from (psip)
#endif
            endif

            write(6, *) ' input sum(|psip|) ',sum(abs(psip(1:nelup_mat*ipc,1:nelup_mat)))

!           do i = 1, nelup_mat
!               do j = 1, nelup_mat
!                   write(6, *) i, j, psip(i, j)
!               enddo
!           enddo

            psisn = 0
            singdet(walker) = .true.
            wconfn = 0.d0
            psidetsn = 0.d0
            iflag = 1
            return
        else

            if(ipc.eq.2) then
                call evaldet_complex(psip, nelup_mat, nelup_mat, ipsip, psidetlnt, psidetlnt(2))
                if(ipf.eq.2) then
                    call  zsktri('U', nelup_mat, psip, nelup_mat, ainv, nelup_mat, ipsip, scratch_getri, info)
                else
       call zgetri_(nelup_mat,psip,nelup_mat,ipsip,scratch_getri,nbdgetri,info)
                endif
            else
                call evaldet(psip, nelup_mat, nelup_mat, ipsip, psidetlnt, psidetsn)
                if(ipf.eq.2) then
                    call  dsktri('U', nelup_mat, psip, nelup_mat, ainv, nelup_mat, ipsip, scratch_getri, info)
                else
       call dgetri_(nelup_mat,psip,nelup_mat,ipsip,scratch_getri,nbdgetri,info)
                endif
            endif
            if(info.ne.0) then
                write(6, *) 'ERROR IN Z/DGETRI info =', info
                psisn = 0
                wconfn = 0.d0
                iflag = 1
                return
            endif

            if(ipf.ne.2) then
#ifdef _CUSOLVER
!$omp target teams distribute parallel do collapse(2)
#endif
                do j = 1, nelup_mat
                    do i = 1, ipc * nelup_mat
                        ainv(i, j) = psip(i, j)
                    enddo
                enddo
            endif
            ! antisymmetrization ainv
            !       write(6,*) ' Output matrix '
            !       cost=0.d0
            !       do i=1,nelup_mat
            !          write(6,*) i,i,ainv(i,i)
            !          cost=cost+abs(ainv(i,i))
            !          do j=i+1,nelup_mat
            !          write(6,*) i,j,ainv(i,j),ainv(j,i)
            !          cost=cost+abs(ainv(i,j)+ainv(j,i))
            !          enddo
            !       enddo
            !       write(6,*) ' ERROR asym old =',cost
            !       call mpi_finalize(ierr)
            !       stop
            if(ipf.eq.2) then
                if(ipc.eq.1) then
                    !       cost=0.d0
                    do i = 1, nelup_mat
                        do j = 1, nelup_mat
                            if(j.gt.i) then
                                ainv(i, j) = 0.5d0 * (ainv(i, j) - ainv(j, i))
                                ainv(j, i) = -ainv(i, j)
                            endif
                        enddo
                        ainv(i, i) = 0.d0
                    enddo

                    !      write(6,*) ' Loss symmetry start ',cost
                else
                    do i = 1, nelup_mat
                        do j = 1, nelup_mat
                            if(j.gt.i) then
                                ainv(2 * i - 1, j) = 0.5d0 * (ainv(2 * i - 1, j) - ainv(2 * j - 1, i))
                                ainv(2 * i, j) = 0.5d0 * (ainv(2 * i, j) - ainv(2 * j, i))
                                ainv(2 * j - 1, i) = -ainv(2 * i - 1, j)
                                ainv(2 * j, i) = -ainv(2 * i, j)
                            endif
                        enddo
                        ainv(2 * i - 1, i) = 0.d0
                        ainv(2 * i, i) = 0.d0
                    enddo
                endif
            endif

            


#ifdef _CUSOLVER
            if(ipf.eq.1) then
            psidetln = 1.d0 /sqrt(dnrmsq_(ipc * nelup_mat * nelup_mat, ainv, 1))
            else
            psidetln = 1.d0 /sqrt(dnrmsq(ipc * nelup_mat * nelup_mat, ainv, 1))
            endif
#else
            psidetln = 1.d0 /sqrt(dnrmsq(ipc * nelup_mat * nelup_mat, ainv, 1))
#endif
            if(.not.yesfn.and..not.test_aad.and.epscuttype.gt.0) then
                if(epscuttype.eq.2) then
#ifdef _CUSOLVER
!$omp target teams distribute parallel do collapse(2) if(ipf.eq.1)
#endif
                    do j = 1, nelup_mat
                        do i = 1, ipc * nelup_mat
                            agpn(i, j) = agp(i, j, walker)
                        enddo
                    enddo
                endif
#ifdef _CUSOLVER
             if(ipf.eq.1) then
             call psireg(ainv, agpn, psip, nelup_mat, neldo, psidetln,.true.)
             else
             call psireg(ainv, agpn, psip, nelup_mat, neldo, psidetln,.false.)
             endif
#else
             call psireg(ainv, agpn, psip, nelup_mat, neldo, psidetln,.false.)
#endif
                !         write(6,*) ' old new psidetln =',psidetln_old,psidetln
                !        For epscuttype<=0 psidetln is not updated
                cost = abs(psidetln_old - psidetln)
                ! When psidetln_old=0 (determinant=1) we require a relative  accuracy of 10^-6
                if(cost.gt.min(1d-2, sqrt(abs(psidetln_old))) * max(abs(psidetln_old), 1d-2).and.signold.ne.0) then
                    write(6, *) ' ERROR updating regularization =', walker, cost, psidetln, psidetln_old
                endif
            endif

        endif
#ifdef _OFFLOAD
#ifdef _CUSOLVER
!$omp target update to (ainv) if(ipf.eq.2)
!$omp target update to (agpn,agp(:,:,walker:walker)) if(yes_ontarget.and.ipf.eq.2) 
#else
!$omp target update to (ainv) 
!$omp target update to (agpn,agp(:,:,walker:walker)) if(yes_ontarget) 
#endif
#endif
    end subroutine task8

    subroutine task10(ainv, winvbar, kel, rion, ainvup, ainvdo, winvup, winvdo, winvfn, psip)
        implicit none
        real*8, intent(in) :: ainv(ipc * nelup_mat, nelup_mat), winvbar(ipc * ipf * nelorbh, nel_mat)&
                &, kel(3, nel, 0:indt), rion(3, nion)
        real*8, intent(inout) :: winvup(ipc * nelup, indt + 4), winvdo(max(ipc * neldo, 1), indt + 4)&
                &, ainvup(ipc * nelup, nelorbh), ainvdo(max(ipc * neldo, 1), nelorbh), psip(ipc * nelorbh,0: *)&
                &, winvfn(*)
        integer jr, ji, k, kr, ki, jelf,jj

        !   Computing winvup winvdo to be used later for local eloc and logpsi
        !   input ainv, winvbar, kel(3,nel,indt)
        !   output winvup,winvdo

        call dscalzero_(size(ainvup),0.d0,ainvup,1) 
        call dscalzero_(size(ainvdo),0.d0,ainvdo,1)
        if(ipc.eq.2) then
            if(ipf.eq.2) then
                call zgemm_('T', 'T', nelup, nelorbh, nelup_mat, zone, ainv, nelup_mat &
                        &, winvbar, 2 * nelorbh, zzero, ainvup, nelup)
            else
                call zgemm_('T', 'T', nelup, nelorbh, nelup, zone, ainv, nelup &
                        &, winvbar(1, nelup + 1), nelorbh, zzero, ainvup, nelup)
            endif
        else
            if(ipf.eq.2) then
                call dgemm_('T', 'T', nelup, nelorbh, nelup_mat, 1.d0, ainv, nelup_mat &
                        &, winvbar, 2 * nelorbh, 0.d0, ainvup, nelup)
            else
                call dgemm_('T', 'T', nelup, nelorbh, nelup, 1.d0, ainv, nelup &
                        &, winvbar(1, nelup + 1), nelorbh, 0.d0, ainvup, nelup)
            endif
        endif
#ifdef _OFFLOAD
!$omp target update from (ainvup,ainvdo)
#endif

        winvup = 0.d0
        winvdo = 0.d0

        if(.not.yesfn.and.ncore.le.0.and..not.lrdmc_deru) then

            ! VMC and SR case without pseudo
            do j = 1, nelup

                if((yesupwf.or.indt4.eq.0).and..not.membigcpu) then
                    call upnewwf(indt, 0, 0, 0, nshell, ioptorb, ioccup, kel(1, j, 0)   &
                            &, nel, r, rmu, dd, zeta, rion, psip(1, indt + ip5), psip, nelorbh, nion, kion         &
                            &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                            &, indpar_tab, indorb_tab, indshell_tab, .true.)

                    if(indt4.ne.0)  winv(:, indt + 1:indt + 4, j) = psip(:, indt + 1:indt + 4)
                else
                    if(indt4.eq.0) then
                    psip(:, indt + 1:indt + 4) = winv_big(:, indt + 1:indt + 4, j)
                    else
                    psip(:, indt + 1:indt + 4) = winv(:, indt + 1:indt + 4, j)
                    endif
                endif


                !  write(6,*) 'psip in task10',sum(abs(psip(1:ipc*nelorbh,1:indt+4)))

                if(ipc.eq.2) then
         call zgemv('T',nelorbh,4,zone,psip(1,indt+1),nelorbh,ainvup(2*j-1,1),nelup,zone,winvup(2*j-1,indt+1),nelup)
!                    jr = 2 * (j - 1) + 1
!                    ji = 2 * j
!!$omp parallel do default(shared) private(i,k,kr,ki) reduction(+:winvup)
!                    do k = 1, nelorbh
!                        do i = 1, 4
!                            kr = 2 * (k - 1) + 1
!                            ki = 2 * k
!                            winvup(jr, indt + i) = winvup(jr, indt + i) + ainvup(jr, k) * psip(kr, indt + i) - &
!                                    & ainvup(ji, k) * psip(ki, indt + i)
!                            winvup(ji, indt + i) = winvup(ji, indt + i) + ainvup(jr, k) * psip(ki, indt + i) + &
!                                    & ainvup(ji, k) * psip(kr, indt + i)
!                        enddo
!                    enddo
                else
         call dgemv('T',nelorbh,4,1.d0,psip(1,indt+1),nelorbh,ainvup(j,1),nelup,1.d0,winvup(j,indt+1),nelup)
!                 do k=1,nelorbh
!                   do i = 1, 4
!        winvup(j, indt + i) =winvup(j,indt+i)+ainvup(j,k)*psip(k, indt + i)
!                   enddo
!                 enddo
                endif
            enddo


        else ! pseudopotential or FN-DMC

            do j = 1, nelup

                ! compute basis set in the LRDMC lattice positions and
                ! basis gradients and laplacian
                if((yesupwf.or.indt4.eq.0).and..not.membigcpu) then
                    psip(1:ipc * nelorb, 1:indt) = 0.d0
                    call upnewwf(indt, 0, indtm(j), 0, nshellh, ioptorb, ioccup, kel(1, j, 0)   &
                            &, nel, r, rmu, dd, zeta, rion, psip(1, indt + ip5), psip, nelorbh, nion, kion      &
                            &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                            &, indpar_tab, indorb_tab, indshell_tab, .true.)
                    if(indt4.gt.0) winv(:, 1:indt + 4, j) = psip(:, 1:indt + 4)
                else
                    if(indt4.eq.0) then
                    psip(:, 1:indt + 4) = winv_big(:, 1:indt + 4, j)
                    else
                    psip(:, 1:indt + 4) = winv(:, 1:indt + 4, j)
                    endif
                endif



                ! fill winfn table for FN-DMC
                if(yesfn.and.yesupwf) then
#ifdef _OFFLOAD
!$omp target update to (psip(:,1:indt+4))
#endif
                    if(contraction.ne.0) then
                        jelf = nmolipf * (indt + ip4) * (j - 1) + 1
                        if(ipc.eq.2) then
                            call zgemm_('T', 'N', nmolipf, indt + ip4, nelorbh, zone, mu_c(firstmmu), ipf * nelorbh&
                                    &, psip(1,1), nelorb, zzero, winvfn(2 * jelf - 1), nmolipf)
                        else
                            call dgemm_('T', 'N', nmolipf, indt + ip4, nelorbh, 1.d0, mu_c(firstmmu), ipf * nelorbh&
                                    &, psip(1,1), nelorb, 0.d0, winvfn(jelf), nmolipf)
                        endif
                    else
                     jelf = ipc * nelorbh * (indt + ip4) * (j - 1) + 1
                     call dcopy_vec_(ipc*nelorbh*(indt+ip4),psip(1,1),winvfn(jelf))
!                        do i = 1, indt + ip4
!                            winvfn(ipc * nelorbh * (i - 1) + jelf:ipc * nelorbh
!* (i - 1) + ipc * nelorbh + jelf - 1) = &
!                                    psip(1:ipc * nelorbh, i)
!                        enddo



                    endif
                endif

                ! compute gradients/laplacian necessary for energy computation
                if(ipc.eq.2) then
         call zgemv('T',nelorbh,indt+4,zone,psip(1,1),nelorbh,ainvup(2*j-1,1),nelup,zone,winvup(2*j-1,1),nelup)
!                    jr = 2 * (j - 1) + 1
!                    ji = 2 * j
!!$omp parallel do default(shared) private(i,k,kr,ki) reduction(+:winvup)
!                  do k = 1, nelorbh
!                    do i = indt + 1, indt + 4
!                            kr = 2 * (k - 1) + 1
!                            ki = 2 * k
!                            winvup(jr, i) = winvup(jr, i) + ainvup(jr, k) * psip(kr, i) - ainvup(ji, k) * psip(ki, i)
!                            winvup(ji, i) = winvup(ji, i) + ainvup(ji, k) * psip(kr, i) + ainvup(jr, k) * psip(ki, i)
!                    enddo
!                    do i = 1, indtm(j)
!                            kr = 2 * (k - 1) + 1
!                            ki = 2 * k
!                            winvup(jr, i) = winvup(jr, i) + ainvup(jr, k) * psip(kr, i) - ainvup(ji, k) * psip(ki, i)
!                            winvup(ji, i) = winvup(ji, i) + ainvup(ji, k) * psip(kr, i) + ainvup(jr, k) * psip(ki, i)
!                    enddo
!                  enddo
                else
!!$omp parallel do default(shared) private(i,k) reduction(+:winvup)
         call dgemv('T',nelorbh,indt+4,1.d0,psip(1,1),nelorbh,ainvup(j,1),nelup,1.d0,winvup(j,1),nelup)
!                 do k=1,nelorbh
!                   do i = indt + 1, indt + 4
!                       winvup(j, i) = winvup(j,i)+ainvup(j, k) * psip(k, i)
!                   enddo
!                   do i = 1, indtm(j)
!                       winvup(j, i) = winvup(j,i)+ainvup(j, k) * psip(k, i)
!                   enddo
!                 enddo
                endif
            enddo
        endif


        ! spin down update
        if(.not.yesfn.and.ncore.le.0.and.neldo.gt.0.and..not.lrdmc_deru) then
            ! you are doing VMC or SR without pseudo
            ainvdo = 0.d0

            if(ipc.eq.2) then
                if(ipf.eq.2) then
                    call zgemm_('T', 'T', neldo, nelorbh, nelup_mat, zone, ainv(1, nelup1), nelup_mat &
                            &, winvbar(1 + 2 * nelorbh, 1), 2 * nelorbh, zzero, ainvdo, neldo)
                else
                    call zgemm_('N', 'T', neldo, nelorbh, nelup, zone, ainv, nelup         &
                            &, winvbar, nelorbh, zzero, ainvdo, neldo)
                endif
            else

                if(ipf.eq.2) then
                    call dgemm_('T', 'T', neldo, nelorbh, nelup_mat, 1.d0, ainv(1, nelup1), nelup_mat &
                            &, winvbar(1 + nelorbh, 1), 2 * nelorbh, 0.d0, ainvdo, neldo)
                else
                    call dgemm_('N', 'T', neldo, nelorbh, nelup, 1.d0, ainv, nelup         &
                            &, winvbar, nelorbh, 0.d0, ainvdo, neldo)
                endif
            endif
#ifdef  _OFFLOAD
!$omp target update from (ainvup,ainvdo)
#endif

            do j = 1, neldo
                if((yesupwf.or.indt4.eq.0).and..not.membigcpu) then
                    call upnewwf(indt, 0, 0, 0, nshellh, ioptorb, ioccup&
                            &, kel(1, j + nelup, 0), nel, r, rmu, dd, zeta, rion, psip(1, indt + ip5), psip, nelorbh&
                            &, nion, kion, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                            &, indpar_tab, indorb_tab, indshell_tab, .false.)
                    if(indt4.ne.0) winv(:, indt + 1:indt + 4, j + nelup) = psip(:, indt + 1:indt + 4)
                else
                    if(indt4.eq.0) then
                    psip(:, indt + 1:indt + 4) = winv_big(:, indt + 1:indt + 4, j + nelup)
                    else
                    psip(:, indt + 1:indt + 4) = winv(:, indt + 1:indt + 4, j + nelup)
                    endif
                endif

                if(ipc.eq.2) then
         call zgemv('T',nelorbh,4,zone,psip(1,indt+1),nelorbh,ainvdo(2*j-1,1),neldo,zone,winvdo(2*j-1,indt+1),neldo)
!                    jr = 2 * (j - 1) + 1
!                    ji = 2 * j
!!$omp parallel do default(shared) private(i,k,kr,ki) reduction(+:winvdo)
!                   do k = 1, nelorbh
!                        do i = 1, 4
!                            kr = 2 * (k - 1) + 1
!                            ki = 2 * k
!                            winvdo(jr, i + indt) = winvdo(jr, i + indt) + ainvdo(jr, k) * psip(kr, indt + i) - ainvdo(ji, k) * psip(ki, indt + i)
!                            winvdo(ji, i + indt) = winvdo(ji, i + indt) + ainvdo(ji, k) * psip(kr, indt + i) + ainvdo(jr, k) * psip(ki, indt + i)
!                        enddo
!                    enddo
                else
         call dgemv('T',nelorbh,4,1.d0,psip(1,indt+1),nelorbh,ainvdo(j,1),neldo,1.d0,winvdo(j,indt+1),neldo)
!!$omp parallel do default(shared) private(i,k) reduction(+:winvdo)
!                   do k = 1, nelorbh
!                    do i = 1, 4
!        winvdo(j, i + indt) = winvdo(j,i+indt)+ainvdo(j, k) * psip(k, indt + i)
!                    enddo
!                   enddo
                endif
            enddo

        elseif(neldo.gt.0) then
            !
            ! you are not doing standard VMC or SR
            call dscalzero_(size(ainvdo),0.d0,ainvdo,1)
!           ainvdo = 0.d0
            if(ipc.eq.2) then
                if(ipf.eq.2) then
                    call zgemm_('T', 'T', neldo, nelorbh, nelup_mat, zone, ainv(1, nelup + 1), nelup_mat &
                            &, winvbar(1 + 2 * nelorbh, 1), 2 * nelorbh, zzero, ainvdo, neldo)
                else
                    call zgemm_('N', 'T', neldo, nelorbh, nelup, zone, ainv, nelup         &
                            &, winvbar, nelorbh, zzero, ainvdo, neldo)
                endif
            else
                if(ipf.eq.2) then
                    call dgemm_('T', 'T', neldo, nelorbh, nelup_mat, 1.d0, ainv(1, nelup + 1), nelup_mat &
                            &, winvbar(1 + nelorbh, 1), 2 * nelorbh, 0.d0, ainvdo, neldo)
                else
                    call dgemm_('N', 'T', neldo, nelorbh, nelup, 1.d0, ainv, nelup         &
                            &, winvbar, nelorbh, 0.d0, ainvdo, neldo)
                endif
            endif
#ifdef _OFFLOAD
!$omp target update from (ainvup,ainvdo)
#endif
            do j = 1, neldo
                if((yesupwf.or.indt4.eq.0).and..not.membigcpu) then
                    psip(1:ipc * nelorb, 1:indt) = 0.d0
                    call upnewwf(indt, 0, indtm(j + nelup), 0, nshellh, ioptorb, ioccup&
                            &, kel(1, j + nelup, 0), nel, r, rmu, dd, zeta, rion, psip(1, indt + ip5), psip, nelorbh&
                            &, nion, kion, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                            &, indpar_tab, indorb_tab, indshell_tab, .false.)
                    if(indt4.ne.0) winv(:, 1:indt + 4, j + nelup) = psip(:, 1:indt + 4)

                else
                    if(indt4.eq.0) then
                    psip(:, 1:indt + 4) = winv_big(:, 1:indt + 4, j + nelup)
                    else
                    psip(:, 1:indt + 4) = winv(:, 1:indt + 4, j + nelup)
                    endif
                endif
                ! fill winvfn table for FN-DMC
                if(yesfn.and.yesupwf) then
#ifdef _OFFLOAD
!$omp target update to (psip(:,1:indt+4))
#endif
                    if(contraction.ne.0) then
                        jelf = nmolipf * (indt + ip4) * (j + nelup - 1) + 1

                        if(ipc.eq.2) then
                            call zgemm_('T', 'N', nmolipf, indt + ip4, nelorbh, zone, mu_c(firstmmu + nmolshift), ipf * nelorbh&
                                    &, psip(1,1), nelorb, zzero, winvfn(2 * jelf - 1), nmolipf)
                            !   call zgemm_('T','N',nmol/ipf,indt+ip4,nelorbh,zone,mu_c(firstmmu+(ipf-1)*(2*nmol+2)*nelorbh),ipf*nelorbh&
                            !               &,psip,nelorb,zzero,winvfn(2*jelf-1),nmol/ipf)
                        else
                            call dgemm_('T', 'N', nmolipf, indt + ip4, nelorbh, 1.d0, mu_c(firstmmu + nmolshift), ipf * nelorbh&
                                    &, psip(1,1), nelorb, 0.d0, winvfn(jelf), nmolipf)
                            !   call dgemm_('T','N',nmol/ipf,indt+ip4,nelorbh,1.d0,mu_c(firstmmu+(ipf-1)*(nmol+1)*nelorbh),ipf*nelorbh&
                            !               &,psip,nelorb,0.d0,winvfn(jelf),nmol/ipf)
                        endif

                    else
                        jelf = ipc * nelorbh * (indt + ip4) * (j + nelup - 1) + 1
                        !            if(ipf.eq.2) then
                     call dcopy_vec_(ipc*nelorbh*(indt+ip4),psip(1,1),winvfn(jelf))
!                        do i = 1, indt + ip4
!                            winvfn(ipc * nelorbh * (i - 1) + jelf:ipc * nelorbh
!* i + jelf - 1) = psip(1:ipc * nelorbh, i)
!                        enddo


                        !            else
                        !            do i=1,indt+ip4
                        !            winvfn(ipc*nelorbh*(i-1)+jelf:ipc*nelorbh*i+jelf-1)=psip(1:ipc*nelorbh,i)
                        !            enddo
                        !            endif
                    endif
                endif

                if(ipc.eq.2) then
         call zgemv('T',nelorbh,indt+4,zone,psip(1,1),nelorbh,ainvdo(2*j-1,1),neldo,zone,winvdo(2*j-1,1),neldo)
!                    jr = 2 * (j - 1) + 1
!                    ji = 2 * j
!!$omp parallel do default(shared) private(i,k,kr,ki) reduction(+:winvdo)
!                    do k = 1, nelorbh
!                        do i = indt + 1, indt + 4
!                            kr = 2 * (k - 1) + 1
!                            ki = 2 * k
!                            winvdo(jr, i) = winvdo(jr, i) + ainvdo(jr, k) * psip(kr, i) - ainvdo(ji, k) * psip(ki, i)
!                            winvdo(ji, i) = winvdo(ji, i) + ainvdo(ji, k) * psip(kr, i) + ainvdo(jr, k) * psip(ki, i)
!                        enddo
!                        do i = 1, indtm(j + nelup)
!                            kr = 2 * (k - 1) + 1
!                            ki = 2 * k
!                            winvdo(jr, i) = winvdo(jr, i) + ainvdo(jr, k) * psip(kr, i) - ainvdo(ji, k) * psip(ki, i)
!                            winvdo(ji, i) = winvdo(ji, i) + ainvdo(ji, k) * psip(kr, i) + ainvdo(jr, k) * psip(ki, i)
!                        enddo
!                    enddo

                else ! ipc=2
         call dgemv('T',nelorbh,indt+4,1.d0,psip(1,1),nelorbh,ainvdo(j,1),neldo,1.d0,winvdo(j,1),neldo)
!!$omp parallel do default(shared) private(i,k) reduction(+:winvdo)
!                    do k=1,nelorbh
!                        do i = indt + 1, indt + 4
!                        winvdo(j, i) = winvdo(j,i)+ainvdo(j,k)*psip(k, i)
!                        enddo
!                        do i = 1, indtm(j + nelup)
!                        winvdo(j, i) = winvdo(j,i)+ainvdo(j,k)*psip(k,i)
!                        enddo
!                    enddo
                endif
            enddo ! do j=1,neldo

        endif ! neldo gt 0



        return

    end subroutine task10

    subroutine task11(winvup, winvdo, tabpip, tmu, vpot, eloc)
        implicit none
        real*8, intent(in) :: winvup(ipc * nelup, indt + 4), winvdo(max(ipc * neldo, 1), indt + 4)&
                &, tabpip(nel, indt + 4), tmu(nelup + neldo, max(indt, 1)), vpot
        real*8, intent(inout) :: eloc(ipc)
        real*8 distsquare_node,weight_wagner
        integer ir
        real*8, external:: f_wagner

        !     compute local energy given winvup winvdo tabpip
        !     input winvup,winvdo tabpip tmu vpot
        !     output Local energy eloc

!       weight_wagner=1.d0
!       if(true_wagner) then
!        if(cutweight.gt.0.d0) then
!         distsquare_node=0.d0
!         do i=1,ipc*nelup
!         distsquare_node=distsquare_node+sum(winvup(i,indt+1:indt+3)**2)
!         enddo
!         do i=1,ipc*neldo
!         distsquare_node=distsquare_node+sum(winvdo(i,indt+1:indt+3)**2)
!         enddo
!         distsquare_node=1.d0/distsquare_node/cutweight**2 
!         weight_wagner=f_wagner(distsquare_node)
!        endif
!       endif

        call subener(indt, nelup, neldo, winvup, winvdo, tabpip, eloc)
        eloc(1) = eloc(1) + vpot
        !     contribution due to the non local pseudo

        if(ncore.gt.0.and..not.lrdmc_deru) then
            if(ipc.eq.1) then
                do jj = istart, indt ! to compute
                    do i = 1, nelup
                        eloc(1) = eloc(1) - tmu(i, jj) * tabpip(i, jj) * winvup(i, jj)
                    enddo

                    do i = 1, neldo
                        eloc(1) = eloc(1) - tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(i, jj)
                    enddo
                enddo
            else
                do jj = istart, indt ! to compute
                    do i = 1, nelup
                        ir = 2 * i - 1
                        eloc(1) = eloc(1) - tmu(i, jj) * tabpip(i, jj) * winvup(ir, jj)
                        eloc(2) = eloc(2) - tmu(i, jj) * tabpip(i, jj) * winvup(ir + 1, jj)
                    enddo
                    do i = 1, neldo
                        ir = 2 * i - 1
                        eloc(1) = eloc(1) - tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir, jj)
                        eloc(2) = eloc(2) - tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir + 1, jj)
                    enddo
                enddo
            endif
        endif
        if(lrdmc_deru) then

            if(ipc.eq.1) then
                do jj = 1, istart - 1 ! to compute
                    do i = 1, nelup
                        if(-tmu(i, jj) * winvup(i, jj).gt.0) then
                            eloc(1) = eloc(1) - gamma * tmu(i, jj) * tabpip(i, jj) * winvup(i, jj)
                        else
                            eloc(1) = eloc(1) + tmu(i, jj) * tabpip(i, jj) * winvup(i, jj)
                        endif
                        eloc(1) = eloc(1) - tmu(i, jj) !  the diagonal term of the discretized kinetic energy
                    enddo
                    do i = 1, neldo
                        if(-tmu(i + nelup, jj) * winvdo(i, jj).gt.0) then
                            eloc(1) = eloc(1) - gamma * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(i, jj)
                        else
                            eloc(1) = eloc(1) + tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(i, jj)
                        endif
                        eloc(1) = eloc(1) - tmu(i + nelup, jj) !  the diagonal term of the discretized kinetic energy
                    enddo
                enddo
                do jj = istart, indt ! to compute
                    do i = 1, nelup
                        if(-tmu(i, jj) * winvup(i, jj).gt.0) then
                            eloc(1) = eloc(1) - (1.d0 + gamma) * tmu(i, jj) * tabpip(i, jj) * winvup(i, jj)
                        elseif(npow.gt.0.d0) then
                            eloc(1) = eloc(1) - npow*(1.d0 + gamma) * tmu(i, jj) * tabpip(i, jj) * winvup(i, jj)
                        endif
                    enddo
                    do i = 1, neldo
                        if(-tmu(i + nelup, jj) * winvdo(i, jj).gt.0) then
                            eloc(1) = eloc(1) - (1.d0 + gamma) * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(i, jj)
                        elseif(npow.gt.0.d0) then
                            eloc(1) = eloc(1) - npow*(1.d0 + gamma) * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(i, jj)
                        endif
                    enddo
                enddo
            else
                do jj = 1, istart - 1 ! to compute
                    do i = 1, nelup
                        ir = 2 * i - 1
                        if(-tmu(i, jj) * winvup(ir, jj).gt.0) then
                            eloc(1) = eloc(1) - gamma * tmu(i, jj) * tabpip(i, jj) * winvup(ir, jj)
                            eloc(2) = eloc(2) - gamma * tmu(i, jj) * tabpip(i, jj) * winvup(ir + 1, jj)
                        else
                            eloc(1) = eloc(1) + tmu(i, jj) * tabpip(i, jj) * winvup(ir, jj)
                            eloc(2) = eloc(2) + tmu(i, jj) * tabpip(i, jj) * winvup(ir + 1, jj)
                        endif
                        eloc(1) = eloc(1) - tmu(i, jj) !  the diagonal term of the discretized kinetic energy
                    enddo
                    do i = 1, neldo
                        ir = 2 * i - 1
                        if(-tmu(i + nelup, jj) * winvdo(ir, jj).gt.0) then
                            eloc(1) = eloc(1) - gamma * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir, jj)
                            eloc(2) = eloc(2) - gamma * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir + 1, jj)
                        else
                            eloc(1) = eloc(1) + tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir, jj)
                            eloc(2) = eloc(2) + tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir + 1, jj)
                        endif
                        eloc(1) = eloc(1) - tmu(i + nelup, jj) !  the diagonal term of the discretized kinetic energy
                    enddo
                enddo
                do jj = istart, indt ! to compute
                    do i = 1, nelup
                        ir = 2 * i - 1
                        if(-tmu(i, jj) * winvup(ir, jj).gt.0) then
                            eloc(1) = eloc(1) - (1.d0 + gamma) * tmu(i, jj) * tabpip(i, jj) * winvup(ir, jj)
                            eloc(2) = eloc(2) - (1.d0 + gamma) * tmu(i, jj) * tabpip(i, jj) * winvup(ir + 1, jj)
                        elseif(npow.gt.0.d0) then
                            eloc(1) = eloc(1) - npow*(1.d0 + gamma) * tmu(i, jj) * tabpip(i, jj) * winvup(ir, jj)
                            eloc(2) = eloc(2) - npow*(1.d0 + gamma) * tmu(i, jj) * tabpip(i, jj) * winvup(ir + 1, jj)
                        endif
                    enddo
                    do i = 1, neldo
                        ir = 2 * i - 1
                        if(-tmu(i + nelup, jj) * winvdo(ir, jj).gt.0) then
                            eloc(1) = eloc(1) - (1.d0 + gamma) * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir, jj)
                            eloc(2) = eloc(2) - (1.d0 + gamma) * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir + 1, jj)
                        elseif(npow.gt.0.d0) then
                            eloc(1) = eloc(1) - npow*(1.d0 + gamma) * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir, jj)
                            eloc(2) = eloc(2) - npow*(1.d0 + gamma) * tmu(i + nelup, jj) * tabpip(i + nelup, jj) * winvdo(ir + 1, jj)
                        endif
                    enddo
                enddo
            endif
        endif

!   if(lrdmc_deru) eloc(1)=eloc(1)*weight_wagner

    end subroutine task11

    subroutine task3(rion, iond)
        !     input rion(3,nion)  output iond(nion,nion)  ion-ion distances
        !     ion_cart,psip  not used later.
        implicit none
        real*8, intent(in) :: rion(3, nion)
        real*8, intent(inout) :: iond(nion, nion)
        call eval_iond(iond, rion, nion, LBox, psip, iond_cart)
    end subroutine task3

    subroutine  task6(kelind, rion, iond, vpot)
        implicit none
        real*8, intent(in) :: kelind(3, nel), rion(3, nion), iond(nion, nion)
        real*8, intent(inout) :: vpot
        !   evaluate classical potential energy
        !   input kelind(3,nel), (NB the mesh is not used)
        !   input rion(3,nion)
        !   input iond(nion,nion)
        !   output vpot (classical Coulomb energy).
        ! compute the potential energy
        call upvpotdiag(kelind, nel, zeta, rion, iond, vpot, nion, oldkappa, LBox&
         &, vpotreg, cutreg, costz, costz3,vpotsav_ee(1,ind_vpotsav),yesfn)
        return
    end subroutine task6

    subroutine  task7(kelind, rion, vpseudolocal, vpot)
        implicit none
        real*8, intent(in) :: kelind(3, nel), rion(3, nion), vpseudolocal
        real*8, intent(inout) :: vpot

        !!********* Add Ewald Sums *******************
        ! if PBC on correct vpot for Ewald contribution.
        !  input rion(3,nion), kelind(3,nel)
        !  output vpot, the rest not used.
        if(LBox.gt.0.d0) then
            call EwaldSum(kelind, rion, ukwald, vpotreg, walker)
            !write(6,*) 'vpot',vpot,ukwald,sum(vpotreg(1:2,1:nel))
            vpot = vpot + ukwald
        endif
        !!********* Add Ewald Sums *******************
        ! adding to vpot vpseudolocal previously computed.
        if(ncore.gt.0) then
            vpot = vpot + 2.d0 * vpseudolocal
            do j = 1, nel
                vpotreg(1:2, j) = vpotreg(1:2, j) + 2.d0 * pseudolocal(j)
            enddo
        endif

        return
    end subroutine task7

end subroutine compute_eloc_logpsi
    function f_wagner(dist)
    real*8 f_wagner,dist
    if(dist.le.1) then
    f_wagner=7*dist*dist*dist-15*dist*dist+9*dist
    else
    f_wagner=1.d0
    endif
    end function f_wagner
   
    subroutine  fwagner_b(dist,distb,f_wagnerb)
    use allio, only: true_wagner
    real*8 f_wagnerb,dist,distb
    if(dist.le.1.and.true_wagner.gt.0) then
!   f_wagner=7*dist**3-15*dist**2+9*dist
    distb=distb+f_wagnerb*(21*dist*dist-30*dist+9)
    endif
    f_wagnerb=0.d0
    end subroutine fwagner_b
    subroutine find_j1j2(n,i,j1,j2)
         implicit none
         integer*8 i,j,jd,maxi
         integer n,j1,j2
         real*8 delta
         maxi=n
         maxi=(maxi*(n-1))/2
         j=1+4*(2*maxi-2*i)
         delta=sqrt(dble(j))
         jd=nint(delta)
         if(jd*jd.eq.j) then
         j=(1+jd)/2
         else
         j=(1.d0+delta)/2
         endif
         j1=n-j
         j2=j1+(maxi-(j*(j-1))/2)-i+1
     end  subroutine find_j1j2
