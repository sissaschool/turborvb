!TL off
SUBROUTINE UPNEWWF_B(indt, i0, indtm, typecomp, nshell, ioptorb, iocc, x&
        &, xb, nel, r, rb, rmu, rmub, dd, ddb, zeta, rion, rionb, distp, distpb&
        &, z, zb, nelskip, nion, kion, iflagnorm, cnorm, lbox, rmucos, rmucosb, &
        &  rmusin, rmusinb, mindist, cellscale, cellscaleb, s2r, s2rb&
        &, indpar_tab, indorb_tab, indshell_tab, adr_nion, ind_nion)
    use allio, only :yesupel,rank,ikshift,kgrid,kgrid_atom,rmunew,rknew,npip&
            !                 &,rmunewb,rknewb,iespbc,yesupel,nelorbh,indt
     &,rmunewb, rknewb, iespbc,yesupel,yes_crystalj,lepsbas&
     &,yes_scemama,novec_loop1,slaterorb_read,nshell_det
    use Cell, ONLY : cellpi, cellscalep, rphasep, rphase, cosphase, cosphaseb, sinphase, &
            &  sinphaseb, gamma_point, phase2pi, phase2pi_down, at, s2rp, car2cry
    use Constants, only : ipc, zzero, pi, ip4

    IMPLICIT NONE
    INTEGER :: nel, indorb, indpar, indshell, iocc(*), ioptorb(*), j, i, &
            &  indt, nelskip, nshell, indtmin, nion, k, kion(*), iflagnorm, indtm, i0&
            &, kk, ii, typecomp, kpip(3), ll, jj, iii, jjj, ishift&
            &, dimp, indp1, indp2, indp3,i0u,mshift
    REAL*8 :: dd(*), ddb(*), rmu(3, 0:indtm,nion), rion(3, nion)&
            &, distp((indtm + 1) * 23 * nshell + nelskip * (indt + 4 - i0 + 1))&
            &, cnorm(*), zeta(*), r(0:indtm,nion), rcos_sav(3&
            &), mindist, rcos(3), rsin(3), rcosb(3), rsinb(3), rphs(0:indtm)
    INTEGER :: indpar_tab(*), indorb_tab(*), indshell_tab(*), adr_nion(*), ind_nion(*)
    REAL*8 :: rmub(3, 0:indtm,nion), rionb(3, nion)&
            &, distpb((indtm + 1) * 23 * nshell + nelskip * (indt + 4 - i0 + 1))&
            &, zb(nelskip * min(ipc, indpar_tab(1) + 2), i0:*)&
            &, z(nelskip * min(ipc, indpar_tab(1) + 2), i0:*)&
            &, rb(0:indtm,nion), rmusav(3), cellscale(3), cellscaleb(3)&
            &, rmup(3), s2r(3, 3), s2rb(3, 3)
    REAL*8 :: x(3, nel, 0:indtm)
    REAL*8 :: xb(3, nel, 0:indtm)
    REAL*8 :: lbox, lbox_sav, rmucos(3, 0:indtm,nion), rmusin(3, 0:indtm,nion)
    REAL*8 :: rmucosb(3, 0:indtm,nion), rmusinb(3, 0:indtm,nion)
    REAL*8 :: rphaseb(3), vecscra(3)
    !  real*8, dimension(:,:), allocatable:: psip_r,psip_rb
    !  real*8 :: psip_rb(nelorbh,0:indt+ip4),psip_r(nelorbh,0:indt+ip4)
    DOUBLE PRECISION :: x1, x2, cellpib(3)
    DOUBLE PRECISION :: temp2b0
    logical iesjas,gammaorj
    real*8 phs(3)
    integer case_if,case_upz
!   logical, external:: slaterorb
    logical do_makefun
    logical yeszero_z(0:indtm)

   if(indtm.eq.0.and.i0.eq.0.and.novec_loop1) then


      CALL UPNEWWF0_B(indt, typecomp, nshell, ioptorb, iocc, x&
        &, xb, nel, r, rb, rmu, rmub, dd, ddb, zeta, rion, rionb, distp, distpb&
        &, z, zb, nelskip, nion, kion, iflagnorm, cnorm, lbox, rmucos, rmucosb, &
        &  rmusin, rmusinb, mindist, cellscale, cellscaleb, s2r, s2rb&
        &, indpar_tab, indorb_tab, indshell_tab, adr_nion, ind_nion)

   else
    if(typecomp.eq.1) then
    i0u=i0
    else
    i0u=0
    endif


    if(indpar_tab(1).eq.-1) then
        iesjas = .true.
        !  indpar_tab(1)=0
        mshift=nshell_det
    else
        iesjas = .false.
        mshift=0
    endif

    LBox_sav = LBox
    ishift = 0
    dimp = (indtm + 1) * 20
    if(LBox.eq.3) then
        if(yesupel) then
            phs(:) = phase2pi(:)
        else
            phs(:) = phase2pi_down(:)
        endif
    else
        phs = 0.d0
    endif
    if(abs(LBox).eq.3.and.iesjas) then
        ! always use the old periodic basis for the Jastrow
        if(.not.yes_crystalj) then
            if(iespbc) then
                !      iflagnorm=2
                LBox = 1.d0
            else
                LBox = -1.d0
            endif
        else
            ishift = ikshift
        endif
        phs(:) = 0.d0
        !  elseif(abs(LBox).eq.3) then
        !  allocate(psip_r(nelskip,i0:indt+4),psip_rb(nelskip,i0:indt+4))
    endif

            if(sum(abs(phs(:))).eq.0.d0) then
            gammaorj=.true.
            else
            gammaorj=.false.
            endif
            if(ipc.eq.2.and..not.gammaorj) then ! complex case
            case_if=1
            elseif(.not.gammaorj.and.ipc.eq.1) then
            case_if=2
            else if(ipc.eq.2.and..not.iesjas) then
            case_if=3
            else
            case_if=4
            endif

            if(ipc.eq.2.and..not.iesjas) then ! complex case
                if(typecomp.eq.1) then
                case_upz=1
                else
                case_upz=2
                endif
            else ! real case or  Jastrow
                if(typecomp.eq.1) then
                case_upz=3
                else
                case_upz=4
                endif
            endif

    !  if(ishift.ne.0) write(6,*) ' dimension work3 inside =',nelskip*(indt+5+(20*(indt+1))/nelskip-i0+1)

    indshell = 0
    indorb = 0
    indpar = 0
    indtmin = 0
    cellpib = 0.d0
    rphaseb = 0.d0
    ! used only here

    if(abs(LBox).ne.3) then
        DO i = indtmin, indtm
            DO k = 1, nion
                rmu(1, i,k) = x(1, 1, i) - rion(1, k)
                rmu(2, i,k) = x(2, 1, i) - rion(2, k)
                rmu(3, i,k) = x(3, 1, i) - rion(3, k)
            END DO
        END DO
        ! ****** Periodic Systems *************
        IF (lbox .eq. 1.d0) THEN
            IF (.NOT.gamma_point.and..not.iesjas) THEN
                DO j = indtmin, indtm
                    DO k = 1, nion
                        rcos_sav(1) = DCOS(rphase(1) * rmu(1, j,k))
                        rcos_sav(2) = DCOS(rphase(2) * rmu(2, j,k))
                        rcos_sav(3) = DCOS(rphase(3) * rmu(3, j,k))
                        cosphase(j,k) = rcos_sav(1) * rcos_sav(2) * rcos_sav(3)
                        sinphase(1, j,k) = DSIN(rphase(1) * rmu(1, j,k)) * rcos_sav(2)&
                                & * rcos_sav(3)
                        sinphase(2, j,k) = DSIN(rphase(2) * rmu(2, j,k)) * rcos_sav(1)&
                                & * rcos_sav(3)
                        sinphase(3, j,k) = DSIN(rphase(3) * rmu(3, j,k)) * rcos_sav(1)&
                                & * rcos_sav(2)
                    END DO
                END DO
            END IF
            DO j = indtmin, indtm
                DO k = 1, nion
                    rmu(1, j,k) = rmu(1, j,k) / cellpi(1)
                    rmu(2, j,k) = rmu(2, j,k) / cellpi(2)
                    rmu(3, j,k) = rmu(3, j,k) / cellpi(3)
                    rmucos(1, j,k) = DCOS(rmu(1, j,k))
                    rmucos(2, j,k) = DCOS(rmu(2, j,k))
                    rmucos(3, j,k) = DCOS(rmu(3, j,k))
                    rmusin(1, j,k) = DSIN(rmu(1, j,k))
                    rmusin(2, j,k) = DSIN(rmu(2, j,k))
                    rmusin(3, j,k) = DSIN(rmu(3, j,k))
                    rmu(1, j,k) = cellpi(1) * rmusin(1, j,k)
                    rmu(2, j,k) = cellpi(2) * rmusin(2, j,k)
                    rmu(3, j,k) = cellpi(3) * rmusin(3, j,k)
                    x1 = DSQRT(rmu(1, j,k)**2 + rmu(2, j,k)**2 + rmu(3, j,k)**2&
                            &)
                    IF (x1 .LT. mindist) THEN
                        r(j,k) = mindist
                    ELSE
                        r(j,k) = x1
                    END IF
                END DO
            END DO
        ELSEIF(LBox.eq.2) then
            cellpib = 0.d0
            rphaseb = 0.d0
            if(gamma_point.or.iesjas) then
                do i = indtmin, indtm
                    do k = 1, nion
                        kpip(1) = anint(rmu(1, i,k) / cellscale(1))
                        kpip(2) = anint(rmu(2, i,k) / cellscale(2))
                        kpip(3) = anint(rmu(3, i,k) / cellscale(3))
                        rmu(1, i,k) = rmu(1, i,k) - kpip(1) * cellscale(1)
                        rmu(2, i,k) = rmu(2, i,k) - kpip(2) * cellscale(2)
                        rmu(3, i,k) = rmu(3, i,k) - kpip(3) * cellscale(3)
                        r(i,k) = dsqrt(rmu(1, i,k)**2 + rmu(2, i,k)**2 + rmu(3, i,k)**2)
                    enddo
                enddo
            else
                do i = indtmin, indtm
                    do k = 1, nion
                        kpip(1) = anint(rmu(1, i,k) / cellscale(1))
                        kpip(2) = anint(rmu(2, i,k) / cellscale(2))
                        kpip(3) = anint(rmu(3, i,k) / cellscale(3))
                        rmu(1,i,k) = rmu(1,i,k) - kpip(1) * cellscale(1)
                        rmu(2,i,k) = rmu(2,i,k) - kpip(2) * cellscale(2)
                        rmu(3,i,k) = rmu(3,i,k) - kpip(3) * cellscale(3)
                        r(i,k) = dsqrt(rmu(1, i,k)**2 + rmu(2, i,k)**2 + rmu(3, i,k)**2)
                        do kk = 1, 3
                            if(rphase(kk).ne.0.d0.and.2 * (kpip(kk) / 2).ne.kpip(kk)) r(i,k) = -r(i,k)
                        enddo
                    enddo
                enddo
            endif
            ! ********** Periodic Systems End *******
        ELSEIF(LBox.eq.-2.d0) then
            do i = indtmin, indtm
                do k = 1, nion
                r(i,k) = dsqrt(rmu(1, i,k)**2 + rmu(2, i,k)**2 + rmu(3, i,k)**2)
                enddo
            enddo
        ELSE
            DO i = indtmin, indtm
                DO k = 1, nion
                    x2 = DSQRT(rmu(1, i,k)**2 + rmu(2, i,k)**2 + rmu(3, i,k)**2)
                    IF (x2 .LT. mindist) THEN
                        r(i,k) = mindist
                    ELSE
                        r(i,k) = x2
                    END IF
                END DO
            END DO
        END IF
    END IF   ! endif iflagnorm  and |LBox| =/ 3

    if(abs(LBox).eq.3.d0) then
        ! 1 --> pointer to scratch  for makefun
        indp1 = nshell * (indtm + 1) * 20 + 1 ! pointer to cphs
        indp2 = indp1 + 2 * nshell * (indtm + 1) ! pointer to rphs
        indp3 = indp2 + nshell * (indtm + 1)  !pointer to real distp
        if(iespbc) then
!$omp parallel do default(shared) private(i,k,vecscra)
            do k = 1, nion
                 do i = indtmin, indtm
                    rmu(:, i,k) = x(:, 1, i) - rion(:, k)
                    !             vecscra(:)=rmu(:,k,i)
                    !             call CartesianToCrystal(vecscra,1)
                    vecscra(:) = &
                            & car2cry(:, 1) * rmu(1, i,k) + car2cry(:, 2) * rmu(2, i,k) + car2cry(:, 3) * rmu(3, i,k)
                    vecscra(1) = anint(vecscra(1) / cellscale(1))
                    vecscra(2) = anint(vecscra(2) / cellscale(2))
                    vecscra(3) = anint(vecscra(3) / cellscale(3))
                    npip(i,:, k) = vecscra(:)
                    rmu(:, i,k) = rmu(:, i,k)&
                            & - s2r(:, 1) * vecscra(1) - s2r(:, 2) * vecscra(2) - s2r(:, 3) * vecscra(3)
                    !     call dgemv('N',3,3,-1.d0,s2r,3,vecscra,1,1.d0,rmu(1,k,i),1)
                    !             rmu(1,k,i)  = rmu(1,k,i)-npip(1,k,i)*cellscale(1)
                    !             rmu(2,k,i)  = rmu(2,k,i)-npip(2,k,i)*cellscale(2)
                    !             rmu(3,k,i)  = rmu(3,k,i)-npip(3,k,i)*cellscale(3)
                    r(i,k) = max(dsqrt(rmu(1, i,k)**2 + rmu(2, i,k)**2 + &
                            rmu(3, i,k)**2), mindist)
                enddo
            enddo
!$omp end parallel do
        else
!$omp parallel do default(shared) private(i,k)
            do k = 1, nion
                do i = indtmin, indtm
                    rmu(1, i,k) = x(1, 1, i) - rion(1, k)
                    rmu(2, i,k) = x(2, 1, i) - rion(2, k)
                    rmu(3, i,k) = x(3, 1, i) - rion(3, k)
                    r(i,k) = max(dsqrt(rmu(1, i,k)**2 + rmu(2, i,k)**2 + &
                            rmu(3, i,k)**2), mindist)
                enddo
            enddo
!$omp end parallel do
        endif
    endif





    ! SECOND PART TRIVIAL

    IF (lbox .eq. 1.d0) THEN
        rb = 0.d0
        rmub = 0.d0
        rmucosb = 0.d0
        rmusinb = 0.d0
        sinphaseb(:, indtmin:indtm, 1:nion) = 0.d0
        cosphaseb(indtmin:indtm,1:nion) = 0.d0
        cellscalep = 0.d0
        rphasep = 0.d0
!$omp parallel do default(shared) private(i,j,ii,indpar)
        DO j = 1, nion
            DO ii = adr_nion(j), adr_nion(j + 1) - 1
                i = ind_nion(ii)
                if(i.le.nshell) then
                    indpar = max(indpar_tab(i), 0)
             CALL MAKEFUN_PBC_B(ioptorb(i), iocc, indt, i0, indtmin, indtm, &
         & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip, z, zb, &
         &   dd, ddb, r(0,kion(i)), rb(0,kion(i)), rmu(1,0,kion(i)),&
         &  rmub(1,0,kion(i)), distp(dimp * (i - 1) + 1),&
         &  distpb(dimp * (i - 1) + 1), &
         & iflagnorm, cnorm(i), rmucos(1,0,kion(i)),rmucosb(1,0,kion(i)),&
         & rmusin(1,0,kion(i)), rmusinb(1,0,kion(i)),&
         & sinphase(1,0,kion(i)), sinphaseb(1,0,kion(i)),&
         & cosphase(0,kion(i)), cosphaseb(0,kion(i)), cellscale,cellscalep(1, j),rphase, rphasep(1, j))
                endif
            enddo
        END DO
!$omp end parallel do
        do j = 1, nion
            cellscaleb(1) = cellscaleb(1) + cellscalep(1, j)
            cellscaleb(2) = cellscaleb(2) + cellscalep(2, j)
            cellscaleb(3) = cellscaleb(3) + cellscalep(3, j)
            rphaseb(1) = rphaseb(1) + rphasep(1, j)
            rphaseb(2) = rphaseb(2) + rphasep(2, j)
            rphaseb(3) = rphaseb(3) + rphasep(3, j)
        enddo
    ELSEIF(abs(LBox).eq.2.d0) then
        rmub = 0.d0
        rb = 0.d0
!$omp parallel do default(shared) private(i,ii,j,indpar)
        do j = 1, nion
            DO ii = adr_nion(j), adr_nion(j + 1) - 1
                i = ind_nion(ii)
                if(i.le.nshell) then
                    indpar = max(indpar_tab(i), 0)
             CALL MAKEFUN_BUMP_B(ioptorb(i), iocc, indt, i0, indtmin, indtm, &
     & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip, z, zb, dd, &
     & ddb, r(0,kion(i)), rb(0,kion(i)), rmu(1,0,kion(i)), rmub(1,0,kion(i)),&
     & distp(dimp * (i - 1) + 1), distpb(dimp * (i - 1) + 1), &
     & iflagnorm, cnorm(i))
                endif
            END DO
        ENDDO
!$omp end parallel do
    ELSEIF(abs(LBox).eq.3) then


        call makefun_grid_b(distp,distpb,distp(indp1)&
                &, distp(indp2),distp(indp3),distpb(indp3))

    ELSE
        rmub = 0.d0
        rb = 0.d0


!$omp parallel do default(shared) private(i,j,ii,ll,indpar,do_makefun,yeszero_z,jjj)
        DO j = 1, nion
            DO ii = adr_nion(j), adr_nion(j + 1) - 1
                i = ind_nion(ii)
                !   updates rb rmub
                    indpar = max(indpar_tab(i), 0)
          do_makefun=.true.
          if(yes_scemama.and.ioptorb(i).ne.200) then
            jjj=kion(i)
!           do_makefun=.false.
            if(slaterorb_read(i+mshift)) then
              do ll=i0u,indtm
              if(dd(indpar+1)*r(ll,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z(ll)=.false.
              else
              yeszero_z(ll)=.true.
              endif
              enddo
            else
              do ll=i0u,indtm
       if(dd(indpar+1)*r(ll,jjj)*r(ll,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z(ll)=.false.
              else
              yeszero_z(ll)=.true.
              endif
              enddo
            endif
            do_makefun=.not.all(yeszero_z(i0u:indtm))
          endif
               if(do_makefun) then
                    if(yes_scemama.and.indtm.gt.0.and.ioptorb(i).ne.200) then
                     do ll=i0,indtm
                      if(yeszero_z(ll)) then
                      zb(indorb_tab(i)+1:indorb_tab(i+1),ll)=0.d0
                      endif
                     enddo
                     if(yeszero_z(0).and.typecomp.ne.1) then
                     do ll=indt+1,indt+4
                     zb(indorb_tab(i)+1:indorb_tab(i+1),ll)=0.d0
                     enddo
                     endif
                    endif
                    CALL MAKEFUN_B(ioptorb(i), indt, i0, indtmin, indtm, &
           & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip, z, zb&
           &, dd, ddb, zeta, r(0,kion(i)), rb(0,kion(i)),&
           & rmu(1,0,kion(i)), rmub(1,0,kion(i)), distp(dimp * (i - 1) + 1), &
           &  distpb(dimp * (i - 1) + 1),iflagnorm, cnorm(i))
                endif
            enddo
            !  update only indorb indpar indshell
            !      CALL MAKEFUN(ioptorb(i), iocc, indt, 1, 1, 0, &
            !&                1, indpar, indorb, indshell, nelskip, z, dd, zeta&
            !&                , r, rmu, distp, kion(i), nion, iflagnorm , cnorm(i))
        END DO
!$omp end parallel do
    END IF

    ! LEFT this instruction for optimization
!$omp barrier

    if(abs(LBox).eq.3) then

        if(iespbc) then
            do i = indtmin, indtm
                do k = 1, nion
                    rmub(1, i,k) = rmub(1, i,k) + rb(i,k) * rmu(1, i,k) / r(i,k)
                    rmub(2, i,k) = rmub(2, i,k) + rb(i,k) * rmu(2, i,k) / r(i,k)
                    rmub(3, i,k) = rmub(3, i,k) + rb(i,k) * rmu(3, i,k) / r(i,k)
                    rb(i,k) = 0.d0
                    !             The original npip has to be used
                    !            npip(:,k,i) = anint(rmu(:,k,i)/cellscale(:))
                    !    reverse  of
                    !    call dgemv('N',3,3,-1.d0,s2r,3,npip(1,k,i),1,1.d0,rmu(1,k,i),1)
!                   vecscra(:) = npip(:, k, i)
                    !              call dger(3,3,-1.d0,rmub(1,k,i),1,vecscra,1,s2rb,3)
                    do kk=1,3
                    s2rb(kk, 1) = s2rb(kk, 1) - rmub(kk, i,k) * npip(i,1,k)
                    s2rb(kk, 2) = s2rb(kk, 2) - rmub(kk, i,k) * npip(i,2,k)
                    s2rb(kk, 3) = s2rb(kk, 3) - rmub(kk, i,k) * npip(i,3,k)
                    enddo

                    xb(1, 1, i) = xb(1, 1, i) + rmub(1, i,k)
                    xb(2, 1, i) = xb(2, 1, i) + rmub(2, i,k)
                    xb(3, 1, i) = xb(3, 1, i) + rmub(3, i,k)
                    rionb(1, k) = rionb(1, k) - rmub(1, i,k)
                    rionb(2, k) = rionb(2, k) - rmub(2, i,k)
                    rionb(3, k) = rionb(3, k) - rmub(3, i,k)
                    rmub(1, i,k) = 0.d0
                    rmub(2, i,k) = 0.d0
                    rmub(3, i,k) = 0.d0
                enddo
            enddo
        else
            do i = indtmin, indtm
                do k = 1, nion
                    rmub(1, i,k) = rmub(1, i,k) + rb(i,k) * rmu(1, i,k) / r(i,k)
                    rmub(2, i,k) = rmub(2, i,k) + rb(i,k) * rmu(2, i,k) / r(i,k)
                    rmub(3, i,k) = rmub(3, i,k) + rb(i,k) * rmu(3, i,k) / r(i,k)
                    rb(i,k) = 0.d0
                    xb(1, 1, i) = xb(1, 1, i) + rmub(1, i,k)
                    xb(2, 1, i) = xb(2, 1, i) + rmub(2, i,k)
                    xb(3, 1, i) = xb(3, 1, i) + rmub(3, i,k)
                    rionb(1, k) = rionb(1, k) - rmub(1, i,k)
                    rionb(2, k) = rionb(2, k) - rmub(2, i,k)
                    rionb(3, k) = rionb(3, k) - rmub(3, i,k)
                    rmub(1, i,k) = 0.d0
                    rmub(2, i,k) = 0.d0
                    rmub(3, i,k) = 0.d0
                enddo
            enddo
        endif
    else ! Lbox eq 3

        IF(LBox.ne.1.d0) then
            ! given rmub rb r rb  update xb and rionb

            !       do i=indtmin,indtm
            !         do k=1,nion
            !         r(k,i)=max(dsqrt(rmu(1,k,i)**2+rmu(2,k,i)**2+rmu(3,k,i)**2),mindist)
            !         enddo
            !       enddo
            if(abs(LBox).ne.2.d0) then
                DO i = indtmin, indtm
                    DO k = 1, nion
                        IF (r(i,k).gt.mindist) THEN
                            temp2b0 = rb(i,k) / r(i,k)
                            rmub(1, i,k) = rmub(1, i,k) + temp2b0 * rmu(1, i,k)
                            rmub(2, i,k) = rmub(2, i,k) + temp2b0 * rmu(2, i,k)
                            rmub(3, i,k) = rmub(3, i,k) + temp2b0 * rmu(3, i,k)
                        END IF
                    END DO
                END DO
            else

                DO i = indtmin, indtm
                    DO k = 1, nion
                        IF (r(i,k).ne.0.d0) THEN
                            !  Takes into account the change of sign in r for LBox=2 apbc
                            temp2b0 = rb(i,k) / abs(r(i,k))
                            rmub(1, i,k) = rmub(1, i,k) + temp2b0 * rmu(1, i,k)
                            rmub(2, i,k) = rmub(2, i,k) + temp2b0 * rmu(2, i,k)
                            rmub(3, i,k) = rmub(3, i,k) + temp2b0 * rmu(3, i,k)
                        END IF
                        If(LBox.eq.2) then
                            !     recompute npip
                            rmup(1) = x(1, 1, i) - rion(1, k)
                            kpip(1) = anint(rmup(1) / cellscale(1))
                            cellscaleb(1) = cellscaleb(1) - kpip(1) * rmub(1, i,k)

                            rmup(2) = x(2, 1, i) - rion(2, k)
                            kpip(2) = anint(rmup(2) / cellscale(2))
                            cellscaleb(2) = cellscaleb(2) - kpip(2) * rmub(2, i,k)

                            rmup(3) = x(3, 1, i) - rion(3, k)
                            kpip(3) = anint(rmup(3) / cellscale(3))
                            cellscaleb(3) = cellscaleb(3) - kpip(3) * rmub(3, i,k)
                        endif
                    END DO
                END DO

            endif

        ELSE ! LBox

            ! Input sinphaseb cosphaseb rmucosb rmusinb rmub rb output rmub
            do j = indtmin, indtm
                do k = 1, nion
                    if(r(j,k).gt.mindist) rmub(1:3, j,k) = rmub(1:3, j,k)&
                            & + rb(j,k) * rmu(1:3, j,k) / r(j,k)

                    rmusinb(1:3, j,k) = rmusinb(1:3, j,k) + cellpi(1:3) * rmub(1:3, j,k)

                    cellpib(1:3) = cellpib(1:3) + rmub(1:3, j,k) * rmusin(1:3, j,k)

                    rmub(1, j,k) = rmusinb(1, j,k) * rmucos(1, j,k)
                    rmub(2, j,k) = rmusinb(2, j,k) * rmucos(2, j,k)
                    rmub(3, j,k) = rmusinb(3, j,k) * rmucos(3, j,k)
                    rmub(1, j,k) = rmub(1, j,k) - rmucosb(1, j,k) * rmusin(1, j,k)
                    rmub(2, j,k) = rmub(2, j,k) - rmucosb(2, j,k) * rmusin(2, j,k)
                    rmub(3, j,k) = rmub(3, j,k) - rmucosb(3, j,k) * rmusin(3, j,k)

                    rmub(1, j,k) = rmub(1, j,k) / cellpi(1)
                    rmub(2, j,k) = rmub(2, j,k) / cellpi(2)
                    rmub(3, j,k) = rmub(3, j,k) / cellpi(3)
                    !      restoring rmusav
                    rmusav(1) = x(1, 1, j) - rion(1, k)
                    rmusav(2) = x(2, 1, j) - rion(2, k)
                    rmusav(3) = x(3, 1, j) - rion(3, k)
                    cellpib(1) = cellpib(1) - rmub(1, j,k) * rmusav(1) / cellpi(1)
                    cellpib(2) = cellpib(2) - rmub(2, j,k) * rmusav(2) / cellpi(2)
                    cellpib(3) = cellpib(3) - rmub(3, j,k) * rmusav(3) / cellpi(3)
                enddo
            enddo

            if(.not.gamma_point.and..not.iesjas) then
                !   propagation sinphaseb
                ! FIRST PART
                !
                do j = indtmin, indtm
                    do k = 1, nion
                        !          restoring old rmu
                        rmusav(:) = x(:, 1, j) - rion(:, k)

                        rcos(1) = dcos(rphase(1) * rmusav(1))
                        rcos(2) = dcos(rphase(2) * rmusav(2))
                        rcos(3) = dcos(rphase(3) * rmusav(3))
                        rsin(1) = dsin(rphase(1) * rmusav(1))
                        rsin(2) = dsin(rphase(2) * rmusav(2))
                        rsin(3) = dsin(rphase(3) * rmusav(3))

                        rmub(1, j,k) = rmub(1, j,k) + sinphaseb(1, j,k) * (rcos(2) * rcos(3) * &
                                &rcos(1) * rphase(1))&
                                & + sinphaseb(2, j,k) * (-rsin(1) * rsin(2) * rcos(3) * rphase(1))&
                                & + sinphaseb(3, j,k) * (-rsin(1) * rcos(2) * rsin(3) * rphase(1))


                        !   propagation sinphase
                        rsinb(1) = sinphaseb(1, j,k) * rcos(2) * rcos(3)
                        rsinb(2) = sinphaseb(2, j,k) * rcos(1) * rcos(3)
                        rsinb(3) = sinphaseb(3, j,k) * rcos(1) * rcos(2)

                        rcosb(1) = cosphaseb(j,k) * rcos(2) * rcos(3) + sinphaseb(2, j,k) * rsin(2) * rcos(3)&
                                & + sinphaseb(3, j,k) * rsin(3) * rcos(2)
                        !   propagation cosphase

                        rcosb(2) = sinphaseb(3, j,k) * rsin(3) * rcos(1) + sinphaseb(1, j,k) * rsin(1) * rcos(3)&
                                & + cosphaseb(j,k) * rcos(1) * rcos(3)

                        rcosb(3) = sinphaseb(1, j,k) * rsin(1) * rcos(2) + sinphaseb(2, j,k) * rsin(2) * rcos(1)&
                                & + cosphaseb(j,k) * rcos(1) * rcos(2)

                        rmub(2, j,k) = rmub(2, j,k) + sinphaseb(1, j,k) * &
                                & (-rsin(2) * rcos(3) * rsin(1) * rphase(2))&
                                & + sinphaseb(2, j,k) * (rcos(1) * rcos(2) * rcos(3) * rphase(2))&
                                & + sinphaseb(3, j,k) * (-rcos(1) * rsin(2) * rsin(3) * rphase(2))

                        rmub(3, j,k) = rmub(3, j,k) + sinphaseb(1, j,k) * &
                                &(-rcos(2) * rsin(3) * rsin(1) * rphase(3))&
                                & + sinphaseb(2, j,k) * (-rcos(1) * rsin(2) * rsin(3) * rphase(3))&
                                & + sinphaseb(3, j,k) * (rcos(1) * rcos(2) * rcos(3) * rphase(3))

                        !   propagation cosphaseb

                        rmub(1, j,k) = rmub(1, j,k) - cosphaseb(j,k) * sinphase(1, j,k) * rphase(1)
                        rmub(2, j,k) = rmub(2, j,k) - cosphaseb(j,k) * sinphase(2, j,k) * rphase(2)
                        rmub(3, j,k) = rmub(3, j,k) - cosphaseb(j,k) * sinphase(3, j,k) * rphase(3)
                        !propagation cellscaleb given the implicit dependence rphase=cost(:)/cellscale

                        rphaseb(1) = rphaseb(1) - rcosb(1) * rmusav(1) * rsin(1) + rsinb(1) * rmusav(1) * rcos(1)
                        rphaseb(2) = rphaseb(2) - rcosb(2) * rmusav(2) * rsin(2) + rsinb(2) * rmusav(2) * rcos(2)
                        rphaseb(3) = rphaseb(3) - rcosb(3) * rmusav(3) * rsin(3) + rsinb(3) * rmusav(3) * rcos(3)

                    enddo
                enddo

                !  reverse of rphase=2 pi phase/cellscale

                cellscaleb(1) = cellscaleb(1) - rphaseb(1) * rphase(1) / cellscale(1)
                cellscaleb(2) = cellscaleb(2) - rphaseb(2) * rphase(2) / cellscale(2)
                cellscaleb(3) = cellscaleb(3) - rphaseb(3) * rphase(3) / cellscale(3)

            endif ! not gamma_point
            !  reverse of cellpi=cellscale/Pi

            cellscaleb(1) = cellscaleb(1) + cellpib(1) / pi
            cellscaleb(2) = cellscaleb(2) + cellpib(2) / pi
            cellscaleb(3) = cellscaleb(3) + cellpib(3) / pi

        ENDIF ! PBC


        ! Last step common for both

        DO i = indtmin, indtm
            DO k = 1, nion
                xb(3, 1, i) = xb(3, 1, i) + rmub(3, i,k)
                rionb(3, k) = rionb(3, k) - rmub(3, i,k)
                rmub(3, i,k) = 0.0_8
                xb(2, 1, i) = xb(2, 1, i) + rmub(2, i,k)
                rionb(2, k) = rionb(2, k) - rmub(2, i,k)
                rmub(2, i,k) = 0.0_8
                xb(1, 1, i) = xb(1, 1, i) + rmub(1, i,k)
                rionb(1, k) = rionb(1, k) - rmub(1, i,k)
                rmub(1, i,k) = 0.0_8
            END DO
        END DO

    endif ! New case LBox=3
!   The action below has to be  done always
    zb(:, i0:indtm) = 0.d0
    if(typecomp.ne.1) zb(:, indt + 1:indt + 4) = 0.d0

    ! if(iflagnorm.lt.0) then
    !   iflagnorm=-iflagnorm
    !   if(iflagnorm.ne.1) iflagnorm=2
    ! endif

    ! Restore that the tables are for the Jastrow
    !  if(abs(LBox).eq.3) deallocate(psip_r,psip_rb)
    if(iesjas) then
        ! indpar_tab(1)=-1
        LBox = LBox_sav
    endif

  endif ! indtm=0

    RETURN

contains

    subroutine  makefun_grid_b(distp_true, distp_trueb, cphs_r&
            &, rphs, distp, distpb)
        implicit none
        real*8 distp_true(0:indtm, 20, nshell), distp(nelskip, i0:indt + 4)
        real*8 distp_trueb(0:indtm, 20, nshell), distpb(nelskip, i0:indt + 4)
        real*8  cphs_r(2,0:indtm, *)
        real*8 rphs(0:indtm, nshell)
        real*8 cphs
        integer nshift,mshift
        logical yeszero_z(0:indtm)
        rb = 0.d0
        rmub = 0.d0
        rmunewb = 0.d0
        rknewb = 0.d0
        distpb = 0.d0  ! Initialize to zero shared variable.
        if(iespbc)  s2rp = 0.d0

            if(iesjas) then
            nshift=nion
            mshift=nshell_det
            else
            nshift=0
            mshift=0
            endif


!$omp parallel do default(shared)  private(iii,jjj,i,j,ii,jj,kk,ll,indpar,indorb,indshell,kpip,cphs,do_makefun,yeszero_z)
        do jjj = 1, nion
            ! write(6,*) ' rmnew before mekefun -2=',jjj,sum(rmunewb(:,jjj,indtmin:indtm)),sum(rknewb(jjj,indtmin:indtm))
          do ii=1,kgrid_atom(jjj+nshift)%dimshell
            kpip(1) = kgrid_atom(jjj + nshift)%kpip(1, ii)
            kpip(2) = kgrid_atom(jjj + nshift)%kpip(2, ii)
            kpip(3) = kgrid_atom(jjj + nshift)%kpip(3, ii)


            do j = indtmin, indtm
            rmunew(1, j, jjj) = rmu(1, j,jjj) + &
          &s2r(1, 1) * kpip(1) + s2r(1, 2) * kpip(2) + s2r(1, 3) * kpip(3)
            rmunew(2, j, jjj) = rmu(2, j,jjj) + &
          &s2r(2, 1) * kpip(1) + s2r(2, 2) * kpip(2) + s2r(2, 3) * kpip(3)
            rmunew(3, j, jjj) = rmu(3, j,jjj) + &
          &s2r(3, 1) * kpip(1) + s2r(3, 2) * kpip(2) + s2r(3, 3) * kpip(3)
    rknew(j,jjj)=max(dsqrt(rmunew(1,j,jjj)**2+rmunew(2,j,jjj)**2+rmunew(3,j,jjj)**2),mindist)
            enddo

            select case(case_if)

            case(1)
!            if(ipc.eq.2.and..not.iesjas) then ! complex case
                    do ll = i0u, indtm
                        cphs = -(phs(1) * (kpip(1) - npip(ll,1,jjj)) + &
                                phs(2) * (kpip(2) - npip(ll,2,jjj)) + &
                                phs(3) * (kpip(3) - npip(ll,3,jjj)))
                        cphs_r(1,ll,jjj) =dcos(cphs)
                        cphs_r(2,ll,jjj) =-dsin(cphs)
                    enddo
            case(2)
                    do ll = i0u, indtm
            rphs(ll, jjj) = dcos(phs(1) * (kpip(1) - npip(ll,1,jjj)) + &
                                phs(2) * (kpip(2) - npip(ll,2,jjj)) + &
                                phs(3) * (kpip(3) - npip(ll,3,jjj)))
                    enddo
            case(3)
              cphs_r(1,i0u:indtm,jjj) = 1.d0
              cphs_r(2,i0u:indtm,jjj) = 0.d0
            case(4)
              rphs(i0u:indtm,jjj)=1.d0
            end select


            DO iii = adr_nion(jjj), adr_nion(jjj + 1) - 1
                i = ind_nion(iii)
                ! write(6,*) ' rmnew before mekefun -1=',iii,sum(rmunewb(:,jjj,indtmin:indtm)),sum(rknewb(jjj,indtmin:indtm))
!               if(i.le.nshell) then
!                   do ii = 1, kgrid(i + ishift)%dimshell
                 if(kgrid(i+ishift)%tobedone(ii))  then
                    

                        indpar = max(indpar_tab(i), 0)
                        indorb = indorb_tab(i + 1)
                        indshell = indshell_tab(i)

          do_makefun=.true.
          if(yes_scemama.and.ioptorb(i).ne.200) then
!           do_makefun=.false.
            if(slaterorb_read(i+mshift)) then
              do ll=i0u,indtm
              if(dd(indpar+1)*rknew(ll,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z(ll)=.false.
              else
              yeszero_z(ll)=.true.
              endif
              enddo
            else
              do ll=i0u,indtm
       if(dd(indpar+1)*rknew(ll,jjj)*rknew(ll,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z(ll)=.false.
              else
              yeszero_z(ll)=.true.
              endif
              enddo
            endif
            do_makefun=.not.all(yeszero_z(i0u:indtm))
          endif
           if(do_makefun) then 
                    select case(case_upz) 
                    case(1) 
                                do ll = i0, indtm
        call zaxpyr(indorb_tab(i),indorb,cphs_r(1,ll,jjj),zb(1,ll),distpb(1,ll))
                                enddo
                     case(2)
                                do ll = i0, indtm
        call zaxpyr(indorb_tab(i),indorb,cphs_r(1,ll,jjj),zb(1,ll),distpb(1,ll))
                                enddo
                                do ll = indt + 1, indt + 4
        call zaxpyr(indorb_tab(i),indorb,cphs_r(1,0,jjj),zb(1,ll),distpb(1,ll))
                                enddo
                           case(3)
                                do ll = i0, indtm
        call daxrpy(indorb_tab(i),indorb,rphs(ll,jjj),zb(1,ll),distpb(1,ll))
                                enddo
                           case(4)
                                do ll = i0, indtm
        call daxrpy(indorb_tab(i),indorb,rphs(ll,jjj),zb(1,ll),distpb(1,ll))
                                enddo
                                do ll = indt + 1, indt + 4
        call daxrpy(indorb_tab(i),indorb,rphs(0,jjj),zb(1,ll),distpb(1,ll))
                                enddo
                        end select
           if(yes_scemama.and.indtm.gt.0.and.ioptorb(i).ne.200) then
             do ll=i0,indtm
              if(yeszero_z(ll)) then
              distpb(indorb_tab(i)+1:indorb_tab(i+1),ll)=0.d0
              endif
             enddo
             if(yeszero_z(0).and.typecomp.ne.1) then
             do ll=indt+1,indt+4
             distpb(indorb_tab(i)+1:indorb_tab(i+1),ll)=0.d0
             enddo
             endif
            endif
                        CALL MAKEFUN_B(ioptorb(i), indt, i0, indtmin, indtm, &
              & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip&
              &, distp, distpb, dd, ddb, zeta, rknew(0,jjj), rknewb(0,jjj)&
              &, rmunew(1,0,jjj), rmunewb(1,0,jjj), distp_true(0, 1, i)&
              &, distp_trueb(0, 1, i), iflagnorm, cnorm(i))

                  endif ! endif do_makefun

                 endif ! endif to be done

                enddo  ! enddo over the shell of each atom

!  reverse of 
!            do j = indtmin, indtm
!            rmunew(:, j, i_ion) = rmu(:, j,i_ion) + &
!          &s2r(:, 1) * kpip(1) + s2r(:, 2) * kpip(2) + s2r(:, 3) * kpip(3)
!    rknew(j,i_ion)=max(dsqrt(rmunew(1,j,i_ion)**2+rmunew(2,j,i_ion)**2+rmunew(3,j,i_ion)**2),mindist)
!            enddo

! Below we do not  put the if rknew>mindist because it is negligible case that
! occurs  rarely and in case  the  correction is negligible instead of being zero.
                        do j = indtmin, indtm

                            rmunewb(1, j,jjj) = rmunewb(1, j,jjj) + rknewb(j,jjj) * &
                                    &rmunew(1, j,jjj) / rknew(j,jjj)
                            rmunewb(2, j,jjj) = rmunewb(2, j,jjj) + rknewb(j,jjj) * &
                                    &rmunew(2, j,jjj) / rknew(j,jjj)
                            rmunewb(3, j,jjj) = rmunewb(3, j,jjj) + rknewb(j,jjj) * &
                                    &rmunew(3, j,jjj) / rknew(j,jjj)
                            rknewb(j,jjj) = 0.d0

! Unrolled by  hand.

             s2rp(1, 1, jjj) = s2rp(1, 1, jjj) + rmunewb(1, j,jjj) * kpip(1)
             s2rp(1, 2, jjj) = s2rp(1, 2, jjj) + rmunewb(1, j,jjj) * kpip(2)
             s2rp(1, 3, jjj) = s2rp(1, 3, jjj) + rmunewb(1, j,jjj) * kpip(3)

             s2rp(2, 1, jjj) = s2rp(2, 1, jjj) + rmunewb(2, j,jjj) * kpip(1)
             s2rp(2, 2, jjj) = s2rp(2, 2, jjj) + rmunewb(2, j,jjj) * kpip(2)
             s2rp(2, 3, jjj) = s2rp(2, 3, jjj) + rmunewb(2, j,jjj) * kpip(3)

             s2rp(3, 1, jjj) = s2rp(3, 1, jjj) + rmunewb(3, j,jjj) * kpip(1)
             s2rp(3, 2, jjj) = s2rp(3, 2, jjj) + rmunewb(3, j,jjj) * kpip(2)
             s2rp(3, 3, jjj) = s2rp(3, 3, jjj) + rmunewb(3, j,jjj) * kpip(3)

                            rmub(1, j,jjj) = rmub(1, j,jjj) + rmunewb(1, j,jjj)
                            rmub(2, j,jjj) = rmub(2, j,jjj) + rmunewb(2, j,jjj)
                            rmub(3, j,jjj) = rmub(3, j,jjj) + rmunewb(3, j,jjj)
                            rmunewb(1, j,jjj) = 0.d0
                            rmunewb(2, j,jjj) = 0.d0
                            rmunewb(3, j,jjj) = 0.d0
                           enddo ! end  j
               enddo ! enddo over the grid ii 
        enddo  ! enddo over jjj
!$omp end parallel do
        !  To avoid an inefficient reduction.
        if(iespbc)  then
            do j = 1, nion
                s2rb(:, :) = s2rb(:, :) + s2rp(:, :, j)
            enddo
        endif
    end subroutine makefun_grid_b
END SUBROUTINE UPNEWWF_B

SUBROUTINE UPNEWWF0_B(indt, typecomp, nshell, ioptorb, iocc, x&
        &, xb, nel, r, rb, rmu, rmub, dd, ddb, zeta, rion, rionb, distp, distpb&
        &, z, zb, nelskip, nion, kion, iflagnorm, cnorm, lbox, rmucos, rmucosb, &
        &  rmusin, rmusinb, mindist, cellscale, cellscaleb, s2r, s2rb&
        &, indpar_tab, indorb_tab, indshell_tab, adr_nion, ind_nion)
    use allio, only :yesupel,rank,ikshift,kgrid,kgrid_atom,rmunew,rknew,npip&
            !                 &,rmunewb,rknewb,iespbc,yesupel,nelorbh,indt
    &,rmunewb, rknewb, iespbc, yesupel, yes_crystalj,lepsbas,yes_scemama,slaterorb_read,nshell_det
    use Cell, ONLY : cellpi, cellscalep, rphasep, rphase, cosphase, cosphaseb, sinphase, &
            &  sinphaseb, gamma_point, phase2pi, phase2pi_down, at, s2rp, car2cry
    use Constants, only : ipc, zzero, pi, ip4

    IMPLICIT NONE
    INTEGER :: nel, indorb, indpar, indshell, iocc(*), ioptorb(*), j, i &
            &, indt, nelskip, nshell, nion, k, kion(*), iflagnorm &
            &, kk, ii, typecomp, kpip(3), ll, jj, iii, jjj, ishift&
            &, dimp, indp1, indp2, indp3,mshift
               REAL*8 :: dd(*), ddb(*), rmu(3, 0:0,nion), rion(3, nion)&
            &, distp(23 * nshell + nelskip * (indt + 5 ))&
            &, cnorm(*), zeta(*), r(0:0,nion), rcos_sav(3) &
            &, mindist, rcos(3), rsin(3), rcosb(3), rsinb(3), rphs(0:0)
    INTEGER :: indpar_tab(*), indorb_tab(*), indshell_tab(*), adr_nion(*), ind_nion(*)
    REAL*8 :: rmub(3, 0:0,nion), rionb(3, nion)&
            &, distpb(23 * nshell + nelskip * (indt + 5 ))&
            &, zb(nelskip * min(ipc, indpar_tab(1) + 2), 0:*)&
            &, z(nelskip * min(ipc, indpar_tab(1) + 2), 0:*)&
            &, rb(0:0,nion), rmusav(3), cellscale(3), cellscaleb(3)&
            &, rmup(3), s2r(3, 3), s2rb(3, 3)
    REAL*8 :: x(3, nel, 0:0)
    REAL*8 :: xb(3, nel, 0:0)
    REAL*8 :: lbox, lbox_sav, rmucos(3, 0:0,nion), rmusin(3, 0:0,nion)
    REAL*8 :: rmucosb(3, 0:0,nion), rmusinb(3, 0:0,nion)
    REAL*8 :: rphaseb(3), vecscra(3)
    !  real*8, dimension(:,:), allocatable:: psip_r,psip_rb
    !  real*8 :: psip_rb(nelorbh,0:indt+ip4),psip_r(nelorbh,0:indt+ip4)
    DOUBLE PRECISION :: x1, x2, cellpib(3)
    DOUBLE PRECISION :: temp2b0
    logical iesjas,gammaorj
    real*8 phs(3)
    integer case_if,case_upz
    logical do_makefun
    logical yeszero_z


    if(indpar_tab(1).eq.-1) then
        iesjas = .true.
        !  indpar_tab(1)=0
        mshift=nshell_det
    else
        iesjas = .false.
        mshift=0
    endif

    LBox_sav = LBox
    ishift = 0
    dimp = 20
    if(LBox.eq.3) then
        if(yesupel) then
            phs(:) = phase2pi(:)
        else
            phs(:) = phase2pi_down(:)
        endif
    else
        phs = 0.d0
    endif
    if(abs(LBox).eq.3.and.iesjas) then
        ! always use the old periodic basis for the Jastrow
        if(.not.yes_crystalj) then
            if(iespbc) then
                !      iflagnorm=2
                LBox = 1.d0
            else
                LBox = -1.d0
            endif
        else
            ishift = ikshift
        endif
        phs(:) = 0.d0
        !  elseif(abs(LBox).eq.3) then
        !  allocate(psip_r(nelskip,i0:indt+4),psip_rb(nelskip,i0:indt+4))
    endif

            if(sum(abs(phs(:))).eq.0.d0) then
            gammaorj=.true.
            else
            gammaorj=.false.
            endif
            if(ipc.eq.2.and..not.gammaorj) then ! complex case
            case_if=1
            elseif(.not.gammaorj.and.ipc.eq.1) then
            case_if=2
            else if(ipc.eq.2.and..not.iesjas) then
            case_if=3
            else
            case_if=4
            endif

            if(ipc.eq.2.and..not.iesjas) then ! complex case
                if(typecomp.eq.1) then
                case_upz=1
                else
                case_upz=2
                endif
            else ! real case or  Jastrow
                if(typecomp.eq.1) then
                case_upz=3
                else
                case_upz=4
                endif
            endif

    !  if(ishift.ne.0) write(6,*) ' dimension work3 inside =',nelskip*(indt+5+(20*(indt+1))/nelskip-i0+1)

    indshell = 0
    indorb = 0
    indpar = 0
    cellpib = 0.d0
    rphaseb = 0.d0
    ! used only here

    if(abs(LBox).ne.3) then
            DO k = 1, nion
                rmu(1, 0,k) = x(1, 1, 0) - rion(1, k)
                rmu(2, 0,k) = x(2, 1, 0) - rion(2, k)
                rmu(3, 0,k) = x(3, 1, 0) - rion(3, k)
            END DO
        ! ****** Periodic Systems *************
        IF (lbox .eq. 1.d0) THEN
            IF (.NOT.gamma_point.and..not.iesjas) THEN
                    DO k = 1, nion
                        rcos_sav(1) = DCOS(rphase(1) * rmu(1, 0,k))
                        rcos_sav(2) = DCOS(rphase(2) * rmu(2, 0,k))
                        rcos_sav(3) = DCOS(rphase(3) * rmu(3, 0,k))
                        cosphase(0,k) = rcos_sav(1) * rcos_sav(2) * rcos_sav(3)
                        sinphase(1, 0,k) = DSIN(rphase(1) * rmu(1, 0,k)) * rcos_sav(2)&
                                & * rcos_sav(3)
                        sinphase(2, 0,k) = DSIN(rphase(2) * rmu(2, 0,k)) * rcos_sav(1)&
                                & * rcos_sav(3)
                        sinphase(3, 0,k) = DSIN(rphase(3) * rmu(3, 0,k)) * rcos_sav(1)&
                                & * rcos_sav(2)
                    END DO
            END IF
                DO k = 1, nion
                    rmu(1, 0,k) = rmu(1, 0,k) / cellpi(1)
                    rmu(2, 0,k) = rmu(2, 0,k) / cellpi(2)
                    rmu(3, 0,k) = rmu(3, 0,k) / cellpi(3)
                    rmucos(1, 0,k) = DCOS(rmu(1, 0,k))
                    rmucos(2, 0,k) = DCOS(rmu(2, 0,k))
                    rmucos(3, 0,k) = DCOS(rmu(3, 0,k))
                    rmusin(1, 0,k) = DSIN(rmu(1, 0,k))
                    rmusin(2, 0,k) = DSIN(rmu(2, 0,k))
                    rmusin(3, 0,k) = DSIN(rmu(3, 0,k))
                    rmu(1, 0,k) = cellpi(1) * rmusin(1, 0,k)
                    rmu(2, 0,k) = cellpi(2) * rmusin(2, 0,k)
                    rmu(3, 0,k) = cellpi(3) * rmusin(3, 0,k)
                    x1 = DSQRT(rmu(1, 0,k)**2 + rmu(2, 0,k)**2 + rmu(3, 0,k)**2&
                            &)
                    IF (x1 .LT. mindist) THEN
                        r(0,k) = mindist
                    ELSE
                        r(0,k) = x1
                    END IF
                END DO
        ELSEIF(LBox.eq.2) then
            cellpib = 0.d0
            rphaseb = 0.d0
            if(gamma_point.or.iesjas) then
                    do k = 1, nion
                        kpip(1) = anint(rmu(1, 0,k) / cellscale(1))
                        kpip(2) = anint(rmu(2, 0,k) / cellscale(2))
                        kpip(3) = anint(rmu(3, 0,k) / cellscale(3))
                        rmu(1, 0,k) = rmu(1, 0,k) - kpip(1) * cellscale(1)
                        rmu(2, 0,k) = rmu(2, 0,k) - kpip(2) * cellscale(2)
                        rmu(3, 0,k) = rmu(3, 0,k) - kpip(3) * cellscale(3)
                        r(0,k) = dsqrt(rmu(1, 0,k)**2 + rmu(2, 0,k)**2 + rmu(3, 0,k)**2)
                    enddo
            else
                    do k = 1, nion
                        kpip(1) = anint(rmu(1, 0,k) / cellscale(1))
                        kpip(2) = anint(rmu(2, 0,k) / cellscale(2))
                        kpip(3) = anint(rmu(3, 0,k) / cellscale(3))
                        rmu(1, 0,k) = rmu(1, 0,k) - kpip(1) * cellscale(1)
                        rmu(2, 0,k) = rmu(2, 0,k) - kpip(2) * cellscale(2)
                        rmu(3, 0,k) = rmu(3, 0,k) - kpip(3) * cellscale(3)
                        r(0,k) = dsqrt(rmu(1, 0,k)**2 + rmu(2, 0,k)**2 + rmu(3, i,k)**2)
                        do kk = 1, 3
                            if(rphase(kk).ne.0.d0.and.2 * (kpip(kk) / 2).ne.kpip(kk)) r(0,k) = -r(0,k)
                        enddo
                    enddo
            endif
            ! ********** Periodic Systems End *******
        ELSEIF(LBox.eq.-2.d0) then
                do k = 1, nion
                r(0,k) = dsqrt(rmu(1, 0,k)**2 + rmu(2, 0,k)**2 + rmu(3, 0,k)**2)
                enddo
        ELSE
                DO k = 1, nion
                    x2 = DSQRT(rmu(1, 0,k)**2 + rmu(2, 0,k)**2 + rmu(3, 0,k)**2)
                    IF (x2 .LT. mindist) THEN
                        r(0,k) = mindist
                    ELSE
                        r(0,k) = x2
                    END IF
                END DO
        END IF
    END IF   ! endif iflagnorm  and |LBox| =/ 3

    if(abs(LBox).eq.3.d0) then
        ! 1 --> pointer to scratch  for makefun
        indp1 = nshell * 20 + 1 ! pointer to cphs
        indp2 = indp1 + 2 * nshell  ! pointer to rphs
        indp3 = indp2 + nshell  !pointer to real distp
        if(iespbc) then
!$omp parallel do default(shared) private(k,vecscra)
            do k = 1, nion
                    rmu(:, 0,k) = x(:, 1, 0) - rion(:, k)
                    !             vecscra(:)=rmu(:,k,i)
                    !             call CartesianToCrystal(vecscra,1)
                    vecscra(:) = &
                            & car2cry(:, 1) * rmu(1, 0,k) + car2cry(:, 2) * rmu(2, 0,k) + car2cry(:, 3) * rmu(3, 0,k)
                    vecscra(1) = anint(vecscra(1) / cellscale(1))
                    vecscra(2) = anint(vecscra(2) / cellscale(2))
                    vecscra(3) = anint(vecscra(3) / cellscale(3))
                    npip(0,:, k) = vecscra(:)
                    rmu(:, 0,k) = rmu(:, 0,k)&
                            & - s2r(:, 1) * vecscra(1) - s2r(:, 2) * vecscra(2) - s2r(:, 3) * vecscra(3)
                    !     call dgemv('N',3,3,-1.d0,s2r,3,vecscra,1,1.d0,rmu(1,k,i),1)
                    !             rmu(1,k,i)  = rmu(1,k,i)-npip(1,k,i)*cellscale(1)
                    !             rmu(2,k,i)  = rmu(2,k,i)-npip(2,k,i)*cellscale(2)
                    !             rmu(3,k,i)  = rmu(3,k,i)-npip(3,k,i)*cellscale(3)
                    r(0,k) = max(dsqrt(rmu(1, 0,k)**2 + rmu(2, 0,k)**2 + &
                            rmu(3, 0,k)**2), mindist)
            enddo
!$omp end parallel do
        else
!$omp parallel do default(shared) private(k)
            do k = 1, nion
                    rmu(1, 0,k) = x(1, 1, 0) - rion(1, k)
                    rmu(2, 0,k) = x(2, 1, 0) - rion(2, k)
                    rmu(3, 0,k) = x(3, 1, 0) - rion(3, k)
                    r(0,k) = max(dsqrt(rmu(1, 0,k)**2 + rmu(2, 0,k)**2 + &
                            rmu(3, 0,k)**2), mindist)
            enddo
!$omp end parallel do
        endif
    endif





    ! SECOND PART TRIVIAL

    IF (lbox .eq. 1.d0) THEN
        rb = 0.d0
        rmub = 0.d0
        rmucosb = 0.d0
        rmusinb = 0.d0
        sinphaseb(:,0, 1:nion) = 0.d0
        cosphaseb(0,1:nion) = 0.d0
        cellscalep = 0.d0
        rphasep = 0.d0
!$omp parallel do default(shared) private(i,j,ii,indpar)
        DO j = 1, nion
            DO ii = adr_nion(j), adr_nion(j + 1) - 1
                i = ind_nion(ii)
                if(i.le.nshell) then
                    indpar = max(indpar_tab(i), 0)
             CALL MAKEFUN0_PBC_B(ioptorb(i), iocc, indt,&
         & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip, z, zb, &
         &   dd, ddb, r(0,kion(i)), rb(0,kion(i)), rmu(1,0,kion(i)),&
         &  rmub(1,0,kion(i)), distp(dimp * (i - 1) + 1),&
         &  distpb(dimp * (i - 1) + 1),&
         & iflagnorm, cnorm(i), rmucos(1,0,kion(i)),rmucosb(1,0,kion(i)),&
         & rmusin(1,0,kion(i)), rmusinb(1,0,kion(i)),sinphase(1,0,kion(i)), sinphaseb(1,0,kion(i)),&
         & cosphase(0,kion(i)),cosphaseb(0,kion(i)),cellscale,cellscalep(1, j),rphase,rphasep(1, j))
                endif
            enddo
        END DO
!$omp end parallel do
        do j = 1, nion
            cellscaleb(1) = cellscaleb(1) + cellscalep(1, j)
            cellscaleb(2) = cellscaleb(2) + cellscalep(2, j)
            cellscaleb(3) = cellscaleb(3) + cellscalep(3, j)
            rphaseb(1) = rphaseb(1) + rphasep(1, j)
            rphaseb(2) = rphaseb(2) + rphasep(2, j)
            rphaseb(3) = rphaseb(3) + rphasep(3, j)
        enddo
    ELSEIF(abs(LBox).eq.2.d0) then
        rmub = 0.d0
        rb = 0.d0
!$omp parallel do default(shared) private(i,ii,j,indpar)
        do j = 1, nion
            DO ii = adr_nion(j), adr_nion(j + 1) - 1
                i = ind_nion(ii)
                if(i.le.nshell) then
                    indpar = max(indpar_tab(i), 0)
             CALL MAKEFUN0_BUMP_B(ioptorb(i), iocc, indt, &
     & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip, z, zb, dd, &
     & ddb, r(0,kion(i)), rb(0,kion(i)), rmu(1,0,kion(i)), rmub(1,0,kion(i)),&
     & distp(dimp * (i - 1) + 1), distpb(dimp * (i - 1) + 1), &
     & iflagnorm, cnorm(i))
                endif
            END DO
        ENDDO
!$omp end parallel do
    ELSEIF(abs(LBox).eq.3) then


        call makefun_grid_b(distp,distpb,distp(indp1)&
                &, distp(indp2),distp(indp3),distpb(indp3))

    ELSE
        rmub = 0.d0
        rb = 0.d0


!$omp parallel do default(shared) private(i,j,ii,indpar,do_makefun,yeszero_z,jjj)
        DO j = 1, nion
            DO ii = adr_nion(j), adr_nion(j + 1) - 1
                i = ind_nion(ii)
                !   updates rb rmub
                    indpar = max(indpar_tab(i), 0)
          do_makefun=.true.
          if(yes_scemama.and.ioptorb(i).ne.200) then
            jjj=kion(i)
!           do_makefun=.false.
            if(slaterorb_read(i+mshift)) then
              if(dd(indpar+1)*r(0,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z=.false.
              else
              yeszero_z=.true.
              endif
            else
       if(dd(indpar+1)*r(0,jjj)*r(0,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z=.false.
              else
              yeszero_z=.true.
              endif
            endif
            do_makefun=.not.yeszero_z
          endif
               if(do_makefun) then
                    CALL MAKEFUN0_B(ioptorb(i), indt,&
           & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip, z, zb&
           &, dd, ddb, zeta, r(0,kion(i)), rb(0,kion(i)),&
           & rmu(1,0,kion(i)), rmub(1,0,kion(i)), distp(dimp * (i - 1) + 1), &
           &  distpb(dimp * (i - 1) + 1),iflagnorm, cnorm(i))
                endif
            enddo
            !  update only indorb indpar indshell
            !      CALL MAKEFUN(ioptorb(i), iocc, indt, 1, 1, 0, &
            !&                1, indpar, indorb, indshell, nelskip, z, dd, zeta&
            !&                , r, rmu, distp, kion(i), nion, iflagnorm , cnorm(i))
        END DO
!$omp end parallel do
    END IF

    ! LEFT this instruction for optimization
!$omp barrier

    if(abs(LBox).eq.3) then

        if(iespbc) then
                do k = 1, nion
                    rmub(1, 0,k) = rmub(1, 0,k) + rb(0,k) * rmu(1, 0,k) / r(0,k)
                    rmub(2, 0,k) = rmub(2, 0,k) + rb(0,k) * rmu(2, 0,k) / r(0,k)
                    rmub(3, 0,k) = rmub(3, 0,k) + rb(0,k) * rmu(3, 0,k) / r(0,k)
                    rb(0,k) = 0.d0
                    !             The original npip has to be used
                    !            npip(:,k,i) = anint(rmu(:,k,i)/cellscale(:))
                    !    reverse  of
                    !    call dgemv('N',3,3,-1.d0,s2r,3,npip(1,k,i),1,1.d0,rmu(1,k,i),1)
!                   vecscra(:) = npip(:, k, i)
                    !              call dger(3,3,-1.d0,rmub(1,k,i),1,vecscra,1,s2rb,3)
                    do kk=1,3
                    s2rb(kk, 1) = s2rb(kk, 1) - rmub(kk, 0,k) * npip(0,1,k)
                    s2rb(kk, 2) = s2rb(kk, 2) - rmub(kk, 0,k) * npip(0,2,k)
                    s2rb(kk, 3) = s2rb(kk, 3) - rmub(kk, 0,k) * npip(0,3,k)
                    enddo

                    xb(1, 1, 0) = xb(1, 1, 0) + rmub(1, 0,k)
                    xb(2, 1, 0) = xb(2, 1, 0) + rmub(2, 0,k)
                    xb(3, 1, 0) = xb(3, 1, 0) + rmub(3, 0,k)
                    rionb(1, k) = rionb(1, k) - rmub(1, 0,k)
                    rionb(2, k) = rionb(2, k) - rmub(2, 0,k)
                    rionb(3, k) = rionb(3, k) - rmub(3, 0,k)
                    rmub(1, 0,k) = 0.d0
                    rmub(2, 0,k) = 0.d0
                    rmub(3, 0,k) = 0.d0
                enddo
        else
                do k = 1, nion
                    rmub(1, 0,k) = rmub(1, 0,k) + rb(0,k) * rmu(1, 0,k) / r(0,k)
                    rmub(2, 0,k) = rmub(2, 0,k) + rb(0,k) * rmu(2, 0,k) / r(0,k)
                    rmub(3, 0,k) = rmub(3, 0,k) + rb(0,k) * rmu(3, 0,k) / r(0,k)
                    rb(i,k) = 0.d0
                    xb(1, 1, 0) = xb(1, 1, 0) + rmub(1, 0,k)
                    xb(2, 1, 0) = xb(2, 1, 0) + rmub(2, 0,k)
                    xb(3, 1, 0) = xb(3, 1, 0) + rmub(3, 0,k)
                    rionb(1, k) = rionb(1, k) - rmub(1, 0,k)
                    rionb(2, k) = rionb(2, k) - rmub(2, 0,k)
                    rionb(3, k) = rionb(3, k) - rmub(3, 0,k)
                    rmub(1, 0,k) = 0.d0
                    rmub(2, 0,k) = 0.d0
                    rmub(3, 0,k) = 0.d0
                enddo
        endif
    else ! Lbox eq 3

        IF(LBox.ne.1.d0) then
            ! given rmub rb r rb  update xb and rionb

            !       do i=indtmin,indtm
            !         do k=1,nion
            !         r(k,i)=max(dsqrt(rmu(1,k,i)**2+rmu(2,k,i)**2+rmu(3,k,i)**2),mindist)
            !         enddo
            !       enddo
            if(abs(LBox).ne.2.d0) then
                    DO k = 1, nion
                        IF (r(0,k).gt.mindist) THEN
                            temp2b0 = rb(0,k) / r(0,k)
                            rmub(1, 0,k) = rmub(1, 0,k) + temp2b0 * rmu(1, 0,k)
                            rmub(2, 0,k) = rmub(2, 0,k) + temp2b0 * rmu(2, 0,k)
                            rmub(3, 0,k) = rmub(3, 0,k) + temp2b0 * rmu(3, 0,k)
                        END IF
                    END DO
            else

                    DO k = 1, nion
                        IF (r(0,k).ne.0.d0) THEN
                            !  Takes into account the change of sign in r for LBox=2 apbc
                            temp2b0 = rb(0,k) / abs(r(0,k))
                            rmub(1, 0,k) = rmub(1, 0,k) + temp2b0 * rmu(1, 0,k)
                            rmub(2, 0,k) = rmub(2, 0,k) + temp2b0 * rmu(2, 0,k)
                            rmub(3, 0,k) = rmub(3, 0,k) + temp2b0 * rmu(3, 0,k)
                        END IF
                        If(LBox.eq.2) then
                            !     recompute npip
                            rmup(1) = x(1, 1, 0) - rion(1, k)
                            kpip(1) = anint(rmup(1) / cellscale(1))
                            cellscaleb(1) = cellscaleb(1) - kpip(1) * rmub(1, 0,k)

                            rmup(2) = x(2, 1, 0) - rion(2, k)
                            kpip(2) = anint(rmup(2) / cellscale(2))
                            cellscaleb(2) = cellscaleb(2) - kpip(2) * rmub(2, 0,k)

                            rmup(3) = x(3, 1, 0) - rion(3, k)
                            kpip(3) = anint(rmup(3) / cellscale(3))
                            cellscaleb(3) = cellscaleb(3) - kpip(3) * rmub(3, 0,k)
                        endif
                    END DO

            endif

        ELSE ! LBox

            ! Input sinphaseb cosphaseb rmucosb rmusinb rmub rb output rmub
                do k = 1, nion
                    if(r(0,k).gt.mindist) rmub(1:3, 0,k) = rmub(1:3, 0,k)&
                            & + rb(0,k) * rmu(1:3, 0,k) / r(0,k)

                    rmusinb(1:3, 0,k) = rmusinb(1:3, 0,k) + cellpi(1:3) * rmub(1:3, 0,k)

                    cellpib(1:3) = cellpib(1:3) + rmub(1:3, 0,k) * rmusin(1:3, 0,k)

                    rmub(1, 0,k) = rmusinb(1, 0,k) * rmucos(1, 0,k)
                    rmub(2, 0,k) = rmusinb(2, 0,k) * rmucos(2, 0,k)
                    rmub(3, 0,k) = rmusinb(3, 0,k) * rmucos(3, 0,k)
                    rmub(1, 0,k) = rmub(1, 0,k) - rmucosb(1, 0,k) * rmusin(1, 0,k)
                    rmub(2, 0,k) = rmub(2, 0,k) - rmucosb(2, 0,k) * rmusin(2, 0,k)
                    rmub(3, 0,k) = rmub(3, 0,k) - rmucosb(3, 0,k) * rmusin(3, 0,k)

                    rmub(1, 0,k) = rmub(1, 0,k) / cellpi(1)
                    rmub(2, 0,k) = rmub(2, 0,k) / cellpi(2)
                    rmub(3, 0,k) = rmub(3, 0,k) / cellpi(3)
                    !      restoring rmusav
                    rmusav(1) = x(1, 1, 0) - rion(1, k)
                    rmusav(2) = x(2, 1, 0) - rion(2, k)
                    rmusav(3) = x(3, 1, 0) - rion(3, k)
                    cellpib(1) = cellpib(1) - rmub(1, 0,k) * rmusav(1) / cellpi(1)
                    cellpib(2) = cellpib(2) - rmub(2, 0,k) * rmusav(2) / cellpi(2)
                    cellpib(3) = cellpib(3) - rmub(3, 0,k) * rmusav(3) / cellpi(3)
                enddo

            if(.not.gamma_point.and..not.iesjas) then
                !   propagation sinphaseb
                ! FIRST PART
                !
                    do k = 1, nion
                        !          restoring old rmu
                        rmusav(:) = x(:, 1, 0) - rion(:, k)

                        rcos(1) = dcos(rphase(1) * rmusav(1))
                        rcos(2) = dcos(rphase(2) * rmusav(2))
                        rcos(3) = dcos(rphase(3) * rmusav(3))
                        rsin(1) = dsin(rphase(1) * rmusav(1))
                        rsin(2) = dsin(rphase(2) * rmusav(2))
                        rsin(3) = dsin(rphase(3) * rmusav(3))

                        rmub(1, 0,k) = rmub(1, 0,k) + sinphaseb(1, 0,k) * (rcos(2) * rcos(3) * &
                                &rcos(1) * rphase(1))&
                                & + sinphaseb(2, 0,k) * (-rsin(1) * rsin(2) * rcos(3) * rphase(1))&
                                & + sinphaseb(3, 0,k) * (-rsin(1) * rcos(2) * rsin(3) * rphase(1))


                        !   propagation sinphase
                        rsinb(1) = sinphaseb(1, 0,k) * rcos(2) * rcos(3)
                        rsinb(2) = sinphaseb(2, 0,k) * rcos(1) * rcos(3)
                        rsinb(3) = sinphaseb(3, 0,k) * rcos(1) * rcos(2)

                        rcosb(1) = cosphaseb(0,k) * rcos(2) * rcos(3) + sinphaseb(2, 0,k) * rsin(2) * rcos(3)&
                                & + sinphaseb(3, 0,k) * rsin(3) * rcos(2)
                        !   propagation cosphase

                        rcosb(2) = sinphaseb(3, 0,k) * rsin(3) * rcos(1) + sinphaseb(1, 0,k) * rsin(1) * rcos(3)&
                                & + cosphaseb(0,k) * rcos(1) * rcos(3)

                        rcosb(3) = sinphaseb(1, 0,k) * rsin(1) * rcos(2) + sinphaseb(2, 0,k) * rsin(2) * rcos(1)&
                                & + cosphaseb(0,k) * rcos(1) * rcos(2)

                        rmub(2, 0,k) = rmub(2, 0,k) + sinphaseb(1, 0,k) * &
                                & (-rsin(2) * rcos(3) * rsin(1) * rphase(2))&
                                & + sinphaseb(2, 0,k) * (rcos(1) * rcos(2) * rcos(3) * rphase(2))&
                                & + sinphaseb(3, 0,k) * (-rcos(1) * rsin(2) * rsin(3) * rphase(2))

                        rmub(3, 0,k) = rmub(3, 0,k) + sinphaseb(1, 0,k) * &
                                &(-rcos(2) * rsin(3) * rsin(1) * rphase(3))&
                                & + sinphaseb(2, 0,k) * (-rcos(1) * rsin(2) * rsin(3) * rphase(3))&
                                & + sinphaseb(3, 0,k) * (rcos(1) * rcos(2) * rcos(3) * rphase(3))

                        !   propagation cosphaseb

                        rmub(1, 0,k) = rmub(1, 0,k) - cosphaseb(0,k) * sinphase(1, 0,k) * rphase(1)
                        rmub(2, 0,k) = rmub(2, 0,k) - cosphaseb(0,k) * sinphase(2, 0,k) * rphase(2)
                        rmub(3, 0,k) = rmub(3, 0,k) - cosphaseb(0,k) * sinphase(3, 0,k) * rphase(3)
                        !propagation cellscaleb given the implicit dependence rphase=cost(:)/cellscale

                        rphaseb(1) = rphaseb(1) - rcosb(1) * rmusav(1) * rsin(1) + rsinb(1) * rmusav(1) * rcos(1)
                        rphaseb(2) = rphaseb(2) - rcosb(2) * rmusav(2) * rsin(2) + rsinb(2) * rmusav(2) * rcos(2)
                        rphaseb(3) = rphaseb(3) - rcosb(3) * rmusav(3) * rsin(3) + rsinb(3) * rmusav(3) * rcos(3)

                    enddo

                !  reverse of rphase=2 pi phase/cellscale

                cellscaleb(1) = cellscaleb(1) - rphaseb(1) * rphase(1) / cellscale(1)
                cellscaleb(2) = cellscaleb(2) - rphaseb(2) * rphase(2) / cellscale(2)
                cellscaleb(3) = cellscaleb(3) - rphaseb(3) * rphase(3) / cellscale(3)

            endif ! not gamma_point
            !  reverse of cellpi=cellscale/Pi

            cellscaleb(1) = cellscaleb(1) + cellpib(1) / pi
            cellscaleb(2) = cellscaleb(2) + cellpib(2) / pi
            cellscaleb(3) = cellscaleb(3) + cellpib(3) / pi

        ENDIF ! PBC


        ! Last step common for both

            DO k = 1, nion
                xb(3, 1, 0) = xb(3, 1, 0) + rmub(3, 0,k)
                rionb(3, k) = rionb(3, k) - rmub(3, 0,k)
                rmub(3, 0,k) = 0.0_8
                xb(2, 1, 0) = xb(2, 1, 0) + rmub(2, 0,k)
                rionb(2, k) = rionb(2, k) - rmub(2, 0,k)
                rmub(2, 0,k) = 0.0_8
                xb(1, 1, 0) = xb(1, 1, 0) + rmub(1, 0,k)
                rionb(1, k) = rionb(1, k) - rmub(1, 0,k)
                rmub(1, 0,k) = 0.0_8
            END DO

    endif ! New case LBox=3
!   The action below has to be  done always
    zb(:, 0) = 0.d0
    if(typecomp.ne.1) zb(:, indt + 1:indt + 4) = 0.d0

    ! if(iflagnorm.lt.0) then
    !   iflagnorm=-iflagnorm
    !   if(iflagnorm.ne.1) iflagnorm=2
    ! endif

    ! Restore that the tables are for the Jastrow
    !  if(abs(LBox).eq.3) deallocate(psip_r,psip_rb)
    if(iesjas) then
        ! indpar_tab(1)=-1
        LBox = LBox_sav
    endif

    RETURN

contains

    subroutine  makefun_grid_b(distp_true, distp_trueb, cphs_r&
            &, rphs, distp, distpb)
        implicit none
        real*8 distp_true(0:0, 20, nshell), distp(nelskip, 0:indt + 4)
        real*8 distp_trueb(0:0, 20, nshell), distpb(nelskip, 0:indt + 4)
        real*8  cphs_r(2,0:0, *)
        real*8 rphs(0:0, nshell)
        real*8 cphs
        integer nshift,mshift
        logical yeszero_z
        rb = 0.d0
        rmub = 0.d0
        rmunewb = 0.d0
        rknewb = 0.d0
        distpb = 0.d0  ! Initialize to zero shared variable.
        if(iespbc)  s2rp = 0.d0

            if(iesjas) then
            nshift=nion
            mshift=nshell_det
            else
            nshift=0
            mshift=0
            endif


!$omp parallel do default(shared)  private(iii,jjj,i,j,ii,jj,kk,ll,indpar,indorb,indshell,kpip,cphs,do_makefun,yeszero_z)
        do jjj = 1, nion
            ! write(6,*) ' rmnew before mekefun -2=',jjj,sum(rmunewb(:,jjj,indtmin:indtm)),sum(rknewb(jjj,indtmin:indtm))
          do ii=1,kgrid_atom(jjj+nshift)%dimshell
            kpip(1) = kgrid_atom(jjj + nshift)%kpip(1, ii)
            kpip(2) = kgrid_atom(jjj + nshift)%kpip(2, ii)
            kpip(3) = kgrid_atom(jjj + nshift)%kpip(3, ii)


            rmunew(:, 0, jjj) = rmu(:, 0,jjj) + &
          &s2r(:, 1) * kpip(1) + s2r(:, 2) * kpip(2) + s2r(:, 3) * kpip(3)
    rknew(0,jjj)=max(dsqrt(rmunew(1,0,jjj)**2+rmunew(2,0,jjj)**2+rmunew(3,0,jjj)**2),mindist)

            select case(case_if)

            case(1)
!            if(ipc.eq.2.and..not.iesjas) then ! complex case
                        cphs = -(phs(1) * (kpip(1) - npip(0,1,jjj)) + &
                                phs(2) * (kpip(2) - npip(0,2,jjj)) + &
                                phs(3) * (kpip(3) - npip(0,3,jjj)))
                        cphs_r(1,0, jjj) = dcos(cphs)
                        cphs_r(2,0, jjj) = -dsin(cphs)
            case(2)
            rphs(0, jjj) = dcos(phs(1) * (kpip(1) - npip(0,1,jjj)) + &
                                phs(2) * (kpip(2) - npip(0,2,jjj)) + &
                                phs(3) * (kpip(3) - npip(0,3,jjj)))
            case(3)
              cphs_r(1,0, jjj) = 1.d0
              cphs_r(2,0, jjj) = 0.d0
            case(4)
              rphs(0,jjj)=1.d0
            end select


            DO iii = adr_nion(jjj), adr_nion(jjj + 1) - 1
                i = ind_nion(iii)
                ! write(6,*) ' rmnew before mekefun -1=',iii,sum(rmunewb(:,jjj,indtmin:indtm)),sum(rknewb(jjj,indtmin:indtm))
!               if(i.le.nshell) then
!                   do ii = 1, kgrid(i + ishift)%dimshell
                 if(kgrid(i+ishift)%tobedone(ii))  then
                        indpar = max(indpar_tab(i), 0)
                        indorb = indorb_tab(i + 1)
                        indshell = indshell_tab(i)
          do_makefun=.true.
          if(yes_scemama.and.ioptorb(i).ne.200) then
!           do_makefun=.false.
            if(slaterorb_read(i+mshift)) then
              if(dd(indpar+1)*rknew(0,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z=.false.
              else
              yeszero_z=.true.
              endif
            else
       if(dd(indpar+1)*rknew(0,jjj)*rknew(0,jjj).lt.lepsbas) then 
!             do_makefun=.true.
              yeszero_z=.false.
              else
              yeszero_z=.true.
              endif
            endif
            do_makefun=.not.yeszero_z
          endif
           if(do_makefun) then 
                    select case(case_upz) 
                    case(1) 
        call zaxpyr(indorb_tab(i),indorb,cphs_r(1,0,jjj),zb(1,0),distpb(1,0))
                    case(2)
        call zaxpyr(indorb_tab(i),indorb,cphs_r(1,0,jjj),zb(1,0),distpb(1,0))
                                do ll = indt + 1, indt + 4
        call zaxpyr(indorb_tab(i),indorb,cphs_r(1,0,jjj),zb(1,ll),distpb(1,ll))
                                enddo
                    case(3)
        call daxrpy(indorb_tab(i),indorb,rphs(0,jjj),zb(1,0),distpb(1,0))
                    case(4)
        call daxrpy(indorb_tab(i),indorb,rphs(0,jjj),zb(1,0),distpb(1,0))
                         do ll = indt + 1, indt + 4
        call daxrpy(indorb_tab(i),indorb,rphs(0,jjj),zb(1,ll),distpb(1,ll))
                         enddo
                    end select
                        CALL MAKEFUN0_B(ioptorb(i), indt,&
              & typecomp, indpar, indorb_tab(i), indshell_tab(i), nelskip&
              &, distp, distpb, dd, ddb, zeta, rknew(0,jjj), rknewb(0,jjj)&
              &, rmunew(1,0,jjj), rmunewb(1,0,jjj), distp_true(0, 1, i)&
              &, distp_trueb(0, 1, i), iflagnorm, cnorm(i))

                  endif ! endif do_makefun

                 endif ! endif to be done

                enddo  ! enddo over the shell of each atom

!  reverse of 
!            do j = indtmin, indtm
!            rmunew(:, j, i_ion) = rmu(:, j,i_ion) + &
!          &s2r(:, 1) * kpip(1) + s2r(:, 2) * kpip(2) + s2r(:, 3) * kpip(3)
!    rknew(j,i_ion)=max(dsqrt(rmunew(1,j,i_ion)**2+rmunew(2,j,i_ion)**2+rmunew(3,j,i_ion)**2),mindist)
!            enddo

! Below we do not  put the if rknew>mindist because it is negligible case that
! occurs  rarely and in case  the  correction is negligible instead of being zero.

                            rmunewb(1, 0,jjj) = rmunewb(1, 0,jjj) + rknewb(0,jjj) * &
                                    &rmunew(1, 0,jjj) / rknew(0,jjj)
                            rmunewb(2, 0,jjj) = rmunewb(2, 0,jjj) + rknewb(0,jjj) * &
                                    &rmunew(2, 0,jjj) / rknew(0,jjj)
                            rmunewb(3, 0,jjj) = rmunewb(3, 0,jjj) + rknewb(0,jjj) * &
                                    &rmunew(3, 0,jjj) / rknew(0,jjj)
                            rknewb(0,jjj) = 0.d0
                                do kk=1,3
                                s2rp(kk, 1, jjj) = s2rp(kk, 1, jjj) + rmunewb(kk, 0,jjj) * kpip(1)
                                s2rp(kk, 2, jjj) = s2rp(kk, 2, jjj) + rmunewb(kk, 0,jjj) * kpip(2)
                                s2rp(kk, 3, jjj) = s2rp(kk, 3, jjj) + rmunewb(kk, 0,jjj) * kpip(3)
                                enddo
                            rmub(1, 0,jjj) = rmub(1, 0,jjj) + rmunewb(1, 0,jjj)
                            rmub(2, 0,jjj) = rmub(2, 0,jjj) + rmunewb(2, 0,jjj)
                            rmub(3, 0,jjj) = rmub(3, 0,jjj) + rmunewb(3, 0,jjj)
                            rmunewb(1, 0,jjj) = 0.d0
                            rmunewb(2, 0,jjj) = 0.d0
                            rmunewb(3, 0,jjj) = 0.d0
               enddo ! enddo over the grid ii 
        enddo  ! enddo over jjj
!$omp end parallel do
        !  To avoid an inefficient reduction.
        if(iespbc)  then
            do j = 1, nion
                s2rb(:, :) = s2rb(:, :) + s2rp(:, :, j)
            enddo
        endif
    end subroutine makefun_grid_b
END SUBROUTINE UPNEWWF0_B
subroutine  zaxpyr(n,m,a,x,z)
implicit none
integer n,m,l,i
complex*16  x(*),a
real*8 z(*)
  l=m-n
  select case(l)
  case(1)
  z(n+1)=z(n+1)+a*x(n+1)
  case default
  do i=n+1,m
  z(i)=z(i)+a*x(i)
  enddo
 end select
return 
end
