!TL off
SUBROUTINE PSEUDOSET_B(jel, kel, kelb, ivic, ivicb, prefactor, &
        &  prefactorb, pseudolocal, pseudolocalb, nion, nintpseudo, dist, distb, rcutoff, &
        &  kindion, rion, rionb, pshell, nparpshell, parshell, lmax, wpseudo, &
        &  wpseudob, legendre, legendreb, versor, wintpseudo, jpseudo, ncore, &
        &  skip, indtm, iopt, angle, psip, psipb, lbox, mindist, iflagpseudo&
        &, cellscaleb, s2rb, car2cryb, metricb)
    USE ALLIO, only : norm_metric
    USE CELL, only : cellscale, at, s2r, Lmin, car2cry, metric, yes_tilted
    !  USE DIFFSIZES
    !  Hint: ISIZE1OFpsip should be the size of dimension 1 of array psip
    IMPLICIT NONE
    INTEGER :: jel, core, ncore, kindion(*), pshell(*), lmax, nparpshell(&
            &  lmax, *), lang, jj, kk, nion, nintpseudo, jpseudo(lmax, *), skip, indtm&
            &, npseudo, iflagpseudo, kimages, maximages, flagplus
    REAL*8 :: ivic(3, *), prefactor(skip, *), dist(nion, *), rcutoff(*), &
            &  distmin, parshell(3, *), wpseudo(*), legendre(lmax - 1, *), &
            &  versor(3, *), rion(3, *), wintpseudo(*), PSEUDOFUN, pseudolocal, scal&
            &, LEGFUN, rver(3), r(3), rsav(3), angle(3, 3), psip(*), lbox, mindist
    REAL*8 :: ivicb(3, *), prefactorb(skip, *), distb(nion, *), distminb, &
            &wpseudob(*), legendreb(lmax - 1, *), rionb(3, *), scalb, rb(3), rsavb(3), &
            &psipb(*), pseudolocalb, tempb, cellscaleb(3), s2rb(3, 3), car2cryb(3, 3), pseudolocali
    real*8  npip(3), rsav_before(3), r_before(3), rsav_beforeb(3), metricb(3, 3)
    !      Direct algorithm:
    !      input kel,rion,dist
    !      output prefactor (passed ) ,legendre,wpseudo, ivic (passed)
    !      pseudolocal (passed)

    LOGICAL :: iopt
    REAL*8 :: kel(3)
    REAL*8 :: kelb(3)

    if(iflagpseudo.ne.0.or.nintpseudo.le.0) return

    npseudo = 0
    do core = 1, ncore
        distmin = dist(kindion(core), jel)

        if(distmin.le.rcutoff(core)) then

            if(distmin.le.mindist) distmin = mindist

            do kk = 1, 3
                rsav(kk) = kel(kk) - rion(kk, kindion(core))
            enddo
            if(Lbox.gt.0.d0) then
                !   Find the image ion rion(:,kindion(core))-rdiff closest to kel
                !   NB I am not using PBC on electron positions
                !   here so it works also for APBC or whatever boundary condition.
                !                    call ApplyPBC(rsav,1)  ! the  output will be rsav-rion(:,..)+rdiff
                !                    if(yes_tilted) then
                rsav_before = rsav
                rsav(:) = car2cry(:, 1) * rsav(1) + car2cry(:, 2) * rsav(2) + car2cry(:, 3) * rsav(3)
                npip(1:3) = anint(rsav(1:3) / cellscale(1:3))
                rsav(1:3) = rsav(1:3) - npip(1:3) * cellscale(1:3)
                !                    else
                !                    npip(1:3)=anint(rsav(1:3)/cellscale(1:3))
                !                    rsav(1:3)=rsav(1:3)-npip(1:3)*cellscale(1:3)
                !                    endif
                if(rcutoff(core).gt.Lmin / 2.d0) then
                    maximages = 4
                else
                    maximages = 1
                endif
            else
                maximages = 1
            endif

            rsavb = 0.d0

            do kimages = 1, maximages

                pseudolocali = pseudolocalb
                !           redefine all for the images
                r = rsav
                if(kimages.eq.2) then
                    !              first image over x
                    if(r(1).gt.0.d0) then
                        r(1) = rsav(1) - cellscale(1)
                        flagplus = -1
                    else
                        r(1) = rsav(1) + cellscale(1)
                        flagplus = 1
                    endif
                    distmin = norm_metric(r, metric)
                elseif(kimages.eq.3) then
                    !              second image over y
                    if(r(2).gt.0.d0) then
                        r(2) = rsav(2) - cellscale(2)
                        flagplus = -1
                    else
                        r(2) = rsav(2) + cellscale(2)
                        flagplus = 1
                    endif
                    distmin = norm_metric(r, metric)
                elseif(kimages.eq.4) then
                    !              third image over z
                    if(r(3).gt.0.d0) then
                        r(3) = rsav(3) - cellscale(3)
                        flagplus = -1
                    else
                        r(3) = rsav(3) + cellscale(3)
                        flagplus = 1
                    endif
                    distmin = norm_metric(r, metric)
                endif

                if(LBox.gt.0.d0) then
                    !          go back in Cartesian Coordinates
                    r(1:3) = r(1:3) / cellscale(1:3)
                    r_before = r
                    r(:) = s2r(:, 1) * r(1) + s2r(:, 2) * r(2) + s2r(:, 3) * r(3)
                endif


                !       defined rb distmin
                distminb = 0.d0
                rb = 0.d0

                if(distmin.le.rcutoff(core)) then

                    !             indtm=indtm+nintpseudo  no update of indtm in reverse mode

                    npseudo = npseudo + 1

                    if(npseudo.gt.skip / nintpseudo) then
                        write(6, *) ' ERROR nuclei too close for pseudo !!! '
                        write(6, *) 'increase max # of nuclei', skip / nintpseudo
                        iflagpseudo = 1
                        return
                    endif

                    do lang = 1, pshell(core) - 1
                        wpseudo(lang) = pseudofun(nparpshell(lang, core)         &
                                &, distmin, parshell(1, jpseudo(lang, core)), psip)
                        !           Initializing wpseudob
                        wpseudob(lang) = 0.d0
                    enddo
                    !           pseudolocal=pseudolocal                                     &
                    !    &      +pseudofun(nparpshell(pshell(core),core)              &
                    !    &      ,distmin,parshell(1,jpseudo(pshell(core),core)),psip)

                    do jj = (npseudo - 1) * nintpseudo + 1, npseudo * nintpseudo

                        !              call rotation(angle
                        !    &         ,versor(1,jj-(npseudo-1)*nintpseudo),rver)
                        call dgemv('N', 3, 3, 1.d0, angle, 3                          &
                                &, versor(1, jj - (npseudo - 1) * nintpseudo), 1, 0.d0, rver, 1)

                        scal = 0.d0
                        scalb = 0.d0

                        do kk = 1, 3
                            scal = scal + rver(kk) * r(kk)
                            !                      ivic(kk,jj)=rver(kk)*distmin-r(kk)
                        enddo

                        scal = scal / distmin

                        do lang = 1, pshell(core) - 1
                            legendre(lang, jj - (npseudo - 1) * nintpseudo) = &
                                    &            legfun(lang - 1, scal)
                            legendreb(lang, jj - (npseudo - 1) * nintpseudo) = 0.d0
                        enddo

                        !                       do lang=1,pshell(core)-1
                        !                          prefactor(jj,jel)=prefactor(jj,jel)+                  &
                        !             &            (2.d0*lang-1.d0)*wpseudo(lang)                        &
                        !             &            *legendre(lang,jj-(npseudo-1)*nintpseudo)
                        !                       enddo
                        !                       prefactor(jj,jel)=prefactor(jj,jel)                      &
                        !             &         *2.d0*wintpseudo(jj-(npseudo-1)*nintpseudo)


                        !      REVERSE starts here

                        prefactorb(jj, jel) = prefactorb(jj, jel)&
                                & * 2.d0 * wintpseudo(jj - (npseudo - 1) * nintpseudo)

                        do lang = 1, pshell(core) - 1
                            tempb = (lang * 2 - 1) * prefactorb(jj, jel)
                            legendreb(lang, jj - (npseudo - 1) * nintpseudo) = &
                                    &         legendreb(lang, jj - (npseudo - 1) * nintpseudo) + &
                                    &         tempb * wpseudo(lang)
                            wpseudob(lang) = wpseudob(lang) + tempb&
                                    & * legendre(lang, jj - (npseudo - 1) * nintpseudo)
                        enddo

                        do lang = 1, pshell(core) - 1
                            call legfun_b(lang - 1, scal, scalb, legendreb(lang, jj - (npseudo - 1) * nintpseudo))
                            !                 legendre(lang,jj-(npseudo-1)*nintpseudo)=             &
                            !    &            legfun(lang-1,scal)
                        enddo
                        !   reverse of:           scal=scal/distmin
                        !           scal used below is the one before the rescaling.

                        distminb = distminb - scalb * scal / distmin
                        scalb = scalb / distmin

                        do kk = 1, 3
                            !    reverse of  scal=scal+rver(kk)*r(kk)
                            rb(kk) = rb(kk) + scalb * rver(kk)

                            !           reverse of ivic(kk,jj)=rver(kk)*distmin-r(kk)
                            distminb = distminb + ivicb(kk, jj) * rver(kk)
                            rb(kk) = rb(kk) - ivicb(kk, jj)
                            ivicb(kk, jj) = 0.d0
                        enddo

                    enddo

                    !            pseudolocal=pseudolocal                                     &
                    !     &      +pseudofun(nparpshell(pshell(core),core)              &
                    !     &      ,distmin,parshell(1,jpseudo(pshell(core),core)),psip)
                    !              mesh electronic point in the pseudo integration.
                    call pseudofun_b(nparpshell(pshell(core), core), distmin, distminb&
                            &, parshell(1, jpseudo(pshell(core), core)), psip, psipb, pseudolocali)

                    do lang = 1, pshell(core) - 1
                        call pseudofun_b(nparpshell(lang, core), distmin, distminb&
                                &, parshell(1, jpseudo(lang, core)), psip, psipb, wpseudob(lang))
                        !              wpseudo(lang)=pseudofun(nparpshell(lang,core)         &
                        !    &         ,distmin,parshell(1,jpseudo(lang,core)),psip)
                    enddo

                endif ! endif check distance images

                if(Lbox.gt.0) then
                    !  reverse of
                    !         r(:)=s2r(:,1)*r(1)+s2r(:,2)*r(2)+s2r(:,3)*r(3)

                    !  rb_sav=0.d0
                    !  call dgemv_b('N',3,3,1.d0,s2r,3,s2rb,3,r_before,1,rb_sav,1,0.d0,rb,1)
                    !  rb=rb_sav

                    s2rb(:, 1) = s2rb(:, 1) + rb(:) * r_before(1)
                    s2rb(:, 2) = s2rb(:, 2) + rb(:) * r_before(2)
                    s2rb(:, 3) = s2rb(:, 3) + rb(:) * r_before(3)

                    rb(:) = s2r(1, :) * rb(1) + s2r(2, :) * rb(2) + s2r(3, :) * rb(3)

                    !reverse of            r(1:3)=r(1:3)/cellscale(1:3)
                    cellscaleb(1:3) = cellscaleb(1:3) - r_before(1:3) / cellscale(1:3) * rb(1:3)
                    rb(1:3) = rb(1:3) / cellscale(1:3)

                endif
                if(kimages.ne.1) then

                    r_before(1:3) = r_before(1:3) * cellscale(1:3)
                    call norm_metric_b(r_before, rb, metric, metricb, distminb)


                    !  reverse of  r(kimages-1)=rsav(kimages-1)-cellscale(kimages-1)
                    cellscaleb(kimages - 1) = cellscaleb(kimages - 1) + rb(kimages - 1) * flagplus
                    rsavb(kimages - 1) = rsavb(kimages - 1) + rb(kimages - 1)
                    rb(kimages - 1) = 0.d0

                else
                    !  reverse of     distmin=dist(jel,kindion(core))
                    !   Note that  only  when kimages=1 the input dist is used otherwise distmin is
                    !    recomputed by scratch with norm_metric

                    if(distmin.gt.mindist) distb(kindion(core), jel) = distb(kindion(core), jel) + distminb
                    distminb = 0.d0

                endif

                !      reverse of  r=rsav
                rsavb = rsavb + rb

            enddo  ! enddo kimages

            if(Lbox.gt.0) then
                !            if(yes_tilted) then
                !         reverse of
                !         rsav(:)=rsav(kk)-cellscale(:)*anint(rsav(:)/cellscale(:))

                cellscaleb(:) = cellscaleb(:) - npip(:) * rsavb(:)


                !          reverse of

                !rsav(:)=car2cry(:,1)*rsav(1)+car2cry(:,2)*rsav(2)+car2cry(:,3)*rsav(3)
                !        rsav_beforeb=0.d0
                !         call dgemv_b('N',3,3,1.d0,car2cry,3,car2cryb,3,rsav_before,1,rsav_beforeb,1,0.d0,rsavb,1)
                !        rsavb=rsav_beforeb
                car2cryb(:, 1) = car2cryb(:, 1) + rsavb(:) * rsav_before(1)
                car2cryb(:, 2) = car2cryb(:, 2) + rsavb(:) * rsav_before(2)
                car2cryb(:, 3) = car2cryb(:, 3) + rsavb(:) * rsav_before(3)
                rsavb(:) = car2cry(1, :) * rsavb(1) + car2cry(2, :) * rsavb(2) + car2cry(3, :) * rsavb(3)

                !             call dger(3,3,-1.d0,rsavb,1,npip,1,s2rb,3)
                !            else
                !            cellscaleb(:)=cellscaleb(:)-npip(:)*rsavb(:)
                !            endif
            endif

            kelb(:) = kelb(:) + rsavb(:)
            rionb(:, kindion(core)) = rionb(:, kindion(core)) - rsavb(:)

        endif  ! endif check distance

    enddo

    pseudolocalb = 0.d0

END SUBROUTINE PSEUDOSET_B
