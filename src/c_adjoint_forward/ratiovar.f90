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

subroutine ratiovar(jel, nelorb, nelorbh, nelup, neldo, table &
        &, ainv, kel, rcart, iesdr, vj, dup, zeta, rion, psip &
        &, ioccup, ioccdo, ioptorb, nshell, nshelldo, ainvs, winvs, r &
        &, rmu, nion, kion, ioccj, kionj, vju, nelorbj, nelorbjh &
        &, ioptorbj, nshellj, winvsj, winvbar, winvjbar, winvjbarsz, iflagnorm &
        &, cnorm, iflagvar, ratiodetr, costz, costz3, iessz, LBox, rmucos, rmusin &
        &, timewf, rcne, jastrowall_ee, jasnew_ee, jastrowall_ei, jasnew_ei &
        &, n_body_on, niesd, nshelltot, epscuttype &
        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

    use Constants, only: ipc, ipj, zmone, zone, zzero, ipf
    use Cell, only: cellscale, metric, car2cry, CartesianToCrystal
    use allio, only: agp, agpn, jtype2, itest, pointvj, nelup_mat, epst, epsvar, yes_crystalj, iespbc, dim_ratiovar, yes_ontarget, norm_metric

    implicit none

    integer nelup, neldo, nel, i, jel, niesd, iflagvar                    &
            &, j, iesd, iesdr, nelorb, ioccup(*), ioccdo(*), ioptorb(*), nshell        &
            &, nshelldo, ispin, ispinj, nion, kion(*)         &
            &, ioccj(*), kionj(*), nelorbj, ioptorbj(*), nshellj, iflagnorm          &
            &, nelorbh, nelorbjh, nelcdo, iesdr2iesd            &
            &, nshelltot, n_body_on, indt4, epscuttype, indexi

    integer :: indpar_tab(*), indorb_tab(*), indshell_tab(*), indparj_tab(*)&
            &, indorbj_tab(*), indshellj_tab(*)

    real(8) rcart(3), rcartpip(3, 0:1), kel(3, nelup + neldo)

    real(8) rion(3, nion), jastrow_ee, jastrow_ei, psip(dim_ratiovar)   &
            &, rold, vj(max(1, niesd)), dup(*), zeta(*), r(nion, 0:1), rmu(3, nion, 1), r0       &
            &, winvjbar(max(ipj*nelorbjh, 1), *), winvjbarsz(*), rc(3), rcn(3) &
            &, winvsj(max(nelorbjh, 1)), vju(*), cost, cnorm(*), weight&
            &, psiav, costz(nion), costzj       &
            &, costz3(nion), costz0, signszup, signszdo, jastrowall_ee(nelup + neldo)     &
            &, jastrowall_ei(nion), jasnew_ee(nelup + neldo), jasnew_ei(nion), diffj3body&
            &, psidetln, psidetsn, winv(2)

    ! real/complex quantities

    real(8) ratiodetr, ratiodet(ipc), g(ipc), ainv(ipc*nelup_mat, *), table(ipc)&
            &, ainvs(ipc*nelup_mat, 2 + min(epscuttype, 1)*(nelup_mat - 2))&
            &, winvs(ipc*nelorbh), winvbar(ipc*ipf*nelorbh, *), g1

    real(8), external :: ddot
    complex(8), external :: zdotu, zdotc
    logical iessz, iesspin, reject
    real(8) cclock
    real(8) LBox, rmucos(*), rmusin(*), rcne(3, nelup + neldo), timewf, timep

    reject = .false.

    if (iflagvar .ne. 0) return

    iesd = iesdr2iesd(iesdr)

    !     iesd=iesdr
    !     if(iesdr.eq.-5) iesd=1
    !     if(iesdr.eq.-6) iesd=4
    !     if(iesdr.eq.-7) iesd=4

    if (iesdr .eq. -7 .or. iesd .eq. 2 .or. iesd .lt. 0) then
        iesspin = .true.
    else
        iesspin = .false.
    end if

    nel = nelup + neldo
    table(1) = 1.d0
    if (ipc .eq. 2) table(2) = 0.d0

    if (jel .le. nelup) then ! move an electron with spin up

        timep = cclock()

        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, rcart, 1, r, rmu          &
                &, dup, zeta, rion, psip, winvs, nelorbh, nion, kion, iflagnorm, cnorm &
                &, LBox, rmucos, rmusin, 1d-9, indpar_tab, indorb_tab, indshell_tab, .true.)
#ifdef _OFFLOAD
!     if(yes_ontarget) then
!$omp target  update to (winvs) if(yes_ontarget)
!     endif
#endif

        if (ipc .eq. 1) then ! real case

            timewf = timewf + cclock() - timep
            if (ipf .eq. 2) then
                call dgemv_('T', nelorbh, nelup_mat, 1.d0, winvbar, 2*nelorbh, winvs, 1, 0.d0, ainvs, 1, yes_ontarget)
                ainvs(jel, 1) = 0.d0
#ifdef _OFFLOAD
!$omp target update to(ainvs(jel:jel,1:1)) if(yes_ontarget)
#endif

                call dgemv_('T', nelup_mat, nelup_mat, 1.d0, ainv, nelup_mat, ainvs, 1, 0.d0, ainvs(1, 2), 1, yes_ontarget)

            else

                call dgemv_('T', nelorbh, nelup, 1.d0, winvbar(1, nelup + 1), nelorbh, winvs, 1, 0.d0, ainvs, 1, yes_ontarget)
                call dgemv_('T', nelup, nelup, 1.d0, ainv, nelup, ainvs, 1, 0.d0, ainvs(1, 2), 1, yes_ontarget)

            end if

#ifdef _OFFLOAD
!$omp target update from(ainvs(jel:jel,2:2)) if(yes_ontarget)
#endif
            ratiodet(1) = ainvs(jel, 2)
            ainvs(jel, 2) = ainvs(jel, 2) - 1.d0
#ifdef _OFFLOAD
!$omp target update to (ainvs(jel:jel,2:2)) if(yes_ontarget)
#endif

            table(1) = table(1)*ratiodet(1)

            if (epscuttype .eq. 2) then
                !          agpn(:,:)=agp(:,:,jtype2)
                call dcopy_vec(nelup_mat*nelup_mat, agp(1, 1, jtype2), agpn)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup_mat
                    agpn(jel, i) = ainvs(i, 1)
                    if (ipf .eq. 2) agpn(i, jel) = -ainvs(i, 1)
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
            end if
            if (epscuttype .gt. 0) then
                !         redefinition of ratiodet
                if (epsvar .ne. 0.d0 .and. (abs(ratiodet(1)) .lt. epsvar .or. abs(ratiodet(1)) .gt. 1.d0/epsvar)) then
                    ratiodet(1) = 0.d0
                    reject = .true.
                else
                    if (ipf .eq. 1) then
                        g1 = -1.d0/ratiodet(1)
                        call dcopy_vec(nelup, ainvs(1, 2), psip)
                        call dcopy_vec(nelup*nelup, ainv, ainvs)
                        call dger_(nelup, nelup, g1, ainv(1, jel), 1, psip, 1, ainvs, nelup)
                    else
                        call dcopy_vec(2*nelup_mat, ainvs, psip)
                        call dcopy_vec(nelup_mat*nelup_mat, ainv, ainvs)
                        call upinvhop(1, 1, jel, ainvs, winv, psidetln, &
                                      psidetsn, epst, psip(2*nelup_mat + 1), nelup, neldo, psip, winvs, .false.)
                    end if
                end if
                !          call dcopy(nelup*nelup,ainv,1,psip,1)
                !          call dger(nelup,nelup,g,ainv(1,jel),1,ainvs(1,2),1,psip,nelup)
            end if
            ratiodetr = ratiodet(1)
        else ! complex case

            timewf = timewf + cclock() - timep

            if (ipf .eq. 2) then
                !SSS controllare se vengono chiamati con la giusta dimensione
                call zgemv_('T', nelorbh, nelup_mat, zone, winvbar, 2*nelorbh, winvs, 1, zzero, ainvs, 1, yes_ontarget)

                ainvs(2*jel - 1:2*jel, 1) = 0.d0
#ifdef _OFFLOAD
!$omp target update to(ainvs(2*jel-1:2*jel,1)) if(yes_ontarget)
#endif
                call zgemv_('T', nelup_mat, nelup_mat, zone, ainv, nelup_mat, ainvs, 1, zzero, ainvs(1, 2), 1, yes_ontarget)
            else

                call zgemv_('T', nelorbh, nelup, zone, winvbar(1, nelup + 1), nelorbh, winvs, 1, zzero, ainvs, 1, yes_ontarget)
                call zgemv_('T', nelup, nelup, zone, ainv, nelup, ainvs, 1, zzero, ainvs(1, 2), 1, yes_ontarget)

            end if

            if (epscuttype .eq. 2) then
                !          agpn(:,:)=agp(:,:,jtype2)
                call dcopy_vec(2*nelup_mat*nelup_mat, agp(1, 1, jtype2), agpn)
                !          call zcopy(nelup_mat,ainvs,1,agpn(2*jel-1,1),nelup_mat)
                !          if(ipf.eq.2) agpn(:,jel)=-ainvs(1:2*nelup_mat,1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup_mat
                    agpn(2*jel - 1, i) = ainvs(2*i - 1, 1)
                    agpn(2*jel, i) = ainvs(2*i, 1)
                    if (ipf .eq. 2) then
                        agpn(2*i - 1, jel) = -ainvs(2*i - 1, 1)
                        agpn(2*i, jel) = -ainvs(2*i, 1)
                    end if
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
            end if

            call update_ratio_complex(jel, ratiodet, ainvs, ainv, nelup_mat, table, g, psip, reject)

            ratiodetr = sqrt(sum(ratiodet(1:ipc)**2))
            if (ratiodet(1) .lt. 0) ratiodetr = -ratiodetr

        end if ! endif ipc.eq.1

    else ! move an electron with spin down

        nelcdo = jel - nelup
        timep = cclock()

        call upnewwf(0, 0, 0, 1, nshelldo, ioptorb, ioccdo, rcart, 1, r, rmu         &
                &, dup, zeta, rion, psip, winvs, nelorbh, nion, kion, iflagnorm, cnorm  &
                &, LBox, rmucos, rmusin, 1d-9, indpar_tab, indorb_tab, indshell_tab, .false.)
#ifdef _OFFLOAD
!     if(yes_ontarget) then
!$omp target  update to (winvs) if(yes_ontarget)
!     endif
#endif

        if (ipc .eq. 1) then ! real case

            timewf = timewf + cclock() - timep

            !    down Non ho ben capito quale debba essere la differenza con il caso jel> nelup,
            !unica modifica fatta è nell'ultima dgemv 'T'-> 'N' per analogia

            ! controllare se è giusta la versione 'T' 'T' o la 'T' 'N' (secondome sempre la 'T' 'T')
            if (ipf .eq. 2) then
              call dgemv_('T', nelorbh, nelup_mat, 1.d0&
                  &, winvbar(1 + nelorbh, 1), 2*nelorbh, winvs, 1, 0.d0, ainvs&
                  &, 1, yes_ontarget)

                ainvs(jel, 1) = 0.d0
#ifdef _OFFLOAD
!$omp target update to(ainvs(jel:jel,1:1)) if(yes_ontarget)
#endif
                call dgemv_('T', nelup_mat, nelup_mat, 1.d0, ainv, nelup_mat, ainvs, 1, 0.d0, ainvs(1, 2), 1, yes_ontarget)
                !         call dgemv('T',nelorbh,2*nelup,1.d0,winvbar(nelorbh+1,1),2*nelorbh,winvs,1,0.d0,ainvs,1)
                !         call dgemv('N',2*nelup,2*nelup,1.d0,ainv,2*nelup,ainvs,1,0.d0,ainvs(1,2),1)
            else

                ! call dgemv('T',nelorbh,nelup,1.d0,winvbar(1,nelup+1),nelorbh,winvs,1,0.d0,ainvs,1)
                ! call dgemv('T',nelup,nelup,1.d0,ainv,nelup,ainvs,1,0.d0,ainvs(1,2),1)

                call dgemv_('T', nelorbh, nelup, 1.d0, winvbar, nelorbh, winvs, 1, 0.d0, ainvs, 1, yes_ontarget)
                call dgemv_('N', nelup, nelup, 1.d0, ainv, nelup, ainvs, 1, 0.d0, ainvs(1, 2), 1, yes_ontarget)
            end if

            if (epscuttype .eq. 2) then
                call dcopy_vec(nelup_mat*nelup_mat, agp(1, 1, jtype2), agpn)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup_mat
                    if (ipf .eq. 2) then
                        agpn(jel, i) = ainvs(i, 1)
                        agpn(i, jel) = -ainvs(i, 1)
                    else
                        agpn(i, nelcdo) = ainvs(i, 1)
                    end if
                end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif
            end if

            if (ipf .eq. 1) then
#ifdef _OFFLOAD
!$omp target update from(ainvs(nelcdo:nelcdo,2:2))  if(yes_ontarget)
#endif
                ratiodet(1) = ainvs(nelcdo, 2)
                ainvs(nelcdo, 2) = ainvs(nelcdo, 2) - 1.d0
                g1 = -1.d0/ratiodet(1)
#ifdef _OFFLOAD
!$omp target update to(ainvs(nelcdo:nelcdo,2:2))  if(yes_ontarget)
#endif
            else
#ifdef _OFFLOAD
!$omp target update from(ainvs(jel:jel,2:2))  if(yes_ontarget)
#endif
                ratiodet(1) = ainvs(jel, 2)
                ainvs(jel, 2) = ainvs(jel, 2) - 1.d0
#ifdef _OFFLOAD
!$omp target update to(ainvs(jel:jel,2:2))  if(yes_ontarget)
#endif
            end if

            table(1) = table(1)*ratiodet(1)

            if (epscuttype .gt. 0) then
                !         redefinition of ratiodet
                if (epsvar .ne. 0.d0 .and. (abs(ratiodet(1)) .lt. epsvar .or. abs(ratiodet(1)) .gt. 1.d0/epsvar)) then
                    ratiodet(1) = 0.d0
                    reject = .true.
                else
                    if (ipf .eq. 1) then
                        !          call dcopy(nelup*nelup,ainv,1,psip,1)
                        call dcopy_vec(nelup, ainvs(1, 2), psip)
                        call dcopy_vec(nelup*nelup, ainv, ainvs)
                        !          call dger(nelup,nelup,g,ainvs(1,2),1,ainv(nelcdo,1),nelup,ainvs,nelup)
                        call dcopy_my(nelup, ainv(nelcdo, 1), nelup, psip(nelup + 1), 1)
                        call dger_(nelup, nelup, g1, psip, 1, psip(nelup + 1), 1, ainvs, nelup)
                    else
                        call dcopy_vec(2*nelup_mat, ainvs, psip)
                        call dcopy_vec(nelup_mat*nelup_mat, ainv, ainvs)
                        call upinvhop(1, 1, jel, ainvs, winv, psidetln, &
                                      psidetsn, epst, psip(2*nelup_mat + 1), nelup, neldo, psip, winvs, .false.)
                    end if
                end if
            end if

            ratiodetr = ratiodet(1)

        else ! complex case

            timewf = timewf + cclock() - timep

            if (ipf .eq. 2) then
           call zgemv_('T', nelorbh, nelup_mat, zone, winvbar(1 + 2*nelorbh, 1)&
               &, 2*nelorbh, winvs, 1, zzero, ainvs, 1, yes_ontarget)
                ainvs(2*jel - 1:2*jel, 1) = 0.d0
#ifdef _OFFLOAD
!$omp target update from (ainvs(2*jel-1:2*jel,1:1)) if(yes_ontarget)
#endif
                call zgemv_('T', nelup_mat, nelup_mat, zone, ainv, nelup_mat, ainvs, 1, zzero, ainvs(1, 2), 1, yes_ontarget)
                !           call zgemv('T',nelorbh,2*nelup,zone,winvbar(1+2*nelorbh,1),2*nelorbh,winvs,1,zzero,ainvs,1)
                !           call zgemv('N',2*nelup,2*nelup,zone,ainv,2*nelup,ainvs,1,zzero,ainvs(1,2),1)
            else
                call zgemv_('T', nelorbh, nelup, zone, winvbar, nelorbh, winvs, 1, zzero, ainvs, 1, yes_ontarget)
                call zgemv_('N', nelup, nelup, zone, ainv, nelup, ainvs, 1, zzero, ainvs(1, 2), 1, yes_ontarget)
            end if
            if (epscuttype .eq. 2) then
                call dcopy_vec(2*nelup_mat*nelup_mat, agp(1, 1, jtype2), agpn)
                !          call zcopy(nelup_mat,ainvs,1,agpn(2*jel-1,1),nelup_mat)
                !          if(ipf.eq.2) agpn(:,jel)=-ainvs(1:2*nelup_mat,1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup_mat
                    if (ipf .eq. 2) then
                        agpn(2*jel - 1, i) = ainvs(2*i - 1, 1)
                        agpn(2*jel, i) = ainvs(2*i, 1)
                        agpn(2*i - 1, jel) = -ainvs(2*i - 1, 1)
                        agpn(2*i, jel) = -ainvs(2*i, 1)
                    else
                        agpn(2*i - 1, nelcdo) = ainvs(2*i - 1, 1)
                        agpn(2*i, nelcdo) = ainvs(2*i, 1)
                    end if
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
            end if
            call update_ratio_complex(jel, ratiodet, ainvs, ainv, nelup_mat, table, g, psip, reject)
            ratiodetr = sqrt(sum(ratiodet(1:ipc)**2))
            if (ratiodet(1) .lt. 0) ratiodetr = -ratiodetr
        end if

    end if !end of deteminant change

    if (epscuttype .gt. 0 .and. .not. reject) &
        call psireg(ainvs, agpn, psip, nelup_mat, neldo, ratiodetr, yes_ontarget)

    costz0 = 0.d0

    if (nshellj .ne. 0) then
        timep = cclock()
        call upnewwf(0, 0, 0, 1, nshellj, ioptorbj, ioccj, rcart, 1, r, rmu      &
                &, vju, zeta, rion, psip, winvsj, nelorbjh, nion, kionj, iflagnorm          &
                &, cnorm(nshelltot + 1), LBox, rmucos, rmusin, 1d-9&
                &, indparj_tab, indorbj_tab, indshellj_tab, .true.)
        timewf = timewf + cclock() - timep
#ifdef _OFFLOAD
!     if(yes_ontarget) then
!$omp target  update to (winvsj) if(yes_ontarget)
!     endif
#endif
        if (jel .le. nelup .or. ipj .eq. 1) then
            call dgemv_('T', nelorbjh, nel, 1.d0, winvjbar, ipj*nelorbjh, winvsj, 1, 0.d0, psip, 1, yes_ontarget)
        else
            call dgemv_('T', nelorbjh, nel, 1.d0, winvjbar(1 + nelorbjh, 1), 2*nelorbjh, winvsj, 1, 0.d0, psip, 1, yes_ontarget)
        end if
#ifdef _OFFLOAD
!     if(yes_ontarget) then
!$omp target  update from (psip(1:nel)) if(yes_ontarget)
!     endif
#endif

        if (iessz) then

            if (jel .le. nelup) then
                signszup = 1.d0
                signszdo = -1.d0
            else
                signszup = -1.d0
                signszdo = 1.d0
            end if

            call dgemv('T', nelorbjh, nelup, signszup, winvjbarsz, nelorbjh, winvsj, 1, 1.d0, psip, 1)
            if(neldo.gt.0) call dgemv('T', nelorbjh, neldo, signszdo, winvjbarsz(1 + nelup * nelorbjh)&
                &, nelorbjh, winvsj, 1, 1.d0, psip(nelup + 1), 1)

        end if

    end if

    if (iesspin) then
        if (jel .le. nelup) then
            ispinj = 1
        else
            ispinj = -1
        end if
    end if

    do i = 1, nel
        rcne(1, i) = rcart(1) - kel(1, i)
        rcne(2, i) = rcart(2) - kel(2, i)
        rcne(3, i) = rcart(3) - kel(3, i)
    end do
    if (LBox .gt. 0.d0) then
        call CartesianToCrystal(rcne, nel)
        call scalevect(nel, cellscale, rcne)
    end if

    if (nshellj .gt. 0) then
!$omp parallel do default(shared) reduction(+:costz0) &
!$omp private(i,ispin)
        do i = 1, nel
            if (iesspin) then
                if (i .le. nelup) then
                    ispin = ispinj
                else
                    ispin = -ispinj
                end if
            else
                ispin = -1
            end if
            if (i .ne. jel) then
                jasnew_ee(i) = jastrow_ee(rcne(1, i), vj, iesd, ispin) + psip(i)
                costz0 = costz0 + jasnew_ee(i) - jastrowall_ee(i)
            else
                jasnew_ee(i) = 0.d0
            end if
        end do
!$omp end parallel do
    else

!$omp parallel do default(shared) reduction(+:costz0) &
!$omp private(i,ispin)
        do i = 1, nel
            if (iesspin) then
                if (i .le. nelup) then
                    ispin = ispinj
                else
                    ispin = -ispinj
                end if
            else
                ispin = -1
            end if
            if (i .ne. jel) then
                jasnew_ee(i) = jastrow_ee(rcne(1, i), vj, iesd, ispin)
                costz0 = costz0 + jasnew_ee(i) - jastrowall_ee(i)
            else
                jasnew_ee(i) = 0.d0
            end if
        end do
!$omp end parallel do
    end if

    if (n_body_on .ne. 0) then
        if (yes_crystalj .or. nshellj .eq. 0) then
            !     recomputation rmu
            do i = 1, nion
                rmu(1:3, i, 1) = rcart(1:3) - rion(1:3, i)
            end do
            if (iespbc) then
                call CartesianToCrystal(rmu, nion)
                call scalevect(nion, cellscale, rmu)
            end if
        end if
        if (itest .eq. 2) then
            costzj = 0.d0 !  Vanish this in the device
!$omp  parallel do default(shared) reduction(+:costzj) &
!$omp  private(i,r0)
            do i = 1, nion
                if (iespbc) then
                    r0 = norm_metric(rmu(1, i, 1), metric)*costz(i)
                else
                    r0 = dsqrt(sum(rmu(:, i, 1)**2))*costz(i)
                end if
                costzj = costzj + costz3(i)*jastrow_ei(r0, vj(pointvj(1, i)), pointvj(2, i))
            end do
!$omp end parallel do
            costz0 = costz0 - costzj + jastrowall_ei(1)
            jasnew_ei(1) = costzj
        else
!$omp  parallel do default(shared) reduction(+:costz0) &
!$omp  private(i,r0)
            do i = 1, nion
                if (iespbc) then
                    r0 = norm_metric(rmu(1, i, 1), metric)*costz(i)
                else
                    r0 = dsqrt(sum(rmu(:, i, 1)**2))*costz(i)
                end if
                jasnew_ei(i) = jastrow_ei(r0, vj(pointvj(1, i)), pointvj(2, i))
                costz0 = costz0 + &
                        &costz3(i)*(jastrowall_ei(i) - jasnew_ei(i))
            end do
!$omp end parallel do
        end if
    else
        jasnew_ei(1) = 0.d0
    end if
    table(1:ipc) = table(1:ipc)*dexp(costz0)

    return

end subroutine ratiovar
!
real(8) function dnrmsq(n, dx, incx)
    implicit none
    integer n, i, incx
    real(8) dx(n*incx)
    dnrmsq = 0.d0
!$omp parallel do default(shared) private(i) reduction(+:dnrmsq)
    do i = 1, n*incx, incx
        dnrmsq = dnrmsq + dx(i)*dx(i)
    end do
    return
end function dnrmsq
!
real(8) function dnrmsq_(n, dx, incx)
    implicit none
    integer n, i, incx
    real(8) dx(n*incx)
    dnrmsq_ = 0.d0
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do reduction(+:dnrmsq_)
#endif
    do i = 1, n*incx, incx
        dnrmsq_ = dnrmsq_ + dx(i)*dx(i)
    end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif
    return
end function dnrmsq_
!
!
subroutine measure_dis(kel, nel, rcart, jel)
    implicit none
    real*8 kel(3, *), dismax, rcart(3), dnrmsq
    integer jel, nel, ii, theone
    dismax = dnrmsq(3, rcart, 1)
    theone = jel
    write (6, *) 'jel-th electron at dis ', dsqrt(dismax)
    do ii = 1, nel
        if (ii .ne. jel .and. dnrmsq(3, kel(1, ii), 1) > dismax) then
            dismax = dnrmsq(3, kel(1, ii), 1)
            theone = ii
        end if
    end do
    write (6, *) 'farthest electron ', theone, ' at dis ', dsqrt(dismax)
end subroutine measure_dis

! this subroutine implementes the regularization scheme for
! OPEN SYSTEMS based on a redefinition of A^-1
! Best choice for open systems: epscuttype=3
! Best choice for PBC systems: epscuttype=1
! followed by complex version
subroutine psireg(ainv, agpn, psip, nelup, neldo, ratiodet, yes_ontarget)
    use allio, only: epscuttype, theta_reg
    use constants, only: ipc, ipf, safemin
    !        USE ieee_arithmetic
    implicit none
    logical yes_ontarget
    integer nelup, neldo, i, j
    real(8) ainv(ipc*nelup, nelup), agpn(ipc*nelup, nelup), ratiodet, scaledet, cost
    real(8), external :: dnrmsq, dnrmsq_, dnrm2
    real(8) psip(1:4*nelup)
    !       variables for epscuttup=2
    real(8) csum
    !       variables for epscuttup=3

    if (epscuttype .lt. 0 .or. epscuttype .gt. 3) return

    if (epscuttype .le. 2 .and. epscuttype .ge. 0) then
        if (yes_ontarget) then
            ratiodet = dnrmsq_(ipc*nelup*nelup, ainv, 1) ! simple regularization for PBC systems
        else
            ratiodet = dnrmsq(ipc*nelup*nelup, ainv, 1)
        end if
    end if

    if (yes_ontarget) then

        if (epscuttype .eq. 2) then
            scaledet = 0.d0
#ifdef _OFFLOAD
!$omp target teams distribute parallel do
#endif
            do i = 1, 2*nelup
                psip(i) = 0.d0
            end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
            if (ipc .eq. 1) then
#ifdef _OFFLOAD
!$omp target teams distribute private(csum)
#else
!$omp parallel do  default(shared) private(i,j,csum)
#endif
                do i = 1, nelup
                    csum = psip(i)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do j = 1, nelup
                        csum = csum + agpn(j, i)**2
                    end do
                    psip(i) = csum
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute
#else
!$omp end parallel do
#endif
                if (ipf .eq. 1) then
#ifdef _OFFLOAD
!$omp target teams distribute  private(csum)
#else
!$omp parallel do default(shared) private(i,j,csum)
#endif
                    do i = 1, nelup
                        csum = psip(i + nelup)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                        do j = 1, nelup
                            csum = csum + agpn(i, j)**2
                        end do
                        psip(i + nelup) = csum
                    end do
#ifdef _OFFLOAD
!$omp end target teams distribute
#else
!$omp end parallel do
#endif
                end if
            else
#ifdef _OFFLOAD
!$omp target teams distribute private(csum)
#else
!$omp parallel do default(shared) private(i,j,csum)
#endif
                do i = 1, nelup
                    csum = psip(i)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do j = 1, 2*nelup
                        csum = csum + agpn(j, i)**2
                    end do
                    psip(i) = csum
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute
#else
!$omp end parallel do
#endif

                if (ipf .eq. 1) then
#ifdef _OFFLOAD
!$omp target teams distribute private(csum)
#else
!$omp parallel do default(shared) private(i,j,csum)
#endif
                    do i = 1, nelup
                        csum = psip(i + nelup)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                        do j = 1, nelup
                            csum = csum + agpn(2*i - 1, j)**2 + agpn(2*i, j)**2
                        end do
                        psip(i + nelup) = csum
                    end do

#ifdef _OFFLOAD
!$omp end target teams distribute
#else
!$omp end parallel do
#endif

                end if
            end if
#ifdef _OFFLOAD
!$omp target teams distribute reduction(max:scaledet) private(cost)
#endif
            do i = 1, nelup
                if (psip(i) .gt. safemin) then
                    cost = 1.d0/psip(i)
                else
                    cost = 1.d0/safemin
                end if
                scaledet = max(scaledet, cost)
                if (ipf .eq. 1) then
                    if (psip(i + nelup) .gt. safemin) then
                        cost = 1.d0/psip(i + nelup)
                    else
                        cost = 1.d0/safemin
                    end if
                    scaledet = max(scaledet, cost)
                end if
            end do
#ifdef _OFFLOAD
!$omp end target teams distribute
#endif

            ratiodet = ratiodet/scaledet
        end if
    else
!   repeated the same with  no OFFLOAD  instructions
        if (epscuttype .eq. 2) then
            scaledet = 0.d0
            do i = 1, 2*nelup
                psip(i) = 0.d0
            end do
            if (ipc .eq. 1) then
!$omp parallel do  default(shared) private(i,j,csum)
                do i = 1, nelup
                    csum = psip(i)
                    do j = 1, nelup
                        csum = csum + agpn(j, i)**2
                    end do
                    psip(i) = csum
                end do
                if (ipf .eq. 1) then
!$omp parallel do default(shared) private(i,j)
                    do i = 1, nelup
                        do j = 1, nelup
                            psip(i + nelup) = psip(i + nelup) + agpn(i, j)**2
                        end do
                    end do
                end if
            else
!$omp parallel do default(shared) private(i,j,csum)
                do i = 1, nelup
                    csum = psip(i)
                    do j = 1, 2*nelup
                        csum = csum + agpn(j, i)**2
                    end do
                    psip(i) = csum
                end do
                if (ipf .eq. 1) then
!$omp parallel do default(shared) private(i,j)
                    do i = 1, nelup
                        do j = 1, nelup
                            psip(i + nelup) = psip(i + nelup) + agpn(2*i - 1, j)**2 + agpn(2*i, j)**2
                        end do
                    end do
                end if
            end if
            do i = 1, nelup
                if (psip(i) .gt. safemin) then
                    cost = 1.d0/psip(i)
                else
                    cost = 1.d0/safemin
                end if
                scaledet = max(scaledet, cost)
                if (ipf .eq. 1) then
                    if (psip(i + nelup) .gt. safemin) then
                        cost = 1.d0/psip(i + nelup)
                    else
                        cost = 1.d0/safemin
                    end if
                    scaledet = max(scaledet, cost)
                end if
            end do
            ratiodet = ratiodet/scaledet
        end if
    end if

    ! final ratio regularization
    ratiodet = 1.d0/ratiodet**theta_reg
    return
end subroutine psireg

subroutine update_ratio_complex(jel, ratiodet, ainvs, ainv, nelup, table, g, psip, reject)

    use allio, only: epscuttype, nelup_mat, epst, epsvar, yes_ontarget
    use constants, only: zmone, zone, ipf
    implicit none

    integer, intent(in) :: jel, nelup
    complex(8) g, ratiodet, table, ainvs(nelup_mat, 2 + min(epscuttype, 1)&
       & * (nelup_mat - 2)), psip(6 * nelup_mat), ainv(nelup_mat, nelup_mat)
    integer indel, i
    real*8 psidetln, psidetsn, winv(2), winvs(2)
    logical reject
    if (jel .le. nelup .or. ipf .eq. 2) then
#ifdef _OFFLOAD
!$omp target update from(ainvs(jel:jel,2:2)) if(yes_ontarget)
#endif
        ratiodet = ainvs(jel, 2)
        ainvs(jel, 2) = ainvs(jel, 2) - zone
#ifdef _OFFLOAD
!$omp target update to(ainvs(jel:jel,2:2)) if(yes_ontarget)
#endif
        table = table*ratiodet
        if (epscuttype .gt. 0) then
            !         redefinition of ratiodet
            if (epsvar .ne. 0.d0 .and. (abs(ratiodet) .lt. epsvar .or. abs(ratiodet) .gt. 1.d0/epsvar)) then
                ratiodet = (0.d0, 0.d0)
                reject = .true.
            else
                if (ipf .eq. 2) then
                    !    call zcopy(2*nelup_mat,ainvs,1,psip,1)
                    call dcopy_vec(4*nelup_mat, ainvs, psip)
                    !    call zcopy(nelup_mat*nelup_mat,ainv,1,ainvs,1)
                    call dcopy_vec(2*nelup_mat*nelup_mat, ainv, ainvs)
                    !   neldo, winv, winvs are not used with .false.
                    call upinvhop(1, 1, jel, ainvs, winv, psidetln, &
                                  psidetsn, epst, psip(2*nelup_mat + 1), nelup, nelup, psip, winvs, .false.)
                else
                    g = zmone/ratiodet
                    !    call zcopy(nelup,ainvs(1,2),1,psip,1)
                    call dcopy_vec(2*nelup, ainvs(1, 2), psip)
                    !    call zcopy(nelup*nelup,ainv,1,ainvs,1)
                    call dcopy_vec(2*nelup*nelup, ainv, ainvs)
                    !   call zcopy(nelup*nelup,ainv,1,psip,1)
                    !   call zgeru(nelup,nelup,g,ainv(1,jel),1,ainvs(1,2),1,psip,nelup)
                    call zgeru_(nelup, nelup, g, ainv(1, jel), 1, psip, 1, ainvs, nelup)
                end if
            end if
        end if
    else
        indel = jel - nelup
#ifdef _OFFLOAD
!$omp target update from(ainvs(indel:indel,2:2)) if(yes_ontarget)
#endif
        ratiodet = ainvs(indel, 2)
        ainvs(indel, 2) = ainvs(indel, 2) - zone
#ifdef _OFFLOAD
!$omp target update to(ainvs(indel:indel,2:2)) if(yes_ontarget)
#endif
        table = table*ratiodet
        if (epscuttype .gt. 0) then
            ! redefinition of ratiodet
            if (epsvar .ne. 0.d0 .and. (abs(ratiodet) .lt. epsvar .or. abs(ratiodet) .gt. 1.d0/epsvar)) then
                ratiodet = (0.d0, 0.d0)
                reject = .true.
            else
                g = zmone/ratiodet
                !     call zcopy(nelup,ainvs(1,2),1,psip,1)
                call dcopy_vec(2*nelup, ainvs(1, 2), psip)
                !     call zcopy(nelup*nelup,ainv,1,ainvs,1)
                call dcopy_vec(2*nelup*nelup, ainv, ainvs)
                !     call zcopy_my(nelup,ainv(indel,1),nelup,psip(nelup+1),1)
                call zcopy_my(nelup, ainv(indel, 1), nelup, psip(nelup + 1), 1)
                call zgeru_(nelup, nelup, g, psip, 1, psip(nelup + 1), 1, ainvs, nelup)
                !     call zcopy(nelup*nelup,ainv,1,psip,1)
                !     call zgeru(nelup,nelup,g,ainvs(1,2),1,ainv(indel,1),nelup,psip,nelup)
            end if
        end if
    end if
    return
end subroutine update_ratio_complex
