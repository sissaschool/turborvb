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

subroutine upinvhop_fnf(nelorb, nelc, indt, ainv               &
        &, winv, winvup, winvdo, psiln, psisn, epst, psip, nelup, neldo             &
        &, ainvs, winvs, psinew, leadpsi, detmat, indtm, nelorbh&
        &, mu_c, nelorb_c, firstmol, nmol, yesfast, winvbarfn, winvfn&
        &, detmat_c, contraction, indt4r, no_dgemv)

    use Constants, only: ip4, ipc, ipf, zzero, zone, zimg, ip_reshuff
    use allio, only: rank, nelup_mat, nel_mat, ndiff, nosingledet, molecular, iscramax, yes_ontarget

    implicit none

    integer Nel, j, kk, nelc, nelcdo, nelorb, nelup, neldo, leadpsi           &
            &, nelorb4, k, indt, nelorbh, yesfast, indt4, indt4r          &
            &, indt5, i, indtm(nelup + neldo), firstmol, nmol, nelorb_c, contraction, nmoltest, nmolipf, nmolshift
    real*8 winv(ipc*nelorb, 0:indt4r, *), ainv(ipc*nelup_mat, nelup_mat)&
            &, epst, psip(ipc*nmol, 6), ainvs(ipc*nelup_mat), winvs(ipc*(indt4r + 1)*nelorb) &
            !      &     ,winvup(ipc*nelup,*),winvdo(max(ipc*neldo,ipc),*),psinew(ipc*nmol/ipf,0:*)     &
            &, winvup(ipc*nelup, indt + ip4), winvdo(max(ipc*neldo, ipc), indt + ip4)&
            &, psinew(leadpsi, 0:indt + ip4), detmat(ipc*ipf*nelorbh, *), psiln&
            &, mu_c(ipc*ipf*nelorbh, *), winvbarfn(ipc*nmol, nel_mat)               &
            &, winvfn(leadpsi, indt + ip4, nelup + neldo), detmat_c(ipc*nelorb_c, *), psisn, costpip(2)
    real*8 gr
    complex*16 gc
    logical no_dgemv

    ! The scratch psip has to have at least dimension: 5*nelorb+max(nelorb,nelup*(indt+5))
    ! Calculation prefactor:
    !   nelc is the index of the electrons that has done the move.
    !   nelc<=nelup up spin, otherwise down spin
    ! Here winv(nelorb,nel) are the wavefunction orbitals

    ! Compute the wavefunction in the new position r and update new functions
    if (ipf .eq. 2 .and. molecular .gt. 0) then
        nmolipf = nmol
        nmolshift = 0
    else
        nmolipf = nmol/ipf
        nmolshift = (ipf - 1)*nmol/2
    end if

    nel = neldo + nelup
    !          nelpair=nelorb-neldiff
    indt4 = indt + ip4
    indt5 = indt4 + 1
    nelorb4 = nmol*indt4
#ifdef _OFFLOAD
!$omp target data map(winvup,winvdo) if(yes_ontarget)
#endif
    !          firstmol=nelorb_c-molecular+1
    ! psip(1:ipc*nmol,1:6)=0.d0  ! clean up the work space
    call dscalzero__(6*ipc*nmol, 0.d0, psip, 1)
    if (nelc .le. nelup .or. ipf .eq. 2) then ! up-spin electrons
        ! ainvs= \Phi_i
        if (ipc .eq. 1) then

            call dgemv_('T', nelup_mat, nelup_mat, 1.d0, ainv, nelup_mat, ainvs, 1 &
                    &, 0.d0, psip(1, 2), 1, yes_ontarget)
#ifdef _OFFLOAD
!$omp target update from (psip(nelc:nelc,2:2))  if(yes_ontarget)
#endif
            gr = psip(nelc, 2)
            psip(nelc, 2) = psip(nelc, 2) - 1.d0

            if (gr .ne. 0.d0) then
                psisn = psisn*int(sign(1.d0, gr))
                psiln = psiln + log(abs(gr))
                gr = -1.d0/gr
            end if

#ifdef _OFFLOAD
!$omp target update to (psip(nelc:nelc,2:2)) if(yes_ontarget)
#endif

        else

            call zgemv_('T', nelup_mat, nelup_mat, zone, ainv, nelup_mat, ainvs, 1 &
                    &, zzero, psip(1, 2), 1, yes_ontarget)
#ifdef _OFFLOAD
!$omp target update from(psip(2*nelc-1:2*nelc,2:2)) if(yes_ontarget)
#endif
            gc = dcmplx(psip(2*nelc - 1, 2), psip(2*nelc, 2))
            psip(2*nelc - 1, 2) = psip(2*nelc - 1, 2) - 1.d0
            call makeg_complex(gc, psiln, psisn)
#ifdef _OFFLOAD
!$omp target update to(psip(2*nelc-1:2*nelc-1,2:2)) if(yes_ontarget)
#endif

        end if

        !    psip(1:ipc*nelup_mat,1)=ainv(1:ipc*nelup_mat,nelc)
        call dcopy_vec(ipc*nelup_mat, ainv(1, nelc), psip)

        if (ipf .eq. 2) then
            !    vector a
            !    vector b parbcs.pdf
            if (ipc .eq. 2) then
                gc = -gc
                !      call zscal(nelup_mat,g,psip,1)
                call scalepsip(nelup_mat, gc, psip)
                !      psip(1:2*nelup_mat,2)=-psip(1:2*nelup_mat,2)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, 2*nelup_mat
                    psip(i, 2) = -psip(i, 2)
                end do

            else
                !      psip(1:nelup_mat,1)=-g(1)*psip(1:nelup_mat,1)
                !      psip(1:nelup_mat,2)=-psip(1:nelup_mat,2)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup_mat
                    psip(i, 1) = -gr*psip(i, 1)
                    psip(i, 2) = -psip(i, 2)
                end do
            end if
            ! Rank-2 update
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
            do i = 1, ipc*nelup_mat
                psip(i, 3) = psip(i, 1)
                psip(i, 4) = -psip(i, 2)
                psip(i, 5) = psip(i, 2)
                psip(i, 6) = psip(i, 1)
                ainvs(i) = -ainv(i, nelc)
            end do
            if (ipc .eq. 1) then
                !      call dgemm('N','T',nelup_mat,nelup_mat,2,1.d0,psip(1,3),nmol,psip(1,5),nmol,1.d0,ainv,nelup_mat)
                call dger2(nelup_mat, nelup_mat, 1.d0, psip(1, 3), nmol, psip(1, 5), nmol, ainv, nelup_mat)
                !      definition of abar bbar postponed after updating winvbarfn
            else
                !      call zgemm('N','T',nelup_mat,nelup_mat,2,zone,psip(1,3),nmol,psip(1,5),nmol,zone,ainv,nelup_mat)
                call zger2(nelup_mat, nelup_mat, zone, psip(1, 3), nmol, psip(1, 5), nmol, ainv, nelup_mat)
            end if

        else

            if (ipc .eq. 1) then

                call dger_(nelup, nelup, gr, psip, 1, psip(1, 2), 1, ainv, nelup)

                call dgemv_('N', nmol, nelup, 1.d0, winvbarfn(1, nelup + 1), nmol   &
                        &, psip, 1, 0.d0, psip(1, 3), 1, yes_ontarget)

                ! now update v (psip(*,4)) and w (psip(*,5)) for the down
                call dgemv_('N', nmol, nelup, 1.d0, winvbarfn, nmol        &
                        &, psip(1, 2), 1, 0.d0, psip(1, 4), 1, yes_ontarget)

            else

                call zgeru_(nelup, nelup, gc, psip, 1, psip(1, 2), 1, ainv, nelup)

                call zgemv_('N', nmol, nelup, zone, winvbarfn(1, nelup + 1), nmol   &
                        &, psip, 1, zzero, psip(1, 3), 1, yes_ontarget)

                call zgemv_('N', nmol, nelup, zone, winvbarfn, nmol     &
                        &, psip(1, 2), 1, zzero, psip(1, 4), 1, yes_ontarget)

            end if

        end if ! ned ipf=2

        ! update winvbar here
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
        do j = 1, nmol*ipc
            psip(j, 5) = -winvbarfn(j, nelc)
        end do
        if (ipc .eq. 1) then
            if (yesfast .ne. 0) then
                if (ipf .eq. 2) then
                    if (nelc .le. nelup) then
                        call dgemv_('N', nmol, nmolipf, 1.d0, detmat_c(firstmol, firstmol), nelorb_c&
                                &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                    else
                        call dgemv_('N', nmol, nmolipf, 1.d0, detmat_c(firstmol, firstmol + nmolshift), nelorb_c&
                                &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                    end if
                else
                    call dgemv_('T', nmol, nmol, 1.d0, detmat_c(firstmol, firstmol), nelorb_c&
                            &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                end if
            else
                if (contraction .ne. 0) then
                    if (ipf .eq. 2) then
                        if (nelc .le. nelup) then
                            call dgemv_('N', nmol, nmolipf, 1.d0, detmat_c(firstmol, firstmol), nelorb_c&
                                    &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                        else
                            call dgemv_('N', nmol, nmolipf, 1.d0, detmat_c(firstmol, firstmol + nmolshift), nelorb_c&
                                    &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                        end if
                    else
                        call dgemv_('T', nmol, nmol, 1.d0, detmat_c(firstmol, firstmol), nelorb_c&
                                &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                    end if
                else
                    if (ipf .eq. 2) then
                        if (nelc .le. nelup) then
                            call dgemv_('N', 2*nelorbh, nelorbh, 1.d0, detmat, 2*nelorbh            &
                                    &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                        else
                            call dgemv_('N', 2*nelorbh, nelorbh, 1.d0, detmat(1, nelorbh + 1), 2*nelorbh&
                                    &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                        end if
                    else
                        call dgemv_('T', nelorbh, nelorbh, 1.d0, detmat, nelorbh&
                                &, winvs, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                    end if
                end if
            end if

        else

            if (yesfast .ne. 0) then
                if (ipf .eq. 2) then
                    if (nelc .le. nelup) then
                        call zgemv_('N', nmol, nmolipf, zone, detmat_c(2*firstmol - 1, firstmol), nelorb_c&
                                &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                    else
                        call zgemv_('N', nmol, nmolipf, zone, detmat_c(2*firstmol - 1, firstmol + nmolshift), nelorb_c&
                                &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                    end if
                else
                    call zgemv_('T', nmol, nmol, zone, detmat_c(2*firstmol - 1, firstmol), nelorb_c&
                            &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                end if
            else
                if (contraction .ne. 0) then
                    if (ipf .eq. 2) then
                        if (nelc .le. nelup) then
                            call zgemv_('N', nmol, nmolipf, zone, detmat_c(2*firstmol - 1, firstmol), nelorb_c&
                                    &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                        else
                            call zgemv_('N', nmol, nmolipf, zone, detmat_c(2*firstmol - 1, firstmol + nmolshift), nelorb_c&
                                    &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                        end if
                    else
                        call zgemv_('T', nmol, nmol, zone, detmat_c(2*firstmol - 1, firstmol), nelorb_c&
                                &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                    end if
                else
                    if (ipf .eq. 2) then
                        if (nelc .le. nelup) then
                            call zgemv_('N', 2*nelorbh, nelorbh, zone, detmat, 2*nelorbh            &
                                    &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                        else
                            call zgemv_('N', 2*nelorbh, nelorbh, zone, detmat(1, nelorbh + 1), 2*nelorbh            &
                                    &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                        end if
                    else
                        call zgemv_('T', nelorbh, nelorbh, zone, detmat, nelorbh            &
                                &, winvs, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                    end if
                end if
            end if
        end if

        if (ipf .eq. 1) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
            do j = 1, ipc*nmol
                psip(j, 5) = psip(j, 5) + winvbarfn(j, nelc)
            end do
        end if

        !         now update winv (all if indt4>0)
        !    call dcopy(ipc*nelorb*(indt4r+1),winvs,1,winv(1,0,nelc),1)
        call dcopy_vec(ipc*nelorb*(indt4r + 1), winvs, winv(1, 0, nelc))

        !    winvfn(:,1:indt4,nelc)=psinew(:,1:indt4)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
        do i = 1, indt4
            do j = 1, leadpsi
                winvfn(j, i, nelc) = psinew(j, i)
            end do
        end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
        if (ipf .eq. 2) then
            !    compute abar bbar Eq.K12 parbcs
            if (ipc .eq. 1) then
                call dgemv_('N', nmol, nel_mat, 1.d0, winvbarfn, nmol   &
                        &, psip, 1, 0.d0, psip(1, 3), 1, yes_ontarget) ! abar
                call dgemv_('N', nmol, nel_mat, 1.d0, winvbarfn, nmol        &
                        &, psip(1, 2), 1, 0.d0, psip(1, 4), 1, yes_ontarget) !bbar
            else
                call zgemv_('N', nmol, nel_mat, zone, winvbarfn, nmol   &
                        &, psip, 1, zzero, psip(1, 3), 1, yes_ontarget) ! abar
                call zgemv_('N', nmol, nel_mat, zone, winvbarfn, nmol     &
                        &, psip(1, 2), 1, zzero, psip(1, 4), 1, yes_ontarget) !bbar
            end if
#ifdef _OFFLOAD
!$omp target teams distribute parallel do  if(yes_ontarget)
#endif
            do j = 1, ipc*nmol
                psip(j, 5) = psip(j, 5) + winvbarfn(j, nelc)
            end do

            if (ipc .eq. 1) then
                !!       contract O_mu winvfn with abar put in psip(1,6)
                !        call  dgemvtry(nmol,nel,indt,ip4,indtm,winvfn(ip_reshuff,1,1)&
                !    &,psip(ip_reshuff,3),ip_reshuff,psip(1,6))
                !  call upwinv_pfaff(nelup,neldo,indt4,winvup,winvdo,psip(ip_reshuff,2),psip(1,6),.true.)

                !!       contract O_mu winvfn with bbar put in psip(1,6)
                !        call  dgemvtry(nmol,nel,indt,ip4,indtm,winvfn(ip_reshuff,1,1)&
                !    &,psip(ip_reshuff,4),ip_reshuff,psip(1,6))

                !  call upwinv_pfaff(nelup,neldo,indt4,winvup,winvdo,psip(ip_reshuff,1),psip(1,6),.false.)
                !
                !!       contract O_mu winvfn with psip(:,5)  put in psip(1,6)
                !        call  dgemvtry(nmol,nel,indt,ip4,indtm,winvfn(ip_reshuff,1,1)&
                !    &,psip(ip_reshuff,5),ip_reshuff,psip(1,6))
                !
                !  call upwinv_pfaff(nelup,neldo,indt4,winvup,winvdo,ainvs(ip_reshuff),psip(1,6),.true.)
                !
                call dgemvtry_pfaff(nmol, nmolipf, nmolshift, nelup, neldo, nel, indt, ip4, indtm, winvfn&
                        &, psip(1, 3), psip, ainvs, winvup, winvdo)

                !   computing ainvfn(nelc) put in psip(1,6)
                call dgemv_('N', nmol, nelup_mat, 1.d0, winvbarfn, nmol&
                    &, ainv(1, nelc), 1, 0.d0, psip(1, 6), 1, yes_ontarget)
                call upwinvp_pfaff(nelc, nelup, neldo, nmol, nmolipf, nmolshift&
                    &, indt4, winvup, winvdo, winvfn(1, 1, nelc), psip(1, 6))
            else
                !!       contract O_mu winvfn with abar put in psip(1,6)
                !        call  zgemvtry(nmol,nel,indt,ip4,indtm,winvfn(2*ip_reshuff-1,1,1)&
                !     &,psip(2*ip_reshuff-1,3),ip_reshuff,psip(1,6))
                !  call upwinv_pfaff_complex(nelup,neldo,indt4,winvup,winvdo,psip(2*ip_reshuff-1,2),psip(1,6),.true.)
                !
                !!       contract O_mu winvfn with bbar put in psip(1,6)
                !        call  zgemvtry(nmol,nel,indt,ip4,indtm,winvfn(2*ip_reshuff-1,1,1)&
                !     &,psip(2*ip_reshuff-1,4),ip_reshuff,psip(1,6))
                !  call upwinv_pfaff_complex(nelup,neldo,indt4,winvup,winvdo,psip(2*ip_reshuff-1,1),psip(1,6),.false.)
                !!       contract O_mu winvfn with psip(:,5)  put in psip(1,6)
                !        call  zgemvtry(nmol,nel,indt,ip4,indtm,winvfn(2*ip_reshuff-1,1,1)&
                !    &,psip(2*ip_reshuff-1,5),ip_reshuff,psip(1,6))
                !  call upwinv_pfaff(nelup,neldo,indt4,winvup,winvdo,ainvs(2*ip_reshuff-1),psip(1,6),.true.)

 call zgemvtry_pfaff(nmol, nmolipf, nmolshift, nelup, neldo, nel, indt, ip4&
     &, indtm, winvfn, psip(1, 3), psip, ainvs, winvup, winvdo)
                !   computing ainvfn (nelc)
                call zgemv_('N', nmol, nelup_mat, zone, winvbarfn, nmol, ainv(1, nelc), 1, zzero, psip(1, 6), 1, yes_ontarget)
     call upwinvp_pfaff_complex(nelc, nelup, neldo, nmol, nmolipf, nmolshift&
         &, indt4, winvup, winvdo, winvfn(1, 1, nelc), psip(1, 6))
                !
            end if !end ipc=2
        else ! ipf=2
            !         now compute psi_j^nu for winvup
            if (ipc .eq. 1) then
                call dgemvtry(nmol, nelup, indt, ip4, indtm, winvfn(ip_reshuff, 1, 1)&
                        &, psip(ip_reshuff, 3), ip_reshuff, psip(1, 6))
            else
                call zgemvtry(nmol, nelup, indt, ip4, indtm, winvfn(2*ip_reshuff - 1, 1, 1)&
                        &, psip(2*ip_reshuff - 1, 3), ip_reshuff, psip(1, 6))
            end if

            ! v  for the up
            if (ipc .eq. 1) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do  if(yes_ontarget)
#endif
                do i = 1, nelup
                    psip(i, 2) = gr*psip(i, 2)
                end do
                call upwinv(nelup, nelc, indt4, nelup, winvup, psip(1, 2), psip(1, 6))
            else
                call scalepsip(nelup, gc, psip(1, 2))
                call upwinv_complex(nelup, nelc, indt4, nelup, winvup, psip(1, 2), psip(1, 6))
            end if

            if (nosingledet) then ! simple Slater determinant not necessary to update below

                !         define v for the down
                !    psip(1:nelup*ipc,3)=ainv(1:nelup*ipc,nelc)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, ipc*nelup
                    psip(i, 3) = ainv(i, nelc)
                end do

                if (ipc .eq. 1) then
                    !       psip(1:nelup,1)=g(1)*psip(1:nelup,1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                    do i = 1, nelup
                        psip(i, 1) = gr*psip(i, 1)
                    end do
                else
                    call scalepsip(nelup, gc, psip)
                end if

                !         now compute new psi and psi bar for winvdo
                !       Here  for ip_reshuff=2   psip(:,4) has odd  non zero components (even are zero)
                if (ipc .eq. 1) then
                    call dgemmtry(neldo, indt, ip4, indtm(nelup + 1), nmol&
      &, winvfn(1, 1, nelup + 1), psip(1, 4), ip_reshuff, psip(1, 6), no_dgemv)
                    call upwinvp(neldo, indt4, winvdo, psip, psip(1, 3), psip(1, 6))
                else
                    call zgemmtry(neldo, indt, ip4, indtm(nelup + 1), nmol            &
                            &, winvfn(1, 1, nelup + 1), psip(1, 4), ip_reshuff, psip(1, 6), no_dgemv)
                    call upwinvp_complex(neldo, indt4, winvdo, psip, psip(1, 3), psip(1, 6))
                end if
            end if ! neldo ne nmol

        end if ! endif ipf=2

    else ! down-spin electrons and ipf=1

        nelcdo = nelc - nelup
        if (ipc .eq. 1) then

            call dgemv_('N', nelup, nelup, 1.d0, ainv, nelup, ainvs, 1, 0.d0       &
                    &, psip(1, 2), 1, yes_ontarget)
#ifdef _OFFLOAD
!$omp target update from (psip(nelcdo:nelcdo,2:2)) if(yes_ontarget)
#endif
            gr = psip(nelcdo, 2)
            psip(nelcdo, 2) = psip(nelcdo, 2) - 1.d0
            ! protect from non accurate machines
            if (gr .ne. 0.d0) then
                psisn = psisn*int(sign(1.d0, gr))
                psiln = psiln + log(abs(gr))
                gr = -1.d0/gr
            end if

#ifdef _OFFLOAD
!$omp target update to (psip(nelcdo:nelcdo,2:2)) if(yes_ontarget)
#endif

        else

            call zgemv_('N', nelup, nelup, zone, ainv, nelup, ainvs, 1, zzero       &
                    &, psip(1, 2), 1, yes_ontarget)

#ifdef _OFFLOAD
!$omp target update from(psip(2*nelcdo-1:2*nelcdo,2:2)) if(yes_ontarget)
#endif
            gc = dcmplx(psip(2*nelcdo - 1, 2), psip(2*nelcdo, 2))
            psip(2*nelcdo - 1, 2) = psip(2*nelcdo - 1, 2) - 1.d0

            call makeg_complex(gc, psiln, psisn)
#ifdef _OFFLOAD
!$omp target update to (psip(2*nelcdo-1:2*nelcdo-1,2:2)) if(yes_ontarget)
#endif

        end if

        if (ipc .eq. 1) then

            !       psip(1:nelup,1)=ainv(nelcdo,1:nelup)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do  if(yes_ontarget)
#endif
            do i = 1, nelup
                psip(i, 1) = ainv(nelcdo, i)
            end do

            call dger_(nelup, nelup, gr, psip(1, 2), 1, psip, 1, ainv, nelup)

            !         now update w (psip(*,3)  and v psip(*,2))  for down
            call dgemv_('N', nmol, nelup, 1.d0, winvbarfn, nmol              &
                    &, psip, 1, 0.d0, psip(1, 3), 1, yes_ontarget)

            call dgemv_('N', nmol, nelup, 1.d0, winvbarfn(1, nelup + 1), nmol  &
                    &, psip(1, 2), 1, 0.d0, psip(1, 4), 1, yes_ontarget)

        else
            call up_ainv_down(psip, ainv, nelup, nelcdo)
            call zgeru_(nelup, nelup, gc, psip(1, 2), 1, psip, 1, ainv, nelup)
            !         now update w (psip(*,3)  and v psip(*,2))  for down
            call zgemv_('N', nmol, nelup, zone, winvbarfn, nmol              &
                    &, psip, 1, zzero, psip(1, 3), 1, yes_ontarget)
            !         psip(1,3)= w
            !    now update v (psip(*,5))  and w =bar bar \Phi   (psip(*,4)) for the
            call zgemv_('N', nmol, nelup, zone, winvbarfn(1, nelup + 1), nmol  &
                    &, psip(1, 2), 1, zzero, psip(1, 4), 1, yes_ontarget)
        end if

        ! update winvbar here
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
        do j = 1, ipc*nmol
            psip(j, 5) = -winvbarfn(j, nelc)
        end do

        if (ipc .eq. 1) then

            if (yesfast .ne. 0) then
                call dgemv_('N', nmol, nmol, 1.d0, detmat_c(firstmol, firstmol), nelorb_c&
                        &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)

            else
                if (contraction .ne. 0) then
                    call dgemv_('N', nmol, nmol, 1.d0, detmat_c(firstmol, firstmol), nelorb_c&
                            &, psinew, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                else
                    call dgemv_('N', nelorbh, nelorbh, 1.d0, detmat, nelorbh            &
                            &, winvs, 1, 0.d0, winvbarfn(1, nelc), 1, yes_ontarget)
                end if
            end if

        else

            if (yesfast .ne. 0) then
                call zgemv_('N', nmol, nmol, zone, detmat_c(2*firstmol - 1, firstmol), nelorb_c&
                        &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)

            else
                if (contraction .ne. 0) then
                    call zgemv_('N', nmol, nmol, zone, detmat_c(2*firstmol - 1, firstmol), nelorb_c&
                            &, psinew, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                else
                    call zgemv_('N', nelorbh, nelorbh, zone, detmat, nelorbh            &
                            &, winvs, 1, zzero, winvbarfn(1, nelc), 1, yes_ontarget)
                end if
            end if

        end if

#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
        do j = 1, ipc*nmol
            psip(j, 5) = psip(j, 5) + winvbarfn(j, nelc)
        end do

        ! now update winv

        !    call dcopy(ipc*nelorb*(indt4r+1),winvs,1,winv(1,0,nelc),1)
        call dcopy_vec(ipc*nelorb*(indt4r + 1), winvs, winv(1, 0, nelc))
        !    winvfn(1:nmol*ipc,1:indt4,nelc)=psinew(1:nmol*ipc,1:indt4)

#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
        do i = 1, nmol*ipc
            winvfn(i, 1:indt4, nelc) = psinew(i, 1:indt4)
        end do

        if (neldo .gt. 0) then
            ! now compute psi_j^nu for winvdo
            if (ipc .eq. 1) then
                !       Here for ip_reshuff=2 psip(:,3) has only odd components (even are zero)
                call dgemvtry(nmol, neldo, indt, ip4, indtm(nelup + 1)           &
                        &, winvfn(1, 1, nelup + 1), psip(1, 3), ip_reshuff, psip(1, 6))
                !       psip(1:nelup,2)=g(1)*psip(1:nelup,2)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup
                    psip(i, 2) = gr*psip(i, 2)
                end do

                call upwinv(neldo, nelcdo, indt4, neldo, winvdo, psip(1, 2)         &
                        &, psip(1, 6))
            else
                call zgemvtry(nmol, neldo, indt, ip4, indtm(nelup + 1)           &
                        &, winvfn(1, 1, nelup + 1), psip(1, 3), ip_reshuff, psip(1, 6))
                call scalepsip(nelup, gc, psip(1, 2))
                call upwinv_complex(neldo, nelcdo, indt4, neldo, winvdo, psip(1, 2)         &
                        &, psip(1, 6))
            end if
        end if

        if (nosingledet) then ! Slater determinant not necessary steps below

            if (ipc .eq. 1) then
                !       psip(1:nelup,3)=ainv(nelcdo,1:nelup)
                !       psip(1:nelup,1)=g(1)*psip(1:nelup,1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup
                    psip(i, 3) = ainv(nelcdo, i)
                    psip(i, 1) = gr*psip(i, 1)
                end do
            else
                call up_ainv_down(psip(1, 3), ainv, nelup, nelcdo)
                call scalepsip(nelup, gc, psip)
            end if

            ! now compute new psi and psi bar for winvup
            !       Here for ip_reshuff=2 psip(:,4) has only even components (odd are zero)
            if (ipc .eq. 1) then
                call dgemmtry(nelup, indt, ip4, indtm, nmol&
     &, winvfn(ip_reshuff, 1, 1), psip(ip_reshuff, 4), ip_reshuff, psip(1, 6), no_dgemv)
                call upwinvp(nelup, indt4, winvup, psip, psip(1, 3), psip(1, 6))
            else
                call zgemmtry(nelup, indt, ip4, indtm, nmol&
                        &, winvfn(2*ip_reshuff - 1, 1, 1), psip(2*ip_reshuff - 1, 4), ip_reshuff, psip(1, 6), no_dgemv)
                call upwinvp_complex(nelup, indt4, winvup, psip, psip(1, 3), psip(1, 6))
            end if

        end if ! Slater determinant if

    end if
#ifdef _OFFLOAD
!$omp end target data
#endif
    return
end subroutine upinvhop_fnf

subroutine dgemmtry(neldo, indt, ip4, indtm, nelorb              &
        &, winv, psiin, ip_reshuff, psiout, no_dgemv)
    use constants, only: yes_ontarget
    implicit none
    integer neldo, indt, indtm(neldo), nelorb, i, j, k, ip4, ip_reshuff
    real*8 winv(nelorb, indt + ip4, neldo), psiin(nelorb, 2)              &
            &, psiout(indt + ip4, neldo, 2)
    real*8 csum1, csum2
    logical no_dgemv
    integer indta
    if (.not. no_dgemv .and. ip_reshuff .eq. 1) then
        indta = (indt + ip4)*neldo
        if (yes_ontarget) then
            call dgemm_('T', 'N', indta, 2, nelorb, 1.d0, winv, nelorb, psiin, nelorb, 0.d0, psiout, indta)
        else
            call dgemm('T', 'N', indta, 2, nelorb, 1.d0, winv, nelorb, psiin, nelorb, 0.d0, psiout, indta)
        end if
    else
        if (yes_ontarget) then
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
            do i = 1, neldo
                do j = 1, indt + ip4
                    psiout(j, i, 1) = 0.d0
                    psiout(j, i, 2) = 0.d0
                end do
            end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(csum1,csum2)
#else
!$omp parallel do default(shared) private(i,k,j,csum1,csum2)
#endif
            do i = 1, neldo
                do j = 1, indt + ip4
                    if (j .le. indtm(i) .or. j .gt. indt) then
                        csum1 = psiout(j, i, 1)
                        csum2 = psiout(j, i, 2)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum1,csum2)
#endif
                        do k = 1, nelorb, ip_reshuff
                            csum1 = csum1 + winv(k, j, i)*psiin(k, 1)
                            csum2 = csum2 + winv(k, j, i)*psiin(k, 2)
                        end do
                        psiout(j, i, 1) = csum1
                        psiout(j, i, 2) = csum2
                    end if
                end do
            end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif
        else
            do i = 1, neldo
                do j = 1, indt + ip4
                    psiout(j, i, 1) = 0.d0
                    psiout(j, i, 2) = 0.d0
                end do
            end do
!$omp parallel do default(shared) private(i,k,j,csum1,csum2)
            do i = 1, neldo
                do j = 1, indt + ip4
                    if (j .le. indtm(i) .or. j .gt. indt) then
                        csum1 = psiout(j, i, 1)
                        csum2 = psiout(j, i, 2)
                        do k = 1, nelorb, ip_reshuff
                            csum1 = csum1 + winv(k, j, i)*psiin(k, 1)
                            csum2 = csum2 + winv(k, j, i)*psiin(k, 2)
                        end do
                        psiout(j, i, 1) = csum1
                        psiout(j, i, 2) = csum2
                    end if
                end do
            end do
!$omp end parallel do
        end if
    end if ! no no_dgemv
    return
end subroutine dgemmtry
subroutine zgemmtry(neldo, indt, ip4, indtm, nelorb              &
        &, winv, psiin, ip_reshuff, psiout, no_dgemv)
    use constants, only: yes_ontarget, zone, zzero
    implicit none
    integer neldo, indt, indtm(neldo), nelorb, i, j, k, ip4, ip_reshuff
    complex*16 winv(nelorb, indt + ip4, neldo), psiin(nelorb, 2)              &
            &, psiout(indt + ip4, neldo, 2)
    complex*16 csum1, csum2
    integer indta
    logical no_dgemv
    if (.not. no_dgemv .and. ip_reshuff .eq. 1) then
        indta = (indt + ip4)*neldo
        if (yes_ontarget) then
            call zgemm_('T', 'N', indta, 2, nelorb, zone, winv, nelorb, psiin, nelorb, zzero, psiout, indta)
        else
            call zgemm('T', 'N', indta, 2, nelorb, zone, winv, nelorb, psiin, nelorb, zzero, psiout, indta)
        end if
    else
        if (yes_ontarget) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
            do i = 1, neldo
                do j = 1, indt + ip4
                    psiout(j, i, 1) = (0.d0, 0.d0)
                    psiout(j, i, 2) = (0.d0, 0.d0)
                end do
            end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(csum1,csum2)
#else
!$omp parallel do default(shared) private(i,k,j,csum1,csum2)
#endif
            do i = 1, neldo
                do j = 1, indt + ip4
                    if (j .le. indtm(i) .or. j .gt. indt) then
                        csum1 = psiout(j, i, 1)
                        csum2 = psiout(j, i, 2)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum1,csum2)
#endif
                        do k = 1, nelorb, ip_reshuff
                            csum1 = csum1 + winv(k, j, i)*psiin(k, 1)
                            csum2 = csum2 + winv(k, j, i)*psiin(k, 2)
                        end do
                        psiout(j, i, 1) = csum1
                        psiout(j, i, 2) = csum2
                    end if
                end do
            end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif
        else
            do i = 1, neldo
                do j = 1, indt + ip4
                    psiout(j, i, 1) = (0.d0, 0.d0)
                    psiout(j, i, 2) = (0.d0, 0.d0)
                end do
            end do

!$omp parallel do default(shared) private(i,k,j,csum1,csum2)
            do i = 1, neldo
                do j = 1, indt + ip4
                    if (j .le. indtm(i) .or. j .gt. indt) then
                        csum1 = psiout(j, i, 1)
                        csum2 = psiout(j, i, 2)
!$omp parallel do reduction(+:csum1,csum2)
                        do k = 1, nelorb, ip_reshuff
                            csum1 = csum1 + winv(k, j, i)*psiin(k, 1)
                            csum2 = csum2 + winv(k, j, i)*psiin(k, 2)
                        end do
                        psiout(j, i, 1) = csum1
                        psiout(j, i, 2) = csum2
                    end if
                end do
            end do
!$omp end parallel do
        end if
    end if
    return
end subroutine zgemmtry
subroutine dgemvtry(nelorb, nelup, indt, ip4, indtm, winv, psiin, ip_reshuff, psiout)
    use constants, only: yes_ontarget
    implicit none
    integer nelorb, nelup, indt, indtm(nelup), i, j, k, ip4, ip_reshuff
    real*8 winv(nelorb, 1:indt + ip4, nelup), psiin(nelorb)                  &
            &, psiout(indt + ip4, nelup)
    real*8 csum
    if (yes_ontarget) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
        do i = 1, nelup
            do j = 1, indt + ip4
                psiout(j, i) = 0.d0
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(csum)
#else
!$omp parallel do default(shared) private(i,k,j,csum)
#endif
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiout(j, i)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do k = 1, nelorb, ip_reshuff
                        csum = csum + winv(k, j, i)*psiin(k)
                    end do
                    psiout(j, i) = csum
                end if
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif
    else
        do i = 1, nelup
            do j = 1, indt + ip4
                psiout(j, i) = 0.d0
            end do
        end do
!$omp parallel do default(shared) private(i,k,j,csum)
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiout(j, i)
                    do k = 1, nelorb, ip_reshuff
                        csum = csum + winv(k, j, i)*psiin(k)
                    end do
                    psiout(j, i) = csum
                end if
            end do
        end do
!$omp end parallel do
    end if
    return
end subroutine dgemvtry

subroutine dgemvtry_pfaff(nelorb, nmolipf, nmolshift                     &
           &, nelup, neldo, nel, indt, ip4, indtm, winv                  &
           &, barpsiin, psiin, ainvs, psiup, psido)
    use constants, only: yes_ontarget
    implicit none
    integer nelorb, nmolipf, nmolshift, nelup, neldo, nel, indt, indtm(nel), i, j, k, ip, ip4
    real*8 winv(nmolipf, 1:indt + ip4, nel), barpsiin(nelorb, 3), psiin(nelorb, 2)&
            &, psiup(nelup, indt + ip4), psido(max(neldo, 1), indt + ip4), ainvs(nel)
    real*8 csum
    if (yes_ontarget) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
        do i = 1, nelup
            do j = 1, indt
                if (j .gt. indtm(i)) then
                    psiup(i, j) = 0.d0
                end if
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) private(ip)
#endif
        do i = 1, neldo
            do j = 1, indt
                ip = i + nelup
                if (j .gt. indtm(ip)) psido(i, j) = 0.d0
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(csum)
#else
!$omp parallel do default(shared) private(i,k,j,csum)
#endif
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiup(i, j)
#ifdef  _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do k = 1, nmolipf
                        csum = csum + winv(k, j, i)*(psiin(i, 2)*barpsiin(k, 1)&
                           & - psiin(i, 1)*barpsiin(k, 2) + ainvs(i)*barpsiin(k, 3))
                    end do
                    psiup(i, j) = csum
                end if
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(ip,csum)
#else
!$omp parallel do default(shared) private(i,ip,k,j,csum)
#endif
        do i = 1, neldo
            do j = 1, indt + ip4
                ip = i + nelup
                if (j .le. indtm(ip) .or. j .gt. indt) then
                    csum = psido(i, j)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do k = nmolshift + 1, nelorb
                        csum = csum + winv(k - nmolshift, j, ip)*(psiin(ip, 2)&
                               &* barpsiin(k, 1) - psiin(ip, 1)*barpsiin(k, 2)&
                               & + ainvs(ip)*barpsiin(k, 3))
                    end do
                    psido(i, j) = csum
                end if
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif
    else
        do i = 1, nelup
            do j = 1, indt
                if (j .gt. indtm(i)) then
                    psiup(i, j) = 0.d0
                end if
            end do
        end do

        do i = 1, neldo
            do j = 1, indt
                ip = i + nelup
                if (j .gt. indtm(ip)) psido(i, j) = 0.d0
            end do
        end do

!$omp parallel do default(shared) private(i,k,j,csum)
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiup(i, j)
                    do k = 1, nmolipf
                        csum = csum + winv(k, j, i)*(psiin(i, 2)*barpsiin(k, 1)&
                           & - psiin(i, 1) * barpsiin(k, 2) + ainvs(i)*barpsiin(k, 3))
                    end do
                    psiup(i, j) = csum
                end if
            end do
        end do
!$omp end parallel do

!$omp parallel do default(shared) private(i,ip,k,j,csum)
        do i = 1, neldo
            do j = 1, indt + ip4
                ip = i + nelup
                if (j .le. indtm(ip) .or. j .gt. indt) then
                    csum = psido(i, j)
                    do k = nmolshift + 1, nelorb
                        csum = csum + winv(k - nmolshift, j, ip)*(psiin(ip, 2)&
                            & * barpsiin(k, 1) - psiin(ip, 1) * barpsiin(k, 2)&
                            & + ainvs(ip)*barpsiin(k, 3))
                    end do
                    psido(i, j) = csum
                end if
            end do
        end do
!$omp end parallel do
    end if
    return
end subroutine dgemvtry_pfaff
subroutine zgemvtry_pfaff(nelorb, nmolipf, nmolshift                     &
           &, nelup, neldo, nel, indt, ip4, indtm, winv                  &
           &, barpsiin, psiin, ainvs, psiup, psido)
    use constants, only: yes_ontarget
    implicit none
    integer nelorb, nelup, neldo, nel, indt, indtm(nel), i, j, k, ip, ip4, nmolipf, nmolshift
    complex*16 winv(nmolipf, 1:indt + ip4, nel), barpsiin(nelorb, 3), psiin(nelorb, 2)&
            &, psiup(nelup, indt + ip4), psido(max(neldo, 1), indt + ip4), ainvs(nel)
    complex*16 csum
    if (yes_ontarget) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
        do i = 1, nelup
            do j = 1, indt
                if (j .gt. indtm(i)) psiup(i, j) = (0.d0, 0.d0)
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) private(ip)
#endif
        do i = 1, neldo
            do j = 1, indt
                ip = i + nelup
                if (j .gt. indtm(ip)) psido(i, j) = (0.d0, 0.d0)
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(csum)
#else
!$omp parallel do default(shared) private(i,k,j,csum)
#endif
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiup(i, j)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do k = 1, nmolipf
                     csum = csum + winv(k, j, i)*(psiin(i, 2)*barpsiin(k, 1)&
                        & - psiin(i, 1)*barpsiin(k, 2) + ainvs(i)*barpsiin(k, 3))
                    end do
                    psiup(i, j) = csum
                end if
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(ip,csum)
#else
!$omp parallel do default(shared) private(i,ip,k,j,csum)
#endif
        do i = 1, neldo
            do j = 1, indt + ip4
                ip = i + nelup
                if (j .le. indtm(ip) .or. j .gt. indt) then
                    csum = psido(i, j)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do k = nmolshift + 1, nelorb
                        csum = csum + winv(k - nmolshift, j, ip)*(psiin(ip, 2)&
                            & * barpsiin(k, 1) - psiin(ip, 1)*barpsiin(k, 2)&
                            & + ainvs(ip)*barpsiin(k, 3))
                    end do
                    psido(i, j) = csum
                end if
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif
    else
        do i = 1, nelup
            do j = 1, indt
                if (j .gt. indtm(i)) psiup(i, j) = (0.d0, 0.d0)
            end do
        end do

        do i = 1, neldo
            do j = 1, indt
                ip = i + nelup
                if (j .gt. indtm(ip)) psido(i, j) = (0.d0, 0.d0)
            end do
        end do

!$omp parallel do default(shared) private(i,k,j,csum)
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiup(i, j)
                    do k = 1, nmolipf
                        csum = csum + winv(k, j, i)*(psiin(i, 2)*barpsiin(k, 1)&
                           & - psiin(i, 1)*barpsiin(k, 2) + ainvs(i)*barpsiin(k, 3))
                    end do
                    psiup(i, j) = csum
                end if
            end do
        end do
!$omp end parallel do

!$omp parallel do default(shared) private(i,ip,k,j,csum)
        do i = 1, neldo
            do j = 1, indt + ip4
                ip = i + nelup
                if (j .le. indtm(ip) .or. j .gt. indt) then
                    csum = psido(i, j)
                    do k = nmolshift + 1, nelorb
                        csum = csum + winv(k - nmolshift, j, ip)*(psiin(ip, 2)&
                            & * barpsiin(k, 1) - psiin(ip, 1)*barpsiin(k, 2) + ainvs(ip)*barpsiin(k, 3))
                    end do
                    psido(i, j) = csum
                end if
            end do
        end do
!$omp end parallel do
    end if
    return
end subroutine zgemvtry_pfaff

subroutine zgemvtry(nelorb, nelup, indt, ip4, indtm, winv, psiin, ip_reshuff, psiout)
    use constants, only: zzero, yes_ontarget
    implicit none
    integer nelorb, nelup, indt, indtm(nelup), i, j, k, ip4, ip_reshuff
    complex*16 winv(nelorb, 1:indt + ip4, nelup), psiin(nelorb)&
            &, psiout(indt + ip4, nelup)
    complex*16 csum
    if (yes_ontarget) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
        do i = 1, nelup
            do j = 1, indt + ip4
                psiout(j, i) = (0.d0, 0.d0)
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute parallel do
#endif

#ifdef _OFFLOAD
!$omp target teams distribute collapse(2) private(csum)
#else
!$omp parallel do default(shared) private(i,k,j,csum)
#endif
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiout(j, i)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                    do k = 1, nelorb, ip_reshuff
                        csum = csum + winv(k, j, i)*psiin(k)
                    end do
                    psiout(j, i) = csum
                end if
            end do
        end do
#ifdef _OFFLOAD
!$omp end  target teams distribute
#else
!$omp end parallel do
#endif
    else
        do i = 1, nelup
            do j = 1, indt + ip4
                psiout(j, i) = (0.d0, 0.d0)
            end do
        end do
!$omp parallel do default(shared) private(i,k,j,csum)
        do i = 1, nelup
            do j = 1, indt + ip4
                if (j .le. indtm(i) .or. j .gt. indt) then
                    csum = psiout(j, i)
                    do k = 1, nelorb, ip_reshuff
                        csum = csum + winv(k, j, i)*psiin(k)
                    end do
                    psiout(j, i) = csum
                end if
            end do
        end do
!$omp end parallel do
    end if
    return
end subroutine zgemvtry

! support routines to deal with complex matrices
subroutine makeg_complex(g, psiln, psisn)
    implicit none
    complex*16 g, ph
    real*8 psiln, psisn
    if (abs(g) .ne. 0.d0) then
        ph = log(g)
        psiln = psiln + real(ph)
        psisn = psisn + aimag(ph)
        g = -(1.d0, 0.d0)/g
    end if
    return
end subroutine makeg_complex
subroutine scalepsip(n, g, psip)
#ifdef _OFFLOAD
    use constants, only: yes_ontarget
#endif
    implicit none
    integer n, j
    complex*16 g, psip(n)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do j = 1, n
        psip(j) = psip(j)*g
    end do
    return
end subroutine scalepsip
subroutine up_ainv_down(psip, ainv, nelup, nelcdo)
#ifdef _OFFLOAD
    use constants, only: yes_ontarget
#endif
    implicit none
    integer i
    integer, intent(in) :: nelup, nelcdo
    complex(8), intent(in) :: ainv(nelup, nelup)
    complex(8), intent(out) :: psip(nelup)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do i = 1, nelup
        psip(i) = ainv(nelcdo, i)
    end do
    return
end subroutine
subroutine dcopy_vec_(n, a, b)
    implicit none
    integer n, j
    real*8 a(n), b(n)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do j = 1, n
        b(j) = a(j)
    end do
    return
end subroutine dcopy_vec_
subroutine dcopy_vec(n, a, b)
#ifdef _OFFLOAD
    use constants, only: yes_ontarget
#endif
    implicit none
    integer n, j
    real*8 a(n), b(n)
!    if(n.le.16384) then
!#ifdef _OFFLOAD
!!$omp target update from (a)  if(yes_ontarget)
!#endif
!    do j = 1, n
!        b(j) = a(j)
!    enddo
!#ifdef _OFFLOAD
!!$omp target update to (b)  if(yes_ontarget)
!#endif
!    else
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do j = 1, n
        b(j) = a(j)
    end do
!    endif
    return
end subroutine dcopy_vec
subroutine dcopy_my(n, a, inca, b, incb)
#ifdef _OFFLOAD
    use constants, only: yes_ontarget
#endif
    implicit none
    integer n, j, inca, incb
    real*8 a(inca*(n - 1) + 1), b(incb*(n - 1) + 1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do j = 1, n
        b((j - 1)*incb + 1) = a((j - 1)*inca + 1)
    end do
    return
end subroutine dcopy_my
subroutine dcopy_my_(n, a, inca, b, incb)
    implicit none
    integer n, j, inca, incb
    real*8 a(inca*(n - 1) + 1), b(incb*(n - 1) + 1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do j = 1, n
        b((j - 1)*incb + 1) = a((j - 1)*inca + 1)
    end do
    return
end subroutine dcopy_my_

subroutine zcopy_my(n, a, inca, b, incb)
#ifdef _OFFLOAD
    use constants, only: yes_ontarget
#endif
    implicit none
    integer n, j, inca, incb
    complex*16 a(inca*(n - 1) + 1), b(incb*(n - 1) + 1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do j = 1, n
        b((j - 1)*incb + 1) = a((j - 1)*inca + 1)
    end do
    return
end subroutine zcopy_my
subroutine zcopy_my_(n, a, inca, b, incb)
    implicit none
    integer n, j, inca, incb
    complex*16 a(inca*(n - 1) + 1), b(incb*(n - 1) + 1)
#ifdef _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do j = 1, n
        b((j - 1)*incb + 1) = a((j - 1)*inca + 1)
    end do
    return
end subroutine zcopy_my_
