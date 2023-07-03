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

subroutine scratchdet(winv, winvbar, ainv, psidetln, psip, scratch_getri, ipsip, info)
    use constants, only: ipc, ipf, zone, zzero, nbdgetri, yes_ontarget
    use allio, only: nelorbh, nelorb, nelup, nel, indt4, neldo, nelup_mat, nel_mat, ndiff, npar_eagp, eagp_pfaff, agpn, epscuttype

    implicit none
    integer nelorb5, i, j, k, psidetsn, ipsip(nelup_mat), info
    !  integer,intent(in) :: nelorb,nelup,nel,indt4
    real*8 winv(ipc*nelorb, 0:indt4, nel), winvbar(ipc*ipf*nelorb, nel_mat), &
        ainv(ipc*nelup_mat, nelup_mat), psip(ipc*nelup_mat, nelup_mat), psidetlnt, &
        psidetln, psidetsnc, scratch_getri(ipc*nbdgetri)

    !      This subroutine recompute by scratch the inverse matrix ainv that
    !      can deteriorate with many rank-1 updates.
#ifdef _OFFLOAD
!$omp target update to (winv,winvbar) if(.not.yes_ontarget)
#endif

    nelorb5 = nelorb*(indt4 + 1)

    if (ipf .eq. 2) then
        if (ndiff .ne. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
            do j = 1, nelup_mat
                do i = 1, ipc*nelup_mat
                    psip(i, j) = 0.d0
                end do
            end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
        end if
        if (ipc .eq. 2) then
            call zgemm_('T', 'N', nelup, nelup_mat, nelorbh, zone, winv               &
                    &, nelorb5, winvbar, 2*nelorbh, zzero, psip, nelup_mat)
            call zgemm_('T', 'N', neldo, nelup_mat, nelorbh, zone, winv(1, 0, nelup + 1)&
                    &, nelorb5, winvbar(2*nelorbh + 1, 1), 2*nelorbh, zzero, psip(2*nelup + 1, 1), nelup_mat)

            if (ndiff .ne. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
                do k = 1, ndiff
                    do j = 1, nelup_mat - ndiff
                        psip(2*nelup_mat - 2*k + 1, j) = -psip(2*j - 1, nelup_mat - k + 1)
                        psip(2*nelup_mat - 2*k + 2, j) = -psip(2*j, nelup_mat - k + 1)
                    end do
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif

                if (npar_eagp .gt. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
                    do k = 1, ndiff
                        do j = 1, ndiff
                            psip(2*nelup_mat - 2*ndiff + 2*k - 1, nelup_mat - ndiff + j) = eagp_pfaff(2*k - 1, j)
                            psip(2*nelup_mat - 2*ndiff + 2*k, nelup_mat - ndiff + j) = eagp_pfaff(2*k, j)
                        end do
                    end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
                end if

            end if
!           agpn = psip
            if (epscuttype .gt. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
                do j = 1, nelup_mat
                    do i = 1, ipc*nelup_mat
                        agpn(i, j) = psip(i, j)
                    end do
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
            end if

#ifdef _OFFLOAD
!$omp target update from (psip)
#endif
            call zsktrf('U', 'N', nelup_mat, psip, nelup_mat, ipsip, scratch_getri, nbdgetri, info)
            !    call zgetrf(nelup_mat,nelup_mat,psip,nelup_mat,ipsip,info)
        else
            call dgemm_('T', 'N', nelup, nelup_mat, nelorbh, 1.d0, winv               &
                    &, nelorb5, winvbar, 2*nelorbh, 0.d0, psip, nelup_mat)
            call dgemm_('T', 'N', neldo, nelup_mat, nelorbh, 1.d0, winv(1, 0, nelup + 1)&
                    &, nelorb5, winvbar(nelorbh + 1, 1), 2*nelorbh, 0.d0, psip(nelup + 1, 1), nelup_mat)

            if (ndiff .ne. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do
#endif
                do k = 1, ndiff
                    psip(nelup_mat - k + 1, 1:nelup_mat - ndiff) = -psip(1:nelup_mat - ndiff, nelup_mat - k + 1)
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
                if (npar_eagp .gt. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
                    do k = 1, ndiff
                        do j = 1, ndiff
                            psip(nelup_mat - ndiff + k, nelup_mat - ndiff + j) = eagp_pfaff(k, j)
                        end do
                    end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
                end if
            end if
!           agpn = psip
            if (epscuttype .gt. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
                do j = 1, nelup_mat
                    do i = 1, ipc*nelup_mat
                        agpn(i, j) = psip(i, j)
                    end do
                end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
            end if
#ifdef _OFFLOAD
!$omp target update from (psip)
#endif
            call dsktrf('U', 'N', nelup_mat, psip, nelup_mat, ipsip, scratch_getri, nbdgetri, info)

            !    call dgetrf(nelup_mat,nelup_mat,psip,nelup_mat,ipsip,info)
        end if

    else

        if (ipc .eq. 2) then

            call zgemm_('T', 'N', nelup, nelup, nelorbh, zone, winv, nelorb5, &
                        winvbar(1, nelup + 1), nelorbh, zzero, psip, nelup)

!           agpn = psip
            if (epscuttype .gt. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
                do j = 1, nelup_mat
                    do i = 1, ipc*nelup_mat
                        agpn(i, j) = psip(i, j)
                    end do
                end do
            end if
#ifdef _OFFLOAD
#ifndef _CUSOLVER
!$omp target update from (psip)
#endif
#endif
            call zgetrf_(nelup, nelup, psip, nelup, ipsip, info)

        else

            call dgemm_('T', 'N', nelup, nelup, nelorbh, 1.d0, winv, nelorb5, &
                        winvbar(1, nelup + 1), nelorbh, 0.d0, psip, nelup)

!           agpn = psip
            if (epscuttype .gt. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2)
#endif
                do j = 1, nelup_mat
                    do i = 1, ipc*nelup_mat
                        agpn(i, j) = psip(i, j)
                    end do
                end do
            end if
#ifdef _OFFLOAD
#ifndef _CUSOLVER
!$omp target update from (psip)
#endif
#endif
            call dgetrf_(nelup, nelup, psip, nelup, ipsip, info)
        end if

    end if

    if (info .ne. 0 .and. ipf .eq. 1) then

        write (6, *) 'ERROR IN Z/DGETRF info  =', info
        call dscalzero(ipc*nelup_mat*nelup_mat, 0.d0, psip, 1)
#ifdef _OFFLOAD
!$omp target update from  (winv,winvbar) if(yes_ontarget)
#endif
        if (ipf .eq. 2) then
            if (ndiff .ne. 0) psip(:, 1:nelup_mat) = 0.d0
            if (ipc .eq. 2) then
                call zgemm('T', 'N', nelup, nelup_mat, nelorbh, zone, winv               &
                        &, nelorb5, winvbar, 2*nelorbh, zzero, psip, nelup_mat)
                call zgemm('T', 'N', neldo, nelup_mat, nelorbh, zone, winv(1, 0, nelup + 1)&
                        &, nelorb5, winvbar(2*nelorbh + 1, 1), 2*nelorbh, zzero, psip(2*nelup + 1, 1), nelup_mat)
                if (ndiff .ne. 0) then
                    do k = 1, ndiff
                        do j = 1, nelup_mat - ndiff
                            psip(2*nelup_mat - 2*k + 1, j) = -psip(2*j - 1, nelup_mat - k + 1)
                            psip(2*nelup_mat - 2*k + 2, j) = -psip(2*j, nelup_mat - k + 1)
                        end do
                    end do
                    if (npar_eagp .gt. 0) then
                        do k = 1, ndiff
                            do j = 1, ndiff
                                psip(2*nelup_mat - 2*ndiff + 2*k - 1, nelup_mat - ndiff + j) = eagp_pfaff(2*k - 1, j)
                                psip(2*nelup_mat - 2*ndiff + 2*k, nelup_mat - ndiff + j) = eagp_pfaff(2*k, j)
                            end do
                        end do
                    end if

                end if
            else
                call dgemm('T', 'N', nelup, nelup_mat, nelorbh, 1.d0, winv               &
                        &, nelorb5, winvbar, 2*nelorbh, 0.d0, psip, nelup_mat)
                call dgemm('T', 'N', neldo, nelup_mat, nelorbh, 1.d0, winv(1, 0, nelup + 1)&
                        &, nelorb5, winvbar(nelorbh + 1, 1), 2*nelorbh, 0.d0, psip(nelup + 1, 1), nelup_mat)
                if (ndiff .ne. 0) then
                    do k = 1, ndiff
                        psip(nelup_mat - k + 1, 1:nelup_mat - ndiff) = -psip(1:nelup_mat - ndiff, nelup_mat - k + 1)
                    end do
                    if (npar_eagp .gt. 0) then
                        do k = 1, ndiff
                            do j = 1, ndiff
                                psip(nelup_mat - ndiff + k, nelup_mat - ndiff + j) = eagp_pfaff(k, j)
                            end do
                        end do
                    end if
                end if
            end if
        else
            if (ipc .eq. 2) then
                call zgemm('T', 'N', nelup, nelup, nelorbh, zone, winv               &
                        &, nelorb5, winvbar(1, nelup + 1), nelorbh, zzero, psip, nelup)
            else
                call dgemm('T', 'N', nelup, nelup, nelorbh, 1.d0, winv               &
                        &, nelorb5, winvbar(1, nelup + 1), nelorbh, 0.d0, psip, nelup)
            end if
        end if
        if (ipc .eq. 2) then
            write (6, *) ' input  matrix element DET AGP (real/imag) '
            do i = 1, nelup_mat
                do j = 1, nelup_mat
                    write (6, *) i, j, psip(2*i - 1, j), psip(2*i, j)
                end do
            end do
        else
            write (6, *) ' input  matrix element DET AGP  '
            do i = 1, nelup_mat
                do j = 1, nelup_mat
                    write (6, *) i, j, psip(i, j)
                end do
            end do
        end if
        return

    else

        if (ipc .eq. 1) then
            call evaldet(psip, nelup_mat, nelup_mat, ipsip, psidetlnt, psidetsn)
            if (ipf .eq. 2) then
! HERE psip is already in the CPU after dsktrf
                call dsktri('U', nelup_mat, psip, nelup_mat, ainv, nelup_mat, ipsip, scratch_getri, info)
            else
! HERE Same as above
                call dgetri_(nelup_mat, psip, nelup_mat, ipsip, scratch_getri, nbdgetri, info)
            end if
        else
            call evaldet_complex(psip, nelup_mat, nelup_mat, ipsip, psidetlnt, psidetsnc)
            if (ipf .eq. 2) then
! HERE Same as above
                call zsktri('U', nelup_mat, psip, nelup_mat, ainv, nelup_mat, ipsip, scratch_getri, info)
            else
!HERE as  above
                call zgetri_(nelup_mat, psip, nelup_mat, ipsip, scratch_getri, nbdgetri, info)
            end if
        end if

        if (info .ne. 0) then
            write (6, *) 'ERROR IN Z/DGETRI info =', info
            return
        end if
        if (ipf .ne. 2) then
! Offload  HERE only if _CUSOLVER psip is in the GPU
#ifdef _CUSOLVER
!$omp target teams distribute parallel do collapse(2)
#endif
            do i = 1, nelup_mat
                do j = 1, ipc*nelup_mat
                    ainv(j, i) = psip(j, i)
                end do
            end do
        end if
        ! antisymmetrization ainv
        if (ipf .eq. 2) then
            if (ipc .eq. 1) then
!HERE ainv is in the CPU
                do i = 1, nelup_mat
                    do j = 1, nelup_mat
                        if (j .gt. i) then
                            ainv(i, j) = 0.5d0*(ainv(i, j) - ainv(j, i))
                            ainv(j, i) = -ainv(i, j)
                        end if
                    end do
                    ainv(i, i) = 0.d0
                end do
            else
!HERE ainv is in the CPU
                do i = 1, nelup_mat
                    do j = 1, nelup_mat
                        if (j .gt. i) then
                            ainv(2*i - 1, j) = 0.5d0*(ainv(2*i - 1, j) - ainv(2*j - 1, i))
                            ainv(2*i, j) = 0.5d0*(ainv(2*i, j) - ainv(2*j, i))
                            ainv(2*j - 1, i) = -ainv(2*i - 1, j)
                            ainv(2*j, i) = -ainv(2*i, j)
                        end if
                    end do
                    ainv(2*i - 1, i) = 0.d0
                    ainv(2*i, i) = 0.d0
                end do
            end if
        end if
    end if
    psidetln = psidetlnt
!   OFFLOAD array distribution
!  agpn is always in the target regardless yes_ontarget
!  ainv is in the target only if _CUSOLVER and ipf==1,  otherwise in  the  cpu
! HERE only if _CUSOLVER ainv is in the GPU
#ifdef _OFFLOAD
!$omp target update from (agpn) if(.not.yes_ontarget)
#ifdef _CUSOLVER
!$omp target update from (ainv) if(.not.yes_ontarget.and.ipf.eq.1)
#else
!$omp target update to (ainv) if(yes_ontarget)
#endif
#endif
    return
end subroutine scratchdet
