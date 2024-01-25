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

subroutine upnewwf(indt, i0, indtm, typecomp, nshell, ioptorb, iocc, kel, nel, r, rmu, dd, zeta, rion, distp, z, nelskip, &
                   nion, kion, iflagnorm, cnorm, LBox, rmucos, rmusin, mindist, indpar_tab, indorb_tab, indshell_tab, yesupel)

    use allio, only: ikshift, iespbc, rank, gamma_point, yes_crystalj&
        &, yes_scemama, lepsbas, novec_loop1, slaterorb_read, nshell_det&
        &, use_qmckl, qmckl_ctx
    use Cell, only: cellscale, cellpi, rphase, phase2pi, phase2pi_down, sinphase, cosphase, s2r, car2cry
    use Constants, only: ipc
    use qmckl
#ifdef _QMCKL_GPU
    !use qmckl_gpu
#endif

    implicit none
    ! input
    integer nel, mshift, indorb, indpar, indshell, dimp &
        , j, i, i_ion, indt, nelskip, nshell, indtmin, nion, k, kk &
        , iflagnorm, indtm, i0, typecomp, ii, jj, ll, ierr, ishift, indp1, indp2, indp3 &
        , indp4, indp5, dimx, dimy, i0u
    real*8 rmu(3, 0:indtm, nion), rion(3, nion), distp((indtm + 1)*27*max(nshell, nion) +&
            &nelskip*(indt + 4 - i0 + 1)), cnorm(nshell), zeta(nion), r(0:indtm, nion), rcos_sav(3)&
            &, mindist, rsign(3)
    ! specify the phase choice for crystal basis, ignored for Jastrow
    logical, intent(in) :: yesupel
    logical do_makefun, yeszero_z(0:indtm)
!   logical, external:: slaterorb
    ! real(8), dimension(:,:,:), allocatable::  distp_true

    ! local variables
    real*8 kpip(3)
    real*8 kel(3, nel, 0:indtm)
    real*8 LBox, Lbox_sav
    ! in case abs(LBox).eq.3.d0, rmucos and rmusin are used as npip and rmunew
    real*8 rmucos(3, 0:indtm, nion), rmusin(3, 0:indtm, nion)
    ! real(8), dimension(:,:), allocatable :: psip_r
    logical iesjas
    integer :: ioptorb(nshell), kion(nshell)&
            &, indpar_tab(nshell + 1), indorb_tab(nshell + 1), indshell_tab(nshell + 1)
    real*8, intent(out) :: z(min(ipc, indpar_tab(1) + 2)*nelskip, i0:&
            &(1 - typecomp)*(indt + 4) + typecomp*indtm)
    integer :: iocc(indshell_tab(nshell + 1))
    real*8 dd(indpar_tab(nshell + 1))
    real*8 phs(3) ! scratch phase for Lbox=3
    !
    ! QMCKL
    integer*4                      :: rc
    integer*8, save                :: ao_num=0, npoints_qmckl=0
    double precision, allocatable  :: ao_vgl_qmckl(:,:), ao_value_qmckl(:,:)
    !
    integer, external :: omp_get_num_threads
    ! ---------------------------------------------------------------------
    ! orbital types are defined by the variable Lbox:
    ! Lbox = -3.d0 --> Complex open
    ! Lbox = -2.d0 --> open boundary conditions with bump orbitals
    ! Lbox = -1.d0 --> open boundary conditions
    ! Lbox =  1.d0 --> PBC/APBC old version (real w.f.)
    ! Lbox =  2.d0 --> PBC/APBC with bump orbitals
    ! Lbox =  3.d0 --> Crystal basis set (complex w.f. for PBC)
    ! ---------------------------------------------------------------------
    ! initialization

    if (indtm .eq. 0 .and. i0 .eq. 0 .and. novec_loop1) then
        call upnewwf0(indt, typecomp, nshell, ioptorb, iocc, kel, nel, r, rmu, dd, zeta, rion, distp, z, nelskip, &
                      nion, kion, iflagnorm, cnorm, LBox, rmucos, rmusin, mindist, indpar_tab, indorb_tab, indshell_tab, yesupel)
    else

        indshell = 0
        indorb = 0
        indpar = 0
        indtmin = 0
        ! rphs = 0.d0
        ! allocate(distp_true(0:indtm,20,nshell))
        if (abs(LBox) .eq. 3.d0) then
            indp5 = nion*(indtm + 1)*27 + 1
            indp1 = 2*nion*(indtm + 1) + 1
            indp2 = indp1 + nion*(indtm + 1)
            indp3 = indp2 + nion*(indtm + 1)
            indp4 = indp3 + 3*nion*(indtm + 1)
            dimx = min(ipc, indpar_tab(1) + 2)*nelskip
            dimy = (1 - typecomp)*(indt + 4) + typecomp*indtm
            z = 0.d0
        elseif (yes_scemama) then
            z = 0.d0
        end if

        if (indpar_tab(1) .eq. -1) then
            iesjas = .true.
            !    indpar_tab(1)=0
            mshift = nshell_det
        else
            iesjas = .false.
            mshift = 0
        end if

        LBox_sav = LBox
        ! Crystal basis set for periodic systems
        phs(:) = 0.d0

        dimp = (indtm + 1)*20

        if (LBox .eq. 3.d0) then
            !
            ! initialize different phase for spin up/spin down electrons
            !
            if (yesupel) then
                phs(:) = phase2pi(:)
            else
                phs(:) = phase2pi_down(:)
            end if
        elseif (LBox .eq. -3.d0) then
            phs = 0.d0
        end if
        ishift = 0
        if (abs(LBox) .eq. 3.d0 .and. iesjas) then
            !
            ! always use the old periodic basis for the Jastrow
            !
            if (.not. yes_crystalj) then
                if (iespbc) then
                    !       iflagnorm=2 ! always recompute distances if Jastrow is considered with Lbox.eq.3!!!
                    ! REMEMBER: with Jastrow it pass from Lbox.eq.3 to Lbox.eq.1, distances
                    ! to be computed are different!!!
                    LBox = 1.d0 ! always use the old periodic basis for the Jastrow
                else
                    Lbox = -1.d0
                end if
            else
                ishift = ikshift
            end if
            phs = 0.d0
            ! elseif(abs(Lbox).eq.3.d0) then
            ! allocate(psip_r(nelskip,i0:indt+4))
        end if

        ! electron-ion distances

        ! if((iflagnorm.gt.0.or.(LBox.eq.2.and..not.gamma_point)).and.abs(LBox).ne.3) then
        if (abs(LBox) .ne. 3) then
            if (typecomp .eq. 1) then
                i0u = i0
            else
                i0u = 0
            end if

            do k = 1, nion
                do i = indtmin, indtm
                    rmu(1, i, k) = kel(1, 1, i) - rion(1, k)
                    rmu(2, i, k) = kel(2, 1, i) - rion(2, k)
                    rmu(3, i, k) = kel(3, 1, i) - rion(3, k)
                end do
            end do

            ! ****** Periodic Systems *************
            if (LBox .eq. 1.d0) then

                if (.not. gamma_point .and. .not. iesjas) then

                    do j = indtmin, indtm
                        do k = 1, nion
                            rcos_sav(1) = dcos(rphase(1)*rmu(1, j, k))
                            rcos_sav(2) = dcos(rphase(2)*rmu(2, j, k))
                            rcos_sav(3) = dcos(rphase(3)*rmu(3, j, k))
                            cosphase(j, k) = rcos_sav(1)*rcos_sav(2)*rcos_sav(3)
                            sinphase(1, j, k) = dsin(rphase(1)*rmu(1, j, k))*rcos_sav(2)*rcos_sav(3)
                            sinphase(2, j, k) = dsin(rphase(2)*rmu(2, j, k))*rcos_sav(1)*rcos_sav(3)
                            sinphase(3, j, k) = dsin(rphase(3)*rmu(3, j, k))*rcos_sav(1)*rcos_sav(2)
                        end do
                    end do

                end if

                do j = indtmin, indtm
                    do k = 1, nion
                        rmu(1, j, k) = rmu(1, j, k)/cellpi(1)
                        rmu(2, j, k) = rmu(2, j, k)/cellpi(2)
                        rmu(3, j, k) = rmu(3, j, k)/cellpi(3)
                        rmucos(1, j, k) = dcos(rmu(1, j, k))
                        rmucos(2, j, k) = dcos(rmu(2, j, k))
                        rmucos(3, j, k) = dcos(rmu(3, j, k))
                        rmusin(1, j, k) = dsin(rmu(1, j, k))
                        rmusin(2, j, k) = dsin(rmu(2, j, k))
                        rmusin(3, j, k) = dsin(rmu(3, j, k))
                        rmu(1, j, k) = cellpi(1)*rmusin(1, j, k)
                        rmu(2, j, k) = cellpi(2)*rmusin(2, j, k)
                        rmu(3, j, k) = cellpi(3)*rmusin(3, j, k)
                        r(j, k) = max(dsqrt(rmu(1, j, k)**2 + rmu(2, j, k)**2 + rmu(3, j, k)**2), mindist)
                    end do
                end do

            elseif (LBox .eq. 2.d0) then

                if (iesjas) then
                    do i = indtmin, indtm
                        do k = 1, nion
                            kpip(1) = anint(rmu(1, i, k)/cellscale(1))
                            kpip(2) = anint(rmu(2, i, k)/cellscale(2))
                            kpip(3) = anint(rmu(3, i, k)/cellscale(3))
                            rmu(1, i, k) = rmu(1, i, k) - kpip(1)*cellscale(1)
                            rmu(2, i, k) = rmu(2, i, k) - kpip(2)*cellscale(2)
                            rmu(3, i, k) = rmu(3, i, k) - kpip(3)*cellscale(3)
                            r(i, k) = dsqrt(rmu(1, i, k)**2 + rmu(2, i, k)**2 + rmu(3, i, k)**2)
                        end do
                    end do
                else
                    do i = indtmin, indtm
                        do k = 1, nion
                            kpip(1) = anint(rmu(1, i, k)/cellscale(1))
                            kpip(2) = anint(rmu(2, i, k)/cellscale(2))
                            kpip(3) = anint(rmu(3, i, k)/cellscale(3))
                            rmu(1, i, k) = rmu(1, i, k) - kpip(1)*cellscale(1)
                            rmu(2, i, k) = rmu(2, i, k) - kpip(2)*cellscale(2)
                            rmu(3, i, k) = rmu(3, i, k) - kpip(3)*cellscale(3)
                            r(i, k) = dsqrt(rmu(1, i, k)**2 + rmu(2, i, k)**2 + rmu(3, i, k)**2)
                            do kk = 1, 3
                                if (phase2pi(kk) .ne. 0.d0 .and. 2*(kpip(kk)/2) .ne. kpip(kk)) r(i, k) = -r(i, k)
                            end do
                        end do
                    end do

                end if

            elseif (LBox .eq. -2.d0) then
                ! No need of the cutoff mindist
                do i = indtmin, indtm
                    do k = 1, nion
                        r(i, k) = dsqrt(rmu(1, i, k)**2 + rmu(2, i, k)**2 + rmu(3, i, k)**2)
                    end do
                end do

            else

                do k = 1, nion
                    do i = indtmin, indtm
                        r(i, k) = max(dsqrt(rmu(1, i, k)**2 + rmu(2, i, k)**2 + rmu(3, i, k)**2), mindist)
                    end do
                end do

            end if

        end if

        ! if(abs(LBox).eq.3.d0.and.iflagnorm.gt.0) then
        if (abs(LBox) .eq. 3.d0) then
            !if(iflagnorm.gt.0) then
            if (iespbc) then
!$omp parallel do default(shared) private(i,k)
                do k = 1, nion
                    do i = indtmin, indtm
                        rmu(:, i, k) = kel(:, 1, i) - rion(:, k)
                        rmucos(:, i, k) = &
                                &car2cry(:, 1)*rmu(1, i, k) + car2cry(:, 2)*rmu(2, i, k) + car2cry(:, 3)*rmu(3, i, k)
                        rmucos(1, i, k) = anint(rmucos(1, i, k)/cellscale(1))
                        rmucos(2, i, k) = anint(rmucos(2, i, k)/cellscale(2))
                        rmucos(3, i, k) = anint(rmucos(3, i, k)/cellscale(3))
                        rmu(:, i, k) = rmu(:, i, k)&
                                & - s2r(:, 1)*rmucos(1, i, k) - s2r(:, 2)*rmucos(2, i, k) - s2r(:, 3)*rmucos(3, i, k)
                        r(i, k) = max(dsqrt(rmu(1, i, k)**2 + rmu(2, i, k)**2 + &
                                            rmu(3, i, k)**2), mindist)
                    end do
                end do
!$omp end parallel do

            else
!$omp parallel do default(shared) private(i,k)
                do k = 1, nion
                    do i = indtmin, indtm
                        rmu(1, i, k) = kel(1, 1, i) - rion(1, k)
                        rmu(2, i, k) = kel(2, 1, i) - rion(2, k)
                        rmu(3, i, k) = kel(3, 1, i) - rion(3, k)
                        r(i, k) = max(dsqrt(rmu(1, i, k)**2 + rmu(2, i, k)**2 + &
                                            rmu(3, i, k)**2), mindist)
                    end do
                end do
!$omp end parallel do
            end if
        end if

        !
        if (LBox .eq. 1.d0) then
!$omp parallel do default(shared)&
!$omp private(i,indpar,indorb,indshell)
            do i = 1, nshell
                indpar = max(indpar_tab(i), 0)
                indorb = indorb_tab(i)
                indshell = indshell_tab(i)
                call makefun_pbc(ioptorb(i), iocc, indt, i0, indtmin, indtm, typecomp &
                            &, indpar, indorb, indshell, nelskip, z, dd&
                            &, r(0, kion(i)), rmu(1, 0, kion(i)), distp(dimp*(i - 1) + 1) &
                            &, iflagnorm, cnorm(i), rmucos(1, 0, kion(i))&
                            &, rmusin(1, 0, kion(i)), sinphase(1, 0, kion(i))&
                            &, cosphase(0, kion(i)))
            end do
!$omp end parallel do
        elseif (abs(LBox) .eq. 2.d0) then
!$omp parallel do default(shared)&
!$omp  private(i,indorb,indshell,indpar)
            do i = 1, nshell
                indpar = max(indpar_tab(i), 0)
                indorb = indorb_tab(i)
                indshell = indshell_tab(i)
                call makefun_bump(ioptorb(i), iocc, indt, i0, indtmin, indtm, typecomp &
                        , indpar, indorb, indshell, nelskip, z, dd, r(0, kion(i)), rmu(1, 0, kion(i))&
                        &, distp(dimp*(i - 1) + 1), iflagnorm, cnorm(i))
            end do
!$omp end parallel do

        elseif (abs(LBox) .eq. 3.d0) then
            call makefun_grid(distp(indp4), distp(indp5), distp(indp3), distp(indp2)&
                    &, distp(indp1), distp&
                    &, nshell, nion, i0, indtmin, indtm, indt, typecomp, ishift, iflagnorm&
                    &, kion, ioptorb, indpar_tab, indorb_tab, indshell_tab, mindist, phs, dd, z, dimx, dimy&
                    &, nelskip, iesjas, rmu, rion, cnorm, zeta, rmucos)
        else
            if (.not.use_qmckl) then
!$omp parallel do default(shared)&
!$omp private(i,ll,i_ion,indpar,indorb,indshell,do_makefun,yeszero_z)
                do i = 1, nshell
                    indpar = max(indpar_tab(i), 0)
                    indorb = indorb_tab(i)
                    indshell = indshell_tab(i)
                    do_makefun = .true.
                    if (yes_scemama .and. ioptorb(i) .ne. 200) then
                        i_ion = kion(i)
                        if (slaterorb_read(i + mshift)) then
                            do ll = i0u, indtm
                                if (dd(indpar + 1)*r(ll, i_ion) .lt. lepsbas) then
                                    yeszero_z(ll) = .false.
                                else
                                    yeszero_z(ll) = .true.
                                end if
                            end do
                        else
                            do ll = i0u, indtm
                                if (dd(indpar + 1)*r(ll, i_ion)*r(ll, i_ion) .lt. lepsbas) then
                                    yeszero_z(ll) = .false.
                                else
                                    yeszero_z(ll) = .true.
                                end if
                            end do
                        end if
                        do_makefun = .not. all(yeszero_z(i0u:indtm))
                    end if
                    if (do_makefun) then
                        call makefun(ioptorb(i), indt, i0, indtmin, indtm, typecomp&
                               &, indpar, indorb, indshell, nelskip, z, dd, zeta(kion(i)), r(0, kion(i)), rmu(1, 0, kion(i))&
                               &, distp(dimp*(i - 1) + 1), iflagnorm, cnorm(i))
                        if (yes_scemama .and. indtm .gt. 0 .and. ioptorb(i) .ne. 200) then
                            do ll = i0, indtm
                                if (yeszero_z(ll)) then
                                    z(indorb_tab(i) + 1:indorb_tab(i + 1), ll) = 0.d0
                                end if
                            end do
                            if (yeszero_z(0) .and. typecomp .ne. 1) then
                                do ll = indt + 1, indt + 4
                                    z(indorb_tab(i) + 1:indorb_tab(i + 1), ll) = 0.d0
                                end do
                            end if
                        end if
                    end if
                end do
!$omp end parallel do
#ifdef _QMCKL
            else
                if (npoints_qmckl == 0) then
                    rc = qmckl_get_ao_basis_ao_num(qmckl_ctx, ao_num)
                    if (rc /= QMCKL_SUCCESS) then
                        write(0,*) 'Error getting ao_num'
                        call abort()
                    end if
                end if

                if (typecomp.eq.1) then   ! Only values
                    npoints_qmckl = (indtm-i0+1)*1_8
                    allocate(ao_value_qmckl(ao_num, i0:indtm))
                    rc = qmckl_set_point(qmckl_ctx, 'N', npoints_qmckl, kel(1:3,1,i0:indtm), 3_8*npoints_qmckl)

                    if (rc /= QMCKL_SUCCESS) then
                        write(0,*) 'Error setting electron coords'
                        call abort()
                    end if

                    rc = qmckl_get_ao_basis_ao_value_inplace(        &
                            qmckl_ctx,                               &
                            ao_value_qmckl,                          &
                            ao_num*npoints_qmckl)

                    if (rc /= QMCKL_SUCCESS) then
                        write(0,*) 'Error getting AOs from QMCkl'
                        call abort()
                    end if

                    do jj=i0,indtm
                        do ii=1,ao_num
                            z(ii,jj) = ao_value_qmckl(ii,jj)
                        end do
                    end do

                    deallocate(ao_value_qmckl)
                else
                    npoints_qmckl = (indtm-i0)*1_8
                    allocate(ao_vgl_qmckl(ao_num, 5))
                    allocate(ao_value_qmckl(ao_num, i0+1:indtm))
                    rc = qmckl_set_point(qmckl_ctx, 'N', npoints_qmckl, kel(1:3,1,i0+1:indtm), 3_8*npoints_qmckl)

                    if (rc /= QMCKL_SUCCESS) then
                        write(0,*) 'Error setting electron coords'
                        call abort()
                    end if

                    rc = qmckl_get_ao_basis_ao_value_inplace(        &
                            qmckl_ctx,                               &
                            ao_value_qmckl,                          &
                            ao_num*npoints_qmckl)

                    if (rc /= QMCKL_SUCCESS) then
                        write(0,*) 'Error getting AOs from QMCkl'
                        call abort()
                    end if

                    rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, kel(1:3,1,i0), 3_8)

                    if (rc /= QMCKL_SUCCESS) then
                        write(0,*) 'Error setting electron coords 2'
                        call abort()
                    end if

                    rc = qmckl_get_ao_basis_ao_vgl_inplace(          &
                            qmckl_ctx,                               &
                            ao_vgl_qmckl,                            &
                            ao_num*5_8)

                    if (rc /= QMCKL_SUCCESS) then
                        write(0,*) 'Error getting AOs from QMCkl 2'
                        call abort()
                    end if

                    do jj=i0+1,indtm
                        do ii=1,ao_num
                            z(ii,jj) = ao_value_qmckl(ii,jj)
                        end do
                    end do

                    do ii=1,ao_num
                        z(ii,i0)     = ao_vgl_qmckl(ii,1)
                        z(ii,indt+1) = ao_vgl_qmckl(ii,2)
                        z(ii,indt+2) = ao_vgl_qmckl(ii,3)
                        z(ii,indt+3) = ao_vgl_qmckl(ii,4)
                        z(ii,indt+4) = ao_vgl_qmckl(ii,5)
                    end do

                    deallocate(ao_value_qmckl)
                    deallocate(ao_vgl_qmckl)
                end if
#endif
            end if
        end if

        if (iesjas) then
            LBox = LBox_sav
        end if
    end if
    return
end subroutine upnewwf
subroutine makefun_grid(distp_true_, distp_, rmusin, r, rphs, cphs_r&
        &, nshell, nion, i0, indtmin, indtm, indt, typecomp, ishift, iflagnorm&
        &, kion, ioptorb, indpar_tab, indorb_tab, indshell_tab, mindist, phs, dd, z, dimx, dimy&
        &, nelskip, iesjas, rmu, rion, cnorm, zeta, rmucos)
!   use allio, only : allgrid, point_shell
    use allio, only: kgrid, kgrid_atom, adr_nion, adrj_nion, ind_nion, indj_nion, lepsbas, yes_scemama, slaterorb_read, nshell_det
    use Cell, only: s2r
    use Constants, only: ipc, zimg
    implicit none
    integer, intent(in) :: nelskip, i0, nshell, indtm, indt, nion, indtmin, ishift&
            &, dimx, dimy, typecomp
    real*8, intent(in) :: mindist, phs(3), rion(3, nion), cnorm(nshell), zeta(nion)&
            &, rmu(3, 0:indtm, nion)
    real*8 distp_true_(0:indtm, 20, nion), distp_(nelskip, i0:*)
    integer i, j, ii, jj, kk, ll, indpar, indorb, indshell
    integer, intent(in) :: kion(nshell), ioptorb(nshell), indpar_tab(nshell + 1)&
            &, indorb_tab(nshell + 1), indshell_tab(nshell + 1)
    real*8, intent(in) :: dd(indpar_tab(nshell + 1)), rmucos(3, 0:indtm, nion)
    real*8 kpip(3)
    complex*16 cphs_r(0:indtm, nion)
    real*8 rphs(0:indtm, nion)
    real*8 rmusin(3, 0:indtm, nion), r(0:indtm, nion)
    real*8, intent(out) :: z(dimx, i0:dimy)
    integer iflagnorm, indfirst, indlast, i_ion, j_ion, nshift, mshift, i0u, case_if, case_upz
    real*8 cphs
    logical iesjas, gammaorj, do_makefun
    logical yeszero_z(0:indtm)
    if (iesjas) then
        nshift = nion
        mshift = nshell_det
    else
        nshift = 0
        mshift = 0
    end if
    if (typecomp .eq. 1) then
        i0u = i0
    else
        i0u = 0
    end if

    if (sum(abs(phs(:))) .eq. 0.d0) then
        gammaorj = .true.
    else
        gammaorj = .false.
    end if
    if (ipc .eq. 2 .and. .not. gammaorj) then ! complex case
        case_if = 1
    elseif (.not. gammaorj .and. ipc .eq. 1) then
        case_if = 2
    else if (ipc .eq. 2 .and. .not. iesjas) then
        case_if = 3
    else
        case_if = 4
    end if
    if (ipc .eq. 2 .and. .not. iesjas) then ! complex case
        if (typecomp .eq. 1) then
            case_upz = 1
        else
            case_upz = 2
        end if
    else ! real case
        if (typecomp .eq. 1) then
            case_upz = 3
        else
            case_upz = 4
        end if
    end if

!$omp parallel do default(shared)&
!$omp   private(i_ion,j_ion,indfirst,indlast,i,j,ii,jj,kk,ll,indpar,indorb,indshell,kpip,cphs,do_makefun,yeszero_z)
    do i_ion = 1, nion
        if (iesjas) then
            indfirst = adrj_nion(i_ion)
            indlast = adrj_nion(i_ion + 1) - 1
        else
            indfirst = adr_nion(i_ion)
            indlast = adr_nion(i_ion + 1) - 1
        end if
        do ii = 1, kgrid_atom(i_ion + nshift)%dimshell
            kpip(1) = kgrid_atom(i_ion + nshift)%kpip(1, ii)
            kpip(2) = kgrid_atom(i_ion + nshift)%kpip(2, ii)
            kpip(3) = kgrid_atom(i_ion + nshift)%kpip(3, ii)
! Vectorized version, almost optimal
            do j = indtmin, indtm
                rmusin(1, j, i_ion) = rmu(1, j, i_ion) + &
            &s2r(1, 1)*kpip(1) + s2r(1, 2)*kpip(2) + s2r(1, 3)*kpip(3)
                rmusin(2, j, i_ion) = rmu(2, j, i_ion) + &
            &s2r(2, 1)*kpip(1) + s2r(2, 2)*kpip(2) + s2r(2, 3)*kpip(3)
                rmusin(3, j, i_ion) = rmu(3, j, i_ion) + &
            &s2r(3, 1)*kpip(1) + s2r(3, 2)*kpip(2) + s2r(3, 3)*kpip(3)
                r(j, i_ion) = max(dsqrt(rmusin(1, j, i_ion)**2 + rmusin(2, j, i_ion)**2 + rmusin(3, j, i_ion)**2), mindist)
            end do
! Slower version, too much  memory  access.
!            do j = indtmin, indtm
!            rmusin(:, j, i_ion) = rmu(:, j,i_ion) + &
!          &s2r(:, 1) * kpip(1) + s2r(:, 2) * kpip(2) + s2r(:, 3) * kpip(3)
!            r(j,i_ion)=rmusin(1,j,i_ion)**2+rmusin(2,j,i_ion)**2+rmusin(3,j,i_ion)**2
!            enddo
!            do j=indtmin,indtm
!            r(j,i_ion)=max(dsqrt(r(j,i_ion)),mindist)
!            enddo
            select case (case_if)

            case (1)
                do ll = i0u, indtm
                    cphs = -(phs(1)*(kpip(1) - rmucos(1, ll, i_ion)) + &
                             phs(2)*(kpip(2) - rmucos(2, ll, i_ion)) + &
                             phs(3)*(kpip(3) - rmucos(3, ll, i_ion)))
                    cphs_r(ll, i_ion) = dcmplx(dcos(cphs), dsin(cphs))
                end do
            case (2)
                do ll = i0u, indtm
                    rphs(ll, i_ion) = dcos(phs(1)*(kpip(1) - rmucos(1, ll, i_ion)) + &
                                           phs(2)*(kpip(2) - rmucos(2, ll, i_ion)) + &
                                           phs(3)*(kpip(3) - rmucos(3, ll, i_ion)))
                end do
            case (3)
                cphs_r(i0u:indtm, i_ion) = dcmplx(1.d0, 0.d0)
            case (4)
                rphs(i0u:indtm, i_ion) = 1.d0
            end select
            do j_ion = indfirst, indlast
                if (iesjas) then
                    i = indj_nion(j_ion)
                else
                    i = ind_nion(j_ion)
                end if

                if (kgrid(i + ishift)%tobedone(ii)) then

                    indpar = max(indpar_tab(i), 0)
                    indorb = indorb_tab(i)
                    indshell = indshell_tab(i)
                    do_makefun = .true.
                    if (yes_scemama .and. ioptorb(i) .ne. 200) then
                        if (slaterorb_read(i + mshift)) then
                            do ll = i0u, indtm
                                if (dd(indpar + 1)*r(ll, i_ion) .lt. lepsbas) then
                                    yeszero_z(ll) = .false.
                                else
                                    yeszero_z(ll) = .true.
                                end if
                            end do
                        else
                            do ll = i0u, indtm
                                if (dd(indpar + 1)*r(ll, i_ion)*r(ll, i_ion) .lt. lepsbas) then
                                    yeszero_z(ll) = .false.
                                else
                                    yeszero_z(ll) = .true.
                                end if
                            end do
                        end if
                        do_makefun = .not. all(yeszero_z(i0u:indtm))
                    end if
                    if (do_makefun) then
                        call makefun(ioptorb(i), indt, i0, indtmin, indtm, typecomp&
                                &, indpar, indorb, indshell, nelskip, distp_, dd, zeta(i_ion), r(0, i_ion), rmusin(1, 0, i_ion)&
                                &, distp_true_(0, 1, i_ion), iflagnorm, cnorm(i))
                        if (yes_scemama .and. indtm .gt. 0 .and. ioptorb(i) .ne. 200) then
                            do ll = i0, indtm
                                if (yeszero_z(ll)) then
                                    distp_(indorb_tab(i) + 1:indorb_tab(i + 1), ll) = 0.d0
                                end if
                            end do
                            if (yeszero_z(0) .and. typecomp .ne. 1) then
                            do ll = indt + 1, indt + 4
                                distp_(indorb_tab(i) + 1:indorb_tab(i + 1), ll) = 0.d0
                            end do
                            end if
                        end if

                        select case (case_upz)
                        case (1)
                            do ll = i0, indtm
                                call zaxrpy(indorb_tab(i), indorb, cphs_r(ll, i_ion), distp_(1, ll), z(1, ll))
                            end do
                        case (2)
                            do ll = i0, indtm
                                call zaxrpy(indorb_tab(i), indorb, cphs_r(ll, i_ion), distp_(1, ll), z(1, ll))
                            end do
                            do ll = indt + 1, indt + 4
                                call zaxrpy(indorb_tab(i), indorb, cphs_r(0, i_ion), distp_(1, ll), z(1, ll))
                            end do
                        case (3)
                            do ll = i0, indtm
                                call daxrpy(indorb_tab(i), indorb, rphs(ll, i_ion), distp_(1, ll), z(1, ll))
                            end do
                        case (4)
                            do ll = i0, indtm
                                call daxrpy(indorb_tab(i), indorb, rphs(ll, i_ion), distp_(1, ll), z(1, ll))
                            end do
                            do ll = indt + 1, indt + 4
                                call daxrpy(indorb_tab(i), indorb, rphs(0, i_ion), distp_(1, ll), z(1, ll))
                            end do
                        end select
                    end if ! endif do_makefun
                end if
            end do
        end do ! enddo over j_ion
    end do ! enddo over i_ion
!$omp end parallel do
end subroutine makefun_grid

subroutine upnewwf0(indt, typecomp, nshell, ioptorb, iocc, kel, nel, r, rmu, dd, zeta, rion, distp, z, nelskip, &
                    nion, kion, iflagnorm, cnorm, LBox, rmucos, rmusin, mindist, indpar_tab, indorb_tab, indshell_tab, yesupel)

    use allio, only: ikshift, iespbc, rank, gamma_point, yes_crystalj, yes_scemama, lepsbas, slaterorb_read, nshell_det, use_qmckl, qmckl_ctx
    use Cell, only: cellscale, cellpi, rphase, phase2pi, phase2pi_down, sinphase, cosphase, s2r, car2cry
    use Constants, only: ipc
    use qmckl
#ifdef _QMCKL_GPU
    !use qmckl_gpu
#endif

    implicit none
    ! input
    integer nel, indorb, indpar, indshell, dimp, mshift &
        , j, i, i_ion, indt, nelskip, nshell, indtmin, nion, k, kk &
        , iflagnorm, typecomp, ii, jj, ll, ierr, ishift, indp1, indp2, indp3 &
        , indp4, indp5, dimx, dimy
    real*8 rmu(3, 0:0, nion), rion(3, nion), distp(27*max(nshell, nion) +&
            &nelskip*(indt + 4 + 1)), cnorm(nshell), zeta(nion), r(0:0, nion), rcos_sav(3)&
            &, mindist, rsign(3)
    ! specify the phase choice for crystal basis, ignored for Jastrow
    logical, intent(in) :: yesupel
    logical do_makefun, yeszero_z

    ! local variables
    real*8 kpip(3)
    real*8 kel(3, nel, 0:0)
    real*8 LBox, Lbox_sav
    ! in case abs(LBox).eq.3.d0, rmucos and rmusin are used as npip and rmunew
    real*8 rmucos(3, 0:0, nion), rmusin(3, 0:0, nion)
    ! real(8), dimension(:,:), allocatable :: psip_r
    logical iesjas
    integer :: ioptorb(nshell), kion(nshell)&
            &, indpar_tab(nshell + 1), indorb_tab(nshell + 1), indshell_tab(nshell + 1)
    real*8, intent(out) :: z(min(ipc, indpar_tab(1) + 2)*nelskip, 0:&
            &(1 - typecomp)*(indt + 4))
    integer :: iocc(indshell_tab(nshell + 1))
    real*8 dd(indpar_tab(nshell + 1))
    real*8 phs(3) ! scratch phase for Lbox=3
    !
    ! QMCKL
    integer*4                      :: rc
    integer*8, save                :: ao_num=0, npoints_qmckl=0
    double precision, allocatable  :: ao_vgl_qmckl(:,:), ao_value_qmckl(:,:)
    !
    integer, external :: omp_get_num_threads
    ! ---------------------------------------------------------------------
    ! orbital types are defined by the variable Lbox:
    ! Lbox = -3.d0 --> Complex open
    ! Lbox = -2.d0 --> open boundary conditions with bump orbitals
    ! Lbox = -1.d0 --> open boundary conditions
    ! Lbox =  1.d0 --> PBC/APBC old version (real w.f.)
    ! Lbox =  2.d0 --> PBC/APBC with bump orbitals
    ! Lbox =  3.d0 --> Crystal basis set (complex w.f. for PBC)
    ! ---------------------------------------------------------------------
    ! initialization
    indshell = 0
    indorb = 0
    indpar = 0
    indtmin = 0
    ! rphs = 0.d0
    ! allocate(distp_true(0:indtm,20,nshell))
    if (abs(LBox) .eq. 3.d0) then
        indp5 = nion*27 + 1
        indp1 = 2*nion + 1
        indp2 = indp1 + nion
        indp3 = indp2 + nion
        indp4 = indp3 + 3*nion
        dimx = min(ipc, indpar_tab(1) + 2)*nelskip
        dimy = (1 - typecomp)*(indt + 4)
        z = 0.d0
    elseif (yes_scemama) then
        z = 0.d0
    end if

    if (indpar_tab(1) .eq. -1) then
        iesjas = .true.
        !    indpar_tab(1)=0
        mshift = nshell_det
    else
        mshift = 0
        iesjas = .false.
    end if

    LBox_sav = LBox
    ! Crystal basis set for periodic systems
    phs(:) = 0.d0

    dimp = 20

    if (LBox .eq. 3.d0) then
        !
        ! initialize different phase for spin up/spin down electrons
        !
        if (yesupel) then
            phs(:) = phase2pi(:)
        else
            phs(:) = phase2pi_down(:)
        end if
    elseif (LBox .eq. -3.d0) then
        phs = 0.d0
    end if
    ishift = 0
    if (abs(LBox) .eq. 3.d0 .and. iesjas) then
        !
        ! always use the old periodic basis for the Jastrow
        !
        if (.not. yes_crystalj) then
            if (iespbc) then
                !       iflagnorm=2 ! always recompute distances if Jastrow is considered with Lbox.eq.3!!!
                ! REMEMBER: with Jastrow it pass from Lbox.eq.3 to Lbox.eq.1, distances
                ! to be computed are different!!!
                LBox = 1.d0 ! always use the old periodic basis for the Jastrow
            else
                Lbox = -1.d0
            end if
        else
            ishift = ikshift
        end if
        phs = 0.d0
        ! elseif(abs(Lbox).eq.3.d0) then
        ! allocate(psip_r(nelskip,i0:indt+4))
    end if

    ! electron-ion distances

    ! if((iflagnorm.gt.0.or.(LBox.eq.2.and..not.gamma_point)).and.abs(LBox).ne.3) then
    if (abs(LBox) .ne. 3) then

        do k = 1, nion
            rmu(1, 0, k) = kel(1, 1, 0) - rion(1, k)
            rmu(2, 0, k) = kel(2, 1, 0) - rion(2, k)
            rmu(3, 0, k) = kel(3, 1, 0) - rion(3, k)
        end do

        ! ****** Periodic Systems *************
        if (LBox .eq. 1.d0) then

            if (.not. gamma_point .and. .not. iesjas) then

                do k = 1, nion
                    rcos_sav(1) = dcos(rphase(1)*rmu(1, 0, k))
                    rcos_sav(2) = dcos(rphase(2)*rmu(2, 0, k))
                    rcos_sav(3) = dcos(rphase(3)*rmu(3, 0, k))
                    cosphase(0, k) = rcos_sav(1)*rcos_sav(2)*rcos_sav(3)
                    sinphase(1, 0, k) = dsin(rphase(1)*rmu(1, 0, k))*rcos_sav(2)*rcos_sav(3)
                    sinphase(2, 0, k) = dsin(rphase(2)*rmu(2, 0, k))*rcos_sav(1)*rcos_sav(3)
                    sinphase(3, 0, k) = dsin(rphase(3)*rmu(3, 0, k))*rcos_sav(1)*rcos_sav(2)
                end do

            end if

            do k = 1, nion
                rmu(1, 0, k) = rmu(1, 0, k)/cellpi(1)
                rmu(2, 0, k) = rmu(2, 0, k)/cellpi(2)
                rmu(3, 0, k) = rmu(3, 0, k)/cellpi(3)
                rmucos(1, 0, k) = dcos(rmu(1, 0, k))
                rmucos(2, 0, k) = dcos(rmu(2, 0, k))
                rmucos(3, 0, k) = dcos(rmu(3, 0, k))
                rmusin(1, 0, k) = dsin(rmu(1, 0, k))
                rmusin(2, 0, k) = dsin(rmu(2, 0, k))
                rmusin(3, 0, k) = dsin(rmu(3, 0, k))
                rmu(1, 0, k) = cellpi(1)*rmusin(1, 0, k)
                rmu(2, 0, k) = cellpi(2)*rmusin(2, 0, k)
                rmu(3, 0, k) = cellpi(3)*rmusin(3, 0, k)
                r(0, k) = max(dsqrt(rmu(1, 0, k)**2 + rmu(2, 0, k)**2 + rmu(3, 0, k)**2), mindist)
            end do

        elseif (LBox .eq. 2.d0) then

            if (iesjas) then
                do k = 1, nion
                    kpip(1) = anint(rmu(1, 0, k)/cellscale(1))
                    kpip(2) = anint(rmu(2, 0, k)/cellscale(2))
                    kpip(3) = anint(rmu(3, 0, k)/cellscale(3))
                    rmu(1, 0, k) = rmu(1, 0, k) - kpip(1)*cellscale(1)
                    rmu(2, 0, k) = rmu(2, 0, k) - kpip(2)*cellscale(2)
                    rmu(3, 0, k) = rmu(3, 0, k) - kpip(3)*cellscale(3)
                    r(0, k) = dsqrt(rmu(1, 0, k)**2 + rmu(2, 0, k)**2 + rmu(3, 0, k)**2)
                end do
            else
                do k = 1, nion
                    kpip(1) = anint(rmu(1, 0, k)/cellscale(1))
                    kpip(2) = anint(rmu(2, 0, k)/cellscale(2))
                    kpip(3) = anint(rmu(3, 0, k)/cellscale(3))
                    rmu(1, 0, k) = rmu(1, 0, k) - kpip(1)*cellscale(1)
                    rmu(2, 0, k) = rmu(2, 0, k) - kpip(2)*cellscale(2)
                    rmu(3, 0, k) = rmu(3, 0, k) - kpip(3)*cellscale(3)
                    r(0, k) = dsqrt(rmu(1, 0, k)**2 + rmu(2, 0, k)**2 + rmu(3, 0, k)**2)
                    do kk = 1, 3
                        if (phase2pi(kk) .ne. 0.d0 .and. 2*(kpip(kk)/2) .ne. kpip(kk)) r(0, k) = -r(0, k)
                    end do
                end do

            end if

        elseif (LBox .eq. -2.d0) then
            ! No need of the cutoff mindist
            do k = 1, nion
                r(0, k) = dsqrt(rmu(1, 0, k)**2 + rmu(2, 0, k)**2 + rmu(3, 0, k)**2)
            end do

        else

            do k = 1, nion
                r(0, k) = max(dsqrt(rmu(1, 0, k)**2 + rmu(2, 0, k)**2 + rmu(3, 0, k)**2), mindist)
            end do

        end if

    end if

    ! if(abs(LBox).eq.3.d0.and.iflagnorm.gt.0) then
    if (abs(LBox) .eq. 3.d0) then
        !if(iflagnorm.gt.0) then
        if (iespbc) then
!$omp parallel do default(shared) private(k)
            do k = 1, nion
                rmu(:, 0, k) = kel(:, 1, 0) - rion(:, k)
                rmucos(:, 0, k) =                                        &
                       &car2cry(:, 1)*rmu(1, 0, k) + car2cry(:, 2)*rmu(2, 0, k) + car2cry(:, 3)*rmu(3, 0, k)
                rmucos(1, 0, k) = anint(rmucos(1, 0, k)/cellscale(1))
                rmucos(2, 0, k) = anint(rmucos(2, 0, k)/cellscale(2))
                rmucos(3, 0, k) = anint(rmucos(3, 0, k)/cellscale(3))
                rmu(:, 0, k) = rmu(:, 0, k)                              &
                       & - s2r(:, 1)*rmucos(1, 0, k) - s2r(:, 2)*rmucos(2, 0, k) - s2r(:, 3)*rmucos(3, 0, k)
                r(0, k) = max(dsqrt(rmu(1, 0, k)**2 + rmu(2, 0, k)**2 +  &
                       &rmu(3, 0, k)**2), mindist)
            end do
!$omp end parallel do

        else
!$omp parallel do default(shared) private(k)
            do k = 1, nion
                rmu(1, 0, k) = kel(1, 1, 0) - rion(1, k)
                rmu(2, 0, k) = kel(2, 1, 0) - rion(2, k)
                rmu(3, 0, k) = kel(3, 1, 0) - rion(3, k)
                r(0, k) = max(dsqrt(rmu(1, 0, k)**2 + rmu(2, 0, k)**2 +  &
                       &rmu(3, 0, k)**2), mindist)
            end do
!$omp end parallel do
        end if
    end if

    !
    if (LBox .eq. 1.d0) then
!$omp parallel do default(shared)&
!$omp private(i,indpar,indorb,indshell)
        do i = 1, nshell
            indpar = max(indpar_tab(i), 0)
            indorb = indorb_tab(i)
            indshell = indshell_tab(i)
            call makefun0_pbc(ioptorb(i), iocc, indt, typecomp           &
                   &, indpar, indorb, indshell, nelskip                  &
                   &, z, dd, r(0, kion(i)), rmu(1, 0, kion(i))           &
                   &, distp(dimp*(i - 1) + 1), iflagnorm, cnorm(i)       &
                   &, rmucos(1, 0, kion(i)), rmusin(1, 0, kion(i))       &
                   &, sinphase(1, 0, kion(i)), cosphase(0, kion(i)))
        end do
!$omp end parallel do
    elseif (abs(LBox) .eq. 2.d0) then
!$omp parallel do default(shared)&
!$omp  private(i,indorb,indshell,indpar)
        do i = 1, nshell
            indpar = max(indpar_tab(i), 0)
            indorb = indorb_tab(i)
            indshell = indshell_tab(i)
            call makefun0_bump(ioptorb(i), iocc, indt, typecomp          &
                   &, indpar, indorb, indshell, nelskip, z, dd, r(0, kion(i)), rmu(1, 0, kion(i))&
                   &, distp(dimp*(i - 1) + 1), iflagnorm, cnorm(i))
        end do
!$omp end parallel do

    elseif (abs(LBox) .eq. 3.d0) then
        call makefun_grid0(distp(indp4), distp(indp5), distp(indp3), distp(indp2)&
               &, distp(indp1), distp                                    &
               &, nshell, nion, indt, typecomp, ishift, iflagnorm        &
               &, kion, ioptorb, indpar_tab, indorb_tab, indshell_tab, mindist, phs, dd, z, dimx, dimy&
               &, nelskip, iesjas, rmu, rion, cnorm, zeta, rmucos)
    else
        if (.not.use_qmckl) then
!$omp parallel do default(shared)&
!$omp private(i,ll,i_ion,indpar,indorb,indshell,do_makefun,yeszero_z)
            do i = 1, nshell
                indpar = max(indpar_tab(i), 0)
                indorb = indorb_tab(i)
                indshell = indshell_tab(i)
                do_makefun = .true.
                if (yes_scemama .and. ioptorb(i) .ne. 200) then
                    i_ion = kion(i)
                    if (slaterorb_read(i + mshift)) then
                        if (dd(indpar + 1)*r(0, i_ion) .lt. lepsbas) then
                            yeszero_z = .false.
                        else
                            yeszero_z = .true.
                        end if
                    else
                        if (dd(indpar + 1)*r(0, i_ion)*r(0, i_ion) .lt. lepsbas) then
                            yeszero_z = .false.
                        else
                            yeszero_z = .true.
                        end if
                    end if
                    do_makefun = .not. yeszero_z
                end if
                if (do_makefun) then
                    call makefun0(ioptorb(i), indt, typecomp             &
                           &, indpar, indorb, indshell, nelskip, z, dd, zeta(kion(i)), r(0, kion(i)), rmu(1, 0, kion(i))&
                           &, distp(dimp*(i - 1) + 1), iflagnorm, cnorm(i))
                end if
            end do
!$omp end parallel do
#ifdef _QMCKL
        else
            if (ao_num == 0) then
                rc = qmckl_get_ao_basis_ao_num(qmckl_ctx, ao_num)
                if (rc /= QMCKL_SUCCESS) then
                    print *, 'Error getting ao_num', rc, qmckl_ctx, ao_num
                    call abort()
                end if
            end if

            if (typecomp.eq.1) then   ! Only values
                rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, kel(1:3,1,0), 3_8)

                if (rc /= QMCKL_SUCCESS) then
                    print *, 'Error setting electron coords'
                    call abort()
                end if

                rc = qmckl_get_ao_basis_ao_value_inplace(                &
                       &qmckl_ctx,                                       &
                       &z(1,0),                                          &
                       &ao_num)

                if (rc /= QMCKL_SUCCESS) then
                    print *, 'Error getting AOs from QMCkl'
                    call abort()
                end if

            else
                allocate(ao_vgl_qmckl(ao_num, 5))

                rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, kel(1:3,1,0), 3_8)

                if (rc /= QMCKL_SUCCESS) then
                    print *, 'Error setting electron coords 2'
                    call abort()
                end if

                rc = qmckl_get_ao_basis_ao_vgl_inplace(                  &
                       &qmckl_ctx,                                       &
                       &ao_vgl_qmckl,                                    &
                       &ao_num*5_8)

                if (rc /= QMCKL_SUCCESS) then
                    print *, 'Error getting AOs from QMCkl 2'
                    call abort()
                end if

                do ii=1,ao_num
                    z(ii,0)      = ao_vgl_qmckl(ii,1)
                    z(ii,indt+1) = ao_vgl_qmckl(ii,2)
                    z(ii,indt+2) = ao_vgl_qmckl(ii,3)
                    z(ii,indt+3) = ao_vgl_qmckl(ii,4)
                    z(ii,indt+4) = ao_vgl_qmckl(ii,5)
                end do

                deallocate(ao_vgl_qmckl)
            end if
#endif
        end if
    end if

    if (iesjas) then
        LBox = LBox_sav
    end if

    return
end subroutine upnewwf0

subroutine makefun_grid0(distp_true_, distp_, rmusin, r, rphs, cphs_r&
        &, nshell, nion, indt, typecomp, ishift, iflagnorm&
        &, kion, ioptorb, indpar_tab, indorb_tab, indshell_tab, mindist, phs, dd, z, dimx, dimy&
        &, nelskip, iesjas, rmu, rion, cnorm, zeta, rmucos)
!   use allio, only : allgrid, point_shell
    use allio, only: kgrid, kgrid_atom, adr_nion, adrj_nion, ind_nion, indj_nion, lepsbas, yes_scemama, slaterorb_read, nshell_det
    use Cell, only: s2r
    use Constants, only: ipc, zimg
    implicit none
    integer, intent(in) :: nelskip, nshell, indt, nion, ishift&
            &, dimx, dimy, typecomp
    real*8, intent(in) :: mindist, phs(3), rion(3, nion), cnorm(nshell), zeta(nion)&
            &, rmu(3, 0:0, nion)
    real*8 distp_true_(0:0, 20, nion), distp_(nelskip, 0:*)
    integer i, j, ii, jj, kk, ll, indpar, indorb, indshell
    integer, intent(in) :: kion(nshell), ioptorb(nshell), indpar_tab(nshell + 1)&
            &, indorb_tab(nshell + 1), indshell_tab(nshell + 1)
    real*8, intent(in) :: dd(indpar_tab(nshell + 1)), rmucos(3, 0:0, nion)
    real*8 kpip(3)
    complex*16 cphs_r(0:0, nion)
    real*8 rphs(0:0, nion)
    real*8 rmusin(3, 0:0, nion), r(0:0, nion)
    real*8, intent(out) :: z(dimx, 0:dimy)
    integer iflagnorm, indfirst, indlast, i_ion, j_ion, nshift, mshift, case_if, case_upz
    real*8 cphs
    logical iesjas, gammaorj, do_makefun
    logical yeszero_z
    if (iesjas) then
        nshift = nion
        mshift = nshell_det
    else
        nshift = 0
        mshift = 0
    end if

    if (sum(abs(phs(:))) .eq. 0.d0) then
        gammaorj = .true.
    else
        gammaorj = .false.
    end if
    if (ipc .eq. 2 .and. .not. gammaorj) then ! complex case
        case_if = 1
    elseif (.not. gammaorj .and. ipc .eq. 1) then
        case_if = 2
    else if (ipc .eq. 2 .and. .not. iesjas) then
        case_if = 3
    else
        case_if = 4
    end if
    if (ipc .eq. 2 .and. .not. iesjas) then ! complex case
        if (typecomp .eq. 1) then
            case_upz = 1
        else
            case_upz = 2
        end if
    else ! real case
        if (typecomp .eq. 1) then
            case_upz = 3
        else
            case_upz = 4
        end if
    end if

!$omp parallel do default(shared)&
!$omp   private(i_ion,j_ion,indfirst,indlast,i,j,ii,jj,kk,ll,indpar,indorb,indshell,kpip,cphs,do_makefun,yeszero_z)
    do i_ion = 1, nion
        if (iesjas) then
            indfirst = adrj_nion(i_ion)
            indlast = adrj_nion(i_ion + 1) - 1
        else
            indfirst = adr_nion(i_ion)
            indlast = adr_nion(i_ion + 1) - 1
        end if
        do ii = 1, kgrid_atom(i_ion + nshift)%dimshell
            kpip(1) = kgrid_atom(i_ion + nshift)%kpip(1, ii)
            kpip(2) = kgrid_atom(i_ion + nshift)%kpip(2, ii)
            kpip(3) = kgrid_atom(i_ion + nshift)%kpip(3, ii)
            rmusin(:, 0, i_ion) = rmu(:, 0, i_ion) + &
          &s2r(:, 1)*kpip(1) + s2r(:, 2)*kpip(2) + s2r(:, 3)*kpip(3)
            r(0, i_ion) = max(dsqrt(rmusin(1, 0, i_ion)**2 + rmusin(2, 0, i_ion)**2 + rmusin(3, 0, i_ion)**2), mindist)
            select case (case_if)

            case (1)
                cphs = -(phs(1)*(kpip(1) - rmucos(1, 0, i_ion)) + &
                         phs(2)*(kpip(2) - rmucos(2, 0, i_ion)) + &
                         phs(3)*(kpip(3) - rmucos(3, 0, i_ion)))
                cphs_r(0, i_ion) = dcmplx(dcos(cphs), dsin(cphs))
            case (2)
                rphs(0, i_ion) = dcos(phs(1)*(kpip(1) - rmucos(1, 0, i_ion)) + &
                                      phs(2)*(kpip(2) - rmucos(2, 0, i_ion)) + &
                                      phs(3)*(kpip(3) - rmucos(3, 0, i_ion)))
            case (3)
                cphs_r(0, i_ion) = dcmplx(1.d0, 0.d0)
            case (4)
                rphs(0, i_ion) = 1.d0
            end select
            do j_ion = indfirst, indlast
                if (iesjas) then
                    i = indj_nion(j_ion)
                else
                    i = ind_nion(j_ion)
                end if

                if (kgrid(i + ishift)%tobedone(ii)) then

                    indpar = max(indpar_tab(i), 0)
                    indorb = indorb_tab(i)
                    indshell = indshell_tab(i)
                    do_makefun = .true.
                    if (yes_scemama .and. ioptorb(i) .ne. 200) then
!           do_makefun=.false.
                        if (slaterorb_read(i + mshift)) then
                            if (dd(indpar + 1)*r(0, i_ion) .lt. lepsbas) then
!              do_makefun=.true.
                                yeszero_z = .false.
                            else
                                yeszero_z = .true.
                            end if
                        else
                            if (dd(indpar + 1)*r(0, i_ion)*r(0, i_ion) .lt. lepsbas) then
!             do_makefun=.true.
                                yeszero_z = .false.
                            else
                                yeszero_z = .true.
                            end if
                        end if
                        do_makefun = .not. yeszero_z
                    end if
                    if (do_makefun) then
                        call makefun0(ioptorb(i), indt, typecomp&
                                &, indpar, indorb, indshell, nelskip, distp_, dd, zeta(i_ion), r(0, i_ion), rmusin(1, 0, i_ion)&
                                &, distp_true_(0, 1, i_ion), iflagnorm, cnorm(i))

                        select case (case_upz)
!           if(ipc.eq.2.and..not.iesjas) then ! complex case
!               if(typecomp.eq.1) then
                        case (1)
                            call zaxrpy(indorb_tab(i), indorb, cphs_r(0, i_ion), distp_(1, 0), z(1, 0))
                        case (2)
                            call zaxrpy(indorb_tab(i), indorb, cphs_r(0, i_ion), distp_(1, 0), z(1, 0))
                            do ll = indt + 1, indt + 4
                                call zaxrpy(indorb_tab(i), indorb, cphs_r(0, i_ion), distp_(1, ll), z(1, ll))
                            end do
                        case (3)
                            call daxrpy(indorb_tab(i), indorb, rphs(0, i_ion), distp_(1, 0), z(1, 0))
                        case (4)
                            call daxrpy(indorb_tab(i), indorb, rphs(0, i_ion), distp_(1, 0), z(1, 0))
                            do ll = indt + 1, indt + 4
                                call daxrpy(indorb_tab(i), indorb, rphs(0, i_ion), distp_(1, ll), z(1, ll))
                            end do
                        end select
!               endif
!           endif
                    end if ! endif do_makefun
                end if
            end do
        end do ! enddo over j_ion
    end do ! enddo over i_ion
!$omp end parallel do
end subroutine makefun_grid0

subroutine zaxrpy(n, m, a, x, z)
    implicit none
    integer n, m, l, i
    complex*16 z(*), a
    real*8 x(*)
    l = m - n
    select case (l)
    case (1)
        z(n + 1) = z(n + 1) + a*x(n + 1)
    case default
        do i = n + 1, m
            z(i) = z(i) + a*x(i)
        end do
    end select
    return
end
subroutine daxrpy(n, m, a, x, z)
    implicit none
    integer n, m, l, i
    real*8 z(*), a
    real*8 x(*)
    l = m - n
    select case (l)
    case (1)
        z(n + 1) = z(n + 1) + a*x(n + 1)
    case default
        do i = n + 1, m
            z(i) = z(i) + a*x(i)
        end do
    end select
    return
end
