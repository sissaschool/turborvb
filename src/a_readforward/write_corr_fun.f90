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
!

subroutine write_corr_fun(ebin, ebin2, wbin, ibin, ibinit, nind, maxf, ddim, ell, nel, nelup &
                          , nrhoind, dxil, ind_offset, psip, read_start, ncell &
                          , ifrho, ifspin, ifkspin &
                          , iespbc, atom_number, datagrid &
                          , npairind, dxil_p, ngrid_p, ind_offset_pair, ifpair, iffluct, drmax &
                          , r_offset, vdim, nskind, ifsofk, fermi_flag, rkcomp, ndim, simap, ifcorrs &
                          , sphere_radius, sec_spc, nshellsp, iioniond, allshells)

    use constants
    use grid_module
    use Assar_module
    use qpwf_module, only: ifqpwf, ifqpwf_h, ifqpwf_k, ifqpwf_extr, decouple_k, decouple_files
    use Spin2, only: ifspin2
    use allio, only: mult, dists, nion, rionsav, rion, chara, rank, wkp, commcolrep_mpi, nproc, nk, nw
    use cell, only: s2r, yes_tilted, recip
    use Rho_corr_module, only: ifrho_corr
    use Dipole_module, only: ifdipole, ndipole, nquad, quad, quad_diag, diag_quad, nuclear_dipole
    use berry_phase, only: ifberry

    implicit none
    integer i1, j1, ibin, ibinit, maxf, nmis, ddim, i, k, kk, l, shift, nel &
        , jj, ix(3), ind_offset(*), nrhoind, j, read_start, nind, ncell(*) &
        , npairind, ind_offset_pair(*), vdim(*), ii, nelup, nbin &
        , ind, ind_last, ngrid_p, sec_spc(*), sec, nshellsp, kmom, iioniond(*), ll, m, mm
    logical blank, ifrho, iespbc, ifpair, iffluct, ifsofk, fermi_flag, ifspin, ifkspin, ifcorrs &
        , allow_write, allshells
    real*8 psip(0:maxf, *), ebin(0:maxf, *), ebin2(0:maxf, *), wbin(0:*) &
        , rx(3), rx_old(3), dxil(*), dv, eps, ell(*), vol, origin(3), cell_loc(3), cellscale(3) &
        , atom_number(nion), datagrid(*), dxil_p(*), dvp, r_offset(3) &
        , rkcomp(ndim, *), antot, signk, anorm(4), av(4), err(4), vgrid, pair_norm, rho, drmax &
        , prefac_g, vold, vnew, rv, vnorm, sang, sphere_radius, tmp_norm, vecscra(3), s2r_loc(3, 3) &
        , knorm, knorm_cry
    integer :: mesh(3), nskind, ndim, simap(*), mesh_plot_3d(3)
    character(14) :: word
    character(1) :: cartesian(3)
    real(8), dimension(:), allocatable :: check

    allocate (check(0:maxf))

    eps = 1d-7

    cartesian(1) = 'x'
    cartesian(2) = 'y'
    cartesian(3) = 'z'

    ! averages
    nmis = ibin - ibinit + 1

    do i = 0, maxf
        do kk = 1, nind
            psip(i, kk) = ebin(i, kk)/wbin(i)
            psip(i, nind + kk) = dsqrt(dabs(ebin2(i, kk)/wbin(i) - psip(i, kk)**2))
            psip(i, nind + kk) = psip(i, nind + kk)/dsqrt(dfloat(nmis))
            !      write(*,*) psip(i,kk),psip(i,nind+kk)
        end do
    end do

    allow_write = .true.
    if (rank .ne. 0) allow_write = .false.
    !average psip over processors using collrep_mpi communicator
    kmom = (rank + nproc/nk)/(nproc/nk)
    !     if(ifberry.and.decouple_k.and..not.decouple_files) then
    !      do i=1,nind,2
    !      !The x components are substitute with the modulus and the error is propagated n=sqrt(x^2+yˆ2), dn^2= (xdx/n)^2+(ydy/n)^2
    !!      check(0:maxf)=sqrt(psip(0:maxf,i)**2+psip(0:maxf,i+1)**2) ! this is n
    !!      psip(0:maxf,nind+i)=(wkp(kmom)*&
    !     &(psip(0:maxf,i)*psip(0:maxf,nind+i)/check(0:maxf)))**2+&
    !     &(wkp(kmom)*(psip(0:maxf,i+1)*psip(0:maxf,nind+i+1)/check(0:maxf)))**2
    !      !the pure immaginary part is unchanged
    !      psip(0:maxf,nind+i+1)=(wkp(kmom)*psip(0:maxf, nind+i+1))**2
    !      psip(0:maxf,i)=check(0:maxf)*wkp(kmom)
    !      psip(0:maxf,i+1)=psip(0:maxf,i+1)*wkp(kmom)
    !      Both real and  immaginary parts are unchanged
    !       psip(0:maxf,i)=psip(0:maxf,i)*wkp(kmom)
    !       psip(0:maxf,nind+i)=(wkp(kmom)*psip(0:maxf, nind+i))**2
    !       psip(0:maxf,i+1)=psip(0:maxf,i+1)*wkp(kmom)
    !       psip(0:maxf,nind+i+1)=(wkp(kmom)*psip(0:maxf, nind+i+1))**2
    !      enddo
    !
    if (decouple_k .and. .not. decouple_files) then
        psip(0:maxf, 1:nind) = psip(0:maxf, 1:nind)*wkp(kmom)
        psip(0:maxf, nind + 1:2*nind) = (psip(0:maxf, nind + 1:2*nind)*wkp(kmom))**2
#ifdef PARALLEL
        call reduce_base_real((maxf + 1)*2*nind, psip, commcolrep_mpi, -1)
#endif
        psip(0:maxf, nind + 1:2*nind) = dsqrt(psip(0:maxf, nind + 1:2*nind))
    end if

    if (ifrho .or. ifspin) then
        do i = 1, 3
            mesh(i) = nint(ell(vdim(i))*dxil(i))
            mesh_plot_3d(i) = mesh(i) + 1
            cell_loc(i) = ell(vdim(i))
            cellscale(i) = ell(vdim(i))*ncell(vdim(i))
            if (yes_tilted) s2r_loc(:, i) = s2r(:, i)/ncell(vdim(i))
        end do

        origin = 0.d0

        vol = 1.d0
        dv = 1.d0
        do i = 1, ddim
            ii = vdim(i)
            vol = vol*ell(ii)*ncell(ii)
            dv = dv*dxil(ii)
        end do

        vgrid = 1./dv

        rho = dble(nel)/vol
        pair_norm = 1.d0/nel/rho
        dvp = dv
        dvp = dvp*pair_norm

        ! valid for 2d and 3d
        sang = (ndim - 1)*2.d0*pi
        vnorm = sang/ndim/drmax**ndim

        if (ifspin2) dv = 1.d0
    end if

    ! Here output the measured properties. Each property has its own subroutine
    ! which contains open, rewind, write and close.
    ! The output order should exactly follow the collecting order in readforward
    ! because the order affects the shift.
    ! Hints: to have a check, simply uncomment the lines which output the
    !         shift here and those with 'check xxx index' in readforward

    shift = 0

    if (ifrho .and. allow_write) then
        call write_rho
    end if
    !write(6,*) 'check write rho index',shift

    if (ifspin .and. allow_write) then
        call write_spin
    end if
    !write(6,*) 'check write spin index',shift

    if (ifsofk .and. allow_write) then
        call write_sofk
    end if
    !write(6,*) 'check write sofk index',shift

    ! NB notice here allow_write is not present because averages and error bars
    ! are computed inside write_corrs
    if (ifcorrs) then
        call write_corrs
    end if
    !write(6,*) 'check write corrs index',shift

    if (ifrho_assar .and. allow_write) then
        call write_rho_assar
    end if

    if (ifspin2 .and. allow_write) then
        call write_spin2
    end if
    !write(6,*) 'check write spin2 index',shift

    if (ifqpwf .and. allow_write) then
        call write_qpwf
    end if
    !write(6,*) 'check write stm index',shift

    if (ifrho_corr .and. allow_write) then
        call write_rho_corr
    end if

    if (ifdipole .and. allow_write) then
        call write_dipole
    end if
    !write(6,*) 'check write spin index',shift
    if (ifberry) then
        call write_berry
    end if

    do kk = 1, 2*nind
        do i = 0, maxf
            psip(i, kk) = 0.d0
        end do
    end do

    deallocate (check)

contains

    subroutine write_rho
        implicit none
        integer ip, jp, kp

        if (ifpair) then
            if (decouple_files) then
                open (read_start + 86, file='rhorhoK'//trim(chara), form='formatted', status='unknown')
            else
                open (read_start + 86, file='rhorho.dat', form='formatted', status='unknown')
            end if
            rewind (read_start + 86)
            if (ngrid_p .ne. 0) then
                if (decouple_files) then
                    open (read_start + 87, file='rhorho_radialK'//trim(chara), form='formatted', status='unknown')
                else
                    open (read_start + 87, file='rhorho_radial.dat', form='formatted', status='unknown')
                end if
                rewind (read_start + 87)
            end if
            if (sphere_radius .ne. 0.d0) then
                if (decouple_files) then
                    open (read_start + 85, file='rho_ionK'//trim(chara), form='formatted', status='unknown')
                else
                    open (read_start + 85, file='rho_ion.dat', form='formatted', status='unknown')
                end if
                rewind (read_start + 85)
            else
                if (decouple_files) then
                    open (read_start + 85, file='rhoK'//trim(chara), form='formatted', status='unknown')
                else
                    open (read_start + 85, file='rho.dat', form='formatted', status='unknown')
                end if
                rewind (read_start + 85)
            end if
        else
            if (decouple_files) then
                open (read_start + 85, file='rhoK'//trim(chara), form='formatted', status='unknown')
            else
                open (read_start + 85, file='rho.dat', form='formatted', status='unknown')
            end if
            rewind (read_start + 85)
        end if

        if (sphere_radius .eq. 0.d0) then

            rx_old = 0.d0
            check = 0.d0

            write (read_start + 85, *) ' Volume unit cell grid ', vgrid, nrhoind

            do k = 1, nrhoind

                jj = k
                blank = .true.
                do j = ddim, 1, -1
                    ix(j) = (jj - 1)/ind_offset(j) + 1
!                    rx(j) = (ix(j) - 0.5d0) / dxil(j) + r_offset(j)
                    rx(j) = (ix(j) - 0.5d0)/dxil(j)
                    if (abs(rx(j) - rx_old(j)) .lt. eps) blank = .false.
                    jj = jj - (ix(j) - 1)*ind_offset(j)
                end do

                if (yes_tilted) then
                    vecscra(:) = 0.d0
                    do j = ddim, 1, -1
                        vecscra(:) = vecscra(:) + s2r(:, vdim(j))*rx(j)/cellscale(j)
                    end do
                    vecscra(:) = vecscra(:) - r_offset(:)
                else
                    do j = ddim, 1, -1
                        vecscra(j) = rx(j) - r_offset(vdim(j))
                    end do
                end if

                check(0:maxf) = check(0:maxf) + psip(0:maxf, k + shift)

                if (blank .and. k .ne. 1 .and. ddim .ne. 1) write (read_start + 85, *)

#ifdef __KCOMP
                if (yes_tilted) then
                    write (read_start + 85, '(32767f14.8)') (vecscra(j), j=1, 3) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                else
                    write (read_start + 85, '(32767f14.8)') (vecscra(j), j=1, ddim) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                end if
#else
                if (yes_tilted) then
                    write (read_start + 85, '(1000000f14.8)') (vecscra(j), j=1, 3) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                else
                    write (read_start + 85, '(1000000f14.8)') (vecscra(j), j=1, ddim) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                end if
#endif
                rx_old = rx

            end do

            if (ddim .eq. 3) then
                if (decouple_files) then
                    word = 'chargeK'//trim(chara)//'M'
                else
                    word = 'charge'
                end if
                datagrid(1:mesh_plot_3d(1)*mesh_plot_3d(2)*mesh_plot_3d(3)) = 0.d0

                do i = maxf, 0, -1
                    j = maxf - i + 1
                    do k = 1, nrhoind

                        jj = k
                        do kk = ddim, 1, -1
                            ix(kk) = (jj - 1)/ind_offset(kk) + 1
                            jj = jj - (ix(kk) - 1)*ind_offset(kk)
                        end do

                        ind = ix(1) + (ix(2) - 1)*mesh_plot_3d(1) + (ix(3) - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                        datagrid(ind) = dv*psip(i, k + shift)

                    end do

                    ! periodic condition on the grid
                    do ii = 1, mesh(1) + 1
                        do jj = 1, mesh(2) + 1
                            do kk = 1, mesh(3) + 1

                                ip = ii
                                if (ip .gt. mesh(1)) ip = 1
                                jp = jj
                                if (jp .gt. mesh(2)) jp = 1
                                kp = kk
                                if (kp .gt. mesh(3)) kp = 1

                                ind = ip + (jp - 1)*mesh_plot_3d(1) + (kp - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)
                                ind_last = ii + (jj - 1)*mesh_plot_3d(1) + (kk - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                                if (ind_last .ne. ind) datagrid(ind_last) = datagrid(ind)

                            end do
                        end do
                    end do

                    if (yes_tilted) then
                        call plot_3d_data_tilted(1, s2r_loc, s2r, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                             &, origin, datagrid, j, word)
                    else
                        call plot_3d_data(1, cell_loc, cellscale, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                             &, origin, datagrid, j, word)
                    end if

                end do

                if (iffluct) then
                    ! plotting charge fluctuations
                    if (decouple_files) then
                        word = 'chargeflK'//trim(chara)//'M'
                    else
                        word = 'chargefluct'
                    end if
                    datagrid(1:mesh_plot_3d(1)*mesh_plot_3d(2)*mesh_plot_3d(3)) = 0.d0

                    do i = maxf, 0, -1
                        j = maxf - i + 1
                        do k = 1, nrhoind

                            jj = k
                            do kk = ddim, 1, -1
                                ix(kk) = (jj - 1)/ind_offset(kk) + 1
                                jj = jj - (ix(kk) - 1)*ind_offset(kk)
                            end do

                            ind = ix(1) + (ix(2) - 1)*mesh_plot_3d(1) + (ix(3) - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                            ! fluctuations assuming uncorrelated bins and walkers
                            datagrid(ind) = dv*psip(i, k + shift + nind)*dsqrt(dfloat(nmis))*dsqrt(dfloat(nw))

                        end do

                        ! periodic condition on the grid
                        do ii = 1, mesh(1) + 1
                            do jj = 1, mesh(2) + 1
                                do kk = 1, mesh(3) + 1

                                    ip = ii
                                    if (ip .gt. mesh(1)) ip = 1
                                    jp = jj
                                    if (jp .gt. mesh(2)) jp = 1
                                    kp = kk
                                    if (kp .gt. mesh(3)) kp = 1

                                    ind = ip + (jp - 1)*mesh_plot_3d(1) + (kp - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)
                                    ind_last = ii + (jj - 1)*mesh_plot_3d(1) + (kk - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                                    if (ind_last .ne. ind) datagrid(ind_last) = datagrid(ind)

                                end do
                            end do
                        end do

                        if (yes_tilted) then
                            call plot_3d_data_tilted(1, s2r_loc, s2r, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                                 &, origin, datagrid, j, word)
                        else
                            call plot_3d_data(1, cell_loc, cellscale, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                                 &, origin, datagrid, j, word)
                        end if

                    end do

                end if

            end if

            shift = shift + nrhoind

            if (rank .eq. 0) write (6, *) 'total charge in the system', check(maxf)

        else

            check = 0.d0

            write (read_start + 85, *) ' Number of ion sites ', nion

            do k = 1, nion

                check(0:maxf) = check(0:maxf) + psip(0:maxf, k + shift)

#ifdef __KCOMP
                write (read_start + 85, '(32767f14.8)') (rion(j, k), j=1, 3) &
                    , (psip(i, k + shift), psip(i, k + shift + nind), i=maxf, 0, -1)
#else
                write (read_start + 85, '(1000000f14.8)') (rion(j, k), j=1, 3) &
                    , (psip(i, k + shift), psip(i, k + shift + nind), i=maxf, 0, -1)
#endif

            end do

            shift = shift + nrhoind

            if (rank .eq. 0) write (6, *) 'total charge within spheres', check(maxf)
            ! compute <n_i>^2 averaged over sites
            check = check/dble(nion)
            check = check**2

        end if

        if (ifpair) then

            if (sphere_radius .eq. 0.d0) then

                do k = 1, npairind

                    jj = k
                    do j = ddim, 1, -1
                        ix(j) = (jj - 1)/ind_offset_pair(j) + 1
                        rx(j) = (ix(j) - 0.5d0)/dxil_p(j)
                        jj = jj - (ix(j) - 1)*ind_offset_pair(j)
                    end do

#ifdef __KCOMP
                    write (read_start + 86, '(32767f20.8)') &
#else
                        write (read_start + 86, '(1000000f22.8)') &
#endif
                        (rx(j), j=1, ddim), (dvp*psip(i, k + shift), dvp*psip(i, k + shift + nind), i=maxf, 0, -1)

                end do

                shift = shift + npairind

                ! radial g
                vold = 0.d0

                do k = 1, ngrid_p

                    rv = (dble(k) - 0.5d0)/drmax
                    vnew = vnorm*dble(k)**ndim
                    prefac_g = pair_norm/(vnew - vold)

#ifdef __KCOMP
                    write (read_start + 87, '(32767f22.8)') rv &
                        , (prefac_g*psip(i, k + shift), prefac_g*psip(i, k + shift + nind), i=maxf, 0, -1)
#else
                    write (read_start + 87, '(100000f22.8)') rv &
                        , (prefac_g*psip(i, k + shift), prefac_g*psip(i, k + shift + nind), i=maxf, 0, -1)
#endif
                    vold = vnew

                end do

                shift = shift + ngrid_p

            else

                do k = 1, npairind - nshellsp

                    tmp_norm = 1.d0/mult(k)

                    ! <n_i n_j>
                    if (allshells) then

                        j1 = (iioniond(k) - 1)/nion + 1
                        i1 = iioniond(k) - (j1 - 1)*nion
                        if (i1 .le. j1) then
#ifdef __KCOMP
                            write (read_start + 86, '(2I9,32767f22.8)') &
#else
                                write (read_start + 86, '(2I9,1000000f22.8)') &
#endif
                                ! <n_i n_j> - <n_i> <n_j>
                                i1, j1, dists(k) &
                                , (tmp_norm*psip(i, k + shift), tmp_norm*psip(i, k + shift + nind), &
                                   tmp_norm*psip(i, k + shift) - psip(i, shift - nrhoind + i1)*psip(i, shift - nrhoind + j1) &
                                   , tmp_norm*psip(i, k + shift + nind), i=maxf, 0, -1)
                        end if
                    else
#ifdef __KCOMP
                        write (read_start + 86, '(32767f22.8)') &
#else
                            write (read_start + 86, '(1000000f22.8)') &
#endif
                            dists(k), (tmp_norm*psip(i, k + shift), tmp_norm*psip(i, k + shift + nind), &
                                       ! <n_i n_j> - <n_i> <n_j>
                                       tmp_norm*psip(i, k + shift) - check(i), tmp_norm*psip(i, k + shift + nind), i=maxf, 0, -1)
                    end if
                end do

                sec = 0
                sec = sec_spc(npairind - nshellsp + 1)
                write (read_start + 86, *) '  '
                write (read_start + 86, *) '  '
                write (read_start + 86, *) ' species-species sector:', sec

                do k = npairind - nshellsp + 1, npairind

                    tmp_norm = 1.d0/mult(k)

                    if (sec_spc(k) .ne. sec) then
                        sec = sec_spc(k)
                        write (read_start + 86, *) '  '
                        write (read_start + 86, *) '  '
                        write (read_start + 86, *) ' species-species sector:', sec
                    end if

#ifdef __KCOMP
                    write (read_start + 86, '(32767f20.8)') &
#else
                        write (read_start + 86, '(1000000f22.8)') &
#endif
                        ! <n_i n_j>
                        dists(k), (tmp_norm*psip(i, k + shift), tmp_norm*psip(i, k + shift + nind), i=maxf, 0, -1)

                end do

                shift = shift + npairind

            end if

        end if

        close (read_start + 85)
        if (ifpair) then
            close (read_start + 86)
            if (ngrid_p .ne. 0) close (read_start + 87)
        end if

    end subroutine write_rho

    subroutine write_spin
        implicit none
        integer ip, jp, kp
        if (ifpair) then
            if (decouple_files) then
                open (read_start + 89, file='spinspinK'//trim(chara), form='formatted', status='unknown')
            else
                open (read_start + 89, file='spinspin.dat', form='formatted', status='unknown')
            end if
            rewind (read_start + 89)
            if (ngrid_p .ne. 0) then
                if (decouple_files) then
                    open (read_start + 90, file='spinspin_radialK'//trim(chara), form='formatted', status='unknown')
                else
                    open (read_start + 90, file='spinspin_radial.dat', form='formatted', status='unknown')
                end if
                rewind (read_start + 90)
            end if
            if (sphere_radius .ne. 0.d0) then
                if (decouple_files) then
                    open (read_start + 88, file='spin_ionK'//trim(chara), form='formatted', status='unknown')
                else
                    open (read_start + 88, file='spin_ion.dat', form='formatted', status='unknown')
                end if
                rewind (read_start + 88)
            else
                if (decouple_files) then
                    open (read_start + 88, file='spinK'//trim(chara), form='formatted', status='unknown')
                else
                    open (read_start + 88, file='spin.dat', form='formatted', status='unknown')
                end if
                rewind (read_start + 88)
            end if
        else
            if (decouple_files) then
                open (read_start + 88, file='spinK'//trim(chara), form='formatted', status='unknown')
            else
                open (read_start + 88, file='spin.dat', form='formatted', status='unknown')
            end if
            rewind (read_start + 88)
        end if

        if (sphere_radius .eq. 0.d0) then

            rx_old = 0.d0
            check = 0.d0
            write (read_start + 88, *) ' Volume unit cell grid ', vgrid, nrhoind

            do k = 1, nrhoind

                jj = k
                blank = .true.
                do j = ddim, 1, -1
                    ix(j) = (jj - 1)/ind_offset(j) + 1
!                    rx(j) = (ix(j) - 0.5d0) / dxil(j) + r_offset(j)
                    rx(j) = (ix(j) - 0.5d0)/dxil(j)
                    if (abs(rx(j) - rx_old(j)) .lt. eps) blank = .false.
                    jj = jj - (ix(j) - 1)*ind_offset(j)
                end do

                if (yes_tilted) then
                    vecscra(:) = 0.d0
                    do j = ddim, 1, -1
                        vecscra(:) = vecscra(:) + s2r(:, vdim(j))*rx(j)/cellscale(j)
                    end do
                    vecscra(:) = vecscra(:) - r_offset(:)
                else
                    do j = ddim, 1, -1
                        vecscra(j) = rx(j) - r_offset(vdim(j))
                    end do
                end if

                check(0:maxf) = check(0:maxf) + psip(maxf, k + shift)

                if (blank .and. k .ne. 1 .and. ddim .ne. 1) write (read_start + 88, *)
#ifdef __KCOMP

!          write(read_start+88,'(32767f14.8)') (rx(j)-0.5d0*ell(vdim(j)),j=1,ddim)&
                !,(dv*psip(i,k+shift),dv*psip(i,k+shift+nind),i=maxf,0,-1)
                ! MC 13/05/2022 to make it consistent with the charge density
                if (yes_tilted) then
                    write (read_start + 88, '(32767f14.8)') (vecscra(j), j=1, 3) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                else
                    write (read_start + 88, '(32767f14.8)') (vecscra(j), j=1, ddim) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                end if
#else
!                write(read_start + 88, '(1000000f14.8)') (rx(j) - 0.5d0 * ell(vdim(j)), j=1, ddim)&
!               ,(dv*psip(i, k+shift), dv*psip(i, k+shift+nind), i = maxf, 0, -1)
! MC 13/05/2022 to make it consistent with the charge density
                if (yes_tilted) then
                    write (read_start + 88, '(1000000f14.8)') (vecscra(j), j=1, 3) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                else
                    write (read_start + 88, '(1000000f14.8)') (vecscra(j), j=1, ddim) &
                        , (dv*psip(i, k + shift), dv*psip(i, k + shift + nind), i=maxf, 0, -1)
                end if
#endif

                rx_old = rx

            end do

            if (ifkspin) then

                write (read_start + 88, *) ' Total structure factor  real (mu_B) ='
#ifdef __KCOMP
                write (read_start + 88, '(32767f14.8)') (psip(i, shift + nrhoind + 1)&
                    &, psip(i, shift + nrhoind + 1 + nind), i=maxf, 0, -1)
#else
                write (read_start + 88, '(1000000f14.8)') (psip(i, shift + nrhoind + 1)&
                        &, psip(i, shift + nrhoind + 1 + nind), i=maxf, 0, -1)
#endif
                write (read_start + 88, *) ' Total structure factor  imaginary (mu_B) ='
#ifdef __KCOMP
                write (read_start + 88, '(32767f14.8)') (psip(i, shift + nrhoind + 2)&
                    &, psip(i, shift + nrhoind + 2 + nind), i=maxf, 0, -1)
#else
                write (read_start + 88, '(1000000f14.8)') (psip(i, shift + nrhoind + 2)&
                        &, psip(i, shift + nrhoind + 2 + nind), i=maxf, 0, -1)
#endif

            else

                !  Computing total magnetization on the grid. The wider the grid the more
                !  reliable is the calculation.
                do i = 0, maxf
                    psip(i, shift + nrhoind + 1) = sum(abs(psip(i, shift + 1:shift + nrhoind)))
                    psip(i, shift + nrhoind + 1 + nind) = &
                            & dsqrt(sum(psip(i, shift + nind + 1:shift + nind + nrhoind)**2))
                end do

                write (read_start + 88, *) ' Total magnetization (2 mu_B) ='
                do i = maxf, 0, -1
                    write (read_start + 88, '(I12,2f14.8)') maxf - i, psip(i, shift + nrhoind + 1)&
                            &, psip(i, shift + nrhoind + 1 + nind)
                end do

            end if

            if (ddim .eq. 3) then

                datagrid(1:mesh_plot_3d(1)*mesh_plot_3d(2)*mesh_plot_3d(3)) = 0.d0

                if (decouple_files) then
                    word = 'spinK'//trim(chara)//'M'
                else
                    word = 'spin'
                end if
                do i = maxf, 0, -1
                    j = maxf - i + 1
                    do k = 1, nrhoind

                        jj = k
                        do kk = ddim, 1, -1
                            ix(kk) = (jj - 1)/ind_offset(kk) + 1
                            jj = jj - (ix(kk) - 1)*ind_offset(kk)
                        end do

                        ind = ix(1) + (ix(2) - 1)*mesh_plot_3d(1) + (ix(3) - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                        datagrid(ind) = dv*psip(i, k + shift)

                    end do

                    ! periodic condition on the grid
                    do ii = 1, mesh(1) + 1
                        do jj = 1, mesh(2) + 1
                            do kk = 1, mesh(3) + 1

                                ip = ii
                                if (ip .gt. mesh(1)) ip = 1
                                jp = jj
                                if (jp .gt. mesh(2)) jp = 1
                                kp = kk
                                if (kp .gt. mesh(3)) kp = 1

                                ind = ip + (jp - 1)*mesh_plot_3d(1) + (kp - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)
                                ind_last = ii + (jj - 1)*mesh_plot_3d(1) + (kk - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                                if (ind_last .ne. ind) datagrid(ind_last) = datagrid(ind)

                            end do
                        end do
                    end do

                    if (yes_tilted) then
                        call plot_3d_data_tilted(1, s2r_loc, s2r, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                            &, origin, datagrid, j, word)
                    else
                        call plot_3d_data(1, cell_loc, cellscale, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                            &, origin, datagrid, j, word)
                    end if

                end do

                if (iffluct) then
                    ! plotting spin fluctuations
                    if (decouple_files) then
                        word = 'spinflK'//trim(chara)//'M'
                    else
                        word = 'spinfluct'
                    end if
                    do i = maxf, 0, -1
                        j = maxf - i + 1
                        do k = 1, nrhoind

                            jj = k
                            do kk = ddim, 1, -1
                                ix(kk) = (jj - 1)/ind_offset(kk) + 1
                                jj = jj - (ix(kk) - 1)*ind_offset(kk)
                            end do

                            ind = ix(1) + (ix(2) - 1)*mesh_plot_3d(1) + (ix(3) - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                            datagrid(ind) = dv*psip(i, k + shift + nind)*dsqrt(dfloat(nmis))*dsqrt(dfloat(nw))

                        end do

                        ! periodic condition on the grid
                        do ii = 1, mesh(1) + 1
                            do jj = 1, mesh(2) + 1
                                do kk = 1, mesh(3) + 1

                                    ip = ii
                                    if (ip .gt. mesh(1)) ip = 1
                                    jp = jj
                                    if (jp .gt. mesh(2)) jp = 1
                                    kp = kk
                                    if (kp .gt. mesh(3)) kp = 1

                                    ind = ip + (jp - 1)*mesh_plot_3d(1) + (kp - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)
                                    ind_last = ii + (jj - 1)*mesh_plot_3d(1) + (kk - 1)*mesh_plot_3d(1)*mesh_plot_3d(2)

                                    if (ind_last .ne. ind) datagrid(ind_last) = datagrid(ind)

                                end do
                            end do
                        end do

                        if (yes_tilted) then
                            call plot_3d_data_tilted(1, s2r_loc, s2r, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                                &, origin, datagrid, j, word)
                        else
                            call plot_3d_data(1, cell_loc, cellscale, nion, rionsav, atom_number, iespbc, mesh_plot_3d&
                                &, origin, datagrid, j, word)
                        end if

                    end do

                end if

            end if

            if (ifkspin) then
                shift = shift + nrhoind + 2
            else
                shift = shift + nrhoind + 1
            end if

            if (rank .eq. 0) write (6, *) 'total spin in the system', check(maxf)

        else

            check = 0.d0

            write (read_start + 88, *) ' Number of ion sites ', nion

            do k = 1, nion

                check(0:maxf) = check(0:maxf) + psip(maxf, k + shift)

#ifdef __KCOMP
                write (read_start + 88, '(32767f14.8)') (rion(j, k), j=1, 3) &
                    , (psip(i, k + shift), psip(i, k + shift + nind), i=maxf, 0, -1)
#else
                write (read_start + 88, '(1000000f14.8)') (rion(j, k), j=1, 3) &
                    , (psip(i, k + shift), psip(i, k + shift + nind), i=maxf, 0, -1)
#endif

            end do

            !       write(read_start+88,*) ' Total magnetization (2 mu_B) ='
            !#ifdef __KCOMP
            !       write(read_start+88,'(32767f14.8)') (psip(i,shift+nrhoind+1)&
            !               &,psip(i,shift+nrhoind+1+nind),i=maxf,0,-1)
            !#else
            !       write(read_start+88,'(1000000f14.8)') (psip(i,shift+nrhoind+1)&
            !            &,psip(i,shift+nrhoind+1+nind),i=maxf,0,-1)
            !#endif

            shift = shift + nrhoind + 1

            if (rank .eq. 0) write (6, *) 'total spin within spheres', check(maxf)
            ! compute <sigma_i>^2 averaged over sites
            check = check/dble(nion)
            check = check**2

        end if

        if (ifpair) then

            if (sphere_radius .eq. 0.d0) then

                do k = 1, npairind

                    jj = k
                    do j = ddim, 1, -1
                        ix(j) = (jj - 1)/ind_offset_pair(j) + 1
                        rx(j) = (ix(j) - 0.5d0)/dxil_p(j)
                        jj = jj - (ix(j) - 1)*ind_offset_pair(j)
                    end do

#ifdef __KCOMP
                    write (read_start + 89, '(32767f20.8)') &
#else
                        write (read_start + 89, '(1000000f22.8)') &
#endif
                        (rx(j) - 0.5d0*ell(vdim(j)), j=1, ddim) &
                        , (dvp*psip(i, k + shift), dvp*psip(i, k + shift + nind), i=maxf, 0, -1)

                end do

                shift = shift + npairind

                ! radial g
                vold = 0.d0

                do k = 1, ngrid_p

                    rv = (dble(k) - 0.5d0)/drmax
                    vnew = vnorm*dble(k)**ndim
                    prefac_g = pair_norm/(vnew - vold)
#ifdef __KCOMP
                    write (read_start + 90, '(32767f22.8)') rv &
                        , (prefac_g*psip(i, k + shift), prefac_g*psip(i, k + shift + nind), i=maxf, 0, -1)
#else
                    write (read_start + 90, '(100000f22.8)') rv &
                        , (prefac_g*psip(i, k + shift), prefac_g*psip(i, k + shift + nind), i=maxf, 0, -1)
#endif
                    vold = vnew

                end do

                shift = shift + ngrid_p

            else

                do k = 1, npairind - nshellsp

                    tmp_norm = 1.d0/mult(k)

                    if (allshells) then
                        j1 = (iioniond(k) - 1)/nion + 1
                        i1 = iioniond(k) - (j1 - 1)*nion
                        if (i1 .le. j1) then
#ifdef __KCOMP
                            write (read_start + 89, '(2I9,32767f22.8)') &
#else
                                write (read_start + 89, '(2I9,1000000f22.8)') &
#endif
                                ! <sigma_i sigma_j>, <sigma_i sigma_j> - <sigma_i> <sigma_j>
                                i1, j1, dists(k) &
                                , (tmp_norm*psip(i, k + shift), tmp_norm*psip(i, k + shift + nind), &
                                   tmp_norm*psip(i, k + shift) &
                                   - psip(i, shift - nrhoind - 1 + i1)*psip(i, shift - nrhoind - 1 + j1) &
                                   , tmp_norm*psip(i, k + shift + nind), i=maxf, 0, -1)
                        end if
                    else

#ifdef __KCOMP
                        write (read_start + 89, '(32767f20.8)') &
#else
                            write (read_start + 89, '(1000000f22.8)') &
#endif
                            ! <sigma_i sigma_j>
                            dists(k), (tmp_norm*psip(i, k + shift), tmp_norm*psip(i, k + shift + nind), &
                                       tmp_norm*psip(i, k + shift) - check(i), tmp_norm*psip(i, k + shift + nind), i=maxf, 0, -1)
                        ! <sigma_i sigma_j> - <sigma_i> <sigma_j>
                    end if

                end do

                sec = sec_spc(npairind - nshellsp + 1)
                write (read_start + 89, *) '  '
                write (read_start + 89, *) '  '
                write (read_start + 89, *) ' species-species sector:', sec

                do k = npairind - nshellsp + 1, npairind

                    tmp_norm = 1.d0/mult(k)

                    if (sec_spc(k) .ne. sec) then
                        sec = sec_spc(k)
                        write (read_start + 89, *) '  '
                        write (read_start + 89, *) '  '
                        write (read_start + 89, *) ' species-species sector:', sec
                    end if

#ifdef __KCOMP
                    write (read_start + 89, '(32767f20.8)') &
#else
                        write (read_start + 89, '(1000000f22.8)') &
#endif
                        ! <n_i n_j>
                        dists(k), (tmp_norm*psip(i, k + shift), tmp_norm*psip(i, k + shift + nind), i=maxf, 0, -1)

                end do

                shift = shift + npairind

            end if

        end if

        close (read_start + 88)
        if (ifpair) then
            close (read_start + 89)
            if (ngrid_p .ne. 0) close (read_start + 90)
        end if
    end subroutine write_spin

    subroutine write_sofk
        implicit none
        integer neldown

        if (decouple_files) then
            open (read_start + 75, file='skvecK'//trim(chara), form='formatted', status='unknown')
        else
            open (read_start + 75, file='skvec.dat', form='formatted', status='unknown')
        end if
        rewind (read_start + 75)
        if (fermi_flag) then
            if (decouple_files) then
                open (read_start + 76, file='skvec_ssK'//trim(chara), form='formatted', status='unknown')
                open (read_start + 77, file='skvec_uuK'//trim(chara), form='formatted', status='unknown')
                open (read_start + 78, file='skvec_ddK'//trim(chara), form='formatted', status='unknown')
                open (read_start + 79, file='skvec_udK'//trim(chara), form='formatted', status='unknown')
            else
                open (read_start + 76, file='skvec_ss.dat', form='formatted', status='unknown')
                open (read_start + 77, file='skvec_uu.dat', form='formatted', status='unknown')
                open (read_start + 78, file='skvec_dd.dat', form='formatted', status='unknown')
                open (read_start + 79, file='skvec_ud.dat', form='formatted', status='unknown')
            end if
            rewind (read_start + 76)
            rewind (read_start + 77)
            rewind (read_start + 78)
        end if

        neldown = nel - nelup

        antot = 1.d0/dble(nel)
        anorm(1) = 1.d0/dble(nel)
        anorm(2) = 1.d0/dble(nelup)
        anorm(3) = 1.d0/dble(neldown)
        anorm(4) = 1.d0/sqrt(dble(nelup*neldown))

        do i = 1, 3
            cellscale(i) = ell(vdim(i))*ncell(vdim(i))
        end do

        rx_old = 0.d0

        do k = 1, 2*nskind

            jj = abs(simap(k))
            signk = simap(k)/jj

            knorm_cry = 0.d0
            do l = 1, ddim
                knorm_cry = knorm_cry + rkcomp(l, jj)**2
            end do
            knorm_cry = sqrt(knorm_cry)

            if (yes_tilted) then
                vecscra(:) = 0.d0
                do j = 1, 3
                    vecscra(:) = vecscra(:) + recip(j, :)*signk*rkcomp(j, jj)*cellscale(j)/2.0/Pi
                end do
                knorm = 0.d0
                do l = 1, ddim
                    knorm = knorm + vecscra(l)**2
                end do
                knorm = sqrt(knorm)
            end if

            blank = .true.
            do l = 1, ddim
                if (dabs(signk*rkcomp(l, jj) - rx_old(l)) .lt. eps) blank = .false.
            end do

            if (blank .and. k .ne. 1 .and. ddim .ne. 1) write (read_start + 75, *)
#ifdef __KCOMP
            if (yes_tilted) then
                write (read_start + 75, '(32767f14.8)') (signk*rkcomp(l, jj), l=1, ddim), knorm_cry, vecscra(:), knorm &
                    , (antot*psip(i, jj + shift), antot*psip(i, jj + nind + shift), i=maxf, 0, -1)
            else
                write (read_start + 75, '(32767f14.8)') (signk*rkcomp(l, jj), l=1, ddim), knorm_cry &
                    , (antot*psip(i, jj + shift), antot*psip(i, jj + nind + shift), i=maxf, 0, -1)
            end if
#else
            if (yes_tilted) then
                write (read_start + 75, '(1000000f14.8)') (signk*rkcomp(l, jj), l=1, ddim), knorm_cry, vecscra(:), knorm &
                    , (antot*psip(i, jj + shift), antot*psip(i, jj + nind + shift), i=maxf, 0, -1)
            else
                write (read_start + 75, '(1000000f14.8)') (signk*rkcomp(l, jj), l=1, ddim), knorm_cry &
                    , (antot*psip(i, jj + shift), antot*psip(i, jj + nind + shift), i=maxf, 0, -1)
            end if
#endif
            if (fermi_flag) then
                do kk = 1, 4
                    if (blank .and. k .ne. 1 .and. ddim .ne. 1) write (read_start + 75 + kk, *)
#ifdef __KCOMP
                    if (yes_tilted) then
                        write (read_start + 75 + kk, '(32767f14.8)') (signk*rkcomp(l, jj), l=1, ddim) &
                            , knorm_cry, vecscra(:), knorm, &
                            (anorm(kk)*psip(i, jj + kk*nskind + shift) &
                             , anorm(kk)*psip(i, jj + nind + kk*nskind + shift), i=maxf, 0, -1)
                    else
                        write (read_start + 75 + kk, '(32767f14.8)') (signk*rkcomp(l, jj), l=1, ddim) &
                            , knorm_cry, &
                            (anorm(kk)*psip(i, jj + kk*nskind + shift) &
                             , anorm(kk)*psip(i, jj + nind + kk*nskind + shift), i=maxf, 0, -1)
                    end if
#else
                    if (yes_tilted) then
                        write (read_start + 75 + kk, '(1000000f14.8)') (signk*rkcomp(l, jj), l=1, ddim) &
                            , knorm_cry, vecscra(:), knorm, &
                            (anorm(kk)*psip(i, jj + kk*nskind + shift) &
                             , anorm(kk)*psip(i, jj + nind + kk*nskind + shift), i=maxf, 0, -1)
                    else
                        write (read_start + 75 + kk, '(1000000f14.8)') (signk*rkcomp(l, jj), l=1, ddim) &
                            , knorm_cry, &
                            (anorm(kk)*psip(i, jj + kk*nskind + shift) &
                             , anorm(kk)*psip(i, jj + nind + kk*nskind + shift), i=maxf, 0, -1)
                    end if
#endif
                end do
            end if

            do l = 1, ddim
                rx_old(l) = signk*rkcomp(l, jj)
            end do

        end do

        shift = shift + nskind
        if (fermi_flag) shift = shift + 4*nskind

        close (read_start + 75)
        if (fermi_flag) then
            close (read_start + 76)
            close (read_start + 77)
            close (read_start + 78)
            close (read_start + 79)
        end if
    end subroutine write_sofk

    subroutine write_corrs
        implicit none

        if (decouple_files) then
            open (read_start + 92, file='corrsamplingK'//trim(chara), form='formatted', status='unknown')
        else
            if (allow_write) &
                open (read_start + 92, file='corrsampling.dat', form='formatted', status='unknown')
        end if
        !  open(read_start+90,file='corrsampling_bootstrap.dat',form='formatted',status='unknown')
        if (allow_write) rewind (read_start + 92)

        !   write(read_start+89,'(a50,2(e,2x))') 'reweighting factor <(fort.10_corr/fort.10)^2>' &
        !   ,psip(0,shift+3),psip(0,nind+shift+3)
        !   write(read_start+89,'(a50,2(e,2x))') 'reference energy: E(fort.10)',psip(0,shift+1),psip(0,nind+shift+1)
        !   write(read_start+89,'(a50,2(e,2x))') 'reweighted energy: E(fort.10_corr)',psip(0,shift+5),psip(0,nind+shift+5)
        !   write(read_start+89,'(a50,2(e,2x))') 'reweighted difference: &
        !  &E(fort.10)-E(fort.10_corr) ',psip(0,shift+4),psip(0,nind+shift+4)

        shift = shift + 5

        call bootstrap(14, av, err, nbin)

        if (decouple_k .and. .not. decouple_files) then
#ifdef PARALLEL
            kmom = (rank + nproc/nk)/(nproc/nk)
            av(:) = av(:)*wkp(kmom)
            err(:) = (err(:)*wkp(kmom))**2
            call reduce_base_real(4, av, commcolrep_mpi, -1)
            call reduce_base_real(4, err, commcolrep_mpi, -1)
            err(1:4) = dsqrt(err(1:4))
#endif
        end if

        if (allow_write) then

            write (read_start + 92, '(a50,I18)') 'Number of bins ', nbin
            write (read_start + 92, '(a50,2(e18.9,2x))') 'reference energy: E(fort.10)', av(1), err(1)
            write (read_start + 92, '(a50,2(e18.9,2x))') 'reweighted energy: E(fort.10_corr)', av(2), err(2)
            write (read_start + 92, '(a50,2(e18.9,2x))') 'reweighted difference: &
                    &E(fort.10)-E(fort.10_corr) ', av(3), err(3)
            write (read_start + 92, '(a50,2(e18.9,2x))') ' Overlap square : &
                    &(fort.10,fort.10_corr) ', av(4), err(4)

            close (read_start + 92)
        end if
    end subroutine write_corrs

    subroutine write_rho_assar
        implicit none
        ! (Matteo) Calculating Assaraf density
        if (decouple_files) then
            open (read_start + 85, file='rho_assarK'//trim(chara), form='formatted', status='unknown')
        else
            open (read_start + 85, file='rho_assar.dat', form='formatted', status='unknown')
        end if
        rewind (read_start + 85)
        do k = 1, nrhoind
            write (read_start + 85, '(5e26.16)') (out_grid(i, k), i=1, 3) &
                , (psip(i, k + shift), psip(i, k + shift + nind), i=maxf, 0, -1)
        end do
        shift = shift + nrhoind
        close (read_start + 85)
    end subroutine write_rho_assar

    subroutine write_rho_corr
        implicit none
        ! (Matteo) Calculating Assaraf density
        if (decouple_files) then
            open (read_start + 85, file='rho_corrK'//trim(chara), form='formatted', status='unknown')
        else
            open (read_start + 85, file='rho_corr.dat', form='formatted', status='unknown')
        end if
        rewind (read_start + 85)
        do k = 1, nrhoind
            write (read_start + 85, '(5e26.16)') (out_grid(i, k), i=1, 3) &
                , (psip(i, k + shift), psip(i, k + shift + nind), i=maxf, 0, -1)
        end do
        shift = shift + nrhoind
        close (read_start + 85)
    end subroutine write_rho_corr

    subroutine write_spin2
        implicit none

        if (decouple_files) then
            open (read_start + 91, file='spinsquareK'//trim(chara), form='formatted', status='unknown')
        else
            open (read_start + 91, file='spinsquare.dat', form='formatted', status='unknown')
        end if
        rewind (read_start + 91)
        do i = maxf, 0, -1
#ifdef __KCOMP
            write (read_start + 91, '(I8,32767f14.8)') maxf - i, psip(i, 1 + shift), psip(i, 1 + shift + nind)
#else
            write (read_start + 91, '(I8,1000000f14.8)') maxf - i, psip(i, 1 + shift), psip(i, 1 + shift + nind)
#endif
        end do
        shift = shift + 1
        close (read_start + 91)
    end subroutine write_spin2

    subroutine write_qpwf
        implicit none

        if (ifqpwf_h) then
            dv = 1/(da(1)*da(2)*da(3))
        else
            dv = 1
        end if

        if (decouple_files) then
            open (read_start + 95, file='qpwfK'//trim(chara), form='formatted', status='unknown')
        else
            open (read_start + 95, file='qpwf.dat', form='formatted', status='unknown')
        end if
        rewind (read_start + 95)

        if (.not. ifqpwf_extr) then
            if (ipc .eq. 1 .and. .not. ifqpwf_k) then
                do k = 1, nrhoind
                    write (read_start + 95, '(5e26.16)') (out_grid(i, k), i=1, 3) &
                        , (psip(i, k + shift)*dv, psip(i, k + shift + nind)*dv, i=maxf, 0, -1)
                end do
            else
                do k = 1, nrhoind/2
                    write (read_start + 95, '(7e26.16)') (out_grid(i, k), i=1, 3), &
                        (psip(i, k + shift)/(2*pi)**(1.5), psip(i, k + shift + nind)/(2*pi)**(1.5), &
                         psip(i, k + shift + nrhoind/2)/(2*pi)**(1.5) &
                         , psip(i, k + shift + nind + nrhoind/2)/(2*pi)**(1.5), i=maxf, 0, -1)
                end do
            end if
        else
            if (ipc .eq. 1) then
                do k = 1, nrhoind
                    dv = 1/((da(1)/nrhoind)*k)**3
                    write (read_start + 95, '(3e26.12)') (da(1)/nrhoind)*k &
                        , (psip(i, k + shift)*dv, psip(i, k + shift + nind)*dv, i=maxf, 0, -1)
                end do
            else
                do k = 1, nrhoind
                    dv = 1/((da(1)/nrhoind)*k)**3
                    write (read_start + 95, '(5e26.12)') (da(1)/nrhoind)*k, &
                        (psip(i, k + 2*shift - 1)*dv, psip(i, k + 2*shift - 1 + nind)*dv, &
                         psip(i, k + 2*shift)*dv, psip(i, k + 2*shift + nind)*dv, i=maxf, 0, -1)
                end do
            end if
        end if
        shift = shift + ipc*nrhoind
        close (read_start + 95)

    end subroutine write_qpwf

    subroutine write_berry
        implicit none
        if (decouple_files) then
            open (read_start + 96, file='BerryK'//trim(chara), form='formatted', status='unknown')
        else
            open (read_start + 96, file='Berry.dat', form='formatted', status='unknown')
        end if
        rewind (read_start + 96)
        do i = maxf, 0, -1
            write (read_start + 96, '(I6,12f15.6)') maxf - i, (psip(i, k + shift), psip(i, k + shift + nind), k=1, 6)
        end do
        shift = shift + 6
    end subroutine write_berry

    subroutine write_dipole
        implicit none
        real(8) :: d_mod
        real(8) :: d_mod_err
        real(8) :: q_mod
        real(8) :: q_mod_err
        !    real(8), dimension(:,:), allocatable :: quad
        !    real(8), dimension(:), allocatable :: quad_diag
        integer shift_start

        d_mod = 0
        d_mod_err = 0
        q_mod = 0
        q_mod_err = 0

        if (decouple_files) then
            open (read_start + 93, file='dipoleK'//trim(chara), form='formatted', status='unknown')
        else
            open (read_start + 93, file='dipole.dat', form='formatted', status='unknown')
        end if
        rewind (read_start + 93)

        !    allocate(quad(nquad,nquad),quad_diag(nquad))
        ! quad and quad_diag already allocated in readforward

        shift_start = shift

        do i = maxf, 0, -1
            shift = shift_start
            d_mod = 0.d0
            d_mod_err = 0.d0
            do k = 1, ndipole
                d_mod = d_mod + psip(i, k + shift)**2
                d_mod_err = d_mod_err + psip(i, k + shift)**2*psip(i, k + shift + nind)**2
            end do
            write (read_start + 93, *) 'Dipole moment(a.u.)', maxf - i
            write (read_start + 93, '("|d| ",2f21.8)') sqrt(d_mod), sqrt(d_mod_err/d_mod)
            do l = 1, ndipole
                ll = vdim(l)
                shift = shift + 1
                write (read_start + 93, '(a3,2f21.8)') 'd_'//trim(cartesian(ll)), psip(i, shift), psip(i, shift + nind)
            end do
            write (read_start + 93, *) 'Nuclear contribution to Dipole moment(a.u.)'
            do l = 1, ndipole
                ll = vdim(l)
                write (read_start + 93, '(a3,2f21.8)') 'd_'//trim(cartesian(ll)), nuclear_dipole(l)
            end do
        end do

        shift_start = shift

        do i = maxf, 0, -1
            shift = shift_start
            !convert psip 1d array into 3*3 matrix
            !      quad=reshape(psip(i,shift+1:shift+nquad**2),[nquad,nquad])
            call dcopy(nquad**2, psip(i, shift + 1), (maxf + 1), quad, 1)
            call diag_quad
            q_mod = 0.d0
            q_mod_err = 0.d0
            do k = 1, nquad
                q_mod = q_mod + quad_diag(k)**2
                q_mod_err = q_mod_err + quad_diag(k)**2*psip(i, k*(nquad + 1) - nquad + shift + nind)**2
            end do

            write (read_start + 93, *) 'Effective quadrupole moment q (a.u.) =  sqrt[2/3(Qxx^2 + Qyy^2 + Qzz^2)]', maxf - i
            write (read_start + 93, '("|q|  ",2f21.8)') sqrt(2.d0/3.d0*q_mod), sqrt(q_mod_err/q_mod)
            do l = 1, nquad
                ll = vdim(l)
                do m = 1, nquad
                    mm = vdim(m)
                    shift = shift + 1

                    write (read_start + 93, '(a4,2f21.8)') 'q_'//trim(cartesian(mm))//trim(cartesian(ll)), &
                        psip(i, shift), psip(i, shift + nind)
                    !      write(read_start+93,'("q_xy ",2f14.8)') psip(i,5+shift),psip(i,5+shift+nind)
                    !      write(read_start+93,'("q_xz ",2f14.8)') psip(i,6+shift),psip(i,6+shift+nind)
                    !      write(read_start+93,'("q_yx ",2f14.8)') psip(i,7+shift),psip(i,7+shift+nind)
                    !      write(read_start+93,'("q_yy ",2f14.8)') psip(i,8+shift),psip(i,8+shift+nind)
                    !      write(read_start+93,'("q_yz ",2f14.8)') psip(i,9+shift),psip(i,9+shift+nind)
                    !      write(read_start+93,'("q_zx ",2f14.8)') psip(i,10+shift),psip(i,10+shift+nind)
                    !      write(read_start+93,'("q_zy ",2f14.8)') psip(i,11+shift),psip(i,11+shift+nind)
                    !      write(read_start+93,'("q_zz ",2f14.8)') psip(i,12+shift),psip(i,12+shift+nind)
                    !      write(read_start+93,'("diag_xx",2f14.8)') quad_diag_c(1)
                    !      write(read_start+93,'("diag_yy",2f14.8)') quad_diag_c(2)
                    !      write(read_start+93,'("diag_zz",2f14.8)') quad_diag_c(3)

                end do
            end do
        end do

        !    deallocate(quad,quad_diag)
        close (read_start + 93)

    end subroutine write_dipole

end subroutine write_corr_fun

subroutine bootstrap(unit, av, err, nbin)
    use constants, only: ipc
    use allio, only: rank
    implicit none

    integer unit, i, nbin, kmain, nmis, k, j
    real(8) av(4), err(4)
    real(8), dimension(:), allocatable :: w, ws, e, es, o, oi
    real(8) :: wt, wts, et, ets, ot, oti, ots, diff, drand1

    rewind (unit)
    nbin = 0
    do while (nbin .ge. 0)
        read (unit, *, end=500)
        nbin = nbin + 1
    end do
500 continue

    if (rank .eq. 0) write (6, *) ' number of bins read in file corrsampling_bin.dat =', nbin

    allocate (e(nbin), es(nbin), w(nbin), ws(nbin), o(nbin))
    if (ipc .eq. 2) allocate (oi(nbin))

    rewind (unit)
    if (ipc .eq. 1) then
        do i = 1, nbin
            read (unit, *) e(i), w(i), es(i), ws(i), o(i)
        end do
    else
        do i = 1, nbin
            read (unit, *) e(i), w(i), es(i), ws(i), o(i), oi(i)
        end do
    end if

    nmis = 1000

    do i = 1, 4
        av(i) = 0.d0
        err(i) = 0.d0
    end do

    do kmain = 1, nmis

        wt = 0.d0
        et = 0.d0
        ets = 0.d0
        wts = 0.d0
        ot = 0.d0
        oti = 0.d0
        ! bootstrap estimate sample
        do k = 1, nbin
            j = drand1()*nbin + 1
            !       j=k
            if (w(j) .ne. 0.d0 .and. ws(j) .ne. 0.d0) then ! e,es may be NaN
                wt = wt + w(j)
                wts = wts + ws(j)
                et = et + e(j)*w(j)
                ot = ot + o(j)*w(j)
                if (ipc .eq. 2) oti = oti + oi(j)*w(j)
                ets = ets + es(j)*ws(j)
            end if
        end do

        et = et/wt
        ets = ets/wts
        ot = ot/wt
        oti = oti/wt
        ots = wts/wt

        ! now average correlation function

        diff = et - ets

        av(3) = av(3) + diff
        err(3) = err(3) + diff**2

        av(1) = av(1) + et
        err(1) = err(1) + et**2

        av(2) = av(2) + ets
        err(2) = err(2) + ets**2

        diff = (ot**2 + oti**2)/ots

        av(4) = av(4) + diff
        err(4) = err(4) + diff**2

    end do

    err = err/nmis
    av = av/nmis

    err(:) = dsqrt(max(err(:) - av(:)**2, 1.d-15))

    ! replace the averages with the true averages without bootstrap

    av(1) = sum(e(:)*w(:))/sum(w(:))
    av(2) = sum(es(:)*ws(:))/sum(ws(:))
    av(3) = av(1) - av(2)

    if (ipc .eq. 1) then
        av(4) = (sum(o(:)*w(:))/sum(w(:)))**2/(sum(ws(:))/sum(w(:)))
    else
        av(4) = (sum(o(:)*w(:))**2 + sum(oi(:)*w(:))**2)/sum(w(:))**2/(sum(ws(:))/sum(w(:)))
    end if
    deallocate (e, es, w, ws, o)

end subroutine bootstrap
