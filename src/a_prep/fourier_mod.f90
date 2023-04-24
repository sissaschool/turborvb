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

module fourier_module

    !
    ! This modules contains all variables and procedures
    ! related to the DFT transform in the code.
    ! The subroutines contained are:
    ! evalvhartreeq = evaluation of the Hartree potential (part independent from the density)
    ! update_vhartree = density-dependent part of the Hartree potential (called in uphamilt_new)
    ! precond = preconditioning of the density with mixingegder (tfcut option in the nml DFT)
    ! evalgrad = compute the gradient of the density in Fourier space and back.
    !            Useful for PBE functional, NOT USED in the present version.
    !

    use allio, only: nprocrep, commrep_mpi, rankrep, rank, iespbc, &
                     nx, ny, nz, ax, ay, az, rion, nion, ipc, zetar, norm_metric
    use Ewald, only: kappa
    use constants, only: pi, zzero
    use fft_scalar, only: cfft3d, cft_1z, cft_2xy
    use setup, only: double_mesh, l0_at, scale_z, nx0, ny0, nz0, weightx, weighty&
            &, weightz, meshproc, linear, sub_comm_fft, nproc_fft, bufbuf, nelorb3, mesh&
            &, do_hartree, vmax0, from_ions, minz_at, meshproc_tot, commensurate_lattice
    use cell, only: metric, yes_tilted, recip, at, ApplyPBC, dist_shift, x_neigh, neigh
    implicit none

    public
    ! Hartree potentials
    integer, dimension(:, :), allocatable :: ncub_min, ncub_max
    integer, dimension(:), allocatable :: ind_init
    integer*8, dimension(:), allocatable :: nx8, nxny8
    integer, dimension(:), allocatable :: nx_proc, nxny_proc
    real(8), dimension(:), allocatable :: vhartree, vhartreeq
    real(8), dimension(:, :, :), allocatable :: vhartree_global
    real*8, dimension(:, :), allocatable :: buffer_fft
    real*8, dimension(:), allocatable :: buffer_local
    complex(8), dimension(:), allocatable :: fp, sndbuf, rcvbuf
    ! local variables
    integer :: dim_fp1, dim_fp2, nxl, nzl, nxl2, nzl2, nxu, nyu, nzu, dim_fft_vh, dim_fft_den
    real(8) :: fx, fy, fz, axu, ayu, azu

contains

    subroutine upload_fp(indmesh, nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, add)
        use setup, only: dent
        implicit none
        integer nxi, nxf, nyi, nyf, nzi, nzf, i, j, k, indproc, indtot, ind, nbufrep, ierr, indmesh&
                &, ii, jj, indi, indj, indk, indi_proc, indj_proc, indk_proc
        real*8 ax, ay, az
        logical add, doextra
        integer, external :: mod_true
#ifdef PARALLEL
        include 'mpif.h'
#endif
        indtot = 0
        indproc = 0
        ind = 0

        if (add) then
            allocate (buffer_fft(dim_fft_den, nprocrep), buffer_local(dim_fft_den))
            buffer_fft = 0.d0
            buffer_local = 0.d0
            nbufrep = nprocrep*dim_fft_den
            indk = nzi
            indj = nyi
            indi = nxi - 1
            do k = nzi, nzf
                do j = nyi, nyf
                    do i = nxi, nxf
                        indtot = indtot + 1
                        if (rankrep .eq. indproc) then
                            indmesh = indmesh + 1
                            ind = ind + 1
                            buffer_local(ind) = dent(indmesh)
                        end if
                        if (indtot .eq. nbufrep .or. (i .eq. nxf .and. j .eq. nyf .and. k .eq. nzf)) then
#ifdef PARALLEL
                            call mpi_allgather(buffer_local, dim_fft_den, MPI_DOUBLE_PRECISION, buffer_fft&
                               &, dim_fft_den, MPI_DOUBLE_PRECISION, commrep_mpi, ierr)
#else
                            buffer_fft(:, 1) = buffer_local(:)
#endif

                            do ii = 1, dim_fft_den
                                do jj = 1, nprocrep
                                    ! Update counter i,j,k in the right sequence
                                    indi = indi + 1
                                    if (indi .gt. nxf) then
                                        indi = nxi
                                        indj = indj + 1
                                        if (indj .gt. nyf) then
                                            indj = nyi
                                            indk = indk + 1
                                        end if
                                    end if
                                    if (sub_comm_fft%yesin) then
                                        if (mod_true(indi - 1, nproc_fft) .eq. sub_comm_fft%rank .and. &
                                            indi .le. nxf .and. indj .le. nyf .and. indk .le. nzf) then
                                            if (iespbc) then
                                                indi_proc = mod_true(indi - 1, nxu) + 1
                                                indj_proc = mod_true(indj - 1, nyu) + 1
                                                indk_proc = mod_true(indk - 1, nzu) + 1
                                            else
                                                indi_proc = indi
                                                indj_proc = indj
                                                indk_proc = indk
                                                if (indi .lt. 1) indi = 1
                                                if (indi .gt. 2*nxu) indi = 2*nxu
                                                if (indj .lt. 1) indj = 1
                                                if (indj .gt. 2*nyu) indj = 2*nyu
                                                if (indk .lt. 1) indk = 1
                                                if (indk .gt. 2*nzu) indk = 2*nzu
                                            end if
                                            indi_proc = (indi_proc - 1)/nproc_fft + 1
                                            if (iespbc) then
                                                fp(indk_proc + nzu*(indi_proc - 1) + nzu*nxl*(indj_proc - 1)) &
                                                    = buffer_fft(ii, jj)
                                            else
                                                fp(indk_proc + 2*nzu*(indi_proc - 1) + 2*nzu*nxl2*(indj_proc - 1)) &
                                                    = buffer_fft(ii, jj)
                                            end if
                                        end if ! sub_comm
                                    end if ! mod_true
                                end do
                            end do
                            buffer_local = 0.d0
                            ind = 0
                            indtot = 0
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                    end do
                end do
            end do
            deallocate (buffer_fft, buffer_local)
        else
            ! do nothing just update the counter indmesh
            do k = nzi, nzf
                do j = nyi, nyf
                    do i = nxi, nxf
                        if (rankrep .eq. indproc) indmesh = indmesh + 1
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                    end do
                end do
            end do
        end if
        return
    end subroutine upload_fp

    subroutine initialize_fourier
        use constants, only: zzero
        use setup, only: do_hartree, double_mesh, scale_z
        implicit none
        integer scale_here
        if (do_hartree .and. double_mesh .and. scale_z .ne. 1) then
            nxl = (nx*scale_z)/nproc_fft
            if (nxl*nproc_fft .ne. nx*scale_z .and. rankrep .lt. nx*scale_z - nproc_fft*nxl) nxl = nxl + 1
            nzl = (nz*scale_z)/nproc_fft
            if (nzl*nproc_fft .ne. nz*scale_z .and. rankrep .lt. nz*scale_z - nproc_fft*nzl) nzl = nzl + 1
            nxl2 = (2*nx*scale_z)/nproc_fft
            nzl2 = (2*nz*scale_z)/nproc_fft
            scale_here = scale_z
        else
            nxl = nx/nproc_fft
            if (nxl*nproc_fft .ne. nx .and. rankrep .lt. nx - nproc_fft*nxl) nxl = nxl + 1
            nzl = nz/nproc_fft
            if (nzl*nproc_fft .ne. nz .and. rankrep .lt. nz - nproc_fft*nzl) nzl = nzl + 1
            nxl2 = (2*nx)/nproc_fft
            nzl2 = (2*nz)/nproc_fft
            scale_here = 1
        end if
        nxu = nx*scale_here
        nyu = ny*scale_here
        nzu = nz*scale_here
        axu = ax/scale_here
        ayu = ay/scale_here
        azu = az/scale_here
        ! dimensions of Fourier transforms
        if (iespbc) then
            dim_fp1 = nxl*nyu*nzu
            dim_fp2 = nxu*nyu*nzl
        else
            dim_fp1 = 4*nxl2*nyu*nzu
            dim_fp2 = 4*nxu*nyu*nzl2
        end if
        !
        allocate (vhartree(meshproc))
        vhartree = 0.d0
        allocate (fp(dim_fp1), vhartreeq(dim_fp2))
        fp = 0.d0
        vhartreeq = 0.d0
#ifdef PARALLEL
        allocate (sndbuf(dim_fp1), rcvbuf(dim_fp2))
        sndbuf = zzero
        rcvbuf = zzero
#endif
        !  allocate with minimal dimension to avoid overallocation.
        !  It should be not too small otherwise latency time will be important
        dim_fft_vh = max(nint(dble(bufbuf)/nproc_fft), 128)
        dim_fft_den = max(nint(dble(bufbuf)/nprocrep), 128)

        if (double_mesh) then
            if (from_ions) then
                allocate (ncub_min(3, 2*nion), ncub_max(3, 2*nion))
                allocate (nx8(2*nion), nxny8(2*nion), nx_proc(2*nion), nxny_proc(2*nion))
                allocate (ind_init(2*nion))
            else
                allocate (nx8(2), nxny8(2), nx_proc(2), nxny_proc(2))
                allocate (ncub_min(3, 2), ncub_max(3, 2))
                allocate (ind_init(2))
            end if
            ind_init = 0
            nx8 = 0
            nxny8 = 0
            nx_proc = 0
            nxny_proc = 0
            ncub_min = 0
            ncub_max = 0
        end if
        return
    end subroutine initialize_fourier

    subroutine evalvhartreeq

        use setup, only: voltot, voltot_double, vmax0, volmesh, double_mesh, time_fft, time_uploadfft, vmax0_in
        use dielectric
        implicit none

        real(8) :: q, q2, vecscra(3), qx, qy, qz, distel, vgauss, dist_kel(3), timep
        integer :: i, j, k, iq, ixl, izl, ii, kk, ip, jj, ll, indmesh, ierr
        real(8), external :: cclock
        real(8) derfc

#ifdef PARALLEL
        include 'mpif.h'
#endif

        vgauss = 0.d0
        vmax0 = 0.d0
        timep = cclock()

#ifdef PARALLEL

        fp = zzero
        sndbuf = zzero
        rcvbuf = zzero
        vhartreeq = 0.d0

        if (sub_comm_fft%yesin) then

            if (iespbc) then
                do i = 1, nxl
                    ii = sub_comm_fft%rank + 1 + (i - 1)*nproc_fft
!       dist_kel(:)=((ii-1)-nxu/2)*axu
                    do j = 1, nyu
!          dist_kel(2)=((j-1)-nyu/2)*ayu
                        do k = 1, nzu
!             dist_kel(3)=((k-1)-nzu/2)*azu
                            dist_kel(:) = ((ii - 1) - nxu/2)*axu*at(:, 1) + ((j - 1) - nyu/2)*ayu*at(:, 2) +&
                        & ((k - 1) - nzu/2)*azu*at(:, 3)
                            call ApplyPBC(dist_kel, 1)
!             distel=sum(dist_kel(1:3)**2)

#ifdef _SIMD
!$omp simd
#endif
                            do ll = 1, neigh
                                dist_shift(ll) = dsqrt((dist_kel(1) + x_neigh(ll, 1))**2 &
                                                       + (dist_kel(2) + x_neigh(ll, 2))**2 &
                                                       + (dist_kel(3) + x_neigh(ll, 3))**2)
                                if (dist_shift(ll) .gt. 0.d0) dist_shift(ll) = rep_erfc(dist_shift(ll), kappa)
                            end do
                            vmax0 = sum(dist_shift(1:neigh))
                            fp(k + nzu*(i - 1) + nzu*nxl*(j - 1)) = dcmplx(vmax0)
                        end do
                    end do
                end do

                call cft_1z(fp, nxl*nyu, nzu, nzu, -1) ! Fourier transform of the Hartree potential

                do i = 1, nxl
                    do j = 1, nyu
                        do ip = 1, nproc_fft
                            do k = 1 + nzl*(ip - 1), nzl*ip
                                sndbuf(i + nxl*(j - 1) + nxl*nyu*(k - 1)) = fp(k + nzu*(i - 1) + nzu*nxl*(j - 1))
                            end do
                        end do
                    end do
                end do
                call MPI_ALLTOALL(sndbuf, nxl*nyu*nzl, MPI_DOUBLE_COMPLEX, &
                                  rcvbuf, nxl*nyu*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
                do ip = 1, nproc_fft
                    do j = 1, nyu
                        do k = 1, nzl
                            jj = nxu*(j - 1) + nxu*nyu*(k - 1)
                            kk = nxl*(j - 1) + nxl*nyu*(k - 1 + (ip - 1)*nzl)
                            do i = 1, nxl
                                fp(ip + nproc_fft*(i - 1) + jj) = rcvbuf(i + kk)
                            end do
                        end do
                    end do
                end do

                call cft_2xy(fp, nzl, nxu, nyu, nxu, nyu, -1)

                vhartreeq = abs(fp)*voltot
                if (vmax0_in .eq. -100.d0) then
!     if(corr_hartree) then
!     vmax0=0.d0
!     else
                    vmax0 = -vhartreeq(1) ! q=0 component of the Ewald potential
                    call mpi_bcast(vmax0, 1, MPI_DOUBLE_PRECISION, 0, sub_comm_fft%comm, ierr)
                    vmax0 = vmax0 + vq0_diel
!     endif
                else
                    vmax0 = vmax0_in*ax**4
                end if
                ! shift to have the correct q=0 limit
                vhartreeq(:) = vhartreeq(:) + vmax0

                do i = 1, nxu
                    if (i > (nxu + 1)/2) then
                        qx = fx*(i - nxu - 1)
                    else
                        qx = fx*(i - 1)
                    end if
                    do j = 1, nyu
                        if (j > (nyu + 1)/2) then
                            qy = fy*(j - nyu - 1)
                        else
                            qy = fy*(j - 1)
                        end if
                        do kk = 1, nzl
                            k = kk + nzl*sub_comm_fft%rank
                            if (k > (nzu + 1)/2) then
                                qz = fz*(k - nzu - 1)
                            else
                                qz = fz*(k - 1)
                            end if
                            if (yes_tilted) then
                                vecscra(:) = recip(1, :)*qx + recip(2, :)*qy + recip(3, :)*qz
                                q2 = sum(vecscra(:)**2)
                            else
                                q2 = qx*qx + qy*qy + qz*qz
                            end if
                            iq = i + nxu*(j - 1) + nxu*nyu*(kk - 1)
                            if (q2 .gt. 0.d0) then
!             if(corr_hartree) then
!             vhartreeq(iq)=4.d0*pi/q2+vmax0 ! The simplest one
!             else
                                vhartreeq(iq) = vhartreeq(iq) + 4.d0*epsilon0*pi/q2*dexp(-0.25d0*q2/kappa**2)
!             endif
                            else
                                ! Already included in the single particle eigenvalues
                                ! to get them independent of kappa
                                vhartreeq(iq) = 0.d0
                            end if
                        end do
                    end do
                end do

            else ! open systems case

                indmesh = 0
                do k = 1, 2*nzu
                    dist_kel(3) = ((k - 1) - nzu)*azu
                    do j = 1, 2*nyu
                        dist_kel(2) = ((j - 1) - nyu)*ayu
                        do i = sub_comm_fft%rank + 1, 2*nxu, nproc_fft
                            indmesh = indmesh + 1
                            dist_kel(1) = ((i - 1) - nxu)*axu
                            distel = sum(dist_kel(1:3)**2)
                            if (distel .gt. 0.d0) then
                                distel = dsqrt(distel)
                                sndbuf(indmesh) = veps(distel)
                                vgauss = vgauss + volmesh*dexp(-distel**2/2.d0)*veps(distel)
                            else
                                sndbuf(indmesh) = dcmplx(vmax0)
                            end if
                        end do
                    end do
                end do
                if (vmax0_in .eq. -100.d0) then
                    vmax0 = vgauss
                    call mpi_allreduce(vmax0, vgauss, 1, MPI_DOUBLE_PRECISION, MPI_SUM, sub_comm_fft%comm, ierr)

                end if

                do i = 1, nxl2
                    do j = 1, 2*nyu
                        do k = 1, 2*nzu
                            fp(k + 2*nzu*(i - 1) + 2*nzu*nxl2*(j - 1)) = sndbuf(i + nxl2*(j - 1) + 2*nxl2*nyu*(k - 1))
                        end do
                    end do
                end do

                call cft_1z(fp, 2*nxl2*nyu, 2*nzu, 2*nzu, -1)

                do i = 1, nxl2
                    do j = 1, 2*nyu
                        do ip = 1, nproc_fft
                            do k = 1 + nzl2*(ip - 1), nzl2*ip
                                sndbuf(i + nxl2*(j - 1) + 2*nxl2*nyu*(k - 1)) = fp(k + 2*nzu*(i - 1) + 2*nzu*nxl2*(j - 1))
                            end do
                        end do
                    end do
                end do
                call MPI_ALLTOALL(sndbuf, 2*nxl2*nyu*nzl2, MPI_DOUBLE_COMPLEX, &
                                  rcvbuf, 2*nxl2*nyu*nzl2, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
                do ip = 1, nproc_fft
                    do j = 1, 2*nyu
                        do k = 1, nzl2
                            jj = 2*nxu*(j - 1) + 4*nxu*nyu*(k - 1)
                            kk = nxl2*(j - 1) + 2*nxl2*nyu*(k - 1 + (ip - 1)*nzl2)
                            do i = 1, nxl2
                                fp(ip + nproc_fft*(i - 1) + jj) = rcvbuf(i + kk)
                            end do
                        end do
                    end do
                end do

                call cft_2xy(fp, nzl2, 2*nxu, 2*nyu, 2*nxu, 2*nyu, -1)

                vhartreeq = 8.d0*abs(fp)*voltot

                !  We use that the integral of 1/r exp(-r^2/2)= 4 pi in D=3
                !  So in order to have this result exactly one has to add a delta function
                !  contribution for r=0 equal to the relation below.
                if (vmax0_in .eq. -100.d0) then
!     if(corr_hartree) then
!     vmax0=0.d0
!     else
                    vmax0 = vgauss_diel - vgauss
!     endif
                else
                    vmax0 = vmax0_in*ax**4
                end if
                vhartreeq(:) = vhartreeq(:) + vmax0

            end if
        end if
        if (nprocrep .gt. nproc_fft) call bcast_real(vmax0, 1, 0, commrep_mpi)

#else

        ! scalar algorithm
        if (iespbc) then

            indmesh = 0
            do k = 1, nzu
                !       dist_kel(3)=((k-1)-nzu/2)*azu
                do j = 1, nyu
                    !          dist_kel(2)=((j-1)-nyu/2)*ayu
                    do i = 1, nxu
                        !             dist_kel(1)=((i-1)-nxu/2)*axu
                        !             if(yes_tilted) then
                        !             distel=norm_metric(dist_kel,metric)**2
                        !             else
                        dist_kel(:) = ((i - 1) - nxu/2)*axu*at(:, 1) + ((j - 1) - nyu/2)*ayu*at(:, 2) + &
                                & ((k - 1) - nzu/2)*azu*at(:, 3)
                        call ApplyPBC(dist_kel, 1)
!                       distel = sum(dist_kel(1:3)**2)
                        !             endif
                        indmesh = indmesh + 1
!                       if(distel.gt.0.d0) then
!                           distel = dsqrt(distel)
!                           fp(indmesh) = derfc(kappa * distel) / distel
!$omp simd
                        do ll = 1, neigh
                            dist_shift(ll) = dsqrt((dist_kel(1) + x_neigh(ll, 1))**2 &
                                                   + (dist_kel(2) + x_neigh(ll, 2))**2 &
                                                   + (dist_kel(3) + x_neigh(ll, 3))**2)
                            if (dist_shift(ll) .gt. 0.d0) dist_shift(ll) = rep_erfc(dist_shift(ll), kappa)
                        end do
!             vmax0=sum(dist_shift(1:neigh))
                        fp(indmesh) = sum(dist_shift(1:neigh))
                    end do
                end do
            end do
            call cfft3d(fp, nxu, nyu, nzu, nxu, nyu, nzu, -1)
            vhartreeq = abs(fp)*voltot

            if (vmax0_in .eq. -100.d0) then
                !     if(corr_hartree) then
                !     vmax0=0.d0
                !     else
                vmax0 = -vhartreeq(1)
                vmax0 = vmax0 + vq0_diel
!               vmax0 = vmax0 + pi / kappa**2
                !     endif
            else
                vmax0 = vmax0_in*ax**4
            end if
            vhartreeq(:) = vhartreeq(:) + vmax0

            do i = 1, nxu
                if (i > (nxu + 1)/2) then
                    qx = fx*(i - nxu - 1)
                else
                    qx = fx*(i - 1)
                end if

                do j = 1, nyu
                    if (j > (nyu + 1)/2) then
                        qy = fy*(j - nyu - 1)
                    else
                        qy = fy*(j - 1)
                    end if
                    do k = 1, nzu
                        if (k > (nzu + 1)/2) then
                            qz = fz*(k - nzu - 1)
                        else
                            qz = fz*(k - 1)
                        end if
                        if (yes_tilted) then
                            vecscra(:) = recip(1, :)*qx + recip(2, :)*qy + recip(3, :)*qz
                            q2 = sum(vecscra(:)**2)
                        else
                            q2 = qx*qx + qy*qy + qz*qz
                        end if
                        iq = i + nxu*(j - 1) + nxu*nyu*(k - 1)
                        if (q2 .gt. 0.d0) then
                            !        if(corr_hartree) then
                            !        vhartreeq(iq) = 4.d0*pi/q2+vmax0
                            !        else
                            vhartreeq(iq) = vhartreeq(iq) + 4.d0*epsilon0*pi/q2*dexp(-0.25d0*q2/kappa**2)
                            !        endif
                        else
                            !          Already included in the single particle eigenvalues
                            !          to get them independent of kappa
                            vhartreeq(iq) = 0.d0
                        end if
                    end do
                end do
            end do

        else
            !
            !  We use that the integral of 1/r exp(-r^2/2)= 4 pi in D=3
            !  So in order to have this result exactly one has to add a delta function
            !  contribution for r=0 equal to the relation below.
            !
            indmesh = 0
            do k = 1, 2*nzu
                dist_kel(3) = ((k - 1) - nzu)*azu
                do j = 1, 2*nyu
                    dist_kel(2) = ((j - 1) - nyu)*ayu
                    do i = 1, 2*nxu
                        dist_kel(1) = ((i - 1) - nxu)*axu
                        distel = sum(dist_kel(1:3)**2)
                        indmesh = indmesh + 1
                        if (distel .gt. 0.d0) then
                            distel = dsqrt(distel)
                            fp(indmesh) = veps(distel)
                            vgauss = vgauss + volmesh*dexp(-distel**2/2.d0)*veps(distel)
                        else
                            fp(indmesh) = vmax0
                        end if
                    end do
                end do
            end do

            call cfft3d(fp, 2*nxu, 2*nyu, 2*nzu, 2*nxu, 2*nyu, 2*nzu, -1)

            vhartreeq(:) = abs(fp)*voltot*8.d0

            if (vmax0_in .eq. -100.d0) then
                !     if(corr_hartree) then
                !     vmax0=0.d0
                !     else
                vmax0 = vgauss_diel - vgauss
!               vmax0 = 4.d0 * pi - vgauss
                !     endif
            else
                vmax0 = vmax0_in*ax**4
            end if

            vhartreeq(:) = vhartreeq(:) + vmax0

        end if

#endif

        if (rank .eq. 0) write (6, *) ' Correction Coulomb =', vmax0

        time_fft = time_fft + cclock() - timep

        return

    end subroutine evalvhartreeq

    subroutine update_vhartree

        use allio, only: zetar
        use setup, only: dent, voltot, rion_shift, rion_ref, eh_ew, vh_test&
                &, vh_att, corr_hartree, meshproc, meshproc_tot, minz_at, do_hartree&
                &, from_ions, rion_from, nx_at, ny_at, nz_at, time_fft, dens0, scale_hartree &
                &, scale_hartreen, weightvh, rion_from, time_uploadfft, volmesh_proc
        use constants, only: TWO_PI, PI
        implicit none
        real(8) :: q, q2, qx, qy, qz, scal, rion_upload(3), axn, ayn, azn, timep, costvhq &
                   , cost_att, timepp, vshould, vhq0, denstot, cost, vecscra(3)
        integer :: i, j, k, ii, jj, kk, iq, ixl, izl, ip, indfp, indmesh, ierr&
                &, nx0n, ny0n, nz0n, nxr, nyr, nzr, nxi, nxf, nyi, nyf, nzi, nzf, scalea
        complex(8) :: ccost
        real*8, external :: attenuate, cclock
#ifdef PARALLEL
        include 'mpif.h'
#endif
        !
        ! computation of Hartree potential in real space
        !
        timep = cclock()

        if (weightvh .eq. 0.d0) return
#ifdef PARALLEL

        vh_test = 0.d0
        vh_att = 0.d0
        eh_ew = 0.d0
        fp = zzero
        sndbuf = zzero
        rcvbuf = zzero
        vhartree = 0.d0
        if (.not. commensurate_lattice) call upload_fp_fromdent

!   cost=sum(abs(dent(:)))
!   call reduce_base_real(1,cost,commrep_mpi,-1)
!   if(rank.eq.0)  write(6,*) ' Input dent =',cost
!   cost=sum(abs(vhartreeq(:)))
!   call reduce_base_real(1,cost,commrep_mpi,-1)
!   if(rank.eq.0)  write(6,*) ' Input vhartreeq =',cost

        if (sub_comm_fft%yesin) then

            if (iespbc) then
                if (commensurate_lattice) then
                do i = 1, nxl
                    do j = 1, nyu
                        do k = 1, nzu
                            fp(k + nzu*(i - 1) + nzu*nxl*(j - 1)) = dent(i + nxl*(j - 1) + nxl*nyu*(k - 1))
                        end do
                    end do
                end do
                end if
                time_uploadfft = time_uploadfft + cclock() - timep
                call cft_1z(fp, nxl*nyu, nzu, nzu, -1)

                do i = 1, nxl
                    do j = 1, nyu
                        do ip = 1, nproc_fft
                            do k = 1 + nzl*(ip - 1), nzl*ip
                                sndbuf(i + nxl*(j - 1) + nxl*nyu*(k - 1)) = fp(k + nzu*(i - 1) + nzu*nxl*(j - 1))
                            end do
                        end do
                    end do
                end do

                call MPI_ALLTOALL(sndbuf, nxl*nyu*nzl, MPI_DOUBLE_COMPLEX, &
                                  rcvbuf, nxl*nyu*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)

                do ip = 1, nproc_fft
                    do j = 1, nyu
                        do k = 1, nzl
                            jj = nxu*(j - 1) + nxu*nyu*(k - 1)
                            kk = nxl*(j - 1) + nxl*nyu*(k - 1 + (ip - 1)*nzl)
                            do i = 1, nxl
                                fp(ip + nproc_fft*(i - 1) + jj) = rcvbuf(i + kk)
                            end do
                        end do
                    end do
                end do

                call cft_2xy(fp, nzl, nxu, nyu, nxu, nyu, -1)

                do i = 1, nxu
                    if (i > (nxu + 1)/2) then
                        qx = fx*(i - nxu - 1)
                    else
                        qx = fx*(i - 1)
                    end if
                    do j = 1, nyu
                        if (j > (nyu + 1)/2) then
                            qy = fy*(j - nyu - 1)
                        else
                            qy = fy*(j - 1)
                        end if
                        do kk = 1, nzl
                            k = kk + nzl*sub_comm_fft%rank
                            if (k > (nzu + 1)/2) then
                                qz = fz*(k - nzu - 1)
                            else
                                qz = fz*(k - 1)
                            end if

                            if (yes_tilted) then
                                vecscra(:) = recip(1, :)*qx + recip(2, :)*qy + recip(3, :)*qz
                                q2 = sum(vecscra(:)**2)
                            else
                                q2 = qx*qx + qy*qy + qz*qz
                            end if

                            iq = i + nxu*(j - 1) + nxu*nyu*(kk - 1)

                            cost_att = 1.d0
                            if (scale_hartree .gt. 0.d0) cost_att = cost_att*scale_hartreen

                            if (q2 .gt. 0) then

                                vh_test = vh_test + 0.5d0*vhartreeq(iq)*fp(iq)*dconjg(fp(iq))*voltot*cost_att**2
                                eh_ew = eh_ew + TWO_PI/q2*fp(iq)*dconjg(fp(iq))*voltot*cost_att**2
                                vh_att = vh_att + TWO_PI/((2.d0 - 2.d0*cos(qx*axu))/axu**2 &
                              &+ (2 - 2.d0*cos(qy*ayu))/ayu**2 + (2.d0 - 2.d0*cos(qz*azu))/azu**2)*fp(iq)*dconjg(fp(iq))*voltot
                                fp(iq) = fp(iq)*vhartreeq(iq)*cost_att

                            else
                                fp(iq) = (0.d0, 0.d0)
                            end if
                        end do
                    end do
                end do

                call cft_2xy(fp, nzl, nxu, nyu, nxu, nyu, 1)

                do ip = 1, nproc_fft
                    do j = 1, nyu
                        do k = 1, nzl
                            jj = nxu*(j - 1) + nxu*nyu*(k - 1)
                            kk = nxl*(j - 1) + nxl*nyu*(k - 1 + (ip - 1)*nzl)
                            do i = 1, nxl
                                rcvbuf(i + kk) = fp((ip + nproc_fft*(i - 1)) + jj)
                            end do
                        end do
                    end do
                end do
                call MPI_ALLTOALL(rcvbuf, nxl*nyu*nzl, MPI_DOUBLE_COMPLEX, &
                                  sndbuf, nxl*nyu*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
                do i = 1, nxl
                    do j = 1, nyu
                        do ip = 1, nproc_fft
                            do k = 1 + nzl*(ip - 1), nzl*ip
                                fp(k + nzu*(i - 1) + nzu*nxl*(j - 1)) = sndbuf(i + nxl*(j - 1) + nxl*nyu*(k - 1))
                            end do
                        end do
                    end do
                end do
                call cft_1z(fp, nxl*nyu, nzu, nzu, 1)
                timepp = cclock()
                if (commensurate_lattice) then
                do i = 1, nxl
                    do j = 1, nyu
                        do k = 1, nzu
                            vhartree(i + nxl*(j - 1) + nxl*nyu*(k - 1)) = fp(k + nzu*(i - 1) + nzu*nxl*(j - 1))
                        end do
                    end do
                end do
                end if
                call reduce_base_real(1, vh_test, sub_comm_fft%comm, -1)
                call reduce_base_real(1, vh_att, sub_comm_fft%comm, -1)
                call reduce_base_real(1, eh_ew, sub_comm_fft%comm, -1)

            else ! open system

                if (commensurate_lattice) then
                    fp = zzero
                    do i = 1, nxl
                        do j = 1, nyu
                            do k = 1, nzu
                                fp(k + 2*nzu*(i - 1) + 2*nzu*nxl2*(j - 1)) = dent(i + nxl*(j - 1) + nxl*nyu*(k - 1))
                            end do
                        end do
                    end do
                end if
!   cost=sum(abs(fp(:)))
!   call reduce_base_real(1,cost,sub_comm_fft%comm,-1)
!   if(rank.eq.0)  write(6,*) ' Input fp =',cost
                time_uploadfft = time_uploadfft + cclock() - timep
                call cft_1z(fp, 2*nxl2*nyu, 2*nzu, 2*nzu, -1)

                do i = 1, nxl2
                    do j = 1, 2*nyu
                        do ip = 1, nproc_fft
                            do k = 1 + nzl2*(ip - 1), nzl2*ip
                                sndbuf(i + nxl2*(j - 1) + 2*nxl2*nyu*(k - 1)) = fp(k + 2*nzu*(i - 1) + 2*nzu*nxl2*(j - 1))
                            end do
                        end do
                    end do
                end do
                call MPI_ALLTOALL(sndbuf, 2*nxl2*nyu*nzl2, MPI_DOUBLE_COMPLEX, &
                                  rcvbuf, 2*nxl2*nyu*nzl2, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
                do ip = 1, nproc_fft
                    do j = 1, 2*nyu
                        do k = 1, nzl2
                            jj = 2*nxu*(j - 1) + 4*nxu*nyu*(k - 1)
                            kk = nxl2*(j - 1) + 2*nxl2*nyu*(k - 1 + (ip - 1)*nzl2)
                            do i = 1, nxl2
                                fp(ip + nproc_fft*(i - 1) + jj) = rcvbuf(i + kk)
                            end do
                        end do
                    end do
                end do

                call cft_2xy(fp, nzl2, 2*nxu, 2*nyu, 2*nxu, 2*nyu, -1)

                do i = 1, 2*nxu
                    if (i > (2*nxu + 1)/2) then
                        qx = 0.5d0*fx*(i - 2*nxu - 1)
                    else
                        qx = 0.5d0*fx*(i - 1)
                    end if
                    do j = 1, 2*nyu
                        if (j > (2*nyu + 1)/2) then
                            qy = 0.5d0*fy*(j - 2*nyu - 1)
                        else
                            qy = 0.5d0*fy*(j - 1)
                        end if
                        do kk = 1, nzl2
                            k = kk + nzl2*rank
                            if (k > (2*nzu + 1)/2) then
                                qz = 0.5d0*fz*(k - 2*nzu - 1)
                            else
                                qz = 0.5d0*fz*(k - 1)
                            end if
                            q2 = qx*qx + qy*qy + qz*qz
                            iq = i + 2*nxu*(j - 1) + 4*nxu*nyu*(kk - 1)

                            cost_att = 1.d0
                            if (scale_hartree .gt. 0.d0) cost_att = cost_att*scale_hartreen

                            vh_test = vh_test + 4.d0*vhartreeq(iq)*fp(iq)*dconjg(fp(iq))*voltot*cost_att**2
                            fp(iq) = fp(iq)*vhartreeq(iq)*cost_att

                        end do
                    end do
                end do

                call cft_2xy(fp, nzl2, 2*nxu, 2*nyu, 2*nxu, 2*nyu, 1)

                do ip = 1, nproc_fft
                    do j = 1, 2*nyu
                        do k = 1, nzl2
                            jj = 2*nxu*(j - 1) + 4*nxu*nyu*(k - 1)
                            kk = nxl2*(j - 1) + 2*nxl2*nyu*(k - 1 + (ip - 1)*nzl2)
                            do i = 1, nxl2
                                rcvbuf(i + kk) = fp((ip + nproc_fft*(i - 1)) + jj)
                            end do
                        end do
                    end do
                end do
                call MPI_ALLTOALL(rcvbuf, 2*nxl2*nyu*nzl2, MPI_DOUBLE_COMPLEX, &
                                  sndbuf, 2*nxl2*nyu*nzl2, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
                do i = 1, nxl2
                    do j = 1, 2*nyu
                        do ip = 1, nproc_fft
                            do k = 1 + nzl2*(ip - 1), nzl2*ip
                                fp(k + 2*nzu*(i - 1) + 2*nzu*nxl2*(j - 1)) = sndbuf(i + nxl2*(j - 1) + 2*nxl2*nyu*(k - 1))
                            end do
                        end do
                    end do
                end do
                call cft_1z(fp, 2*nxl2*nyu, 2*nzu, 2*nzu, 1)
                timepp = cclock()
                if (commensurate_lattice) then
                do i = 1, nxl
                    do j = 1, nyu
                        do k = 1, nzu
                            vhartree(i + nxl*(j - 1) + nxl*nyu*(k - 1)) = fp(k + 2*nzu*(i - 1) + 2*nzu*nxl2*(j - 1))
                        end do
                    end do
                end do
                end if
!   cost=sum(abs(fp(:)))
!   call reduce_base_real(1,cost,sub_comm_fft%comm,-1)
!   if(rank.eq.0)  write(6,*) ' Output fp =',cost
                call reduce_base_real(1, vh_test, sub_comm_fft%comm, -1)
            end if
        end if

        if (.not. commensurate_lattice) then
            if (iespbc) then
                call bcast_real(vh_test, 1, 0, commrep_mpi)
                call bcast_real(vh_att, 1, 0, commrep_mpi)
                call bcast_real(eh_ew, 1, 0, commrep_mpi)
            else
                call bcast_real(vh_test, 1, 0, commrep_mpi)
            end if
            call upload_vhartree_fromfp
        end if
!     time_uploadfft=time_uploadfft+cclock()-timepp
#else

        ! scalar algorithm
        vh_test = 0.d0
        vh_att = 0.d0
        eh_ew = 0.d0

        if (iespbc) then

            fp = dcmplx(dent)

            call cfft3d(fp, nxu, nyu, nzu, nxu, nyu, nzu, -1)

            do i = 1, nxu
                if (i > (nxu + 1)/2) then
                    qx = fx*(i - nxu - 1)
                else
                    qx = fx*(i - 1)
                end if
                do j = 1, nyu
                    if (j > (nyu + 1)/2) then
                        qy = fy*(j - nyu - 1)
                    else
                        qy = fy*(j - 1)
                    end if
                    do k = 1, nzu
                        if (k > (nzu + 1)/2) then
                            qz = fz*(k - nzu - 1)
                        else
                            qz = fz*(k - 1)
                        end if
                        if (yes_tilted) then
                            vecscra(:) = recip(1, :)*qx + recip(2, :)*qy + recip(3, :)*qz
                            q2 = sum(vecscra(:)**2)
                        else
                            q2 = qx*qx + qy*qy + qz*qz
                        end if

                        iq = i + nxu*(j - 1) + nxu*nyu*(k - 1)

                        cost_att = 1.d0
                        if (scale_hartree .gt. 0.d0) cost_att = cost_att*scale_hartreen

                        if (q2 .gt. 0) then
                            vh_test = vh_test + 0.5d0*vhartreeq(iq)*fp(iq)*dconjg(fp(iq))*voltot*cost_att**2
                            vh_att = vh_att + TWO_PI/((2.d0 - 2.d0*cos(qx*axu))/axu**2 &
                                                      + (2 - 2.d0*cos(qy*ayu))/ayu**2 &
                                                      + (2.d0 - 2.d0*cos(qz*azu))/azu**2)*fp(iq)*dconjg(fp(iq))*voltot
                            eh_ew = eh_ew + TWO_PI/q2*fp(iq)*dconjg(fp(iq))*voltot*cost_att**2
                            fp(iq) = fp(iq)*vhartreeq(iq)*cost_att
                        else
                            fp(iq) = (0.d0, 0.d0)
                        end if
                    end do
                end do
            end do

            call cfft3d(fp, nxu, nyu, nzu, nxu, nyu, nzu, +1)

            vhartree = fp

        else ! open system
            !               loading fp
            indmesh = 0
            !               write(6,*) ' Input dent =',sum(dent(1:nx*ny*nz))
            fp = (0.d0, 0.d0)
            do k = 1, nzu
                do j = 1, nyu
                    do i = 1, nxu
                        indmesh = indmesh + 1
                        indfp = 4*nxu*nyu*(k - 1) + 2*nxu*(j - 1) + i
                        fp(indfp) = dent(indmesh)
                    end do
                end do
            end do
            !                write(6,*) ' fp after loading ',sum(abs(fp(1:8*nx*ny*nz)))
            call cfft3d(fp, 2*nxu, 2*nyu, 2*nzu, 2*nxu, 2*nyu, 2*nzu, -1)
            !                write(6,*) ' fp after fft  ',sum(abs(fp(1:8*nx*ny*nz)))
            do i = 1, 2*nxu
                if (i > (2*nxu + 1)/2) then
                    qx = 0.5d0*fx*(i - 2*nxu - 1)
                else
                    qx = 0.5d0*fx*(i - 1)
                end if
                do j = 1, 2*nyu
                    if (j > (2*nyu + 1)/2) then
                        qy = 0.5d0*fy*(j - 2*nyu - 1)
                    else
                        qy = 0.5d0*fy*(j - 1)
                    end if
                    do k = 1, 2*nzu
                        if (k > (2*nzu + 1)/2) then
                            qz = 0.5d0*fz*(k - 2*nzu - 1)
                        else
                            qz = 0.5d0*fz*(k - 1)
                        end if
                        iq = i + 2*nxu*(j - 1) + 4*nxu*nyu*(k - 1)
                        q2 = qx*qx + qy*qy + qz*qz
                        cost_att = 1.d0
                        if (scale_hartree .gt. 0.d0) cost_att = cost_att*scale_hartreen
                        vh_test = vh_test + 4.d0*vhartreeq(iq)*fp(iq)*dconjg(fp(iq))*voltot*cost_att**2
                        fp(iq) = fp(iq)*vhartreeq(iq)*cost_att
                    end do
                end do
            end do

            call cfft3d(fp, 2*nxu, 2*nyu, 2*nzu, 2*nxu, 2*nyu, 2*nzu, +1)

            ! restoring vhartree in real space
            indmesh = 0
            do k = 1, nzu
                do j = 1, nyu
                    do i = 1, nxu
                        indmesh = indmesh + 1
                        indfp = 4*nxu*nyu*(k - 1) + 2*nxu*(j - 1) + i
                        vhartree(indmesh) = fp(indfp)
                    end do
                end do
            end do

        end if

#endif

        if (double_mesh) then

            rion_upload(:) = rion_ref(:) - (nx + 1)/2.d0*ax*at(:, 1) - (ny + 1)/2.d0*ay*at(:, 2)&
                    & - (nz + 1)/2.d0*az*at(:, 3)

            !     initialize counter
            indmesh = meshproc
            if (corr_hartree) then
                if (from_ions) then
                    do ii = 1, nion
                        !        locate atom
                        if (zetar(ii) .ge. minz_at) then
                            scalea = scale_z
                            if (nx_at .gt. 0) then
                                nx0 = nx_at
                                ny0 = ny_at
                                nz0 = nz_at
                            else
                                nx0 = (2*l0_at)/ax + 1
                                ny0 = (2*l0_at)/ay + 1
                                nz0 = (2*l0_at)/az + 1
                            end if

                            allocate (weightx(nx0 + 2), weighty(ny0 + 2), weightz(nz0 + 2))
                            weightx = 0.d0
                            weighty = 0.d0
                            weightz = 0.d0

                            allocate (vhartree_global(nx0 + 2, ny0 + 2, nz0 + 2))
                            vhartree_global = 0.d0

                            call set_interval(scalea, nx0, ny0, nz0, rion(1, ii), rion_upload, ax, ay, az &
                                    , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .true.&
                                    &, weightx, weighty, weightz)

                            call upload_vector_global(vhartree, vhartree_global, nxi, nxf, nyi, nyf, nzi, nzf)

                            call upload_vhartree(indmesh, scalea, nxi, nxf, nyi, nyf, nzi, nzf, .false.)

                            deallocate (weightx, weighty, weightz)

                            axn = ax
                            ayn = ay
                            azn = az

                            nx0n = scalea*(nx0 - 1) + 3
                            ny0n = scalea*(ny0 - 1) + 3
                            nz0n = scalea*(nz0 - 1) + 3

                            allocate (weightx(nx0n), weighty(ny0n), weightz(nz0n))
                            weightx = 0.d0
                            weighty = 0.d0
                            weightz = 0.d0

                            nx0n = nx0
                            ny0n = ny0
                            nz0n = nz0

                            call set_interval(scalea, nx0n, ny0n, nz0n, rion(1, ii)&
                                    &, rion_upload, axn, ayn, azn, nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr&
                                    &, .false., weightx, weighty, weightz)

                            call upload_vhartree(indmesh, scalea, nxi, nxf, nyi, nyf, nzi, nzf, .true.)

                            deallocate (weightx, weighty, weightz, vhartree_global)
                        end if
                    end do
                elseif (nx_at .gt. 0) then
                    scalea = scale_z
                    nx0 = nx_at
                    ny0 = ny_at
                    nz0 = nz_at

                    allocate (weightx(nx0 + 2), weighty(ny0 + 2), weightz(nz0 + 2))
                    weightx = 0.d0
                    weighty = 0.d0
                    weightz = 0.d0

                    allocate (vhartree_global(nx0 + 2, ny0 + 2, nz0 + 2))
                    vhartree_global = 0.d0

                    call set_interval(scalea, nx0, ny0, nz0, rion_from, rion_upload, ax, ay, az &
                            , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .true.&
                            &, weightx, weighty, weightz)

                    call upload_vector_global(vhartree, vhartree_global, nxi, nxf, nyi, nyf, nzi, nzf)

                    call upload_vhartree(indmesh, scalea, nxi, nxf, nyi, nyf, nzi, nzf, .false.)

                    deallocate (weightx, weighty, weightz)

                    axn = ax
                    ayn = ay
                    azn = az

                    nx0n = scalea*(nx0 - 1) + 3
                    ny0n = scalea*(ny0 - 1) + 3
                    nz0n = scalea*(nz0 - 1) + 3

                    allocate (weightx(nx0n), weighty(ny0n), weightz(nz0n))
                    weightx = 0.d0
                    weighty = 0.d0
                    weightz = 0.d0

                    nx0n = nx0
                    ny0n = ny0
                    nz0n = nz0

                    call set_interval(scalea, nx0n, ny0n, nz0n, rion_from&
                            &, rion_upload, axn, ayn, azn, nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr&
                            &, .false., weightx, weighty, weightz)

                    call upload_vhartree(indmesh, scalea, nxi, nxf, nyi, nyf, nzi, nzf, .true.)

                    deallocate (weightx, weighty, weightz, vhartree_global)
                end if
                !     Add the delta also to the first part
                !  vhartree(1:meshproc)=vhartree(1:meshproc)+vmax0*(dent(1:meshproc)-dens0)
            elseif (.not. do_hartree) then ! already defined by the fft
                vhartree(meshproc + 1:meshproc_tot) = 0.d0
            end if
        end if
        time_fft = time_fft + cclock() - timep
        !   cost=sum(abs(vhartree(:)))
        !   call reduce_base_real(1,cost,commrep_mpi,-1)
        !   if(rank.eq.0) write(6,*) ' Output vhartree =',cost

        return

    end subroutine update_vhartree

    subroutine upload_vector_global(vector, vector_global, nxi, nxf, nyi, nyf, nzi, nzf)
        implicit none
        integer indmesh, i, j, k, ip, jp, kp, indproc, nxi, nxf, nyi, nyf, nzi, nzf&
                &, nxiu, nyiu, nziu, nxfu, nyfu, nzfu, signi, signj, signk, ierr&
                &, nx0, ny0, nz0
        real*8 vector_global(nxf - nxi + 1, nyf - nyi + 1, nzf - nzi + 1), vector(*), countm
        integer, external :: mod_true
        vector_global = 0.d0
        indmesh = 0
        indproc = 0
        nx0 = nxf - nxi + 1
        ny0 = nyf - nyi + 1
        nz0 = nzf - nzi + 1
        if (iespbc) then
            countm = 0.d0
            signi = 1
            signj = 1
            signk = 1
            nxiu = mod_true(nxi - 1, nx) + 1

            if (nx0 .le. nx) then
                nxfu = mod_true(nxf - 1, nx) + 1
            else
                nxfu = nxi + nx - 1
                nxfu = mod_true(nxfu - 1, nx) + 1
            end if

            if (nxfu .lt. nxiu) signi = -1

            nyiu = mod_true(nyi - 1, ny) + 1
            if (ny0 .le. ny) then
                nyfu = mod_true(nyf - 1, ny) + 1
            else
                nyfu = nyi + ny - 1
                nyfu = mod_true(nyfu - 1, ny) + 1
            end if

            if (nyfu .lt. nyiu) signj = -1

            nziu = mod_true(nzi - 1, nz) + 1
            if (nz0 .le. nz) then
                nzfu = mod_true(nzf - 1, nz) + 1
            else
                nzfu = nzi + nz - 1
                nzfu = mod_true(nzfu - 1, nz) + 1
            end if

            if (nzfu .lt. nziu) signk = -1

            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        if (indproc .eq. rankrep) then
                            indmesh = indmesh + 1
                            if (((signi .eq. 1 .and. i .ge. nxiu .and. i .le. nxfu) .or. &
                                 (signi .eq. -1 .and. (i .ge. nxiu .or. i .le. nxfu))) .and. &
                                ((signj .eq. 1 .and. j .ge. nyiu .and. j .le. nyfu) .or. &
                                 (signj .eq. -1 .and. (j .ge. nyiu .or. j .le. nyfu))) .and. &
                                ((signk .eq. 1 .and. k .ge. nziu .and. k .le. nzfu) .or. &
                                 (signk .eq. -1 .and. (k .ge. nziu .or. k .le. nzfu)))) then
                                ip = i
                                if (signi .lt. 0 .and. ip - nxiu + 1 .le. 0) ip = ip + nx
                                jp = j
                                if (signj .lt. 0 .and. jp - nyiu + 1 .le. 0) jp = jp + ny
                                kp = k
                                if (signk .lt. 0 .and. kp - nziu + 1 .le. 0) kp = kp + nz
                                countm = countm + 1.d0
                                vector_global(ip - nxiu + 1, jp - nyiu + 1, kp - nziu + 1) = vector(indmesh)
                            end if
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                    end do
                end do
            end do
#ifdef  PARALLEL
            call reduce_base_real(1, countm, commrep_mpi, -1)
#endif

            if (countm .ne. size(vector_global) .and. nx0 .le. nx .and. ny0 .le. ny .and. nz0 .le. nz)&
                    &call error(' upload_vector_global  ', ' not read all mesh points as it should', 1, rankrep)
            ! Apply PBC to load the unknown values
            if (nx0 .gt. nx .or. ny0 .gt. ny .or. nz0 .gt. nz) then
                do k = nzi, nzf
                    do j = nyi, nyf
                        do i = nxi, nxf
                            vector_global(i - nxi + 1, j - nyi + 1, k - nzi + 1) = &
                                    &vector_global(mod_true(i - nxi, nx) + 1&
                                    &, mod_true(j - nyi, ny) + 1, mod_true(k - nzi, nz) + 1)
                        end do
                    end do
                end do
            end if

        else ! if iespbc
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        if (indproc .eq. rankrep) then
                            indmesh = indmesh + 1
                            if (i .ge. nxi .and. i .le. nxf&
                                    &.and. j .ge. nyi .and. j .le. nyf&
                                    &.and. k .ge. nzi .and. k .le. nzf) then
                                vector_global(i - nxi + 1, j - nyi + 1, k - nzi + 1) = vector(indmesh)
                            end if
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                    end do
                end do
            end do
        end if
#ifdef PARALLEL
        call reduce_base_real(size(vector_global), vector_global, commrep_mpi, -1)
#endif
    end subroutine upload_vector_global

    subroutine upload_vhartree(indmesh, scalea, nxi, nxf, nyi, nyf, nzi, nzf, add)
        implicit none
        integer indproc, i, j, k, ii, jj, kk, nxi, nxf, nyi, nyf, nzi, nzf, scalea, indmesh, irest
        integer i1, i2, i3, ind1, ind2, ind3, ierr
        logical add
        real*8 px, pxm, py, pym, pz, pzm, ppx(4), ppy(4), ppz(4)
        indproc = 0

        do k = nzi, nzf
            do j = nyi, nyf
                do i = nxi, nxf
                    if (indproc .eq. rankrep) then
                        indmesh = indmesh + 1
                        if (add) then
                            ! When you add you do not consider the first and last points.
                            ! Here we may add some interpolation formula of the potential in a finer grid.
                            !   Linear interpolation in the given mesh
                            if (linear) then
                                if (i .eq. nxi) then
                                    ii = 1
                                    irest = scalea - 1
                                elseif (i .eq. nxf) then
                                    ii = nx0 + 1
                                    irest = 1
                                else
                                    ii = (i - nxi - 1)/scalea + 2
                                    irest = (i - nxi - 1) - (ii - 2)*scalea
                                end if
                                px = dble(irest)/dble(scalea)
                                pxm = 1.d0 - px

                                if (j .eq. nyi) then
                                    jj = 1
                                    irest = scalea - 1
                                elseif (j .eq. nyf) then
                                    jj = ny0 + 1
                                    irest = 1
                                else
                                    jj = (j - nyi - 1)/scalea + 2
                                    irest = (j - nyi - 1) - (jj - 2)*scalea
                                end if
                                py = dble(irest)/dble(scalea)
                                pym = 1.d0 - py

                                if (k .eq. nzi) then
                                    kk = 1
                                    irest = scalea - 1
                                elseif (k .eq. nzf) then
                                    kk = nz0 + 1
                                    irest = 1
                                else
                                    kk = (k - nzi - 1)/scalea + 2
                                    irest = (k - nzi - 1) - (kk - 2)*scalea
                                end if
                                pz = dble(irest)/dble(scalea)
                                pzm = 1.d0 - pz

                                vhartree(indmesh) = 0.d0

                                vhartree(indmesh) = &
                                    px*py*pz*vhartree_global(ii + 1, jj + 1, kk + 1) &
                                    + px*py*pzm*vhartree_global(ii + 1, jj + 1, kk) &
                                    + px*pym*pz*vhartree_global(ii + 1, jj, kk + 1) &
                                    + px*pym*pzm*vhartree_global(ii + 1, jj, kk) &
                                    + pxm*py*pz*vhartree_global(ii, jj + 1, kk + 1) &
                                    + pxm*py*pzm*vhartree_global(ii, jj + 1, kk) &
                                    + pxm*pym*pz*vhartree_global(ii, jj, kk + 1) &
                                    + pxm*pym*pzm*vhartree_global(ii, jj, kk)

                            else

                                !  ii is the point where interpolation is done 1<=ii<=nx0+2
                                !  px is the value of the interpolation value x--> f(x)
                                !  in unit of ax (the original not scaled)

                                if (i .eq. nxi) then
                                    ii = 2
                                    px = -1.d0/scalea
                                elseif (i .eq. nxf) then
                                    ii = nx0
                                    px = 1.d0 + 1.d0/scalea
                                else
                                    ii = (i - nxi - 1)/scalea + 2
                                    irest = (i - nxi - 1) - (ii - 2)*scalea
                                    px = dble(irest)/scalea
                                end if

                                if (j .eq. nyi) then
                                    jj = 2
                                    py = -1.d0/scalea
                                elseif (j .eq. nyf) then
                                    jj = ny0
                                    py = 1.d0 + 1.d0/scalea
                                else
                                    jj = (j - nyi - 1)/scalea + 2
                                    irest = (j - nyi - 1) - (jj - 2)*scalea
                                    py = dble(irest)/scalea
                                end if

                                if (k .eq. nzi) then
                                    kk = 2
                                    pz = -1.d0/scalea
                                elseif (k .eq. nzf) then
                                    kk = nz0
                                    pz = 1.d0 + 1.d0/scalea
                                else
                                    kk = (k - nzi - 1)/scalea + 2
                                    irest = (k - nzi - 1) - (kk - 2)*scalea
                                    pz = dble(irest)/scalea
                                end if

                                call cubic(px, ppx)
                                call cubic(py, ppy)
                                call cubic(pz, ppz)
                                vhartree(indmesh) = 0.d0
                                do i3 = 1, 4
                                    ind3 = kk - 2 + i3
                                    do i2 = 1, 4
                                        ind2 = jj - 2 + i2
                                        do i1 = 1, 4
                                            ind1 = ii - 2 + i1
                                            if (ppx(i1)*ppy(i2)*ppz(i3) .ne. 0.d0) then
                                                vhartree(indmesh) = vhartree(indmesh) + ppx(i1)*ppy(i2)*ppz(i3)* &
                                                        & vhartree_global(ind1, ind2, ind3)
                                            end if
                                        end do
                                    end do
                                end do
                            end if
                        else
                            vhartree(indmesh) = vhartree_global(i - nxi + 1, j - nyi + 1, k - nzi + 1)
                        end if
                        !             vhartree(indmesh)=vhartree(indmesh)+vmax0*(dent(indmesh)-dens0)
                    end if
                    indproc = indproc + 1
                    if (indproc .eq. nprocrep) indproc = 0
                end do
            end do
        end do

    end subroutine upload_vhartree

    !
    ! Input:  dent, dent unchanged
    ! Output: gradt convoluted with mixingder q^2/ ( tfcut+ q^2)
    ! This subroutine compute the gradient of the electronic density.
    ! It is supposed to be used in the evaluation of the GGA XC potential
    ! which is not implemented yet in the current versione of the code.
    ! Therefore this subroutine is UNUSED for the moment.
    !

    subroutine evalgrad

        use setup, only: dent, gradt, voltot, volmesh
        implicit none

        real(8) :: q, q2, fac, qx, qy, qz
        integer :: iq, ixl, izl, ierr, kk, ip, jj, i, j, k

#ifdef PARALLEL
        include 'mpif.h'
        ! x-grad
        gradt = 0.d0
        if (sub_comm_fft%yesin) then

            do i = 1, nxl
                do j = 1, ny
                    do k = 1, nz
                        fp(k + nz*(i - 1) + nz*nxl*(j - 1)) = dent(i + nxl*(j - 1) + nxl*ny*(k - 1))
                    end do
                end do
            end do

            call cft_1z(fp, nxl*ny, nz, nz, -1)

            do i = 1, nxl
                do j = 1, ny
                    do ip = 1, nproc_fft
                        do k = 1 + nzl*(ip - 1), nzl*ip
                            sndbuf(i + nxl*(j - 1) + nxl*ny*(k - 1)) = fp(k + nz*(i - 1) + nz*nxl*(j - 1))
                        end do
                    end do
                end do
            end do
            call MPI_ALLTOALL(sndbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, &
                              rcvbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
            do ip = 1, nproc_fft
                do j = 1, ny
                    do k = 1, nzl
                        jj = nx*(j - 1) + nx*ny*(k - 1)
                        kk = nxl*(j - 1) + nxl*ny*(k - 1 + (ip - 1)*nzl)
                        do i = 1, nxl
                            fp(ip + nproc_fft*(i - 1) + jj) = rcvbuf(i + kk)
                        end do
                    end do
                end do
            end do

            call cft_2xy(fp, nzl, nx, ny, nx, ny, -1)

            do i = 1, nx
                if (i > (nx + 1)/2) then
                    qx = fx*(i - nx - 1)
                else
                    qx = fx*(i - 1)
                end if
                !             qx = qx*qx
                !             qx=dsin(qx*ax)/ax
                do j = 1, ny
                    if (j > (ny + 1)/2) then
                        qy = fy*(j - ny - 1)
                    else
                        qy = fy*(j - 1)
                    end if
                    !                qy = qy*qy
                    !                if(qy.ne.0.d0) qy = 2.d0/ay**2*(1.d0-dcos(qy*ay))

                    do kk = 1, nzl
                        k = kk + nzl*sub_comm_fft%rank
                        if (k > (nz + 1)/2) then
                            qz = fz*(k - nz - 1)
                        else
                            qz = fz*(k - 1)
                        end if
                        !                if(qz.ne.0.d0) qz = 2.d0/az**2*(1.d0-dcos(qz*az))
                        !                   q2 = qx + qy + qz*qz
                        iq = i + nx*(j - 1) + nx*ny*(kk - 1)
                        fp(iq) = fp(iq)*dcmplx(0.d0, qx)
                    end do
                end do
            end do

            call cft_2xy(fp, nzl, nx, ny, nx, ny, 1)

            do ip = 1, nproc_fft
                do j = 1, ny
                    do k = 1, nzl
                        jj = nx*(j - 1) + nx*ny*(k - 1)
                        kk = nxl*(j - 1) + nxl*ny*(k - 1 + (ip - 1)*nzl)
                        do i = 1, nxl
                            rcvbuf(i + kk) = fp((ip + nproc_fft*(i - 1)) + jj)
                        end do
                    end do
                end do
            end do
            call MPI_ALLTOALL(rcvbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, &
                              sndbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
            do i = 1, nxl
                do j = 1, ny
                    do ip = 1, nproc_fft
                        do k = 1 + nzl*(ip - 1), nzl*ip
                            fp(k + nz*(i - 1) + nz*nxl*(j - 1)) = sndbuf(i + nxl*(j - 1) + nxl*ny*(k - 1))
                        end do
                    end do
                end do
            end do
            call cft_1z(fp, nxl*ny, nz, nz, 1)
            do i = 1, nxl
                do j = 1, ny
                    do k = 1, nz
                        gradt(1, i + nxl*(j - 1) + nxl*ny*(k - 1)) = dreal(fp(k + nz*(i - 1) + nz*nxl*(j - 1)))
                    end do
                end do
            end do

            !          y-grad
            do i = 1, nxl
                do j = 1, ny
                    do k = 1, nz
                        fp(k + nz*(i - 1) + nz*nxl*(j - 1)) = dent(i + nxl*(j - 1) + nxl*ny*(k - 1))
                    end do
                end do
            end do

            call cft_1z(fp, nxl*ny, nz, nz, -1)

            do i = 1, nxl
                do j = 1, ny
                    do ip = 1, nproc_fft
                        do k = 1 + nzl*(ip - 1), nzl*ip
                            sndbuf(i + nxl*(j - 1) + nxl*ny*(k - 1)) = fp(k + nz*(i - 1) + nz*nxl*(j - 1))
                        end do
                    end do
                end do
            end do
            call MPI_ALLTOALL(sndbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, &
                              rcvbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
            do ip = 1, nproc_fft
                do j = 1, ny
                    do k = 1, nzl
                        jj = nx*(j - 1) + nx*ny*(k - 1)
                        kk = nxl*(j - 1) + nxl*ny*(k - 1 + (ip - 1)*nzl)
                        do i = 1, nxl
                            fp(ip + nprocrep*(i - 1) + jj) = rcvbuf(i + kk)
                        end do
                    end do
                end do
            end do

            call cft_2xy(fp, nzl, nx, ny, nx, ny, -1)

            do i = 1, nx
                if (i > (nx + 1)/2) then
                    qx = fx*(i - nx - 1)
                else
                    qx = fx*(i - 1)
                end if
                !             qx = qx*qx
                do j = 1, ny
                    if (j > (ny + 1)/2) then
                        qy = fy*(j - ny - 1)
                    else
                        qy = fy*(j - 1)
                    end if
                    !                qy = qy*qy
                    !                if(qy.ne.0.d0) qy = 2.d0/ay**2*(1.d0-dcos(qy*ay))
                    !                qy=dsin(qy*ay)/ay

                    do kk = 1, nzl
                        k = kk + nzl*sub_comm_fft%rank
                        if (k > (nz + 1)/2) then
                            qz = fz*(k - nz - 1)
                        else
                            qz = fz*(k - 1)
                        end if

                        iq = i + nx*(j - 1) + nx*ny*(kk - 1)

                        !                if(qz.ne.0.d0) qz = 2.d0/az**2*(1.d0-dcos(qz*az))
                        !                   q2 = qx + qy + qz*qz
                        fp(iq) = fp(iq)*dcmplx(0.d0, qy)
                    end do
                end do
            end do

            call cft_2xy(fp, nzl, nx, ny, nx, ny, 1)

            do ip = 1, nproc_fft
                do j = 1, ny
                    do k = 1, nzl
                        jj = nx*(j - 1) + nx*ny*(k - 1)
                        kk = nxl*(j - 1) + nxl*ny*(k - 1 + (ip - 1)*nzl)
                        do i = 1, nxl
                            rcvbuf(i + kk) = fp((ip + nproc_fft*(i - 1)) + jj)
                        end do
                    end do
                end do
            end do
            call MPI_ALLTOALL(rcvbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, &
                              sndbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
            do i = 1, nxl
                do j = 1, ny
                    do ip = 1, nprocrep
                        do k = 1 + nzl*(ip - 1), nzl*ip
                            fp(k + nz*(i - 1) + nz*nxl*(j - 1)) = sndbuf(i + nxl*(j - 1) + nxl*ny*(k - 1))
                        end do
                    end do
                end do
            end do
            call cft_1z(fp, nxl*ny, nz, nz, 1)
            do i = 1, nxl
                do j = 1, ny
                    do k = 1, nz
                        gradt(2, i + nxl*(j - 1) + nxl*ny*(k - 1)) = dreal(fp(k + nz*(i - 1) + nz*nxl*(j - 1)))
                    end do
                end do
            end do

            !          z-grad
            do i = 1, nxl
                do j = 1, ny
                    do k = 1, nz
                        fp(k + nz*(i - 1) + nz*nxl*(j - 1)) = dent(i + nxl*(j - 1) + nxl*ny*(k - 1))
                    end do
                end do
            end do

            call cft_1z(fp, nxl*ny, nz, nz, -1)

            do i = 1, nxl
                do j = 1, ny
                    do ip = 1, nproc_fft
                        do k = 1 + nzl*(ip - 1), nzl*ip
                            sndbuf(i + nxl*(j - 1) + nxl*ny*(k - 1)) = fp(k + nz*(i - 1) + nz*nxl*(j - 1))
                        end do
                    end do
                end do
            end do
            call MPI_ALLTOALL(sndbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, &
                              rcvbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
            do ip = 1, nproc_fft
                do j = 1, ny
                    do k = 1, nzl
                        jj = nx*(j - 1) + nx*ny*(k - 1)
                        kk = nxl*(j - 1) + nxl*ny*(k - 1 + (ip - 1)*nzl)
                        do i = 1, nxl
                            fp(ip + nprocrep*(i - 1) + jj) = rcvbuf(i + kk)
                        end do
                    end do
                end do
            end do

            call cft_2xy(fp, nzl, nx, ny, nx, ny, -1)

            do i = 1, nx
                if (i > (nx + 1)/2) then
                    qx = fx*(i - nx - 1)
                else
                    qx = fx*(i - 1)
                end if
                do j = 1, ny
                    if (j > (ny + 1)/2) then
                        qy = fy*(j - ny - 1)
                    else
                        qy = fy*(j - 1)
                    end if
                    !                qy = qy*qy
                    !                if(qy.ne.0.d0) qy = 2.d0/ay**2*(1.d0-dcos(qy*ay))

                    do kk = 1, nzl
                        k = kk + nzl*sub_comm_fft%rank
                        if (k > (nz + 1)/2) then
                            qz = fz*(k - nz - 1)
                        else
                            qz = fz*(k - 1)
                        end if
                        !                   qz=dsin(az*qz)/az
                        iq = i + nx*(j - 1) + nx*ny*(kk - 1)
                        fp(iq) = fp(iq)*dcmplx(0.d0, qz)
                    end do
                end do
            end do

            call cft_2xy(fp, nzl, nx, ny, nx, ny, 1)

            do ip = 1, nproc_fft
                do j = 1, ny
                    do k = 1, nzl
                        jj = nx*(j - 1) + nx*ny*(k - 1)
                        kk = nxl*(j - 1) + nxl*ny*(k - 1 + (ip - 1)*nzl)
                        do i = 1, nxl
                            rcvbuf(i + kk) = fp((ip + nproc_fft*(i - 1)) + jj)
                        end do
                    end do
                end do
            end do
            call MPI_ALLTOALL(rcvbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, &
                              sndbuf, nxl*ny*nzl, MPI_DOUBLE_COMPLEX, sub_comm_fft%comm, ierr)
            do i = 1, nxl
                do j = 1, ny
                    do ip = 1, nproc_fft
                        do k = 1 + nzl*(ip - 1), nzl*ip
                            fp(k + nz*(i - 1) + nz*nxl*(j - 1)) = sndbuf(i + nxl*(j - 1) + nxl*ny*(k - 1))
                        end do
                    end do
                end do
            end do
            call cft_1z(fp, nxl*ny, nz, nz, 1)
            do i = 1, nxl
                do j = 1, ny
                    do k = 1, nz
                        gradt(3, i + nxl*(j - 1) + nxl*ny*(k - 1)) = dreal(fp(k + nz*(i - 1) + nz*nxl*(j - 1)))
                    end do
                end do
            end do
        end if

#else

        !  x-grad
        fp = dent

        call cfft3d(fp, nx, ny, nz, nx, ny, nz, -1)

        do i = 1, nx
            if (i > (nx + 1)/2) then
                qx = fx*(i - nx - 1)
            else
                qx = fx*(i - 1)
            end if
            !             qx = dsin(qx*ax)/ax

            do j = 1, ny
                if (j > (ny + 1)/2) then
                    qy = fy*(j - ny - 1)
                else
                    qy = fy*(j - 1)
                end if
                !                qy=qy*qy
                !                if(qy.ne.0.d0) qy = 2.d0/ay**2*(1.d0-dcos(qy*ay))
                do k = 1, nz
                    if (k > (nz + 1)/2) then
                        qz = fz*(k - nz - 1)
                    else
                        qz = fz*(k - 1)
                    end if

                    !                   if(qz.ne.0.d0) qz = 2.d0/az**2*(1.d0-dcos(qz*az))
                    !                   q2 = qx + qy + qz*qz
                    iq = i + nx*(j - 1) + nx*ny*(k - 1)
                    fp(iq) = fp(iq)*dcmplx(0.d0, qx)
                end do
            end do
        end do

        call cfft3d(fp, nx, ny, nz, nx, ny, nz, +1)

        gradt(1, :) = dreal(fp(:))

        !  y-grad

        fp = dent

        call cfft3d(fp, nx, ny, nz, nx, ny, nz, -1)

        do i = 1, nx
            if (i > (nx + 1)/2) then
                qx = fx*(i - nx - 1)
            else
                qx = fx*(i - 1)
            end if

            do j = 1, ny
                if (j > (ny + 1)/2) then
                    qy = fy*(j - ny - 1)
                else
                    qy = fy*(j - 1)
                end if
                !                qy=qy*qy
                !                if(qy.ne.0.d0) qy = 2.d0/ay**2*(1.d0-dcos(qy*ay))
                !             qy = dsin(qy*ay)/ay
                do k = 1, nz
                    if (k > (nz + 1)/2) then
                        qz = fz*(k - nz - 1)
                    else
                        qz = fz*(k - 1)
                    end if

                    !                   if(qz.ne.0.d0) qz = 2.d0/az**2*(1.d0-dcos(qz*az))
                    !                   q2 = qx + qy + qz*qz
                    iq = i + nx*(j - 1) + nx*ny*(k - 1)
                    fp(iq) = fp(iq)*dcmplx(0.d0, qy)
                end do
            end do
        end do

        call cfft3d(fp, nx, ny, nz, nx, ny, nz, +1)

        gradt(2, :) = dreal(fp(:))

        !  z-grad

        fp = dent

        call cfft3d(fp, nx, ny, nz, nx, ny, nz, -1)

        do i = 1, nx
            if (i > (nx + 1)/2) then
                qx = fx*(i - nx - 1)
            else
                qx = fx*(i - 1)
            end if

            do j = 1, ny
                if (j > (ny + 1)/2) then
                    qy = fy*(j - ny - 1)
                else
                    qy = fy*(j - 1)
                end if
                !                qy=qy*qy
                !                if(qy.ne.0.d0) qy = 2.d0/ay**2*(1.d0-dcos(qy*ay))
                do k = 1, nz
                    if (k > (nz + 1)/2) then
                        qz = fz*(k - nz - 1)
                    else
                        qz = fz*(k - 1)
                    end if

                    !                   if(qz.ne.0.d0) qz = 2.d0/az**2*(1.d0-dcos(qz*az))
                    !                   q2 = qx + qy + qz*qz
                    !                   qz = dsin(qz*az)/az
                    iq = i + nx*(j - 1) + nx*ny*(k - 1)
                    fp(iq) = fp(iq)*dcmplx(0.d0, qz)
                end do
            end do
        end do

        call cfft3d(fp, nx, ny, nz, nx, ny, nz, +1)

        gradt(3, :) = dreal(fp(:))

#endif

        return

    end subroutine evalgrad

    subroutine upload_vhartree_fromfp
        implicit none
        integer i, i1, j, k, ii, jj, kk, ind, indx, indmesh, indtot, indproc, nbufrep, ierr&
                &, indproc_vh, indsi, indsj, indsk, indsiu, indsju, indsku
        integer indi, indj, indk, indiu, nionu, ind_hartree, ind_upper
        integer, dimension(:), allocatable :: indju, indku
        integer indtot_proc, nxupo
        integer*8 indtot_mesh, i8
        logical, external :: condp
        !logical, dimension(:), allocatable:: yestouched
#ifdef PARALLEL
        include 'mpif.h'
#endif
        !      allocate(yestouched(meshproc_tot))

        !      yestouched=.false.
        if (iespbc) then
            nxupo = nxu
        else
            nxupo = nxu*2
        end if
        allocate (buffer_fft(dim_fft_vh, nproc_fft), buffer_local(dim_fft_vh))
        buffer_fft = 0.d0
        buffer_local = 0.d0
        indmesh = 0
        indtot = 0
        nbufrep = nproc_fft*dim_fft_vh
        ind = 0
        indproc_vh = 0
        indi = 0
        indj = 1
        indk = 1
        indproc = 0

        do k = 1, nzu
            do j = 1, nyu
                indx = 0
                do i = 1, nxupo
                    indtot = indtot + 1
                    if (sub_comm_fft%yesin) then
                        if (indproc .eq. sub_comm_fft%rank) then
                            ind = ind + 1
                            indx = indx + 1
                            if (iespbc) then
                                buffer_local(ind) = fp(k + nzu*(indx - 1) + nzu*nxl*(j - 1))
                            else
                                buffer_local(ind) = fp(k + 2*nzu*(indx - 1) + 2*nzu*nxl2*(j - 1))
                            end if
                        end if
                    end if
                    indproc = indproc + 1
                    if (indproc .eq. nproc_fft) indproc = 0
                    if (indtot .eq. nbufrep .or. (i .eq. nxupo .and. j .eq. nyu .and. k .eq. nzu)) then
#ifdef PARALLEL
                        if (sub_comm_fft%yesin) then
                            call mpi_gather(buffer_local, dim_fft_vh, MPI_DOUBLE_PRECISION, buffer_fft&
                            &, dim_fft_vh, MPI_DOUBLE_PRECISION, 0, sub_comm_fft%comm, ierr)
                        end if
                        call bcast_real(buffer_fft, nbufrep, 0, commrep_mpi)
#else
                        buffer_fft(:, 1) = buffer_local(:)
#endif
                        do ii = 1, dim_fft_vh
                            do jj = 1, nproc_fft
                                ! Update counter i,j,k in the right sequence
                                indi = indi + 1
                                if (indi .gt. nxupo) then
                                    indi = 1
                                    indj = indj + 1
                                    if (indj .gt. nyu) then
                                        indj = 1
                                        indk = indk + 1
                                    end if
                                end if
                                if (indk .le. nzu .and. indi .le. nxu) then
                                    ! a point of the sparse lattice as it is in the same order indmesh should work
                                    if (do_hartree) then
                                        ! To be done
                                    else

                                        if (indproc_vh .eq. rankrep) then
                                            indmesh = indmesh + 1
                                            vhartree(indmesh) = buffer_fft(ii, jj)
                                            !               yestouched(indmesh)=.true.
                                        end if
                                        indproc_vh = indproc_vh + 1
                                        if (indproc_vh .eq. nprocrep) indproc_vh = 0
                                    end if
                                end if ! endif indk.le.nzu
                            end do
                        end do
                        buffer_local = 0.d0
                        ind = 0
                        indtot = 0
                    end if
                end do
            end do
        end do

        deallocate (buffer_local, buffer_fft)

    end subroutine upload_vhartree_fromfp

    subroutine upload_fp_fromdent
        use setup, only: dent, do_hartree
        use allio, only: zetar
        use constants, only: zzero
        implicit none
        integer i, j, k, ii, jj, kk, ind, indmesh, indtot, indproc, nbufrep, indi, indj, indk&
                &, indi_proc, ierr, nxi, nxf, nyi, nyf, nzi, nzf, nx0, ny0, nz0, i_fit, j_fit, k_fit&
                &, indfp, i1, i2, i3, ind1, ind2, ind3, scalea, nxr, nyr, nzr, nx0n, ny0n, nz0n
        real*8 px, py, pz, pxm, pym, pzm, ppx(4), ppy(4), ppz(4), axn, ayn, azn
        real*8, dimension(:, :), allocatable :: buffer_fft
        real*8, dimension(:), allocatable :: buffer_local
#ifdef PARALLEL
        include 'mpif.h'
#endif
        fp = zzero
        if (do_hartree) then
            ! To be done
        else ! do_hartree
            allocate (buffer_fft(dim_fft_den, nprocrep), buffer_local(dim_fft_den))
            buffer_fft = 0.d0
            buffer_local = 0.d0
            indproc = 0
            indmesh = 0
            indtot = 0
            nbufrep = nprocrep*dim_fft_den
            indi = 0
            indj = 1
            indk = 1
            ind = 0

            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        indtot = indtot + 1
                        if (indproc .eq. rankrep) then
                            indmesh = indmesh + 1
                            ind = ind + 1
                            buffer_local(ind) = dent(indmesh)
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0
                        if (indtot .eq. nbufrep .or. (i .eq. nx .and. j .eq. ny .and. k .eq. nz)) then
#ifdef PARALLEL
                            call mpi_allgather(buffer_local, dim_fft_den, MPI_DOUBLE_PRECISION, buffer_fft&
                            &, dim_fft_den, MPI_DOUBLE_PRECISION, commrep_mpi, ierr)
#else
                            buffer_fft(:, 1) = buffer_local(:)
#endif
                            do ii = 1, dim_fft_den
                                do jj = 1, nprocrep
                                    ! Update counter i,j,k in the right sequence
                                    indi = indi + 1
                                    if (indi .gt. nx) then
                                        indi = 1
                                        indj = indj + 1
                                        if (indj .gt. ny) then
                                            indj = 1
                                            indk = indk + 1
                                        end if
                                    end if
                                    if (sub_comm_fft%yesin) then
                                        if (mod(indi - 1, nproc_fft) .eq. &
                                            sub_comm_fft%rank .and. indi .le. nx .and. indj .le. ny .and. indk .le. nz) then
                                            indi_proc = (indi - 1)/nproc_fft + 1
                                            if (iespbc) then
                                                fp(indk + nzu*(indi_proc - 1) + nzu*nxl*(indj - 1)) = buffer_fft(ii, jj)
                                            else
                                                fp(indk + 2*nzu*(indi_proc - 1) + 2*nzu*nxl2*(indj - 1)) = buffer_fft(ii, jj)
                                            end if
                                        end if
                                    end if
                                end do
                            end do
                            buffer_local = 0.d0
                            ind = 0
                            indtot = 0
                        end if
                    end do
                end do
            end do
            deallocate (buffer_local, buffer_fft)
        end if
    end subroutine upload_fp_fromdent

end module fourier_module
subroutine cubic(x, p)
    implicit none
    real*8 x ! input x in unit of lattice space a   0<=x<=1
    !      xfit=x2+x a
    !      input data y1 y2 y3 y4
    !      NB for x=0 p(2)=1 and all  the other are zero. This will be used
    !      to avoid to compute a mesh point that is not sampled.
    real*8 p(4), x1, x2, x3, x4
    x1 = -1.d0
    x2 = 0.d0
    x3 = 1.d0
    x4 = 2.d0
    p(1) = (x - x2)/(x1 - x2)*(x - x3)/(x1 - x3)*(x - x4)/(x1 - x4)
    p(2) = (x - x1)/(x2 - x1)*(x - x3)/(x2 - x3)*(x - x4)/(x2 - x4)
    p(3) = (x - x1)/(x3 - x1)*(x - x2)/(x3 - x2)*(x - x4)/(x3 - x4)
    p(4) = (x - x1)/(x4 - x1)*(x - x2)/(x4 - x2)*(x - x3)/(x4 - x3)
    !      By Lagrange polynomial
    !      cubic= y(1)*(x-x2)/(x1-x2)*(x-x3)/(x1-x3)*(x-x4)/(x1-x4)+&
    !           & y(2)*(x-x1)/(x2-x1)*(x-x3)/(x2-x3)*(x-x4)/(x2-x4)+&
    !           & y(3)*(x-x1)/(x3-x1)*(x-x2)/(x3-x2)*(x-x4)/(x3-x4)+&
    !           & y(4)*(x-x1)/(x4-x1)*(x-x2)/(x4-x2)*(x-x3)/(x4-x3)
    return
end

function attenuate(x)
    use constants, only: Pi
    use setup, only: linear
    implicit none
    real*8 attenuate, x, z
    z = Pi*x
    if (linear) then
        !       linear fit
        if (abs(x) .gt. 1d-8) then
            attenuate = (sin(z)/z)**2
        else
            attenuate = 1.d0
        end if
    elseif (.not. linear) then
        !       cubic fit
        if (abs(x) .gt. 1d-8) then
            attenuate = (sin(z)/z)**4*(1.d0 + 0.6666666666667d0*z**2)
        else
            attenuate = 1.d0
        end if
    else
        !      Highest polynomial fit
        if (abs(x) .gt. 1d-8) then
            attenuate = (sin(z)/z)**8*(1.d0 + 1.33333333333333d0*z**2 + .933333333333333d0*z**4 + .4571428571428571d0*z**6)
        else
            attenuate = 1.d0
        end if
    end if
    return
end
function mod_true(i, n)
    implicit none
    integer mod_true, i, n
    mod_true = mod(i, n)
    if (mod_true .lt. 0) mod_true = mod_true + n
    return
end
function condp(i, mini, maxi, iespbc)
    logical condp, iespbc
    integer i, mini, maxi
    if (iespbc) then
        if (mini .le. maxi) then
            if (i .le. maxi .and. i .ge. mini) then
                condp = .true.
            else
                condp = .false.
            end if
        else
            if (i .le. maxi .or. i .ge. mini) then
                condp = .true.
            else
                condp = .false.
            end if
        end if
    else
        if (i .le. maxi .and. i .ge. mini) then
            condp = .true.
        else
            condp = .false.
        end if
    end if
    return
end

