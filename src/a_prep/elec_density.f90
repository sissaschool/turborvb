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

subroutine updenorb_new

    ! This subroutine updates charge and spin density
    ! after diagonalization of KS hamiltonian.
    ! Compute charge (dent) and spin (spint) density
    ! on the mesh from occupations of DFT orbitals and
    ! k-points weights if using k-points.
    ! This subroutine is general for both real and complex
    ! algorithms.

    use kpoints_mod, only: wkp, wkp_down, sum_kpoints_array_real8, decoupled_run
    use constants, only: ipc, zzero, zone
    use cell, only: cellscale, map, metric, car2cry, at, unit_volume, yes_tilted
    use symmetries, only: isymm, nrot
    use parallel_module, only: old_threads
    use buffers

    use allio, only: nshell, ioptorb, ioccup, dupr, zetar, rank, &
                     rion, psip, nion, kion, iflagnorm, cnorm, indt, &
                     LBox, indpar_tab, indorb_tab, &
                     indshell_tab, nel, nprocrep, rankrep, commrep_mpi, &
                     iespbc, nelorbh, nelorb, costz, costz3, niesd, mu_c, &
                     vj, n_body_on, scale_one_body, same_phase, &
                     skip_equivalence, commcolrep_mpi, pointvj, norm_metric

    use setup, only: molecorb, molecorbdo, wf, occupations, occupationdo

    use setup, only: dent, spint, bufbuf, nelorb3, memlarge, nx, ny, nz, ax, ay, &
            zgemm_time, zgemm_timep, dgemm_time, dgemm_timep, az, &
            yeslsda, bands, meshproc, indk, contracted_on, rion_ref, &
            bands, bandsdo, nelorbu, mindist, double_overs, &
            buffer_grid, fix_density, volmesh, double_mesh, scale_z, &
            & l0_at, nx0, ny0, nz0, weightx, weighty, weightz, minz_at&
            &, from_ions, rion_from, nx_at, ny_at, nz_at, rion_upload

    implicit none
    ! local
    integer :: i, j, k, ii, jj, kk, indmax, indmesh, ibuf, nbas_1, nbas_tot &
               , buf_dim, ierr, ind, wf_true_dim, indproc, indtot, nbufrep &
               , nxr, nyr, nzr, nxi, nxf, nyi, nyf, nzi, nzf, scalea, nx0n &
               , ny0n, nz0n
    real(8) :: dens, axn, ayn, azn
    real(8), external :: cclock
    integer :: tid, thread_active

    ! local buffers
    real(8), dimension(:, :), allocatable :: distp_scratch, r_scratch
    real(8), dimension(:, :), allocatable :: rmu_scratch, rmusin_scratch, rmucos_scratch
    real(8), dimension(:, :), allocatable :: x_buffer

#ifdef PARALLEL
    include "mpif.h"
#endif

#if defined(_OPENMP)
    integer, external :: omp_get_max_threads, omp_get_thread_num
    thread_active = omp_get_max_threads()
#if defined(__NOOMPDFT)
    if (.not. memlarge) then
        call omp_set_num_threads(1) ! scalar code
        thread_active = 1
    end if
#endif
#else
    thread_active = 1
#endif

    ! initialize variables
    ind = 0
    indmesh = 0
    ! define indices of the buffers
    nbas_1 = nelorbu ! # of basis elements (contracted or uncontracted)
    nbas_tot = nelorb ! # of uncontracted basis elements >= nbas_1
    buf_dim = nelorb3 ! leading dimension of the buffer
    ! define index for the second
    ! part of buffer/buffer_on used as scratch
    if (ipc .eq. 1) then
        ibuf = nbas_1 + 1
        wf_true_dim = nbas_1
    else
        ibuf = 2*nbas_1 + 1
        if (double_overs) then
            wf_true_dim = 2*nbas_1
        else
            wf_true_dim = nbas_1
        end if
    end if

    ! allocate scratch vectors and buffers
    allocate (distp_scratch(nbas_tot*(indt + 5) + 27*(indt + 1)*max(nshell, nion), thread_active))
    allocate (r_scratch(nion, thread_active))
    allocate (rmu_scratch(3*nion, thread_active))
    allocate (rmusin_scratch(3*nion, thread_active))
    allocate (rmucos_scratch(3*nion, thread_active))
    allocate (x_buffer(3, bufbuf))
    distp_scratch = 0.d0
    r_scratch = 0.d0
    rmu_scratch = 0.d0
    rmusin_scratch = 0.d0
    rmucos_scratch = 0.d0
    call allocate_buffers(bufbuf, nbas_tot, buf_dim, nbas_tot, nbas_tot, nbas_tot, thread_active)

    nbufrep = nprocrep*bufbuf
    indmesh = 0
    dent = 0.d0
    if (yeslsda) spint = 0.d0

    call upload_dens(1, nx, 1, ny, 1, nz, ax, ay, az, .true.)

    if (double_mesh) then

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
                    call set_interval(scalea, nx0, ny0, nz0, rion(1, ii), rion_upload, ax, ay, az &
                            , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .true.&
                            &, weightx, weighty, weightz)
                    deallocate (weightx, weighty, weightz)

                    call upload_dens(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, .false.)

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

                    call set_interval(scalea, nx0n, ny0n, nz0n, rion(1, ii), rion_upload, axn, ayn, azn &
                            , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .false.&
                            &, weightx, weighty, weightz)
                    deallocate (weightx, weighty, weightz)

                    call upload_dens(nxi, nxf, nyi, nyf, nzi, nzf, axn, ayn, azn, .true.)
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
            call set_interval(scalea, nx0, ny0, nz0, rion_from, rion_upload, ax, ay, az &
                    , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .true.&
                    &, weightx, weighty, weightz)
            deallocate (weightx, weighty, weightz)

            call upload_dens(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, .false.)

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

            call set_interval(scalea, nx0n, ny0n, nz0n, rion_from, rion_upload, axn, ayn, azn &
                    , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .false.&
                    &, weightx, weighty, weightz)
            deallocate (weightx, weighty, weightz)

            call upload_dens(nxi, nxf, nyi, nyf, nzi, nzf, axn, ayn, azn, .true.)
        end if
        !   restore previous value of volmesh

        volmesh = ax*ay*az*unit_volume

    end if

    !
    ! deallocate buffers and scratch vectors
    !
    call deallocate_buffers()
    !
    deallocate (distp_scratch)
    deallocate (r_scratch, rmu_scratch, rmusin_scratch, rmucos_scratch)
    deallocate (x_buffer)

#if defined (_OPENMP) && defined (__NOOMPDFT)
    if (.not. memlarge) call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

#ifdef PARALLEL
    call mpi_allreduce(ind, indmax, 1, MPI_INTEGER, MPI_MAX, commrep_mpi, ierr)
#else
    indmax = ind
#endif

    if (indmax .ne. 0) &
        call error(' updenorb ', ' check input nbufd and/or code ', 1, rank)

    ! Do not average the electronic density over k-points in the case of decoupled
    ! calculations, except if the flag fix_density is active. In this case the k-points
    ! are evolved independently using the averaged density.
    if (.not. decoupled_run .or. fix_density) then

        ! collect the density from the pools
        call sum_kpoints_array_real8(dent, meshproc, commcolrep_mpi, -1)
        if (yeslsda) then
            call sum_kpoints_array_real8(spint, meshproc, commcolrep_mpi, -1)
        end if

        ! symmetrize density using the Bravais lattice symmetries
        if (.not. skip_equivalence .and. .not. yes_tilted) then
            if (.not. allocated(buffer_grid)) then
                allocate (buffer_grid(meshproc, nprocrep))
                buffer_grid = 0.d0
            end if
            call symmetrize_density(dent, buffer_grid, meshproc, nx, ny, nz, isymm, nrot)
            if (yeslsda) call symmetrize_density(spint, buffer_grid, meshproc, nx, ny, nz, isymm, nrot)
            deallocate (buffer_grid)
        end if

    end if

#ifdef PARALLEL
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier

    return

contains

    subroutine upload_dens(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, add)
        implicit none
        integer nxi, ii, nxf, nyi, nyf, nzi, nzf, nx, ny, nz, indtot, bufmax
        logical add
        real*8 ax, ay, az

        if (memlarge) then
            !

            indproc = 0
            indtot = 0
            ind = 0

            !
            do k = nzi, nzf
                do j = nyi, nyf
                    do i = nxi, nxf
                        indtot = indtot + 1

                        if (indproc .eq. rankrep) then
                            indmesh = indmesh + 1
                            ind = ind + 1
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0

                        if (indtot .eq. nbufrep .or. (i .eq. nxf .and. j .eq. nyf .and. k .eq. nzf)) then
#ifdef _OFFLOAD
!$omp target data map(buffer(:,1:ind)) &
!$omp& map(to:molecorb(:,1:bands),wf(:,indmesh-ind+1:indmesh))
#endif
                            if (ipc .eq. 1) then

                                dgemm_timep = cclock()
                                call dgemm_('T', 'N', bands, ind, nbas_1, 1.d0, molecorb, nbas_1, &
                                            wf(1, indmesh - ind + 1), wf_true_dim, 0.d0, buffer, buf_dim)
                                dgemm_time = dgemm_time + cclock() - dgemm_timep

                            else

                                zgemm_timep = cclock()
                                call zgemm_('T', 'N', bands, ind, nbas_1, zone, molecorb, nbas_1, &
                                            wf(1, indmesh - ind + 1), wf_true_dim, zzero, buffer, buf_dim)
                                zgemm_time = zgemm_time + cclock() - zgemm_timep

                            end if
#ifdef _OFFLOAD
!$omp end target data
#endif

!$omp parallel do default(shared) private(ii,jj,dens) schedule(static)
                            do ii = 1, ind
                                dens = 0.d0
                                do jj = 1, bands
                                    if (ipc .eq. 1) then
                                        dens = dens + occupations(jj)*buffer(jj, ii)**2
                                    else
                                        dens = dens + occupations(jj)* &
                                               ((buffer(2*jj - 1, ii))**2 + (buffer(2*jj, ii))**2)
                                    end if
                                end do
                                dent(indmesh - ind + ii) = dens
                            end do
!$omp end parallel do

                            if (yeslsda .or. ipc .eq. 2) then
                                ! in the case of different boundaries for up/down spin or for LSDA calculations
                                ! the basis set is stored in the array: wf(nelorb+1:2*nelorb,:)
                                if (ipc .eq. 1) then
#ifdef _OFFLOAD
!$omp target data map(buffer(:,1:ind)) &
!$omp& map(to:molecorbdo(:,1:bands),wf(:,indmesh-ind+1:indmesh))
#endif

                                    dgemm_timep = cclock()
                                    call dgemm_('T', 'N', bands, ind, nbas_1, 1.d0, molecorbdo, nbas_1&
                                            &, wf(1, indmesh - ind + 1), wf_true_dim, 0.d0, buffer, buf_dim)
                                    dgemm_time = zgemm_time + cclock() - dgemm_timep
#ifdef _OFFLOAD
!$omp end target data
#endif

                                else

                                    zgemm_timep = cclock()
                                    if (double_overs) then
#ifdef _OFFLOAD
!$omp target data map(buffer(:,1:ind)) &
!$omp& map(to:molecorbdo(:,1:bands),wf(:,indmesh-ind+1:indmesh))
#endif
                                        call zgemm_('T', 'N', bands, ind, nbas_1, zone, molecorbdo, nbas_1&
                                                &, wf(2*nbas_1 + 1, indmesh - ind + 1), wf_true_dim, zzero, buffer, buf_dim)
#ifdef _OFFLOAD
!$omp end target data
#endif
                                    else
                                        ! if same phase the buffer does not change
                                        ! if opposite phase, use the complex conjugate of the wf array computed before
                                        call fill_phase_wf(molecorbdo, nbas_1, wf(1, indmesh - ind + 1) &
                                                           , wf_true_dim, buffer, buf_dim, same_phase, ind)
                                    end if
                                    zgemm_time = zgemm_time + cclock() - zgemm_timep

                                end if

!$omp parallel do default(shared) private(ii,jj,dens) schedule(static)
                                do ii = 1, ind
                                    dens = 0.d0
                                    do jj = 1, bands
                                        if (ipc .eq. 1) then
                                            dens = dens + occupationdo(jj)*buffer(jj, ii)**2
                                        else
                                            dens = dens + occupationdo(jj)*((buffer(2*jj - 1, ii))**2 + &
                                                                            (buffer(2*jj, ii))**2)
                                        end if
                                    end do
                                    ! spin density = density_up - density_down
                                    ! compute the spin density only in the case of a TRUE LSDA functional
                                    ! and not for the "fake" LSDA used when there's a complex WF with LDA.
                                    if (yeslsda) spint(indmesh - ind + ii) = dent(indmesh - ind + ii) - dens
                                    dent(indmesh - ind + ii) = dent(indmesh - ind + ii) + dens
                                end do
!$omp end parallel do

                            end if
                            ind = 0
                            indtot = 0
                        end if
                    end do
                end do
            end do

        else ! memlarge=.false. (default value)
            !
            ind = 0
            indproc = 0
            indtot = 0
            !
            do k = nzi, nzf
                do j = nyi, nyf
                    do i = nxi, nxf
                        ! grid point buffer counter
                        indtot = indtot + 1
                        if (indproc .eq. rankrep) then
                            indmesh = indmesh + 1
                            ind = ind + 1
                            ! buffering the grid point coordinates
                            x_buffer(:, ind) = i*ax*at(:, 1) + j*ay*at(:, 2) + k*az*at(:, 3) + rion_upload(:)
                            !               x_buffer(:,ind)=i*ax+rion_upload(1)
                            !               x_buffer(2,ind)=j*ay+rion_upload(2)
                            !               x_buffer(3,ind)=k*az+rion_upload(3)
                        end if
                        indproc = indproc + 1
                        if (indproc .eq. nprocrep) indproc = 0

                        ! exceeding buffer dimension: buffering
                        if (indtot .eq. nbufrep .or. (i .eq. nxf .and. j .eq. nyf .and. k .eq. nzf)) then
!$omp parallel default(shared) private(tid,ii)
#if defined(_OPENMP)
                            tid = 1 + omp_get_thread_num()
#else
                            tid = 1
#endif

!$omp do
                            do ii = 1, ind
                                call compute_one_grid_point(x_buffer(1, ii), tid, ii, indmesh - ind + ii)
                            end do
!$omp end do nowait
!$omp end parallel
#if defined (_OPENMP) && defined (__NOOMPDFT)
                            call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
#ifdef _OFFLOAD
!$omp target data map(from:buffer(:,1:ind)) &
!$omp& map(to:mu_c(:,1:nbas_1),molecorb(:,1:bands)&
!$omp& ,buffer_on(:,1:ind)) if(contracted_on)
!$omp target data map(buffer(:,1:ind)) &
!$omp& map(to:molecorb(:,1:bands)) if(.not.contracted_on)
#endif
                            if (ipc .eq. 1) then
                                dgemm_timep = cclock()
                                if (contracted_on) then
                                    call dgemm_('T', 'N', nbas_1, ind, nelorbh, 1.d0, mu_c, nelorbh, &
                                                buffer_on, nbas_tot, 0.d0, buffer, buf_dim)
                                end if
                                call dgemm_('T', 'N', bands, ind, nbas_1, 1.d0, molecorb, nbas_1, &
                                            buffer, buf_dim, 0.d0, buffer(ibuf, 1), buf_dim)
                                dgemm_time = dgemm_time + cclock() - dgemm_timep
                            else
                                zgemm_timep = cclock()
                                if (contracted_on) then
                                    call zgemm_('T', 'N', nbas_1, ind, nelorbh, zone, mu_c, nelorbh, &
                                                buffer_on, nbas_tot, zzero, buffer, buf_dim)
                                end if
                                call zgemm_('T', 'N', bands, ind, nbas_1, zone, molecorb, nbas_1, &
                                            buffer, buf_dim, zzero, buffer(ibuf, 1), buf_dim)
                                zgemm_time = zgemm_time + cclock() - zgemm_timep
                            end if
#ifdef  _OFFLOAD
!$omp end target data
!$omp end target data
#endif
!$omp parallel do default(shared) private(ii,jj,dens) schedule(static)
                            do ii = 1, ind
                                dens = 0.d0
                                do jj = 1, bands
                                    if (ipc .eq. 1) then
                                        dens = dens + occupations(jj)*buffer(nbas_1 + jj, ii)**2
                                    else
                                        dens = dens + occupations(jj)* &
                                               ((buffer(2*(nbas_1 + jj) - 1, ii))**2 + (buffer(2*(nbas_1 + jj), ii))**2)
                                    end if
                                end do
                                dent(indmesh - ind + ii) = dens
                            end do
!$omp end parallel do

                            if (yeslsda .or. ipc .eq. 2) then ! in case of spin calculations or complex wave function
#ifdef _OFFLOAD
                                if (ipc .eq. 2 .and. double_overs) then
                                    bufmax = bufbuf + ind
                                else
                                    bufmax = ind
                                end if
!$omp target data map(from:buffer(:,1:bufmax))&
!$omp& map(to:mu_c(:,1:nbas_1),molecorbdo(:,1:bands)&
!$omp& ,buffer_on(:,1:bufmax)) if(contracted_on.and.ipc.eq.2)
!$omp target data map(buffer(:,1:bufmax))&
!$omp& map(to:molecorbdo(:,1:bands)) if(.not.contracted_on.or.ipc.eq.1)
#endif

                                if (ipc .eq. 1) then
                                    dgemm_timep = cclock()
                                    call dgemm_('T', 'N', bands, ind, nbas_1, 1.d0, molecorbdo, nbas_1, &
                                                buffer, buf_dim, 0.d0, buffer(ibuf, 1), buf_dim)
                                    dgemm_time = dgemm_time + cclock() - dgemm_timep
                                else
                                    zgemm_timep = cclock()

                                    if (double_overs) then
                                        if (contracted_on) then
                                            call zgemm_('T', 'N', nbas_1, ind, nelorbh, zone, mu_c, nelorbh, &
                                                        buffer_on(1, bufbuf + 1), nbas_tot, zzero, buffer(1, bufbuf + 1), buf_dim)
                                        end if
                                        call zgemm_('T', 'N', bands, ind, nbas_1, zone, molecorbdo, nbas_1, &
                                                    buffer(1, bufbuf + 1), buf_dim, zzero, buffer(ibuf, 1), buf_dim)
                                    else
                                        if (contracted_on) then
                                            if (.not. same_phase) &
                                                call conjmat_(nbas_tot, ind, buffer_on, nbas_tot)
                                            call zgemm_('T', 'N', nbas_1, ind, nelorbh, zone, mu_c, nelorbh, &
                                                        buffer_on, nbas_tot, zzero, buffer, buf_dim)
                                        else
                                            if (.not. same_phase) &
                                                call conjmat_(buf_dim, ind, buffer, buf_dim)
                                        end if
                                        call zgemm_('T', 'N', bands, ind, nbas_1, zone, molecorbdo, nbas_1, &
                                                    buffer, buf_dim, zzero, buffer(ibuf, 1), buf_dim)
                                    end if
                                    zgemm_time = zgemm_time + cclock() - zgemm_timep

                                end if
#ifdef _OFFLOAD
!$omp end target data
!$omp end target data
#endif
!$omp parallel do default(shared) private(ii,jj,dens) schedule(static)
                                do ii = 1, ind
                                    dens = 0.d0
                                    do jj = 1, bands
                                        if (ipc .eq. 1) then
                                            dens = dens + occupationdo(jj)*buffer(nbas_1 + jj, ii)**2
                                        else
                                            dens = dens + occupationdo(jj)* &
                                                   ((buffer(2*(nbas_1 + jj) - 1, ii))**2 + (buffer(2*(nbas_1 + jj), ii))**2)
                                        end if
                                    end do

                                    !  write(6,*) 'density:',rank,dens,sum(buffer(:,:))
                                    !  call mpi_finalize(ierr)
                                    !  stop
                                    ! compute the spin density only in the case of a TRUE LSDA functional
                                    ! and not for the "fake" LSDA used when there's a complex wf.
                                    if (yeslsda) spint(indmesh - ind + ii) = dent(indmesh - ind + ii) - dens
                                    dent(indmesh - ind + ii) = dent(indmesh - ind + ii) + dens
                                end do
!$omp end parallel do
                            end if

#if defined (_OPENMP) && defined (__NOOMPDFT)
                            call omp_set_num_threads(1) ! restore the scalar code
#endif
                            ind = 0
                            indtot = 0
                        end if
                    end do
                end do
            end do
        end if
    end subroutine upload_dens

    subroutine compute_one_grid_point(x, tid, ind, indmesh)
        ! it is thread-safe and intended to be called from threaded region.
        ! computing the wave function for up/down spin electrons
        implicit none
        ! arguments
        ! single electron positions
        real(8), intent(in) :: x(3)
        ! index to the scratch space
        integer, intent(in) :: tid
        ! index in the grid buffer
        integer, intent(in) :: ind
        ! global index of the current grid point
        integer, intent(in) :: indmesh
        !
        ! local variables
        real(8) :: jas_1body, r0, rc(3)
        real(8), external :: jastrow_ei
        integer :: jj, kk

        ! computing the wave function for up/down spin electrons
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, x, 1, r_scratch(1, tid), rmu_scratch(1, tid), &
                     dupr, zetar, rion, distp_scratch(1, tid), wf_threading_scratch(1, tid), nbas_tot &
                     , nion, kion, iflagnorm, cnorm, LBox, rmucos_scratch(1, tid), rmusin_scratch(1, tid) &
                     , mindist, indpar_tab, indorb_tab, indshell_tab, .true.)
        if (double_overs) then
            call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, x, 1, r_scratch(1, tid), rmu_scratch(1, tid), &
                         dupr, zetar, rion, distp_scratch(1, tid), wf_threading_scratch_down(1, tid), nbas_tot &
                         , nion, kion, iflagnorm, cnorm, LBox, rmucos_scratch(1, tid), rmusin_scratch(1, tid) &
                         , mindist, indpar_tab, indorb_tab, indshell_tab, .false.)
        end if

        ! adding 1b Jastrow if any
        if (n_body_on .ne. 0) then
            jas_1body = -scale_one_body
            do jj = 1, nion
                if (iespbc) then
                    rc(1) = x(1) - rion(1, jj)
                    rc(2) = x(2) - rion(2, jj)
                    rc(3) = x(3) - rion(3, jj)
!                   call CartesianToCrystal(rc, 1)
                    rc(:) = car2cry(:, 1)*rc(1) + car2cry(:, 2)*rc(2) + car2cry(:, 3)*rc(3)

                    do kk = 1, 3
                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                    end do
                    r0 = norm_metric(rc, metric)
                else
                    rc(1) = (x(1) - rion(1, jj))*costz(jj)
                    rc(2) = (x(2) - rion(2, jj))*costz(jj)
                    rc(3) = (x(3) - rion(3, jj))*costz(jj)
                    r0 = dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2)
                end if
                jas_1body = jas_1body - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
            end do
            jas_1body = dexp(jas_1body)
            wf_threading_scratch(1:ipc*nbas_tot, tid) = wf_threading_scratch(1:ipc*nbas_tot, tid)*jas_1body
            if (double_overs) &
                wf_threading_scratch_down(1:ipc*nbas_tot, tid) = wf_threading_scratch_down(1:ipc*nbas_tot, tid)*jas_1body
        end if

        ! load the buffer
        if (contracted_on) then
            buffer_on(1:ipc*nbas_tot, ind) = wf_threading_scratch(1:ipc*nbas_tot, tid)
            if (double_overs) then
                buffer_on(1:ipc*nbas_tot, bufbuf + ind) = wf_threading_scratch_down(1:ipc*nbas_tot, tid)
            end if
        else
            buffer(1:ipc*nbas_tot, ind) = wf_threading_scratch(1:ipc*nbas_tot, tid)
            if (double_overs) then
                buffer(1:ipc*nbas_tot, bufbuf + ind) = wf_threading_scratch_down(1:ipc*nbas_tot, tid)
            end if
        end if
    end subroutine compute_one_grid_point

end subroutine updenorb_new

! This subrutine integrates the charge/spin
! over the whole mesh and computes the total charge.
! It is called at the end of each DFT cycle to check
! if eigenvectors are correct: /int_mesh n(r) dr = nel

subroutine cutdens

    use allio, only: nelup, neldo, nel, rank, commrep_mpi, rankcolrep, rankrep

    use setup, only: dent, spint, volmesh, meshproc, meshproc_tot, yeslsda, &
                     denstot, spintot, spin2tot, spingrid, voltot, &
                     occupations, occupationdo, dens0, double_mesh, volmesh_proc, &
                     scale_hartreen, corr_hartree, scale_hartree

    use fourier_module, only: vhartree

    implicit none

    integer i, ierr
    real(8) cost, densproc, densproc_cut, minden, mindenproc, spinproc, spin2proc, &
        minspinproc, minspin, spingridproc, vhq0, denstot_cut

#ifdef PARALLEL
    include 'mpif.h'
#endif

    mindenproc = dent(1)
    do i = 2, meshproc_tot
        if (dent(i) .lt. mindenproc) mindenproc = dent(i)
    end do
    if (double_mesh) then
        densproc_cut = sum(dent(1:meshproc)*volmesh_proc(1:meshproc))
        densproc = sum(dent(1:meshproc_tot)*volmesh_proc(1:meshproc_tot))
        vhq0 = sum(vhartree(1:meshproc_tot)*volmesh_proc(1:meshproc_tot))
    else
        densproc = sum(dent(1:meshproc))*volmesh
    end if

    if (yeslsda) then
        minspinproc = dent(1) - abs(spint(1))
        do i = 2, meshproc_tot
            cost = dent(i) - abs(spint(i))
            if (cost .lt. minspinproc) minspinproc = cost
        end do
        if (double_mesh) then
            spinproc = sum(spint(1:meshproc_tot)*volmesh_proc(1:meshproc_tot))
            spin2proc = sum(spint(1:meshproc_tot)**2*volmesh_proc(1:meshproc_tot))
            spingridproc = sum(abs(spint(1:meshproc_tot))*volmesh_proc(1:meshproc_tot))

        else
            spinproc = sum(spint(1:meshproc))*volmesh
            spin2proc = sum(spint(1:meshproc)**2)*volmesh
            spingridproc = 0.d0
            do i = 1, meshproc
                spingridproc = spingridproc + volmesh*abs(spint(i))
            end do
        end if

    else
        minspinproc = 0.d0
        spinproc = 0.d0
        spin2proc = 0.d0
        spingridproc = 0.d0
    end if
#ifdef PARALLEL
    call reduce_base_real_to(1, densproc, denstot, commrep_mpi, -1)
    call reduce_base_real(1, densproc_cut, commrep_mpi, -1)
    call reduce_base_real(1, vhq0, commrep_mpi, -1)
    if (yeslsda) then
        call reduce_base_real_to(1, spinproc, spintot, commrep_mpi, -1)
        call reduce_base_real_to(1, spin2proc, spin2tot, commrep_mpi, -1)
        call reduce_base_real_to(1, spingridproc, spingrid, commrep_mpi, -1)
        call mpi_allreduce(minspinproc, minspin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, commrep_mpi, ierr)
    end if
    call mpi_allreduce(mindenproc, minden, 1, MPI_DOUBLE_PRECISION, MPI_MIN, commrep_mpi, ierr)
    spin2tot = spin2tot*voltot - spintot**2
#else

    denstot = densproc ! total charge
    spintot = spinproc ! S_tot
    spingrid = spingridproc ! spin charge on the mesh
    spin2tot = spin2proc*voltot ! S^2
    spin2tot = spin2tot - spintot**2 ! spin fluctuations
    minspin = minspinproc ! minimum value of the spin density
    minden = mindenproc ! minimum value of the spin density

#endif
    if (double_mesh) then
        densproc_cut = densproc_cut/denstot
        !    if(scale_hartree.gt.0.d0) then
        !     scale_hartreen=(1.d0/densproc_cut)**scale_hartree
        !    else
        !     scale_hartreen=1.d0
        !    endif
        if (rank .eq. 0) write (6, *) ' Old/new Total charge ratio ', densproc_cut
    end if

    if (rank .eq. 0) then
        !      if(rankrep.eq.0) then
        if (yeslsda) then
            write (6, '(a,2F18.6)') ' Total charge/spin: ', denstot, spintot
        else
            write (6, '(a,F18.6)') ' Total charge: ', denstot
            !          write(6,'(a,F10.6,I4)') ' Total charge: ',denstot,rankcolrep
            if (double_mesh) write (6, *) ' V Hartree q=0: ', vhq0
        end if
        if (minden .lt. 0.d0) write (6, *) ' Warning negative charge detected: ', minden/dens0
        if (abs(denstot - nel) .gt. 1d-4) then
            write (6, *) ' Warning total charge not conserved: ', denstot - nel
            !          write(6,'(a,F10.6,I4)') ' Warning total charge not conserved: ',denstot-nel,rankcolrep
        end if
        if (yeslsda) then
            if (minspin .lt. 0.d0) write (6, *) ' Warning negative spin charge detected: ', minspin/dens0
            if (abs(spintot - nelup + neldo) .gt. 1d-4) then
                write (6, '(a,F18.6)') ' Warning total spin not conserved: ', spintot - nelup + neldo
                write (6, '(a,F18.6)') ' Total spin occupations found: ', sum(occupations(:) - occupationdo(:))
            end if
        end if
    end if

    return

end subroutine cutdens

!
! This subroutine symmetrize the charge/spin densities according
! to the symmetry operation of the Bravais lattice.
! Adapted from QuantumESPRESSO
!
subroutine symmetrize_density(rho, buffer_grid, meshproc, nx, ny, nz, isymm, nsym)

    use symmetries, only: transform_point
    use allio, only: rankrep, commrep_mpi, nprocrep, rank, rankcolrep
    use setup, only: symtime
    use parallel_module, only: gather_from_grid, scatter_to_grid

    implicit none
#ifdef PARALLEL
    include "mpif.h"
#endif
    ! input
    integer, intent(in) :: nsym, isymm(3, 3, 48), meshproc, nx, ny, nz
    real(8), intent(inout) :: rho(meshproc), buffer_grid(meshproc, nprocrep)
    ! local
    real(8) :: dens
    logical, dimension(:, :, :), allocatable :: symflag
    real(8), dimension(:, :, :), allocatable :: dens_grid
    ! map for the transformed coordinates in the real space grid
    integer :: ri(48), rj(48), rk(48), isym, ierr
    integer :: i, j, k, i_, j_, k_

    real(8), external :: cclock
    real(8) :: symtimep

    if (rank .eq. 0) write (6, *) ' Warning: symmetrisation of the charge density! '

    allocate (dens_grid(nx, ny, nz), symflag(nx, ny, nz))
    dens_grid = 0.d0
    symflag = .true.

    if (nsym .eq. 1) return

    symtimep = cclock()

    call gather_from_grid(rho, buffer_grid, dens_grid, meshproc, nx, ny, nz)

    if (rankrep .eq. 0) then

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    if (symflag(i, j, k)) then
                        !
                        ! compute all the points in the grid which are
                        ! equivalent given the symmetries of the system
                        !
                        dens = 0.d0
                        do isym = 1, nsym
                            call transform_point(isymm(:, :, isym), i, j, k, &
                                                 nx, ny, nz, ri(isym), rj(isym), rk(isym))
                            dens = dens + dens_grid(ri(isym), rj(isym), rk(isym))
                        end do
                        dens = dens/nsym
                        !
                        ! sum contains the symmetrized charge density at point r.
                        ! now fill the star of r with this sum.
                        !
                        do isym = 1, nsym
                            dens_grid(ri(isym), rj(isym), rk(isym)) = dens
                            symflag(ri(isym), rj(isym), rk(isym)) = .false.
                        end do

                    end if
                end do
            end do
        end do

    end if

    call scatter_to_grid(dens_grid, buffer_grid, rho, meshproc, nx, ny, nz)

    deallocate (dens_grid)

    symtime = symtime + (cclock() - symtimep)

    return

end subroutine symmetrize_density

subroutine fill_phase_wf(molecorb, nbas, wf, wf_dim, buffer, buf_dim, same_phase, bufbuf)

    use constants, only: zzero, zone
    use setup, only: bands

    implicit none
    !
    ! input
    integer, intent(in) :: nbas, wf_dim, buf_dim, bufbuf
    complex(8), intent(inout) :: wf(wf_dim, bufbuf)
    complex(8), intent(in) :: molecorb(nbas, *)
    complex(8), intent(inout) :: buffer(buf_dim, *)
    logical same_phase
#ifdef _OFFLOAD
!$omp target data map(to:molecorb(:,1:bands),wf) map(from:buffer(:,1:bufbuf))
#endif
    if (.not. same_phase) call conjmat_(wf_dim, bufbuf, wf, wf_dim)
    call zgemm_('T', 'N', bands, bufbuf, nbas, zone, molecorb, nbas&
            &, wf, wf_dim, zzero, buffer, buf_dim)
    if (.not. same_phase) call conjmat_(wf_dim, bufbuf, wf, wf_dim)
#ifdef _OFFLOAD
!$omp end target data
#endif

    return

end subroutine fill_phase_wf

! end module density
