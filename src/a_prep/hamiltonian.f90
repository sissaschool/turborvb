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

subroutine uphamilt_new

    use cell, only: cellscale, car2cry, map, metric, unit_volume
    use constants, only: zzero, zone, Pi, ipc, safemin
    use allio
    use buffers
    use setup, only: overham, overhamdo, hamilt, hamiltdo, wf, &
            ax, ay, az, nx, ny, nz, nelorb, contracted_on, bufbuf, memlarge, yeslsda, &
            typedft, dent, spint, voltot, gradt, double_overs, &
            weightvh, weightcorr, weightxc, totvpot, exchange, ecorr, ehartree, volmesh, &
            eh_ew, zgemm_time, mindist, h_field, h_charge, nelorb3, &
            gridnospin, dgemm_timep, zgemm_timep, dgemm_time, rion_ref, &
            gridspin, nelorbu, vpotaa, gridcharge, gridnocharge, &
            & double_mesh, scale_z, l0_at, meshproc, nx0, ny0, nz0&
            &, weightx, weighty, weightz, minz_at, from_ions, rion_from&
            &, nx_at, ny_at, nz_at, rion_upload, scale_hartree, scale_hartreen, corr_hartree&
            &, mixing, vh_test

    use parallel_module, only: old_threads
    use fourier_module, only: update_vhartree, vhartree

#if defined __SCALAPACK
    use descriptors
    use setup, only: overhaml, overhamldo, hamiltl, hamiltldo, np_ortho, leg_ortho, desch
#endif

    implicit none

    ! local variables
    integer :: i, j, k, ii, jj, kk, ip, ind, max_index, indmesh, indtotal, nbufrep, indproc
    integer :: nbas_tot, nbas_1, nbas_2, nbas_3, buf_dim, buf_num, nxr, nyr, nzr, nxi&
            &, nxf, nyi, nyf, nzi, nzf, scalea, nx0n, ny0n, nz0n, buf_dim_conj

    real(8) :: scal, axn, ayn, azn
    real(8) :: totvpot_local, exchange_local, ehartree_local, ecorr_local
    real(8), external :: dlamch, cclock
#ifdef __SCALAPACK
    integer, allocatable :: irc_ip(:), nrc_ip(:), rank_ip(:, :)
    integer :: myid
#endif
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

    safemin = dlamch('S')
#ifdef __SCALAPACK
    allocate (irc_ip(np_ortho(1)))
    allocate (nrc_ip(np_ortho(1)))
    allocate (rank_ip(np_ortho(1), np_ortho(2)))
    irc_ip = 0
    nrc_ip = 0
    rank_ip = 0

    do j = 0, descla(la_npc_) - 1
        call descla_local_dims(irc_ip(j + 1), nrc_ip(j + 1), descla(la_n_), &
                               descla(la_nx_), np_ortho(1), j)
        do i = 0, descla(la_npr_) - 1
            call GRID2D_RANK('R', descla(la_npr_), descla(la_npc_), i, j, myid)
            rank_ip(i + 1, j + 1) = myid*leg_ortho
        end do
    end do
    hamiltl = 0.d0
    if (yeslsda .or. ipc .eq. 2) hamiltldo = 0.d0
#else
    hamilt = 0.d0
    if (yeslsda .or. ipc .eq. 2) hamiltdo = 0.d0
#endif

    ! mesh indices
    indmesh = 0

    ! indices of the buffers
    nbas_1 = nelorbu ! # of basis elements (contracted or uncontracted)
    nbas_2 = 2*nelorbu ! spin up part of the density-dependent hamiltonian
    nbas_3 = 3*nelorbu ! spin down part of the density-dependent hamiltonian
    nbas_tot = nelorb ! # of uncontracted basis elements >= nbas_1
    buf_dim = nelorb3 ! leading dimension of the buffer
    buf_dim_conj = nbas_1

    ! allocate buffers
    call allocate_buffers(bufbuf, nbas_tot, buf_dim, nbas_tot, buf_dim_conj, nbas_tot, thread_active)

    nbufrep = nprocrep*bufbuf

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

    ! update hartree potential in real space
    call update_vhartree

    ecorr = 0.d0
    exchange = 0.d0
    ehartree = 0.d0
    totvpot = 0.d0
    allocate (weightx(nx), weighty(ny), weightz(nz))
    weightx = 1.d0
    weighty = 1.d0
    weightz = 1.d0
    call upload_hamilt(1, nx, 1, ny, 1, nz, ax, ay, az, .true., memlarge, .true.)
    deallocate (weightx, weighty, weightz)

    if (double_mesh) then

        if (from_ions) then
            do ii = 1, nion
                !      locate atom
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

                    call upload_hamilt(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, .false., memlarge, .false.)

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

                    call set_interval(scalea, nx0n, ny0n, nz0n, rion(1, ii), rion_upload, axn, ayn, azn &
                            , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .false.&
                            &, weightx, weighty, weightz)

                    call upload_hamilt(nxi, nxf, nyi, nyf, nzi, nzf, axn, ayn, azn, .true., memlarge, .false.)

                    deallocate (weightx, weighty, weightz)
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

            call upload_hamilt(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, .false., memlarge, .false.)

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

            call set_interval(scalea, nx0n, ny0n, nz0n, rion_from, rion_upload, axn, ayn, azn &
                    , nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, .false.&
                    &, weightx, weighty, weightz)

            call upload_hamilt(nxi, nxf, nyi, nyf, nzi, nzf, axn, ayn, azn, .true., memlarge, .false.)

            deallocate (weightx, weighty, weightz)
        end if
        !    Restore original value of volmesh
        volmesh = ax*ay*az*unit_volume
    end if
    !
    ! collect the potentials from the grid
    !
#ifdef PARALLEL
    call reduce_base_real(1, totvpot, commrep_mpi, -1)
    call reduce_base_real(1, ecorr, commrep_mpi, -1)
    call reduce_base_real(1, exchange, commrep_mpi, -1)
    call reduce_base_real(1, ehartree, commrep_mpi, -1)
#endif
    totvpot = totvpot + vpotaa/2.d0

#ifdef __SCALAPACK
    ! symmetrize hamiltonian to reduce round off
    if (ipc .eq. 1) then
        call symmetrize_mat(nbas_1, size(hamiltl, 1), hamiltl)
        if (yeslsda) call symmetrize_mat(nbas_1, size(hamiltldo, 1), hamiltldo)
    else
        if (.not. double_overs .and. .not. yeslsda) then
            hamiltldo = hamiltl
            if (.not. same_phase) &
                call conjmat(size(hamiltldo, 1)/2, size(hamiltldo, 2), hamiltldo, size(hamiltldo, 1)/2)
        end if
        ! now symmetrize both matrices
        call symmetrize_mat_complex(nbas_1, size(hamiltl, 1)/2, hamiltl)
        call symmetrize_mat_complex(nbas_1, size(hamiltldo, 1)/2, hamiltldo)
    end if

#else

    ! collect hamilt/hamiltdo distributed on the grid and symmetrize
    if (ipc .eq. 1) then
#ifdef PARALLEL
        call reduce_base_real(size(hamilt), hamilt, commrep_mpi, -1)
        if (yeslsda) call reduce_base_real(size(hamiltdo), hamiltdo, commrep_mpi, -1)
#endif
    else
#ifdef PARALLEL
        call reduce_base_real(size(hamilt), hamilt, commrep_mpi, -1)
#endif
        if (double_overs .or. yeslsda) then
#ifdef PARALLEL
            call reduce_base_real(size(hamiltdo), hamiltdo, commrep_mpi, -1)
#endif
        else
            hamiltdo = hamilt
            if (.not. same_phase) &
                call conjmat(size(hamiltdo, 1)/2, size(hamiltdo, 2), hamiltdo, size(hamiltdo, 1)/2)
        end if
    end if

    ! symmetrize to reduce round-off
    if (ipc .eq. 1) then
        call symmetrize_mat(nbas_1, size(hamilt, 1), hamilt)
        if (yeslsda) call symmetrize_mat(nbas_1, size(hamiltdo, 1), hamiltdo)
    else
        call symmetrize_mat_complex(nbas_1, size(hamilt, 1)/2, hamilt)
        call symmetrize_mat_complex(nbas_1, size(hamiltdo, 1)/2, hamiltdo)
    end if

#endif
    !
    ! adding kinetic energy part, it doesn't change throughout
    ! the calculation since it's indipendent from the density
    !
#ifdef __SCALAPACK

    if (descla(lambda_node_) > 0) then
        if (ipc .eq. 1) then
            hamiltl = hamiltl + overhaml
            if (yeslsda) hamiltldo = hamiltldo + overhaml
        else
            hamiltl = hamiltl + overhaml
            hamiltldo = hamiltldo + overhamldo
        end if
    end if

#else
    !
    ! add the density-independent part of the Hamiltonian
    !
    if (ipc .eq. 1) then
        hamilt = hamilt + overham
        if (yeslsda) hamiltdo = hamiltdo + overham
    else
        hamilt = hamilt + overham
        hamiltdo = hamiltdo + overhamdo
    end if

#endif

#ifdef DEBUG
    if (rank .eq. 0) then
        write (6, *)
        write (6, *) ' Checking energy contributions:'
        write (6, *) ' total pot/exchange/correlation/hartree'
        write (6, *) totvpot, exchange, ecorr, ehartree
        write (6, *)
    end if
#endif

    if (corr_hartree .and. scale_hartree .gt. 0 .and. vh_test .gt. 0.d0) then
        scale_hartreen = scale_hartreen + 0.5d0*mixing*(ehartree/vh_test - 1.d0)
        if (rank .eq. 0) write (6, *) ' New trial value of scale_hartree =', scale_hartreen
    end if

    ! deallocate buffers and scratch vectors
    !
    call deallocate_buffers()
    !
    deallocate (distp_scratch)
    deallocate (r_scratch, rmu_scratch, rmusin_scratch, rmucos_scratch)
    deallocate (x_buffer)
    !
#ifdef __SCALAPACK
    deallocate (irc_ip, nrc_ip, rank_ip)
#endif

#if defined (_OPENMP) && defined (__NOOMPDFT)
    call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

#ifdef PARALLEL
    call mpi_barrier(MPI_COMM_WORLD, ierr) ! make sure all the nodes have released the memory
!$omp barrier
#endif

#ifdef PARALLEL
    call mpi_allreduce(ind, max_index, 1, MPI_INTEGER, MPI_MAX, commrep_mpi, ierr)
#else
    max_index = ind
#endif
    if (max_index .ne. 0) &
        call error(' uphamilt ', ' check input nbufd and/or code ', 1, rank)

    return

contains

    subroutine upload_hamilt(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, add, memlarge, do_field)
        implicit none
        integer nxi, nxf, nyi, nyf, nzi, nzf, ii, bufmax
        real*8 ax, ay, az, volmeshu, r_eion
        logical add, memlarge, do_field
        volmesh = ax*ay*az*unit_volume

        ind = 0
        indproc = 0
        indtotal = 0
        do k = nzi, nzf
            do j = nyi, nyf
                do i = nxi, nxf
                    indtotal = indtotal + 1

                    if (indproc .eq. rankrep) then
                        ! grid point buffer counter
                        indmesh = indmesh + 1
                        ind = ind + 1
                        weight_buff(ind) = weightx(i - nxi + 1)*weighty(j - nyi + 1)*weightz(k - nzi + 1)
                        if (.not. add) weight_buff(ind) = -weight_buff(ind)
                        if (.not. memlarge) then ! in this case no optimization of the overlaps is done
                            ! buffering the grid point coordinates
                            x_buffer(:, ind) = i*ax*at(:, 1) + j*ay*at(:, 2) + k*az*at(:, 3) + rion_upload(:)

                            !             x_buffer(1,ind)=i*ax+rion_upload(1)
                            !             x_buffer(2,ind)=j*ay+rion_upload(2)
                            !             x_buffer(3,ind)=k*az+rion_upload(3)
                        end if ! endif memlarge
                    end if
                    indproc = indproc + 1
                    if (indproc .eq. nprocrep) indproc = 0

                    if (indtotal .eq. nbufrep .or. (i .eq. nxf .and. j .eq. nyf .and. k .eq. nzf)) then

!$omp parallel default(shared) private(tid,ii)
#if defined(_OPENMP)
                        tid = 1 + omp_get_thread_num()
#else
                        tid = 1
#endif
!$omp do
                        do ii = 1, ind
                            if (memlarge) then ! in this case no optimization of the overlaps is done
                                buffer(1:ipc*nbas_1, ii) = wf(1:ipc*nbas_1, indmesh - ind + ii) ! up spin
                                if (ipc .eq. 2) then
                                    if (double_overs) then
                                        buffer(1:ipc*nbas_1, bufbuf + ii) &
                                            = wf(ipc*nbas_1 + 1:ipc*nbas_2, indmesh - ind + ii) ! down spin
                                    elseif (yeslsda) then
                                        if (same_phase) then
                                            buffer(1:ipc*nbas_1, bufbuf + ii) &
                                                = wf(1:ipc*nbas_1, indmesh - ind + ii) ! down spin
                                        else
                                            ! The complex conjugate
                                            buffer(1:ipc*nbas_1:2, bufbuf + ii) &
                                                = wf(1:ipc*nbas_1:2, indmesh - ind + ii) ! down spin
                                            buffer(2:ipc*nbas_1:2, bufbuf + ii) &
                                                = -wf(2:ipc*nbas_1:2, indmesh - ind + ii) ! down spin
                                        end if
                                    end if
                                end if

                            else
                                call compute_wf_one_grid_point(x_buffer(1, ii), tid, ii, indmesh - ind + ii)
                            end if ! endif memlarge
                        end do
!$omp end do nowait
!$omp end parallel

                        if (contracted_on .and. .not. memlarge) then
                            if (ipc .eq. 1) then
#ifdef _OFFLOAD
!$omp target data map(from:buffer(:,1:ind)) map(to:buffer_on(:,1:ind),mu_c(:,1:nbas_1))
#endif
                                call dgemm_('T', 'N', nbas_1, ind, nbas_tot, 1.d0, mu_c, &
                                            nbas_tot, buffer_on, nbas_tot, 0.d0, buffer, buf_dim)
#ifdef _OFFLOAD
!$omp end target data
#endif
                            else
#ifdef _OFFLOAD
                                if (double_overs .or. yeslsda) then
                                    bufmax = bufbuf + ind
                                else
                                    bufmax = ind
                                end if
!$omp target data map(from:buffer(:,1:bufmax)) map(to:buffer_on(:,1:bufmax),mu_c(:,1:nbas_1))
#endif
                                call zgemm_('T', 'N', nbas_1, ind, nbas_tot, zone, mu_c, &
                                            nbas_tot, buffer_on, nbas_tot, zzero, buffer, buf_dim)
                                if (double_overs .or. yeslsda) then
                                    call zgemm_('T', 'N', nbas_1, ind, nbas_tot, zone, mu_c, &
                                                nbas_tot, buffer_on(1, bufbuf + 1), nbas_tot, zzero &
                                                , buffer(1, bufbuf + 1), buf_dim)
                                end if
#ifdef _OFFLOAD
!$omp end target data
#endif
                            end if
                        end if
                        !
                        ! exchange correlation does not depend on the k-point
                        ! computed locally on the pool
                        !
!$omp parallel do default(shared) reduction(+:totvpot,exchange,ehartree,ecorr) &
!$omp private(ii,totvpot_local,exchange_local,ehartree_local,ecorr_local,volmeshu)
                        do ii = 1, ind
                            volmeshu = volmesh*weight_buff(ii)
                            call compute_xc_one_grid_point(ii, indmesh - ind + ii &
                                                           , totvpot_local, exchange_local, ehartree_local, ecorr_local &
                                                           , volmeshu, do_field)
                            totvpot = totvpot + totvpot_local
                            exchange = exchange + exchange_local
                            ehartree = ehartree + ehartree_local
                            ecorr = ecorr + ecorr_local
                        end do ! enddo bufbuf
!$omp end parallel do

                        ! avoid any large temporary allocation in zgemm
                        if (ipc .eq. 2) then
                            do ii = 1, ind
                                buf_conj(1:2*buf_dim_conj, ii) = buffer(1:2*buf_dim_conj, ii)*weight_buff(ii)
                            end do
                            call conjmat(buf_dim_conj, ind, buf_conj, buf_dim_conj)
                            if (double_overs .or. yeslsda) then
                                do ii = 1, ind
                                    buf_conj(1:2*buf_dim_conj, bufbuf + ii) = buffer(1:2*buf_dim_conj, bufbuf + ii)*weight_buff(ii)
                                end do
                                call conjmat(buf_dim_conj, ind, buf_conj(1, bufbuf + 1), buf_dim_conj)
                            end if
                        else
                            do ii = 1, ind
                                buf_conj(1:buf_dim_conj, ii) = buffer(1:buf_dim_conj, ii)*weight_buff(ii)
                            end do
                        end if

                        ! now compute < \Psi_l | V
#ifdef __SCALAPACK

#if defined (_OPENMP) && defined (__NOOMPDFT)
                        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

                        if (ipc .eq. 1) then
                            call fill_matrix_scalapack(hamiltl, buffer, buf_dim, buf_conj, buf_dim_conj, ind &
                                                       , irc_ip, nrc_ip, rank_ip, descla(nlax_), nbas_1)
                            if (yeslsda) &
                                call fill_matrix_scalapack(hamiltldo, buffer, buf_dim, buf_conj, buf_dim_conj, ind &
                                                           , irc_ip, nrc_ip, rank_ip, descla(nlax_), nbas_2)
                        else
                            call fill_matrix_scalapack_complex(hamiltl, buffer, buf_dim, buf_conj, buf_dim_conj, ind &
                                                               , irc_ip, nrc_ip, rank_ip, descla(nlax_), nbas_1)

                            if (double_overs .or. yeslsda) then
                                call fill_matrix_scalapack_complex(hamiltldo, buffer, buf_dim &
                                                                   , buf_conj(1, bufbuf + 1), buf_dim_conj, ind &
                                                                   , irc_ip, nrc_ip, rank_ip, descla(nlax_), nbas_2)
                            end if
                        end if
#else
#ifdef _OFFLOAD
                        if ((double_overs .or. yeslsda) .and. ipc .eq. 2) then
                            bufmax = bufbuf + ind
                        else
                            bufmax = ind
                        end if
!$omp target data map(hamilt,hamiltdo) &
!$omp map(to:buf_conj(:,1:bufmax),buffer(:,1:ind))
#endif

                        if (ipc .eq. 1) then
                            dgemm_timep = cclock()
                            call dgemm_('N', 'T', nbas_1, nbas_1, ind, volmesh, buf_conj, &
                                        buf_dim_conj, buffer(nbas_1 + 1, 1), buf_dim, 1.d0, hamilt, nbas_1)
                            if (yeslsda) &
                                call dgemm_('N', 'T', nbas_1, nbas_1, ind, volmesh, buf_conj, &
                                            buf_dim_conj, buffer(nbas_2 + 1, 1), buf_dim, 1.d0, hamiltdo, nbas_1)
                            dgemm_time = dgemm_time + cclock() - dgemm_timep
                        else
                            zgemm_timep = cclock()
                            ! up spin
                            call zgemm_('N', 'T', nbas_1, nbas_1, ind, dcmplx(volmesh), buf_conj, &
                                        buf_dim_conj, buffer(2*nbas_1 + 1, 1), buf_dim, zone, hamilt, nbas_1)
                            ! down spin in any case for a complex wave function
                            if (double_overs .or. yeslsda) then
                                call zgemm_('N', 'T', nbas_1, nbas_1, ind, dcmplx(volmesh), buf_conj(1, bufbuf + 1), &
                                            buf_dim_conj, buffer(2*nbas_2 + 1, 1), buf_dim, zone, hamiltdo, nbas_1)
                            end if
                            zgemm_time = zgemm_time + cclock() - zgemm_timep
                        end if
#ifdef _OFFLOAD
!$omp end target data
#endif

#endif

#if defined (_OPENMP) && defined (__NOOMPDFT)
                        call omp_set_num_threads(1) ! restore the scalar code
#endif

                        ind = 0
                        indtotal = 0

                    end if ! if bufbuf
                end do ! nx
            end do ! ny
        end do ! nz
    end subroutine upload_hamilt

    subroutine compute_wf_one_grid_point(x, tid, ind, indmesh)
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

        ! computing the wave function
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, x, 1, r_scratch(1, tid), rmu_scratch(1, tid), &
                     dupr, zetar, rion, distp_scratch(1, tid), wf_threading_scratch(1, tid), nbas_tot, nion, kion, &
                     iflagnorm, cnorm, LBox, rmucos_scratch(1, tid), rmusin_scratch(1, tid), mindist, &
                     indpar_tab, indorb_tab, indshell_tab, .true.)
        if (double_overs) then
            call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, x, 1, r_scratch(1, tid), rmu_scratch(1, tid), &
                         dupr, zetar, rion, distp_scratch(1, tid), wf_threading_scratch_down(1, tid), nbas_tot, nion, kion, &
                         iflagnorm, cnorm, LBox, rmucos_scratch(1, tid), rmusin_scratch(1, tid), mindist, &
                         indpar_tab, indorb_tab, indshell_tab, .false.)
        end if
        if (n_body_on .ne. 0) then
            jas_1body = -scale_one_body
            do jj = 1, nion
                if (iespbc) then
                    rc(1) = x(1) - rion(1, jj)
                    rc(2) = x(2) - rion(2, jj)
                    rc(3) = x(3) - rion(3, jj)
!                   call  CartesianToCrystal(rc, 1)
                    rc(:) = car2cry(:, 1)*rc(1) + car2cry(:, 2)*rc(2) + car2cry(:, 3)*rc(3)
                    do kk = 1, 3
                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                    end do
                    r0 = norm_metric(rc, metric)
                else
                    rc(1) = (x(1) - rion(1, jj))*costz(jj)
                    rc(2) = (x(2) - rion(2, jj))*costz(jj)
                    rc(3) = (x(3) - rion(3, jj))*costz(jj)
                    r0 = dsqrt(sum(rc(:)**2))
                end if
                jas_1body = jas_1body - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
            end do
            jas_1body = dexp(jas_1body)
            do jj = 1, ipc*nbas_tot
                wf_threading_scratch(jj, tid) = wf_threading_scratch(jj, tid)*jas_1body
                if (double_overs) wf_threading_scratch_down(jj, tid) = wf_threading_scratch_down(jj, tid)*jas_1body
            end do
        end if

        if (contracted_on) then

            buffer_on(1:ipc*nbas_tot, ind) = wf_threading_scratch(1:ipc*nbas_tot, tid)
            if (ipc .eq. 2) then
                if (double_overs) then
                    buffer_on(1:ipc*nbas_tot, bufbuf + ind) = wf_threading_scratch_down(1:ipc*nbas_tot, tid)
                elseif (yeslsda) then
                    if (same_phase) then
                        buffer_on(1:ipc*nbas_tot, bufbuf + ind) = wf_threading_scratch(1:ipc*nbas_tot, tid)
                    else
                        buffer_on(1:ipc*nbas_tot:2, bufbuf + ind) = wf_threading_scratch(1:ipc*nbas_tot:2, tid)
                        buffer_on(2:ipc*nbas_tot:2, bufbuf + ind) = -wf_threading_scratch(2:ipc*nbas_tot:2, tid)
                    end if
                end if
            end if

        else
            buffer(1:ipc*nbas_tot, ind) = wf_threading_scratch(1:ipc*nbas_tot, tid)
            if (ipc .eq. 2) then
                if (double_overs) then
                    buffer(1:ipc*nbas_tot, bufbuf + ind) = wf_threading_scratch_down(1:ipc*nbas_tot, tid)
                elseif (yeslsda) then
                    if (same_phase) then
                        buffer(1:ipc*nbas_tot, bufbuf + ind) = wf_threading_scratch(1:ipc*nbas_tot, tid)
                    else
                        buffer(1:ipc*nbas_tot:2, bufbuf + ind) = wf_threading_scratch(1:ipc*nbas_tot:2, tid)
                        buffer(2:ipc*nbas_tot:2, bufbuf + ind) = -wf_threading_scratch(2:ipc*nbas_tot:2, tid)
                    end if
                end if
            end if
        end if
    end subroutine compute_wf_one_grid_point

    subroutine compute_xc_one_grid_point(ind, indmesh, totvpot, exchange, ehartree, ecorr, volmesh, do_field)
        ! it is thread-safe and intended to be called from threaded region.
        ! local density on the real space mesh
        implicit none
        ! arguments
        ! index in the grid buffer
        integer, intent(in) :: ind
        ! global index of the current grid point
        integer, intent(in) :: indmesh
        ! exchange and corr energy
        real(8), intent(out) :: totvpot, exchange, ehartree, ecorr

        ! local variables
        real(8) :: dens_true, rsdft, densup, densdo
        real(8) :: vx, vc, ec, ex, exup, exdo, vxup, vxdo, rsup, rsdo, vcup, vcdo, zetaxi
        real(8) :: vpotden, vh, v_field, v_charge, vpotup, vpotdo, volmesh
        logical do_field

        dens_true = max(dent(indmesh), safemin)

        if (typedft .ne. 0) then
            ! find rs
            rsdft = (3.d0/(4.d0*Pi*dens_true))**(1.d0/3.d0)
            ! compute DFT exchange and correlation energies and potentials
            if (abs(typedft) .le. 2) then ! LDA
                call slater(rsdft, ex, vx)
                call pz(rsdft, typedft, ec, vc)

            elseif (abs(typedft) .eq. 3) then ! LDA KZK
                call slaterKZK(rsdft, ex, vx, voltot, typedft)
                call pzKZK(rsdft, ec, vc, voltot)

            elseif (abs(typedft) .eq. 4) then ! SLDA
                densup = dens_true + spint(indmesh)
                densdo = dens_true - spint(indmesh)

                if (densup .gt. 0.d0) then
                    rsup = (3.d0/(4.d0*Pi*densup))**(1.d0/3.d0)
                    call slater(rsup, exup, vxup)
                    exup = exup/2.d0
                else
                    exup = 0.d0
                    vxup = 0.d0
                end if

                if (densdo .gt. 0.d0) then
                    rsdo = (3.d0/(4.d0*Pi*densdo))**(1.d0/3.d0)
                    call slater(rsdo, exdo, vxdo)
                    exdo = exdo/2.d0
                else
                    exdo = 0.d0
                    vxdo = 0.d0
                end if

                zetaxi = spint(indmesh)/dens_true
                call pz_spin(rsdft, zetaxi, ec, vcup, vcdo)

            elseif (abs(typedft) .eq. 5) then ! SLDA KZK

                densup = max(dens_true + spint(indmesh), safemin)
                densdo = max(dens_true - spint(indmesh), safemin)
                rsup = (3.d0/(4.d0*Pi*densup))**(1.d0/3.d0)
                call slaterKZK(rsup, exup, vxup, voltot, typedft)
                exup = exup/2.d0
                rsdo = (3.d0/(4.d0*Pi*densdo))**(1.d0/3.d0)
                call slaterKZK(rsdo, exdo, vxdo, voltot, typedft)
                exdo = exdo/2.d0
                zetaxi = spint(indmesh)/dens_true
                call pz_spinKZK(rsdft, zetaxi, ec, vcup, vcdo, voltot)

            end if

        else ! typedft==0, no exchange-correalation (Hartree)

            vx = 0.d0
            vc = 0.d0
            ec = 0.d0
            ex = 0.d0
            exup = 0.d0
            exdo = 0.d0
            vxup = 0.d0
            vxdo = 0.d0
            vcup = 0.d0
            vcdo = 0.d0

        end if
        !
        ! compute < \Psi_l | V | \Psi_m > = hamilt(l,m) for the buffer
        ! V is a LOCAL potential, terms are additive
        !
        ! first perform : V | \Psi_m >
        !
        !      if(indmesh.le.meshproc) then
        vh = vhartree(indmesh)
        v_charge = 0.d0
        if (h_charge .ne. 0.d0 .and. do_field) then
            if (.not. gridnocharge(indmesh)) then
                if (gridcharge(indmesh)) then
                    v_charge = -h_charge
                else
                    v_charge = h_charge
                end if
            end if
        end if
        !      else
        !       vh=0.d0
        !       v_charge=0.d0
        !      endif

        if (yeslsda) then

            if (indmesh .le. meshproc .and. do_field .and. h_field .ne. 0.d0) then
                if (h_field .eq. 0.d0 .or. gridnospin(indmesh)) then
                    v_field = 0.d0
                else
                    if (gridspin(indmesh)) then
                        v_field = -h_field
                    else
                        v_field = h_field
                    end if
                end if
            else
                v_field = 0.d0
            end if

            !
            ! Addressing of the buffer:
            !   buffer(1:ipc*nbas_1,1:bufbuf)          = w.f. spin up electrons
            !   buffer(1:ipc*nbas_1,bufbuf+1:2*bufbuf) = w.f. spin down electrons
            !   buffer(nbas_1+1:nbas_2) = LDA potential  (spin up electrons)
            !   buffer(nbas_2+1:nbas_3) = LSDA potential (spin down electrons)
            !
            ! up spin
            vpotup = vh*weightvh + vxup*weightxc + vcup*weightcorr + v_field + v_charge
            buffer(ipc*nbas_1 + 1:ipc*nbas_2, ind) = vpotup*buffer(1:ipc*nbas_1, ind)
            ! down spin
            vpotdo = vh*weightvh + vxdo*weightxc + vcdo*weightcorr - v_field + v_charge
            if (ipc .eq. 2) then
                buffer(2*nbas_2 + 1:2*nbas_3, ind) = vpotdo*buffer(1:2*nbas_1, bufbuf + ind)
            else
                buffer(ipc*nbas_2 + 1:ipc*nbas_3, ind) = vpotdo*buffer(1:ipc*nbas_1, ind)
            end if
            !
            ! total potential up+down
            !
            totvpot = volmesh* &
                    &((dens_true*ec - 0.5d0*densup*vcup - 0.5d0*densdo*vcdo)*weightcorr&
                            & + weightxc*(exup*densup + exdo*densdo - 0.5d0*densup*vxup - 0.5d0*densdo*vxdo)&
                            & - 0.5d0*weightvh*vh*dens_true)
            ! exchange part
            exchange = volmesh*(exup*densup + exdo*densdo)

        else

            vpotden = vh*weightvh + vx*weightxc + vc*weightcorr + v_charge
            buffer(ipc*nbas_1 + 1:ipc*nbas_2, ind) = vpotden*buffer(1:ipc*nbas_1, ind)
            if (double_overs .and. ipc .eq. 2) then
                buffer(2*nbas_2 + 1:2*nbas_3, ind) = vpotden*buffer(1:2*nbas_1, bufbuf + ind)
            end if

            totvpot = volmesh*dens_true*((ec - vc)*weightcorr&
                    & + weightxc*(ex - vx) - 0.5d0*weightvh*vh)
            exchange = volmesh*ex*dens_true

        end if

        ehartree = 0.5d0*volmesh*vh*dens_true
        ecorr = volmesh*ec*dens_true

    end subroutine compute_xc_one_grid_point
end subroutine uphamilt_new

!
! INPUT
! molecorb_old -> molecular orthonormal orbitals
! dent (distributed over the grid) density corresponding to the molecular orbitals
! hamiltl (hamilt) matrix elements Kohn-Sham hamiltonian (not distributed)
! All these quantities are not changed in output.
!
! OUTPUT
! evardft -> variational DFT energy corresponding to the
!            density and the corresponding molecular orbitals.
! If dent != sum(occupations(:) molecorbold^T  w.f. ), then the output
! edftvar is not variational
! Used molecorb, changed on output.
!

subroutine evalevar

    use constants, only: ipc, zone, zzero
    use allio, only: commrep_mpi, commcolrep_mpi
    use kpoints_mod, only: sum_kpoints_scalar_real8
    use setup, only: molecorb, molecorbdo, molecorb_old, molecorbdo_old, hamilt, hamiltdo
    use setup, only: bands, edftvar, nelorb, nelorbu, occupations, nelocc, &
                     occupationdo, yeslsda, indk, neloccdo, totvpot

    use parallel_module, only: old_threads
#ifdef __SCALAPACK
    use descriptors
    use setup, only: molecorbl, molecorbldo, hamiltl, hamiltldo, desch
#endif
    implicit none
    ! local
    real(8), allocatable :: eigv_scratch(:, :)
    real(8), external :: ddot
    complex(8), external :: zdotc_
    integer ii
#ifdef PARALLEL
    include 'mpif.h'
#endif
#ifdef __SCALAPACK
    integer :: i, j
    integer :: irow, icol
#endif
#if defined (_OPENMP) && defined (__NOOMPDFT)
    integer, external :: omp_get_max_threads
    call omp_set_num_threads(1) ! scalar code
#endif
    !
    !  Calculation variational DFT energy with the molecorb_old molecular orbitals.
    !
#ifdef __SCALAPACK

    molecorb = 0.d0
    if (descla(lambda_node_) > 0) then
        allocate (eigv_scratch(ipc*descla(nlax_), descla(nlax_)))
        eigv_scratch = 0.d0
        ! load molecorb in molecorbl
        irow = descla(ilar_)
        icol = descla(ilac_)
        do j = 1, descla(nlac_)
            do i = 1, descla(nlar_)
                if ((j + icol - 1) <= bands) then
                    if (ipc .eq. 1) then
                        molecorbl(i, j) = molecorb_old((i + irow - 1), (j + icol - 1))
                    else
                        molecorbl(2*i - 1:2*i, j) = &
                            molecorb_old(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1))
                    end if
                end if
            end do
        end do

#if defined (_OPENMP) && defined (__NOOMPDFT)
        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

        if (ipc .eq. 1) then
            call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.d0, hamiltl, 1, 1,&
                 & desch, molecorbl, 1, 1, desch, 0.d0, eigv_scratch, 1, 1, desch)
        else
            call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, zone, hamiltl, 1, 1,&
                 & desch, molecorbl, 1, 1, desch, zzero, eigv_scratch, 1, 1, desch)
        end if
#if defined (_OPENMP) && defined (__NOOMPDFT)
        call omp_set_num_threads(1) ! restore the scalar
#endif

        irow = descla(ilar_)
        icol = descla(ilac_)
        do j = 1, descla(nlac_)
            do i = 1, descla(nlar_)
                if ((j + icol - 1) <= bands) then
                    if (ipc .eq. 1) then
                        molecorb((i + irow - 1), (j + icol - 1)) = eigv_scratch(i, j)
                    else
                        molecorb(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1)) = eigv_scratch(2*i - 1:2*i, j)
                    end if
                end if
            end do
        end do
        deallocate (eigv_scratch)
    end if
#ifdef PARALLEL
    call reduce_base_real(size(molecorb), molecorb, commrep_mpi, -1)
#endif

#else

    if (ipc .eq. 1) then
        call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, hamilt, nelorbu&
                &, molecorb_old, nelorb, 0.d0, molecorb, nelorbu)
    else
        call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, hamilt, nelorbu&
                &, molecorb_old, nelorb, zzero, molecorb, nelorbu)
    end if
#endif
! compute DFT energy
    do ii = 1, nelocc
        molecorb(:, ii) = occupations(ii)*molecorb(:, ii)
    end do
    edftvar = 0.d0
    do ii = 1, nelocc
        if (ipc .eq. 1) then
            edftvar = edftvar + ddot(nelorbu, molecorb(1, ii), 1, molecorb_old(1, ii), 1)
        else
            edftvar = edftvar + zdotc_(nelorbu, molecorb(1, ii), 1, molecorb_old(1, ii), 1)
        end if
    end do

    if (yeslsda .or. ipc .eq. 2) then

#ifdef __SCALAPACK

        molecorbdo = 0.d0
        if (descla(lambda_node_) > 0) then
            allocate (eigv_scratch(ipc*descla(nlax_), descla(nlax_)))
            eigv_scratch = 0.d0
            !       load molecorb in molecorbl
            irow = descla(ilar_)
            icol = descla(ilac_)
            do j = 1, descla(nlac_)
                do i = 1, descla(nlar_)
                    if ((j + icol - 1) <= bands) then
                        if (ipc .eq. 1) then
                            molecorbldo(i, j) = molecorbdo_old((i + irow - 1), (j + icol - 1))
                        else
                            molecorbldo(2*i - 1:2*i, j) = &
                                molecorbdo_old(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1))
                        end if
                    end if
                end do
            end do

#if defined (_OPENMP) && defined (__NOOMPDFT)
            call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

            if (ipc .eq. 1) then
                call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.d0, hamiltldo, 1, 1,&
                     & desch, molecorbldo, 1, 1, desch, 0.d0, eigv_scratch, 1, 1, desch)
            else
                call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, zone, hamiltldo, 1, 1,&
                     & desch, molecorbldo, 1, 1, desch, zzero, eigv_scratch, 1, 1, desch)
            end if

#if defined (_OPENMP) && defined (__NOOMPDFT)
            call omp_set_num_threads(1) ! restore the scalar
#endif

            irow = descla(ilar_)
            icol = descla(ilac_)
            do j = 1, descla(nlac_)
                do i = 1, descla(nlar_)
                    if ((j + icol - 1) <= bands) then
                        if (ipc .eq. 1) then
                            molecorbdo((i + irow - 1), (j + icol - 1)) = eigv_scratch(i, j)
                        else
                            molecorbdo(2*(i + irow - 1) - 1:2*(i + irow - 1), (j + icol - 1)) = eigv_scratch(2*i - 1:2*i, j)
                        end if
                    end if
                end do
            end do
            deallocate (eigv_scratch)
        end if
#ifdef PARALLEL
        call reduce_base_real(size(molecorbdo), molecorbdo, commrep_mpi, -1)
#endif

#else
        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, hamiltdo, nelorbu&
                    &, molecorbdo_old, nelorb, 0.d0, molecorbdo, nelorbu)
        else
            call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, hamiltdo, nelorbu&
                    &, molecorbdo_old, nelorb, zzero, molecorbdo, nelorbu)
        end if
#endif

        do ii = 1, neloccdo
            molecorbdo(:, ii) = occupationdo(ii)*molecorbdo(:, ii)
        end do

        do ii = 1, neloccdo
            if (ipc .eq. 1) then
                edftvar = edftvar + ddot(nelorbu, molecorbdo(1, ii), 1, molecorbdo_old(1, ii), 1)
            else
                edftvar = edftvar + zdotc_(nelorbu, molecorbdo(1, ii), 1, molecorbdo_old(1, ii), 1)
            end if
        end do

    end if
! sum over k-points
    call sum_kpoints_scalar_real8(edftvar, commcolrep_mpi, -1)
    edftvar = edftvar + totvpot

#if defined (_OPENMP) && defined (__NOOMPDFT)
    call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

    return

end subroutine evalevar
