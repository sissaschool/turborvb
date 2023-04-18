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

! This subroutine initializes the main matrices used in the SC run. It is
! written to be generally used for both real and complex matrices.
! NB: ipc.eq.1 -> real case
!     ipc.eq.2 -> complex case, all matrices and buffer/scratch arrays
!                 are defined with double dimension.

subroutine initialize_mats_new

    ! use cell,            only: ApplyPBC,CartesianToCrystal,map,metric,unit_volume
    ! use ewald,           only: ewaldion1b,ewaldel1b
    ! use constants,       only: zzero,zone,Pi,ipc
    use allio
    use parallel_module, only: old_threads
    use buffers
    use setup, only: overs, oversdo, overham, overhamdo, wf, write_matrix
    use setup, only: ax, ay, az, nx, ny, nz, nelorb, contracted_on, bufbuf, memlarge, yeslsda, &
            weightvh, volmesh, zgemm_time, mindist, nelorb3, &
            dgemm_timep, zgemm_timep, dgemm_time, meshproc_wf, noread, rion_ref, &
            indk, vlocaltot, vq0tots, unit_scratch_distributed, &
            meshproc, nelorbu, double_overs, init_time, wf_dim, l0_at, double_mesh, &
            &scale_z, mesh, nx0, ny0, nz0, weightx, weighty, weightz, meshproc_tot, &
            &volmesh_proc, minz_at, from_ions, rion_from, nx_at, ny_at, nz_at, rion_upload
    use fourier_module, only: ncub_min, ncub_max, nx_proc, nxny_proc, nx8, nxny8&
            &, ind_init
#ifdef __SCALAPACK
    use descriptors
    use setup, only: np_ortho, leg_ortho, desch, oversl, oversldo, overhaml, overhamldo
#endif

    implicit none
    ! local variables
    integer :: nbas_1, nbas_tot, buf_dim, nbuf, i, ii, jj, j, k, ind, indmesh, indmesh_local, &
            & max_index, buf_num, iii, mpi_ierr, indproc, nxr, nyr, nzr, nxi, nxf, nyi, nyf, &
            & nzi, nzf, scalea, nx0n, ny0n, nz0n, buf_dim_conj
    real(8) :: timep_init, timep_buf, single_buffer_time, costq0, costq0_one, vlocaltot_one, scratch_mem, axn, ayn, azn
    integer, external :: iesdr1iesd, mod_true
    real(8), external :: cclock
    integer :: tid, thread_active
    logical :: to_print

    ! local buffers
    real(8), dimension(:, :), allocatable :: wpseudo_and_psip_scratch, tcost_scratch
    real(8), dimension(:, :), allocatable :: prefactor_scratch
    real(8), dimension(:, :, :), allocatable :: ivic_scratch
    real(8), dimension(:, :), allocatable :: distp_scratch, r_scratch, dist_scratch, dist_shift_scratch
    real(8), dimension(:, :), allocatable :: rmu_scratch, rmusin_scratch, rmucos_scratch
    real(8), dimension(:, :, :), allocatable :: x_buffer
    real(8), dimension(:, :), allocatable :: angle_buffer
    double complex, allocatable :: pshfac1(:), pshfac2(:), pshfac3(:)
    double complex twiddlep

#ifdef __SCALAPACK
    integer, allocatable :: irc_ip(:), nrc_ip(:), rank_ip(:, :)
    integer :: myid, nrr, nc, ir, ic, indpc, indpr, root
#endif

#ifdef PARALLEL
    include "mpif.h"
#endif
#if defined(_OPENMP)
    integer, external :: omp_get_max_threads, omp_get_thread_num
    thread_active = omp_get_max_threads()
#if defined(__NOOMPDFT)
    call omp_set_num_threads(1) ! scalar code
    thread_active = 1
#endif
#else
    thread_active = 1
#endif
    !
    ! initialize indices and processors of the SCALAPACK grid
    !
#ifdef __SCALAPACK
    ! irc_ip(*)    = the indices of the first local element
    ! nrc_ip(*)    = the total dimensions (# of elements) of the local block
    ! rank_ip(i,j) = rank of the MPI task associated with element (i,j) of the grid
    allocate (irc_ip(np_ortho(1)))
    allocate (nrc_ip(np_ortho(1)))
    allocate (rank_ip(np_ortho(1), np_ortho(2)))
    irc_ip = 0
    nrc_ip = 0
    rank_ip = 0

    do j = 0, descla(la_npc_) - 1
        call descla_local_dims(irc_ip(j + 1), nrc_ip(j + 1), descla(la_n_), descla(la_nx_), np_ortho(1), j)
        do i = 0, descla(la_npr_) - 1
            ! compute the rank of the processors whose carthesian coordinates
            ! are "i" and "j" and put in "myid".
            call GRID2D_RANK('R', descla(la_npr_), descla(la_npc_), i, j, myid)
            rank_ip(i + 1, j + 1) = myid*leg_ortho
        end do
    end do

    to_print = .false.
    if (rank .eq. 0 .and. .not. compute_bands) then
        to_print = .true.
    end if

    overhaml = 0.d0
    oversl = 0.d0
    if (ipc .eq. 2) then
        overhamldo = 0.d0
        oversldo = 0.d0
    end if

#else

    overs = 0.d0
    overham = 0.d0
    if (ipc .eq. 2) then
        oversdo = 0.d0
        overhamdo = 0.d0
    end if

    to_print = .false.
    if (rank .eq. 0 .and. .not. compute_bands) then
        to_print = .true.
    end if

#endif

    timep_init = cclock()
    ! initialize variables
    ! define indices of the buffers
    nbas_1 = nelorbu ! # of basis elements (contracted or uncontracted)
    nbas_tot = nelorb ! # of uncontracted basis elements >= nbas_1
    buf_dim = nelorb3 ! leading dimension of the buffer
    indmesh = 0 ! global counter mesh of processor
    buf_dim_conj = nbas_1

    ! allocate buffers for threading
    allocate (wpseudo_and_psip_scratch(lmax*2 + nparshellmax, thread_active))
    wpseudo_and_psip_scratch = 0.d0
    allocate (ivic_scratch(3, max(indt, 1), thread_active))
    ivic_scratch = 0
    allocate (prefactor_scratch((indt - istart + 1), thread_active))
    prefactor_scratch = 0.d0
    allocate (tcost_scratch(max(1, indt), thread_active))
    tcost_scratch = 0.d0
    allocate (distp_scratch(nbas_tot*(indt + 5) + 27*(indt + 1)*max(nshell, nion), thread_active))
    distp_scratch = 0.d0
    allocate (dist_scratch(nion, thread_active), dist_shift_scratch(neigh, thread_active))
    dist_scratch = 0.d0
    dist_shift_scratch = 0.d0
    allocate (r_scratch(nion*(max(indt, 1) + 1), thread_active))
    r_scratch = 0.d0
    allocate (rmu_scratch(3*nion*(max(indt, 1) + 1), thread_active))
    rmu_scratch = 0.d0
    allocate (rmusin_scratch(3*nion*(max(indt, 1) + 1), thread_active))
    rmusin_scratch = 0.d0
    allocate (rmucos_scratch(3*nion*(max(indt, 1) + 1), thread_active))
    rmucos_scratch = 0.d0
    ! buffers required in the size of nbufd
    allocate (x_buffer(3, 0:indt, bufbuf))
    x_buffer = 0.d0
    allocate (angle_buffer(18, bufbuf))
    angle_buffer = 0.d0
    allocate (pshfac1(-nr(1):nr(1)), pshfac2(-nr(2):nr(2))&
            &, pshfac3(0:nr(3)))

    ! total scratch space
    scratch_mem = (1.d0*thread_active*(lmax*2 + nparshellmax + 3*max(indt, 1) &
            & + (indt - istart + 1)*thread_active + max(1, indt) &
            & + nbas_tot*(indt + 5) + 20*(indt + 1)*nshell + nion &
            & + nion*(max(indt, 1) + 1)*10) &
            & + 1.0*bufbuf*(3*(indt + 1) + 10))*8
    if (rank .eq. 0 .and. .not. compute_bands) then
        write (6, '(E12.5, A)') scratch_mem/1d9, " Gbyte per MPI task for threading in initialize_mats_new!"
    end if

    ! allocate buffers
    buf_dim_conj = nbas_1
    call allocate_buffers(bufbuf, nbas_tot, buf_dim, 2*nbas_tot, buf_dim_conj, nbas_tot*(indt + 5), thread_active, .true.)

    ! 1-body Jastrow factor

    ! double_overs = .true. -> compute up/down spin overlaps from scratch. Only when abs(phase_up) != abs(phase_down)
    ! same_phase = .true.   -> up/down spin phases are the same, otherwise they are opposite
    if (to_print) then
        write (6, '(A)') ' Computing overlap matrices:'
        write (6, '(A,L3,L3)') ' double_overs/same_phase: ', double_overs, same_phase
        write (6, *) ' scale one body =', scale_one_body
    end if
    costq0 = 0.d0
    rion_upload(1) = rion_ref(1) - (nx + 1)/2.d0*ax
    rion_upload(2) = rion_ref(2) - (ny + 1)/2.d0*ay
    rion_upload(3) = rion_ref(3) - (nz + 1)/2.d0*az

    if (noread) then

        ! of the orbitals.
        allocate (weightx(nx), weighty(ny), weightz(nz))
        weightx = 1.d0
        weighty = 1.d0
        weightz = 1.d0
        nbuf = 0

        call upload_mat(1, nx, 1, ny, 1, nz, ax, ay, az, .true., memlarge)

        deallocate (weightx, weighty, weightz)

        meshproc = indmesh
        mesh = nx
        mesh = mesh*ny
        mesh = mesh*nz

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

                        if (iespbc) then
                            ncub_min(1, ii + nion) = mod_true(nxi - 1, nx) + 1
                            ncub_max(1, ii + nion) = mod_true(nxf - 1, nx) + 1
                        else
                            ncub_min(1, ii + nion) = nxi
                            ncub_max(1, ii + nion) = nxf
                        end if

                        if (iespbc) then
                            ncub_min(2, ii + nion) = mod_true(nyi - 1, ny) + 1
                            ncub_max(2, ii + nion) = mod_true(nyf - 1, ny) + 1
                        else
                            ncub_min(2, ii + nion) = nyi
                            ncub_max(2, ii + nion) = nyf
                        end if

                        if (iespbc) then
                            ncub_min(3, ii + nion) = mod_true(nzi - 1, nz) + 1
                            ncub_max(3, ii + nion) = mod_true(nzf - 1, nz) + 1
                        else
                            ncub_min(3, ii + nion) = nzi
                            ncub_max(3, ii + nion) = nzf
                        end if

                        nx_proc(ii + nion) = mod(nxf - nxi + 1, nprocrep)
                        nxny_proc(ii + nion) = mod((nyf - nyi + 1)*nx_proc(ii + nion), nprocrep)
                        nx8(ii + nion) = nxf - nxi + 1
                        nxny8(ii + nion) = nx8(ii + nion)*(nyf - nyi + 1)
                        ind_init(ii + nion) = indmesh

                        call upload_mat(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, .false., memlarge)

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

                        if (iespbc) then
                            ncub_min(1, ii) = mod_true(nxi - 1, nx*scalea) + 1
                            ncub_max(1, ii) = mod_true(nxf - 1, nx*scalea) + 1
                        else
                            ncub_min(1, ii) = nxi
                            ncub_max(1, ii) = nxf
                        end if

                        if (iespbc) then
                            ncub_min(2, ii) = mod_true(nyi - 1, ny*scalea) + 1
                            ncub_max(2, ii) = mod_true(nyf - 1, ny*scalea) + 1
                        else
                            ncub_min(2, ii) = nyi
                            ncub_max(2, ii) = nyf
                        end if

                        if (iespbc) then
                            ncub_min(3, ii) = mod_true(nzi - 1, nz*scalea) + 1
                            ncub_max(3, ii) = mod_true(nzf - 1, nz*scalea) + 1
                        else
                            ncub_min(3, ii) = nzi
                            ncub_max(3, ii) = nzf
                        end if

                        nx_proc(ii) = mod(nxf - nxi + 1, nprocrep)
                        nxny_proc(ii) = mod((nyf - nyi + 1)*nx_proc(ii), nprocrep)
                        nx8(ii) = nxf - nxi + 1
                        nxny8(ii) = nx8(ii)*(nyf - nyi + 1)
                        ind_init(ii) = indmesh

                        call upload_mat(nxi, nxf, nyi, nyf, nzi, nzf, axn, ayn, azn, .true., memlarge)

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

                if (iespbc) then
                    ncub_min(1, 2) = mod_true(nxi - 1, nx) + 1
                    ncub_max(1, 2) = mod_true(nxf - 1, nx) + 1
                else
                    ncub_min(1, 2) = nxi
                    ncub_max(1, 2) = nxf
                end if

                if (iespbc) then
                    ncub_min(2, 2) = mod_true(nyi - 1, ny) + 1
                    ncub_max(2, 2) = mod_true(nyf - 1, ny) + 1
                else
                    ncub_min(2, 2) = nyi
                    ncub_max(2, 2) = nyf
                end if

                if (iespbc) then
                    ncub_min(3, 2) = mod_true(nzi - 1, nz) + 1
                    ncub_max(3, 2) = mod_true(nzf - 1, nz) + 1
                else
                    ncub_min(3, 2) = nzi
                    ncub_max(3, 2) = nzf
                end if

                nx_proc(2) = mod(nxf - nxi + 1, nprocrep)
                nxny_proc(2) = mod((nyf - nyi + 1)*nx_proc(2), nprocrep)
                nx8(2) = nxf - nxi + 1
                nxny8(2) = nx8(2)*(nyf - nyi + 1)
                ind_init(2) = indmesh
                call upload_mat(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, .false., memlarge)

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

                if (iespbc) then
                    ncub_min(1, 1) = mod_true(nxi - 1, nx*scalea) + 1
                    ncub_max(1, 1) = mod_true(nxf - 1, nx*scalea) + 1
                else
                    ncub_min(1, 1) = nxi
                    ncub_max(1, 1) = nxf
                end if

                if (iespbc) then
                    ncub_min(2, 1) = mod_true(nyi - 1, ny*scalea) + 1
                    ncub_max(2, 1) = mod_true(nyf - 1, ny*scalea) + 1
                else
                    ncub_min(2, 1) = nyi
                    ncub_max(2, 1) = nyf
                end if

                if (iespbc) then
                    ncub_min(3, 1) = mod_true(nzi - 1, nz*scalea) + 1
                    ncub_max(3, 1) = mod_true(nzf - 1, nz*scalea) + 1
                else
                    ncub_min(3, 1) = nzi
                    ncub_max(3, 1) = nzf
                end if

                nx_proc(1) = mod(nxf - nxi + 1, nprocrep)
                nxny_proc(1) = mod((nyf - nyi + 1)*nx_proc(1), nprocrep)
                nx8(1) = nxf - nxi + 1
                nxny8(1) = nx8(1)*(nyf - nyi + 1)
                ind_init(1) = indmesh

                call upload_mat(nxi, nxf, nyi, nyf, nzi, nzf, axn, ayn, azn, .true., memlarge)

                deallocate (weightx, weighty, weightz)

            end if

            !    Restore original value of volmesh
            volmesh = ax*ay*az*unit_volume
        end if
        meshproc_tot = indmesh

        if (to_print) then
            write (6, *) 'overlap/hamiltonian matrix elements computed'
        end if
        call deallocate_buffers()

        deallocate (wpseudo_and_psip_scratch)
        deallocate (prefactor_scratch, ivic_scratch)
        deallocate (tcost_scratch)
        deallocate (dist_scratch, dist_shift_scratch, distp_scratch)
        deallocate (r_scratch, rmu_scratch, rmusin_scratch, rmucos_scratch)
        deallocate (x_buffer)
        deallocate (angle_buffer)
        deallocate (pshfac1, pshfac2, pshfac3)

        !
        ! Computed quantities :
        !   overs   = <psi_i | psi_j >
        !   overham = <psi_i | H | psi_j >

#ifdef PARALLEL
        ! find the maximum ind among all processors
        call mpi_allreduce(ind, max_index, 1, MPI_INTEGER, MPI_MAX, commrep_mpi, mpi_ierr)
#else
        max_index = ind
#endif
        if (max_index .ne. 0) &
            call error(' initialize_mats_complex ', ' check input nbufd and/or code ', 1, rank)

#ifdef PARALLEL
        call reduce_base_real(1, vlocaltot, commrep_mpi, -1)
        call reduce_base_real(1, costq0, commrep_mpi, -1)
#endif
        vq0tots = 0.5d0*vq0tots ! Unit Hartree
        if (to_print) then
            write (6, *) ' q=0 contribution Ewald (H)  =', &
                -costq0/2.d0, dble(nel)*pi/kappa**2
        end if

#ifdef __SCALAPACK

        ! symmetrize matrices, important to remove roundoff from
        ! block-distributed matrices and improve diagonalization procedure
        if (descla(lambda_node_) > 0) then
            !
            if (ipc .eq. 1) then
                call symmetrize_mat(nbas_1, size(oversl, 1), oversl)
                call symmetrize_mat(nbas_1, size(overhaml, 1), overhaml)
            else
                call symmetrize_mat_complex(nbas_1, size(oversl, 1)/2, oversl)
                call symmetrize_mat_complex(nbas_1, size(overhaml, 1)/2, overhaml)
                if (double_overs) then
                    call symmetrize_mat_complex(nbas_1, size(oversldo, 1)/2, oversldo)
                    call symmetrize_mat_complex(nbas_1, size(overhamldo, 1)/2, overhamldo)
                else
                    oversldo = oversl
                    overhamldo = overhaml
                    if (.not. same_phase) then
                        call conjmat(size(oversldo, 1)/2, size(oversldo, 2), oversldo, size(oversldo, 1)/2)
                        call conjmat(size(overhamldo, 1)/2, size(overhamldo, 2), overhamldo, size(overhamldo, 1)/2)
                    end if
                end if
            end if
        end if

#else

        ! collect matrices distributed over the real-space grid
#ifdef PARALLEL
        call reduce_base_real(size(overs), overs, commrep_mpi, -1)
        call reduce_base_real(size(overham), overham, commrep_mpi, -1)
        if (ipc .eq. 2) then
            call reduce_base_real(size(oversdo), oversdo, commrep_mpi, -1)
            call reduce_base_real(size(overhamdo), overhamdo, commrep_mpi, -1)
        end if
#endif
        ! symmetrization
        if (ipc .eq. 1) then
            call symmetrize_mat(nbas_1, size(overs, 1), overs)
            call symmetrize_mat(nbas_1, size(overham, 1), overham)
        else
            call symmetrize_mat_complex(nbas_1, size(overs, 1)/2, overs)
            call symmetrize_mat_complex(nbas_1, size(overham, 1)/2, overham)
            if (double_overs) then
                call symmetrize_mat_complex(nbas_1, size(oversdo, 1)/2, oversdo)
                call symmetrize_mat_complex(nbas_1, size(overhamdo, 1)/2, overhamdo)
            else
                oversdo = overs
                overhamdo = overham
                if (.not. same_phase) then
                    call conjmat(nbas_1, nbas_1, oversdo, nbas_1)
                    call conjmat(nbas_1, nbas_1, overhamdo, nbas_1)
                end if
            end if
        end if

#endif

        ! write(6,*) 'CHECK oversl sum',rank,sum(oversl(:,:))
        ! write(6,*) 'CHECK oversldo sum',rank,sum(oversldo(:,:))
        ! write(6,*) 'CHECK overhaml sum',rank,sum(overhaml(:,:))
        ! write(6,*) 'CHECK overhamldo sum',rank,sum(overhamldo(:,:))

        ! ------------------- DEBUG ---------------------------

#ifdef DEBUG

#ifdef __SCALAPACK

        if (rankrep .eq. 0) then
            write (6, *) '# Overlap matrix SCALAPACK ', size(oversl), rankcolrep, cell_phase(:)
            do j = 1, descla(nlax_)
                do i = 1, descla(nlax_)
                    write (6, *) j, i, oversl(i, j)
                end do
            end do
            write (6, *) '# Hamiltonian matrix SCALAPACK ', size(overhaml), rankcolrep, cell_phase(:)
            do j = 1, descla(nlax_)
                do i = 1, descla(nlax_)
                    write (6, *) j, i, overhaml(i, j)
                end do
            end do
        end if

#else

        if (rankrep .eq. 0) then
            write (6, *) '# Hamiltonian matrix ', rankcolrep, cell_phase(:)
            do i = 1, nbas_1
                do j = 1, nbas_1
                    write (6, *) i, j, overham(ipc*(i - 1) + 1:ipc*i, j)
                end do
            end do
            write (6, *) '# Overlap matrix ', rankcolrep, cell_phase(:)
            do i = 1, nbas_1
                do j = 1, nbas_1
                    write (6, *) i, j, overs(ipc*(i - 1) + 1:ipc*i, j)
                end do
            end do
        end if

#endif

#endif

        ! ------------------ END DEBUG ------------------

#ifndef __SCALAPACK
        if (write_matrix .and. rank .eq. 0) then
            write (6, *) '# Hamiltonian/Overlap matrix H_0 '
            do i = 1, nbas_1
                do j = 1, nbas_1
                    if (ipc .eq. 2) then
                        write (6, '(2I8,2e14.6,2X,2e14.6)') &
                            i, j, overham(ipc*(i - 1) + 1:ipc*i, j), overs(ipc*(i - 1) + 1:ipc*i, j)
                    else
                        write (6, '(2I8,1e14.6,2X,1e14.6)') &
                            i, j, overham(ipc*(i - 1) + 1:ipc*i, j), overs(ipc*(i - 1) + 1:ipc*i, j)
                    end if
                end do
            end do
        end if
#endif

#ifdef PARALLEL
#ifdef __SCALAPACK

        if (rankrep .eq. 0) then
            if (to_print) write (6, *) ' Raws/Column processors =', descla(la_npr_), descla(la_npc_)
            do indpc = 1, descla(la_npc_) !  loop on column procs
                nc = nrc_ip(indpc)
                nrr = nrc_ip(indpc)
            end do
        end if
        deallocate (irc_ip, nrc_ip, rank_ip)

#endif
        ! make sure all the nodes have released the memory
        call mpi_barrier(MPI_COMM_WORLD, mpi_ierr)
!$omp barrier
#endif

        if (writescratch .eq. 0) then ! do not write anything in case of NonSC run
            !
            ! write overham,wf,overs on scratch files
            ! to speed up continuation. Every processor belonging
            ! to a pool will read/write the corresponding matrix elements.
            !
#ifdef __SCALAPACK
            if (memlarge) then
!    write(6,*) ' rank meshproc write =',rank,size(overhaml),size(oversl),nbas_1,meshproc,meshproc_tot
                write (unit_scratch_distributed) meshproc, meshproc_tot
                write (unit_scratch_distributed) overhaml, oversl, wf(1:ipc*nbas_1, 1:meshproc_tot)

                if (ipc .eq. 2) then
                    if (double_overs) then
                        write (unit_scratch_distributed) overhamldo, oversldo, wf(2*nbas_1 + 1:4*nbas_1, 1:meshproc_tot)
                    else
                        write (unit_scratch_distributed) overhamldo, oversldo
                    end if
                end if
            else
                write (unit_scratch_distributed) overhaml, oversl
                if (ipc .eq. 2) write (unit_scratch_distributed) overhamldo, oversldo
            end if
#else

            if (memlarge) then
                write (unit_scratch_distributed) meshproc, meshproc_tot
                write (unit_scratch_distributed) overham, overs, wf(1:ipc*nbas_1, 1:meshproc_tot)
                if (ipc .eq. 2) then
                    if (double_overs) then
                        write (unit_scratch_distributed) overhamdo, oversdo, wf(2*nbas_1 + 1:4*nbas_1, 1:meshproc_tot)
                    else
                        write (unit_scratch_distributed) overhamdo, oversdo
                    end if
                end if
            else
                write (unit_scratch_distributed) overham, overs
                if (ipc .eq. 2) write (unit_scratch_distributed) overhamdo, oversdo
            end if
#endif

            if (double_mesh) write (unit_scratch_distributed) volmesh_proc(1:meshproc_tot)
        end if

    else ! if(noread) : reading the matrix elements from scratch files

        ! read overham wf overs from files for continuation
#ifdef __SCALAPACK

        if (memlarge) then
            read (unit_scratch_distributed) meshproc, meshproc_tot
!     write(6,*) ' rank meshproc read =',rank,size(overhaml),size(oversl),nbas_1,meshproc,meshproc_tot
            read (unit_scratch_distributed) overhaml, oversl, wf(1:ipc*nbas_1, 1:meshproc_tot)

            if (ipc .eq. 2) then
                if (double_overs) then
                    read (unit_scratch_distributed) overhamldo, oversldo, wf(2*nbas_1 + 1:4*nbas_1, 1:meshproc_tot)
                else
                    read (unit_scratch_distributed) overhamldo, oversldo
                end if
            end if

        else
            read (unit_scratch_distributed) overhaml, oversl
            if (ipc .eq. 2) read (unit_scratch_distributed) overhamldo, oversldo
        end if
#else

        if (memlarge) then
            read (unit_scratch_distributed) meshproc, meshproc_tot
            read (unit_scratch_distributed) overham, overs, wf(1:ipc*nbas_1, 1:meshproc)
            if (ipc .eq. 2) then
                if (double_overs) then
                    read (unit_scratch_distributed) overhamdo, oversdo, wf(2*nbas_1 + 1:4*nbas_1, 1:meshproc_tot)
                else
                    read (unit_scratch_distributed) overhamdo, oversdo
                end if
            else
                read (unit_scratch_distributed) overham, overs
                if (ipc .eq. 2) read (unit_scratch_distributed) overhamdo, oversdo
            end if
        end if
#endif
        if (double_mesh) read (unit_scratch_distributed) volmesh_proc(1:meshproc_tot)
    end if

#if defined (_OPENMP) && defined (__NOOMPDFT)
    call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

#ifdef PARALLEL
    call mpi_barrier(MPI_COMM_WORLD, mpi_ierr) ! syncronize all the groups
!$omp barrier
#endif

    init_time = cclock() - timep_init
    if (to_print) write (6, *) 'Total time initialization =', init_time
    if (allocated(buffer)) call deallocate_buffers()

    return

contains

    subroutine upload_mat(nxi, nxf, nyi, nyf, nzi, nzf, ax, ay, az, add, memlarge)
        implicit none
        integer nxi, ii, jj, kk, i, j, k, nxf, nyi, nyf, nzi, nzf, nx, ny, nz, indtot, nbufrep
        logical add, memlarge
        real*8 ax, ay, az, volmeshu, r_eion
        volmesh = ax*ay*az*unit_volume
        iflagnorm = 3 ! always compute atomic distances and normalization coefficients
        timep_buf = cclock()
        indproc = 0
        ind = 0
        indtot = 0
        nbufrep = nprocrep*bufbuf

        do k = nzi, nzf
            do j = nyi, nyf
                do i = nxi, nxf

                    indtot = indtot + 1

                    if (indproc .eq. rankrep) then

                        ! grid point buffer counter
                        ind = ind + 1
                        indmesh = indmesh + 1
                        if (double_mesh) then
                            weight_buff(ind) = weightx(i - nxi + 1)*weighty(j - nyi + 1)*weightz(k - nzi + 1)
                            if (.not. add) weight_buff(ind) = -weight_buff(ind)
                            volmesh_proc(indmesh) = volmesh*weight_buff(ind)
                        else
                            weight_buff(ind) = 1.d0
                        end if

                        ! buffering the grid point coordinates
                        x_buffer(:, 0, ind) = i*ax*at(:, 1) + j*ay*at(:, 2) + k*az*at(:, 3) + rion_upload(:)

                        !             x_buffer(1,0,ind)=i*ax+rion_upload(1)
                        !             x_buffer(2,0,ind)=j*ay+rion_upload(2)
                        !             x_buffer(3,0,ind)=k*az+rion_upload(3)
                    end if ! endif indproc
                    indproc = indproc + 1
                    if (indproc .eq. nprocrep) indproc = 0

                    if ((i .eq. nxf .and. j .eq. nyf .and. k .eq. nzf) .or. indtot .eq. nbufrep) then

                        if (pseudologic) then
                            do ii = 1, ind
                                call fillmatrix(angle_buffer(1, ii))
                            end do
                        else
                            do ii = 1, ind
                                angle_buffer(:, ii) = angle(:, 1)
                            end do
                        end if
!$omp parallel default(shared) private(tid,ii,costq0_one,vlocaltot_one,volmeshu) reduction(+:costq0,vlocaltot)
#if defined(_OPENMP)
                        tid = 1 + omp_get_thread_num()
#else
                        tid = 1
#endif
                        !OK
!$omp do
                        do ii = 1, ind
                            volmeshu = volmesh*weight_buff(ii)
                            call compute_one_grid_point(x_buffer(1, 0, ii), tid, ii, indmesh - ind + ii &
                                                        , costq0_one, vlocaltot_one, memlarge, volmeshu)
                            costq0 = costq0 + costq0_one
                            vlocaltot = vlocaltot + vlocaltot_one
                        end do
!$omp end do nowait
!$omp end parallel

                        if (contracted_on) then

                            if (ipc .eq. 1) then

                                dgemm_timep = cclock()

                                call dgemm('T', 'N', nbas_1, ind, nbas_tot, 1.d0, mu_c, &
                                           nbas_tot, buffer_on, 2*nbas_tot, 0.d0, buffer, buf_dim)

                                call dgemm('T', 'N', nbas_1, ind, nbas_tot, 1.d0, mu_c, &
                                           nbas_tot, buffer_on(nbas_tot + 1, 1), 2*nbas_tot, 0.d0, buffer(nbas_1 + 1, 1), buf_dim)

                                dgemm_time = dgemm_time + cclock() - dgemm_timep

                            else

                                zgemm_timep = cclock()

                                ! up spin
                                call zgemm('T', 'N', nbas_1, ind, nbas_tot, zone, mu_c, &
                                           nbas_tot, buffer_on, 2*nbas_tot, zzero, buffer, buf_dim)

                                call zgemm('T', 'N', nbas_1, ind, nbas_tot, zone, mu_c, &
                                           nbas_tot, buffer_on(2*nbas_tot + 1, 1), 2*nbas_tot &
                                           , zzero, buffer(2*nbas_1 + 1, 1), buf_dim)

                                ! down spin
                                if (double_overs) then

                                    call zgemm('T', 'N', nbas_1, ind, nbas_tot, zone, mu_c, &
                                               nbas_tot, buffer_on(1, bufbuf + 1), 2*nbas_tot &
                                               , zzero, buffer(1, bufbuf + 1), buf_dim)

                                    call zgemm('T', 'N', nbas_1, ind, nbas_tot, zone, mu_c, &
                                               nbas_tot, buffer_on(2*nbas_tot + 1, bufbuf + 1), 2*nbas_tot, zzero, &
                                               buffer(2*nbas_1 + 1, bufbuf + 1), buf_dim)
                                end if

                                zgemm_time = zgemm_time + cclock() - zgemm_timep

                            end if

                            if (memlarge) call update_wf_memlarge(nbas_1, wf, wf_dim, ind, indmesh)

                        end if
                        !
                        ! To take into account the sign only the first matrix elements are multiplied by it.
                        ! and fill conjugate buffer needed by the matrix elements computation
                        !
                        if (ipc .eq. 2) then
                            do ii = 1, ind
                                buf_conj(1:2*buf_dim_conj, ii) = buffer(1:2*buf_dim_conj, ii)*weight_buff(ii)
                            end do
                            call conjmat(buf_dim_conj, ind, buf_conj, buf_dim_conj)
                            if (double_overs) then
                                do ii = 1, ind
                                    buf_conj(1:2*buf_dim_conj, ii + bufbuf) = buffer(1:2*buf_dim_conj, ii + bufbuf)*weight_buff(ii)
                                end do
                                call conjmat(buf_dim_conj, ind, buf_conj(1, bufbuf + 1), buf_dim_conj)
                            end if
                        else
                            do ii = 1, ind
                                buf_conj(1:buf_dim_conj, ii) = buffer(1:buf_dim_conj, ii)*weight_buff(ii)
                            end do
                        end if

#if defined (_OPENMP) && defined (__NOOMPDFT)
                        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

#ifdef __SCALAPACK

                        if (ipc .eq. 1) then
                            call fill_matrix_scalapack(oversl, buffer, buf_dim, buf_conj&
                      &, buf_dim_conj, ind, irc_ip, nrc_ip, rank_ip, descla(nlax_), 0)
                            call fill_matrix_scalapack(overhaml, buffer, buf_dim, buf_conj&
                      &, buf_dim_conj, ind, irc_ip, nrc_ip, rank_ip, descla(nlax_), nbas_1)
                        else
                            ! up spin
                            call fill_matrix_scalapack_complex(oversl, buffer, buf_dim, buf_conj, buf_dim_conj, &
                                                               ind, irc_ip, nrc_ip, rank_ip, descla(nlax_), 0)
                            call fill_matrix_scalapack_complex(overhaml, buffer, buf_dim, buf_conj, buf_dim_conj, &
                                                               ind, irc_ip, nrc_ip, rank_ip, descla(nlax_), nbas_1)

                            ! down spin
                            if (double_overs) then
                                call fill_matrix_scalapack_complex(oversldo, buffer(1, bufbuf + 1) &
                                                                   , buf_dim, buf_conj(1, bufbuf + 1), buf_dim_conj, &
                                                                   ind, irc_ip, nrc_ip, rank_ip, descla(nlax_), 0)
                                call fill_matrix_scalapack_complex(overhamldo, buffer(1, bufbuf + 1) &
                                                                   , buf_dim, buf_conj(1, bufbuf + 1), buf_dim_conj, &
                                                                   ind, irc_ip, nrc_ip, rank_ip, descla(nlax_), nbas_1)
                            end if
                        end if
#else
                        ! Computing overlaps: overs/oversdo and overham/overhamdo
                        if (ipc .eq. 1) then
                            dgemm_timep = cclock()
                            call dgemm('N', 'T', nbas_1, nbas_1, ind, volmesh, buf_conj, &
                                       buf_dim_conj, buffer, buf_dim, 1.d0, overs, nbas_1)
                            call dgemm('N', 'T', nbas_1, nbas_1, ind, volmesh, buf_conj, &
                                       buf_dim_conj, buffer(nbas_1 + 1, 1), buf_dim, 1.d0, overham, nbas_1)
                            dgemm_time = dgemm_time + cclock() - dgemm_timep
                        else
                            ! up spin
                            zgemm_timep = cclock()
                            call zgemm('N', 'T', nbas_1, nbas_1, ind, dcmplx(volmesh), buf_conj, &
                                       buf_dim_conj, buffer, buf_dim, zone, overs, nbas_1)
                            call zgemm('N', 'T', nbas_1, nbas_1, ind, dcmplx(volmesh), buf_conj, &
                                       buf_dim_conj, buffer(2*nbas_1 + 1, 1), buf_dim, zone, overham, nbas_1)
                            ! down spin
                            if (double_overs) then
                                call zgemm('N', 'T', nbas_1, nbas_1, ind, dcmplx(volmesh), buf_conj(1, bufbuf + 1), &
                                           buf_dim_conj, buffer(1, bufbuf + 1), buf_dim, zone, oversdo, nbas_1)
                                call zgemm('N', 'T', nbas_1, nbas_1, ind, dcmplx(volmesh), buf_conj(1, bufbuf + 1), &
                                           buf_dim_conj, buffer(2*nbas_1 + 1, bufbuf + 1), buf_dim, zone, overhamdo, nbas_1)
                            end if
                            zgemm_time = zgemm_time + cclock() - zgemm_timep
                        end if
#endif
                        nbuf = nbuf + 1
                        single_buffer_time = cclock() - timep_buf

                        if (to_print) write (6, '(A, I5, A, F8.2, A)') ' Buf number = ', &
                            nbuf, ' taking ', single_buffer_time, ' sec.'

                        timep_buf = timep_buf + single_buffer_time
                        ind = 0
                        indtot = 0
                    end if ! endif for ind global
                end do
            end do
        end do ! enddo over the mesh

    end subroutine upload_mat

    subroutine compute_one_grid_point(x, tid, ind, indmesh, costq0, vlocaltot &
                                      , memlarge, volmesh)
        ! This routine computes the WF on a single grid point.
        ! it is thread-safe and intended to be called from threaded region.
        use allio, only: kappa, iespbc, nion, rion, zetar, neigh, x_neigh, iflagerr, metric, applypbc
        use dielectric
        implicit none
        ! arguments
        ! single electron positions
        real(8), intent(inout) :: x(3, 0:indt)
        ! index to the scratch space
        integer, intent(in) :: tid
        ! index in the grid buffer
        integer, intent(in) :: ind
        ! global index of the current grid point
        integer, intent(in) :: indmesh
        real(8), intent(in) :: volmesh
        ! return values
        real(8), intent(out) :: costq0, vlocaltot
        !
        ! local variables
        real(8) :: pseudolocalnew
        real(8) :: jastrow, jas_1body, r0, grad(3), lap, tgrad(3), tlap, rc(3), hess(3)
        real(8), external :: jastrow_ei
        integer :: indtm, ii, jj
        logical memlarge
        real(8) vec(3), vpots, const, dist_kel(3), ewald_energy, derfc
        integer kk
        vec = 0.d0 !
        vpots = 0.d0
        costq0 = 0.d0
        vlocaltot = 0.d0
        if (iespbc) then
            do kk = 1, nion
                dist_kel(1) = x(1, 0) - rion(1, kk)
                dist_kel(2) = x(2, 0) - rion(2, kk)
                dist_kel(3) = x(3, 0) - rion(3, kk)
                call ApplyPBC(dist_kel, 1)
                dist_scratch(kk, tid) = max(dsqrt(dist_kel(1)**2 + dist_kel(2)**2 + dist_kel(3)**2), mindist)
                dist_shift_scratch(1, tid) = -2.d0*zetar(kk)*rep_erfc(dist_scratch(kk, tid), kappa)
#ifdef _SIMD
!$omp simd
#endif
                do ii = 2, neigh
!               dist_kel(1) = dist_kel(1)+x_neigh(ii,1)
!               dist_kel(2) = dist_kel(2)+x_neigh(ii,2)
!               dist_kel(3) = dist_kel(3)+x_neigh(ii,3)

                    dist_shift_scratch(ii, tid) = max(dsqrt((dist_kel(1) + x_neigh(ii, 1))**2 &
                                                            + (dist_kel(2) + x_neigh(ii, 2))**2 &
                                                            + (dist_kel(3) + x_neigh(ii, 3))**2), mindist)

!               if(dist_scratch(kk, tid).lt.mindist) dist_scratch(kk, tid) = mindist
!               const = -2.d0 * zetar(kk) * derfc(kappa * dist_scratch(kk, tid)) / dist_scratch(kk, tid)
!               const =
                    dist_shift_scratch(ii, tid) = -2.d0*zetar(kk)*rep_erfc(dist_shift_scratch(ii, tid), kappa)
                end do
                vpots = vpots + sum(dist_shift_scratch(1:neigh, tid))
            end do
            ! Ewald contribution after appropriate initialization
            costq0 = costq0 + volmesh*vpots
            call Ewaldup1b(x, tid, ewald_energy)
            vpots = vpots + ewald_energy + ewaldion1b + (1.d0 - weightvh)*ewaldel1b
        else
            do kk = 1, nion
                dist_scratch(kk, tid) = dsqrt((x(1, 0) - rion(1, kk))**2 + &
                                              (x(2, 0) - rion(2, kk))**2 + &
                                              (x(3, 0) - rion(3, kk))**2)
                if (dist_scratch(kk, tid) .lt. mindist) dist_scratch(kk, tid) = mindist
!               const = -2.d0 * zetar(kk) / dist_scratch(kk, tid)
                const = -2.d0*zetar(kk)*veps(dist_scratch(kk, tid))
                vpots = vpots + const
            end do
        end if
        !
        ! evaluate pseudopotential contribution
        !
        if (npsa .gt. 0) then

            indtm = istart - 1

            call pseudoset(tid, x(1, 0), ivic_scratch(1, istart, tid), prefactor_scratch, pseudolocalnew &
                    &, thread_active, dist_scratch(1, tid), rion, wpseudo_and_psip_scratch(nparshellmax + 1, tid) &
                    &, npsa, indt - istart + 1, indtm, .false., angle_buffer(1, ind) &
                    &, wpseudo_and_psip_scratch(1, tid), Lbox, mindist, iflagerr, vec)

            vpots = vpots + 2.d0*pseudolocalnew

            vlocaltot = vlocaltot + pseudolocalnew*volmesh

            !if(rank.eq.0) write(6,'(A, L2, 2I10, 3F15.10)') "ivic_debug "&
            !, pseudologic, indmesh, tid, ivic(1,istart+1,tid), ivic(2,istart+1,tid), ivic(3,istart+1,tid)
            ! definition mesh for pseudo
            do ii = 1, indtm
                do kk = 1, 3
                    x(kk, ii) = x(kk, 0) + ivic_scratch(kk, ii, tid)
                end do
            end do
            ! computation coefficients put in tcost_scratch
            do kk = 1, indtm
                !          tcost_scratch(kk,tid)=t(x(1,0),nion,rion,zetar,kk,alat&
                !               &,ivic_scratch,plat,prefactor_scratch(1,tid),itest,tid,indt-istart+1,LBox&
                !               &,istart)
                tcost_scratch(kk, tid) = -prefactor_scratch(kk - istart + 1, tid)
            end do

            !if(rank.eq.0) write(6,'(A, I10, 5F15.10)') "tcost_debug "&
            !, indmesh, tcost_scratch(1,tid), tcost_scratch(2,tid), sum(tcost_scratch(1:indtm,tid)), vpots, pseudolocalnew
        else
            indtm = 0
        end if
        !
        ! evaluating basis set for each electron position on the mesh
        !
        call upnewwf(indt, 0, indtm, 0, nshell, ioptorb, ioccup, x, 1, r_scratch(1, tid), rmu_scratch(1, tid), &
                     dupr, zetar, rion, distp_scratch(1, tid), wf_threading_scratch(1, tid), nbas_tot, nion, kion &
                     , iflagnorm, cnorm, LBox, rmucos_scratch(1, tid), rmusin_scratch(1, tid), mindist, indpar_tab &
                     , indorb_tab, indshell_tab, .true.)
        if (double_overs) then
            call upnewwf(indt, 0, indtm, 0, nshell, ioptorb, ioccup, x, 1, r_scratch(1, tid), rmu_scratch(1, tid), &
                         dupr, zetar, rion, distp_scratch(1, tid), wf_threading_scratch_down(1, tid), nbas_tot, nion, kion &
                         , iflagnorm, cnorm, LBox, rmucos_scratch(1, tid), rmusin_scratch(1, tid), mindist, indpar_tab &
                         , indorb_tab, indshell_tab, .false.)
        end if

        !1body Jastrow
        if (n_body_on .ne. 0) then
            ! updating also the pseudo if present
!$omp parallel do default(shared) private(ii,jj,kk,jas_1body,rc,r0) schedule(static)
            do ii = 1, indtm
                jas_1body = -scale_one_body
                do jj = 1, nion
                    if (iespbc) then

                        !rc(:)=x(:,ii)-rion(:,jj)
                        rc(1) = x(1, ii) - rion(1, jj)
                        rc(2) = x(2, ii) - rion(2, jj)
                        rc(3) = x(3, ii) - rion(3, jj)

!                       call CartesianToCrystal(rc, 1)
                        rc(:) = car2cry(:, 1)*rc(1) + car2cry(:, 2)*rc(2) + car2cry(:, 3)*rc(3)

                        do kk = 1, 3
                            rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                        end do
                        r0 = norm_metric(rc, metric)

                    else

                        !                  rc(:)=(x(:,ii)-rion(:,jj))*costz(jj)
                        rc(1) = (x(1, ii) - rion(1, jj))*costz(jj)
                        rc(2) = (x(2, ii) - rion(2, jj))*costz(jj)
                        rc(3) = (x(3, ii) - rion(3, jj))*costz(jj)
                        r0 = dsqrt(sum(rc(:)**2))

                    end if
                    jas_1body = jas_1body - &
                            &jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                end do
                jas_1body = dexp(jas_1body)
                do iii = ipc*nbas_tot*ii + 1, ipc*nbas_tot*(ii + 1)
                    wf_threading_scratch(iii, tid) = wf_threading_scratch(iii, tid)*jas_1body
                    if (double_overs) wf_threading_scratch_down(iii, tid) = wf_threading_scratch_down(iii, tid)*jas_1body
                end do
            end do
!$omp end parallel do

            ! updating gradients/u  and laplacian/u for 1b Jastrow
            tlap = 0.d0
            tgrad = 0.d0
            jas_1body = -scale_one_body
            do ii = 1, nion
                if (iespbc) then

                    !rc(:)=x(:,0)-rion(:,ii)
                    rc(1) = x(1, 0) - rion(1, ii)
                    rc(2) = x(2, 0) - rion(2, ii)
                    rc(3) = x(3, 0) - rion(3, ii)

                    call jastrowgrad_pbc(rc, vj(pointvj(1, ii))&
                            &, pointvj(2, ii), jastrow, grad, lap, -1, costz(ii))
                else

                    !               rc(:)=(x(:,0)-rion(:,ii))*costz(ii)
                    rc(1) = (x(1, 0) - rion(1, ii))*costz(ii)
                    rc(2) = (x(2, 0) - rion(2, ii))*costz(ii)
                    rc(3) = (x(3, 0) - rion(3, ii))*costz(ii)

                    call jastrowgrad(rc, vj(pointvj(1, ii)), pointvj(2, ii), jastrow, grad, lap, -1)
                end if
                ! to avoid overflow or too small numbers
                jas_1body = jas_1body - jastrow*costz3(ii)
                tlap = tlap - lap*costz(ii)**2*costz3(ii)
                tgrad(1) = tgrad(1) - grad(1)*costz(ii)*costz3(ii)
                tgrad(2) = tgrad(2) - grad(2)*costz(ii)*costz3(ii)
                tgrad(3) = tgrad(3) - grad(3)*costz(ii)*costz3(ii)

            end do
            jas_1body = dexp(jas_1body)

            if (ipc .eq. 1) then
                call up_1body(wf_threading_scratch(1, tid), nbas_tot, indt, tlap, tgrad, jas_1body)
            else
                call up_1body_complex(wf_threading_scratch(1, tid), nbas_tot, indt, tlap, tgrad, jas_1body) ! up spin
                if (double_overs) &
                    call up_1body_complex(wf_threading_scratch_down(1, tid), nbas_tot, indt, tlap, tgrad, jas_1body) ! down spin
            end if
        end if ! end if(n_body_on.ne.0)
        !
        ! defining buffers for the |psi> and H|psi> :
        ! buffer = for uncontracted basis
        ! buffer_on = for contracted basis
        ! indexing for each electron type : 1:nbas_tot = w.f.
        !                                   nbas_1+1:2*nbas_tot = hamiltonian
        ! buffer distribution: (:,1:bufbuf) = up spin electrons
        !                      (:,bufbuf+1,2*bufbuf) = down spin electrons
        ! nbas_1+1=nelorb+1
        !
        !OK
        if (contracted_on) then

            if (ipc .eq. 1) then

                call dcopy(nbas_tot, wf_threading_scratch(1, tid), 1, buffer_on(1, ind), 1)
                call fillbuff(nbas_tot, indt, istart, indtm, tcost_scratch(1, tid), vpots, wf_threading_scratch(1, tid)&
                        &, buffer_on(nbas_tot + 1, ind))

            else

                ! fill buffer for up spin electrons
                call zcopy(nbas_tot, wf_threading_scratch(1, tid), 1, buffer_on(1, ind), 1)
                call fillbuff_complex(nbas_tot, indt, istart, indtm, tcost_scratch(1, tid), vpots, wf_threading_scratch(1, tid)&
                        &, buffer_on(2*nbas_tot + 1, ind))

                ! fill buffer for down spin electrons and optimize if same or opposite phase
                if (double_overs) then

                    call zcopy(nbas_tot, wf_threading_scratch_down(1, tid), 1, buffer_on(1, bufbuf + ind), 1)
                    call fillbuff_complex(nbas_tot, indt, istart, indtm, tcost_scratch(1, tid), vpots &
                                          , wf_threading_scratch_down(1, tid) &
                                          , buffer_on(2*nbas_tot + 1, bufbuf + ind))

                end if
            end if

        else

            if (ipc .eq. 1) then

                call dcopy(nbas_tot, wf_threading_scratch(1, tid), 1, buffer(1, ind), 1)
                call fillbuff(nbas_tot, indt, istart, indtm, tcost_scratch(1, tid), vpots, wf_threading_scratch(1, tid), &
                              buffer(nbas_1 + 1, ind))
            else

                call zcopy(nbas_tot, wf_threading_scratch(1, tid), 1, buffer(1, ind), 1)
                call fillbuff_complex(nbas_tot, indt, istart, indtm, &
                                      tcost_scratch(1, tid), vpots, wf_threading_scratch(1, tid), buffer(2*nbas_1 + 1, ind))
                !
                ! fill the down spin part of the overlap and optimize in case of same or opposite phase
                !
                if (double_overs) then
                    call zcopy(nbas_tot, wf_threading_scratch_down(1, tid), 1, buffer(1, bufbuf + ind), 1)
                    call fillbuff_complex(nbas_tot, indt, istart, indtm, tcost_scratch(1, tid), &
                                          vpots, wf_threading_scratch_down(1, tid), buffer(2*nbas_1 + 1, bufbuf + ind))
                end if

            end if
            if (memlarge) call update_wf_memlarge(nbas_tot, wf, wf_dim, ind, indmesh, tid)

        end if

        return
    end subroutine compute_one_grid_point

end subroutine initialize_mats_new
subroutine set_interval(scalea, nx0r, ny0r, nz0r, rion, rion_upload&
        &, ax, ay, az, nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr, changer, weightx&
        &, weighty, weightz)
    use allio, only: iespbc, nx, ny, nz
    use cell, only: CartesianToCrystal
    implicit none
    integer scalea, i, nx0, ny0, nz0, nxi, nxf, nyi, nyf, nzi, nzf, nxr, nyr, nzr&
            &, nx0r, ny0r, nz0r
    real*8 rion(3), rion_upload(3), rion_diff(3), ax, ay, az
    real*8 weightx(*), weighty(*), weightz(*)
    logical changer, nofit(3)

    nx0 = nx0r
    ny0 = ny0r
    nz0 = nz0r

    rion_diff(:) = rion(:) - rion_upload(:)
    if (iespbc) call CartesianToCrystal(rion_diff, 1)

    if (changer) then
        nxr = nint(rion_diff(1)/ax)
        nyr = nint(rion_diff(2)/ay)
        nzr = nint(rion_diff(3)/az)
        !   The total number of points in the mesh including the extra ones for interpolation at the boarder.
        nx0 = nx0 + 2
        ny0 = ny0 + 2
        nz0 = nz0 + 2
    else
        nxr = nxr*scalea
        nyr = nyr*scalea
        nzr = nzr*scalea

        ax = ax/scalea
        ay = ay/scalea
        az = az/scalea

        nx0 = (nx0 - 1)*scalea + 3
        ny0 = (ny0 - 1)*scalea + 3
        nz0 = (nz0 - 1)*scalea + 3

    end if

    if (changer) then
        if (nx0 .gt. nx) then
            nxi = 1
            nxf = nx
            nofit(1) = .true.
        else
            nofit(1) = .false.
            if (mod(nx0, 2) .eq. 0) then
                if (rion_diff(1) - nxr*ax .gt. 0) then
                    nxi = nxr - nx0/2 + 1
                    nxf = nxr + nx0/2
                else
                    nxi = nxr - nx0/2
                    nxf = nxr + nx0/2 - 1
                end if
            else
                nxi = nxr - nx0/2
                nxf = nxr + nx0/2
            end if
            if (.not. iespbc .and. nxi .lt. 1) nxi = 1
            if (.not. iespbc .and. nxf .gt. nx) nxf = nx
        end if
    else
        if (nxf - nxi + 1 .eq. nx) then
            nofit(1) = .true.
            nxi = 1
            nxf = scalea*nx
        else
            nofit(1) = .false.
            nxi = (nxi + 1)*scalea - 1
            nxf = (nxf - 1)*scalea + 1
        end if
    end if
    if (nofit(1)) then
        weightx(nxi:nxf) = 1.d0
    else
        call prepw(nx0, weightx)
    end if

    if (changer) then
        if (ny0 .gt. ny) then
            nofit(2) = .true.
            nyi = 1
            nyf = ny
        else
            nofit(2) = .false.
            if (mod(ny0, 2) .eq. 0) then
                if (rion_diff(2) - nyr*ay .gt. 0) then
                    nyi = nxr - ny0/2 + 1
                    nyf = nxr + ny0/2
                else
                    nyi = nyr - ny0/2
                    nyf = nyr + ny0/2 - 1
                end if
            else
                nyi = nyr - ny0/2
                nyf = nyr + ny0/2
            end if
            if (.not. iespbc .and. nyi .lt. 1) nyi = 1
            if (.not. iespbc .and. nyf .gt. ny) nyf = ny
        end if
    else
        if (nyf - nyi + 1 .eq. ny) then
            nofit(2) = .true.
            nyi = 1
            nyf = scalea*ny
        else
            nofit(2) = .false.
            nyi = (nyi + 1)*scalea - 1
            nyf = (nyf - 1)*scalea + 1
        end if
    end if
    if (nofit(2)) then
        weighty(nyi:nyf) = 1.d0
    else
        call prepw(ny0, weighty)
    end if
    if (changer) then
        if (nz0 .gt. nz) then
            nofit(3) = .true.
            nzi = 1
            nzf = nz
        else
            nofit(3) = .false.
            if (mod(nz0, 2) .eq. 0) then
                if (rion_diff(3) - nzr*az .gt. 0) then
                    nzi = nzr - nz0/2 + 1
                    nzf = nzr + nz0/2
                else
                    nzi = nzr - nz0/2
                    nzf = nzr + nz0/2 - 1
                end if
            else
                nzi = nzr - nz0/2
                nzf = nzr + nz0/2
            end if
            if (.not. iespbc .and. nzi .lt. 1) nzi = 1
            if (.not. iespbc .and. nzf .gt. nz) nzf = nz
        end if
    else
        if (nzf - nzi + 1 .eq. nz) then
            nofit(3) = .true.
            nzi = 1
            nzf = scalea*nz
        else
            nofit(3) = .false.
            nzi = (nzi + 1)*scalea - 1
            nzf = (nzf - 1)*scalea + 1
        end if
    end if
    if (nofit(3)) then
        weightz(nzi:nzf) = 1.d0
    else
        call prepw(nz0, weightz)
    end if
    return
end
subroutine prepw(n, w)
    integer n, i
    real*8 w(n), fb(4), fe(4)
    fb(1) = -1.d0/24.d0
    fb(2) = 13.d0/24.d0
    fb(3) = 13.d0/24.d0
    fb(4) = -1.d0/24.d0
    fe(1) = 9.d0/24d0
    fe(2) = 19.d0/24d0
    fe(3) = -5.d0/24.d0
    fe(4) = 1.d0/24.d0

    !     if(add) then
    !      if(n.eq.1) then
    !      w(1)=1.d0
    !      elseif(n.eq.2) then
    !      w(1:2)=0.5d0
    !      elseif(n.eq.3) then
    !      w(1)=1.d0/3.d0
    !      w(2)=4.d0/3.d0
    !      w(3)=1.d0/3.d0
    !      elseif(n.eq.4) then
    !      w(1)=3.d0/8.d0
    !      w(4)=3.d0/8.d0
    !      w(2)=9.d0/8.d0
    !      w(3)=9.d0/8.d0
    !      elseif(n.eq.5) then
    !      w(1)=fe(1)+fb(1)
    !      w(2)=fe(2)+fb(2)+fb(1)+fe(4)
    !      w(3)=fe(3)+fb(3)+fb(2)+fe(3)
    !      w(4)=fe(4)+fb(4)+fb(3)+fe(2)
    !      w(5)=fb(4)+fe(1)
    !      elseif(n.eq.6) then
    !      w(1)=fe(1)+fb(1)
    !      w(2)=fe(2)+fb(2)+fb(1)
    !      w(3)=fe(3)+fb(3)+fb(2)+fb(1)+fe(4)
    !      w(4)=fe(4)+fb(4)+fb(3)+fb(2)+fe(3)
    !      w(5)=fb(4)+fb(3)+fe(2)
    !      w(6)=fb(4)+fe(1)
    !      elseif(n.eq.7) then
    !      w(1)=fe(1)+fb(1)
    !      w(2)=fe(2)+fb(2)+fb(1)
    !      w(3)=fe(3)+fb(3)+fb(2)+fb(1)
    !      w(4)=fe(4)+fb(4)+fb(3)+fb(2)+fb(1)+fe(4)
    !      w(5)=fb(4)+fb(3)+fb(2)+fe(3)
    !      w(6)=fb(4)+fb(3)+fe(2)
    !      w(7)=fb(4)+fe(1)
    !      elseif(n.gt.7) then
    !      w(1)=fe(1)+fb(1)
    !      w(2)=fe(2)+fb(2)+fb(1)
    !      w(3)=fe(3)+fb(3)+fb(2)+fb(1)
    !      w(4)=fe(4)+fb(4)+fb(3)+fb(2)+fb(1)
    !      do i=5,n-4
    !      w(i)=fb(4)+fb(3)+fb(2)+fb(1)
    !      enddo
    !      w(n)=fb(4)+fe(1)
    !      w(n-1)=fb(4)+fb(3)+fe(2)
    !      w(n-2)=fb(4)+fb(3)+fb(2)+fe(3)
    !      w(n-3)=fb(4)+fb(3)+fb(2)+fb(1)+fe(4)
    !      endif
    !     else
    if (n .eq. 1) then
        w(1) = 1.d0
    elseif (n .eq. 2) then
        w(1:2) = 0.5d0
    elseif (n .eq. 3) then
        w(1) = 1.d0/3.d0
        w(2) = 4.d0/3.d0
        w(3) = 1.d0/3.d0
    elseif (n .eq. 4) then
        w(1:4) = fb(1:4)
    elseif (n .eq. 5) then
        w(1) = fb(1)
        w(2) = fb(2) + fb(1)
        w(3) = fb(3) + fb(2)
        w(4) = fb(4) + fb(3)
        w(5) = fb(4)
    elseif (n .eq. 6) then
        w(1) = fb(1)
        w(2) = fb(2) + fb(1)
        w(3) = fb(3) + fb(2) + fb(1)
        w(4) = fb(4) + fb(3) + fb(2)
        w(5) = fb(4) + fb(3)
        w(6) = fb(4)
    elseif (n .ge. 7) then
        w(1) = fb(1)
        w(2) = fb(2) + fb(1)
        w(3) = fb(3) + fb(2) + fb(1)
        do i = 4, n - 3
            w(i) = fb(4) + fb(3) + fb(2) + fb(1)
        end do
        w(n - 2) = fb(4) + fb(3) + fb(2)
        w(n - 1) = fb(4) + fb(3)
        w(n) = fb(4)
    end if
    !     endif

    return
end
