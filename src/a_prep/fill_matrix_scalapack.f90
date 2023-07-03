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

! These subroutines update a block matrices in the SCALAPACK formalism
! after one buffering over the MOs grid. Real version is followed by
! the complex version.
! They are used in the routines: initialize_mats_complex, uphamilt_kpoints

subroutine fill_matrix_scalapack(mat_in, buffer, nelorb, buffer_weight, nelorb_w, bufbuf, &
                                 irc_ip, nrc_ip, rank_ip, nlax, adr_buf)

#ifdef __SCALAPACK

    use constants, only: zzero
    use descriptors
    use allio, only: ortho_cntx, me_blacs, np_ortho, me_ortho, &
                     ortho_comm, ortho_comm_id, rank, commrep_mpi
    use setup, only: volmesh
    implicit none
    ! input
    integer, intent(in) :: nelorb, bufbuf, irc_ip(np_ortho(1)), nrc_ip(np_ortho(1)), &
                           rank_ip(np_ortho(1), np_ortho(2)), adr_buf, nlax, nelorb_w

    ! REAL VERSION
    ! nelorb  = leading dimension of the buffer as it has been allocated in the calling routine
    ! descla( nlax_ ) = leading dimension of the block matrix
    ! adr_buf = index to be added to the address of the buffers buffer/buf_conj
    !           in order to access its different parts
    ! This trick must be done since the BLAS on the block matrix needs to know
    ! the exact leading dimension of the buffer.

    real(8), intent(in) :: buffer(nelorb, *), buffer_weight(nelorb_w, *)
    real(8), intent(out) :: mat_in(nlax, nlax)
    ! local
    integer :: irow, icol, nr, ir, root, ierr, ic, nc, ii
    real(8), dimension(:, :), allocatable :: work, work2
    !
#ifdef PARALLEL
    include "mpif.h"
#endif
    !
    ! allocate workspace
    allocate (work(nlax, nlax), work2(nlax, nlax))
    work = 0.d0
    work2 = 0.d0
    !
#ifdef _OFFLOAD
!$omp target data map(to:buffer(:,1:bufbuf)&
!$omp& ,buffer_weight(:,1:bufbuf)) map(to:work)
#endif
    do icol = 1, descla(la_npc_) !  loop on column procs
        nc = nrc_ip(icol)
        ic = irc_ip(icol) + adr_buf
        do irow = 1, descla(la_npr_)
            nr = nrc_ip(irow)
            ir = irc_ip(irow)
            !  rank of the processor for which this block (irow,icol) is destinated
            root = rank_ip(irow, icol)
            ! use blas subs. on the matrix block
            if (nr .gt. 0 .and. nc .gt. 0) then
                call dgemm_('N', 'T', nr, nc, bufbuf, volmesh, &
                            buffer_weight(ir, 1), nelorb_w, buffer(ic, 1), nelorb, 0.d0, work, nlax)
#ifdef _OFFLOAD
!$omp target update from(work)
#endif
            end if
            ! accumulate result on dm of root proc.
#ifdef UNREL
            call reduce_base_real_to(size(work), work, work2, commrep_mpi, -1)
            if (rank .eq. root) mat_in = work2 + mat_in
#else
            call reduce_base_real_to(size(work), work, work2, commrep_mpi, root)
#endif
        end do
    end do
#ifdef _OFFLOAD
!$omp end target data
#endif

#ifndef UNREL
    if (descla(lambda_node_) > 0) then
        mat_in = work2 + mat_in
    end if
#endif
    deallocate (work, work2)

    return

#else
    return
#endif

end subroutine fill_matrix_scalapack

subroutine fill_matrix_scalapack_complex(mat_in, buffer, nelorb, buf_conj, nelorb_conj, &
                                         bufbuf, irc_ip, nrc_ip, rank_ip, nlax, adr_buf)

#ifdef __SCALAPACK

    use constants, only: zzero
    use descriptors
    use allio, only: ortho_cntx, me_blacs, np_ortho, me_ortho, &
                     ortho_comm, ortho_comm_id, commrep_mpi, rankrep, rank
    use setup, only: volmesh
    implicit none
    ! input
    integer, intent(in) :: nelorb, bufbuf, irc_ip(np_ortho(1)), nrc_ip(np_ortho(1)), &
                           rank_ip(np_ortho(1), np_ortho(2)), adr_buf, nlax, nelorb_conj

    ! COMPLEX VERSION
    ! nelorb  = leading dimension of the buffer as it has been allocated in the calling routine
    ! descla( nlax_ ) = leading dimension of the block matrix
    ! adr_buf = index to be added to the address of the buffers buffer/buf_conj
    !           in order to access its different parts
    ! This trick must be done since the BLAS on the block matrix needs to know
    ! the exact leading dimension of the buffer.

    complex(8), intent(in) :: buffer(nelorb, *), buf_conj(nelorb_conj, *)
    complex(8), intent(out) :: mat_in(nlax, nlax)
    ! local
    integer :: irow, icol, nr, ir, root, ierr, ic, nc, ii
    complex(8), dimension(:, :), allocatable :: work, work2
    !

#ifdef PARALLEL
    include "mpif.h"
#endif
    !
    ! allocate workspace
    allocate (work(nlax, nlax), work2(nlax, nlax))
    work = zzero
    work2 = zzero
    !
#ifdef _OFFLOAD
!$omp target data map(to:buffer(:,1:bufbuf),buf_conj(:,1:bufbuf)) map(to:work)
#endif
    do icol = 1, descla(la_npc_) !  loop on column procs
        nc = nrc_ip(icol)
        ic = irc_ip(icol) + adr_buf
        do irow = 1, descla(la_npr_)
            nr = nrc_ip(irow)
            ir = irc_ip(irow)
            !  rank of the processor for which this block (irow,icol) is destinated
            root = rank_ip(irow, icol)
            ! use blas subs. on the matrix block
            if (nr .gt. 0 .and. nc .gt. 0) then
                call zgemm_('N', 'T', nr, nc, bufbuf, dcmplx(volmesh), &
                            buf_conj(ir, 1), nelorb_conj, buffer(ic, 1), nelorb, zzero, work, nlax)
#ifdef _OFFLOAD
!$omp target update from(work)
#endif
            end if
            ! accumulate result on dm of root proc.
#ifdef UNREL
            call reduce_base_complex_to(size(work), work, work2, commrep_mpi, -1)
            if (rankrep .eq. root) mat_in = work2 + mat_in
#else
            call reduce_base_complex_to(size(work), work, work2, commrep_mpi, root)
#endif
        end do
    end do
#ifdef _OFFLOAD
!$omp end target data
#endif

#ifndef UNREL
    if (descla(lambda_node_) > 0) then
        mat_in = work2 + mat_in
    end if
#endif
    deallocate (work, work2)

    return

#else
    return
#endif

end subroutine fill_matrix_scalapack_complex

