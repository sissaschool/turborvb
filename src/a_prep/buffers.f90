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

module buffers

    ! This module contains a couple of subroutines
    ! for allocating and deallocating the buffers needed
    ! when dealing with real space grid calculations (subroutines
    ! initialize_matrices, upham_new, elec_dens).

    use constants, only: ipc, zone, zzero
    use allio, only: rank, compute_bands
    use setup, only: contracted_on, double_overs, yeslsda

    real(8), dimension(:, :), allocatable :: wf_threading_scratch, wf_threading_scratch_down
    real(8), dimension(:, :), allocatable :: buffer, buffer_on, buf_conj
    real(8), dimension(:), allocatable :: weight_buff
contains

    subroutine allocate_buffers(bufbuf, & ! number of grid points to be buffered
                                nbas_tot, & ! total atomic basis set dimension
                                buf_dim, & ! buffer leading dimension
                                buf_dim_contr, & ! buffer_on leading dimension
                                buf_dim_conj, & ! buf_conj leading dimension
                                wf_dim, & ! dimension of the wf_threading_scratch buffer
                                thread_active, & ! the number of active threads
                                print_size) ! output the final allocated size
        implicit none

        ! input dimensions
        integer, intent(in) :: bufbuf, nbas_tot, buf_dim, buf_dim_contr, &
                               buf_dim_conj, wf_dim, thread_active
        ! print switch
        logical, intent(in), optional :: print_size

        ! local variables
        integer :: buf_num
        ! final allocated size
        real(8) :: allocated_size_buf, allocated_size_thread

        allocated_size_buf = 0
        allocated_size_thread = 0

        if (ipc .eq. 2 .and. (double_overs .or. yeslsda)) then
            buf_num = 2*bufbuf
        else
            buf_num = bufbuf
        end if
        ! allocate buffers
        allocate (buffer(ipc*buf_dim, buf_num))
        buffer = 0.d0
        allocate (weight_buff(bufbuf))
        weight_buff = 1.d0
        allocated_size_buf = allocated_size_buf + ipc*buf_dim*buf_num
        if (contracted_on) then
            allocate (buffer_on(ipc*buf_dim_contr, buf_num))
            allocated_size_buf = allocated_size_buf + ipc*buf_dim_contr*buf_num
            buffer_on = 0.d0
        end if
        if (allocated(buf_conj)) deallocate (buf_conj)
        allocate (buf_conj(ipc*buf_dim_conj, buf_num))
        buf_conj = 0.d0
        allocated_size_buf = allocated_size_buf + ipc*buf_dim_conj*buf_num
        ! allocate scratch vectors
        allocate (wf_threading_scratch(ipc*wf_dim, thread_active))
        allocated_size_thread = allocated_size_thread + ipc*wf_dim*thread_active
        wf_threading_scratch = 0.d0
        if (double_overs) then
            if (allocated(wf_threading_scratch_down)) deallocate (wf_threading_scratch_down)
            allocate (wf_threading_scratch_down(ipc*wf_dim, thread_active))
            allocated_size_thread = allocated_size_thread + ipc*wf_dim*thread_active
            wf_threading_scratch_down = 0.d0
        end if
        if (present(print_size)) then
            if (print_size .and. rank .eq. 0) then
                if (.not. compute_bands) then
                    write (6, '(E12.5, A)') 8.d0*allocated_size_buf/1d9 &
                        , " Gbyte per MPI task for the buffer, proportional to nbufd!"
                    write (6, '(E12.5, A)') 8.d0*allocated_size_thread/1d9 &
                        , " Gbyte per MPI task for the buffer, proportional to threads!"
                end if
            end if
        end if
        return

    end subroutine allocate_buffers
    !
    ! fill the support array wf(:,indmesh) with the wavefunction
    ! values evaluated at one grid point
    ! If the thread_id (tid) argument is present, fill only one grid
    ! point, otherwise fill the full buffer.
    !
    subroutine update_wf_memlarge(nbas, & ! true basis set dimension
                                  wf, & ! array containing the wave function for memlarge option
                                  wf_dim, & ! leading dimension of wf
                                  bufbuf, & ! number of buffered points in the grid
                                  indmesh, & ! last index on the grid
                                  tid) ! thread identification

        implicit none
        !
        ! input
        integer, intent(in) :: nbas, indmesh, bufbuf, wf_dim
        real(8), intent(inout) :: wf(wf_dim, *)
        integer, intent(in), optional :: tid
        !
        ! local
        integer :: i

        !      write(6,*) 'debugging: ',tid,indmesh,wf_dim,nbas,size(wf_threading_scratch,1)

        if (ipc .eq. 1) then
            if (present(tid)) then ! one point
                wf(1:nbas, indmesh) = wf_threading_scratch(1:nbas, tid)
            else ! full buffer
                do i = 1, bufbuf
                    wf(1:nbas, indmesh - bufbuf + i) = buffer(1:nbas, i)
                end do
            end if
        else
            if (present(tid)) then ! one point
                wf(1:2*nbas, indmesh) = wf_threading_scratch(1:2*nbas, tid)
                if (double_overs) then
                    wf(2*nbas + 1:4*nbas, indmesh) = wf_threading_scratch_down(1:2*nbas, tid)
                end if
            else ! full buffer
                do i = 1, bufbuf
                    wf(1:2*nbas, indmesh - bufbuf + i) = buffer(1:2*nbas, i)
                    if (double_overs) then
                        wf(2*nbas + 1:4*nbas, indmesh - bufbuf + i) = buffer(1:2*nbas, bufbuf + i)
                    end if
                end do
            end if
        end if

        return
    end subroutine update_wf_memlarge

    subroutine deallocate_buffers()

        implicit none

        deallocate (buffer)
        if (allocated(wf_threading_scratch)) deallocate (wf_threading_scratch)
        if (contracted_on) then
            deallocate (buffer_on)
        end if
        if (allocated(buf_conj)) deallocate (buf_conj)
        if (allocated(wf_threading_scratch_down)) deallocate (wf_threading_scratch_down)
        deallocate (weight_buff)

        return

    end subroutine deallocate_buffers

end module buffers
