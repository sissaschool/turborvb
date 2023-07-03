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

program offload_pointer_transcription
    implicit none

    real :: A(10)
    integer :: b = 0

    A = 0

!#ifndef _OFFLOAD
!#error OFFLOAD is not supported, turn on _OFFLOAD flag
!#endif
!$omp target data map(alloc:A)
#ifdef RISC
    call transcribe_(A, b)
#else
    call transcribe(A, b)
#endif
!$omp end target data

    if (b .eq. 0) then
        stop 0
    else
        stop 1
    end if

end program offload_pointer_transcription

