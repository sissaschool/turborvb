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

module handle
#ifdef _CUSOLVE
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), target :: cusolver_handle
    integer :: cusolver_error
    interface

        integer(c_int) function cusolver_handle_init(handle) bind(c, name="cusolverDnCreate")
            use, intrinsic :: iso_c_binding
            type(c_ptr) :: handle
        end function cusolver_handle_init

        integer(c_int) function cusolver_handle_destroy(handle) bind(c, name="cusolverDnDestroy")
            use, intrinsic :: iso_c_binding
            type(c_ptr), value :: handle
        end function cusolver_handle_destroy

    end interface

contains

    subroutine init_handle
        implicit none
        cusolver_error = cusolver_handle_init(cusolver_handle)
    end subroutine init_handle

    subroutine finalize_handle
        implicit none
        cusolver_error = cusolver_handle_destroy(cusolver_handle)
    end subroutine finalize_handle
#endif
end module handle
