! MIT License
!
! Copyright (c) 2019-2021 stdlib contributors
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

! Fortran-lang/stdlib v0.8.1

!> Version: experimental
!>
!> The specification of this module is available [here](../page/specs/stdlib_kinds.html).
module stdlib_kinds
  use iso_fortran_env, only: int8, int16, int32, int64
  use iso_c_binding, only: c_bool, c_char
  implicit none
  private
  public :: sp, dp, xdp, qp, int8, int16, int32, int64, lk, c_bool, c_char

  !> Single precision real numbers
  integer, parameter :: sp = selected_real_kind(6)

  !> Double precision real numbers
  integer, parameter :: dp = selected_real_kind(15)

  !> Extended double precision real numbers
  integer, parameter :: xdp = -1

  !> Quadruple precision real numbers
  integer, parameter :: qp = -1

  !> Default logical kind parameter
  integer, parameter :: lk = kind(.true.)

end module stdlib_kinds
