#ifdef _QMCKL
#include "qmckl_f.F90"
#else
module qmckl
  ! This module defines the dummy Fortran interface to the QMCkl library.
  use, intrinsic :: iso_c_binding
  integer  , parameter :: qmckl_context = c_int64_t
  integer*8, parameter :: QMCKL_NULL_CONTEXT = 0
  integer  , parameter :: qmckl_exit_code = c_int32_t
  integer(qmckl_exit_code), parameter :: QMCKL_SUCCESS = 0
end module
#endif
