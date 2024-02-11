#ifdef _TREXIO
#include "trexio_f.f90"
#else

module trexio
 use, intrinsic :: iso_c_binding
 implicit none

 integer, parameter :: trexio_exit_code  = c_int32_t
 integer, parameter :: trexio_back_end_t = c_int32_t
 integer, parameter :: trexio_t          = c_size_t

 integer(trexio_exit_code), parameter :: TREXIO_SUCCESS = 0

end
#endif
