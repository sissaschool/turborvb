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

!> @file buffer.f90
!> @brief Module providing buffer arrays for reduction operations in TurboRVB
!> @author TurboRVB group
!> @date 2022
!> @version 1.0
!> @details
!> This module defines buffer arrays used for MPI reduction operations and
!> temporary storage in TurboRVB calculations. The buffer size is configurable
!> based on compilation flags to optimize memory usage for different scenarios.
!>
!> The module provides:
!> - Integer buffer for index-based reductions
!> - Double precision buffer for real number reductions  
!> - Complex buffer for complex number reductions
!> - Configurable buffer size based on compilation flags
!>
!> @note Buffer size is determined by UNREL_DIAG compilation flag
!> @warning Buffer arrays are shared across all routines using this module

#if UNREL_DIAG
#define Buffersize 100000
#else
#define Buffersize 1000000
#endif

!> @brief Module containing buffer arrays for reduction operations
!> @details Provides pre-allocated buffer arrays for efficient MPI reduction
!> operations and temporary storage in quantum Monte Carlo calculations.
!> The buffer size is automatically determined based on compilation flags:
!> - UNREL_DIAG defined: 100,000 elements (smaller memory footprint)
!> - UNREL_DIAG undefined: 1,000,000 elements (larger capacity)
!>
!> These buffers are used throughout TurboRVB for:
!> - MPI collective operations (reduce, gather, scatter)
!> - Temporary storage during matrix operations
!> - Accumulation of partial results in parallel calculations
!> - Buffer space for communication in distributed memory computations
!>
!> @note All arrays are statically allocated for performance
!> @warning Buffer arrays are shared - avoid concurrent writes to same elements
module buffer
    implicit none
    
    !> @var bufdim Buffer dimension parameter
    !> @details The actual buffer size used in calculations.
    !> This parameter is set to the value of Buffersize macro.
    integer bufdim
    parameter(bufdim=Buffersize)
    
    !> @var buffi_reduce Integer buffer for reduction operations
    !> @details Array of integers used for index-based reductions and
    !> temporary storage of integer data in MPI operations.
    !> Size: Buffersize elements
    integer buffi_reduce(Buffersize)
    
    !> @var buff_reduce Double precision buffer for reduction operations
    !> @details Array of double precision real numbers used for
    !> accumulation and reduction of real-valued data in parallel
    !> quantum Monte Carlo calculations.
    !> Size: Buffersize elements
    real(8) buff_reduce(Buffersize)
    
    !> @var buff_reducec Complex buffer for reduction operations
    !> @details Array of double precision complex numbers used for
    !> accumulation and reduction of complex-valued data, particularly
    !> in calculations involving complex wave functions and operators.
    !> Size: Buffersize elements
    complex(8) buff_reducec(Buffersize)
    
end module buffer

#undef Buffersize
