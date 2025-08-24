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

!> @brief Set specific elements of a vector to a constant value
!>
!> This subroutine sets every m-th element of a vector to a specified value.
!> It is used for initializing arrays with a specific pattern, particularly
!> useful in quantum Monte Carlo calculations for setting up initial conditions
!> or clearing specific array elements.
!>
!> @param[in] n Number of elements to set
!> @param[in] zero Value to assign to the selected elements
!> @param[in,out] vet Vector to be modified (m*(n-1)+1)
!> @param[in] m Stride between elements to be set
!>
!> @details
!> The subroutine sets vet(m*(i-1)+1) = zero for i = 1 to n.
!> This effectively sets every m-th element starting from the first element.
!> For m=1, all n elements are set to zero.
!> For m=2, elements 1, 3, 5, ... are set to zero.
!>
!> @note This is a serial version without OpenMP parallelization
subroutine dscalzero(n, zero, vet, m)
    implicit none
    integer n, m, i
    real*8 zero, vet(m*(n - 1) + 1)
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
    return
end

!> @brief Set specific elements of a vector to a constant value (OpenMP version)
!>
!> This subroutine sets every m-th element of a vector to a specified value
!> using OpenMP parallelization for improved performance on multi-core systems.
!>
!> @param[in] n Number of elements to set
!> @param[in] zero Value to assign to the selected elements
!> @param[in,out] vet Vector to be modified (m*(n-1)+1)
!> @param[in] m Stride between elements to be set
!>
!> @details
!> The subroutine uses OpenMP target teams distribute parallel do directive
!> when _OFFLOAD is defined, enabling parallel execution on accelerators.
!> Sets vet(m*(i-1)+1) = zero for i = 1 to n in parallel.
!>
!> @note This version includes OpenMP offloading support
subroutine dscalzero_(n, zero, vet, m)
    implicit none
    integer n, m, i
    real*8 zero, vet(m*(n - 1) + 1)
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end

!> @brief Set specific elements of a vector to a constant value (conditional OpenMP version)
!>
!> This subroutine sets every m-th element of a vector to a specified value
!> using conditional OpenMP parallelization based on the yes_ontarget flag.
!>
!> @param[in] n Number of elements to set
!> @param[in] zero Value to assign to the selected elements
!> @param[in,out] vet Vector to be modified (m*(n-1)+1)
!> @param[in] m Stride between elements to be set
!>
!> @details
!> The subroutine uses conditional OpenMP target teams distribute parallel do
!> directive when _OFFLOAD is defined and yes_ontarget is true.
!> This allows runtime control over whether to use accelerator offloading.
!> Sets vet(m*(i-1)+1) = zero for i = 1 to n.
!>
!> @note This version includes conditional OpenMP offloading support
subroutine dscalzero__(n, zero, vet, m)
    use constants, only: yes_ontarget
    implicit none
    integer n, m, i
    real*8 zero, vet(m*(n - 1) + 1)
!    if(n.le.16384.and.m.eq.1) then
!    do i = 1, n
!    vet(m*(i-1)+1) = zero
!    enddo
!#ifdef  _OFFLOAD
!!$omp target update to (vet)  if(yes_ontarget)
!#endif
!    else
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end

!> @brief Set specific elements of an integer vector to a constant value (OpenMP version)
!>
!> This subroutine sets every m-th element of an integer vector to a specified value
!> using OpenMP parallelization for improved performance on multi-core systems.
!>
!> @param[in] n Number of elements to set
!> @param[in] zero Integer value to assign to the selected elements
!> @param[in,out] vet Integer vector to be modified (m*(n-1)+1)
!> @param[in] m Stride between elements to be set
!>
!> @details
!> The subroutine uses OpenMP target teams distribute parallel do directive
!> when _OFFLOAD is defined, enabling parallel execution on accelerators.
!> Sets vet(m*(i-1)+1) = zero for i = 1 to n in parallel.
!>
!> @note This version handles integer arrays with OpenMP offloading support
subroutine iscalzero_(n, zero, vet, m)
    implicit none
    integer n, m, i
    integer zero, vet(m*(n - 1) + 1)
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end

!> @brief Set specific elements of a complex vector to a constant value (OpenMP version)
!>
!> This subroutine sets every m-th element of a complex vector to a specified value
!> using OpenMP parallelization for improved performance on multi-core systems.
!>
!> @param[in] n Number of elements to set
!> @param[in] zero Complex value to assign to the selected elements
!> @param[in,out] vet Complex vector to be modified (m*(n-1)+1)
!> @param[in] m Stride between elements to be set
!>
!> @details
!> The subroutine uses OpenMP target teams distribute parallel do directive
!> when _OFFLOAD is defined, enabling parallel execution on accelerators.
!> Sets vet(m*(i-1)+1) = zero for i = 1 to n in parallel.
!>
!> @note This version handles complex arrays with OpenMP offloading support
subroutine zscalzero_(n, zero, vet, m)
    implicit none
    integer n, m, i
    complex*16 zero, vet(m*(n - 1) + 1)
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end
