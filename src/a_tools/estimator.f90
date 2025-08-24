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

!> @brief Statistical estimator module for TurboRVB
!> @details This module provides data structures and procedures for computing
!> statistical estimators including averages, variances, and correlation functions.
!> It supports both scalar and array data types with different levels of
!> statistical analysis.
!>
!> The module defines several derived types:
!> - avg_scalar/avg_array: For simple averaging
!> - esti_scalar/esti_array: For statistical analysis with variance
!> - corr_basic/corr_advanced: For correlation functions
!>
!> Each type has associated procedures for allocation, reset, data pushing,
!> calculation, and cleanup.
!>
!> @author TurboRVB group
!> @date 2022
module estimator
    implicit none

    !> @brief Scalar average estimator type
    !> @details Contains fields for number of samples, summation, and average
    type :: avg_scalar
        integer :: num        !< Number of samples
        double precision :: summation  !< Sum of all values
        double precision :: average    !< Computed average
    end type avg_scalar

    !> @brief Array average estimator type
    !> @details Contains fields for number of samples, dimension, and arrays
    !> for summation and average values
    type :: avg_array
        integer :: num        !< Number of samples
        integer :: dimen      !< Dimension of the array
        double precision, allocatable :: summation(:)  !< Sum of all values
        double precision, allocatable :: average(:)    !< Computed averages
    end type avg_array

    !> @brief Scalar statistical estimator type
    !> @details Contains fields for statistical analysis including variance
    !> and standard deviation
    type :: esti_scalar
        integer :: num        !< Number of samples
        double precision :: summation  !< Sum of all values
        double precision :: sumsqr     !< Sum of squared values
        double precision :: average    !< Computed average
        double precision :: vari       !< Computed variance
        double precision :: deviation  !< Computed standard deviation
    end type esti_scalar

    !> @brief Array statistical estimator type
    !> @details Contains fields for statistical analysis of array data
    type :: esti_array
        integer :: num        !< Number of samples
        integer :: dimen      !< Dimension of the array
        double precision, allocatable :: summation(:)  !< Sum of all values
        double precision, allocatable :: sumsqr(:)     !< Sum of squared values
        double precision, allocatable :: average(:)    !< Computed averages
        double precision, allocatable :: vari(:)       !< Computed variances
        double precision, allocatable :: deviation(:)  !< Computed standard deviations
    end type esti_array

    !> @brief Basic correlation function estimator type
    !> @details Uses average estimators for correlation analysis
    type :: corr_basic
        integer :: num        !< Number of samples
        type(avg_array) :: estia  !< Estimator for first variable
        type(avg_array) :: estib  !< Estimator for second variable
        double precision, allocatable :: corrsum(:)  !< Sum of correlations
        double precision, allocatable :: corrfun(:)  !< Computed correlation function
    end type corr_basic

    !> @brief Advanced correlation function estimator type
    !> @details Uses statistical estimators for correlation analysis
    type :: corr_advanced
        integer :: num        !< Number of samples
        type(esti_array) :: estia  !< Estimator for first variable
        type(esti_array) :: estib  !< Estimator for second variable
        double precision, allocatable :: corrsum(:)  !< Sum of correlations
        double precision, allocatable :: corrfun(:)  !< Computed correlation function
    end type corr_advanced

    !> @brief Generic interface for allocation procedures
    interface alloc
        module procedure allocate_esti_array
        module procedure allocate_avg_array
        module procedure allocate_corr_advanced_array
        module procedure allocate_corr_basic_array
    end interface alloc

    !> @brief Generic interface for deallocation procedures
    interface free
        module procedure free_esti_array
        module procedure free_avg_array
        module procedure free_corr_advanced_array
        module procedure free_corr_basic_array
    end interface free

    !> @brief Generic interface for reset procedures
    interface reset
        module procedure reset_esti_scalar, reset_esti_array
        module procedure reset_avg_scalar, reset_avg_array
        module procedure reset_corr_advanced_array
        module procedure reset_corr_basic_array
    end interface reset

    !> @brief Generic interface for data pushing procedures
    interface push
        module procedure push_esti_scalar, push_esti_array
        module procedure push_avg_scalar, push_avg_array
        module procedure push_corr_advanced_array
        module procedure push_corr_basic_array
    end interface push

    !> @brief Generic interface for calculation procedures
    interface calc
        module procedure calc_esti_scalar, calc_esti_array
        module procedure calc_avg_scalar, calc_avg_array
        module procedure calc_corr_advanced_array
        module procedure calc_corr_basic_array
    end interface calc

contains

    !!!avg!!!
    !> @brief Reset scalar average estimator
    !> @param estimator The average estimator to reset
    subroutine reset_avg_scalar(estimator)
        implicit none
        type(avg_scalar) :: estimator

        estimator%num = 0
        estimator%summation = 0.d0
        estimator%average = 0.d0
    end subroutine reset_avg_scalar

    !> @brief Reset array average estimator
    !> @param estimator The average estimator to reset
    subroutine reset_avg_array(estimator)
        implicit none
        type(avg_array) :: estimator

#ifdef _DEBUG
        if (.not. allocated(estimator%summation)) then
            write (*, *) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
            stop
        end if
#endif
        estimator%num = 0
        estimator%summation = 0.d0
        estimator%average = 0.d0
    end subroutine reset_avg_array

    !> @brief Allocate array average estimator
    !> @param estimator The average estimator to allocate
    !> @param ndim Dimension of the array
    subroutine allocate_avg_array(estimator, ndim)
        implicit none
        type(avg_array) :: estimator
        integer, intent(in) :: ndim

        if (.not. allocated(estimator%summation)) then
            estimator%dimen = ndim
            allocate (estimator%summation(ndim), estimator%average(ndim))
            estimator%num = 0
            estimator%summation = 0.d0
            estimator%average = 0.d0
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't allocate the same array twice!"
#endif
        end if
    end subroutine allocate_avg_array

    !> @brief Push data to scalar average estimator
    !> @param estimator The average estimator
    !> @param element The value to add
    subroutine push_avg_scalar(estimator, element)
        implicit none
        type(avg_scalar) :: estimator
        double precision :: element

        estimator%num = estimator%num + 1
        estimator%summation = estimator%summation + element
    end subroutine push_avg_scalar

    !> @brief Push data to array average estimator
    !> @param estimator The average estimator
    !> @param element The array of values to add
    subroutine push_avg_array(estimator, element)
        implicit none
        type(avg_array) :: estimator
        double precision, intent(in) :: element(:)

        estimator%num = estimator%num + 1
        estimator%summation = estimator%summation + element
    end subroutine push_avg_array

    !> @brief Calculate average for scalar estimator
    !> @param estimator The average estimator
    subroutine calc_avg_scalar(estimator)
        implicit none
        type(avg_scalar) :: estimator

        if (estimator%num <= 0) return
        estimator%average = estimator%summation/estimator%num
    end subroutine calc_avg_scalar

    !> @brief Calculate average for array estimator
    !> @param estimator The average estimator
    subroutine calc_avg_array(estimator)
        implicit none
        type(avg_array) :: estimator

        if (estimator%num <= 0) return
        estimator%average = estimator%summation/estimator%num
    end subroutine calc_avg_array

    !> @brief Free array average estimator memory
    !> @param estimator The average estimator to deallocate
    subroutine free_avg_array(estimator)
        implicit none
        type(avg_array) :: estimator

        if (allocated(estimator%summation)) then
            deallocate (estimator%summation, estimator%average)
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't free the same array twice!"
#endif
        end if
    end subroutine free_avg_array

    !!!esti!!!
    !> @brief Reset scalar statistical estimator
    !> @param estimator The statistical estimator to reset
    subroutine reset_esti_scalar(estimator)
        implicit none
        type(esti_scalar) :: estimator

        estimator%num = 0
        estimator%summation = 0.d0
        estimator%sumsqr = 0.d0
        estimator%average = 0.d0
        estimator%vari = 0.d0
        estimator%deviation = 0.d0

    end subroutine reset_esti_scalar

    !> @brief Allocate array statistical estimator
    !> @param estimator The statistical estimator to allocate
    !> @param ndim Dimension of the array
    subroutine allocate_esti_array(estimator, ndim)
        implicit none
        type(esti_array) :: estimator
        integer, intent(in) :: ndim
        if (.not. allocated(estimator%summation)) then
            estimator%dimen = ndim
            allocate (estimator%summation(ndim), estimator%average(ndim), &
                      estimator%sumsqr(ndim), estimator%vari(ndim), estimator%deviation(ndim))
            estimator%num = 0
            estimator%summation = 0.d0
            estimator%sumsqr = 0.d0
            estimator%average = 0.d0
            estimator%vari = 0.d0
            estimator%deviation = 0.d0
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't allocate the same array twice!"
#endif
        end if
    end subroutine allocate_esti_array

    !> @brief Reset array statistical estimator
    !> @param estimator The statistical estimator to reset
    subroutine reset_esti_array(estimator)
        implicit none
        type(esti_array) :: estimator

#ifdef _DEBUG
        if (.not. allocated(estimator%summation)) then
            write (*, *) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
            stop
        end if
#endif
        estimator%num = 0
        estimator%summation = 0.d0
        estimator%sumsqr = 0.d0
        estimator%average = 0.d0
        estimator%vari = 0.d0
        estimator%deviation = 0.d0

    end subroutine reset_esti_array

    !> @brief Push data to scalar statistical estimator
    !> @param estimator The statistical estimator
    !> @param element The value to add
    subroutine push_esti_scalar(estimator, element)
        implicit none
        type(esti_scalar) :: estimator
        double precision :: element

        estimator%num = estimator%num + 1
        estimator%summation = estimator%summation + element
        estimator%sumsqr = estimator%sumsqr + element*element
    end subroutine push_esti_scalar

    !> @brief Push data to array statistical estimator
    !> @param estimator The statistical estimator
    !> @param element The array of values to add
    subroutine push_esti_array(estimator, element)
        implicit none
        type(esti_array) :: estimator
        double precision, intent(in) :: element(:)

        estimator%num = estimator%num + 1
        estimator%summation = estimator%summation + element
        estimator%sumsqr = estimator%sumsqr + element*element
    end subroutine push_esti_array

    !> @brief Calculate statistical analysis for scalar estimator
    !> @param estimator The statistical estimator
    subroutine calc_esti_scalar(estimator)
        implicit none
        type(esti_scalar) :: estimator
        integer :: ii

        if (estimator%num <= 0) return
        estimator%average = estimator%summation/estimator%num
        estimator%vari = estimator%sumsqr/estimator%num - estimator%average*estimator%average
        estimator%deviation = dsqrt(estimator%vari)
        if (estimator%vari < 0.d0) write (6, *) "Warning : negative variation due to machine precision!"
    end subroutine calc_esti_scalar

    !> @brief Calculate statistical analysis for array estimator
    !> @param estimator The statistical estimator
    subroutine calc_esti_array(estimator)
        implicit none
        type(esti_array) :: estimator
        integer :: ii

        if (estimator%num <= 0) return
        estimator%average = estimator%summation/estimator%num
        estimator%vari = estimator%sumsqr/estimator%num - estimator%average*estimator%average
        estimator%deviation = dsqrt(estimator%vari)
        do ii = 1, estimator%dimen
            if (estimator%vari(ii) < 0.d0) write (6, *) "Warning : negative variation due to machine precision!"
        end do
    end subroutine calc_esti_array

    !> @brief Free array statistical estimator memory
    !> @param estimator The statistical estimator to deallocate
    subroutine free_esti_array(estimator)
        implicit none
        type(esti_array) :: estimator

        if (allocated(estimator%summation)) then
            deallocate (estimator%summation, estimator%average, estimator%sumsqr, estimator%vari, estimator%deviation)
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't free the same array twice!"
#endif
        end if
    end subroutine free_esti_array

    !!!corr_basic!!!
    !> @brief Allocate basic correlation function estimator
    !> @param corr The correlation estimator to allocate
    !> @param ndima Dimension of the first variable
    !> @param ndimb Dimension of the second variable
    subroutine allocate_corr_basic_array(corr, ndima, ndimb)
        implicit none
        type(corr_basic) :: corr
        integer, intent(in) :: ndima, ndimb

        if (.not. allocated(corr%corrsum)) then
            call alloc(corr%estia, ndima)
            call alloc(corr%estib, ndimb)
            allocate (corr%corrsum(ndima*ndimb), corr%corrfun(ndima*ndimb))
            corr%num = 0
            corr%corrsum = 0.d0
            corr%corrfun = 0.d0
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't allocate the same array twice!"
#endif
        end if
    end subroutine allocate_corr_basic_array

    !> @brief Reset basic correlation function estimator
    !> @param corr The correlation estimator to reset
    subroutine reset_corr_basic_array(corr)
        implicit none
        type(corr_basic) :: corr

#ifdef _DEBUG
        if (.not. allocated(corr%corrsum)) then
            write (*, *) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
            stop
        end if
#endif
        call reset(corr%estia)
        call reset(corr%estib)
        corr%num = 0
        corr%corrsum = 0.d0
        corr%corrfun = 0.d0
    end subroutine reset_corr_basic_array

    !> @brief Push data to basic correlation function estimator
    !> @param corr The correlation estimator
    !> @param elementa The first variable array
    !> @param elementb The second variable array
    subroutine push_corr_basic_array(corr, elementa, elementb)
        implicit none
        type(corr_basic) :: corr
        double precision, intent(in) :: elementa(:), elementb(:)
        integer :: ii, jj

        corr%num = corr%num + 1
        call push(corr%estia, elementa)
        call push(corr%estib, elementb)
        call dgemm('N', 'N', corr%estia%dimen, corr%estib%dimen, 1, 1.d0, &
                   elementa, corr%estia%dimen, elementb, 1, 1.d0, corr%corrsum, corr%estia%dimen)
    end subroutine push_corr_basic_array

    !> @brief Calculate correlation function for basic correlation function estimator
    !> @param corr The correlation estimator
    subroutine calc_corr_basic_array(corr)
        implicit none
        type(corr_basic) :: corr
        integer :: ii, jj

        if (corr%num <= 0) return
        call calc(corr%estia)
        call calc(corr%estib)
        corr%corrfun = corr%corrsum/corr%num
        call dgemm('N', 'N', corr%estia%dimen, corr%estib%dimen, 1, -1.d0, &
                   corr%estia%average, corr%estia%dimen, corr%estib%average, 1, 1.d0, corr%corrfun, corr%estia%dimen)
    end subroutine calc_corr_basic_array

    !> @brief Free basic correlation function estimator memory
    !> @param corr The correlation estimator to deallocate
    subroutine free_corr_basic_array(corr)
        implicit none
        type(corr_basic) :: corr

        if (allocated(corr%corrsum)) then
            call free_avg_array(corr%estia)
            call free_avg_array(corr%estib)
            deallocate (corr%corrsum, corr%corrfun)
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't free the same array twice!"
#endif
        end if
    end subroutine free_corr_basic_array

    !!!corr_advanced!!!
    !> @brief Allocate advanced correlation function estimator
    !> @param corr The correlation estimator to allocate
    !> @param ndima Dimension of the first variable
    !> @param ndimb Dimension of the second variable
    subroutine allocate_corr_advanced_array(corr, ndima, ndimb)
        implicit none
        type(corr_advanced) :: corr
        integer, intent(in) :: ndima, ndimb

        if (.not. allocated(corr%corrsum)) then
            call alloc(corr%estia, ndima)
            call alloc(corr%estib, ndimb)
            allocate (corr%corrsum(ndima*ndimb), corr%corrfun(ndima*ndimb))
            corr%num = 0
            corr%corrsum = 0.d0
            corr%corrfun = 0.d0
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't allocate the same array twice!"
#endif
        end if
    end subroutine allocate_corr_advanced_array

    !> @brief Reset advanced correlation function estimator
    !> @param corr The correlation estimator to reset
    subroutine reset_corr_advanced_array(corr)
        implicit none
        type(corr_advanced) :: corr

#ifdef _DEBUG
        if (.not. allocated(corr%corrsum)) then
            write (*, *) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
            stop
        end if
#endif
        call reset(corr%estia)
        call reset(corr%estib)
        corr%num = 0
        corr%corrsum = 0.d0
        corr%corrfun = 0.d0
    end subroutine reset_corr_advanced_array

    !> @brief Push data to advanced correlation function estimator
    !> @param corr The correlation estimator
    !> @param elementa The first variable array
    !> @param elementb The second variable array
    subroutine push_corr_advanced_array(corr, elementa, elementb)
        implicit none
        type(corr_advanced) :: corr
        double precision, intent(in) :: elementa(:), elementb(:)
        integer :: ii, jj

        corr%num = corr%num + 1
        call push(corr%estia, elementa)
        call push(corr%estib, elementb)
        call dgemm('N', 'N', corr%estia%dimen, corr%estib%dimen, 1, 1.d0, &
                   elementa, corr%estia%dimen, elementb, 1, 1.d0, corr%corrsum, corr%estia%dimen)
    end subroutine push_corr_advanced_array

    !> @brief Calculate correlation function for advanced correlation function estimator
    !> @param corr The correlation estimator
    subroutine calc_corr_advanced_array(corr)
        implicit none
        type(corr_advanced) :: corr
        integer :: ii, jj

        if (corr%num <= 0) return
        call calc(corr%estia)
        call calc(corr%estib)
        call dgemm('N', 'N', corr%estia%dimen, corr%estib%dimen, 1, -1.d0, &
                   corr%estia%average, corr%estia%dimen, corr%estib%average, 1, 1.d0, corr%corrfun, corr%estia%dimen)
    end subroutine calc_corr_advanced_array

    !> @brief Free advanced correlation function estimator memory
    !> @param corr The correlation estimator to deallocate
    subroutine free_corr_advanced_array(corr)
        implicit none
        type(corr_advanced) :: corr

        if (allocated(corr%corrsum)) then
            call free_esti_array(corr%estia)
            call free_esti_array(corr%estib)
            deallocate (corr%corrsum, corr%corrfun)
#ifdef _DEBUG
        else
            write (*, *) "Warning! you can't free the same array twice!"
#endif
        end if
    end subroutine free_corr_advanced_array

end module estimator
