#ifdef _QMCKL
#include "qmckl_f.F90"

subroutine setup_qmckl_ctx(&
                &  number_of_atoms&
                &, number_of_shells&
                &, atomic_numbers&
                &, atomic_coordinates&
                &, shell_to_ion&
                &, shell_types&
                &, shell_exponents&
                &, context)

    use qmckl
    use constants, only: pi

    implicit none
    integer*4, intent(in)          :: number_of_atoms
    integer*4, intent(in)          :: number_of_shells
    integer(kind=qmckl_context), intent(out) :: context
    integer*4, intent(in)          :: shell_to_ion(number_of_shells)
    integer*4, intent(in)          :: shell_types(number_of_shells)

    real*8, intent(in)             :: shell_exponents(number_of_shells)
    real*8, intent(in)             :: atomic_numbers(number_of_atoms)
    real*8, intent(in)             :: atomic_coordinates(3,number_of_atoms)

    integer(kind=qmckl_exit_code)  :: rc

    integer*4                      :: ii, jj, kk, ll, mm, count

    integer*4                      :: shell_types_(number_of_shells)
    integer*8                      :: primitive_numbers_(number_of_shells)
    integer*8                      :: primitive_to_shell_(number_of_shells)
    integer*8                      :: ao_factor_number_

    real*8                         :: shell_factor_(number_of_shells)
    real*8                         :: shell_exponents_(number_of_shells)
    real*8                         :: shell_coefficients_(number_of_shells)
    real*8                         :: prim_factor_(number_of_shells)
    real*8, dimension(:), allocatable :: ao_factor_

    ! Number of shell for each atom
    integer*8                      :: nucleus_shell_num(number_of_atoms)
    integer*8                      :: nucleus_shell_index(number_of_atoms)

    !#####################################################################
    !#                                                                   #
    !#  This subroutine sets up the QMCkl context.                       #
    !#  Setup nuclei information.                                        #
    !#                                                                   #
    !#####################################################################

#ifdef _DEBUG
    write(*,*) "QMCKL: Creating context"
#endif

    context = qmckl_context_create()
    if (context.eq.QMCKL_NULL_CONTEXT) then
        write(0,*) "Error: qmckl_context_create"
        stop 1
    end if

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up nuclei number: ",I4)') number_of_atoms
#endif

    rc = qmckl_set_nucleus_num(context, 1_8 * number_of_atoms)
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_nucleus_num"
        stop 1
    end if

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up nuclei charge: ",F4.2)') atomic_numbers
#endif

    rc = qmckl_set_nucleus_charge(context, 1.0_8 * idnint(atomic_numbers), 1_8 * size(atomic_numbers))
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_nucleus_charge"
        stop 1
    end if

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up nuclei coord: ")')
    do ii = 1, number_of_atoms
        write(*,'(3F8.3)') atomic_coordinates(:,ii)
    end do
#endif

    rc = qmckl_set_nucleus_coord(context, "N", atomic_coordinates, 3_8 * number_of_atoms)
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_nucleus_coord"
        stop 1
    end if

    !#####################################################################
    !#                                                                   #
    !#  Setup the basis set.                                             #
    !#                                                                   #
    !#####################################################################

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up basis set type: G")')
#endif
    
    rc = qmckl_set_ao_basis_type(context, "G")
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_type"
        stop 1
    end if

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up basis set shell number: ",I4)') number_of_shells
#endif

    rc = qmckl_set_ao_basis_shell_num(context, 1_8 * number_of_shells)
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_shell_num"
        stop 1
    end if

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up basis set prim number:")')
    do ii = 1, number_of_shells
        write(*,'(I4)') shell_to_ion(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_prim_num(context, 1_8 * number_of_shells)
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_prim_num"
        stop 1
    end if

    nucleus_shell_num = 0_8
    do ii = 1, size(shell_to_ion)
        nucleus_shell_num(shell_to_ion(ii)) = nucleus_shell_num(shell_to_ion(ii)) + 1
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up basis set nucleus shell number: ")')
    do ii = 1, number_of_atoms
        write(*,'(I4)') nucleus_shell_num(ii)
    end do
#endif
    
    rc = qmckl_set_ao_basis_nucleus_shell_num(context, nucleus_shell_num, 1_8 * number_of_atoms)

    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_nucleus_shell_num"
        stop 1
    end if

    do ii = 1, number_of_atoms
        nucleus_shell_index(ii) = sum(nucleus_shell_num(1:ii-1))
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up basis set nucleus shell index: ")')
    do ii = 1, number_of_atoms
        write(*,'(I4)') nucleus_shell_index(ii)
    end do
#endif

    rc =  qmckl_set_ao_basis_nucleus_index(context, nucleus_shell_index, 1_8 * number_of_atoms)

    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_get_ao_basis_nucleus_index"
        stop 1
    end if

#ifdef _DEBUG
    print *, "QMCKL: Setting up basis set shell type (reading): "
    do ii = 1, number_of_shells
        write(*,'(I4)') shell_types(ii)
    end do
#endif

    do ii = 1, number_of_shells
        if (shell_types(ii).gt.99) then
            write(0,*) "Error: shell_indeces(ii) > 99"
            stop 1
        end if
        if (shell_types(ii).lt.90) then
            write(0,*) "Error: shell_indeces(ii) < 90"
            stop 1
        end if
        shell_types_(ii) = shell_types(ii) - 90
    end do

#ifdef _DEBUG
    print *, "QMCKL: Setting up basis set shell type: "
    do ii = 1, number_of_shells
        write(*,'(I4)') shell_types_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_shell_ang_mom(context, shell_types_, 1_8 * number_of_shells)

    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_shell_ang_mom"
        stop 1
    end if

    do ii = 1, number_of_shells
        primitive_numbers_(ii) = 1_8
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up primitive numbers for shells:")')
    do ii = 1, number_of_shells
        write(*,'(I4)') primitive_numbers_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_shell_prim_num(context&
                &, primitive_numbers_&
                &, 1_8 * number_of_shells)

    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_shell_prim_num"
        stop 1
    end if

    do ii = 1, number_of_shells
        primitive_to_shell_(ii) = ii - 1
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up primitive to shell map:")')
    do ii = 1, number_of_shells
        write(*,'(I4)') primitive_to_shell_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_shell_prim_index(context&
                &, primitive_to_shell_&
                &, 1_8 * number_of_shells)

    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_nucleus_shell_num"
        stop 1
    end if

    do ii = 1, number_of_shells
        shell_factor_(ii) = 1.0_8
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up shell factor:")')
    do ii = 1, number_of_shells
        write(*,'(F8.3)') shell_factor_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_shell_factor(context&
            &, shell_factor_&
            &, 1_8 * number_of_shells)

    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_shell_factor"
        stop 1
    end if

    do ii = 1, number_of_shells
        shell_exponents_(ii) = shell_exponents(ii)
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up shell exponents:")')
    do ii = 1, number_of_shells
        write(*,'(F8.3)') shell_exponents_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_exponent(context&
            &, shell_exponents_&
            &, 1_8 * number_of_shells)

    if (rc /= 0) then
        write(0,*) "qmckl_set_ao_basis_exponent failed"
        stop 1
    end if

    do ii = 1, number_of_shells
        shell_coefficients_(ii) = 1.0_8
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up shell coefficients:")')
    do ii = 1, number_of_shells
        write(*,'(F8.3)') shell_coefficients_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_coefficient(context&
            &, shell_coefficients_&
            &, 1_8 * number_of_shells)

    if (rc /= 0) then
        write(0,*) "qmckl_set_ao_basis_coefficient failed"
        stop 1
    end if

    do ii = 1, number_of_shells
        prim_factor_(ii) = (2.0_8 * shell_exponents_(ii) / pi)**(3.0_8/4.0_8)
        prim_factor_(ii) = prim_factor_(ii)&
                       & * (8.0_8 * shell_exponents_(ii))**(1.0_8 * shell_types_(ii) * 0.5_8)
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up primitive factor:")')
    do ii = 1, number_of_shells
        write(*,'(F8.3)') prim_factor_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_prim_factor(context&
            &, prim_factor_&
            &, 1_8 * number_of_shells)

    if (rc /= 0) then
        write(0,*) "qmckl_set_ao_basis_prim_factor failed"
        stop 1
    end if

    ao_factor_number_ = 0_8
    do ii = 1, number_of_shells
        ao_factor_number_ = ao_factor_number_&
                         &+ (shell_types_(ii) + 1_8) * (shell_types_(ii) + 2_8) / 2_8
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up ao factor number: ",I4)') ao_factor_number_
#endif

    rc = qmckl_set_ao_basis_ao_num(context&
            &, ao_factor_number_)

    if (rc /= 0) then
        write(0,*) "qmckl_set_ao_basis_ao_num failed"
        stop 1
    end if

    if (allocated(ao_factor_)) deallocate(ao_factor_)
    allocate(ao_factor_(ao_factor_number_))

    count = 1
    do ll = 1, number_of_shells
        do ii = shell_types_(ll), 0, -1
            do jj = shell_types_(ll) - ii, 0, -1
                kk = shell_types_(ll) - ii - jj
                ao_factor_(count) = 1.0_8
                do mm = ii + 1, 2 * ii
                   ao_factor_(count)  = ao_factor_(count) / mm
                end do
                do mm = jj + 1, 2 * jj
                   ao_factor_(count)  = ao_factor_(count) / mm
                end do
                do mm = kk + 1, 2 * kk
                   ao_factor_(count)  = ao_factor_(count) / mm
                end do
                ao_factor_(count) = dsqrt(ao_factor_(count))
                count = count + 1
            end do
        end do
    end do

#ifdef _DEBUG
    write(*,'("QMCKL: Setting up ao factor: ")')
    do ii = 1, ao_factor_number_
        write(*,'(F8.3)') ao_factor_(ii)
    end do
#endif

    rc = qmckl_set_ao_basis_ao_factor(context&
            &, ao_factor_&
            &, 1_8 * ao_factor_number_)

    if (rc /= 0) then
        write(0,*) "qmckl_set_ao_basis_ao_factor failed"
        stop 1
    end if

    if (allocated(ao_factor_)) deallocate(ao_factor_)

end subroutine setup_qmckl_ctx

#else
module qmckl
    ! This module defines the dummy Fortran interface to the QMCkl library.
    use, intrinsic                 :: iso_c_binding
    integer  , parameter           :: qmckl_context = c_int64_t
    integer*8, parameter           :: QMCKL_NULL_CONTEXT = 0
    integer  , parameter           :: qmckl_exit_code = c_int32_t
    integer(qmckl_exit_code), parameter :: QMCKL_SUCCESS = 0
    
end module
#endif
