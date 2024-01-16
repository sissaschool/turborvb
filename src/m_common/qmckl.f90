#ifdef _QMCKL
#include "qmckl_f.F90"

subroutine setup_qmckl_ctx(                                              &
           &atomic_numbers                                               &
           &, atomic_coordinates                                         &
           &, number_of_atoms                                            &
           &, shell_indeces                                              &
           &, shell_number                                               &
           &, context)

    use qmckl

    implicit none
    integer*4, intent(in)          :: number_of_atoms, shell_number
    integer(kind=qmckl_context), intent(out) :: context
    integer*4, intent(in)          :: shell_indeces(shell_number)

    real*8, intent(in)             :: atomic_numbers(number_of_atoms)
    real*8, intent(in)             :: atomic_coordinates(3,number_of_atoms)

    integer(kind=qmckl_exit_code)  :: rc

    integer*4                      :: ii

    integer*8                      :: nucleus_indeces_(number_of_atoms)

    !#####################################################################
    !#                                                                   #
    !#  This subroutine sets up the QMCkl context.                       #
    !#  Setup nuclei information.                                        #
    !#                                                                   #
    !#####################################################################

    context = qmckl_context_create()
    if (context.eq.QMCKL_NULL_CONTEXT) then
        write(0,*) "Error: qmckl_context_create"
        stop 1
    end if

    rc = qmckl_set_nucleus_num(context, 1_8 * number_of_atoms)
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_nucleus_num"
        stop 1
    end if

    rc = qmckl_set_nucleus_charge(context, 1.0_8 * idnint(atomic_numbers), 1_8 * size(atomic_numbers))
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_nucleus_charge"
        stop 1
    end if

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

    
    rc = qmckl_set_ao_basis_type(context, "G")
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_type"
        stop 1
    end if

    rc = qmckl_set_ao_basis_shell_num(context, 1_8 * shell_number)
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_shell_num"
        stop 1
    end if

    rc = qmckl_set_ao_basis_prim_num(context, 1_8 * shell_number)
    if (rc.ne.QMCKL_SUCCESS) then
        write(0,*) "Error: qmckl_set_ao_basis_prim_num"
        stop 1
    end if

    do ii = 1, shell_number
        if (shell_indeces(ii).gt.99) then
            write(0,*) "Error: shell_indeces(ii) > 99"
            stop 1
        end if
        if (shell_indeces(ii).lt.90) then
            write(0,*) "Error: shell_indeces(ii) < 90"
            stop 1
        end if
        nucleus_indeces_(ii) = shell_indeces(ii) - 90
    end do

    !rc = qmckl_set_ao_basis_nucleus_shell_num(context, nucleus_shell_num, 1_8 * shell_number)

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
