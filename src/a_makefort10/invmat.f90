! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2004 PWSCF-CP-FPMD group
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

!
#include "f_defs.h"
!
subroutine invmat(n, a, a_inv, da)
    !-----------------------------------------------------------------------
    ! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
    ! if the matrix is dimensioned 3x3, it also computes determinant "da"
    ! matrix "a" is unchanged on output - LAPACK
    !
    use kinds
    implicit none
    integer :: n
    real(DP), dimension(n, n) :: a, a_inv
    real(DP) :: da
    !
    integer :: info, lda, lwork, ipiv(n)
    ! info=0: inversion was successful
    ! lda   : leading dimension (the same as n)
    ! ipiv  : work space for pivoting (assumed of length lwork=n)
    real(DP) :: work(n)
    ! more work space
    !
    lda = n
    lwork = n
    !
    a_inv(:, :) = a(:, :)
    !
    call DGETRF(n, n, a_inv, lda, ipiv, info)
    if (info .ne. 0) call errore('invmat', 'error in DGETRF', abs(info))
    call DGETRI(n, a_inv, lda, ipiv, work, lwork, info)
    if (info .ne. 0) call errore('invmat', 'error in DGETRI', abs(info))
    !
    if (n == 3) then
        da = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) + &
             a(1, 2)*(a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)) + &
             a(1, 3)*(a(2, 1)*a(3, 2) - a(3, 1)*a(2, 2))
        if (abs(da) < 1.d-10) call errore(' invmat ', ' singular matrix ', 1)
    else
        da = 0.d0
    end if

    return
end subroutine invmat
