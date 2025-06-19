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

/**
 * @file molec_pfaff.f90
 * @brief Pfaffian molecular orbital calculations for quantum Monte Carlo
 *
 * This file contains subroutines for calculating Pfaffian molecular orbitals
 * in TurboRVB. Pfaffians are used to represent antisymmetric wave functions
 * for systems with odd numbers of electrons or for certain types of
 * correlated wave functions.
 *
 * @details
 * The main algorithm involves:
 * 1. Tridiagonalization of the skew-symmetric matrix using Pfaffian decomposition
 * 2. Transformation to real symmetric form for eigenvalue calculation
 * 3. Diagonalization to obtain eigenvalues and eigenvectors
 * 4. Back-transformation to obtain the final Pfaffian molecular orbitals
 *
 * Key subroutines:
 * - pfaffian_mo: Main driver for Pfaffian molecular orbital calculation
 * - pfatriag: Tridiagonalization using Pfaffian decomposition
 * - symmtriang: Transformation to real symmetric form
 * - finalize_mopfaff: Back-transformation to final orbitals
 *
 * @note
 * - Supports both real and complex matrices (ipc=1,2)
 * - Uses LAPACK routines for eigenvalue calculations
 * - Includes debug options for verification
 * - Handles both even and odd numbers of orbitals
 *
 * @author TurboRVB group
 * @date 2022
 */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutines for the calculation of the pfaffians
!molecular orbitals            C.G.
!
!Output: outvl (nelorb_c/2) eigenvalues sorted by the
!        value (only the real part)
!Output: outvct(ipc*nelorb_c,nelorb_c) eigenvectors
!        ordered as egvl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Calculate Pfaffian molecular orbitals from a skew-symmetric matrix
 *
 * This subroutine computes the Pfaffian molecular orbitals by diagonalizing
 * a skew-symmetric matrix. The algorithm involves tridiagonalization,
 * transformation to real symmetric form, eigenvalue calculation, and
 * back-transformation to obtain the final orbitals.
 *
 * @param[in] lda Leading dimension of the input matrix
 * @param[in] nelorb_c Number of orbitals
 * @param[in] ipc Complex flag (1=real, 2=complex)
 * @param[in] detmat_c Input skew-symmetric matrix
 * @param[out] outvl Eigenvalues of the Pfaffian (nelorb_c/2 values)
 * @param[out] outvct Eigenvectors/molecular orbitals (ipc*lda x nelorb_c)
 *
 * @details
 * The algorithm proceeds as follows:
 * 1. Copy input matrix to U1 for tridiagonalization
 * 2. Apply Pfaffian tridiagonalization (pfatriag)
 * 3. Transform to real symmetric form (symmtriang)
 * 4. Diagonalize using DSTEVX (LAPACK)
 * 5. Apply gauge fixing to eigenvectors
 * 6. Back-transform to final orbitals (finalize_mopfaff)
 *
 * The eigenvalues are sorted by magnitude and only the positive
 * eigenvalues are returned (skew-symmetric matrices have paired
 * eigenvalues ±λ).
 *
 * @note
 * - Uses LAPACK DSTEVX for eigenvalue calculation
 * - Includes gauge fixing for consistent eigenvector signs
 * - Supports both real and complex matrices
 * - Memory is allocated and deallocated within the subroutine
 *
 * @see pfatriag(), symmtriang(), finalize_mopfaff(), orb_max()
 *
 * @author C.G.
 * @date 2022
 */
subroutine pfaffian_mo(lda, nelorb_c, ipc, detmat_c, outvl, outvct)
    implicit none
    integer :: nelorb_c, ipc, lda
    integer :: i, j, ierr ! auxiliary variables
    real(8) :: detmat_c(ipc*lda, nelorb_c), outvl(nelorb_c), outvct(ipc*lda, nelorb_c)
    real(8) :: auxr, orbmax
    real(8), external :: dlamch, orb_max
    !U2 matrix is no more necessary (ZSKTRD gives a real matrix as output)
    !U1, U2, and U3 are the operations required to convert the
    !detmat_c in the real symmetric  matrix whose elements are
    !stored in lambdars (the first nelorb_c values are 0.d0 and
    !correspond to the diagonal elements, the second nelorb_c
    !values are the off-diagonal ones)(it has
    !to be used in the DSTEVX routine)
    !vector with information necessary for the DSTEVX routine
    !work, iwork, ifail are vectors for the dstevx routine,
    !detmattr is the tridiagonal matrix needed during the calculations
    real(8), allocatable :: U1(:, :), U3(:, :), work(:), lambdars(:), eigvalues(:), &
                            eigvect(:, :), detmattr(:, :), auxmat(:, :), auxmat2(:, :), auxmat1(:, :)
    integer, allocatable :: iwork(:), ifail(:)
    complex(8) :: zzero, zone

    zzero = (0.d0, 0.d0)
    zone = (1.d0, 0.d0)

    allocate (U3(2*nelorb_c, nelorb_c), lambdars(2*nelorb_c - 1), eigvalues(nelorb_c), &
              eigvect(nelorb_c, nelorb_c), work(5*nelorb_c), iwork(5*nelorb_c), ifail(nelorb_c), &
              detmattr(ipc*nelorb_c, nelorb_c))
    allocate (U1(ipc*nelorb_c, nelorb_c))
    U1(1:ipc*nelorb_c, 1:nelorb_c) = detmat_c(1:ipc*nelorb_c, 1:nelorb_c)

    call pfatriag(nelorb_c, ipc, detmattr, U1)

    !The U3 matrix and lambdars are calculated
    call symmtriang(nelorb_c, ipc, detmattr, lambdars, U3)

    !Double precision symmetric triangular matrix diagonalization
    call dstevx("V", "A", nelorb_c, lambdars, lambdars(nelorb_c + 1), auxr, auxr, j, j, &
                2*dlamch('S'), i, eigvalues, eigvect, nelorb_c, work, iwork, ifail, ierr)
    !  choose a gauge
    do i = 1, nelorb_c
        orbmax = orb_max(nelorb_c, eigvect(1, i))
        if (orbmax .lt. 0.d0) eigvect(:, i) = -eigvect(:, i)
    end do

    if (ierr .gt. 0) then
        write (6, *) "ERROR DSTEVX: the eigenvectors:", ifail(1:ierr), "did not converged!"
    else if (ierr .lt. 0) then
        write (6, *) "ERROR DSTEVX: the parameter:", ierr, "has an illegal value!"
    end if

    call finalize_mopfaff(lda, nelorb_c, ipc, U1, U3, eigvalues, eigvect, outvl, outvct, detmat_c)

    deallocate (U1, U3, lambdars, eigvalues, eigvect, work, iwork, ifail, detmattr)
end subroutine pfaffian_mo

/**
 * @brief Find the maximum element in a vector for gauge fixing
 *
 * This function finds the first element in a vector that exceeds a threshold
 * value, used for gauge fixing of eigenvectors in Pfaffian calculations.
 * The threshold is set to avoid numerical issues with very small elements.
 *
 * @param[in] n Size of the vector
 * @param[in] vect Input vector of real values
 * @return Real value of the first element above threshold, or error if none found
 *
 * @details
 * The function searches for the first element in vect that satisfies
 * vect(i)^2 > safemin, where safemin = 0.5773576451/n.
 * This threshold is chosen to avoid numerical issues with poorly
 * normalized eigenvectors.
 *
 * @note
 * - Used for gauge fixing in pfaffian_mo
 * - Returns error message if no suitable element is found
 * - Threshold depends on vector size n
 * - Assumes vector should be properly normalized
 *
 * @see pfaffian_mo()
 *
 * @author C.G.
 * @date 2022
 */
function orb_max(n, vect)
    integer n, i
    real(8) orb_max, safemin
    real(8) vect(n)
    safemin = 0.5773576451d0/n ! To avoid particular cases.
    !  The gauge is defined by  using the first element > threshold in  the chosen space
    i = 1
    do while (i .le. n .and. vect(i)**2 .le. safemin)
        i = i + 1
    end do
    if (i .le. n) then
        orb_max = vect(i)
    else
        !  if all the elements satisfies the dowhile ineq. the normalization of vect < 0.57.. not possible
        write (6, *) ' ERROR check normalization dstevx in molec_pfaff '
    end if
    return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! It cretes the operator U_3 and if #debug also the
! checks the calculation maps the hermitian
! matrix \lambda_{iH} in \lambda_R
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Transform skew-symmetric matrix to real symmetric form
 *
 * This subroutine creates the transformation operator U3 and transforms
 * the skew-symmetric tridiagonal matrix to real symmetric form for
 * eigenvalue calculation. The transformation maps the Hermitian matrix
 * λ_{iH} to the real symmetric matrix λ_R.
 *
 * @param[in] nelorb_c Number of orbitals
 * @param[in] ipc Complex flag (1=real, 2=complex)
 * @param[in] detmattr Tridiagonal skew-symmetric matrix
 * @param[out] lambdars Real symmetric matrix in packed storage format
 * @param[out] U3 Transformation matrix for the conversion
 *
 * @details
 * The subroutine performs the following operations:
 * 1. Constructs the U3 transformation matrix with alternating signs
 * 2. Extracts off-diagonal elements from detmattr to form lambdars
 * 3. For real matrices (ipc=1): copies upper diagonal elements
 * 4. For complex matrices (ipc=2): copies real parts of upper diagonal
 * 5. In debug mode: verifies the transformation is correct
 *
 * The U3 matrix has a specific pattern:
 * - Elements (2i-1,i) = 1 for i mod 4 = 1
 * - Elements (2i,i) = 1 for i mod 4 = 2  
 * - Elements (2i-1,i) = -1 for i mod 4 = 3
 * - Elements (2i,i) = -1 for i mod 4 = 0
 *
 * @note
 * - lambdars is stored in packed format for LAPACK DSTEVX
 * - Debug mode includes verification of the transformation
 * - Supports both real and complex input matrices
 * - U3 matrix is used in the final back-transformation
 *
 * @see pfatriag(), finalize_mopfaff()
 *
 * @author C.G.
 * @date 2022
 */
subroutine symmtriang(nelorb_c, ipc, detmattr, lambdars, U3)
    implicit none
    integer :: nelorb_c, ipc
    integer :: i, j ! auxiliary variables
    real(8) :: detmattr(ipc*nelorb_c, nelorb_c), lambdars(2*nelorb_c - 1)
    real(8) :: U3(2*nelorb_c, nelorb_c)
    !the lambda*_test are matrices for the debug version of the code
    real(8), allocatable :: lambdars_test(:, :), lambdaih_test(:, :), auxmat(:, :)
    complex(8) :: zzero, zone

    zzero = (0.d0, 0.d0)
    zone = (1.d0, 0.d0)
    !building U3
    U3 = 0.d0
    do i = 1, nelorb_c
        if (mod(i, 4) .eq. 1) then
            U3(2*i - 1, i) = 1.d0
        else if (mod(i, 4) .eq. 2) then
            U3(2*i, i) = 1.d0
        else if (mod(i, 4) .eq. 3) then
            U3(2*i - 1, i) = -1.d0
        else
            U3(2*i, i) = -1.d0
        end if
    end do

    lambdars = 0.d0
    !building lambdars ipc.eq.1
    if (ipc .eq. 1) then
        do i = 1, nelorb_c - 1
            lambdars(nelorb_c + i) = detmattr(i, i + 1)
        end do
    else !building lambdars ipc.eq.2
        do i = 1, nelorb_c - 1
            lambdars(nelorb_c + i) = detmattr(2*i - 1, i + 1)
        end do
    end if

#ifdef DEBUG
    allocate (lambdars_test(2*nelorb_c, nelorb_c), lambdaih_test(2*nelorb_c, nelorb_c), &
              auxmat(2*nelorb_c, nelorb_c))
    lambdars_test = 0.d0
    lambdaih_test = 0.d0
    do i = 1, nelorb_c - 1
        lambdars_test(2*i - 1, i + 1) = lambdars(i + nelorb_c)
        lambdars_test(2*i + 1, i) = lambdars(i + nelorb_c)
        lambdaih_test(2*i, i + 1) = lambdars(i + nelorb_c)
        lambdaih_test(2*(i + 1), i) = -lambdars(i + nelorb_c)
    end do

    call ZGEMM("C", "N", nelorb_c, nelorb_c, nelorb_c, zone, U3, nelorb_c, &
               lambdars_test, nelorb_c, zzero, auxmat, nelorb_c)
    call ZGEMM("N", "N", nelorb_c, nelorb_c, nelorb_c, zone, auxmat, nelorb_c, &
               U3, nelorb_c, zzero, lambdars_test, nelorb_c)

    !call print_matrix(nelorb_c, 2, lambdars_test)

    lambdaih_test = lambdaih_test - lambdars_test
    write (6, *) "If everything is correct no output before <Check U3 Completed>"
    call print_matrix(nelorb_c, nelorb_c, 2, lambdaih_test)
    write (6, *) "Check U3 Completed"

    deallocate (lambdars_test, lambdaih_test, auxmat)

#endif

end subroutine symmtriang

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Using the pfapack library to tridiagonalize the
!matrix detmattr and to calculate U1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Tridiagonalize skew-symmetric matrix using Pfaffian decomposition
 *
 * This subroutine uses the Pfaffian decomposition algorithm to tridiagonalize
 * a skew-symmetric matrix. It employs the Pfaffian library routines (DSKTRD/ZSKTRD)
 * to reduce the matrix to tridiagonal form and compute the transformation matrix U1.
 *
 * @param[in] nelorb_c Number of orbitals
 * @param[in] ipc Complex flag (1=real, 2=complex)
 * @param[out] detmattr Tridiagonal skew-symmetric matrix
 * @param[in,out] U1 Transformation matrix (input: original matrix, output: orthogonal matrix)
 *
 * @details
 * The subroutine performs the following steps:
 * 1. Allocates work arrays for the tridiagonalization
 * 2. For real matrices (ipc=1):
 *    - Calls DSKTRD to tridiagonalize the matrix
 *    - Extracts off-diagonal elements to detmattr
 *    - Calls DORGTR to compute the orthogonal transformation matrix
 * 3. For complex matrices (ipc=2):
 *    - Calls ZSKTRD to tridiagonalize the matrix
 *    - Extracts real and imaginary parts of off-diagonal elements
 *    - Calls ZUNGTR to compute the unitary transformation matrix
 * 4. In debug mode: verifies the tridiagonalization is correct
 *
 * The tridiagonal matrix detmattr has the form:
 * - Diagonal elements are zero (skew-symmetric property)
 * - Off-diagonal elements contain the tridiagonal values
 * - For complex matrices, real and imaginary parts are stored separately
 *
 * @note
 * - Uses Pfaffian library routines (DSKTRD/ZSKTRD, DORGTR/ZUNGTR)
 * - Work array size is determined by workspace query
 * - Debug mode includes verification of the decomposition
 * - U1 matrix is used in the final back-transformation
 * - Supports both real and complex matrices
 *
 * @see symmtriang(), finalize_mopfaff()
 *
 * @author C.G.
 * @date 2022
 */
subroutine pfatriag(nelorb_c, ipc, detmattr, U1)
    implicit none
    integer :: nelorb_c, ipc
    real(8) :: detmattr(ipc*nelorb_c, nelorb_c), U1(ipc*nelorb_c, nelorb_c)
    integer :: i, lwork, info
    real(8) :: testr
    complex(8) :: testc
    real(8), allocatable :: aux(:), tau(:), work(:), auxmat(:, :)
    complex(8) :: zzero, zone
    zzero = (0.d0, 0.d0)
    zone = (1.d0, 0.d0)

    allocate (aux(ipc*(nelorb_c - 1)), tau(ipc*(nelorb_c - 1)))
    detmattr = 0.d0
    lwork = -1

#ifdef DEBUG
    allocate (auxmat(ipc*nelorb_c, nelorb_c))
    auxmat = U1
#endif

    if (ipc .eq. 1) then
        !Look for the documentation, this thing is a real mess
        call dsktrd("U", "N", nelorb_c, U1, nelorb_c, aux, tau, testr, lwork, info)
        lwork = idnint(testr)
        allocate (work(lwork))
        call dsktrd("U", "N", nelorb_c, U1, nelorb_c, aux, tau, work, lwork, info)
        if (info .ne. 0) write (6, *) "Parameter n", info, "of dsktrd is incorrect"

        do i = 1, nelorb_c - 1
            detmattr(i, i + 1) = U1(i, i + 1)
            detmattr(i + 1, i) = -U1(i, i + 1)
            U1(i, i + 1) = 1.d0
        end do
        call dorgtr("U", nelorb_c, U1, nelorb_c, tau, work, lwork, info)
        if (info .ne. 0) write (6, *) "Parameter n", info, "of dorgtr is incorrect"
    else
        !Look for the documentation, this thing is a real mess
        call zsktrd("U", "N", nelorb_c, U1, nelorb_c, aux, tau, testc, lwork, info)
        lwork = idnint(dreal(testc))
        allocate (work(lwork*ipc))
        call zsktrd("U", "N", nelorb_c, U1, nelorb_c, aux, tau, work, lwork, info)
        if (info .ne. 0) write (6, *) "Parameter n", info, "of dsktrd is incorrect"

        !     stop
        do i = 1, nelorb_c - 1
            detmattr(2*i - 1, i + 1) = U1(2*i - 1, i + 1)
            detmattr(2*i, i + 1) = U1(2*i, i + 1)
            detmattr(2*(i + 1) - 1, i) = -U1(2*i - 1, i + 1)
            detmattr(2*(i + 1), i) = -U1(2*i, i + 1)
            U1(2*i - 1, i + 1) = 1.d0
            U1(2*i, i + 1) = 0.d0
        end do

        call zungtr("U", nelorb_c, U1, nelorb_c, tau, work, lwork, info)
        if (info .ne. 0) write (6, *) "Parameter n", info, "of dorgtr is incorrect"
        !     call print_matrix(nelorb_c,2,U1)
    end if

#ifdef  DEBUG
    deallocate (aux)
    allocate (aux(ipc*nelorb_c*nelorb_c))
    if (ipc .eq. 1) then
        call DGEMM("N", "N", nelorb_c, nelorb_c, nelorb_c, 1.d0, U1, nelorb_c, &
                   detmattr, nelorb_c, 0.d0, aux, nelorb_c)
        call DGEMM("N", "T", nelorb_c, nelorb_c, nelorb_c, -1.d0, aux, nelorb_c, &
                   U1, nelorb_c, 1.d0, auxmat, nelorb_c)
    else
        call ZGEMM("N", "N", nelorb_c, nelorb_c, nelorb_c, zone, U1, nelorb_c, &
                   detmattr, nelorb_c, zzero, aux, nelorb_c)
        call ZGEMM("N", "T", nelorb_c, nelorb_c, nelorb_c, -zone, aux, nelorb_c, &
                   U1, nelorb_c, zone, auxmat, nelorb_c)
    end if
    write (6, *) "If everything is correct no output before <Check U1 Completed>"
    call print_matrix(nelorb_c, nelorb_c, ipc, auxmat)
    write (6, *) "Check U1 Completed"
    deallocate (auxmat)
#endif
    deallocate (aux, tau, work)
end subroutine pfatriag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Applying the transformations to prepare the output
!eigenvectors and eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Apply final transformations to obtain Pfaffian molecular orbitals
 *
 * This subroutine applies the final transformations to convert the diagonalized
 * eigenvectors back to the original basis, producing the final Pfaffian molecular
 * orbitals. It handles the pairing of eigenvalues and the construction of
 * proper molecular orbitals from the transformed eigenvectors.
 *
 * @param[in] lda Leading dimension of the output matrix
 * @param[in] nelorb_c Number of orbitals
 * @param[in] ipc Complex flag (1=real, 2=complex)
 * @param[in] U1 Orthogonal transformation matrix from tridiagonalization
 * @param[in] U3 Transformation matrix for real symmetric form
 * @param[in] eigvalues Eigenvalues from diagonalization
 * @param[in] eigvect Eigenvectors from diagonalization
 * @param[out] outvl Final eigenvalues (nelorb_c/2 values)
 * @param[out] outvct Final molecular orbitals (ipc*lda x nelorb_c)
 * @param[in] detmat_c Original skew-symmetric matrix (for debug verification)
 *
 * @details
 * The subroutine performs the following transformations:
 * 1. Applies U3^dagger to the eigenvectors to convert back from real symmetric form
 * 2. Handles eigenvalue pairing and selection:
 *    - Skew-symmetric matrices have paired eigenvalues ±λ
 *    - Only positive eigenvalues are returned in outvl
 *    - Eigenvectors are paired accordingly
 * 3. Applies U1 transformation to convert back to original basis
 * 4. Handles special cases for odd numbers of orbitals
 * 5. In debug mode: verifies the final result matches the original matrix
 *
 * The eigenvalue pairing logic:
 * - Normal pairs: eigenvalues with significant magnitude difference
 * - Singular pairs: eigenvalues corresponding to zero eigenvalues
 * - Single orbital: for odd numbers of orbitals, one unpaired orbital
 *
 * @note
 * - Uses LAPACK BLAS routines (DGEMM/ZGEMM) for matrix multiplications
 * - Handles both real and complex matrices
 * - Includes numerical stability checks for eigenvalue pairing
 * - Debug mode verifies the complete transformation
 * - Supports both even and odd numbers of orbitals
 *
 * @see pfaffian_mo(), pfatriag(), symmtriang()
 *
 * @author C.G.
 * @date 2022
 */
subroutine finalize_mopfaff(lda, nelorb_c, ipc, U1, U3, eigvalues, eigvect, outvl, outvct, detmat_c)
    implicit none
    integer :: nelorb_c, ipc, lda
    integer :: i, ind_even
    real(8) :: U1(ipc*nelorb_c, nelorb_c), U3(2*nelorb_c, nelorb_c), &
               eigvect(nelorb_c, nelorb_c), eigvalues(nelorb_c), outvct(ipc*lda, nelorb_c), &
               outvl(nelorb_c), detmat_c(ipc*lda, nelorb_c)
    real(8), allocatable :: auxvect(:, :), auxmat(:, :), auxvect2(:, :)
    complex(8) :: zzero, zone
    real(8) :: sqrt2
    real(8), external :: dlamch
    zzero = (0.d0, 0.d0)
    zone = (1.d0, 0.d0)
    allocate (auxvect(2*nelorb_c, nelorb_c), auxmat(2*nelorb_c, nelorb_c))

    outvl = 0.d0

    auxvect(1:2*nelorb_c - 1:2, 1:nelorb_c) = eigvect(1:nelorb_c, 1:nelorb_c)
    auxvect(2:2*nelorb_c:2, 1:nelorb_c) = 0.d0

    call ZGEMM("C", "N", nelorb_c, nelorb_c, nelorb_c, zone, U3, nelorb_c, &
               auxvect, nelorb_c, zzero, auxmat, nelorb_c)

    !Creating the vector needed for the tridiagonal rapresentation
    outvct = 0.d0
    sqrt2 = 2.d0**0.5d0
    if (ipc .eq. 1) then
        allocate (auxvect2(nelorb_c, nelorb_c))

        do i = 1, (nelorb_c + 1)/2
            ind_even = 2*(i - mod(nelorb_c, 2))
            outvl(i) = eigvalues(nelorb_c/2 + i)
            if (abs(outvl(i)) .gt. dlamch('S')*1000000 .and. .not. (i .eq. 1 .and. mod(nelorb_c, 2) .ne. 0)) then
                !      if it  is not  zero within the max cond number it should be condidered.
                !   It is only a poorly converged non zero eigenvalue. Hopefully never happens
                if (abs(abs((eigvalues(nelorb_c/2 + i) - eigvalues((nelorb_c + 1)/2 - i + 1))/ &
                        (2*eigvalues(nelorb_c/2 + i))) - 1.d0) .lt. dlamch('eps')*1000000 .or.&
                        &abs(eigvalues(nelorb_c/2 + i)/eigvalues(nelorb_c)) .gt. dlamch('eps')) then
                    auxvect2(1:nelorb_c, 2*i - 1) = auxmat(1:2*nelorb_c - 1:2, i + nelorb_c/2)*sqrt2
                    auxvect2(1:nelorb_c, ind_even) = auxmat(2:2*nelorb_c:2, i + nelorb_c/2)*sqrt2
                    !           write(6,*) ' Normal pair  put in/from = ',2*i-1,ind_even,i+nelorb_c/2
                else
                    !  consider this pair  singular eigenvectors corresponding to zero eigenvalue
                    auxvect2(:, 2*i - 1) = auxmat(1:2*nelorb_c - 1:2, i + nelorb_c/2) + &
                                           auxmat(2:2*nelorb_c:2, i + nelorb_c/2)
                    auxvect2(:, ind_even) = auxmat(1:2*nelorb_c - 1:2, 1 - i + (nelorb_c + 1)/2) + &
                                            auxmat(2:2*nelorb_c:2, 1 - i + (nelorb_c + 1)/2)
                    !           write(6,*) ' Singular pair put in/from = ',2*i-1,ind_even,i+nelorb_c/2,1-i+(nelorb_c+1)/2
                end if
            else
                if (i .eq. 1 .and. mod(nelorb_c, 2) .ne. 0) then
                    auxvect2(:, 2*i - 1) = auxmat(1:2*nelorb_c - 1:2, i + nelorb_c/2) + &
                                           auxmat(2:2*nelorb_c:2, i + nelorb_c/2)
                    !           write(6,*) ' Singular single put in/from = ',2*i-1,i+nelorb_c/2
                else
                    auxvect2(:, 2*i - 1) = auxmat(1:2*nelorb_c - 1:2, i + nelorb_c/2) + &
                                           auxmat(2:2*nelorb_c:2, i + nelorb_c/2)
                    auxvect2(:, ind_even) = auxmat(1:2*nelorb_c - 1:2, 1 - i + (nelorb_c + 1)/2) + &
                                            auxmat(2:2*nelorb_c:2, 1 - i + (nelorb_c + 1)/2)
                    !           write(6,*) ' Singular pair II put in/from = ',2*i-1,ind_even,i+nelorb_c/2,1-i+(nelorb_c+1)/2
                end if
            end if
        end do

        call DGEMM("N", "N", nelorb_c, nelorb_c, nelorb_c, 1.d0, U1, nelorb_c, &
                   auxvect2, nelorb_c, 0.d0, outvct, lda)
        deallocate (auxvect2)
    else
        do i = 1, (nelorb_c + 1)/2
            ind_even = 2*(i - mod(nelorb_c, 2))
            outvl(i) = eigvalues(nelorb_c/2 + i)
            if (abs(outvl(i)) .gt. dlamch('S')*1000000 .and. .not. (i .eq. 1 .and. mod(nelorb_c, 2) .ne. 0)) then
                !      if it  is not  zero within the max cond number it should be condidered.
                !   It  is only a poorly converged non zero eigenvalue. Hopefully never happens
                if (abs(abs((eigvalues(nelorb_c/2 + i) - eigvalues((nelorb_c + 1)/2 - i + 1))/ &
                        (2*eigvalues(nelorb_c/2 + i))) - 1.d0) .lt. dlamch('eps')*1000000 .or.&
                        &abs(eigvalues(nelorb_c/2 + i)/eigvalues(nelorb_c)) .gt. dlamch('eps')) then
                    auxvect(1:2*nelorb_c - 1:2, 2*i - 1) = auxmat(1:2*nelorb_c - 1:2, i + nelorb_c/2)*sqrt2
                    auxvect(1:2*nelorb_c - 1:2, ind_even) = auxmat(2:2*nelorb_c:2, i + nelorb_c/2)*sqrt2

                else
                    auxvect(1:2*nelorb_c:2, 2*i - 1)&
                       & = auxmat(1:2*nelorb_c:2, i + nelorb_c/2)&
                       & + auxmat(2:2*nelorb_c:2, i + nelorb_c/2)
                    auxvect(1:2*nelorb_c:2, ind_even)&
                       & = auxmat(1:2*nelorb_c:2, 1 - i + (nelorb_c + 1)/2)&
                       & + auxmat(2:2*nelorb_c:2, 1 - i + (nelorb_c + 1)/2)
                end if
            else
                if (i .eq. 1 .and. mod(nelorb_c, 2) .ne. 0) then
                    auxvect(1:2*nelorb_c:2, 2*i - 1)&
                       & = auxmat(1:2*nelorb_c:2, i + nelorb_c/2)&
                       & + auxmat(2:2*nelorb_c:2, i + nelorb_c/2)
                else
                    auxvect(1:2*nelorb_c:2, 2*i - 1)&
                       & = auxmat(1:2*nelorb_c:2, i + nelorb_c/2)&
                       & + auxmat(2:2*nelorb_c:2, i + nelorb_c/2)
                    auxvect(1:2*nelorb_c:2, ind_even)&
                       & = auxmat(1:2*nelorb_c:2, 1 - i + (nelorb_c + 1)/2)&
                       & + auxmat(2:2*nelorb_c:2, 1 - i + (nelorb_c + 1)/2)
                end if
            end if
        end do

        call ZGEMM("N", "N", nelorb_c, nelorb_c, nelorb_c, zone, U1, nelorb_c, &
                   auxvect, nelorb_c, zzero, outvct, lda)

    end if

#ifdef DEBUG
    allocate (auxvect2(ipc*nelorb_c, nelorb_c))
    !Testing if the matrix has been built correctly
    deallocate (auxmat)
    allocate (auxmat(ipc*lda, nelorb_c))
    if (ipc .eq. 1) then
        auxmat(:, :) = 0.d0
        do i = 1, nelorb_c/2
            auxmat(2*i - 1, 2*i) = -outvl(i)
            auxmat(2*i, 2*i - 1) = outvl(i)
        end do
        call DGEMM("N", "N", nelorb_c, nelorb_c, nelorb_c, 1.d0, outvct, lda, &
                   auxmat, lda, 0.d0, auxvect2, nelorb_c)
        call DGEMM("N", "T", nelorb_c, nelorb_c, nelorb_c, 1.d0, auxvect2, nelorb_c, &
                   outvct, lda, 0.d0, auxmat, lda)
    else
        auxmat(:, :) = 0.d0
        do i = 1, nelorb_c/2
            auxmat(2*(2*i - 1) - 1, 2*i) = -outvl(i)
            auxmat(2*(2*i) - 1, 2*i - 1) = outvl(i)
        end do
        call ZGEMM("N", "N", nelorb_c, nelorb_c, nelorb_c, zone, outvct, lda, &
                   auxmat, lda, zzero, auxvect2, nelorb_c)
        call ZGEMM("N", "T", nelorb_c, nelorb_c, nelorb_c, zone, auxvect2, nelorb_c, &
                   outvct, lda, zzero, auxmat, lda)
    end if

    auxmat = auxmat - detmat_c

    write (6, *) "If everything is correct no output before <Choice Completed>"
    call print_matrix(lda, nelorb_c, ipc, auxmat)
    write (6, *) "<Choice Completed>"
    deallocate (auxvect2)
#endif

    deallocate (auxvect, auxmat)
end subroutine finalize_mopfaff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Auxiliary subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Filling a skew symmetric tridiagonal matrix for tests
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Fill a matrix with random skew-symmetric values for testing
 *
 * This subroutine fills a matrix with random skew-symmetric values for
 * testing and debugging purposes. It generates random numbers and
 * constructs a skew-symmetric matrix with the property A(i,j) = -A(j,i).
 *
 * @param[in] nelorb_c Number of orbitals
 * @param[in] ipc Complex flag (1=real, 2=complex)
 * @param[out] detmat_c Skew-symmetric matrix filled with random values
 *
 * @details
 * The subroutine generates random values in the range [-1, 1] and
 * constructs a skew-symmetric matrix:
 * - For real matrices (ipc=1): A(i,j) = random value, A(j,i) = -A(i,j)
 * - For complex matrices (ipc=2): real and imaginary parts are generated separately
 * - Upper triangular elements are filled with random values
 * - Lower triangular elements are set to negative of corresponding upper elements
 * - Diagonal elements remain zero (skew-symmetric property)
 *
 * @note
 * - Uses Fortran random_number() for random value generation
 * - Sets seed to 4 for reproducible results
 * - Only fills elements where i+j > nelorb_c for tridiagonal-like structure
 * - Used primarily for testing and debugging
 *
 * @see fill_tridiag(), print_matrix()
 *
 * @author C.G.
 * @date 2022
 */
subroutine fill_skw(nelorb_c, ipc, detmat_c)
    implicit none
    integer :: nelorb_c, ipc
    integer :: i, j, seed !Auxiliary variables
    real(8) :: x, detmat_c(ipc*nelorb_c, nelorb_c)
    seed = 4
    call random_seed(seed)
    detmat_c = 0.d0
    if (ipc .eq. 1) then
        do i = 1, nelorb_c
            do j = i + 1, nelorb_c
                call random_number(x)
                if (i + j .gt. nelorb_c) x = 0.d0
                detmat_c(j, i) = (x - 0.5)*2
                detmat_c(i, j) = -detmat_c(j, i)
            end do
        end do
    else
        do i = 1, nelorb_c
            do j = i + 1, nelorb_c
                call random_number(x)
                if (i + j .gt. nelorb_c) x = 0.d0
                detmat_c(2*j - 1, i) = (x - 0.5)*2
                detmat_c(2*i - 1, j) = -detmat_c(2*j - 1, i)
                call random_number(x)
                if (i + j .gt. nelorb_c) x = 0.d0
                detmat_c(2*j, i) = (x - 0.5)*2
                detmat_c(2*i, j) = -detmat_c(2*j, i)
            end do
        end do

    end if
end subroutine fill_skw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Filling a skew symmetric tridiagonal matrix for tests
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Fill a matrix with random tridiagonal skew-symmetric values for testing
 *
 * This subroutine fills a matrix with random tridiagonal skew-symmetric values
 * for testing and debugging purposes. It generates a matrix with non-zero
 * elements only on the first superdiagonal and subdiagonal.
 *
 * @param[in] nelorb_c Number of orbitals
 * @param[in] ipc Complex flag (1=real, 2=complex)
 * @param[out] detmattr Tridiagonal skew-symmetric matrix filled with random values
 *
 * @details
 * The subroutine constructs a tridiagonal skew-symmetric matrix:
 * - Diagonal elements are zero (skew-symmetric property)
 * - Superdiagonal elements (i, i+1) are filled with random values
 * - Subdiagonal elements (i+1, i) are set to negative of corresponding superdiagonal
 * - For complex matrices (ipc=2): real and imaginary parts are generated separately
 * - All other elements remain zero
 *
 * @note
 * - Uses Fortran random_number() for random value generation
 * - Sets seed to 4 for reproducible results
 * - Creates a proper tridiagonal structure
 * - Used primarily for testing tridiagonalization algorithms
 *
 * @see fill_skw(), print_matrix()
 *
 * @author C.G.
 * @date 2022
 */
subroutine fill_tridiag(nelorb_c, ipc, detmattr)
    implicit none
    integer :: nelorb_c, ipc
    integer :: i, j, seed !Auxiliary variables
    real(8) :: x, detmattr(ipc*nelorb_c, nelorb_c)
    seed = 4
    call random_seed(seed)

    do i = 1, nelorb_c
        do j = 1, ipc*nelorb_c
            detmattr(j, i) = 0.d0
        end do
        j = i + 1
        if (j .le. nelorb_c) then
            if (ipc .eq. 2) then
                call random_number(x)
                detmattr(2*j - 1, i) = (x - 0.5)*2
            end if
            call random_number(x)
            detmattr(ipc*j, i) = (x - 0.5)*2
        end if
        j = i - 1
        if (j .ne. 0) then
            detmattr(ipc*j, i) = -detmattr(ipc*i, j)
            if (ipc .eq. 2) detmattr(ipc*j - 1, i) = -detmattr(ipc*i - 1, j)
        end if
    end do

end subroutine fill_tridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine that  matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Print matrix elements above a threshold for debugging
 *
 * This subroutine prints matrix elements that exceed a specified threshold
 * value, useful for debugging and verification of matrix operations.
 * Only elements with magnitude greater than the threshold are printed.
 *
 * @param[in] lda Leading dimension of the matrix
 * @param[in] nelorb_c Number of orbitals (columns)
 * @param[in] ipc Complex flag (1=real, 2=complex)
 * @param[in] detmat_c Matrix to be printed
 *
 * @details
 * The subroutine prints matrix elements in the format:
 * - For real matrices (ipc=1): "j i value"
 * - For complex matrices (ipc=2): "j i real_part imaginary_part"
 * - Only elements with magnitude > prec are printed
 * - Default precision threshold is 1e-7
 *
 * @note
 * - Used primarily for debugging matrix operations
 * - Helps verify matrix structure and values
 * - Suppresses printing of very small elements
 * - Output format is suitable for manual inspection
 *
 * @see fill_skw(), fill_tridiag()
 *
 * @author C.G.
 * @date 2022
 */
subroutine print_matrix(lda, nelorb_c, ipc, detmat_c)
    implicit none
    integer :: nelorb_c, lda, ipc
    integer :: i, j !Auxiliary variables
    real(8) :: detmat_c(ipc*lda, nelorb_c)
    real(8) :: prec

    prec = 1d-7

    do i = 1, nelorb_c
        do j = 1, nelorb_c
            if (ipc .eq. 2) then
                if (abs(detmat_c(2*j, i)) + abs(detmat_c(2*j - 1, i)) .gt. prec) &
                    write (6, *) j, i, detmat_c(2*j - 1, i), detmat_c(2*j, i)
            else
                if (abs(detmat_c(j, i)) .gt. prec) &
                    write (6, *) j, i, detmat_c(j, i)

            end if
        end do
    end do
end subroutine print_matrix
