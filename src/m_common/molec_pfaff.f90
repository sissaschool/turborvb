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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutines for the calculation of the pfaffians
!molecular orbitals            C.G.
!
!Output: outvl (nelorb_c/2) eigenvalues sorted by the
!        value (only the real part)
!Output: outvct(ipc*nelorb_c,nelorb_c) eigenvectors
!        ordered as egvl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine that  matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
