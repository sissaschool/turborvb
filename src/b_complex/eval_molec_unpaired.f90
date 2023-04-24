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

subroutine eval_molec_unpaired(nelorb_c, overs, mat_in        &
        &, molecorb, eig, lda, eps, nproc, rank, rank_mpi, comm_mpi, ndiff, psi_unp, orthoyes, symmagp)
    use constants, only: ipc, ipf, zzero, zone, zmone
    use allio, only: yes_hermite, symmetrize_agp&
            &, opposite_phase, same_phase, real_contracted, gauge_fixing
    implicit none
    integer nelorb_c, i, j, nproc, rank, info, n, lda, mine, neig, lwork, ndiff&
            &, minemax, dimo, mino, dimn, dimorb, ierr, miner, dimor, il, iu, indi, indj&
            &, indunp, comm_mpi, rank_mpi, nshift
    real*8 overs(ipc*lda, *), mat_in(ipc*lda, *)                       &
            &, molecorb(ipc*lda, *), eig(nelorb_c), psi_unp(ipc*nelorb_c, ndiff)     &
            &, abstol, eps, orbmax, dlamch, ddot, eigmax, cost, dnrm2, vl, vu
    real(8), dimension(:, :), allocatable :: mat, umat, umatright
    real(8), dimension(:), allocatable :: eig_sav, work, eigmat, eigright, eigmr, psint
    integer, dimension(:), allocatable :: iwork
    logical orthoyes, symmagp, eqover, yes_constrainm, yesh
#ifdef PARALLEL
    include 'mpif.h'
#ifdef __TEST
    if (nproc .gt. 1) then
        call bcast_real(overs, ipc*lda*nelorb_c, 0, comm_mpi)
        call bcast_real(mat_in, ipc*lda*nelorb_c, 0, comm_mpi)
        call bcast_real(psi_unp, ipc*ndiff*nelorb_c, 0, comm_mpi)
    end if
#endif
#endif
    !        integer iwork(5*nelorb_c),ifail(nelorb_c)

    !        if(rank_mpi.eq.0.or.nproc.gt.1) then

    lwork = 6*nelorb_c
    allocate (mat(ipc*nelorb_c, nelorb_c), umat(ipc*nelorb_c, nelorb_c), &
            &eig_sav(nelorb_c), work(lwork)&
            &, eigmat(nelorb_c), iwork(5*nelorb_c), psint(nelorb_c))
    vl = 0.d0
    vu = 0.d0
    il = 0
    iu = 0
    if (ipc .eq. 2 .and. yes_hermite .and. symmagp .and. symmetrize_agp .and. ipf .eq. 1) then
        yes_constrainm = .true.
    else
        yes_constrainm = .false.
    end if

    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then
        eqover = .true.
        if (ipc .eq. 1) then
            do i = 1, nelorb_c
                do j = 1, nelorb_c
                    if (overs(i, j) .ne. overs(i, j + nelorb_c)) eqover = .false.
                end do
            end do
        else
            eqover = .false.
        end if
        allocate (umatright(ipc*nelorb_c, nelorb_c))
        allocate (eigright(ipc*nelorb_c), eigmr(nelorb_c))
    end if

    ! destroy invo
    !        calculation overlap old
    !        calculation normalized orbitals and eigenvectors               RR
    ! maximum accuracy
    abstol = 0.d0
    if (eps .eq. 0.d0) eps = 10.d0*dlamch('E')

    !    This code compute the generalized eigenvalue equation:
    !      mat_in x over_s  y = lambda y
    !      by stable preconditioning and pre-diagonalization of the matrix over_s

    !         There is numerical instability to use premat ne 1
    !         the orthogonality of the eigenvectors
    !         with respect to the metric given by overs  is poorly satisfied
    !         if premat ne 1 is used.

    !         write(6,*) ' premat chosen '
    !         do i=1,nelorb_c
    !         if(overs(i,i).ne.0.d0) then
    !         premat(i)=1.d0/dsqrt(abs(overs(i,i)))
    !         else
    !         if(rank.eq.0) write(6,*) ' Singular norm for orbital ',i,overs(i,i)
    !         premat(i)=1.d0
    !         stop
    !         endif
    !         write(6,*) i,premat(i)
    !         enddo

    n = nelorb_c
    !        if(rank.eq.0) then
    !         write(6,*) ' Input matrices ',n
    !         do i=1,n
    !           do j=i,n
    !           write(6,*) i,j,mat_in(i,j),overs(i,j)
    !           enddo
    !         enddo
    !        endif
    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then
        !        ARROCCO input
        do i = 1, n
            umatright(1:n*ipc, i) = overs(1:n*ipc, i)
            !        In the following the spin up (left) orbitals are in overs(:,n+1:2n)
            !                         the spin-down (right) ==        in overs(:,1:n)
            overs(1:n*ipc, i) = overs(1:n*ipc, n + i)
            overs(1:n*ipc, n + i) = umatright(1:n*ipc, i)
        end do
    end if

    if (ndiff .gt. 0) then
        if (ipc .eq. 1) then
            psint(:) = psi_unp(:, 1)**2
            do j = 2, ndiff
                psint(:) = psint(:) + psi_unp(:, j)**2
            end do
        else
            do i = 1, n
                psint(i) = psi_unp(2*i - 1, 1)**2 + psi_unp(2*i, 1)**2
            end do
            do j = 2, ndiff
                do i = 1, n
                    psint(i) = psint(i) + psi_unp(2*i - 1, j)**2 + psi_unp(2*i, j)**2
                end do
            end do
        end if
    else
        psint(:) = 0.d0
    end if

    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        if (ipf .eq. 2) then
            indi = n
            indj = n
            eig_sav = 0.d0
            do i = 1, n
                do j = i, n
                    if (mat_in(i, j) .ne. 0.d0) then
                        eig_sav(i) = 1.d0
                        eig_sav(j) = 1.d0
                    end if
                end do
            end do

            do j = 1, n
                do i = 1, n
                    if (i .ne. j) then
                        overs((i - 1)*ipc + 1:ipc*i, j) = overs((i - 1)*ipc + 1:ipc*i, j)*eig_sav(i)*eig_sav(j)
                    elseif (eig_sav(i) .eq. 0.d0) then
                        overs((i - 1)*ipc + 1:ipc*i, j) = 1.d0
                        if (ipc .eq. 2) overs(2*i, j) = 0.d0
                    end if
                end do
            end do
            eig_sav = 1.d0
        else
            eig_sav = 0.d0
            do i = 1, n
                do j = i, n
                    if (mat_in(i, j) .ne. 0.d0 .or. psint(i) .ne. 0.d0 .or. psint(j) .ne. 0.d0) then
                        eig_sav(i) = 1.d0
                        eig_sav(j) = 1.d0
                    end if
                end do
            end do
            indi = 0
            do i = 1, n
                if (eig_sav(i) .eq. 1.d0) then
                    indi = indi + 1
                    indj = 0
                    do j = 1, n
                        if (eig_sav(j) .eq. 1.d0) then
                            indj = indj + 1
                            overs(indi, indj) = overs(i, j)
                        end if
                    end do
                end if
            end do
        end if ! ipf=2

    else

        !     Determine the non zero left  basis elements for the up spin
        eig_sav = 0.d0
        do i = 1, n
            do j = 1, n
                if (sum(abs(mat_in(ipc*(i - 1) + 1:ipc*i, j))) .ne. 0.d0 .or. psint(i) .ne. 0.d0) then
                    eig_sav(i) = 1.d0
                end if
            end do
        end do
        indi = 0
        do i = 1, n
            if (eig_sav(i) .eq. 1.d0) then
                indi = indi + 1
                indj = 0
                do j = 1, n
                    if (eig_sav(j) .eq. 1.d0) then
                        indj = indj + 1
                        if (ipc .eq. 1) then
                            overs(indi, indj + n) = overs(i, j + n)
                        else
                            overs(2*indi - 1, indj + n) = overs(2*i - 1, j + n)
                            overs(2*indi, indj + n) = overs(2*i, j + n)
                        end if
                    end if
                end do
            end if
        end do
    end if

    !        over_sav(1:n,1:n)=overs(1:n,1:n)
    if (ipc .eq. 2 .and. ipf .eq. 1) then
        call ZSYEV_MY('V', 'L', indi, overs(1, n + 1), LDA, eig, INFO&
                &, nproc, rank_mpi, comm_mpi)
        umat(1:2*indi, 1:indi) = overs(1:2*indi, n + 1:n + indi)
    else
        if (symmagp .or. ipf .eq. 2) then
            if (ipc .eq. 1) then
                call DSYEV_MY('V', 'L', indi, overs, LDA, eig, INFO&
                        &, nproc, rank_mpi, comm_mpi)
            else ! below only ipf=2 is possible
                call ZSYEV_MY('V', 'L', indi, overs, LDA, eig, INFO&
                        &, nproc, rank_mpi, comm_mpi)
            end if
            umat(1:ipc*indi, 1:indi) = overs(1:ipc*indi, 1:indi)
        else
            call DSYEV_MY('V', 'L', indi, overs(1, n + 1), LDA, eig, INFO&
                    &, nproc, rank_mpi, comm_mpi)
            umat(1:indi, 1:indi) = overs(1:indi, n + 1:n + indi)
        end if
    end if
    work(1:indj) = eig(1:indj)
    eig = -1.d0
    eig(n - indj + 1:n) = work(1:indj)

    !         Reordering  without loosing information
    do j = indj, 1, -1
        work(1:ipc*indj) = umat(1:ipc*indj, j)
        umat(:, n - indj + j) = 0.d0
        indi = 0
        do i = 1, n
            if (eig_sav(i) .eq. 1.d0) then
                indi = indi + 1
                umat(ipc*(i - 1) + 1:ipc*i, n - indj + j) = work(ipc*(indi - 1) + 1:ipc*indi)
            end if
        end do
    end do
    do i = 1, n - indj
        umat(1:ipc*n, i) = 0.d0
    end do

    mine = 1

    !         write(6,*) ' Eigenvalues overlap matrix '

    do i = 1, n
        !         write(6,*) i,eig(i)
        if (eig(i)/eig(n) .gt. eps) then ! the condition number criterium
            eigmat(i) = dsqrt(1.d0/eig(i))
        else
            if (ipf .eq. 1) mine = i + 1
            if (ipf .eq. 2) then
                eigmat(i) = dsqrt(1.d0/(eps*eig(n)))
            else
                eigmat(i) = 0.d0
            end if
            !        write(6,*) ' warning zero eigenvalue !!! '
        end if
    end do
    if (info .ne. 0 .and. rank .eq. 0) write (6, *)                         &
            &' info > 0 in dsyevx !!! ', info

    if (mine .ne. 1 .and. rank .eq. 0)                                     &
            &write (6, *) ' disregarded coll. =', mine - 1

    minemax = max(mine - 1, ndiff)

    !      Calculation unpaired orbitals in the orthogonal basis diagonalizing overs

    mat = 0.d0

    !  mat contains the coefficient of the unpaired in the orthogonal left basis
    !  defined by umat(:,i)*sqrt(lambda_i), where lambda_i are the eigenvalues of
    !  the overlap matrix
    if (ipc .eq. 2) then
        call zgemm_my('C', 'N', n, ndiff, n, zone, umat, n, psi_unp, n, zzero, mat&
                &, n, nproc, rank_mpi, comm_mpi)
    else
        call dgemm_my('T', 'N', n, ndiff, n, 1.d0, umat, n, psi_unp, n, 0.d0, mat&
                &, n, nproc, rank_mpi, comm_mpi)
    end if
    dimo = n - mine + 1

    do j = 1, ndiff
        do i = 1, n
            if (eigmat(i) .ne. 0.d0) then
                mat(ipc*(i - 1) + 1:i*ipc, j) = mat(ipc*(i - 1) + 1:i*ipc, j)/eigmat(i)
            else
                mat(ipc*(i - 1) + 1:ipc*i, j) = 0.d0
            end if
        end do
    end do

    if (.not. orthoyes .or. ndiff .eq. 1) then
        !       only normalization
        do j = 1, ndiff
            cost = dnrm2(ipc*n, psi_unp(:, j), 1)
            if (cost .gt. 0.d0) psi_unp(:, j) = psi_unp(:, j)/cost
        end do
    end if

    !       orthogonalize this ndiff  vectors only (numerically stable)

    if (ipc .eq. 1) then
        call grahamo(mat, work(ndiff + 1), work, n, ndiff, indunp)
    else
        call grahamo_complex(mat, work(2*ndiff + 1), work, n, ndiff, indunp)
    end if
    if (indunp .ne. 0 .and. rank .eq. 0) write (6, *) &
            &'Warning dependency, unpaired not changed', indunp
    !       changing psi_unp back in the original space
    do i = 1, ndiff
        do j = 1, n
            molecorb(ipc*(j - 1) + 1:ipc*j, i) = mat(ipc*(j - 1) + 1:ipc*j, i)*eigmat(j)
        end do
    end do

    if (orthoyes .and. indunp .eq. 0) then
        if (ipc .eq. 1) then
            call dgemm_my('N', 'N', n, ndiff, dimo, 1.d0, umat(1, mine)&
                    &, n, molecorb(mine, 1), lda, 0.d0, psi_unp, n, nproc, rank_mpi, comm_mpi)
        else
            call zgemm_my('N', 'N', n, ndiff, dimo, zone, umat(1, mine)&
                    &, n, molecorb(2*mine - 1, 1), lda, zzero, psi_unp, n, nproc, rank_mpi, comm_mpi)
        end if
    end if

    !       For test removing this part
    !       mat=0.d0
    !       do i=1,n
    !       mat(2*i-1,i)=1.d0
    !       enddo
    !   added old I

    !      if(.not.symmagp.or.ipf.eq.2)  then ! This algorithm destroy symmetry if imposed

    molecorb(:, 1:n) = 0.d0
    do j = mine, n
        molecorb(ipc*(j - 1) + 1, j) = -1.d0
    end do

    !       Now define a fictitious matrix with only ndiff eigenvalues
    !       equal to -2 the other -1 (orthogonal)

    if (indunp .eq. 0) then
        if (ipc .eq. 1) then
            call dgemm_my('N', 'T', dimo, dimo, ndiff, -1.d0, mat(mine, 1), n&
                    &, mat(mine, 1), n, 1.d0, molecorb(mine, mine), lda, nproc, rank_mpi, comm_mpi)
        else
            call zgemm_my('N', 'C', dimo, dimo, ndiff, zmone, mat(2*mine - 1, 1), n&
                    &, mat(2*mine - 1, 1), n, zone, molecorb(2*mine - 1, mine), lda, nproc, rank_mpi, comm_mpi)
        end if
    end if

    mat = 0.d0

    if (ipc .eq. 1) then
        call DSYEV_MY('V', 'L', dimo, molecorb(mine, mine), LDA, eig, INFO&
                &, nproc, rank_mpi, comm_mpi)
    else
        call ZSYEV_MY('V', 'L', dimo, molecorb(2*mine - 1, mine), LDA, eig, INFO&
                &, nproc, rank_mpi, comm_mpi)
    end if
    !  Now the highest  eigenvectors are automatically orthogonal to the unpaired
    mat(ipc*(mine - 1) + 1:ipc*(mine + dimo - 1), mine:mine + dimo - 1) = &
            &molecorb(ipc*(mine - 1) + 1:ipc*(mine + dimo - 1), mine:mine + dimo - 1)

    !       else
    !       mat=0.d0
    !       do i=1,n
    !       mat(ipc*(i-1)+1,i)=1.d0
    !       enddo
    !       endif

    !       end added old I

    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then
        if (eqover) then
            umatright = umat
            eigright = eig
            eigmr = eigmat
            miner = mine
        else
            !  The right spin-down  eigenvectors
            work(1:n) = eig_sav(1:n)
            eig_sav = 0.d0
            do i = 1, n
                do j = 1, n
                    if (sum(abs(mat_in(ipc*(j - 1) + 1:ipc*j, i))) .ne. 0.d0) then
                        eig_sav(i) = 1.d0
                    end if
                end do
            end do
            indi = 0
            do i = 1, n
                if (eig_sav(i) .eq. 1.d0) then
                    indi = indi + 1
                    indj = 0
                    do j = 1, n
                        if (eig_sav(j) .eq. 1.d0) then
                            indj = indj + 1
                            if (ipc .eq. 1) then
                                overs(indi, indj) = overs(i, j)
                            else
                                overs(2*indi - 1, indj) = overs(2*i - 1, j)
                                overs(2*indi, indj) = -overs(2*i, j) ! the right have complex conjugate overlap
                            end if
                        end if
                    end do
                end if
            end do
            if (yes_constrainm) then
                yesh = .true.
                do i = 1, n
                    if (eig_sav(i) .ne. work(i)) yesh = .false.
                end do
            end if
            if (rank .eq. 0) then
                if (yes_constrainm) then
                    write (6, *) ' Warning constraining molecular orbitals'
                else
                    write (6, *) ' Warning NO constraining molecular orbitals'
                end if
            end if
            if (indi .gt. 0 .and.&
                    &((.not. opposite_phase .and. .not. same_phase) .or. .not. real_contracted)) then
                ! protection from fully polarized state.
                if (ipc .eq. 2) then
                    call ZSYEV_MY('V', 'L', indi, overs, LDA, eigright, INFO&
                            &, nproc, rank_mpi, comm_mpi)
                else
                    call DSYEV_MY('V', 'L', indi, overs, LDA, eigright, INFO&
                            &, nproc, rank_mpi, comm_mpi)
                end if
                !    right should be changed later
                do i = 1, indi
                    umatright(1:ipc*indi, i) = overs(1:ipc*indi, i)
                end do
            end if

            if ((.not. opposite_phase .and. .not. same_phase) .or. .not. real_contracted) then
                if (indj .gt. 0) then
                    work(1:indj) = eigright(1:indj)
                    eigright = -1.d0
                    eigright(n - indj + 1:n) = work(1:indj)
                else
                    eigright = -1.d0
                end if

                do j = indj, 1, -1
                    work(1:ipc*indj) = umatright(1:ipc*indj, j)
                    umatright(:, n - indj + j) = 0.d0
                    indi = 0
                    do i = 1, n
                        if (eig_sav(i) .eq. 1.d0) then
                            indi = indi + 1
                            if (ipc .eq. 1) then
                                umatright(i, n - indj + j) = work(indi)
                            else
                                umatright(2*i - 1, n - indj + j) = work(2*indi - 1)
                                umatright(2*i, n - indj + j) = work(2*indi)
                            end if
                        end if
                    end do
                end do
                do i = 1, n - indj
                    umatright(1:ipc*n, i) = 0.d0
                end do

                miner = 1

                do i = 1, n
                    if (eigright(i)/eigright(n) .gt. eps) then ! the condition number criterium

                        eigmr(i) = dsqrt(1.d0/eigright(i))
                    else
                        miner = i + 1
                        !        write(6,*) ' warning zero eigenvalue !!! '
                        eigmr(i) = 0.d0
                    end if
                end do
                if (info .ne. 0 .and. rank .eq. 0) write (6, *)                         &
                        &' info > 0 in dsyevx !!! ', info
            else
                if (opposite_phase .or. ipc .eq. 1) then
                    umatright = umat
                else
                    do i = 1, n
                        do j = 1, n
                            umatright(2*j - 1, i) = umat(2*j - 1, i)
                            umatright(2*j, i) = -umat(2*j, i)
                        end do
                    end do
                end if
                miner = mine
                eigmr = eigmat
                eigright = eig
            end if
        end if ! endif eqover
    end if ! endif symmagp

    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        miner = mine
        dimor = dimo
    else
        dimor = n - miner + 1
    end if

    !        destroy overs mat_in

    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        if (ipc .eq. 2) then
            if (ipf .eq. 2) call conjmat(n, n, umat, n)
            call zgemm_my('N', 'N', n, n, n, zone, mat_in, lda, umat, n, zzero, overs&
                    &, lda, nproc, rank_mpi, comm_mpi)
            if (ipf .eq. 2) call conjmat(n, n, umat, n)
        else
            call dgemm_my('N', 'N', n, n, n, 1.d0, mat_in, lda, umat, n, 0.d0, overs&
                    &, lda, nproc, rank_mpi, comm_mpi)
        end if
    else
        if (ipc .eq. 2) then
            call zgemm_my('N', 'N', n, n, n, zone, mat_in, lda, umatright, n, zzero, overs&
                    &, lda, nproc, rank_mpi, comm_mpi)
        else
            call dgemm_my('N', 'N', n, n, n, 1.d0, mat_in, lda, umatright, n, 0.d0, overs&
                    &, lda, nproc, rank_mpi, comm_mpi)
        end if
    end if
    if (ipc .eq. 2) then
        call zgemm_my('C', 'N', n, n, n, zone, umat, n, overs, lda, zzero, mat_in&
                &, lda, nproc, rank_mpi, comm_mpi)
    else
        call dgemm_my('T', 'N', n, n, n, 1.d0, umat, n, overs, lda, 0.d0, mat_in&
                &, lda, nproc, rank_mpi, comm_mpi)
    end if

    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        do i = 1, n
            do j = 1, n
                if (eigmat(j)*eigmat(i) .ne. 0.d0) then
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = mat_in(ipc*(i - 1) + 1:ipc*i, j)/eigmat(i)/eigmat(j)
                else
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = 0.d0
                end if
            end do
        end do
    else
        do i = 1, n
            do j = 1, n
                if (eigmr(j)*eigmat(i) .ne. 0.d0) then
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = mat_in(ipc*(i - 1) + 1:ipc*i, j)/eigmat(i)/eigmr(j)
                else
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = 0.d0
                end if
            end do
        end do
    end if
    ! addded back II

    !        destroy overs mat_in

    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        if (ipc .eq. 1) then
            call dgemm_my('N', 'N', n, n, n, 1.d0, mat_in, lda, mat, n, 0.d0, overs&
                    &, lda, nproc, rank_mpi, comm_mpi)
        else
            call conjmat(n, n, mat, n)
            call zgemm_my('N', 'N', n, n, n, zone, mat_in, lda, mat, n, zzero, overs&
                    &, lda, nproc, rank_mpi, comm_mpi)
            call conjmat(n, n, mat, n)
        end if
    else
        overs(1:n*ipc, 1:n) = mat_in(1:n*ipc, 1:n)
    end if
    ! Further changing basis to take into account the orthogonalization with
    ! the unpaired
    if (ipc .eq. 2) then
        call zgemm_my('C', 'N', n, n, n, zone, mat, n, overs, lda, zzero, mat_in&
                &, lda, nproc, rank_mpi, comm_mpi)
    else
        call dgemm_my('T', 'N', n, n, n, 1.d0, mat, n, overs, lda, 0.d0, mat_in&
                &, lda, nproc, rank_mpi, comm_mpi)
    end if

    ! end added back II

    if (indunp .ne. 0) then
        minemax = mine - 1
    else
        if (.not. symmagp .or. ipf .eq. 2) then
            minemax = mine - 1 + ndiff
        else
            minemax = mine - 1
        end if
    end if

    dimo = n - minemax
    mino = minemax + 1
    if ((.not. symmagp .or. ipc .eq. 2 .or. ipf .eq. 2) .and. minemax .gt. 0) then

        if (yes_constrainm .or. ipf .eq. 2) then
            mat_in(1:ipc*minemax, 1:n) = 0.d0
            mat_in(1:ipc*n, 1:minemax) = 0.d0
        else
            mat_in(1:ipc*minemax, 1:n) = 0.d0
        end if
    end if

    if (symmagp .and. ipc .eq. 1 .or. ipf .eq. 2) then
        molecorb(:, 1:n) = 0.d0
        if (mino .gt. 1) eig(1:minemax) = 0.d0

        if (ipf .eq. 2) then

            call pfaffian_mo(LDA, nelorb_c, ipc, mat_in, eig, molecorb)

        else
            call DSYEV_MY('V', 'L', dimo, mat_in(mino, mino), LDA, eig(mino), INFO&
                    &, nproc, rank_mpi, comm_mpi)
            molecorb(mino:mino + dimo - 1, mino:mino + dimo - 1) = &
                    & mat_in(mino:mino + dimo - 1, mino:mino + dimo - 1)
        end if

    else

        molecorb(1:ipc*n, n + 1:2*n) = mat_in(1:ipc*n, 1:n)
        if (ipc .eq. 2) then
            if (yes_constrainm .and. yesh) then
                call ZSYEV_MY('V', 'L', N, molecorb(1, n + 1), LDA, eig, INFO&
                        &, nproc, rank_mpi, comm_mpi)
                !        eig(:)=eig(:)**2
                do i = 1, n
                    do j = 1, n
                        molecorb(2*i - 1, j) = molecorb(2*j - 1, i + n)
                        molecorb(2*i, j) = molecorb(2*j, i + n)
                    end do
                end do
            else

                call ZGESVD_MY('V', N, N, molecorb(1, n + 1), LDA, molecorb, eig, nproc, rank_mpi, comm_mpi)
                if (gauge_fixing) then
                    do i = 1, N
                        call gauge_fix(N, molecorb(1, n + i), 1, molecorb(2*i - 1, 1), LDA)
                    end do
                end if

            end if
        else ! ipc=2
            if (nproc .gt. 1) then
                call DGESVD_MY('V', N, N, molecorb(1, n + 1), LDA, molecorb, eig, nproc, rank_mpi, comm_mpi)
            else
                call dgesvd('A', 'A', n, n, mat_in, lda&
                        &, eig, molecorb(1, n + 1), lda, molecorb, lda, work, lwork, info)
            end if
            if (gauge_fixing) then
                do i = 1, N
                    call gauge_fixr(N, molecorb(1, n + i), 1, molecorb(i, 1), LDA)
                end do
            end if
        end if ! ipc=2
    end if ! symmagp

    !      now go back to the original representation changing umat
    !      First the left spin up eigenvectors
    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then

        !    added old III
        overs(1:ipc*n, 1:n) = molecorb(1:ipc*n, 1:n)
        if (ipf .eq. 2) then
            if (ipc .eq. 1) then
                call dgemm_my('N', 'N', n, n, n, 1.d0, mat, n               &
                        &, overs, lda, 0.d0, molecorb, lda, nproc, rank_mpi, comm_mpi)
            else
                call zgemm_my('N', 'N', n, n, n, zone, mat, n               &
                        &, overs, lda, zzero, molecorb, lda, nproc, rank_mpi, comm_mpi)
            end if
        else
            call dgemm_my('N', 'N', n, dimo, dimo, 1.d0, mat(1, mino), n               &
                    &, overs(mino, mino), lda, 0.d0, molecorb(1, mino), lda, nproc, rank_mpi, comm_mpi)
        end if
        !     end added III
        do i = 1, n
            molecorb(ipc*(i - 1) + 1:ipc*i, 1:n) = molecorb(ipc*(i - 1) + 1:ipc*i, 1:n)*eigmat(i)
        end do
        if (ipc .eq. 2) then
            call zgemm_my('N', 'N', n, n, n, zone, umat, n               &
                    &, molecorb, lda, zzero, overs, lda, nproc, rank_mpi, comm_mpi)
        else
            call dgemm_my('N', 'N', n, n, n, 1.d0, umat, n               &
                    &, molecorb, lda, 0.d0, overs, lda, nproc, rank_mpi, comm_mpi)
        end if
    else
        !  added
        overs(1:ipc*n, n + 1:2*n) = molecorb(1:ipc*n, n + 1:2*n)
        !end
        if (ipc .eq. 2) then
            ! added
            call zgemm_my('N', 'N', n, n, n, zone, mat, n               &
                    &, overs(1, n + 1), lda, zzero, molecorb(1, n + 1), lda, nproc, rank_mpi, comm_mpi)
            !end
            do i = 1, n
                molecorb(2*i - 1:2*i, n + 1:2*n) = molecorb(2*i - 1:2*i, n + 1:2*n)*eigmat(i)
            end do
            call zgemm_my('N', 'N', n, n, n, zone, umat, n               &
                    &, molecorb(1, n + 1), lda, zzero, overs(1, n + 1), lda, nproc, rank_mpi, comm_mpi)
        else

            !    added
            call dgemm_my('N', 'N', n, n, n, 1.d0, mat, n               &
                    &, overs(1, n + 1), lda, 0.d0, molecorb(1, n + 1), lda, nproc, rank_mpi, comm_mpi)
            !    end added
            do i = 1, n
                molecorb(i, n + 1:2*n) = molecorb(i, n + 1:2*n)*eigmat(i)
            end do
            call dgemm_my('N', 'N', n, n, n, 1.d0, umat, n               &
                    &, molecorb(1, n + 1), lda, 0.d0, overs(1, n + 1), lda, nproc, rank_mpi, comm_mpi)
        end if
    end if

    !     We need the complex conjugate of the right eigenvectors
    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then
        do i = 1, n
            umatright(1:n*ipc, i) = umatright(1:n*ipc, i)*eigmr(i)
        end do
        if (ipc .eq. 2) then
            call zgemm_my('N', 'T', n, n, n, zone, umatright, n&
                    &, molecorb, lda, zzero, overs, lda, nproc, rank_mpi, comm_mpi)
            call conjmat(n, n, overs, lda) ! the down orbitals have to be conjugated
        else
            call dgemm_my('N', 'T', n, n, n, 1.d0, umatright, n&
                    &, molecorb, lda, 0.d0, overs, lda, nproc, rank_mpi, comm_mpi)
        end if
    end if

    if (ipf .eq. 2) then
        nshift = mod(n, 2)
    else
        nshift = 0
    end if

    !         now sort according to the norm of eig(i)
    do i = 1, n/ipf
        work(i) = abs(eig(i + nshift))
    end do

    call dsortx(work, 1, n/ipf, iwork)

    eig_sav = eig

    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        if (ipf .eq. 2) then
            do i = 1, nelorb_c/2
                eig(i) = eig_sav(iwork(i) + nshift)
                molecorb(:, 2*i - 1 + nshift) = overs(:, 2*iwork(i) - 1 + nshift)
                molecorb(:, 2*i + nshift) = overs(:, 2*iwork(i) + nshift)
            end do
            if (nshift .eq. 1) molecorb(:, 1) = overs(:, 1)
        else
            do i = 1, nelorb_c
                eig(i) = eig_sav(iwork(i))
                molecorb(:, i) = overs(:, iwork(i))
            end do
        end if
    else

        if (yes_constrainm) then
            do i = 1, nelorb_c
                eig(i) = eig_sav(iwork(i))
                molecorb(:, 2*i - 1) = overs(:, iwork(i))
                do j = 1, nelorb_c
                    molecorb(2*j - 1, 2*i) = overs(2*j - 1, iwork(i))
                    molecorb(2*j, 2*i) = -overs(2*j, iwork(i))
                end do
            end do
        else
            do i = 1, nelorb_c
                eig(i) = eig_sav(iwork(i))
                !        The even are the left eigenvectors spin up
                molecorb(:, 2*i) = overs(:, iwork(i) + n)
                !        The odd  are the right eigenvectors spin down
                molecorb(:, 2*i - 1) = overs(:, iwork(i))
            end do
        end if

    end if

    if (ipf .eq. 2) then
        !        choose a gauge (the gauge is already fixed in molec_pfaff)
        !         do i=1,n/2
        !         orbmax=sum(molecorb(1:ipc*nelorb_c,2*i-1))
        !         if(orbmax.lt.0) then
        !            molecorb(:,2*i-1)=-molecorb(:,2*i-1)
        !            molecorb(:,2*i)=-molecorb(:,2*i)
        !         endif
        !         enddo
        !    change even with odd for consistency of the sign
        overs(:, 1:n) = molecorb(:, 1:n)

        do i = 1, n/2
            molecorb(:, 2*i - 1 + nshift) = overs(:, 2*i + nshift)
            molecorb(:, 2*i + nshift) = overs(:, 2*i - 1 + nshift)
        end do

    elseif (symmagp .and. ipc .eq. 1) then
        do i = 1, nelorb_c
            orbmax = sum(molecorb(1:nelorb_c, i))
            !           do j=2,nelorb_c
            !           if(abs(molecorb(j,i)).gt.abs(orbmax)) orbmax=molecorb(j,i)
            !           enddo
            if (orbmax .lt. 0) molecorb(:, i) = -molecorb(:, i)
        end do
    else
        do i = 1, n
            orbmax = sum(molecorb(1:n*ipc, 2*i - 1))
            if (orbmax .lt. 0) then
                molecorb(:, 2*i - 1) = -molecorb(:, 2*i - 1)
                molecorb(:, 2*i) = -molecorb(:, 2*i)
            end if
        end do
    end if

    !         check orthogonality
    !        if(ipc.eq.2) then
    !        call zgemm_my('N','N',n,n,n,zone,over_sav,n,molecorb,lda,zzero,overs,lda,nproc,rank)
    !        else
    !        call dgemm_my('N','N',n,n,n,1.d0,over_sav,n,molecorb,lda,0.d0,overs,lda,nproc,rank)
    !        endif
    !        write(6,*) ' Check orthogonality '
    !        do i=mino,n
    !          do j=i,n
    !          if(ipc.eq.2) then
    !          write(6,*) i,j,zdotc_(n,molecorb(1,i),1,overs(1,j),1)
    !          else
    !          write(6,*) i,j,ddot(n,molecorb(1,i),1,overs(1,j),1)
    !          endif
    !          enddo
    !        enddo

    if (rank .eq. 0) then
        write (6, *) ' Eigenvalues hamiltonian '
        do i = 1, nelorb_c/ipf
            write (6, *) i, eig(i)
        end do
    end if
    deallocate (mat, umat, eig_sav, work, eigmat, iwork, psint)
    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) deallocate (umatright, eigmr, eigright)

#ifdef  PARALLEL
!        bcast just the relevant info to the nodes
    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        dimorb = lda*(nelorb_c - 1) + nelorb_c
    else
        dimorb = ipc*(lda*(2*nelorb_c - 1) + nelorb_c)
    end if
    call bcast_real(molecorb, dimorb, 0, comm_mpi)
    if (ndiff .gt. 0) then
        dimorb = ndiff*nelorb_c*ipc
        call bcast_real(psi_unp, dimorb, 0, comm_mpi)
    end if
    call bcast_real(eig, nelorb_c, 0, comm_mpi)
#endif
    return
end
