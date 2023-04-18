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

subroutine eval_molec_epsdgel(nelorb_c, overs, mat_in               &
        &, molecorb, eig, lda, eps, nproc, rank, rank_mpi, comm_mpi, optprint, symmagp)
    use constants, only: ipc, ipf, zone, zzero
    use allio, only: molecular, yes_hermite, symmetrize_agp&
            &, opposite_phase, same_phase, real_contracted, gauge_fixing, printoverlap
    implicit none
    integer nelorb_c, i, j, k, rank, info, n, lda, mine, dimo, neig, lwork&
            &, optprint, dimorb, ierr, il, iu, indi, indj, nproc, old_threads, comm_mpi&
            &, rank_mpi, nshift
    real*8 overs(ipc*lda, *), mat_in(ipc*lda, *)                       &
            &, molecorb(ipc*lda, *), eig(nelorb_c)
    real*8 abstol, eps, orbmax, dlamch, vl, vu, cost, maxerr, timep
    real*8, external :: cclock
    real*8, dimension(:, :), allocatable :: mat, umat, umatleft
    real*8, dimension(:), allocatable :: eig_sav, work, eigmat, eigleft, eigml
    integer, dimension(:), allocatable :: iwork
    integer, external :: omp_get_max_threads
    logical symmagp, eqover, yes_constrainm, yesh
#ifdef PARALLEL
    include 'mpif.h'
#ifdef __TEST
    integer dimmat, dimmat2
    if (nproc .gt. 1) then
        dimmat = ipc*((nelorb_c - 1)*lda + nelorb_c)
        if (symmagp .and. ipc .eq. 1) then
            allocate (mat(ipc*lda, nelorb_c))
            mat(1:ipc*nelorb_c, 1:nelorb_c) = overs(1:ipc*nelorb_c, 1:nelorb_c)
            call bcast_real(overs, dimmat, 0, comm_mpi)
            if (sum(abs(mat(1:ipc*nelorb_c, 1:nelorb_c) - overs(1:ipc*nelorb_c, 1:nelorb_c))) .ne. 0.d0)&
              & write (6, *) ' Inconsistent input ', rank

        else
            allocate (mat(ipc*lda, 2*nelorb_c))
            mat(1:ipc*nelorb_c, 1:2*nelorb_c) = overs(1:ipc*nelorb_c, 1:2*nelorb_c)
            dimmat2 = ((2*nelorb_c - 1)*lda + nelorb_c)*ipc
            call bcast_real(overs, dimmat2, 0, comm_mpi)
            if (sum(abs(mat(1:ipc*nelorb_c, 1:2*nelorb_c) - overs(1:ipc*nelorb_c, 1:2*nelorb_c)))&
           &.ne. 0.d0) write (6, *) ' Inconsistent input ', rank
        end if
        mat(1:ipc*nelorb_c, 1:nelorb_c) = mat_in(1:ipc*nelorb_c, 1:nelorb_c)
        call bcast_real(mat_in, dimmat, 0, comm_mpi)
        if (sum(abs(mat(1:ipc*nelorb_c, 1:nelorb_c) - mat_in(1:ipc*nelorb_c, 1:nelorb_c))) .ne. 0.d0)&
           & write (6, *) ' Inconsistent input ', rank
        deallocate (mat)
    end if
#endif
#endif

    !       if symmagp=.false. molecorb is defined ldax 2 x nelorb_c
    !       the first nelorb_c are the right eigenvectors of the SDV and the second
    !       ones are the left.

    vl = 0.d0
    vu = 0.d0
    il = 0
    iu = 0

    !      if(rank_mpi.eq.0.or.nproc.gt.1) then

    if (ipc .eq. 2 .and. yes_hermite .and. symmagp .and. symmetrize_agp .and. ipf .eq. 1) then
        yes_constrainm = .true.
    else
        yes_constrainm = .false.
    end if

    lwork = 5*nelorb_c ! memory required by dgesvd

    allocate (mat(ipc*nelorb_c, nelorb_c), umat(ipc*nelorb_c, nelorb_c), &
            &eig_sav(nelorb_c), work(lwork)&
            &, eigmat(nelorb_c), iwork(5*nelorb_c))
    n = nelorb_c

    !    write(6,*) ' Input matrix ',sum(abs(mat_in(1:ipc*n,1:n))),sum(mat_in(1:ipc*n,1:n))
    !        do j=1,n
    !           do k=1,n
    !           write(6,*) j,k,mat_in(ipc*(j-1)+1:ipc*j,k)
    !           enddo
    !        enddo

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

        allocate (umatleft(ipc*nelorb_c, nelorb_c))
        allocate (eigleft(nelorb_c), eigml(nelorb_c))
        !        ARROCCO input
        ! OK
        !        Here we have the spin up left in overs(:,n+1:2*n), the right are in overs(:,1:n)
        do i = 1, n
            umat(1:n*ipc, i) = overs(1:n*ipc, i)
            overs(1:n*ipc, i) = overs(1:n*ipc, i + n)
            overs(1:n*ipc, n + i) = umat(1:n*ipc, i)
        end do
    end if

    ! destroy invo
    !        calculation overlap old
    !        calculation normalized orbitals and eigenvectors
    ! maximum accuracy
    abstol = 0.d0
    if (eps .eq. 0.d0) eps = 10.d0*dlamch('E')

    !    This code compute the generalized eigenvalue equation:
    !      mat_in x over_s  y = lambda y
    !      by stable preconditioning and pre-diagonalization of the matrix over_s

    !        put to zero all matrix elements not contained in mat_in
    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        if (ipf .eq. 2) then
            indj = n
            indi = n
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
            !!       CCC bisogna eliminare
            eig_sav = 0.d0
            do i = 1, n
                do j = i, n
                    if (mat_in(i, j) .ne. 0.d0) then
                        eig_sav(i) = 1.d0
                        eig_sav(j) = 1.d0
                    end if
                end do
            end do
            indi = 0
            indj = 0
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
        end if

    else

        !       Determine the non zero columns to be excluded as the matrix has no connection with them
        eig_sav = 0.d0
        do i = 1, n
            if (ipc .eq. 1) then
                do j = 1, n
                    if (mat_in(j, i) .ne. 0.d0) then
                        eig_sav(i) = 1.d0
                    end if
                end do
            else
                do j = 1, n
                    if (mat_in(2*j - 1, i) .ne. 0.d0 .or. mat_in(2*j, i) .ne. 0.d0) then
                        eig_sav(i) = 1.d0
                    end if
                end do
            end if
        end do
        indi = 0
        indj = 0
        do i = 1, n
            if (eig_sav(i) .eq. 1.d0) then
                indi = indi + 1
                indj = 0
                if (ipc .eq. 1) then
                    do j = 1, n
                        if (eig_sav(j) .eq. 1.d0) then
                            indj = indj + 1
                            overs(indi, indj) = overs(i, j)
                        end if
                    end do
                else
                    do j = 1, n
                        if (eig_sav(j) .eq. 1.d0) then
                            indj = indj + 1
                            overs(2*indi - 1, indj) = overs(2*i - 1, j)
                            overs(2*indi, indj) = -overs(2*i, j) ! the right have the complex conjugate
                        end if
                    end do
                end if
            end if
        end do
    end if

    umat = 0.d0
    !  Build right orthogonal basis for the spin down, put in umat the unitary transform

    if (ipc .eq. 1) then

        if (molecular .eq. 0 .and. ipf .eq. 2) then
            call DSYEV_MY('V', 'L', indi/2, overs, LDA, eig, INFO&
                    &, nproc, rank_mpi, comm_mpi)
            call DSYEV_MY('V', 'L', indi/2, overs(indi/2 + 1, indi/2 + 1), LDA&
                    &, eig(indi/2 + 1), INFO, nproc, rank_mpi, comm_mpi)
        else
            call DSYEV_MY('V', 'L', indi, overs, LDA, eig, INFO&
                    &, nproc, rank_mpi, comm_mpi)
        end if

    else
        if (molecular .eq. 0 .and. ipf .eq. 2) then
            call ZSYEV_MY('V', 'L', indi/2, overs, LDA, eig, INFO&
                    &, nproc, rank_mpi, comm_mpi)
            call ZSYEV_MY('V', 'L', indi/2, overs(indi + 1, indi/2 + 1), LDA&
                    &, eig(indi/2 + 1), INFO, nproc, rank_mpi, comm_mpi)
        else
            call ZSYEV_MY('V', 'L', indi, overs, LDA, eig, INFO&
                    &, nproc, rank_mpi, comm_mpi)
        end if
    end if
#ifdef DEBUG
    if (ipc .eq. 1) then
        call dgemm_my('T', 'N', indi, indi, indi, 1.d0, overs, LDA, overs, LDA, 0.d0, umat&
       &, nelorb_c, nproc, rank_mpi, comm_mpi)
    else
        call zgemm_my('C', 'N', indi, indi, indi, zone, overs, LDA, overs, LDA, zzero, umat&
       &, nelorb_c, nproc, rank_mpi, comm_mpi)
    end if
    maxerr = 0.d0
    do i = 1, indi
        do j = i + 1, indi
            cost = sum(abs(umat(ipc*(i - 1) + 1:ipc*i, j)))
            if (cost .gt. maxerr) maxerr = cost
        end do
    end do
    if (rank .eq. 0 .and. maxerr .gt. 1.d-6) &
   &write (6, *) ' Warning precision orthogonality eigenvectors =', maxerr
    umat = 0.d0
#endif

    umat(1:ipc*indi, 1:indi) = overs(1:ipc*indi, 1:indi)

    ! Reorder the considered taken directions as the largest ones
    if (indj .gt. 0) then
        work(1:indj) = eig(1:indj)
        eig = -1.d0
        eig(n - indj + 1:n) = work(1:indj)
    end if
    do j = indj, 1, -1
        work(1:ipc*indj) = umat(1:ipc*indj, j)
        umat(:, n - indj + j) = 0.d0
        indi = 0
        do i = 1, n
            if (eig_sav(i) .eq. 1.d0) then
                indi = indi + 1
                if (ipc .eq. 1) then
                    umat(i, n - indj + j) = work(indi)
                else
                    umat(2*i - 1, n - indj + j) = work(2*indi - 1)
                    umat(2*i, n - indj + j) = work(2*indi)
                end if
            end if
        end do
    end do

    do i = 1, n - indj
        umat(1:ipc*n, i) = 0.d0
    end do

    mine = 1

    if (rank .eq. 0 .and. optprint .ne. 0) then
        write (6, *) ' min/max Eigenvalues overlap matrix ', indj, eig(min(n - indj + 1, n)), eig(n)
        !        Extended output
        !        do i=1,indj
        !        write(6,*) i,eig(n-indj+i)
        !        enddo
    end if
    do i = 1, n
        if (eig(i)/eig(n) .gt. eps) then ! the condition number criterium
            eigmat(i) = dsqrt(1.d0/eig(i))
        else
            mine = i + 1
            !        write(6,*) ' warning zero eigenvalue !!! '
            if (ipf .eq. 2) then
                eigmat(i) = dsqrt(1.d0/(eps*eig(n)))
            else
                eigmat(i) = 0.d0
            end if
        end if
    end do
    if (info .ne. 0 .and. rank .eq. 0) write (6, *)                         &
            &' info > 0 in dsyevx !!! ', info

    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then

        if (eqover) then
            !        umatleft here is indeed the unitary matrix corresponding to right basis for the spin down
            umatleft = umat
            eigleft = eig
            eigml = eigmat
        else

            work(1:n) = eig_sav(1:n)
            eig_sav = 0.d0
            !        Determine the non zero raws for the left eigenvectors
            do i = 1, n
                do j = 1, n
                    if (sum(abs(mat_in(ipc*(i - 1) + 1:ipc*i, j))) .ne. 0.d0) then
                        eig_sav(i) = 1.d0
                    end if
                end do
            end do
            indi = 0
            indj = 0
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

            if ((.not. opposite_phase .and. .not. same_phase) .or. .not. real_contracted) then
                if (ipc .eq. 1) then
                    call DSYEV_MY('V', 'L', indi, overs(1, n + 1), LDA, eigleft, INFO&
                            &, nproc, rank_mpi, comm_mpi)
                else
                    call ZSYEV_MY('V', 'L', indi, overs(1, n + 1), LDA, eigleft, INFO&
                            &, nproc, rank_mpi, comm_mpi)
                end if
                umatleft(1:ipc*indi, 1:indi) = overs(1:ipc*indi, n + 1:n + indi)

                !   Reorder the eigenvectors

                if (indj .gt. 0) then
                    work(1:indj) = eigleft(1:indj)
                    eigleft = -1.d0
                    eigleft(n - indj + 1:n) = work(1:indj)
                else
                    eigleft = -1.d0
                end if
                do j = indj, 1, -1
                    work(1:ipc*indj) = umatleft(1:ipc*indj, j)
                    umatleft(:, n - indj + j) = 0.d0
                    indi = 0
                    do i = 1, n
                        if (eig_sav(i) .eq. 1.d0) then
                            indi = indi + 1
                            if (ipc .eq. 1) then
                                umatleft(i, n - indj + j) = work(indi)
                            else
                                umatleft(2*i - 1, n - indj + j) = work(2*indi - 1)
                                umatleft(2*i, n - indj + j) = work(2*indi)
                            end if
                        end if
                    end do
                end do

                do i = 1, n - indj
                    umatleft(1:ipc*n, i) = 0.d0
                end do

                if (rank .eq. 0 .and. optprint .ne. 0) then
                    write (6, *) 'min/max Eigenvalues overlap matrix left eigenvectors ', &
                        indj, eigleft(min(n - indj + 1, n)), eigleft(n)
                    !     Extended output
                    !     do i=1,indj
                    !     write(6,*) i,eigleft(n-indj+i)
                    !     enddo
                end if

                do i = 1, n
                    if (eigleft(i)/eigleft(n) .gt. eps) then ! the condition number criterium
                        eigml(i) = dsqrt(1.d0/eigleft(i))
                    else
                        !        write(6,*) ' warning zero eigenvalue !!! '
                        eigml(i) = 0.d0
                    end if
                end do
                if (info .ne. 0 .and. rank .eq. 0) write (6, *)                         &
                        &' info > 0 in dsyevx !!! ', info
            else ! real_contracted
                if (opposite_phase .or. ipc .eq. 1) then
                    umatleft = umat
                else
                    do i = 1, n
                        do j = 1, n
                            umatleft(2*j - 1, i) = umat(2*j - 1, i)
                            umatleft(2*j, i) = -umat(2*j, i)
                        end do
                    end do
                end if
                eigleft = eig
                eigml = eigmat
            end if
        end if ! endif eqover
    end if ! endif not symmagp

    if (mine .ne. 1 .and. rank .eq. 0 .and. optprint .ne. 0)&
            &write (6, *) ' disregarded coll. =', mine - 1

    dimo = n - mine + 1

    !        Now modifying the matrix detmat_c (and forgetting it)
    !        destroy overs mat_in
    if (ipc .eq. 1) then
        call dgemm_my('N', 'N', n, n, n, 1.d0, mat_in, lda, umat, n, 0.d0, overs, lda&
                &, nproc, rank_mpi, comm_mpi)
        if (symmagp .or. ipf .eq. 2) then
            call dgemm_my('T', 'N', n, n, n, 1.d0, umat, n, overs, lda, 0.d0, mat_in, lda&
                    &, nproc, rank_mpi, comm_mpi)
        else
            call dgemm_my('T', 'N', n, n, n, 1.d0, umatleft, n, overs, lda, 0.d0, mat_in&
                    &, lda, nproc, rank_mpi, comm_mpi)
        end if
    else
        if (ipf .eq. 2) then
            call conjmat(n, n, umat, n)
            call zgemm_my('N', 'N', n, n, n, zone, mat_in, lda, umat, n, zzero, overs, lda&
                    &, nproc, rank_mpi, comm_mpi)
            call conjmat(n, n, umat, n)
            call zgemm_my('C', 'N', n, n, n, zone, umat, n, overs, lda, zzero, mat_in, lda&
                    &, nproc, rank_mpi, comm_mpi)
        else
            call zgemm_my('N', 'N', n, n, n, zone, mat_in, lda, umat, n, zzero, overs, lda&
                    &, nproc, rank_mpi, comm_mpi)
            call zgemm_my('C', 'N', n, n, n, zone, umatleft, n, overs, lda, zzero, mat_in&
                    &, lda, nproc, rank_mpi, comm_mpi)
        end if
    end if

    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        do i = 1, n
            do j = 1, n
                if (eigmat(i)*eigmat(j) .ne. 0.d0) then
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = mat_in(ipc*(i - 1) + 1:ipc*i, j)/eigmat(i)/eigmat(j)
                else
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = 0.d0
                end if
            end do
        end do
    else
        do i = 1, n
            do j = 1, n
                if (eigml(i)*eigmat(j) .ne. 0.d0) then
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = mat_in(ipc*(i - 1) + 1:ipc*i, j)/eigml(i)/eigmat(j)
                else
                    mat_in(ipc*(i - 1) + 1:ipc*i, j) = 0.d0
                end if
            end do
        end do
    end if

    if (symmagp .and. ipc .eq. 1 .or. ipf .eq. 2) then
        if (mine .gt. 1) eig(1:mine - 1) = 0.d0
        molecorb(1:lda, 1:n) = 0.d0
        !      if(rank.eq.0) write(6,*) ' dimo mine inside main ',dimo,mine,eig(1:mine-1)
        if (ipf .eq. 2) then
            !      CCC tua routine

            !          do i=1, nelorb_c
            !             do j=1, nelorb_c
            !         if (abs(mat_in(ipc*(i-1)+1,j)+mat_in(ipc*(i-1)+2,j)+mat_in(ipc*(j-1)+1,i)+mat_in(ipc*(j-1)+2,i)).gt.1d-8) &
            !                                        write (6,*) i,j,mat_in(ipc*(i-1)+1,j), mat_in(ipc*(i-1)+2,j)
            !                   if (ipc.eq.1) then
            !                      if (abs(mat_in(i,j)+mat_in(j,i)).gt.1d-8) write (6,*) i,j, mat_in(i,j)
            !                   end if
            !                end do
            !             end do

            call pfaffian_mo(LDA, nelorb_c, ipc, mat_in, eig, molecorb)
            !         if (rank.eq.0) write (6,*) "eig", eig
            !         if (rank.eq.0) write (6,*) "Stopping after pfaffian mo"

            if (rank .eq. 0 .and. printoverlap) then
                write (6, *) '  Normalization up/down molecular orbitals '
                if (mod(nelorb_c, 2) .eq. 0) then
                    do i = 1, nelorb_c
                        write (6, *) i, eig((i + 1)/2), sum(molecorb(1:ipc*nelorb_c/2, i)**2)&
                                &, sum(molecorb(ipc*nelorb_c/2 + 1:ipc*nelorb_c, i)**2)
                    end do
                else
                    do i = 1, nelorb_c
                        write (6, *) i, eig(i/2 + 1), sum(molecorb(1:ipc*nelorb_c/2, i)**2)&
                                &, sum(molecorb(ipc*nelorb_c/2 + 1:ipc*nelorb_c, i)**2)
                    end do
                end if
            end if

            !call mpi_finalize(ierr)

        else
            if (dimo .gt. 0) then
                call DSYEV_MY('V', 'L', dimo, mat_in(mine, mine), LDA, eig(mine), INFO&
                        &, nproc, rank_mpi, comm_mpi)
                molecorb(mine:mine + dimo - 1, mine:mine + dimo - 1) = &
                        & mat_in(mine:mine + dimo - 1, mine:mine + dimo - 1)
            end if
        end if

    else
        timep = cclock()
        if ((nproc .gt. 1 .or. ipc .eq. 2) .and. ipf .eq. 1) then
            molecorb(1:ipc*n, n + 1:2*n) = mat_in(1:ipc*n, 1:n)
            if (ipc .eq. 1) then
                call DGESVD_MY('V', N, N, molecorb(1, n + 1), LDA, molecorb, eig, nproc, rank_mpi, comm_mpi)
            else

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

            end if ! ipc
        else ! if nproc>1 or ipc = 2
            call dgesvd('A', 'A', n, n, mat_in, lda&
                    &, eig, molecorb(1, n + 1), lda, molecorb, lda, work, lwork, info)
        end if ! if nproc>1 or ipc = 2
        if (ipc .eq. 1 .and. gauge_fixing) then
            do i = 1, N
                call gauge_fixr(N, molecorb(1, n + i), 1, molecorb(i, 1), LDA)
            end do
        end if
        timep = cclock() - timep
        if (rank .eq. 0) write (6, *) ' Time SVD =', timep

    end if ! endif ipc.eq.1.and.symmagp

    !CCC Finisce l'if
    !
    !      now go back to the original representation changing umat
    do i = 1, n
        do j = 1, ipc*n
            umat(j, i) = umat(j, i)*eigmat(i)
        end do
    end do

    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then

        do i = 1, n
            do j = 1, ipc*n
                umatleft(j, i) = umatleft(j, i)*eigml(i)
            end do
        end do
        !   The right (down) have to be complex conjugated at the end
        if (ipc .eq. 2) then
            call zgemm_my('N', 'T', n, n, n, zone, umat, n, molecorb, lda, zzero&
                    &, overs, lda, nproc, rank_mpi, comm_mpi)
        else
            call dgemm_my('N', 'T', n, n, n, 1.d0, umat, n, molecorb, lda, 0.d0&
                    &, overs, lda, nproc, rank_mpi, comm_mpi)
        end if
    else

        if (ipc .eq. 2) then
            call zgemm_my('N', 'N', n, n, n, zone, umat, n, molecorb, lda, zzero&
                    &, overs, lda, nproc, rank_mpi, comm_mpi)
        else
            call dgemm_my('N', 'N', n, n, n, 1.d0, umat, n, molecorb, lda, 0.d0&
                    &, overs, lda, nproc, rank_mpi, comm_mpi)
        end if
    end if

    !        now sort according to the norm of eig(i)

    if (ipf .eq. 2) then
        nshift = mod(n, 2)
    else
        nshift = 0
    end if

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

            !    Check identity
            !      if(ipc.eq.2) then
            !      call  check_complex(n,eig,molecorb,mat)
            !      else
            !      mat=0.d0
            !        do i=1,nelorb_c/2
            !          do j=1,nelorb_c
            !            do k=1,nelorb_c
            !            mat(j,k)=mat(j,k)-eig(i)*(molecorb(j,2*i-1)*molecorb(k,2*i)-&
            !           &molecorb(k,2*i-1)*molecorb(j,2*i))
            !            enddo
            !          enddo
            !        enddo
            !       endif

            !       write(6,*) ' Output matrix ',sum(abs(mat(:,:))),sum(mat(:,:))
            !         do j=1,n
            !            do k=1,n
            !            write(6,*) j,k,mat(ipc*(j-1):ipc*j,k)
            !            enddo
            !         enddo
            !       write(6,*) ' Final molecular inside eval '
            !       do i=1,n
            !       write(6,*) i,sum(molecorb(:,i))
            !       enddo
            !       stop

        else
            do i = 1, nelorb_c
                eig(i) = eig_sav(iwork(i))
                molecorb(:, i) = overs(:, iwork(i))
            end do
        end if
    else

        !     Below are the even
        if (ipc .eq. 2) then
            call zgemm_my('N', 'N', n, n, n, zone, umatleft, n, molecorb(1, n + 1), lda, zzero&
                    &, overs(1, n + 1), lda, nproc, rank_mpi, comm_mpi)
            !     The odd are the right eigenvectors spin down , complex conjugated at the end
            call conjmat(n, n, overs, lda)
        else
            call dgemm_my('N', 'N', n, n, n, 1.d0, umatleft, n, molecorb(1, n + 1), lda, 0.d0&
                    &, overs(1, n + 1), lda, nproc, rank_mpi, comm_mpi)
        end if
        !        write(6,*) ' Output svd '
        !        do i=1,n
        !        write(6,*) ' Eigenvalue ',i,eig(i)
        !        write(6,*) ' Sorted index =',i,iwork(i)
        !        write(6,*) ' Left  eigenvector spin-up =',overs(1:n,n+i)
        !        write(6,*) ' Right eigenvector spin-down=',overs(1:n,i)
        !        enddo
        !        The even  are the left  eigenvectors spin up
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
                molecorb(:, 2*i - 1) = overs(:, iwork(i))
                molecorb(:, 2*i) = overs(:, iwork(i) + n)
            end do
        end if
    end if

    if (ipf .eq. 2) then
        !        gauge  chosen in molec_pfaff
        !     do i=1,n/2
        !       orbmax=sum(molecorb(1:ipc*nelorb_c,2*i-1))
        !         if(orbmax.lt.0) then
        !            molecorb(:,2*i-1)=-molecorb(:,2*i-1)
        !            molecorb(:,2*i)=-molecorb(:,2*i)
        !         endif
        !      enddo
        !    change even with odd for consistency of the sign
        overs(:, 1:n) = molecorb(:, 1:n)
        do i = 1, n/2
            molecorb(:, 2*i - 1 + nshift) = overs(:, 2*i + nshift)
            molecorb(:, 2*i + nshift) = overs(:, 2*i - 1 + nshift)
        end do

    elseif (symmagp .and. ipc .eq. 1) then
        do i = 1, nelorb_c
            orbmax = sum(molecorb(1:nelorb_c, i))
            if (orbmax .lt. 0) then
                molecorb(:, i) = -molecorb(:, i)
            end if
        end do
    else
        do i = 1, nelorb_c
            orbmax = sum(molecorb(1:ipc*nelorb_c, 2*i - 1))
            if (orbmax .lt. 0) then
                molecorb(:, 2*i - 1) = -molecorb(:, 2*i - 1)
                molecorb(:, 2*i) = -molecorb(:, 2*i)
            end if
        end do
    end if
    !          if(rank.eq.0) then
    !          write(6,*) ' Output svd '
    !          do i=1,n
    !          write(6,*) ' Eigenvalue ',i,eig(i)
    !          write(6,*) ' Left  eigenvector =',sum(molecorb(1:n*ipc,2*i)**2)
    !          write(6,*) ' Right eigenvector =',sum(molecorb(1:ipc*n,2*i-1)**2)
    !          enddo
    !          endif

    !        if(info.eq.0) then
    if (optprint .ne. 0 .and. rank .eq. 0) then
        write (6, *) ' Eigenvalues Det '
        do i = 1, nelorb_c/ipf
            write (6, *) i, eig(i)
        end do
    end if

    deallocate (mat, umat, eig_sav, work, eigmat, iwork)
    if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) deallocate (umatleft, eigml, eigleft)
    !       endif ! endif rank.eq.0
#ifdef  PARALLEL
    if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
        dimorb = ipc*(lda*(nelorb_c - 1) + nelorb_c)
    else
        dimorb = ipc*(lda*(2*nelorb_c - 1) + nelorb_c)
    end if
    call bcast_real(molecorb, dimorb, 0, comm_mpi)
    call bcast_real(eig, nelorb_c, 0, comm_mpi)
#endif
    return

end subroutine eval_molec_epsdgel

subroutine gauge_fix(N, molecorbup, ldup, molecorbdown, lddown)
    implicit none
    integer n, ldup, lddown
    complex*16 molecorbup(ldup, *), molecorbdown(lddown, *), over
    over = sum(molecorbup(1, 1:n)*conjg(molecorbdown(1, 1:n)))
    over = over/abs(over)
    molecorbdown(1, 1:n) = molecorbdown(1, 1:n)*over
    return
end
subroutine gauge_fixr(N, molecorbup, ldup, molecorbdown, lddown)
    implicit none
    integer n, ldup, lddown
    real*8 molecorbup(ldup, *), molecorbdown(lddown, *), over
    over = sum(molecorbup(1, 1:n)*molecorbdown(1, 1:n))
    if (over .lt. 0) molecorbdown(1, 1:n) = -molecorbdown(1, 1:n)
    return
end

subroutine check_complex(n, eig, molecorb, mat)
    implicit none
    real*8 eig(*)
    integer n, i, j, k
    complex*16 molecorb(n, *), mat(n, n)
    mat = (0.d0, 0.d0)
    do i = 1, n/2
        do j = 1, n
            do k = 1, n
                mat(j, k) = mat(j, k) - eig(i)*(molecorb(j, 2*i - 1)*molecorb(k, 2*i) - &
                        &molecorb(k, 2*i - 1)*molecorb(j, 2*i))
            end do
        end do
    end do

    write (6, *) ' Output complex matrix '
    do j = 1, n
        do k = 1, n
            write (6, *) j, k, mat(j, k)
        end do
    end do
    return
end
