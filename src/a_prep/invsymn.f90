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

subroutine invsymn(n, a, lda, info, epsr, mine, rank, epssr, bands, nlax)

#ifdef __SCALAPACK
    use allio, only: ortho_cntx, me_blacs, np_ortho, me_ortho, ortho_comm, ortho_comm_id
    use descriptors
    use dspev_module
#endif
    implicit none

#ifdef PARALLEL
    include 'mpif.h'
#endif
    !        The input  a is assumed to be a positive definite matrix.
    !        The output a is the inverse of this matrix, without singular
    !        direction with too small eigenvalues, so as to have
    !        condition number > epsr.
    !        Cholesky is not used for the inverse since for matrices with
    !        very small condition number slightly negative eigenvalues can appear
    !        and Choleski crash.
    integer n, lda, info, i, j, lwork, mine, neig, rank, ierr, neginfo, bands, nlax
    !        integer iwork(5*n),ifail(n)
    !        real*8 a(lda,n),work(n*8),eig(n),eps,b(n,n),umat(n,n)
    !    1,eigmat(n),premat(n),abstol,dlamch
    real*8 eps, epsr, abstol, dlamch, epssr, condnumber, cost
    real*8, dimension(:), allocatable :: work, eig, eigmat, premat
    integer, dimension(:), allocatable :: iwork, ifail
    real*8, allocatable :: umatl(:, :)
#ifdef __SCALAPACK
    real*8, allocatable :: diag(:, :), vv(:, :), overlap(:, :)
    real*8, allocatable :: overs(:, :)
    real*8 a(nlax, nlax)
#else
    real*8 a(lda, n)
    real*8, dimension(:, :), allocatable :: b, umat
    allocate (b(n, n), work(30*n), eig(n), iwork(5*n), ifail(n))
    allocate (umat(n, n), premat(n), eigmat(n))
    b = 0.d0
    work = 0.d0
    eig = 0.d0
    iwork = 0
    ifail = 0
    umat = 0.d0
    premat = 0.d0
    eigmat = 0.d0
#endif

#ifdef __SCALAPACK
    integer :: desch(20)
    integer :: desc(descla_siz_)
    integer :: LIWORK, ic, ir, mm, nz, ii, jj, ip, jp, nrlx
    real*8 :: rtmp(10)
    integer :: itmp(10)
    real*8 :: PDLAMCH
    integer :: INDXG2L, INDXG2P, INDXL2G
#endif

#ifdef __SCALAPACK
    allocate (eig(n), eigmat(n), premat(n))
    mine = 1
    eig = 0.d0
    eigmat = 0.d0
    premat = 1.d0

    call descla_init(desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id)

    if (desc(lambda_node_) > 0) then

        allocate (umatl(desc(nlax_), desc(nlax_)))
        allocate (overlap(desc(nlax_), desc(nlax_)))
        allocate (overs(desc(nlax_), desc(nlax_)))
        umatl = 0.d0
        overlap = 0.d0
        overs = 0.d0

        ir = desc(ilar_)
        ic = desc(ilac_)

        ABSTOL = 2.0d0*PDLAMCH(ortho_cntx, 'U')

        eps = epsr

        if (eps .eq. 0.d0) eps = abstol

        call descinit(desch, n, n, desc(nlax_), desc(nlax_), 0, 0, ortho_cntx, size(umatl, 1), info)

        if (info /= 0) call errore(' cdiaghg ', ' desckinit ', abs(info))

        do j = 1, desc(nlac_)
            do i = 1, desc(nlar_)
                overs(i, j) = a(i, j)
                if (epssr .gt. 0.d0 .and. (i + ir .eq. j + ic)) overs(i, j) = overs(i, j)*(1.d0 + epssr)
            end do
        end do

        nrlx = desc(la_nrlx_)

        allocate (diag(nrlx, n), vv(nrlx, n))
        diag = 0.d0
        vv = 0.d0
        !
        call blk2cyc_redist(n, diag, nrlx, overs, size(overs, 1), desc(1))
        !
        call pdspev_drv_ss('V', diag, nrlx, eig, vv, nrlx, desc(la_nrl_), n, &
      &desc(la_npr_)*desc(la_npc_), desc(la_me_), desc(la_comm_), eps, 0)
        !
        !  Redistribute matrix "vv" into "s"
        !  matrix "s" is block distributed
        !  across 2D processors grid ( ortho_comm )
        !
        call cyc2blk_redist(n, vv, nrlx, umatl, size(umatl, 1), desc(1))
        !
        deallocate (diag, vv)

!       IF( info /= 0 ) CALL errore( ' D 02 ', ' PDSYEVD ', info )

        if (rank .eq. 0) write (6, *) ' Eigenvalues inside  inverse '
        mine = 1
        do i = 1, n
            if (eig(i)/eig(n) .gt. eps) then ! the condition number criterium
                if (rank .eq. 0) write (6, *) i, eig(i)
                eigmat(i) = dsqrt(1.d0/eig(i))
            else
                mine = mine + 1
                eigmat(i) = 0.d0
            end if
        end do

        overs = 0.d0
        do j = 1, desc(nlac_)
            if (eigmat(j + ic - 1) .ne. 0.d0) then
                do i = 1, desc(nlar_)
                    overs(i, j) = umatl(i, j)*premat(i + ir - 1)*eigmat(j + ic - 1)
                end do
            end if
        end do

!        it should be authomatically symmetric
        call pdgemm('N', 'T', n, n, n, 1.d0, overs, 1, 1, desch, overs, 1, 1, desch, 0.d0, a, 1, 1, desch)

        if (mine .ne. 1 .and. rank .eq. 0) then
            write (6, *) ' Warning neglecting', mine - 1, 'directions in SDV'
        end if
        if (rank .eq. 0) then
            condnumber = 1.d0
            do i = 1, n
                if (eigmat(i) .ne. 0.d0) then
                    cost = abs(eig(i)/eig(n))
                    if (cost .lt. condnumber) condnumber = cost
                end if
            end do
            write (6, *) ' Condition number basis set ', condnumber
        end if

        deallocate (overs)
        deallocate (umatl)
        deallocate (overlap)

    end if

    deallocate (eig, eigmat, premat)

#else

    b = 0.d0
    work = 0.d0
    umat = 0.d0
    premat = 1.d0
    !         There is numerical instability to use premat ne 1
    !         the orthogonality of the eigenvectors
    !         with respect to the metric given by overs  is poorly satisfied
    !         if premat ne 1 is used.
    eigmat = 0.d0
    eig = 0.d0
    iwork = 0
    ifail = 0

    !       Now the inverse of this matrix using SDV diagonalization
    !       by eliminating the irrelevant directions
    lwork = 30*n

    eps = abs(epsr)
    ! maximum accuracy
    abstol = 2.d0*dlamch('S')

    if (eps .eq. 0.d0) eps = abstol

    do i = 1, n
        do j = 1, n
            a(i, j) = premat(i)*premat(j)*a(i, j)
        end do
    end do
    if (epssr .gt. 0.d0) then
        do i = 1, n
            a(i, i) = a(i, i)*(1.d0 + epssr)
        end do
    end if

    call dsyevx('V', 'A', 'L', n, a, lda, 0.d0, 0.d0, 1, 1, abstol&
            &, neig, eig, umat, n, work, lwork, iwork, ifail, info)

    if (rank .eq. 0) write (6, *) ' Eigenvalues inside  inverse '
    mine = 1
    do i = 1, n
        if (eig(i)/eig(n) .gt. eps) then ! the condition number criterium
            if (rank .eq. 0) write (6, *) i, eig(i)
            eigmat(i) = dsqrt(1.d0/eig(i))
        else
            mine = mine + 1
            eigmat(i) = 0.d0
        end if
    end do
    if (info .ne. 0 .and. rank .eq. 0)                                    &
            & write (6, *) ' info > 0 in dsyevx !!! ', info

    do i = 1, info
        if (eigmat(ifail(i)) .ne. 0.d0) then
            mine = mine + 1
            eigmat(ifail(i)) = 0.d0
        end if
    end do

    !       positive eigenvalues
    do j = 1, n
        if (eigmat(j) .ne. 0.d0) then
            do i = 1, n
                b(i, j) = premat(i)*umat(i, j)*eigmat(j)
            end do
        else
            b(:, j) = 0.d0
        end if
    end do

    !        it should be authomatically symmetric
    call dgemm('N', 'T', n, n, n, 1.d0, b, n, b, n, 0.d0, a, lda)

    if (mine .ne. 1 .and. rank .eq. 0) write (6, *) ' Warning neglecting', mine - 1, 'directions in SDV'
    if (rank .eq. 0) then
        condnumber = 1.d0
        do i = 1, n
            if (eigmat(i) .ne. 0.d0) then
                cost = abs(eig(i)/eig(n))
                if (cost .lt. condnumber) condnumber = cost
            end if
        end do
        write (6, *) ' Condition number basis set ', condnumber
    end if

    deallocate (b, work, eig, iwork, ifail, umat, premat, eigmat)
#endif

    return
end
