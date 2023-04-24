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

subroutine eval_hamilt(nelorb_c, oversav, matsav&
        &, molecorb, eig, lda, epsr, eps_mach, rank, optprint, lworkr, iopt&
        &, premat, eigmat, info, bands, nlax, mincond)

    !      This subroutine was written by Sandro Sorella and Carlo Cavazzoni
    !      14 November 2008. Revised 20 January 2008
    !
    !      The task of this subroutine is to compute in a stable way the first
    !      bands= (~# electrons/2)  eigenvectors of an Hamiltonian defined in
    !      a non orthogonal basis of linear dimension nelorb_c>>bands.
    !      The overlap matrix input matrix elements
    !      being oversav and the corresponding hamiltonian matrix elements matsav.
    !      These matrices are not destroyed on output.
    !      The problem is that the matrix oversav may have a very small condition
    !      number and straighforward diagonalization using scalapack  PDSYGVX
    !      leads to inaccurate and often garbage eigenvectors.
    !      In parallel computation there is also the further complication that
    !      it is difficult to preserve orthogonality (and accuracy) of eigenvectors
    !      without having a huge buffer memory   at disposal, as no memory
    !      distribution  is possible in order to orthogonalize.
    !      The parallel code is therefore slightly different from the scalar one.
    ! #########################################################################
    !      Scalar algorithm:
    !      1)    Diagonalize the full overlap matrix with maximum accuracy
    !   (ABSTOL=2.*dlamch('S')), compute eigenvalies eig(i) and eigenvectors umat).
    !       After removing the singular eigenvalues recompute and reorthoganilize
    !       them exactly as is done in the parallel algorithm.
    !      2)  Recompute the matrix elements of the Hamiltonian in the new
    !          orthogonal basis by taking only the vectors {i}  with eigenvalues
    !          eig(i)/eig_max < epsr (input). So that the new overlap matrix
    !          has by definition condition number > eps.
    !      3) Now the diagonalization of the matrix in the new basis can be done
    !         because it is a well conditioned  matrix as
    !         the eigenvalues are the physical  eigenvalues of the Hamiltonian.
    !      4)  Go back in the original basis for the eigenvectors computed.
    ! #########################################################################
    !      PARALLEL  algorithm (Using SCALAPACK version...):
    !      1) Apply Cavazzoni's diagonalization  to the overlap matrix.
    !         Disregard eigenvalues  of the overlap matrix with
    !         eig(i)/eig_max < eps.
    !         This is justified from the fact that an eigenvector of the
    !         overlap matrix with nearly zero eigenvalue corresponds  to
    !         a linear combination of states in the original basis with
    !         vanishing norm, i.e. the non orthogonal basis is redundant and
    !         this direction can be eliminated within numerical accuracy<sqrt(eps).
    !          In this step we do not require neither orthogonality between
    !          the selected eigenvectors nor maximum accuracy.
    !          However the accuracy has to be sufficient to guarantee that
    !          the eigenvectors are at least linearly independent.
    !          In a second step we diagonalize again the overlap matrix in
    !          this basis without the singular vectors corresponding to the
    !          disregarded eigenvalues. Since in this basis
    !          the matrix has very good condition number, we can assume that
    !          the Cavazzoni's diagonalization routine provides a very accurate
    !          orthonormal basis set. Then we store in umatl the transformation
    !          between the original basis and the orthogonal one.
    !          This step is done only at the inizialization (iopt=1) and umatl is
    !          saved, and used  for all further iterations (iopt=0).
    !      2)  At each iteration recompute the hamiltonian matrix H
    !          in this new basis.
    !      3)  Then  apply again the Cavazzoni's diagonalization routine.
    !          Computation and memory requirements are minimal.
    !      4)  Go back to  the original basis for the eigenvectors computed at
    !          each iteration.

#ifdef __SCALAPACK
    use allio, only: ortho_cntx, me_blacs, np_ortho, me_ortho, ortho_comm, ortho_comm_id, commrep_mpi
    use descriptors
    use dspev_module
#else
    use allio, only: rankrep, commrep_mpi
#endif
    implicit none
    integer, intent(in) :: nlax
    integer nelorb_c, i, j, rank, info, n, lda, mine, dimo, neig, lwork, optprint&
            &, lworkr, iopt, countexc, bands, dimorb, ierr, mincond
    real*8 eig(nelorb_c), premat(nelorb_c)&
            &, eigmat(*), abstol, eps, epsr, orbmax, dlamch, minover&
            &, condnumber, cost, maxeig, eps_mach
    real(8), dimension(:, :), allocatable :: mat_in, overs
    real(8), dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: iwork, ifail
#ifdef __SCALAPACK
    real*8 molecorb(nlax, nlax), oversav(nlax, nlax), matsav(nlax, nlax)
#else
    real*8 molecorb(lda, nelorb_c), oversav(lda, nelorb_c), matsav(lda, nelorb_c)
#endif

#ifdef PARALLEL
    include 'mpif.h'
#endif

#ifdef __SCALAPACK
    integer :: desch(20)
    integer :: desc(descla_siz_)
    integer :: LIWORK, ic, ir, mm, nz, ii, jj, ip, jp, nrlx
    real*8 :: PDLAMCH
    integer :: INDXG2L, INDXG2P, INDXL2G
#endif

    real*8, allocatable, save :: umatl(:, :)
    real*8, allocatable :: diag(:, :), vv(:, :), overlap(:, :)

    premat = 1.d0

#ifdef __SCALAPACK

    n = nelorb_c
    mine = 1

    call descla_init(desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id)

    if (desc(lambda_node_) > 0) then

        allocate (mat_in(desc(nlax_), desc(nlax_)))
        mat_in = 0.d0
        if (.not. allocated(umatl)) then
            allocate (umatl(desc(nlax_), desc(nlax_)))
            umatl = 0.d0
        end if
        allocate (overs(desc(nlax_), desc(nlax_)))
        overs = 0.d0

        ir = desc(ilar_)
        ic = desc(ilac_)

!    if(rank.eq.0) write(6,*) ' abstol =',abstol

        eps = epsr

        if (eps .eq. 0.d0) eps = abs(eps_mach)

!    For the first part we do not assume to have accurate eigenvectors.

        call descinit(desch, n, n, desc(nlax_), desc(nlax_), 0, 0, ortho_cntx, size(umatl, 1), info)

        if (info /= 0) call errore(' cdiaghg ', ' descinit ', abs(info))

!    do j=1,desc( nlac_ )
!       do i=1,desc( nlar_ )
!           overs(i,j)=oversav(i+ir-1,j+ic-1)
!       enddo
!    enddo
!    do j=1,desc( nlac_ )
!       do i=1,desc( nlar_ )
!           mat_in(i,j)=matsav(i+ir-1,j+ic-1)
!       enddo
!    enddo
!    To be sure there is no garbage.
        overs = oversav
        mat_in = matsav

!    call pdgeadd('N',n,n,1.d0,oversav,1,1,desch,0.d0,overs,1,1,desch )
!    call pdgeadd('N',n,n,1.d0,matsav,1,1,desch,0.d0,mat_in,1,1,desch )

        if (iopt == 1) then

            allocate (overlap(desc(nlax_), desc(nlax_)))
            overlap = 0.d0

            do j = 1, desc(nlac_)
                do i = 1, desc(nlar_)
                    overs(i, j) = overs(i, j)*premat(i + ir - 1)*premat(j + ic - 1)
                end do
            end do

            nrlx = desc(la_nrlx_)

            allocate (diag(nrlx, n), vv(nrlx, n))

            diag = 0.d0
            vv = 0.d0
            !
            call blk2cyc_redist(n, diag, nrlx, overs, size(overs, 1), desc(1))
            !
!    iopt=0, the standard one
            call pdspev_drv_ss('V', diag, nrlx, eig, vv, nrlx, desc(la_nrl_), n,&
         &desc(la_npr_)*desc(la_npc_), desc(la_me_), desc(la_comm_), eps_mach, 1)
            !
            !  Redistribute matrix "vv" into "s"
            !  matrix "s" is block distributed
            !  across 2D processors grid ( ortho_comm )
            !
            call cyc2blk_redist(n, vv, nrlx, umatl, size(umatl, 1), desc(1))
            !
            deallocate (diag, vv)

            if (rank == 0) then
                write (6, *) ' Lowest/Max  eigenvalue overlap mat =', eig(n), eig(1)
                condnumber = abs(eig(n)/eig(1))
                if (eig(1) .lt. 0.d0) condnumber = condnumber/100.d0
                write (6, *) 1, eig(1), 1.d0
                do i = 2, n
                    cost = abs(eig(i)/eig(1))
                    if (eig(i) .lt. 0.d0) cost = cost/100.d0
                    write (6, *) i, eig(i), cost
                    if (cost .lt. condnumber) condnumber = cost
                end do
                write (6, *) ' Inverse Condition Number basis set =', condnumber
            end if

            do i = 1, n
                if ((eig(i)/eig(1) .gt. abs(eps) .or. (eig(i) .lt. 0 .and. -eig(i)/eig(1)&
               &.gt. abs(eps*100.d0))) .and. i .le. n - mincond + 1) then ! the condition number criterium
                    eigmat(i) = dsqrt(1.d0/abs(eig(i)))
                else
                    mine = mine + 1
                    eigmat(i) = 0.d0
                end if
            end do

            !       we assume here that the garbage eigenvectors are the ones
            !       close to zero eigenvalue.

            if (mine .ne. 1 .and. rank .eq. 0 .and. optprint .ne. 0) write (6, *) ' disregarded coll. =', mine - 1

            !       we assume here that the garbage eigenvectors are the ones
            !       close to zero eigenvalue.

!        rescale umatl
!    the matrix used  for the transformation in the new basis
            do j = 1, desc(nlac_)
                do i = 1, desc(nlar_)
                    umatl(i, j) = umatl(i, j)*premat(i + ir - 1)*eigmat(j + ic - 1)
                end do
            end do

            call PDGEMM('N', 'N', n, n, n, 1.0d0, oversav, 1, 1, desch &
                        , umatl(1, 1), 1, 1, desch, 0.0d0, overs, 1, 1, desch)

            call PDGEMM('T', 'N', n, n, n, 1.0d0, umatl(1, 1), 1, 1, desch &
                        , overs(1, 1), 1, 1, desch, 0.0d0, overlap(1, 1), 1, 1, desch)

!      do j=1,desc( nlac_ )
!         do i=1,desc( nlar_ )
!            overlap(i,j)=overlap(i,j)*eigmat(i+ir-1)*eigmat(j+ic-1)
!         enddo
!      enddo

            maxeig = 100

!    diagonalize again overlap

            do i = 1, n
                ii = INDXG2L(i, desc(nlax_), me_ortho(1), 0, np_ortho(1))
                jj = INDXG2L(i, desc(nlax_), me_ortho(2), 0, np_ortho(2))
                ip = INDXG2P(i, desc(nlax_), me_ortho(1), 0, np_ortho(1))
                jp = INDXG2P(i, desc(nlax_), me_ortho(2), 0, np_ortho(2))
                if (ip == me_ortho(1) .and. jp == me_ortho(2) .and. eigmat(i) == 0.d0) then
                    overlap(ii, jj) = maxeig + dble(i)
                end if
            end do

!        if(rank.eq.0) then
!         write(6,*) ' Output matrix overlap '
!         do i=1,n
!          write(6,*) i,overlap(i,i)
!         enddo
!        endif

            molecorb = 0.d0

            nrlx = desc(la_nrlx_)

            allocate (diag(nrlx, n), vv(nrlx, n))
            diag = 0.d0
            vv = 0.d0
            !
            call blk2cyc_redist(n, diag, nrlx, overlap, size(overlap, 1), desc(1))
            !
!! maxeig should be calculated here as max_i   sum_j |mat_ij| and then
!        replaced in the diagonal elements with eigmat=0.

!     Here iopt=0  the standard one.

            call pdspev_drv_ss('V', diag, nrlx, eig, vv, nrlx, desc(la_nrl_), n,&
        & desc(la_npr_)*desc(la_npc_), desc(la_me_), desc(la_comm_), eps_mach, 0)
            !
            !  Redistribute matrix "vv" into "s"
            !  matrix "s" is block distributed
            !  across 2D processors grid ( ortho_comm )
            !
            call cyc2blk_redist(n, vv, nrlx, molecorb, size(molecorb, 1), desc(1))
            !
            deallocate (diag, vv, overlap)

!         if(rank.eq.0) then
!         write(6,*) ' Eigenvalues new overlap matrix ~ 1  '
!         do i=1,n
!         write(6,*) i,eig(i)
!         enddo
!         endif

            do i = 1, n
            if (eig(i) .le. eps) then
                if (rank .eq. 0 .and. eig(i) .gt. -1.d0) write (6, *) ' Further singular eigenvalue ', i
                eigmat(i) = 0.d0
            end if
            end do

            !    the matrix used  to go back to the original representation
            do j = 1, desc(nlac_)
                do i = 1, desc(nlar_)
                    overs(i, j) = umatl(i, j)
                end do
            end do

            do j = 1, desc(nlac_)
            if (eig(j + ic - 1) .gt. eps) then
                do i = 1, desc(nlar_)
                    molecorb(i, j) = molecorb(i, j)/dsqrt(eig(j + ic - 1))
                end do
            else
                do i = 1, desc(nlar_)
                    molecorb(i, j) = 0.d0
                end do
            end if
            end do
!    FINAL  PDGEMM
            call PDGEMM('N', 'N', n, n, n, 1.0d0, overs(1, 1), 1, 1 &
                        , desch, molecorb(1, 1), 1, 1, desch, 0.0d0, umatl(1, 1), 1, 1, desch)

        else

!   Calculation new lwork for next iterations
            do i = 1, n
                if (eigmat(i) .eq. 0.d0) mine = mine + 1
            end do
            dimo = n - mine + 1

        end if

        !        destroy overs mat_in
        call PDGEMM('N', 'N', n, n, n, 1.0d0, mat_in(1, 1), 1, 1 &
                    , desch, umatl(1, 1), 1, 1, desch, 0.0d0, &
                    overs(1, 1), 1, 1, desch)
        call PDGEMM('T', 'N', n, n, n, 1.0d0, umatl(1, 1), 1, 1 &
                    , desch, overs(1, 1), 1, 1, desch, 0.0d0, &
                    mat_in(1, 1), 1, 1, desch)

        do j = 1, desc(nlac_)
            do i = 1, desc(nlar_)
                if (eigmat(i + ir - 1) .eq. 0.d0 .or. eigmat(j + ic - 1) .eq. 0.d0) mat_in(i, j) = 0.d0
            end do
        end do

        maxeig = max(1000.d0, 10.d0*n)

        do i = 1, n
            ii = INDXG2L(i, desc(nlax_), me_ortho(1), 0, np_ortho(1))
            jj = INDXG2L(i, desc(nlax_), me_ortho(2), 0, np_ortho(2))
            ip = INDXG2P(i, desc(nlax_), me_ortho(1), 0, np_ortho(1))
            jp = INDXG2P(i, desc(nlax_), me_ortho(2), 0, np_ortho(2))
            if (ip == me_ortho(1) .and. jp == me_ortho(2) .and. eigmat(i) == 0.d0) then
                mat_in(ii, jj) = maxeig + dble(10*i)
            end if
        end do

        molecorb = 0.d0

        nrlx = desc(la_nrlx_)

        allocate (diag(nrlx, n), vv(nrlx, n))
        diag = 0.d0
        vv = 0.d0
        !
        call blk2cyc_redist(n, diag, nrlx, mat_in, size(mat_in, 1), desc(1))
        !
!!   maxeig should be calculated here as max_i   sum_j |mat_ij| and then
!        replaced in the diagonal elements with eigmat=0.

!    Here the standard one.
        call pdspev_drv_ss('V', diag, nrlx, eig, vv, nrlx, desc(la_nrl_), n, &
       &desc(la_npr_)*desc(la_npc_), desc(la_me_), desc(la_comm_), eps_mach, 0)
        !
        !  Redistribute matrix "vv" into "s"
        !  matrix "s" is block distributed
        !  across 2D processors grid ( ortho_comm )
        !
        call cyc2blk_redist(n, vv, nrlx, molecorb, size(molecorb, 1), desc(1))
        !
        deallocate (diag, vv)

        overs = molecorb

!    call pdgeadd('N',n,n,1.d0,molecorb,1,1,desch,0.d0,overs,1,1,desch )

        call PDGEMM('N', 'N', n, n, n, 1.0d0, umatl(1, 1), 1, 1 &
                    , desch, overs(1, 1), 1, 1, desch, 0.0d0, molecorb(1, 1), 1, 1, desch)

        !  if(rank.eq.0.and.optprint.ne.0) then
        !    do i=1,bands
        !        write(6,*) i,eig(i)
        !    enddo
        !  endif

        deallocate (overs, mat_in)

    end if

! In case not all processors are involved in the diagonalization
#ifdef PARALLEL
! call mpi_bcast(eig,nelorb_c,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call bcast_real(eig, nelorb_c, 0, commrep_mpi)
#endif

#else

    if (rankrep .eq. 0) then

        if (.not. allocated(umatl)) then
            allocate (umatl(nelorb_c, nelorb_c))
            umatl = 0.d0
        end if

        allocate (mat_in(nelorb_c, nelorb_c)&
                &, work(max(30*nelorb_c, lworkr)), overs(lda, nelorb_c)&
                &, iwork(5*nelorb_c), ifail(nelorb_c))
        mat_in = 0.d0
        work = 0.d0
        overs = 0.d0
        iwork = 0
        ifail = 0
        ! destroy invo
        !        calculation overlap old
        !        calculation normalized orbitals and eigenvectors
        ! maximum accuracy
        abstol = 2.d0*dlamch('S')
        eps = epsr

        if (eps .eq. 0.d0) eps = eps_mach

        n = nelorb_c
        lwork = 30*n
        mine = 1

        !    This code compute the generalized eigenvalue equation:
        !      mat_in x   y = over_s lambda y
        !      by stable preconditioning and pre-diagonalization of the matrix over_s

        overs = oversav
        minover = -1.d0

        if (iopt .eq. 1) then

            !         do i=1,nelorb_c
            !         if(overs(i,i).gt.0.d0) then
            !          if(overs(i,i).lt.minover.or.minover.eq.-1.d0) minover=overs(i,i)
            !          premat(i)=1.d0/dsqrt(overs(i,i))
            !         else
            !          if(rank.eq.0) write(6,*) ' Singular norm for orbital ',i,overs(i,i)
            !          premat(i)=-1.d0
            !         endif
            !         enddo
            !         do i=1,nelorb_c
            !         if(premat(i).eq.-1.d0) premat(i)=1.d0/dsqrt(minover)
            !         enddo

            !         There is numerical instability to use premat ne 1
            !         the orthogonality of the eigenvectors
            !         with respect to the metric given by overs  is poorly satisfied
            !         if premat ne 1 is used.
            premat = 1.d0

            !         reduced matrix
            do i = 1, nelorb_c
                do j = 1, nelorb_c
                    overs(i, j) = overs(i, j)*premat(i)*premat(j)
                end do
            end do

            call dsyevx('V', 'A', 'L', n, overs, lda, 0.d0, 0.d0, 1, 1, abstol       &
                    &, neig, eig, umatl, nelorb_c, work, lwork, iwork, ifail, info)

            if (rank .eq. 0) then
                write (6, *) ' Lowest/Max  eigenvalue overlap mat =', eig(1), eig(n)
                condnumber = abs(eig(1)/eig(n))
                if (eig(1) .lt. 0) condnumber = condnumber/100.d0
                do i = 2, n
                    cost = abs(eig(i)/eig(n))
                    if (eig(i) .lt. 0.d0) cost = cost/100.d0
                    if (cost .lt. condnumber) condnumber = cost
                end do
                write (6, *) ' Condition number basis set =', condnumber
            end if

            do i = 1, n
                if ((eig(i)/eig(n) .gt. abs(eps) .or. (eig(i) .lt. 0 .and. -eig(i)/eig(n)&
                        &.gt. abs(eps*100.d0))) .and. i .ge. mincond) then ! the condition number criterium
                    eigmat(i) = dsqrt(1.d0/eig(i))
                else
                    mine = mine + 1
                    eigmat(i) = 0.d0
                end if
            end do

            if (info .gt. 0 .and. rank .eq. 0) write (6, *)                         &
                    &' info > 0 in dsyevx !!! ', info

            !       we assume here that the garbage eigenvectors are the ones
            !       close to zero eigenvalue.
            do i = 1, info
                if (eigmat(ifail(i)) .ne. 0.d0) then
                    if (ifail(i) .ge. mine) mine = ifail(i) + 1
                    eigmat(1:ifail(i)) = 0.d0
                end if
            end do

            if (mine .ne. 1 .and. rank .eq. 0 .and. optprint .ne. 0)                   &
                    &write (6, *) ' disregarded coll. =', mine - 1

            !       first transformation  umatl

            do i = 1, n
                umatl(:, i) = umatl(:, i)*eigmat(i)
            end do

            !       now compute the overlap matrix in the new basis
            call dgemm('N', 'N', n, n, n, 1.d0, oversav, lda, umatl, n, 0.d0, mat_in, n)
            call dgemm('T', 'N', n, n, n, 1.d0, umatl, n, mat_in, lda, 0.d0, overs, lda)

            do i = 1, n
                if (eigmat(i) .eq. 0) overs(i, i) = -1000.d0 - dble(i)
            end do

            call dsyevx('V', 'A', 'L', n, overs, lda, 0.d0, 0.d0, 1, 1, abstol       &
                    &, neig, eig, mat_in, n, work, lwork, iwork, ifail, info)

            do i = 1, n
                if (eig(i) .gt. eps) then
                    mat_in(:, i) = mat_in(:, i)/dsqrt(eig(i))
                else
                    if (rankrep .eq. 0 .and. eig(i) .gt. -1.d0) write (6, *) ' Further singular eigenvalue ', i
                    eigmat(i) = 0.d0
                    mat_in(:, i) = 0.d0
                end if
            end do
            overs(1:n, 1:n) = umatl(1:n, 1:n)

            !       Final dgemm
            call dgemm('N', 'N', n, n, n, 1.d0, overs, lda, mat_in, n, 0.d0, umatl, n)

        else

            do i = 1, n
                if (eigmat(i) .eq. 0.d0) mine = mine + 1
            end do

        end if

        !        Now modifying the matrix detmat_c (and forgetting it)
        mat_in(1:n, 1:n) = matsav(1:n, 1:n)
        do i = 1, n
            do j = 1, n
                mat_in(i, j) = mat_in(i, j)*premat(i)*premat(j)
            end do
        end do

        !        destroy overs mat_in
        call dgemm('N', 'N', n, n, n, 1.d0, mat_in, n, umatl, n, 0.d0, overs, lda)
        call dgemm('T', 'N', n, n, n, 1.d0, umatl, n, overs, lda, 0.d0, mat_in, n)

        !        do i=1,n
        !          do j=1,n
        !          mat_in(i,j)=mat_in(i,j)*eigmat(i)*eigmat(j)
        !          enddo
        !        enddo

        if (lworkr .eq. 0) then
            lwork = 30*n
        else
            lwork = lworkr
        end if

        dimo = n - mine + 1

        molecorb = 0.d0
        !      avoid to include the garbage matrix elements.

        if (bands .lt. n) then

            call dsyevx('V', 'I', 'L', dimo, mat_in(mine, mine), n, 0.d0, 0.d0, 1, bands &
                    &, abstol, neig, eig, molecorb(mine, 1), lda, work, lwork, iwork   &
                    &, ifail, info)

        else

            call dsyevx('V', 'A', 'L', dimo, mat_in(mine, mine), n, 0.d0, 0.d0, 1, 1 &
                    &, abstol, neig, eig, molecorb(mine, 1), lda, work, lwork, iwork   &
                    &, ifail, info)

        end if

        if (lworkr .eq. 0) then
            if (info .eq. 0) then
                lworkr = work(1)
                if (rank .eq. 0) write (6, *) ' Optimal lwork  found =', lworkr
            else
                lworkr = 30*n
            end if
        end if
        if (rank .eq. 0 .and. info .ne. 0) &
                &write (6, *) ' Warning info ne 0 in dsyevx ', info

        overs = molecorb

        call dgemm('N', 'N', n, bands, n, 1.d0, umatl, n               &
                &, overs, lda, 0.d0, molecorb, lda)

        !        choose a gauge
        do i = 1, bands
            orbmax = sum(molecorb(1:nelorb_c, i))
            !           do j=2,nelorb_c
            !           if(abs(molecorb(j,i)).gt.abs(orbmax)) orbmax=molecorb(j,i)
            !           enddo
            if (orbmax .lt. 0) molecorb(:, i) = -molecorb(:, i)
        end do

        !        if(info.eq.0) then
        !if(rank.eq.0.and.optprint.ne.0) then
        ! do i=1,bands
        ! write(6,*) i,eig(i)
        ! enddo
        !endif

        !         write(6,*) ' Eigenvectors  '
        !         do i=1,nelorb_c
        !         write(6,*) i,(molecorb(j,i),j=1,nelorb_c)
        !         enddo
        !        endif
        deallocate (mat_in, work, iwork, ifail, overs)
    end if

#ifdef PARALLEL
!        bcast just the relevant info to the nodes
    dimorb = lda*(nelorb_c - 1) + nelorb_c
    call bcast_real(molecorb, dimorb, 0, commrep_mpi)
    call bcast_real(eig, nelorb_c, 0, commrep_mpi)
#endif

#endif

    !call mpi_finalize(ierr)
    !stop

    return
end
