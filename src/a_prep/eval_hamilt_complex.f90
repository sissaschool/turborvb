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

!
! complex version of eval_hamilt.f90
! It uses the same algorithm for diagonalization on a non-orthogonal basis.
!
subroutine eval_hamilt_complex(nelorb_c, oversav, matsav, molecorb, umatl, eig, lda, epsr, &
                               eps_mach, rank, optprint, lworkr, iopt, premat, eigmat, &
                               info, bands, nlax, mincond)

#ifdef __SCALAPACK
    use allio, only: ortho_cntx, me_blacs, np_ortho, me_ortho, ortho_comm, &
                     ortho_comm_id
    use descriptors
    use dspev_module
#endif
    use kpoints_mod, only: xkp, nk, kaverage, compute_bands
    use allio, only: commrep_mpi, commcolrep_mpi, rankrep, rankcolrep
    use allio, only: rank_print => rank ! for printing eigenvalues with k-points
    use constants
    use setup, only: indk

    implicit none

    integer, intent(in) :: nlax
    integer nelorb_c, i, j, rank, info, n, m, lda, mine, dimo, neig, lwork, optprint&
            &, lworkr, iopt, countexc, bands, dimorb, ierr, mincond
    real(8) eig(nelorb_c), premat(nelorb_c)&
            &, eigmat(*), abstol, eps, epsr, orbmax, dlamch, minover&
            &, condnumber, cost, maxeig, eps_mach, trace, print_eig(nelorb_c)
    integer, dimension(:), allocatable :: iwork, ifail
    ! auxiliary variables for printing
    real(8), dimension(:), allocatable :: eigov
    real(8), dimension(:, :), allocatable :: eigham
    !
#ifdef __SCALAPACK
    complex(8) molecorb(nlax, nlax), oversav(nlax, nlax), matsav(nlax, nlax)
    complex(8), dimension(:, :), allocatable :: mat_in, overs
    complex(8), dimension(:), allocatable :: work
    !
    integer :: desch(20)
    integer :: desc(descla_siz_)
    integer :: LIWORK, ic, ir, mm, nz, ii, jj, ip, jp, nrlx
    real(8) :: PDLAMCH
    integer :: INDXG2L, INDXG2P, INDXL2G
    complex(8) :: umatl(nlax, nlax)
    complex(8), allocatable :: diag(:, :), vv(:, :), overlap(:, :)
    real(8), dimension(:), allocatable :: rwork
#else
    complex(8) molecorb(lda, nelorb_c), oversav(lda, nelorb_c), matsav(lda, nelorb_c)
    complex(8), dimension(:, :), allocatable :: mat_in, overs
    complex(8), dimension(:), allocatable :: work
    complex(8) :: umatl(nelorb_c, nelorb_c)
    real(8), dimension(:), allocatable :: rwork
#endif

#ifdef PARALLEL
    include 'mpif.h'
#endif

    premat = 1.d0
    allocate (eigov(nelorb_c), eigham(nelorb_c, nk))
    eigov = 0.d0
    eigham = 0.d0

#ifdef __SCALAPACK

    n = nelorb_c
    mine = 1

    call descla_init(desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id)

    if (desc(lambda_node_) > 0) then

        allocate (mat_in(desc(nlax_), desc(nlax_)))
        mat_in = zzero
        !IF( .NOT. ALLOCATED( umatl ) ) then
        !   ALLOCATE( umatl( desc( nlax_ ), desc( nlax_ ) ) )
        !   umatl=zzero
        !ENDIF

        allocate (overs(desc(nlax_), desc(nlax_)))
        overs = zzero
        ir = desc(ilar_)
        ic = desc(ilac_)

        !    if(rank.eq.0) write(6,*) ' abstol =',abstol

        eps = epsr

        if (eps .eq. 0.d0) eps = abs(eps_mach)

        !    For the first part we do not assume to have accurate eigenvectors.

        call descinit(desch, n, n, desc(nlax_), desc(nlax_), 0, 0, ortho_cntx, size(umatl, 1), info)

        if (info /= 0) call errore(' eval_hamilt_complex ', ' descinit ', abs(info))

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
            overlap = zzero
!       trace=0.d0

            do j = 1, desc(nlac_)
                do i = 1, desc(nlar_)
                    overs(i, j) = overs(i, j)*premat(i + ir - 1)*premat(j + ic - 1)
!             overs(i,j)=real(overs(i,j))
!          trace=trace+abs(dimag(overs(i,j)))
                end do
            end do
!       write(6,*) ' Input imag =',trace

            nrlx = desc(la_nrlx_)

            allocate (diag(nrlx, n), vv(nrlx, n))

            diag = zzero
            vv = zzero
            !
            call blk2cyc_zredist(n, diag, nrlx, overs, size(overs, 1), desc(1))
            !
            !    iopt=0, the standard one
            call zdspev_drv_ss('V', diag, nrlx, eig, vv, nrlx, desc(la_nrl_), n,&
                 &desc(la_npr_)*desc(la_npc_), desc(la_me_), desc(la_comm_), eps_mach, 1)

            !  Redistribute matrix "vv" into "s"
            !  matrix "s" is block distributed
            !  across 2D processors grid ( ortho_comm )
            !
            call cyc2blk_zredist(n, vv, nrlx, umatl, size(umatl, 1), desc(1))
            !
            deallocate (diag, vv)

            if (rank_print .eq. 0) then
!          write(6,*) ' Lowest/Max  eigenvalue overlap mat =',eig(n),eig(1)
                condnumber = abs(eig(n)/eig(1))
!          write(6,*) 1,eig(1),1.d0
                if (eig(1) .lt. 0.d0) condnumber = condnumber/100.d0
                do i = 2, n
                    cost = abs(eig(i)/eig(1))
                    if (eig(i) .lt. 0.d0) cost = cost/100.d0
                    if (cost .lt. condnumber) condnumber = cost
!          write(6,*) i,eig(i),cost
                end do
!          write(6,*) ' Inverse Condition number basis set =',condnumber
            end if

            do i = 1, n
                eigov(i) = eig(n - i + 1) ! save overlap eigenvalues for printing
            end do

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
            if (mine .ne. 1 .and. rank_print .eq. 0 .and. optprint .ne. 0) write (6, *) ' disregarded coll. =', mine - 1

            !       we assume here that the garbage eigenvectors are the ones
            !       close to zero eigenvalue.

            !        rescale umatl
            !    the matrix used  for the transformation in the new basis
            do j = 1, desc(nlac_)
                do i = 1, desc(nlar_)
                    umatl(i, j) = umatl(i, j)*premat(i + ir - 1)*eigmat(j + ic - 1)
                end do
            end do

            call PZGEMM('N', 'N', n, n, n, (1.0d0, 0.d0), oversav, 1, 1, desch, &
                        umatl(1, 1), 1, 1, desch, (0.0d0, 0.d0), overs, 1, 1, desch)

            call PZGEMM('C', 'N', n, n, n, (1.0d0, 0.d0), umatl(1, 1), 1, 1, desch, &
                        overs(1, 1), 1, 1, desch, (0.0d0, 0.d0), overlap(1, 1), 1, 1, desch)

            maxeig = 100

            ! diagonalize again overlap
            do i = 1, n
                ii = INDXG2L(i, desc(nlax_), me_ortho(1), 0, np_ortho(1))
                jj = INDXG2L(i, desc(nlax_), me_ortho(2), 0, np_ortho(2))
                ip = INDXG2P(i, desc(nlax_), me_ortho(1), 0, np_ortho(1))
                jp = INDXG2P(i, desc(nlax_), me_ortho(2), 0, np_ortho(2))
                if (ip == me_ortho(1) .and. jp == me_ortho(2) .and. eigmat(i) == 0.d0) then
                    overlap(ii, jj) = maxeig + dble(i)
                end if
            end do

!       allocate(vv(n,n))
!       if(rank_print.eq.0) then
!
!        trace=0.d0
!        vv(1:n,1:n)=real(overlap(1:n,1:n))
!        vv(1:n,1:n)=overlap(1:n,1:n)
!        write(6,*) ' Output matrix overlap '
!        do i=1,n
!         do j=1,n
!         write(6,*) i,overlap(i,j)
!         overlap(i,j)=vv(n-i+1,n-j+1)
!         overlap(i,j)=vv(i,j)
!         enddo
!         trace=trace+vv(i,i)
!        enddo
!        write(6,*) ' Input trace =',trace
!       endif
!       deallocate(vv)
!     abstol=2.d0*dlamch('S')

!     allocate(work(max(30*nelorb_c,lworkr)),&
!              iwork(5*nelorb_c),ifail(nelorb_c),rwork(7*nelorb_c))
!        call zheevx('V','A','L',n,overlap,lda,0.d0,0.d0,1,1,abstol       &
!             &,neig,eig,mat_in,n,work,lwork,rwork,iwork,ifail,info)

!        go to 23

            molecorb = zzero

            nrlx = desc(la_nrlx_)

            allocate (diag(nrlx, n), vv(nrlx, n))
            diag = 0.d0
            vv = 0.d0
            !
            call blk2cyc_zredist(n, diag, nrlx, overlap, size(overlap, 1), desc(1))
            !
        !! maxeig should be calculated here as max_i   sum_j |mat_ij| and then
            !        replaced in the diagonal elements with eigmat=0.

            !     Here iopt=0  the standard one.

            call zdspev_drv_ss('V', diag, nrlx, eig, vv, nrlx, desc(la_nrl_), n,&
       & desc(la_npr_)*desc(la_npc_), desc(la_me_), desc(la_comm_), eps_mach, 0)
            !
            !  Redistribute matrix "vv" into "s"
            !  matrix "s" is block distributed
            !  across 2D processors grid ( ortho_comm )
            !
            call cyc2blk_zredist(n, vv, nrlx, molecorb, size(molecorb, 1), desc(1))
            !
            deallocate (diag, vv, overlap)

!23 continue
!        allocate(vv(n,n))
!        vv(1:n,1:n)=molecorb(1:n,1:n)
!        do i=1,n
!        molecorb(n-i+1,:)=vv(i,:)
!        enddo
!        deallocate(vv)

#ifdef DEBUG
            if (rank .eq. 0) then
                write (6, *) ' Eigenvalues new overlap matrix ~ 1  '
                do i = 1, n
                    write (6, *) i, eig(i)
                end do
            end if
#endif

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
                        molecorb(i, j) = zzero
                    end do
                end if
            end do
            !    FINAL  PDGEMM
            call PZGEMM('N', 'N', n, n, n, (1.0d0, 0.d0), overs(1, 1), 1, 1, desch, &
                        molecorb(1, 1), 1, 1, desch, (0.0d0, 0.d0), umatl(1, 1), 1, 1, desch)

        else

            ! Calculation new lwork for next iterations
            do i = 1, n
                if (eigmat(i) .eq. 0.d0) mine = mine + 1
            end do
            dimo = n - mine + 1

        end if

        !write(6,*) rankcolrep,xkp(:,indk),sum(mat_in(:,:)),rankrep
        !
        !if(rankcolrep.eq.0) then
        !   do i=1,desc(nlax_)
        !           do j=1,desc(nlax_)
        !             write(6,*) i,j,mat_in(i,j)
        !      enddo
        !   enddo
        !endif

        !        destroy overs mat_in
        call PZGEMM('N', 'N', n, n, n, (1.0d0, 0.d0), mat_in(1, 1), 1, 1, desch, umatl(1, 1), 1, 1, desch, (0.0d0, 0.d0), &
                    overs(1, 1), 1, 1, desch)
        call PZGEMM('C', 'N', n, n, n, (1.0d0, 0.d0), umatl(1, 1), 1, 1, desch, overs(1, 1), 1, 1, desch, (0.0d0, 0.d0), &
                    mat_in(1, 1), 1, 1, desch)

        do j = 1, desc(nlac_)
            do i = 1, desc(nlar_)
                if (eigmat(i + ir - 1) .eq. 0.d0 .or. eigmat(j + ic - 1) .eq. 0.d0) mat_in(i, j) = zzero
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

        molecorb = zzero

        nrlx = desc(la_nrlx_)

        allocate (diag(nrlx, n), vv(nrlx, n))
        diag = 0.d0
        vv = 0.d0
        call blk2cyc_zredist(n, diag, nrlx, mat_in, size(mat_in, 1), desc(1))
        !
     !!   maxeig should be calculated here as max_i   sum_j |mat_ij| and then
        !        replaced in the diagonal elements with eigmat=0.

        !    Here the standard one.
        call zdspev_drv_ss('V', diag, nrlx, eig, vv, nrlx, desc(la_nrl_), n, &
             &desc(la_npr_)*desc(la_npc_), desc(la_me_), desc(la_comm_), eps_mach, 0)
        !
        !  Redistribute matrix "vv" into "s"
        !  matrix "s" is block distributed
        !  across 2D processors grid ( ortho_comm )
        !
        call cyc2blk_zredist(n, vv, nrlx, molecorb, size(molecorb, 1), desc(1))
        !
        deallocate (diag, vv)

        overs = molecorb

        call PZGEMM('N', 'N', n, n, n, (1.0d0, 0.d0), umatl(1, 1), 1, 1 &
                    , desch, overs(1, 1), 1, 1, desch, (0.0d0, 0.d0), molecorb(1, 1), 1, 1, desch)

#ifdef DEBUG
        if (rank .eq. 0) then
            write (6, *) ' # Eigenvectors Hamiltonian matrix'
            do i = 1, nlax
                write (6, *) 'Eigenvector # i mod, real phase', i
                do j = 1, nlax
                    write (6, *) j, real(abs(molecorb(j, i))), &
                        real(molecorb(j, i)/abs(molecorb(j, i)))
                end do
            end do
        end if
#endif

        deallocate (overs, mat_in)

    end if

    ! In case not all processors are involved in the diagonalization
#ifdef PARALLEL
    call bcast_real(eig, nelorb_c, 0, commrep_mpi)
#endif

    eigham(:, indk) = eig ! save hamiltonian eigenvalues for printing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !if(rankcolrep.eq.1) then
    !   write(6,*) 'eig',xkp(:,rankcolrep+1)
    !   do i=1,nelorb_c
    !      write(6,*) i,eig(i)
    !   enddo
    !endif
    !write(6,*) rankcolrep,sum(molecorb(:,:)),rankrep
    !call mpi_finalize(ierr)
    !stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#else

    if (rankrep .eq. 0) then

        !if(.not.allocated(umatl)) allocate(umatl(nelorb_c,nelorb_c))

        allocate (mat_in(nelorb_c, nelorb_c), work(max(30*nelorb_c, lworkr)), overs(lda, nelorb_c), &
                  iwork(5*nelorb_c), ifail(nelorb_c), rwork(7*nelorb_c))
        mat_in = 0.d0
        work = 0.d0
        overs = 0.d0
        iwork = 0
        ifail = 0
        rwork = 0.d0

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
        !    mat_in x   y = over_s lambda y
        !    by stable preconditioning and pre-diagonalization of the matrix over_s

        overs = oversav
        minover = -1.d0

        if (iopt .eq. 1) then

            !         There is numerical instability to use premat ne 1
            !         the orthogonality of the eigenvectors
            !         with respect to the metric given by overs is poorly satisfied
            !         if premat ne 1 is used.

            premat = 1.d0
            ! reduced matrix
            do i = 1, nelorb_c
                do j = 1, nelorb_c
                    overs(i, j) = overs(i, j)*premat(i)*premat(j)
                end do
            end do

            call zheevx('V', 'A', 'L', n, overs, lda, 0.d0, 0.d0, 1, 1, abstol   &
                    &, neig, eig, umatl, nelorb_c, work, lwork, rwork, iwork, ifail, info)
#ifdef DEBUG
            if (rank .eq. 0) then
                write (6, *) '# Eigenvalues Overlap matrix'
                do i = 1, n
                    write (6, *) i, eig(i)
                end do
            end if
#endif

            eigov = eig ! save overlap eigenvalues for printing

            condnumber = abs(eig(1)/eig(n))
            if (eig(1) .lt. 0) condnumber = condnumber/100.d0
            do i = 2, n
                cost = abs(eig(i)/eig(n))
                if (eig(i) .lt. 0.d0) cost = cost/100.d0
                if (cost .lt. condnumber) condnumber = cost
            end do

            do i = 1, n
                if ((eig(i)/eig(n) .gt. abs(eps) .or. (eig(i) .lt. 0 .and. -eig(i)/eig(n)&
                        &.gt. abs(eps*100.d0))) .and. i .ge. mincond) then ! the condition number criterium
                    eigmat(i) = dsqrt(1.d0/eig(i))
                else
                    mine = mine + 1
                    eigmat(i) = 0.d0
                end if
            end do

            if (info .gt. 0 .and. rank_print .eq. 0) write (6, *) ' info > 0 in zhpevx !!! ', info
            ! we assume here that the garbage eigenvectors are the ones
            ! close to zero eigenvalue.
            do i = 1, info
                if (eigmat(ifail(i)) .ne. 0.d0) then
                    if (ifail(i) .ge. mine) mine = ifail(i) + 1
                    eigmat(1:ifail(i)) = 0.d0
                end if
            end do

            if (mine .ne. 1 .and. optprint .ne. 0 .and. rank_print .eq. 0)                   &
                    &write (6, *) ' disregarded coll. =', mine - 1

            ! first transformation  umatl
            do i = 1, n
                umatl(:, i) = umatl(:, i)*dcmplx(eigmat(i))
            end do

            ! now compute the overlap matrix in the new basis
            call zgemm('N', 'N', n, n, n, (1.d0, 0.d0), oversav, lda, umatl, n, (0.d0, 0.d0), mat_in, n)
            call zgemm('C', 'N', n, n, n, (1.d0, 0.d0), umatl, n, mat_in, lda, (0.d0, 0.d0), overs, lda)

#ifdef DEBUG
            if (rank .eq. 0) then
                write (6, *) ' Overlap matrix '
                do i = 1, n
                    do j = i, n
                        write (6, *) i, j, overs(i, j)
                    end do
                end do
            end if
#endif

            do i = 1, n
                if (eigmat(i) .eq. 0) overs(i, i) = -1000.d0 - dcmplx(i)
            end do

            call zheevx('V', 'A', 'L', n, overs, lda, 0.d0, 0.d0, 1, 1, abstol       &
                    &, neig, eig, mat_in, n, work, lwork, rwork, iwork, ifail, info)

            do i = 1, n
                if (eig(i) .gt. eps) then
                    mat_in(:, i) = mat_in(:, i)/(zone*dsqrt(eig(i)))
                else
                    if (rank_print .eq. 0 .and. (eig(i) .gt. -1.d0)) write (6, *) ' Further singular eigenvalue ', i
                    eigmat(i) = 0.d0
                    mat_in(:, i) = zzero
                end if
            end do
            overs(1:n, 1:n) = umatl(1:n, 1:n)

            ! final zgemm
            call zgemm('N', 'N', n, n, n, (1.d0, 0.d0), overs, lda, mat_in, n, (0.d0, 0.d0), umatl, n)

            !       if(rank.eq.0) then
            !          write(6,*) '# Eigenvalues Overlap matrix ~ 1'
            !          do i=1,n
            !             write(6,*) i,eig(i)
            !          enddo
            !       endif

        else

            do i = 1, n
                if (eigmat(i) .eq. 0.d0) mine = mine + 1
            end do

        end if

        ! Now modifying the matrix detmat_c (and forgetting it)
        mat_in(1:n, 1:n) = matsav(1:n, 1:n)
        do i = 1, n
            do j = 1, n
                mat_in(i, j) = mat_in(i, j)*premat(i)*premat(j)
            end do
        end do

        ! destroy overs mat_in
        call zgemm('N', 'N', n, n, n, (1.d0, 0.d0), mat_in, n, umatl, n, (0.d0, 0.d0), overs, lda)
        call zgemm('C', 'N', n, n, n, (1.d0, 0.d0), umatl, n, overs, lda, (0.d0, 0.d0), mat_in, n)

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
        ! avoid to include the garbage matrix elements
        if (bands .lt. n) then

            call zheevx('V', 'I', 'L', dimo, mat_in(mine, mine), n, 0.d0, 0.d0, 1, bands &
                    &, abstol, neig, eig, molecorb(mine, 1), lda, work, lwork, rwork, iwork   &
                    &, ifail, info)

        else

            call zheevx('V', 'A', 'L', dimo, mat_in(mine, mine), n, 0.d0, 0.d0, 1, 1 &
                    &, abstol, neig, eig, molecorb(mine, 1), lda, work, lwork, rwork, iwork   &
                    &, ifail, info)

        end if

        if (lworkr .eq. 0) then
            if (info .eq. 0) then
                lworkr = work(1)
                if (rank_print .eq. 0) write (6, *) ' Optimal lwork found =', lworkr
            else
                lworkr = 30*n
            end if
        end if
        if (info .ne. 0 .and. rank_print .eq. 0) &
            write (6, *) ' Warning info ne 0 in zhpevx ', info

        overs = molecorb

        call zgemm('N', 'N', n, bands, n, (1.d0, 0.d0), umatl, n, overs, lda, (0.d0, 0.d0), molecorb, lda)

        ! choose a gauge
        do i = 1, bands
            orbmax = sum(real(molecorb(1:nelorb_c, i)))
            !           do j=2,nelorb_c
            !           if(abs(molecorb(j,i)).gt.abs(orbmax)) orbmax=molecorb(j,i)
            !           enddo
            if (orbmax .lt. 0) molecorb(:, i) = zmone*molecorb(:, i)
        end do

        eigham(:, indk) = eig ! save hamiltonian eigenvalues for printing

#ifdef DEBUG
        write (6, *) ' # Eigenvectors Hamiltonian matrix'
        do i = 1, nelorb_c
            write (6, *) i, molecorb(1:nelorb_c, i)
        end do
#endif
        deallocate (mat_in, work, iwork, ifail, overs, rwork)
    end if

#ifdef PARALLEL
    ! bcast just the relevant info to the nodes
    dimorb = lda*(nelorb_c - 1) + nelorb_c

    call bcast_complex(molecorb, dimorb, 0, commrep_mpi)
    call bcast_real(eig, nelorb_c, 0, commrep_mpi)

#endif

#endif
!
! print eigenvalues. If kaverage, it prints the eigenvalues for
! all k-points.
!
#ifdef PARALLEL
    if (kaverage) &
        call reduce_base_real(nk*nelorb_c, eigham, commcolrep_mpi, -1)
#endif

    if (optprint .ne. 0 .and. rank_print .eq. 0 .and. (.not. compute_bands)) &
        call print_eigenvalues(eigov, eigham, nelorb_c, bands, condnumber, iopt)
!
    deallocate (eigov, eigham)

    return

end subroutine eval_hamilt_complex
!
! print eigenvalues depending on optprint option
!
subroutine print_eigenvalues(eigov, eigham, dim, bands, condnumber, iopt)

    use allio, only: commcolrep_mpi, rank, xkp
    use freeelmod_complex, only: nk, indk

    implicit none

    integer, intent(in) :: iopt, bands, dim
    real(8), intent(in) :: eigov(*), eigham(dim, *), condnumber
    !
    integer i, j
    real(8) cost

    ! overlap eigenvalues
    if (iopt .eq. 1) then
        write (6, *) ' Lowest/Max  eigenvalue overlap mat =', eigov(1), eigov(dim)
        write (6, *) ' Eigenvalues Overlap '
        do i = 1, dim
            cost = abs(eigov(i)/eigov(dim))
            if (eigov(i) .lt. 0.d0) cost = cost/100.d0
            write (6, *) i, eigov(i), cost
        end do
        write (6, *) ' Inverse Condition Number basis set =', condnumber
    end if
#ifdef DEBUG
    ! hamiltoninan eigenvalues
    write (6, *) ' Eigenvalues Hamiltonian '
    do i = 1, nk
        do j = 1, bands
            write (6, *) j, eigham(j, i)
        end do
    end do
#endif

    return

end subroutine print_eigenvalues

