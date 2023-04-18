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

#include "f_defs.h"

module dspev_module

    implicit none
    save
    private

    public :: pdspev_drv, pdspev_drv_ss, zdspev_drv_ss, dspev_drv

contains

    subroutine ptredv(a, lda, d, e, v, ldv, nrl, n, nproc, me, comm, machp)
        !
        !     Parallel version of the famous HOUSEHOLDER tridiagonalization
        !     Algorithm for simmetric matrix.
        !
        !     AUTHOR : Carlo Cavazzoni - SISSA in 1997
        !     Modified: S. Sorella -SISSA in 2008
        !
        ! REFERENCES :
        !
        ! NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
        ! W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
        ! CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
        !
        ! PARALLEL NUMERICAL ALGORITHMS,
        ! T.L. FREEMAN AND C.PHILLIPS,
        ! PRENTICE HALL INTERNATIONAL (1992).
        !
        !
        !
        !     INPUTS :
        !
        !     A(NRL,N) Local part of the global matrix A(N,N) to be reduced,
        !              only the upper triangle is needed.
        !              The rows of the matrix are distributed among processors
        !              with blocking factor 1.
        !              Example for NPROC = 4 :
        !              ROW | PE
        !              1   | 0
        !              2   | 1
        !              3   | 2
        !              4   | 3
        !              5   | 0
        !              6   | 1
        !              ..  | ..
        !
        !     LDA      LEADING DIMENSION OF MATRIX A.
        !
        !     LDV      LEADING DIMENSION OF MATRIX V.
        !
        !     NRL      NUMBER OF ROWS BELONGING TO THE LOCAL PROCESSOR.
        !
        !     N        DIMENSION OF THE GLOBAL MATRIX.
        !
        !     NPROC    NUMBER OF PROCESSORS.
        !
        !     ME       INDEX OF THE LOCAL PROCESSOR (Starting from 0).
        !
        !
        !     OUTPUTS :
        !
        !     V(NRL,N) Orthogonal transformation that tridiagonalize A,
        !              this matrix is distributed among processor
        !              in the same way as A.
        !
        !     D(N)     Diagonal elements of the tridiagonal matrix
        !              this vector is equal on all processors.
        !
        !     E(N)     Subdiagonal elements of the tridiagonal matrix
        !              this vector is equal on all processors. E(1)=0
        !
        !
        !     machprec Input relative precision of eigenvalues if >=0. Otherwise
        !     lapack machine relative precision obtained by DLAMCH('E') is used.
        !     It may be useful to have less converged eigenvectors/eigenvalues
        !     with a cheap computation for large matrices.

        implicit none
#ifdef PARALLEL
    include 'mpif.h'
#endif
        integer, intent(in) :: N, NRL, LDA, LDV
        integer, intent(in) :: NPROC, ME, comm
        real*8 :: A(LDA, N), D(N), E(N), V(LDV, N)
        real*8, intent(in) :: machp
        !
        real*8, external :: ddot
        !
        real*8 :: g, scalef, sigma, kappa, f, h, tmp, machprec, safemin, maxmat&
                &, safemin2
        real*8, external :: DLAMCH
        real*8, ALLOCATABLE :: u(:)
        real*8, ALLOCATABLE :: p(:)
        real*8, ALLOCATABLE :: vtmp(:)

        real*8 :: tu, tp, one_over_h
        real*8 :: one_over_scale
        real*8 :: redin(3), redout(3)
        real*8, allocatable :: ul(:)
        real*8, allocatable :: asav(:)
        integer :: l, i, j, k, t, tl, ierr, countstrange
        integer :: kl, jl, ks, lloc
        integer, allocatable :: is(:)
        integer, allocatable :: ri(:)
#if defined (_OPENMP) && defined (__NOOMP)
        integer, external :: omp_get_max_threads
        integer old_threads
        old_threads = omp_get_max_threads()
        call omp_set_num_threads(1)   ! scalar code
#endif



        !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........

        if(N == 0) then
            return
        end if

        countstrange = 0

        !     The allocation of p below is necessary when nproc> n otherwise some
        !     call (DGER below) uses undefined memory (possible abort).
        allocate(u(n + 2), p(max(n + 1, nproc)), vtmp(n + 2), ul(n)&
                &, is(n), ri(n), Asav(lda))
        u = 0.d0
        p = 0.d0
        vtmp = 0.d0
        ul = 0.d0
        is = 0
        ri = 0
        asav = 0.d0

        do I = N, 1, -1
            IS(I) = (I - 1) / NPROC
            RI(I) = MOD((I - 1), NPROC)  !  owner of I-th row
            if(ME .le. RI(I)) then
                IS(I) = IS(I) + 1
            end if
        end do
#ifdef PARALLEL
      tmp=0.d0
      do i=1,N
        do k=1,is(i)
        tmp=max(tmp,abs(A(k,i)))
        end do
      end do
      call mpi_allreduce(tmp,maxmat,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM,ierr)
#else
        maxmat = 0.d0
        do i = 1, N
            do k = 1, i
                maxmat = max(maxmat, abs(A(k, i)))
            end do
        end do
#endif
        if(maxmat.eq.0.d0) then
            D = 0.d0
            E = 0.d0
            return
        end if

        machprec = abs(machp)
        safemin = machprec**2 * maxmat
        safemin2 = dlamch('S') * 1000000.d0 * N  ! safer minimum than lapack

        do I = N, 2, -1

            L = I - 1         ! first element
            H = 0.0d0

            if (L > 1) THEN

                SCALEF = 0.0d0
                do K = 1, is(l)
                    SCALEF = SCALEF + DABS(A(K, I))
                end do

#if defined PARALLEL
                tmp = scalef
                call mpi_allreduce(tmp, scalef, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, ierr)
#endif

                if (SCALEF .le. SAFEMIN) then
                    !
                    if (RI(L).eq.ME) then
                        tmp = A(is(L), I)
                    else
                        tmp = 0.d0
                    end if
#ifdef PARALLEL
     call mpi_allreduce(tmp,E(I),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM,ierr)
#else
                    E(I) = tmp
#endif
                ELSE

                    !  ......  CALCULATION OF SIGMA AND H

                    ONE_OVER_SCALE = 1.0d0 / SCALEF
                    SIGMA = 0.0d0
                    do k = 1, is(L)
                        asav(k) = A(k, i)
                        A(k, I) = A(k, I) * ONE_OVER_SCALE
                        SIGMA = SIGMA + A(k, I)**2
                    end do

                    if(ri(l) .eq. me) then
                        F = A(is(l), i)
                    else
                        F = 0.0d0
                    end if

                    !  CONSTRUCTION OF VECTOR U

                    vtmp(1:l) = 0.0d0

                    k = ME + 1
                    do kl = 1, is(l)
                        vtmp(k) = A(kl, I)
                        k = k + NPROC
                    end do

                    do kl = 1, is(l)
                        UL(kl) = A(kl, I)
                    end do

#if defined PARALLEL
                    vtmp(l + 1) = sigma
                    vtmp(l + 2) = f
                    call reduce_base_real_to(L + 2, vtmp, u, comm, -1)
                    sigma = u(l + 1)
                    f = u(l + 2)
#else
             u(1:l) = vtmp(1:l)
#endif

                    g = -sign(sqrt(sigma), f)
                    h = sigma - f * g

                end if

                if(H.lt.SAFEMIN2) then ! S.S. Added another control: in principle
                    ! H can be very small so that 1/H may overflow.
                    countstrange = countstrange + 1

                    !  restore the previous unscaled A
                    do k = 1, is(L)
                        A(k, I) = Asav(k)
                    end do
                    !            I am not sure E(I) is global in this case so mpi_allreduce needed
                    if (RI(L).eq.ME) then
                        tmp = A(is(L), I)
                    else
                        tmp = 0.d0
                    end if
#ifdef PARALLEL
     call mpi_allreduce(tmp,E(I),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM,ierr)
#else
                    E(I) = tmp
#endif

                    !            H=0.d0  do not make useless approximations.

                else

                    ONE_OVER_H = 1.0d0 / H
                    E(I) = SCALEF * G

                    U(L) = F - G

                    if(RI(L) ==  ME) THEN
                        UL(is(l)) = F - G
                        A(is(l), I) = F - G
                    end if

                    !  CONSTRUCTION OF VECTOR P
#ifdef __DOOMP
!$omp parallel do shared(vtmp,L,IS,A,UL,RI,U,ME,ONE_OVER_H)  private(J,KL,K)
#endif
                    do J = 1, L
                        vtmp(j) = 0.0d0

                        do KL = 1, IS(J)
                            vtmp(J) = vtmp(J) + A(KL, J) * UL(KL)
                        end do

                        if(L > J .and. ME == RI(J)) then
                            do K = J + 1, L
                                vtmp(J) = vtmp(J) + A(IS(J), K) * U(K)
                            end do
                        end if

                        vtmp(J) = vtmp(J) * ONE_OVER_H

                    end do
#ifdef __DOOMP
!$omp end parallel do
#endif

#ifdef PARALLEL
             call reduce_base_real_to( L, vtmp, p, comm, -1 )
#else
                    p(1:l) = vtmp(1:l)
#endif
                    KAPPA = 0.5d0 * ONE_OVER_H * ddot(l, p, 1, u, 1)

                    p(1:l) = p(1:l) - kappa * u(1:l)

                    call DGER(is(l), l, -1.0d0, ul, 1, p, 1, a, lda)
                    call DGER(is(l), l, -1.0d0, p(me + 1), nproc, u, 1, a, lda)
                    !#if defined (_OPENMP) && defined (__NOOMP)
                    !         call omp_set_num_threads(1)  ! restore the previous threads
                    !#endif

                end if  ! endif on H>safemin
            else
                if(RI(L).eq.ME) then
                    G = A(is(l), I)
                end if

#if defined PARALLEL
                call mpi_bcast(g, 1, MPI_DOUBLE_PRECISION, ri(L), comm, ierr)
#endif
                E(I) = G

            end if

            D(I) = H

        end do

        E(1) = 0.0d0
        D(1) = 0.0d0

        do J = 1, N
            V(1:nrl, J) = 0.0d0
            if(RI(J).eq.ME) THEN
                V(IS(J), J) = 1.0d0
            end if
        end do

        !#if defined (_OPENMP) && defined (__NOOMP)
        !         call omp_set_num_threads(old_threads)  ! restore the previous threads
        !#endif
        do I = 2, N
            L = I - 1
            LLOC = IS(L)
            !
            ! here D(i) if  =/ 0 is such that  |D(i)|> safemin2
            if(D(I).ge.SAFEMIN2) then
                ONE_OVER_H = 1.0d0 / D(I)
                !
                if(lloc > 0) then
                    call DGEMV('T', lloc, l, 1.0d0, v(1, 1), ldv, a(1, i), 1, 0.0d0, p(1), 1)
                else
                    P(1:l) = 0.0d0
                end if
#if defined PARALLEL
                call reduce_base_real_to(L, p, vtmp, comm, -1)
#else
          vtmp(1:l) = p(1:l)
#endif

                if(lloc > 0) then
                    call DGER(lloc, l, -ONE_OVER_H, a(1, i), 1, vtmp, 1, v, ldv)
                end if

            end if

        end do
        !#if defined (_OPENMP) && defined (__NOOMP)
        !         call omp_set_num_threads(1)  ! restore the serial
        !#endif

        do I = 1, N
            U(I) = 0.0d0
            if(RI(I).eq.ME) then
                U(I) = A(IS(I), I)
            end if
        end do
#if defined PARALLEL
        call reduce_base_real_to(n, u, d, comm, -1)
#else
      D(1:N) = U(1:N)
#endif
        deallocate(u, p, vtmp, ul, is, ri, Asav)
        !    if(countstrange.ne.0.and.me.eq.0) &
        !   & write(6,*) ' Warning strange iterations with H <<1 !!! ',countstrange
#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads)  ! restore the previous threads
#endif
        return
    end subroutine ptredv
    !==----------------------------------------------==!
    !==----------------------------------------------==!
    subroutine ptqliv(d, e, n, z, ldz, nrl, mpime, comm, machp)
        !
        ! Modified QL algorithm for CRAY T3E PARALLEL MACHINE
        ! calculate the eigenvectors and eigenvalues of a matrix reduced to
        ! tridiagonal form by PTREDV.
        !
        ! AUTHOR : Carlo Cavazzoni - SISSA 1997
        !          comments and suggestions to : carlo.cavazzoni@cineca.it
        ! Modified: S. Sorella -SISSA in 2008
        ! REFERENCES :
        !
        ! NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
        ! W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
        ! CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
        !
        ! PARALLEL NUMERICAL ALGORITHMS,
        ! T.L. FREEMAN AND C.PHILLIPS,
        ! PRENTICE HALL INTERNATIONAL (1992).
        !
        ! NOTE : the algorithm that finds the eigenvalues is not parallelized
        !        ( it scales as O(N^2) ), I preferred to parallelize only the
        !        updating of the eigenvectors because it is the most costly
        !        part of the algorithm ( it scales as O(N^3) ).
        !        For large matrix in practice all the time is spent in the updating
        !        that in this routine scales linearly with the number of processors,
        !        in fact there is no communication at all.
        !
        !
        !     INPUTS :
        !
        !     D(N)     Diagonal elements of the tridiagonal matrix
        !              this vector is equal on all processors.
        !
        !     E(N)     Subdiagonal elements of the tridiagonal matrix
        !              this vector is equal on all processors.
        !
        !     N        DIMENSION OF THE GLOBAL MATRIX.
        !
        !     NRL      NUMBER OF ROWS OF Z BELONGING TO THE LOCAL PROCESSOR.
        !
        !     LDZ      LEADING DIMENSION OF MATRIX Z.
        !
        !     Z(LDZ,N) Orthogonal transformation that tridiagonalizes the original
        !              matrix A.
        !              The rows of the matrix are distributed among processors
        !              with blocking factor 1.
        !              Example for NPROC = 4 :
        !              ROW | PE
        !              1   | 0
        !              2   | 1
        !              3   | 2
        !              4   | 3
        !              5   | 0
        !              6   | 1
        !              ..  | ..
        !     MACHPREC    : RELATIVE/ABSOLUTE  PRECISION CALCULATION
        !     If MACHPREC>0 ABSOLUTE PRECISION, THE PROGRAM STOPS WHEN
        !     SUBDIAGONAL ELEMENTS ARE LESS THAN MACHPREC X MAXMAT WHERE
        !     MAXMAT= MAX(D(:),E(:)).
        !
        !     IF MACHPREC=0: ORIGINAL NUMERICAL RECIPES ALGORITHM BASED ON
        !                    RELATIVE MACHINE PRECISION (DISCORAGED-->UNSTABLE).
        !     IF MACHPREC<0: THE PROGRAM STOPS WHEN
        !                     E(I) <= |MACHPREC|X ( |D(I)|+|D(I+1)|) IS REACHED
        !                    IT MAY BE STABLE FOR MACHPREC>> RELATIVE MACHINE PRECISION.
        !
        !
        !     OUTPUTS :
        !
        !     Z(LDZ,N) EIGENVECTORS OF THE ORIGINAL MATRIX.
        !              THE Jth COLUMN of Z contains the eigenvectors associated
        !              with the jth eigenvalue.
        !              The eigenvectors are scattered among processors (4PE examp. )
        !              eigenvector | PE
        !               elements   |
        !                 V(1)     | 0
        !                 V(2)     | 1
        !                 V(3)     | 2
        !                 V(4)     | 3
        !                 V(5)     | 0
        !                 V(6)     | 1
        !                 ....       ..
        !
        !     D(N)     Eigenvalues of the original matrix,
        !              this vector is equal on all processors.
        !
        !
        !
        !

        implicit none
#ifdef PARALLEL
    include 'mpif.h'
#endif
        integer, intent(in) :: n, nrl, ldz, mpime, comm
        real*8, intent(in) :: machp
        real*8 :: d(n), e(n)
        real*8 :: z(ldz, n)
        logical :: nogoto, nogotomain

        integer :: i, lm, mm, iter, itersav, mk, k, l, m, ierr, maxit, countuflow&
                &, intbuf(3)
        real*8 :: maxmat, machprec, safemin, maxitf
#ifdef __DDIAG
        real*16  :: b,dd,f,g,p,r,c,s,rs,rn,one
#else
        real*8 :: b, dd, f, g, p, r, c, s, rs, rn, one
#endif
        real*8, external :: DLAMCH
        !     Not to be allocated in the stack
        real*8, allocatable :: cv(:, :)
        real*8, allocatable :: fv1(:)
        real*8, allocatable :: fv2(:)
#if defined (_OPENMP) && defined (__NOOMP)
        integer, external :: omp_get_max_threads
        integer old_threads
        old_threads = omp_get_max_threads()
        call omp_set_num_threads(1)   ! scalar code
#endif


        allocate(cv(2, n))
        allocate(fv1(nrl))
        allocate(fv2(nrl))

        one = 1.0
        maxmat = 0.d0
        do i = 1, n
            maxmat = max(maxmat, abs(d(i)), abs(e(i)))
        end do

        machprec = abs(machp)

        if(machp.gt.0.d0) then
            safemin = machprec**2 * maxmat ! for absolute precision
        else
            safemin = 0.d0
        end if

        maxit = 2000  ! no more than n^3 parallel operations S. Sorella change

        maxitf = 0.d0

        countuflow = 0
        do l = 2, n
            e(l - 1) = e(l)
        end do

        do l = 1, n
            iter = 0
            nogotomain = .true.
            do while(nogotomain) ! QL iterations
#ifdef UNREL_DIAG
! In unreliable network  the precision can be different for different
! processors and it is better that only the master make the nasty calculation.
       if ( mpime == 0 ) then
#endif
                itersav = iter
                m = l
                if(machp.gt.0.d0) then
                    do while(m.lt.n.and.abs(e(m)).gt.safemin)
                        !       absolute precision is more stable
                        !       because it eventually leads to a diagonal  matrix + non diagonal
                        !       with error bounded by machprec.
                        !       Relative precision (as below) does not lead to a predictable end.
                        !       In principle it may never stop (as happens for large matrices).
                        m = m + 1
                    end do
                elseif(machp.eq.0.d0) then  ! Original Numerical Recipes algorithm
                    dd = abs(d(m)) + abs(d(m + 1))
                    do while(m.lt.n.and.abs(e(m)) + dd.ne.dd)
                        m = m + 1
                        dd = abs(d(m)) + abs(d(m + 1))
                    end do
                else  ! relative precision equal to machprec in input.
                    dd = abs(d(m)) + abs(d(m + 1))
                    do while(m.lt.n.and.abs(e(m)).gt.machprec * dd)
                        m = m + 1
                        dd = abs(d(m)) + abs(d(m + 1))
                    end do
                end if

                if (m /= l .and.iter <= maxit) then   ! S. Sorella change
                    if (iter == maxit) then

                        write(6, *) ' Warning too many iterations !!! ', iter

                    end if
                    iter = iter + 1
                    !
                    ! iteration is performed on one processor and results broadcast
                    ! to all others to prevent potential problems if all processors
                    ! do not behave in exactly the same way (even with the same data!)
                    !
                    g = (d(l + 1) - d(l)) / (2.0d0 * e(l))
                    r = pythag(g, one)
                    g = d(m) - d(l) + e(l) / (g + sign(r, g))
                    s = 1.0
                    c = 1.0
                    p = 0.0
                    nogoto = .true.
                    i = m - 1
                    do while(nogoto.and.i.ge.l)
                        f = s * e(i)
                        b = c * e(i)
                        !               r=pythag(f,g)
                        call pythags(f, g, rs, rn)
                        r = rs * rn
                        e(i + 1) = r
                        if(r.gt.safemin) then  ! standard Numerical Recipes protection
                            c = (g / rn) / rs
                            s = (f / rn) / rs
                        else
                            countuflow = countuflow + 1
                            !         continue anyway with the identity
                            !         the identity provides the minimum possible accumulation of roundoff
                            c = 1.0
                            s = 0.0
                        end if

                        g = d(i + 1) - p
                        r = (d(i) - g) * s + 2.0 * c * b
                        p = s * r
                        d(i + 1) = g + p
                        g = c * r - b
                        cv(1, i - l + 1) = c
                        cv(2, i - l + 1) = s
                        i = i - 1
                    end do  ! enddo first do while
                    !
                    d(l) = d(l) - p
                    e(l) = g
                    e(m) = 0.0
                end if  ! endif m =/ l
#ifdef UNREL_DIAG
          end if  ! if mpime eq 0
#endif

#if defined PARALLEL
#ifdef UNREL_DIAG
#ifdef __DOOMP
!$omp barrier
#endif
           if(mpime.eq.0) then
           intbuf(1)=m
           intbuf(2)=l
           intbuf(3)=itersav
           end if
           call mpi_bcast(intbuf,3,MPI_INTEGER,0,comm,ierr)
           mm=intbuf(1)
           lm=intbuf(2)
           itersav=intbuf(3)
           call mpi_barrier(comm,ierr)
#ifdef __DOOMP
!$omp barrier
#endif
           call bcast_real(cv,2*(mm-lm),0,comm)
           call bcast_real(d(lm),mm-lm+1,0,comm)
           call bcast_real(e(lm),mm-lm+1,0,comm)
           call mpi_barrier(comm,ierr)
#ifdef __DOOMP
!$omp barrier
#endif
#else
                mm = m
                lm = l
                call mpi_barrier(comm, ierr)
#ifdef __DOOMP
!$omp barrier
#endif
#endif
#endif
                nogotomain = .false.
                if(itersav.le.maxit.and.mm.ne.lm) then
                    do i = mm - 1, lm, -1
                        do k = 1, nrl
                            fv2(k) = z(k, i + 1)
                        end do
                        do k = 1, nrl
                            fv1(k) = z(k, i)
                        end do
                        c = cv(1, i - lm + 1)
                        s = cv(2, i - lm + 1)
                        do k = 1, nrl
                            z(k, i + 1) = s * fv1(k) + c * fv2(k)
                            z(k, i) = c * fv1(k) - s * fv2(k)
                        end do
                    end do
                    nogotomain = .true.
                end if
            end do  ! enddo main do while
            maxitf = maxitf + iter
        end do

        if(mpime.eq.0) then

            maxitf = maxitf / n

            if(maxitf.ge.30) write(6, *) &
                    &'Warning av. # of QL iterations/eigenvector >= 30, check roundoff ', maxitf

            !     if(countuflow.ne.0) write(6,*) &
            !    &' Warning recovered by underflow ',countuflow,'times'

        end if

        deallocate(cv)
        deallocate(fv1)
        deallocate(fv2)

#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads)  ! restore the previous threads
#endif

        return
    end subroutine ptqliv

    !==----------------------------------------------==!


    subroutine peigsrtv(d, v, ldv, n, nrl, iopt)

        !
        !     This routine sorts eigenvalues and eigenvectors
        !     generated by PTREDV and PTQLIV.
        !
        !     AUTHOR : Carlo Cavazzoni - SISSA 1997
        !              comments and suggestions to : carlo.cavazzoni@cineca.it
        !     Modified: S. Sorella -SISSA in 2008
        !

        implicit none
        integer, intent (in) :: n, ldv, nrl, iopt
        real*8, intent(inout) :: d(n), v(ldv, n)

        integer :: i, j, k
        real*8 :: p

        do 13 i = 1, n - 1
            k = i
            p = d(i)
            if(mod(iopt, 2).eq.0) then
                do j = i + 1, n
                    if(d(j).le.p)then
                        k = j
                        p = d(j)
                    end if
                end do
            else
                do j = i + 1, n
                    if(d(j).ge.p)then
                        k = j
                        p = d(j)
                    end if
                end do
            end if
            if(k.ne.i)then
                d(k) = d(i)
                d(i) = p
                !
                !         Exchange local elements of eigenvectors.
                !
                do j = 1, nrl
                    p = v(j, i)
                    v(j, i) = v(j, k)
                    v(j, k) = p
                end do

            end if
        13    continue
        return
    end subroutine peigsrtv

    !
    !-------------------------------------------------------------------------
    function pythag(a, b)
        implicit none
#ifdef __DDIAG
        real*16 :: a, b, pythag
        real*16 :: absa, absb
#else
        real*8 :: a, b, pythag
        real*8 :: absa, absb
#endif
        absa = abs(a)
        absb = abs(b)
        if(absa.gt.absb)then
            pythag = absa * sqrt(1.0 + (absb / absa)**2)
        else
            if(absb.eq.0.0)then
                pythag = 0.0
            else
                pythag = absb * sqrt(1.0 + (absa / absb)**2)
            end if
        end if
        return
    end function pythag
    !

    subroutine  pythags(a, b, rs, rn)
        implicit none
        !    This subroutine mimics the (numerically unstable) function pythag
        !    pythag(a,b)=rs*rn, but rs is given even when rn=0.d0
#ifdef __DDIAG
      real*16 :: a, b, rs, rn
      real*16 :: absa, absb, root2
      parameter(root2=1.4142135623730950488016887242097_16)
#else
        real*8 :: a, b, rs, rn
        real*8 :: absa, absb, root2
        parameter(root2 = 1.4142135623730950D0)
#endif
        !     1 <=  rs <= sqrt(2)

        absa = abs(a)
        absb = abs(b)
        if(absa.gt.absb)then
            rn = absa
            rs = sqrt(1.0 + (absb / absa)**2)
        else
            rn = absb
            if(absb.eq.0.0)then
                rs = root2  ! assumed absa=absb-->0
            else
                rs = sqrt(1.0 + (absa / absb)**2)
            end if
        end if
        !==----------------------------------------------==!
        return

    end subroutine pythags


    subroutine pdspev_drv(jobz, ap, lda, w, z, ldz, &
            nrl, n, nproc, mpime, comm)
        implicit none
        character, intent(in) :: jobz
        integer, intent(in) :: lda, ldz, nrl, n, nproc, mpime
        INTEGER, INTENT(IN) :: comm
        integer i, j
        real*8 :: ap(lda, *), w(*), z(ldz, *)
        !    Not to be allocated in the stack
        real*8, allocatable :: sd(:)
        !
        if(n < 1) return
        !
        allocate (sd (n))
        call ptredv(ap, lda, w, sd, z, ldz, nrl, n, nproc, mpime, comm, 0.d0)

        call ptqliv(w, sd, n, z, ldz, nrl, mpime, comm, 0.d0)

        deallocate (sd)
        call peigsrtv(w, z, ldz, n, nrl, 0)
        return
    end subroutine pdspev_drv

    subroutine zdspev_drv_ss(jobz, ap, lda, w, z, ldz, &
            nrl, n, nproc, mpime, comm, machp, iopt)
        implicit none
        character, intent(in) :: jobz
        real*8, intent(in) :: machp
        real*8, external :: dlamch
        real*8 :: machprec, safm
        complex*16 cdet, cdetm, cdetn
        integer, intent(in) :: lda, ldz, nrl, n, nproc, mpime
        integer, intent(in) :: comm, iopt
        integer i, j, ind
        complex*16 :: ap(lda, *), z(ldz, *)
        real*8 :: w(*)
        complex*16 :: cost, costtot
        !    Not to be allocated in the stack
        complex*16, allocatable :: zd(:)
        real*8, allocatable :: sd(:)
        !    iopt>1 long output
        !

        if(machp.eq.0.d0) then
            machprec = DLAMCH('E')  ! relative machine precision ~ 1d-16
        else
            machprec = machp
        end if
        safm = dlamch('S')

        !     if(mpime.eq.0) write(6,*) ' Input machine precision safe minimum ',machprec,safm
        !
        if(n < 1) return
        !
        allocate (zd (n), sd(n))
        call ztredv(ap, lda, w, zd, z, ldz, nrl, n, nproc, mpime, comm, machprec)

        !      if(mpime.eq.0) then ! long output
        if(mpime.eq.0.and.iopt.gt.1) then ! long output
            write(6, *) ' after COMPLEX HOUSHOLDER '
            do i = 1, n
                write(6, *) i, w(i), zd(i)
            end do
            write(6, *) ' Trace =', sum(w(1:n))
            cdet = 1.d0
            cdetm = 0.d0
            do i = 1, n
                cdetn = cdet * w(i) - abs(zd(i))**2 * cdetm
                cdetm = cdet
                cdet = cdetn
            end do
            write(6, *) ' Determinant ', cdet
        end if

        cost = (1.d0, 0.d0)
        zd(1) = cost
        sd(1) = 0.d0
        do i = 2, n
            sd(i) = abs(zd(i))
            if(sd(i).gt.safm) then
                zd(i) = cost * zd(i) / sd(i)
            else
                zd(i) = cost
            end if
            cost = zd(i) / abs(zd(i)) ! for numerical stability remains a pure phase
            zd(i) = cost
        end do
        !     scaling eigenvector with appropriate phase, i.e. apply the transform
        !     that leads to a real tridiagonal matrix

        !   Proof: The Housolder gives Q H Q^+  tridiagonal complex H'
        !   Then we applied the transformation U dag H' U^dag , with U = zd(j) delta_ij
        !   The new eigenvectors will be written in terms of Q^+ U^dag.
        do j = 1, n
            z(:, j) = conjg(zd(j)) * z(:, j)
        end do

        call ztqliv(w, sd, n, z, ldz, nrl, mpime, comm, machprec)

        !      if(mpime.eq.0) then ! long output
        if(mpime.eq.0.and.iopt.gt.1) then ! long output
            write(6, *) ' after ztqliv '
            do i = 1, n
                write(6, *) i, w(i), sd(i)
            end do
        end if
        deallocate (sd, zd)

        call zeigsrtv(w, z, ldz, n, nrl, iopt)
        if(mpime.eq.0.and.iopt.gt.1) then ! long output
            write(6, *) ' after zeigsrtv '
            do i = 1, n
                write(6, *) i, w(i)
            end do
        end if
        if(iopt.gt.1) then
            costtot = (0.d0, 0.d0)
            do i = 1, n
                do j = 1, i
                    cost = sum(conjg(z(1:nrl, i)) * z(1:nrl, j))
#ifdef PARALLEL
      call reduce_base_complex(1,cost,comm,0)
#endif
                    costtot = cost + costtot
                    !     if(mpime.eq.0) write(6,*) i,j,cost
                end do
            end do
            if(mpime.eq.0) write(6, *) ' Check orthogonality after zeigsrtv', costtot - n
            !    write(6,*) ' Last eigenvector same phase ?  '
            !    do i=1,n
            !    write(6,*) i,real(abs(z(i,n))),z(i,n)/z(3,n)
            !    enddo
        end if

        return
    end subroutine zdspev_drv_ss

    !==----------------------------------------------==!

    subroutine pdspev_drv_ss(jobz, ap, lda, w, z, ldz, &
            nrl, n, nproc, mpime, comm, machp, iopt)
        implicit none
        character, intent(in) :: jobz
        real*8, intent(in) :: machp
        real*8, external :: dlamch
        real*8 :: machprec, cdet, cdetn, cdetm
        integer, intent(in) :: lda, ldz, nrl, n, nproc, mpime
        integer, intent(in) :: comm, iopt
        integer i, j
        real*8 :: ap(lda, *), w(*), z(ldz, *)
        !    Not to be allocated in the stack
        real*8, allocatable :: sd(:)
        !    iopt>1 long output
        !

        if(machp.eq.0.d0) then
            machprec = DLAMCH('E')  ! relative machine precision ~ 1d-16
        else
            machprec = machp
        end if

        !
        if(n < 1) return
        !
        allocate (sd (n))
        call ptredv(ap, lda, w, sd, z, ldz, nrl, n, nproc, mpime, comm, machprec)

        !      if(mpime.eq.0) then ! long output
        if(mpime.eq.0.and.iopt.gt.1) then ! long output
            write(6, *) ' after HOUSHOLDER '
            do i = 1, n
                write(6, *) i, w(i), sd(i)
            end do
            write(6, *) ' Trace =', sum(w(1:n))
            cdet = 1.d0
            cdetm = 0.d0
            do i = 1, n
                cdetn = cdet * w(i) - abs(sd(i))**2 * cdetm
                cdetm = cdet
                cdet = cdetn
            end do
            write(6, *) ' Determinant ', cdet
        end if

        call ptqliv(w, sd, n, z, ldz, nrl, mpime, comm, machprec)

        !      if(mpime.eq.0) then ! long output
        if(mpime.eq.0.and.iopt.gt.1) then ! long output
            write(6, *) ' after ptqliv '
            do i = 1, n
                write(6, *) i, w(i), sd(i)
            end do
        end if

        deallocate (sd)
        call peigsrtv(w, z, ldz, n, nrl, iopt)
        if(mpime.eq.0.and.iopt.gt.1) then ! long output
            write(6, *) ' after peigsrtv '
            cdetn = 1.d0
            do i = 1, n
                write(6, *) i, w(i)
                cdetn = cdetn * w(i)
            end do
            write(6, *) ' Determinant new =', cdetn, cdet
            !      write(6,*) 'Last eigenvector '
            !      do i=1,n
            !      write(6,*) i,z(i,20)
            !      enddo
        end if
        RETURN
    end subroutine pdspev_drv_ss

    subroutine dspev_drv(jobz, uplo, n, ap, w, z, ldz)
        implicit none
        character :: jobz, uplo
        integer :: iopt, info, ldz, n
        real*8 :: ap(*), w(*), z(ldz, *)
        real*8, allocatable :: work(:)

        if(n < 1) return

        allocate(work(3 * n))

#if defined __ESSL
        IOPT = 0
        if((JOBZ .eq. 'V') .or. (JOBZ .eq. 'v')) iopt = iopt + 1
        if((UPLO .eq. 'U') .or. (UPLO .eq. 'u')) iopt = iopt + 20
        call dspev(iopt, ap, w, z, ldz, n, work, 3 * n)
#else
          call dspev(jobz, uplo, n, ap(1), w(1), z(1,1), ldz, work, info)
          if( info .ne. 0 ) THEN
            call errore( ' dspev_drv ', ' diagonalization failed ',info )
          end if
#endif

        deallocate(work)

        return
    end subroutine dspev_drv


    subroutine ztredv(a, lda, d, e, v, ldv, nrl, n, nproc, me, comm, machp)
        !
        !     Parallel version of the famous HOUSEHOLDER tridiagonalization
        !     Algorithm for simmetric matrix.
        !
        !     AUTHOR : Carlo Cavazzoni - SISSA 1997
        !     modified by S. Sorella on January 2008
        !              comments and suggestions to : carlo.cavazzoni@cineca.it
        !
        ! REFERENCES :
        !
        ! NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
        ! W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
        ! CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
        !
        ! PARALLEL NUMERICAL ALGORITHMS,
        ! T.L. FREEMAN AND C.PHILLIPS,
        ! PRENTICE HALL INTERNATIONAL (1992).
        !
        !
        !
        !     INPUTS :
        !
        !     A(NRL,N) Local part of the global matrix A(N,N) to be reduced,
        !              only the upper triangle is needed.
        !              The rows of the matrix are distributed among processors
        !              with blocking factor 1.
        !              Example for NPROC = 4 :
        !              ROW | PE
        !              1   | 0
        !              2   | 1
        !              3   | 2
        !              4   | 3
        !              5   | 0
        !              6   | 1
        !              ..  | ..
        !
        !     LDA      LEADING DIMENSION OF MATRIX A.
        !
        !     LDV      LEADING DIMENSION OF MATRIX V.
        !
        !     NRL      NUMBER OF ROWS BELONGING TO THE LOCAL PROCESSOR.
        !
        !     N        DIMENSION OF THE GLOBAL MATRIX.
        !
        !     NPROC    NUMBER OF PROCESSORS.
        !
        !     ME       INDEX OF THE LOCAL PROCESSOR (Starting from 0).
    !
    !
    !     OUTPUTS :
    !
    !     V(NRL,N) Orthogonal transformation that tridiagonalize A,
    !              this matrix is distributed among processor
    !              in the same way as A.
    !
    !     D(N)     Diagonal elements of the tridiagonal matrix
    !              this vector is equal on all processors.
    !
    !     E(N)     Subdiagonal elements of the tridiagonal matrix
    !              this vector is equal on all processors.
    !
    !
    !     machprec Input relative precision of eigenvalues if >=0. Otherwise
    !     lapack machine relative precision obtained by DLAMCH('E') is used.
    !     It may be useful to have less converged eigenvectors/eigenvalues
    !     with a cheap computation for large matrices.

    implicit none
#ifdef PARALLEL
include 'mpif.h'
#endif
    integer, intent(in) :: N, NRL, LDA, LDV
    integer, intent(in) :: NPROC, ME, comm
    complex*16 :: A(LDA, N), E(N), V(LDV, N)
    real*8 :: D(N)
    real*8, intent(in) :: machp
    !
    complex*16, external :: zdotc_
    !
    real*8 :: scalef, sigma, h, tmp, machprec, safemin, maxmat&
                &, safemin2, safm, kappa
        complex*16 :: f, g, ctmp, czero, cone, trace
        real*8, external :: DLAMCH
        complex*16, allocatable :: u(:)
        complex*16, allocatable :: p(:)
        complex*16, allocatable :: vtmp(:)

        real*8 :: one_over_h
        real*8 :: one_over_scale, cost
        complex*16, allocatable :: ul(:)
        complex*16, allocatable :: asav(:)
        integer :: l, i, j, ii, k, t, tl, ierr, countstrange
        integer :: kl, jl, ks, lloc
        integer, allocatable :: is(:)
        integer, allocatable :: ri(:)
#if defined (_OPENMP) && defined (__NOOMP)
        integer, external :: omp_get_max_threads
        integer old_threads
        old_threads = omp_get_max_threads()
        call omp_set_num_threads(1)   ! scalar code
#endif

        czero = (0.d0, 0.d0)
        cone = (1.d0, 0.d0)

        !     a(:,:)=real(a(:,:))
        !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........

        if(N == 0) then
            return
        end if

        countstrange = 0

        !     The allocation of p below is necessary when nproc> n otherwise some
        !     call (DGER below) uses undefined memory (possible abort).
        allocate(u(n + 2), p(max(n + 1, nproc)), vtmp(n + 2), ul(n)&
                &, is(n), ri(n), Asav(lda))
        u = czero
        p = czero
        vtmp = czero
        ul = czero
        is = 0
        ri = 0
        asav = czero

        do I = N, 1, -1
            IS(I) = (I - 1) / NPROC
            RI(I) = MOD((I - 1), NPROC)  !  owner of I-th row
            if(ME .le. RI(I)) then
                IS(I) = IS(I) + 1
            end if
        end do
#ifdef PARALLEL
      tmp=0.d0
      do i=1,N
          do k=1,is(i)
              tmp=max(tmp,abs(A(k,i)))
          end do
      end do
!     The diagonal part should be real for stability/consistency reasons
!     ctmp=(0.d0,0.d0)
      do i=1,N
      if(me.eq.ri(i)) then
!     ctmp=ctmp+A(is(i),i)
      A(is(i),i)=real(A(is(i),i))
      end if
      end do

      call mpi_allreduce(tmp,maxmat,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM,ierr)
!     call mpi_allreduce(ctmp,trace,1,MPI_COMPLEX16,MPI_SUM,COMM,ierr)
!     if(me.eq.0) write(6,*) ' Trace input ',trace
#else
        maxmat = 0.d0
        do i = 1, N
            do k = 1, i
                maxmat = max(maxmat, abs(A(k, i)))
            end do
            a(i, i) = real(a(i, i))
        end do
#endif

        if(maxmat.eq.0.d0) then
            D = 0.d0
            E = czero
            return
        end if
        machprec = abs(machp)
        safemin = machprec**2 * maxmat
        safm = dlamch('S')
        safemin2 = dlamch('S') * 1000000.d0 * N  ! safer minimum than lapack

        do I = N, 2, -1

            L = I - 1         ! first element
            H = 0.0d0

            if (L > 1) then

                SCALEF = 0.0d0
                do K = 1, is(l)
                    SCALEF = SCALEF + ABS(A(K, I))
                end do

#if defined PARALLEL
                tmp = scalef
                call mpi_allreduce(tmp, scalef, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, ierr)
#endif

                if (SCALEF .le. SAFEMIN)  then
                    !
                    if (RI(L).eq.ME) then
                        ctmp = A(is(L), I)
                    else
                        ctmp = czero
                    end if
#ifdef PARALLEL
     call mpi_allreduce(ctmp,E(I),1,MPI_COMPLEX16,MPI_SUM,COMM,ierr)
#else
                    E(I) = ctmp
#endif
                else

                    !  ......  CALCULATION OF SIGMA AND H

                    ONE_OVER_SCALE = 1.0d0 / SCALEF
                    SIGMA = 0.0d0
                    do k = 1, is(L)
                        asav(k) = A(k, i)
                        A(k, I) = A(k, I) * ONE_OVER_SCALE
                        SIGMA = SIGMA + A(k, I) * DCONJG(A(k, I))
                    end do

                    if(ri(l) .eq. me) then
                        F = A(is(l), i)
                    else
                        F = czero
                    end if

                    !  CONSTRUCTION OF VECTOR U

                    vtmp(1:l) = czero

                    k = ME + 1
                    do kl = 1, is(l)
                        vtmp(k) = A(kl, I)
                        k = k + NPROC
                    end do

                    do kl = 1, is(l)
                        UL(kl) = A(kl, I)
                    end do

#if defined PARALLEL
                    vtmp(l + 1) = sigma
                    vtmp(l + 2) = f
                    call reduce_base_complex_to(L + 2, vtmp, u, comm, -1)
                    sigma = u(l + 1)
                    f = u(l + 2)
#else
             u(1:l) = vtmp(1:l)
#endif

                    !             G          = -SIGN(SQRT(SIGMA),F)

                    cost = abs(f)
                    if(cost.gt.safm) then
                        G = -SQRT(SIGMA) * (F / cost)
                        H = SIGMA - DCONJG(F) * G
                    else
                        G = -SQRT(SIGMA)
                        H = SIGMA
                    end if

                end if

                if(H.lt.SAFEMIN2) then ! S.S. Added another control: in principle
                    ! H can be very small so that 1/H may overflow.
                    countstrange = countstrange + 1

                    !  restore the previous unscaled A
                    do k = 1, is(L)
                        A(k, I) = Asav(k)
                    end do
                    !            I am not sure E(I) is global in this case so mpi_allreduce needed
                    if (RI(L).eq.ME) THEN
                        ctmp = A(is(L), I)
                    else
                        ctmp = czero
                    end if
#ifdef PARALLEL
     call mpi_allreduce(ctmp,E(I),1,MPI_COMPLEX16,MPI_SUM,COMM,ierr)
#else
                    E(I) = ctmp
#endif

                    !            H=0.d0  do not make useless approximations.

                else
                    ONE_OVER_H = 1.0d0 / H
                    E(I) = SCALEF * G
                    U(L) = F - G
                    if(RI(L) ==  ME) then
                        UL(is(l)) = F - G
                        A(is(l), I) = F - G
                    end if
                    !  CONSTRUCTION OF VECTOR P
#ifdef __DOOMP
!$omp parallel do shared(vtmp,L,IS,A,UL,RI,U,ME,ONE_OVER_H)  private(J,KL,K)
#endif
                    do J = 1, L
                        vtmp(j) = CZERO
                        do KL = 1, IS(J)
                            vtmp(J) = vtmp(J) + dconjg(A(KL, J)) * UL(KL)
                        end do
                        if(L > J .and. ME == RI(J)) then
                            do K = J + 1, L
                                vtmp(J) = vtmp(J) + A(IS(J), K) * U(K)
                            end do
                        end if
                        vtmp(J) = vtmp(J) * ONE_OVER_H
                    end do
#ifdef __DOOMP
!$omp end parallel do
#endif

#ifdef PARALLEL
             call reduce_base_complex_to( L, vtmp, p, comm, -1 )
#else
                    p(1:l) = vtmp(1:l)
#endif
                    !            KAPPA = 0.5d0 * ONE_OVER_H * ddot( l, p, 1, u, 1 )
                    !            KAPPA is real
                    KAPPA = 0.5d0 * ONE_OVER_H * zdotc_(l, u, 1, p, 1)

                    p(1:l) = p(1:l) - kappa * u(1:l)

                    call ZGERC(is(l), l, -CONE, ul, 1, p, 1, a, lda)
                    call ZGERC(is(l), l, -CONE, p(me + 1), nproc, u, 1, a, lda)

                    !            CALL DGER( is(l), l, -1.0d0, ul, 1, p, 1, a, lda )
                    !            CALL DGER( is(l), l, -1.0d0, p( me + 1 ), nproc, u, 1, a, lda )
                    !#if defined (_OPENMP) && defined (__NOOMP)
                    !         call omp_set_num_threads(1)  ! restore the previous threads
                    !#endif

                end if  ! endif on H>safemin

            else

                if(RI(L).eq.ME) THEN
                    G = A(is(l), I)
                end if

#if defined PARALLEL
                call mpi_bcast(g, 1, MPI_COMPLEX16, ri(L), comm, ierr)
#endif
                E(I) = G

            end if

            D(I) = H

        end do

        E(1) = CZERO
        D(1) = 0.0d0

        !     The identity

        do J = 1, N
            V(1:nrl, J) = CZERO
            if(RI(J).eq.ME) THEN
                V(IS(J), J) = CONE
            end if
        end do

        !#if defined (_OPENMP) && defined (__NOOMP)
        !         call omp_set_num_threads(old_threads)  ! restore the previous threads
        !#endif
        do I = 2, N
            !      a(*,i) contains the vector used in the Householder elimination

            L = I - 1
            LLOC = IS(L)
            !
            ! here D(i) if  =/ 0 is such that  |D(i)|> safemin2
            if(D(I).ge.SAFEMIN2) then
                ONE_OVER_H = 1.0d0 / D(I)
                !
                if(lloc > 0) then
                    !     CALL DGEMV( 'T', lloc, l, 1.0d0, v(1,1), ldv, a(1,i), 1, 0.0d0, p(1), 1 )
                    call ZGEMV('C', lloc, l, CONE, v, ldv, a(1, i), 1, CZERO, p(1), 1)
                else
                    P(1:l) = CZERO
                end if
#if defined PARALLEL
                call reduce_base_complex_to(L, p, vtmp, comm, -1)
#else
          vtmp(1:l) = p(1:l)
#endif

                if(lloc > 0) then
                    !         CALL DGER( lloc, l, -ONE_OVER_H, a(1,i), 1, vtmp, 1, v, ldv )
                    ctmp = -ONE_OVER_H
                    call ZGERC(lloc, l, ctmp, a(1, i), 1, vtmp, 1, v, ldv)
                end if
            end if

        end do
        !#if defined (_OPENMP) && defined (__NOOMP)
        !         call omp_set_num_threads(1)  ! restore the serial
        !#endif

        do I = 1, N
            D(I) = 0.d0
            if(RI(I).eq.ME) then
                D(I) = A(IS(I), I)
            end if
        end do
#if defined PARALLEL
        call reduce_base_real(n, d, comm, -1)
#endif
        deallocate(u, p, vtmp, ul, is, ri, asav)
        !    if(countstrange.ne.0.and.me.eq.0) &
        !   & write(6,*) ' Warning strange iterations with H <<1 !!! ',countstrange
#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads)  ! restore the previous threads
#endif
        return
    end subroutine ztredv

    !==----------------------------------------------==!

    subroutine zeigsrtv(d, v, ldv, n, nrl, iopt)

        !
        !     This routine sorts eigenvalues and eigenvectors
        !     generated by PTREDV and PTQLIV.
        !
        !     AUTHOR : Carlo Cavazzoni - SISSA 1997
        !              comments and suggestions to : carlo.cavazzoni@cineca.it
        !     Modified: S. Sorella -SISSA in 2008
        !     if iopt is odd they are ordered from maximum (first element) to minimum
        !         =      even  they are ordered in the standard from minimum to max.

        implicit none
        integer, intent (in) :: n, ldv, nrl, iopt
        real*8, intent(inout) :: d(n)
        complex*16 v(ldv, n)

        integer :: i, j, k
        real*8 :: p
        complex*16 :: cp

        do 13 i = 1, n - 1
            k = i
            p = d(i)
            if(mod(iopt, 2).eq.0) then
                do j = i + 1, n
                    if(d(j).le.p)then
                        k = j
                        p = d(j)
                    end if
                end do
            else
                do j = i + 1, n
                    if(d(j).ge.p)then
                        k = j
                        p = d(j)
                    end if
                end do
            end if
            if(k.ne.i)then
                d(k) = d(i)
                d(i) = p
                !
                !         Exchange local elements of eigenvectors.
                !
                do j = 1, nrl
                    cp = v(j, i)
                    v(j, i) = v(j, k)
                    v(j, k) = cp
                end do
            end if
        13    continue
        return
    end subroutine zeigsrtv


    subroutine ztqliv(d, e, n, z, ldz, nrl, mpime, comm, machp)
        !
        ! Modified QL algorithm for CRAY T3E PARALLEL MACHINE
        ! calculate the eigenvectors and eigenvalues of a matrix reduced to
        ! tridiagonal form by PTREDV.
        !
        ! AUTHOR : Carlo Cavazzoni - SISSA 1997
        !          comments and suggestions to : carlo.cavazzoni@cineca.it
        ! Modified: S. Sorella -SISSA in 2008
        ! REFERENCES :
        !
        ! NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
        ! W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
        ! CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
        !
        ! PARALLEL NUMERICAL ALGORITHMS,
        ! T.L. FREEMAN AND C.PHILLIPS,
        ! PRENTICE HALL INTERNATIONAL (1992).
        !
        ! NOTE : the algorithm that finds the eigenvalues is not parallelized
        !        ( it scales as O(N^2) ), I preferred to parallelize only the
        !        updating of the eigenvectors because it is the most costly
        !        part of the algorithm ( it scales as O(N^3) ).
        !        For large matrix in practice all the time is spent in the updating
        !        that in this routine scales linearly with the number of processors,
        !        in fact there is no communication at all.
        !
        !
        !     INPUTS :
        !
        !     D(N)     Diagonal elements of the tridiagonal matrix
        !              this vector is equal on all processors.
        !
        !     E(N)     Subdiagonal elements of the tridiagonal matrix
        !              this vector is equal on all processors.
        !
        !     N        DIMENSION OF THE GLOBAL MATRIX.
        !
        !     NRL      NUMBER OF ROWS OF Z BELONGING TO THE LOCAL PROCESSOR.
        !
        !     LDZ      LEADING DIMENSION OF MATRIX Z.
        !
        !     Z(LDZ,N) Orthogonal transformation that tridiagonalizes the original
        !              matrix A.
        !              The rows of the matrix are distributed among processors
        !              with blocking factor 1.
        !              Example for NPROC = 4 :
        !              ROW | PE
        !              1   | 0
        !              2   | 1
        !              3   | 2
        !              4   | 3
        !              5   | 0
        !              6   | 1
        !              ..  | ..
        !     MACHPREC    : RELATIVE/ABSOLUTE  PRECISION CALCULATION
        !     If MACHPREC>0 ABSOLUTE PRECISION, THE PROGRAM STOPS WHEN
        !     SUBDIAGONAL ELEMENTS ARE LESS THAN MACHPREC X MAXMAT WHERE
        !     MAXMAT= MAX(D(:),E(:)).
        !
        !     IF MACHPREC=0: ORIGINAL NUMERICAL RECIPES ALGORITHM BASED ON
        !                    RELATIVE MACHINE PRECISION (DISCORAGED-->UNSTABLE).
        !     IF MACHPREC<0: THE PROGRAM STOPS WHEN
        !                     E(I) <= |MACHPREC|X ( |D(I)|+|D(I+1)|) IS REACHED
        !                    IT MAY BE STABLE FOR MACHPREC>> RELATIVE MACHINE PRECISION.
        !
        !
        !     OUTPUTS :
        !
        !     Z(LDZ,N) EIGENVECTORS OF THE ORIGINAL MATRIX.
        !              THE Jth COLUMN of Z contains the eigenvectors associated
        !              with the jth eigenvalue.
        !              The eigenvectors are scattered among processors (4PE examp. )
        !              eigenvector | PE
        !               elements   |
        !                 V(1)     | 0
        !                 V(2)     | 1
        !                 V(3)     | 2
        !                 V(4)     | 3
        !                 V(5)     | 0
        !                 V(6)     | 1
        !                 ....       ..
        !
        !     D(N)     Eigenvalues of the original matrix,
        !              this vector is equal on all processors.
        !
        !
        !
        !

        implicit none
#ifdef PARALLEL
    include 'mpif.h'
#endif
        integer, intent(in) :: n, nrl, ldz, mpime, comm
        real*8, intent(in) :: machp
        real*8 :: d(n), e(n)
        complex*16 :: z(ldz, n)
        logical :: nogoto, nogotomain

        integer :: i, lm, mm, iter, itersav, mk, k, l, m, ierr, maxit, countuflow&
                &, intbuf(3)
        real*8 :: maxmat, machprec, safemin, maxitf
#ifdef __DDIAG
        real*16  :: b,dd,f,g,p,r,c,s,rs,rn,one
#else
        real*8 :: b, dd, f, g, p, r, c, s, rs, rn, one
#endif
        real*8, external :: DLAMCH
        !     Not to be allocated in the stack
        real*8, allocatable :: cv(:, :)
        complex*16, allocatable :: fv1(:)
        complex*16, allocatable :: fv2(:)
#if defined (_OPENMP) && defined (__NOOMP)
        integer, external :: omp_get_max_threads
        integer old_threads
        old_threads = omp_get_max_threads()
        call omp_set_num_threads(1)   ! scalar code
#endif


        allocate(cv(2, N))
        allocate(fv1(nrl))
        allocate(fv2(nrl))

        one = 1.0
        maxmat = 0.d0
        do i = 1, n
            maxmat = max(maxmat, abs(d(i)), abs(e(i)))
        end do

        machprec = abs(machp)

        if(machp.gt.0.d0) then
            safemin = machprec**2 * maxmat ! for absolute precision
        else
            safemin = 0.d0
        end if

        maxit = 2000  ! no more than n^3 parallel operations S. Sorella change

        maxitf = 0.d0

        countuflow = 0
        do l = 2, n
            e(l - 1) = e(l)
        end do

        do l = 1, n
            iter = 0
            nogotomain = .true.
            do while(nogotomain) ! QL iterations
#ifdef UNREL_DIAG
! In unreliable network  the precision can be different for different
! processors and it is better that only the master make the nasty calculation.
       if ( mpime == 0 ) then
#endif
                itersav = iter
                m = l
                if(machp.gt.0.d0) then
                    do while(m.lt.n.and.abs(e(m)).gt.safemin)
                        !       absolute precision is more stable
                        !       because it eventually leads to a diagonal  matrix + non diagonal
                        !       with error bounded by machprec.
                        !       Relative precision (as below) does not lead to a predictable end.
                        !       In principle it may never stop (as happens for large matrices).
                        m = m + 1
                    end do
                elseif(machp.eq.0.d0) then  ! Original Numerical Recipes algorithm
                    dd = abs(d(m)) + abs(d(m + 1))
                    do while(m.lt.n.and.abs(e(m)) + dd.ne.dd)
                        m = m + 1
                        dd = abs(d(m)) + abs(d(m + 1))
                    end do
                else  ! relative precision equal to machprec in input.
                    dd = abs(d(m)) + abs(d(m + 1))
                    do while(m.lt.n.and.abs(e(m)).gt.machprec * dd)
                        m = m + 1
                        dd = abs(d(m)) + abs(d(m + 1))
                    end do
                end if

                if (m /= l .and.iter <= maxit) then   ! S. Sorella change
                    if (iter == maxit) then

                        write(6, *) ' Warning too many iterations !!! ', iter

                    end if
                    iter = iter + 1
                    !
                    ! iteration is performed on one processor and results broadcast
                    ! to all others to prevent potential problems if all processors
                    ! do not behave in exactly the same way (even with the same data!)
                    !
                    g = (d(l + 1) - d(l)) / (2.0d0 * e(l))
                    r = pythag(g, one)
                    g = d(m) - d(l) + e(l) / (g + sign(r, g))
                    s = 1.0
                    c = 1.0
                    p = 0.0
                    nogoto = .true.
                    i = m - 1
                    do while(nogoto.and.i.ge.l)
                        f = s * e(i)
                        b = c * e(i)
                        !               r=pythag(f,g)
                        call pythags(f, g, rs, rn)
                        r = rs * rn
                        e(i + 1) = r
                        if(r.gt.safemin) then  ! standard Numerical Recipes protection
                            c = (g / rn) / rs
                            s = (f / rn) / rs
                        else
                            countuflow = countuflow + 1
                            !         continue anyway with the identity
                            !         the identity provides the minimum possible accumulation of roundoff
                            c = 1.0
                            s = 0.0
                        end if

                        g = d(i + 1) - p
                        r = (d(i) - g) * s + 2.0 * c * b
                        p = s * r
                        d(i + 1) = g + p
                        g = c * r - b
                        cv(1, i - l + 1) = c
                        cv(2, i - l + 1) = s
                        i = i - 1
                    end do  ! enddo first do while
                    !
                    d(l) = d(l) - p
                    e(l) = g
                    e(m) = 0.0
                end if  ! endif m =/ l
#ifdef UNREL_DIAG
          end if  ! if mpime eq 0
#endif

#if defined PARALLEL
#ifdef UNREL_DIAG
#ifdef __DOOMP
!$omp barrier
#endif
           if(mpime.eq.0) then
           intbuf(1)=m
           intbuf(2)=l
           intbuf(3)=itersav
           end if
           call mpi_bcast(intbuf,3,MPI_INTEGER,0,comm,ierr)
           mm=intbuf(1)
           lm=intbuf(2)
           itersav=intbuf(3)
           call mpi_barrier(comm,ierr)
#ifdef __DOOMP
!$omp barrier
#endif
           call bcast_real(cv,2*(mm-lm),0,comm)
           call bcast_real(d(lm),mm-lm+1,0,comm)
           call bcast_real(e(lm),mm-lm+1,0,comm)
           call mpi_barrier(comm,ierr)
#ifdef __DOOMP
!$omp barrier
#endif
#else
                mm = m
                lm = l
                call mpi_barrier(comm, ierr)
#ifdef __DOOMP
!$omp barrier
#endif
#endif
#else
    mm=m
    lm=l
#endif


                nogotomain = .false.
                if(itersav.le.maxit.and.mm.ne.lm) then

                    do i = mm - 1, lm, -1
                        do k = 1, nrl
                            fv2(k) = z(k, i + 1)
                        end do
                        do k = 1, nrl
                            fv1(k) = z(k, i)
                        end do
                        c = cv(1, i - lm + 1)
                        s = cv(2, i - lm + 1)
                        do k = 1, nrl
                            z(k, i + 1) = s * fv1(k) + c * fv2(k)
                            z(k, i) = c * fv1(k) - s * fv2(k)
                        end do
                    end do
                    nogotomain = .true.
                end if
            end do  ! enddo main do while
            maxitf = maxitf + iter
        end do

        if(mpime.eq.0) then

            maxitf = maxitf / n

            if(maxitf.ge.30) write(6, *) &
                    &'Warning av. # of QL iterations/eigenvector >= 30, check roundoff ', maxitf

            !     if(countuflow.ne.0) write(6,*) &
            !    &' Warning recovered by underflow ',countuflow,'times'

        end if

        deallocate(cv)
        deallocate(fv1)
        deallocate(fv2)

#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads)  ! restore the previous threads
#endif

        return
    end subroutine ztqliv

end module dspev_module
