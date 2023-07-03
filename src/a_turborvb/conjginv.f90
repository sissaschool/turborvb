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

subroutine conjginv(np, kp_complex, symmagp, nbin, rank_t, rank, comm_mpi, mat, forza&
        &, g, h, psip, maxit, eps, epsdgel, fkav, x, parcut, eps_umrigar)
    implicit none
    integer np, npk, i, ii, j, k, nl, npp, nbin, maxit, iter, rank, info, np3, np4&
            &, npm, ierr, countzero, comm_mpi, kpr, kpc, kp_complex, kpcp, npt, kpp, rank_t
    real*8 mat(nbin, np), h(np), g(np), psip(*), eps                       &
            &, epsdgel, sv, eigmin, z, error, dnrm2, fkav(*), forza(*), cost            &
            &, lambda, ng2, ddot, gamma, Fun, errorg, x(np), parcut, maxsr, eps_umrigar&
            &, cost_umrigar
    logical symmagp, yes_umrigar
#ifdef PARALLEL
    include 'mpif.h'
    call mpi_barrier(comm_mpi, ierr)
#endif
    !      input mat(i,j)=  \sqrt( w_j)  (O_i(x_j)-\bar O_i)/sqrt(sum_k w_k)
    !      modified (normalized) on output
    !      input forza(i) i the force acting on the parameter i
    !      fkav(i)=sqrt( s_{i,i}) where:
    !      s_{i,j} = sum_k  mat_{i,k} mat_{j,k}
    !      Output  x  solution of the linear system  s x =force
    !      within a tollerance eps
    !      For QMC convenience (avoid too much noise in the inversion)
    !      The matrix s is regularized by s_{i,i} ---> s_{i,i} (1 + epsdgel)
    !----------------------------------------------------------------------------
    !
    !       COMPLEX CASE
    !       In this case the complex parameters are at the beginning.
    !      For each variational parameter (both real and imaginary part if complex)
    !      we have the real and imaginary part corresponding to the operator.
    !   e.g. alpha_k (always real, even when it correspond to an imaginary part)
    !        --> Ok = dlog Psi^*/d alpha_k
    !      There are two cases:
    !   1) symmagp = .true.
    !     The parameter value (real or imaginary) alpha_k Ok
    !       occupy only odd position h(1),h(3)... the other assumed zero
    !       the operator value (in the matrix) both.
    !      the SR matrix in this case is    Re (< O_k^* O_k'>)
    !   2) symmagp=.false.
    !      The real part of the parameter occupy odd position h(1),h(3),...
    !      The imaginary part the even position h(2),h(4),....
    !      The operator O(1) O(2) means the real and -imaginary part of alpha_1 (real)
    !  O = dlog Psi^*/dalphak  alpha_k (complex)= Re alpha_k, Im alpha_k both real
    !        O(1)=dlog Re Psi/dRe alpha_1 , O(2)=-dIm log Psi/dRe alpha_1
    !       Im  alpha_1 -->    dlog Re Psi/dIm alpha_1 , -dIm log Psi/dIm alpha_1
    !      By the Cauchy relations the ones corresponding to the imaginary part are
    !         dlog Re log Psi/dIm alpha= -dlog Im log Psi/dRe alpha
    !         dlog Im log Psi/dIm alpha= dlog Re log Psi/dRe alpha
    !      we get:
    !      Im alpha -->  O(2)-i O(1)
    !       Then the matrix SR between complex parameters read:
    !        real real:
    !        SR = Re  <O(2k-1)+i O(2k)| O(2k'-1) +i O(2k') >=
    !        < O(2k-1) O(2k'-1)>+ <O(2k) O(2k')>
    !        real imag:
    !        SR = Re  < O(2k-1)+i O(2k) | O(2k')-i O(2k'-1)>
    !        < O(2k-1) O(2k')>- <O(2k) O(2k'-1)>
    !        imag real:
    !        SR = Re  < O(2k)-i O(2k-1) | O(2k'-1)+i O(2k')>
    !        < O(2k) O(2k'-1)>- <O(2k-1) O(2k')>
    !        imag imag:
    !        SR = Re  < O(2k-1)-i O(2k) | O(2k'-1)-i O(2k')>
    !        < O(2k-1) O(2k'-1)>+ <O(2k) O(2k')>

    !        Instead between a complex and real
    !        real of complex, real
    !        SR = Re < O(2k-1) +i O(2k) | O(k')> = < O(2k-1) | O(k')>
    !        Imag of complex, real
    !        SR = Re < O(2k)+i O(2k-1) | O(k')> = < O(2k)| O(k')>
    !        real, real of complex
    !        SR = Re < O(k) | O(2k'-1)+iO(2k')> = < O(k)| O(2k'-1)>
    !        real, imag of complex
    !        SR = Re < O(k) | O(2k')+iO(2k'-1)> = < O(k)| O(2k')>
    !        nothing change between real and real SR=  < O(k)| O(k')>.
    !  Now h has three parts hr (h(1),h(3)...), hi (h(2), h(4),...) and hn
    !   real parameters for i>kp_complex
    !   The loaded matrix is
    !          | Srr Sri Srn |
    !          | Sir Sii Sin |   = mat^T x mat    (is stored in compact form in mat)
    !          | Snr Sni Snn |
    !    And one has to apply a matrix vector of the type (according to above):
    !          | Srr+Sii  Sri-Sir Srn | hr | =
    !          | Sir-Sri  Sii+Srr Sin | hi | =
    !          | Snr      Sni    Snn | hn | =
    !          | Srr Sri Srn | hr |   | Sii  -Sir 0| hr |
    !          | Sir Sii Sin | hi | + | -Sri  Srr 0| hi |
    !          | Snr Sni Snn | hn |   |  0     0  0| 0  |

    !                 = mat^T (mat h) + | -v_2  |
    !                                   | v_1  |
    !                                   |  0   |
    !                 where
    !                | v_1 |    =   | Srr   Sri | hi |=  Srr hi - Sri hr
    !                | v_2 |        | Sir   Sii | -hr |=  Sir hi - Sii hr
    !This is an ARROCCO of the stored matrix limited only to the complex parameters.
    !     ARROCCO=| 0 I|
    !             |-I 0|
    !      So that the added contribution is written as ARROCCO^T  mat^T mat ARROCCO
    !      that remains symmetric as it should.
    npk = np - kp_complex ! the number of  real parameter
    kpr = kp_complex + 1 ! the first real parameter
    kpc = kp_complex/2
    do i = 1, kp_complex, 2
        cost = sum(mat(1:nbin, i)**2) + sum(mat(1:nbin, i + 1)**2)
        if (fkav(i) .ne. 0.d0) fkav(i) = cost
        if (fkav(i + 1) .ne. 0.d0) fkav(i + 1) = cost
    end do
    do i = kpr, np
        if (fkav(i) .ne. 0.d0) fkav(i) = sum(mat(1:nbin, i)**2)
    end do

#ifdef PARALLEL
    call reduce_base_real(np, fkav, comm_mpi, -1)
#endif
    !       if(rank.eq.0) then
    !       write(6,*) ' fkav inside '
    !       do i=1,np
    !       write(6,*) i,fkav(i)
    !       enddo
    !       endif
    if (eps_umrigar .ne. 0.d0) then
        yes_umrigar = .true.
        cost_umrigar = eps_umrigar
        do i = 1, np
            if (fkav(i) .gt. 0.d0) then
                fkav(i) = dsqrt(fkav(i) + eps_umrigar)
            else
                fkav(i) = 0.d0
            end if
        end do
    else
        cost_umrigar = 0.d0
        yes_umrigar = .false.
        do i = 1, np
            if (fkav(i) .gt. 0.d0) then
                fkav(i) = dsqrt(fkav(i))
            else
                fkav(i) = 0.d0
            end if
        end do
    end if
    !      call mpi_finalize(ierr)
    !      stop
    !      rescaling mat
#ifdef UNREL_DIAG
! Different processor may have different fkav due to different
! accuracy in sqrt
    call bcast_real(fkav, np, 0, comm_mpi)
#endif

    iter = 0
    error = 2.d0*eps
    npp = np + 1
    if (kpc .ne. 0) then
        kpp = kp_complex + 1
        npt = npk + kpc
        kpcp = npp + kpc
    end if

    maxsr = 0.d0
    do i = 1, kp_complex, 2
        if (fkav(i) .gt. maxsr) maxsr = fkav(i) ! Only the real part matters
    end do
    do i = kp_complex + 1, np
        if (fkav(i) .gt. maxsr) maxsr = fkav(i)
    end do

    !      do i=1,np
    !      if(fkav(i).gt.parcut*maxsr) then
    !      x(i)=forza(i)/fkav(i)
    !      else
    !      x(i)=0.d0
    !      endif
    !      enddo

    countzero = 0
    do i = 1, kp_complex, 2
        if (fkav(i) .gt. parcut*maxsr .or. (yes_umrigar .and. fkav(i) .gt. 0.d0)) then
            do j = 1, nbin
                mat(j, i) = mat(j, i)/fkav(i)
                mat(j, i + 1) = mat(j, i + 1)/fkav(i)
            end do
            x(i) = forza(i)/fkav(i)
            if (fkav(i + 1) .ne. 0.d0) then
                x(i + 1) = forza(i + 1)/fkav(i)
            else
                x(i + 1) = 0.d0
            end if
        else
            if (rank_t .eq. 0 .and. fkav(i) .ne. 0.d0) write (6, *) ' Warning sr parameter too small !!! '
            countzero = countzero + 1
            do j = 1, nbin
                mat(j, i) = 0.d0
                mat(j, i + 1) = 0.d0
            end do
            fkav(i) = 0.d0
            fkav(i + 1) = 0.d0
            x(i) = 0.d0
            x(i + 1) = 0.d0
        end if
    end do
    do i = kp_complex + 1, np
        if (fkav(i) .gt. parcut*maxsr .or. (yes_umrigar .and. fkav(i) .ne. 0.d0)) then
            do j = 1, nbin
                mat(j, i) = mat(j, i)/fkav(i)
            end do
            x(i) = forza(i)/fkav(i)
        else
            if (rank_t .eq. 0 .and. fkav(i) .ne. 0.d0) &
                    &  write (6, *) ' Warning sr parameter too small !!! '
            countzero = countzero + 1
            do j = 1, nbin
                mat(j, i) = 0.d0
            end do
            fkav(i) = 0.d0
            x(i) = 0.d0
        end if
    end do

#if defined UNREL_DIAG && defined PARALLEL
    ! Different processor may have different countzero due to different
    ! accuracy in sqrt
    !  Just to be sure because also fkav was bcasted for safety.
    call mpi_bcast(countzero, 1, MPI_INTEGER, 0, comm_mpi, ierr)
#endif

    if (np - countzero .gt. 1) then

        do while (iter .lt. maxit .and. error .gt. eps)

            if (iter .eq. 0) h = x

            !      now psip= H h

#ifdef PARALLEL
#ifdef UNREL_DIAG
            call bcast_real(h, np, 0, comm_mpi)
#endif
#endif
            if (kpc .ne. 0) then
                if (symmagp) then
                    if (npk .gt. 0) then
                        call dgemv('N', nbin, npk, 1.d0, mat(1, kpr), nbin, h(kpr), 1, 0.d0, psip(npp), 1)
                    else
                        psip(npp:npp + nbin - 1) = 0.d0
                    end if

                    psip(1:kp_complex) = 0.d0
                    !      real-complex real
                    call dgemv('N', nbin, kpc, 1.d0, mat, 2*nbin, h, 2, 1.d0, psip(npp), 1)
                    call dgemv('T', nbin, kpc, 1.d0, mat, 2*nbin, psip(npp), 1, 1.d0, psip, 2)
                    if (npk .gt. 0) &
                        call dgemv('T', nbin, npk, 1.d0, mat(1, kpr), nbin, psip(npp), 1, 0.d0, psip(kpr), 1)
                    !      Adding remaining part coming form imag matrix.
                    call dgemv('N', nbin, kpc, 1.d0, mat(1, 2), 2*nbin, h, 2, 0.d0, psip(npp), 1)
                    call dgemv('T', nbin, kpc, 1.d0, mat(1, 2), 2*nbin, psip(npp), 1, 1.d0, psip, 2)
                else
                    !   Making the ARROCCO v_1,v_2 as above
                    do i = 1, kp_complex, 2
                        psip(i) = h(i + 1)
                        psip(i + 1) = -h(i)
                    end do
                    call dgemv('N', nbin, kp_complex, 1.d0, mat, nbin, psip, 1, 0.d0, psip(npp), 1)
                    !   Now psip contains mat x  |hi|
                    !                            |-hr|
                    !                            |0|
                    call dgemv('T', nbin, kp_complex, 1.d0, mat, nbin, psip(npp), 1, 0.d0&
                            &, psip, 1)
                    !     ARROCCO^T  changing -v_2 with v_1 as above comments
                    do i = 1, kp_complex, 2
                        cost = psip(i)
                        psip(i) = -psip(i + 1)
                        psip(i + 1) = cost
                    end do

                    if (npk .gt. 0) psip(kpr:np) = 0.d0 ! Reinizializing psip to zero real parameters.
                    !   Now adding the total stored matrix  to psip
                    call dgemv('N', nbin, np, 1.d0, mat, nbin, h, 1, 0.d0, psip(npp), 1)
                    !   Collecting mat^T x mat and add to psip
                    call dgemv('T', nbin, np, 1.d0, mat, nbin, psip(npp), 1, 1.d0, psip, 1)
                end if
            else
                call dgemv('N', nbin, npk, 1.d0, mat(1, kpr), nbin, h(kpr), 1, 0.d0, psip(npp), 1)
                call dgemv('T', nbin, npk, 1.d0, mat(1, kpr), nbin, psip(npp), 1, 0.d0, psip(kpr), 1)
            end if
            !      Projecting
            do i = 1, np
                if (fkav(i) .eq. 0.d0) psip(i) = 0.d0
            end do

#ifdef  PARALLEL
            if (kpc .gt. 0 .and. symmagp) then
!       To save a little communication.
                if (npk .gt. 0) call dcopy(npk, psip(kpp), 1, psip(kpcp), 1)
                call dcopy(kpc, psip, 2, psip(npp), 1)
                call reduce_base_real(npt, psip(npp), comm_mpi, -1)
                psip(1:kp_complex) = 0.d0
                call dcopy(kpc, psip(npp), 1, psip, 2)
                if (npk .gt. 0) call dcopy(npk, psip(kpcp), 1, psip(kpp), 1)
            else
                call reduce_base_real(np, psip, comm_mpi, -1)
            end if
#endif
            !       regularization
            if (yes_umrigar) then
                do ii = 1, np
                    if (fkav(ii) .ne. 0.d0) psip(ii) = psip(ii) + cost_umrigar*h(ii)/fkav(ii)**2
                end do
            end if
            call daxpy(np, epsdgel, h, 1, psip, 1)

            if (iter .ne. 0) then

                cost = ddot(np, h, 1, psip, 1)
                !       gauge invariant error
#if defined  UNREL_DIAG && defined PARALLEL
                call mpi_bcast(cost, 1, MPI_DOUBLE_PRECISION, 0, comm_mpi, ierr)
#endif
                if (cost .gt. 0.d0) then
                    lambda = -ddot(np, h, 1, g, 1)/cost
                    error = abs(lambda)*dsqrt(cost)
                else
                    if (rank .eq. 0)                                                   &
                            & write (6, *) ' ERROR  Non positive  definite matrix !!!'

#ifdef PARALLEL
                    call mpi_finalize(ierr)
#endif
                    stop

                end if

                ng2 = dnrm2(np, g, 1)**2
                !       write(6,*) ' scalar g A h, h A h ng2  '
                !    1,ddot(np,g,1,psip,1),ddot(np,h,1,psip,1),ng2,ddot(np,h,1,g,1)
                !       write(6,*) ' Overlap g g_i-1 ',ng2+lambda*ddot(np,psip,1,g,1)
                !       change g

                call daxpy(np, lambda, psip, 1, g, 1)

                gamma = lambda*ddot(np, psip, 1, g, 1)/ng2

                !       change  x  the solution
                lambda = -lambda
                call daxpy(np, lambda, h, 1, x, 1)

                !       change h
                do i = 1, np
                    h(i) = gamma*h(i) + g(i)
                end do

                !       compute function for check
                !      call dgemv('T',np,nbin,1.d0,mat,np,x,1,0.d0,psip(npp),1)
                !      call dgemv('N',np,nbin,1.d0,mat,np,psip(npp),1,0.d0,psip,1)
                !       call daxpy(np,epsdgel,x,1,psip,1)
                !        Fun=0.5d0*ddot(np,psip,1,x,1)
                !        write(6,*) ' New gradient at new position '
                !        do i=1,np
                !        if(fkav(i).gt.0.d0) then
                !        Fun=Fun-x(i)*forza(i)/fkav(i)
                !        write(6,*) i,g(i),-psip(i)+forza(i)/fkav(i)
                !        endif
                !        enddo

                !       error=abs(lambda*dnrm2(np,h,1))
                !       errorg=dnrm2(np,g,1)
                !       Value of the function
            else
                do i = 1, np
                    if (fkav(i) .gt. 0.d0) then
                        h(i) = -psip(i) + forza(i)/fkav(i)
                    else
                        h(i) = 0.d0
                    end if
                end do
                call dcopy(np, h, 1, g, 1)
                error = dnrm2(np, g, 1)
            end if
#ifdef UNREL_DIAG
!   The master force the calculation to be the same for all.
            call mpi_bcast(error, 1, MPI_DOUBLE_PRECISION, 0, comm_mpi, ierr)
#endif
            !       error=max(error,errorg)

            !      write(iopen) h,g

            iter = iter + 1

        end do
        if (rank .eq. 0) write (6, *) ' Output  cg  iter error   =', iter, error
        ! endif countzero
    end if
    do i = 1, np
        if (fkav(i) .ne. 0.d0) then
            x(i) = x(i)/fkav(i)
        else
            x(i) = 0.d0
        end if
    end do
    !      write(6,*) ' Check orthogonality gs '

    !     do j=1,iter
    !      rewind(iopen)
    !       do k=1,j-1
    !       read(iopen)
    !       enddo
    !       read(iopen) h,g
    !       call dcopy(np,g,1,psip,1)

    !      rewind(iopen)
    !      do i=1,iter
    !      read(iopen) h,g
    !       write(6,*) i,j,ddot(np,psip,1,g,1)
    !      enddo

    !     enddo
#ifdef UNREL_DIAG
!   The master force the calculation to be the same for all.
    call bcast_real(x, np, 0, comm_mpi)
#endif

    return
end
