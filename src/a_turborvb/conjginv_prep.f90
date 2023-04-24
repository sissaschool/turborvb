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

subroutine conjginv_prep(np, prep, nprepm, kp_complex, symmagp, nbin, rank_t, rank&
        &, comm_mpi, comm_raw, comm_col, mat, forza, g, h, psip, maxit, eps, epsdgel, fkav&
        &, x, parcut, eps_umrigar, yes_ontarget)
    implicit none
    integer np, npk, i, ii, j, k, nl, npp, nbin, maxit, iter, rank, info, np3, np4         &
            &, npm, ierr, countzero, comm_mpi, kpr, kpc, kp_complex, kpcp, npt, kpp, mini, maxi&
            &, npu, kpu, prep, nprepm, ind, comm_raw, comm_col, dimp, nptot, nproc, myid, rank_t
    real*8 mat(nprepm, nbin), h(np), g(np), psip(*), eps                       &
            &, epsdgel, sv, eigmin, z, error, dnrm2, fkav(*), forza(*), cost, eps_umrigar&
            &, cost_umrigar&
            &, lambda, ng2, ddot, gamma, Fun, errorg, x(np), parcut, maxsr, parcut_max, cost_mpi(2)
    logical symmagp, yes_umrigar, yes_ontarget
#ifdef PARALLEL
    include 'mpif.h'
    call mpi_barrier(comm_mpi, ierr)
#endif
    !  This code is roundoff free, as it has a master-slave paradigm and
    !  only the master compute the quantities necessary to the iteration
    !  with its own roundoff, not assumed equal to the slaves.
    !  For this reason at each iteration the h and g vectors are set consistent
    !  between different slaves as multiplication by a constant may have
    !  a slightly different value in different slaves.
    !
    !      mat is distributed in  raw groups with prep processor each.
    !      ind= mod(rank,prep)*nprepm+1
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
!   For GPU acceleration  everything  is assumed to be defined in the CPU
    npk = np - kp_complex ! the number of  real parameter
    kpr = kp_complex + 1 ! the first real parameter
    kpc = kp_complex/2
    mini = mod(rank, prep)*nprepm + 1
    maxi = min(mini + nprepm - 1, np)

    dimp = maxi - mini + 1
    nptot = nprepm*prep
    npp = nptot + 1
    do i = 1, dimp
        ind = mod(rank, prep)*nprepm + i
        if (ind .ge. kpr .and. fkav(ind) .ne. 0.d0) fkav(ind) = sum(mat(i, 1:nbin)**2)
    end do
    if (kpc .gt. 0) then
        do i = 1, dimp, 2
            ind = mod(rank, prep)*nprepm + i
            if (ind .lt. kpr) then
                cost = sum(mat(i, 1:nbin)**2) + sum(mat(i + 1, 1:nbin)**2)
                if (fkav(ind) .ne. 0.d0) fkav(ind) = cost
                if (fkav(ind + 1) .ne. 0.d0) fkav(ind + 1) = cost
            end if
        end do
    end if
#ifdef PARALLEL
    psip(1:nptot) = 0.d0
    if (maxi .ge. mini) then
        psip(mini:maxi) = fkav(mini:maxi)
    end if
    call reduce_base_real(nprepm, psip(mini), comm_col, -1)
    psip(npp:npp + nprepm - 1) = 0.d0
    do i = mini, maxi
        if (psip(i) .gt. 0.d0) then
            psip(i - mini + npp) = dsqrt(psip(i))
        else
            psip(i - mini + npp) = 0.d0
        end if
    end do

    call mpi_allgather(psip(npp), nprepm, mpi_double_precision&
   &, psip, nprepm, mpi_double_precision, comm_raw, ierr)
    fkav(1:np) = psip(1:np)

#else
    do i = 1, np
        if (fkav(i) .gt. 0.d0) then
            fkav(i) = dsqrt(fkav(i))
        else
            fkav(i) = 0.d0
        end if
    end do
#endif

    if (eps_umrigar .ne. 0.d0) then
        yes_umrigar = .true.
        cost_umrigar = eps_umrigar
        do i = 1, np
            if (fkav(i) .gt. 0.d0) then
                fkav(i) = sqrt(fkav(i)**2 + eps_umrigar)
            else
                fkav(i) = 0.d0
            end if
        end do
    else
        yes_umrigar = .false.
        cost_umrigar = 0.d0
    end if

    !       if(rank.eq.0) then
    !       write(6,*) ' fkav inside wrong ',sum(fkav(1:np)),sum(forza(1:np))
    !       do i=1,np
    !       write(6,*) i,fkav(i),forza(i)
    !       enddo
    !       endif

    iter = 0
    error = 2.d0*eps
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
    parcut_max = parcut*maxsr

#ifdef PARALLEL
    call bcast_real(parcut_max, 1, 0, comm_mpi)
#endif

    !      if(rank.eq.0) write(6,*) ' parcut_max = ',parcut_max

    !      if(rank.eq.0) then
    !      write(6,*) ' fkav inside '
    !      do i=1,np
    !      write(6,*) i,fkav(i)**2
    !      enddo
    !      endif
    !      call mpi_finalize(ierr)
    !      stop

    countzero = 0
    do i = 1, kp_complex, 2
        ind = i - mini + 1
        if (fkav(i) .gt. parcut_max .or. (yes_umrigar .and. fkav(i) .gt. 0.d0)) then
            if (i .le. maxi .and. i .ge. mini) then
                do j = 1, nbin
                    mat(ind, j) = mat(ind, j)/fkav(i)
                end do
            end if
            if (i + 1 .le. maxi .and. i + 1 .ge. mini) then
                do j = 1, nbin
                    mat(ind + 1, j) = mat(ind + 1, j)/fkav(i)
                end do
            end if
            x(i) = forza(i)/fkav(i)
            if (fkav(i + 1) .ne. 0.d0) then
                x(i + 1) = forza(i + 1)/fkav(i)
            else
                x(i + 1) = 0.d0
            end if
        else
            if (rank_t .eq. 0 .and. fkav(i) .ne. 0.d0) &
                    &  write (6, *) ' Warning sr parameter too small !!! '
            countzero = countzero + 1
            if (i .le. maxi .and. i .ge. mini) then
                do j = 1, nbin
                    mat(ind, j) = 0.d0
                end do
            end if
            fkav(i) = 0.d0
            x(i) = 0.d0
            if (i + 1 .ge. mini .and. i + 1 .le. maxi) then
                do j = 1, nbin
                    mat(ind + 1, j) = 0.d0
                end do
            end if
            fkav(i + 1) = 0.d0
            x(i + 1) = 0.d0
        end if
    end do
    do i = kp_complex + 1, np
        ind = i - mini + 1
        if (fkav(i) .gt. parcut_max .or. (yes_umrigar .and. fkav(i) .gt. 0.d0)) then
            if (i .ge. mini .and. i .le. maxi) then
                do j = 1, nbin
                    mat(ind, j) = mat(ind, j)/fkav(i)
                end do
            end if
            x(i) = forza(i)/fkav(i)
        else
            if (rank_t .eq. 0 .and. fkav(i) .ne. 0.d0) &
                    &  write (6, *) ' Warning sr parameter too small !!! '
            countzero = countzero + 1
            if (i .ge. mini .and. i .le. maxi) then
                do j = 1, nbin
                    mat(ind, j) = 0.d0
                end do
            end if
            fkav(i) = 0.d0
            x(i) = 0.d0
        end if
    end do
#ifdef _OFFLOAD
    h = 0.d0 ! inizialize h to avoid NaN
!$omp target data map(to:mat,h) map(alloc:psip(1:2*nptot+nbin)) if(yes_ontarget)
#endif

    if (np - countzero .gt. 1) then

        do while (iter .lt. maxit .and. error .gt. eps)

            if (iter .eq. 0) then
                h = x

                !      now psip= H h

#ifdef PARALLEL
!       At the beginning set consistent h among different processors
                psip(npp:npp + nprepm - 1) = 0.d0
                psip(npp:npp + maxi - mini) = h(mini:maxi)
                call bcast_real(psip, nprepm, 0, comm_col)
                h(mini:maxi) = psip(npp:npp + maxi - mini)
#endif
            end if
            psip(1:nptot) = 0.d0
            if (kpc .ne. 0) then
                if (symmagp) then
                    if (npk .gt. 0) then
                        if (maxi .ge. kpr) then
                            kpu = max(kpr, mini) - mini + 1
                            npu = maxi - max(kpr, mini) + 1
                            if (npu .gt. 0) then
                                call dgemv__('T', npu, nbin, 1.d0, mat(kpu, 1), nprepm, h(max(kpr, mini))&
                                &, 1, 0.d0, psip(npp), 1, yes_ontarget)
                            else
                                psip(npp:npp + nbin - 1) = 0.d0
                            end if
                        else
                            psip(npp:npp + nbin - 1) = 0.d0
                        end if
                    else
                        psip(npp:npp + nbin - 1) = 0.d0
                    end if

                    !      real-complex real
                    if (mini .le. kp_complex) then
                        npu = min(kp_complex, maxi) - mini + 1
                        psip(npp + nbin:npp + nbin + npu - 1) = 0.d0
                        call dcopy(npu/2, h(mini), 2, psip(npp + nbin), 2)
                        call dgemv__('T', npu, nbin, 1.d0, mat, nprepm, psip(npp + nbin), 1, 1.d0, psip(npp), 1, yes_ontarget)
                    end if
#ifdef PARALLEL
                    call reduce_base_real(nbin, psip(npp), comm_raw, -1)
#endif
                    if (mini .le. kp_complex .and. npu .gt. 0) then
                        !      Do not change the imaginary part
                        call dcopy(npu/2, psip(mini + 1), 2, psip(npp + nbin), 1)
                        call dgemv__('N', npu, nbin, 1.d0, mat, nprepm, psip(npp), 1, 1.d0, psip(mini), 1, yes_ontarget)
                        call dcopy(npu/2, psip(npp + nbin), 1, psip(mini + 1), 2)

                    end if
                    if (npk .gt. 0 .and. maxi .ge. kpr) then
                        kpu = max(kpr, mini) - mini + 1
                        npu = maxi - max(kpr, mini) + 1
                        if (npu .gt. 0) then
                            call dgemv__('N', npu, nbin, 1.d0, mat(kpu, 1), nprepm, psip(npp), 1, 0.d0&
                                                     &, psip(max(kpr, mini)), 1, yes_ontarget)
                        end if
                    end if
                    !      Adding remaining part coming form imag matrix.
                    if (mini .le. kp_complex) then
                        kpu = 2
                        npu = min(kp_complex, maxi) - mini + 1
                        if (npu .gt. 0) then
                            psip(npp + nbin:npp + nbin + npu - 1) = 0.d0
                            call dcopy(npu/2, h(mini), 2, psip(npp + nbin + 1), 2)
                            call dgemv__('T', npu, nbin, 1.d0, mat, nprepm, psip(npp + nbin), 1, 0.d0, psip(npp), 1, yes_ontarget)
                        else
                            psip(npp:npp + nbin - 1) = 0.d0
                        end if
                    else
                        psip(npp:npp + nbin - 1) = 0.d0
                    end if
#ifdef PARALLEL
                    call reduce_base_real(nbin, psip(npp), comm_raw, -1)
#endif
                    if (mini .le. kp_complex .and. npu .gt. 0) then
                        call dcopy(npu/2, psip(mini + 1), 2, psip(npp + nbin), 1)
                        call dcopy(npu/2, psip(mini), 2, psip(mini + 1), 2)
                        call dgemv__('N', npu, nbin, 1.d0, mat, nprepm, psip(npp), 1, 1.d0, psip(mini), 1, yes_ontarget)
                        !  shift to the real part

                        call dcopy(npu/2, psip(mini + 1), 2, psip(mini), 2)

                        call dcopy(npu/2, psip(npp + nbin), 1, psip(mini + 1), 2)
                    end if

                else

                    !   Making the ARROCCO v_1,v_2 as above
                    do i = mini, min(maxi, kp_complex), 2
                        psip(i) = h(i + 1)
                        psip(i + 1) = -h(i)
                    end do
                    if (mini .le. kp_complex) then
                        kpu = mini
                        npu = min(kp_complex, maxi) - mini + 1
                        if (npu .gt. 0) then
                            call dgemv__('T', npu, nbin, 1.d0, mat, nprepm, psip(kpu), 1, 0.d0, psip(npp), 1, yes_ontarget)
                        else
                            psip(npp:npp + nbin - 1) = 0.d0
                        end if
                    else
                        psip(npp:npp + nbin - 1) = 0.d0
                    end if

#ifdef PARALLEL
                    call reduce_base_real(nbin, psip(npp), comm_raw, -1)
#endif
                    !   Now psip contains mat x  |hi|
                    !                            |-hr|
                    !                            |0|
                    if (mini .le. kp_complex .and. npu .gt. 0) then
                        call dgemv__('N', npu, nbin, 1.d0, mat, nprepm, psip(npp), 1, 0.d0&
                                &, psip(kpu), 1, yes_ontarget)
                    end if
                    !     ARROCCO^T  changing -v_2 with v_1 as above comments
                    do i = mini, min(maxi, kp_complex), 2
                        cost = psip(i)
                        psip(i) = -psip(i + 1)
                        psip(i + 1) = cost
                    end do

                    if (npk .gt. 0) psip(kpr:nptot) = 0.d0 ! Reinizializing psip to zero real parameters.
                    !   Now adding the total stored matrix  to psip
                    npu = dimp
                    if (npu .gt. 0) then
                        call dgemv__('T', npu, nbin, 1.d0, mat, nprepm, h(mini), 1, 0.d0, psip(npp), 1, yes_ontarget)
                    else
                        psip(npp:npp + nbin - 1) = 0.d0
                    end if
                    !   Collecting mat^T x mat and add to psip
#ifdef PARALLEL
                    call reduce_base_real(nbin, psip(npp), comm_raw, -1)
#endif
                    if (npu .gt. 0) then
                        call dgemv__('N', npu, nbin, 1.d0, mat, nprepm, psip(npp), 1, 1.d0, psip(mini), 1, yes_ontarget)
                    end if
                end if
            else
                if (maxi .ge. kpr) then
                    kpu = max(kpr, mini) - mini + 1
                    npu = maxi - max(kpr, mini) + 1
                    if (npu .gt. 0) then
                        call dgemv__('T', npu, nbin, 1.d0, mat(kpu, 1), nprepm, h(max(kpr, mini))&
                        &, 1, 0.d0, psip(npp), 1, yes_ontarget)
                    else
                        psip(npp:npp + nbin - 1) = 0.d0
                    end if
                else
                    psip(npp:npp + nbin - 1) = 0.d0
                end if
#ifdef PARALLEL
                call reduce_base_real(nbin, psip(npp), comm_raw, -1)
#endif
                if (maxi .ge. kpr .and. npu .gt. 0) then
                    call dgemv__('N', npu, nbin, 1.d0, mat(kpu, 1), nprepm, psip(npp), 1, 0.d0&
                    &, psip(max(kpr, mini)), 1, yes_ontarget)
                end if
            end if
            !      Projecting
            do i = mini, maxi
                if (fkav(i) .eq. 0.d0) psip(i) = 0.d0
            end do
#ifdef PARALLEL
!  Collecting all the small pieces
            call reduce_base_real(nprepm, psip(mini), comm_col, -1)
#endif
            !        write(6,*) ' First step =',rank,sum(psip(mini:maxi))
            !        psip(1:mini-1)=0.d0
            !        psip(maxi+1:np)=0.d0
            !        call reduce_base_real(np,psip,comm_raw,-1)
            !        if(rank.eq.0) write(6,*) ' Total = ',sum(psip(1:np))
            !        call mpi_finalize(ierr)
            !        stop
            !       regularization
            if (dimp .gt. 0) then
                if (yes_umrigar) then
                    do ii = mini, mini + dimp - 1
                        if (fkav(ii) .ne. 0.d0) psip(ii) = psip(ii) + cost_umrigar*h(ii)/fkav(ii)**2
                    end do
                end if
                call daxpy(dimp, epsdgel, h(mini), 1, psip(mini), 1)
            end if
            !       call daxpy(np,epsdgel,h,1,psip,1)

            if (iter .ne. 0) then

                if (dimp .gt. 0) then
                    cost = ddot(dimp, h(mini), 1, psip(mini), 1)
                else
                    cost = 0.d0
                end if
#ifdef PARALLEL
                cost_mpi(1) = cost
                if (dimp .gt. 0) then
                    cost_mpi(2) = -ddot(dimp, h(mini), 1, g(mini), 1)
                else
                    cost_mpi(2) = 0.d0
                end if
                call reduce_base_real(2, cost_mpi, comm_raw, -1)
                call bcast_real(cost_mpi, 2, 0, comm_col)
                cost = cost_mpi(1)
#endif

                if (cost .gt. 0.d0) then
#ifdef PARALLEL
                    lambda = cost_mpi(2)/cost
                    error = abs(lambda)*dsqrt(abs(cost))
#else
                    lambda = -ddot(np, h, 1, g, 1)/cost
                    error = abs(lambda)*dsqrt(abs(cost))
#endif
                else
                    if (rank .eq. 0)&
                            & write (6, *) ' ERROR  Non positive  definite matrix !!!', cost
#ifdef PARALLEL
                    call mpi_finalize(ierr)
#endif
                    stop

                end if

                if (dimp .gt. 0) then
                    ng2 = dnrm2(dimp, g(mini), 1)**2
                else
                    ng2 = 0.d0
                end if
                !        ng2=dnrm2(np,g,1)**2

                if (dimp .gt. 0) call daxpy(dimp, lambda, psip(mini), 1, g(mini), 1)
                !        call daxpy(np,lambda,psip,1,g,1)

                if (dimp .gt. 0) then
                    gamma = ddot(dimp, psip(mini), 1, g(mini), 1)
                else
                    gamma = 0.d0
                end if
#ifdef PARALLEL
                cost_mpi(1) = ng2
                cost_mpi(2) = gamma
                call reduce_base_real(2, cost_mpi, comm_raw, -1)
                ng2 = cost_mpi(1)
                gamma = cost_mpi(2)
                gamma = gamma*lambda/ng2
                cost_mpi(2) = gamma
                call bcast_real(cost_mpi, 2, 0, comm_col)
                gamma = cost_mpi(2)
                ng2 = cost_mpi(1)
#endif
                !        gamma=lambda*ddot(np,psip,1,g,1)/ng2

                !       change  x  the solution
                lambda = -lambda
                if (dimp .gt. 0) call daxpy(dimp, lambda, h(mini), 1, x(mini), 1)
                !       call daxpy(np,lambda,h,1,x,1)

                !       change h
                do i = mini, maxi
                    !       do i=1,np
                    h(i) = gamma*h(i) + g(i)
                end do

            else
                !        do i=1,np
                do i = mini, maxi
                    if (fkav(i) .gt. 0.d0) then
                        h(i) = -psip(i) + forza(i)/fkav(i)
                    else
                        h(i) = 0.d0
                    end if
                end do
                lambda = 0.d0
                cost = 0.d0
                gamma = 0.d0
                ng2 = 0.d0
                if (dimp .gt. 0) then
                    call dcopy(dimp, h(mini), 1, g(mini), 1)
                    error = sum(g(mini:maxi)**2)
                else
                    error = 0.d0
                end if
#ifdef PARALLEL
                call reduce_base_real(1, error, comm_raw, -1)
                error = dsqrt(error)
                call bcast_real(error, 1, 0, comm_col)
#else
                error = dsqrt(error)
#endif
                !        call dcopy(np,h,1,g,1)
                !        error=dsqrt(sum(g(1:np)**2))

            end if

#ifdef PARALLEL
!   set consistent g and h
            psip(npp:npp + 2*nprepm - 1) = 0.d0
            psip(npp:npp + maxi - mini) = h(mini:maxi)
            psip(npp + nprepm:npp + nprepm + maxi - mini) = g(mini:maxi)
            call bcast_real(psip, 2*nprepm, 0, comm_col)
            h(mini:maxi) = psip(npp:npp + maxi - mini)
            g(mini:maxi) = psip(npp + nprepm:npp + nprepm + maxi - mini)
#endif

            iter = iter + 1
        end do
        if (rank .eq. 0) write (6, *) ' Output  cg  iter error   =', iter, error
        ! endif countzero
    end if
#ifdef _OFFLOAD
!$omp end target data
#endif
    do i = mini, maxi
        if (fkav(i) .ne. 0.d0) then
            x(i) = x(i)/fkav(i)
        else
            x(i) = 0.d0
        end if
    end do
#ifdef PARALLEL
    psip(npp:npp + nprepm - 1) = 0.d0
    if (maxi .ge. mini) then
        psip(npp:npp + maxi - mini) = x(mini:maxi)
    end if
!      To set all equal and avoid a final bcast global
    call bcast_real(psip(npp), nprepm, 0, comm_col)
    call mpi_allgather(psip(npp), nprepm, mpi_double_precision, psip, nprepm, mpi_double_precision, comm_raw, ierr)
    x(1:np) = psip(1:np)
#endif
    return
end
