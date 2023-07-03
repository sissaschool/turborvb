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

subroutine conjginvs(np, npp, nbin, rank_out, rank, comm_mpi, mat, forza&
        &, g, h, psip, maxit, eps, epsdgel, fkav, x, parcut, nproc, eps_umrigar)
    implicit none
    integer np, npp, i, ii, j, k, nl, nbin, maxit, iter, rank, info, np3, np4, rank_out&
            &, npm, indj, itert, ierr, jblock, countzero, nproc, nbinr, nppp, comm_mpi
    real*8 mat(npp, nbin), h(np), g(np), psip(*), eps                      &
            &, epsdgel, sv, eigmin, z, error, dnrm2, fkav(*), forza(*), cost, cost_umrigar &
            &, lambda, ng2, ddot, gamma, Fun, errorg, x(np), parcut, maxsr, eps_umrigar
    logical yes_umrigar
#ifdef PARALLEL
    include 'mpif.h'
    call mpi_barrier(comm_mpi, ierr)
#endif
    !      input mat(i,j)= reduced SR matrix
    !      input forza(i) i the force acting on the parameter i
    !      fkav(i)=s_{i,i} where:
    !      s_{i,j} = reduce SR matrix
    !      Output  x  solution of the linear system  s x =forza
    !      within a tollerance eps
    !      For QMC convenience (avoid too much noise in the inversion)
    !      The matrix s is regularized by s_{i,i} ---> s_{i,i} (1 + epsdgel)
    !
    !      notice that in the last block for jblock=(nproc-1)*nbin+1, nproc *n
    !      i.e. may be jblock>npp. But it is never used this info.

    !      recomputing fkav inside for stability reasons. if zero do not
    !      change
    nbinr = np - rank*nbin
    if (nbinr .gt. nbin) nbinr = nbin
    if (nbinr .gt. 0) then
        if (rank .ne. 0) fkav(1:rank*nbin) = 0.d0
        if (rank*nbin + nbinr .lt. np) fkav(rank*nbin + nbinr + 1:np) = 0.d0
        do j = 1, nbinr
            indj = rank*nbin + j
            if (fkav(indj) .ne. 0.d0) fkav(indj) = mat(indj, j)
        end do
    else
        fkav(1:np) = 0.d0
    end if
#ifdef PARALLEL
    call reduce_base_real(np, fkav, comm_mpi, -1)
#endif
    if (eps_umrigar .ne. 0.d0) then
        yes_umrigar = .true.
        cost_umrigar = eps_umrigar
        do i = 1, np
            if (fkav(i) .gt. 0.d0) then
                fkav(i) = fkav(i) + eps_umrigar
            else
                fkav(i) = 0.d0
            end if
        end do
    else
        yes_umrigar = .false.
        cost_umrigar = 0.d0
    end if

    iter = 0
    itert = 0
    error = 2.d0*eps

    !      write(6,*) ' parcut inside =',rank,parcut,parcut2

    nppp = nbin*nproc + 1
    !      write(6,*) ' Initial fkav , v1 '
    countzero = 0
    maxsr = 0.d0
    do i = 1, np
        if (fkav(i) .gt. maxsr) maxsr = fkav(i)
    end do

    do i = 1, np
        if (fkav(i) .gt. parcut*maxsr .or. (yes_umrigar .and. fkav(i) .gt. 0.d0)) then
            fkav(i) = dsqrt(fkav(i))
            x(i) = forza(i)/fkav(i)
        else
            countzero = countzero + 1
            x(i) = 0.d0
            fkav(i) = 0.d0
        end if
        !      write(6,*) i,x(i),fkav(i)
    end do

#ifdef UNREL_DIAG
! Different processor may have different fkav due to different
! accuracy in sqrt
    call bcast_real(fkav, np, 0, comm_mpi)
#endif
    !     preconditioning the matrix
    !
    do i = 1, np
        if (fkav(i) .gt. 0.d0) then
            do j = 1, nbinr
                if (fkav(rank*nbin + j) .gt. 0.d0) then
                    mat(i, j) = mat(i, j)/(fkav(i)*fkav(rank*nbin + j))
                else
                    mat(i, j) = 0.d0
                end if
            end do
        else
            do j = 1, nbin
                mat(i, j) = 0.d0
            end do
        end if
    end do

#ifdef UNREL_DIAG
! Different processor may have different countzero due to different
! accuracy in sqrt
    call mpi_bcast(countzero, 1, MPI_INTEGER, 0, comm_mpi, ierr)
#endif

    if (np - countzero .gt. 1) then

        do while (iter .lt. maxit .and. error .gt. eps)

            if (itert .eq. 0) h = x

            !      now psip= H h

#ifdef PARALLEL
#ifdef UNREL_DIAG
!   The master force the calculation to be the same for all.
            call bcast_real(h, np, 0, comm_mpi)
#endif

            jblock = rank*nbin + 1
            if (nbinr .gt. 0) then
                call dgemv('N', np, nbinr, 1.d0, mat, npp, h(jblock), 1              &
           &, 0.d0, psip(nppp), 1)
            else
                call dscalzero(np, 0.d0, psip(nppp), 1)
            end if
!       call dgemv('T',np,nbin,1.d0,mat,npp,h,1,0.d0,psip(nppp),1)
!       call mpi_allgather(psip(nppp),nbin,MPI_DOUBLE_PRECISION
!    1,psip(jblock),nbin,MPI_DOUBLE_PRECISION,comm_mpi,ierr)
            call reduce_base_real_to(np, psip(nppp), psip, comm_mpi, -1)
#else
            call dgemv('T', np, nbin, 1.d0, mat, npp, h, 1, 0.d0, psip, 1)
#endif

            !       write(6,*) ' rank =',rank,psip(1)
            !       regularization
            if (yes_umrigar) then
                do ii = 1, np
                    psip(ii) = psip(ii) + cost_umrigar*h(ii)/fkav(ii)**2
                end do
            end if
            call daxpy(np, epsdgel, h, 1, psip, 1)

            if (itert .ne. 0) then

                cost = ddot(np, h, 1, psip, 1)
                !       gauge invariant error
#ifdef  UNREL_DIAG
!  Only the master has the correct value of cost
                call mpi_bcast(cost, 1, MPI_DOUBLE_PRECISION, 0, comm_mpi, ierr)
#endif
                if (cost .gt. 0.d0) then
                    lambda = -ddot(np, h, 1, g, 1)/cost
                    error = abs(lambda)*dsqrt(cost)
                else
                    if (rank .eq. 0)                                                   &
                            &  write (6, *) ' ERROR  Non positive  definite matrix !!!'
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

            !       error=max(error,errorg)

            !      write(iopen) h,g

            iter = iter + 1
            itert = itert + 1

            !      if(mod(itert,10).eq.0) itert=0
            !      if(rank.eq.0) write(6,*) ' Output  cg  iter error   =',iter,error

#ifdef  UNREL_DIAG
            call mpi_bcast(error, 1, MPI_DOUBLE_PRECISION, 0, comm_mpi, ierr)
#endif

        end do

        ! endif countzero
    end if
    if (rank_out .eq. 0) write (6, *) ' Output  cg  iter error   =', iter, error
    !      write(6,*) ' Check orthogonality gs '

    do i = 1, np
        if (fkav(i) .ne. 0.d0) then
            x(i) = x(i)/fkav(i)
        else
            x(i) = 0.d0
        end if
    end do

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
