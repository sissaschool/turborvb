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

subroutine bootparameter(nbin, efenergy, e, ndimj, ndimp, eta       &
        &, derr, eall, esav, efall, fk, dimfk, okav, weightall_c         &
        &, min_par, max_par, rank, comm_c, nproc_c, comm_f, nproc_f, nproc)
    !     use allio, only: col_id
    use constants, only: ipc
    use allio, only: ndims, ndimsp, ndimjp, ndimiesup, symmagp, yes_hermite, yes_correct, kaverage, lrdmc_der
    !,iesuptrans,iesupr_2&
    !     &,contraction,nozero_c,nelorb_c,jbradet,nnozero_c
    implicit none

    ! *** Does BootStrap to evalute force and its error for each bin
    !  energy !  derrivative of the ener! log derivative wf! energy x log de

    integer mdone, ix, iy, k, j, jj, dimfk, nbin, kmain, kk, i, np, ndimp, rank, nproc_c&
            &, ierr, comm_c, comm_f, min_par, max_par, nproc_f, nproc, ndimj, nstart
    real*8 e(ndimp, 2, *), givenmeaskk, enerboot(1:ipc), summpi(2*ipc + 1), outmpi(2*ipc + 1)&
            &, esav_c(1:ipc), esav_f(1:ipc), derr(ndimp), eall(2, ndimp), efall_c(ipc), efall_f(ipc)&
            &, eta(ndimp), efenergy(1 + ipc, *), weightall_f, givenmeaski&
            &, cost, weight, weightall_c, fk(dimfk, *), okav(*)&
            &, nsamp_c, nsamp_f, ekk(2), eki(2), esav(ipc), efall(ipc), weightall, nsamp&
            &, ekt(2), eti(2)
    !      Any processor has nbin samples
#ifdef PARALLEL
    include 'mpif.h'
#endif
    !       initialization opara e opars
    np = max_par - min_par + 1
    if (ipc .eq. 1) then
        mdone = max_par
    else
        if (yes_correct) then
            mdone = min(max_par, ndimj)
            nstart = ndimjp
        else
            mdone = max_par
            nstart = max_par + 1
        end if
    end if

    nsamp_c = dble(nbin)*dble(nproc_c)
    nsamp_f = dble(nbin)*dble(nproc_f)
    nsamp = dble(nbin)*dble(nproc)

    if (nsamp_c .le. 1 .or. nsamp_f .le. 1 .or. max_par .gt. ndimp .or. min_par .lt. 1) then
        if (rank .eq. 0 .and. (nsamp_c .le. 1 .or. nsamp_f .le. 1)) then
            write (6, *)                                                     &
                    &' You should have more than 1 bin to estimate error bars !!!'
            write (6, *) ' nbin , nproc_c, nproc_f =', nbin, nproc_c, nproc_f
        end if
        if (max_par .gt. ndimp .or. min_par .lt. 1) then
            if (rank .eq. 0) write (6, *) 'ERROR outside range 1<= # <= ndimp'
        end if
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif
        stop
    end if

    !      calculation eall once for all
    efall_c = 0.d0
    esav_c = 0.d0
    weightall_c = 0.d0
    do i = 1, nbin
        weightall_c = weightall_c + efenergy(1 + ipc, i)
        efall_c(1:ipc) = efall_c(1:ipc) + efenergy(1:ipc, i)
        esav_c(1:ipc) = esav_c(1:ipc) + efenergy(1:ipc, i)**2/efenergy(1 + ipc, i)
    end do
    weightall = weightall_c
    efall = efall_c
    esav = esav_c
    weightall_f = weightall_c
    efall_f = efall_c
    esav_f = esav_c
#ifdef PARALLEL
    if (nproc_c .gt. 1) then
        summpi(1) = weightall_c
        summpi(2:2 + ipc - 1) = esav_c(1:ipc)
        summpi(2 + ipc:2*ipc + 1) = efall_c(1:ipc)
        call mpi_allreduce(summpi, outmpi, 2*ipc + 1, MPI_DOUBLE_PRECISION     &
    &  , MPI_SUM, comm_c, ierr)
        weightall_c = outmpi(1)
        esav_c(1:ipc) = outmpi(2:2 + ipc - 1)
        efall_c(1:ipc) = outmpi(2 + ipc:2*ipc + 1)
        weightall = weightall_c
        esav = esav_c
        efall = efall_c
    end if
    if (nproc_f .gt. 1 .and. nproc_c .ne. nproc_f) then
        summpi(1) = weightall_f
        summpi(2:2 + ipc - 1) = esav_f(1:ipc)
        summpi(2 + ipc:2*ipc + 1) = efall_f(1:ipc)
        call mpi_allreduce(summpi, outmpi, 2*ipc + 1, MPI_DOUBLE_PRECISION     &
    &  , MPI_SUM, comm_f, ierr)
        weightall_f = outmpi(1)
        esav_f = outmpi(2:2 + ipc - 1)
        efall_f = outmpi(2 + ipc:2*ipc + 1)
    elseif (nproc_c .eq. nproc_f) then
        weightall_f = weightall_c
        esav_f = esav_c
        efall_f = efall_c
    end if
    if (.not. kaverage) then
        if (nproc .ne. nproc_c .and. nproc .ne. nproc_f) then
            summpi(1) = weightall
            summpi(2:2 + ipc - 1) = esav(1:ipc)
            summpi(2 + ipc:2*ipc + 1) = efall(1:ipc)
            call mpi_allreduce(summpi, outmpi, 2*ipc + 1, MPI_DOUBLE_PRECISION     &
        &  , MPI_SUM, MPI_COMM_WORLD, ierr)
            weightall = outmpi(1)
            esav = outmpi(2:2 + ipc - 1)
            efall = outmpi(2 + ipc:2*ipc + 1)
        elseif (nproc .eq. nproc_c) then
            weightall = weightall_c
            esav = esav_c
            efall = efall_c
        elseif (nproc .eq. nproc_f) then
            weightall = weightall_f
            esav = esav_f
            efall = efall_f
        end if
    end if
#endif
    esav_c = esav_c/weightall_c
    efall_c = efall_c/weightall_c
    esav_c = dsqrt(dabs(esav_c - efall_c**2))/dsqrt(nsamp_c - 1.d0)

    esav_f = esav_f/weightall_f
    efall_f = efall_f/weightall_f
    esav_f = dsqrt(dabs(esav_f - efall_f**2))/dsqrt(nsamp_f - 1.d0)

    esav = esav/weightall
    efall = efall/weightall
    if (kaverage) then
        esav = dsqrt(dabs(esav - efall**2))/dsqrt(nsamp_c - 1.d0)
    else
        esav = dsqrt(dabs(esav - efall**2))/dsqrt(nsamp - 1.d0)
    end if
    eall(:, min_par:max_par) = 0.d0
    do kk = min_par, max_par
        do k = 1, nbin
            do j = 1, 2
                eall(j, kk) = eall(j, kk) + e(kk, j, k)
            end do
        end do
    end do ! enddo kk
    eall(:, min_par:max_par) = eall(:, min_par:max_par)/weightall_c
#ifdef PARALLEL
    if (nproc_c .gt. 1) then
        call reduce_base_real(2*np, eall(1, min_par), comm_c, -1)
    end if
#endif
    !       calculation  average <O_k>
    do kk = min_par, max_par
        okav(kk) = eall(1, kk)
    end do

    !      calculation fk
    if (nproc_c .eq. nproc_f) then
        eta(min_par:max_par) = 0.d0
        derr(min_par:max_par) = 0.d0
        do kk = min_par, mdone
            do kmain = 1, nbin
                weight = efenergy(1 + ipc, kmain)
                enerboot(1) = (efall_c(1)*weightall_c - efenergy(1, kmain))/(weightall_c - weight)
                do jj = 1, 2
                    ekk(jj) = (eall(jj, kk)*weightall_c - e(kk, jj, kmain))                 &
                            & /(weightall_c - weight) !wt(jj)
                end do
                !       now average correlation function on the given bootstrap
                if (lrdmc_der) then
                    givenmeaskk = -ekk(1)
                else
                    givenmeaskk = 2.d0*(-ekk(2) + enerboot(1)*ekk(1))
                end if
                eta(kk) = eta(kk) + givenmeaskk
                derr(kk) = derr(kk) + givenmeaskk**2
                fk(kmain, kk) = givenmeaskk
            end do
        end do
        if (ipc .eq. 2) then

            do kk = nstart, max_par, 2
                do kmain = 1, nbin
                    weight = efenergy(1 + ipc, kmain)
                    enerboot(1:2) = (efall_c(1:2)*weightall_c - efenergy(1:2, kmain))/(weightall_c - weight)
                    do jj = 1, 2
                        ekk(jj) = (eall(jj, kk)*weightall_c - e(kk, jj, kmain))                 &
                                & /(weightall_c - weight) !wt(jj)
                    end do
                    do jj = 1, 2
                        eki(jj) = (eall(jj, kk + 1)*weightall_c - e(kk + 1, jj, kmain))                 &
                                & /(weightall_c - weight) !wt(jj)
                    end do
                    !       now average correlation function on the given bootstrap
                    if (lrdmc_der) then
                        givenmeaskk = -ekk(1)
                    else
                        givenmeaskk = 2.d0*(-ekk(2) + enerboot(1)*ekk(1) - enerboot(2)*eki(1))
                    end if
                    eta(kk) = eta(kk) + givenmeaskk
                    derr(kk) = derr(kk) + givenmeaskk**2
                    fk(kmain, kk) = givenmeaskk
                    if (lrdmc_der) then
                        givenmeaski = -eki(1)
                    else
                        givenmeaski = 2.d0*(-eki(2) + enerboot(2)*ekk(1) + enerboot(1)*eki(1))
                    end if
                    eta(kk + 1) = eta(kk + 1) + givenmeaski
                    derr(kk + 1) = derr(kk + 1) + givenmeaski**2
                    fk(kmain, kk + 1) = givenmeaski
                end do
            end do
        end if

        eta(min_par:max_par) = eta(min_par:max_par)/nsamp_c
        derr(min_par:max_par) = derr(min_par:max_par)/nsamp_c
#ifdef PARALLEL
        if (nproc_c .gt. 1) then
            call reduce_base_real(np, eta(min_par), comm_c, -1)
            call reduce_base_real(np, derr(min_par), comm_c, -1)
        end if
#endif
    elseif (nproc_c .ge. 1) then
        eta(min_par:max_par) = 0.d0
        do kk = min_par, mdone
            do kmain = 1, nbin
                weight = efenergy(1 + ipc, kmain)
                enerboot(1) = (efall_c(1)*weightall_c - efenergy(1, kmain))/(weightall_c - weight)
                do jj = 1, 2
                    ekk(jj) = (eall(jj, kk)*weightall_c - e(kk, jj, kmain))                 &
                            & /(weightall_c - weight) !wt(jj)
                end do
                !       now average correlation function on the given bootstrap
                if (lrdmc_der) then
                    givenmeaskk = -ekk(1)
                else
                    givenmeaskk = 2.d0*(-ekk(2) + enerboot(1)*ekk(1))
                end if
                eta(kk) = eta(kk) + givenmeaskk
                fk(kmain, kk) = givenmeaskk
            end do
        end do
        if (ipc .eq. 2) then
            do kk = nstart, max_par, 2
                do kmain = 1, nbin
                    weight = efenergy(1 + ipc, kmain)
                    enerboot(1:2) = (efall_c(1:2)*weightall_c - efenergy(1:2, kmain))/(weightall_c - weight)
                    do jj = 1, 2
                        ekk(jj) = (eall(jj, kk)*weightall_c - e(kk, jj, kmain))&
                                & /(weightall_c - weight) !wt(jj)
                    end do
                    do jj = 1, 2
                        eki(jj) = (eall(jj, kk + 1)*weightall_c - e(kk + 1, jj, kmain))&
                                & /(weightall_c - weight) !wt(jj)
                    end do
                    !       now average correlation function on the given bootstrap
                    if (lrdmc_der) then
                        givenmeaskk = -ekk(1)
                    else
                        givenmeaskk = 2.d0*(-ekk(2) + enerboot(1)*ekk(1) - enerboot(2)*eki(1))
                    end if
                    eta(kk) = eta(kk) + givenmeaskk
                    fk(kmain, kk) = givenmeaskk
                    if (lrdmc_der) then
                        givenmeaski = -eki(1)
                    else
                        givenmeaski = 2.d0*(-eki(2) + enerboot(2)*ekk(1) + enerboot(1)*eki(1))
                    end if
                    eta(kk + 1) = eta(kk + 1) + givenmeaski
                    fk(kmain, kk + 1) = givenmeaski
                end do
            end do

        end if
        eta(min_par:max_par) = eta(min_par:max_par)/nsamp_c
#ifdef PARALLEL
        if (nproc_c .gt. 1) then
            call reduce_base_real(np, eta(min_par), comm_c, -1)
        end if
#endif
    end if
    !      Now calculation matrix cov using the prefactor of the
    !      Jackknife
    cost = dsqrt((nsamp_c - 1.d0)/nsamp_c)
    do j = min_par, max_par
        do i = 1, nbin
            fk(i, j) = cost*(fk(i, j) - eta(j))
        end do
    end do

    if (nproc_c .ne. nproc_f) then
        eall(:, min_par:max_par) = 0.d0
        do kk = min_par, max_par
            do k = 1, nbin
                do j = 1, 2
                    eall(j, kk) = eall(j, kk) + e(kk, j, k)
                end do
            end do
        end do ! enddo kk
        eall(:, min_par:max_par) = eall(:, min_par:max_par)/weightall_f
#ifdef PARALLEL
        if (nproc_f .gt. 1) then
            call reduce_base_real(2*np, eall(1, min_par), comm_f, -1)
        end if
#endif

        eta(min_par:max_par) = 0.d0
        derr(min_par:max_par) = 0.d0
        do kk = min_par, mdone
            do kmain = 1, nbin
                weight = efenergy(1 + ipc, kmain)
                enerboot(1) = (efall_f(1)*weightall_f - efenergy(1, kmain))/(weightall_f - weight)
                do jj = 1, 2
                    ekk(jj) = (eall(jj, kk)*weightall_f - e(kk, jj, kmain))                 &
                            & /(weightall_f - weight) !wt(jj)
                end do
                !       now average correlation function on the given bootstrap
                if (lrdmc_der) then
                    givenmeaskk = -ekk(1)
                else
                    givenmeaskk = 2.d0*(-ekk(2) + enerboot(1)*ekk(1))
                end if
                eta(kk) = eta(kk) + givenmeaskk
                derr(kk) = derr(kk) + givenmeaskk**2
            end do
        end do
        if (ipc .eq. 2) then
            do kk = nstart, max_par, 2
                do kmain = 1, nbin
                    weight = efenergy(1 + ipc, kmain)
                    enerboot(1:ipc) = (efall_f(1:ipc)*weightall_f - efenergy(1:ipc, kmain))/(weightall_f - weight)
                    do jj = 1, 2
                        ekk(jj) = (eall(jj, kk)*weightall_f - e(kk, jj, kmain))                 &
                                & /(weightall_f - weight) !wt(jj)
                    end do
                    do jj = 1, 2
                        eki(jj) = (eall(jj, kk + 1)*weightall_f - e(kk + 1, jj, kmain))                 &
                                & /(weightall_f - weight) !wt(jj)
                    end do
                    !       now average correlation function on the given bootstrap
                    if (lrdmc_der) then
                        givenmeaskk = -ekk(1)
                    else
                        givenmeaskk = 2.d0*(-ekk(2) + enerboot(1)*ekk(1) - enerboot(2)*eki(1))
                    end if
                    eta(kk) = eta(kk) + givenmeaskk
                    derr(kk) = derr(kk) + givenmeaskk**2
                    if (lrdmc_der) then
                        givenmeaski = -eki(1)
                    else
                        givenmeaski = 2.d0*(-eki(2) + enerboot(2)*ekk(1) + enerboot(1)*eki(1))
                    end if
                    eta(kk + 1) = eta(kk + 1) + givenmeaskk
                    derr(kk + 1) = derr(kk + 1) + givenmeaskk**2
                end do
            end do
        end if

        eta(min_par:max_par) = eta(min_par:max_par)/nsamp_f
        derr(min_par:max_par) = derr(min_par:max_par)/nsamp_f
#ifdef PARALLEL
        if (nproc_f .gt. 1) then
            call reduce_base_real(np, eta(min_par), comm_f, -1)
            call reduce_base_real(np, derr(min_par), comm_f, -1)
        end if
#endif
    end if

    cost = nsamp_f - 1.d0
    do kk = min_par, max_par
        derr(kk) = dsqrt(cost*(derr(kk) - eta(kk)**2))
    end do

    if (lrdmc_der) then
        do kk = min_par, mdone
            eta(kk) = -eall(1, kk)
        end do
    else
        do kk = min_par, mdone
            eta(kk) = 2.d0*(-eall(2, kk) + efall_f(1)*eall(1, kk))
        end do
    end if
    if (ipc .eq. 2) then
        if (lrdmc_der) then
            do kk = nstart, max_par, 2
                eta(kk) = -eall(1, kk)
                eta(kk + 1) = -eall(1, kk + 1)
            end do
        else
            do kk = nstart, max_par, 2
                eta(kk) = 2.d0*(-eall(2, kk) + efall_f(1)*eall(1, kk) - efall_f(2)*eall(1, kk + 1))
                eta(kk + 1) = 2.d0*(-eall(2, kk + 1) + efall_f(2)*eall(1, kk) + efall_f(1)*eall(1, kk + 1))
            end do
        end if
    end if
    !!    Vanishing imaginary part in this case.
    !              if(yes_hermite.or.contraction.eq.0) then
    !               do kk=ndimsp,ndimiesup,2
    !               eta(kk+1)=0.d0
    !!               enddo
    !              else
    !!         Vanish Imaginary part of Z forces only
    !               do i=1,iesupr_2
    !                if(iesuptrans(i).ne.0) then
    !!                eta(2*iesuptrans(i)+ndims)=0.d0
    !                endif
    !               enddo
    !              endif
    !          else
    !!         Vanish Imaginary part of Z forces only
    !          if(contraction.ne.0) then
    !!          do i=1,iesupr_2
    !             if(iesuptrans(i).ne.0) then
    !             eta(2*iesuptrans(i)+ndims)=0.d0
    !             endif
    !          enddo
    !          else
    !               do kk=ndimsp,ndimiesup,2
    !               eta(kk+1)=0.d0
    !               enddo
    !!          endif
    !         endif ! symmagp
    !         if(yes_hermite.and.symmagp) then
    !!  Vanishing real part of diagonal terms forces
    !          do i=1,nnozero_c
    !          j=abs(jbradet(i))
    !          iy=(nozero_c(i)-1)/nelorb_c+1
    !          ix=nozero_c(i)-(iy-1)*nelorb_c
    !           if(ix.eq.iy.and.j.ne.0) then
    !           eta(4*j-2+ndimj)=0.d0
    !           endif
    !          enddo
    !         endif
    !       endif ! ipc=2
    !!      write(6,*) ' Output fk=',rank,sum(fk(:,min_par:max_par)**2)
    !!      write(6,*) ' Output eta=',rank,sum(eta(min_par:max_par)**2)
    !!      call mpi_finalize(ierr)
    !!      stop

    return
end
