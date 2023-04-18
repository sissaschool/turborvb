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

subroutine bootforcecov(nbin, efenergy, ndimp, ef, ieskin, eta       &
        &, derr, eall, fk, dimfk, okav, rank, comm_f, nproc_f, comm_c, nproc_c)
    !     use allio, only: col_id
    use constants, only: ipc
    use allio, only: normcorr
    implicit none

    ! *** Does BootStrap to evalute force and its error for each bin
    !  energy !  derrivative of the ener! log derivative wf! energy x log de

    integer k, j, jj, dimfk, ieskin, nbin, kmain, kk, i&
            &, comm_c, np, ndimp, rank&
            &, comm_f, ierr, nproc_f, nproc_c, nmax, nmin
    real*8 ef(ieskin, 3, *), ekk(4), givenmeaskk     &
            &, enerboot, summpi(3), outmpi(3)             &
            &, derr(ieskin + ndimp), eall(3, ieskin + ndimp)              &
            &, eta(ieskin + ndimp), efenergy(1 + ipc, *), cost, weight &
            &, fk(dimfk, *), okav(*), nsamp_c, nsamp_f&
            &, weightall_f, esav_f, efall_f, efall_c, esav_c, weightall_c
#ifdef PARALLEL
    include 'mpif.h'
#endif

    !       initialization opara e opars

    np = ieskin
    nmax = ndimp + ieskin
    nmin = ndimp + 1
    nsamp_f = dble(nbin)*dble(nproc_f)
    nsamp_c = dble(nbin)*dble(nproc_c)
    if (nsamp_c .le. 1 .or. nsamp_f .le. 1) then
        if (rank .eq. 0) then
            write (6, *)                                                     &
                    &' You should have more than 1 bin to estimate error bars !!!'
            write (6, *) ' nbin , nproc =', nbin, nproc_c, nproc_f
        end if
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif
        stop
    end if

    !      calculation eall once for all
    weightall_c = 0.d0
    efall_c = 0.d0
    esav_c = 0.d0
    do i = 1, nbin
        weightall_c = weightall_c + efenergy(1 + ipc, i)
        efall_c = efall_c + efenergy(1, i)
        esav_c = esav_c + efenergy(1, i)**2/efenergy(1 + ipc, i)
    end do
    weightall_f = weightall_c
    esav_f = esav_c
    efall_f = efall_c

#ifdef PARALLEL
    if (nproc_f .gt. 1) then
        summpi(1) = weightall_f
        summpi(2) = esav_f
        summpi(3) = efall_f
        call mpi_allreduce(summpi, outmpi, 3, MPI_DOUBLE_PRECISION     &
    &  , MPI_SUM, comm_f, ierr)
        weightall_f = outmpi(1)
        esav_f = outmpi(2)
        efall_f = outmpi(3)
!        Average for the variational parameters always assumed for all
!        processors
    end if
    if (nproc_c .gt. 1 .and. nproc_c .ne. nproc_f) then
        summpi(1) = weightall_c
        summpi(2) = esav_c
        summpi(3) = efall_c
        call mpi_allreduce(summpi, outmpi, 3, MPI_DOUBLE_PRECISION     &
    &  , MPI_SUM, comm_c, ierr)
        weightall_c = outmpi(1)
        esav_c = outmpi(2)
        efall_c = outmpi(3)
    elseif (nproc_c .eq. nproc_f .and. nproc_c .gt. 1) then
        weightall_c = weightall_f
        esav_c = esav_f
        efall_c = efall_f
    end if
#endif

    esav_f = esav_f/weightall_f
    efall_f = efall_f/weightall_f

    esav_c = esav_c/weightall_c
    efall_c = efall_c/weightall_c

    esav_f = dsqrt(esav_f - efall_f**2)/dsqrt(nsamp_f - 1.d0)
    esav_c = dsqrt(esav_c - efall_c**2)/dsqrt(nsamp_c - 1.d0)

    eall(:, nmin:nmax) = 0.d0
    do kk = nmin, nmax
        do k = 1, nbin
            do j = 1, 3
                eall(j, kk) = eall(j, kk) + ef(kk - ndimp, j, k)
            end do
        end do
    end do
    eall(:, nmin:nmax) = eall(:, nmin:nmax)/weightall_c
#ifdef PARALLEL
    if (nproc_c .gt. 1) then
        call reduce_base_real(3*np, eall(:, nmin), comm_c, -1)
    end if
#endif

    !       calculation normalization for preconditioning
    !            (diagonal elements of SR matrix)
    do kk = nmin, nmax
        okav(kk) = eall(2, kk)
    end do
    !      calculation fk
    if (nproc_c .eq. nproc_f) then
        eta(nmin:nmax) = 0.d0
        derr(nmin:nmax) = 0.d0
        do kk = nmin, nmax
            do kmain = 1, nbin
                weight = efenergy(1 + ipc, kmain)
                enerboot = (efall_c*weightall_c - efenergy(1, kmain))/(weightall_c - weight)
                do jj = 1, 3
                    ekk(jj) = (eall(jj, kk)*weightall_c - ef(kk - ndimp, jj, kmain))          &
                            & /(weightall_c - weight) !wt(jj)
                end do
                !       now average correlation function on the given bootstrap
                givenmeaskk = -ekk(1) + 2.d0*(-ekk(3) + enerboot*ekk(2))
                if (normcorr .ge. 0) then
                    fk(kmain, kk) = givenmeaskk
                else
                    !       neglect the non estensive fluctuations of Pulay
                    fk(kmain, kk) = -ekk(1)
                end if
                eta(kk) = eta(kk) + givenmeaskk
                derr(kk) = derr(kk) + givenmeaskk**2
            end do
        end do
        eta(nmin:nmax) = eta(nmin:nmax)/nsamp_c
        derr(nmin:nmax) = derr(nmin:nmax)/nsamp_c
#ifdef PARALLEL
        if (nproc_c .gt. 1) then
            call reduce_base_real(ieskin, eta(nmin), comm_c, -1)
            call reduce_base_real(ieskin, derr(nmin), comm_c, -1)
        end if
#endif
    else
        eta(nmin:nmax) = 0.d0
        do kk = nmin, nmax
            do kmain = 1, nbin
                weight = efenergy(1 + ipc, kmain)
                enerboot = (efall_c*weightall_c - efenergy(1, kmain))/(weightall_c - weight)
                do jj = 1, 3
                    ekk(jj) = (eall(jj, kk)*weightall_c - ef(kk - ndimp, jj, kmain))          &
                            & /(weightall_c - weight) !wt(jj)
                end do
                !       now average correlation function on the given bootstrap
                givenmeaskk = -ekk(1) + 2.d0*(-ekk(3) + enerboot*ekk(2))
                if (normcorr .ge. 0.d0) then
                    fk(kmain, kk) = givenmeaskk
                else
                    fk(kmain, kk) = -ekk(1)
                end if
                eta(kk) = eta(kk) + givenmeaskk
            end do
        end do
        eta(nmin:nmax) = eta(nmin:nmax)/nsamp_c
#ifdef PARALLEL
        if (nproc_c .gt. 1) then
            call reduce_base_real(ieskin, eta(nmin), comm_c, -1)
        end if
#endif
    end if
    !      Now calculation matrix cov using the prefactor of the
    !      Jackknife

    cost = dsqrt((nsamp_c - 1.d0)/nsamp_c)

    if (normcorr .ge. 0.d0) then
        do j = nmin, nmax
            do i = 1, nbin
                fk(i, j) = cost*(fk(i, j) - eta(j))
            end do
        end do
    else
        do j = nmin, nmax
            do i = 1, nbin
                fk(i, j) = cost*(fk(i, j) + eall(1, j))
            end do
        end do
    end if

    if (nproc_f .ne. nproc_c) then
        !     Recomputing eall
        eall(:, nmin:nmax) = 0.d0
        do kk = nmin, nmax
            do k = 1, nbin
                do j = 1, 3
                    eall(j, kk) = eall(j, kk) + ef(kk - ndimp, j, k)
                end do
            end do
        end do ! enddo kk
        eall(:, nmin:nmax) = eall(:, nmin:nmax)/weightall_f
#ifdef PARALLEL
        if (nproc_f .gt. 1) then
            call reduce_base_real(3*np, eall(:, nmin), comm_f, -1)
        end if
#endif
        eta(nmin:nmax) = 0.d0
        derr(nmin:nmax) = 0.d0
        do kk = nmin, nmax
            do kmain = 1, nbin
                weight = efenergy(1 + ipc, kmain)
                enerboot = (efall_f*weightall_f - efenergy(1, kmain))/(weightall_f - weight)
                do jj = 1, 3
                    ekk(jj) = (eall(jj, kk)*weightall_f - ef(kk - ndimp, jj, kmain))          &
                            & /(weightall_f - weight) !wt(jj)
                end do
                !       now average correlation function on the given bootstrap
                givenmeaskk = -ekk(1) + 2.d0*(-ekk(3) + enerboot*ekk(2))
                eta(kk) = eta(kk) + givenmeaskk
                derr(kk) = derr(kk) + givenmeaskk**2
            end do
        end do
        eta(nmin:nmax) = eta(nmin:nmax)/nsamp_f
        derr(nmin:nmax) = derr(nmin:nmax)/nsamp_f
#ifdef PARALLEL
        if (nproc_f .gt. 1) then
            call reduce_base_real(ieskin, eta(nmin), comm_f, -1)
            call reduce_base_real(ieskin, derr(nmin), comm_f, -1)
        end if
#endif
    end if

    cost = nsamp_f - 1.d0
    do kk = nmin, nmax
        derr(kk) = dsqrt(cost*(derr(kk) - eta(kk)**2))
    end do

    do kk = nmin, nmax
        eta(kk) = -eall(1, kk) + 2.d0*(-eall(3, kk) + efall_f*eall(2, kk))
    end do

    return
end
