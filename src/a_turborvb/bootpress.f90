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

subroutine bootpress(nbin, efenergy, e, eta, derr                     &
        &, etanp, derrnp, LBox, omega, ris, rank, nproc)
    use constants, only: ipc
    implicit none

    ! *** Does BootStrap to evalute force and its error for each bin
    !  energy !  derrivative of the ener! log derivative wf! energy x log de

    !******* etanp= Pressure without pulay contributes
    !
    integer k, j, jj, kk, ieskin, nbin, kmain, rank, nproc, nbins, ierr
    real*8 e(3, *), et(4), eall(4), eta, givenmeas, derr, omega, ris(4)      &
            &, efenergy(1 + ipc, *), weight, weightall, LBox, etanp, derrnp          &
            &, summpi(5), outmpi(5)

#ifdef PARALLEL
    include 'mpif.h'
#endif

    !       initialization opara e opars

    eta = 0.d0
    derr = 0.d0
    etanp = 0.d0
    derrnp = 0.d0
    call dscalzero(4, 0.d0, eall, 1)

    weightall = 0.d0

    do k = 1, nbin
        do j = 2, 4
            eall(j) = eall(j) + e(j - 1, k)
        end do
        ! BUG in scaling pressure, the Pulay was transfomed in Hartree while the HF contr. was in
        ! Ry. Then both contributions were divided by the same factor. Thus the
        ! Pulay in pressure before 3 Sept. 2014 was underestimated by a factor 2
        ! always. That's why there were always almost the same.
        ! The calculation was instead OK when computing pressure at fixed conf.
        ! with forcevmc.sh or forcefn.sh
        !        eall(1)=eall(1)+efenergy(1,k)*ris(2)
        eall(1) = eall(1) + efenergy(1, k)
        weightall = weightall + efenergy(1 + ipc, k)
    end do

#ifdef PARALLEL
    summpi(1:4) = eall(1:4)
    summpi(5) = weightall
    call mpi_allreduce(summpi, outmpi, 5, MPI_DOUBLE_PRECISION          &
&  , MPI_SUM, MPI_COMM_WORLD, ierr)
    eall(1:4) = outmpi(1:4)
    weightall = outmpi(5)
#endif

    do j = 1, 4
        eall(j) = eall(j)/weightall
    end do

    !       write(*,*) 'eall = ',eall(1),eall(2),eall(3),eall(4)

    !       if(LBox.ne.1.d0) write(6,*) ' eall 2 scaled =',eall(2)*LBox

    !       Jet Lag evaluation of error bar with bin length lbin=100

    do kmain = 1, nbin

        do jj = 2, 4
            !w(jj,j)*e(jj,j)
            et(jj) = e(jj - 1, kmain)
        end do
        !        BUG ???   Why the following scaling from Ry-H while everything
        !        is done in Ry?
        et(1) = efenergy(1, kmain)
        weight = efenergy(1 + ipc, kmain)

        !        write(6,*) ' et read ',et(1)/weight,et(2)/weight,et(3)/weight
        !    1,et(4)/weight

        do jj = 1, 4
            !wt(jj)
            et(jj) = (eall(jj)*weightall - et(jj))/(weightall - weight)
        end do

        !       now average correlation function on the given bootstrap

        givenmeas = -(et(2) + 2.d0*(et(4) - et(1)*et(3)))

        !       write(6,*) ' givenmeas =',givenmeas,-et(2)

        eta = eta + givenmeas
        derr = derr + givenmeas**2

        etanp = etanp - et(2)
        derrnp = derrnp + et(2)**2

    end do

#ifdef PARALLEL
    summpi(1) = eta
    summpi(2) = derr
    summpi(3) = etanp
    summpi(4) = derrnp
    call mpi_allreduce(summpi, outmpi, 4, MPI_DOUBLE_PRECISION           &
&  , MPI_SUM, MPI_COMM_WORLD, ierr)
    eta = outmpi(1)
    derr = outmpi(2)
    etanp = outmpi(3)
    derrnp = outmpi(4)
#endif
    nbins = nbin*nproc

    eta = eta/nbins
    derr = derr/nbins

    etanp = etanp/nbins
    derrnp = derrnp/nbins

    !      write(6,*) 'etanp derrnp =',eta,derr,etanp,derrnp

    derrnp = dsqrt(dble(nbins - 1)**2/dble(nbins)                        &
                  & *(derrnp - etanp**2))

    derr = dsqrt(dble(nbins - 1)**2/dble(nbins)*(derr - eta**2))
    eta = -(eall(2) + 2.d0*(eall(4) - eall(1)*eall(3)))
    etanp = -eall(2)

    !      scaling for unit Hartree at the end for all (after 3 Sept. 2014).
    eta = eta*ris(2)
    etanp = etanp*ris(2)
    derr = derr*ris(2)
    derrnp = derrnp*ris(2)

    return
end
