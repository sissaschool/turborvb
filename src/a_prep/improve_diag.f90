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
! INPUT: molecorb, hamiltl, oversl = approximate eigenvector and eigenvalues
! of the KS hamiltonian.
! OUTPUT: molecorb, eigmol = orthogonal eigenvectors and more accurate eigenvalues.
! It is used also molecorb_part dynamically.
!

subroutine improvediag

    use constants, only: ipc, zzero, zone
    use allio, only: commrep_mpi, rankrep, rank, psip
    use setup, only: overs, oversdo, molecorb, molecorbdo, hamilt, hamiltdo, &
                     nelorbu, nelocc, neloccdo, bands, yeslsda, eigmol, eigmoldo, orthodiag
#ifdef __SCALAPACK
    use descriptors
    use setup, only: desch, oversl, oversldo, molecorbl, molecorbldo, hamiltl, hamiltldo
#endif

    implicit none

    integer :: ir, ic, ierr, info, i, j
    real(8), dimension(:, :), allocatable :: molecorb_part
    real(8), external :: ddot
    complex(8), external :: zdotc_, zdotu
#ifdef PARALLEL
    include "mpif.h"
#endif

    info = 0 ! error flag for orthogonalization
    if (rank .eq. 0) write (6, '(a)') ' # WARNING: orthogonalization of the orbitals! '

#ifdef __SCALAPACK
    if (ipc .eq. 1) then
        call graham_scalapack(molecorb, oversl, psip, nelorbu, nelorbu, nelocc, &
                              info, descla, desch, size(oversl, 1), rankrep)
        if (yeslsda) &
            call graham_scalapack(molecorbdo, oversl, psip, nelorbu, nelorbu, neloccdo, &
                                  info, descla, desch, size(oversl, 1), rankrep)
    else
        call graham_scalapack_complex(molecorb, oversl, psip, nelorbu, nelorbu, nelocc, &
                                      info, descla, desch, size(oversl, 1)/2, rankrep)
        call graham_scalapack_complex(molecorbdo, oversldo, psip, nelorbu, nelorbu, neloccdo, &
                                      info, descla, desch, size(oversldo, 1)/2, rankrep)
    end if
#else
    if (ipc .eq. 1) then
        call graham(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
        if (yeslsda) &
            call graham(molecorbdo, overs, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
    else
        call graham_complex(molecorb, overs, nelorbu, psip, nelorbu, nelorbu, nelocc, info)
        call graham_complex(molecorbdo, oversdo, nelorbu, psip, nelorbu, nelorbu, neloccdo, info)
    end if
#endif

#ifdef PARALLEL
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
!$omp barrier

#ifdef DEBUG
    if (rank .eq. 0) write (6, *) ' Check orthogonality eigenvectors: '
    allocate (molecorb_part(ipc*nelorbu, nelocc + neloccdo))
    molecorb_part = 0.d0

    if (ipc .eq. 1) then
        call dgemm('N', 'N', nelorbu, nelocc, nelorbu, 1.d0, overs &
                   , nelorbu, molecorb, nelorbu, 0.d0, molecorb_part, nelorbu)
        if (yeslsda) then
            call dgemm('N', 'N', nelorbu, neloccdo, nelorbu, 1.d0, overs &
                       , nelorbu, molecorbdo, nelorbu, 0.d0, molecorb_part(1, nelocc + 1), nelorbu)
        end if
    else
        call zgemm('N', 'N', nelorbu, nelocc, nelorbu, zone, overs &
                   , nelorbu, molecorb, nelorbu, zzero, molecorb_part, nelorbu)
        call zgemm('N', 'N', nelorbu, neloccdo, nelorbu, zone, oversdo &
                   , nelorbu, molecorbdo, nelorbu, zzero, molecorb_part(1, nelocc + 1), nelorbu)
    end if

    do i = 1, nelocc
        do j = i, nelocc
            if (ipc .eq. 1) then
                if (rank .eq. 0) then
                    if (yeslsda) then
                        write (6, *) i, j, dabs(ddot(nelorbu, molecorb(1, i), 1, molecorb_part(1, j), 1)), &
                            dabs(ddot(nelorbu, molecorbdo(1, i), 1, molecorb_part(1, j + nelocc), 1))
                    else
                        write (6, *) i, j, dabs(ddot(nelorbu, molecorb(1, i), 1, molecorb_part(1, j), 1))
                    end if
                end if
            else
                if (rank .eq. 0) write (6, *) i, j, abs(zdotc_(nelorbu, molecorb(1, i), 1, molecorb_part(1, j), 1)), &
                    abs(zdotc_(nelorbu, molecorbdo(1, i), 1, molecorb_part(1, j + nelocc), 1))
            end if
        end do
    end do
    deallocate (molecorb_part)
#endif

#ifdef __SCALAPACK

    eigmol(1:bands) = 0.d0

    if (descla(lambda_node_) > 0) then

        allocate (molecorb_part(ipc*descla(nlax_), descla(nlax_)))
        molecorb_part = 0.d0
        !
        ! copy the orthogonalized eigenvectors in the distributed matrix molecorbl
        !
        ir = descla(ilar_)
        ic = descla(ilac_)
        do j = 1, descla(nlac_)
            do i = 1, descla(nlar_)
                if ((j + ic - 1) <= bands) then
                    if (ipc .eq. 1) then
                        molecorbl(i, j) = molecorb((i + ir - 1), (j + ic - 1))
                    else
                        molecorbl(2*i - 1:2*i, j) = molecorb(2*(i + ir - 1) - 1:2*(i + ir - 1), (j + ic - 1))
                    end if
                end if
            end do
        end do
        !
        ! compute H |\psi> = E |\psi> and put in the scratch distributed matrix molecorb_part
        !
        if (ipc .eq. 1) then
            call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.0d0, hamiltl, 1, 1,&
                 & desch, molecorbl, 1, 1, desch, 0.0d0, molecorb_part, 1, 1, desch)
        else
            call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, zone, hamiltl, 1, 1,&
                 & desch, molecorbl, 1, 1, desch, zzero, molecorb_part, 1, 1, desch)
        end if
        !
        ! recompute the eigenvalues improved after orthogonalization
        !
        ir = descla(ilar_)
        ic = descla(ilac_)
        do j = 1, descla(nlac_)
            if ((j + ic - 1) <= bands) then
                do i = 1, descla(nlar_)
                    if (ipc .eq. 1) then
                        eigmol(j + ic - 1) = eigmol(j + ic - 1) + molecorb_part(i, j)*molecorb(i + ir - 1, j + ic - 1)
                    else
                        ! molecorb * dconj(molecorb)
                        eigmol(j + ic - 1) = eigmol(j + ic - 1) + &
                                             ((molecorb_part(2*i - 1, j)*molecorb(2*(i + ir - 1) - 1, j + ic - 1)) + & ! Re
                                              (molecorb_part(2*i, j)*molecorb(2*(i + ir - 1), j + ic - 1))) ! Im
                    end if
                end do
            end if
        end do

        !ir = descla( ilar_ )
        !ic = descla( ilac_ )
        !do j=1,descla( nlac_ )
        !   do i=1,descla( nlar_ )
        !      IF( (j+ic-1) <= bands ) THEN
        !         if(ipc.eq.1) then
        !            molecorbl( i, j ) =  molecorb( (i+ir-1), (j+ic-1) )
        !         else
        !            molecorbl( 2*i-1:2*i, j ) =  molecorb( 2*(i+ir-1)-1:2*(i+ir-1), (j+ic-1) )
        !         endif
        !      END IF
        !   enddo
        !enddo

        !ir = descla( ilar_ )
        !ic = descla( ilac_ )
        !do j=1,descla( nlac_ )
        !   IF( (j+ic-1) <= bands ) THEN
        !      do i=1,descla( nlar_ )
        !         eigmol(j+ic-1)=eigmol(j+ic-1)+molecorb_part(i,j)*molecorb(i+ir-1,j+ic-1)
        !      enddo
        !   END IF
        !enddo

        deallocate (molecorb_part)
    end if
#ifdef PARALLEL
    call reduce_base_real(bands, eigmol, commrep_mpi, -1)
#endif
    !
    ! same story for the spin down part
    !
    if (yeslsda .or. ipc .eq. 2) then

        eigmoldo(1:bands) = 0.d0

        if (descla(lambda_node_) > 0) then

            allocate (molecorb_part(ipc*descla(nlax_), descla(nlax_)))
            molecorb_part = 0.d0

            ir = descla(ilar_)
            ic = descla(ilac_)
            do j = 1, descla(nlac_)
                do i = 1, descla(nlar_)
                    if ((j + ic - 1) <= bands) then
                        if (ipc .eq. 1) then
                            molecorbldo(i, j) = molecorbdo((i + ir - 1), (j + ic - 1))
                        else
                            molecorbldo(2*i - 1:2*i, j) = molecorbdo(2*(i + ir - 1) - 1:2*(i + ir - 1), (j + ic - 1))
                        end if
                    end if
                end do
            end do

            if (ipc .eq. 1) then
                call PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.0d0, hamiltldo, 1, 1,&
                     & desch, molecorbldo, 1, 1, desch, 0.0d0, molecorb_part, 1, 1, desch)
            else
                call PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, zone, hamiltldo, 1, 1,&
                     & desch, molecorbldo, 1, 1, desch, zzero, molecorb_part, 1, 1, desch)
            end if

            ir = descla(ilar_)
            ic = descla(ilac_)
            do j = 1, descla(nlac_)
                if ((j + ic - 1) <= bands) then
                    do i = 1, descla(nlar_)
                        if (ipc .eq. 1) then
                            eigmoldo(j + ic - 1) = eigmoldo(j + ic - 1) + molecorb_part(i, j)*molecorbdo(i + ir - 1, j + ic - 1)
                        else
                            ! molecorb * dconj(molecorb)
                            eigmoldo(j + ic - 1) = eigmoldo(j + ic - 1) + &
                                                   ((molecorb_part(2*i - 1, j)*molecorbdo(2*(i + ir - 1) - 1, j + ic - 1)) + & ! Re
                                                    (molecorb_part(2*i, j)*molecorbdo(2*(i + ir - 1), j + ic - 1))) ! Im
                        end if
                    end do
                end if
            end do

            !if(ipc.eq.1) then
            !   CALL PDGEMM('N', 'N', nelorbu, nelorbu, nelorbu, 1.0d0, hamiltldo, 1, 1,&
            !        & desch, molecorbldo, 1, 1, desch, 0.0d0, molecorb_part,1,1,desch)
            !else
            !   CALL PZGEMM('N', 'N', nelorbu, nelorbu, nelorbu, zone, hamiltldo, 1, 1,&
            !        & desch, molecorbldo, 1, 1, desch, zzero, molecorb_part,1,1,desch)
            !endif

            !ir = descla( ilar_ )
            !ic = descla( ilac_ )
            !do j=1,descla( nlac_ )
            !   IF( (j+ic-1) <= bands ) THEN
            !      do i=1,descla( nlar_ )
            !         eigmoldo(j+ic-1)=eigmoldo(j+ic-1)+molecorb_part(i,j)*molecorbdo(i+ir-1,j+ic-1)
            !      enddo
            !   END IF
            !enddo

            deallocate (molecorb_part)

        end if
#ifdef PARALLEL
        call reduce_base_real(bands, eigmoldo, commrep_mpi, -1)
#endif
    end if

#else

    allocate (molecorb_part(ipc*nelorbu, bands))

    if (ipc .eq. 1) then
        call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, hamilt, nelorbu&
                &, molecorb, nelorbu, 0.d0, molecorb_part, nelorbu)
        do i = 1, bands
            eigmol(i) = ddot(nelorbu, molecorb(1, i), 1, molecorb_part(1, i), 1)
        end do
    else
        call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, hamilt, nelorbu&
                &, molecorb, nelorbu, zzero, molecorb_part, nelorbu)
        do i = 1, bands
            eigmol(i) = zdotc_(nelorbu, molecorb(1, i), 1, molecorb_part(1, i), 1)
        end do

    end if

    if (yeslsda .or. ipc .eq. 2) then
        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbu, bands, nelorbu, 1.d0, hamiltdo, nelorbu&
                    &, molecorbdo, nelorbu, 0.d0, molecorb_part, nelorbu)
            do i = 1, bands
                eigmoldo(i) = ddot(nelorbu, molecorbdo(1, i), 1, molecorb_part(1, i), 1)
            end do
        else
            call zgemm('N', 'N', nelorbu, bands, nelorbu, zone, hamiltdo, nelorbu&
                    &, molecorbdo, nelorbu, zzero, molecorb_part, nelorbu)
            do i = 1, bands
                eigmoldo(i) = zdotc_(nelorbu, molecorbdo(1, i), 1, molecorb_part(1, i), 1)
            end do
        end if
    end if
    deallocate (molecorb_part)

#endif

    return

end subroutine improvediag
