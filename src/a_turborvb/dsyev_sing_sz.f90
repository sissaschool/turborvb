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

subroutine dsyev_sing(nbin, cov, eig, tcin, psip, lwork           &
        &, iwork, ilwork, info)
    implicit none
    integer ipiv, nbin, lwork, info, i, j, nbin2, nbinm, iwork(*), ilwork
    real*8 cov(nbin, *), eig(*), psip(*), nmax, dnrm2                  &
            &, ddot, cost, sum, tcin(nbin, *), arg

    !         this subroutine returns eigenvector and eigenvalues
    !         of a matrix cov, satisfying by definition
    !         sum_j cov(i,j)=0., i.e. has a singulare eigenvector
    !         which is constant.

    arg = 4.d0*dasin(1.d0)/nbin
    cost = 1.d0/dsqrt(dble(nbin))
    call dscalzero(nbin, cost, tcin, 1)

    nbin2 = nbin/2

    if (2*nbin2 .eq. nbin) then
        do i = 2, nbin2
            do j = 1, nbin
                tcin(j, i) = dcos(arg*(i - 1)*j)
                tcin(j, i + nbin2) = dsin(arg*(i - 1)*j)
            end do
            cost = 1.d0/dnrm2(nbin, tcin(1, i), 1)
            call dscal(nbin, cost, tcin(1, i), 1)
            cost = 1.d0/dnrm2(nbin, tcin(1, i + nbin2), 1)
            call dscal(nbin, cost, tcin(1, i + nbin2), 1)
        end do
        do j = 1, nbin
            tcin(j, nbin2 + 1) = dcos(arg*nbin2*j)
        end do
        cost = 1.d0/dnrm2(nbin, tcin(1, nbin2 + 1), 1)
        call dscal(nbin, cost, tcin(1, nbin2 + 1), 1)

    else

        do i = 2, nbin2 + 1
            do j = 1, nbin
                tcin(j, i) = dcos(arg*(i - 1)*j)
                tcin(j, i + nbin2) = dsin(arg*(i - 1)*j)
            end do
            cost = 1.d0/dnrm2(nbin, tcin(1, i), 1)
            call dscal(nbin, cost, tcin(1, i), 1)
            cost = 1.d0/dnrm2(nbin, tcin(1, i + nbin2), 1)
            call dscal(nbin, cost, tcin(1, i + nbin2), 1)
        end do
    end if

    !         preconditioning before
    !          do j=1,nbin
    !             do i=1,nbin
    !             tcin(j,i)=fkav(j)*tcin(j,i)
    !             enddo
    !          enddo

    !          write(6,*) ' zero root before dsyev ',ipiv

    !          write(6,*) ' raw/column ipip '
    !          do i=1,nbin
    !          write(6,*) i,cov(i,ipiv),cov(ipiv,i)
    !          enddo

    !      Input L =cov , Cholesky S= L L^+
    !          now transform the matrix L^+ T

    call dgemm('T', 'N', nbin, nbin, nbin, 1.d0, cov, nbin, tcin, nbin     &
            &, 0.d0, psip, nbin)

    call dgemm('T', 'N', nbin, nbin, nbin, 1.d0, psip, nbin, psip, nbin    &
            &, 0.d0, cov, nbin)

    !          this transformed matrix should have the first
    !          index=0

    !      write(6,*) ' are really zero ?? '
    !      do i=1,nbin
    !      write(6,*) i,cov(i,1),cov(1,i)
    !      enddo

    !      write(6,*) ' full transformed matrix '
    !      do i=1,nbin
    !         do j=i,nbin
    !         write(6,*)  i,j,cov(i,j)
    !         enddo
    !      enddo

    nbinm = nbin - 1

    call dsyevd('V', 'U', nbinm, cov(2, 2), nbin, eig(2), psip, lwork        &
            &, iwork, ilwork, info)

    if (eig(2) .le. 0.) then

        !c         preconditioning after
        do j = 1, nbin
            do i = 1, nbin
                write (6, *) i, j, cov(i, j), cov(j, i)
                !             tcin(j,i)=tcin(j,i)/fkav(j)
            end do
        end do

    end if

    !          now transform back the eigenvectors
    call dgemm('N', 'N', nbin, nbinm, nbinm, 1.d0, tcin(1, 2), nbin          &
            &, cov(2, 2), nbin, 0.d0, psip, nbin)

    call dcopy(nbin*nbinm, psip, 1, cov(1, 2), 1)

    eig(1) = 0.d0
    cost = 1.d0/dsqrt(dble(nbin))
    call dscalzero(nbin, cost, cov, 1)

    !          write(6,*) ' check orthogonality '
    !          do i=1,nbin
    !            do j=i,nbin
    !            write(6,*) i,j,ddot(nbin,cov(1,i),1,cov(1,j),1)
    !            enddo
    !          enddo
    !          stop
    return
end
