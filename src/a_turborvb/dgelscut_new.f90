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

subroutine dgelscut(sov, nmat, npm, psip, epsdgel, ipsip, np, info    &
        &, scalpar, lwork, epsmach, indexpar, rank, nprocopt, rankopt, commopt_mpi)
    use allio, only: nproc_diag
    implicit none
    integer ip, nmat, i, j, npm, n1, n2, n3, n4, n5, info, rank, ipsip(*)    &
            &, ind, lwork, iscra, k, ierr, imax, imaxg, indi, indj, np &
            &, indexpar(*), indscali, nprocopt, rankopt, commopt_mpi
    real*8 sov(npm, npm, *), psip(*), epsdgel, cost, epsmach, ratoff      &
            &, costmax, costmaxg, scalpar(*), eigmin, eigs, eigb, eigs_min, eigb_min, maxsr

    !        epsmach=1d-12
    n1 = nmat + 1
    n2 = n1 + nmat
    n3 = n2 + nmat
    n4 = n3 + nmat
    n5 = n4 + nmat
    !        preconditioning the matrix
    maxsr = 0.d0
    do i = 1, nmat
        costmax = abs(sov(i, i, 1)*scalpar(i))
        if (costmax .gt. maxsr) maxsr = costmax
    end do
    indj = 0
    do i = 1, nmat
        if (abs(sov(i, i, 1)*scalpar(i)) .gt. epsmach*maxsr) then
            psip(n1 + i - 1) = 1.d0/dsqrt(dabs(sov(i, i, 1)))
        else
            !        write(6,*) ' Warning eliminated =',i,sov(i,i,1)
            psip(n1 + i - 1) = 0.d0
            ipsip(i) = 0
            indj = indj + 1
        end if
    end do

    indi = 0
    do i = 1, nmat
        if (ipsip(i) .eq. 1) then
            indi = indi + 1
            indj = 0
            do j = 1, nmat
                if (ipsip(j) .eq. 1) then
                    indj = indj + 1
                    sov(indi, indj, 2) = sov(i, j, 1)*psip(n1 + i - 1)*psip(n1 + j - 1)
                end if
            end do
        end if
    end do

    !        first step diagonalize the matrix
    if (indi .gt. 1) then
        if (nproc_diag .gt. 1) then
            call dsyev_my('N', 'L', indi, sov(1, 1, 2), npm, psip, info     &
                    &, nprocopt, rankopt, commopt_mpi)
        else
            call dsyev('N', 'L', indi, sov(1, 1, 2), npm, psip, psip(n2), lwork, info)
        end if
        if (rank .eq. 0) write (6, *) ' Lowest/Max  eigenvalues SR mat ', psip(1), psip(indi)

        eigs = psip(1)
        eigb = psip(2)
        eigs_min = psip(1)
        eigb_min = psip(2)

        eigmin = psip(1)

    else
        if (rank .eq. 0) write (6, *) ' Warning no par to eliminate =', indi
        eigmin = 1.d0

    end if

    if (eigmin .lt. epsdgel .and. indi .gt. 1) then
        do while (eigmin .lt. epsdgel)

            !        definition matrix for the chosen raw
            indi = 0
            do i = 1, nmat
                if (ipsip(i) .eq. 1) then
                    indi = indi + 1
                    indj = 0
                    do j = 1, nmat
                        if (ipsip(j) .eq. 1) then
                            indj = indj + 1
                            sov(indi, indj, 2) = sov(i, j, 1)*psip(n1 + i - 1)*psip(n1 + j - 1)
                        end if
                    end do
                end if
            end do
            do j = 1, indi
                sov(j, j, 2) = 1.d0
            end do

            if (nproc_diag .gt. 1) then
                call dsyev_my('V', 'L', indi, sov(1, 1, 2), npm, psip, info &
                        &, nprocopt, rankopt, commopt_mpi)
            else
                call dsyev('V', 'L', indi, sov(1, 1, 2), npm, psip, psip(n2), lwork, info)
            end if

            eigs_min = psip(1)
            eigb_min = psip(2)

            !      write(6,*) ' before findzero '

            call findzero(eigs_min, eigb_min, npm, indi, sov(1, 1, 2), psip        &
                    &, psip(n2), psip(n3), psip(n4), psip(n5))

            !       write(6,*) ' after findzero '
            costmax = eigs_min
            imax = 0
            indj = 0
            !         write(6,*) ' Eigenvalues inside new '

            !    In principle we can eliminate all parameters < epsdgel
            !    with the constraint that they are independent
            do k = 1, nmat
                if (ipsip(k) .eq. 1) then
                    indj = indj + 1
                    if (psip(n5 + indj - 1) .ge. costmax) then
                        imax = k
                        costmax = psip(n5 + indj - 1)
                    end if
                    !         write(6,*) indj,psip(n5+indj-1)
                end if
            end do

            indscali = indexpar(imax)

            if (rank .eq. 0) then
                if (indscali .eq. 0) then
                    write (6, *) 'eliminated collective par', imax
                else
                    write (6, *) 'eliminated normal  par', indscali
                end if
            end if

            if (imax .eq. 0) then
                if (rank .eq. 0) write (6, *) ' Some error in findzero ', eigs_min, eigb_min
#ifdef PARALLEL
                call mpi_finalize(ierr)
#endif
                stop
            else
                ipsip(imax) = 0
                eigmin = costmax
            end if

        end do

        if (rank .eq. 0) write (6, *) ' New  Lowest eigenval ', eigmin
    end if
    if (rank .eq. 0) write (6, *) ' Fixed parameters '
    ind = 0
    do i = 1, nmat
        ipsip(i + np) = ipsip(i)
        if (ipsip(i) .eq. 0) then
            ind = ind + 1
            indscali = indexpar(i)
            if (rank .eq. 0) then
                if (indscali .eq. 0) then
                    write (6, *) ind, i
                else
                    write (6, *) ind, indscali
                end if
            end if
        end if
        ! only the regul
        if (scalpar(i) .gt. 0.d0 .and. ipsip(i) .ne. 0) then
            ipsip(i + np) = 2
        end if
    end do

    return
end
