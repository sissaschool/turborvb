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

subroutine initconf(nel, nelup, psiln, psisn, kel, wconf, rion, dist  &
        &, indt, nion, zeta, in1, rank, ierr, LBox, alat, rion_ref, iesrandoma, itest)
    use Cell
    implicit none

    integer nhalf, nelup, nel, i, j, ii, kk, ih, indt, nion, ion, in1, itest
    real*8 kel(3, nel, 0:indt, *)
    real*8 b, c, zeta(nion), LBox, dst
    integer cup, cdown
    real*8 wconf(*), drand1, rion(3, *)
    real*8 dist(nion, nel, *)
    integer occflag, nelupeff, ztot
    integer ist, ien, rank, ierr, d
    real*8 psiln(*), dist_kel(3), scart(3), rion_ref(3), alat
    real*8 psisn(*)
    integer, dimension(:), allocatable :: occup, occdo

    real*8 dist_max, dist_try
    integer i_seq, i_max, i_previous
    integer, dimension(:), allocatable :: ion_seq
    logical, dimension(:), allocatable :: ion_sel
    logical iesrandoma, yeslat, checkdiff

    allocate (occup(nion), occdo(nion))
    occup = 0
    occdo = 0

    if (alat .ne. 0.d0 .and. .not. iesrandoma .and. itest .ne. 2) then
        yeslat = .true.
    else
        yeslat = .false.
    end if

    if (yeslat) then
        dst = 2.d0*abs(alat)*dble(nelup)**(1.d0/3.d0) ! to avoid to generate the same conf
    else
        dst = 2.d0
    end if
    ! find a new sequence of ions in order to improve fairness
    ! of the initial el configuration
    ! the sequence maximizes the distance from one ion
    ! to the next one.
    allocate (ion_sel(nion), ion_seq(nion))
    ! ion_sel(i) true i-th ion not included yet
    ion_sel = .true.
    ! start new sequence with the same element
    ion_seq(1) = 1
    ! previous element in the new list
    i_previous = 1
    ! the first is included
    ion_sel(i_previous) = .false.
    ! number of elements already included
    i_seq = 1
    do while (i_seq .lt. nion)
        ! first value for distance
        dist_max = 0.d0
        do i = 1, nion
            ! search among those not included yet
            if (ion_sel(i)) then
                dist_try = 0.d0
                do ii = 1, 3
                    dist_try = dist_try + (rion(ii, i_previous) - rion(ii, i))**2
                end do
                ! check for the farest
                if (dist_try .ge. dist_max) then
                    dist_max = dist_try
                    i_max = i
                end if
            end if
        end do
        i_seq = i_seq + 1
        ion_seq(i_seq) = i_max
        ion_sel(i_max) = .false.
        i_previous = i_max
    end do

    ztot = 0
    do i = 1, nion
        ztot = ztot + zeta(i)
    end do

    if (nel .gt. ztot) then
        nelupeff = nel - ztot
    else
        nelupeff = 0
    end if
    ! if nelupeff gt 0 --> electron affinity calculation

    if (rank .eq. 0) then
        write (6, *) '****************************************'
        write (6, *) '****** INITIALIZATION ******************'
    end if

    ! the constraint here is that the spin up  alone are forbidden

    ! find the rion in the simulations coordinates
    !      lion=1

    do kk = 1, in1
        ih = 0
        wconf(kk + rank*in1) = 1.d0
        psisn(kk) = 0
        psiln(kk) = 0.d0

        do i = 1, nion
            occup(i) = 0
        end do

        cup = 0
        cdown = 0

        j = 0
        do i = 1, nel - nelup

            occflag = 0

            do while (occflag .eq. 0)
                ion = ion_seq(mod(j, nion) + 1)
                if (occup(ion) .lt. zeta(ion)/2 .or. nelupeff .gt. 1) then
                    kel(1, nelup + i, 0, kk) = dst*(drand1() - 0.5d0) + rion(1, ion)
                    kel(2, nelup + i, 0, kk) = dst*(drand1() - 0.5d0) + rion(2, ion)
                    kel(3, nelup + i, 0, kk) = dst*(drand1() - 0.5d0) + rion(3, ion)

                    checkdiff = .false.
                    do while (.not. checkdiff)
                        checkdiff = .true.
                        ii = 1
                        do while (checkdiff .and. ii .lt. i)
                            if (sum(abs(kel(:, nelup + ii, 0, kk) - kel(:, nelup + i, 0, kk))) .lt. 1d-6) checkdiff = .false.
                            ii = ii + 1
                        end do
                        if (.not. checkdiff) then
                            kel(1, nelup + i, 0, kk) = dst*(drand1() - 0.5d0) + rion(1, ion)
                            kel(2, nelup + i, 0, kk) = dst*(drand1() - 0.5d0) + rion(2, ion)
                            kel(3, nelup + i, 0, kk) = dst*(drand1() - 0.5d0) + rion(3, ion)
                            if (yeslat) then
                                kel(1, nelup + i, 0, kk) = alat*nint((kel(1, nelup + i, 0, kk) - rion_ref(1))/alat)&
                                        & + rion_ref(1)
                                kel(2, nelup + i, 0, kk) = alat*nint((kel(2, nelup + i, 0, kk) - rion_ref(2))/alat)&
                                        & + rion_ref(2)
                                kel(3, nelup + i, 0, kk) = alat*nint((kel(3, nelup + i, 0, kk) - rion_ref(3))/alat)&
                                        & + rion_ref(3)
                            end if
                        end if
                    end do

                    !*********** PERIODIC SYSTEMS *********************
                    if (LBox .gt. 0.d0) call ApplyPBC(kel(1, nelup + i, 0, kk), 1)
                    !**************************************************

                    occup(ion) = occup(ion) + 1
                    cdown = cdown + 1
                    occflag = 1
                end if

                j = j + 1

            end do

        end do

        do i = 1, nion
            occdo(i) = occup(i)
        end do

        j = 0
        if (nelupeff .le. 1) then
            nhalf = nelup - nelupeff
        else
            nhalf = nelup
        end if

        do i = 1, nhalf

            occflag = 0

            do while (occflag .eq. 0)

                ion = ion_seq(mod(j, nion) + 1)
                if (occup(ion) .lt. zeta(ion) .or. nelupeff .gt. 1) then
                    kel(1, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(1, ion)
                    kel(2, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(2, ion)
                    kel(3, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(3, ion)
                    checkdiff = .false.
                    do while (.not. checkdiff)
                        checkdiff = .true.
                        ii = 1
                        do while (checkdiff .and. ii .lt. i)
                            if (sum(abs(kel(:, ii, 0, kk) - kel(:, i, 0, kk))) .lt. 1d-6) checkdiff = .false.
                            ii = ii + 1
                        end do
                        if (.not. checkdiff) then
                            kel(1, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(1, ion)
                            kel(2, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(2, ion)
                            kel(3, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(3, ion)
                            if (yeslat) then
                                kel(1, i, 0, kk) = alat*nint((kel(1, i, 0, kk) - rion_ref(1))/alat) + rion_ref(1)
                                kel(2, i, 0, kk) = alat*nint((kel(2, i, 0, kk) - rion_ref(2))/alat) + rion_ref(2)
                                kel(3, i, 0, kk) = alat*nint((kel(3, i, 0, kk) - rion_ref(3))/alat) + rion_ref(3)
                            end if
                        end if
                    end do

                    !*********** PERIODIC SYSTEMS ********************
                    if (LBox .gt. 0.d0) call ApplyPBC(kel(1, i, 0, kk), 1)
                    !*************************************************

                    occup(ion) = occup(ion) + 1
                    cup = cup + 1
                    occflag = 1
                end if

                j = j + 1

            end do

        end do

        if (nelupeff .le. 1) then

            do i = nelup - nelupeff + 1, nelup

                ion = ion_seq(int(drand1()*nion) + 1)
                kel(1, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(1, ion)
                kel(2, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(2, ion)
                kel(3, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(3, ion)

                checkdiff = .false.
                do while (.not. checkdiff)
                    checkdiff = .true.
                    ii = 1
                    do while (checkdiff .and. ii .lt. i)
                        if (sum(abs(kel(:, ii, 0, kk) - kel(:, i, 0, kk))) .lt. 1d-6) checkdiff = .false.
                        ii = ii + 1
                    end do
                    if (.not. checkdiff) then
                        kel(1, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(1, ion)
                        kel(2, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(2, ion)
                        kel(3, i, 0, kk) = dst*(drand1() - 0.5d0) + rion(3, ion)
                        if (yeslat) then
                            kel(1, i, 0, kk) = alat*nint((kel(1, i, 0, kk) - rion_ref(1))/alat) + rion_ref(1)
                            kel(2, i, 0, kk) = alat*nint((kel(2, i, 0, kk) - rion_ref(2))/alat) + rion_ref(2)
                            kel(3, i, 0, kk) = alat*nint((kel(3, i, 0, kk) - rion_ref(3))/alat) + rion_ref(3)
                        end if
                    end if
                end do

                !*********** PERIODIC SYSTEMS ********************
                if (LBox .gt. 0.d0) call ApplyPBC(kel(1, i, 0, kk), 1)

                occup(ion) = occup(ion) + 1
                cup = cup + 1

            end do

        end if

        if (rank .eq. 0 .and. kk .eq. 1) then
            do i = 1, nion
                write (6, *) 'zeta, up and down ', i, zeta(i), occup(i) - occdo(i), occdo(i)
            end do

        end if

        if (nelup .ne. cup) then
            if (rank .eq. 0)                                                  &
                    &   write (6, *) 'error in up nesting', nelup, cup
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        if (nel .ne. cup + cdown) then
            if (rank .eq. 0)                                                  &
                    &   write (6, *) 'error in electrons nesting', nel, cup + cdown
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        if (LBox .le. 0.d0) then
            do ih = 1, nel
                do j = 1, nion
                    dist(j, ih, kk) = dsqrt((kel(1, ih, 0, kk) - rion(1, j))**2 + &
                            &   (kel(2, ih, 0, kk) - rion(2, j))**2 + (kel(3, ih, 0, kk) - rion(3, j))**2)
                end do
            end do
        else
            do ih = 1, nel
                do j = 1, nion
                    dist_kel(:) = kel(:, ih, 0, kk) - rion(:, j)
                    call ApplyPBC(dist_kel, 1)
                    dist(j, ih, kk) = dsqrt(sum(dist_kel(:)**2))
                end do
            end do
        end if
    end do

    ! test
    !         kel(1,1,0,1)=0.5d0
    !         kel(2,1,0,1)=0.25d0
    !         kel(3,1,0,1)=-0.3d0
    !         kel(1,2,0,1)=0.34d0
    !         kel(2,2,0,1)=0.25d0
    !         kel(3,2,0,1)=0.55d0
    !         kel(1,3,0,1)=-0.2d0
    !         kel(2,3,0,1)=0.5d0
    !         kel(3,3,0,1)=-0.2d0
    !         kel(1,4,0,1)=0.2d0
    !         kel(2,4,0,1)=0.5d0
    !         kel(3,4,0,1)=0.2d0
    !         kel(1,5,0,1)=0.1d0
    !         kel(2,5,0,1)=0.7d0
    !         kel(3,5,0,1)=0.2d0
    !         kel(1,6,0,1)=-0.1d0
    !         kel(2,6,0,1)=-0.62d0
    !         kel(3,6,0,1)=-0.2d0
    !         kel(1,7,0,1)=0.1d0
    !         kel(2,7,0,1)=-1.2d0
    !         kel(3,7,0,1)=-0.2d0
    !         kel(1,8,0,1)=0.1d0
    !         kel(2,8,0,1)=-0.7d0
    !         kel(3,8,0,1)=-0.2d0

    !     write(6,*) ' Initial conf '
    !     do kk=ist,ien
    !        do ih=1,nel
    !           do ii=1,3
    !              write(6,*) ii,ih,kel(ii,ih,0,kk)
    !           enddo
    !        enddo
    !     enddo

    deallocate (occup, occdo, ion_seq, ion_sel)

    return
end
