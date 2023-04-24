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

program fort10toxyz
    use allio
    use atom_names
    implicit none
    integer c, n, i, j, countpar
    integer, allocatable :: intzeta(:)
    double precision, allocatable :: new_rion(:, :)
    double precision, allocatable :: distanza(:), tmp(:)
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'plotlambda'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added
    write (*, *) " * * * READ fort.10 and make xyz file * * * "
    open (unit=10, file='fort.10', form='formatted', status='unknown')
    call default_allocate
    call read_fort10(10)

    call eval_iond(iond, rion, nion, LBox, psip, iond_cart)

    open (unit=11, file="lambda.dat", form="formatted", status="unknown")
    open (unit=12, file="lambdaj.dat", form="formatted", status="unknown")
    if (iessz)                                                         &
            &open (unit=13, file="lambdajsz.dat", form="formatted"                &
            &, status="unknown")
    allocate (distanza(nnozero_c))

    distanza(:) = -1.d0

    do i = 1, nnozero_c
        iy = (nozero_c(i) - 1)/nelorb_c + 1
        ix = nozero_c(i) - (iy - 1)*nelorb_c
        if (iy .le. nelorb_c .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
            distanza(i) = iond(kiontot(iy) + (kiontot(ix) - 1)*nion)
        end if
    end do

    call dsortx(distanza, 1, nnozero_c, ipsip)

    countpar = 0
    do i = 1, nnozero_c
        if (distanza(i) .ge. 0.d0) then
            countpar = countpar + 1
            iy = (nozero_c(ipsip(i)) - 1)/nelorb_c + 1
            ix = nozero_c(ipsip(i)) - (iy - 1)*nelorb_c
            if (ipc .eq. 1) then
                write (11, '(2f12.6,3i12)') real(distanza(i)), real(scale_c(ipsip(i))), ix, iy, countpar
            else
                write (11, '(3f12.6,3i12)') real(distanza(i)), real(scale_c(2*ipsip(i) - 1)), &
                        & real(scale_c(2*ipsip(i))), ix, iy, countpar
            end if
        end if
    end do

    deallocate (distanza)

    allocate (distanza(nnozeroj_c))
    distanza = -1.d0

    do i = 1, nnozeroj_c
        iy = (nozeroj_c(i) - 1)/(ipj*nelorbj_c) + 1
        ix = nozeroj_c(i) - (iy - 1)*nelorbj_c*ipj
        if (.not. orbcostn(ix) .and. .not. orbcostn(iy)) then
            distanza(i) = iond(kiontotj(iy) + (kiontotj(ix) - 1)*nion)
        end if
    end do

    call dsortx(distanza, 1, nnozeroj_c, ipsip)

    countpar = 0
    do i = 1, nnozeroj_c
        if (distanza(i) .ge. 0.d0) then
            countpar = countpar + 1
            iy = (nozeroj_c(ipsip(i)) - 1)/(ipj*nelorbj_c) + 1
            ix = nozeroj_c(ipsip(i)) - (iy - 1)*ipj*nelorbj_c
            write (12, *) real(distanza(i)), real(scalej_c(ipsip(i))), ix, iy, countpar
        end if
    end do

    if (iessz) then
        distanza = -1.d0
        do i = 1, nnozeroj_c
            iy = (nozeroj_c(i) - 1)/nelorbj_c + 1
            ix = nozeroj_c(i) - (iy - 1)*nelorbj_c
            if (.not. orbcostn(ix) .and. .not. orbcostn(iy)) then
                distanza(i) = iond(kiontotj(iy) + (kiontotj(ix) - 1)*nion)
            end if
        end do

        call dsortx(distanza, 1, nnozeroj_c, ipsip)

        do i = 1, nnozeroj_c
            if (distanza(i) .ge. 0.d0) then
                iy = (nozeroj_c(ipsip(i)) - 1)/nelorbj_c + 1
                ix = nozeroj_c(ipsip(i)) - (iy - 1)*nelorbj_c
                write (13, *) real(distanza(i)), real(scalejsz_c(ipsip(i))), ix, iy
            end if
        end do

    end if

    deallocate (distanza)

    if (iessz) close (13)
    close (12)
    close (11)
    close (10)
end
