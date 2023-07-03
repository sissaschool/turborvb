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

program zerojasfort10
    use allio
    implicit none
    real(8) Rx, Ry, Rz
    integer i, j, inputyes
    logical leave1body

    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'zerojasfort10'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added

    leave1body = .false.
    if (trim(str) .eq. "leave1body") leave1body = .true.

    call default_allocate
    open (unit=10, file='fort.10', status='old', form='formatted')
    call read_fort10(10)
    close (10)

    if (iessz) then
        write (6, *) ' Do you want to vanish also the charge Jastrow? yes/no 1/0'
        read (5, *) inputyes
    else
        inputyes = 1
    end if

    if (inputyes .eq. 1) then
        if (contractionj .ne. 0) then
            do i = 1, ipj*nelorbj_c
                do j = 1, ipj*nelorbj_c
                    if (.not. leave1body .or. (.not. orbcostn(i) .and. .not. orbcostn(j))) then
                        jasmat_c((j - 1)*ipj*nelorbj_c + i) = 0.d0
                    end if
                end do
            end do
        else
            do i = 1, ipj*nelorbj
                do j = 1, ipj*nelorbj
                    if (.not. leave1body .or. (.not. orbcostn(i) .and. .not. orbcostn(j))) then
                        jasmat((j - 1)*ipj*nelorbj + i) = 0.d0
                    end if
                end do
            end do
        end if
    end if

    if (iessz) then
        if (contractionj .ne. 0) then
            do i = 1, nelorbj_c
                do j = 1, nelorbj_c
                    if (.not. leave1body .or. (.not. orbcostn(i) .and. .not. orbcostn(j))) then
                        jasmatsz_c((j - 1)*nelorbj_c + i) = 0.d0
                    end if
                end do
            end do
        else
            do i = 1, nelorbj
                do j = 1, nelorbj
                    if (.not. leave1body .or. (.not. orbcostn(i) .and. .not. orbcostn(j))) then
                        jasmatsz((j - 1)*nelorbj + i) = 0.d0
                    end if
                end do
            end do
        end if
    end if

    open (unit=10, file='fort.10_new', status='unknown'                &
            &, form='formatted')

    call write_fort10(10)

    close (10)

end
