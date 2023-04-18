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

program find_kpoints

    ! This tool is useful to determine the total number of k-points
    ! given a certain input file. It requires the wave function fort.10
    ! and the input file and it prints out all the k-points associated
    ! to the calculation.

    use allio
    use kpoints_mod
    use cell

    implicit none
    integer, parameter :: ufort10 = 10
    integer :: i, j
    !
    character(100) name_tool
    character(20) str
    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then
        ! input the name of the file exactly as it is in /doc
        name_tool = 'find_kpoints'
        call help_online(name_tool)
        stop
    end if
    !
    ! read the input wave function in fort.10
    !
    call default_allocate
    open (unit=ufort10, file='fort.10', status='old', form='formatted', err=101)
    call read_fort10(ufort10)
    ! some useful warnings
    if (.not. yes_complex) then
        write (6, *)
        write (6, *) 'WARNING: You need a complex wave function in order to run a k-points calculations!'
        stop
    end if
    if (.not. iespbc) then
        write (6, *)
        write (6, *) 'ERROR: You need a PBC wave function in order to run a k-points calculations!'
        stop
    end if
    ! the cell has been already initialized with read_fort10
    kp_type = -1
    k1 = 0
    k2 = 0
    k3 = 0
    nk = 1
    nk1 = -1
    nk2 = -1
    nk3 = -1
    time_reversal = .false.
    skip_equivalence = .true.
    double_kpgrid = .false.
    ! k-point card must be given in input
    iflagerr = 1
    read (5, nml=kpoints, err=102)
    iflagerr = 0
    if (nk2 .eq. -1) nk2 = nk1
    if (nk3 .eq. -1) nk3 = nk1
    !
    ! find total number of k-points and
    ! allocate k-points and weigths
    !
    select case (kp_type)
    case (0)
        nk = 1
    case (1, -1)
        nk = nk1*nk2*nk3
        if (nk .le. 0) go to 103
    case (2, -2)
        nk = nk1
        if (nk .le. 0) go to 104
    case (3, -3)
        nk = (nk1 - 1)*nk2 + 1
        if (nk .le. 0) go to 105
    case (4, -4)
        nk = nk1
        if (nk .le. 0) go to 106
    case (5, -5)
        nk = (nk1*nk2*nk3)**2
        if (.not. double_kpgrid) double_kpgrid = .true.
        if (nk .le. 0) go to 107
    case default
        go to 108
    end select

    allocate (xkp(3, nk), wkp(nk))
    allocate (xkp_down(3, nk), wkp_down(nk))
    !
    ! initialization main quantities related to k-points
    !
    iflagerr = 0
    call get_kpoints(iflagerr, nel, nion, rion, atom_number, rs)

    if (iflagerr .ne. 0) go to 102

    if (.not. double_kpgrid) then
        xkp_down(:, :) = xkp(:, :)
        wkp_down(:) = wkp(:)
    end if

    tot_wt = sum(wkp(:))
    if (abs(tot_wt - 1.d0) .gt. 1d-6) then
        wkp(:) = 1.d0/nk
        tot_wt = 1.d0
    end if

    tot_wt_down = sum(wkp_down(:))
    if (abs(tot_wt - 1.d0) .gt. 1d-6) then
        wkp_down(:) = 1.d0/nk
        tot_wt_down = 1.d0
    end if

    write (6, *)
    write (6, *) ' type k-points/# k-points/total weight ', kp_type, nk, tot_wt
    write (6, *) ' k-points up/weights:'
    do i = 1, nk
        write (6, 90) i, xkp(1, i), xkp(2, i), xkp(3, i), wkp(i)
    end do
    write (6, *)
    write (6, *) ' k-points down/weights:'
    do i = 1, nk
        write (6, 90) i, xkp_down(1, i), xkp_down(2, i), xkp_down(3, i), wkp_down(i)
    end do
90  format(3x, I6, X, 4f11.7)

    if (allocated(xkp)) deallocate (xkp)
    if (allocated(xkp_down)) deallocate (xkp_down)
    if (allocated(wkp)) deallocate (wkp)
    if (allocated(wkp_down)) deallocate (wkp_down)
    call deallocate_all
    close (ufort10)

    stop

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! ERRORS !!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

101 write (6, *) 'ERROR: input fort.10 not found!'
    stop
102 write (6, *) 'ERROR: reading/generating k-points'
    stop
103 write (6, *) 'ERROR: check k-points grid !!'
    stop
104 write (6, *) 'ERROR: specify nk1 greater than 0 if using kp_type=2 !!'
    stop
105 write (6, *) 'ERROR: specify nk1 and nk2 greater than 0 if using kp_type=3 !!'
    stop
106 write (6, *) 'ERROR: specify nk1 greater than 0 if using kp_type=4 !!'
    stop
107 write (6, *) 'ERROR: check k-points grid for kp_type=5 !!'
    stop
108 write (6, *) 'ERROR: kp_type not recognized !!'
    stop

end program find_kpoints
