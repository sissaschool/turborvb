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

program copyjas

    use allio
    implicit none
    real(8), dimension(:), allocatable :: vj_sav, vju_sav, jasmat_sav, jasmatsz_sav, atom_number_sav, dup_c_store
    real(8), dimension(:, :), allocatable :: rion_store
    real(8) rs_store, celldm_store(6), rdiff(3), rtry(3), distnew, distr, center(3), &
        center_store(3)
    integer, dimension(:), allocatable :: multj_sav, nparamj_sav, &
                                          ioptorbj_sav, kionj_sav, kionj_new, ioccj_sav, nozeroj_sav, jbraj_sav, &
                                          jbrajn_sav, jbraiesm_savsav, map_c, mapj
    integer :: j, occj_sav, nnozeroj_sav, contractionj_sav, nshellj_sav, &
               dimjbraj, dimjbrajn, dimjbraiesm, dimvj, dimvju, iesmind_sav, iesfreer_sav, &
               iesdrr_sav, npar3body_sav, nelorbj_sav, nion_sav, iesdr_sav, niesd_sav, &
               nshelljm, nion_new, i, ipj_sav, stat, dimdup, nmoltot, nmoldiff, nelorb_old, ipj_old
    logical :: iessz_sav, compatibility, yesghost, found, yeslbox_store, iessz_old
    logical, dimension(:), allocatable :: occupied

    integer, parameter :: start_index = 100
    integer :: nkpoints, ikp, nproc_in, using_kcomp, index_file
    logical :: do_kpoints, eof, copyrion, copydup, copy_jas, copy_out, copy_in
    integer nel_sav, case_map_sav
    logical yes_crystalj_sav, chosen_map_sav, nomul

    !   AAA    Lines to be added just after all definitions of variables.
    character(100) :: name_tool
    character(20) :: str
    real*8 rcn(3)

    call get_command_argument(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then
        ! Input the name of the file exactly as it is in /doc
        name_tool = 'copyjas'
        call help_online(name_tool)
        stop
    end if
    !    AAA   end lines to be added

    ! k-points option as input line command
    do_kpoints = .false.
    nkpoints = 0
    nproc_in = 1
    copyrion = .false.
    copydup = .false.
    copy_jas = .true.
    copy_out = .false.
    copy_in = .false.

    if (trim(str) .eq. "kpoints") then
        do_kpoints = .true.
        using_kcomp = 0
    elseif (trim(str) .eq. "kpointsR") then
        do_kpoints = .true.
        using_kcomp = 0
        copyrion = .true.
        write (6, *) ' Warning copying ion positions '
    elseif (trim(str) .eq. "copyR") then
        do_kpoints = .false.
        using_kcomp = 1
        copyrion = .true.
        copy_jas = .false.
        write (6, *) ' Warning copying ONLY ion positions '
    elseif (trim(str) .eq. "kpointsRd") then
        do_kpoints = .true.
        using_kcomp = 0
        copyrion = .true.
        copydup = .true.
        write (6, *) ' Warning copying also contracted and ion positions '
    elseif (trim(str) .eq. "kpointsK") then
        using_kcomp = 1
        write (6, *) ' Insert number of processors and number of k-points. '
        read (5, *) nproc_in, nkpoints
    elseif (trim(str) .eq. "copy_out") then
        copy_out = .true.
    elseif (trim(str) .eq. "copy_in") then
        copy_in = .true.
    end if
    ! read the wavefunction with the new Jastrow
    call default_allocate
    open (unit=10, file='fort.10_new', status='old', form='formatted', err=101)
    call read_fort10(10)

    close (10)
    ! E perche' ???
    ! if( do_kpoints .and. (.not. yes_complex) ) go to 102
    !   ndiff=nelup-neldo
    nmoldiff = nmol + ndiff
    nmoltot = nmoldiff
    if (.not. symmagp .or. ipc .eq. 2) nmoltot = nmoltot + nmol

    ipj_sav = ipj
    nion_sav = 0
    nel_sav = nel
    do i = 1, nion
        if (atom_number(i) .gt. 0) nion_sav = nion_sav + 1
    end do
    nshelljm = max(nshellj_c, 1)

    allocate (kionj_sav(nshelljm), nparamj_sav(nshelljm), &
              ioptorbj_sav(nshelljm), multj_sav(nshelljm), &
              ioccj_sav(max(occj_c, 1)), atom_number_sav(nion), rion_store(3, nion))

    rion_store(:, :) = rion(:, :)
    rs_store = rs
    yeslbox_store = yeslbox
    celldm_store(:) = celldm(:)

    iesdrr_sav = iesdrr
    iesdr_sav = iesdr
    niesd_sav = niesd
    npar3body_sav = npar3body_c
    atom_number_sav = atom_number

    dimvj = size(vj)
    dimdup = size(dup_c) - 2*ipc*nelorbh*nmoltot ! do not change the molecular orbitals

    allocate (vj_sav(dimvj))
    allocate (dup_c_store(dimdup))
    vj_sav = vj
    dimvju = size(vju_c)
    allocate (vju_sav(dimvju))
    vju_sav = vju_c
    dup_c_store(1:dimdup) = dup_c(1:dimdup)

    nshellj_sav = nshellj_c
    occj_sav = occj_c
    yes_crystalj_sav = yes_crystalj
    chosen_map_sav = chosen_map
    case_map_sav = case_map

    kionj_sav = kionj_c
    nparamj_sav = nparamj_c
    ioptorbj_sav = ioptorbj_c
    multj_sav = multj_c
    ioccj_sav = ioccj_c

    contractionj_sav = contractionj
    iessz_sav = iessz
    nnozeroj_sav = nnozeroj_c
    if (contractionj .gt. 0) then
        nelorbj_sav = nelorbj_c
        allocate (jasmat_sav(max(ipj*ipj*nelorbj_c*nelorbj_c, 1)))
        jasmat_sav = jasmat_c
        allocate (nozeroj_sav(max(nnozeroj_sav, 1)))
        nozeroj_sav = nozeroj_c
        if (iessz) then
            allocate (jasmatsz_sav(max(nelorbj_c*nelorbj_c, 1)))
            jasmatsz_sav = jasmatsz_c
        end if
    else
        nelorbj_sav = nelorbjh
        allocate (jasmat_sav(max(ipj*ipj*nelorbjh*nelorbjh, 1)))
        nnozeroj_sav = nnozeroj
        jasmat_sav = jasmat
        allocate (nozeroj_sav(max(nnozeroj_sav, 1)))
        nozeroj_sav = nozeroj
        if (iessz) then
            allocate (jasmatsz_sav(max(nelorbjh*nelorbjh, 1)))
            jasmatsz_sav = jasmatsz
        end if
    end if

    write (6, *) ' nelorbj_sav here =', nelorbj_sav

    iesfreer_sav = iesfreer
    iesmind_sav = iesmind

    dimjbraj = size(jbraj)
    allocate (jbraj_sav(dimjbraj))
    jbraj_sav = jbraj

    dimjbrajn = size(jbrajn)
    allocate (jbrajn_sav(dimjbrajn))
    jbrajn_sav = jbrajn

    dimjbraiesm = size(jbraiesm_sav)
    allocate (jbraiesm_savsav(dimjbraiesm))
    jbraiesm_savsav = jbraiesm_sav

    ! read the wavefunction
    if (do_kpoints) then
        ! search for fort.10_000xxx files
        call system('ls fort.10_* > temporary_list')
        open (unit=11, file='temporary_list', status='unknown', form='formatted')
        rewind (11)

        ! search the number of k-points automatically if
        ! not using K computer
        if (using_kcomp .eq. 0) then
            eof = .false.
            str = ""
            do while (.not. eof)
                read (11, *, end=104) str
                if (trim(str) .ne. "fort.10_new") nkpoints = nkpoints + 1
            end do
104         write (6, '(A,I5/)') " Number of k-points found: ", nkpoints
            nproc_in = nkpoints
        end if

        ! open all fort.10_000xxx files related to different k-points.
        ! Needed for writing afterwards.
        rewind (11)
        index_file = 0
        write (6, *) ' Wavefunctions found '
        open (unit=start_index, file='fort.10', status='unknown', form='formatted', err=103)
        write (6, '(A,I6,2A)') '  Wavefunction', 0, ' = ', 'fort.10'
        do i = 1, nproc_in
            read (11, *) str
            if (mod((i - 1), nproc_in/nkpoints) .eq. 0) then
                index_file = index_file + 1
                write (6, '(A,I6,2A)') '  Wavefunction', i, ' = ', trim(str)
                open (unit=start_index + index_file, file=trim(str), status='old', form='formatted')
            end if
        end do
        close (11, status='delete')

    else

        open (unit=start_index, file='fort.10', status='unknown', form='formatted', err=103)

    end if

    dimjbraiesm = size(jbraiesm_sav)

    !   do ikp = 0,nkpoints
    !
    !      write(6,*) start_index+ikp
    !      read(start_index+ikp,*) str
    !      read(start_index+ikp,*) str
    !
    !      write(6,*) str
    !
    !   enddo
    !
    !   stop

    do ikp = 0, nkpoints

        !      write(6,*) ' Copying Jastrow to wavefunction: ',start_index+ikp

        call deallocate_all
        call default_allocate
        call read_fort10(start_index + ikp)

        nelorb_old = nelorbj_c*ipj

        if (copy_in) then
            allocate (mapj(nelorb_old), kionj_new(nshellj_sav))
            mapj = 0
            kionj_new = kionj_sav
!           checking the atom order

            do j = 1, nshellj_sav
                found = .false.
                do i = 1, nshellj_c
                    rcn(:) = rion_store(:, kionj_new(j)) - rion(:, kionj_c(i))
                    if (iespbc) call ApplyPBC(rcn, 1)
                    if (((ioptorbj_sav(j) .eq. 200 .and. ioptorbj_c(i) .eq. 200)&
                   &.or. (sum(rcn(:)**2) .lt. 1d-8 .and. ioptorbj_sav(j) .ne. 200&
                   &.and. ioptorbj_sav(j) .eq. ioptorbj_c(i))) .and. .not. found) then
                        kionj_sav(j) = kionj_c(i)
                        found = .true.
                    end if
                end do
                if (.not. found) then
                    write (6, *) ' fort.10 and fort.10_new should have the same ion positions '
                    stop
                end if
!         write(6,*) ' reshuff kionj_sav =',j,kionj_sav(j)
            end do
            deallocate (kionj_new)

            if (contractionj .eq. 0 .and. contractionj_sav .eq. 0) then

                call mappingu(nshellj_sav, multj_sav, ioccj_sav, kionj_sav&
               &, ioptorbj_sav, nion_sav, vju_sav, nshellj_c, multj_c, ioccj_c&
               &, kionj_c, ioptorbj_c, nion, vju_c, mapj)

            else

                call mapping(nshellj_sav, multj_sav, ioccj_sav, kionj_sav &
                             , ioptorbj_sav, nion_sav, nshellj_c, multj_c &
                             , ioccj_c, kionj_c, ioptorbj_c, nion, mapj)

            end if
            if (ipj_sav .eq. 2) then
                if (ipj .eq. 2) then
                    do i = 1, nelorbj_c
                        if (mapj(i) .gt. 0) then
                            mapj(nelorbj_c + i) = mapj(i) + nelorbj_sav
                        elseif (mapj(i) .lt. 0) then
                            mapj(nelorbj_c + i) = mapj(i) - nelorbj_sav
                        end if
                    end do
                else
                    do i = 1, nelorbj_c
                        mapj(i + nelorbj_c) = mapj(i)
                    end do
                end if
            elseif (ipj .eq. 2) then
                write (6, *) ' ERROR new and old Jastrow not compatible with ipj=2 ! '
                stop
            end if
            write (6, *) ' mapping found '
            do i = 1, nelorbj_c*ipj
                write (6, *) i, mapj(i)
            end do
            nelorb_old = nelorbj_sav*ipj_sav
        end if
        if (copy_out) then
            allocate (mapj(ipj_sav*nelorbj_sav), kionj_new(nshellj_c))
            mapj = 0
            kionj_new = kionj_c
!           checking the atom order
            do j = 1, nshellj_c
                found = .false.
                do i = 1, nshellj_sav
                    rcn(:) = rion_store(:, kionj_new(j)) - rion(:, kionj_sav(i))
                    if (iespbc) call ApplyPBC(rcn, 1)
                    if (((ioptorbj_sav(j) .eq. 200 .and. ioptorbj_c(i) .eq. 200)&
                   &.or. (sum(rcn(:)**2) .lt. 1d-8 .and. ioptorbj_sav(j) .ne. 200&
                   &.and. ioptorbj_sav(j) .eq. ioptorbj_c(i))) .and. .not. found) then
                        kionj_c(j) = kionj_sav(i)
                        found = .true.
                    end if
                end do
                if (.not. found) then
                    write (6, *) ' fort.10 and fort.10_new should have the same ion positions '
                    stop
                end if
            end do
            deallocate (kionj_new)
            if (contractionj .eq. 0 .and. contractionj_sav .eq. 0) then
                call mappingu(nshellj_c, multj_c, ioccj_c, kionj_c, ioptorbj_c&
               &, nion, vju_c, nshellj_sav, multj_sav, ioccj_sav, kionj_sav&
               &, ioptorbj_sav, nion_sav, vju_sav, mapj)
            else
                call mapping(nshellj_c, multj_c, ioccj_c, kionj_c, ioptorbj_c &
                             , nion_sav, nshellj_sav, multj_sav, ioccj_sav, kionj_sav &
                             , ioptorbj_sav, nion, mapj)
            end if
            if (ipj_sav .eq. 2) then
                if (ipj .eq. 2) then
                    do i = 1, nelorbj_sav
                        if (mapj(i) .gt. 0) then
                            mapj(nelorbj_sav + i) = mapj(i) + nelorbj_c
                        elseif (mapj(i) .lt. 0) then
                            mapj(nelorbj_sav + i) = mapj(i) - nelorbj_c
                        end if
                    end do
                else
                    do i = 1, nelorbj_sav
                        mapj(i + nelorbj_sav) = mapj(i)
                    end do
                end if
            elseif (ipj .eq. 2) then
                write (6, *) ' ERROR new and old Jastrow not compatible with ipj=2 ! '
                stop
            end if
            write (6, *) ' mapping found '
            do i = 1, nelorbj_sav*ipj_sav
                write (6, *) i, mapj(i)
            end do
        end if

        nion_new = 0
        yesghost = .false.
        do i = 1, nion
            if (atom_number(i) .gt. 0) then
                nion_new = nion_new + 1
            else
                yesghost = .true.
            end if
        end do
        nomul = .false.
        if (nion_sav .eq. nion_new) nomul = .true. ! no copy on different subsystems

        if (nion_sav .ne. nion_new .and. .not. copy_in) go to 105

        if (copy_jas) then
            !replace the jastrow read in fort.10_new
            ipj_old = ipj

            if (.not. copy_in) then
                nshellj_c = nshellj_sav
                case_map = case_map_sav
                chosen_map = chosen_map_sav
                yes_crystalj = yes_crystalj_sav
                ipj = ipj_sav

                !       allocate(kionj_sav(nshellj_c),nparamj_sav(nshellj_c)&
                !     &,ioptorbj_sav(nshellj_c),multj_sav(nshellj_c),ioccj_sav(occj_c))

                if (copy_out) then
                    !     do not change one-two body Jastrow corresponding to fort.10 read
                    niesd_sav = niesd
                    dimvj = size(vj)
                    if (dimvj .gt. size(vj_sav)) then
                        deallocate (vj_sav)
                        allocate (vj_sav(dimvj))
                    end if
                    vj_sav(1:dimvj) = vj(1:dimvj)
                end if
                npar3body_c = npar3body_sav
            end if ! endif copy_in
            deallocate (vj)
            allocate (vj(dimvj))

            iesdrr = iesdrr_sav
            iesdr = iesdr_sav
            niesd = niesd_sav
            vj = vj_sav
            if (.not. copy_in) then
                deallocate (vju_c)
                allocate (vju_c(dimvju))
                vju_c = vju_sav
            end if
        end if ! endif copy_jas

        if (copyrion) then
            !      center_store(:)=0.d0
            !      center(:)=0.d0
            !      do j=1,nion
            !      center_store(:)=center_store(:)+rion_store(:,j)
            !      center(:)=center(:)+rion(:,j)
            !      enddo
            !      center(:)=center(:)/nion
            !      center_store(:)=center_store(:)/nion
            !
            center(:) = rion(:, 1)
            center_store(:) = rion_store(:, 1)

            do j = 1, nion
                rion(:, j) = rion(:, j) + center_store(:) - center(:)
            end do

            do j = 1, nion
                rdiff(:) = rion_store(:, 1) - rion(:, j)
                if (iespbc) call ApplyPBC(rdiff, 1)
                distr = sum(rdiff(:)**2)
                rtry(:) = rion_store(:, 1)
                do i = 2, nion
                    rdiff(:) = rion_store(:, i) - rion(:, j)
                    if (iespbc) call ApplyPBC(rdiff, 1)
                    distnew = sum(rdiff(:)**2)
                    if (distnew .lt. distr) then
                        rtry(:) = rion_store(:, i)
                        distr = distnew
                    end if
                end do
                write (6, *) ' min dist =', j, sqrt(distr)
                rion(:, j) = rtry(:)
            end do
            celldm = celldm_store
            rs = rs_store
            yeslbox = yeslbox_store
        end if
        if (copydup) dup_c(1:dimdup) = dup_c_store(1:dimdup)

        if (copy_jas) then
            if (.not. copy_in) then
                nshellj_c = nshellj_sav
                case_map = case_map_sav
                chosen_map = chosen_map_sav
                yes_crystalj = yes_crystalj_sav
                occj_c = occj_sav
                deallocate (kionj_c, nparamj_c, ioptorbj_c, multj_c, ioccj_c)

                allocate (kionj_c(nshelljm), nparamj_c(nshelljm), ioptorbj_c(nshelljm), &
                          multj_c(nshelljm), ioccj_c(max(occj_c, 1)), occupied(nion), map_c(nion_new))

                if (yesghost) then
                    occupied = .true.
                    do i = 1, nion_new
                        j = 1
                        found = .true.
                        do while (j .le. nion .and. found)
                            if (atom_number_sav(i) .eq. -atom_number(j) .and. occupied(j)) then
                                found = .false.
                                occupied(j) = .false.
                                map_c(i) = j
                            end if
                            j = j + 1
                        end do
                        if (found) then
                            j = 1
                            found = .true.
                            do while (j .le. nion .and. found)
                                if (atom_number_sav(i) .eq. atom_number(j) .and. occupied(j)) then
                                    found = .false.
                                    occupied(j) = .false.
                                    map_c(i) = j
                                end if
                                j = j + 1
                            end do
                        else
                            do j = 1, nion
                                if (atom_number(j) .eq. atom_number_sav(i)) occupied(j) = .false.
                            end do
                        end if
                        if (found) go to 106

                    end do

                    do i = 1, nshelljm
                        kionj_c(i) = map_c(kionj_sav(i))
                    end do
                else
                    kionj_c = kionj_sav
                end if
                multj_c = multj_sav
                ioccj_c = ioccj_sav
                nparamj_c = nparamj_sav
                ioptorbj_c = ioptorbj_sav

                contractionj = contractionj_sav
                iessz_old = iessz
                iessz = iessz_sav

                nnozeroj_c = nnozeroj_sav
            end if

            if (contractionj .gt. 0) then
                if (copy_out) then
                    deallocate (jasmat_sav)
                    allocate (jasmat_sav(max(nelorb_old*nelorb_old, 1)))
                    jasmat_sav = jasmat_c
                end if
                if (.not. copy_in) then
                    nelorbj_c = nelorbj_sav
                    deallocate (jasmat_c)
                    allocate (jasmat_c(max(ipj*ipj*nelorbj_c*nelorbj_c, 1)))
                end if
                if (copy_out .or. copy_in) then
                    if (ipj_old .eq. 1 .and. ipj .eq. 2) then
                        call copyin2to1(jasmat_sav, nelorb_old, jasmat_c, nelorbj_c*ipj, mapj, nelup, neldo, nomul)
                    else
                        call copyin(jasmat_sav, nelorb_old, nel_sav, jasmat_c, nelorbj_c*ipj, nel, mapj, nomul)
                    end if
                    if (ipj .eq. 2 .and. iessz_old) then
                        call copyinsz(jasmatsz, nelorb_old, jasmat_c, nelorbj_c*ipj, mapj, nomul)
                    end if
                else
                    jasmat_c = jasmat_sav
                end if
                if (.not. copy_in) then
                    deallocate (nozeroj_c)
                    allocate (nozeroj_c(max(nnozeroj_c, 1)))
                    nozeroj_c = nozeroj_sav
                end if
                if (iessz) then
                    if (allocated(jasmatsz_c) .and. .not. copy_in) then
                        if (copy_out) then
                            deallocate (jasmat_sav)
                            allocate (jasmat_sav(max(nelorb_old*nelorb_old, 1)))
                            jasmat_sav = jasmatsz_c
                        end if
                        deallocate (jasmatsz_c)
                    end if
                    if (.not. copy_in) allocate (jasmatsz_c(max(nelorbj_c*nelorbj_c, 1)))
                    if (copy_out .or. copy_in) then
                        call copyin(jasmat_sav, nelorb_old, nel_sav, jasmatsz_c, nelorbj_c, nel, mapj, nomul)
                    else
                        jasmatsz_c = jasmatsz_sav
                    end if
                end if

            else
                if (copy_out) then
                    deallocate (jasmat_sav)
                    allocate (jasmat_sav(max(nelorb_old*nelorb_old, 1)))
                    jasmat_sav = jasmat
                end if
                if (.not. copy_in) then
                    nelorbjh = nelorbj_sav
                    nelorbj_c = nelorbj_sav
                    deallocate (jasmat)
                    allocate (jasmat(max(ipj*ipj*nelorbjh*nelorbjh, 1)))
                    nnozeroj = nnozeroj_sav
                end if
                if (copy_out .or. copy_in) then
                    if (ipj_old .eq. 1 .and. ipj .eq. 2) then
                        call copyin2to1(jasmat_sav, nelorb_old, jasmat, nelorbjh*ipj, mapj, nelup, neldo, nomul)
                    else
!  passes here in the normal case
                        call copyin(jasmat_sav, nelorb_old, nel_sav, jasmat, nelorbjh*ipj, nel, mapj, nomul)
                    end if
                    if (ipj .eq. 2 .and. iessz_old) then
                        call copyinsz(jasmatsz, nelorb_old, jasmat, nelorbjh*ipj, mapj, nomul)
                    end if
                else
                    jasmat = jasmat_sav
                end if

                if (.not. copy_in) then
                    deallocate (nozeroj, nozeroj_c)
                    allocate (nozeroj(max(nnozeroj, 1)), nozeroj_c(max(nnozeroj, 1)))
                    nozeroj = nozeroj_sav
                    nozeroj_c = nozeroj_sav
                end if
                if (iessz) then
                    if (allocated(jasmatsz) .and. .not. copy_in) then
                        deallocate (jasmat_sav)
                        if (copy_out) then
                            allocate (jasmat_sav(max(nelorb_old*nelorb_old, 1)))
                            jasmat_sav = jasmatsz
                        end if
                        deallocate (jasmatsz)
                    end if
                    if (.not. copy_in) allocate (jasmatsz(max(nelorbjh*nelorbjh, 1)))
                    if (copy_out .or. copy_in) then
                        call copyin(jasmat_sav, nelorb_old, nel_sav, jasmatsz, nelorbjh, nel, mapj, nomul)
                    else
                        jasmatsz = jasmatsz_sav
                    end if
                end if
            end if

            if (.not. copy_in) then
                iesfreer = iesfreer_sav
                iesmind = iesmind_sav

                deallocate (jbraj)
                allocate (jbraj(dimjbraj))
                jbraj = jbraj_sav

                deallocate (jbrajn)
                allocate (jbrajn(dimjbrajn))
                jbrajn = jbrajn_sav

                deallocate (jbraiesm_sav)
                allocate (jbraiesm_sav(dimjbraiesm))

                jbraiesm_sav = jbraiesm_savsav
            end if

            if (.not. copy_in) deallocate (occupied, map_c)
            if (copy_out .or. copy_in) deallocate (mapj)

        end if ! endif copy_jas

        call write_fort10(start_index + ikp)

        close (start_index + ikp)

    end do
    if (allocated(kionj_sav)) deallocate (kionj_sav, nparamj_sav&
    &, ioptorbj_sav, multj_sav, ioccj_sav, atom_number_sav, rion_store)
    if (allocated(mapj)) deallocate (mapj)
    stop

101 write (6, *) ' ERROR: wavefunction fort.10_new with the new Jastrow not found or wrong! '
    stop
102 write (6, *) ' ERROR: the wavefunction must be complex in the case of k-points!'
    stop
103 write (6, *) ' ERROR: wavefunction fort.10 not found or wrong! '
    stop
105 write (6, *) ' ERROR: fort.10_new has a different number of ion, fort.10 unchanged!'
    stop
106 write (6, *) ' ERROR: inconsistency found in the definition of ghost atoms!'
    stop

end program copyjas
subroutine mappingu(nshell_c, mult_c, ioccup_c, kion_c, ioptorb_c, nion_c &
                    , vju_c, nshell, mult, ioccup, kion, ioptorb, nion, vju &
                    , map)
    implicit none
    integer nshell, nshell_c, i, j, k, l, indorb, indorbnew, ind, indnew&
            &, nion_c, nion
    integer mult_c(*), ioccup_c(*), kion_c(*), mult(*), ioccup(*), kion(*), map(*), ioptorb_c(*), ioptorb(*)
    real*8 vju_c(*), vju(*)
    integer indp, indp_c, indpu, indpu_c
    indorbnew = 0
    indnew = 0
    do k = 1, nshell
        do l = 1, mult(k)
            indnew = indnew + 1
            if (ioccup(indnew) .eq. 1) then
                indorbnew = indorbnew + 1
                map(indorbnew) = 0
            end if
        end do
    end do
    indorb = 0
    ind = 0
    indp_c = 0
    do i = 1, nshell_c
        if (ioptorb_c(i) .ne. 200) indp_c = indp_c + 1
        indpu_c = max(indp_c, 1) ! to avoid the first orbital is 200
        do j = 1, mult_c(i)
            ind = ind + 1
            if (ioccup_c(ind) .eq. 1) then
                indorb = indorb + 1
                indorbnew = 0
                indnew = 0
                indp = 0
                do k = 1, nshell
                    if (ioptorb(k) .ne. 200) indp = indp + 1
                    indpu = max(indp, 1) ! to avoid the first orbital is 200
                    do l = 1, mult(k)
                        indnew = indnew + 1
                        if (ioccup(indnew) .eq. 1) then
                            indorbnew = indorbnew + 1
                            if ((l .eq. j .and. kion(k) .eq. kion_c(i) .and.&
                           &(ioptorb(k) .eq. ioptorb_c(i)) .and. abs(vju(indpu) - vju_c(indpu_c)) .le. 1d-4)&
                           &.or. (ioptorb(k) .eq. 200 .and. ioptorb_c(i) .eq. 200)) then
                                map(indorbnew) = indorb
                                if (ioptorb(k) .eq. 200) map(indorbnew) = -map(indorbnew)
!                                                write(6,*) ' Match found =',kion(k),kion_c(i),indorb,indorbnew,adr,adrnew
                            end if
                        end if
                    end do
                end do
            end if
        end do
    end do
    return
end
subroutine mapping(nshell_c, mult_c, ioccup_c, kion_c, ioptorb_c, nion_c, nshell, mult, ioccup, kion, ioptorb, nion, map)
    implicit none
    integer nshell, nshell_c, i, j, k, l, indorb, indorbnew, ind, indnew, adr, adrnew&
            &, nion_c, nion
    integer mult_c(*), ioccup_c(*), kion_c(*), mult(*), ioccup(*), kion(*), map(*), ioptorb_c(*), ioptorb(*)
    indorbnew = 0
    indnew = 0
    do k = 1, nshell
        if (k .eq. 1) then
            adr = 0
        else
            if (kion(k - 1) .ne. kion(k) .or. ioptorb(k - 1) .ne. ioptorb(k)) adr = 0
        end if
        do l = 1, mult(k)
            indnew = indnew + 1
            if (ioccup(indnew) .eq. 1) then
                adr = adr + 1
                indorbnew = indorbnew + 1
                map(indorbnew) = 0
            end if
        end do
    end do
    write (6, *) ' New basis inside =', indorbnew
    indorb = 0
    ind = 0
    do i = 1, nshell_c
        if (i .eq. 1) then
            adr = 0
        else
!    if i-term points to different atom or orbital compared to the previous one
            if (kion_c(i - 1) .ne. kion_c(i) .or. ioptorb_c(i - 1) .ne. ioptorb_c(i)) adr = 0
        end if
        do j = 1, mult_c(i)
            ind = ind + 1
            if (ioccup_c(ind) .eq. 1) then
                indorb = indorb + 1
                adr = adr + 1
                indorbnew = 0
                indnew = 0
                do k = 1, nshell
                    if (k .eq. 1) then
                        adrnew = 0
                    else
                        if (kion(k) .ne. kion(k - 1) .or. ioptorb(k) .ne. ioptorb(k - 1)) adrnew = 0
                    end if
                    do l = 1, mult(k)
                        indnew = indnew + 1
                        if (ioccup(indnew) .eq. 1) then
                            adrnew = adrnew + 1
                            indorbnew = indorbnew + 1
                            if ((kion(k) .eq. kion_c(i) .and. (adrnew .eq. adr .and. ioptorb(k) .eq. ioptorb_c(i))) &
                                .or. (ioptorb(k) .eq. 200 .and. ioptorb_c(i) .eq. 200)) then
                                map(indorbnew) = indorb
                                if (ioptorb(k) .eq. 200) map(indorbnew) = -map(indorbnew)
!                                                write(6,*) ' Match found =',kion(k),kion_c(i),indorb,indorbnew,adr,adrnew
                            end if
                        end if
                    end do
                end do
            end if
        end do
    end do
    return
end
subroutine copyin(jastrow_old, nelorbj_old, nel_old &
                  , jastrow_new, nelorbj_new, nel_new &
                  , mapj, nomul)
    implicit none
    integer nelorbj_old, nelorbj_new, i, j, nel_new, nel_old, nelorbu
    real*8 jastrow_old(nelorbj_old, *), jastrow_new(nelorbj_new, nelorbj_new)
    real*8 cost
    integer mapj(*)
    logical nomul, yescost
    nelorbu = nelorbj_old
    yescost = .false.
    do i = 1, nelorbj_new
        if (mapj(i) .lt. 0) yescost = .true.
    end do
    if (yescost) nelorbu = nelorbu - 1
    cost = dble(nel_old - 1)/dble(nel_new - 1)
    jastrow_new = 0.d0
    do i = 1, nelorbj_new
        do j = 1, nelorbj_new
            if (mapj(i) .ne. 0 .and. mapj(j) .ne. 0) then
                if (mapj(i) .lt. 0 .or. mapj(j) .lt. 0) then
                    jastrow_new(i, j) = jastrow_new(i, j) + jastrow_old(abs(mapj(i)), abs(mapj(j)))*cost
                else
!                if i and j refer to the same subsystem
                    if ((i - 1)/nelorbu .eq. (j - 1)/nelorbu .or. nomul) then
                        jastrow_new(i, j) = jastrow_old(abs(mapj(i)), abs(mapj(j)))
                    else
                        jastrow_new(i, j) = 0.d0
                    end if
                end if
            else
                jastrow_new(i, j) = 0.d0
            end if
        end do
    end do
    return
end
subroutine copyin2to1(jastrow_old, nelorbj_old, jastrow_new, nelorbj_new, mapj, nelup, neldo, nomul)
    implicit none
    integer nelorbj_old, nelorbj_new, i, j, nelup, neldo, nelorbh, nelorbu
    real*8 jastrow_old(nelorbj_old, *), jastrow_new(nelorbj_new, *)
    integer mapj(*)
    real*8 cost(2, 2)
    logical nomul, yescost
    nelorbu = nelorbj_old
    yescost = .false.
    do i = 1, nelorbj_new
        if (mapj(i) .lt. 0) yescost = .true.
    end do
    if (yescost) nelorbu = nelorbu - 1

    cost = 0.d0
    cost(1, 2) = dble(nelup + neldo - 1)/dble(nelup)
    if (nelup .gt. 1) cost(1, 1) = dble(nelup + neldo - 1)/dble(nelup - 1)
    if (neldo .gt. 1) cost(2, 2) = dble(nelup + neldo - 1)/dble(neldo - 1)
    if (neldo .gt. 0) cost(2, 1) = dble(nelup + neldo - 1)/dble(neldo)
    nelorbh = nelorbj_new/2
    do i = 1, nelorbj_new
        do j = 1, nelorbj_new
!                if i and j refer to the same subsystem
            if (mapj(i) .ne. 0 .and. mapj(j) .ne. 0) then
                if (mapj(i) .lt. 0 .or. mapj(j) .lt. 0) then
                    jastrow_new(i, j) = jastrow_old(abs(mapj(i)), abs(mapj(j)))
                    if (mapj(i) .lt. 0) then
                        jastrow_new(i, j) = jastrow_new(i, j)*cost((i - 1)/nelorbh + 1, (j - 1)/nelorbh + 1)
                    elseif (mapj(j) .lt. 0) then
                        jastrow_new(i, j) = jastrow_new(i, j)*cost((j - 1)/nelorbh + 1, (i - 1)/nelorbh + 1)
                    end if
                else
!                if i and j refer to the same subsystem
                    if ((i - 1)/nelorbu .eq. (j - 1)/nelorbu .or. nomul) then
                        jastrow_new(i, j) = jastrow_old(abs(mapj(i)), abs(mapj(j)))
                    else
                        jastrow_new(i, j) = 0.d0
                    end if
                end if
            else
                jastrow_new(i, j) = 0.d0
            end if
        end do
    end do
    return
end

subroutine copyinsz(jastrow_old, nelorbj_old, jastrow_new, nelorbj_new, mapj, nomul)
    implicit none
    integer nelorbj_old, nelorbj_new, i, j, nelorbjh, nelorbu
    real*8 jastrow_old(nelorbj_old, *), jastrow_new(nelorbj_new, *)
    integer mapj(*)
    logical nomul, yescost
    nelorbu = nelorbj_old
    yescost = .false.
    do i = 1, nelorbj_new
        if (mapj(i) .lt. 0) yescost = .true.
    end do
    if (yescost) nelorbu = nelorbu - 1

    nelorbjh = nelorbj_new/2
    do i = 1, nelorbj_new
        do j = 1, nelorbj_new
            if (mapj(i) .ne. 0 .and. mapj(j) .ne. 0 .and. ((i - 1)/nelorbu .eq. (j - 1)/nelorbu .or. nomul)) then
                if (i .gt. nelorbjh .and. j .le. nelorbjh .or.&
                        &   i .le. nelorbjh .and. j .gt. nelorbjh) then
                    jastrow_new(i, j) = jastrow_new(i, j) - jastrow_old(abs(mapj(i)), abs(mapj(j)))
                else
                    jastrow_new(i, j) = jastrow_new(i, j) + jastrow_old(abs(mapj(i)), abs(mapj(j)))
                end if
                !            else
                !             jastrow_new(i,j)=0.d0
            end if
        end do
    end do
    return
end
