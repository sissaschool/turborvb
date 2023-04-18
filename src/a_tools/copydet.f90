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

program copydet

    use allio
    implicit none
    real(8), dimension(:), allocatable :: vj_sav, vju_sav, jasmat_sav, jasmatsz_sav, atom_number_sav, dup_c_store
    real(8), dimension(:, :), allocatable :: rion_store
    real(8) rs_store, celldm_store(6), rdiff(3), rtry(3), distnew, distr, center(3), &
        center_store(3), cellscale_sav(3)
    integer, dimension(:), allocatable :: multj_sav, nparamj_sav, &
                                          ioptorbj_sav, kionj_sav, ioccj_sav, nozeroj_sav, jbraj_sav, &
                                          jbrajn_sav, jbraiesm_savsav, mapj, mapion, kiontot_sav, atbasis
    integer :: j, occj_sav, nnozeroj_sav, contractionj_sav, nshellj_sav, &
            dimjbraj, dimjbrajn, dimjbraiesm, dimvj, dimvju, iesmind_sav, iesfreer_sav, &
            iesdrr_sav, npar3body_sav, nelorbj_sav, nion_sav, iesdr_sav, niesd_sav, &
            nshelljm, nion_new, i, ipf_sav, stat, dimdup, nmoltot, nmoldiff&
            &, nelorb_old, ipf_old, imin, k, imap(3), indpar, jj, countat_old, countat_new&
            &, souncazzoio, ix_old, iy_old, ind, ion_ix, ion_iy
    logical :: iessz_sav, compatibility, yesghost, found, yeslbox_store, iessz_old
    logical, dimension(:), allocatable :: occupied

    integer, parameter :: start_index = 100
    integer :: nkpoints, ikp, nproc_in, using_kcomp, index_file
    logical :: do_kpoints, eof, copyrion, copydup, copy_jas, copy_out

    !   AAA    Lines to be added just after all definitions of variables.
    character(100) :: name_tool
    character(20) :: str

    call get_command_argument(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then
        ! Input the name of the file exactly as it is in /doc
        name_tool = 'copydet'
        call help_online(name_tool)
        stop
    end if
    ! read the wavefunction with the new Jastrow
    call default_allocate
    open (unit=10, file='fort.10_new', status='old', form='formatted', err=101)
    call read_fort10(10)
    close (10)
    nmoldiff = nmol + ndiff
    nmoltot = nmoldiff
    if (.not. symmagp .or. ipc .eq. 2) nmoltot = nmoltot + nmol

    ipf_sav = ipf
    nion_sav = 0
    do i = 1, nion
        if (atom_number(i) .gt. 0) nion_sav = nion_sav + 1
    end do
    nshelljm = max(nshell_c, 1)

    allocate (kionj_sav(nshelljm), nparamj_sav(nshelljm), &
              ioptorbj_sav(nshelljm), multj_sav(nshelljm), &
              ioccj_sav(size(ioccup_c)), atom_number_sav(nion), rion_store(3, nion))

    rion_store(:, :) = rion(:, :)
    rs_store = rs
    yeslbox_store = yeslbox
    celldm_store(:) = celldm(:)
    cellscale_sav(1:3) = cellscale(1:3)

    npar3body_sav = iesupr_c
    atom_number_sav = atom_number

    dimdup = size(dup_c) - 2*ipc*nelorbh*nmoltot ! do not change the molecular orbitals

    allocate (kiontot_sav(size(kiontot)))
    kiontot_sav = kiontot

    allocate (dup_c_store(dimdup))
    dup_c_store(1:dimdup) = dup_c(1:dimdup)

    nshellj_sav = nshell_c
    occj_sav = occ_c

    kionj_sav = kion_c
    nparamj_sav = nparam_c
    ioptorbj_sav = ioptorb_c
    multj_sav = mult_c
    ioccj_sav = ioccup_c

    contractionj_sav = contraction
    nnozeroj_sav = nnozero_c
    if (contraction .gt. 0) then
        nelorbj_sav = nelorb_c
        allocate (jasmat_sav(max(ipc*nelorb_c*nelcol_c, 1)))
        write (6, *) ' size match ? =', size(detmat_c), size(jasmat_sav)
        jasmat_sav = detmat_c
        allocate (nozeroj_sav(max(nnozeroj_sav, 1)))
        nozeroj_sav = nozero_c
    else
        nelorbj_sav = nelorbh
        allocate (jasmat_sav(max(ipc*ipf*nelorbh*nelcolh, 1)))
        nnozeroj_sav = nnozero
        jasmat_sav = detmat
        allocate (nozeroj_sav(max(nnozeroj_sav, 1)))
        nozeroj_sav = nozero
    end if

    iesfreer_sav = iesswr
    iesmind_sav = iesupind
    open (unit=10, file='fort.10', status='unknown', form='formatted', err=103)
    !    write(6,*) ' Copying Jastrow to wavefunction: ',start_index+ikp
    call deallocate_all
    call default_allocate
    call read_fort10(10)
    nelorb_old = nelorb_c
    allocate (mapion(nion))

    write (6, *) ' cellscale new/old =', (cellscale(i), cellscale_sav(i), i=1, 3)
    do i = 1, nion
        ! Find the closest atom to the old
        do j = 1, nion_sav
            rdiff(:) = rion(:, i) - rion_store(:, j)
            rdiff(1) = rdiff(1) - cellscale_sav(1)*anint(rdiff(1)/cellscale_sav(1))
            rdiff(2) = rdiff(2) - cellscale_sav(2)*anint(rdiff(2)/cellscale_sav(2))
            rdiff(3) = rdiff(3) - cellscale_sav(3)*anint(rdiff(3)/cellscale_sav(3))
            cost = sqrt(sum(rdiff(1:3)**2))
            if (j .eq. 1 .or. cost .lt. distnew) then
                rdiff(:) = rion(:, i) - rion_store(:, j)
                imap(:) = anint(rdiff(:)/cellscale_sav(:))
                imin = j
                distnew = cost
            end if
        end do
        mapion(i) = imin
        rion(:, i) = rion_store(:, imin) + imap(:)*cellscale_sav(:)
        write (6, *) ' mapion distmin', i, mapion(i), distnew
    end do
    do i = 1, nion_sav
        k = 0
        do j = 1, nion
            if (mapion(j) .eq. i) k = k + 1
        end do
        write (6, *) ' Atom i/mult =', i, k
    end do
    write (6, *) ' Shell read ='
    !    allocate(atbasis(nion_sav))
    !    atbasis=0
    !    do j=1,nion_sav
    !     do i=1,nshellj_sav
    !     if(kionj_sav(i).eq.j.and.ioptorbj_sav(i).lt.900000) atbasis(j)=atbasis(j)+multj_sav(i)
    !     enddo
    !    enddo
    indpar = 0
    !    dup_c=0.d0
    do i = 1, nshellj_sav
        write (6, *) i, kionj_sav(i), (dup_c_store(indpar + j), j=1, ipc*nparamj_sav(i))
        if (i .eq. 1) then
            countat_old = 0
        elseif (kionj_sav(i) .ne. kionj_sav(i - 1)) then
            countat_old = 0
        end if
        k = 0
        do j = 1, nshell_c
            if (j .eq. 1) then
                countat_new = 0
            elseif (kion_c(j) .ne. kion_c(j - 1)) then
                countat_new = 0
            end if
            if (mapion(kion_c(j)) .eq. kionj_sav(i)) then
                if (countat_new .eq. countat_old .and. ioptorbj_sav(i) .eq. ioptorb_c(j)) then
                    !        write(6,*) ' Found match ',j,k,kionj_sav(i),kion_c(j)
                    !
                    !         if(ioptorb_c(j).eq.16) write(6,*) ' Matching s ',j
                    !         if(ioptorb_c(j).eq.36) write(6,*) ' Matching p ',j
                    !         if(ioptorb_c(j).eq.68) write(6,*) ' Matching d ',j
                    if (ioptorbj_sav(i) .eq. 900000) then
                        !         Leave it unchanged
                        !         souncazzoio=atbasis(i)*(kion_c(j)-kionj_sav(i))
                        !         do jj=1,nparamj_sav(i)/2
                        !         dup_c(k+ipc*(jj-1)+1)=dup_c_store(indpar+ipc*(jj-1)+1)+souncazzoio
                        !         enddo
                        do jj = nparamj_sav(i)/2 + 1, nparamj_sav(i)
                            if (ipc .eq. 1) then
                                dup_c(k + jj) = dup_c_store(indpar + jj)
                            else
                                dup_c(k + ipc*(jj - 1) + 1:k + ipc*jj) = dup_c_store(indpar + ipc*(jj - 1) + 1:indpar + ipc*jj)
                            end if
                        end do
                    else
                        do jj = 1, ipc*nparamj_sav(i)
                            dup_c(k + jj) = dup_c_store(indpar + jj)
                        end do
                    end if
                end if
                countat_new = countat_new + ipc*nparamj_sav(i)
            end if
            k = k + ipc*nparam_c(j)
        end do
        indpar = indpar + ipc*nparamj_sav(i)
        countat_old = countat_old + ipc*nparamj_sav(i)
    end do

    write (6, *) ' Final new dup_c ', k, iesupr, iesup_c, size(dup_c)
    do i = 1, k
        write (6, *) i, dup_c(i)
    end do

    allocate (atbasis(nion_sav))
    atbasis = 0
    do j = 1, nion_sav
        do i = 1, nshellj_sav
            if (kionj_sav(i) .eq. j .and. ioptorbj_sav(i) .eq. 900000) atbasis(j) = atbasis(j) + multj_sav(i)
        end do
        if (atbasis(j) .ne. 0) write (6, *) ' atbasis hybrid =', j, atbasis(j)
    end do
    if (sum(atbasis(:)) .eq. 0) then
        atbasis = 0
        do j = 1, nion_sav
            do i = 1, nshellj_sav
                if (kionj_sav(i) .eq. j) atbasis(j) = atbasis(j) + multj_sav(i)
            end do
            write (6, *) ' atbasis standard =', j, atbasis(j)
        end do
    end if

    !     write(6,*) ' Input detmat ',nnozeroj_sav,nelorbj_sav
    !       do j=1,nnozeroj_sav
    !       ind=nozeroj_sav(j)
    !       iy_old=(ind-1)/nelorbj_sav+1
    !       ix_old=ind-(iy_old-1)*nelorbj_sav
    !       write(6,*) ix_old,iy_old,kiontot_sav(ix_old),kiontot_sav(iy_old),&
    !                  jasmat_sav(ipc*(nozeroj_sav(j)-1)+1:ipc*nozeroj_sav(j))
    !       enddo

    do i = 1, nnozero_c
        if (mod(i, 1000) .eq. 0) write (6, *) 'done =', i, 'over', nnozero_c
        ind = nozero_c(i)
        iy = (ind - 1)/nelorb_c + 1
        ix = ind - (iy - 1)*nelorb_c
        rdiff(:) = rion(:, kiontot(ix)) - rion(:, kiontot(iy))
        rdiff(1) = rdiff(1) - cellscale(1)*anint(rdiff(1)/cellscale(1))
        rdiff(2) = rdiff(2) - cellscale(2)*anint(rdiff(2)/cellscale(2))
        rdiff(3) = rdiff(3) - cellscale(3)*anint(rdiff(3)/cellscale(3))
        do j = 1, nnozeroj_sav
            ind = nozeroj_sav(j)
            iy_old = (ind - 1)/nelorbj_sav + 1
            ix_old = ind - (iy_old - 1)*nelorbj_sav
            ion_ix = kiontot_sav(ix_old)
            ion_iy = kiontot_sav(iy_old)
            rtry(:) = rion_store(:, ion_ix) - rion_store(:, ion_iy)
            rtry(1) = rtry(1) - cellscale_sav(1)*anint(rtry(1)/cellscale_sav(1))
            rtry(2) = rtry(2) - cellscale_sav(2)*anint(rtry(2)/cellscale_sav(2))
            rtry(3) = rtry(3) - cellscale_sav(3)*anint(rtry(3)/cellscale_sav(3))

            if (abs(rtry(1) - rdiff(1)) .lt. deps .and. abs(rtry(2) - rdiff(2)) .lt. deps .and.&
                    &abs(rtry(3) - rdiff(3)) .lt. deps .and. mod(ix, atbasis(ion_ix)) .eq. mod(ix_old, atbasis(ion_ix))&
                    &.and. mod(iy, atbasis(ion_iy)) .eq. mod(iy_old, atbasis(ion_iy))) then
                if (contraction .ne. 0) then
                    detmat_c(ipc*(nozero_c(i) - 1) + 1:ipc*nozero_c(i)) = &
                        jasmat_sav(ipc*(nozeroj_sav(j) - 1) + 1:ipc*nozeroj_sav(j))
                else
                    detmat(ipc*(nozero_c(i) - 1) + 1:ipc*nozero_c(i)) = &
                        jasmat_sav(ipc*(nozeroj_sav(j) - 1) + 1:ipc*nozeroj_sav(j))
                end if
            elseif (symmagp .and. abs(rtry(1) + rdiff(1)) .lt. deps .and. abs(rtry(2) + rdiff(2)) .lt. deps .and.&
                    &abs(rtry(3) + rdiff(3)) .lt. deps .and. mod(ix, atbasis(ion_iy)) .eq. mod(iy_old, atbasis(ion_iy))&
                    &.and. mod(iy, atbasis(ion_ix)) .eq. mod(ix_old, atbasis(ion_ix))) then
                ind = nelorbj_sav*(ix_old - 1) + iy_old
                if (contraction .ne. 0) then
                    detmat_c(ipc*(nozero_c(i) - 1) + 1:ipc*nozero_c(i)) = &
                        jasmat_sav(ipc*(ind - 1) + 1:ipc*ind)
                else
                    detmat(ipc*(nozero_c(i) - 1) + 1:ipc*nozero_c(i)) = &
                        jasmat_sav(ipc*(ind - 1) + 1:ipc*ind)
                end if
            end if
        end do
    end do

    call write_fort10(10)

    close (10)

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

end program copydet
subroutine mapping(nshell_c, mult_c, ioccup_c, kion_c, ioptorb_c, nshell, mult, ioccup, kion, ioptorb, map)
    implicit none
    integer nshell, nshell_c, i, j, k, l, indorb, indorbnew, ind, indnew, adr, adrnew
    integer mult_c(*), ioccup_c(*), kion_c(*), mult(*), ioccup(*), kion(*), map(*), ioptorb_c(*), ioptorb(*)
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
    write (6, *) ' New basis inside =', indorbnew
    indorb = 0
    ind = 0
    do i = 1, nshell_c
        do j = 1, mult_c(i)
            ind = ind + 1
            if (ioccup_c(ind) .eq. 1) then
                indorb = indorb + 1
                if (i .eq. 1) then
                    adr = 0
                else
                    if (kion_c(i - 1) .ne. kion_c(i)) adr = 0
                end if
                adr = adr + 1
                indorbnew = 0
                indnew = 0
                do k = 1, nshell
                    do l = 1, mult(k)
                        indnew = indnew + 1
                        if (ioccup(indnew) .eq. 1) then
                            if (k .eq. 1) then
                                adrnew = 0
                            else
                                if (kion(k) .ne. kion(k - 1)) adrnew = 0
                            end if
                            adrnew = adrnew + 1
                            indorbnew = indorbnew + 1
                            if (kion(k) .eq. kion_c(i) &
                                .and. ((adrnew .eq. adr .and. ioptorb(k) .eq. ioptorb_c(i)) &
                                       .or. (ioptorb(k) .eq. 200 .and. ioptorb_c(i) .eq. 200))) then
                                map(indorbnew) = indorb
                                if (ioptorb(k) .eq. 200) map(indorbnew) = -map(indorbnew)
                                ! write(6,*) ' Match found =',kion(k),kion_c(i),indorb,indorbnew,adr,adrnew
                            end if
                        end if
                    end do
                end do
            end if
        end do
    end do
    return
end
subroutine copyin(jastrow_old, nelorbj_old, jastrow_new, nelorbj_new, mapj)
    implicit none
    integer nelorbj_old, nelorbj_new, i, j
    real*8 jastrow_old(nelorbj_old, *), jastrow_new(nelorbj_new, *)
    integer mapj(*)
    do i = 1, nelorbj_new
        do j = 1, nelorbj_new
            if (mapj(i) .ne. 0 .and. mapj(j) .ne. 0) then
                if (mapj(i) .lt. 0 .or. mapj(j) .lt. 0) then
                    jastrow_new(i, j) = jastrow_old(abs(mapj(i)), abs(mapj(j)))
                else
                    jastrow_new(i, j) = jastrow_old(abs(mapj(i)), abs(mapj(j)))
                end if
            else
                jastrow_new(i, j) = 0.d0
            end if
        end do
    end do
    return
end
subroutine copyin2to1(jastrow_old, nelorbj_old, jastrow_new, nelorbj_new, mapj, nelup, neldo)
    implicit none
    integer nelorbj_old, nelorbj_new, i, j, nelup, neldo, nelorbh
    real*8 jastrow_old(nelorbj_old, *), jastrow_new(nelorbj_new, *)
    integer mapj(*)
    real*8 cost(2, 2)
    cost = 0.d0
    cost(1, 2) = dble(nelup + neldo - 1)/dble(nelup)
    if (nelup .gt. 1) cost(1, 1) = dble(nelup + neldo - 1)/dble(nelup - 1)
    if (neldo .gt. 1) cost(2, 2) = dble(nelup + neldo - 1)/dble(neldo - 1)
    if (neldo .gt. 0) cost(2, 1) = dble(nelup + neldo - 1)/dble(neldo)
    nelorbh = nelorbj_new/2
    do i = 1, nelorbj_new
        do j = 1, nelorbj_new
            if (mapj(i) .ne. 0 .and. mapj(j) .ne. 0) then
                if (mapj(i) .lt. 0 .or. mapj(j) .lt. 0) then
                    jastrow_new(i, j) = jastrow_old(abs(mapj(i)), abs(mapj(j)))
                    if (mapj(i) .lt. 0) then
                        jastrow_new(i, j) = jastrow_new(i, j)*cost((i - 1)/nelorbh + 1, (j - 1)/nelorbh + 1)
                    elseif (mapj(j) .lt. 0) then
                        jastrow_new(i, j) = jastrow_new(i, j)*cost((j - 1)/nelorbh + 1, (i - 1)/nelorbh + 1)
                    end if
                else
                    jastrow_new(i, j) = jastrow_old(abs(mapj(i)), abs(mapj(j)))
                end if
            else
                jastrow_new(i, j) = 0.d0
            end if
        end do
    end do
    return
end

subroutine copyinsz(jastrow_old, nelorbj_old, jastrow_new, nelorbj_new, mapj)
    implicit none
    integer nelorbj_old, nelorbj_new, i, j, nelorbjh
    real*8 jastrow_old(nelorbj_old, *), jastrow_new(nelorbj_new, *)
    integer mapj(*)
    nelorbjh = nelorbj_new/2
    do i = 1, nelorbj_new
        do j = 1, nelorbj_new
            if (mapj(i) .ne. 0 .and. mapj(j) .ne. 0) then
                if (i .gt. nelorbjh .and. j .le. nelorbjh .or.&
                        &   i .le. nelorbjh .and. j .gt. nelorbjh) then
                    jastrow_new(i, j) = jastrow_new(i, j) - jastrow_old(abs(mapj(i)), abs(mapj(j)))
                else
                    jastrow_new(i, j) = jastrow_new(i, j) + jastrow_old(abs(mapj(i)), abs(mapj(j)))
                end if
            else
                jastrow_new(i, j) = 0.d0
            end if
        end do
    end do
    return
end
