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

subroutine makelambda(norb, symrot, nrot, symtra, ntra, occupied, recordsym &
                      , lenrec, nrec, symmagp, ipsip, rion, ntotatoms, orb2atom &
                      , cellscale, deps, yes_hermite, opposite_phase, zeta, J3_off &
                      , check_same_atom, no_4body_jas, cut_hybrid)
    use symm_data, only: CartesianToCrystal
    implicit none
    integer norb, nrot, ntra, symrot(norb, nrot), symtra(norb, ntra)
    logical occupied(norb, norb), symmagp, accept, yes_hermite, J3_off(ntotatoms), check_same_atom, yes_cut
    integer recordsym(2, ntra*nrot, *), lenrec(*), i, j, k, l, jfirst, imap, jmap&
            &, signi, signj, imapr, jmapr, signir, signjr, ipsip(2, ntra*nrot), indsym, countrec&
            &, nrec, count, kk, ntotatoms
    real*8 rion(3, ntotatoms), rdiff(3), rdiff_new(3), rdiff_rot(3), cellscale(3), cost, deps, zeta(2, ntotatoms)
    integer orb2atom(*), posi, posj, posi_new, posj_new
    integer cut_hybrid(ntotatoms)
    logical phase_locked, phase_bound, yesfirst, opposite_phase, no_4body_jas
    integer, dimension(:), allocatable :: atbasis, refat
    logical, external :: check_longbond

    !  Input symrot(norb,nrot)
    !   j=symrot(i,is)   is=1,nrot are the possible rotation symmetries
    !  |j|>0  is the orbital obtained after application of the symmetry is
    !  to the orbital i. if symrot(i,is)=0  it is not possible to apply this
    !  symmetry to the orbital i (generates more than a single orbital).
    !  the sign of j is the sign obtained after this transformation to the
    !  orbital including the sign of the boundary conditions.

    !   symtra(norb,ntra) as above for the translation symmetries.
    !   assumed that the translation symmetries are always allowed (never symtra(i,j)=0.
    !   Output occupied=.true. for an allowed matrix element of the
    !   lambda_ij
    !   output record(2,1:lenrec(k),k) k<=nrec, if nrec large enough all matrix
    !   elements are classified otherwise a Warning appears and nrec is the number
    !   of records that are actually required.
    !   lenrec(k)  is the length of the record #k

    occupied = .false.

    phase_locked = .false.
    countrec = 0
    yes_cut = .false.
    do i = 1, ntotatoms
        if (cut_hybrid(i) .ne. -1) yes_cut = .true.
    end do

    if (yes_cut .or. yes_hermite) then
        allocate (atbasis(ntotatoms), refat(ntotatoms))
        atbasis = 0
        do i = 1, norb
            atbasis(orb2atom(i)) = atbasis(orb2atom(i)) + 1
        end do
        ! computation cumulative index
        refat = atbasis
        do i = 2, ntotatoms
            refat(i) = refat(i) + refat(i - 1)
        end do
        do i = 1, ntotatoms
            refat(i) = refat(i) - atbasis(i)
        end do
        ! compute the first reference
        !  write(6,*) ' HERE refat '
        !  do i=1,ntotatoms
        !  write(6,*) i,refat(i)
        !  enddo
    end if

    !  write(6,*) ' HERE orb2atom ',yes_hermite
    !  do i=1,norb
    !  write(6,*) i,orb2atom(i)
    !  enddo

    do i = 1, norb
        if (symmagp) then
            jfirst = i
        else
            jfirst = 1
        end if
        do j = jfirst, norb
            posi = orb2atom(i)
            posj = orb2atom(j)
            if (yes_hermite) then
                rdiff(:) = rion(:, posi) - rion(:, posj)
                call CartesianToCrystal(rdiff, 1)
                call makeimage(rdiff, cellscale, deps)
                phase_locked = .false.
                phase_bound = .false.
            end if

            if (.not. occupied(i, j)) then
                if (check_same_atom .and. orb2atom(i) .eq. orb2atom(j)) then
                    ! avoid J1
                    if (orb2atom(j) .gt. 0) then
                        if (J3_off(orb2atom(i))) then
                            !write(6,*) "debug XXX ", i, j, orb2atom(i), J3_off(orb2atom(i))
                            cycle
                        end if
                    end if
                end if
                indsym = 0
                do k = 1, ntra
                    signi = 1
                    signj = 1
                    imap = symtra(i, k)
                    jmap = symtra(j, k)
                    if (imap .lt. 0) then
                        imap = -imap
                        signi = -signi
                    end if
                    if (jmap .lt. 0) then
                        jmap = -jmap
                        signj = -signj
                    end if

                    if (imap .eq. 0 .or. jmap .eq. 0) then
                        write (6, *) ' Error in symtra !!! ', imap, jmap
                        stop
                    end if
                    do l = 1, nrot
                        imapr = symrot(imap, l)
                        jmapr = symrot(jmap, l)

                        if (imapr .lt. 0) then
                            imapr = -imapr
                            signir = -signi
                        else
                            signir = signi
                        end if
                        if (jmapr .lt. 0) then
                            jmapr = -jmapr
                            signjr = -signj
                        else
                            signjr = signj
                        end if

                        if (imapr .ne. 0 .and. jmapr .ne. 0) then ! if it is an allowed symmetry
                            if (yes_hermite) then
                                if (symmagp) then
                                    posi = orb2atom(abs(imapr))
                                    posj = orb2atom(abs(jmapr))
                                    ! rdiff_rot remains n.n. image, apart for boundaries.
                                    rdiff_rot(:) = rion(:, posi) - rion(:, posj)
                                    call CartesianToCrystal(rdiff_rot, 1)
                                    call makeimage(rdiff_rot, cellscale, deps)
                                    phase_bound = check_longbond(rdiff_rot, cellscale, deps)
                                    !          Only for the same type of orbital and the same type of ion
                                    if (opposite_phase) then
                                        if (phase_bound .and. (abs(imapr) - refat(posi)) .eq. (abs(jmapr) - refat(posj))&
                                                &.and. opposite_phase .and. zeta(1, orb2atom(i)) .eq. zeta(1, orb2atom(j))&
                                                &.and. posi .ne. posj) then
                                            phase_locked = .true.
                                        end if
                                    end if
                                end if

                            end if ! endif conj
                            indsym = indsym + 1
                            ipsip(1, indsym) = imapr*signjr*signir
                            ipsip(2, indsym) = jmapr
                        end if
                    end do
                end do
                !       check if it is a compatible record
                if (indsym .gt. 0) then
                    accept = .true.
                    do k = 1, indsym
                        occupied(abs(ipsip(1, k)), abs(ipsip(2, k))) = .true.
                    end do
                    if (symmagp) then
                        do k = 1, indsym
                            occupied(abs(ipsip(2, k)), abs(ipsip(1, k))) = .true.
                        end do
                    end if
                else
                    accept = .false.
                end if
                if (no_4body_jas) then
                    if (orb2atom(i) .ne. orb2atom(j) .and. orb2atom(i) .ne. 0 .and. orb2atom(j) .ne. 0) accept = .false.
                end if
                if (yes_cut) then
                    if (posi .eq. posj .and. cut_hybrid(posi) .ne. -1) then
                        if (abs(imapr) - refat(posi) .gt. cut_hybrid(posi) &
                            .or. abs(jmapr) - refat(posj) .gt. cut_hybrid(posj)) then
                            accept = .false.
                        end if
                    end if
                    if (posi .ne. posj .and. cut_hybrid(posi) .ne. -1 .and. cut_hybrid(posj) .ne. -1) then
                        if (abs(jmapr) - refat(posj) .le. cut_hybrid(posj)) accept = .false.
                        if (abs(imapr) - refat(posi) .le. cut_hybrid(posi)) accept = .false.
                    end if
                    if (posi .ne. posj .and. cut_hybrid(posi) .ne. -1 .and. cut_hybrid(posj) .eq. -1) then
                        if (abs(imapr) - refat(posi) .le. cut_hybrid(posi)) accept = .false.
                    end if
                    if (posi .ne. posj .and. cut_hybrid(posi) .eq. -1 .and. cut_hybrid(posj) .ne. -1) then
                        if (abs(jmapr) - refat(posj) .le. cut_hybrid(posj)) accept = .false.
                    end if
                end if

                !           if(i.eq.1.and.j.eq.54) then
                !            write(6,*) ' Path 1-54 ',phase_bound,rdiff(1:3)
                !            do k=1,indsym
                !            write(6,*) k,ipsip(1,k),ipsip(2,k)
                !            enddo
                !           endif

                !          if(.not.yes_hermite) then  ! no parameter eliminated in the complex case.
                if (symmagp) then
                    do k = 1, indsym
                        if ((abs(ipsip(2, k)) .eq. i .and. ipsip(1, k) .eq. -j) .or. (ipsip(1, k) .eq. -i &
                                                                                      .and. abs(ipsip(2, k)) .eq. j)) then
                            accept = .false.
                        end if
                    end do
                else
                    do k = 1, indsym
                        if (abs(ipsip(2, k)) .eq. j .and. ipsip(1, k) .eq. -i) accept = .false.
                    end do
                end if

                if (accept) then
                    countrec = countrec + 1
                    if (countrec .le. nrec) then
                        !    group toghether the various pairs
                        lenrec(countrec) = 1
                        recordsym(1:2, 1, countrec) = ipsip(1:2, 1)
                        do k = 2, indsym
                            accept = .true.
                            if (symmagp) then
                                do l = 1, k - 1
                                    if ((abs(ipsip(1, k)) .eq. abs(ipsip(1, l)) .and. abs(ipsip(2, k)) .eq. abs(ipsip(2, l))) .or.&
                                      &(abs(ipsip(1, k)) .eq. abs(ipsip(2, l)) .and. abs(ipsip(2, k)) .eq. abs(ipsip(1, l)))) then
                                        accept = .false.
                                    end if
                                end do
                            else
                                do l = 1, k - 1
                                    if ((abs(ipsip(1, k)) .eq. abs(ipsip(1, l)) .and. abs(ipsip(2, k)) .eq. abs(ipsip(2, l)))) then
                                        accept = .false.
                                    end if
                                end do
                            end if
                            if (accept) then
                                lenrec(countrec) = lenrec(countrec) + 1
                                recordsym(1:2, lenrec(countrec), countrec) = ipsip(1:2, k)
                            end if
                        end do
                        if (yes_hermite) then
                            do l = 1, lenrec(countrec)
                                if (recordsym(2, l, countrec) .lt. 0) then
                                    ! Forget about the phase as it is taken into account by the attaching flux.
                                    recordsym(1, l, countrec) = -abs(recordsym(1, l, countrec))
                                    recordsym(2, l, countrec) = abs(recordsym(2, l, countrec))
                                else
                                    recordsym(1:2, l, countrec) = abs(recordsym(1:2, l, countrec))
                                end if
                            end do
                            if (phase_locked .and. opposite_phase) then
                                recordsym(2, 1:lenrec(countrec), countrec) = -abs(recordsym(2, 1:lenrec(countrec), countrec))
                            end if
                        end if
                    end if ! if countrec<=nrec
                end if ! endif accepted
            end if ! endif occupied
        end do
    end do

    if (countrec .gt. nrec) then
        write (6, *) ' Warning # record too small ', countrec, nrec
    else

        occupied = .false.
        count = 0
        do i = 1, min(nrec, countrec)
            do j = 1, lenrec(i)

                occupied(abs(recordsym(1, j, i)), abs(recordsym(2, j, i))) = .true.

                do k = 1, lenrec(i)

                    if (abs(recordsym(1, j, i)) .eq. abs(recordsym(1, k, i)) .and. abs(recordsym(2, j, i)) .eq.&
                            &abs(recordsym(2, k, i)) .and. k .ne. j) then

                        write (6, *) ' Repetition !!! j,k, record  ', j, k, i
                        write (6, *) ' first pair --> ', recordsym(1, j, i), recordsym(2, j, i)
                        write (6, *) ' second pair --> ', recordsym(1, k, i), recordsym(2, k, i)
                        write (6, *) lenrec(i), (recordsym(1, kk, i), recordsym(2, kk, i), kk=1, lenrec(i))
                        stop

                    end if

                end do

            end do
            count = count + lenrec(i)
        end do
        write (6, *) ' Matrix elements inside makelambda =', count
        count = 0
        do i = 1, norb
            do j = 1, norb
                if (occupied(i, j)) count = count + 1
            end do
        end do

        write (6, *) ' Count check inside makelambda ', count

    end if
    nrec = countrec
    if (yes_cut .or. yes_hermite) deallocate (atbasis, refat)

    return
end subroutine makelambda

function check_longbond(rdiff, cellscale, deps)
    implicit none
    logical check_longbond
    integer i, check
    real*8 rdiff(3), cellscale(3), deps
    check = 0
    do i = 1, 3
        if (abs(2*abs(rdiff(i)) - cellscale(i)) .lt. deps) then
            check = check - 1
        elseif (abs(rdiff(i)) .gt. deps) then
            check = check + 4
        end if
    end do
    check_longbond = .false.
    if (check .lt. 0) check_longbond = .true.
end function check_longbond

subroutine makeimageo_fake(rdiff, cellscale, deps)
    implicit none
    integer kk, npip, m
    real*8 cost, rdiff(3), cellscale(3), deps
    do kk = 1, 3
        cost = rdiff(kk)/cellscale(kk)
        npip = nint(2*cost)
        if ((npip/2)*2 .ne. npip .and. abs(cost - npip/2.d0) .lt. deps) then ! Boarder cases
            m = (npip + 1)/2
        else
            m = nint(cost)
        end if
        rdiff(kk) = rdiff(kk) - cellscale(kk)*m
    end do
    return
end
