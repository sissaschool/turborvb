! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2004 PWSCF-CP-FPMD group
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
! It builds the matrix for the pfaffian combining the makelambda imported
! from a previour run (in the makefort10) and uses a similar procedure for
! the skew symmetric part, then it combines the two matrices obtained to
! create the wf
!

!
! The first part of the code builds the skew symmetric block of the matrix
! Run 2 times and for the second execution allocate recordsym
!

subroutine makeskew(norb, symrot, nrot, symtra, ntra, occupied, recordsym, &
                    lenrec, nrec, ipsip, rion, ntotatoms, orb2atom, cellscale, deps, zeta, yes_hermite)
    implicit none

    integer norb, nrot, ntra, symrot(norb, nrot), symtra(norb, ntra)
    logical occupied(norb, norb), accept, sameorb, yes_hermite
    integer recordsym(2, ntra*nrot, *), lenrec(*), i, j, k, l, jfirst, imap, jmap&
            &, signi, signj, imapr, jmapr, signir, signjr, ipsip(2, ntra*nrot), indsym, countrec&
            &, nrec, count, kk, ntotatoms, check
    real*8 rion(3, ntotatoms), rdiff(3), rdiff_new(3), rdiff_rot(3), rdiff_roti(3), cellscale(3), cost, deps, zeta(2, ntotatoms)
    integer orb2atom(*), posi, posj, posi_new, posj_new
    integer, dimension(:), allocatable :: atbasis, refat

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
    !   elements are classified otherwise a Warning appears and rec is the number
    !   of records that are actually required.
    !   lenrec(k)  is the length of the record #k

    occupied = .false.
    countrec = 0

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

    do i = 1, norb
        jfirst = i + 1
        do j = jfirst, norb
            posi = orb2atom(i)
            posj = orb2atom(j)
            rdiff(:) = rion(:, posi) - rion(:, posj)
            call makeimage(rdiff, cellscale, deps)
            accept = .false.
            if (.not. occupied(i, j)) then
                accept = .true.
                indsym = 0
                do k = 1, ntra
                    signi = 1
                    signj = 1
                    imap = symtra(i, k)
                    jmap = symtra(j, k)
                    !Check with Sandro
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

                        !Only for the same type of orbital and the same type of ion
                        !Checking if the two atoms are at distance cellscale/2 and with the same orbital
                        if (imapr .ne. 0 .and. jmapr .ne. 0) then
                            if (yes_hermite) then ! if it is an allowed symmetry
                                posi = orb2atom(abs(imapr))
                                posj = orb2atom(abs(jmapr))
                                rdiff(:) = rion(:, posi) - rion(:, posj)
                                sameorb = .false.
                                call makeimage(rdiff, cellscale, deps)

                                if ((abs(imapr) - refat(posi)) .eq. (abs(jmapr) - refat(posj)) &
                                    .and. zeta(1, orb2atom(i)) .eq. zeta(1, orb2atom(j)) .and. posi .ne. posj) sameorb = .true.

                                !If the distance of the atoms has only component L/2 and 0 (at least one has to be L/2)
                                if (sameorb) then
                                    check = 0
                                    do kk = 1, 3
                                        if (abs(2*abs(rdiff(kk)) - cellscale(kk)) .lt. deps) then
                                            check = check - 1
                                        else if (abs(2*abs(rdiff(kk))) .gt. deps) then
                                            check = check + 4
                                        end if
                                    end do
                                    if (check .lt. 0) accept = .false.
                                end if
                            end if

                            if (accept) then
                                indsym = indsym + 1
                                ipsip(1, indsym) = imapr*signjr*signir
                                ipsip(2, indsym) = jmapr
                            end if
                        end if
                    end do
                end do

                !       check if it is a compatible record
                if (indsym .gt. 0 .and. accept) then
                    accept = .true.
                    do k = 1, indsym
                        occupied(abs(ipsip(1, k)), abs(ipsip(2, k))) = .true.
                        occupied(abs(ipsip(2, k)), abs(ipsip(1, k))) = .true.
                    end do
                else
                    accept = .false.
                end if

                !Checking if some elements cannot be accepted because they are zero by definition
                do k = 1, indsym
                    if ((abs(ipsip(2, k)) .eq. i .and. ipsip(1, k) .eq. j) .or. &
                        (ipsip(1, k) .eq. -i .and. abs(ipsip(2, k)) .eq. j)) then
                        accept = .false.
                    end if
                end do

                !Checking if it is already occupied
                if (accept) then
                    countrec = countrec + 1
                    if (countrec .le. nrec) then
                        !    group toghether the various pairs

                        lenrec(countrec) = 1
                        recordsym(1:2, 1, countrec) = ipsip(1:2, 1)
                        do k = 2, indsym
                            accept = .true.

                            do l = 1, k - 1
                                if ((abs(ipsip(1, k)) .eq. abs(ipsip(1, l)) .and. abs(ipsip(2, k)) .eq. abs(ipsip(2, l))) .or.&
                                        &(abs(ipsip(1, k)) .eq. abs(ipsip(2, l)) .and. abs(ipsip(2, k)) .eq. abs(ipsip(1, l)))) &
                                        & accept = .false.
                            end do

                            if (accept) then
                                lenrec(countrec) = lenrec(countrec) + 1
                                recordsym(1:2, lenrec(countrec), countrec) = ipsip(1:2, k)
                            end if
                        end do

                    end if ! if countrec<=nrec
                end if ! endif accepted
            end if
        end do

    end do
    deallocate (atbasis, refat)

    write (6, *) "SKEW"
    do i = 1, nrec
        write (6, *) "step", i, "lenrec(step)", lenrec(i)
        do j = 1, lenrec(i)
            write (6, *) recordsym(1, j, i), recordsym(2, j, i)
        end do
    end do

    nrec = countrec

end subroutine makeskew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Assembling the pfaffian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine makepfaff(norb, nrot, ntra, occupied, recordsymagp, lenrecagp, recordsymskw, lenrecskw, recordsym, &
                     lenrec, nrec, nrecagp, nrecskw, nouppfaff, nodownpfaff, opposite_phase)
    implicit none
    integer :: norb, nrot, ntra, nrec, nrecagp, nrecskw, recordymmagp, count, i, j, k, kk
    logical :: occupied(2*norb, 2*norb), nodownpfaff, nouppfaff, opposite_phase
    integer :: recordsym(2, ntra*nrot, nrec), recordsymagp(2, ntra*nrot, nrecagp), recordsymskw(2, ntra*nrot, nrecskw), &
               lenrec(nrec), lenrecagp(nrecagp), lenrecskw(nrecskw)

    recordsym = 0

    !Building a new and larger recordsym and lenrec
    if (nodownpfaff .and. .not. nouppfaff) then
        recordsym(:, 1:nrot*ntra, 1:nrecskw) = recordsymskw(:, 1:nrot*ntra, 1:nrecskw)
        recordsym(1, 1:nrot*ntra, nrecskw + 1:nrecskw + nrecagp) = recordsymagp(1, 1:nrot*ntra, 1:nrecagp)
        recordsym(2, 1:nrot*ntra, nrecskw + 1:nrecskw + nrecagp) = recordsymagp(2, 1:nrot*ntra, 1:nrecagp) &
                                                                   + sign(1, recordsymagp(2, 1:nrot*ntra, 1:nrecagp))*norb
        lenrec(1:nrecskw) = lenrecskw(1:nrecskw)
        lenrec(1 + nrecskw:nrecagp + nrecskw) = lenrecagp(1:nrecagp)
        nrec = nrecagp + nrecskw

    else if (nouppfaff .and. .not. nodownpfaff) then
        recordsym(:, 1:nrot*ntra, 1:nrecskw) = recordsymskw(:, 1:nrot*ntra, 1:nrecskw) &
                                               + norb*sign(1, recordsymskw(:, 1:nrot*ntra, 1:nrecskw))
        recordsym(1, 1:nrot*ntra, nrecskw + 1:nrecskw + nrecagp) = recordsymagp(1, 1:nrot*ntra, 1:nrecagp)
        recordsym(2, 1:nrot*ntra, nrecskw + 1:nrecskw + nrecagp) = recordsymagp(2, 1:nrot*ntra, 1:nrecagp) &
                                                                   + norb*sign(1, recordsymagp(2, 1:nrot*ntra, 1:nrecagp))
        lenrec(1:nrecskw) = lenrecskw(1:nrecskw)
        lenrec(1 + nrecskw:nrecagp + nrecskw) = lenrecagp(1:nrecagp)
        nrec = nrecagp + nrecskw
        !  else if (symmagp) then
        !     write (6,*) "Symmagp is active with same/opposite  phase"
        !     do i=1, nrecskw
        !       recordsym(:, 1:lenrecskw(i) , i )=recordsymskw(:, 1:lenrecskw(i), i)
        !        recordsym(:, 1+lenrecskw(i):2*lenrecskw(i), i)= recordsymskw(:, 1:lenrecskw(i), i)&
        !        +norb*sign(1,recordsymskw(:, 1:lenrecskw(i), i))
        !     end do
        !     recordsym(1,1:nrot*ntra ,nrecskw+1:nrecskw+nrecagp)=recordsymagp(1,1:nrot*ntra,1:nrecagp)
        !     recordsym(2,1:nrot*ntra ,nrecskw+1:nrecskw+nrecagp)=recordsymagp(2,1:nrot*ntra,1:nrecagp)&
        !      + norb*sign(1,recordsymagp(2,1:nrot*ntra,1:nrecagp))
        !      lenrec(1:nrecskw)=lenrecskw(1:nrecskw)
        !     lenrec(1+nrecskw:nrecagp+nrecskw)=lenrecagp(1:nrecagp)
        !     nrec=nrecagp+nrecskw
    elseif (nouppfaff .and. nodownpfaff) then
        recordsym(1, 1:nrot*ntra, 1:nrecagp) = recordsymagp(1, 1:nrot*ntra, 1:nrecagp)
        recordsym(2, 1:nrot*ntra, 1:nrecagp) = recordsymagp(2, 1:nrot*ntra, 1:nrecagp) &
                                               + norb*sign(1, recordsymagp(2, 1:nrot*ntra, 1:nrecagp))
        lenrec(1:nrecagp) = lenrecagp(1:nrecagp)
        nrec = nrecagp
    else
        recordsym(:, 1:nrot*ntra, 1:nrecskw) = recordsymskw(:, 1:nrot*ntra, 1:nrecskw)
        recordsym(:, 1:nrot*ntra, nrecskw + 1:nrecskw*2) = recordsymskw(:, 1:nrot*ntra, 1:nrecskw) &
                                                           + norb*sign(1, recordsymskw(:, 1:nrot*ntra, 1:nrecskw))
        recordsym(1, 1:nrot*ntra, nrecskw*2 + 1:nrecskw*2 + nrecagp) = recordsymagp(1, 1:nrot*ntra, 1:nrecagp)
        recordsym(2, 1:nrot*ntra, nrecskw*2 + 1:nrecskw*2 + nrecagp) = recordsymagp(2, 1:nrot*ntra, 1:nrecagp) &
                                                                       + norb*sign(1, recordsymagp(2, 1:nrot*ntra, 1:nrecagp))
        lenrec(1:nrecskw) = lenrecskw(1:nrecskw)
        lenrec(1 + nrecskw:nrecskw*2) = lenrecskw(1:nrecskw)
        lenrec(1 + nrecskw*2:nrecagp + nrecskw*2) = lenrecagp(1:nrecagp)
        nrec = nrecagp + 2*nrecskw
    end if

    !   write (6,*) "Inside makepfaff"

    !   do i=1, nrec
    !      write (6,*) "step", i  ,"lenrec(step)",lenrec(i)
    !     do j=1, lenrec(i)
    !         write (6,*)  recordsym(1,j,i), recordsym(2,j,i)
    !      end do
    !   end do

    occupied = .false.
    count = 0
    do i = 1, nrec
        !    write (6,*) i,"lenrec", lenrec(i)
        do j = 1, lenrec(i)
            occupied(abs(recordsym(1, j, i)), abs(recordsym(2, j, i))) = .true.
            !       occupied(abs(recordsym(2,j,i)),abs(recordsym(1,j,i)))=.true.

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
        !   write (6,*) "Passi pfaff 2"
    end do
    count = 0
    do i = 1, 2*norb
        do j = 1, 2*norb
            if (occupied(i, j)) count = count + 1
        end do
    end do

    write (6, *) ' Count check inside makepfaff ', count
    write (6, *) 'record found  for the agp part of the pfaffian=', nrecagp
    write (6, *) 'record found  for the skew part of the pfaffian=', nrecskw
    write (6, *) ' record found =', nrec
    return
end subroutine makepfaff
