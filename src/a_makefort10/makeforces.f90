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

subroutine makeforces(ntotatoms, symrot, nrot, nrotu, symtra, ntra, ntrau, occupied, recordsym&
        &, lenrec, nrec, ipsip)
    implicit none
    integer norb, ntotatoms, nrot, ntra, symrot(3*ntotatoms, nrot)&
            &, symtra(3*ntotatoms, ntra), ntrau, nrotu
    logical occupied(3*ntotatoms), accept
    integer recordsym(2, ntra*nrot, *), lenrec(*), i, j, k, l, jfirst, imap, jmap&
            &, signi, signj, imapr, jmapr, signir, signjr, ipsip(2, ntra*nrot), indsym, countrec&
            &, nrec, count, kk, iatom, icomp

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

    countrec = 0

    norb = 3*ntotatoms

    do i = 1, norb
        if (.not. occupied(i)) then
            indsym = 0
            do k = 1, ntrau
                signi = 1
                imap = symtra(i, k)
                if (imap .lt. 0) then
                    imap = -imap
                    signi = -signi
                end if
                if (imap .eq. 0) then
                    write (6, *) ' Error in symtra !!! ', imap, jmap
                    stop
                end if
                do l = 1, nrotu
                    imapr = symrot(imap, l)
                    if (imapr .lt. 0) then
                        imapr = -imapr
                        signir = -signi
                    else
                        signir = signi
                    end if
                    if (imapr .ne. 0) then ! if it is an allowed symmetry
                        indsym = indsym + 1
                        iatom = (imapr - 1)/3 + 1
                        icomp = imapr - (iatom - 1)*3
                        ipsip(1, indsym) = iatom*signir
                        ipsip(2, indsym) = icomp

                    end if
                end do
            end do
            !       check if it is a compatible record

            if (indsym .gt. 0) then

                accept = .true.
                do k = 1, indsym
                    occupied(3*(abs(ipsip(1, k)) - 1) + ipsip(2, k)) = .true.
                end do
                !       Reject record if I can go back to the same component changing sign.
                do k = 2, indsym
                    if (ipsip(1, k) .eq. -ipsip(1, 1) .and. ipsip(2, k) .eq. ipsip(2, 1)) accept = .false.
                end do
            else
                accept = .false.
            end if

            !        do k=2,indsym
            !    if(ipsip(1,k).eq.-ipsip(1,1).and.ipsip(2,k).eq.ipsip(2,1)) accept=.false.
            !        enddo

            if (accept) then
                countrec = countrec + 1
                if (countrec .le. nrec) then
                    !    group toghether the various pairs
                    lenrec(countrec) = 1
                    recordsym(1:2, 1, countrec) = ipsip(1:2, 1)
                    do k = 2, indsym
                        accept = .true.
                        do l = 1, k - 1
                            if (abs(ipsip(1, l)) .eq. abs(ipsip(1, k)) .and. ipsip(2, l) .eq. ipsip(2, k))&
                                    & accept = .false.
                        end do
                        if (accept) then
                            lenrec(countrec) = lenrec(countrec) + 1
                            recordsym(1:2, lenrec(countrec), countrec) = ipsip(1:2, k)
                        end if
                    end do
                end if ! if countrec<=nrec
            end if ! endif accepted
        end if ! endif occupied
    end do

    nrec = countrec

    return

end

