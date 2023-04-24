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

! This subroutine is called at the beginning of the self-consistent cycle
! and it reads the molecular orbitals from the input real or complex wave
! function.
!
! Remember: ipc = 1 --> real wave function, real eigenvectors molecorb
!           ipc = 2 --> complex wave function, complex eigenvectors molecorb

subroutine read_molecular

    use allio
    use constants, only: ipc
    use setup, only: yeslsda, newocc, bands, &
                     molecorb_old, molecorbdo_old
    use freeelmod_complex, only: nozeroc_in, nnozeroc_in, jbradetc_in, jbradetnc_in, nelcolc_in, &
                                 nelorbc_in, niesd_in, iesdrr_in, contraction_in, nmoltot, count2n, &
                                 countalln, costalln

    implicit none
    ! local
    integer :: indpar, i, j, k, i_new, j_new, maxup, nunpn, ind, indx, indy
    logical :: done

    ! by definition the unpaired are the last molecular orbitals
    indpar = iesup_c - 2*molecular*nelorbh
    ! read also occupations from detmat_c
    done = .true.
    newocc = 0.d0
    count2n = 0
    countalln = 0
    ind = 0
    !
    ! save initial w.f. important parameters
    !
    niesd_in = niesd ! 1/2 Jastrow parameters
    iesdrr_in = iesdrr ! 1/2 Jastrow type
    contraction_in = contraction ! number of contracted orbitals (atomic+molecular)
    nelcolc_in = nelcol_c ! # of contracted parameters
    nelorbc_in = nelorb_c
    allocate (nozeroc_in(nnozero_c), jbradetc_in(nnozero_c), jbradetnc_in(3*nnozero_c))
    if (contraction .gt. 0) then
        nnozeroc_in = nnozero_c ! # of non-zero elements in the determinant
        nozeroc_in = nozero_c ! addresses of non-zero elements of the determinant
        jbradetc_in = jbradet ! symmetries of the determinant
        jbradetnc_in = jbradetn ! parameters of the determinant locked by symmetry
    else
        nnozeroc_in = nnozero_c
        nozeroc_in = nozero
        jbradetc_in = jbradet
        jbradetnc_in = jbradetn
    end if

    molecorb_old = 0.d0
    if (yeslsda .or. ipc .eq. 2) molecorbdo_old = 0.d0
    if (rank .eq. 0) write (6, '(a,5I6)') ' Read molecular orbitals', nmol, nmoltot, nelorb_c, nelorbh, bands

    do i = 1, nmol
        if (i .le. bands) then
            if (symmagp .and. ipc .eq. 1) then
                indx = nelorb_c - nmoltot + i
                indy = nelorb_c - nmoltot + i
            else
                indx = nelorb_c - nmoltot + 2*i
                indy = indx - 1
            end if
            ind = nelorb_c*(indy - 1) + indx
            if (ipc .eq. 2) ind = 2*ind - 1 ! Real part of detmat_c
            newocc(i) = detmat_c(ind)

            if (detmat_c(ind) .ne. 0.d0) countalln = countalln + 1
            if (detmat_c(ind) .eq. 1.d0) count2n = count2n + 1
            !
            ! real case with symmetric AGP
            !
            if (symmagp .and. ipc .eq. 1) then
                do j = nelorbh + 1, 2*nelorbh
                    molecorb_old(j - nelorbh, i) = dup_c(indpar + j)
                end do
            end if
            !
            ! non-symmetric AGP for real case (spin calculation) or complex case
            !
            if (.not. symmagp .or. ipc .eq. 2) then
                if (ipc .eq. 1) then ! real case: read eigenvectors for down-spin electrons only
                    ! in the case of a LSDA calculation
                    if (yeslsda) then
                        do j = nelorbh + 1, 2*nelorbh
                            molecorbdo_old(j - nelorbh, i) = dup_c(indpar + j)
                        end do
                        do j = nelorbh + 1, 2*nelorbh
                            molecorb_old(j - nelorbh, i) = dup_c(indpar + j + 2*nelorbh)
                        end do
                    else
                        do j = nelorbh + 1, 2*nelorbh
                            molecorb_old(j - nelorbh, i) = dup_c(indpar + j)
                        end do
                    end if
                else ! complex case: always read different eigenvectors
                    ! for up/down spin electrons.
                    do j = nelorbh + 1, 2*nelorbh ! down
                        molecorbdo_old(2*(j - nelorbh) - 1, i) = dup_c(2*(indpar + j) - 1) ! Re
                        molecorbdo_old(2*(j - nelorbh), i) = dup_c(2*(indpar + j)) ! Im
                    end do
                    do j = nelorbh + 1, 2*nelorbh ! up
                        molecorb_old(2*(j - nelorbh) - 1, i) = dup_c(2*(indpar + j) - 1 + 4*nelorbh) ! Re
                        molecorb_old(2*(j - nelorbh), i) = dup_c(2*(indpar + j) + 4*nelorbh) ! Im
                    end do
                end if
            end if
        else
            if (rank .eq. 0 .and. done) &
                write (6, '(a)') ' Warning not all orbitals are read from fort10'
            done = .false.
        end if
        indpar = indpar + 2*nelorbh
        if (.not. symmagp .or. ipc .eq. 2) indpar = indpar + 2*nelorbh

    end do

    if (ndiff .gt. 0 .and. rank .eq. 0) write (6, '(a)') ' Read unpaired molecular orbitals'
    do i = nmol + 1, nmol + ndiff
        i_new = i - nmol
        j_new = nelorb_c - ndiff + i_new
        ind = nelorb_c*(nelorb_c + i_new - 1) + j_new
        if (ipc .eq. 1) then
            newocc(bands + i_new) = -abs(detmat_c(ind))
        else
            newocc(bands + i_new) = -abs(detmat_c(2*ind - 1))
        end if
        ind = i - nmol + bands - ndiff
        do j = nelorbh + 1, 2*nelorbh
            if (ipc .eq. 1) then
                molecorb_old(j - nelorbh, i) = dup_c(indpar + j)
            else
                molecorb_old(2*(j - nelorbh) - 1, i) = dup_c(2*indpar + 2*j - 1) ! Re
                molecorb_old(2*(j - nelorbh), i) = dup_c(2*indpar + 2*j) ! Im
            end if
        end do
!     if(rank.eq.0)  write(6,*) i,sum(dup_c(indpar+nelorbh+1:indpar+2*nelorbh))
        indpar = indpar + 2*nelorbh
    end do
    !
    ! definition of costalln
    !
    maxup = (countalln - count2n) - (neldo - count2n)
    if (countalln .gt. count2n) then
        if (maxup .gt. ndiff) then
            nunpn = 0
            costalln = dble(neldo + nelup - 2*count2n)/dble(countalln - count2n)
        else
            nunpn = ndiff - maxup
            costalln = dble(neldo + nelup - 2*count2n - nunpn)/dble(countalln - count2n)
        end if
    else
        if (maxup .gt. ndiff) then
            nunpn = 0
            costalln = 0
        else
            nunpn = ndiff - maxup
            costalln = 0
        end if
    end if

    return

end subroutine read_molecular

! This subroutine is called at the end of the self-consistent
! cycle. It updates the wave function (written in the file fort.10_new)
! with the new molecular orbitals, eigenvectors of the converged Kohn-Sham hamiltonian.
! It also recomputes all symmetries of the orbitals and update several variables
! needed by the subroutine write_fort10().
!
! Remember: ipc = 1 --> real wave function, real eigenvectors molecorb
!           ipc = 2 --> complex wave function, complex eigenvectors molecorb

subroutine update_fort10

    use allio
    use constants, only: ipc
    use setup, only: molecorb_old, molecorbdo_old
    use setup, only: bands, yeslsda, occupations, occupationdo
    use freeelmod_complex, only: contraction_in, nelcolc_in, nelorbc_in, yest, &
                                 jbradetc_in, jbradetnc_in, nnozeroc_in, nozeroc_in, &
                                 countall, optocc, iesswrc_in, iesup_cmol, nmoltot, &
                                 iesupindmol, iesupr_cmol, indnn_in, inddiff, inds, &
                                 molecsw, nmol_in, nelocc, neloccdo, nshellmol, occ_cmol

    implicit none
    ! local
    integer :: i, j, k, ii, jj, kk, dimvect, icount1, icount2, indfirst, &
               indlast, ind, indpar, indx, indy, occmax
    integer, dimension(:), allocatable :: indexo, indexn
    real(8), dimension(:), allocatable :: sortvect
    logical found1

    ! nmoltot -> total number of molecular orbitals + unpaired electrons
    nmoltot = nmol + ndiff
    if (.not. symmagp .or. ipc .eq. 2) nmoltot = nmoltot + nmol
    !
!    deallocate(ipsip)
!    allocate(ipsip(nshell_c))
!    ipsip = 0
!    nshell_c = nshell_c - molecular
!    nshellmol = nshell_c + nmoltot

    ! mult_c -> orbitals multiplicity
!    ipsip = mult_c
!    deallocate(mult_c)
!    allocate(mult_c(nshellmol))
!    mult_c(1:nshell_c) = ipsip(1:nshell_c)
!    mult_c(nshell_c + 1:nshellmol) = 1

    ! nparam_c -> basis set parameters
!    ipsip = nparam_c
!    deallocate(nparam_c)
!    allocate(nparam_c(nshellmol))
!    nparam_c(1:nshell_c) = ipsip(1:nshell_c)
!    nparam_c(nshell_c + 1:nshellmol) = 2 * nelorbh

    ! ioptorb_c -> orbitals type
    ! ioptorb_c = 1000000 -> molecular orbital
!    ipsip = ioptorb_c
!    deallocate(ioptorb_c)
!    allocate(ioptorb_c(nshellmol))
!    ioptorb_c(1:nshell_c) = ipsip(1:nshell_c)
!    ioptorb_c(nshell_c + 1:nshellmol) = 1000000

    ! kion_c -> ion type in a particular shell
!    ipsip = kion_c
!    deallocate(kion_c)
!    allocate(kion_c(nshellmol))
!    kion_c(1:nshell_c) = ipsip(1:nshell_c)
!    kion_c(nshell_c + 1:nshellmol) = 1

    ! redefinition of matrix dup_c with the
    ! new MOs coefficients taken from the KS cycle.
!    if(allocated(psip)) deallocate(psip)
!    allocate(psip(ipc * iesupr_c))
!    psip(1:ipc * iesupr_c) = dup_c(1:ipc * iesupr_c)
!    deallocate(dup_c)
!    if(molecular.eq.0) then
!        iesup_cmol = iesup_c + 2 * nmoltot * nelorbh
!    else
!        iesup_cmol = iesup_c + 2 * (nmoltot - molecular) * nelorbh
!    endif
!    iesupr_cmol = iesup_cmol
!    if(yeszagp) iesupr_cmol = 2 * iesupr_cmol
!    allocate(dup_c(ipc * iesupr_cmol))
!    dup_c = 0.d0
    ! select the MOs part of dup_c matrix only and
    ! do not change the atomic one.
    indpar = iesup_c - 2*molecular*nelorbh
!    dup_c(1:ipc * indpar) = psip(1:ipc * indpar)
    !
    ! NB: in the case of a complex wave function, always consider the AGP
    ! non symmetric (as symmagp=.false.) and therefore I have different
    ! eigenvectors for spin up/down electrons.
    !
    do i = 1, nmol

        do j = 1, nelorbh
            if (ipc .eq. 1) then
                dup_c(indpar + j) = j
            else
                dup_c(2*(indpar + j) - 1) = j
                dup_c(2*(indpar + j)) = 0.d0
            end if
        end do

        if (i .le. bands) then ! fill only MOs smaller than bands value

            if (.not. yeslsda .and. ipc .eq. 1) then

                do j = nelorbh + 1, 2*nelorbh
                    dup_c(indpar + j) = molecorb_old(j - nelorbh, i)
                end do
                ! for real functions the spin-down part for non-symmetric lambda
                ! is equal to the spin-up part if no LSDA functional is used
                if (.not. symmagp) then
                    indpar = indpar + 2*nelorbh
                    do j = 1, nelorbh
                        dup_c(indpar + j) = j
                    end do
                    do j = nelorbh + 1, 2*nelorbh
                        dup_c(indpar + j) = molecorb_old(j - nelorbh, i)
                    end do
                end if

            else ! LSDA calculation or complex wave function, different eigenvectors for up/down electrons

                if (ipc .eq. 1) then
                    do j = nelorbh + 1, 2*nelorbh
                        dup_c(indpar + j) = molecorbdo_old(j - nelorbh, i)
                    end do
                else
                    do j = nelorbh + 1, 2*nelorbh
                        dup_c(2*(indpar + j) - 1) = molecorbdo_old(2*(j - nelorbh) - 1, i) ! Re
                        dup_c(2*(indpar + j)) = molecorbdo_old(2*(j - nelorbh), i) ! Im
                    end do
                end if

                indpar = indpar + 2*nelorbh
                if (ipc .eq. 1) then
                    do j = 1, nelorbh
                        dup_c(indpar + j) = j
                    end do
                    do j = nelorbh + 1, 2*nelorbh
                        dup_c(indpar + j) = molecorb_old(j - nelorbh, i)
                    end do
                else
                    do j = 1, nelorbh
                        dup_c(2*(indpar + j) - 1) = j
                        dup_c(2*(indpar + j)) = 0.d0
                    end do
                    do j = nelorbh + 1, 2*nelorbh
                        dup_c(2*(indpar + j) - 1) = molecorb_old(2*(j - nelorbh) - 1, i) ! Re
                        dup_c(2*(indpar + j)) = molecorb_old(2*(j - nelorbh), i) ! Im
                    end do
                end if

            end if

        else ! MOs index greater than number of bands required in input

            if (ipc .eq. 1) then
                do j = 1, nelorbh
                    dup_c(indpar + j) = j
                end do
                do j = nelorbh + 1, 2*nelorbh
                    dup_c(indpar + j) = 0.d0
                end do
            else
                do j = 1, nelorbh
                    dup_c(2*(indpar + j) - 1) = j
                    dup_c(2*(indpar + j)) = 0.d0
                end do
                do j = 2*nelorbh + 1, 4*nelorbh
                    dup_c(2*indpar + j) = 0.d0
                end do
            end if

            if (.not. symmagp .or. ipc .eq. 2) then
                indpar = indpar + 2*nelorbh
                if (ipc .eq. 1) then
                    do j = 1, nelorbh
                        dup_c(indpar + j) = j
                    end do
                    do j = nelorbh + 1, 2*nelorbh
                        dup_c(indpar + j) = 0.d0
                    end do
                else
                    do j = 1, nelorbh
                        dup_c(2*(indpar + j) - 1) = j
                        dup_c(2*(indpar + j)) = 0.d0
                    end do
                    do j = 2*nelorbh + 1, 4*nelorbh
                        dup_c(2*indpar + j) = 0.d0
                    end do
                end if
            end if

        end if ! endif i.le.bands
        indpar = indpar + 2*nelorbh
    end do
    !
    ! now filling the unpaired electrons in the last positions of dup_c
    !
    do i = nmol + 1, nmol + ndiff

        if (optocc .ne. 1) then
            ! In any case the last ones
            ind = nelocc + i - nmol - ndiff
        else
            ind = neldo + i - nmol
        end if
        if (ipc .eq. 1) then
            do j = 1, nelorbh
                dup_c(indpar + j) = j
            end do
            do j = nelorbh + 1, 2*nelorbh
                dup_c(indpar + j) = molecorb_old(j - nelorbh, ind)
            end do
        else
            do j = 1, nelorbh
                dup_c(2*(indpar + j) - 1) = j ! Re
                dup_c(2*(indpar + j)) = 0.d0 ! Im
            end do
            do j = nelorbh + 1, 2*nelorbh
                dup_c(2*(indpar + j) - 1) = molecorb_old(2*(j - nelorbh) - 1, ind) ! Re
                dup_c(2*(indpar + j)) = molecorb_old(2*(j - nelorbh), ind) ! Im
            end do
        end if
        indpar = indpar + 2*nelorbh

    end do ! enddo i=nmol+1,nmol+ndiff
    !
    ! if optimization of Z, required double dimension of dup_c matrix
    ! necessary for write_fort10()
    !
!    if(yeszagp) &
!            dup_c(ipc * iesup_cmol + 1:ipc * iesupr_cmol) = dup_c(1:ipc * iesup_cmol)

!    deallocate(ipsip)
!    allocate(ipsip(occ_c))
!    ipsip = ioccup_c
!    deallocate(ioccup_c)
!    occ_cmol = occ_c + nmoltot - molecular
!    allocate(ioccup_c(occ_cmol))
!    if(occ_cmol.gt.occ_c) then
!        ioccup_c(1:occ_c) = ipsip(1:occ_c)
!        ioccup_c(occ_c + 1:occ_cmol) = 1
!    else
!        ioccup_c(1:occ_cmol) = ipsip(1:occ_cmol)
!    endif

!    contraction = contraction + nmoltot - molecular

!    if(allocated(detmat_c)) deallocate(detmat_c)
!    if(allocated(nozero_c)) deallocate(nozero_c)
!    deallocate(jbradet, jbradetn)

!    if(molecular.gt.0) nelorb_c = nelorb_c - molecular
!    nelorb_c = nelorb_c + nmoltot
!    nelcol_c = nelorb_c + ndiff
!    allocate(detmat_c(ipc * nelorb_c * nelcol_c))
    detmat_c = 0.d0
    ind = nnozero_c - nmoltot
    if (.not. symmagp .or. ipc .eq. 2) ind = ind + nmol

!    if(molecular.eq.0) then
!        nnozero_c = nnozeroc_in + nmol + ndiff
!        iesswrc_in = iesswr
!        iesswr = iesswr + nmol + ndiff
!        molecsw = 0
!    else
!        if(symmagp.and.ipc.eq.1) then
!            nmol_in = molecular - ndiff
!            molecsw = molecular
!            nnozeroc_in = nnozeroc_in - molecular
!!            iesswrc_in = iesswr - molecsw
!        else
!            nmol_in = (molecular - ndiff) / 2
!            molecsw = molecular - nmol_in
!            nnozeroc_in = nnozeroc_in - molecsw
!            !nnozeroc_in=nnozeroc_in-molecular
!            iesswrc_in = iesswr - molecsw
!        endif
!        nnozero_c = nnozeroc_in + nmol + ndiff
!        iesswr = iesswrc_in + nmol + ndiff
!    endif

!    allocate(nozero_c(nnozero_c))
!    allocate(jbradet(nnozero_c), jbradetn(3 * nnozero_c))
!    nozero_c = 0
!    jbradet = 0
!    jbradetn = 0
!    if(allocated(ipsip)) deallocate(ipsip)
!    allocate(ipsip(3 * nnozero_c))
!    ipsip = 0
!    dimvect = nnozeroc_in + molecsw

!    allocate(indexo(dimvect))
!    indexo = 0
!    allocate(sortvect(dimvect + 1))
!    sortvect = 0
!    sortvect(1:dimvect) = abs(jbradetc_in(1:dimvect))
!    call dsortx(sortvect, 1, dimvect, indexo)
!
!    jbradet = 0
!    jbradetn = 0
!    icount1 = 1
!    ! individuation indnn
!    indnn = 0
!    indnn_in = 0
!
    !   if(rank.eq.0) then
    !      write(6,*) 'INFO:',iesswrc_in,iesswr,molecsw,nnozeroc_in
    !      write(6,*) 'INFO:',molecular,dimvect,nmol,nmol_in
    !      write(6,*) 'INFO:',iesup_cmol,iesupr_c,iesupr_cmol,iesupind
    !      do i = 1,nnozeroc_in
    !       write(6,*) '#',i,nozeroc_in(i)
    !      enddo
    !   endif

!    do i = 1, iesswrc_in + molecsw
!        ind = 0
!        do while(sortvect(icount1).lt.i.and.icount1.le.dimvect)
!            icount1 = icount1 + 1
!        enddo
!        ! check the very last
!        if(icount1.le.dimvect) then
!            do while(sortvect(icount1).eq.i.and.icount1.le.dimvect)
!                ind = ind + 1
!                icount1 = icount1 + 1
!            enddo
!            icount1 = icount1 - 1
!        endif
!
!        if(ind.eq.0) then
!            indnn = indnn + 1
!            ! load ipsip
!            ipsip(1) = jbradetnc_in(indnn)
!            ii = ipsip(1)
!            yest = .true.
!            do j = 1, -2 * ii, 2
!                ipsip(j + 1) = jbradetnc_in(j + indnn)
!                ipsip(j + 2) = jbradetnc_in(j + 1 + indnn)
!                if((abs(ipsip(j + 1)).gt.nelorbc_in - molecular&
!                        &.and.abs(ipsip(j + 1)).le.nelorbc_in)&
!                        &.or.(abs(ipsip(j + 2)).gt.nelorbc_in - molecular.and.&
!                                &abs(ipsip(j + 2)).le.nelorbc_in)) yest = .false.
!            enddo
!            if(yest) then
!                do j = 1, -2 * ii + 1
!                    if(ipsip(j).gt.nelorbc_in.and.j.ne.1) ipsip(j) = ipsip(j) + nelorb_c - nelorbc_in
!                    jbradetn(indnn_in + j) = ipsip(j)
!                enddo
!                indnn_in = indnn_in + 1 - 2 * ii
!            endif
!            indnn = indnn - 2 * ii
!        endif
!    enddo

!    deallocate(indexo, sortvect)

!    if(molecular.ne.0) then
!        inddiff = (indnn - indnn_in) / 3
!    else
!        inddiff = 1
!    endif

!    ind = 0
!
!    do i = 1, nnozeroc_in + molecsw
!
!        indy = (nozeroc_in(i) - 1) / nelorbc_in + 1
!        indx = nozeroc_in(i) - (indy - 1) * nelorbc_in
!
!        if(indx.le.nelorbc_in - molecular.and.indy.le.nelorbc_in - molecular) then
!            ind = ind + 1
!            jbradet(ind) = jbradetc_in(i)
!            nozero_c(ind) = (indy - 1) * nelorb_c + indx
!        elseif(indx.le.nelorbc_in - molecular.and.indy.gt.nelorbc_in) then
!            ind = ind + 1
!            jbradet(ind) = jbradetc_in(i)
!            nozero_c(ind) = (indy + nelorb_c - nelorbc_in - 1) * nelorb_c + indx
!        endif
!    enddo

    ! now sort jbradet
!    dimvect = ind
!    allocate(indexo(dimvect), indexn(dimvect))
!    indexo = 0
!    indexn = 0
!    allocate(sortvect(dimvect + 1)) ! Important ovedundant allocation
!    sortvect = 0
!    sortvect(1:dimvect) = abs(jbradet(1:dimvect))
!    indexn(1:dimvect) = jbradet(1:dimvect)  ! save the value...
!    jbradet = 0
!    inds = 0
!    i = 0
!    call dsortx(sortvect, 1, dimvect, indexo)
!    icount1 = 1
!
!    do while(inds.lt.iesswrc_in.and.i.lt.ind)
!        i = i + 1
!        indj = 0
!        ! Go to the first |jbradet|<i
!        do while(sortvect(icount1).lt.i.and.icount1.le.dimvect)
!            icount1 = icount1 + 1
!        enddo
!        ! If you are not at the end count how many |jbradet| =i
!        if(icount1.le.dimvect) then
!            do while(sortvect(icount1).eq.i.and.icount1.le.dimvect)
!                indj = indj + 1
!                if(indj.eq.1) indfirst = icount1
!                icount1 = icount1 + 1
!            enddo
!            ! At the end of the loop the condition |jbradet| /=i and icount1=icount1-1
!            icount1 = icount1 - 1
!            indlast = icount1
!            if(indj.ne.0) then
!                inds = inds + 1
!                do icount2 = indfirst, indlast
!                    j = indexo(icount2)
!                    if(indexn(j).gt.0) then
!                        jbradet(j) = inds
!                    elseif(indexn(j).lt.0) then
!                        jbradet(j) = -inds
!                    endif
!                enddo
!            endif
!        endif
!    enddo
!
!    deallocate(sortvect, indexo, indexn)

!    ind = nnozeroc_in
!    inds = iesswrc_in
    ! restored indnn here due to error in compiler ifort.10
!   indnn = indnn_in
    occmax = 0
    found1 = .false.
    do i = 1, min(nmol, bands)
        if (yeslsda .or. ipc .eq. 2) then
            cost = (occupations(i) + occupationdo(i))/2.d0
            if (cost .gt. 0.9d0) cost = 1.d0
        else
            cost = occupations(i)/2.d0
            if (cost .gt. 0.9d0) cost = 1.d0
        end if
        if (cost .eq. 1.d0) found1 = .true.
        if (cost .lt. 1.d0 .and. cost .gt. 1.d0/bands) then
            cost = 0.001d0
        elseif (cost .lt. 1.d0 .and. cost .le. 1.d0/bands) then
            cost = 0.d0
        end if
        if (cost .eq. 0 .and. .not. found1) cost = 0.001d0
        if (cost .ne. 0.d0) occmax = occmax + 1
    end do

    ii = 0
    found1 = .false.
    do i = 1, nmol

        if (symmagp .and. ipc .eq. 1) then
            indx = nelorb_c - nmoltot + i
            indy = nelorb_c - nmoltot + i
        else
            indx = nelorb_c - nmoltot + 2*i
            indy = indx - 1
        end if

        ind = ind + 1
!       nozero_c(ind) = nelorb_c * (indy - 1) + indx
        if (i .le. bands) then
            if (occupations(i) .ne. 0.d0) ii = ii + 1
            if (ii .le. countall .or. optocc .eq. 1) then

                if (yeslsda .or. ipc .eq. 2) then
                    cost = (occupations(i) + occupationdo(i))/2.d0
                    if (cost .gt. 0.9d0) cost = 1.d0
                else
                    cost = occupations(i)/2.d0
                    if (cost .gt. 0.9d0) cost = 1.d0
                end if
                ! degenerate shell, small coefficient
                if (cost .eq. 1.d0) found1 = .true.

                if (cost .lt. 1.d0 .and. cost .gt. 1.d0/bands) then
                    cost = 0.001d0
                elseif (cost .lt. 1.d0 .and. cost .le. 1.d0/bands) then
                    cost = 0.d0
                end if
                if (cost .eq. 0 .and. .not. found1) cost = 0.001d0
                !          For a single determinant wf does not make sense to put 0.001d0
                if (nmol .eq. neldo .or. (occmax .le. neldo .and. i .le. neldo)) cost = 1.d0
                if (ipc .eq. 1) then
                    detmat_c(nozero_c(ind)) = cost
                else
                    detmat_c(2*nozero_c(ind) - 1) = cost
                end if

            else
                if (ipc .eq. 1) then
                    detmat_c(nozero_c(ind)) = 0.d0
                else
                    detmat_c(2*nozero_c(ind) - 1) = 0.d0
                end if
            end if
        else
            if (ipc .eq. 1) then
                detmat_c(nozero_c(ind)) = 0.d0
            else
                detmat_c(2*nozero_c(ind) - 1) = 0.d0
            end if
        end if
!       indnn = indnn + 1
!       jbradetn(indnn) = -1
!       jbradetn(indnn + 1) = indx
!       jbradetn(indnn + 2) = indy
!       indnn = indnn + 2
    end do

    do i = 1, ndiff
        ind = ind + 1
        j = nelorb_c - ndiff + i ! the unpaired are always the last ones
!       nozero_c(ind) = nelorb_c * (nelorb_c + i - 1) + j
        if (ipc .eq. 1) then
            detmat_c(nozero_c(ind)) = 1.d0
        else
            detmat_c(2*nozero_c(ind) - 1) = 1.d0
        end if
!       indnn = indnn + 1
!       jbradetn(indnn) = -1
!       jbradetn(indnn + 1) = j
!       jbradetn(indnn + 2) = nelorb_c + i
!       indnn = indnn + 2
    end do

    do i = nnozero_c - nmol - ndiff + 1, nnozero_c
        if (ipc .eq. 1) then
            scale_c(i) = detmat_c(nozero_c(i))
        else
            scale_c(2*i - 1) = detmat_c(2*nozero_c(i) - 1)
            scale_c(2*i) = detmat_c(2*nozero_c(i))
        end if
    end do

!    deallocate(ipsip)
!    allocate(ipsip(iesup_c + iesupind))
!    iesupindmol = iesupind
!    ipsip = jbraiesup_sav
!    deallocate(jbraiesup_sav)
!    allocate(jbraiesup_sav(iesup_cmol + iesupindmol))
!!    jbraiesup_sav = 0
!    ! calculation of ind
!    ind = 0
!    do i = 1, iesupind
!        ind = ind + 1
!        ii = abs(ipsip(ind))
!        ind = ind + ii
!    enddo
!    jbraiesup_sav(1:ind) = ipsip(1:ind)
!    if(molopt.ge.2.and.contraction_in.ne.0) then
!        ! do not optimize the coefficients of the contracted orbitals
    ! But only if molopt>0 otherwise optimize also the coefficients
!        ind = 0
!        do i = 1, iesupind
!            ind = ind + 1
!            ii = jbraiesup_sav(ind)
!            do j = 1, abs(ii)
!                if(iesuptransb(jbraiesup_sav(ind + j)).eq.0.and.ii.gt.0)&
!                        jbraiesup_sav(ind) = -abs(jbraiesup_sav(ind))
!            enddo
!            ind = ind + abs(ii)
!        enddo
!    endif

!    indpar = iesup_c
!    nshell_c = nshellmol
!    iesupr_c = iesupr_cmol
!    iesup_c = iesup_cmol
!    occ_c = occ_cmol
!    iesupind = iesupindmol

    return

end subroutine update_fort10
