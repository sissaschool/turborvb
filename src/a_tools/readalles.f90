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

program er0read
    use allio
    implicit none
    integer ngenr, maxk, k, kt, i, j, maxj, ibin, nmis &
        , nbuf, kk, indpp, shellj, nbuf_dim &
        , dmat, jmat, ntot, nread, jas &
        , maxp, nnozero_read, nnozeroj_read, indswc, ii &
        , ntmp, ntmp1, iflagread, wfort10, drawpar, file_unit, fcol &
        , lcol, query, npipsz, ind, nnozero_fake, nnozeroj_fake, nwrite, ieskin_sav &
        , iesm_short, iesupr_short, iesinv_short, iesfree_short, nnozero_short, jj &
        , nnozeror, iesm_sav, iesd_sav, iesup_sav, iesfree_sav, iessw_sav, iesinv_sav &
        , iesking_sav, nnozero_ipc, iesupr_ipc, nnozero_eagp_short, ind_ghost

    real*8 wsk, tmp, tmp1
    real*8, dimension(:), allocatable :: wbuf, wkr, wbin, eskip
    real*8, dimension(:, :), allocatable :: ebuf, ek, ebin, ebin2
    logical, dimension(:), allocatable :: yespar

    character(64) command_file
    parameter(nbuf=1, nbuf_dim=0)
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'readalles'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added

    !real*8  ebuf(nbuf,npm),wbuf(nbuf),ek(0:nbuf,npm),wkr(0:nbuf),wsk&
    !,ebin(0:nbuf,npm),wbin(0:nbuf),ebin2(0:nbuf,npm),eskip(npm)&
    !,tbra,etry,tmp,tmp1,npow

    write (6, *) 'bin length, ibinit, write fort.10 (0/1), draw (0/1) ? '
    read (5, *) lbin, ibinit, wfort10, drawpar
    write (6, *) 'number of generations from standard input? (1  yes, 0 no) '
    read (5, *) query

    if (query .eq. 1) then
        write (6, *) 'ngen'
        read (5, *) ngen
    end if

    open (unit=10, file='fort.10', form='formatted', status='old')
    open (unit=23, file='parminimized.d', form='formatted', status='old')
    open (unit=11, file='fort.11', form='unformatted', status='old')
    open (unit=12, file='fort.12', form='unformatted', status='old')

    open (unit=20, file='Average_parameters.dat', form='formatted', status='unknown')
    open (unit=22, file='story.d', form='formatted', status='unknown')

    rewind (11)
    read (11) ngenr, nbra, npow, etry, nrest, nw, tbra, np, npbra, nnozero_read &
        , nnozeroj_read, yesmin

    !write(6,*) ' yesminjsz = ',yesminjsz

    !write(6,*) ' nnozero inside readalles =',nnozero
    !write(6,*) 'record 11 '
    !write(6,*) ngenr,nbra,npow,etry,nrest,nw,tbra,np,npbra,nnozero

    !           np=np-1
    np = np + 1

    rewind (23)
    read (23, *, end=133) iesinv, iesm, iesd, iesfree, iessw, iesup, ieskin, jj, isfix &
        , ieser, ieser, iesking

    iesinv_sav = iesinv
    iesm_sav = iesm
    iesd_sav = iesd
    iesfree_sav = iesfree
    iessw_sav = iessw
    iesup_sav = iesup
    ieskin_sav = ieskin
    iesking_sav = iesking

    !mod by Kosuke Nakano on 29th May 2019
    !this arises inconsistency between drawpar=1 and 0
    !if(drawpar.ge.1) then
    write (6, *) 'max number of ind par  for each part of the wf  '
    read (5, *) maxp
    !else
    !maxp=1
    !endif
    !mod by Kosuke Nakano on 29th May 2019

133 continue
    if (ieskin + iesking .gt. 0) then
        ieskin = 1 ! default value
    else
        ieskin = 0
    end if

    call read_fort10_fast ! read and correct parameters if they are not accurate

    iesm = iesm_sav
    iesd = iesd_sav
    iesfree = iesfree_sav
    iessw = iessw_sav
    iesup = iesup_sav
    iesinv = iesinv_sav
    ieskin = ieskin_sav
    iesking = iesking_sav

    nnozero_fake = nnozero + nnozero_eagp
    nnozeroj_fake = nnozeroj
    nnozero = nnozero_read + nnozero_eagp ! restore the value read in fort.11
    nnozeroj = nnozeroj_read

    iesinv = abs(iesinv)
    iesfree = abs(iesfree)
    iessw = abs(iessw)

    ieskint = ieskin + iesking

    if (yesmin .ne. 0) then
        iesup_read = iesup_atom*ipc
        iesupr = iesup_atom
    else
        iesup_read = iesupr*ipc
    end if

    if (ieskint .gt. ieskinr) then
        ntot = ipc*nnozero + iesinv + iesfree + iesm + iesup_read + 3*nion + iesd + 3

        iessw = 1 ! read in any event the matrix
        nnozeror = nnozero
    elseif (ieskint .gt. 0) then
        ntot = ipc*nnozero + iesinv + iesfree + iesm + iesup_read + 3*nion + iesd
        iessw = 1 ! read in any event the matrix
        nnozeror = nnozero
    else
        ntot = ipc*nnozero + iesinv + iesfree + iesm + iesup_read + iesd
        ! minimize read
        if (iesup .eq. 0) ntot = ntot - ipc*iesupr
        nnozeror = nnozero
        if (iessw .eq. 0) then
            ntot = ntot - ipc*nnozero
            nnozeror = 0
        end if
    end if

    indswc = iesinv + iesfree + iesm + iesd

    write (6, *) ' # words read from unit 11 ', ntot
    write (6, *) ' maxp= ', maxp
    rewind (12)

    !           stop

    !real*8  ebuf(nbuf,npm),wbuf(nbuf),ek(0:nbuf,npm),wkr(0:nbuf),wsk&
    !,ebin(0:nbuf,npm),wbin(0:nbuf),ebin2(0:nbuf,npm),eskip(npm)&
    !,tbra,etry,tmp,tmp1,npow
    allocate (wbuf(nbuf), wkr(0:nbuf_dim), wbin(0:nbuf_dim), eskip(ntot), psip(ntot))
    allocate (ebuf(nbuf, ntot), ek(0:nbuf_dim, ntot), ebin(0:nbuf_dim, ntot), ebin2(0:nbuf_dim, ntot))
    allocate (yespar(ntot), ipsip(ntot))

    yespar = .true.
    iesm_short = iesm
    iesupr_short = iesupr
    iesinv_short = iesinv
    iesfree_short = iesfree
    nnozero_short = nnozeror

    np = ntot
    irej = 0

    do i = 0, nbuf_dim
        wkr(i) = 0.d0
        wbin(i) = 0.d0
        do k = 1, np
            ebin2(i, k) = 0.d0
            ebin(i, k) = 0.d0
            ek(i, k) = 0.d0
        end do
    end do

    if (query .eq. 0) then
        ngen = 0
        do while (ngen .ge. 0)
            read (12, end=500)
            ngen = ngen + 1
            if (mod(ngen, 1) .eq. 0) write (6, *) ' record read =', ngen
        end do
500     continue
        rewind (12)
        if (ngen .ne. ngenr .and. ngenr .ne. 0) then
            write (6, *) ' partial file read , ngen =', ngen
            !     ngen=ngen-1
        end if
    end if

    maxk = 0
    ng = ngen
    nbin = ngen/lbin

    rewind (12)
    do i = 1, ng - 1
        read (12)
    end do
    read (12) wkr(0), tmp, (ek(0, k), k=1, ntot)

    !      write(*,*) wkr(0),tmp,(ek(0,k),k=1,ntot)

    !     write(6,*) ' final scale read '
    !     do ii=1,nnozero
    !     write(6,*) ii,ek(0,indswc+ii)
    !     enddo

    nnozero_ipc = ipc*nnozeror
    iesupr_ipc = ipc*iesupr

    call defyespar(maxk, nbuf_dim, ek, wkr, psip, iesinv, iesm, iesd, iesfree&
            &, nnozero_ipc, iesupr_ipc, nion, ieskint, ieskinr, ind, yespar, ipsip, iesm_short&
            &, iesupr_short, iesfree_short, iesinv_short, nnozero_short, maxp, iesup, iessw&
            &, ieskin, nnozero_eagp, nnozero_eagp_short)

    rewind (12)

    write (20, *) ' Number  of generations  read =', ngen
    write (20, *) ' # of bins =', nbin

    ibin = 0
    icount = 0
    !          check

    maxj = 1

    do kt = 1, ng
        if (ibin .ge. ibinit - 1) then
            do j = 1, maxj
                read (12) wbuf(j), tmp, (eskip(k), k=1, ntot)
                !      write(6,*) ' scale read ',kt
                !      do ii=1,nnozero
                !      write(6,*) ii,eskip(indswc+ii)
                !      enddo
                if (wbuf(j) .eq. 0.d0) irej = irej + 1
                do k = 1, ntot
                    ebuf(j, k) = eskip(k)
                end do
            end do
        else
            do j = 1, maxj
                read (12)
                wbuf(j) = 0.d0
                ebuf(j, 1:ntot) = 0.d0
            end do
        end if

        do j = 1, maxj
            icount = icount + 1
            do kk = 1, ntot
                ek(0, kk) = ek(0, kk) + ebuf(j, kk)*wbuf(j)
            end do
            wkr(0) = wkr(0) + wbuf(j)

            if (mod(icount, lbin) .eq. 0) then

                ibin = ibin + 1

                if (ibin .ge. ibinit) then
                    if (wkr(0) .ne. 0.d0) then
                        do k = 0, maxk
                            wbin(k) = wbin(k) + wkr(k)
                            do kk = 1, ntot
                                ebin(k, kk) = ebin(k, kk) + ek(k, kk)
                                ebin2(k, kk) = ebin2(k, kk) + ek(k, kk)**2/wkr(k)
                            end do
                        end do

                        nnozero_ipc = ipc*nnozeror
                        iesupr_ipc = ipc*iesupr

                        call sortpsip(maxk, nbuf_dim, ek, wkr, psip, iesinv, iesm, iesd, iesfree&
                                &, nnozero_ipc, iesupr_ipc, nion, ieskint, ieskinr, nwrite, yespar, iesup&
                                &, iessw, ieskint, nnozero_eagp)

                        write (22, 123) ibin, (psip(kk), kk=1, nwrite)

                    end if
                end if
                do k = 0, maxk
                    do kk = 1, np
                        ek(k, kk) = 0.d0
                    end do
                    wkr(k) = 0.d0
                    !                wkr(k,2)=0.d0
                end do
            end if
        end do
    end do
    !            calculation error bars

    if (ibinit .gt. 0) then
        nmis = ibin - ibinit + 1
    else
        nmis = ibin
    end if

    write (6, *) ' number of measures done =', nmis
    write (6, *) ' Rejected measures =', irej
    write (6, *) ' Rejection ratio =', dble(irej)/icount

    do i = 0, maxk
        do kk = 1, ntot
            ebin(i, kk) = ebin(i, kk)/wbin(i)
            ebin2(i, kk) = dsqrt(dabs(ebin2(i, kk)/wbin(i) - ebin(i, kk)**2))
            if (nmis .gt. 1) ebin2(i, kk) = ebin2(i, kk)/dsqrt(dfloat(nmis - 1))
        end do
    end do
    write (20, *) ' Independent bins ', nmis, 'of length ', lbin
    i = 0
    if (ieskint .gt. ieskinr) then
        write (20, *) ' Atomic position '
        do kk = ntot - 3*nion - 2, ntot - 3
            write (20, 123) kk - ntot + 3*nion + 3, ebin(i, kk), ebin2(i, kk)
        end do
        write (20, *) ' Cell derivatives '
        do kk = ntot - 2, ntot
            write (20, 123) kk - ntot + 3, ebin(i, kk), ebin2(i, kk)
        end do
        write (20, *) ' wf  parameters '
        do kk = 1, ntot - 3*nion - 3
            write (20, 123) kk, ebin(i, kk), ebin2(i, kk)
        end do
    elseif (ieskint .gt. 0) then
        write (20, *) ' Atomic position '
        do kk = ntot - 3*nion + 1, ntot
            write (20, 123) kk - ntot + 3*nion, ebin(i, kk), ebin2(i, kk)
        end do
        write (20, *) ' wf  parameters '
        do kk = 1, ntot - 3*nion
            write (20, 123) kk, ebin(i, kk), ebin2(i, kk)
        end do
    else
        do kk = 1, ntot
            write (20, 123) kk, ebin(i, kk), ebin2(i, kk)
        end do
    end if

    write (20, *) ' with no error bars '
    write (6, *)
    if (wfort10 .eq. 1) then
        close (10)
        open (unit=10, file='fort.10', form='formatted', status='old', position='append')
        !   kk=0
        !   do while(kk.ge.0)
        !   kk=kk+1
        !   read(10,*,END=1000)
        !   enddo
        !
        !   stop

        write (10, *) '# new parameters'
        do i = 0, maxk
#ifdef __KCOMP
            write (10, *) (ebin(i, kk), kk=1, ntot)
#else
            write (10, 125) (ebin(i, kk), kk=1, ntot)
#endif
        end do
    end if
#ifdef __KCOMP
123 format(i15, 32767e24.15)
#else
123 format(i15, 1000000000e24.15)
125 format(1000000000e24.15)
#endif
    close (11)
    close (12)
    close (20)
    close (10)

    if (drawpar .ge. 1) then

        lcol = ntot
        indpp = 0
        write (*, *) ' Draw the story of parameters'
        file_unit = 30
        command_file = 'commands.story'
        open (unit=file_unit, file=command_file, form='formatted', status='unknown')
        write (file_unit, '(a)') 'set title "Fort.22 Par"'
        write (file_unit, '(a)') 'set xlabel "x"'
        do i = 1, iesd
            write (file_unit, '(a,i15,a)') 'set ylabel  "Jastrow Two Body"'
            write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', i + 1
            write (file_unit, '(a)') 'pause -1'
            indpp = indpp + 1
        end do

        if (ieskint .gt. ieskinr) then
            ind = 0
            do i = iesd + 1, iesd + 3
                ind = ind + 1
                indpp = indpp + 1
                if (ind .eq. 1) write (file_unit, '(a,i15,a)') 'set ylabel  "rs"'
                if (ind .eq. 2) write (file_unit, '(a,i15,a)') 'set ylabel  "b/a"'
                if (ind .eq. 3) write (file_unit, '(a,i15,a)') 'set ylabel  "c/a"'
                write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', indpp + 1
                write (file_unit, '(a)') 'pause -1'
            end do
        end if
        ind = indpp
        do i = ind + 1, ind + iesm_short
            indpp = indpp + 1
            write (file_unit, '(a,i15,a)') 'set ylabel "Three Body Z Parameter ', i - ind, '"'
            write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', indpp + 1
            write (file_unit, '(a)') 'pause -1'
        end do

        ind = indpp
        do i = ind + 1, ind + iesupr_short
            indpp = indpp + 1
            write (file_unit, '(a,i15,a)') 'set ylabel  "Orbital Z Parameter ', i - ind, '"'
            write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', indpp + 1
            write (file_unit, '(a)') 'pause -1'
        end do
        if (ieskint .ne. 0) then
            ind = indpp
            do i = ind + 1, ind + 3*nion
                indpp = indpp + 1
                write (file_unit, '(a,i15,a)') 'set ylabel  "Nuclear Position ', i - ind, '"'
                write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', indpp + 1
                write (file_unit, '(a)') 'pause -1'
            end do
        end if
        ind = indpp
        do i = ind + 1, ind + iesfree_short
            indpp = indpp + 1
            write (file_unit, '(a,i15,a)') 'set ylabel "Three Body Parameter Matrix ', i - ind, '"'
            write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', indpp + 1
            write (file_unit, '(a)') 'pause -1'
        end do
        ind = indpp
        do i = ind + 1, ind + iesinv_short
            indpp = indpp + 1
            write (file_unit, '(a,i15,a)') 'set ylabel "Three Body-Sz Parameter Matrix', i - ind, '"'
            write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', indpp + 1
            write (file_unit, '(a)') 'pause -1'
        end do
        ind = indpp
        ind_ghost = ind
        do i = ind + 1, ind + nnozero_short
            indpp = indpp + 1
            if (i - ind .le. nnozero_eagp_short) then
                write (file_unit, '(a,i15,a)') 'set ylabel  "Ghost coeff   ', i - ind, '"'
                ind_ghost = i
            else
                write (file_unit, '(a,i15,a)') 'set ylabel  "Det Matrix ', i - ind_ghost, '"'
            end if
            write (file_unit, '(a,i15)') 'plot "story.d"  using 1:', indpp + 1
            write (file_unit, '(a)') 'pause -1'
        end do

        write (file_unit, '(a)') 'q'
        close (unit=file_unit)
        close (unit=22)
        if (drawpar .eq. 1) call run_gnuplot(command_file)
    end if

    stop
end

subroutine run_gnuplot(command_line)
    implicit none
    character(50) command
    character(64) command_line

    write (command, *) 'gnuplot commands.story'
    call system(trim(command))

    return
end
subroutine sortpsip(maxk, nbuf, ek, wkr, psip, nfreesz, diesm, djas, jorb, nnozero &
                    , iesupr, nion, dieskin, ieskinr, ind, yespar, iesup, iessw, ieskin, nnozero_eagp)
    implicit none
    integer maxk, nbuf, nfreesz, diesm, djas, jorb, nnozero, iesupr, nion &
        , dieskin, ieskinr, i, ind, iesup, iessw, ieskin, nnozero_eagp
    !       this subroutine put the long output nnozero nfreesz at the end
    real*8 ek(0:nbuf, *), wkr(0:nbuf), psip(*)
    logical yespar(*)
    ind = 0
    do i = nfreesz + diesm + 1, nfreesz + diesm + djas
        ind = ind + 1
        psip(ind) = ek(maxk, i)/wkr(maxk)
    end do
    if (dieskin .gt. ieskinr .and. ieskin .ne. 0) then
        do i = nfreesz + diesm + djas + jorb + nnozero + iesupr + 3*nion + 1 &
               , nfreesz + diesm + djas + jorb + nnozero + iesupr + 3*nion + 3
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
    end if

    do i = nfreesz + 1, nfreesz + diesm
        if (yespar(i)) then
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end if
    end do
    if (iesup .gt. 0) then
        do i = nfreesz + 1 + diesm + djas + jorb + nnozero &
               , nfreesz + diesm + djas + jorb + nnozero + iesupr
            if (yespar(i)) then
                ind = ind + 1
                psip(ind) = ek(maxk, i)/wkr(maxk)
            end if
        end do
    end if
    if (dieskin .ne. 0 .and. ieskin .gt. 0) then
        do i = nfreesz + 1 + diesm + djas + jorb + nnozero + iesupr &
               , nfreesz + diesm + djas + jorb + nnozero + iesupr + 3*nion
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
    end if
    do i = nfreesz + 1 + diesm + djas, nfreesz + diesm + djas + jorb
        if (yespar(i)) then
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end if
    end do
    do i = 1, nfreesz
        if (yespar(i)) then
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end if
    end do
    if (iessw .ne. 0) then
        do i = nfreesz + diesm + djas + jorb + nnozero - nnozero_eagp + 1 &
               , nfreesz + diesm + djas + jorb + nnozero
            if (yespar(i)) then
                ind = ind + 1
                psip(ind) = ek(maxk, i)/wkr(maxk)
            end if
        end do
        do i = nfreesz + 1 + diesm + djas + jorb &
               , nfreesz + diesm + djas + jorb + nnozero - nnozero_eagp
            if (yespar(i)) then
                ind = ind + 1
                psip(ind) = ek(maxk, i)/wkr(maxk)
            end if
        end do
    end if
    return
end

subroutine defyespar(maxk, nbuf, ek, wkr, psip, nfreesz, diesm, djas, jorb, nnozero&
        &, iesupr, nion, dieskin, ieskinr, ind, yespar, ipsip, diesm_short, iesupr_short&
        &, jorb_short, nfreesz_short, nnozero_short, maxp, iesup, iessw, ieskin&
        &, nnozero_eagp, nnozero_eagp_short)
    implicit none
    integer maxk, nbuf, nfreesz, diesm, djas, jorb, nnozero, iesupr, nion &
            , dieskin, ieskinr, i, ind, inds, ipsip(*), ind_short, indp &
            , diesm_short, iesupr_short, jorb_short, nfreesz_short, nnozero_short, maxp&
            &, iesup, iessw, ieskin, nnozero_eagp, nnozero_eagp_short
    !       this subroutine put the long output nnozero nfreesz at the end
    real*8 ek(0:nbuf, *), wkr(0:nbuf), psip(*), epsmach, epsmachd
    logical yespar(*)

    epsmach = 1d-12
    epsmachd = 1d-8
    ind = 0
    do i = nfreesz + diesm + 1, nfreesz + diesm + djas
        ind = ind + 1
        yespar(i) = .true.
        psip(ind) = ek(maxk, i)/wkr(maxk)
    end do
    if (dieskin .gt. ieskinr .and. ieskin .gt. 0) then
        do i = nfreesz + diesm + djas + jorb + nnozero + iesupr + 3*nion + 1 &
               , nfreesz + diesm + djas + jorb + nnozero + iesupr + 3*nion + 3
            yespar(i) = .true.
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
    end if

    inds = ind
    diesm_short = 0
    if (diesm .ne. 0) then
        do i = nfreesz + 1, nfreesz + diesm
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
        call dsortx(psip(inds + 1), 1, diesm, ipsip)
        ind = inds
        indp = nfreesz
        do i = nfreesz + 1, nfreesz + diesm
            ind = ind + 1
            if (ind .eq. inds + 1) then
                yespar(indp + ipsip(1)) = .true.
                diesm_short = diesm_short + 1
            elseif (abs(psip(ind) - psip(ind - 1)) .gt. epsmach .and. nint(psip(ind)) .ne. psip(ind)) then
                yespar(indp + ipsip(ind - inds)) = .true.
                diesm_short = diesm_short + 1
            else
                yespar(indp + ipsip(ind - inds)) = .false.
            end if
        end do
        if (diesm_short .gt. maxp) then
            ind = inds
            diesm_short = 0
            do i = nfreesz + 1, nfreesz + diesm
                ind = ind + 1
                if (yespar(i) .and. diesm_short .lt. maxp) then
                    diesm_short = diesm_short + 1
                else
                    yespar(i) = .false.
                end if
            end do
        end if

    end if

    inds = ind
    iesupr_short = 0
    if (iesupr .ne. 0 .and. iesup .ne. 0) then

        do i = nfreesz + 1 + diesm + djas + jorb + nnozero, nfreesz + diesm + djas + jorb + nnozero + iesupr
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
        call dsortx(psip(inds + 1), 1, iesupr, ipsip)

        ind = inds
        indp = nfreesz + diesm + djas + jorb + nnozero
        do i = nfreesz + 1 + diesm + djas + jorb + nnozero, nfreesz + diesm + djas + jorb + nnozero + iesupr
            ind = ind + 1
            if (ind .eq. inds + 1) then
                yespar(indp + ipsip(1)) = .true.
                iesupr_short = iesupr_short + 1
            elseif (abs(psip(ind) - psip(ind - 1)) .gt. epsmach .and. nint(psip(ind)) .ne. psip(ind)) then
                yespar(indp + ipsip(ind - inds)) = .true.
                iesupr_short = iesupr_short + 1
            else
                yespar(indp + ipsip(ind - inds)) = .false.
            end if
        end do

        if (iesupr_short .gt. maxp) then
            ind = inds
            iesupr_short = 0
            do i = nfreesz + 1 + diesm + djas + jorb + nnozero, nfreesz + diesm + djas + jorb + nnozero + iesupr
                ind = ind + 1
                if (yespar(i) .and. iesupr_short .lt. maxp) then
                    iesupr_short = iesupr_short + 1
                else
                    yespar(i) = .false.
                end if
            end do
        end if

    end if

    if (dieskin .ne. 0 .and. ieskin .gt. 0) then
        do i = nfreesz + 1 + diesm + djas + jorb + nnozero + iesupr &
               , nfreesz + diesm + djas + jorb + nnozero + iesupr + 3*nion
            ind = ind + 1
            yespar(i) = .true.
        end do
    end if
    inds = ind
    jorb_short = 0
    if (jorb .gt. 0) then
        do i = nfreesz + 1 + diesm + djas, nfreesz + diesm + djas + jorb
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
        call dsortx(psip(inds + 1), 1, jorb, ipsip)
        ind = inds
        indp = nfreesz + diesm + djas
        do i = nfreesz + 1 + diesm + djas, nfreesz + diesm + djas + jorb
            ind = ind + 1
            if (ind .eq. inds + 1) then
                yespar(indp + ipsip(1)) = .true.
                jorb_short = jorb_short + 1
            elseif (abs(psip(ind) - psip(ind - 1)) .gt. epsmachd .and. nint(psip(ind)) .ne. psip(ind)) then
                yespar(indp + ipsip(ind - inds)) = .true.
                jorb_short = jorb_short + 1
            else
                yespar(indp + ipsip(ind - inds)) = .false.
            end if
        end do
        if (jorb_short .gt. maxp) then
            ind = inds
            jorb_short = 0
            do i = nfreesz + 1 + diesm + djas, nfreesz + diesm + djas + jorb
                ind = ind + 1
                if (yespar(i) .and. jorb_short .lt. maxp) then
                    jorb_short = jorb_short + 1
                else
                    yespar(i) = .false.
                end if
            end do
        end if
    end if

    inds = ind
    nfreesz_short = 0
    if (nfreesz .gt. 0) then
        do i = 1, nfreesz
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
        call dsortx(psip(inds + 1), 1, nfreesz, ipsip)
        ind = inds
        indp = 0
        do i = 1, nfreesz
            ind = ind + 1
            if (ind .eq. inds + 1) then
                yespar(indp + ipsip(1)) = .true.
                nfreesz_short = nfreesz_short + 1
            elseif (abs(psip(ind) - psip(ind - 1)) .gt. epsmachd .and. nint(psip(ind)) .ne. psip(ind)) then
                yespar(indp + ipsip(ind - inds)) = .true.
                nfreesz_short = nfreesz_short + 1
            else
                yespar(indp + ipsip(ind - inds)) = .false.
            end if
        end do
        if (nfreesz_short .gt. maxp) then
            ind = inds
            nfreesz_short = 0
            do i = 1, nfreesz
                ind = ind + 1
                if (yespar(i) .and. nfreesz_short .lt. maxp) then
                    nfreesz_short = nfreesz_short + 1
                else
                    yespar(i) = .false.
                end if
            end do
        end if
    end if

    nnozero_short = 0
    if (nnozero .gt. 0 .and. iessw .ne. 0) then
        inds = ind
        do i = nfreesz + 1 + diesm + djas + jorb, nfreesz + diesm + djas + jorb + nnozero
            ind = ind + 1
            psip(ind) = ek(maxk, i)/wkr(maxk)
        end do
        call dsortx(psip(inds + 1), 1, nnozero, ipsip)

        ind = inds
        indp = nfreesz + diesm + djas + jorb
        do i = nfreesz + 1 + diesm + djas + jorb, nfreesz + diesm + djas + jorb + nnozero
            ind = ind + 1
            if (ind .eq. inds + 1) then
                yespar(indp + ipsip(1)) = .true.
                nnozero_short = nnozero_short + 1
            elseif (abs(psip(ind) - psip(ind - 1)) .gt. epsmachd .and. nint(psip(ind)) .ne. psip(ind)) then
                yespar(indp + ipsip(ind - inds)) = .true.
                nnozero_short = nnozero_short + 1
            else
                yespar(indp + ipsip(ind - inds)) = .false.
            end if
        end do
        ! if(nnozero_short.gt.maxp) then
        ind = inds
        nnozero_short = 0
        nnozero_eagp_short = 0
        do i = nfreesz + diesm + djas + jorb + nnozero - nnozero_eagp + 1, nfreesz + diesm + djas + jorb + nnozero
            ind = ind + 1
            if (yespar(i) .and. nnozero_short .lt. maxp) then
                nnozero_eagp_short = nnozero_eagp_short + 1
                nnozero_short = nnozero_short + 1
            else
                yespar(i) = .false.
            end if
        end do
        do i = nfreesz + 1 + diesm + djas + jorb, nfreesz + diesm + djas + jorb + nnozero - nnozero_eagp
            ind = ind + 1
            if (yespar(i) .and. nnozero_short .lt. maxp) then
                nnozero_short = nnozero_short + 1
            else
                yespar(i) = .false.
            end if
        end do
        ! endif
    end if
    return
end

