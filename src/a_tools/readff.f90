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

program erread
    implicit none
    integer(4) ngen, maxk, k, kt, i, j, ng, nrest, nbin, npm         &
            &, iskip, ibinit, icount, ibin, nmis, lbin, nbuf, nbufp, iskipr, countread&
            &, ntime, countmax, iseed
    real(8) wsk

    real(8), dimension(:), allocatable :: ebuf, wbuf, ek, wk          &
            &, wbufm, eskip, wtot
    real*8, dimension(:), allocatable :: eav, wav
    real*8, dimension(:, :), allocatable :: ekw, wkw
    real*8, dimension(:), allocatable :: gapjack, gap, gapav, gaps&
            &, gapt, gapst, error, errors

    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'readff'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added

    write (6, *) ' max k corrections, bin length, ibinit,iskip  ?'
    read (5, *) maxk, lbin, ibinit, iskipr

    if (iskipr .ge. 0) then
        open (unit=12, file='fort.12', form='unformatted', status='old')
        iskip = iskipr
    else
        iskip = -iskipr - 1
        open (unit=12, file='fort.12.fn', form='unformatted', status='old')
    end if

    open (unit=20, file='fort.20', form='formatted', status='unknown')
    open (unit=21, file='fort.21', form='formatted', status='unknown')
    !           rewind (11)
    !           read(11) ngen
    rewind (12)

    nbuf = max(maxk, 1)

    !      if(nbuf.ne.1) then
    !      nbuf=10*((nbuf-1)/10+1)
    !      endif

    npm = iskip + 1

    !      write(6,*) ' nbuf/npm chosen =',nbuf,npm

    allocate (ebuf(nbuf), wbuf(nbuf), ek(0:nbuf), wk(0:nbuf)             &
            &, wbufm(nbuf), eskip(npm), wtot(nbuf))

    iseed = 123754123
    call rand_init(iseed)

    !          Initialization
    wbufm(:) = 1.d0
    wbuf(:) = 1.d0

    ebuf = 0.d0
    eskip = 0.d0
    wtot = 0.d0

    nbufp = nbuf + 1
    do i = 0, nbuf
        ek(i) = 0.d0
        wk(i) = 0.d0
    end do
    i = 0
    do while (i .ge. 0)
        read (12, end=123, err=123)
        i = i + 1
    end do
123 countmax = (i/lbin)*lbin
    rewind (12)

    nbin = countmax/lbin - ibinit + 1
    ntime = maxk + 1

    if (nbin .le. 0 .or. maxk .lt. 0) stop

    allocate (ekw(ntime, nbin), wkw(ntime, nbin))
    allocate (eav(ntime), wav(ntime), gapjack(ntime), gap(ntime), errors(ntime)&
            &, gaps(ntime), gapav(ntime), gapt(ntime), gapst(ntime), error(ntime))

    icount = 0
    ibin = 0
    !          check
    kt = 0
    countread = 1
    do while (kt .ge. 0 .and. icount .lt. countmax)
        kt = kt + 1
        do i = 1, nbuf
            wbufm(i) = wbuf(i)
        end do
        !             call dcopy(nbuf,wbuf,1,wbufm,1)

        !             write(6,*) ' buf number ',kt
        countread = 0
        do j = 1, nbuf
            if (iskip .ne. 0) then
                read (12, end=100, err=100) wbuf(j), wtot(j), (eskip(k), k=1, iskip), ebuf(j)
            else
                read (12, end=100, err=100) wbuf(j), wtot(j), ebuf(j)
            end if
            countread = countread + 1
        end do
100     continue

        if (kt .eq. 1 .and. countread .ne. 0) then
            do j = 1, countread
                icount = icount + 1
                do k = 0, maxk
                    if (k .eq. 0) then
                        wsk = wtot(j)
                    else
                        if (j - k .gt. 0) then
                            wsk = wsk*wbuf(j - k)
                        else
                            wsk = wsk*wbufm(j - k + nbuf)
                        end if
                    end if
                    ek(k) = ek(k) + ebuf(j)*wsk
                    wk(k) = wk(k) + wsk
                end do

                !        write(6,*) ' countread wsk ',icount,wsk/wtot(j),wtot(j),ebuf(j)

                if (mod(icount, lbin) .eq. 0) then
                    ibin = ibin + 1
                    if (ibin .ge. ibinit) then
                        ekw(1:ntime, ibin - ibinit + 1) = ek(0:maxk)
                        wkw(1:ntime, ibin - ibinit + 1) = wk(0:maxk)
                        if (wk(maxk) .ne. 0.d0) then
                            write (21, *) ek(maxk)/wk(maxk), wk(maxk), ibin

                        else
                            write (21, *) 0., wk(maxk), ibin
                        end if

                    end if
                    do k = 0, maxk
                        ek(k) = 0.d0
                        wk(k) = 0.d0
                    end do
                end if

                ! end j
            end do

        elseif (countread .ne. 0) then

            do j = 1, countread
                icount = icount + 1
                do k = 0, maxk
                    if (k .eq. 0) then
                        wsk = wtot(j)
                    else
                        if (j - k .ge. 1) then
                            wsk = wsk*wbuf(j - k)
                        else
                            wsk = wsk*wbufm(j - k + nbuf)
                        end if
                    end if
                    ek(k) = ek(k) + ebuf(j)*wsk
                    wk(k) = wk(k) + wsk
                    ! end do k
                end do

                !        write(6,*) ' countread wsk ',icount,wsk/wtot(j),wtot(j),ebuf(j)

                if (mod(icount, lbin) .eq. 0) then
                    ibin = ibin + 1
                    if (ibin .ge. ibinit) then
                        ekw(1:ntime, ibin - ibinit + 1) = ek(0:maxk)
                        wkw(1:ntime, ibin - ibinit + 1) = wk(0:maxk)
                        if (wk(maxk) .ne. 0.d0) then
                            write (21, *) ek(maxk)/wk(maxk), wk(maxk), ibin
                        else
                            write (21, *) 0., wk(maxk), ibin
                        end if
                    end if
                    do k = 0, maxk
                        ek(k) = 0.d0
                        wk(k) = 0.d0
                    end do
                end if
            end do
        end if
    end do
    !            calculation error bars
    nmis = ibin - ibinit + 1
    call boot
    write (6, *) ' number of measures done =', nmis
    write (20, *) ' Independent bins ', nmis, 'of length ', lbin
    write (20, *)
    write (20, *) ' Energy , error, derror, # of bias  correcting factor '
    do i = 0, maxk
        write (20, '(3e20.10,I6)') gapjack(i + 1), gaps(i + 1), error(i + 1), i
    end do
    deallocate (ebuf, wbuf, ek, wk, wbufm, eskip, wtot)
    deallocate (ekw, wkw)
    deallocate (eav, wav, gapjack, gap, gaps, gapav, gapt, gapst, error)
    stop
contains
    subroutine boot
        implicit none
        integer nmis, kmain, k
        real*8, external :: drand1
        nmis = 1000
        eav = 0.d0
        wav = 0.d0
        do i = 1, nbin
            eav(:) = eav(:) + ekw(:, i)
            wav(:) = wav(:) + wkw(:, i)
        end do
        eav(:) = eav(:)/wav(:)

        gapjack(:) = eav(:)

        gaps = 0.d0
        gapav = 0.d0
        error = 0.d0
        errors = 0.d0
        do kmain = 1, nmis
            gapt = 0.d0
            gapst = 0.d0
            wav = 0.d0
            do k = 1, nbin
                i = drand1()*nbin + 1
                !               do j=1,ntime
                !               gapt(j)=(eav(j)*wav(j)-ekw(j,i))/(wav(j)-wkw(j,i))
                !               enddo
                do j = 1, ntime
                    wav(j) = wav(j) + wkw(j, i)
                    gapst(j) = gapst(j) + ekw(j, i)**2/wkw(j, i)
                    gapt(j) = gapt(j) + ekw(j, i)
                end do
            end do
            gapt(:) = gapt(:)/wav(:)
            gapst(:) = gapst(:)/wav(:)

            gapav(:) = gapav(:) + gapt(:)
            gaps(:) = gaps(:) + gapt(:)**2

            gap(:) = dsqrt(max(gapst(:) - gapt(:)**2, 0.d0)/(nbin - 1))

            error(:) = error(:) + gap(:)
            errors(:) = errors(:) + gap(:)**2
        end do
        gapav = gapav/nmis
        gaps = gaps/nmis
        error = error/nmis
        errors = errors/nmis

        gaps(:) = dsqrt(max(gaps(:) - gapav(:)**2, 0.d0))

        error(:) = dsqrt(max(errors(:) - error(:)**2, 0.d0))

    end subroutine boot

end
