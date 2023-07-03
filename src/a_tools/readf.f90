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
    integer(4) k, kt, i, iread, j, ng, nrest, nbin, npm         &
            &, iskip, icount, ibin, nmis, nbuf, nbufp, iskipr, countread&
            &, ntime, countmax, iflag, iflag_ave, rest, ikp, krep
    real(8) wsk, count_mes

    integer(4) :: maxk ! number of correcting factors
    integer(4) :: ibinit, ibinitr ! initial bin for averages
    integer(4) :: lbin ! bin length
    integer :: kp_read ! current replica/k-point index
    integer :: kave_read ! current replica/k-point index
    integer :: nkps ! total number of k-points/replicas
    real(8) :: wkps, kps(3) ! auxiliary variables

    real(8), dimension(:), allocatable :: wbuf, wk, wbufm, eskip, wtot
    real(8), dimension(:, :), allocatable :: ebuf, ek, eav
    real*8, dimension(:), allocatable :: wav, cost_mes
    real*8, dimension(:, :, :), allocatable :: ekw
    real*8, dimension(:, :), allocatable :: wkw
    real*8, dimension(:, :), allocatable :: gapjack, gap, gapav, gaps

    ! AAA  Lines to be added just after all definitions of variables.
    character(100) :: name_tool
    character(20) :: str
    logical :: yes_read, yes_weight

    call get_command_argument(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then
        ! Input the name of the file exactly as it is in /doc
        name_tool = 'readf'
        call help_online(name_tool)
        stop
    end if

    ! read the index of the considered k-point and
    ! put to zero if no index is specified.
    kp_read = -1
    iflag = 1
    read (str, *, err=101, end=101) kp_read
    iflag = 0
101 if (iflag .eq. 1) kp_read = 1

    write (6, *) ' max k corrections, bin length, ibinit,iskip  ?'
    read (5, *) maxk, lbin, ibinitr, iskipr

    yes_weight = .false.
    if (maxk .lt. 0) then
        yes_weight = .true.
        maxk = 0
    end if

    if (iskipr .ge. 0) then
        open (unit=12, file='fort.12', form='unformatted', status='old', err=102)
        iskip = iskipr
    else
        iskip = -iskipr - 1
        open (unit=12, file='fort.12.fn', form='unformatted', status='old', err=103)
    end if
    if (ibinitr .lt. 0) then
        krep = iskip
        ibinit = -ibinitr
    else
        krep = 1
        ibinit = ibinitr
    end if

    if (iflag .eq. 0) open (unit=37, file='kp_info.dat', form='formatted', status='unknown', err=104)

    yes_read = .false.
    if (lbin .lt. 0) then
        lbin = -lbin
        write (6, *) ' ngen max read ? '
        read (5, *) iread
        yes_read = .true.
    end if

    open (unit=20, file='fort.20', form='formatted', status='unknown')
    open (unit=21, file='fort.21', form='formatted', status='unknown')

    rewind (12)

    if (iflag .eq. 1) then
        nkps = 1
    else
        rewind (37)
        read (37, *) nkps
        if (nkps .eq. 1) kp_read = 1
    end if

    kave_read = 1
    if (kp_read .gt. nkps) then
        kave_read = (kp_read + nkps - 1)/nkps
        kp_read = kp_read - (kave_read - 1)*nkps
        write (6, *) ' kp_read kave_read =', kp_read, kave_read
    end if
    if (nkps .gt. 1) then
        write (6, '(2X,A,X,2I5)') ' number of k-points / index read ', nkps, kp_read
        do i = 1, kp_read
            read (37, *)
        end do
        read (37, *) i, kps(:), wkps
        write (6, '(2X,A,X,4F12.5)') 'kpoint/weight chosen: ', kps(:), wkps
    end if
    close (37)

    nbuf = max(maxk, 1) ! number of buffers for k corrections
    nbufp = nbuf + 1
    npm = iskip + 1 ! number of correlation functions to skip for computing the average

    allocate (ebuf(nbuf, krep), wbuf(nbuf), ek(0:nbuf, krep), wk(0:nbuf), &
              wbufm(nbuf), eskip(npm), wtot(nbuf), cost_mes(krep))

    ! Initialization
    wbufm(:) = 1.d0
    wbuf(:) = 1.d0

    ebuf = 0.d0
    eskip = 0.d0
    wtot = 0.d0

    do i = 0, nbuf
        ek(i, :) = 0.d0
        wk(i) = 0.d0
    end do
    if (yes_read) then
        i = iread
    else
        i = 0
        do while (i .ge. 0)
            read (12, end=123, err=123)
            i = i + 1
        end do
    end if
123 if (nkps .le. 1) then
        countmax = (i/lbin)*lbin
    else
        rest = mod(i, nkps)
        if (rest .ne. 0) go to 107
        countmax = i/nkps
        countmax = (countmax/lbin)*lbin
    end if
    rewind (12)

    nbin = countmax/lbin - ibinit + 1
    ntime = maxk + 1

    if (nbin .le. 0 .or. maxk .lt. 0) go to 108

    allocate (ekw(ntime, nbin, krep), wkw(ntime, nbin))
    allocate (eav(ntime, krep), wav(ntime), gapjack(ntime, krep), gap(ntime, krep)&
            &, gaps(ntime, krep), gapav(ntime, krep))

    icount = 0
    ibin = 0
    ! check
    kt = 0
    countread = 1
    do while (kt .ge. 0 .and. icount .lt. countmax)
        kt = kt + 1
        do i = 1, nbuf
            wbufm(i) = wbuf(i)
        end do

        countread = 0
        do j = 1, nbuf
            if (kave_read .eq. 1) then
                do ikp = 1, nkps
                    if (ikp .eq. kp_read) then
                        if (krep .eq. 1) then
                            if (iskip .ne. 0) then
                                read (12, end=109, err=109) wbuf(j), wtot(j), (eskip(k), k=1, iskip), ebuf(j, 1)
                            else
                                read (12, end=109, err=109) wbuf(j), wtot(j), ebuf(j, 1)
                            end if
                            countread = countread + 1
                        else
                            read (12, end=109, err=109) wbuf(j), wtot(j), eskip(1), (ebuf(j, k), k=1, iskip)
                            countread = countread + 1
                        end if
                        if (yes_weight) then
                            wbuf(j) = 1.d0
                            wtot(j) = 1.d0
                        end if
                    else
                        read (12, end=109, err=109)
                    end if
                end do
            else
                count_mes = 0.d0
                ebuf(j, :) = 0.d0
                do ikp = 1, nkps
                    if (ikp .ge. kp_read .and. ikp .lt. kp_read + kave_read) then
                        if (krep .eq. 1) then
                            if (iskip .ne. 0) then
                                read (12, end=109, err=109) wbuf(j), wtot(j), (eskip(k), k=1, iskip), cost_mes(1)
                            else
                                read (12, end=109, err=109) wbuf(j), wtot(j), cost_mes(1)
                            end if
                        else
                            read (12, end=109, err=109) wbuf(j), wtot(j), eskip(1), (cost_mes(k), k=1, krep)
                        end if

                        if (yes_weight) then
                            wbuf(j) = 1.d0
                            wtot(j) = 1.d0
                        end if

                        ebuf(j, :) = ebuf(j, :) + cost_mes(:)
                        count_mes = count_mes + 1.d0
                        if (ikp .eq. kp_read) countread = countread + 1
                    else
                        read (12, end=109, err=109)
                    end if
                end do
                ebuf(j, :) = ebuf(j, :)/count_mes
            end if
        end do
109     continue

        !      write(6,*) 'countread ',countread,kt,nbuf

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
                    ek(k, :) = ek(k, :) + ebuf(j, :)*wsk
                    wk(k) = wk(k) + wsk
                end do

                !                  write(6,*) ' countread wsk ',icount,wsk/wtot(j),wtot(j)

                if (mod(icount, lbin) .eq. 0) then
                    ibin = ibin + 1
                    if (ibin .ge. ibinit) then
                        do k = 1, krep
                            ekw(1:ntime, ibin - ibinit + 1, k) = ek(0:maxk, k)
                        end do
                        wkw(1:ntime, ibin - ibinit + 1) = wk(0:maxk)
                        if (wk(maxk) .ne. 0.d0) then
#ifdef __KCOMP
                            write (21, '(32767e20.12)') (ek(maxk, k)/wk(maxk), k=1, krep), wk(maxk), dble(ibin)
#else
                            write (21, '(1000000e20.12)') (ek(maxk, k)/wk(maxk), k=1, krep), wk(maxk), dble(ibin)
#endif
                        else
                            write (21, *) 0., wk(maxk), ibin
                        end if
                    end if
                    do k = 0, maxk
                        ek(k, :) = 0.d0
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
                    ek(k, :) = ek(k, :) + ebuf(j, :)*wsk
                    wk(k) = wk(k) + wsk
                    ! end do k
                end do

                !        write(6,*) ' countread wsk ',icount,wsk/wtot(j),wtot(j),ebuf(j)

                if (mod(icount, lbin) .eq. 0) then
                    ibin = ibin + 1
                    if (ibin .ge. ibinit) then
                        do k = 1, krep
                            ekw(1:ntime, ibin - ibinit + 1, k) = ek(0:maxk, k)
                        end do
                        wkw(1:ntime, ibin - ibinit + 1) = wk(0:maxk)
                        if (wk(maxk) .ne. 0.d0) then
#ifdef __KCOMP
                            write (21, '(32767e20.12)') (ek(maxk, k)/wk(maxk), k=1, krep), wk(maxk), dble(ibin)

#else
                            write (21, '(1000000e20.12)') (ek(maxk, k)/wk(maxk), k=1, krep), wk(maxk), dble(ibin)
#endif
                        else
                            write (21, *) 0., wk(maxk), ibin
                        end if
                    end if
                    do k = 0, maxk
                        ek(k, :) = 0.d0
                        wk(k) = 0.d0
                    end do
                end if
            end do
        end if
    end do
    ! Calculation error bars
    nmis = ibin - ibinit + 1
    call jack
    write (6, *) ' number of measures done =', nmis
    write (20, *) ' Independent bins ', nmis, 'of length ', lbin
    write (20, *)
    write (20, *) ' Energy , error, # of bias correcting factor '
    do i = 0, maxk
        write (20, *) (gapjack(i + 1, k), gaps(i + 1, k), k=1, krep), i
    end do
    deallocate (ebuf, wbuf, ek, wk, wbufm, eskip, wtot)
    deallocate (ekw, wkw)
    deallocate (eav, wav, gapjack, gap, gaps, gapav)
    stop

    !!!!!!!!!!! ERRORS !!!!!!!!!!!

102 write (6, *) 'ERROR: file fort.12 not found or empty!'
    stop
103 write (6, *) 'ERROR: file fort.12.fn not found or empty!'
    stop
104 write (6, *) 'ERROR: file kp_info.dat not found or empty!'
    stop
106 write (6, *) 'ERROR: k-point index out of bounds!'
    stop
107 write (6, *) 'ERROR: number of measures not multiple of number of k-points!'
    stop
108 write (6, *) 'ERROR: wrong number of bins or k correcting factors!'
    stop

contains

    subroutine jack
        implicit none
        integer i, k
        eav = 0.d0
        wav = 0.d0
        do i = 1, nbin
            do k = 1, krep
                eav(:, k) = eav(:, k) + ekw(:, i, k)
            end do
            wav(:) = wav(:) + wkw(:, i)
        end do
        do k = 1, krep
            eav(:, k) = eav(:, k)/wav(:)
        end do
        gapjack(:, :) = eav(:, :)
        gaps = 0.d0
        gapav = 0.d0
        do i = 1, nbin
            do j = 1, ntime
                do k = 1, krep
                    gap(j, k) = (eav(j, k)*wav(j) - ekw(j, i, k))/(wav(j) - wkw(j, i))
                end do
            end do
            gapav(:, :) = gapav(:, :) + gap(:, :)
        end do
        gapav = gapav/nbin
        do i = 1, nbin
            do j = 1, ntime
                do k = 1, krep
                    gap(j, k) = (eav(j, k)*wav(j) - ekw(j, i, k))/(wav(j) - wkw(j, i))
                end do
            end do
            gaps(:, :) = gaps(:, :) + (gap(:, :) - gapav(:, :))**2
        end do
        gaps = gaps/nbin
        gaps(:, :) = dsqrt(gaps(:, :))*dsqrt(dble(nbin - 1))
    end subroutine jack

end program erread
