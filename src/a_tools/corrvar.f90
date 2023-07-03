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

! Given in input the block averages of energy and variance, this program
! computes the global average and standard deviation of these quantities
! using bootstrap resampling technique.

program bootback

    implicit none

    integer :: iseed, nbin, n, nbinm, k, j, i, nmis, kmain, nwalk, bin_length
    real(8), dimension(:), allocatable :: w, e, es
    real(8) :: wpip, drand1
    real(8) :: eps, wt, et, ets, etsv, variance
    real(8) :: errn, givenmeas, givenmeas2
    real(8) :: eta, eav, esav, eav0, err
    real(8) :: eta_, eav_, esav_, err_
    logical file_exists
    iseed = 12375731
    call rand_init(iseed)

    ! open files with block-average of the energy fort.21
    ! In case of k-points calculations it already contains
    ! the block-average of the variance.
    open (unit=11, form='formatted', status='old', file='fort.21', err=101)
    open (unit=12, form='formatted', status='old', file='fort.22', err=102)
    open (unit=13, form='formatted', status='old', file='parminimized.d', err=100)

    !       inquire(file='kp_info.dat',exist=file_exists)
    !       if(file_exists) then
    !       open(unit=14,form='formatted',status='old',file='kp_info.dat')
    !       read(14,*) nk
    !       else
    !       nk=1
    !       endif

    read (13, *)
    read (13, *) nwalk
    close (13)

    ! nwalk=nwalk/nk

    read (5, *) bin_length

    ! inquire if k-points are present

    n = 1
    eps = 1d-6

    ! reading datas
    nbin = 0
    do while (nbin .ge. 0)
        read (11, *, err=101, end=500)
        read (12, *, err=102, end=500)
        nbin = nbin + 1
    end do
500 continue
    rewind (11)
    rewind (12)
    allocate (w(nbin), e(nbin), es(nbin))
    do i = 1, nbin
        read (11, *, err=101, end=500) e(i), w(i)
        read (12, *, err=102, end=500) es(i), wpip
        if (wpip .ne. w(i)) then
            write (6, *) ' warning the measures are not correlated !!!! '
        end if
    end do
    close (11)
    close (12)

    !  write(6,*) ' Measures read '
    !  do i=1,nbin
    !  write(6,*) i,e(i),es(i),sqrt(es(i)-e(i)**2)
    !  enddo

    eav0 = sum(e(:)*w(:))/sum(w(:))
    write (6, *) ' number of bins read =', nbin

    nmis = 1000
    err = 0.d0
    eta = 0.d0
    eav = 0.d0
    esav = 0.d0
    err_ = 0.d0
    eta_ = 0.d0
    eav_ = 0.d0
    esav_ = 0.d0

    do kmain = 1, nmis
        !
        wt = 0.d0
        et = 0.d0
        ets = 0.d0
        etsv = 0.d0
        !       bootstrap estimate sample
        do k = 1, nbin
            j = drand1()*nbin + 1
            wt = wt + w(j)
            et = et + w(j)*e(j)
            ets = ets + w(j)*es(j) ! correlated estimate of <E**2>
            etsv = etsv + w(j)*e(j)**2 ! uncorrelated estimate of <E**2>
        end do

        !       now average correlation functions on the given bootstrap
        et = et/wt
        ets = ets/wt
        etsv = etsv/wt
        errn = dsqrt(abs(etsv - et**2)/nbin)
        givenmeas = ets - et**2
        eta = eta + givenmeas
        err = err + givenmeas**2
        eav = eav + et
        esav = esav + et**2
        variance = ets - et**2
        givenmeas2 = errn
        eta_ = eta_ + givenmeas2
        err_ = err_ + givenmeas2**2
        givenmeas2 = nbin*bin_length*errn**2/variance
        eav_ = eav_ + givenmeas2
        esav_ = esav_ + givenmeas2**2

    end do

    !      normalize measures
    eta = eta/nmis
    err = err/nmis
    esav = esav/nmis
    eav = eav/nmis
    err = sqrt(err - eta**2)
    esav = sqrt(esav - eav**2)

    eta_ = eta_/nmis
    err_ = err_/nmis
    esav_ = esav_/nmis
    eav_ = eav_/nmis
    err_ = dsqrt(err_ - eta_**2)
    esav_ = dsqrt(esav_ - eav_**2)

    write (6, *) ' Energy =', dble(eav0), dble(esav)
    write (6, *) ' Variance square =', dble(eta), dble(err)
    write (6, *) ' Est. energy error bar =', eta_, err_
    write (6, *) ' Est. corr. time  =', eav_*nwalk, esav_*nwalk

    deallocate (e, es, w)
    stop

    !!!!!!!!!!! ERRORS !!!!!!!!!!!
100 write (6, *) 'ERROR: file parminimized.d not found or wrong!'
    stop
101 write (6, *) 'ERROR: file fort.21 not found or wrong!'
    stop
102 write (6, *) 'ERROR: file fort.22 not found or wrong!'
    stop
103 write (6, *) 'ERROR: file kp_info.dat empty or wrong!'
    stop

end program bootback

