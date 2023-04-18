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

program bootback
    implicit none
    integer nh, nbinm, nm, nmis, n, nbin, iseed, i, j, jj, k, kmain, nel
    parameter(nh=500, nbinm=20000, nm=1)
    real*8 err, eta, esav, eav, et(4), wt(4), par(10), alatb, alat, eps&
            &, givenmeas, drand1, w(nbinm, 4), e(nbinm, 4), deltab, den, num&
            &, val(4), errval(4), eta0(5), errav, errerr, errv, scalepulay
    logical ifneg

    iseed = 12375731
    call rand_init(iseed)

    !        y =  a*x*2+ b*x + c
    !        legge z,y da for011.dat
    eps = 1d-6

    open (unit=11, form='formatted', status='old', file='fort.21.1')
    open (unit=12, form='formatted', status='old', file='fort.21.2')
    open (unit=13, form='formatted', status='old', file='fort.21.3')
    open (unit=14, form='formatted', status='old', file='fort.21.4')

    read (5, *) scalepulay
    if (scalepulay .lt. 0.d0) then
        ifneg = .true.
        if (scalepulay .ge. -100.d0) then
            scalepulay = -scalepulay
        else
            !        really negative value
            scalepulay = scalepulay + 100.d0
        end if
    else
        ifneg = .false.
    end if

    rewind (20)
    n = 1
    !        lettura dati

    nbin = 0
    do while (nbin .lt. nbinm)
        nbin = nbin + 1
        read (11, *, end=500) e(nbin, 1), w(nbin, 1) !  energy
        read (12, *, end=500) e(nbin, 2), w(nbin, 2) !  derrivative of the ener
        read (13, *, end=500) e(nbin, 3), w(nbin, 3) ! log derivative wf
        read (14, *, end=500) e(nbin, 4), w(nbin, 4) ! energy x log der. wf.
    end do
500 continue

    if (nbin .ne. nbinm) then
        nbin = nbin - 1
    else
        write (6, *) ' Warning maximum number of bins exceeded !!!'
    end if

    !      write(6,*) ' number of bins read =',nbin

    nmis = 200

    !       initialization opara e opars
    err = 0.d0
    eta = 0.d0
    eav = 0.d0
    esav = 0.d0
    errav = 0.d0
    errerr = 0.d0

    val(1) = 0.d0
    val(2) = 0.d0
    val(3) = 0.d0
    val(4) = 0.d0

    errval(1) = 0.d0
    errval(2) = 0.d0
    errval(3) = 0.d0
    errval(4) = 0.d0

    wt = 0.d0
    et = 0.d0
    do j = 1, nbin
        do jj = 1, 4
            wt(jj) = wt(jj) + w(j, jj)
            et(jj) = et(jj) + w(j, jj)*e(j, jj)
        end do
    end do

    do jj = 1, 4
        et(jj) = et(jj)/wt(jj)
    end do

    if (ifneg) then
        eta0(1) = -et(2)*scalepulay + 2.d0*(-et(4) + et(1)*et(3))
    else
        eta0(1) = -et(2) + 2.d0*scalepulay*(-et(4) + et(1)*et(3))
    end if

    eta0(2) = -et(2)
    eta0(3) = -et(4)
    eta0(4) = et(1)*et(3)
    eta0(5) = 2.d0*(-et(4) + et(1)*et(3))

    do kmain = 1, nmis

        !
        !       call dscal(4,0.d0,wt,1)
        !       call dscal(4,0.d0,et,1)
        wt = 0.d0
        et = 0.d0

        !       bootstrap estimate sample

        eav = 0.d0
        errv = 0.d0
        do k = 1, nbin
            j = drand1()*nbin + 1
            !       j=k
            do jj = 1, 4
                wt(jj) = wt(jj) + w(j, jj)
                et(jj) = et(jj) + w(j, jj)*e(j, jj)
            end do
            if (ifneg) then
                givenmeas = -e(j, 2)*scalepulay + 2.d0*(-e(j, 4) + e(j, 1)*e(j, 3))
            else
                givenmeas = -e(j, 2) + 2.d0*scalepulay*(-e(j, 4) + e(j, 1)*e(j, 3))
            end if
            eav = eav + givenmeas
            errv = errv + givenmeas**2
        end do
        errv = errv/nbin
        eav = eav/nbin
        errv = dsqrt((errv - eav**2)/nbin)

        errav = errav + errv
        errerr = errerr + errv**2

        do jj = 1, 4
            et(jj) = et(jj)/wt(jj)
        end do

        !       now average correlation function on the given bootstrap

        if (ifneg) then
            givenmeas = -et(2)*scalepulay + 2.d0*(-et(4) + et(1)*et(3))
        else
            givenmeas = -et(2) + 2.d0*scalepulay*(-et(4) + et(1)*et(3))
        end if
        val(1) = val(1) - et(2)
        val(2) = val(2) - et(4)
        val(3) = val(3) + et(1)*et(3)
        val(4) = val(4) + 2.d0*(-et(4) + et(1)*et(3))

        eta = eta + givenmeas
        err = err + givenmeas**2

        errval(1) = errval(1) + et(2)**2
        errval(2) = errval(2) + et(4)**2
        errval(3) = errval(3) + (et(1)*et(3))**2
        errval(4) = errval(4) + (2.d0*(-et(4) + et(1)*et(3)))**2

    end do

    eta = eta/nmis
    err = err/nmis

    val(1) = val(1)/nmis
    errval(1) = errval(1)/nmis
    errval(1) = dsqrt(errval(1) - val(1)**2)

    val(2) = val(2)/nmis
    errval(2) = errval(2)/nmis
    errval(2) = dsqrt(errval(2) - val(2)**2)

    val(3) = val(3)/nmis
    errval(3) = errval(3)/nmis
    errval(3) = dsqrt(errval(3) - val(3)**2)

    val(4) = val(4)/nmis
    errval(4) = errval(4)/nmis
    errval(4) = dsqrt(errval(4) - val(4)**2)

    err = dsqrt(err - eta**2)

    errav = errav/nmis
    errerr = errerr/nmis
    errerr = dsqrt(errerr - errav**2)

    write (6, *) '  Force   =', eta0(1), err, errerr
    write (6, *) '  Der Eloc =', eta0(2), errval(1)
    write (6, *) '  <OH> =', eta0(3), errval(2)
    write (6, *) '  <O><H> =', eta0(4), errval(3)
    write (6, *) '2*(<OH> - <O><H>) =', eta0(5), errval(4)
    !       if(abs(eta0(5)/errval(4)).gt.2) write(6,*) ' Warning Pulay large !!! '

    close (11)
    close (12)
    close (13)
    close (14)

    stop
end

