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
    parameter(nh=500, nbinm=100000, nm=1)
    real*8 err, eta, esav, eav, et(3), wt(3), alatb, alat, eps&
            &, givenmeas, drand1, w(nbinm, 3), e(nbinm, 3), deltab, den, num, eta0, errav, errerr&
            &, errv
    iseed = 12375731
    call rand_init(iseed)

    !        y =  a*x*2+ b*x + c
    !        legge z,y da for011.dat
    eps = 1d-6

    open (unit=11, form='formatted', status='old', file='fort.21.1')
    open (unit=12, form='formatted', status='old', file='fort.21.2')
    open (unit=13, form='formatted', status='old', file='fort.21.3')

    rewind (20)
    n = 1
    !        lettura dati

    nbin = 0
    do while (nbin .lt. nbinm)
        nbin = nbin + 1
        read (11, *, end=500) e(nbin, 1), w(nbin, 1) !  energy
        read (12, *, end=500) e(nbin, 2), w(nbin, 2) ! log  derivative wf
        read (13, *, end=500) e(nbin, 3), w(nbin, 3) ! log derivative wf x local energy
    end do
500 continue

    if (nbin .ne. nbinm) then
        nbin = nbin - 1
    else
        write (6, *) ' Warning maximum number of bins exceeded !!!'
    end if

    write (6, *) ' number of bins read =', nbin

    nmis = 200

    !       initialization opara e opars
    err = 0.d0
    eta = 0.d0
    eav = 0.d0
    esav = 0.d0

    wt = 0.d0
    et = 0.d0
    do j = 1, nbin
        do jj = 1, 3
            wt(jj) = wt(jj) + w(j, jj)
            et(jj) = et(jj) + w(j, jj)*e(j, jj)
        end do
    end do

    do jj = 1, 3
        et(jj) = et(jj)/wt(jj)
    end do
    eta0 = 2.d0*(-et(3) + et(1)*et(2))

    errav = 0.d0
    errerr = 0.d0

    do kmain = 1, nmis
        !
        !       call dscal(3,0.d0,wt,1)
        !       call dscal(3,0.d0,et,1)
        wt = 0.d0
        et = 0.d0

        !       bootstrap estimate sample
        eav = 0.d0
        errv = 0.d0

        do k = 1, nbin
            j = drand1()*nbin + 1
            do jj = 1, 3
                wt(jj) = wt(jj) + w(j, jj)
                et(jj) = et(jj) + w(j, jj)*e(j, jj)
            end do
            givenmeas = 2.d0*(-e(j, 3) + e(j, 1)*e(j, 2))
            eav = eav + givenmeas
            errv = errv + givenmeas**2
        end do

        errv = errv/nbin
        eav = eav/nbin
        errv = dsqrt((errv - eav**2)/nbin)

        errav = errav + errv
        errerr = errerr + errv**2

        do jj = 1, 3
            et(jj) = et(jj)/wt(jj)
        end do

        !       now average correlation function on the given bootstrap

        givenmeas = 2.d0*(-et(3) + et(1)*et(2))

        eta = eta + givenmeas
        err = err + givenmeas**2

    end do

    eta = eta/nmis
    err = err/nmis
    errav = errav/nmis
    errerr = errerr/nmis

    errerr = dsqrt(errerr - errav**2)

    err = dsqrt(err - eta**2)

    if (abs(eta0) .gt. 3*err) then
        write (6, *) '  Force   =', eta0, err, errerr&
                &, ' Warning large fluctuation=', abs(eta/err)
    else
        write (6, *) '  Force   =', eta0, err, errerr
    end if

    close (11)
    close (12)
    close (13)
    close (14)

    stop
end

