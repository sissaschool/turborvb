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
    real*8 err(2), eta(2), esav(2), eav(2), et(6), wt(6)&
            &, givenmeas(2), drand1, eta0(2), errav(2), errerr(2), errv(2)
    real*8, dimension(:, :), allocatable :: w, e
    iseed = 12375731
    call rand_init(iseed)

    !        y =  a*x*2+ b*x + c
    !        legge z,y da for011.dat

    open (unit=11, form='formatted', status='old', file='fort.21.1')
    open (unit=12, form='formatted', status='old', file='fort.21.2')
    open (unit=13, form='formatted', status='old', file='fort.21.3')
    open (unit=14, form='formatted', status='old', file='fort.21.4')
    open (unit=15, form='formatted', status='old', file='fort.21.5')
    open (unit=16, form='formatted', status='old', file='fort.21.1i')

    n = 1
    !        lettura dati

    nbin = 0
    do while (nbin .ge. 0)
        read (11, *, end=500) !  energy
        read (16, *, end=500) ! imaginary part local energy
        read (12, *, end=500) !  derivative of the ener
        read (13, *, end=500) ! log derivative wf x local energy
        read (14, *, end=500) !  derivative of the ener  Im
        read (15, *, end=500) ! log derivative wf x local energy Im
        nbin = nbin + 1
    end do
500 continue
    allocate (e(nbin, 6), w(nbin, 6))
    rewind (11)
    rewind (12)
    rewind (13)
    rewind (14)
    rewind (15)
    rewind (16)

    do i = 1, nbin
        read (11, *) e(i, 1), w(i, 1) !  energy
        read (16, *) e(i, 6), w(i, 6) ! imaginary part local energy
        read (12, *) e(i, 2), w(i, 2) !  derivative of the ener
        read (13, *) e(i, 3), w(i, 3) ! log derivative wf x local energy
        read (14, *) e(i, 4), w(i, 4) !  derivative of the ener  Im
        read (15, *) e(i, 5), w(i, 5) ! log derivative wf x local energy Im
    end do

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
        do jj = 1, 6
            wt(jj) = wt(jj) + w(j, jj)
            et(jj) = et(jj) + w(j, jj)*e(j, jj)
        end do
    end do

    do jj = 1, 6
        et(jj) = et(jj)/wt(jj)
    end do
    eta0(1) = 2.d0*(-et(3) + et(1)*et(2) - et(4)*et(6))
    eta0(2) = 2.d0*(-et(5) + et(6)*et(2) + et(4)*et(1))
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
            do jj = 1, 6
                wt(jj) = wt(jj) + w(j, jj)
                et(jj) = et(jj) + w(j, jj)*e(j, jj)
            end do
            givenmeas(1) = 2.d0*(-e(j, 3) + e(j, 1)*e(j, 2) - e(j, 4)*e(j, 6))
            givenmeas(2) = 2.d0*(-e(j, 5) + e(j, 1)*e(j, 4) + e(j, 2)*e(j, 6))
            eav = eav + givenmeas
            errv = errv + givenmeas**2
        end do

        errv = errv/nbin
        eav = eav/nbin
        errv = dsqrt((errv - eav**2)/nbin)

        errav = errav + errv
        errerr = errerr + errv**2

        do jj = 1, 6
            et(jj) = et(jj)/wt(jj)
        end do

        !       now average correlation function on the given bootstrap

        givenmeas(1) = 2.d0*(-et(3) + et(1)*et(2) - et(4)*et(6))
        givenmeas(2) = 2.d0*(-et(5) + et(1)*et(4) + et(2)*et(6))

        eta = eta + givenmeas
        err = err + givenmeas**2

    end do

    eta = eta/nmis
    err = err/nmis
    errav = errav/nmis
    errerr = errerr/nmis

    errerr = dsqrt(errerr - errav**2)

    err = dsqrt(err - eta**2)

    write (6, *) '  Real Force   =', eta0(1), err(1), errerr(1)
    if (abs(eta(1)) .gt. 3*err(1)) then
        write (6, *) ' Warning,  this parameter is not at minimum !!! '&
                &, abs(eta(1)/err(1))
    end if
    write (6, *) '  Imag Force   =', eta0(2), err(2), errerr(2)
    if (abs(eta(2)) .gt. 3*err(2)) then
        write (6, *) ' Warning,  this parameter is not at minimum !!! '&
                &, abs(eta(2)/err(2))
    end if

    close (11)
    close (12)
    close (13)
    close (14)
    close (15)
    close (16)
    deallocate (e, w)

    stop
end

