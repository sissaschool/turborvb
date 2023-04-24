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
    integer nh, npar, nm, nmis, n, nbin, iseed, i, j, jj, kk, k, kmain, nel&
            &, lwork, info
    parameter(nh=500, nm=1)
    real*8 err, esav, eav, par(10), et(4), wt(4), alatb, alat, eps&
            &, givenmeas, drand1, deltab, den, num, cost, scalef, maxsn, scalepulay&
            &, val(4), errval(4), eta0(5), errav, errerr, errv, weight, weightall
    integer ir(7)
    real*8, dimension(:, :), allocatable :: cov, fk, w, e, eall
    real*8, dimension(:, :, :), allocatable :: ef, wf
    real*8, dimension(:), allocatable :: fkav, psip, eig, eta

    iseed = 12375731
    call rand_init(iseed)

    !        y =  a*x*2+ b*x + c
    !        legge z,y da for011.dat
    eps = 1d-6

    open (unit=11, form='formatted', status='old', file='fort.21.1')
    open (unit=12, form='formatted', status='old', file='fort.21.2')
    open (unit=13, form='formatted', status='old', file='fort.21.3')
    open (unit=14, form='formatted', status='old', file='fort.21.4')
    open (unit=15, form='formatted', status='old', file='parminimized.d')

    read (15, *) ir(1:7)

    npar = ir(7)
    lwork = 3*npar

    read (5, *) scalepulay

    rewind (20)
    n = 1
    !        lettura dati

    nbin = 0
    do while (nbin .ge. 0)
        read (11, *, end=500) et(1:2) !  energy
        nbin = nbin + 1
        !       read(12,*,end=500) e(nbin,2),w(nbin,2)  ! derivative of the ener
        !       read(13,*,end=500) e(nbin,3),w(nbin,3)  ! log derivative wf
        !       read(14,*,end=500) e(nbin,4),w(nbin,4)  ! energy x log der. wf.
    end do
500 continue

    allocate (fk(nbin, npar), ef(nbin, npar, 3), wf(nbin, npar, 3))
    allocate (e(nbin, 1), w(nbin, 1))
    allocate (fkav(npar), eta(npar), eall(4, npar), cov(npar, npar))
    allocate (eig(npar), psip(lwork))
    rewind (11)
    do i = 1, nbin
        read (11, *) e(i, 1), w(i, 1)
        read (12, *) (ef(i, k, 1), wf(i, k, 1), ir(1), k=1, npar)
        read (13, *) (ef(i, k, 2), wf(i, k, 2), ir(1), k=1, npar)
        read (14, *) (ef(i, k, 3), wf(i, k, 3), ir(1), k=1, npar)

    end do

    write (6, *) ' number of bins read =', nbin

    nmis = nbin
    if (nbin .lt. 10*npar) then
        write (6, *) ' Warning too small number of bins !!!  '
    end if

    !       initialization opara e opars
    err = 0.d0
    eav = 0.d0
    esav = 0.d0

    val = 0.d0
    errval = 0.d0
    wt = 0.d0
    et = 0.d0
    do j = 1, nbin
        wt(1) = wt(1) + w(j, 1)
        et(1) = et(1) + e(j, 1)*w(j, 1)
    end do
    weightall = wt(1)
    et(1) = et(1)/wt(1)

    do kk = 1, npar
        wt(2:4) = 0.d0
        et(2:4) = 0.d0
        do j = 1, nbin
            do jj = 1, 3
                wt(jj + 1) = wt(jj + 1) + wf(j, kk, jj)
                et(jj + 1) = et(jj + 1) + wf(j, kk, jj)*ef(j, kk, jj)
            end do
        end do

        do jj = 2, 4
            et(jj) = et(jj)/wt(jj)
        end do
        eall(1:4, kk) = et(1:4)

        eta0(1) = -et(2) + 2.d0*scalepulay*(-et(4) + et(1)*et(3))

        eta0(2) = -et(2)
        eta0(3) = -et(4)
        eta0(4) = et(1)*et(3)
        eta0(5) = 2.d0*(-et(4) + et(1)*et(3))
        fkav(kk) = eta0(1)
    end do

    eta = 0.d0 ! Jackknife average.

    do kk = 1, npar

        do kmain = 1, nmis

            weight = w(kmain, 1)
            !       jacknife  estimate sample

            et(1) = (eall(1, kk)*weightall - weight*e(kmain, 1))/(weightall - weight)
            do k = 1, 3
                et(k + 1) = (eall(k + 1, kk)*weightall - weight*ef(kmain, kk, k))/(weightall - weight)
            end do

            givenmeas = -et(2) + 2.d0*scalepulay*(-et(4) + et(1)*et(3))
            fk(kmain, kk) = givenmeas
            eta(kk) = eta(kk) + givenmeas

        end do
    end do

    cost = dsqrt(dble(max(nbin - 1, 1))/dble(max(nbin, 1)))

    eta = eta/nmis

    do j = 1, npar
        do i = 1, nbin
            fk(i, j) = cost*(fk(i, j) - eta(j))
        end do
    end do
    write (6, *) ' Force before covariance '
    do i = 1, npar, 3
        write (6, *) (i - 1)/3 + 1, (fkav(i + k), k=0, 2)
    end do

    call dgemm('T', 'N', npar, npar, nbin, 1.d0, fk, nbin, fk, nbin, 0.d0, cov, npar)

    write (6, *) ' Covariance matrix '
    do i = 1, npar
        do j = i, npar
            write (6, *) i, j, cov(i, j)
        end do
    end do

    call dsyev('V', 'L', npar, cov, npar, eig, psip, lwork, info)
    if (info .ne. 0) write (6, *) ' ERROR in diagonalization '
    write (6, *) ' Eigenvalues covariance matrix ', eig(1)/eig(npar)
    do i = 1, npar
        write (6, *) i, eig(i)
    end do
    call dgemv('T', npar, npar, 1.d0, cov, npar, fkav, 1, 0.d0, psip, 1)
    maxsn = 0.d0
    do i = 1, npar
        if (i .le. 3) then
            psip(i) = 0.d0
        else
            psip(i) = psip(i)/eig(i)
            maxsn = maxsn + psip(i)**2*eig(i)
        end if
    end do
    if (maxsn .gt. 0) maxsn = dsqrt(maxsn)
    call dgemv('N', npar, npar, 1.d0, cov, npar, psip, 1, 0.d0, fkav, 1)

    write (6, *) ' Direction maximum signal/noise ratio, value = ', maxsn
    scalef = 0.d0
    do i = 1, npar
        scalef = scalef + eig(i)
    end do
    scalef = scalef/npar

    write (6, *) ' force after covariance '
    do i = 1, npar, 3
        write (6, *) (i - 1)/3 + 1, (fkav(i + k)*scalef, k=0, 2)
    end do

    close (11)
    close (12)
    close (13)
    close (14)
    close (15)

    stop
end

