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

program fitxyz
    implicit none
    real(8), dimension(:), allocatable :: x, erx, y, z, er, ys, xs, xt, yn&
            &, ern, yns, alat
    real(8) eps, pi, am, cm, delm, sam, scm, sdelm, a, c, an, cn, stot, as, cs&
            &, r1, r2, zeta(2), del(2), xmax, drand1, era, erc, yp, cost, cmat(2, 2)&
            &, eig(2), psip(20), anm, cnm, sanm, scnm, stotn, ans, cns, eran, ercn
    integer i, n, nh, nmis, imis, nmi, iopt, iter, imax, iseed, lwork, info
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'evsa'
        !           call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added

    open (unit=11, file='evsa.d', form='formatted', status='old')
    open (unit=20, file='evsa.fitl', form='formatted', status='unknown')
    !        y = a*z +  c
    !        legge z,y da for011.dat
    eps = 1e-12
    lwork = 20
    !        lettura dati
    rewind (11)
    n = 0
    do while (n .ge. 0)
        read (11, *, end=500) psip(1:9)
        n = n + 1
    end do
500 continue
    allocate (x(n), erx(n), y(n), z(n), er(n), ys(n), xs(n), xt(n), yn(n)&
            &, ern(n), yns(n), alat(n))

    nh = n

    rewind (11)
    do n = 1, nh
        read (11, *) alat(n), y(n), er(n), yn(n), ern(n), a, c, xt(n), erx(n)
        write (6, *) xt(n), erx(n), y(n), er(n), yn(n), ern(n)
    end do
    n = nh

    do i = 1, n
        if (er(i) .ne. 0) then
            er(i) = 1.0/er(i)**2
        else
            er(i) = 1.d0/eps
        end if
        if (ern(i) .ne. 0) then
            ern(i) = 1.0/ern(i)**2
        else
            ern(i) = 1.d0/eps
        end if
        if (erx(i) .ne. 0) then
            erx(i) = 1.0/erx(i)**2
        else
            erx(i) = 1.d0/eps
        end if
    end do

    nmis = 10000

    iseed = 3523
    call rand_init(iseed) ! ibm
    pi = 2.d0*dasin(1.d0)
    am = 0.
    cm = 0.
    anm = 0.
    cnm = 0.
    delm = 0.

    sam = 0.
    scm = 0.
    sanm = 0.
    scnm = 0.
    sdelm = 0.

    nmi = 0
    iopt = 0
    do imis = 1, nmis
        iter = 0
        call fun(xt, y, er, n, stot, a, c)
        call fun(xt, yn, ern, n, stotn, an, cn)

        if (imis .eq. 1) then
            !        update ys
            as = a
            cs = c
            ans = an
            cns = cn
            do i = 1, n
                ys(i) = y(i)
                yns(i) = yn(i)
                xs(i) = xt(i)
            end do

            write (20, *)
            write (20, *) 'Linear fit'
            write (20, *) ' sum weighted squares = ', stot, stotn

        else
            !        generate new set of y random distributed
            do i = 1, n

                !        calculation 2x2 covariance matrix
                cmat(1, 1) = 1.d0/er(i)
                cmat(2, 2) = 1.d0/ern(i)
                cmat(1, 2) = (cmat(1, 1) + cmat(2, 2) - 1.d0/erx(i))/2.d0
                cmat(2, 1) = cmat(1, 2)
                !        diagonalize the matrix
                call dsyev('V', 'L', 2, cmat, 2, eig, psip, lwork, info)

                r1 = drand1()
                r2 = drand1()
                zeta(1) = dsqrt(-2.*dlog(r1)*eig(1))*dcos(2.*pi*r2)
                zeta(2) = dsqrt(-2.*dlog(r1)*eig(2))*dsin(2.*pi*r2)

                call dgemv('N', 2, 2, 1.d0, cmat, 2, zeta, 1, 0.d0, del, 1)

                y(i) = ys(i) + del(1)
                yn(i) = yns(i) + del(2)
                xt(i) = xs(i) + del(2) - del(1)

                !         write(6,*) ' y new =',y(i),' del =',del
            end do
        end if
        !        WRITE(6,*) A,B,C
        am = am + a
        cm = cm + c
        sam = sam + a**2
        scm = scm + c**2
        anm = anm + an
        cnm = cnm + cn
        sanm = sanm + an**2
        scnm = scnm + cn**2
        nmi = nmi + 1

    end do ! enddo main if

    cost = 1.d0/dble(nmi)
    am = am*cost
    cm = cm*cost
    sam = sam*cost
    scm = scm*cost
    era = sqrt(abs(am**2 - sam))
    erc = sqrt(abs(cm**2 - scm))

    anm = anm*cost
    cnm = cnm*cost
    sanm = sanm*cost
    scnm = scnm*cost
    eran = sqrt(abs(anm**2 - sanm))
    ercn = sqrt(abs(cnm**2 - scnm))

    !
    write (20, *) 'EXTRAPOLATION a-->0 '
    write (20, *) 'Slope (MA)  = ', as, ' +/- ', era
    write (20, *) 'Energy (MA) = ', cs, ' +/- ', erc
    write (20, *) 'Slope (H^a)  = ', ans, ' +/- ', eran
    write (20, *) 'Energy (H^a)  = ', cns, ' +/- ', ercn

    !       imax=1
    !       xmax=0.
    !       do i=1,n
    !       if(xs(i).gt.xmax) then
    !       imax=i
    !       xmax=xs(i)
    !       endif
    !       yp=as*xs(i) +cs
    !       write(20,*) ' predicted --> ',yp, ' measured --> ',ys(i)
    !       enddo

    !       yp=as*xs(imax)  +cs

    !       write(21,*) 0.,cs,erc
    !       write(21,*) xs(imax),yp,0.

    stop
end

subroutine fun(z, y, er, n, stot, a, c)
    implicit none
    real(8) y(*), z(*), er(*), somz, somy, somzy, somzq, wtot, stot, yp, a, c&
            &, Szy, Szz
    integer i, n
    somz = 0.0
    somy = 0.0
    somzy = 0.0
    somzq = 0.0

    wtot = 0.d0

    do i = 1, n
        somy = somy + y(i)*er(i)
        somz = somz + z(i)*er(i)
        somzy = somzy + z(i)*y(i)*er(i)
        somzq = somzq + z(i)**2*er(i)
        wtot = wtot + er(i)
    end do
    somy = somy/wtot
    somz = somz/wtot
    somzy = somzy/wtot
    somzq = somzq/wtot

    Szy = somzy - somz*somy
    Szz = somzq - somz**2

    !  ONLY TWO PARAMETER FIT    ....
    a = Szy/Szz
    c = somy - a*somz
    stot = 0.0
    do i = 1, n
        yp = a*z(i) + c
        stot = stot + (yp - y(i))**2*er(i)
    end do

    return
end

