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

program fitvsa
    implicit none
    real*8 res, errscale, r1, r2, zeta, pi, drand1, xk, xmax, power&
            &, ermax, xw, evidencem, evidences
    integer order, orderp, N, iopt, i, j, iseed, itry, nbin, nmis, Ntrue

    real(8), dimension(:, :), allocatable :: psip
    integer, dimension(:), allocatable :: ipsip
    real(8), dimension(:), allocatable :: x, y, wy, ypred, dypred, wx, coeffa&
            &, coefferr, yt, xt, fxt, fxa, fxerr, ypreda, yprederr, dypreda, dyprederr
#ifdef __DDIAG
    real*16, dimension(:), allocatable :: coeff, coefft
#else
    real*8, dimension(:), allocatable :: coeff, coefft
#endif

    logical checkconv, fixbeta
    real*16 xd, xtr, cost, dcost
    double precision, external :: dlamch
    double precision map, beta, gamma, evidence
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'funvsa'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added

    write (6, *) ' Program with precision =', dlamch('E')

    map = 0.d0
    beta = 1.d0
    fixbeta = .false.
    read (5, *) order, N, iopt, power
    if (N .lt. 0) then
        N = -N
        map = 1.d0
        if (order .lt. 0) then
            fixbeta = .true.
            order = -order
        end if
    end if
    orderp = order + 1
    allocate (x(N), y(N), yt(N), wy(N), coeff(orderp), psip(orderp, orderp + 1)&
            &, ipsip(orderp), coefft(orderp), coeffa(orderp), coefferr(orderp)&
            &, ypred(N), ypreda(N), yprederr(N), dypred(N), dypreda(N), dyprederr(N))

    if (iopt .eq. 2) then
        read (5, *) errscale
    end if

    if (iopt .eq. 0) then
        do i = 1, N
            read (5, *) x(i), y(i), wy(i)
        end do
    elseif (iopt .eq. 1) then
        allocate (wx(N), xt(N))
        do i = 1, N
            read (5, *) x(i), wx(i), y(i), wy(i)
        end do
    elseif (iopt .eq. 3) then
        do i = 1, N
            read (5, *) x(i), y(i), wy(i)
            x(i) = dsqrt(x(i))
        end do
    elseif (iopt .eq. 4) then
        do i = 1, N
            read (5, *) x(i), y(i), wy(i)
            x(i) = x(i)**2
        end do
    elseif (iopt .eq. 5) then
        do i = 1, N
            read (5, *) x(i), y(i), wy(i)
            x(i) = 1/dsqrt(x(i))
        end do
    elseif (iopt .eq. 6) then
        do i = 1, N
            read (5, *) x(i), y(i), wy(i)
            x(i) = 1/x(i)
        end do
    elseif (iopt .eq. 7) then
        do i = 1, N
            read (5, *) x(i), y(i), wy(i)
            x(i) = 1/x(i)**2
        end do
    elseif (iopt .eq. 8) then
        do i = 1, N
            read (5, *) x(i), y(i), wy(i)
            x(i) = 1/log(x(i))
        end do
    else
        do i = 1, N
            read (5, *) x(i), y(i)
            !     x(i)=x(i)**power
            wy(i) = 1.d0
        end do
        if (iopt .eq. 2) then
            do i = 1, N
                wy(i) = errscale
            end do
        end if
    end if
    write (6, *) ' Read data OK '
    xmax = x(1)
    do i = 2, N
        if (x(i) .gt. xmax) xmax = x(i)
    end do

    iseed = 3523
    nbin = 200

    call rand_init(iseed) ! ibm
    pi = 2.d0*dasin(1.d0)

    Ntrue = 0
    do i = 1, N
        if (wy(i) .ne. 0.d0) Ntrue = Ntrue + 1
    end do

    call funmult(order, N, Ntrue, x, y, wy, ypred, dypred, coeff, res, psip, ipsip, power, map, beta, gamma, evidence, fixbeta)
    if (map .ne. 0.d0) write (6, *) ' Machine learning alpha,beta =', map, beta
    if (Ntrue .gt. orderp) then
        write (6, *) ' Reduced chi^2  =', res/dble(Ntrue - orderp)
    else
        write (6, *) ' Warning, no degrees of freedom '
    end if
    !     estimation error bars
    coeffa = 0.d0
    coefferr = 0.d0
    ypreda = 0.d0
    yprederr = 0.d0
    dypreda = 0.d0
    dyprederr = 0.d0
    evidencem = 0.d0
    evidences = 0.d0
    nmis = 0
    itry = 0
    do while (nmis .lt. nbin)
        itry = itry + 1
        do i = 1, N
            r1 = drand1()
            r2 = drand1()
            zeta = dsqrt(-2.*dlog(1.d0 - r1))*dcos(2.*pi*r2)
            yt(i) = y(i) + zeta*wy(i)
        end do
        if (iopt .eq. 1) then
            do i = 1, N
                r1 = drand1()
                r2 = drand1()
                zeta = dsqrt(-2.*dlog(1.d0 - r1))*dcos(2.*pi*r2)
                xt(i) = x(i) + zeta*wx(i)
            end do
            call funmult(order, N, Ntrue, xt, yt, wy, ypred, dypred, coefft &
                         , res, psip, ipsip, power, map, beta, gamma, evidence, fixbeta)
        else
            call funmult(order, N, Ntrue, x, yt, wy, ypred, dypred, coefft &
                         , res, psip, ipsip, power, map, beta, gamma, evidence, fixbeta)
        end if
        nmis = nmis + 1
        coeffa = coeffa + coefft
        coefferr = coefferr + coefft**2
        ypreda = ypreda + ypred
        yprederr = yprederr + ypred**2
        dypreda = dypreda + dypred
        dyprederr = dyprederr + dypred**2
        evidencem = evidencem + evidence
        evidences = evidences + evidence**2
    end do

    coefferr = coefferr/nmis
    coeffa = coeffa/nmis
    yprederr = yprederr/nmis
    ypreda = ypreda/nmis
    dyprederr = dyprederr/nmis
    dypreda = dypreda/nmis
    evidencem = evidencem/nmis
    evidences = evidences/nmis

    coefferr = dsqrt(max(coefferr - coeffa**2, 0.d0))
    yprederr = dsqrt(max(yprederr - ypreda**2, 0.d0))
    dyprederr = dsqrt(max(dyprederr - dypreda**2, 0.d0))
    evidences = dsqrt(max(evidences - evidencem**2, 0.d0))
    if (map .ne. 0.d0) write (6, *) ' Machine learning Evidence =', real(evidencem), '+/-', real(evidences)
    write (6, *) ' Coefficient found '
    do i = 1, orderp
        write (6, *) i, coeff(i), coefferr(i)
    end do

    write (6, *) ' Predicted error /measured error '

    ermax = 0.d0
    do i = 1, N
        cost = coeff(1)
        dcost = 0.d0
        xd = x(i)
        if (power .eq. 1.d0) then
            xtr = 1.d0
        else
            xtr = xd**(power - 1.d0)
        end if
        do j = 2, orderp
            dcost = dcost + (j - 2 + power)*(xtr*coeff(j))
            xtr = xtr*xd
            cost = cost + coeff(j)*xtr
        end do
        !      Go back to the original definition of x
        if (iopt .eq. 3) then
            xw = xd**2
        elseif (iopt .eq. 0 .or. iopt .eq. 2) then
            xw = xd
        elseif (iopt .eq. 4) then
            xw = sqrt(xd)
        elseif (iopt .eq. 5) then
            xw = 1.d0/xd**2
        elseif (iopt .eq. 6) then
            xw = 1.d0/xd
        elseif (iopt .eq. 7) then
            xw = 1.d0/sqrt(xd)
        elseif (iopt .eq. 8) then
            xw = exp(1.d0/xd)
        end if

        if (wy(i) .ne. 0.d0) then
            if (abs(cost - y(i)) .gt. ermax) ermax = abs(cost - y(i))
            write (6, 123) i, xw, cost, yprederr(i), dcost, dyprederr(i), y(i), wy(i)
        else
            write (6, 123) i, xw, cost, yprederr(i), dcost, dyprederr(i)
        end if
    end do
    write (6, *) ' Max error in fit =', ermax

123 format(I4, 7f15.7)

    stop
end

subroutine funmult(order, N, Ntrue, x, y, wy, ypred, dypred, coeff, res, psip, ipsip, power&
        &, map, beta, gamma, evidence, fixbeta)
    implicit none
    integer i, j, k, l, order, orderp, orderpp, orderp3, N, ipsip(*), info, iter, maxit, Ntrue
#ifdef __DDIAG
    real*16 coeff(order + 1)
#else
    real*8 coeff(order + 1)
#endif
    real*8 y(N), x(N), wy(N), psip(order + 1, *), res, power&
            &, ypred(N), dypred(N), error, eps
    double precision, dimension(:, :), allocatable :: mat
    double precision, dimension(:), allocatable :: eig, work, vec
    real*16, dimension(:), allocatable :: coeffd, vecd
    real*16, dimension(:, :), allocatable :: psipd
    real*16 xd, xk, fx, dfx, resd
    double precision, parameter :: TWO_PI = 6.28318530717958647692d0
    double precision, parameter :: one = 1.d0
    double precision, parameter :: zero = 0.d0

    logical fixbeta
    double precision map, beta, gamma, evidence, eig_min
    double precision, external :: dlamch

    orderp = order + 1
    orderpp = order + 2
    !  Implementing quadruple precision in the O(N^2) operations (load matrix)
    allocate (coeffd(orderp), vecd(orderp), psipd(orderp, orderp))

    coeffd = 0
    vecd = 0
    psipd = 0
    do i = 1, orderp
        do j = 1, orderp
            psip(i, j) = 0.d0
        end do
        coeff(i) = 0.d0
    end do

    do i = 1, N
        if (wy(i) .ne. 0.d0) then
            xd = x(i)
            if (power .gt. 1) then
                vecd(1) = xd**(power - 1)
                !                psip(1,orderpp)=x(i)**(power-1)
            elseif (abs(power) .eq. 1) then
                vecd(1) = 1
                !                psip(1,orderpp)=1.d0
            else
                vecd(1) = xd**(power + 1)
                !                psip(1,orderpp)=x(i)**(power+1)
            end if
            if (power .gt. 0) then
                do k = 2, orderp
                    !                psip(k,orderpp)=psip(k-1,orderpp)*x(i)
                    vecd(k) = vecd(k - 1)*xd
                end do
            else
                do k = 2, orderp
                    !                psip(k,orderpp)=psip(k-1,orderpp)/x(i)
                    vecd(k) = vecd(k - 1)/xd
                end do
            end if
            !                psip(1,orderpp)=1.d0
            vecd(1) = 1
            do k = 1, orderp
                vecd(k) = vecd(k)/wy(i)
            end do
            !        update the matrix
            do k = 1, orderp
                do l = 1, orderp
                    !                 psip(k,l)=psip(k,l)+vecd(k)*psip(l,orderpp)/wy(i)**2
                    !                 psip(k,l)=psip(k,l)+vecd(k)*vecd(l)
                    psipd(k, l) = psipd(k, l) + vecd(k)*vecd(l)
                end do
                !                 coeff(k)=coeff(k)+y(i)*psip(k,orderpp)/wy(i)**2
                coeffd(k) = coeffd(k) + y(i)*(vecd(k)/wy(i))
            end do
        end if
    end do
    do k = 1, orderp
        coeff(k) = coeffd(k)
    end do
    psip(1:orderp, 1:orderp) = psipd(1:orderp, 1:orderp)

    if (map .ne. 0.d0) then
        allocate (mat(orderp, orderp), eig(orderp), vec(orderp), work(3*orderp))
        mat(1:orderp, 1:orderp) = psip(1:orderp, 1:orderp)
        call dsyev('V', 'L', orderp, mat, orderp, eig, work, 3*orderp, info)

        iter = 0
        error = 1.d0
        eps = 1d-8
        maxit = 100
        if (eig(1) .le. 0.d0) write (6, *) &
                &' Warning accuracy diag/cond numb. ', eig(1), abs(eig(orderp)/eig(1))
        eig_min = eig(orderp)*dlamch('E')
        do i = 1, orderp
            if (eig(i) .lt. eig_min) eig(i) = eig_min
        end do
        call dgemv('T', orderp, orderp, one, mat, orderp, coeff, 1, zero, vec, 1)
        do while (iter .le. maxit .and. error .gt. eps)
            iter = iter + 1
            gamma = 0.d0
            do i = 1, orderp
                gamma = gamma + beta*eig(i)/(beta*eig(i) + map)
            end do
            do i = 1, orderp
                work(i) = beta*vec(i)/(map + beta*eig(i))
            end do
            call dgemv('N', orderp, orderp, one, mat, orderp, work, 1, zero, coeff, 1)
            error = map
            map = gamma/sum(coeff(1:orderp)**2)
            error = abs(map - error)
            resd = 0.d0
            do i = 1, N
                fx = coeff(1)
                dfx = 0.d0
                xd = x(i)
                if (abs(power) .eq. 1.d0) then
                    xk = 1.d0
                elseif (power .gt. 0) then
                    xk = xd**(power - 1)
                else
                    xk = xd**(power + 1)
                end if
                if (power .gt. 0) then
                    do k = 2, orderp
                        dfx = dfx + (k - 2 + power)*(xk*coeff(k))
                        xk = xk*x(i)
                        fx = fx + coeff(k)*xk
                    end do
                else
                    do k = 2, orderp
                        xk = xk/x(i)
                        fx = fx + coeff(k)*xk
                    end do
                end if
                ypred(i) = fx
                dypred(i) = dfx
                if (wy(i) .ne. 0.d0) resd = resd + ((fx - y(i))/wy(i))**2
            end do
            res = resd
            if (.not. fixbeta) beta = (Ntrue - gamma)/res
            if (map .gt. 0.d0 .and. beta .gt. 0.d0) then
                evidence = orderp*log(map) + Ntrue*log(beta/TWO_PI) - (beta*res + gamma)
            end if
            do i = 1, orderp
                if (beta*eig(i) + map .gt. 0.d0) then
                    evidence = evidence - log(beta*eig(i) + map)
                end if
            end do
            evidence = evidence/2.d0
            !         write(6,*) ' New value of alpha,beta =',iter,map,beta,res,error
        end do
        if (error .gt. eps .or. iter .ge. maxit) &
                &write (6, *) ' ERROR not converged', error, iter
        deallocate (mat, eig, work, vec)
        !        stop
    end if

    !        Now solve the linear system

    if (map .eq. 0.d0) then
#ifdef __DDIAG
        call dgetrf(orderp, orderp, psipd, orderp, ipsip, info)
#else
        call dgetrf(orderp, orderp, psip, orderp, ipsip, info)
#endif
        if (info .eq. 0) then
#ifdef __DDIAG
            call dgetrs('N', orderp, 1, psipd, orderp, ipsip, coeffd, orderp, info)
            coeff(1:orderp) = coeffd(1:orderp)
#else
            call dgetrs('N', orderp, 1, psip, orderp, ipsip, coeff, orderp, info)
#endif
        else
            write (6, *) ' There is depdendency !!! '
        end if
        resd = 0.d0
        do i = 1, N
#ifdef __DDIAG
            fx = coeffd(1)
#else
            fx = coeff(1)
#endif
            dfx = 0.d0
            xd = x(i)
            if (abs(power) .eq. 1.d0) then
                xk = 1.d0
            elseif (power .gt. 0) then
                xk = xd**(power - 1)
            else
                xk = xd**(power + 1)
            end if

            if (power .gt. 0) then
                do k = 2, orderp
#ifdef __DDIAG
                    dfx = dfx + (k - 2 + power)*coeffd(k)*xk
#else
                    dfx = dfx + (k - 2 + power)*coeff(k)*xk
#endif
                    xk = xk*x(i)
#ifdef __DDIAG
                    fx = fx + coeffd(k)*xk
#else
                    fx = fx + coeff(k)*xk
#endif
                end do
            else
                do k = 2, orderp
                    xk = xk/x(i)
#ifdef __DDIAG
                    fx = fx + coeffd(k)*xk
#else
                    fx = fx + coeff(k)*xk
#endif
                end do
            end if
            ypred(i) = fx
            dypred(i) = dfx
            if (wy(i) .ne. 0.d0) resd = resd + ((fx - y(i))/wy(i))**2
        end do
        res = resd
    end if
    deallocate (coeffd, vecd, psipd)
    return
end

