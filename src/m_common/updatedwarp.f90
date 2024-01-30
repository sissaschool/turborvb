subroutine updatedwarp(nel, nion, kel, rion, rmu, r, kelelocb, kellogb&
        &, rionelocb, rionlogb, cellb, celllb, force, pulay&
        &, iespbc, warp, power, atom_number, warpmat, niong)
    use allio, only: norm_metric
    use Cell, only: CartesianToCrystal, map, cellscale, metric, dmap, car2cry
    implicit none

    ! argument variables
    logical, intent(in) :: iespbc, warp
    integer, intent(in) :: nel, nion, niong
    real*8, intent(in) :: power, kel(3, nel), rion(3, nion), atom_number(*)
    real*8, intent(inout) :: rmu(3, max(nel, nion), nion), r(nel, nion), &
                            &kelelocb(3, nel), kellogb(3, nel), rionelocb(3, nion), &
                            &rionlogb(3, nion), force(3, nion), pulay(3, nion), &
                            &cellb(3), celllb(3), warpmat(nion - niong, *)

    ! local variables
    integer i, j, k, n, irefg, ireft
    real*8 sumw, sumdw, wfunc, wderiv, derpul, wder, xmu(3)

    if (iespbc) then
        do k = 1, nion
            if (atom_number(k) .gt. 0) then
                do i = 1, nel
                    rmu(1, i, k) = kel(1, i) - rion(1, k)
                    rmu(2, i, k) = kel(2, i) - rion(2, k)
                    rmu(3, i, k) = kel(3, i) - rion(3, k)
                end do
                call CartesianToCrystal(rmu(1, 1, k), nel)
                do i = 1, nel
                    rmu(1, i, k) = map(rmu(1, i, k), cellscale(1))
                    rmu(2, i, k) = map(rmu(2, i, k), cellscale(2))
                    rmu(3, i, k) = map(rmu(3, i, k), cellscale(3))
                end do
            else
                do i = 1, nion
                    rmu(1, i, k) = rion(1, k) - rion(1, i)
                    rmu(2, i, k) = rion(2, k) - rion(2, i)
                    rmu(3, i, k) = rion(3, k) - rion(3, i)
                end do
                call CartesianToCrystal(rmu(1, 1, k), nion)
                do i = 1, nion
                    rmu(1, i, k) = map(rmu(1, i, k), cellscale(1))
                    rmu(2, i, k) = map(rmu(2, i, k), cellscale(2))
                    rmu(3, i, k) = map(rmu(3, i, k), cellscale(3))
                end do
            end if
        end do
    else
        do k = 1, nion
            if (atom_number(k) .gt. 0) then
                do i = 1, nel
                    rmu(1, i, k) = kel(1, i) - rion(1, k)
                    rmu(2, i, k) = kel(2, i) - rion(2, k)
                    rmu(3, i, k) = kel(3, i) - rion(3, k)
                end do
            else
                do i = 1, nion
                    rmu(1, i, k) = rion(1, k) - rion(1, i)
                    rmu(2, i, k) = rion(2, k) - rion(2, i)
                    rmu(3, i, k) = rion(3, k) - rion(3, i)
                end do
            end if
        end do
    end if

    if (iespbc) then
        do k = 1, nion
            do i = 1, nel
                r(i, k) = max(norm_metric(rmu(1, i, k), metric), 1d-9)
            end do
        end do
    else
        do k = 1, nion
            do i = 1, nel
                r(i, k) = max(dsqrt(sum(rmu(:, i, k)**2)), 1d-9)
            end do
        end do
    end if

    if (iespbc) then
        do k = 1, nion
            if (atom_number(k) .gt. 0) then
                do i = 1, nel
                    rmu(1, i, k) = kel(1, i) - rion(1, k)
                    rmu(2, i, k) = kel(2, i) - rion(2, k)
                    rmu(3, i, k) = kel(3, i) - rion(3, k)
                end do
                call CartesianToCrystal(rmu(1, 1, k), nel)
                do i = 1, nel
                    xmu(1) = map(rmu(1, i, k), cellscale(1))
                    xmu(2) = map(rmu(2, i, k), cellscale(2))
                    xmu(3) = map(rmu(3, i, k), cellscale(3))
                    rmu(1, i, k) = (metric(1, 1)*xmu(1) + metric(1, 2)*xmu(2) + metric(1, 3)*xmu(3)) &
                                   *dmap(rmu(1, i, k), cellscale(1))
                    rmu(2, i, k) = (metric(2, 1)*xmu(1) + metric(2, 2)*xmu(2) + metric(2, 3)*xmu(3)) &
                                   *dmap(rmu(2, i, k), cellscale(2))
                    rmu(3, i, k) = (metric(3, 1)*xmu(1) + metric(3, 2)*xmu(2) + metric(3, 3)*xmu(3)) &
                                   *dmap(rmu(3, i, k), cellscale(3))
!  HERE rmu is the derivative of r  vs r_cell times r
!    chain rule for r_cell = car2cry x r_physical
!  dr/dr_phisical = dr/dr_cell x dr_cell/dr_physical
                    xmu(:) = rmu(:, i, k)
                    rmu(1, i, k) = xmu(1)*car2cry(1, 1) + xmu(2)*car2cry(2, 1) + xmu(3)*car2cry(3, 1)
                    rmu(2, i, k) = xmu(1)*car2cry(1, 2) + xmu(2)*car2cry(2, 2) + xmu(3)*car2cry(3, 2)
                    rmu(3, i, k) = xmu(1)*car2cry(1, 3) + xmu(2)*car2cry(2, 3) + xmu(3)*car2cry(3, 3)
                end do
            else
                do i = 1, nion
                    rmu(1, i, k) = rion(1, k) - rion(1, i)
                    rmu(2, i, k) = rion(2, k) - rion(2, i)
                    rmu(3, i, k) = rion(3, k) - rion(3, i)
                end do
                call CartesianToCrystal(rmu(1, 1, k), nion)
                do i = 1, nion
                    xmu(1) = map(rmu(1, i, k), cellscale(1))
                    xmu(2) = map(rmu(2, i, k), cellscale(2))
                    xmu(3) = map(rmu(3, i, k), cellscale(3))
                    rmu(1, i, k) = (metric(1, 1)*xmu(1) + metric(1, 2)*xmu(2) + metric(1, 3)*xmu(3)) &
                                   *dmap(rmu(1, i, k), cellscale(1))
                    rmu(2, i, k) = (metric(2, 1)*xmu(1) + metric(2, 2)*xmu(2) + metric(2, 3)*xmu(3)) &
                                   *dmap(rmu(2, i, k), cellscale(2))
                    rmu(3, i, k) = (metric(3, 1)*xmu(1) + metric(3, 2)*xmu(2) + metric(3, 3)*xmu(3)) &
                                   *dmap(rmu(3, i, k), cellscale(3))
!  HERE rmu is the derivative of r  vs r_cell times r
!    chain rule for r_cell = car2cry x r_physical
!  dr/dr_phisical = dr/dr_cell x dr_cell/dr_physical
                    xmu(:) = rmu(:, i, k)
                    rmu(1, i, k) = xmu(1)*car2cry(1, 1) + xmu(2)*car2cry(2, 1) + xmu(3)*car2cry(3, 1)
                    rmu(2, i, k) = xmu(1)*car2cry(1, 2) + xmu(2)*car2cry(2, 2) + xmu(3)*car2cry(3, 2)
                    rmu(3, i, k) = xmu(1)*car2cry(1, 3) + xmu(2)*car2cry(2, 3) + xmu(3)*car2cry(3, 3)
                end do
            end if
        end do
    end if
    ! NB Obviously in the calculation of pressures all electron and ion
    ! positions have to be put in the same box with coordinates
    ! | r_i | < cellscale(i)/2, i=1,2,3.

    force = rionelocb
    pulay = rionlogb

    if (warp) then
        do n = 1, nel
            sumw = 0.d0
            if (power .ne. 0.d0) then
                do k = 1, nion
                    if (atom_number(k) .gt. 0.d0) then
                        wfunc = 1.d0/r(n, k)**power
                        sumw = sumw + wfunc
                    end if
                end do
            end if
            do j = 1, 3
                sumdw = 0.d0
                if (power .ne. 0.d0) then
                    do k = 1, nion
                        !          taken a factor 1/2 of the Jacobian into account
                        if (atom_number(k) .gt. 0.d0) then
                            wderiv = -power/2.d0/r(n, k)**(power + 2)*rmu(j, n, k)
                            sumdw = sumdw + wderiv
                        end if
                    end do
                end if
                do i = 1, nion
                    !          taken a factor 1/2 of the Jacobian into account
                    !          wderiv=-power/2.d0/r(n,i)**(power+2)*rmu(j,n,i)
                    !          wfunc=1.d0/r(n,i)**power

                    if (atom_number(i) .gt. 0.d0) then

                        if (power .eq. 0.d0) then
                            wder = 1.d0
                            derpul = 0.d0
                        else
                            wderiv = -(power/2.d0)/r(n, i)**(power + 2)*rmu(j, n, i)
                            wfunc = 1.d0/r(n, i)**power
                            !          Jacobian contribution
                            derpul = wderiv/sumw - wfunc/sumw**2*sumdw
                            wder = wfunc/sumw
                        end if
                        !          warp contribution
                        pulay(j, i) = pulay(j, i) + wder*kellogb(j, n) + derpul
                        force(j, i) = force(j, i) + wder*kelelocb(j, n)
                    end if
                end do
            end do
        end do
        irefg = 0
        do n = 1, nion
            if (atom_number(n) .le. 0) then
                irefg = irefg + 1
                sumw = 0.d0
                if (power .ne. 0.d0) then
                    do k = 1, nion
                        if (atom_number(k) .gt. 0.d0) then
                            wfunc = 1.d0/r(k, n)**power
                            sumw = sumw + wfunc
                        end if
                    end do
                end if
                do j = 1, 3
                    ireft = 0
                    do i = 1, nion
                        !          taken a factor 1/2 of the Jacobian into account
                        !          wderiv=-power/2.d0/r(n,i)**(power+2)*rmu(j,n,i)
                        !          wfunc=1.d0/r(n,i)**power

                        if (atom_number(i) .gt. 0.d0) then
                            ireft = ireft + 1
                            if (power .ne. 0.d0) then
                                wfunc = 1.d0/r(i, n)**power
                                wder = wfunc/sumw
                            else
                                wder = 1.d0
                            end if
                            !          warp contribution
                            pulay(j, i) = pulay(j, i) + wder*rionlogb(j, n)
                            force(j, i) = force(j, i) + wder*rionelocb(j, n)
                            warpmat(ireft, irefg) = wder
                        end if
                    end do
                end do
            end if
        end do
    end if ! endif warp
    if (iespbc) then
        !     call ApplyPBC(rmu,nion*nel)
        !     press=cellb
        !     press_pulay=celllb
        rmu(:, 1:nion, 1) = rion(:, 1:nion)
        call CartesianToCrystal(rmu, nion)
!       call ApplyPBC(rmu, nion)

        do k = 1, nion
            cellb(1) = cellb(1) + rionelocb(1, k)*rmu(1, k, 1)/cellscale(1)
            celllb(1) = celllb(1) + rionlogb(1, k)*rmu(1, k, 1)/cellscale(1)

            cellb(2) = cellb(2) + rionelocb(2, k)*rmu(2, k, 1)/cellscale(2)
            celllb(2) = celllb(2) + rionlogb(2, k)*rmu(2, k, 1)/cellscale(2)

            cellb(3) = cellb(3) + rionelocb(3, k)*rmu(3, k, 1)/cellscale(3)
            celllb(3) = celllb(3) + rionlogb(3, k)*rmu(3, k, 1)/cellscale(3)
        end do
        rmu(:, 1:nel, 1) = kel(:, 1:nel)
!       call ApplyPBC(rmu, nel)
        call CartesianToCrystal(rmu, nel)

        do k = 1, nel
            cellb(1) = cellb(1) + kelelocb(1, k)*rmu(1, k, 1)/cellscale(1)
            celllb(1) = celllb(1) + kellogb(1, k)*rmu(1, k, 1)/cellscale(1)

            cellb(2) = cellb(2) + kelelocb(2, k)*rmu(2, k, 1)/cellscale(2)
            celllb(2) = celllb(2) + kellogb(2, k)*rmu(2, k, 1)/cellscale(2)

            cellb(3) = cellb(3) + kelelocb(3, k)*rmu(3, k, 1)/cellscale(3)
            celllb(3) = celllb(3) + kellogb(3, k)*rmu(3, k, 1)/cellscale(3)
        end do
    end if

    return
end
