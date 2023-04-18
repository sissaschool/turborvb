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

!******************************
! by E.Coccia (4/1/11)        !
!******************************

! Subroutines:
!
!   - forces_interpolate
!   - forces_evaluate
!   - ext_force
!   - deallocate_forces
!   - vdw_force
!   - force_capping
!   - force_angle
!   - force_dihed
!   - force_improper
!   - cos_mult_der
!   - rf_bond
!   - rf_angle
!   - rf_dihe
!   - rf_dimp

! SUBROUTINE FORCES_INTERPOLATE:
! 3-dimensional (x, y e z) interpolation
! of the external potential for the calculation
! of the forces

subroutine forces_interpolate()

    use ext_forces
    use bspline
    use allio, only: rank, nion, ieskin, nel
    use extpot, only: n_x, n_y, n_z, pot, xdata, ydata, zdata, link_atom
    use splines, only: nxcoef, nycoef, nzcoef
    use van_der_waals, only: vdw

    implicit none

    f_kx = 5
    f_ky = 5
    f_kz = 5

    f_nx = n_x + f_kx
    f_ny = n_y + f_ky
    f_nz = n_z + f_kz

    allocate (f_coef(n_x, n_y, n_z))
    allocate (f_xknot(f_nx), f_yknot(f_ny), f_zknot(f_nz))
    allocate (forcext(3, nion), forcext_el(3, nel))
    if (vdw) allocate (force_vdw(3, nion))

    !generate knots
    call dbsnak(n_x, xdata, f_kx, f_xknot)
    call dbsnak(n_y, ydata, f_ky, f_yknot)
    call dbsnak(n_z, zdata, f_kz, f_zknot)

    !interpolate
    if (rank .eq. 0) then
        write (*, *)
        write (*, *) '|******************************************|'
        write (*, *) '|     FORCES FROM THE QMC/MM POTENTIAL     |'
        write (*, *) '|******************************************|'
        write (*, *)
        write (*, '(a50,3i4)') &
                &  ' Interpolating with 3D spline of order (x,y,z) : '  &
                &, f_kx, f_ky, f_kz
        write (*, *) ''
    end if

    ! 3D interpolation
    call dbs3in(n_x, xdata, n_y, ydata, n_z, zdata, pot, &
            &            n_x, n_y, f_kx, f_ky, f_kz, f_xknot, f_yknot, f_zknot, &
            &            f_coef)

    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) ' Interpolation done.'
        write (*, *) ''
        write (*, *) ''
        write (*, *) '|******************************************|'
        write (*, *) '|     FORCES FROM THE QMC/MM POTENTIAL     |'
        write (*, *) '|******************************************|'
        write (*, *)
    end if
    nxcoef = n_x
    nycoef = n_y
    nzcoef = n_z

    ! by E. Coccia (17/5/11): forces from QM/MM region
    if (link_atom) then
        allocate (mm_f_theta(3, nion), mm_f_dihed(3, nion), mm_f_impr(3, nion))
    end if

    return

end subroutine forces_interpolate

! SUBROUTINE FORCES_EVALUATE:
! evaluation of the x-(y-, z-) derivative
! of the external potential

subroutine forces_evaluate(nxvec, nyvec, nzvec, xvec, yvec, zvec, dedx, dedy, dedz)

    use ext_forces
    use splines, only: nxcoef, nycoef, nzcoef
    use bspline

    implicit none

    integer :: nxvec, nyvec, nzvec
    real*8 :: xvec(nxvec), yvec(nyvec), zvec(nzvec)
    real*8 :: dedx(nxvec, nyvec, nzvec), dedy(nxvec, nyvec, nzvec), dedz(nxvec, nyvec, nzvec)

    ! derivative with respect to x
    call dbs3gd(1, 0, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec, &
            &            f_kx, f_ky, f_kz, f_xknot, f_yknot, f_zknot, nxcoef, &
            &            nycoef, nzcoef, f_coef, dedx, nxvec, nyvec)
    ! derivative with respect to y
    call dbs3gd(0, 1, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec, &
            &            f_kx, f_ky, f_kz, f_xknot, f_yknot, f_zknot, nxcoef, &
            &            nycoef, nzcoef, f_coef, dedy, nxvec, nyvec)
    ! derivative with respect to z
    call dbs3gd(0, 0, 1, nxvec, xvec, nyvec, yvec, nzvec, zvec, &
            &            f_kx, f_ky, f_kz, f_xknot, f_yknot, f_zknot, nxcoef, &
            &            nycoef, nzcoef, f_coef, dedz, nxvec, nyvec)

    return

end subroutine forces_evaluate

!SUBROUTINE EXT_FORCE
! evaluation of the forces
! due to the external potential

subroutine ext_force(nion, rion, zetar, nel, rel)

    use ext_forces, only: forcext, forcext_el
    use extpot, only: x0, delta, n_x, n_y, n_z

    implicit none

    integer, intent(in) :: nion, nel
    real*8, intent(in) :: rion(3, nion), rel(3, nel)
    real*8, intent(in) :: zetar(nion)

    integer :: i, j
    real*8 :: point(3), x1(3)
    real*8 :: dedx(1, 1, 1), dedy(1, 1, 1), dedz(1, 1, 1)

    ! nuclear forces
    do i = 1, nion
        point(:) = rion(:, i)
        call forces_evaluate(1, 1, 1, point(1), point(2), point(3), dedx, dedy, dedz)

        forcext(1, i) = dedx(1, 1, 1)*zetar(i)
        forcext(2, i) = dedy(1, 1, 1)*zetar(i)
        forcext(3, i) = dedz(1, 1, 1)*zetar(i)

    end do

    ! electronic forces
    x1(1) = x0(1) + (n_x - 1)*delta(1)
    x1(2) = x0(2) + (n_y - 1)*delta(2)
    x1(3) = x0(3) + (n_z - 1)*delta(3)
    do i = 1, nel
        point(:) = rel(:, i)
        ! check if the point is out of the potential box
        if (point(1) .gt. x0(1) .and. point(1) .lt. x1(1) .and. &
                &           point(2) .gt. x0(2) .and. point(2) .lt. x1(2) .and. &
                &           point(3) .gt. x0(3) .and. point(3) .lt. x1(3)) then
            call forces_evaluate(1, 1, 1, point(1), point(2), point(3), dedx, dedy, dedz)

            forcext_el(1, i) = -dedx(1, 1, 1)
            forcext_el(2, i) = -dedy(1, 1, 1)
            forcext_el(3, i) = -dedz(1, 1, 1)
        else
            forcext_el(1, i) = 0.d0
            forcext_el(2, i) = 0.d0
            forcext_el(3, i) = 0.d0
        end if
    end do

    return

end subroutine ext_force

!SUBROUTINE DEALLOCATE_FORCES
! deallocate the arrays involved
! in the evaulation of the forces

subroutine deallocate_forces()

    use ext_forces
    use allio, only: rank
    use van_der_waals, only: vdw
    use extpot, only: link_atom

    implicit none

    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) '                    END OF                 '
        write (*, *) ''
        write (*, *) '|******************************************|'
        write (*, *) '|     FORCES FROM THE QMC/MM POTENTIAL     |'
        write (*, *) '|******************************************|'
        write (*, *)
    end if

    deallocate (f_coef)
    deallocate (f_xknot, f_yknot, f_zknot)
    deallocate (forcext, forcext_el)
    if (vdw) deallocate (force_vdw)
    if (link_atom) then
        deallocate (mm_f_theta, mm_f_dihed, mm_f_impr)
    end if

    return

end subroutine deallocate_forces

! Calculation of the vdw forces
subroutine vdw_force(nion, rion)

    use van_der_waals
    use ext_forces, only: force_vdw
    use exc_list
    use extpot, only: link_atom

    implicit none

    integer, intent(in) :: nion
    real*8, intent(in) :: rion(3, nion)

    integer :: i, ii, jj, j, iii, jjj
    real*8 :: xr, yr, zr, r13, r7, flj, tmp, r1, rr
    logical :: compute_vdw, compute_14

    force_vdw = 0.d0
    ! Exclusion list in the case of link atoms
    if (link_atom) then
        do i = 1, nat_tot
            do j = i + 1, nat_tot
                if (grom(i)%qm .and. .not. grom(j)%qm) then
                    if (compute_vdw(i, j)) then
                        iii = grom(i)%it
                        jjj = grom(j)%it
                        xr = rion(1, iii) - coord_nn(1, jjj)
                        yr = rion(2, iii) - coord_nn(2, jjj)
                        zr = rion(3, iii) - coord_nn(3, jjj)
                        rr = sqrt(xr**2 + yr**2 + zr**2)
                        r1 = rr**(-1)
                        r7 = rr**(-7)
                        r13 = rr**(-13)
                        ii = qmc_vdw(iii)
                        jj = nn_vdw(jjj)
                        if (compute_14(i, j)) then
                            tmp = (12.d0*r13*cs12(ii, jj) - 6.d0*r7*cs6(ii, jj))*r1
                        else
                            tmp = (12.d0*r13*c12(ii, jj) - 6.d0*r7*c6(ii, jj))*r1
                        end if
                        flj = tmp*xr
                        force_vdw(1, iii) = force_vdw(1, iii) + flj
                        flj = tmp*yr
                        force_vdw(2, iii) = force_vdw(2, iii) + flj
                        flj = tmp*zr
                        force_vdw(3, iii) = force_vdw(3, iii) + flj
                    end if
                elseif (grom(j)%qm .and. .not. grom(i)%qm) then
                    if (compute_vdw(i, j)) then
                        iii = grom(i)%it
                        jjj = grom(j)%it
                        xr = rion(1, jjj) - coord_nn(1, iii)
                        yr = rion(2, jjj) - coord_nn(2, iii)
                        zr = rion(3, jjj) - coord_nn(3, iii)
                        rr = sqrt(xr**2 + yr**2 + zr**2)
                        r1 = rr**(-1)
                        r7 = rr**(-7)
                        r13 = rr**(-13)
                        ii = qmc_vdw(jjj)
                        jj = nn_vdw(iii)
                        if (compute_14(i, j)) then
                            tmp = (12.d0*r13*cs12(ii, jj) - 6.d0*r7*cs6(ii, jj))*r1
                        else
                            tmp = (12.d0*r13*c12(ii, jj) - 6.d0*r7*c6(ii, jj))*r1
                        end if
                        flj = tmp*xr
                        force_vdw(1, jjj) = force_vdw(1, jjj) + flj
                        flj = tmp*yr
                        force_vdw(2, jjj) = force_vdw(2, jjj) + flj
                        flj = tmp*zr
                        force_vdw(3, jjj) = force_vdw(3, jjj) + flj
                    end if
                end if
            end do
        end do
    else
        do j = 1, nat_nn
            do i = 1, nion
                xr = rion(1, i) - coord_nn(1, j)
                yr = rion(2, i) - coord_nn(2, j)
                zr = rion(3, i) - coord_nn(3, j)
                rr = sqrt(xr**2 + yr**2 + zr**2)
                r1 = rr**(-1)
                r7 = rr**(-7)
                r13 = rr**(-13)
                ii = qmc_vdw(i)
                jj = nn_vdw(j)
                tmp = (12.d0*r13*c12(ii, jj) - 6.d0*r7*c6(ii, jj))*r1
                flj = tmp*xr
                force_vdw(1, i) = force_vdw(1, i) + flj
                flj = tmp*yr
                force_vdw(2, i) = force_vdw(2, i) + flj
                flj = tmp*zr
                force_vdw(3, i) = force_vdw(3, i) + flj
            end do
        end do
    end if

    return

end subroutine vdw_force

!**************************************************************************
subroutine force_capping(rionb)
    !**************************************************************************

    use extpot, only: latoms
    use link_atoms, only: capping, calpha, qm, prt, cap
    use allio, only: nion
    use van_der_waals, only: grom

    implicit none

    real*8, intent(inout), dimension(3, nion) :: rionb
    integer :: i, j, k, m, jj, kk
    real*8 :: factor

    ! rionb(:,:) contains the gradients of the energy
    ! with respect to the atomic positions.
    ! With the definition of R_capping (subroutine r_capping)
    ! the gradient for the capping atom is:
    ! dE/dX_capping = 1/(1-alpha) * dE/dX_QM

    do i = 1, latoms
        !j = capping(i)%cap
        !k = capping(i)%qm
        !jj = grom(j)%is
        !kk = grom(k)%is
        !m = capping(i)%mm

        jj = cap(i)
        kk = qm(i)
        factor = 1.d0 - calpha(i)
        factor = 1.d0/factor
        rionb(:, jj) = factor*rionb(:, kk)
    end do

    return

end subroutine force_capping

!**************************************************
subroutine force_angle()
    !**************************************************

    use allio, only: rion, nion
    use extpot, only: latoms
    use van_der_waals, only: grom, coord_nn, nat_tot
    use link_angle, only: linangle, ntheta
    use ext_forces, only: mm_f_theta

    implicit none

    integer :: i, j, k, ii, jj, kk, iatom
    real*8 :: rijsq, rkjsq
    real*8 :: costh, dvdcos, dvdx, dvdy, dvdz
    real*8 :: dot, ix, iy, iz, kx, ky, kz, tmp
    real*8, dimension(3) :: ai, aj, ak, rij, rkj
    logical :: qm_atom(3)

    ! Derivative of V(theta) with respect to
    ! the QM atoms
    ! See the GROMOS manual (II-17)

    mm_f_theta = 0.d0
    do i = 1, latoms
        do j = 1, ntheta(i)
            iatom = 0
            qm_atom = .false.
            ! BOND ANGLE for I J and K
            !
            !   I     K
            !    \   /
            !     \ /
            !      J
            do k = 1, nat_tot
                ! I ATOM
                if (grom(k)%ind .eq. linangle(i, j)%i) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ai(:) = rion(:, grom(k)%it)
                        qm_atom(1) = .true.
                        ii = grom(k)%it
                    else
                        ! MM atom
                        ai(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! J ATOM
                elseif (grom(k)%ind .eq. linangle(i, j)%j) then
                    if (grom(k)%qm) then
                        ! QM atom
                        aj(:) = rion(:, grom(k)%it)
                        qm_atom(2) = .true.
                        jj = grom(k)%it
                    else
                        ! MM atom
                        aj(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! K ATOM
                elseif (grom(k)%ind .eq. linangle(i, j)%k) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ak(:) = rion(:, grom(k)%it)
                        qm_atom(3) = .true.
                        kk = grom(k)%it
                    else
                        ! MM atom
                        ak(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                end if
                if (iatom .eq. 3) exit
            end do

            ! I-J distance
            rij(:) = ai(:) - aj(:)
            rijsq = sqrt(dot_product(rij, rij))

            ! K-J distance
            rkj(:) = ak(:) - aj(:)
            rkjsq = sqrt(dot_product(rkj, rkj))

            ! COS(THETA)
            dot = dot_product(rij, rkj)
            tmp = rijsq*rkjsq
            costh = dot/tmp
            ! Derivative of V with respect to COS(THETA)
            dvdcos = -linangle(i, j)%ktheta*(costh - linangle(i, j)%thetaeq)

            ix = (rkj(1)/rkjsq - (rij(1)/rijsq)*costh)/rijsq
            iy = (rkj(2)/rkjsq - (rij(2)/rijsq)*costh)/rijsq
            iz = (rkj(3)/rkjsq - (rij(3)/rijsq)*costh)/rijsq
            kx = (rij(1)/rijsq - (rkj(1)/rkjsq)*costh)/rkjsq
            ky = (rij(2)/rijsq - (rkj(2)/rkjsq)*costh)/rkjsq
            kz = (rij(3)/rijsq - (rkj(3)/rkjsq)*costh)/rkjsq
            ! ATOM I is QM
            if (qm_atom(1)) then
                mm_f_theta(1, ii) = mm_f_theta(1, ii) + dvdcos*ix
                mm_f_theta(2, ii) = mm_f_theta(2, ii) + dvdcos*iy
                mm_f_theta(3, ii) = mm_f_theta(3, ii) + dvdcos*iz
            end if
            ! ATOM J IS QM
            if (qm_atom(2)) then
                dvdx = ix + kx
                dvdy = iy + ky
                dvdz = iz + kz
                mm_f_theta(1, jj) = mm_f_theta(1, jj) - dvdcos*dvdx
                mm_f_theta(2, jj) = mm_f_theta(2, jj) - dvdcos*dvdy
                mm_f_theta(3, jj) = mm_f_theta(3, jj) - dvdcos*dvdz
            end if
            ! ATOM K IS QM
            if (qm_atom(3)) then
                mm_f_theta(1, kk) = mm_f_theta(1, kk) + dvdcos*kx
                mm_f_theta(2, kk) = mm_f_theta(2, kk) + dvdcos*ky
                mm_f_theta(3, kk) = mm_f_theta(3, kk) + dvdcos*kz
            end if
        end do
    end do

    return

end subroutine force_angle

!***************************************************************
subroutine force_dihed()
    !***************************************************************

    use allio, only: rion, nion
    use extpot, only: latoms
    use van_der_waals, only: grom, coord_nn, nat_tot
    use link_angle, only: lindhd, nphi, pi
    use ext_forces, only: mm_f_dihed

    implicit none

    integer :: i, j, k, ii, jj, kk, ll, iatom
    real*8 :: phi, cosphi, dvdphi, dvdx, dvdy, dvdz, der
    real*8 :: dum, dot, dot1, app, rkj2, rim(3), rln(3), app1(3), app4
    real*8 :: ix, iy, iz, lx, ly, lz
    real*8, dimension(3) :: ai, aj, ak, al
    real*8 :: rij(3), rkj(3), rkl(3), rmj2, rnk2, rnk(3), rmj(3), rlnsq, rimsq, app2, app3
    logical :: qm_atom(4)

    ! Derivative of V(phi) with respect to
    ! the QM atoms

    mm_f_dihed = 0.d0
    do i = 1, latoms
        do j = 1, nphi(i)
            iatom = 0
            qm_atom = .false.
            do k = 1, nat_tot
                ! I ATOM
                if (grom(k)%ind .eq. lindhd(i, j)%i) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ai(:) = rion(:, grom(k)%it)
                        qm_atom(1) = .true.
                        ii = grom(k)%it
                    else
                        ! MM atom
                        ai(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! J ATOM
                elseif (grom(k)%ind .eq. lindhd(i, j)%j) then
                    if (grom(k)%qm) then
                        ! QM atom
                        aj(:) = rion(:, grom(k)%it)
                        qm_atom(2) = .true.
                        jj = grom(k)%it
                    else
                        ! MM atom
                        aj(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! K ATOM
                elseif (grom(k)%ind .eq. lindhd(i, j)%k) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ak(:) = rion(:, grom(k)%it)
                        qm_atom(3) = .true.
                        kk = grom(k)%it
                    else
                        ! MM atom
                        ak(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! L ATOM
                elseif (grom(k)%ind .eq. lindhd(i, j)%l) then
                    if (grom(k)%qm) then
                        ! QM atom
                        al(:) = rion(:, grom(k)%it)
                        qm_atom(4) = .true.
                        ll = grom(k)%it
                    else
                        ! MM atom
                        al(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                end if
                if (iatom .eq. 4) exit
            end do

            call compute_dihed(ai, aj, ak, al, dum, dum, rij, rkj, rkl, rmj2, rnk2, rnk, rmj)
            !See the GROMOS manual (II-26)
            rkj2 = dot_product(rkj, rkj)
            dot = dot_product(rij, rkj)
            dot1 = dot_product(rkl, rkj)
            app2 = dot/rkj2 - 1.d0
            app3 = dot1/rkj2
            app4 = app2 + 1.d0
            rim(:) = rij(:) - (dot/rkj2)*rkj(:)
            rln(:) = -rkl(:) + (dot1/rkj2)*rkj(:)

            app = dot_product(rim, rln)
            app = app/sqrt(dot_product(rim, rim)*dot_product(rln, rln))
            phi = dacos(app)
            call cross(rkj, rkl, app1)
            phi = sign(1.d0, dot_product(rij, app1))*phi
            phi = modulo(phi, 2.d0*pi)

            ! dV/dcos(m*phi) * dcos(m*phi)/dcos(phi)
            call cos_mult_der(phi, lindhd(i, j)%mult, der)
            dvdphi = -lindhd(i, j)%kphi*lindhd(i, j)%pcos*der

            cosphi = cos(phi)
            rlnsq = sqrt(dot_product(rln, rln))
            rimsq = sqrt(dot_product(rim, rim))
            rlnsq = 1.d0/rlnsq
            rimsq = 1.d0/rimsq

            ix = (rln(1)*rlnsq - rim(1)*rimsq*cosphi)*rimsq
            iy = (rln(2)*rlnsq - rim(2)*rimsq*cosphi)*rimsq
            iz = (rln(3)*rlnsq - rim(3)*rimsq*cosphi)*rimsq
            lx = (rim(1)*rimsq - rln(1)*rlnsq*cosphi)*rlnsq
            ly = (rim(2)*rimsq - rln(2)*rlnsq*cosphi)*rlnsq
            lz = (rim(3)*rimsq - rln(3)*rlnsq*cosphi)*rlnsq

            ! ATOM I is QM
            if (qm_atom(1)) then

                mm_f_dihed(1, ii) = mm_f_dihed(1, ii) + dvdphi*ix
                mm_f_dihed(2, ii) = mm_f_dihed(2, ii) + dvdphi*iy
                mm_f_dihed(3, ii) = mm_f_dihed(3, ii) + dvdphi*iz
            end if
            ! ATOM J IS QM
            if (qm_atom(2)) then

                dvdx = app2*ix - app3*lx
                dvdy = app2*iy - app3*ly
                dvdz = app2*iz - app3*lz

                mm_f_dihed(1, jj) = mm_f_dihed(1, jj) + dvdphi*dvdx
                mm_f_dihed(2, jj) = mm_f_dihed(2, jj) + dvdphi*dvdy
                mm_f_dihed(3, jj) = mm_f_dihed(3, jj) + dvdphi*dvdz
            end if
            ! ATOM K IS QM
            if (qm_atom(3)) then

                dvdx = app4*ix - (app3 - 1.d0)*lx
                dvdy = app4*iy - (app3 - 1.d0)*ly
                dvdz = app4*iz - (app3 - 1.d0)*lz

                mm_f_dihed(1, kk) = mm_f_dihed(1, kk) - dvdphi*dvdx
                mm_f_dihed(2, kk) = mm_f_dihed(2, kk) - dvdphi*dvdy
                mm_f_dihed(3, kk) = mm_f_dihed(3, kk) - dvdphi*dvdz
            end if
            ! ATOM L IS QM
            if (qm_atom(4)) then

                mm_f_dihed(1, ll) = mm_f_dihed(1, ll) + dvdphi*lx
                mm_f_dihed(2, ll) = mm_f_dihed(2, ll) + dvdphi*ly
                mm_f_dihed(3, ll) = mm_f_dihed(3, ll) + dvdphi*lz
            end if
        end do
    end do

    return

end subroutine force_dihed

!***************************************************************
subroutine force_improper()
    !***************************************************************

    use allio, only: rion, nion
    use extpot, only: latoms
    use van_der_waals, only: grom, coord_nn, nat_tot
    use link_angle, only: linimp, nimp, pi
    use ext_forces, only: mm_f_impr

    implicit none

    integer :: i, j, k, ii, jj, kk, ll, iatom
    real*8 :: phi, cosphi, dvdphi, dvdx, dvdy, dvdz, rkj2sq
    real*8 :: rkj2, rmj2, tmp, dot, rnk2, tmp1, dot1
    real*8 :: ix, iy, iz, lx, ly, lz
    real*8, dimension(3) :: ai, aj, ak, al
    real*8 :: rij(3), rkj(3), rkl(3), dp1, rnk(3), rmj(3)
    logical :: qm_atom(4)

    ! Derivative of V(phi)[IMPROPER] with respect to
    ! the QM atoms
    ! See the GROMOS manual (II-22)

    mm_f_impr = 0.d0
    do i = 1, latoms
        do j = 1, nimp(i)
            iatom = 0
            qm_atom = .false.
            do k = 1, nat_tot
                ! I ATOM
                if (grom(k)%ind .eq. linimp(i, j)%i) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ai(:) = rion(:, grom(k)%it)
                        qm_atom(1) = .true.
                        ii = grom(k)%it
                    else
                        ! MM atom
                        ai(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! J ATOM
                elseif (grom(k)%ind .eq. linimp(i, j)%j) then
                    if (grom(k)%qm) then
                        ! QM atom
                        aj(:) = rion(:, grom(k)%it)
                        qm_atom(2) = .true.
                        jj = grom(k)%it
                    else
                        ! MM atom
                        aj(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! K ATOM
                elseif (grom(k)%ind .eq. linimp(i, j)%k) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ak(:) = rion(:, grom(k)%it)
                        qm_atom(3) = .true.
                        kk = grom(k)%it
                    else
                        ! MM atom
                        ak(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! L ATOM
                elseif (grom(k)%ind .eq. linimp(i, j)%l) then
                    if (grom(k)%qm) then
                        ! QM atom
                        al(:) = rion(:, grom(k)%it)
                        qm_atom(4) = .true.
                        ll = grom(k)%it
                    else
                        ! MM atom
                        al(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                end if
                if (iatom .eq. 4) exit
            end do

            call compute_dihed(ai, aj, ak, al, phi, cosphi, rij, rkj, rkl, rmj2, rnk2, rnk, rmj)
            rkj2 = dot_product(rkj, rkj)
            rkj2sq = sqrt(rkj2)
            tmp = rkj2sq/rmj2
            tmp1 = rkj2sq/rnk2
            dot = dot_product(rij, rkj)
            dot1 = dot_product(rkl, rkj)

            ix = tmp*rmj(1)
            iy = tmp*rmj(2)
            iz = tmp*rmj(3)
            lx = tmp1*rnk(1)
            ly = tmp1*rnk(2)
            lz = tmp1*rnk(3)

            ! Derivative of V with respect to PHI
            dvdphi = -linimp(i, j)%kqhi*modulo(phi - linimp(i, j)%qcos, 2.d0*pi)

            ! ATOM I is QM
            if (qm_atom(1)) then

                mm_f_impr(1, ii) = mm_f_impr(1, ii) + dvdphi*ix
                mm_f_impr(2, ii) = mm_f_impr(2, ii) + dvdphi*iy
                mm_f_impr(3, ii) = mm_f_impr(3, ii) + dvdphi*iz
            end if
            ! ATOM J IS QM
            if (qm_atom(2)) then

                dvdx = (dot/rkj2 - 1.d0)*ix - (dot1/rkj2)*lx
                dvdy = (dot/rkj2 - 1.d0)*iy - (dot1/rkj2)*ly
                dvdz = (dot/rkj2 - 1.d0)*iz - (dot1/rkj2)*lz

                mm_f_impr(1, jj) = mm_f_impr(1, jj) + dvdphi*dvdx
                mm_f_impr(2, jj) = mm_f_impr(2, jj) + dvdphi*dvdy
                mm_f_impr(3, jj) = mm_f_impr(3, jj) + dvdphi*dvdz
            end if
            ! ATOM K IS QM
            if (qm_atom(3)) then

                dvdx = (dot/rkj2)*ix - (dot1/rkj2 - 1.d0)*lx
                dvdy = (dot/rkj2)*iy - (dot1/rkj2 - 1.d0)*ly
                dvdz = (dot/rkj2)*iz - (dot1/rkj2 - 1.d0)*lz

                mm_f_impr(1, kk) = mm_f_impr(1, kk) - dvdphi*dvdx
                mm_f_impr(2, kk) = mm_f_impr(2, kk) - dvdphi*dvdy
                mm_f_impr(3, kk) = mm_f_impr(3, kk) - dvdphi*dvdz
            end if
            ! ATOM L IS QM
            if (qm_atom(4)) then

                mm_f_impr(1, ll) = mm_f_impr(1, ll) + dvdphi*lx
                mm_f_impr(2, ll) = mm_f_impr(2, ll) + dvdphi*ly
                mm_f_impr(3, ll) = mm_f_impr(3, ll) + dvdphi*lz
            end if

        end do
    end do

    return

end subroutine force_improper

!*****************************************************************
subroutine cos_mult_der(phi, mult, der)
    !*****************************************************************

    implicit none

    real*8, intent(in) :: phi
    integer, intent(in) :: mult
    real*8, intent(out) :: der

    ! Calculation of dcos(mult*phi)/dcos(phi)
    if (mult .eq. 0) then
        der = 0.d0
    elseif (mult .eq. 1) then
        der = 1.d0
    elseif (mult .eq. 2) then
        der = 4.d0*cos(phi)
    elseif (mult .eq. 3) then
        der = 12.d0*cos(phi)**2 - 3.d0
    elseif (mult .eq. 4) then
        der = 32.d0*cos(phi)**3 - 16.d0*cos(phi)
    elseif (mult .eq. 5) then
        der = 80.d0*cos(phi)**4 - 60.d0*cos(phi)**2 + 5.d0
    elseif (mult .eq. 6) then
        der = 192.d0*cos(phi)**5 - 192.d0*cos(phi)**3 + 36.d0*cos(phi)
    end if

    return

end subroutine cos_mult_der

!******************************************************************
subroutine rf_bond()
    !******************************************************************

    ! Forces for the harmonic distance restraint

    use allio, only: rion
    use cl_restr, only: mm_fact, nbonds, restr_f_bond
    use link_angle, only: cl_bond

    implicit none

    integer :: i, j, kk
    real*8 :: x, y, z, rr, dvdr

    restr_f_bond = 0.d0
    do kk = 1, nbonds

        i = cl_bond(kk)%i
        j = cl_bond(kk)%j

        x = rion(1, i) - rion(1, j)
        y = rion(2, i) - rion(2, j)
        z = rion(3, i) - rion(3, j)

        rr = sqrt(x**2 + y**2 + z**2)

        dvdr = -cl_bond(kk)%kbond*(rr - cl_bond(kk)%req)

        !ATOM I
        restr_f_bond(1, i) = restr_f_bond(1, i) + dvdr*x/rr*mm_fact
        restr_f_bond(2, i) = restr_f_bond(2, i) + dvdr*y/rr*mm_fact
        restr_f_bond(3, i) = restr_f_bond(3, i) + dvdr*z/rr*mm_fact

        !ATOM J
        restr_f_bond(1, j) = -restr_f_bond(1, i)
        restr_f_bond(2, j) = -restr_f_bond(2, i)
        restr_f_bond(3, j) = -restr_f_bond(3, i)

    end do

    return

end subroutine rf_bond

!****************************************************************
subroutine rf_angle()
    !****************************************************************

    ! Forces for the harmonic bond angle restraint

    use allio, only: rion
    use cl_restr, only: mm_fact, nth, restr_f_angle
    use link_angle, only: cl_angle

    implicit none

    integer :: i, j, k, kk
    real*8 :: rijsq, rkjsq
    real*8 :: costh, dvdcos, dvdx, dvdy, dvdz
    real*8 :: dot, ix, iy, iz, kx, ky, kz, tmp
    real*8, dimension(3) :: rij, rkj

    restr_f_angle = 0.d0
    do kk = 1, nth

        i = cl_angle(kk)%i
        j = cl_angle(kk)%j
        k = cl_angle(kk)%k

        ! I-J distance
        rij(:) = rion(:, i) - rion(:, j)
        rijsq = sqrt(dot_product(rij, rij))

        ! K-J distance
        rkj(:) = rion(:, k) - rion(:, j)
        rkjsq = sqrt(dot_product(rkj, rkj))

        ! COS(THETA)
        dot = dot_product(rij, rkj)
        tmp = rijsq*rkjsq
        costh = dot/tmp
        ! Derivative of V with respect to COS(THETA)
        dvdcos = -cl_angle(kk)%ktheta*(costh - cl_angle(kk)%thetaeq)

        ix = (rkj(1)/rkjsq - (rij(1)/rijsq)*costh)/rijsq
        iy = (rkj(2)/rkjsq - (rij(2)/rijsq)*costh)/rijsq
        iz = (rkj(3)/rkjsq - (rij(3)/rijsq)*costh)/rijsq
        kx = (rij(1)/rijsq - (rkj(1)/rkjsq)*costh)/rkjsq
        ky = (rij(2)/rijsq - (rkj(2)/rkjsq)*costh)/rkjsq
        kz = (rij(3)/rijsq - (rkj(3)/rkjsq)*costh)/rkjsq
        ! ATOM I
        restr_f_angle(1, i) = restr_f_angle(1, i) + dvdcos*ix*mm_fact
        restr_f_angle(2, i) = restr_f_angle(2, i) + dvdcos*iy*mm_fact
        restr_f_angle(3, i) = restr_f_angle(3, i) + dvdcos*iz*mm_fact
        ! ATOM J
        dvdx = ix + kx
        dvdy = iy + ky
        dvdz = iz + kz
        restr_f_angle(1, j) = restr_f_angle(1, j) - dvdcos*dvdx*mm_fact
        restr_f_angle(2, j) = restr_f_angle(2, j) - dvdcos*dvdy*mm_fact
        restr_f_angle(3, j) = restr_f_angle(3, j) - dvdcos*dvdz*mm_fact
        ! ATOM K
        restr_f_angle(1, k) = restr_f_angle(1, k) + dvdcos*kx*mm_fact
        restr_f_angle(2, k) = restr_f_angle(2, k) + dvdcos*ky*mm_fact
        restr_f_angle(3, k) = restr_f_angle(3, k) + dvdcos*kz*mm_fact
    end do

    return

end subroutine rf_angle

!******************************************************************
subroutine rf_dihe()
    !******************************************************************

    ! Forces for the harmonic restraint around a proper dihedral

    use allio, only: rion, rank
    use cl_restr, only: mm_fact, ndihe, restr_f_dihe
    use link_angle, only: cl_dihe, pi

    implicit none

    integer :: i, j, k, l, mm
    real*8 :: phi, cosphi, dvdphi, dvdx, dvdy, dvdz, der
    real*8 :: dum, dot, dot1, app, rkj2, rim(3), rln(3), app1(3), app4
    real*8 :: ix, iy, iz, lx, ly, lz
    real*8 :: rij(3), rkj(3), rkl(3), rmj2, rnk2, rnk(3), rmj(3), rlnsq, rimsq, app2, app3

    restr_f_dihe = 0.d0

    do mm = 1, ndihe

        i = cl_dihe(mm)%i
        j = cl_dihe(mm)%j
        k = cl_dihe(mm)%k
        l = cl_dihe(mm)%l

        call compute_dihed(rion(:, i), rion(:, j), rion(:, k), rion(:, l), dum, dum, rij, rkj, rkl, rmj2, rnk2, rnk, rmj)
        !See the GROMOS manual (II-26)
        rkj2 = dot_product(rkj, rkj)
        dot = dot_product(rij, rkj)
        dot1 = dot_product(rkl, rkj)
        app2 = dot/rkj2 - 1.d0
        app3 = dot1/rkj2
        app4 = app2 + 1.d0
        rim(:) = rij(:) - (dot/rkj2)*rkj(:)
        rln(:) = -rkl(:) + (dot1/rkj2)*rkj(:)

        app = dot_product(rim, rln)
        app = app/sqrt(dot_product(rim, rim)*dot_product(rln, rln))
        phi = dacos(app)
        call cross(rkj, rkl, app1)
        phi = sign(1.d0, dot_product(rij, app1))*phi
        phi = modulo(phi, 2.d0*pi)
        cosphi = dcos(phi)
        ! dV/dcos(m*phi) * dcos(m*phi)/dcos(phi)
        !call cos_mult_der(phi, cl_dihe(mm)%mult,der)
        dvdphi = -cl_dihe(mm)%kqhi*(cosphi - cl_dihe(mm)%qcos)

        !cosphi = cos(phi)
        rlnsq = sqrt(dot_product(rln, rln))
        rimsq = sqrt(dot_product(rim, rim))
        rlnsq = 1.d0/rlnsq
        rimsq = 1.d0/rimsq

        ix = (rln(1)*rlnsq - rim(1)*rimsq*cosphi)*rimsq
        iy = (rln(2)*rlnsq - rim(2)*rimsq*cosphi)*rimsq
        iz = (rln(3)*rlnsq - rim(3)*rimsq*cosphi)*rimsq
        lx = (rim(1)*rimsq - rln(1)*rlnsq*cosphi)*rlnsq
        ly = (rim(2)*rimsq - rln(2)*rlnsq*cosphi)*rlnsq
        lz = (rim(3)*rimsq - rln(3)*rlnsq*cosphi)*rlnsq

        ! ATOM I
        restr_f_dihe(1, i) = restr_f_dihe(1, i) + dvdphi*ix*mm_fact
        restr_f_dihe(2, i) = restr_f_dihe(2, i) + dvdphi*iy*mm_fact
        restr_f_dihe(3, i) = restr_f_dihe(3, i) + dvdphi*iz*mm_fact

        ! ATOM J
        dvdx = app2*ix - app3*lx
        dvdy = app2*iy - app3*ly
        dvdz = app2*iz - app3*lz

        restr_f_dihe(1, j) = restr_f_dihe(1, j) + dvdphi*dvdx*mm_fact
        restr_f_dihe(2, j) = restr_f_dihe(2, j) + dvdphi*dvdy*mm_fact
        restr_f_dihe(3, j) = restr_f_dihe(3, j) + dvdphi*dvdz*mm_fact

        ! ATOM K
        dvdx = app4*ix - (app3 - 1.d0)*lx
        dvdy = app4*iy - (app3 - 1.d0)*ly
        dvdz = app4*iz - (app3 - 1.d0)*lz

        restr_f_dihe(1, k) = restr_f_dihe(1, k) - dvdphi*dvdx*mm_fact
        restr_f_dihe(2, k) = restr_f_dihe(2, k) - dvdphi*dvdy*mm_fact
        restr_f_dihe(3, k) = restr_f_dihe(3, k) - dvdphi*dvdz*mm_fact

        ! ATOM L

        restr_f_dihe(1, l) = restr_f_dihe(1, l) + dvdphi*lx*mm_fact
        restr_f_dihe(2, l) = restr_f_dihe(2, l) + dvdphi*ly*mm_fact
        restr_f_dihe(3, l) = restr_f_dihe(3, l) + dvdphi*lz*mm_fact

    end do

    return

end subroutine rf_dihe

!***********************************************************************
subroutine rf_dimp()
    !***********************************************************************

    !  Forces for the harmonic restraint around an improper dihedral

    use allio, only: rion
    use cl_restr, only: mm_fact, ndimp, restr_f_dimp
    use link_angle, only: cl_dimp, pi

    implicit none

    integer :: i, j, k, l, mm
    real*8 :: phi, cosphi, dvdphi, dvdx, dvdy, dvdz, rkj2sq
    real*8 :: rkj2, rmj2, tmp, dot, rnk2, tmp1, dot1
    real*8 :: ix, iy, iz, lx, ly, lz
    real*8 :: rij(3), rkj(3), rkl(3), dp1, rnk(3), rmj(3)

    restr_f_dimp = 0.d0

    do mm = 1, ndimp

        i = cl_dimp(mm)%i
        j = cl_dimp(mm)%j
        k = cl_dimp(mm)%k
        l = cl_dimp(mm)%l

        call compute_dihed(rion(:, i), rion(:, j), rion(:, k), rion(:, l), phi, cosphi, rij, rkj, rkl, rmj2, rnk2, rnk, rmj)
        rkj2 = dot_product(rkj, rkj)
        rkj2sq = sqrt(rkj2)
        tmp = rkj2sq/rmj2
        tmp1 = rkj2sq/rnk2
        dot = dot_product(rij, rkj)
        dot1 = dot_product(rkl, rkj)

        ix = tmp*rmj(1)
        iy = tmp*rmj(2)
        iz = tmp*rmj(3)
        lx = tmp1*rnk(1)
        ly = tmp1*rnk(2)
        lz = tmp1*rnk(3)

        ! Derivative of V with respect to PHI
        dvdphi = -cl_dimp(mm)%kqhi*modulo(phi - cl_dimp(mm)%qcos, 2.d0*pi)

        ! ATOM I
        restr_f_dimp(1, i) = restr_f_dimp(1, i) + dvdphi*ix*mm_fact
        restr_f_dimp(2, i) = restr_f_dimp(2, i) + dvdphi*iy*mm_fact
        restr_f_dimp(3, i) = restr_f_dimp(3, i) + dvdphi*iz*mm_fact

        ! ATOM J
        dvdx = (dot/rkj2 - 1.d0)*ix - (dot1/rkj2)*lx
        dvdy = (dot/rkj2 - 1.d0)*iy - (dot1/rkj2)*ly
        dvdz = (dot/rkj2 - 1.d0)*iz - (dot1/rkj2)*lz

        restr_f_dimp(1, j) = restr_f_dimp(1, j) + dvdphi*dvdx*mm_fact
        restr_f_dimp(2, j) = restr_f_dimp(2, j) + dvdphi*dvdy*mm_fact
        restr_f_dimp(3, j) = restr_f_dimp(3, j) + dvdphi*dvdz*mm_fact

        ! ATOM K
        dvdx = (dot/rkj2)*ix - (dot1/rkj2 - 1.d0)*lx
        dvdy = (dot/rkj2)*iy - (dot1/rkj2 - 1.d0)*ly
        dvdz = (dot/rkj2)*iz - (dot1/rkj2 - 1.d0)*lz

        restr_f_dimp(1, k) = restr_f_dimp(1, k) - dvdphi*dvdx*mm_fact
        restr_f_dimp(2, k) = restr_f_dimp(2, k) - dvdphi*dvdy*mm_fact
        restr_f_dimp(3, k) = restr_f_dimp(3, k) - dvdphi*dvdz*mm_fact

        ! ATOM L
        restr_f_dimp(1, l) = restr_f_dimp(1, l) + dvdphi*lx*mm_fact
        restr_f_dimp(2, l) = restr_f_dimp(2, l) + dvdphi*ly*mm_fact
        restr_f_dimp(3, l) = restr_f_dimp(3, l) + dvdphi*lz*mm_fact

    end do

    return

end subroutine rf_dimp
