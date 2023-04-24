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

!****************************************
! by E. Coccia and L. Guidoni (13/1/11) !
!****************************************

! Subroutines:
!
! - extpot_read
! - interpolate
! - evaluate
! - extpot_ene
! - extpot_final
! - deallocate_extpot
! - ion_ene
! - ion_final
! - vdw_read
! - deallocate_vdw
! - deallocate_link
! - vdw_ene
! - vdw_final
! - r_capping
! - link_read
! - mm_angle
! - mm_dihed
! - cos_mult
! - cross
! - mm_improper
! - compute_dihed
! - exclusion_list
! - compute_vdw
! - compute_14
! - elist
! - restr_read
! - deallocate_restr
! - r_bond
! - r_angle
! - r_dihe
! - r_dimp
! - write_traj

!*********************************************************************
subroutine extpot_read
    !*********************************************************************

    use extpot
    use allio, only: rank

    implicit none
    integer :: i, ix, iy, iz
    real*8 :: tmp(3)

    filename_cube = 'potential.cube'
    open (900, file=filename_cube, form='formatted', status='old')
    read (900, '(a80)') title_cube(1)
    read (900, '(a80)') title_cube(2)
    read (900, *) n_atoms, x0
    read (900, *) n_x, delta(1), tmp(2), tmp(3)
    read (900, *) n_y, tmp(1), delta(2), tmp(3)
    read (900, *) n_z, tmp(1), tmp(2), delta(3)
    allocate (x_atom(n_atoms, 3), id_atom(n_atoms))
    allocate (chrg_atom(n_atoms))
    allocate (pot(n_x, n_y, n_z))
    do i = 1, n_atoms
        read (900, *) id_atom(i), chrg_atom(i), x_atom(i, :)
    end do
    do ix = 1, n_x
        do iy = 1, n_y
            read (900, '(6e13.5)') pot(ix, iy, :)
        end do
    end do
    allocate (xdata(n_x), ydata(n_y), zdata(n_z))
    do i = 1, n_x
        xdata(i) = x0(1) + dble(i - 1)*delta(1)
    end do
    do i = 1, n_y
        ydata(i) = x0(2) + dble(i - 1)*delta(2)
    end do
    do i = 1, n_z
        zdata(i) = x0(3) + dble(i - 1)*delta(3)
    end do

    do i = 1, n_x
        iy = (n_y - 1)/2 + 1
        iz = (n_z - 1)/2 + 1
    end do
    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) '|******************************************|'
        write (*, *) '|       EXTERNAL QMC/MM POTENTIAL          |'
        write (*, *) '|******************************************|'
        write (*, *)
        write (*, *) " Reading cube file ...", filename_cube
        write (*, '(a15,3i7)') '  Dimensions  :', n_x, n_y, n_z
        write (*, '(a10,3f12.6)') '  Mesh  : ', delta
        write (*, *) ''
    end if

    call interpolate

    ave = 0.d0; ave2 = 0.d0
    t_ave = 0.d0; t_ave2 = 0.d0
    ave_ion = 0.d0; ave2_ion = 0.d0
    t_ave_ion = 0.d0; t_ave2_ion = 0.d0

    return
end subroutine extpot_read

!*********************************************************************
subroutine interpolate
    !*********************************************************************
    ! set up the interpolation parameters of the splines

    use extpot
    use splines
    use bspline
    use allio, only: rank

    implicit none

    kxord = 5
    kyord = 5
    kzord = 5

    nxknot = n_x + kxord
    nyknot = n_y + kyord
    nzknot = n_z + kzord

    allocate (bscoef(n_x, n_y, n_z))
    allocate (xknot(nxknot), yknot(nyknot), zknot(nzknot))

    !        generate knots

    call dbsnak(n_x, xdata, kxord, xknot)
    call dbsnak(n_y, ydata, kyord, yknot)
    call dbsnak(n_z, zdata, kzord, zknot)

    !        interpolate
    if (rank .eq. 0) then
        write (*, *)
        write (*, '(a50,3i4)') &
                &  ' Interpolating with 3D spline of order (x,y,z) : '  &
                &, kxord, kyord, kzord
        write (*, *) ''
    end if

    call dbs3in(n_x, xdata, n_y, ydata, n_z, zdata, pot, &
            &            n_x, n_y, kxord, kyord, kzord, xknot, yknot, zknot, &
            &            bscoef)

    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) ' Interpolation done.'
        write (*, *) ''
        write (*, *) '|******************************************|'
        write (*, *) '|       EXTERNAL QMC/MM POTENTIAL          |'
        write (*, *) '|******************************************|'
        write (*, *)
    end if
    nxcoef = n_x
    nycoef = n_y
    nzcoef = n_z

    return
end subroutine interpolate

!*********************************************************************
subroutine evaluate(nxvec, nyvec, nzvec, xvec, yvec, zvec, value)
    !*********************************************************************
    ! evaluate the function on the points

    use extpot
    use splines
    use bspline

    implicit none
    integer :: nxvec, nyvec, nzvec
    integer :: i, j, k
    real*8 :: f
    real*8 :: xvec(nxvec), yvec(nyvec), zvec(nzvec)
    real*8 :: value(nxvec, nyvec, nzvec)

    call dbs3gd(0, 0, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec, &
            &            kxord, kyord, kzord, xknot, yknot, zknot, nxcoef, &
            &            nycoef, nzcoef, bscoef, value, nxvec, nyvec)

    return
end subroutine evaluate

!*********************************************************************
subroutine extpot_ene(coord, nelec, ext_ene, sum_ext)
    !*********************************************************************

    use extpot
    use splines
    use vector
    use bspline

    implicit none

    integer :: nelec
    real*8 :: coord(3, nelec), ext_ene(nelec), sum_ext
    integer :: n
    real*8 :: point(3), valuep(1, 1, 1), x1(3)

    x1(1) = x0(1) + (n_x - 1)*delta(1)
    x1(2) = x0(2) + (n_y - 1)*delta(2)
    x1(3) = x0(3) + (n_z - 1)*delta(3)
    ext_ene = 0.d0
    do n = 1, nelec
        point(:) = coord(:, n)
        ! check if the point is out of the potential box
        if (point(1) .gt. x0(1) .and. point(1) .lt. x1(1) .and. &
                &           point(2) .gt. x0(2) .and. point(2) .lt. x1(2) .and. &
                &           point(3) .gt. x0(3) .and. point(3) .lt. x1(3)) then
            call evaluate(1, 1, 1, point(1), point(2), point(3), valuep)
            ext_ene(n) = valuep(1, 1, 1)
        else
            nout = nout + 1
        end if
    end do

    sum_ext = sum(ext_ene)
    ncount = ncount + 1
    ave = ave + sum_ext
    ave2 = ave2 + sum_ext*sum_ext

    return
end subroutine extpot_ene

!*********************************************************************
subroutine extpot_final(nelec)
    !*********************************************************************

    use extpot

    implicit none

    integer, intent(in) :: nelec

    ! by E. Coccia (29/11/10): t_* variables for MPI runs
    t_ave = t_ave/t_ncount
    err_el = dsqrt((t_ave2/t_ncount - t_ave**2)/t_ncount)

    write (*, *) ''
    write (*, *) '|******************************************|'
    write (*, *) '|       EXTERNAL QMC/MM  POTENTIAL         |'
    write (*, *) '|******************************************|'
    write (*, *)
    write (*, *) ' Total number of N-el evaluations:', t_ncount
    write (*, 110) ' Rate of out of box electrons:', &
            &             t_nout/dble(t_ncount*nelec)
    write (*, 120) ' Average QMC/MM electronic energy:', t_ave, '(', err_el, ')'
    write (*, *)

110 format(a31, 1x, 1f10.6)
120 format(a35, 1x, 1f10.6, 2x, 1a, 1f10.6, 1a)

    return
end subroutine extpot_final

!*******************************
subroutine deallocate_extpot()
    !*******************************

    use extpot
    use splines, only: bscoef, xknot, yknot, zknot

    ! by E. Coccia (30/12/10): deallocate arrays (all the procs!)
    deallocate (bscoef)
    deallocate (xknot, yknot, zknot)

    deallocate (pot, xdata, ydata, zdata)
    deallocate (x_atom, id_atom, chrg_atom)

    return

end subroutine deallocate_extpot

!**************************************************
subroutine ion_ene(rion, zeta, nion, ext_ion, sum_ion)
    !**************************************************

    use extpot
    use splines
    use vector
    use bspline

    implicit none

    integer :: nion
    real*8 :: rion(3, nion), ext_ion(nion), zeta(nion), sum_ion
    integer :: n, j
    real*8 :: point(3), valuep(1, 1, 1), x1(3)

    ext_ion = 0.d0
    do n = 1, nion
        point(:) = rion(:, n)
        call evaluate(1, 1, 1, point(1), point(2), point(3), valuep)
        ext_ion(n) = -valuep(1, 1, 1)*zeta(n)
    end do

    ! Accumulating the nuclear potential
    sum_ion = sum(ext_ion)
    ncount_ion = ncount_ion + 1
    ave_ion = ave_ion + sum_ion
    ave2_ion = ave2_ion + sum_ion*sum_ion

    return

end subroutine ion_ene

!*******************************************************************
subroutine ion_final(nion)
    !*********************************************************************

    use extpot
    use allio, only: idyn
    use van_der_waals, only: vdw

    implicit none

    integer, intent(in) :: nion

    real*8 :: app

    ! t_* variables for MPI runs
    t_ave_ion = t_ave_ion/t_ncount_ion
    ! by E. Coccia (14/1/11): the QMC/MM nuclear term is just
    ! a constant in the case of wf optimization (or VMC and DMC)
    if (idyn .gt. 0) then
        err_ion = dsqrt((t_ave2_ion/t_ncount_ion - t_ave_ion**2)/t_ncount_ion)
    else
        err_ion = 0.d0
    end if

    ! by E. Coccia (19/1/11): to avoid NaN (in general)
    app = t_ave2_ion/t_ncount_ion - t_ave_ion**2
    if (app .le. 0.d0) then
        err_ion = 0.d0
    end if

    write (*, *)
    write (*, *) ' Total number of evaluations:', t_ncount_ion
    write (*, 120) 'Average QMC/MM nuclear energy:', t_ave_ion, '(', err_ion, ')'
    write (*, *)
    total_ave = t_ave + t_ave_ion
    total_err = dsqrt(err_el**2 + err_ion**2)

    if (.not. vdw) then
        write (*, 130) ' Total (N+e) QMC/MM energy:', total_ave, '(', total_err, ')'
        write (*, *)
        write (*, *) '|******************************************|'
        write (*, *) '|         EXTERNAL QMC/MM POTENTIAL        |'
        write (*, *) '|******************************************|'
        write (*, *) ''
    end if

120 format(a33, 1x, 1f10.6, 2x, 1a, 1f10.6, 1a)
130 format(a28, 1x, 1f10.6, 2x, 1a, 1f10.6, 1a)

    return
end subroutine ion_final

!*******************************************************************
subroutine vdw_read()
    !*******************************************************************

    use van_der_waals
    use allio, only: rank, nion, rion
    use extpot, only: link_atom

    implicit none
    integer :: idum, i, j, jdum, i_err, qmcount, mmcount
    character(3) :: cdum, cdum1
    real*8 :: dum1, dum2, dum3, qdum

    filename_vdw = 'vdw.dat'
    open (901, file=filename_vdw, form='formatted', status='old')

    ! nratt = number of Gromos vdw types
    ! idum  = number of QM atoms (must be equal to nion!)
    read (901, *) nratt, idum
    if (idum .ne. nion) then
        if (rank .eq. 0) then
            write (*, *) 'ERROR: the number QM atoms read  must be equal to nion!', idum, '.ne.', nion
        end if
#ifdef PARALLEL
        call mpi_finalize(i_err)
        stop
#else
        stop
#endif
    end if

    ! read c12 and c6 parameters
    ! They depend on the Gromos vdw type
    allocate (c12(nratt, nratt), c6(nratt, nratt))
    ! vdW parameters for 14 interactions
    if (link_atom) then
        allocate (cs12(nratt, nratt), cs6(nratt, nratt))
        do i = 1, nratt
            do j = 1, i
                read (901, *) idum, idum, c12(j, i), c6(j, i), cs12(j, i), cs6(j, i)
                c12(i, j) = c12(j, i)
                c6(i, j) = c6(j, i)
                cs12(i, j) = cs12(j, i)
                cs6(i, j) = cs6(j, i)
            end do
        end do
        read (901, *)
    else
        do i = 1, nratt
            do j = 1, i
                read (901, *) idum, idum, c12(j, i), c6(j, i), dum1, dum1
                c12(i, j) = c12(j, i)
                c6(i, j) = c6(j, i)
            end do
        end do
        read (901, *)
    end if

    nat_tot = 0
    do
        read (901, *, end=600) cdum
        nat_tot = nat_tot + 1
    end do
    ! nat_tot = number of QM + NN atoms
    ! nat_nn  = number of NN atoms
600 nat_nn = nat_tot - nion

    allocate (coord_nn(3, nat_nn), nn_vdw(nat_nn))

    rewind (901)
    read (901, *) idum, idum
    do i = 1, nratt
        do j = 1, i
            read (901, *) idum, idum, dum1, dum1
        end do
    end do

    allocate (qmc_vdw(nion), grom(nat_tot), cpmd(nat_tot))
    grom(:)%qm = .false.
    coord_nn = 0.d0
    qmcount = 0
    mmcount = 0
    ! Dummy array for the atom type ("QUA" or "PRT")
    ! qmc_vdw -> associates the Gromos vdw type to the QMC atom
    ! nn_vdw -> associates the Gromos vdw type to the NN atom
    ! grom%ind -> Gromos index for QM or NN atoms
    ! cpmd%ind -> CPMD index for QM or NN atoms

    do i = 1, nat_tot
        !read(901,*) cdum, grom(i)%ind, cpmd(i), qmc_vdw(qmcount), dum1, dum2, dum3, qdum
        read (901, *) cdum, cdum1, grom(i)%ind, cpmd(i)%ind, idum, dum1, dum2, dum3, qdum
        if (cdum .eq. 'QUA') then
            qmcount = qmcount + 1
            qmc_vdw(qmcount) = idum
            grom(i)%qm = .true.
            grom(i)%it = qmcount
            cpmd(i)%it = qmcount
        elseif (cdum .eq. 'PRT') then
            mmcount = mmcount + 1
            nn_vdw(mmcount) = idum
            grom(i)%it = mmcount
            cpmd(i)%it = mmcount
            coord_nn(1, mmcount) = dum1
            coord_nn(2, mmcount) = dum2
            coord_nn(3, mmcount) = dum3
        end if
    end do

    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) 'Adding the van der Waals contribution'
        write (*, *) 'Number of QM atoms:', nion
        write (*, *) 'Number of NN atoms:', nat_nn
        write (*, *) ''
    end if

    close (901)

    return

end subroutine vdw_read

!*******************************************************************
subroutine deallocate_vdw()
    !*******************************************************************

    use van_der_waals
    use extpot, only: link_atom

    deallocate (coord_nn)
    deallocate (c12, c6)
    deallocate (qmc_vdw, nn_vdw)
    deallocate (grom, cpmd)
    if (link_atom) deallocate (cs12, cs6)

    return

end subroutine deallocate_vdw

!*****************************************************************
subroutine deallocate_link()
    !*****************************************************************

    use link_atoms
    use link_angle
    use exc_list

    deallocate (capping)
    deallocate (log_cap)
    deallocate (linangle, ntheta)
    deallocate (lindhd, nphi)
    deallocate (linimp, nimp)
    deallocate (exc, exc_at)
    deallocate (exc14, exc_14)
    deallocate (qm, cap, prt)

    return

end subroutine deallocate_link

!*******************************************************************
subroutine vdw_ene(rion, nion, sum_vdw)
    !*******************************************************************
    ! Calculation of the vdW interaction between QMC and NN (classical) atoms
    ! Lennard-Jones formula is used:
    ! V_QMC_NN = \sum_QMC \sum_NN [ C12_QMC_NN/(R_QMC_NN)**12  -  C6_QMC_NN/(R_QMC_NN)**6 ]

    use van_der_waals
    use exc_list
    use extpot, only: link_atom

    implicit none

    integer, intent(in) :: nion
    real*8, intent(in) :: rion(3, nion)
    real*8, intent(inout) :: sum_vdw

    integer i, j, ii, jj, iii, jjj, k
    real*8 :: rr, xr, yr, zr, lj, r6, r12
    logical :: compute_vdw, compute_14

    sum_vdw = 0.d0
    ! Exclusion list in the case of link atoms
    if (link_atom) then
        do i = 1, nat_tot
            do j = 1 + i, nat_tot
                if (grom(i)%qm .and. .not. grom(j)%qm) then
                    if (compute_vdw(i, j)) then
                        iii = grom(i)%it
                        jjj = grom(j)%it
                        xr = (rion(1, iii) - coord_nn(1, jjj))**2
                        yr = (rion(2, iii) - coord_nn(2, jjj))**2
                        zr = (rion(3, iii) - coord_nn(3, jjj))**2
                        rr = sqrt(xr + yr + zr)
                        r6 = rr**(-6)
                        r12 = rr**(-12)
                        ii = qmc_vdw(iii)
                        jj = nn_vdw(jjj)
                        if (compute_14(i, j)) then
                            lj = r12*cs12(ii, jj) - r6*cs6(ii, jj)
                        else
                            lj = r12*c12(ii, jj) - r6*c6(ii, jj)
                        end if
                        sum_vdw = sum_vdw + lj
                    end if
                elseif (grom(j)%qm .and. .not. grom(i)%qm) then
                    if (compute_vdw(i, j)) then
                        iii = grom(i)%it
                        jjj = grom(j)%it
                        xr = (rion(1, jjj) - coord_nn(1, iii))**2
                        yr = (rion(2, jjj) - coord_nn(2, iii))**2
                        zr = (rion(3, jjj) - coord_nn(3, iii))**2
                        rr = sqrt(xr + yr + zr)
                        r6 = rr**(-6)
                        r12 = rr**(-12)
                        ii = qmc_vdw(jjj)
                        jj = nn_vdw(iii)
                        if (compute_14(i, j)) then
                            lj = r12*cs12(ii, jj) - r6*cs6(ii, jj)
                        else
                            lj = r12*c12(ii, jj) - r6*c6(ii, jj)
                        end if
                        sum_vdw = sum_vdw + lj
                    end if
                end if
            end do
        end do
    else
        do i = 1, nion
            do j = 1, nat_nn
                xr = (rion(1, i) - coord_nn(1, j))**2
                yr = (rion(2, i) - coord_nn(2, j))**2
                zr = (rion(3, i) - coord_nn(3, j))**2
                rr = sqrt(xr + yr + zr)
                r6 = rr**(-6)
                r12 = rr**(-12)
                ii = qmc_vdw(i)
                jj = nn_vdw(j)
                lj = r12*c12(ii, jj) - r6*c6(ii, jj)
                sum_vdw = sum_vdw + lj
            end do
        end do
    end if

    ! Accumulating the vdw potential
    ncount_vdw = ncount_vdw + 1
    ave_vdw = ave_vdw + sum_vdw
    ave2_vdw = ave2_vdw + sum_vdw*sum_vdw

    return

end subroutine vdw_ene

!*******************************************************************
subroutine vdw_final()
    !*******************************************************************

    use van_der_waals
    use extpot, only: mm_restr, total_ave, err_el, err_ion, link_atom
    use allio, only: idyn
    use tot_angle

    implicit none

    real*8 :: app

    if (.not. mm_restr) then
        ! t_* variables for MPI runs
        t_ave_vdw = t_ave_vdw/t_ncount_vdw

        ! No error in the case of wf optimization (or VMC and DMC)
        if (idyn .gt. 0) then
            err_vdw = dsqrt((t_ave2_vdw/t_ncount_vdw - t_ave_vdw**2)/t_ncount_vdw)
        else
            err_vdw = 0.d0
        end if

        ! by E. Coccia (2/3/11): to avoid NaN (in general)
        app = t_ave2_vdw/t_ncount_vdw - t_ave_vdw**2
        if (app .le. 0.d0) then
            err_vdw = 0.d0
        end if
    end if

    if (link_atom .or. mm_restr) then
        ! t_* variables for MPI runs
        if (mm_restr) t_ave_bond = t_ave_bond/t_ncount_bond
        t_ave_angle = t_ave_angle/t_ncount_angle
        t_ave_dihed = t_ave_dihed/t_ncount_dihed
        t_ave_impr = t_ave_impr/t_ncount_impr

        ! No error in the case of wf optimization (or VMC and DMC)
        if (idyn .gt. 0) then
            if (mm_restr) err_bond = dsqrt((t_ave2_bond/t_ncount_bond - t_ave_bond**2)/t_ncount_bond)
            err_angle = dsqrt((t_ave2_angle/t_ncount_angle - t_ave_angle**2)/t_ncount_angle)
            err_dihed = dsqrt((t_ave2_dihed/t_ncount_dihed - t_ave_dihed**2)/t_ncount_dihed)
            err_impr = dsqrt((t_ave2_impr/t_ncount_impr - t_ave_impr**2)/t_ncount_impr)
        else
            if (mm_restr) err_bond = 0.d0
            err_angle = 0.d0
            err_dihed = 0.d0
            err_impr = 0.d0
        end if

        ! To avoid NaN (in general)
        if (mm_restr) then
            app = t_ave2_bond/t_ncount_bond - t_ave_bond**2
            if (app .le. 0.d0) then
                err_bond = 0.d0
            end if
        end if
        app = t_ave2_angle/t_ncount_angle - t_ave_angle**2
        if (app .le. 0.d0) then
            err_angle = 0.d0
        end if
        app = t_ave2_dihed/t_ncount_dihed - t_ave_dihed**2
        if (app .le. 0.d0) then
            err_dihed = 0.d0
        end if
        app = t_ave2_impr/t_ncount_impr - t_ave_impr**2
        if (app .le. 0.d0) then
            err_impr = 0.d0
        end if
    end if

    if (.not. mm_restr) then
        write (*, *)
        write (*, *) ' Total number of vdW evaluations:', t_ncount_vdw
        write (*, 120) ' Average QMC/MM vdW energy:', t_ave_vdw, '(', err_vdw, ')'
        write (*, *)
    end if
    if (link_atom .or. mm_restr) then
        if (mm_restr) then
            write (*, *) ' Total number of bond distance evaluations:', t_ncount_bond
            write (*, 120) ' Average bond distance energy:', t_ave_bond, '(', err_bond, ')'
            write (*, *)
        end if
        write (*, *) ' Total number of bond angle evaluations:', t_ncount_angle
        write (*, 120) ' Average bond angle energy:', t_ave_angle, '(', err_angle, ')'
        write (*, *)
        write (*, *) ' Total number of proper dihedral evaluations:', t_ncount_dihed
        write (*, 120) ' Average proper dihedral energy:', t_ave_dihed, '(', err_dihed, ')'
        write (*, *)
        write (*, *) ' Total number of improper dihedral evaluations:', t_ncount_impr
        write (*, 120) ' Average improper dihedral energy:', t_ave_impr, '(', err_impr, ')'
        write (*, *)
        if (link_atom) then
            sum_pot = total_ave + t_ave_vdw + t_ave_angle + t_ave_dihed + t_ave_impr
            sum_err = dsqrt(err_el**2 + err_ion**2 + err_vdw**2 + err_angle**2 + err_dihed**2 + err_impr**2)
        elseif (mm_restr) then
            sum_pot = t_ave_bond + t_ave_angle + t_ave_dihed + t_ave_impr
            sum_err = dsqrt(err_bond + err_angle**2 + err_dihed**2 + err_impr**2)
        end if
        if (link_atom) then
            write (*, 130) ' Total QMC/MM (electrostatic + vdW + MM) energy:', sum_pot, '(', sum_err, ')'
            write (*, *)
        elseif (mm_restr) then
            write (*, 130) ' Total QMC/MM (restraints) energy:', sum_pot, '(', sum_err, ')'
        end if

    elseif (vdw .and. .not. link_atom) then
        sum_pot = total_ave + t_ave_vdw
        sum_err = dsqrt(err_el**2 + err_ion**2 + err_vdw**2)
        write (*, 130) ' Total QMC/MM (electrostatic + vdW) energy:', sum_pot, '(', sum_err, ')'
        write (*, *)
    end if
    if (.not. mm_restr) then
        write (*, *) '|******************************************|'
        write (*, *) '|         EXTERNAL QMC/MM POTENTIAL        |'
        write (*, *) '|******************************************|'
        write (*, *) ''
    end if

120 format(a34, 1x, 1f12.6, 2x, 1a, 1f12.6, 1a)
130 format(a48, 1x, 1f12.6, 2x, 1a, 1f12.6, 1a)

    return

end subroutine vdw_final

!*****************************************************************************
subroutine r_capping()
    !*****************************************************************************

    use extpot, only: latoms
    use link_atoms, only: capping, calpha, qm, cap, prt
    use allio, only: rion
    use van_der_waals, only: coord_nn, nat_tot, grom

    implicit none

    integer :: i, j, k, m, jj, kk, ll
    ! The capping atom is not an independent variable
    ! But it is to be on the QM-link/MM-link bond
    ! From Theor. Chem. Acc., vol. 100, 307 (1998)
    ! j -> capping atom
    ! k -> QM-link atom
    ! m -> MM-link atom
    ! R_cap(:) = R_QM(:) + alpha*(R_MM(:) - R_QM(:))
    ! The orginal implementation for alpha is
    ! alpha = (||R_MM(:) - R_QM(:)|| - deltaR)/||R_MM(:) - R_QM(:)||
    ! This gives ||R_MM(:) - R_QM(:)|| = ||R_capping(:) - R_QM(:)|| + DeltaR (with DeltaR constant)
    ! I use a constant value for alhpa:
    !  ||R_MM(:) - R_QM(:)|| = alpha * ||R_capping(:) - R_QM(:)||
    do i = 1, latoms

        !j = capping(i)%cap
        !k = capping(i)%qm
        !m = capping(i)%mm
        !jj = grom(j)%is
        !kk =  grom(k)%is
        !ll = grom(m)%is

        jj = cap(i)
        kk = qm(i)
        ll = prt(i)

        rion(:, jj) = rion(:, kk) + calpha(i)*(coord_nn(:, ll) - rion(:, kk))

    end do

    return

end subroutine r_capping

!********************************************************
subroutine link_read()
    !********************************************************

    use extpot, only: latoms
    use link_atoms, only: capping, calpha
    use link_angle
    use van_der_waals, only: grom, cpmd, nat_tot
    use allio, only: rank

    implicit none

    integer :: i, j, k, idum, kk, icap, iqm, imm
    character(100) :: cdum

    filename_link = 'link.dat'
    open (902, file=filename_link, form='formatted', status='old')

    !Number of capping atoms
    read (902, *) cdum
    read (902, *) latoms
    if (rank .eq. 0) then
        write (*, *) ''
        write (6, *) '|****************************************|'
        write (6, *) '| Link atoms in the QMC/MM framework     |'
        write (6, *) '|****************************************|'
        write (*, *) ''
        write (*, *) 'Number of link atoms:', latoms
        write (*, *) ''
    end if

    allocate (linangle(latoms, maxth), ntheta(latoms))
    allocate (capping(latoms))

    capping(:)%cap = 1
    capping(:)%qm = 1
    capping(:)%mm = 1

    linangle(:, :)%i = 1
    linangle(:, :)%j = 1
    linangle(:, :)%k = 1
    linangle(:, :)%ktheta = 0.d0
    linangle(:, :)%thetaeq = 0.d0

    ! Header of link.dat containing the (Gromos) list
    ! of the capping, QM-link and MM-link atoms
    read (902, *) cdum

    do i = 1, latoms
        read (902, *) icap, iqm, imm
        do j = 1, nat_tot
            if (icap .eq. grom(j)%ind) then
                capping(i)%cap = cpmd(j)%ind
            end if
            if (iqm .eq. grom(j)%ind) then
                capping(i)%qm = cpmd(j)%ind
            end if
            if (imm .eq. grom(j)%ind) then
                capping(i)%mm = cpmd(j)%ind
            end if
        end do
    end do

    if (rank .eq. 0) then
        write (*, *) 'Capping  QM-link  MM-link (Turbo/CPMD list)'
        do i = 1, latoms
            write (*, '(3I8)') capping(i)%cap, capping(i)%qm, capping(i)%mm
        end do
        write (*, *) ''
        write (6, *) '|****************************************|'
        write (6, *) '| Link atoms in the QMC/MM framework     |'
        write (6, *) '|****************************************|'
        write (*, *) ''
        write (*, *)
    end if

    read (902, *) cdum

    ! ANGULAR CONTRIBUTION
    ! For each link atom, there exists a term
    ! 1/2 * K_theta * (cos(theta) - cos(theta_eq))**2
    ! The first part of link.dat is (for each link atom!)
    ! NTHETA
    ! I   J  K    KTHETA    COS(THETA_EQ)
    !   ...        ...     ...     ...              ...
    ! I J and K are Gromos indexes
    ! KTHETA is the force constant
    ! COS(THETA_EQ) is an adimensiona factor

    do i = 1, latoms
        read (902, *) ntheta(i)
        do j = 1, ntheta(i)
            read (902, *) linangle(i, j)%i, linangle(i, j)%j, linangle(i, j)%k, linangle(i, j)%ktheta, linangle(i, j)%thetaeq
        end do
    end do

    ! PROPER DIHEDRAL CONTRIBUTION

    allocate (lindhd(latoms, maxphi), nphi(latoms))

    ! NPHI
    ! I J K L MULT KPHI COS(PHI_EQ)
    ! .. .. .. .. .. .. ..
    ! I J K and L are Gromos indexes
    ! MULT is an integer
    ! KPHI is the force constant
    ! COS(PHI_EQ) is an adimensional  factor
    ! V_dihed = sum_dihed KPHI*(1.d0 + COS(PHI_EQ)*COS(MULT*PHI(I,J,K,L)))

    lindhd(:, :)%i = 1
    lindhd(:, :)%j = 1
    lindhd(:, :)%k = 1
    lindhd(:, :)%l = 1
    lindhd(:, :)%mult = 0
    lindhd(:, :)%kphi = 0.d0
    lindhd(:, :)%pcos = 0.d0

    read (902, *) cdum
    do i = 1, latoms
        read (902, *) nphi(i)
        do j = 1, nphi(i)
            read (902, *) lindhd(i, j)%i, lindhd(i, j)%j, &
                lindhd(i, j)%k, lindhd(i, j)%l, lindhd(i, j)%mult, &
                lindhd(i, j)%kphi, lindhd(i, j)%pcos
        end do
    end do

    ! IMPROPER DIHEDRAL CONTRIBUTION
    ! NQHI
    ! I J K L KQHI COS(QHI_EQ)
    ! .. .. .. .. .. ..
    ! I J K and L are Gromos indexes
    ! KQHI is the force constant
    ! COS(QHI_EQ) is an adimensional factor
    ! V_imp = sum_imp 0.5*KQHI*(COS(QHI) - COS(QHI_EQ) )**2

    allocate (linimp(latoms, maxqhi), nimp(latoms))

    linimp(:, :)%i = 1
    linimp(:, :)%j = 1
    linimp(:, :)%k = 1
    linimp(:, :)%l = 1
    linimp(:, :)%kqhi = 0.d0
    linimp(:, :)%qcos = 0.d0

    read (902, *) cdum
    do i = 1, latoms
        read (902, *) nimp(i)
        do j = 1, nimp(i)
            read (902, *) linimp(i, j)%i, linimp(i, j)%j, linimp(i, j)%k, linimp(i, j)%l, linimp(i, j)%kqhi, linimp(i, j)%qcos
        end do
    end do
    close (902)

    ! Logical array for recognizing capping atoms
    ! and other important definitions regarding capping atoms
    call define_link()

    return

end subroutine link_read

!******************************************************************
subroutine define_link()
    !******************************************************************

    use van_der_waals, only: grom, cpmd, nat_tot
    use link_atoms, only: capping, log_cap, cap, qm, prt
    use allio, only: nion
    use extpot, only: latoms

    implicit none

    integer :: i, iqua, j

    allocate (log_cap(nion))

    log_cap = .false.
    iqua = 0
    do i = 1, nat_tot
        if (grom(i)%qm) then
            iqua = iqua + 1
            do j = 1, latoms
                if (cpmd(i)%ind .eq. capping(j)%cap) then
                    log_cap(grom(i)%it) = .true.
                    exit
                end if
            end do
            if (iqua .eq. nion) exit
        end if
    end do

    ! index for capping anf link atoms
    allocate (cap(latoms), qm(latoms), prt(latoms))
    do j = 1, latoms
        do i = 1, nat_tot
            if (cpmd(i)%ind .eq. capping(j)%cap) then
                cap(j) = cpmd(i)%it
            elseif (cpmd(i)%ind .eq. capping(j)%qm) then
                qm(j) = cpmd(i)%it
            elseif (cpmd(i)%ind .eq. capping(j)%mm) then
                prt(j) = cpmd(i)%it
            end if
        end do
    end do

    return

end subroutine define_link

!******************************************************************
subroutine mm_angle()
    !******************************************************************

    use allio, only: rion, nion
    use extpot, only: latoms
    use van_der_waals, only: coord_nn, grom, nat_tot
    use link_angle, only: mm_pot_theta, linangle, ntheta
    use tot_angle, only: ncount_angle, ave_angle, ave2_angle

    implicit none

    integer :: i, j, k, m, n, iatom
    real*8 :: x, y, z, bij, bjk, x1, y1, z1
    real*8 :: costh
    real*8, dimension(3) :: ai, aj, ak

    ! For each triplet ( QM-link -- MM-link -- MM) or (QM -- QM-link -- MM-link)
    ! the routine calculates
    ! 1/2 * k_theta * (cos(theta) - cos(theta_eq))**2
    ! See the Gromos manual (II-17)

    mm_pot_theta = 0.d0

    do i = 1, latoms
        do j = 1, ntheta(i)
            iatom = 0
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
                    else
                        ! MM atom
                        aj(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                elseif (grom(k)%ind .eq. linangle(i, j)%k) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ak(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        ak(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                end if
                if (iatom .eq. 3) exit
            end do

            ! I-J distance
            x = ai(1) - aj(1)
            y = ai(2) - aj(2)
            z = ai(3) - aj(3)
            bij = sqrt(x**2 + y**2 + z**2)

            ! K-J distance
            x1 = ak(1) - aj(1)
            y1 = ak(2) - aj(2)
            z1 = ak(3) - aj(3)
            bjk = sqrt(x1**2 + y1**2 + z1**2)
            ! COS(THETA)
            costh = (x*x1 + y*y1 + z*z1)/(bij*bjk)
            mm_pot_theta = mm_pot_theta + linangle(i, j)%ktheta*(costh - linangle(i, j)%thetaeq)**2
        end do
    end do
    mm_pot_theta = mm_pot_theta*0.5d0

    !Accumulating the bond angle potential
    ncount_angle = ncount_angle + 1
    ave_angle = ave_angle + mm_pot_theta
    ave2_angle = ave2_angle + mm_pot_theta*mm_pot_theta

    return

end subroutine mm_angle

!******************************************************************
subroutine mm_dihed()
    !******************************************************************

    use extpot, only: latoms
    use allio, only: rion, nion
    use van_der_waals, only: coord_nn, grom, nat_tot
    use link_angle, only: mm_pot_dihed, nphi, lindhd, pi
    use tot_angle, only: ncount_dihed, ave_dihed, ave2_dihed

    implicit none

    integer :: i, j, k, iatom
    real*8, dimension(3) :: ai, aj, ak, al
    real*8 :: phi, cosphi, rij(3), rkj(3), rkl(3), rmj2, rnk2, rnk(3), rmj(3)
    real*8 :: dot, dot1, rkj2, app, app1(3), cosine, rim(3), rln(3), dum

    mm_pot_dihed = 0.d0
    do i = 1, latoms
        do j = 1, nphi(i)
            iatom = 0
            ! PROPER DIHEDRAL ANGLE for I J K and L
            !    I\
            !      \J____K
            !             \
            !              \L
            do k = 1, nat_tot
                !I-th atom
                if (grom(k)%ind .eq. lindhd(i, j)%i) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ai(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        ai(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! J-th atom
                elseif (grom(k)%ind .eq. lindhd(i, j)%j) then
                    if (grom(k)%qm) then
                        ! QM atom
                        aj(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        aj(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! K-th atom
                elseif (grom(k)%ind .eq. lindhd(i, j)%k) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ak(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        ak(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! L-th atom
                elseif (grom(k)%ind .eq. lindhd(i, j)%l) then
                    if (grom(k)%qm) then
                        ! QM atom
                        al(:) = rion(:, grom(k)%it)
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
            rim(:) = rij(:) - (dot/rkj2)*rkj(:)
            rln(:) = -rkl(:) + (dot1/rkj2)*rkj(:)

            app = dot_product(rim, rln)
            app = app/sqrt(dot_product(rim, rim)*dot_product(rln, rln))
            phi = dacos(app)
            call cross(rkj, rkl, app1)
            phi = sign(1.d0, dot_product(rij, app1))*phi
            !phi = modulo(phi,2.d0*pi)

            call cos_mult(phi, lindhd(i, j)%mult, cosine)
            mm_pot_dihed = mm_pot_dihed + lindhd(i, j)%kphi*(1.d0 + lindhd(i, j)%pcos*cosine)

        end do
    end do

    !Accumulating the proper dihedral potential
    ncount_dihed = ncount_dihed + 1
    ave_dihed = ave_dihed + mm_pot_dihed
    ave2_dihed = ave2_dihed + mm_pot_dihed*mm_pot_dihed

    return

end subroutine mm_dihed

!****************************************************************
subroutine cos_mult(phi, mult, cosine)

    implicit none

    real*8, intent(in) :: phi
    integer, intent(in) :: mult
    real*8, intent(out) :: cosine

    ! Calculation of cos(mult*phi)
    if (mult .eq. 0) then
        cosine = 1.d0
    elseif (mult .eq. 1) then
        cosine = cos(phi)
    elseif (mult .eq. 2) then
        cosine = 2.d0*cos(phi)**2 - 1.d0
    elseif (mult .eq. 3) then
        cosine = 4.d0*cos(phi)**3 - 3.d0*cos(phi)
    elseif (mult .eq. 4) then
        cosine = 8.d0*cos(phi)**4 - 8.d0*cos(phi)**2 + 1.d0
    elseif (mult .eq. 5) then
        cosine = 16.d0*cos(phi)**5 - 20.d0*cos(phi)**3 + 5.d0*cos(phi)
    elseif (mult .eq. 6) then
        cosine = 32.d0*cos(phi)**6 - 48.d0*cos(phi)**4 + 18.d0*cos(phi)**2 - 1.d0
    end if

    return

end subroutine cos_mult
!****************************************************************
subroutine cross(a, b, c)
    !****************************************************************

    implicit none

    real*8, intent(in), dimension(3) :: a, b
    real*8, intent(out), dimension(3) :: c

    c(1) = (a(2)*b(3) - a(3)*b(2))
    c(2) = (a(3)*b(1) - a(1)*b(3))
    c(3) = (a(1)*b(2) - a(2)*b(1))

    return

end subroutine cross

!*************************************************************
subroutine mm_improper()
    !*************************************************************

    use extpot, only: latoms
    use allio, only: rion, nion
    use van_der_waals, only: coord_nn, grom, nat_tot
    use link_angle, only: mm_pot_impr, nimp, linimp, pi
    use tot_angle, only: ave_impr, ave2_impr, ncount_impr

    implicit none

    integer :: i, j, k, iatom
    real*8, dimension(3) :: ai, aj, ak, al
    real*8 :: cosphi, phi, rij(3), rkj(3), rkl(3), rmj2, rnk2, rnk(3), rmj(3)

    mm_pot_impr = 0.d0
    do i = 1, latoms
        do j = 1, nimp(i)
            iatom = 0
            do k = 1, nat_tot
                ! I-th atom
                if (grom(k)%ind .eq. linimp(i, j)%i) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ai(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        ai(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! J-th atom
                elseif (grom(k)%ind .eq. linimp(i, j)%j) then
                    if (grom(k)%qm) then
                        ! QM atom
                        aj(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        aj(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! K-th atom
                elseif (grom(k)%ind .eq. linimp(i, j)%k) then
                    if (grom(k)%qm) then
                        ! QM atom
                        ak(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        ak(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                    ! L-th atom
                elseif (grom(k)%ind .eq. linimp(i, j)%l) then
                    if (grom(k)%qm) then
                        ! QM atom
                        al(:) = rion(:, grom(k)%it)
                    else
                        ! MM atom
                        al(:) = coord_nn(:, grom(k)%it)
                    end if
                    iatom = iatom + 1
                end if
                if (iatom .eq. 4) exit
            end do

            call compute_dihed(ai, aj, ak, al, phi, cosphi, rij, rkj, rkl, rmj2, rnk2, rnk, rmj)
            mm_pot_impr = mm_pot_impr + linimp(i, j)%kqhi*modulo(phi - linimp(i, j)%qcos, 2.d0*pi)**2
        end do
    end do
    mm_pot_impr = mm_pot_impr*0.5d0

    !Accumulating the improper dihedral potential
    ncount_impr = ncount_impr + 1
    ave_impr = ave_impr + mm_pot_impr
    ave2_impr = ave2_impr + mm_pot_impr*mm_pot_impr

    return

end subroutine mm_improper

!***************************************************************
subroutine compute_dihed(ai, aj, ak, al, phi, cosphi, rij, rkj, rkl, rmj2, rnk2, rnk, rmj)
    !***************************************************************

    implicit none

    real*8, intent(in) :: ai(3), aj(3), ak(3), al(3)
    real*8, intent(out) :: cosphi, phi, rmj2, rnk2
    real*8, intent(out), dimension(3) :: rij, rkj, rkl, rmj, rnk

    real*8 :: dp1, dp4
    integer :: i

    !See the Gromos manual (II-21)
    do i = 1, 3
        rij(i) = ai(i) - aj(i)
        rkj(i) = ak(i) - aj(i)
        rkl(i) = ak(i) - al(i)
    end do

    call cross(rij, rkj, rmj)
    call cross(rkj, rkl, rnk)

    dp1 = dot_product(rmj, rnk)
    rmj2 = dot_product(rmj, rmj)
    rnk2 = dot_product(rnk, rnk)
    cosphi = dp1/sqrt(rmj2*rnk2)
    dp4 = dot_product(rij, rnk)
    phi = sign(1.d0, dp4)*dacos(cosphi)

    return

end subroutine compute_dihed

!************************************************************
subroutine exclusion_list()
    !************************************************************
    use exc_list
    use van_der_waals, only: nat_tot, grom
    use allio, only: rank

    implicit none

    integer :: i, j, ii, idum1, idum2, jj
    character(50) :: cdum
    logical :: elist

    ! First and second neighbours

    allocate (exc(nat_tot))
    allocate (exc_at(nat_tot, nmax))
    exc_at = 0.d0

    open (702, file='exclusion.dat', form='formatted', status='old')
    read (702, *) cdum
    read (702, *) cdum
    read (702, *) cdum
    ii = 0

    do
        read (702, *, end=600) idum1, idum2
        !read(702,'(2I5)',end=600) idum1, idum2
        ! Atom idum1 belongs to the NN list
        if (elist(idum1)) then
            ii = ii + 1
            exc(ii)%ref = idum1
            exc(ii)%n = idum2
            if (idum2 .ne. 0) then
                read (702, *, end=600) (exc_at(ii, jj), jj=1, exc(ii)%n)
                !read(702,'(20I7)',end=600) (exc_at(ii,jj), jj=1,exc(ii)%n)
            end if
        else
            if (idum2 .ne. 0) read (702, *, end=600) idum1
            !if (idum2.ne.0) read(702,'(20I7)',end=600) idum1
        end if
    end do
600 close (702)

    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) 'Exclusion list for first and second neighbours (GROMOS)'
        write (*, *) 'ATOM I   NEXCL'
        write (*, *) 'J1  J2  J3  J4 ...JNEXCL'
        do i = 1, nat_tot
            write (*, *) exc(i)%ref, exc(i)%n
            if (exc(i)%n .ne. 0) write (*, *) (exc_at(i, j), j=1, exc(i)%n)
        end do
        write (*, *) ''
    end if

    ! 1-4 neighbours
    allocate (exc14(nat_tot))
    allocate (exc_14(nat_tot, nmax14))
    exc_14 = 0.d0

    open (703, file='14.dat', form='formatted', status='old')
    read (703, *) cdum
    read (703, *) cdum
    read (703, *) cdum
    ii = 0

    do
        read (703, *, end=601) idum1, idum2
        !read(703,'(2I7)',end=601) idum1, idum2
        ! Atom idum1 belongs to the NN list
        if (elist(idum1)) then
            ii = ii + 1
            exc14(ii)%ref = idum1
            exc14(ii)%n = idum2
            if (idum2 .ne. 0) read (703, *, end=601) (exc_14(ii, jj), jj=1, exc14(ii)%n)
            !if (idum2.ne.0) read(703,'(40I7)',end=601) (exc_14(ii,jj), jj=1,exc14(ii)%n)
        else
            if (idum2 .ne. 0) read (703, *, end=601) idum1
            !if (idum2.ne.0) read(703,'(20I7)',end=601) idum1
        end if
    end do
601 close (703)

    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) 'Exclusion list for 1-4 neighbours (GROMOS)'
        write (*, *) 'ATOM I   NEXCL'
        write (*, *) 'J1  J2  J3  J4 ...JNEXCL'
        do i = 1, nat_tot
            write (*, *) exc14(i)%ref, exc14(i)%n
            if (exc14(i)%n .ne. 0) write (*, *) (exc_14(i, j), j=1, exc14(i)%n)
        end do
        write (*, *) ''
    end if

    return

end subroutine exclusion_list

!**********************************************************
logical function compute_vdw(i, j)
    !**********************************************************

    use exc_list
    use van_der_waals, only: grom

    implicit none

    integer, intent(in) :: i, j
    integer :: k

    compute_vdw = .true.
    do k = 1, exc(i)%n
        if (grom(j)%ind .eq. exc_at(i, k)) then
            compute_vdw = .false.
            exit
        end if
    end do

    return

end function compute_vdw

!**********************************************************
logical function compute_14(i, j)
    !**********************************************************

    use exc_list
    use van_der_waals, only: grom

    implicit none

    integer, intent(in) :: i, j
    integer :: k

    compute_14 = .false.
    do k = 1, exc14(i)%n
        if (grom(j)%ind .eq. exc_14(i, k)) then
            compute_14 = .true.
            exit
        end if
    end do

    return

end function compute_14

!**********************************************************
logical function elist(dum)
    !**********************************************************

    use van_der_waals, only: grom, nat_tot

    implicit none

    integer, intent(in) :: dum
    integer :: k

    elist = .false.
    do k = 1, nat_tot
        if (grom(k)%ind .eq. dum) then
            elist = .true.
            exit
        end if
    end do

    return

end function elist

!***********************************************************
subroutine restr_read()
    !***********************************************************

    ! Red classical restraints for restrained Wf and GEO optimizations
    ! if mm_restr=.true.

    use allio, only: rank, nion
    use cl_restr
    use link_angle, only: cl_bond, cl_angle, cl_dimp, cl_dihe, pi

    character*50 :: cdum
    integer :: i

    open (913, file='restr.dat')

    if (rank .eq. 0) then
        write (*, *) ''
        write (*, *) 'Classical (MM) restraints'
        write (*, *) ''
    end if

    !Multiplicative factor
    read (913, *) mm_fact

    !Default values
    nbonds = 0
    nth = 0
    ndihe = 0
    ndimp = 0

    !Bonds
    read (913, *) cdum
    read (913, *) nbonds
    allocate (cl_bond(nbonds))
    do i = 1, nbonds
        read (913, *) cl_bond(i)%i, cl_bond(i)%j, cl_bond(i)%kbond, cl_bond(i)%req
    end do

    if (rank .eq. 0) then
        write (*, *) nbonds, 'bond restraints'
        write (*, *) ''
    end if

    ! Angles
    read (913, *) cdum
    read (913, *) nth
    allocate (cl_angle(nth))
    do i = 1, nth
        read (913, *) cl_angle(i)%i, cl_angle(i)%j, cl_angle(i)%k, cl_angle(i)%ktheta, cl_angle(i)%thetaeq
        cl_angle(i)%thetaeq = pi/180.*cl_angle(i)%thetaeq
        cl_angle(i)%thetaeq = dcos(cl_angle(i)%thetaeq)
    end do

    if (rank .eq. 0) then
        write (*, *) nth, 'angle restraints'
        write (*, *) ''
    end if

    ! Dihedrals
    read (913, *) cdum
    read (913, *) ndihe
    allocate (cl_dihe(ndihe))
    do i = 1, ndihe
        read (913, *) cl_dihe(i)%i, cl_dihe(i)%j, cl_dihe(i)%k, cl_dihe(i)%l, cl_dihe(i)%kqhi, cl_dihe(i)%qcos
        cl_dihe(i)%qcos = pi/180.*cl_dihe(i)%qcos
        cl_dihe(i)%qcos = dcos(cl_dihe(i)%qcos)
    end do

    if (rank .eq. 0) then
        write (*, *) ndihe, 'dihedral restraints'
        write (*, *) ''
    end if

    ! Improper dihedrals
    read (913, *) cdum
    read (913, *) ndimp
    allocate (cl_dimp(ndimp))
    do i = 1, ndimp
        read (913, *) cl_dimp(i)%i, cl_dimp(i)%j, cl_dimp(i)%k, cl_dimp(i)%l, cl_dimp(i)%kqhi, cl_dimp(i)%qcos
        cl_dimp(i)%qcos = pi/180.*cl_dimp(i)%qcos
        cl_dimp(i)%qcos = dcos(cl_dimp(i)%qcos)
    end do

    if (rank .eq. 0) then
        write (*, *) ndimp, 'improper dihedral restraints'
        write (*, *) ''
    end if

    ! Allocate arrays for forces
    allocate (restr_f_bond(3, nion), restr_f_angle(3, nion), restr_f_dihe(3, nion), restr_f_dimp(3, nion))

    close (913)

    return

end subroutine restr_read

!**********************************************************************
subroutine deallocate_restr()

    use cl_restr
    use link_angle, only: cl_bond, cl_angle, cl_dimp, cl_dihe

    deallocate (cl_bond, restr_f_bond)
    deallocate (cl_angle, restr_f_angle)
    deallocate (cl_dimp, restr_f_dimp)
    deallocate (cl_dihe, restr_f_dihe)

    return

end subroutine deallocate_restr

!****************************************************************
subroutine r_bond()
    !****************************************************************

    ! Simple harmonic restraint around the distance req

    use allio, only: rion
    use cl_restr, only: mm_fact, nbonds, restr_bond
    use link_angle, only: cl_bond
    use tot_angle, only: ncount_bond, ave_bond, ave2_bond

    implicit none

    integer :: i, j, kk
    real*8 :: rr, x, y, z

    restr_bond = 0.d0
    do kk = 1, nbonds
        i = cl_bond(kk)%i
        j = cl_bond(kk)%j
        x = rion(1, i) - rion(1, j)
        y = rion(2, i) - rion(2, j)
        z = rion(3, i) - rion(3, j)
        rr = sqrt(x**2 + y**2 + z**2)
        restr_bond = restr_bond + cl_bond(kk)%kbond*(rr - cl_bond(kk)%req)**2
    end do
    restr_bond = 0.5d0*mm_fact*restr_bond

    ! Accumulating the bond potential
    ncount_bond = ncount_bond + 1
    ave_bond = ave_bond + restr_bond
    ave2_bond = ave2_bond + restr_bond*restr_bond

    return

end subroutine r_bond

!********************************************************************
subroutine r_angle()
    !********************************************************************

    ! Simple harmonic restraint aroun the bond angle

    use allio, only: rion
    use cl_restr, only: mm_fact, nth, restr_angle
    use link_angle, only: cl_angle
    use tot_angle, only: ncount_angle, ave_angle, ave2_angle

    implicit none

    integer :: i, j, k, ll
    real*8 :: bij, x, y, z, x1, y1, z1, bjk, costh

    restr_angle = 0.d0

    do ll = 1, nth
        i = cl_angle(ll)%i
        j = cl_angle(ll)%j
        k = cl_angle(ll)%k

        ! I-J DISTANCE
        x = rion(1, i) - rion(1, j)
        y = rion(2, i) - rion(2, j)
        z = rion(3, i) - rion(3, j)
        bij = sqrt(x**2 + y**2 + z**2)

        ! K-J distance
        x1 = rion(1, k) - rion(1, j)
        y1 = rion(2, k) - rion(2, j)
        z1 = rion(3, k) - rion(3, j)
        bjk = sqrt(x1**2 + y1**2 + z1**2)
        ! COS(THETA)
        costh = (x*x1 + y*y1 + z*z1)/(bij*bjk)
        restr_angle = restr_angle + cl_angle(ll)%ktheta*(costh - cl_angle(ll)%thetaeq)**2
    end do

    restr_angle = 0.5d0*mm_fact*restr_angle

    !Accumulating the bond angle potential
    ncount_angle = ncount_angle + 1
    ave_angle = ave_angle + restr_angle
    ave2_angle = ave2_angle + restr_angle*restr_angle

    return

end subroutine r_angle

!**********************************************************************************
subroutine r_dihe()
    !**********************************************************************************

    ! Simple harmonic restraint around a proper dihedral

    use allio, only: rion
    use cl_restr, only: mm_fact, ndihe, restr_dihe
    use link_angle, only: cl_dihe
    use tot_angle, only: ncount_dihed, ave_dihed, ave2_dihed
    use constants, only: pi

    implicit none

    integer :: i, j, k, l, mm
    real*8 :: phi, cosphi, rij(3), rkj(3), rkl(3), rmj2, rnk2, rnk(3), rmj(3)
    real*8 :: dot, dot1, rkj2, app, app1(3), cosine, rim(3), rln(3), dum

    restr_dihe = 0.d0

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
        rim(:) = rij(:) - (dot/rkj2)*rkj(:)
        rln(:) = -rkl(:) + (dot1/rkj2)*rkj(:)

        app = dot_product(rim, rln)
        app = app/sqrt(dot_product(rim, rim)*dot_product(rln, rln))
        phi = dacos(app)
        call cross(rkj, rkl, app1)
        phi = sign(1.d0, dot_product(rij, app1))*phi
        !phi = modulo(phi,2.d0*pi)
        cosine = dcos(phi)
        !call cos_mult(phi,cl_dihe(mm)%mult,cosine)
        restr_dihe = restr_dihe + cl_dihe(mm)%kqhi*(cosine - cl_dihe(mm)%qcos)**2

    end do

    restr_dihe = 0.5d0*mm_fact*restr_dihe

    !Accumulating the proper dihedral potential
    ncount_dihed = ncount_dihed + 1
    ave_dihed = ave_dihed + restr_dihe
    ave2_dihed = ave2_dihed + restr_dihe*restr_dihe

    return

end subroutine r_dihe

!*****************************************************************************
subroutine r_dimp()
    !*****************************************************************************

    ! Simple harmonic restraint for improper dihedrals

    use allio, only: rion
    use cl_restr, only: mm_fact, ndimp, restr_dimp
    use link_angle, only: cl_dimp
    use constants, only: pi
    use tot_angle, only: ncount_impr, ave_impr, ave2_impr

    implicit none

    integer :: i, j, k, l, mm
    real*8 :: cosphi, phi, rij(3), rkj(3), rkl(3), rmj2, rnk2, rnk(3), rmj(3)

    restr_dimp = 0.d0

    do mm = 1, ndimp

        i = cl_dimp(mm)%i
        j = cl_dimp(mm)%j
        k = cl_dimp(mm)%k
        l = cl_dimp(mm)%l

        call compute_dihed(rion(:, i), rion(:, j), rion(:, k), rion(:, l), phi, cosphi, rij, rkj, rkl, rmj2, rnk2, rnk, rmj)
        restr_dimp = restr_dimp + cl_dimp(mm)%kqhi*modulo(phi - cl_dimp(mm)%qcos, 2.d0*pi)**2
    end do

    restr_dimp = 0.5d0*mm_fact*restr_dimp

    !Accumulating the improper dihedral potential
    ncount_impr = ncount_impr + 1
    ave_impr = ave_impr + restr_dimp
    ave2_impr = ave2_impr + restr_dimp*restr_dimp

    return

end subroutine r_dimp

!***********************************************************
subroutine write_traj(nel, nelup, kel, eloc)
    !***********************************************************

    implicit none

    integer, intent(IN) :: nel, nelup
    double precision, intent(IN) :: kel(3, nel), eloc

    integer :: i
    double precision :: autoang

    autoang = 0.529177

    ! only the first walker (the subroutine is called only if rank.eq.0)
    ! Random walk in Angstroms
    write (671, '(I4)') nel
    write (671, *) 'Random walk:', nelup, 'up', nel - nelup, 'down, Eloc = ', eloc*0.5d0
    ! nelup electrons
    do i = 1, nelup
        write (671, '(1A,1X,3F16.8)') 'N', kel(1, i)*autoang, kel(2, i)*autoang, kel(3, i)*autoang
    end do
    ! neldo electrons
    do i = nelup + 1, nel
        write (671, '(1A,1X,3F16.8)') 'O', kel(1, i)*autoang, kel(2, i)*autoang, kel(3, i)*autoang
    end do

    return

end subroutine write_traj

