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

subroutine read_corr_fun(nel, nelup, nion, iespbc, celldm, rs, cellscale &
                         , ngen, ell, ncell, nbias, maxf, ibinit, lbin, iskip, ddim, cutk, vdim, ngrid_l, r_offset, ngrid_p &
                         , ifrho, ifspin, ifkspin, kspin, ifpair, iffluct, ifsofk &
                         , fermi_flag, err_stop, ifcorrs, shiftlog, ioptread, longio, noeloc &
                         , sphere_radius, allshells, outofplane, rdf_for_atom)
    ! read fort.10 and readforward.input (standard input) with namelist

    use constants
    use Assar_module, only: ifrho_assar, kswitch, nswitch, assar_parr, assar_cut
    use grid_module, only: ext_grid, grid_points, da, center, grid_start
    use qpwf_module, only: ifqpwf, ifqpwf_h, ifqpwf_k, ifqpwf_extr, &
                           n_extr_points, extr_point, decouple_k, decouple_files
    use Rho_corr_module, only: ifrho_corr
    use Spin2, only: ifspin2
    use Dipole_module, only: ifdipole
    use allio, only: rion
    use cell, only: givens2r, s2r, yes_tilted
    use berry_phase, only: ifberry
    implicit none
    integer nion, corr_factors, fwd_propagations, j &
        , bin_length, initial_bin, fwd_skip, corrfun_dim, i, ngen, nbias, maxf, lbin, iskip, ddim, ibinit, iopt, ioptread
    integer ncell(3), ngrid(3), cart_axes(3), vdim(3), ngrid_l(3), nel, nelup, ngrid_p, radial_grid
    integer :: info, ipiv(3)
    real*8 ell(3), k_cutoff, cutk, offset(3), r_offset(3) &
        , cellscale(3), omega, celldm(*), rs, shiftlog, kspin(3), sphere_radius
    logical err_stop, fermi_flag, assar_density, new_density, charge_density, spin_density, pair_corr_fun, structure_factor &
        , allshells, ifrho, ifspin, ifpair, ifkspin, ifsofk, correlated_samp, ifcorrs, iespbc, fluctuations, iffluct
    logical longio
    character(2) cartesian(3)
    ! (Matteo change) Added logical variables assar_density and input_grid
    logical noeloc
    ! (Matteo) Logical for computing STM
    logical quasi_particle, quasi_hole, quasi_hole_k, quasi_hole_extr
    integer center_type
    double precision :: starting_point(3)
    double precision :: matscra(3, 3), vecscra(3)
    real(8) outofplane
    ! (Ye) an logical variable for computing spin^2.
    logical compute_spin2
    logical dipole_moment
    ! (Kosuke Nakano) an logical variable for computing Radial distribution function of atom
    logical rdf_for_atom
    !namelist /simulation/ ngen,nw,number_files,wherescratch
    namelist /simulation/ ngen, iopt, longio, decouple_k, decouple_files
    namelist /system/ ncell, ell, center_type, starting_point, da, extr_point, n_extr_points
    namelist /corrfun/ corr_factors, fwd_propagations, bin_length, initial_bin, fwd_skip, corrfun_dim &
        , k_cutoff, assar_density, grid_points, assar_parr, charge_density, new_density, spin_density, pair_corr_fun, &
        structure_factor, ngrid, cart_axes, offset, correlated_samp, shiftlog, ifkspin, kspin, assar_cut, noeloc, &
        kswitch, nswitch, compute_spin2, radial_grid, sphere_radius, quasi_particle, quasi_hole, quasi_hole_k, &
        quasi_hole_extr, dipole_moment, allshells, fluctuations, ifberry, outofplane, rdf_for_atom

    err_stop = .false.
    cartesian(1) = ' x'
    cartesian(2) = ' y'
    cartesian(3) = ' z'

    open (55, file='readforward.input', status='old', form='formatted')

    ! read commands of readforward
    decouple_k = .true.
    decouple_files = .false.
    ngen = 0
    iopt = 1 ! default value of iopt
    read (55, nml=simulation)
    ioptread = iopt

    ell = 0.d0
    ncell = 1

    n_extr_points = 0
    extr_point(:) = 0
    center_type = 1
    starting_point(:) = 0
    da(:) = 0

    read (55, nml=system)

    if (iespbc) then
        if (.not. givens2r) then
            omega = celldm(2)*celldm(3)
            cellscale(1) = (PI*nel*4.d0/3.d0/omega)**(1.d0/3.d0)*rs
            cellscale(2) = cellscale(1)*celldm(2)
            cellscale(3) = cellscale(1)*celldm(3)
        else
            matscra = s2r
            call dgetrf(3, 3, matscra, 3, ipiv, info)
            omega = matscra(1, 1)
            do i = 2, 3
                omega = omega*matscra(i, i)
            end do
            omega = abs(omega)
            do i = 1, 3
                cellscale(i) = sqrt(s2r(1, i)**2 + s2r(2, i)**2 + s2r(3, i)**2)
            end do
        end if
    else
        omega = 1.d0
        celldm(2) = 1.d0
        celldm(3) = 1.d0
        rs = 1.d0
        if (ell(1)*ell(2)*ell(3) .eq. 0.d0) then
            write (6, *) ' Default box 10x10x10 a.u.. Change ell in input to modify it! '
            cellscale(1) = 10.
            cellscale(2) = 10.
            cellscale(3) = 10.
            ell(:) = cellscale(:)
        else
            cellscale(:) = ell(:)
            !        write(6,*) 'cellscale inside',cellscale(:)
        end if
    end if

    ! change cellscale
    cellscale(:) = cellscale(:)/ncell(:)

    if (ell(1) .eq. 0.d0 .and. ell(2) .eq. 0.d0 .and. ell(3) .eq. 0.d0) then
        ell(:) = cellscale(:)
    end if

    write (6, *) 'SYSTEM PARAMETERS'

    if (nel - nelup .eq. 0) then
        fermi_flag = .false.
        write (6, *) 'Bose particles, or spinless fermions'
    else
        fermi_flag = .true.
        write (6, *) 'Fermi particles with spin'
    end if

    new_density = .false. ! (Matteo) Added default for flag assaraf density
    assar_density = .false. ! (Matteo) Added default for flag assaraf density
    assar_parr = 1.d0
    assar_cut = 1000.d0 ! (Matteo) default per il cut off
    longio = .false.
    noeloc = .false.
    charge_density = .false.
    rdf_for_atom = .false. ! (Kosuke Nakano) computing rdf of an atom
    spin_density = .false.
    pair_corr_fun = .false.
    fluctuations = .false.
    structure_factor = .false.
    correlated_samp = .false.
    shiftlog = 0.d0
    compute_spin2 = .false. ! (Ye) default computing spin^2 switched off
    dipole_moment = .false. ! (Shibing) no dipole moment by default

    quasi_particle = .false.
    quasi_hole = .false.
    quasi_hole_k = .false.
    quasi_hole_extr = .false.

    grid_points = 0
    ext_grid = .false.

    corr_factors = 0
    fwd_propagations = 0
    bin_length = ngen/100
    initial_bin = 2
    fwd_skip = 1
    corrfun_dim = 3
    k_cutoff = 3.d0
    shiftlog = 0.d0

    cart_axes(1) = 1
    cart_axes(2) = 2
    cart_axes(3) = 3

    offset(1) = 0.d0
    offset(2) = 0.d0
    offset(3) = 0.d0

    ngrid(1) = 0
    ngrid(2) = 0
    ngrid(3) = 0
    kspin = 0.d0
    ifkspin = .false.
    ! by E. Coccia (19/7/11): default value for Chernomor formulas
    kswitch = 0.1d0
    nswitch = 8

    radial_grid = 0
    sphere_radius = 0.d0
    allshells = .false.
    ifberry = .false.

    outofplane = 0.d0

    read (55, nml=corrfun)

    write (6, *) 'CORRELATION FUNCTION PARAMETERS'
    write (6, *) 'corrective factors (corr_factors) =', corr_factors
    write (6, *) 'forward propagations (fwd_propagations) =', fwd_propagations
    write (6, *) 'skips forward (fwd_skip) =', fwd_skip
    write (6, *) 'bin length (bin_length) =', bin_length
    write (6, *) 'initial bin (initial_bin) =', initial_bin
    write (6, *) 'active dimension (corrfun_dim) =', corrfun_dim
    write (6, *) 'cartesian coordinates (cart_axes) =', (cartesian(cart_axes(i)), i=1, corrfun_dim)
    write (6, *) 'offset =', (offset(i), i=1, corrfun_dim)
    write (6, *) 'cutoff log wave function (shiftlog) =', shiftlog

    if (sum(abs(offset(:))) .eq. 0.d0) then
        if (ngrid(1) .ge. 1 .and. ngrid(2) .ge. 1 .and. ngrid(3) .ge. 1) then
            write (6, *) 'Default shift offset for centering xcrysden -0.5,-0.5,-0.5 in mesh grid units'
            if (.not. yes_tilted) then
                offset(:) = offset(:) - 0.5d0*cellscale(:)/ngrid(:)
            else
                vecscra(:) = 0.d0
                do j = 1, 3
                    vecscra(:) = vecscra(:) - s2r(:, j)*0.5d0/ngrid(j)
                end do
                offset(:) = offset(:) + vecscra(:)
            end if
        end if
    end if
    !  offset(:)=offset(:)*cellscale(:)/ngrid(:)

    nbias = corr_factors
    maxf = fwd_propagations
    lbin = bin_length
    ibinit = initial_bin
    iskip = fwd_skip
    ddim = corrfun_dim
    cutk = k_cutoff

    if (iskip .le. 0) then
        write (6, *) 'input fatal error'
        write (6, *) 'fwd_skip must be positive!'
        err_stop = .true.
    end if

    ! added by Kosuke Nakano on 29 May 2019
    if (rdf_for_atom .and. ddim .ne. 1) then
        write (6, *) 'input fatal error'
        write (6, *) 'rdf_for_atom is valid only for corrfun_dim=1'
        err_stop = .true.
    end if

    ifrho_corr = new_density ! (Matteo) Added default for flag assaraf density
    ifrho_assar = assar_density !(Matteo)

    ifrho = charge_density
    ifspin = spin_density

    iffluct = fluctuations

    if ((.not. ifrho) .and. (.not. ifspin) .and. fluctuations) then
        write (6, *) 'warning: to compute local charge or spin fluctuations, you must switch on the charge or spin flag!'
        iffluct = .false.
    end if

    ifqpwf = quasi_particle
    ifqpwf_h = quasi_hole
    ifqpwf_k = quasi_hole_k
    ifqpwf_extr = quasi_hole_extr

    if (ifqpwf_k .or. ifqpwf_h .or. ifqpwf_extr) ifqpwf = .true.

    if (ifqpwf_extr) then
        ifqpwf_h = .false.
        ifqpwf_k = .false.
    end if

    ifpair = pair_corr_fun
    ifsofk = structure_factor
    ifcorrs = correlated_samp
    ifspin2 = compute_spin2 !(Ye)

    ifdipole = dipole_moment !(Shibing)

    center = center_type
    grid_start(:) = starting_point(:)

    vdim = cart_axes
    if (iespbc) then
        r_offset = offset
    else
        ! added by Kosuke Nakano 30 May 2019
!        if(.not.rdf_for_atom.and..not.ifdipole) then
!            if(corrfun_dim.eq.3) then
!               r_offset(:) = offset(:) - cellscale(:) / 2.d0
!            endif
!        else
!            r_offset = offset
!        endif
        if (.not. rdf_for_atom .and. .not. ifdipole .and. corrfun_dim .eq. 3) then
            r_offset(:) = offset(:) - cellscale(:)/2.d0
        else
            r_offset = offset
        end if

        !     write(6,*) 'offset',offset(:)
        !     write(6,*) 'cellscale',cellscale(:)
        !     write(6,*) 'r_offset',r_offset(:)
    end if
    ngrid_l = ngrid
    ngrid_p = radial_grid

    if (ifrho_assar) then
        write (6, *) 'warning: calculating only the charge distribution with the &
                &method proposed \n by Assaraf the rest is switched off'
        charge_density = .false.
        spin_density = .false.
        pair_corr_fun = .false.
        quasi_particle = .false.
        quasi_hole = .false.
        quasi_hole_k = .false.
        structure_factor = .false.
        correlated_samp = .false.
        compute_spin2 = .false.
        ifrho = charge_density
        ifspin = spin_density
        ifpair = pair_corr_fun
        ifsofk = structure_factor
        ifcorrs = correlated_samp
        ifspin2 = compute_spin2 !(Ye)
        da(:) = 0
        if (grid_points .gt. 0) then
            ext_grid = .true.
        else
            grid_points = ngrid(1)*ngrid(2)*ngrid(3)
        end if
        write (*, *) 'Switch value for the grid =', kswitch
        write (*, *) 'Steepness n value for Ks and Ka (see O. Chernomor''s thesis) =', nswitch
    end if

    ! (Matteo) Added for STM to fix the grid type.
    if ((ifqpwf .and. .not. ifqpwf_extr) .or. ifrho_corr) then
        if (grid_points .gt. 0) then
            ext_grid = .true.
        else
            grid_points = ngrid(1)*ngrid(2)*ngrid(3)
        end if
    end if

    if (ifpair .or. ifspin .or. ifrho .or. (grid_points .eq. 0 .and. ifrho_assar)) then
        if (sphere_radius .eq. 0.d0) then
            write (6, *) 'ngrid (ngrid) =', (ngrid(i), i=1, corrfun_dim)
        else
            write (6, *) 'Warning: site occupation instead of standard density!!!'
            if (ddim .eq. 2) then
                write (6, *) 'site defined by cylinder radius of', sphere_radius, 'a_0'
                if (outofplane .ne. 0.d0) write (6, *) 'take only out-of-plane configurations >', outofplane, 'a_0'
            else
                write (6, *) 'site defined by sphere radius of', sphere_radius, 'a_0'
            end if
            if (allshells) write (6, *) ' Warning: All bonds considered without averaging according to the distance '
            ngrid_l = 0
            ngrid_p = 0
            if (ifkspin) then
                write (6, *) 'warning: site occupation not compatible with SDW yet !!!'
                write (6, *) 'ifkspin set to .false.'
                ifkspin = .false.
            end if
            if (iffluct) then
                write (6, *) 'warning: in the site approach, fluctuations are already computed as on-site pair correlations!'
                write (6, *) 'fluctuations option has been disabled'
                iffluct = .false.
            end if
        end if
    end if

    if (ifsofk) then
        write (6, *) 'cutoff in k space for S(k)', cutk
    end if

    if ((ifpair .and. .not. ifspin .and. .not. ifrho)) then
        ifrho = .true.
        ifspin = .true.
        write (6, *) 'warning: charge_density and spin_density switched on'
    end if

    if (ifpair .and. .not. iespbc) then
        write (6, *) 'warning: to compute pair_corr_fun we assume translational invariance!'
    end if

    if (ifsofk .and. .not. iespbc) then
        write (6, *) 'warning: to compute the structure factor we assume translational invariance!'
    end if

    if ((ifpair .or. ifspin .or. ifrho .or. (grid_points .eq. 0 .and. ifrho_assar)) .and. sphere_radius .eq. 0.d0) then
        if (ngrid(1) .eq. 0 .and. ngrid(2) .eq. 0 .and. ngrid(3) .eq. 0) then
            write (6, *) 'input fatal error'
            write (6, *) 'you must specify ngrid in each direction to compute density and pair_corr_fun'
            err_stop = .true.
        end if
    end if

    if (ifpair .or. ifspin .or. ifrho .or. ifsofk .or. (grid_points .eq. 0 .and. ifrho_assar)) then
        write (6, *) 'dimension of the cell for correlation functions'
        write (6, *) 'ell(x) =', ell(1)
        write (6, *) 'ell(y) =', ell(2)
        write (6, *) 'ell(z) =', ell(3)
    end if

    !

    close (55)

end subroutine read_corr_fun
