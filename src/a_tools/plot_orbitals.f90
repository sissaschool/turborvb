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

program plot_orbitals
    use allio
    ! use constants, only: pi,ipc,ipf

    implicit none
    integer, parameter :: ufort10 = 10
    integer :: iorb, nfil, ntot, indorb, imu
    integer :: i1, i, kk, ii, iii, mesh(3), ind, j, imol, lq, lq2, ind_mol, chosen_mol, k, up_down
    real(8) :: rpoint(3), step(3), origin(3), center(3), cell_loc(3)
    complex(8) :: totmagc
    real(8) :: totmag, magorb
    real(8), allocatable :: datagrid(:, :, :), datagrid_sum(:, :, :), datagrid_prod(:, :, :)
    real(8), allocatable :: buffer_winv(:), norm(:), writefile(:), psi(:)
    integer, allocatable :: imap_loc(:)
    real(8) :: r1, r2, kspin(3)
    character(15) :: question, question2
    logical, allocatable :: control(:)
    logical :: spinon, chargeon, nochsp
    integer :: lower, upper, choice, iesdvj, jj, iesdr1iesd
    real(8), allocatable :: psir(:, :, :, :)
    real(8) :: vol, costocc, sumval, r0, rc(3), psilg
    complex(8) :: psi_c

    real(8), parameter :: minimum_distance = 1.0d-9 ! minimum el/ion distance accepted

    ! scalar products
    real(8), external :: ddot, jastrow_ei
    complex(8), external :: zdotu

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(60) :: name_tool
    character(20) :: str
    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'plot_orbitals'
        call help_online(name_tool)
        stop
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open (unit=ufort10, file='fort.10', form='formatted', status='unknown')
    call default_allocate
    yesfast = 0 ! it works only with allocation of detmat
    rank = 0
    nw = 1
    nproc = 1
    in1 = 1

    call read_fort10_fast
    npsa = npsar

    pseudofile = "pseudo.dat"
    iflagerrall = 0

    if (npsa .gt. 0) nintpsa = 6
    call read_pseudo
    call read_fort10(ufort10)

    iesdvj = iesdr1iesd(iesdr)

    if (.not. molyes) call error('plot_orbitals', ' This tool work only with molecular orbitals', 1, 0)
    write (*, *) ' Number of molecular orbitals : ', molecular

    do i = 1, 3
        center(i) = sum(rion(i, 1:nion))/nion
    end do

    if (.not. iespbc) then
        write (*, *) ' Choose box size (x,y,z) '
        read (*, *) cell_loc(:)
        write (*, *) cell_loc(:)
        cellscale(1:3) = cell_loc(1:3)
        do i = 1, nion
            rion(1:3, i) = rion(1:3, i) - center(1:3) + cell_loc(1:3)/2.d0
        end do
        ! for open system the center of the box is the (0,0,0) point
        origin = 0.d0
    else
        !call ApplyPBC(rion,nion)
        !do i=1,nion
        !rion(1:3,i)=rion(1:3,i)+cellscale(1:3)/2.d0
        !enddo
        cell_loc(1:3) = cellscale(1:3)
        write (*, *) ' Choose shift reference unit cell PBC : '
        read (*, *) origin(:)
        origin(1:3) = origin(1:3)*cellscale(1:3)
    end if

    write (*, *) ' Choose number of mesh points (x,y,z) : '
    read (*, *) mesh(:)
    write (*, *) mesh(:)

    ! To match with output plot in xcrysden
    if (iespbc) mesh(:) = mesh(:) + 1

    write (*, *) ' Choose orbitals to tabulate (possible answers all, partial, charge, spin) : '
    read (*, *) question
    lq = index(question, ' ') - 1
    write (*, *) question(1:lq)

    spinon = .false.
    chargeon = .false.
    nochsp = .true.
    lower = 1
    upper = molecular

    if (question(1:lq) .eq. 'partial') then
        write (*, *) ' Please give a range between 1 and ', molecular
        read (*, *) lower, upper
        write (*, *) lower, upper
    end if

    if (question(1:lq) .eq. 'charge') then
        write (*, *) ' Please give  the lowest molecular orbital within 1 and ', molecular
        read (*, *) lower
        upper = -1
        write (6, *) 'Number of fully occupied molecular orbital/total number occupied by up and down ?'
        read (5, *) nfil, ntot

    end if

    if (question(1:lq) .eq. 'spin') then
        write (*, *) ' Please give  the lowest molecular orbital within 1 and ', molecular
        read (*, *) lower
        upper = -2
        write (6, *) 'Number of fully occupied molecular orbital/total number occupied by up and down ?'
        read (5, *) nfil, ntot

        write (6, *) ' Momentum magnetization ? (unit 2pi/cellscale) '
        read (5, *) kspin(1:3)
        kspin(1:3) = 2.d0*pi*kspin(1:3)/cellscale(1:3)

        write (6, *) ' K rescaled =', kspin(1:3)

    end if

    if (upper .eq. -1) then
        chargeon = .true.
        nochsp = .false.
        if (symmagp) then
            if (nelup - neldo .ne. 0) then
                upper = molecular
            else
                upper = ntot
            end if
        else
            if (nelup - neldo .ne. 0) then
                upper = molecular
            else
                upper = 2*ntot
            end if
        end if
        write (6, *) ' Warning computing charge density '
    elseif (upper .eq. -2) then
        spinon = .true.
        nochsp = .false.
        if (nelup - neldo .ne. 0) then
            upper = molecular
        else
            upper = 2*ntot
        end if
        if (symmagp) call error("plot_orbitals", " Spin density is not possible for symmagp=.true.", 1, 0)
        write (6, *) ' Warning computing spin density '
    end if

    allocate (imap_loc(molecular), control(molecular), psi(ipc))

    control = .false.
    indorb = 0
    imol = 0
    ind = 0
    ind_mol = 0

    do i = 1, nshell_c
        if (ioptorb_c(i) .eq. 1000000) then
            ind = ind + 1
            if (ioccup_c(ind) .eq. 1) then
                indorb = indorb + 1
                imol = imol + 1
                imap_loc(imol) = indorb
                if (imol .ge. lower .and. imol .le. upper) then
                    control(imol) = .true.
                    ind_mol = ind_mol + 1
                end if
            end if
        else
            do j = 1, mult_c(i)
                ind = ind + 1
                if (ioccup_c(ind) .eq. 1) indorb = indorb + ipf
            end do
        end if
    end do

    chosen_mol = ind_mol
    write (*, *) ' # of orbitals written =', chosen_mol

    allocate (norm(molecular + 1))
    norm = 0.d0

    if (iespbc) then
        step(:) = cell_loc(:)/dble(mesh(:) - 1)
    else
        step(:) = cell_loc(:)/dble(mesh(:))
    end if
    vol = step(1)*step(2)*step(3)

    if (chargeon .or. spinon) then
        allocate (datagrid(mesh(1), mesh(2), mesh(3)))
    else
        allocate (datagrid(ipc*mesh(1), mesh(2), mesh(3)))
        if (ipf .eq. 2) allocate (datagrid_sum(ipc*mesh(1), mesh(2), mesh(3)), datagrid_prod(ipc*mesh(1), mesh(2), mesh(3)))
    end if

    allocate (buffer_winv(ipc*nelorb))
    iflagnorm = 3

    ! compute psi^2 normalization for each orbital
    if (nochsp) then

        do j = 1, molecular
            norm(j) = 1.d0
        end do

    else
        if (symmagp) then
            norm = 0.d0
            norm(1:nfil) = 2.d0

            if (nfil .lt. ntot) then
                costocc = dble(nel - 2*nfil - nelup + neldo)/(ntot - nfil)
                norm(nfil + 1:ntot) = costocc
            end if

            if (nelup .gt. neldo) norm(molecular - nelup + neldo + 1:molecular) = 1.d0

        else
            norm(1:2*nfil) = 1.d0

            if (nfil .lt. ntot) then
                costocc = dble(nel - 2*nfil - nelup + neldo)/(2*ntot - 2*nfil)
                norm(2*nfil + 1:2*ntot) = costocc
            elseif (nelup .gt. neldo) then
                norm(molecular - nelup + neldo + 1:molecular) = 1.d0
            end if

        end if
    end if

    datagrid = 0.d0
    magorb = 0.d0

    ind_mol = 0

    psip = 0.d0

    do j = 1, molecular

        !!It does not do anything if ipf=1, otherwise it allows you to print
        !!the up and the down part of the molecular orbitals
        do up_down = 1, ipf
            imu = (up_down - 1)*nelorbh + 1

            if (control(j)) then
                !
                !
                if (up_down .eq. 1) ind_mol = ind_mol + 1
                do iii = 1, mesh(3)
                    do ii = 1, mesh(2)
                        do i = 1, mesh(1)
                            !
                            !                   rpoint(1)=(i-1)*step(1)-origin(1)
                            !                   rpoint(2)=(ii-1)*step(2)-origin(2)
                            !                   rpoint(3)=(iii-1)*step(3)-origin(3)

                            rpoint(:) = (i - 1)*step(1)*at(:, 1) + (ii - 1)*step(2)*at(:, 2) + (iii - 1)*step(3)*at(:, 3)

                            if (ipc .eq. 2 .and. mod(j, 2) .eq. 1) then
                                call upnewwf(1, 0, 0, 1, nshell, ioptorb, ioccup, rpoint, 1, r, rmu, &
                                             dupr, zetar, rion, psip, buffer_winv, nelorb, nion, kion, &
                                             iflagnorm, cnorm, LBox, rmucos, rmusin, minimum_distance, &
                                             indpar_tab, indorb_tab, indshell_tab, .false.)
                            else
                                call upnewwf(1, 0, 0, 1, nshell, ioptorb, ioccup, rpoint, 1, r, rmu, &
                                             dupr, zetar, rion, psip, buffer_winv, nelorb, nion, kion, &
                                             iflagnorm, cnorm, LBox, rmucos, rmusin, minimum_distance, &
                                             indpar_tab, indorb_tab, indshell_tab, .true.)
                            end if

                            !                  if(add_onebody2det) then
                            psilg = -scale_one_body
                            do jj = 1, nion
                                if (iespbc) then
                                    rc(:) = rpoint(:) - rion(:, jj)
                                    call CartesianToCrystal(rc, 1)
                                    do kk = 1, 3
                                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                                    end do
                                    r0 = norm_metric(rc, metric)
                                else
                                    rc(:) = (rpoint(:) - rion(:, jj))*costz(jj)
                                    r0 = dsqrt(sum(rc(:)**2))
                                end if
                                psilg = psilg - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                            end do

                            buffer_winv(1:ipc*nelorb) = buffer_winv(1:ipc*nelorb)*dexp(psilg)

                            !                  endif

                            if (ipc .eq. 1) then
                                psi(1) = ddot(nelorbh, buffer_winv, 1, mu_c(imu, imap_loc(j)), 1)
                            else
                                psi_c = zdotu(nelorbh, buffer_winv, 1, mu_c(2*imu - 1, imap_loc(j)), 1)
                                psi(1) = real(psi_c)
                                psi(2) = aimag(psi_c)
                            end if

                            if (nochsp) then
                                datagrid(ipc*(i - 1) + 1:ipc*i, ii, iii) = psi(:)
                            elseif (chargeon) then
                                datagrid(i, ii, iii) = datagrid(i, ii, iii) + norm(j)*sum(psi(:)**2)
                            elseif (spinon) then
                                if (mod(j, 2) .eq. 0 .or. j .gt. 2*nfil) then
                                    datagrid(i, ii, iii) = datagrid(i, ii, iii) + norm(j)*sum(psi(:)**2)
                                else
                                    datagrid(i, ii, iii) = datagrid(i, ii, iii) - norm(j)*sum(psi(:)**2)
                                end if
                            end if
                        end do
                    end do
                end do

                if (spinon) then
                    if (iespbc) then
                        totmag = sum(abs(datagrid(1:mesh(1) - 1, 1:mesh(2) - 1, 1:mesh(3) - 1)))*vol
                    else
                        totmag = sum(abs(datagrid(:, :, :)))*vol
                    end if
                    write (6, *) ' Mag up to orbital =', ind_mol, imap_loc(j), totmag
                else
                    if (ipf .eq. 1) then
                        write (6, *) ' Molecular =', j, ind_mol
                    else
                        if (up_down .eq. 1) then
                            write (6, *) ' Up pfaffian Molecular =', j, ind_mol
                        else
                            write (6, *) ' Down pfaffian Molecular =', j, ind_mol
                        end if
                    end if
                end if
                !
                if (nochsp) then
                    if (ipf .eq. 1) then
                        datagrid = datagrid*norm(ind_mol)
                        call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                          , atom_number, iespbc, mesh, origin, datagrid, j, 'orbital       ')
                        datagrid = datagrid**2
                        call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                          , atom_number, iespbc, mesh, origin, datagrid, j, 'orbsqrd       ')
                    else
                        if (up_down .eq. 1) then
                            datagrid = datagrid*norm(ind_mol)
                            call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid, j, 'orbital_up       ')
                            datagrid_prod = datagrid
                            datagrid = datagrid**2
                            datagrid_sum = datagrid
                            call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid, j, 'orbsqrd_up       ')
                        else
                            datagrid = datagrid*norm(ind_mol)
                            datagrid_prod = datagrid_prod*datagrid
                            call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid, j, 'orbital_down       ')
                            call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid_prod, j, 'prod_updown       ')
                            datagrid = datagrid**2
                            datagrid_sum = datagrid_sum + datagrid
                            call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid, j, 'orbsqrd_down       ')
                            call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid_sum, j, 'orbsqrd_sum       ')
                            datagrid_sum = datagrid_sum - 2.0*datagrid
                            call plot_3d_data(ipc, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid_sum, j, 'orbsqrd_diff       ')
                        end if
                    end if

                elseif (j .eq. upper) then
                    if (spinon) then
                        if (iespbc) then
                            totmag = sum(abs(datagrid(1:mesh(1) - 1, 1:mesh(2) - 1, 1:mesh(3) - 1)))*vol
                        else
                            totmag = sum(abs(datagrid(:, :, :)))*vol
                        end if
                        totmagc = (0.d0, 0.d0)
                        do iii = 1, mesh(3) - 1
                            do ii = 1, mesh(2) - 1
                                do i = 1, mesh(1) - 1
                                    rpoint(1) = (i - 1)*step(1) - origin(1)
                                    rpoint(2) = (ii - 1)*step(2) - origin(2)
                                    rpoint(3) = (iii - 1)*step(3) - origin(3)
                                    sumval = sum(kspin(1:3)*rpoint(1:3))
                                    totmagc = totmagc + exp(dcmplx(0.d0, sumval))*datagrid(i, ii, iii)*vol
                                end do
                            end do
                        end do

                        write (6, *) ' Total magnetization in the cell (Bohr)', totmag
                        write (6, *) ' Total square root structure factor at K (Bohr)', abs(totmagc)

                    elseif (chargeon) then
                        if (iespbc) then
                            totmag = sum(datagrid(1:mesh(1) - 1, 1:mesh(2) - 1, 1:mesh(3) - 1))*vol
                        else
                            totmag = sum(datagrid(:, :, :))*vol
                        end if
                        write (6, *) ' Total charge in the cell ', totmag
                    end if
                    if (chargeon) then
                        if (ipf .eq. 1) then
                            call plot_3d_data(1, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid, 0, 'charge        ')
                        else
                            if (up_down .eq. 1) then
                                call plot_3d_data(1, cell_loc, cellscale, nion, rion &
                                                  , atom_number, iespbc, mesh, origin, datagrid, 0, 'charge_up        ')
                            else
                                call plot_3d_data(1, cell_loc, cellscale, nion, rion &
                                                  , atom_number, iespbc, mesh, origin, datagrid, 0, 'charge_down        ')
                            end if
                        end if
                    elseif (spinon) then
                        !     datagrid=datagrid/2.d0  Bohr magn. units.
                        if (ipf .eq. 1) then
                            call plot_3d_data(1, cell_loc, cellscale, nion, rion &
                                              , atom_number, iespbc, mesh, origin, datagrid, 0, 'spin          ')
                        else
                            if (up_down .eq. 1) then
                                call plot_3d_data(1, cell_loc, cellscale, nion, rion &
                                                  , atom_number, iespbc, mesh, origin, datagrid, 0, 'spin_up          ')
                            else
                                call plot_3d_data(1, cell_loc, cellscale, nion, rion &
                                                  , atom_number, iespbc, mesh, origin, datagrid, 0, 'spin_down          ')
                            end if
                        end if
                    end if !
                end if ! endif nocshp
            end if ! endif control
        end do
        !
        !
    end do

    call deallocate_all
    close (ufort10)

    stop

end program plot_orbitals
