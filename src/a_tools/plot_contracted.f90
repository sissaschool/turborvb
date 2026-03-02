!TL off
program plot_contracted

    ! this tool plots contracted orbitals if present in the wave function
    ! one needs to select the atom the contracted are centered on,
    ! and the number of total contracted plotted
    ! the other input parameters are similar to the ones of plot_orbitals
    !
    ! when hubbardU_logic is set to .true. the code computes the spread of the
    ! contracted orbitals and the Coulomb integrals involving them, instead of
    ! computing the orbitals on grid for xcrysden plotting

    use allio
    use constants, only : pi, ipc
    use logger_io, only: log_info, log_warning

    implicit none

    integer, parameter :: ufort10 = 10
    real(8), parameter :: minimum_distance = 1.0d-9  ! minimum el/ion distance accepted

    integer :: iorb, nfil, ntot, indorb, chosen_ion, imu
    integer :: i1, i, ii, kk, iii, mesh(3), ind, j, imol, ind_mol, chosen_mol, k
    integer, allocatable :: nhybrids(:)
    double precision :: rpoint(3), step(3), origin(3), center(3), cell_loc(3)
    double precision, allocatable :: datagrid(:, :, :)
    double precision, allocatable :: buffer_winv(:), norm(:), psi(:)
    integer, allocatable :: imap_loc(:)
    double precision :: r1, r2, kspin(3)
    integer :: nradial, iteri, iterj, ind_ri
    integer :: nshift
    double precision :: U, U_ex, integ1, integ2, psi2
    logical, allocatable :: control(:)
    logical :: spinon, chargeon, nochsp, hubbardU_logic, analytic
    integer :: lower, upper, choice, jj
    real(8), allocatable :: psiri(:, :), ri(:, :), psirj(:, :), rj(:, :)
    real(8), allocatable :: HubbardU(:, :), HubbardJ(:, :)
    real(8), allocatable :: orb_center(:, :)
    real(8) :: vol, r0, rc(3), psilg
    real(8) :: ravx, ravy, ravz, rav2, rvec(3), r2vec, spread, psi2_norm
    real(8) :: analytic_integral
    complex(8) psi_c

    real(8), external :: ddot
    complex(8), external :: zdotu
    real(8), external :: jastrow_ei
    integer, external :: ioptorbcontr

    character(100) :: name_tool
    CHARACTER(20) :: str

    CALL getarg(1, str)
    IF (str.eq."--help".or.str.eq."-help".or.str.eq."help") THEN

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'plot_contracted'
        call help_online(name_tool)

        stop
    ENDIF

    open(unit = ufort10, file = 'fort.10', form = 'formatted', status = 'unknown')
    call default_allocate
    yesfast = 0             ! it works only with allocation of detmat
    rank = 0
    nw = 1
    nproc = 1
    in1 = 1

    call read_fort10_fast
    npsa = npsar

    pseudofile = "pseudo.dat"
    iflagerrall = 0

    if(npsa.gt.0) nintpsa = 6

    call read_pseudo

    call read_fort10(ufort10)

    if(contraction.eq.0) then
        call log_info(' You should have contracted orbitals ! ')
    endif

    allocate(nhybrids(nion), psi(ipc))

    nhybrids = 0
    ntot = 0
    ind = 0

    do i = 1, nshell_c
        do j = 1, mult_c(i)
            ind = ind + 1
            if(ioptorb_c(i).ne.ioptorbcontr(ioptorb_c(i), Lbox, 1).and.ioccup_c(ind).ne.0) then
                nhybrids(kion_c(i)) = nhybrids(kion_c(i)) + ipf
                ntot = ntot + ipf
            endif
        enddo
    enddo

    call log_info('Total number of occupied contracted orbitals : ', sum(nhybrids(:)))

    call log_info('Reporting number of occupied contracted orbitals per atom')
    do i = 1, nion
        call log_info('atom : ', i, ' # of occupied contracted : ', nhybrids(i))
    enddo
    call log_info()

    do i = 1, 3
        center(i) = sum(rion(i, 1:nion)) / nion
    enddo

    if(.not.iespbc) then
        call log_info(' Choose box size (x,y,z) ')
        read(*, *) cell_loc(:)
        call log_info(cell_loc(1), cell_loc(2), cell_loc(3))
        cellscale(1:3) = cell_loc(1:3)
        do i = 1, nion
            rion(1:3, i) = rion(1:3, i) - center(1:3) + cell_loc(1:3) / 2.d0
        enddo
        ! for open system the center of the box is the (0,0,0) point
        origin = 0.d0
    else
        !call ApplyPBC(rion,nion)
        !do i=1,nion
        !rion(1:3,i)=rion(1:3,i)+cellscale(1:3)/2.d0
        !enddo
        cell_loc(1:3) = cellscale(1:3)
        call log_info(' Choose shift reference unit cell PBC : ')
        read(*, *) origin(:)
        origin(1:3) = origin(1:3) * cellscale(1:3)
    endif

    call log_info(' Choose number of mesh points (x,y,z) : ')
    read(*, *) mesh(:)
    call log_info(mesh(1), mesh(2), mesh(3))

    call log_info(' Choose ionic center between 1 and', nion)
    read(*, *) chosen_ion
    call log_info(chosen_ion)


    !spinon=.false.
    !chargeon=.false.
    !nochsp=.true.

    lower = 1
    upper = nhybrids(chosen_ion)
    molecular = nhybrids(chosen_ion)

    call log_info(' Please give a range between 1 and ', nhybrids(chosen_ion))
    read(*, *) lower, upper
    call log_info(lower, upper)

    call log_info(' Compute orbital spread and Hubbard U matrix elements')
    read(*, *) hubbardU_logic

    if(hubbardU_logic) then
        open(113, file = 'HubbardU_matrix.dat', form = 'formatted', status = 'unknown')
        call log_warning(' Warning the order of d-orbitals is the following: 3z^2-r^2,x^2-y^2,xy,yz,xz ')
    endif
    ! To match with output plot in xcrysden
    if(iespbc.and..not.hubbardU_logic) mesh(:) = mesh(:) + 1

    nradial = mesh(1) * mesh(2) * mesh(3)

    allocate(imap_loc(molecular), control(molecular), orb_center(3, molecular))

    control = .false.

    indorb = 0
    imol = 0
    ind = 0
    ind_mol = 0
    do i = 1, nshell_c
        if(ioptorb_c(i).ne.ioptorbcontr(ioptorb_c(i), Lbox, 1).and.kion_c(i).eq.chosen_ion) then
            !
            !         write(*,*) ' imol = ',imol,ind,occ_c
            !         write(*,*) ' iopcc ',ioccup_c(ind)
            !
            do j = 1, mult_c(i)
                ind = ind + 1
                if(ioccup_c(ind).eq.1) then
                    indorb = indorb + 1
                    imol = imol + 1
                    imap_loc(imol) = indorb
                    orb_center(:, imol) = rion(:, kion_c(i))
                    if(imol.ge.lower.and.imol.le.upper) then
                        control(imol) = .true.
                        ind_mol = ind_mol + 1
                    endif
                    if(ipf.eq.2) then
                        imol = imol + 1
                        imap_loc(imol) = indorb + nelorb_c / 2
                        orb_center(:, imol) = rion(:, kion_c(i))
                        if(imol.ge.lower.and.imol.le.upper) then
                            control(imol) = .true.
                            ind_mol = ind_mol + 1
                        endif
                    endif
                endif
            enddo
        else
            do j = 1, mult_c(i)
                ind = ind + 1
                if(ioccup_c(ind).eq.1) indorb = indorb + 1
            enddo
        endif
    enddo

    chosen_mol = ind_mol
    call log_info(' # of orbitals written =', chosen_mol)

    allocate(norm(chosen_mol))
    norm = 1.d0

    if(iespbc.and..not.hubbardU_logic) then
        step(:) = cell_loc(:) / dble(mesh(:) - 1)
    else
        step(:) = cell_loc(:) / dble(mesh(:))
    endif
    vol = step(1) * step(2) * step(3)

    if(hubbardU_logic.and.(step(1).ne.step(2).or.step(2).ne.step(3).or.step(1).ne.step(3))) then
        analytic = .false.
        call log_info('Computing spread and local Coulomb repulsion')
        call log_info('Integration algorithm for the U matrix: shifted grids')
        nshift = 2
    elseif(hubbardU_logic) then
        analytic = .true.
        call log_info('Computing spread and local Coulomb repulsion')
        call log_info('Integration algorithm for the U matrix: analytic integration of divergence')
    else
        call log_info('Plotting orbitals in xcrysden format file')
    endif

    if(hubbardU_logic) then
        allocate(psiri(ipc * nradial, chosen_mol), ri(3, nradial))
        psiri = 0.d0
        ri = 0.d0
        if(.not.analytic) then
            allocate(psirj(ipc * nradial, chosen_mol), rj(3, nradial))
            psirj = 0.d0
            rj = 0.d0
        endif
    else
        allocate(datagrid(ipc * mesh(1), mesh(2), mesh(3)))
        datagrid = 0.d0
    endif

    allocate(buffer_winv(ipc * nelorb))
    iflagnorm = 3
    ind_mol = 0

    do j = 1, molecular

        if(control(j)) then

            rav2 = 0.d0
            ravx = 0.d0
            ravy = 0.d0
            ravz = 0.d0
            psi2_norm = 0.d0

            ind_mol = ind_mol + 1
            ind_ri = 1

            do iii = 1, mesh(3)
                do ii = 1, mesh(2)
                    do i = 1, mesh(1)

                        !                rpoint(1)=(i-1)*step(1)-origin(1)
                        !                rpoint(2)=(ii-1)*step(2)-origin(2)
                        !                rpoint(3)=(iii-1)*step(3)-origin(3)

                        rpoint(:) = (i - 1) * step(1) * at(:, 1) + (ii - 1) * step(2) * at(:, 2) + (iii - 1) * step(3) * at(:, 3)


                        ! only the up
                        call  upnewwf(1, 0, 0, 1, nshell, ioptorb, ioccup, rpoint, 1, r, rmu, &
                                dupr, zetar, rion, psip, buffer_winv, nelorb, nion, kion, &
                                iflagnorm, cnorm, LBox, rmucos, rmusin, minimum_distance, &
                                indpar_tab, indorb_tab, indshell_tab, .true.)

                        if(n_body_on.ne.0) then
                            psilg = -scale_one_body
                            do jj = 1, nion
                                if(iespbc) then
                                    rc(:) = rpoint(:) - rion(:, jj)
                                    call CartesianToCrystal(rc, 1)
                                    do kk = 1, 3
                                        rc(kk) = costz(jj) * map(rc(kk), cellscale(kk))
                                    enddo
                                    r0 = norm_metric(rc, metric)
                                else
                                    rc(:) = (rpoint(:) - rion(:, jj)) * costz(jj)
                                    r0 = dsqrt(sum(rc(:)**2))
                                endif
                                psilg = psilg - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj)) * costz3(jj)
                            enddo
                            buffer_winv(1:ipc * nelorbh) = buffer_winv(1:ipc * nelorbh) * dexp(psilg)
                        endif

                        if(ipf.eq.1) then
                            imu = 1
                        else
                            if(imap_loc(j).le.nelorb_c / 2) then
                                imu = 1
                            else
                                imu = nelorbh + 1
                            endif
                        endif
                        if(ipc.eq.1) then
                            psi(1) = ddot(nelorbh, buffer_winv, 1, mu_c(imu, imap_loc(j)), 1)
                        else
                            psi_c = zdotu(nelorbh, buffer_winv, 1, mu_c(2 * imu - 1, imap_loc(j)), 1)
                            psi(1) = dreal(psi_c)
                            psi(2) = aimag(psi_c)
                        endif

                        if(hubbardU_logic) then

                            rvec(:) = rpoint(:) - orb_center(:, j)
                            if(iespbc) call ApplyPBC(rvec, 1)
                            r2vec = sum(rvec(:)**2)

                            psiri(ipc * (ind_ri - 1) + 1:ipc * ind_ri, ind_mol) = psi(:)
                            ri(1, ind_ri) = rvec(1)
                            ri(2, ind_ri) = rvec(2)
                            ri(3, ind_ri) = rvec(3)
                            ind_ri = ind_ri + 1
                            if(ipc.eq.1) then
                                psi2 = psi(1)**2
                            else
                                psi2 = psi_c * dconjg(psi_c)
                            endif

                            rav2 = rav2 + r2vec * psi2
                            ravx = ravx + rvec(1) * psi2
                            ravy = ravy + rvec(2) * psi2
                            ravz = ravz + rvec(3) * psi2
                            psi2_norm = psi2_norm + psi2

                        else

                            datagrid(ipc * (i - 1) + 1:ipc * i, ii, iii) = psi(:)

                        endif

                    enddo
                enddo
            enddo

            if(.not.hubbardU_logic) then
                datagrid = datagrid * norm(ind_mol)
                call plot_3d_data(ipc, cell_loc, cellscale, nion, rion, atom_number, iespbc, mesh, origin, datagrid, j, 'orbital       ')
                if(ipc.eq.2) then
                    call plot_3d_data(ipc, cell_loc, cellscale, nion, rion, atom_number, iespbc, mesh, origin, datagrid(2, 1, 1), j, 'orbitalimag   ')
                endif
                do i = 1, mesh(1)
                    datagrid(ipc * (i - 1) + 1, :, :) = datagrid(ipc * (i - 1) + 1, :, :)**2
                    if(ipc.eq.2) datagrid(ipc * (i - 1) + 1, :, :) = datagrid(ipc * (i - 1) + 1, :, :) + &
                            &datagrid(2 * i, :, :)**2
                enddo
                call plot_3d_data(ipc, cell_loc, cellscale, nion, rion, atom_number, iespbc, mesh, origin, datagrid, j, 'orbsqrd       ')
            else
                ravx = ravx / psi2_norm
                ravy = ravy / psi2_norm
                ravz = ravz / psi2_norm

                spread = rav2 / psi2_norm - ravx**2 - ravy**2 - ravz**2
                psi2_norm = vol * psi2_norm

                norm(ind_mol) = psi2_norm

                call log_info('contracted orbital #', ind_mol, 'norm (a_0**3)=', psi2_norm, 'spread (a_0**2)=', spread)
            endif

        endif

    enddo

    if(hubbardU_logic) then
        !!!!!!!!!!!!
        !integral Hubbard U

        allocate(HubbardU(chosen_mol, chosen_mol), HubbardJ(chosen_mol, chosen_mol))
        HubbardU = 0.d0
        HubbardJ = 0.d0

        if(analytic) then
            ! analytic integration of Integrate[1/Sqrt[x^2 + y^2 + z^2], {x,-1,1}, {y,-1,1}, {z,-1,1}]
            analytic_integral = -2.d0 * (pi + log(64.d0) - 12.d0 * log(1.d0 + sqrt(3.d0)))

            do i = 1, chosen_mol

                do ii = i, chosen_mol

                    U = 0.d0
                    U_ex = 0.d0

                    do iteri = 1, nradial
                        integ2 = psiri(iteri, i)**2 * psiri(iteri, ii)**2 * vol**2 / step(1) / 4.d0 * analytic_integral
                        U = U + integ2
                        U_ex = U_ex + integ2
                    enddo

                    do iteri = 1, nradial
                        do iterj = 1, iteri - 1
                            integ1 = vol**2 / sqrt((ri(1, iteri) - ri(1, iterj))**2 + (ri(2, iteri) - ri(2, iterj))**2 + (ri(3, iteri) - ri(3, iterj))**2)
                            integ2 = sum(psiri(ipc * (iteri - 1) + 1:ipc * iteri, i)**2 * psiri(ipc * (iterj - 1) + 1:&
                                    ipc * iterj, ii)**2) * integ1
                            U = U + integ2
                            integ2 = sum(psiri(ipc * (iteri - 1) + 1:ipc * iteri, i) * psiri(ipc * (iteri - 1) + 1:&
                                    &ipc * iteri, ii) * psiri(ipc * (iterj - 1) + 1:ipc * iterj, i) * &
                                    &psiri(ipc * (iterj - 1) + 1:ipc * iterj, ii)) * integ1
                            U_ex = U_ex + integ2
                        enddo

                        do iterj = iteri + 1, nradial
                            integ1 = vol**2 / sqrt((ri(1, iteri) - ri(1, iterj))**2 + (ri(2, iteri) - ri(2, iterj))**2 + (ri(3, iteri) - ri(3, iterj))**2)
                            integ2 = sum(psiri(ipc * (iteri - 1) + 1:ipc * iteri, i)**2 * psiri(ipc * (iterj - 1) + 1:&
                                    ipc * iterj, ii)**2) * integ1
                            U = U + integ2
                            integ2 = sum(psiri(ipc * (iteri - 1) + 1:ipc * iteri, i) * psiri(ipc * (iteri - 1) + 1:&
                                    &ipc * iteri, ii) * psiri(ipc * (iterj - 1) + 1:ipc * iterj, i) * &
                                    &psiri(ipc * (iterj - 1) + 1:ipc * iterj, ii)) * integ1
                            U_ex = U_ex + integ2

                            !           integ2=psiri(iteri,i)**2*psiri(iterj,ii)**2*integ1
                            !           U=U+integ2
                            !            integ2=psiri(iteri,i)*psiri(iteri,ii)*psiri(iterj,i)*psiri(iterj,ii)*integ1
                            !            U_ex=U_ex+integ2
                        enddo
                    enddo

                    HubbardU(i, ii) = U / norm(i) / norm(ii)
                    HubbardJ(i, ii) = U_ex / norm(i) / norm(ii)
                    HubbardU(ii, i) = HubbardU(i, ii)
                    HubbardJ(ii, i) = HubbardJ(i, ii)

                enddo

            enddo

        else

            ! compute orbitals on shifted grid
            ind_mol = 0

            do j = 1, molecular
                !
                if(control(j)) then

                    ind_ri = 1
                    ind_mol = ind_mol + 1

                    do iii = 1, mesh(3)
                        do ii = 1, mesh(2)
                            do i = 1, mesh(1)
                                !
                                rpoint(1) = (i - 1) * step(1) - origin(1) + step(1) / nshift
                                rpoint(2) = (ii - 1) * step(2) - origin(2) + step(2) / nshift
                                rpoint(3) = (iii - 1) * step(3) - origin(3) + step(3) / nshift
                                !
                                call  upnewwf(1, 0, 0, 1, nshell, ioptorb, ioccup, rpoint, 1, r, rmu, &
                                        dupr, zetar, rion, psip, buffer_winv, nelorb, nion, kion, &
                                        iflagnorm, cnorm, LBox, rmucos, rmusin, minimum_distance, &
                                        indpar_tab, indorb_tab, indshell_tab, .true.)

                                if(n_body_on.ne.0) then
                                    psilg = -scale_one_body
                                    do jj = 1, nion
                                        if(iespbc) then
                                            rc(:) = rpoint(:) - rion(:, jj)
                                            call CartesianToCrystal(rc, 1)
                                            do kk = 1, 3
                                                rc(kk) = costz(jj) * map(rc(kk), cellscale(kk))
                                            enddo
                                            r0 = norm_metric(rc, metric)
                                        else
                                            rc(:) = (rpoint(:) - rion(:, jj)) * costz(jj)
                                            r0 = dsqrt(sum(rc(:)**2))
                                        endif
                                        psilg = psilg - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj)) * costz3(jj)
                                    enddo
                                    buffer_winv(1:ipc * nelorbh) = buffer_winv(1:ipc * nelorbh) * dexp(psilg)
                                endif
                                if(ipf.eq.1) then
                                    imu = 1
                                else
                                    if(imap_loc(j).le.nelorb_c / 2) then
                                        imu = 1
                                    else
                                        imu = nelorbh + 1
                                    endif
                                endif

                                if(ipc.eq.1) then
                                    psi(1) = ddot(nelorbh, buffer_winv, 1, mu_c(imu, imap_loc(j)), 1)
                                else
                                    psi_c = zdotu(nelorbh, buffer_winv, 1, mu_c(2 * imu - 1, imap_loc(j)), 1)
                                    psi(1) = dreal(psi_c)
                                    psi(2) = aimag(psi_c)
                                endif

                                rvec(:) = rpoint(:) - orb_center(:, j)
                                if(iespbc) call ApplyPBC(rvec, 1)

                                psirj(ipc * (ind_ri - 1) + 1:ipc * ind_ri, ind_mol) = psi(:)
                                rj(1, ind_ri) = rvec(1)
                                rj(2, ind_ri) = rvec(2)
                                rj(3, ind_ri) = rvec(3)
                                ind_ri = ind_ri + 1
                            enddo
                        enddo
                    enddo

                endif
            enddo

            do i = 1, chosen_mol

                do ii = i, chosen_mol

                    U = 0.
                    U_ex = 0.d0
                    do iteri = 1, nradial
                        do iterj = 1, nradial
                            integ1 = 1.d0 / sqrt((ri(1, iteri) - rj(1, iterj))**2 + (ri(2, iteri) - rj(2, iterj))**2 + (ri(3, iteri) - rj(3, iterj))**2)
                            integ2 = sum(psiri(ipc * (iteri - 1) + 1:ipc * iteri, i)**2 * psirj(ipc * (iterj - 1) + 1:&
                                    ipc * iterj, ii)**2) * integ1
                            U = U + integ2
                            integ2 = sum(psiri(ipc * (iteri - 1) + 1:ipc * iteri, i) * psiri(ipc * (iteri - 1) + 1:&
                                    &ipc * iteri, ii) * psirj(ipc * (iterj - 1) + 1:ipc * iterj, i) * &
                                    &psirj(ipc * (iterj - 1) + 1:ipc * iterj, ii)) * integ1
                            U_ex = U_ex + integ2
                        enddo
                    enddo

                    HubbardU(i, ii) = U * vol**2 / norm(i) / norm(ii)
                    HubbardJ(i, ii) = U_ex * vol**2 / norm(i) / norm(ii)
                    HubbardU(ii, i) = HubbardU(i, ii)
                    HubbardJ(ii, i) = HubbardJ(i, ii)

                enddo

            enddo

        endif
        ! endif analytic

        HubbardU = HubbardU * 27.211396132
        HubbardJ = HubbardJ * 27.211396132

        write(113, *) 'U matrix'
        do i = 1, chosen_mol
            write(113, '(1000(f8.4,1x))') (HubbardU(ii, i), ii = 1, chosen_mol)
        enddo
        write(113, *)
        write(113, *) 'J matrix'
        do i = 1, chosen_mol
            write(113, '(1000(f8.4,1x))') (HubbardJ(ii, i), ii = 1, chosen_mol)
        enddo
        write(113, *)
        write(113, *) 'U-J matrix'
        do i = 1, chosen_mol
            write(113, '(1000(f8.4,1x))') ((HubbardU(ii, i) - HubbardJ(ii, i)), ii = 1, chosen_mol)
        enddo

    endif

    call deallocate_all
    close(ufort10)

    stop

end program plot_contracted
