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

subroutine compute_corr_fun(rnew, nspacedim, nparts, ell, nel_read, nelup_read&
        &, ifrho, ifspin, ifkspin, kspin, nrhoind, ddim, dxil&
        &, ind_offset, density_c, sdensity_c, mesh, vell, ifpair, npairind, ngrid_p, drmax&
        &, paircorr, spaircorr, r_offset, r_spin2, vdim, ifsofk&
        &, nvects, rhok, rhoks, pwmat, kmult, rkcomp, sofk, ifcorrs, corrsamp&
        &, seconf, spsiln, signpsi, sangle&
        &, shiftlog, logsamp, noeloc, ipart, sphere_radius, adr_spc, multcell, outofplane, rdf_for_atom)

    use constants, only: zimg, ipc
    use allio
    use ewald
    use Spin2
    use grid_module
    use Assar_module
    use Rho_corr_module, only: ifrho_corr, &
                               compute_rho_corr
    use qpwf_module, only: ifqpwf, ifqpwf_h, ifqpwf_k, ifqpwf_extr, &
                           compute_elec_delete, &
                           compute_elec_add, &
                           compute_QPWF_k_space, &
                           compute_extrapolate
    use Dipole_module, only: ifdipole, &
                             calc_dipole, &
                             calc_quadrupole
    use berry_phase
    ! calculates pair properties  of rnew
    implicit none

    integer i, nspacedim, nparts, ddim, ind, l, nrhoind, npairind, ind_offset(*), j, ii, jj, k, jp, ngrid_p, nspec
    integer nvects, kmult(0:*), ipart(*), adr_spc(*)
    real*8 rnew(nspacedim, *), ell(*), r_offset(*), r_spin2(*), vell(*)
    real*8 dxil(*), density_c(*), sdensity_c(*), paircorr(*), spaircorr(*), corrsamp(*)
    integer mesh(3), indcheck, iix(3), vdim(*), point(3), nelorbjmax, indtmax, nshelljmax
    real*8 drmax
    complex*16 phase_
    logical ifrho, ifspin, ifkspin, control, ifpair, ifsofk, ifcorrs, noeloc
    real*8 spin, rx(3), ry(3), dx(3), shiftlog, logsamp, kspin(3), scalcost, rmod
    real*8 zero(3), sofk(nvects, *), rhok(*), rhoks(2*nvects, *), pwmat(2*nvects, *), rkcomp(nspacedim, *)
    integer kp, nppss(2)
    real*8 spsiln, seconf
    real*4 sangle(18, *), signpsi
    real*8 vcutav(2), voffpseudo, psi_new, sphere_radius
    real*8, external :: enercont
    logical :: yesprojm, check_inside
    real*4 dx_test(3)
    real*8 eps, elocm(2), logpsim(2), outofplane
    !allind::he first column contains the box in wich the electron is, the second one is
    !the number of the electron (necessary for the dsort output understanding)

    !listf::each column corresponds to an occupied box (dimention is dimlist that is
    !nel (or nparts) if the calculation is on the mesh, else nion), the first number
    !is the index relative to the box, the second is the number of up electron inside
    !the box, the third is the number of down electron in the box, starting from
    !the fourth we have the electrons
    integer :: dimlist, electron, newbox, elec
    real(8), allocatable :: allind(:, :), test1(:), test2(:)
    integer, allocatable :: listf(:, :)
    integer nel_read, nelup_read, multcell
    ! (Kosuke Nakano) an logical variable for computing Radial distribution function of atom
    logical rdf_for_atom
    real*8 r_dist

    parameter(eps=1d-5)
    !    defines once for all the parameters for upscratch
    iflagnorm = 3
    if (allocated(psiln)) psiln = 0.d0
    if (allocated(psiln)) singdet = .true.
    test_aad = .true.
    pseudologic = .false.
    iesrandoml = .false.
    itestrfn = 2
    itest = 2
    pseudologic = .false.
    iesrandoml = .false.
    nshelljmax = max(nshellj, 1)
    neldomax = max(neldo, 1)
    indtmax = max(indt, 1)
    nelorbjmax = max(nelorbj, 1)
    if (contraction .ne. 0) then
        yesprojm = .true.
    else
        yesprojm = .false.
    end if
    zero = 0.d0

    nppss(1) = nelup
    nppss(2) = nparts - nelup
    if (ifberry) then
        call berry_phase_update(rnew, nparts)

    end if

    if (ifrho) then

        do i = 1, nrhoind
            density_c(i) = 0.d0
        end do

        if (sphere_radius .eq. 0.d0) then

            do i = 1, nparts

                ind = 1

                check_inside = .true.
                do l = 1, ddim
                    ll = vdim(l)
                    if (.not. rdf_for_atom) then
                        indcheck = int(rnew(ll, i)*dxil(l))
                    else
                        r_dist = sqrt(sum(rnew(:, i)**2)) ! newly added to compute rdf
                        indcheck = int(r_dist*dxil(l)) ! newly added to compute rdf
                    end if
                    ! protection against roundoff misplacements
                    if (indcheck .lt. 0) then
                        indcheck = mesh(l) - 1
                        check_inside = .false.
                    elseif (indcheck .gt. mesh(l) - 1) then
                        indcheck = 0
                        check_inside = .false.
                    end if
                    ind = ind + indcheck*ind_offset(l)
                end do

                if (check_inside) then
                    density_c(ind) = density_c(ind) + 1.d0/multcell
                end if
            end do

        else

            !This function calculates the staggered dipole in a chain along Z,
            !The staggered dipole is the dypole calculated in every site times a
            ! (-1)**k
            if (ddim == 1) then
                call el_ion_distance_for_dipole(dists_kel, rnew, nparts, rion, nion, iespbc) !La subroutine porta il nome del file
                do i = 1, nparts
                    do k = 1, nion
                        !It makes the calculation only if it is in the considered region
                        if (abs(dists_kel(k, i)) .lt. abs(sphere_radius)) then
                            if (sphere_radius .lt. 0.d0) then
                                density_c(k) = density_c(k) + dists_kel(k, i)/multcell
                            else
                                density_c(k) = density_c(k) + 1.d0/multcell
                            end if
                        end if
                    end do
                end do

            else
                !! Here we start calculating density of charge

                if (ddim .eq. 2) then
                    call el_ion_distance_2d(dists_kel, rnew, nparts, rion, nion, iespbc, outofplane)
                else
                    call el_ion_distance(dists_kel, rnew, nparts, rion, nion, iespbc)
                end if
                do i = 1, nparts
                    do k = 1, nion
                        if (dists_kel(k, i) .lt. abs(sphere_radius)) then
                            density_c(k) = density_c(k) + 1.d0/multcell
                        end if
                    end do
                end do
            end if
        end if

        !This seems to be the pair correlation function for the ions

        if (ifpair) then

            do i = 1, npairind + ngrid_p
                paircorr(i) = 0.d0
            end do

            if (sphere_radius .eq. 0.d0) then

                ! compute pair_corr_fun
                do i = 1, nparts
                    do k = i, nparts

                        ! r(i)-r(k) in PBC
                        do l = 1, ddim
                            dx_test(l) = rnew(l, i) - rnew(l, k)
                        end do
                        call MyApplyPBC(dx_test, 1, ell, r_offset, iespbc)

                        ! find ind from dx_test
                        ind = 1
                        check_inside = .true.
                        do l = 1, ddim
                            ll = vdim(l)
                            indcheck = int(dble(dx_test(ll))*dxil(l)) + 1
                            ! protection against roundoff misplacements
                            if (indcheck .gt. mesh(l)) then
                                check_inside = .false.
                                indcheck = 1
                            elseif (indcheck .lt. 1) then
                                indcheck = mesh(l)
                                check_inside = .false.
                            end if
                            ind = ind + (indcheck - 1)*ind_offset(l)
                        end do
                        if (check_inside) then

                            if (i .ne. k) then
                                paircorr(ind) = paircorr(ind) + 2.d0
                            else
                                paircorr(ind) = paircorr(ind) + 1.d0
                            end if

                            ! dx_test modulus
                            rmod = 0.d0
                            do l = 1, ddim
                                ll = vdim(l)
                                rmod = rmod + dble(dx_test(ll))**2
                            end do
                            rmod = sqrt(rmod)

                            ! find ind from rmod
                            ind = int(rmod*drmax) + 1
                            !              ! protection against roundoff misplacements
                            !              if(ind.gt.ngrid_p) then
                            !                 ind=ngrid_p
                            !              elseif(ind.lt.1) then
                            !                 ind=1
                            !              endif

                            ind = ind + npairind

                            if (i .ne. k) then
                                paircorr(ind) = paircorr(ind) + 2.d0
                                !              else
                                !                 paircorr(ind)=paircorr(ind)+1.d0
                            end if
                        end if

                    end do
                end do

            else

                do i = 1, nion
                    do k = 1, nion
                        !
                        ! general pair correlation function. No distinction between
                        ! atomic species, contained in paircorr(1:nshell)
                        ind = k + (i - 1)*nion
                        ind = ipip_adr(ind)
                        paircorr(ind) = paircorr(ind) + density_c(i)*density_c(k)
                        !
                        ! pair correlation function for each species sector
                        ! contained in paircorr(nshell+1:nshell+nshellsp)
                        ind = k + (i - 1)*nion
                        ind = adr_spc(ind)
                        paircorr(ind) = paircorr(ind) + density_c(i)*density_c(k)

                    end do
                end do
                !
                !If ddim=1 the 0 term is equal to the dipole dipole correlation
                if (ddim == 1) then
                    paircorr(1) = 0
                    do i = 1, nion
                        paircorr(1) = paircorr(1) + (-1)**(i)*density_c(i)
                    end do
                    paircorr(1) = (paircorr(1)/nion)**2*mult(1)
                end if
            end if

        end if

    end if

    !if(rank.eq.0) write(6,*) ' after rho '

    if (ifspin) then

        !This calculate the spin2 and all the related quantities
        !or inside the sphere_radius the actual on a mesh
        !calculation for the spin calculations follow this if.
        if (ifspin2) then

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Preliminary part that calculates the ratiospin matrix necessary for the
            !following calculation
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, rnew, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, contraction, elocm, logpsim&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            call compute_spin2(kel, ainv, winv, winvj, winvbar, winvjbar, winvjbarsz, detmat, projm, mu_c, nmolfn)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do i = 1, nrhoind + 1
                sdensity_c(i) = 0.d0
            end do

            if (sphere_radius .eq. 0.d0) then
                dimlist = nparts
            else
                dimlist = nion
            end if

            allocate (allind(nparts, 2), listf(nparts + 3, dimlist))
            allind = 0.d0

            !Calculating allind (on the mesh here, with sphere_radius after)
            if (sphere_radius .eq. 0.d0) then

                call DDMyApplyPBC(rnew, nparts, cellscale, r_spin2, iespbc)

                do i = 1, nparts
                    ind = 1
                    check_inside = .true.
                    do l = 1, ddim
                        ll = vdim(l)
                        indcheck = int(rnew(ll, i)*dxil(l))
                        ! protection against roundoff misplacements
                        if (indcheck .lt. 0) then
                            check_inside = .false.
                            indcheck = mesh(l) - 1
                        elseif (indcheck .gt. mesh(l) - 1) then
                            check_inside = .false.
                            indcheck = 0
                        end if
                        ind = ind + indcheck*ind_offset(l)
                    end do

                    if (check_inside) then
                        allind(i, 2) = dble(i)
                        allind(i, 1) = dble(ind)
                    end if
                end do
            else
                !Starting the calculation of allind on the spheres
                if (.not. ifrho) then
                    if (ddim == 1) then
                        call el_ion_distance_for_dipole(dists_kel, rnew, nparts, rion, nion, iespbc) !
                    elseif (ddim .eq. 2) then
                        call el_ion_distance_2d(dists_kel, rnew, nparts, rion, nion, iespbc, outofplane)
                    else
                        call el_ion_distance(dists_kel, rnew, nparts, rion, nion, iespbc)
                    end if
                end if
                do i = 1, nparts
                    do k = 1, nion
                        if (abs(dists_kel(k, i)) .lt. abs(sphere_radius)) then
                            allind(i, 1) = dble(k)
                            allind(i, 2) = dble(i)
                        end if
                    end do
                end do

            end if
            !        do i=1,nparts
            !           write(6,*) "allind", i,allind(i,1),allind(i,2)
            !        enddo
            !Once we have allind we can create listf
            call dsort(allind(1, 1), allind(1, 2), nparts, 2)

            listf = 0
            !In the case in which no electron is in the boxes
            newbox = 0
            !It creates the list only if there is at least one electron inside the boxes
            if (int(allind(nparts, 2)) .ne. 0) then
                !New box is the last box in which we put an electron
                newbox = 1
                do i = 1, nparts
                    elec = int(allind(i, 2))
                    ind = int(allind(i, 1))
                    if (elec .ne. 0) then
                        !I change the column if I have finished with that box
                        if (listf(1, newbox) .ne. ind) then
                            if (listf(1, newbox) .eq. 0) then
                                listf(1, newbox) = ind
                            else
                                newbox = newbox + 1
                                listf(1, newbox) = ind
                            end if
                        end if
                        !For the moment I put all the down electrons at the end of the row
                        !then I will move then after the up ones (in this way I can save
                        !time loading from the memory still having all the electron in the
                        !right order!)
                        if (elec .gt. nelup) then
                            listf(3, newbox) = listf(3, newbox) + 1
                            listf(4 + nparts - listf(3, newbox), newbox) = elec
                        else
                            listf(2, newbox) = listf(2, newbox) + 1
                            listf(3 + listf(2, newbox), newbox) = elec
                        end if
                    end if
                end do
                !Putting the down electrons in their positions
                do i = 1, newbox
                    if (listf(3, i) .ne. 0 .and. listf(3, i) + listf(2, i) .ne. nparts) then
                        do l = 1, listf(3, i)
                            listf(3 + listf(2, i) + l, i) = listf(3 + nparts - listf(3, i) + l, i)
                            listf(3 + nparts - listf(3, i) + l, i) = 0
                        end do
                    end if
                end do
            end if

            !        do i=1,nparts
            !           write(6,*) "after allind", i,allind(i,1),allind(i,2)
            !        enddo
            !        write (6,*) "nelup=", nelup
            !        write (6,*) "listf", listf(1:3,1)
            !        write (6,*) "listf", listf(4:3+nparts,1)
            !        write (6,*) "listf", listf(1:3,2)
            !        write (6,*) "listf", listf(4:3+nel,2)
            !        call mpi_finalize(ierr)
            !        stop

            !Now listf is ready and we can use it
            if (newbox .gt. 0) then
                do i = 1, newbox
                    sdensity_c(listf(1, i)) = 0.5d0*dble(listf(2, i) + listf(3, i)) + 0.25d0*dble(listf(2, i) - listf(3, i))**2
                    if (listf(2, i) .ne. 0 .and. listf(3, i) .ne. 0) then
                        do l = 1, listf(3, i)
                            do ll = 1, listf(2, i)
                                sdensity_c(listf(1, i)) = sdensity_c(listf(1, i)) &
                                                          + ratiospin(listf(3 + ll, i), listf(3 + listf(2, i) + l, i) - nelup)
                            end do
                        end do
                    end if
                end do
            end if
            deallocate (allind, listf)

            !Starting the regular spin calculation
        else
            do i = 1, nrhoind + 1
                sdensity_c(i) = 0.d0
            end do
            if (ifkspin) sdensity_c(nrhoind + 2) = 0.d0

            if (sphere_radius .eq. 0.d0) then
                do i = 1, nelup
                    ind = 1
                    check_inside = .true.
                    do l = 1, ddim
                        ll = vdim(l)
                        indcheck = int(rnew(ll, i)*dxil(l))
                        ! protection against roundoff misplacements
                        if (indcheck .lt. 0) then
                            check_inside = .false.
                            indcheck = mesh(l) - 1
                        elseif (indcheck .gt. mesh(l) - 1) then
                            check_inside = .false.
                            indcheck = 0
                        end if
                        ind = ind + indcheck*ind_offset(l)
                    end do

                    if (check_inside) then
                        sdensity_c(ind) = sdensity_c(ind) + 0.5d0/multcell
                        if (ifkspin) then
                            scalcost = sum(kspin(1:3)*rnew(1:3, i))
                            sdensity_c(nrhoind + 1) = sdensity_c(nrhoind + 1) + dcos(scalcost)
                            sdensity_c(nrhoind + 2) = sdensity_c(nrhoind + 2) + dsin(scalcost)
                        end if
                    end if
                end do

                do i = nelup + 1, nparts
                    control = .true.
                    ind = 1
                    check_inside = .true.
                    do l = 1, ddim
                        ll = vdim(l)
                        indcheck = int(rnew(ll, i)*dxil(l))
                        ! protection against roundoff misplacements
                        if (indcheck .lt. 0) then
                            indcheck = mesh(l) - 1
                            check_inside = .false.
                        elseif (indcheck .gt. mesh(l) - 1) then
                            indcheck = 0
                            check_inside = .false.
                        end if
                        ind = ind + indcheck*ind_offset(l)
                    end do
                    if (check_inside) then
                        sdensity_c(ind) = sdensity_c(ind) - 0.5d0/multcell
                        if (ifkspin) then
                            scalcost = sum(kspin(1:3)*rnew(1:3, i))
                            sdensity_c(nrhoind + 1) = sdensity_c(nrhoind + 1) - dcos(scalcost)
                            sdensity_c(nrhoind + 2) = sdensity_c(nrhoind + 2) - dsin(scalcost)
                        end if
                    end if
                end do

            else

                if (.not. ifrho) then
                    if (ddim == 1) then
                        call el_ion_distance_for_dipole(dists_kel, rnew, nparts, rion, nion, iespbc) !
                    elseif (ddim .eq. 2) then
                        call el_ion_distance_2d(dists_kel, rnew, nparts, rion, nion, iespbc, outofplane)
                    else
                        call el_ion_distance(dists_kel, rnew, nparts, rion, nion, iespbc)
                    end if
                end if
                if (ddim == 1) then
                    do i = 1, nelup
                        do k = 1, nion
                            if (abs(dists_kel(k, i)) .lt. abs(sphere_radius)) then
                                if (sphere_radius .lt. 0.d0) then
                                    sdensity_c(k) = sdensity_c(k) + 0.5d0*dists_kel(k, i)/multcell
                                else
                                    sdensity_c(k) = sdensity_c(k) + 0.5d0/multcell
                                end if
                            end if
                        end do
                    end do

                    do i = nelup + 1, nparts
                        do k = 1, nion
                            if (abs(dists_kel(k, i)) .lt. abs(sphere_radius)) then
                                if (sphere_radius .lt. 0.d0) then
                                    sdensity_c(k) = sdensity_c(k) - dists_kel(k, i)*0.5d0/multcell
                                else
                                    sdensity_c(k) = sdensity_c(k) - 0.5d0/multcell

                                end if
                            end if
                        end do
                    end do

                else

                    do i = 1, nelup
                        do k = 1, nion
                            if (dists_kel(k, i) .lt. abs(sphere_radius)) then
                                sdensity_c(k) = sdensity_c(k) + 0.5d0/multcell
                            end if
                        end do
                    end do

                    do i = nelup + 1, nparts
                        do k = 1, nion
                            if (dists_kel(k, i) .lt. abs(sphere_radius)) then
                                sdensity_c(k) = sdensity_c(k) - 0.5d0/multcell
                            end if
                        end do
                    end do
                end if

            end if

            if (ifpair) then

                do i = 1, npairind + ngrid_p
                    spaircorr(i) = 0.d0
                end do

                if (sphere_radius .eq. 0.d0) then

                    ipart(nelup + 1:nparts) = -1
                    ipart(1:nelup) = 1

                    do i = 1, nparts
                        do k = i, nparts
                            ! r(i)-r(k)
                            do l = 1, ddim
                                dx_test(l) = rnew(l, i) - rnew(l, k)
                            end do
                            call MyApplyPBC(dx_test, 1, ell, r_offset, iespbc)

                            ! find ind from dx_test
                            ind = 1
                            check_inside = .true.
                            do l = 1, ddim
                                ll = vdim(l)
                                indcheck = int(dble(dx_test(ll))*dxil(l)) + 1
                                ! protection against roundoff misplacements
                                if (indcheck .gt. mesh(l)) then
                                    indcheck = 1
                                    check_inside = .false.
                                elseif (indcheck .lt. 1) then
                                    indcheck = mesh(l)
                                    check_inside = .false.
                                end if
                                ind = ind + (indcheck - 1)*ind_offset(l)
                            end do

                            if (check_inside) then

                                if (i .ne. k) then
                                    spaircorr(ind) = spaircorr(ind) + 2.d0*sign(1, ipart(k))*sign(1, ipart(i))
                                else
                                    spaircorr(ind) = spaircorr(ind) + sign(1, ipart(k))*sign(1, ipart(i))
                                end if

                                ! dx_test modulus
                                rmod = 0.d0
                                do l = 1, ddim
                                    ll = vdim(l)
                                    rmod = rmod + dble(dx_test(ll))**2
                                end do
                                rmod = sqrt(rmod)

                                ! find ind from rmod
                                ind = int(rmod*drmax) + 1
                                !              ! protection against roundoff misplacements
                                !              if(ind.gt.ngrid_p) then
                                !                 ind=ngrid_p
                                !              elseif(ind.lt.1) then
                                !                 ind=1
                                !              endif

                                ind = ind + npairind
                                if (i .ne. k) then
                                    spaircorr(ind) = spaircorr(ind) + 2.d0*sign(1, ipart(k))*sign(1, ipart(i))
                                    !              else
                                    !                 spaircorr(ind)=spaircorr(ind)+sign(1,ipart(k))*sign(1,ipart(i))
                                end if
                            end if

                        end do
                    end do

                else
                    !Beginning of the calculation for the pair correlation function
                    do i = 1, nion
                        do k = 1, nion

                            ind = k + (i - 1)*nion
                            ind = ipip_adr(ind)
                            spaircorr(ind) = spaircorr(ind) + sdensity_c(i)*sdensity_c(k)
                            !
                            ind = k + (i - 1)*nion
                            ind = adr_spc(ind)
                            spaircorr(ind) = spaircorr(ind) + sdensity_c(i)*sdensity_c(k)

                        end do
                    end do
                end if
            end if
        end if

    end if

    !write(6,*) ' if so ',ifsofk

    if (ifsofk) then
        call CompPW(rnew, nspacedim, nparts, nvects, rkcomp, pwmat, ddim, vdim)
        call CompRhok(rhok, rhoks, nparts, nvects, pwmat, nppss)
        do k = 1, nvects
            do i = 1, 5
                sofk(k, i) = 0.d0
            end do
            do i = 1, 2 ! i=1 --> cos, i=2 --> sin
                kp = 2*(k - 1) + i
                sofk(k, 1) = sofk(k, 1) + rhok(kp)*rhok(kp) ! charge-charge
                sofk(k, 2) = sofk(k, 2) + (rhoks(kp, 1) - rhoks(kp, 2))*(rhoks(kp, 1) - rhoks(kp, 2)) ! spin-spin
                sofk(k, 3) = sofk(k, 3) + rhoks(kp, 1)*rhoks(kp, 1) !up,up
                sofk(k, 4) = sofk(k, 4) + rhoks(kp, 2)*rhoks(kp, 2) !down,down
                sofk(k, 5) = sofk(k, 5) + rhoks(kp, 1)*rhoks(kp, 2) !up,down
            end do
        end do
    end if

    if (ifrho_assar) then !(Matteo) calculating observabls for Assaraf charge density_c

        !    do i=1,nel
        !       do k=1,nspacedim
        !          kel(k,i)=rnew(k,i)
        !       enddo
        !    enddo

        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                &, tabpip, rnew, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                &, nintpseudo, prefactor, rcutoff, parshell                            &
                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                &, jpseudo, pseudolocal, istart, costz, costz3         &
                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                &, psidetln, 1, nelorbh, nelorbjh                         &
                &, niesd&
                &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, contraction, elocm, logpsim&
                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

        call compute_rho_assar(dist, rion, nion, kel, nelup, neldo, indt, tabpip, winvup, winvdo, density_c, zetaq)
    end if

    if (ifrho_corr) then
        call compute_rho_corr(rnew, nel_read, nelup_read, density_c, signpsi, spsiln, shiftlog)
    end if

    !(Matteo) for STM
    ! Calculate the STM contribution from the creation operator.
    ! In this case the loop is over the grid points.
    ! There is no error due to the discretization effect.
    if (ifqpwf .and. .not. ifqpwf_h .and. .not. ifqpwf_k .and. .not. ifqpwf_extr) then
        call compute_elec_add(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog)
        ! Calculate the STM contribution from the destruction operator.
        ! In this case the loop is over the electrons and the value of the grid
        ! must be allocated.
        ! There is an error due to the discretization of the space.
    elseif (ifqpwf .and. ifqpwf_h) then
        call compute_elec_delete(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog, mesh, vell)
    elseif (ifqpwf .and. ifqpwf_k) then
        call compute_QPWF_k_space(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog)
    elseif (ifqpwf .and. ifqpwf_extr) then
        call compute_extrapolate(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog)
    end if

    if (ifcorrs) then
        !
        if (pseudorandom) then
            do i = 1, nel
                do k = 1, 18
                    angle(k, i) = sangle(k, i)
                end do
            end do

            do j = 1, nel
                do i = 1, istart - 1
                    call dgemv('N', 3, 3, 1.d0, angle(1, j), 3, versoralat(1, i), 1, 0.d0, ivic(1, i, j), 1)
                end do
            end do
        end if

        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                &, tabpip, rnew, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                &, nintpseudo, prefactor, rcutoff, parshell                            &
                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                &, jpseudo, pseudolocal, istart, costz, costz3         &
                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                &, psidetln, 1, nelorbh, nelorbjh                         &
                &, niesd&
                &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, contraction, elocm, logpsim&
                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

        !    call  upscratch(indt,indt4,indt4j,nelorb,nelup,neldo             &
        !         &,tabpip,diag,kel,winv             &
        !         &,winvup,winvdo,ainv,ainvup,ainvdo                                 &
        !         &,psip,ipsip,wconf,wconfn,psiln,psisn,iesdr,vj,dupr                &
        !         &,zetar,rion,dist,ioccup,ioccdo,ioptorb,nshell                     &
        !         &,nshelldo,ivic,alat,plat,vpot                                     &
        !         &,tmu,itest,bcost,ccost,nion,r,rmu,kion,iond                       &
        !         &,winvj,ioccj,kionj,vjur,nelorbj,ioptorbj,nshellj                  &
        !         &,winvbar,detmat,winvjbar,winvjbarsz                               &
        !         &,jasmat,jasmatsz,iessz,iflagnorm,cnorm,iflagerr,kl,parcut            &
        !         &,npsa,lmax,nintpseudo,prefactor                                  &
        !         &,rcutoff,parshell,nparpshell,kindion,pshell,wpseudo,legendre      &
        !         &,versor,wintpseudo,jpseudo,pseudolocal,coeff,istart               &
        !         &,costz,costz3,.false.,angle,indtm                     &
        !         &,LBox,rmucos,rmusin,kappa,vpotreg,cutreg                       &
        !         &,selfsum,psidetln,1                                               &
        !         &,nelorbh,nelorbjh                                                 &
        !         &,jastrowall_ee,jastrowall_ei,niesd,.false.,versoralat,iond_cart   &
        !         &,mu_c,projm,nelorb_c,firstmol,nmolfn,yesfast&
        !         &,winvbarfn,winvfn,contraction,detmat_c                            &
        !         &,indpar_tab,indorb_tab,indshell_tab,indparj_tab,indorbj_tab,indshellj_tab)
        if (iflagerr .eq. 0) then

            econtnew = elocm(1)
            psi_new = logpsim(1)

            corrsamp(2) = seconf*0.5d0

            if (noeloc) then
                ! By neglecting the variation of the local energy one makes a negligible error.
                !
                corrsamp(3) = seconf*dexp(2.d0*(psi_new - spsiln - shiftlog))*0.5d0
            else
                corrsamp(3) = econtnew*dexp(2.d0*(psi_new - spsiln - shiftlog))*0.5d0
            end if
            logsamp = psi_new - spsiln
            if (ipc .eq. 2) then
                phase_ = exp(zimg*(psisn(1) - signpsi))
                corrsamp(1) = dexp(psi_new - spsiln - shiftlog)*real(phase_)
            else
                corrsamp(1) = dexp(psi_new - spsiln - shiftlog)*psisn(1)
            end if
            if (ipc .eq. 2) then
                corrsamp(5) = dexp(psi_new - spsiln - shiftlog)*aimag(phase_)
                !    add the correction due to the sign read
            else
                corrsamp(1) = corrsamp(1)*signpsi
            end if
            corrsamp(4) = corrsamp(1)**2
            if (ipc .eq. 2) corrsamp(4) = corrsamp(4) + corrsamp(5)**2
        else
            iflagerr = 0 ! to continue the run
            corrsamp(1:5) = 0.d0
        end if
    end if

    if (ifspin2 .and. .not. ifspin) then

        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                &, tabpip, rnew, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                &, nintpseudo, prefactor, rcutoff, parshell                            &
                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                &, jpseudo, pseudolocal, istart, costz, costz3         &
                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                &, psidetln, 1, nelorbh, nelorbjh                         &
                &, niesd&
                &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, contraction, elocm, logpsim&
                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

        call compute_spin2(kel, ainv, winv, winvj, winvbar, winvjbar, winvjbarsz, detmat, projm, mu_c, nmolfn)

    end if

    if (ifdipole) then
        call calc_dipole(rnew, vdim, ell)
        call calc_quadrupole(rnew, vdim, ell)
    end if
    return
end subroutine compute_corr_fun
