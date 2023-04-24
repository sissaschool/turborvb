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

module qpwf_module
    use constants, only: zimg

    logical :: ifqpwf, ifqpwf_h, ifqpwf_k, ifqpwf_extr, decouple_k, decouple_files
    integer :: n_extr_points
    real*8 :: extr_point(3)
    real*8, dimension(:), allocatable :: qpwf_image
    logical, private :: yesprojm
    integer, private :: nelorbjmax, indtmax, neldomax, nshelljmax
    real*8, private :: logpsim(2), elocm(2)

contains

    subroutine compute_elec_add(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog)
        use allio
        use grid_module
        implicit none

        integer grid_point
        double precision rnew(3, *)
        double precision shiftlog, spsiln

        complex*16 phase_

        integer nel_read, nelup_read

        real*4 signpsi
        integer i, j

        !    defines once for all the parameters for upscratch
        iflagnorm = 3
        test_AAD = .false.
        pseudologic = .false.
        iesrandoml = .false.
        itestrfn = 2
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

        itest = 2

        if (nelup_read .lt. nelup) then
            do i = 1, nelup_read
                kel(:, i) = rnew(:, i)
            end do
            kel(:, nelup) = out_grid(:, 1)
            do i = nelup + 1, nel
                kel(:, i) = rnew(:, i - 1)
            end do
        else
            do i = 1, nel_read
                kel(:, i) = rnew(:, i)
            end do
            kel(:, nel) = out_grid(:, 1)
        end if

        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                &, tabpip, kel, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                &, nintpseudo, prefactor, rcutoff, parshell                            &
                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                &, jpseudo, pseudolocal, istart, costz, costz3         &
                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                &, psidetln, 1, nelorbh, nelorbjh                         &
                &, niesd&
                &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

        if (ipc .eq. 1) then
            qpwf_image(1) = signpsi*dexp(psiln(1) - spsiln - shiftlog)*psisn(1)
        else
            phase_ = exp(zimg*(psisn(1) - signpsi))
            qpwf_image(1) = dexp(psiln(1) - spsiln - shiftlog)*real(phase_)
            qpwf_image(2) = dexp(psiln(1) - spsiln - shiftlog)*aimag(phase_)
        end if

        !    write(*,*) psiln(1)-spsiln

        do grid_point = 2, grid_points
            if (nelup_read .lt. nelup) then
                rcart(:) = out_grid(:, grid_point)
                iout = nelup
            else
                rcart(:) = out_grid(:, grid_point)
                iout = nel
            end if

            call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                    &, ainv, kel, rcart, iesdr, vj, dupr&
                    &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr&
                    &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur&
                    &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar&
                    &, winvjbar, winvjbarsz, iflagnorm, cnorm, iflagerr&
                    &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne&
                    &, jastrowall_ee(1, iout, 0, 1), jasnew_ee, jastrowall_ei(1, iout, 1)&
                    &, jasnew_ei, n_body_on, niesd, nshell, 0&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            if (ipc .eq. 1) then
                qpwf_image(grid_point) = qpwf_image(1)*ratior(1)
            else
                qpwf_image(2*grid_point - 1) = qpwf_image(1)*ratior(1) - qpwf_image(2)*ratior(2)
                qpwf_image(2*grid_point) = qpwf_image(1)*ratior(2) + qpwf_image(2)*ratior(1)
            end if
        end do

    end subroutine compute_elec_add

    subroutine compute_extrapolate(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog)
        use allio
        use grid_module, only: da
        implicit none

        logical ele_first
        double precision rnew(3, *)
        double precision shiftlog, spsiln
        double precision stm_image_first
        double precision delta(3)
        real*4 signpsi
        integer nel_read, nelup_read
        integer :: i, j

        delta(:) = 0.d0

        psisn = 0
        itest = 2
        iflagnorm = 3
        !    psiln=0.d0
        singdet = .true.
        test_AAD = .false.
        pseudologic = .false.
        iesrandoml = .false.
        itestrfn = 2
        nshelljmax = max(nshellj, 1)
        neldomax = max(neldo, 1)
        indtmax = max(indt, 1)
        nelorbjmax = max(nelorbj, 1)
        if (contraction .ne. 0) then
            yesprojm = .true.
        else
            yesprojm = .false.
        end if

        ele_first = .true.

        do i = 1, n_extr_points
            qpwf_image(i) = 0.0
        end do

        if (nelup_read .gt. nelup) then

            do i = 1, nelup_read

                delta(1) = abs(rnew(1, i) - extr_point(1))
                delta(2) = abs(rnew(2, i) - extr_point(2))
                delta(3) = abs(rnew(3, i) - extr_point(3))

                if (delta(1) .le. da(1)/2 .and. delta(2) .le. da(2)/2 .and. delta(3) .le. da(3)/2) then

                    if (ele_first) then

                        ele_first = .false.

                        rcart(:) = rnew(:, i)
                        do j = 1, i - 1
                            kel(:, j) = rnew(:, j)
                        end do
                        do j = i + 1, nel_read
                            kel(:, j - 1) = rnew(:, j)
                        end do

                        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo, tabpip, kel, &
                                                 kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip, &
                                                 ipsip, wconfn, psisn, iesdr, vj, dupr, zetar, rion, dist, &
                                                 ioccup, ioccdo, ioptorb, nshell, nshelldo, ivic, alat, &
                                                 plat, vpot, tmu, nion, r, rmu, kion, iond, winvj, ioccj, kionj, &
                                                 vjur, nelorbj, ioptorbj, nshellj, winvbar, detmat, winvjbar, &
                                                 winvjbarsz, jasmat, jasmatsz, muj_c, jasmat_c, jasmatsz_c, &
                                                 contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax, &
                                                 nintpseudo, prefactor, rcutoff, parshell, nparpshell, &
                                                 kindion, pshell, wpseudo, legendre, versor, wintpseudo, &
                                                 jpseudo, pseudolocal, istart, costz, costz3, angle, indtm, &
                                                 LBox, rmucos, rmusin, kappa, vpotreg, cutreg, psidetln, 1, &
                                                 nelorbh, nelorbjh, niesd, iond_cart, mu_c, detmat_c, projm, &
                                                 yesprojm, nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim, &
                                                 nelorbjmax, neldomax, indtmax, nshelljmax, cellscale, &
                                                 indpar_tab, indorb_tab, indshell_tab, indparj_tab, &
                                                 indorbj_tab, indshellj_tab)

                        stm_image_first = signpsi*dexp(psiln(1) - spsiln - shiftlog)*psisn(1)*(-1)**(i - nelup_read)
                        qpwf_image(n_extr_points) = stm_image_first
                        do j = 1, n_extr_points - 1
                            if (delta(1) .le. (da(1) - da(1)/n_extr_points*j)/2 .and. &
                                delta(2) .le. (da(2) - da(2)/n_extr_points*j)/2 .and. &
                                delta(3) .le. (da(3) - da(3)/n_extr_points*j)/2) then
                                qpwf_image(n_extr_points - j) = qpwf_image(n_extr_points - j) + stm_image_first
                            else
                                exit
                            end if
                        end do
                    else
                        iout = i - 1

                        call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                                &, ainv, kel, rcart, iesdr, vj, dupr&
                                &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr&
                                &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur&
                                &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar&
                                &, winvjbar, winvjbarsz, iflagnorm, cnorm, iflagerr&
                                &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne&
                                &, jastrowall_ee(1, iout, 0, 1), jasnew_ee, jastrowall_ei(1, iout, 1)&
                                &, jasnew_ei, n_body_on, niesd, nshell, 0&
                                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                        qpwf_image(n_extr_points) = qpwf_image(n_extr_points) + stm_image_first*ratior(1)*(-1)
                        do j = 1, n_extr_points - 1
                            if (delta(1) .le. (da(1) - da(1)/n_extr_points*j)/2 .and. &
                                delta(2) .le. (da(2) - da(2)/n_extr_points*j)/2 .and. &
                                delta(3) .le. (da(3) - da(3)/n_extr_points*j)/2) then
                                qpwf_image(n_extr_points - j) = qpwf_image(n_extr_points - j) + stm_image_first*ratior(1)*(-1)
                            else
                                exit
                            end if
                        end do
                    end if
                end if
            end do
        else
            do i = nelup_read + 1, nel_read

                delta(1) = abs(rnew(1, i) - extr_point(1))
                delta(2) = abs(rnew(2, i) - extr_point(2))
                delta(3) = abs(rnew(3, i) - extr_point(3))

                if (delta(1) .le. da(1)/2 .and. delta(2) .le. da(2)/2 .and. delta(3) .le. da(3)/2) then

                    if (ele_first) then

                        ele_first = .false.
                        rcart(:) = rnew(:, i)

                        do j = 1, i - 1
                            kel(:, j) = rnew(:, j)
                        end do
                        do j = i + 1, nel_read
                            kel(:, j - 1) = rnew(:, j)
                        end do

                        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                                &, tabpip, kel, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                                &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                                &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                                &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                                &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                                &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm&
                                &, iflagerr, npsa, lmax &
                                &, nintpseudo, prefactor, rcutoff, parshell                            &
                                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                                &, jpseudo, pseudolocal, istart, costz, costz3         &
                                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                                &, psidetln, 1, nelorbh, nelorbjh                         &
                                &, niesd&
                                &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast&
                                &, elocm, logpsim&
                                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                        stm_image_first = signpsi*dexp(psiln(1) - spsiln - shiftlog)*psisn(1)*(-1)**(i - nel_read)
                        qpwf_image(n_extr_points) = stm_image_first
                        do j = 1, n_extr_points - 1
                            if (delta(1) .le. (da(1) - da(1)/n_extr_points*j)/2 .and. &
                                delta(2) .le. (da(2) - da(2)/n_extr_points*j)/2 .and. &
                                delta(3) .le. (da(3) - da(3)/n_extr_points*j)/2) then
                                qpwf_image(n_extr_points - j) = qpwf_image(n_extr_points - j) + stm_image_first
                            else
                                exit
                            end if
                        end do
                    else
                        iout = i - 1

                        call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                                &, ainv, kel, rcart, iesdr, vj, dupr&
                                &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr&
                                &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur&
                                &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar&
                                &, winvjbar, winvjbarsz, iflagnorm, cnorm, iflagerr&
                                &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne&
                                &, jastrowall_ee(1, iout, 0, 1), jasnew_ee, jastrowall_ei(1, iout, 1)&
                                &, jasnew_ei, n_body_on, niesd, nshell, 0&
                                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                        qpwf_image(n_extr_points) = qpwf_image(n_extr_points) + stm_image_first*ratior(1)*(-1)

                        do j = 1, n_extr_points - 1
                            if (delta(1) .le. (da(1) - da(1)/n_extr_points*j)/2 .and. &
                                delta(2) .le. (da(2) - da(2)/n_extr_points*j)/2 .and. &
                                delta(3) .le. (da(3) - da(3)/n_extr_points*j)/2) then
                                qpwf_image(n_extr_points - j) = qpwf_image(n_extr_points - j) + stm_image_first*ratior(1)*(-1)
                            else
                                exit
                            end if
                        end do
                    end if
                end if
            end do
        end if

    end subroutine compute_extrapolate

    subroutine compute_elec_delete(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog, mesh, vell)
        use allio
        use grid_module
        implicit none

        integer first_grid_point
        integer grid_point
        double precision vell(3)
        double precision rnew(3, *)
        double precision shiftlog, spsiln
        double precision delta(3)
        logical :: ele_first
        real*4 signpsi
        integer nel_read, nelup_read
        integer mesh(3)
        integer point(3)
        integer :: i, j

        delta(:) = 0.d0

        psisn = 0
        itest = 2
        iflagnorm = 3
        !    psiln=0.d0
        singdet = .true.
        test_AAD = .false.
        pseudologic = .false.
        iesrandoml = .false.
        itestrfn = 2
        nshelljmax = max(nshellj, 1)
        neldomax = max(neldo, 1)
        indtmax = max(indt, 1)
        nelorbjmax = max(nelorbj, 1)
        if (contraction .ne. 0) then
            yesprojm = .true.
        else
            yesprojm = .false.
        end if

        ele_first = .true.

        do i = 1, grid_points
            qpwf_image(i) = 0.0
        end do

        if (nelup_read .gt. nelup) then

            do i = 1, nelup_read

                !        i=nelup_read

                !        point(:)=int((rnew(:,i)-grid_start(:)+vell(:))/vell(:))
                point(1) = int((rnew(1, i) - grid_start(1) + vell(1))/vell(1))
                point(2) = int((rnew(2, i) - grid_start(2) + vell(2))/vell(2))
                point(3) = int((rnew(3, i) - grid_start(3) + vell(3))/vell(3))

                if ((point(1) .gt. 0 .and. point(1) .le. mesh(1)) .and.&
                        &(point(2) .gt. 0 .and. point(2) .le. mesh(2)) .and.&
                        &(point(3) .gt. 0 .and. point(3) .le. mesh(3))) then

                    delta(1) = abs(rnew(1, i) - (grid_start(1) + point(1)*vell(1) - vell(1)/2.0))
                    delta(2) = abs(rnew(2, i) - (grid_start(2) + point(2)*vell(2) - vell(2)/2.0))
                    delta(3) = abs(rnew(3, i) - (grid_start(3) + point(3)*vell(3) - vell(3)/2.0))

                    if (delta(1) .le. da(1)/2 .and. delta(2) .le. da(2)/2 .and. delta(3) .le. da(3)/2) then

                        grid_point = ((point(3) - 1.0)*mesh(2) + (point(2) - 1.0))*mesh(1) + point(1)

                        if (ele_first) then

                            ele_first = .false.

                            first_grid_point = grid_point

                            rcart(:) = rnew(:, i)

                            do j = 1, i - 1
                                kel(:, j) = rnew(:, j)
                            end do
                            do j = i + 1, nel_read
                                kel(:, j - 1) = rnew(:, j)
                            end do

                            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                                    &, tabpip, kel, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz &
                                    &, cnorm, iflagerr, npsa, lmax &
                                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                                    &, psidetln, 1, nelorbh, nelorbjh                         &
                                    &, niesd&
                                    &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast&
                                    &, elocm, logpsim&
                                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                            qpwf_image(first_grid_point) &
                                = signpsi*dexp(psiln(1) - spsiln - shiftlog)*psisn(1)*(-1)**(i - nelup_read)

                        else
                            iout = i - 1

                            call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                                    &, ainv, kel, rcart, iesdr, vj, dupr&
                                    &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr&
                                    &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur&
                                    &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar&
                                    &, winvjbar, winvjbarsz, iflagnorm, cnorm, iflagerr&
                                    &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne&
                                    &, jastrowall_ee(1, iout, 0, 1), jasnew_ee, jastrowall_ei(1, iout, 1)&
                                    &, jasnew_ei, n_body_on, niesd, nshell, 0&
                                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                            qpwf_image(grid_point) = qpwf_image(grid_point) + qpwf_image(first_grid_point)*ratior(1)*(-1)

                        end if
                    end if
                end if
            end do
        else
            do i = nelup_read + 1, nel_read
                !        i=nel_read
                point(1) = int((rnew(1, i) - grid_start(1) + vell(1))/vell(1))
                point(2) = int((rnew(2, i) - grid_start(2) + vell(2))/vell(2))
                point(3) = int((rnew(3, i) - grid_start(3) + vell(3))/vell(3))

                if ((point(1) .gt. 0 .and. point(1) .le. mesh(1)) .and.&
                        &(point(2) .gt. 0 .and. point(2) .le. mesh(2)) .and.&
                        &(point(3) .gt. 0 .and. point(3) .le. mesh(3))) then

                    delta(1) = abs(rnew(1, i) - (grid_start(1) + point(1)*vell(1) - vell(1)/2.0))
                    delta(2) = abs(rnew(2, i) - (grid_start(2) + point(2)*vell(2) - vell(2)/2.0))
                    delta(3) = abs(rnew(3, i) - (grid_start(3) + point(3)*vell(3) - vell(3)/2.0))

                    if (delta(1) .le. da(1)/2 .and. delta(2) .le. da(2)/2 .and. delta(3) .le. da(3)/2) then

                        grid_point = ((point(3) - 1)*mesh(2) + (point(2) - 1))*mesh(1) + point(1)

                        if (ele_first) then

                            ele_first = .false.
                            first_grid_point = grid_point
                            rcart(:) = rnew(:, i)

                            do j = 1, i - 1
                                kel(:, j) = rnew(:, j)
                            end do
                            do j = i + 1, nel_read
                                kel(:, j - 1) = rnew(:, j)
                            end do

                            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                                    &, tabpip, kel, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz &
                                    &, cnorm, iflagerr, npsa, lmax &
                                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                                    &, psidetln, 1, nelorbh, nelorbjh                         &
                                    &, niesd&
                                    &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn&
                                    &, yesfast, elocm, logpsim&
                                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                            qpwf_image(first_grid_point) = signpsi*dexp(psiln(1) - spsiln - shiftlog)*psisn(1)*(-1)**(i - nel_read)

                        else
                            iout = i - 1

                            call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                                    &, ainv, kel, rcart, iesdr, vj, dupr&
                                    &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr&
                                    &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur&
                                    &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar&
                                    &, winvjbar, winvjbarsz, iflagnorm, cnorm, iflagerr&
                                    &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne&
                                    &, jastrowall_ee(1, iout, 0, 1), jasnew_ee, jastrowall_ei(1, iout, 1)&
                                    &, jasnew_ei, n_body_on, niesd, nshell, 0&
                                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                            qpwf_image(grid_point) = qpwf_image(grid_point) + qpwf_image(first_grid_point)*ratior(1)*(-1)
                        end if
                    end if
                end if
            end do
        end if
    end subroutine compute_elec_delete

    subroutine compute_QPWF_k_space(rnew, nel_read, nelup_read, signpsi, spsiln, shiftlog)
        use allio
        use grid_module

        implicit none

        double precision stm_init
        double precision rnew(3, *)
        double precision shiftlog, spsiln

        real*4 signpsi

        integer nel_read, nelup_read

        integer :: i, j

        !    if(allocated(psiln)) psiln=0.d0
        !    if(allocated(psiln)) singdet=.true.
        itest = 2
        iflagnorm = 3
        test_AAD = .false.
        pseudologic = .false.
        iesrandoml = .false.
        itestrfn = 2
        nshelljmax = max(nshellj, 1)
        neldomax = max(neldo, 1)
        indtmax = max(indt, 1)
        nelorbjmax = max(nelorbj, 1)

        if (contraction .ne. 0) then
            yesprojm = .true.
        else
            yesprojm = .false.
        end if

        if (nelup_read .gt. nelup) then

            rcart(:) = rnew(:, 1)

            do j = 2, nel_read
                kel(:, j - 1) = rnew(:, j)
            end do

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kel, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            stm_init = dexp(psiln(1) - spsiln - shiftlog)*signpsi*psisn(1)*(-1)**(nelup_read - 1)

            do i = 1, grid_points
                if (active_points(i)) then
                    qpwf_image(i) &
                        = stm_init*cos(out_grid(1, i)*rnew(1, 1) + out_grid(2, i)*rnew(2, 1) + out_grid(3, i)*rnew(3, 1))
                    qpwf_image(i + grid_points) &
                        = stm_init*sin(out_grid(1, i)*rnew(1, 1) + out_grid(2, i)*rnew(2, 1) + out_grid(3, i)*rnew(3, 1))
                else
                    qpwf_image(i) = 0
                    qpwf_image(i + grid_points) = 0
                end if
            end do

            do j = 2, nelup_read
                iout = j - 1
                call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                        &, ainv, kel, rcart, iesdr, vj, dupr&
                        &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr&
                        &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur&
                        &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar&
                        &, winvjbar, winvjbarsz, iflagnorm, cnorm, iflagerr&
                        &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne&
                        &, jastrowall_ee(1, iout, 0, 1), jasnew_ee, jastrowall_ei(1, iout, 1)&
                        &, jasnew_ei, n_body_on, niesd, nshell, 0&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                do i = 1, grid_points
                    if (active_points(i)) then
                        qpwf_image(i) = qpwf_image(i) + stm_init*ratior(1)*(-1) &
                                        *cos(out_grid(1, i)*rnew(1, j) + out_grid(2, i)*rnew(2, j) + out_grid(3, i)*rnew(3, j))
                        qpwf_image(i + grid_points) &
                            = qpwf_image(i + grid_points) + stm_init*ratior(1)*(-1) &
                              *sin(out_grid(1, i)*rnew(1, j) + out_grid(2, i)*rnew(2, j) + out_grid(3, i)*rnew(3, j))
                    end if
                end do
            end do

        else

            rcart(:) = rnew(:, nelup_read + 1)

            do j = 1, nelup_read
                kel(:, j) = rnew(:, j)
            end do
            do j = nelup_read + 2, nel_read
                kel(:, j - 1) = rnew(:, j)
            end do

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kel, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            stm_init = dexp(psiln(1) - spsiln - shiftlog)*signpsi*psisn(1)*(-1)**(nel_read - nelup_read - 1)
            do i = 1, grid_points
                if (active_points(i)) then
                    qpwf_image(i) &
                        = stm_init*cos(out_grid(1, i)*rnew(1, nelup_read + 1) &
                                       + out_grid(2, i)*rnew(2, nelup_read + 1) + out_grid(3, i)*rnew(3, nelup_read + 1))
                    qpwf_image(i + grid_points) &
                        = stm_init*sin(out_grid(1, i)*rnew(1, nelup_read + 1) &
                                       + out_grid(2, i)*rnew(2, nelup_read + 1) + out_grid(3, i)*rnew(3, nelup_read + 1))
                else
                    qpwf_image(i) = 0
                    qpwf_image(i + grid_points) = 0
                end if

            end do

            do j = nelup_read + 2, nel_read
                iout = j - 1
                call ratiovar(iout, nelorb, nelorbh, nelup, neldo, ratior            &
                        &, ainv, kel, rcart, iesdr, vj, dupr&
                        &, zetar, rion, psip, ioccup, ioccdo, ioptorb, nshellr, nshellr&
                        &, ainvs, winvs, r, rmu, nion, kion, ioccj, kionj, vjur&
                        &, nelorbj, nelorbjh, ioptorbj, nshelljr, winvsj, winvbar&
                        &, winvjbar, winvjbarsz, iflagnorm, cnorm, iflagerr&
                        &, ratiodet, costz, costz3, iessz, LBox, rmucos, rmusin, timewf, rcne&
                        &, jastrowall_ee(1, iout, 0, 1), jasnew_ee, jastrowall_ei(1, iout, 1)&
                        &, jasnew_ei, n_body_on, niesd, nshell, 0&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                do i = 1, grid_points
                    if (active_points(i)) then
                        qpwf_image(i) &
                            = qpwf_image(i) + stm_init*ratior(1)*(-1) &
                              *cos(out_grid(1, i)*rnew(1, j) + out_grid(2, i)*rnew(2, j) + out_grid(3, i)*rnew(3, j))
                        qpwf_image(i + grid_points) &
                            = qpwf_image(i + grid_points) + stm_init*ratior(1)*(-1) &
                              *sin(out_grid(1, i)*rnew(1, j) + out_grid(2, i)*rnew(2, j) + out_grid(3, i)*rnew(3, j))
                    end if
                end do
            end do

            !      do i=1,grid_points
            !        qpwf_image(i)=qpwf_image(i) !*4*pi !/(sqrt(2*pi))**3
            !      enddo

        end if
    end subroutine compute_QPWF_k_space

end module qpwf_module

