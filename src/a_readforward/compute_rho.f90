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

module Rho_corr_module

    logical :: ifrho_corr
    logical, private :: yesprojm
    integer, private :: nelorbjmax, indtmax, neldomax, nshelljmax
    real*8, private :: logpsim(2), elocm(2)

contains

    subroutine compute_rho_corr(rnew, nel_read, nelup_read, density_c, signpsi, spsiln, shiftlog)
        use allio
        use grid_module
        implicit none

        integer nel_read, nelup_read

        double precision rnew(3, *)
        double precision density_c(*)
        double precision shiftlog, spsiln

        real signpsi

        integer i, j
        !    defines once for all the parameters for upscratch
        iflagnorm = 3
        !     psisn=0
        singdet = .true.
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
        itest = 2

        if (contraction .ne. 0) then
            yesprojm = .true.
        else
            yesprojm = .false.
        end if

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
                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                &, nintpseudo, prefactor, rcutoff, parshell                            &
                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                &, jpseudo, pseudolocal, istart, costz, costz3         &
                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                &, psidetln, 1, nelorbh, nelorbjh                         &
                &, niesd&
                &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscale&
                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

        density_c(1) = dexp(2.d0*(psiln(1) - spsiln - shiftlog))
        write (20, *) density_c(1), signpsi*dexp(psiln(1) - spsiln - shiftlog)*psisn(1)

        do i = 2, grid_points
            if (nelup_read .lt. nelup) then
                rcart(:) = out_grid(:, i)
                iout = nelup
            else
                rcart(:) = out_grid(:, i)
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

            density_c(i) = density_c(1)*ratior(1)**2.d0
        end do

    end subroutine compute_rho_corr

end module Rho_corr_module

