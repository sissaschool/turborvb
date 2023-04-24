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

module convertmod
    use allio
    implicit none
    real(dp), allocatable, dimension(:, :) :: overs
    integer nummol, nmolmatdo
    integer, allocatable, dimension(:) :: wheremol
    real(dp), private, allocatable, dimension(:, :) :: &
            &buffer, buffer_c, oversl, oversj, oversjsz, molecorb, molecorb_c
    !     logical, allocatable, dimension(:) :: occorb
    integer, private :: i, j, k, ind, nbuf, nbufu, nleft               &
            &, nelorbc_in, nelorbjc_in, niesd_in, iesdrr_in, imax, molecsw, molecsym  &
            &, nelcolc_in, nshellmol, iesupr_cmol, occ_cmol, iesup_cmol, nmol_in           &
            &, iesupindmol, nmoltot, i_max, iyt, ii, jj, kk, nnozeroc_in, nnozerojc_in    &
            &, iesswrc_in, indnn_in, inds, nelorb2, iesfreerc_in, iesupc_in&
            &, contraction_in, contractionj_in, indunp, nelorb_diag&
            &, nelorbl2, nmoldiff, ntotpsi, indref, indleft, imin, nelorb_diagu, nelorbpf
    integer*8, private :: indr, mesh
    real(dp), private :: x(3), volmesh, tracem, overmax, dnrm2, maxpsi, rannum
    real*8 rion_ref(3)
    real(dp), public, allocatable, dimension(:) :: eigmol
    real(dp), private, allocatable, dimension(:) :: eigmolsz&
            &, psipu, contrnorm, agp_kp, agpo_kp
    real(dp), private, allocatable, dimension(:, :) :: mat, oversnl, psi_out&
            &, psi_bar, over_sav, psi_sav
    complex*16 volmeshc
    logical occopt
    logical, private :: iessz_in, yest
    logical, private, allocatable, dimension(:) :: occorb, sjbradet_in
    integer, private, allocatable, dimension(:) :: ind_unpair         &
            &, nozeroc_in, jbradetc_in, jbrajc_in, nozerojc_in, jbradetnc_in, jbrajnc_in&
            &, multpointer, address

    integer, private, allocatable, dimension(:, :) :: mupointer
    integer, private, allocatable, dimension(:) :: imin_kp
    integer, private :: indorb, indpar
contains

    subroutine shift_originref
        implicit none
        real*8 mind(3)
        real*8, allocatable :: rion_sav(:, :)
#ifdef PARALLEL
        include 'mpif.h'
#endif
        !     the reference is the average ion position
        if (.not. iespbc) then
            !     the reference is the average ion position
            do j = 1, 3
                rion_ref(j) = sum(rion(j, :))/nion
            end do
            rion_ref(1) = rion_ref(1) - (nx - 1)/2.d0*ax
            rion_ref(2) = rion_ref(2) - (ny - 1)/2.d0*ay
            rion_ref(3) = rion_ref(3) - (nz - 1)/2.d0*az
        else
            rion_ref = 0.d0
        end if

        if (shiftx) rion_ref(:) = rion_ref(:) + ax/2.d0*at(:, 1)
        if (shifty) rion_ref(:) = rion_ref(:) + ay/2.d0*at(:, 2)
        if (shiftz) rion_ref(:) = rion_ref(:) + az/2.d0*at(:, 3)

        if (shift_origin) then
            allocate (rion_sav(3, nion))
            rion_sav = rion
            if (iespbc) then
                call CartesianToCrystal(rion_sav, nion)
                call CartesianToCrystal(rion_ref, 1)
            end if
            call findrionfref(nion, ax, rion_sav, rion_ref, mind(1))
            call findrionfref(nion, ay, rion_sav(2, 1), rion_ref(2), mind(2))
            call findrionfref(nion, az, rion_sav(3, 1), rion_ref(3), mind(3))
            if (iespbc) rion_ref(:) = rion_ref(1)*at(:, 1) + rion_ref(2)*at(:, 2) + rion_ref(3)*at(:, 3)
            deallocate (rion_sav)
        end if

#ifdef PARALLEL
        !  Just to be consistent with all processors.
        call mpi_bcast(rion_ref, 3, MPI_DOUBLE_PRECISION, 0, commopt_mpi, ierr)
#endif

        rion_ref(:) = rion_ref(:) + (nx - 1)/2.d0*ax*at(:, 1)
        rion_ref(:) = rion_ref(:) + (ny - 1)/2.d0*ay*at(:, 2)
        rion_ref(:) = rion_ref(:) + (nz - 1)/2.d0*az*at(:, 3)

    end subroutine shift_originref

    subroutine convertmol_fast
        use allio, only: norm_metric
        implicit none
        real*8 ddot, dnrm2, overmax, cost, costn, r0, psiln, rc(3), jastrow_ei, &
                &overlapsquare, overo, timechange, sumdet
        real*8, external :: cclock
        real ran
        integer iseed, dimvect, icount, iicount, indfirst, indlast, maxdimmol&
                &, mine, dimo, inddetc, nbufp, ind_c, shiftnelorb_diag, ind_t&
                &, ind, nmol_add
        logical symmagp_in
        integer, dimension(:), allocatable :: indexo, indexn
        real*8, dimension(:), allocatable :: sortvect, eigmat, umat, mat, mats
        real*8, dimension(:, :), allocatable :: inv_sav
        real*8, dimension(:, :), allocatable :: oeff
        real*8, external :: tracemat, tracemat2
#ifdef PARALLEL
        include 'mpif.h'
#endif

        !---------------------------------------------------------------------------------
        !       This (complicated) subroutine appends or replaces molecular orbitals at the
        !       end of the fort.10, the wf in TurboRVB notations.
        !       main input parameter: molopt
        !       In this definition we pay also attention to the coefficients of the
        !       contracted orbitals that do not correspond to molecular orbitals
        !       in other words to the atomic contracted orbitals.
        !       molopt=+/-1   do not care about the coefficients of the contracted
        !       molopt=+/-2   evaluate also the coefficients of the AGP with DMRG
        !       molopt=+/-3   evaluate also the coefficients of the AGP and J with DMRG
        !       molopt >= 2 do not optimize these coefficients
        !       molopt < 2  optimize always the coefficients of the contracted orbitals
        !--------------------------------------------------------------------------------

        !#ifdef __CASO
        !    nprocu=nprocopt
        !#else
        !    nprocu=1
        !#endif
        timechange = cclock()

        allocate (nozeroc_in(nnozero_c), jbradetc_in(nnozero_c + npar_eagp)      &
                &, jbradetnc_in(3*nnozero_c), sjbradet_in(nnozero_c))

        if (nnozeroj_c .gt. 0) then
            allocate (nozerojc_in(nnozeroj_c), jbrajc_in(nnozeroj_c) &
                      , jbrajnc_in(3*nnozeroj_c))
        end if

        iessz_in = iessz
        niesd_in = niesd
        iesdrr_in = iesdrr
        contraction_in = contraction

        symmagp_in = symmagp

        iesupc_in = iesup_c

        nelcolc_in = nelcol_c
        nelorbc_in = nelorb_c
        nelorbjc_in = nelorbj_c

        if (contraction .gt. 0) then
            nnozeroc_in = nnozero_c
            nozeroc_in = nozero_c
            jbradetc_in = jbradet
            jbradetnc_in = jbradetn
        else
            nnozeroc_in = nnozero_c
            nozeroc_in = nozero
            jbradetc_in = jbradet
            jbradetnc_in = jbradetn
        end if
        sjbradet_in = sjbradet

        nelcolc_in = nelcol_c
        nelorbc_in = nelorb_c

        mesh = nx
        mesh = mesh*ny
        mesh = mesh*nz

        volmesh = ax*ay*az*unit_volume
        volmeshc = dcmplx(volmesh)

        call shift_originref

        x = 0.d0

        !       first part calculation overlap determinant
        ! just the minimum allocation

        nmoldiff = nmol + ndiff
        nmoltot = nmoldiff
        if ((.not. symmagp .or. ipc .eq. 2) .or. ipf .eq. 2) nmoltot = nmoltot + nmol

        nelorb_at = nelorb_c - molecular
        !     count the number of matrix elements not concerning the atomic basis
        molecsw = 0
        molecsym = 0
        if (molecular .gt. 0) then
            allocate (occorb(iesswr))
            occorb = .false.
            do i = 1, nnozeroc_in
                iy = (nozero_c(i) - 1)/nelorb_c + 1
                ix = nozero_c(i) - (iy - 1)*nelorb_c
                !      whatever matrix element coupling a molecular orbital to be removed
                if (ix .gt. nelorb_at .or. (iy .gt. nelorb_at .and. iy .le. nelorb_c)) then
                    molecsw = molecsw + 1
                    if (jbradet(i) .ne. 0) then
                        if (.not. occorb(abs(jbradet(i)))) molecsym = molecsym + 1
                    end if
                end if
                if (jbradet(i) .ne. 0) then
                    occorb(abs(jbradet(i))) = .true.
                end if
            end do
            !     now count how many are left
            indleft = 0
            do i = 1, iesswr
                if (.not. occorb(i)) indleft = indleft + 1
            end do

            !     sroll jbradetn (the matrix element not optimized)
            indnn = 0
            do i = 1, indleft
                indnn = indnn + 1
                ii = jbradetn(indnn)
                ix = abs(jbradetn(indnn + 1))
                iy = abs(jbradetn(indnn + 2))
                if (ix .gt. nelorb_at .or. (iy .gt. nelorb_at .and. iy .le. nelorb_c))&
                        & molecsym = molecsym + 1
                indnn = indnn - 2*ii
            end do

            deallocate (occorb)

        end if

        nelorb_diag = nelorb_c
        nelorbpf = ipf*nelorbh
        nelorb_diagu = nelorb_diag
        if (npar_eagp .gt. 0) then
            nelorbpf = nelorbpf + ndiff
            nelorb_diagu = nelorb_diagu + ndiff
        end if

        if (symmagp .and. ipc .eq. 1 .or. ipf .eq. 2) then
            nelorb2 = nelorb_diagu*nelorb_diagu
        else
            nelorb2 = ipc*2*nelorb_diagu*nelorb_diagu
        end if

        nelorbl2 = nelorbpf*nelorb_diagu
        if (ipf .eq. 1 .and. (.not. symmagp .or. ipc .eq. 2)) nelorbl2 = 2*nelorbl2

        nbufu = nbufd
        if (mesh .lt. nbufd) nbufu = mesh
        if (ipc .eq. 1) then
            allocate (buffer(ipf*nelorbh, nbufu))
        else
            if (ipf .eq. 2) then
                allocate (buffer(4*nelorbh, nbufu))
            else
                allocate (buffer(2*nelorbh, 2*nbufu))
            end if
        end if
        nbufp = nbufu + 1
        if (symmagp .and. ipc .eq. 1) then
            allocate (overs(nelorb_diagu, nelorb_diagu))
            allocate (oversl(nelorbpf, nelorb_diagu))
        else
            if (ipf .eq. 2) then
                allocate (overs(ipc*nelorb_diagu, nelorb_diagu))
                allocate (oversl(ipc*nelorbpf, nelorb_diagu))
            else
                allocate (overs(ipc*nelorb_diag, 2*nelorb_diag))
                allocate (oversl(ipc*nelorbh, 2*nelorb_diag))
            end if
        end if
        overs = 0.d0
        oversl = 0.d0
        allocate (eigmol(max(nelorb_diagu, nmoltot)))
        eigmol = 0.d0
        if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
            maxdimmol = max(nelorbl2, nelorb_diagu*nelorb_diagu, nmoltot*nelorbpf&
                    &, nelorb_diag*nbufu)
            maxdimmol = maxdimmol/(ipf*nelorbh) + 1
            allocate (molecorb(ipc*nelorbpf, maxdimmol))
        else
            maxdimmol = max(nelorbl2, 2*nelorb_diag*nelorb_diag, nmoltot*ipf*nelorbh&
                    &, 2*nelorb_diag*nbufu)
            maxdimmol = maxdimmol/(ipf*nelorbh) + 1
            ! complex MOs coefficients
            allocate (molecorb(ipc*ipf*nelorbh, maxdimmol))
        end if
        !ccc
        !    deallocate(molecorb)
        !    allocate(molecorb(nelorb_diag*ipc, nbufu))

        if (.not. allocated(contrnorm)) allocate (contrnorm(nelorb_diag))

        contrnorm = 1.d0
        buffer = 0.d0
        eigmol = 0.d0
        molecorb = 0.d0
        overs = 0.d0
        oversl = 0.d0

        ! compute overlap matrices and diagonalize
        ! atomic detmat in order to evaluate molecular orbitals
        ! (subroutine eval_molec_epsdgel)

        !CCC I didn't touch anything before!
        if (contraction .eq. 0 .and. .not. detc_proj) then
            if (allocated(mu_c)) deallocate (mu_c, detmat_c)
            allocate (mu_c(ipc*nelorbh*ipf, nelorbh*ipf), detmat_c(ipf*ipc*nelorbh*(ipf*nelorbh + ndiff)))
            mu_c = 0.d0
            do i = 1, ipf*nelorbh
                mu_c(ipc*(i - 1) + 1, i) = 1.d0
            end do
            detmat_c = detmat
        end if

        if (epsdgm .ge. 0.d0) then
            !     The odd molecular orbitals are already conjugated

            overs = 0.d0
            oversl = 0.d0

            iflagnorm = 3
            indr = 0
            ind = 0
            nbuf = 0

            do k = 1, nz
                !         x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
                do j = 1, ny
                    !            x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                    do i = 1, nx
                        indr = indr + 1

                        if (indr - (indr/nprocopt)*nprocopt .eq. rankopt) then
                            ind = ind + 1
                            !                  x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)
                            x(:) = rion_ref(:) + (-(nz + 1)/2.d0 + k)*az*at(:, 3) + &
                                   (-(ny + 1)/2.d0 + j)*ay*at(:, 2) + &
                                   (-(nx + 1)/2.d0 + i)*ax*at(:, 1)

                            call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                    &, dupr, zetar, rion, psip, buffer(1, ind), nelorbh, nion, kion             &
                                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                    &, indpar_tab, indorb_tab, indshell_tab, .true.)
                            if (ipc .eq. 1 .and. ipf .eq. 2) then
                                buffer(nelorbh + 1:2*nelorbh, ind) = buffer(1:nelorbh, ind)
                            end if

                            if (ipc .eq. 2) then
                                if (ipf .eq. 2) then
                                    call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu &
                                            &, dupr, zetar, rion, psip, buffer(2*nelorbh + 1, ind), nelorbh, nion, kion &
                                            &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                            &, indpar_tab, indorb_tab, indshell_tab, .false.)
                                else
                                    call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu &
                                            &, dupr, zetar, rion, psip, buffer(1, ind + nbufu), nelorbh, nion, kion &
                                            &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                            &, indpar_tab, indorb_tab, indshell_tab, .false.)
                                end if
                            end if
                            !  REMOVED. It is overflow unstable for small vj and the convergence in the mesh is
                            !  slower due to implicit cusps in the orbitals added by the one body Jastrow.
                            !  With effective Hamiltonian it is important to include these terms to be
                            !  consistent with DFT and with the fact that we are interested in the
                            !  matrix elements of the orbitals that are indeed used in the determinantal
                            !  part, including the one body term that allows to satisfy the el-ion cusp.

                            if (add_onebody2det) then
                                psiln = -scale_one_body
                                do jj = 1, nion
                                    if (iespbc) then
                                        rc(:) = x(:) - rion(:, jj)
                                        call CartesianToCrystal(rc, 1)
                                        do kk = 1, 3
                                            rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                                        end do
                                        r0 = norm_metric(rc, metric)
                                    else
                                        rc(:) = (x(:) - rion(:, jj))*costz(jj)
                                        r0 = dsqrt(sum(rc(:)**2))
                                    end if
                                    psiln = psiln - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                                end do
                                buffer(1:ipf*ipc*nelorbh, ind) = buffer(1:ipf*ipc*nelorbh, ind)*dexp(psiln)
                                if (ipc .eq. 2 .and. ipf .eq. 1) &
                                        &buffer(1:2*nelorbh, nbufu + ind) = buffer(1:2*nelorbh, ind + nbufu)*dexp(psiln)
                            end if

                            !CCC qui si calcolano le matrici di overlap, credo che i conti con la base scontratta
                            !vadano divisi in due come nel compute_fast per il caso contratto, mentre per il caso
                            !scontratto dovrebbe essere deallocata la matrice di overlap e riallocata di dimensione
                            !doppia mettendo la parte up e down nei due quadranti diagonali, non lo faccio in attesa
                            !di conferma.
                            !Non ho capito perché questa parte è messa dentro i cicli!

                            if (mod(ind, nbufu) .eq. 0) then
                                if (contraction .eq. 0) then
                                    if (ipc .eq. 1) then
                                        call dgemm('N', 'T', nelorbh, nelorb_diag/ipf, nbufu, volmesh, buffer&
                                                &, ipf*nelorbh, buffer, ipf*nelorbh, 1.d0, oversl, nelorbpf)
                                    else
                                        call zgemm('N', 'C', nelorbh, nelorb_diag/ipf, nbufu, volmeshc, buffer&
                                                &, ipf*nelorbh, buffer, ipf*nelorbh, zone, oversl, nelorbpf)
                                        if (ipf .eq. 2) then
                                            call zgemm('N', 'C', nelorbh, nelorb_diag/2, nbufu, volmeshc&
                                                    &, buffer(2*nelorbh + 1, 1), 2*nelorbh, buffer(2*nelorbh + 1, 1), 2*nelorbh&
                                                    &, zone, oversl(2*nelorbh + 1, nelorb_diag/2 + 1), nelorbpf)
                                        end if

                                    end if
                                else
                                    if (ipc .eq. 1) then
                                        if (ipf .eq. 1) then
                                            call dgemm('T', 'N', nelorb_diag, nbufu, nelorbh, 1.d0, mu_c, nelorbh&
                                                    &, buffer, nelorbh, 0.d0, molecorb, nelorb_diag)
                                            call dgemm('N', 'T', nelorbh, nelorb_diag, nbufu, volmesh, buffer&
                                                    &, nelorbh, molecorb, nelorb_diag, 1.d0, oversl, nelorbh)
                                        else

                                            !             Up overlap
                                            call dgemm('T', 'N', nelorb_diag, nbufu, nelorbh, 1.d0, mu_c, 2*nelorbh&
                                                    &, buffer, 2*nelorbh, 0.d0, molecorb, nelorb_diag)

                                            call dgemm('N', 'T', nelorbh, nelorb_diag, nbufu, volmesh, buffer&
                                                    &, 2*nelorbh, molecorb, nelorb_diag, 1.d0, oversl, nelorbpf)
                                            !             Down  overlap
                                            call dgemm('T', 'N', nelorb_diag, nbufu, nelorbh, 1.d0, mu_c(nelorbh + 1, 1)&
                                                    &, 2*nelorbh, buffer(nelorbh + 1, 1), 2*nelorbh, 0.d0, molecorb, nelorb_diag)
                                            call dgemm('N', 'T', nelorbh, nelorb_diag, nbufu, volmesh, buffer(nelorbh + 1, 1)&
                                                    &, 2*nelorbh, molecorb, nelorb_diag, 1.d0, oversl(nelorbh + 1, 1), nelorbpf)
                                        end if
                                    else
                                        if (ipf .eq. 1) then
                                            call zgemm('T', 'N', nelorb_diag, nbufu, nelorbh, zone, mu_c, nelorbh&
                                                    &, buffer, nelorbh, zzero, molecorb, nelorb_diag)
                                            call zgemm('N', 'C', nelorbh, nelorb_diag, nbufu, volmeshc, buffer&
                                                    &, nelorbh, molecorb, nelorb_diag, zone, oversl, nelorbh)
                                        else
                                            ! Up overlap
                                            call zgemm('T', 'N', nelorb_diag, nbufu, nelorbh, zone, mu_c, 2*nelorbh&
                                                    &, buffer, 2*nelorbh, zzero, molecorb, nelorb_diag)

                                            call zgemm('N', 'C', nelorbh, nelorb_diag, nbufu, volmeshc, buffer&
                                                    &, 2*nelorbh, molecorb, nelorb_diag, zone, oversl, nelorbpf)
                                            ! Down overlap
                                            call zgemm('T', 'N', nelorb_diag, nbufu, nelorbh, zone, mu_c(2*nelorbh + 1, 1)&
                                                        &, 2*nelorbh, buffer(2*nelorbh + 1, 1), 2*nelorbh, zzero, molecorb&
                                                        &, nelorb_diag)

                                            call zgemm('N', 'C', nelorbh, nelorb_diag, nbufu, volmeshc, buffer(2*nelorbh + 1, 1)&
                                                    &, 2*nelorbh, molecorb, nelorb_diag, zone, oversl(2*nelorbh + 1, 1), nelorbpf)

                                        end if
                                    end if
                                end if
                                if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then
                                    if (contraction .eq. 0) then
                                        if (ipc .eq. 1) then
                                            call dgemm('N', 'T', nelorbh, nelorb_diag, nbufu, volmesh, buffer&
                                                    &, nelorbh, buffer, nelorbh, 1.d0, oversl(1, nelorb_diag + 1), nelorbh)
                                        else
                                            call zgemm('N', 'C', nelorbh, nelorb_diag, nbufu, volmeshc, buffer(1, nbufp)&
                                                   &, nelorbh, buffer(1, nbufp), nelorbh, zone, oversl(1, nelorb_diag + 1)&
                                                   &, nelorbh)
                                        end if
                                    else
                                        if (ipc .eq. 1) then
                                            call dgemm('T', 'N', nelorb_diag, nbufu, nelorbh, 1.d0, mu_c, nelorbh&
                                                    &, buffer, nelorbh, 0.d0, molecorb, nelorb_diag)
                                            call dgemm('N', 'T', nelorbh, nelorb_diag, nbufu, volmesh, buffer&
                                                    &, nelorbh, molecorb, nelorb_diag, 1.d0, oversl(1, nelorb_diag + 1), nelorbh)
                                        else
                                            call zgemm('T', 'N', nelorb_diag, nbufu, nelorbh, zone, mu_c, nelorbh&
                                                    &, buffer(1, nbufp), nelorbh, zzero, molecorb, nelorb_diag)
                                            call zgemm('N', 'C', nelorbh, nelorb_diag, nbufu, volmeshc, buffer(1, nbufp)&
                                                    &, nelorbh, molecorb, nelorb_diag, zone, oversl(1, nelorb_diag + 1), nelorbh)
                                        end if
                                    end if
                                end if

                                nbuf = nbuf + 1
                                ind = 0

                            end if
                        end if
                    end do
                end do
            end do

            if (mod(ind, nbufu) .ne. 0) then
                nleft = ind
                if (contraction .eq. 0) then
                    if (ipc .eq. 1) then
                        call dgemm('N', 'T', nelorbh, nelorb_diag/ipf, nleft, volmesh, buffer&
                                &, ipf*nelorbh, buffer, ipf*nelorbh, 1.d0, oversl, nelorbpf)
                    else
                        call zgemm('N', 'C', nelorbh, nelorb_diag/ipf, nleft, volmeshc, buffer&
                                &, ipf*nelorbh, buffer, ipf*nelorbh, zone, oversl, nelorbpf)
                        if (ipf .eq. 2) then
                            call zgemm('N', 'C', nelorbh, nelorb_diag/2, nleft, volmeshc&
                                    &, buffer(2*nelorbh + 1, 1), 2*nelorbh, buffer(2*nelorbh + 1, 1), 2*nelorbh&
                                    &, zone, oversl(2*nelorbh + 1, nelorb_diag/2 + 1), nelorbpf)
                        end if

                    end if
                else
                    if (ipc .eq. 1) then
                        if (ipf .eq. 1) then
                            call dgemm('T', 'N', nelorb_diag, nleft, nelorbh, 1.d0, mu_c, nelorbh&
                                    &, buffer, nelorbh, 0.d0, molecorb, nelorb_diag)
                            call dgemm('N', 'T', nelorbh, nelorb_diag, nleft, volmesh, buffer&
                                    &, nelorbh, molecorb, nelorb_diag, 1.d0, oversl, nelorbh)
                        else
                            ! Overlap up
                            call dgemm('T', 'N', nelorb_diag, nleft, nelorbh, 1.d0, mu_c, 2*nelorbh&
                                    &, buffer, 2*nelorbh, 0.d0, molecorb, nelorb_diag)
                            call dgemm('N', 'T', nelorbh, nelorb_diag, nleft, volmesh, buffer&
                                    &, 2*nelorbh, molecorb, nelorb_diag, 1.d0, oversl, nelorbpf)
                            ! Overlap down
                            call dgemm('T', 'N', nelorb_diag, nleft, nelorbh, 1.d0, mu_c(nelorbh + 1, 1), 2*nelorbh&
                                    &, buffer(nelorbh + 1, 1), 2*nelorbh, 0.d0, molecorb, nelorb_diag)
                            call dgemm('N', 'T', nelorbh, nelorb_diag, nleft, volmesh, buffer(nelorbh + 1, 1)&
                                    &, 2*nelorbh, molecorb, nelorb_diag, 1.d0, oversl(nelorbh + 1, 1), nelorbpf)
                        end if
                    else
                        if (ipf .eq. 1) then
                            call zgemm('T', 'N', nelorb_diag, nleft, nelorbh, zone, mu_c, nelorbh&
                                    &, buffer, nelorbh, zzero, molecorb, nelorb_diag)
                            call zgemm('N', 'C', nelorbh, nelorb_diag, nleft, volmeshc, buffer&
                                    &, nelorbh, molecorb, nelorb_diag, zone, oversl, nelorbh)
                        else
                            !  Overlap up
                            call zgemm('T', 'N', nelorb_diag, nleft, nelorbh, zone, mu_c, 2*nelorbh&
                                    &, buffer, 2*nelorbh, zzero, molecorb, nelorb_diag)
                            call zgemm('N', 'C', nelorbh, nelorb_diag, nleft, volmeshc, buffer&
                                    &, 2*nelorbh, molecorb, nelorb_diag, zone, oversl, nelorbpf)
                            !  Overlap down
                            call zgemm('T', 'N', nelorb_diag, nleft, nelorbh, zone, mu_c(2*nelorbh + 1, 1), 2*nelorbh&
                                    &, buffer(2*nelorbh + 1, 1), 2*nelorbh, zzero, molecorb, nelorb_diag)
                            call zgemm('N', 'C', nelorbh, nelorb_diag, nleft, volmeshc, buffer(2*nelorbh + 1, 1)&
                                    &, 2*nelorbh, molecorb, nelorb_diag, zone, oversl(2*nelorbh + 1, 1), nelorbpf)
                        end if
                    end if
                end if
                if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then
                    if (contraction .eq. 0) then
                        if (ipc .eq. 1) then
                            call dgemm('N', 'T', nelorbh, nelorb_diag, nleft, volmesh, buffer&
                                    &, nelorbh, buffer, nelorbh, 1.d0, oversl(1, nelorb_diag + 1), nelorbh)
                        else
                            call zgemm('N', 'C', nelorbh, nelorb_diag, nleft, volmeshc, buffer(1, nbufp)&
                                    &, nelorbh, buffer(1, nbufp), nelorbh, zone, oversl(1, nelorb_diag + 1), nelorbh)
                        end if
                    else
                        if (ipc .eq. 1) then
                            call dgemm('T', 'N', nelorb_diag, nleft, nelorbh, 1.d0, mu_c, nelorbh&
                                    &, buffer, nelorbh, 0.d0, molecorb, nelorb_diag)
                            call dgemm('N', 'T', nelorbh, nelorb_diag, nleft, volmesh, buffer&
                                    &, nelorbh, molecorb, nelorb_diag, 1.d0, oversl(1, nelorb_diag + 1), nelorbh)
                        else
                            call zgemm('T', 'N', nelorb_diag, nleft, nelorbh, zone, mu_c, nelorbh&
                                    &, buffer(1, nbufp), nelorbh, zzero, molecorb, nelorb_diag)
                            call zgemm('N', 'C', nelorbh, nelorb_diag, nleft, volmeshc, buffer(1, nbufp)&
                                    &, nelorbh, molecorb, nelorb_diag, zone, oversl(1, nelorb_diag + 1), nelorbh)
                        end if
                    end if
                end if
                nbuf = nbuf + 1
            end if
#ifdef PARALLEL
            call reduce_base_real(size(oversl), oversl, commopt_mpi, -1)
#endif
            if (contraction .eq. 0 .and. ipf .eq. 2 .and. ipc .eq. 1) then
                oversl(nelorbh + 1:2*nelorbh, nelorb_diag/2 + 1:nelorb_diag) = &
                        &oversl(1:nelorbh, 1:nelorb_diag/2)
            end if

            if (ipc .eq. 2) then
                if (ipf .eq. 1) then
                    call conjmat(nelorbh, 2*nelorb_diag, oversl, nelorbh) ! Both conjugated
                else
                    call conjmat(nelorbpf, nelorb_diagu, oversl, nelorbpf) ! Both conjugated
                end if
            end if
            !       Computing the small matrix overs

            !Qua viene calcolato overs, ma se oversl viene calcolato prima non devo
            if (contraction .eq. 0) then
                overs = oversl
            else
                if (ipc .eq. 1) then
                    call dgemm_my('T', 'N', nelorb_diag, nelorb_diag, ipf*nelorbh, 1.d0, mu_c, ipf*nelorbh&
                            &, oversl, nelorbpf, 0.d0, overs, nelorb_diagu, nprocu, rankopt, commopt_mpi)

                else
                    call zgemm_my('C', 'N', nelorb_diag, nelorb_diag, ipf*nelorbh, zone, mu_c, ipf*nelorbh&
                            &, oversl, nelorbpf, zzero, overs, nelorb_diagu, nprocu, rankopt, commopt_mpi)
                end if
                if ((.not. symmagp .or. ipc .eq. 2) .and. ipf .eq. 1) then
                    if (ipc .eq. 1) then
                        call dgemm_my('T', 'N', nelorb_diag, nelorb_diag, nelorbh, 1.d0, mu_c, nelorbh&
                                &, oversl(1, nelorb_diag + 1), nelorbh, 0.d0, overs(1, nelorb_diag + 1), nelorb_diag&
                                &, nprocu, rankopt, commopt_mpi)
                    else
                        call zgemm_my('C', 'N', nelorb_diag, nelorb_diag, nelorbh, zone, mu_c, nelorbh&
                                &, oversl(1, nelorb_diag + 1), nelorbh, zzero, overs(1, nelorb_diag + 1), nelorb_diag&
                                &, nprocu, rankopt, commopt_mpi)
                    end if
                end if
            end if
            !CCC controllare se nelorb_at non va molt. per ipf (e' ridefinito sopra, discuterne con Sandro)

            sumdet = 0.d0
            do i = 1, nelorb_at
                do j = 1, ipc*nelorb_at
                    sumdet = sumdet + abs(detmat_c(ipc*nelorb_c*(i - 1) + j))
                end do
            end do

            ! The part below (until 'endif big if ') is used only when optimization
            ! with molecular orbitals is on. Thus it does  not need to be touched
            ! for converfort10mol use.
            ! This part translates an agp or a  pfaffian written with molecular
            ! orbitals in an agp or a pfaffian  written in the corresponding atomic
            ! basis. A much simpler and stable algorithm could be written if we
            ! would not have the problem that the molecular orbitals in Turbo
            ! are always written in the uncontracted basis in input (fort.10), so that
            ! one has to find the appropriate coefficients in the contracted basis.
            ! If the overlap is not 1, the wf is not consistent but the best
            ! approximation with the given contracted orbitals is done.

            if (detc_proj .and. sumdet .eq. 0.d0) then
                !  Translates molecular orbitals in atomic orbitals.
                dimo = nmoltot - ndiff
                allocate (umat(ipc*nelorb_at*nelorb_at), eigmat(nelorb_at)&
                        &, mat(ipc*dimo*max(nelorb_at, dimo)))
                if (.not. symmagp .or. ipc .eq. 2 .or. ipf .eq. 2) allocate (mats(ipc*nelorb_at*nelorb_at))
                if (ipc .eq. 2) then
                    allocate (inv_sav(2*nelorb_at, 2*nelorb_at))
                else
                    allocate (inv_sav(nelorb_at, nelorb_at))
                end if
                !         Compute both inverse at the beginning, in order to use umat later.
                do i = 1, ipc*nelorb_at
                    do j = 1, nelorb_at
                        inv_sav(i, j) = overs(i, j)
                    end do
                end do
                call invsymeps(ipc, nelorb_at, inv_sav&
                        &, nelorb_at, info, epsdgm, mine, umat, eigmat, nprocu, rankopt, commopt_mpi)
                if (ipf .eq. 2) then
                    shiftnelorb_diag = 0
                else
                    shiftnelorb_diag = nelorb_diag
                end if
                if (ipc .eq. 2) then

                    do i = 1, ipc*nelorb_at
                        do j = 1, nelorb_at
                            inv_sav(i, j + nelorb_at) = overs(i, j + shiftnelorb_diag)
                        end do
                    end do
                    call invsymeps(ipc, nelorb_at, inv_sav(1, nelorb_at + 1)&
                            &, nelorb_at, info, epsdgm, mine, umat, eigmat, nprocu, rankopt, commopt_mpi)
                end if
                inddetc = ipc*(nelorb_diag*nelorb_at + nelorb_at) + 1
                if (ipc .eq. 1) then
                    call dgemm_my('N', 'N', dimo, dimo, dimo, 1.d0, overs(nelorb_at + 1, nelorb_at + 1), nelorb_diag&
                            &, detmat_c(inddetc), nelorb_diag, 0.d0, molecorb, dimo, nprocu, rankopt, commopt_mpi)
                else
                    call conjmat(dimo, dimo, detmat_c(inddetc), nelorb_diag)
                    !  S_L^T \lambda^*
                    call zgemm_my('T', 'N', dimo, dimo, dimo, zone, overs(2*nelorb_at + 1, nelorb_at + 1), nelorb_diag&
                            &, detmat_c(inddetc), nelorb_diag, zzero, molecorb, dimo, nprocu, rankopt, commopt_mpi)
                    call conjmat(dimo, dimo, detmat_c(inddetc), nelorb_diag)
                end if

                if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                    overo = tracemat(dimo, molecorb, dimo)
                else
                    if (ipc .eq. 1) then
                        call dgemm_my('N', 'T', dimo, dimo, dimo, 1.d0, overs(nelorb_at + 1, shiftnelorb_diag + nelorb_at + 1)&
                                                &, nelorb_diag, detmat_c(inddetc), nelorb_diag, 0.d0, mat, dimo, nprocu, rankopt&
                                                &, commopt_mpi)
                    else
                        ! S_R \lambda^T
                        call zgemm_my('N', 'T', dimo, dimo, dimo, zone, overs(2*nelorb_at + 1, shiftnelorb_diag + nelorb_at + 1)&
                                                &, nelorb_diag, detmat_c(inddetc), nelorb_diag, zzero, mat, dimo, nprocu, rankopt&
                                                &, commopt_mpi)
                    end if
                    overo = tracemat2(dimo, dimo, molecorb, dimo, mat, dimo)
                end if !
                if (ipc .eq. 1) then

                    call dgemm_my('N', 'N', dimo, nelorb_at, nelorb_at, 1.d0, overs(nelorb_at + 1, 1)&
                            &, nelorb_diag, inv_sav, nelorb_at, 0.d0, molecorb, dimo, nprocu, rankopt, commopt_mpi)
                    !
                    call dgemm_my('N', 'N', dimo, nelorb_at, dimo, 1.d0&
                            &, detmat_c(inddetc), nelorb_diag, overs(nelorb_at + 1, 1), nelorb_diag, 0.d0, mat, dimo&
                            &, nprocu, rankopt, commopt_mpi)
                    call dgemm_my('T', 'N', nelorb_at, nelorb_at, dimo, 1.d0, molecorb&
                            &, dimo, mat, dimo, 0.d0, umat, nelorb_at, nprocu, rankopt, commopt_mpi)

                else

                    ! \bar S_L (S^\prime_L)^-1
                    call zgemm_my('N', 'N', dimo, nelorb_at, nelorb_at, zone, overs(2*nelorb_at + 1, 1)&
                            &, nelorb_diag, inv_sav, nelorb_at, zzero, molecorb, dimo, nprocu, rankopt, commopt_mpi)
                    !
                    !   lambda \bar S_R^*
                    call conjmat(dimo, nelorb_at, overs(2*nelorb_at + 1, shiftnelorb_diag + 1), nelorb_diag)
                    call zgemm_my('N', 'N', dimo, nelorb_at, dimo, zone, detmat_c(inddetc), nelorb_diag&
                            &, overs(2*nelorb_at + 1, shiftnelorb_diag + 1), nelorb_diag, zzero, mat, dimo&
                            &, nprocu, rankopt, commopt_mpi)
                    call conjmat(dimo, nelorb_at, overs(2*nelorb_at + 1, shiftnelorb_diag + 1), nelorb_diag)
                    !  umat=A = molecorb^dag mat= (S^\prime_L)^-1 \bar S_L^dag lambda \bar S_R^*
                    call zgemm_my('C', 'N', nelorb_at, nelorb_at, dimo, zone, molecorb&
                            &, dimo, mat, dimo, zzero, umat, nelorb_at, nprocu, rankopt, commopt_mpi)
                end if
                !

                !
                if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                    !
                    overlapsquare = tracemat(nelorb_at, umat, nelorb_at)/overo

                else
                    !

                    if (ipc .eq. 1) then
                        !   call dgemm_my('T','T',nelorb_at,dimo,dimo,1.d0,overs(nelorb_at+1,1),nelorb_diag&
                        !                 &,detmat_c(inddetc),nelorb_diag,0.d0,mat,dimo,nprocu,rankopt,commopt_mpi)
                        !          call dgemm_my('N','N',nelorb_at,dimo,nelorb_at,1.d0,inv_sav,nelorb_at&
                        !     &,mat,nelorb_at,0.d0,molecorb,nelorb_at,nprocu,rankopt,commopt_mpi)
                        !          call dgemm_my('N','N',nelorb_at,nelorb_at,dimo,1.d0,molecorb&
                        !      &,nelorb_at,overs(nelorb_at+1,1),nelorb_diag,0.d0,mats,nelorb_at,nprocu,rankopt,commopt_mpi)
                        !       In case the s R and L are the same the simplification below holds.

                        call dgemm_my('T', 'N', dimo, nelorb_at, dimo, 1.d0&
                                    &, detmat_c(inddetc), nelorb_diag, overs(nelorb_at + 1, 1)&
                                    &, nelorb_diag, 0.d0, mat, dimo, nprocu, rankopt, commopt_mpi)
                        call dgemm_my('T', 'N', nelorb_at, nelorb_at, dimo, 1.d0, molecorb&
                                &, dimo, mat, dimo, 0.d0, mats, nelorb_at, nprocu, rankopt, commopt_mpi)

                    else
                        !! S_R \lambda^T
                        !    call zgemm_my('N','T',dimo,dimo,dimo,zone,overs(2*nelorb_at+1,nelorb_diag+nelorb_at+1),nelorb_diag&
                        !             &,detmat_c(inddetc),nelorb_diag,zzero,mat,dimo,nprocu,rankopt,commopt_mpi)
                        !            compute inverse right. Here we cannot assume R and L are the same.
                        !   Definition of mat=\bar S_R^T \lambda^dag
                        call zgemm_my('T', 'C', nelorb_at, dimo, dimo, zone, overs(2*nelorb_at + 1, shiftnelorb_diag + 1)&
                                    &, nelorb_diag, detmat_c(inddetc), nelorb_diag, zzero, mat, nelorb_at, nprocu&
                                    &, rankopt, commopt_mpi)
                        !    B= mats=S_R^*-1 \bar S_R^T \lambda^\dag  \bar S_L
                        call conjmat(nelorb_at, nelorb_at, inv_sav(1, nelorb_at + 1), nelorb_at)
                        call zgemm_my('N', 'N', nelorb_at, dimo, nelorb_at, zone, inv_sav(1, nelorb_at + 1), nelorb_at&
                                &, mat, nelorb_at, zzero, molecorb, nelorb_at, nprocu, rankopt, commopt_mpi)
                        call conjmat(nelorb_at, nelorb_at, inv_sav(1, nelorb_at + 1), nelorb_at)

                        call zgemm_my('N', 'N', nelorb_at, nelorb_at, dimo, zone, molecorb&
                            &, nelorb_at, overs(2*nelorb_at + 1, 1), nelorb_diag, zzero, mats, nelorb_at, nprocu&
                            &, rankopt, commopt_mpi)

                    end if

                    overlapsquare = tracemat2(nelorb_at, nelorb_at, umat, nelorb_at, mats, nelorb_at)/overo

                end if
                if (kaverage) then
#ifdef PARALLEL
                    allocate (agp_kp(nk))
                    if (rankrep .eq. 0) call mpi_gather(overlapsquare, 1, MPI_DOUBLE_PRECISION, agp_kp,&
                        &1, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
                    if (rank .eq. 0) then
                        do ii = 1, nk
                            write (6, *) ' K#/Overlap square  ', ii, agp_kp(ii)
                            if (agp_kp(ii) .lt. 0.99999) write (6, *) ' ERROR Overlap K=', ii
                        end do
                    end if
                    deallocate (agp_kp)
#endif
                else
                    if (rank .eq. 0) write (6, *) ' Overlap square atomic vs molec =', overlapsquare
                end if

                if (ipc .eq. 1) then
                    call dgemm_my('N', 'N', nelorb_at, nelorb_at, nelorb_at, 1.d0, umat&
                            &, nelorb_at, inv_sav, nelorb_at, 0.d0, detmat_proj, nelorb_c, nprocu, rankopt, commopt_mpi)
                else
                    call zgemm_my('N', 'T', nelorb_at, nelorb_at, nelorb_at, zone, umat&
                        &, nelorb_at, inv_sav(1, nelorb_at + 1), nelorb_at, zzero, detmat_proj, nelorb_c, nprocu&
                        &, rankopt, commopt_mpi)
                end if

                !            if(rank.eq.0) then
                !             write(6,*) ' Output detmat_proj after convertmol_c '
                !             cost=0.d0
                !              do i=1,nelorb_at
                !               do j=i,nelorb_at
                !                write(6,*) i,j,detmat_proj(nelorb_c*(j-1)+i),detmat_proj(nelorb_c*(i-1)+j)
                !                write(6,*) i,j,detmat_proj(2*nelorb_c*(j-1)+2*i-1:2*nelorb_c*(j-1)+2*i),&
                !                           detmat_proj(2*nelorb_c*(i-1)+2*j-1:2*nelorb_c*(i-1)+2*j)
                !               cost=cost+abs(detmat_proj(2*nelorb_c*(j-1)+2*i-1)-detmat_proj(2*nelorb_c*(i-1)+2*j-1))
                !              enddo
                !             enddo
                !             write(6,*) ' Total asymmetry =',cost
                !             endif

                !
                if (ndiff .gt. 0) then

                    dimo = nmoltot
                    inddetc = ipc*(nelorb_diag*nelorb_diag + nelorb_at) + 1

                    if (ipc .eq. 1) then

                        call dgemm_my('N', 'N', nelorb_at, ndiff, dimo, 1.d0&
                            &, overs(1, nelorb_at + 1), nelorb_diag, detmat_c(inddetc), nelorb_diag, 0.d0, mat, nelorb_at&
                            &, nprocu, rankopt, commopt_mpi)

                        call dgemm_my('N', 'N', nelorb_at, ndiff, nelorb_at, 1.d0&
                            &, inv_sav, nelorb_at, mat, nelorb_at, 0.d0, detmat_proj(nelorb_c*nelorb_c + 1), nelorb_c&
                            &, nprocu, rankopt, commopt_mpi)

                        call dgemm_my('N', 'N', dimo, ndiff, dimo, 1.d0&
                            &, overs(nelorb_at + 1, nelorb_at + 1), nelorb_diag, detmat_c(inddetc), nelorb_diag&
                            &, 0.d0, umat, dimo, nprocu, rankopt, commopt_mpi)

                        !         Normalization

                    else

                        call zgemm_my('N', 'N', nelorb_at, ndiff, dimo, zone&
                            &, overs(1, nelorb_at + 1), nelorb_diag, detmat_c(inddetc), nelorb_diag, zzero, mat, nelorb_at&
                            &, nprocu, rankopt, commopt_mpi)

                        call zgemm_my('N', 'N', nelorb_at, ndiff, nelorb_at, zone&
                            &, inv_sav, nelorb_at, mat, nelorb_at, zzero, detmat_proj(2*nelorb_c*nelorb_c + 1), nelorb_c&
                            &, nprocu, rankopt, commopt_mpi)

                        call zgemm_my('N', 'N', dimo, ndiff, dimo, zone&
                            &, overs(2*nelorb_at + 1, nelorb_at + 1), nelorb_diag, detmat_c(inddetc), nelorb_diag&
                            &, zzero, umat, dimo, nprocu, rankopt, commopt_mpi)

                    end if

                    cost = 2.d0
                    imin = 1
                    do i = 1, ndiff
                        overo = ddot(ipc*dimo, umat(ipc*dimo*(i - 1) + 1), 1&
                                &, detmat_c(ipc*(nelorb_diag*(nelorb_diag + i - 1) + nelorb_at) + 1), 1)
                        costn = ddot(ipc*nelorb_at, mat(ipc*nelorb_at*(i - 1) + 1), 1, &
                                &detmat_proj(ipc*nelorb_c*(nelorb_c + i - 1) + 1), 1)/overo
                        if (costn .le. cost) then
                            cost = costn
                            imin = i
                        end if
                    end do

                    if (kaverage) then
#ifdef PARALLEL
                        allocate (agp_kp(nk), imin_kp(nk))
                        if (rankrep .eq. 0) then
                            call mpi_gather(cost, 1, MPI_DOUBLE_PRECISION, agp_kp,&
                                 &1, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
                            call mpi_gather(imin, 1, MPI_INTEGER, imin_kp,&
                                 &1, MPI_INTEGER, 0, commcolrep_mpi, ierr)
                        end if
                        if (rank .eq. 0) then
                            do ii = 1, nk
                                write (6, *) ' K#/ Min Overlap square unpaired   ', ii, imin_kp(ii), agp_kp(ii)
                                if (agp_kp(ii) .lt. 0.99999) write (6, *) ' ERROR Overlap K=', ii
                            end do
                        end if
                        deallocate (agp_kp, imin_kp)
#endif
                    else

                        if (rank .eq. 0) write (6, *) 'Minimum Overlap unpaired=', imin, cost

                    end if

                end if ! endif unpaired

                deallocate (eigmat, umat, mat, inv_sav)
                if (.not. symmagp .or. ipc .eq. 2 .or. ipf .eq. 2) deallocate (mats)

            elseif (detc_proj) then
                if (rank .eq. 0) write (6, *) ' Warning atomic orbitals found assuming&
                        & molecular  empty', nelorb_c, nelorb_diag
                sumdet = 0.d0
                do i = nelorb_at + 1, nelorb_c
                    do j = ipc*nelorb_at + 1, ipc*nelorb_c
                        sumdet = sumdet + abs(detmat_c(ipc*nelorb_c*(i - 1) + j))
                    end do
                end do
                if (sumdet .ne. 0.d0) then
                    if (rankopt .eq. 0) write (6, *) ' ERROR molecular not empty ', rank
#ifdef PARALLEL
                    call mpi_finalize(ierr)
#endif
                    stop
                end if
                detmat_proj = detmat_c

            end if ! endif big if detc_proj and sumdet=0

            if (detc_proj) then
                !       In this way we avoid a useless diagonalization that can
                !       affect the accuracy of the initial wf.
                deallocate (nozeroc_in, jbradetc_in, jbradetnc_in, sjbradet_in)
                if (allocated(nozerojc_in)) deallocate (nozerojc_in, jbrajc_in, jbrajnc_in)
                deallocate (buffer, overs, oversl, eigmol, molecorb)
                if (allocated(contrnorm)) deallocate (contrnorm)
                return
            end if

            !       write(6,*) ' Normalization basis contracted '
            do i = 1, nelorb_diag
                cost = overs(ipc*(i - 1) + 1, i)

                if (cost .gt. 0.d0) then
                    contrnorm(i) = cost
                    !        write(6,*) i,contrnorm(i)
                else

                    if (i .le. nelorb_at) then
                        if (rank .eq. 0) then
                            write (6, *) ' Warning atomic contracted #', i, ' singular norm ', overs(ipc*(i - 1) + 1, i)
                        end if
                    elseif (i .gt. nelorb_at .and. i .le. nelorb_at + nmolmax) then
                        if (rank .eq. 0) then
                            write (6, *) ' Warning molecular orbital  #', i, ' zero or singular ', overs(ipc*(i - 1) + 1, i)
                        end if
                    end if
                    contrnorm(i) = 0.d0
                end if
            end do

            !       write(6,*) ' Atomic normalizazion '
            !       do i=1,nelorb_at
            !       write(6,*) i,overs(i,i),contrnorm(i)
            !       enddo

            !       write(6,*) ' Input matrices overs  ',size(overs),nelorb_diag**2,nelorb2
            !       write(6,*) ' Input matrices detmat_c ',size(detmat_c),nelorb_diag**2

            !      define the metric of the Ghost's states
            do i = nelorb_diag + 1, nelorb_diagu
                overs(ipc*(i - 1) + 1, i) = 1.d0
                oversl(ipc*(i - nelorb_diag + ipf*nelorbh - 1) + 1, i) = 1.d0
            end do

            if (printoverlap .and. rank .eq. 0) then
                write (6, *) ' Overlap matrix ', nmoltot, nmolmax, nmol
                do i = 1, nelorb_diagu
                    do j = i, nelorb_diagu
                        write (6, *) i, j, overs(ipc*(i - 1) + 1:ipc*i, j)
                    end do
                end do
                write (6, *) ' Diagonal part '
                do i = 1, nelorb_diagu
                    write (6, *) i, overs(ipc*(i - 1) + 1:ipc*i, i)
                end do
            end if

            if (rank .eq. 0.) write (6, *) ' Eigenvalue det 1 '

            if (ndiff .eq. 0 .or. npar_eagp .gt. 0) then

                !              if(rank.eq.0.and.ipc.eq.2) then
                !              do i=1,nelorb_diag
                !                do j=1,nelorb_diag
                !                write(6,*) i,j,overs(2*j-1:2*j,i+nelorb_diag),overs(2*i-1:2*i,j+nelorb_diag)
                !                enddo
                !              enddo
                !            write(6,*) ' Detmat matrix '
                !              do i=1,nelorb_diag
                !                do j=1,nelorb_diag
                !                write(6,*) i,j,detmat_c(2*nelorb_diag*(j-1)+2*i-1:2*nelorb_diag*(j-1)+2*i)
                !                enddo
                !              enddo
                !             endif

                if (npar_eagp .gt. 0) then
                    call dcopy(ipc*nelorb_c*nelcol_c, detmat_c, 1, molecorb, 1)
                    !   TEST

                    !  write(6,*) ' Input detmat_c unpaired =', &
                    !  sum(abs(detmat_c(ipc*nelorb_c*(nelorb_c-1)+1:ipc*nelorb_c*nelcol_c)))
                    !
                    !             if(allocated(detmat)) deallocate(detmat)
                    !             allocate(detmat(ipc*ipf*nelorbh*(ipf*nelorbh+ndiff)))
                    !             detmat=0.d0
                    !         if(rank.eq.0) write(6,*) ' Input  scontract = ',nelorb,nelorbh,nelorbh+ndiff,nelcol&
                    !    &,nelcol_c,size(psip),ipf*nelorbh*nelcol_c
                    !        call scontract_mat_det(nelorbh,nelorb,nelcol,nelorb_c&
                    !       &,nelcol_c,detmat,detmat_c,mu_c,psip)
                    !        if(rank.eq.0) then
                    !        do i=1,ndiff
                    !          cost=0.d0
                    !          do k=1,ipf*nelorbh
                    !          cost=cost+detmat(ipf*nelorbh*ipc*(ipf*nelorbh+i-1)+k)
                    !          enddo
                    !        write(6,*) 'Unpaired correct Unscaled = ',i,cost
                    !        enddo
                    !        cost=0.d0
                    !        do i=1,ipf*nelorbh
                    !          do k=i+1,ipf*nelorbh
                    !          cost=cost+detmat(ipf*nelorbh*ipc*(i-1)+k)
                    !          enddo
                    !        enddo
                    !        write(6,*) 'Geminal  correct Unscaled = ',cost
                    !        endif
                    !   END TEST

                    !       call mpi_finalize(ierr)
                    !       stop

                    deallocate (detmat_c)
                    allocate (detmat_c(ipc*nelcol_c*nelcol_c))
                    call copy_eagp(.true., ipc, nelorb_c, nelcol_c, molecorb, eagp_pfaff, detmat_c)
                end if
                !                          if(rank.eq.0) then
                !                     write(6,*) ' Overlap matrix input ',nelorb_diagu,nelcol_c
                !                         do i=1,nelorb_diagu
                !                           do j=1,nelorb_diagu
                !                         write(6,*) i,j,overs(i,j),detmat_c(nelcol_c*(j-1)+i)
                !                           enddo
                !                         enddo
                !                         endif

                call eval_molec_epsdgel(nelorb_diagu, overs, detmat_c, molecorb&
                        &, eigmol, nelorb_diagu, epsdgm, nprocu, rank, rankopt, commopt_mpi, 1, symmagp)
                call dcopy(nelorb2, molecorb, 1, overs, 1)

                !             if(rank.eq.0) write(6,*) ' Output overs =',sum(abs(overs(:,:)))

                if (printoverlap .and. rank .eq. 0) then
                    if (ipf .eq. 2 .or. (ipc .eq. 1 .and. symmagp)) then
                        write (6, *) ' Molecular orbitals in the contracted basis '
                        do i = 1, nelorb_diagu
                            write (6, *) ' Eigenvalue = ', i, eigmol((i - 1)/ipf + 1)
                            do j = 1, nelorb_diagu
                                write (6, *) j, overs(ipc*(j - 1) + 1:ipc*j, i)
                            end do
                        end do
                    else
                        write (6, *) ' Left/right  Molecular orbitals in the contracted basis '
                        do i = 1, nelorb_diag
                            write (6, *) ' Eigenvalue = ', eigmol(i)
                            do j = 1, nelorb_diag
                                write (6, *) j, overs(ipc*(j - 1) + 1:ipc*j, i), overs(ipc*(j - 1) + 1:ipc*j, i + nelorb_diag)
                            end do
                        end do
                    end if
                end if

                !                 write(6,*) ' Eigenvectors ',nelorb_diag,nelorb_c
                !                 do i=1,nelorb_diag
                !                 write(6,*) i,sum(overs(1:nelorb_diag,i))
                !                 enddo

                !       scontract molecular orbitals ipf=2 is similar to symmagp=true
                !   since no distintion between left and right eigenvect. is possible
                if (ipf .eq. 2) then

                    if (ipc .eq. 1) then

                        call dgemm_my('N', 'N', 2*nelorbh, nelorb_diagu, nelorb_diag, 1.d0, mu_c, 2*nelorbh&
                                &, overs, nelorb_diagu, 0.d0, molecorb, nelorbpf, nprocu, rankopt, commopt_mpi)

                    else

                        call zgemm_my('N', 'N', 2*nelorbh, nelorb_diagu, nelorb_diag, zone, mu_c, 2*nelorbh&
                                &, overs, nelorb_diagu, zzero, molecorb, nelorbpf, nprocu, rankopt, commopt_mpi)

                    end if

                    if (npar_eagp .gt. 0) then
                        do i = 1, nelorb_diagu
                            do j = nelorb_diag + 1, nelorb_diagu
                                ind = ipc*(j - nelorb_diag + ipf*nelorbh - 1) + 1
                                molecorb(ind:ind + ipc - 1, i) = overs(ipc*(j - 1) + 1:ipc*j, i)
                            end do
                        end do
                    end if

                elseif (symmagp .and. ipc .eq. 1) then
                    call dgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, 1.d0, mu_c, nelorbh&
                            &, overs, nelorb_diag, 0.d0, molecorb, nelorbh, nprocu, rankopt, commopt_mpi)
                else
                    !       Odd eigenvectors are right eigenvectors

                    !       write(6,*) ' Input molecorb '
                    !       do i=1,nelorb_diag
                    !       write(6,*) ' Left eigenvectors =',overs(1:nelorb_diag,2*i)
                    !       write(6,*) ' Right eigenvectors =',overs(1:nelorb_diag,2*i-1)
                    !       enddo

                    if (ipc .eq. 1) then
                        call dgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, 1.d0, mu_c, nelorbh&
                                &, overs, 2*nelorb_diag, 0.d0, molecorb, 2*nelorbh, nprocu, rankopt, commopt_mpi)
                        !       Even  eigenvectors are left eigenvectors
                        call dgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, 1.d0, mu_c, nelorbh&
                                &, overs(1, 2), 2*nelorb_diag, 0.d0, molecorb(1, 2), 2*nelorbh, nprocu, rankopt, commopt_mpi)
                    else
                        call zgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, zone, mu_c, nelorbh&
                                &, overs, 2*nelorb_diag, zzero, molecorb, 2*nelorbh, nprocu, rankopt, commopt_mpi)
                        !       Even  eigenvectors are left eigenvectors
                        call zgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, zone, mu_c, nelorbh&
                                &, overs(1, 2), 2*nelorb_diag, zzero, molecorb(1, 2), 2*nelorbh, nprocu, rankopt, commopt_mpi)
                    end if

                    !       write(6,*) ' Output  molecorb scontracted  '
                    !       do i=1,nelorb_diag
                    !       write(6,*) ' Eigenvalue ',i,eigmol(i)
                    !       write(6,*) ' Left eigenvectors =',sum(molecorb(1:nelorbh,2*i))
                    !       write(6,*) ' Right eigenvectors =',sum(molecorb(1:nelorbh,2*i-1))
                    !       enddo

                end if

                if (nelorb_diagu .ne. nelorb_c) then
                    deallocate (detmat_c)
                    allocate (detmat_c(ipc*nelcol_c*nelorb_c))
                    detmat_c = 0.d0
                    indunp = ipc*nelorb_c*nelorb_c + 1
                end if

            else ! ndiff=0, ndiff > 0 below
                !  othorgonalize the AGP part with respect to the unpaired
                !  then orthogonalize all and fill properly the unpaired
                !       for test

                indunp = ipc*nelorb_diag*nelorb_diag + 1

                if (symmagp .and. ipc .eq. 1 .or. ipf .eq. 2) then
                    if (allocated(molecorb)) deallocate (molecorb)
                    allocate (molecorb(ipc*nelorb_diag, nelorb_diag))
                else
                    if (allocated(molecorb)) deallocate (molecorb)
                    allocate (molecorb(ipc*nelorb_diag, 2*nelorb_diag))
                end if
                molecorb = 0.d0

                !       if the input is a det the output will be a det with orth orbitals.

                call eval_molec_unpaired(nelorb_diag, overs, detmat_c, molecorb&
                        &, eigmol, nelorb_diag, epsdgm, nprocu, rank, rankopt, commopt_mpi, ndiff&
                        &, detmat_c(indunp), orthoyes, symmagp)

                overs = molecorb

                deallocate (molecorb)
                if (symmagp .and. ipc .eq. 1 .or. ipf .eq. 2) then
                    allocate (molecorb(ipf*nelorbh, nelorb_diag))
                else
                    allocate (molecorb(ipc*nelorbh, 2*nelorb_diag))
                end if
                !       scontract molecular orbitals

                if (ipf .eq. 2) then
                    if (ipc .eq. 1) then

                        call dgemm_my('N', 'N', 2*nelorbh, nelorb_diag, nelorb_diag, 1.d0, mu_c, 2*nelorbh&
                                &, overs, nelorb_diag, 0.d0, molecorb, 2*nelorbh, nprocu, rankopt, commopt_mpi)

                    else

                        call zgemm_my('N', 'N', 2*nelorbh, nelorb_diag, nelorb_diag, zone, mu_c, 2*nelorbh&
                                &, overs, nelorb_diag, zzero, molecorb, 2*nelorbh, nprocu, rankopt, commopt_mpi)

                    end if

                elseif (symmagp .and. ipc .eq. 1) then

                    call dgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, 1.d0, mu_c, nelorbh&
                            &, overs, nelorb_diag, 0.d0, molecorb, nelorbh, nprocu, rankopt, commopt_mpi)

                else

                    if (ipc .eq. 1) then
                        !       Odd eigenvectors are right eigenvectors
                        call dgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, 1.d0, mu_c, nelorbh&
                                &, overs, 2*nelorb_diag, 0.d0, molecorb, 2*nelorbh, nprocu, rankopt, commopt_mpi)
                        !       Even  eigenvectors are left eigenvectors
                        call dgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, 1.d0, mu_c, nelorbh&
                                &, overs(1, 2), 2*nelorb_diag, 0.d0, molecorb(1, 2), 2*nelorbh, nprocu, rankopt, commopt_mpi)
                    else
                        !       Odd eigenvectors are right eigenvectors
                        call zgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, zone, mu_c, nelorbh&
                                &, overs, 2*nelorb_diag, zzero, molecorb, 2*nelorbh, nprocu, rankopt, commopt_mpi)
                        !       Even  eigenvectors are left eigenvectors
                        call zgemm_my('N', 'N', nelorbh, nelorb_diag, nelorb_diag, zone, mu_c, nelorbh&
                                &, overs(1, 2), 2*nelorb_diag, zzero, molecorb(1, 2), 2*nelorbh, nprocu, rankopt, commopt_mpi)

                    end if

                end if

                !  Scontracting the unpaired stored in detmat_c

                if (ipf .eq. 2) then

                    if (ipc .eq. 1) then
                        call dgemm_my('N', 'N', 2*nelorbh, ndiff, nelorb_diag, 1.d0, mu_c, 2*nelorbh&
                                &, detmat_c(indunp), nelorb_diag, 0.d0, molecorb(1, nelorb_diag - nmoltot + 1)&
                                &, 2*nelorbh, nprocu, rankopt, commopt_mpi)

                    else

                        call zgemm_my('N', 'N', 2*nelorbh, ndiff, nelorb_diag, zone, mu_c, 2*nelorbh&
                                &, detmat_c(indunp), nelorb_diag, zzero, molecorb(1, nelorb_diag - nmoltot + 1)&
                                &, 2*nelorbh, nprocu, rankopt, commopt_mpi)

                    end if

                elseif (symmagp .and. ipc .eq. 1) then
                    call dgemm_my('N', 'N', nelorbh, ndiff, nelorb_diag, 1.d0, mu_c, nelorbh&
                            &, detmat_c(indunp), nelorb_diag, 0.d0, molecorb(1, nelorb_diag - nmoltot + 1)&
                            &, nelorbh, nprocu, rankopt, commopt_mpi)
                else
                    if (ipc .eq. 1) then
                        call dgemm_my('N', 'N', nelorbh, ndiff, nelorb_diag, 1.d0, mu_c, nelorbh&
                                &, detmat_c(indunp), nelorb_diag, 0.d0, molecorb(1, 2*nelorb_diag - nmoltot + 1)&
                                &, nelorbh, nprocu, rankopt, commopt_mpi)
                    else
                        call zgemm_my('N', 'N', nelorbh, ndiff, nelorb_diag, zone, mu_c, nelorbh&
                                &, detmat_c(indunp), nelorb_diag, zzero, molecorb(1, 2*nelorb_diag - nmoltot + 1)&
                                &, nelorbh, nprocu, rankopt, commopt_mpi)
                    end if
                end if

                if (nelorb_c .ne. nelorb_diag) then
                    deallocate (detmat_c)
                    allocate (detmat_c(ipc*nelcol_c*nelorb_c))
                    detmat_c = 0.d0
                    indunp = ipc*nelorb_c*nelorb_c + 1
                end if
            end if ! end if ndiff

            !CCC Fine di ndiff.eq.0

            if (ndiff .gt. 0 .and. npar_eagp .eq. 0) then
                do i = 1, (ndiff + ipf - 1)/ipf
                    if (ipf .eq. 2 .and. nelorb_diag/2 + i - nmol - (ndiff + 1)/2 .gt. 0) then
                        eigmol(nelorb_diag/2 - nmol - (ndiff + ipf - 1)/ipf + i) = 0.d0
                        !         write(6,*) ' VANISHED TERM = ',nelorb_diag,nmol,nelorb_diag/2-nmol-(ndiff+ipf-1)/ipf+i
                    elseif (nelorb_diag + i - nmol - ndiff .gt. 0) then
                        eigmol(nelorb_diag - nmol - ndiff + i) = 0.d0

                    else
                        if (rank .eq. 0) write (6, *) ' ERROR too many molecular orbitals, pls. decrease nmol  !!! '
#ifdef  PARALLEL
                        call mpi_finalize(ierr)
#endif
                        stop
                    end if
                end do
            end if ! ndiff > 0
        else ! epsdgm.lt.0 --> no diagonalization nor evaluation of overlaps

            eigmol = 1.d0

            ! generation random molecular orbitals below

            if ((ipc .eq. 2 .or. .not. symmagp) .and. ipf .eq. 1) then
                do j = max(2*nelorb_diag - nmoltot + 1, 1), max(2*nelorb_diag, nmoltot)
                    if (.not. yes_complex) then
                        do i = 1, nelorbh
                            molecorb(i, j) = ran(iseed) - 0.5d0
                        end do
                    else
                        do i = 1, 2*nelorbh, 2
                            molecorb(i, j) = ran(iseed) - 0.5d0
                            !molecorb(i+1,j)=rannum+0.5d0
                            molecorb(i + 1, j) = 0.d0
                        end do
                    end if
                end do
            else
                do j = max(nelorb_diag - nmoltot + 1, 1), max(nelorb_diag, nmoltot)
                    do i = 1, ipf*nelorbh
                        molecorb(i, j) = ran(iseed) - 0.5d0
                    end do
                end do
            end if
        end if

        ! * * * *  end evaluation of MOs  * * * * !

        ! Now updating the wave function using same conventions as in fort10_io

        if (allocated(ipsip)) deallocate (ipsip)
        allocate (ipsip(nshell_c))
        ipsip = 0

        nshell_c = nshell_c - molecular ! atomic contracted shells
        nshellmol = nshell_c + nmoltot ! atomic + molecular contracted shells

        ipsip = mult_c ! mult_c = array which contains multiplicty of all basis orbitals
        deallocate (mult_c)
        allocate (mult_c(nshellmol))
        mult_c(1:nshell_c) = ipsip(1:nshell_c)
        mult_c(nshell_c + 1:nshellmol) = 1 ! multiplicity of MOs is always 1
        ipsip = nparam_c
        deallocate (nparam_c)
        allocate (nparam_c(nshellmol))
        nparam_c(1:nshell_c) = ipsip(1:nshell_c)

        nparam_c(nshell_c + 1:nshellmol) = 2*nelorbh*ipf ! # of parameters in the MOs basis
        ! first set of params = pointers to basis orbitals
        ! second set of params = coefficients
        ipsip = ioptorb_c
        deallocate (ioptorb_c)
        allocate (ioptorb_c(nshellmol))
        ioptorb_c(1:nshell_c) = ipsip(1:nshell_c)
        ioptorb_c(nshell_c + 1:nshellmol) = 1000000
        ipsip = kion_c
        deallocate (kion_c)
        allocate (kion_c(nshellmol))
        kion_c(1:nshell_c) = ipsip(1:nshell_c)
        kion_c(nshell_c + 1:nshellmol) = 1
        ! reallocating dup_c to include molecular orbitals
        ! MOs are considered as contracted orbitals

        ! ****
        ! CONVENTION for complex algorithm :
        ! index are left unchanged while allocation is doubled!!!
        ! ****

        ! updating matrix of coefficients dup_c
        deallocate (psip)
        if (.not. yes_complex) then
            allocate (psip(iesupr_c))
            psip(1:iesupr_c) = dup_c(1:iesupr_c)
        else
            allocate (psip(2*iesupr_c))
            psip(1:2*iesupr_c) = dup_c(1:2*iesupr_c)
        end if
        deallocate (dup_c)

        if (molecular .eq. 0) then
            iesup_cmol = iesup_c + 2*nmoltot*nelorbh*ipf
        else
            iesup_cmol = iesup_c + 2*(nmoltot - molecular)*nelorbh*ipf
        end if

        ! following same convention as in fort10_io
        iesupr_cmol = iesup_cmol
        !
        if (.not. yes_complex) then
            allocate (dup_c(iesupr_cmol))
        else
            allocate (dup_c(2*iesupr_cmol))
        end if

        call update_dup_c

        ! ----- end update coefficients matrix dup_c ----- !

        deallocate (ipsip)
        allocate (ipsip(occ_c))
        ipsip = ioccup_c
        deallocate (ioccup_c)
        occ_cmol = occ_c + nmoltot - molecular
        allocate (ioccup_c(occ_cmol))
        if (occ_cmol .gt. occ_c) then
            if (only_molecular) then
                ioccup_c(1:occ_c) = 0
            else
                ioccup_c(1:occ_c) = ipsip(1:occ_c)
            end if
            ioccup_c(occ_c + 1:occ_cmol) = 1
        else
            ioccup_c(1:occ_cmol) = ipsip(1:occ_cmol)
        end if

        contraction = contraction + nmoltot - molecular
        ! change here

        if (allocated(detmat_c)) deallocate (detmat_c)
        if (allocated(nozero_c)) deallocate (nozero_c)
        deallocate (jbradetn, sjbradet)
        if (only_molecular) then
            nelorb_c = 0
            nelorb_at = 0
        else
            if (molecular .gt. 0) nelorb_c = nelorb_c - molecular
        end if
        nelorb_c = nelorb_c + nmoltot
        nelcol_c = nelorb_c + ndiff

        ! new allocation of detmat_c to store MOs
        if (.not. yes_complex) then
            allocate (detmat_c(nelorb_c*nelcol_c))
        else
            allocate (detmat_c(2*nelorb_c*nelcol_c))
        end if
        detmat_c = 0.d0

        !
        if (npar_eagp .gt. 0) then
            allocate (nozero_c(npar_eagp))
            nozero_c(1:npar_eagp) = jbradet(nnozeroc_in + 1:nnozeroc_in + npar_eagp)
        end if
        deallocate (jbradet)
        !
        if (only_molecular) then
            nnozeroc_in = 0
            iesswrc_in = 0
            molecsw = 0
            molecsym = 0
        else
            nnozeroc_in = nnozeroc_in - molecsw
            iesswrc_in = iesswr - molecsym
        end if

        nnozero_c = nnozeroc_in + nmol + ndiff

        iesswr = iesswrc_in + nmol + ndiff
        if (add_offmol .and. nmolmax .gt. nmolmin) then
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                nmol_add = ((nmolmax - nmolmin + 1)*(nmolmax - nmolmin))/2
            else
                if (ipf .eq. 2) then
                    nmol_add = (nmolmax - nmolmin)*(2*(nmolmax - nmolmin) - 1)
                else
                    if (.not. symmagp) then
                        nmol_add = (nmolmax - nmolmin + 1)*(nmolmax - nmolmin)
                    else
                        nmol_add = ((nmolmax - nmolmin + 1)*(nmolmax - nmolmin))/2
                    end if
                end if
            end if
            nnozero_c = nnozero_c + nmol_add
            iesswr = iesswr + nmol_add
        end if
        allocate (jbradet(nnozero_c + npar_eagp), sjbradet(nnozero_c)&
                &, jbradetn(3*nnozero_c))

        if (npar_eagp .gt. 0) then
            jbradet(nnozero_c + 1:nnozero_c + npar_eagp) = nozero_c(1:npar_eagp)
            deallocate (nozero_c)
        end if

        allocate (nozero_c(nnozero_c))
        nozero_c = 0
        if (allocated(ipsip)) deallocate (ipsip)
        allocate (ipsip(3*nnozero_c))
        sjbradet = .false.
        jbradet = 0
        jbradetn = 0

        if (.not. only_molecular) then
            dimvect = nnozeroc_in + molecsw
            allocate (indexo(dimvect))
            allocate (sortvect(dimvect + 1))
            sortvect = 0
            sortvect(1:dimvect) = abs(jbradetc_in(1:dimvect))
            call dsortx(sortvect, 1, dimvect, indexo)
            icount = 1
            ! individuation indnn
            indnn = 0
            indnn_in = 0
            do i = 1, iesswrc_in + molecsym
                ind = 0
                ! To be optimized sort jbradetc_in
                ! skip all not equal
                do while (sortvect(icount) .lt. i .and. icount .le. dimvect)
                    icount = icount + 1
                end do
                ! check the very last
                if (icount .le. dimvect) then
                    do while (sortvect(icount) .eq. i .and. icount .le. dimvect)
                        ind = ind + 1
                        icount = icount + 1
                    end do
                    icount = icount - 1
                end if

                if (ind .eq. 0) then
                    indnn = indnn + 1
                    ! load ipsip
                    ipsip(1) = jbradetnc_in(indnn)
                    ii = ipsip(1)
                    yest = .true.
                    do j = 1, -2*ii, 2
                        ipsip(j + 1) = jbradetnc_in(j + indnn)
                        ipsip(j + 2) = jbradetnc_in(j + 1 + indnn)
                        if ((abs(ipsip(j + 1)) .gt. nelorbc_in - molecular&
                                &.and. abs(ipsip(j + 1)) .le. nelorbc_in)&
                                &.or. (abs(ipsip(j + 2)) .gt. nelorbc_in - molecular .and.&
                                        &abs(ipsip(j + 2)) .le. nelorbc_in)) yest = .false.
                    end do
                    if (yest) then
                        do j = 1, -2*ii + 1
                            if (ipsip(j) .gt. nelorbc_in .and. j .ne. 1) ipsip(j) = ipsip(j) + nelorb_c - nelorbc_in
                            jbradetn(indnn_in + j) = ipsip(j)
                        end do
                        indnn_in = indnn_in + 1 - 2*ii
                    end if
                    indnn = indnn - 2*ii
                end if
            end do

            deallocate (indexo, sortvect)

            ind = 0
            do i = 1, nnozeroc_in + molecsw
                iy = (nozeroc_in(i) - 1)/nelorbc_in + 1
                ix = nozeroc_in(i) - (iy - 1)*nelorbc_in
                if (ix .le. nelorbc_in - molecular .and. iy .le. nelorbc_in - molecular) then
                    ind = ind + 1
                    jbradet(ind) = jbradetc_in(i)
                    sjbradet(ind) = sjbradet_in(i)
                    nozero_c(ind) = (iy - 1)*nelorb_c + ix
                elseif (ix .le. nelorbc_in - molecular .and. iy .gt. nelorbc_in) then
                    ind = ind + 1
                    jbradet(ind) = jbradetc_in(i)
                    sjbradet(ind) = sjbradet_in(i)
                    nozero_c(ind) = (iy + nelorb_c - nelorbc_in - 1)*nelorb_c + ix
                end if
            end do

            ! now sort jbradet
            dimvect = ind
            allocate (indexo(dimvect), indexn(dimvect))
            allocate (sortvect(dimvect + 1)) ! Important ovedundant allocation
            sortvect = 0
            sortvect(1:dimvect) = abs(jbradet(1:dimvect))
            indexn(1:dimvect) = jbradet(1:dimvect) ! save the value...
            jbradet = 0
            inds = 0
            i = 0
            call dsortx(sortvect, 1, dimvect, indexo)
            icount = 1
            do while (inds .lt. iesswrc_in .and. i .lt. ind)
                i = i + 1
                indj = 0
                ! Go to the first |jbradet|<i
                do while (sortvect(icount) .lt. i .and. icount .le. dimvect)
                    icount = icount + 1
                end do
                ! If you are not at the end count how many |jbradet| =i
                if (icount .le. dimvect) then
                    do while (sortvect(icount) .eq. i .and. icount .le. dimvect)
                        indj = indj + 1
                        if (indj .eq. 1) indfirst = icount
                        icount = icount + 1
                    end do
                    ! At the end of the loop the condition |jbradet| /=i and icount=icount-1
                    icount = icount - 1
                    indlast = icount
                    if (indj .ne. 0) then
                        inds = inds + 1
                        do iicount = indfirst, indlast
                            j = indexo(iicount)
                            if (indexn(j) .gt. 0) then
                                jbradet(j) = inds
                            elseif (indexn(j) .lt. 0) then
                                jbradet(j) = -inds
                            end if
                        end do
                    end if
                end if
            end do

            deallocate (sortvect, indexo, indexn)

        end if

        if (only_molecular) then
            ind = 0
            inds = 0
            indnn = 0
        else
            ind = nnozeroc_in ! # of non zero values of detmat_c
            ! already filled by contracted atomic orbitals
            inds = iesswrc_in
            indnn = indnn_in
            ! updating detmat_c
        end if
        do i = 1, nmol
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                ix = nelorb_c - nmoltot + i
                iy = ix
            else
                ix = nelorb_c - nmoltot + 2*i
                iy = ix - 1
            end if
            ind = ind + 1
            nozero_c(ind) = nelorb_c*(iy - 1) + ix
            ! index to distinguish between real and complex wfs
            ind_c = nozero_c(ind)
            ind_t = nelorb_c*(ix - 1) + iy

            if (i .le. nmolmin) then ! the case for a perfect SD with nmol_min=nmol_max=N/2
                if (.not. yes_complex) then
!                   detmat_c(ind_c) = dble(nmolmin - i + 1)
                    detmat_c(ind_c) = 1.d0
                else
!                   detmat_c(2 * ind_c - 1) = dble(nmolmin - i + 1)
                    detmat_c(2*ind_c - 1) = 1.d0
                    detmat_c(2*ind_c) = 0.d0
                end if
            elseif (i .le. nmolmaxw) then

                if (power .eq. 1.d0) then
                    if (.not. yes_complex) then
                        detmat_c(ind_c) = eigmol(nelorb_diagu/ipf - i + 1)/eigmol(nelorb_diagu/ipf - nmolmin + 1)
                    else
                        detmat_c(2*ind_c - 1) = eigmol(nelorb_diagu/ipf - i + 1) &
                                                /eigmol(nelorb_diagu/ipf - nmolmin + 1)
                        detmat_c(2*ind_c) = 0.d0
                    end if
                else
                    cost = eigmol(nelorb_diagu/ipf - i + 1)/eigmol(nelorb_diagu/ipf - nmolmin + 1)
                    if (.not. yes_complex) then
                        detmat_c(ind_c) = abs(cost)**power*sign(1.d0, cost)
                    else
                        detmat_c(2*ind_c - 1) = abs(cost)**power*sign(1.d0, cost)
                        detmat_c(2*ind_c) = 0.d0
                    end if
                end if
            else
                if (.not. yes_complex) then
                    detmat_c(ind_c) = 0.d0
                else
                    detmat_c(2*ind_c - 1) = 0.d0
                    detmat_c(2*ind_c) = 0.d0
                end if
            end if
            !      Fullfill antisymmetry
            if (ipf .eq. 2) then
                if (ipc .eq. 1) then
                    detmat_c(ind_t) = -detmat_c(ind_c)
                else
                    detmat_c(2*ind_t - 1) = -detmat_c(2*ind_c - 1)
                    detmat_c(2*ind_t) = -detmat_c(2*ind_c)
                end if
            end if
            indnn = indnn + 1
            jbradetn(indnn) = -1
            jbradetn(indnn + 1) = ix
            jbradetn(indnn + 2) = iy
            indnn = indnn + 2
            !      psip(i)=detmat_c(ipc*(ind_c-1)+1)
        end do
        if (add_offmol .and. nmolmax .gt. nmolmin) then
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                do i = nmolmin, nmolmax
                    do j = i + 1, nmolmax
                        ix = nelorb_c - nmoltot + j
                        iy = nelorb_c - nmoltot + i ! iy<ix
                        ind = ind + 1
                        nozero_c(ind) = nelorb_c*(iy - 1) + ix
                        ! index to distinguish between real and complex wfs
                        indnn = indnn + 1
                        jbradetn(indnn) = -1
                        jbradetn(indnn + 1) = ix
                        jbradetn(indnn + 2) = iy
                        indnn = indnn + 2
                    end do
                end do
            else
                if (ipf .eq. 2) then
                do i = 2*nmolmin, 2*nmolmax
                    do j = i + 2, 2*nmolmax
                        ix = nelorb_c - nmoltot + j
                        iy = nelorb_c - nmoltot + i ! iy< ix-1
                        ind = ind + 1
                        nozero_c(ind) = nelorb_c*(iy - 1) + ix
                        ! index to distinguish between real and complex wfs
                        indnn = indnn + 1
                        jbradetn(indnn) = -1
                        jbradetn(indnn + 1) = ix
                        jbradetn(indnn + 2) = iy
                        indnn = indnn + 2
                    end do
                end do
                else
                if (.not. symmagp) then
                do i = nmolmin, nmolmax
                    do j = nmolmin, nmolmax
                    if (i .ne. j) then
                        ix = nelorb_c - nmoltot + 2*i
                        iy = nelorb_c - nmoltot + 2*j - 1
                        ind = ind + 1
                        nozero_c(ind) = nelorb_c*(iy - 1) + ix
                        ! index to distinguish between real and complex wfs
                        indnn = indnn + 1
                        jbradetn(indnn) = -1
                        jbradetn(indnn + 1) = ix
                        jbradetn(indnn + 2) = iy
                        indnn = indnn + 2
                    end if
                    end do
                end do
                else
                do i = nmolmin, nmolmax
                    do j = i + 1, nmolmax
                        ix = nelorb_c - nmoltot + 2*j
                        iy = nelorb_c - nmoltot + 2*i - 1
                        ind = ind + 1
                        nozero_c(ind) = nelorb_c*(iy - 1) + ix
                        ! index to distinguish between real and complex wfs
                        indnn = indnn + 1
                        jbradetn(indnn) = -1
                        jbradetn(indnn + 1) = ix
                        jbradetn(indnn + 2) = iy
                        indnn = indnn + 2
                    end do
                end do
                end if
                end if
            end if
        end if
        !  TEST uncomment also psip line above.
        !    allocate(detmat_proj(ipc*nelorb_c*nelorb_c))
        ! !       Compute detmat_proj
        !    detmat_proj=0.d0
        !     if(rank.eq.0) write(6,*) ' Before detmat_proj  overs =',nelorb_diag,sum(abs(overs(:,:)))
        !     if(symmagp.and.ipc.eq.1) then
        !        do i=1,nmol
        !           do j=1,nelorb_at
        !              molecorb(j,i)=molecorb(j,i)*psip(i)
        !           enddo
        !        enddo
        ! !
        !        call dgemm_my('N','T',nelorb_at,nelorb_at,nmol,1.d0,molecorb,nelorb_at&
        !             &,molecorb,nelorbh,0.d0,detmat_proj,nelorb_c,nprocu,rankopt,commopt_mpi)
        !     else
        !        do i=1,nmol
        !           if(rank.eq.0) write(6,*) 'psip =',i,psip(i)
        !           do j=1,ipc*nelorb_at
        !           molecorb(j,i)=overs(j,2*nelorb_diag-2*i+2)*psip(i)
        !          molecorb(j,i+nmol)=overs(j,2*nelorb_diag-2*i+1)
        !           enddo
        !        enddo
        !        if(ipc.eq.2) then
        ! !      Here I am consistent with the definition without the complex conjugate
        !        call zgemm_my('N','T',nelorb_at,nelorb_at,nmol,zone,molecorb,nelorbh&
        !  &,molecorb(1,nmol+1),nelorbh,zzero,detmat_proj,nelorb_c,nprocu,rankopt,commopt_mpi)
        !        else
        !        call dgemm_my('N','T',nelorb_at,nelorb_at,nmol,1.d0,molecorb,nelorbh&
        !  &,molecorb,nelorbh,0.d0,detmat_proj,nelorb_c,nprocu,rankopt,commopt_mpi)
        !        endif
        !     endif
        !! end first part of  TEST

        ! unpaired orbitals
        do i = 1, ndiff
            ind = ind + 1
            j = nelorb_c - ndiff + i
            nozero_c(ind) = nelorb_c*(nelorb_c + i - 1) + j
            if (.not. yes_complex) then
                detmat_c(nozero_c(ind)) = 1.d0
            else
                detmat_c(2*nozero_c(ind) - 1) = 1.d0
                detmat_c(2*nozero_c(ind)) = 0.d0
            end if
            indnn = indnn + 1
            jbradetn(indnn) = -1
            jbradetn(indnn + 1) = j
            jbradetn(indnn + 2) = nelorb_c + i
            indnn = indnn + 2
        end do

        deallocate (ipsip)

        !  Do not change jbradet from the chosen input value.
        if (size(jbradet) .eq. size(jbradetc_in) .and. epsdgm .ge. 0.d0) then
            !       if(rankopt.eq.0) write(6,*) ' Warning jbradet,nozero unchanged '
            dimvect = min(size(jbradetn), size(jbradetnc_in))
            jbradet = jbradetc_in
            sjbradet = sjbradet_in
            jbradetn(1:dimvect) = jbradetnc_in(1:dimvect)
            if (contraction .eq. 0) then
                nozero = nozeroc_in
            else
                nozero_c = nozeroc_in
            end if
        end if

        allocate (ipsip(iesup_c + iesupind))
        iesupindmol = iesupind
        ipsip = jbraiesup_sav
        deallocate (jbraiesup_sav)
        allocate (jbraiesup_sav(iesup_cmol + iesupindmol))

        ! calculation of ind
        ind = 0
        do i = 1, iesupind
            ind = ind + 1
            ii = abs(ipsip(ind))
            ind = ind + ii
        end do
        jbraiesup_sav(1:ind) = ipsip(1:ind)
        if (molopt .ge. 2 .and. contraction_in .ne. 0 .and. .not. detc_proj) then
            ! do not optimize the coefficients of the contracted orbitals
            ! But only if molopt>0 otherwise optimize also the coefficients
            ind = 0
            do i = 1, iesupind
                ind = ind + 1
                ii = jbraiesup_sav(ind)
                do j = 1, abs(ii)
                    if (iesuptransb(jbraiesup_sav(ind + j)) .eq. 0 .and. ii .gt. 0) &
                            & jbraiesup_sav(ind) = -abs(jbraiesup_sav(ind))
                end do
                ind = ind + abs(ii)
            end do
        end if

        indpar = iesup_c

        nshell_c = nshellmol
        iesupr_c = iesupr_cmol
        iesup_c = iesup_cmol
        occ_c = occ_cmol
        iesupind = iesupindmol

        ! definition molecular orbitals
        deallocate (overs, oversl, molecorb, eigmol)
        ! recompute and reallocate mu_c in ANY case
        ! recompute also multranspip and transpip

        call upmuctranspip

        deallocate (contrnorm)

        ! ---------- END molecular orbitals on the determinant --------- !

        !!!!   nshell_c --> nshell_c+nmol   redefine mu_c nelorb_c
        !!!!   ioptorb_c, dup_c nparam_c contraction>0 nozero_c(i),detmat_c
        !!!    iesswr jbradet jbradetn
        !!!    jbraiesup_sav
        !!!    the unpaired substituted with conventional molecular orbital

        deallocate (nozeroc_in, jbradetc_in, jbradetnc_in, sjbradet_in)
        deallocate (buffer)

        if (nmoltot .gt. molecular) then
            iscraipsip = max(iscraipsip, ipc*iesupr_c + ipc*iesupind, npar3bodyr_c + iesmind)
            molecular = nmoltot !  redefine the number or molecular
        end if

        timechange = cclock() - timechange

        if (rank .eq. 0) write (6, *) ' Time change fort.10 =', timechange
        !!  TEST
        !!  From real to effective
        !    if(rank.eq.0) write(6,*) ' before attach ',nelorb_c,2*nnozero_c,size(dsw),size(psip),size(ipsip)
        !    call update_kiontot
        !    deallocate(dsw)
        !    allocate(dsw(2*nnozero_c))
        !    deallocate(psip,ipsip)
        !    allocate(psip(2*nnozero_c),ipsip(2*nnozero_c))
        !    if(allowed_averagek) call attach_phase2det(.false.,detmat_proj)
        !    if(rank.eq.0) write(6,*) ' after  attach '
        !    if(rank.eq.0) write(6,*) ' Output detmat_proj '
        !
        !       cost=0.d0
        !       do i=1,nelorb_at
        !        do j=i,nelorb_at
        ! if(rank.eq.0)  write(6,*) i,j,detmat_proj(ipc*nelorb_c*(j-1)+2*i-1:ipc*nelorb_c*(j-1)+2*i),&
        !        &detmat_proj(ipc*nelorb_c*(i-1)+2*j-1:ipc*nelorb_c*(i-1)+2*j)
        !        cost=cost+sum(abs(detmat_proj(ipc*nelorb_c*(j-1)+2*i-1:ipc*nelorb_c*(j-1)+2*i)))+&
        !        &sum(abs(detmat_proj(ipc*nelorb_c*(i-1)+2*j-1:ipc*nelorb_c*(i-1)+2*j)))
        !        enddo
        !       enddo
        !       if(rank.eq.0) write(6,*) ' sum rule ',cost
        !           call  symmetrizeagp(nnozero_c,nozero_c,jbradet,sjbradet,jbradetn,dsw&
        !          &,iessw0,psip,ipsip,detmat_proj,nelorb_c,nelorb_at,nelcol_c&
        !          &,symmagp,yes_hermite)
        !!
        !    if(rank.eq.0) write(6,*) ' after  symmetrize_agp '
        !    if(rank.eq.0) write(6,*) ' Output detmat_proj '
        !!
        !       cost=0.d0
        !       do i=1,nelorb_at
        !        do j=i,nelorb_at
        !        if(rank.eq.0) write(6,*) i,j,detmat_proj(ipc*nelorb_c*(j-1)+2*i-1:ipc*nelorb_c*(j-1)+2*i),&
        !        &detmat_proj(ipc*nelorb_c*(i-1)+2*j-1:ipc*nelorb_c*(i-1)+2*j)
        !        cost=cost+sum(abs(detmat_proj(ipc*nelorb_c*(j-1)+2*i-1:ipc*nelorb_c*(j-1)+2*i)))+&
        !        &sum(abs(detmat_proj(ipc*nelorb_c*(i-1)+2*j-1:ipc*nelorb_c*(i-1)+2*j)))
        !        enddo
        !       enddo
        !       if(rank.eq.0) write(6,*) ' sum rule ',cost
        !
        !
        !    deallocate(detmat_proj)
        !
        !   call mpi_finalize(ierr)
        !   stop
        !!   END TEST
        if (allocated(ipsip)) deallocate (ipsip)
        allocate (ipsip(iscraipsip))
        ipsip = 0
        if (allocated(psip)) deallocate (psip)
        allocate (psip(iscramax))
        psip = 0.d0

    end subroutine convertmol_fast

    subroutine convertmol_c
        use allio, only: norm_metric
        implicit none
        real*8 ddot, dnrm2, overmax, cost, r0, psiln, rc(3)&
                &, jastrow_ei, timechange, eigref, LBox_sav
        real*8, external :: cclock
        real ran
        integer iseed, dimvect, icount, iicount, indfirst, ind_c, ind_t&
                &, indlast, ndiffdim_c, i, ix, iy, mine, maxe, nbufp, nmolu, shiftnelorb_at
        logical flagdo_down

#ifdef PARALLEL
        include 'mpif.h'
#endif
        !#ifdef __CASO
        !    nprocu=nprocopt
        !#else
        !    nprocu=1
        !#endif
        timechange = cclock()

        ! This subroutine changes detmat_c after an optimization step.
        if ((ipc .eq. 2 .and. (opposite_phase .or. same_phase) .and. real_contracted) .or. ipc .eq. 1) then
            flagdo_down = .false.
        else
            flagdo_down = .true.
        end if
        mesh = nx
        mesh = mesh*ny
        mesh = mesh*nz

        LBox_sav = LBox

        volmesh = ax*ay*az*unit_volume
        volmeshc = dcmplx(volmesh)

        !     the reference is the average ion position
        call shift_originref

        x = 0.d0

        !     nelorb_at=nelorb_c-molecular
        !     ndiff unchanged in input
        nmolu = nmol
        nmoldiff = nmolu + ndiff
        nmoltot = nmoldiff
        ndiffdim_c = ndiff*nelorb_at
        if ((.not. symmagp .or. ipc .eq. 2) .or. ipf .eq. 2) nmoltot = nmoltot + nmolu
        ! first part calculation overlap determinant

        if (ipc .eq. 2) then

            if (mesh .gt. nbufd) then
                allocate (buffer(2*nelorbh, 2*nbufd))
                allocate (buffer_c(2*nelorb_at, nbufd))
                nbufp = nbufd + 1
                nbufu = nbufd
            else
                allocate (buffer(2*nelorbh, 2*mesh))
                allocate (buffer_c(2*nelorb_at, mesh))
                nbufp = mesh + 1
                nbufu = mesh
            end if

        else

            if (mesh .gt. nbufd) then
                allocate (buffer(nelorbh, nbufd))
                allocate (buffer_c(nelorb_at, nbufd))
                nbufu = nbufd
                nbufp = nbufd + 1
            else
                allocate (buffer(nelorbh, mesh))
                allocate (buffer_c(nelorb_at, mesh))
                nbufp = mesh + 1
                nbufu = mesh
            end if

        end if

        if (.not. allocated(eigmol)) then
            allocate (eigmol(nelorb_at))
            eigmol = 0.d0
        end if

        if (ipf .eq. 2) then
            shiftnelorb_at = 0
        else
            shiftnelorb_at = nelorb_at
        end if

        ! allocate (oversl(nelorbh,nelorb_at))
        if (symmagp .and. ipc .eq. 1 .or. ipf .eq. 2) then
            allocate (overs(ipc*nelorb_at, nelorb_at))
            allocate (molecorb(ipf*ipc*nelorbh, nmoltot))
            allocate (molecorb_c(ipc*nelorb_at, nelorb_at))
        else
            allocate (overs(ipc*nelorb_at, 2*nelorb_at))
            allocate (molecorb(ipc*nelorbh, nmoltot))
            allocate (molecorb_c(ipc*nelorb_at, 2*nelorb_at))
        end if
        buffer = 0.d0
        buffer_c = 0.d0
        eigmol = 0.d0
        molecorb = 0.d0

        overs = 0.d0
        iflagnorm = 3
        indr = 0
        ind = 0
        nbuf = 0
#ifdef _OFFLOAD
! In any case take what is defined in the CPU for mu_c
!$omp target update to (mu_c)
!$omp target data map(overs) map(alloc:buffer,buffer_c)
#endif
        do k = 1, nz
            !          x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
                !             x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
                    indr = indr + 1

                    if (indr - (indr/nprocopt)*nprocopt .eq. rankopt) then
                        ind = ind + 1
                        !                   x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)

                        x(:) = rion_ref(:) + (-(nz + 1)/2.d0 + k)*az*at(:, 3) + &
                                &(-(ny + 1)/2.d0 + j)*ay*at(:, 2) + (-(nx + 1)/2.d0 + i)*ax*at(:, 1)

                        call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu &
                                &, dupr, zetar, rion, psip, buffer(1, ind), nelorbh, nion, kion &
                                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                &, indpar_tab, indorb_tab, indshell_tab, .true.)

                        if (flagdo_down) then
                            call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu &
                                    &, dupr, zetar, rion, psip, buffer(1, ind + nbufu), nelorbh, nion, kion &
                                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                    &, indpar_tab, indorb_tab, indshell_tab, .false.)
                        end if

                        !  DEPRECATED. It is overflow unstable for small vj and the convergence in the mesh is
                        !  slower due to implicit cusps in the orbitals added by the one body Jastrow.
                        if (add_onebody2det) then
                            psiln = -scale_one_body
                            do jj = 1, nion
                                if (iespbc) then
                                    rc(:) = x(:) - rion(:, jj)
                                    call CartesianToCrystal(rc, 1)
                                    do kk = 1, 3
                                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                                    end do
                                    r0 = norm_metric(rc, metric)
                                else
                                    rc(:) = (x(:) - rion(:, jj))*costz(jj)
                                    r0 = dsqrt(sum(rc(:)**2))
                                end if
                                psiln = psiln - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                            end do
                            buffer(1:ipc*nelorbh, ind) = buffer(1:ipc*nelorbh, ind)*dexp(psiln)
                            if (flagdo_down) &
                                    &buffer(1:2*nelorbh, ind + nbufu) = buffer(1:2*nelorbh, ind + nbufu)*dexp(psiln)
                        end if

                        if (mod(ind, nbufu) .eq. 0) then
#ifdef _OFFLOAD
!$omp target update to (buffer)
#endif

                            if (ipc .eq. 1) then
                                call dgemm_('T', 'N', nelorb_at, nbufu, nelorbh, 1.d0, mu_c, ipf*nelorbh, buffer&
                                        &, nelorbh, 0.d0, buffer_c, nelorb_at)
                                call dgemm_('N', 'T', nelorb_at, nelorb_at, nbufu, volmesh, buffer_c&
                                        &, nelorb_at, buffer_c, nelorb_at, 1.d0, overs, nelorb_at)
                                if (ipf .eq. 2) then

                                    call dgemm_('T', 'N', nelorb_at, nbufu, nelorbh, 1.d0, mu_c(nelorbh + 1, 1)&
                                            &, ipf*nelorbh, buffer, nelorbh, 0.d0, buffer_c, nelorb_at)

                                    call dgemm_('N', 'T', nelorb_at, nelorb_at, nbufu, volmesh, buffer_c&
                                            &, nelorb_at, buffer_c, nelorb_at, 1.d0, overs, nelorb_at)

                                end if

                            else
                                call zgemm_('T', 'N', nelorb_at, nbufu, nelorbh, zone, mu_c, ipf*nelorbh, buffer&
                                        &, nelorbh, zzero, buffer_c, nelorb_at)
                                call zgemm_('N', 'C', nelorb_at, nelorb_at, nbufu, volmeshc, buffer_c&
                                        &, nelorb_at, buffer_c, nelorb_at, zone, overs, nelorb_at)
                                if (flagdo_down) then

                                    call zgemm_('T', 'N', nelorb_at, nbufu, nelorbh, zone&
                                            &, mu_c(2*(ipf - 1)*nelorbh + 1, 1), ipf*nelorbh, buffer(1, nbufp)&
                                            &, nelorbh, zzero, buffer_c, nelorb_at)
                                    call zgemm_('N', 'C', nelorb_at, nelorb_at, nbufu, volmeshc, buffer_c&
                                            &, nelorb_at, buffer_c, nelorb_at, zone, overs(1, shiftnelorb_at + 1), nelorb_at)
                                elseif (ipf .eq. 2) then
                                    if (opposite_phase) call conjmat_(nelorbh, nbufu, buffer, nelorbh)
                                    call zgemm_('T', 'N', nelorb_at, nbufu, nelorbh, zone&
                                            &, mu_c(2*nelorbh + 1, 1), 2*nelorbh, buffer, nelorbh, zzero, buffer_c, nelorb_at)
                                    if (opposite_phase) call conjmat_(nelorbh, nbufu, buffer, nelorbh)
                                    call zgemm_('N', 'C', nelorb_at, nelorb_at, nbufu, volmeshc, buffer_c&
                                            &, nelorb_at, buffer_c, nelorb_at, zone, overs, nelorb_at)
                                end if
                            end if
                            !            call dgemm('N','T',nelorbh,nelorb_at,nbufd,volmesh,buffer&
                            !    &,nelorbh,buffer_c,nelorb_at,1.d0,oversl,nelorbh)
                            nbuf = nbuf + 1
                            ind = 0
                        end if
                    end if
                end do
            end do
        end do

        if (mod(ind, nbufu) .ne. 0) then
#ifdef _OFFLOAD
!$omp target update to (buffer)
#endif
            nleft = ind
            if (ipc .eq. 1) then
                call dgemm_('T', 'N', nelorb_at, nleft, nelorbh, 1.d0, mu_c, ipf*nelorbh, buffer&
                        &, nelorbh, 0.d0, buffer_c, nelorb_at)
                call dgemm_('N', 'T', nelorb_at, nelorb_at, nleft, volmesh, buffer_c&
                        &, nelorb_at, buffer_c, nelorb_at, 1.d0, overs, nelorb_at)
                if (ipf .eq. 2) then
                    call dgemm_('T', 'N', nelorb_at, nleft, nelorbh, 1.d0, mu_c(nelorbh + 1, 1)&
                            &, 2*nelorbh, buffer, nelorbh, 0.d0, buffer_c, nelorb_at)
                    call dgemm_('N', 'T', nelorb_at, nelorb_at, nleft, volmesh, buffer_c&
                            &, nelorb_at, buffer_c, nelorb_at, 1.d0, overs, nelorb_at)
                end if
            else
                call zgemm_('T', 'N', nelorb_at, nleft, nelorbh, zone, mu_c, ipf*nelorbh, buffer&
                        &, nelorbh, zzero, buffer_c, nelorb_at)
                call zgemm_('N', 'C', nelorb_at, nelorb_at, nleft, volmeshc, buffer_c&
                        &, nelorb_at, buffer_c, nelorb_at, zone, overs, nelorb_at)
                if (flagdo_down) then
                    call zgemm_('T', 'N', nelorb_at, nleft, nelorbh, zone&
                            &, mu_c(2*(ipf - 1)*nelorbh + 1, 1), ipf*nelorbh, buffer(1, nbufp)&
                            &, nelorbh, zzero, buffer_c, nelorb_at)
                    call zgemm_('N', 'C', nelorb_at, nelorb_at, nleft, volmeshc, buffer_c&
                            &, nelorb_at, buffer_c, nelorb_at, zone, overs(1, shiftnelorb_at + 1), nelorb_at)
                elseif (ipf .eq. 2) then
                    if (opposite_phase) call conjmat_(nelorbh, nleft, buffer, nelorbh)
                    call zgemm_('T', 'N', nelorb_at, nleft, nelorbh, zone&
                            &, mu_c(2*nelorbh + 1, 1), 2*nelorbh, buffer, nelorbh, zzero, buffer_c, nelorb_at)
                    if (opposite_phase) call conjmat_(nelorbh, nleft, buffer, nelorbh)
                    call zgemm_('N', 'C', nelorb_at, nelorb_at, nleft, volmeshc, buffer_c&
                            &, nelorb_at, buffer_c, nelorb_at, zone, overs, nelorb_at)
                end if
            end if
            !            call dgemm('N','T',nelorbh,nelorb_at,nleft,volmesh,buffer   &
            !    &,nelorbh,buffer_c,nelorb_at,1.d0,oversl,nelorbh)
            nbuf = nbuf + 1
        end if
#ifdef _OFFLOAD
!$omp end target data
#endif

#ifdef PARALLEL
        !     nelorb2=nelorbh*nelorb_at
        !     CALL reduce_base_real(nelorb2,oversl,commopt_mpi,-1)
        if (ipc .eq. 2 .and. ipf .eq. 1) then
            nelorb2 = 4*nelorb_at*nelorb_at
        else
            nelorb2 = ipc*nelorb_at*nelorb_at
        end if

        call reduce_base_real(nelorb2, overs, commopt_mpi, -1)

#endif
        if (ipc .eq. 2 .and. ipf .eq. 1) then
            if (.not. flagdo_down) then
                overs(1:2*nelorb_at, nelorb_at + 1:2*nelorb_at) = overs(1:2*nelorb_at, 1:nelorb_at)
                if (opposite_phase) call conjmat(nelorb_at, nelorb_at, overs(1, nelorb_at + 1), nelorb_at)
            end if
            call conjmat(nelorb_at, 2*nelorb_at, overs, nelorb_at) ! The complex conjugate in standard metric
        end if

        !   NB for ipc=2 overlap for the down is explicitly computed and generally different
        if (.not. symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) overs(1:ipc*nelorb_at, nelorb_at + 1:2*nelorb_at) = &
                &overs(1:ipc*nelorb_at, 1:nelorb_at)

        if (yesmin .ne. 0 .or. read_molecul) then
            do i = 1, nelorb_at
                psip(ipc*nelorb_at*(i - 1) + 1:ipc*nelorb_at*i) = &
                        &detmat_proj(ipc*nelorb_c*(i - 1) + 1:ipc*nelorb_c*(i - 1) + ipc*nelorb_at)&
                        & + detmat_c(ipc*nelorb_c*(i - 1) + 1:ipc*nelorb_c*(i - 1) + ipc*nelorb_at)
            end do
            do i = 1, ndiff
                psip(ipc*nelorb_at*(nelorb_at + i - 1) + 1:ipc*nelorb_at*(nelorb_at + i)) = &
                        detmat_proj(ipc*nelorb_c*(nelorb_c + i - 1) + 1:ipc*nelorb_c*(nelorb_c + i - 1) + ipc*nelorb_at) + &
                                &detmat_c(ipc*nelorb_c*(nelorb_c + i - 1) + 1:ipc*nelorb_c*(nelorb_c + i - 1) + ipc*nelorb_at)
            end do
            detmat_c(1:ipc*nelorb_at*(nelorb_at + ndiff)) = psip(1:ipc*nelorb_at*(nelorb_at + ndiff))

            if (.not. allocated(projmat_c)) then
                if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                    allocate (projmat_c(ipc*nelorb_at, 2*nmolmat))
                else
                    allocate (projmat_c(ipc*nelorb_at, 4*nmolmat))
                end if
                projmat_c = 0.d0
            else
                call projectmat(nelorb_at, nelorb_at, nmolmatdo&
                        &, nmolmat, projmat_c, nmolmat, detmat_c, psip, symmagp, yesmin)
            end if

            if (symmetrize_agp) then

                psip(1:ipc*nelorb_at*(nelorb_at + ndiff)) = detmat_c(1:ipc*nelorb_at*(nelorb_at + ndiff))
                detmat_c = 0.d0
                do i = 1, nelorb_at
                    detmat_c(ipc*nelorb_c*(i - 1) + 1:ipc*nelorb_c*(i - 1) + ipc*nelorb_at) = &
                        psip(ipc*nelorb_at*(i - 1) + 1:ipc*nelorb_at*i)
                end do
                do i = 1, ndiff
                    detmat_c(ipc*nelorb_c*(nelorb_c + i - 1) + 1:ipc*nelorb_c*(nelorb_c + i - 1) + ipc*nelorb_at) = &
                            &psip(ipc*nelorb_at*(nelorb_at + i - 1) + 1:ipc*nelorb_at*(nelorb_at + i))
                end do

                if (rankrep .eq. 0) then
                    !  write(6,*) ' before  symmetrize agp '
                    cost = 0.d0
                    do i = 1, nelorb_at
                        do j = 1, nelorb_at
                            !     write(6,*) i,j,detmat_c(2*nelorb_c*(j-1)+2*i-1),detmat_c(2*nelorb_c*(i-1)+2*j-1)
                            cost = cost + sum(abs(detmat_c(ipc*nelorb_c*(j - 1) + ipc*(i - 1) + 1:&
                                    &ipc*nelorb_c*(j - 1) + ipc*i)))
                        end do
                    end do
                end if
                if (kaverage) then
#ifdef  PARALLEL
                    allocate (agpo_kp(nk), agp_kp(nk))
                    if (rankrep .eq. 0) call mpi_gather(cost, 1, MPI_DOUBLE_PRECISION, agpo_kp,&
                   &1, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
#endif
                else
                    if (rank .eq. 0) write (6, *) ' before  symmetrize agp sum rule= ', cost
                end if

                !   dsw can be used as is recomputed outside
                ! it is done before changing detmat_proj because we want that detmat_c and
                ! detmat_proj are exactly consistent after the convertmol_c.
                ! In this way we accept a small roundoff non symmetric component given by the
                ! diagonalization. It is clear that if the number of molecular orbitals is not
                ! consistent with the symmetry chosen, there may be some problem to satisfy symmetries.
                !          From real to effective
                if (allowed_averagek) call attach_phase2det(.false., detmat_c)

                call symmetrizeagp(nnozero_c, nozero_c, jbradet, jbradetn, dsw&
                        &, iessw0, psip, ipsip, detmat_c, nelorb_c, nelorb_at, nelcol_c&
                        &, symmagp, yes_hermite)

                !          From effective to real
                if (allowed_averagek) call attach_phase2det(.true., detmat_c)
                if (rankrep .eq. 0) then
                    !   write(6,*) ' after symmetrize agp rank  ='
                    cost = 0.d0
                    do i = 1, nelorb_at
                        do j = 1, nelorb_at
                            !   write(6,*) i,j,detmat_c(2*nelorb_c*(j-1)+2*i-1),detmat_c(2*nelorb_c*(i-1)+2*j-1)
                            cost = cost + sum(abs(detmat_c(ipc*nelorb_c*(j - 1) + ipc*(i - 1) + 1:&
                                    &ipc*nelorb_c*(j - 1) + ipc*i)))
                        end do
                    end do
                    !  write(6,*) ' new sum rule =',cost
                end if

                if (kaverage) then
#ifdef  PARALLEL
                    if (rankrep .eq. 0) call mpi_gather(cost, 1, MPI_DOUBLE_PRECISION, agp_kp,&
                    &1, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
                    if (rank .eq. 0) then
                    do ii = 1, nk
                    if (abs((agpo_kp(ii) - agp_kp(ii))/agpo_kp(ii)) .gt. 1d-6) then
                        write (6, *) ' ERROR agp symmetrize in momentum K= before/after'&
                        &, ii, agpo_kp(ii), agp_kp(ii)
                    end if
                    end do
                    write (6, *) ' after symmetrize agp relative error='&
                  &, sum(abs(agpo_kp(:) - agp_kp(:)))/sum(abs(agpo_kp(:)))
                    end if
                    deallocate (agp_kp, agpo_kp)
#endif
                else
                    if (rank .eq. 0) write (6, *) ' after  symmetrize agp sum rule= ', cost
                end if

                do i = 1, nelorb_at
                    psip(ipc*nelorb_at*(i - 1) + 1:ipc*nelorb_at*i) = &
                        detmat_c(ipc*nelorb_c*(i - 1) + 1:ipc*nelorb_c*(i - 1) + ipc*nelorb_at)
                end do
                do i = 1, ndiff
                    psip(ipc*nelorb_at*(nelorb_at + i - 1) + 1:ipc*nelorb_at*(nelorb_at + i)) = &
                        detmat_c(ipc*nelorb_c*(nelorb_c + i - 1) + 1:ipc*nelorb_c*(nelorb_c + i - 1) + ipc*nelorb_at)
                end do
                detmat_c(1:ipc*nelorb_at*(nelorb_at + ndiff)) = psip(1:ipc*nelorb_at*(nelorb_at + ndiff))
            end if ! endif symmetrize_agp
            ! save overlap matrix
            if (ipc .eq. 1 .and. symmagp .or. ipf .eq. 2) then
                allocate (over_sav(ipc*nelorb_at, nelorb_at))
                over_sav(1:ipc*nelorb_at, 1:nelorb_at) = overs(1:ipc*nelorb_at, 1:nelorb_at)
            else
                allocate (over_sav(ipc*nelorb_at, 2*nelorb_at))
                over_sav(1:ipc*nelorb_at, 1:2*nelorb_at) = overs(1:ipc*nelorb_at, 1:2*nelorb_at)
            end if
        end if ! endif yesmin (basically always)

        if (ndiff .eq. 0) then
            if (rank .eq. 0) write (6, *) ' Eigenvalue det 2 '
            call eval_molec_epsdgel(nelorb_at, overs, detmat_c, molecorb_c&
                    &, eigmol, nelorb_at, epsdgm, nprocu, rank, rankopt, commopt_mpi, 1, symmagp)

            overs = molecorb_c
            psip(1:nelorb_at) = eigmol(1:nelorb_at)
            if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
                do i = 1, nelorb_at
                    molecorb_c(1:ipc*nelorb_at, i) = overs(1:ipc*nelorb_at, nelorb_at - i + 1)
                end do
                do i = 1, nelorb_at/ipf
                    eigmol(i) = psip(nelorb_at/ipf - i + 1)
                end do
            else
                do i = 1, nelorb_at
                    molecorb_c(1:ipc*nelorb_at, 2*i - 1) = overs(1:ipc*nelorb_at, 2*nelorb_at - 2*i + 1)
                    molecorb_c(1:ipc*nelorb_at, 2*i) = overs(1:ipc*nelorb_at, 2*nelorb_at - 2*i + 2)
                    eigmol(i) = psip(nelorb_at - i + 1)
                end do
            end if
        else
            !  othorgonalize the AGP part with respect to the unpaired
            !  then orthogonalize all and fill properly the unpaired
            !       for test

            if (.not. allocated(over_sav)) then

                allocate (over_sav(ipc*nelorb_at, nelorb_at))
                over_sav(1:ipc*nelorb_at, 1:nelorb_at) = overs(1:ipc*nelorb_at, 1:nelorb_at)

            end if

            indunp = ipc*nelorb_at*nelorb_at + 1

            call eval_molec_unpaired(nelorb_at, overs, detmat_c, molecorb_c&
                    &, eigmol, nelorb_at, epsdgm, nprocu, rank, rankopt, commopt_mpi, ndiff, detmat_c(indunp)&
                    &, orthoyes, symmagp)
            psip(1:nelorb_at) = eigmol(1:nelorb_at)
            do i = 1, nelorb_at/ipf
                eigmol(i) = psip(nelorb_at/ipf - i + 1)
            end do
            if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
                do i = 1, nelorb_at
                    overs(1:ipc*nelorb_at, i) = molecorb_c(1:ipc*nelorb_at, nelorb_at - i + 1)
                end do
            else
                do i = 1, nelorb_at
                    overs(1:ipc*nelorb_at, i) = molecorb_c(1:ipc*nelorb_at, 2*nelorb_at - 2*i + 2)
                end do
                do i = 1, nelorb_at
                    overs(1:ipc*nelorb_at, nelorb_at + i) = molecorb_c(1:ipc*nelorb_at, 2*nelorb_at - 2*i + 1)
                end do
            end if

            !       putting in the same order
            if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
                do i = 1, nelorb_at
                    molecorb_c(1:ipc*nelorb_at, i) = overs(1:ipc*nelorb_at, i)
                end do
            else
                do i = 1, nelorb_at
                    molecorb_c(1:ipc*nelorb_at, 2*i) = overs(1:ipc*nelorb_at, i)
                end do
                do i = 1, nelorb_at
                    molecorb_c(1:ipc*nelorb_at, 2*i - 1) = overs(1:ipc*nelorb_at, nelorb_at + i)
                end do
            end if

            if (gramyes .and. nmoltot .eq. nelup) then
                !   Copy the unpaired in detmat_c in overs

                do i = 1, ndiff
                    call dcopy(ipc*nelorb_at, detmat_c(indunp + (i - 1)*(ipc*nelorb_at)), 1&
                            &, overs(1, nmoltot - ndiff + i), 1)
                end do

                if (rankopt .eq. 0) then
                    if (ipc .eq. 1) then
                        call graham(overs, over_sav, nelorb_at, psip, nelorb_at&
                                &, nelorb_at, nelup, info)
                    else
                        call graham_complex(overs, over_sav, nelorb_at, psip, nelorb_at&
                                &, nelorb_at, nelup, info)
                    end if
                end if
#ifdef PARALLEL
                call bcast_real(overs, ipc*nelorb_at*nelup, 0, commopt_mpi)
                call mpi_bcast(info, 1, MPI_INTEGER, 0, commopt_mpi, ierr)
#endif

                if (info .ne. 0) then
                    if (rank .eq. 0) write (6, *) ' Error dependency in unpaired orbitals ', info
#ifdef PARALLEL
                    call mpi_finalize(ierr)
#endif
                    stop
                end if
                !   Copy back detmat_c now orthogonal to the paired orbitals.
                do i = 1, ndiff
                    call dcopy(ipc*nelorb_at, overs(1, nmoltot - ndiff + i), 1&
                            &, detmat_c(indunp + (i - 1)*(ipc*nelorb_at)), 1)
                end do

            end if

        end if ! endif ndiff

        if (yesmin .ne. 0) then
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                call dgemm_my('N', 'N', nelorb_at, nmolmat, nelorb_at, 1.d0 &
                        &, over_sav, nelorb_at, molecorb_c, nelorb_at, 0.d0&
                        &, projmat_c(1, nmolmat + 1), nelorb_at, nprocu, rankopt, commopt_mpi)
                do ii = 1, nmolmat
                    projmat_c(1:nelorb_at, ii) = molecorb_c(1:nelorb_at, ii)
                end do
            else
                if (ipc .eq. 1) then
                    if (ipf .eq. 2) then
                        call dgemm_my('N', 'N', nelorb_at, nmolmat, nelorb_at, 1.d0&
                                &, over_sav, nelorb_at, molecorb_c, nelorb_at, 0.d0&
                                &, projmat_c(1, 3*nmolmat + 1), nelorb_at, nprocu, rankopt, commopt_mpi)
                        projmat_c(:, 2*nmolmat + 1:3*nmolmat) = projmat_c(:, 3*nmolmat + 1:4*nmolmat)
                        do ii = 1, nmolmat
                            projmat_c(1:nelorb_at, ii) = molecorb_c(1:nelorb_at, ii)
                            projmat_c(1:nelorb_at, ii + nmolmat) = molecorb_c(1:nelorb_at, ii)
                        end do
                    else
                        call dgemm_my('N', 'N', nelorb_at, 2*nmolmat, nelorb_at, 1.d0&
                                &, over_sav, nelorb_at, molecorb_c, nelorb_at, 0.d0&
                                &, projmat_c, nelorb_at, nprocu, rankopt, commopt_mpi)
                    end if
                else
                    !  The left and right have different overlap matrix
                    !            Left eigenvectors with over_sav NO ARROCCO here

                    if (ipf .eq. 2) then
                        call zgemm_my('N', 'N', nelorb_at, nmolmat, nelorb_at, zone&
                                &, over_sav, nelorb_at, molecorb_c, nelorb_at, zzero&
                                &, projmat_c(1, 3*nmolmat + 1), nelorb_at, nprocu, rankopt, commopt_mpi)
                        projmat_c(:, 2*nmolmat + 1:3*nmolmat) = projmat_c(:, 3*nmolmat + 1:4*nmolmat)
                        do ii = 1, nmolmat
                            projmat_c(1:ipc*nelorb_at, ii) = molecorb_c(1:ipc*nelorb_at, ii)
                            projmat_c(1:ipc*nelorb_at, ii + nmolmat) = molecorb_c(1:ipc*nelorb_at, ii)
                        end do
                        call conjmat(nelorb_at, nmolmat, projmat_c(1, 2*nmolmat + 1), nelorb_at)
                        call conjmat(nelorb_at, nmolmat, projmat_c, nelorb_at)

                    else

                        call zgemm_my('N', 'N', nelorb_at, nmolmat, nelorb_at, zone&
                                &, over_sav, nelorb_at, molecorb_c(1, 2), 2*nelorb_at, zzero&
                                &, projmat_c(1, 2), 2*nelorb_at, nprocu, rankopt, commopt_mpi)

                        !  The right eigenvectors are the complex conjugate of the output
                        call conjmat(nelorb_at, nmolmat, molecorb_c, 2*nelorb_at)
                        call conjmat(nelorb_at, nelorb_at, over_sav(1, nelorb_at + 1), nelorb_at)
                        call zgemm_my('N', 'N', nelorb_at, nmolmat, nelorb_at, zone&
                                &, over_sav(1, nelorb_at + 1), nelorb_at, molecorb_c, 2*nelorb_at, zzero&
                                &, projmat_c, 2*nelorb_at, nprocu, rankopt, commopt_mpi)
                        ! Restoring the output molecorb
                        call conjmat(nelorb_at, nmolmat, molecorb_c, 2*nelorb_at)
                        call conjmat(nelorb_at, nelorb_at, over_sav(1, nelorb_at + 1), nelorb_at)
                    end if ! end ipf
                end if ! endif ipc

                if (ipf .ne. 2) then
                    do ii = 1, nmolmat
                        projmat_c(1:ipc*nelorb_at, 2*nmolmat + ii) = projmat_c(1:ipc*nelorb_at, 2*ii - 1)
                        projmat_c(1:ipc*nelorb_at, 3*nmolmat + ii) = projmat_c(1:ipc*nelorb_at, 2*ii)
                    end do
                    !  The right eigenvectors are the complex conjugate of the output
                    if (ipc .eq. 2) call conjmat(nelorb_at, nmolmat, molecorb_c, 2*nelorb_at)
                    do ii = 1, nmolmat
                        projmat_c(1:ipc*nelorb_at, ii) = molecorb_c(1:ipc*nelorb_at, 2*ii - 1)
                        projmat_c(1:ipc*nelorb_at, ii + nmolmat) = molecorb_c(1:ipc*nelorb_at, 2*ii)
                    end do
                    !   Restoring the output
                    if (ipc .eq. 2) call conjmat(nelorb_at, nmolmat, molecorb_c, 2*nelorb_at)
                    !          if(ipc.eq.2.and.symmagp.and.same_phase.and.ndiff.gt.0) then
                    !           do ii=1,nmolmat
                    !           projmat_c(1:2*nelorb_at,ii)=projmat_c(1:2*nelorb_at,ii+nmolmat)
                    !           projmat_c(1:2*nelorb_at,2*nmolmat+ii)=projmat_c(1:2*nelorb_at,ii+3*nmolmat)
                    !           enddo
                    !          endif
                    !          if(ipc.eq.2.and.symmagp.and.opposite_phase.and.ndiff.gt.0) then
                    !           do ii=1,nmolmat
                    !           projmat_c(1:2*nelorb_at:2,ii)=projmat_c(1:2*nelorb_at:2,ii+nmolmat)
                    !           projmat_c(2:2*nelorb_at:2,ii)=-projmat_c(2:2*nelorb_at:2,ii+nmolmat)
                    !           projmat_c(1:2*nelorb_at:2,2*nmolmat+ii)=projmat_c(1:2*nelorb_at:2,ii+3*nmolmat)
                    !           projmat_c(2:2*nelorb_at:2,2*nmolmat+ii)=-projmat_c(2:2*nelorb_at:2,ii+3*nmolmat)
                    !           enddo
                    !          endif
                end if ! if ipf =/2
            end if ! symmagp.and.ipc.eq.1.and.ipf.eq.1
            if (allocated(over_sav)) deallocate (over_sav)
        end if ! endif yesmin

        do i = 1, ndiff
            call dcopy(ipc*nelorb_at, detmat_c(indunp + (i - 1)*(ipc*nelorb_at)), 1&
                    &, molecorb_c(1, nmoltot - ndiff + i), 1)
        end do

! the matrix-matrix below are huge as are prop. to the scontracted basis nelorb
! Thus  I put it in the GPU
        if (ipc .eq. 1) then
#ifdef  _OFFLOAD
!$omp target data map(from:molecorb) map(to:molecorb_c)
            call dgemm_('N', 'N', ipf*nelorbh, nmoltot, nelorb_at, 1.d0&
   &, mu_c, ipf*nelorbh, molecorb_c, nelorb_at, 0.d0, molecorb, ipf*nelorbh)
!$omp end target data
#else
            call dgemm_my('N', 'N', ipf*nelorbh, nmoltot, nelorb_at, 1.d0, mu_c, ipf*nelorbh&
                    &, molecorb_c, nelorb_at, 0.d0, molecorb, ipf*nelorbh, nprocu, rankopt, commopt_mpi)
#endif
        else
#ifdef  _OFFLOAD
!$omp target data map(from:molecorb) map(to:molecorb_c)
            call zgemm_('N', 'N', ipf*nelorbh, nmoltot, nelorb_at, zone&
   &, mu_c, ipf*nelorbh, molecorb_c, nelorb_at, zzero, molecorb, ipf*nelorbh)
!$omp end target data
#else
            call zgemm_my('N', 'N', ipf*nelorbh, nmoltot, nelorb_at, zone, mu_c, ipf*nelorbh&
                    &, molecorb_c, nelorb_at, zzero, molecorb, ipf*nelorbh, nprocu, rankopt, commopt_mpi)
#endif
        end if

        indpar = iesup_c - 2*nmoltot*nelorbh*ipf
        !   if(ipc.eq.1) then
        do i = 1, nmoltot
            do j = 1, ipc*nelorbh*ipf, ipc
                dup_c(ipc*indpar + j) = (j - 1)/ipc + 1
                if (ipc .eq. 2) dup_c(2*indpar + j + 1) = 0.d0
            end do
            do j = nelorbh*ipf*ipc + 1, 2*ipf*ipc*nelorbh
                dup_c(indpar*ipc + j) = molecorb(j - ipf*ipc*nelorbh, i)
            end do
            indpar = indpar + 2*nelorbh*ipf
        end do

        ind = nnozero_c - nmol - ndiff
        detmat_c = 0

        do i = 1, nmol
            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
                ix = nelorb_c - nmoltot + i
                iy = ix
            else
                ix = nelorb_c - nmoltot + 2*i
                iy = ix - 1
            end if
            ind = ind + 1
            ind_c = nozero_c(ind)
            ind_t = nelorb_c*(ix - 1) + iy
            if (i .le. nmolmin) then
                detmat_c(ipc*(nozero_c(ind) - 1) + 1) = dble(nmolmin - i + 1)
                if (ipc .eq. 2) detmat_c(2*nozero_c(ind)) = 0.d0
            elseif (i .le. nmolmaxw) then
                if (power .eq. 1.d0) then
                    detmat_c(ipc*(nozero_c(ind) - 1) + 1) = eigmol(i)/eigmol(nmolmin)
                    if (ipc .eq. 2) detmat_c(2*nozero_c(ind)) = 0.d0
                else
                    cost = eigmol(i)/eigmol(nmolmin)
                    detmat_c(ipc*(nozero_c(ind) - 1) + 1) = abs(cost)**power*sign(1.d0, cost)
                    if (ipc .eq. 2) detmat_c(2*nozero_c(ind)) = 0.d0
                end if
            else
                detmat_c(ipc*(nozero_c(ind) - 1) + 1) = 0.d0
                if (ipc .eq. 2) detmat_c(2*nozero_c(ind)) = 0.d0
            end if
            psip(i) = detmat_c(ipc*(nozero_c(ind) - 1) + 1)
            !      Fullfill antisymmetry
            if (ipf .eq. 2) then
                if (ipc .eq. 1) then
                    detmat_c(ind_t) = -detmat_c(ind_c)
                else
                    detmat_c(2*ind_t - 1) = -detmat_c(2*ind_c - 1)
                    detmat_c(2*ind_t) = -detmat_c(2*ind_c)
                end if
            end if
        end do

        !       Compute detmat_proj
        detmat_proj = 0.d0
        if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then
            do i = 1, nmolu
                do j = 1, nelorb_at
                    molecorb(j, i) = molecorb_c(j, i)*psip(i)
                end do
            end do
            call dgemm_my('N', 'T', nelorb_at, nelorb_at, nmolu, 1.d0, molecorb_c, nelorb_at&
                    &, molecorb, nelorbh, 0.d0, detmat_proj, nelorb_c, nprocu, rankopt, commopt_mpi)
        else
            do i = 1, nmolu
                do j = 1, ipc*nelorb_at
                    molecorb(j, i) = molecorb_c(j, 2*i)*psip(i)
                end do
            end do
            if (ipc .eq. 2) then
                !      Here I am consistent with the definition without the complex conjugate
                call zgemm_my('N', 'T', nelorb_at, nelorb_at, nmolu, zone, molecorb, ipf*nelorbh&
                        &, molecorb_c, 2*nelorb_at, zzero, detmat_proj, nelorb_c, nprocu, rankopt, commopt_mpi)
            else
                call dgemm_my('N', 'T', nelorb_at, nelorb_at, nmolu, 1.d0, molecorb, ipf*nelorbh&
                        &, molecorb_c, 2*nelorb_at, 0.d0, detmat_proj, nelorb_c, nprocu, rankopt, commopt_mpi)
            end if
            if (ipf .eq. 2) then
                ! Antisymmetrize the pfaffian
                do i = 1, nelorb_at
                    detmat_proj(nelorb_c*ipc*(i - 1) + ipc*(i - 1) + 1:nelorb_c*ipc*(i - 1) + ipc*i) = 0.d0
                    if (ipc .eq. 1) then
                        do j = i + 1, nelorb_at
                            cost = detmat_proj(i + (j - 1)*nelorb_c) - detmat_proj(j + (i - 1)*nelorb_c)
                            detmat_proj(i + (j - 1)*nelorb_c) = cost
                            detmat_proj(j + (i - 1)*nelorb_c) = -cost
                        end do
                    else
                        do j = i + 1, nelorb_at
                            cost = detmat_proj(2*i - 1 + (j - 1)*2*nelorb_c) - detmat_proj(2*j - 1 + (i - 1)*2*nelorb_c)
                            detmat_proj(2*i - 1 + (j - 1)*2*nelorb_c) = cost
                            detmat_proj(2*j - 1 + (i - 1)*2*nelorb_c) = -cost
                            cost = detmat_proj(2*i + (j - 1)*2*nelorb_c) - detmat_proj(2*j + (i - 1)*2*nelorb_c)
                            detmat_proj(2*i + (j - 1)*2*nelorb_c) = cost
                            detmat_proj(2*j + (i - 1)*2*nelorb_c) = -cost
                        end do
                    end if
                end do
            end if
        end if

        do i = 1, ndiff
            detmat_proj(ipc*nelorb_c*(nelorb_c + i - 1) + 1:&
                    &ipc*nelorb_c*(nelorb_c + i - 1) + ipc*nelorb_at) = &
                    &molecorb_c(1:ipc*nelorb_at, nmoltot - ndiff + i)
        end do

        do i = 1, ndiff
            ind = ind + 1
            j = nelorb_c - ndiff + i
            detmat_c(ipc*(nozero_c(ind) - 1) + 1) = 1.d0
            if (ipc .eq. 2) detmat_c(2*nozero_c(ind)) = 0.d0
        end do

        deallocate (molecorb, molecorb_c, overs)
        deallocate (eigmol)

        if (ipc .eq. 1) then
            do jj = 1, iesup_c
                do kk = 1, multranspip(jj)
                    iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                    ix = transpip(kk)%col(jj) - (iy - 1)*nelorbh*ipf
                    mu_c(ix, iy) = dup_c(jj)
                end do
            end do
        else
            do jj = 1, iesup_c
                do kk = 1, multranspip(jj)
                    iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                    ix = transpip(kk)%col(jj) - (iy - 1)*nelorbh*ipf
                    mu_c(2*ix - 1:2*ix, iy) = dup_c(2*jj - 1:2*jj)
                end do
            end do
        end if

        if (allocated(buffer)) deallocate (buffer)
        if (allocated(buffer_c)) deallocate (buffer_c)

        if (rank .eq. 0) write (6, *) ' Time change fort.10 =', cclock() - timechange
        LBox = LBox_sav

    end subroutine convertmol_c

    subroutine ortho_fast
        implicit none
        real*8 ddot, dnrm2, overmax, cost, r0, psiln, rc(3)
        real ran
        integer iseed, dimvect, icount, iicount, indfirst, indlast, maxdimmol, dimpsip&
                &, nbufp
        real*8, dimension(:), allocatable :: psip_loc
#ifdef PARALLEL
        include 'mpif.h'
#endif
        !#ifdef __CASO
        !    nprocu=nprocopt
        !#else
        !    nprocu=1
        !#endif

        mesh = nx
        mesh = mesh*ny
        mesh = mesh*nz

        volmesh = ax*ay*az*unit_volume
        volmeshc = dcmplx(volmesh)
        !     the reference is the average ion position

        call shift_originref

        if (.not. iespbc .and. rank .eq. 0) write (6, *) ' Updated Center of mesh =', rion_ref(:)

        x = 0.d0

        !       first part calculation overlap determinant
        ! just the minimum allocation

        nelorb_at = nelorb_c - molecular

        nelorb_diag = nelorb_c
        nelorb2 = nelorb_diag*nelorb_diag

        nbufu = min(nbufd, nelorb_diag)
        if (mesh .lt. nbufu) nbufu = mesh

        if (ipc .eq. 1) then
            allocate (buffer(nelorbh, nbufu))
            allocate (overs(nelorb_diag, nelorb_diag))
        else
            allocate (buffer(2*nelorbh, 2*nbufu))
            allocate (overs(2*nelorb_diag, 2*nelorb_diag))
        end if
        nbufp = nbufu + 1

        overs = 0.d0
        maxdimmol = max(nelorb2, nelorb_diag*nbufu)
        maxdimmol = maxdimmol/nelorbh + 1
        allocate (molecorb(ipc*nelorbh, maxdimmol))
        molecorb = 0.d0
        dimpsip = max((indt + 5)*10, nummol)

        allocate (psip_loc(ipc*dimpsip))

        buffer = 0.d0
        molecorb = 0.d0
        overs = 0.d0

        overs = 0.d0

        iflagnorm = 3
        indr = 0
        ind = 0
        nbuf = 0
        do k = 1, nz
            !      x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
                !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
                    indr = indr + 1

                    if (indr - (indr/nprocopt)*nprocopt .eq. rankopt) then
                        ind = ind + 1
                        !               x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)
                        x(:) = rion_ref(:) + (-(nz + 1)/2.d0 + k)*az*at(:, 3) + &
                                &(-(ny + 1)/2.d0 + j)*ay*at(:, 2) + (-(nx + 1)/2.d0 + i)*ax*at(:, 1)

                        call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                &, dupr, zetar, rion, psip_loc, buffer(1, ind), nelorbh, nion, kion             &
                                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                &, indpar_tab, indorb_tab, indshell_tab, .true.)
                        if (ipc .eq. 2) then
                            call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                    &, dupr, zetar, rion, psip_loc, buffer(1, ind + nbufu), nelorbh, nion, kion             &
                                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                    &, indpar_tab, indorb_tab, indshell_tab, .false.)
                        end if

                        if (mod(ind, nbufu) .eq. 0) then
                            if (ipc .eq. 1) then
                                call dgemm('T', 'N', nelorb_diag, nbufu, nelorbh, 1.d0, mu_c, nelorbh&
                                        &, buffer, nelorbh, 0.d0, molecorb, nelorb_diag)
                                call dgemm('N', 'T', nelorb_diag, nelorb_diag, nbufu, volmesh, molecorb&
                                        &, nelorb_diag, molecorb, nelorb_diag, 1.d0, overs, nelorb_diag)
                            else
                                call zgemm('T', 'N', nelorb_diag, nbufu, nelorbh, zone, mu_c, nelorbh&
                                        &, buffer, nelorbh, zzero, molecorb, nelorb_diag)
                                call zgemm('N', 'C', nelorb_diag, nelorb_diag, nbufu, volmeshc, molecorb&
                                        &, nelorb_diag, molecorb, nelorb_diag, zone, overs, nelorb_diag)
                                call zgemm('T', 'N', nelorb_diag, nbufu, nelorbh, zone, mu_c, nelorbh&
                                        &, buffer(1, nbufp), nelorbh, zzero, molecorb, nelorb_diag)
                                call zgemm('N', 'C', nelorb_diag, nelorb_diag, nbufu, volmeshc, molecorb&
                                        &, nelorb_diag, molecorb, nelorb_diag, zone, overs(1, nelorb_diag + 1), nelorb_diag)
                            end if

                            nbuf = nbuf + 1
                            ind = 0
                        end if
                    end if
                end do
            end do
        end do

        if (mod(ind, nbufu) .ne. 0) then
            nleft = ind
            if (ipc .eq. 1) then
                call dgemm('T', 'N', nelorb_diag, nleft, nelorbh, 1.d0, mu_c, nelorbh&
                        &, buffer, nelorbh, 0.d0, molecorb, nelorb_diag)
                call dgemm('N', 'T', nelorb_diag, nelorb_diag, nleft, volmesh, molecorb&
                        &, nelorb_diag, molecorb, nelorb_diag, 1.d0, overs, nelorb_diag)
            else
                call zgemm('T', 'N', nelorb_diag, nleft, nelorbh, zone, mu_c, nelorbh&
                        &, buffer, nelorbh, zzero, molecorb, nelorb_diag)
                call zgemm('N', 'C', nelorb_diag, nelorb_diag, nleft, volmeshc, molecorb&
                        &, nelorb_diag, molecorb, nelorb_diag, zone, overs, nelorb_diag)
                call zgemm('T', 'N', nelorb_diag, nleft, nelorbh, zone, mu_c, nelorbh&
                        &, buffer(1, nbufp), nelorbh, zzero, molecorb, nelorb_diag)
                call zgemm('N', 'C', nelorb_diag, nelorb_diag, nleft, volmeshc, molecorb&
                        &, nelorb_diag, molecorb, nelorb_diag, zone, overs(1, nelorb_diag + 1), nelorb_diag)
            end if

            nbuf = nbuf + 1
        end if

#ifdef PARALLEL
        call reduce_base_real(size(overs), overs, commopt_mpi, -1)
#endif

        if (ipc .eq. 2) call conjmat(nelorb_diag, 2*nelorb_diag, overs, nelorb_diag) ! The complex conjugate in standard metric

        !       Input wheremol(1:nummol),nummol

        if (rank .eq. 0) then
            write (6, *) ' Overlap matrix between molecular orbitals  '
            do i = 1, nummol
                do j = i, nummol
                    write (6, *) i, wheremol(i), j, wheremol(j), overs(ipc*(wheremol(i) - 1) + 1:&
                            &ipc*wheremol(i), wheremol(j))
                end do
            end do
        end if

        ntotpsi = nelorb_c - nelorb_at

        allocate (psi_sav(ipc*nelorb_c, ntotpsi))

        psi_sav = 0.d0
        do i = 1, nummol
            psi_sav(ipc*(wheremol(i) - 1) + 1, i) = 1.d0
        end do
        if (ipc .eq. 1) then
            call graham(psi_sav(nelorb_at + 1, 1), overs(nelorb_at + 1, nelorb_at + 1), nelorb_c&
                    &, psip_loc, nelorb_c, ntotpsi, nummol, info)
        else
            call graham_complex(psi_sav(2*nelorb_at + 1, 1), overs(2*nelorb_at + 1, nelorb_at + 1), nelorb_c&
                    &, psip_loc, nelorb_c, ntotpsi, nummol, info)
        end if

        overs(:, 1:nummol) = psi_sav(:, 1:nummol)

        psi_sav = 0.d0

        do ii = nelorb_at + 1, nelorb_c
            psi_sav(ipc*(ii - 1) + 1, ii - nelorb_at) = 1.d0
        end do
        !        replace the orthogonalized molecular orbitals
        do i = 1, nummol
            ii = wheremol(i) - nelorb_at
            psi_sav(:, ii) = overs(:, i)
        end do

        overs = 0.d0
        do ii = 1, ntotpsi
            overs(1:ipc*nelorb_c, ii) = psi_sav(1:ipc*nelorb_c, ii)
        end do

        deallocate (psi_sav)

        deallocate (molecorb)
        allocate (molecorb(nelorbh, ntotpsi))
        molecorb = 0.d0
        if (ipc .eq. 1) then
            call dgemm_my('N', 'N', nelorbh, ntotpsi, nelorb_diag, 1.d0, mu_c, nelorbh&
                    &, overs, nelorb_diag, 0.d0, molecorb, nelorbh, nprocu, rankopt, commopt_mpi)
        else
            call zgemm_my('N', 'N', nelorbh, ntotpsi, nelorb_diag, zone, mu_c, nelorbh&
                    &, overs, nelorb_diag, zzero, molecorb, nelorbh, nprocu, rankopt, commopt_mpi)
        end if

        indpar = iesup_c - 2*molecular*nelorbh

        do i = 1, ntotpsi
            do j = 1, nelorbh
                dup_c(ipc*(indpar + j - 1) + 1) = j
                if (ipc .eq. 2) dup_c(2*(indpar + j)) = 0.d0
            end do
            do j = nelorbh + 1, 2*nelorbh
                dup_c(ipc*(indpar + j - 1) + 1:ipc*(indpar + j)) = &
                        &molecorb(ipc*(j - nelorbh - 1) + 1:ipc*(j - nelorbh), i)
            end do
            indpar = indpar + 2*nelorbh
        end do

        deallocate (molecorb, overs, buffer, psip_loc)

    end subroutine ortho_fast

    subroutine evalovers
        use allio, only: norm_metric
        implicit none
        real*8 ddot, dnrm2, overmax, cost, r0, psiln, rc(3), jastrow_ei
        real ran
        integer iseed, nbufp

#ifdef PARALLEL
        include 'mpif.h'
#endif

        mesh = nx
        mesh = mesh*ny
        mesh = mesh*nz

        volmesh = ax*ay*az*unit_volume
        volmeshc = dcmplx(volmesh)

        !     the reference is the average ion position
        call shift_originref

        !       first part calculation overlap determinant

        if (mesh .gt. nbufd) then
            allocate (buffer(ipc*nelorbh, ipc*nbufd))
            nbufu = nbufd
        else
            allocate (buffer(ipc*nelorbh, ipc*mesh))
            nbufu = mesh
        end if
        nbufp = nbufu + 1
        allocate (overs(ipc*nelorbh, ipc*nelorbh))
        buffer = 0.d0
        overs = 0.d0

        iflagnorm = 3
        indr = 0
        ind = 0
        nbuf = 0
        do k = 1, nz
            !      x(3)=(-(nz+1)/2.d0+k)*az+rion_ref(3)
            do j = 1, ny
                !         x(2)=(-(ny+1)/2.d0+j)*ay+rion_ref(2)
                do i = 1, nx
                    indr = indr + 1

                    if (indr - (indr/nprocopt)*nprocopt .eq. rankopt) then
                        ind = ind + 1
                        !               x(1)=(-(nx+1)/2.d0+i)*ax+rion_ref(1)

                        x(:) = rion_ref(:) + (-(nz + 1)/2.d0 + k)*az*at(:, 3) + &
                                &(-(ny + 1)/2.d0 + j)*ay*at(:, 2) + (-(nx + 1)/2.d0 + i)*ax*at(:, 1)
                        call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                &, dupr, zetar, rion, psip, buffer(1, ind), nelorbh, nion, kion             &
                                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                &, indpar_tab, indorb_tab, indshell_tab, .true.)
                        if (ipc .eq. 2) then
                            call upnewwf(0, 0, 0, 1, nshellr, ioptorb, ioccup, x, 1, r, rmu        &
                                    &, dupr, zetar, rion, psip, buffer(1, ind + nbufu), nelorbh, nion, kion             &
                                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                                    &, indpar_tab, indorb_tab, indshell_tab, .false.)
                        end if

                        if (add_onebody2det) then
                            psiln = -scale_one_body
                            do jj = 1, nion
                                if (iespbc) then
                                    rc(:) = x(:) - rion(:, jj)
                                    call CartesianToCrystal(rc, 1)
                                    do kk = 1, 3
                                        rc(kk) = costz(jj)*map(rc(kk), cellscale(kk))
                                    end do
                                    r0 = norm_metric(rc, metric)
                                else
                                    rc(:) = (x(:) - rion(:, jj))*costz(jj)
                                    r0 = dsqrt(sum(rc(:)**2))
                                end if
                                psiln = psiln - jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj))*costz3(jj)
                            end do
                            buffer(1:ipc*nelorbh, ind) = buffer(1:ipc*nelorbh, ind)*dexp(psiln)
                            if (ipc .eq. 2) buffer(1:2*nelorbh, nbufu + ind) = buffer(1:2*nelorbh, nbufu + ind)*dexp(psiln)
                        end if

                        if (mod(ind, nbufu) .eq. 0) then
                            if (ipc .eq. 1) then
                                call dgemm('N', 'T', nelorbh, nelorbh, nbufu, volmesh, buffer    &
                                        &, nelorbh, buffer, nelorbh, 1.d0, overs, nelorbh)
                            else
                                call zgemm('N', 'C', nelorbh, nelorbh, nbufu, volmeshc, buffer    &
                                        &, nelorbh, buffer, nelorbh, zone, overs, nelorbh)
                                call zgemm('N', 'C', nelorbh, nelorbh, nbufu, volmeshc, buffer(1, nbufp)&
                                        &, nelorbh, buffer(1, nbufp), nelorbh, zone, overs(1, nelorbh + 1), nelorbh)
                            end if
                            nbuf = nbuf + 1
                            ind = 0
                        end if

                    end if
                end do
            end do
        end do

        if (mod(ind, nbufu) .ne. 0) then
            nleft = ind
            if (ipc .eq. 1) then
                call dgemm('N', 'T', nelorbh, nelorbh, nleft, volmesh, buffer    &
                        &, nelorbh, buffer, nelorbh, 1.d0, overs, nelorbh)
            else
                call zgemm('N', 'C', nelorbh, nelorbh, nleft, volmeshc, buffer    &
                        &, nelorbh, buffer, nelorbh, zone, overs, nelorbh)
                call zgemm('N', 'C', nelorbh, nelorbh, nleft, volmeshc, buffer(1, nbufp)    &
                        &, nelorbh, buffer(1, nbufp), nelorbh, zone, overs(1, nelorbh + 1), nelorbh)
            end if
            nbuf = nbuf + 1
        end if

#ifdef PARALLEL
        call reduce_base_real(size(overs), overs, commopt_mpi, -1)
#endif

        if (ipc .eq. 2) call conjmat(nelorbh, 2*nelorbh, overs, nelorbh)

        deallocate (buffer)

    end subroutine evalovers

    subroutine upmuctranspip

        real*8, dimension(:, :), allocatable :: mu_sav
        integer ind, jj, i1, iesupw

        ! transpip and multranspip are not changed in case of complex wfs
        ! mu_c is doubled instead

        ! reset mu_c in ANY case after reading MOs
        allocate (mu_sav(ipc*ipf*nelorbh, nelorbc_in))
        if (nelorbc_in .ne. nelorb_c .or. only_molecular) then
            ! redefine multranspip transpip
            if (.not. only_molecular) then
                deallocate (ipsip)
                allocate (ipsip(iesupc_in))
                ipsip(1:iesupc_in) = multranspip(1:iesupc_in)
            end if
            if (allocated(multranspip)) deallocate (multranspip)
            allocate (multranspip(iesup_c))
            if (.not. only_molecular) then
                multranspip(1:min(iesupc_in, iesup_c)) = ipsip(1:min(iesupc_in, iesup_c))
                if (iesup_c .gt. iesupc_in) multranspip(iesupc_in + 1:iesup_c) = 0
            else
                multranspip(1:iesup_c) = 0
            end if

            deallocate (ipsip)
            allocate (ipsip(iesup_atom*maxshell))

            ipsip = 0

            if (allocated(transpip)) then
                do i = 1, iesup_atom
                    do j = 1, multranspip(i)
                        ipsip(maxshell*(i - 1) + j) = transpip(j)%col(i)
                    end do
                end do
                deallocate (transpip)
            end if

            allocate (transpip(maxshell))
            do i1 = 1, maxshell
                if (i1 .eq. 1) then
                    allocate (transpip(i1)%col(iesup_c))
                    transpip(i1)%col = 0
                elseif (.not. only_molecular) then
                    allocate (transpip(i1)%col(iesup_atom))
                    transpip(i1)%col = 0
                end if
            end do
            if (.not. only_molecular) then
            do i = 1, iesup_atom
                do j = 1, multranspip(i)
                    transpip(j)%col(i) = ipsip(maxshell*(i - 1) + j)
                end do
            end do
            end if
            if (allocated(projm)) then
                deallocate (projm)
                allocate (projm(ipf*ipc*nelorbh*nelcol_c))
                projm = 0.d0
            end if
        end if

        if (contraction_in .gt. 0 .and. .not. only_molecular) then

            if (.not. yes_complex) then
                mu_sav(1:ipf*nelorbh, 1:nelorbc_in) = mu_c(1:ipf*nelorbh, 1:nelorbc_in)
            else
                mu_sav(1:2*ipf*nelorbh, 1:nelorbc_in) = mu_c(1:2*ipf*nelorbh, 1:nelorbc_in)
            end if
            deallocate (mu_c)
            ! it works only with a contracted basis
            if (.not. yes_complex) then
                allocate (mu_c(ipf*nelorbh, nelorb_c))
                mu_c(1:ipf*nelorbh, 1:min(nelorbc_in, nelorb_c)) = &
                    mu_sav(1:ipf*nelorbh, 1:min(nelorbc_in, nelorb_c))
            else
                allocate (mu_c(2*ipf*nelorbh, nelorb_c))
                mu_c(1:2*ipf*nelorbh, 1:min(nelorbc_in, nelorb_c)) = &
                    mu_sav(1:2*ipf*nelorbh, 1:min(nelorbc_in, nelorb_c))
            end if
        else

            if (allocated(mu_c)) deallocate (mu_c)
            if (.not. yes_complex) then
                allocate (mu_c(ipf*nelorbh, nelorb_c))
            else
                allocate (mu_c(2*ipf*nelorbh, nelorb_c))
            end if
            mu_c = 0.d0

        end if

        deallocate (mu_sav)

        ! reupdate transpip and multitranspip
        indorb = 0
        indpar = 0
        ind = 0
        do i = 1, nshell_c
            if (ioptorb_c(i) .eq. 1000000) then
                ind = ind + 1
                if (ioccup_c(ind) .eq. 1) then
                    indorb = indorb + 1
                    ! definition mu_c
                    do kk = 1, nparam_c(i)/2
                        ix = kk
                        indn = indpar + kk + nparam_c(i)/2
                        if (ipf .eq. 1 .or. only_molecular) then
                            multranspip(indn) = 1
                            transpip(1)%col(indn) = ipf*nelorbh*(indorb - 1) + ix
                        else
                            multranspip(indn) = 1
                            transpip(1)%col(indn) = nelorbh*ipf*(indorb + nelorb_at/2 - 1) + ix
                        end if
                    end do
                    indpar = indpar + nparam_c(i)
                end if
            else
                do j = 1, mult_c(i)
                    ind = ind + 1
                    if (ioccup_c(ind) .eq. 1 .and. .not. only_molecular) indorb = indorb + 1
                end do
                indpar = indpar + nparam_c(i)
            end if
        end do

        !   ! reupdate mu_c
        !   if(allocated(muc_np)) mu_c=muc_np

        do jj = 1, iesup_c
            do kk = 1, multranspip(jj)
                iy = (transpip(kk)%col(jj) - 1)/(ipf*nelorbh) + 1
                ix = transpip(kk)%col(jj) - (iy - 1)*nelorbh*ipf
                if (.not. yes_complex) then
                    mu_c(ix, iy) = dup_c(jj)
                else
                    mu_c(2*ix - 1, iy) = dup_c(2*jj - 1)
                    mu_c(2*ix, iy) = dup_c(2*jj)
                end if
            end do
        end do
        !   if(allocated(muc_np)) muc_np=mu_c

    end subroutine upmuctranspip

    ! ---------------------------------------------------------------
    ! this routine fills the matrix dup_c with the molecular orbitals
    ! coefficients contained in the matrix molecorb(:,:) for both
    ! real and complex wfs.
    !----------------------------------------------------------------

    subroutine update_dup_c

        implicit none
        integer i, j, k, ind
        integer maxref
        real*8, allocatable :: molecorb_unpaired(:, :), eigmolu(:)

        if (epsdgm .lt. 0.d0) then
            if (symmagp .and. ipc .eq. 1 .or. ipf .eq. 2) then
                maxref = max(nelorb_diag, nmoltot)
            else
                maxref = max(2*nelorb_diag, nmoltot)
            end if
        else
            if ((symmagp .and. ipc .eq. 1) .or. ipf .eq. 2) then
                if (ipf .eq. 2) then
                    !  if nelorb_diagu is odd the last eigenvector does not matter
                    maxref = nelorb_diagu
                else
                    maxref = nelorb_diag
                end if
            else
                maxref = 2*nelorb_diag
            end if
        end if

        ! iesup_c = total number of determinant parameters
        indpar = iesup_c - 2*molecular*nelorbh*ipf ! first indpar entries of dup_c are related

        if (npar_eagp .gt. 0) then
            !   Compute here unpaired molecular orbitals and eagp_pfaff according to nmolmax    chosen
            allocate (molecorb_unpaired(ipc*nelorbh*ipf, ndiff), eigmolu(nelorb_diagu/2))

            do i = nelorb_diagu/2 - nmolmin + 1, nelorb_diagu/2
                eigmolu(i) = dble(i - nelorb_diag + nmolmin)
            end do
            do i = 1, nelorb_diagu/2 - nmolmin
                eigmolu(i) = eigmol(i)/eigmol(nelorb_diagu/2 - nmolmin + 1)
                if (power .ne. 1.d0) eigmolu(i) = abs(eigmolu(i))**power*sign(1.d0, eigmolu(i))
            end do

            molecorb_unpaired = 0.d0
            eagp_pfaff = 0.d0
            if (ipc .eq. 1) then
                do j = 1, ndiff
                    do i = 1, ndiff
                        do k = 1, nelorb_diagu/2
                            eagp_pfaff(i, j) = eagp_pfaff(i, j)&
                                    & - molecorb(ipf*nelorbh + i, 2*k - 1)*molecorb(ipf*nelorbh + j, 2*k)*eigmolu(k)&
                                    & + molecorb(ipf*nelorbh + i, 2*k)*molecorb(ipf*nelorbh + j, 2*k - 1)*eigmolu(k)
                        end do
                    end do
                    do k = 1, nelorb_diagu/2
                        molecorb_unpaired(:, j) = molecorb_unpaired(:, j)&
                                & - molecorb(1:ipf*nelorbh, 2*k - 1)*molecorb(ipf*nelorbh + j, 2*k)*eigmolu(k)&
                                & + molecorb(1:ipf*nelorbh, 2*k)*molecorb(ipf*nelorbh + j, 2*k - 1)*eigmolu(k)
                    end do
                end do
            else

                call update_eagp(ndiff, nelorb_diagu, nelorbh, nelorbpf, eigmolu, molecorb&
                        &, molecorb_unpaired, eagp_pfaff)

            end if
        end if

        ! to atomic contracted orbitals
        if (.not. yes_complex) then
            !------------------------------------------------
            !                 real w.f.
            !------------------------------------------------
            dup_c(1:indpar) = psip(1:indpar)

            do i = 1, nmoltot - ndiff
                ! filling 1st set: pointers to atomic orbitals
                do j = 1, nelorbh*ipf
                    dup_c(indpar + j) = j
                end do
                ! filling 2nd set: coefficients
                do j = ipf*nelorbh + 1, 2*ipf*nelorbh
                    if (symmagp .and. ipf .eq. 1) then
                        if (i .le. maxref) dup_c(indpar + j) = molecorb(j - ipf*nelorbh, maxref - i + 1)
                    else
                        ! The even remain even orbitals and the odd remain odd orbitals.
                        if (mod(maxref, 2) .eq. 0) then
                            if (mod(i, 2) .eq. 1) then
                                if (i .lt. maxref) dup_c(indpar + j) = molecorb(j - ipf*nelorbh, maxref - i)
                            else
                                if (i .le. maxref) dup_c(indpar + j) = molecorb(j - ipf*nelorbh, maxref - i + 2)
                            end if
                        else
                            if (i .le. maxref) dup_c(indpar + j) = molecorb(j - ipf*nelorbh, maxref - i + 1)
                        end if

                    end if
                end do
                indpar = indpar + 2*ipf*nelorbh
            end do ! enddo over nmoltot
            ! unpaired orbitals
            do i = 1, ndiff
                do j = 1, ipf*nelorbh
                    dup_c(indpar + j) = j
                end do
                if (symmagp .or. ipf .eq. 2) then
                    if (npar_eagp .gt. 0) then
                        do j = ipf*nelorbh + 1, 2*ipf*nelorbh
                            dup_c(indpar + j) = molecorb_unpaired(j - ipf*nelorbh, i)
                        end do
                    else
                        do j = ipf*nelorbh + 1, 2*ipf*nelorbh
                            if (maxref - nmoltot + i .gt. 0) dup_c(indpar + j) = molecorb(j - ipf*nelorbh, maxref - nmoltot + i)
                        end do
                    end if
                else
                    do j = nelorbh + 1, 2*nelorbh
                        if (maxref - nmoltot + i .gt. 0) dup_c(indpar + j) = molecorb(j - nelorbh, maxref - nmoltot + i)
                    end do
                end if
                indpar = indpar + 2*ipf*nelorbh
            end do
            !------------------------------------------------
            !                 complex w.f.
            !------------------------------------------------
        else

            dup_c(1:2*indpar) = psip(1:2*indpar)

            do i = 1, nmoltot - ndiff

                do j = 1, 2*ipf*nelorbh, 2
                    dup_c(2*indpar + j) = j/2 + 1
                    dup_c(2*indpar + j + 1) = 0.d0 ! put imaginary part to zero
                end do

                do j = 2*ipf*nelorbh + 1, 4*ipf*nelorbh, 2

                    if (mod(maxref, 2) .eq. 0) then
                        if (mod(i, 2) .eq. 1) then
                            if (i .lt. maxref) then
                                dup_c(2*indpar + j) = molecorb(j - 2*ipf*nelorbh, maxref - i)
                                dup_c(2*indpar + j + 1) = molecorb(j + 1 - 2*ipf*nelorbh, maxref - i)
                            end if
                        else
                            if (i .le. maxref) then
                                dup_c(2*indpar + j) = molecorb(j - 2*ipf*nelorbh, maxref - i + 2)
                                dup_c(2*indpar + j + 1) = molecorb(j - 2*ipf*nelorbh + 1, maxref - i + 2)
                            end if
                        end if
                    else
                        if (i .le. maxref) then
                            dup_c(2*indpar + j) = molecorb(j - 2*ipf*nelorbh, maxref - i + 1)
                            dup_c(2*indpar + j + 1) = molecorb(j + 1 - 2*ipf*nelorbh, maxref - i + 1)
                        end if
                    end if
                end do
                indpar = indpar + 2*ipf*nelorbh
            end do

            do i = 1, ndiff
                do j = 1, 2*ipf*nelorbh, 2
                    dup_c(2*indpar + j) = (j + 1)/2
                    dup_c(2*indpar + j + 1) = 0.d0
                end do
                if (npar_eagp .gt. 0) then
                    do j = 2*ipf*nelorbh + 1, 4*ipf*nelorbh, 2
                        dup_c(2*indpar + j) = molecorb_unpaired(j - 2*ipf*nelorbh, i)
                        dup_c(2*indpar + j + 1) = molecorb_unpaired(j + 1 - 2*ipf*nelorbh, i)
                    end do
                else
                    do j = 2*ipf*nelorbh + 1, 4*ipf*nelorbh, 2
                        if (maxref - nmoltot + i .gt. 0) then
                            dup_c(2*indpar + j) = molecorb(j - 2*ipf*nelorbh, maxref - nmoltot + i)
                            dup_c(2*indpar + j + 1) = molecorb(j + 1 - 2*ipf*nelorbh, maxref - nmoltot + i)
                        end if
                    end do
                end if
                indpar = indpar + 2*ipf*nelorbh
            end do

        end if

        if (allocated(molecorb_unpaired)) deallocate (molecorb_unpaired, eigmolu)
        ! COMPLEX DEB
        !write(6,*) ' Final dup_c ',size(dup_c)
        !k=iesup_c-2*molecular*nelorbh
        !do i=1,nmoltot-ndiff
        !   if(yes_complex) then
        !      do j=1,4*nelorbh,2
        !         write(6,*) i,j/2+1,dup_c(2*k+j),dup_c(2*k+j+1)
        !      enddo
        !   else
        !      do j=1,2*nelorbh
        !         write(6,*) i,j,dup_c(k+j)
        !      enddo
        !   endif
        !   k=k+2*nelorbh
        !enddo
        !stop

    end subroutine update_dup_c

    subroutine projectder(nelorbh, nelorb, nmolmatdo, nmol, projmat&
            &, nmolmat, dermat, psip, yesmin, symmagp)
        use constants, only: ipc, ipf
        implicit none
        integer nelorb, nelorbh, nmolmatdo, i, j, nmol, nmolmat, yesmin
        real*8 dermat(ipc*nelorb, nelorb), psip(*), projmat(ipc*nelorbh, nmolmat, *)&
                &, cost
        logical symmagp
        real*8, dimension(:, :), allocatable :: dersav

        !   Test, avoiding the projection means to symmetrize by multiplying the off by 2.
        !          ! symmetrization
        !           do i=1,nelorbh
        !              do j=i+1,nelorbh
        !                 cost=dermat(i,j)+dermat(j,i)
        !                 dermat(i,j)=cost
        !                 dermat(j,i)=cost
        !              enddo
        !           enddo
        !        return

        if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then

            ! Lambda is assumed not symmetric in the input (given by AAD) and symmetrized at the
            ! end to avoid stupid mistakes.
            ! Be input Lambda=dermat
            ! The projector is written as P_1 [projmat(:,:1)] x P_2 [projmat(:,:,2)]^dag
            ! The projected derivative without any constraint on the AGP can be written as
            !  D= P_2 P_1^dag Lambda + Lambda P_1 P_2^dag - P_2 P_1^dag Lambda P_1 P_2^dag
            !  Then the constrained derivative with symmetric AGP is given by:
            !  D--> D+D^dag - delta_ij D_ii   as the diagonal elements are not multiplied by two.
            !  Now D+D^dag can be written as 2 (C+C^dag):
            !  C= A -1/2 B P_2^dag.
            !  where A=1/2 P_2 P_1^dag (Lambda + Lambda^dag)
            !  B= A P_1
            ! The above is equivalent to:
            !  Lambda = 2 ( C+ C^dag)-2 delta_ij C_ii

            !     The algorithm.
            do i = 1, nelorbh
                do j = i + 1, nelorbh
                    !         cost=(dermat(i,j)+dermat(j,i))*0.5d0
                    dermat(i, j) = dermat(i, j)/2.d0
                    dermat(j, i) = dermat(j, i)/2.d0
                end do
            end do
#ifdef _OFFLOAD
!$omp target data map(dermat) map(to:projmat(:,: ,1:2))
#endif
            !
            !   A=  1/2 P_2 P_1^dag (Lambda + Lambda^dag)--> dermat
            call dgemm_('T', 'N', nmolmat, nelorbh, nelorbh, 1.d0, projmat, nelorbh&
                    &, dermat, nelorb, 0.d0, psip, nmolmat)
            call dgemm_('N', 'N', nelorbh, nelorbh, nmolmat, 1.d0, projmat(1, 1, 2)&
                    &, nelorbh, psip, nmolmat, 0.d0, dermat, nelorb)
            !   B= A P_1--> psip
            call dgemm_('N', 'N', nelorbh, nmolmat, nelorbh, 1.d0, dermat, nelorb   &
                    &, projmat, nelorbh, 0.d0, psip, nelorbh)
            !   C--> A -1/2 B P_2^dag --> dermat
            call dgemm_('N', 'T', nelorbh, nelorbh, nmolmat, -0.5d0, psip, nelorbh   &
                    &, projmat(1, 1, 2), nelorbh, 1.d0, dermat, nelorb)
#ifdef _OFFLOAD
!$omp end target data
#endif
            do i = 1, nelorbh
                ! symmetrization
                !       Lambda = 2*(C+C^dag)-2*delta_ij C_ii
                dermat(i, i) = 2.d0*dermat(i, i) ! NB we do not multiply by 4 but by 2
                do j = i + 1, nelorbh
                    cost = 2.d0*(dermat(i, j) + dermat(j, i))
                    dermat(i, j) = cost
                    dermat(j, i) = cost
                end do
            end do
        else

            ! Here see ../doc/parbcs.tex last part for the algorithm.
            ! First we emply a projected  derivative, without assuming any symmetry for the AGP.
            ! Then we symmetrize according to logical variables: symmagp, yes_hermite.
            allocate (dersav(ipc*nelorbh, nelorbh))
            dersav(1:ipc*nelorbh, 1:nelorbh) = dermat(1:ipc*nelorbh, 1:nelorbh)
#ifdef _OFFLOAD
!$omp target data map(dermat) map(to:projmat(:,: ,1:4))
#endif

            if (ipc .eq. 1) then
                !      compute P^L lambda -lambda
                call dgemm_('T', 'N', nmolmat, nelorbh, nelorbh, 1.d0, projmat(1, 1, 2), nelorbh&
                        &, dermat, nelorb, 0.d0, psip, nmolmat)
                call dgemm_('N', 'N', nelorbh, nelorbh, nmolmat, 1.d0, projmat(1, 1, 4)&
                        &, nelorbh, psip, nmolmat, -1.d0, dermat, nelorb)

                !        ( P^L lambda -lambda) -( P^L lambda -lambda) P^R
                call dgemm_('N', 'N', nelorbh, nmolmatdo, nelorbh, 1.d0, dermat, nelorb   &
                        &, projmat, nelorbh, 0.d0, psip, nelorbh)
                call dgemm_('N', 'T', nelorbh, nelorbh, nmolmatdo, -1.d0, psip, nelorbh   &
                        &, projmat(1, 1, 3), nelorbh, 1.d0, dermat, nelorb)

            else

                !      compute P^L^+ lambda -lambda
                call zgemm_('C', 'N', nmolmat, nelorbh, nelorbh, zone, projmat(1, 1, 2), nelorbh&
                        &, dermat, nelorb, zzero, psip, nmolmat)

                call zgemm_('N', 'N', nelorbh, nelorbh, nmolmat, zone, projmat(1, 1, 4)&
                        &, nelorbh, psip, nmolmat, zmone, dermat, nelorb)
                !        ( P^L lambda -lambda) -( P^L lambda -lambda) P^R
                call zgemm_('N', 'N', nelorbh, nmolmatdo, nelorbh, zone, dermat, nelorb   &
                        &, projmat, nelorbh, zzero, psip, nelorbh)
                call zgemm_('N', 'C', nelorbh, nelorbh, nmolmatdo, zmone, psip, nelorbh   &
                        &, projmat(1, 1, 3), nelorbh, zone, dermat, nelorb)

            end if
#ifdef _OFFLOAD
!$omp end target data
#endif
            !        adding back the original lambda
            do i = 1, ipc*nelorbh
                do j = 1, nelorbh
                    dermat(i, j) = dermat(i, j) + dersav(i, j)
                end do
            end do
            deallocate (dersav)
        end if
    end subroutine projectder

    subroutine projectmat(nelorbh, nelorb, nmolmatdo, nmol, projmat, nmolmat&
            &, dermat, psip, symmagp, yesmin)
        use constants, only: ipc, ipf, zone, zmone, zzero
        implicit none
        integer nelorb, nelorbh, nmolmatdo, nmax, i, j, nelorbp, yesmin, nmolmat&
                &, nmol
        real*8 dermat(ipc*nelorbh, *), psip(*), projmat(ipc*nelorbh, nmolmat, *)&
                &, cost
        logical symmagp
        real*8, dimension(:, :), allocatable :: dermat_sav
        !     max memory psip nmolmax/nmol*nelorbh for yesmin=1/2
        !     yesmin =1  Projection on rankopt nmolmax SDV:
        !M--> (I-P^nmolmnin_L) M (I -P^nmolmin_R)- ( I -P^nmolmax_L) M (I -P^nmolmax_R)
        !         symmagp=.true.
        !     P_L^n = sum i<=n  phi1_L^i  phi2_L^i dag
        !     P_R^n = sum i<=n  phi2 R^i  phi1_R^i dag
        !     phi2= S phi1  where S is the overlap matrix and phi1 are S-orthogonal
        !         symmagp=.false.
        !         P_L=   phi2 phi4^+
        !         P_R=   phi3 phi1^+
        !         phi4 and phi3 contain the left and the right overlap matrices, respectively.
        !     vectors
        !#ifdef __CASO
        !    nprocu=nprocopt
        !#else
        !    nprocu=1
        !#endif

        if (yesmin .eq. 1) then

            allocate (dermat_sav(ipc*nelorbh, nelorbh))
            dermat_sav(1:ipc*nelorbh, 1:nelorbh) = dermat(1:ipc*nelorbh, 1:nelorbh)

            if (symmagp .and. ipc .eq. 1 .and. ipf .eq. 1) then

                call dgemm_my('T', 'N', nmolmat, nelorbh, nelorbh, 1.d0, projmat(1, 1, 2), nelorbh&
                        &, dermat, nelorbh, 0.d0, psip, nmolmat, nprocu, rankopt, commopt_mpi)
                call dgemm_my('N', 'N', nelorbh, nelorbh, nmolmat, -1.d0, projmat, nelorbh&
                        &, psip, nmolmat, 1.d0, dermat, nelorbh, nprocu, rankopt, commopt_mpi)

                call dgemm_my('N', 'N', nelorbh, nmolmat, nelorbh, 1.d0, dermat, nelorbh&
                        &, projmat(1, 1, 2), nelorbh, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
                call dgemm_my('N', 'T', nelorbh, nelorbh, nmolmat, -1.d0, psip, nelorbh   &
                        &, projmat, nelorbh, 1.d0, dermat, nelorbh, nprocu, rankopt, commopt_mpi)

            else

                if (ipc .eq. 1) then

                    call dgemm_my('T', 'N', nmolmat, nelorbh, nelorbh, 1.d0, projmat(1, 1, 4), nelorbh&
                            &, dermat, nelorbh, 0.d0, psip, nmolmat, nprocu, rankopt, commopt_mpi)
                    call dgemm_my('N', 'N', nelorbh, nelorbh, nmolmat, -1.d0, projmat(1, 1, 2), nelorbh&
                            &, psip, nmolmat, 1.d0, dermat, nelorbh, nprocu, rankopt, commopt_mpi)

                    call dgemm_my('N', 'N', nelorbh, nmolmatdo, nelorbh, 1.d0, dermat, nelorbh   &
                            &, projmat(1, 1, 3), nelorbh, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
                    call dgemm_my('N', 'T', nelorbh, nelorbh, nmolmatdo, -1.d0, psip, nelorbh   &
                            &, projmat, nelorbh, 1.d0, dermat, nelorbh, nprocu, rankopt, commopt_mpi)

                else

                    !  Here P_L = psi psi^dag
                    !        (I-P_L) D

                    call zgemm_my('C', 'N', nmolmat, nelorbh, nelorbh, zone, projmat(1, 1, 4), nelorbh&
                            &, dermat, nelorbh, zzero, psip, nmolmat, nprocu, rankopt, commopt_mpi)
                    call zgemm_my('N', 'N', nelorbh, nelorbh, nmolmat, zmone, projmat(1, 1, 2), nelorbh&
                            &, psip, nmolmat, zone, dermat, nelorbh, nprocu, rankopt, commopt_mpi)
                    !     (I-P_L) D (I-P_R)   here P_R =\bar psi \bar psi^dag

                    call zgemm_my('N', 'N', nelorbh, nmolmatdo, nelorbh, zone, dermat, nelorbh   &
                            &, projmat(1, 1, 3), nelorbh, zzero, psip, nelorbh, nprocu, rankopt, commopt_mpi)
                    call zgemm_my('N', 'C', nelorbh, nelorbh, nmolmatdo, zmone, psip, nelorbh   &
                            &, projmat, nelorbh, zone, dermat, nelorbh, nprocu, rankopt, commopt_mpi)
                end if

            end if

        end if

        if (yesmin .eq. 1) then
            dermat(1:ipc*nelorbh, 1:nelorbh) = &
                    &dermat_sav(1:ipc*nelorbh, 1:nelorbh) - dermat(1:ipc*nelorbh, 1:nelorbh)
            deallocate (dermat_sav)
        end if

    end subroutine projectmat

end module convertmod

function tracem(n, a, lda)
    implicit none
    integer i, j, n, lda
    real*8 tracem, a(lda, *)
    tracem = 0.d0
    do i = 1, n
        do j = 1, n
            tracem = tracem + a(i, j)*a(j, i)
        end do
    end do
    return
end

subroutine invsym(n, a, lda, info)
    use constants, only: ipc
    implicit none
    integer n, lda, info, i, j
    real*8 a(ipc*lda, n)
    !       Now the inverse of this matrix using Cholesky
    if (ipc .eq. 1) then
        call dpotrf('L', n, a, lda, info)
        call dpotri('L', n, a, lda, info)
    else
        call zpotrf('L', n, a, lda, info)
        call zpotri('L', n, a, lda, info)
    end if
    do i = 1, n
        do j = 1, i
            a(ipc*(j - 1) + 1:ipc*j, i) = a(ipc*(i - 1) + 1:ipc*i, j)
        end do
    end do
    return
end subroutine invsym

subroutine update_eagp(ndiff, nelorb_diagu, nelorbh, nelorbpf, eigmol, molecorb&
        &, molecorb_unpaired, eagp_pfaff)
    implicit none
    integer ndiff, nelorbh, nelorb_diagu, nelorbpf, i, j, k
    real*8 eigmol(*)
    complex*16 molecorb(nelorbpf, *), molecorb_unpaired(2*nelorbh, ndiff)
    complex*16 eagp_pfaff(ndiff, ndiff)
    molecorb_unpaired = (0.d0, 0.d0)
    eagp_pfaff = (0.d0, 0.d0)
    do j = 1, ndiff
        do i = 1, ndiff
            do k = 1, nelorb_diagu/2
                eagp_pfaff(i, j) = eagp_pfaff(i, j)&
                        & - molecorb(2*nelorbh + i, 2*k - 1)*molecorb(2*nelorbh + j, 2*k)*eigmol(k)&
                        & + molecorb(2*nelorbh + i, 2*k)*molecorb(2*nelorbh + j, 2*k - 1)*eigmol(k)
            end do
        end do
        do k = 1, nelorb_diagu/2
            molecorb_unpaired(:, j) = molecorb_unpaired(:, j)&
                    & - molecorb(1:2*nelorbh, 2*k - 1)*molecorb(2*nelorbh + j, 2*k)*eigmol(k)&
                    & + molecorb(1:2*nelorbh, 2*k)*molecorb(2*nelorbh + j, 2*k - 1)*eigmol(k)
        end do
    end do
    return
end
