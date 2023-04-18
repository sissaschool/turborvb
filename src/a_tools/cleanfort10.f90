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

program join

    implicit none
    integer i, j, k, ii, kk, nelup_1, nel_1, nion_1, nelup_2, nel_2, nion_2    &
            &, nelup, nel, nion, niesd, iesdr_1, iesdr_2, iesdim, iesupr_1, iesupr_2    &
            &, npar3body, npar3body_1, npar3body_2, nshell, nshell_1, nshell_2       &
            &, nshellj, nshellj_1, nshellj_2, iesdr, iesupr, nnozero, nnozero_1       &
            &, nnozero_2, nnozeroj, nnozeroj_1, nnozeroj_2, iesupind_1, iesupind_2   &
            &, iesupind, iesmind_1, iesmind_2, iesmind, iesfreer_1, iesfreer_2       &
            &, iesfreer, iesswr_1, iesswr_2, iesswr, ieskinr_1, ieskinr_2, ieskinr    &
            &, indpar, occ_1, occ_2, occtot_1, occtotj, occtotj_2, ind, indj, iflag     &
            &, neldo_1, neldo_2, occtot_2, occj_1, occj_2, iflagdetpp, iflagjaspp     &
            &, iflagjasss, iflagdetss, iflagdetsp, iflagjassp, neldiff_1, neldiff_2  &
            &, itest, iflagdetdd, iflagjasdd, iflagdetda, iflagjasda, maxnozero      &
            &, iflagdetspa, iflagjasspa, iflagdetdaa, iflagjasdaa, iflagr, inds      &
            &, maxnozeroj, indpar_1, indpar_2, ind_1, indj_1, fill_1, fill_2, jj       &
            &, iflag_1, nrep, norb, indorb, neldiff, pipsz, iesdrr, norbj              &
            &, count_s, count_p, count_d, count_f, halfpar, iesupind_new, iii         &
            &, npar3body_new, count_sj, count_pj, count_dj, count_fj, halfparj       &
            &, iesmind_new, iesup_new, menouno, sumnapar, sumnaparj, maxic           &
            &, indtot, indtot0, indorb0, indorbi, indorbj, irankdet, info, ind0        &
            &, icost, ixt, iyt, inddet, ipc, indpar_at, ipj, nelorb_at, nozero_at&
            &, count_at, ipf, ndimh, jmax, numpaired, npar_eagp, ix, iy, iesswr_eagp&
            &, nnozero_eagp, case_map

    real(8) scale, range, rangej, distmin, dist, amax, cost, costmax, sign       &
            &, rs1, Ly1, Lz1, phase(3), phase_down(3), rtarget, u(2), v(2), costc(2), celldm1(9)
    real(8), dimension(:), allocatable :: zetar_1       &
            &, dupr_1, vj, vjur_1, detmat_1, overlap_1           &
            &, jasmat_1, scaleorb, jasmatsz_1                                     &
            &, eig, psip, atomic_number, sortvect
    real(8), dimension(:, :), allocatable :: eagp_pfaff
    integer, dimension(:), allocatable :: nparam_1      &
            &, ic_1, kion_1, mult_1, ioptorb_1        &
            &, kionj_1, multj_1, ioptorbj_1            &
            &, ioccup_1, nparamj_1, ix_1, iy_1       &
            &, ioccupj_1, icd_1, icj_1, ixj_1, iyj_1                      &
            &, icp_1, icq_1, ioptorbja_1, ioptorba                        &
            &, map, mapj, occ_shell, occj_shell, ind_old                            &
            &, ipiv, ns, halfp, cek_w, cek_wj                                       &
            &, nsj, halfpj, mult_new, mult_newj, ipipsav, indexo, opt_par
    real(8), dimension(:, :), allocatable :: rion_1              &
            &, UU
    real(8), dimension(:, :), allocatable :: lambda_mat, lambda_matj       &
            &, lambda_matszj
    real(8), dimension(:, :), allocatable :: psi_in, psi_out, umat, umatj
    integer, dimension(:, :), allocatable :: ipsip_1      &
            &, ipsipd_1, ipsipj_1, ipsipp_1                                       &
            &, ipsipq_1, address, addressorb, addressorbtot
    character(2) chars
    character(3) checkpbc
    character(5) checkpbc_c
    logical iessz, iespbc, check, iesdet, yesbump, mol_off, molj_off&
            &, mol_ies, symmagp, yes9, yes9j, yes_complex, yes_crystal, forcecomplex&
            &, forcesymm, yes_hermite, yesalldet, yesfixphase, pfaffup, yes_tilted&
            &, yes_crystalj, chosen_map
    logical, dimension(:), allocatable :: cek_occ, done, signpar, signparj, maxdupr

    !    ZZZ   lines added by Zen
    integer flagconstrains
    !    ZZZ   end lines added by Zen

    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'cleanfort10'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added

    !    ZZZ   lines added by Zen
    !          to include a flag that removes the constrains
    !          on the contracted coefficients
    !            if flagconstrains is
    !            0: constrains both on DET. and JAS. contracted coeff.
    !            1: no constrains in DET. contracted coeff.
    !            2: no constrains in JAS. contracted coeff.
    !            3: no constrains in all contracted coeff.
    flagconstrains = 0
    if (str .eq. "noconstr" .or. str .eq. "noconstrains" .or. &
            & str .eq. "NOCONSTR" .or. str .eq. "NOCONSTRAINS") then
        flagconstrains = 3
        write (6, *) ' Warning: read flag for NO-constrains on contracted orbitals! '
    end if
    if (str .eq. "noconstrjas" .or. str .eq. "noconstrainsjas" .or. &
            & str .eq. "NOCONSTRJAS" .or. str .eq. "NOCONSTRAINSJAS") then
        flagconstrains = 2
        write (6, *) ' Warning: read flag for NO-constrains on contracted orbitals of Jastrow! '
    end if
    if (str .eq. "noconstrdet" .or. str .eq. "noconstrainsdet" .or. &
            & str .eq. "NOCONSTRDET" .or. str .eq. "NOCONSTRAINSDET") then
        flagconstrains = 1
        write (6, *) ' Warning: read flag for NO-constrains on contracted orbitals of det. part! '
    end if
    !    ZZZ   end lines added by Zen

    open (unit=10, file='fort.10', status='old', form='formatted')
    open (unit=11, file='fort.10_clean', status='unknown'           &
            &, form='formatted')

    yes_crystal = .false.
    yes_complex = .false.
    yes_tilted = .false.

    rewind (10)
    read (10, *) chars, checkpbc_c
    if (checkpbc_c .eq. 'PBC_T') then
        iespbc = .true.
        yes_crystal = .true.
        yes_tilted = .true.
    elseif (checkpbc_c .eq. 'PBC_C') then
        iespbc = .true.
        yes_crystal = .true.
    else
        rewind (10)
        read (10, *) chars, checkpbc
        if (checkpbc .eq. 'PBC') then
            iespbc = .true.
        else
            iespbc = .false.
        end if
    end if
    rewind (10)
    if (iespbc) then
        read (10, *)
        if (.not. yes_crystal) then
            read (10, *) rs1, Ly1, Lz1, phase(:)
        elseif (yes_tilted) then
            read (10, *) celldm1(1:9), phase(:), phase_down(:)
        else
            read (10, *) rs1, Ly1, Lz1, phase(:), phase_down(:)
        end if
    end if
    if (yes_crystal) then
        yes_hermite = .true.
    else
        yes_hermite = .false.
    end if

    if (iespbc) then
        do i = 1, 3
            if (abs(phase(i) - nint(phase(i))) .ne. 0.5d0 .and. abs(phase(i) - nint(phase(i))) .ne. 0.d0) yes_complex = .true.
            ! treat gamma and other real boundaries, like (0.5,0.5,0.5), without complex w.f.
            ! unless complex algorithm is forced on input (nshell<0)
            if (yes_crystal) then
                if (abs(phase_down(i) - nint(phase_down(i))) .ne. 0.5d0 .and. abs(phase_down(i) - nint(phase_down(i)))&
                        &.ne. 0.d0) yes_complex = .true.
                if ((abs(phase(i) - nint(phase(i))) .ne. 0.5 .and.&
                        &(phase_down(i) - nint(phase_down(i)) .ne. nint(phase(i)) - phase(i))) .or.&
                        &(abs(phase(i) - nint(phase(i))) .eq. 0.5 .and.&
                                &abs(phase_down(i) - nint(phase_down(i))) .ne. 0.5)) yes_hermite = .false.
            end if
        end do
    end if

    rewind (11)

    read (10, *)
    read (10, *) nelup_1, nel_1, nion_1
    if (nelup_1 .le. 0) then
        pfaffup = .true.
        nelup_1 = -nelup_1
    else
        pfaffup = .false.
    end if

    if (nel_1 .lt. 0) then
        ipf = 2
        nel_1 = -nel_1
    else
        ipf = 1
    end if

    if (mod(nelup_1, nel_1) .ne. 0) then
        numpaired = 2*(nelup_1/nel_1)
        nelup_1 = mod(nelup_1, nel_1)
    else
        numpaired = 2*(nelup_1/nel_1 - 1)
        nelup_1 = nel_1
    end if
    if (nion_1 .lt. 0) then
        yesbump = .true.
        nion_1 = -nion_1
    else
        yesbump = .false.
    end if
    neldo_1 = nel_1 - nelup_1

    if (ipf .eq. 2) then
        neldiff_1 = mod(nel_1, 2) + numpaired
        npar_eagp = (neldiff_1*(neldiff_1 - 1))/2
    else
        npar_eagp = 0
        neldiff_1 = nelup_1 - neldo_1
    end if

    nelup = nelup_1
    nel = nel_1
    nion = nion_1

    read (10, *)
    read (10, *) nshell_1, nshellj_1

    forcecomplex = .false.
    if (nshell_1 .lt. 0) then
        yes_complex = .true.
        nshell_1 = -nshell_1
        forcecomplex = .true.
    end if
    yes_crystalj = .false.
    if (nshellj_1 .lt. 0) then
        yes_crystalj = .true.
        nshellj_1 = -nshellj_1
    end if

    if (yes_complex) then
        ipc = 2
    else
        ipc = 1
    end if

    nshell = nshell_1
    nshellj = nshellj_1

    read (10, *)
    if (npar_eagp .gt. 0) then
        read (10, *) iesdrr, iesupr_1, npar3body_1, iesswr_eagp, nnozero_eagp
    else
        iesswr_eagp = 0
        nnozero_eagp = 0
        read (10, *) iesdrr, iesupr_1, npar3body_1
    end if

    chosen_map = .false.

    if (abs(iesdrr) .ge. 100) then
        chosen_map = .true.
        case_map = abs(iesdrr)/100
        if (iesdrr .gt. 0) then
            iesdrr = iesdrr - 100*case_map
        else
            iesdrr = iesdrr + 100*case_map
        end if
    end if

    if (iesdrr .eq. -8 .or. iesdrr .eq. -9 .or. iesdrr .eq. -18 .or. iesdrr .eq. -19&
            &.or. iesdrr .eq. -16 .or. iesdrr .eq. -28 .or. iesdrr .eq. -29) then
        iesdr_1 = -7
        iessz = .true.
    else
        iesdr_1 = iesdrr
        iessz = .false.
    end if

    if (iesdrr .eq. -12 .or. iesdrr .eq. -22 .or. iesdrr .eq. -26 .or. iesdrr .eq. -27 &
        .or. iesdrr .eq. -30 .or. iesdrr .eq. -31) then
        ipj = 2
    else
        ipj = 1
    end if

    iesdr = iesdr_1

    iesdim = max(abs(iesdr_1), 1)

    iesupr = iesupr_1
    npar3body = npar3body_1

    yesalldet = .false.
    yesfixphase = .false.

    read (10, *)
    read (10, *) nnozero_1, nnozeroj_1

    if (nnozero_1 .lt. 0) then
        forcesymm = .true.
        nnozero_1 = -nnozero_1
        yes_hermite = .false.
    else
        forcesymm = .false.
    end if

    nnozero = nnozero_1
    nnozeroj = nnozeroj_1

    read (10, *)
    read (10, *) iesupind_1, iesmind_1

    iesupind = iesupind_1
    iesmind = iesmind_1

    read (10, *)
    read (10, *) iesfreer_1, iesswr_1, ieskinr_1, iflag_1

    ieskinr = ieskinr_1

    iesfreer = iesfreer_1
    iesswr = iesswr_1

    allocate (rion_1(3, nion_1))
    allocate (ipsip_1(ieskinr_1, 6*nion_1))

    allocate (zetar_1(nion_1), dupr_1(ipc*iesupr_1)     &
            &, vjur_1(npar3body_1), ic_1(ieskinr_1), vj(iesdim)                   &
            &, kion_1(nshell_1), mult_1(nshell_1), ioptorb_1(nshell_1)            &
            &, nparam_1(nshell_1), occ_shell(nshell_1), occj_shell(nshellj_1)     &
            &, kionj_1(nshellj_1), multj_1(nshellj_1), ioptorbj_1(nshellj_1)      &
            &, nparamj_1(nshellj_1)                                             &
            &, ix_1(nnozero_1 + nnozero_eagp), iy_1(nnozero_1 + nnozero_eagp)             &
            &, detmat_1(ipc*nnozero_1), atomic_number(nion_1))
    allocate (signpar(iesupr_1), signparj(npar3body_1))
    if (npar_eagp .gt. 0) then
        allocate (eagp_pfaff(ipc*neldiff_1, neldiff_1))
    end if

    read (10, *)
    do j = 1, nion_1
        read (10, *) zetar_1(j), atomic_number(j), (rion_1(i, j), i=1, 3)
    end do

    read (10, *)
    do i = 1, ieskinr_1
        read (10, *) ic_1(i), (ipsip_1(i, j), j=1, 2*abs(ic_1(i)))
    end do

    read (10, *)

    if (iesdr_1 .ne. 0) then
        if (abs(iesdr_1) .le. 3) then
            symmagp = .true.
            read (10, *) (vj(i), i=1, abs(iesdr_1))
            niesd = abs(iesdr_1)
        else
            read (10, *) niesd, (vj(i), i=1, abs(niesd))
            if (niesd .lt. 0) then
                symmagp = .false.
                niesd = -niesd
            else
                symmagp = .true.
            end if
        end if
    else
        read (10, *) niesd
        if (niesd .lt. 0) then
            symmagp = .false.
            niesd = 0
        else
            symmagp = .true.
        end if
    end if

    if (.not. symmagp .and. yes_hermite) yes_hermite = .false.

    indpar_at = 0
    dupr_1 = 0.d0
    read (10, *)
    indpar = 0
    nshell_2 = 0
    yes9 = .false.
    allocate (maxdupr(iesupr_1))
    maxdupr = .false.

    do i = 1, nshell_1
        read (10, *) mult_1(i), nparam_1(i), ioptorb_1(i)
        if (ioptorb_1(i) .eq. 900000) yes9 = .true.
        nshell_2 = nshell_2 + 1
        if (yes_complex .and. nparam_1(i) .gt. 1) then
            read (10, *) kion_1(i), (dupr_1(2*(indpar + j) - 1), j=1, nparam_1(i)/2), & ! reading exponents
                    &(dupr_1(2*indpar + nparam_1(i) + j), j=1, nparam_1(i)) ! reading complex coefficients
            if (ipf .eq. 2 .and. ioptorb_1(i) .ge. 900000) then
                costmax = -1.d0
                do j = 1, nparam_1(i)/2, 2
                    cost = dupr_1(2*indpar + nparam_1(i) + j)**2 + dupr_1(2*indpar + nparam_1(i) + j + 1)**2
                    if (cost .gt. costmax) then
                        costmax = cost
                        jmax = j
                    end if
                end do
                maxdupr(indpar + nparam_1(i)/2 + (jmax + 1)/2) = .true.
                costmax = -1.d0
                do j = nparam_1(i)/2 + 1, nparam_1(i), 2
                    cost = dupr_1(2*indpar + nparam_1(i) + j)**2 + dupr_1(2*indpar + nparam_1(i) + j + 1)**2
                    if (cost .gt. costmax) then
                        costmax = cost
                        jmax = j
                    end if
                end do
                maxdupr(indpar + nparam_1(i)/2 + (jmax + 1)/2) = .true.
            else
                costmax = -1.d0
                do j = 1, nparam_1(i), 2
                    cost = dupr_1(2*indpar + nparam_1(i) + j)**2 + dupr_1(2*indpar + nparam_1(i) + j + 1)**2
                    if (cost .gt. costmax) then
                        costmax = cost
                        jmax = j
                    end if
                end do
                maxdupr(indpar + nparam_1(i)/2 + (jmax + 1)/2) = .true.
            end if
        elseif (yes_complex .and. nparam_1(i) .eq. 1) then ! complex w.f. with uncontracted orb
            read (10, *) kion_1(i), dupr_1(2*indpar + 1)
            dupr_1(2*indpar + 2) = 0.d0 ! exponents always real
        else
            read (10, *) kion_1(i), (dupr_1(indpar + j), j=1, nparam_1(i))

            if (nparam_1(i) .gt. 1) then
                if (ipf .eq. 2 .and. ioptorb_1(i) .ge. 900000) then
                    costmax = -1.d0
                    do j = 1, nparam_1(i)/4
                        cost = abs(dupr_1(indpar + nparam_1(i)/2 + j))
                        if (cost .gt. costmax) then
                            costmax = cost
                            jmax = j
                        end if
                    end do
                    maxdupr(indpar + nparam_1(i)/2 + jmax) = .true.
                    costmax = -1.d0
                    do j = nparam_1(i)/4 + 1, nparam_1(i)/2
                        cost = abs(dupr_1(indpar + nparam_1(i)/2 + j))
                        if (cost .gt. costmax) then
                            costmax = cost
                            jmax = j
                        end if
                    end do
                    maxdupr(indpar + nparam_1(i)/2 + jmax) = .true.
                else
                    costmax = -1.d0
                    do j = 1, nparam_1(i)/2
                        cost = abs(dupr_1(indpar + nparam_1(i)/2 + j))
                        if (cost .gt. costmax) then
                            costmax = cost
                            jmax = j
                        end if
                    end do
                    maxdupr(indpar + nparam_1(i)/2 + jmax) = .true.
                end if
            end if
        end if
        indpar = indpar + nparam_1(i)
        if (ioptorb_1(i) .ne. 1000000) indpar_at = indpar_at + nparam_1(i)
    end do

    if (yes_hermite .and. ipc .eq. 2) then
        !         Set to zero imaginary part of contracted, but the molecular orbitals
        call dscalzero(indpar_at, 0.d0, dupr_1(2), 2)
    end if

    read (10, *)

    indpar = 0
    do i = 1, nshellj_1
        read (10, *) multj_1(i), nparamj_1(i), ioptorbj_1(i)
        if (ioptorbj_1(i) .eq. 900000) yes9j = .true.
        read (10, *) kionj_1(i), (vjur_1(indpar + j), j=1, nparamj_1(i))
        indpar = indpar + nparamj_1(i)
    end do

    !         count the shell according to the type....

    occ_1 = 0
    do i = 1, nshell_1
        occ_1 = occ_1 + mult_1(i)
    end do

    allocate (ioccup_1(occ_1))

    read (10, *)
    do i = 1, occ_1
        read (10, *) ioccup_1(i)
    end do

    occtot_1 = 0
    do i = 1, occ_1
        if (ioccup_1(i) .ne. 0) occtot_1 = occtot_1 + 1
    end do

    allocate (ioptorba(ipf*occtot_1 + neldiff_1))
    ioptorba = 0
    ind = 0
    indj = 0
    mol_ies = .false.
    nelorb_at = 0
    occtot_1 = 0
    do i = 1, nshell_1
        if (ioptorb_1(i) .eq. 1000000) mol_ies = .true.
        do j = 1, mult_1(i)
            ind = ind + 1
            if (ioccup_1(ind) .eq. 1) then
                indj = indj + 1
                ioptorba(indj) = ioptorb_1(i)
                if (ioptorb_1(i) .ne. 1000000) then
                    nelorb_at = nelorb_at + ipf
                    occtot_1 = occtot_1 + ipf
                else
                    occtot_1 = occtot_1 + 1
                end if
            end if
        end do
    end do
    if (ipf .eq. 2) then
        !         first put in order the molecular orbitals
        if (occtot_1 .gt. nelorb_at) ioptorba(occtot_1 - nelorb_at + 1:occtot_1) = ioptorba(nelorb_at/2 + 1:indj)
        ioptorba(nelorb_at/2 + 1:nelorb_at) = ioptorba(1:nelorb_at/2)
    end if

    !         mol_ies --> .true. if there are molecular orbitals.

    occj_1 = 0
    do i = 1, nshellj_1
        occj_1 = occj_1 + multj_1(i)
    end do

    !         DEALLOCATE(ioccup_1)
    allocate (ioccupj_1(occj_1), ioptorbja_1(occj_1*ipj))

    read (10, *)
    do i = 1, occj_1
        read (10, *) ioccupj_1(i)
    end do

    occj_2 = 0
    do i = 1, nshellj_1
        occj_2 = occj_2 + multj_1(i)
    end do

    ind = 0
    indj = 0
    do i = 1, nshellj_1
        do j = 1, multj_1(i)
            ind = ind + 1
            if (ioccupj_1(ind) .eq. 1) then
                indj = indj + 1
                ioptorbja_1(indj) = ioptorbj_1(i)
            end if
        end do
    end do
    occtotj = 0
    do i = 1, occj_1
        if (ioccupj_1(i) .ne. 0) occtotj = occtotj + 1
    end do
    if (ipj .eq. 2) then
        do i = 1, occtotj
            ioptorbja_1(i + occtotj) = ioptorbja_1(i)
        end do
    end if
    indj = 1

    if (ipf .eq. 2) then

        if (symmagp) then
            nozero_at = ((nelorb_at/2 + 1)*nelorb_at)/2
            nozero_at = nozero_at + ((nelorb_at/2 - 1)*nelorb_at)/2
        else
            nozero_at = ((nelorb_at - 1)*nelorb_at)/2
            if (pfaffup) nozero_at = nozero_at - ((nelorb_at/2 - 1)*nelorb_at/2)/2
        end if

    else
        if (symmagp) then
            nozero_at = ((nelorb_at + 1)*nelorb_at)/2
        else
            nozero_at = nelorb_at*nelorb_at
        end if
    end if
    !------------------------------UMAT by neek
    !            BUG
    !            allocate( lambda_mat(occ_1,occ_1))
    allocate (lambda_mat(ipc*occtot_1, occtot_1 + neldiff_1))
    lambda_mat = 0

    !-----------------------------
    read (10, *)
    mol_off = .false.
    count_at = 0
    do i = 1, nnozero_1
        read (10, *) ix_1(i), iy_1(i), detmat_1(ipc*(i - 1) + 1:ipc*i)
        if (symmagp .and. ipc .eq. 1) then
            if (ix_1(i) .ne. iy_1(i) .and. iy_1(i) .le. occtot_1 .and.&
                    &(ioptorba(ix_1(i)) .eq. 1000000 .or. ioptorba(iy_1(i)) .eq. 1000000)&
                    &.and. detmat_1(i) .ne. 0.d0) mol_off = .true.
        elseif (ipc .eq. 2) then
            if (ix_1(i) .ne. iy_1(i) + 1 .and. iy_1(i) .le. occtot_1 .and.&
                    &(ioptorba(ix_1(i)) .eq. 1000000 .or. ioptorba(iy_1(i)) .eq. 1000000)&
                    &.and. (detmat_1(2*i - 1) .ne. 0.d0 .or. detmat_1(2*i) .ne. 0.d0)) mol_off = .true.
        else
            if (ix_1(i) .ne. iy_1(i) + 1 .and. iy_1(i) .le. occtot_1 .and.&
                    &(ioptorba(ix_1(i)) .eq. 1000000 .or. ioptorba(iy_1(i)) .eq. 1000000)&
                    &.and. detmat_1(i) .ne. 0.d0) mol_off = .true.
        end if
        if (ioptorba(ix_1(i)) .ne. 1000000 .and. ioptorba(iy_1(i)) .ne. 1000000) &
            count_at = count_at + 1
    end do
    if (npar_eagp .ne. 0) then
        allocate (psip(ipc*nnozero_eagp))
        if (ipc .eq. 1) then
            do i = 1, nnozero_eagp
                read (10, *) ix_1(nnozero_1 + i), iy_1(nnozero_1 + i), psip(i)
            end do
            do i = 1, nnozero_eagp
                eagp_pfaff(ix_1(nnozero_1 + i), iy_1(nnozero_1 + i)) = psip(i)
            end do
        else
            do i = 1, nnozero_eagp
                read (10, *) ix_1(nnozero_1 + i), iy_1(nnozero_1 + i), psip(2*i - 1), psip(2*i)
            end do
            do i = 1, nnozero_eagp
                eagp_pfaff(2*ix_1(nnozero_1 + i) - 1, iy_1(nnozero_1 + i)) = psip(2*i - 1)
                eagp_pfaff(2*ix_1(nnozero_1 + i), iy_1(nnozero_1 + i)) = psip(2*i)
            end do
        end if
        deallocate (psip)
    end if
    !    mol_off = .false. --> standard fort.10 with molecular orbitals e.g. output
    !    DFT. mol_off =.true. --> non trivial molecular orbitals matrix elements.
    !    if mol_ies is present.

    if (count_at .ge. nozero_at .and. ipf .eq. 1) yesalldet = .true.

    maxnozero = nnozero_1
    maxnozeroj = nnozeroj_1
    maxnozero = 2*maxnozero

    allocate (ipipsav(2*(nnozero_1 + nnozero_eagp)), icd_1(iesswr_1 + iesswr_eagp))

    maxic = 0
    ind = 0
    read (10, *)
    do i = 1, iesswr_1 + iesswr_eagp
        read (10, *) icd_1(i), (ipipsav(j + ind), j=1, 2*abs(icd_1(i)))
        if (ipc .eq. 2) then
            do j = 2, 2*abs(icd_1(i)), 2
                if (ipipsav(j + ind) .lt. 0) yesfixphase = .true.
            end do
        end if
        if (abs(icd_1(i)) .gt. maxic) maxic = abs(icd_1(i))
        ind = ind + 2*abs(icd_1(i))
    end do
    if (yesfixphase) yesalldet = .false. ! To be safe.

    if (yesalldet) write (6, *) ' Warning changing contracted and tranforming the Det matrix '

    allocate (ipsipd_1(iesswr_1 + iesswr_eagp, 2*maxic))

    ind = 0
    do i = 1, iesswr_1 + iesswr_eagp
        do j = 1, 2*abs(icd_1(i))
            ipsipd_1(i, j) = ipipsav(j + ind)
        end do
        ind = ind + 2*abs(icd_1(i))
    end do

    deallocate (ipipsav)

    allocate (ixj_1(nnozeroj_1), iyj_1(nnozeroj_1)               &
            &, jasmat_1(nnozeroj_1), jasmatsz_1(nnozeroj_1))

    scale = dble(nel_1 - 1)/dble(nel - 1)

    !     new lambdas
    allocate (lambda_matj(ipj*occtotj, ipj*occtotj))
    allocate (lambda_matszj(occtotj, occtotj))
    lambda_matj = 0.d0
    lambda_matszj = 0.d0
    jasmat_1 = 0.d0
    jasmatsz_1 = 0.d0

    read (10, *)
    molj_off = .false.
    do i = 1, nnozeroj_1
        read (10, *) ixj_1(i), iyj_1(i), jasmat_1(i)
        lambda_matj(ixj_1(i), iyj_1(i)) = jasmat_1(i)
        lambda_matj(iyj_1(i), ixj_1(i)) = jasmat_1(i)
        if (ixj_1(i) .ne. iyj_1(i) .and. (ioptorbja_1(ixj_1(i)) .eq. 1000000&
                &.or. ioptorbja_1(iyj_1(i)) .eq. 1000000)) molj_off = .true.
    end do
    if (iessz) then
        read (10, *)
        do i = 1, nnozeroj_1
            read (10, *) ixj_1(i), iyj_1(i), jasmatsz_1(i)
            lambda_matszj(ixj_1(i), iyj_1(i)) = jasmatsz_1(i)
            lambda_matszj(iyj_1(i), ixj_1(i)) = jasmatsz_1(i)
        end do
    end if

    allocate (icj_1(iesfreer_1))

    allocate (ipipsav(2*nnozeroj_1))

    maxic = 1
    ind = 0
    read (10, *)
    do i = 1, iesfreer_1
        !            read(10,*) icj_1(i),(ipsipj_1(i,j),j=1,2*abs(icj_1(i)))
        read (10, *) icj_1(i), (ipipsav(j + ind), j=1, 2*abs(icj_1(i)))
        ind = ind + 2*abs(icj_1(i))
        if (abs(icj_1(i)) .gt. maxic) maxic = abs(icj_1(i))
    end do
    allocate (ipsipj_1(max(iesfreer_1, 1), 2*maxic))
    ind = 0
    do i = 1, iesfreer_1
        do j = 1, 2*abs(icj_1(i))
            ipsipj_1(i, j) = ipipsav(j + ind)
        end do
        ind = ind + 2*abs(icj_1(i))
    end do

    deallocate (ipipsav)

    allocate (cek_w(iesupr_1))

    !            turn on  the constant

    allocate (cek_wj(npar3body_1))

    cek_wj = 0
    cek_w = 0

    read (10, *)

    maxic = 0
    ind = 0
    allocate (ipipsav(iesupr_1), icp_1(iesupind_1))

    do i = 1, iesupind_1
        read (10, *) icp_1(i), (ipipsav(ind + j), j=1, abs(icp_1(i)))
        icost = abs(icp_1(i))
        if (icost .gt. maxic) maxic = icost
        icp_1(i) = abs(icp_1(i))
        if (.not. yesalldet) then
            check = .false.
            do j = 1, abs(icp_1(i))
                if (maxdupr(abs(ipipsav(ind + j)))) check = .true.
            end do
            if (check) icp_1(i) = -icp_1(i)
        end if
        ind = ind + icost
    end do

    allocate (ipsipp_1(iesupind_1, maxic))

    signpar = .true.
    ind = 0
    do i = 1, iesupind_1
        icost = abs(icp_1(i))
        do j = 1, icost
            ipsipp_1(i, j) = ipipsav(ind + j)
            if (ipsipp_1(i, j) .lt. 0) signpar(-ipsipp_1(i, j)) = .false.
        end do
        do j = 1, abs(icp_1(i))
            cek_w(abs(ipsipp_1(i, j))) = i
        end do
        ind = ind + icost
    end do
    deallocate (ipipsav)

    !           write(6,*) ' signpar determined '
    !           do i=1,iesupr_1
    !           write(6,*) i,signpar(i)
    !           enddo

    allocate (icq_1(iesmind_1), ipipsav(npar3body_1))
    read (10, *)
    ind = 0
    maxic = 0
    do i = 1, iesmind_1
        read (10, *) icq_1(i), (ipipsav(ind + j), j=1, abs(icq_1(i)))
        icost = abs(icq_1(i))
        icq_1(i) = icost
        ind = ind + icost
        if (icost .gt. maxic) maxic = icost
    end do

    allocate (ipsipq_1(iesmind_1, maxic))
    signparj = .true.
    ind = 0
    do i = 1, iesmind_1
        do j = 1, abs(icq_1(i))
            ipsipq_1(i, j) = ipipsav(ind + j)
            if (ipsipq_1(i, j) .lt. 0) signparj(-ipsipq_1(i, j)) = .false.
        end do
        do j = 1, abs(icq_1(i))
            cek_wj(abs(ipsipq_1(i, j))) = i
        end do
        ind = ind + icq_1(i)
    end do

    deallocate (ipipsav)

    !   ################################################## Counting AGP...
    !   start here

    !------------------------------------   iesupind_new
    !                nshell_1>number of each shell so :

    allocate (ns(nshell_1), address(nshell_1, nshell_1))
    allocate (halfp(nshell_1), cek_occ(nshell_1)                        &
            &, addressorb(nshell_1, nshell_1), addressorbtot(nshell_1, nshell_1))
    allocate (mult_new(nshell_1))

    cek_occ = .false.

    k = 0
    ns = 0
    indtot0 = 0
    indorb0 = 0
    ind0 = 0

    !      enumeration shell to be touched

    do i = 1, nshell_1

        if (.not. cek_occ(i) .and. (nparam_1(i) .ge. 2 .and.&
                &(mol_off .or. ioptorb_1(i) .ne. 1000000))) then

            k = k + 1

            mult_new(k) = mult_1(i)
            halfp(k) = nparam_1(i)/2
            ind = 0
            indorb = 0
            indtot = 0
            do j = 1, nshell_1

                if ((kion_1(i) .eq. kion_1(j) .or. ioptorb_1(i) .ge. 1000000) .and.       &
                        &          ioptorb_1(i) .eq. ioptorb_1(j)                            &
                        &      .and. nparam_1(i) .eq. nparam_1(j)) then

                    cek_occ(j) = .true.

                    !   check also that the occupations are exactly the same

                    do kk = 1, mult_1(i)
                        if (ioccup_1(indtot + kk) .ne. ioccup_1(indtot0 + kk)) cek_occ(j)        &
                                & = .false.
                    end do

                    !  check also that also all the exponent Z are the same and in the
                    !  same order

                    do kk = 1, halfp(k)
                        if (dupr_1(ipc*(ind + kk - 1) + 1) .ne. dupr_1(ipc*(ind0 + kk - 1) + 1)) cek_occ(j) = .false.
                    end do

                    if (cek_occ(j)) then
                        !        Update addresses
                        ns(k) = ns(k) + 1
                        address(ns(k), k) = ind
                        addressorb(ns(k), k) = indorb
                        addressorbtot(ns(k), k) = indtot
                    end if
                end if
                do kk = 1, mult_1(j)
                    if (ioccup_1(indtot + kk) .eq. 1) indorb = indorb + 1
                end do
                indtot = indtot + mult_1(j)
                ind = ind + nparam_1(j)
                ! enddo j
            end do

            ! endif occupation
        end if

        do kk = 1, mult_1(i)
            if (ioccup_1(indtot0 + kk) .eq. 1) indorb0 = indorb0 + 1
        end do

        indtot0 = indtot0 + mult_1(i)

        ind0 = ind0 + nparam_1(i)

        ! enddo i
    end do

    !      write(6,*) ' number of contracted groups found =',k
    !      do i=1,k
    !       write(6,*) ' multeplicity shell =',ns(i)
    !       write(6,*) ' multeplicity type  =',mult_new(i)
    !
    !       write(6,*) ' Address found '
    !       do j=1,ns(i)
    !       write(6,*) j,address(j,i),addressorb(j,i),addressorbtot(j,i)
    !       enddo
    !
    !      enddo
    !
    !      stop

    !ccccccccccccccccccccccccccccc UMAT part
    !ccccccccccccccccccccccccccccccccccccccccccc
    allocate (UU(ipc*occtot_1, occtot_1))

    UU = 0.d0
    do i = 1, occtot_1
        UU(ipc*(i - 1) + 1, i) = 1.d0
    end do

    if (yesalldet) then

        do i = 1, k

            allocate (psi_in(ipc*halfp(i), ns(i)))
            allocate (psi_out(ipc*halfp(i), ns(i)))
            allocate (ipiv(ns(i)))
            allocate (umat(ipc*ns(i), ns(i)))

            do kk = 1, ns(i)
                do j = 1, halfp(i)
                    psi_in(ipc*(j - 1) + 1:ipc*j, kk) = dupr_1(ipc*(address(kk, i) + j + halfp(i) - 1) &
                                                               + 1:ipc*(address(kk, i) + j + halfp(i)))
                end do
            end do

            !     write(6,*) ' Input ns(i) =',ns(i)

            if (ipc .eq. 1) then
                call independent(ns(i), halfp(i), psi_in, psi_out, ipiv)
            else
                call independent_complex(ns(i), halfp(i), psi_in, psi_out, ipiv, yesfixphase)
            end if
            !     choose the gauge
            do kk = 1, ns(i)
                if (.not. signpar(address(kk, i) + ipiv(kk) + halfp(i))) then
                    psi_out(1:ipc*halfp(i), kk) = -psi_out(1:ipc*halfp(i), kk)
                end if
            end do

            do kk = 1, ns(i)
                do j = 1, halfp(i)
                    dupr_1(ipc*(address(kk, i) + j + halfp(i) - 1) + 1:ipc*(address(kk, i) + j + halfp(i))) &
                        = psi_out(ipc*(j - 1) + 1:ipc*j, kk)
                end do
            end do
            !---------------------------Umat
            !       write(*,*)ns(i),halfp(i)
            !      makeumat(n,ncoeff,psi_in,psi_out,umat)
            !      error is here !

            call makeumat(ipc, ns(i), halfp(i), psi_in, psi_out, umat)
            !      write(6,*) ' Transformation found ',ns(i)
            !      do ii=1,ns(i)
            !         do jj=1,ns(i)
            !         write(6,*) ii,jj,umat(ipc*(ii-1)+1:ipc*ii,jj)
            !         enddo
            !      enddo

            !       write(6,*) ' full matrix ',ns(i),mult_new(i)
            do ii = 1, ns(i)
                do jj = 1, ns(i)
                    indorbi = addressorb(ii, i)
                    indorbj = addressorb(jj, i)
                    do kk = 1, mult_new(i)
                        !          write(6,*) ' addressorbtot =',addressorbtot(jj,ii)
                        if (ioccup_1(addressorbtot(jj, i) + kk) .eq. 1) then
                            indorbj = indorbj + 1
                            indorbi = indorbi + 1
                            UU(ipc*(indorbi - 1) + 1:ipc*indorbi, indorbj) = umat(ipc*(ii - 1) + 1:ipc*ii, jj)
                            !          write(6,*) indorbi,indorbj,umat(ii,jj)
                        end if
                    end do
                end do
            end do

            do kk = 1, ns(i)
                do jj = 1, ns(i)

                    !      write(6,*) kk,address(jj,i),halfp(i),ipiv(kk)
                    if (cek_w(address(jj, i) + halfp(i) + ipiv(kk)) .ne. 0) then
                        icp_1(cek_w(address(jj, i) + halfp(i) + ipiv(kk))) = &
                                &     -abs(icp_1(cek_w(address(jj, i) + halfp(i) + ipiv(kk))))
                        !      else
                        !      icp_1(cek_w(address(jj,i)+halfp(i)+ipiv(kk)))=                   &
                        !    &     abs(icp_1(cek_w(address(jj,i)+halfp(i)+ipiv(kk))))
                    end if
                end do
            end do

            deallocate (psi_in)
            deallocate (psi_out)
            deallocate (ipiv)
            deallocate (umat)

        end do

    end if ! endif yesalldet
    !          write(6,*) ' Final umat before '
    !          do i=1,occ_1
    !            do j=1,occ_1
    !            write(6,*) i,j,uu(ipc*(i-1)+1:ipc*i,j)
    !            enddo
    !          enddo

    !         lambda has to be symmetrized  only for symmagp=.true.

    ndimh = nelorb_at/2

    do i = 1, nnozero_1
        lambda_mat(ipc*(ix_1(i) - 1) + 1:ipc*ix_1(i), iy_1(i)) = detmat_1(ipc*(i - 1) + 1:ipc*i)
        if (ipf .eq. 2 .and. iy_1(i) .le. occtot_1) then
            lambda_mat(ipc*(iy_1(i) - 1) + 1:ipc*iy_1(i), ix_1(i)) = -detmat_1(ipc*(i - 1) + 1:ipc*i)
        end if
        if (symmagp .and. ipf .eq. 2 .and. iy_1(i) .le. nelorb_at .and. ix_1(i) .le. nelorb_at) then
            if (ipc .eq. 1) then
                if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                    lambda_mat(ix_1(i) + ndimh, iy_1(i) + ndimh) = detmat_1(i)
                    lambda_mat(iy_1(i) + ndimh, ix_1(i) + ndimh) = -detmat_1(i)
                elseif (iy_1(i) .gt. ndimh .and. iy_1(i) .le. nelorb_at .and. ix_1(i) .le. ndimh) then
                    lambda_mat(iy_1(i) - ndimh, ix_1(i) + ndimh) = detmat_1(i)
                    lambda_mat(ix_1(i) + ndimh, iy_1(i) - ndimh) = -detmat_1(i)
                end if
            else
                if (yes_hermite) then
                    if (iy_1(i) .le. ndimh .and. ix_1(i) .le. ndimh .and. .not. pfaffup) then
                        lambda_mat(2*(ix_1(i) + ndimh) - 1, iy_1(i) + ndimh) = detmat_1(2*i - 1)
                        lambda_mat(2*(ix_1(i) + ndimh), iy_1(i) + ndimh) = -detmat_1(2*i)
                        lambda_mat(2*(iy_1(i) + ndimh) - 1, ix_1(i) + ndimh) = -detmat_1(2*i - 1)
                        lambda_mat(2*(iy_1(i) + ndimh), ix_1(i) + ndimh) = detmat_1(2*i)
                    elseif (iy_1(i) .gt. ndimh .and. ix_1(i) .le. ndimh) then
                        lambda_mat(2*(iy_1(i) - ndimh) - 1, ix_1(i) + ndimh) = detmat_1(2*i - 1)
                        lambda_mat(2*(iy_1(i) - ndimh), ix_1(i) + ndimh) = -detmat_1(2*i)
                        lambda_mat(2*(ix_1(i) + ndimh) - 1, iy_1(i) - ndimh) = -detmat_1(2*i - 1)
                        lambda_mat(2*(ix_1(i) + ndimh), iy_1(i) - ndimh) = detmat_1(2*i)
                    end if
                else
                    if (iy_1(i) .le. ndimh .and. ix_1(i) .le. ndimh .and. .not. pfaffup) then
                        lambda_mat(2*(ix_1(i) + ndimh) - 1:2*(ix_1(i) + ndimh), iy_1(i) + ndimh) = detmat_1(2*i - 1:2*i)
                        lambda_mat(2*(iy_1(i) + ndimh) - 1:2*(iy_1(i) + ndimh), ix_1(i) + ndimh) = -detmat_1(2*i - 1:2*i)
                    elseif (iy_1(i) .gt. ndimh .and. ix_1(i) .le. ndimh) then
                        lambda_mat(2*(iy_1(i) - ndimh) - 1:2*(iy_1(i) - ndimh), ix_1(i) + ndimh) = detmat_1(2*i - 1:2*i)
                        lambda_mat(2*(ix_1(i) + ndimh) - 1:2*(ix_1(i) + ndimh), iy_1(i) - ndimh) = -detmat_1(2*i - 1:2*i)
                    end if
                end if
            end if
        end if
        if (ipf .eq. 1 .and. iy_1(i) .le. nelorb_at .and. ix_1(i) .le. nelorb_at .and. symmagp) then
            if (yes_hermite .and. ipc .eq. 2) then
                lambda_mat(2*iy_1(i) - 1, ix_1(i)) = detmat_1(2*i - 1)
                lambda_mat(2*iy_1(i), ix_1(i)) = -detmat_1(2*i)
            else
                lambda_mat(ipc*(iy_1(i) - 1) + 1:ipc*iy_1(i), ix_1(i)) = detmat_1(ipc*(i - 1) + 1:ipc*i)
            end if
        end if
    end do

    call lambda_transf(ipc, occtot_1, occtot_1 + neldiff_1, lambda_mat, uu)

    !         copy lambda_mat in detmat_1
    !          write(6,*) ' input detmat '
    do i = 1, nnozero_1
        detmat_1(ipc*(i - 1) + 1:ipc*i) = lambda_mat(ipc*(ix_1(i) - 1) + 1:ipc*ix_1(i), iy_1(i))
        !          write(6,*) i,detmat_1(ipc*(i-1)+1:ipc*i)
    end do

    if (ipc .eq. 1 .and. symmagp .and. ipf .eq. 1) then
        allocate (eig(occtot_1), psip(3*occtot_1))
        call dsyev('N', 'L', occtot_1, lambda_mat, occtot_1, eig, psip, 3*occtot_1, info)
    else
        if (yes_complex) then
            allocate (eig(occtot_1), psip(11*occtot_1))
            call zgesvd('N', 'N', occtot_1, occtot_1, lambda_mat, occtot_1, eig, u, 1, v, 1, psip, 3*occtot_1&
                    &, psip(6*occtot_1 + 1), info)
            !      write(6,*) ' output eigenvalues '
            !      do i=1,occtot_1
            !      write(6,*) i,eig(i)
            !      enddo
            !      stop
        else
            allocate (eig(occtot_1), psip(5*occtot_1))
            call dgesvd('N', 'N', occtot_1, occtot_1, lambda_mat, occtot_1, eig, u, 1, v, 1, psip, 5*occtot_1, info)
        end if
    end if

    irankdet = 0
    do i = 1, occtot_1
        if (abs(eig(i)) .gt. 1d-11) then
            irankdet = irankdet + 1
        end if
    end do

    if (ipf .eq. 2) then
        if (irankdet .eq. neldo_1 + nelup_1) then
            if (occtot_1 .eq. neldo_1 + nelup_1) then
                iesdet = .true.
            else
                iesdet = .false.
            end if
        elseif (irankdet .lt. neldo_1 + nelup_1) then
            write (6, *) ' Singular fort.10 !!!! '
            iesdet = .false.
        else
            iesdet = .false.
        end if
    else
        if (irankdet .eq. neldo_1) then
            if (occtot_1 .eq. neldo_1) then
                iesdet = .true.
            else
                iesdet = .false.
            end if
        elseif (irankdet .lt. neldo_1) then
            write (6, *) ' Singular fort.10 !!!! '
            iesdet = .false.
        else
            iesdet = .false.
        end if
    end if

    deallocate (eig, psip)

    !         write(6,*) ' output  lambda_mat '
    !         do j=1,occ_1+neldiff_1
    !           write(6,*) ' raw ',j
    !           do i=1,occ_1
    !           write(6,*) i,lambda_mat(i,j)
    !           enddo
    !         enddo

    deallocate (UU, address, addressorb, addressorbtot, cek_occ)

    !cccccccccccccccccccccccccccccccccccccccccccccccccc
    indpar = 0
    k = 0
    !                nshell_1>number of each shell so :

    allocate (nsj(nshellj_1))
    allocate (halfpj(nshellj_1), cek_occ(nshellj_1)                     &
            &, address(nshellj_1, nshellj_1), addressorb(nshellj_1, nshellj_1)     &
            &, addressorbtot(nshellj_1, nshellj_1))
    allocate (mult_newj(nshellj_1))

    !      enumeration shell to be touched

    cek_occ = .false.

    k = 0
    nsj = 0
    indtot0 = 0
    indorb0 = 0
    ind0 = 0

    do i = 1, nshellj_1

        if (.not. cek_occ(i) .and. (nparamj_1(i) .ge. 2 .and. (molj_off .or. ioptorbj_1(i) .ne. 1000000))) then
            k = k + 1
            mult_newj(k) = multj_1(i)
            halfpj(k) = nparamj_1(i)/2
            ind = 0
            indorb = 0
            indtot = 0
            do j = 1, nshellj_1

                if ((kionj_1(i) .eq. kionj_1(j) .or. ioptorbj_1(i) .ge. 1000000) .and.&
                        &ioptorbj_1(i) .eq. ioptorbj_1(j) .and. nparamj_1(i) .eq. nparamj_1(j)) then

                    cek_occ(j) = .true.

                    !   check also that the occupations are exactly the same

                    do kk = 1, multj_1(i)
                        if (ioccupj_1(indtot + kk) .ne. ioccupj_1(indtot0 + kk))                 &
                                &  cek_occ(j) = .false.
                    end do

                    !  check also that also all the exponent Z are the same and in the
                    !  same order

                    do kk = 1, halfpj(k)
                        if (vjur_1(ind + kk) .ne. vjur_1(ind0 + kk)) cek_occ(j) = .false.
                    end do

                    if (cek_occ(j)) then
                        !        Update addresses
                        nsj(k) = nsj(k) + 1
                        address(nsj(k), k) = ind
                        addressorb(nsj(k), k) = indorb
                        addressorbtot(nsj(k), k) = indtot
                    end if
                end if
                do kk = 1, multj_1(j)
                    if (ioccupj_1(indtot + kk) .eq. 1) indorb = indorb + 1
                end do
                indtot = indtot + multj_1(j)
                ind = ind + nparamj_1(j)
                ! enddo j
            end do

            ! endif occupation
        end if

        do kk = 1, multj_1(i)
            if (ioccupj_1(indtot0 + kk) .eq. 1) indorb0 = indorb0 + 1
        end do

        indtot0 = indtot0 + multj_1(i)
        ind0 = ind0 + nparamj_1(i)

        ! enddo i
    end do

    !      write(6,*) ' number of contracted groups found =',k
    !      do i=1,k
    !       write(6,*) ' multeplicity shell =',ns(i)
    !       write(6,*) ' multeplicity type  =',mult_new(i)

    !       write(6,*) ' Address found '
    !       do j=1,ns(i)
    !       write(6,*) j,address(j,i),addressorb(j,i),addressorbtot(j,i)
    !       enddo
    !      enddo

    allocate (UU(ipj*occtotj, ipj*occtotj))

    !      The identity for all untouched orbitals

    UU = 0.d0
    do i = 1, ipj*occtotj
        UU(i, i) = 1.d0
    end do

    !     write(6,*) ' before loop jas '
    !.......................................
    do i = 1, k

        allocate (psi_in(halfpj(i), nsj(i)))
        allocate (psi_out(halfpj(i), nsj(i)))
        allocate (ipiv(nsj(i)))
        allocate (umatj(nsj(i), nsj(i)))

        do kk = 1, nsj(i)
            do j = 1, halfpj(i)
                psi_in(j, kk) = vjur_1(address(kk, i) + j + halfpj(i))
            end do
        end do

        call independent(nsj(i), halfpj(i), psi_in, psi_out, ipiv)

        !   write(6,*) choose the gauge
        do kk = 1, nsj(i)
            if (.not. signparj(address(kk, i) + ipiv(kk) + halfpj(i))) psi_out(1:halfpj(i), kk) = -psi_out(1:halfpj(i), kk)
        end do

        do kk = 1, nsj(i)
            do j = 1, halfpj(i)
                vjur_1(address(kk, i) + j + halfpj(i)) = psi_out(j, kk)
            end do
        end do
        !              write(*,*)'mehdiam',i,nsj(i)

        !cccccccccccccc

        call makeumat(1, nsj(i), halfpj(i), psi_in, psi_out, umatj)

        do ii = 1, nsj(i)
            do jj = 1, nsj(i)
                indorbj = addressorb(jj, i)
                indorbi = addressorb(ii, i)
                do kk = 1, mult_newj(i)
                    if (ioccupj_1(addressorbtot(jj, i) + kk) .eq. 1) then
                        indorbj = indorbj + 1
                        indorbi = indorbi + 1
                        UU(indorbi, indorbj) = umatj(ii, jj)
                    end if
                end do
            end do
        end do

        do kk = 1, nsj(i)
            do jj = 1, nsj(i)
                if (cek_wj(address(jj, i) + halfpj(i) + ipiv(kk)) .ne. 0) then
                    icq_1(cek_wj(address(jj, i) + halfpj(i) + ipiv(kk))) = &
                            &        -abs(icq_1(cek_wj(address(jj, i) + halfpj(i) + ipiv(kk))))
                    !         else
                    !      icq_1(cek_wj(address(jj,i)+halfpj(i)+ipiv(kk)))=                 &
                    !    &        abs(icq_1(cek_wj(address(jj,i)+halfpj(i)+ipiv(kk))))
                end if
            end do
        end do

        deallocate (psi_in)
        deallocate (psi_out)
        deallocate (ipiv)
        deallocate (umatj)

        ! do i=1,k
    end do
    if (ipj .eq. 2) then
        do i = 1, occtotj
            do j = 1, occtotj
                UU(i + occtotj, j + occtotj) = UU(i, j)
            end do
        end do
    end if

    if (occj_1 .gt. 0) then
        call lambda_transf(1, ipj*occtotj, ipj*occtotj, lambda_matj, uu)
        do i = 1, nnozeroj_1
            jasmat_1(i) = lambda_matj(ixj_1(i), iyj_1(i))
        end do

        if (iessz) then
            call lambda_transf(1, occtotj, occtotj, lambda_matszj, uu)
            do i = 1, nnozeroj_1
                jasmatsz_1(i) = lambda_matszj(ixj_1(i), iyj_1(i))
            end do
        end if

    end if

    deallocate (uu)

    if (iespbc) then
        if (yes_crystal) then
            if (yes_tilted) then
                write (11, *) '# PBC_T a b c phase up phase down '
                write (11, '(9f18.12,6f12.8)') celldm1(1:9), phase(:), phase_down(:)
            else
                write (11, *) '# PBC_C rs, Ly/Lx, Lz/Lx '
                write (11, '(3f18.12,6f12.8)') rs1, Ly1, lz1, phase(:), phase_down(:)
            end if
        else
            write (11, *) '# PBC rs, Ly/Lx, Lz/Lx '
            write (11, '(3f18.12,3f12.8)') rs1, Ly1, lz1, phase(:)
        end if
    end if

    write (11, *) '# Nelup  #Nel  # Ion '
    nelup_1 = nelup_1 + numpaired/2*nel_1
    if (pfaffup) nelup_1 = -nelup_1
    if (ipf .eq. 2) nel_1 = -nel_1
    if (yesbump) then
        write (11, *) nelup_1, nel_1, -nion_1
    else
        write (11, *) nelup_1, nel_1, nion_1
    end if
    if (ipf .eq. 2) nel_1 = -nel_1
    if (pfaffup) nelup_1 = -nelup_1
    write (11, *) '# Shell Det.   # Shell Jas. '
    if (yes_crystalj) nshellj_1 = -nshellj_1
    if (forcecomplex) then
        write (11, *) - nshell_2, nshellj_1
    else
        write (11, *) nshell_2, nshellj_1
    end if
    if (yes_crystalj) nshellj_1 = -nshellj_1

    if (ipf .eq. 2) then
        neldiff = mod(abs(nel_1), 2) + numpaired
    else
        neldiff = nelup_1 - (nel_1 - nelup_1)
    end if
    !         calculation of iesupr_2

    iesupr_2 = 0
    do i = 1, nshell_1
        iesupr_2 = iesupr_2 + nparam_1(i)
    end do

    if (chosen_map) then
        if (iesdrr .ge. 0) then
            iesdrr = iesdrr + 100*case_map
        else
            iesdrr = iesdrr - 100*case_map
        end if
    end if

    if (npar_eagp .ne. 0) then
        write (11, *) '# Jas 2body  # Det   #  3 body atomic par.  iessw_eagp #eagp '
        write (11, *) iesdrr, iesupr_1, npar3body_1, iesswr_eagp, nnozero_eagp
    else
        write (11, *) '# Jas 2body  # Det   #  3 body atomic par.  '
        write (11, *) iesdrr, iesupr_1, npar3body_1
    end if

    write (11, *) '# Det mat. =/0  # Jas mat. =/0  '
    if (forcesymm) then
        write (11, *) - nnozero_1, nnozeroj_1
    else
        write (11, *) nnozero_1, nnozeroj_1
    end if

    write (11, *) ' # Eq. Det atomic par. # Eq. 3 body atomic. par. '

    write (11, *) iesupind_1, iesmind_1

    !          we dont count jastrow for orthogonality
    !          write(11,*)  iesupind_new,iesmind_new

    write (11, *) '# unconstrained iesfree,iessw,ieskinr,I/O flag '
    write (11, *) iesfreer_1, iesswr_1, ieskinr_1, iflag_1

    write (11, *) '# Ion coordinates '
    do j = 1, nion_1
        write (11, *) zetar_1(j), atomic_number(j), (rion_1(i, j), i=1, 3)
    end do

    write (11, *) '#  Constraints for forces: ion - coordinate'
    do i = 1, ieskinr_1
        write (11, *) ic_1(i), (ipsip_1(i, j), j=1, 2*abs(ic_1(i)))
    end do

    write (11, *) '#          Parameters Jastrow two body'
    if (iesdr_1 .eq. 0) then
        if (symmagp) then
            write (11, *) 0
        else
            write (11, *) - 1
        end if
    elseif (abs(iesdr_1) .le. 3) then
        write (11, *) (vj(i), i=1, abs(iesdr_1))
    else
        if (symmagp) then
            write (11, *) abs(niesd), (vj(i), i=1, abs(niesd))
        else
            write (11, *) - abs(niesd), (vj(i), i=1, abs(niesd))
        end if
    end if

    !         if(abs(iesdr_1).ge.4) then
    !         write(11,*) niesd,(vj(i),i=1,abs(niesd))
    !         else
    !         write(11,*) (vj(i),i=1,abs(iesdr_1))
    !         endif

    write (11, *) '#          Parameters atomic wf'
    !         symmetrization  par. before writing it

    do i = 1, iesupind_1
        costc = 0.d0
        do j = 1, abs(icp_1(i))
            if (ipsipp_1(i, j) .gt. 0) then
                costc(1:ipc) = costc(1:ipc) + dupr_1(ipc*(ipsipp_1(i, j) - 1) + 1:ipc*ipsipp_1(i, j))
            elseif (ipsipp_1(i, j) .lt. 0) then
                costc(1:ipc) = costc(1:ipc) - dupr_1(ipc*(-ipsipp_1(i, j) - 1) + 1:ipc*(-ipsipp_1(i, j)))
            end if
        end do
        costc(1:ipc) = costc(1:ipc)/abs(icp_1(i))
        do j = 1, abs(icp_1(i))
            if (ipsipp_1(i, j) .gt. 0) then
                dupr_1(ipc*(ipsipp_1(i, j) - 1) + 1:ipc*ipsipp_1(i, j)) = costc(1:ipc)
            else
                dupr_1(ipc*(-ipsipp_1(i, j) - 1) + 1:ipc*(-ipsipp_1(i, j))) = -costc(1:ipc)
            end if
        end do
    end do

    indpar = 0
    do i = 1, nshell_1

        write (11, *) mult_1(i), nparam_1(i), ioptorb_1(i)
        if (ioptorb_1(i) .eq. 1000000 .or. ioptorb_1(i) .eq. 900000) then
            write (11, *) kion_1(i)                                        &
                    &, (nint(dupr_1(ipc*(indpar + j - 1) + 1)), j=1, nparam_1(i)/2)                       &
                    &, (dupr_1(ipc*(indpar + j - 1) + 1:ipc*(indpar + j)), j=nparam_1(i)/2 + 1, nparam_1(i))
        else
            if (yes_complex) then
                if (nparam_1(i) .gt. 1) then
                    write (11, *) kion_1(i), (dupr_1(2*(indpar + j - 1) + 1), j=1, nparam_1(i)/2)&
                            &, (dupr_1(nparam_1(i) + 2*indpar + j), j=1, nparam_1(i))
                else
                    write (11, *) kion_1(i), dupr_1(2*indpar + 1)
                end if
            else
                write (11, *) kion_1(i), (dupr_1(indpar + j), j=1, nparam_1(i))
            end if
        end if
        indpar = indpar + nparam_1(i)
    end do

    write (11, *) '#  Parameters atomic Jastrow wf '

    !         symmetrization parameters according to sym table
    !            do i=1,iesmind_1
    !          write(11,*) icq_1(i),(ipsipq_1(i,j),j=1,abs(icq_1(i)))
    do i = 1, iesmind_1
        cost = 0.d0
        do j = 1, abs(icq_1(i))
            if (ipsipq_1(i, j) .gt. 0) then
                cost = cost + vjur_1(ipsipq_1(i, j))
            elseif (ipsipq_1(i, j) .lt. 0) then
                cost = cost - vjur_1(-ipsipq_1(i, j))
            end if
        end do
        cost = cost/abs(icq_1(i))
        do j = 1, abs(icq_1(i))
            if (ipsipq_1(i, j) .gt. 0) then
                vjur_1(ipsipq_1(i, j)) = cost
            else
                vjur_1(-ipsipq_1(i, j)) = -cost
            end if
        end do
    end do

    indpar = 0

    do i = 1, nshellj_1
        write (11, *) multj_1(i), nparamj_1(i), ioptorbj_1(i)
        if (ioptorbj_1(i) .eq. 1000000 .or. ioptorbj_1(i) .eq. 900000) then
            write (11, *) kionj_1(i)                                       &
                    &, (nint(vjur_1(indpar + j)), j=1, nparamj_1(i)/2)                      &
                    &, (vjur_1(indpar + j), j=nparamj_1(i)/2 + 1, nparamj_1(i))
        else
            write (11, *) kionj_1(i), (vjur_1(indpar + j), j=1, nparamj_1(i))
        end if
        !           endif
        !          write(11,*) kionj_1(i),(vjur_1(indpar+j),j=1,nparamj_1(i))
        indpar = indpar + nparamj_1(i)
    end do
    !cccccccccccccccccccccccccccccccccccccccc

    write (11, *) '#  Occupation atomic orbitals  '
    norb = 0
    do i = 1, occ_1
        write (11, *) ioccup_1(i)
        norb = norb + ioccup_1(i)
    end do
    norb = norb*ipf

    allocate (scaleorb(norb + neldiff), done(norb + neldiff_1))
    ind = 0
    indorb = 0

    scaleorb = 0.d0

    write (11, *) '#  Occupation atomic orbitals  Jastrow  '
    do i = 1, occj_1
        write (11, *) ioccupj_1(i)
    end do

    norbj = 0
    do i = 1, occj_1
        norbj = norbj + ioccupj_1(i)
    end do

    ind = 0
    indorb = 0

    !            find the maximum detmax
    !            to normalize the maximum geminal det to 1
    !          for each unpaired orbital scale the
    !          orbital in a way that the maximum component is set to 1
    !           symmetrize detmat_1 according to input

    allocate (sortvect(max(nnozero_1, nnozero_eagp))&
            &, indexo(max(nnozero_1, nnozero_eagp)))

    do kk = 1, nnozero_1
        if (ix_1(kk) .le. iy_1(kk) .or. .not. symmagp) then
            sortvect(kk) = (iy_1(kk) - 1)*norb + ix_1(kk)
        else
            sortvect(kk) = (ix_1(kk) - 1)*norb + iy_1(kk)
        end if
    end do

    call dsortx(sortvect, 1, nnozero_1, indexo)

    do i = 1, iesswr_1
        if (abs(icd_1(i)) .ne. 1) then
            costc = 0.d0
            do k = 1, abs(icd_1(i))
                ixt = abs(ipsipd_1(i, 2*k - 1))
                iyt = abs(ipsipd_1(i, 2*k))
                if (ixt .le. iyt .or. .not. symmagp) then
                    rtarget = (iyt - 1)*norb + ixt
                else
                    rtarget = (ixt - 1)*norb + iyt
                end if
                call search(rtarget, sortvect, nnozero_1, kk)
                inddet = indexo(kk)
                if (ipsipd_1(i, 2*k - 1) .lt. 0) then
                    if (ipc .eq. 2 .and. yes_crystal) then
                        costc(2) = costc(2) - detmat_1(2*inddet)
                        costc(1) = costc(1) + detmat_1(2*inddet - 1)
                    else
                        costc(1:ipc) = costc(1:ipc) - detmat_1(ipc*(inddet - 1) + 1:ipc*inddet)
                    end if
                else
                    costc(1:ipc) = costc(1:ipc) + detmat_1(ipc*(inddet - 1) + 1:ipc*inddet)
                end if
            end do

            costc(1:ipc) = costc(1:ipc)/abs(icd_1(i))
            do k = 1, abs(icd_1(i))
                ixt = abs(ipsipd_1(i, 2*k - 1))
                iyt = abs(ipsipd_1(i, 2*k))
                if (ixt .le. iyt .or. .not. symmagp) then
                    rtarget = (iyt - 1)*norb + ixt
                else
                    rtarget = (ixt - 1)*norb + iyt
                end if
                call search(rtarget, sortvect, nnozero_1, kk)
                inddet = indexo(kk)
                if (ipsipd_1(i, 2*k - 1) .lt. 0) then
                    if (yes_crystal .and. ipc .eq. 2) then
                        detmat_1(2*inddet - 1) = costc(1)
                        if (ipsipd_1(i, 2*k) .lt. 0) then
                            detmat_1(2*inddet) = 0.d0
                        else
                            detmat_1(2*inddet) = -costc(2)
                        end if
                    else
                        detmat_1(ipc*(inddet - 1) + 1:ipc*inddet) = -costc(1:ipc)
                    end if
                else
                    if (yes_crystal .and. ipc .eq. 2 .and. ipsipd_1(i, 2*k) .lt. 0) then
                        detmat_1(2*inddet - 1) = costc(1)
                        detmat_1(2*inddet) = 0.d0
                    else
                        detmat_1(ipc*(inddet - 1) + 1:ipc*inddet) = costc(1:ipc)
                    end if
                end if
            end do
        end if
    end do

    amax = 0.d0
    do i = 1, nnozero_1
        if (iy_1(i) .le. norb .or. npar_eagp .ne. 0) then
            costc(1:ipc) = detmat_1(ipc*(i - 1) + 1:ipc*i)
            cost = costc(1)
            if (ipc .eq. 2 .and. abs(costc(2)) .gt. abs(costc(1))) cost = costc(2)
            if (abs(cost) .gt. abs(amax)) amax = cost
        else
            costc(1:ipc) = detmat_1(ipc*(i - 1) + 1:ipc*i)
            cost = costc(1)
            if (ipc .eq. 2 .and. abs(costc(2)) .gt. abs(costc(1))) cost = costc(2)
            if (abs(cost) .gt. abs(scaleorb(iy_1(i))))                   &
                    &        scaleorb(iy_1(i)) = cost
        end if
    end do

    if (npar_eagp .gt. 0) then

        do i = iesswr_1 + 1, iesswr_1 + iesswr_eagp
            if (abs(icd_1(i)) .ne. 1) then
                costc = 0.d0
                do k = 1, abs(icd_1(i))
                    ixt = abs(ipsipd_1(i, 2*k - 1))
                    iyt = abs(ipsipd_1(i, 2*k))
                    if (ipsipd_1(i, 2*k - 1) .lt. 0) then
                        costc(1:ipc) = costc(1:ipc) - eagp_pfaff(ipc*(ixt - 1) + 1:ipc*ixt, iyt)
                    else
                        costc(1:ipc) = costc(1:ipc) + eagp_pfaff(ipc*(ixt - 1) + 1:ipc*ixt, iyt)
                    end if
                end do

                costc(1:ipc) = costc(1:ipc)/abs(icd_1(i))
                do k = 1, abs(icd_1(i))
                    ixt = abs(ipsipd_1(i, 2*k - 1))
                    iyt = abs(ipsipd_1(i, 2*k))
                    if (ipsipd_1(i, 2*k - 1) .lt. 0) then
                        eagp_pfaff(ipc*(ixt - 1) + 1:ipc*ixt, iyt) = -costc(1:ipc)
                    else
                        eagp_pfaff(ipc*(ixt - 1) + 1:ipc*ixt, iyt) = costc(1:ipc)
                    end if
                end do
            end if
        end do

        do i = 1, nnozero_eagp
            !            do iy=1,neldiff_1
            !             do ix=iy+1,neldiff_1
            ix = ix_1(i + nnozero_1)
            iy = iy_1(i + nnozero_1)
            costc(1:ipc) = eagp_pfaff(ipc*(ix - 1) + 1:ipc*ix, iy)
            cost = costc(1)
            if (ipc .eq. 2 .and. abs(costc(2)) .gt. abs(costc(1))) cost = costc(2)
            if (abs(cost) .gt. abs(amax)) amax = cost
        end do

    end if

    write (11, *) ' #          Nonzero values of  detmat '
    do i = 1, nnozero_1
        if (iy_1(i) .le. norb .or. npar_eagp .ne. 0) then
            detmat_1(ipc*(i - 1) + 1:ipc*i) = detmat_1(ipc*(i - 1) + 1:ipc*i)/amax
            if (yes_complex .and. yes_hermite .and. iy_1(i) .eq. ix_1(i)) detmat_1(2*i) = 0.d0
            write (11, *) ix_1(i), iy_1(i), detmat_1(ipc*(i - 1) + 1:ipc*i)
        else
            detmat_1(ipc*(i - 1) + 1:ipc*i) = detmat_1(ipc*(i - 1) + 1:ipc*i)/scaleorb(iy_1(i))
            if (yes_complex) then
                write (11, '(2I9,2e18.10)') ix_1(i), iy_1(i), detmat_1(2*i - 1), detmat_1(2*i)
            else
                write (11, *) ix_1(i), iy_1(i), detmat_1(i)
            end if
        end if
    end do

    do i = 1, nnozero_eagp
        ix = ix_1(i + nnozero_1)
        iy = iy_1(i + nnozero_1)
        eagp_pfaff(ipc*(ix - 1) + 1:ipc*ix, iy) = eagp_pfaff(ipc*(ix - 1) + 1:ipc*ix, iy)/amax
        write (11, '(2I9,2e18.10)') ix, iy, eagp_pfaff(ipc*(ix - 1) + 1:ipc*ix, iy)
    end do

    if (iesdet .and. .not. mol_off .and. mol_ies) iesdet = .false.
    if (npar_eagp .gt. 0) iesdet = .false.

    do kk = 1, nnozero_1
        sortvect(kk) = (iy_1(kk) - 1)*norb + ix_1(kk)
    end do

    call dsortx(sortvect, 1, nnozero_1, indexo)

    write (11, *) '#   Grouped par.  in the chosen ordered basis'
    ! done(:) is used to avoid to repeat by chance the fixing (two values ar
    done = .true.
    do i = 1, iesswr_1
        icd_1(i) = abs(icd_1(i))
        check = .false.
        do j = 1, icd_1(i)
            ixt = abs(ipsipd_1(i, 2*j - 1))
            iyt = abs(ipsipd_1(i, 2*j))
            rtarget = (iyt - 1)*norb + ixt
            call search(rtarget, sortvect, nnozero_1, k)
            inddet = indexo(k)
            if (ipc .eq. 2) then
                if (detmat_1(2*inddet - 1) .eq. 1.d0) check = .true.
                if (detmat_1(2*inddet) .eq. 1.d0) check = .true.
            else
                if (detmat_1(inddet) .eq. 1.d0) check = .true.
            end if
        end do

        if (check .and. done(abs(ipsipd_1(i, 2)))) then

            if (.not. mol_off .and. mol_ies) then

                check = .false.
                do j = 1, 2*abs(icd_1(i))
                    if (ioptorba(abs(ipsipd_1(i, j))) .eq. 1000000) check = .true.
                end do
                if (check) then
                    write (11, *) - abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
                else
                    write (11, *) abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
                end if
            else
                write (11, *) - abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
            end if
            if (abs(ipsipd_1(i, 2)) .le. norb) then
                done = .false.
            else
                do j = 1, icd_1(i)
                    done(abs(ipsipd_1(i, 2*j))) = .false.
                end do
            end if
        else
            if (iesdet) then
                write (11, *) - abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
            elseif (.not. mol_off .and. (ioptorba(abs(ipsipd_1(i, 1))) .eq. 1000000 .or.&
                    &ioptorba(abs(ipsipd_1(i, 2))) .eq. 1000000)) then
                write (11, *) - abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
            else
                write (11, *) abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
            end if
        end if

    end do

    if (npar_eagp .ne. 0) then
        do kk = nnozero_1 + 1, nnozero_1 + nnozero_eagp
            sortvect(kk - nnozero_1) = (iy_1(kk) - 1)*neldiff_1 + ix_1(kk)
        end do
        call dsortx(sortvect, 1, nnozero_eagp, indexo)
        do i = iesswr_1 + 1, iesswr_1 + iesswr_eagp
            icd_1(i) = abs(icd_1(i))
            check = .false.
            do j = 1, icd_1(i)
                ixt = abs(ipsipd_1(i, 2*j - 1))
                iyt = abs(ipsipd_1(i, 2*j))
                rtarget = (iyt - 1)*neldiff_1 + ixt
                call search(rtarget, sortvect, nnozero_eagp, k)
                inddet = indexo(k)
                iyt = (inddet - 1)/neldiff_1 + 1
                ixt = inddet - (iyt - 1)*neldiff_1
                if (ipc .eq. 2) then
                    if (eagp_pfaff(2*ixt - 1, iyt) .eq. 1.d0) check = .true.
                    if (eagp_pfaff(2*ixt, iyt) .eq. 1.d0) check = .true.
                else
                    if (eagp_pfaff(ixt, iyt) .eq. 1.d0) check = .true.
                end if
            end do

            if (check .and. done(norb + abs(ipsipd_1(i, 2)))) then

                if (.not. mol_off .and. mol_ies) then

                    check = .false.
                    do j = 1, 2*abs(icd_1(i))
                        if (ioptorba(abs(ipsipd_1(i, j))) .eq. 1000000) check = .true.
                    end do
                    if (check) then
                        write (11, *) - abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
                    else
                        write (11, *) abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
                    end if
                else
                    write (11, *) - abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
                end if
                if (abs(ipsipd_1(i, 2)) .le. neldiff_1) then
                    done = .false.
                else
                    do j = 1, icd_1(i)
                        done(norb + abs(ipsipd_1(i, 2*j))) = .false.
                    end do
                end if
            else
                if (.not. mol_off .and. (ioptorba(abs(ipsipd_1(i, 1))) .eq. 1000000 .or.&
                        &ioptorba(abs(ipsipd_1(i, 2))) .eq. 1000000)) then
                    write (11, *) - abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
                else
                    write (11, *) abs(icd_1(i)), (ipsipd_1(i, j), j=1, 2*abs(icd_1(i)))
                end if
            end if

        end do
    end if

    deallocate (sortvect, indexo)

    if (nnozeroj_1 .gt. 0) then

        allocate (sortvect(nnozeroj_1), indexo(nnozeroj_1))

        do kk = 1, nnozeroj_1
            if (ixj_1(kk) .le. iyj_1(kk)) then
                sortvect(kk) = (iyj_1(kk) - 1)*norbj*ipj + ixj_1(kk)
            else
                sortvect(kk) = (ixj_1(kk) - 1)*norbj*ipj + iyj_1(kk)
            end if
        end do

        call dsortx(sortvect, 1, nnozeroj_1, indexo)

    end if

    !
    do i = 1, iesfreer_1
        if (abs(icj_1(i)) .ne. 1) then
            cost = 0.d0
            do k = 1, abs(icj_1(i))
                ixt = abs(ipsipj_1(i, 2*k - 1))
                iyt = abs(ipsipj_1(i, 2*k))
                if (ixt .le. iyt) then
                    rtarget = (iyt - 1)*norbj*ipj + ixt
                else
                    rtarget = (ixt - 1)*norbj*ipj + iyt
                end if
                call search(rtarget, sortvect, nnozeroj_1, kk)
                inddet = indexo(kk)

                if (ipsipj_1(i, 2*k) .lt. 0 .or. ipsipj_1(i, 2*k - 1) .lt. 0) then
                    cost = cost - jasmat_1(inddet)
                else
                    cost = cost + jasmat_1(inddet)
                end if
            end do

            cost = cost/abs(icj_1(i))
            do k = 1, abs(icj_1(i))
                ixt = abs(ipsipj_1(i, 2*k - 1))
                iyt = abs(ipsipj_1(i, 2*k))
                if (ixt .le. iyt) then
                    rtarget = (iyt - 1)*norbj*ipj + ixt
                else
                    rtarget = (ixt - 1)*norbj*ipj + iyt
                end if
                call search(rtarget, sortvect, nnozeroj_1, kk)

                inddet = indexo(kk)

                if (ipsipj_1(i, 2*k) .lt. 0 .or. ipsipj_1(i, 2*k - 1) .lt. 0) then
                    jasmat_1(inddet) = -cost
                else
                    jasmat_1(inddet) = cost
                end if
            end do
        end if

    end do

    write (11, *) ' #          Nonzero values of  jasmat '
    do i = 1, nnozeroj_1
        if (ioptorbja_1(ixj_1(i)) .eq. 200 .and. ioptorbja_1(iyj_1(i)) .eq. 200) then
            write (11, *) ixj_1(i), iyj_1(i), 0.d0
        else
            write (11, *) ixj_1(i), iyj_1(i), jasmat_1(i)
        end if
    end do

    if (iessz) then
        do i = 1, iesfreer_1
            if (abs(icj_1(i)) .ne. 1) then
                cost = 0.d0
                do k = 1, abs(icj_1(i))
                    ixt = abs(ipsipj_1(i, 2*k - 1))
                    iyt = abs(ipsipj_1(i, 2*k))
                    if (ixt .le. iyt) then
                        rtarget = (iyt - 1)*norbj + ixt
                    else
                        rtarget = (ixt - 1)*norbj + iyt
                    end if
                    call search(rtarget, sortvect, nnozeroj_1, kk)
                    inddet = indexo(kk)

                    if (ipsipj_1(i, 2*k) .lt. 0 .or. ipsipj_1(i, 2*k - 1) .lt. 0) then
                        cost = cost - jasmatsz_1(inddet)
                    else
                        cost = cost + jasmatsz_1(inddet)
                    end if
                end do

                cost = cost/abs(icj_1(i))
                do k = 1, abs(icj_1(i))
                    ixt = abs(ipsipj_1(i, 2*k - 1))
                    iyt = abs(ipsipj_1(i, 2*k))
                    if (ixt .le. iyt) then
                        rtarget = (iyt - 1)*norbj + ixt
                    else
                        rtarget = (ixt - 1)*norbj + iyt
                    end if
                    call search(rtarget, sortvect, nnozeroj_1, kk)
                    inddet = indexo(kk)

                    if (ipsipj_1(i, 2*k) .lt. 0 .or. ipsipj_1(i, 2*k - 1) .lt. 0) then
                        jasmatsz_1(inddet) = -cost
                    else
                        jasmatsz_1(inddet) = cost
                    end if
                end do
            end if
        end do
        write (11, *) ' #          Nonzero values of SZ jasmat  '
        do i = 1, nnozeroj_1
            if (ioptorbja_1(ixj_1(i)) .eq. 200 .and. ioptorbja_1(iyj_1(i)) .eq. 200) then
                write (11, *) ixj_1(i), iyj_1(i), 0.d0
            else
                write (11, *) ixj_1(i), iyj_1(i), jasmatsz_1(i)
            end if
        end do
    end if

    write (11, *) '# Eq. par. in the 3-body Jastrow                &
            & in the chosen basis '
    do i = 1, iesfreer_1
        if (ioptorbja_1(abs(ipsipj_1(i, 1))) .eq. 200                  &
                &.and. ioptorbja_1(abs(ipsipj_1(i, 2))) .eq. 200) then
            write (11, *) - abs(icj_1(i)), (ipsipj_1(i, j), j=1, 2*abs(icj_1(i)))
        elseif (nel - nelup .eq. 0 .and. ipj .eq. 2 &
                .and. (abs(ipsipj_1(i, 1)) .gt. norbj .or. abs(ipsipj_1(i, 2)) .gt. norbj)) then
            write (11, *) - abs(icj_1(i)), (ipsipj_1(i, j), j=1, 2*abs(icj_1(i)))
            !         elseif(nel-nelup.ge.1.and.ipj.eq.2.and.&
            ! &((ioptorbja_1(abs(ipsipj_1(i,1))).eq.200.and.abs(ipsipj_1(i,1)).gt.norbj)&
            ! &.or.&
            ! &((ioptorbja_1(abs(ipsipj_1(i,2))).eq.200.and.abs(ipsipj_1(i,2)).gt.norbj)&
            ! &)))  then
            !! &.or.(abs(ipsipj_1(i,1)).gt.norbj.and.abs(ipsipj_1(i,2)).le.norbj))) then
            !         write(11,*) -abs(icj_1(i)),(ipsipj_1(i,j),j=1,2*abs(icj_1(i)))
        elseif (.not. molj_off .and. (ioptorbja_1(abs(ipsipj_1(i, 1))) .eq. 1000000 .or.&
                &ioptorbja_1(abs(ipsipj_1(i, 2))) .eq. 1000000)) then
            write (11, *) - abs(icj_1(i)), (ipsipj_1(i, j), j=1, 2*abs(icj_1(i)))
        else
            write (11, *) abs(icj_1(i)), (ipsipj_1(i, j), j=1, 2*abs(icj_1(i)))
        end if
    end do

    if (nnozeroj_1 .gt. 0) deallocate (sortvect, indexo)

    !           ALLOCATE(cek_write(iesupind_1))
    !           cek_write=.true.

    !=========================================================
    write (11, *) '# Eq. par.in the atomic Det par. in the chosen basis '
    do i = 1, iesupind_1
        !    ZZZ   lines modified by Zen
        if (flagconstrains .eq. 1 .or. flagconstrains .eq. 3) then
            write (11, *) abs(icp_1(i)), (ipsipp_1(i, j), j=1, abs(icp_1(i)))
        else
            write (11, *) icp_1(i), (ipsipp_1(i, j), j=1, abs(icp_1(i)))
        end if
        !    ZZZ   end lines modified by Zen
    end do

    !           write(*,*)'deallocate cek_write'
    !==================================================

    !        write(*,*)'allocate cek_writej'

    !           ALLOCATE(cek_writej(iesmind_1))
    !           cek_writej=.true.

    write (11, *) &
        '# Eq. par. in the atomic 3-body  par. in the chosen basis '

    do i = 1, iesmind_1

        !              do ii=1,abs(icq_1(i))
        !               if(cek_zeroj(ipsipq_1(i,ii) ))  cek_writej(i)=.false.
        !              enddo

        !              if (cek_writej(i)) then
        !            this is not necessary
        !              iesmind_new=iesmind_new+1

        !    ZZZ   lines modified by Zen
        if (flagconstrains .eq. 2 .or. flagconstrains .eq. 3) then
            write (11, *) abs(icq_1(i)), (ipsipq_1(i, j), j=1, abs(icq_1(i)))
        else
            write (11, *) icq_1(i), (ipsipq_1(i, j), j=1, abs(icq_1(i)))
        end if
        !    ZZZ   end lines modified added by Zen

        !              endif

    end do

    !           DEALLOCATE(cek_writej)!,cek_zero)

    !
    !====================================================
    stop

end
subroutine independent_complex(n, ncoeff, psi_in, psi_out, ipiv, yesfixphase)
    implicit none
    integer n, i, ipiv(n), jmax, j, ncoeff
    complex*16 psi_in(ncoeff, n), psi_out(ncoeff, n), cost
    real*8 maxv, costr
    logical yesfixphase
    ipiv = 0
    do i = 1, n
        !       first step set to zero all coefficient

        psi_out(:, i) = psi_in(:, i)

        if (.not. yesfixphase) then
            do j = 1, i - 1
                if (ipiv(j) .ne. 0) then
                    cost = psi_out(ipiv(j), i)
                    psi_out(:, i) = psi_out(:, i) - cost*psi_out(:, j)
                end if
            end do
        end if
        !         Now find the maximum

        maxv = 0.d0
        jmax = 0
        do j = 1, ncoeff
            costr = abs(psi_out(j, i))
            if (costr .gt. maxv) then
                jmax = j
                maxv = costr
            end if
        end do

        ipiv(i) = jmax

        if (jmax .gt. 0) then

            cost = psi_out(jmax, i)

            psi_out(:, i) = psi_out(:, i)/cost
            ! to avoid roundoff
            psi_out(jmax, i) = (1.d0, 0.d0)

            if (.not. yesfixphase) then
                do j = 1, i - 1
                    cost = psi_out(jmax, j)
                    psi_out(:, j) = psi_out(:, j) - cost*psi_out(:, i)
                    !         cost=psi_out(ipiv(j),j)
                    !         psi_out(:,j)=psi_out(:,j)/cost
                end do
            end if
        else
            !         It may happen that the vector is vanishing
            psi_out(:, i) = (0.d0, 0.d0)
        end if

    end do
    return
end
!
subroutine independent(n, ncoeff, psi_in, psi_out, ipiv)
    implicit none
    integer n, i, ipiv(n), jmax, j, ncoeff
    real*8 psi_in(ncoeff, n), psi_out(ncoeff, n), cost, maxv
    ipiv = 0
    do i = 1, n
        !       first step set to zero all coefficient

        psi_out(:, i) = psi_in(:, i)

        do j = 1, i - 1
            if (ipiv(j) .ne. 0) then
                cost = psi_out(ipiv(j), i)
                psi_out(:, i) = psi_out(:, i) - cost*psi_out(:, j)
            end if
        end do
        !         Now find the maximum

        maxv = 0.d0
        jmax = 0
        do j = 1, ncoeff
            cost = abs(psi_out(j, i))
            if (cost .gt. maxv) then
                jmax = j
                maxv = cost
            end if
        end do

        ipiv(i) = jmax

        if (jmax .gt. 0) then

            cost = psi_out(jmax, i)

            psi_out(:, i) = psi_out(:, i)/cost
            ! to avoid roundoff
            psi_out(jmax, i) = 1.d0

            do j = 1, i - 1
                cost = psi_out(jmax, j)
                psi_out(:, j) = psi_out(:, j) - cost*psi_out(:, i)
                !         cost=psi_out(ipiv(j),j)
                !         psi_out(:,j)=psi_out(:,j)/cost
            end do
        else
            !         It may happen that the vector is vanishing
            psi_out(:, i) = 0.d0
        end if

    end do
    return
end
!================
subroutine makeumat(ipc, n, ncoeff, psi_in, psi_out, umat)
    use constants, only: zone, zzero
    implicit none
    integer n, i, j, ncoeff, info, lwork, ipc
    real*8 psi_in(ncoeff*ipc, n), psi_out(ncoeff*ipc, n), umat(ipc*n, n)
    real*8, dimension(:, :), allocatable :: smat, psip
    integer, dimension(:), allocatable :: ipsip

    if (ncoeff .le. 0 .or. n .le. 0) return

    allocate (smat(ipc*n, n), psip(ipc*n, n), ipsip(n))

    !      The matrix psip has to have n^2 elemnts
    !      the psi_out due to roundoff are likely to be independent
    !      even when psi_in are dependent

    if (ipc .eq. 1) then
        call dgemm('T', 'N', n, n, ncoeff, 1.d0, psi_out, ncoeff, psi_out, ncoeff&
                &, 0.d0, smat, n)
        call dgetrf(n, n, smat, n, ipsip, info)
    else
        call zgemm('C', 'N', n, n, ncoeff, zone, psi_out, ncoeff, psi_out, ncoeff&
                &, zzero, smat, n)
        call zgetrf(n, n, smat, n, ipsip, info)
    end if

    !       write(6,*) ' Output dgetrf '

    !       do i=1,n
    !       write(6,*) i,smat(i,i)
    !       enddo

    if (info .eq. 0) then
        lwork = n*n
        if (ipc .eq. 1) then

            call dgetri(n, smat, n, ipsip, psip, lwork, info)

            call dgemm('T', 'N', n, n, ncoeff, 1.d0, psi_out, ncoeff, psi_in, ncoeff &
                    &, 0.d0, psip, n)

            call dgemm('T', 'N', n, n, n, 1.d0, psip, n, smat, n                     &
                    &, 0.d0, umat, n)

        else
            call zgetri(n, smat, n, ipsip, psip, lwork, info)

            call zgemm('C', 'N', n, n, ncoeff, zone, psi_out, ncoeff, psi_in, ncoeff &
                    &, zzero, psip, n)

            call zgemm('C', 'N', n, n, n, zone, psip, n, smat, n                     &
                    &, zzero, umat, n)
            call conjmat(n, n, umat, n)
        end if

    else
        write (6, *) ' Orbitals are not independent, ERROR in dgetrf !!!! ', info
        do i = 1, n
            write (6, *) i, smat(ipc*(i - 1) + 1, ipc*(i - 1) + 1)
        end do
        !       stop
    end if

    deallocate (smat, psip, ipsip)
    return
end
!...................
subroutine lambda_transf(ipc, n, m, lambda, umat)

    use constants, only: zone, zzero
    implicit none

    integer n, m, neldiff, i, j, ipc
    real*8 lambda(ipc*n, m), umat(ipc*n, n)
    real*8, dimension(:, :), allocatable :: psip
    allocate (psip(ipc*n, m))

    psip = 0.d0 ! to avoid NaN

    !          New  lambda = emme lambda emme^{T}

    if (m .gt. n) then
        neldiff = m - n
        call dcopy(ipc*neldiff*n, lambda(1, n + 1), 1, psip(1, n + 1), 1)
    end if

    if (ipc .eq. 1) then

        call dgemm('N', 'N', n, n, n, 1.d0, lambda, n, &
                &    umat, n, 0.d0, psip, n)

        call dgemm('T', 'N', n, m, n, 1.d0, umat, n, &
                &    psip, n, 0.d0, lambda, n)

    else

        call zgemm('N', 'N', n, n, n, zone, lambda, n, &
                &    umat, n, zzero, psip, n)

        call zgemm('T', 'N', n, m, n, zone, umat, n, &
                &     psip, n, zzero, lambda, n)
    end if

    deallocate (psip)
    return
end

subroutine search(target, sortvect, n, index)
    implicit none
    real*8 target, sortvect(*)
    integer n, index, i, k
    i = 1
    index = n
    do
        k = (i + index)/2
        if (target .lt. sortvect(k)) then
            index = k
        else
            i = k
        end if
        if (i + 1 .ge. index) then
            if (sortvect(index) .ne. target) index = i
            return
        end if
    end do
    return
end

