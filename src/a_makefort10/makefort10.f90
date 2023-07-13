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

program makefort10
    use symm_data
    use constants
    use mod_orbital, only: parsymm, print_orbitals &
            &, apply_symm_to_orbitals, apply_symm_to_forces&
            &, orbmap, generate_orbidx, par_symm, read_orbitals &
            &, orbital, lsym_type, read_atoms, atomstypes, yesmolat &
            &, yesmolatj, yesalloc_jas, atomic_jasmat, real_contracted

    implicit none

    real(8) :: at(3, 3), bg(3, 3), alphap, betap, gammap, unit_volume, omega, volume_new
    integer :: i1, i2, i3, i4, j1, ix, iy, is, ione, itwo, ionpos, i2map, ncellmap, signmap, nrec&
            &, njastot, ndettot, nrecj, nrecj_ok, nrecf, nrecagp, nrecskw, ind, indmap, imap&
            &, i, j, k, jj, nsymu, ntrau, nrecj_save, index1, index2, niesd, const_term, npar_eagp, tf, case_map

    !character(30) :: gformat
    character(100) :: name_tool
    character(20) :: str

    integer :: ntra, ntraind, vecpbc(3)
    integer, allocatable :: detorbidx(:), detorbidxz(:) ! index that map orbitals of equal atoms
    integer, allocatable :: jasorbidx(:), jasorbidxz(:), orb2atom(:)
    integer, allocatable :: atoms_maprot(:, :), atoms_maptra(:, :)

    type(parsymm), allocatable :: eqdet(:), eqjas(:)
    double precision, allocatable :: detmat(:, :), jasmat(:, :)
    logical, allocatable :: jasyes(:, :), jasyes_1(:, :), jasyes_2(:, :), detyes(:, :), forceyes(:)
    type(orbmap) :: det_map, jas_map
    real*8, dimension(:, :), allocatable :: trasl

    integer :: natoms, ntotatoms, nel_read, nbas
    integer :: ntyp ! Number of atomic types
    integer, allocatable :: mytype(:) ! the type of each atom
    integer :: nxyz(3), ncell, axyz(3, 3)
    real(8) :: celldm(6), smallcell(3), newcell(3, 3), cellscale(3), rs_read, rs_try, L_read
    character(len=10) :: posunits
    character(len=12) :: unit_crystal
    real(8), allocatable :: rion(:, :)
    real(8), allocatable :: zeta(:, :), tmp_zeta(:, :)
    logical, allocatable :: J3_off(:)
    type(orbital), pointer :: detorb(:), jasorb(:), detorbat(:), jasorbat(:)
    integer, allocatable :: symrot(:, :), symtra(:, :), ipsip(:)&
  &, recordsymj(:, :, :), lenrecj(:), recordsym(:, :, :), lenrec(:), indion(:)&
  &, recordsymf(:, :, :), lenrecf(:), orbf_map(:, :), recordsymj_1(:, :, :)&
  &, recordsymj_2(:, :, :), lenrecj_1(:), lenrecj_2(:), lenrecagp(:)&
  &, lenrecskw(:), recordsymskw(:, :, :), recordsymagp(:, :, :), indionz(:)
    logical :: found, is_ok, read_vecpbc, yes_hermite, accept
    integer :: nel, nelup, neldiff, numpaired, twobody, ndetpar, njaspar
    ! double precision :: twobodypar, onebodypar
    integer :: n_onebody, n_twobody !KN
    double precision, allocatable :: twobodypar(:) !KN
    double precision, allocatable :: onebodypar(:) !KN
    double precision :: phase(3), phasedo(3)
    logical :: apbc(3), apbcmap(3), pbcfort10, noonebody, complexfort10, complexfort10_sav, yes_trivial, no_4body_jas, nopseudo
    character(20) :: filling
    type(atomstypes), allocatable :: atypes(:) ! Atoms types
    character(20) :: orbtype, jorbtype
    integer :: nshelldet, nshelljas, counter
    integer :: norb_unpaired, nlambda_unp
    integer, allocatable :: orb_unpaired(:, :), orb_numbers(:), cut_hybrid(:), cut_hybridj(:)
    logical :: unpaired, nosym_contrj

    logical :: no_orthorombic
    integer :: max_factor(3)
    real(8), allocatable :: new_rion(:, :)
    real(8) :: test_rion(3), test_rion_before(3)
    integer :: max_ions
    ! integer, external:: makenpipo
    integer, external :: makenpip
    logical :: readatoms ! read atomic wave-function from file
    logical :: rot_det ! if false exclude rotations from the determinant
    logical :: rot_jas ! if false exclude rotations from the determinant
    logical :: rot_pfaff ! if false exclude rotations from the triplet component of the pfaffian (usually odd)
    !logical :: makes2s1p ! make symmetries as in makes2s1p, not implemented yet!
    logical :: symmagp, yesbump, yes_crystal, yes_crystalj, genjason, opposite_phase&
            &, yes_tilted

    double precision, pointer :: atomic_detmat(:, :)
    integer :: njasorb, ndetorb, ndetorbat, njasorbat
    integer :: neqdet, neqjas, shiftbeta

    logical :: onlycontrdet ! put minus sign to zeta in the determinant
    ! for the optimization of the contraction coefficients
    ! only

    logical :: yes_pfaff, nodownpfaff, nodownpfaff_read, nouppfaff ! yes_pfaff creates a pfaffian wf
    !that can be complete or without the up up or the down down term depending on the
    !other two variables

    logical :: onlycontrjas ! put minus sign to zeta in the jastrow
    ! for the optimization of the contraction coefficients
    ! only

    logical :: readunpaired, forcesymm
    real*8 scale_jasfat

    ! output version information
    call print_version

    namelist /system/ natoms, posunits, nxyz, celldm, at, phase, phasedo&
            &, pbcfort10, complexfort10, real_contracted&
            &, ntyp, rs_read, write_log, axyz, nel_read, L_read, yes_pfaff&
            &, nodownpfaff, nouppfaff, yes_tilted, unit_crystal

#ifdef __PORT
    namelist /electrons/ twobody, twobodypar, filling, noonebody, readatoms&
            &, orbtype, nel, neldiff, numpaired, jorbtype, onlycontrdet&
            &, onlycontrjas, shiftbeta, readunpaired, vecpbc, yesbump&
            &, niesd, yes_crystal, yes_crystalj, no_4body_jas&
            &, nopseudo, scale_jasfat
    write (*, *)
    write (*, *) ' ! ! ! WARNING DUE TO POSSIBLE BUG IN NVIDIA COMPILER IS NOT POSSIBLE TO SET ONEBODY PAR ! ! ! '
    write (*, *)
    write (*, *)
#else
    namelist /electrons/ twobody, twobodypar, filling, noonebody, readatoms&
            &, orbtype, nel, neldiff, numpaired, jorbtype, onlycontrdet&
            &, onlycontrjas, shiftbeta, readunpaired, vecpbc, yesbump&
            &, onebodypar, niesd, yes_crystal, yes_crystalj, no_4body_jas&
            &, nopseudo, scale_jasfat
#endif

    namelist /symmetries/ nosym, notra, forces_sym, notra_forces, nosym_forces&
            &, eqatoms, eq_intatoms, rot_det, rot_jas, rot_pfaff, nosym_contr&
            &, nosym_contrj, symmagp, forcesymm

    !   AAA    Lines to be added just after all definitions of variables.

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'makefort10'
        call help_online(name_tool)
        stop
    end if

    write (*, *)
    write (*, *) ' * * * Creates fort.10 for solids and open systems * * * '
    write (*, *)
    write (*, *)

    ! reading input file and generating orthorombic
    ! cell from non-orthorombic one
    call read_input(5)

    orbtype = trim(orbtype)
    jorbtype = trim(jorbtype)
    tf = twobody

    if (abs(twobody) .ge. 100) then
        case_map = abs(twobody)/100
        if (twobody .gt. 0) then
            tf = tf - 100*case_map
        else
            tf = tf + 100*case_map
        end if
    end if

    if (tf .eq. -12 .or. tf .eq. -22 .or. tf .eq. -26 .or. tf .eq. -27&
            &.or. tf .eq. -30 .or. tf .eq. -31) then
        genjason = .true.
        ipj = 2
    else
        genjason = .false.
        ipj = 1
    end if

    ! set n_twobody added by Kosuke Nakano on 18.Feb.
    if (tf .eq. -26 .or. tf .eq. -27 .or. tf .eq. -30 .or. tf .eq. -31) then
        n_twobody = 2
    else
        n_twobody = 1
        if (twobodypar(2) .ne. 1.0d0) then
            write (6, *) ' Warning twobodypar(2) is ignored.'
        end if
    end if
    write (6, *) ' n_twobody is set', n_twobody

    if (complexfort10_sav .neqv. complexfort10) then
        if (complexfort10) write (6, *) ' Warning the wavefunction is considered complex'
        if (.not. complexfort10) write (6, *) ' Warning the wavefunction is considered real'
    elseif (complexfort10) then
        write (6, *) ' Creating complex wave function...'
    end if

    if (nel == -1 .and. rs_read == -1 .and. L_read == -1.d0) then
        call errore("makefort10", " Error you have to set nel , rs_read or L_read!! ", 1)
    end if

    if (filling == "undefined") filling = "diagonal" ! default filling for jastrow and determinant matrices

    if (nshelljas .eq. -1) then
        if (readatoms) then
            nshelljas = 1
        else
            nshelljas = 0
        end if
    end if

    if (nshelljas == 0) noonebody = .true.

    ! print characteristics of all orbitals in the log file
    call print_orbitals(detorb, ndetorb, "DETERMINANT")
    if (nshelljas .ne. 0) call print_orbitals(jasorb, njasorb, "JASTROW")

    ! find number of symmetries of the Bravais lattice (only orthorombic)
    call find_rotations
    ! ------------------

    allocate (detorbidx(ndetorb), detorbidxz(ndetorb))
    allocate (atoms_maprot(ntotatoms, nsym), atoms_maptra(ntotatoms, ntraind))
    !write(6,*) ' symmetry map '
    !do is=1,nsym
    !write(6,*) ' Symmetry det =',is
    !  do i1=1,ndetorb
    !  write(6,*) i1,det_map%iorb(i1,is),det_map%comp(i1,is)
    !  enddo
    !enddo

    nrot = 48

    at_qe(:, 1) = newcell(:, 1)/dsqrt(sum(newcell(:, 1)**2))
    at_qe(:, 2) = newcell(:, 2)/dsqrt(sum(newcell(:, 2)**2))
    at_qe(:, 3) = newcell(:, 3)/dsqrt(sum(newcell(:, 3)**2))

    call init_cell

    call update_atoms_map(apbc)

    if (vecpbc(1) .eq. -1) then
        read_vecpbc = .false.

        do is = 1, 3
            if (apbc(is)) then
                vecpbc(is) = 1
            else
                vecpbc(is) = 0
            end if
        end do

    else
        read_vecpbc = .true.
    end if

    ! check Bravais lattice symmetries allowed by the atomic structure
    call check_symm(vecpbc, rot_det, yes_hermite)
    ! ------------------------------

    call update_atoms_map(apbc)

    call apply_symm_to_orbitals(detorb, ndetorb, det_map)

    ! --------- CHECK ------------
    ! write(6,*) ' symmetry rotation '
    ! do is=1,nsym
    !   write(6,*) ' Symmetry det =',is
    !   do i1=1,ntotatoms
    !   write(6,*) i1,atoms_maprot(i1,is)
    !   enddSo
    ! enddo

    ndettot = ndetorb*ncell

    write (6, *) 'ndetorb', ndetorb, ndettot, ncell

    allocate (orb2atom(ndettot))

    if (.not. yesmolat) then
        allocate (symrot(ndettot, nsym), symtra(ndettot, ntraind))
        call update_symrot(nsym, symrot, ndettot, atoms_maprot, detorb, ndetorb, det_map)
        write (6, *) ' Orb to atom mapping '
        do i1 = 1, ndetorb
            do i2 = 1, ncell
                ionpos = detorb(i1)%kion + (i2 - 1)*natoms
                i3 = i1 + (i2 - 1)*ndetorb
                orb2atom(i3) = ionpos
            end do
        end do
        do i1 = 1, ndettot
            write (6, *) i1, orb2atom(i1)
        end do
    end if

    !write(6,*) ' Symmetry rot '
    !do is=1,nsym
    ! write(6,*) ' Sym =',is
    ! do i1=1,ndettot
    ! write(6,*) i1,symrot(i1,is)
    ! enddo
    !enddo

    if (nosym) then
        nsym = 1 ! Use only the identity
        write (6, *) ' Warning forced not to use all point symmetries nsym=', nsym
    end if

    allocate (orbf_map(3, nsym))

    call apply_symm_to_forces(orbf_map)

    write (6, *) ' Mapping force components '

    do i = 1, nsym
        write (6, *) i, orbf_map(1:3, i)
    end do

    if (.not. yesmolat) then
        allocate (detyes(ndettot, ndettot)) ! logical matrix which keeps
        ! track of detmat occupation
        ! allocation of geminal matrix
        ! double allocation in the case of complex wave function
        if (.not. complexfort10) then
            allocate (detmat(ndetorb*ncell, ndetorb*ncell))
        else
            allocate (detmat(2*ndetorb*ncell, ndetorb*ncell))
        end if

        if (complexfort10 .and. .not. yes_hermite) then
            ntra = 1
        else
            ntra = ntraind
        end if
        call update_symtra(ntra, symtra, ndettot, atoms_maptra, detorb, ndetorb)
    end if
    call generate_indion
    write (6, *) ' Independent atoms in the orthorombic cell '
    do i1 = 1, natoms
        write (6, *) i1, indion(i1), indionz(i1)
    end do

    call generate_orbidx(zeta, natoms, detorbidx, detorb, ndetorb, indion)
    if (.not. eqatoms) then
        detorbidxz = detorbidx
    else
        call generate_orbidx(zeta, natoms, detorbidxz, detorb, ndetorb, indionz)
    end if
    if (complexfort10 .and. .not. yes_hermite) then
        ntra = 1 ! only the identity
    end if

    if (.not. yesmolat) then

        nrec = 1
        allocate (ipsip(2*nsym*ntra))
        allocate (recordsym(2, nsym*ntra, nrec), lenrec(nrec))
        nrec = 0

        !  first run count the number of records

        write (6, *) ' before makelambda det  nsym ntra=', nsym, ntra
        call makelambda(ndettot, symrot, nsym, symtra, ntra, detyes, recordsym, lenrec, nrec, &
                        symmagp, ipsip, rion, ntotatoms, orb2atom, cellscale, deps, yes_hermite, &
                        opposite_phase, zeta, J3_off, .false., .false., cut_hybrid)
        write (6, *) ' Output records det =', nrec
        deallocate (recordsym, lenrec)
        allocate (recordsym(2, nsym*ntra, nrec), lenrec(nrec))
        call makelambda(ndettot, symrot, nsym, symtra, ntra, detyes, recordsym, lenrec, nrec, &
                        symmagp, ipsip, rion, ntotatoms, orb2atom, cellscale, deps, yes_hermite, &
                        opposite_phase, zeta, J3_off, .false., .false., cut_hybrid)

        !If ipf2 It uses the lambda for the up-down part of the matrix and calculates the
        !skew symmetric quater to fill the pfaffian matrix
        if (ipf .eq. 2) then
            nbas = ndettot
            nrecagp = nrec
            allocate (recordsymagp(2, ntra*nsym, nrecagp), lenrecagp(nrecagp))
            recordsymagp(:, :, :) = recordsym(:, :, :)
            lenrecagp(:) = lenrec(:)

            nrecskw = 1

            if (.not. rot_pfaff) then
                do i1 = 1, nsym
                    do i2 = 1, nbas
                        symrot(i2, i1) = i2
                    end do
                end do
            end if

            allocate (recordsymskw(2, nsym*ntra, nrecskw), lenrecskw(nrecskw))
            call makeskew(nbas, symrot, nsym, symtra, ntra, detyes, recordsymskw, &
                          lenrecskw, nrecskw, ipsip, rion, ntotatoms, orb2atom, cellscale, &
                          deps, zeta, yes_hermite)
            deallocate (recordsymskw, lenrecskw)
            allocate (recordsymskw(2, ntra*nsym, nrecskw), lenrecskw(nrecskw))
            call makeskew(nbas, symrot, nsym, symtra, ntra, detyes, recordsymskw, &
                          lenrecskw, nrecskw, ipsip, rion, ntotatoms, orb2atom, cellscale, &
                          deps, zeta, yes_hermite)

            if ((nouppfaff .and. .not. nodownpfaff) .or. (nodownpfaff .and. .not. nouppfaff)) then
                nrec = nrecagp + nrecskw
            elseif (.not. nouppfaff .and. .not. nodownpfaff) then
                nrec = nrecagp + 2*nrecskw
            else
                nrec = nrecagp
            end if

            deallocate (lenrec, recordsym, detyes)
            allocate (lenrec(nrec), recordsym(2, ntra*nsym, nrec), detyes(2*nbas, 2*nbas))
            call makepfaff(nbas, nsym, ntra, detyes, recordsymagp, lenrecagp, recordsymskw, &
                           lenrecskw, recordsym, lenrec, nrec, nrecagp, nrecskw, nouppfaff, nodownpfaff, &
                           opposite_phase)
            deallocate (lenrecskw, lenrecagp, recordsymagp, recordsymskw)
            deallocate (detmat)
            nbas = nbas*ipf
            ndettot = nbas
            if (.not. complexfort10) then
                allocate (detmat(nbas, nbas))
            else

                allocate (detmat(2*nbas, nbas))
            end if
        end if
    else

        ! first count the basis used.
        nbas = 0
        do i2 = 1, ndetorb
            if (detorb(i2)%ioptorb .eq. 900000) nbas = nbas + 1
        end do
        ndetorbat = nbas

        allocate (detorbat(ndetorbat))
        nbas = 0
        do i2 = 1, ndetorb
            if (detorb(i2)%ioptorb .eq. 900000) then
                nbas = nbas + 1
                detorbat(nbas) = detorb(i2)
            end if
        end do

        nbas = nbas*ncell
        write (6, *) ' n basis found =', nbas

        do i2 = 1, ndetorbat
            write (6, *) i2, detorbat(i2)%ioptorb, detorbat(i2)%kion
        end do
        !    No rotation symmetry imposed.
        nsym = 1
        if (complexfort10 .and. .not. yes_hermite) then
            ntra = 1
        else
            ntra = ntraind
        end if
        allocate (symrot(nbas, nsym), symtra(nbas, ntra))
        !    Identity for symtra
        do i2 = 1, nbas
            symrot(i2, 1) = i2
        end do

        call update_symtra(ntra, symtra, nbas, atoms_maptra, detorbat, ndetorbat)
        write (6, *) ' Symmetry tra ', ntra
        do is = 1, ntra
            write (6, *) ' Sym =', is
            do i1 = 1, nbas
                write (6, *) i1, symtra(i1, is)
            end do
        end do
        allocate (detyes(nbas, nbas))
        if (.not. complexfort10) then
            allocate (detmat(nbas, nbas))
        else
            allocate (detmat(2*nbas, nbas))
        end if
        ndettot = nbas

        write (6, *) ' Orb to atom mapping '
        do i1 = 1, ndetorbat
            do i2 = 1, ncell
                ionpos = detorbat(i1)%kion + (i2 - 1)*natoms
                i3 = i1 + (i2 - 1)*ndetorbat
                orb2atom(i3) = ionpos
            end do
        end do
        do i1 = 1, ndettot
            write (6, *) i1, orb2atom(i1)
        end do

        nrec = 1
        allocate (ipsip(2*nsym*ntra))
        allocate (recordsym(2, nsym*ntra, nrec), lenrec(nrec))
        nrec = 0

        !    first run count the number of records

        !     write(6,*) ' before makelambda =',nsym
        call makelambda(nbas, symrot, nsym, symtra, ntra, detyes, recordsym, lenrec, nrec, &
                        symmagp, ipsip, rion, ntotatoms, orb2atom, cellscale, deps, yes_hermite, &
                        opposite_phase, zeta, J3_off, .false., .false., cut_hybrid)

        write (6, *) ' Output records det =', nrec

        deallocate (recordsym, lenrec)

        allocate (recordsym(2, nsym*ntra, nrec), lenrec(nrec))

        call makelambda(nbas, symrot, nsym, symtra, ntra, detyes, recordsym, lenrec, nrec, &
                        symmagp, ipsip, rion, ntotatoms, orb2atom, cellscale, deps, yes_hermite, &
                        opposite_phase, zeta, J3_off, .false., .false., cut_hybrid)

        !If ipf2 It uses the lambda for the up-down part of the matrix and calculates the
        !skew symmetric quater to fill the pfaffian matrix
        if (ipf .eq. 2) then
            nbas = ndettot
            nrecagp = nrec
            allocate (recordsymagp(2, ntra*nsym, nrecagp), lenrecagp(nrecagp))
            recordsymagp(:, :, :) = recordsym(:, :, :)
            lenrecagp(:) = lenrec(:)
            if (.not. rot_pfaff) then
                do i1 = 1, nsym
                    do i2 = 1, nbas
                        symrot(i2, i1) = i2
                    end do
                end do
            end if

            !       write(6,*) ' Output record agp '
            !       do i1=1,nrecagp
            !       write(6,*) i1,(recordsymagp(1,i2,i1),recordsymagp(2,i2,i1),i2=1,lenrecagp(i1))
            !       enddo
            nrecskw = 1
            allocate (recordsymskw(2, nsym*ntra, nrecskw), lenrecskw(nrecskw))
            call makeskew(nbas, symrot, nsym, symtra, ntra, detyes, recordsymskw, &
                          lenrecskw, nrecskw, ipsip, rion, ntotatoms, orb2atom, cellscale, &
                          deps, zeta, yes_hermite)
            deallocate (recordsymskw, lenrecskw)
            allocate (recordsymskw(2, ntra*nsym, nrecskw), lenrecskw(nrecskw))
            call makeskew(nbas, symrot, nsym, symtra, ntra, detyes, recordsymskw, &
                          lenrecskw, nrecskw, ipsip, rion, ntotatoms, orb2atom, cellscale, &
                          deps, zeta, yes_hermite)

            !       write(6,*) ' Output record skew '
            !       do i1=1,nrecskw
            !       write(6,*) i1,(recordsymskw(1,i2,i1),recordsymskw(2,i2,i1),i2=1,lenrecskw(i1))
            !       enddo

            if (nouppfaff .or. nodownpfaff .or. symmagp) then
                nrec = nrecagp + nrecskw
            else
                nrec = nrecagp + 2*nrecskw
            end if

            deallocate (lenrec, recordsym, detyes)
            allocate (lenrec(nrec), recordsym(2, ntra*nsym, nrec), detyes(2*nbas, 2*nbas))

            call makepfaff(nbas, nsym, ntra, detyes, recordsymagp, lenrecagp, recordsymskw, &
                           lenrecskw, recordsym, lenrec, nrec, nrecagp, nrecskw, nouppfaff, nodownpfaff, &
                           opposite_phase)
            !       write(6,*) ' Output record tot ',nouppfaff,nodownpfaff,nrec,nrecskw,nrecagp
            !       do i1=1,nrec
            !       write(6,*) i1,(recordsym(1,i2,i1),recordsym(2,i2,i1),i2=1,lenrec(i1))
            !       enddo
            !       stop

            deallocate (lenrecskw, lenrecagp, recordsymagp, recordsymskw)
            deallocate (detmat)
            nbas = nbas*ipf
            ndettot = nbas
            if (.not. complexfort10) then
                allocate (detmat(nbas, nbas))
            else
                allocate (detmat(2*nbas, nbas))
            end if

        end if
        !          write(6,*) ' record found =',nrec
        !          do i=1,nrec
        !          write(6,*) i,recordsym(1:2,1:lenrec(i),i)
        !          enddo

        counter = 0
        do i = 1, nbas
            do j = 1, nbas
                !    write(6,*) i,j,detyes(i,j)
                if (detyes(i, j)) counter = counter + 1
            end do
        end do
        write (6, *) ' Total number of non zero det =', counter

    end if

    !  write(6,*) ' after  makelambda =',nsym
    !  countlines=0
    !  do i1=1,ndettot
    !     do i2=1,ndettot
    !     if(detyes(i1,i2)) countlines=countlines+1
    !     enddo
    !  enddo
    !
    !  write(6,*) ' After  makelambda ',countlines

    ntra = ntraind
    deallocate (symrot, symtra, ipsip)
    allocate (symrot(3*ntotatoms, nsym), symtra(3*ntotatoms, ntra)&
            &, forceyes(3*ntotatoms), ipsip(2*nsym*ntra))

    do i = 1, 3
        do j = 1, ntotatoms
            ind = (j - 1)*3 + i
            do k = 1, ntra
                jj = abs(atoms_maptra(j, k))
                indmap = (jj - 1)*3 + i
                symtra(ind, k) = indmap
            end do
            do k = 1, nsym
                jj = abs(atoms_maprot(j, k))
                imap = orbf_map(i, k)
                indmap = (jj - 1)*3 + abs(imap)
                if (imap .gt. 0) then
                    symrot(ind, k) = indmap
                else
                    symrot(ind, k) = -indmap
                end if
            end do
        end do
    end do

    ! write(6,*) ' Rotation symmetry'
    ! do j=1,nsym
    !  write(6,*) ' Symmetry # ',j
    !  do i=1,3*ntotatoms
    !  write(6,*) i,symrot(i,j)
    !  enddo
    ! enddo

    ! write(6,*) ' Translation symmetry'
    ! do j=1,ntra
    !  write(6,*) ' Symmetry # ',j
    !  do i=1,3*ntotatoms
    !  write(6,*) i,symtra(i,j)
    !  enddo
    ! enddo

    if (notra_forces) then
        ntrau = 1
    else
        ntrau = ntra
    end if

    if (nosym_forces) then
        nsymu = 1
    else
        nsymu = nsym
    end if

    nrecf = 1
    allocate (recordsymf(2, nsym*ntra, nrecf), lenrecf(nrecf))
    nrecf = 0

    call makeforces(ntotatoms, symrot, nsym, nsymu, symtra, ntra, ntrau, forceyes, recordsymf, lenrecf, nrecf, ipsip)

    deallocate (recordsymf, lenrecf)

    allocate (recordsymf(2, nsym*ntra, nrecf), lenrecf(nrecf))

    call makeforces(ntotatoms, symrot, nsym, nsymu, symtra, ntra, ntrau, forceyes, recordsymf, lenrecf, nrecf, ipsip)

    deallocate (ipsip, symrot, symtra, forceyes)

    !Filling the detmat
    if (.not. complexfort10) then
        call fill_detmat
    else
        call fill_detmat_complex
    end if
    if (yesalloc_jas) then
        if (noonebody) then
            allocate (jasmat(ipj*njasorb*ncell, ipj*njasorb*ncell))
        else
            allocate (jasmat(ipj*(njasorb*ncell + 1), ipj*(njasorb*ncell + 1)))
        end if
        jasmat = 0.d0
    end if
    !  write(6,*) ' Output record det '
    !  do i1=1,nrec
    !  write(6,*) lenrec(i1),recordsym(1:2,1:lenrec(i1),i1)
    !  enddo
    ! write(6,*) ' sumdet I  =',ndetorb*ncell,sum(detmat(:,:))

    !********* Symmetries of the parameters ***************

    allocate (eqdet(ndetpar))
    call par_symm(eqdet, detorb, detorbidx, detorbidxz, neqdet, ndetpar, ncell, ndetorb, onlycontrdet, ipf, symmagp)

    nrecj = 0
    njastot = 0

    if (nshelljas .ne. 0) then
        !
        ! Always use the rotations
        !
        nsym = nsym_full
        isymm(:, :, 1:nsym_full) = isymm_full(:, :, 1:nsym_full)

        apbcmap = .false.
        call update_atoms_map(apbcmap)

        if (.not. read_vecpbc) then
            vecpbc(:) = 0
        end if

        call check_symm(vecpbc, rot_jas, .false.)

        call update_atoms_map(apbcmap)

        call apply_symm_to_orbitals(jasorb, njasorb, jas_map)

        ! do is=1,nsym
        !   write(6,*) ' Symmetry =',is
        !   do i1=1,ntotatoms
        !   write(6,*) i1,atoms_maprot(i1,is)
        !   enddo
        ! enddo

        ! call generate_cellmap(cellmap_jas,atoms_map_jas)
        !
        ! call find_additional_rotations(map_itself_jas)
        !
        allocate (jasorbidx(njasorb), jasorbidxz(njasorb))
        if (write_log) write (lunit, *) ' Jastrow Data '

        if (noonebody) then
            njastot = njasorb*ncell
            allocate (symrot(njastot, nsym), symtra(njastot, ntra))
            call update_symrot(nsym, symrot, njastot, atoms_maprot, jasorb, njasorb, jas_map)
        else
            njastot = njasorb*ncell + 1
            allocate (symrot(njastot, nsym), symtra(njastot, ntra))
            call update_symrot(nsym, symrot, njastot, atoms_maprot, jasorb, njasorb, jas_map)
        end if

        ! prepare orb2atom for Jastrow orbitals
        deallocate (orb2atom)
        allocate (orb2atom(njastot))

        if (.not. yesmolatj) then
            write (6, *) ' Jastrow Orb to atom mapping '
            do i1 = 1, njasorb
                do i2 = 1, ncell
                    ionpos = jasorb(i1)%kion + (i2 - 1)*natoms
                    i3 = i1 + (i2 - 1)*njasorb
                    orb2atom(i3) = ionpos
                end do
            end do
            do i1 = 1, njasorb
                write (6, *) i1, orb2atom(i1)
            end do
            if (.not. noonebody) orb2atom(njastot) = 0
        end if

        !  write(6,*) ' Output symrot '
        !  do is=1,nsym
        !   write(6,*) ' Symmetry considered =',is
        !   do i1=1,njasorb*ncell
        !   write(6,*) i1,symrot(i1,is)
        !   enddo
        !  enddo

        !  update translation symmetries to be done
        !  ntra=1  ! only the identity
        !  allocate(symtra(njasorb*ncell+1,ntra))
        !  do is=1,njasorb*ncell+1
        !  symtra(is,1)=is
        !  enddo

        if (nosym) nsym = 1 ! Use only the identity

        call update_symtra(ntra, symtra, njastot, atoms_maptra, jasorb, njasorb)

        ! Recompute independent ions with Jastrow (in principle may be different)

        call generate_indion

        if (.not. yesmolatj) then

            const_term = njastot

            !  write(6,*) ' Input jasorb no hybrid ',njasorb
            !  do i2=1,njasorb
            !  write(6,*) i2,jasorb(i2)%ioptorb
            !  enddo

            nrecj = 1
            allocate (ipsip(2*nsym*ntra))
            allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))
            allocate (jasyes(njastot, njastot))

            nrecj = 0
            !  first run count the number of records
            write (6, *) ' before makelambda Jas  nsym ntra=', nsym, ntra
            call makelambda(njastot, symrot, nsym, symtra, ntra, jasyes &
                            , recordsymj, lenrecj, nrecj, .true., ipsip, rion, ntotatoms &
                            , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                            , .true., no_4body_jas, cut_hybridj)

            write (6, *) ' Output records Jastrow =', nrecj

            deallocate (recordsymj, lenrecj)

            allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))

            call makelambda(njastot, symrot, nsym, symtra, ntra, jasyes &
                            , recordsymj, lenrecj, nrecj, .true., ipsip, rion, ntotatoms &
                            , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                            , .true., no_4body_jas, cut_hybridj)

            if (genjason) then
                allocate (recordsymj_1(2, nsym*ntra, nrecj), lenrecj_1(nrecj))
                allocate (jasyes_1(njastot, njastot))
                allocate (jasyes_2(njastot, njastot))

                recordsymj_1 = recordsymj !checkit
                jasyes_1 = jasyes !cehckit
                lenrecj_1 = lenrecj

                nrecj_save = nrecj

                deallocate (recordsymj, jasyes, lenrecj)
                allocate (jasyes(njastot, njastot))

                nrecj = 1
                allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))

                nrecj = 0
                !  first run count the number of records
                call makelambda(njastot, symrot, nsym, symtra, ntra, jasyes &
                                , recordsymj, lenrecj, nrecj, .false., ipsip, rion, ntotatoms &
                                , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                                , .true., no_4body_jas, cut_hybridj)

                write (6, *) ' Output records Jastrow =', nrecj

                deallocate (recordsymj, lenrecj)

                allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))

                call makelambda(njastot, symrot, nsym, symtra, ntra, jasyes &
                                , recordsymj, lenrecj, nrecj, .false., ipsip, rion, ntotatoms &
                                , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                                , .true., no_4body_jas, cut_hybridj)

                allocate (recordsymj_2(2, nsym*ntra, nrecj), lenrecj_2(nrecj))

                recordsymj_2 = recordsymj
                lenrecj_2 = lenrecj
                jasyes_2 = jasyes

                deallocate (recordsymj, jasyes, lenrecj)

                !       call update_genjas

            end if ! endif of jasgenon

        else

            const_term = 1
            !  first count the basis used.
            allocate (ipsip(2*nsym*ntra))
            njasorbat = 0
            do i2 = 1, njasorb
                if (jasorb(i2)%ioptorb .eq. 900000) njasorbat = njasorbat + 1
            end do
            nbas = njasorbat
            nbas = nbas*ncell
            if (.not. noonebody) nbas = nbas + 1

            allocate (jasorbat(njasorbat))

            i3 = 0
            !  write(6,*) ' Input jasorb ',njasorb
            do i2 = 1, njasorb
                !  write(6,*) i2,jasorb(i2)%ioptorb
                if (jasorb(i2)%ioptorb .eq. 900000) then
                    i3 = i3 + 1
                    jasorbat(i3) = jasorb(i2)
                end if
            end do

            write (6, *) ' n basis Jastrow found =', nbas

            do i2 = 1, njasorbat
                write (6, *) i2, jasorbat(i2)%ioptorb, jasorbat(i2)%kion
            end do
            !  No rotation symmetry imposed.
            nsym = 1
            ntra = ntraind
            if (allocated(symrot)) deallocate (symrot)
            if (allocated(symtra)) deallocate (symtra)

            allocate (symrot(nbas, nsym), symtra(nbas, ntra))
            !  Identity for symrot

            do i2 = 1, nbas
                symrot(i2, 1) = i2
            end do

            call update_symtraj(ntra, symtra, nbas, atoms_maptra, jasorbat, njasorbat)
            !  if(allocated(detyes)) deallocate(detyes)

            if (allocated(jasyes)) deallocate (jasyes)
            allocate (jasyes(nbas, nbas))
            njastot = nbas

            nrecj = 1
            if (allocated(ipsip)) deallocate (ipsip)
            allocate (ipsip(2*nsym*ntra))
            if (allocated(recordsymj)) deallocate (recordsymj)
            if (allocated(lenrecj)) deallocate (lenrecj)

            allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))
            nrecj = 0
            call makelambda(nbas, symrot, nsym, symtra, ntra, jasyes &
                            , recordsymj, lenrecj, nrecj, .true., ipsip, rion, ntotatoms &
                            , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                            , .false., no_4body_jas, cut_hybridj)
            write (6, *) ' Output records jas =', nrecj
            deallocate (recordsymj, lenrecj)
            allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))
            call makelambda(nbas, symrot, nsym, symtra, ntra, jasyes &
                            , recordsymj, lenrecj, nrecj, .true., ipsip, rion, ntotatoms &
                            , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                            , .false., no_4body_jas, cut_hybridj)

            if (genjason) then
                allocate (recordsymj_1(2, nsym*ntra, nrecj), lenrecj_1(nrecj))
                allocate (jasyes_1(njastot, njastot))
                allocate (jasyes_2(njastot, njastot))

                recordsymj_1 = recordsymj !checkit
                jasyes_1 = jasyes !cehckit
                lenrecj_1 = lenrecj

                nrecj_save = nrecj

                deallocate (recordsymj, jasyes, lenrecj)
                allocate (jasyes(njastot, njastot))

                nrecj = 1
                allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))

                nrecj = 0
                !  first run count the number of records
                call makelambda(nbas, symrot, nsym, symtra, ntra, jasyes &
                                , recordsymj, lenrecj, nrecj, .false., ipsip, rion, ntotatoms &
                                , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                                , .false., no_4body_jas, cut_hybridj)

                write (6, *) ' Output records Jastrow =', nrecj

                deallocate (recordsymj, lenrecj)

                allocate (recordsymj(2, nsym*ntra, nrecj), lenrecj(nrecj))

                call makelambda(nbas, symrot, nsym, symtra, ntra, jasyes &
                                , recordsymj, lenrecj, nrecj, .false., ipsip, rion, ntotatoms &
                                , orb2atom, cellscale, deps, .false., .false., zeta, J3_off &
                                , .false., no_4body_jas, cut_hybridj)

                allocate (recordsymj_2(2, nsym*ntra, nrecj), lenrecj_2(nrecj))

                recordsymj_2 = recordsymj
                lenrecj_2 = lenrecj
                jasyes_2 = jasyes

                deallocate (recordsymj, jasyes, lenrecj)

                !       call update_genjas

            end if ! endif of jasgenon

        end if
        deallocate (ipsip, atoms_maprot, atoms_maptra, symrot, symtra)

        call generate_orbidx(zeta, natoms, jasorbidx, jasorb, njasorb, indion)
        if (.not. eqatoms) then
            jasorbidxz = jasorbidx
        else
            call generate_orbidx(zeta, natoms, jasorbidxz, jasorb, njasorb, indionz)
        end if
        allocate (eqjas(njaspar))
        nosym_contr = nosym_contrj
        call par_symm(eqjas, jasorb, jasorbidx, jasorbidxz, neqjas, njaspar &
                      , ncell, njasorb, onlycontrjas, 1, symmagp)
    end if

    call write_xsf

    call write_fort10

    if (write_log) close (lunit)

    if (allocated(J3_off)) deallocate (J3_off)
    if (allocated(cut_hybrid)) deallocate (cut_hybrid, cut_hybridj)
    if (allocated(jasmat)) deallocate (jasmat, atomic_jasmat)
    !================= END makefort10 ===================!

contains

    subroutine write_xsf
        implicit none
        integer, parameter :: xsfunit = 21
        integer j
        open (unit=xsfunit, file="structure.xsf", form="formatted", status="unknown")
        if (pbcfort10) then
            write (xsfunit, *) 'CRYSTAL '
            write (xsfunit, *) 'PRIMVEC '
            if (yes_tilted) then
                write (xsfunit, *) (nxyz(1)*newcell(j, 1)*length_unit, j=1, 3)
                write (xsfunit, *) (nxyz(2)*newcell(j, 2)*length_unit, j=1, 3)
                write (xsfunit, *) (nxyz(3)*newcell(j, 3)*length_unit, j=1, 3)
            else
                write (xsfunit, *) cellscale(1)*length_unit, 0.d0, 0.d0
                write (xsfunit, *) 0.d0, cellscale(2)*length_unit, 0.d0
                write (xsfunit, *) 0.d0, 0.d0, cellscale(3)*length_unit
            end if
            write (xsfunit, *) '# empty line'
            write (xsfunit, *) 'PRIMCOORD ', 1
            write (xsfunit, *) natoms*ncell, 1
        else
            write (xsfunit, *) 'ATOMS '
        end if
        do i1 = 1, ntotatoms
            write (xsfunit, '(i3,6(F16.8,4X))') int(zeta(1, i1)), rion(:, i1)*length_unit
        end do
        close (xsfunit)
    end subroutine write_xsf

    !****** Find rotations using PWSCF subroutines ************

    subroutine find_rotations
        implicit none
        double precision, parameter :: eps = 0.01

        at(:, :) = 0.d0
        at(1, 1) = 1.d0
        at(2, 2) = 1.d0
        at(3, 3) = 1.d0

        bg = 0.d0
        bg(1, 1) = 1.d0/at(1, 1)
        bg(2, 2) = 1.d0/at(2, 2)
        bg(3, 3) = 1.d0/at(3, 3)

        nrot = 48

        call cubicsym(at, isymm, isname, nrot)

        write (*, *) 'Number of Cell Symmetries : ', nrot

        do i1 = 1, nrot
            write (6, *) ' Sym op 3x3 matrix ', i1
            do is = 1, 3
                write (6, *) isymm(1:3, is, i1)
            end do
        end do

        nsym = nrot

        write (*, *) 'Number of allowed symmetries: ', nsym
        if (write_log) write (lunit, *) 'Number of allowed symmetries: ', nsym

        open (file="symmetries.dat", unit=15, status="unknown", form="formatted")
        write (15, *) nsym
        do i1 = 1, nsym
            write (15, *) i1, ') ', isname(i1)
            write (15, '(3i3,a,3i3)') (isymm(j, 1, i1), j=1, 3)
            write (15, '(3i3)') (isymm(j, 2, i1), j=1, 3)
            write (15, '(3i3)') (isymm(j, 3, i1), j=1, 3)
            if (write_log) then
                write (lunit, *) i1, ')', isname(i1)
                write (lunit, '(3i3,a,3i3)') (isymm(j, 1, i1), j=1, 3)
                write (lunit, '(3i3)') (isymm(j, 2, i1), j=1, 3)
                write (lunit, '(3i3)') (isymm(j, 3, i1), j=1, 3)
            end if
        end do

        close (15)

        nsym_full = nsym
        isymm_full(:, :, 1:nsym_full) = isymm(:, :, 1:nsym_full)

    end subroutine find_rotations

    subroutine update_atoms_map(apbcmap)
        implicit none
        integer apbcs, j, is, i2, npip(3)
        logical found, apbcmap(3)
        real*8 rionmap(3), dist(3), check1, check2
        logical, dimension(:), allocatable :: yessite
        npip = 0
        atoms_maprot = 0
        atoms_maptra = 0

        write (6, *) 'Atom coordinates'
        do j = 1, ntotatoms
            write (6, *) j, rion(:, j)
            dist(:) = rion(:, j) - rion(:, 1)
            call CartesianToCrystal(dist, 1)
        end do

        do j = 1, ntotatoms
            do is = 1, nsym
                rionmap(1:3) = matmul(isymm(:, :, is), rion(:, j))
                found = .true.
                i2 = 1

                do while (i2 .le. ntotatoms .and. found)
                    dist(:) = rionmap(:) - rion(:, i2)
                    call CartesianToCrystal(dist, 1)
                    if (pbcfort10) then
                        apbcs = 1
                        do i3 = 1, 3
                            !                  itest=makenpipo(dist(i3),cellscale(i3),deps)
                            npip(i3) = makenpip(dist(i3), cellscale(i3), deps)
                            !                  itest=nint(dist(i3)/cellscale(i3))
                            if ((npip(i3)/2)*2 .ne. npip(i3) .and. apbcmap(i3)) then
                                apbcs = -apbcs
                            end if
                        end do
                        !               call makeimageo(dist,cellscale,deps)
                        call makeimage(dist, cellscale, deps)
                        !               dist(:)=dist(:)-cellscale(:)*anint(dist(:)/cellscale(:))
                    else
                        apbcs = 1
                    end if
                    check1 = zeta(1, i2)
                    if (eq_intatoms) check1 = int(check1)
                    check2 = zeta(1, j)
                    if (eq_intatoms) check2 = int(check2)
                    if (sum(abs(dist(:))) .lt. deps .and. check1 .eq. check2) then
                        atoms_maprot(j, is) = apbcs*i2
                        found = .false.
                    end if
                    i2 = i2 + 1
                end do
                !              if(found) then
                !              write(6,*) ' Map not found rotation site/transform  !!! ',j,is
                !              stop
                !              endif
            end do
        end do

        !    check surrettivity
        allocate (yessite(ntotatoms))
        !   Identity
        do j = 1, ntotatoms
            atoms_maprot(j, 1) = j
        end do

        do is = 1, nsym
            yessite = .false.
            do j = 1, ntotatoms
                i2 = abs(atoms_maprot(j, is))
                if (i2 .ne. 0) yessite(i2) = .true.
            end do
            found = .true.
            do j = 1, ntotatoms
                if (.not. yessite(j)) found = .false.
            end do
            if (.not. found) then
                ! Eliminate this symmetry from the chosen ones.
                atoms_maprot(:, is) = 0
            end if
        end do
        deallocate (yessite)

        do j = 1, ntotatoms
            do is = 1, ntraind
                rionmap(1:3) = rion(:, j) + trasl(:, is)
                found = .true.
                i2 = 1

                do while (i2 .le. ntotatoms .and. found)
                    dist(:) = rionmap(:) - rion(:, i2)
                    call CartesianToCrystal(dist, 1)
                    do i3 = 1, 3
                    end do
                    if (pbcfort10) then
                        apbcs = 1
                        do i3 = 1, 3
                            !                  itest=nint(dist(i3)/cellscale(i3))
                            !                  itest=makenpipo(dist(i3),cellscale(i3),deps)
                            npip(i3) = makenpip(dist(i3), cellscale(i3), deps)
                            if ((npip(i3)/2)*2 .ne. npip(i3) .and. apbcmap(i3)) then
                                apbcs = -apbcs
                            end if
                        end do
                        !               call makeimageo(dist,cellscale,deps)
                        call makeimage(dist, cellscale, deps)
                        !               dist(:)=dist(:)-cellscale(:)*anint(dist(:)/cellscale(:))
                    else
                        apbcs = 1
                    end if
                    check1 = zeta(1, i2)
                    if (eq_intatoms) check1 = int(check1)
                    check2 = zeta(1, j)
                    if (eq_intatoms) check2 = int(check2)
                    if (sum(abs(dist(:))) .lt. deps .and. check1 .eq. check2) then
                        atoms_maptra(j, is) = apbcs*i2
                        found = .false.
                    end if
                    i2 = i2 + 1
                end do
                !     if(found) then
                !     write(6,*) ' Map not found translation  site/transform  !!! ',j,is
                !     stop
                !     endif
            end do
        end do

    end subroutine update_atoms_map

    subroutine generate_indion
        implicit none
        integer i, j, k, ind, imap, jmap
        real*8 check1, check2
        if (.not. allocated(indion)) allocate (indion(natoms), indionz(natoms))
        !  Using all possible symmetries of the lattice compute the independent
        !  atoms

        do k = 1, natoms
            indion(k) = k
            do i = 1, ntra
                imap = abs(atoms_maptra(k, i))
                do j = 1, nsym
                    jmap = abs(atoms_maprot(imap, j))
                    if (jmap .ne. 0 .and. jmap .lt. indion(k)) indion(k) = jmap
                end do
            end do
        end do
        if (.not. eqatoms) then
            indionz = indion
        else
            indionz = 0
            ind = 0
            do k = 1, natoms
                if (indionz(k) .eq. 0) then
                    indionz(k) = k
                    do j = 1, natoms
                        check1 = zeta(1, j)
                        if (eq_intatoms) check1 = int(check1)
                        check2 = zeta(1, k)
                        if (eq_intatoms) check2 = int(check2)
                        if (check1 .eq. check2) indionz(j) = k
                    end do
                end if
            end do
        end if

    end subroutine generate_indion

    subroutine write_fort10
        implicit none
        integer, parameter :: ufort10 = 10
        double precision :: rs, zetaw
        integer :: i1, i2, i3, iskipcell, itothyb
        integer :: iesfree, iessw
        integer :: good, gooddet
        logical, dimension(:), allocatable :: yeskion
        logical cmplx_coeff
        character(len=30) gformat

        allocate (yeskion(natoms))
        ! volume_new=abs(newcell(1,1)*(newcell(2,2)*newcell(3,3)-newcell(2,3)*newcell(3,2)) &
        ! -newcell(1,2)*(newcell(2,1)*newcell(3,3)-newcell(2,3)*newcell(3,1))&
        ! +newcell(1,3)*(newcell(2,1)*newcell(3,2)-newcell(2,2)*newcell(3,1)))

        rs = (3.d0*volume_new*product(nxyz(:))/(4.d0*sum(zeta(2, :)*PI)))**(1.d0/3.d0)
        write (*, *) ' Effective rs = ', rs

        iessw = nrec
        gooddet = 0
        do i1 = 1, ndettot
            do i2 = 1, ndettot
                if (detyes(i1, i2)) gooddet = gooddet + 1
            end do
        end do

        !  given at compute alphap betap and gammap as defined in crystallography

        !   alphap = angle  from b and c
        alphap = dacos(sum(newcell(:, 2)*newcell(:, 3))/smallcell(2)/smallcell(3))
        !   betap = angle from c and a
        betap = dacos(sum(newcell(:, 3)*newcell(:, 1))/smallcell(3)/smallcell(1))
        !   gammap = angle between  a  and b
        gammap = dacos(sum(newcell(:, 1)*newcell(:, 2))/smallcell(1)/smallcell(2))

        open (unit=ufort10, file="fort.10_new", status="unknown", form="formatted")
        if (pbcfort10) then
            if (complexfort10 .or. yes_crystal) then
                if (yes_tilted) then
                    write (ufort10, *) '# PBC_T a b c phase up phase down'
                    write (ufort10, '(9f18.12,6f12.8)') (nxyz(j)*newcell(:, j), j=1, 3) &
                            &, (phase(j), j=1, 3), (phasedo(j), j=1, 3)
                else
                    write (ufort10, *) '# PBC_C rs, Ly/Lx, Lz/Lx'
                    write (ufort10, '(3f18.12,6f12.8)') rs, celldm(2)*dble(nxyz(2))/dble(nxyz(1)) &
                        , celldm(3)*dble(nxyz(3))/dble(nxyz(1)), (phase(j), j=1, 3), (phasedo(j), j=1, 3)
                end if
            else
                write (ufort10, *) '# PBC rs, Ly/Lx, Lz/Lx'
                write (ufort10, '(3f18.12,3f12.8)') rs, celldm(2)*dble(nxyz(2))/dble(nxyz(1)) &
                    , celldm(3)*dble(nxyz(3))/dble(nxyz(1)), (phase(j), j=1, 3)
            end if
        end if

        write (ufort10, *) '# Nelup  #Nel  # Ion'

        nel = nel*ncell

        if (neldiff == -1) then
            nelup = (nel + 1)/2
            if (ipf .eq. 2) then
                neldiff = 2*numpaired
                if (mod(nel, 2) .ne. 0) neldiff = 1 + 2*numpaired
            else
                if (nelup*2 .ne. nel) then
                    neldiff = 1
                else
                    neldiff = 0
                end if
            end if
        else
            if (mod(nel, 2) .eq. 0) then
                nelup = nel/2
                if (mod(neldiff, 2) .eq. 0) then
                    nelup = nelup + neldiff/2
                else
                    write (6, *) ' Warning Total charge not neutral '
                    nel = nel + 1
                    !         now nel is odd
                    nelup = nelup + (neldiff + 1)/2
                end if
                if (ipf .eq. 2) neldiff = 2*numpaired
            else
                nelup = nel/2
                if (mod(neldiff, 2) .ne. 0) then
                    nelup = nelup + (neldiff + 1)/2
                else
                    write (6, *) ' Total charge not neutral '
                    nel = nel + 1
                    !         now nel is even
                    nelup = nel/2 + neldiff/2
                end if
                if (ipf .eq. 2) neldiff = 1 + 2*numpaired
            end if
        end if

        if (ipf .eq. 2 .and. neldiff .gt. 1) then
            npar_eagp = (neldiff*(neldiff - 1))/2
        else
            npar_eagp = 0
        end if

        if (genjason .and. njastot .ne. 0) call update_genjas

        iesfree = nrecj
        good = 0
        if (iesfree .ne. 0) then
            do i1 = 1, njastot
                do i2 = 1, njastot
                    if (jasyes(i1, i2)) good = good + 1
                end do
            end do
        end if

        !If ipf=2 then nel=-nel in the fort.10
        if (ipf .eq. 2) then
            nelup = numpaired*nel + nelup
            if (nodownpfaff_read .and. symmagp) nelup = -nelup
            if (yesbump) then
                write (ufort10, *) nelup, -nel, -natoms*ncell
            else
                write (ufort10, *) nelup, -nel, natoms*ncell
            end if
            if (nodownpfaff_read .and. symmagp) nelup = -nelup
        else
            if (yesbump) then
                write (ufort10, *) nelup, nel, -natoms*ncell
            else
                write (ufort10, *) nelup, nel, natoms*ncell
            end if
        end if

        if (yesalloc_jas) call fill_jasmat

        write (ufort10, *) '# Shell Det.   # Shell Jas. '
        if (complexfort10) then
            if (noonebody) then
                if (yes_crystalj) then
                    write (ufort10, *) - nshelldet*ncell, -nshelljas*ncell
                else
                    write (ufort10, *) - nshelldet*ncell, nshelljas*ncell
                end if
            else
                if (yes_crystalj) then
                    write (ufort10, *) - nshelldet*ncell, -(nshelljas*ncell + 1)
                else
                    write (ufort10, *) - nshelldet*ncell, nshelljas*ncell + 1
                end if
            end if
        else
            if (noonebody) then
                if (yes_crystalj) then
                    write (ufort10, *) nshelldet*ncell, -nshelljas*ncell
                else
                    write (ufort10, *) nshelldet*ncell, nshelljas*ncell
                end if
            else
                if (yes_crystalj) then
                    write (ufort10, *) nshelldet*ncell, -(nshelljas*ncell + 1)
                else
                    write (ufort10, *) nshelldet*ncell, nshelljas*ncell + 1
                end if
            end if
        end if
        if (npar_eagp .gt. 0) then
            write (ufort10, *) '# Jas 2body  # Det   #  3 body atomic par.  iessw_eagp #npar '
            write (ufort10, *) twobody, ndetpar*ncell, njaspar*ncell, npar_eagp, npar_eagp
        else
            write (ufort10, *) '# Jas 2body  # Det   #  3 body atomic par. '
            write (ufort10, *) twobody, ndetpar*ncell, njaspar*ncell
        end if

        write (ufort10, *) '# Det mat. =/0  # Jas mat. =/0'

        if (readunpaired) then
            nlambda_unp = neldiff*norb_unpaired
        else
            ! add the counting of unpaired hybrid orbitals
            nlambda_unp = ndettot*neldiff
            !Check unpaired with Sandro
        end if

        if (forcesymm) then
            write (ufort10, *) - gooddet - nlambda_unp, good
        else
            write (ufort10, *) gooddet + nlambda_unp, good
        end if
        write (ufort10, *) '# Eq. Det atomic par. # Eq. 3 body atomic. par.'
        write (ufort10, *) neqdet, neqjas
        write (ufort10, *) '# unconstrained iesfree,iessw,ieskinr,I/O flag'

        write (ufort10, *) iesfree, iessw + nlambda_unp, nrecf, 0
        write (ufort10, *) '# Ion coordinates'
        do i1 = 1, ntotatoms
            zetaw = zeta(1, i1)
            if ((int(zeta(1, i1)) .ne. 1 .and. int(zeta(1, i1)) .ne. 2 &
                 .and. int(zeta(1, i1)) .ne. 3 .and. int(zeta(1, i1)) .ne. 4) &
                .or. nopseudo) then
                zetaw = int(zeta(1, i1))
            end if
            write (ufort10, '(1f6.0,1f8.2,3f22.14)') zeta(2, i1), zetaw, (rion(j, i1), j=1, 3)
        end do
        write (ufort10, *) '#  Constraints for forces: ion - coordinate'
        do i1 = 1, nrecf
            write (ufort10, *) lenrecf(i1), (recordsymf(1, i2, i1), recordsymf(2, i2, i1)&
                    &, i2=1, lenrecf(i1))
        end do
        write (ufort10, *) '#          Parameters Jastrow two body'

        ! added by Kosuke Nakano on 18 Feb. 2021.
        ! check if niesd is consistent with what the code expects if it is explicitly defined.
        if (niesd .ne. -1) then
            if ((n_onebody + n_twobody) .eq. niesd) then
                write (6, *) " niesd is consistent with what the code expects"
                write (6, *) " niesd=", niesd, "expected=", n_onebody + n_twobody
                niesd = n_onebody + n_twobody
            else
                write (6, *) " Error!! niesd is inconsistent with what the code expects"
                write (6, *) " niesd=", niesd, "expected=", n_onebody + n_twobody
                stop
            end if
        else
            niesd = n_onebody + n_twobody
        end if

        if (symmagp) then
            if (tf == 0) then
                write (ufort10, *) 0
            else
                ! added by K.N. on 18.Feb
                if (n_onebody .eq. 0) then
                    write (ufort10, *) niesd, (twobodypar(itwo), itwo=1, n_twobody)
                else
                    write (ufort10, *) niesd, (twobodypar(itwo), itwo=1, n_twobody), (onebodypar(ione), ione=1, n_onebody)
                end if
                !if(onebodypar.eq.0.d0) then
                !    write(ufort10, *) 1, twobodypar
                !else
                !    write(ufort10, *) niesd, twobodypar, (onebodypar, i1 = 1, niesd - 1)
                !endif
            end if
        else
            if (tf == 0) then
                write (ufort10, *) - 1
            else
                if (n_onebody .eq. 0) then
                    write (ufort10, *) - niesd, (twobodypar(itwo), itwo=1, n_twobody)
                else
                    write (ufort10, *) - niesd, (twobodypar(itwo), itwo=1, n_twobody), (onebodypar(ione), ione=1, n_onebody)
                end if
                !if(onebodypar.eq.0.d0) then
                !    write(ufort10, *) -1, twobodypar
                !else
                !    write(ufort10, *) -niesd, twobodypar, (onebodypar, i1 = 1, niesd - 1)
                !endif
            end if
        end if

#ifdef __KCOMP
        gformat = "(i12,32767e26.18)"
#else
        gformat = "(i12,1000000000e26.18)"
#endif

        ! writing basis set for each atom
        write (ufort10, *) '#          Parameters atomic wf'
        cmplx_coeff = .false.
        do i2 = 1, ncell
            if (detorb(1)%ioptorb .ne. 900000) then
                if (complexfort10 .and. detorb(1)%nparm .gt. 1) cmplx_coeff = .true.
                write (ufort10, *) detorb(1)%itype, detorb(1)%nparm, detorb(1)%ioptorb
                if (.not. cmplx_coeff) then
                    write (ufort10, gformat) detorb(1)%kion + (i2 - 1)*natoms, (detorb(1)%parms(i3), i3=1, detorb(1)%nparm)
                else
                    write (ufort10, gformat) detorb(1)%kion + (i2 - 1)*natoms, &
                        (detorb(1)%parms(i3), i3=1, 3*detorb(1)%nparm/2)
                    !        (detorb(1)%parms(i3),0.,i3=detorb(1)%nparm/2+1,detorb(1)%nparm)
                end if
                !#ifdef __KCOMP
                !write(ufort10,'(i3,32767e26.18)') detorb(1)%kion+(i2-1)*natoms,(detorb(1)%parms(i3),i3=1,detorb(1)%nparm)
                !#else
                !write(ufort10,'(i3,1000000e26.18)') detorb(1)%kion+(i2-1)*natoms,(detorb(1)%parms(i3),i3=1,detorb(1)%nparm)
                !#endif
            end if
            do i1 = 2, ndetorb
                cmplx_coeff = .false.
                if (detorb(i1)%shell .ne. detorb(i1 - 1)%shell) then
                    if (detorb(i1)%ioptorb .ne. 900000) then
                        if (complexfort10 .and. detorb(i1)%nparm .gt. 1) cmplx_coeff = .true.
                        write (ufort10, *) detorb(i1)%itype, detorb(i1)%nparm, detorb(i1)%ioptorb
                        if (.not. cmplx_coeff) then
                            write (ufort10, gformat) detorb(i1)%kion + (i2 - 1)*natoms &
                                , (detorb(i1)%parms(i3), i3=1, detorb(i1)%nparm)
                        else
                            write (ufort10, gformat) detorb(i1)%kion + (i2 - 1)*natoms, &
                                (detorb(i1)%parms(i3), i3=1, 3*detorb(i1)%nparm/2)
                            !         &(detorb(i1)%parms(i3),0.,i3=detorb(i1)%nparm/2+1,detorb(i1)%nparm)
                        end if
                        !#ifdef __KCOMP
                        !write(ufort10,'(i3,32767e26.18)') detorb(i1)%kion+(i2-1)*natoms &
                        !#else
                        !write(ufort10,'(i3,1000000e26.18)') detorb(i1)%kion+(i2-1)*natoms &
                        !#endif
                        !& ,(detorb(i1)%parms(i3),i3=1,detorb(i1)%nparm)
                    end if
                end if
            end do
        end do

        ! writing hybrid orbitals
        ! COMPLEX DEB bisogna generalizzare al complesso anche gli orbitali ibridi????
        iskipcell = 0
        do i2 = 1, ncell
            itothyb = 0
            yeskion = .true.
            if (detorb(1)%ioptorb .eq. 900000) then
                write (ufort10, *) detorb(1)%itype, detorb(1)%nparm, detorb(1)%ioptorb

                if (complexfort10) then
                    write (ufort10, *) detorb(1)%kion + (i2 - 1)*natoms, &
                            &(nint(detorb(1)%parms(i3)) + iskipcell, i3=1, detorb(1)%nparm/2),&
                    &(real(detorb(1)%parms(i3)), i3=detorb(1)%nparm/2 + 1, 3*detorb(1)%nparm/2)
                else
                    write (ufort10, *) detorb(1)%kion + (i2 - 1)*natoms, &
                            &(nint(detorb(1)%parms(i3)) + iskipcell, i3=1, detorb(1)%nparm/2),&
                    &(real(detorb(1)%parms(i3)), i3=detorb(1)%nparm/2 + 1, detorb(1)%nparm)
                end if
                !      On a given atom the dimension of the space is given by the length of the
                !      hybrid atomic orbital
                if (yeskion(detorb(1)%kion)) then
                    itothyb = itothyb + detorb(1)%nparm/2
                    yeskion(detorb(1)%kion) = .false.
                end if
            end if
            do i1 = 2, ndetorb
                if (detorb(i1)%shell .ne. detorb(i1 - 1)%shell) then
                    if (detorb(i1)%ioptorb .eq. 900000) then
                        write (ufort10, *) detorb(i1)%itype, detorb(i1)%nparm, detorb(i1)%ioptorb
                        if (complexfort10) then
                            write (ufort10, *) detorb(i1)%kion + (i2 - 1)*natoms &
                                    &, (nint(detorb(i1)%parms(i3)) + iskipcell, i3=1, detorb(i1)%nparm/2)&
                                    &, (real(detorb(i1)%parms(i3)), i3=detorb(i1)%nparm/2 + 1, 3*detorb(i1)%nparm/2)
                        else
                            write (ufort10, *) detorb(i1)%kion + (i2 - 1)*natoms &
                                    &, (nint(detorb(i1)%parms(i3)) + iskipcell, i3=1, detorb(i1)%nparm/2)&
                                    &, (real(detorb(i1)%parms(i3)), i3=detorb(i1)%nparm/2 + 1, detorb(i1)%nparm)
                        end if
                        if (yeskion(detorb(i1)%kion)) then
                            itothyb = itothyb + detorb(i1)%nparm/2
                            yeskion(detorb(i1)%kion) = .false.
                        end if
                    end if
                end if
            end do
            iskipcell = iskipcell + itothyb/ipf
        end do

        i2 = 0
        do i1 = 1, natoms
            if (.not. yeskion(i1)) i2 = i2 + 1
        end do
        if (i2 .ne. 0 .and. i2 .ne. natoms) then
            write (6, *) ' Warning  all atoms should have hybrid orbitals!!!'
        end if

        write (ufort10, *) '#  Parameters atomic Jastrow wf'
        if (nshelljas .ne. 0) then
            !
            do i2 = 1, ncell
                if (jasorb(1)%ioptorb .ne. 900000) then
                    write (ufort10, *) jasorb(1)%itype, jasorb(1)%nparm, jasorb(1)%ioptorb
#ifdef __KCOMP
                    write (ufort10, '(i12,32767e26.18)') jasorb(1)%kion + (i2 - 1)*natoms &
                        , (jasorb(1)%parms(i3), i3=1, jasorb(1)%nparm)
#else
                    write (ufort10, '(i12,1000000e26.18)') jasorb(1)%kion + (i2 - 1)*natoms &
                        , (jasorb(1)%parms(i3), i3=1, jasorb(1)%nparm)
#endif
                end if
                do i1 = 2, njasorb
                    if (jasorb(i1)%shell .ne. jasorb(i1 - 1)%shell) then
                        if (jasorb(i1)%ioptorb .ne. 900000) then
                            write (ufort10, *) jasorb(i1)%itype, jasorb(i1)%nparm, jasorb(i1)%ioptorb
#ifdef __KCOMP
                            write (ufort10, '(i12,32767e26.18)') jasorb(i1)%kion + (i2 - 1)*natoms &
                            & , (jasorb(i1)%parms(i3), i3=1, jasorb(i1)%nparm)
#else
                            write (ufort10, '(i12,1000000e26.18)') jasorb(i1)%kion + (i2 - 1)*natoms &
                            & , (jasorb(i1)%parms(i3), i3=1, jasorb(i1)%nparm)
#endif
                        end if
                    end if
                end do
            end do
            ! One-Body Term
            if (.not. noonebody) then
                write (ufort10, '(3i4)') 1, 0, 200
                write (ufort10, '(i3)') 1
            end if
            iskipcell = 0
            do i2 = 1, ncell
                itothyb = 0
                yeskion = .true.
                if (jasorb(1)%ioptorb .eq. 900000) then
                    write (ufort10, *) jasorb(1)%itype, jasorb(1)%nparm, jasorb(1)%ioptorb
#ifdef __KCOMP
                    write (ufort10, '(i12,32767e26.18)') jasorb(1)%kion + (i2 - 1)*natoms,&
                         &(nint(jasorb(1)%parms(i3)) + iskipcell, i3=1, jasorb(1)%nparm/2),&
                         &(real(jasorb(1)%parms(i3)), i3=jasorb(1)%nparm/2 + 1, jasorb(1)%nparm)
#else
                    write (ufort10, '(i12,1000000e26.18)') jasorb(1)%kion + (i2 - 1)*natoms, &
                            &(nint(jasorb(1)%parms(i3)) + iskipcell, i3=1, jasorb(1)%nparm/2),&
                    &(real(jasorb(1)%parms(i3)), i3=jasorb(1)%nparm/2 + 1, jasorb(1)%nparm)
#endif
                    if (yeskion(jasorb(1)%kion)) then
                        itothyb = itothyb + jasorb(1)%nparm/2
                        yeskion(jasorb(1)%kion) = .false.
                    end if
                end if
                do i1 = 2, njasorb
                    if (jasorb(i1)%shell .ne. jasorb(i1 - 1)%shell) then
                        if (jasorb(i1)%ioptorb .eq. 900000) then
                            write (ufort10, *) jasorb(i1)%itype, jasorb(i1)%nparm, jasorb(i1)%ioptorb
                            write (ufort10, *) jasorb(i1)%kion + (i2 - 1)*natoms &
                                    &, (nint(jasorb(i1)%parms(i3)) + iskipcell, i3=1, jasorb(i1)%nparm/2) &
                                    &, (real(jasorb(i1)%parms(i3)), i3=jasorb(i1)%nparm/2 + 1, jasorb(i1)%nparm)
                            if (yeskion(jasorb(i1)%kion)) then
                                itothyb = itothyb + jasorb(i1)%nparm/2
                                yeskion(jasorb(i1)%kion) = .false.
                            end if
                        end if
                    end if
                end do
                iskipcell = iskipcell + itothyb
            end do
            !
            !
        end if
        write (ufort10, *) '#  Occupation atomic orbitals'
        if (yesmolat) then
            do i2 = 1, ncell
                do i1 = 1, ndetorb
                    if (detorb(i1)%ioptorb .ne. 900000) write (ufort10, *) 0
                end do
            end do
            do i2 = 1, ncell
                do i1 = 1, ndetorb
                    if (detorb(i1)%ioptorb .eq. 900000) write (ufort10, *) 1
                end do
            end do
        else
            do i1 = 1, ndetorb*ncell
                write (ufort10, *) 1
            end do
        end if

        write (ufort10, *) '#  Occupation atomic orbitals  Jastrow'

        if (yesmolatj) then
            do i2 = 1, ncell
                do i1 = 1, njasorb
                    if (jasorb(i1)%ioptorb .ne. 900000 .and. jasorb(i1)%ioptorb .ne. 200) &
                            & write (ufort10, *) 0
                end do
            end do
            do i2 = 1, ncell
                do i1 = 1, njasorb
                    if (jasorb(i1)%ioptorb .eq. 900000) write (ufort10, *) 1
                end do
            end do
        else
            do i1 = 1, njasorb*ncell
                write (ufort10, *) 1
            end do
        end if

        write (*, *) ' Number orb det ', ndetorb*ncell
        write (*, *) ' Number orb jas ', njasorb*ncell

        ! One-Body Term
        if (.not. noonebody) write (ufort10, *) 1

        ! write(6,*) ' sumdet III  =',ndetorb*ncell,sum(detmat(:,:))

        write (ufort10, *) '#          Nonzero values of  detmat'
        ! if complex w.f. write real part and imaginary part of each element of detmat
        if (.not. complexfort10) then
            do i1 = 1, ndettot
                do i2 = 1, ndettot
                    if (detyes(i1, i2)) write (ufort10, *) i1, i2, detmat(i1, i2)
                end do
            end do
        else
            do i1 = 1, 2*ndettot, 2
                do i2 = 1, ndettot
                    if (detyes(i1/2 + 1, i2)) write (ufort10, *) i1/2 + 1, i2, detmat(i1, i2), detmat(i1 + 1, i2)
                end do
            end do

        end if
        ! COMPLEX DEB to modify
        ! complex w.f. for unpaired orbitals is not implemented yet
        if (readunpaired) then
            do i2 = 1, ncell
                do i1 = 1, norb_unpaired
                    write (ufort10, *) orb_numbers(i1) + (i2 - 1)*ndetorb, ncell*ndetorb + 1, 1.d0/dble(norb_unpaired*ncell)
                end do
            end do
        elseif (neldiff .gt. 0) then
            ! add the unpaired hybrid orbital parameters
            do i2 = 1, neldiff
                if (complexfort10) then
                    do i1 = 1, ndettot
                        if (i2 == i1) then
                            write (ufort10, *) i1, ndettot + i2, 1.d0, 0.d0
                        else
                            write (ufort10, *) i1, ndettot + i2, 0.d0, 0.d0
                        end if
                    end do
                else
                    do i1 = 1, ndettot
                        if (i2 == i1) then
                            write (ufort10, *) i1, ndettot + i2, 1.d0
                        else
                            write (ufort10, *) i1, ndettot + i2, 0.d0
                        end if
                    end do
                end if
            end do
            write (6, *) neldiff, ' Unpaired orbital parameters are added with no symmetry !'
        end if

        if (npar_eagp .ne. 0) then
            do iy = 1, neldiff
                do ix = iy + 1, neldiff
                    if (complexfort10) then
                        write (ufort10, *) ix, iy, 0.d0, 0.d0
                    else
                        write (ufort10, *) ix, iy, 0.d0
                    end if
                end do
            end do
        end if

        write (ufort10, *) '#  Grouped par. in the chosen ordered basis'
!       if(ndetorb * ncell.lt.100) then
!           gformat = "(1000i5)"
!       elseif(ndetorb * ncell.lt.1000) then
!           gformat = "(1000i6)"
!       elseif(ndetorb * ncell.lt.10000) then
!           gformat = "(10000i7)"
!       else
#ifdef __KCOMP
        gformat = "(32767i12)"
#else
        gformat = "(1000000000i12)"
#endif
!       endif

        do i1 = 1, nrec
            write (ufort10, gformat) lenrec(i1), recordsym(1:2, 1:lenrec(i1), i1)
        end do
        if (readunpaired) then
            do i2 = 1, ncell
                do i1 = 1, norb_unpaired
                    write (ufort10, gformat) 1, orb_numbers(i1) + (i2 - 1)*ndetorb, ncell*ndetorb + 1
                end do
            end do
        elseif (neldiff .gt. 0) then
            ! add the unpaired hybrid orbital symmetry tables
            do i2 = 1, neldiff
                do i1 = 1, ndettot
                    write (ufort10, *) 1, i1, ndettot + i2
                end do
            end do
            write (6, *) neldiff, ' Unpaired  orbital symmetry tables are added with no symmetry!'
            write (6, *) ' Use readunpaired in case  you need a more efficient fort.10 !'
        end if
        if (npar_eagp .ne. 0) then
            do iy = 1, neldiff
                do ix = iy + 1, neldiff
                    write (ufort10, *) - 1, ix, iy
                end do
            end do
        end if
        write (ufort10, *) ' #          Nonzero values of  jasmat'
        do i1 = 1, njastot
            do i2 = 1, njastot
                if (jasyes(i1, i2)) then
                    if (yesalloc_jas) then
                        write (ufort10, *) i1, i2, jasmat(i1, i2)
                    else
                        write (ufort10, *) i1, i2, 0.d0
                    end if
                end if
            end do
        end do

        if (tf .eq. -9 .or. tf .eq. -8 .or. tf .eq. -19 .or. tf .eq. -18 .or. tf .eq. -16 &
            .or. tf .eq. -28 .or. tf .eq. -29) then
            write (*, *) " Spin Jastrow Added. "
            write (ufort10, *) "  #       Nonzero values of  jasmat Sz"
            do i1 = 1, njastot
                do i2 = 1, njastot
                    if (jasyes(i1, i2)) write (ufort10, *) i1, i2, 0.d0
                end do
            end do
        end if

        write (ufort10, *) '#  Eq. par.in the 3-body Jastrow in the chosen basis'
        do i1 = 1, iesfree
            write (ufort10, *) lenrecj(i1), recordsymj(1:2, 1:lenrecj(i1), i1)
        end do

        write (ufort10, *) '# Eq. par. in the atomic Det par.in the chosen basis'
        do i1 = 1, ndetpar
            if (eqdet(i1)%neq .ne. 0) &
                    & write (ufort10, gformat) eqdet(i1)%neq, (eqdet(i1)%idx(i2), i2=1, iabs(eqdet(i1)%neq))
        end do
        write (ufort10, *) '# Eq. par. in the atomic 3-body  par. in the chosen basis'
        do i1 = 1, njaspar
            if (eqjas(i1)%neq .ne. 0) &
                    & write (ufort10, gformat) eqjas(i1)%neq, (eqjas(i1)%idx(i2), i2=1, iabs(eqjas(i1)%neq))
        end do
        close (ufort10)
        deallocate (yeskion)
    end subroutine write_fort10

    subroutine fill_jasmat
        implicit none
        double precision :: r1
        integer :: i1, i3, i4
        integer i_atom, i_mol, itotup, icostup
        real*8 fatjas
        i_mol = ipj*ncell*njasorb
        jasmat = 0.d0
        do i1 = 1, ncell
            do i3 = 1, njasorb
                do i4 = i3, njasorb
                    jasmat(i3 + (i1 - 1)*njasorb, i4 + (i1 - 1)*njasorb) = atomic_jasmat(i3, i4)
                end do
            end do
        end do
        if (ipj .eq. 2) then
            itotup = ncell*njasorb
            if (.not. noonebody) itotup = itotup + 1
! down-down
            do i1 = 1, ncell
                do i3 = 1, njasorb
                    do i4 = i3, njasorb
                        jasmat(i3 + (i1 - 1)*njasorb + itotup, i4 + (i1 - 1)*njasorb + itotup) = atomic_jasmat(i3, i4)
                    end do
                end do
            end do
! up-down
            do i1 = 1, ncell
                do i3 = 1, njasorb
                    do i4 = 1, njasorb
                        jasmat(i3 + (i1 - 1)*njasorb, i4 + (i1 - 1)*njasorb + itotup) = atomic_jasmat(i3, i4)
                    end do
                end do
            end do
        end if
        if (.not. noonebody) then
            i_atom = njasorb + 1
            i_mol = ipj*(ncell*njasorb + 1)
            if (ipj .eq. 2) then
!  coupling with  the constant orbitals corresponding to down electrons
                fatjas = dble(nel - 1)/dble(nel - nelup)
            else
                fatjas = 1.d0
            end if
            do i1 = 1, ncell
                do i3 = 1, njasorb
                    jasmat(i3 + (i1 - 1)*njasorb, i_mol) = atomic_jasmat(i3, i_atom)*fatjas
                end do
            end do
            if (ipj .eq. 2) then
!  coupling with  the constant orbitals corresponding to up electrons
                fatjas = dble(nel - 1)/dble(nelup)
                icostup = ncell*njasorb + 1
                do i1 = 1, ncell
                    do i3 = 1, njasorb
                        jasmat(i3 + (i1 - 1)*njasorb + itotup, icostup) = atomic_jasmat(i3, i_atom)*fatjas
                    end do
                end do
            end if
        end if
!        symmetrize just for generality
        do i1 = 1, i_mol
            do i3 = i1 + 1, i_mol
                jasmat(i3, i1) = jasmat(i1, i3)
            end do
        end do
    end subroutine fill_jasmat

    subroutine fill_detmat
        implicit none
        double precision :: r1
        integer :: i1, i3, i4

        detmat = 0.d0
        if (ipf .eq. 2) then
            do i1 = 1, ndettot/2
                detmat(i1, ndettot/2 + i1) = 1.d0
                detmat(i1 + ndettot/2, i1) = 1.d0
            end do
        elseif (filling == "diagonal") then
            do i1 = 1, ndettot
                detmat(i1, i1) = 1.d0
            end do
        elseif (filling == "random") then
            call random_seed
            do i1 = 1, ndettot
                do i2 = i1, ndettot
                    call random_number(r1)
                    detmat(i1, i2) = r1
                end do
            end do
        elseif (filling == "semidiagonal") then
            call random_seed
            do i1 = 1, ndettot
                do i2 = i1 + 1, ndettot
                    call random_number(r1)
                    detmat(i1, i2) = r1*0.1
                end do
            end do
            do i1 = 1, ndettot
                detmat(i1, i1) = 1.d0
            end do
        elseif (filling == "atomic") then
            do i1 = 1, ncell
                do i3 = 1, ndetorb
                    do i4 = i3, ndetorb
                        detmat(i3 + (i1 - 1)*ndetorb, i4 + (i1 - 1)*ndetorb) = atomic_detmat(i3, i4)
                    end do
                end do
            end do
        end if
    end subroutine fill_detmat

    ! suboutine to fill detmat for complex wave function
    ! detmat always remains real in main subroutines
    subroutine fill_detmat_complex
        implicit none
        double precision :: r1, r2
        integer :: i1, i2, i3, i4

        detmat = 0.d0

        if (ipf .eq. 2) then
            do i1 = 1, ndettot - 1
                detmat(2*i1 - 1, i1 + 1) = 1.d0
                detmat(2*i1 + 1, i1) = -1.d0
            end do
        elseif (filling == "diagonal") then
            do i1 = 1, 2*ndettot, 2
                detmat(i1, i1/2 + 1) = 1.d0
                detmat(i1 + 1, i1/2 + 1) = 0.d0
            end do
        elseif (filling == "random") then
            call random_seed
            do i1 = 1, 2*ndettot, 2
                do i2 = (i1 - i1/2), ndettot
                    call random_number(r1)
                    detmat(i1, i2) = r1 ! real part
                    detmat(i1 + 1, i2) = 0 ! imaginary part
                end do
            end do
        elseif (filling == "semidiagonal") then
            call random_seed
            do i1 = 1, 2*ndettot, 2
                do i2 = i1 - (i1/2) + 1, ndettot
                    call random_number(r1)
                    detmat(i1, i2) = r1*0.1 ! real part
                    detmat(i1 + 1, i2) = 0 ! imaginary part
                end do
            end do
            do i1 = 1, 2*ndettot, 2
                detmat(i1, i1/2 + 1) = 1.d0
                detmat(i1 + 1, i1/2 + 1) = 0
            end do
            !
            ! atomic detmat not implemented yet for complex wave function
            !elseif(filling=="atomic") then
            !   do i1=1,ncell
            !      do i2=i1,ncell
            !         if(i1==i2) then
            !            do i3=1,ndetorb
            !               do i4=i3,ndetorb
            !                  detmat(i3+(i1-1)*ndetorb,i4+(i2-1)*ndetorb)=atomic_detmat(i3,i4)
            !               enddo
            !            enddo
            !         endif
            !      enddo
            !   enddo
            !
        end if
    end subroutine fill_detmat_complex

    subroutine check_symm(vecpbc, rot_det, yes_hermite)
        implicit none
        integer, intent(in) :: vecpbc(3)
        integer :: vecmap(3), trace, i
        logical ok_pbc, ok_try, yes_hermite
        logical, intent(in) :: rot_det
        logical allreal
        real*8 vecmul(3), phased(3), one3(3)
        nrot = 48
        nsym = 0
        one3 = 1.d0
        phased = phase
        call makeimage(phased, one3, deps)
        allreal = .true.
        do is = 1, nrot
            ok_pbc = .true.

            ! Exclude rotation upon request.
            if (.not. rot_det) then
                trace = 0
                do i1 = 1, 3
                    trace = trace + abs(isymm(i1, i1, is))
                end do
                if (trace .ne. 3) ok_pbc = .false.
            end if

            found = .true.
            do i1 = 1, ntotatoms
                if (atoms_maprot(i1, is) .eq. 0) found = .false.
            end do
            if (found .and. ok_pbc) then
                nsym = nsym + 1
                write (6, *) ' Rotation symmetry accepted =', is
                if (nsym .ne. is) then
                    atoms_maprot(1:ntotatoms, nsym) = atoms_maprot(1:ntotatoms, is)
                    isymm(1:3, 1:3, nsym) = isymm(1:3, 1:3, is)
                end if
            end if
        end do

        write (6, *) ' Rotation symmetries found =', nsym

        ntra = 0
        do is = 1, ntraind
            found = .true.
            do i1 = 1, ntotatoms
                if (atoms_maptra(i1, is) .eq. 0) found = .false.
            end do
            if (found) then
                ntra = ntra + 1
                if (ntra .ne. is) atoms_maptra(1:ntotatoms, ntra) = atoms_maptra(1:ntotatoms, is)
            end if
        end do

        write (6, *) ' Translation symmetries found =', ntra, ntraind
        if (ntra .ne. ntraind) then
            write (6, *) ' ERROR you should check your primitive cell '
            stop
        end if

    end subroutine check_symm

    subroutine read_input(funit)
        implicit none
        integer :: i1, i2, i3, i4, i5, nskip, ind, ref_atom
        integer :: c_onebodypar, count_non_zero_onebodypar ! KN
        integer, intent(in) :: funit
        real(8) :: volume, scalfact, vol_ratio

        scalfact = 1.d0
        volume = 1.d0
        vol_ratio = 1.d0

        nel = -1
        neldiff = -1
        numpaired = 0
        rs_read = -1.d0
        L_read = -1.d0
        unit_crystal = 'default'
        celldm(1) = -1.d0
        celldm(2) = 1.d0
        celldm(3) = 1.d0
        celldm(4:6) = Pi/2.d0

        ntraind = 1

        at = 0.0

        axyz = 0
        do i = 1, 3
            axyz(i, i) = 1.d0
        end do
        phase(:) = 0.d0
        phasedo(:) = -100.d0
        nxyz(:) = 1
        posunits = "bohr"
        pbcfort10 = .true.
        complexfort10 = .false. ! flag for creating a complex wave function
        real_contracted = .false.
        ! only detmat and contracted coefficients are modified

        ntyp = 0 ! not used if fort.10 are not read
        write_log = .false.

        twobody = -6
        !twobodypar = 1. !KN
        !onebodypar = 0. !KN
        filling = "undefined"
        noonebody = .false.
        readatoms = .false.
        orbtype = "normal"
        jorbtype = "normal"
        niesd = -1
        n_onebody = -1 !KN
        n_twobody = -1 !KN
        c_onebodypar = 0 !KN
        count_non_zero_onebodypar = 0 !KN
        symmagp = .true.
        onlycontrdet = .false.
        onlycontrjas = .false.
        yes_pfaff = .false.
        nodownpfaff = .false.
        nouppfaff = .false.
        eqatoms = .true.
        eq_intatoms = .true.
        nosym = .false.
        nosym_contr = .false.
        nosym_contrj = .false.
        notra = .false.
        rot_det = .true.
        rot_jas = .true.
        rot_pfaff = .false.
        readunpaired = .false.

        nshelljas = -1
        nshelldet = 0

        vecpbc(1:3) = -1

        yes_pfaff = .false.
        yes_tilted = .false.

        read (funit, nml=system, err=100, end=100)

        if (.not. pbcfort10 .and. celldm(1) .eq. -1) then
            write (6, *) ' Warning celldm(1) set to 50 ! '
            celldm(1) = 50.d0
        end if
        if (yes_pfaff) then
            ipf = 2
        else
            ipf = 1
        end if
        nodownpfaff_read = nodownpfaff

        if (write_log) &
                &   open (file="makefort10.log", status="unknown", form="formatted", unit=lunit)

        if (celldm(1) /= -1.d0 .and. sum(abs(at(:, :))) /= 0.d0) then
            call errore("read_input", " Use celldm(:) or at(:,:) but not both!!", 1)
        end if

        if (unit_crystal == "default") then
!           if(yes_tilted.and.product(nxyz(:)).eq.1) then
!               unit_crystal = 'primitive'
!               write(6, *) ' Unit L_read = primitive (at) '
!           else
            unit_crystal = 'conventional'
            write (6, *) ' Unit L_read = conventional (axyz X at) '
!           endif
        end if

        if (posunits == "angstrom") then
            if (celldm(1) .ne. -1.d0) celldm(1) = celldm(1)/length_unit
            if (rs_read .ne. -1.d0) rs_read = rs_read/length_unit
            if (L_read .ne. -1.d0) L_read = L_read/length_unit
        end if

        if ((nxyz(1) .ne. 1 .or. nxyz(2) .ne. 1 .or. nxyz(3) .ne. 1 &
             .or. (sum(abs(at(:, :))) .ne. 0.d0) .and. celldm(1) .eq. -1.d0 .and. pbcfort10)) then
            yes_trivial = .false.
        else
            yes_trivial = .true.
        end if

        complexfort10_sav = complexfort10

        if (.not. pbcfort10) then
            phase = 0.d0
            phasedo = 0.d0
        end if

        do i1 = 1, 3
            if (phasedo(i1) .eq. -100.d0) then
                if (phase(i1) - nint(phase(i1)) .ne. 0 .and. abs(phase(i1) - nint(phase(i1))) .ne. 0.5) then
                    phasedo(i1) = -phase(i1)
                else
                    phasedo(i1) = phase(i1)
                end if
            end if
        end do
        do i1 = 1, 3
            if (abs(phase(i1) - nint(phase(i1))) .ne. 0.d0 .and. abs(phase(i1) - nint(phase(i1))) .ne. 0.5d0) then
                complexfort10 = .true.
            end if
            if (abs(phasedo(i1) - nint(phasedo(i1))) .ne. 0.d0 .and. abs(phasedo(i1) - nint(phasedo(i1))) .ne. 0.5d0) then
                complexfort10 = .true.
            end if

            !  Real case with different boundary conditions not implemented
            if (phasedo(i1) - nint(phasedo(i1)) .ne. phase(i1) - nint(phase(i1)) .and. .not. &
                (abs(phasedo(i1) - nint(phasedo(i1))) .eq. 0.5 .and. abs(phase(i1) - nint(phase(i1))) .eq. 0.5)) then
                complexfort10 = .true.
            end if
        end do

        if (sum(abs(at(:, :))) /= 0.d0 .and. sum(abs(axyz(:, :))) /= 0.d0) then

            !
            newcell = 0.d0
            no_orthorombic = .true.
            !
            newcell = matmul(at, axyz)
            !
            ! Orthorombic Check
            ! Check off-diagonal elements
            !
            is_ok = .true.
            !
            write (6, *) ' newcell ='

            do i1 = 1, 3
                write (6, *) i1, newcell(:, i1)
                do i2 = i1 + 1, 3
                    if (abs(newcell(i1, i2)) .gt. 1.e-7 .or. abs(newcell(i2, i1)) .gt. 1e-7) is_ok = .false.
                end do
            end do
            !
            if (.not. is_ok .and. .not. yes_tilted) &
                    &   call errore("read_input", " Only Orthorombic cells supported! check axyz(:,:) ! ", 1)
            !
            celldm(1) = dsqrt(sum(newcell(:, 1)**2))
            celldm(2) = dsqrt(sum(newcell(:, 2)**2))/celldm(1)
            celldm(3) = dsqrt(sum(newcell(:, 3)**2))/celldm(1)
            write (*, '(a,3f12.6)') ' New orthorombic cell : ', celldm(1), celldm(2:3)*celldm(1)
            smallcell(1) = celldm(1)
            smallcell(2) = celldm(2)*celldm(1)
            smallcell(3) = celldm(3)*celldm(1)
            if (yes_tilted) then
                !   alphap angle between b and c
                celldm(4) = dacos(sum(newcell(:, 2)*newcell(:, 3))/smallcell(2)/smallcell(3))
                !   betap = angle from c and a
                celldm(5) = dacos(sum(newcell(:, 3)*newcell(:, 1))/smallcell(3)/smallcell(1))
                !   gammap = angle between  a  and b
                celldm(6) = dacos(sum(newcell(:, 1)*newcell(:, 2))/smallcell(1)/smallcell(2))

            else
                celldm(4:6) = Pi/2.d0
            end if
            !
            !
        elseif (sum(abs(at(:, :))) /= 0.d0) then
            newcell = at
            no_orthorombic = .false.
            celldm(1) = dsqrt(sum(newcell(:, 1)**2))
            celldm(2) = dsqrt(sum(newcell(:, 2)**2))/celldm(1)
            celldm(3) = dsqrt(sum(newcell(:, 3)**2))/celldm(1)
            smallcell(1) = celldm(1)
            smallcell(2) = celldm(2)*celldm(1)
            smallcell(3) = celldm(3)*celldm(1)
            write (*, '(a,3f12.6)') ' New cell : ', celldm(1), celldm(2:3)*celldm(1)
            if (yes_tilted) then
                !   alphap angle between b and c
                celldm(4) = dacos(sum(newcell(:, 2)*newcell(:, 3))/smallcell(2)/smallcell(3))
                !   betap = angle from c and a
                celldm(5) = dacos(sum(newcell(:, 3)*newcell(:, 1))/smallcell(3)/smallcell(1))
                !   gammap = angle between  a  and b
                celldm(6) = dacos(sum(newcell(:, 1)*newcell(:, 2))/smallcell(1)/smallcell(2))

            else
                celldm(4:6) = Pi/2.d0
            end if
        elseif (celldm(1) .gt. 0.d0) then

            smallcell(1) = celldm(1)
            smallcell(2) = celldm(2)*celldm(1)
            smallcell(3) = celldm(3)*celldm(1)
            write (*, '(a,3f12.6)') ' New cell : ', celldm(1), celldm(2:3)*celldm(1)
            if (yes_tilted) then
                !   alphap angle between b and c
                alphap = celldm(4)
                betap = celldm(5)
                gammap = celldm(6)

                omega = product(smallcell(:))*dsqrt(abs(&
                        & 1.d0 - dcos(alphap)**2.d0 - dcos(betap)**2.d0 - dcos(gammap)**2.d0 + &
                                & 2.d0*dcos(alphap)*dcos(betap)*dcos(gammap)))

                newcell(1, 1) = smallcell(1)
                newcell(2, 1) = 0.d0
                newcell(3, 1) = 0.d0
                newcell(1, 2) = smallcell(2)*dcos(gammap)
                newcell(2, 2) = smallcell(2)*dsin(gammap)
                newcell(3, 2) = 0.d0
                newcell(1, 3) = smallcell(3)*dcos(betap)
                newcell(2, 3) = smallcell(3)*(dcos(alphap) - dcos(betap)*dcos(gammap))       &
                        & /dsin(gammap)
                newcell(3, 3) = omega/(smallcell(1)*smallcell(2)*dsin(gammap))

            else
                newcell = 0.d0
                do i = 1, 3
                    newcell(i, i) = smallcell(i)
                end do
            end if
            at = newcell
            no_orthorombic = .false.
        else

            no_orthorombic = .false.

        end if

        shiftbeta = 1 ! default value for shiftbeta

        yesbump = .false. ! default value for yesbump
        yes_crystal = .true. ! default value for PBC
        yes_crystalj = .false. ! false in any case
        no_4body_jas = .false.
        nopseudo = .false.

        !! added by K.Nakano on 18 Feb. for generalizing onebody and twobody params.
        allocate (twobodypar(2))
        allocate (onebodypar(ntyp))
        twobodypar = 1.0d0 ! 1.0d0 for all element
        onebodypar = 0.0d0 ! 0.0d0 for all element
        scale_jasfat = 0.d0 ! default do not use scale_jasfat

        read (funit, nml=electrons, err=101, end=101)

        ! count non zero onebodypar added by KN on 18.Feb.
        count_non_zero_onebodypar = 0
        do c_onebodypar = 1, size(onebodypar)
            if (onebodypar(c_onebodypar) .ne. 0.0d0) then
                count_non_zero_onebodypar = count_non_zero_onebodypar + 1
            end if
        end do

        ! set n_onebody added by KN on 18.Feb.
        ! n_twobody is defined in the main routine.
        if (n_onebody .eq. -1) then
            if (twobody .eq. 0) then
                n_onebody = 0
            else
                if (count_non_zero_onebodypar .eq. 0) then
                    n_onebody = 0
                elseif (count_non_zero_onebodypar .eq. 1) then
                    n_onebody = 1
                else
                    n_onebody = ntyp
                end if
            end if
        end if

        write (6, *) ' n_onebody is set', n_onebody

        ! commented out by K.Nakano on 18 Feb.
        !if(niesd.eq.-1) then
        !    if(twobody.eq.0) then
        !        niesd = 0
        !    elseif(onebodypar.eq.0) then
        !        niesd = 1
        !    else
        !        niesd = 2
        !    endif
        !endif

        if (complexfort10 .and. pbcfort10 .and. .not. yes_crystal) then
            write (6, *) ' Warning complex and pbc only with yes_crystal=.true. basis, &
                    & Changing to yes_crystal=.true. !!! '
            yes_crystal = .true.
        end if
        if (.not. yes_crystal .and. yes_crystalj) then
            write (6, *) ' Warning yes_crystalj=.true. possible only when yes_crystal=.true. yes_crystalj changed to false '
            yes_crystalj = .false.
        end if

        if (ntyp .eq. 0) call errore("read_input", " Warning  ntyp == 0 ! ", -1)

        if (rs_read .eq. -1.d0 .and. L_read .eq. -1.d0) then
            !
            smallcell(1) = celldm(1)
            smallcell(2) = celldm(2)*celldm(1)
            smallcell(3) = celldm(3)*celldm(1)

            ! volume_new=abs(newcell(1,1)*(newcell(2,2)*newcell(3,3)-newcell(2,3)*newcell(3,2))&
            ! -newcell(1,2)*(newcell(2,1)*newcell(3,3)-newcell(2,3)*newcell(3,1))&
            ! +newcell(1,3)*(newcell(2,1)*newcell(3,2)-newcell(2,2)*newcell(3,1)))

            !        unit_volume=volume_new/sqrt(sum(newcell(:,1)**2))**3

            !
        else
            !
            !  if(no_orthorombic) call errore("read_input",'rs_try and not orthorombic still not implemented!',1)

            if (no_orthorombic) then

                volume = abs(at(1, 1)*(at(2, 2)*at(3, 3) - at(2, 3)*at(3, 2)) - at(1, 2)*(at(2, 1)*at(3, 3) - &
                        &at(2, 3)*at(3, 1)) + at(1, 3)*(at(2, 1)*at(3, 2) - at(2, 2)*at(3, 1)))

                volume_new = abs(newcell(1, 1)*(newcell(2, 2)*newcell(3, 3) - newcell(2, 3)*newcell(3, 2)) &
                                 - newcell(1, 2)*(newcell(2, 1)*newcell(3, 3) - newcell(2, 3)*newcell(3, 1)) &
                                 + newcell(1, 3)*(newcell(2, 1)*newcell(3, 2) - newcell(2, 2)*newcell(3, 1)))

                !        unit_volume=volume_new/sqrt(sum(newcell(:,1)**2))**3
                !         volume=abs(newcell(1,1)*newcell(2,2)*newcell(3,3))/volume

                vol_ratio = volume_new/volume

                if (yes_tilted) then

                    write (6, *) ' Ratio unit supercell volume / unit cell volume =', vol_ratio

                else

                    write (6, *) ' Ratio ortho cell volume / unit cell volume =', vol_ratio

                end if

                smallcell(1) = celldm(1)
                smallcell(2) = celldm(2)*celldm(1)
                smallcell(3) = celldm(3)*celldm(1)

            else

                !
                !

                smallcell(1) = celldm(1)
                smallcell(2) = celldm(2)*celldm(1)
                smallcell(3) = celldm(3)*celldm(1)

            end if
            !
        end if

        cellscale(:) = smallcell(:)

        if (write_log) then
            write (lunit, *) ' Primitive Cell '
            write (lunit, *) smallcell(:)
        end if

        volume_new = abs(newcell(1, 1)*(newcell(2, 2)*newcell(3, 3) - newcell(2, 3)*newcell(3, 2)) &
                         - newcell(1, 2)*(newcell(2, 1)*newcell(3, 3) - newcell(2, 3)*newcell(3, 1)) &
                         + newcell(1, 3)*(newcell(2, 1)*newcell(3, 2) - newcell(2, 2)*newcell(3, 1)))

        unit_volume = volume_new/sqrt(sum(newcell(:, 1)**2))**3

        forces_sym = .false.
        notra_forces = .false.
        nosym_forces = .false.
        forcesymm = .false.

        read (funit, nml=symmetries, err=102, end=102)

        !   If different phases we do not have the code for the real
        !   also if symmagp=.true. in input should be changed to false
        yes_hermite = .true.

        if (symmagp .and. yes_pfaff) then
            if (.not. nouppfaff .and. .not. nodownpfaff) then
                write (6, *) ' Warning set nodownpfaff true '
                nodownpfaff = .true.
                !      elseif(nouppfaff.and.nodownpfaff) then
                !      write(6,*) ' Warning set nouppfaff false '
                !      nouppfaff=.false.
            end if
        end if
        if (forcesymm) then
            opposite_phase = .false.
            phasedo(:) = phase(:)
        else
            opposite_phase = .true.
        end if
        do i1 = 1, 3
            if ((abs(phase(i1) - nint(phase(i1))) .ne. 0.5 .and.&
                    &phasedo(i1) - nint(phasedo(i1)) .ne. nint(phase(i1)) - phase(i1)) .or.&
                    &(abs(phase(i1) - nint(phase(i1))) .eq. 0.5 .and. abs(phasedo(i1) - nint(phasedo(i1))) .ne. .5)) then
                opposite_phase = .false.
            end if
        end do
        if (.not. opposite_phase) then
            do i1 = 1, 3
                if (phasedo(i1) - nint(phasedo(i1)) .ne. phase(i1) - nint(phase(i1))) then
                    if (symmagp) then
                        symmagp = .false.
                        write (6, *) ' Warning different phase up =/down, symmagp forced to false '
                    end if
                end if
            end do
        end if

        write (6, *) ' opposite_phase =', opposite_phase

        if (.not. yes_crystal .and. yes_hermite) yes_hermite = .false.
        !
        ! Build a new orthorombic cell from a non-orthorombic one
        if ((.not. yes_crystal .or. (sum(abs(phase(:) - nint(phase(:)))) .eq. 0.d0 .and. .not. complexfort10)) &
            .and. yes_hermite) then
            yes_hermite = .false.
        end if
        !
        if (yes_hermite) write (6, *) ' Warning boundary conditions with flux attaching: different algorithm '

        if (.not. forces_sym) then
            nosym_forces = nosym
            notra_forces = notra
        end if

        if (.not. no_orthorombic) then
            ncell = product(nxyz)
            ntotatoms = ncell*natoms
        else
            ntotatoms = natoms
            ncell = 1
        end if

        allocate (zeta(2, ntotatoms), rion(3, ntotatoms), J3_off(ntotatoms)&
                &, cut_hybrid(ntotatoms), cut_hybridj(ntotatoms))
        cut_hybrid = -1 ! No cut by default

        cut_hybridj = -1 ! No cut in Jastrow

        call findsection(funit, "ATOMIC_POSITIONS")
        do i1 = 1, natoms
            ! same notations as in fort.10 numbers are more general than letters
            read (funit, *, err=103, end=103) zeta(2, i1), zeta(1, i1), rion(:, i1)
            write (6, *) ' read zeta =', zeta(1, i1)
        end do

        nel_read = sum(zeta(2, 1:natoms))

        if (celldm(1) .eq. -1 .and. .not. no_orthorombic) then
            if (rs_read .ne. -1.d0) then
                volume = 4.d0/3.d0*pi*rs_read**3*nel_read
                !          celldm(1)=volume**(1.d0/3.d0)
                smallcell(1:3) = volume**(1.d0/3.d0)

                newcell = 0.d0
                do i = 1, 3
                    newcell(i, i) = smallcell(i)
                end do

            elseif (L_read .eq. -1.d0) then
                !          celldm(1)=L_read
                volume = L_read**3
                smallcell(1:3) = L_read
                newcell = 0.d0
                do i = 1, 3
                    newcell(i, i) = smallcell(i)
                end do
            end if
        end if

        write (6, *) ' nel_read found =', nel_read

        if (posunits == "crystal" .and. .not. pbcfort10) &
                & call errore("read_input", 'You cannot use "crystal" units with pbcfort10=.false.', 1)

        if (posunits == "angstrom") then
            !
            rion(:, :) = rion(:, :)/length_unit
            !
        elseif (posunits == "crystal") then

            do i1 = 1, natoms
                rion(:, i1) = matmul(at, rion(:, i1))
            end do

            !
        elseif (posunits == "bohr") then
            !
        else
            !
            call errore("read_input", "Unknown units", 1)
            !
        end if

        if (no_orthorombic) then ! creating a new orthorombic cell
            !
            max_factor(1) = maxval(abs(axyz(:, 1))) + 1
            max_factor(2) = maxval(abs(axyz(:, 2))) + 1
            max_factor(3) = maxval(abs(axyz(:, 3))) + 1
            !
            allocate (new_rion(3, product(max_factor(:)*natoms)))
            allocate (tmp_zeta(2, product(max_factor(:)*natoms)))
            !
            max_ions = 0
            do i1 = 1, 3
                at_qe(:, i1) = newcell(:, i1)/dsqrt(sum(newcell(:, i1)**2))
            end do
            call init_cell

            !
            ! Find the equivalent atoms sequentially so that the independent translations
            ! are easy to find.
            do i4 = 1, natoms
                do i1 = -max_factor(1), max_factor(1)
                    do i2 = -max_factor(2), max_factor(2)
                        do i3 = -max_factor(3), max_factor(3)
                            test_rion(:) = rion(:, i4) + at(:, 1)*dble(i1) + at(:, 2)*dble(i2) + at(:, 3)*dble(i3)
                            test_rion_before = test_rion

                            call CartesianToCrystal(test_rion, 1)
                            if (test_rion(1) >= -smallcell(1)/2.d0 - deps .and. test_rion(1) < smallcell(1)/2.d0 - deps .and. &
                                test_rion(2) >= -smallcell(2)/2.d0 - deps .and. test_rion(2) < smallcell(2)/2.d0 - deps .and. &
                                test_rion(3) >= -smallcell(3)/2.d0 - deps .and. test_rion(3) < smallcell(3)/2.d0 - deps) then
                                !
                                found = .false.

                                !
                                i5 = 1
                                do while (.not. found .and. i5 .lt. max_ions)
                                    if (sum(abs(test_rion_before(:) - new_rion(:, i5))) .lt. deps&
                                            &.and. abs(int(zeta(1, i4)) - int(tmp_zeta(1, i5))) .lt. deps) found = .true.

                                    i5 = i5 + 1
                                end do
                                !
                                if (.not. found) then
                                    max_ions = max_ions + 1
                                    new_rion(:, max_ions) = test_rion_before(:)
                                    tmp_zeta(:, max_ions) = zeta(:, i4)
                                end if
                                !
                            end if
                        end do
                    end do
                end do
            end do
            !
            deallocate (rion, zeta, J3_off)
            !
            ntraind = max_ions/natoms
            if (ntraind*natoms .ne. max_ions) then
                write (6, *) ' There should be some error !!! ', max_ions, ntraind, natoms
                stop
            end if

            ! reading Wigner-Seize radius
            if (rs_read .ne. -1.d0) then
                if (celldm(1) .ne. -1.d0) then
                    rs_try = (3.d0/4.d0*volume_new/pi/nel_read)**(1.d0/3d0)
                    celldm(1) = celldm(1)*rs_read/rs_try
                else
                    celldm(1) = smallcell(1)*rs_read/rs_try
                end if
                !
                scalfact = celldm(1)/smallcell(1)
                smallcell(:) = smallcell(:)*scalfact
                new_rion(1:3, 1:max_ions) = new_rion(1:3, 1:max_ions)*scalfact
                ! reading cell parameter
            elseif (L_read .ne. -1.d0) then

                if (trim(unit_crystal) == "primitive") then
                    !        L_read refers to the primitive cell
                    scalfact = L_read/sqrt(sum(at(:, 1)**2))
                    celldm(1) = smallcell(1)*scalfact
                    smallcell(:) = smallcell(:)*scalfact
                    new_rion(1:3, 1:max_ions) = new_rion(1:3, 1:max_ions)*scalfact
                    newcell = newcell*scalfact
                    at = at*scalfact
                else
                    celldm(1) = L_read
                    scalfact = L_read/smallcell(1)
                    smallcell(:) = smallcell(:)*scalfact
                    new_rion(1:3, 1:max_ions) = new_rion(1:3, 1:max_ions)*scalfact
                    newcell = newcell*scalfact
                    at = at*scalfact
                end if
            end if

            natoms = max_ions
            ncell = product(nxyz)
            ntotatoms = ncell*natoms
            !
            write (*, *) 'Number of atoms in the orthorombic cell : ', max_ions
            write (*, *) 'Total number of atoms : ', ntotatoms
            !
            allocate (zeta(2, ntotatoms), rion(3, ntotatoms)&
                    &, J3_off(ntotatoms))
            if (allocated(cut_hybrid)) deallocate (cut_hybrid, cut_hybridj)
            allocate (cut_hybrid(ntotatoms), cut_hybridj(ntotatoms))
            cut_hybridj = -1 ! no cut on Jastrow
            cut_hybrid = -1 ! no cut by default
            !
            zeta(1:2, 1:natoms) = tmp_zeta(1:2, 1:natoms)
            rion(1:3, 1:natoms) = new_rion(1:3, 1:natoms)
            !
            !

        else ! cell is already orthorombic
            write (6, *) ' HERE L_read, at =', rs_read, L_read, at(:, 1)

            ! reading Wigner-Seize radius
            if (rs_read .ne. -1.d0) then

                if (celldm(1) .ne. -1.d0) then
                    rs_try = (3.d0/4.d0*volume_new/pi/nel_read)**(1.d0/3d0)

                    celldm(1) = celldm(1)*rs_read/rs_try
                else
                    celldm(1) = smallcell(1)
                end if

                !
                scalfact = celldm(1)/smallcell(1)

                smallcell(:) = smallcell(:)*scalfact
!               new_rion(1:3, 1:max_ions) = new_rion(1:3, 1:max_ions) * scalfact
                ! reading cell parameter
            elseif (L_read .ne. -1.d0) then

                if (trim(unit_crystal) == "primitive") then
                    !        L_read refers to the primitive cell
                    scalfact = L_read/sqrt(sum(at(:, 1)**2))
                    celldm(1) = smallcell(1)*scalfact
                    smallcell(:) = smallcell(:)*scalfact
!                   new_rion(1:3, 1:max_ions) = new_rion(1:3, 1:max_ions) * scalfact
                    newcell = newcell*scalfact
                    at = at*scalfact
                else
                    celldm(1) = L_read
                    scalfact = L_read/smallcell(1)
                    smallcell(:) = smallcell(:)*scalfact
!                   new_rion(1:3, 1:max_ions) = new_rion(1:3, 1:max_ions) * scalfact
                    newcell = newcell*scalfact
                    at = at*scalfact
                end if
            end if

        end if ! end if(no_orthorombic)

        cellscale(1:3) = smallcell(1:3)*dble(nxyz(1:3))

        write (*, *) 'Atoms positions '
        do i1 = 1, natoms
            write (*, '(1f6.0,1f7.1,3f12.6)') zeta(:, i1), rion(:, i1)
        end do

        write (*, *) 'Unit cell volume (a_0^3) =', smallcell(1)**3*unit_volume
        write (*, *) 'Unit cell volume (A^3)   =', smallcell(1)**3*unit_volume*length_unit**3
        write (*, *) 'Supercell volume (a_0^3) =', smallcell(1)**3*product(nxyz(:))*unit_volume
        write (*, *) 'Supercell volume (A^3)   =', smallcell(1)**3*product(nxyz(:))*unit_volume*length_unit**3

        !    redefine newcell
        do j = 1, 3
            newcell(:, j) = newcell(:, j)*smallcell(j)/sqrt(sum(newcell(:, j)**2))
            write (6, *) j, newcell(:, j)
        end do
        ! ********** Find the type of each atom *********

        volume_new = abs(newcell(1, 1)*(newcell(2, 2)*newcell(3, 3) - newcell(2, 3)*newcell(3, 2)) &
                         - newcell(1, 2)*(newcell(2, 1)*newcell(3, 3) - newcell(2, 3)*newcell(3, 1)) &
                         + newcell(1, 3)*(newcell(2, 1)*newcell(3, 2) - newcell(2, 2)*newcell(3, 1)))

        if (readatoms) then
            call findsection(funit, "ATOMIC_SPECIES")
            allocate (atypes(ntyp))
            do i1 = 1, ntyp
                read (funit, *, err=105, end=105) atypes(i1)%atomic_number, atypes(i1)%wf_name
            end do
        end if

        if (readatoms) then
            allocate (mytype(natoms))
            do i1 = 1, natoms
                found = .false.
                i2 = 1
                do while (i2 <= ntyp .and. .not. found)
                    if (int(zeta(1, i1)) == atypes(i2)%atomic_number) then
                        found = .true.
                        mytype(i1) = i2
                    end if
                    i2 = i2 + 1
                end do
                if (.not. found) call errore("read_input", " Atom type not found! ", 1)
            end do
        end if

        if (nel == -1) nel = sum(zeta(2, 1:natoms))
        ! neldiff is a switch for putting the unpaired orbitals in the det.
        ! The case for Cr, 10 nelup and 4 neldo

        ! For the odd number of electron system is handled by UNPAIRED section
        ! if(mod(nel*ncell,2).ne.0) neldiff=1

        ! ****** reading the basis set for DET and JAS ******

        if (readatoms) then ! reading wave function from external file
            !
            if (filling == "undefined") filling = "atomic"
            !
            call read_atoms(natoms, zeta, atomic_detmat, detorb, ndetorb, nshelldet, ndetpar &
           &, jasorb, njasorb, nshelljas, njaspar, mytype, atypes, ntyp)
            !
        else ! reading wave function from input

            call read_orbitals(natoms, nel, noonebody, zeta, J3_off, cut_hybrid &
                               , detorb, ndetorb, nshelldet, ndetpar, jasorb, njasorb, nshelljas &
                               , njaspar, orbtype, jorbtype, funit, shiftbeta, complexfort10 &
                               , symmagp, ncell, scale_jasfat)
            !
            ! Put one to the first coefficient if it is zero
            !
            !  if(detorb(1)%nparm.gt.1) &
            !&  detorb(1)%parms(detorb(1)%nparm/2+1)=1.d0
            !  do i1=2,ndetorb
            !    if(detorb(i1)%shell.ne.detorb(i1-1)%shell) then
            !      if(detorb(i1)%parms(detorb(i1)%nparm/2+1)==0.d0) &
            !&            detorb(i1)%parms(detorb(i1)%nparm/2+1)=1.d0
            !    endif
            !  enddo
            !
            if (nshelljas .ne. 0) then
                !
                !    if(jasorb(1)%nparm.gt.1) &
                ! &   jasorb(1)%parms(jasorb(1)%nparm/2+1)=1.d0
                !    do i1=2,njasorb
                !      if(jasorb(i1)%shell.ne.jasorb(i1-1)%shell) then
                !        if(jasorb(i1)%parms(jasorb(i1)%nparm/2+1)==0.d0) &
                !&         jasorb(i1)%parms(jasorb(i1)%nparm/2+1)=1.d0
                !      endif
                !    enddo
                !
            else
                !
                njaspar = 0
                njasorb = 0
                neqjas = 0
                !
            end if
            !
        end if
        ! ***** end reading basis set *****

        if (rs_read .ne. -1.d0) then
            !if(nel_read.ne.nel) then
            if ((.not. yes_trivial .and. nel_read .ne. nel)&
                    &.or. (yes_trivial .and. nel_read .ne. nel)) then
                write (6, *) ' Error in input  nel_read does not match nel '&
                        &, volume, nel_read, nel
                stop
            end if
        end if

        write (*, *) ' # Orbitals in the Determinant : ', ndetorb
        write (*, *) ' # Parameters in the Determinant : ', ndetpar
        write (*, *) ' # Orbitals in Jastrow : ', njasorb
        write (*, *) ' # Parameters in Jastrow : ', njaspar

        unpaired = .false.
        norb_unpaired = 0

        if (readunpaired) then
            !
            write (*, *) ' Odd Number of Electrons... reading UNPAIRED section '
            call findsection(funit, "UNPAIRED")
            unpaired = .true.
            read (funit, *, err=106, end=106) norb_unpaired
            allocate (orb_unpaired(2, norb_unpaired))
            read (funit, *, err=106, end=106) (orb_unpaired(:, i1), i1=1, norb_unpaired)
            allocate (orb_numbers(norb_unpaired))
            orb_numbers = 0
            do i1 = 1, norb_unpaired
                i2 = 1
                found = .false.
                do while (i2 .le. ndetorb .and. .not. found)
                    !
                    if (orb_unpaired(1, i1) == detorb(i2)%kion .and. orb_unpaired(2, i1) == detorb(i2)%norbatom) then
                        orb_numbers(i1) = i2
                        found = .true.
                    end if
                    i2 = i2 + 1
                    !
                end do
            end do
            !
            if (any(orb_numbers == 0)) &
                    & call errore("makefort10", " Error in UNPAIRED section check orbital numbers! ", 1)
        end if

        ! ********* Set APBC flag for each side ********
        apbc = .false.
        if (.not. complexfort10 .and. .not. yes_crystal) then
            do i1 = 1, 3
                if (abs(phase(i1) - nint(phase(i1))) .eq. 0.5d0) then
                    apbc(i1) = .true.
                else
                    apbc(i1) = .false.
                    if (.not. complexfort10 .and. abs(phase(i1) - nint(phase(i1))) .ne. 0.d0) then
                        call errore("makefort10", &
                                    " Symmetries with phase != 0 or 0.5 not implemented for real wave function! &
&                                    Use PBC_C instead", 1)
                    end if
                end if
            end do
        end if

        !    if(yes_crystal.or.complexfort10) then
        !    do i1=1,3
        !     if(apbc(i1)) then
        !     write(6,*) ' ERROR !!! in this case translation invariance assumed'
        !     write(6,*) ' Please check makefort10 input...'
        !     stop
        !     endif
        !    enddo
        !    endif
        do i1 = 2, ncell
            do i2 = 1, natoms
                zeta(:, (i1 - 1)*natoms + i2) = zeta(:, i2)
                J3_off((i1 - 1)*natoms + i2) = J3_off(i2)
                cut_hybrid((i1 - 1)*natoms + i2) = cut_hybrid(i2)
            end do
        end do

        counter = 0
        write (*, *) ' Total number of cells/further ind. transl.  = ', ncell, ntraind

        if (.not. notra) then
            nskip = ntraind
            ntraind = ntraind*ncell
        else
            ntraind = ncell
            nskip = 1
        end if

        write (6, *) ' Number of independent translation =', ntraind
        allocate (trasl(3, ntraind))
        ind = 0
        write (6, *) ' Independent translation ='
!       do i1 = 1, 3
!           write(6, *) i1, newcell(:, i1)
!       enddo
        do i1 = 0, nxyz(1) - 1
            do i2 = 0, nxyz(2) - 1
                do i3 = 0, nxyz(3) - 1
                    ref_atom = 1

                    do i4 = 1, natoms
                        counter = counter + 1
                        rion(:, counter) = -i1*newcell(:, 1) - i2*newcell(:, 2) - i3*newcell(:, 3) + rion(:, i4)
                        !               rion(1,counter)=-i1*smallcell(1)+rion(1,i4)
                        !               rion(2,counter)=-i2*smallcell(2)+rion(2,i4)
                        !               rion(3,counter)=-i3*smallcell(3)+rion(3,i4)
                        if (i4 .le. nskip) then
                            ind = ind + 1
                            trasl(:, ind) = rion(:, counter) - rion(:, ref_atom)
                            write (6, *) ind, trasl(1:3, ind)
                        end if
                    end do
                end do
            end do
        end do
        if (yes_pfaff) write (6, *) 'Pfaffian wave function'
        !write(6,*) ' final counter =',counter

        !ntraind=1

        !do i1=1,ntraind
        !trasl(:,i1)=trasl(:,2*i1-1)
        !enddo

        !stop
        !close(funit)

        return
100     call errore("read_input", "Error reading 'system' namelist! ", 1)
101     call errore("read_input", "Error reading 'electrons' namelist! ", 1)
102     call errore("read_input", "Error reading 'symmetries' namelist! ", 1)
103     call errore("read_input", "Error reading ATOMIC_POSITIONS! ", 1)
104     call errore("read_input", "Error reading ATOMIC_WF! ", 1)
105     call errore("read_input", "Error reading ATOMIC_SPECIES! ", 1)
106     call errore("read_input", "Error reading UNPAIRED! ", 1)
    end subroutine read_input

    subroutine update_symrot(nsym, symmat, lda, atoms_map, orb, norb, orb_map)
        implicit none
        integer, intent(in) :: norb, lda, nsym, atoms_map(ntotatoms, nsym)
        type(orbmap), intent(in) :: orb_map
        type(orbital), intent(in) :: orb(norb)
        integer, intent(out) :: symmat(lda, nsym)
        integer ordorb, ordorbtry
        symmat = 0
        !write(6,*)  ' mapping output '
        do is = 1, nsym
            do i1 = 1, norb
                do i2 = 1, ncell
                    ionpos = orb(i1)%kion + (i2 - 1)*natoms
                    i2map = atoms_map(ionpos, is)
                    signmap = 1
                    if (i2map .lt. 0) then
                        signmap = -1
                        i2map = -i2map
                    end if

                    ncellmap = (i2map - 1)/natoms + 1
                    i2map = i2map - (ncellmap - 1)*natoms
                    if (orb_map%comp(i1, is)) then
                        if (orb_map%iorb(i1, is) .lt. 0) then
                            signmap = -signmap
                            indmap = -orb_map%iorb(i1, is)
                        else
                            indmap = orb_map%iorb(i1, is)
                        end if

                        if (i2map .eq. orb(i1)%kion) then
                            symmat(i1 + (i2 - 1)*norb, is) = signmap*(indmap + (ncellmap - 1)*norb)
                        else
                            !  if there are more than one orbital of the same type acting on the same
                            !  atom determine its position ordorb (ordorb=1 typically)
                            i3 = 0
                            ordorb = 0
                            do while (i3 .lt. i1)
                                i3 = i3 + 1
                                if (orb(i3)%kion .eq. orb(i1)%kion .and. orb(i3)%comp .eq. orb(i1)%comp&
                                        &.and. orb(i3)%ioptorb .eq. orb(i1)%ioptorb) then
                                    ordorb = ordorb + 1
                                end if
                            end do

                            found = .false.
                            i3 = 0
                            ordorbtry = 0
                            do while (i3 .lt. norb .and. .not. found)
                                i3 = i3 + 1
                                if (orb(i3)%kion .eq. i2map .and. orb(i3)%comp .eq. orb(indmap)%comp&
                                        &.and. orb(i3)%ioptorb .eq. orb(indmap)%ioptorb) then
                                    ordorbtry = ordorbtry + 1
                                    if (ordorbtry .eq. ordorb) then
                                        found = .true.
                                        symmat(i1 + (i2 - 1)*norb, is) = signmap*(i3 + (ncellmap - 1)*norb)
                                    end if
                                end if
                            end do
                        end if
                    end if

                end do
            end do
        end do
        ! the constant for the Jastrow
        if (lda .eq. norb*ncell + 1) then
            do is = 1, nsym
                symmat(lda, is) = lda
            end do
        end if
    end subroutine update_symrot

    subroutine update_symtraj(nsym, symmat, lda, atoms_map, orb, norb)
        implicit none
        integer, intent(in) :: norb, lda, nsym, atoms_map(ntotatoms, nsym)
        type(orbital), intent(in) :: orb(norb)
        integer, intent(out) :: symmat(lda, nsym)
        integer ordorb, ordorbtry
        integer shiftorb
        shiftorb = 0
        if (lda .eq. norb*ncell + 1) shiftorb = 1
        !  Here the first is the constant
        symmat = 0
        do is = 1, nsym
            do i1 = 1, norb
                do i2 = 1, ncell
                    ionpos = orb(i1)%kion + (i2 - 1)*natoms
                    i2map = atoms_map(ionpos, is)
                    signmap = 1
                    if (i2map .lt. 0) then
                        signmap = -1
                        i2map = -i2map
                    end if
                    ncellmap = (i2map - 1)/natoms + 1
                    i2map = i2map - (ncellmap - 1)*natoms

                    indmap = i1 ! translation does not change type of orbital

                    if (i2map .eq. orb(i1)%kion) then
                        symmat(i1 + (i2 - 1)*norb + shiftorb, is) = signmap*(indmap + (ncellmap - 1)*norb + shiftorb)
                    else
                        !  if there are more than one orbital of the same type acting on the same
                        !  atom determine its position ordorb (ordorb=1 typically)
                        i3 = 0
                        ordorb = 0
                        do while (i3 .lt. i1)
                            i3 = i3 + 1
                            if (orb(i3)%kion .eq. orb(i1)%kion .and. orb(i3)%comp .eq. orb(i1)%comp&
                                    &.and. orb(i3)%ioptorb .eq. orb(i1)%ioptorb) then
                                ordorb = ordorb + 1
                            end if
                        end do

                        found = .false.
                        i3 = 0
                        ordorbtry = 0
                        do while (i3 .lt. norb .and. .not. found)
                            i3 = i3 + 1
                            if (orb(i3)%kion .eq. i2map .and. orb(i3)%comp .eq. orb(indmap)%comp&
                                    &.and. orb(i3)%ioptorb .eq. orb(indmap)%ioptorb) then
                                ordorbtry = ordorbtry + 1
                                if (ordorbtry .eq. ordorb) then
                                    found = .true.
                                    symmat(i1 + (i2 - 1)*norb + shiftorb, is) = signmap*(i3 + (ncellmap - 1)*norb + shiftorb)
                                end if
                            end if
                        end do
                    end if

                end do
            end do
        end do
        ! the constant for the Jastrow
        if (lda .eq. norb*ncell + 1) then
            do is = 1, nsym
                symmat(1, is) = 1
            end do
        end if
    end subroutine update_symtraj

    subroutine update_symtra(nsym, symmat, lda, atoms_map, orb, norb)
        implicit none
        integer, intent(in) :: norb, lda, nsym, atoms_map(ntotatoms, nsym)
        type(orbital), intent(in) :: orb(norb)
        integer, intent(out) :: symmat(lda, nsym)
        integer ordorb, ordorbtry
        symmat = 0
        do is = 1, nsym
            do i1 = 1, norb
                do i2 = 1, ncell
                    ionpos = orb(i1)%kion + (i2 - 1)*natoms
                    i2map = atoms_map(ionpos, is)
                    signmap = 1
                    if (i2map .lt. 0) then
                        signmap = -1
                        i2map = -i2map
                    end if
                    ncellmap = (i2map - 1)/natoms + 1
                    i2map = i2map - (ncellmap - 1)*natoms

                    indmap = i1 ! translation does not change type of orbital

                    if (i2map .eq. orb(i1)%kion) then
                        symmat(i1 + (i2 - 1)*norb, is) = signmap*(indmap + (ncellmap - 1)*norb)
                    else
                        !  if there are more than one orbital of the same type acting on the same
                        !  atom determine its position ordorb (ordorb=1 typically)
                        i3 = 0
                        ordorb = 0
                        do while (i3 .lt. i1)
                            i3 = i3 + 1
                            if (orb(i3)%kion .eq. orb(i1)%kion .and. orb(i3)%comp .eq. orb(i1)%comp&
                                    &.and. orb(i3)%ioptorb .eq. orb(i1)%ioptorb) then
                                ordorb = ordorb + 1
                            end if
                        end do

                        found = .false.
                        i3 = 0
                        ordorbtry = 0
                        do while (i3 .lt. norb .and. .not. found)
                            i3 = i3 + 1
                            if (orb(i3)%kion .eq. i2map .and. orb(i3)%comp .eq. orb(indmap)%comp&
                                    &.and. orb(i3)%ioptorb .eq. orb(indmap)%ioptorb) then
                                ordorbtry = ordorbtry + 1
                                if (ordorbtry .eq. ordorb) then
                                    found = .true.
                                    symmat(i1 + (i2 - 1)*norb, is) = signmap*(i3 + (ncellmap - 1)*norb)
                                end if
                            end if
                        end do
                    end if

                end do
            end do
        end do
        ! the constant for the Jastrow
        if (lda .eq. norb*ncell + 1) then
            do is = 1, nsym
                symmat(lda, is) = lda
            end do
        end if
    end subroutine update_symtra
    subroutine update_genjas
        implicit none
        allocate (recordsymj(2, nsym*ntra, nrecj + 2*nrecj_save), lenrecj(nrecj + 2*nrecj_save))
        allocate (jasyes(2*njastot, 2*njastot))
        jasyes(:, :) = .false.
        nrecj_ok = 0
        write (6, *) ' nelup nel inside update =', nelup, nel, const_term
        !  No constant term and nothing if nelup<2
        do j = 1, nrecj_save
            if (nelup .ge. 2 .and. abs(recordsymj_1(1, 1, j)) .ne. const_term&
                    &.and. abs(recordsymj_1(2, 1, j)) .ne. const_term) then
                nrecj_ok = nrecj_ok + 1
                lenrecj(nrecj_ok) = lenrecj_1(j)
                do k = 1, lenrecj_1(j)
                    do i = 1, 2
                        recordsymj(i, k, nrecj_ok) = recordsymj_1(i, k, j)
                    end do
                end do
            end if
        end do
        !  No constant term and nothing if neldo<2
        do j = 1, nrecj_save
            if (nel - nelup .ge. 2 .and. abs(recordsymj_1(1, 1, j)) .ne. const_term&
                    &.and. abs(recordsymj_1(2, 1, j)) .ne. const_term) then
                nrecj_ok = nrecj_ok + 1
                lenrecj(nrecj_ok) = lenrecj_1(j)
                do k = 1, lenrecj_1(j)
                    do i = 1, 2
                        recordsymj(i, k, nrecj_ok) = sign(1, recordsymj_1(i, k, j))*(abs(recordsymj_1(i, k, j)) + njastot)
                    end do
                end do
            end if
        end do

        do i = 1, nrecj
            accept = .false.
            if (nelup .ge. 1 .and. nel - nelup .ge. 1 .and. &
                !       accept always but the constant
                (abs(recordsymj_2(1, 1, i)) .ne. const_term .or. abs(recordsymj_2(2, 1, i)) .ne. const_term)) then
                accept = .true.
            elseif (nelup .eq. 0) then
                !  Only the constant term with the down
                if (abs(recordsymj_2(1, 1, i)) .eq. const_term .and. abs(recordsymj_2(2, 1, i)) .ne. const_term) accept = .true.
            elseif (nel - nelup .eq. 0) then
                !  Only the constant term with the up
                if (abs(recordsymj_2(2, 1, i)) .eq. const_term .and. abs(recordsymj_2(1, 1, i)) .ne. const_term) accept = .true.
            end if
            if (accept) then
                nrecj_ok = nrecj_ok + 1
                lenrecj(nrecj_ok) = lenrecj_2(i)
                do k = 1, lenrecj_2(i)
                    recordsymj(1, k, nrecj_ok) = recordsymj_2(1, k, i)
                    recordsymj(2, k, nrecj_ok) = sign(1, recordsymj_2(2, k, i))*(abs(recordsymj_2(2, k, i)) + njastot)
                end do
            end if
        end do

        do i = 1, nrecj_ok
            do k = 1, lenrecj(i)
                ix = abs(recordsymj(1, k, i))
                iy = abs(recordsymj(2, k, i))
                jasyes(ix, iy) = .true.
            end do
        end do
        njastot = 2*njastot
        nrecj = nrecj_ok
        deallocate (recordsymj_1, recordsymj_2, lenrecj_1, lenrecj_2, jasyes_1, jasyes_2)
    end subroutine update_genjas

end program makefort10

