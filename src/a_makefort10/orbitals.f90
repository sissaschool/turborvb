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

module mod_orbital
    use symm_data, only: nosym, nosym_contr, isymm, nsym, lunit, eps, nsym, write_log
    use constants, only: ipf
    implicit none
    logical yesmolat, yesmolatj, yesjas_atom, yesalloc_jas, real_contracted
    integer nnozero
    type atomstypes
        character(60) :: pseudo_name
        character(60) :: wf_name
        integer :: atomic_number
    end type
    double precision, allocatable :: atomic_jasmat(:, :)
    ! Orbital type
    type orbital
        integer ioptorb
        integer itype ! 1=s  3=p  5=d 7=f
        integer comp ! e.g. x,y,z if itype=3
        integer nparm ! # parameters
        integer kion ! # atom
        integer shell ! shell as in fort.10
        integer shellatom ! the shell on each atom
        integer norbatom ! nomber of the orbital on its own atom
        double precision, allocatable :: parms(:) ! orbital parameters
        double precision, allocatable :: pvec(:) !  used to rotate p or d orbital
    end type orbital

    type lsym_type
        integer :: neq
        integer, allocatable :: idx(:, :)
        integer :: goodsym
    end type

    type parsymm
        integer :: neq
        integer, allocatable :: idx(:)
        logical :: zpar
        integer :: ioptorb
        integer :: kion
        integer :: mult
        double precision :: value
    end type

    ! map an orbital in another using a given symmetry
    type orbmap
        ! iorb(i1,isym) = number of the orbital that corresponds to the orbital i1 with
        !                  symmetry isym, it contain the sign
        integer, allocatable :: iorb(:, :)

        ! flag to check if the symmetry is compatible with the orbital
        logical, allocatable :: comp(:, :)
    end type

    character(20) :: section_name

    integer :: nshelldet, nshelljas, njas_hyb, ndet_hyb, ncut_hyb
    logical :: no_3body_jas

    namelist /shells/ nshelljas, nshelldet, njas_hyb, ndet_hyb, no_3body_jas, ncut_hyb

contains

    subroutine read_orbitals(natoms, nel, noonebody, zeta, J3_off, cut_hybrid, detorb, ndetorb&
            &, totshelldet, ndetpar, jasorb, njasorb, totshelljas, njaspar, orbtype, jorbtype&
            &, funit, shiftbeta, complexfort10, symmagp, ncell, scale_jasfat)
        implicit none
        integer, intent(in) :: funit, natoms, nel, shiftbeta, ncell
        real(8), intent(in) :: zeta(2, natoms)
        logical, intent(out) :: J3_off(natoms)
        integer, intent(out) :: cut_hybrid(natoms)
        integer, intent(out) :: ndetorb, totshelldet, ndetpar, njasorb, totshelljas, njaspar
        type(orbital), pointer :: detorb(:), jasorb(:)
        character(20), intent(in) :: orbtype, jorbtype
        integer :: i1, i, j, ii, ix, iy, ishelljas, ishelldet
        integer :: iorbd, iorbj, iorbj_before, orbj_atom, nunc, nuncj, ncost
        character(6) :: intc
        logical :: complexfort10, symmagp
        logical, intent(in) :: noonebody
        real*8 jascost, jasfat, scale_jasfat
        yesmolat = .false.
        yesmolatj = .false.
        ndetpar = 0
        njaspar = 0
        ndetorb = 0
        njasorb = 0
        totshelljas = 0
        totshelldet = 0
        yesalloc_jas = .false.

        do i1 = 1, natoms
            write (6, *) ' atom number =', nint(zeta(1, i1)), i1, orbtype
            !
            nshelldet = 0
            nshelljas = 0
            njas_hyb = 0
            ndet_hyb = 0
            ncut_hyb = -1
            no_3body_jas = .false.
            !
            section_name = "ATOM_"//trim(intc(zeta(1, i1)))
            call findsection(funit, section_name)
            read (funit, nml=shells)
            if (nshelljas .lt. 0) then
                nshelljas = -nshelljas
                yesalloc_jas = .true.
            end if
            totshelldet = totshelldet + nshelldet + ndet_hyb
            totshelljas = totshelljas + nshelljas + njas_hyb
            J3_off(i1) = no_3body_jas
            cut_hybrid(i1) = ncut_hyb
            !  if(njas_hyb.gt.0) call errore('read_orbitals',"the automatical &
            !    &generation of hybrid jastrow orbitals hasn't been implemented!",1)

            ! computing dimension of the basis for DET and JAS
            call counting_orb(funit, nshelldet, ndet_hyb, orbtype, ndetorb, ndetpar, complexfort10, ipf)

            if (nshelljas .gt. 0) then
                read (funit, *) ! read a comment line between detorb and jasorb
                write (6, *) ' nshelljas, njas_hyb: read =', nshelljas, njas_hyb
                call counting_orb(funit, nshelljas, njas_hyb, jorbtype, njasorb, njaspar, .false., 1)
            else
                if (njas_hyb .gt. 0) call errore('read_orbitals', 'No jastrow shell, no jastrow hybrid!', 1)
            end if
            !
        end do

        write (*, *) ' # Det Shells : ', totshelldet
        write (*, *) ' # Jas Shells : ', totshelljas
        write (*, *) ' # Det Orbitals : ', ndetorb
        write (*, *) ' # Jas Orbitals : ', njasorb
        write (*, *) ' # Det Parameters : ', ndetpar
        write (*, *) ' # Jas Parameters : ', njaspar

        allocate (detorb(ndetorb))
        if (totshelljas .gt. 0) allocate (jasorb(njasorb))
        if (yesalloc_jas) then
            if (noonebody) then
                allocate (atomic_jasmat(njasorb, njasorb))
            else
                allocate (atomic_jasmat(njasorb + 1, njasorb + 1))
            end if
            atomic_jasmat = 0.d0
        end if

        ishelldet = 0
        ishelljas = 0
        iorbd = 1
        iorbj = 1
        nunc = 0
        nuncj = 0
        do i1 = 1, natoms
            write (6, *) ' atom number =', i1
            !
            nshelldet = 0
            nshelljas = 0
            njas_hyb = 0
            ndet_hyb = 0
            no_3body_jas = .false.
            yesjas_atom = .false.
            !
            section_name = "ATOM_"//trim(intc(zeta(1, i1)))
            call findsection(funit, section_name)
            read (funit, nml=shells, err=105, end=105)
            if (nshelljas .lt. 0) then
                yesjas_atom = .true.
                yesalloc_jas = .true.
                nshelljas = -nshelljas
            end if
            ! reading the basis set and filling detorb and jasorb
            call read_shells(detorb, iorbd, nunc, ishelldet, nshelldet, ndet_hyb, orbtype, i1, funit&
                    &, shiftbeta, yesmolat, complexfort10, ipf, symmagp)
            if (nshelljas .gt. 0) then
                if (scale_jasfat .ne. 0.d0) then
                    jasfat = (scale_jasfat - 1.d0)/dble(nel*ncell - 1) ! For H we take the H2 Jastrow

                else
                    jasfat = max((zeta(2, i1) - 1.d0), 1.d0)/dble(nel*ncell - 1) ! For H we take the H2 Jastrow
                end if

                iorbj_before = iorbj - 1
                read (funit, *, err=105, end=105) ! read a comment line between detorb and jasorb
                call read_shells(jasorb, iorbj, nuncj, ishelljas, nshelljas, njas_hyb, jorbtype, i1, funit&
                        &, shiftbeta, yesmolatj, .false., 1, symmagp)
                if (yesjas_atom) then
                    read (funit, *, err=105, end=105) ! read a comment line
                    read (funit, *) nnozero, ncost
                    if (ncost .eq. -1) then ! all one body term
                        do ii = 1, nnozero
                            read (funit, *) ix, jascost
                            atomic_jasmat(ix + iorbj_before, njasorb + 1) = jascost*jasfat
                            atomic_jasmat(njasorb + 1, ix + iorbj_before) = jascost*jasfat
                        end do
                    elseif (ncost .eq. 0) then ! all two body
                        do ii = 1, nnozero
                            read (funit, *) ix, iy, jascost
                            atomic_jasmat(ix + iorbj_before, iy + iorbj_before) = jascost
                        end do
                    else ! one body and two body
                        do ii = 1, nnozero
                            read (funit, *) ix, iy, jascost
                            if (iy .ne. ncost) then
                                atomic_jasmat(ix + iorbj_before, iy + iorbj_before) = jascost
                                atomic_jasmat(iy + iorbj_before, ix + iorbj_before) = jascost
                            else
                                atomic_jasmat(ix + iorbj_before, njasorb + 1) = jascost*jasfat
                                atomic_jasmat(njasorb + 1, ix + iorbj_before) = jascost*jasfat
                            end if
                        end do
                    end if
!    HERE load atomic_jasmat with the right scaling
                end if
            end if
            !
            ishelldet = ishelldet + nshelldet + ndet_hyb
            ishelljas = ishelljas + nshelljas + njas_hyb
            !
        end do
        if (ipf .eq. 2) then
            do i = 1, ndetorb
                if (detorb(i)%ioptorb .eq. 900000) then
                    do j = detorb(i)%nparm/4 + 1, detorb(i)%nparm/2
                        detorb(i)%parms(j) = detorb(i)%parms(j) + nunc*ncell
                    end do
                end if
            end do
        end if

        return
105     call errore("read_orbitals", " Error reading orbitals ! ", 1)
    end subroutine read_orbitals

    !-----------------------------------------------------------------
    ! This routine reads from standard input the basis set describing
    ! to each atom and computes the total # of parameters.
    !-----------------------------------------------------------------

    subroutine counting_orb(funit, nshell, nhyb, orbtype, norb, npar, complexfort10, ipf)
        implicit none
        integer, intent(in) :: funit, nshell, nhyb, ipf
        integer, intent(inout) :: norb, npar
        character(20), intent(in) :: orbtype
        !
        integer :: atomicnparm, nparm, nparm2, i2, i3, multi, ioptorb
        character(200) :: linedata
        !
        integer :: tempint
        double precision :: tempdb
        integer, external :: check_multioptorb
        logical complexfort10, cmplx_coeff

        ! COMPLEX DEB
        ! complex coefficients implemented only for 'normal' orbtype, so far...

        atomicnparm = 0
        cmplx_coeff = .false.
        !
        do i2 = 1, nshell
            read (funit, "(a)", err=100) linedata
            read (linedata, *, err=100) multi, nparm, ioptorb ! reading multiplicity and orbital type
            nparm2 = 0
            ! check to have a correct number of parameters
            if (.not. (mod(nparm, 2) .eq. 0 .or. nparm .eq. 1)) then
                call errore("counting_orb", " ERROR number of parameters (nparm) must be even ! ", 1)
            end if
            if (complexfort10 .and. nparm .gt. 1) cmplx_coeff = .true.
            !
            if (orbtype == "mixed" .and. ioptorb .ne. 900000) then
                read (linedata, *, err=100) multi, nparm, ioptorb, nparm2
                if (check_multioptorb(multi, ioptorb, 'W') .eq. 1)&
                        &write (6, *) ' Pay attention to the multiplicity!'
                write (6, *) ' shell number =', i2, multi, nparm, nparm2
            else
                if (check_multioptorb(multi, ioptorb, 'W') .eq. 1)&
                        &write (6, *) ' Pay attention to the multiplicity!'
                write (6, *) ' shell number =', i2, multi, nparm
            end if

            if (orbtype == "tempered" .and. ioptorb .ne. 900000) then
                read (funit, *, err=100) tempint, tempdb, tempdb
            elseif (orbtype == "mixed" .and. nparm2 .gt. 0 .and. ioptorb .ne. 900000) then
                read (funit, *, err=100) tempint, (tempdb, i3=1, nparm + 2)
            else
                ! orbtype=='normal'
                if (.not. cmplx_coeff .or. real_contracted) then
                    read (funit, *, err=100) tempint, (tempdb, i3=1, nparm)
                else ! contracted coefficients are considered complex
                    ! while exponents remain real.
                    ! in the case of uncontracted basis only detmat is complex
                    read (funit, *, err=100) tempint, (tempdb, i3=1, nparm*3/2)
                end if
            end if
            if (ioptorb .ne. 900000) atomicnparm = atomicnparm + (max(nparm/2, 1)*2 + nparm2)*multi
            if (ioptorb .eq. 900000 .and. ipf .eq. 2) then
                npar = npar + nparm
            end if
            nparm = nparm + nparm2
            npar = npar + nparm
            norb = norb + multi
        end do
        !
        ! adding hybrid orbitals

        npar = npar + atomicnparm*nhyb*ipf
        norb = norb + nhyb
        !
        return
100     call errore("counting_orb", " ERROR reading orbitals ! ", 1)
    end subroutine counting_orb

    ! -----------------------------------------------------------------------
    ! This routine stores all information related to orbitals of determinant
    ! and Jastrow, such as orbitals types, # of parameters and so on.
    ! -----------------------------------------------------------------------

    subroutine read_shells(orbitals, iorb, nunc, ishell, nshell, nhyb, orbtype, kion, funit&
            &, shiftbeta, yesmolat, complexfort10, ipf, symmagp)
        implicit none
        integer, intent(in) :: funit, nshell, nhyb, kion, ishell, shiftbeta, ipf
        character(20), intent(in) :: orbtype
        logical :: yesmolat, yeshyb
        integer :: iorb, indocc, nuncr, nunc
        type(orbital) :: orbitals(*)
        integer :: i1, i2, nparm, c1
        integer :: shellatom, inorbatom
        integer :: nparm2
        double precision :: alpha, beta
        double precision, allocatable :: tpar(:)
        character(200) :: linedata
        integer, external :: check_multioptorb
        logical complexfort10, cmplx_coeff, symmagp

        nuncr = nunc
        shellatom = 0
        inorbatom = 0
        cmplx_coeff = .false.
        yeshyb = .false.
        c1 = iorb

        do i1 = 1, nshell + nhyb
            if (i1 <= nshell) then
                read (funit, "(a)", err=101) linedata
                read (linedata, *, err=101) orbitals(c1)%itype, nparm, orbitals(c1)%ioptorb
                ! first line of each orbital section:
                ! itype = s,p,d... (multiplicity)
                ! nparm = # of parameters: exponents + contracted coefficients
                ! ioptorb = functional form of the orbital
                nparm2 = 0
                if (orbtype == "mixed" .and. orbitals(c1)%ioptorb .ne. 900000) then
                    read (linedata, *, err=101) orbitals(c1)%itype, nparm, orbitals(c1)%ioptorb, nparm2
                end if

                if (check_multioptorb(orbitals(c1)%itype, orbitals(c1)%ioptorb, 'W') .eq. 1)&
                        &write (6, *) ' Pay attention to the multiplicity!'

                if (.not. (mod(nparm, 2) == 0 .or. nparm == 1)) &
                        & call errore("read_shells", " Error nparm must be even!", 1)
                if (mod(nparm2, 2) .ne. 0) &
                        & call errore("read_shells", " Error when orbtype=mixed nparm2 must be even!", 1)

                if (complexfort10 .and. nparm .gt. 1) cmplx_coeff = .true.
                if (orbitals(c1)%ioptorb .ne. 900000 .or. ipf .eq. 1) then
                    if (.not. cmplx_coeff) then
                        allocate (orbitals(c1)%parms(nparm + nparm2))
                    else
                        allocate (orbitals(c1)%parms(nparm*3/2 + nparm2))
                    end if
                else ! below ipf=2  & orbital=900000
                    if (complexfort10) then
                        allocate (orbitals(c1)%parms(3*nparm))
                    else
                        allocate (orbitals(c1)%parms(2*nparm))
                    end if
                end if
                !
                !
                orbitals(c1)%parms(:) = 0.d0
                if (orbitals(c1)%ioptorb .ne. 900000) then
                    nunc = nunc + (max(nparm/2, 1) + nparm2/2)*orbitals(c1)%itype ! gaussians exponents
                end if

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! COMPLEX MODIF : complex coefficients implemented only for orbital type 'normal'
                ! so far. Necessary extension to 'tempered' and 'mixed' types. Eventually to hybrids too!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (orbtype == "tempered" .and. orbitals(c1)%ioptorb .ne. 900000) then
                    !
                    read (funit, *, err=101) orbitals(c1)%kion, alpha, beta
                    call adjustbeta(nparm, alpha, beta, shiftbeta)
                    call defz(nparm/2, alpha, beta, orbitals(c1)%parms)

                    if (cmplx_coeff) then
                        indocc = mod(2*(i1 - 1), nparm/2) + 1
                    else
                        indocc = mod(i1 - 1, nparm/2) + 1
                    end if

                    orbitals(c1)%parms(nparm/2 + indocc) = 1.d0
                elseif (orbtype == "mixed" .and. nparm2 .gt. 0 .and. orbitals(c1)%ioptorb .ne. 900000) then
                    allocate (tpar(nparm))
                    read (funit, *, err=101) orbitals(c1)%kion, (tpar(i2), i2=1, nparm), alpha, beta
                    call adjustbeta(nparm2, alpha, beta, shiftbeta)
                    call defz(nparm2/2, alpha, beta, orbitals(c1)%parms(nparm/2 + 1))

                    do i2 = 1, nparm/2
                        orbitals(c1)%parms(i2) = tpar(i2)
                        orbitals(c1)%parms((nparm + nparm2)/2 + i2) = tpar(nparm/2 + i2)
                    end do

                    if (nparm .eq. 0) then
                        indocc = mod(i1 - 1, nparm2/2) + 1
                        orbitals(c1)%parms(nparm2/2 + indocc) = 1.d0
                    elseif (sum(abs(tpar(nparm/2 + 1:nparm))) .eq. 0) then
                        indocc = mod(i1 - 1, nparm2/2) + 1
                        orbitals(c1)%parms(nparm2/2 + nparm + indocc) = 1.d0
                    end if

                    nparm = nparm + nparm2
                    deallocate (tpar)
                else ! orbtype='normal'
                    if (orbitals(c1)%ioptorb .ne. 900000 .or. ipf .eq. 1) then
                        if (.not. cmplx_coeff) then
                            read (funit, *, err=101) orbitals(c1)%kion, (orbitals(c1)%parms(i2), i2=1, nparm)
                        elseif (real_contracted) then
                            read (funit, *, err=101) orbitals(c1)%kion, (orbitals(c1)%parms(i2), i2=1, nparm/2) &
                                , (orbitals(c1)%parms(i2), i2=nparm/2 + 1, nparm*3/2, 2)
                            do i2 = nparm/2 + 2, nparm*3/2, 2
                                orbitals(c1)%parms(i2) = 0.d0
                            end do
                        else
                            read (funit, *, err=101) orbitals(c1)%kion, (orbitals(c1)%parms(i2), i2=1, nparm*3/2)
                        end if
                    else ! below ipf=2 and hybrid orbital
                        if (.not. cmplx_coeff) then
                            read (funit, *, err=101) orbitals(c1)%kion, (orbitals(c1)%parms(i2), i2=1, nparm/2) &
                                , (orbitals(c1)%parms(i2), i2=nparm + 1, nparm + nparm/2)
                            do i2 = 1, nparm/2
                                orbitals(c1)%parms(nparm/2 + i2) = orbitals(c1)%parms(i2)
                                orbitals(c1)%parms(nparm + nparm/2 + i2) = orbitals(c1)%parms(i2 + nparm)
                            end do
                        else
                            if (real_contracted) then
                                read (funit, *, err=101) orbitals(c1)%kion, (orbitals(c1)%parms(i2), i2=1, nparm/2) &
                                    , (orbitals(c1)%parms(i2), i2=nparm + 1, 2*nparm, 2)
                                do i2 = nparm + 2, 2*nparm, 2
                                    orbitals(c1)%parms(i2) = 0.d0
                                end do
                            else
                                read (funit, *, err=101) orbitals(c1)%kion, (orbitals(c1)%parms(i2), i2=1, nparm/2) &
                                    , (orbitals(c1)%parms(i2), i2=nparm + 1, 2*nparm)
                            end if
                            do i2 = 1, nparm/2
                                orbitals(c1)%parms(nparm/2 + i2) = orbitals(c1)%parms(i2)
                            end do
                            do i2 = 1, nparm
                                orbitals(c1)%parms(2*nparm + i2) = orbitals(c1)%parms(i2 + nparm)
                            end do
                        end if
                    end if
                end if
                !
                if (orbitals(c1)%ioptorb .eq. 900000) then
                    yeshyb = .true.
                end if
            else ! automatically generate hybrid orbitals
                orbitals(c1)%itype = 1
                orbitals(c1)%ioptorb = 900000
                orbitals(c1)%kion = 1
                nparm = (nunc - nuncr)*2 ! sick definition for hybrid jastrow
                if (complexfort10) then
                    allocate (orbitals(c1)%parms(3*ipf*(nparm/2)))
                else
                    allocate (orbitals(c1)%parms(ipf*nparm))
                end if
                orbitals(c1)%parms(:) = 0.d0
                do i2 = 1, nparm/2
                    orbitals(c1)%parms(i2) = i2
                end do
                if (ipf .eq. 2) then
                    do i2 = nparm/2 + 1, nparm
                        orbitals(c1)%parms(i2) = i2 - nparm/2
                    end do
                end if
                if (ipf .eq. 2 .and. symmagp) then
                    indocc = mod(i1 - nshell - 1, nparm/2) + 1
                else
                    indocc = mod(i1 - nshell - 1, ipf*nparm/2) + 1
                end if
                if (complexfort10) then
                    orbitals(c1)%parms(ipf*nparm/2 + 2*indocc - 1) = 1.d0
                    if (symmagp .and. ipf .eq. 2) orbitals(c1)%parms(nparm*3/2 + 2*indocc - 1) = 1.d0
                else
                    orbitals(c1)%parms(ipf*nparm/2 + indocc) = 1.d0
                    if (symmagp .and. ipf .eq. 2) orbitals(c1)%parms(nparm*3/2 + indocc) = 1.d0
                end if

                !    write(6,*) ' symmagp here =',symmagp

                !   if(orbitals(c1)%ioptorb.eq.900000) then
                !   do i2=1,ipf*(nparm/2)
                !     orbitals(c1)%parms(i2)=orbitals(c1)%parms(i2)+nuncr
                !   enddo
                !   !!write(6,*) "hyb", orbitals(c1)%parms(1:nparm/2)
                !   endif
                !   redefinition nparam
                !   nparm=ipf*nparm
            end if
            !
            !
            if (orbitals(c1)%ioptorb .eq. 900000) then
                do i2 = 1, ipf*(nparm/2)
                    orbitals(c1)%parms(i2) = orbitals(c1)%parms(i2) + nuncr
                end do
                !!write(6,*) "hyb", orbitals(c1)%parms(1:nparm/2)
                nparm = ipf*nparm
            end if
            !   redefinition nparam

            orbitals(c1)%kion = kion ! set the ions automaticaly
            shellatom = shellatom + 1
            orbitals(c1)%shellatom = shellatom
            inorbatom = inorbatom + 1
            orbitals(c1)%norbatom = inorbatom
            orbitals(c1)%nparm = nparm
            orbitals(c1)%shell = i1 + ishell
            !
            if (orbitals(c1)%itype == 1) then ! nothing to allocate for s
                !
                orbitals(c1)%comp = 1
                c1 = c1 + 1
                !
            else
                !
                allocate (orbitals(c1)%pvec(orbitals(c1)%itype)) ! the p-d-f-g vector
                orbitals(c1)%pvec(:) = 0.d0
                orbitals(c1)%pvec(1) = 1.d0
                ! **** Create the orbital for py and pz
                orbitals(c1)%comp = 1
                do i2 = 1, orbitals(c1)%itype - 1
                    allocate (orbitals(c1 + i2)%pvec(orbitals(c1)%itype))
                    orbitals(c1 + i2) = orbitals(c1)
                    orbitals(c1 + i2)%pvec(:) = 0.d0
                    orbitals(c1 + i2)%pvec(i2 + 1) = 1.d0
                    orbitals(c1 + i2)%shell = i1 + ishell
                    orbitals(c1 + i2)%kion = orbitals(c1)%kion
                    orbitals(c1 + i2)%shellatom = shellatom
                    orbitals(c1 + i2)%comp = i2 + 1
                    inorbatom = inorbatom + 1
                    orbitals(c1 + i2)%norbatom = inorbatom
                end do
                c1 = c1 + orbitals(c1)%itype
                !
            end if
        end do

        if (nhyb .gt. 0 .and. yeshyb) &
            write (6, *) "Warning: hybrid orbitals added both by hand and by automatic generation!"
        if (yesmolat .and. (nhyb .eq. 0 .and. .not. yeshyb)) then
            write (6, *) 'Warning hybrid orbitals have to be defined for all atoms or none!!!'
            !stop
        end if

        if (nhyb .gt. 0 .or. yeshyb) yesmolat = .true.

        iorb = c1

        return
101     call errore("read_shells", " Error reading atomic shells! ", 1)
    end subroutine read_shells

    subroutine adjustbeta(nparm, alpha, beta, shiftbeta)
        implicit none
        integer, intent(in) :: nparm, shiftbeta
        double precision :: alpha, beta
        !  update beta
        if (beta .lt. 0.d0) then
            if (nparm - 2*shiftbeta .ge. 2) then
                beta = abs(beta/alpha)**(1.d0/dble(nparm/2 - shiftbeta))
                if (shiftbeta .eq. -1) alpha = alpha*beta ! replacing the minimum alpha
            else
                beta = 1.d0
            end if
        elseif (beta .eq. 0.d0) then
            beta = 1.d0
        end if
    end subroutine adjustbeta

    subroutine read_atoms(natoms, izeta, atomic_detmat, detorb, ndetorb, nshelldet&
            &, ndetpar, jasorb, njasorb, nshelljas, njaspar, mytype, atypes, ntyp)

        use allio, only: zetar, atom_number, nshell_c, mult_c, nparam_c &
                &, nshellj_c, ioptorbj_c, dup_c, nnozero_c, iy, ix, occ_c &
                &, iesupind, jbraiesup_sav, jbraiesm_sav, iesmind, npar3body_c &
                &, nparamj_c, vju_c, iesup_c, multj_c, nozero_c, nelorb_c, scale_c &
                &, ioptorb_c, occj_c, nnozeroj_c, nelorbj_c, nozeroj_c, scalej_c
        implicit none
        integer, intent(out) :: ndetorb, nshelldet, ndetpar, njasorb, njaspar
        integer, intent(inout) :: nshelljas
        integer, intent(in) :: natoms
        real(8), intent(in) :: izeta(2, natoms)
        integer, intent(in) :: ntyp, mytype(natoms)
        type(atomstypes) :: atypes(ntyp)
        type(orbital), pointer :: detorb(:), jasorb(:)
        double precision, pointer :: atomic_detmat(:, :)
        integer, parameter :: aunit = 10
        logical :: readjastrow
        integer :: i1, i2, i3, iorb, ijorb, inorbatom
        integer :: ishell, ishellj, dpos, norbatom, indpar
        integer :: ppos

        yesmolat = .false.

        readjastrow = .false.
        if (nshelljas .ne. 0) readjastrow = .true.

        ndetorb = 0
        njasorb = 0
        nshelljas = 0
        nshelldet = 0
        ndetpar = 0
        njaspar = 0

        if (write_log) then
            write (lunit, *)
            write (lunit, *) ' Reading Atoms '
        end if

        do i1 = 1, natoms
            !
            if (nint(izeta(2, i1)) == 0) then
                if (write_log) write (lunit, *) ' Wave-function for atom ', i1, ' skipped '
            else
                !
                if (write_log) write (lunit, *) 'For Atom ', i1, ' reading file ', atypes(mytype(i1))%wf_name
                open (file=trim(atypes(mytype(i1))%wf_name), unit=aunit, status="old", form="formatted")
                !
                call read_fort10
                !
                if (nint(zetar(1)) .ne. nint(izeta(2, i1)) .or. nint(atom_number(1)) .ne. nint(izeta(1, i1))) then
                    call errore("read_atoms", &
                                "Error! Zeta in "//trim(atypes(mytype(i1))%wf_name)//" different from the one in the input ", 1)
                end if
                nshelldet = nshelldet + nshell_c
                do i2 = 1, nshell_c
                    ndetorb = ndetorb + mult_c(i2)
                    ndetpar = ndetpar + nparam_c(i2)
                end do
                !
                if (readjastrow .and. nshellj_c .ne. 0) then
                    if (any(ioptorbj_c(:) == 200)) then
                        nshelljas = nshelljas + nshellj_c - 1
                    else
                        nshelljas = nshelljas + nshellj_c
                    end if
                    do i2 = 1, nshellj_c
                        if (ioptorbj_c(i2) .ne. 200) then
                            njasorb = njasorb + multj_c(i2)
                            njaspar = njaspar + nparamj_c(i2)
                        end if
                    end do
                end if
                !
                call deallocate_all
                close (aunit)
                !
            end if
            !
        end do
        !
        ! add one orbital for the constant one (ioptorb=200)
        !  if(readjastrow.and.njasorb.ne.0) njasorb=njasorb+1
        !
        allocate (detorb(ndetorb), atomic_detmat(ndetorb, ndetorb))
        atomic_detmat = 0.d0
        if (njasorb .ne. 0) then
            allocate (jasorb(njasorb))
        end if

        iorb = 1
        ijorb = 1
        ishell = 1
        ishellj = 1

        dpos = 0
        ppos = 0

        do i1 = 1, natoms
            !
            if (nint(izeta(2, i1)) .ne. 0) then
                !
                open (file=trim(atypes(mytype(i1))%wf_name), unit=aunit, status="old", form="formatted")
                !
                call read_fort10
                !
                indpar = 0
                norbatom = 0
                inorbatom = 1
                !
                do i2 = 1, nshell_c
                    detorb(iorb)%norbatom = inorbatom
                    detorb(iorb)%itype = mult_c(i2)
                    detorb(iorb)%nparm = nparam_c(i2)
                    detorb(iorb)%ioptorb = ioptorb_c(i2)
                    detorb(iorb)%kion = i1
                    detorb(iorb)%shell = ishell
                    detorb(iorb)%shellatom = i2
                    allocate (detorb(iorb)%parms(nparam_c(i2)))
                    allocate (detorb(iorb)%pvec(detorb(iorb)%itype))
                    detorb(iorb)%pvec = 0
                    detorb(iorb)%pvec(1) = 1
                    !
                    do i3 = 1, nparam_c(i2)
                        detorb(iorb)%parms(i3) = dup_c(indpar + i3)
                    end do
                    !
                    do i3 = 2, mult_c(i2), 1
                        inorbatom = inorbatom + 1
                        detorb(iorb + i3 - 1)%norbatom = inorbatom
                        detorb(iorb + i3 - 1)%itype = detorb(iorb)%itype
                        detorb(iorb + i3 - 1)%nparm = detorb(iorb)%nparm
                        detorb(iorb + i3 - 1)%ioptorb = detorb(iorb)%ioptorb
                        detorb(iorb + i3 - 1)%kion = detorb(iorb)%kion
                        detorb(iorb + i3 - 1)%shellatom = i2
                        detorb(iorb + i3 - 1)%shell = ishell
                        allocate (detorb(iorb + i3 - 1)%parms(nparam_c(i2)))
                        detorb(iorb + i3 - 1)%parms = detorb(iorb)%parms
                        allocate (detorb(iorb + i3 - 1)%pvec(detorb(iorb + i3 - 1)%itype))
                        detorb(iorb + i3 - 1)%pvec = 0
                        detorb(iorb + i3 - 1)%pvec(i3) = 1
                    end do
                    iorb = iorb + mult_c(i2)
                    norbatom = norbatom + mult_c(i2)
                    indpar = indpar + nparam_c(i2)
                    ishell = ishell + 1
                    inorbatom = inorbatom + 1
                    !
                end do
                !
                do i2 = 1, nnozero_c
                    iy = (nozero_c(i2) - 1)/nelorb_c + 1
                    ix = nozero_c(i2) - (iy - 1)*nelorb_c
                    if (iy .gt. norbatom) then
                        if (write_log) write (lunit, *) ' WARNING! Skipping detmat for unpaired electrons! '
                    else
                        atomic_detmat(dpos + ix, dpos + iy) = scale_c(i2)
                    end if
                end do
                !
                dpos = dpos + occ_c
                ppos = ppos + iesup_c ! position of Z indices
                !
                call deallocate_all
                close (aunit)
                !
            end if
            !
        end do
        !
        if (readjastrow) then
            !
            ! dpos=0

            jasorb(:)%ioptorb = 0

            do i1 = 1, natoms
                !
                if (nint(izeta(2, i1)) .ne. 0) then
                    !
                    open (file=trim(atypes(mytype(i1))%wf_name), unit=aunit, status="old", form="formatted")
                    !
                    call read_fort10

                    indpar = 0
                    norbatom = 0
                    inorbatom = 1
                    !
                    do i2 = 1, nshellj_c
                        ! do not read two time the constant orbital 200
                        if (ioptorbj_c(i2) == 200) cycle
                        !
                        write (*, *) njasorb, ijorb, nshellj_c
                        !
                        jasorb(ijorb)%itype = multj_c(i2)
                        jasorb(ijorb)%nparm = nparamj_c(i2)
                        jasorb(ijorb)%ioptorb = ioptorbj_c(i2)
                        jasorb(ijorb)%kion = i1
                        jasorb(ijorb)%shell = ishellj
                        jasorb(ijorb)%shellatom = i2
                        allocate (jasorb(ijorb)%parms(nparamj_c(i2)))
                        allocate (jasorb(ijorb)%pvec(jasorb(ijorb)%itype))
                        jasorb(ijorb)%norbatom = inorbatom
                        jasorb(ijorb)%pvec = 0
                        jasorb(ijorb)%pvec(1) = 1
                        !
                        do i3 = 1, nparamj_c(i2)
                            jasorb(ijorb)%parms(i3) = vju_c(indpar + i3)
                        end do
                        !
                        do i3 = 2, multj_c(i2), 1
                            inorbatom = inorbatom + 1
                            jasorb(ijorb + i3 - 1)%norbatom = inorbatom
                            jasorb(ijorb + i3 - 1)%itype = jasorb(ijorb)%itype
                            jasorb(ijorb + i3 - 1)%nparm = jasorb(ijorb)%nparm
                            jasorb(ijorb + i3 - 1)%ioptorb = jasorb(ijorb)%ioptorb
                            jasorb(ijorb + i3 - 1)%kion = jasorb(ijorb)%kion
                            jasorb(ijorb + i3 - 1)%shellatom = i2
                            jasorb(ijorb + i3 - 1)%shell = ishellj
                            allocate (jasorb(ijorb + i3 - 1)%parms(nparamj_c(i2)))
                            jasorb(ijorb + i3 - 1)%parms = jasorb(ijorb)%parms
                            allocate (jasorb(ijorb + i3 - 1)%pvec(3))
                            jasorb(ijorb + i3 - 1)%pvec = 0
                            jasorb(ijorb + i3 - 1)%pvec(i3) = 1
                        end do
                        ijorb = ijorb + multj_c(i2)
                        norbatom = norbatom + multj_c(i2)
                        indpar = indpar + nparamj_c(i2)
                        ishellj = ishellj + 1
                        inorbatom = inorbatom + 1
                        !
                    end do
                    !
                    ! atomic_jasmat = 0.
                    !    maxy=0
                    !    do i2=1,nnozeroj_c
                    !      if(maxy<(nozeroj_c(i2)-1)/nelorbj_c+1) &
                    !        & maxy=(nozeroj_c(i2)-1)/nelorbj_c+1
                    !    enddo
                    !    !
                    !    do i2=1,nnozeroj_c
                    !      iy=(nozeroj_c(i2)-1)/nelorbj_c+1
                    !      ix=nozeroj_c(i2)-(iy-1)*nelorbj_c
                    !      if(iy.eq.ix) then
                    !              atomic_jasmat(dpos+ix,1)=scalej_c(i2)
                    !       endif
                    !      if(iy==maxy) then
                    !         atomic_jasmat(dpos+ix,2)=scalej_c(i2)
                    !      endif
                    !    enddo
                    !    !
                    !    dpos=dpos+occj_c-1
                    !
                    call deallocate_all
                    close (aunit)
                    !
                end if
                !
            end do
        end if

        write (*, *)
        write (*, *)
        write (*, *) ' Number of determinant shells : ', nshelldet
        write (*, *) ' Number of determinant orbitals : ', ndetorb
        write (*, *) ' Number of parameters in determinant : ', ndetpar

        write (*, *) ' Number of Jastrow shells : ', nshelljas
        write (*, *) ' Number of Jastrow orbitals : ', njasorb
        write (*, *) ' Number of parameters in Jastrow : ', njaspar

        if (write_log) then
            write (lunit, *) ' Atomic Lambda Matrix '
            do i1 = 1, ndetorb
                do i2 = i1, ndetorb
                    write (lunit, '(2i3,f12.6)') i1, i2, atomic_detmat(i1, i2)
                end do
            end do
        end if
    end subroutine read_atoms

    subroutine print_orbitals(orbitals, norbitals, orbname)
        implicit none
        integer, intent(in) :: norbitals
        character(*), intent(in) :: orbname
        type(orbital), intent(in) :: orbitals(norbitals)
        integer :: i1

        if (write_log) then
            write (lunit, *) '# read orbitals ', orbname
            do i1 = 1, norbitals
                write (lunit, '(a,i5)') "---> number ", i1
                write (lunit, '(a,i5)') "     type ", orbitals(i1)%itype
                write (lunit, '(a,i5)') "     norbatom ", orbitals(i1)%norbatom
                write (lunit, '(a,i5)') "     shell ", orbitals(i1)%shell
                write (lunit, '(a,i5)') "     atom-shell ", orbitals(i1)%shellatom
                write (lunit, '(a,i5)') "     kion ", orbitals(i1)%kion
                if (orbitals(i1)%itype .eq. 3) write (lunit, '(a,3f12.4)') "     p-vec ", orbitals(i1)%pvec
            end do
            write (lunit, *)
        end if
    end subroutine print_orbitals

    subroutine apply_symm_to_orbitals(orb, norb, orb_map)
        implicit none
        integer, intent(in) :: norb
        type(orbmap), intent(out) :: orb_map
        type(orbital), intent(in) :: orb(norb)
        integer :: i1, i2, is, j, jj, k
        double precision :: dvec(9), mat(3, 3), ut(3, 3), matf(3, 3, 3), matg(3, 3, 3, 3)
        double precision :: matblok(3, 3), matscra(3, 3, 3), vecscra(3), vecscra2(3)
        double precision :: matscrag(3, 3, 3, 3), matconv(3, 3)
        logical :: found

        ! map an orbital in another using a given symmetry
        ! orb_map%iorb(i1,isym) = number of the orbital that corresponds
        !                         to the orbital i1 with  symmetry isym,
        !                         it contain the sign
        ! orb_map%comp(i1,isym) = flag to check if the symmetry is compatible
        !                         with the orbital, namely if under the symmetry
        !                         isymm(:,:,is) an orbital is transformed in another
        !                         one (ex. px-> py), TurboRVB do not support yet
        !                         linear combination of orbitals in the symmetries

        allocate (orb_map%iorb(norb, nsym), orb_map%comp(norb, nsym))

        !write(6,*) ' nsym found by pw ',nsym

        do is = 1, nsym
            !  write(6,*) ' sym matrix =',is
            !  do i1=1,3
            !  write(6,*) isymm(:,i1,is)
            !  enddo
            matconv(:, :) = isymm(:, :, is)
            do i1 = 1, norb
                found = .false.
                i2 = 1
                dvec = 0.d0
                !
                !******** Rotate p and d orbitals *************
                if (orb(i1)%itype == 3) then ! p-orbital
                    dvec(1:3) = matmul(isymm(:, :, is), orb(i1)%pvec(1:3))
                elseif (orb(i1)%itype == 5) then ! d-orbital
                    call mat2d(mat, orb(i1)%pvec, -1)
                    call dgemm('N', 'N', 3, 3, 3, 1.d0, matconv, 3, mat, 3, 0.d0, ut, 3)
                    call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, matconv, 3, 0.d0, mat, 3)
                    call mat2d(mat, dvec, 1)
                elseif (orb(i1)%itype == 7) then ! f-orbital
                    call mat3f(matf, orb(i1)%pvec, -1)
                    do j = 1, 3
                        matblok(:, :) = matf(j, :, :)
                        call dgemm('N', 'N', 3, 3, 3, 1.d0, matconv, 3, matblok, 3, 0.d0, ut, 3)
                        call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, matconv, 3, 0.d0, matblok, 3)
                        matscra(j, :, :) = matblok(:, :)
                    end do
                    !
                    do j = 1, 3
                        do jj = 1, 3
                            vecscra(:) = matscra(:, j, jj)
                            call dgemv('N', 3, 3, 1.d0, matconv, 3, vecscra(:), 1, 0.d0, vecscra2, 1)
                            matf(:, j, jj) = vecscra2(:)
                        end do
                    end do
                    call mat3f(matf, dvec, 1)
                    !
                elseif (orb(i1)%itype == 9) then ! g-orbital
                    call mat2g(matg, orb(i1)%pvec, -1)
                    do j = 1, 3
                        do k = 1, 3
                            matblok(:, :) = matg(j, k, :, :)
                            call dgemm('N', 'N', 3, 3, 3, 1.d0, matconv, 3, matblok, 3, 0.d0, ut, 3)
                            call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, matconv, 3, 0.d0, matblok, 3)
                            matscrag(j, k, :, :) = matblok(:, :)
                        end do
                    end do
                    !
                    !
                    do j = 1, 3
                        do k = 1, 3
                            matblok(:, :) = matscrag(:, :, j, k)
                            call dgemm('N', 'N', 3, 3, 3, 1.d0, matconv, 3, matblok, 3, 0.d0, ut, 3)
                            call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, matconv, 3, 0.d0, matblok, 3)
                            matscrag(:, :, j, k) = matblok(:, :)
                        end do
                    end do
                    !
                    !         Cristhmas BUG (2008) by Claudio
                    !         call mat2g(matscra,dvec,1)
                    call mat2g(matscrag, dvec, 1)
                    !
                end if
                !
                do while (i2 .le. norb)
                    if (orb(i1)%shellatom == orb(i2)%shellatom .and. orb(i2)%kion == orb(i1)%kion) then
                        !
                        if (orb(i1)%itype == 1) then ! s-orbital
                            orb_map%iorb(i1, is) = i2
                            orb_map%comp(i1, is) = .true.
                            found = .true.
                        elseif (orb(i1)%itype .gt. 1 .and. orb(i1)%itype .le. 9) then ! p-d-f-g-orbital
                            if (sum(abs(orb(i2)%pvec(:) - dvec(1:orb(i1)%itype))) .lt. eps) then
                                orb_map%iorb(i1, is) = i2
                                orb_map%comp(i1, is) = .true.
                                found = .true.
                            elseif (sum(abs(orb(i2)%pvec(:) + dvec(1:orb(i1)%itype))) .lt. eps) then
                                orb_map%iorb(i1, is) = -i2
                                orb_map%comp(i1, is) = .true.
                                found = .true.
                            end if
                        end if
                    end if
                    i2 = i2 + 1
                end do
                ! symmetry transformation not campatible with this orbital
                if (.not. found) then
                    orb_map%iorb(i1, is) = i1
                    orb_map%comp(i1, is) = .false.
                end if
            end do
        end do

    end subroutine apply_symm_to_orbitals

    subroutine apply_symm_to_forces(orb_map)
        implicit none
        integer, intent(out) :: orb_map(3, nsym)
        integer :: i1, i2, is, j, jj, k
        double precision :: dvec(3), dvecf(3), mat(3, 3), ut(3, 3)
        logical :: found

        orb_map = 0

        do is = 1, nsym
            mat(:, :) = isymm(:, :, is)
            do i1 = 1, 3
                found = .false.
                i2 = 1
                dvec = 0.d0
                dvecf = 0.d0
                dvecf(i1) = 1.d0
                !
                !******** Rotate component i1
                dvec(1:3) = matmul(mat, dvecf)
                !
                do i2 = 1, 3
                    if (abs(abs(dvec(i2)) - 1.d0) .lt. 1d-7) then
                        found = .true.
                        if (dvec(i2) .gt. 0.d0) then
                            orb_map(i1, is) = i2
                        else
                            orb_map(i1, is) = -i2
                        end if
                    end if
                end do
                ! symmetry transformation not campatible with this orbital
            end do
        end do

    end subroutine apply_symm_to_forces

    subroutine generate_orbidx(zeta, natoms, newidx, orb, norb, indion)
        implicit none
        integer, intent(in) :: norb, natoms, indion(natoms)
        real(8), intent(in) :: zeta(2, natoms)
        integer, intent(out) :: newidx(norb)
        type(orbital) :: orb(norb)
        integer :: i1, i2
        logical :: found

        ! For a given parameter  select the one compatible with simmetry with the
        ! largest parameter number < given parameter.
        do i1 = 1, norb
            found = .false.
            i2 = 1
            do while (i2 .le. i1 .and. .not. found)
                if (all(nint(zeta(:, orb(i1)%kion)) == nint(zeta(:, orb(i2)%kion))) .and.&
                        & orb(i1)%norbatom == orb(i2)%norbatom .and.&
                        &indion(orb(i1)%kion) .eq. indion(orb(i2)%kion)) then
                    found = .true.
                    newidx(i1) = i2
                end if
                !
                i2 = i2 + 1
            end do
            !
            if (.not. found) newidx(i1) = i1
            !
        end do

        if (write_log) then
            write (lunit, *)
            write (lunit, *) ' Orbital idx '

            do i1 = 1, norb
                write (lunit, *) i1, ' --> ', newidx(i1)
            end do
        end if
    end subroutine generate_orbidx

    subroutine par_symm(eqpar, orb, orbidx, orbidxz, neq, ntotpar, ncell, norb, onlycontr, ipf, symmagp)
        implicit none
        integer, intent(in) :: ncell, norb, ntotpar, ipf
        type(orbital), intent(in) :: orb(norb)
        integer, intent(in) :: orbidx(norb), orbidxz(norb)
        logical, intent(in) :: onlycontr
        integer, intent(out) :: neq
        type(parsymm), intent(out) :: eqpar(ntotpar)
        integer :: maxpar, ntothyb, ntotnohyb
        integer :: i1, i2, ip, ic, i3, ind, i1old, iold
        integer :: tneq
        logical :: found
        logical :: yes9
        logical :: symmagp
        logical :: dosym
        integer :: neqsign
        integer, dimension(:), allocatable :: pos, poseq, tmp

        allocate (pos(norb), poseq(norb), tmp(ntotpar*ncell))

        maxpar = maxval(orb(:)%nparm)

        !write(6,*) ' symmagp inside =',symmagp,dosym

        !write(6,*) ' Input  ntotpar ',ntotpar
        !write(6,*) ' symmagp inside =',symmagp,dosym
        !write(6,*) ' Input orb '
        !do i1=1,norb
        !write(6,*) i1,orb(i1)%ioptorb,orb(i1)%nparm
        !enddo

        eqpar(:)%zpar = .false.

        neq = 0
        ind = 0
        pos = 0
        i1old = 0
        i1 = 0
        ! Take the first non-hybrid orbital
        do while (i1old .eq. 0 .and. i1 .lt. norb)
            i1 = i1 + 1
            if (orb(i1)%ioptorb .ne. 900000) i1old = i1
        end do
        i1 = 0
        yes9 = .false.
        iold = 0
        do while (.not. yes9 .and. i1 .lt. norb)
            i1 = i1 + 1
            if (orb(i1)%ioptorb .eq. 900000) then
                yes9 = .true.
                iold = i1
            end if
        end do

        if (symmagp .and. ipf .eq. 2 .and. .not. nosym_contr .and. yes9) then

            dosym = .true.

        else

            dosym = .false.

        end if

        ! First the non hybrid
        do i1 = i1old + 1, norb
            if (orb(i1)%ioptorb .ne. 900000) then
                if (orb(i1)%shell .ne. orb(i1old)%shell) then
                    pos(i1) = pos(i1old) + orb(i1old)%nparm
                else
                    pos(i1) = pos(i1old)
                end if
                i1old = i1
            end if
        end do
        ! Then the hybrid
        do i1 = 1, norb
            if (orb(i1)%ioptorb .eq. 900000) then
                if (orb(i1)%shell .ne. orb(i1old)%shell) then
                    pos(i1) = pos(i1old) + orb(i1old)%nparm
                else
                    pos(i1) = pos(i1old)
                end if
                i1old = i1
            end if
        end do

        poseq = pos

        if (yes9) then
            ntotnohyb = pos(iold)
            ntothyb = ntotpar - pos(iold)

            do i1 = 1, norb
                if (orb(i1)%ioptorb .eq. 900000) pos(i1) = pos(i1) + (ncell - 1)*ntotnohyb
            end do
        else
            ntothyb = ntotpar
            ntotnohyb = ntotpar
        end if
        !write(6,*) ' Position orbitals found ',yes9,ntothyb,ntotnohyb
        !do i1=1,norb
        !write(6,*) i1,orb(i1)%ioptorb,pos(i1)
        !enddo

        eqpar(:)%neq = 0

        if (dosym) then
            do i1 = 1, ntotpar
                allocate (eqpar(i1)%idx(2*ncell*ntotpar))
            end do
        else
            do i1 = 1, ntotpar
                allocate (eqpar(i1)%idx(ncell*ntotpar))
            end do
        end if

        !
        !********* Put the paramters of the first orbital
        do ip = 1, orb(1)%nparm
            if (ip .le. (orb(1)%nparm + 1)/2 .and. orb(1)%ioptorb .ne. 900000) then
                eqpar(ip + poseq(1))%zpar = .true.
                eqpar(ip + poseq(1))%ioptorb = orb(1)%ioptorb
                eqpar(ip + poseq(1))%value = orb(1)%parms(ip)
                eqpar(ip + poseq(1))%kion = orb(1)%kion
                eqpar(ip + poseq(1))%mult = orb(1)%itype
            end if

            if (orb(1)%ioptorb .ne. 900000 .or. ip .gt. (orb(1)%nparm + 1)/2) then
                if (dosym .and. orb(1)%ioptorb .eq. 900000 .and. ip .le. orb(1)%nparm*3/4) then
                    eqpar(ip + poseq(1))%neq = 2*ncell
                    do ic = 1, ncell
                        eqpar(ip + poseq(1))%idx(ic) = ip + pos(1) + (ic - 1)*ntothyb
                    end do
                    do ic = 1, ncell
                        eqpar(ip + poseq(1))%idx(ic + ncell) = ip + pos(1) + orb(1)%nparm/4 + (ic - 1)*ntothyb
                    end do
                elseif (.not. dosym .or. orb(1)%ioptorb .ne. 900000) then
                    eqpar(ip + poseq(1))%neq = ncell
                    if (orb(1)%ioptorb .ne. 900000) then
                        do ic = 1, ncell
                            eqpar(ip + poseq(1))%idx(ic) = ip + pos(1) + (ic - 1)*ntotnohyb
                        end do
                    else
                        do ic = 1, ncell
                            eqpar(ip + poseq(1))%idx(ic) = ip + pos(1) + (ic - 1)*ntothyb
                        end do
                    end if
                end if
            end if
        end do
        if (orb(1)%ioptorb .ne. 900000) then
            neq = neq + orb(1)%nparm
        else
            if (dosym) then
                neq = neq + orb(1)%nparm/4
            else
                neq = neq + orb(1)%nparm/2
            end if
        end if

        !
        !********* Loop on the other orbitals to check symmetries
        !        First the exponents
        !
        do i1 = 2, norb
            !
            if (orb(i1)%shell .ne. orb(i1 - 1)%shell) then
                !
                i2 = iabs(orbidxz(i1))
                if (i2 .lt. i1) then
                    !      Update

                    !
                    do ip = 1, orb(i1)%nparm
                        if (ip .le. (orb(i1)%nparm + 1)/2 .and. orb(i1)%ioptorb .ne. 900000) then
                            eqpar(ip + poseq(i2))%zpar = .true.
                            eqpar(ip + poseq(i2))%ioptorb = orb(i1)%ioptorb
                            eqpar(ip + poseq(i2))%value = orb(i1)%parms(ip)
                            eqpar(ip + poseq(i2))%kion = orb(i1)%kion
                            eqpar(ip + poseq(i2))%mult = orb(i1)%itype
                        end if
                        if (orb(i1)%ioptorb .ne. 900000 .and. ip .le. (orb(i1)%nparm + 1)/2) then
                            do ic = 1, ncell
                                eqpar(ip + poseq(i2))%idx(eqpar(ip + poseq(i2))%neq + ic) = ip + pos(i1) + (ic - 1)*ntotnohyb
                            end do
                            eqpar(ip + poseq(i2))%neq = eqpar(ip + poseq(i2))%neq + ncell
                        end if
                    end do
                    !
                else
                    !
                    !      Begin new symmetry parameter record

                    do ip = 1, orb(i1)%nparm
                        if (ip .le. (orb(i1)%nparm + 1)/2 .and. orb(i1)%ioptorb .ne. 900000) then
                            eqpar(ip + poseq(i1))%zpar = .true.
                            eqpar(ip + poseq(i1))%ioptorb = orb(i1)%ioptorb
                            eqpar(ip + poseq(i1))%value = orb(i1)%parms(ip)
                            eqpar(ip + poseq(i1))%kion = orb(i1)%kion
                            eqpar(ip + poseq(i1))%mult = orb(i1)%itype
                        end if
                        if (orb(i1)%ioptorb .ne. 900000 .and. ip .le. (orb(i1)%nparm + 1)/2) then
                            eqpar(ip + poseq(i1))%neq = ncell
                            do ic = 1, ncell
                                eqpar(ip + poseq(i1))%idx(ic) = ip + pos(i1) + (ic - 1)*ntotnohyb
                            end do
                        end if
                    end do
                    ! this number is #Eq. Det Atomic or #Eq. Det. Jastrow
                    if (orb(i1)%ioptorb .ne. 900000) then
                        neq = neq + orb(i1)%nparm
                    else
                        neq = neq + orb(i1)%nparm/2
                    end if
                end if
                !
            end if
            !
        end do

        !********* Loop on the other orbitals to check symmetries

        !    Then the coefficients
        !
        do i1 = 2, norb
            !
            if (orb(i1)%shell .ne. orb(i1 - 1)%shell) then
                !
                if (nosym_contr) then
                    i2 = i1
                else
                    i2 = iabs(orbidx(i1))
                end if
                if (i2 .lt. i1) then
                    !      Update
                    do ip = 1, orb(i1)%nparm

                        if (ip .gt. (orb(i1)%nparm + 1)/2) then
                            if (orb(i1)%ioptorb .ne. 900000) then
                                do ic = 1, ncell
                                    eqpar(ip + poseq(i2))%idx(eqpar(ip + poseq(i2))%neq + ic) = ip + pos(i1) + (ic - 1)*ntotnohyb
                                end do
                                eqpar(ip + poseq(i2))%neq = eqpar(ip + poseq(i2))%neq + ncell
                            else
                                if (dosym .and. ip .le. orb(i1)%nparm*3/4) then
                                    do ic = 1, ncell
                                        eqpar(ip + poseq(i2))%idx(eqpar(ip + poseq(i2))%neq + ic) = ip + pos(i1) + (ic - 1)*ntothyb
                                    end do
                                    do ic = 1, ncell
                                        eqpar(ip + poseq(i2))%idx(eqpar(ip + poseq(i2))%neq + ic + ncell) &
                                            = ip + pos(i1) + orb(i1)%nparm/4 + (ic - 1)*ntothyb
                                    end do
                                    eqpar(ip + poseq(i2))%neq = eqpar(ip + poseq(i2))%neq + 2*ncell

                                elseif (.not. dosym) then
                                    do ic = 1, ncell
                                        eqpar(ip + poseq(i2))%idx(eqpar(ip + poseq(i2))%neq + ic) &
                                            = ip + pos(i1) + (ic - 1)*ntothyb
                                    end do
                                    eqpar(ip + poseq(i2))%neq = eqpar(ip + poseq(i2))%neq + ncell
                                end if
                            end if
                        end if
                    end do
                    !
                else
                    !
                    !      Begin new symmetry parameter record

                    do ip = 1, orb(i1)%nparm
                        if (ip .gt. (orb(i1)%nparm + 1)/2) then
                            if (orb(i1)%ioptorb .ne. 900000) then
                                do ic = 1, ncell
                                    eqpar(ip + poseq(i1))%idx(ic) = ip + pos(i1) + (ic - 1)*ntotnohyb
                                end do
                                eqpar(ip + poseq(i1))%neq = ncell
                            else
                                if (dosym .and. ip .le. orb(i1)%nparm*3/4) then
                                    do ic = 1, ncell
                                        eqpar(ip + poseq(i1))%idx(ic) = ip + pos(i1) + (ic - 1)*ntothyb
                                    end do
                                    do ic = 1, ncell
                                        eqpar(ip + poseq(i1))%idx(ic + ncell) = ip + pos(i1) + orb(i1)%nparm/4 + (ic - 1)*ntothyb
                                    end do
                                    eqpar(ip + poseq(i1))%neq = 2*ncell
                                elseif (.not. dosym) then
                                    do ic = 1, ncell
                                        eqpar(ip + poseq(i1))%idx(ic) = ip + pos(i1) + (ic - 1)*ntothyb
                                    end do
                                    eqpar(ip + poseq(i1))%neq = ncell
                                end if
                            end if
                        end if
                    end do
                end if
                !
            end if
            !
        end do
        !

        ! Put equal Zeta equal by symmetry
        do i1 = 1, ntotpar
            i2 = 1
            found = .false.
            !
            if (eqpar(i1)%zpar .and. eqpar(i1)%neq /= 0) then
                do while (i2 .lt. i1 .and. .not. found)
                    if (eqpar(i2)%zpar .and. eqpar(i2)%neq /= 0 .and. eqpar(i1)%kion == eqpar(i2)%kion&
                            &.and. eqpar(i1)%value == eqpar(i2)%value&
                            &.and. eqpar(i1)%mult == eqpar(i2)%mult) then
                        found = .true.
                        neqsign = isign(1, eqpar(i2)%neq)
                        tneq = iabs(eqpar(i2)%neq)
                        tmp(1:tneq) = eqpar(i2)%idx(1:tneq)
                        deallocate (eqpar(i2)%idx)
                        eqpar(i1)%neq = iabs(eqpar(i1)%neq)
                        eqpar(i2)%neq = tneq + eqpar(i1)%neq
                        allocate (eqpar(i2)%idx(eqpar(i2)%neq))
                        eqpar(i2)%idx(1:tneq) = tmp(1:tneq)
                        eqpar(i2)%idx(tneq + 1:tneq + eqpar(i1)%neq) = neqsign*eqpar(i1)%idx(1:eqpar(i1)%neq)
                        eqpar(i1)%neq = 0
                    end if
                    i2 = i2 + 1
                end do
            end if
        end do

        ! Count Parameters and put minus signs

        neq = 0
        do i2 = 1, orb(1)%nparm
            if (eqpar(poseq(1) + i2)%neq .ne. 0) then
                neq = neq + 1
            end if
        end do

        if (orb(1)%nparm .gt. 1) then
            eqpar(poseq(1) + orb(1)%nparm/2 + 1)%neq = -iabs(eqpar(poseq(1) + orb(1)%nparm/2 + 1)%neq)
            if (ipf .eq. 2 .and. orb(1)%ioptorb .eq. 900000) &
                    & eqpar(poseq(1) + (orb(1)%nparm*3)/4 + 1)%neq = -iabs(eqpar(poseq(1) + (orb(1)%nparm*3)/4 + 1)%neq)
            if (onlycontr .or. (yes9 .and. orb(1)%ioptorb .ne. 900000)) then
                do i3 = orb(1)%nparm/2 + 2, orb(1)%nparm
                    eqpar(poseq(1) + i3)%neq = -iabs(eqpar(poseq(1) + i3)%neq)
                end do
            end if
        end if

        do i1 = 2, norb
            if (orb(i1)%shell .ne. orb(i1 - 1)%shell) then
                do i2 = 1, orb(i1)%nparm
                    if (eqpar(poseq(i1) + i2)%neq .ne. 0) then
                        neq = neq + 1
                    end if
                end do
                if (orb(i1)%nparm .gt. 1) then
                    eqpar(poseq(i1) + orb(i1)%nparm/2 + 1)%neq = -iabs(eqpar(poseq(i1) + orb(i1)%nparm/2 + 1)%neq)
                    if (ipf .eq. 2 .and. orb(i1)%ioptorb .eq. 900000) then
                        eqpar(poseq(i1) + (orb(i1)%nparm*3)/4 + 1)%neq &
                            = -iabs(eqpar(poseq(i1) + (orb(i1)%nparm*3)/4 + 1)%neq)
                    end if
                end if
                if (onlycontr .or. (yes9 .and. orb(i1)%ioptorb .ne. 900000)) then
                    do i3 = orb(i1)%nparm/2 + 2, orb(i1)%nparm
                        eqpar(poseq(i1) + i3)%neq = -iabs(eqpar(poseq(i1) + i3)%neq)
                    end do
                end if
            end if
        end do

        deallocate (pos, poseq, tmp)
    end subroutine par_symm

    subroutine defz(n, alpha, beta, z)
        implicit none
        integer, intent(in) :: n
        double precision, intent(out) :: z(n)
        double precision, intent(in) :: alpha, beta
        integer i

        do i = 1, n
            z(i) = alpha*beta**(i - 1)
        end do
    end subroutine defz

end module mod_orbital
