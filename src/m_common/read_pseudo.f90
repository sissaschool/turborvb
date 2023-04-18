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

subroutine read_pseudo
    use allio
    implicit none
    integer :: i, j

#ifdef PARALLEL
    include 'mpif.h'
#endif

    if (npsa .le. 0) then

        nparshellmax = 0
        npseudopar = 1
        lmax = 1

        allocate (nparpshell(lmax, 1), jpseudo(lmax, 1), kindion(1)     &
                &, pshell(1), parshell(3, 1), rcutoff(1))
        allocate (wpseudo(2*lmax), legendre(lmax - 1, 1), versor(3, 1)    &
                &, wintpseudo(1))

        !            Initialize to zero
        nparpshell = 0
        jpseudo = 0
        kindion = 0
        pshell = 0
        parshell = 0
        rcutoff = 0
        wpseudo = 0
        legendre = 0
        versor = 0.d0
        wintpseudo = 0.d0

        if (npow .ne. 0.d0) then
            call error("main", " alpha .ne. 0 only with pseudo !! ", 1, rank)
        end if

    else
        ! pseudopotential calculation
        ! open file containing pseudo data
        ! npsa is the number of pseudo ion
        if (rank .eq. 0) then
            open (8, file=trim(pseudofile), status='unknown', form='formatted')
            read (8, *, iostat=iflagerr) pseudoname
        end if
        call checkiflagerr(iflagerr, rank, 'ERROR reading pseudo name')
#ifdef PARALLEL
        call mpi_bcast(pseudoname, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
        softcusp = .false.
        if (pseudoname .eq. 'ECP' .or. pseudoname .eq. 'ECS') then

            if (pseudoname .eq. 'ECS') then
                softcusp = .true.
                if (rank .eq. 0) write (6, *) ' Warning use of soft cusp condition, new implementation '
            else
                softcusp = .false.
            end if

            if (rank .eq. 0) then
                write (6, *) '#############################################'
                write (6, *) '    EFFECTIVE CORE POTENTIAL CALCULATION     '
                write (6, *) '#############################################'
            end if

            if ((itestrfn .eq. 6 .or. itestrfn .eq. 1) .and. npow .eq. 0.d0 .and. rank .eq. 0) then
                write (*, *) '    DMC FIXED NODE WITH LOCALITY APPROXIMATION   '
                if (itestrfn .eq. 1) write (6, *) ' Umrigar 93 algorithm '
                if (itestrfn .eq. 6) write (6, *) ' Improved local move integration '

            elseif ((itestrfn .eq. -2 .or. itestrfn .eq. -3) .and. npow .eq. 0.d0 .and. rank .eq. 0) then
                write (*, *) '    DMC FIXED NODE WITH OFF DIAGONAL PSEUDO    '
                if (itestrfn .eq. -2) write (6, *) ' heat bath after all electron diffusion '
                if (itestrfn .eq. -3) write (6, *) ' heat bath after single electron diffusion '

            elseif (itestrfn .eq. -6 .and. npow&
                    &.eq. 0.d0 .and. rank .eq. 0) then
                write (*, *) '    LRDMC FIXED NODE WITH OFF DIAGONAL PSEUDO     '

            elseif (itestrfn .eq. -6 .or. itestrfn .eq. -2&
                    &.or. itestrfn .eq. -3 .and. npow .gt. 0.d0) then
                if (itestrfn .eq. -6 .and. rank .eq. 0) then
                    write (*, *) ' LRDMC INTERPOLATION BETWEEN LOCALITY AND OFF DIAG '
                    write (*, *) ' Interpolating parameter', npow
                elseif (((itestrfn .eq. -2) .or. itestrfn .eq. -3)&
                        &.and. rank .eq. 0) then
                    write (*, *) ' DMC INTERPOLATION BETWEEN LOCALITY AND OFF DIAG '
                    write (*, *) ' Interpolating parameter', npow
                    if (itestrfn .eq. -2) write (6, *) ' heat bath after all electron diffusion '
                    if (itestrfn .eq. -3) write (6, *) ' heat bath after single electron diffusion '
                end if
                if (gamma .gt. 1.d0/npow - 1.d0) then
                    write (errmsg, *) 'gamma violates the sign! choose gamma <', 1.d0/npow - 1.d0
                    call checkiflagerr(1, rank, errmsg)
                end if
            end if
            nintpseudo = nintpsa

            if (rank .eq. 0) write (*, *) ' ************ nintpseudo read '      &
                    &, nintpseudo

            allocate (kindion(npsa + 1), pshell(npsa), rcutoff(npsa))
            allocate (versor(3, nintpseudo), wintpseudo(nintpseudo))

            kindion = 0
            pshell = 0
            rcutoff = 0.d0
            versor = 0.d0
            wintpseudo = 0.d0

            npseudopar = 1

            lmax = 1

            do ion = 1, npsa

                if (rank .eq. 0) then
                    read (8, *, iostat=iflagerr) kindion(ion), rcutoff(ion), pshell(ion)
                    iflagerrall = iflagerrall + iflagerr
                end if
                call checkiflagerr(iflagerr, rank, ' ERROR reading pseudo kindion')
#ifdef PARALLEL
                call mpi_bcast(kindion(ion), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(pshell(ion), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call mpi_bcast(rcutoff(ion), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

                allocate (ipsip(pshell(ion)))
                if (rank .eq. 0) then
                    read (8, *, iostat=iflagerr) (ipsip(j), j=1, pshell(ion))
                    iflagerrall = iflagerrall + iflagerr
                end if
                call checkiflagerr(iflagerr, rank, ' ERROR reading pseudo ipsip')

#ifdef PARALLEL
                call bcast_integer(ipsip, pshell(ion), 0, MPI_COMM_WORLD)
#endif

                ! kindion: ion type the ion-th pseudo refers to
                ! pshell: # of pseudo shell (L+1)
                ! nparpshell: # of gaussian per shell
                ! read for each ion all those parameters
                !            if(pshell(ion).gt.lmax) then
                !          if(rank.eq.0) write(6,*) '# of pseudo shell .gt.',lmax
                !ifdef PARALLEL
                !         call mpi_finalize(ierr)
                !endif
                !         stop
                !           endif

                if (pshell(ion) .gt. lmax) lmax = pshell(ion)

                i = 0
                do j = 1, pshell(ion)
                    !             jpseudo(j,ion)=npseudopar
                    npseudopar = npseudopar + ipsip(j)
                    i = i + ipsip(j)
                end do

                ! skip the pseudo parameters before their allocation
                if (rank .eq. 0) then
                    do j = 1, i
                        iflagerr = 1
                        read (8, *, err=103, end=103)
                        iflagerr = 0
                    end do
103                 continue
                    iflagerrall = iflagerrall + iflagerr
                end if
                call checkiflagerr(iflagerr, rank, ' ERROR skipping pseudo parameters')

                deallocate (ipsip)
            end do !npsa

            ! allocate the correct number of pseudo  parameters
            npseudopar = npseudopar - 1

            allocate (wpseudo(2*lmax), legendre(lmax - 1, nintpseudo), &
                    &   nparpshell(lmax, npsa), jpseudo(lmax, npsa))

            wpseudo = 0.d0
            legendre = 0.d0
            nparpshell = 0
            jpseudo = 0

            ! define versor and legendre
            iflagerr = 0
            call definition(nintpseudo, wintpseudo, lmax, versor, legendre, rank, iflagerr)
            call checkiflagerr(iflagerr, rank, ' ERROR checking pseudo definition')

            allocate (parshell(3, npseudopar))

            if (rank .eq. 0) then
                write (6, *) '------  parameters for pseudopotentials -------'
                write (6, *) 'Max angular momentum pseudo ', lmax
                write (6, *) '# of quadrature points in the projector', nintpsa
                write (6, *) '# of pseudo atoms', npsa
                write (6, *) '# of gaussian for pseudo', npseudopar
                write (6, *) '-----------------------------------------------'
            end if

            ! now read the pseudo parameters (parshell)
            npseudoparn = 1
            if (rank .eq. 0) then
                rewind (8)
                iflagerr = 1
                read (8, *, err=104, end=104)
                iflagerr = 0
104             continue
                iflagerrall = iflagerrall + iflagerr
            end if
            call checkiflagerr(iflagerr, rank, ' ERROR reading pseudo name')

            do ion = 1, npsa
                if (rank .eq. 0) then
                    iflagerr = 1
                    read (8, *, err=105, end=105)
                    read (8, *, err=105, end=105) (nparpshell(j, ion), j=1, pshell(ion))
                    iflagerr = 0
105                 continue
                    iflagerrall = iflagerrall + iflagerr
                end if
                call checkiflagerr(iflagerr, rank, ' ERROR reading pseudo kindion')
#ifdef PARALLEL
                call bcast_integer(nparpshell(1, ion), pshell(ion), 0, MPI_COMM_WORLD)
#endif

                do j = 1, pshell(ion)
                    jpseudo(j, ion) = npseudoparn
                    npseudoparn = npseudoparn + nparpshell(j, ion)
                end do

                indmax = jpseudo(pshell(ion), ion) + nparpshell(pshell(ion), ion) - 1
                do i = jpseudo(1, ion), indmax
                    if (rank .eq. 0) then
                        iflagerr = 1
                        read (8, *, end=106, err=106) (parshell(j, i), j=1, 3)
                        iflagerr = 0
106                     continue
                        iflagerrall = iflagerrall + iflagerr
                    end if
                    call checkiflagerr(iflagerr, rank, ' ERROR reading pseudo parameters')
#ifdef PARALLEL
                    call mpi_bcast(parshell(1, i), 3, MPI_DOUBLE_PRECISION, 0            &
                         &, MPI_COMM_WORLD, ierr)
#endif
                end do
            end do

        end if ! endif for pseudoname

        nparshellmax = 0
        do ion = 1, npsa
            do j = 1, pshell(ion)
                if (nparpshell(j, ion) .gt. nparshellmax)&
                        & nparshellmax = nparpshell(j, ion)
            end do
        end do
    end if ! endif for  npsa

    call checkiflagerr(iflagerrall, rank, ' at least one ERROR in reading pseudo')
end subroutine read_pseudo
