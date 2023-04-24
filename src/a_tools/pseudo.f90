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

program pseudo
    implicit none
    integer npsa, ion, i, j, npseudopar, npseudoparn, lmax, indmax
    real*8 r, rcut, fun, toll, pseudofun
    real(8), dimension(:), allocatable :: rcutoff, psip
    real(8), dimension(:, :), allocatable :: parshell
    integer, dimension(:, :), allocatable :: nparpshell, jpseudo
    integer, dimension(:), allocatable :: kindion, pshell, ipsip
    character(3) pseudoname
    !   AAA    Lines to be added just after all definitions of variables.
    character(100) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'traslafort10'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added
    open (8, file='pseudo.dat', status='unknown', form='formatted')
    open (10, file='pseudo.plot', status='unknown', form='formatted')
    npsa = 1
    read (8, *) pseudoname
    allocate (kindion(npsa + 1), pshell(npsa), rcutoff(npsa))
    write (6, *) ' Tollerance pseudo energy/pseudo '
    read (5, *) toll

    npseudopar = 1
    lmax = 1

    do ion = 1, npsa
        read (8, *) kindion(ion), rcutoff(ion), pshell(ion)
        allocate (ipsip(pshell(ion)))
        read (8, *) (ipsip(j), j=1, pshell(ion))

        if (pshell(ion) .gt. lmax) lmax = pshell(ion)

        i = 0
        do j = 1, pshell(ion)
            !              jpseudo(j,ion)=npseudopar
            npseudopar = npseudopar + ipsip(j)
            i = i + ipsip(j)
        end do

        do j = 1, i
            read (8, *)
        end do
        deallocate (ipsip)

    end do
    ! allocate the correct number of pseudo  parameters
    npseudopar = npseudopar - 1

    allocate (nparpshell(lmax, npsa), jpseudo(lmax, npsa))

    ! define versor and legendre
    allocate (parshell(3, npseudopar), psip(npseudopar))

    rewind (8)
    read (8, *)
    npseudoparn = 1
    do ion = 1, npsa
        read (8, *)
        read (8, *) (nparpshell(j, ion), j=1, pshell(ion))

        do j = 1, pshell(ion)
            jpseudo(j, ion) = npseudoparn
            npseudoparn = npseudoparn + nparpshell(j, ion)
        end do

        indmax = jpseudo(pshell(ion), ion) + nparpshell(pshell(ion), ion) - 1
        write (6, *) ' indmax found =', indmax
        write (6, *) ' last jpseudo  =', jpseudo(pshell(ion), ion)
        write (6, *) ' last  nparpshell  =', nparpshell(pshell(ion), ion)

        do i = jpseudo(1, ion), indmax
            read (8, *) (parshell(j, i), j=1, 3)
        end do
    end do

    rcut = 0.d0

    do j = 1, pshell(1)
        write (6, *) ' Angular component =', j
        do i = 0, 1000
            r = i*0.01d0
            fun = pseudofun(nparpshell(j, 1), r                               &
                    &, parshell(1, jpseudo(j, 1)), psip)

            if (abs(fun) .gt. toll .and. r .gt. rcut) rcut = r

        end do
    end do

    write (6, *) ' Suggested cut-off pseudo  =', rcut

    do i = 1, 1000
        r = i*0.01d0
        write (10, 100) r, (pseudofun(nparpshell(j, 1), r                  &
                &, parshell(1, jpseudo(j, 1)), psip), j=1, pshell(1))
    end do
100 format(100e17.8e3)
!    100    format(100f17.8)

    if (npsa .eq. 1) then
        rewind (8)
        write (8, '(3A)') 'ECP'
        npseudoparn = 1
        do ion = 1, npsa
            write (8, '(I6,1f11.6,I6)') 1, rcut, pshell(ion)
            write (8, '(1000I6)') (nparpshell(j, ion), j=1, pshell(ion))

            do j = 1, pshell(ion)
                jpseudo(j, ion) = npseudoparn
                npseudoparn = npseudoparn + nparpshell(j, ion)
            end do

            indmax = jpseudo(pshell(ion), ion) + nparpshell(pshell(ion), ion) - 1

            do i = jpseudo(1, ion), indmax
                j = nint(parshell(2, i))
                write (8, '(1f19.14,I6,1f19.14)') parshell(1, i), j, parshell(3, i)
            end do
        end do
    end if
    stop
end

function pseudofun(nmax, r, param, psip)
    !
    implicit none
    real*8 pseudofun, param(3, *), r, logr, r2, psip(*)
    integer nmax, i
    !
    if (r .lt. 1d-9) r = 1d-9
    r2 = r**2
    logr = dlog(r)

    do i = 1, nmax
        psip(i) = dexp(-param(3, i)*r2 + logr*param(2, i))
    end do
    !        call my_dexp(nmax,psip,psip)
    !

    !        do i=1,nmax
    !        psip(i)=param(1,i)*psip(i)
    !        enddo

    pseudofun = 0.d0
    do i = 1, nmax
        !        pseudofun=pseudofun+param(1,i)*r**param(2,i)
        !    &               *dexp(-param(3,i)*r2)
        pseudofun = pseudofun + psip(i)*param(1, i)
    end do
    pseudofun = pseudofun/r2
    !
    return
end
