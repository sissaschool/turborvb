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

subroutine write_output
    use allio
    implicit none
    real*8 cclock
    integer i
    logical tst
    write (6, *) ' All files written correctly ...'

    if (kaverage) then
        sumdiff = sumdiff/dble((nw/nk))/dble(i_main - iend)
    else
        sumdiff = sumdiff/dble(nw)/dble(i_main - iend)
    end if
    if (rata .gt. 0) rata = nacc/rata
    if (signflip .ne. 0. .and. tbra .ne. 0.) then
        signflip = signflip/(dble(i_main - iend)*dble(nw)*tbra)
        write (6, *) ' Sign flips per unit time =', signflip
    end if

    !       time=cclock()-time
    write (6, *) ' #####   TurboRVB  profiling (sec.) #####   '
    write (6, *) ' Time initialization =', timeinit
    write (6, *) ' Total time with no initialization =', time
    write (6, *) ' Total time with no measures ', timemc
    write (6, *) ' Time measures =', time_meas
    write (6, *) ' Time main =', time_main
    if (itestr .eq. -5) then
        write (6, *) ' Time optimization part =', timeopt
    end if
    if (iesbra) write (6, *) ' Time branching  =', time_branch
    write (6, *) ' Tracing the qmc update  move '
    write (6, *) ' Time ratiovar=', time_ratiovar
    write (6, *) ' Time uptabtot=', time_uptabtot
    write (6, *) ' Tracing the main routines '
    write (6, *) ' Time uptabpip  in uptabtot=', timepip
    write (6, *) ' Time upnewwf in ratiovar/uptabtot=', timewf
    write (6, *) ' Time upscratch  =', timescra
    write (6, *) ' accept. rate off diagonal moves =', rata
    if (itest .eq. 2) then
        rata = timemc/(time_meas + timeopt)
        write (6, *) ' Optimal nbra suggested =', max(1, nint(nbra/rata))
    end if

    write (6, *) ' Average time for 1000 generations ', time/(i_main - iend)*1000.d0
    psiav = psiav/counttot
    psisav = psisav/counttot
    countreg = countreg/counttot

    if (psisav .gt. psiav**2) then
        psisav = dsqrt(psisav - psiav**2)
    else
        psisav = 0.d0
    end if

    if (hopfraction .ne. 0.d0 .and. abs(parcut) .ne. 100.d0 .and. nacc .gt. 0) then
        write (6, *) ' Fraction large hopping moves ='                    &
                &, (nion - 1)*acclarge/nacc
    end if
    if (epscuttype .eq. -100 .and. nacc .gt. 0) then
        write (6, *) ' Fraction large lattice a  moves ='                  &
                &, acclarge/nacc
    end if
    if (npsa .gt. 0 .and. itest .eq. 1) then
        write (6, *) '# pseudo off diag moves per generation per walker =', &
                &naccpseudo/(dble(i_main - iend)*dble(nw))
        inquire (8, opened=tst)
        if (tst) close (8)
    end if
    if (epscuttype .gt. 0) then
        write (6, *) '# non trivial accepted/tried  =', &
                &nontr/(dble(i_main - iend)*dble(nw)*dble(nbra))
    end if
    if (epstldmc .ne. 0.d0 .and. nacc .gt. 0) then
        write (6, *) ' Number of rejected moves inside dmc cutoff ', nint(countcut)
        write (6, *) ' Corresponding to a fraction of ', countcut/nacc
    end if

    if (iesbra) write (6, *) 'Av. num.  of survived walkers/ # walkers&
            & in the branching', sumdiff
    if (itest .eq. 2) write (6, *) &
            &'Average inverse A  wf =', psiav, '+/-', psisav
    if (itest .ne. 2) then
        if (mod(typereg, 2) .eq. 0) then
            write (6, *) ' Average log det =', psiav, '+/-', psisav
        else
            write (6, *) ' Average log Psi  =', psiav, '+/-', psisav
        end if
        if (yesnleft) then
            write (6, *) ' Final value of av.energy  =', -lambda*ris(2)
        end if
        write (6, *) ' Average time x branching =', tave_cyrus/tcount_cyrus/ris(2)

    end if
    if (change_parr .and. itestr .eq. -5) then
        write (6, *) ' Number of times with parr <=', parr, ' amounts to ', iesconv
        write (6, *) ' Warning you can safely average with the previous', iesconv, ' Optimization steps'
    end if

    if (parcutg .ne. 0) then
        write (6, *) ' Fraction nodal surface/vol = ', countreg
    end if
    if (change_tstep .and. itest .eq. 2) then
        write (6, *) ' Final tstep found ', tstep
    end if
#ifdef _TIME
    if (yes_crystalj) then
        write (6, *) ' Fraction zero wfs (det & J)  =  ', count_zerowf/count_allwf
    else
        write (6, *) ' Fraction zero wfs (det only)  =  ', count_zerowf/count_allwf
    end if
    do i = 1, 11
        write (6, *) ' Time task direct/adjoint # ', i, timings(i), timingsb(i)
    end do
    write (6, *) ' Total  time  compute_fast/compute_fast_b ', sum(timings(1:11)), sum(timingsb(1:11))
#endif
end subroutine write_output
