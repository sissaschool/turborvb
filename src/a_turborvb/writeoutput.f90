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
    use logger_io, only: log_info, log_warning
    implicit none
    real*8 cclock
    integer i
    logical tst
    call log_info(' All files written correctly ...')

    if (kaverage) then
        sumdiff = sumdiff/dble((nw/nk))/dble(i_main - iend)
    else
        sumdiff = sumdiff/dble(nw)/dble(i_main - iend)
    end if
    if (rata .gt. 0) rata = nacc/rata
    if (signflip .ne. 0. .and. tbra .ne. 0.) then
        signflip = signflip/(dble(i_main - iend)*dble(nw)*tbra)
        call log_info(' Sign flips per unit time =', signflip)
    end if

    !       time=cclock()-time
    call log_info(' #####   TurboRVB  profiling (sec.) #####   ')
    call log_info(' Time initialization =', timeinit)
    call log_info(' Total time with no initialization =', time)
    call log_info(' Total time with no measures ', timemc)
    call log_info(' Time measures =', time_meas)
    call log_info(' Time main =', time_main)
    if (itestr .eq. -5) then
        call log_info(' Time optimization part =', timeopt)
    end if
    if (iesbra) call log_info(' Time branching  =', time_branch)
    call log_info(' Tracing the qmc update  move ')
    call log_info(' Time ratiovar=', time_ratiovar)
    call log_info(' Time uptabtot=', time_uptabtot)
    call log_info(' Tracing the main routines ')
    call log_info(' Time uptabpip  in uptabtot=', timepip)
    call log_info(' Time upnewwf in ratiovar/uptabtot=', timewf)
    call log_info(' Time upscratch  =', timescra)
    call log_info(' accept. rate off diagonal moves =', rata)
    if (itest .eq. 2) then
        rata = timemc/(time_meas + timeopt)
        call log_info(' Optimal nbra suggested =', max(1, nint(nbra/rata)))
    end if

    call log_info(' Average time for 1000 generations ', time/(i_main - iend)*1000.d0)
    psiav = psiav/counttot
    psisav = psisav/counttot
    countreg = countreg/counttot

    if (psisav .gt. psiav**2) then
        psisav = dsqrt(psisav - psiav**2)
    else
        psisav = 0.d0
    end if

    if (hopfraction .ne. 0.d0 .and. abs(parcut) .ne. 100.d0 .and. nacc .gt. 0) then
        call log_info(' Fraction large hopping moves =', (nion - 1)*acclarge/nacc)
    end if
    if (epscuttype .eq. -100 .and. nacc .gt. 0) then
        call log_info(' Fraction large lattice a  moves =', acclarge/nacc)
    end if
    if (npsa .gt. 0 .and. itest .eq. 1) then
        call log_info('# pseudo off diag moves per generation per walker =', naccpseudo/(dble(i_main - iend)*dble(nw)))
        inquire (8, opened=tst)
        if (tst) close (8)
    end if
    if (epscuttype .gt. 0) then
        call log_info('# non trivial accepted/tried  =', nontr/(dble(i_main - iend)*dble(nw)*dble(nbra)))
    end if
    if (epstldmc .ne. 0.d0 .and. nacc .gt. 0) then
        call log_info(' Number of rejected moves inside dmc cutoff ', nint(countcut))
        call log_info(' Corresponding to a fraction of ', countcut/nacc)
    end if

    if (iesbra) call log_info('Av. num.  of survived walkers/ # walkers in the branching', sumdiff)
    if (itest .eq. 2) call log_info('Average inverse A  wf =', psiav, '+/-', psisav)
    if (itest .ne. 2) then
        if (mod(typereg, 2) .eq. 0) then
            call log_info(' Average log det =', psiav, '+/-', psisav)
        else
            call log_info(' Average log Psi  =', psiav, '+/-', psisav)
        end if
        if (yesnleft) then
            call log_info(' Final value of av.energy  =', -lambda*ris(2))
        end if
        call log_info(' Average time x branching =', tave_cyrus/tcount_cyrus/ris(2))

    end if
    if (change_parr .and. itestr .eq. -5) then
        call log_info(' Number of times with parr <=', parr, ' amounts to ', iesconv)
        call log_warning(' Warning you can safely average with the previous', iesconv, ' Optimization steps')
    end if

    if (parcutg .ne. 0) then
        call log_info(' Fraction nodal surface/vol = ', countreg)
    end if
    if (change_tstep .and. itest .eq. 2) then
        call log_info(' Final tstep found ', tstep)
    end if
#ifdef _TIME
    if (yes_crystalj) then
        call log_info(' Fraction zero wfs (det & J)  =  ', count_zerowf/count_allwf)
    else
        call log_info(' Fraction zero wfs (det only)  =  ', count_zerowf/count_allwf)
    end if
    do i = 1, 11
        call log_info(' Time task direct/adjoint # ', i, timings(i), timingsb(i))
    end do
    call log_info(' Total  time  compute_fast/compute_fast_b ', sum(timings(1:11)), sum(timingsb(1:11)))
#endif
end subroutine write_output
