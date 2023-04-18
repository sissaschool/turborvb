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

subroutine energy(Lz, table, tabler, diag, wtot, veff, veffright, enerdiff&
        &, itest, npow, gamma, nel, istart, indtm, epscutdmc)
    use allio, only: typereg
    implicit none
    integer Lz, i, itest, index, nel, istart, ipart, indtm(*)
    real*8 tabler(*), table(*), diag, wtot, veff, gamma, npow, veffright&
            &, epscutdmc, costnpow, enerdiff

    !     wtot is the sum over x' of  | G_x',x | where it is assumed that
    !     G_x,x (the diagonal part) is given by the lattice FN aproximation:
    !     G_x,x= diag+ (1+gamma) veff
    !     Here |G_x',x| =table(i)          if   table > 0 and x' ne x
    !     Here |G_x',x| =|gamma table(i)|  if   table < 0 and x' ne x
    !        veff = sum_x' G_{x'.x} < 0 and x' ne x.

    enerdiff = 0.d0

!   if(itest.eq.1.or.itest.eq.3.or.itest.eq.-2.or.itest.eq.-3) then
    wtot = diag
    veff = 0.d0
    veffright = 0.d0

    if (epscutdmc .eq. 0.d0 .and. typereg .lt. 6) then
        do i = 1, Lz
            index = (i - 1)/nel + 1
            ipart = i - (index - 1)*nel
            wtot = wtot + table(i)
            if (index .ge. istart .and. index .le. indtm(ipart)) then
                if (table(i) .lt. 0.) then
                    veff = veff + table(i)
                else
                    veffright = veffright + table(i)
                end if
            else
                if (table(i) .lt. 0.) veff = veff + table(i)
            end if
        end do

    else
        costnpow = npow*(1.d0 + gamma)
        !  The Guiding function ME are in table. The Vsf is defined in term
        !  of a different wf with the same nodes, the corresponding ME are
        !  written in tabler
        do i = 1, Lz
            index = (i - 1)/nel + 1
            ipart = i - (index - 1)*nel
            if (index .ge. istart .and. index .le. indtm(ipart)) then
                if (table(i) .lt. 0) then
                    veff = veff + tabler(i)
                    enerdiff = enerdiff - (1.d0 + gamma)*(tabler(i) - table(i))
                    wtot = wtot - gamma*table(i) + (1.d0 + gamma)*tabler(i)
                else
                    veffright = veffright + tabler(i)
                    wtot = wtot + table(i)*(1.d0 - costnpow) + costnpow*tabler(i)
                    enerdiff = enerdiff - costnpow*(tabler(i) - table(i))
                end if
            else
                if (table(i) .lt. 0) then
                    enerdiff = enerdiff - (1.d0 + gamma)*(tabler(i) - table(i))
                    veff = veff + tabler(i)
                    wtot = wtot - gamma*table(i) + (1.d0 + gamma)*tabler(i)
                else
                    wtot = wtot + table(i)
                end if
            end if
        end do
    end if

!   else
!       veff = 0.d0
!       veffright = 0.d0
!       wtot = 0.d0
!       costnpow = npow * (1.d0 + gamma)
!       do i = 1, Lz
!           index = (i - 1) / nel + 1
!           if(npow.ne.0.d0.and.table(i).gt.0.d0.and.index.ge.istart) then
!               wtot = wtot + table(i) * (1.d0 - costnpow)
!           elseif(table(i).gt.0.d0) then
!               wtot = wtot + table(i)
!           else
!               wtot = wtot - table(i) * gamma
!           endif
!       enddo
!   endif

    return
end subroutine energy

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!

subroutine energy_complex(Lz, table, tabler, diag, wtot, veff, veffright &
                          , enerdiff, itest, npow, gamma, nel, istart, indtm, epscutdmc)

    use allio, only: typereg
    implicit none
    integer Lz, i, itest, index, nel, istart, ipart, indtm(*)
    real(8) tabler(*), diag, wtot, veff, gamma, npow, veffright &
        , epscutdmc, costnpow, enerdiff
    !
    complex(8) table(*)

    enerdiff = 0.d0

!   if(itest.eq.1.or.itest.eq.3.or.itest.eq.-2.or.itest.eq.-3) then
    ! for dmclrdmc
    wtot = diag
    veff = 0.d0
    veffright = 0.d0
    if (epscutdmc .eq. 0.d0 .and. typereg .lt. 6) then
        do i = 1, Lz
            index = (i - 1)/nel + 1
            ipart = i - (index - 1)*nel
            wtot = wtot + table(i)
            if (index .ge. istart .and. index .le. indtm(ipart)) then
                if (dreal(table(i)) .lt. 0.d0) then
                    veff = veff + dreal(table(i))
                else
                    veffright = veffright + dreal(table(i))
                end if
            else
                if (dreal(table(i)) .lt. 0.d0) veff = veff + dreal(table(i))
            end if
        end do
    else
        costnpow = npow*(1.d0 + gamma)
        do i = 1, Lz
            !  The Guiding function ME are in table. The Vsf is defined in term
            !  of a different wf with the same nodes, the corresponding ME are
            !  written in tabler
            index = (i - 1)/nel + 1
            ipart = i - (index - 1)*nel
            if (index .ge. istart .and. index .le. indtm(ipart)) then
                if (tabler(i) .lt. 0) then
                    veff = veff + tabler(i)
                    enerdiff = enerdiff - (1.d0 + gamma)*(tabler(i) - table(i))
                    wtot = wtot - gamma*table(i) + (1.d0 + gamma)*tabler(i)
                else
                    veffright = veffright + tabler(i)
                    enerdiff = enerdiff - costnpow*(tabler(i) - table(i))
                    wtot = wtot + table(i)*(1.d0 - costnpow) + costnpow*tabler(i)
                end if
            else
                if (tabler(i) .lt. 0) then
                    enerdiff = enerdiff - (1.d0 + gamma)*(tabler(i) - table(i))
                    veff = veff + tabler(i)
                    wtot = wtot - gamma*table(i) + (1.d0 + gamma)*tabler(i)
                else
                    wtot = wtot + table(i)
                end if
            end if
        end do
    end if
!   else ! itest
!       ! vmc case
!       veff = 0.d0
!       veffright = 0.d0
!       wtot = 0.d0
!       costnpow = npow * (1.d0 + gamma)
!       do i = 1, Lz
!           index = (i - 1) / nel + 1
!           if(npow.ne.0.d0.and.dreal(table(i)).gt.0.d0.and.index.ge.istart) then
!               wtot = wtot + dreal(table(i)) * (1.d0 - costnpow)  ! weight wtot remains real also in the case of complex w.f.
!           elseif(dreal(table(i)).gt.0.d0) then
!               wtot = wtot + dreal(table(i))
!           else
!               wtot = wtot - dreal(table(i)) * gamma
!           endif
!       enddo
!   endif
    return

end subroutine energy_complex
