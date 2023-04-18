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

subroutine uptable_complex(nelup, neldo, indt, table, tabler, winvup &
                           , winvdo, tabpip, tmu, epscut, psiln, psidetln, costa, parcutg, istart, typereg)
    implicit none

    integer nelup, neldo, i, jel, indt, jel2, istart, parcutg, typereg
    real(8) tabpip(nelup + neldo, *) &
        , tabler(nelup + neldo, indt) &
        , tmu(nelup + neldo, *), epscut, psidetln, psiln, costa
    logical flagsign
    ! complex variables
    complex(8) winvup(nelup, *), winvdo(neldo, *), table(nelup + neldo, indt), ratio

    ! Table(nel,indt) contains off diagonal pseudopotential terms of the
    ! global w.f.

    ! Here the bosonic regularizing wfb is Sqrt( Det^2+exp(2 epscut))
    ! The guiding is the fermionic one wf = Sign Det x Sqrt( Det^2+exp(2 epscut))
    ! abs(table) are the ratio of the bosonic wf.
    ! costa on input=plat(2)/alat**2 is biased

    if (indt .le. 0) return

    if (epscut .eq. 0.d0 .and. typereg .lt. 6) then
        do i = 1, indt
            do jel = 1, nelup
                table(jel, i) = tmu(jel, i)*tabpip(jel, i)*winvup(jel, i)
            end do
            do jel = 1, neldo
                jel2 = jel + nelup
                table(jel2, i) = tmu(jel2, i)*tabpip(jel2, i)*winvdo(jel, i)
            end do
        end do
        tabler(:, :) = real(table(:, :))
    else
        if (mod(typereg, 2) .eq. 0) then
            do i = 1, indt
                do jel = 1, nelup
                    ratio = winvup(jel, i)
                    table(jel, i) = tmu(jel, i)*tabpip(jel, i)*ratio
                    if (typereg .eq. 6) then
                        if (i .ge. istart) then
                            tabler(jel, i) = tmu(jel, i)*ratio
                        else
                            tabler(jel, i) = table(jel, i)
                        end if
                    else
                        call ratiofn_psi(epscut, psidetln, abs(ratio), tabler(jel, i), typereg)
                    end if
                end do
                do jel = 1, neldo
                    jel2 = jel + nelup
                    ratio = winvdo(jel, i)
                    table(jel2, i) = tmu(jel2, i)*tabpip(jel2, i)*ratio
                    if (typereg .eq. 6) then
                        if (i .ge. istart) then
                            tabler(jel2, i) = tmu(jel2, i)*ratio
                        else
                            tabler(jel2, i) = table(jel2, i)
                        end if
                    else
                        call ratiofn_psi(epscut, psidetln, abs(ratio), tabler(jel2, i), typereg)
                    end if
                end do
            end do
        else
            do i = 1, indt
                do jel = 1, nelup
                    ratio = winvup(jel, i)*tabpip(jel, i)
                    table(jel, i) = tmu(jel, i)*ratio
                    call ratiofn_psi(epscut, psiln, abs(ratio), tabler(jel, i), typereg)
                end do
                do jel = 1, neldo
                    jel2 = jel + nelup
                    ratio = winvdo(jel, i)*tabpip(jel2, i)
                    table(jel2, i) = tmu(jel2, i)*ratio
                    call ratiofn_psi(epscut, psiln, abs(ratio), tabler(jel2, i), typereg)
                end do
            end do
        end if
        if (typereg .ne. 6) then
            ! Table remains untouched tabler is the wf with the guiding
            do i = 1, indt
                do jel = 1, nelup + neldo
                    if (tabler(jel, i) .ne. 0.d0) then
                        tabler(jel, i) = real(table(jel, i))/tabler(jel, i)
                    end if
                end do
            end do
!       else
!       tabler=table
        end if
    end if ! end if over epscut.eq.0
    return

end subroutine uptable_complex
