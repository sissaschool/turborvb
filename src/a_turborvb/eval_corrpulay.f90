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

subroutine eval_corrpulay(jel, k_ion, qver, rad, kel, keln, versor, indt, rion &
                          , nion, nelup, neldo, nelorb, nelorbh, nelorbj, nelorbjh &
                          , winv, winvj, ainvup, ainvdo, winvjbar, winvjbarsz &
                          , iesdr, vj, costz, costz3, LBox &
                          , iessz, add_pulay, psip, pulay, corrpulay)

    use Cell

    implicit none

    integer k_ion, jel, indt, i, j, d, qver, nel, nelup, neldo &
        , nelorb, nelorbj, nelorbh, nelorbjh, add_pulay, nion &
        , iesd, iesd1, iesdr, ispin, ispinjel, iesdr2iesd, iesdr1iesd

    real*8 versor(3, *), rad, scal, psip(nelorbjh, *) &
        , winv(nelorb, 0:indt + 4, *), winvj(nelorbj, 0:indt + 4, *) &
        , rion(3, *), corrpulay(3, *), kel(3, nelup + neldo, 0:*) &
        , ainvup(nelup, *), ainvdo(neldo, *) &
        , winvjbar(nelorbjh, *), winvjbarsz(nelorbjh, *) &
        , jsign, pulay(3, *), rc(3), hess(3), tmp, grad(3), LBox, vj(*) &
        , costz(*), costz3(*), keln(3, *), gradver(3)

    logical iessz, iesspin

    ! PURPOSE
    ! this subroutine computes the correction to the gradient at the mesh points
    ! due to the dipendence of the rotated electron on the nuclear position

    ! INPUT
    ! jel rotated electron
    ! qver versor index

    if (add_pulay .eq. 0) return

    if (add_pulay .ge. 1) then

        nel = nelup + neldo

        gradver(:) = rion(:, k_ion) - kel(:, jel, 0)
        if (Lbox .gt. 0.d0) call ApplyPBC(gradver, 1)
        gradver(:) = gradver(:)/rad

        ! define a vector for (rion(d,k_ion)-kel(d,jel,0))/rad !!!!

        iesd = iesdr2iesd(iesdr)
        iesd1 = iesdr1iesd(iesdr)
        !     if(iesdr.eq.-5) iesd=1
        !     if(iesdr.eq.-6) iesd=4
        !     if(iesdr.eq.-7) iesd=-1

        if (iesdr .eq. -7 .or. iesd .eq. 2 .or. iesd .lt. 0) then
            iesspin = .true.
        else
            iesspin = .false.
        end if

        if (jel .le. nelup) then
            ispinjel = 1
        else
            ispinjel = -1
        end if

        if (LBox .le. 0.d0) then

            ! one body
            if (iesdr .le. -5) then

                do i = 1, nion

                    rc(:) = costz(i)*(keln(:, 1) - rion(:, i))
                    call jastrowgrad(rc, vj, iesd1, tmp, grad, tmp, -1)

                    corrpulay(:, 1) = grad(:)*costz(i)*costz3(i)

                    scal = 0.d0
                    do d = 1, 3
                        scal = scal + corrpulay(d, 1)*versor(d, qver)
                    end do
                    do d = 1, 3
                        corrpulay(d, 2) = corrpulay(d, 1) + scal*gradver(d)
                    end do

                    pulay(:, k_ion) = pulay(:, k_ion) - corrpulay(:, 2)

                end do

            end if

            ! two body
            do i = 1, nel

                if (i .ne. jel) then

                    if (iesspin) then
                        if (i .le. nelup) then
                            ! parallel spin
                            ispin = ispinjel
                        else
                            ! antiparallel spin
                            ispin = -ispinjel
                        end if
                    else
                        ! spin independent 2b jastrow
                        ispin = -1
                    end if

                    rc(:) = keln(:, 1) - kel(:, i, 0)
                    call jastrowgrad(rc, vj, iesd, tmp, grad, tmp, ispin)

                    scal = 0.d0
                    do d = 1, 3
                        scal = scal + grad(d)*versor(d, qver)
                    end do
                    do d = 1, 3
                        corrpulay(d, 1) = grad(d) + scal*gradver(d)
                    end do

                    pulay(:, k_ion) = pulay(:, k_ion) + corrpulay(:, 1)

                end if ! i.ne.jel

            end do

        else
            ! PBC
            if (iesdr .le. -5) then

                do i = 1, nion

                    rc(:) = keln(:, 1) - rion(:, i)
                    if (yes_tilted) then
                        call jastrowgrad_pbc(rc, vj, iesd1, tmp, grad, tmp, -1, costz(i), hess)
                    else
                        call jastrowgrad_pbc_ortho(rc, vj, iesd1, tmp, grad, tmp, -1, costz(i), Lbox, hess)
                    end if

                    corrpulay(:, 1) = grad(:)*costz(i)*costz3(i)

                    scal = 0.d0
                    do d = 1, 3
                        scal = scal + corrpulay(d, 1)*versor(d, qver)
                    end do
                    do d = 1, 3
                        corrpulay(d, 2) = corrpulay(d, 1) + scal*gradver(d)
                    end do

                    pulay(:, k_ion) = pulay(:, k_ion) - corrpulay(:, 2)

                end do

            end if

            ! two body
            do i = 1, nel

                if (i .ne. jel) then

                    if (iesspin) then
                        if (i .le. nelup) then
                            ! parallel spin
                            ispin = ispinjel
                        else
                            ! antiparallel spin
                            ispin = -ispinjel
                        end if
                    else
                        ! spin independent 2b jastrow
                        ispin = -1
                    end if

                    rc(:) = keln(:, 1) - kel(:, i, 0)
                    if (yes_tilted) then
                        call jastrowgrad_pbc(rc, vj, iesd, tmp, grad, tmp, ispin, 1.d0, hess)
                    else
                        call jastrowgrad_pbc_ortho(rc, vj, iesd, tmp, grad, tmp, ispin, 1.d0, Lbox, hess)
                    end if

                    scal = 0.d0
                    do d = 1, 3
                        scal = scal + grad(d)*versor(d, qver)
                    end do
                    do d = 1, 3
                        corrpulay(d, 1) = grad(d) + scal*gradver(d)
                    end do

                    pulay(:, k_ion) = pulay(:, k_ion) + corrpulay(:, 1)

                end if

            end do

        end if

        ! prepare caculation for three body
        ! quite expensive part (use blas ?)
        do i = 1, nelorbjh
            psip(i, 1) = 0.d0
        end do
        do j = 1, nel
            do i = 1, nelorbjh
                psip(i, 1) = psip(i, 1) + winvjbar(i, j)
            end do
        end do
        do i = 1, nelorbjh
            psip(i, 1) = psip(i, 1) - winvjbar(i, jel)
        end do

        if (iessz) then

            if (jel .le. nelup) then
                jsign = 1.d0
            else
                jsign = -1.d0
            end if

            ! quite expensive part
            do i = 1, nelorbjh
                psip(i, 2) = 0.d0
            end do
            do j = 1, nelup
                do i = 1, nelorbjh
                    psip(i, 2) = psip(i, 2) + winvjbarsz(i, j)
                end do
            end do
            do j = nelup + 1, nel
                do i = 1, nelorbjh
                    psip(i, 2) = psip(i, 2) - winvjbarsz(i, j)
                end do
            end do
            do i = 1, nelorbjh
                psip(i, 2) = psip(i, 2) - jsign*winvjbarsz(i, jel)
            end do

        end if

        ! three body
        do i = 1, nelorbjh

            scal = 0.d0
            do d = 1, 3
                scal = scal + winvj(i, indt + d, jel)*versor(d, qver)
            end do
            do d = 1, 3
                corrpulay(d, i) = winvj(i, indt + d, jel) + scal*gradver(d)
            end do

            pulay(:, k_ion) = pulay(:, k_ion) + corrpulay(:, i)*psip(i, 1)

        end do

        ! spin dependent three body
        if (iessz) then

            do i = 1, nelorbjh
                pulay(:, k_ion) = pulay(:, k_ion) + jsign*corrpulay(:, i)*psip(i, 2)
            end do

        end if

        ! determinantal part
        if (add_pulay .ge. 2) then

            do i = 1, nelorbh

                scal = 0.d0
                do d = 1, 3
                    scal = scal + winv(i, indt + d, jel)*versor(d, qver)
                end do
                do d = 1, 3
                    corrpulay(d, i) = winv(i, indt + d, jel) + scal*gradver(d)
                end do

                if (jel .le. nelup) then
                    pulay(:, k_ion) = pulay(:, k_ion) + corrpulay(:, i)*ainvup(jel, i)
                else
                    pulay(:, k_ion) = pulay(:, k_ion) + corrpulay(:, i)*ainvdo(jel - nelup, i)
                end if

            end do

        end if

    end if

    !      do i=1,nion
    !         write(6,*) 'pulay',i,(pulay(d,i),d=1,3)
    !      enddo

    return
end subroutine eval_corrpulay
