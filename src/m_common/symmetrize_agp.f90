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

subroutine symmetrizeagp(nnozero_c, nozero_c, jbradet, jbradetn, dsw&
        &, iessw0, scale_c, itouch, detmat_c, nelorb_c, nelorb_at, nelcol_c&
        &, symmagp, yes_hermite)
    use constants, only: ipc, ipf, deps
    use cell, only: cellscale, phase2pi
    use allio, only: rank, rion, yes_crystal, real_agp, molyes, pfaffup, kiontot
    implicit none

    ! argument variables
    integer, intent(in) :: nnozero_c, iessw0, nelorb_c, nelcol_c, nelorb_at, &
                        &  nozero_c(*), jbradet(*), jbradetn(*)
    integer, intent(inout) :: itouch(*)
    real*8, intent(inout) :: dsw(*)
    real*8, intent(inout) :: detmat_c(ipc*nelorb_c*nelcol_c), scale_c(*)
    logical, intent(in) :: symmagp, yes_hermite

    ! local variables
    integer ii, jj, kk, ind, indc, iessw, iesswt, ix, iy, iyr, ierr, ndim
    real*8 dsw0, dswc(2), riondiff(3), cosphihalf, sinphihalf, cost, detr, deti

    !     Input detmat_c output detmat_c symmetrized
    do ii = 1, nnozero_c
        scale_c(ipc*(ii - 1) + 1:ipc*ii) = detmat_c(ipc*(nozero_c(ii) - 1) + 1:ipc*nozero_c(ii))

    end do
    if (ipc .eq. 2 .and. real_agp .and. .not. molyes) then
        do ii = 1, nnozero_c
            scale_c(2*ii) = 0.d0
        end do
    end if

    !     Remove not allowed matrix elements
    detmat_c = 0.d0
    do ii = 1, nnozero_c
        detmat_c(ipc*(nozero_c(ii) - 1) + 1:ipc*nozero_c(ii)) = scale_c(ipc*(ii - 1) + 1:ipc*ii)
    end do
    iessw = 0
    ndim = nelorb_at/2

    if (ipc .eq. 1) then
        do ii = 1, nnozero_c
            jj = abs(jbradet(ii))
            if (jj .gt. iessw) iessw = jj
            if (jj .ne. 0) then
                dsw(jj) = 0.d0
                itouch(jj) = 0
            end if
        end do
        do ii = 1, nnozero_c
            if (jbradet(ii) .gt. 0) then
                dsw(jbradet(ii)) = dsw(jbradet(ii)) + scale_c(ii)
                itouch(jbradet(ii)) = itouch(jbradet(ii)) + 1
            elseif (jbradet(ii) .lt. 0) then
                dsw(-jbradet(ii)) = dsw(-jbradet(ii)) - scale_c(ii)
                itouch(-jbradet(ii)) = itouch(-jbradet(ii)) + 1
            end if
        end do
    else
        do ii = 1, nnozero_c
            jj = abs(jbradet(ii))
            if (jj .ne. 0) then
                dsw(2*jj - 1:2*jj) = 0.d0
                itouch(jj) = 0
            end if
            if (jj .gt. iessw) iessw = jj
        end do
        do ii = 1, nnozero_c
            iy = (nozero_c(ii) - 1)/nelorb_c + 1
            ix = nozero_c(ii) - (iy - 1)*nelorb_c
            ind = abs(jbradet(ii))

            if (jbradet(ii) .gt. 0) then
                dsw(2*jbradet(ii) - 1) = dsw(2*jbradet(ii) - 1) + scale_c(2*ii - 1)
                dsw(2*jbradet(ii)) = dsw(2*jbradet(ii)) + scale_c(2*ii)
                itouch(jbradet(ii)) = itouch(jbradet(ii)) + 1
            elseif (jbradet(ii) .lt. 0) then
                dsw(-2*jbradet(ii) - 1) = dsw(-2*jbradet(ii) - 1) - scale_c(2*ii - 1)
                dsw(-2*jbradet(ii)) = dsw(-2*jbradet(ii)) - scale_c(2*ii)
                itouch(-jbradet(ii)) = itouch(-jbradet(ii)) + 1
            end if
        end do
    end if

    if (ipc .eq. 1) then

        ind = 1
        do ii = 1, iessw0
            kk = -jbradetn(ind)
            dsw0 = 0.d0
            do jj = 1, kk
                ix = jbradetn(ind + 1)
                iy = jbradetn(ind + 2)
                if (ix .gt. 0 .and. iy .gt. 0) then
                    dsw0 = dsw0 + detmat_c(nelorb_c*(iy - 1) + ix)
                elseif (abs(ix) .gt. 0 .and. abs(iy) .gt. 0) then
                    ix = abs(ix)
                    iy = abs(iy)
                    dsw0 = dsw0 - detmat_c(nelorb_c*(iy - 1) + ix)
                end if
                ind = ind + 2
            end do
            !      replace detmat_c symmetrized
            dsw0 = dsw0/kk
            ind = ind - 2*kk
            do jj = 1, kk
                ix = jbradetn(ind + 1)
                iy = jbradetn(ind + 2)
                if (ix .gt. 0 .and. iy .gt. 0) then
                    detmat_c(nelorb_c*(iy - 1) + ix) = dsw0
                    if (ipf .eq. 2 .and. iy .le. nelorb_c) then
                        detmat_c(nelorb_c*(ix - 1) + iy) = -dsw0
                        if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                            if (ix .le. ndim .and. iy .le. ndim .and. .not. pfaffup) then
                                detmat_c(nelorb_c*(iy + ndim - 1) + ix + ndim) = dsw0
                                detmat_c(nelorb_c*(ix + ndim - 1) + iy + ndim) = -dsw0
                            elseif (iy .gt. ndim .and. ix .le. ndim) then
                                detmat_c(nelorb_c*(ix + ndim - 1) + iy - ndim) = dsw0
                                detmat_c(nelorb_c*(iy - ndim - 1) + ix + ndim) = -dsw0
                            end if
                        end if
                    elseif (ipf .eq. 1) then
                        if (symmagp .and. iy .lt. nelorb_c) detmat_c(nelorb_c*(ix - 1) + iy) = dsw0
                    end if
                elseif (abs(ix) .gt. 0 .and. abs(iy) .gt. 0) then
                    ix = abs(ix)
                    iy = abs(iy)
                    detmat_c(nelorb_c*(iy - 1) + ix) = -dsw0
                    if (ipf .eq. 2 .and. iy .le. nelorb_c) then
                        detmat_c(nelorb_c*(ix - 1) + iy) = dsw0
                        if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                            if (ix .le. ndim .and. iy .le. ndim .and. .not. pfaffup) then
                                detmat_c(nelorb_c*(iy + ndim - 1) + ix + ndim) = -dsw0
                                detmat_c(nelorb_c*(ix + ndim - 1) + iy + ndim) = dsw0
                            elseif (iy .gt. ndim .and. ix .le. ndim) then
                                detmat_c(nelorb_c*(iy - ndim - 1) + ix + ndim) = dsw0
                                detmat_c(nelorb_c*(ix + ndim - 1) + iy - ndim) = -dsw0
                            end if
                        end if
                    elseif (ipf .eq. 1) then
                        if (symmagp .and. iy .lt. nelorb_c) detmat_c(nelorb_c*(ix - 1) + iy) = -dsw0
                    end if
                end if
                ind = ind + 2
            end do
            ind = ind + 1
        end do

    else

        ind = 1
        do ii = 1, iessw0
            kk = -jbradetn(ind)
            dswc = 0.d0
            do jj = 1, kk
                ix = jbradetn(ind + 1)
                iy = jbradetn(ind + 2)
                if (ix .gt. 0) then
                    iyr = iy
                    iy = abs(iy)
                    indc = nelorb_c*(iy - 1) + ix
                    dswc(1) = dswc(1) + detmat_c(2*indc - 1)
                    dswc(2) = dswc(2) + detmat_c(2*indc)
                elseif (abs(ix) .gt. 0 .and. abs(iy) .gt. 0) then
                    ix = abs(ix)
                    iyr = iy
                    iy = abs(iy)
                    indc = nelorb_c*(iy - 1) + ix
                    dswc(1) = dswc(1) - detmat_c(2*indc - 1)
                    dswc(2) = dswc(2) - detmat_c(2*indc)
                end if
                ind = ind + 2
            end do
            !      replace detmat_c symmetrized
            dswc = dswc/kk
            ind = ind - 2*kk
            do jj = 1, kk
                ix = jbradetn(ind + 1)
                iy = jbradetn(ind + 2)
                if (ix .gt. 0) then
                    iyr = iy
                    iy = abs(iy)
                    indc = nelorb_c*(iy - 1) + ix
                    detmat_c(2*indc - 1:2*indc) = dswc(:)
                    if (ipf .eq. 2 .and. iy .le. nelorb_c) then
                        detr = detmat_c(2*indc - 1)
                        deti = detmat_c(2*indc)
                        indc = nelorb_c*(ix - 1) + iy
                        detmat_c(2*indc - 1) = -detr
                        detmat_c(2*indc) = -deti
                        if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                            if (ix .le. ndim .and. iy .le. ndim .and. .not. pfaffup) then
                                detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim) - 1) = detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim)) = -deti
                                else
                                    detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim)) = deti
                                end if
                                detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim) - 1) = -detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim)) = deti
                                else
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim)) = -deti
                                end if
                            elseif (iy .gt. ndim .and. ix .le. ndim) then
                                detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim) - 1) = -detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim)) = deti
                                else
                                    detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim)) = -deti
                                end if
                                detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim) - 1) = detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim)) = -deti
                                else
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim)) = deti
                                end if
                            end if
                        end if
                    elseif (ipf .eq. 1) then
                        if (symmagp .and. ix .le. nelorb_at .and. iy .le. nelorb_at) then
                            detr = detmat_c(2*indc - 1)
                            deti = detmat_c(2*indc)
                            indc = nelorb_c*(ix - 1) + iy
                            detmat_c(2*indc - 1) = detr
                            if (yes_hermite .and. ix .eq. iy) then
                                detmat_c(2*indc) = 0.d0
                            elseif (yes_hermite) then
                                detmat_c(2*indc) = -deti
                            else
                                detmat_c(2*indc) = deti
                            end if
                        end if
                    end if
                elseif (abs(ix) .gt. 0 .and. abs(iy) .gt. 0) then
                    ix = abs(ix)
                    iyr = iy
                    iy = abs(iy)
                    indc = nelorb_c*(iy - 1) + ix
                    detmat_c(2*indc - 1) = -dswc(1)
                    detmat_c(2*indc) = -dswc(2)
                    if (ipf .eq. 2 .and. iy .le. nelorb_c) then
                        detr = detmat_c(2*indc - 1)
                        deti = detmat_c(2*indc)
                        indc = nelorb_c*(ix - 1) + iy
                        detmat_c(2*indc - 1) = -detr
                        detmat_c(2*indc) = -deti
                        if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                            if (ix .le. ndim .and. iy .le. ndim .and. .not. pfaffup) then
                                detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim) - 1) = detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim)) = -deti
                                else
                                    detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim)) = deti
                                end if
                                detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim) - 1) = -detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim)) = deti
                                else
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim)) = -deti
                                end if
                            elseif (iy .gt. ndim .and. ix .le. ndim) then
                                detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim) - 1) = -detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim)) = deti
                                else
                                    detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim)) = -deti
                                end if
                                detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim) - 1) = detr
                                if (yes_hermite) then
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim)) = -deti
                                else
                                    detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim)) = deti
                                end if
                            end if
                        end if
                    elseif (ipf .eq. 1) then
                        if (symmagp .and. ix .le. nelorb_at .and. iy .le. nelorb_at) then
                            detr = detmat_c(2*indc - 1)
                            deti = detmat_c(2*indc)
                            indc = nelorb_c*(ix - 1) + iy
                            if (yes_hermite .and. ix .eq. iy) then
                                detmat_c(2*indc) = 0.d0
                                detmat_c(2*indc - 1) = detr
                            elseif (yes_hermite) then
                                detmat_c(2*indc - 1) = detr
                                detmat_c(2*indc) = -deti
                            else
                                detmat_c(2*indc - 1) = detr
                                detmat_c(2*indc) = deti
                            end if
                        end if
                    end if
                end if
                ind = ind + 2
            end do
            ind = ind + 1
        end do
    end if

    if (ipc .eq. 1) then
        do ii = 1, iessw
            dsw(ii) = dsw(ii)/itouch(ii)
        end do
    else
        do ii = 1, iessw
            dsw(2*ii - 1:2*ii) = dsw(2*ii - 1:2*ii)/itouch(ii)
        end do
    end if

    if (ipc .eq. 1) then
        do ii = 1, nnozero_c
            iy = (nozero_c(ii) - 1)/nelorb_c + 1
            ix = nozero_c(ii) - (iy - 1)*nelorb_c
            if (jbradet(ii) .gt. 0) then
                detmat_c(nozero_c(ii)) = dsw(jbradet(ii))
            elseif (jbradet(ii) .lt. 0) then
                detmat_c(nozero_c(ii)) = -dsw(-jbradet(ii))
            end if

            if (ipf .eq. 2 .and. iy .le. nelorb_c) then
                indc = nelorb_c*(ix - 1) + iy
                detmat_c(indc) = -detmat_c(nozero_c(ii))
                if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                    if (ix .le. ndim .and. iy .le. ndim .and. .not. pfaffup) then
                        detmat_c(nelorb_c*(iy + ndim - 1) + ix + ndim) = detmat_c(nozero_c(ii))
                        detmat_c(nelorb_c*(ix + ndim - 1) + iy + ndim) = -detmat_c(nozero_c(ii))
                    elseif (iy .gt. ndim .and. ix .le. ndim) then
                        detmat_c(nelorb_c*(ix + ndim - 1) + iy - ndim) = detmat_c(nozero_c(ii))
                        detmat_c(nelorb_c*(iy - ndim - 1) + ix + ndim) = -detmat_c(nozero_c(ii))
                    end if
                end if
            elseif (ipf .eq. 1) then
                if (symmagp .and. iy .le. nelorb_c .and. ix .ne. iy) then
                    indc = nelorb_c*(ix - 1) + iy
                    detmat_c(indc) = detmat_c(nozero_c(ii))
                end if
            end if
        end do
    else
        do ii = 1, nnozero_c
            iy = (nozero_c(ii) - 1)/nelorb_c + 1
            ix = nozero_c(ii) - (iy - 1)*nelorb_c
            if (jbradet(ii) .gt. 0) then
                detmat_c(2*nozero_c(ii) - 1:2*nozero_c(ii)) = dsw(2*jbradet(ii) - 1:2*jbradet(ii))
            elseif (jbradet(ii) .lt. 0) then
                detmat_c(2*nozero_c(ii) - 1:2*nozero_c(ii)) = -dsw(-2*jbradet(ii) - 1:-2*jbradet(ii))
            end if
            if (ipf .eq. 2 .and. iy .le. nelorb_c) then
                detr = detmat_c(2*nozero_c(ii) - 1)
                deti = detmat_c(2*nozero_c(ii))
                indc = nelorb_c*(ix - 1) + iy
                detmat_c(2*indc - 1) = -detr
                detmat_c(2*indc) = -deti
                if (symmagp .and. kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
                    if (ix .le. ndim .and. iy .le. ndim .and. .not. pfaffup) then
                        detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim) - 1) = detr
                        if (yes_hermite) then
                            detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim)) = -deti
                        else
                            detmat_c(2*nelorb_c*(iy + ndim - 1) + 2*(ix + ndim)) = deti
                        end if
                        detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim) - 1) = -detr
                        if (yes_hermite) then
                            detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim)) = deti
                        else
                            detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy + ndim)) = -deti
                        end if
                    elseif (iy .gt. ndim .and. ix .le. ndim) then
                        detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim) - 1) = -detr
                        if (yes_hermite) then
                            detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim)) = deti
                        else
                            detmat_c(2*nelorb_c*(iy - ndim - 1) + 2*(ix + ndim)) = -deti
                        end if
                        detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim) - 1) = detr
                        if (yes_hermite) then
                            detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim)) = -deti
                        else
                            detmat_c(2*nelorb_c*(ix + ndim - 1) + 2*(iy - ndim)) = deti
                        end if
                    end if
                end if
            elseif (ipf .eq. 1) then
                !      The effective lambda is not hermitian.
                if (symmagp .and. ix .le. nelorb_at .and. iy .le. nelorb_at) then
                    indc = nelorb_c*(ix - 1) + iy
                    detmat_c(2*indc - 1) = detmat_c(2*nozero_c(ii) - 1)
                    if (yes_hermite .and. ix .eq. iy) then
                        detmat_c(2*indc) = 0.d0
                    elseif (yes_hermite) then
                        detmat_c(2*indc) = -detmat_c(2*nozero_c(ii))
                    else
                        detmat_c(2*indc) = detmat_c(2*nozero_c(ii))
                    end if
                end if
            end if
        end do
    end if
    return
end
