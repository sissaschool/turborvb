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

subroutine pareff(npar, initpower, nparsw, initpowersw, nparinv     &
        &, initpowerinv, endinv, nlead, kp_ion, rpar, reduce, jas_invariant&
        &, adrlambda, nmax_ion, type_atom, allfit, orbps)
    implicit none
    integer npar, nparsw, initpowersw, ncg, kp_ion, i, j, initpower        &
            &, nparinv, initpowerinv, endinv, nlead, nfun, ii, icek, powerexp, ntpar
    real*8 rpar(*), reduce(nlead, *), funloc, jas_invariant(4, *)
    real(8), dimension(:), allocatable :: rmax, rmin, rmaxj, rminj, rind &
            &, rmaxinv, rmininv
    real*8, external :: damping
    integer powermin, powermax, nmax_ion, ifirst, adr_ion, k1, k2, kk
    integer adrlambda(2, *), type_atom(*)
    integer, external :: findindex
    integer, dimension(:), allocatable :: imap
    logical allfit
    logical orbps(*)
    !        rpar=0  untouched variational parameters
    !        rpar>0 distance   lambda_{i,j}  Jastrow
    !        rpar<0 distance lambda_{i,j}   Det
    !        prepare the vectors consecutevely on reduce
    !        special correlation for lambdas

    !        parametrizzazione costante per almeno un tipo di lambda
    if (initpowersw .eq. -1 .or. initpower .eq. -1 .or. initpowerinv .eq. -1) &
            &   allocate (rind(kp_ion), imap(kp_ion))

    !#########################################################
    !        parametrizzazione costante per il determinante
    if (nparsw .gt. 0 .and. initpowersw .eq. -1) then
        allocate (rmax(nparsw), rmin(nparsw))
        rmin = 0.d0
        rmax = 0.d0
        !        define rmax,rmin interval,in terms of rpar
        call preprminmax(kp_ion, nparsw, rpar, rmin, rmax, rind, imap, 0, endinv)
    end if

    !##########################################################
    !        parametrizzazione costante per il jastrow
    if (npar .gt. 0 .and. initpower .eq. -1) then
        allocate (rmaxj(npar), rminj(npar))
        rminj = 0.d0
        rmaxj = 0.d0
        !        define rmax,rmin interval,in terms of rpar
        call preprminmax(kp_ion, npar, rpar, rminj, rmaxj, rind, imap, 1, endinv)
    end if

    !############################################################
    !       parametrizzazione costante per il jastrow Sz
    if (nparinv .gt. 0 .and. initpowerinv .eq. -1) then
        allocate (rmaxinv(nparinv), rmininv(nparinv))
        rmininv = 0.d0
        rmaxinv = 0.d0
        !        define rmax,rmin interval,in terms of rpar
        call preprminmax(kp_ion, nparinv, rpar, rmininv, rmaxinv             &
                &, rind, imap, 2, endinv)
    end if
    !#############################################################

    !        npar=sum npar_type
    ncg = npar + nparsw + nparinv

    !        Initialize to zero reduce
    do i = 1, ncg
        do j = 1, kp_ion
            reduce(i, j) = 0.d0
        end do
    end do

    !#########################################################
    !Firstly the parametrization of the Spin Jastrow
    if (initpowerinv .gt. 0 .and. nparinv .ne. 0) then

        do i = ncg - nparinv + 1, ncg

            !         indice che conta i collettivi

            do j = 1, endinv
                if (rpar(j) .gt. 0.d0)                                          &
                        &     reduce(i, j) = rpar(j)**(-dble(i - ncg + nparinv - 1 + initpowerinv))

            end do
        end do
    elseif (initpowerinv .eq. 0 .and. nparinv .ne. 0) then

        do j = 1, endinv
            if (rpar(j) .gt. 0.d0) reduce(ncg - nparinv + 1, j) = dlog(rpar(j))
        end do

        do i = ncg - nparinv + 2, ncg
            do j = 1, endinv
                if (rpar(j) .gt. 0.d0) then

                    reduce(i, j) = rpar(j)**(-dble(i - ncg + nparinv - 1))

                end if
            end do
        end do
    elseif (nparinv .ne. 0 .and. initpowerinv .eq. -1) then
        do i = ncg - nparinv + 1, ncg
            nfun = i + nparinv - ncg
            do j = 1, endinv
                if (rpar(j) .gt. 0.d0)                                          &
                        & reduce(i, j) = funloc(nfun, rpar(j), rmaxinv, rmininv)
            end do
        end do
        deallocate (rmaxinv, rmininv)
    elseif (initpowerinv .lt. -1 .and. nparinv .ne. 0) then
        ! parametrizzazione
        if (allfit) then
            ntpar = nmax_ion*nmax_ion
        else
            ntpar = ((nmax_ion + 1)*nmax_ion)/2
        end if
        powermin = -2 - initpowerinv
        powermax = powermin + nparinv/ntpar/4 - 1
        !        do  i=ncg-nparinv+1,ncg,3
        do powerexp = powermin, powermax
            ifirst = ncg - nparinv + 1 + (powerexp - powermin)*4*ntpar
            do j = 1, endinv
                if (rpar(j) .gt. 0.d0) then
                    k1 = type_atom(adrlambda(1, j))
                    k2 = type_atom(adrlambda(2, j))
                    if (orbps(j)) then
                        if (k1 .lt. k2) then
                            kk = k2
                            k2 = k1
                            k1 = kk
                        end if
                    else
                        if (k1 .gt. k2) then
                            kk = k2
                            k2 = k1
                            k1 = kk
                        end if
                    end if
                    adr_ion = findindex(k1, k2, nmax_ion, allfit)
                    i = ifirst + 4*(adr_ion - 1)
                    !           powerexp=(i-ncg+nparinv-1)/3-2-initpowerinv    !
                    if (powerexp .eq. 0) then
                        reduce(i, j) = jas_invariant(1, j)*dlog(rpar(j))*damping(rpar(j), 0)
                        reduce(i + 1, j) = jas_invariant(2, j)/rpar(j)*damping(rpar(j), 1)
                        reduce(i + 2, j) = jas_invariant(3, j)/rpar(j)**2*damping(rpar(j), 2)
                        reduce(i + 3, j) = jas_invariant(4, j)/rpar(j)**2*damping(rpar(j), 2)
                    else
                        reduce(i, j) = jas_invariant(1, j)*rpar(j)**(-powerexp)&
                                      & *damping(rpar(j), powerexp)
                        reduce(i + 1, j) = jas_invariant(2, j)*rpar(j)**(-powerexp - 1)&
                                          & *damping(rpar(j), powerexp + 1)
                        reduce(i + 2, j) = jas_invariant(3, j)*rpar(j)**(-powerexp - 2)&
                                          & *damping(rpar(j), powerexp + 2)
                        reduce(i + 3, j) = jas_invariant(4, j)*rpar(j)**(-powerexp - 2)&
                                          & *damping(rpar(j), powerexp + 2)
                    end if
                end if
            end do
        end do

    end if

    if (initpower .gt. 0 .and. npar .ne. 0) then
        !        s-wave
        do i = ncg - npar - nparinv + 1, ncg - nparinv
            do j = endinv + 1, kp_ion
                if (rpar(j) .gt. 0.d0)                                          &
                        &     reduce(i, j) = rpar(j)**(-dble(i - ncg + npar + nparinv - 1 + initpower))
            end do
        end do

    elseif (initpower .lt. -1 .and. npar .ne. 0) then
        !        Parametrization  density-density Jastrow
        ! parametrizzazione
        if (allfit) then
            ntpar = nmax_ion*nmax_ion
        else
            ntpar = ((nmax_ion + 1)*nmax_ion)/2
        end if
        powermin = -2 - initpower
        powermax = powermin + npar/ntpar/4 - 1
        do powerexp = powermin, powermax
            ifirst = ncg - npar - nparinv + 1 + 4*ntpar*(powerexp - powermin)
            do j = endinv + 1, kp_ion
                if (rpar(j) .gt. 0.d0) then
                    !           powerexp=(i-ncg+npar+nparinv-1)/3-2-initpower    !
                    k1 = type_atom(adrlambda(1, j))
                    k2 = type_atom(adrlambda(2, j))
                    !     The p(ix)-s(iy)  matrix element connecting an heavier p with a
                    !     lighter s is put in the upper part of the max_ion x max_ion matrix
                    !     as well as  s(ix)-p(iy) connecting a lighter s with an hevier p.
                    !     In the lower part of the matrix we put a lighter p with an heavier  s

                    if (orbps(j)) then
                        !            p -s
                        if (k1 .lt. k2) then
                            !            k1 is lighter atom than k2
                            kk = k2
                            k2 = k1
                            k1 = kk
                        end if
                    else
                        !            s-p or anything else
                        if (k1 .gt. k2) then
                            !            k1 is heavier than k2
                            kk = k2
                            k2 = k1
                            k1 = kk
                        end if
                    end if
                    adr_ion = findindex(k1, k2, nmax_ion, allfit)
                    i = ifirst + 4*(adr_ion - 1)
                    if (powerexp .eq. 0) then
                        reduce(i, j) = jas_invariant(1, j)*dlog(rpar(j))*damping(rpar(j), 0)
                        reduce(i + 1, j) = jas_invariant(2, j)/rpar(j)*damping(rpar(j), 1)
                        reduce(i + 2, j) = jas_invariant(3, j)/rpar(j)**2*damping(rpar(j), 2)
                        reduce(i + 3, j) = jas_invariant(4, j)/rpar(j)**2*damping(rpar(j), 2)
                    else
                        reduce(i, j) = jas_invariant(1, j)*rpar(j)**(-powerexp)&
                                      & *damping(rpar(j), powerexp)
                        reduce(i + 1, j) = jas_invariant(2, j)*rpar(j)**(-powerexp - 1)&
                                          & *damping(rpar(j), powerexp + 1)
                        reduce(i + 2, j) = jas_invariant(3, j)*rpar(j)**(-powerexp - 2)&
                                          & *damping(rpar(j), powerexp + 2)
                        reduce(i + 3, j) = jas_invariant(4, j)*rpar(j)**(-powerexp - 2)&
                                          & *damping(rpar(j), powerexp + 2)
                    end if
                end if
            end do
        end do

    elseif (initpower .eq. 0 .and. npar .ne. 0) then
        do j = endinv + 1, kp_ion
            if (rpar(j) .gt. 0.d0) reduce(ncg - npar - nparinv + 1, j) = dlog(rpar(j))
        end do
        do i = ncg - npar - nparinv + 2, ncg - nparinv

            do j = endinv + 1, kp_ion
                if (rpar(j) .gt. 0.d0) reduce(i, j) = &
                        &    rpar(j)**(-dble(i - ncg + npar + nparinv - 1))
            end do
        end do

    elseif (npar .ne. 0 .and. initpower .eq. -1) then
        do i = ncg - nparinv - npar + 1, ncg - nparinv

            nfun = i + nparinv + npar - ncg
            do j = endinv + 1, kp_ion
                if (rpar(j) .gt. 0.d0)                                          &
                        & reduce(i, j) = funloc(nfun, rpar(j), rmaxj, rminj)
            end do
        end do
        deallocate (rmaxj, rminj)
    end if

    !         output reduce  -->
    !          le ultime npar direzioni

    !#######      End parametrization for the 3Body Jastrow        ##########
    !########################################################################

    if (initpowersw .gt. 0 .and. nparsw .ne. 0) then
        do i = ncg - nparsw - npar - nparinv + 1, ncg - npar - nparinv
            !        sostituire direzioni con una parametrizzazione prefissata
            !        a npar parametri
            !        direzioni d labda_{i,j}
            !        introdurre un modulo che definisce
            do j = 1, kp_ion
                if (rpar(j) .lt. 0.d0)                                          &
                        & reduce(i, j) = (-rpar(j))                                           &
                        &**(-dble(i - ncg + npar + nparsw + nparinv - 1 + initpowersw))
            end do
        end do
    elseif (initpowersw .eq. 0 .and. nparsw .ne. 0) then
        do j = 1, kp_ion
            if (rpar(j) .lt. 0.d0) reduce(ncg - npar - nparsw - nparinv + 1, j)        &
                    & = dlog(-rpar(j))
        end do
        do i = ncg - npar - nparsw - nparinv + 2, ncg - npar - nparinv
            !        sostituire direzioni con una parametrizzazione prefissata
            !        a npar parametri
            !        direzioni d labda_{i,j}
            !        introdurre un modulo che definisce
            do j = 1, kp_ion
                if (rpar(j) .lt. 0.d0) reduce(i, j) = (-rpar(j))**                  &
                        & (-dble(i - ncg + npar + nparsw + nparinv - 1))
            end do
        end do
    elseif (nparsw .ne. 0 .and. initpowersw .eq. -1) then
        do i = ncg - nparsw - npar - nparinv + 1, ncg - npar - nparinv
            !        sostituire direzioni con una parametrizzazione prefissata
            !        a npar parametri
            !        direzioni d labda_{i,j}
            !        introdurre un modulo che definisce
            nfun = i + nparsw + npar + nparinv - ncg
            do j = 1, kp_ion
                if (rpar(j) .lt. 0.d0)                                          &
                        & reduce(i, j) = funloc(nfun, -rpar(j), rmax, rmin)
            end do
        end do
        deallocate (rmax, rmin)
    end if
    if (initpowersw .eq. -1 .or. initpower .eq. -1 .or. initpowerinv .eq. -1) &
            &   deallocate (rind, imap)
    return
end

subroutine preprpar(rpar, kp_ion, iond, nion, kiontotj            &
        &, nozeroj_c, nnozeroj_c, nelorbj_c, jbraj, iesfree, indfree             &
        &, jbrajsz, iesinv, indinv, orbcost                                    &
        &, kiontot, nozero_c, nnozero_c, nelorb_c, jbra, iessw, indsw             &
        &, adrlambda, whereiesm, iesm, whereiesup, iesup, &
        &iond_cart, typeorb, nshellj_c, multij_c, ioccj_c, occj_c&
        &, jas_invariant, orbps)
    use Constants, only: ipc, ipj
    use allio, only: symmagp, yes_correct, yes_sparse
    implicit none
    integer kp_ion, i, ii, iii, j, ix, iy, jbraj(*), nozeroj_c(*)    &
            &, nelorbj_c, kiontotj(*), iesfree, indfree, nion, nnozeroj_c            &
            &, nozero_c(*), iessw, indsw, nnozero_c, nelorb_c, jbra(*), kiontot(*)    &
            &, iesinv, indinv, jbrajsz(*), whereiesup(*), iesup, adrlambda(2, *)      &
            &, whereiesm(*), iesm, typeorb(*), multij_c(*)                         &
            &, nshellj_c, ioccj_c(*), occj_c, icek, icek2, jcek
    real*8 iond(nion, nion), rpar(*), jas_invariant(4, *), &
            & iond_cart(3, nion, nion), pip, vec_r(3)
    logical orbcost(*), orbps(*)
    integer k1, k2
    integer*8 ind
    !          indfree=primo indirizzo da riempire-1

    call dscalzero(kp_ion, 0.d0, rpar, 1)

    if (iesfree .ne. 0) then

        do i = 1, nnozeroj_c
            j = abs(jbraj(i))
            if (j .ne. 0) then
                if (yes_sparse) then
                    ix = nozeroj_c(i)
                    iy = nozeroj_c(i + nnozeroj_c)
                    ind = ix + int(ipj*nelorbj_c, 8)*(iy - 1)
                else
                    ind = nozeroj_c(i)
                    iy = (ind - 1)/(ipj*nelorbj_c) + 1
                    ix = ind - (iy - 1)*nelorbj_c*ipj
                end if
                if (ix .ne. iy .and. iy .le. ipj*nelorbj_c .and. .not. orbcost(ix)       &
                        &.and. .not. orbcost(iy)) then
                    adrlambda(1, j + indfree) = kiontotj(ix)
                    adrlambda(2, j + indfree) = kiontotj(iy)

                    rpar(j + indfree) = iond(kiontotj(ix), kiontotj(iy))

                    if (kiontotj(ix) .ne. kiontotj(iy)) then
                        k1 = typeorb(ix)
                        k2 = typeorb(iy)
                        ! In order to compute the invariant 2 in a way that is invariant
                        ! we have to define the direction from p to s, so change sign if k1<k2
                        ! (s->p).
                        if (k1 .gt. k2) then
                            do icek = 1, 3
                                vec_r(icek) = iond_cart(icek, kiontotj(ix), kiontotj(iy))&
                                        & /iond(kiontotj(ix), kiontotj(iy)) ! the unit vecors obviously
                            end do
                        else
                            do icek = 1, 3
                                vec_r(icek) = iond_cart(icek, kiontotj(iy), kiontotj(ix))&
                                        & /iond(kiontotj(iy), kiontotj(ix)) ! the unit vecors obviously
                            end do
                        end if
                        call evaluate_invariant(vec_r, ix, iy, typeorb, &
                                &  jas_invariant(1, j + indfree), orbps(j + indfree))
                    end if

                else
                    if (orbcost(ix) .and. orbcost(iy)) then
                        adrlambda(1, j + indfree) = 0
                        adrlambda(2, j + indfree) = 0
                    elseif (orbcost(ix)) then
                        adrlambda(1, j + indfree) = kiontotj(iy)
                        adrlambda(2, j + indfree) = kiontotj(iy)
                    elseif (orbcost(iy)) then
                        adrlambda(1, j + indfree) = kiontotj(ix)
                        adrlambda(2, j + indfree) = kiontotj(ix)
                    elseif (ix .eq. iy) then
                        adrlambda(1, j + indfree) = kiontotj(ix)
                        adrlambda(2, j + indfree) = kiontotj(ix)
                    end if
                end if
            end if
        end do

    end if

    !############################################################

    if (iesinv .ne. 0) then
        do i = 1, nnozeroj_c
            j = abs(jbrajsz(i))
            if (j .ne. 0) then
                ind = nozeroj_c(i)
                iy = (ind - 1)/nelorbj_c + 1
                ix = ind - (iy - 1)*nelorbj_c
                if (ix .ne. iy .and. iy .le. nelorbj_c .and. .not. orbcost(ix)       &
                        &.and. .not. orbcost(iy)) then
                    rpar(j + indinv) = iond(kiontotj(ix), kiontotj(iy))
                    adrlambda(1, j + indinv) = kiontotj(ix)
                    adrlambda(2, j + indinv) = kiontotj(iy)
                    if (kiontotj(ix) .ne. kiontotj(iy)) then
                        k1 = typeorb(ix)
                        k2 = typeorb(iy)
                        ! In order to compute the invariant 2 in a way that is invariant
                        ! we have to define the direction from p to s, so vec_r changes its sign if k1<k2
                        ! (s->p).
                        if (k1 .gt. k2) then
                            do icek = 1, 3
                                vec_r(icek) = iond_cart(icek, kiontotj(ix), kiontotj(iy))&
                                        & /iond(kiontotj(ix), kiontotj(iy)) ! the unit vecors obviously
                            end do
                        else
                            do icek = 1, 3
                                vec_r(icek) = iond_cart(icek, kiontotj(iy), kiontotj(ix))&
                                        & /iond(kiontotj(iy), kiontotj(ix)) ! the unit vecors obviously
                            end do
                        end if
                        call evaluate_invariant(vec_r, ix, iy, typeorb, &
                                &  jas_invariant(1, j + indinv), orbps(j + indinv))
                    end if
                else
                    if (orbcost(ix) .and. orbcost(iy)) then
                        adrlambda(1, j + indinv) = 0
                        adrlambda(2, j + indinv) = 0
                    elseif (orbcost(ix)) then
                        adrlambda(1, j + indinv) = kiontotj(iy)
                        adrlambda(2, j + indinv) = kiontotj(iy)
                    elseif (orbcost(iy)) then
                        adrlambda(1, j + indinv) = kiontotj(ix)
                        adrlambda(2, j + indinv) = kiontotj(ix)
                    elseif (ix .eq. iy) then
                        adrlambda(1, j + indinv) = kiontotj(ix)
                        adrlambda(2, j + indinv) = kiontotj(ix)
                    end if
                end if
            end if
        end do

    end if

    if (iessw .ne. 0) then
        if (ipc .eq. 1) then
            do i = 1, nnozero_c
                j = abs(jbra(i))
                if (j .ne. 0) then
                    ind = nozero_c(i)
                    iy = (ind - 1)/nelorb_c + 1
                    ix = ind - (iy - 1)*nelorb_c
                    if (ix .ne. iy .and. iy .le. nelorb_c) then
                        if (kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) rpar(j + indsw) = -iond(kiontot(ix), kiontot(iy))
                        adrlambda(1, j + indsw) = kiontot(ix)
                        adrlambda(2, j + indsw) = kiontot(iy)
                    else
                        adrlambda(1, j + indsw) = kiontot(ix)
                        adrlambda(2, j + indsw) = kiontot(ix)
                    end if
                end if
            end do
        elseif (symmagp .and. yes_correct) then
            do i = 1, nnozero_c
                j = abs(jbra(i))
                if (j .ne. 0) then
                    ind = nozero_c(i)
                    iy = (ind - 1)/nelorb_c + 1
                    ix = ind - (iy - 1)*nelorb_c
                    if (ix .ne. iy .and. iy .le. nelorb_c) then
                        if (kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0)&
                           & rpar(4*j - 3 + indsw:4*j + indsw) = -iond(kiontot(ix), kiontot(iy))
                        adrlambda(1, 4*j - 3 + indsw:4*j + indsw) = kiontot(ix)
                        adrlambda(2, 4*j - 3 + indsw:4*j + indsw) = kiontot(iy)
                    else
                        adrlambda(1, 4*j - 3 + indsw:4*j + indsw) = kiontot(ix)
                        adrlambda(2, 4*j - 3 + indsw:4*j + indsw) = kiontot(ix)
                    end if
                end if
            end do
        else
            do i = 1, nnozero_c
                j = abs(jbra(i))
                if (j .ne. 0) then
                    ind = nozero_c(i)
                    iy = (ind - 1)/nelorb_c + 1
                    ix = ind - (iy - 1)*nelorb_c
                    if (ix .ne. iy .and. iy .le. nelorb_c) then
                        if (kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) rpar(2*j - 1 + indsw) = -iond(kiontot(ix), kiontot(iy))
                        if (kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) rpar(2*j + indsw) = -iond(kiontot(ix), kiontot(iy))
                        adrlambda(1, 2*j - 1 + indsw) = kiontot(ix)
                        adrlambda(1, 2*j + indsw) = kiontot(ix)
                        adrlambda(2, 2*j - 1 + indsw) = kiontot(iy)
                        adrlambda(2, 2*j + indsw) = kiontot(iy)
                    else
                        adrlambda(1, 2*j - 1 + indsw) = kiontot(ix)
                        adrlambda(1, 2*j + indsw) = kiontot(ix)
                        adrlambda(2, 2*j - 1 + indsw) = kiontot(ix)
                        adrlambda(2, 2*j + indsw) = kiontot(ix)
                    end if
                end if
            end do
        end if
    end if
    if (iesup .ne. 0) then
        if (ipc .eq. 2) then
            do i = 1, iesup/2
                adrlambda(1, 2*i - 1 + indsw + iessw) = whereiesup(i)
                adrlambda(1, 2*i + indsw + iessw) = whereiesup(i)
                adrlambda(2, 2*i - 1 + indsw + iessw) = whereiesup(i)
                adrlambda(2, 2*i + indsw + iessw) = whereiesup(i)
            end do
        else
            do i = 1, iesup
                adrlambda(1, i + indsw + iessw) = whereiesup(i)
                adrlambda(2, i + indsw + iessw) = whereiesup(i)
            end do
        end if
    end if
    if (iesm .ne. 0) then
        do i = 1, iesm
            adrlambda(1, i + iesinv) = whereiesm(i)
            adrlambda(2, i + iesinv) = whereiesm(i)
        end do
    end if
    return
end

function funloc(n, r, rmax, rmin)
    implicit none
    integer n
    real*8 r, rmax(*), rmin(*), funloc
    if (r .lt. rmax(n) .and. r .ge. rmin(n)) then
        funloc = 1.d0
    else
        funloc = 0.d0
    end if
    return
end
subroutine preprminmax(kp, npar, rpar, rmin, rmax, rind, imap          &
        &, iopt, endinv)
    implicit none
    integer kp, npar, ndist, ngroup, nrest, ind, i, imap(*)                  &
            &, iopt, endinv
    real*8 rpar(*), rmin(*), rmax(*), rind(*), eps, dist, rmint, rmaxt, delta
    !     compute how many independent distances are present
    eps = 1d-6
    ndist = 0
    if (iopt .eq. 0) then
        do i = 1, kp
            if (rpar(i) .lt. 0) then
                ndist = ndist + 1
                rind(ndist) = -rpar(i)
            end if
        end do
    elseif (iopt .eq. 1) then
        do i = endinv + 1, kp
            if (rpar(i) .gt. 0) then
                ndist = ndist + 1
                rind(ndist) = rpar(i)
            end if
        end do
    elseif (iopt .eq. 2) then
        do i = 1, endinv
            if (rpar(i) .gt. 0) then
                ndist = ndist + 1
                rind(ndist) = rpar(i)
            end if
        end do
    end if

    !     dividing the full interval in npar subintervals

    call dsortx(rind, 1, ndist, imap)
    rmint = rind(1) - eps
    rmaxt = rind(ndist) + eps
    delta = (rmaxt - rmint)/npar

    do i = 1, npar
        rmin(i) = rmint + (i - 1)*delta
        rmax(i) = rmint + i*delta
    end do
    return
end
function findindex(iar, ibr, nmax, allfit)
    implicit none
    integer ia, ib, iar, ibr, nmax, findindex
    logical allfit
    if (allfit) then
        findindex = nmax*(ibr - 1) + iar
    else
        if (iar .gt. ibr) then
            ia = iar
            ib = ibr
        else
            ia = ibr
            ib = iar
        end if
        findindex = (ia - ib + 1) + nmax*(ib - 1) - ((ib - 2)*(ib - 1))/2
    end if
    return
end

subroutine attach_phase2det(minus, detmat_c)
    use allio, only: rion, kiontot, symmagp, yes_hermite, sjbradet&
            &, jbradet, nnozero_c, nozero_c, nelorb_c, nelorb_at, opposite_phase, rank, ndiff, pfaffup
    use cell, only: cellscale, phase2pi, phase2pi_down, car2cry, CartesianToCrystal
    use constants, only: ipc, deps, ipf
    implicit none

    !nelorb_cm is half of nelorb_c-ndiff if ipf.eq.2
    integer npip(3), npipn(3), ix, iy, ind, i, k, nelorb_cm
    complex*16 exp_phase
    real*8 detmat_c(ipc*nelorb_c, nelorb_c)
    real*8 riondiff(3), cost, cost_old, rdiff(3), phase2pi1(3), phase2pi2(3)
    logical minus, dolong, dodolong
    integer, external :: makenpip

    !       minus=.true.  From effective parametrization to real determinant
    !       minus=.false. From real to effective parametrization
    !       NB In te case symmagp=.false. also the hermitian conjugate is an
    !       independent transformation.
    if (ipf .eq. 2) nelorb_cm = (nelorb_c)/2
    if (ipf .eq. 1) nelorb_cm = nelorb_c

    !CCC questo è uguale se ipf=2, no?
    if (symmagp .and. opposite_phase .and. yes_hermite) then
        dolong = .true.
    else
        dolong = .false.
    end if

    do i = 1, nnozero_c
        ind = nozero_c(i)
        iy = (ind - 1)/nelorb_c + 1
        ix = ind - (iy - 1)*nelorb_c
        !checking if the first phase is given by an up electron
        if (ipf .eq. 1 .or. ix .le. nelorb_cm) then
            phase2pi1(:) = phase2pi(:)
            if (ipf .eq. 1 .or. iy .gt. nelorb_cm) then
                phase2pi2(:) = phase2pi_down(:)
                dodolong = .true.
            else
                phase2pi2(:) = phase2pi(:)
                dodolong = .false.
            end if
        else
            phase2pi1(:) = phase2pi_down(:)
            if (iy .gt. nelorb_cm) then
                phase2pi2(:) = phase2pi_down(:)
                dodolong = .false.
            else
                phase2pi2(:) = phase2pi(:)
                dodolong = .true.
            end if
        end if

        if (.not. dolong) dodolong = .false.

        if (kiontot(ix) .ne. 0 .and. kiontot(iy) .ne. 0) then
            if (ix .ne. iy) then
                riondiff(:) = rion(:, kiontot(ix)) - rion(:, kiontot(iy))
                !   L_mu of the notes (parbcs.tex) corresponds to -npip(:)*cellscale(:)
                call CartesianToCrystal(riondiff, 1)
                do k = 1, 3
                    npip(k) = makenpip(riondiff(k), cellscale(k), deps)
                end do
                if (dodolong) then
                    rdiff(:) = 2*riondiff(:)
                    call makeimagep(rdiff, npipn, cellscale, deps)
                end if
                if (kiontot(ix) .ne. kiontot(iy) .and. abs(sum(rdiff(:)**2)) .lt. deps**2 .and. dodolong) then
                    if (ipc .eq. 2) then
                        if (minus) then
                            cost = sum(npipn(:)/4.d0*(phase2pi1(:) - phase2pi2(:)))
                        else
                            cost = -sum(npipn(:)/4.d0*(phase2pi1(:) - phase2pi2(:)))
                        end if
                    else
                        cost = 0.d0 ! do nothing
                    end if
                else
                    !        npip * cellscale == - L_nu in the notes parbcs.tex because determined
                    !        by      R(ix)-cellscale(:)*npip is inside the reference < cellscale/2
                    if (minus) then
                        cost = sum(npip(:)/2.d0*(phase2pi1(:) - phase2pi2(:)))
                    else
                        cost = -sum(npip(:)/2.d0*(phase2pi1(:) - phase2pi2(:)))
                    end if
                end if
                exp_phase = dcmplx(dcos(cost), dsin(cost))
            else ! ix=iy
                exp_phase = (1.d0, 0.d0)
            end if ! ix=iy

            if (minus) then
                cost = sum((rion(:, kiontot(iy)) + rion(:, kiontot(ix)))*(phase2pi1(:) + phase2pi2(:))/cellscale(1:3))/2.d0
            else
                cost = -sum((rion(:, kiontot(iy)) + rion(:, kiontot(ix)))*(phase2pi1(:) + phase2pi2(:))/cellscale(1:3))/2.d0
            end if
            exp_phase = exp_phase*dcmplx(dcos(cost), dsin(cost))

            if (ipc .eq. 2) then
                call zscal(1, exp_phase, detmat_c(2*ix - 1, iy), 1)
            else
                call dscal(1, dreal(exp_phase), detmat_c(ix, iy), 1)
            end if

            if (ipf .eq. 2 .and. ix .ne. iy) then
                if (minus) then
                    if (ix .le. nelorb_cm) then
                        if (iy .le. nelorb_cm) then
                            if (ipc .eq. 2) then
                                if (symmagp .and. yes_hermite) then
                                    detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                                    if (.not. pfaffup) then
                                        detmat_c(2*(ix + nelorb_cm) - 1, iy + nelorb_cm) = detmat_c(2*ix - 1, iy)
                                        detmat_c(2*(ix + nelorb_cm), iy + nelorb_cm) = -detmat_c(2*ix, iy)
                                        detmat_c(2*(iy + nelorb_cm) - 1, ix + nelorb_cm) = -detmat_c(2*ix - 1, iy)
                                        detmat_c(2*(iy + nelorb_cm), ix + nelorb_cm) = detmat_c(2*ix, iy)
                                    end if
                                else if (symmagp) then
                                    detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                                    if (.not. pfaffup) then
                                        detmat_c(2*(ix + nelorb_cm) - 1, iy + nelorb_cm) = detmat_c(2*ix - 1, iy)
                                        detmat_c(2*(ix + nelorb_cm), iy + nelorb_cm) = detmat_c(2*ix, iy)
                                        detmat_c(2*(iy + nelorb_cm) - 1, ix + nelorb_cm) = -detmat_c(2*ix - 1, iy)
                                        detmat_c(2*(iy + nelorb_cm), ix + nelorb_cm) = -detmat_c(2*ix, iy)
                                    end if
                                else
                                    detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                                end if
                            else
                                if (symmagp) then
                                    detmat_c(iy, ix) = -detmat_c(ix, iy)
                                    if (.not. pfaffup) then
                                        detmat_c(ix + nelorb_cm, iy + nelorb_cm) = detmat_c(ix, iy)
                                        detmat_c(iy + nelorb_cm, ix + nelorb_cm) = -detmat_c(ix, iy)
                                    end if
                                else
                                    detmat_c(iy, ix) = -detmat_c(ix, iy)
                                end if
                            end if
                        else
                            if (ipc .eq. 2) then
                                if (symmagp .and. yes_hermite) then
                                    detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                                    detmat_c(2*(iy - nelorb_cm) - 1, ix + nelorb_cm) = detmat_c(2*ix - 1, iy)
                                    detmat_c(2*(iy - nelorb_cm), ix + nelorb_cm) = -detmat_c(2*ix, iy)
                                    detmat_c(2*(ix + nelorb_cm) - 1, iy - nelorb_cm) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*(ix + nelorb_cm), iy - nelorb_cm) = detmat_c(2*ix, iy)
                                else if (symmagp) then
                                    detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                                    detmat_c(2*(iy - nelorb_cm) - 1, ix + nelorb_cm) = detmat_c(2*ix - 1, iy)
                                    detmat_c(2*(iy - nelorb_cm), ix + nelorb_cm) = detmat_c(2*ix, iy)
                                    detmat_c(2*(ix + nelorb_cm) - 1, iy - nelorb_cm) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*(ix + nelorb_cm), iy - nelorb_cm) = -detmat_c(2*ix, iy)
                                else
                                    detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                                    detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                                end if
                            else
                                if (symmagp) then
                                    detmat_c(iy, ix) = -detmat_c(ix, iy)
                                    detmat_c(iy - nelorb_cm, ix + nelorb_cm) = detmat_c(ix, iy)
                                    detmat_c(ix + nelorb_cm, iy - nelorb_cm) = -detmat_c(ix, iy)
                                else
                                    detmat_c(iy, ix) = -detmat_c(ix, iy)
                                end if
                            end if
                        end if
                    else
                        if (ipc .eq. 2) then
                            detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                            detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                        else
                            detmat_c(iy, ix) = -detmat_c(ix, iy)
                        end if
                    end if
                else
                    !Inizia il minus dello pfaff
                    if (symmagp) then
                        if (opposite_phase) exp_phase = dconjg(exp_phase)

                        !If symm agp the first here we have only up-up and
                        !up-down
                        if (iy .gt. nelorb_cm) then
                            if (ipc .eq. 2) then
                                call zscal(1, exp_phase, detmat_c(2*(iy - nelorb_cm) - 1, ix + nelorb_cm), 1)
                                detmat_c(2*(ix + nelorb_cm) - 1, iy - nelorb_cm) =&
                                   & -detmat_c(2*(iy - nelorb_cm) - 1, ix + nelorb_cm)
                                detmat_c(2*(ix + nelorb_cm), iy - nelorb_cm) = -detmat_c(2*(iy - nelorb_cm), ix + nelorb_cm)
                            else
                                call dscal(1, dreal(exp_phase), detmat_c(iy - nelorb_cm, ix + nelorb_cm), 1)
                                detmat_c(ix + nelorb_cm, iy - nelorb_cm) = -detmat_c(iy - nelorb_cm, ix + nelorb_cm)
                            end if
                        else
                            if (.not. pfaffup) then
                                if (ipc .eq. 2) then
                                    call zscal(1, exp_phase, detmat_c(2*(iy + nelorb_cm) - 1, ix + nelorb_cm), 1)
                                    detmat_c(2*(ix + nelorb_cm) - 1, iy + nelorb_cm) =&
                                       & -detmat_c(2*(iy + nelorb_cm) - 1, ix + nelorb_cm)
                                    detmat_c(2*(ix + nelorb_cm), iy + nelorb_cm) = -detmat_c(2*(iy + nelorb_cm), ix + nelorb_cm)
                                else
                                    call dscal(1, dreal(exp_phase), detmat_c(iy + nelorb_cm, ix + nelorb_cm), 1)
                                    detmat_c(ix + nelorb_cm, iy + nelorb_cm) = -detmat_c(iy + nelorb_cm, ix + nelorb_cm)
                                end if
                            end if
                        end if
                    end if
                    if (ipc .eq. 1) then
                        detmat_c(iy, ix) = -detmat_c(ix, iy)
                    else
                        detmat_c(2*iy - 1, ix) = -detmat_c(2*ix - 1, iy)
                        detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                    end if
                end if
                !Finisce il minus dello pfaff
            elseif (symmagp .and. ix .ne. iy) then
                !                 From effective to real
                if (minus) then

                    if (ipc .eq. 2) then
                        if (yes_hermite) then
                            detmat_c(2*iy - 1, ix) = detmat_c(2*ix - 1, iy)
                            detmat_c(2*iy, ix) = -detmat_c(2*ix, iy)
                        else
                            detmat_c(2*iy - 1, ix) = detmat_c(2*ix - 1, iy)
                            detmat_c(2*iy, ix) = detmat_c(2*ix, iy)
                        end if
                    else
                        detmat_c(iy, ix) = detmat_c(ix, iy)
                    end if
                else
                    ! From real to effective: NB we cannot assume that detmat_c is hermitian or
                    ! symmetric as it is applied to derivative matrices after AAD (computeb_global)
                    if (opposite_phase) exp_phase = dconjg(exp_phase)
                    if (ipc .eq. 2) then
                        call zscal(1, exp_phase, detmat_c(2*iy - 1, ix), 1)
                    else
                        call dscal(1, dreal(exp_phase), detmat_c(iy, ix), 1)
                    end if
                end if
            end if
        end if
    end do
    return
end subroutine attach_phase2det

function makenpip_fake(rdiff, cellscale, deps)
    implicit none
    integer kk, makenpip_fake, npip
    real*8 cost, deps, rdiff, cellscale
    !   L/2 and -L/2 are not equivalent in this routine.
    cost = rdiff/cellscale
    npip = nint(2*cost)
    if ((npip/2)*2 .ne. npip .and. abs(cost - npip/2.d0) .lt. deps) then
        if (npip .eq. -1) then
            makenpip_fake = 0
        elseif (npip .gt. 0) then
            makenpip_fake = (npip - 1)/2
        else
            makenpip_fake = (npip + 1)/2
        end if
    else
        makenpip_fake = nint(cost)
    end if
    return
end function makenpip_fake

function makenpip(rdiff, cellscale, eps)
    implicit none
    integer kk, makenpip
    real*8 cut, cost, eps, rdiff, cellscale
    !       rdiff -cellscale* makenpip is inside the reference -L/2<= rdiff < L/2.
    !       NB L/2 is mapped to -L/2 in this unique definition.
    cut = 0.5d0 - eps
    cost = rdiff/cellscale
    !  if cost = n+0.5-small  makenpip--> n , cost  -npip = 0.5-small -->
    !   makenpip--mekenpip+1. cost=n+0.5+small, makenpip -> n+1, cost-makenpip = -0.5-small
    makenpip = nint(cost)
    if (cost - makenpip .gt. cut) makenpip = makenpip + 1
    return
end function makenpip

subroutine makeimage(rdiff, cellscale, eps)
    !   L/2 and -L/2 are equivalent in this routine.
    implicit none
    integer kk, npip
    real*8 cut, cost, eps, rdiff(3), cellscale(3)
    cut = 0.5d0 - eps
    do kk = 1, 3
        cost = rdiff(kk)/cellscale(kk)
        npip = nint(cost)
        if (cost - npip .gt. cut) npip = npip + 1
        rdiff(kk) = rdiff(kk) - cellscale(kk)*npip
    end do
    return
end

subroutine makeimagep(rdiff, npip, cellscale, eps)
    !   L/2 and -L/2 are equivalent in this routine. Namely L/2 -> -L/2
    implicit none
    integer kk, npip(3)
    real*8 cut, cost, eps, rdiff(3), cellscale(3)
    cut = 0.5d0 - eps
    do kk = 1, 3
        cost = rdiff(kk)/cellscale(kk)
        npip(kk) = nint(cost)
        if (cost - npip(kk) .gt. cut) npip(kk) = npip(kk) + 1
        rdiff(kk) = rdiff(kk) - cellscale(kk)*npip(kk)
        ! BUG 1/8/2016  with the line below we do not put any phase to the
        ! long distance bonds (see parbcs.tex).
        !     if(abs(rdiff(kk)+cellscale(kk)/2).gt.eps) rdiff(kk)=0.d0
    end do
    return
end

subroutine makeimage_fake(rdiff, cellscale, deps)
    implicit none
    !   L/2 and -L/2 are not equivalent in this routine.
    integer kk, npip, m
    real*8 cost, rdiff(3), cellscale(3), deps
    do kk = 1, 3
        cost = rdiff(kk)/cellscale(kk)
        npip = nint(2*cost)
        if ((npip/2)*2 .ne. npip .and. abs(cost - npip/2.d0) .lt. deps) then ! Boarder cases
            if (npip .eq. -1) then
                m = 0
            elseif (npip .gt. 0) then
                m = (npip - 1)/2
            else
                m = (npip + 1)/2
            end if
        else
            m = nint(cost)
        end if
        rdiff(kk) = rdiff(kk) - cellscale(kk)*m
    end do
    return
end
