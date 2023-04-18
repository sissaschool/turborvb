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

subroutine uptabpip(indt, nel, nelup, jel                           &
        &, rkel, rkeln, tabpip, tcost, iesdr, vj, winvsj, winvj, nelorbj            &
        &, psip, psip_store, winvjbar, winvjbarsz, winvjbarn, winvjbarszn, rion, nion         &
        &, costz, costz3, iessz, LBox, n_body_on                                &
      &, jastrowall_ee, jasnew_ei, rmu, niesd, indtm, nelorbjh, psiln, no_dgemv)

    use Cell, only: map, cellscale, metric, car2cry
    use Constants, only: ipj
    use allio, only: pointvj, yes_ontarget, norm_metric
    implicit none

    integer nel, ii, jel, indt, i, j, k, iesd, iesdr, nelup, ispin, ispinj   &
     &, nelorbj5, nelorbj, indtold, indtnew, indtoldsz, indtnewsz, indt5, nion&
            &, neldo, indjel, iesdr2iesd, kk, niesd, indtm(nel), indtp, nelorbjh, n_body_on&
            &, indtmn, indt4n
    real*8 tabpip(nel, indt + 4), cost, vj(max(1, niesd)), tcost(indt)           &
            &, zj, xj, yj, zjo, yjo, xjo, xjn, yjn, zjn, xi, yi, zi, n1, n2, n3, n4, jastrow_ee &
            &, jastrow1, jastrow2, grad1(3), grad2(3), grad12, grad22, rc(3), rcn(3)   &
            &, winvsj(max(nelorbj, 1), 0:indt + 4), winvj(max(nelorbj, 1), 0:indt + 4, nel)&
            &, psip(2*indt + 9, nel), psip_store(indt + 5, nel), costj, winvjbar(max(ipj*nelorbjh, 1), nel)&
            &, winvjbarn(max(ipj*nelorbjh, 1)), rion(3, nion), hess(3)&
            &, costz(nion), costz3(nion), costz0, jastrow_ei, winvjbarsz(*), hess2(3)&
            &, winvjbarszn(*), signszup, signszdo, costjnew, r0, psiln, rt(3), csum, costtab(4)&
            &, csum1, csum2, csum3, csum4
    real*8 rkel(3, nel, 0:indt), rkeln(3, 0:indt), LBox
    real*8 jastrowall_ee(nel, nel, 0:indt + 4), jasnew_ei(nion)
    real*8 rmu(3, 0:indt, nion)
    logical iessz, iesspin, no_dgemv
    integer lpsip, indt_all

!   if(yes_ontarget) then

    lpsip = 2*indt + 9

    iesd = iesdr2iesd(iesdr)

    if (iesdr .eq. -7 .or. iesd .eq. 2 .or. iesd .lt. 0) then
        iesspin = .true.
    else
        iesspin = .false.
    end if

    if (iessz) then
        indjel = (jel - 1)*nelorbjh + 1
    else
        indjel = 1
    end if

    !         winvjbar is here before  the updating
    !         winvbarjn  is the new raw
    indtmn = indtm(jel)
    indt4n = indtmn + 4

    if (nelorbj .ne. 0) then

        do k = 1, 4
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
            do i = 1, nelorbj
                winvsj(i, indtmn + k) = winvsj(i, indt + k)
            end do
        end do

        neldo = nel - nelup
        nelorbj5 = nelorbj*(indt + 5)

        !         new terms
        !          call dgemm('T','N',nel,indt4n,nelorbjh,1.d0,winvjbar      &
        !     &,nelorbjh,winvsj(1,1),nelorbj,0.d0,psip,nel)

        if (yes_ontarget) then
            if ((jel <= nelup) .or. (ipj == 1)) then
                call dgemm_tn(indt4n, nel, nelorbjh, 1.d0, winvsj(1, 1), nelorbj, winvjbar &
                          &, ipj*nelorbjh, 0.d0, psip, lpsip)
            else
                call dgemm_tn(indt4n, nel, nelorbjh, 1.d0, winvsj(1, 1), nelorbj&
               &, winvjbar(1 + nelorbjh, 1), 2*nelorbjh, 0.d0, psip, lpsip)
            end if
        else
            if ((jel <= nelup) .or. (ipj == 1)) then
                call dgemm('T', 'N', indt4n, nel, nelorbjh, 1.d0, winvsj(1, 1), nelorbj, winvjbar &
                          &, ipj*nelorbjh, 0.d0, psip, lpsip)
            else
                call dgemm('T', 'N', indt4n, nel, nelorbjh, 1.d0, winvsj(1, 1), nelorbj&
               &, winvjbar(1 + nelorbjh, 1), 2*nelorbjh, 0.d0, psip, lpsip)
            end if
        end if

        indt5 = indt + 5
        indtnew = indt5

        if (no_dgemv) then
            if (ipj == 1) then
#ifdef _OFFLOAD
                if (yes_ontarget) then
!$omp target teams distribute collapse(2) private(csum)
                    do i = 1, nel
                        do k = 0, indt + 4
                            if (k .le. indtm(i) .or. k .gt. indt) then
                                csum = 0.d0
!$omp parallel do reduction(+:csum)
                                do j = 1, nelorbjh
                                    csum = csum + winvj(j, k, i)*winvjbarn(j)
                                end do
                                psip(indtnew + k, i) = csum
                            end if
                        end do
                    end do
!$omp end target teams distribute
                else
!$omp  parallel do default(shared) private(i,j,k,csum)
                    do i = 1, nel
                        do k = 0, indt + 4
                            if (k .le. indtm(i) .or. k .gt. indt) then
                                csum = 0.d0
                                do j = 1, nelorbjh
                                    csum = csum + winvj(j, k, i)*winvjbarn(j)
                                end do
                                psip(indtnew + k, i) = csum
                            end if
                        end do
                    end do
!$omp end parallel do
                end if
#else
!$omp  parallel do default(shared) private(i,j,k,csum)
                do i = 1, nel
                    do k = 0, indt + 4
                        if (k .le. indtm(i) .or. k .gt. indt) then
                            csum = 0.d0
                            do j = 1, nelorbjh
                                csum = csum + winvj(j, k, i)*winvjbarn(j)
                            end do
                            psip(indtnew + k, i) = csum
                        end if
                    end do
                end do
!$omp end parallel do
#endif

                !        do i=indt+1,indt+4
                !       call dgemv('T',nelorbjh,nel,1.d0,winvj(1,i,1),nelorbj5           &
                !     &,winvjbarn,1,0.d0,psip(1,indtnew+i),1)
                !        enddo
            else
#ifdef _OFFLOAD
                if (yes_ontarget) then
!$omp target teams distribute collapse(2) private(csum)
                    do i = 1, nelup
                        do k = 0, indt + 4
                            if (k .le. indtm(i) .or. k .gt. indt) then
                                csum = 0.d0
!$omp parallel do reduction(+:csum)
                                do j = 1, nelorbjh
                                    csum = csum + winvj(j, k, i)*winvjbarn(j)
                                end do
                                psip(indtnew + k, i) = csum
                            end if
                        end do
                    end do
!$omp end target teams distribute
                else
!$omp  parallel do default(shared) private(i,j,k,csum)
                    do i = 1, nelup
                        do k = 0, indt + 4
                            if (k .le. indtm(i) .or. k .gt. indt) then
                                csum = 0.d0
                                do j = 1, nelorbjh
                                    csum = csum + winvj(j, k, i)*winvjbarn(j)
                                end do
                                psip(indtnew + k, i) = csum
                            end if
                        end do
                    end do
!$omp end parallel do
                end if
#else
!$omp  parallel do default(shared) private(i,j,k,csum)
                do i = 1, nelup
                    do k = 0, indt + 4
                        if (k .le. indtm(i) .or. k .gt. indt) then
                            csum = 0.d0
                            do j = 1, nelorbjh
                                csum = csum + winvj(j, k, i)*winvjbarn(j)
                            end do
                            psip(indtnew + k, i) = csum
                        end if
                    end do
                end do
!$omp end parallel do
#endif
#ifdef _OFFLOAD
                if (yes_ontarget) then
!$omp target teams distribute collapse(2) private(csum)
                    do i = nelup + 1, nel
                        do k = 0, indt + 4
                            if (k .le. indtm(i) .or. k .gt. indt) then
                                csum = 0.d0
!$omp parallel do reduction(+:csum)
                                do j = 1, nelorbjh
                                    csum = csum + winvj(j, k, i)*winvjbarn(nelorbjh + j)
                                end do
                                psip(indtnew + k, i) = csum
                            end if
                        end do
                    end do
!$omp end target teams distribute
                else
!$omp  parallel do default(shared) private(i,j,k,csum)
                    do i = nelup + 1, nel
                        do k = 0, indt + 4
                            if (k .le. indtm(i) .or. k .gt. indt) then
                                csum = 0.d0
                                do j = 1, nelorbjh
                                    csum = csum + winvj(j, k, i)*winvjbarn(nelorbjh + j)
                                end do
                                psip(indtnew + k, i) = csum
                            end if
                        end do
                    end do
!$omp end parallel do
                end if
#else
!$omp  parallel do default(shared) private(i,j,k,csum)
                do i = nelup + 1, nel
                    do k = 0, indt + 4
                        if (k .le. indtm(i) .or. k .gt. indt) then
                            csum = 0.d0
                            do j = 1, nelorbjh
                                csum = csum + winvj(j, k, i)*winvjbarn(nelorbjh + j)
                            end do
                            psip(indtnew + k, i) = csum
                        end if
                    end do
                end do
!$omp end parallel do
#endif
                !       do i=indt+1,indt+4
                !      call dgemv('T',nelorbjh,nelup,1.d0,winvj(1,i,1),nelorbj5           &
                !    &,winvjbarn,1,0.d0,psip(1,indtnew+i),1)
                !      call dgemv('T',nelorbjh,neldo,1.d0,winvj(1,i,nelup+1),nelorbj5     &
                !    &,winvjbarn(1+nelorbjh),1,0.d0,psip(nelup+1,indtnew+i),1)
                !       enddo
            end if ! endif ipj
        else
            if (ipj .eq. 1) then
                call dgemv_('T', nelorbjh, nel*indt5, 1.d0, winvj(1, 0, 1), nelorbjh&
                      &, winvjbarn, 1, 0.d0, psip_store, 1, yes_ontarget)
                call copy_simple(indt5, lpsip, nel, psip_store, psip(indtnew, 1), yes_ontarget)
            else
                call dgemv_('T', nelorbjh, nelup*indt5, 1.d0, winvj(1, 0, 1), nelorbjh&
                      &, winvjbarn, 1, 0.d0, psip_store, 1, yes_ontarget)
                call dgemv_('T', nelorbjh, neldo*indt5, 1.d0, winvj(1, 0, nelup + 1)&
            &, nelorbjh, winvjbarn(nelorbjh + 1), 1, 0.d0, psip_store(1, nelup + 1), 1, yes_ontarget)
                call copy_simple(indt5, lpsip, nel, psip_store, psip(indtnew, 1), yes_ontarget)
            end if
        end if

        if (iessz) then
#ifdef _OFFLOAD
!$omp target update from (psip) if(yes_ontarget)
#endif
            if (jel .le. nelup) then
                signszup = 1.d0
                signszdo = -1.d0
            else
                signszup = -1.d0
                signszdo = 1.d0
            end if

            !      new terms
            call dgemm('T', 'N', indt4n, nelup, nelorbjh, signszup, winvsj(1, 1), nelorbj&
           &, winvjbarsz, nelorbjh, 1.d0, psip, lpsip)
            call dgemm('T', 'N', indt4n, neldo, nelorbjh, signszdo, winvsj(1, 1), nelorbj&
                    &, winvjbarsz(nelup*nelorbjh + 1), nelorbjh        &
                    &, 1.d0, psip(1, nelup + 1), lpsip)
            !      end new terms

            if (indt .ne. 0) then

                do i = 1, nelup
                    call dgemv('T', nelorbjh, indtm(i) + 1, signszup, winvj(1, 0, i), nelorbj &
                            &, winvjbarszn, 1, 1.d0, psip(indtnew, i), 1)
                end do
                do i = nelup + 1, nel
                    call dgemv('T', nelorbjh, indtm(i) + 1, signszdo, winvj(1, 0, i), nelorbj  &
                            &, winvjbarszn, 1, 1.d0, psip(indtnew, i), 1)
                end do
                do i = indt + 1, indt + 4
                    call dgemv('T', nelorbjh, nelup, signszup, winvj(1, i, 1), nelorbj5      &
                            &, winvjbarszn, 1, 1.d0, psip(indtnew + i, 1), lpsip)
                    call dgemv('T', nelorbjh, neldo, signszdo, winvj(1, i, nelup + 1) &
                 &, nelorbj5, winvjbarszn, 1, 1.d0, psip(indtnew + i, nelup + 1), lpsip)
                end do
            else
                do i = 0, indt + 4
                    call dgemv('T', nelorbjh, nelup, signszup, winvj(1, i, 1), nelorbj5      &
                            &, winvjbarszn, 1, 1.d0, psip(indtnew + i, 1), lpsip)
                    call dgemv('T', nelorbjh, neldo, signszdo, winvj(1, i, nelup + 1)&
               &, nelorbj5, winvjbarszn, 1, 1.d0, psip(indtnew + i, nelup + 1), lpsip)
                end do
            end if

        else

#ifdef _OFFLOAD
!$omp target update from (psip) if(yes_ontarget)
#endif
            ! endif iessz
        end if
        !    Restoring original winvsj
        if (indt .ne. indtmn) then
            do k = 4, 1, -1
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelorbj
                    winvsj(i, indt + k) = winvsj(i, indtmn + k)
                end do
            end do
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
            do k = indtmn + 1, indt
                do i = 1, nelorbj
                    winvsj(i, k) = 0.d0
                end do
            end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
        end if

        ! endif nelorbj
    end if

    zj = rkeln(3, 0)
    yj = rkeln(2, 0)
    xj = rkeln(1, 0)

    zjo = rkel(3, jel, 0)
    yjo = rkel(2, jel, 0)
    xjo = rkel(1, jel, 0)

    if (jel .le. nelup) then
        ispinj = 1
    else
        ispinj = -1
    end if
    do j = 1, nel
        do k = 1, indt
            if (k .gt. indtm(j)) tabpip(j, k) = 0.d0
        end do
    end do

!$omp parallel do default(shared) reduction(+:psiln) &
!$omp private(i,k,ispin,rc,grad1,xi,yi,zi,jastrow1   &
!$omp ,jastrow2,grad12,grad2,grad22,cost,costj,costjnew)

    do i = 1, nel

        if (i .ne. jel) then

            if (iesspin) then
                if (i .le. nelup) then
                    ispin = ispinj
                else
                    ispin = -ispinj
                end if
            else
                ispin = -1
            end if

            zi = rkel(3, i, 0)
            yi = rkel(2, i, 0)
            xi = rkel(1, i, 0)

            rc(1) = xi - xj
            rc(2) = yi - yj
            rc(3) = zi - zj

            if (LBox .gt. 0.d0) then
                call jastrowgrad_pbc(rc, vj, iesd, jastrow1, grad1, grad12, ispin, 1.d0)
            else
                call jastrowgrad(rc, vj, iesd, jastrow1, grad1, grad12, ispin)
            end if

            !     old 3-body grad-lap
            if (nelorbj .ne. 0) then
                !     new 3-body grad-lap
                jastrow1 = jastrow1 + psip(indt5, i)
                !       Update jastrowall the symmetric
                jastrowall_ee(i, jel, 0) = jastrow1
                do k = 1, 3
                    jastrowall_ee(i, jel, indt + k) = -grad1(k) + psip(indtmn + k, i)
                end do
                jastrowall_ee(i, jel, indt + 4) = grad12 + psip(indtmn + 4, i)

                do k = 1, 3
                    grad1(k) = grad1(k) + psip(indtnew + indt + k, i)
                end do
                grad12 = grad12 + psip(indtnew + indt + 4, i)
            else
                !       Update jastrowall the symmetric
                jastrowall_ee(i, jel, 0) = jastrow1
                jastrowall_ee(i, jel, indt + 1:indt + 3) = -grad1(1:3)
                jastrowall_ee(i, jel, indt + 4) = grad12
            end if

            jastrow2 = jastrowall_ee(jel, i, 0)
            grad2(:) = jastrowall_ee(jel, i, indt + 1:indt + 3)
            grad22 = jastrowall_ee(jel, i, indt + 4)
            !       Update jastrowall
            jastrowall_ee(jel, i, 0) = jastrow1
            jastrowall_ee(jel, i, indt + 1:indt + 3) = grad1(1:3)
            jastrowall_ee(jel, i, indt + 4) = grad12

            cost = dexp(jastrow2 - jastrow1)

            psiln = psiln + jastrow1 - jastrow2

            tabpip(i, 1 + indt) = tabpip(i, 1 + indt) + grad1(1)                     &
                    & - grad2(1)
            tabpip(i, 2 + indt) = tabpip(i, 2 + indt) + grad1(2)                     &
                    & - grad2(2)
            tabpip(i, 3 + indt) = tabpip(i, 3 + indt) + grad1(3)                     &
                    & - grad2(3)

            tabpip(i, 4 + indt) = tabpip(i, 4 + indt) + grad12 - grad22

            do ii = 1, indtm(i)
                if (tabpip(i, ii) .ne. 0.d0) then
                    zi = rkel(3, i, ii)
                    yi = rkel(2, i, ii)
                    xi = rkel(1, i, ii)

                    if (LBox .gt. 0.d0) then
                        rc(1) = xi - xj
                        rc(2) = yi - yj
                        rc(3) = zi - zj
                        !              call CartesianToCrystal(rc,1)
                        rc(:) = car2cry(:, 1)*rc(1) + car2cry(:, 2)*rc(2) + car2cry(:, 3)*rc(3)
                        rc(1) = map(rc(1), cellscale(1))
                        rc(2) = map(rc(2), cellscale(2))
                        rc(3) = map(rc(3), cellscale(3))
                    else
                        rc(1) = xi - xj
                        rc(2) = yi - yj
                        rc(3) = zi - zj
                    end if

                    !     old 3-body grad-lap
                    if (nelorbj .ne. 0) then
                        costj = psip(indtnew + ii, i)
                    else
                        costj = 0.d0
                    end if
                    costjnew = jastrow_ee(rc, vj, iesd, ispin) + costj

                    tabpip(i, ii) = tabpip(i, ii)*dexp(costjnew - &
                            & jastrowall_ee(jel, i, ii))*cost
                    !         update jastrowall
                    jastrowall_ee(jel, i, ii) = costjnew

                end if
            end do
            ! i ne jel
        end if
    end do
    !         update winvjbar
!$omp end parallel do
    !        compute by scratch  grad e grad^2
    i = jel
    do k = 1, 4
        costtab(k) = 0.d0
    end do
    if (nelorbj .ne. 0) then
        do j = 1, nel
            if (j .ne. i) then
                costtab(1:3) = costtab(1:3) &
                        & - jastrowall_ee(jel, j, indt + 1:indt + 3)                            &
                        & + psip(indtnew + indt + 1:indtnew + indt + 3, j) + psip(indtmn + 1:indtmn + 3, j)
                costtab(4) = costtab(4) &
                        & + jastrowall_ee(jel, j, indt + 4) - psip(indtnew + indt + 4, j) + psip(indtmn + 4, j)
            end if
        end do
    else
        do j = 1, nel
            if (j .ne. i) then
                costtab(1:3) = costtab(1:3) &
                        & - jastrowall_ee(jel, j, indt + 1:indt + 3)
                costtab(4) = costtab(4) + jastrowall_ee(jel, j, indt + 4)
            end if
        end do
    end if
    do k = 1, 4
        tabpip(i, indt + k) = costtab(k)
    end do

    if (n_body_on .ne. 0) then
        do k = 1, 4
            costtab(k) = tabpip(i, k + indt)
        end do
!$omp parallel do reduction(+:costtab,psiln) default(shared) &
!$omp  private(j,zi,yi,xi,rc,costz0,jastrow1,grad1,grad12)
        do j = 1, nion
            zi = rion(3, j)
            yi = rion(2, j)
            xi = rion(1, j)

            rc(1) = xj - xi
            rc(2) = yj - yi
            rc(3) = zj - zi

            costz0 = costz(j)*costz3(j)

            if (Lbox .gt. 0.d0) then
                call jastrowgrad_pbc(rc, vj(pointvj(1, j)), pointvj(2, j) &
                                     , jastrow1, grad1, grad12, -1, costz(j))
                psiln = psiln - costz3(j)*(jastrow1 - jasnew_ei(j))
                jasnew_ei(j) = jastrow1

            else
                rc(:) = costz(j)*rc(:)
                call jastrowgrad(rc, vj(pointvj(1, j)), pointvj(2, j) &
                                 , jastrow1, grad1, grad12, -1)
                psiln = psiln - costz3(j)*(jastrow1 - jasnew_ei(j))
                jasnew_ei(j) = jastrow1
            end if

            !        tabpip(i,1+indt)=tabpip(i,1+indt)-grad1(1)*costz0
            !        tabpip(i,2+indt)=tabpip(i,2+indt)-grad1(2)*costz0
            !        tabpip(i,3+indt)=tabpip(i,3+indt)-grad1(3)*costz0
            !        tabpip(i,4+indt)=tabpip(i,4+indt)-grad12*costz0*costz(j)
            costtab(1) = costtab(1) - grad1(1)*costz0
            costtab(2) = costtab(2) - grad1(2)*costz0
            costtab(3) = costtab(3) - grad1(3)*costz0
            costtab(4) = costtab(4) - grad12*costz0*costz(j)
        end do
!$omp end parallel do
        do k = 1, 4
            tabpip(i, k + indt) = costtab(k)
        end do

    end if

!$omp parallel do default(shared) &
!$omp  private(ii,kk,j,zjn,yjn,xjn,zi,yi,xi,ispin,rcn,costz0,costjnew,r0)

    do ii = 1, indtm(jel)
        tabpip(i, ii) = tcost(ii)

        if (tcost(ii) .ne. 0.d0) then
            zjn = rkeln(3, ii)
            yjn = rkeln(2, ii)
            xjn = rkeln(1, ii)

            costz0 = 0.d0
            do j = 1, nel
                if (j .ne. i) then
                    if (iesspin) then
                        if (j .le. nelup) then
                            ispin = ispinj
                        else
                            ispin = -ispinj
                        end if
                    else
                        ispin = -1
                    end if
                    zi = rkel(3, j, 0)
                    yi = rkel(2, j, 0)
                    xi = rkel(1, j, 0)

                    if (LBox .gt. 0.d0) then
                        rcn(1) = xi - xjn
                        rcn(2) = yi - yjn
                        rcn(3) = zi - zjn
                        !                call CartesianToCrystal(rcn,1)
                        rcn(:) = car2cry(:, 1)*rcn(1) + car2cry(:, 2)*rcn(2) + car2cry(:, 3)*rcn(3)
                        rcn(1) = map(rcn(1), cellscale(1))
                        rcn(2) = map(rcn(2), cellscale(2))
                        rcn(3) = map(rcn(3), cellscale(3))

                    else
                        rcn(1) = xi - xjn
                        rcn(2) = yi - yjn
                        rcn(3) = zi - zjn
                    end if

                    if (nelorbj .ne. 0) then
                        costjnew = jastrow_ee(rcn, vj, iesd, ispin) + psip(ii, j)
                        costz0 = costz0 + costjnew - jastrowall_ee(j, jel, 0)

                        jastrowall_ee(j, jel, ii) = costjnew
                    else
                        costjnew = jastrow_ee(rcn, vj, iesd, ispin)
                        costz0 = costz0 + costjnew - jastrowall_ee(j, jel, 0)
                        jastrowall_ee(j, jel, ii) = costjnew
                    end if
                end if
            end do
            tabpip(i, ii) = tabpip(i, ii)*dexp(costz0)

            if (n_body_on .ne. 0) then

                costz0 = 0.d0
                do j = 1, nion

                    if (Lbox .gt. 0.d0) then
                        rcn(:) = rkeln(:, ii) - rion(:, j)
!                       call CartesianToCrystal(rcn, 1)
                        rcn(:) = car2cry(:, 1)*rcn(1) + car2cry(:, 2)*rcn(2) + car2cry(:, 3)*rcn(3)
                        do kk = 1, 3
                            rcn(kk) = costz(j)*map(rcn(kk), cellscale(kk))
                        end do
                        r0 = norm_metric(rcn, metric)
                    else
                        r0 = costz(j)*dsqrt(sum(rmu(:, ii, j)**2))
                    end if

                    costz0 = costz0 + &
                            &costz3(j)*(jasnew_ei(j) - jastrow_ei(r0, vj(pointvj(1, j)), pointvj(2, j)))
                end do
                tabpip(i, ii) = tabpip(i, ii)*dexp(costz0)
                ! endif iesdr<= -5
            end if
            ! endif tcost=0
        end if
        ! enddo on ii=1,indt
    end do
!$omp end parallel do

    if (nelorbj .ne. 0) then
#ifdef _OFFLOAD
!$omp target teams distribute  parallel do if(yes_ontarget)
#endif
        !      winvjbar(1:nelorbjh*ipj,jel)=winvjbarn(1:nelorbjh*ipj)
        do k = 1, nelorbjh*ipj
            winvjbar(k, jel) = winvjbarn(k)
        end do

        if (iessz) winvjbarsz(indjel:indjel + nelorbjh - 1) = winvjbarszn(1:nelorbjh)
    end if
    !       write(*,*) " winvjbarn in uptappip: ",sum(winvjbarn(1:nelorbjh*ipj)),sum(winvjbar(:,:))
    do i = 1, nel
        do j = indtm(i) + 1, indt
            tabpip(i, j) = 0.d0
        end do
    end do

    return
end
subroutine copy_simple(n, n_tra, m, mat, mat_tra, yes_ontarget)
    implicit none
    integer n, n_tra, m, i, j
    real*8 mat(n, m), mat_tra(n_tra, m)
    logical yes_ontarget
#ifdef _OFFLOAD
    if (yes_ontarget) then
!$omp target teams distribute parallel do collapse(2)
        do j = 1, m
            do i = 1, n
                mat_tra(i, j) = mat(i, j)
            end do
        end do
    else
        do j = 1, m
            do i = 1, n
                mat_tra(i, j) = mat(i, j)
            end do
        end do
    end if
#else
    do j = 1, m
        do i = 1, n
            mat_tra(i, j) = mat(i, j)
        end do
    end do
#endif
end subroutine copy_simple
