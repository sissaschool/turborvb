!TL off
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

subroutine upvpotdiag(rkel, nel, zeta, rion, iond, vpot, nion           &
        &, kappa, LBox, vpotreg, cutreg, costz, costz3,vpotsav_ee,yes_fn)

    ! by E. Coccia (30/11/10)
    use extpot, only : ext_pot, link_atom, mm_restr
    ! by E. Coccia (4/2/11)
    use van_der_waals, only : vdw
    ! by E. Coccia (12/5/11)
    use link_angle, only : mm_pot_theta, mm_pot_dihed, mm_pot_impr
    ! by E. Coccia (7/12/11)
    use cl_restr, only : restr_bond, restr_angle, restr_dihe, restr_dimp
    use dielectric
    use cell, only: x_neigh,neigh,dist_shift,distreg_shift
    use allio, only: iond_cart,pot_aasunel
    implicit none
    integer nel,ii, jj, i, j, k, nion
    integer*8 maxll,ll
    real(8) rkel(3, *)
    real(8) zeta(nion), vpot, rion(3, *), iond(nion, nion)   &
            &, cost, tmp, kappa, LBox, pot_ee, pot_ea, pot_aa, derfc, max, costz(*)&
            &, distb, vpotreg(2, nel), potreg_ea, costz3(*), cutreg, Zreg2, distbreg
    real(8), external :: ngivej
    real(8) x_shift(3),cost_z
    logical yes_fn
    real(8) vpotsav_ee(nel,nel)

    ! by E. Coccia (10/11/10)
    real*8 ext_ene(nel)
    real*8 sum_ext
    ! by E. Coccia (12/1/11)
    real*8 ext_ion(nion)
    real*8 sum_ion
    ! by E. Coccia (4/2/11)
    real*8 sum_vdw

    vpot = 0.d0
    pot_aa = 0.d0
    vpotreg = 0.d0
    ! by E. Coccia
    ext_ene = 0.d0
    ext_ion = 0.d0

    if(LBox.le.0.d0) then
!$omp parallel do default(shared)  private(k,i,Zreg2,distb) reduction(+:vpotreg)
        do k = 1, nion
            Zreg2 = costz(k) * costz3(k)
            do i = 1, nel
                distb = ngivej(rion(1, k), rkel(1, i), LBox)
                vpotreg(1, i) = vpotreg(1, i) - 2.d0 * zeta(k) *veps(distb)
                !        pot_ea=pot_ea-2.d0*zeta(k)/distb
                !        pot_ea=pot_ea-2.d0*zeta(k)/distb
                !        potreg_ea=potreg_ea+Zreg2/distb-Zreg2/max(distb,cutreg)
                vpotreg(2, i) = vpotreg(2, i) + (Zreg2 - 2.d0 * zeta(k)) *veps(distb) - Zreg2 * veps(max(distb, cutreg))
            enddo
        enddo

!$omp parallel do default(shared)  private(j,i,cost) reduction(+:vpotreg)
        do i = 1, nel
            do j = i+1, nel
                    cost =veps( ngivej(rkel(1, i), rkel(1, j), LBox))
                    !         pot_ee=pot_ee+cost
                    vpotreg(1, j) = vpotreg(1, j) + cost
                    vpotreg(2, j) = vpotreg(2, j) + cost
                    vpotreg(1, i) = vpotreg(1, i) + cost
                    vpotreg(2, i) = vpotreg(2, i) + cost
                    if(yes_fn) then
                    vpotsav_ee(i,j)=2.d0*cost
                    vpotsav_ee(j,i)=2.d0*cost
                    endif
            enddo
        enddo

!$omp parallel do default(shared) private(i,j) reduction(+:pot_aa)
        do i = 1, nion
            do j = i + 1, nion
                if(zeta(i) * zeta(j).ne.0.d0) pot_aa = pot_aa + 2.d0 * zeta(i) * zeta(j) * veps(iond(i, j))
            enddo
        enddo

        ! ********  EWALD REAL PART *************
    else

        ! write(6,*) cutreg
        ! write(6,*) costz(1:nion)
        ! write(6,*) costz3(1:nion)
        ! write(6,*) zeta(1:nion)

!$omp parallel do default(shared)  private(k,i,ii,Zreg2,x_shift,cost,dist_shift,distreg_shift) reduction(+:vpotreg)
        do k = 1, nion
            Zreg2 = costz(k) * costz3(k)  ! 2 Z_singular
            if(Zreg2 * cutreg.gt.0) then
                do i = 1, nel
                   call dgivej(rion(1, k), rkel(1, i), LBox,x_shift)
#ifdef _SIMD
!$omp simd
#endif
                   do ii=1,neigh  
                    dist_shift(ii) = max(dsqrt((x_shift(1)+x_neigh(ii,1))**2+&
      & (x_shift(2)+x_neigh(ii,2))**2+(x_shift(3)+x_neigh(ii,3))**2),1d-9) 
                    distreg_shift(ii) = max(dist_shift(ii), cutreg)
!                   cost = derfc(kappa * distb) / distb
                    dist_shift(ii)=rep_erfc(dist_shift(ii),kappa)
                    !        pot_ea=pot_ea-2.d0*zeta(k)*cost
                    distreg_shift(ii)=rep_erfc(distreg_shift(ii),kappa)
                   enddo
                   do ii=1,neigh
               vpotreg(1, i) = vpotreg(1, i) - 2.d0 * zeta(k) * dist_shift(ii) 
                    !        potreg_ea=potreg_ea+Zreg2*(cost-derfc(kappa*distbreg)/distbreg)
      vpotreg(2, i) = vpotreg(2, i) + (Zreg2 - 2.d0 * zeta(k)) * dist_shift(ii)&
!                           & - Zreg2 * derfc(kappa * distbreg) / distbreg
                            & - Zreg2 * distreg_shift(ii)

                   enddo
                enddo
            else
                do i = 1, nel
                   call dgivej(rion(1, k), rkel(1, i), LBox,x_shift)
#ifdef _SIMD
!$omp simd
#endif
                   do ii=1,neigh
                    dist_shift(ii) = max(dsqrt((x_shift(1)+x_neigh(ii,1))**2+&
      & (x_shift(2)+x_neigh(ii,2))**2+(x_shift(3)+x_neigh(ii,3))**2),1d-9) 
!                   distb = ngivej(rion(1, k), rkel(1, i), LBox)
!                   cost = -2.d0 * zeta(k) * derfc(kappa * distb) / distb
                    dist_shift(ii)=-2.d0*zeta(k)*rep_erfc(dist_shift(ii),kappa)
                    !         pot_ea=pot_ea+cost
                   enddo
                   cost=sum(dist_shift(1:neigh))
                   vpotreg(1, i) = vpotreg(1, i) + cost 
                   vpotreg(2, i) = vpotreg(2, i) + cost 
                enddo
            endif
        enddo

        !      potreg_ea=pot_ea+potreg_ea
        maxll=nel
        maxll=(maxll*(nel-1))/2
!$omp parallel do default(shared)  private(ll,j,i,ii,cost,dist_shift,x_shift) reduction(+:vpotreg)
        do ll=1,maxll
        call find_j1j2(nel,ll,i,j)
!       do i = 1, nel
!           do j = i+1, nel
!               if(j.ne.i) then
!                   tmp = ngivej(rkel(1, i), rkel(1, j), LBox)
                    call dgivej(rkel(1, i), rkel(1, j), LBox,x_shift)
!                   cost = derfc(tmp * kappa) / tmp
#ifdef _SIMD
!$omp simd
#endif
                    do ii=1,neigh
                    dist_shift(ii) = max(dsqrt((x_shift(1)+x_neigh(ii,1))**2+&
      & (x_shift(2)+x_neigh(ii,2))**2+(x_shift(3)+x_neigh(ii,3))**2),1d-9) 

                    dist_shift(ii)=rep_erfc(dist_shift(ii),kappa)
                    !          pot_ee=pot_ee+cost
!                   vpotreg(1, j) = vpotreg(1, j) + cost
!                   vpotreg(2, j) = vpotreg(2, j) + cost
                    enddo
                    cost=sum(dist_shift(1:neigh))
                    vpotreg(1,j)=vpotreg(1,j)+cost
                    vpotreg(2,j)=vpotreg(2,j)+cost
                    vpotreg(1,i)=vpotreg(1,i)+cost
                    vpotreg(2,i)=vpotreg(2,i)+cost
                    if(yes_fn) then
                    vpotsav_ee(i,j)=2.d0*cost
                    vpotsav_ee(j,i)=2.d0*cost
                    endif
!               endif
!           enddo
        enddo
        maxll=nion
        maxll=(maxll*(nion-1))/2
!$omp parallel do default(shared)  private(ll,j,i,ii,jj,cost_z,dist_shift) reduction(+:pot_aa)
        do ll=1,maxll
        call find_j1j2(nion,ll,i,j)
!       do i = 1, nion
!           do j = i + 1, nion
            jj=nion*(j-1)+i
            cost_z=2.d0*zeta(i)*zeta(j)
                if(cost_z.ne.0.d0)   then 
#ifdef _SIMD
!$omp simd
#endif
                   do ii=1,neigh
   dist_shift(ii)=dsqrt((iond_cart(1,jj)+x_neigh(ii,1))**2+(iond_cart(2,jj)+x_neigh(ii,2))**2+(iond_cart(3,jj)+x_neigh(ii,3))**2)
                dist_shift(ii)=cost_z*rep_erfc(dist_shift(ii),kappa)
!                   pot_aa = pot_aa + cost_z*rep_erfc(dist_shift(ii),kappa)
                   enddo
                pot_aa=pot_aa+sum(dist_shift(1:neigh))
                endif
            enddo
!       enddo
    endif

    pot_aasunel=pot_aa/nel
    vpot = 0.d0
    do i = 1, nel
        vpotreg(1:2, i) = pot_aa / nel + vpotreg(1:2, i)
        vpot = vpot + vpotreg(1, i)
    enddo

    ! by E. Coccia (10/11/10): add an external potential (electronic)
    if (ext_pot) then
        call extpot_ene(rkel, nel, ext_ene, sum_ext)
        do i = 1, nel
            vpotreg(:, i) = vpotreg(:, i) + 2.d0 * ext_ene(i)
        enddo
        vpot = vpot + 2.d0 * sum_ext
        ! by E. Coccia (12/1/11): add an external potential (nuclear)
        call ion_ene(rion, zeta, nion, ext_ion, sum_ion)
        do i = 1, nel
            vpotreg(:, i) = vpotreg(:, i) + 2.d0 * sum_ion / real(nel)
        enddo
        vpot = vpot + 2.d0 * sum_ion
        ! by E. Coccia (4/2/11): add the vdw term
        if (vdw) then
            call vdw_ene(rion, nion, sum_vdw)
            do i = 1, nel
                vpotreg(:, i) = vpotreg(:, i) + 2.d0 * sum_vdw / real(nel)
            enddo
            vpot = vpot + 2.d0 * sum_vdw
        endif
        ! by E. Coccia (11/5/11): add the mm angular term in case of QMC/MM boundary
        if (link_atom) then
            call mm_angle()
            do i = 1, nel
                vpotreg(:, i) = vpotreg(:, i) + 2.d0 * mm_pot_theta / real(nel)
            enddo
            vpot = vpot + 2.d0 * mm_pot_theta
            call mm_dihed()
            do i = 1, nel
                vpotreg(:, i) = vpotreg(:, i) + 2.d0 * mm_pot_dihed / real(nel)
            enddo
            vpot = vpot + 2.d0 * mm_pot_dihed
            call mm_improper()
            do i = 1, nel
                vpotreg(:, i) = vpotreg(:, i) + 2.d0 * mm_pot_impr / real(nel)
            enddo
            vpot = vpot + 2.d0 * mm_pot_impr
        endif
    endif

    ! by E. Coccia (7/12/11): classical restraints
    if (mm_restr) then
        call r_bond()
        do i = 1, nel
            vpotreg(:, i) = vpotreg(:, i) + 2.d0 * restr_bond / real(nel)
        enddo
        vpot = vpot + 2.d0 * restr_bond
        call r_angle()
        do i = 1, nel
            vpotreg(:, i) = vpotreg(:, i) + 2.d0 * restr_angle / real(nel)
        enddo
        vpot = vpot + 2.d0 * restr_angle
        call r_dihe()
        do i = 1, nel
            vpotreg(:, i) = vpotreg(:, i) + 2.d0 * restr_dihe / real(nel)
        enddo
        vpot = vpot + 2.d0 * restr_dihe
        call r_dimp()
        do i = 1, nel
            vpotreg(:, i) = vpotreg(:, i) + 2.d0 * restr_dimp / real(nel)
        enddo
        vpot = vpot + 2.d0 * restr_dimp

    endif

    return

END subroutine upvpotdiag
