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

!=======================================================================
! This module contains all the data structures and subrotuines
! needed for performing the Ewald summation of the long-range Coulomb
! potential. Its need information about the supercell, and particles.
!=======================================================================

module Ewald
    use Constants
    use Cell, only: cellscale, recip, lmin, omega, celldm
    use dielectric, only: epsilon0
    implicit none
    double precision ksq
    double precision, allocatable :: q(:)
    integer kmax
    integer natoms, nelectrons
    ! tot num of charges including virtual ones for threading
    integer :: ncharges
    ! tot num of physical charges including ions and electrons
    integer :: n_physical_charges
    double precision :: kappa ! cut off in real space
    double precision :: ecut ! cut off in rec. space
    double precision :: gcut ! |G| max
    double precision :: eself, eself1b, ewaldion1b, ewaldel1b ! Ewald self-energy
    integer :: n_gvec ! number of plane waves
    integer :: nwalkers
    double precision, dimension(:, :), allocatable :: g ! G vectors
    double precision, dimension(:), allocatable :: gg ! |G|^2
    double precision, dimension(:), allocatable :: factor
    ! exp(-G^2/4 kappa^2)/G^2
    double complex, dimension(:, :), allocatable :: phsfac1, phsfac2, phsfac3
    double complex, dimension(:, :), allocatable :: twiddle
    double precision, dimension(:, :, :), allocatable :: psip_sin, psip_cos
    integer, dimension(3) :: nr
    integer, dimension(:, :), allocatable :: ind_g ! (h,k,l)
    double complex, dimension(:), allocatable :: qphase
    double precision, dimension(:), allocatable :: q_cos_gr, q_sin_gr&
            &, sum_q_cos, sum_q_sin
    double precision, dimension(:, :), allocatable :: s_cord
    double precision, dimension(:, :), allocatable :: sum_q_cos_gr&
            &, sum_q_sin_gr

    integer movedions
    ! flag to not recalculate the ionic part of ewald
contains

    subroutine InitEwald(nion, zetar, nel, nw, nthreads)
        implicit none
        integer, intent(in) :: nion, nel, nw
        double precision, intent(in) :: zetar(nion)
        ! should be effective only when nthreads > nel
        integer, intent(in), optional :: nthreads
        real*8 sab, sac, sbc, cond
        real*8, external :: cond_find
        if (allocated(q)) then
            deallocate (sum_q_cos_gr, sum_q_sin_gr &
                    &, q, s_cord, sum_q_cos, sum_q_sin)
        end if
        if (ksq .eq. 0.5d0) then
            gcut = dsqrt(2.d0*ksq)*kappa**2*lmin
            ecut = gcut**2
        else
            sab = sum(recip(:, 1)*recip(:, 2))/dsqrt(sum(recip(:, 1)**2)*sum(recip(:, 2)**2))
            sbc = sum(recip(:, 3)*recip(:, 2))/dsqrt(sum(recip(:, 3)**2)*sum(recip(:, 2)**2))
            sac = sum(recip(:, 3)*recip(:, 1))/dsqrt(sum(recip(:, 3)**2)*sum(recip(:, 1)**2))
            cond = cond_find(sab, sac, sbc)

!       sab=dsqrt(abs(1.d0-sab**2))
!       sbc=dsqrt(abs(1.d0-sbc**2))
!       sac=dsqrt(abs(1.d0-sac**2))

            ecut = -4.d0*kappa**2*log(ksq)
!       gcut=dsqrt(ecut)/cond
!  S. Sorella,  the above conditions are sufficient to have
!  a sphere of radius sqrt(ecut) inside the cell.
            gcut = dsqrt(ecut)/cond

        end if

        natoms = nion
        nelectrons = nel
        nwalkers = nw
        ncharges = natoms + nelectrons
        n_physical_charges = ncharges
        if (present(nthreads)) then
            if (nthreads > nel) then
                ncharges = ncharges + nthreads - nel
            end if
        end if

        allocate (q(ncharges), s_cord(3, ncharges))
        q(1:natoms) = zetar(1:natoms)
        q(natoms + 1:ncharges) = -1.d0
        call Ggen
        allocate (sum_q_cos_gr(n_gvec + 1, nwalkers), sum_q_sin_gr(n_gvec, nwalkers)&
                &, sum_q_cos(n_gvec), sum_q_sin(n_gvec))

        call EwaldSelf

    end subroutine InitEwald

    subroutine Ggen_samegrid
        double precision :: b, gv(3), g2
        integer :: i, j, k, ig
        g = 0.d0
        gg = 0.d0
        factor = 0.d0
        do ig = 1, n_gvec
            i = ind_g(ig, 1); j = ind_g(ig, 2); k = ind_g(ig, 3)

            gv(:) = dble(i)*recip(1, :) + dble(j)*recip(2, :) + dble(k)*recip(3, :)
            g2 = dot_product(gv, gv)
            g(:, ig) = gv
            gg(ig) = g2
            factor(ig) = epsilon0*dexp(-g2/(4.d0*kappa**2.d0))/g2
        end do
    end subroutine Ggen_samegrid

    subroutine Ggen
        double precision :: b, gv(3), g2
        integer :: i, j, k, j0, i0
        integer nr_max(3)

        do i = 1, 3
            b = dsqrt(sum(recip(i, :)**2.d0))
            if (b .ne. 0.d0) then
                nr(i) = int(gcut/b) + 1
            else
                nr(i) = 0.d0
            end if
        end do

        if (allocated(qphase)) deallocate (qphase, q_cos_gr, q_sin_gr, twiddle&
                &, phsfac1, phsfac2, phsfac3)

        allocate (twiddle(ncharges, 3))
        allocate (qphase(ncharges), q_cos_gr(ncharges), q_sin_gr(ncharges))

        ! count G vectors
        n_gvec = 0
        nr_max = 0
        do k = 0, nr(3)
            j0 = -nr(2); if (k == 0) j0 = 0
            do j = j0, nr(2)
                i0 = -nr(1); if ((k == 0) .and. (j == 0)) i0 = 1
                do i = i0, nr(1)
                    gv(:) = dble(i)*recip(1, :) + dble(j)*recip(2, :) + dble(k)*recip(3, :)

                    g2 = dot_product(gv, gv)
                    if (g2 < ecut) then
                        n_gvec = n_gvec + 1
                        nr_max(1) = max(nr_max(1), abs(i))
                        nr_max(2) = max(nr_max(2), abs(j))
                        nr_max(3) = max(nr_max(3), k)
                    end if
                end do
            end do
        end do
        nr(:) = nr_max(:)

        allocate (phsfac1(-nr(1):nr(1), ncharges))
        allocate (phsfac2(-nr(2):nr(2), ncharges))
        allocate (phsfac3(0:nr(3), ncharges)) ! see in TrigRecur

        if (n_gvec .eq. 0) n_gvec = 1
        if (allocated(g)) deallocate (g, gg, factor, ind_g)
        allocate (g(3, n_gvec), gg(n_gvec), factor(n_gvec), ind_g(n_gvec, 3))
        n_gvec = 0
        g = 0.d0
        gg = 0.d0
        factor = 0.d0
        ind_g = 0
        do k = 0, nr(3)
            j0 = -nr(2); if (k == 0) j0 = 0
            do j = j0, nr(2)
                i0 = -nr(1); if ((k == 0) .and. (j == 0)) i0 = 1
                do i = i0, nr(1)
                    gv(:) = dble(i)*recip(1, :) + dble(j)*recip(2, :) + dble(k)*recip(3, :)
                    g2 = dot_product(gv, gv)
                    if (g2 < ecut) then
                        n_gvec = n_gvec + 1
                        g(:, n_gvec) = gv
                        ind_g(n_gvec, 1) = i
                        ind_g(n_gvec, 2) = j
                        ind_g(n_gvec, 3) = k
                        gg(n_gvec) = g2
                        if (i .ne. 0 .or. j .ne. 0 .or. k .ne. 0) then
                            factor(n_gvec) = epsilon0*dexp(-g2/(4.d0*kappa**2))/g2
                        else
                            factor(n_gvec) = 0.d0
                        end if
                    end if
                end do
            end do
        end do
        if (n_gvec .eq. 0) n_gvec = 1
    end subroutine Ggen

    !=====================================================================
    ! Calculate Ewald self-energy
    !=====================================================================
    subroutine EwaldSelf
        implicit none
        eself = -2.d0*epsilon0*kappa/dsqrt(PI)*sum(q(1:n_physical_charges)**2.d0)
        eself1b = -2.d0*epsilon0*kappa/dsqrt(PI)*sum(q(1:natoms)**2)/nelectrons
    end subroutine EwaldSelf

    !=====================================================================
    ! Calculate by recursion, the sin and cos of G.r, where G is a
    ! reciprocal lattice vector, and r is an atomic position
    !=====================================================================

    subroutine TrigRecur(particle_start, particle_end, s)
        implicit none
        integer, intent(in) :: particle_start, particle_end
        double precision, intent(in) :: s(3, particle_end)
        double precision :: g_dot_r
        !           complex*16 :: twiddle
        integer :: icharge, j

        if (particle_start .eq. particle_end .and. old_threads .gt. 1) then
            phsfac1(0, particle_start) = (1.d0, 0.d0) ! h = 0
            g_dot_r = sum(recip(1, :)*s(:, particle_start))
            twiddle(particle_start, 1) = dcmplx(dcos(g_dot_r), dsin(g_dot_r))

            phsfac2(0, particle_start) = (1.d0, 0.d0) ! k = 0
            g_dot_r = sum(recip(2, :)*s(:, particle_start))
            twiddle(particle_start, 2) = dcmplx(dcos(g_dot_r), dsin(g_dot_r))

            phsfac3(0, particle_start) = (1.d0, 0.d0) ! l = 0
            g_dot_r = sum(recip(3, :)*s(:, particle_start))
            twiddle(particle_start, 3) = dcmplx(dcos(g_dot_r), dsin(g_dot_r))
!$omp parallel do default(shared) private(j)
            do j = 1, nr(1)
                phsfac1(j, particle_start) = twiddle(particle_start, 1)**j
                phsfac1(-j, particle_start) = dconjg(phsfac1(j, particle_start))
            end do
!$omp parallel do default(shared) private(j)
            do j = 1, nr(2)
                phsfac2(j, particle_start) = twiddle(particle_start, 2)**j
                phsfac2(-j, particle_start) = dconjg(phsfac2(j, particle_start))
            end do
            ! In order to save memory when doing long cells in the z direction
            ! here only for l >= 0
!$omp parallel do default(shared) private(j)
            do j = 1, nr(3)
                phsfac3(j, particle_start) = twiddle(particle_start, 3)**j
            end do
        else

!$omp parallel do default(shared) private(icharge,j,g_dot_r)
            do icharge = particle_start, particle_end
                phsfac1(0, icharge) = (1.d0, 0.d0) ! h = 0

                g_dot_r = sum(recip(1, :)*s(:, icharge))
                twiddle(icharge, 1) = dcmplx(dcos(g_dot_r), dsin(g_dot_r))
                phsfac1(1, icharge) = twiddle(icharge, 1) ! h = 1
                phsfac1(-1, icharge) = dconjg(twiddle(icharge, 1)) ! h =-1

                do j = 2, nr(1)
                    phsfac1(j, icharge) = phsfac1(j - 1, icharge)*twiddle(icharge, 1)
                    phsfac1(-j, icharge) = dconjg(phsfac1(j, icharge))
                end do

                phsfac2(0, icharge) = (1.d0, 0.d0) ! k = 0
                g_dot_r = sum(recip(2, :)*s(:, icharge))
                twiddle(icharge, 2) = dcmplx(dcos(g_dot_r), dsin(g_dot_r))
                phsfac2(1, icharge) = twiddle(icharge, 2) ! k = 1
                phsfac2(-1, icharge) = dconjg(twiddle(icharge, 2)) ! k =-1

                do j = 2, nr(2)
                    phsfac2(j, icharge) = phsfac2(j - 1, icharge)*twiddle(icharge, 2)
                    phsfac2(-j, icharge) = dconjg(phsfac2(j, icharge))
                end do

                ! In order to save memory when doing long cells in the z direction
                ! here only for l >= 0
                phsfac3(0, icharge) = (1.d0, 0.d0) ! l = 0
                g_dot_r = sum(recip(3, :)*s(:, icharge))
                twiddle(icharge, 3) = dcmplx(dcos(g_dot_r), dsin(g_dot_r))
                phsfac3(1, icharge) = twiddle(icharge, 3) ! l = 1
                do j = 2, nr(3)
                    phsfac3(j, icharge) = phsfac3(j - 1, icharge)*twiddle(icharge, 3)
                end do
            end do
        end if
    end subroutine TrigRecur

    !=====================================================================
    ! This compute the reciprocal part of the Ewald summation without
    ! computation of stress
    !=====================================================================

    subroutine EwaldSum(kel, rion, energy, vpotreg, iw)
        implicit none
        double precision, intent(in) :: rion(3, natoms)
        double precision, intent(in) :: kel(3, nelectrons)
        double precision, intent(out) :: energy
        double precision :: vpotreg(2, nelectrons), sumce, sumci, sumse&
                &, sumsi, cost, costii, fac0, sumcos, sumsin
        double complex :: cphase
        integer, intent(in) :: iw
        integer :: ig, h, k, l, i, p, j, start

        !           call CartesianToCrystal(rion,s_cord,natoms)
        !           call CartesianToCrystal(kel,s_cord(1,natoms+1),nelectrons)
        do i = 1, natoms
            s_cord(:, i) = rion(:, i)
        end do
        do i = natoms + 1, n_physical_charges
            s_cord(:, i) = kel(:, i - natoms)
        end do
        start = 1
        if (movedions == 0) start = natoms + 1

        !           do p=start,n_physical_charges
        call TrigRecur(start, n_physical_charges, s_cord)
        !           enddo

        !           fac0=8*pi/omega
        energy = 0.d0
        ! with half g-space, G=0 is not in the list
!$omp parallel do default(shared) reduction(+:energy) private(ig,i,h,k,l,sumcos,sumsin,cphase)
        do ig = 1, n_gvec
            h = ind_g(ig, 1); k = ind_g(ig, 2); l = ind_g(ig, 3)
            !       qphase(:) = q(:) * phsfac1(:,h)* phsfac2(:,k) * phsfac3(:,l)

            sumcos = 0.d0
            sumsin = 0.d0
            do i = 1, ncharges
                cphase = q(i)*phsfac1(h, i)*phsfac2(k, i)*phsfac3(l, i)
                sumcos = sumcos + dreal(cphase)
                sumsin = sumsin + dimag(cphase)
            end do
            !      sumce=sum(q_cos_gr(natoms+1:n_physical_charges))
            !      sumci=sum(q_cos_gr(1:natoms))

            !      sum_q_cos_gr(ig,iw) = sumce+sumci

            !      sumse=sum(q_sin_gr(natoms+1:n_physical_charges))
            !      sumsi=sum(q_sin_gr(1:natoms))

            !      sum_q_sin_gr(ig,iw) = sumse+sumsi

            sum_q_cos_gr(ig, iw) = sumcos
            sum_q_sin_gr(ig, iw) = sumsin

            energy = energy + factor(ig)*(sum_q_cos_gr(ig, iw)**2 + sum_q_sin_gr(ig, iw)**2)

            !     costii=factor(ig)*(sumsi**2+sumci**2)/nelectrons

            !     do j=1,nelectrons
            !     cost=costii+factor(ig)*(q_sin_gr(j+natoms)*(2.d0*sumsi+sumse)&
            !    &+q_cos_gr(j+natoms)*(2.d0*sumci+sumce))
            !     vpotreg(1:2,j)=vpotreg(1:2,j)+cost*fac0
            !     enddo

        end do

        !     fac0=eself/nelectrons
        !     do j=1,nelectrons
        !     vpotreg(:,j)=vpotreg(:,j)+fac0
        !     enddo

        energy = energy*8.d0*Pi/omega
        sum_q_cos_gr(n_gvec + 1, iw) = energy
        energy = energy + eself

        !   The long range Ewald contribution is non singular and is averaged over all
        !   the electrons in the regularized potentials
        do j = 1, nelectrons
            vpotreg(:, j) = vpotreg(:, j) + energy/nelectrons
        end do
    end subroutine EwaldSum

    !=====================================================================
    ! This compute the reciprocal part of the Ewald summation without
    ! computation of stress
    !=====================================================================
    subroutine EwaldSum_b(kel, kelb, rion, rionb, energy, energyb&
            &, cellscale, cellscaleb, recipb, omegab, iw)
        implicit none
        double precision, intent(in) :: rion(3, natoms)
        double precision, intent(in) :: kel(3, nelectrons)
        double precision, intent(out) :: energy
        double precision :: rionb(3, natoms), kelb(3, nelectrons), energyb, cost, sumcos, sumsin, q_cos, q_sin
        real*8 :: cellscale(3), cellscaleb(3), ck2, cost2, recipb(3, 3)
        real*8 :: gvb(3), omegab
        double complex :: cphase

        !           reverse algorithm  of EwaldSum
        !           direct
        !        see Claudio's thesis Eq.4.17 pag.53

        integer, intent(in) :: iw
        integer :: ig, h, k, l, p, start, i

        !           call CartesianToCrystal(rion,s_cord,natoms)
        !           call CartesianToCrystal(kel,s_cord(1,natoms+1),nelectrons)
        !#ifdef _OFFLOAD
        !!$omp target data map(alloc:phsfac1,phsfac2,phsfac3,q_cos_gr,q_sin_gr,qphase)
        !#endif

        s_cord(:, 1:natoms) = rion(:, 1:natoms)
        s_cord(:, natoms + 1:n_physical_charges) = kel(:, 1:nelectrons)

        start = 1
        if (movedions == 0) start = natoms + 1

        !           do p=start,n_physical_charges
        call TrigRecur(start, n_physical_charges, s_cord)
        !           enddo

        ck2 = 1.d0/(4.d0*kappa**2)

        energy = 0.d0
        ! with half g-space, G=0 is not in the list
!$omp parallel do default(shared) private(ig,h,k,l,i,cost&
!$omp& ,cost2,gvb,cphase,sumcos,sumsin,q_cos,q_sin) &
!$omp& reduction(+:energy,rionb,kelb,recipb)
        do ig = 1, n_gvec
            h = ind_g(ig, 1); k = ind_g(ig, 2); l = ind_g(ig, 3)
!           qphase(:) = q(:) * phsfac1(h, :) * phsfac2(k, :) * phsfac3(l, :)

            cost = factor(ig)*energyb*16.d0*pi/omega
            sumcos = 0.d0
            sumsin = 0.d0
            do i = 1, ncharges
                cphase = q(i)*phsfac1(h, i)*phsfac2(k, i)*phsfac3(l, i)
                sumcos = sumcos + dreal(cphase)
                sumsin = sumsin + dimag(cphase)
            end do

            sum_q_cos_gr(ig, iw) = sumcos
            sum_q_sin_gr(ig, iw) = sumsin

!           q_cos_gr = dreal(qphase)
!           q_sin_gr = dimag(qphase)

!           sum_q_cos_gr(ig, iw) = sum(q_cos_gr)
!           sum_q_sin_gr(ig, iw) = sum(q_sin_gr)

            cost2 = sum_q_cos_gr(ig, iw)**2 + sum_q_sin_gr(ig, iw)**2

            energy = energy + factor(ig)*cost2

            !     derivative of the energy with respect to input rion,kel

            !         derivative of factor
            gvb(:) = -cost*cost2*(ck2 + 1.d0/gg(ig))*g(:, ig)
            recipb(1, :) = recipb(1, :) + gvb(:)*h
            recipb(2, :) = recipb(2, :) + gvb(:)*k
            recipb(3, :) = recipb(3, :) + gvb(:)*l
            !   derivative of the sum_cos and sum_sin
            do p = 1, natoms
                cphase = q(p)*phsfac1(h, p)*phsfac2(k, p)*phsfac3(l, p)
                q_cos = dreal(cphase)
                q_sin = dimag(cphase)
                rionb(:, p) = rionb(:, p) + cost*(sum_q_sin_gr(ig, iw)*q_cos - &
                        & sum_q_cos_gr(ig, iw)*q_sin)*g(:, ig)
                recipb(1, :) = recipb(1, :) - cost*h*rion(:, p)* &
                        &(sum_q_cos_gr(ig, iw)*q_sin - sum_q_sin_gr(ig, iw)*q_cos)
                recipb(2, :) = recipb(2, :) - cost*k*rion(:, p)* &
                        &(sum_q_cos_gr(ig, iw)*q_sin - sum_q_sin_gr(ig, iw)*q_cos)
                recipb(3, :) = recipb(3, :) - cost*l*rion(:, p)* &
                        &(sum_q_cos_gr(ig, iw)*q_sin - sum_q_sin_gr(ig, iw)*q_cos)
            end do
            do p = natoms + 1, n_physical_charges
                cphase = q(p)*phsfac1(h, p)*phsfac2(k, p)*phsfac3(l, p)
                q_cos = dreal(cphase)
                q_sin = dimag(cphase)
                kelb(:, p - natoms) = kelb(:, p - natoms) + cost*(sum_q_sin_gr(ig, iw)*q_cos - &
                          & sum_q_cos_gr(ig, iw)*q_sin)*g(:, ig)
                recipb(1, :) = recipb(1, :) - cost*h*kel(:, p - natoms)* &
                        &(sum_q_cos_gr(ig, iw)*q_sin - sum_q_sin_gr(ig, iw)*q_cos)
                recipb(2, :) = recipb(2, :) - cost*k*kel(:, p - natoms)* &
                        &(sum_q_cos_gr(ig, iw)*q_sin - sum_q_sin_gr(ig, iw)*q_cos)
                recipb(3, :) = recipb(3, :) - cost*l*kel(:, p - natoms)* &
                        &(sum_q_cos_gr(ig, iw)*q_sin - sum_q_sin_gr(ig, iw)*q_cos)
            end do

            !  We neglect the implicit dependence of kappa on cellscale. This is OK,
            !  since ewald sum do not depend on the kappa chosen, when the sum are
            !  converged as we suppose is the case.

        end do
!$omp end parallel do
        !  The statement below should be at the beginning but it can be commuted.
        !       reverse of
        !        energy=energy*8.d0*Pi/omega+eself
        omegab = omegab - energy*8.d0*PI/omega**2*energyb

        energyb = 0_8
        !#ifdef _OFFLOAD
        !!$omp end target data
        !#endif

    end subroutine EwaldSum_b

    !=====================================================================
    ! This compute the reciprocal part of the Ewald summation without
    ! computation of stress
    !=====================================================================
    subroutine EwaldUpdate(kel, keln, jel, energy, denergy, iw)
        implicit none
        double precision, intent(in) :: kel(3), keln(3)
        double complex :: qp
        integer jelnatoms
        integer, intent(in) :: jel, iw
        double precision, intent(out) :: energy, denergy
        integer :: ig, h, k, l
        double precision :: old_energy, new_energy
        new_energy = 0.d0
        !           old_energy = sum_q_cos_gr(n_gvec+1,iw)
        !           energy=0.d0
        !           energy=-old_energy

        jelnatoms = jel + natoms

        !     call CartesianToCrystal(kel(1,jel),s_cord(1,jelnatoms),1)
        !     s_cord(:,jelnatoms)=kel(:,jel)
        !     call TrigRecur(jelnatoms,jelnatoms,s_cord)
        !     s_cord(:,jelnatoms)=kel(:,jel)
        call TrigRecur(1, 1, kel)
!$omp parallel do default(shared) private(ig,qp,h,k,l)
        do ig = 1, n_gvec
            h = ind_g(ig, 1)
            k = ind_g(ig, 2)
            l = ind_g(ig, 3)
            qp = q(jelnatoms)*phsfac1(h, 1)*phsfac2(k, 1)*phsfac3(l, 1)
            sum_q_cos_gr(ig, iw) = sum_q_cos_gr(ig, iw) - dreal(qp)
            sum_q_sin_gr(ig, iw) = sum_q_sin_gr(ig, iw) - dimag(qp)
        end do

        !     call CartesianToCrystal(keln,s_cord(1,jelnatoms),1)
        !     s_cord(:,jelnatoms)=keln(:)
        !     call TrigRecur(jelnatoms,jelnatoms,s_cord)
        call TrigRecur(1, 1, keln)

!$omp parallel do default(shared)  private(ig,qp,h,k,l)
        do ig = 1, n_gvec
            h = ind_g(ig, 1)
            k = ind_g(ig, 2)
            l = ind_g(ig, 3)
            qp = q(jelnatoms)*phsfac1(h, 1)*phsfac2(k, 1)*phsfac3(l, 1)
            sum_q_cos_gr(ig, iw) = sum_q_cos_gr(ig, iw) + dreal(qp)
            sum_q_sin_gr(ig, iw) = sum_q_sin_gr(ig, iw) + dimag(qp)
        end do
!$omp parallel do default(shared) private(ig) reduction(+: new_energy)
        do ig = 1, n_gvec
            new_energy = new_energy + factor(ig)*(sum_q_cos_gr(ig, iw)**2.d0 + sum_q_sin_gr(ig, iw)**2.d0)
        end do
        new_energy = new_energy*8.d0*Pi/omega
        denergy = -sum_q_cos_gr(n_gvec + 1, iw) + new_energy ! The output is the energy  change
        energy = new_energy ! The output is the new energy
        sum_q_cos_gr(n_gvec + 1, iw) = new_energy !  Replace the old with new energy
    end subroutine EwaldUpdate

    subroutine Ewaldup1b(keln, jel, energy)
        implicit none
        double precision :: keln(3)
        double complex :: qp
        integer :: jelnatoms
        integer, intent(in) :: jel
        double precision, intent(out) :: energy
        integer :: ig, h, k, l
        ! with half g-space, G=0 is not in the list
        !       In input it is assumed that sum_q_cos and sum_q_sin do not
        !       contain the cos and sin corresponding to the electron coordinates
        !        namely instead of EwaldSum this arrays have to be initialized
        !        with EwaldSum1b
        jelnatoms = jel + natoms
        s_cord(:, jelnatoms) = keln(:)
        call TrigRecur(jelnatoms, jelnatoms, s_cord)

        !      call TrigRecur(1,1,keln)
        !      call TrigRecur_DFT(keln)

        energy = 0.d0
        !#ifdef __REDUCTION
        !!$omp parallel do shared(ind_g,n_gvec,q,phsfac1,phsfac2, &
        !! &phsfac3,factor,sum_q_cos,sum_q_sin,jelnatoms) private(ig,qp,h,k,l) reduction(+: energy)
        !#endif
        do ig = 1, n_gvec
            h = ind_g(ig, 1)
            k = ind_g(ig, 2)
            l = ind_g(ig, 3)
            qp = q(jelnatoms)*phsfac1(h, jelnatoms)*phsfac2(k, jelnatoms)*phsfac3(l, jelnatoms)
            !     energy=energy+(sum_q_cos(ig)*dreal(qp)+sum_q_sin(ig)*dimag(qp))/gg(ig)
            energy = energy + factor(ig)*(sum_q_cos(ig)*dreal(qp) + sum_q_sin(ig)*dimag(qp))
        end do
        !#ifdef __REDUCTION
        !!$omp end parallel do
        !#endif
        !       Adding the q=0 contribution of the short range ewald
        !       energy=energy+0.25d0/kappa**2*q(jelnatoms)*dble(natoms)

        energy = energy*16.d0*Pi/omega ! a factor two was missing
    end subroutine Ewaldup1b

    subroutine EwaldSum1b(rion)
        implicit none
        double precision, intent(in) :: rion(3, natoms)
        integer :: ig, h, k, l, p
        !           call CartesianToCrystal(rion,s_cord,natoms)

        s_cord(:, 1:natoms) = rion(:, 1:natoms)
        !           do p=1,natoms
        call TrigRecur(1, natoms, s_cord)
        !           enddo

        ! with half g-space, G=0 is not in the list
        ewaldion1b = 0.d0
        ewaldel1b = 0.d0
        do ig = 1, n_gvec
            h = ind_g(ig, 1); k = ind_g(ig, 2); l = ind_g(ig, 3)
            qphase(1:natoms) = q(1:natoms)*phsfac1(h, 1:natoms)* &
                    & phsfac2(k, 1:natoms)*phsfac3(l, 1:natoms)

            q_cos_gr(1:natoms) = dreal(qphase(1:natoms))
            q_sin_gr(1:natoms) = dimag(qphase(1:natoms))

            sum_q_cos(ig) = sum(q_cos_gr(1:natoms))
            sum_q_sin(ig) = sum(q_sin_gr(1:natoms))

            ewaldion1b = ewaldion1b + factor(ig)*(sum_q_cos(ig)**2 + sum_q_sin(ig)**2)
            ewaldel1b = ewaldel1b + factor(ig)
        end do
        ewaldion1b = ewaldion1b*8.d0*pi/omega/nelectrons
        ewaldel1b = ewaldel1b*8.d0*pi/omega - 2.d0*epsilon0*kappa/dsqrt(PI)

    end subroutine EwaldSum1b
end module Ewald
function cond_find(sab, sac, sbc)
    implicit none
    real*8 cond_find, sab, sac, sbc, cond_try, checkcond
!       a b c
    cond_try = (1.d0 - sac**2) - (sbc - sab*sac)**2/(1.d0 - sab**2)
    if (cond_try .gt. 0.d0) then
        cond_try = dsqrt(cond_try)
    else
        cond_try = dsqrt(abs(1.d0 - sbc**2)) ! 2d condition
    end if
    cond_find = cond_try

!       b a c
!       cond_try=(1.d0-sac**2)-(sbc-sab*sac)**2/(1.d0-sab**2)
    cond_try = (1.d0 - sbc**2) - (sac - sab*sbc)**2/(1.d0 - sab**2)

    if (cond_try .gt. 0.d0) then
        cond_try = dsqrt(cond_try)
    else
        cond_try = dsqrt(abs(1.d0 - sab**2))
    end if

    if (cond_try .lt. cond_find) cond_find = cond_try

!       c b a
!       cond_try=(1.d0-sac**2)-(sbc-sab*sac)**2/(1.d0-sab**2)
    cond_try = (1.d0 - sac**2) - (sab - sbc*sac)**2/(1.d0 - sbc**2)
    if (cond_try .gt. 0.d0) then
        cond_try = dsqrt(cond_try)
    else
        cond_try = dsqrt(abs(1.d0 - sac**2))
    end if

    if (cond_try .lt. cond_find) cond_find = cond_try

!       c a b
!       cond_try=(1.d0-sac**2)-(sbc-sab*sac)**2/(1.d0-sab**2)
    cond_try = (1.d0 - sbc**2) - (sab - sac*sbc)**2/(1.d0 - sac**2)
    if (cond_try .gt. 0.d0) then
        cond_try = dsqrt(cond_try)
    else
        cond_try = dsqrt(abs(1.d0 - sac**2))
    end if

    if (cond_try .lt. cond_find) cond_find = cond_try

!       a c b
!       cond_try=(1.d0-sac**2)-(sbc-sab*sac)**2/(1.d0-sab**2)
    cond_try = (1.d0 - sab**2) - (sbc - sac*sab)**2/(1.d0 - sac**2)
    if (cond_try .gt. 0.d0) then
        cond_try = dsqrt(cond_try)
    else
        cond_try = dsqrt(abs(1.d0 - sac**2))
    end if

    if (cond_try .lt. cond_find) cond_find = cond_try

!       b c a
!       cond_try=(1.d0-sac**2)-(sbc-sab*sac)**2/(1.d0-sab**2)
    cond_try = (1.d0 - sab**2) - (sac - sbc*sab)**2/(1.d0 - sbc**2)
    if (cond_try .gt. 0.d0) then
        cond_try = dsqrt(cond_try)
    else
        cond_try = dsqrt(abs(1.d0 - sac**2))
    end if

    if (cond_try .lt. cond_find) cond_find = cond_try
    checkcond = 1.d0 - sac**2 - sbc**2 - sab**2 + 2.d0*sab*sbc*sac
    if (checkcond .lt. 0.d0) cond_find = cond_find*2.d0
    ! checkcond = det(s_ij)>0 because s_ij= <v_i | v_j>, v_1=a/|a|, v_2=b/|b|, v_3=c/|c|. Thus checkcond<0 is never met.

    return
end
