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

module compute_efermi

    use allio, only: rank, qp
    use kpoints_mod, only: nk, wkp, wkp_down

    implicit none

    public
    integer :: smear_type ! smearing technique
    ! current choice: smear_type = 1 : Fermi-Dirac

    real(qp) :: efermi, efermido, ef_up, ef_do, metS, metSdo
    integer, dimension(:), allocatable :: number_particles, number_particles_do

    ! prcedure interfaces for the calculation of Fermi energy in the case of:
    ! !decoupled_run -> *_kaverage
    ! decoupled_run  -> *_kindependent
    interface efermiFix
        module procedure efermiFix_kaverage
        module procedure efermiFix_kindependent
    end interface efermiFix

    interface efermiMet
        module procedure efermiMet_kaverage
        module procedure efermiMet_kindependent
    end interface efermiMet

contains
    !
    ! Fermi energy in the case of fixed occupations
    subroutine efermiFix_kaverage(eig, occupations, bands, ef)
        implicit none
        integer, intent(in) :: bands
        real(8), intent(in) :: eig(bands, *), occupations(bands, *)
        real(qp), intent(out) :: ef
        integer :: ibnd, i, ikp
        real(qp) :: ef_tmp
        ! assuming that each k-point has the same (fixed) occupation
        ibnd = 0
        do i = 1, bands
            if (occupations(i, 1) .ne. 0.d0) ibnd = ibnd + 1
        end do
        ef = 0.d0
        if (ibnd .gt. 0) then
            if (ibnd + 1 .le. bands) then
                do ikp = 1, nk
                    ! take Efermi in the middle of the band gap
                    ef_tmp = (eig(ibnd, ikp) + eig(ibnd + 1, ikp))/2.d0
                    ef = ef + ef_tmp*wkp(ikp)
                end do
            else
                do ikp = 1, nk
                    ! take Efermi just above the HOMO level
                    ef_tmp = eig(ibnd, ikp) + 1d-6
                    ef = ef + ef_tmp*wkp(ikp)
                end do
            end if
        end if
        return
    end subroutine efermiFix_kaverage

    subroutine efermiFix_kindependent(eig, occupations, bands, ef)
        implicit none
        integer, intent(in) :: bands
        real(8), intent(in) :: eig(bands), occupations(bands)
        real(qp), intent(out) :: ef
        integer :: ibnd, i, ikp
        real(qp) :: ef_tmp
        ! assuming that each k-point has the same (fixed) occupation
        ibnd = 0
        ef = 0.d0
        do i = 1, bands
            if (occupations(i) .ne. 0.d0) ibnd = ibnd + 1
        end do
        if (ibnd .gt. 0) then ! the occupation can be empty
            if (ibnd + 1 .le. bands) then
                ef = (eig(ibnd) + eig(ibnd + 1))/2.d0
            else
                ef = eig(ibnd) + 1d-6
            end if
        end if
        return
    end subroutine efermiFix_kindependent
    !
    ! Fermi energy in case of smearing
    !
    subroutine efermiMet_kaverage(eig, bands, nel, dt, ef)

        implicit none
        ! input
        integer, intent(in) :: bands, nel
        real(8), intent(in) :: eig(bands, *)
        real(qp), intent(in) :: dt
        real(qp), intent(inout) :: ef
        ! local
        integer :: iter, ibnd, ikp
        real(qp) :: efmin, efmax, efsum, ntry, error
        integer, parameter :: maxit = 100
        real(qp), parameter :: eps = 1.d-14

        ! find energy bounds among k-points
        efmin = eig(1, 1)
        efmax = eig(bands, 1)
        do ikp = 1, nk
            if (eig(1, ikp) .lt. efmin) efmin = eig(1, ikp)
            if (eig(bands, ikp) .gt. efmax) efmax = eig(bands, ikp)
        end do
        efmin = efmin - dt*10.0_qp
        efmax = efmax + dt*10.0_qp
        error = (efmax - efmin)/2.0_qp
        iter = 0
        ! use bisection method to evaluate Fermi energy
        do while (error .gt. eps .and. iter .le. maxit)
            ntry = 0.0_qp
            ef = (efmin + efmax)/2.0_qp
            do ikp = 1, nk
                efsum = 0.0_qp
                do ibnd = 1, bands
                    efsum = efsum + smear(eig(ibnd, ikp), dt, ef, smear_type)
                end do
                ntry = ntry + wkp(ikp)*efsum
            end do
            if (ntry .gt. nel) then
                efmax = ef
            elseif (ntry .lt. nel) then
                efmin = ef
            end if
            error = max(efmax - efmin, 0.0_qp)
            iter = iter + 1
        end do
        return

    end subroutine efermiMet_kaverage

    subroutine efermiMet_kindependent(eig, bands, nel, dt, ef)

        implicit none
        ! input
        integer, intent(in) :: bands, nel
        real(8), intent(in) :: eig(bands)
        real(qp), intent(in) :: dt
        real(qp), intent(inout) :: ef
        ! local
        integer :: iter, ibnd, ikp
        real(qp) :: efmin, efmax, efsum, ntry, error
        integer, parameter :: maxit = 100
        real(qp), parameter :: eps = 1.d-14

        ! find energy bounds among k-points
        efmin = eig(1) - 10.0_qp*dt
        efmax = eig(bands) + 10.0_qp*dt
        error = (efmax - efmin)/2.0_qp
        iter = 0
        ! use bisection method to evaluate Fermi energy
        do while (error .gt. eps .and. iter .le. maxit)
            ntry = 0.0_qp
            ef = (efmin + efmax)/2.0_qp
            efsum = 0.0_qp
            do ibnd = 1, bands
                efsum = efsum + smear(eig(ibnd), dt, ef, smear_type)
            end do
            ntry = ntry + efsum
            if (ntry .gt. nel) then
                efmax = ef
            elseif (ntry .lt. nel) then
                efmin = ef
            end if
            error = max(efmax - efmin, 0.0_qp)
            iter = iter + 1
        end do
        return

    end subroutine efermiMet_kindependent

    ! -------------------
    ! Smearing functions
    !--------------------

    real(qp) function smear(e, dt, ef, ityp)

        implicit none

        real(8), intent(in) :: e
        real(qp), intent(in) :: dt, ef
        integer, intent(in) :: ityp
        real(qp) cost, cost2

        smear = 0.0_qp
        if (ityp .eq. 1) then
            if (e .gt. ef) then ! protection from overflow
                cost = exp(-(e - ef)/dt)
                cost2 = cost/(1.0_qp + cost)
            else
                cost2 = 1.0_qp/(exp((e - ef)/dt) + 1.0_qp)
            end if
        end if
        smear = cost2
        return

    end function smear

    ! -------------------
    ! Entropy functions
    !--------------------

    real(8) function smearS(e, dt, ef, ityp)

        implicit none

        real(8), intent(in) :: e
        real(qp), intent(in) :: dt, ef
        integer, intent(in) :: ityp
        real(qp) :: cost, mcost
        real(qp), parameter :: safemin = 1.0d-14

        smearS = 0.0_qp
        if (ityp .eq. 1) then
            !  Real protection from overflow
            if (e .lt. ef) then
                cost = exp((e - ef)/dt)
                cost = 1.0_qp/(1.0_qp + cost)
            else
                cost = exp(-(e - ef)/dt)
                cost = cost/(1.0_qp + cost)
            end if
            mcost = 1.0_qp - cost
            if ((abs(cost) .lt. safemin) .or. (abs(mcost) .lt. safemin)) then
                smearS = 0.0_qp
            else
                smearS = cost*log(cost) + mcost*log(mcost)
            end if
        end if
        return

    end function smearS

    ! -------------------
    ! Delta function
    !--------------------

    real(qp) function smearD(e, dt, ef, ityp)

        implicit none

        real(8), intent(in) :: e
        real(qp), intent(in) :: dt, ef
        integer, intent(in) :: ityp
        real(qp) :: cost, costexp

        smearD = 0.0_qp
        if (ityp .eq. 1) then
            cost = (e - ef)/dt
            ! protection from overflow
            costexp = exp(-abs(cost))
            smearD = costexp/(2.0_qp*costexp + 1.0_qp + costexp**2)
            !         if((cost.gt.0.)) then
            !            costexp=exp(-cost)
            !            smearD = costexp/(2.0_qp*costexp+1.0_qp+costexp**2)
            !!  multiplying num and den by costexp
            !!           smearD = 1.0_qp / ( 2.0_qp + exp(cost) + exp(-cost) )
            !         else
            !            costexp=exp(cost)
            !!  multiplying num and den by costexp
            !            smearD=costexp/(2.0_qp*costexp+costexp**2+1.0_qp)
            !!           smearD = 1.0_qp / ( 2.0_qp + exp(cost) + exp(-cost) )
            !         endif
        end if
        smearD = smearD/dt
        return

    end function smearD

end module compute_efermi
