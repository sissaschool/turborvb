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

subroutine upinvhop(nelorb, nelorb5, nelc, ainv, winv &
        &, psiln, psisn, epst, psip, nelup, neldo, ainvs, winvs, update_psi)

    use Constants, only: ipc, zmone, zone, ipf, yes_ontarget
    use allio, only: nelup_mat

    implicit none
    !SSS attenzione alla dichiarazione di ainv!
    integer nel, j, i, nelc, nelcdo, nelorb, nelorb5, nelup, neldo, ierr
    real(8) epst, psiln, psisn
    real(8) ainv(nelup_mat*ipc, nelup_mat), ainvs(nelup_mat*ipc, 2), &
        winvs(*), winv(nelorb5*ipc, *)
    real(8) gr
    ! scratch vectors for complex wf
    real(8) :: psip(ipc*(nelup_mat + 3*(ipf - 1)*nelup_mat))
    complex(8) gc, phs
    logical update_psi

    ! calculation prefactor
    ! nelc is the # electrons that moves nelc<=nelup up spin,
    ! otherwise down spin
    ! jold is its position
    ! indtc is the direction where it moves according to table ivi
    ! here winv(nelorb5,nel) are the wavefunction orbitals

    ! compute the wavefunction in the new position r
    ! update new functions

    nel = neldo + nelup
    if (ipf .eq. 2) then
        !In the pfaffian case the update is done using a matrix matrix product
        !of an auxiliary vector by itself, recovering the expression in the
        !book (Becca, Sorella)
        !    if(update_psi) winv(1:nelorb*ipc,nelc)=winvs(1:nelorb*ipc)
        if (update_psi) call dcopy_vec(nelorb*ipc, winvs, winv(1, nelc))
        if (ipc .eq. 2) then
#ifdef _OFFLOAD
!$omp target update from(ainvs(2*nelc-1:2*nelc,2)) if(yes_ontarget)
#endif
            gc = dcmplx(ainvs(nelc*2 - 1, 2) + 1.d0, ainvs(nelc*2, 2))
            if (update_psi) then
                phs = log(gc)
                psisn = psisn + aimag(phs)
                psiln = psiln + real(phs)
            end if
        else
#ifdef _OFFLOAD
!$omp target update from(ainvs(nelc:nelc,2)) if(yes_ontarget)
#endif
            gr = ainvs(nelc, 2) + 1.d0
            if (update_psi) then
                psisn = psisn*sign(1.d0, gr)
                psiln = psiln + log(abs(gr))
            end if
        end if
        !    psip(1:nelup_mat*ipc)=ainv(1:nelup_mat*ipc,nelc)
        !    psip(nelup_mat*ipc+1:2*nelup_mat*ipc)=-ainvs(1:nelup_mat*ipc,2)
        !    psip(nelup_mat*ipc*2+1:3*nelup_mat*ipc)=-ainvs(1:nelup_mat*ipc,2)
        !    psip(nelup_mat*ipc*3+1:4*nelup_mat*ipc)=-ainv(1:nelup_mat*ipc,nelc)

#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
        do i = 1, ipc*nelup_mat
            psip(i) = ainv(i, nelc)
            psip(nelup_mat*ipc + i) = -ainvs(i, 2)
            psip(nelup_mat*ipc*2 + i) = -ainvs(i, 2)
            psip(nelup_mat*ipc*3 + i) = -ainv(i, nelc)
        end do

        if (ipc .eq. 2) then
            gc = zone/gc
            !call zgemm ('N','T',nelup_mat, nelup_mat,2, gc, psip, nelup_mat, psip(nelup_mat*ipc*2+1)&
            ! &, nelup_mat, zone, ainv, nelup_mat)
            call zger2(nelup_mat, nelup_mat, gc, psip, nelup_mat, psip(4*nelup_mat + 1), nelup_mat, ainv, nelup_mat)
        else
            gr = 1.d0/gr
            !       Symmetric update deprecated as the antisymmetrization may be slow:
            !       unthreaded, unvectorized, cash limited...
            !        call dger(nelup_mat,nelup_mat,2*g(1),psip,1,psip(2*nelup_mat+1),1,ainv,nelup_mat)
            !         do i=1,nelup_mat
            !           ainv(i,i)=0.d0
            !           do j=i+1,nelup_mat
            !           ainv(i,j)=0.5d0*(ainv(i,j)-ainv(j,i))
            !           ainv(j,i)=-ainv(i,j)
            !           enddo
            !         enddo
            !        call dgemm ('N','T',nelup_mat, nelup_mat, 2, g(1), psip&
            ! &, nelup_mat, psip(nelup_mat*2+1), nelup_mat, 1.d0, ainv, nelup_mat)
            call dger2(nelup_mat, nelup_mat, gr, psip, nelup_mat, psip(2*nelup_mat + 1), nelup_mat, ainv, nelup_mat)
        end if

    else
        if (nelc .le. nelup) then ! up spin electrons
            ! update basis orbitals
            !       if(update_psi) winv(1:nelorb*ipc,nelc)=winvs(1:nelorb*ipc)
            if (update_psi) call dcopy_vec(nelorb*ipc, winvs, winv(1, nelc))
            ! update ratio

            if (ipc .eq. 2) then
#ifdef  _OFFLOAD
!$omp target update from(ainvs(2*nelc-1:2*nelc,2)) if(yes_ontarget)
#endif
                gc = dcmplx(ainvs(nelc*2 - 1, 2) + 1.d0, ainvs(nelc*2, 2))
                if (update_psi) then
                    phs = log(gc)
                    psisn = psisn + aimag(phs)
                    psiln = psiln + real(phs)
                end if
            else
#ifdef  _OFFLOAD
!$omp target update from(ainvs(nelc:nelc,2)) if(yes_ontarget)
#endif
                gr = ainvs(nelc, 2) + 1.d0
                if (update_psi) then
                    psisn = psisn*sign(1.d0, gr)
                    psiln = psiln + log(abs(gr))
                end if
            end if
            !       psip(1:nelup*ipc)=ainv(1:nelup*ipc,nelc)
            call dcopy_vec(nelup*ipc, ainv(1, nelc), psip)
            if (ipc .eq. 2) then
                gc = zmone/gc
                call zgeru_(nelup, nelup, gc, psip, 1, ainvs(1, 2), 1, ainv, nelup)
            else
                gr = -1.d0/gr
                call dger_(nelup, nelup, gr, psip, 1, ainvs(1, 2), 1, ainv, nelup)
            end if
        else ! down spin electrons
            !       if(update_psi) winv(1:nelorb*ipc,nelc)=winvs(1:nelorb*ipc)
            if (update_psi) call dcopy_vec(nelorb*ipc, winvs, winv(1, nelc))
            nelcdo = nelc - nelup
            if (ipc .eq. 2) then
#ifdef  _OFFLOAD
!$omp target update from(ainvs(2*nelcdo-1:2*nelcdo,2)) if(yes_ontarget)
#endif
                gc = dcmplx(ainvs(2*nelcdo - 1, 2) + 1.d0, ainvs(2*nelcdo, 2))
                if (update_psi) then
                    phs = log(gc)
                    psisn = psisn + aimag(phs)
                    psiln = psiln + real(phs)
                end if
            else
#ifdef  _OFFLOAD
!$omp target update from(ainvs(nelcdo:nelcdo,2)) if(yes_ontarget)
#endif
                gr = ainvs(nelcdo, 2) + 1.d0
                if (update_psi) then
                    psisn = psisn*sign(1.d0, gr)
                    psiln = psiln + log(abs(gr))
                end if
            end if

            if (ipc .eq. 2) then
                call up_ainv_down(psip, ainv, nelup, nelcdo)
                gc = zmone/gc
                call zgeru_(nelup, nelup, gc, ainvs(1, 2), 1, psip, 1, ainv, nelup)
            else
#ifdef  _OFFLOAD
!$omp target teams  distribute parallel do if(yes_ontarget)
#endif
                do i = 1, nelup
                    psip(i) = ainv(nelcdo, i)
                end do
                gr = -1.d0/gr
                call dger_(nelup, nelup, gr, ainvs(1, 2), 1, psip, 1, ainv, nelup)
            end if
        end if
    end if ! endif ipf=2
    return
end subroutine upinvhop
