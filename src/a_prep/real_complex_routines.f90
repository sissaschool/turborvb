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

! Subroutines to update the one-body part of the wave function.
! They are called in the routine "initialize_mats_new" only if
! the flag nbody_on is set to .true., if a one-body Jastrow is
! attached in front of the DFT wave function.

subroutine up_1body_complex(wf_scr, nbas_tot, indt, tlap, tgrad, jas_1body)
    use constants, only: zone
    implicit none
    integer, intent(in) :: indt, nbas_tot
    real(8), intent(in) :: tlap, tgrad(3), jas_1body
    complex(8), intent(out) :: wf_scr(nbas_tot*(indt + 5))
    integer :: ii
    wf_scr(1:nbas_tot) = wf_scr(1:nbas_tot)*dcmplx(jas_1body)
    wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot) = wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot)* &
                                                          dcmplx(jas_1body) + dcmplx(tlap)*wf_scr(1:nbas_tot)
    do ii = indt + 2, indt + 4
        wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot) = &
            wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot) &
            + (2.d0*zone*wf_scr(nbas_tot*(ii - 1) + 1:ii*nbas_tot) &
               *dcmplx(jas_1body)*dcmplx(tgrad(ii - indt - 1)) &
               + dcmplx(tgrad(ii - indt - 1))**2*wf_scr(1:nbas_tot))
    end do
    do ii = indt + 2, indt + 4
        wf_scr(nbas_tot*(ii - 1) + 1:nbas_tot*ii) = wf_scr(nbas_tot*(ii - 1) + 1:nbas_tot*ii)* &
                &dcmplx(jas_1body) + wf_scr(1:nbas_tot)*dcmplx(tgrad(ii - indt - 1))
    end do
    return
end subroutine up_1body_complex

subroutine up_1body(wf_scr, nbas_tot, indt, tlap, tgrad, jas_1body)
    implicit none
    integer, intent(in) :: indt, nbas_tot
    real(8), intent(in) :: tlap, tgrad(3), jas_1body
    real(8), intent(out) :: wf_scr(nbas_tot*(indt + 5))
    integer ii
    wf_scr(1:nbas_tot) = wf_scr(1:nbas_tot)*jas_1body
    wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot) = wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot)* &
                                                          jas_1body + tlap*wf_scr(1:nbas_tot)
    do ii = indt + 2, indt + 4
        wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot) = &
            wf_scr(nbas_tot*(indt + 4) + 1:(indt + 5)*nbas_tot) &
            + (2.d0*wf_scr(nbas_tot*(ii - 1) + 1:ii*nbas_tot) &
               *jas_1body*tgrad(ii - indt - 1) &
               + tgrad(ii - indt - 1)**2*wf_scr(1:nbas_tot))
    end do
    do ii = indt + 2, indt + 4
        wf_scr(nbas_tot*(ii - 1) + 1:nbas_tot*ii) = &
            wf_scr(nbas_tot*(ii - 1) + 1:nbas_tot*ii)*jas_1body + wf_scr(1:nbas_tot)*tgrad(ii - indt - 1)
    end do
    return
end subroutine up_1body

! Subroutines for matrix symmetrization. They are used in the subroutines
! "initialize_mats_new" and "uphamilt_new" in order to symmetrize overlap
! and Hamiltonian matrix. Useful to reduce roundoff especially in the
! __SCALAPACK version.

subroutine symmetrize_mat_complex(n, ld, mat_in)
#ifdef __SCALAPACK
    use descriptors
    use setup, only: desch
#endif
    implicit none

    integer, intent(in) :: ld, n
    complex(8) :: mat_in(ld, ld)
#ifdef __SCALAPACK
    complex(8), dimension(:, :), allocatable :: work
#endif
    integer :: i, j

#ifdef __SCALAPACK
    allocate (work(ld, ld))
    work = 0.d0
    if (descla(lambda_node_) > 0) then
        ! check dimensions
        if (n /= descla(la_n_)) &
            call errore(" symmetrize_mat_complex ", " wrong global dim n ", n)
        if (ld /= descla(nlax_)) &
            call errore(" symmetrize_mat_complex ", " wrong leading dim lda ", ld)
        ! symmetrize
        work = dconjg(mat_in)
        call pzgeadd('T', n, n, (0.5d0, 0.d0), work, 1, 1, desch, &
                     (0.5d0, 0.d0), mat_in, 1, 1, desch)
    end if
    deallocate (work)
#else
    do i = 1, n
        do j = i, n
            mat_in(i, j) = 0.5d0*(mat_in(i, j) + dconjg(mat_in(j, i)))
            mat_in(j, i) = dconjg(mat_in(i, j))
        end do
    end do
#endif
    return
end subroutine symmetrize_mat_complex

subroutine symmetrize_mat(n, ld, mat_in)
#ifdef __SCALAPACK
    use descriptors
    use setup, only: desch
#endif
    implicit none

    integer, intent(in) :: ld, n
    real(8) :: mat_in(ld, ld)
#ifdef __SCALAPACK
    real(8), dimension(:, :), allocatable :: work
#endif
    integer :: i, j

#ifdef __SCALAPACK
    allocate (work(ld, ld))
    work = 0.d0
    if (descla(lambda_node_) > 0) then
        ! check dimensions
        if (n /= descla(la_n_)) &
            call errore(" symmetrize_mat ", " wrong global dim n ", n)
        if (ld /= descla(nlax_)) &
            call errore(" symmetrize_mat ", " wrong leading dim lda ", ld)
        ! symmetrize
        work = mat_in
        call pdgeadd('T', n, n, 0.5d0, work, 1, 1, desch, &
                     0.5d0, mat_in, 1, 1, desch)
    end if
    deallocate (work)
#else
    do i = 1, n
        do j = i + 1, n
            mat_in(i, j) = 0.5d0*(mat_in(i, j) + mat_in(j, i))
            mat_in(j, i) = mat_in(i, j)
        end do
    end do
#endif
    return
end subroutine symmetrize_mat
