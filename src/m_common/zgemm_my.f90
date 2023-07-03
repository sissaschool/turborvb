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

subroutine ZGEMM_MY(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA&
        &, B, LDB, BETA, C, LDC, nproc, rank, comm_mpi)
    implicit none
#ifdef PARALLEL
    include 'mpif.h'
    integer ndim2, i, nu, nm, indr
#endif
    character*1 TRANSA, TRANSB
    integer M, N, K, LDA, LDB, LDC, nproc, rank, comm_mpi
    complex*16 ALPHA, BETA
    complex*16 A(LDA, *), B(LDB, *), C(LDC, *)
#ifdef PARALLEL
#ifdef __TEST
    integer dima, dimb, dimc
!imposing consistent input
    if (nproc .gt. 1) then
        if (transa .eq. 'N' .or. transa .eq. 'n') then
            dima = LDA*(K - 1) + M
        else
            dima = LDA*(M - 1) + K
        end if
        if (transb .eq. 'N' .or. transb .eq. 'n') then
            dimb = LDB*(N - 1) + K
        else
            dimb = LDB*(K - 1) + N
        end if
        dimc = LDC*(N - 1) + M
        call bcast_real(a, dima, 0, comm_mpi)
        call bcast_real(b, dimb, 0, comm_mpi)
        call bcast_real(c, dimc, 0, comm_mpi)
    end if
#endif
    ndim2 = 2*(LDC*(N - 1) + M)
    if (nproc .gt. 1) then
!             allocate(psip(ldc,N))
        nm = k/nproc
        if (nm*nproc .ne. k) nm = nm + 1
        indr = rank*nm + 1
        nu = nm
        if (indr + nm - 1 .gt. k) nu = k - indr + 1

        if (nu .gt. 0) then
            if ((transa .eq. 'N' .or. transa .eq. 'n') .and. (transb .eq. 'N' .or. transb .eq. 'n')) then
                if (rank .ne. 0) then
                    call zgemm('N', 'N', m, n, nu, alpha, a(1, indr), lda&
         &, b(indr, 1), ldb, (0.d0, 0.d0), c, ldc)
                else
                    call zgemm('N', 'N', m, n, nu, alpha, a(1, indr), lda&
         &, b(indr, 1), ldb, beta, c, ldc)
                end if
            elseif ((transa .eq. 'N' .or. transa .eq. 'n') .and.&
      &(transb .eq. 'T' .or. transb .eq. 't' .or. transb .eq. 'C' .or. transb .eq. 'c')) then
                if (rank .ne. 0) then
                    call zgemm('N', transb, m, n, nu, alpha, a(1, indr), lda&
         &, b(1, indr), ldb, (0.d0, 0.d0), c, ldc)
                else
                    call zgemm('N', transb, m, n, nu, alpha, a(1, indr), lda&
         &, b(1, indr), ldb, beta, c, ldc)
                end if
            elseif ((transa .eq. 'T' .or. transa .eq. 't' .or. transa .eq. 'c' .or&
      &. transa .eq. 'C') .and. (transb .eq. 'N' .or. transb .eq. 'n')) then
                if (rank .ne. 0) then
                    call zgemm(transa, 'N', m, n, nu, alpha, a(indr, 1), lda&
         &, b(indr, 1), ldb, (0.d0, 0.d0), c, ldc)
                else
                    call zgemm(transa, 'N', m, n, nu, alpha, a(indr, 1), lda&
         &, b(indr, 1), ldb, beta, c, ldc)
                end if
            elseif ((transa .eq. 'T' .or. transa .eq. 't' .or. transa .eq. 'c' .or&
      &. transa .eq. 'C') .and. (transb .eq. 'T' .or. transb .eq. 't' .or&
      &. transb .eq. 'c' .or. transb .eq. 'C')) then
                if (rank .ne. 0) then
                    call zgemm(transa, transb, m, n, nu, alpha, a(indr, 1), lda&
         &, b(1, indr), ldb, (0.d0, 0.d0), c, ldc)
                else
                    call zgemm(transa, transb, m, n, nu, alpha, a(indr, 1), lda&
         &, b(1, indr), ldb, beta, c, ldc)
                end if
            end if
        elseif (rank .ne. 0) then
            do i = 1, n
                c(1:m, i) = 0.d0
            end do
        end if
        if (rank .ne. 0 .and. ldc .gt. m) then
            do i = 1, n - 1
                c(m + 1:ldc, i) = 0.d0
            end do
        end if
        call reduce_base_real(ndim2, c, comm_mpi, -1)
    else
#endif
        call ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA&
                &, B, LDB, BETA, C, LDC)
#ifdef PARALLEL
    end if
#endif
    return
end
