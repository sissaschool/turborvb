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

!> @brief Generic matrix-matrix multiplication supporting mixed real/complex operations
!> @details This subroutine performs matrix-matrix multiplication C = alpha*A*B + beta*C
!> with support for various combinations of real and complex matrices. The operation
!> type is specified by TRANSA and TRANSB parameters which can be:
!> - 'NR': Real normal (no transpose)
!> - 'NC': Complex normal (no transpose)
!> - 'TR': Real transpose
!> - 'TC': Complex transpose (no conjugate)
!> - 'CC': Complex adjoint (transpose + conjugate)
!> - 'SC': Simple conjugate (no transpose)
!> The subroutine handles mixed cases where one matrix is real and the other is complex
!> by decomposing the operation into separate real matrix multiplications.
!> @param[in] TRANSA Operation on matrix A ('NR', 'NC', 'TR', 'TC', 'CC', 'SC')
!> @param[in] TRANSB Operation on matrix B ('NR', 'NC', 'TR', 'TC', 'CC', 'SC')
!> @param[in] M Number of rows in matrix C
!> @param[in] N Number of columns in matrix C
!> @param[in] K Number of columns in A and rows in B
!> @param[in] ALPHA Complex scaling factor for A*B
!> @param[in] A Input matrix A
!> @param[in] LDA Leading dimension of matrix A
!> @param[in] B Input matrix B
!> @param[in] LDB Leading dimension of matrix B
!> @param[in] BETA Complex scaling factor for C
!> @param[in,out] C Input/output matrix C
!> @param[in] LDC Leading dimension of matrix C
subroutine gemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    use allio, only: rank
    implicit none
    complex*16 alpha, beta
    integer M, N, K, LDA, LDB, LDC, i, j
    real*8 A(LDA, *), B(LDB, *), C(LDC, *)
    real*8, dimension(:, :), allocatable :: Creal, Ccomp
    character*2 TRANSA, TRANSB
    character*1 TRANSAB, TRANSBB
    logical yesr

    !        Generic matrix matrix multiplication not dependent on the type
    !        complex or real of the matrices.
    !        Always capital letter assumed for simplicity
    !          Trans= 'NR'       Real normal
    !                 'NC'       Complex normal
    !          Trans= 'TR'       Real Transpose
    !          Trans= 'TC'       Complex Transpose (no Conjugate)
    !          Trans= 'CC'       Complex Adjoint (Transpose+Conjugate)
    !          Trans= 'SC'       Simple  Conjugate (no Transpose)
    !         alpha and beta are complex constants
    !         NB in case A,B,C are complex matrices the leading dimension
    !         in input has to be obviously equal to twice the number of complex
    !         components in each column of the matrix.
    ! Does not work for mixed cases!

    ! if in input they are real some simplification occurs.
    if (aimag(alpha) .eq. 0.d0 .and. aimag(beta) .eq. 0.d0) then
        yesr = .true.
    else
        yesr = .false.
    end if

    if (TRANSA .eq. 'NR' .and. TRANSB .eq. 'NR') then
        if (yesr) then
            call dgemm('N', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
        else
            !          C is complex but AxB is real
            allocate (creal(M, N))
            call dgemm('N', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Creal, M)
            call putcrealinc(M, N, LDC, alpha, beta, C, Creal)
            deallocate (creal)
        end if
    elseif (TRANSA .eq. 'NR' .and. TRANSB .eq. 'TR') then
        if (yesr) then
            call dgemm('N', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
        else
            !          C is complex but AxB is real
            allocate (creal(M, N))
            call dgemm('N', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Creal, M)
            call putcrealinc(M, N, LDC, alpha, beta, C, Creal)
            deallocate (creal)
        end if
    elseif (TRANSA .eq. 'TR' .and. TRANSB .eq. 'NR') then
        if (yesr) then
            call dgemm('T', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
        else
            !          C is complex but AxB is real
            allocate (creal(M, N))
            call dgemm('T', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Creal, M)
            call putcrealinc(M, N, LDC, alpha, beta, C, Creal)
            deallocate (creal)
        end if
    elseif (TRANSA .eq. 'TR' .and. TRANSB .eq. 'TR') then
        if (yesr) then
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
        else
            !          C is complex but AxB is real
            allocate (creal(M, N))
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Creal, M)
            call putcrealinc(M, N, LDC, alpha, beta, C, Creal)
            deallocate (creal)
        end if
    elseif (TRANSA .eq. 'NC' .and. TRANSB .eq. 'NC') then
        call zgemm('N', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'NC' .and. TRANSB .eq. 'TC') then
        call zgemm('N', 'T', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'NC' .and. TRANSB .eq. 'CC') then
        call zgemm('N', 'C', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'NC' .and. TRANSB .eq. 'SC') then
        ! We use that two times the conjugate is the identity
        call conjb
        call zgemm('N', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        call conjb
    elseif (TRANSA .eq. 'TC' .and. TRANSB .eq. 'NC') then
        call zgemm('T', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'TC' .and. TRANSB .eq. 'TC') then
        call zgemm('T', 'T', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'TC' .and. TRANSB .eq. 'CC') then
        call zgemm('T', 'C', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'TC' .and. TRANSB .eq. 'SC') then
        ! We use that two times the conjugate is the identity
        call conjb
        call zgemm('T', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        call conjb
    elseif (TRANSA .eq. 'CC' .and. TRANSB .eq. 'NC') then
        call zgemm('C', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'CC' .and. TRANSB .eq. 'TC') then
        call zgemm('C', 'T', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
    elseif (TRANSA .eq. 'CC' .and. TRANSB .eq. 'CC') then
        call zgemm('C', 'C', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        !          Now non trivial cases : A is real and  B is complex
    elseif (TRANSA .eq. 'CC' .and. TRANSB .eq. 'SC') then
        ! We use that two times the conjugate is the identity
        call conjb
        call zgemm('C', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        call conjb
    elseif (TRANSA .eq. 'SC' .and. TRANSB .eq. 'NC') then
        call conja
        call zgemm('N', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        call conja
    elseif (TRANSA .eq. 'SC' .and. TRANSB .eq. 'TC') then
        call conja
        call zgemm('N', 'T', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        call conja
    elseif (TRANSA .eq. 'SC' .and. TRANSB .eq. 'CC') then
        call conja
        call zgemm('N', 'C', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        call conja
    elseif (TRANSA .eq. 'SC' .and. TRANSB .eq. 'SC') then
        ! We use that two times the conjugate is the identity
        call conja
        call conjb
        call zgemm('N', 'N', M, N, K, alpha, A, LDA/2, B, LDB/2, &
                &                   beta, C, LDC/2)
        call conjb
        call conja
    elseif (TRANSA .eq. 'NR' .and. TRANSB .eq. 'NC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        if (yesr) then
            !           real part
            call dgemm('N', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, 1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSA .eq. 'NR' .and. TRANSB .eq. 'TC') then
        if (yesr) then
            !           real part
            call dgemm('N', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   beta, C, LDC)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   beta, C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, 1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSA .eq. 'NR' .and. TRANSB .eq. 'CC') then
        if (yesr) then
            !           real part
            call dgemm('N', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, -real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, -1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSA .eq. 'NR' .and. TRANSB .eq. 'SC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        call conjb
        if (yesr) then
            !           real part
            call dgemm('N', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, 1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
        call conjb
    elseif (TRANSA .eq. 'TR' .and. TRANSB .eq. 'NC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        if (yesr) then
            !           real part
            call dgemm('T', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, 1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSA .eq. 'TR' .and. TRANSB .eq. 'TC') then
        if (yesr) then
            !           real part
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSA .eq. 'TR' .and. TRANSB .eq. 'CC') then
        if (yesr) then
            !           real part
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, -real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, -1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSA .eq. 'TR' .and. TRANSB .eq. 'TC') then
        if (yesr) then
            !           real part
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSA .eq. 'TR' .and. TRANSB .eq. 'SC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        call conjb
        if (yesr) then
            !           real part
            call dgemm('T', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, real(alpha), A, LDA, B(2, 1), LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, 1.d0, A, LDA, B(2, 1), LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
        call conjb
    elseif (TRANSB .eq. 'NR' .and. TRANSA .eq. 'NC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        if (yesr) then
            !           real part
            call dgemm('N', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, 1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSB .eq. 'NR' .and. TRANSA .eq. 'TC') then
        if (yesr) then
            !           real part
            call dgemm('T', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, 1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSB .eq. 'NR' .and. TRANSA .eq. 'CC') then
        if (yesr) then
            !           real part
            call dgemm('T', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, -real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'N', M, N, K, -1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSB .eq. 'NR' .and. TRANSA .eq. 'SC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        call conja
        if (yesr) then
            !           real part
            call dgemm('N', 'N', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'N', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'N', M, N, K, 1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
        call conja
    elseif (TRANSB .eq. 'TR' .and. TRANSA .eq. 'NC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        if (yesr) then
            !           real part
            call dgemm('N', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, 1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSB .eq. 'TR' .and. TRANSA .eq. 'TC') then
        if (yesr) then
            !           real part
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, 1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSB .eq. 'TR' .and. TRANSA .eq. 'CC') then
        if (yesr) then
            !           real part
            call dgemm('T', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, -real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !           real part
            allocate (ccomp(2*M, N))
            call dgemm('T', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('T', 'T', M, N, K, -1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
    elseif (TRANSB .eq. 'TR' .and. TRANSA .eq. 'SC') then
        ! Here I use that the product of a real times a complex matrix
        ! is equivalent to two realxreal  matrix products.
        call conja
        if (yesr) then
            !           real part
            call dgemm('N', 'T', M, N, K, real(alpha), A, LDA, B, LDB, &
                    &                   real(beta), C, LDC)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, real(alpha), A(2, 1), LDA, B, LDB, &
                    &                   real(beta), C(2, 1), LDC)
        else
            !          real part
            allocate (ccomp(2*M, N))
            call dgemm('N', 'T', M, N, K, 1.d0, A, LDA, B, LDB, &
                    &                   0.d0, Ccomp, 2*M)
            !           imaginary  part
            call dgemm('N', 'T', M, N, K, 1.d0, A(2, 1), LDA, B, LDB, &
                    &                   0.d0, Ccomp(2, 1), 2*M)
            call putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
            deallocate (ccomp)
        end if
        call conja
    else
        write (6, *) ' Case not treated in gemm !!! ', TRANSA, TRANSB
    end if
    return
contains
    !> @brief Conjugate matrix B by negating imaginary parts
    !> @details This subroutine performs complex conjugation of matrix B by negating
    !> the imaginary components (even-indexed rows) while keeping real parts unchanged.
    !> The matrix is modified in-place.
    !> @note Matrix B is assumed to be stored in interleaved real/imaginary format
    !> where even rows contain imaginary parts and odd rows contain real parts.
    subroutine conjb
        implicit none
        integer i, j
        do j = 1, k
            do i = 1, n
                b(2*i, j) = -b(2*i, j)
            end do
        end do
    end subroutine conjb

    !> @brief Conjugate matrix A by negating imaginary parts  
    !> @details This subroutine performs complex conjugation of matrix A by negating
    !> the imaginary components (even-indexed rows) while keeping real parts unchanged.
    !> The matrix is modified in-place.
    !> @note Matrix A is assumed to be stored in interleaved real/imaginary format
    !> where even rows contain imaginary parts and odd rows contain real parts.
    subroutine conja
        implicit none
        integer i, j
        do j = 1, k
            do i = 1, m
                a(2*i, j) = -a(2*i, j)
            end do
        end do
    end subroutine conja
end subroutine gemm

!> @brief Put real matrix result into complex matrix with scaling
!> @details This helper subroutine combines a real matrix result with a complex matrix
!> using the formula: C = beta*C + alpha*Creal. It is used when the matrix multiplication
!> result is real but needs to be stored in a complex matrix format.
!> @param[in] M Number of rows
!> @param[in] N Number of columns
!> @param[in] LDC Leading dimension of matrix C
!> @param[in] alpha Complex scaling factor
!> @param[in] beta Complex scaling factor
!> @param[in,out] C Complex output matrix
!> @param[in] Creal Real input matrix
subroutine putcrealinc(M, N, LDC, alpha, beta, C, Creal)
    implicit none
    integer M, N, LDC
    complex*16 C(LDC/2, *), alpha, beta
    real*8 Creal(M, N)
    C(1:M, 1:N) = beta*C(1:M, 1:N) + alpha*Creal(1:M, 1:N)
    return
end subroutine putcrealinc

!> @brief Put complex matrix result into complex matrix with scaling
!> @details This helper subroutine combines two complex matrices using the formula:
!> C = beta*C + alpha*Ccomp. It is used when both input and output matrices are complex.
!> @param[in] M Number of rows
!> @param[in] N Number of columns
!> @param[in] LDC Leading dimension of matrix C
!> @param[in] alpha Complex scaling factor
!> @param[in] beta Complex scaling factor
!> @param[in,out] C Complex output matrix
!> @param[in] Ccomp Complex input matrix
subroutine putccompinc(M, N, LDC, alpha, beta, C, Ccomp)
    implicit none
    integer M, N, LDC
    complex*16 C(LDC/2, *), Ccomp(M, *), alpha, beta
    C(1:M, 1:N) = beta*C(1:M, 1:N) + alpha*Ccomp(1:M, 1:N)
    return
end subroutine putccompinc
